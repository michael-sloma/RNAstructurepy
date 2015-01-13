#include <stdio.h>
#include <stdlib.h>
#include "TurboFold_object.h"
#include "../RNA_class/RNA.h"
#include <string>
#include "../src/phmm/phmm_aln.h"
#include "../src/phmm/structure/structure_object.h"
#include "../src/phmm/utils/xmath/matrix/matrix.h"

#ifdef COMPILE_SMP
	#include "../src/phmm/utils/mutex/ansi_mutex.h"
	#include "../src/phmm/utils/ansi_thread/ansi_thread.h"
#endif

#include "TurboFold_thread.h"

#ifdef __WINDOWS
	#include <float.h>
	#define isnan _isnan
#endif

TurboFold::TurboFold(const char fasta_fp[])
{
	this->err_code = 0;

	this->saves = NULL;

	progress=NULL;

	this->sequences = new vector<t_structure*>();
	this->folders = new vector<RNA*>();
	char _fasta_fp[1000];
	strcpy(_fasta_fp, fasta_fp);
	vector<t_structure*>* new_sequences = t_structure::read_multi_seq(_fasta_fp);

	for(int i_str = 0; i_str < new_sequences->size(); i_str++)
	{
		sequences->push_back(new t_structure(new_sequences->at(i_str)));
		RNA* new_folder = new RNA(&(new_sequences->at(i_str)->nucs[1]));

		if(new_folder->GetErrorCode() != 0)
		{
			this->err_code = RNA_LIB_RNA_CONSTRUCTOR_ERROR;
			return;
		}

		folders->push_back(new_folder);
                //printf("Added %s (%d nucs)\n", new_sequences->at(i_str)->ctlabel, folders->at(i_str)->GetSequenceLength());
		delete(new_sequences->at(i_str));
	} 
	new_sequences->clear();
	delete(new_sequences);

	if(sequences->size() == 0)
	{
		printf("Need at least 1 sequence to predict structure for.\n");
                return;
	}

	// Allocate extrinsic info.
        this->extrinsic_info_list = new vector<t_matrix*>();

        for(int i_seq = 0; i_seq < sequences->size(); i_seq++)
        {
                this->extrinsic_info_list->push_back(new t_matrix(this->sequences->at(i_seq)->numofbases, this->sequences->at(i_seq)->numofbases, true));
	} // i_seq loop. 

	this->refolding_threads = NULL;
	this->turbofold_threads_mutex = NULL;

	if(this->initialize_alignment_information() != 0)
	{
                this->err_code = ALIGNMENT_INFO_COMPUTATION_ERROR;
                return;
	}
}

TurboFold::TurboFold(vector<string>* _sequences, vector<string>* _saves)
{
	this->err_code = 0;

	// Do a simple sanity check.
	if(_sequences->size() !=  _saves->size())
	{
		this->err_code = CONSTRUCTOR_ERROR;
		return;
	}

	this->sequences = new vector<t_structure*>();
	this->folders = new vector<RNA*>();

	progress=NULL;

	// Allocate saves if it is supplied in the arguments.
	if(_saves != NULL)
	{
		this->saves = new vector<char*>();
	}
	else
	{
		this->saves = NULL;
	}

	this->extrinsic_info_list = new vector<t_matrix*>();
	
	for(int i_seq = 0; i_seq < _sequences->size(); i_seq++)
	{
		char cur_seq_fp[1000];
		strcpy(cur_seq_fp, _sequences->at(i_seq).c_str());
		this->sequences->push_back(new t_structure(cur_seq_fp));
		this->folders->push_back(new RNA(cur_seq_fp, 2)); 

		if(_saves != NULL)
		{
			char* cur_save_fp = (char*)malloc(sizeof(char) * (_saves->at(i_seq).length() + 2));	
			strcpy(cur_save_fp, _saves->at(i_seq).c_str());
			this->saves->push_back(cur_save_fp);
		}

		// Push back a symmetric matrix to store the extrinsic information for current sequence.
		this->extrinsic_info_list->push_back(new t_matrix(this->sequences->at(i_seq)->numofbases, this->sequences->at(i_seq)->numofbases, true));
	} // i_seq loop.

/*
        for(int i_seq = 0; i_seq < _sequences->size(); i_seq++)
        {
		printf("%d\n", this->folders->at(i_seq)->GetSequenceLength());
		for(int i_nuc = 1; i_nuc <= this->folders->at(i_seq)->GetSequenceLength(); i_nuc++)
		{
			printf("%c", this->folders->at(i_seq)->GetNucleotide(i_nuc));
		}
		printf("\n");
	}
*/

        this->refolding_threads = NULL;
        this->turbofold_threads_mutex = NULL;

        if(this->initialize_alignment_information() != 0)
        {
                this->err_code = ALIGNMENT_INFO_COMPUTATION_ERROR;
                return;
        }
}

TurboFold::~TurboFold()
{
	int n_seq = this->sequences->size();

	if(this->saves != NULL)
	{
		for(int i_seq = 0; i_seq < n_seq; i_seq++)
		{
			free(this->saves->at(i_seq));
		} // i_seq
		this->saves->clear();
		delete(this->saves);
	}

	for(int i_seq = 0; i_seq < n_seq; i_seq++)
	{
		delete(this->folders->at(i_seq));
		delete(this->extrinsic_info_list->at(i_seq));
	} // i_seq

	this->folders->clear();
	this->extrinsic_info_list->clear();
	delete(this->folders);
	delete(this->extrinsic_info_list);
	
	// Free the alignment information.
        for(int i_seq1 = 0; i_seq1 < n_seq; i_seq1++)
        {
                free(this->similarities[i_seq1]);
        } // i_seq1 loop
	free(this->similarities);

        for(unsigned int i_seq1 = 0; i_seq1 < n_seq; i_seq1++)
        {
                for(unsigned int i_seq2 = i_seq1+1; i_seq2 < n_seq; i_seq2++)
                {
                        if(i_seq1 != i_seq2)
                        {
                                t_phmm_aln* phmm_aln = new t_phmm_aln(sequences->at(i_seq1), sequences->at(i_seq2));

                                // Free alignment mapping probabilities.
                                //this->aln_mapping_probs[i_seq1][i_seq2] = (double**)malloc(sizeof(double*) * (this->sequences->at(i_seq1)->numofbases + 2));                          
                                for(int i = 1; i <= this->sequences->at(i_seq1)->numofbases; i++)
                                {
                                        int min_k = aln_env_results[i_seq1][i_seq2]->low_limits[i];
                                        int max_k = aln_env_results[i_seq1][i_seq2]->high_limits[i];
                                        ///this->aln_mapping_probs[i_seq1][i_seq2][i] = (double*)malloc(sizeof(double) * (max_k - min_k + 2));
                                        //this->aln_mapping_probs[i_seq1][i_seq2][i] -= min_k;
                                        this->aln_mapping_probs[i_seq1][i_seq2][i] += min_k;
                                        free(this->aln_mapping_probs[i_seq1][i_seq2][i]);
                                } // i loop
                                free(this->aln_mapping_probs[i_seq1][i_seq2]);

                                // Free alignment envelope result.
                                phmm_aln->free_aln_env_result(this->aln_env_results[i_seq1][i_seq2]);
                                delete(phmm_aln);
                        }
                        else
                        {
                        }
                } // i_seq2 loop.

                aln_env_results[i_seq1] += i_seq1;
                free(aln_env_results[i_seq1]);

                aln_mapping_probs[i_seq1] += i_seq1;
                free(this->aln_mapping_probs[i_seq1]);
        } // i_seq1 loop.
	

	// sequences are used in memory free'ing. Must make sure it is free'ed at the very last step.
        for(int i_seq = 0; i_seq < n_seq; i_seq++)
        {
                delete(this->sequences->at(i_seq));
        } // i_seq

        this->sequences->clear();
        delete(this->sequences);
}


//Set the maximum pairing distance for all sequences.
int TurboFold::SetMaxPairingDistance(int distance) {
	int error;

	for(int i_seq = 0; i_seq < this->sequences->size(); i_seq++) {

		error = folders->at(i_seq)->ForceMaximumPairingDistance(distance);
		if (error!=0) {
			err_code = RNALIB_SETDISTANCE_ERROR;
			return RNALIB_SETDISTANCE_ERROR;
		}

	}
	return 0;



}




int TurboFold::fold(double gamma, int _n_iterations, int _n_parallel_pfunctions)
{
#ifdef COMPILE_SMP
	// Allocate the mutex.
	this->turbofold_threads_mutex = new t_ansi_mutex();
	this->n_parallel_pfunctions = _n_parallel_pfunctions;
#else
	this->turbofold_threads_mutex = NULL;
	this->n_parallel_pfunctions = 0;
#endif



	// Generate the alignment information.
	//
	this->n_iterations = _n_iterations;
	
	// Initialize the loops.
	for(int i_iter = 0; i_iter <= n_iterations; i_iter++)
	{

		//Add a coarse update of progress:
		if (progress!=NULL) {

			progress->update((int)((100.0*((double) i_iter))/((double) n_iterations+1)));
		}

		// Set the extrinsic information for each sequence.
		if(i_iter == 0)
		{
			// Initialize the extrinsic information to 1.0 for all the sequences.
			for(int i_seq = 0; i_seq < this->sequences->size(); i_seq++)
			{
				for(int i = 1; i <= this->sequences->at(i_seq)->numofbases; i++)
				{
					for(int j = i+1; j <= this->sequences->at(i_seq)->numofbases; j++)
					{
						//this->SetExtrinsic(i_seq, i, j, 1.0f); 
						this->folders->at(i_seq)->SetExtrinsic(i, j, 1.0f);
					} // j loop
				} // i loop 		
			} // i_Seq loop	
		}
		else
		{
			// Compute the extrinsic information using the base pairing probabilities just computed for each sequence.
			if(this->generate_extrinsic_information(i_iter, gamma) != 0)
			{
				return(this->err_code);
			}	 	
		
			// Set ExtrinsicInformation(onst int i, const int j, double ext_info), for each base pair in each sequence.	
                        for(int i_seq = 0; i_seq < this->sequences->size(); i_seq++)
                        {
                                for(int i = 1; i <= this->sequences->at(i_seq)->numofbases; i++)
                                {
                                        for(int j = i+1; j <= this->sequences->at(i_seq)->numofbases; j++)
                                        {
                                                //this->SetExtrinsic(i_seq, i, j, this->extrinsic_info_list->at(i_seq)->x(i, j)); 
						this->folders->at(i_seq)->SetExtrinsic(i, j, this->extrinsic_info_list->at(i_seq)->x(i, j));
                                        } // j loop
                                } // i loop             
                        } // i_Seq loop 			 
		}

		/*
		The list of threads that are refolding the respective sequence:
		refolding_threads[i_seq] is the pointer to the thread that refolds sequence at i_seq.
		*/
		this->refolding_threads = (t_turbofold_thread**)malloc(sizeof(t_turbofold_thread*) * (this->sequences->size() + 1));
		for(int _i_seq = 0; _i_seq < this->sequences->size(); _i_seq++)
		{
			refolding_threads[_i_seq] = NULL; // NULL refers to a thread not being assigned to this sequence, yet.
		} // _i_seq loop for initing folding thread pointers.

#ifdef COMPILE_SMP
		// If there are more processors available than are threads, start (K-1) threads and Kth sequence will be folded by the main thread.
		int n_threads_to_init = (this->sequences->size() > this->n_parallel_pfunctions)?(this->n_parallel_pfunctions-1):(this->sequences->size()-1);
		for(int _i_seq = 0; _i_seq < n_threads_to_init; _i_seq++)
		{
			this->start_next_refolding_thread(false, i_iter); // Start children threads for computing.

			if(this->err_code != 0)
			{
				return(this->err_code);
			}
		}
#endif

		// Infinite loop until all sequences are refolded.
		while(1)
		{
			int i_seq = this->start_next_refolding_thread(true, i_iter); // This is the main thread.

                        if(this->err_code != 0)
                        {
                                return(this->err_code);
                        }

			// If there is no sequence to refold for main thread, exit from the infinite loop.
			if(i_seq == this->sequences->size())
			{
				break;
			}

			//printf("Refolding %s with main thread.\n", sequences->at(i_seq)->ctlabel);

                        if(i_iter == n_iterations &&
                                this->saves != NULL)
                        {
                                this->folders->at(i_seq)->PartitionFunction(this->saves->at(i_seq));
                        }
                        else
                        {
                                int ret = this->folders->at(i_seq)->PartitionFunction();
				if(ret != 0)
				{
					//printf("PartitionFunction() returned %d\n", ret);
					this->err_code = RNALIB_PARTITIONFUNCTION_ERROR;
					return(this->err_code);
				}
                        }

			//t_matrix* pp_matrix = new t_matrix(new_energy_loop->coinc_posteriors->v_coinc_info);
			//char cur_pp_matrix_fp[500];
			//sprintf(cur_pp_matrix_fp, "%s_pp_%d.txt", sequences->at(i_seq)->ctlabel, i_iter);
			//pp_matrix->dump_matrix(cur_pp_matrix_fp);
			//delete(pp_matrix);

			// Insert this to energy_loops, remove the element before that, which was NULL.
			//printf("Refolded %s\n", sequences->at(i_seq)->ctlabel);
		} // i_seq loop.

#ifdef COMPILE_SMP
		// Wait for other threads which may still be refolding sequences.
		for(int _i_seq = 0; _i_seq < this->sequences->size(); _i_seq++)
		{
			// A sanity check.
			if(this->refolding_threads[_i_seq] == NULL)
			{
				//printf("Sequence %s is not folded, yet, while main thread is not refolding anything.\n", this->sequences->at(_i_seq)->ctlabel);
				this->err_code = THREAD_SCHEDULE_ERROR;
				return(this->err_code);
			}

			// Wait for this thread: Do not wait for threads that were actually main thread.
			if(this->refolding_threads[_i_seq]->computation_thread != NULL)
			{
				// Wait for this thread.
				//printf("Waiting for thread %d to finish.\n", this->refolding_threads[_i_seq]);
				this->refolding_threads[_i_seq]->computation_thread->wait_thread();

				// Check the local signaling, based on the class member thread_state.
				if(this->refolding_threads[_i_seq]->computation_thread->thread_state != THREAD_TERMINAL)
				{
					//printf("Thread %d not finished @ %s(%d)\n", this->refolding_threads[_i_seq]->computation_thread,
					//	__FILE__, __LINE__);
					this->err_code = UNFINISHED_THREAD_ERROR;
					return(this->err_code);						
				}
			}
		} // cur_thread_i loop.
#endif // COMPILE_SMP

/*
		for(int _i_seq = 0; _i_seq < this->sequences->size(); _i_seq++)
		{
			printf("%s: Thread %d\n", this->sequences->at(_i_seq)->ctlabel, this->refolding_threads[_i_seq]);
		} // _i_seq loop.
*/
		// Free the threads now.
		for(int _i_seq = 0; _i_seq < this->sequences->size(); _i_seq++)
		{
			// Delete the threads, make sure the destructors of threads do not free the energy loops.
			delete(this->refolding_threads[_i_seq]);

			// Set all the thread pointers to NULL.
			this->refolding_threads[_i_seq] = NULL;
		} // cur_thread_i loop.

		// The list of refolding thread pointers is reallocated above, should be freed here.
		free(this->refolding_threads);
	} // i_iter loop.

	//Add a coarse update of progress:
	if (progress!=NULL) {

		progress->update(100);
	}

	return(0);
}

int TurboFold::generate_extrinsic_information(int i_iter, const double gamma)
{
	// Reset the extrinsic information.
	for(unsigned int i_seq = 0; i_seq < this->sequences->size(); i_seq++)
	{
		for(int i = 1; i <= this->sequences->at(i_seq)->numofbases; i++)
		{
			for(int j = i+1; j <= this->sequences->at(i_seq)->numofbases; j++)
			{
				this->extrinsic_info_list->at(i_seq)->x(i,j) = 0.0f;
			} // j loop
		} // i loop 
	} // i_seq loop

       // Main computation loop.
        for(unsigned int i_seq1 = 0; i_seq1 < this->sequences->size(); i_seq1++)
        {
                //t_energy_loops* seq1_energy_loops = energy_loops->at(i_seq1);
		RNA* RNA1 = this->folders->at(i_seq1);

                for(unsigned int i_seq2 = i_seq1+1; i_seq2 < this->sequences->size(); i_seq2++)
                {       
                        //t_energy_loops* seq2_energy_loops = energy_loops->at(i_seq2);
			RNA* RNA2 = this->folders->at(i_seq2);
                 
                        for(int i = 1; i <= this->sequences->at(i_seq1)->numofbases; i++)
                        {
                                for(int j = i+1; j <= this->sequences->at(i_seq1)->numofbases; j++)
                                {
                                        int low_k = max(1, this->aln_env_results[i_seq1][i_seq2]->low_limits[i]);
                                        int high_k = this->aln_env_results[i_seq1][i_seq2]->high_limits[i];

                                        for(int k = low_k; 
                                                k <= high_k;
                                                k++)
                                        {
                                                int low_l = max(k+1, this->aln_env_results[i_seq1][i_seq2]->low_limits[j]);
                                                int high_l = this->aln_env_results[i_seq1][i_seq2]->high_limits[j];

                                                for(int l = low_l; l <= high_l; l++)
                                                {
                                                        double seq1_seq2_mapping_probability = this->aln_mapping_probs[i_seq1][i_seq2][i][k] * this->aln_mapping_probs[i_seq1][i_seq2][j][l];
                                                        double seq1_seq2_seq_similarity_weight = (1.0f - this->similarities[i_seq1][i_seq2]);

							//printf("%d, %d: %lf\n", k, l, this->extrinsic_info_list->at(i_seq2)->x(k,l));

                                                        this->extrinsic_info_list->at(i_seq2)->x(k,l) += seq1_seq2_seq_similarity_weight * seq1_seq2_mapping_probability * RNA1->GetPairProbability(i,j);

							this->extrinsic_info_list->at(i_seq1)->x(i,j) += seq1_seq2_seq_similarity_weight * seq1_seq2_mapping_probability * RNA2->GetPairProbability(k,l);

                                                        if(RNA1->GetErrorCode() != 0)
                                                        {
                                                                printf("Problem getting pairing probability for (%d, %d) in sequence %d\n", i, j, i_seq1);
								this->err_code = RNALIB_GETPAIRPROBABILITY_ERROR;
								return(this->err_code);
                                                        }

                                                        if(RNA2->GetErrorCode() != 0)
                                                        {
                                                                printf("Problem getting pairing probability for (%d, %d) in sequence %d\n", k, l, i_seq2);
                                                                this->err_code = RNALIB_GETPAIRPROBABILITY_ERROR;
                                                                return(this->err_code);
                                                        }

                                                        double seq2_seq1_mapping_probability = seq1_seq2_mapping_probability;
                                                        double seq2_seq1_seq_similarity_weight = (1.0f - this->similarities[i_seq1][i_seq2]);
                                                } // l loop
                                        } // k loop
                                } // j loop
                        } // i loop.
                } // i_seq2 loop.
        } // i_seq1 loop.


        // Post process. For normalization and powerizing.
        for(unsigned int i_seq = 0; i_seq < this->sequences->size(); i_seq++)
        {
                //char cur_ext_info_fp[1000];
                //sprintf(cur_ext_info_fp, "ext_info_%d_iter_%d.txt", i_seq, i_iter);
                //this->extrinsic_info_list->at(i_seq)->dump_matrix(cur_ext_info_fp);

                // Now normalize everything. The case where maximum is 0 is handled inside normalize_by_max function.
                this->extrinsic_info_list->at(i_seq)->normalize_by_max();
                this->extrinsic_info_list->at(i_seq)->powerize_each_element(gamma);

        } // i_seq loop.

	return(0);
}

int TurboFold::initialize_alignment_information()
{
        // Allocate similarities.
        this->similarities = (double**)malloc(sizeof(double*) * (this->sequences->size() + 2));
        for(int i_seq1 = 0; i_seq1 < this->sequences->size(); i_seq1++)
        {
                this->similarities[i_seq1] = (double*)malloc(sizeof(double) * (this->sequences->size() + 2));

                for(int i_seq2 = 0; i_seq2 < this->sequences->size(); i_seq2++)
                {
                        this->similarities[i_seq1][i_seq2] = 0.0f;
                } // i_seq2 loop
        } // i_seq1 loop

        // Compute the alignment probabilities.
        //this->pp_results = (t_pp_result***)malloc(sizeof(t_pp_result**) * (sequences->size() + 1));
        this->aln_env_results = (t_aln_env_result***)malloc(sizeof(t_aln_env_result**) * (sequences->size() + 1));
        this->aln_mapping_probs = (double****)malloc(sizeof(double***) * (sequences->size() + 1));

        double alignment_info_mem_size = 0.0f;

        for(unsigned int i_seq1 = 0; i_seq1 < sequences->size(); i_seq1++)
        {
                aln_env_results[i_seq1] = (t_aln_env_result**)malloc((sequences->size() + 2) * sizeof(t_aln_env_result*));

                this->aln_mapping_probs[i_seq1] = (double***)malloc((sequences->size() + 2) * sizeof(double**));
                alignment_info_mem_size += ((sequences->size() + 2) * sizeof(double**));

                // Move the pointers.
                this->aln_env_results[i_seq1] -= i_seq1;
                this->aln_mapping_probs[i_seq1] -= i_seq1;

                //pp_results[i_seq1] = (t_pp_result**)malloc((sequences->size() + 2) * sizeof(t_pp_result*));
                for(unsigned int i_seq2 = i_seq1+1; i_seq2 < sequences->size(); i_seq2++)
                {
                        if(i_seq1 != i_seq2)
                        {
                                t_phmm_aln* phmm_aln = new t_phmm_aln(sequences->at(i_seq1), sequences->at(i_seq2));

                                t_pp_result* cur_pp_results = phmm_aln->compute_posterior_probs();
                                //pp_results[i_seq1][i_seq2] = cur_pp_results;

                                t_aln_env_result* cur_aln_env_result = phmm_aln->compute_alignment_envelope(PROB_ALN_ENV, cur_pp_results, cur_pp_results->fam_threshold, 7);

                                // Copy the ML similarity.
                                this->similarities[i_seq1][i_seq2] = cur_pp_results->ml_similarity;

				//printf("%d-%d similarity: %lf\n", i_seq1, i_seq2, cur_pp_results->ml_similarity);

                                aln_env_results[i_seq1][i_seq2] = cur_aln_env_result;

                                // Allocate the alignment mapping probability matrix for current sequences. Note that these are 1 based.
                                this->aln_mapping_probs[i_seq1][i_seq2] = (double**)malloc(sizeof(double*) * (this->sequences->at(i_seq1)->numofbases + 2));  
                                alignment_info_mem_size += (sizeof(double*) * (this->sequences->at(i_seq1)->numofbases + 2));

                                for(int i = 1; i <= this->sequences->at(i_seq1)->numofbases; i++)
                                {
                                        int min_k = cur_aln_env_result->low_limits[i];
                                        int max_k = cur_aln_env_result->high_limits[i];
                                        this->aln_mapping_probs[i_seq1][i_seq2][i] = (double*)malloc(sizeof(double) * (max_k - min_k + 2));
                                        alignment_info_mem_size += (sizeof(double) * (max_k - min_k + 2));

                                        this->aln_mapping_probs[i_seq1][i_seq2][i] -= min_k;
                                        for(int k = min_k; k <= max_k; k++)
                                        {
                                                // Mapping probability is coincidence probability of the nucleotides.
                                                double aln_prob = exp(cur_pp_results->aln_probs[i][k]);
                                                double ins1_prob = exp(cur_pp_results->ins1_probs[i][k]);
                                                double ins2_prob = exp(cur_pp_results->ins2_probs[i][k]);

                                                double mapping_probability = (aln_prob + ins1_prob + ins2_prob);

                                                this->aln_mapping_probs[i_seq1][i_seq2][i][k] = mapping_probability;
                                        } // k loop
                                } // i loop

/*
				char cur_map_fp[1000];
				sprintf(cur_map_fp, "aln_mapping_probs_%d_%d.txt", i_seq1, i_seq2);
				FILE* f_aln_map_prob = fopen(cur_map_fp, "w");
				for(int i = 1; i <= this->sequences->at(i_seq1)->numofbases; i++)
				{
					for(int k = 1; k <= this->sequences->at(i_seq2)->numofbases; k++)
					{
						if(k >= cur_aln_env_result->low_limits[i] && k <= cur_aln_env_result->high_limits[i])
						{
							fprintf(f_aln_map_prob, "%lf ", this->aln_mapping_probs[i_seq1][i_seq2][i][k]);
						}
						else
						{
							fprintf(f_aln_map_prob, "0.0 ");
						}
					} // k loop

					fprintf(f_aln_map_prob, "\n");
				} // i loop
				fclose(f_aln_map_prob);
*/

                                // Free current pp result.
                                phmm_aln->free_pp_result(cur_pp_results);
                                delete(phmm_aln);
                        }
                        else
                        {
                                //pp_results[i_seq1][i_seq2] = NULL;
                                aln_env_results[i_seq1][i_seq2] = NULL;
                        }
                } // i_seq2 loop.
        } // i_seq1 loop.

	return(0);
}

int TurboFold::ProbKnot(const int i_seq, const int n_iterations, const int min_helix_length)
{
        if(i_seq > this->GetNumberSequences())
        {
                this->err_code = ISEQ_ARGUMENT_OVERFLOW_ERROR;

        }
        else
        {

		int ret = this->folders->at(i_seq-1)->ProbKnot(n_iterations, min_helix_length);
	
		if(ret != 0)
		{
			this->err_code = RNALIB_PROBKNOT_ERROR;
		}
	        else
	        {
	                this->err_code = 0;
	        }
	}

        return(this->err_code);
}

int TurboFold::PredictProbablePairs(const int i_seq, const float probability)
{
        if(i_seq > this->GetNumberSequences())
        {
                this->err_code = ISEQ_ARGUMENT_OVERFLOW_ERROR;

        }
        else
        {
	        int ret = this->folders->at(i_seq-1)->PredictProbablePairs(probability);
	
	        if(ret != 0)
	        {
	                this->err_code = RNALIB_THRESHOLDING_ERROR;
	        }
	        else
	        {
			this->err_code = 0;
	        }
	}
	return(this->err_code);

}

int TurboFold::MaximizeExpectedAccuracy(const int i_seq, const double maxPercent, const int maxStructures, const int window, const double gamma)
{
        if(i_seq > this->GetNumberSequences())
        {
		this->err_code = ISEQ_ARGUMENT_OVERFLOW_ERROR;
		
        }
	else
	{
	        int ret = this->folders->at(i_seq-1)->MaximizeExpectedAccuracy(maxPercent, maxStructures, window, gamma);
	
	        if(ret != 0)
	        {
	                this->err_code = RNALIB_MEA_ERROR;
	        }
	        else
	        {
	                this->err_code = 0;
	        }
	}
        return(this->err_code);

}

int TurboFold::GetPair(const int i_seq, const int i, const int structurenumber)
{
        if(i_seq > this->GetNumberSequences())
        {
                this->err_code = ISEQ_ARGUMENT_OVERFLOW_ERROR;
        }
        else
        {
		int ret = this->folders->at(i_seq-1)->GetPair(i, structurenumber);

		if(ret != 0)
		{
			this->err_code = RNALIB_GETPAIR_ERROR;
		}
		else
		{
			this->err_code = 0;
		}
	}

	return(this->err_code);
}

int TurboFold::WriteCt(const int i_seq, const char fp[])
{
        if(i_seq > this->GetNumberSequences())
        {
                this->err_code = ISEQ_ARGUMENT_OVERFLOW_ERROR;
        }
        else
        {
                int ret = this->folders->at(i_seq-1)->WriteCt(fp);

                if(ret != 0)
                {
                        this->err_code = RNALIB_WRITECT_ERROR;
                }
                else
                {
                        this->err_code = 0;
                }
        }       

        return(this->err_code);

}

double TurboFold::GetPairProbability(const int i_seq, const int i, const int j)
{
        if(i_seq > this->GetNumberSequences())
        {
                this->err_code = ISEQ_ARGUMENT_OVERFLOW_ERROR;
		return(0.0f);
        }
        else
        {
                double pair_prob = this->folders->at(i_seq-1)->GetPairProbability(i,j);
		int ret = this->folders->at(i_seq-1)->GetErrorCode();

                if(ret != 0)
                {
                        this->err_code = RNALIB_GETPAIRPROBABILITY_ERROR;
			return(0.0f);
                }
                else
                {
                        this->err_code = 0;
			return(pair_prob);
                }
        }
}

int TurboFold::ReadSHAPE(const int i_seq, const char fp[], const double par1, const double par2)
{
        if(i_seq > this->GetNumberSequences())
        {
                this->err_code = ISEQ_ARGUMENT_OVERFLOW_ERROR;
        }
        else
        {
				//Note that 0.0 and 0.0 are hard-wired because these are single-stranded terms that do not apply using current practices.
                int ret = this->folders->at(i_seq-1)->ReadSHAPE(fp, par1, par2, 0.0, 0.0);

                if(ret != 0)
                {
                        this->err_code = RNALIB_READSHAPE_ERROR;
                }
                else
                {
                        this->err_code = 0;
                }
        }

        return(this->err_code);

}

int TurboFold::GetNumberSequences()
{
	return(this->sequences->size());
}

int TurboFold::SetTemperature(const double temp)
{
	this->err_code = 0;

	for(int i_seq = 1; i_seq <= this->sequences->size(); i_seq++)
	{
                int ret = this->folders->at(i_seq-1)->SetTemperature(temp);

                if(ret != 0)
                {
                        this->err_code = RNALIB_SETTEMPERATURE_ERROR;
                }
                else
                {
                        this->err_code = 0;
                }
        }

        return(this->err_code);
}

/*
Starts a new refolding thread. A slot is always dedicated to the main thread.
*/
int TurboFold::start_next_refolding_thread(bool is_main_thread, int _i_iter)
{
#ifdef COMPILE_SMP
	// MUST LOCK THE COMPUTATIONS AT THIS POINT!
	this->turbofold_threads_mutex->lock_mutex();
#endif // COMPILE_SMP

	if(is_main_thread)
	{
		// Check if there is a sequence to start folding.
		for(int _i_seq = 0; _i_seq < this->sequences->size(); _i_seq++)
		{
			// Check if this sequence is being refolded.
			if(this->refolding_threads[_i_seq] == NULL)
			{
				// Assign sequence at _i_seq to main thread.
				t_turbofold_thread* cur_refolding_thread = new t_turbofold_thread(this, _i_seq, _i_iter);
				this->refolding_threads[_i_seq] = cur_refolding_thread;

				//printf("Assigned %s to main thread.\n", this->sequences->at(_i_seq)->ctlabel);

#ifdef COMPILE_SMP
				// RELEASE THE COMPUTATION LOCK AND RETURN.
				this->turbofold_threads_mutex->release_mutex();
#endif // COMPILE_SMP
				return(_i_seq);
			}
		} // _i_seq loop.

		// If there are no more sequences to refold, return an indicator, which is K (Number of sequences) in this case.
#ifdef COMPILE_SMP 
		this->turbofold_threads_mutex->release_mutex();
#endif // COMPILE_SMP
		return(this->sequences->size());
	}
	else
	{
#ifndef COMPILE_SMP
		return(this->sequences->size());
#endif		
	}

#ifdef COMPILE_SMP
	// Following is for children threads.
	// Check the number of threads alive.
	int n_alive_children_threads = 0;

	// Count the number of children threads that are alive, note that main thread is always running. 
	for(int _i_seq = 0; _i_seq < this->sequences->size(); _i_seq++)
	{
		// Check if this sequence is being refolded.
		if(this->refolding_threads[_i_seq] != NULL)
		{
			// Is there a child thread that is refolding this sequence? Following checks if a child is folding the sequence.
			if(this->refolding_threads[_i_seq]->computation_thread != NULL &&
				this->refolding_threads[_i_seq]->computation_thread->thread_state != THREAD_TERMINAL)
			{
				n_alive_children_threads++;
			}
		} 
	} // _i_seq loop.

	if(n_alive_children_threads < this->n_parallel_pfunctions - 1)
	{
		// Check if there is a sequence to start folding.
		for(int _i_seq = 0; _i_seq < this->sequences->size(); _i_seq++)
		{
			// Check if this sequence is being refolded.
			if(this->refolding_threads[_i_seq] == NULL)
			{
				// Can assign this _i_seq to a new thread.
				t_turbofold_thread* cur_refolder_thread = new t_turbofold_thread(this, _i_seq, _i_iter);
				if(cur_refolder_thread->run_refolding_thread() != 0)
				{
					this->err_code = RUN_REFOLDING_THREAD_ERROR;
					return(err_code);
				}

				this->refolding_threads[_i_seq] = cur_refolder_thread;

				//printf("Assigned %s to thread %d.\n", this->sequences->at(_i_seq)->ctlabel, cur_refolder_thread);

				// RELEASE THE COMPUTATION LOCK AND RETURN.
				this->turbofold_threads_mutex->release_mutex();
				return(_i_seq);
			}
		} // _i_seq loop.

		// At this point, there are no more sequences left to refold, release computation lock and return.
		// RELEASE THE COMPUTATION LOCK AND RETURN.
		this->turbofold_threads_mutex->release_mutex();
		return(this->sequences->size());
	}
	else
	{
		// Cannot start a new thread at this point, must wait until a child thread finishes.

		// RELEASE THE COMPUTATION LOCK AND RETURN.
		this->turbofold_threads_mutex->release_mutex();
		return(this->sequences->size());
	}
#endif // COMPILE_SMP
}

int TurboFold::GetErrorCode()
{
	return(this->err_code);
}

char* TurboFold::GetErrorMessage(const int err_code)
{
	return(err_strings[err_code]);
}

string TurboFold::GetErrorString(const int err_code)
{
	return(string(err_strings[err_code]));
}

//Provide a TProgressDialog for following calculation progress.
//A TProgressDialog class has a public function void update(int percent) that indicates the progress of a long calculation.
void TurboFold::SetProgress(TProgressDialog& Progress) {


	progress = &Progress;

	return;

}


//Provide a means to stop using a TProgressDialog.
//StopProgress tells the RNA class to no longer follow progress.  This should be called if the TProgressDialog is deleted, so that this class does not make reference to it.
void TurboFold::StopProgress() {

	progress=NULL;
	return;

}
