
//Class RNA -- Wrapper for RNAstructure for use in object-oriented applications

#include "RNA.h"

#include "../src/structure.h"
#include "../src/arrayclass.h"
#include "../src/dotarray.h"
#include "../src/forceclass.h"
#include "../src/rna_library.h"
#include "../src/stackclass.h"
#include "../src/stackstruct.h"



#include "../src/algorithm.h"
#include "../src/outputconstraints.h"
#include "../src/pfunction.h"
#include "../src/boltzmann.h"
#include "../src/alltrace.h"
#include "../src/stochastic.h"
#include "../src/MaxExpect.h"
#include "../src/probknot.h"

#ifdef _CUDA_CALC_
#include "../partition-smp/param.h"
#include "../partition-smp/prna.h"
#include "../partition-smp/util.h"
#include "../partition-smp/base.h"
#endif//_CUDA_CALC

#include <iostream>

const float epsilon = 1e-6; // a small number for a tolerance in comparing floats


//constructor where user provides a string with the sequence
RNA::RNA(const char sequence[], const bool IsRNA):Thermodynamics(IsRNA) {
	int i;



	//allocate ct
	ct = new structure();

	//Specify the sequence length based on the string length
	
	//ct->numofbases = (short) strlen(sequence);
	ct->allocate((int) strlen(sequence));//allocate the space required in the arrays for the sequence


	//Now store the sequence information
	for (i=1;i<=ct->GetSequenceLength();i++) {
		if (sequence[i-1]=='A'||sequence[i-1]=='a') ct->numseq[i]=1;
		else if (sequence[i-1]=='C'||sequence[i-1]=='c') ct->numseq[i]=2;
		else if (sequence[i-1]=='G'||sequence[i-1]=='g') ct->numseq[i]=3;
		else if (sequence[i-1]=='U'||sequence[i-1]=='u'||sequence[i-1]=='T'||sequence[i-1]=='t') ct->numseq[i]=4;
		else ct->numseq[i]=0;

		ct->nucs[i] = sequence[i-1];
		ct->hnumber[i]=i;
	}



	//These should not be located here: (DHM 5/3/2011)
	// Experimental bonuses array
   	//ct->EX = new double *[ct->GetSequenceLength()+1];
   	//for(i=0;i<ct->GetSequenceLength()+1;i++) {
	//  ct->EX[i] = new double[ct->GetSequenceLength()+1];
   	//}
	//for (i=0;i<ct->GetSequenceLength()+1;i++){
	//  for (j=0;j<ct->GetSequenceLength()+1;j++){
	 //   ct->EX[i][j] = 0.0;
	 // }
     //   }


	//Indicate that the partition function calculation has not been performed.
	partitionfunctionallocated = false;

	//Indicate that the energy data is not read.
	energyallocated = false;

	//Drawing coordinates have not been determined.
	drawallocated = false;

	//set error status to zero
	ErrorCode=0;

	//Do not report progress by default:
	progress=NULL;



}

//constructor where user provides a string with a filename
//	The flag type indicates the file type: type = 1 => ct file, type = 2 => .seq file, type = 3 => .pfs file.
//	The fconstructor saves an error code in ErrorCode:
//	0 = no error, 1 = file not found
//	2 = error opening file.
//  GetErrorCode provides public access to this errorcode.
RNA::RNA(const char filename[], const int type, const bool IsRNA ):Thermodynamics(IsRNA) {


	//allocate ct
	ct = new structure();

	//Indicate that the partition function calculation has not been performed.
	partitionfunctionallocated = false;

	//Indicate that the energy data is not (yet) read.
	energyallocated = false;

	//Drawing coordinates have not been determined.
	drawallocated = false;

	//Do not report progress by default:
	progress=NULL;

	ErrorCode = FileReader(filename, type);


	return;


}

//Default constructor.
RNA::RNA(const bool IsRNA):Thermodynamics(IsRNA) {





	//allocate the underlying structure class and nothing more.
	//User is then required to propogate the sequence and structural information manually.
	ct = new structure();

	//set the number of nucleotides to zero because no sequence has been read
	//ct->numofbases = 0;

	//Indicate that the partition function calculation has not been performed.
	partitionfunctionallocated = false;

	//Indicate that the energy data is not (yet) read.
	energyallocated = false;

	//Drawing coordinates have not been determined.
	drawallocated = false;

	//set error status to zero
	ErrorCode=0;

	//Do not report progress by default:
	progress=NULL;



}


//Return the value of ErrorCode
int RNA::GetErrorCode() {
	return ErrorCode;
}

//Return a c string that describes errors from GetErrorCode and other errors.
char* RNA::GetErrorMessage(const int error) {

	if (error==0) return "No Error.\n";
	else if (error==1) return "Input file not found.\n";
	else if (error==2) return "Error opening file.\n";
	else if (error==3) return "Structure number out of range.\n";
	else if (error==4) return "Nucleotide number out of range.\n";
	else if (error==5) return "Error reading thermodynamic parameters.\nPlease set environment variable $DATAPATH to the location of the thermodynamic parameters.\n";
	else if (error==6) return "This would form a pseudoknot and is not allowed.\n";
	else if (error==7) return "This pair is non-canonical and is therefore not allowed.\n";
	else if (error==8) return "Too many restraints specified.\n";
	else if (error==9) return "This nucleotide already under a conflicting constraint.\n";
	else if (error==10) return "There are no structures to write to file.\n";
	else if (error==11) return "Nucleotide is not a U.\n";
	else if (error==12) return "Maximum pairing distance is too short.\n";
	else if (error==13) return "Error reading constraint file.\n";
	else if (error==14) return "A traceback error occurred.\n";
	else if (error==15) return "No partition function data is available.\n";
	else if (error==16) return "Wrong save file version used or file format not recognized.\n";
	else if (error==17) return "This function cannot be performed unless a save file (.sav) was correctly loaded by the RNA constructor.\n";
	else if (error==18) return "This threshold is too low to generate valide secondary structures.\n";
	else if (error==19) return "The structure coordinates have not been determined, use DetermineDrawingCoordinates() to calculate the coordinates.\n";
	else if (error==20) return "No sequence has been read.\n";
	else if (error==21) return "Probabilities summed to greater than 1 in stochastic traceback.\n";
	else if (error==22) return "Programming error.  Incorrect file type passed to constructor.\n";
	else if (error==23) return "There are no structures present.\n";
	else if (error==24) return "Too few iterations.  There must be at least one iteration.\n";
	else if (error==25) return "Index is not a multiple of 10.\n";
	else if (error==26) return "k, the equilibrium constant, needs to be greater than or equal to 0.\n";
	else return "Unknown Error\n";


}

//Return a string that describes errors from GetErrorCode and other errors.
//This uses GetErrorMessage to acrually get the errors.
std::string RNA::GetErrorMessageString(const int error) {
	std::string temp;

	temp = GetErrorMessage(error);

	return temp;


}


//User specifies a base pair between i and j in structure # structurenumber.
int RNA::SpecifyPair(const int i, const int j, const int structurenumber) {
	


	//start with error checking:
	if (i<0||i>ct->GetSequenceLength()||j<0||j>ct->GetSequenceLength()) return 4;
	else if (structurenumber<1) return 3;

	//also keep track of the maximum number of structures for which pairs have been specified
	if (structurenumber>ct->GetNumberofStructures()) {
		

		//Add one structure for each position between those available and the one being specified:
		for (int index=ct->GetNumberofStructures()+1;index<=structurenumber;++index) ct->AddStructure();
			

	}

	//now register the pair:
	ct->SetPair(i,j,structurenumber);


	return 0;

}
//Break a pair that i is involved in
int RNA::RemoveBasePair(const int i, const int structurenumber) {

	//start with error checking:
	if (i<0||i>ct->GetSequenceLength()) return 4;
	else if (structurenumber<1||structurenumber>ct->GetNumberofStructures()) return 3;

	//Call the function for this in the underlying structure class.
	ct->RemovePair(i, structurenumber);


	//return that there was no error
	return 0;

}

//remove all pairs in structure # structurenumber.
//Also, roll back the number of specified structures if this is the last specified structure.
int RNA::RemovePairs(const int structurenumber) {
	
	//do some error checking
	if (structurenumber>ct->GetNumberofStructures()||structurenumber<1) return 5;


	//decrement the number of structures, if appropriate, i.e. this is the last structure
	if (structurenumber==ct->GetNumberofStructures()) {
		ct->RemoveLastStructure();
		return 0;
	}

	//otherwise, clean the selected structure of pairs:
	ct->CleanStructure(structurenumber);
	return 0;

}


//Calculate and return the folding free energy change for structure number structurenumber.
double RNA::CalculateFreeEnergy(const int structurenumber, const bool UseSimpleMBLoopRules) {
	//Do some simple error checking
	if (structurenumber<1||structurenumber>ct->GetNumberofStructures()) return 0.0;

	if (!energyread) {
		//The thermodynamic data tables have not yet been read
		if (ReadThermodynamic()!=0) {
			ErrorCode = 5;//record an error
			return 0.0;//return 0.0 if a problem occurs

		}
		else ErrorCode=0;//record that there was no error.
	}
	else ErrorCode=0;//Set the error code to zero because no errors were encountered.

	efn2(data,ct,structurenumber,UseSimpleMBLoopRules);

	//conversion factor is set in defines.h.  Free energies are multiplied by this factor internally so that integer math can be used.
	return (((double) ct->GetEnergy(structurenumber)/conversionfactor));


}

//Write the details on the energy caclulation for all structures.
int RNA::WriteThermodynamicDetails(const char filename[], const bool UseSimpleMBLoopRules) {

	if (!energyread) {
		//The thermodynamic data tables have not yet been read
		if (ReadThermodynamic()!=0) return 5;//return non-zero if a problem occurs
	}
	efn2(data,ct,0,UseSimpleMBLoopRules,filename);
	return 0;


}

#ifndef DYNALIGN_II
//Predict the secondary structure by free energy minimization.
//Also generate subooptimal solutions using a heuristic.
int RNA::FoldSingleStrand(const float percent, const int maximumstructures, const int window, const char savefile[], const int maxinternalloopsize, bool mfeonly) {
	char *savefilename;
	int percenti;
	int tracebackstatus;

	//check to make sure that a sequence has been read
	if (ct->GetSequenceLength()==0) return 20;


	if (!energyread) {
		//The thermodynamic data tables have not been read and need to be read now.
		if (ReadThermodynamic()!=0) return 5;//return non-zero if a problem occurs

	}

	//savefile will be passed to the function dynamic for structure prediction.
	//Dynamic expects a null pointer if no file is to be created, so set savefilename to null if savefile is an empty string.
	if (!strcmp(savefile,"")) savefilename = NULL;
	else {
		savefilename=new char[((int) strlen(savefile))+1];
		strcpy(savefilename,savefile);
	}

	//right now, dynamic requires an integer specification of percent energy change.
	//FoldSingleStrand takes this is a float to provide the opportunity to reform this in the future.
	//For now, cast percent as an integer to pass to dynamic.
	percenti = (int) percent;

	//Predict the secondary structures.
	tracebackstatus = dynamic(ct, data, maximumstructures, percenti, window, progress, false, savefilename, maxinternalloopsize,mfeonly);

	//Clean up the memory use.
	delete[] savefilename;

	if(tracebackstatus!=0) return 14;//This indicates a traceback error.
	else return 0;



}

#else
#endif
// Predict the lowest free energy secondary structure and generate all suboptimal structures.
int RNA::GenerateAllSuboptimalStructures(const float percent, const double deltaG) {

	//check to make sure that a sequence has been read
	if (ct->GetSequenceLength()==0) return 20;


	if (!energyread) {
		//The thermodynamic data tables have not been read and need to be read now.
		if (ReadThermodynamic()!=0) return 5;//return non-zero if a problem occurs

	}


	//Call the alltrace function to do the work:
	alltrace(ct,data, ((short) percent), ((short) (deltaG*conversionfactor)),progress,NULL);

	return 0;


}

// Predict the structure with maximum expected accuracy and suboptimal structures.
int RNA::MaximizeExpectedAccuracy(const double maxPercent, const int maxStructures, const int window, const double gamma) {

	//first trap some possible errors
	if (!partitionfunctionallocated) {
		//There is no partition function data available.
		return 15;
	}


	//Past error trapping
	MaxExpectFill(ct, v, w5, pfdata, lfce, mod, fce, maxPercent, maxStructures, window, gamma, progress);

	return 0;//no error return functionality right now

}

// This function predicts structures composed of probable base pairs.
int RNA::PredictProbablePairs(const float probability) {
	int i,j,count;
	char thresh[8];
	string label;//A string for making ct file labels

	//first trap some possible errors
	if (probability > epsilon && probability < 0.500-epsilon) {
		//The threshold is too low to be valie and not low enough that it will be considered zero, a the default
		return 18;
	}

	if (!partitionfunctionallocated) {
		//There is no partition function data available.
		return 15;
	}


	//Past error trapping

	
	
	

	
	if (probability>epsilon) {
		//The user specified a threshold, so use that and generate one structure
		

		//Get one clean structure
		if (ct->GetNumberofStructures()>0) {
			ct->CleanStructure(1);
			for (i=ct->GetNumberofStructures();i>1;--i) {
				ct->RemoveLastStructure();
			}
		}
		else {

			ct->AddStructure();
		}



		for (i=1;i<ct->GetSequenceLength();i++) {
			for (j=i+1;j<=ct->GetSequenceLength();j++) {

				if (calculateprobability(i,j,v,w5,ct,pfdata,lfce,mod,pfdata->scaling,fce) > probability) {
					//This pair exceeded the threshold, so add it to the list
					ct->SetPair(i,j);
					

				}

			}
		}

		//put the threshold in the ctlabel, which will appear in the header of a ct file
		sprintf(thresh,"%1.5f",probability);
		
		
		//Insert a label at the beging of ctlabel to provide the threshold.  This will be written in an output ct file.
		
		label = " >";
		label+=thresh;
		label+=" pairing probability; ";
		label+=ct->GetCtLabel(1);

		ct->SetCtLabel(label,1);
		
	}
	else {
		//The default threshold was specified, so create 8 structures, with thresholds of >=0.99, >=0.97, >=0.95, >=0.90, >=0.80, >=0.70, >=0.60, >0.50.

		//Set up 8 clean structures:
		
		//Get 8 clean structures
		if (ct->GetNumberofStructures()>8) {
			for (i=ct->GetNumberofStructures();i>8;--i) {
				ct->RemoveLastStructure();
			}
			for (i=1;i<=8;++i) ct->CleanStructure(i);
		}
		else {
			for (i=1;i<=ct->GetNumberofStructures();++i) ct->CleanStructure(i);
			for (i=ct->GetNumberofStructures();i<8;++i) ct->AddStructure();
		}
		

		//loop over the structures and thresholds
		for (count=1;count<=8;count++) {

			for (i=1;i<ct->GetSequenceLength();i++) {
				for (j=i+1;j<=ct->GetSequenceLength();j++) {

					if (count==1) {
						if (calculateprobability(i,j,v,w5,ct,pfdata,lfce,mod,pfdata->scaling,fce)>=.99) {
							
							//set this pair because it meets the threshold
							ct->SetPair(i,j,count);
							
						}
					}
					else if (count==2) {
						if (calculateprobability(i,j,v,w5,ct,pfdata,lfce,mod,pfdata->scaling,fce)>=.97) {
							
							//set this pair because it meets the threshold
							ct->SetPair(i,j,count);
						}
					}
					else if (count==3) {
						if (calculateprobability(i,j,v,w5,ct,pfdata,lfce,mod,pfdata->scaling,fce)>=.95) {
							
							//set this pair because it meets the threshold
							ct->SetPair(i,j,count);
						}
					}
					else if (count==4) {
						if (calculateprobability(i,j,v,w5,ct,pfdata,lfce,mod,pfdata->scaling,fce)>=.90) {
							
							//set this pair because it meets the threshold
							ct->SetPair(i,j,count);
						}
					}
					else if (count==5) {
						if (calculateprobability(i,j,v,w5,ct,pfdata,lfce,mod,pfdata->scaling,fce)>=.80) {
							
							//set this pair because it meets the threshold
							ct->SetPair(i,j,count);
						}
					}
					else if (count==6) {
						if (calculateprobability(i,j,v,w5,ct,pfdata,lfce,mod,pfdata->scaling,fce)>=.70) {
							
							//set this pair because it meets the threshold
							ct->SetPair(i,j,count);
						}
					}
					else if (count==7) {
						if (calculateprobability(i,j,v,w5,ct,pfdata,lfce,mod,pfdata->scaling,fce)>=.60) {
							
							//set this pair because it meets the threshold
							ct->SetPair(i,j,count);
						}
					}
					else if (count==8) {
						if (calculateprobability(i,j,v,w5,ct,pfdata,lfce,mod,pfdata->scaling,fce)>.50) {
							
							//set this pair because it meets the threshold
							ct->SetPair(i,j,count);
						}
					}
				}
			}

		}

		//add labels that would appear in a ct file header
		
		label = " >=97% probable pairs ";
		label+= ct->GetCtLabel(1);
		ct->SetCtLabel(label,2);
		
		label = " >=95% probable pairs ";
		label+= ct->GetCtLabel(1);
		ct->SetCtLabel(label,3);

		label = " >=90% probable pairs ";
		label+= ct->GetCtLabel(1);
		ct->SetCtLabel(label,4);

		label = " >=80% probable pairs ";
		label+= ct->GetCtLabel(1);
		ct->SetCtLabel(label,5);

		label = " >=70% probable pairs ";
		label+= ct->GetCtLabel(1);
		ct->SetCtLabel(label,6);

		label = " >=60% probable pairs ";
		label+= ct->GetCtLabel(1);
		ct->SetCtLabel(label,7);

		label = " >50% probable pairs ";
		label+= ct->GetCtLabel(1);
		ct->SetCtLabel(label,8);

		label = " >=99% probable pairs ";
		label+= ct->GetCtLabel(1);
		ct->SetCtLabel(label,1);

	
			

	}

	return 0;



}

//Calculate the partition function for the current sequence.
int RNA::PartitionFunction(const char savefile[], double temperature) {
	int i,j;
	char *savefilename;
	//check to make sure that a sequence has been read
	if (ct->GetSequenceLength()==0) return 20;
	if (!energyread) {
		//The thermodynamic data tables have not been read and need to be read now.
		if (ReadThermodynamic()!=0) return 5;//return non-zero if a problem occurs

	}

	//savefile will be passed to the function dynamic for structure prediction.
	//Dynamic expects a null pointer if no file is to be created, so set savefilename to null if savefile is an empty string.
	if (!strcmp(savefile,"")) savefilename = NULL;
	else {
		savefilename=new char[((int) strlen(savefile))+1];
		strcpy(savefilename,savefile);
	}

	if (partitionfunctionallocated) {
		delete v;
		delete w;
		delete wmb;
		delete wl;
		delete wmbl;
		delete wcoax;
		delete fce;
		delete[] lfce;
		delete[] mod;
		delete[] w3;
		delete[] w5;
		delete pfdata;
	}
	//Allocate the memory needed (only if this is the first call to pfunction):
	//indicate that the memory has been allocated so that the destructor will delete it.
	partitionfunctionallocated = true;

	//allocate space for the v and w arrays:
	w = new pfunctionclass(ct->GetSequenceLength());
	v = new pfunctionclass(ct->GetSequenceLength());
	wmb = new pfunctionclass(ct->GetSequenceLength());
	wl = new pfunctionclass(ct->GetSequenceLength());
	wmbl = new pfunctionclass(ct->GetSequenceLength());
	wcoax = new pfunctionclass(ct->GetSequenceLength());
	fce = new forceclass(ct->GetSequenceLength());
	lfce = new bool [2*ct->GetSequenceLength()+1];
	mod = new bool [2*ct->GetSequenceLength()+1];

	for (i=0;i<=2*ct->GetSequenceLength();i++) {
		lfce[i] = false;
		mod[i] = false;
	}

	for (i=0;i<ct->GetNumberofModified();i++) {

		if (ct->GetModified(i)!=1&&ct->GetModified(i)!=ct->GetSequenceLength()) {
			mod[ct->GetModified(i)]=true;
			mod[ct->GetModified(i)+ct->GetSequenceLength()]=true;
		}
	}

	w5 = new PFPRECISION [ct->GetSequenceLength()+1];
	w3 = new PFPRECISION [ct->GetSequenceLength()+2];

	if (ct->intermolecular) {
		//take advantage of templating to prevent intramolecular base pairs

		ct->allocatetem();//This allocates a bool array that indicates what pairs are allowed.  It is initilaized to true, i.e. all nucs can pair with each other.
		for (i=1;i<ct->inter[0];i++) {
			for (j=i+1;j<=ct->inter[2];j++) {

				//Set intermolecular pairs to false, i.e. not allowed.  Note the indexing with the high index and then the low index number.
				ct->tem[j][i]=false;

			}
		}
		for (i=ct->inter[2]+1;i<ct->GetSequenceLength();i++) {
			for (j=i+1;j<=ct->GetSequenceLength();j++) {

				//Another set of intermolecular pairs to forbid.
				ct->tem[j][i]=false;

			}
		}
	}

	//Initialize the partition function datatable:
		//Ignore the setting of parameter temperature if it is less than zero.
		//Generally, this parameter should be left at the default.
	if (temperature<0) pfdata=new pfdatatable(data,scalingdefinition,temp);
	else pfdata=new pfdatatable(data,scalingdefinition,temperature);

	//This code converts the SHAPE array of data to equilibrium constants.  This is
	//needed for the partition function.  NOTE, however, there is no going back so a structure
	//that has been used for partition functions cannot then be used to predict a structure
	//by free energy minimization. This is a compromise for efficiency, but clearly something
	//that could cause a problem for programmers.
	if (ct->shaped) {
	  for (i=1;i<=2*ct->GetSequenceLength();i++) ct->SHAPE[i]=boltzman( ct->SHAPE[i], pfdata->temp);
	}

	if (ct->experimentalPairBonusExists){
	  //symmetrize as well...
	  for (i = 1; i <= 2*ct->GetSequenceLength(); i++){
	    for(j = i; j <= 2*ct->GetSequenceLength(); j++){
	      double avg = (double)( 0.5*(ct->EX[i][j] + ct->EX[j][i]) );
	      ct->EX[i][j] = boltzman( avg, pfdata->temp);
	      ct->EX[j][i] = ct->EX[i][j];
	    }
	  }
	}


	//The next section handles the case where base pairs are not
			//not allowed to form between nucs more distant
			//than ct->GetPairingDistanceLimit()
	if (ct->DistanceLimited()) {
		//This allocates a bool array that indicates what pairs are allowed.  It is initilaized to true, i.e. all nucs can pair with each other.
		if (!ct->templated) ct->allocatetem();
		for (j=minloop+2;j<=ct->GetSequenceLength();j++) {
			for (i=1;i<j;i++) {
				//Set distant pairs to false, i.e. not allowed.  Note the indexing with the high index and then the low index number.
				if (j-i>=ct->GetPairingDistanceLimit()) ct->tem[j][i]=false;
			}
		}
	}

#ifndef _CUDA_CALC_
	//default behavior: calculate the partition function on the CPU
	calculatepfunction(ct,pfdata,progress,savefilename,false,&Q,w,v,wmb,wl,wmbl,wcoax,fce,w5,w3,mod,lfce);
#else //ifdef _CUDA_CALC_
	//if cuda flag is set, calculate on GPU
	//this requires compilation with nvcc
	//it will work for pair probabilities but it will break stochastic traceback

//read nearest neighbor parameters at desired temperature
	const char *path = getenv("DATAPATH");
	if (!path)
		die("%s: need to set environment variable $DATAPATH", "Could not find nearest neighbor parameters");
	struct param par;
	param_read_alternative_temperature(path, &par, pfdata->temp, !isrna);
//copy the sequence, correcting for 1-indexing
	char* buf = new char[ct->GetSequenceLength()+1];
	for(int i=0;i<GetSequenceLength();i++){
		buf[i] = ct->nucs[i+1];
	}
	buf[ct->GetSequenceLength()] = '\0'; //make it null-terminated
//calculate the partition function on the GPU
	prna_t p = prna_new(buf, &par);
//transfer over the partition function from the GPU data structure
//because the GPU partition function is calculated in log space, we do not directly transfer the partition function
//instead, we put the pair probs, calculated in log space in V, 1 everywhere else
//so RNA::GetPairProbability(i,j) will return V(i,j)*1/1, which will be the pairing probability
	for(int i=1;i<=ct->GetSequenceLength();i++){
		for(int j=i+1;j<=ct->GetSequenceLength();j++){
			v->f(i,j) = (PFPRECISION) probability_of_pair(p,i-1,j-1); //V(i,j) = P(i,j)
			v->f(j,i+ct->GetSequenceLength()) = 1.0; //V'(i,j) = 1
		}
		w5[i] = 1.0/(pfdata->scaling*pfdata->scaling);//have to divide by scaling*2
		w3[i] = 1.0/(pfdata->scaling*pfdata->scaling);//because RNA::calculateprobability expects scaling values
	}
	delete[] buf;
	prna_delete(p);
#endif
	if (savefilename!=NULL) {
		writepfsave(savefilename,ct,w5,w3,v,w,wmb,wl,wmbl,wcoax,fce,mod,lfce,pfdata);

		//clean up some memory use:
		delete[] savefilename;
	}
	return 0;
}


//Predict maximum expected accuracy structures that contain pseudoknots from either a sequence or a partition function save file.
int RNA::ProbKnot(int iterations, int MinHelixLength) {

	//first trap some possible errors
	if (!partitionfunctionallocated) {
		//There is no partition function data available.
		return 15;
	}

	if (iterations < 1) {
		//there can't be fewer than one iteration
		return 24;

	}



	//Past error trapping
	//Call the ProbKnot Program:
	return ProbKnotAssemble(v, w5, ct, pfdata, lfce, mod, pfdata->scaling, fce, iterations, MinHelixLength );


}

//Predict maximum expected accuracy structures that contain pseudoknots from a file containing ensemble of structures.
int RNA::ProbKnotFromSample(int iterations, int MinHelixLength) {

    if (iterations < 1) {
		//there can't be fewer than one iteration
		return 24;

	}

	//Past error trapping
	//Call the ProbKnot Program:
	return ProbKnotAssemble( ct, iterations, MinHelixLength );

}


//Refold a sequence using data from a save file.
int RNA::ReFoldSingleStrand(const float percent, const int maximumstructures, const int window) {

	if (!energyallocated) {
		//A .sav file was not read by the constructor.  Therefore, this function cannot be performed.
		return 17;

	}

	//Now do the refolding.
	return traceback(ct, data, ev, ew, ewmb, ew2, ewmb2, ew3, ew5, fce, lfce, vmin, maximumstructures, (int) percent, window,mod);

}


//Sample structures from the Boltzman ensemble.
int RNA::Stochastic(const int structures, const int seed) {

	if (!partitionfunctionallocated) {
		//There is no partition function data available.
		return 15;
	}


	//Past error trapping, call the stochastic traceback function
	return stochastictraceback(w,wmb,wmbl,wcoax,wl,v,
		fce, w3,w5,pfdata->scaling, lfce, mod, pfdata, structures,
		ct, seed, progress);


}


//Force a nucleotide to be double stranded (base paired).
//Return an integer that indicates an error code (0 = no error, 4 = nucleotide out of range, 8 = too many restraints specified, 9 = same nucleotide in conflicting restraint).
int RNA::ForceDoubleStranded(const int i) {
	int index;

	//check to make sure that a sequence has been read
	if (ct->GetSequenceLength()==0) return 20;

	//Check that nucleotide is valid
	if (i<1||i>ct->GetSequenceLength()) return 4;//i is out of range

	//Check for conflicting restraints (anything specifying i to be unpaired):
	for (index=0;index<ct->GetNumberofSingles();index++) {

		if(ct->GetSingle(index)==i)	return 9;
	}

	
	ct->AddDouble(i);
	
	return 0;

	


}

//Function to specify a nucleotide, i, that is a U in a GU pair
//Returns an integer that indicates an error code (0 = no error, 4 = nucleotide out of range, 8 = too many restraints specified, 9 = same nucleotide in conflicting restraint, 11 = nucleotide not U).
int RNA::ForceFMNCleavage(const int i) {
	int index;

	//check to make sure that a sequence has been read
	if (ct->GetSequenceLength()==0) return 20;

	//Check that nucleotide is valid
	if (i<1||i>ct->GetSequenceLength()) return 4;//i is out of range

	//Check to make sure the nucleotide is U.
	if (ct->numseq[i]!=4) return 11;

	//Check for conflicting restraints (anything specifying i to be unpaired):
	for (index=0;index<ct->GetNumberofSingles();index++) {

		if(ct->GetSingle(index)==i) return 9;
	}

	//Check for a nucleotide already forced to be in a pair that is not a GU pair.
	for (index=0;index<ct->GetNumberofPairs();index++) {

		if (i==ct->GetPair5(index)&&ct->numseq[ct->GetPair3(index)]!=3) return 9;
		else if (i==ct->GetPair3(index)&&ct->numseq[ct->GetPair5(index)]!=3) return 9;

	}


	ct->AddGUPair(i);
	
	return 0;

	


}

//Specify the maximum distance allowed between paired nucleotides in subsequent structure prediction.
//return An integer that indicates an error code (0 = no error, 12 = too long or too short).
int RNA::ForceMaximumPairingDistance(const int distance) {

	//check to make sure that a sequence has been read
	if (ct->GetSequenceLength()==0) return 20;

	if (distance < minloop+1) return 12;
	else {
		
		ct->SetPairingDistance(distance);
		return 0;
	}


}


//Indicate a nucleotide that is accessible to chemical modification.
//Returns an integer that indicates an error code (0 = no error, 4 = nucleotide out of range, 8 = too many restraints specified).
int RNA::ForceModification(const int i) {

	//check to make sure that a sequence has been read
	if (ct->GetSequenceLength()==0) return 20;

	//Check that nucleotide is valid
	if (i<1||i>ct->GetSequenceLength()) return 4;//i is out of range


	


	
	//Go ahead and record the constraint.
	ct->AddModified(i);
	
	return 0;
	


}

//Force a base pair between nucleotides i and j.
//Returns an error code: (0 = no error, 4 = nucleotide out of range, 6 = pseudoknot formation, 7 = non-canonical pair, 8 = too many restraints specified, 9 = same nucleotide in conflicting restraint).
int RNA::ForcePair(const int i, const int j) {
	bool allowedpairs[6][6]={{false,false,false,false,false,false},{false,false,false,false,true,false},{false,false,false,true,false,false},{false,false,true,false,true,false},
	{false,true,false,true,false,false},{false,false,false,false,false,false}};
	int index;
	int locali,localj;

	//First perform the error checking:

	//check to make sure that a sequence has been read
	if (ct->GetSequenceLength()==0) return 20;

	//Note: In structure, forced pairs run between index of 1 and a maximum of maxforce-1.

	//if (ct->npair==(maxforce-1)) return 8;//This means there are too many pair constraints.

	if (i<1||i>ct->GetSequenceLength()) return 4;//i is out of range
	if (j<1||j>ct->GetSequenceLength()) return 4;//j is out of range

	if (!allowedpairs[ct->numseq[i]][ct->numseq[j]]) return 7;//non-canonical pair

	//sort indexes from 5' to 3':
	if (i>j) {
		locali=j;
		localj=i;
	}
	else {
		locali=i;
		localj=j;
	}

	//check for pseudoknots with any other forced pair or the same nucleotide forced into two pairs:
	for (index=0;index<ct->GetNumberofPairs();index++) {
		if (locali<ct->GetPair5(index)&&ct->GetPair5(index)<localj&&localj<ct->GetPair3(index)) return 6;//a pseudoknot

		if (locali==ct->GetPair5(index)||locali==ct->GetPair3(index)||localj==ct->GetPair5(index)||localj==ct->GetPair3(index)) return 9;//i or j is in another forced pair

	}

	//now check for other conflicting restraints:
	for (index=0;index<ct->GetNumberofForbiddenPairs();index++) {
		if(ct->GetForbiddenPair5(index)==locali && ct->GetForbiddenPair3(index)==localj ) return 9;//The pair was forbidden.

	}
	for (index=0;index<ct->GetNumberofSingles();index++) {

		if(ct->GetSingle(index)==locali||ct->GetSingle(index)==localj)	return 9;//i or j was previously forced single-stranded.
	}

	//Now register the restraint because the error checking was clear or errors.
	ct->AddPair(locali,localj);

	return 0;

}

//Prohibit a pair between two nucleotides in subsequent structure prediction.
//Returns an integer that indicates an error code (0 = no error, 4 = nucleotide out of range, 8 = too many restraints specified, 9 = nucleotide in conflicting restraint).
int RNA::ForceProhibitPair(const int i, const int j) {
	int index,locali,localj;

	//First perform the error checking:

	//Note: In structure, forced pairs run between index of 0 and a maximum of maxforce-1.

	//check to make sure that a sequence has been read
	if (ct->GetSequenceLength()==0) return 20;


	if (i<1||i>ct->GetSequenceLength()) return 4;//i is out of range
	if (j<1||j>ct->GetSequenceLength()) return 4;//j is out of range


	//sort indexes from 5' to 3':
	if (i>j) {
		locali = j;
		localj = i;
	}
	else {
		locali=i;
		localj=j;
	}

	//check to make sure this pair hasn't been forced:
	for (index=0;index<ct->GetNumberofPairs();index++) {

		if (locali==ct->GetPair5(index)&&localj==ct->GetPair3(index)) return 9;//i or j is in a forced pair

	}


	//Now register the restraint because the error checking was clear or errors.
	ct->AddForbiddenPair(locali,localj);

	return 0;

}

//Force a nucleotide to be single stranded in subsequent structure prediction.
//An integer that indicates an error code (0 = no error, 4 = nucleotide out of range, 8 = too many restraints specified, 9 = same nucleotide in conflicting restraint).
int RNA::ForceSingleStranded(const int i) {
	int index;

	//check to make sure that a sequence has been read
	if (ct->GetSequenceLength()==0) return 20;

	//Check that nucleotide is valid
	if (i<1||i>ct->GetSequenceLength()) return 4;//i is out of range

	//Check for conflicting constraints; anything forcing a nucleotide to be paired.

	for (index=0;index<ct->GetNumberofPairs();index++) {//check all the forced pairs
		if (i==ct->GetPair5(index)||i==ct->GetPair3(index)) return 9;//i is in a forced pair
	}
	for (index=0;index<ct->GetNumberofDoubles();index++) {//check all the force doubles
		if(ct->GetDouble(index)==i)	return 9;
	}
	for (index=0;index<ct->GetNumberofGU();index++) {//check all the force FMN
		if(ct->GetGUpair(index)==i)	return 9;
	}


	//Register the constraint:
	ct->AddSingle(i);

	return 0;

}

//Return a nucleotide that is forced double stranded.
int RNA::GetForcedDoubleStranded(const int constraintnumber) {

	//First make sure the constraintnumber is valid.
	if (constraintnumber<0||constraintnumber>=ct->GetNumberofDoubles()) return 0;

	//Now return the constraint.
	return ct->GetDouble(constraintnumber);

}

//Return a nucleotide that is accessible to FMN cleavage.
int RNA::GetForcedFMNCleavage(const int constraintnumber) {

	//First make sure the constraintnumber is valid.
	if (constraintnumber<0||constraintnumber>=ct->GetNumberofGU()) return 0;

	//Now return the constraint.
	return ct->GetGUpair(constraintnumber);

}

//Return a nucleotide that is accessible to modification.
int RNA::GetForcedModification(const int constraintnumber) {

	//First make sure the constraintnumber is valid.
	if (constraintnumber<0||constraintnumber>=ct->GetNumberofModified()) return 0;

	//Now return the constraint.
	return ct->GetModified(constraintnumber);//note that the underlying ct indexes from 1 to ndbl.

}

//Return a nucleotide in a forced pair.
//fiveprime determines if the nucleotide is the five prime or the three prime nucleotide in the constraint.  true = five prime nucleotide.
int RNA::GetForcedPair(const int constraintnumber, const bool fiveprime) {

	//First make sure the constraintnumber is valid.
	if (constraintnumber<0||constraintnumber>=ct->GetNumberofPairs()) return 0;

	//Now return the constraint.
	if (fiveprime) return ct->GetPair5(constraintnumber);
	else return ct->GetPair3(constraintnumber);

}

//Return a nucleotide in a prohibited pair.
//fiveprime determines if the nucleotide is the five prime or the three prime nucleotide in the constraint.  true = five prime nucleotide.
int RNA::GetForcedProhibitedPair(const int constraintnumber, const bool fiveprime) {

	//First make sure the constraintnumber is valid.
	if (constraintnumber<0||constraintnumber>=ct->GetNumberofForbiddenPairs()) return 0;

	//Now return the constraint.
	if (fiveprime) return ct->GetForbiddenPair5(constraintnumber);
	else return ct->GetForbiddenPair3(constraintnumber);

}

//Return a nucleotide that is forced single stranded.
int RNA::GetForcedSingleStranded(const int constraintnumber) {

	//First make sure the constraintnumber is valid.
	if (constraintnumber<0||constraintnumber>=ct->GetNumberofSingles()) return 0;

	//Now return the constraint.
	return ct->GetSingle(constraintnumber);//note that the underlying ct indexes from 1 to ndbl.


}


//Return the maximum pairing distance.
//Return an integer that indicates the maximum distance allowed between paired nucleotides, where -1 indicates that the maximum distance is not set.
int RNA::GetMaximumPairingDistance() {

	if (ct->DistanceLimited()) return ct->GetPairingDistanceLimit();
	else return -1;

}


//Return the number of nucletides forced to be paired.
int RNA::GetNumberOfForcedDoubleStranded() {

	return ct->GetNumberofDoubles();

}

// Add an experimental bonus to a pair of nucleotides
//void RNA::SetExperimentalBonus(const int i, const int j, const double bonus){

//	ct->EX[i][j] = bonus;

//}


//!Return the number of nucleotides accessible to FMN cleavage.
int RNA::GetNumberOfForcedFMNCleavages() {

	return ct->GetNumberofGU();

}

//!Return the number of nucleotides accessible to chemical modification.
int RNA::GetNumberOfForcedModifications() {

	return ct->GetNumberofModified();

}

//!Return the number of forced base pairs.
int RNA::GetNumberOfForcedPairs() {

	return ct->GetNumberofPairs();

}

//!Return the number of prohibited base pairs.
int RNA::GetNumberOfForcedProhibitedPairs() {

	return ct->GetNumberofForbiddenPairs();

}

//!Return the number of nucleotides that are not allowed to pair.
int RNA::GetNumberOfForcedSingleStranded() {

	return ct->GetNumberofSingles();

}

//Read a set of folding constraints to disk in a plain text file.
//filename is a c string that is the file name to be read.
//Returns an integer that indicates an error code (0 = no error, 1 = file not found, 13 = error reading constraint file).
int RNA::ReadConstraints(const char filename[]) {
	FILE *check;

	//check that the file exists.
	if ((check = fopen(filename, "r"))== NULL) {
		//the file is not found
		fclose(check);
		return 1;
	}
	fclose(check);

	//Now read the constraints
	if (readconstraints(filename, ct)) return 0;
	else return 13;



}

//Read SHAPE data to constrain structure prediction on subsequent structure predictions.
//filename is a c string that indicates a file that contains SHAPE data.
//IsPseudoEnergy indicates whether this is the pseudo folding free energy constraint (the preferred method).  This defaults to true.
//parameter1 is the slope when IsPseudoEnergy=true and is a threshold above which nucleotides are forced single stranded otherwise.
//parameter2 is the intercept when IsPseudoEnergy=true and is a threshold above which a nucleotide is considered chemically modified otherwise.
//modifier is the type of chemical modification probe that was used (currently accepted values are SHAPE, diffSHAPE, DMS, and CMCT). Defaults to SHAPE.
//Returns an integer that indicates an error code (0 = no error, 1 = input file not found).
int RNA::ReadSHAPE(const char filename[], const double parameter1, const double parameter2, std::string modifier, const bool IsPseudoEnergy) {
	FILE *check;

	//check that the SHAPE input file exists
	if ((check = fopen(filename, "r"))== NULL) {
		//the file is not found
		fclose(check);
		return 1;
	}

	fclose(check);

	if (IsPseudoEnergy) {
		//This is the pseudo energy version

		ct->SHAPEslope=parameter1*conversionfactor;//register the slope in tenths of kcal/mol
		ct->SHAPEintercept=parameter2*conversionfactor;//register the intercept in tenths of a kcal/mol
		ct->ReadSHAPE(filename, modifier);//call ReadSHAPE() to read the file and determine pseudo energies

	}
	else {
		ct->ReadSHAPE(filename,(float) parameter1,(float) parameter2);//call ReadSHAPE() with parameters to parse thresholds
	}
	return 0;


}


//Read SHAPE data to constrain structure prediction on subsequent structure predictions.
//filename is a c string that indicates a file that contains SHAPE data.
//parameter1 is the slope.
//parameter2 is the intercept.
//Returns an integer that indicates an error code (0 = no error, 1 = input file not found).
int RNA::ReadExperimentalPairBonus(const char filename[], double const experimentalOffset, double const experimentalScaling ) {
	FILE *check;

	//check that the SHAPE input file exists
	if ( strlen( filename ) > 0  ) {
	  if ( (check = fopen(filename, "r"))== NULL) {
	    //the file is not found
	    fclose(check);
	    return 1;
	  }

	  fclose(check);
	}

	ct->ReadExperimentalPairBonus(filename, experimentalOffset, experimentalScaling );

	return 0;

}


//Read SHAPE data to constrain structure prediction on subsequent structure predictions - overloaded version for including single-stranded SHAPE.
//filename is a c string that indicates a file that contains SHAPE data.
//parameter1 is the double-stranded slope.
//parameter2 is the double-stranded intercept.
//modifier is the type of chemical modification probe that was used (currently accepted values are SHAPE, DMS, and CMCT). Defaults to SHAPE.
//ssm is the single-stranded slope.
//ssb in the single-stranded intercept.
//Returns an integer that indicates an error code (0 = no error, 1 = input file not found).
int RNA::ReadSHAPE(const char filename[], const double parameter1, const double parameter2, const double ssm, const double ssb, std::string modifier) {
	FILE *check;

	//check that the SHAPE input file exists
	if ((check = fopen(filename, "r"))== NULL) {
		//the file is not found
		fclose(check);
		return 1;
	}

	fclose(check);



	ct->SHAPEslope=parameter1*conversionfactor;//register the slope in tenths of kcal/mol
	ct->SHAPEintercept=parameter2*conversionfactor;//register the intercept in tenths of a kcal/mol
	ct->SHAPEslope_ss=ssm*conversionfactor;//register the slope in tenths of kcal/mol
	ct->SHAPEintercept_ss=ssb*conversionfactor;//register the intercept in tenths of a kcal/mol
	ct->ReadSHAPE(filename, modifier);//call ReadSHAPE() to read the file and determine pseudo energies


	return 0;


}

//Read Double Strand Offset
int RNA::ReadDSO(const char filename[]) {
	FILE *check;

	//check that the SHAPE input file exists
	if ((check = fopen(filename, "r"))== NULL) {
		//the file is not found
		fclose(check);
		return 1;
	}

	ct->ReadOffset(NULL,filename);

	return 0;

}

//Read Single Strand Offset
int RNA::ReadSSO(const char filename[]) {
	FILE *check;

	//check that the SHAPE input file exists
	if ((check = fopen(filename, "r"))== NULL) {
		//the file is not found
		fclose(check);
		return 1;
	}

	ct->ReadOffset(filename,NULL);

	return 0;

}


//Remove all previously defined constraints.
void RNA::RemoveConstraints() {

	ct->RemoveConstraints();


	ct->min_gu=0;
	ct->min_g_or_u=0;
	ct->nneighbors=0;
	ct->nregion=0;
	ct->nmicroarray=0;

}

//Add extrinsic restraints for partition function calculations.
int RNA::SetExtrinsic(int i, int j, double k) {
	int locali, localj;


	//First do the error checking:
	//check the indices
	if (i<1||i>ct->GetSequenceLength()||j<1||j>ct->GetSequenceLength()) return 4;

	//make sure the equilibrium constant is not less than zero.
	if (k<0) return 26;

	//Now past error checking:

	//sort indexes from 5' to 3':
	if (i>j) {
		locali = j;
		localj = i;
	}
	else {
		locali=i;
		localj=j;
	}



	if (ct->constant==NULL) {
		//allocate the space needed in structure to store the constants
		ct->allocateconstant();
	}


	ct->constant[localj][locali] = k;

	return 0;
}


//Write the set of folding constraints to disk.
int RNA::WriteConstraints(const char filename[]) {

	outputconstraints(filename,ct);
	return 0;

}



//Specify a comment for inclusion in subsequently written .ct files.
int RNA::AddComment(const char comment[], const int structurenumber) {
	string label;

	//start with error checking:
	if (structurenumber<1||structurenumber>ct->GetNumberofStructures()) return 3;

	//now register the comment (at the end of existing string, but before the newline there by default):
	label = ct->GetCtLabel(structurenumber);
	
	//Remove the existing newline character, if it exists:
	if (label.length()>0) if (label[label.length()-1]=='\n') label.erase(label.length()-1);

	//Add the comment
	label+=comment;

	//Add back the newline
	label+="\n";

	//add a newline at the end of the comment -- required when ct is written
	ct->SetCtLabel(label,structurenumber);

	return 0;

}

//Write a ct file of the structures
int RNA::WriteCt(const char filename[], bool append) {
	if (ct->GetNumberofStructures()>0) {
		ct->ctout(filename,append);
		return 0;
	}
	else return 10; //an error code


}

//Write a ct file of the structures
int RNA::WriteDotBracket(const char filename[]) {
	if (ct->GetNumberofStructures()>0) {
		ct->writedotbracket(filename);
		return 0;
	}
	else return 10; //an error code


}


//Break any pseudoknots that might be in a structure.
int RNA::BreakPseudoknot(const bool minimum_energy, const int structurenumber) {

	int i,j,structures;
	structure *tempct;


	double **bpProbArray; //contains the raw bp probabilities for each bp
	double *bpSSProbArray; //contains the raw single strand probability for a base
	double **vwArray;  //contains v and w recursion values
	double **vwPArray; //the v' and w' recursion values

	double *w3Array=0;//w3Array[i] is the maximum score from nucletides i to ct->GetSequenceLength()
	double *w5Array=0;//w5Array[i] is the maximum score from nucleotides 1 to i
							

	int Length;
	int start,stop;


	double sumPij; //holds the sum of probabilities for base pairs based on a specific i over js


	//Make sure there are structures...
	if (ct->GetNumberofStructures()<=0) {
		return 23;
	}

	//If a specific structure is requested, i.e. structurenumber!=0, make sure the structure is valid.
	if (structurenumber!=0) {
		if (structurenumber<0||structurenumber>ct->GetNumberofStructures()) return 3;

	}

	if (minimum_energy) {
		//mnimum_energy==true, so minimize free energy in the pseudoknot free structure


		if (!energyread) {
			//The thermodynamic data tables have not been read and need to be read now.
			if (ReadThermodynamic()!=0) return 5;//return non-zero if a problem occurs

		}

		//allocate tempct:
		tempct = new structure(2);
		tempct->allocate(ct->GetSequenceLength());
		tempct->SetSequenceLabel("temp\n");
		//tempct->numofbases = ct->GetSequenceLength();

		for (i=1;i<=ct->GetSequenceLength();i++) {
			tempct->numseq[i]=ct->numseq[i];

		}
		//allocate a template of allowed pairs
		tempct->allocatetem();

		//Break pairs for each structure (first check if only one structure should be treated)

		if (structurenumber!=0) {
			start = structurenumber;
			stop = structurenumber;

		}
		else {
			start = 1;
			stop = ct->GetNumberofStructures();
		}
		for (structures=start;structures<=stop;structures++) {

			//initialize all base pairs as unallowed:
			for (i=0;i<=ct->GetSequenceLength();i++) {
				for (j=i+1;j<=ct->GetSequenceLength();j++) {
    				tempct->tem[j][i] = false;
				}
		   }

			//now allow the pairs that are in the loaded ct:
			for (i=1;i<=ct->GetSequenceLength();i++) {
				if (ct->GetPair(i,structures)>i) {
					tempct->tem[ct->GetPair(i,structures)][i] = true;
				}
			}

			if (tempct->GetNumberofStructures() > 0) {
				//strip the structure from the ct at this point
				for (int index=tempct->GetNumberofStructures();index>0;--index) tempct->RemoveLastStructure();

			}

			//Predict the secondary structures.
			alltrace(tempct,data,0,0,progress,NULL,false);//
			//dynamic(tempct, data, 1, 0, 0, progress, false, NULL, 30);



			//copy the pairs back to ct
			ct->CleanStructure(structures);
			for (i=1;i<=ct->GetSequenceLength();i++) {
				if (tempct->GetPair(i)>i) ct->SetPair(i,tempct->GetPair(i),structures);

			}

			//Also copy the energy back to ct
			ct->SetEnergy(structures,tempct->GetEnergy(1));



		}
		delete tempct;

	}//end of minimum_energy == true
	else {
		//This is minimum_energy == false, so maximize pairs

		//allocate tempct:
		tempct = new structure(2);
		tempct->allocate(ct->GetSequenceLength());
		tempct->SetSequenceLabel("temp\n");
		//tempct->numofbases = ct->GetSequenceLength();

		for (i=1;i<=ct->GetSequenceLength();i++) {
			tempct->numseq[i]=ct->numseq[i];

		}


		//loop over all structures unless the programmer requested a specific structure
		if (structurenumber!=0) {
			start = structurenumber;
			stop = structurenumber;

		}
		else {
			start = 1;
			stop = ct->GetNumberofStructures();
		}
		for (structures=start;structures<=stop;structures++) {



			 
			//strip the structure from the ct at this point
			if (tempct->GetNumberofStructures()>0) {
				for(int index=tempct->GetNumberofStructures();index>0;index--) tempct->RemoveLastStructure();
			}



			//allocate main arrays and initialize the Arrays to defaults
			bpProbArray = new double *[ct->GetSequenceLength()+1];
			bpSSProbArray = new double [ct->GetSequenceLength()+1];
			vwArray = new double *[ct->GetSequenceLength()+1];
			vwPArray = new double *[ct->GetSequenceLength()+1];


			sumPij = 0;

			for (i=0;i<=ct->GetSequenceLength();i++) {
				bpProbArray[i] = new double [ct->GetSequenceLength()+1];
				vwArray[i] = new double [ct->GetSequenceLength()+1];
				vwPArray[i] = new double [ct->GetSequenceLength()+1];

				//tempct->basepr[1][i] = 0;
				bpSSProbArray[i] = 0;

				for (j=0;j<=ct->GetSequenceLength();j++) {
					bpProbArray[i][j]=0;
					vwArray[i][j]=-0;
					vwPArray[i][j]=-0;
				}
			}




			tempct->nucs[0] = ' ';

			// Recursion rules investigate for vwArray
			//    1)  if the base pair (BP) probability value is 0, skip that pair
			//    2)  hairpin turns at 5 BP
			//    3)  stack/internal/bulge pairing at 7 BPs
			//    4)  multibranching at 12 BPs (2 hairpins and a stack)
			// Because of the rules for hairpin, the probabilities
			//    can be taken from the bpProbArray directly until j-i > 5

			// Calculate the single stranded probabilities for each base
			// Pi = 1 - (for all j, sum(Pij)
			// fill in w for the diagonal for the Pi,i
			for (i=1; i<=ct->GetSequenceLength(); i++)
			{


				if (ct->GetPair(i,structures)!=0) bpSSProbArray[i] = 0.0;
				else bpSSProbArray[i] = 1.0;



				vwArray[i][i] = bpSSProbArray[i];
			} // end loop over each base pair


			//Calculate the base pair probabilities to start...
			for (Length = 2; Length <=ct->GetSequenceLength(); Length++)
			{

				//begin populating v and w along the diagonal starting with the
				//   shortest hairpin length
				for (i=1, j=i+Length-1; j<=ct->GetSequenceLength(); i++, j++)
				{
					if (ct->GetPair(i,structures)==j) {
						bpProbArray[j][i]=1.0;
					}
					else {
						bpProbArray[j][i]=-1.0;
					}

				}
			}

			//Call the MEAFill routine.
				//Note the false at the end "allows" non-canonical pairs.  This is required so that
				//non-canonical pairs aren't spuriosly broken
			MEAFill(tempct, bpProbArray, bpSSProbArray, vwArray, vwPArray, w5Array, w3Array, 1.0, 0, progress,false);



			// start traceback
			trace(tempct, vwArray, vwPArray, bpProbArray, 1.0, 0, 1, 0);






			// Deallocate memory for the MaxExpect calculation
			//Arrays with functionality in the fill step
			for (i=0;i<=ct->GetSequenceLength();i++) {
				delete[] bpProbArray[i];
			}
			delete[] bpProbArray;

			delete[] bpSSProbArray;

			for (i=0; i<=ct->GetSequenceLength(); i++) {
				delete[] vwArray[i];
				delete[] vwPArray[i];
			}
			delete[] vwArray;
			delete[] vwPArray;




			//copy the pairs back to ct

			//clean the structure of pairs first
			ct->CleanStructure(structures);

			//loop over all bases to find pairs
			for (i=1;i<=ct->GetSequenceLength();i++) {
				if(tempct->GetPair(i,1)>i) ct->SetPair(i,tempct->GetPair(i,1),structures);

			}

			//Also copy the energy back to ct
			ct->SetEnergy(structures,tempct->GetEnergy(1));



		}
		delete tempct;



	}
	return 0;





}


// Report if there are any pseudoknots in a structure.


bool RNA::ContainsPseudoknot(const int structurenumber) {
	int i,j;

	//make sure structurenumber is a valid structure
	if (structurenumber<1||structurenumber>ct->GetNumberofStructures()) {
		ErrorCode=3;
		return false;
	}
	else {
		//passed error trapping:

		//check all nucs for pairs
		for (i=1;i<ct->GetSequenceLength();i++) {

			if (ct->GetPair(i,structurenumber)>i) {
				//found pair, check for crossing pairs
				for (j=i+1;j<ct->GetPair(i,structurenumber);j++) {

					if (ct->GetPair(j,structurenumber)>j) {

						if (ct->GetPair(j,structurenumber)>ct->GetPair(i,structurenumber)) {

							return true;

						}

					}
				}
			}
		}

		return false;//no pseudoknot was found
	}
}


//Get the ensemble folding free energy change as determined by the partition function.
double RNA::GetEnsembleEnergy() {

	//check to see if partition function data is present.
	if (!partitionfunctionallocated) {
		ErrorCode = 15;
		return 0.0;
	}

	//past error trapping, so set no error found
	ErrorCode = 0;

	//calculate the ensemble folding free energy change, -RT ln (Q).
	//One needs to also account for the fact that a scaling is applied per nucleotide.
	return (double) ((-RKC*temp)*(log(w5[ct->GetSequenceLength()])-ct->GetSequenceLength()*log(pfdata->scaling)));

}

//Get the folding free energy change for a predicted structure.
double RNA::GetFreeEnergy(const int structurenumber) {

	//make sure structurenumber is a valid structure
	if (structurenumber<1||structurenumber>ct->GetNumberofStructures()) {
		ErrorCode=3;
		return 0.0;
	}
	else {
		//error trapping complete
		return (((double) ct->GetEnergy(structurenumber))/((double) conversionfactor));

	}

}

// Get the nucleotide to which the specified nucleotide is paired.
int RNA::GetPair(const int i, const int structurenumber) {

	//check to make sure i is a valid nucleotide
	if (i<1||i>ct->GetSequenceLength()) {
		ErrorCode=4;
		return 0;
	}
	//now make sure structurenumber is a valid structure
	else if (structurenumber<1||structurenumber>ct->GetNumberofStructures()) {
		ErrorCode=3;
		return 0;
	}
	else {
		//error trapping complete
		return ct->GetPair(i,structurenumber);

	}

}

//Extract the lowest free energy for a structure that contains the i-j pair using data from a save file (.sav).
double RNA::GetPairEnergy(const int i, const int j) {
	int locali, localj;

	//check whether the a save file was read (by the constructor)
	if (!energyallocated) {
		ErrorCode = 17;
		return 0.0;
	}


	if (i<1||i>ct->GetSequenceLength()) {
		//i is out of range
		ErrorCode = 4;
		return 0.0;
	}
	if (j<1||j>ct->GetSequenceLength()) {
		//j is out of range
		ErrorCode = 4;
		return 0.0;
	}


	//sort indexes from 5' to 3':
	if (i>j) {
		locali = j;
		localj = i;
	}
	else {
		locali=i;
		localj=j;
	}

	//No error;
	ErrorCode = 0;

	//calculate and return the energy
	return ((((double) (ev->f(locali,localj)+ev->f(localj,locali+ct->GetSequenceLength())))/conversionfactor));



}

//Get the total number of specified structures
int RNA::GetStructureNumber() {

		return ct->GetNumberofStructures();

}

//return the base pairing probability between nucleotides i and j
double RNA::GetPairProbability(const int i, const int j) {

	//check to see if partition function data is present.
	if (!partitionfunctionallocated) {
		ErrorCode = 15;
		return 0.0;
	}

	//check that the nucleotides are in the correct range
	if (i<1||j>ct->GetSequenceLength()||j<0||j>ct->GetSequenceLength()) {
		ErrorCode = 4;
		return 0.0;
	}

	//past error trapping, so set no error found
	ErrorCode = 0;

	//calculate the base pair probability
	return (double) calculateprobability(i,j,v,w5,ct,pfdata,lfce,mod,pfdata->scaling,fce);


}

//Determine the coordinates for drawing a secondary structure.
int RNA::DetermineDrawingCoordinates(const int height, const int width, const int structurenumber) {

	//check to make sure that a sequence has been read
	if (ct->GetSequenceLength()==0) return 20;

	//First check for errors in the specification of the structure number:
	if (structurenumber<0||structurenumber>ct->GetNumberofStructures()) return 3;

	if (!drawallocated) {
		//If the memory has not been allocated for coordinates, go ahead and allocate it now
		structurecoordinates = new coordinates(ct->GetSequenceLength());
		drawallocated = true;
	}

	//now perform the calculation for coordinate determination:
	place(structurenumber, ct, structurecoordinates, height, width);

	//no errors, return 0
	return 0;


}

//Provide the comment from the ct file as a string.
std::string RNA::GetCommentString(const int structurenumber) {
	std::string temp;


	//Add some code for backwards compatibility:
	//In the past, all labels were assocaited with structures.  Now, labels can come with the sequence.
	//If there are no structures, return the sequence label:

	if (ct->GetNumberofStructures()==0) {

		temp = ct->GetSequenceLabel();
		return temp;
	}



	//start with error checking:
	if (structurenumber<1||structurenumber>ct->GetNumberofStructures()){
		//The request is for a structure that is out of range
		ErrorCode = 3;
		temp = "";

		return temp;

	}

	temp = ct->GetCtLabel(structurenumber);

	return temp;

}

//Get the X coordinate for nucleotide i for drawing a structure.
int RNA::GetNucleotideXCoordinate(const int i) {

	if (!drawallocated) {
		//The drawing coordinates were not pre-determined by DetermineDrawingCoordinates(), so indicate an error.
		ErrorCode = 19;
		return 0;
	}

	if (i<0||i>ct->GetSequenceLength()) {
		//The nucleotide is invalid, indicate an error.
		ErrorCode = 4;
		return 0;
	}

	//Fetch the coordinate from the coordinates structure, structurecoordinates.
	return structurecoordinates->x[i];


}

//Get the Y coordinate for nucleotide i for drawing a structure.
int RNA::GetNucleotideYCoordinate(const int i) {

	if (!drawallocated) {
		//The drawing coordinates were not pre-determined by DetermineDrawingCoordinates(), so indicate an error.
		ErrorCode = 19;
		return 0;
	}

	if (i<0||i>ct->GetSequenceLength()) {
		//The nucleotide is invalid, indicate an error.
		ErrorCode = 4;
		return 0;
	}

	//Fetch the coordinate from the coordinates structure, structurecoordinates.
	return structurecoordinates->y[i];


}

// Get the X coordinate for placing the nucleotide index label specified by i.
int RNA::GetLabelXCoordinate(const int i) {

	if (!drawallocated) {
		//The drawing coordinates were not pre-determined by DetermineDrawingCoordinates(), so indicate an error.
		ErrorCode = 19;
		return 0;
	}

	if (i<0||i>ct->GetSequenceLength()) {
		//The nucleotide is invalid, indicate an error.
		ErrorCode = 4;
		return 0;
	}

	if (i%10!=0) {
		//The nucleotide index is not a multiple of 10, return an error
		ErrorCode = 25;
		return 0;

	}

	//Fetch the coordinate from the coordinates structure, structurecoordinates.
	return structurecoordinates->num[i/10][0];


}

// Get the Y coordinate for placing the nucleotide index label specified by i.
int RNA::GetLabelYCoordinate(const int i) {

	if (!drawallocated) {
		//The drawing coordinates were not pre-determined by DetermineDrawingCoordinates(), so indicate an error.
		ErrorCode = 19;
		return 0;
	}

	if (i<0||i>ct->GetSequenceLength()) {
		//The nucleotide is invalid, indicate an error.
		ErrorCode = 4;
		return 0;
	}

	if (i%10!=0) {
		//The nucleotide index is not a multiple of 10, return an error
		ErrorCode = 25;
		return 0;

	}

	//Fetch the coordinate from the coordinates structure, structurecoordinates.
	return structurecoordinates->num[i/10][1];

}

//Get the identity of nucleotide i.
char RNA::GetNucleotide(const int i) {

	//check to make sure that a sequence has been read
	if (ct->GetSequenceLength()==0) {
		ErrorCode = 20;
		return '-';
	}

	//Check that nucleotide is valid
	if (i<1||i>ct->GetSequenceLength()) {
		ErrorCode = 4;//i is out of range
		return '-';
	}

	return ct->nucs[i];

}

//Get the total length of the sequence
int RNA::GetSequenceLength() {

	return ct->GetSequenceLength();

}

//Return the type of backbone (true = RNA, false = DNA).
bool RNA::GetBackboneType() {

	return isrna;

}






//Access the underlying structure class.
structure *RNA::GetStructure() {

	return ct;

}


//Provide a TProgressDialog for following calculation progress.
//A TProgressDialog class has a public function void update(int percent) that indicates the progress of a long calculation.
void RNA::SetProgress(TProgressDialog& Progress) {


	progress = &Progress;

	return;

}


//Provide a means to stop using a TProgressDialog.
//StopProgress tells the RNA class to no longer follow progress.  This should be called if the TProgressDialog is deleted, so that this class does not make reference to it.
void RNA::StopProgress() {

	progress=NULL;
	return;

}

TProgressDialog* RNA::GetProgress() {

		return progress;

}


RNA::~RNA() {






	if (partitionfunctionallocated) {
		//The partition function calculation was performed, so the memory allocated for the partition function needs to be deleted.
		delete[] lfce;
		delete[] mod;
		delete[] w5;
		delete[] w3;
		delete v;
		delete w;
		delete wmb;
		delete wl;
		delete wmbl;
		delete wcoax;
		delete fce;
		delete pfdata;

	}

	if (energyallocated) {
		//A folding save file was opened, so clean up the memory use.

		delete[] lfce;
		delete[] mod;



		delete[] ew5;
		delete[] ew3;


		if (ct->intermolecular) {
			delete ew2;
			delete ewmb2;
		}

		delete ev;
		delete ew;
		delete ewmb;
		delete fce;

	}

	if (drawallocated) {
		//The drawing coordiantes have been determined, so the memory needs to be cleaned up.
		delete structurecoordinates;

	}

	delete ct;//delete the structure


}



//This is a protected function for handling file input.
int RNA::FileReader(const char filename[], const int type) {
	FILE *check;
	short vers;
	//double scaling;
	int i;



	if ((check = fopen(filename, "r"))== NULL) {
		//the file is not found
		return 1;

	}
	else {
		fclose(check);
		//open the file based on type:
		if (type==1) {
			//type indicates a ct file
			long linenumber = ct->openct(filename);

			//if linenumber==0, there was no error
			if (linenumber==0) return 0;
			//Otherwise, there was an error opening the file
			else return 2;

		}
		else if (type==2) {
			//type indicates a .seq file

			//ct->numofstructures=0;
			if (ct->openseq(filename)==1) return 0;//no error
			else return 2;//File open error



		}
		else if (type==3) {
			//type indicates a partition function save file

			//allocate the ct file by reading the save file to get the sequence length:
			std::ifstream sav(filename,std::ios::binary);

			read(&sav,&(vers));//read the version of the save file

			if (vers!=pfsaveversion) {
				//Wrong version!
				sav.close();
				return 16;
			}

			int sequencelength;
			//read the length of the sequence
			read(&sav,&(sequencelength));
			sav.close();

			//allocate everything
			ct->allocate(sequencelength);

			w = new pfunctionclass(ct->GetSequenceLength());
			v = new pfunctionclass(ct->GetSequenceLength());
			wmb = new pfunctionclass(ct->GetSequenceLength());
			wmbl = new pfunctionclass(ct->GetSequenceLength());
			wcoax = new pfunctionclass(ct->GetSequenceLength());
			wl = new pfunctionclass(ct->GetSequenceLength());
			fce = new forceclass(ct->GetSequenceLength());

			w5 = new PFPRECISION [ct->GetSequenceLength()+1];
			w3 = new PFPRECISION [ct->GetSequenceLength()+2];

			lfce = new bool [2*ct->GetSequenceLength()+1];
			mod = new bool [2*ct->GetSequenceLength()+1];

			pfdata = new pfdatatable();

			//indicate that the memory has been allocated so that the destructor will delete it.
			partitionfunctionallocated = true;

			//load all the data from the pfsavefile:
			readpfsave(filename, ct, w5, w3,v, w, wmb,wl, wmbl, wcoax, fce,&pfdata->scaling,mod,lfce,pfdata);
			return 0;

		}
		else if (type==4) {
			//Type indicates restoration of a folding save file.


			//peek at the length of the sequence and whether the folding is intermolecular to allocate arrays:
			std::ifstream sav(filename,std::ios::binary);

			//Read the save file version information.
			read(&sav,&vers);

			//check the version
			if (vers!=safiversion) {
				//Wrong version!
				sav.close();
				return 16;
			}

			int sequencelength;
			read(&sav,&(sequencelength));
			read(&sav,&(ct->intermolecular));
			sav.close();

			//indicate that everything is allocated and needs to be deleted in the destructor
			energyallocated = true;

			//allocate everything
			ct->allocate(sequencelength);

			ew = new arrayclass(ct->GetSequenceLength());
			ev = new arrayclass(ct->GetSequenceLength());
			ewmb = new arrayclass(ct->GetSequenceLength());
			fce = new forceclass(ct->GetSequenceLength());


			lfce = new bool [2*ct->GetSequenceLength()+1];
			mod = new bool [2*ct->GetSequenceLength()+1];

			ew5 = new integersize [ct->GetSequenceLength()+1];
			ew3 = new integersize [ct->GetSequenceLength()+2];

			if (ct->intermolecular) {
				ew2 = new arrayclass(ct->GetSequenceLength());
				ewmb2 = new arrayclass(ct->GetSequenceLength());

				for (i=0;i<3;i++) read(&sav,&(ct->inter[i]));

			}
			else {
				ew2 = NULL;
				ewmb2 = NULL;
			}

			//indicate that the thermodynamic parameters are read and available (and need to be deleted).
			energyread = true;
			data = new datatable();

			//now read the file.
			readsav(filename, ct, ew2, ewmb2, ew5, ew3, lfce, mod, data,
					 ev, ew, ewmb, fce, &vmin);

			//set error status to zero
			 return 0;


		}
		return 22;

	}


}




//The following should not be included for compilations for Windows:
#ifndef _WINDOWS_GUI

//A global function for error reporting
void errmsg(int err,int erri) {

if (err==30) {
	std::cout << "End Reached at traceback #"<<erri<<"\n";
   return;
}
if (err==100) {
	std::cout << "error # "<<erri;
   return;
}
switch (err) {
	case 1:
   	std::cout << "Could not allocate enough memory";
      break;
   case 2:
   	std::cout << "Too many possible base pairs";
      break;
   case 3:
   	std::cout << "Too many helixes in multibranch loop";
   case 4:
   	std::cout << "Too many structures in CT file";
   default:
   	std::cout << "Unknown error";
}
//std::cin >> err;
return;

}

#endif
