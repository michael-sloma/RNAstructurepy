#ifndef _TURBOFOLD_OBJECT_
#define _TURBOFOLD_OBJECT_

#include <vector>
#include <string>
using namespace std;

//Include the correct flavor of TProgressDialog.h
//This is used to provide progress information back to the interface.
#ifdef _WINDOWS_GUI
#include "../RNAstructure_windows_interface/TProgressDialog.h"
#else
#ifdef _JAVA_GUI
#include "../RNAstructure_java_interface/SWIG/TProgressDialog.h"
#else
#include "../src/TProgressDialog.h"
#endif
#endif

class RNA;
class t_matrix;
struct t_aln_env_result;
class t_structure;
class t_ansi_mutex;
class t_turbofold_thread;

enum{NO_ERRORS,
CONSTRUCTOR_ERROR, 
ISEQ_ARGUMENT_OVERFLOW_ERROR,
THREAD_SCHEDULE_ERROR,
UNFINISHED_THREAD_ERROR,
RUN_REFOLDING_THREAD_ERROR,
ALIGNMENT_INFO_COMPUTATION_ERROR,
RNALIB_PROBKNOT_ERROR, 
RNALIB_THRESHOLDING_ERROR, 
RNALIB_MEA_ERROR,
RNALIB_PARTITIONFUNCTION_ERROR,
RNALIB_GETPAIR_ERROR,
RNALIB_WRITECT_ERROR,
RNALIB_GETPAIRPROBABILITY_ERROR,
RNALIB_READSHAPE_ERROR,
RNALIB_SETTEMPERATURE_ERROR,
RNA_LIB_RNA_CONSTRUCTOR_ERROR,
RNALIB_SETDISTANCE_ERROR,
N_TF_ERRORS};

static char err_strings[N_TF_ERRORS+1][100] = { "No error (Cross fingers).",//NO_ERRORS,
"Constructor Error.", //CONSTRUCTOR_ERROR,
"Sequence index argument error.", //ISEQ_ARGUMENT_OVERFLOW_ERROR,
"Thread schedule error.", //THREAD_SCHEDULE_ERROR,
"Child thread could not finish folding.", //UNFINISHED_THREAD_ERROR,
"Could not start a refolding thread.", //RUN_REFOLDING_THREAD_ERROR,
"Alignment information computation failed.", //ALIGNMENT_INFO_COMPUTATION_ERROR,
"ProbKnot error.", //RNALIB_PROBKNOT_ERROR,
"Thresholding error.", //RNALIB_THRESHOLDING_ERROR,
"MEA computation error.", //RNALIB_MEA_ERROR,
"Partition function error.", //RNALIB_PARTITIONFUNCTION_ERROR,
"GetPair error.", //RNALIB_GETPAIR_ERROR,
"WriteCt error.", //RNALIB_WRITECT_ERROR,
"GetPairProbability error.", //RNALIB_GETPAIRPROBABILITY_ERROR,
"ReadSHAPE error.", //RNALIB_READSHAPE_ERROR,
"SetTemperature error.", //RNALIB_SETTEMPERATURE_ERROR,
"RNA class constructor failed." //RNA_LIB_RNA_CONSTRUCTOR_ERROR,
"Error in setting maximum distance at level of an RNA class."//RNALIB_SETDISTANCE_ERROR
}; //N_TF_ERRORS


//! TurboFold Class.
/*!
	The TurboFold class provides an entry point for the TurboFold algorithm in RNAstructure.
*/

//Note the stylized comments provide facility for automatic documentation via doxygen.

class TurboFold
{
public:

	//******************************
	//Constructors:
	//******************************

	//!Constructor - user provides a filename for a FASTA file.

	//!	Input file should contain two or more FASTA sequences.
	//! The constructor reads parameter files from disk.  The location should be specified in the DATAPATH environment variable.
	//! If DATAPATH is undefined, the program will attempt to load the files from the present working directory.
	//!	This constructor generates internal error codes that can be accessed by GetErrorCode() after the constructor is called.  0 = no error.
	//! The errorcode can be resolved to a c string using GetErrorMessage.
	//!	\param fasta_fp is a NULL terminated c string that give a filename.  This must be 1000 or fewer characters.
	TurboFold(const char fasta_fp[]);

	//!Constructor - user provides a vector array of strings that provide input and output file names.

	//!	The output files are partition function save files, which can be read by the RNA class to determine ppair probabilities.
	//! There needs to be one output file per input sequence.
	//! The constructor reads parameter files from disk.  The location should be specified in the DATAPATH environment variable.
	//! If DATAPATH is undefined, the program will attempt to load the files from the present working directory.
	//!	This constructor generates internal error codes that can be accessed by GetErrorCode() after the constructor is called.  0 = no error.
	//! The errorcode can be resolved to a c string using GetErrorMessage.
	//!	\param sequences is a vector of strings that provide sequence file names.  These files need to be either FASTA or .seq.
	//! \param saves is a vector of strings that provide file names for output partition function save files.  This array needs to have exactly the same number of elements as sequences.
	TurboFold(vector<string>* sequences, vector<string>* saves);	
	
	//******************************
	//Destructor:
	//******************************

	//!Destructor
	~TurboFold();

	//******************************
	// Internal error code for this TurboFold object:
	//******************************
	
	//! Get an integer that reports the current error status of the class.

	//! Functions generate internal errors that can be accessed using this function.
	//! An error code of zero is no error.
	//! A non-zero error code can be resolved to a cstring or string using GetErrorMessage() or GetErrorString().
	//! \return An integer that indicates error status.
	int GetErrorCode();

	//!	Return error messages based on code from GetErrorCode or function-returned error codes.		

	//!\param err_code is the integer error code provided by GetErrorCode().
	//!\return A pointer to a c string that provides an error message.
	char* GetErrorMessage(int err_code);

	//!	Return error messages based on code from GetErrorCode and other error codes.		

	//!\param err_code is the integer error code provided by GetErrorCode() or from other functions that return integer error codes.
	//!\return A string that provides an error message.
	string GetErrorString(int err_code); 

	//******************************
	// Functions to apply restraints:
	//******************************

	//! Set a maximum distance between nucleotides that can pair.

	//! This function must be called before fold() and will limit the distance between nucleotides that can pair.
	//! \param distance is an integer that specifies the maximum distance between nucleotides that can pair, , i.e. |j-i| < distance for nucleotide i to pair to j.
	//!	\return An integer error code that can be resolved to an error message using GetErrorMessage() or GetErrorString().  0 is no error. 
	int SetMaxPairingDistance(int distance);

	//! Read and apply SHAPE mapping data to a specific sequence.

	//! The pseudofree energy approach will be used to apply SHAPE data to restrain structure prediction.  Where DG(per stack with ith nucleotide) = slope (SHAPE on ith nucleotide) + intercept.
	//! This function must be called before fold().
	//! \param i_seq is the sequence number to which the restraint should be applied, where the number starts at 1.
	//! \param fp[] is a cstring that provides the name of the file that contains the normalized SHAPE mapping data.
	//! \param par1 is the slope in kcal/mol.
	//! \param par2 is the intercept in kcal/mol.
	//!	\return An integer error code that can be resolved to an error message using GetErrorMessage() or GetErrorString().  0 is no error. 
	int ReadSHAPE(const int i_seq, const char fp[], const double par1, const double par2);

	//******************************
	// Functions to set the temperature:
	//******************************

	//! Set the folding temperature.

	//! This function must be called before fold(). 
	//! If this function is not called, the default temperature of 310.15 K (37 degrees C) is used.
	//! \param temp is the temperature in Kelvin.
	//!	\return An integer error code that can be resolved to an error message using GetErrorMessage() or GetErrorString().  0 is no error. 
	int SetTemperature(double temp);

	//******************************
	// Main TurboFold algorithm:
	//******************************

	//! The main TurboFold algorithm.

	//! This function accomplishes the task of determining the pair probabilities.
	//! This function must be called before any of the structure prediction methods can be used.
	//! \param gamma is the weight of the extrinsic information.  Larger gamma will result in more consistent structures.  The default is 0.3 and this provided a good structure prediction accuracy in benchmarks.
	//! \param n_iterations is the number of iterations that should be performed to converge the base pairing probabilities.  The default is 3 because benchmarks showed only marginal improvement with further iterations.
	//! \param _n_parallel_pfunctions is the number of threads to use.  For code compiled in serial, this must be 1, which is the default.  Define COMPILE_SMP to build for multithreading.
	//!	\return An integer error code that can be resolved to an error message using GetErrorMessage() or GetErrorString().  0 is no error. 
	int fold(double gamma = 0.3, int n_iterations = 3, int _n_parallel_pfunctions = 1); 
	
	//******************************
	// Functions that predict structures:
	//****************************** 

	//! Use the ProbKnot algorithm to predict a structure for a sequence.

	//! This function can predict pseudoknots.
	//! This function can only be called after fold() is called.
	//! \param i_seq is the sequence number, where the number starts at 1.
	//! \param n_iterations is the number of ProbKnot iterations.
	//! \param minhelixlength is the length of the shortest helix allowed.
	//!	\return An integer error code that can be resolved to an error message using GetErrorMessage() or GetErrorString().  0 is no error.
	int ProbKnot(const int i_seq, const int n_iterations, const int minhelixlength);


	//! Predict a structure for a sequence that is composed of highly probably pairs.

	//! This function can only be called after fold() is called.
	//! \param i_seq is the sequence number, where the number starts at 1.
	//! \param probability is the pairing probability threshold, where pairs will be predicted if they have a higher probability.  Note that a value of less than 0.5 (50%), will cause an error.  The default value of zero will trigger the creation of 8 structures, with thresholds of >=0.99, >=0.97, >=0.95, >=0.90, >=0.80, >=0.70, >=0.60, >0.50.
	//!	\return An integer error code that can be resolved to an error message using GetErrorMessage() or GetErrorString().  0 is no error.
	int PredictProbablePairs(const int i_seq, const float probability);


	//! Predict maximum expected accuracy structures for a sequence.

	//! This function can only be called after fold() is called.
	//! The expectd accuracy score for a structure is = gamma * 2 * (sum of pairing probabilities for pairs) + (sum of unpairing probabilities for single stranded nucleotides).
	//! \param i_seq is the sequence number, where the number starts at 1.
	//! \param maxPercent is the maximum difference in score allowed for generation of suboptimal structures.
	//! \param maxStructures is the maximum number of suboptimal structures allowed.
	//! \param window is the window parameter that controls what suboptimal structures can be included.  0 is the minimum and the higher the window, the more different suboptimal structures must be from each other.
	//! \param gamma is the weight on base pairs.  The default of 1.0 works well based on benchmarks on single sequence calculations.
	//!	\return An integer error code that can be resolved to an error message using GetErrorMessage() or GetErrorString().  0 is no error.
	int MaximizeExpectedAccuracy(const int i_seq, const double maxPercent, const int maxStructures, const int window, const double gamma=1.0);

	//******************************
	// Functions that provide structure information:
	//****************************** 
	
	//! Provide pairing information.

	//! This function can only be called after one of the structure prediction methods is a called.
	//!	This function generates internal error codes that can be accessed by GetErrorCode() after the constructor is called.  0 = no error.
	//! The errorcode can be resolved to a c string using GetErrorMessage.
	//! \param i_seq is the sequence number, where the number starts at 1.
	//! \param i is the nucleotide.
	//! \param structurenumber is the structure number.  This can be used to specify a suboptimal structure, but defaults to 1.
	//!	\return The nucleotide to which i is paired in sequence i_seq and suboptimal structure number structurenumber.  Zero indicates that the nucleotide is unpaired.
	int GetPair(const int i_seq, const int i, const int structurenumber=1);
	
	//! Provide pairing probability information.

	//! This function can only be called after fold() is a called.
	//!	This function generates internal error codes that can be accessed by GetErrorCode() after the constructor is called.  0 = no error.
	//! The errorcode can be resolved to a c string using GetErrorMessage.
	//! \param i_seq is the sequence number, where the number starts at 1.
	//! \param i is the 5' nucleotide in a pair.
	//! \param j is the 3' nucleotide in a pair.
	//!	\return The probability that i is paired to j in sequence i_seq.
	double GetPairProbability(const int i_seq, const int i, const int j);
	

	//******************************
	// Function that provides information about the number of sequences:
	//****************************** 

	//! Provide the number of sequences used in the calculation.

	//! \return The number of sequences used in the calculation.
	int GetNumberSequences();


	//******************************
	// Function that writes structure output:
	//****************************** 

	//! Write the predicted structures for a specific sequence to a ct file.

	//! This function can only be called after one of the structure prediction methods is a called.
	//! \param i_seq is the sequence number, where the number starts at 1.
	//! \param fp is a cstring that gives the filename to which the ct table is to be written.
	//!	\return An integer error code that can be resolved to an error message using GetErrorMessage() or GetErrorString().  0 is no error.
	int WriteCt(const int i_seq, const char fp[]);


	//*****************************
	// Functions for providing progress information
	//*****************************
	//!Provide a TProgressDialog for following calculation progress.
	//!A TProgressDialog class has a public function void update(int percent) that indicates the progress of a long calculation.
	//!\param Progress is a TProgressDialog class.
	void SetProgress(TProgressDialog& Progress);


	//!Provide a means to stop using a TProgressDialog.
	//!StopProgress tells the RNA class to no longer follow progress.  This should be called if the TProgressDialog is deleted, so that this class does not make reference to it.
	void StopProgress();

	

	//Multithreading support variables and functions
	t_turbofold_thread** refolding_threads;
	t_ansi_mutex* turbofold_threads_mutex;
	int n_parallel_pfunctions;
	int n_iterations;
	int start_next_refolding_thread(bool is_main_thread, int _i_iter);

	vector<t_structure*>* sequences;
	vector<RNA*>* folders; // These are the folders, although they are RNA objects, their main job is to  
	vector<char*>* saves;
	

private:

	int generate_extrinsic_information(int i_iter, const double gamma); 
	int initialize_alignment_information();


	int err_code;//Internal error code

	TProgressDialog *progress;//Used to provide progress information to the interface

	

	double** similarities; // Similarities between the sequences, between 0 and 1.
	double**** aln_mapping_probs; // Basically, the co-incidence probabilities. 
	vector<t_matrix*>* extrinsic_info_list; // Current set of extrinsic for each sequence.

	t_aln_env_result*** aln_env_results;
};

#endif // _TURBOFOLD_OBJECT_


