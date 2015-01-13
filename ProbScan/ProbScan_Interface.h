/*
 * A program that calculates the probscan function for a strand of nucleic acids.
 * This strand of nucleic acids can be composed of either DNA or RNA.
 *
 * (c) 2009 Mathews Lab, University of Rochester Medical Center.
 * Written by Jessica S. Reuter
 */

#ifndef PROBSCANINTERFACE
#define PROBSCANINTERFACE

#include "../RNA_class/RNA.h"
#include "../RNA_class/ProbScan.h"
#include "../src/ErrorChecker.h"
#include "../src/ParseCommandLine.h"

class probscanInterface {
 public:
	// Public constructor and methods.

	/*
	 * Name:        Constructor.
	 * Description: Initializes all private variables.
	 */
	probscanInterface();

	/*
	 * Name:        parse
	 * Description: Parses command line arguments to determine what options are required for a particular calculation.
	 * Arguments:
	 *     1.   The number of command line arguments.
	 *     2.   The command line arguments themselves.
	 * Returns:
	 *     True if parsing completed without errors, false if not.
	 */
	bool parse( int argc, char** argv );
    void processLine(string input);

	/*
	 * Name:        run
	 * Description: Run calculations.
	 */
	void run();

 private:
	// Private variables.

	// Description of the calculation type.
	string calcType;

    //nuc indices for scan
    int i,j,k,l;
    vector<int>mb;

	// Input and output file names.
	string seqFile;          // The input sequence file.
	string pfsFile;          // The output probscan function save file.
	string inputFile;
    string loop_file;   //input file with multibranch loops
	string constraintFile;   // The constraints file.
	string doubleOffsetFile; // The optional double strand offset file.
	string experimentalFile; // The experimental pair bonus file.
	string SHAPEFile;        // The SHAPE constraints file.

	// The experimental pair bonus offset.
	double experimentalOffset;

	// The experimental pair bonus scaling.
	double experimentalScaling;

	// The intercept for SHAPE constraints.
	double intercept;

	// Flag signifying if calculation handles RNA (true) or DNA (false).
	bool isRNA;
    //Flad signifiying if we're working with muiltibranch loops
    bool multibranch;

	bool fromSequence;

	// The maximum pairing distance.
	int maxDistance;

	// The slope for SHAPE constraints.
	double slope;

	// The temperature at which calculation occurs.
	double temperature;
};

#endif /* PROBSCANINTERFACE */
