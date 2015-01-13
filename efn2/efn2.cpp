/*
 * A program that calculates the free energy of a strand of nucleic acids.
 * This strand of nucleic acids can be composed of either DNA or RNA.
 *
 * (c) 2008 Mathews Lab, University of Rochester Medical Center.
 * Written by Jessica S. Reuter
 */
#include "efn2.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
efn2Interface::efn2Interface() {

	// Initialize the calculation type description.
	calcType = "efn2";

	// Initialize the nucleic acid type.
	isRNA = true;

	// Initialize the SHAPE intercept.
	intercept = -0.6;

	// Initialize the SHAPE slope.
	slope = 1.8;

	// Initialize flag that signifies streaming to standard output.
	stdPrint = false;

	// Initialize the calculation temperature.
	temperature = 310.15;

	// Initialize the flag that signifies writing a thermodynamic details file.
	writeTherm = false;

	// Initialize the simple energy function.
	simple = false;
}

///////////////////////////////////////////////////////////////////////////////
// Parse the command line arguments.
///////////////////////////////////////////////////////////////////////////////
bool efn2Interface::parse( int argc, char** argv ) {

	// Create the command line parser and build in its required parameters.
	ParseCommandLine* parser = new ParseCommandLine( "efn2" );
	parser->addParameterDescription( "ct file", "The name of a file containing structure CT data." );
	parser->addParameterDescription( "output file", "The energy file to which output is written. Depending on the options selected, this may be one of the following file types. 1) Simple list (Lists free energy for each structure, lowest first). 2) Thermodynamic details (Writes details of every substructure in each structure, and the corresponding free energy of each)." );

	// Add the DNA option.
	vector<string> dnaOptions;
	dnaOptions.push_back( "-d" );
	dnaOptions.push_back( "-D" );
	dnaOptions.push_back( "--DNA" );
	parser->addOptionFlagsNoParameters( dnaOptions, "Specify that the sequence is DNA, and DNA parameters are to be used. Default is to use RNA parameters." );

	// Add the simple eneergy function option.
	vector<string> simpleOptions;
	simpleOptions.push_back( "-s" );
	simpleOptions.push_back( "-S" );
	simpleOptions.push_back( "--simple" );
	parser->addOptionFlagsNoParameters( simpleOptions, "Specify the simple energy function for multibranch loops, used by the dynamic programming algorithms (Fold, partition, stochastic, AllSub, etc.), should be used. If this is not specified, an more sophisticated energy function is used, and the energies might not match those estimated for structures during structure prediction." );

	// Add the print option.
	vector<string> printOptions;
	printOptions.push_back( "-p" );
	printOptions.push_back( "-P" );
	printOptions.push_back( "--print" );
	parser->addOptionFlagsNoParameters( printOptions, "Print the simple list file to standard output. This won't override default behavior of writing to a file. Thermodynamic files (if written) are not piped." );

	// Add the SHAPE option.
	vector<string> shapeOptions;
	shapeOptions.push_back( "-sh" );
	shapeOptions.push_back( "-SH" );
	shapeOptions.push_back( "--SHAPE" );
	parser->addOptionFlagsWithParameters( shapeOptions, "Specify a SHAPE constraints file to be applied. These constraints are pseudoenergy constraints. Default is to have no constraints applied." );

	// Add the SHAPE intercept option.
	vector<string> shapeInterceptOptions;
	shapeInterceptOptions.push_back( "-si" );
	shapeInterceptOptions.push_back( "-SI" );
	shapeInterceptOptions.push_back( "--SHAPEintercept" );
	parser->addOptionFlagsWithParameters( shapeInterceptOptions, "Specify an intercept used with SHAPE constraints. Default is -0.6 kcal/mol." );

	// Add the SHAPE slope option.
	vector<string> shapeSlopeOptions;
	shapeSlopeOptions.push_back( "-sm" );
	shapeSlopeOptions.push_back( "-SM" );
	shapeSlopeOptions.push_back( "--SHAPEslope" );
	parser->addOptionFlagsWithParameters( shapeSlopeOptions, "Specify a slope used with SHAPE constraints. Default is 1.8 kcal/mol." );

	// Add the temperature option.
	vector<string> tempOptions;
	tempOptions.push_back( "-t" );
	tempOptions.push_back( "-T" );
	tempOptions.push_back( "--temperature" );
	parser->addOptionFlagsWithParameters( tempOptions, "Specify the temperature at which calculation takes place in Kelvin. Default is 310.15 K, which is 37 degrees C." );

	// Add the details option.
	vector<string> detailsOptions;
	detailsOptions.push_back( "-w" );
	detailsOptions.push_back( "-W" );
	detailsOptions.push_back( "--writedetails" );
	parser->addOptionFlagsNoParameters( detailsOptions, "Write a thermodynamic details file. The thermodynamic details file replaces the list file that is outputted by default." );

	// Parse the command line into pieces.
	parser->parseLine( argc, argv );

	// Get required parameters from the parser.
	if( !parser->isError() ) {
		ctFile = parser->getParameter( 1 );
		outFile = parser->getParameter( 2 );
	}

	// Get the DNA option.
	if( !parser->isError() ) { isRNA = !parser->contains( dnaOptions ); }

	// Get the simple energy rule option.
	if (!parser->isError() ) { simple = parser->contains( simpleOptions); }

	// Get the print option.
	if( !parser->isError() ) { stdPrint = parser->contains( printOptions ); }

	// Get the SHAPE data and options.
	if( !parser->isError() ) {
		SHAPEFile = parser->getOptionString( shapeOptions );
		if( !parser->isError() ) { parser->setOptionDouble( shapeInterceptOptions, intercept ); }
		if( !parser->isError() ) { parser->setOptionDouble( shapeSlopeOptions, slope ); }
	}

	// Get the temperature option.
	if( !parser->isError() ) {
		parser->setOptionDouble( tempOptions, temperature );
		if( temperature < 0 ) { parser->setError( "temperature" ); }
	}

	// Get the write thermodynamic details option.
	if( !parser->isError() ) { writeTherm = parser->contains( detailsOptions ); }

	// Delete the parser and return whether the parser encountered an error.
	bool noError = ( parser->isError() == false );
	delete parser;
	return noError;
}

///////////////////////////////////////////////////////////////////////////////
// Run calculations.
///////////////////////////////////////////////////////////////////////////////
void efn2Interface::run() {

	// Create a variable that handles errors.
	int error = 0;

	/*
	 * Use the constructor for RNA that specifies a filename.
	 * Specify type = 1 (CT file).
	 * isRNA identifies whether the strand is RNA (true) or DNA (false).
	 *
	 * After construction of the strand data structure, create the error checker which monitors for errors.  
	 * Throughout, the error status of the calculation is checked with a variant of the isErrorStatus method, which returns 0 if no error occurred.
	 * The calculation proceeds as long as error = 0.
	 */
	cout << "Initializing nucleic acids..." << flush;
	RNA* strand = new RNA( ctFile.c_str(), 1, isRNA );
	ErrorChecker<RNA>* checker = new ErrorChecker<RNA>( strand );
	error = checker->isErrorStatus();
	if( error == 0 ) { cout << "done." << endl; }

	/*
	 * Check the strand for pseudoknots by looking through each structure with the ContainsPseudoknot method.
	 * efn2 cannot handle pseudoknots, so if the strand contains pseudoknots, this is considered an error.
	 */
	if( error == 0 ) {
		structures = strand->GetStructureNumber();
		for( int i = 1; i <= structures; i++ ) {
			if( strand->ContainsPseudoknot( i ) ) {
				cerr << "Nucleic acids contain pseudoknots; cannot proceed." << endl;
				error = -1;
			}
		}
	}

	/*
	 * Set the temperature using the SetTemperature method.
	 * Only set the temperature if a given temperature doesn't equal the default.
	 * If the temperature does need to be set, use the error checker's isErrorStatus method to check for errors.
	 */
	if( ( error == 0 ) && ( temperature != 310.15 ) ) {

		// Show a message saying that the temperature is being set.
		cout << "Setting temperature..." << flush;

		// Set the temperature and check for errors.
		int tempError = strand->SetTemperature( temperature );
		error = checker->isErrorStatus( tempError );

		// If no error occurred, print a message saying that temperature is set.
		if( error == 0 ) { cout << "done." << endl; }
	}

	/*
	 * Read SHAPE constraints, if applicable.
	 * When reading SHAPE constraints, use the ReadSHAPE method.
	 * After constraints are read, use the error checker's isErrorStatus method to check for errors.
	 */
	if( error == 0 && SHAPEFile != "" ) {

		// Show a message saying that SHAPE constraints are being read.
		cout << "Applying SHAPE constraints..." << flush;

		// Initialize the single stranded SHAPE slope and intercept.
		// For now, these are hard-coded as 0.
		double slopeSingle = 0;
		double interceptSingle = 0;

		// Read SHAPE constraints and check for errors.
		int constraintError = strand->ReadSHAPE( SHAPEFile.c_str(), slope, intercept, slopeSingle, interceptSingle );
		error = checker->isErrorStatus( constraintError );

		// If no error occurred, print a message saying that SHAPE constraints are set.
		if( error == 0 ) { cout << "done." << endl; }
	}

	/*
	 * Do the efn2 calculation.
	 * If the user wants a simple output file, get free energies for each structure using the CalculateFreeEnergy method.
	 * If the user wants a thermodynamic details file, write the file with the WriteThemodynamicDetails method.
	 */
	if( error == 0 ) {

		// Write a thermodynamic details file, if asked.
		if( writeTherm ) {

			// Show a message saying that the details file is being written.
			cout << "Writing thermodynamic details file..." << flush;

			// Write the thermodynamic details file and check for errors.
			int thermError = strand->WriteThermodynamicDetails( outFile.c_str() );
			error = checker->isErrorStatus( thermError );

			// Print a message saying that the details file has been written.
			if( error == 0 ) { cout << "done." << endl; }
		}

		// Write a simple list file, if asked.
		else {

			// Show a message saying that the list file is being written.
			cout << "Writing free energy list file..." << flush;

			// For each structure, calculate its energy and push it into the energies vector.
			// Changed this 5-14-2013, no longer pushes into a vector. Instead use calculatefreeenergy method from RNA class.
			// This is so it works with openMP -- vector would be filled out of order when loop runs in parallel and the indexing would be messed up.
			// If an error occurs, set it.
			//vector<double> energies;
			strand->CalculateFreeEnergy( 1 , simple );			//calculate the free energy of the first structure first (this is to make it play well with openMP)			
			#ifdef SMP
			#pragma omp parallel for
			#endif
			for( int i = 2; i <= structures; i++ ) { 
				strand->CalculateFreeEnergy( i , simple ); 		//now calculate the free energy of the remaining structures
				error = checker->isErrorStatus();
				//if( error == 0 ) { energies.push_back( energy ); }
				//else break;
			}

			// If all free energies were calculated correctly, write the output file.
			if( error == 0 ) {
				ofstream out( outFile.c_str() );
				for( int i = 1; i <= structures; i++ ) {
					int index = i - 1;
					out << "Structure: " << i << "   Energy = " << fixed << setprecision( 1 ) << strand->GetFreeEnergy(i)  << endl; //changed to GetFreeEnergy instead of writing the energies vector
				}
				out.close();
			}

			// Print a message saying that the list file has been written.
			if( error == 0 ) { cout << "done." << endl; }

			// If the output should be piped to standard output, then pipe it.
                        if( ( error == 0 ) && ( stdPrint ) ) {
                          cout << endl << "Generated output file: " << outFile << endl << endl;
                          for( int i = 1; i <= structures; i++ ) {
	                          cout << "Structure: " << i << "   Energy = " << fixed << setprecision( 1 ) << strand->GetFreeEnergy(i) << endl;
			  }
                          cout << endl;
                        }
		}
	}

	// Delete the error checker and data structure.
	delete checker;
	delete strand;

	// Print confirmation of run finishing.
	if( error == 0 ) { cout << calcType << " complete." << endl; }
	else { cerr << calcType << " complete with errors." << endl; }
}

///////////////////////////////////////////////////////////////////////////////
// Main method to run the program.
///////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] ) {

	efn2Interface* runner = new efn2Interface();
	bool parseable = runner->parse( argc, argv );
	if( parseable == true ) { runner->run(); }
	delete runner;
	return 0;
}

