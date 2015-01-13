/*
 * A program that predicts structures using the TurboFold algorithm.
 *
 * (c) 2010 Mathews Lab, University of Rochester Medical Center.
 * Written by Jessica S. Reuter
 */

#include "TurboFold_Interface.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
TurboFold_Interface::TurboFold_Interface() {

	// Initialize the calculation type description.
	calcType = "TurboFold";

	// Initialize default number of sequences and number of processors to 1.
	seqNumber = 1;
	processors = 1;

	// Initialize the default mode.
	mode = "MEA";

	// Initialize general options, which are applied independent of mode.
	{
		// Initialize the maximum pairing distance to 0, which indicates no limit.
		distance = 0;

		// Initialize the SHAPE intercept.
		intercept = -0.6;

		// Initialize the SHAPE slope.
		slope = 1.8;

		// Initialize the default temperature.
		temperature = 310.15;

		// Initialize the TurboFold gamma.
		turboGamma = 0.3;

		// Initialize the TurboFold iterations.
		turboIterations = 3;
	}

	// Initialize maximum expected accuracy mode parameters.
	{
		// Initialize the maximum number of structures.
		maxStructures = 1000;

		// Initialize the maximum expected accuracy mode gamma.
		meaGamma = 1;

		// Initialize the maximum percent energy difference.
		percent = 50;

		// Initialize the window size.
		windowSize = 5;
	}

	// Initialize ProbKnot mode parameters.
	{
		// Initialize the number of ProbKnot iterations.
		pkIterations = 1;

		// Initialize the minimum helix length for a pseudoknot.
		minHelixLength = 3;
	}

	// Initialize Threshold mode parameters.
	{
		// Initialize the default probability threshold.
		// The default of 0 means that structures are generated at multiple thresholds.
		threshold = 0;
	}
}

///////////////////////////////////////////////////////////////////////////////
// Parse the command line arguments.
///////////////////////////////////////////////////////////////////////////////
bool TurboFold_Interface::parse( int argc, char** argv ) {

	// Determine the proper executable name, depending on if the executable is serial or SMP.
#ifndef COMPILE_SMP
	string type = "TurboFold";
#else
	string type = "TurboFold-smp";
#endif

	// Create the command line parser and build in its required parameters.
	ParseCommandLine* parser = new ParseCommandLine( type );
	parser->addParameterDescription( "configuration file", "The name of a file containing configuration data." );

	// Tell the parser that a special usage message is needed for this interface.
	parser->setSpecializedUsage();

	// Parse the command line into pieces.
	parser->parseLine( argc, argv );

	// If the specialized usage has been asked for, print it out.
	// An error is set here only for the purpose of preventing further parsing and calculation; it's not a real error.
	if( parser->isSpecializedUsage() ) {
		usage( parser );
		parser->setError();
		return false;
	}

	// Open the config file.
	// If the file isn't valid, delete the parser and return false.
	ConfigFile file( parser->getParameter( 1 ) );
	if( file.isValid() == false ) {
		delete parser;
		return false;
	}

	// Check if the config file contains the Mode flag, which is always necessary for reading data.
	// If it exists, read it in.
	// If the flag isn't present, set an error, delete the parser, and return false.
	bool isReadable = file.contains( "Mode" );
	if( isReadable ) { mode = file.getOption<string>( "Mode" ); }
	else {
		parser->setErrorSpecialized( "Configuration file must contain the Mode flag." );
		delete parser;
		return false;
	}

	// If SMP calculations are being done, make sure the processors flag has been specified.
	// If it has, set the number of processors if possible.
	// If it hasn't, show an error, delete the parser, and return false.
#ifdef COMPILE_SMP
	if( file.contains( "Processors" ) ) {
		processors = file.getOption<int>( "Processors" );
		if( processors < 0 ) {
			parser->setError( "number of processors" );
			delete parser;
			return false;
		}
	} else {
		parser->setErrorSpecialized( "Configuration file must contain the Processors flag." );
		delete parser;
		return false;
	}
#endif

	// Initialize the files vector.
	// Index 0: Sequence files
	// Index 1: CT files
	// Index 2: Optional partition function file names
	// Index 3: Optional SHAPE constraint files
	for( int i = 1; i <= 4; i++ ) {
		vector<string> row;
		files.push_back( row );
	}

	// Check whether sequence files can be read singly.
	// If they can, attempt to read them into the 2D file vector.
	// An error occurs doing this, set an error, delete the parser, and return false.
	bool isReadableSingly =
		file.contains( "Seq1" ) &&
		file.contains( "CT1" ) &&
		file.contains( "SequenceNumber" );
	if( isReadableSingly ) {

		// Set the sequence number.
		// If the sequence number is less than or equal to 0, set an error, delete the parser, and return false.
		seqNumber = file.getOption<int>( "SequenceNumber" );
		if( seqNumber <= 0 ) {
			parser->setErrorSpecialized( "Sequence number must be greater than 0." );
			delete parser;
			return false;
		}

		// Read in all required files.
		// For each possible sequence number, read in the corresponding data file.
		for( int i = 1; i <= seqNumber; i++ ) {

			// Get the next sequence number as a string.
			stringstream numStream( stringstream::in | stringstream::out );
			numStream << i;
			string num = numStream.str();

			// Read in the next sequence file.
			// If an error occurs, set an error, delete the parser, and return false.
			string seqString = "Seq" + num;
			if( file.contains( seqString ) ) { files[0].push_back( file.getOption<string>( seqString ) ); }
			else {
				parser->setErrorSpecialized( "The number of sequence files specified must be equal to SequenceNumber." );
				delete parser;
				return false;
			}

			// Read in the next CT file.
			// If an error occurs, set an error, delete the parser, and return false.
			string ctString = "CT" + num;
			if( file.contains( ctString ) ) { files[1].push_back( file.getOption<string>( ctString ) ); }
			else {
				parser->setErrorSpecialized( "The number of CT files specified must be equal to SequenceNumber." );
				delete parser;
				return false;
			}
		}
	}
	// Check whether sequence files can be read in groups.
	// If they can, attempt to read them into the 2D file vector.
	// If an error occurs doing this, set an error, delete the parser, and return false.
	bool isReadableGroups =
		file.contains( "InSeq" ) &&
		file.contains( "OutCT" );
	if( isReadableGroups ) {

		// Get the sequences.
		string seqData = file.getOption<string>( "InSeq" );
		unsigned int seqLast = seqData.length() - 1;
		if( seqData[0] == '{' && seqData[seqLast] == '}' ) {
			seqData = seqData.erase( 0, 1 );
			seqData = seqData.erase( seqLast - 1, 1 );
			seqLast = seqData.length() - 1;
			if( seqData[seqLast] == ';' ) { seqData = seqData.erase( seqLast, 1 ); }
			stringstream seqStr( seqData );
			string seqFile;
			while( seqStr.good() ) {
				getline( seqStr, seqFile, ';' );
				if( seqFile != "" ) { files[0].push_back( seqFile ); }
			}
		} else {
			if( seqData[0] != '{' ) { parser->setErrorSpecialized( "Sequence group has no start bracket." ); }
			else { parser->setErrorSpecialized( "Sequence group has no end bracket." ); }
			delete parser;
			return false;
		}

		// Get the CT files.
		string ctData = file.getOption<string>( "OutCT" );
		unsigned int ctLast = ctData.length() - 1;
		if( ctData[0] == '{' && ctData[ctLast] == '}' ) {
			ctData = ctData.erase( 0, 1 );
			ctData = ctData.erase( ctLast - 1, 1 );
			ctLast = ctData.length() - 1;
			if( ctData[ctLast] == ';' ) { ctData = ctData.erase( ctLast, 1 ); }
			stringstream ctStr( ctData );
			string ctFile;
			while( ctStr.good() ) {
				getline( ctStr, ctFile, ';' );
				if( ctFile != "" ) { files[1].push_back( ctFile ); }
			}
		} else {
			if( ctData[0] != '{' ) { parser->setErrorSpecialized( "Sequence group has no start bracket." ); }
			else { parser->setErrorSpecialized( "Sequence group has no end bracket." ); }
			delete parser;
			return false;
		}

		// If the number of sequence files equals the number of CT files, set that number as the number of sequences.
		// If not, set an error, delete the parser, and return.
		if( files[0].size() == files[1].size() ) { seqNumber = files[0].size(); }
		else {
			parser->setErrorSpecialized( "Number of sequence files does not equal number of CT files." );
			delete parser;
			return false;
		}
	}

	// If the configuration file isn't set up properly to read sequence files singly or in groups, set an error, delete the parser, and return false.
	if( ( isReadableSingly == false ) && ( isReadableGroups == false ) ) {
		parser->setErrorSpecialized( "File groups are not specified properly. Please run \"" + type + " -h\" for help." );
		delete parser;
		return false;
	}

	// Check to see if optional save files or SHAPE files were specified.
	stringstream optionalStream( stringstream::in | stringstream::out );
	string saveString, shapeString;
	hasSHAPE = false;
	for( int i = 1; i <= seqNumber; i++ ) {

		// If save file i is specified, get the name.
		optionalStream << "Save" << i;
		saveString = optionalStream.str();
		if( file.contains( saveString ) ) { files[2].push_back( file.getOption<string>( saveString ) ); }
		else { files[2].push_back( "" ); }
		optionalStream.str( "" );

		// If SHAPE file i is specified, get the name.
		optionalStream << "SHAPE" << i;
		shapeString = optionalStream.str();
		if( file.contains( shapeString ) ) {
			files[3].push_back( file.getOption<string>( shapeString ) );
			hasSHAPE = true;
		} else { files[3].push_back( "" ); }
		optionalStream.str( "" );
	}

	// Set the general options.
	if( true ) {

		// Get the TurboFold gamma.
		if( !parser->isError() ) {
			if( file.contains( "Gamma" ) ) {
				turboGamma = file.getOption<double>( "Gamma" );
				if( turboGamma < 0.0 ) { parser->setError( "TurboFold gamma" ); }
			}
		}

		// Get the TurboFold iterations.
		if( !parser->isError() ) {
			if( file.contains( "Iterations" ) ) {
				turboIterations = file.getOption<int>( "Iterations" );
				if( turboIterations < 0 ) { parser->setError( "TurboFold iterations" ); }
			}
		}

		// Get the maximum pairing distance.
		if( !parser->isError() ) {
			if( file.contains( "MaximumPairingDistance" ) ) {
				distance = file.getOption<int>( "MaximumPairingDistance" );
				if( distance < 0 ) { parser->setError( "maximum pairing distance" ); }
			}
		}

		// Get the SHAPE intercept.
		if( !parser->isError() ) {
			if( hasSHAPE ) {
				if( file.contains( "SHAPEintercept" ) ) { intercept = file.getOption<double>( "SHAPEintercept" ); }
			}
		}

		// Get the SHAPE slope.
		if( !parser->isError() ) {
			if( hasSHAPE ) {
				if( file.contains( "SHAPEslope" ) ) { slope = file.getOption<double>( "SHAPEslope" ); }
			}
		}

		// Get the temperature.
		if( !parser->isError() ) {
			if( file.contains( "Temperature" ) ) {
				temperature = file.getOption<double>( "Temperature" );
				if( temperature < 0.0 ) { parser->setError( "temperature" ); }
			}
		}
	}

	// Set the MEA mode options, if applicable.
	if( mode == "MEA" ) {

		// Get the maximum percent energy difference.
		if( !parser->isError() ) {
			if( file.contains( "MaxPercent" ) ) {
				percent = file.getOption<double>( "MaxPercent" );
				if( percent < 0.0 ) { parser->setError( "maximum percent energy difference" ); }
			}
		}

		// Get the maximum number of structures.
		if( !parser->isError() ) {
			if( file.contains( "MaxStructures" ) ) {
				maxStructures = file.getOption<int>( "MaxStructures" );
				if( maxStructures < 0 ) { parser->setError( "maximum number of structures" ); }
			}
		}

                // Get the MEA gamma.
                if( !parser->isError() ) {
		  if( file.contains( "MeaGamma" ) ) { meaGamma = file.getOption<double>( "MeaGamma" ); }
                }

		// Get the window size.
		if( !parser->isError() ) {
			if( file.contains( "Window" ) ) {
				windowSize = file.getOption<int>( "Window" );
				if( windowSize < 0 ) { parser->setError( "window size" ); }
			}
		}
	}

	// Set the ProbKnot mode options, if applicable.
	else if( mode == "ProbKnot" ) {

		// Get the number of ProbKnot iterations.
		if( !parser->isError() ) {
			if( file.contains( "PkIterations" ) ) {
				pkIterations = file.getOption<int>( "PkIterations" );
				if( pkIterations < 0 ) { parser->setError( "ProbKnot iterations" ); }
			}
		}

		// Get the minimum helix length, if applicable.
		if( !parser->isError() ) {
			if( file.contains( "MinHelixLength" ) ) {
				minHelixLength = file.getOption<int>( "MinHelixLength" );
				if( minHelixLength < 0 ) { parser->setError( "minimum helix length" ); }
			}
		}
	}

	// Set the Threshold mode options, if applicable.
	else if( mode == "Threshold" ) {

		// Get the threshold for pairs.
		if( !parser->isError() ) {
			if( file.contains( "Threshold" ) ) {
				threshold = file.getOption<double>( "Threshold" );
				if( threshold < 0.0 ) { parser->setError( "pairing threshold" ); }
			}
		}
	}

	// If an invalid mode was specified, show an error message.
	else {
		parser->setErrorSpecialized( "Invalid mode given; mode must be 'MEA', 'ProbKnot', or 'TurboFold'." );
	}

	// Delete the parser and return whether the parser encountered an error.
	bool noError = ( parser->isError() == false );
	delete parser;
	return noError;
}

///////////////////////////////////////////////////////////////////////////////
// Run calculations.
///////////////////////////////////////////////////////////////////////////////
void TurboFold_Interface::run() {

	// Create a variable that handles errors.
	int error = 0;

	/*
	 * Use the constructor for TurboFold_object that specifies vectors of file names.
	 * This allows for many varied sequence files to be used as input.
	 *
	 * After construction of the data structure, create the error checker which monitors for errors.
	 * Throughout, the error status of the calculation is checked with a variant of the isErrorStatus method, which returns 0 if no error occurred.
	 * The calculation proceeds as long as error = 0.
	 */
	cout << "Initializing nucleic acids..." << flush;
	TurboFold* object = new TurboFold( &files[0], &files[2] );
	ErrorChecker<TurboFold>* checker = new ErrorChecker<TurboFold>( object );
	error = checker->isErrorStatus();
	if( error == 0 ) { cout << "done." << endl; }

	/*
	 * Set the temperature using the SetTemperature method.
	 * Only set the temperature if a given temperature doesn't equal the default.
	 * If the temperature does need to be set, use the error checker's isErrorStatus method to check for errors.
	 */
	if( ( error == 0 ) && ( temperature != 310.15 ) ) {

		// Show a message saying that the temperature is being set.
		cout << "Setting temperature..." << flush;

		// Set the temperature and check for errors.
		int tempError = object->SetTemperature( temperature );
		error = checker->isErrorStatus( tempError );

		// If no error occurred, print a message saying temperature is set.
		if( error == 0 ) { cout << "done." << endl; }
	}

	/*
	 * Set the maximum pairing distance, if required.
	 */
	if( ( error == 0 ) && ( distance > 0 ) ) {

		// Show a message saying that maximum pairing distance is being set.
		cout << "Setting maximum pairing distance..." << flush;

		// Set the maximum pairing distance and check for errors.
		int distError = object->SetMaxPairingDistance( distance );
		error = checker->isErrorStatus( distError );

		// If no error occurred, print a message saying that maximum pairing
		// distance is set.
		if( error == 0 ) { cout << "done." << endl; }
	}

	/*
	 * Read SHAPE data where necessary for specific sequence files.
	 * Since SHAPE data must be read individually for each sequence, use the ReadSHAPE method.
	 */
	if( error == 0 && hasSHAPE ) {

		// Show a message saying SHAPE data is being read.
		cout << "Reading SHAPE data..." << flush;

		// For each sequence, read in SHAPE data, if applicable.
		for( int i = 1; i <= seqNumber; i++ ) {
			if( files[3][i-1] != "" ) {
				int shapeError = object->ReadSHAPE( i, files[3][i-1].c_str(), slope, intercept );
				if( ( error = checker->isErrorStatus( shapeError ) ) == true ) {
					i += seqNumber;
				}
			}
		}

		// Show a message saying SHAPE data reading is done.
		if( error == 0 ) { cout << "done." << endl; }
	}

	/*
	 * Run the TurboFold algorithm using the fold method.
	 * During calculation, monitor progress using the TProgressDialog class and the Start/StopProgress methods of the RNA class.
	 * Neither of these methods require any error checking.
	 * After the main calculation is complete, use the error checker's isErrorStatus method to check for errors.
	 */
	if( error == 0 ) {

		// Show a message saying that the main calculation has started.
		cout << "Folding nucleic acids..." << endl;

		// Create the progress monitor.
		TProgressDialog* progress = new TProgressDialog();
		object->SetProgress( *progress );

		// Run the TurboFold algorithm and check for errors.
		int mainCalcError = object->fold( turboGamma, turboIterations, processors );
		error = checker->isErrorStatus( mainCalcError );

		// Delete the progress monitor.
		object->StopProgress();
		delete progress;

		// If no error occurred, print message that main calculation is done.
		if( error == 0 ) { cout << "done." << endl; }
	}

	/*
	 * Resolve the structures generated by TurboFold using specific methods, depending on the mode selected for TurboFold.
	 * In MEA mode, use the MaximizeExpectedAccuracy method.
	 * In ProbKnot mode, use the ProbKnot method.
	 * In Threshold mode, use the PredictProbablePairs method.
	 * During calculation, monitor progress using the TProgressDialog class and the Start/StopProgress methods of the RNA class.
	 * Neither of these methods require any error checking.
	 * After the resolution calculation is complete, use the error checker's isErrorStatus method to check for errors.
	 */
	if( error == 0 ) {

		// Show a message saying that the resolving calculation has started.
		if( mode == "MEA" ) {
			cout << "Calculating maximum expected accuracy structures..." << flush;
		} else if( mode == "ProbKnot" ) {
			cout << "Predicting pseudoknots..." << flush;
		} else {
			cout << "Calculating probable pairs..." << flush;
		}

		// Resolve the structures and check for errors.
		for( int i = 1; i <= seqNumber; i++ ) {
			int resolveError = 0;
			if( mode == "MEA" ) {
				resolveError = object->MaximizeExpectedAccuracy( i, percent, maxStructures, windowSize, meaGamma );
			} else if( mode == "ProbKnot" ) {
				resolveError = object->ProbKnot( i, pkIterations, minHelixLength );
			} else {
				resolveError = object->PredictProbablePairs( i, threshold );
			}
			error = checker->isErrorStatus( resolveError );

			if( error != 0 ) { i += seqNumber; }
		}

		// If no error occurred, print message that resolving is done.
		if( error == 0 ) { cout << "done." << endl; }
	}

	/*
	 * Write CT output files using the WriteCt method.
	 * After writing is complete, use the error checker's isErrorStatus method to check for errors.
	 */
	if( error == 0 ) {

		// Show a message saying that the CT files are being written.
		cout << "Writing output ct files..." << flush;

		// Write the CT files and check for errors.
		for( int i = 1; i <= seqNumber; i++ ) {
			int writeError = object->WriteCt( i, files[1][i-1].c_str() );
			error = checker->isErrorStatus( writeError );
			if( error != 0 ) { i += seqNumber; }
		}

		// If no errors occurred, show a CT files writing completion message.
		if( error == 0 ) { cout << "done." << endl; }
	}

	// Delete the error checker and data structure.
	delete checker;
	delete object;

	// Print confirmation of run finishing.
	if( error == 0 ) { cout << calcType << " complete." << endl; }
	else { cerr << calcType << " complete with errors." << endl; }
}

///////////////////////////////////////////////////////////////////////////////
// Print out a special usage message for this interface.
///////////////////////////////////////////////////////////////////////////////
void TurboFold_Interface::usage( ParseCommandLine* parser ) {

	cout << "===================================" << endl
	     << "==== Configuration File Format ====" << endl
	     << "===================================" << endl
	     << "Note that configuration options are not case-sensitive." << endl
	     << "Any unrecognized options are ignored." << endl << endl
	     << "Configuration file line format:" << endl
	     << "<Option> = <Value>" << endl << endl;

	cout << "Required input when specifying file groups" << endl
	     << "------------------------------------------" << endl
	     << "The sequence and CT file groups must each specify the same number of files." << endl << endl;
	cout << "InSeq" << endl;
	parser->wrapString( "Flag that can be used to specify a group of sequence files, from Seq1 to Seq<SequenceNumber>. Only one sequence file group can be specified without overwriting files. Group format: {seq1File;seq2File;seq3File;}" );
	cout << "OutCT" << endl;
	parser->wrapString( "Flag that can be used to specify a group of CT files, from CT1 to CT<SequenceNumber>. Only one CT file group can be specified without overwriting files. Group format: {ct1File;ct2File;ct3File;}" );
	cout << "Mode" << endl;
	parser->wrapString( "The mode in which TurboFold is run. A mode can be specified in the following ways. 1) MEA (Maximum expected accuracy). 2) ProbKnot (Pseudoknots). 3) Threshold (Probable pairs)." );

	cout << "Required input when specifying files singly" << endl
	     << "-------------------------------------------" << endl
	     << "Every specified sequence file must have a corresponding CT file specified." << endl << endl;
	cout << "SequenceNumber" << endl;
	parser->wrapString( "The number of sequences given as input." );
	cout << "Seq1 ... Seq<SequenceNumber>" << endl;
	parser->wrapString( "Names of sequence files used as input, from 1 to SequenceNumber." );
	cout << "CT1 ... CT<SequenceNumber>" << endl;
	parser->wrapString( "Names of CT files written as output, from 1 to SequenceNumber." );
	cout << "Mode" << endl;
	parser->wrapString( "The mode in which TurboFold is run. A mode can be specified in the following ways. 1) MEA (Maximum expected accuracy). 2) ProbKnot (Pseudoknots). 3) Threshold (Probable pairs)." );

	cout << "General options" << endl
	     << "---------------" << endl;
	cout << "Gamma" << endl;
	parser->wrapString( "The TurboFold gamma. Default is 0.3." );
	cout << "Iterations" << endl;
	parser->wrapString( "The number of iterations TurboFold goes through. Default is 3 iterations." );
	cout << "MaximumPairingDistance" << endl;
	parser->wrapString( "The maximum distance between nucleotides that can pair. For nucleotide i to pair with j, [i - j| < MaximumPairingDistance. This applies to each sequence. Default is no limit." );
#ifdef COMPILE_SMP
	cout << "Processors" << endl;
	parser->wrapString( "The number of processors on which the calculation runs. Default is 1." );
#endif
	cout << "Save1 ... Save<SequenceNumber>" << endl;
	parser->wrapString( "Names of save files written to by TurboFold, from 1 to SequenceNumber. The number at the end of the flag (1 to SequenceNumber) identifies which sequence the save file will be written for." );
	cout << "SHAPE1 ... SHAPE<SequenceNumber>" << endl;
	parser->wrapString( "Names of SHAPE constraint files. The number at the end of the flag (1 to SequenceNumber) identifies which sequence the constraints will be applied to." );
	cout << "SHAPEintercept" << endl;
	parser->wrapString( "The SHAPE intercept. This value is only used when at least one SHAPE constraint file is specified. Default is 1.8 kcal/mol." );
	cout << "SHAPEslope" << endl;
	parser->wrapString( "The SHAPE slope. This value is only used when at least one SHAPE constraint file is specified. Default is -0.6 kcal/mol." );
	cout << "Temperature" << endl;
	parser->wrapString( "The temperature at which calculations are run, in Kelvin. Default is 310.15 K, which is 37 degrees C." );

	cout << "Maximum expected accuracy (MEA) mode options" << endl
	     << "--------------------------------------------" << endl;
	cout << "MaxPercent" << endl;
	parser->wrapString( "The maximum percent energy difference. Default is 50 percent (Specified as 50, not 0.5)." );
	cout << "MaxStructures" << endl;
	parser->wrapString( "The maximum number of structures. Default is 1000 structures." );
	cout << "MeaGamma" << endl;
	parser->wrapString( "The weight given to pairs. Default is 1.0." );
	cout << "Window" << endl;
	parser->wrapString( "The window size. Default is 5 nucleotides." );

	cout << "Pseudoknot (ProbKnot) mode options" << endl
	     << "----------------------------------" << endl;
	cout << "MinHelixLength" << endl;
	parser->wrapString( "The minimum helix length. Default is 3 nucleotides." );
	cout << "PkIterations" << endl;
	parser->wrapString( "The number of iterations. Default is 1 iteration." );

	cout << "Probable pairs (Threshold) mode options" << endl
	     << "---------------------------------------" << endl;
	cout << "Threshold" << endl;
	parser->wrapString( "The threshold at which pairs should be included in a structure. This should be expressed as a number: 0.5 <= x <= 1.0. Default is 0, which signifies that structures should be generated at multiple thresholds: >= 0.99, >= 0.97, >= 0.95, >= 0.90, >= 0.80, >= 0.70, >= 0.60, and >= 0.50." );
}

///////////////////////////////////////////////////////////////////////////////
// Main method to run the program.
///////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] ) {

	TurboFold_Interface* runner = new TurboFold_Interface();
	bool parseable = runner->parse( argc, argv );
	if( parseable == true ) { runner->run(); }
	delete runner;
	return 0;
}

