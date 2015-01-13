/*
 * An implementation file for a class that holds structure methods for
 * Java drawing.
 *
 * (c) 2012 Mathews Lab, University of Rochester Medical Center.
 * Written by Jessica S. Reuter.
 */

#include "StructureBackend.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
StructureBackend::StructureBackend() { structures = 0; }

///////////////////////////////////////////////////////////////////////////////
// Add probability annotation to all structures.
///////////////////////////////////////////////////////////////////////////////
bool StructureBackend::addAnnotationProbability( string file ) {

	for( int i = 1; i <= structures; i++ ) {
		string result = structureHandlers[i-1].addAnnotationProbability( file );
		if( result != "" ) { return false; }
	}
	return true;
}

///////////////////////////////////////////////////////////////////////////////
// Add SHAPE annotation to all structures.
///////////////////////////////////////////////////////////////////////////////
bool StructureBackend::addAnnotationSHAPE( string file ) {

	for( int i = 1; i <= structures; i++ ) {
		string result = structureHandlers[i-1].addAnnotationSHAPE( file );
		if( result != "" ) { return false; }
	}
	return true;
}

///////////////////////////////////////////////////////////////////////////////
// Flip all structures horizontally.
///////////////////////////////////////////////////////////////////////////////
void StructureBackend::flip() {

	for( int i = 1; i <= structures; i++ ) {
		structureHandlers[i-1].flipHorizontally();
	}
}

///////////////////////////////////////////////////////////////////////////////
// Get the data for a particular structure.
///////////////////////////////////////////////////////////////////////////////
string StructureBackend::getStructureData( int number ) {

	// If the structure handler for this data piece doesn't exist, create it.
	// If an error occurred, return the empty string.
	if( structureHandlers[number-1].toString() == "" ) {
		string result =
			structureHandlers[number-1].readRadial( structureDataFile, number );
		if( result != "" ) { return ""; }
	}

	// Build and return the structure string.
	stringstream stream( stringstream::in | stringstream::out );
	stream << "Structure " << number << " of " << structures << endl
	       << structureHandlers[number-1].toString();
	return stream.str();
}

///////////////////////////////////////////////////////////////////////////////
// Read structure data.
///////////////////////////////////////////////////////////////////////////////
bool StructureBackend::readStructureData( string file ) {

	// Determine how many structures are in the file.
	// If an error occurred while trying to do this, return false.
	RNA* rna = new RNA( file.c_str(), 1 );
	ErrorChecker<RNA>* rnaChecker = new ErrorChecker<RNA>( rna );
	structures = rna->GetStructureNumber();
	bool isConstructorError = ( rnaChecker->returnError() != "" );
	delete rnaChecker;
	delete rna;
	if( isConstructorError ) { return false; }

	// Save the name of the structures data file.
	structureDataFile = file;

	// Reserve space for all possible structures.
	structureHandlers.reserve( structures );

	// If the number of structures is greater than 20, only set the possible
	// number of structures to 20 for now.
	int possibleStructures = structures;
	if( structures > 20 ) { possibleStructures = 20; }

	// Fully read in structures up to the initial reading limit.
	// If an error occurs at any point, return false.
	for( int i = 1; i <= possibleStructures; i++ ) {
		StructureImageHandler handler;
		string result = handler.readRadial( file, i );
		if( result != "" ) { return false; }
		structureHandlers.push_back( handler );
	}

	// If any other structures should exist over the initial limit, initialize
	// them without fully reading their data.
	if( structures > possibleStructures ) {
		for( int i = possibleStructures + 1; i <= structures; i++ ) {
			StructureImageHandler handler;
			structureHandlers.push_back( handler );
		} 
	}

	// Return true; if the method got to the end, data was read correctly.
	return true;
}

///////////////////////////////////////////////////////////////////////////////
// Remove annotation from all structures.
///////////////////////////////////////////////////////////////////////////////
void StructureBackend::removeAnnotation() {
	for( int i = 1; i <= structures; i++ ) {
		structureHandlers[i-1].removeAnnotation();
	}
}

///////////////////////////////////////////////////////////////////////////////
// Set all nucleotides circled or uncircled.
///////////////////////////////////////////////////////////////////////////////
void StructureBackend::setNucleotidesCircled( bool circled ) {
	for( int i = 1; i <= structures; i++ ) {
		structureHandlers[i-1].setNucleotidesCircled( circled );
	}
}

///////////////////////////////////////////////////////////////////////////////
// Write a dot bracket file.
///////////////////////////////////////////////////////////////////////////////
string StructureBackend::writeDotBracketFile( string inFile, string outFile ) {

	// Create the data strand and its error checker.
	RNA* rna = new RNA( inFile.c_str(), 1 );
	ErrorChecker<RNA>* rnaChecker = new ErrorChecker<RNA>( rna );

	// Check for errors.
	string complete = rnaChecker->returnError();

	// If no error occurred, write the file and check for errors afterward.
	if( complete == "" ) {
		int writeError = rna->WriteDotBracket( outFile.c_str() );
		complete = rnaChecker->returnError( writeError );
	}

	// Delete the data strand and error checker.
	delete rna;
	delete rnaChecker;

	// Return the completion string.
	return complete;
}

///////////////////////////////////////////////////////////////////////////////
// Write a helix file.
///////////////////////////////////////////////////////////////////////////////
string StructureBackend::writeHelixFile(
	string inFile, string outFile, int number ) {

	// Create the data strand and its error checker.
	RNA* rna = new RNA( inFile.c_str(), 1 );
	ErrorChecker<RNA>* rnaChecker = new ErrorChecker<RNA>( rna );

	// Check for errors.
	string complete = rnaChecker->returnError();

	// If no error occurred, write the file and check for errors afterward.
	if( complete == "" ) {

		// Open the output stream.
		ofstream out( outFile.c_str() );

		// Determine helices and write them out.
		int length = rna->GetSequenceLength();
		int counter = 1;
		while( counter <= length ) {

			// Get whether the next nucleotide is paired.
			int pair = rna->GetPair( counter, number );

			// If the nucleotide is paired, figure out the extent of the helix
			// and write it out.
			if( pair > counter ) {
				int count = 1;

				while( ( rna->GetPair( counter + 1, number ) ) ==
				       ( rna->GetPair( counter ) - 1 ) ) {
					counter++;
					count++;
				}

				out << ( counter - count + 1 ) << " " << pair << " " << count
				    << endl;
			}

			// Increment the nucleotide counter by 1.
			counter++;
		}

		// Close the output stream.
		out.close();
	}

	// Delete the data strand and error checker.
	delete rna;
	delete rnaChecker;

	// Return the completion string.
	return complete;
}

///////////////////////////////////////////////////////////////////////////////
// Write a Postscript file.
///////////////////////////////////////////////////////////////////////////////
void StructureBackend::writePostscriptFile( string outFile, int number ) {

	structureHandlers[number-1].writePostscript( outFile, false );
}

///////////////////////////////////////////////////////////////////////////////
// Write an SVG file.
///////////////////////////////////////////////////////////////////////////////
void StructureBackend::writeSVGFile( string outFile, int number ) {

	structureHandlers[number-1].writeSVG( outFile );
}
