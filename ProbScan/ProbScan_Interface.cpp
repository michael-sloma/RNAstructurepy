/*
 * A program that calculates the partition function for a strand of nucleic acids.
 * This strand of nucleic acids can be composed of either DNA or RNA.
 *
 * (c) 2009 Mathews Lab, University of Rochester Medical Center.
 * Written by Jessica S. Reuter
 */

#include "ProbScan_Interface.h"
#include <vector>
#include <string>
#include <algorithm>
///////////////////////////////////////////////////////////////////////////////
// Constructor.
//////////////////////////////////////////////////////////////////////////////
probscanInterface::probscanInterface() {

	// Initialize the "experimental" offset.
	experimentalOffset = 0.0;

	// Initialize the "experimental" scaling.
	experimentalScaling = 1.0;

	// Initialize the SHAPE intercept.
	intercept = -0.6;

	// Initialize the nucleic acid type.
	isRNA = true;

	// Initialize the maximum pairing distance between nucleotides.
	maxDistance = -1;

	// Initialize the SHAPE slope.
	slope = 1.8;

	// Initialize the calculation temperature.
	temperature = 310.15;

    //initialize the index variables
    i=j=k=l=1;
}

///////////////////////////////////////////////////////////////////////////////
// Parse the command line arguments.
///////////////////////////////////////////////////////////////////////////////
bool probscanInterface::parse( int argc, char** argv ) {

	// Create the command line parser and build in its required parameters.
	ParseCommandLine* parser = new ParseCommandLine( "ProbScan" );
	parser->addParameterDescription( "input file", "The name of a file containing a partition function or input sequence." );

	// Add the fromSequence option.
	vector<string> seqFileOptions;
	seqFileOptions.push_back( "-s" );
	seqFileOptions.push_back( "-S" );
	seqFileOptions.push_back( "--sequence" );
	parser->addOptionFlagsNoParameters( seqFileOptions, "Provide RNA from sequence file. Partition function will be calculated (may take a while)." );

	// Add the fromSequence option.
	vector<string> multibranchOptions;
	multibranchOptions.push_back( "-m" );
	multibranchOptions.push_back( "-M" );
	multibranchOptions.push_back( "--multibranch" );
	parser->addOptionFlagsWithParameters( multibranchOptions, "Provide a file with multibranch loops. These multibranch loops' probabilities will be checked." );

	// Add the maximum pairing distance option.
	vector<string> distanceOptions;
	distanceOptions.push_back( "-md" );
	distanceOptions.push_back( "-MD" );
	distanceOptions.push_back( "--maxdistance" );
	parser->addOptionFlagsWithParameters( distanceOptions, "Specify a maximum pairing distance between nucleotides. Default is no restriction on distance between pairs." );
    // Add the DNA option.
    vector<string> dnaOptions;
    dnaOptions.push_back( "-d" );
    dnaOptions.push_back( "-D" );
    dnaOptions.push_back( "--DNA" );
    parser->addOptionFlagsNoParameters( dnaOptions, "Specify that the sequence is DNA, and DNA parameters are to be used. Default is to use RNA parameters." );

	// Parse the command line into pieces.
	parser->parseLine( argc, argv );

	// Get required parameters from the parser.
	if( !parser->isError() ) {
		inputFile = parser->getParameter( 1 );
	}
	// Get the sequence file option.
	if( !parser->isError() ) { fromSequence = parser->contains( seqFileOptions ); }

	// Get the multibranch file option.
	if( !parser->isError() ) { 
        multibranch = parser->contains( multibranchOptions ); 
        loop_file = parser->getOptionString( multibranchOptions);
        //cout<<"loop file is "<<loop_file<<endl;
    }

    // Get the DNA option.
    if( !parser->isError() ) { isRNA = !parser->contains( dnaOptions ); }

	// Get the maximum distance option.
	if( !parser->isError() ) {
		parser->setOptionInteger( distanceOptions, maxDistance );
		bool badDistance =
		  ( maxDistance < 0 ) &&
		  ( maxDistance != -1 );
		if( badDistance ) { parser->setError( "maximum pairing distance" ); }
	}

	// Delete the parser and return whether the parser encountered an error.
	bool noError = ( parser->isError() == false );
	delete parser;
	return noError;
}

void probscanInterface::processLine(string input)
{
    vector<char> buff(input.size() + 1); 
    copy(input.begin(), input.end(), buff.begin());
    char* s = &buff[0];

    calcType = strtok(s," ,.\t");
    i = atoi(strtok(NULL," ,.\t"));
    j = atoi(strtok(NULL," ,.\t"));
    if (calcType=="internal"){
        k = atoi(strtok(NULL," ,.\t"));
        l = atoi(strtok(NULL," ,.\t"));
    }
}
//functions to split an input string, for multibranch calculation IO
std::vector<std::string> &spl(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}
//split string s in delimiter delim and return a vector of strings with the result
std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    spl(s, delim, elems);
    return elems;
}
//convert string of the form "#-#" and returns a std::pair<int,int>
//eg "10-20" -> std::pair(10,20)
std::pair<int,int> stringtopair(const std::string s){
    std::vector<std::string> str_indices = split(s,'-');
//    cout<<s<<endl;
//    cout<<str_indices[0]<<endl;
//    cout<<str_indices[1]<<endl;
//    cout<<atoi(str_indices[0].c_str())<<endl;
    std::pair<int,int> p = std::make_pair(atoi(str_indices[0].c_str()),
                              atoi(str_indices[1].c_str()));
    return p;
}
//takes a string of form "#-#\t#-#\t ... #-#" and returns a vector of std::pair<int,int>
multibranch_loop_t stringtombl(const std::string s){
    std::vector<std::string> pairs = split(s,'\t');
    multibranch_loop_t mb;
    for(std::vector<std::string>::iterator it=pairs.begin()+1;it!=pairs.end();++it){
        mb.branches.push_back(stringtopair(*it));
    }
    return mb;
}


///////////////////////////////////////////////////////////////////////////////
// Run calculations.
///////////////////////////////////////////////////////////////////////////////
void probscanInterface::run() {

	// Create a variable that handles errors.
	int error = 0;

	/*
	 * Use the constructor for RNA that specifies a filename.
	 * Specify type = 2 (sequence file).
	 * isRNA identifies whether the strand is RNA (true) or DNA (false).
	 *
	 * After construction of the strand data structure, create the error checker which monitors for errors.  
	 * Throughout, the error status of the calculation is checked with a variant of the isErrorStatus method, which returns 0 if no error occurred.
	 * The calculation proceeds as long as error = 0.
	 */
    try{
//	cout << "Initializing nucleic acids..." << flush;
	ProbScan* ps = new ProbScan(inputFile.c_str(),fromSequence?2:3,isRNA);
	ErrorChecker<RNA>* checker = new ErrorChecker<RNA>( ps );
	error = checker->isErrorStatus();
//	if( error == 0 ) { cout << "done." << endl; }
/*    string input;
    cout << ">> ";
    while(getline(cin,input)){//horrible hacked together interface for testing. needs ot be made user friendly!
        processLine(input);
        if (calcType=="hairpin")
            cout<< ps->probability_of_individual_hairpin(i,j)<<endl;
        if (calcType=="internal"&&(i<k<l<j))
            cout<< ps->probability_of_internal_loop(i,j,k,l)<<endl;
        cout << ">> ";
    }
*/
    if(!multibranch){
        cout << ps->GetStructure()->GetSequenceLabel()<<endl;
        show_hairpins(ps->probability_of_all_hairpins(3,ps->GetSequenceLength()-2,0.01));
        show_internal_loops(ps->probability_of_all_internal_loops(0.01));
    }

//    ps->probscan_test();
    if(multibranch){
        ifstream infile(loop_file.c_str());
        if (!infile.good()) throw "failed to open file\n";
        std::string line;
        while(getline(infile,line)){
            multibranch_loop_t mbl = stringtombl(line);
//            for(multibranch_loop_t::iterator it = mbl.begin();it!=mbl.end();++it)
//                cout<<it->first<<" "<<it->second<<"\n";
            mbl.probability = ps->probability_of_multibranch_loop(mbl);
            show_mbl(mbl);
        }
    }
//    ps->probscan_test();
	// Delete the error checker and data structure.
	delete checker;
	delete ps;

	// Print confirmation of run finishing.
//	if( error == 0 ) { cout << calcType << " complete." << endl; }
//	else { cerr << calcType << " complete with errors." << endl; }
    }
    catch(const char* oops){
	cout<<oops;
    }

}

///////////////////////////////////////////////////////////////////////////////
// Main method to run the program.
///////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] ) {

	probscanInterface* runner = new probscanInterface();
	bool parseable = runner->parse( argc, argv );
	if( parseable == true ) { runner->run(); }
	delete runner;
	return 0;
}
