/*
 * ShapeKnots, a program that predicts RNA secondary structures with pseudoknots.
 *
 * (c) 2013 
 * Mathews Lab, University of Rochester Medical Center
 * Weeks Lab, The University at North Carolina at Chapel Hill
 * Code contributors: Wayne Higgins, Stanislav Bellaousov, David H. Mathews
 */

#include "ShapeKnots_Interface.h"
#define OUTPUT_TO_SCREEN
//#undef OUTPUT_TO_SCREEN

///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
ShapeKnots_Interface::ShapeKnots_Interface() {

	//Initialize the differntial SHAPE slope.
	Dslope = 2.11;
	
	// Initialize the maximum number of possible pseudoknotted helices (used in function convertToHunderHelices())
	finallistSize=100; 

	// Initialize the maximum number of structures to be generated internally.
	InMaxStructures=100;

	// Initialize the maximum percent energy difference for the generated structures.
	InPercent=20;

	// Initialize the SHAPE intercept (kcal/mol).
	intercept=-0.6;

	// Initialize the folding window size or how different internal suboptimal structures can be.
	InWindowSize=0;

	// Initialize the maximum number of structures to be outputted.
	OutMaxStructures=20;

	// Initialize the maximum percent energy difference for the outputted suboptimal structures.
	OutPercent=10;

	// Initialize the folding window size or hwo different internal suboptimal structures can be.
	// This value will be changed based on the length of the sequence.
	OutWindowSize=0;

	// Initialize pseudoknot energy model parameters (kcal/mol).
	P1=0.35;

	// Initialize pseudoknot energy model parameters (kcal/mol).
	P2=0.65;

	// Initialize the SHAPE slope (kcal/mol).
	slope=1.8;


}

bool ShapeKnots_Interface::Parse(int argc, char** argv){
	
	// Create the command line parser and build in its required parameters.
	ParseCommandLine* parser = new ParseCommandLine( "ShapeKnots" );
	parser->addParameterDescription( "seq file", "The name of a sequence file containing input data. Note that lowercase nucleotides are forced single-stranded in structure prediction." );
	parser->addParameterDescription( "ct file", "The name of a CT file to which output will be written." );

    // Add the constraint file option.
	vector<string> constraintOptions;
	constraintOptions.push_back( "-c" );
	constraintOptions.push_back( "-C" );
	constraintOptions.push_back( "--constraint" );
	parser->addOptionFlagsWithParameters( constraintOptions, "Specify a constraints file to be applied. Default is to have no constraints applied." );

    // Add the DMS option.
	vector<string> dmsOptions;
	dmsOptions.push_back( "-dms" );
	dmsOptions.push_back( "-DMS" );
	dmsOptions.push_back( "--DMS" );
	parser->addOptionFlagsWithParameters( dmsOptions, "Specify a DMS constraints file to be applied. These constraints are pseudoenergy constraints. Default is to have no constraints applied." );

	// Add the differential SHAPE restraint file option.
	vector<string> DshapeOptions;
	DshapeOptions.push_back( "-dsh" );
	DshapeOptions.push_back( "-DSH" );
	DshapeOptions.push_back( "--DSHAPE" );
	parser->addOptionFlagsWithParameters(DshapeOptions, "Specify a differential SHAPE restraints file to be applied. These constraints are pseudoenergy restraints. Default is to have no restraints applied." );

	// Add the differential SHAPEslope option.
	vector<string> DshapeSlopeOptions;
	DshapeSlopeOptions.push_back( "-dsm" );
	DshapeSlopeOptions.push_back( "-DSM" );
	DshapeSlopeOptions.push_back( "--DSHAPEslope" );
	parser->addOptionFlagsWithParameters( DshapeSlopeOptions, "Specify a slope used with differential SHAPE restraints. Default is 2.11 kcal/mol." );

	// Add the double offset file option.
	vector<string> doubleOffsetOptions;
	doubleOffsetOptions.push_back( "-dso" );
	doubleOffsetOptions.push_back( "-DSO" );
	doubleOffsetOptions.push_back( "--doubleOffset" );
	parser->addOptionFlagsWithParameters( doubleOffsetOptions, "Specify a double-stranded offset file, which adds energy bonuses to particular double-stranded nucleotides. Default is to have no double-stranded offset specified." );

	// Add the maximum number of structures option.                                                                                                                                                        
	vector<string> INmaxStructuresOptions;
	INmaxStructuresOptions.push_back( "-im" );
	INmaxStructuresOptions.push_back( "-IM" );
	INmaxStructuresOptions.push_back( "--IMaximum" );
	parser->addOptionFlagsWithParameters( INmaxStructuresOptions, "Specify a maximum number of internally generated structures for each call of the dynamic programming algorithm. Note that suboptimal structures are generated until either the maximum number of structures is reached or the maximum percent difference is reached (below).  This is not the maximum number of structures provided to the user, which is controlled by –m, -M, --maximum. Default is 100 structures." );
	
	// Add the percent energy difference option.
	vector<string> INpercentOptions;
	INpercentOptions.push_back( "-ip" );
	INpercentOptions.push_back( "-IP" );
	INpercentOptions.push_back( "--IPercent" );
	parser->addOptionFlagsWithParameters( INpercentOptions, "Specify a maximum percent difference in folding free energy change for internally generated suboptimal structures for each call of the dynamic programming algorithm. For example, 20 would indicate 20%. This is not the maximum percent difference in energy for structures provided to the user, which is controlled by –p, -P, --percent. Default is 20%." );
	
	// Add the window size option.                                                                                                                                                                         
	vector<string> INwindowOptions;
	INwindowOptions.push_back( "-iw" );
	INwindowOptions.push_back( "-IW" );
	INwindowOptions.push_back( "--IWindow" );
	parser->addOptionFlagsWithParameters( INwindowOptions, "Specify a window size for the internally generated suboptimal structures for each call of the dynamic programming algorithm.  This is not the window for structures provided to the user, which is controlled by –w, -W, --window. Default is determined by the length of the sequence." );

	// Add the maximum number of structures option.                                                                                                                                                        
	vector<string> OUTmaxStructuresOptions;
	OUTmaxStructuresOptions.push_back( "-m" );
	OUTmaxStructuresOptions.push_back( "-M" );
	OUTmaxStructuresOptions.push_back( "--maximum" );
	parser->addOptionFlagsWithParameters( OUTmaxStructuresOptions, "Specify a maximum number of structures to be outputted. Note that suboptimal structures are generated until either the maximum number of structures is reached or the maximum percent difference is reached (below). Default is 20 structures." );
	
	// Add the percent energy difference option.
	vector<string> OUTpercentOptions;
	OUTpercentOptions.push_back( "-p" );
	OUTpercentOptions.push_back( "-P" );
	OUTpercentOptions.push_back( "--percent" );
	parser->addOptionFlagsWithParameters( OUTpercentOptions, "Specify a maximum percent difference in folding free energy change for generating suboptimal structures in the output. For example, 10 would indicate 10%. Default is 10%." );

	// Add the Penalty1 option.
	vector<string> Penalty1Options;
	Penalty1Options.push_back( "-p1" );
	Penalty1Options.push_back( "-P1" );
	Penalty1Options.push_back( "--Penalty1" );
	parser->addOptionFlagsWithParameters( Penalty1Options, "Specify a pseudoknot penalty P1. Default is 0.35 kcal/mol." );

	// Add the Penalty2 option.
	vector<string> Penalty2Options;
	Penalty2Options.push_back( "-p2" );
	Penalty2Options.push_back( "-P2" );
	Penalty2Options.push_back( "--Penalty2" );
	parser->addOptionFlagsWithParameters( Penalty2Options, "Specify a pseudoknot penalty P2. Default is 0.65 kcal/mol." );

	// Add the finallistSize option.
	vector<string> finallistSizeOptions;
	finallistSizeOptions.push_back( "-ph" );
	finallistSizeOptions.push_back( "-PH" );
	finallistSizeOptions.push_back( "--PseudoknottedHelices" );
	parser->addOptionFlagsWithParameters( finallistSizeOptions, "Specify maximum number of helices to be processed. Default is 100 helices." );

	// Add the SHAPE restraint file option.
	vector<string> shapeOptions;
	shapeOptions.push_back( "-sh" );
	shapeOptions.push_back( "-SH" );
	shapeOptions.push_back( "--SHAPE" );
	parser->addOptionFlagsWithParameters(shapeOptions, "Specify a SHAPE restraints file to be applied. These restraints specifically use SHAPE pseudoenergy restraints. Default is no SHAPE restraint file specified." );

	// Add the SHAPEintercept option.
	vector<string> shapeInterceptOptions;
	shapeInterceptOptions.push_back( "-si" );
	shapeInterceptOptions.push_back( "-SI" );
	shapeInterceptOptions.push_back( "--SHAPEintercept" );
	parser->addOptionFlagsWithParameters( shapeInterceptOptions, "Specify an intercept used with SHAPE restraints. Default is -0.6 kcal/mol." );

	// Add the SHAPEslope option.
	vector<string> shapeSlopeOptions;
	shapeSlopeOptions.push_back( "-sm" );
	shapeSlopeOptions.push_back( "-SM" );
	shapeSlopeOptions.push_back( "--SHAPEslope" );
	parser->addOptionFlagsWithParameters( shapeSlopeOptions, "Specify an slope used with SHAPE restraints. Default is 1.8 kcal/mol." );

	// Add the single offset file option.
	vector<string> singleOffsetOptions;
	singleOffsetOptions.push_back( "-sso" );
	singleOffsetOptions.push_back( "-SSO" );
	singleOffsetOptions.push_back( "--singleOffset" );
	parser->addOptionFlagsWithParameters( singleOffsetOptions, "Specify a single-stranded offset file, which adds energy bonuses to particular single-stranded nucleotides. Default is to have no single-stranded offset specified." );
	
	// Add the window size option.                                                                                                                                                                         
	vector<string> OUTwindowOptions;
	OUTwindowOptions.push_back( "-w" );
	OUTwindowOptions.push_back( "-W" );
	OUTwindowOptions.push_back( "--window" );
	parser->addOptionFlagsWithParameters( OUTwindowOptions, "Specify a window size for outputted suboptimal structures. Default is determined by the length of the sequence." );

	// Parse the command line into pieces.
	parser->parseLine( argc, argv );

	// Get required parameters from the parser.                                                                                                                                                            
	if( !parser->isError() ) {
		seqFile = parser->getParameter( 1 );
		ctFile = parser->getParameter( 2 );
	}

	// Get the constraint file option.
	if( !parser->isError() ) { constraintFile = parser->getOptionString( constraintOptions, true ); }

	// Get the differential SHAPE slope option.
	if( !parser->isError() ) {
		parser->setOptionDouble( DshapeSlopeOptions, Dslope );
	}

	// Get the double strand offset file option.
	if( !parser->isError() ) { doubleOffsetFile = parser->getOptionString( doubleOffsetOptions, true ); }

	// Get the maximum number of structures to be internally generated.
	if( !parser->isError() ) {
		parser->setOptionInteger( INmaxStructuresOptions, InMaxStructures );
		if( OutMaxStructures <= 0 ) { parser->setError( "maximum number of structures" ); }
	}

	// Get the percent energy difference option for internally generated structures.
	if( !parser->isError() ) {
		parser->setOptionInteger( INpercentOptions, InPercent );
		if( OutPercent < 0 ) { parser->setError( "percent energy difference" ); }
	}
	
	// Get the window size option for internally generated structures.
	if( !parser->isError() ) {
		parser->setOptionInteger( INwindowOptions, InWindowSize );
		if( OutWindowSize < 0 ) { parser->setError( "window size" ); }
	}

	// Get the maximum number of structures to be outputted option.
	if( !parser->isError() ) {
		parser->setOptionInteger( OUTmaxStructuresOptions, OutMaxStructures );
		if( OutMaxStructures <= 0 ) { parser->setError( "maximum number of structures" ); }
	}

	// Set modifier type
	if( !parser->isError() ) {
		if(parser->contains(dmsOptions))
		    modifier = "DMS"; 
		if(parser->contains(shapeOptions))
		    modifier = "SHAPE"; 
	}

	// Get the SHAPE, and DMS data and options.
	if( !parser->isError() ) {
		SHAPEFile = parser->getOptionString( shapeOptions );
		DSHAPEFile = parser->getOptionString( DshapeOptions );
        DMSFile = parser->getOptionString( dmsOptions ); 
        if( !parser->isError() ) { parser->setOptionDouble( shapeInterceptOptions, intercept ); }
        if( !parser->isError() ) { parser->setOptionDouble( shapeSlopeOptions, slope ); }
		if( !parser->isError() ) { parser->setOptionDouble( DshapeSlopeOptions, Dslope ); }
	}
    
	// Get the percent energy difference option for outputted suboptimal structures.
    if( !parser->isError() ) {
        parser->setOptionInteger( OUTpercentOptions, OutPercent );
        if( OutPercent < 0 ) { parser->setError( "percent energy difference" ); }
	}

	// Get the pseudoknot penalty 1.
	if( !parser->isError() ) {
		parser->setOptionDouble( Penalty1Options, P1 );
	}

	// Get the pseudoknot penalty 2.
	if( !parser->isError() ) {
		parser->setOptionDouble( Penalty2Options, P2 );
	}

	// Get the finallistSize option.
	if( !parser->isError() ) {
		parser->setOptionInteger( finallistSizeOptions, finallistSize );
		if( finallistSize < 1 ) { parser->setError( "maximum number of possible pseudoknotted helices" ); }
	}

	// Get the single strand offset file option.
	if( !parser->isError() ) { singleOffsetFile = parser->getOptionString( singleOffsetOptions, true ); }
	
	// Get the window size option for outputted structures.
	if( !parser->isError() ) {
		parser->setOptionInteger( OUTwindowOptions, OutWindowSize );
		if( OutWindowSize < 0 ) { parser->setError( "window size" ); }
	}
	//Record if the windowOption was set
	ifWindowOptions=parser->contains(OUTwindowOptions);
	
	// Delete the parser and return whether the parser encountered an error.
	bool noError = ( parser->isError() == false );
	delete parser;
	return noError;
}


///////////////////////////////////////////////////////////////////////////////
// Run calculations.
///////////////////////////////////////////////////////////////////////////////
void ShapeKnots_Interface::run() {
	
	char loop[maxfil],stackf[maxfil],tstackh[maxfil],tstacki[maxfil],
		tloop[maxfil],miscloop[maxfil],danglef[maxfil],int22[maxfil],
		int21[maxfil],coax[maxfil],tstackcoax[maxfil],coaxstack[maxfil],tstack[maxfil],tstackm[maxfil],triloop[maxfil],int11[maxfil],hexaloop[maxfil],
		tstacki23[maxfil], tstacki1n[maxfil],datapath[maxfil],*pointer;
	
#ifdef OUTPUT_TO_SCREEN
	cout << "############################################\n";
	cout << "         ShapeKnots Folding Started\n";
	cout << "############################################\n";
	cout << "Reading Files..." << flush ;
#endif
	
	datatable *data=new datatable;
	pointer = getenv("DATAPATH");
	if (pointer!=NULL){
		strcpy(datapath,pointer);//open the data files -- must reside in the DATAPATH
		strcat(datapath,"/");//The path to DATAPATH has to end with '/'
	}		
	else strcpy(datapath,"");//open the data files -- must reside in the same directory as the executable
	//open the thermodynamic data tables
	getdat (loop, stackf, tstackh, tstacki,tloop, miscloop, danglef, int22,
	        int21,coax, tstackcoax,coaxstack, tstack, tstackm, triloop,
	        int11, hexaloop, tstacki23, tstacki1n, datapath, true);//the true indicates RNA parameters
	
	if (opendat (loop,stackf,tstackh,tstacki,tloop,miscloop,danglef,int22,int21,
	             coax,tstackcoax,coaxstack,tstack,tstackm,triloop,int11,hexaloop,tstacki23,tstacki1n,data)==0) {
		
		cerr << "A data file was lost";
	}
	
	//Check if the files exist:
	
	//Check if the sequence file exists
	ifstream inseq(seqFile.c_str());
	if(!inseq) {
		cerr << " SEQUENCE INPUT FILE NOT FOUND\n";
		delete data;
		inseq.close();
		return;
	}
	
    RNA *rnaCT=new RNA;//Initialize RNA class
    structure *ct=rnaCT->GetStructure();//Create a 'structure' pointer pointing to rnaCT.
    
	//Store the path to output CT file
	ct->openseq(seqFile.c_str());//open the sequence

	//If the user has specified a SHAPE restraints file, read the file.
	if(!SHAPEFile.empty()){
		//Read in the shape data file and convert the data to a linear penalty 
		ct->SHAPEslope = slope*conversionfactor;//read the slope
		ct->SHAPEintercept = intercept*conversionfactor;//read the intercept
		ct->ReadSHAPE(SHAPEFile.c_str());
	}

    //If the user has specified a DSHAPE restraints file, read the file.
	if(!DSHAPEFile.empty()){
		//Read in the shape data file and convert the data to a linear penalty 
		ct->SHAPEslope = Dslope*conversionfactor;//read the slope
		ct->SHAPEintercept = 0;//read the intercept
		ct->ReadSHAPE(DSHAPEFile.c_str(),"diffSHAPE");
	}

    //If the user has specified a DMS restraints file, read the file.
	if(!DMSFile.empty()){
		//Read in the dms data file and convert the data to a linear penalty 
        ct->ReadSHAPE(DMSFile.c_str(),"DMS");
	}

 	//Read the Double-strand offset data into 'ct'
	if(!doubleOffsetFile.empty()) ct->ReadOffset(NULL,doubleOffsetFile.c_str());//read the offset data

	//Read the Single-strand offset data into 'ct'
	if(!singleOffsetFile.empty()) ct->ReadOffset(NULL,singleOffsetFile.c_str());//read the offset data
	
	//The rest of the code is located in the 'pseudoknot' function.
	pseudoknot(rnaCT, data, InMaxStructures, InPercent, InWindowSize, ctFile, P1, P2, slope, intercept, DMSFile, SHAPEFile, Dslope, DSHAPEFile, doubleOffsetFile, OutPercent, OutWindowSize, OutMaxStructures, ifWindowOptions, finallistSize);
	
	delete rnaCT;
	delete data;
}

///////////////////////////////////////////////////////////////////////////////
// Main method to run the program.
///////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] ) {
	
	ShapeKnots_Interface* runner = new ShapeKnots_Interface();
	bool parseable = runner->Parse( argc, argv );
	if( parseable == true ) { runner->run(); }
	delete runner;
	return 0;
}

