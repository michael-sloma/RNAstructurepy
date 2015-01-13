/*
 * NAPSS, a program that predicts RNA secondary structures with pseudoknots with the aid of NMR constraints data.
 *
 * (c) 2012 Mathews Lab, University of Rochester Medical Center.
 * Written by James Hart, David Mathews, Jon Chen, Stanislav Bellaousov
 *
 * Modified from the command-line interface to Dynalign, 
 * written by David Mathews; Copyright 2002, 2003, 2004, 2005, 2006
 *
 * Modified to parse command line interface by Jessica Reuter and Stanislav Bellaousov 2012
 *
 * Modified to read Triplet Constraints by Stanislav Bellaousov in August 2012
 *
 */

#include "napss.h"

// Flags for text output
#undef VERBOSE_MODE
//#define VERBOSE_MODE

// Flags for output associated with Triplet Constraints
#undef TRIPLET_VERBOSE_MODE
//#define TRIPLET_VERBOSE_MODE

// Flags for debug features
#undef DEBUG_MODE
//#define DEBUG_MODE
#undef ALREADY_USED_CHECK
//#define ALREADY_USED_CHECK


// The main entry point for NAPSS using a text interface:
int main(int argc, char* argv[]) {

	char loop[maxfil],stackf[maxfil],tstackh[maxfil],tstacki[maxfil],
		tloop[maxfil],miscloop[maxfil],danglef[maxfil],int22[maxfil],
		int21[maxfil],coax[maxfil],tstackcoax[maxfil],
		coaxstack[maxfil],tstack[maxfil],tstackm[maxfil],triloop[maxfil],int11[maxfil],hexaloop[maxfil],
		tstacki23[maxfil], tstacki1n[maxfil],datapath[maxfil],*pointer;
	structure *ct;
	datatable data;
	int i,j,iter,iter2;
	string a;

	//Parse the commans line inputs
	bool parsable = Parse(argc, argv);
	if(parsable==false)	return 0;

	cout << "####################################################\n";
	cout << "                  NAPSS Started\n";
	cout << "####################################################\n\n";
	cout << "Reading files..." << flush ;

	// Get the location of the data files
	pointer = getenv("DATAPATH");
	if (pointer!=NULL) {
		strcpy(datapath,pointer);
		strcat(datapath,"/");
	}
	else strcpy(datapath,"");//if the DATAPATH is not specified, look in current directory

	// Open the thermodynamic data tables
	getdat (loop, stackf, tstackh, tstacki,tloop, miscloop, danglef, int22,
	        int21,coax, tstackcoax,coaxstack, tstack, tstackm, triloop,
	        int11, hexaloop, tstacki23, tstacki1n, datapath, true);//the true indicates RNA parameters
	if (opendat (loop, stackf, tstackh, tstacki,tloop, miscloop, danglef, int22,
	             int21,coax, tstackcoax,coaxstack, tstack, tstackm, triloop,
	             int11,hexaloop,tstacki23, tstacki1n, &data)==0) cerr << "A data file was lost\n";
	
	//Initialize RNA class with RNA sequence.
	RNA *rnaCT=new RNA(inseq.c_str(),2);
	//Check for errors
	if (rnaCT->GetErrorCode()!=0) {
		//If there was an error output the error and exit
		cerr << rnaCT->GetErrorMessage(rnaCT->GetErrorCode())<<"\n";
		delete rnaCT;
		//delete ct;
		return 1;
	}
	
	//Read the SHAPE data from the disk
	if(!inSHAPEfile.empty()) rnaCT->ReadSHAPE(inSHAPEfile.c_str(),slope,intercept);

	//Read the constraint data from the disk
	if(!constraintFile.empty()) rnaCT->ReadConstraints(constraintFile.c_str());

	//Set pointer 'ct' to point to 'structure' in 'rna'
	ct=rnaCT->GetStructure();

#if defined(VERBOSE_MODE)
	// Display basic info
	cout << "\n\n";
	for (i = 1; i <= ct->numofbases; i++) cout << ct->nucs[i];
	cout << "\n";
	for (i = 1; i <= ct->numofbases; i++) cout << ct->numseq[i];
	cout << "\nLength of sequence: " << ct->numofbases << "\n";
#endif

	// Create a basepair type lookup matrix (AU = 5, GC = 6, GU = 7)
	short bpLookup[5][5] = {{0,0,0,0,0},
	                        {0,0,0,0,5},
	                        {0,0,0,6,0},
	                        {0,0,6,0,7},
	                        {0,5,0,7,0}};

	// Create 2D matrices to store converted dotplot, DeltaG, and vmb/vext values
	short** convertedDotplot;
	short** dgArray;
	short** mbDotplot;
	convertedDotplot = new short* [(ct->GetSequenceLength())+1];
	dgArray = new short* [(ct->GetSequenceLength())+1];
	mbDotplot = new short* [(ct->GetSequenceLength())+1];	
	for(i = 0; i <= ct->GetSequenceLength(); i++) {
		*(convertedDotplot + i) = new short[(ct->GetSequenceLength())+1];
		*(dgArray + i) = new short[(ct->GetSequenceLength())+1];
		*(mbDotplot + i) = new short[(ct->GetSequenceLength())+1];
	}
	
	// Clear the values stored in the matrices
	for (i = 0; i <= ct->GetSequenceLength(); i++) {
		for (j = 0; j <= ct->GetSequenceLength(); j++) {
			convertedDotplot[i][j] = 0;
			dgArray[i][j] = 0;
			mbDotplot[i][j] = 0;
		}
	}

	//allocate space for v, vmb, and vext arrays:
	arrayclass v(ct->GetSequenceLength());
	arrayclass vmb(ct->GetSequenceLength());
	arrayclass vext(ct->GetSequenceLength());
	double pairs;
	double mbpairs;
	double extpairs;
	int dgMin = 0;

	// Load NMR constraints
	short maxConLength = 0, numOfCon = 0, totalConLength = 0;
	
	ifstream inNMRConFile;
	inNMRConFile.open(inNMRconstraints.c_str());
	if (!inNMRConFile) {
		cerr << "Unable to open constraints file";
		cin >> i;
		exit(1);
	}

	string s;

	// Loop through each line of the file to get count of constraints and maximum length of a constraint
	while (getline(inNMRConFile,s)) {
		if(!s.empty()&&s.length() < 2){cerr<<"Constraints cannot be less than two basepairs in length!\n";cin>>i;return -1;}
		if(!s.empty()){//make sure not to read empty lines
			numOfCon++;
			//Counting the constraints omitting the constraints in parentheses
			const char *TempRead = s.c_str();//Read the string into character array
			short TempConLength=0;//Re-Set the constraint length counter
			while(*TempRead != '\0'){//Check if it is the end of array. If not, execute while loop
				if(*TempRead++ != '(')++TempConLength;//if not the parenthesis, count the constraints
				else{//if find parenthesis
					while(*TempRead++ != ')');//skip till the closing parenthesis
				}
			}
#if defined(TRIPLET_VERBOSE_MODE)
			cerr << "-------------\nTempConLength = " << TempConLength << endl;	//debug	
#endif
			totalConLength+=TempConLength;//add to the total constraint length counter
			if(TempConLength > maxConLength) maxConLength = TempConLength;//find the longest constraint
			
		}//END:Making sure not to read empty lines
	}//END:while loop that reads file line by line
	inNMRConFile.close();
	inNMRConFile.clear();
	
#if defined(VERBOSE_MODE)
	cout << numOfCon << " NMR constraints\n";
	cout << "Total length of NMR constraints: " << totalConLength << "\n";
	cout << "Length of largest NMR constraint: " << maxConLength << "\n\n";
#endif


	// Allocate 2D array for storing constraints (first column contains length of each constraint)
	// Also allocate two more copies of this array which will be used later for storing temporary
	// dotplot coordinates
	char**** tripletArray;//Initialize a 4D array to hold the triplet constraints. 1D = non-triplet constraint number;
	                      //2D = constraint position; 3D = # of triplet constraint (could be 1 or 2);
	                      //4D = the triplet constraint, ex: +RGY
	short** conArray;
	short** xCoords;
	short** yCoords;
	tripletArray = new char*** [numOfCon+1];
	conArray = new short* [numOfCon+1];
	xCoords = new short* [numOfCon+1];
	yCoords = new short* [numOfCon+1];
	for(i = 0; i <= numOfCon; i++){
		tripletArray[i] = new char**[(maxConLength+1)];
		for(int j=0;j<=maxConLength;j++){
			tripletArray[i][j] = new char*[2];
			for(int k=0;k<2;k++){
				tripletArray[i][j][k] = new char[4];
				for(int m=0;m<4;m++) tripletArray[i][j][k][m] = 0;
			}
		}
		*(conArray + i) = new short[(maxConLength+1)];
		*(xCoords + i) = new short[(maxConLength+1)];
		*(yCoords + i) = new short[(maxConLength+1)];
	}
	// Re-open NMR constraints file and fill the 2D arrays (only fill Coords arrays with 0's)
	i = 0;
	inNMRConFile.open(inNMRconstraints.c_str());
	while (getline(inNMRConFile,s)) {
		if(!s.empty()){//make sure that the line is not empty
			while(s.find('(') != string::npos){//check if there are any open parenthesis in the string
				for(int tr=0;tr<4;tr++){//for the four positions of the triplet constraint
					tripletArray[i][s.find('(')][0][tr]=s.at(s.find('(')+tr+1);//store the triplet constraint in the 'tripletArray' array
				}
				if(s.at(s.find('(')+5)=='/'){//check if there is a second constraint separated by '/'
					for(int tr=0;tr<4;tr++){//for the four positions of the triplet constraint
						tripletArray[i][s.find('(')][1][tr]=s.at(s.find('(')+tr+6);//store the triplet constraint in the 'tripletArray' array
					}
					s.erase(s.find('('),s.find(')')-s.find('(')+1);
				}//END: if there is a second constraint
				else(s.erase(s.find('('),s.find(')')-s.find('(')+1));//erase the triplet constraint from the string
			}//END: while loop that is looking for triplet constraints
		
			conArray[i][0] = s.size();
			xCoords[i][0] = 0;
			yCoords[i][0] = 0;
			
			for (j = 1; j <= s.size(); j++) {
				conArray[i][j] = atoi(s.substr(j-1,1).c_str());
				xCoords[i][j] = 0;
				yCoords[i][j] = 0;
			}
			i++;		
		}
	}
	inNMRConFile.close();
	inNMRConFile.clear();

#if defined(VERBOSE_MODE)
	// Display the constraints array
	cout << "NMR Constraints:\n";
	for (int im = 0; im < numOfCon; im++) {
		for (int jm = 1; jm <= conArray[im][0]; jm++) {
			cout << conArray[im][jm];
			for (int km = 0; km <2; km++){
				for(int m=0;m<4;m++) cout << tripletArray[im][jm][km][m];
			}
		}
		cout << "\n";
	}
#endif

	// Create 1D bool array for storing the nucleotides that are being considered for a match
	bool* alreadyUsed = new bool[ct->GetSequenceLength()+1];
	for (i = 0; i <= ct->GetSequenceLength(); i++) {
		alreadyUsed[i] = false;
	}

	// Create storage container for matches
	vector<conMatch> matchVector;

	// Create 1D bool array for storing information about which constraints are symmetric
	bool* isSymmetric = new bool[numOfCon];
	for (i = 0; i < numOfCon; i++) {
		isSymmetric[i] = false;
	}

	// Test to see which constraints are symmetric
	bool tempEquality;
	for (i = 0; i < numOfCon; i++) {
		j = 0;
		tempEquality = true;
		while (j < conArray[i][0] && tempEquality) {
			if (conArray[i][conArray[i][0]-j] != conArray[i][j+1]) {
				tempEquality = false;
				continue;
			}
			j++;
			if (j == conArray[i][0]) {
				isSymmetric[i] = true;

#if defined(VERBOSE_MODE)
				cout << "Symmetric constraint: " << i << "\n";
#endif

			}
		}
	}
	short currConNum = 0;
	short currConPos = 1;

	cout << "\t\t\t\tDONE\nRunning dynamic programming algorithm..." << flush;

	// Fill arrays 'v', 'vmb', and 'vext' with thermodynamic data
	dynamic(ct,&data,maxtracebacks,99,windowsize,&v,&vmb,&vext);

	cout << "\tDONE\nGenerating dotplot..." << flush;

	// Calculate the dot plot information that tracks MB loops and exterior loops
	for (i=1;i<=ct->GetSequenceLength();i++) {
		for (j=(i+1);j<=ct->GetSequenceLength();j++) {
			if ((v.f(i,j)+v.f(j,i+ct->GetSequenceLength()))<0) {
				pairs= ((double)(v.f(i,j)+v.f(j,i+ct->GetSequenceLength()))/10.0);
				mbpairs=min( ((double)(vmb.f(i,j)+v.f(j,i+ct->GetSequenceLength()))/10.0), ((double)(v.f(i,j)+vmb.f(j,i+ct->GetSequenceLength()))/10.0));
				extpairs= ((double)(v.f(i,j)+vext.f(j,i+ct->GetSequenceLength()))/10.0);
				//cerr << i << "\t" << j << "\t" << pairs<<"\t"<< mbpairs<< "\t"<<extpairs<<"\n"; 
				//             if (instream >> x         >> y            >> z         >> vmb           >> vext) {
				convertedDotplot[i][j] = bpLookup[ct->numseq[i]][ct->numseq[j]];
				dgArray[i][j] = short(pairs*10);
				mbDotplot[i][j] = min(short(mbpairs*10),short(extpairs*10));
				// Copy these values to the other half of the matrices
				convertedDotplot[j][i] = bpLookup[ct->numseq[i]][ct->numseq[j]];
				dgArray[j][i] = short(pairs*10);
				mbDotplot[j][i] = min(short(mbpairs*10),short(extpairs*10));
				if(pairs*10 < dgMin) dgMin = short(pairs*10); // Finds and stores minimum free energy value from dotplot
				//cerr << i << "\t" << j << "\t" << dgArray[i][j] << "\t" << mbDotplot[i][j] << "\t" << convertedDotplot[i][j] << endl;
			}
		}
	}
	// Store cutoff DG value
	int dgCutoff = dgMin*(100-percent)/100;
	// Loop back through convertedDotplot and remove those values that have a DG value greater than dgCutoff
	for (i = 0; i <= ct->GetSequenceLength(); i++) {
		for (j = 0; j <= ct->GetSequenceLength(); j++) {
			if (dgArray[i][j] > dgCutoff) convertedDotplot[i][j] = 0;
		}
	}

	cout << "\t\t\t\tDONE\nLooking for Matches: " << flush;
	
#if defined(VERBOSE_MODE)
	// Display the converted dotplot
	cout << "\nConverted Dotplot:\n";
	for(i = 0; i <= ct->numofbases; i++) {
		for(j = 0; j <= ct->numofbases; j++) cout << convertedDotplot[i][j];
		cout << "\n";
	}
	cout << "\n";
#endif

	// Ready to begin recursive constraint matching
	firstDotMatch(conArray,alreadyUsed,isSymmetric,totalConLength,numOfCon,xCoords,yCoords,&matchVector,
	              &currConNum,&currConPos,convertedDotplot,mbDotplot,dgCutoff,ct->GetSequenceLength(),ct,tripletArray);

	// When finished with constraint matching
	cout << "\rLooking for matches: " << matchVector.size() << " matches found\tDONE\n" << flush;

#if defined(VERBOSE_MODE)
	// Display the matches
	conMatch tempMatchVector;
	for (i = 0; i < matchVector.size(); i++) {
		tempMatchVector = matchVector[i];
		for (j = 0; j < totalConLength; j++) {
			//			cout << tempMatchVector.xCoords[j*2] << "," << tempMatchVector.xCoords[j*2+1] << " ";
		}
		cout << "\n";
	}
#endif

	// Determine if refolding should occur
	if (matchVector.size() == 0) {
		cout << "No possible constraint matches - NAPSS will terminate.\n";
		return(0);
	}
	
	// Initialize ct's that will be used for storing all refolded structures and for calculating pseudoknot energies
	structure *ct2;
	structure *ctbroken;

	//Initialize RNA class with RNA sequence.
	RNA *rnaCT2=new RNA(inseq.c_str(),2);
	RNA *rnaCTBROKEN=new RNA(inseq.c_str(),2);
	//Check for errors
	if (rnaCT2->GetErrorCode()!=0||rnaCTBROKEN->GetErrorCode()!=0) {
		//If there was an error output the error and exit
		cerr << rnaCT2->GetErrorMessage(rnaCT2->GetErrorCode())<<"\n";
		cerr << rnaCTBROKEN->GetErrorMessage(rnaCTBROKEN->GetErrorCode())<<"\n";
		delete rnaCT,ct,rnaCT2,ct2,rnaCTBROKEN,ctbroken;
		return 1;
	}

	//Read the SHAPE data from the disk
	if(!inSHAPEfile.empty()) rnaCT2->ReadSHAPE(inSHAPEfile.c_str(),slope,intercept);
	if(!inSHAPEfile.empty()) rnaCTBROKEN->ReadSHAPE(inSHAPEfile.c_str(),slope,intercept);
	//Read the constraint data from the disk
	if(!constraintFile.empty()) rnaCT2->ReadConstraints(constraintFile.c_str());
	if(!constraintFile.empty()) rnaCTBROKEN->ReadConstraints(constraintFile.c_str());

	//Set pointer 'ct2' to point to 'structure' in 'rnaCT2'
	ct2=rnaCT2->GetStructure();
	//Set pointer 'ctbroken' to point to 'structure' in 'rnaCTBROKEN'
	ctbroken=rnaCTBROKEN->GetStructure();

	//ct2->numofstructures = 0;
	//ctbroken->numofstructures = 0;

	int start = 0;
	int count = 0;
	// Also create a 1D array for temporarily storing possible extensions and a control switch for that loop
	int* helixExtend = new int[ct->GetSequenceLength()+1];
	bool extensionAdded = true;

	// Create 2D boolean array (triangular matrix with sides of length ct->numofbases+1) to store dotplot data
	ct->allocatetem();

	// Loop over all matches
	for (iter = 0; iter < matchVector.size(); iter++) {

		cout << "\rFolding structure " << iter+1 << " of " << matchVector.size() << flush;

		if(!pseudoknotFree){// If pseudoknot-free mode is NOT specified (default behavior)

			// Initialize all positions to false:
			for (i = 0; i < ct->GetSequenceLength(); i++) {
				for (j = i+1; j <= ct->GetSequenceLength(); j++) {
					ct->tem[j][i] = false;
				}
			}
			
#if defined(VERBOSE_MODE)
			cout << "Non-conflicting match " << iter+1 << ":\n";
#endif

			// Copy original dotplot into tem** - change all non-zero convertedDotplot values to true
			for (i = 0; i < ct->GetSequenceLength(); i++) {
				for (j = i+1; j <= ct->GetSequenceLength(); j++) {
					if (convertedDotplot[j][i] != 0) ct->tem[j][i] = true;
				}
			}

			// Create new array to temporarily hold matched basepairs (initialize all to 0)
			int* tempbasepr = new int[ct->GetSequenceLength()+1];
			for (i=0; i<=ct->GetSequenceLength(); i++) tempbasepr[i]=0;
		
			// Trim the dotplot according to the base pairs of the current match
			int trimming_i, trimming_j;

			// Loop over all base pairs in the match combination
			for (iter2 = 0; iter2 < totalConLength; iter2++) {
				trimming_j = matchVector[iter][(2*iter2)];
				trimming_i = matchVector[iter][(2*iter2)+1];

				// Store coordinates in temporary basepair array
				tempbasepr[trimming_i] = trimming_j;
				tempbasepr[trimming_j] = trimming_i;
			
				// Change all tem values to false where (i or j) = (trimming_i or trimming_j)
				for (j = 0; j <= trimming_i-1; j++) {
					ct->tem[trimming_i][j] = false;
				}
				for (i = trimming_i + 1; i <= ct->GetSequenceLength(); i++) {
					ct->tem[i][trimming_i] = false;
				}
				for (j = 0; j <= trimming_j-1; j++) {
					ct->tem[trimming_j][j] = false;
				}
				for (i = trimming_j + 1; i <=ct->GetSequenceLength(); i++) {
					ct->tem[i][trimming_j] = false;
				}
			}

			// Search for regions of dotplot that can be excluded
			pseudodptrim(ct, tempbasepr, &count);
			delete[] tempbasepr;

			// Do the structure prediction considering only the allowable basepairs from the trimmed dotplot
			dynamic(ct,&data,maxtracebacks,99,windowsize);

#if defined(VERBOSE_MODE)
			cout << "Refolding yields " << ct->numofstructures << " structures.\n\n";
#endif
			// Reinsert the basepairs that match the NMR constraints
			for (i = 1; i <= ct->GetNumberofStructures(); i++) {
				for (iter2 = 0; iter2 < totalConLength; iter2++) {
					trimming_j = matchVector[iter][(2*iter2)];
					trimming_i = matchVector[iter][(2*iter2)+1];
					ct->SetPair(trimming_j,trimming_i,i);
					
				}
			}

			// Check for helical extension possibilities
			HelicalExtension(ct, convertedDotplot, helixExtend, dgArray);
		}

		else{// If pseudoknot-free mode is specified
#if defined(VERBOSE_MODE)//output forsed paires
		cerr << "\nForced Pairs in Pseudoknot-Free mode:\n";
#endif
			// Loop over all base pairs in the match combination
			for (iter2 = 0; iter2 < totalConLength; iter2++) {
				if(matchVector[iter][(2*iter2)]<matchVector[iter][(2*iter2)+1]){//if the 2*iter2 position is a 5' position
					ct->AddPair(matchVector[iter][(2*iter2)],matchVector[iter][(2*iter2)+1]);
					
				}
				else{//else if the 2*iter2 position is a 3' position
					ct->AddPair(matchVector[iter][(2*iter2)],matchVector[iter][(2*iter2)+1]);
					
				}
#if defined(VERBOSE_MODE)//output forsed paires
				cerr << "pair\t" << matchVector[iter][(2*iter2)] << "-" << matchVector[iter][(2*iter2)+1] << endl;
#endif
				
			}
#if defined(VERBOSE_MODE)//output forsed paires
			cerr << "------------------------------------------\n";
#endif

			dynamic(ct,&data,maxtracebacks,99,windowsize);
		}

		// Copy current ct to new ct file that is used for appending the results from all matches
		for (i = 1; i <= ct->GetNumberofStructures(); i++) {
			//ct2->checknumberofstructures();
			ct2->AddStructure();
			ctbroken->AddStructure();
			
						
			for (j = 1; j <= ct->GetSequenceLength(); j++) {
				if (ct->GetPair(j,i)>j) ct2->SetPair(j,ct->GetPair(j,i),i+start);
				
			}
			ct2->SetCtLabel(ct->GetCtLabel(i),i+start);
			
		}
		start += ct->GetNumberofStructures();
	}
	
	cout << "\t\t\tDONE\n" << flush;

#if defined(VERBOSE_MODE)
	cout << "Refolding yields " << ct2->numofstructures << " total structures.\n\n" << flush;
#endif
	//QUESTIONS QUESTIONS:	
	delete[] helixExtend;

	rnaCT2->WriteCt("ct2.ct");
	rnaCTBROKEN->WriteCt("ctBROKEN.ct");

	// Re-calculate free energy with modified efn2
	for (i=1;i<=ct2->GetNumberofStructures();i++){
		
		//if(i % 100 == 0 && i < 11440||i==1){
			cout << "\rRecalculating energies " << i << " of " << ct2->GetNumberofStructures() << flush;
			//}
			//		else if(i>=11440){
			//			cout << "\rRecalculating energies " << i << " of " << ct2->numofstructures << flush;
			//		}
		
			//rnaCTBROKEN->WriteCt("ctbroken");

		//cerr << "\n" << rnaCT2->GetFreeEnergy(i) << " ###############\n";
		
		efn2mod(&data, ct2, i, false, ctbroken);

#if defined(VERBOSE_MODE)
		if (i % 100 == 0 && i > 0) cout << i << " free energies calculated so far.\n";
#endif

		// Re-insert pseudoknot stems that were broken for energy calculation
		//for (j=1; j<=ct2->GetSequenceLength(); j++){
		//	if (ctbroken->basepr[0][j]!=0) ct2->basepr[i][j] = ctbroken->basepr[0][j];
		//}
		//This is now moved inside efn2mod, and is not needed here.
	}
	cout << "\t\tDONE\nSorting structures..." << flush;
	// Re-sort the structures according to efn2mod-calculated energy
	ct2->sort();
	
	cout << "\t\t\t\tDONE\n" << flush;

	// Calculate value of cutoff for output structures
	if (cutoff !=0) cutoff = (100-cutoff)*ct2->GetEnergy(1)/100;

	// Re-output ct file
	cout << ctout2(ct2, cutoff, outct.c_str()) << " structures are within the specified cutoff percentage.\n";

	//  Optional: output paired-positions text file (for structure viewing in PseudoViewer3)
	//	Note: this outputs one concatenated file, individual structures must be manually cut from this file
	//  and placed in a new text file for PseudoViewer3 to read it.
	if(outpairs != "") pairout(ct2, cutoff, outpairs.c_str());

	// Insert steps to clean up allocated memory
	matchVector.clear();
	delete[] alreadyUsed;
	delete[] isSymmetric;
	
	for (i=0;i<=ct->GetSequenceLength();i++){
		delete[]convertedDotplot[i];
		delete[]dgArray[i];
		delete[]mbDotplot[i];
	}
	delete[]convertedDotplot;
	delete[]dgArray;
	delete[]mbDotplot;

	for(i=0;i<=numOfCon;i++){
		for(j=0;j<=maxConLength;j++){
			for(int k=0;k<2;k++){
				delete[]tripletArray[i][j][k];
			}
		}
	}
 	for(i=0;i<=numOfCon;i++){
		for(j=0;j<=maxConLength;j++){
			delete[]tripletArray[i][j];
		}
	}
	for(i=0;i<=numOfCon;i++){
		delete[]tripletArray[i];
	}
	delete[]tripletArray;

	for (i=0;i<=numOfCon;i++){
		delete[]conArray[i];
		delete[]xCoords[i];
		delete[]yCoords[i];
	}
	delete[]conArray;
	delete[]xCoords;
	delete[]yCoords;
	delete rnaCT;
	delete rnaCT2;
	delete rnaCTBROKEN;

	cout << "\n####################################################\n";
	cout << "                        DONE\n";
	cout << "####################################################\n" << flush;

	return 0;
}

void getdat(char *loop, char *stackf, char *tstackh, char *tstacki,
            char *tloop, char *miscloop, char *danglef, char *int22,
            char *int21,char *coax, char *tstackcoax,
            char *coaxstack, char *tstack, char *tstackm, char *triloop,
            char *int11, char *hexaloop, char *tstacki23, char *tstacki1n,
            char *datapath, bool isRNA)
{
	if( !isRNA) strcat( datapath,"dna");
	strcpy (loop,datapath);
	strcpy (stackf,datapath);
	strcpy (tstackh,datapath);
	strcpy (tstacki,datapath);
	strcpy (tloop,datapath);
	strcpy (miscloop,datapath);
	strcpy (danglef,datapath);
	strcpy (int22,datapath);
	strcpy (int21,datapath);
	strcpy (triloop,datapath);
	strcpy (coax,datapath);
	strcpy (tstackcoax,datapath);
	strcpy (coaxstack,datapath);
	strcpy (tstack,datapath);
	strcpy (tstackm,datapath);
	strcpy (int11,datapath);
	strcpy (hexaloop,datapath);
	strcpy (tstacki23,datapath);
	strcpy (tstacki1n,datapath);

	strcat (loop,"loop.dat");
	strcat (stackf,"stack.dat");
	strcat (tstackh,"tstackh.dat");
	strcat (tstacki,"tstacki.dat");
	strcat (tloop,"tloop.dat");
	strcat (miscloop,"miscloop.dat");
	strcat (danglef,"dangle.dat");
	strcat (int22,"int22.dat");
	strcat (int21,"int21.dat");
	strcat (triloop,"triloop.dat");
	strcat (coax,"coaxial.dat");
	strcat (tstackcoax,"tstackcoax.dat");
	strcat (coaxstack,"coaxstack.dat");
	strcat (tstack,"tstack.dat");
	strcat (tstackm,"tstackm.dat");
	strcat (int11,"int11.dat");
	strcat (hexaloop,"hexaloop.dat");
	strcat (tstacki23,"tstacki23.dat");
	strcat (tstacki1n,"tstacki1n.dat");
}

void errmsg2(int err,int erri) {

	if (err==30) {
		cout << "End Reached at traceback #"<<erri<<"\n";
		exit(1);
	}
	if (err==100) {
		cout << "error # "<<erri;
		exit(1);
	}
	switch (err) {
	case 1:
		cout << "Could not allocate enough memory";
		break;
	case 2:
		cout << "Too many possible base pairs";
		break;
	case 3:
		cout << "Too many helixes in multibranch loop";
	case 4:
		cout << "Too many structures in CT file";
	default:
		cout << "Unknown error";
	}
	cin >> err;
	exit(1);
	return;
}


void firstDotMatch(short** conArray, bool* alreadyUsed, bool* isSymmetric, short totalConLength,
                   short numOfCon, short** xCoords, short** yCoords, vector<conMatch>* matchVector,
                   short* pCurrConNum, short* pCurrConPos, short** convertedDotplot, short** mbDotplot,
                   int dgCutoff, short sequenceLength,structure *ct,char**** tripletArray) {

	short i, j;
	short currConNum = *(pCurrConNum);
	short currConPos = *(pCurrConPos);

#if defined(DEBUG_MODE)
	cout << "currConNum="<< currConNum << " currConPos=" << currConPos << "\n";
	cin >> i;
#endif

	if(matchVector->size()>warningLimit && ifwarningMessage==false){
		cout << "\r                                                   ";
		cout << "\n!! There are too many matches. This will make the program run very long time.\n"<< flush;
		cout << "!! To speed up the calculation try lowering the -d and -m parameters (see help),\n"<< flush;
		cout << "!! or try removing shortest constraints from the constraint file.\n\n" << flush;
		ifwarningMessage=true;
	}
	else if(matchVector->size()%1000 == 0){
		cout << "\rLooking for matches: " << matchVector->size() << " matches found" << flush;
	}

	// Loop through x,y of dotplot to find potential match for first position of the current constraint
	// Skip all values that have already been used by matches for previous constraints
	for (j = 1; j < sequenceLength+1; j++) {
		if (alreadyUsed[j] == true) continue;

		for (i = 1; i < sequenceLength+1; i++) {
			if (alreadyUsed[i] == true) continue;

			// If constraint is symmetric, only consider matches from upper-right half of dotplot
			if (isSymmetric[currConNum]==false || (isSymmetric[currConNum]==true && i<j)) {

				// Check the first value only here - if this value matches,
				// then call the recursive method, passing it the proper parameters
				if (conArray[currConNum][1] == convertedDotplot[i][j]) {
					xCoords[currConNum][1] = i;
					yCoords[currConNum][1] = j;
					alreadyUsed[i] = true;
					alreadyUsed[j] = true;
					*(pCurrConPos) = 2;
					recursiveMatch(conArray, alreadyUsed, totalConLength, numOfCon, xCoords, yCoords, 
					               matchVector, pCurrConNum, pCurrConPos, convertedDotplot, mbDotplot, dgCutoff, sequenceLength, 
					               isSymmetric,ct,tripletArray);
					currConNum = *(pCurrConNum);
					currConPos = *(pCurrConPos);
				}
			} 
		}
	}
	// Once we've looped through all coordinates for this base pair, step back values for constraint 
	// number and position, flip switches for last base pair, and return
	currConNum--;
	if (currConNum > -1) {
		currConPos = conArray[currConNum][0];
		alreadyUsed[xCoords[currConNum][currConPos]] = false;
		alreadyUsed[yCoords[currConNum][currConPos]] = false;
	}
	*(pCurrConNum) = currConNum;
	*(pCurrConPos) = currConPos;
	return;
}

void recursiveMatch(short** conArray, bool* alreadyUsed, short totalConLength, short numOfCon,
                    short** xCoords, short** yCoords, vector<conMatch>* matchVector, 
                    short* pCurrConNum, short* pCurrConPos, short** convertedDotplot, short** mbDotplot,
                    int dgCutoff, short sequenceLength, bool* isSymmetric,structure *ct,char**** tripletArray) {

	short i,j,k,xIter,yIter,start,stop;
	short currConNum = *(pCurrConNum);
	short currConPos = *(pCurrConPos);
	//conMatch* tempMatch;
	string a;
	bool searchX;  // Switch to indicate in which direction we're searching
	bool searchExtCoax; // Switch to indicate whether we're currently searching for a coaxial stack in 
	// the external direction
	
 subroutine: // Had to use this so that part of the code could return here from a nested "for" loop

	// This function will loop until we reach the end of the current constraint, 
	// or until no more matches are found using the current nucleotides for the
	// first base pair of this constraint
	while (currConPos > 1) {
		
#if defined(DEBUG_MODE)
		cout << currConNum << " " << currConPos << " " << xCoords[currConNum][currConPos-1] << 
			" " << yCoords[currConNum][currConPos-1] << " " << xCoords[currConNum][currConPos] <<
			" " << yCoords[currConNum][currConPos] << " First\n";
		cin >> i;
#endif

		// Check if we've matched all the base pairs for the current constraint
		if (currConPos == conArray[currConNum][0]+1) {
			// Check if this is the final constraint
			if (currConNum == numOfCon-1) {
				// If so, we've completed a match - compile 1D short array of dotplot 
				// coordinates and append this to matchVector

				//				conMatch tempMatch;
				/*
				tempMatch = new conMatch;
				tempMatch->coords = new short[totalConLength*2];
				*/
				//				tempMatch.coords=vector<short>(totalConLength*2);
				conMatch tempMatch(totalConLength*2);

				i = 0; // current base pair to be written
				j = 0; // takes place of currConNum
				k = 1; // takes place of currConPos
				
				while (i < totalConLength) {
					if (k > conArray[j][0]) {j++; k=1; continue;} // reached end of current constraint
					tempMatch[i*2] = xCoords[j][k];
					tempMatch[i*2+1] = yCoords[j][k];
					i++;
					k++;
				}

#if defined(VERBOSE_MODE)
				cout << "Match found!: ";
				for (i = 0; i < totalConLength*2; i++) {
					cout << tempMatch[i] << " ";
				}
				cout << "\n";
#endif

#if defined(ALREADY_USED_CHECK)
				// Verify that the number of "true" alreadyUsed values = twice the number of basepairs in all constraints
				j = 0;
				for (i = 1; i<=sequenceLength; i++){
					if (alreadyUsed[i]) j++;
				}
				if (j!=totalConLength*2){cerr<<"ERROR: Discrepancy in alreadyUsed array.\n";}
#endif

				matchVector->push_back(tempMatch);

#if defined(VERBOSE_MODE)
				if (matchVector->size() % 50000 == 0) {
					// Display a current count of matches
					cout<<matchVector->size()<<" ";
					for (i=0; i<numOfCon; i++){
						cout <<xCoords[i][1]<<","<<yCoords[i][1]<<" ";
					}
					cout << "\n";
				}
#endif
				// Now step back two base pairs, clear final position in the Coords arrays,
				// flip the alreadyUsed values of last 2 positions, and continue searching
				currConPos=currConPos-2;
				alreadyUsed[xCoords[currConNum][currConPos+1]] = false;
				alreadyUsed[yCoords[currConNum][currConPos+1]] = false;
				xCoords[currConNum][currConPos+1] = 0;
				yCoords[currConNum][currConPos+1] = 0;
				alreadyUsed[xCoords[currConNum][currConPos]] = false;
				alreadyUsed[yCoords[currConNum][currConPos]] = false;
				continue;
			}

			// If this is not the last constraint, we need to move on and match the first position 
			// of the next constraint
			currConNum++;
			currConPos = 1;
			*(pCurrConNum) = currConNum;
			*(pCurrConPos) = currConPos;
			firstDotMatch(conArray,alreadyUsed,isSymmetric,totalConLength,numOfCon,xCoords,yCoords,matchVector,
			              pCurrConNum,pCurrConPos,convertedDotplot,mbDotplot,dgCutoff,sequenceLength,ct,tripletArray);
			// Once this returns, need to continue the search for matches to this constraint
			currConNum = *(pCurrConNum);
			currConPos = *(pCurrConPos);
			continue;
		}
		
		// Check if we're still in-bounds on the next position across and down.
		// If we're not, step back, flip switches, and continue
		if (xCoords[currConNum][currConPos-1]-1 < 1 || yCoords[currConNum][currConPos-1]+1 > sequenceLength) {
			currConPos--;
			alreadyUsed[xCoords[currConNum][currConPos]] = false;
			alreadyUsed[yCoords[currConNum][currConPos]] = false;
			continue;
		}
		// If this is the first or last step, or if the previous basepair is unstable as a closing pair for 
		// external and mb loops, only look for non-bulged/non-stacked match on current step
		if (currConPos==2 || currConPos==(conArray[currConNum][0]) || 
		    mbDotplot[xCoords[currConNum][currConPos-1]][yCoords[currConNum][currConPos-1]] >= dgCutoff /*||
			                            xCoords[currConNum][currConPos-2]-xCoords[currConNum][currConPos-1]>1 ||
			                            yCoords[currConNum][currConPos-1]-yCoords[currConNum][currConPos-2]>1*/) {
			if (xCoords[currConNum][currConPos] != 0 || yCoords[currConNum][currConPos] !=0) {
				// This would indicate that this basepair has already been checked
				// Reset coords, step back, flip switches, and continue
				xCoords[currConNum][currConPos] = 0;
				yCoords[currConNum][currConPos] = 0;
				currConPos--;
				alreadyUsed[xCoords[currConNum][currConPos]] = false;
				alreadyUsed[yCoords[currConNum][currConPos]] = false;
				continue;
			}

			if (conArray[currConNum][currConPos]==convertedDotplot[xCoords[currConNum][currConPos-1]-1]
			    [yCoords[currConNum][currConPos-1]+1] 
			    && alreadyUsed[xCoords[currConNum][currConPos-1]-1] == false
			    && alreadyUsed[yCoords[currConNum][currConPos-1]+1] == false) {

				// Match found at lastX-1,lastY+1
				xCoords[currConNum][currConPos] = xCoords[currConNum][currConPos-1]-1;
				yCoords[currConNum][currConPos] = yCoords[currConNum][currConPos-1]+1;

#if defined(TRIPLET_VERBOSE_MODE)
				cerr << "Location:A Constraint_Number:"<< currConNum << " Constraints_From/To:" << conArray[currConNum][currConPos-1] << "->" << conArray[currConNum][currConPos] 
				     << " Constraint_Position_From/To:" << currConPos-1 << "->" << currConPos 
				     << " Seqence_Position_From/To: xCoords[" << xCoords[currConNum][currConPos-1]<< "->" << xCoords[currConNum][currConPos] 
				     << "] Ycoord[" << yCoords[currConNum][currConPos-1]<< "->" << yCoords[currConNum][currConPos] << "]\n";
#endif
				//Run the 'TripletMatch' function. This function test if there is a triplet constraint for the current location on the sequence.
				//It also checks if there are any bulges or loops between the last 3 matched positions.
				//If there is a triplet constraint and there are no loops or bulges separating the constraint, it checks if the triplet
				//constraint matches the current matched helix, and returns 'true' for a MATCH, and 'false' for NON-MATCH.
				//Otherwise it by default returns 'true'. 
				if(TripletMatch(ct,tripletArray,xCoords,yCoords,currConNum,currConPos)){
					alreadyUsed[xCoords[currConNum][currConPos-1]-1] = true;
					alreadyUsed[yCoords[currConNum][currConPos-1]+1] = true;
					currConPos++;
					continue;
				}
				//If the triplet constraint hasn't matched, return the current position xCoords and yCoords to previous state	
				else{
					// If not a match, step back, flip switches, and continue
					// Match found at lastX-1,lastY+1
					xCoords[currConNum][currConPos] = 0;
					yCoords[currConNum][currConPos] = 0;
					currConPos--;
					alreadyUsed[xCoords[currConNum][currConPos]] = false;
					alreadyUsed[yCoords[currConNum][currConPos]] = false;
					continue;
				}
			}	
			// If not a match, step back, flip switches, and continue
			currConPos--;
			alreadyUsed[xCoords[currConNum][currConPos]] = false;
			alreadyUsed[yCoords[currConNum][currConPos]] = false;
			continue;
		}

		// Next step is in-bounds - find out where the search left off
		if (xCoords[currConNum][currConPos] == 0 && yCoords[currConNum][currConPos] == 0){
			// This is the first time searching forward from this position
			// Set coords values to -1,-1 to act as placeholder
			xCoords[currConNum][currConPos] = -1;
			yCoords[currConNum][currConPos] = -1;
			searchX = true;
			searchExtCoax = false;
			start = xCoords[currConNum][currConPos-1] - 1;
		}
	
		else if (yCoords[currConNum][currConPos] == yCoords[currConNum][currConPos-1]+1 &&
		         xCoords[currConNum][currConPos] <= xCoords[currConNum][currConPos-1]-1) {
			// We left off searching in the X direction
			if (xCoords[currConNum][currConPos-1] -1 < 1){ 
				// Last found a match at left edge of dot plot - start searching in Y direction
				xCoords[currConNum][currConPos] = -1;
				yCoords[currConNum][currConPos] = -1;
				continue;
			}

			// Otherwise, continue searching in X direction
			searchX = true;
			searchExtCoax = false;
			start = xCoords[currConNum][currConPos] - 1;
		}

		else if (yCoords[currConNum][currConPos] == -1) {
			// Search in X direction was completed
			if (yCoords[currConNum][currConPos-1] + 2 > sequenceLength) {
				// Cannot search in Y direction because last match is near bottom edge of dot plot
				// Time to search for external coaxial stacks
				yCoords[currConNum][currConPos] = yCoords[currConNum][currConPos-1]+1;
				xCoords[currConNum][currConPos] = sequenceLength+1;
				continue;
			}
			
			// Otherwise, start at +2 because +1 was already covered in searchX
			searchX = false;
			searchExtCoax = false;
			start = yCoords[currConNum][currConPos-1] + 2;
		}

		else if (xCoords[currConNum][currConPos] == xCoords[currConNum][currConPos-1]-1 &&
		         yCoords[currConNum][currConPos] >= yCoords[currConNum][currConPos-1]+1) {
			// We left off searching in the Y direction
			if (yCoords[currConNum][currConPos] == sequenceLength){ 
				// Last found a match at bottom edge of dot plot - start searching for external coaxial stacks
				yCoords[currConNum][currConPos] = yCoords[currConNum][currConPos-1]+1;
				xCoords[currConNum][currConPos] = sequenceLength+1;
				continue;
			}

			// Otherwise, continue searching in Y direction
			searchX = false;
			searchExtCoax = false;
			start = yCoords[currConNum][currConPos] + 1;
		}

		else if (yCoords[currConNum][currConPos] == yCoords[currConNum][currConPos-1]+1 &&
		         xCoords[currConNum][currConPos] > yCoords[currConNum][currConPos]) {
			// We're currently searching for a external coaxial stack in the X direction
			searchX = true;
			searchExtCoax = true;
			start = xCoords[currConNum][currConPos] - 1;
		}

		else if (xCoords[currConNum][currConPos] == xCoords[currConNum][currConPos-1]-1 &&
		         yCoords[currConNum][currConPos] != -1 && 
		         yCoords[currConNum][currConPos] < xCoords[currConNum][currConPos]) {
			// We're currently searching for a external coaxial stack in the Y direction
			searchX = false;
			searchExtCoax = true;
			start = yCoords[currConNum][currConPos] + 1;
		}


		else {
			cerr<<"ERROR: Unclassified looping condition in constraint matching.\n"
			    <<xCoords[currConNum][currConPos]<<" "<<yCoords[currConNum][currConPos]<<"\n";			
			cin>>i;break;
		}
		// Enter this loop if we need to search in the X direction
		if (searchX && !searchExtCoax) {
			if (alreadyUsed[yCoords[currConNum][currConPos-1]+1]) {
				// Y coordinate has already been used, skip searching in X direction
				xCoords[currConNum][currConPos] = -1;
				yCoords[currConNum][currConPos] = -1;
				continue;
			}
			
			// If we're in the upper-right half of the dot plot...
			// Search the entire subspace of dotplot where next coordinate is 
			// x = {start,...,lastY + 1}, y = lastY+ 1
			if (xCoords[currConNum][currConPos-1] > yCoords[currConNum][currConPos-1]) {
				stop = (yCoords[currConNum][currConPos-1]);
			}

			// If we're in the lower-left half of the dot plot...
			// Search the entire subspace of dotplot where next coordinate is 
			// x = {start,...,1}, y = lastY + 1
			else stop = 1;

#if defined(DEBUG_MODE)
			cout << "start="<< start << " stop=" << stop << " Second\n";
			cin >> i;
#endif
			
			for (xIter=0;(start-xIter) >= stop; xIter++) {
				// Don't include certain values of xCoord
				if (alreadyUsed[start-xIter] ||
				    (xCoords[currConNum][currConPos-1]-(start-xIter) > 2 &&
				     xCoords[currConNum][currConPos-1]-(start-xIter) < 8)) continue;

#if defined(DEBUG_MODE)
				cout << "start-xIter=" << start-xIter << " Third";
#endif

				if (conArray[currConNum][currConPos]==convertedDotplot[start-xIter][yCoords[currConNum][currConPos-1]+1] && 
				    (start-xIter == xCoords[currConNum][currConPos-1] - 1 ||
				     mbDotplot[start-xIter][yCoords[currConNum][currConPos-1]+1] <= dgCutoff)) {
					// Match found: store coords, flip switches, step forward, and continue
					xCoords[currConNum][currConPos] = start-xIter;
					yCoords[currConNum][currConPos] = yCoords[currConNum][currConPos-1]+1;

#if defined(TRIPLET_VERBOSE_MODE)
					cerr << "Location:B Constraint_Number:"<< currConNum << " Constraints_From/To:" << conArray[currConNum][currConPos-1] << "->" << conArray[currConNum][currConPos] 
					     << " Constraint_Position_From/To:" << currConPos-1 << "->" << currConPos 
					     << " Seqence_Position_From/To: xCoords[" << xCoords[currConNum][currConPos-1]<< "->" << xCoords[currConNum][currConPos] 
					     << "] Ycoord[" << yCoords[currConNum][currConPos-1]<< "->" << yCoords[currConNum][currConPos] << "]\n";
#endif
					//Run the 'TripletMatch' function. This function test if there is a triplet constraint for the current location on the sequence.
					//It also checks if there are any bulges or loops between the last 3 matched positions.
					//If there is a triplet constraint and there are no loops or bulges separating the constraint, it checks if the triplet
					//constraint matches the current matched helix, and returns 'true' for a MATCH, and 'false' for NON-MATCH.
					//Otherwise it by default returns 'true'. 
					
					if(TripletMatch(ct,tripletArray,xCoords,yCoords,currConNum,currConPos)){
						alreadyUsed[start-xIter] = true;
						alreadyUsed[yCoords[currConNum][currConPos-1]+1] = true;
						currConPos++;
						// Can't simply continue the "while" loop because we're in a nested "for" loop
						goto subroutine;
					}
				}
			}

			// If we reach this point, we've iterated through all X values at this Y coordinate
			// Time to continue the search in the Y direction
			xCoords[currConNum][currConPos] = -1;
			yCoords[currConNum][currConPos] = -1;
			continue;
		}

		// Enter this loop if we need to search in the y direction 
		else if (!searchX && !searchExtCoax) {
			if (alreadyUsed[xCoords[currConNum][currConPos-1]-1]) {
				// X coordinate has already been used
				
				// Check if we're in the lower-left half of the dot plot - if so, need to search for external coaxial stack
				if (xCoords[currConNum][currConPos-1] < yCoords[currConNum][currConPos-1]) {
					yCoords[currConNum][currConPos] = yCoords[currConNum][currConPos-1] + 1;
					xCoords[currConNum][currConPos] = sequenceLength + 1;
					continue;
				}
				
				// We're in the upper-right half of the dot plot - reset the coords, step back, flip switches, and continue
				xCoords[currConNum][currConPos] = 0;
				yCoords[currConNum][currConPos] = 0;
				currConPos--;
				alreadyUsed[xCoords[currConNum][currConPos]] = false;
				alreadyUsed[yCoords[currConNum][currConPos]] = false;
				continue;
			}

			// If we're in the upper-right half of the dot plot...
			// Search the entire subspace of dotplot where next coordinate is 
			// x = lastX - 1, y = {start,...,lastX - 1}
			if (xCoords[currConNum][currConPos-1] > yCoords[currConNum][currConPos-1]) {
				stop = (xCoords[currConNum][currConPos-1]-1);
			}

			// If we're in the lower-left half of the dot plot...
			// Search the entire subspace of dotplot where next coordinate is 
			// x = lastX - 1, y = {start,...,sequenceLength}
			else stop = sequenceLength;

#if defined(DEBUG_MODE)
			cout << start << " " << stop << " Fourth\n";
			cin >> i;
#endif

			for (yIter=0;(start+yIter) <= stop; yIter++) {
				// Don't include certain values of yCoord
				if (alreadyUsed[start+yIter] ||
				    ((start+yIter)-yCoords[currConNum][currConPos-1] > 2 &&
				     (start+yIter)-yCoords[currConNum][currConPos-1] < 8)) continue;

#if defined(DEBUG_MODE)
				cout << "start+yIter=" << start+yIter << " Fifth";
#endif

				if (conArray[currConNum][currConPos]==
				    convertedDotplot[xCoords[currConNum][currConPos-1]-1][start+yIter] &&
				    mbDotplot[xCoords[currConNum][currConPos-1]-1][start+yIter] <= dgCutoff) {
					
					// Match found: store coords, flip switches, step forward, and continue
					xCoords[currConNum][currConPos] = xCoords[currConNum][currConPos-1]-1;
					yCoords[currConNum][currConPos] = start+yIter;

#if defined(TRIPLET_VERBOSE_MODE)
					cerr << "Location:C Constraint_Number:"<< currConNum << " Constraints_From/To:" << conArray[currConNum][currConPos-1] << "->" << conArray[currConNum][currConPos] 
					     << " Constraint_Position_From/To:" << currConPos-1 << "->" << currConPos 
					     << " Seqence_Position_From/To: xCoords[" << xCoords[currConNum][currConPos-1]<< "->" << xCoords[currConNum][currConPos] 
					     << "] Ycoord[" << yCoords[currConNum][currConPos-1]<< "->" << yCoords[currConNum][currConPos] << "]\n";
#endif
					//Run the 'TripletMatch' function. This function test if there is a triplet constraint for the current location on the sequence.
					//It also checks if there are any bulges or loops between the last 3 matched positions.
					//If there is a triplet constraint and there are no loops or bulges separating the constraint, it checks if the triplet
					//constraint matches the current matched helix, and returns 'true' for a MATCH, and 'false' for NON-MATCH.
					//Otherwise it by default returns 'true'. 
					
					if(TripletMatch(ct,tripletArray,xCoords,yCoords,currConNum,currConPos)){
						alreadyUsed[xCoords[currConNum][currConPos-1]-1] = true;
						alreadyUsed[start+yIter] = true;
						currConPos++;
						// Can't simply continue the "while" loop because we're in a nested "for" loop
						goto subroutine;
					}
				}

			}

			// If we reach this point, we've iterated through all x and y values for this base pair

			// Check if we're in the lower-left half of the dot plot - if so, need to search for external coaxial stack
			if (xCoords[currConNum][currConPos-1] < yCoords[currConNum][currConPos-1]) {
				yCoords[currConNum][currConPos] = yCoords[currConNum][currConPos-1] + 1;
				xCoords[currConNum][currConPos] = sequenceLength + 1;
				continue;
			}

			// We're in the upper-right half of the dot plot - reset the coords, step back, flip switches, and continue
			xCoords[currConNum][currConPos] = 0;
			yCoords[currConNum][currConPos] = 0;
			currConPos--;
			alreadyUsed[xCoords[currConNum][currConPos]] = false;
			alreadyUsed[yCoords[currConNum][currConPos]] = false;
			continue;
		}

		// Enter this loop if we're searching for an external coaxial stack in the X direction
		else if (searchX && searchExtCoax) {
			if (alreadyUsed[yCoords[currConNum][currConPos-1]+1]) {
				// Y coordinate has already been used, skip searching in X direction
				xCoords[currConNum][currConPos] = xCoords[currConNum][currConPos-1]-1;
				yCoords[currConNum][currConPos] = 0;
				continue;
			}
			
			stop = (yCoords[currConNum][currConPos-1]+1);
			
#if defined(DEBUG_MODE)
			cout << "start=" << start << " stop=" << stop << " Sixth\n";
			cin >> i;
#endif
			
			for (xIter=0;(start-xIter) >= stop; xIter++) {
				// Don't include certain values of xCoord
				if (alreadyUsed[start-xIter]) continue;

#if defined(DEBUG_MODE)
				cout << "start-xIter=" << start-xIter << " Seventh";
#endif

				if (conArray[currConNum][currConPos]==
				    convertedDotplot[start-xIter][yCoords[currConNum][currConPos-1]+1] && 
				    mbDotplot[start-xIter][yCoords[currConNum][currConPos-1]+1] <= dgCutoff) {
					// Match found: store coords, flip switches, step forward, and continue
					xCoords[currConNum][currConPos] = start-xIter;
					yCoords[currConNum][currConPos] = yCoords[currConNum][currConPos-1]+1;

#if defined(TRIPLET_VERBOSE_MODE)
					cerr << "Location:D Constraint_Number:"<< currConNum << " Constraints_From/To:" << conArray[currConNum][currConPos-1] << "->" << conArray[currConNum][currConPos] 
					     << " Constraint_Position_From/To:" << currConPos-1 << "->" << currConPos 
					     << " Seqence_Position_From/To: xCoords[" << xCoords[currConNum][currConPos-1]<< "->" << xCoords[currConNum][currConPos] 
					     << "] Ycoord[" << yCoords[currConNum][currConPos-1]<< "->" << yCoords[currConNum][currConPos] << "]\n";
#endif
					//Run the 'TripletMatch' function. This function test if there is a triplet constraint for the current location on the sequence.
					//It also checks if there are any bulges or loops between the last 3 matched positions.
					//If there is a triplet constraint and there are no loops or bulges separating the constraint, it checks if the triplet
					//constraint matches the current matched helix, and returns 'true' for a MATCH, and 'false' for NON-MATCH.
					//Otherwise it by default returns 'true'. 
					
					if(TripletMatch(ct,tripletArray,xCoords,yCoords,currConNum,currConPos)){
						alreadyUsed[start-xIter] = true;
						alreadyUsed[yCoords[currConNum][currConPos-1]+1] = true;
						currConPos++;
						// Can't simply continue the "while" loop because we're in a nested "for" loop
						goto subroutine;
					}
				}
			}

			// If we reach this point, we've iterated through all x values at this y coordinate
			// Time to continue the external coaxial stack search in the y direction
			xCoords[currConNum][currConPos] = xCoords[currConNum][currConPos-1]-1;
			yCoords[currConNum][currConPos] = 0;
			continue;
		}

		// Enter this loop if we need to search for external coaxial stacks in the y direction 
		else if (!searchX && searchExtCoax) {

			if (alreadyUsed[xCoords[currConNum][currConPos-1]-1]) {
				// X coordinate has already been used, time to reset, step back, flip, and continue 
				xCoords[currConNum][currConPos] = 0;
				yCoords[currConNum][currConPos] = 0;
				currConPos--;
				alreadyUsed[xCoords[currConNum][currConPos]] = false;
				alreadyUsed[yCoords[currConNum][currConPos]] = false;
				continue;
			}
			
			stop = (xCoords[currConNum][currConPos-1]-1);

#if defined(DEBUG_MODE)
			cout << "start=" << start << " stop=" << stop << " Eighths\n";
			cin >> i;
#endif

			for (yIter=0;(start+yIter) <= stop; yIter++) {
				// Don't include certain values of yCoord
				if (alreadyUsed[start+yIter]) continue;

#if defined(DEBUG_MODE)
				cout << "start+yIter="<< start+yIter << " Nineth";
#endif

				if (conArray[currConNum][currConPos]==
				    convertedDotplot[xCoords[currConNum][currConPos-1]-1][start+yIter] &&
				    mbDotplot[xCoords[currConNum][currConPos-1]-1][start+yIter] <= dgCutoff) {

					// Match found: store coords, flip switches, step forward, and continue
					xCoords[currConNum][currConPos] = xCoords[currConNum][currConPos-1]-1;
					yCoords[currConNum][currConPos] = start+yIter;

#if defined(TRIPLET_VERBOSE_MODE)
					cerr << "Location:E Constraint_Number:"<< currConNum << " Constraints_From/To:" << conArray[currConNum][currConPos-1] << "->" << conArray[currConNum][currConPos] 
					     << " Constraint_Position_From/To:" << currConPos-1 << "->" << currConPos 
					     << " Seqence_Position_From/To: xCoords[" << xCoords[currConNum][currConPos-1]<< "->" << xCoords[currConNum][currConPos] 
					     << "] Ycoord[" << yCoords[currConNum][currConPos-1]<< "->" << yCoords[currConNum][currConPos] << "]\n";
#endif
					//Run the 'TripletMatch' function. This function test if there is a triplet constraint for the current location on the sequence.
					//It also checks if there are any bulges or loops between the last 3 matched positions.
					//If there is a triplet constraint and there are no loops or bulges separating the constraint, it checks if the triplet
					//constraint matches the current matched helix, and returns 'true' for a MATCH, and 'false' for NON-MATCH.
					//Otherwise it by default returns 'true'. 
					
					if(TripletMatch(ct,tripletArray,xCoords,yCoords,currConNum,currConPos)){
						alreadyUsed[xCoords[currConNum][currConPos-1]-1] = true;
						alreadyUsed[start+yIter] = true;
						currConPos++;
						// Can't simply continue the "while" loop because we're in a nested "for" loop
						goto subroutine;
					}
				}
			}

			// If we reach this point, we've iterated through all x and y values for this base pair
			// Time to reset the coords, step back, flip switches, and continue
			xCoords[currConNum][currConPos] = 0;
			yCoords[currConNum][currConPos] = 0;
			currConPos--;
			alreadyUsed[xCoords[currConNum][currConPos]] = false;
			alreadyUsed[yCoords[currConNum][currConPos]] = false;
			continue;
		}
	}
	// The current position has gone back to the first base pair of a constraint
	// Flip switches and return
	alreadyUsed[xCoords[currConNum][currConPos]] = false;
	alreadyUsed[yCoords[currConNum][currConPos]] = false;
	*(pCurrConNum) = currConNum;
	*(pCurrConPos) = currConPos;
	return;
}

// Search dotplot for complicated pseudoknot folds to exclude
void pseudodptrim(structure *ct, int *tempbp, int * count) {
	int i, iter, iter2, iter3, init=1, a=0, b=0, c=0, d=0, e=0, f=0, j=0;
	bool rightSingle, findCF, leftSingle, findAD;
	while(init <= ct->GetSequenceLength()) { 
		if(tempbp[init] > init) {
			for(i = 1; i < init; i++) {
				if(tempbp[i] > init && tempbp[i] < tempbp[init]) {
					// Define stems of pseudoknotted region
					j = tempbp[init];
					e = tempbp[i];
					b = tempbp[j];

					rightSingle = false;
					findCF = false;
					f = e+1; // Initial guess
					while( findCF == false && rightSingle == false) {
						if( tempbp[f] != 0 && tempbp[f] < e) {
							c = tempbp[f];
							findCF = true;
						}
						else {
							f++;
							if( f >= j) {
								rightSingle = true;
								f = j;
								c = b;
							}
						}
					}

					leftSingle = false;
					findAD = false;
					a = b - 1; // Initial guess
					while( findAD == false && leftSingle == false) {
						if( tempbp[a] != 0 && tempbp[a] >  c) {
							d = tempbp[a];			
							findAD = true;
						}
						else {
							a--;
							if( a <= i) {
								leftSingle = true;
								a = i;
								d = e;
							}
						}
					}
					
#if defined(VERBOSE_MODE)
					cout << "Pseudoknot detected: " << 
						i << "-" << e << "," << a << "-" << d << " ; " <<
						c << "-" << f << "," << b << "-" << j << "\n";
#endif


					// Remove all possible base pairs that would involved pairing between the 2 gap regions
					int count = 0;
					for (iter = a+1; iter < b; iter++) {
						for (iter2 = e+1; iter2 < f; iter2++) {
							if (ct->tem[iter2][iter] == true) count++;
							ct->tem[iter2][iter] = false;
						}
					}

#if defined(VERBOSE_MODE)
					if (count > 0) cout << "Excluded pairings between gaps: " << count << "\n";
#endif
				
					// Scan gaps for contained secondary structures - forbid this region from pairing outside the gap
					// a-b gap:
					count = 0;
					iter = a+1;
					while (iter < b ) {
						if (tempbp[iter] > iter && tempbp[iter] > a && tempbp[iter] < b) { 
							// Gap has self-contained secondary structure
							// First condition is to prevent endless looping scenario
							// Forbid all pairings for iter-tempbp[iter] to anything outside the gap
							for (iter2 = iter; iter2 <= tempbp[iter]; iter2++) {
								for (iter3 = 1; iter3 <= a; iter3++) {
									if (ct->tem[iter2][iter3] == true) count++;
									ct->tem[iter2][iter3] = false;
								}
								for (iter3 = b; iter3 <= ct->GetSequenceLength(); iter3++) {
									if (ct->tem[iter3][iter2] == true) count++;
									ct->tem[iter3][iter2] = false;
								}
							}
							// Continue loop at end of contained secondary structure
							iter = (tempbp[iter] + 1);
							continue;
						}
						iter++;
						continue;
					}
#if defined(VERBOSE_MODE)
					if (count > 0) cout << "Pairings forbidden from structures in a-b gap: " << count << "\n";
#endif
					
					// c-d gap:
					count = 0;
					iter = c+1;
					while (iter < d ) {
						if (tempbp[iter] > iter && tempbp[iter] > c && tempbp[iter] < d) { 
							// Gap has self-contained secondary structure
							// First condition is to prevent endless looping scenario
							// Forbid all pairings for iter-tempbp[iter] to anything outside the gap
							for (iter2 = iter; iter2 <= tempbp[iter]; iter2++) {
								for (iter3 = 1; iter3 <= c; iter3++) {
									if (ct->tem[iter2][iter3] == true) count++;
									ct->tem[iter2][iter3] = false;
								}
								for (iter3 = d; iter3 <= ct->GetSequenceLength(); iter3++) {
									if (ct->tem[iter3][iter2] == true) count++;
									ct->tem[iter3][iter2] = false;
								}
							}
							// Continue loop at end of contained secondary structure
							iter = (tempbp[iter] + 1);
							continue;
						}
						iter++;
						continue;
					}
#if defined(VERBOSE_MODE)
					if (count > 0) cout << "Pairings forbidden from structures in c-d gap: " << count << "\n";
#endif
				
					// e-f gap:
					count = 0;
					iter = e+1;
					while (iter < f ) {
						if (tempbp[iter] > iter && tempbp[iter] > e && tempbp[iter] < f) { 
							// Gap has self-contained secondary structure
							// First condition is to prevent endless looping scenario
							// Forbid all pairings for iter-tempbp[iter] to anything outside the gap
							for (iter2 = iter; iter2 <= tempbp[iter]; iter2++) {
								for (iter3 = 1; iter3 <= e; iter3++) {
									if (ct->tem[iter2][iter3] == true) count++;
									ct->tem[iter2][iter3] = false;
								}
								for (iter3 = f; iter3 <= ct->GetSequenceLength(); iter3++) {
									if (ct->tem[iter3][iter2] == true) count++;
									ct->tem[iter3][iter2] = false;
								}
							}
							// Continue loop at end of contained secondary structure
							iter = (tempbp[iter] + 1);
							continue;
						}
						iter++;
						continue;
					}
					
#if defined(VERBOSE_MODE)
					if (count > 0) cout << "Pairings forbidden from structures in e-f gap: " << count << "\n";
#endif
					
					// Continue loop at 3' end of pseudoknot
					init = j;
					continue;
				}
			}				
		}
		// If no pseudoknot detected, step up one position and recheck
		init++;
		continue;
	}
	return;
}
/*	Function efn2mod - pseudoknot-capable version of efn2

	Calculates the free energy of each structure in a structure called structure.

	Structures can contain certain types of pseudoknots (those allowed in the Dirks and Pierce algorithm)

	Structnum indicates which structure to calculate the free energy
	the default, 0, indicates "all" structures

	simplemb = true means that the energy function for mb loops is linear,
	so identical to the dynamic programming algorithms.
	The default is false, so logarithmic in the number of unpaired nucleotides.

	This new version of efn2 is now stacking aware -- if ct->stacking==true, it
	will use the stacking information stored in ct, rather than inferring stacking
	by free energy minimization.
*/
void efn2mod(datatable *data,structure *ct, int structnum, bool simplemb, structure *ctbroken) {
	int i,j,k,open,null,stz,count,sum,ip,jp;
	stackstruct stack;
	forceclass fce(ct->GetSequenceLength());
	int start, stop;

	ofstream out;

	//build an array to store broken base pairs"
	int *brokenpairs;
	brokenpairs = new int[ctbroken->GetSequenceLength()+1];

#if defined(debugmode)
	char filename[maxfil],temp[maxfil];
	int tempsum;
#endif
	
	/*	stack = a place to keep track of where efn2 is located in a structure
		inter = indicates whether there is an intermolecular interaction involved in
		a multi-branch loop
	*/

	stack.sp = 0;  //set stack counter
	
	if (ct->intermolecular) {//this indicates an intermolecular folding
		for (i=0;i<3;i++) {
			forceinterefn(ct->inter[i],ct,&fce);
		}
	}
	
	if (structnum!=0) {
		start = structnum;
		stop = structnum;
	}
	else {
		start = 1;
		stop = ct->GetNumberofStructures();
	}
	
	for (count=start;count<=stop;count++){//one structure at a time

#if defined(debugmode) //open an output file for debugging
		strcpy(filename,"efn2dump");
		itoa(count,temp,10);
		strcat(filename,temp);
		strcat(filename,".out");
		out.open(filename);
#endif

		// Reset the zero-order ctbroken basepr array for storing pairing information of all broken pseudoknot stems
		for (j = 1; j <=ctbroken->GetSequenceLength(); j++){
			brokenpairs[j] = 0;
		}


		// Scan the structure to determine if pseudoknot(s) are present
		bool nopseudos = false;
		int init = 1;
		ct->SetEnergy(count,0);
		while (!nopseudos) {
			ct->SetEnergy(count,ct->GetEnergy(count)+pseudoenergy(data, ct, count, &nopseudos, ctbroken, &init, brokenpairs));
		}

		// Original efn2 takes over from here
		ct->SetEnergy(count,ct->GetEnergy(count)+ergexterior(count, ct, data));
		

#if defined(debugmode)
		gcvt((float (ct->energy[count]))/conversionfactor,6,temp);
		out << "Exterior loop = "<<temp<<"\n";	
#endif
		
		i=1;
		while (i<ct->GetSequenceLength()) { 	
			if (ct->GetPair(i,count)!=0) {
				push(&stack,i,ct->GetPair(i,count),1,0);
				i = ct->GetPair(i,count);
			}
			i++;
		}

	subroutine://loop starts here (I chose the goto statement to make the code more readable)
		pull(&stack,&i,&j,&open,&null,&stz);//take a substructure off the stack

		while (stz!=1){
			while (ct->GetPair(i,count)==j) { //are i and j paired?

#if defined(debugmode) 
				tempsum=0;
#endif

				while (ct->GetPair(i+1,count)==j-1) {//are i,j and i+1,j-1 stacked?
					ct->SetEnergy(count,ct->GetEnergy(count) + erg1(i,j,i+1,j-1,ct,data));

#if defined(debugmode)
					tempsum = tempsum + erg1(i,j,i+1,j-1,ct,data); 
					gcvt((float (erg1(i,j,i+1,j-1,ct,data)))/conversionfactor,6,temp);
					out << "\tStack = "<<temp<<"  for stack of "<<i<<"-"<<j<<"\n";
#endif

					i++;
					j--;
				}

#if defined(debugmode)
				gcvt((float (tempsum))/conversionfactor,6,temp);
				out << "Helix total = "<<temp<<"\n";	
#endif

				sum = 0;
				k = i + 1;

				// now efn2 is past the paired region, so define the intervening non-paired segment
				while (k<j) {
					if (ct->GetPair(k,count)>k)	{
						sum++;
						ip = k;
						k = ct->GetPair(k,count) + 1;
						jp = k-1;
					}
					else if (ct->GetPair(k,count)==0) k++;
				}

				if (sum==0) {//hairpin loop
					ct->SetEnergy(count,ct->GetEnergy(count) + erg3(i,j,ct,data,fce.f(i,j-i)));

#if defined(debugmode)
					gcvt((float (erg3(i,j,ct,data,fce.f(i,j-i))))/conversionfactor,6,temp);
					out << "Hairpin = "<<temp<<"  for closure of "<<i<<"-"<<j<<"\n";	
#endif
					
					goto subroutine;
				}

				else if (sum==1) {//bulge/internal loop
					ct->SetEnergy(count, ct->GetEnergy(count) +	erg2(i,j,ip,jp,ct,data,fce.f(i,ip-i),fce.f(jp,j-jp)));

#if defined(debugmode)
					gcvt((float (erg2(i,j,ip,jp,ct,data,fce.f(i,ip-i),fce.f(jp,j-jp))))/conversionfactor,6,temp);
					out << "Internal/bulge = "<<temp<<"  for closure of "<<i<<"-"<<j<<"\n";	
#endif

					i = ip;
					j = jp;
				}

				else {//multi-branch loop
					ct->SetEnergy(count, ct->GetEnergy(count)+ ergmulti(count, i, ct, data, simplemb));
					
#if defined(debugmode)
					gcvt((float (ergmulti(count, i, ct, data)))/conversionfactor,6,temp);
					out << "Multiloop = "<<temp<<"  for closure of "<<i<<"-"<<j<<"\n";	
#endif

					//put the exiting helices on the stack:
					sum++;//total helixes = sum + 1
					i++;
					for (k=1;k<sum;k++) {
						while (ct->GetPair(i,count)==0) i++;
						push (&stack,i,ct->GetPair(i,count),1,0);
						i = ct->GetPair(i,count)+1;
					}
					
					goto subroutine;
				}
			}    
		}
		
#if defined(debugmode)
		gcvt((float (ct->energy[count]))/conversionfactor,6,temp);
		out << "\n\nTotal energy = "<<temp<<"\n";
		out.close();
#endif

		//return;-Removed by DHM - no need to return here
		//Restore broken base pairs:
		for (k=1;k<=ctbroken->GetSequenceLength();++k) {

			if (brokenpairs[k]>k) ctbroken->SetPair(k,brokenpairs[k],count);
		}
	}

	
	delete[] brokenpairs;
}//end efn2mod

// Scan structure ct[count] for pseudoknots, returning the energy for breaking pseudoknot stem(s) * 10
//Store any pairs broken in structure, count, in the brokenpairs array.
int pseudoenergy(datatable *data, structure *ct, int count, bool *nopseudo, structure *ctbroken, int *pinit, int *brokenpairs) {
	int i, iter, a=0, b=0, c=0, d=0, e=0, f=0, j=0;
	int stem1bps = 0, stem2bps = 0;
	int init = *(pinit);
	bool rightSingle, findCF, leftSingle, findAD;

	for(init; init <= ct->GetSequenceLength(); init++) {
		if(ct->GetPair(init,count) > init) {
			for(i = 1; i < init; i++) {
				if(ct->GetPair(i,count) > init && ct->GetPair(i,count) < ct->GetPair(init,count)) {
					// Update pointer so that we remember where the pseudoknot search left off
					*(pinit) = init;
					// Define stems of pseudoknotted region
					j = ct->GetPair(init,count);
					e = ct->GetPair(i,count);
					b = ct->GetPair(j,count);

					rightSingle = false;
					findCF = false;
					f = e+1; // Initial guess
					while( findCF == false && rightSingle == false) {
						if( ct->GetPair(f,count) != 0 && ct->GetPair(f,count) < e) {
							c = ct->GetPair(f,count);
							findCF = true;
						}
						else {
							f++;
							if( f >= j) {
								rightSingle = true;
								f = j;
								c= b;
							}
						}
					}

					leftSingle = false;
					findAD = false;
					a = b - 1; // Initial guess
					while( findAD == false && leftSingle == false) {
						if( ct->GetPair(a,count) != 0 && ct->GetPair(a,count) >  c) {
							d = ct->GetPair(a,count);			
							findAD = true;
						}
						else {
							a--;
							if( a <= i) {
								leftSingle = true;
								a = i;
								d = e;
							}
						}
					}
#if defined(VERBOSE_MODE)
					cout << "Pseudoknot detected in structure " << count << ": " << 
						i << "-" << e << "," << a << "-" << d << " ; " <<
						c << "-" << f << "," << b << "-" << j << "\n";
#endif
					//QUESTION QUESTION: should this be j1 instead of j?
					// Reset the current ctbroken basepr array for storing pairing information of the latest 
					// broken pseudoknot stem
					/*for (j = 1; j <= ctbroken->numofbases; j++){
						ctbroken->basepr[count][j] = 0;
					}*/
					for (int j1 = 1; j1 <= ctbroken->GetSequenceLength(); j1++){
						//ctbroken->basepr[count][j1] = 0;
						ctbroken->RemovePair(j1,count);
					}
					// Determine which pseudoknot formation penalty to use
					int beta1 = 0;
					bool pseudoOrMulti = false;
					for (iter = 1; iter <= i; iter++) {
						if (ct->GetPair(iter,count) > j) pseudoOrMulti = true;
					}
					if (pseudoOrMulti) beta1 = BETA_1PM;
					else beta1 = BETA_1;

					// Calculate RNAstructure hairpin energy of each pseudoknot stem
					
					ctbroken->SetPair(a,d,count);
					
					efn2(data, ctbroken, count, false);
					int hairpin1Energy = ctbroken->GetEnergy(count);
					int stackingEnergy1 = data->tstack[ctbroken->numseq[d]][ctbroken->numseq[a]]
						[ctbroken->numseq[d+1]][ctbroken->numseq[a-1]];
					
					ctbroken->RemovePair(a,count);

					
					ctbroken->SetPair(c,f,count);
					efn2(data, ctbroken, count, false);
					int hairpin2Energy = ctbroken->GetEnergy(count);					
					int stackingEnergy2 = data->tstack[ctbroken->numseq[f]][ctbroken->numseq[c]]
						[ctbroken->numseq[f+1]][ctbroken->numseq[c-1]];
					
					ctbroken->RemovePair(c, count);

#if defined(VERBOSE_MODE)
					cout<<"Hairpin Energies: "<<hairpin1Energy<<" "<<hairpin2Energy<<"\n";
#endif

					// Calculate paired and unpaired nucleotide penalties for the interior of the pseudoknot
					int beta2 = 0, beta3 = 0;

					// Add 2 beta2 penalties for a-d and c-f pairs
					beta2 += 2*(BETA_2);

					// Walk through gap region 1
					iter = a+1;
					while (iter < b){
						// Check for closing bp of a contained secondary structure in the gap
						if (ct->GetPair(iter,count) > iter && ct->GetPair(iter,count) < b){
							beta2 += BETA_2;
							iter = ct->GetPair(iter,count) + 1;
						}
						else if (ct->GetPair(iter,count) != 0 && 
						         (ct->GetPair(iter,count) <= a || ct->GetPair(iter,count) >= b)) {
							// If nt(iter) pairs outside the gap, add large penalty to the structure
							beta2 += 1000;
							iter++;
						}
						else {
							// Handles cases where nt(iter) is unpaired
							beta3 += BETA_3;
							iter++;
						}
					}

					// Walk through gap region 2
					iter = c+1;
					while (iter < d){
						// Check for closing bp of a contained secondary structure in the gap
						if (ct->GetPair(iter,count) > iter && ct->GetPair(iter,count) < d){
							beta2 += BETA_2;
							iter = ct->GetPair(iter,count) + 1;
						}
						else if (ct->GetPair(iter,count) != 0 && 
						         (ct->GetPair(iter,count) <= c || ct->GetPair(iter,count) >= d)) {
							// If nt(iter) pairs outside the gap, add large penalty to the structure
							beta2 += 1000;
							iter++;
						}
						else {
							// Handles cases where nt(iter) is unpaired
							beta3 += BETA_3;
							iter++;
						}
					}

					// Walk through gap region 3
					iter = e+1;
					while (iter < f){
						// Check for closing bp of a contained secondary structure in the gap
						if (ct->GetPair(iter,count) > iter && ct->GetPair(iter,count) < f){
							beta2 += BETA_2;
							iter = ct->GetPair(iter,count) + 1;
						}
						else if (ct->GetPair(iter,count) != 0 && 
						         (ct->GetPair(iter,count) <= e || ct->GetPair(iter,count) >= f)) {
							// If nt(iter) pairs outside the gap, add large penalty to the structure
							beta2 += 1000;
							iter++;
						}
						else {
							// Handles cases where nt(iter) is unpaired or paired to something outside the gap
							beta3 += BETA_3;
							iter++;
						}
					}

					// Scans both stems to determine which one is simpler to break
					for (iter = i; iter <= a; iter++) {
						if (ct->GetPair(iter,count) > iter) stem1bps++;
					}
					for (iter = d; iter <= e; iter++) {
						if (ct->GetPair(iter,count) > iter) stem1bps++;
					}
					for (iter = b; iter <= c; iter++) {
						if (ct->GetPair(iter,count) > iter) stem2bps++;
					}
					for (iter = f; iter < j; iter++) {//j now is iqual to numofbases+1, so must use 'iter<j' instead of 'iter<=j'
						if (ct->GetPair(iter,count) > iter) stem2bps++;
					}

					if (stem1bps <= stem2bps) {
#if defined(VERBOSE_MODE)
						cout << "Breaking stem 1...";
#endif

						for (iter = i; iter <= a; iter++) {
							if (ct->GetPair(iter,count) >= d && ct->GetPair(iter,count)<= e) {
#if defined(VERBOSE_MODE)
								cout << iter << "-" << ct->basepr[count][iter] << " ";
#endif

								ctbroken->SetPair(iter, ct->GetPair(iter,count),count);
								//ctbroken->basepr[count][ct->basepr[count][iter]] = 
								//	ct->basepr[count][ct->basepr[count][iter]];
								ct->RemovePair(iter,count);
								//ct->basepr[count][ct->basepr[count][iter]] = 0;
								//ct->basepr[count][iter] = 0;
							}
						}
#if defined(VERBOSE_MODE)
						cout << "\n";
#endif
					}

					else {
#if defined(VERBOSE_MODE)
						cout << "Breaking stem 2...\n";
#endif

						for (iter = b; iter <= c; iter++) {
							if (ct->GetPair(iter,count) >= f && ct->GetPair(iter,count) <= j) {
#if defined(VERBOSE_MODE)
								cout << iter << "-" << ct->basepr[count][iter] << " ";
#endif

								ctbroken->SetPair(iter, ct->GetPair(iter,count), count);
								//ctbroken->basepr[count][ct->basepr[count][iter]] = ct->basepr[count][ct->basepr[count][iter]];
								ct->RemovePair(iter,count);
								//ct->basepr[count][ct->basepr[count][iter]] = 0;
								//ct->basepr[count][iter] = 0;
							}
						}
#if defined(VERBOSE_MODE)
						cout << "\n";
#endif
					}
					
					// **Calculate energy of broken stem
					efn2(data, ctbroken, count, false);

#if defined(VERBOSE_MODE)
					cout << "Energy from broken stem = " << ctbroken->energy[count] << "\n";
#endif

					int uncorrectedEnergy = beta1+beta2+beta3+ctbroken->GetEnergy(count);

					// Update zero-order ctbroken basepr array with the latest broken pseudoknot stem
					for (j = 1; j <=ctbroken->GetSequenceLength(); j++){
						if (ctbroken->GetPair(j,count)!=0) brokenpairs[j] = ctbroken->GetPair(j,count);
					}

					return (uncorrectedEnergy-(hairpin1Energy-stackingEnergy1)-(hairpin2Energy-stackingEnergy2));
				}
			}				
		}
	}
	// no more pseudoknots found
	*(nopseudo) = true;
	return 0;
}

int ctout2 (structure *ct,int cutoff,const char *ctoutfile) {
	int count,i,j=0;
	char line[2*ctheaderlength],number[2*numlen];

	FILE *ctfile;
	ctfile=fopen(ctoutfile,"w");

	for (count=1;count<=(ct->GetNumberofStructures());count++) {

		if (cutoff !=0 && ct->GetEnergy(count) > cutoff) break;
		j++;

		strcpy(line,"");
		sprintf(line,"%5i",ct->GetSequenceLength());

		
		if (ct->GetEnergy(count)!=0) {
			strcat(line,"  ENERGY = ");

			//gcvt((float (ct->energy[count]))/conversionfactor,6,number);
			if (conversionfactor==10)
				sprintf(number,"%.1f",(float (ct->GetEnergy(count)))/conversionfactor);
			else if (conversionfactor==100)
				sprintf(number,"%.2f",(float (ct->GetEnergy(count)))/conversionfactor);
			else sprintf(number,"%f",(float (ct->GetEnergy(count)))/conversionfactor);
	
			strcat(line,number);
			strcat(line,"  ");
		}
		else strcat(line,"  ");
		strcat(line,ct->GetCtLabel(count).c_str());
		fputs (line,ctfile);
		for (i=1;i<ct->GetSequenceLength();i++) {
			if (ct->stacking) {
				sprintf(line,"%5i%2c%8i%5i%5i%5i%5i\n",
				        i,ct->nucs[i],(i-1),(i+1),ct->GetPair(i,count),ct->hnumber[i],ct->GetPair(i+ct->GetSequenceLength(),count));
			}
			else {
				sprintf(line,"%5i%2c%8i%5i%5i%5i\n",
				        i,ct->nucs[i],(i-1),(i+1),ct->GetPair(i,count),ct->hnumber[i]);
			}
			fputs(line,ctfile);
		}
		
   
		//last nucleotide not connected--
		i = ct->GetSequenceLength();
		if (ct->stacking) {
			sprintf(line,"%5i %1c%8i%5i%5i%5i%5i\n",
			        i,ct->nucs[i],(i-1),0,ct->GetPair(i,count),ct->hnumber[i],ct->GetPair(i+ct->GetSequenceLength(),count));
		}
		else {
			sprintf(line,"%5i %1c%8i%5i%5i%5i\n",
			        i,ct->nucs[i],(i-1),0,ct->GetPair(i,count),ct->hnumber[i]);
		}
		fputs(line,ctfile);

	}

	fclose (ctfile);
	return(j);
}

void pairout(structure *ct, int cutoff, const char* pairsoutfile) {
	int iter, iter2, a, b, c, d, bpCount;
	ofstream outFile(pairsoutfile);
	if (outFile.bad()) {
		cerr << "Error opening file " << pairsoutfile << "\n";
		return; 
	}
	for (iter = 1; iter <= ct->GetNumberofStructures(); iter++) {
		if (cutoff !=0 && ct->GetEnergy(iter) > cutoff) break;
		outFile << "#  Structure #" << iter << "; ENERGY = " << double(ct->GetEnergy(iter))/10 << "; " << 
			ct->GetCtLabel(iter) << endl;
		outFile << "1" << endl;
		for (iter2 = 1; iter2 <= ct->GetSequenceLength(); iter2++) outFile << ct->nucs[iter2];
		outFile << endl << endl;
		outFile << "Positions paired." << endl;
		iter2 = 1;
		while (iter2 <= ct->GetSequenceLength()) {
			if (ct->GetPair(iter2,iter) < iter2) {
				iter2++;
				continue;
			}
			else {
				a = iter2;
				b = ct->GetPair(iter2,iter);
				bpCount = 1;
				while (ct->GetPair(iter2+bpCount,iter) == (b-bpCount)) bpCount++;
				c = iter2 + bpCount - 1;
				d = ct->GetPair(iter2+(bpCount-1),iter);
				outFile << a << "-" << c << "; " << d << "-" << b << endl;
				iter2 = c + 1;
				continue;
			}
		}
		outFile << endl << "---------------------------------------" << endl << endl;
	}
	if (outFile.bad()) {
		cerr << "Error writing to file " << pairsoutfile << "\n";
		return; 
	}
	outFile.close();
	return;
}

bool TripletMatch(structure *ct, char**** tripletArray, short** xCoords, short** yCoords, int currConNum, int currConPos) {
	
	bool TripConCheck=true;
		
	//Check if the nucleotides in the constraint are not separated by bulges are loops 
	if(currConPos>2&&xCoords[currConNum][currConPos-2]-xCoords[currConNum][currConPos]==2&&yCoords[currConNum][currConPos]-yCoords[currConNum][currConPos-2]==2){
		//if there is a triplet constraint in the PREVIOUS position, this means that there is an active constraint
		if(tripletArray[currConNum][currConPos-1][0][0]!=0){
			//if the constraint is on XCoord side
			if(NucCompare(tripletArray[currConNum][currConPos-1][0][2],ct->numseq[xCoords[currConNum][currConPos-1]],0)){
#if defined(TRIPLET_VERBOSE_MODE)
				cerr << "    Triplet Constraint found on xCoord:\n    ";
#endif
				//if the nucleotides matchs the constraint
				if((NucCompare(tripletArray[currConNum][currConPos-1][0][1],ct->numseq[xCoords[currConNum][currConPos]],tripletArray[currConNum][currConPos-1][0][0])) &&
				   (NucCompare(tripletArray[currConNum][currConPos-1][0][3],ct->numseq[xCoords[currConNum][currConPos-2]],tripletArray[currConNum][currConPos-1][0][0]))){
#if defined(TRIPLET_VERBOSE_MODE)
					cerr << "\n       Match at (x+1)[" << xCoords[currConNum][currConPos] << "]=" << ct->numseq[xCoords[currConNum][currConPos]] << " and (x-1)[" << xCoords[currConNum][currConPos-2] << "]=" << ct->numseq[xCoords[currConNum][currConPos-2]] << endl;
#endif
					//QUESTION: temporary line
					cerr << " foo " << tripletArray[currConNum][currConPos-1][0][0] << " " << xCoords[currConNum][currConPos] << "-" << xCoords[currConNum][currConPos-2] << endl;
					
					//if there is a second constraint for that match
					if(tripletArray[currConNum][currConPos-1][1][0]!=0){
						//if the second constraint doesnt' match
						if(!(NucCompare(tripletArray[currConNum][currConPos-1][1][1],ct->numseq[xCoords[currConNum][currConPos]],tripletArray[currConNum][currConPos-1][1][0])) &&
						   !(NucCompare(tripletArray[currConNum][currConPos-1][1][3],ct->numseq[xCoords[currConNum][currConPos-2]],tripletArray[currConNum][currConPos-1][1][0]))){
#if defined(TRIPLET_VERBOSE_MODE)
							cerr << "\n             !!! Second Constraint Did Not Match !!!";
#endif
							TripConCheck=false;//set TripConCheck to FALSE if the second constraint didn't match
						}//END: if the second constraint doesn't match
					}//END: if there is a second constraint
				}
				else{
#if defined(TRIPLET_VERBOSE_MODE)
					cerr << "\n                     !!! First Constraint Did Not Match !!!";
#endif
					TripConCheck=false;//set TripConCheck to false if the first constraint didn't match
				}//END: if the nucleotide in the third position matches the constraint
			}//END: if the constraint is on xCoord.

			//if the constraint is on yCoord side
			else if(NucCompare(tripletArray[currConNum][currConPos-1][0][2],ct->numseq[yCoords[currConNum][currConPos-1]],0)){
#if defined(TRIPLET_VERBOSE_MODE)
				cerr << "    Triplet Constraint found on yCoord:\n    ";
#endif
				//if the nucleotide matches the constraint
				if((NucCompare(tripletArray[currConNum][currConPos-1][0][3],ct->numseq[yCoords[currConNum][currConPos]],tripletArray[currConNum][currConPos-1][0][0])) &&
				   (NucCompare(tripletArray[currConNum][currConPos-1][0][1],ct->numseq[yCoords[currConNum][currConPos-2]],tripletArray[currConNum][currConPos-1][0][0]))){
#if defined(TRIPLET_VERBOSE_MODE)
					cerr << "\n      Match at (y+1)[" << yCoords[currConNum][currConPos] << "]=" << ct->numseq[yCoords[currConNum][currConPos]] << " and (y-1)[" << yCoords[currConNum][currConPos-2] << "]=" << ct->numseq[yCoords[currConNum][currConPos-2]] << endl;
#endif
					//QUESTION: temporary lines
					cerr << " foo " << tripletArray[currConNum][currConPos-1][0][0] << " " << yCoords[currConNum][currConPos] << "-" << yCoords[currConNum][currConPos-2] << endl;

					//if there is a second constraint for that match
					if(tripletArray[currConNum][currConPos-1][1][0]!=0){
						//if the second constraint doesnt' match
						if(!(NucCompare(tripletArray[currConNum][currConPos-1][1][3],ct->numseq[yCoords[currConNum][currConPos]],tripletArray[currConNum][currConPos-1][1][0])) &&
						   !(NucCompare(tripletArray[currConNum][currConPos-1][1][1],ct->numseq[yCoords[currConNum][currConPos-2]],tripletArray[currConNum][currConPos-1][1][0]))){
#if defined(TRIPLET_VERBOSE_MODE)
							cerr << "\n             !!! Second Constraint Did Not Match !!!";
#endif
							TripConCheck=false;//set TripConCheck to FALSE if the second constraint didn't match
						}//END: if the second constraint doesn't match
					}//END: if there is a second constraint
				}
				else{
#if defined(TRIPLET_VERBOSE_MODE)
					cerr << "\n                     !!! First Constraint Did Not Match !!!";
#endif
					TripConCheck=false;//set TripConCheck to false if the first constraint didn't match
				}//END: if the nucleotide in the third position matches the constraint
			}//END: if the constraint is on xCoord.

#if defined(TRIPLET_VERBOSE_MODE)
			if(TripConCheck==true)cerr << "\n**********MATCH*********\n";
			else cerr << "\n___________NO MATCH__________\n";
#endif
			
		}//END: if there is a triplet constraint for the previous position
	}//END: if the nucleotides are not separated by loops or bulges

	return TripConCheck;
	/*Return TRUE if either: 1. There were no triplet constraint for either the current or previous positions
	  2. There was a bulge or a loop in the middle of the triplet constraint
	  3. There were constraints and they all matched
	
	  Return FALSE if:       1. The constraint(s) didn't match
	*/
}

bool NucCompare(char ConNuc,int SeqNuc, char ConSign){//Given a constraint and a nucleotide, and a constraint sigh, the program returns 'true' if there is a match and 'false' if not
	
	/* Structure class nomenclature:
	     A = 1 (R)
	     C = 2 (Y)
	     G = 3 (R)
	     U = 4 (Y)
	   NAPSS Triplet constraint nomenclature
	     Y = Pyrimidines
	     R = Purines
	*/
#if defined(TRIPLET_VERBOSE_MODE)
	cerr << "  Comparing Nucleotides: Constraint=" << ConNuc << " Sequence=" << SeqNuc << " Constraint_Sign=" << ConSign ;
#endif

	//If the input for the constraint sign was 0, this lets the program know that it is comparing nucleotides, and not constraints
	//This operation is performed when the program tries to determine where is (xCoords or yCoords) 'A' or 'G' nucleotide of the triplet constraint located
	if(ConSign==0){
		if(ConNuc=='A' && SeqNuc==1 || ConNuc=='G' && SeqNuc==3) return true;
		else return false;
	}
	//If the comparison sigh is a '+', this means that it is a positive constraint (Can be ONLY this, and nothing else).
	else if(ConSign=='+'){
		if(ConNuc=='Y' && (SeqNuc==2 || SeqNuc==4))	return true;//Check if the nucleotide is a pyrimidine
		else if(ConNuc=='R' && (SeqNuc==1 || SeqNuc==3)) return true;//Check if the nucleotide is a purine
		else if(ConNuc=='X') return true;//Check if there is an 'X' in the triplet constraint. This means that there is no data.
		else return false;
	}
	//If the comparison sigh is a '-', this means that it is a negotive constraint (Can be anything BUT this)
	else if(ConSign=='-'){
		if(ConNuc=='Y' && (SeqNuc!=2 || SeqNuc!=4)) return true;//Check if the nucleotide is NOT a pyrimidine
		else if(ConNuc=='R' && (SeqNuc!=1 || SeqNuc!=3)) return true;//Check if the nucleotide is NOT a purine
		else if(ConNuc=='X') return true;//Check if there is an 'X' in the triplet constraint. This means that there is no data.
		else return false;
	}


}

// Check for helical extension possibilities
void HelicalExtension(structure* ct, short** convertedDotplot, int* helixExtend, short** dgArray){
	for (int i = 1; i <= ct->GetNumberofStructures(); i++) {
		bool extensionAdded = true;
		while (extensionAdded) {
			extensionAdded = false;
			// Initialize the temporary array
			for (int iter2 = 0; iter2 <= ct->GetSequenceLength(); iter2++) {
				helixExtend[iter2] = 0;
			}
			
			int iter2 = 1;
			while (iter2 <= ct->GetSequenceLength()) {
				if (ct->GetPair(iter2,i) < iter2) {
					iter2++;
					continue;
				}
				else {
					// Check to see if the nucleotides downstream can pair with each other
					if (iter2 > 1 && ct->GetPair(iter2,i) < ct->GetSequenceLength() &&
					    ct->GetPair(iter2-1,i) == 0 && ct->GetPair(ct->GetPair(iter2,i)+1,i) == 0) { 
						// Nucleotides are currently unpaired 
						if (convertedDotplot[iter2-1][ct->GetPair(iter2,i)+1] !=0) { 
							// They could form a valid, stable pair
							if (helixExtend[iter2-1] == 0 && helixExtend[ct->GetPair(iter2,i)+1] == 0) {
								// There are no other possible extensions to these nucleotides yet - store this one
								helixExtend[iter2-1] = ct->GetPair(iter2,i)+1;
								helixExtend[ct->GetPair(iter2,i)+1] = iter2-1;
							}
							else {
								// One or both nucleotides already have a possible extension - if both do, 
								// make no changes because one new pair shouldn't break up two other possible pairs
								if (helixExtend[iter2-1] != 0 && helixExtend[ct->GetPair(iter2,i)+1] == 0) {
									// iter2 already has a possible helical extension to it - determine if the latest
									// pairing is more favorable in the original dotplot
									if (dgArray[iter2-1][helixExtend[iter2-1]] < 
									    dgArray[iter2-1][ct->GetPair(iter2,i)+1]) {
										helixExtend[helixExtend[iter2-1]] = 0;
										helixExtend[iter2-1] = ct->GetPair(iter2,i)+1;
										helixExtend[ct->GetPair(iter2,i)+1] = iter2-1;
									}
								}
								if (helixExtend[ct->GetPair(iter2,i)+1] != 0 && helixExtend[iter2-1] == 0) {
									// the nucleotide opposite iter2 already has a possible helical extension to it - 
									// determine if the latest pairing is more favorable in the original dotplot
									if (dgArray[ct->GetPair(iter2,i)+1][helixExtend[ct->GetPair(iter2,i)+1]] < 
									    dgArray[ct->GetPair(iter2,i)+1][iter2-1]) {
										helixExtend[helixExtend[ct->GetPair(iter2,i)+1]] = 0;
										helixExtend[ct->GetPair(iter2,i)+1] = iter2-1;
										helixExtend[iter2-1] = ct->GetPair(iter2,i)+1;
									}
								}
							}
						}
					}

					// Check to see if the nucleotides upstream can pair with each other
					if (ct->GetPair(iter2+1,i) == 0 && ct->GetPair(ct->GetPair(iter2,i)-1,i) == 0) { 
						// Nucleotides are currently unpaired 
						if (convertedDotplot[iter2+1][ct->GetPair(iter2,i)-1] !=0) { 
							// They could form a valid, stable pair
							if (helixExtend[iter2+1] == 0 && helixExtend[ct->GetPair(iter2,i)-1] == 0) {
								// There are no other possible extensions to these nucleotides yet - store this one
								helixExtend[iter2+1] = ct->GetPair(iter2,i)-1;
								helixExtend[ct->GetPair(iter2,i)-1] = iter2+1;
							}
							else {
								// One or both nucleotides already have a possible extension - if both do, 
								// make no changes because one new pair shouldn't break up two other possible pairs
								if (helixExtend[iter2+1] != 0 && helixExtend[ct->GetPair(iter2,i)-1] == 0) {
									// iter2 already has a possible helical extension to it - determine if the latest
									// pairing is more favorable in the original dotplot
									if (dgArray[iter2+1][helixExtend[iter2+1]] < 
									    dgArray[iter2+1][ct->GetPair(iter2,i)-1]) {
										helixExtend[helixExtend[iter2+1]] = 0;
										helixExtend[iter2+1] = ct->GetPair(iter2,i)-1;
										helixExtend[ct->GetPair(iter2,i)-1] = iter2+1;
									}
								}
								if (helixExtend[ct->GetPair(iter2,i)-1] != 0 && helixExtend[iter2+1] == 0) {
									// the nucleotide opposite iter2 already has a possible helical extension to it - 
									// determine if the latest pairing is more favorable in the original dotplot
									if (dgArray[ct->GetPair(iter2,i)-1][helixExtend[ct->GetPair(iter2,i)-1]] < 
									    dgArray[ct->GetPair(iter2,i)-1][iter2+1]) {
										helixExtend[helixExtend[ct->GetPair(iter2,i)-1]] = 0;
										helixExtend[ct->GetPair(iter2,i)-1] = iter2+1;
										helixExtend[iter2+1] = ct->GetPair(iter2,i)-1;
									}
								}
							}
						}
					}
					iter2++;
					continue;
				}
			}

			// Now loop through the temporary array and add any remaining potential pairs to the final structure
			for (int iter2 = 1; iter2 <= ct->GetSequenceLength(); iter2++) {
				if (helixExtend[iter2] != 0) {
					// Adding a pair - flip the switch
					extensionAdded = true;
					ct->SetPair(iter2,helixExtend[iter2],i);
				}
			}
		}
	}
}

bool Parse(int argc, char* argv[]){
	
	// Create the command line parser and build in its required parameters.
	ParseCommandLine* parser = new ParseCommandLine( "NAPSS" );
	parser->addParameterDescription( "seq file", "The name of a file containing an input sequence." );
	parser->addParameterDescription( "NMR file", "The name of an NMR file with constraints." );
	parser->addParameterDescription( "ct file", "The name of a CT file to which output will be written." );

	// Add the constraint file option.
	vector<string> constraintOptions;
	constraintOptions.push_back( "-c" );
	constraintOptions.push_back( "-C" );
	constraintOptions.push_back( "--constraint" );
	parser->addOptionFlagsWithParameters( constraintOptions, "Specify a constraints file to be applied. Default is to have no constraints applied." );

	// Add the maximum percent energy difference to consider in the dotplot option.
	vector<string> DotPlotPercentOptions;
	DotPlotPercentOptions.push_back( "-d" );
	DotPlotPercentOptions.push_back( "-D" );
	DotPlotPercentOptions.push_back( "--DotPercent" );
	parser->addOptionFlagsWithParameters( DotPlotPercentOptions, "Specify a maximum percent energy difference to consider in the dotplot. Default is 5 percent." );
	
	// Add the maximum number of structures option.                                                                                                                                                        
	vector<string> maxStructuresOptions;
	maxStructuresOptions.push_back( "-m" );
	maxStructuresOptions.push_back( "-M" );
	maxStructuresOptions.push_back( "--maximum" );
	parser->addOptionFlagsWithParameters( maxStructuresOptions, "Specify a maximum number of structures per matched constraint set. Default is 100 structures." );
	
	// Add the percent energy difference option.
	vector<string> percentOptions;
	percentOptions.push_back( "-p" );
	percentOptions.push_back( "-P" );
	percentOptions.push_back( "--percent" );
	parser->addOptionFlagsWithParameters( percentOptions, "Specify a maximum percent energy difference. Default is 0 which means that all structures are outputted." );

	// Add the positions paired output file option.
	vector<string> posPairsOptions;
	posPairsOptions.push_back( "-pp" );
	posPairsOptions.push_back( "-PP" );
	posPairsOptions.push_back( "--posPaired" );
	parser->addOptionFlagsWithParameters( posPairsOptions, "Specify the name of the positions paired style output file. Default is to have no file specified." );

	// Add the pseudoknot-free output option.
	vector<string> pseudoknotOptions;
	pseudoknotOptions.push_back( "-pf" );
	pseudoknotOptions.push_back( "-PF" );
	pseudoknotOptions.push_back( "--pseudoknotFree" );
	parser->addOptionFlagsNoParameters( pseudoknotOptions, "Specify pseudoknot-free prediction mode. Default is to predict pseudoknots." );

	// Add the SHAPE option.
	vector<string> shapeOptions;
	shapeOptions.push_back( "-sh" );
	shapeOptions.push_back( "-SH" );
	shapeOptions.push_back( "--SHAPE" );
	parser->addOptionFlagsWithParameters( shapeOptions, "Specify a SHAPE date file to be used to generate restraints. These restraints specifically use SHAPE pseudoenergy restraints. Default is to have no SHAPE data file specified." );

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

	// Add the window size option.                                                                                                                                                                         
	vector<string> windowOptions;
	windowOptions.push_back( "-w" );
	windowOptions.push_back( "-W" );
	windowOptions.push_back( "--window" );
	parser->addOptionFlagsWithParameters( windowOptions, "Specify a window size. Default is 0 nucleotides." );

	// Parse the command line into pieces.
	parser->parseLine( argc, argv );

	// Get required parameters from the parser.                                                                                                                                                            
	if( !parser->isError() ) {
		inseq = parser->getParameter( 1 );
		inNMRconstraints = parser->getParameter( 2 );
		outct = parser->getParameter( 3 );
	}

	// Get the constraint file option.
	if( !parser->isError() ) { constraintFile = parser->getOptionString( constraintOptions, true ); }

	// Get the maximum number of structures option.
	if( !parser->isError() ) {
		parser->setOptionInteger( maxStructuresOptions, maxtracebacks );
		if( maxtracebacks <= 0 ) { parser->setError( "maximum number of structures" ); }
	}
	
	// Get the percent energy difference option.
	if( !parser->isError() ) {
		parser->setOptionInteger( percentOptions, cutoff );
		if( cutoff < 0 ) { parser->setError( "percent energy difference" ); }
	}
	
	// Get the window size option.
	if( !parser->isError() ) {
		parser->setOptionInteger( windowOptions, windowsize );
		if( windowsize < 0 ) { parser->setError( "window size" ); }
	}

	// Get the DotPLotPercent energy difference option.
	if( !parser->isError() ) {
		parser->setOptionInteger( DotPlotPercentOptions, percent );
		if( percent < 0 ) { parser->setError( "DotPlot Percent energy difference" ); }
	}

	// Get the pseudoknot-free output option.
	pseudoknotFree = parser->contains( pseudoknotOptions);

	// Get the SHAPE data and options.
	if( !parser->isError() ) {
		inSHAPEfile = parser->getOptionString( shapeOptions );
		if( !parser->isError() ) { parser->setOptionDouble( shapeInterceptOptions, intercept ); }
		if( !parser->isError() ) { parser->setOptionDouble( shapeSlopeOptions, slope ); }
	}

	// Get the outpairs file option.
	if( !parser->isError() ) { outpairs = parser->getOptionString( posPairsOptions, false ); }
	
	// Delete the parser and return whether the parser encountered an error.
	bool noError = ( parser->isError() == false );
	delete parser;
	return noError;
}
