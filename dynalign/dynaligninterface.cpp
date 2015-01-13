/*
 * This is the command-line interface to Dynalign, 
 * written by David Mathews; Copyright 2002, 2003, 2004, 2005, 2006
 *
 * Contributors:
 *  Chris Connett and Andrew Yohn, 2006
 *  Josh Keegan, 2006
 *  Arif Harmanci, 2006
 *
 * Dynalign is described in:
 * Mathews & Turner, JMB, 317:191-203 (2002).
 *
 *
 * 
 *----------------------------------------------------------------
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 * 02111-1307, USA.
 */

#include <fstream>
#include <iostream>
#include <cstring>
#include <cstdlib>

#include "../src/configfile.h"
#include "../src/defines.h"
#include "../src/dynalign.h"
#include "../src/dynalignheap.h"
#include "../src/outputconstraints.h"
#include "../src/platform.h"
#include "../src/rna_library.h"
#include "../src/structure.h"
#include "../src/TProgressDialog.h"
#include "../src/phmm/phmm_aln.h"
#include "../src/phmm/structure/structure_object.h"

//TOLERANCE is the maximum deviation from 310.15 K before which the enthalpy parameters 
//are read from disk to adjust the free energyy changes from 310.15 K.
#define TOLERANCE 0.01 

#ifdef COMPILE_SMP
#include "../src/observingtextprogressbar.h"
#endif

using namespace std;



/*	Function getdat

	Function gets the names of data files to open

*/

void GetDat(char *loop, char *stackf, char *tstackh, char *tstacki,
	char *tloop, char *miscloop, char *danglef, char *int22,
	char *int21,char *coax, char *tstackcoax,
    char *coaxstack, char *tstack, char *tstackm, char *triloop,
    char *int11, char *hexaloop, char *tstacki23, char *tstacki1n,
	char *datapath, bool isRNA, bool isEnthalpy=false)

{

 
  //copy the path to each name
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


  if( !isRNA) {
	  //these are dna parameters and so they need to start with "dna"
	strcat (loop,"dna");
	strcat (stackf,"dna");
	  strcat (tstackh,"dna");
	  strcat (tstacki,"dna");
	  strcat (tloop,"dna");
	  strcat (miscloop,"dna");
	  strcat (danglef,"dna");
	  strcat (int22,"dna");
	  strcat (int21,"dna");
	  strcat (triloop,"dna");
	  strcat (coax,"dna");
	  strcat (tstackcoax,"dna");
	  strcat (coaxstack,"dna");
	  strcat (tstack,"dna");
	  strcat (tstackm,"dna");
	  strcat (int11,"dna");
	  strcat (hexaloop,"dna");
	  strcat (tstacki23,"dna");
	  strcat (tstacki1n,"dna");  
	 
  }
	 
	//append the actual file name
  strcat (loop,"loop.");
  strcat (stackf,"stack.");
  strcat (tstackh,"tstackh.");
  strcat (tstacki,"tstacki.");
  strcat (tloop,"tloop.");
  strcat (miscloop,"miscloop.");
  strcat (danglef,"dangle.");
  strcat (int22,"int22.");
  strcat (int21,"int21.");
  strcat (triloop,"triloop.");
  strcat (coax,"coaxial.");
  strcat (tstackcoax,"tstackcoax.");
  strcat (coaxstack,"coaxstack.");
  strcat (tstack,"tstack.");
  strcat (tstackm,"tstackm.");
  strcat (int11,"int11.");
  strcat (hexaloop,"hexaloop.");
  strcat (tstacki23,"tstacki23.");
  strcat (tstacki1n,"tstacki1n.");
  
  if (isEnthalpy) {
	  //thse are enthalpy partameters so they need to end in .dh
	strcat (loop,"dh");
	strcat (stackf,"dh");
	  strcat (tstackh,"dh");
	  strcat (tstacki,"dh");
	  strcat (tloop,"dh");
	  strcat (miscloop,"dh");
	  strcat (danglef,"dh");
	  strcat (int22,"dh");
	  strcat (int21,"dh");
	  strcat (triloop,"dh");
	  strcat (coax,"dh");
	  strcat (tstackcoax,"dh");
	  strcat (coaxstack,"dh");
	  strcat (tstack,"dh");
	  strcat (tstackm,"dh");
	  strcat (int11,"dh");
	  strcat (hexaloop,"dh");
	  strcat (tstacki23,"dh");
	  strcat (tstacki1n,"dh");  
  }
  else {
	  //these are free energy parameters and the files end in .dat
	strcat (loop,"dat");
	strcat (stackf,"dat");
	  strcat (tstackh,"dat");
	  strcat (tstacki,"dat");
	  strcat (tloop,"dat");
	  strcat (miscloop,"dat");
	  strcat (danglef,"dat");
	  strcat (int22,"dat");
	  strcat (int21,"dat");
	  strcat (triloop,"dat");
	  strcat (coax,"dat");
	  strcat (tstackcoax,"dat");
	  strcat (coaxstack,"dat");
	  strcat (tstack,"dat");
	  strcat (tstackm,"dat");
	  strcat (int11,"dat");
	  strcat (hexaloop,"dat");
	  strcat (tstacki23,"dat");
	  strcat (tstacki1n,"dat");  


  }


}

// Zane: to get seq filename to pass to sanity check
#ifdef CHECK_ARRAY
string seq1;
string seq2;
#endif

//The main entry point for Dynalign using a text interface:
int main(int argc, char* argv[]) {
  string inseq1;
  string inseq2;
  string outct;
  string outct2;
  string aout;
  int imaxseparation;
#ifdef DYNALIGN_II
  double slope;
  double intercept;
  int max_elongation;
#else
#endif
  double fgap;
  int maxpairs=-1;

#ifdef DYNALIGN_II  
  int islope;
  int iintercept;
#else
#endif
  int igapincrease;
  int maxtrace;
  int percent;
  int singlefold_subopt_percent,maximumpairingdistance1,maximumpairingdistance2;
  int bpwin;
  int awin;
#ifdef DYNALIGN_II
#else
  bool insert;
#endif
  string savefile;
  string constraint1;
  string constraint2;
  string constrainta;
  string shape1;
  string shape2;
  bool dsv_templated, ct_templated;
  string dsvtemplatename;
	float maxdsvchange;
  int numProcessors = 1;
  bool **allowed_alignments,optimalonly;
  bool local;
  FILE *check;//a c FILE for checking that specified files exist
  
	structure ct1,ct2;
	int i,query,index;
	datatable data;
	char loop[maxfil],stackf[maxfil],tstackh[maxfil],tstacki[maxfil],
		tloop[maxfil],miscloop[maxfil],danglef[maxfil],int22[maxfil],
    int21[maxfil],coax[maxfil],tstackcoax[maxfil],
    coaxstack[maxfil],tstack[maxfil],tstackm[maxfil],triloop[maxfil],int11[maxfil],
	  datapath[maxfil],hexaloop[maxfil],tstacki1n[maxfil],tstacki23[maxfil],*pointer;
	//int dsvtemplatenumber;
	short **align;

	double shapeslope1,shapeintercept1,shapeslope2,shapeintercept2;

	double Temperature;
	bool DNA;


	TProgressDialog *progress = NULL;
#ifndef COMPILE_SMP
  progress = new TProgressDialog();
#else
  progress = new ObservingTextProgressBar();
#endif
  
	bool constrained;
	short **forcealign;

	
	
  
  // Check if we have a command-line parameter; if not, fall back to
  // interactive mode.
	if( argc == 2 ) {

    ConfigFile config(argv[1]);
    
    // Check for mandatory config file settings.
    bool valid_config =
      config.contains("inseq1") &&
      config.contains("inseq2") &&
      config.contains("outct") &&
      config.contains("outct2") &&
      config.contains("aout");
      //(config.contains("imaxseparation")||config.contains("align_threshold")) &&
      //config.contains("fgap") &&
      //config.contains("maxtrace") &&
      //config.contains("percent") &&
      //config.contains("singlefold_subopt_percent") &&
      //config.contains("bpwin") &&
      //config.contains("awin") &&
      //config.contains("insert") &&
      //config.contains("dsv_templated");

#ifdef COMPILE_SMP
    valid_config = valid_config && config.contains("num_processors");
#endif
    
    if (valid_config) {
      // Read all mandatory config file settings.
      inseq1         = config.getOption<string>("inseq1");
      inseq2         = config.getOption<string>("inseq2");
      outct          = config.getOption<string>("outct");
      outct2         = config.getOption<string>("outct2");
      aout           = config.getOption<string>("aout");
      
	  

#ifdef COMPILE_SMP
      numProcessors = config.getOption<int>("num_processors");
#endif
    }

    // Read optional config file settings
    if (valid_config) {

		//check to see if a ct file is specified for sequence 1, i.e. ct_templated = true.
	   if (config.contains("ct_templated")) ct_templated = config.getOption<bool>("ct_templated");
	   else ct_templated = false;
	  if (config.contains("imaxseparation")) imaxseparation = config.getOption<int>("imaxseparation");
	  else imaxseparation = -99;//set a flag to indicate the M is not being used

#ifdef DYNALIGN_II
          if (config.contains("slope")) slope           = config.getOption<double>("slope");
	  else slope = 0.1;
          if (config.contains("intercept")) intercept           = config.getOption<double>("intercept");
	  else intercept = 0.5;
          if (config.contains("max_elongation")) max_elongation           = config.getOption<int>("max_elongation");
	  else max_elongation = 5;
#else
#endif
          if (config.contains("fgap")) fgap           = config.getOption<double>("fgap");
	  else fgap = 0.4;
      if (config.contains("maxtrace")) maxtrace       = config.getOption<int>("maxtrace");
	  else maxtrace = 750;
      if (config.contains("percent")) percent        = config.getOption<int>("percent");
	  else percent =20;
      if (config.contains("singlefold_subopt_percent")) singlefold_subopt_percent =
                       config.getOption<int>("singlefold_subopt_percent");
	  else singlefold_subopt_percent=30;
      if (config.contains("bpwin")) bpwin          = config.getOption<int>("bpwin");
	  else bpwin = 2;
      if (config.contains("awin")) awin           = config.getOption<int>("awin");
	  else awin = 1;
#ifdef DYNALIGN_II
#else
      if (config.contains("insert")) insert         = config.getOption<bool>("insert");
	  else insert = true;
#endif
      if (config.contains("dsv_templated")) dsv_templated  = config.getOption<bool>("dsv_templated");
	  else dsv_templated = false;
      if (config.contains("savefile")) {
        savefile = config.getOption<string>("savefile");
      }
      if (config.contains("maxpairs")) maxpairs = config.getOption<int>("maxpairs");
		else maxpairs = -1;
      if (config.contains("constraint_1_file")) {
        constraint1 = config.getOption<string>("constraint_1_file");
      }

      if (config.contains("constraint_2_file")) {
        constraint2 = config.getOption<string>("constraint_2_file");
      }

      if (config.contains("shape_1_file")) {
        shape1 = config.getOption<string>("shape_1_file");
      }

	  if (config.contains("shape_2_file")) {
        shape2 = config.getOption<string>("shape_2_file");
      }

	  if (config.contains("constraint_align_file")) {
        constrainta = config.getOption<string>("constraint_align_file");
      }

	  if (config.contains("optimal_only")) {
		optimalonly = config.getOption<bool>("optimal_only");
	  }
	  //by default, setup to determine suboptimal structures
	  else optimalonly=false;
	  if (config.contains("local")) {
		local = config.getOption<bool>("local");
	  }
	  //by default, run global alignment
	  else local=false;

	  //Read SHAPE parameters or set defaults
	  if (config.contains("shapeslope1")) shapeslope1 = config.getOption<double>("shapeslope1");
	  else shapeslope1 = 1.8;
	  if (config.contains("shapeslope2")) shapeslope2 = config.getOption<double>("shapeslope2");
	  else shapeslope2 = 1.8;
	  if (config.contains("shapeintercept1")) shapeintercept1 = config.getOption<double>("shapeintercept1");
	  else shapeintercept1 = -0.6;
	  if (config.contains("shapeintercept2")) shapeintercept2 = config.getOption<double>("shapeintercept2");
	  else shapeintercept2 = -0.6;

	  if (config.contains("maximumpairingdistance1")) maximumpairingdistance1 = config.getOption<int>("maximumpairingdistance1");
	  else maximumpairingdistance1 = 0;
	  if (config.contains("maximumpairingdistance2")) maximumpairingdistance2 = config.getOption<int>("maximumpairingdistance2");
	  else maximumpairingdistance2 = 0;
	  if (config.contains("temperature")) Temperature = config.getOption<double>("temperature");
	  else Temperature = 310.15;
	  if (config.contains("DNA")) DNA = config.getOption<bool>("DNA");
	  else DNA = false;

    }

    // Read settings dependent on previous settings.
    if (valid_config && dsv_templated) {
      valid_config =
        config.contains("dsvtemplatename") &&
        config.contains("maxdsvchange");
      if (valid_config) {
        dsvtemplatename = config.getOption<string>("dsvtemplatename");
        maxdsvchange    = config.getOption<float>("maxdsvchange");
      }
    }
    
    if (!valid_config) {
      cerr << "ERROR: At least one parameter could not be read from the "
           << "configuration file. Aborting." << endl;
      exit(1);
    }
      
  } else {
	  
	  //set ct_templated to false by default:
	  ct_templated = false;
      
    // Interactive mode:
#ifdef DYNALIGN_II
	cout << "Usage: dynalign_ii [config file]\nNo config file specified, using interactive mode:\nNote that most, but not all, features are available through the interactive mode.\n\n";
#else
	cout << "Usage: dynalign [config file]\nNo config file specified, using interactive mode:\nNote that most, but not all, features are available through the interactive mode.\n\n";
#endif      
    cout << "Enter the name of the first sequence: ";
    cin >> inseq1;
	    
    cout << "Enter the name of the second sequence: ";
    cin >> inseq2;
	    
    cout << "Enter the name of the first output ct file: ";
    cin >> outct;
	    
    cout << "Enter the name of the second output ct file: ";
    cin >> outct2;
	    
    cout << "Enter name for the output of the alignment: ";
    cin >> aout;
	    
    cout << "Enter the max separation (M) <Recommend enter -99 for constraint by HMM forward-backward>: ";
    cin >> imaxseparation;
#ifdef DYNALIGN_II
    cout << "Enter the slope of domain insertion penalty: <Recommend 0.1>";
    cin >> slope;

    cout << "Enter the intercept of domain insertion penalty: <Recommend 0.5>";
    cin >> intercept;


    cout << "Enter the maximum value for helix elongation: <Recommend 5>";
    cin >> max_elongation;
#else
#endif
    cout << "Enter the gap penalty: <Recommend 0.4>";
    cin >> fgap;
	    
    cout << "Enter the maximum number of structures: ";
    cin >> maxtrace;
	    
    cout << "Enter the maximum percent energy difference "
         << "(where entering 20 is 20%): ";
    cin >> percent;

    cout << "Enter the maximum percent energy difference for keeping "
         << "suboptimal structures from the single-sequence folding "
         << "calculations: <Recommend 30>";
    cin >> singlefold_subopt_percent;
	    
    cout << "Enter the base pair window size: <Recommend 2>";
    cin >> bpwin;
	    
    cout << "Enter the alignment window size: <Recommend 1>";
    cin >> awin;
#ifdef DYNALIGN_II
#else    
    cout << "Allow single BP inserts into one sequence? (1/0) <Recommend 1>";
    cin >> insert;
#endif	    
    cout << "Write save file (needed for refold or dot plots)? (1/0) ";
    cin >> query;
	    
    if (query) {
      cout << "Enter the save file name: ";
      cin >> savefile;
    }
	    
    cout << "Input constraints for sequence 1? (1/0) ";
    cin >> query;
	    
    if (query) {
      cout << "Enter the constraint file name: ";
      cin >> constraint1;
    }
	    
    cout << "Input constraints for sequence 2? (1/0) ";
    cin >> query;
	    
    if (query) {
      cout << "Enter the constraint file name: ";
      cin >> constraint2;
    }
	    
    cout << "Input constraints for alignment? (1/0) ";
    cin >> query;
	    
    if (query) {
      cout << "Enter the constraint file name: ";
      cin >> constrainta;
    }
    
    cout << "Dsv template to be used (i.e. a progressive calculation) "
         << "(1/0): ";
    cin >> dsv_templated;
    if(dsv_templated)
      {
        cout << "Specify dsv template file: ";
        cin >> dsvtemplatename;
        cout << "Percent energy window for bases to keep from dsv template "
             << "(where 20 = 20% difference from optimal): ";
        cin >> maxdsvchange;
		cout << "Number of base pairs to allow for template <recommend -1, meaning the length of the 1st sequence: ";
		cin >> maxpairs;
      }
#ifdef COMPILE_SMP
    cout << "Enter the number of processors on this machine: ";
    cin >> numProcessors;
#endif
  }

  // Validate settings
#ifdef DYNALIGN_II
  islope = (int)(slope * 10.0);
  iintercept = (int)(intercept * 10.0);
#else
#endif

  igapincrease = (int)(fgap * 10.0);

#ifdef COMPILE_SMP
  if (numProcessors < 1) {
    cerr << "Warning: invalid number of processors specified; "
         << "defaulting to 1." << endl;
    numProcessors = 1;
  }
#endif

	//open the sequences
  if (ct_templated) {
	//This means the first sequence should be read from a ct file, not a .seq file.
	ct1.openct(inseq1.c_str());

	//also append a needed \n to the description string
	//try skipping this
	//strcat(ct1.ctlabel[1],"\n");


  }
	else if(!ct1.openseq(inseq1.c_str())) {cerr << "ERROR: Could not open sequence file "<<inseq1<<"\n"; return -1; }
	
	//open the sequence file for sequence 2.
	if(!ct2.openseq(inseq2.c_str())) {cerr << "ERROR: Could not open sequence file "<<inseq2<<"\n"; return -1; }
	
	
	
	
	//The next section was commented out by DHM when DNA and temperature were added
	//get the location of the data files
/*	pointer = getenv("DATAPATH");
	if (pointer!=NULL) {
		strcpy(datapath,pointer);
		strcat(datapath,"/");
	}
	else strcpy(datapath,"");


	//open the thermodynamic data tables
	getdat (loop, stackf, tstackh, tstacki,tloop, miscloop, danglef, int22,
          int21,coax, tstackcoax,coaxstack, tstack, tstackm, triloop,
          int11, hexaloop, tstacki23, tstacki1n, datapath, true);//the true indicates RNA parameters
	if (opendat (loop, stackf, tstackh, tstacki,tloop, miscloop, danglef, int22,
               int21,coax, tstackcoax,coaxstack, tstack, tstackm, triloop,
               int11,hexaloop,tstacki23, tstacki1n, &data)==0) {
    cerr << "A data file was lost\n";
	cerr << "Please set the environment variable DATAPATH to indicate the location of the data_tables/ directory.\n";
	cerr << "The current value of DATAPATH is: "<<datapath<<"\n";
    exit(1);
  }*/

	
	//Get the information from $DATAPATH, if available
	pointer = getenv("DATAPATH");
	if (pointer!=NULL) {
		strcpy(datapath,pointer);
		strcat(datapath,"/");
	}
	else strcpy(datapath,"");
	
	//open the data files -- must reside in pwd or $DATAPATH.
	//open the thermodynamic data tables
	GetDat(loop, stackf, tstackh, tstacki,tloop, miscloop, danglef, int22,
          int21,coax, tstackcoax,coaxstack, tstack, tstackm, triloop,
          int11, hexaloop, tstacki23, tstacki1n, datapath, (!DNA));//the false for DNA indicates RNA parameters
	if (opendat (loop,stackf,tstackh,tstacki,tloop,miscloop,danglef,int22,int21,
   		coax,tstackcoax,coaxstack,tstack,tstackm,triloop,int11,hexaloop,tstacki23, tstacki1n,&data)==0) {
      	
			cerr << "A data file was lost\n";
			cerr << "Please set the environment variable DATAPATH to indicate the location of the data_tables/ directory.\n";
			cerr << "The current value of DATAPATH is: "<<datapath<<"\n";
			return 0;//an error code
		
	}
	
		
	//now check to see if the temperature is other than 310.15 K:
	if (Temperature>(310.15+TOLERANCE)||Temperature<(310.15-TOLERANCE)) {

		//temperature is altered, so the enthalpy tables need to be read:
		datatable *localenthalpy;

		//get the names of the enthalpy files
		GetDat(loop, stackf, tstackh, tstacki,
			tloop, miscloop, danglef, int22,
			int21,coax, tstackcoax,
			coaxstack, tstack, tstackm, triloop,
			int11, hexaloop, tstacki23, tstacki1n, datapath, (!DNA),true);//rtue means this is enthalpy files

		//allocate a table to storte enthalpy parameters
		localenthalpy = new datatable();

		//open the enthlpy parameters and check for errors
		if (opendat(loop, stackf, tstackh, tstacki,
			tloop, miscloop, danglef, int22,
			int21,coax, tstackcoax,
			coaxstack, tstack, tstackm, triloop,
			int11,hexaloop,tstacki23, tstacki1n, localenthalpy)==0) {

			delete localenthalpy;

			cerr << "A data file was lost\n";
			cerr << "Please set the environment variable DATAPATH to indicate the location of the data_tables/ directory.\n";
			cerr << "The current value of DATAPATH is: "<<datapath<<"\n";
			return 0;//an error code

		}

		//using the enthalpy parameters and the folding free energy changes at 37 degrees C,
		//	calculate folding free energy changes at temp, the current temperature, and overwrite the 
		//	folding free energy changes at 37 degrees C
		dG_T((float)Temperature,data,*localenthalpy,data);

		//remove the enthalpy parameters, there is no need to keep them
		delete localenthalpy;

	}

#ifdef CHECK_ARRAY
        seq1 = inseq1;
        seq2 = inseq2;
#endif
	//allocate space for the alignment
	align = new short *[maxtrace];//maximum number of tracebacks and next line and below at delete
	for (i=0;i<maxtrace;i++)  align[i] = new short [ct1.GetSequenceLength()+1];

	//do the dynalign calculation

	constrained = false;
	if (constraint1 != "") {
		//check that the file exists.
		if ((check = fopen(constraint1.c_str(), "r"))== NULL) {
			//the file is not found
			cerr << "constraint_1_file, "<<constraint1<<" not found.\n";
			fclose(check);
			return 1;		
		}
		fclose(check);
		constrained = true;
		readconstraints(constraint1.c_str(),&ct1);
	}
	if (constraint2 != "") {
		//check that the file exists.
		if ((check = fopen(constraint2.c_str(), "r"))== NULL) {
			//the file is not found
			cerr << "constraint_2_file, "<<constraint2<<" not found.\n";
			fclose(check);
			return 1;		
		}
		fclose(check);
		constrained = true;
		readconstraints(constraint2.c_str(),&ct2);
	}
	
	if (shape1 != "") {

		//check that the file exists.
		if ((check = fopen(shape1.c_str(), "r"))== NULL) {
			//the file is not found
			cerr << "shape_1_file, "<<shape1<<" not found.\n";
			fclose(check);
			return 1;		
		}
		fclose(check);

		ct1.SHAPEslope = shapeslope1*conversionfactor;
		ct1.SHAPEintercept = shapeintercept1*conversionfactor;

		ct1.ReadSHAPE(shape1.c_str());
	}
	if (shape2 != "") {
		//check that the file exists.
		if ((check = fopen(shape2.c_str(), "r"))== NULL) {
			//the file is not found
			cerr << "shape_2_file, "<<shape2<<" not found.\n";
			fclose(check);
			return 1;		
		}
		fclose(check);

		ct2.SHAPEslope = shapeslope2*conversionfactor;
		ct2.SHAPEintercept = shapeintercept2*conversionfactor;

		ct2.ReadSHAPE(shape2.c_str());
	}
	if (constrainta != "") {
		//Apply alignment constraints

		//check that the file exists.
		if ((check = fopen(constrainta.c_str(), "r"))== NULL) {
			//the file is not found
			cerr << "constraint_align_file, "<<constrainta<<" not found.\n";
			fclose(check);
			return 1;		
		}
		fclose(check);
		constrained = true;
		forcealign=new short *[2];
		
		forcealign[0]=new short [ct1.GetSequenceLength()+1];
		forcealign[1]=new short [ct2.GetSequenceLength()+1];
		for (index=1;index<=ct1.GetSequenceLength();index++) {
			forcealign[0][index]=0;
		}
		for (index=1;index<=ct2.GetSequenceLength();index++) {
			forcealign[1][index]=0;
		}
		readalignmentconstraints(constrainta.c_str(),forcealign,&ct1,&ct2);
	} 
	else {
		forcealign = NULL;
	}
  
	//check if lowercase nucleotides were entered in either sequences - these will be single-stranded
	if (ct1.GetNumberofSingles()>0||ct2.GetNumberofSingles()>0) constrained=true;
	
	
	
	//This section folds thie individual sequences to find insignificant pairs to be ingnored by Dynalign.
	ct1.allocatetem();
	ct2.allocatetem();

	

	//dsv_templated is helpful for progressive calculations from multiple dynalign calculations
	if(dsv_templated) templatefromdsv( &ct1, dsvtemplatename.c_str(), maxdsvchange, maxpairs);

	//ct_templted is helpful if the structure is known for sequence 1
	else if (ct_templated) templatefromct(&ct1);

	//otherwise, fold the sequence to determine which pairs can result in low free energy structures
	else { 

		//if maximumpairingdistance1 is > 0 (the default), apply this to the templatefromfold calulation
		if (maximumpairingdistance1>0) {
			ct1.SetPairingDistance(maximumpairingdistance1);
		}
		templatefromfold(&ct1, &data, singlefold_subopt_percent);

	}

	//if maximumpairingdistance2 is > 0 (the default), apply this to the templatefromfold calulation
	if (maximumpairingdistance2>0) {
		
		ct2.SetPairingDistance(maximumpairingdistance2);
	}
	templatefromfold( &ct2, &data, singlefold_subopt_percent );

	//This next section determined the allowed nucleotide alignments if the HMM forward-backward is used:
	if (imaxseparation == -99) {
		//allocate space in allowed_alignments
		allowed_alignments = new bool *[ct1.GetSequenceLength()+1];
		for (i=0;i<=ct1.GetSequenceLength();i++) {
			allowed_alignments[i] = new bool [ct2.GetSequenceLength()+1];	
	
		}

		// Needed for having nucleotide sequences as c strings.
		ct1.nucs[ct1.GetSequenceLength() + 1] = 0;
		ct2.nucs[ct2.GetSequenceLength() + 1] = 0;


		//calculate_aln_probs_env(&ct1, &ct2, NULL, allowed_alignments, align_threshold);
		calculate_coinc_probs_env(&ct1, &ct2, allowed_alignments, forcealign);
		
	}
	else allowed_alignments = NULL;
#ifdef DYNALIGN_II
if (dynalign(&ct1, &ct2, align, imaxseparation, islope, iintercept, igapincrease, &data,
                 maxtrace, bpwin, awin, percent, forcealign, max_elongation, allowed_alignments, progress,
           (savefile != "") ? savefile.c_str() : NULL, optimalonly, local,
		   /*force =*/ constrained, numProcessors)==14)
#else
    if (dynalign(&ct1, &ct2, align, imaxseparation, igapincrease, &data,
           insert, maxtrace, bpwin, awin, percent, forcealign, allowed_alignments, progress,
           (savefile != "") ? savefile.c_str() : NULL, optimalonly, local,
		   /*force =*/ constrained, numProcessors)==14)
#endif

 {
		
		   //dynalign returns an int that indicates an error.  14 is currently the only possible error, a traceback error
		   cerr << "Dynalign encountered a traceback error.  Please report this possible bug to David_Mathews@urmc.rochester.edu\n";



	}
  
	//output the structures
	ct1.ctout(outct.c_str());
	ct2.ctout(outct2.c_str());

	//output the alignment
	alignout(align,aout.c_str(),&ct1,&ct2);



	//clean up memory allocation
	for (i=0;i<maxtrace;i++)  delete[] align[i];
	delete[] align;

	if (constrainta != "") {
		delete[] forcealign[0];
		delete[] forcealign[1];
		delete[] forcealign;

	}

  delete progress;


	if (imaxseparation == -99) {
		//delete space in allowed_alignments
		
		for (i=0;i<=ct1.GetSequenceLength();i++) {
			delete[] allowed_alignments[i];	
	
		}

		delete[] allowed_alignments;

	}
  
	return 0;
}


