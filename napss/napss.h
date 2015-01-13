/*
 * NAPSS, a program that predicts RNA secondary structures with pseudoknots with the aid of NMR constraints data.
 *
 * (c) 2012 Mathews Lab, University of Rochester Medical Center.
 * Written by Stanislav Bellaousov
 */

#ifndef NAPSS_H
#define NAPSS_H


#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <cstring>
#include <cstdlib>

//#include "stdafx.h"
#include "../RNA_class/RNA.h"
#include "../src/configfile.h"
#include "../src/defines.h"
#include "../src/platform.h"
#include "../src/rna_library.h"
#include "../src/structure.h"
#include "../src/algorithm.h"
#include "../src/ParseCommandLine.h"

using namespace std;

#define BETA_1 96
#define BETA_1PM 150
#define BETA_2 1
#define BETA_3 1

typedef vector<short> conMatch;

void getdat(char *loop, char *stackf, char *tstackh, char *tstacki,
            char *tloop, char *miscloop, char *danglef, char *int22,
            char *int21,char *coax, char *tstackcoax,
            char *coaxstack, char *tstack, char *tstackm, char *triloop,
            char *int11, char *hexaloop, char *tstacki23, char *tstacki1n,
            char *datapath, bool isRNA);

void firstDotMatch(short**,bool*,bool*,short,short,short**,short**,vector<conMatch>*,short*,short*,short**,short**,int,short, structure *ct,char**** tripletArray);
// Loops through x,y of entire dotplot searching for matches to the first basepair of each constraint

void recursiveMatch(short**,bool*,short,short,short**,short**,vector<conMatch>*,short*,short*,short**,short**,int,short,bool*,structure *ct,char**** tripletArray);
// Takes over for firstDotMatch to find the matches to the remaining basepairs in a constraint

void pseudodptrim(structure*,int*, int*); // Searches dotplot for complicated pseudoknot folds to exclude

void efn2mod(datatable *data,structure *ct, int structnum, bool simplemb, structure *ctbroken);
// Pseudoknot-capable energy calculation

int pseudoenergy(datatable *data, structure *ct, int count, bool *nopseudos, structure *ctbroken, int *init, int *brokenpairs);
// Searches ct for pseudoknots, returns energy for breaking apart the pseudoknot (including Beta penalties)

int ctout2(structure *ct, int cutoff, const char *ctoutfile); 
// Modified version of ctout that can truncate unstable structures; returns number of output structures

void pairout(structure *ct, int cutoff, const char* pairsoutfile); // Outputs PseudoViewer3 structural data

bool TripletMatch(structure *ct, char**** tripletArray, short** xCoords, short** yCoords, int currConNum, int currConPos);//Match triplet constraints

bool NucCompare(char ConNuc, int SeqNuc, char ConSign=0);//Given a constraint and a nucleotide the program returns 'true' if there is a match and 'false' if not

void HelicalExtension(structure *ct, short** convertedDotplot, int* helixExtended, short** dgArray);//Check for helical extension possibilities

bool Parse(int argc, char* argv[]);//Parses the data from the command line and returns true is parsing works properly, and false if it doesn't

string inseq;//address for the input sequence file
string inNMRconstraints;//address for the NMR constraints file
string outct;//address for the output ct file
string inSHAPEfile;//address for the input SHAPE file
string outpairs;//address for the output file with pairs
string constraintFile;//address for the structural constraints file
double slope=1.8,intercept=-0.6;//shape slope and shape intercept
int maxtracebacks=100,percent=5,windowsize=0,cutoff=0;
bool pseudoknotFree=false;//default is to have pseudoknot prediction mode enabled
int warningLimit=50000;//number of matches (matchVector->size) before the warning message is printed
bool ifwarningMessage=false;//keeps track if the warning message was printed. false - warning message was not printed yet, so should be printed when warning limit is reached; true - warning message has already been printed, so no need to print it again.

#endif /* NAPSS_H */
