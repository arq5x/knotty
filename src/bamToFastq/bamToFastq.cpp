/*
  ***************************************************************************
   bamToFastqMain.cpp (c) 2009 Aaron Quinlan

   Hall Lab
   Department of Biochemistry and Molecular Genetics
   University of Virginia

   All rights reserved.
 ***************************************************************************
*/

#define PROGRAM_NAME "bamToFastq"
// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
using namespace std;

// BamTools includes
#include "version.h"
#include "sequenceUtils.h"
#include "api/BamReader.h"
using namespace BamTools;



// function declarations
void ShowHelp(void);
	

int main(int argc, char* argv[]) {

	// our configuration variables
	bool showHelp = false;

	bool haveInBam     = true;
	bool haveFastq1    = false;
	bool haveFastq2    = false;
		
	// input files
	string inBamFile = "stdin";
	
	//output files
	string fastq1, fastq2;

	// check to see if we should print out some help
	if(argc <= 1) showHelp = true;

	for(int i = 1; i < argc; i++) {
		int parameterLength = (int)strlen(argv[i]);

		if((PARAMETER_CHECK("-h", 2, parameterLength)) || 
		(PARAMETER_CHECK("--help", 5, parameterLength))) {
			showHelp = true;
		}
	}

	if(showHelp) ShowHelp();


	// do some parsing (all of these parameters require 2 strings)
	for(int i = 1; i < argc; i++) {

		int parameterLength = (int)strlen(argv[i]);

 		if (PARAMETER_CHECK("-i", 2, parameterLength)) {
			if ((i+1) < argc) {
				inBamFile = argv[i + 1];
				i++;
			}
		}
		else if (PARAMETER_CHECK("-fq1", 4, parameterLength)) {
			if ((i+1) < argc) {
				haveFastq1 = true;
				fastq1 = argv[i + 1];
				i++;
			}
		}
		else if (PARAMETER_CHECK("-fq2", 4, parameterLength)) {
			if ((i+1) < argc) {
				haveFastq2 = true;
				fastq2 = argv[i + 1];
				i++;
			}
		}					
		else {
		  cerr << endl << "*****ERROR: Unrecognized parameter: " << argv[i] << " *****" << endl << endl;
			showHelp = true;
		}		
	}

	// make sure we have a BAM file
	if (!haveInBam) {
	  cerr << endl << "*****" << endl << "*****ERROR: Need -bam. " << endl << "*****" << endl;
	  showHelp = true;
	}
	// make sure we have an end1 FASTQ file
	if (!haveFastq1) {
	  cerr << endl << "*****" << endl << "*****ERROR: Need -fq1. " << endl << "*****" << endl;
	  showHelp = true;
	}
	// make sure we have an end2 FASTQ file
	if (!haveFastq2) {
	  cerr << endl << "*****" << endl << "*****ERROR: Need -fq2. " << endl << "*****" << endl;
	  showHelp = true;
	}

	
	// let 'er rip.
	if (!showHelp) {
		
		// open the 1st fastq file for writing
    	ofstream fq1(fastq1.c_str(), ios::out);
    	if ( !fq1 ) {
    		cerr << "Error: The first fastq file (" << fastq1 << ") could not be opened.  Exiting!" << endl;
    		exit (1);
    	}
    	// open the 2nd fastq file for writing
    	ofstream fq2(fastq2.c_str(), ios::out);
    	if ( !fq2 ) {
    		cerr << "Error: The second fastq file (" << fastq2 << ") could not be opened.  Exiting!" << endl;
    		exit (1);
    	}
    	
    	// open the BAM file
    	BamReader reader;	
    	reader.Open(inBamFile);
	
    	BamAlignment bam1, bam2;
    	while (reader.GetNextAlignment(bam1)) {		
            
            reader.GetNextAlignment(bam2);

            // get the seqs and quals
    		string seq1  = bam1.QueryBases;
    		string qual1 = bam1.Qualities;
    		string seq2  = bam2.QueryBases;
    		string qual2 = bam2.Qualities;

            // revcomp the seq and rev the quals si necessaire
    		if (bam1.IsReverseStrand() == true) {
    			reverseComplement(seq1);
    			reverseSequence(qual1);
    		}
    		if (bam2.IsReverseStrand() == true) {
    			reverseComplement(seq2);
    			reverseSequence(qual2);
    		}
    		
			// write the ends to their correct FASTQ output files.
			if (bam1.IsFirstMate() == true) {
				// bam1 is mate1
				fq1 << "@" << bam1.Name << "/1" << endl;
				fq1 << seq1 << endl;
				fq1 << "+" << endl;
				fq1 << qual1 << endl;
				// bam2 is mate2
				fq2 << "@" << bam2.Name << "/2" << endl;
				fq2 << seq2 << endl;
				fq2 << "+" << endl;
				fq2 << qual2 << endl;			
			}
			else {
				// bam2 is mate1
				fq1 << "@" << bam2.Name << "/1" << endl;
				fq1 << seq2 << endl;
				fq1 << "+" << endl;
				fq1 << qual2 << endl;
				// bam1 is mate2
				fq2 << "@" << bam1.Name << "/2" << endl;
				fq2 << seq1 << endl;
				fq2 << "+" << endl;
				fq2 << qual1 << endl;					
			}
		}
        reader.Close();
	}
	else {
		ShowHelp();
	}
}


void ShowHelp(void) {

	cerr << endl << "Program: " << PROGRAM_NAME << " (v" << VERSION << ")" << endl;
	
	cerr << "Author:  Aaron Quinlan (aaronquinlan@gmail.com)" << endl;

	cerr << "Summary: Creates paired-end FASTQ files from a *name sorted* BAM file" << endl << endl;

	cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -i <BAM> -fq1 <FQ> -fq2 <FQ>" << endl << endl;

	exit(1);	
}

