/*****************************************************************************
  clipper.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "version.h"
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/BamAux.h"
#include "bedFile.h"
#include "Align.h"
using namespace BamTools;

#include <vector>
#include <algorithm>    // for swap()
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <assert.h>

using namespace std;


// define our program name
#define PROGRAM_NAME "merger"

// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)


// function declarations
void ShowHelp(void);

void MergePairs(const string &bamFile, int minOverlap, bool verbose);
string makeOverlapConsensus (const string aln1, const string aln2, const string qual1, const string qual2);

int main(int argc, char* argv[]) {

    // our configuration variables
    bool showHelp = false;

    string bamFile = "stdin";
    bool haveBam    = true;
    int  minOverlap = 10;
    bool verbose    = false;


    // check to see if we should print out some help
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

        if(PARAMETER_CHECK("-i", 2, parameterLength)) {
            if ((i+1) < argc) {
                bamFile = argv[i + 1];
                i++;
            }
        }
        else if(PARAMETER_CHECK("-mo", 3, parameterLength)) {
            if ((i+1) < argc) {
                minOverlap = atoi(argv[i + 1]);
                i++;
            }
        }
        else if(PARAMETER_CHECK("-verbose", 8, parameterLength)) {
            verbose = true;
        }
        else {
            cerr << endl << "*****ERROR: Unrecognized parameter: " << argv[i] << " *****" << endl << endl;
            showHelp = true;
        }
    }

    // make sure we have an input files
    if (haveBam == false) {
        cerr << endl << "*****" << endl << "*****ERROR: Need -i (BAM) file. " << endl << "*****" << endl;
        showHelp = true;
    }
    // if there are no problems, let's convert BAM to BED or BEDPE
    if (!showHelp) {
        MergePairs(bamFile, minOverlap, verbose);
    }
    else {
        ShowHelp();
    }
}


void ShowHelp(void) {

    cerr << endl << "Program: " << PROGRAM_NAME << " (v" << VERSION << ")" << endl;

    cerr << "Author:  Aaron Quinlan (aaronquinlan@gmail.com)" << endl;

    cerr << "Summary: Converts BAM alignments to BED6 or BEDPE format." << endl << endl;

    cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -i <bam> " << endl << endl;

    cerr << "Options: " << endl;
    cerr << "\t-mo\t"       << "Minimum number of base pairs of overlap req'd for merging (def. 10)" << endl << endl;
    
    cerr << "NOTE:  Requires a name-sorted BAM file." << endl << endl;
    // end the program here
    exit(1);
}


void MergePairs(const string &bamFile, int minOverlap, bool verbose = false) {
    // open the BAM file
    BamReader reader;
    reader.Open(bamFile);

    // get header & reference information
    string bamHeader  = reader.GetHeaderText();
    RefVector refs    = reader.GetReferenceData();

    // rip through the BAM file and convert each mapped entry to BEDPE
    BamAlignment bam1, bam2;
    while (reader.GetNextAlignment(bam1)) {
        // slurp the second alignment for the pair.
        reader.GetNextAlignment(bam2);        
        // require that the alignments are from the same query
        if (bam1.Name == bam2.Name) {
            
            // align end1 and end2 and see if we can merge them.
            Align aln(bam1.QueryBases,bam2.QueryBases,false); // false == local alignment

            // store the important results from the alignment to avoid
            // recalling Align() methods.
            int start1 = aln.getSeq1Start();
            int start2 = aln.getSeq2Start();
            int end1   = aln.getSeq1End();
            int end2   = aln.getSeq2End();
            int len1   = end1 - start1 + 1;
            int len2   = end2 - start2 + 1;
            string aln1 = aln.getSeq1Alignment();
            string aln2 = aln.getSeq2Alignment();
            
            if (len1 >= minOverlap && start2 <= minOverlap) {

                if (verbose) {
                    printf("%s\n", bam1.QueryBases.c_str());
                    printf("%s\n", bam2.QueryBases.c_str());
                    printf("%s\n", bam1.Qualities.c_str());
                    printf("%s\n", bam2.Qualities.c_str());
                    printf("READ1: %s\t%d\t%d\n", aln.getSeq1Alignment(), start1, end1);
                    printf("       %s\n", aln.getMatchString());
                    printf("READ2: %s\t%d\t%d\n", aln.getSeq2Alignment(), start2, end2);
        
                    string pad1(start1, ' ');
                    string pad2(start1-start2, ' ');
                    cout << bam1.QueryBases.substr(1,start1) 
                         << aln1 << endl
                         << pad1 << aln.getMatchString() << endl
                         << pad2 << bam2.QueryBases.substr(1,start2)
                         << aln2 << bam2.QueryBases.substr(start2,(end2-start2+1))
                         << endl;            
                }

                // derive a consensus sequence from the overlap of end1 and end2.
                // if bases differ, we choose the one with the higher quality score.
                string overlap_consensus = makeOverlapConsensus(aln1,
                                                                aln2,
                                                                bam1.Qualities.substr(start1, len1),
                                                                bam2.Qualities.substr(start2, len2));

                // make a FASTA entry from 
                cout << ">" << bam1.Name << "-consensus" << endl;
                cout << bam1.QueryBases.substr(1,start1) 
                     << overlap_consensus 
                     << bam2.QueryBases.substr(start2,(end2-start2+1))
                     << endl;
                
                if (verbose) { printf("\n\n"); }
            }
            else {
                string mate1, mate2;
                mate1 = mate2 = "/1";
                if (!bam1.IsFirstMate()) mate1 = "/2";
                if (!bam2.IsFirstMate()) mate2 = "/2";

                cout << ">" << bam1.Name << mate1 << endl << bam1.QueryBases << endl;
                cout << ">" << bam2.Name << mate2 << endl << bam2.QueryBases << endl;
            }
        }
    }
}

string makeOverlapConsensus (const string aln1, const string aln2, const string qual1, const string qual2) {
    ostringstream overlap_consensus;
    for (size_t i = 0; i < aln1.size(); ++i) {
        if (aln1[i] == aln2[i]) {
            overlap_consensus << aln1[i];
        }
        else {
            if (qual1[i] > qual2[i]) {
                overlap_consensus << aln1[i];
            }
            else {
                overlap_consensus << aln2[i];
            }
        }
    }
    return overlap_consensus.str();
}







