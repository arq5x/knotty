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
#include "api/BamConstants.h"
#include "bedFile.h"
using namespace BamTools;

#include <vector>
#include <algorithm>    // for swap()
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <assert.h>

using namespace std;


// define our program name
#define PROGRAM_NAME "clipper"

// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)


// function declarations
void ShowHelp(void);


bool GetSoftClips(const vector<CigarOp> &cigar, int alignmentStart,
                                                vector<int> &clipSizes,
                                                vector<int> &readPositions,
                                                vector<int> &genomePositions);
int  GetSoftClipLocation(const BamAlignment &bam);
void ExtractSplitsAndDisordants(const string &bamFile, int minClipSize);

int main(int argc, char* argv[]) {

    // our configuration variables
    bool showHelp = false;

    // input files
    string bamFile  = "stdin";
    int minClipSize = 10;
    bool haveBam    = true;


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
        else if(PARAMETER_CHECK("-mc", 3, parameterLength)) {
            if ((i+1) < argc) {
                minClipSize = atoi(argv[i + 1]);
                i++;
            }
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
        ExtractSplitsAndDisordants(bamFile, minClipSize);
    }
    else {
        ShowHelp();
    }
}


void ShowHelp(void) {

    cerr << endl << "Program: " << PROGRAM_NAME << " (v" << VERSION << ")" << endl;

    cerr << "Author:  Aaron Quinlan (aaronquinlan@gmail.com)" << endl;

    cerr << "Summary: Extracts discordant and clipped BAM alignments." << endl << endl;

    cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -i <bam> " << endl << endl;

    cerr << "Options: " << endl;
    cerr << "\t-mc\t"   << "Minimum number of base pairs of clipping req'd (def. 10)" << endl << endl;

    cerr << "NOTE:  Requires a name-sorted BAM file." << endl << endl;
    // end the program here
    exit(1);
}


void ExtractSplitsAndDisordants(const string &bamFile, int minClipSize) {
    // open the BAM file
    BamReader reader;
    reader.Open(bamFile);

    // get header & reference information
    string bamHeader  = reader.GetHeaderText();
    RefVector refs    = reader.GetReferenceData();

    // create output BAMs for each case
    BamWriter discordants;
    BamWriter softclips;

    // name the output BAMs
    string discFile  = bamFile + ".disc.bam";
    string clipFile  = bamFile + ".clip.bam";

    // open our BAM writer
    discordants.Open(discFile, bamHeader, refs);
    softclips.Open(clipFile, bamHeader, refs);

    // NOTE: Assumes input BAM is sorted by query
    BamAlignment bam1, bam2;
    while (reader.GetNextAlignment(bam1)) {
        reader.GetNextAlignment(bam2);
        
        if (bam1.Name != bam2.Name) {
            while (bam1.Name != bam2.Name)
            {
                bam1 = bam2;
                reader.GetNextAlignment(bam2);
            }
        }
            
        // exlude duplicate reads.
        if (!bam1.IsDuplicate() && !bam2.IsDuplicate()) {
            // count the number of soft clips for each end
            vector<int> clipSizes1, readPositions1, genomePositions1;
            vector<int> clipSizes2, readPositions2, genomePositions2;
            bool hasSoftClip1 = bam1.GetSoftClips(clipSizes1, readPositions1, genomePositions1);
            bool hasSoftClip2 = bam2.GetSoftClips(clipSizes2, readPositions2, genomePositions2);

            // plain old proper pair without a clip.  move on.
            if (bam1.IsProperPair() && !hasSoftClip1 && !hasSoftClip2) 
                continue;
            
            // the end has a single soft clipping and it 
            // affects enough base pairs in the read (i.e., it's likely not just the last base)
            if (bam1.IsProperPair() && (hasSoftClip1 || hasSoftClip2)) {
                // keep this alignment if at least one 
                // of the soft clips was big enough.
                int maxClip1 = 0;
                int maxClip2 = 0;
                if (hasSoftClip1) {
                    vector<int>::const_iterator it1 = std::max_element(clipSizes1.begin(), clipSizes1.end());
                    maxClip1 = *it1;
                }
                if (hasSoftClip2) {
                    vector<int>::const_iterator it2 = std::max_element(clipSizes2.begin(), clipSizes2.end());
                    maxClip2 = *it2;
                }
                                
                if (maxClip1 >= minClipSize || maxClip2 >= minClipSize) {
                    softclips.SaveAlignment(bam1);
                    softclips.SaveAlignment(bam2);
                }
            }
            // discordant or orphaned
            else if (!bam1.IsProperPair()) {
                discordants.SaveAlignment(bam1);
                discordants.SaveAlignment(bam2);
            }
        }
    }
    reader.Close();
    softclips.Close();
    discordants.Close();
}


int GetSoftClipLocation(const BamAlignment &bam) {

    // this assumes that there is only one soft clip.
    // either as the first (left clip) or last (right clip)
    // CIGAR operation
    int location = 1;
    if (bam.CigarData[0].Type == 'S') // left clip, the bp is at the alignment start
    {
        location = bam.Position;
    }
    else if (bam.CigarData[bam.CigarData.size() - 1].Type == 'S') // right clip, the bp is at the alignment end
    {
        location =bam.GetEndPosition(false, false);
    }
    return location;
}






