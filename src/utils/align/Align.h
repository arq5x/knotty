#include <string>
#include <stdlib.h>
#include "stdaln.h"
#include <iostream>


using namespace std;

class Align {

public:
    // ctor. alignment of seq1 and seq2 is done 
    // automatically as part of construction
    Align(const string &seq1, const string &seq2, bool global);
    ~Align(void);
    
    // access the aligned portions of each sequence. example below is global
    // ---------------------------------------------------------------------
    char* getSeq1Alignment(void);   // aatctg-   // aligned portion of seq1
    char* getMatchString(void);     // |||| |    // "|" for matches
    char* getSeq2Alignment(void);   // aatcgga   // aligned portion of seq2
    
    int getScore(void);     // return the score for the alignment
    
    int getSeq1Start(void); // return seq1's alignment start coordinate
    int getSeq2Start(void); // return seq2's alignment start coordinate
    int getSeq1End(void);   // return seq1's alignment end coordinate
    int getSeq2End(void);   // return seq2's alignment end coordinate
    void getCigar(void);     // return CIGAR
    
private:
    void alignSequences(void);
    
private:
    AlnAln   *_alignment;
    //AlnParam *_parameters;
    string _seq1;
    string _seq2;
    int _len1;
    int _len2;
    bool _global;
};
