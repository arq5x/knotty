#include "Align.h"

// ctor
Align::Align(const string &seq1, const string &seq2, bool global)
: _seq1(seq1),
  _seq2(seq2),
  _len1(seq1.size()),
  _len2(seq2.size()),
  _global(global)
{
    // trigger the alignment of se1 and seq2 right away.
    // other methods are used to retrieve specific results
    // from the alignment.
    alignSequences();
}

// dtor
Align::~Align(void)
{
    delete _alignment;
    //delete _parameters;
}

// align seq1 to seq2 either locally (Smith-Waterman) or globally (Needleman-Wunsch)
void Align::alignSequences(void) 
{
    if (_global)
        _alignment = aln_stdaln_aux(_seq1.c_str(), _seq2.c_str(), &aln_param_nt2nt, ALN_TYPE_GLOBAL, 1, _len1, _len2);
    else
        _alignment = aln_stdaln_aux(_seq1.c_str(), _seq2.c_str(), &aln_param_nt2nt, ALN_TYPE_LOCAL, 1, _len1, _len2);
}

char* Align::getSeq1Alignment(void)  { return _alignment->out1; }

char* Align::getSeq2Alignment(void)  { return _alignment->out2; }

char* Align::getMatchString(void)    { return _alignment->outm; }

int   Align::getScore(void)          { return _alignment->score; }

int   Align::getSeq1Start(void)      { return _alignment->start1; }

int   Align::getSeq2Start(void)      { return _alignment->start2; }

int   Align::getSeq1End(void)        { return _alignment->end1; }

int   Align::getSeq2End(void)        { return _alignment->end2; }

void   Align::getCigar(void) {
    for (int i = 0; i < _alignment->n_cigar; ++i) {
        printf("%d%c", _alignment->cigar32[i]>>4, "MID"[_alignment->cigar32[i]&0xf]);
    }
    printf("\n");
}
