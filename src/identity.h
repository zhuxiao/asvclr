#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <cmath>  
#include <limits.h>
#include <unordered_set>
#include <cstdlib>
#include <algorithm>
#include <fstream>
#include <unordered_map>
#include "util.h"

using namespace std;

#define WINDOWSIZE          10
#define KMERSIZE            15

#define MAX_SEQ_LEN			5000

#define MATCHSCORE			2
#define	MISMATCHSCORE		-1
#define GAPPENALTY			-1

#define RELIEF_FACTOR		0.1


struct Minimizer {
    size_t position;
    size_t hashValue;
    string sequence;
    bool matchflag;

    Minimizer(size_t pos, size_t hash, const string& seq, bool match)
        : position(pos), hashValue(hash), sequence(seq), matchflag(match) {}
};

int max(int a, int b, int c);
double computeVarseqIdentity(const string& seq, const string& seq1);
void needleman_wunschOp(const string& seq1, const string& seq2, int32_t match_score, int32_t mismatch_score, int32_t gap_penalty, string& seq1_new, string& seq2_new);
double calculate_identity(const string& seq1, const string& seq2);
vector<Minimizer> findMinimizers(const string& sequence, size_t windowSize, size_t kmerSize);
size_t customHashFunction(const string& kmer, size_t kmerSize);
void retainCommonHashValuesMinimizers(vector<Minimizer>& containerA, vector<Minimizer>& containerB);
void findSimilarityPosMinimizers(vector<Minimizer>& minimizers, vector<Minimizer>& minimizers1);
void FilterDissimilaPosMinimizers(vector<Minimizer>& minimizers, vector<Minimizer>& minimizers1);
void LongSequenceSplitAlignment(vector<Minimizer>& minimizers, vector<Minimizer>& minimizers1, int k, string& AlignSeq, string& AlignSeq1, string sequence, string sequence1);
//void Readfile(const string file_path, string& sequence);
