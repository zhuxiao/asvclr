#ifndef SRC_UTIL_H_
#define SRC_UTIL_H_

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sys/stat.h>
#include <htslib/sam.h>

#include "Region.h"
#include "clipReg.h"

using namespace std;

vector<string> split(const  string& s, const string& delim);
bool isExistStr(string &str, vector<string> &vec);
int copySingleFile(string &infilename, ofstream &outfile);
char getBase(string &seq, size_t pos, size_t orient);
void reverseSeq(string &seq);
void reverseComplement(string &seq);
void upperSeq(string &seq);
size_t getCtgCount(string &contigfilename);
reg_t* findVarvecItem(size_t startPos, size_t endPos, vector<reg_t*> &varVec);
reg_t* findVarvecItemExtSize(size_t startRefPos, size_t endRefPos, vector<reg_t*> &varVec, int32_t leftExtSize, int32_t rightExtSize);
int32_t getVectorIdx(reg_t *reg, vector<reg_t*> &varVec);
bool isOverlappedReg(reg_t* reg1, reg_t* reg2);
bool isOverlappedPos(size_t startPos1, size_t endPos1, size_t startPos2, size_t endPos2);
int32_t getOverlapSize(size_t startPos1, size_t endPos1, size_t startPos2, size_t endPos2);
bool isAdjacent(size_t startPos1, size_t endPos1, size_t startPos2, size_t endPos2, int32_t dist_thres);
bool isOverlappedMateClipReg(mateClipReg_t *mate_clip_reg1, mateClipReg_t *mate_clip_reg2);
mateClipReg_t* getOverlappedMateClipReg(mateClipReg_t *mate_clip_reg_given, vector<mateClipReg_t*> &mateClipRegVec);
bam_hdr_t* loadSamHeader(string &inBamFile);
bool isInReg(size_t pos, vector<reg_t*> &vec);


class Time{
	private:
		time_t timestamp;
		struct tm *time_tm;
		struct timeval time_val;

	public:
		Time();
		virtual ~Time();
		string getTime();
		void printTime();
};

#endif /* SRC_UTIL_H_ */
