#ifndef SRC_UTIL_H_
#define SRC_UTIL_H_

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sys/stat.h>
#include <htslib/sam.h>

#include "Base.h"
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
reg_t* findVarvecItem(int32_t startPos, int32_t endPos, vector<reg_t*> &varVec);
vector<reg_t*> findVarvecItemAll(int32_t startPos, int32_t endPos, vector<reg_t*> &varVec);
reg_t* findVarvecItemExtSize(int32_t startRefPos, int32_t endRefPos, vector<reg_t*> &varVec, int32_t leftExtSize, int32_t rightExtSize);
int32_t getVectorIdx(reg_t *reg, vector<reg_t*> &varVec);
reg_t* getOverlappedReg(reg_t *reg, vector<reg_t*> &varVec);
int32_t getOverlappedRegIdx(reg_t *reg, vector<reg_t*> &varVec);
bool isOverlappedReg(reg_t* reg1, reg_t* reg2);
bool isOverlappedPos(size_t startPos1, size_t endPos1, size_t startPos2, size_t endPos2);
int32_t getOverlapSize(size_t startPos1, size_t endPos1, size_t startPos2, size_t endPos2);
bool isAdjacent(size_t startPos1, size_t endPos1, size_t startPos2, size_t endPos2, int32_t dist_thres);
bool isOverlappedMateClipReg(mateClipReg_t *mate_clip_reg1, mateClipReg_t *mate_clip_reg2);
mateClipReg_t* getOverlappedMateClipReg(mateClipReg_t *mate_clip_reg_given, vector<mateClipReg_t*> &mateClipRegVec);
bam_hdr_t* loadSamHeader(string &inBamFile);
bool isInReg(int32_t pos, vector<reg_t*> &vec);
int32_t computeDisagreeNum(Base *baseArray, int32_t arr_size);
void mergeOverlappedReg(vector<reg_t*> &regVector);
void updateReg(reg_t* reg1, reg_t* reg2);
void mergeAdjacentReg(vector<reg_t*> &regVec, size_t dist_thres);
void printRegVec(vector<reg_t*> &regVec, string header);
void printMateClipReg(mateClipReg_t *mate_clip_reg);
string preprocessPipeChar(string &cmd_str);
bool isFileExist(string &filename);
void removeRedundantItems(vector<reg_t*> &reg_vec);
int32_t getLineCount(string &filename);
bool isBaseMatch(char ctgBase, char refBase);
bool isRegValid(reg_t *reg);
void exchangeRegLoc(reg_t *reg);
void blatAln(string &alnfilename, string &contigfilename, string &refseqfilename);

class Time{
	private:
		time_t timestamp, start, end, start_all, end_all;
		struct tm *time_tm;
		struct timeval time_val;

	public:
		Time();
		virtual ~Time();
		string getTime();
		void printTime();
		void setStartTime();
		void printSubCmdElapsedTime();
		void printOverallElapsedTime();
};

#endif /* SRC_UTIL_H_ */
