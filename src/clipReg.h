#ifndef SRC_CLIPREG_H_
#define SRC_CLIPREG_H_

#include <iostream>
#include <string>
#include <vector>
#include <limits.h>
#include <htslib/sam.h>

#include "Region.h"
#include "alnDataLoader.h"
#include "clipAlnDataLoader.h"
#include "covLoader.h"
#include "varCand.h"

using namespace std;


#define MAX_CLIP_REG_SIZE					10000  // to be parameterized

#define CLIP_END_EXTEND_SIZE				100
#define MIN_CLIP_END_SIZE					50

#define MIN_CLIP_DIST_THRES					1000

#define LONGER_NON_SA_SEG_NUM_THRES			3  // to be parameterized
#define LONGER_NON_SA_SEG_RATIO_THRES		(0.1)
//#define MIN_VALID_CLIP_POS_RATIO			0.7

//#define FAKE_CLIP_RATIO_THRES				(0.1)

#define LEFT_END							0
#define RIGHT_END							1


typedef struct{
	string chrname;
	size_t clipRefPos, clipLocalRefPos, clipQueryPos, aln_orient;
	bool same_orient_flag;  // true: ++, --; false: -+, +-
}clipPos_t;

typedef struct{
	reg_t *leftClipReg, *rightClipReg;
	size_t leftClipPosNum, rightClipPosNum;
	size_t leftMeanClipPos, rightMeanClipPos, sv_type, dup_num;
	bool reg_mated_flag, valid_flag;
	varCand *var_cand, *left_var_cand_tra, *right_var_cand_tra;  // TRA
	int32_t leftClipPosTra1, rightClipPosTra1, leftClipPosTra2, rightClipPosTra2;
	string refseq_tra, altseq_tra;
}mateClipReg_t;

class clipReg {
	public:
		string chrname, inBamFile;
		size_t chrlen, startRefPos, endRefPos;
		faidx_t *fai;

		mateClipReg_t mate_clip_reg;

		vector<clipAlnData_t*> clipAlnDataVector;
		vector<clipPos_t> leftClipPosVector, rightClipPosVector;

	public:
		clipReg(string &chrname, size_t startRefPos, size_t endRefPos, size_t chrlen, string &inBamFile, faidx_t *fai);
		virtual ~clipReg();
		void computeMateClipReg();

	private:
		void destroyClipAlnDataVector();
		void fillClipAlnDataVectorWithSATag();
		void removeNonclipItems();
		bool isValidClipReg();
		void computeMateAlnClipReg();
		void extractClipPosVec();
		vector<clipAlnData_t*> getQueryClipAlnSegs(string &queryname, vector<clipAlnData_t*> &clipAlnDataVector);
		vector<int32_t> getAdjacentClipAlnSeg(int32_t arr_idx, int32_t clip_end_flag, vector<clipAlnData_t*> &query_aln_segs);
		void sortClipPos();
		void sortClipPosSingleVec(vector<clipPos_t> &leftClipPosVector);
		void removeFakeClips();
		void removeFakeClipsDifferentChrSingleVec(vector<clipPos_t> &clipPosVector);
		void removeFakeClipsLongDistSingleVec(vector<clipPos_t> &clipPosVector);
		void computeClipRegs();
		reg_t* computeClipRegSingleVec(vector<clipPos_t> &clipPosVector);
		int32_t getItemIdxClipPosVec(clipPos_t &item, vector<clipPos_t> &vec);
		void removeFalseOverlappedMateClipReg();
		size_t computeMeanClipPos(vector<clipPos_t> &clipPosVector);
		void computeVarTypeClipReg(mateClipReg_t &mate_clip_reg, string &inBamFile);
		void resetClipCheckFlag(vector<clipAlnData_t*> &clipAlnDataVector);
		bool isSameOrient(vector<clipAlnData_t*> &query_aln_segs);
		bool containCompleteDup(vector<clipAlnData_t*> &query_aln_segs, mateClipReg_t &mate_clip_reg);
		size_t extractVarType(size_t dup_type_num, size_t inv_type_num, size_t tra_type_num);
		size_t computeDupNumClipReg(vector<size_t> &dup_num_vec);
		void printClipVec();
		void printResultClipRegs();
};

#endif /* SRC_CLIPREG_H_ */
