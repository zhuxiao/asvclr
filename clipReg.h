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


#define CLIP_END_EXTEND_SIZE				100
#define MIN_CLIP_END_SIZE					200	// 50

#define MIN_CLIP_DIST_THRES					1000

#define LONGER_NON_SA_SEG_NUM_THRES			3  // to be parameterized
#define LONGER_NON_SA_SEG_RATIO_THRES		(0.1f)
//#define MIN_VALID_CLIP_POS_RATIO			0.7

//#define FAKE_CLIP_RATIO_THRES				(0.2f)

#define CLIP_DIFF_LEN_RATIO_SV				(0.05f)

#define LEFT_END							0
#define RIGHT_END							1

#define MIN_QUERY_SELF_OVERLAP_SIZE			100		// 100
//#define MAX_SAME_REG_THRES_SAME_ORIENT		5000	// 500, use paras.maxClipRegSize instead

#define MAX_DIST_SAME_CLIP_END				100000


typedef struct{
	string chrname;
	size_t clipRefPos, clipLocalRefPos, clipQueryPos, aln_orient;
	bool same_orient_flag;  // true: ++, --; false: -+, +-
}clipPos_t;

//typedef struct{
//	reg_t *leftClipReg, *rightClipReg;
//	size_t leftClipPosNum, rightClipPosNum;
//	size_t leftMeanClipPos, rightMeanClipPos;
//}clipRegItem_t;

typedef struct{
	reg_t *leftClipReg, *rightClipReg;
	size_t leftClipPosNum, rightClipPosNum;
	size_t leftMeanClipPos, rightMeanClipPos;

	reg_t *leftClipReg2, *rightClipReg2;
	size_t leftClipPosNum2, rightClipPosNum2;
	size_t leftMeanClipPos2, rightMeanClipPos2;

	size_t leftClipRegNum, rightClipRegNum;

	size_t sv_type, dup_num;
	bool reg_mated_flag, valid_flag, call_success_flag;
	varCand *var_cand, *left_var_cand_tra, *right_var_cand_tra;  // TRA
	string chrname_leftTra1, chrname_rightTra1, chrname_leftTra2, chrname_rightTra2;
	int32_t leftClipPosTra1, rightClipPosTra1, leftClipPosTra2, rightClipPosTra2;
	string refseq_tra, altseq_tra, refseq_tra2, altseq_tra2;
}mateClipReg_t;

class clipReg {
	public:
		Paras *paras;
		string chrname, inBamFile;
		size_t chrlen, startRefPos, endRefPos;
		faidx_t *fai;

		mateClipReg_t mate_clip_reg;

		vector<clipAlnData_t*> clipAlnDataVector;
		vector<clipPos_t> leftClipPosVector, rightClipPosVector;
		vector<clipPos_t> leftClipPosVector2, rightClipPosVector2;

	public:
		clipReg(string &chrname, size_t startRefPos, size_t endRefPos, size_t chrlen, string &inBamFile, faidx_t *fai, Paras *paras);
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
		void splitClipPosVec();
		void removeFakeClips();
		void removeFakeClipsDifferentChrSingleVec(vector<clipPos_t> &clipPosVector);
		void removeFakeClipsLongDistSameOrientSingleVec(vector<clipPos_t> &clipPosVector);
		void removeFakeClipsLowCov(vector<clipPos_t> &clipPosVector, size_t min_clip_reads_num);
		void computeClipRegs();
		reg_t* computeClipRegSingleVec(vector<clipPos_t> &clipPosVector);
		int32_t getItemIdxClipPosVec(clipPos_t &item, vector<clipPos_t> &vec);
		void removeFalseOverlappedMateClipReg();
		size_t computeMeanClipPos(vector<clipPos_t> &clipPosVector);
		void removeFPClipSingleEnd(mateClipReg_t &mate_clip_reg);
		void checkLocOrder(mateClipReg_t &mate_clip_reg);
		void computeVarTypeClipReg(mateClipReg_t &mate_clip_reg, string &inBamFile);
		void resetClipCheckFlag(vector<clipAlnData_t*> &clipAlnDataVector);
		bool isSameChrome(vector<clipAlnData_t*> &query_aln_segs);
		bool isSameOrient(vector<clipAlnData_t*> &query_aln_segs);
		bool isQuerySelfOverlap(vector<clipAlnData_t*> &query_aln_segs);
		bool isSegSelfOverlap(clipAlnData_t *clip_aln1, clipAlnData_t *clip_aln2);
		bool isSameAlnReg(vector<clipAlnData_t*> &query_aln_segs);
		void sortQueryAlnSegs(vector<clipAlnData_t*> &query_aln_segs);
		bool isAdjacentClipAlnSeg(clipAlnData_t *clip_aln1, clipAlnData_t *clip_aln2, size_t dist_thres);
		bool containCompleteDup(vector<clipAlnData_t*> &query_aln_segs, mateClipReg_t &mate_clip_reg);
		size_t extractVarType(size_t dup_type_num, size_t inv_type_num, size_t tra_type_num, size_t min_reads_thres);
		size_t computeDupNumClipReg(vector<size_t> &dup_num_vec);
		void printClipVecs();
		void printSingleClipVec(vector<clipPos_t> &clip_pos_vec);
		void printResultClipRegs();
};

#endif /* SRC_CLIPREG_H_ */
