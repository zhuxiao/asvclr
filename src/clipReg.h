#ifndef SRC_CLIPREG_H_
#define SRC_CLIPREG_H_

#include <iostream>
#include <string>
#include <vector>
#include <limits.h>
#include <htslib/sam.h>

#include "structures.h"
#include "alnDataLoader.h"
#include "varCand.h"
#include "covLoader.h"

using namespace std;


#define CLIP_END_EXTEND_SIZE				300 // 100
#define MIN_CLIP_END_SIZE					200	// 50
#define MIN_ALN_SIZE_SAME_ORIENT			200
#define MIN_ALN_SIZE_DIFF_ORIENT			50

#define MIN_CLIP_DIST_THRES					1000
#define GROUP_DIST_THRES					1000

#define LONGER_NON_SA_SEG_NUM_THRES			3  // to be parameterized
#define LONGER_NON_SA_SEG_RATIO_THRES		(0.1f)
//#define MIN_VALID_CLIP_POS_RATIO			0.7

//#define FAKE_CLIP_RATIO_THRES				(0.2f)

#define CLIP_DIFF_LEN_RATIO_SV				(0.05f)

#define LEFT_END							0
#define RIGHT_END							1

#define MIN_QUERY_SELF_OVERLAP_SIZE			100		// 100
//#define MAX_SAME_REG_THRES_SAME_ORIENT		5000	// 500, use paras.maxVarRegSize instead

#define MAX_DIST_SAME_CLIP_END				100000

#define MAX_ALN_SEG_NUM_PER_READ_TRA		4


typedef struct{
	reg_t *leftClipReg, *leftClipReg2, *rightClipReg, *rightClipReg2;
	int32_t leftClipPosNum, leftClipPosNum2, rightClipPosNum, rightClipPosNum2;
	int64_t leftMeanClipPos, leftMeanClipPos2, rightMeanClipPos, rightMeanClipPos2;
	int8_t leftClipRegNum, rightClipRegNum;

	int32_t sv_type:8, dup_num:24;
	bool reg_mated_flag, valid_flag, call_success_flag, tra_rescue_success_flag;
	varCand *var_cand, *left_var_cand_tra, *right_var_cand_tra;  // TRA
	string chrname_leftTra1, chrname_leftTra2, chrname_rightTra1, chrname_rightTra2;
	int32_t leftClipPosTra1, leftClipPosTra2, rightClipPosTra1, rightClipPosTra2;
	string refseq_tra, altseq_tra, refseq_tra2, altseq_tra2;
	string bnd_mate_reg_strs[4]; // mate strings for BND format: mate_reg_id1|clip_loc1|sup_num1|cov1,mate_reg_id2|clip_loc2|sup_num2|cov2;......
}mateClipReg_t;

class clipReg {
	public:
		Paras *paras;
		string chrname, inBamFile;
		int64_t chrlen, startRefPos, endRefPos;
		int32_t minClipEndSize, maxVarRegSize, minAlnSize_same_orient, minAlnSize_diff_orient;
		faidx_t *fai;

		mateClipReg_t mate_clip_reg;

		vector<clipAlnData_t*> clipAlnDataVector, clipAlnDataVector2, rightClipAlnDataVector, rightClipAlnDataVector2;
		vector<clipPos_t*> leftClipPosVector, rightClipPosVector;
		vector<clipPos_t*> leftClipPosVector2, rightClipPosVector2;
		bool left_part_changed, right_part_changed;

		//string bnd_mate_reg_strs[4]; // mate region strings for BND format

		// deal with ultra-high coverage region


	public:
		clipReg(string &chrname, int64_t startRefPos, int64_t endRefPos, string &inBamFile, faidx_t *fai, Paras *paras);
		virtual ~clipReg();
		void computeMateClipReg();

	private:
		void destroyClipAlnDataVector(vector<clipAlnData_t*> &clipAlnDataVec);
		void destroyClipPosVec(vector<clipPos_t*> &clipPosVec);
		void fillClipAlnDataVectorWithSATag();
		void removeNonclipItems();
		void removeNonclipItemsOp(vector<clipAlnData_t*> &clipAlnDataVector);
		bool isValidClipReg();
		void computeMateAlnClipReg();
		void extractClipPosVec();
		void splitClipPosVec();
		int32_t getClipPosVecId(clipAlnData_t *clip_aln, int32_t clip_end);
		clipPos_t* getClipPosItemFromSingleVec(clipAlnData_t *clip_aln, int32_t clip_end, vector<clipPos_t*> &clip_pos_vec);
		void appendClipPos();
		void appendClipPosSingleVec(vector<clipPos_t*> &clipPosVec, vector<clipAlnData_t*> &clipAlnDataVec);
		vector<string> getClipPosMaxOcc(vector<clipPos_t*> &clipPosVec);
		void addClipPosAtClipEnd(string chrname_clip, int64_t start_pos, int64_t end_pos, int32_t clip_end, vector<clipPos_t*> &clipPosVec, vector<clipAlnData_t*> &clip_aln_data_vec, vector<clipAlnData_t*> &mainClipAlnDataVec);
		clipAlnData_t* getClipAlnItemFromMainVec(clipAlnData_t *clip_aln, vector<clipAlnData_t*> &mainClipAlnDataVec);
		void sortClipPos();
		void sortClipPosSingleVec(vector<clipPos_t*> &clipPosVector);
		void removeFakeClips();
		void removeFakeClipsDifferentChrSingleVec(vector<clipPos_t*> &clipPosVector);
		void removeFakeClipsLongDistSameOrientSingleVec(vector<clipPos_t*> &clipPosVector, string &vec_name);
		void removeFakeClipsLowCov(vector<clipPos_t*> &clipPosVector, int32_t min_clip_reads_num);
		string computeMateSingleRegStrForBND(vector<clipPos_t*> &clip_pos_vec, int32_t vec_id, string &chrname_clip, int64_t meanClipPos);
		int32_t computeCovNumClipPos(string &chrname, int64_t meanClipPos, int32_t clip_end, faidx_t *fai, Paras *paras);
		void computeClipRegs();
		reg_t* computeClipRegSingleVec(vector<clipPos_t*> &clipPosVector);
		int32_t getItemIdxClipPosVec(clipPos_t *item, vector<clipPos_t*> &vec);
		void removeFalseOverlappedMateClipReg();
		size_t computeMeanClipPos(vector<clipPos_t*> &clipPosVector);
		void removeFPClipSingleEnd(mateClipReg_t &mate_clip_reg);
		void checkLocOrder(mateClipReg_t &mate_clip_reg);
		void computeVarTypeClipReg(mateClipReg_t &mate_clip_reg, string &inBamFile);
		void resetClipCheckFlag(vector<clipAlnData_t*> &clipAlnDataVector);
		bool isSameChrome(vector<clipAlnData_t*> &query_aln_segs);
		bool isSameOrient(vector<clipAlnData_t*> &query_aln_segs);
//		bool isQuerySelfOverlap(vector<clipAlnData_t*> &query_aln_segs);
//		bool isSegSelfOverlap(clipAlnData_t *clip_aln1, clipAlnData_t *clip_aln2);
		bool isSameAlnReg(vector<clipAlnData_t*> &query_aln_segs);
		void sortQueryAlnSegs(vector<clipAlnData_t*> &query_aln_segs);
		bool isAdjacentClipAlnSeg(clipAlnData_t *clip_aln1, clipAlnData_t *clip_aln2, size_t dist_thres);
		bool containCompleteDup(vector<clipAlnData_t*> &query_aln_segs, mateClipReg_t &mate_clip_reg);
		size_t extractVarType(size_t dup_type_num, size_t inv_type_num, size_t tra_type_num, size_t min_reads_thres);
		size_t computeDupNumClipReg(vector<size_t> &dup_num_vec);
		void printClipVecs(string head_info);
		void printSingleClipVec(vector<clipPos_t*> &clip_pos_vec);
		void printResultClipRegs();
};

#endif /* SRC_CLIPREG_H_ */
