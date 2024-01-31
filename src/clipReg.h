#ifndef SRC_CLIPREG_H_
#define SRC_CLIPREG_H_

#include <iostream>
#include <string>
#include <vector>
#include <limits.h>
#include <htslib/sam.h>
//#include <pthread.h>

#include "structures.h"
#include "alnDataLoader.h"
#include "varCand.h"
#include "covLoader.h"
#include "Block.h"

using namespace std;


#define CLIP_END_EXTEND_SIZE				50 // 100  300
#define MIN_CLIP_END_SIZE					200	// 50
#define MIN_ALN_SIZE_SAME_ORIENT			200
#define MIN_ALN_SIZE_DIFF_ORIENT			50

#define MIN_CLIP_DIST_THRES					1000
#define GROUP_DIST_THRES					1000

//#define LONGER_NON_SA_SEG_NUM_THRES			5  // to be parameterized
#define LONGER_NON_SA_SEG_RATIO_THRES		(0.1f)
//#define MIN_VALID_CLIP_POS_RATIO			0.7

//#define FAKE_CLIP_RATIO_THRES				(0.2f)

#define CLIP_DIFF_LEN_RATIO_SV				(0.05f)

#define LEFT_END							0
#define RIGHT_END							1

#define MIN_QUERY_SELF_OVERLAP_SIZE			100		// 100
//#define MAX_SAME_REG_THRES_SAME_ORIENT		5000	// 500, use paras.maxVarRegSize instead

#define MAX_DIST_SAME_CLIP_END				100000

#define MAX_ALN_SEG_NUM_PER_READ_TRA		5  // to be parameterized, 4
#define MAX_CLIP_REG_MERGE_DIST				50  // to be parameterized, 50
#define MAX_CLIP_REG_MERGE_DIST_ADJUST		500  //self
#define MAX_INNER_MISSING_IGNORE_SIZE		100

#define MIN_CLIP_REG_MATED_RATIO			(0.1f) // 0.3

#define MAX_REF_DIST_SAME_CHR				500


class clipReg {
	public:
		Paras *paras;
		string chrname, inBamFile;
		int64_t chrlen, startRefPos, endRefPos;
		int32_t minClipEndSize, maxVarRegSize, minAlnSize_same_orient, minAlnSize_diff_orient;
		faidx_t *fai;

		mateClipReg_t mate_clip_reg;

		vector<clipAlnData_t*> clipAlnDataVector, clipAlnDataVector2, rightClipAlnDataVector, rightClipAlnDataVector2;
		vector<clipPos_t*> leftClipPosVector, leftClipPosVector2, rightClipPosVector, rightClipPosVector2;
		bool left_part_changed, right_part_changed, large_indel_flag;

		// large indels
		reg_t *largeIndelClipReg;
		int32_t supp_num_largeIndel, depth_largeIndel;
		vector<clipPos_t*> largeIndelClipPosVector;
//		vector<Block*> *blockVector;
//		pthread_mutex_t *p_mutex_mate_clip_reg;

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
		int32_t getMinDistVecId(clipPos_t *clip_pos_item);
		int32_t getMeanDistSingleVec(clipPos_t *clip_pos_item, vector<clipPos_t*> &clip_pos_vec);
		clipPos_t* getClipPosSameAlnSeg(clipPos_t *clip_pos_item, vector<clipPos_t*> &clip_pos_vec);
		void splitClipPosVec();
		int32_t computeBPConsistencyFlag(clipAlnData_t *clip_aln, int32_t clip_end, clipAlnData_t *adj_aln, int32_t adj_clip_end, vector<clipPos_t*> &clip_pos_vec);
		void AdjustClipPosVecByBPConsistency(vector<clipPos_t*> &clip_pos_vec, int32_t vec_id);
		void AdjustClipPosVecByBPConsistencySingleVec(vector<clipPos_t*> &clip_pos_vec, vector<clipPos_t*> &clip_pos_vec_main);
		vector<clipPos_t*> getClipPosItemsByQueryname(string &queryname, vector<clipPos_t*> &clip_pos_vec);
		clipPos_t* getMinDistClipPosItem(clipPos_t *clip_pos, vector<clipPos_t*> &clip_pos_vec);
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
		void determineLargeIndelReg();
		reg_t* computeClipRegSingleVec(vector<clipPos_t*> &clipPosVector);
		int32_t getItemIdxClipPosVec(clipPos_t *item, vector<clipPos_t*> &vec);
		void removeBNDUnmatedClipRegs();
		string getBNDMateRegStr(int32_t reg_idx, int32_t end_flag, string bnd_mate_reg_strs[4]);
		vector<int32_t> getMateMateRegID(int32_t reg_idx, int32_t end_flag, string bnd_mate_reg_strs[4]);
		void removeFalseOverlappedMateClipReg();
		size_t computeMeanClipPos(vector<clipPos_t*> &clipPosVector);
		void removeFPClipSingleEnd(mateClipReg_t &mate_clip_reg);
		void checkLocOrder(mateClipReg_t &mate_clip_reg);
		void computeVarTypeClipReg(mateClipReg_t &mate_clip_reg);
		void resetClipCheckFlag(vector<clipAlnData_t*> &clipAlnDataVector);
		bool isSameChrome(vector<clipAlnData_t*> &query_aln_segs);
		bool isSameOrient(vector<clipAlnData_t*> &query_aln_segs);
		bool ischeckSameOrient(vector<clipAlnData_t*> &query_aln_segs);//self
		bool Filteredbychrname(vector<clipAlnData_t*> &query_aln_segs);//self
		void FilteredbyAlignmentSegment(vector<clipAlnData_t*> &query_aln_segs);//self
		bool isSameAlnReg(vector<clipAlnData_t*> &query_aln_segs);
		void sortQueryAlnSegs(vector<clipAlnData_t*> &query_aln_segs);
		bool isAdjacentClipAlnSeg(clipAlnData_t *clip_aln1, clipAlnData_t *clip_aln2, size_t dist_thres);
		bool containCompleteDup(vector<clipAlnData_t*> &query_aln_segs, mateClipReg_t &mate_clip_reg);
		size_t extractVarType(mateClipReg_t &mate_clip_reg, size_t dup_type_num, size_t inv_type_num, size_t tra_type_num, size_t min_reads_thres);
		size_t computeDupNumClipReg(vector<size_t> &dup_num_vec);
		void printClipVecs(string head_info);
		void printSingleClipVec(vector<clipPos_t*> &clip_pos_vec);
		void printResultClipRegs();
};

#endif /* SRC_CLIPREG_H_ */
