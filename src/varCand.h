#ifndef SRC_VARCAND_H_
#define SRC_VARCAND_H_

#include <iostream>
#include <cstring>
#include <fstream>
#include <vector>
#include <sys/stat.h>
#include <limits.h>

#include <htslib/faidx.h>

#include "FastaSeqLoader.h"
#include "Region.h"
#include "RefSeqLoader.h"
#include "alnDataLoader.h"
#include "covLoader.h"
//#include "clipReg.h"

using namespace std;

#define BLAT_ALN							1
#define ALIGN_DEBUG							0

#define MIN_VALID_BLAT_SEG_SIZE				200  // the minimal valid blat align segment size
#define MIN_VALID_BLAT_SEG_FRATCION			0.6  // the minimal valid blat align segment fraction
#define INVALID_COV_RATIO_THRES				0.5	 // the minimal coverage ratio to determine a invalid blat alignment
#define MIN_AVER_SIZE_ALN_SEG				1000

#define VAR_ALN_EXTEND_SIZE					200
#define SHORT_VAR_ALN_CHECK_EXTEND_SIZE		50
#define MIN_MISMATCH_NUM_SHORT_VAR_CALL		2
#define MIN_DISAGR_NUM_SHORT_VAR_CALL		1

#define MIN_HIGH_INDEL_BASE_RATIO			0.4

#define MATCH_SCORE							2
#define MISMATCH_SCORE						-2
#define GAP_SCORE							-2
#define GAP_OPEN_SCORE						-4

#define MARGIN_MISMATCH_NUM_IGNORE_THRES	2
#define MIN_MARGIN_LEN_MISMATCH_IGNORE		5

#define MIN_SHORT_DUP_SIZE					30		// to be parameterized
#define MIN_SHORT_INV_SIZE					30		// to be parameterized

#define SIMILARITY_THRES_INV				(0.9)	// to be parameterized

//#define MIN_REF_DIST_DUP_NUM_EST			1000

//#define MAX_SV_SIZE							5000


typedef struct{
	size_t query_start, query_end, subject_start, subject_end;
	float ident_percent;
	size_t aln_orient;
	size_t ref_start, ref_end;  // reference positions
}aln_seg_t;

typedef struct{
	//string query_name, subject_name;
	size_t blat_aln_id, query_len, subject_len, query_id, aln_orient;
	bool head_hanging, tail_hanging, best_aln, valid_aln;  // default: false
	vector<aln_seg_t*> aln_segs;
} blat_aln_t;

typedef struct{
	reg_t *reg, *cand_reg;
	aln_seg_t *aln_seg;
	int32_t startRefPos, endRefPos, startLocalRefPos, endLocalRefPos, startQueryPos, endQueryPos;
	string ctgseq, refseq;
	vector<string> alignResultVec;		// [0]: aligned ctgseq, [1]: match character, [2]: aligned refseq
	int32_t overlapLen, queryLeftShiftLen, queryRightShiftLen, localRefLeftShiftLen, localRefRightShiftLen;
	int32_t start_aln_idx_var, end_aln_idx_var;
}localAln_t;


class varCand {
	public:
		string var_cand_filename, out_dir_call, chrname, inBamFile, misAln_filename;
		faidx_t *fai;

		string refseqfilename, ctgfilename, readsfilename, alnfilename;
		int32_t ref_left_shift_size, ref_right_shift_size, ctg_num;
		vector<reg_t*> varVec, newVarVec;
		bool assem_success, align_success, call_success, clip_reg_flag;  	// default: false
		vector<blat_aln_t*> blat_aln_vec;               	// blat aligned segments

		// clippings
		reg_t *clip_reg; //, *clip_reg_mate;
		size_t leftClipRefPos, rightClipRefPos, sv_type, dup_num; // leftClipQueryPos, rightClipQueryPos, aln_orient;
		bool margin_adjusted_flag;
//		int32_t sv_len;
//		string chrname, refseq_var, altseq_var;

		//mateClipReg_t *mate_clip_reg;
		//int32_t clip_end_flag;  // -1: unassigned; 0: both ends; 1: left end; 2: right end

	public:
		varCand();
		virtual ~varCand();
		void callVariants();
		vector<int32_t> computeDisagreeNumAndHighIndelBaseNum(string &chrname, size_t startRefPos, size_t endRefPos, string &inBamFile, faidx_t *fai);
		void fillVarseq();
		void loadBlatAlnData();

	private:
		void init();
		void alnCtg2Refseq();;
		void blatAln(string &alnfilename, string &contigfilename, string &refseqfilename);
		void assignBlatAlnStatus();
		void callIndelVariants();
		void callClipRegVariants();
		void blatParse();
		int32_t getQueryNameLoc(string query_name, vector<string> query_name_vec);
		void assignBestAlnFlag();
		void assignAlnOrientFlag();
		void blatFilter();
		void blatFilterByShortAlnSegs();
		void blatFilterByGenomeCoverage();
		void blatFilterByQueryCoverage();
		//void removeInvalidAlnSegs();
		float computeRepeatCovRatio(blat_aln_t *blat_aln, int8_t *cov_array, bool query_flag);
		void updateCovArray(blat_aln_t *blat_aln, int8_t *cov_array, bool query_flag);
		void determineIndelType();
		void callShortVariants();
		bool isUnsplitAln(reg_t *reg);
		bool determineQueryidReg(reg_t *reg);
		aln_seg_t *getAlnSeg(reg_t *reg);
		void computeLocalLocsAlnShortVar(localAln_t *local_aln);
		void computeSeqAlignment(localAln_t *local_aln);
		void computeVarLoc(localAln_t *local_aln);
		void confirmShortVar(localAln_t *local_aln);
		int32_t getMismatchNumAln(vector<string> &alignResultVec, int32_t start_check_idx, int32_t end_check_idx);
		int32_t getDisagreeNum(Base *baseArray, int32_t arr_size);
		int32_t computeHighIndelBaseNum(Base *baseArray, int32_t arr_size, float threshold);
		void adjustVarLoc(localAln_t *local_aln);
		int32_t getAdjustedStartAlnIdxVar(localAln_t *local_aln);
		int32_t getAdjustedEndAlnIdxVar(localAln_t *local_aln);
		void adjustLeftLocShortVar(localAln_t *local_aln, int32_t newStartAlnIdx);
		void adjustRightLocShortVar(localAln_t *local_aln, int32_t newEndAlnIdx);
		void processInnerContainedRegs(vector<reg_t*> &foundRegVec, vector<reg_t*> &candRegVec);
		void removeVariantsDuplicatedContigs(vector<reg_t*> &foundRegVec, vector<reg_t*> &candRegVec);
		void mergeNeighboringVariants(vector<reg_t*> &foundRegVec, vector<reg_t*> &candRegVec);
		vector< vector<reg_t*> > dealWithTwoVariantSets(vector<reg_t*> &foundRegVec, vector<reg_t*> &candRegVec);
		void computeVarType(reg_t *reg);
		void updateVarVec(vector<vector<reg_t*>> &regVec, vector<reg_t*> &foundRegVec, vector<reg_t*> &varVec);
		vector<int32_t> computeDisagreeNumAndHighIndelBaseNumAndMarginDist(string &chrname, size_t startRefPos, size_t endRefPos, int32_t query_id, size_t startQueryPos, size_t endQueryPos, size_t aln_orient, string &inBamFile, faidx_t *fai);
		vector<int32_t> computeVarMargins(Base *baseArray, int32_t arr_size, float threshold);
		vector<int32_t> confirmVarMargins(int32_t left_dist, int32_t right_dist, string &chrname, size_t startRefPos, size_t endRefPos, int32_t query_id, size_t startQueryPos, size_t endQueryPos, size_t aln_orient, faidx_t *fai);
		void computeVarRegLoc(reg_t *reg, reg_t *reg_tmp);
		void computeLocalLocsAln(localAln_t *local_aln);
		void distinguishShortDupInvFromIndels();
		double computeSimilaritySeqAln(localAln_t *local_aln);

		// clippings
		void determineClipRegVarType();
		void determineClipRegDupType();
		void determineClipRegInvType();
		int32_t getInvBlatAlnItemIdx(reg_t *reg, vector<blat_aln_t*> &blat_aln_vec, size_t given_idx, aln_seg_t *seg1, aln_seg_t *seg2);
		void updateVarVec();
		reg_t* computeClipPos(blat_aln_t *blat_aln, aln_seg_t *seg1, aln_seg_t *seg2, string &refseq, string &queryseq, size_t var_type);
		vector<size_t> computeLeftShiftSizeDup(reg_t *reg, aln_seg_t *seg1, aln_seg_t *seg2, string &refseq, string &queryseq);
		vector<size_t> computeRightShiftSizeDup(reg_t *reg, aln_seg_t *seg1, aln_seg_t *seg2, string &refseq, string &queryseq);
		size_t computeMismatchNumLocalAln(localAln_t *local_aln);
		vector<size_t> computeQueryClipPosDup(blat_aln_t *blat_aln, size_t leftClipRefPos, string &refseq, string &queryseq);

		// output
		void printSV();
		void outputNewVarcandVec();
		void outputMultiBlatAlnBlocks();
};

#endif /* SRC_VARCAND_H_ */
