#ifndef SRC_VARCAND_H_
#define SRC_VARCAND_H_

#include <iostream>
#include <cstring>
#include <fstream>
#include <vector>
#include <algorithm>
#include <limits.h>
#include <stdexcept>
#include <htslib/faidx.h>

#include "structures.h"
#include "FastaSeqLoader.h"
#include "Region.h"
#include "RefSeqLoader.h"
#include "alnDataLoader.h"
#include "meminfo.h"
#include "genotyping.h"
#include "clipAlnDataLoader.h"
#include "util.h"
#include "consistency.h"


using namespace std;


#define BLAT_ALN							1
#define ALIGN_DEBUG							0

#define MIN_VALID_BLAT_SEG_SIZE				200  // the minimal valid blat align segment size
#define MIN_VALID_BLAT_SEG_FRATCION			0.6  //0.6 the minimal valid blat align segment fraction
#define INVALID_COV_RATIO_THRES				0.5	 // the minimal coverage ratio to determine a invalid blat alignment
#define MIN_AVER_SIZE_ALN_SEG				200 //1000

#define VAR_ALN_EXTEND_SIZE					1000	// 200
#define SHORT_VAR_ALN_CHECK_EXTEND_SIZE		30	// 20
#define MIN_MISMATCH_NUM_SHORT_VAR_CALL		1	// 2
#define MIN_DISAGR_NUM_SHORT_VAR_CALL		1

#define MIN_HIGH_INDEL_BASE_RATIO			0.4
#define IGNORE_POLYMER_RATIO_THRES			0.6	// polymer flag will be ignored if ratio >= threshold

#define MATCH_SCORE							2
#define MISMATCH_SCORE						-2
#define GAP_SCORE							-2
#define GAP_OPEN_SCORE						-4

#define MARGIN_MISMATCH_NUM_IGNORE_THRES	2
#define MIN_MARGIN_LEN_MISMATCH_IGNORE		5

#define SV_MIN_DIST							8

#define MIN_SHORT_DUP_SIZE					30		// to be parameterized
#define MIN_SHORT_INV_SIZE					30		// to be parameterized

#define SIMILARITY_THRES_INV				(0.9f)	// to be parameterized


#define EXT_SIZE_CHK_VAR_LOC				500		// 100

#define MIN_VALID_POLYMER_SIZE				3
#define MAX_SHORT_MIS_REG_SIZE				1

#define MIN_INNER_BLAT_SEG_SIZE				50


class varCand {
	public:
		string var_cand_filename, out_dir_call, chrname, inBamFile, misAln_filename;
		faidx_t *fai;
		vector<simpleReg_t*> sub_limit_reg_vec;
		bool limit_reg_process_flag, limit_reg_delete_flag;

		string refseqfilename, ctgfilename, readsfilename, alnfilename, clusterfilename;
//		string rescue_refseqfilename, rescue_cnsfilename, rescue_readsfilename, rescue_alnfilename;
		int32_t ref_left_shift_size, ref_right_shift_size, ctg_num, min_sv_size, minReadsNumSupportSV, minClipEndSize;
		vector<reg_t*> varVec, newVarVec;
		bool cns_success, align_success, call_success, clip_reg_flag, killed_flag;  	// default: false
		vector<blat_aln_t*> blat_aln_vec;               	// blat aligned segments
		vector<minimap2_aln_t*> minimap2_aln_vec;           // minimap2 aligned segments
		//vector<clipAlnData_t*> clipAlnDataVector;

		// clippings
		reg_t *clip_reg, *clip_reg_allele; //, *clip_reg_mate;
		int64_t leftClipRefPos, rightClipRefPos;
		int32_t sv_type, dup_num, depth_largeIndel; // leftClipQueryPos, rightClipQueryPos, aln_orient;
		string bnd_mate_reg_strs[4];
		bool margin_adjusted_flag, large_indel_flag;
//		int32_t sv_len;
//		string chrname, refseq_var, altseq_var;

		//mateClipReg_t *mate_clip_reg;
		//int32_t clip_end_flag;  // -1: unassigned; 0: both ends; 1: left end; 2: right end

		// blat aligned information
		vector<varCand*> *blat_aligned_info_vec;	// previously blat aligned information

		// minimap2 aligned information
		vector<varCand*> *minimap2_aligned_info_vec;	// previously minimap2 aligned information

		// blat aligned information file
		ofstream *blat_var_cand_file;
		ofstream *minimap2_var_cand_file;

		// pairwise alignment result
		vector<localAln_t*> local_aln_vec;

		// process monitor killed blat work
		int32_t max_proc_running_minutes;
		vector<killedBlatWork_t*> *killed_blat_work_vec;
		ofstream *killed_blat_work_file;
		pthread_mutex_t *mtx_killed_blat_work;

		// process monitor killed minimap2 work
		//int32_t max_proc_running_minutes;
		vector<killedMinimap2Work_t*> *killed_minimap2_work_vec;
		ofstream *killed_minimap2_work_file;
		pthread_mutex_t *mtx_killed_minimap2_work;

		//genotyping parameter
		int32_t gt_min_sig_size, gt_min_sup_num_recover; // not used
		double gt_min_consistency_merge, gt_homo_ratio_thres, gt_hete_ratio_thres;
		double gt_size_ratio_match; // not used

		//SV position correction
		vector<svpos_correction_t*> svpos_correction_vec;

	public:
		varCand();
		virtual ~varCand();
		void callVariants();
		void callVariants02();
		void setBlatVarcandFile(ofstream *blat_var_cand_file, vector<varCand*> *blat_aligned_info_vec);
		void resetBlatVarcandFile();
		void setGtParas(int32_t gt_min_sig_size, double gt_size_ratio_match, double gt_min_consistency_merge, double gt_homo_ratio_thres, double gt_hete_ratio_thres, int32_t min_sup_num_recover);
		vector<int32_t> computeDisagreeNumAndHighIndelBaseNum(string &chrname, size_t startRefPos, size_t endRefPos, string &inBamFile, faidx_t *fai);
		void adjustVarLocSlightly();
		void fillVarseq();
		void alnCtg2Refseq02();
		void alnCtg2Refseq();
		void loadMinimap2AlnData();
		void loadBlatAlnData();
		void destroyVarCand();

	private:
		void init();
		//void alnCtg2Refseq();
		bool getMinimap2WorkKilledFlag();
		bool getBlatWorkKilledFlag();
		void recordMinimap2AlnInfo();
		void recordBlatAlnInfo();
		bool getMinimap2AlnDoneFlag();
		bool getBlatAlnDoneFlag();
		void assignMinimap2AlnStatus();
		void assignBlatAlnStatus();
		void callIndelVariants();
		void callIndelVariants02();
		void callClipRegVariants();
		void callClipRegVariants02();
		void minimap2Parse();
		vector<minimap2_aln_t*> minimap2Parse(string &alnfilename, string &cnsfilename, string &refseqfilename);
		//allocatePafInDelAlnSeg
		void getMissingPafInDelsAtAlnSegEnd(vector<struct pafalnSeg*> &all_pafalnsegs, vector<struct pafalnSeg*> &pafalnsegs, int32_t region_refstart, int32_t region_refend, int32_t querystart, int32_t queryend, int32_t query_len, int32_t subjectlen, int32_t subjectstart, int32_t subjectend, int32_t match_base_num, int32_t match_ref_len, string cons_seq, string ref_seq, int32_t aln_orient);
		void blatParse();
		void assignBestAlnFlag();
		void assignAlnOrientFlag();
		void blatFilter();
		void blatFilterByShortAlnSegs();
		void blatFilterByGenomeCoverage();
		void blatFilterByQueryCoverage();
		float computeRepeatCovRatio(blat_aln_t *blat_aln, int8_t *cov_array, bool query_flag);
		void updateCovArray(blat_aln_t *blat_aln, int8_t *cov_array, bool query_flag);
		void determineIndelType();
		void eraseFalsePositiveVariants();
		void svPosCorrection(reg_t* reg);
		vector<int32_t> computeSuppNumFromRegionAlnSegs(vector<string> &clu_qname_vec, struct pafalnSeg* paf_alnseg, string chrname, int64_t startRefPos_cns, int64_t endRefPos_cns, double size_ratio_match_thres);
		void destoryClipAlnData(vector<clipAlnData_t*> &clipAlnDataVector);
		void destoryPosCorrectionVec();
		void callShortVariants();
		bool isUnsplitAln(reg_t *reg);
		bool determineQueryidReg(reg_t *reg);
		aln_seg_t *getAlnSeg(reg_t *reg, int64_t *blat_aln_id);
		void computeExtendAlnSegs(localAln_t *local_aln);
		void computeLocalLocsAlnShortVar(localAln_t *local_aln);
		void computeSeqAlignment(localAln_t *local_aln);
		void computeSeqAlignmentOp(localAln_t *local_aln);
		void adjustAlnResult(localAln_t *local_aln);
		void computeVarLoc(localAln_t *local_aln);
		void adjustVarLoc(localAln_t *local_aln);
		//void confirmVarLoc(localAln_t *local_aln);
		void confirmShortVar(localAln_t *local_aln);
		vector<int32_t> getMismatchNumAln(string &mid_seq, int32_t start_check_idx, int32_t end_check_idx, vector<mismatchReg_t*> &misReg_vec, int32_t min_match_misReg_size);
		//int32_t getDisagreeNum(Base *baseArray, int32_t arr_size);
		int32_t computeHighIndelBaseNum(Base *baseArray, int32_t arr_size, float threshold, float polymer_ignore_ratio_thres);
		void adjustVarLocShortVar(localAln_t *local_aln);
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
		vector<int32_t> computeVarMargins(Base *baseArray, int32_t arr_size, float threshold, float polymer_ignore_ratio_thres);
		vector<int32_t> confirmVarMargins(int32_t left_dist, int32_t right_dist, string &chrname, size_t startRefPos, size_t endRefPos, int32_t query_id, size_t startQueryPos, size_t endQueryPos, size_t aln_orient, faidx_t *fai);
		int32_t getEndShiftLenFromNumVec(vector<int32_t> &numVec, size_t end_flag);
		void computeVarRegLoc(reg_t *reg, reg_t *reg_tmp);
		void computeLocalLocsAln(localAln_t *local_aln);
		void distinguishShortDupInvFromIndels();
		double computeSimilaritySeqAln(localAln_t *local_aln);
		void callVariantsAlnSegEnd();
		vector<double> computeDisagreeNumAndHighIndelBaseNumAndClipNum(string &chrname, size_t startRefPos, size_t endRefPos, string &inBamFile, faidx_t *fai);

		// clippings
		void computeClipRegVarLoc();
		vector<reg_t*> computeClipRegVarLocOp(string &alnfilename, string &refseqfilename, string &cnsfilename, string &clusterfilename, vector<minimap2_aln_t*> &minimap2_aln_vec, bool rescue_flag);
		void computeGenotypeClipReg(vector<reg_t*> &var_vec);
		vector<reg_t*> rescueDupInvClipReg();
		vector<reg_t*> rescueLargeIndelClipReg();
		vector<vector<string>> getClusterInfo(string &clusterfilename);
		vector<int32_t> getMinimapItemIdVec(vector<string> &queryname_vec, vector<minimap2_aln_t*> &minimap2_aln_vec);
		vector<int32_t> getLeftRightMinimapItemId(int64_t leftClipRefPos, int64_t rightClipRefPos, vector<int32_t> &minimap2_item_id_vec, vector<minimap2_aln_t*> &minimap2_aln_vec);
		vector<int32_t> getLeftRightClipAlnId(int64_t leftClipRefPos, int64_t rightClipRefPos, vector<clipAlnData_t*> &query_aln_segs);
		void determineClipRegVarType();
		void determineClipRegDupType();
		void determineClipRegInvType();
		vector<int32_t> getInvBlatAlnItemIdx(reg_t *reg, vector<blat_aln_t*> &blat_aln_vec, size_t given_idx, aln_seg_t *seg1, aln_seg_t *seg2);
		//reg_t* computeInvReg(reg_t *reg, int32_t blat_aln_idx, vector<int32_t> blat_aln_idx_inv_vec);
		void updateVarVec();
		reg_t* computeClipPos(blat_aln_t *blat_aln, aln_seg_t *seg1, aln_seg_t *seg2, string &refseq, string &queryseq, size_t var_type);
		reg_t* computeClipPos2(vector<blat_aln_t*> &blat_aln_vec, size_t var_type);
		vector<size_t> computeLeftShiftSizeDup(reg_t *reg, aln_seg_t *seg1, aln_seg_t *seg2, string &refseq, string &queryseq);
		vector<size_t> computeRightShiftSizeDup(reg_t *reg, aln_seg_t *seg1, aln_seg_t *seg2, string &refseq, string &queryseq);
		size_t computeMismatchNumLocalAln(localAln_t *local_aln);
		vector<size_t> computeQueryClipPosDup(blat_aln_t *blat_aln, int32_t leftClipRefPos, string &refseq, string &queryseq);

		// local alignment result -- 2021-08-22
		localAln_t* generateNewLocalAlnItem_OnlyAlnInfo(localAln_t *local_aln);
		void copyLocalAlnInInfo(localAln_t *local_aln_dest, localAln_t *local_aln_src);
		bool isLocalAlnInfoComplete(localAln_t *local_aln);
		bool isIdenticalLoclaAlnItems(localAln_t *local_aln1, localAln_t *local_aln2);
		localAln_t* getIdenticalLocalAlnItem(localAln_t *local_aln, vector<localAln_t*> &local_aln_vec);
		void addLocalAlnItemToVec(localAln_t *local_aln, vector<localAln_t*> &local_aln_vec);
		void destroyLocalAlnVec(vector<localAln_t*> &local_aln_vec);
		void destroyMinimap2AlnVec(vector<minimap2_aln_t*> &minimap2_aln_vec);

		// genotyping
		void indelGenotyping();
		void indelGenotyping02();

		// output
		void printSV();
		void outputNewVarcandVec();
		void outputMultiBlatAlnBlocks();
};

#endif /* SRC_VARCAND_H_ */
