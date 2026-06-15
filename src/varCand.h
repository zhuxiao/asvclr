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
#include "clipAlnDataLoader.h"
#include "seqsim.h"
#include "util.h"

using namespace std;

#define ALIGN_DEBUG							0
#define READS_CALL_DEBUG					0
#define RECOVER_OUTPUT_DEBUG				0
#define GT_REFINE_DEBUG						0
#define MINIMAP2_ALN						0


#define MIN_VALID_MINIMAP2_SEG_SIZE			200  // the minimal valid minimap2 align segment size
#define INVALID_COV_RATIO_THRES				0.5	 // the minimal coverage ratio to determine a invalid minimap2 alignment
#define MIN_AVER_SIZE_ALN_SEG				200 //1000

#define VAR_ALN_EXTEND_SIZE					1000	// 200
#define SHORT_VAR_ALN_CHECK_EXTEND_SIZE		30	// 20
#define MIN_MISMATCH_NUM_SHORT_VAR_CALL		1	// 2
#define MIN_DISAGR_NUM_SHORT_VAR_CALL		1

#define MIN_HIGH_INDEL_BASE_RATIO			0.4
#define IGNORE_POLYMER_RATIO_THRES			0.6	// polymer flag will be ignored if ratio >= threshold

#define MARGIN_MISMATCH_NUM_IGNORE_THRES	2
#define MIN_MARGIN_LEN_MISMATCH_IGNORE		5

#define SV_MIN_DIST							8

#define MIN_SHORT_DUP_SIZE					30		// to be parameterized
#define MIN_SHORT_INV_SIZE					30		// to be parameterized

#define SIMILARITY_THRES_INV				(0.9f)	// to be parameterized

#define MIN_ALN_SIZE_RATIO_HIGHSIM			(0.95f) // added on 2026-05-19

#define EXT_SIZE_CHK_VAR_LOC				500		// 100
#define EXT_SIZE_CHK_VAR_LOC_SMALL			300

#define MIN_VALID_POLYMER_SIZE				3
#define MAX_SHORT_MIS_REG_SIZE				1

#define MIN_SPAN_SINGLE_QUERY				2000

#define MAX_INV_ALN_MARGIN_DIST				5000

#define MIN_RATIO_ALLELE_CALL				(0.1f)

#define ULTRA_SHORT_SV_SIZE_FACTOR			(0.8f)  // 0.6f

#define VALID_DUP_INV_SIZE_RATIO_CLIP_REG	(0.75f)	// 0.3 (deleted on 2026-04-14)

#define MIN_SIZE_RATIO_RECOVER				(0.8f)  // added on 2026-05-17
#define MIN_ALN_SIZE_RATIO_RECOVER			(0.95f) // added on 2026-05-19


class varCand {
	public:
		string var_cand_filename, out_dir_call, out_dir_cns, chrname, inBamFile, technology, pg_runid_str;
		faidx_t *fai;
		vector<simpleReg_t*> sub_limit_reg_vec;
		bool limit_reg_process_flag, limit_reg_delete_flag;

		string refseqfilename, ctgfilename, readsfilename, alnfilename, clusterfilename;
		string rescue_refseqfilename, rescue_cnsfilename, rescue_readsfilename, rescue_alnfilename;
		int32_t ref_left_shift_size, ref_right_shift_size, ctg_num, min_sv_size, minReadsNumSupportSV, minClipEndSize, maxVarRegSize, minConReadLen, minMapQ: 10, minHighMapQ: 10, max_seg_num_per_read: 12, min_distance_merge;
		double max_ultra_high_cov, min_seqsim_match, min_seqsim_merge, max_seg_size_ratio, max_seg_nm_ratio, max_absig_density;
		vector<reg_t*> varVec, newVarVec;
		bool cns_success, align_success, call_success, clip_reg_flag, killed_flag;  	// default: false
		vector<minimap2_aln_t*> minimap2_aln_vec; // minimap2 aligned segments

		// clippings
		reg_t *clip_reg, *clip_reg_allele;
		int64_t leftClipRefPos, rightClipRefPos;
		int32_t sv_type, dup_num, supp_num_largeIndel, depth_largeIndel; // leftClipQueryPos, rightClipQueryPos, aln_orient;
		string bnd_mate_reg_strs[4];
		bool margin_adjusted_flag, large_indel_flag;

		// minimap2 aligned information
		vector<varCand*> *minimap2_aligned_info_vec;	// previously minimap2 aligned information

		// minimap2 aligned information file
		ofstream *minimap2_var_cand_file;

		// process monitor work
		int32_t max_proc_running_minutes;

		//genotyping parameter
		int32_t gt_min_sig_size, gt_min_sup_num_recover; // not used
		double gt_min_seqsim_merge, gt_homo_ratio_thres, gt_hete_ratio_thres;
		double gt_size_ratio_match; // not used

	public:
		varCand();
		virtual ~varCand();
		void callVariants();
		void setGtParas(int32_t gt_min_sig_size, double gt_size_ratio_match, double gt_min_seqsim_merge, double gt_homo_ratio_thres, double gt_hete_ratio_thres, int32_t min_sup_num_recover);
		vector<int32_t> computeDisagreeNumAndHighIndelBaseNum(string &chrname, size_t startRefPos, size_t endRefPos, string &inBamFile, faidx_t *fai);
		void adjustVarLocSlightly();
		void alnCtg2Refseq();
		void loadMinimap2AlnData();
		void destroyVarCand();

	private:
		void init();
		void destroyMinimap2AlnVec(vector<minimap2_aln_t*> &minimap2_aln_vec);
		void destoryClipAlnData(vector<clipAlnData_t*> &clipAlnDataVector);
		void recordMinimap2AlnInfo();
		bool getMinimap2AlnDoneFlag();
		void assignMinimap2AlnStatus();
		void callIndelVariants();
		void callClipRegVariants();
		void minimap2Parse();
		vector<minimap2_aln_t*> minimap2Parse(string &alnfilename, string &cnsfilename, string &refseqfilename);
		void getMissingPafInDelsAtAlnSegEnd(vector<struct pafalnSeg*> &all_pafalnsegs, vector<struct pafalnSeg*> &pafalnsegs, int32_t region_refstart, int32_t region_refend, int32_t querystart, int32_t queryend, int32_t query_len, int32_t subjectlen, int32_t subjectstart, int32_t subjectend, int32_t match_base_num, int32_t match_ref_len, string cons_seq, string ref_seq, int32_t aln_orient);
		vector<reg_t*> computeIndelVarLoc(vector<minimap2_aln_t*> &minimap2_aln_vec, string &ctgfilename_para, string &refseqfilename_para, bool rescue_flag);
		void recoverCandVarItems(vector<reg_t*> &var_vec, vector<reg_t*> &varVec, vector<int32_t> &cand_mm2_aln_idx_vec, vector<struct pafalnSeg*> &cand_paf_alnseg_vec, vector<int32_t> &recovered_qid_vec, int32_t mm2_aln_idx, vector<minimap2_aln_t*> &minimap2_aln_vec, string &refseq, int64_t start_var_pos, int64_t end_var_pos, int64_t start_var_pos_pure, int64_t end_var_pos_pure, vector<string> &qname_vec, vector<clipAlnData_t*> &clipAlnDataVector, string &ctgfilename_para, bool rescue_flag);
		void computeIntraMisPafAlnSeg(vector<reg_t*> &var_vec, vector<reg_t*> &varVec, vector<int32_t> &cand_mm2_aln_idx_vec, vector<struct pafalnSeg*> &cand_paf_alnseg_vec, vector<minimap2_aln_t*> &minimap2_aln_vec, string &ctgfilename_para, bool rescue_flag);
		void computeIndelFromSplitSegs(vector<reg_t*> &var_vec, vector<minimap2_aln_t*> &minimap2_aln_vec, string &ctgfilename_para, int64_t startRefPos_cns, int64_t endRefPos_cns, double size_ratio_match_thres, double size_ratio_match_thres_merge, bool rescue_flag);
		vector<reg_t*> computeIndelFromSplitSegsSingleQuery(vector<minimap2_aln_t*> &minimap2_aln_vec, string &ctgfilename_para, int32_t query_id_para, int64_t startRefPos_cns, int64_t endRefPos_cns, double size_ratio_match_thres, double size_ratio_match_thres_merge, bool rescue_flag);

		// merge neighbouring indels
		void mergeNeighbouringVars(vector<reg_t*> &regVector, int32_t min_ref_dist_arbitary_thres, int32_t max_ref_dist_thres, double min_merge_seqsim_thres, double min_valid_sig_size_ratio_thres, faidx_t *fai, string &contigfilename, string &reffilename);

		vector<reg_t*> rescueIndelVarLoc(vector<int32_t> &clusterId_incomplete);
		vector<int32_t> computeSuppNumFromRegionAlnSegs(vector<string> &clu_qname_vec, struct pafalnSeg* paf_alnseg, vector<clipAlnData_t*> &clipAlnDataVector, string &chrname, int64_t startRefPos_cns, int64_t endRefPos_cns, double size_ratio_match_thres, double size_ratio_match_thres_merge);
		vector<int32_t> computeSuppNumFromRegionAlnSegs(vector<string> &clu_qname_vec, int32_t opflag_para, int32_t oplen_para, string &chrname_para, int64_t start_var_pos_para, int64_t end_var_pos_para, int64_t startRefPos_cns, int64_t endRefPos_cns, double size_ratio_match_thres, double size_ratio_match_thres_merge);
		void refineAlleleFreq(vector<reg_t*> &var_vec);
		bool isAlleleFreqOk(vector<reg_t*> &var_vec);
		void refineAlleleFreqOp(vector<reg_t*> &var_vec);
		void computeGenotypeIndelReg(vector<reg_t*> &var_vec);

		// clippings
		vector<reg_t*> computeClipRegVarLoc(string &alnfilename, string &refseqfilename, string &cnsfilename, string &clusterfilename, vector<minimap2_aln_t*> &minimap2_aln_vec, vector<int32_t> *clusterId_incomplete, bool rescue_flag);
		vector<int32_t> getClusterIdCalledIncomplete(vector<reg_t*> &var_vec, string &clusterfilename, vector<minimap2_aln_t*> &minimap2_aln_vec);
		vector<int32_t> getSuppNumCovClipReg(string bnd_mate_reg_strs[4], int32_t supp_num_largeIndel, int32_t depth_largeIndel, vector<string> &target_qname_vec, bool large_indel_flag);
		void computeGenotypeClipReg(vector<reg_t*> &var_vec);
		vector<reg_t*> rescueDupInvClipReg(vector<int32_t> &clusterId_incomplete);
		vector<reg_t*> rescueLargeIndelClipReg(vector<int32_t> &clusterId_incomplete);
		vector<int32_t> getMinimapItemIdVec(vector<string> &queryname_vec, vector<minimap2_aln_t*> &minimap2_aln_vec);
		vector<int32_t> getLeftRightMinimapItemId(int64_t leftClipRefPos, int64_t rightClipRefPos, vector<int32_t> &minimap2_item_id_vec, vector<minimap2_aln_t*> &minimap2_aln_vec);
		vector<int32_t> getLeftRightClipAlnId(int64_t leftClipRefPos, int64_t rightClipRefPos, vector<clipAlnData_t*> &query_aln_segs);
		vector<int32_t> getLeftRightClipAlnIdINV(int64_t leftClipRefPos, int64_t rightClipRefPos, vector<clipAlnData_t*> &query_aln_segs);

		// call indel from reads
		vector<reg_t*> callIndelFromReadsIndelReg(vector<int32_t> &clusterId_incomplete);
		vector<reg_t*> callVarFromReadsClipReg(vector<int32_t> &clusterId_incomplete);
		void updateClusterIdIncomplete(vector<reg_t*> &var_vec, vector<reg_t*> &var_vec_rescue, vector<int32_t> &clusterId_incomplete);

		// output
		void printSV();
		void outputNewVarcandVec();
};

#endif /* SRC_VARCAND_H_ */
