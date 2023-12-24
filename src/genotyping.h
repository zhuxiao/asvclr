#ifndef SRC_GENOTYPING_H_
#define SRC_GENOTYPING_H_

#include <iostream>
#include <cstring>
#include <vector>
#include <algorithm>
#include <htslib/faidx.h>

#include "structures.h"
#include "alnDataLoader.h"
#include "util.h"

using namespace std;

#define GT_SIG_MATCH_SCORE					2
#define GT_SIG_MISMATCH_SCORE				-5
#define GT_SIG_GAP_SCORE					-2
#define GT_SIG_GAP_OPEN_SCORE				-4

#define GT_NOZYGOUS							0
#define GT_HOMOZYGOUS						1	// homozygous
#define GT_HETEROZYGOUS						2	// heterozygous

#define GT_NOZYGOUS_STR						"0/0"
#define GT_HOMOZYGOUS_STR					"1/1"	// homozygous
#define GT_HETEROZYGOUS_STR					"0/1"	// heterozygous

#define GT_HOMO_RATIO_THRES					(0.75f)
#define GT_HETER_RATIO_THRES				(0.2f)

#define GT_STR_DEFAULT						"GT:AD:DP\t./.:.,.:."


class genotyping{
	public:
		faidx_t *fai;
		reg_t *reg;
		string inBamFile;
		int32_t sig_size_thres, clip_size_thres, clip_extend_match_thres, min_sup_num_recover, noProfileQueryNum;
		double size_ratio_match_thres, valid_summed_size_ratio_read, min_alle_ratio_thres, max_alle_ratio_thres;
		queryGtSig_t *seed_gtQuery;
		vector<bool> validGtSigFlagVec;
		vector<string> gt_str_vec;

		vector<bam1_t*> alnDataVector;
		vector<queryGtSig_t*> queryGtSig_vec;
		vector<profile_pat_t*> match_profile_pat_vec;

	public:
		genotyping(reg_t *reg, faidx_t *fai, string &inBamFile, int32_t sig_size_thres, double size_ratio_match_thres, double min_alle_ratio_thres, double max_alle_ratio_thres, int32_t min_sup_num_recover);
		virtual ~genotyping();
		void destroyMatchProfilePatVec(vector<profile_pat_t*> &match_profile_pat_vec);
		void computeGenotype();
		void filterInvalidAlnData(vector<bam1_t*> &alnDataVector, double valid_summed_size_ratio_read);
		vector<bam1_t*> getQueryAlnSegs(vector<bam1_t*> &alnDataVector, string &queryname);
		int32_t getAlnSizeSingleQuery(vector<bam1_t*> &query_aln_segs);
		vector<int32_t> getAlnSizeSingleSeg(bam1_t *b);
		vector<queryGtSig_t*> extractGtSigVec();
		//vector<gtSig_t*> extractGtSigsSingleQuery(bam1_t *b);
		vector<gtSig_t*> extractGtSigsFromAlnSegsSingleQuery(vector<struct alnSeg*> &alnSegs, int64_t startPos, int64_t endPos, string &refseq, int32_t sig_size_thres, int32_t clip_size_thres);
		gtSig_t *allocateGtSigNode(struct alnSeg *aln_seg);
		void destroyQueryGtSigVec(vector<queryGtSig_t*> &gtSig_vec);
		void assignRegContainFlag(vector<queryGtSig_t*> &queryGtSig_vec, reg_t *reg);
		void gtSigMatch(vector<queryGtSig_t*> &gtSig_vec);
		queryGtSig_t *chooseSeedGtQuery(vector<queryGtSig_t*> &queryGtSig_vec);
		void computeGtMatchProfile(vector<queryGtSig_t*> &queryGtSig_vec, queryGtSig_t *seed_gtQuery);
		vector<int8_t> computeGtMatchProfileSingleQuery(queryGtSig_t *queryGtSig, queryGtSig_t *seed_gtQuery);
		bool isGtSigMatch(gtSig_t *gt_sig, gtSig_t *seed_gt_sig, double size_ratio_match_thres);
		vector<int8_t> computeSigMatchProfile(struct alnScoreNode *scoreArr, int32_t rowsNum, int32_t colsNum, queryGtSig_t *queryGtSig, queryGtSig_t *seed_gtQuery);
		vector<string> computeGenotypeString(vector<queryGtSig_t*> &queryGtSig_vec, vector<bool> &validGtSigFlagVec);
		vector<bool> computeValidGtSigFlags(queryGtSig_t *seed_gtQuery, reg_t *reg);
		void recoverVariants(reg_t *reg, vector<profile_pat_t*> &match_profile_pat_vec, queryGtSig_t *seed_gtQuery, vector<queryGtSig_t*> &queryGtSig_vec, vector<bool> &validGtSigFlagVec);
		int32_t getTargetPatSigIDFromProfile(vector<profile_pat_t*> &match_profile_pat_vec);
		int32_t computeSigColID(int32_t pat_sig_idx, vector<bool> &validGtSigFlagVec);
		vector<int64_t> computeAverVarReg(gtSig_t *target_sig, queryGtSig_t *seed_gtQuery, vector<queryGtSig_t*> &queryGtSig_vec);
		int64_t getSigEndPos(gtSig_t *gt_sig);
		int32_t getSigLen(gtSig_t *gt_sig);
		void compareAndUpdateVarReg(reg_t *reg, gtSig_t *target_sig, vector<int64_t> &aver_pos_vec);
};

#endif /* SRC_GENOTYPING_H_ */
