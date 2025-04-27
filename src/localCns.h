#ifndef SRC_LOCALCNS_H_
#define SRC_LOCALCNS_H_

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <limits.h>
#include <sys/stat.h>

#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/faidx.h>

#include "structures.h"
#include "RefSeqLoader.h"
#include "util.h"
#include "indelRegCluster.h"
#include "clipRegCluster.h"


using namespace std;

#define MIN_CANU_VERSION_NO_RAW			"2.0"
#define MIN_CANU_VERSION_HIFI			"1.9"
#define MIN_CANU_VERSION_NO_GNPPLOT		"1.8"

#define MAX_SV_LEN_USING_POA			2000
#define MIN_SEQ_LEN_USING_MINIMIZER		5000 // 2000

class localCns {
	public:
		string chrname, readsfilename, contigfilename, refseqfilename, clusterfilename, tmpdir, inBamFile, technology, canu_version;
		string readsfilename_prefix, readsfilename_suffix;
		vector<string> readsfilename_vec;
		int64_t chrlen, cns_extend_size, startRefPos_cns, endRefPos_cns;
		int32_t num_threads_per_cns_work, minClipEndSize, minConReadLen, min_sv_size, min_supp_num, minMapQ: 16, minHighMapQ: 16, sv_len_est;
		double max_ultra_high_cov, max_seg_size_ratio, min_identity_match, max_seg_nm_ratio;
		bool cns_success_preDone_flag, cns_success_flag, use_poa_flag, clip_reg_flag;
		double min_input_cov_canu;

		vector<struct seqsVec*> seqs_vec;

		// start time and end time
		time_t start_time, end_time;

		// limit process regions
		bool limit_reg_process_flag;
		vector<simpleReg_t*> limit_reg_vec;

		vector<reg_t*> varVec;
		faidx_t *fai;

		// sampling
		size_t ref_seq_size, reads_count_original, total_bases_original, reads_count, total_bases;
		double local_cov_original, sampled_cov, expected_cov, compensation_coefficient, mean_read_len;
		bool sampling_flag, delete_reads_flag, keep_failed_reads_flag;

	private:
		vector<bam1_t*> alnDataVector;
		vector<clipAlnData_t*> clipAlnDataVector;

	public:
		localCns(string &readsfilename, string &contigfilename, string &refseqfilename, string &clusterfilename, string &tmpdir, string &technology, double min_identity_match, int32_t sv_len_est, size_t num_threads_per_cns_work, vector<reg_t*> &varVec, string &chrname, string &inBamFile, faidx_t *fai, size_t cns_extend_size, double expected_cov, double min_input_cov, double max_ultra_high_cov, int32_t minMapQ, int32_t minHighMapQ, bool delete_reads_flag, bool keep_failed_reads_flag, bool clip_reg_flag, int32_t minClipEndSize, int32_t minConReadLen, int32_t min_sv_size, int32_t min_supp_num, double max_seg_size_ratio, double max_seg_nm_ratio);
		virtual ~localCns();
		void extractRefseq();
		void extractReadsDataFromBAM();
		bool localCnsCanu();
		bool localCnsWtdbg2();
		bool cnsByPoa();
		bool localConsensus();
		void recordCnsInfo(ofstream &cns_info_file);
		void setLimitRegs(bool limit_reg_process_flag, vector<simpleReg_t*> limit_reg_vec);

	private:
		void destoryClipAlnData();
		void destorySeqsVec(vector<struct seqsVec*> &seqs_vec);//
		void destoryQueryCluVec(vector<struct querySeqInfoVec*> &query_clu_vec);
		void destoryQueryCluVecClipReg(vector<qcSigListVec_t*> &query_clu_vec_clipReg);
		void saveClusterInfo(string &clusterfilename, vector<struct seqsVec*> &seqs_vec);
		bool updateUsepoaFlag(vector<clipAlnData_t*> &clipAlnDataVector, int32_t maxVarRegSize, double repeat_reads_thres);
		int32_t calculateSignalDistance(alnSeg *query_alnSeg1, alnSeg *query_alnSeg2);
		struct seqsVec* collectQuerySeqDataClipReg(vector<qcSigList_t*> &qcSigList);
		double computeCompensationCoefficient(size_t startRefPos_cns, size_t endRefPos_cns, double mean_read_len);
		double computeLocalCovClipReg(vector<struct fqSeqNode*> &fq_seq_vec, double compensation_coefficient);
		double computeLocalCovIndelReg(vector<struct querySeqInfoNode*> &query_seq_info_all, double compensation_coefficient);
		void samplingReadsClipReg(vector<struct fqSeqNode*> &fq_seq_vec, double expect_cov_val, double compensation_coefficient);
		void samplingReadsClipRegOp(vector<struct fqSeqNode*> &fq_seq_vec, double expect_cov_val);
		void samplingReads(vector<struct querySeqInfoNode*> &query_seq_info_all, double expect_cov_val, double compensation_coefficient);
		void samplingReadsOp(vector<struct querySeqInfoNode*> &query_seq_info_all, double expect_cov_val);
		void saveSampledReads(string &readsfilename, vector<struct fqSeqNode*> &fq_seq_vec);
		void markHardClipSegs(size_t idx, vector<clipAlnData_t*> &query_aln_segs);
		vector<string> getQuerySeqWithoutSoftClipSeqs(clipAlnData_t* clip_aln);
		bool localCnsCanu_IncreaseGenomeSize();
		bool localCnsCanu_DecreaseGenomeSize();
		vector<string> joinQueryAlnSegs(vector<clipAlnData_t*> &query_aln_segs);
		int32_t getLeftMostAlnSeg(vector<clipAlnData_t*> &query_aln_segs);
};

#endif /* SRC_LOCALCNS_H_ */
