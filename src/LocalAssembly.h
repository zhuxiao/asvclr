#ifndef SRC_LOCALASSEMBLY_H_
#define SRC_LOCALASSEMBLY_H_

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
#include "abpoa.h"

using namespace std;

#define MIN_CANU_VERSION_NO_RAW			"2.0"
#define MIN_CANU_VERSION_HIFI			"1.9"
#define MIN_CANU_VERSION_NO_GNPPLOT		"1.8"
#define QC_SIZE_RATIO_MATCH_THRES		 0.7


struct fqSeqNode{
	size_t seq_id;
	string seq_name, seq, qual;
	bool selected_flag;
};

struct querySeqInfoNode{
	//size_t seq_id;
	string qname, seq;
	bool cluster_finished_flag;
	vector<struct alnSeg*> query_alnSegs;
};

struct querySeqInfoVec{
	vector<struct querySeqInfoNode*> query_seq_info_all;
};

struct seqsVec{
	vector<string> qname;
	vector<string> seqs;
};

struct seedQueryInfo{
	int64_t startSpanPos, endSpanPos, span_intersection;
	int32_t id;
};

//for q_cluster init B
typedef struct{
	int32_t num;
	vector<int32_t> pos_id;
	double SR;
}srmin_match_t;

// single signature for cluster
typedef struct qcSigNode{
	int32_t sig_id: 24, cigar_op: 8, cigar_op_len;
	bool reg_contain_flag;
	int64_t ref_pos;
	string chrname;
	struct qcSigNode *mate_qcSig;
}qcSig_t;

// cluster signatures of a query
typedef struct queryCluSigNode{
	//int32_t group_id;
	//bool seed_flag;
	//string queryname;
	//size_t seq_id;
	vector<int8_t> match_profile_vec;
	vector<qcSig_t*> qcSig_vec;
}queryCluSig_t;

class LocalAssembly {
	public:
		string chrname, readsfilename, contigfilename, refseqfilename, tmpdir, inBamFile, technology, canu_version;
		int64_t chrlen, assembly_extend_size, startRefPos_assembly, endRefPos_assembly;
		int32_t num_threads_per_assem_work, minClipEndSize, minConReadLen, min_sv_size;
		double max_seg_size_ratio;
		bool assem_success_preDone_flag, assem_success_flag, use_poa_flag, clip_reg_flag;
		double min_input_cov_canu;

		//AaCcGgTtNn ==> 0,1,2,3,4
		unsigned char nt4_table[256] = {
		    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
		    16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 26, 26, 26, 26, 26,
		    26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26,
		    26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26,
		    26, 0, 5, 1, 6, 7, 8, 2, 9, 10, 11, 12, 13, 14, 4, 15,
		    16, 17, 18, 19, 3, 20, 21, 22, 23, 24, 25, 26, 26, 26, 26, 26,
		    26, 0, 5, 1, 6, 7, 8, 2, 9, 10, 11, 12, 13, 14, 4, 15,
		    16, 17, 18, 19, 3, 20, 21, 22, 23, 24, 25, 26, 26, 26, 26, 26,
		    26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26,
		    26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26,
		    26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26,
		    26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26,
		    26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26,
		    26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26,
		    26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26,
		    26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26};

		vector<struct seqsVec*> seqs_vec;
		//vector<string> seqs;
		//vector<struct alnSeg*> query_alnSegs;

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
		LocalAssembly(string &readsfilename, string &contigfilename, string &refseqfilename, string &tmpdir, string &technology, string &canu_version, size_t num_threads_per_assem_work, vector<reg_t*> &varVec, string &chrname, string &inBamFile, faidx_t *fai, size_t assembly_extend_size, double expected_cov, double min_input_cov, bool delete_reads_flag, bool keep_failed_reads_flag, bool clip_reg_flag, int32_t minClipEndSize, int32_t minConReadLen, int32_t min_sv_size, double max_seg_size_ratio);
		virtual ~LocalAssembly();
		void extractRefseq();
		void extractReadsDataFromBAM();
		bool localAssembleCanu();
		bool cnsByPoa();
		void recordAssemblyInfo(ofstream &assembly_info_file);
		void setLimitRegs(bool limit_reg_process_flag, vector<simpleReg_t*> limit_reg_vec);

	private:
		void destoryAlnData();
		void destoryClipAlnData();
		void destoryFqSeqs(vector<struct fqSeqNode*> &fq_seq_vec);
		void destorySeqsVec(vector<struct seqsVec*> &seqs_vec);//
		void destoryQueryCluVec(vector<struct querySeqInfoVec*> &query_clu_vec);
		void destoryQuerySeqInfoAll(vector<struct querySeqInfoNode*> &query_seq_info_all);
		void destroyQueryQcSig(queryCluSig_t *qc_Sig);
		vector<struct querySeqInfoVec*> queryCluster(vector<struct querySeqInfoNode*> &query_seq_info_all);
		void resucueCluster(vector<struct querySeqInfoNode*> &query_seq_info_all, vector<struct querySeqInfoNode*> &q_cluster_a, vector<struct querySeqInfoNode*> &q_cluster_b);
		double computeScoreRatio(struct querySeqInfoNode *query_seq_info_node, struct querySeqInfoNode *q_cluster_node, int64_t startSpanPos, int64_t endSpanPos);
		vector<int8_t> computeQcMatchProfileSingleQuery(queryCluSig_t *queryCluSig, queryCluSig_t *seed_qcQuery);
		vector<int8_t> qComputeSigMatchProfile(struct alnScoreNode *scoreArr, int32_t rowsNum, int32_t colsNum, queryCluSig_t *queryCluSig, queryCluSig_t *seed_qcQuery);
		bool isQcSigMatch(qcSig_t *qc_sig, qcSig_t *seed_qc_sig, double size_ratio_match_thres);
		struct seedQueryInfo* chooseSeedClusterQuery(struct querySeqInfoNode* query_seq_info_node, vector<struct querySeqInfoNode*> &q_cluster);
		struct seedQueryInfo* chooseSeedClusterQuery02(struct seedQueryInfo* seedinfo, vector<struct querySeqInfoNode*> &q_cluster);
		vector<qcSig_t*> extractQcSigsFromAlnSegsSingleQuery(struct querySeqInfoNode *query_seq_info_node, int64_t startSpanPos, int64_t endSpanPos);
		qcSig_t *allocateQcSigNode(struct alnSeg *aln_seg);
		//void destroyQueryCluSigVec(vector<queryCluSig_t*> &queryCluSig_vec);
		struct seqsVec *smoothQuerySeqData(string &refseq, vector<struct querySeqInfoNode*> &query_seq_info_vec);
		double computeCompensationCoefficient(size_t startRefPos_assembly, size_t endRefPos_assembly, double mean_read_len);
		double computeLocalCov(vector<struct fqSeqNode*> &fq_seq_vec, double compensation_coefficient);
		void samplingReads(vector<struct fqSeqNode*> &fq_seq_vec, double expect_cov_val, double compensation_coefficient);
		void samplingReadsOp(vector<struct fqSeqNode*> &fq_seq_vec, double expect_cov_val);
		void saveSampledReads(string &readsfilename, vector<struct fqSeqNode*> &fq_seq_vec);
		int32_t getNoHardClipAlnItem(vector<clipAlnData_t*> &clipAlnDataVector);
		void markHardClipSegs(size_t idx, vector<clipAlnData_t*> &query_aln_segs);
		vector<string> getQuerySeqWithSoftClipSeqs(clipAlnData_t* clip_aln, string &refseq, int32_t bam_type);
		vector<string> getQuerySeqWithoutSoftClipSeqs(clipAlnData_t* clip_aln);
		vector<int32_t> computeQueryStartEndLocByRefPos(bam1_t *b, int64_t startRefPos, int64_t endRefPos, string &refseq, int32_t bam_type);
		void removeHeadTailAlnSegs(vector<struct alnSeg*> &alnSegs);
		bool localAssembleCanu_IncreaseGenomeSize();
		bool localAssembleCanu_DecreaseGenomeSize();
		vector<string> joinQueryAlnSegs(vector<clipAlnData_t*> &query_aln_segs);
		int32_t getLeftMostAlnSeg(vector<clipAlnData_t*> &query_aln_segs);
};

#endif /* SRC_LOCALASSEMBLY_H_ */
