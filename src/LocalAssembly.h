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

using namespace std;

#define REFSEQ_SIDE_LEN					20000	// 10000
#define ASSEMBLY_SIDE_LEN				10000
#define ASSEMBLY_GENOME_SIZE_INITIAL	30000
#define ASSEMBLY_STEP_SIZE				10000
#define ASSEMBLY_SIDE_EXT_SIZE			5000

struct fqSeqNode{
	size_t seq_id;
	string seq_name, seq, qual;
	bool selected_flag;
};

class LocalAssembly {
	public:
		string chrname, readsfilename, contigfilename, refseqfilename, tmpdir, inBamFile; // canu_version;
		int64_t chrlen, assembly_extend_size, startRefPos_assembly, endRefPos_assembly;
		size_t num_threads_per_assem_work;
		bool assem_success_preDone_flag;

		// limit process regions
		bool limit_reg_process_flag;
		vector<simpleReg_t*> limit_reg_vec;

		vector<reg_t*> varVec;
		faidx_t *fai;

		// sampling
		size_t ref_seq_size, reads_count_original, total_bases_original, reads_count, total_bases;
		double local_cov_original, sampled_cov, expected_cov, compensation_coefficient, mean_read_len;
		bool sampling_flag, delete_reads_flag;

	private:
		vector<bam1_t*> alnDataVector;
		vector<clipAlnData_t*> clipAlnDataVector;

	public:
		LocalAssembly(string &readsfilename, string &contigfilename, string &refseqfilename, string &tmpdir, size_t num_threads_per_assem_work, vector<reg_t*> &varVec, string &chrname, string &inBamFile, faidx_t *fai, size_t assembly_extend_size, double expected_cov, bool delete_reads_flag);
		virtual ~LocalAssembly();
		void extractRefseq();
		void extractReadsDataFromBAM();
		bool localAssembleCanu();
		void recordAssemblyInfo(ofstream &assembly_info_file);
		void setLimitRegs(bool limit_reg_process_flag, vector<simpleReg_t*> limit_reg_vec);

	private:
		void destoryAlnData();
		void destoryClipAlnData();
		void destoryFqSeqs(vector<struct fqSeqNode*> &fq_seq_vec);
		double computeCompensationCoefficient(size_t startRefPos_assembly, size_t endRefPos_assembly, double mean_read_len);
		double computeLocalCov(vector<struct fqSeqNode*> &fq_seq_vec, double compensation_coefficient);
		void samplingReads(vector<struct fqSeqNode*> &fq_seq_vec, double expect_cov_val, double compensation_coefficient);
		void samplingReadsOp(vector<struct fqSeqNode*> &fq_seq_vec, double expect_cov_val);
		void saveSampledReads(string &readsfilename, vector<struct fqSeqNode*> &fq_seq_vec);
		vector<clipAlnData_t*> getQueryClipAlnSegs(string &queryname, vector<clipAlnData_t*> &clipAlnDataVector);
		int32_t getNoHardClipAlnItem(vector<clipAlnData_t*> &clipAlnDataVector);
		void markHardClipSegs(size_t idx, vector<clipAlnData_t*> &query_aln_segs);
		vector<string> getQuerySeqWithSoftClipSeqs(clipAlnData_t* clip_aln);
		vector<string> getQuerySeqWithoutSoftClipSeqs(clipAlnData_t* clip_aln);
		bool localAssembleCanu_IncreaseGenomeSize();
		bool localAssembleCanu_DecreaseGenomeSize();
		vector<string> joinQueryAlnSegs(vector<clipAlnData_t*> &query_aln_segs);
		int32_t getLeftMostAlnSeg(vector<clipAlnData_t*> &query_aln_segs);
};

#endif /* SRC_LOCALASSEMBLY_H_ */
