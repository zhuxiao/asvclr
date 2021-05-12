#ifndef SRC_PARAS_H_
#define SRC_PARAS_H_

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <getopt.h>
#include <unistd.h>

#include "structures.h"

using namespace std;

// program variables
#define PROG_NAME					"ASVCLR"
#define PROG_DESC					"Accurate Structural Variant Caller for Long Reads"
#define PROG_VERSION				"0.10.2"
#define VCF_VERSION					"4.2"

#define CANU_RECOMMEND_VERSION		"2.1"

#define PACBIO_CLR_TECH_STR			"pacbio"
#define PACBIO_CCS_TECH_STR			"pacbio-hifi"
#define NANOPORE_TECH_STR			"nanopore"
#define SEQUENCING_TECH_DEFAULT		PACBIO_CLR_TECH_STR

#define SIZE_EST_OP					0
#define NUM_EST_OP					1
#define SNV_EST_OP					2

#define LARGE_INDEL_SIZE_FACTOR		3

#define MD_MISMATCH   				10
#define BASE_INDEL_CON_UNUSED		(-1)


// default parameter values
#define BLOCKSIZE  					1000000
#define SLIDESIZE  					500

#define MIN_INDEL_EVENT_SIZE		2
#define MIN_INDEL_EVENT_NUM			5
#define MIN_SV_SIZE_USR				2
#define MAX_SV_SIZE_USR				50000

//#define MIN_CLIP_READS_NUM_THRES	7
#define MIN_SUPPORT_READS_NUM		7

#define MAX_VAR_REG_SIZE			50000
#define ASM_CHUNK_SIZE_INDEL		10000

#define SIZE_PERCENTILE_EST			0.95
#define NUM_PERCENTILE_EST			0.99995
#define AUX_ARR_SIZE				1001

#define MAX_REG_SUM_SIZE_EST		500000

#define ASSEMBLY_SUCCESS			"ASM_SUCCESS"
#define ASSEMBLY_FAILURE			"ASM_FAILURE"
#define ALIGN_SUCCESS				"ALN_SUCCESS"
#define ALIGN_FAILURE				"ALN_FAILURE"
#define CALL_SUCCESS				"CALL_SUCCESS"
#define CALL_FAILURE				"CALL_FAILURE"

#define DONE_STR					"DONE"

#define SAMPLED_STR					"SAMPLED"
#define UNSAMPLED_STR				"UNSAMPLED"

#define LIMIT_REG_ALL_STR			"ALL"

#define EXPECTED_COV_ASSEMBLE		30.0f

#define NUM_PARTS_PROGRESS			50
#define NUM_THREADS_PER_ASSEM_WORK	0  // 0: unspecify the limited number of threads for each Canu work

#define OUT_DIR						"output"
#define SAMPLE_DEFAULT				"sample"

#define MAX_ULTRA_HIGH_COV_THRES	300		// maximal coverage threshold for ultra-high coverage

#define MIN_ADJACENT_REG_DIST		50

#define MIN_HIGH_CONSENSUS_INS_RATIO		0.3f
#define MIN_HIGH_CONSENSUS_DEL_RATIO		0.5f	// 0.4


// program parameters
class Paras
{
	public:
		// user/system defined parameters
		string command, refFile, inBamFile, outFilePrefix, sample, pg_cmd_str, canu_version, technology;
		string outDir;
		string out_dir_detect = "1_candidates";    // "1_candidates"
		string out_dir_assemble = "2_assemble";  // "2_assemble"
		string out_dir_call = "3_call";      // "3_call"
		string out_dir_tra = out_dir_call + "/" + "tra";
		string out_dir_result = "4_results";	// "4_results"
		int32_t blockSize, slideSize, min_sv_size_usr, max_sv_size_usr, num_threads, large_indel_size_thres; // , assemSlideSize;
		bool maskMisAlnRegFlag, load_from_file_flag;
		size_t misAlnRegLenSum = 0;
		int64_t minReadsNumSupportSV; //, minClipReadsNumSupportSV;

		// limit SV regions, item format: CHR | CHR:START-END
		vector<simpleReg_t*> limit_reg_vec;
		bool limit_reg_process_flag = false;	// true for limit regions; default is false for disable limit regions (i.e. process all regions)
		string limit_reg_filename = "limit_regions.bed";

		int32_t maxVarRegSize, minClipEndSize, assemChunkSize;

		int64_t mean_read_len, total_read_num_est;
		int32_t reg_sum_size_est, max_reg_sum_size_est;

		// reads sampling parameters
		double expected_cov_assemble;
		bool delete_reads_flag;

		// clipping reads sampling parameters
		double max_ultra_high_cov;

		// estimated parameters
		int32_t min_ins_size_filt, min_del_size_filt, min_clip_size_filt;
		int32_t min_ins_num_filt, min_del_num_filt, min_clip_num_filt;

		// auxiliary data for estimation
		int64_t insSizeEstArr[AUX_ARR_SIZE+1], delSizeEstArr[AUX_ARR_SIZE+1], clipSizeEstArr[AUX_ARR_SIZE+1];
		int64_t insNumEstArr[AUX_ARR_SIZE+1], delNumEstArr[AUX_ARR_SIZE+1], clipNumEstArr[AUX_ARR_SIZE+1];

		// assembly regions for thread pool
		vector<assembleWork_opt*> assem_work_vec;
		int32_t assemble_reg_preDone_num, assemble_reg_work_total, assemble_reg_workDone_num;
		int16_t num_parts_progress, num_threads_per_assem_work;
		pthread_mutex_t mtx_assemble_reg_workDone_num;

	public:
		Paras();
		Paras(int argc, char **argv);
		virtual ~Paras();
		void initEst();
		void estimate(size_t op_est);
		void outputParas();
		void outputEstParas(string info);

	private:
		void init();
		string getCanuVersion();
		bool isRecommendCanuVersion(string &canu_version, const string &recommend_version);
		int checkBamFile();
		int parseParas(int argc, char **argv);
		string getPgCmd(int argc, char **argv);
		void showVersion();
		int parseDetectParas(int argc, char **argv);
		int parseAssembleParas(int argc, char **argv);
		int parseCallParas(int argc, char **argv);
		int parseAllParas(int argc, char **argv);
		void showUsage();
		void showDetectUsage();
		void showAssembleUsage();
		void showCallUsage();
		void showAllUsage();
		int64_t estimateSinglePara(int64_t *arr, int32_t n, double threshold, int32_t min_val);
		int parse_long_opt(int32_t option_index, const char *optarg, const struct option *lopts);
};

#endif /* SRC_PARAS_H_ */
