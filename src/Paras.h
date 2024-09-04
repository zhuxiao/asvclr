#ifndef SRC_PARAS_H_
#define SRC_PARAS_H_

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <getopt.h>
#include <unistd.h>

#include "structures.h"
#include "meminfo.h"


using namespace std;

// program variables
#define PROG_NAME					"ASVCLR"
#define PROG_DESC					"Allele-aware Structural Variant Caller for Long Reads"
#define PROG_VERSION				"1.4.4"
#define VCF_VERSION					"4.2"

#define CMD_DET_STR					"det"
#define CMD_CNS_STR					"cns"
#define CMD_CALL_STR				"call"
#define CMD_ALL_STR					"all"
#define CMD_DET_CNS_STR				"det-cns"

#define CANU_RECOMMEND_VERSION		"2.1"

//#define PACBIO_CLR_TECH_STR			"pacbio"
//#define PACBIO_CCS_TECH_STR			"pacbio-hifi"
//#define NANOPORE_TECH_STR			"nanopore"
#define PACBIO_CCS_TECH_STR			"ccs"
#define PACBIO_RS_TECH_STR			"rs"
#define PACBIO_SQ_TECH_STR			"sq"
#define NANOPORE_TECH_STR			"ont"
#define SEQUENCING_TECH_DEFAULT		PACBIO_CCS_TECH_STR

#define DECOY_PREFIX				"hs37d"
#define DECOY_PREFIX2				"hs38d"

#define CHR_X_STR1					"chrX"
#define CHR_X_STR2					"X"
#define CHR_Y_STR1					"chrY"
#define CHR_Y_STR2					"Y"

#define VAR_DISCOV_L_UNUSED					0
#define VAR_DISCOV_L_CNS_ALN				1
#define VAR_DISCOV_L_RESCUE_CNS_ALN			2
#define VAR_DISCOV_L_READS					3

#define VAR_DISCOV_L_TITLE_STR				"LDISCOV"
#define VAR_DISCOV_L_UNUSED_STR				"UNC"
#define VAR_DISCOV_L_CNS_ALN_STR			"ALN"
#define VAR_DISCOV_L_RESCUE_CNS_ALN_STR		"RESCUE_ALN"
#define VAR_DISCOV_L_READS_STR				"READS"

#define SIZE_EST_OP					0
#define NUM_EST_OP					1
#define SNV_EST_OP					2

#define LARGE_INDEL_SIZE_FACTOR		3

#define MD_MISMATCH   				10
#define BASE_INDEL_CON_UNUSED		(-1)


// default parameter values
#define BLOCKSIZE  					1000000
#define SLIDESIZE  					500//500

#define MIN_INDEL_EVENT_SIZE		2
#define MIN_INDEL_EVENT_NUM			5
#define MIN_SV_SIZE_USR				20
#define MAX_SEG_SIZE_RATIO			0.5f
#define MAX_SV_SIZE_USR				50000

#define MIN_DUP_SIZE				30  // 2021-07-27

//#define MIN_CLIP_READS_NUM_THRES	7
#define MIN_SUPPORT_READS_NUM_EST		-1
#define MIN_SUPPORT_READS_NUM_FACTOR	0.03f	// 0.07f (updated 2024-06-25), 0.05f (updated 2024-09-04)
#define READS_NUM_SUPPORT_FACTOR		0.5f
#define MIN_SUPPORT_READS_NUM_DEFAULT	3

#define MAX_VAR_REG_SIZE			50000
#define CNS_CHUNK_SIZE_INDEL		1000	// 10000
#define CNS_CHUNK_SIZE_EXT_INDEL	1000	//1000
#define CNS_EXT_INDEL_FACTOR_500BP	2	// side extend size factor for mid chunk (> 500bp)
#define CNS_EXT_INDEL_FACTOR_1K		5	// side extend size factor for large chunk (> 1kb)
#define CNS_CHUNK_SIZE_EXT_CLIP		5000	//1000, 10000, 20000
#define CNS_EXT_CLIPREG_FACTOR_1K	(1.5f)	// side extend size factor for mid chunk (> 1kb)
#define CNS_EXT_CLIPREG_FACTOR_2K	2	// side extend size factor for mid chunk (> 2kb)
#define CNS_EXT_CLIPREG_FACTOR_4K	3	// side extend size factor for large chunk (> 4kb)
#define CNS_EXT_CLIPREG_FACTOR_6K	4	// side extend size factor for large chunk (> 6kb)
#define MIN_CONS_READ_LEN			100

#define MAX_REF_DIST_IDENTITY		2000

#define SIZE_PERCENTILE_EST			0.95	//0.95 (d-2024-03-23)
#define NUM_PERCENTILE_EST			0.99995
#define AUX_ARR_SIZE				1001

#define MAX_REG_SUM_SIZE_EST		500000

#define CNS_SUCCESS					"CNS_SUCCESS"
#define CNS_FAILURE					"CNS_FAILURE"
#define ALIGN_SUCCESS				"ALN_SUCCESS"
#define ALIGN_FAILURE				"ALN_FAILURE"
#define CALL_SUCCESS				"CALL_SUCCESS"
#define CALL_FAILURE				"CALL_FAILURE"

#define DONE_STR					"DONE"

#define SAMPLED_STR					"SAMPLED"
#define UNSAMPLED_STR				"UNSAMPLED"

#define LIMIT_REG_ALL_STR			"ALL"

#define MIN_INPUT_COV_CANU			5  		// 2021-08-01
#define EXPECTED_COV_CNS			40.0f	// 30.0f

#define NUM_PARTS_PROGRESS			100
#define NUM_THREADS_PER_CNS_WORK	0  		// 0: unspecify the limited number of threads for each CNS work

#define OUT_DIR						"output"
#define SAMPLE_DEFAULT				"sample"
#define RESULT_PREFIX_DEFAULT		"asvclr"

#define MAX_ULTRA_HIGH_COV_THRES	300		// maximal coverage threshold for ultra-high coverage, 100
#define MIN_MAPQ_THRES				0		// 10
#define MIN_HIGH_MAPQ_THRES			10

#define MIN_ADJACENT_REG_DIST		20		// 50

#define MIN_HIGH_CONSENSUS_INS_RATIO		0.3f
#define MIN_HIGH_CONSENSUS_DEL_RATIO		0.5f	// 0.4

#define MAX_CNS_MINUTES						15
#define MAX_ALN_MINUTES						15

//#define MAX_PROC_RUNNING_MINUTES				120
#define MAX_PROC_RUNNING_MINUTES_CNS			30	//300
#define MAX_PROC_RUNNING_MINUTES_CALL			30  //120
#define MONITOR_WAIT_SECONDS					60
#define ULTRA_LOW_PROC_RUNNING_MINUTES			30

//#define DEFAULT_MONITOR_PROC_NAMES				"overlapInCore,falconsense,blat"
#define DEFAULT_MONITOR_PROC_NAMES_CNS			"overlapInCore,falconsense"
#define DEFAULT_MONITOR_PROC_NAMES_CALL			"blat"

// consensus parameters
#define CNS_GENOME_SIZE_INITIAL				30000
#define CNS_STEP_SIZE						10000
#define CNS_SIDE_EXT_SIZE					5000

#define MAX_RESCUE_VAR_SIZE					40000

// program parameters
class Paras
{
	public:
		// user/system defined parameters
		string command, refFile, inBamFile, outFilePrefix, sample, pg_cmd_str, technology;
		string wtdbg2_version, minimap2_version, abpoa_version, canu_version;
		string outDir;
		string out_dir_detect = "1_candidates";    // "1_candidates"
		string out_dir_cns = "2_cns";		// "2_cns"
		string out_dir_call = "3_call";      // "3_call"
		string out_dir_tra = out_dir_call + "/" + "tra";
		string out_dir_result = "4_results";	// "4_results"
		int32_t blockSize, slideSize, min_sv_size_usr, max_sv_size_usr, num_threads, large_indel_size_thres;
		double max_seg_size_ratio_usr, min_identity_match;
		bool maskMisAlnRegFlag, load_from_file_flag, include_decoy, include_alt;
		size_t misAlnRegLenSum = 0;
		int32_t minReadsNumSupportSV: 29, min_Nsupp_est_flag: 3; //, minClipReadsNumSupportSV; Nsupp_est_flag: 1 for estimated, 0 for user-specified
		int32_t minMapQ: 10, minHighMapQ: 10, max_seg_num_per_read: 12;

		// process monitor
		string monitoring_proc_names_cns, monitoring_proc_names_call;
		int32_t max_proc_running_minutes_cns, max_proc_running_minutes_call;

		// limit SV regions, item format: CHR | CHR:START-END
		vector<simpleReg_t*> limit_reg_vec;
		bool limit_reg_process_flag = false;	// true for limit regions; default is false for disable limit regions (i.e. process all regions)
		string limit_reg_filename = "limit_regions.bed";

		int32_t maxVarRegSize, minClipEndSize, cnsChunkSize, cnsSideExtSize, cnsSideExtSizeClip, minConReadLen;

		int64_t mean_read_len, total_read_num_est;
		//int64_t chr_mean_depth, total_depth, chrome_num;
		vector<int64_t> mean_depth_vec;
		int32_t reg_sum_size_est, max_reg_sum_size_est;

		// reads sampling parameters
		double expected_cov_cns;
		bool delete_reads_flag, keep_failed_reads_flag, recns_failed_work_flag;

		double min_input_cov_canu;
		//int32_t min_overlap_length_canu;

		// clipping reads sampling parameters
		double max_ultra_high_cov;

		// estimated parameters
		int32_t min_ins_size_filt, min_del_size_filt, min_clip_size_filt;
		int32_t min_ins_num_filt, min_del_num_filt, min_clip_num_filt;

		// auxiliary data for estimation
		int64_t insSizeEstArr[AUX_ARR_SIZE+1], delSizeEstArr[AUX_ARR_SIZE+1], clipSizeEstArr[AUX_ARR_SIZE+1];
		int64_t insNumEstArr[AUX_ARR_SIZE+1], delNumEstArr[AUX_ARR_SIZE+1], clipNumEstArr[AUX_ARR_SIZE+1];

		// consensus regions for thread pool
		vector<cnsWork_opt*> cns_work_vec;
		int32_t cns_reg_preDone_num, cns_reg_work_total, cns_reg_workDone_num;
		int16_t num_parts_progress, num_threads_per_cns_work;
		pthread_mutex_t mtx_cns_reg_workDone_num;

		// call works for thread pool
		vector<varCand*> call_work_vec;
		int32_t call_work_num, call_workDone_num;
		pthread_mutex_t mtx_call_workDone_num;

		// process monitor killed blat work
		vector<killedBlatWork_t*> killed_blat_work_vec;
		string killed_blat_work_filename = out_dir_call + "/" + "monitor_killed_blat_work";	// colums: alnfilename,ctgfilename,refseqfilename
		ofstream killed_blat_work_file;
		pthread_mutex_t mtx_killed_blat_work;

		// process monitor killed blat work
		vector<killedMinimap2Work_t*> killed_minimap2_work_vec;
		string killed_minimap2_work_filename = out_dir_call + "/" + "monitor_killed_minimap2_work";	// colums: alnfilename,ctgfilename,refseqfilename
		ofstream killed_minimap2_work_file;
		pthread_mutex_t mtx_killed_minimap2_work;

		//genotyping parameters
		int32_t gt_min_sig_size; // not used
		double gt_min_identity_merge, gt_homo_ratio, gt_hete_ratio;
		double gt_size_ratio_match; // not used

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
		string getProgramVersion(const string &cmd_str);
		bool isRecommendCanuVersion(string &canu_version, const string &recommend_version);
//		string getMinimap2Version();
//		string getAbpoaVersion();
		void initPreset();
		int checkBamFile();
		int parseParas(int argc, char **argv);
		string getPgCmd(int argc, char **argv);
		void show_version();
		int parseDetectParas(int argc, char **argv);
		int parseCnsParas(int argc, char **argv);
		int parseCallParas(int argc, char **argv);
		int parseAllParas(int argc, char **argv, const string &cmd_str);
		int parseDetectCnsParas(int argc, char **argv);
		void showUsage();
		void showDetectUsage();
		void showCnsUsage();
		void showCallUsage();
		void showAllUsage(const string &cmd_str);
		void showDetectCnsUsage();
		int64_t estimateSinglePara(int64_t *arr, int32_t n, double threshold, int32_t min_val);
		int64_t estimateMinReadsNumSupportSV(vector<int64_t> &mean_depth_vec);
		int parse_long_opt(int32_t option_index, const char *optarg, const struct option *lopts);
};

#endif /* SRC_PARAS_H_ */
