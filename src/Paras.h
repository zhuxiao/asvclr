#ifndef SRC_PARAS_H_
#define SRC_PARAS_H_

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <getopt.h>
#include <unistd.h>

#include <htslib/sam.h>

using namespace std;

// program variables
#define PROG_NAME		"ASVCLR"
#define PROG_DESC		"Accurate Structural Variation Caller for Long Reads"
#define PROG_VERSION	"0.5.0"

#define SIZE_EST_OP		0
#define NUM_EST_OP		1
#define SNV_EST_OP		2

#define LARGE_INDEL_SIZE_FACTOR		3


// default parameter values
#define BLOCKSIZE  						1000000
#define SLIDESIZE  						500
#define ASSEM_SLIDE_SIZE				5000

#define MIN_INDEL_EVENT_SIZE			2
#define MIN_INDEL_EVENT_NUM				5
#define MIN_SV_SIZE_USR					2

#define MIN_CLIP_READS_NUM_THRES		7	// to be parameterized

#define MAX_CLIP_REG_SIZE				10000  // to be parameterized

#define SIZE_PERCENTILE_EST				0.95
#define NUM_PERCENTILE_EST				0.99995
#define AUX_ARR_SIZE					1001

#define MAX_REG_SUM_SIZE_EST			500000

#define ASSEMBLY_SUCCESS			"ASS_SUCCESS"
#define ASSEMBLY_FAILURE			"ASS_FAILURE"
#define ALIGN_SUCCESS				"ALN_SUCCESS"
#define ALIGN_FAILURE				"ALN_FAILURE"
#define CALL_SUCCESS				"CALL_SUCCESS"
#define CALL_FAILURE				"CALL_FAILURE"

#define DONE_STR					"DONE"

#define SAMPLED_STR					"SAMPLED"
#define UNSAMPLED_STR				"UNSAMPLED"

#define EXPECTED_COV_ASSEMBLE		40.0f


// program parameters
class Paras
{
	public:
		// user/system defined parameters
		string command, refFile, inBamFile, outFilePrefix; //, canu_version;
		size_t blockSize, slideSize, assemSlideSize, min_sv_size_usr, num_threads, large_indel_size_thres;
		bool maskMisAlnRegFlag, load_from_file_flag;
		size_t misAlnRegLenSum = 0;
		size_t minClipReadsNumSupportSV;

		size_t maxClipRegSize;

		size_t reg_sum_size_est, max_reg_sum_size_est;

		// reads sampling parameters
		double expected_cov_assemble;
		bool delete_reads_flag;

		// estimated parameters
		size_t min_ins_size_filt, min_del_size_filt, min_clip_size_filt;
		size_t min_ins_num_filt, min_del_num_filt, min_clip_num_filt;

		// auxiliary data for estimation
		size_t insSizeEstArr[AUX_ARR_SIZE], delSizeEstArr[AUX_ARR_SIZE], clipSizeEstArr[AUX_ARR_SIZE];
		size_t insNumEstArr[AUX_ARR_SIZE], delNumEstArr[AUX_ARR_SIZE], clipNumEstArr[AUX_ARR_SIZE];

	public:
		Paras();
		Paras(int argc, char **argv);
		void initEst();
		void estimate(size_t op_est);
		void outputParas();
		void outputEstParas(string info);

	private:
		void init();
		//string getCanuVersion();
		int checkBamFile();
		int parseParas(int argc, char **argv);
		int parseDetectParas(int argc, char **argv);
		int parseAssembleParas(int argc, char **argv);
		int parseCallParas(int argc, char **argv);
		int parseAllParas(int argc, char **argv);
		void showUsage();
		void showDetectUsage();
		void showAssembleUsage();
		void showCallUsage();
		void showAllUsage();
		size_t estimateSinglePara(size_t *arr, size_t n, double threshold, size_t min_val);
};

#endif /* SRC_PARAS_H_ */
