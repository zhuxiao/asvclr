#ifndef SRC_CHROME_H_
#define SRC_CHROME_H_

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <iomanip>
#include <htslib/faidx.h>

#include "structures.h"
#include "Paras.h"
#include "Block.h"
#include "varCand.h"
#include "clipReg.h"
#include "Phasing.h"

using namespace std;

#define SNV_DETECT_FLAG		0

#define BLOCK_SIZE_EST		20000
#define DIST_CHR_END		20000
#define MIN_CHR_SIZE_EST	((DIST_CHR_END*2+BLOCK_SIZE_EST)*2)

#define MIN_SIZE_RATIO_REDUNDANT_FILTER		(0.98f)
#define MIN_SEQSIM_REDUNDANT_FILTER			(0.92f)		// (0.92f), updated on 2025-06-02, 0.93 deleted on 2026-05-09

#define REFSEQ_PATTERN		"refseq"
#define CLIPREG_PATTERN		"clipReg_refseq"

#define ULTRA_LARGE_STATUS	"ULarge"
#define NORMAL_STATUS       "Normal"

class Chrome{
	public:
		Paras *paras;
		vector<Chrome*> *chr_vec;
		string chrname;
		int64_t chrlen;
		bool print_flag, decoy_flag, alt_flag, valid_flag;

		// output directory
		string out_dir_detect, out_dir_cns, out_dir_call;
		string out_filename_detect_snv, out_filename_detect_indel, out_filename_detect_clipReg;
		string out_filename_call_snv, out_filename_call_indel, out_filename_call_clipReg;

		vector<mateClipReg_t*> mateClipRegVector;
		vector<varCand*> var_cand_vec;
		vector<varCand*> var_cand_clipReg_vec;

		vector<Block*> blockVector;

		int32_t blockNum, process_block_num;
		string blocks_out_file;
		faidx_t *fai;

		// clip regions
		vector<reg_t*> clipRegVector;  // only for duplication and inversion

		// consensus info
		string var_cand_indel_filename, var_cand_clipReg_filename;
		ofstream var_cand_indel_file, var_cand_clipReg_file;

	public:
		Chrome(string& chrname, int chrlen, faidx_t *fai, Paras *paras, vector<Chrome*> *chr_vec);
		virtual ~Chrome();
		void setOutputDir(string& out_dir_detect_prefix, string& out_dir_cns_prefix, string& out_dir_call_prefix);
		string getVarcandIndelFilename();
		string getVarcandClipregFilename();
		int generateChrBlocks();
		void saveChrBlocksToFile();
		void chrFillDataEst(size_t op_est);
		int chrDetect();
		void removeFPIndelSnvInClipReg(vector<mateClipReg_t*> &mate_clipReg_vec);
		void chrMergeDetectResultToFile();
		void loadPrevConsInfo();
		void chrLoadDataCons();
		int chrGenerateLocalConsWorkOpt();
		void chrResetConsData();
		Block *computeBlocByPos(int64_t begPos, vector<Block*> &block_vec);
		void chrLoadDataCall();
		void chrCollectCallWork();
		void removeVarCandNodeClipReg(varCand *var_cand);
		void chrPhasing();
		void saveCallSV2File();
		void chrLoadMateClipRegData();


	private:
		void init();
		Block *allocateBlock(string& chrname, int64_t begPos, int64_t endPos, faidx_t *fai, bool headIgnFlag, bool tailIgnFlag, bool block_process_flag);
		void destroyBlockVector();
		void destroyVarCandVector(vector<varCand*> &var_cand_vec);
		void destroyClipRegVector();
		void destroyMateClipRegVector();
		void chrLoadIndelDataCons();
		void chrLoadIndelData(bool limit_reg_process_flag, vector<simpleReg_t*> &limit_reg_vec);
		void chrLoadClipRegDataCons();
		int chrDetect_st();
		int chrDetect_mt();
		void removeRedundantIndelDetect();
		void removeRedundantIndelItemDetect(reg_t *reg, int32_t bloc_idx, int32_t indel_vec_idx);
		void chrResetVarCandFiles();
		int32_t computeBlocID(int64_t begPos, vector<Block*> &block_vec);
		int chrGenerateLocalConsWorkOpt_st();
		int chrGenerateLocalConsWorkOpt_mt();
		void loadPrevConsInfo2(bool clipReg_flag, bool limit_reg_process_flag, vector<simpleReg_t*> &limit_reg_vec, vector<cnsWork_opt*> &cns_work_vec);
		string getCnsFileHeaderLine();

		void loadVarCandData();
		void loadVarCandDataFromFile(vector<varCand*> &var_cand_vec, string &var_cand_filename, bool clipReg_flag, bool limit_reg_process_flag, vector<simpleReg_t*> &limit_reg_vec);
		void loadClipRegCandData();
		void sortVarCandData(vector<varCand*> &var_cand_vec);
		bool isVarCandDataSorted(vector<varCand*> &var_cand_vec);

		// DUP, INV, TRA
		void chrComputeMateClipReg();
		void processMateClipRegDetectWork();
		void sortMateClipRegsByWorkId();
		void bubbleSortMateClipRegs();
		void checkSortMateClipRegs();
		//void removeRedundantItemsClipReg(vector<reg_t*> &clipReg_vec, vector<bool> &clip_processed_flag_vec);
		void removeRedundantItemsClipReg(vector<reg_t*> &clipReg_vec);
		void removeFPClipRegsDupInv();
		void chrLoadMateClipRegDataOp(bool limit_reg_process_flag, vector<simpleReg_t*> &limit_reg_vec);
		string getDirnameCall(string &chrname_given);
		mateClipReg_t* getMateClipReg(reg_t *reg1, reg_t *reg2, int32_t *clip_reg_idx_tra, string &chrname_mate_clip_reg);
		mateClipReg_t* getMateClipRegOp(reg_t *reg1, reg_t *reg2, int32_t *clip_reg_idx_tra, vector<mateClipReg_t*> &mate_clipReg_vec);

		// save SV to file
		void saveCallIndelClipReg2File(string &outfilename_indel, string &outfilename_clipReg);
};

#endif /* SRC_CHROME_H_ */
