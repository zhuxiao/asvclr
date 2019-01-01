#ifndef SRC_CHROME_H_
#define SRC_CHROME_H_

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <htslib/faidx.h>

#include "Paras.h"
#include "Block.h"
#include "varCand.h"
#include "clipReg.h"
#include "util.h"

using namespace std;

#define BLOCK_SIZE_EST		20000
#define DIST_CHR_END		20000
#define MIN_CHR_SIZE_EST	(DIST_CHR_END*2+BLOCK_SIZE_EST)


class Chrome{
	public:
		// output directory
		string out_dir_detect, out_dir_assemble, out_dir_call;
		string out_filename_detect_snv, out_filename_detect_indel, out_filename_detect_clipReg;
		string out_filename_call_snv, out_filename_call_indel, out_filename_call_clipReg;

		vector<mateClipReg_t*> mateClipRegVector;

	private:
		Paras *paras;
		string chrname;
		int32_t chrlen, blockNum;
		string blocks_out_file;
		vector<Block*> blockVector;
		faidx_t *fai;

		// clip regions
		vector<reg_t*> clipRegVector;  // only for duplication and inversion

		// call
		vector<varCand*> var_cand_vec;
		vector<reg_t*> mis_aln_vec;
		vector<varCand*> var_cand_clipReg_vec;

		// assembly info and misAln region
		string var_cand_indel_filename, misAln_reg_filename, var_cand_clipReg_filename;
		ofstream var_cand_indel_file, misAln_reg_file, var_cand_clipReg_file;

	public:
		Chrome(string& chrname, int chrlen, faidx_t *fai, Paras *paras);
		virtual ~Chrome();
		void setOutputDir(string& out_dir_detect_prefix, string& out_dir_assemble_prefix, string& out_dir_call_prefix);
		int generateChrBlocks();
		void saveChrBlocksToFile();
		void chrFillDataEst(size_t op_est);
		int chrDetect();
		void chrMergeDetectResultToFile();
		void chrLoadDataAssemble();
		int chrLocalAssemble();
		void chrLoadDataCall();
		int chrCall();
		void chrFillVarseq();
		void saveCallSV2File();

	private:
		void init();
		Block *allocateBlock(string& chrname, size_t chrlen, int begPos, int endPos, faidx_t *fai, bool headIgnFlag, bool tailIgnFlag);
		void destroyBlockVector();
		void destroyVarCandVector(vector<varCand*> &var_cand_vec);
		void destroyMisAlnVector();
		void destroyClipRegVector();
		void destroyMateClipRegVector();
		void chrLoadIndelDataAssemble();
		void chrLoadClipRegDataAssemble();
		int chrDetect_st();
		int chrDetect_mt();
		void chrSetVarCandFiles();
		void chrResetVarCandFiles();
		void chrSetMisAlnRegFile();
		void chrResetMisAlnRegFile();
		int32_t computeBlocID(size_t begPos);
		int chrLocalAssemble_st();
		int chrLocalAssemble_mt();
		void outputAssemDataToFile(string &filename);

		void chrCall_st();
		void chrCall_mt();
		void chrCallVariants(vector<varCand*> &var_cand_vec);
		void loadVarCandData();
		void loadVarCandDataFromFile(vector<varCand*> &var_cand_vec, string &var_cand_filename, bool clipReg_flag);
		void loadClipRegCandData();
		void sortVarCandData(vector<varCand*> &var_cand_vec);
		bool isVarCandDataSorted(vector<varCand*> &var_cand_vec);
		void loadMisAlnRegData();
		void sortMisAlnRegData();
		void chrFillVarseqSingleVec(vector<varCand*> &var_cand_vec);
		void removeRedundantVar();
		void removeRedundantIndel(vector<varCand*> &var_cand_vec);
		//void removeRedundantClipReg(vector<varCand*> &var_cand_clipReg_vec);
		bool isRedundantVarItemBinSearch(reg_t *reg, vector<varCand*> &var_cand_vec);
		void removeRedundantVarItemsInNewCalledVarvec(reg_t *reg_idx, int32_t idx, vector<varCand*> &var_cand_vec);
		void removeFPNewVarVec();
		void removeFPNewVarVecIndel(vector<varCand*> &var_cand_vec);
		void removeFPIndelSnvInClipReg();
		bool isIndelInClipReg(reg_t *reg, vector<mateClipReg_t*> &mate_clipReg_vec);
		bool isSnvInClipReg(size_t pos, vector<mateClipReg_t*> &mate_clipReg_vec);

		// duplications and inversions
		void computeMateClipRegDupInv();
		void precessClipRegs(size_t idx, vector<bool> &clip_processed_flag_vec, mateClipReg_t &mate_clip_reg);
		void removeFPClipRegs();
		vector<size_t> getOverlapClipReg(reg_t *given_reg);
		void chrLoadMateClipRegData();
		mateClipReg_t* getMateClipReg(reg_t *reg1, reg_t *reg2, int32_t *clip_reg_idx_tra);

		// save SV to file
		void saveCallIndelClipReg2File(string &outfilename_indel, string &outfilename_clipReg);
		//void saveCallIndel2File(string &outfilename);
		//void saveCallClipReg2File(string &outfilename);
};

#endif /* SRC_CHROME_H_ */
