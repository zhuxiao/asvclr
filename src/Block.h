#ifndef SRC_BLOCK_H_
#define SRC_BLOCK_H_

#include <iostream>
#include <string>
#include <fstream>
#include <vector>

#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/faidx.h>

#include "structures.h"
#include "Paras.h"
#include "misAlnReg.h"
#include "alnDataLoader.h"
#include "localCns.h"

//#define MIN_ADJACENT_REG_DIST		50

using namespace std;

class Block{
	public:
		Paras *paras;
		string chrname;
		int64_t chrlen, startPos, endPos, winSize;      // 1-based position
		string workdir, outCovFile;
		vector<simpleReg_t*> sub_limit_reg_vec;
		bool process_flag;

		Base *baseArr;
		vector<bam1_t*> alnDataVector;
		faidx_t *fai;

		// SNV and indel
		vector<size_t> snvVector;
		vector<reg_t*> indelVector;
		vector<misAlnReg> misAlnRegVector;
		vector<reg_t*> zeroCovRegVector;	// for long deletions

		// clip regions
		vector<reg_t*> clipRegVector;
		vector<mateClipReg_t*> mateClipRegVector;

		// detect abnormal signatures
		bool headIgnFlag, tailIgnFlag;

		double meanCov;  // excluding the gap regions filled by Ns

		// output file
		string out_dir_detect, out_dir_cns, out_dir_call;
		string snvDetectPrefix_local, indelDetectPrefix_local, clipRegDetectPrefix_local, snvFilenameDetect, indelFilenameDetect, clipRegFilenameDetect;

		string indelCnsPrefix_local, indelFilenameCns;
		ofstream *var_cand_indel_file, *misAln_reg_file, *var_cand_clipReg_file;

	public:
		Block(string chrname, size_t startPos, size_t endPos, faidx_t *fai, Paras *paras);
		virtual ~Block();
		void setOutputDir(string& out_dir_detect_prefix, string& out_dir_cns_prefix, string& out_dir_call_prefix);
		void setLimitRegs(vector<simpleReg_t*> &sub_limit_reg_vec);
		void setProcessFlag(bool process_flag);
		void setRegIngFlag(bool headIgnFlag, bool tailIgnFlag);
		void blockFillDataEst(size_t op_est);
		void blockDetect();
		void blockGenerateLocalConsWorkOpt();
		void setVarCandFiles(ofstream *var_cand_indel_file, ofstream *var_cand_clipReg_file);
		void resetVarCandFiles();
		void setMisAlnRegFile(ofstream *misAln_reg_file);
		void resetMisAlnRegFile();

	private:
		void destroyBaseArray();
		void destoryAlnData();
		void destroySnvVector();
		void destroyIndelVector();
		void destroyClipRegVector();
		void destroyZeroCovRegVector();
		void destroyMisAlnRegVector();
		Base *initBaseArray();
		int loadAlnData();
		int computeBlockBaseInfo();
		void computeBlockMeanCov();
		void fillDataEst(size_t op_est);
		void AddReadLenEstInfo();
		void AddCovDepthEstInfo();
		int outputCovFile();
		void maskMisAlnRegs();
		int computeMisAlnDisagrReg();
		void extractMisAlnRegions();
		void saveMisAlnRegToFile();
		void computeDisagrNumSingleRegion(size_t startRpos, size_t endRPos, size_t regFlag);
		bool isMisAlnReg(Region &reg);
		int computeAbSigs();
		int processSingleRegion(int64_t startRpos, int64_t endRPos, int64_t regFlag);
		void copySVEvents(Region &reg);
		void computeZeroCovReg(Region &reg);
		void updateZeroCovRegUsingIndelReg(vector<reg_t*> &zeroCovRegVec, vector<reg_t*> &indelVec);
		void updateSVRegUsingLongZeroCov();
		void updateClipRegUsingLongZeroCov(vector<reg_t*> &regVec, vector<reg_t*> &zero_cov_vec);
		void updateIndelRegUsingLongZeroCov(vector<reg_t*> &regVec, vector<reg_t*> &zero_cov_vec);
		void removeFalseIndel();
		void removeFalseSNV();
		//void sortRegVec(vector<reg_t*> &regVector);
		void blockGenerateLocalConsWorkOpt_Indel();
		void blockGenerateLocalConsWorkOpt_ClipReg();
		void saveSV2File();
		vector<simpleReg_t*> computeLimitRegsForConsWork(vector<reg_t*> &varVec, bool limit_reg_process_flag, vector<simpleReg_t*> &limit_reg_vec);
		void generateCnsWork(vector<reg_t*> &varVec, bool limit_reg_process_flag, vector<simpleReg_t*> &sub_limit_reg_vec_work, bool clip_reg_flag);

		// duplication and inversion
		void mergeOverlappedClipReg();
};

#endif /* SRC_BLOCK_H_ */
