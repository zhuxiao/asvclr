#ifndef SRC_BLOCK_H_
#define SRC_BLOCK_H_

#include <iostream>
#include <string>
#include <fstream>
#include <vector>

#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/faidx.h>

#include "Region.h"
#include "Paras.h"
#include "misAlnReg.h"
#include "LocalAssembly.h"
#include "alnDataLoader.h"
#include "covLoader.h"
#include "clipReg.h"

#define MIN_ADJACENT_REG_DIST		50

using namespace std;

class Block{
	public:
		Paras *paras;
		string chrname;
		size_t chrlen, startPos, endPos, winSize;      // 1-based position
		string workdir, outCovFile;

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
		string out_dir_detect, out_dir_assemble, out_dir_call;
		string snvDetectPrefix_local, indelDetectPrefix_local, clipRegDetectPrefix_local, snvFilenameDetect, indelFilenameDetect, clipRegFilenameDetect;

		string indelAssemblePrefix_local, indelFilenameAssemble;
		ofstream *var_cand_indel_file, *misAln_reg_file, *var_cand_clipReg_file;

	public:
		Block(string chrname, size_t chrlen, size_t startPos, size_t endPos, faidx_t *fai, Paras *paras);
		virtual ~Block();
		void setOutputDir(string& out_dir_detect_prefix, string& out_dir_assemble_prefix, string& out_dir_call_prefix);
		void setRegIngFlag(bool headIgnFlag, bool tailIgnFlag);
		void blockFillDataEst(size_t op_est);
		void blockDetect();
		void blockLocalAssemble();
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
		void destroyMisAlnRegVector();
		Base *initBaseArray();
		int loadAlnData();
		int computeBlockBaseInfo();
		void computeBlockMeanCov();
		void fillDataEst(size_t op_est);
		int outputCovFile();
		void maskMisAlnRegs();
		int computeMisAlnDisagrReg();
		void extractMisAlnRegions();
		void saveMisAlnRegToFile();
		void computeDisagrNumSingleRegion(size_t startRpos, size_t endRPos, size_t regFlag);
		bool isMisAlnReg(Region &reg);
		int computeAbSigs();
		int processSingleRegion(size_t startRpos, size_t endRPos, size_t regFlag);
		void copySVEvents(Region &reg);
		void computeZeroCovReg(Region &reg);
		void updateZeroCovRegUsingIndelReg(vector<reg_t*> &zeroCovRegVec, vector<reg_t*> &indelVec);
		void updateSVRegUsingLongZeroCov();
		void updateClipRegUsingLongZeroCov(vector<reg_t*> &regVec, vector<reg_t*> &zero_cov_vec);
		void updateIndelRegUsingLongZeroCov(vector<reg_t*> &regVec, vector<reg_t*> &zero_cov_vec);
		void removeFalseIndel();
		void removeFalseSNV();
		void sortRegVec(vector<reg_t*> &regVector);
		void blockLocalAssembleIndel();
		void blockLocalAssembleClipReg();
		void performLocalAssembly(string &readsfilename, string &contigfilename, string &refseqfilename, string &tmpdir, vector<reg_t*> &varVec, string &chrname, string &inBamFile, faidx_t *fai, string &canu_version, ofstream &assembly_info_file);
		void saveSV2File();

		// duplication and inversion
		void mergeOverlappedClipReg();
};

#endif /* SRC_BLOCK_H_ */
