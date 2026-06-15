#ifndef SRC_GENOME_H_
#define SRC_GENOME_H_

#include <iostream>
#include <string>
#include <thread>
#include <htslib/faidx.h>

#include "structures.h"
#include "Paras.h"
#include "Chrome.h"
#include "util.h"
#include "sv_sort.h"

using namespace std;

#define MIN_VALID_TRA_RATIO			(0.95f)

class Genome{
	private:
		Paras *paras;
		vector<Chrome*> chromeVector;
		faidx_t *fai;
		bam_hdr_t *header;
		int64_t genomeSize, ps_num, phased_var_num;

		// output directory and files
		string out_dir, out_dir_detect, out_dir_cns, out_dir_call, out_dir_tra, out_dir_result;

		string out_filename_detect_snv, out_filename_detect_indel, out_filename_detect_clipReg;
		string out_filename_result_snv, out_filename_result_indel, out_filename_result_clipReg, out_filename_result_tra, out_filename_result_vars, out_filename_result_vars_vcf;
		string work_finish_filename_cns, work_finish_filename_call;
		string limit_reg_filename;

	public:
		Genome(Paras *paras);
		virtual ~Genome();
		int generateGenomeBlocks();
		int genomeDetect();
		int genomeLocalCons();
		int genomeCall();
		void estimateSVSizeNum();
		void saveResultVCF();
		void removeTempResults();

	private:
		void init();
		void genomeSetMaxOpenFileNum(size_t chr_num);
		void saveLimitRegsToFile(string &limit_reg_filename, vector<simpleReg_t*> &limit_reg_vec);
		void loadLimitRegs();
		Chrome* allocateChrome(string& chrname, int chrlen, faidx_t *fai);
		void sortChromes(vector<Chrome*> &chr_vec, vector<Chrome*> &chr_vec_tmp);
		void destroyChromeVector();
		int computeCoverage();
		void removeRedundantTra();
		void removeInvalidMateClipItem();
		void removeRedundantMateClipReg();
		void genomeRemoveRedundantClipReg(Chrome *chr, vector<Chrome*> &chr_vec);
		void removeOverlappedIndelFromMateClipReg();
		void genomeRemoveFPIndelSnvInClipReg(Chrome *chr, vector<Chrome*> &chr_vec);
		mateClipReg_t* genomeGetOverlappedMateClipReg(mateClipReg_t *clip_reg_given, vector<Chrome*> &chrome_vec);
		mateClipReg_t* getSameClipRegTRA(mateClipReg_t *clip_reg_given, vector<Chrome*> &chrome_vec);
		void saveDetectResultToFile();
		void mergeDetectResult();
		void genomeLoadDataCons();
		int processConsWork();
		ofstream* getVarcandFile(string &chrname, vector<Chrome*> &chrome_vec, bool clip_reg_flag);
		void generateFile(string &filename);

		void genomeCollectCallWork();
		int processAlnWork();
		int processCallWork();
		void genomeLoadMateClipRegData();
		void genomePhasing();
		void genomeSaveCallSV2File();
		void saveTraCall2File();
		void mergeCallResult();
		void saveResultToVCF(string &in, string &out_vcf);
		void saveVCFHeader(ofstream &fp, string &sample_str);
		void saveVCFContigHeader(ofstream &fp);

		// output statistics
		void computeVarNumStatDetect();
		void computeVarNumStatCons();
		vector<int32_t> getSuccFailNumCns(string &filename);
		void computeVarNumStatCall();
		size_t getSVTypeSingleLine(string &line);
		void sortVarResults(string &infilename, int32_t filetype);
		void computePhaseSetNum(vector<vector<SV_item*>> &subsets);
		vector<int32_t> computePhaseSetNumSubset(vector<SV_item*> &subset);
		string extractPSFromVCFLine(string &line);
		void outputResult(string &outfilename, vector<vector<SV_item*>> &subsets, int32_t filetype);
		void removeTempMonitorFiles();
};


#endif /* SRC_GENOME_H_ */
