#ifndef SRC_GENOME_H_
#define SRC_GENOME_H_

#include <iostream>
#include <string>
#include <thread>

#include "structures.h"
#include "Paras.h"
#include "Chrome.h"
#include "util.h"
#include "blatAlnTra.h"
#include "sv_sort.h"

using namespace std;

#define MIN_VALID_TRA_RATIO			(0.95f)
#define MAX_BED_COLS_NUM			12		//	maximum column number of BED file


class Genome{
	private:
		Paras *paras;
		vector<Chrome*> chromeVector;
		faidx_t *fai;
		bam_hdr_t *header;
		int64_t genomeSize;
//		int32_t minMapQ;

		// output directory and files
		string out_dir, out_dir_detect, out_dir_cns, out_dir_call, out_dir_tra, out_dir_result;

		string out_filename_detect_snv, out_filename_detect_indel, out_filename_detect_clipReg;
		string out_filename_result_snv, out_filename_result_indel, out_filename_result_clipReg, out_filename_result_tra, out_filename_result_vars, out_filename_result_vars_vcf;
		string work_finish_filename;
		string limit_reg_filename;

		//vector<varCand*> var_cand_vec;

		string blat_aln_info_filename_tra;
		vector<varCand*> blat_aligned_tra_varCand_vec;

	public:
		Genome(Paras *paras);
		virtual ~Genome();
		int generateGenomeBlocks();
		void saveGenomeBlocksToFile();
		int genomeDetect();
		int genomeLocalCons();
		int genomeCall();
		void estimateSVSizeNum();
		void saveResultVCF();

	private:
		void init();
		void saveLimitRegsToFile(string &limit_reg_filename, vector<simpleReg_t*> &limit_reg_vec);
		void loadLimitRegs();
		Chrome* allocateChrome(string& chrname, int chrlen, faidx_t *fai);
		void sortChromes(vector<Chrome*> &chr_vec, vector<Chrome*> &chr_vec_tmp);
		int getGenomeSize();
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

		void initMonitorKilledMinimap2WorkMem();
		void initMonitorKilledBlatWorkMem();
		void releaseMonitorKilledMinimap2WorkMem();
		void releaseMonitorKilledBlatWorkMem();
		void genomeCollectCallWork();
		void genomeFinishCallWork();
		int processAlnWork();
		int processBlatAlnWork();
		int processCallWork();
		void genomeLoadMateClipRegData();
		void recallIndelsFromTRA();
		bool isExistSuccessClipReg(varCand *var_cand, mateClipReg_t *clip_reg_exclude);
		void genomeCallTra();
		void generateBlatAlnFilenameTra();
		void genomeCallTraOp();
		void finalizeCallTra();
		vector<blatAlnTra*> loadBlatAlnDataTra();
		void releaseBlatAlnDataTra(vector<blatAlnTra*> &blat_aln_tra_vec);
		void blatAlnTra_st(vector<blatAlnTra*> *blat_aln_tra_vec);
		void blatAlnTra_mt(vector<blatAlnTra*> *blat_aln_tra_vec);
		varCand* constructNewVarCand(varCand *var_cand, varCand *var_cand_tmp);
		vector<int32_t> computeTraLoc(varCand *var_cand, varCand *var_cand_tmp, mateClipReg_t *clip_reg, int32_t round_num);
		vector<int32_t> getMateTraReg(int32_t query_id, int32_t start_query_pos, int32_t end_query_pos, varCand *var_cand_tmp);
		bool isValidTraReg(size_t start_query_pos, size_t end_query_pos, size_t start_mate_query_pos, size_t end_mate_query_pos);
		vector<int32_t> computeTraClippingLoc(size_t query_clip_part_flag, varCand *var_cand, vector<int32_t> &aln_idx_vec, varCand *var_cand_tmp, vector<int32_t> &mate_aln_idx_vec);
		bool computeTraCallSuccessFlag(vector<int32_t> &tra_loc_vec, vector<int32_t> &tra_loc_vec2, varCand *var_cand, varCand *var_cand_tmp,  mateClipReg_t *clip_reg, bool same_orient_flag, int32_t round_num);
		void saveTraLoc2ClipReg(mateClipReg_t *clip_reg, vector<int32_t> &tra_loc_vec, varCand *var_cand, varCand *var_cand_tmp, size_t round_num);
		void saveTraLoc2ClipRegForAlnFailure(mateClipReg_t *clip_reg);
		void mergeCloseRangeTra();
		mateClipReg_t* getNearDistedClipReg(mateClipReg_t *clip_reg, vector<Chrome*> &chrome_vec);
		vector<int32_t> computeDistsTra(mateClipReg_t *clip_reg1, mateClipReg_t *clip_reg2);
		void genomeFillVarseq();
		void genomeFillVarseqTra();
		void fillVarseqSingleMateClipReg(mateClipReg_t *clip_reg, ofstream &cns_info_file);
		void performLocalCnsTra(string &readsfilename, string &contigfilename, string &refseqfilename, string &clusterfilename, string &tmpdir, string &technology, double min_identity_match, int32_t sv_len_est, size_t num_threads_per_cns_work, vector<reg_t*> &varVec, string &chrname, string &inBamFile, faidx_t *fai, size_t cns_extend_size, ofstream &cns_info_file);
		vector<int32_t> getRefShiftSize(string &refseqfilename);
		vector<size_t> computeQueryLocTra(varCand *var_cand, mateClipReg_t *clip_reg, size_t end_flag);
		void genomeSaveCallSV2File();
		void saveTraCall2File();
		//void removeFPs();
		void mergeCallResult();
		void saveResultToVCF(string &in, string &out_vcf);
		//void saveIndelVCF(string &in, string &out_vcf);
		//void saveSnvVCF(string &in, string &out_vcf);
		void saveVCFHeader(ofstream &fp, string &sample_str);
		void saveVCFContigHeader(ofstream &fp);

		// output statistics
		void computeVarNumStatDetect();
		void computeVarNumStatCons();
		vector<int32_t> getSuccFailNumCns(string &filename);
		void computeVarNumStatCall();
		size_t getSVTypeSingleLine(string &line);
		void sortVarResults(string &infilename, int32_t filetype);
		void outputResult(string &outfilename, vector<vector<SV_item*>> &subsets, int32_t filetype);
};


#endif /* SRC_GENOME_H_ */
