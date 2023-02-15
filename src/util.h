#ifndef SRC_UTIL_H_
#define SRC_UTIL_H_

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sys/stat.h>
#include <dirent.h>
#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <signal.h>
#include <thread>

#include "structures.h"
#include "Base.h"
#include "clipReg.h"


using namespace std;

vector<string> split(const  string& s, const string& delim);
bool isExistStr(string &str, vector<string> &vec);
int copySingleFile(string &infilename, ofstream &outfile);
char getBase(string &seq, size_t pos, size_t orient);
void reverseSeq(string &seq);
void reverseComplement(string &seq);
void upperSeq(string &seq);
size_t getCtgCount(string &contigfilename);
reg_t* findVarvecItem(int32_t startPos, int32_t endPos, vector<reg_t*> &varVec);
vector<reg_t*> findVarvecItemAll(int32_t startPos, int32_t endPos, vector<reg_t*> &varVec);
reg_t* findVarvecItemExtSize(int32_t startRefPos, int32_t endRefPos, vector<reg_t*> &varVec, int32_t leftExtSize, int32_t rightExtSize);
vector<reg_t*> findVarvecItemAllExtSize(int32_t startRefPos, int32_t endRefPos, vector<reg_t*> &varVec, int32_t leftExtSize, int32_t rightExtSize);
int32_t getVectorIdx(reg_t *reg, vector<reg_t*> &varVec);
reg_t* getOverlappedRegByCallFlag(reg_t *reg, vector<reg_t*> &varVec);
reg_t* getOverlappedReg(reg_t *reg, vector<reg_t*> &varVec);
int32_t getOverlappedRegIdx(reg_t *reg, vector<reg_t*> &varVec);
int32_t getOverlappedRegIdxByCallFlag(reg_t *reg, vector<reg_t*> &varVec);
bool isOverlappedReg(reg_t* reg1, reg_t* reg2);
bool isOverlappedRegExtSize(reg_t* reg1, reg_t* reg2, int32_t leftExtSize, int32_t rightExtSize);
bool isOverlappedPos(size_t startPos1, size_t endPos1, size_t startPos2, size_t endPos2);
int32_t getOverlapSize(size_t startPos1, size_t endPos1, size_t startPos2, size_t endPos2);
bool isAdjacent(size_t startPos1, size_t endPos1, size_t startPos2, size_t endPos2, int32_t dist_thres);
bool isOverlappedMateClipReg(mateClipReg_t *mate_clip_reg1, mateClipReg_t *mate_clip_reg2);
mateClipReg_t* getOverlappedMateClipReg(mateClipReg_t *mate_clip_reg_given, vector<mateClipReg_t*> &mateClipRegVec);
bam_hdr_t* loadSamHeader(string &inBamFile);
bool isInReg(int32_t pos, vector<reg_t*> &vec);
int32_t computeDisagreeNum(Base *baseArray, int32_t arr_size);
void mergeOverlappedReg(vector<reg_t*> &regVector);
void updateReg(reg_t* reg1, reg_t* reg2);
void mergeAdjacentReg(vector<reg_t*> &regVec, size_t dist_thres);
void printRegVec(vector<reg_t*> &regVec, string header);
void printMateClipReg(mateClipReg_t *mate_clip_reg);
vector<string> getLeftRightPartChrname(mateClipReg_t *mate_clip_reg);
string preprocessPipeChar(string &cmd_str);
bool isFileExist(const string &filename);
void removeRedundantItems(vector<reg_t*> &reg_vec);
int32_t getLineCount(string &filename);
bool isBaseMatch(char ctgBase, char refBase);
bool isRegValid(reg_t *reg, int32_t min_size);
void exchangeRegLoc(reg_t *reg);
int32_t minimap2Aln(string &alnfilename, string &contigfilename, string &refseqfilename, int32_t max_proc_running_minutes);
int32_t blatAln(string &alnfilename, string &contigfilename, string &refseqfilename, int32_t max_proc_running_minutes);
int32_t getSuccessStatusSystemCmd(int32_t status);
bool isBlatAlnResultMatch(string &contigfilename, string &alnfilename);
bool isMinimap2AlnResultMatch(string &contigfilename, string &alnfilename);
int32_t getQueryNameLoc(string &query_name, vector<string> &query_name_vec);
bool isBlatAlnCompleted(string &alnfilename);
void cleanPrevAssembledTmpDir(const string &assem_dir_str, const string &dir_prefix);
string getCallFileHeaderBed(string &sample);
string getCallFileHeaderBedpe(string &sample);
assembleWork_opt* allocateAssemWorkOpt(string &chrname, string &readsfilename, string &contigfilename, string &refseqfilename, string &tmpdir, vector<reg_t*> &varVec, bool clip_reg_flag, bool limit_reg_process_flag, vector<simpleReg_t*> &limit_reg_vec);
void releaseAssemWorkOpt(assembleWork_opt *assem_work_opt);
void destroyAssembleWorkOptVec(vector<assembleWork_opt*> &assem_work_vec);
void deleteItemFromAssemWorkVec(int32_t item_id, vector<assembleWork_opt*> &assem_work_vec);
int32_t getItemIDFromAssemWorkVec(string &contigfilename, vector<assembleWork_opt*> &assem_work_vec);
void *doit_canu(void *arg);
int test_canu(int n, vector<string> &cmd_vec);
void* processSingleAssembleWork(void *arg);
void performLocalAssembly(string &readsfilename, string &contigfilename, string &refseqfilename, string &tmpdir, string &technology, string &canu_version, size_t num_threads_per_assem_work, vector<reg_t*> &varVec, string &chrname, string &inBamFile, faidx_t *fai, int32_t assembly_extend_size, ofstream &assembly_info_file, double expected_cov_assemble, double min_input_cov_canu, bool delete_reads_flag, bool keep_failed_reads_flag, bool clip_reg_flag, int32_t minClipEndSize, int32_t minConReadLen, int32_t min_sv_size, double max_seg_size_ratio, bool limit_reg_process_flag, vector<simpleReg_t*> &limit_reg_vec);
bool isReadableFile(string &filename);
void* processSingleMinimap2AlnWork(void *arg);
void* processSingleBlatAlnWork(void *arg);
void* processSingleCallWork(void *arg);
void outputAssemWorkOptToFile_debug(vector<assembleWork_opt*> &assem_work_opt_vec);
string getOldOutDirname(string &filename, string &sub_work_dir);
string getUpdatedItemFilename(string &filename, string &out_dir, string &old_out_dir);
string getChrnameByFilename(string &filename);
string deleteTailPathChar(string &dirname);
// get overlapped simple regions
vector<simpleReg_t*> getOverlappedSimpleRegsExt(string &chrname, int64_t begPos, int64_t endPos, vector<simpleReg_t*> &limit_reg_vec, int32_t ext_size);
vector<simpleReg_t*> getOverlappedSimpleRegs(string &chrname, int64_t begPos, int64_t endPos, vector<simpleReg_t*> &limit_reg_vec);
vector<simpleReg_t*> getFullyContainedSimpleRegs(string &chrname, int64_t begPos, int64_t endPos, vector<simpleReg_t*> &limit_reg_vec);
bool isFullyContainedReg(string &chrname1, int64_t begPos1, int64_t endPos1, string &chrname2, int64_t startPos2, int64_t endPos2);
vector<simpleReg_t*> getPosContainedSimpleRegs(string &chrname, int64_t begPos, int64_t endPos, vector<simpleReg_t*> &limit_reg_vec);
simpleReg_t* allocateSimpleReg(string &simple_reg_str);
void destroyLimitRegVector(vector<simpleReg_t*> &limit_reg_vec);
void printLimitRegs(vector<simpleReg_t*> &limit_reg_vec, string &description);
void getRegByFilename(simpleReg_t *reg, string &filename, string &pattern_str);
vector<simpleReg_t*> extractSimpleRegsByStr(string &regs_str);
string getLimitRegStr(vector<simpleReg_t*> &limit_reg_vec);
void createDir(string &dirname);
vector<double> getTotalHighIndelClipRatioBaseNum(Base *regBaseArr, int64_t arr_size);
vector<mismatchReg_t*> getMismatchRegVec(localAln_t *local_aln);
vector<mismatchReg_t*> getMismatchRegVecWithoutPos(localAln_t *local_aln);
void removeShortPolymerMismatchRegItems(localAln_t *local_aln, vector<mismatchReg_t*> &misReg_vec, string &inBamFile, faidx_t *fai);
void adjustVarLocByMismatchRegs(reg_t *reg, vector<mismatchReg_t*> &misReg_vec, int32_t start_aln_idx_var, int32_t end_aln_idx_var);
void releaseMismatchRegVec(vector<mismatchReg_t*> &misReg_vec);
mismatchReg_t *getMismatchReg(int32_t aln_idx, vector<mismatchReg_t*> &misReg_vec);
bool isPolymerSeq(string &seq);
vector<clipAlnData_t*> getQueryClipAlnSegs(string &queryname, vector<clipAlnData_t*> &clipAlnDataVector);
vector<clipAlnData_t*> getQueryClipAlnSegsAll(string &queryname, vector<clipAlnData_t*> &clipAlnDataVector);
bool isQuerySelfOverlap(vector<clipAlnData_t*> &query_aln_segs, int32_t maxVarRegSize);
bool isSegSelfOverlap(clipAlnData_t *clip_aln1, clipAlnData_t *clip_aln2, int32_t maxVarRegSize);
vector<int32_t> getAdjacentClipAlnSeg(int32_t arr_idx, int32_t clip_end_flag, vector<clipAlnData_t*> &query_aln_segs, int32_t minClipEndSize);
vector<BND_t*> generateBNDItems(int32_t reg_id, int32_t clip_end, int32_t checked_arr[][2], string &chrname1, string &chrname2, int64_t tra_pos_arr[4], vector<string> &bnd_str_vec, faidx_t *fai);
void checkBNDStrVec(mateClipReg_t &mate_clip_reg);
bool isValidBNDStr(int32_t reg_id, int32_t clip_end, int32_t checked_arr[][2], string &chrname1, string &chrname2, int64_t tra_pos_arr[4], vector<string> &bnd_str_vec);
void checkSuppNum(mateClipReg_t &mate_clip_reg, int32_t support_num_thres);
bool isSizeSatisfied(int64_t ref_dist, int64_t query_dist, int64_t min_sv_size_usr, int64_t max_sv_size_usr);
bool isSizeSatisfied2(int64_t sv_len, int64_t min_sv_size_usr, int64_t max_sv_size_usr);
bool isNotAlreadyExists(vector<reg_t*> &varVec, reg_t *reg);
bool isNotAlreadyExists(vector<reg_t*> &varVec, reg_t *reg, int32_t pos);
bool isDecoyChr(string &chrname);
void removeVarCandNode(varCand *var_cand, vector<varCand*> &var_cand_vec);
void *processSingleMateClipRegDetectWork(void *arg);
void processClipRegs(int32_t work_id, mateClipReg_t &mate_clip_reg, reg_t *reg, vector<mateClipReg_t*> *mateClipRegVector, vector<reg_t*> *clipRegVector, vector<varCand*> &var_cand_clipReg_vec, vector<Block*> &blockVector, Paras *paras, pthread_mutex_t *p_mutex_mate_clip_reg);
Block* computeBlocByPos_util(int64_t begPos, vector<Block*> &block_vec, Paras *paras);
int32_t computeBlocID_util(int64_t begPos, vector<Block*> &block_vec, Paras *paras);
void sortRegVec(vector<reg_t*> &regVector);
bool isRegSorted(vector<reg_t*> &regVector);
void startWorkProcessMonitor(string &work_finish_filename, string &monitoring_proc_names, int32_t max_proc_running_minutes);
void *workProcessMonitor(void *arg);

vector<struct alnSeg*> generateAlnSegs(bam1_t* b);
vector<struct alnSeg*> generateAlnSegs2(bam1_t* b, int64_t startRefPos_paras, int64_t endRefPos_paras);
vector<struct alnSeg*> generateAlnSegs_no_MD(bam1_t* b, Base* BaseArr, int64_t startRefPos_paras, int64_t endRefPos_paras);
vector<struct alnSeg*> generateAlnSegs_no_MD2(bam1_t* b, Base* BaseArr, int64_t startRefPos_paras, int64_t endRefPos_paras);
vector<struct alnSeg*> generateAlnSegs_no_MD(bam1_t* b, string &refseq, int64_t startRefPos_paras, int64_t endRefPos_paras);
vector<struct alnSeg*> generateAlnSegs_no_MD2(bam1_t* b, string &refseq, int64_t startRefPos_paras, int64_t endRefPos_paras);
vector<struct pafalnSeg*> generatePafAlnSegs(minimap2_aln_t* minimap2_aln_item, string cons_seq, string ref_seq);
struct alnSeg* allocateAlnSeg(int64_t startRpos, int64_t startQpos, int32_t seglen, uint32_t opFlag, const string &seg_MD);
struct pafalnSeg* allocatePafAlnSeg(int64_t startRpos, int64_t startQpos, int64_t startSubjectPos, int32_t seglen, uint32_t opflag, string cons_seq, string ref_seq, int32_t aln_orient);
void destroyAlnSegs(vector<struct alnSeg*> &alnSegs);
int32_t getBamType(bam1_t *b);
vector<struct MD_seg*> extractMDSegs(bam1_t* b);
struct MD_seg* allocateMDSeg(string& seg, uint32_t opflag);
void destroyMDSeg(vector<struct MD_seg*> &segs_MD);

// testing functions
void outputAlnSegs(vector<struct alnSeg*> &alnSegs);
void outputMDSeg(vector<struct MD_seg*> &segs_MD);
bool haveOpflagCigar(bam1_t* b, uint32_t opfalg);

void testAlnSegVec(string &inBamFile, faidx_t *fai);
void checkAlnSegVecSingleQuery(vector<struct alnSeg*> &alnSegs);

int32_t getOriginalQueryLen(bam1_t *b);
int64_t getEndRefPosAlnSeg(struct alnSeg *aln_seg);
int64_t getEndRefPosAlnSeg(int64_t startRpos, int32_t opflag, int32_t op_len);


class Time{
	private:
		time_t timestamp, start, end, start_all, end_all;
		struct tm *time_tm;
		struct timeval time_val;

	public:
		Time();
		virtual ~Time();
		string getTime();
		void printTime();
		void setStartTime();
		void printElapsedTime();
		void printSubCmdElapsedTime();
		void printOverallElapsedTime();
};

#endif /* SRC_UTIL_H_ */
