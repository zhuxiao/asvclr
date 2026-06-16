#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <string>
#include <pthread.h>
#include <htslib/thread_pool.h>

#include "Genome.h"
#include "Thread.h"

using namespace std;
extern pthread_mutex_t mutex_fai;

//class MultiThread;

//int32_t sig_num_arr[1001] = {0};

// Constructor with parameters
Genome::Genome(Paras *paras){
	this->paras = paras;
	init();

//	cout << "Checking BAM AlnSegs ..." << endl;
//	testAlnSegVec(paras->inBamFile, fai);
//	cout << "Checking BAM AlnSegs finished." << endl;
}

//Destructor
Genome::~Genome(){
	destroyChromeVector();
	fai_destroy(fai);
	bam_hdr_destroy(header);
}

// initialization
void Genome::init(){
	Chrome *chr;
	vector<Chrome*> chr_vec_tmp;
	string chrname_tmp, result_prefix;
	int64_t chrlen_tmp;
	bool exist_flag;

	ps_num = 0;
	out_dir = paras->outDir;
	if(out_dir.size()>0){
		//mkdir(out_dir.c_str(), S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
		createDir(out_dir);   // create the output directory

		out_dir_detect = out_dir + "/" + paras->out_dir_detect;
		out_dir_cns = out_dir + "/" + paras->out_dir_cns;
		out_dir_call = out_dir + "/" + paras->out_dir_call;
		out_dir_tra = out_dir + "/" + paras->out_dir_tra;
		out_dir_result = out_dir + "/" + paras->out_dir_result;
	}else{
		out_dir_detect = paras->out_dir_detect;
		out_dir_cns = paras->out_dir_cns;
		out_dir_call = paras->out_dir_call;
		out_dir_tra = paras->out_dir_tra;
		out_dir_result = paras->out_dir_result;
	}

	mkdir(out_dir_detect.c_str(), S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);  // create the directory  for detect command
	mkdir(out_dir_cns.c_str(), S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);  // create directory for consensus command
	mkdir(out_dir_call.c_str(), S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);  // create directory for call command
	mkdir(out_dir_tra.c_str(), S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);  // create the directory for TRA

	result_prefix = "";
	if(paras->outFilePrefix.size()) result_prefix = paras->outFilePrefix + "_";
	out_filename_detect_snv = out_dir_detect + "/" + result_prefix + "SNV_candidates";
	out_filename_detect_indel = out_dir_detect + "/" + result_prefix + "INDEL_candidate";
	out_filename_detect_clipReg = out_dir_detect + "/" + result_prefix + "clipReg_candidate";
	out_filename_result_snv = out_dir_result + "/" + result_prefix + "snv";
	out_filename_result_indel = out_dir_result + "/" + result_prefix + "indelReg.bed";
	out_filename_result_clipReg = out_dir_result + "/" + result_prefix + "clipReg.bed";
	out_filename_result_tra = out_dir_result + "/" + result_prefix + "tra.bedpe";
	if(paras->outFilePrefix.compare(RESULT_PREFIX_DEFAULT)==0){
		out_filename_result_vars = out_dir_result + "/" + result_prefix + "sv.bed";
		out_filename_result_vars_vcf = out_dir_result + "/" + result_prefix + "sv.vcf";
	}else{
		out_filename_result_vars = out_dir_result + "/" + result_prefix + "asvclr.bed";
		out_filename_result_vars_vcf = out_dir_result + "/" + result_prefix + "asvclr.vcf";
	}

	work_finish_filename_cns = out_dir + "/" + "cns_work_finished";
	work_finish_filename_call = out_dir + "/" + "call_work_finished";

	limit_reg_filename = out_dir_detect + "/" + paras->limit_reg_filename;

	if(paras->limit_reg_process_flag) saveLimitRegsToFile(limit_reg_filename, paras->limit_reg_vec);
	else {
		if(paras->command.compare(CMD_DET_STR)==0 or paras->command.compare(CMD_ALL_STR)==0 or paras->command.compare(CMD_DET_CNS_STR)==0){
			if(isFileExist(limit_reg_filename)) remove(limit_reg_filename.c_str());
		}else loadLimitRegs();
	}

	// load the fai
	fai = fai_load(paras->refFile.c_str());
	if ( !fai ) {
		cerr << __func__ << ": could not load fai index of " << paras->refFile << endl;;
		exit(1);
	}

	// load the sam/bam header
	header = loadSamHeader(paras->inBamFile);

//	if(header->n_targets!=faidx_nseq(fai))
//		cout << "Warn: the number of sequences in reference is not the same with the number of sequences of BAM header, the same reference is recommended for better performance." << endl;

	// confirm chrlen
	for(int i=0; i<header->n_targets; i++){
		chrname_tmp = header->target_name[i];
		exist_flag = faidx_has_seq(fai, chrname_tmp.c_str());
		if(exist_flag){
			chrlen_tmp = faidx_seq_len64(fai, chrname_tmp.c_str());
			if(chrlen_tmp!=header->target_len[i]){
				cerr << "The sequence length of " << chrname_tmp << " is different between the reference FASTA file and the BAM file, please confirm whether the reference and BAM file are match before running the command." << endl;
				exit(1);
			}
		}
	}

	// allocate each genome
	for(int i=0; i<header->n_targets; i++){
		chrname_tmp = header->target_name[i];
		chr = allocateChrome(chrname_tmp, header->target_len[i], fai);
		chr->setOutputDir(out_dir_detect, out_dir_cns, out_dir_call);
		chr_vec_tmp.push_back(chr);
	}
	chr_vec_tmp.shrink_to_fit();

	// sort chromosomes
	sortChromes(chromeVector, chr_vec_tmp);

	genomeSetMaxOpenFileNum(chromeVector.size());
}

// set the current maximum open files if necessary
void Genome::genomeSetMaxOpenFileNum(size_t chr_num){
	rlim_t currentLimit = getMaxOpenFileNum();
	//cout << ">>>>>>>>>>>>>> Current maximum open files: " << currentLimit << endl;
	if(currentLimit<3*chr_num){
		if(setMaxOpenFileNum(3*chr_num)==false){
			cout << ">>>>>>>>>> Failed to set maximum open files" << endl;
		}
//		currentLimit = getMaxOpenFileNum();
//		cout << ">>>>>>>>>>>>>> Current maximum open files: " << currentLimit << endl;
	}
}

// save limit regions to file
void Genome::saveLimitRegsToFile(string &limit_reg_filename, vector<simpleReg_t*> &limit_reg_vec){
	ofstream outfile;
	simpleReg_t *simple_reg;
	string line;

	outfile.open(limit_reg_filename);
	if(!outfile.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file: " << limit_reg_filename << endl;
		exit(1);
	}

	outfile << "#Chr\tStart\tEnd" << endl;
	for(size_t i=0; i<limit_reg_vec.size(); i++){
		simple_reg = limit_reg_vec.at(i);
		line = simple_reg->chrname;
		if(simple_reg->startPos!=-1 and simple_reg->endPos!=-1) line += "\t" + to_string(simple_reg->startPos) + "\t" + to_string(simple_reg->endPos);
		outfile << line << endl;
	}
	outfile.close();
}

// load limit regions
void Genome::loadLimitRegs(){
	ifstream infile;
	string line, simple_reg_str, desc;
	vector<string> str_vec;
	simpleReg_t *simple_reg;

	// check the file
	if(isFileExist(limit_reg_filename)){
		infile.open(limit_reg_filename);
		if(!infile.is_open()){
			cerr << __func__ << ": cannot open file " << limit_reg_filename << ", error." << endl;
			exit(1);
		}

		while(getline(infile, line)){
			if(line.size()>0 and line.at(0)!='#'){
				str_vec = split(line, "\t");
				simple_reg_str = str_vec.at(0);
				if(str_vec.size()==3) simple_reg_str += ":" + str_vec.at(1) + "-" + str_vec.at(2);
				simple_reg = allocateSimpleReg(simple_reg_str);
				if(simple_reg) paras->limit_reg_vec.push_back(simple_reg);
			}
		}
		infile.close();

		if(paras->limit_reg_vec.size()) {
			paras->limit_reg_process_flag = true;
			desc = "Load limit regions to process:";
			printLimitRegs(paras->limit_reg_vec, desc);
		}
	}
}

// allocate the Chrome node
Chrome* Genome::allocateChrome(string& chrname, int chrlen, faidx_t *fai){
	Chrome *chr_tmp;
	chr_tmp = new Chrome(chrname, chrlen, fai, paras, &chromeVector);
	if(!chr_tmp){ cerr << "Genome: cannot allocate memory" << endl; exit(1); }
	return chr_tmp;
}

// sort chromosomes
void Genome::sortChromes(vector<Chrome*> &chr_vec, vector<Chrome*> &chr_vec_tmp){
	int8_t *selected_flag_array;
	size_t i, j;
	int32_t idx_chr;
	Chrome *chr;
	string chr_str1, chr_str2, head_str, head_str_chr;
	vector<string> chr_str_vec;
	bool have_chr_prefix_flag;

	chr_str1 = "chr1_chr2_chr3_chr4_chr5_chr6_chr7_chr8_chr9_chr10_chr11_chr12_chr13_chr14_chr15_chr16_chr17_chr18_chr19_chr20_chr21_chr22_chrX_chrY_chrM";
	chr_str2 = "1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_X_Y_MT";

	selected_flag_array = (int8_t*) calloc(chr_vec_tmp.size(), sizeof(int8_t));

	// determine the have_chr_prefix_flag
	have_chr_prefix_flag = false;
	for(i=0; i<chr_vec_tmp.size(); i++){
		chr = chr_vec_tmp.at(i);
		if(chr->chrname.compare("chr1")==0){
			have_chr_prefix_flag = true;
			break;
		}
	}

	if(have_chr_prefix_flag) chr_str_vec = split(chr_str1, "_");
	else chr_str_vec = split(chr_str2, "_");

	// select chromosomes by 'chr'
	for(i=0; i<chr_str_vec.size(); i++){
		idx_chr = -1;
		for(j=0; j<chr_vec_tmp.size(); j++){
			chr = chr_vec_tmp.at(j);
			if(chr->chrname.compare(chr_str_vec.at(i))==0){
				idx_chr = j;
				break;
			}
		}
		if(idx_chr!=-1){
			chr = chr_vec_tmp.at(idx_chr);
			chr_vec.push_back(chr);
			selected_flag_array[idx_chr] = 1;
		}
	}

	// select chromosomes by 'chrA_*'
	for(i=0; i<chr_str_vec.size(); i++){
		head_str = chr_str_vec.at(i) + "_";
		for(j=0; j<chr_vec_tmp.size(); j++){
			chr = chr_vec_tmp.at(j);
			head_str_chr = chr->chrname.substr(0, head_str.size());
			if(head_str_chr.compare(head_str)==0){
				chr_vec.push_back(chr);
				selected_flag_array[j] = 1;
			}
		}
	}

	// add unselected items
	for(i=0; i<chr_vec_tmp.size(); i++){
		if(selected_flag_array[i]==0){
			chr = chr_vec_tmp.at(i);
			chr_vec.push_back(chr);
		}
	}

	vector<Chrome*>().swap(chr_vec_tmp);	// delete all items

	free(selected_flag_array);
}

// release each chrome in chromeVector
void Genome::destroyChromeVector(){
	vector<Chrome*>::iterator chr;
	for(chr=chromeVector.begin(); chr!=chromeVector.end(); chr++)
		delete (*chr);   // release each chrome
	vector<Chrome*>().swap(chromeVector);
}

// generate the genome blocks
int Genome::generateGenomeBlocks(){
	vector<Chrome*>::iterator chr;
	for(chr=chromeVector.begin(); chr!=chromeVector.end(); chr++)
		(*chr)->generateChrBlocks();
	return 0;
}

// estimate the insertion/deletion/clipping parameters
void Genome::estimateSVSizeNum(){
	Chrome *chr;
	size_t i;

	cout << "Estimating parameters:" << endl;

	// initialize the data
	paras->initEst();

	// fill the size data using each chromosome
	paras->reg_sum_size_est = 0;
	//paras->total_depth = 0;
	for(i=0; i<chromeVector.size(); i++){
		chr = chromeVector.at(i);
		if(chr->valid_flag){
			chr->chrFillDataEst(SIZE_EST_OP);
			if(paras->reg_sum_size_est>=paras->max_reg_sum_size_est) break;
		}
		//paras->total_depth += paras->chr_mean_depth;
	}
	//paras->chrome_num = chromeVector.size() - 1;

	// size estimate
	paras->estimate(SIZE_EST_OP);

	// fill the num data using each chromosome
	paras->reg_sum_size_est = 0;
	for(i=0; i<chromeVector.size(); i++){
		chr = chromeVector.at(i);
		if(chr->valid_flag){
			chr->chrFillDataEst(NUM_EST_OP);
			if(paras->reg_sum_size_est>=paras->max_reg_sum_size_est) break;
		}
	}
	// num estimate
	paras->estimate(NUM_EST_OP);
}

// detect variants for genome
int Genome::genomeDetect(){
	Chrome *chr;

	for(size_t i=0; i<chromeVector.size(); i++){
		chr = chromeVector.at(i);
		if((chr->decoy_flag==false or paras->include_decoy) and (chr->alt_flag==false or paras->include_alt) and chr->valid_flag){
			//cout << chr->chrname << endl;
			chr->chrDetect();
		}
	}

	Time time;
	cout << "[" << time.getTime() << "]: Finalizing detect ..." << endl;

	// remove redundant translocations
	removeRedundantTra();

	// remove redundant mate clipping regions
	removeRedundantMateClipReg();

	//cout << "[" << time.getTime() << "]: removeOverlappedIndelFromMateClipReg() ..." << endl;
	// remove overlapped indels from mate clipping regions
	removeOverlappedIndelFromMateClipReg();

	//cout << "[" << time.getTime() << "]: saveDetectResultToFile() ..." << endl;
	// save detect result to file for each chrome
	saveDetectResultToFile();

	mergeDetectResult();

	// compute statistics for detect command
	computeVarNumStatDetect();

	return 0;
}

// remove repeatedly detected translocations
void Genome::removeRedundantTra(){
	size_t i, j, clipPosNum, clipPosNum_overlapped;
	Chrome *chr;
	mateClipReg_t *clip_reg, *clip_reg_ret;
	for(i=0; i<chromeVector.size(); i++){
		chr = chromeVector.at(i);
		if(chr->valid_flag==false) continue;
		for(j=0; j<chr->mateClipRegVector.size(); j++){
			clip_reg = chr->mateClipRegVector.at(j);
			if(clip_reg->valid_flag and clip_reg->reg_mated_flag and clip_reg->sv_type==VAR_TRA){
				if(clip_reg->leftClipPosTra1==-1 and clip_reg->leftClipPosTra2==-1 and clip_reg->rightClipPosTra1==-1 and clip_reg->rightClipPosTra2==-1){
					clip_reg_ret = genomeGetOverlappedMateClipReg(clip_reg, chromeVector);
					if(clip_reg_ret){
						clipPosNum = clip_reg->leftClipPosNum + clip_reg->leftClipPosNum2 + clip_reg->rightClipPosNum + clip_reg->rightClipPosNum2;
						clipPosNum_overlapped = clip_reg_ret->leftClipPosNum + clip_reg_ret->leftClipPosNum2 + clip_reg_ret->rightClipPosNum + clip_reg_ret->rightClipPosNum2;
						if(clipPosNum>=clipPosNum_overlapped) clip_reg_ret->valid_flag = false;
						else clip_reg->valid_flag = false;
					}
				}else{
					clip_reg_ret = getSameClipRegTRA(clip_reg, chromeVector);
					if(clip_reg_ret){
						if(clip_reg->leftClipRegNum==1 and clip_reg->rightClipRegNum==1){
							if(clip_reg->leftMeanClipPos<clip_reg_ret->leftMeanClipPos){
								clip_reg->leftClipReg2 = clip_reg_ret->leftClipReg;
								clip_reg->leftMeanClipPos2 = clip_reg_ret->leftMeanClipPos;
							}else{
								clip_reg->leftClipReg2 = clip_reg->leftClipReg;
								clip_reg->leftMeanClipPos2 = clip_reg->leftMeanClipPos;
								clip_reg->leftClipReg = clip_reg_ret->leftClipReg;
								clip_reg->leftMeanClipPos = clip_reg_ret->leftMeanClipPos;
							}
							clip_reg->leftClipRegNum ++;
							clip_reg_ret->leftClipReg = NULL;

							if(clip_reg->rightMeanClipPos<clip_reg_ret->rightMeanClipPos){
								clip_reg->rightClipReg2 = clip_reg_ret->rightClipReg;
								clip_reg->rightMeanClipPos2 = clip_reg_ret->rightMeanClipPos;
							}else{
								clip_reg->rightClipReg2 = clip_reg->rightClipReg;
								clip_reg->rightMeanClipPos2 = clip_reg->rightMeanClipPos;
								clip_reg->rightClipReg = clip_reg_ret->rightClipReg;
								clip_reg->rightMeanClipPos = clip_reg_ret->rightMeanClipPos;
							}
							clip_reg->rightClipRegNum ++;
							clip_reg_ret->rightClipReg = NULL;
							clip_reg_ret->valid_flag = false;
						}
					}
				}
			}
		}
	}

	// remove invalid elements
	removeInvalidMateClipItem();
}

// remove invalid mate clip region items
void Genome::removeInvalidMateClipItem(){
	size_t i, j;
	mateClipReg_t *clip_reg;
	Chrome *chr;

	for(i=0; i<chromeVector.size(); i++){
		chr = chromeVector.at(i);
		if(chr->valid_flag){
			for(j=0; j<chr->mateClipRegVector.size(); ){
				clip_reg = chr->mateClipRegVector.at(j);
				if(clip_reg->valid_flag==false){
					// cout << "\t";
					// if(clip_reg->leftClipReg) { cout << ", leftClipReg=" << clip_reg->leftClipReg->chrname << ":" << clip_reg->leftClipReg->startRefPos << "-" << clip_reg->leftClipReg->endRefPos; }
					// if(clip_reg->leftClipReg2) { cout << ", leftClipReg2=" << clip_reg->leftClipReg2->chrname << ":" << clip_reg->leftClipReg2->startRefPos << "-" << clip_reg->leftClipReg2->endRefPos; }
					// if(clip_reg->rightClipReg) { cout << ", rightClipReg=" << clip_reg->rightClipReg->chrname << ":" << clip_reg->rightClipReg->startRefPos << "-" << clip_reg->rightClipReg->endRefPos; }
					// if(clip_reg->rightClipReg2) { cout << ", rightClipReg2=" << clip_reg->rightClipReg2->chrname << ":" << clip_reg->rightClipReg2->startRefPos << "-" << clip_reg->rightClipReg2->endRefPos; }
					// if(clip_reg->largeIndelClipReg) { cout << ", largeIndelClipReg=" << clip_reg->largeIndelClipReg->chrname << ":" << clip_reg->largeIndelClipReg->startRefPos << "-" << clip_reg->largeIndelClipReg->endRefPos << ", sv_len=" << clip_reg->largeIndelClipReg->sv_len << ", var_type=" << clip_reg->sv_type; }
					// cout << endl;

					if(clip_reg->leftClipReg) { delete clip_reg->leftClipReg; clip_reg->leftClipReg = NULL; }
					if(clip_reg->leftClipReg2) { delete clip_reg->leftClipReg2; clip_reg->leftClipReg2 = NULL; }
					if(clip_reg->rightClipReg) { delete clip_reg->rightClipReg; clip_reg->rightClipReg = NULL; }
					if(clip_reg->rightClipReg2) { delete clip_reg->rightClipReg2; clip_reg->rightClipReg2 = NULL; }
					if(clip_reg->largeIndelClipReg) { delete clip_reg->largeIndelClipReg; clip_reg->largeIndelClipReg = NULL; }
					if(clip_reg->var_cand) { chr->removeVarCandNodeClipReg(clip_reg->var_cand); clip_reg->var_cand = NULL; } // free item
					if(clip_reg->left_var_cand_tra) { chr->removeVarCandNodeClipReg(clip_reg->left_var_cand_tra); clip_reg->left_var_cand_tra = NULL; }  // free item
					if(clip_reg->right_var_cand_tra) { chr->removeVarCandNodeClipReg(clip_reg->right_var_cand_tra); clip_reg->right_var_cand_tra = NULL; }  // free item
					delete clip_reg;
					chr->mateClipRegVector.erase(chr->mateClipRegVector.begin()+j);
				}else j++;
			}
		}
	}
}

// remove redundant mate clipping regions
void Genome::removeRedundantMateClipReg(){
	for(size_t i=0; i<chromeVector.size(); i++){
		if(chromeVector.at(i)->valid_flag)
			genomeRemoveRedundantClipReg(chromeVector.at(i), chromeVector);
	}

	// remove invalid elements
	removeInvalidMateClipItem();
}

// remove redundant mate clipping regions
void Genome::genomeRemoveRedundantClipReg(Chrome *chr, vector<Chrome*> &chr_vec){
	size_t i;
	vector<mateClipReg_t*> mate_clip_reg_vec;
	mateClipReg_t *mate_clip_reg, *mate_clip_reg_ret;
	int32_t clipPosNum, clipPosNum_overlapped;

	mate_clip_reg_vec = chr->mateClipRegVector;
	for(i=0; i<mate_clip_reg_vec.size(); i++){
		mate_clip_reg = mate_clip_reg_vec.at(i);
		if(mate_clip_reg->valid_flag){
			mate_clip_reg_ret = genomeGetOverlappedMateClipReg(mate_clip_reg, chromeVector);
			if(mate_clip_reg_ret){
				if(mate_clip_reg->large_indel_flag==false and mate_clip_reg_ret->large_indel_flag==false){ // clipping regions
					clipPosNum = mate_clip_reg->leftClipPosNum + mate_clip_reg->leftClipPosNum2 + mate_clip_reg->rightClipPosNum + mate_clip_reg->rightClipPosNum2;
					clipPosNum_overlapped = mate_clip_reg_ret->leftClipPosNum + mate_clip_reg_ret->leftClipPosNum2 + mate_clip_reg_ret->rightClipPosNum + mate_clip_reg_ret->rightClipPosNum2;
				}else if(mate_clip_reg->large_indel_flag==true and mate_clip_reg_ret->large_indel_flag==true){ // large indel
					clipPosNum = mate_clip_reg->supp_num_largeIndel;
					clipPosNum_overlapped = mate_clip_reg_ret->supp_num_largeIndel;
				}
				if(clipPosNum>=clipPosNum_overlapped) mate_clip_reg_ret->valid_flag = false;
				else mate_clip_reg->valid_flag = false;
			}
		}
	}
}

void Genome::removeOverlappedIndelFromMateClipReg(){
	for(size_t i=0; i<chromeVector.size(); i++){
		if(chromeVector.at(i)->valid_flag)
			genomeRemoveFPIndelSnvInClipReg(chromeVector.at(i), chromeVector);
	}
}

// remove FP Indels and SNVs in clipping regions
void Genome::genomeRemoveFPIndelSnvInClipReg(Chrome *chr, vector<Chrome*> &chr_vec){
	size_t i, j, k;
	Chrome *chr_tmp;
	vector<mateClipReg_t*> mate_clip_reg_vec;
	mateClipReg_t *mate_clip_reg;
	string chrname_arr[4];
	int64_t pos_arr[4], pos;
	Chrome *chr_arr[4];
	Block *block_arr[4], *block;
	reg_t *reg;
	bool flag;

	//size_t num = 0;
	mate_clip_reg_vec = chr->mateClipRegVector;
	for(i=0; i<mate_clip_reg_vec.size(); i++){
		mate_clip_reg = mate_clip_reg_vec.at(i);
		if(mate_clip_reg->valid_flag){
			for(j=0; j<4; j++) { chrname_arr[j] = ""; pos_arr[j] = 0; block_arr[j] = NULL; chr_arr[j] = NULL; }

			// initialize
			if(mate_clip_reg->leftClipReg) { chrname_arr[0] = mate_clip_reg->leftClipReg->chrname; pos_arr[0] = mate_clip_reg->leftClipReg->startRefPos; }
			if(mate_clip_reg->leftClipReg2) { chrname_arr[1] = mate_clip_reg->leftClipReg2->chrname; pos_arr[1] = mate_clip_reg->leftClipReg2->startRefPos; }
			if(mate_clip_reg->rightClipReg) { chrname_arr[2] = mate_clip_reg->rightClipReg->chrname; pos_arr[2] = mate_clip_reg->rightClipReg->startRefPos; }
			if(mate_clip_reg->rightClipReg2) { chrname_arr[3] = mate_clip_reg->rightClipReg2->chrname; pos_arr[3] = mate_clip_reg->rightClipReg2->startRefPos; }

			// get the chromes
			for(j=0; j<4; j++){
				if(j>0 and chrname_arr[j].compare(chrname_arr[j-1])==0){
					chr_arr[j] = chr_arr[j-1];
				}else{
					for(k=0; k<chr_vec.size(); k++){
						chr_tmp = chr_vec.at(k);
						if(chr!=chr_tmp and chrname_arr[j].size()>0 and chr_tmp->chrname.compare(chrname_arr[j])==0 and chr_tmp->valid_flag){
							chr_arr[j] = chr_tmp;
							break;
						}
					}
				}
			}

			// get the blocks
			for(j=0; j<4; j++){
				chr_tmp = chr_arr[j];
				if(chr_tmp){
					block_arr[j] = chr_tmp->computeBlocByPos(pos_arr[j], chr_tmp->blockVector);
				}
			}

			for(j=0; j<4; j++){
				block = block_arr[j];
				if(block){
					for(k=0; k<block->indelVector.size(); ){
						reg = block->indelVector.at(k);
						flag = isIndelInSingleClipReg(reg, mate_clip_reg);
						if(flag){
							//printRegVec(block->indelVector, "indelVector");
							delete reg;
							block->indelVector.erase(block->indelVector.begin()+k);
							//num ++;
						}else k++;
					}
					for(k=0; k<block->snvVector.size(); ){
						pos = block->snvVector.at(k);
						flag = isSnvInSingleClipReg(block->chrname, pos, mate_clip_reg);
						if(flag) block->snvVector.erase(block->snvVector.begin()+k);
						else k++;
					}
				}
			}
		}
	}

	//cout << "------------- num=" << num << endl;
}

// get overlapped mate clip regions
mateClipReg_t* Genome::genomeGetOverlappedMateClipReg(mateClipReg_t *clip_reg_given, vector<Chrome*> &chrome_vec){
	mateClipReg_t *clip_reg_overlapped = NULL;
	Chrome *chr;
	for(size_t i=0; i<chrome_vec.size(); i++){
		chr = chrome_vec.at(i);
		clip_reg_overlapped = getOverlappedMateClipReg(clip_reg_given, chr->mateClipRegVector);
		if(clip_reg_overlapped)
			break;
	}
	return clip_reg_overlapped;
}

mateClipReg_t* Genome::getSameClipRegTRA(mateClipReg_t *clip_reg_given, vector<Chrome*> &chrome_vec){
	size_t i, j;
	mateClipReg_t *clip_reg, *clip_reg_ret = NULL;
	Chrome *chr;
	bool same_flag = true;

	for(i=0; i<chrome_vec.size(); i++){
		chr = chrome_vec.at(i);
		for(j=0; j<chr->mateClipRegVector.size(); j++){
			clip_reg = chr->mateClipRegVector.at(j);
			if(clip_reg==clip_reg_given) continue;

			if(clip_reg->leftClipRegNum==clip_reg_given->leftClipRegNum and clip_reg->rightClipRegNum==clip_reg_given->rightClipRegNum){
				if(clip_reg->leftClipPosTra1!=-1 and clip_reg_given->leftClipPosTra1!=-1){
					if(clip_reg->chrname_leftTra1.compare(clip_reg_given->chrname_leftTra1)!=0 or clip_reg->leftClipPosTra1<clip_reg_given->leftClipPosTra1-CLIP_END_EXTEND_SIZE or clip_reg->leftClipPosTra1>clip_reg_given->leftClipPosTra1+CLIP_END_EXTEND_SIZE)
						same_flag = false;
				}else same_flag = false;
				if(same_flag and clip_reg->rightClipPosTra1!=-1 and clip_reg_given->rightClipPosTra1!=-1){
					if(clip_reg->chrname_rightTra1.compare(clip_reg_given->chrname_rightTra1)!=0 or clip_reg->rightClipPosTra1<clip_reg_given->rightClipPosTra1-CLIP_END_EXTEND_SIZE or clip_reg->rightClipPosTra1>clip_reg_given->rightClipPosTra1+CLIP_END_EXTEND_SIZE)
						same_flag = false;
				}else same_flag = false;
				if(same_flag and clip_reg->leftClipPosTra2!=-1 and clip_reg_given->leftClipPosTra2!=-1){
					if(clip_reg->chrname_leftTra2.compare(clip_reg_given->chrname_leftTra2)!=0 or clip_reg->leftClipPosTra2<clip_reg_given->leftClipPosTra2-CLIP_END_EXTEND_SIZE or clip_reg->leftClipPosTra2>clip_reg_given->leftClipPosTra2+CLIP_END_EXTEND_SIZE)
						same_flag = false;
				}else same_flag = false;
				if(same_flag and clip_reg->rightClipPosTra2!=-1 and clip_reg_given->rightClipPosTra2!=-1){
					if(clip_reg->chrname_rightTra2.compare(clip_reg_given->chrname_rightTra2)!=0 or clip_reg->rightClipPosTra2<clip_reg_given->rightClipPosTra2-CLIP_END_EXTEND_SIZE or clip_reg->rightClipPosTra2>clip_reg_given->rightClipPosTra2+CLIP_END_EXTEND_SIZE)
						same_flag = false;
				}else same_flag = false;
			}else same_flag = false;

			if(same_flag){
				clip_reg_ret = clip_reg;
				break;
			}
		}
		if(same_flag) break;
	}

	return clip_reg_ret;
}

// save detect result to file for each chrome
void Genome::saveDetectResultToFile(){
	Chrome *chr;
	for(size_t i=0; i<chromeVector.size(); i++){
		chr = chromeVector.at(i);
		if(chr->valid_flag)
			chr->chrMergeDetectResultToFile();
	}
}

// merge detect result to single file
void Genome::mergeDetectResult(){
	Chrome *chr;
	ofstream out_file_snv, out_file_indel, out_file_clipReg;

	out_file_snv.open(out_filename_detect_snv);
	if(!out_file_snv.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << out_filename_detect_snv << endl;
		exit(1);
	}
	out_file_indel.open(out_filename_detect_indel);
	if(!out_file_indel.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << out_filename_detect_indel << endl;
		exit(1);
	}
	out_file_clipReg.open(out_filename_detect_clipReg);
	if(!out_file_clipReg.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << out_filename_detect_clipReg << endl;
		exit(1);
	}

	for(size_t i=0; i<chromeVector.size(); i++){
		chr=chromeVector.at(i);
		if(chr->valid_flag){
			copySingleFile(chr->out_filename_detect_indel, out_file_indel); // indel
			copySingleFile(chr->out_filename_detect_snv, out_file_snv); // snv
			copySingleFile(chr->out_filename_detect_clipReg, out_file_clipReg); // clip regions
		}
	}

	out_file_snv.close();
	out_file_indel.close();
	out_file_clipReg.close();
}

// local consensus for genome
int Genome::genomeLocalCons(){
	Chrome *chr;
	Time time;

	// load consensus data
	genomeLoadDataCons();

	cout << "Number of previously processed regions: " << paras->cns_reg_preDone_num << endl;
	cout << "Number of regions to be processed: " << paras->cns_reg_work_total << endl;

	// invoke the monitor of conseneus work process
	startWorkProcessMonitor(work_finish_filename_cns, paras->monitoring_proc_names_cns, paras->max_proc_running_minutes_cns);

	// begin consensus
	if(!paras->cns_work_vec.empty()) cout << "[" << time.getTime() << "]: start local consensus ..." << endl;
	processConsWork();

	// generate work process finish file
	generateFile(work_finish_filename_cns);

	cout << "[" << time.getTime() << "]: Finalizing consensus ..." << endl;

	computeVarNumStatCons(); // compute statistics for cns command

	if(!paras->cns_work_vec.empty()) destroyConsWorkOptVec(paras->cns_work_vec);

	// reset consensus data
	for(size_t i=0; i<chromeVector.size(); i++){
		chr = chromeVector.at(i);
		if(chr->valid_flag)
			chr->chrResetConsData();
	}

	return 0;
}

void Genome::genomeLoadDataCons(){
	Chrome *chr;

	if(chromeVector.size()==0) return;

	// load X, Y, hs37d5, alt
	for(size_t i=0; i<chromeVector.size(); i++){
		chr = chromeVector.at(i);
		if((chr->decoy_flag==false or paras->include_decoy==false) and (chr->alt_flag==false or paras->include_alt==false) and chr->valid_flag){ // non-decoy, non-alt
			if(chr->chrname.compare(CHR_X_STR1)==0 or chr->chrname.compare(CHR_X_STR2)==0 or chr->chrname.compare(CHR_Y_STR1)==0 or chr->chrname.compare(CHR_Y_STR2)==0){
				chr->chrLoadDataCons();  // load the variant data
				chr->chrGenerateLocalConsWorkOpt();     // generate local consensus work
			}
		}else if(((chr->decoy_flag and paras->include_decoy) or (chr->alt_flag and paras->include_alt)) and chr->valid_flag){ // decoy and alt
			chr->chrLoadDataCons();  // load the variant data
			chr->chrGenerateLocalConsWorkOpt();     // generate local consensus work
		}
	}

	// load other chrs
	for(size_t i=0; i<chromeVector.size(); i++){
		chr = chromeVector.at(i);
		if(chr->decoy_flag==false and chr->alt_flag==false and chr->valid_flag)
		{
			if(chr->chrname.compare(CHR_X_STR1)!=0 and chr->chrname.compare(CHR_X_STR2)!=0 and chr->chrname.compare(CHR_Y_STR1)!=0 and chr->chrname.compare(CHR_Y_STR2)!=0){
				chr->chrLoadDataCons();  // load the variant data
				chr->chrGenerateLocalConsWorkOpt();     // generate local consensus work
			}
		}
	}

	// load previously consensus information
	for(size_t i=0; i<chromeVector.size(); i++){
		chr = chromeVector.at(i);
		if(chr->valid_flag)
			chr->loadPrevConsInfo();
	}
}

// process consensus work using thread pool
int Genome::processConsWork(){
	cnsWork_opt *cns_work_opt;
	cnsWork *cns_work;
	ofstream *var_cand_file;
	size_t num_threads_work, num_work, num_work_percent;
	int64_t sv_len_sum, min_pos, max_pos, dist;
	size_t i, j;

	if(paras->cns_work_vec.empty()) return 0;  // no consensus work, then return directly

	//num_threads_work = (paras->num_threads>=0.5*sysconf(_SC_NPROCESSORS_ONLN)) ? 0.5*sysconf(_SC_NPROCESSORS_ONLN) : paras->num_threads;
	num_threads_work = (paras->num_threads>=sysconf(_SC_NPROCESSORS_ONLN)) ? sysconf(_SC_NPROCESSORS_ONLN) : paras->num_threads;

	if(num_threads_work<1) num_threads_work = 1;

//	if(paras->num_threads_per_cns_work==0)
//		cout << "Local consensus will be processed using " << num_threads_work << " threads, and each work will use the default number of threads" << endl;
//	else
//		cout << "Local consensus will be processed using " << num_threads_work << " concurrent works, and each work will be limited to " << paras->num_threads_per_cns_work << " threads" << endl;

	hts_tpool *p = hts_tpool_init(num_threads_work);
	hts_tpool_process *q = hts_tpool_process_init(p, num_threads_work*2, 1);

	pthread_mutex_init(&paras->mtx_cns_reg_workDone_num, NULL);

	paras->cns_reg_workDone_num = 0;
	num_work = paras->cns_work_vec.size();
	num_work_percent = num_work / (paras->num_parts_progress >> 1);
	if(num_work_percent==0) num_work_percent = 1;
	for(i=0; i<num_work; i++){
		cns_work_opt = paras->cns_work_vec.at(i);
		var_cand_file = getVarcandFile(cns_work_opt->chrname, chromeVector, cns_work_opt->clip_reg_flag);
		if(var_cand_file==NULL){
			cerr << __func__ << ", line=" << __LINE__ << ": cannot get car_cand file for CHR: " << cns_work_opt->chrname << ", error!" << endl;
			exit(1);
		}
		// clipReg_reads_1_26966226-26974816.fq, clipReg_reads_1_21506254-21506762.fq, clipReg_cns_16_209822-210880.fa, clipReg_cns_1_197756788-197757987.fa (INV)
		// cns_chr1_1461692-1461767.fa, cns_chr1_1595082-1595934.fa, cns_chr1_1663381-1663402.fa
		// chr1_hg002_clr_20241218/2_cns/1/cns_1_19387785-19387946.fa
		// chr1_hg002_clr_20241218/2_cns/1/cns_1_19386711-19386935.fa
		// chr1_hg002_clr_20241218/2_cns/1/cns_1_1074689-1076105.fa, cns_1_3097653-3098167.fa, cns_1_670807-676101.fa, cns_1_1223650-1224865.fa, cns_1_7174821-7175786.fa, cns_1_8391251-8391639.fa
		// cns_1_844226-845192.fa, cns_1_46679031-46679316.fa, cns_1_48187816-48187819, cns_1_44605948-44606254.fa, cns_1_30243897-30244578.fa
		// clipReg_cns_chr1_114265459-114272005.fa
//		if(cns_work_opt->contigfilename.compare("debug_tmp/2_cns/X/cns_X_811154-811194.fa")==0){
//			cout << "cnsfilename=" << cns_work_opt->contigfilename << endl;
//		}else continue;

		cns_work = new cnsWork();
		cns_work->cns_work_opt = cns_work_opt;
		cns_work->work_id = i;
		cns_work->num_work = num_work;
		cns_work->num_work_percent = num_work_percent;
		cns_work->p_cns_reg_workDone_num = &(paras->cns_reg_workDone_num);
		cns_work->p_mtx_cns_reg_workDone_num = &(paras->mtx_cns_reg_workDone_num);
		cns_work->num_threads_per_cns_work = paras->num_threads_per_cns_work;
		cns_work->minClipEndSize = paras->minClipEndSize;
		cns_work->maxVarRegSize = paras->maxVarRegSize;
		cns_work->min_seqsim_match = paras->min_seqsim_match;

		sv_len_sum = 0;
		min_pos = INT_MAX, max_pos = 0;
		for(j=0; j<cns_work_opt->arr_size; j++) {
			sv_len_sum += abs(cns_work_opt->var_array[j]->sv_len);
			if(min_pos>cns_work_opt->var_array[j]->startRefPos) min_pos = cns_work_opt->var_array[j]->startRefPos;
			if(max_pos<cns_work_opt->var_array[j]->endRefPos) max_pos = cns_work_opt->var_array[j]->endRefPos;
		}
		dist = max_pos - min_pos + 1;
		if(sv_len_sum<dist) sv_len_sum = dist;

		if(cns_work_opt->clip_reg_flag==false){ // indel
			if(sv_len_sum>paras->cnsChunkSize)
				cns_work->cnsSideExtSize = paras->cnsSideExtSize * CNS_EXT_INDEL_FACTOR_1K + sv_len_sum;
			else if(sv_len_sum>0.5*paras->cnsChunkSize)
				cns_work->cnsSideExtSize = paras->cnsSideExtSize * CNS_EXT_INDEL_FACTOR_500BP + sv_len_sum;
			else
				cns_work->cnsSideExtSize = paras->cnsSideExtSize + sv_len_sum;

//			if(sv_len_sum>2000){
//				cout << "Large indel: arr_size=" << cns_work_opt->arr_size << ", " << cns_work_opt->contigfilename << ", sv_len=" << sv_len_sum << endl;
//			}

		}else{ // clipping
			if(sv_len_sum>1.2*paras->cnsSideExtSizeClip or sv_len_sum==0) // large size or size undetermined
				cns_work->cnsSideExtSize = paras->cnsSideExtSizeClip * CNS_EXT_CLIPREG_FACTOR_6K;
			else if(sv_len_sum>0.8*paras->cnsSideExtSizeClip)
				cns_work->cnsSideExtSize = paras->cnsSideExtSizeClip * CNS_EXT_CLIPREG_FACTOR_4K;
			else if(sv_len_sum>0.4*paras->cnsSideExtSizeClip)
				cns_work->cnsSideExtSize = paras->cnsSideExtSizeClip * CNS_EXT_CLIPREG_FACTOR_2K;
			else if(sv_len_sum>0.2*paras->cnsSideExtSizeClip)
				cns_work->cnsSideExtSize = paras->cnsSideExtSizeClip * CNS_EXT_CLIPREG_FACTOR_1K;
			else
				cns_work->cnsSideExtSize = paras->cnsSideExtSizeClip;
		}
		cns_work->minConReadLen = paras->minConReadLen;
		cns_work->min_sv_size = paras->min_sv_size_usr;
		cns_work->min_supp_num = paras->minReadsNumSupportSV;
		cns_work->sv_len_est = sv_len_sum;
		cns_work->max_seg_size_ratio = paras->max_seg_size_ratio_usr;
		cns_work->max_seg_nm_ratio = paras->max_seg_nm_ratio_usr;
		cns_work->max_absig_density = paras->max_absig_density;
		cns_work->inBamFile = paras->inBamFile;
		cns_work->pg_runid_str = paras->pg_runid_str;
		cns_work->fai = fai;
		cns_work->var_cand_file = var_cand_file;
		cns_work->expected_cov_cns = paras->expected_cov_cns;
		cns_work->min_input_cov_canu = paras->min_input_cov_canu;
		cns_work->max_ultra_high_cov = paras->max_ultra_high_cov;
		cns_work->delete_reads_flag = paras->delete_reads_flag;
		cns_work->keep_failed_reads_flag = paras->keep_failed_reads_flag;
		cns_work->technology = paras->technology;
		//cns_work->canu_version = paras->canu_version;
		cns_work->minMapQ = paras->minMapQ;
		cns_work->minHighMapQ = paras->minHighMapQ;

		hts_tpool_dispatch(p, q, processSingleConsWork, cns_work);
	}

	hts_tpool_process_flush(q);
	hts_tpool_process_destroy(q);
	hts_tpool_destroy(p);

	// create the work finish file
	//generateFile(work_finish_filename);

	return 0;
}

// get car_cand output file according to given 'chrname'
ofstream* Genome::getVarcandFile(string &chrname, vector<Chrome*> &chrome_vec, bool clip_reg_flag){
	ofstream *var_cand_file = NULL;
	Chrome *chr;
	for(size_t i=0; i<chrome_vec.size(); i++){
		chr = chrome_vec.at(i);
		if(chr->chrname.compare(chrname)==0){
			if(clip_reg_flag) var_cand_file = &chr->var_cand_clipReg_file;
			else var_cand_file = &chr->var_cand_indel_file;
			break;
		}
	}
	return var_cand_file;
}

// generate file
void Genome::generateFile(string &filename){
	// create the work finish file
	ofstream out_file(filename);
	if(!out_file.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << filename << endl;
		exit(1);
	}
	out_file << "Work finished" << endl;
	out_file.close();
}

// call variants for genome
int Genome::genomeCall(){
	Time time;
	
	// load clipping region information
	genomeLoadMateClipRegData();

	// collect regions
	cout << "Loading regions ..." << endl;
	genomeCollectCallWork();

	cout << "Number of regions to be processed: " << paras->call_work_num << endl;

	// invoke the monitor of consensus work process
	startWorkProcessMonitor(work_finish_filename_call, paras->monitoring_proc_names_call, paras->max_proc_running_minutes_call);

	// alignment work
	time.setStartTime();
	processAlnWork();
	time.printElapsedTime();

	// process call work
	time.setStartTime();
	processCallWork();
	time.printElapsedTime();

	// finish call work, deleted
	//genomeFinishCallWork();

	// generate work process finish file
	generateFile(work_finish_filename_call);

	// call TRA according to mate clip regions
	//genomeCallTra(); // 2023-12-25
	//genomeCallTra2(); // 2025-09-12

	// fill variant sequence
	//cout << "4444444444444" << endl;
//	cout << "[" << time.getTime() << "]: fill variant sequences... " << endl;
//	genomeFillVarseq();

	// phasing
	if(paras->phasing_flag){
		cout << "[" << time.getTime() << "]: local phasing ..." << endl;
		genomePhasing();
	}

	// save SV to file
	//cout << "5555555555555" << endl;
	cout << "[" << time.getTime() << "]: save call result to file... " << endl;
	genomeSaveCallSV2File();

	// merge call results into single file
	//cout << "6666666666666" << endl;
	cout << "[" << time.getTime() << "]: merge call result... " << endl;
	mergeCallResult();

	// sort variant items in BED format
	//cout << "[" << time.getTime() << "]: sortVarResults 1... " << endl;
	sortVarResults(out_filename_result_vars, 0);

	// save results in VCF file format
	//cout << "[" << time.getTime() << "]: saveResultVCF... " << endl;
	saveResultVCF();

	// sort variant items in VCF format
	//cout << "[" << time.getTime() << "]: sortVarResults 2... " << endl;
	sortVarResults(out_filename_result_vars_vcf, 1);

	// compute statistics for call command
	cout << "[" << time.getTime() << "]: compute variant NUMBER statistics... " << endl;
	computeVarNumStatCall();

	// remove temporary monitor files
	removeTempMonitorFiles();

	return 0;
}

// collect call work
void Genome::genomeCollectCallWork(){
	Chrome *chr;

	if(chromeVector.size()==0) return;

	// load X, Y, hs37d5, alt
	for(size_t i=0; i<chromeVector.size(); i++){
		chr = chromeVector.at(i);
		if((chr->decoy_flag==false or paras->include_decoy==false) and (chr->alt_flag==false or paras->include_alt==false) and chr->valid_flag){ // non-decoy, non-alt
			if(chr->chrname.compare(CHR_X_STR1)==0 or chr->chrname.compare(CHR_X_STR2)==0 or chr->chrname.compare(CHR_Y_STR1)==0 or chr->chrname.compare(CHR_Y_STR2)==0){
				chr->chrLoadDataCall();
				chr->chrCollectCallWork();
			}
		}else if(((chr->decoy_flag and paras->include_decoy) or (chr->alt_flag and paras->include_alt)) and chr->valid_flag){ // decoy or alt
			chr->chrLoadDataCall();
			chr->chrCollectCallWork();
		}
	}

	// load other chrs
	for(size_t i=0; i<chromeVector.size(); i++){
		chr = chromeVector.at(i);
		if(chr->decoy_flag==false and chr->alt_flag==false and chr->valid_flag){ // non-decoy and no alt
			if(chr->chrname.compare(CHR_X_STR1)!=0 and chr->chrname.compare(CHR_X_STR2)!=0 and chr->chrname.compare(CHR_Y_STR1)!=0 and chr->chrname.compare(CHR_Y_STR2)!=0){
				chr->chrLoadDataCall();
				chr->chrCollectCallWork();
			}
		}
	}
}

// call variants using thread pool
int Genome::processAlnWork(){
	callWork_opt *call_work_opt;
	varCand *var_cand;
	size_t num_work, num_work_percent;

	if(paras->call_work_vec.empty()) return 0;  // no align work, then return directly

	cout << "Begin sequence alignment ..." << endl;

	hts_tpool *p = hts_tpool_init(paras->num_threads);
	hts_tpool_process *q = hts_tpool_process_init(p, paras->num_threads*2, 1);

	pthread_mutex_init(&paras->mtx_call_workDone_num, NULL);

	paras->call_workDone_num = 0;
	num_work = paras->call_work_vec.size();
	//num_work_percent = num_work / (paras->num_parts_progress >> 1);
	num_work_percent = num_work / (paras->num_parts_progress / 10);
	if(num_work_percent==0) num_work_percent = 1;
	for(size_t i=0; i<num_work; i++){
		var_cand = paras->call_work_vec.at(i);
		call_work_opt = new callWork_opt();
		call_work_opt->var_cand = var_cand;
		call_work_opt->work_id = i;
		call_work_opt->num_work = num_work;
		call_work_opt->num_work_percent = num_work_percent;
		call_work_opt->p_call_workDone_num = &(paras->call_workDone_num);
		call_work_opt->p_mtx_call_workDone_num = &(paras->mtx_call_workDone_num);

		hts_tpool_dispatch(p, q, processSingleMinimap2AlnWork, call_work_opt);
	}

	hts_tpool_process_flush(q);
	hts_tpool_process_destroy(q);
	hts_tpool_destroy(p);

	return 0;
}

// call variants using thread pool
int Genome::processCallWork(){
	callWork_opt *call_work_opt;
	varCand *var_cand;
	size_t num_work, num_work_percent;

	if(paras->call_work_vec.empty()) return 0;  // no call work, then return directly

	cout << "Begin variants call ..." << endl;

	hts_tpool *p = hts_tpool_init(paras->num_threads);
	hts_tpool_process *q = hts_tpool_process_init(p, paras->num_threads*2, 1);

	pthread_mutex_init(&paras->mtx_call_workDone_num, NULL);

	paras->call_workDone_num = 0;
	num_work = paras->call_work_vec.size();
	//num_work_percent = num_work / (paras->num_parts_progress >> 1); // deleted on 2024-09-04
	num_work_percent = num_work / (paras->num_parts_progress / 10);
	if(num_work_percent==0) num_work_percent = 1;
	for(size_t i=0; i<num_work; i++){
		var_cand = paras->call_work_vec.at(i);

		// DUP not precise (CCS30x): blat_1_2936746-2942685.sim4, blat_contig_1_1180102-1180675.sim4, blat_contig_1_1183812-1185067.sim4, blat_contig_1_1317611-1318285.sim4, blat_1_1860801-1869285.sim4
		// DUP check: blat_1_802003-808928.sim4, blat_1_843249-844905.sim4
		// diploid: blat_1_4480337-4489601.sim4, blat_1_5364079-5371326.sim4, blat_1_5727691-5736300.sim4 (good), blat_1_7613307-7613853.sim4
		// blat_1_19156546-19164246.sim4, blat_1_2415202-2415425.sim4, tra_blat_1_2686251-2691837.sim4
		// blat_contig_chr1_253768-256236.sim4, blat_contig_chr1_1772083-1775128.sim4, blat_contig_chr1_1772083-1775128.sim4, blat_contig_chr1_2068132-2073121.sim4
		// minimap2_1_2213173-2214000.paf, minimap2_contig_1_26966226-26974816.paf, minimap2_contig_1_21506254-21506762.paf, minimap2_contig_1_29382165-29382465.paf
		// minimap2_cns_5_1414495-1414633.paf, minimap2_cns_5_1192021-1192323.paf, minimap2_cns_1_26966226-26974816.paf, minimap2_chr1_143246329-143247038.paf, minimap2_1_649705-650290.paf
		// minimap2_1_2765904-2766241.paf, minimap2_1_5447230-5447236.paf, minimap2_1_3097653-3098167.paf, minimap2_1_2618216-2619007.paf, minimap2_1_1184268-1184983.paf, minimap2_1_5446916-5447358
		// minimap2_1_3560147-3560992.paf, minimap2_1_801940-802476.paf, minimap2_1_3560147-3560992.paf, minimap2_1_3561166-3562123.paf
//		if(var_cand->alnfilename.compare("debug_tmp/3_call/1/minimap2_cns_1_110663085-110663198.paf")==0){
//			cout << var_cand->alnfilename << endl;
//		}else continue;

		call_work_opt = new callWork_opt();
		call_work_opt->var_cand = var_cand;
		call_work_opt->work_id = i;
		call_work_opt->num_work = num_work;
		call_work_opt->num_work_percent = num_work_percent;
		call_work_opt->p_call_workDone_num = &(paras->call_workDone_num);
		call_work_opt->p_mtx_call_workDone_num = &(paras->mtx_call_workDone_num);

		hts_tpool_dispatch(p, q, processSingleCallWork, call_work_opt);
	}

	hts_tpool_process_flush(q);
	hts_tpool_process_destroy(q);
	hts_tpool_destroy(p);

	return 0;
}

// load clipping region information
void Genome::genomeLoadMateClipRegData(){
	Chrome *chr;
	for(size_t i=0; i<chromeVector.size(); i++){
		chr = chromeVector.at(i);
		if(chr->valid_flag)
			chr->chrLoadMateClipRegData();
	}
}

// phasing
void Genome::genomePhasing(){
	Chrome *chr;
	for(size_t i=0; i<chromeVector.size(); i++){
		chr = chromeVector.at(i);
		if(chr->process_block_num and chr->valid_flag)
			chr->chrPhasing();
	}
}

// save variants to file
void Genome::genomeSaveCallSV2File(){
	mkdir(out_dir_result.c_str(), S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);  // create the directory for final results

	Chrome *chr;
	for(size_t i=0; i<chromeVector.size(); i++){
		chr = chromeVector.at(i);
		if(chr->valid_flag)
			chr->saveCallSV2File();
	}

	// save TRA
	saveTraCall2File();
}

// save TRA
void Genome::saveTraCall2File(){
	size_t i, j, m;
	Chrome *chr;
	ofstream outfile_tra;
	vector<mateClipReg_t*> mate_clipReg_vec;
	mateClipReg_t *mate_clip_reg;
	string line, chrname_tmp, header_line_bedpe, mate_reg_str, id_str, sv_type_str, sv_type_str2;
	string gt_seq, gt_header, gt_str, ad_str1, ad_str2, dp_str, id_str_mate, id_vec_mate_str, id_vec_mate_mate_str, sep_str, mate_item_str;
	int32_t k, reg_id, DP, AD, id_vec_mate, clip_end, clip_end2, supp_num, DP_num, supp_num_mate, id_vec_mate_mate;
	reg_t *reg1, *reg2;
	int64_t left_clip_pos, right_clip_pos;
	double AF;
	vector<string> str_vec, str_vec2, str_vec_tmp, str_vec2_tmp, mate_str_vec2_tmp, bnd_mate_reg_strs_new;

	outfile_tra.open(out_filename_result_tra);
	if(!outfile_tra.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << out_filename_result_tra << endl;
		exit(1);
	}
	
	header_line_bedpe = getCallFileHeaderBedpe(paras->sample);
	outfile_tra << header_line_bedpe << endl;

	//num_bnd = num_tra = num_dup = num_inv = num_del = 1;
	for(i=0; i<chromeVector.size(); i++){
		chr = chromeVector.at(i);
		if(chr->valid_flag){
			mate_clipReg_vec = chr->mateClipRegVector;
			for(j=0; j<mate_clipReg_vec.size(); j++){
				mate_clip_reg = mate_clipReg_vec.at(j);
				if(mate_clip_reg->valid_flag){
					if(mate_clip_reg->ultra_large_dist_flag==false and (mate_clip_reg->sv_type==VAR_TRA or mate_clip_reg->sv_type==VAR_BND)){ // valid TRA
						// only process half part
						for(reg_id=0; reg_id<2; reg_id++){
							if(mate_clip_reg->bnd_mate_reg_strs[reg_id].compare("-")!=0){
								str_vec = split(mate_clip_reg->bnd_mate_reg_strs[reg_id], ",");
								for(m=0; m<str_vec.size(); m++){
									if(str_vec.at(m).compare("-")!=0){
										reg1 = reg2 = NULL;
										left_clip_pos = right_clip_pos = -1;
										supp_num = supp_num_mate = 0;

										// main part
										str_vec2_tmp = split(str_vec.at(m), "_");
										for(const auto& item : str_vec2_tmp){
											str_vec2 = split(item, "|");
											id_str_mate = str_vec2.at(0);
											id_vec_mate_str = id_str_mate.at(0);
											id_vec_mate = stoi(id_vec_mate_str);
											if(m==0){
												clip_end = LEFT_END;
												if(id_str_mate.at(1)=='+') clip_end2 = RIGHT_END;
												else clip_end2 = LEFT_END;
											}else{
												clip_end = RIGHT_END;
												if(id_str_mate.at(1)=='+') clip_end2 = LEFT_END;
												else clip_end2 = RIGHT_END;
											}
											left_clip_pos = stoi(str_vec2.at(1));
											supp_num = stoi(str_vec2.at(2));
											DP_num = stoi(str_vec2.at(3));

											// deleted on 2026-02-23
		//									if(clip_end==LEFT_END){
		//										if(mate_clip_reg->leftClipReg) reg1 = mate_clip_reg->leftClipReg;
		//									}else{
		//										if(mate_clip_reg->leftClipReg2) reg1 = mate_clip_reg->leftClipReg2;
		//									}
											if(mate_clip_reg->leftClipReg) reg1 = mate_clip_reg->leftClipReg;
											else if(mate_clip_reg->leftClipReg2) reg1 = mate_clip_reg->leftClipReg2;

											// mate part
											if(id_vec_mate==2 or id_vec_mate==3){
												id_vec_mate_mate = right_clip_pos = supp_num_mate = -1;
												mate_str_vec2_tmp.clear();
												mate_item_str = "";
												if(mate_clip_reg->bnd_mate_reg_strs[id_vec_mate].compare("-")!=0){
													str_vec_tmp = split(mate_clip_reg->bnd_mate_reg_strs[id_vec_mate], ",");
													if(clip_end2==LEFT_END){ // mate at left end
														if(str_vec_tmp.at(0).compare("-")!=0) mate_str_vec2_tmp = split(str_vec_tmp.at(0), "_");
													}else{ // mate at right end
														if(str_vec_tmp.at(1).compare("-")!=0) mate_str_vec2_tmp = split(str_vec_tmp.at(1), "_");
													}
													for(const auto& item2 : mate_str_vec2_tmp){
														str_vec2 = split(item2, "|");
														id_vec_mate_mate_str = str_vec2.at(0).at(0);
														id_vec_mate_mate = stoi(id_vec_mate_mate_str);
														right_clip_pos = stoi(str_vec2.at(1));
														supp_num_mate = stoi(str_vec2.at(2));
														//DP_num_mate = stoi(str_vec2.at(3));
														if(id_vec_mate_mate==reg_id) {
															mate_item_str = item2;
															break;
														}
													}
													if(mate_item_str.size()==0) mate_item_str = "-";

													if(id_vec_mate==2){
														if(mate_clip_reg->rightClipReg) reg2 = mate_clip_reg->rightClipReg;
													}else{
														if(mate_clip_reg->rightClipReg2) reg2 = mate_clip_reg->rightClipReg2;
													}
												}
											}

											if(reg1 and reg2 and left_clip_pos!=-1 and right_clip_pos!=-1 and supp_num>=paras->minReadsNumSupportSV and supp_num_mate>=paras->minReadsNumSupportSV){
												//cout << "sv_type=" << mate_clip_reg->sv_type << ": " << reg1->chrname << ":" << left_clip_pos << " <--> " << reg2->chrname << ":" << right_clip_pos << endl;

												if(mate_clip_reg->sv_type==VAR_TRA){
													sv_type_str = VAR_TRA_STR;
													paras->num_tra ++;
													id_str = "ASVCLR." + sv_type_str + "." + to_string(paras->num_tra);
												}else{
													sv_type_str = VAR_BND_STR;
													paras->num_bnd ++;
													id_str = "ASVCLR." + sv_type_str + "." + to_string(paras->num_bnd);
												}
												sv_type_str2 = "<" + sv_type_str + ">";

												line = reg1->chrname + "\t" + to_string(left_clip_pos) + "\t" + to_string(left_clip_pos) + "\t" + reg2->chrname + "\t" + to_string(right_clip_pos) + "\t" + to_string(right_clip_pos); // chr1, start1, end1, chr2, start2, end2
												line +=  + "\t" + id_str + "\t" + sv_type_str2 + "\t-\t-"; // ID, sv_type, SVLEN1, SVLEN2

												// update BND mate regions
												mate_reg_str = "";
												for(k=0; k<4; k++){
													if(k==0) sep_str = "";
													else sep_str = ";";
													if(k==reg_id){
														if(clip_end==LEFT_END) mate_reg_str += sep_str + item + ",-";
														else mate_reg_str += sep_str + "-," + item;
													}else if(k==id_vec_mate){
														if(clip_end2==LEFT_END) mate_reg_str += sep_str + mate_item_str + ",-";
														else mate_reg_str += sep_str + "-," + mate_item_str;
													}else{
														mate_reg_str += sep_str + "-";
													}
												}

	//											mate_reg_str = mate_clip_reg->bnd_mate_reg_strs[0];
	//											for(k=1; k<4; k++) mate_reg_str += ";" + mate_clip_reg->bnd_mate_reg_strs[k];

												line += "\t" + mate_reg_str;
												line += "\t-\t-\t-\t-"; // Ref1, Alt1, Ref2, Alt2

												// extra information, split with ";": DP, AF, LDISCOV
												AD = supp_num;
												DP = DP_num;
												AF = (double)AD / DP;
												if(AF<0){
													cout << "\t---------- AF=" << AF << ", sv_type=" << mate_clip_reg->sv_type << ": " << reg1->chrname << ":" << left_clip_pos << " <--> " << reg2->chrname << ":" << right_clip_pos << ", error!" << endl;
													exit(1);
												}
												stringstream ss;
												ss << setprecision(3) << AF;
												line += "\tDP=" + to_string(DP) + ";AF=" + ss.str() + ";LDISCOV=READS";

												gt_header = "GT:AD:DP";
												gt_str = "./.";
												ad_str1 = ".";
												ad_str2 = ".";
												dp_str = ".";

												if(AF!=-1){
													if(AF>=paras->gt_homo_ratio){
														gt_str = GT_HOMOZYGOUS_STR;
													}else if(AF>=paras->gt_hete_ratio){
														gt_str = GT_HETEROZYGOUS_STR;
													}else{
														gt_str = GT_NOZYGOUS_STR;
													}
													ad_str1 = to_string(DP-AD);
													ad_str2 = to_string(AD);
													dp_str = to_string(DP);
												}
												gt_seq = gt_header + "\t" + gt_str + ":" + ad_str1 + "," + ad_str2 + ":" + dp_str;
												line += "\t" + gt_seq;

												outfile_tra << line << endl;
											}
										}
									}
								}
							}
						}
					}else if(mate_clip_reg->ultra_large_dist_flag){ // process ultra-large items
						//printMateClipReg(mate_clip_reg);

						for(reg_id=0; reg_id<2; reg_id++){
							if(mate_clip_reg->bnd_mate_reg_strs[reg_id].compare("-")!=0){
								str_vec = split(mate_clip_reg->bnd_mate_reg_strs[reg_id], ",");
								for(m=0; m<str_vec.size(); m++){
									if(str_vec.at(m).compare("-")!=0){
										reg1 = reg2 = NULL;
										left_clip_pos = right_clip_pos = -1;
										supp_num = supp_num_mate = 0;

										str_vec2_tmp = split(str_vec.at(m), "_");
										for(const auto& item : str_vec2_tmp){
											str_vec2 = split(item, "|");
											id_str_mate = str_vec2.at(0);
											id_vec_mate_str = id_str_mate.at(0);
											id_vec_mate = stoi(id_vec_mate_str);

											if(m==0){
												//clip_end = LEFT_END;
												if(id_str_mate.at(1)=='+') clip_end2 = RIGHT_END;
												else clip_end2 = LEFT_END;
											}else{
												//clip_end = RIGHT_END;
												if(id_str_mate.at(1)=='+') clip_end2 = LEFT_END;
												else clip_end2 = RIGHT_END;
											}

											left_clip_pos = stoi(str_vec2.at(1));
											supp_num = stoi(str_vec2.at(2));
											DP_num = stoi(str_vec2.at(3));

											// added on 2026-01-23
											if(mate_clip_reg->leftClipReg) reg1 = mate_clip_reg->leftClipReg;
											else if(mate_clip_reg->leftClipReg2) reg1 = mate_clip_reg->leftClipReg2;

											if(id_vec_mate==2 or id_vec_mate==3){
												id_vec_mate_mate = right_clip_pos = supp_num_mate = -1;
												mate_str_vec2_tmp.clear();
												mate_item_str = "";
												if(mate_clip_reg->bnd_mate_reg_strs[id_vec_mate].compare("-")!=0){
													str_vec_tmp = split(mate_clip_reg->bnd_mate_reg_strs[id_vec_mate], ",");
													if(clip_end2==LEFT_END){
														if(str_vec_tmp.at(0).compare("-")!=0) mate_str_vec2_tmp = split(str_vec_tmp.at(0), "_");
													}else{
														if(str_vec_tmp.at(1).compare("-")!=0) mate_str_vec2_tmp = split(str_vec_tmp.at(1), "_");
													}
													for(const auto& item2 : mate_str_vec2_tmp){
														str_vec2 = split(item2, "|");
														id_vec_mate_mate_str = str_vec2.at(0).at(0);
														id_vec_mate_mate = stoi(id_vec_mate_mate_str);
														right_clip_pos = stoi(str_vec2.at(1));
														supp_num_mate = stoi(str_vec2.at(2));
														if(id_vec_mate_mate==reg_id){
															mate_item_str = item2;
															break;
														}
													}
													if(mate_item_str.size()==0) mate_item_str = "-";

													if(id_vec_mate==2){
														if(mate_clip_reg->rightClipReg) reg2 = mate_clip_reg->rightClipReg;
													}else{
														if(mate_clip_reg->rightClipReg2) reg2 = mate_clip_reg->rightClipReg2;
													}
												}
											}

		//									if(reg1==NULL or reg2==NULL){
		//										printMateClipReg(mate_clip_reg);
		//									}

											if(reg1 and reg2 and left_clip_pos!=-1 and right_clip_pos!=-1 and supp_num>=paras->minReadsNumSupportSV and supp_num_mate>=paras->minReadsNumSupportSV){
												switch(mate_clip_reg->sv_type){
													case VAR_DUP:
														sv_type_str = VAR_DUP_STR; // "DUP"
														paras->num_dup ++;
														id_str = "ASVCLR." + sv_type_str + "." + to_string(paras->num_dup);
														break;
													case VAR_INV:
														sv_type_str = VAR_INV_STR; // "INV"
														paras->num_inv ++;
														id_str = "ASVCLR." + sv_type_str + "." + to_string(paras->num_inv);
														break;
													case VAR_DEL:
														sv_type_str = VAR_DEL_STR; // "DEL"
														paras->num_del ++;
														id_str = "ASVCLR." + sv_type_str + "." + to_string(paras->num_del);
														break;
													case VAR_TRA:
														sv_type_str = VAR_TRA_STR;
														paras->num_tra ++;
														id_str = "ASVCLR." + sv_type_str + "." + to_string(paras->num_tra);
														break;
													case VAR_BND:
														sv_type_str = VAR_BND_STR;
														paras->num_bnd ++;
														id_str = "ASVCLR." + sv_type_str + "." + to_string(paras->num_bnd);
														break;
													default:
														cout << "line=" << __LINE__ << ", sv_type=" << mate_clip_reg->sv_type << ", error!" << endl;
														break;
												}

												sv_type_str2 = "<" + sv_type_str + ">";
												line = reg1->chrname + "\t" + to_string(left_clip_pos) + "\t" + to_string(left_clip_pos) + "\t" + reg2->chrname + "\t" + to_string(right_clip_pos) + "\t" + to_string(right_clip_pos);
												line += "\t" + id_str + "\t" + sv_type_str2 + "\t-\t-";

												// update BND mate regions
												mate_reg_str = "";
												for(k=0; k<4; k++){
													if(k==0) sep_str = "";
													else sep_str = ";";
													if(k==reg_id){
														if(clip_end==LEFT_END) mate_reg_str += sep_str + item + ",-";
														else mate_reg_str += sep_str + "-," + item;
													}else if(k==id_vec_mate){
														if(clip_end2==LEFT_END) mate_reg_str += sep_str + mate_item_str + ",-";
														else mate_reg_str += sep_str + "-," + mate_item_str;
													}else{
														mate_reg_str += sep_str + "-";
													}
												}

	//											mate_reg_str = mate_clip_reg->bnd_mate_reg_strs[0];
	//											for(k=1; k<4; k++) mate_reg_str += ";" + mate_clip_reg->bnd_mate_reg_strs[k];

												line += "\t" + mate_reg_str;
												line += "\t-\t-\t-\t-";

												AD = supp_num;
												DP = DP_num;
												AF = (double)AD / DP;

												stringstream ss;
												ss << setprecision(3) << AF;
												line += "\tDP=" + to_string(DP) + ";AF=" + ss.str() +";ULDist=" + to_string(mate_clip_reg->dist_breakpoint) + ";LDISCOV=READS";

												gt_header = "GT:AD:DP";
												gt_str = "./.";
												ad_str1 = ".";
												ad_str2 = ".";
												dp_str = ".";

												if(AF!=-1){
													if(AF>=paras->gt_homo_ratio){
														gt_str = GT_HOMOZYGOUS_STR;
													}else if(AF>=paras->gt_hete_ratio){
														gt_str = GT_HETEROZYGOUS_STR;
													}else{
														gt_str = GT_NOZYGOUS_STR;
													}
													ad_str1 = to_string(DP-AD);
													ad_str2 = to_string(AD);
													dp_str = to_string(DP);
												}

												gt_seq = gt_header + "\t" + gt_str + ":" + ad_str1 + "," + ad_str2 + ":" + dp_str;
												line += "\t" + gt_seq;

												outfile_tra << line << endl;
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}

	outfile_tra.close();
}

// merge call results into single file
void Genome::mergeCallResult(){
	ofstream out_file_indel, out_file_clipReg, out_file_vars;
	string header_line_bed, header_line_bedpe;
	Chrome *chr;

	out_file_indel.open(out_filename_result_indel);
	if(!out_file_indel.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << out_filename_result_indel << endl;
		exit(1);
	}
	out_file_clipReg.open(out_filename_result_clipReg);
	if(!out_file_indel.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << out_filename_result_clipReg << endl;
		exit(1);
	}
	out_file_vars.open(out_filename_result_vars);
	if(!out_file_vars.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << out_filename_result_vars << endl;
		exit(1);
	}

	// header line
	header_line_bed = getCallFileHeaderBed(paras->sample);
	out_file_indel << header_line_bed << endl;
	out_file_clipReg << header_line_bed << endl;

	for(size_t i=0; i<chromeVector.size(); i++){
		chr = chromeVector.at(i);
		if(chr->valid_flag){
			copySingleFile(chr->out_filename_call_indel, out_file_indel); // indel
			copySingleFile(chr->out_filename_call_clipReg, out_file_clipReg); // clip_reg
		}
	}
	out_file_indel.close();
	out_file_clipReg.close();

	// merge indels, clipping variants, translocations to single file
	header_line_bedpe = getCallFileHeaderBedpe(paras->sample);
	out_file_vars << header_line_bed << endl;
	out_file_vars << header_line_bedpe << endl;
	copySingleFile(out_filename_result_indel, out_file_vars); // indels
	copySingleFile(out_filename_result_clipReg, out_file_vars); // INV, DUP
	copySingleFile(out_filename_result_tra, out_file_vars); // TRA
	out_file_vars.close();
}

// save results in VCF file format for detect command
void Genome::saveResultVCF(){
	// save results
	cout << "Output results to VCF file ..." << endl;
	saveResultToVCF(out_filename_result_vars, out_filename_result_vars_vcf);
}

// save indel in VCF file format for detect command
void Genome::saveResultToVCF(string &in, string &out_vcf){
	string line, line_vcf, chr, start_pos, end_pos, id, ref, alt, qual, filter, info, format, format_val, sample, reg;
	string sv_type, sv_type_tmp, sv_len, ULDist, dup_num, extra_info, seq_str;
	ifstream infile;
	ofstream outfile;
	vector<string> str_vec, SVLEN;

	string bnd_id_str, mate_bnd_id_str, mate_reg_str, mate_bnd_str, mate_bnd_str2, reg_str, mate_chr, chrname1, chrname2, gt_str;
	vector<string> bnd_str_vec, bnd_str_vec2;
	int64_t tra_pos_arr[4];
	int32_t i, j, checked_arr[4][2], seq_len; //, ins_id, del_id, dup_id, inv_id, tra_id;
	vector<BND_t*> bnd_vec, sub_bnd_vec;
	BND_t *bnd_item;
	string format_name_bnd;
	char *seq;
	size_t start_str, end_str, base, support_num, DP;
	double AF;
	// open files
	infile.open(in);
	if(!infile.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << in << endl;
		exit(1);
	}

	outfile.open(out_vcf);
	if(!outfile.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << out_vcf << endl;
		exit(1);
	}

	// save VCF header
	saveVCFHeader(outfile, paras->sample);

	format_name_bnd = "GT:AD:DP";
	// save results
	while(getline(infile, line))
		if(line.size()>0 and line.at(0)!='#'){	// skip header
			str_vec = split(line, "\t");
			if(str_vec.size()<=MAX_BED_COLS_NUM){ // INS, DEL, DUP, INV, MIX, UNC
				chr = str_vec.at(0);		// CHROM
				start_pos = str_vec.at(1);	// POS
				id = str_vec.at(3);			// ID
				ref = str_vec.at(7);		// REF
				alt = str_vec.at(8);		// ALT
				extra_info = str_vec.at(9);	// ExtraInfo
				qual = ".";					// QUAL
				filter = "PASS";			// FILTER

				if(ref.compare("-")==0) ref = ".";
				if(alt.compare("-")==0) alt = ".";

				end_pos = str_vec.at(2);
				sv_type_tmp = str_vec.at(4);
				if(sv_type_tmp.at(0)=='<') {
					sv_type = sv_type_tmp.substr(1, sv_type_tmp.size()-2);
					alt = sv_type_tmp;
				}else sv_type = sv_type_tmp;
				sv_len = str_vec.at(5);
				if(sv_type.compare(VAR_UNC_STR)==0 or sv_type.compare(VAR_UNC_STR1)==0 or sv_type.compare(VAR_UNC_STR2)==0) id = ref = alt = "."; // UNC
				info = "SVTYPE=" + sv_type + ";" + "SVLEN=" + sv_len + ";END=" + end_pos;	// INFO: SVTYPE=sv_type;SVLEN=sv_len;END=end
				if(sv_type.compare(VAR_DUP_STR)==0 or sv_type.compare(VAR_DUP_STR1)==0) { // DUPNUM=dup_num
					dup_num = str_vec.at(6);
					info += ";DUPNUM=" + dup_num;
				}
				info += ";" + extra_info; // DP=xx;AF=yy

				format = str_vec.at(10);			// FORMAT
				format_val = str_vec.at(11);		// gt_value

				line_vcf = chr + "\t" + start_pos + "\t" + id + "\t" + ref + "\t" + alt + "\t" + qual + "\t" + filter + "\t" + info + "\t" + format + "\t" + format_val;
				outfile << line_vcf << endl;

			}else{	// TRA, MIX
				sv_type_tmp = str_vec.at(7);
				if(sv_type_tmp.at(0)=='<') {
					sv_type = sv_type_tmp.substr(1, sv_type_tmp.size()-2);
					alt = sv_type_tmp;
				}else sv_type = sv_type_tmp;

				//if(str_vec.at(7).compare(VAR_TRA_STR)==0){ // TRA, then convert to BND, deleted on 2025-09-19
				if(sv_type.compare(VAR_TRA_STR)==0){ // TRA, then convert to BND

					for(i=0; i<4; i++) { checked_arr[i][0] = checked_arr[i][1] = 0; } // initialize

					chrname1 = str_vec.at(0);		// CHR1
					chrname2 = str_vec.at(3);		// CHR2

					if(str_vec.at(1).compare("-")!=0) tra_pos_arr[0] = stoi(str_vec.at(1));
					else tra_pos_arr[0] = -1;
					if(str_vec.at(2).compare("-")!=0) tra_pos_arr[1] = stoi(str_vec.at(2));
					else tra_pos_arr[1] = -1;
					if(str_vec.at(4).compare("-")!=0) tra_pos_arr[2] = stoi(str_vec.at(4));
					else tra_pos_arr[2] = -1;
					if(str_vec.at(5).compare("-")!=0) tra_pos_arr[3] = stoi(str_vec.at(5));
					else tra_pos_arr[3] = -1;

					bnd_str_vec = split(str_vec.at(10), ";"); // "2+,3-"

					// four regions
					for(i=0; i<4; i++){
						if(bnd_str_vec.at(i).compare("-")!=0){
						//if(bnd_str_vec.at(i).compare("-")!=0 and tra_pos_arr[i]>0){
							// right end
							sub_bnd_vec = generateBNDItems(i, RIGHT_END, checked_arr, chrname1, chrname2, tra_pos_arr, bnd_str_vec, fai, paras->gt_homo_ratio, paras->gt_hete_ratio);
							for(j=0; j<(int32_t)sub_bnd_vec.size(); j++) {
								bnd_item = sub_bnd_vec.at(j);
								//cout << "i=" << i << ", end=" << RIGHT_END << ", reg_id=" << to_string(bnd_item->reg_id) << ", mate_reg_id=" << to_string(bnd_item->mate_reg_id) << endl;
								if(bnd_item->support_num>=paras->minReadsNumSupportSV and ((bnd_item->reg_id<=1 and bnd_item->mate_reg_id>=2) or (bnd_item->reg_id>=2 and bnd_item->mate_reg_id<=1))){
									outfile << bnd_item->vcf_line << endl;
									//cout << bnd_item->vcf_line << endl;
									bnd_vec.push_back(bnd_item);
								}
							}

							// left end
							sub_bnd_vec = generateBNDItems(i, LEFT_END, checked_arr, chrname1, chrname2, tra_pos_arr, bnd_str_vec, fai, paras->gt_homo_ratio, paras->gt_hete_ratio);
							for(j=0; j<(int32_t)sub_bnd_vec.size(); j++) {
								bnd_item = sub_bnd_vec.at(j);
								//cout << "i=" << i << ", end=" << LEFT_END << ", reg_id=" << to_string(bnd_item->reg_id) << ", mate_reg_id=" << to_string(bnd_item->mate_reg_id) << endl;
								if(bnd_item->support_num>=paras->minReadsNumSupportSV and ((bnd_item->reg_id<=1 and bnd_item->mate_reg_id>=2) or (bnd_item->reg_id>=2 and bnd_item->mate_reg_id<=1))){
									outfile << bnd_item->vcf_line << endl;
									//cout << bnd_item->vcf_line << endl;
									bnd_vec.push_back(bnd_item);
								}
							}
						}
					}
				}else{ // MIX
					//cout << __func__ << ": unexpected variant item: " << line << ", skipped." << endl;
					sv_type_tmp = str_vec.at(7);
					if(sv_type_tmp.at(0)=='<') {
						sv_type = sv_type_tmp.substr(1, sv_type_tmp.size()-2);
						alt = sv_type_tmp;
					}else sv_type = sv_type_tmp;
					if(sv_type == VAR_INV_STR || sv_type == VAR_DUP_STR){
						chr = str_vec.at(0);
						id = str_vec.at(6);
						start_str = (str_vec.at(10)).find('|');
						end_str = (str_vec.at(10)).find('|', str_vec.at(10).find(';'));
						start_pos = str_vec.at(10).substr(start_str + 1, str_vec.at(10).find('|', start_str + 1) - (start_str + 1));
						end_pos = str_vec.at(10).substr(end_str + 1, str_vec.at(10).find('|', end_str + 1)   - (end_str + 1));
						base = (stoll(start_pos) < stol(end_pos)) ? start_str : end_str;
						if(stoll(start_pos)>stoll(end_pos)){
							start_pos = str_vec.at(10).substr(end_str + 1, str_vec.at(10).find('|', end_str + 1)   - (end_str + 1));
							end_pos = str_vec.at(10).substr(start_str + 1, str_vec.at(10).find('|', start_str + 1) - (start_str + 1));
						}
						support_num = str_vec.at(10).find('|', base + 1);
						DP = stoi(str_vec.at(10).substr(str_vec.at(10).find('|', support_num + 1) + 1));
						support_num = std::stoi(str_vec.at(10).substr(support_num + 1));
						if((int)support_num>=paras->minReadsNumSupportSV){
							AF = (double)support_num / DP;
							if(AF>=paras->gt_homo_ratio) gt_str = GT_HOMOZYGOUS_STR;
							else if(AF>=paras->gt_hete_ratio) gt_str = GT_HETEROZYGOUS_STR;
							else gt_str = GT_NOZYGOUS_STR;
							stringstream ss;
							ss << fixed << setprecision(3) << AF;

							reg_str = chr + ":" + start_pos + "-" + start_pos;
							pthread_mutex_lock(&mutex_fai);
							seq = fai_fetch(fai, reg_str.c_str(), &seq_len);
							pthread_mutex_unlock(&mutex_fai);
							seq_str = seq;
							free(seq);

							ref = seq_str;
							alt = str_vec.at(7);
							qual = ".";	
							filter = "PASS";

							SVLEN = split(str_vec.at(15), ";");
							sv_len = (SVLEN.at(2)).substr(7);
							format = str_vec.at(16);
							format_val = gt_str + ":" + to_string(DP-support_num) + "," + to_string(support_num) + ":" + to_string(DP);

							line_vcf = chr + "\t" + start_pos + "\t" + id + "\t" + ref + "\t" + alt + "\t" + qual + "\t" + filter + "\t" + "SVTYPE=" + sv_type + ";" + "SVLEN=" + sv_len + ";" + "END=" + end_pos + ";" + "DP=" + to_string(DP) + ";" +"AF=" + ss.str() + ";" + "ULDist=" + sv_len + ";" + "LDISCOV=READS" + "\t" + format + "\t" + format_val;

							outfile << line_vcf << endl;
						}
					}
				}
			}
			//free(seq);
		}
	infile.close();
	outfile.close();

	// release memory
	for(i=0; i<(int32_t)bnd_vec.size(); i++) delete bnd_vec.at(i);
	vector<BND_t*>().swap(bnd_vec);
}

// save VCF header
void Genome::saveVCFHeader(ofstream &fp, string &sample_str){
	time_t rawtime;
	struct tm *timeinfo;
	char buffer [128];

	// get the local time
	time(&rawtime);
	timeinfo = localtime (&rawtime);
	strftime(buffer, sizeof(buffer), "%Y%m%d", timeinfo);

	fp << "##fileformat=VCFv" << VCF_VERSION << endl;
	fp << "##fileDate=" << buffer << endl;
	fp << "##source=" << PROG_NAME << " " << PROG_VERSION << endl;
	fp << "##reference=" << paras->refFile << endl;
	fp << "##PG=\"" << paras->pg_cmd_str << "\"" << endl;

	saveVCFContigHeader(fp);	// contigs information

	//fp << "##phasing=None" << endl;
	if(paras->phasing_flag) fp << "##phasing=Local phased" << endl;
	else fp << "##phasing=None" << endl;
//	fp << "##ALT=<ID=INS,Description=\"Insertion\">" << endl;
//	fp << "##ALT=<ID=DEL,Description=\"Deletion\">" << endl;
//	fp << "##ALT=<ID=DUP,Description=\"Duplication\">" << endl;
//	fp << "##ALT=<ID=INV,Description=\"Inversion\">" << endl;
//	fp << "##ALT=<ID=INVDUP,Description=\"Inverted DUP with unknown boundaries\">" << endl;
//	fp << "##ALT=<ID=TRA,Description=\"Translocation\">" << endl;
//	fp << "##ALT=<ID=BND,Description=\"Breakend\">" << endl;
//	fp << "##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Chromosome for END coordinate in case of a translocation\">" << endl;
	fp << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant described in this record\">" << endl;
	fp << "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">" << endl;
	fp << "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">" << endl;
	fp << "##INFO=<ID=DUPNUM,Number=1,Type=Integer,Description=\"Copy number of DUP\">" << endl;
	fp << "##INFO=<ID=MATEID,Number=.,Type=String,Description=\"ID of mate breakends\">" << endl;
	fp << "##INFO=<ID=MATEDIST,Number=1,Type=Integer,Description=\"Distance to the mate breakend for mates on the same contig\">" << endl;
	fp << "##INFO=<ID=LDISCOV,Number=1,Type=String,Description=\"Variant discover level: ALN for variants are discovered from consensus alignments, RES-ALN for variants are discovered from rescued consensus alignments around variant regions, READS for variants are discovered from reads alignments directly\">" << endl;
	fp << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">" << endl;
	fp << "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">" << endl;
	fp << "##INFO=<ID=ULDist,Number=1,Type=Integer,Description=\"UltraLarge Distance\">" << endl;
//	fp << "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">" << endl;
//	fp << "##FILTER=<ID=q10,Description=\"Quality below 10\">" << endl;
	if(paras->phasing_flag) fp << "##FORMAT=<ID=PS,Number=.,Type=String,Description=\"Phase set\">" << endl;
	fp << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl;
	fp << "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Read depth per allele\">" << endl;
	fp << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth at this position for this sample\">" << endl;
	fp << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << sample_str << endl;
}

// save VCF contig header
void Genome::saveVCFContigHeader(ofstream &fp){
	Chrome *chr;
	for(size_t i=0; i<chromeVector.size(); i++){
		chr = chromeVector.at(i);
		fp << "##contig=<ID=" << chr->chrname << ",length=" << chr->chrlen << ">" << endl;
	}
}

// compute the statistics for 'detect' command
void Genome::computeVarNumStatDetect(){
	size_t total, indel_num, clipReg_num;

	indel_num = getLineCount(out_filename_detect_indel);
	clipReg_num = getLineCount(out_filename_detect_clipReg);

	total = indel_num + clipReg_num;

	cout << "############## Brief statistics for variants detection ##############" << endl;
	cout << "There are " << total << " variant candidates detected in total:" << endl;
	cout << "\t" << indel_num << " candidate indels" << endl;
	cout << "\t" << clipReg_num << " candidate clipping regions" << endl;
}

// compute statistics for 'cns' command
void Genome::computeVarNumStatCons(){
	int32_t total, total_indel, total_clipReg, total_succ_indel, total_fail_indel, total_succ_clipReg, total_fail_clipReg, succ_tmp, fail_tmp;
	vector<int32_t> num_vec;
	Chrome *chr;
	string filename;

	total_indel = total_clipReg = total_succ_indel = total_fail_indel = total_succ_clipReg = total_fail_clipReg = 0;
	for(size_t i=0; i<chromeVector.size(); i++){
		chr = chromeVector.at(i);
		if(chr->valid_flag){
			filename = chr->getVarcandIndelFilename();
			num_vec = getSuccFailNumCns(filename);
			succ_tmp = num_vec.at(0);
			fail_tmp = num_vec.at(1);
			total_succ_indel += succ_tmp;
			total_fail_indel += fail_tmp;
			total_indel += succ_tmp + fail_tmp;

			filename = chr->getVarcandClipregFilename();
			num_vec = getSuccFailNumCns(filename);
			succ_tmp = num_vec.at(0);
			fail_tmp = num_vec.at(1);
			total_succ_clipReg += succ_tmp;
			total_fail_clipReg += fail_tmp;
			total_clipReg += succ_tmp + fail_tmp;
		}
	}
	total = total_indel + total_clipReg;

	cout << "############## Brief statistics for local consensus ##############" << endl;
	cout << "There are " << total << " regions in total:" << endl;
	cout << "\t" << total_indel << " indel regions: " << total_succ_indel << " successful, " << total_fail_indel << " failed" << endl;
	cout << "\t" << total_clipReg << " clipping regions: " << total_succ_clipReg << " successful, " << total_fail_clipReg << " failed" << endl;
}

vector<int32_t> Genome::getSuccFailNumCns(string &filename){
	ifstream infile;
	string line;
	vector<string> str_vec;
	vector<int32_t> succ_fail_vec;
	int32_t num_succ, num_fail;

	if(isFileExist(filename)){
		infile.open(filename);
		if(!infile.is_open()){
			cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << filename << endl;
			exit(1);
		}

		// read each line and check the success and failed flag
		num_succ = num_fail = 0;
		while(getline(infile, line))
			if(line.size()>0 and line.at(0)!='#'){
				str_vec = split(line, "\t");
				if(str_vec.at(6).compare(CNS_SUCCESS)==0) num_succ ++;
				else num_fail ++;
			}

		succ_fail_vec.push_back(num_succ);
		succ_fail_vec.push_back(num_fail);

		infile.close();
	}else{
		succ_fail_vec.push_back(0);
		succ_fail_vec.push_back(0);
	}

	return succ_fail_vec;
}

// compute statistics for 'call' command
void Genome::computeVarNumStatCall(){
	ifstream infile;
	string line;
	vector<string> str_vec;
	vector<int32_t> succ_fail_vec;
	int32_t total, num_ins, num_del, num_inv, num_dup, num_tra, num_unc;
	size_t sv_type;

	infile.open(out_filename_result_vars);
	if(!infile.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << out_filename_result_vars << endl;
		exit(1);
	}

	num_ins = num_del = num_inv = num_dup = num_tra = num_unc = 0;
	while(getline(infile, line)){
		if(line.size()>0 and line.at(0)!='#'){
			sv_type = getSVTypeSingleLine(line);
			switch(sv_type){
				case VAR_UNC: num_unc ++; break;
				case VAR_INS: num_ins ++; break;
				case VAR_DEL: num_del ++; break;
				case VAR_DUP: num_dup ++; break;
				case VAR_INV: num_inv ++; break;
				case VAR_TRA:
				case VAR_BND: num_tra ++; break;
				//case VAR_INV_TRA: num_ins ++; break;
				//case VAR_MIX: num_ins ++; break;
			}
		}
	}
	infile.close();

	total = num_unc + num_ins + num_del + num_inv + num_dup + num_tra;

	cout << "############## Brief statistics for variants call ##############" << endl;
	cout << "There are " << total << " variants in total:" << endl;
	cout << "\t" << "insertions: " << num_ins << endl;
	cout << "\t" << "deletions: " << num_del << endl;
	cout << "\t" << "tandem duplications: " << num_dup << endl;
	cout << "\t" << "inversions: " << num_inv << endl;
	cout << "\t" << "translocations: " << num_tra << endl;
	cout << "\t" << "uncertain: " << num_unc << endl;

	if(paras->phasing_flag){
		cout << "Number of phased sets : " << ps_num << endl;
		cout << "number of phased items: " << phased_var_num << endl;
	}
}

size_t Genome::getSVTypeSingleLine(string &line){
	size_t sv_type = VAR_UNC;
	vector<string> str_vec;
	string str_tmp, str_tmp2;

	str_vec = split(line, "\t");
	str_tmp2 = str_vec.at(4);
	if(str_tmp2.at(0)=='<') str_tmp = str_tmp2.substr(1, str_tmp2.size()-2);
	else str_tmp = str_tmp2;
	if(str_tmp.compare(VAR_UNC_STR)==0 or str_tmp.compare(VAR_UNC_STR1)==0 or str_tmp.compare(VAR_UNC_STR2)==0){
		sv_type = VAR_UNC;
	}else if(str_tmp.compare(VAR_INS_STR)==0 or str_tmp.compare(VAR_INS_STR1)==0 or str_tmp.compare(VAR_INS_STR2)==0){
		sv_type = VAR_INS;
	}else if(str_tmp.compare(VAR_DEL_STR)==0 or str_tmp.compare(VAR_DEL_STR1)==0 or str_tmp.compare(VAR_DEL_STR2)==0){
		sv_type = VAR_DEL;
	}else if(str_tmp.compare(VAR_DUP_STR)==0 or str_tmp.compare(VAR_DUP_STR1)==0 or str_tmp.compare(VAR_DUP_STR2)==0){
		sv_type = VAR_DUP;
	}else if(str_tmp.compare(VAR_INV_STR)==0 or str_tmp.compare(VAR_INV_STR1)==0 or str_tmp.compare(VAR_INV_STR2)==0){
		sv_type = VAR_INV;
	}else{
		str_tmp2 = str_vec.at(7);
		if(str_tmp2.at(0)=='<') str_tmp = str_tmp2.substr(1, str_tmp2.size()-2);
		else str_tmp = str_tmp2;
		if(str_vec.size()>MAX_BED_COLS_NUM){
			if(str_tmp.compare(VAR_INV_STR)==0 or str_tmp.compare(VAR_INV_STR1)==0 or str_tmp.compare(VAR_INV_STR2)==0) sv_type = VAR_INV;
			else if(str_tmp.compare(VAR_TRA_STR)==0 or str_tmp.compare(VAR_TRA_STR1)==0 or str_tmp.compare(VAR_TRA_STR2)==0) sv_type = VAR_TRA;
			else if(str_tmp.compare(VAR_BND_STR)==0) sv_type = VAR_BND;
			else sv_type = VAR_MIX;
		}else{
			sv_type = VAR_MIX;
		}
	}

	return sv_type;
}

void Genome::sortVarResults(string &infilename, int32_t filetype){
	string outfilename;
	vector<SV_item*> sv_vec;
	vector<vector<SV_item*>> subsets;

	outfilename = infilename + ".sotred";

	if(filetype==0) sv_vec = loadDataBED(infilename); // BED format
	else if(filetype==1) sv_vec = loadDataVcf(infilename); // VCF format
	else{
		cerr << __func__ << ", line=" << __LINE__ << ": invalid file type: " << filetype << endl;
		exit(1);
	}

	subsets = constructSubsetByChr(sv_vec, fai);
	sortSVitem(subsets);
	rmRedundantSVitem(subsets, MAX_DIST_MATCH_CLIP_POS, MIN_SIZE_RATIO_REDUNDANT_FILTER, MIN_SEQSIM_REDUNDANT_FILTER, fai); // remove redundant items
	if(filetype==1 and paras->phasing_flag) computePhaseSetNum(subsets);
	outputResult(outfilename, subsets, filetype);

	remove(infilename.c_str());
	rename(outfilename.c_str(), infilename.c_str());

	destroyData(sv_vec);
}

void Genome::computePhaseSetNum(vector<vector<SV_item*>> &subsets){
	ps_num = phased_var_num = 0;
	vector<int32_t> num_vec;

	for(size_t i=0; i<subsets.size(); i++){
		num_vec = computePhaseSetNumSubset(subsets.at(i));
		ps_num += num_vec.at(0);
		phased_var_num += num_vec.at(1);
		//cout << "subset [" << i << "]: ps_num_tmp=" << num_vec.at(0) << ", phased_item_num_tmp=" << num_vec.at(1) << endl;
	}
	//cout << "ps_num=" << ps_num << ", phased_var_num=" << phased_var_num << endl;
}

vector<int32_t> Genome::computePhaseSetNumSubset(vector<SV_item*> &subset){
	SV_item *item;
	vector<int32_t> ret_num_vec;
	int32_t ps_num, item_num;
	map<string, int32_t> ps_id_map;
	string ps_str;
	vector<string> ps_vec;

	item_num = 0;
	for(size_t i=0; i<subset.size(); i++){
		item = subset.at(i);
		if(item->valid_flag){
			ps_str = extractPSFromVCFLine(item->info);
			if(ps_str.size()>0){
				item_num ++;
				ps_vec = split(ps_str, ",");
				for(size_t j=0; j<ps_vec.size(); j++){
					ps_id_map[ps_vec.at(j)] ++;
				}
			}
		}
	}

	ps_num = 0;
	for(const auto &ps : ps_id_map) if(ps.second>=2) ps_num ++;

	ret_num_vec.push_back(ps_num);
	ret_num_vec.push_back(item_num);

	return ret_num_vec;
}

string Genome::extractPSFromVCFLine(string &line){
	string ps_str, str;
	vector<string> str_vec, info_vec;
	uint64_t idx;

	str_vec = split(line, "\t");
	info_vec = split(str_vec.at(7), ";");
	for(size_t i=0; i<info_vec.size(); i++){
		str = info_vec.at(i);
		idx = str.find("PS=");
		if(idx!=string::npos){ // found
			ps_str = str.substr(idx+3);
			break;
		}
	}

	return ps_str;
}

void Genome::outputResult(string &outfilename, vector<vector<SV_item*>> &subsets, int32_t filetype){
	ofstream outfile;
	vector<SV_item*> sv_vec;
	string header_line_bed, header_line_bedpe;

	outfile.open(outfilename);
	if(!outfile.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << outfilename << endl;
		exit(1);
	}

	// header line
	if(filetype==0){ // BED format
		header_line_bed = getCallFileHeaderBed(paras->sample);
		header_line_bedpe = getCallFileHeaderBedpe(paras->sample);
		outfile << header_line_bed << endl;
		outfile << header_line_bedpe << endl;
	}else if(filetype==1){ // VCF format
		// save VCF header
		saveVCFHeader(outfile, paras->sample);
	}else{
		cerr << __func__ << ", line=" << __LINE__ << ": invalid file type: " << filetype << endl;
		exit(1);
	}

	for(size_t i=0; i<subsets.size(); i++){
		sv_vec = subsets.at(i);
		for(size_t j=0; j<sv_vec.size(); j++) if(sv_vec.at(j)->valid_flag) outfile << sv_vec.at(j)->info << endl;
	}
	outfile.close();
}

// remove temporary monitor files
void Genome::removeTempMonitorFiles(){
	if(isFileExist(work_finish_filename_cns)) remove(work_finish_filename_cns.c_str());
	if(isFileExist(work_finish_filename_call)) remove(work_finish_filename_call.c_str());
}

// remove temporary results
void Genome::removeTempResults(){
	string cmd = "rm -rf " + out_dir_detect + " " + out_dir_cns + " " + out_dir_call;
	cout << "Cleaning temporary results:" << endl;
	cout << "\t" << out_dir_detect << endl << "\t" << out_dir_cns << endl << "\t" << out_dir_call << endl;
	cout << "If you want to keep the temporary results, please specify the `--no-clean` option." << endl;
	system(cmd.c_str());
}
