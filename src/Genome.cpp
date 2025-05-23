#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <string>
#include <pthread.h>
#include <htslib/thread_pool.h>

#include "Genome.h"
#include "Thread.h"

using namespace std;

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
	blat_aln_info_filename_tra  = out_dir_tra + "/" + "blat_aln_info_tra";

	work_finish_filename = out_dir + "/" + "work_finished";

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

// get genome size
int Genome::getGenomeSize(){
    genomeSize = 0;
    for(int i=0; i<header->n_targets; i++) genomeSize += header->target_len[i];
	return genomeSize;
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

// save the genome blocks to file
void Genome::saveGenomeBlocksToFile(){
	vector<Chrome*>::iterator chr;
	for(chr=chromeVector.begin(); chr!=chromeVector.end(); chr++)
		(*chr)->saveChrBlocksToFile();
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
		chr->chrFillDataEst(SIZE_EST_OP);
		if(paras->reg_sum_size_est>=paras->max_reg_sum_size_est) break;

		//paras->total_depth += paras->chr_mean_depth;
	}
	//paras->chrome_num = chromeVector.size() - 1;
	// size estimate
	paras->estimate(SIZE_EST_OP);

	// fill the num data using each chromosome
	paras->reg_sum_size_est = 0;
	for(i=0; i<chromeVector.size(); i++){
		chr = chromeVector.at(i);
		chr->chrFillDataEst(NUM_EST_OP);
		if(paras->reg_sum_size_est>=paras->max_reg_sum_size_est) break;
	}
	// num estimate
	paras->estimate(NUM_EST_OP);
}

// detect variants for genome
int Genome::genomeDetect(){
	Chrome *chr;
	for(size_t i=0; i<chromeVector.size(); i++){
		chr = chromeVector.at(i);
		if((chr->decoy_flag==false or paras->include_decoy) and (chr->alt_flag==false or paras->include_alt))
			chr->chrDetect();
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

// remove redundant mate clipping regions
void Genome::removeRedundantMateClipReg(){
	for(size_t i=0; i<chromeVector.size(); i++)
		genomeRemoveRedundantClipReg(chromeVector.at(i), chromeVector);

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
	for(size_t i=0; i<chromeVector.size(); i++)
		genomeRemoveFPIndelSnvInClipReg(chromeVector.at(i), chromeVector);
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
						if(chr!=chr_tmp and chrname_arr[j].size()>0 and chr_tmp->chrname.compare(chrname_arr[j])==0){
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
		//if(chr->process_block_num)
		{
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
	//startWorkProcessMonitor(work_finish_filename, paras->monitoring_proc_names_cns, paras->max_proc_running_minutes_cns);

	// begin consensus
	if(!paras->cns_work_vec.empty()) cout << "[" << time.getTime() << "]: start local consensus ..." << endl;
	processConsWork();

//	for(size_t i=0; i<1001; i++){
//		cout << i << "\t" << sig_num_arr[i] << endl;
//	}

	cout << "[" << time.getTime() << "]: Finalizing consensus ..." << endl;

	computeVarNumStatCons(); // compute statistics for cns command

	if(!paras->cns_work_vec.empty()) destroyConsWorkOptVec(paras->cns_work_vec);

	// reset consensus data
	for(size_t i=0; i<chromeVector.size(); i++){
		chr = chromeVector.at(i);
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
		if((chr->decoy_flag==false or paras->include_decoy==false) and (chr->alt_flag==false or paras->include_alt==false)){ // non-decoy, non-alt
			if(chr->chrname.compare(CHR_X_STR1)==0 or chr->chrname.compare(CHR_X_STR2)==0 or chr->chrname.compare(CHR_Y_STR1)==0 or chr->chrname.compare(CHR_Y_STR2)==0){
				chr->chrLoadDataCons();  // load the variant data
				chr->chrGenerateLocalConsWorkOpt();     // generate local consensus work
			}
		}else if((chr->decoy_flag and paras->include_decoy) or (chr->alt_flag and paras->include_alt)){ // decoy and alt
			chr->chrLoadDataCons();  // load the variant data
			chr->chrGenerateLocalConsWorkOpt();     // generate local consensus work
		}
	}

	// load other chrs
	for(size_t i=0; i<chromeVector.size(); i++){
		chr = chromeVector.at(i);
		if(chr->decoy_flag==false and chr->alt_flag==false)
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
//		if(cns_work_opt->contigfilename.compare("debug_tmp_chr1_1-50m_20250413/2_cns/1/cns_1_9683995-2256585.fa")==0){
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
		cns_work->min_identity_match = paras->min_identity_match;

		sv_len_sum = 0;
		min_pos = INT_MAX, max_pos = 0;
		for(j=0; j<cns_work_opt->arr_size; j++) {
			sv_len_sum += cns_work_opt->var_array[j]->sv_len;
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
		cns_work->inBamFile = paras->inBamFile;
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
	
	// initialize the variables for process monitor killed blat works
	//initMonitorKilledBlatWorkMem();
	//initMonitorKilledMinimap2WorkMem();

	// load clipping region information
	genomeLoadMateClipRegData();

	// collect regions
	cout << "Loading regions ..." << endl;
	genomeCollectCallWork();

	cout << "Number of regions to be processed: " << paras->call_work_num << endl;

	// invoke the monitor of consensus work process
	//startWorkProcessMonitor(work_finish_filename, paras->monitoring_proc_names_call, paras->max_proc_running_minutes_call);

	// blat alignment work
	time.setStartTime();
	processAlnWork();
	time.printElapsedTime();

	// process call work
	time.setStartTime();
	processCallWork();
	time.printElapsedTime();

	// finish call work
	//genomeFinishCallWork();

	// call TRA according to mate clip regions
	//genomeCallTra(); // 2023-12-25

	// fill variant sequence
	//cout << "4444444444444" << endl;
	cout << "[" << time.getTime() << "]: fill variant sequences... " << endl;
	genomeFillVarseq();

	// generate work process finish file
	//generateFile(work_finish_filename);

	// save SV to file
	//cout << "5555555555555" << endl;
	cout << "[" << time.getTime() << "]: save call result to file... " << endl;
	genomeSaveCallSV2File();

	// merge call results into single file
	//cout << "6666666666666" << endl;
	cout << "[" << time.getTime() << "]: merge call result... " << endl;
	mergeCallResult();

	// sort variant items in BED format
	sortVarResults(out_filename_result_vars, 0);

	// save results in VCF file format
	saveResultVCF();

	// sort variant items in VCF format
	sortVarResults(out_filename_result_vars_vcf, 1);

	// compute statistics for call command
	cout << "[" << time.getTime() << "]: compute variant NUMBER statistics... " << endl;
	computeVarNumStatCall();

	//releaseMonitorKilledBlatWorkMem();
	//releaseMonitorKilledMinimap2WorkMem();

	return 0;
}

// initialize variables for monitor killed minimap2 work
void Genome::initMonitorKilledMinimap2WorkMem(){
	string line, alnfilename, contigfilename, refseqfilename, old_out_dir, killed_minimap2_work_filename, tmp_filename;
	ifstream infile;
	killedMinimap2Work_t *killed_minimap2_work;
	vector<string> str_vec;
	bool file_exist_flag;

	pthread_mutex_init(&paras->mtx_killed_minimap2_work, NULL);

	killed_minimap2_work_filename = out_dir + "/" + paras->killed_minimap2_work_filename;
	file_exist_flag = isFileExist(killed_minimap2_work_filename);
	if(file_exist_flag){
		tmp_filename = killed_minimap2_work_filename + "_tmp";
		rename(killed_minimap2_work_filename.c_str(), tmp_filename.c_str());
	}

	paras->killed_minimap2_work_file.open(killed_minimap2_work_filename);
	if(!paras->killed_minimap2_work_file.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << paras->killed_minimap2_work_filename << endl;
		exit(1);
	}

	if(file_exist_flag){
		infile.open(tmp_filename);
		if(!infile.is_open()){
			cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << tmp_filename << endl;
			exit(1);
		}
		// read each line and save to the output file
		while(getline(infile, line)){
			if(line.size()>0 and line.at(0)!='#'){
				str_vec = split(line, "\t");
				if(str_vec.size()==3){
					if(old_out_dir.size()==0) old_out_dir = getOldOutDirname(str_vec.at(0), paras->out_dir_call);
					alnfilename = getUpdatedItemFilename(str_vec.at(0), paras->outDir, old_out_dir);
					contigfilename = getUpdatedItemFilename(str_vec.at(1), paras->outDir, old_out_dir);
					refseqfilename = getUpdatedItemFilename(str_vec.at(2), paras->outDir, old_out_dir);

					killed_minimap2_work = new killedMinimap2Work_t();
					killed_minimap2_work->alnfilename = alnfilename;		// alnfilename;
					killed_minimap2_work->ctgfilename = contigfilename;		// ctgfilename;
					killed_minimap2_work->refseqfilename = refseqfilename;	// refseqfilename;

					paras->killed_minimap2_work_vec.push_back(killed_minimap2_work);
					paras->killed_minimap2_work_file << line << endl;
				}
			}
		}
		infile.close();

		remove(tmp_filename.c_str());
	}
}

// initialize variables for monitor killed blat work
void Genome::initMonitorKilledBlatWorkMem(){
	string line, alnfilename, contigfilename, refseqfilename, old_out_dir, killed_blat_work_filename, tmp_filename;
	ifstream infile;
	killedBlatWork_t *killed_blat_work;
	vector<string> str_vec;
	bool file_exist_flag;

	pthread_mutex_init(&paras->mtx_killed_blat_work, NULL);

	killed_blat_work_filename = out_dir + "/" + paras->killed_blat_work_filename;
	file_exist_flag = isFileExist(killed_blat_work_filename);
	if(file_exist_flag){
		tmp_filename = killed_blat_work_filename + "_tmp";
		rename(killed_blat_work_filename.c_str(), tmp_filename.c_str());
	}

	paras->killed_blat_work_file.open(killed_blat_work_filename);
	if(!paras->killed_blat_work_file.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << paras->killed_blat_work_filename << endl;
		exit(1);
	}

	if(file_exist_flag){
		infile.open(tmp_filename);
		if(!infile.is_open()){
			cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << tmp_filename << endl;
			exit(1);
		}
		// read each line and save to the output file
		while(getline(infile, line)){
			if(line.size()>0 and line.at(0)!='#'){
				str_vec = split(line, "\t");
				if(str_vec.size()==3){
					if(old_out_dir.size()==0) old_out_dir = getOldOutDirname(str_vec.at(0), paras->out_dir_call);
					alnfilename = getUpdatedItemFilename(str_vec.at(0), paras->outDir, old_out_dir);
					contigfilename = getUpdatedItemFilename(str_vec.at(1), paras->outDir, old_out_dir);
					refseqfilename = getUpdatedItemFilename(str_vec.at(2), paras->outDir, old_out_dir);

					killed_blat_work = new killedBlatWork_t();
					killed_blat_work->alnfilename = alnfilename;		// alnfilename;
					killed_blat_work->ctgfilename = contigfilename;		// ctgfilename;
					killed_blat_work->refseqfilename = refseqfilename;	// refseqfilename;

					paras->killed_blat_work_vec.push_back(killed_blat_work);
					paras->killed_blat_work_file << line << endl;
				}
			}
		}
		infile.close();

		remove(tmp_filename.c_str());
	}
}

// initialize variables for monitor killed blat work
void Genome::releaseMonitorKilledMinimap2WorkMem(){
	paras->killed_minimap2_work_file.close();

	for(size_t i=0; i<paras->killed_minimap2_work_vec.size(); i++) delete paras->killed_minimap2_work_vec.at(i);
	vector<killedMinimap2Work_t*>().swap(paras->killed_minimap2_work_vec);
}

// initialize variables for monitor killed blat work
void Genome::releaseMonitorKilledBlatWorkMem(){
	paras->killed_blat_work_file.close();

	for(size_t i=0; i<paras->killed_blat_work_vec.size(); i++) delete paras->killed_blat_work_vec.at(i);
	vector<killedBlatWork_t*>().swap(paras->killed_blat_work_vec);
}

// collect call work
void Genome::genomeCollectCallWork(){
	Chrome *chr;

	if(chromeVector.size()==0) return;

	// load X, Y, hs37d5, alt
	for(size_t i=0; i<chromeVector.size(); i++){
		chr = chromeVector.at(i);
		if((chr->decoy_flag==false or paras->include_decoy==false) and (chr->alt_flag==false or paras->include_alt==false)){ // non-decoy, non-alt
			if(chr->chrname.compare(CHR_X_STR1)==0 or chr->chrname.compare(CHR_X_STR2)==0 or chr->chrname.compare(CHR_Y_STR1)==0 or chr->chrname.compare(CHR_Y_STR2)==0){
				chr->chrLoadDataCall();
				chr->chrCollectCallWork();
			}
		}else if((chr->decoy_flag and paras->include_decoy) or (chr->alt_flag and paras->include_alt)){ // decoy or alt
			chr->chrLoadDataCall();
			chr->chrCollectCallWork();
		}
	}

	// load other chrs
	for(size_t i=0; i<chromeVector.size(); i++){
		chr = chromeVector.at(i);
		if(chr->decoy_flag==false and chr->alt_flag==false){ // non-decoy and no alt
			if(chr->chrname.compare(CHR_X_STR1)!=0 and chr->chrname.compare(CHR_X_STR2)!=0 and chr->chrname.compare(CHR_Y_STR1)!=0 and chr->chrname.compare(CHR_Y_STR2)!=0){
				chr->chrLoadDataCall();
				chr->chrCollectCallWork();
			}
		}
	}
}

void Genome::genomeFinishCallWork(){
	Chrome *chr;

	for(size_t i=0; i<chromeVector.size(); i++){
		chr = chromeVector.at(i);
		chr->resetBlatVarcandFiles();
	}
	vector<varCand*>().swap(paras->call_work_vec);
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

		// DUP not precise (CCS30x): blat_contig_1_1180102-1180675.sim4, blat_contig_1_1183812-1185067.sim4, blat_contig_1_1317611-1318285.sim4
//		if(var_cand->alnfilename.compare("output_ccs_v1.0.1_20210528/3_call/1/blat_1_1860801-1869285.sim4")!=0){
//			continue;
//		}

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
int Genome::processBlatAlnWork(){
	callWork_opt *call_work_opt;
	varCand *var_cand;
	size_t num_work, num_work_percent;

	if(paras->call_work_vec.empty()) return 0;  // no align work, then return directly

	cout << "Begin blat alignment ..." << endl;

	hts_tpool *p = hts_tpool_init(paras->num_threads);
	hts_tpool_process *q = hts_tpool_process_init(p, paras->num_threads*2, 1);

	pthread_mutex_init(&paras->mtx_call_workDone_num, NULL);

	paras->call_workDone_num = 0;
	num_work = paras->call_work_vec.size();
	num_work_percent = num_work / (paras->num_parts_progress >> 1);
	if(num_work_percent==0) num_work_percent = 1;
	for(size_t i=0; i<num_work; i++){
		var_cand = paras->call_work_vec.at(i);

		// DUP not precise (CCS30x): blat_contig_1_1180102-1180675.sim4, blat_contig_1_1183812-1185067.sim4, blat_contig_1_1317611-1318285.sim4
//		if(var_cand->alnfilename.compare("output_ccs_v1.0.1_20210528/3_call/1/blat_1_1860801-1869285.sim4")!=0){
//			continue;
//		}

		call_work_opt = new callWork_opt();
		call_work_opt->var_cand = var_cand;
		call_work_opt->work_id = i;
		call_work_opt->num_work = num_work;
		call_work_opt->num_work_percent = num_work_percent;
		call_work_opt->p_call_workDone_num = &(paras->call_workDone_num);
		call_work_opt->p_mtx_call_workDone_num = &(paras->mtx_call_workDone_num);

		hts_tpool_dispatch(p, q, processSingleBlatAlnWork, call_work_opt);
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
//		 if(var_cand->alnfilename.compare("debug_tmp_chr1_1-50m_20250413/3_call/1/minimap2_1_2256290-2256585.paf")==0){
//		 	cout << var_cand->alnfilename << endl;
//		 }else continue;

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
		chr->chrLoadMateClipRegData();
	}
}

// find indels back from mate clip regions
void Genome::recallIndelsFromTRA(){
	size_t i, j, k, m, n, round_num, success_num, total, total_success;
	Chrome *chr;
	vector<mateClipReg_t*> mate_clipReg_vec;
	mateClipReg_t *clip_reg;
	varCand *var_cand;
	int32_t segIdx_start, segIdx_end, ref_dist, query_dist;
	blat_aln_t *blat_aln;
	aln_seg_t *seg1, *seg2;
	reg_t *reg, *reg1, *reg2, *reg3;
	string chrname_tmp;
	bool overlap_flag, exist_flag;

	total = total_success = 0;
	for(i=0; i<chromeVector.size(); i++){
		success_num = 0;
		chr = chromeVector.at(i);
		mate_clipReg_vec = chr->mateClipRegVector;
		for(j=0; j<mate_clipReg_vec.size(); j++){
			//if(i>=0 and j>=247)
			//	cout << "line=" << __LINE__ << ", i=" << i << ", j=" << j << ", " << chr->chrname << endl;

			total ++;
			clip_reg = chr->mateClipRegVector.at(j);
			reg1 = reg2 = reg3 = NULL;
			if(clip_reg->valid_flag and (clip_reg->call_success_flag==false and clip_reg->tra_rescue_success_flag==false) and clip_reg->sv_type==VAR_TRA){ // valid TRA
				for(round_num=0; round_num<3; round_num++){
					if(round_num==0){
						var_cand = clip_reg->left_var_cand_tra;
						chrname_tmp = clip_reg->leftClipReg ? clip_reg->leftClipReg->chrname : clip_reg->leftClipReg2->chrname;
					}else if(round_num==1){
						var_cand = clip_reg->right_var_cand_tra;
						chrname_tmp = clip_reg->rightClipReg ? clip_reg->rightClipReg->chrname : clip_reg->rightClipReg2->chrname;
					}else{
						var_cand = clip_reg->var_cand;
						if(var_cand) chrname_tmp = var_cand->chrname;
						else chrname_tmp = "";
					}

					if(var_cand){
						var_cand->loadBlatAlnData();

						segIdx_start = segIdx_end = -1;
						for(k=0; k<var_cand->blat_aln_vec.size(); k++){
							blat_aln = var_cand->blat_aln_vec.at(k);
							if(blat_aln->valid_aln){
								for(m=1; m<blat_aln->aln_segs.size(); m++){
									seg1 = blat_aln->aln_segs.at(m-1);
									seg2 = blat_aln->aln_segs.at(m);
									if((seg1->ref_end>=var_cand->leftClipRefPos-CLIP_END_EXTEND_SIZE and seg1->ref_end<=var_cand->rightClipRefPos+CLIP_END_EXTEND_SIZE) and (seg2->ref_start>=var_cand->leftClipRefPos-CLIP_END_EXTEND_SIZE and seg2->ref_start<=var_cand->rightClipRefPos+CLIP_END_EXTEND_SIZE)){
										segIdx_start = m - 1;

										// compute segIdx_end
										segIdx_end = -1;
										for(n=segIdx_start+1; n<blat_aln->aln_segs.size(); n++){
											seg2 = blat_aln->aln_segs.at(n);
											if(seg2->ref_end-seg2->ref_start+1>=MIN_VALID_BLAT_SEG_SIZE){ // ignore short align segments
												segIdx_end = n;
												break;
											}
										}
										if(segIdx_end!=-1) break;
									}
								}

								if(segIdx_start!=-1 and segIdx_end!=-1){
									seg1 = blat_aln->aln_segs.at(segIdx_start);
									seg2 = blat_aln->aln_segs.at(segIdx_end);
									//overlap_flag = isOverlappedPos(var_cand->leftClipRefPos, var_cand->rightClipRefPos, seg1->ref_end, seg2->ref_start);
									overlap_flag = isOverlappedPos(var_cand->leftClipRefPos-CLIP_END_EXTEND_SIZE, var_cand->rightClipRefPos+CLIP_END_EXTEND_SIZE, seg1->ref_end-CLIP_END_EXTEND_SIZE, seg2->ref_start+CLIP_END_EXTEND_SIZE);
									exist_flag = isExistSuccessClipReg(var_cand, clip_reg);
									if(overlap_flag and exist_flag==false){ // indel
										var_cand->call_success = true;
										var_cand->clip_reg_flag = false;

										reg = new reg_t();
										reg->chrname = chrname_tmp;
										reg->startRefPos =  seg1->ref_end;
										reg->endRefPos =  seg2->ref_start;
										reg->startLocalRefPos =  seg1->subject_end;
										reg->endLocalRefPos =  seg2->subject_start;
										reg->startQueryPos =  seg1->query_end;
										reg->endQueryPos =  seg2->query_start;
										reg->aln_orient = blat_aln->aln_orient;
										reg->query_id = blat_aln->query_id;
										//reg->blat_aln_id = i;
										reg->blat_aln_id = k;
										reg->minimap2_aln_id = -1;
										reg->call_success_status = true;
										reg->short_sv_flag = false;
										reg->zero_cov_flag = false;
										reg->aln_seg_end_flag = false;
										reg->query_pos_invalid_flag = false;
										reg->large_indel_flag = false;
										reg->gt_type = -1;
										reg->gt_seq = "";
										reg->AF = 0;
										reg->supp_num = reg->DP = 0;
										reg->discover_level = VAR_DISCOV_L_UNUSED;

										ref_dist = reg->endLocalRefPos - reg->startLocalRefPos + 1;
										query_dist = reg->endQueryPos - reg->startQueryPos + 1;
										reg->sv_len = query_dist - ref_dist + 1;
										if(reg->sv_len>0) reg->var_type = VAR_INS;
										else reg->var_type = VAR_DEL;

										if(round_num==0) reg1 = reg;
										else if(round_num==1) reg2 = reg;
										else reg3 = reg;

										break;
									}else
										segIdx_start = segIdx_end = -1;
								}
							}
						}
					}
				}

				if(reg1){ // left part
					clip_reg->left_var_cand_tra->varVec.push_back(reg1);

					if(clip_reg->leftClipReg) { delete clip_reg->leftClipReg; clip_reg->leftClipReg = NULL; clip_reg->leftMeanClipPos = 0; clip_reg->leftClipPosNum = 0; }
					if(clip_reg->leftClipReg2) { delete clip_reg->leftClipReg2; clip_reg->leftClipReg2 = NULL; clip_reg->leftMeanClipPos2 = 0; clip_reg->leftClipPosNum2 = 0; }
					clip_reg->leftClipRegNum = 0;
					clip_reg->left_var_cand_tra = NULL;
				}
				if(reg2){ // right part
					clip_reg->right_var_cand_tra->varVec.push_back(reg2);

					if(clip_reg->rightClipReg) { delete clip_reg->rightClipReg; clip_reg->rightClipReg = NULL; clip_reg->rightMeanClipPos = 0; clip_reg->rightClipPosNum = 0; }
					if(clip_reg->rightClipReg2) { delete clip_reg->rightClipReg2; clip_reg->rightClipReg2 = NULL; clip_reg->rightMeanClipPos2 = 0; clip_reg->rightClipPosNum2 = 0; }
					clip_reg->rightClipRegNum = 0;
					clip_reg->right_var_cand_tra = NULL;
				}
				if(reg3){ // one region
					clip_reg->var_cand->varVec.push_back(reg3);

					if(clip_reg->leftClipReg) { delete clip_reg->leftClipReg; clip_reg->leftClipReg = NULL; clip_reg->leftMeanClipPos = 0; clip_reg->leftClipPosNum = 0; }
					if(clip_reg->leftClipReg2) { delete clip_reg->leftClipReg2; clip_reg->leftClipReg2 = NULL; clip_reg->leftMeanClipPos2 = 0; clip_reg->leftClipPosNum2 = 0; }
					if(clip_reg->rightClipReg) { delete clip_reg->rightClipReg; clip_reg->rightClipReg = NULL; clip_reg->rightMeanClipPos = 0; clip_reg->rightClipPosNum = 0; }
					if(clip_reg->rightClipReg2) { delete clip_reg->rightClipReg2; clip_reg->rightClipReg2 = NULL; clip_reg->rightMeanClipPos2 = 0; clip_reg->rightClipPosNum2 = 0; }
					clip_reg->leftClipRegNum = 0;
					clip_reg->rightClipRegNum = 0;
					clip_reg->left_var_cand_tra = NULL;
					clip_reg->right_var_cand_tra = NULL;
					clip_reg->var_cand = NULL;
				}

				if(reg1 or reg2 or reg3){
					success_num ++;
					total_success ++;
					clip_reg->valid_flag = false;
				}else{
					//cout << "chr " << i << ", j=" << j << ", recall Indel failed" << endl;
					//printMateClipReg(clip_reg);
				}
			}
		}
		//cout << "chr " << i << ": success_num=" << success_num << endl;
	}
	//cout << "Genome: total=" << total << ", total_success=" << total_success << endl;
}

// determine whether the region is successful clipping regions
bool Genome::isExistSuccessClipReg(varCand *var_cand, mateClipReg_t *clip_reg_exclude){
	size_t i, j;
	Chrome *chr;
	vector<mateClipReg_t*> mate_clipReg_vec;
	mateClipReg_t *clip_reg;
	string chrname_tmp;
	int64_t start_pos, end_pos;
	bool success_flag, overlap_flag;

	success_flag = false;
	for(i=0; i<chromeVector.size(); i++){
		chr = chromeVector.at(i);
		mate_clipReg_vec = chr->mateClipRegVector;
		for(j=0; j<mate_clipReg_vec.size(); j++){
			clip_reg = chr->mateClipRegVector.at(j);
			if(clip_reg!=clip_reg_exclude and (clip_reg->call_success_flag or clip_reg->tra_rescue_success_flag) and clip_reg->sv_type==VAR_TRA){
				// check left part
				if(clip_reg->leftClipRegNum==1){
					if(clip_reg->leftClipReg){
						chrname_tmp = clip_reg->leftClipReg->chrname;
						start_pos = clip_reg->leftClipReg->startRefPos;
						end_pos = clip_reg->leftClipReg->endRefPos;
					}else{
						chrname_tmp = clip_reg->leftClipReg2->chrname;
						start_pos = clip_reg->leftClipReg2->startRefPos;
						end_pos = clip_reg->leftClipReg2->endRefPos;
					}
				}else if(clip_reg->leftClipRegNum==2){
					chrname_tmp = clip_reg->leftClipReg->chrname;
					start_pos = clip_reg->leftClipReg->startRefPos;
					end_pos = clip_reg->leftClipReg2->endRefPos;
				}else{
					cerr << __func__ << ", line=" << __LINE__ << ": invalid left region count: " << clip_reg->leftClipRegNum << ", error!" << endl;
					exit(1);
				}

				if(var_cand->chrname.compare(chrname_tmp)==0){
					overlap_flag = isOverlappedPos(var_cand->leftClipRefPos, var_cand->rightClipRefPos, start_pos, end_pos);
					if(overlap_flag){
						success_flag = true;
						break;
					}
				}

				// check right part
				if(clip_reg->rightClipRegNum==1){
					if(clip_reg->rightClipReg){
						chrname_tmp = clip_reg->rightClipReg->chrname;
						start_pos = clip_reg->rightClipReg->startRefPos;
						end_pos = clip_reg->rightClipReg->endRefPos;
					}else{
						chrname_tmp = clip_reg->rightClipReg2->chrname;
						start_pos = clip_reg->rightClipReg2->startRefPos;
						end_pos = clip_reg->rightClipReg2->endRefPos;
					}
				}else if(clip_reg->rightClipRegNum==2){
					chrname_tmp = clip_reg->rightClipReg->chrname;
					start_pos = clip_reg->rightClipReg->startRefPos;
					end_pos = clip_reg->rightClipReg2->endRefPos;
				}else{
					cerr << __func__ << ", line=" << __LINE__ << ": invalid right region count: " << clip_reg->rightClipRegNum << ", error!" << endl;
					exit(1);
				}

				if(var_cand->chrname.compare(chrname_tmp)==0){
					overlap_flag = isOverlappedPos(var_cand->leftClipRefPos, var_cand->rightClipRefPos, start_pos, end_pos);
					if(overlap_flag){
						success_flag = true;
						break;
					}
				}
			}
		}
		if(success_flag) break;
	}

	return success_flag;
}

// call TRA according to mate clip regions
void Genome::genomeCallTra(){
	Time time;
	vector<blatAlnTra*> blat_aln_tra_vec;

	cout << "[" << time.getTime() << "]: call genomic translocations ..." << endl;

	generateBlatAlnFilenameTra();
	blat_aln_tra_vec = loadBlatAlnDataTra();

	if(paras->num_threads<=1){
		cout << "Blat alignment using single thread ..." << endl;
		blatAlnTra_st(&blat_aln_tra_vec);
		cout << "Blat alignment finished." << endl;
	}else{
		// blat aln using multiple threads
		cout << "Blat alignment using " << paras->num_threads << " threads ..." << endl;
		blatAlnTra_mt(&blat_aln_tra_vec);
		cout << "Blat alignment finished." << endl;
	}
	releaseBlatAlnDataTra(blat_aln_tra_vec);

	// call translocations
	cout << "calling TRAs ..." << endl;
	genomeCallTraOp();

	// finalize TRA calling
	finalizeCallTra();
}

// call TRA according to mate clip regions
void Genome::generateBlatAlnFilenameTra(){
	size_t i, j, k, t, round_num;
	Chrome *chr;
	vector<mateClipReg_t*> mate_clipReg_vec;
	mateClipReg_t *clip_reg;
	varCand *var_cand, *var_cand_tmp;
	blat_aln_t *item_blat;
	vector<int32_t> tra_loc_vec;
	ofstream aln_file_tra;
	string aln_filename_tra;

	aln_filename_tra = out_dir_tra + "/" + "aln_filename_tra";
	aln_file_tra.open(aln_filename_tra);
	if(!aln_file_tra.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << aln_filename_tra << endl;
		exit(1);
	}

	for(i=0; i<chromeVector.size(); i++){
		chr = chromeVector.at(i);
		//if(chr->decoy_flag==false and paras->include_decoy){
		if(chr->decoy_flag==false){
			mate_clipReg_vec = chr->mateClipRegVector;
			for(j=0; j<mate_clipReg_vec.size(); j++){
				clip_reg = chr->mateClipRegVector.at(j);
				clip_reg->call_success_flag = false;
				if(clip_reg->valid_flag and clip_reg->sv_type==VAR_TRA){ // valid TRA
					for(round_num=0; round_num<2; round_num++){
						// construct var_cand
						if(round_num==0){
							var_cand = clip_reg->left_var_cand_tra;
							var_cand_tmp = constructNewVarCand(clip_reg->left_var_cand_tra, clip_reg->right_var_cand_tra);
						}else if(round_num==1){
							var_cand = clip_reg->right_var_cand_tra;
							var_cand_tmp = constructNewVarCand(clip_reg->right_var_cand_tra, clip_reg->left_var_cand_tra);
						}

						if(var_cand and var_cand_tmp){
							// save to file
							aln_file_tra << var_cand->alnfilename << "\t" << var_cand->ctgfilename << "\t" << var_cand->refseqfilename << endl;
							aln_file_tra << var_cand_tmp->alnfilename << "\t" << var_cand_tmp->ctgfilename << "\t" << var_cand_tmp->refseqfilename << endl;

							// release memory
							for(k=0; k<var_cand_tmp->blat_aln_vec.size(); k++){ // free blat_aln_vec
								item_blat = var_cand_tmp->blat_aln_vec.at(k);
								for(t=0; t<item_blat->aln_segs.size(); t++) // free aln_segs
									delete item_blat->aln_segs.at(t);
								delete item_blat;
							}
							delete var_cand_tmp;
						}else{
							//cout << "hhhhhhhhhhhh" << endl;
						}
					}
				}
			}
		}
	}
	aln_file_tra.close();

	//cout << "File generated." << endl;
}

// call TRA according to mate clip regions
void Genome::genomeCallTraOp(){
	size_t i, j, k, t, round_num, success_num, total, fail_reason_id1, fail_reason_id2; // 0- unknown; 1- align failure
	Chrome *chr;
	vector<mateClipReg_t*> mate_clipReg_vec;
	mateClipReg_t *clip_reg;
	varCand *var_cand, *var_cand_tmp;
	blat_aln_t *item_blat;
	vector<int32_t> tra_loc_vec;

	total = 0;
	for(i=0; i<chromeVector.size(); i++){
		fail_reason_id1 = fail_reason_id2 = 0;
		success_num = 0;
		chr = chromeVector.at(i);
		mate_clipReg_vec = chr->mateClipRegVector;
		for(j=0; j<mate_clipReg_vec.size(); j++){
			clip_reg = chr->mateClipRegVector.at(j);

//			if(i>=5 and j>=30)
//				cout << "gggggggggggggggggggg, i=" << i << ", j=" << j << ", " << chr->chrname << endl;
//			else continue;

			clip_reg->call_success_flag = false;
			if(clip_reg->valid_flag and clip_reg->sv_type==VAR_TRA){ // valid TRA
				for(round_num=0; round_num<2; round_num++){
					// construct var_cand
					if(round_num==0){
						var_cand = clip_reg->left_var_cand_tra;
						var_cand_tmp = constructNewVarCand(clip_reg->left_var_cand_tra, clip_reg->right_var_cand_tra);
					}else if(round_num==1){
						var_cand = clip_reg->right_var_cand_tra;
						var_cand_tmp = constructNewVarCand(clip_reg->right_var_cand_tra, clip_reg->left_var_cand_tra);
					}

					if(var_cand and var_cand_tmp){
						// load align data
						var_cand->loadBlatAlnData();
						var_cand_tmp->loadBlatAlnData();

						// compute TRA locations
						tra_loc_vec = computeTraLoc(var_cand, var_cand_tmp, clip_reg, round_num);

						if(tra_loc_vec.size()>0){
							saveTraLoc2ClipReg(clip_reg, tra_loc_vec, var_cand, var_cand_tmp, round_num); // save TRA location to clip region
							clip_reg->call_success_flag = true;
						}else{ // set fail reason
							if(var_cand->align_success==false or var_cand_tmp->align_success==false){
								if(round_num==0) fail_reason_id1 = 1;
								else fail_reason_id2  = 1;
							}
						}

						// release memory
						for(k=0; k<var_cand_tmp->blat_aln_vec.size(); k++){ // free blat_aln_vec
							item_blat = var_cand_tmp->blat_aln_vec.at(k);
							for(t=0; t<item_blat->aln_segs.size(); t++) // free aln_segs
								delete item_blat->aln_segs.at(t);
							delete item_blat;
						}
						delete var_cand_tmp;

						if(clip_reg->call_success_flag) break;
					}else{
						//cout << "hhhhhhhhhhhh" << endl;
					}
				}
			}
			if(clip_reg->call_success_flag){
				//cout << "chr " << i << ", j=" << j << ", TRA success" << endl;
				success_num ++;
				total ++;
			}else{
				if(fail_reason_id1==1 or fail_reason_id2==1){
					//cout << clip_reg->leftClipPosNum << ", " << clip_reg->leftClipPosNum2 << ", " << clip_reg->rightClipPosNum << ", " << clip_reg->rightClipPosNum2 << endl;
					saveTraLoc2ClipRegForAlnFailure(clip_reg);
				}
			}
		}
		//cout << "chr " << i << ": success_num=" << success_num << endl;
	}
	//cout << "Genome: total=" << total << endl;
}

// finalize TRA calling
void Genome::finalizeCallTra(){
	Time time;

	cout << "1111111111111" << endl;
	cout << "[" << time.getTime() << "]: remove redundant TRAs... " << endl;
	removeRedundantTra();

	cout << "2222222222222" << endl;
	cout << "[" << time.getTime() << "]: merge close range TRAs... " << endl;
	mergeCloseRangeTra();

	// find indels back for mate clip regions
	cout << "3333333333333" << endl;
	cout << "[" << time.getTime() << "]: recall indels from TRAs... " << endl;
	recallIndelsFromTRA();
}

vector<blatAlnTra*> Genome::loadBlatAlnDataTra(){
	ifstream aln_file_tra;
	string aln_filename_tra, line;
	vector<blatAlnTra*> blat_aln_tra_vec;
	blatAlnTra *blat_aln_tra;
	vector<string> str_vec;

	aln_filename_tra = out_dir_tra + "/" + "aln_filename_tra";
	aln_file_tra.open(aln_filename_tra);
	if(!aln_file_tra.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << aln_filename_tra << endl;
		exit(1);
	}

	while(getline(aln_file_tra, line))
	{
		if(line.size()){
			str_vec = split(line, "\t");
			blat_aln_tra = new blatAlnTra(str_vec.at(0), str_vec.at(1), str_vec.at(2), paras->max_proc_running_minutes_call);
			blat_aln_tra->killed_blat_work_vec = &paras->killed_blat_work_vec;
			blat_aln_tra->killed_blat_work_file = &paras->killed_blat_work_file;
			blat_aln_tra->mtx_killed_blat_work = &paras->mtx_killed_blat_work;
			blat_aln_tra_vec.push_back(blat_aln_tra);
		}
	}
	blat_aln_tra_vec.shrink_to_fit();
	aln_file_tra.close();

	return blat_aln_tra_vec;
}

void Genome::releaseBlatAlnDataTra(vector<blatAlnTra*> &blat_aln_tra_vec){
	blatAlnTra *blat_aln_tra;
	for(size_t i=0; i<blat_aln_tra_vec.size(); i++){
		blat_aln_tra = blat_aln_tra_vec.at(i);
		delete blat_aln_tra;
	}
	vector<blatAlnTra*>().swap(blat_aln_tra_vec);
}

void Genome::blatAlnTra_st(vector<blatAlnTra*> *blat_aln_tra_vec){
	blatAlnTra *blat_aln_tra;
	for (size_t i=0; i < blat_aln_tra_vec->size(); i++){
		blat_aln_tra = blat_aln_tra_vec->at(i);
		//cout << i << ": " << blat_aln_tra->alnfilename << endl;
		blat_aln_tra->generateBlatResult();
	}
}

void Genome::blatAlnTra_mt(vector<blatAlnTra*> *blat_aln_tra_vec){
	int32_t i;
	MultiThread mt[paras->num_threads];
	for(i=0; i<paras->num_threads; i++){
		mt[i].setNumThreads(paras->num_threads);
		mt[i].setBlatAlnTraVec(blat_aln_tra_vec);
		mt[i].setUserThreadID(i);
		if(!mt[i].startBlatAlnTra()){
			cerr << __func__ << ", line=" << __LINE__ << ": unable to create thread, error!" << endl;
			exit(1);
		}
	}
	for(i=0; i<paras->num_threads; i++){
		if(!mt[i].join()){
			cerr << __func__ << ", line=" << __LINE__ << ": unable to join, error!" << endl;
			exit(1);
		}
	}
}

varCand* Genome::constructNewVarCand(varCand *var_cand, varCand *var_cand_tmp){
	varCand *var_cand_new = NULL;
	if(var_cand and var_cand_tmp){
		var_cand_new = new varCand();
		var_cand_new->chrname = var_cand_tmp->chrname;
		var_cand_new->var_cand_filename = "";
		var_cand_new->out_dir_call = out_dir_call;
		var_cand_new->misAln_filename = "";
		var_cand_new->inBamFile = paras->inBamFile;
		var_cand_new->fai = fai;
		var_cand_new->technology = paras->technology;

		var_cand_new->refseqfilename = var_cand_tmp->refseqfilename;  // refseq file name
		var_cand_new->ctgfilename = var_cand->ctgfilename;  // contig file name
		var_cand_new->readsfilename = var_cand->readsfilename;  // reads file name
		var_cand_new->ref_left_shift_size = var_cand_tmp->ref_left_shift_size;  // ref_left_shift_size
		var_cand_new->ref_right_shift_size = var_cand_tmp->ref_right_shift_size;  // ref_right_shift_size

		var_cand_new->blat_aligned_info_vec = NULL;
		var_cand_new->blat_var_cand_file = NULL;

		var_cand_new->minimap2_aligned_info_vec = NULL;
		var_cand_new->minimap2_var_cand_file = NULL;

		var_cand_new->min_sv_size = paras->min_sv_size_usr;
		var_cand_new->minReadsNumSupportSV = paras->minReadsNumSupportSV;
		var_cand_new->minClipEndSize = paras->minClipEndSize;
		var_cand_new->minConReadLen = paras->minConReadLen;
		var_cand_new->min_identity_match = paras->min_identity_match;
		var_cand_new->min_identity_merge = paras->min_identity_merge;
		var_cand_new->min_distance_merge = paras->min_distance_merge;
		var_cand_new->max_seg_num_per_read = paras->max_seg_num_per_read;
		var_cand_new->minMapQ = paras->minMapQ;
		var_cand_new->minHighMapQ = paras->minHighMapQ;

		var_cand_new->cns_success = var_cand->cns_success;
		var_cand_new->ctg_num = var_cand->ctg_num;

		for(size_t k=0; k<var_cand_tmp->varVec.size(); k++) var_cand_new->varVec.push_back(var_cand_tmp->varVec.at(k));
		var_cand_new->varVec.shrink_to_fit();

		var_cand_new->alnfilename = var_cand->alnfilename.substr(0, var_cand->alnfilename.size()-5) + "_" + var_cand_tmp->chrname + "_" + to_string(var_cand_tmp->leftClipRefPos) + "-" + to_string(var_cand_tmp->rightClipRefPos) + ".sim4";
		var_cand_new->alnfilename = preprocessPipeChar(var_cand_new->alnfilename);
		var_cand_new->align_success = false;
		var_cand_new->clip_reg_flag = true;

		var_cand_new->leftClipRefPos = var_cand_tmp->leftClipRefPos;
		var_cand_new->rightClipRefPos = var_cand_tmp->rightClipRefPos;
		var_cand_new->sv_type = var_cand_tmp->sv_type;
		var_cand_new->dup_num = var_cand_tmp->dup_num;

		//set genotyping parameters
		var_cand_tmp->setGtParas(paras->gt_min_sig_size, paras->gt_size_ratio_match, paras->gt_min_identity_merge, paras->gt_homo_ratio, paras->gt_hete_ratio, paras->minReadsNumSupportSV);

		// process monitor killed blat work
		var_cand_tmp->max_proc_running_minutes = paras->max_proc_running_minutes_call;
		var_cand_tmp->killed_blat_work_vec = &paras->killed_blat_work_vec;
		var_cand_tmp->killed_blat_work_file = &paras->killed_blat_work_file;
		var_cand_tmp->mtx_killed_blat_work = &paras->mtx_killed_blat_work;

		// process monitor killed minimap2 work
		//var_cand_tmp->max_proc_running_minutes = paras->max_proc_running_minutes_call;
		var_cand_tmp->killed_minimap2_work_vec = &paras->killed_minimap2_work_vec;
		var_cand_tmp->killed_minimap2_work_file = &paras->killed_minimap2_work_file;
		var_cand_tmp->mtx_killed_minimap2_work = &paras->mtx_killed_minimap2_work;
	}

	return var_cand_new;
}

// compute variant location
vector<int32_t> Genome::computeTraLoc(varCand *var_cand, varCand *var_cand_tmp, mateClipReg_t *clip_reg, int32_t round_num){
	size_t i, j, aln_orient, missing_size_head, missing_size_inner, missing_size_tail;
	int32_t start_query_pos, end_query_pos;
	int32_t maxValue, maxIdx, secValue, query_dist, query_part_flag_maxValue, query_part_flag_secValue;  // query_part_flag: 0-head, 1-inner, 2-tail
	int32_t tra_blat_aln_idx, start_tra_aln_seg_idx, end_tra_aln_seg_idx, value_tmp, query_part_flag;
	blat_aln_t *blat_aln, *blat_aln_mate;
	aln_seg_t *aln_seg1, *aln_seg2;
	vector<int32_t> aln_idx_vec, mate_aln_idx_vec, tra_loc_vec, tra_loc_vec2;
	bool call_success_flag, same_orient_flag;
	int64_t idx_arr[4], tra_pos_arr[4], minIdx, minValue, dist, tmp_pos;

	call_success_flag = false;
	for(i=0; i<var_cand->blat_aln_vec.size(); i++){
		blat_aln = var_cand->blat_aln_vec.at(i);
		if(blat_aln->valid_aln){
			aln_orient = blat_aln->aln_orient;

			// query head
			if(aln_orient==ALN_PLUS_ORIENT){
				if(blat_aln->aln_segs.at(0)->subject_start>MIN_CLIP_DIST_THRES)
					missing_size_head = blat_aln->aln_segs.at(0)->query_start - 1;
				else
					missing_size_head = 0;
			}else{
				if(blat_aln->subject_len-blat_aln->aln_segs.at(blat_aln->aln_segs.size()-1)->subject_end>MIN_CLIP_DIST_THRES)
					missing_size_head = blat_aln->query_len - blat_aln->aln_segs.at(blat_aln->aln_segs.size()-1)->query_end;
				else
					missing_size_head = 0;
			}

			// query inner part
			maxValue = 0; maxIdx = -1;
			for(j=0; j<blat_aln->aln_segs.size()-1; j++){
				aln_seg1 = blat_aln->aln_segs.at(j);
				aln_seg2 = blat_aln->aln_segs.at(j+1);
				query_dist = aln_seg2->query_start - aln_seg1->query_end;
				if(query_dist>maxValue) {
					maxValue = query_dist;
					maxIdx = j;
				}
			}
			missing_size_inner = maxValue;

			// query tail
			if(aln_orient==ALN_PLUS_ORIENT){
				if(blat_aln->subject_len-blat_aln->aln_segs.at(blat_aln->aln_segs.size()-1)->subject_end>MIN_CLIP_DIST_THRES)
					missing_size_tail = blat_aln->query_len - blat_aln->aln_segs.at(blat_aln->aln_segs.size()-1)->query_end;
				else
					missing_size_tail = 0;
			}else{
				if(blat_aln->aln_segs.at(0)->subject_start>MIN_CLIP_DIST_THRES)
					missing_size_tail = blat_aln->aln_segs.at(0)->query_start - 1;
				else
					missing_size_tail = 0;
			}

			// choose the maximum and the second maximum
			query_part_flag_maxValue = 0; maxValue = missing_size_head; secValue = 0; query_part_flag_secValue = -1;
			if(maxValue<(int32_t)missing_size_inner) { query_part_flag_secValue = query_part_flag_maxValue; secValue = maxValue; query_part_flag_maxValue = 1; maxValue = missing_size_inner; }
			else if(secValue<(int32_t)missing_size_inner) { query_part_flag_secValue = 1; secValue = missing_size_inner; }
			if(maxValue<(int32_t)missing_size_tail) { query_part_flag_secValue = query_part_flag_maxValue; secValue = maxValue; query_part_flag_maxValue = 2; maxValue = missing_size_tail; }
			else if(secValue<(int32_t)missing_size_tail) { query_part_flag_secValue = 2; secValue = missing_size_tail; }

			for(j=0; j<2; j++){
				if(j==0){
					value_tmp = maxValue;
					query_part_flag = query_part_flag_maxValue;
				}else{
					value_tmp = secValue;
					query_part_flag = query_part_flag_secValue;
				}

				if(value_tmp>MIN_CLIP_DIST_THRES){
					tra_blat_aln_idx = i; start_query_pos = end_query_pos = 0;
					if(query_part_flag==0){ // query head
						start_query_pos = 1;
						end_query_pos = missing_size_head;
						start_tra_aln_seg_idx = -1;
						if(aln_orient==ALN_PLUS_ORIENT)  end_tra_aln_seg_idx = 0;
						else end_tra_aln_seg_idx = blat_aln->aln_segs.size() - 1;

						//cout << "head: " << missing_size_head << ", " << var_cand->chrname << ":" << start_query_pos << "-" << end_query_pos << endl;
					}else if(query_part_flag==1){ // query inner
						aln_seg1 = blat_aln->aln_segs.at(maxIdx);
						aln_seg2 = blat_aln->aln_segs.at(maxIdx+1);
						if(aln_orient==ALN_PLUS_ORIENT){ // plus orient
							start_query_pos = aln_seg1->query_end + 1;
							end_query_pos = aln_seg2->query_start - 1;
						}else{ // minus orient
							start_query_pos = blat_aln->query_len - (aln_seg2->query_start - 1) + 1;
							end_query_pos = blat_aln->query_len - (aln_seg1->query_end + 1) + 1;
						}
						if(aln_orient==ALN_PLUS_ORIENT){
							start_tra_aln_seg_idx = maxIdx;
							end_tra_aln_seg_idx = maxIdx + 1;
						}else{
							start_tra_aln_seg_idx = maxIdx + 1;
							end_tra_aln_seg_idx = maxIdx;
						}

						//cout << "inner: " << missing_size_inner << ", " << var_cand->chrname << ":" << start_query_pos << "-" << end_query_pos << endl;
					}else if(query_part_flag==2){ // query tail
						start_query_pos = blat_aln->query_len - missing_size_tail + 1;
						end_query_pos = blat_aln->query_len;

						if(aln_orient==ALN_PLUS_ORIENT) start_tra_aln_seg_idx = blat_aln->aln_segs.size() - 1;
						else start_tra_aln_seg_idx = 0;
						end_tra_aln_seg_idx = -1;

						//cout << "tail: " << missing_size_tail << ", " << var_cand->chrname << ":" << start_query_pos << "-" << end_query_pos << endl;
					}
					aln_idx_vec.push_back(tra_blat_aln_idx);
					aln_idx_vec.push_back(start_tra_aln_seg_idx);
					aln_idx_vec.push_back(end_tra_aln_seg_idx);

					// get the missing region in the mated clipping region
					mate_aln_idx_vec = getMateTraReg(blat_aln->query_id, start_query_pos, end_query_pos, var_cand_tmp);
					if(mate_aln_idx_vec.at(0)!=-1){
						blat_aln_mate = var_cand_tmp->blat_aln_vec.at(mate_aln_idx_vec.at(0));
//						aln_seg1 = blat_aln_mate->aln_segs.at(mate_aln_idx_vec.at(1));
//						aln_seg2 = blat_aln_mate->aln_segs.at(mate_aln_idx_vec.at(2));
//						if(blat_aln_mate->aln_orient==ALN_PLUS_ORIENT){ // plus orient
//							start_mate_query_pos = aln_seg1->query_start;
//							end_mate_query_pos = aln_seg2->query_end;
//						}else{ // minus orient
//							start_mate_query_pos = blat_aln_mate->query_len - aln_seg1->query_end + 1;
//							end_mate_query_pos = blat_aln_mate->query_len - aln_seg2->query_start + 1;
//						}

						//cout << "mate TRA: " << "start_mate_query_pos=" << start_mate_query_pos << ", end_mate_query_pos=" << end_mate_query_pos << endl;

						// compute TRA clipping locations
						tra_loc_vec = computeTraClippingLoc(query_part_flag, var_cand, aln_idx_vec, var_cand_tmp, mate_aln_idx_vec);

						//cout << "TRA location: " << var_cand->chrname << ":" << tra_loc_vec.at(0) << "-" << tra_loc_vec.at(1) << ", " << var_cand_tmp->chrname << ":" << tra_loc_vec.at(2) << "-" << tra_loc_vec.at(3) << endl;

						if(blat_aln->aln_orient==blat_aln_mate->aln_orient) same_orient_flag = true;
						else same_orient_flag = false;
						call_success_flag = computeTraCallSuccessFlag(tra_loc_vec, tra_loc_vec2, var_cand, var_cand_tmp, clip_reg, same_orient_flag, round_num);

						if(call_success_flag) break;
						else if(tra_loc_vec.size()>0){
							tra_loc_vec2.push_back(tra_loc_vec.at(0));
							tra_loc_vec2.push_back(tra_loc_vec.at(1));
							tra_loc_vec2.push_back(tra_loc_vec.at(2));
							tra_loc_vec2.push_back(tra_loc_vec.at(3));
							aln_idx_vec.clear();
							mate_aln_idx_vec.clear();
							tra_loc_vec.clear();
						}
					}
					aln_idx_vec.clear();
				}
			}
			if(call_success_flag) break;
		}
	}

	// adjust the item order
	if(tra_loc_vec.size()){
		if(round_num==0){
			tra_pos_arr[0] = clip_reg->leftMeanClipPos>1 ? clip_reg->leftMeanClipPos-1 : 0;   // 0 for null value
			tra_pos_arr[1] = clip_reg->leftMeanClipPos2>1 ? clip_reg->leftMeanClipPos2+1 : 0;
			tra_pos_arr[2] = clip_reg->rightMeanClipPos>1 ? clip_reg->rightMeanClipPos-1 : 0;
			tra_pos_arr[3] = clip_reg->rightMeanClipPos2>1 ? clip_reg->rightMeanClipPos2+1 : 0;
		}else{
			tra_pos_arr[0] = clip_reg->rightMeanClipPos>1 ? clip_reg->rightMeanClipPos-1 : 0;   // 0 for null value
			tra_pos_arr[1] = clip_reg->rightMeanClipPos2>1 ? clip_reg->rightMeanClipPos2+1 : 0;
			tra_pos_arr[2] = clip_reg->leftMeanClipPos>1 ? clip_reg->leftMeanClipPos-1 : 0;
			tra_pos_arr[3] = clip_reg->leftMeanClipPos2>1 ? clip_reg->leftMeanClipPos2+1 : 0;
		}

		// left part
		for(i=0; i<2; i++){
			idx_arr[i] = -1;
			tmp_pos = tra_loc_vec.at(i);
			if(tmp_pos!=-1){
				minIdx = -1; minValue = INT_MAX;
				for(j=0; j<2; j++){
					if(tra_pos_arr[j]!=0){
						dist = tmp_pos - tra_pos_arr[j];
						if(dist<0) dist = -dist;
						if(minValue>dist) {
							minValue = dist;
							minIdx = j;
						}
					}
				}
				idx_arr[i] = minIdx;
			}
		}
		// right part
		for(i=2; i<4; i++){
			idx_arr[i] = -1;
			tmp_pos = tra_loc_vec.at(i);
			if(tmp_pos!=-1){
				minIdx = -1; minValue = INT_MAX;
				for(j=2; j<4; j++){
					if(tra_pos_arr[j]!=0){
						dist = tmp_pos - tra_pos_arr[j];
						if(dist<0) dist = -dist;
						if(minValue>dist) {
							minValue = dist;
							minIdx = j;
						}
					}
				}
				idx_arr[i] = minIdx;
			}
		}

		for(i=0; i<4; i++){
			if(idx_arr[i]!=-1 and idx_arr[i]!=(int64_t)i){ // swap
				tmp_pos = tra_loc_vec.at(i);
				tra_loc_vec.at(i) = tra_loc_vec.at(idx_arr[i]);
				tra_loc_vec.at(idx_arr[i]) = tmp_pos;
			}
		}
	}

	return tra_loc_vec;
}

// get the missing region in the mated clipping region: [0]-blat_aln_id, [1]-start_aln_seg_id, [2]-end_aln_seg_id
vector<int32_t> Genome::getMateTraReg(int32_t query_id, int32_t start_query_pos, int32_t end_query_pos, varCand *var_cand_tmp){
	vector<int32_t> mate_idx_vec;
	size_t i, aln_orient, start_pos_tmp, end_pos_tmp, start_mate_query_pos, end_mate_query_pos;
	int32_t j, blat_aln_idx, start_seg_idx, end_seg_idx;
	blat_aln_t *blat_aln;
	aln_seg_t *aln_seg;
	bool valid_tra_flag;

	blat_aln_idx = start_seg_idx = end_seg_idx = -1;
	for(i=0; i<var_cand_tmp->blat_aln_vec.size(); i++){
		blat_aln_idx = start_seg_idx = end_seg_idx = -1;
		start_mate_query_pos = end_mate_query_pos = 0;

		blat_aln = var_cand_tmp->blat_aln_vec.at(i);
		if(blat_aln->valid_aln and blat_aln->query_id==query_id){
			aln_orient = blat_aln->aln_orient;
			if(aln_orient==ALN_PLUS_ORIENT){ // plus orient
				// compute the start pos
				for(j=0; j<(int32_t)blat_aln->aln_segs.size(); j++){
					aln_seg = blat_aln->aln_segs.at(j);
					start_pos_tmp = aln_seg->query_start;
					end_pos_tmp = aln_seg->query_end;
					if(isOverlappedPos(start_pos_tmp, end_pos_tmp, start_query_pos, end_query_pos)){
						start_mate_query_pos = start_pos_tmp;
						start_seg_idx = j;
						break;
					}
				}

				// compute the end pos
				if(start_seg_idx!=-1){
					for(j=(int32_t)blat_aln->aln_segs.size()-1; j>=start_seg_idx; j--){
						aln_seg = blat_aln->aln_segs.at(j);
						start_pos_tmp = aln_seg->query_start;
						end_pos_tmp = aln_seg->query_end;
						if(isOverlappedPos(start_pos_tmp, end_pos_tmp, start_query_pos, end_query_pos)){
							end_mate_query_pos = end_pos_tmp;
							end_seg_idx = j;
							break;
						}
					}
				}
			}else{ // minus orient
				// compute the start pos
				for(j=(int32_t)blat_aln->aln_segs.size()-1; j>=0; j--){
					aln_seg = blat_aln->aln_segs.at(j);
					start_pos_tmp = blat_aln->query_len - aln_seg->query_end + 1;
					end_pos_tmp = blat_aln->query_len - aln_seg->query_start + 1;
					if(isOverlappedPos(start_pos_tmp, end_pos_tmp, start_query_pos, end_query_pos)){
						start_mate_query_pos = start_pos_tmp;
						start_seg_idx = j;
						break;
					}
				}

				// compute the end pos
				if(start_seg_idx!=-1){
					for(j=0; j<(int32_t)blat_aln->aln_segs.size(); j++){
						aln_seg = blat_aln->aln_segs.at(j);
						start_pos_tmp = blat_aln->query_len - aln_seg->query_end + 1;
						end_pos_tmp = blat_aln->query_len - aln_seg->query_start + 1;
						if(isOverlappedPos(start_pos_tmp, end_pos_tmp, start_query_pos, end_query_pos)){
							end_mate_query_pos = end_pos_tmp;
							end_seg_idx = j;
							break;
						}
					}
				}
			}

			if(start_seg_idx!=-1 and end_seg_idx!=-1){
				valid_tra_flag = isValidTraReg(start_query_pos, end_query_pos, start_mate_query_pos, end_mate_query_pos);
				if(valid_tra_flag){
					blat_aln_idx = i;
					break;
				}
			}else
				blat_aln_idx = start_seg_idx = end_seg_idx = -1;
		}
	}

	mate_idx_vec.push_back(blat_aln_idx);
	mate_idx_vec.push_back(start_seg_idx);
	mate_idx_vec.push_back(end_seg_idx);

	return mate_idx_vec;
}

// determine whether the TRA region is valid
bool Genome::isValidTraReg(size_t start_query_pos, size_t end_query_pos, size_t start_mate_query_pos, size_t end_mate_query_pos){
	bool valid_flag;
	double shared_ratio, shared_len, total_len;
	size_t minPos, maxPos, min_shared_pos, max_shared_pos;

	// compute the ratio of shared region
	if(start_query_pos<start_mate_query_pos){
		minPos = start_query_pos;
		min_shared_pos = start_mate_query_pos;
	}else{
		minPos = start_mate_query_pos;
		min_shared_pos = start_query_pos;
	}
	if(end_query_pos<end_mate_query_pos){
		maxPos = end_mate_query_pos;
		max_shared_pos = end_query_pos;
	}else{
		maxPos = end_query_pos;
		max_shared_pos = end_mate_query_pos;
	}
	shared_len = max_shared_pos - min_shared_pos + 1;
	total_len = maxPos - minPos + 1;
	shared_ratio = shared_len / total_len;

	//cout << "shared_len=" << shared_len << ", total_len=" << total_len << ", shared_ratio=" << shared_ratio << endl;

	if(shared_ratio>MIN_VALID_TRA_RATIO) valid_flag = true;
	else valid_flag = false;

	return valid_flag;
}

// compute TRA clipping locations
vector<int32_t> Genome::computeTraClippingLoc(size_t query_clip_part_flag, varCand *var_cand, vector<int32_t> &aln_idx_vec, varCand *var_cand_tmp, vector<int32_t> &mate_aln_idx_vec){
	vector<int32_t> tra_loc_vec;
	int32_t left_tra_refPos1, right_tra_refPos1, left_tra_refPos2, right_tra_refPos2;
	blat_aln_t *blat_aln, *blat_aln_tmp;
	aln_seg_t *aln_seg1, *aln_seg2;

	left_tra_refPos1 = left_tra_refPos2 = right_tra_refPos1 = right_tra_refPos2 = -1;
	blat_aln = var_cand->blat_aln_vec.at(aln_idx_vec.at(0));
	blat_aln_tmp = var_cand_tmp->blat_aln_vec.at(mate_aln_idx_vec.at(0));

	if(query_clip_part_flag==0){ // query head
		aln_seg1 = blat_aln->aln_segs.at(aln_idx_vec.at(2));
		left_tra_refPos1 = -1;
		if(blat_aln->aln_orient==ALN_PLUS_ORIENT)
			left_tra_refPos2 = aln_seg1->ref_start;
		else
			left_tra_refPos2 = aln_seg1->ref_end;
	}else if(query_clip_part_flag==1){ // query inner
		aln_seg1 = blat_aln->aln_segs.at(aln_idx_vec.at(1));
		aln_seg2 = blat_aln->aln_segs.at(aln_idx_vec.at(2));
		if(blat_aln->aln_orient==ALN_PLUS_ORIENT){
			left_tra_refPos1 = aln_seg1->ref_end;
			left_tra_refPos2 = aln_seg2->ref_start;
		}else{
			left_tra_refPos1 = aln_seg1->ref_start;
			left_tra_refPos2 = aln_seg2->ref_end;
		}
	}else if(query_clip_part_flag==2){ // query tail
		aln_seg2 = blat_aln->aln_segs.at(aln_idx_vec.at(1));
		if(blat_aln->aln_orient==ALN_PLUS_ORIENT)
			left_tra_refPos1 = aln_seg2->ref_end;
		else
			left_tra_refPos1 = aln_seg2->ref_start;
		left_tra_refPos2 = -1;
	}
	//cout << "left_tra_refPos1=" << left_tra_refPos1 << ", left_tra_refPos2=" << left_tra_refPos2 << endl;

	if(query_clip_part_flag==0){ // query head
		aln_seg2 = blat_aln_tmp->aln_segs.at(mate_aln_idx_vec.at(2));
		right_tra_refPos1 = -1;
		if(blat_aln_tmp->aln_orient==ALN_PLUS_ORIENT){ // plus orient
			right_tra_refPos2 = aln_seg2->ref_end;
		}else{ // minus orient
			right_tra_refPos2 = aln_seg2->ref_start;
		}
	}else if(query_clip_part_flag==1){ // query inner
		aln_seg1 = blat_aln_tmp->aln_segs.at(mate_aln_idx_vec.at(1));
		aln_seg2 = blat_aln_tmp->aln_segs.at(mate_aln_idx_vec.at(2));
		if(blat_aln_tmp->aln_orient==ALN_PLUS_ORIENT){ // plus orient
			right_tra_refPos1 = aln_seg1->ref_start;
			right_tra_refPos2 = aln_seg2->ref_end;
		}else{ // minus orient
			right_tra_refPos1 = aln_seg1->ref_end;
			right_tra_refPos2 = aln_seg2->ref_start;
		}
	}else if(query_clip_part_flag==2){ // query tail
		aln_seg1 = blat_aln_tmp->aln_segs.at(mate_aln_idx_vec.at(1));
		if(blat_aln_tmp->aln_orient==ALN_PLUS_ORIENT){ // plus orient
			right_tra_refPos1 = aln_seg1->ref_start;
		}else{ // minus orient
			right_tra_refPos1 = aln_seg1->ref_end;
		}
		right_tra_refPos2 = -1;
	}
	//cout << "right_tra_refPos1=" << right_tra_refPos1 << ", right_tra_refPos2=" << right_tra_refPos2 << endl;

	tra_loc_vec.push_back(left_tra_refPos1);
	tra_loc_vec.push_back(left_tra_refPos2);
	tra_loc_vec.push_back(right_tra_refPos1);
	tra_loc_vec.push_back(right_tra_refPos2);

	return tra_loc_vec;
}


// determine whether the TRA is called successfully according to TRA location vector and mate cliping information
bool Genome::computeTraCallSuccessFlag(vector<int32_t> &tra_loc_vec, vector<int32_t> &tra_loc_vec2, varCand *var_cand, varCand *var_cand_tmp, mateClipReg_t *clip_reg, bool same_orient_flag, int32_t round_num){
	bool flag = false, overlap_flag1, overlap_flag2, overlap_flag3, overlap_flag4;
	size_t i, tmp;
	reg_t *reg1_clipReg, *reg2_clipReg;
	int32_t ref_dist1, ref_dist2, tra_loc1, tra_loc2, tra_loc3, tra_loc4, tra_loc1_clipReg, tra_loc2_clipReg, tra_loc3_clipReg, tra_loc4_clipReg;
	int32_t near_id1, near_id2, dist1, dist2;
	double ratio;

	if(tra_loc_vec.size()>0){
		if(clip_reg->leftClipRegNum==1 and clip_reg->rightClipRegNum==1){
			if(tra_loc_vec.at(0)!=-1 and tra_loc_vec.at(1)!=-1 and tra_loc_vec.at(2)!=-1 and tra_loc_vec.at(3)!=-1){ // four values in 'tra_loc_vec'
				if(tra_loc_vec.at(0)<tra_loc_vec.at(1)){
					tra_loc1 = tra_loc_vec.at(0);
					tra_loc2 = tra_loc_vec.at(1);
				}else{
					tra_loc1 = tra_loc_vec.at(1);
					tra_loc2 = tra_loc_vec.at(0);
				}
				if(tra_loc_vec.at(2)<tra_loc_vec.at(3)){
					tra_loc3 = tra_loc_vec.at(2);
					tra_loc4 = tra_loc_vec.at(3);
				}else{
					tra_loc3 = tra_loc_vec.at(3);
					tra_loc4 = tra_loc_vec.at(2);
				}
				reg1_clipReg = clip_reg->leftClipReg ? clip_reg->leftClipReg : clip_reg->leftClipReg2;
				reg2_clipReg = clip_reg->rightClipReg ? clip_reg->rightClipReg : clip_reg->rightClipReg2;

				overlap_flag1 = isOverlappedPos(tra_loc1, tra_loc2, reg1_clipReg->startRefPos, reg1_clipReg->endRefPos);
				overlap_flag2 = isOverlappedPos(tra_loc3, tra_loc4, reg2_clipReg->startRefPos, reg2_clipReg->endRefPos);
				overlap_flag3 = isOverlappedPos(tra_loc1, tra_loc2, reg2_clipReg->startRefPos, reg2_clipReg->endRefPos);
				overlap_flag4 = isOverlappedPos(tra_loc3, tra_loc4, reg1_clipReg->startRefPos, reg1_clipReg->endRefPos);

				if((overlap_flag1 and overlap_flag2) or (overlap_flag3 and overlap_flag4))
					flag = true;
			}else if((tra_loc_vec.at(0)!=-1 or tra_loc_vec.at(1)!=-1) and (tra_loc_vec.at(2)!=-1 or tra_loc_vec.at(3)!=-1)){ // missing some values in 'tra_loc_vec'
				tra_loc1 = tra_loc_vec.at(0)!=-1 ? tra_loc_vec.at(0) : tra_loc_vec.at(1);
				tra_loc2 = tra_loc_vec.at(2)!=-1 ? tra_loc_vec.at(2) : tra_loc_vec.at(3);
				reg1_clipReg = clip_reg->leftClipReg ? clip_reg->leftClipReg : clip_reg->leftClipReg2;
				reg2_clipReg = clip_reg->rightClipReg ? clip_reg->rightClipReg : clip_reg->rightClipReg2;
				if(tra_loc1!=-1 and tra_loc2!=-1 and reg1_clipReg and reg2_clipReg){
					if((tra_loc1>=reg1_clipReg->startRefPos and tra_loc1<=reg1_clipReg->endRefPos and tra_loc2>=reg2_clipReg->startRefPos and tra_loc2<=reg2_clipReg->endRefPos)
						or (tra_loc1>=reg2_clipReg->startRefPos and tra_loc1<=reg2_clipReg->endRefPos and tra_loc2>=reg1_clipReg->startRefPos and tra_loc2<=reg1_clipReg->endRefPos)){
						flag = true;
					}
				}
			}
		}else if(clip_reg->leftClipRegNum==2 and clip_reg->rightClipRegNum==2){
			if(tra_loc_vec.at(0)!=-1 and tra_loc_vec.at(1)!=-1 and tra_loc_vec.at(2)!=-1 and tra_loc_vec.at(3)!=-1){ // in correct region
				if(tra_loc_vec.at(0)<tra_loc_vec.at(1)){
					tra_loc1 = tra_loc_vec.at(0);
					tra_loc2 = tra_loc_vec.at(1);
				}else{
					tra_loc1 = tra_loc_vec.at(1);
					tra_loc2 = tra_loc_vec.at(0);
				}
				if(tra_loc_vec.at(2)<tra_loc_vec.at(3)){
					tra_loc3 = tra_loc_vec.at(2);
					tra_loc4 = tra_loc_vec.at(3);
				}else{
					tra_loc3 = tra_loc_vec.at(3);
					tra_loc4 = tra_loc_vec.at(2);
				}

				if(clip_reg->leftMeanClipPos<clip_reg->leftMeanClipPos2){
					tra_loc1_clipReg = clip_reg->leftMeanClipPos;
					tra_loc2_clipReg = clip_reg->leftMeanClipPos2;
				}else{
					tra_loc1_clipReg = clip_reg->leftMeanClipPos2;
					tra_loc2_clipReg = clip_reg->leftMeanClipPos;
				}
				if(clip_reg->rightMeanClipPos<clip_reg->rightMeanClipPos2){
					tra_loc3_clipReg = clip_reg->rightMeanClipPos;
					tra_loc4_clipReg = clip_reg->rightMeanClipPos2;
				}else{
					tra_loc3_clipReg = clip_reg->rightMeanClipPos2;
					tra_loc4_clipReg = clip_reg->rightMeanClipPos;
				}

				overlap_flag1 = isOverlappedPos(tra_loc1, tra_loc2, tra_loc1_clipReg, tra_loc2_clipReg);
				overlap_flag2 = isOverlappedPos(tra_loc3, tra_loc4, tra_loc3_clipReg, tra_loc4_clipReg);
				overlap_flag3 = isOverlappedPos(tra_loc1, tra_loc2, tra_loc3_clipReg, tra_loc4_clipReg);
				overlap_flag4 = isOverlappedPos(tra_loc3, tra_loc4, tra_loc1_clipReg, tra_loc2_clipReg);
				if((overlap_flag1 and overlap_flag2) or (overlap_flag3 and overlap_flag4))
					flag = true;
			}

			// reconstruct TRA locations
			if(flag==false and tra_loc_vec2.size()>0){
				for(i=0; i<tra_loc_vec.size(); i++)
					if(tra_loc_vec.at(i)==-1 and tra_loc_vec2.at(i)!=-1)
						tra_loc_vec.at(i) = tra_loc_vec2.at(i);
				if(tra_loc_vec.at(0)==-1 and tra_loc_vec2.at(0)==-1 and tra_loc_vec.at(1)!=-1 and tra_loc_vec2.at(1)!=-1)
					tra_loc_vec.at(0) = tra_loc_vec2.at(1);
				else if(tra_loc_vec.at(1)==-1 and tra_loc_vec2.at(1)==-1 and tra_loc_vec.at(0)!=-1 and tra_loc_vec2.at(0)!=-1)
					tra_loc_vec.at(1) = tra_loc_vec2.at(0);
				if(tra_loc_vec.at(2)==-1 and tra_loc_vec2.at(2)==-1 and tra_loc_vec.at(3)!=-1 and tra_loc_vec2.at(3)!=-1)
					tra_loc_vec.at(2) = tra_loc_vec2.at(3);
				else if(tra_loc_vec.at(3)==-1 and tra_loc_vec2.at(3)==-1 and tra_loc_vec.at(2)!=-1 and tra_loc_vec2.at(2)!=-1)
					tra_loc_vec.at(3) = tra_loc_vec2.at(2);

				if(tra_loc_vec.at(0)!=-1 and tra_loc_vec.at(1)!=-1 and tra_loc_vec.at(2)!=-1 and tra_loc_vec.at(3)!=-1)
					flag = true;
			}
		}

		// check region size
		if(flag){
			if(tra_loc_vec.at(0)!=-1 and tra_loc_vec.at(1)!=-1 and tra_loc_vec.at(2)!=-1 and tra_loc_vec.at(3)!=-1){
				ref_dist1 = tra_loc_vec.at(1) - tra_loc_vec.at(0) + 1;
				if(ref_dist1<0) ref_dist1 = -ref_dist1;
				ref_dist2 = tra_loc_vec.at(3) - tra_loc_vec.at(2) + 1;
				if(ref_dist2<0) ref_dist2 = -ref_dist2;

				if(ref_dist1==0 or ref_dist2==0) ratio = 0;
				else ratio = (double)ref_dist1 / ref_dist2;
				if(ratio<1.0-CLIP_DIFF_LEN_RATIO_SV or ratio>1.0+CLIP_DIFF_LEN_RATIO_SV){
					if(ref_dist1<ref_dist2) { // same left/right region
						if(round_num==0){
							dist1 = tra_loc_vec.at(0) - clip_reg->leftMeanClipPos;
							dist2 = tra_loc_vec.at(0) - clip_reg->leftMeanClipPos2;
							if(dist1<0) dist1 = -dist1;
							if(dist2<0) dist2 = -dist2;
							if(dist1<dist2) near_id1 = 0;
							else near_id1 = 1;

							dist1 = tra_loc_vec.at(1) - clip_reg->leftMeanClipPos;
							if(dist1<0) dist1 = -dist1;
							dist2 = tra_loc_vec.at(1) - clip_reg->leftMeanClipPos2;
							if(dist2<0) dist2 = -dist2;
							if(dist1<dist2) near_id2 = 0;
							else near_id2 = 1;

							if(near_id1==near_id2){
								if(near_id1==0) tra_loc_vec.at(1) = clip_reg->leftMeanClipPos2;
								else tra_loc_vec.at(0) = clip_reg->leftMeanClipPos;
							}
						}else{
							dist1 = tra_loc_vec.at(0) - clip_reg->rightMeanClipPos;
							dist2 = tra_loc_vec.at(0) - clip_reg->rightMeanClipPos2;
							if(dist1<0) dist1 = -dist1;
							if(dist2<0) dist2 = -dist2;
							if(dist1<dist2) near_id1 = 2;
							else near_id1 = 3;

							dist1 = tra_loc_vec.at(1) - clip_reg->rightMeanClipPos;
							if(dist1<0) dist1 = -dist1;
							dist2 = tra_loc_vec.at(1) - clip_reg->rightMeanClipPos2;
							if(dist2<0) dist2 = -dist2;
							if(dist1<dist2) near_id2 = 2;
							else near_id2 = 3;

							if(near_id1==near_id2){
								if(near_id1==2) tra_loc_vec.at(1) = clip_reg->rightMeanClipPos2;
								else tra_loc_vec.at(0) = clip_reg->rightMeanClipPos;
							}
						}
					}else{ // same right/left region
						if(round_num==0){
							dist1 = tra_loc_vec.at(2) - clip_reg->rightMeanClipPos;
							if(dist1<0) dist1 = -dist1;
							dist2 = tra_loc_vec.at(2) - clip_reg->rightMeanClipPos2;
							if(dist2<0) dist2 = -dist2;
							if(dist1<dist2) near_id1 = 2;
							else near_id1 = 3;

							dist1 = tra_loc_vec.at(3) - clip_reg->rightMeanClipPos;
							if(dist1<0) dist1 = -dist1;
							dist2 = tra_loc_vec.at(3) - clip_reg->rightMeanClipPos2;
							if(dist2<0) dist2 = -dist2;
							if(dist1<dist2) near_id2 = 2;
							else near_id2 = 3;

							if(near_id1==near_id2){
								if(near_id1==2) tra_loc_vec.at(3) = clip_reg->rightMeanClipPos2;
								else tra_loc_vec.at(2) = clip_reg->rightMeanClipPos;
							}
						}else{
							dist1 = tra_loc_vec.at(2) - clip_reg->leftMeanClipPos;
							if(dist1<0) dist1 = -dist1;
							dist2 = tra_loc_vec.at(2) - clip_reg->leftMeanClipPos2;
							if(dist2<0) dist2 = -dist2;
							if(dist1<dist2) near_id1 = 0;
							else near_id1 = 1;

							dist1 = tra_loc_vec.at(3) - clip_reg->leftMeanClipPos;
							if(dist1<0) dist1 = -dist1;
							dist2 = tra_loc_vec.at(3) - clip_reg->leftMeanClipPos2;
							if(dist2<0) dist2 = -dist2;
							if(dist1<dist2) near_id2 = 0;
							else near_id2 = 0;

							if(near_id1==near_id2){
								if(near_id1==0) tra_loc_vec.at(3) = clip_reg->leftMeanClipPos2;
								else tra_loc_vec.at(2) = clip_reg->leftMeanClipPos;
							}
						}
					}
				}

				if(tra_loc_vec.at(0)>tra_loc_vec.at(1)){
					tmp = tra_loc_vec.at(0);
					tra_loc_vec.at(0) = tra_loc_vec.at(1);
					tra_loc_vec.at(1) = tmp;
				}
				if(tra_loc_vec.at(2)>tra_loc_vec.at(3)){
					tmp = tra_loc_vec.at(2);
					tra_loc_vec.at(2) = tra_loc_vec.at(3);
					tra_loc_vec.at(3) = tmp;
				}
				ref_dist1 = tra_loc_vec.at(1) - tra_loc_vec.at(0) + 1;
				ref_dist2 = tra_loc_vec.at(3) - tra_loc_vec.at(2) + 1;
				if(ref_dist1==0 or ref_dist2==0) ratio = 0;
				else ratio = (double)ref_dist1 / ref_dist2;
				if(ratio<1.0-CLIP_DIFF_LEN_RATIO_SV or ratio>1.0+CLIP_DIFF_LEN_RATIO_SV)
					flag = false;
			}
		}
	}
	return flag;
}

// save TRA location to clip region
void Genome::saveTraLoc2ClipReg(mateClipReg_t *clip_reg, vector<int32_t> &tra_loc_vec, varCand *var_cand, varCand *var_cand_tmp, size_t round_num){
//	if(tra_loc_vec.at(0)!=-1 and tra_loc_vec.at(1)!=-1 and tra_loc_vec.at(0)>tra_loc_vec.at(1)){  // INV_TRA
//		if(round_num==0){
//			clip_reg->chrname_leftTra1 = var_cand->chrname;
//			clip_reg->leftClipPosTra1 = tra_loc_vec.at(1);
//			clip_reg->chrname_rightTra1 = var_cand->chrname;
//			clip_reg->rightClipPosTra1 = tra_loc_vec.at(0);
//		}else{
//			clip_reg->chrname_leftTra2 = var_cand->chrname;
//			clip_reg->leftClipPosTra2 = tra_loc_vec.at(1);
//			clip_reg->chrname_rightTra2 = var_cand->chrname;
//			clip_reg->rightClipPosTra2 = tra_loc_vec.at(0);
//		}
//	}else{  // TRA
		if(round_num==0){
			if(tra_loc_vec.at(0)!=-1) clip_reg->chrname_leftTra1 = var_cand->chrname;
			clip_reg->leftClipPosTra1 = tra_loc_vec.at(0);
			if(tra_loc_vec.at(1)!=-1)
				//clip_reg->chrname_rightTra1 = var_cand->chrname;
				clip_reg->chrname_leftTra2 = var_cand->chrname;
			//clip_reg->rightClipPosTra1 = tra_loc_vec.at(1);
			clip_reg->leftClipPosTra2 = tra_loc_vec.at(1);
		}else{
			if(tra_loc_vec.at(0)!=-1)
				//clip_reg->chrname_leftTra2 = var_cand->chrname;
				clip_reg->chrname_rightTra1 = var_cand->chrname;
			//clip_reg->leftClipPosTra2 = tra_loc_vec.at(0);
			clip_reg->rightClipPosTra1 = tra_loc_vec.at(0);
			if(tra_loc_vec.at(1)!=-1) clip_reg->chrname_rightTra2 = var_cand->chrname;
			clip_reg->rightClipPosTra2 = tra_loc_vec.at(1);
		}
//	}
//	if(tra_loc_vec.at(2)!=-1 and tra_loc_vec.at(3)!=-1 and tra_loc_vec.at(2)>tra_loc_vec.at(3)){  // INV_TRA
//		if(round_num==0){
//			clip_reg->chrname_leftTra2 = var_cand_tmp->chrname;
//			clip_reg->leftClipPosTra2 = tra_loc_vec.at(3);
//			clip_reg->chrname_rightTra2 = var_cand_tmp->chrname;
//			clip_reg->rightClipPosTra2 = tra_loc_vec.at(2);
//		}else{
//			clip_reg->chrname_leftTra1 = var_cand_tmp->chrname;
//			clip_reg->leftClipPosTra1 = tra_loc_vec.at(3);
//			clip_reg->chrname_rightTra1 = var_cand_tmp->chrname;
//			clip_reg->rightClipPosTra1 = tra_loc_vec.at(2);
//		}
//	}else{  // TRA
		if(round_num==0){
			if(tra_loc_vec.at(2)!=-1)
				//clip_reg->chrname_leftTra2 = var_cand_tmp->chrname;
				clip_reg->chrname_rightTra1 = var_cand_tmp->chrname;
			//clip_reg->leftClipPosTra2 = tra_loc_vec.at(2);
			clip_reg->rightClipPosTra1 = tra_loc_vec.at(2);
			if(tra_loc_vec.at(3)!=-1) clip_reg->chrname_rightTra2 = var_cand_tmp->chrname;
			clip_reg->rightClipPosTra2 = tra_loc_vec.at(3);
		}else{
			if(tra_loc_vec.at(2)!=-1) clip_reg->chrname_leftTra1 = var_cand_tmp->chrname;
			clip_reg->leftClipPosTra1 = tra_loc_vec.at(2);
			if(tra_loc_vec.at(3)!=-1)
				//clip_reg->chrname_rightTra1 = var_cand_tmp->chrname;
				clip_reg->chrname_leftTra2 = var_cand_tmp->chrname;
			//clip_reg->rightClipPosTra1 = tra_loc_vec.at(3);
			clip_reg->leftClipPosTra2 = tra_loc_vec.at(3);
		}
//	}
}

// save TRA location to clip region for alignment failures
void Genome::saveTraLoc2ClipRegForAlnFailure(mateClipReg_t *clip_reg){
	if(clip_reg->call_success_flag==false){
		if(clip_reg->leftClipRegNum==1 and clip_reg->rightClipRegNum==1 and ((clip_reg->leftClipPosNum>=paras->minReadsNumSupportSV or clip_reg->leftClipPosNum2>=paras->minReadsNumSupportSV) and (clip_reg->rightClipPosNum>=paras->minReadsNumSupportSV or clip_reg->rightClipPosNum2>=paras->minReadsNumSupportSV))){
			if(clip_reg->leftClipReg){
				clip_reg->chrname_leftTra1 = clip_reg->leftClipReg->chrname;
				clip_reg->leftClipPosTra1 = clip_reg->leftMeanClipPos;
			}else if(clip_reg->leftClipReg2){
				//clip_reg->chrname_rightTra1 = clip_reg->leftClipReg2->chrname;
				//clip_reg->rightClipPosTra1 = clip_reg->leftMeanClipPos2;
				clip_reg->chrname_leftTra2 = clip_reg->leftClipReg2->chrname;
				clip_reg->leftClipPosTra2 = clip_reg->leftMeanClipPos2;
			}
			if(clip_reg->rightClipReg){
				//clip_reg->chrname_leftTra2 = clip_reg->rightClipReg->chrname;
				//clip_reg->leftClipPosTra2 = clip_reg->rightMeanClipPos;
				clip_reg->chrname_rightTra1 = clip_reg->rightClipReg->chrname;
				clip_reg->rightClipPosTra1 = clip_reg->rightMeanClipPos;
			}else if(clip_reg->rightClipReg2){
				clip_reg->chrname_rightTra2 = clip_reg->rightClipReg2->chrname;
				clip_reg->rightClipPosTra2 = clip_reg->rightMeanClipPos2;
			}
			clip_reg->tra_rescue_success_flag = true;
		}else if(clip_reg->leftClipRegNum==2 and clip_reg->rightClipRegNum==2 and (clip_reg->leftClipPosNum>=paras->minReadsNumSupportSV and clip_reg->leftClipPosNum2>=paras->minReadsNumSupportSV and clip_reg->rightClipPosNum>=paras->minReadsNumSupportSV and clip_reg->rightClipPosNum2>=paras->minReadsNumSupportSV)){
			if(clip_reg->leftClipReg){
				clip_reg->chrname_leftTra1 = clip_reg->leftClipReg->chrname;
				clip_reg->leftClipPosTra1 = clip_reg->leftMeanClipPos;
			}
			if(clip_reg->leftClipReg2){
				//clip_reg->chrname_rightTra1 = clip_reg->leftClipReg2->chrname;
				//clip_reg->rightClipPosTra1 = clip_reg->leftMeanClipPos2;
				clip_reg->chrname_leftTra2 = clip_reg->leftClipReg2->chrname;
				clip_reg->leftClipPosTra2 = clip_reg->leftMeanClipPos2;
			}
			if(clip_reg->rightClipReg){
				//clip_reg->chrname_leftTra2 = clip_reg->rightClipReg->chrname;
				//clip_reg->leftClipPosTra2 = clip_reg->rightMeanClipPos;
				clip_reg->chrname_rightTra1 = clip_reg->rightClipReg->chrname;
				clip_reg->rightClipPosTra1 = clip_reg->rightMeanClipPos;
			}
			if(clip_reg->rightClipReg2){
				clip_reg->chrname_rightTra2 = clip_reg->rightClipReg2->chrname;
				clip_reg->rightClipPosTra2 = clip_reg->rightMeanClipPos2;
			}
			clip_reg->tra_rescue_success_flag = true;
		}
	}
}

// merge close range TRA
void Genome::mergeCloseRangeTra(){
	size_t i, j, pos1, pos2;
	string chrname1, chrname2;
	Chrome *chr;
	mateClipReg_t *clip_reg, *clip_reg_ret;
	vector<int32_t> dist_vec;

	for(i=0; i<chromeVector.size(); i++){
		chr = chromeVector.at(i);
		for(j=0; j<chr->mateClipRegVector.size(); j++){
			clip_reg = chr->mateClipRegVector.at(j);
			if(clip_reg->valid_flag and clip_reg->sv_type==VAR_TRA and clip_reg->leftClipRegNum==1 and clip_reg->rightClipRegNum==1){
				// get close range
				clip_reg_ret = getNearDistedClipReg(clip_reg, chromeVector);

				// merge near distance TRA
				if(clip_reg_ret){

					dist_vec = computeDistsTra(clip_reg, clip_reg_ret);

					if(dist_vec.at(0)<dist_vec.at(1)){
						if(clip_reg_ret->leftClipPosTra1!=-1) { chrname1 = clip_reg_ret->chrname_leftTra1; pos1 = clip_reg_ret->leftClipPosTra1; }
						else { chrname1 = clip_reg_ret->chrname_rightTra1; pos1 = clip_reg_ret->rightClipPosTra1; }
					}else{
						if(clip_reg_ret->leftClipPosTra2!=-1) { chrname1 = clip_reg_ret->chrname_leftTra2; pos1 = clip_reg_ret->leftClipPosTra2; }
						else { chrname1 = clip_reg_ret->chrname_rightTra2; pos1 = clip_reg_ret->rightClipPosTra2; }
					}
					if(clip_reg->leftClipPosTra1==-1){
						clip_reg->chrname_leftTra1 = chrname1;
						clip_reg->leftClipPosTra1 = pos1;
					}else{
						clip_reg->chrname_rightTra1 = chrname1;
						clip_reg->rightClipPosTra1 = pos1;
					}
					clip_reg->leftClipRegNum ++;

					if(dist_vec.at(2)>dist_vec.at(3)){
						if(clip_reg_ret->leftClipPosTra2!=-1) { chrname2 = clip_reg_ret->chrname_leftTra2; pos2 = clip_reg_ret->leftClipPosTra2; }
						else { chrname2 = clip_reg_ret->chrname_rightTra2; pos2 = clip_reg_ret->rightClipPosTra2; }
					}else{
						if(clip_reg_ret->leftClipPosTra1!=-1) { chrname2 = clip_reg_ret->chrname_leftTra1; pos2 = clip_reg_ret->leftClipPosTra1; }
						else { chrname2 = clip_reg_ret->chrname_rightTra1; pos2 = clip_reg_ret->rightClipPosTra1; }
					}
					if(clip_reg->leftClipPosTra2==-1){
						clip_reg->chrname_leftTra2 = chrname2;
						clip_reg->leftClipPosTra2 = pos2;
					}else{
						clip_reg->chrname_rightTra2 = chrname2;
						clip_reg->rightClipPosTra2 = pos2;
					}
					clip_reg->rightClipRegNum ++;
					clip_reg_ret->valid_flag = false;

					// exchange
					if(clip_reg->leftClipPosTra1>clip_reg->rightClipPosTra1){
						pos1 = clip_reg->leftClipPosTra1;
						clip_reg->leftClipPosTra1 = clip_reg->rightClipPosTra1;
						clip_reg->rightClipPosTra1 = pos1;
					}
					if(clip_reg->leftClipPosTra2>clip_reg->rightClipPosTra2){
						pos2 = clip_reg->leftClipPosTra2;
						clip_reg->leftClipPosTra2 = clip_reg->rightClipPosTra2;
						clip_reg->rightClipPosTra2 = pos2;
					}
				}
			}
		}
	}

	// remove invalid elements
	removeInvalidMateClipItem();
}

mateClipReg_t* Genome::getNearDistedClipReg(mateClipReg_t *clip_reg_given, vector<Chrome*> &chrome_vec){
	size_t i, j;
	Chrome *chr;
	mateClipReg_t *clip_reg, *clip_reg_ret = NULL;
	int32_t dist1, dist2, dist3, dist4, min_dist1, min_dist2;
	vector<int32_t> dist_vec;
	double len_ratio;

	for(i=0; i<chromeVector.size(); i++){
		chr = chromeVector.at(i);
		for(j=0; j<chr->mateClipRegVector.size(); j++){
			clip_reg = chr->mateClipRegVector.at(j);

			if(clip_reg!=clip_reg_given and clip_reg->valid_flag and clip_reg->reg_mated_flag and clip_reg->sv_type==clip_reg_given->sv_type and clip_reg->leftClipRegNum==1 and clip_reg->rightClipRegNum==1){

				dist_vec = computeDistsTra(clip_reg_given, clip_reg);
				dist1 = dist_vec.at(0);
				dist2 = dist_vec.at(1);
				dist3 = dist_vec.at(2);
				dist4 = dist_vec.at(3);

				// chose the minimum
				if(dist1<dist2) min_dist1 = dist1;
				else min_dist1 = dist2;

				if(dist3<dist4) min_dist2 = dist3;
				else min_dist2 = dist4;

				if(min_dist1<paras->maxVarRegSize and min_dist2<paras->maxVarRegSize){
					len_ratio = (double)min_dist1 / min_dist2;
					if(len_ratio>1-CLIP_DIFF_LEN_RATIO_SV and len_ratio<1+CLIP_DIFF_LEN_RATIO_SV){
						clip_reg_ret = clip_reg;
						break;
					}
				}
			}
		}
		if(clip_reg_ret) break;
	}

	return clip_reg_ret;
}

vector<int32_t> Genome::computeDistsTra(mateClipReg_t *clip_reg1, mateClipReg_t *clip_reg2){
	string chrname1, chrname2, chrname3, chrname4;
	int32_t pos1, pos2, pos3, pos4, dist1, dist2, dist3, dist4;
	vector<int32_t> dist_vec;

	dist1 = dist2 = dist3 = dist4 = INT_MAX;
	if(clip_reg1!=clip_reg2 and clip_reg1->valid_flag and clip_reg2->reg_mated_flag and clip_reg2->sv_type==clip_reg1->sv_type and clip_reg1->leftClipRegNum==clip_reg2->leftClipRegNum and clip_reg1->rightClipRegNum==clip_reg2->rightClipRegNum){

		chrname1 = chrname2 = ""; pos1 = pos2 = 0;
		if(clip_reg1->leftClipPosTra1!=-1) { chrname1 = clip_reg1->chrname_leftTra1; pos1 = clip_reg1->leftClipPosTra1; }
		else if(clip_reg1->rightClipPosTra1!=-1) { chrname1 = clip_reg1->chrname_rightTra1; pos1 = clip_reg1->rightClipPosTra1; }
		if(clip_reg1->leftClipPosTra2!=-1) { chrname2 = clip_reg1->chrname_leftTra2; pos2 = clip_reg1->leftClipPosTra2; }
		else if(clip_reg1->rightClipPosTra2!=-1) { chrname2 = clip_reg1->chrname_rightTra2; pos2 = clip_reg1->rightClipPosTra2; }

		chrname1 = chrname2 = ""; pos3 = pos4 = 0;
		if(clip_reg2->leftClipPosTra1!=-1) { chrname3 = clip_reg2->chrname_leftTra1; pos3 = clip_reg2->leftClipPosTra1; }
		else if(clip_reg2->rightClipPosTra1!=-1) { chrname3 = clip_reg2->chrname_rightTra1; pos3 = clip_reg2->rightClipPosTra1; }
		if(clip_reg2->leftClipPosTra2!=-1) { chrname4 = clip_reg2->chrname_leftTra2; pos4 = clip_reg2->leftClipPosTra2; }
		else if(clip_reg2->rightClipPosTra2!=-1) { chrname4 = clip_reg2->chrname_rightTra2; pos4 = clip_reg2->rightClipPosTra2; }

		// compute distances
		if(chrname1.compare(chrname3)==0){
			dist1 = pos1 - pos3;
			if(dist1<0) dist1 = -dist1;
		}
		if(chrname1.compare(chrname4)==0){
			dist2 = pos1 - pos4;
			if(dist2<0) dist2 = -dist2;
		}

		if(chrname2.compare(chrname3)==0){
			dist3 = pos2 - pos3;
			if(dist3<0) dist3 = -dist3;
		}
		if(chrname2.compare(chrname3)==0){
			dist4 = pos2 - pos4;
			if(dist4<0) dist4 = -dist4;
		}
	}

	dist_vec.push_back(dist1);
	dist_vec.push_back(dist2);
	dist_vec.push_back(dist3);
	dist_vec.push_back(dist4);
	return dist_vec;
}

// fill variant sequences
void Genome::genomeFillVarseq(){
//	Chrome *chr;
//	for(size_t i=0; i<chromeVector.size(); i++){
//		chr = chromeVector.at(i);
//		if(chr->process_block_num)
//			chr->chrFillVarseq();
//	}

	// fill sequences for TRA
	genomeFillVarseqTra();
}

void Genome::genomeFillVarseqTra(){
	size_t i, j;
	Chrome *chr;
	mateClipReg_t *clip_reg;
	ofstream cns_info_file;
	string cns_info_filename_tra;

	cns_info_filename_tra = out_dir_tra + "/" + "tra_cns_info";
	cns_info_file.open(cns_info_filename_tra);
	if(!cns_info_file.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << cns_info_filename_tra << endl;
		exit(1);
	}

	for(i=0; i<chromeVector.size(); i++){
		chr = chromeVector.at(i);
		for(j=0; j<chr->mateClipRegVector.size(); j++){
			clip_reg = chr->mateClipRegVector.at(j);
			if(clip_reg->valid_flag and (clip_reg->call_success_flag or clip_reg->tra_rescue_success_flag) and clip_reg->sv_type==VAR_TRA)
				fillVarseqSingleMateClipReg(clip_reg, cns_info_file);
		}
	}
	cns_info_file.close();
}

void Genome::fillVarseqSingleMateClipReg(mateClipReg_t *clip_reg, ofstream &cns_info_file){
	varCand *var_cand_tmp;
	string tmpdir;
	reg_t *reg;
	vector<int32_t> ref_shift_size_vec;
	vector<size_t> query_loc_vec; // [0]-blat_aln_id, [1]-query_id, [2]-start_query_pos, [3]-end_query_pos, [4]-start_local_ref_pos, [5]-end_loal_ref_pos, [6]-start_ref_pos, [7]-end_ref_pos
	blat_aln_t *blat_aln;
	size_t i, cns_extend_size;

	if(clip_reg->valid_flag){
		if(clip_reg->sv_type==VAR_TRA){
			// left part
			//if(clip_reg->leftClipPosTra1!=-1 and clip_reg->rightClipPosTra1!=-1 and clip_reg->chrname_leftTra1.compare(clip_reg->chrname_rightTra1)==0){
			if(clip_reg->leftClipPosTra1!=-1 and clip_reg->leftClipPosTra2!=-1 and clip_reg->chrname_leftTra1.compare(clip_reg->chrname_leftTra2)==0){

				// construct var_cand
				var_cand_tmp = new varCand();

				var_cand_tmp->chrname = clip_reg->chrname_leftTra1;
				var_cand_tmp->var_cand_filename = "";
				var_cand_tmp->out_dir_call = out_dir_tra;
				var_cand_tmp->misAln_filename = "";
				var_cand_tmp->inBamFile = paras->inBamFile;
				var_cand_tmp->fai = fai;
				var_cand_tmp->technology = paras->technology;
				var_cand_tmp->call_success = false;

				var_cand_tmp->blat_aligned_info_vec = NULL;
				var_cand_tmp->blat_var_cand_file = NULL;

				var_cand_tmp->minimap2_aligned_info_vec = NULL;
				var_cand_tmp->minimap2_var_cand_file = NULL;

				var_cand_tmp->min_sv_size = paras->min_sv_size_usr;
				var_cand_tmp->minReadsNumSupportSV = paras->minReadsNumSupportSV;
				var_cand_tmp->minClipEndSize = paras->minClipEndSize;
				var_cand_tmp->minConReadLen = paras->minConReadLen;
				var_cand_tmp->min_identity_match = paras->min_identity_match;
				var_cand_tmp->min_identity_merge = paras->min_identity_merge;
				var_cand_tmp->min_distance_merge = paras->min_distance_merge;
				var_cand_tmp->max_seg_num_per_read = paras->max_seg_num_per_read;
				var_cand_tmp->minMapQ = paras->minMapQ;
				var_cand_tmp->minHighMapQ = paras->minHighMapQ;

				var_cand_tmp->leftClipRefPos = clip_reg->leftClipPosTra1;
				//var_cand_tmp->rightClipRefPos = clip_reg->rightClipPosTra1;
				var_cand_tmp->rightClipRefPos = clip_reg->leftClipPosTra2;

				reg = new reg_t();
				reg->chrname = clip_reg->chrname_leftTra1;
				reg->startRefPos = clip_reg->leftClipPosTra1;
				//reg->endRefPos = clip_reg->rightClipPosTra1;
				reg->endRefPos = clip_reg->leftClipPosTra2;
				reg->var_type = clip_reg->sv_type;
				reg->call_success_status = false;
				reg->short_sv_flag = false;
				reg->zero_cov_flag = false;
				reg->aln_seg_end_flag = false;
				reg->query_pos_invalid_flag = false;
				reg->large_indel_flag = false;
				reg->gt_type = -1;
				reg->gt_seq = "";
				reg->AF = 0;
				reg->supp_num = reg->DP = 0;
				reg->discover_level = VAR_DISCOV_L_UNUSED;

				var_cand_tmp->varVec.push_back(reg);

				// construct the variant region
				var_cand_tmp->readsfilename = out_dir_tra  + "/" + "tra_reads_" + reg->chrname + "_" + to_string(reg->startRefPos) + "-" + to_string(reg->endRefPos) + ".fa";
				var_cand_tmp->ctgfilename = out_dir_tra  + "/" + "tra_contig_" + reg->chrname + "_" + to_string(reg->startRefPos) + "-" + to_string(reg->endRefPos) + ".fa";
				var_cand_tmp->refseqfilename = out_dir_tra  + "/" + "tra_refseq_" + reg->chrname + "_" + to_string(reg->startRefPos) + "-" + to_string(reg->endRefPos) + ".fa";
				var_cand_tmp->alnfilename = out_dir_tra  + "/" + "tra_blat_" + reg->chrname + "_" + to_string(reg->startRefPos) + "-" + to_string(reg->endRefPos) + ".sim4";
				tmpdir = out_dir_tra + "/" + "tmp_tra_" + reg->chrname + "_" + to_string(reg->startRefPos) + "-" + to_string(reg->endRefPos);

				//set genotyping parameters
				var_cand_tmp->setGtParas(paras->gt_min_sig_size, paras->gt_size_ratio_match, paras->gt_min_identity_merge, paras->gt_homo_ratio, paras->gt_hete_ratio, paras->minReadsNumSupportSV);

				// process monitor killed blat work
				var_cand_tmp->max_proc_running_minutes = paras->max_proc_running_minutes_call;
				var_cand_tmp->killed_blat_work_vec = &paras->killed_blat_work_vec;
				var_cand_tmp->killed_blat_work_file = &paras->killed_blat_work_file;
				var_cand_tmp->mtx_killed_blat_work = &paras->mtx_killed_blat_work;

				// process monitor killed blat work
				//var_cand_tmp->max_proc_running_minutes = paras->max_proc_running_minutes_call;
				var_cand_tmp->killed_minimap2_work_vec = &paras->killed_minimap2_work_vec;
				var_cand_tmp->killed_minimap2_work_file = &paras->killed_minimap2_work_file;
				var_cand_tmp->mtx_killed_minimap2_work = &paras->mtx_killed_minimap2_work;

				for(i=0; i<3; i++){
					cns_extend_size = paras->cnsSideExtSizeClip * i;
					// local consensus
					performLocalCnsTra(var_cand_tmp->readsfilename, var_cand_tmp->ctgfilename, var_cand_tmp->refseqfilename, var_cand_tmp->clusterfilename, tmpdir, paras->technology, paras->min_identity_match, reg->endRefPos-reg->startRefPos+1, paras->num_threads_per_cns_work, var_cand_tmp->varVec, reg->chrname, paras->inBamFile, fai, cns_extend_size, cns_info_file);

					ref_shift_size_vec = getRefShiftSize(var_cand_tmp->refseqfilename);
					var_cand_tmp->ref_left_shift_size = ref_shift_size_vec.at(0);
					var_cand_tmp->ref_right_shift_size = ref_shift_size_vec.at(1);

					var_cand_tmp->sv_type = clip_reg->sv_type;
					var_cand_tmp->dup_num = clip_reg->dup_num;
					var_cand_tmp->margin_adjusted_flag = false;
					var_cand_tmp->ctg_num = getCtgCount(var_cand_tmp->ctgfilename);

					// consensus status
					if (isFileExist(var_cand_tmp->ctgfilename)) var_cand_tmp->cns_success = true;
					else var_cand_tmp->cns_success = false;
					if(var_cand_tmp->cns_success){

						var_cand_tmp->loadBlatAlnData(); // BLAT align
						if(var_cand_tmp->align_success){
							// compute the aligned query locations for TRA
							query_loc_vec = computeQueryLocTra(var_cand_tmp, clip_reg, LEFT_END);
							if(query_loc_vec.size()>0){
								clip_reg->leftClipPosTra1 = query_loc_vec.at(6);
								//clip_reg->rightClipPosTra1 = query_loc_vec.at(7);
								clip_reg->leftClipPosTra2 = query_loc_vec.at(7);

								blat_aln = var_cand_tmp->blat_aln_vec.at(query_loc_vec.at(0));
								// get sequences
								FastaSeqLoader refseqloader(var_cand_tmp->refseqfilename);
								clip_reg->refseq_tra = refseqloader.getFastaSeqByPos(0, query_loc_vec.at(4), query_loc_vec.at(5), ALN_PLUS_ORIENT);
								FastaSeqLoader ctgseqloader(var_cand_tmp->ctgfilename);
								clip_reg->altseq_tra = ctgseqloader.getFastaSeqByPos(query_loc_vec.at(1), query_loc_vec.at(2), query_loc_vec.at(3), blat_aln->aln_orient);

								var_cand_tmp->call_success = true;
								break;
							}
						}
					}
				}

				// release memory
				var_cand_tmp->destroyVarCand();
				delete var_cand_tmp;
			}

			// right part
			//if(clip_reg->leftClipPosTra2!=-1 and clip_reg->rightClipPosTra2!=-1 and clip_reg->chrname_leftTra2.compare(clip_reg->chrname_rightTra2)==0){
			if(clip_reg->rightClipPosTra1!=-1 and clip_reg->rightClipPosTra2!=-1 and clip_reg->chrname_rightTra1.compare(clip_reg->chrname_rightTra2)==0){
				// construct var_cand
				var_cand_tmp = new varCand();

				//var_cand_tmp->chrname = clip_reg->chrname_leftTra2;
				var_cand_tmp->chrname = clip_reg->chrname_rightTra1;
				var_cand_tmp->var_cand_filename = "";
				var_cand_tmp->out_dir_call = out_dir_tra;
				var_cand_tmp->misAln_filename = "";
				var_cand_tmp->inBamFile = paras->inBamFile;
				var_cand_tmp->fai = fai;
				var_cand_tmp->technology = paras->technology;
				var_cand_tmp->align_success = false;

				var_cand_tmp->blat_aligned_info_vec = NULL;
				var_cand_tmp->blat_var_cand_file = NULL;

				var_cand_tmp->minimap2_aligned_info_vec = NULL;
				var_cand_tmp->minimap2_var_cand_file = NULL;

				//var_cand_tmp->leftClipRefPos = clip_reg->leftClipPosTra2;
				var_cand_tmp->leftClipRefPos = clip_reg->rightClipPosTra1;
				var_cand_tmp->rightClipRefPos = clip_reg->rightClipPosTra2;

				var_cand_tmp->min_sv_size = paras->min_sv_size_usr;
				var_cand_tmp->minReadsNumSupportSV = paras->minReadsNumSupportSV;
				var_cand_tmp->minClipEndSize = paras->minClipEndSize;
				var_cand_tmp->minConReadLen = paras->minConReadLen;
				var_cand_tmp->min_identity_match = paras->min_identity_match;
				var_cand_tmp->min_identity_merge = paras->min_identity_merge;
				var_cand_tmp->min_distance_merge = paras->min_distance_merge;
				var_cand_tmp->max_seg_num_per_read = paras->max_seg_num_per_read;
				var_cand_tmp->minMapQ = paras->minMapQ;
				var_cand_tmp->minHighMapQ = paras->minHighMapQ;

				reg = new reg_t();
				//reg->chrname = clip_reg->chrname_leftTra2;
				reg->chrname = clip_reg->chrname_rightTra1;
				//reg->startRefPos = clip_reg->leftClipPosTra2;
				reg->startRefPos = clip_reg->rightClipPosTra1;
				reg->endRefPos = clip_reg->rightClipPosTra2;
				reg->var_type = clip_reg->sv_type;
				reg->call_success_status = false;
				reg->short_sv_flag = false;
				reg->zero_cov_flag = false;
				reg->aln_seg_end_flag = false;
				reg->query_pos_invalid_flag = false;
				reg->large_indel_flag = false;
				reg->gt_type = -1;
				reg->gt_seq = "";
				reg->AF = 0;
				reg->supp_num = reg->DP = 0;
				reg->discover_level = VAR_DISCOV_L_UNUSED;

				var_cand_tmp->varVec.push_back(reg);

				// construct the variant region
				var_cand_tmp->readsfilename = out_dir_tra  + "/" + "tra_reads_" + reg->chrname + "_" + to_string(reg->startRefPos) + "-" + to_string(reg->endRefPos) + ".fa";
				var_cand_tmp->ctgfilename = out_dir_tra  + "/" + "tra_contig_" + reg->chrname + "_" + to_string(reg->startRefPos) + "-" + to_string(reg->endRefPos) + ".fa";
				var_cand_tmp->refseqfilename = out_dir_tra  + "/" + "tra_refseq_" + reg->chrname + "_" + to_string(reg->startRefPos) + "-" + to_string(reg->endRefPos) + ".fa";
				var_cand_tmp->alnfilename = out_dir_tra  + "/" + "tra_blat_" + reg->chrname + "_" + to_string(reg->startRefPos) + "-" + to_string(reg->endRefPos) + ".sim4";
				tmpdir = out_dir_tra + "/" + "tmp_tra_" + reg->chrname + "_" + to_string(reg->startRefPos) + "-" + to_string(reg->endRefPos);

				// process monitor killed blat work
				var_cand_tmp->max_proc_running_minutes = paras->max_proc_running_minutes_call;
				var_cand_tmp->killed_blat_work_vec = &paras->killed_blat_work_vec;
				var_cand_tmp->killed_blat_work_file = &paras->killed_blat_work_file;
				var_cand_tmp->mtx_killed_blat_work = &paras->mtx_killed_blat_work;

				// process monitor killed minimap2 work
				//var_cand_tmp->max_proc_running_minutes = paras->max_proc_running_minutes_call;
				var_cand_tmp->killed_minimap2_work_vec = &paras->killed_minimap2_work_vec;
				var_cand_tmp->killed_minimap2_work_file = &paras->killed_minimap2_work_file;
				var_cand_tmp->mtx_killed_minimap2_work = &paras->mtx_killed_minimap2_work;

				//set genotyping parameters
				var_cand_tmp->setGtParas(paras->gt_min_sig_size, paras->gt_size_ratio_match, paras->gt_min_identity_merge, paras->gt_homo_ratio, paras->gt_hete_ratio, paras->minReadsNumSupportSV);

				for(i=0; i<3; i++){
					cns_extend_size = paras->cnsSideExtSizeClip * i;
					// local consensus
					performLocalCnsTra(var_cand_tmp->readsfilename, var_cand_tmp->ctgfilename, var_cand_tmp->refseqfilename, var_cand_tmp->clusterfilename, tmpdir, paras->technology, paras->min_identity_match, reg->endRefPos-reg->startRefPos+1, paras->num_threads_per_cns_work, var_cand_tmp->varVec, reg->chrname, paras->inBamFile, fai, cns_extend_size, cns_info_file);

					ref_shift_size_vec = getRefShiftSize(var_cand_tmp->refseqfilename);
					var_cand_tmp->ref_left_shift_size = ref_shift_size_vec.at(0);
					var_cand_tmp->ref_right_shift_size = ref_shift_size_vec.at(1);

					var_cand_tmp->sv_type = clip_reg->sv_type;
					var_cand_tmp->dup_num = clip_reg->dup_num;
					var_cand_tmp->margin_adjusted_flag = false;
					var_cand_tmp->ctg_num = getCtgCount(var_cand_tmp->ctgfilename);

					// consensus status
					if (isFileExist(var_cand_tmp->ctgfilename)) var_cand_tmp->cns_success = true;
					else var_cand_tmp->cns_success = false;
					if(var_cand_tmp->cns_success){

						var_cand_tmp->loadBlatAlnData(); // BLAT align
						if(var_cand_tmp->align_success){
							// compute the aligned query locations for TRA
							query_loc_vec = computeQueryLocTra(var_cand_tmp, clip_reg, RIGHT_END);
							if(query_loc_vec.size()>0){
								//clip_reg->leftClipPosTra2 = query_loc_vec.at(6);
								clip_reg->rightClipPosTra1 = query_loc_vec.at(6);
								clip_reg->rightClipPosTra2 = query_loc_vec.at(7);

								blat_aln = var_cand_tmp->blat_aln_vec.at(query_loc_vec.at(0));
								// get sequences
								FastaSeqLoader refseqloader(var_cand_tmp->refseqfilename);
								clip_reg->refseq_tra2 = refseqloader.getFastaSeqByPos(0, query_loc_vec.at(4), query_loc_vec.at(5), ALN_PLUS_ORIENT);
								FastaSeqLoader ctgseqloader(var_cand_tmp->ctgfilename);
								clip_reg->altseq_tra2 = ctgseqloader.getFastaSeqByPos(query_loc_vec.at(1), query_loc_vec.at(2), query_loc_vec.at(3), blat_aln->aln_orient);

								var_cand_tmp->call_success = true;
								break;
							}
						}
					}
				}

				// release memory
				var_cand_tmp->destroyVarCand();
				delete var_cand_tmp;
			}
		}
	}
}

// perform local consensus
void Genome::performLocalCnsTra(string &readsfilename, string &contigfilename, string &refseqfilename, string &clusterfilename, string &tmpdir, string &technology, double min_identity_match, int32_t sv_len_est, size_t num_threads_per_cns_work, vector<reg_t*> &varVec, string &chrname, string &inBamFile, faidx_t *fai, size_t cns_extend_size, ofstream &cns_info_file){

	localCns local_cns(readsfilename, contigfilename, refseqfilename, clusterfilename, tmpdir, technology, min_identity_match, sv_len_est, num_threads_per_cns_work, varVec, chrname, inBamFile, fai, cns_extend_size, paras->expected_cov_cns, paras->min_input_cov_canu, paras->max_ultra_high_cov, paras->minMapQ, paras->minHighMapQ, paras->delete_reads_flag, paras->keep_failed_reads_flag, true, paras->minClipEndSize, paras->minConReadLen, paras->min_sv_size_usr, paras->minReadsNumSupportSV, paras->max_seg_size_ratio_usr, paras->max_seg_nm_ratio_usr);

	// extract the corresponding refseq from reference
	local_cns.extractRefseq();

	// extract the reads data from BAM file
	local_cns.extractReadsDataFromBAM();

	// local consensus
//	local_cns.localCnsCanu();
	local_cns.cnsByPoa();

	// record consensus information
	local_cns.recordCnsInfo(cns_info_file);

	// empty the varVec
	//varVec.clear();
}

vector<int32_t> Genome::getRefShiftSize(string &refseqfilename){
	ifstream infile;
	string header_str;
	int32_t left_shift_size, right_shift_size;
	vector<string> str_vec;
	vector<int32_t> shift_size_vec;

	// ref shift size
	infile.open(refseqfilename);
	if(!infile.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file " << refseqfilename << endl;
		exit(1);
	}

	getline(infile, header_str);
	str_vec = split(header_str, "___");
	left_shift_size = stoi(str_vec[1]);
	right_shift_size = stoi(str_vec[2]);

	shift_size_vec.push_back(left_shift_size);
	shift_size_vec.push_back(right_shift_size);

	infile.close();

	return shift_size_vec;
}

// compute the aligned query locations for TRA
vector<size_t> Genome::computeQueryLocTra(varCand *var_cand, mateClipReg_t *clip_reg, size_t end_flag){
	vector<size_t> query_loc_vec;
	size_t i, j, aln_orient, missing_size_head, missing_size_inner, missing_size_tail;
	int32_t maxValue, maxIdx, query_dist, query_part_flag_maxValue;  // query_part_flag: 0-head, 1-inner, 2-tail
	blat_aln_t *blat_aln;
	aln_seg_t *aln_seg1, *aln_seg2;
	int32_t start_query_pos, end_query_pos, start_ref_pos, end_ref_pos, start_local_ref_pos, end_local_ref_pos, start_ref_pos_tra, end_ref_pos_tra;
	int32_t k, segIdx_start, segIdx_end;
	bool valid_tra_flag = false;

	//if(end_flag==LEFT_END) { start_ref_pos_tra = clip_reg->leftClipPosTra1; end_ref_pos_tra = clip_reg->rightClipPosTra1; }
	//else { start_ref_pos_tra = clip_reg->leftClipPosTra2; end_ref_pos_tra = clip_reg->rightClipPosTra2; }
	if(end_flag==LEFT_END) { start_ref_pos_tra = clip_reg->leftClipPosTra1; end_ref_pos_tra = clip_reg->leftClipPosTra2; }
	else { start_ref_pos_tra = clip_reg->rightClipPosTra1; end_ref_pos_tra = clip_reg->rightClipPosTra2; }

	for(i=0; i<var_cand->blat_aln_vec.size(); i++){
		blat_aln = var_cand->blat_aln_vec.at(i);
		if(blat_aln->valid_aln){
			aln_orient = blat_aln->aln_orient;

			// query head
			if(aln_orient==ALN_PLUS_ORIENT){
				if(blat_aln->aln_segs.at(0)->subject_start>MIN_CLIP_DIST_THRES)
					missing_size_head = blat_aln->aln_segs.at(0)->query_start - 1;
				else
					missing_size_head = 0;
			}else{
				if(blat_aln->subject_len-blat_aln->aln_segs.at(blat_aln->aln_segs.size()-1)->subject_end>MIN_CLIP_DIST_THRES)
					missing_size_head = blat_aln->query_len - blat_aln->aln_segs.at(blat_aln->aln_segs.size()-1)->query_end;
				else
					missing_size_head = 0;
			}

			// query inner part
			maxValue = 0; maxIdx = -1;
			for(j=0; j<blat_aln->aln_segs.size()-1; j++){
				aln_seg1 = blat_aln->aln_segs.at(j);
				aln_seg2 = blat_aln->aln_segs.at(j+1);
				query_dist = aln_seg2->query_start - aln_seg1->query_end;
				if(query_dist>maxValue) {
					maxValue = query_dist;
					maxIdx = j;
				}
			}
			missing_size_inner = maxValue;

			// query tail
			if(aln_orient==ALN_PLUS_ORIENT){
				if(blat_aln->subject_len-blat_aln->aln_segs.at(blat_aln->aln_segs.size()-1)->subject_end>MIN_CLIP_DIST_THRES)
					missing_size_tail = blat_aln->query_len - blat_aln->aln_segs.at(blat_aln->aln_segs.size()-1)->query_end;
				else
					missing_size_tail = 0;
			}else{
				if(blat_aln->aln_segs.at(0)->subject_start>MIN_CLIP_DIST_THRES)
					missing_size_tail = blat_aln->aln_segs.at(0)->query_start - 1;
				else
					missing_size_tail = 0;
			}

			// choose the maximum and the second maximum
			query_part_flag_maxValue = 0; maxValue = missing_size_head;
			if(maxValue<(int32_t)missing_size_inner) { query_part_flag_maxValue = 1; maxValue = missing_size_inner; }
			if(maxValue<(int32_t)missing_size_tail) { query_part_flag_maxValue = 2; maxValue = missing_size_tail; }

			aln_seg1 = aln_seg2 = NULL;
			if(maxValue>=MIN_CLIP_DIST_THRES and query_part_flag_maxValue==1){ // inner missing
				aln_seg1 = blat_aln->aln_segs.at(maxIdx);
				aln_seg2 = blat_aln->aln_segs.at(maxIdx+1);

				// compute locations
				start_ref_pos = aln_seg1->ref_end + 1;
				end_ref_pos = aln_seg2->ref_start - 1;
				start_local_ref_pos = aln_seg1->subject_end + 1;
				end_local_ref_pos = aln_seg2->subject_start - 1;
				if(aln_orient==ALN_PLUS_ORIENT){ // plus orient
					start_query_pos = aln_seg1->query_end + 1;
					end_query_pos = aln_seg2->query_start - 1;
				}else{ // minus orient
					start_query_pos = blat_aln->query_len - (aln_seg2->query_start - 1) + 1;
					end_query_pos = blat_aln->query_len - (aln_seg1->query_end + 1) + 1;
				}
				//cout << ">>>>>> missing inner TRA: " << var_cand->chrname << ":" << start_ref_pos << "-" << end_ref_pos << ", " << blat_aln->query_id << ":" << start_query_pos << "-" << end_query_pos << endl;
			}else if(missing_size_inner<MIN_CLIP_DIST_THRES){
				// no missing parts for inner query, and compute the start and end align segments
				segIdx_start = segIdx_end = -1;
				for(j=0; j<blat_aln->aln_segs.size(); j++){
					aln_seg1 = blat_aln->aln_segs[j];
					if(aln_seg1->ref_end-aln_seg1->ref_start+1>=MIN_VALID_BLAT_SEG_SIZE and aln_seg1->ref_start+CLIP_END_EXTEND_SIZE>=start_ref_pos_tra and aln_seg1->ref_start<=end_ref_pos_tra+CLIP_END_EXTEND_SIZE){ // ignore short align segments
						segIdx_start = j;
						break;
					}
				}
				if(segIdx_start!=-1){
					for(k=blat_aln->aln_segs.size()-1; k>=segIdx_start; k--){
						aln_seg2 = blat_aln->aln_segs.at(k);
						if(aln_seg2->ref_end-aln_seg2->ref_start+1>=MIN_VALID_BLAT_SEG_SIZE and aln_seg2->ref_end+CLIP_END_EXTEND_SIZE>=start_ref_pos_tra and aln_seg2->ref_end<=end_ref_pos_tra+CLIP_END_EXTEND_SIZE){ // ignore short align segments
							segIdx_end = k;
							break;
						}
					}
				}
				if(segIdx_start!=-1 and segIdx_end!=-1) {
					aln_seg1 = blat_aln->aln_segs.at(segIdx_start);
					aln_seg2 = blat_aln->aln_segs.at(segIdx_end);

					// compute locations
					start_ref_pos = aln_seg1->ref_start;
					end_ref_pos = aln_seg2->ref_end;
					start_local_ref_pos = aln_seg1->subject_start;
					end_local_ref_pos = aln_seg2->subject_end;
					if(aln_orient==ALN_PLUS_ORIENT){ // plus orient
						start_query_pos = aln_seg1->query_start;
						end_query_pos = aln_seg2->query_end;
					}else{ // minus orient
						start_query_pos = blat_aln->query_len - aln_seg2->query_end + 1;
						end_query_pos = blat_aln->query_len - aln_seg1->query_start + 1;
					}
					//cout << ">>>>>> inner TRA: " << var_cand->chrname << ":" << start_ref_pos << "-" << end_ref_pos << ", " << blat_aln->query_id << ":" << start_query_pos << "-" << end_query_pos << endl;
				}else aln_seg1 = aln_seg2 = NULL;
			}

			if(aln_seg1 and aln_seg2){ // confirm TRA locations
				valid_tra_flag = isValidTraReg(start_ref_pos, end_ref_pos, start_ref_pos_tra, end_ref_pos_tra);
				if(valid_tra_flag){
					query_loc_vec.push_back(blat_aln->blat_aln_id);
					query_loc_vec.push_back(blat_aln->query_id);
					query_loc_vec.push_back(start_query_pos);
					query_loc_vec.push_back(end_query_pos);
					query_loc_vec.push_back(start_local_ref_pos);
					query_loc_vec.push_back(end_local_ref_pos);
					query_loc_vec.push_back(start_ref_pos);
					query_loc_vec.push_back(end_ref_pos);
				}
			}
		}
		if(valid_tra_flag) break;
	}

	return query_loc_vec;
}

// save variants to file
void Genome::genomeSaveCallSV2File(){
	mkdir(out_dir_result.c_str(), S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);  // create the directory for final results

	Chrome *chr;
	for(size_t i=0; i<chromeVector.size(); i++){
		chr = chromeVector.at(i);
		//if(chr->process_block_num)
			//chr->saveCallSV2File();
			chr->saveCallSV2File02();
	}

	// save TRA
	saveTraCall2File();
}

// save TRA
void Genome::saveTraCall2File(){
	size_t i, j, k;
	Chrome *chr;
	ofstream outfile_tra;
	vector<mateClipReg_t*> mate_clipReg_vec;
	mateClipReg_t *clip_reg;
	string line, chrname_tmp, header_line_bedpe, mate_reg_str;
	vector<string> chrname_vec;
	int32_t sv_len, sv_len2;

	outfile_tra.open(out_filename_result_tra);
	if(!outfile_tra.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << out_filename_result_tra << endl;
		exit(1);
	}

	header_line_bedpe = getCallFileHeaderBedpe(paras->sample);
	outfile_tra << header_line_bedpe << endl;

	for(i=0; i<chromeVector.size(); i++){
		chr = chromeVector.at(i);
		mate_clipReg_vec = chr->mateClipRegVector;
		for(j=0; j<mate_clipReg_vec.size(); j++){
			clip_reg = mate_clipReg_vec.at(j);
			if(clip_reg->valid_flag and (clip_reg->call_success_flag or clip_reg->tra_rescue_success_flag) and clip_reg->sv_type==VAR_TRA and (clip_reg->leftClipPosTra1>0 or clip_reg->leftClipPosTra2>0) and (clip_reg->rightClipPosTra1>0 or clip_reg->rightClipPosTra2>0)){
//				if(clip_reg->left_var_cand_tra==NULL or clip_reg->right_var_cand_tra==NULL){
//					cout << "i=" << i << ", j=" << j << endl;
//				}
				//line = clip_reg->left_var_cand_tra->chrname;
				//line = clip_reg->chrname_leftTra1.size()>0 ? clip_reg->chrname_leftTra1 : "-";
				chrname_vec = getLeftRightPartChrname(clip_reg);
				line = chrname_vec.at(0).size()>0 ? chrname_vec.at(0) : "-";
				if(clip_reg->leftClipPosTra1>0) line += "\t" + to_string(clip_reg->leftClipPosTra1);
				else line += "\t-";
				if(clip_reg->leftClipPosTra2>0) line += "\t" + to_string(clip_reg->leftClipPosTra2);
				else line += "\t-";
				chrname_tmp = chrname_vec.at(1).size()>0 ? chrname_vec.at(1) : "-";
				line += "\t" + chrname_tmp;
				if(clip_reg->rightClipPosTra1>0) line += "\t" + to_string(clip_reg->rightClipPosTra1);
				else line += "\t-";
				if(clip_reg->rightClipPosTra2>0) line += "\t" + to_string(clip_reg->rightClipPosTra2);
				else line += "\t-";
				line += "\tTRA";

				sv_len = sv_len2 = -1;
				if(clip_reg->leftClipPosTra1>0 and clip_reg->leftClipPosTra2>0){
					sv_len = clip_reg->leftClipPosTra2 - clip_reg->leftClipPosTra1 + 1;
					line += "\t" + to_string(sv_len);
				}else line += "\t-";
				if(clip_reg->rightClipPosTra1>0 and clip_reg->rightClipPosTra2>0){
					sv_len2 = clip_reg->rightClipPosTra2 - clip_reg->rightClipPosTra1 + 1;
					line += "\t" + to_string(sv_len2);
				}else line += "\t-";

				// BND mate regions
				mate_reg_str = clip_reg->bnd_mate_reg_strs[0];
				for(k=1; k<4; k++) mate_reg_str += ";" + clip_reg->bnd_mate_reg_strs[k];
				line += "\t" + mate_reg_str;

				if(clip_reg->refseq_tra.size()>0) line += "\t" + clip_reg->refseq_tra;
				else line += "\t-";
				if(clip_reg->altseq_tra.size()>0) line += "\t" + clip_reg->altseq_tra;
				else line += "\t-";
				if(clip_reg->refseq_tra2.size()>0) line += "\t" + clip_reg->refseq_tra2;
				else line += "\t-";
				if(clip_reg->altseq_tra2.size()>0) line += "\t" + clip_reg->altseq_tra2;
				else line += "\t-";

//				if(reg->gt_seq.size()==0) reg->gt_seq = GT_STR_DEFAULT;
//				line += "\t" + reg->gt_seq;

				// size select
				//if((sv_len==-1 and sv_len2==-1) or (sv_len!=-1 and sv_len>=paras->min_sv_size_usr and sv_len<=paras->max_sv_size_usr) or (sv_len2!=-1 and sv_len2>=paras->min_sv_size_usr and sv_len2<=paras->max_sv_size_usr)) // deleted on 2025-03-20
				if((sv_len==-1 and sv_len2==-1) or (sv_len!=-1 and sv_len>=paras->min_sv_size_usr_final and sv_len<=paras->max_sv_size_usr) or (sv_len2!=-1 and sv_len2>=paras->min_sv_size_usr_final and sv_len2<=paras->max_sv_size_usr))
					outfile_tra << line << endl;
					//cout << line << endl;
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
		copySingleFile(chr->out_filename_call_indel, out_file_indel); // indel
		copySingleFile(chr->out_filename_call_clipReg, out_file_clipReg); // clip_reg
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

// save results in VCF file format for detect command
//void Genome::saveResultVCFDetect(){
//	string outIndelFile, outSnvFile;
//
//	cout << "Output results to file ..." << endl;
//
//	// initialize the output filenames
//	if(paras->outFilePrefix.size()>0) {
//		outIndelFile = out_dir_detect + "/" + paras->outFilePrefix + "_INDEL.vcf";
//		outSnvFile = out_dir_detect + "/" + paras->outFilePrefix + "_SNV.vcf";
//	}else {
//		outIndelFile = out_dir_detect + "/" + "INDEL.vcf";
//		outSnvFile = out_dir_detect + "/" + "SNV.vcf";
//	}
//
//	// save results
//	saveIndelVCFDetect(out_filename_detect_indel, outIndelFile);
//	saveSnvVCFDetect(out_filename_detect_snv, outSnvFile);
//}

// save indel in VCF file format for detect command
void Genome::saveResultToVCF(string &in, string &out_vcf){
	string line, line_vcf, chr, start_pos, id, ref, alt, qual, filter, info, format, format_val, sample, reg;
	string end_pos, sv_type, sv_len, dup_num, blat_inner, extra_info;
	ifstream infile;
	ofstream outfile;
	vector<string> str_vec;

	string bnd_id_str, mate_bnd_id_str, mate_reg_str, mate_bnd_str, mate_bnd_str2, reg_str, mate_chr, chrname1, chrname2;
	vector<string> bnd_str_vec, bnd_str_vec2;
	int64_t tra_pos_arr[4];
	int32_t i, j, checked_arr[4][2]; //, ins_id, del_id, dup_id, inv_id, tra_id;
	vector<BND_t*> bnd_vec, sub_bnd_vec;
	string format_name_bnd;

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

			//line_vcf = str_vec.at(0) + "\t" + str_vec.at(1) + "\t."; // CHROM, POS, ID
			//reg = str_vec[0] + ":" + str_vec[1] + "-" + str_vec[2];
			//seq = fai_fetch(fai, reg.c_str(), &seq_len);

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
				sv_type = str_vec.at(4);
				sv_len = str_vec.at(5);
				if(sv_type.compare(VAR_UNC_STR)==0 or sv_type.compare(VAR_UNC_STR1)==0 or sv_type.compare(VAR_UNC_STR2)==0) id = ref = alt = "."; // UNC
				info = "SVTYPE=" + sv_type + ";" + "SVLEN=" + sv_len + ";END=" + end_pos;	// INFO: SVTYPE=sv_type;SVLEN=sv_len;END=end
				if(sv_type.compare(VAR_DUP_STR)==0 or sv_type.compare(VAR_DUP_STR1)==0) { // DUPNUM=dup_num
					dup_num = str_vec.at(6);
					info += ";DUPNUM=" + dup_num;
				}
				info += ";" + extra_info; // DP=xx;AF=yy
				if(str_vec.size()==MAX_BED_COLS_NUM){ // convert 'shortSV' to 'BLATINNER' into info
					blat_inner = "BLATINNER";
					info += ";" + blat_inner;
				}

//				format = "GT";		// FORMAT and the values
//				format_val = "./.";

				// id = "";
				// if(sv_type.compare(VAR_INS_STR)==0){
				// 	ins_id++;
				// 	id = sv_type + "." + to_string(ins_id);
				// }else if(sv_type.compare(VAR_DEL_STR)==0){
				// 	del_id++;
				// 	id = sv_type + "." + to_string(del_id);
				// }else if(sv_type.compare(VAR_DUP_STR)==0){
				// 	dup_id++;
				// 	id = sv_type + "." + to_string(dup_id);
				// }else if(sv_type.compare(VAR_INV_STR)==0){
				// 	inv_id++;
				// 	id = sv_type + "." + to_string(inv_id);
				// }else if(sv_type.compare(VAR_TRA_STR)==0){
				// 	tra_id++;
				// 	id = sv_type + "." + to_string(tra_id);
				// }

//				line_vcf = chr + "\t" + start_pos + "\t" + id + "\t" + ref + "\t" + alt + "\t" + qual + "\t" + filter + "\t" + info + "\t" + format + "\t" + format_val;

				format = str_vec.at(10);			// FORMAT
				format_val = str_vec.at(11);		// gt_value

				line_vcf = chr + "\t" + start_pos + "\t" + id + "\t" + ref + "\t" + alt + "\t" + qual + "\t" + filter + "\t" + info + "\t" + format + "\t" + format_val;
				outfile << line_vcf << endl;

			}else{	// TRA, MIX
				if(str_vec.at(7).compare(VAR_TRA_STR)==0){ // TRA, then convert to BND

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
							sub_bnd_vec = generateBNDItems(i, RIGHT_END, checked_arr, chrname1, chrname2, tra_pos_arr, bnd_str_vec, fai);
							for(j=0; j<(int32_t)sub_bnd_vec.size(); j++) {
								outfile << sub_bnd_vec.at(j)->vcf_line << endl;
								//cout << sub_bnd_vec.at(j)->vcf_line << endl;
								bnd_vec.push_back(sub_bnd_vec.at(j));
							}

							// left end
							sub_bnd_vec = generateBNDItems(i, LEFT_END, checked_arr, chrname1, chrname2, tra_pos_arr, bnd_str_vec, fai);
							for(j=0; j<(int32_t)sub_bnd_vec.size(); j++) {
								outfile << sub_bnd_vec.at(j)->vcf_line << endl;
								//cout << sub_bnd_vec.at(j)->vcf_line << endl;
								bnd_vec.push_back(sub_bnd_vec.at(j));
							}
						}
					}
				}else{ // MIX

					cout << __func__ << ": unexpected variant item: " << line << ", skipped." << endl;

//					chr = str_vec.at(0);		// CHROM
//					if(str_vec.at(1).compare("-")==0) start_pos = str_vec.at(2); // POS
//					else start_pos = str_vec.at(1);
//					id = ".";					// ID
//					ref = ".";					// REF
//					alt = ".";					// ALT
//
//					qual = ".";					// QUAL
//					filter = "PASS";			// FILTER
//
//					//chr2 = str_vec.at(3);
//					if(str_vec.at(4).compare("-")==0) end_pos = str_vec.at(5);
//					else end_pos = str_vec.at(4);
//					sv_type = str_vec.at(6);
//					sv_len = "False";
//					info = "CHR2=" + chr2 + ";END=" + end_pos + ";SVTYPE=" + sv_type + ";" + "SVLEN=" + sv_len;	// INFO: CHR2=chr2;END=end;SVTYPE=sv_type;SVLEN=sv_len
//
//					format = "GT";		// FORMAT and the values
//					format_val = "./.";
//
//					line_vcf = chr + "\t" + start_pos + "\t" + id + "\t" + ref + "\t" + alt + "\t" + qual + "\t" + filter + "\t" + info + "\t" + format + "\t" + format_val;
//					outfile << line_vcf << endl;
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

//// save SNV in VCF file format for detect command
//void Genome::saveSnvVCFToFile(string &in, string &out_vcf){
//	string line, data_vcf, ref, reg;
//	ifstream infile;
//	ofstream outfile;
//	vector<string> str_vec;
//	char *seq;
//	int seq_len;
//
//	// open files
//	infile.open(in);
//	if(!infile.is_open()){
//		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << in << endl;
//		exit(1);
//	}
//
//	outfile.open(out_vcf);
//	if(!outfile.is_open()){
//		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << out_vcf << endl;
//		exit(1);
//	}
//
//	// save VCF header
//	saveVCFHeader(outfile);
//
//	// save results
//	while(getline(infile, line))
//		if(line.size()>0){
//			str_vec = split(line, "\t");
//			data_vcf = str_vec[0] + "\t" + str_vec[1] + "\t."; // CHROM, POS, ID
//
//			reg = str_vec[0] + ":" + str_vec[1] + "-" + str_vec[1];
//			seq = fai_fetch(fai, reg.c_str(), &seq_len);
//			data_vcf += "\t" + (string)seq + "\t."; // REF, ALT
//
//			// QUAL, FILTER, INFO, FORMAT ........
//
//
//			outfile << data_vcf << endl;
//
//			free(seq);
//		}
//	infile.close();
//	outfile.close();
//}

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

	fp << "##phasing=None" << endl;
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
	//fp << "##INFO=<ID=BLATINNER,Number=1,Type=Flag,Description=\"variants called from inner parts of BLAT align segment\">" << endl;
//	fp << "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">" << endl;
//	fp << "##FILTER=<ID=q10,Description=\"Quality below 10\">" << endl;
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

	total = num_unc + num_ins + num_del + num_inv + num_dup + num_tra;

	cout << "############## Brief statistics for variants call ##############" << endl;
	cout << "There are " << total << " variants in total:" << endl;
	cout << "\t" << "insertions: " << num_ins << endl;
	cout << "\t" << "deletions: " << num_del << endl;
	cout << "\t" << "tandem duplications: " << num_dup << endl;
	cout << "\t" << "inversions: " << num_inv << endl;
	cout << "\t" << "translocations: " << num_tra << endl;
	cout << "\t" << "uncertain: " << num_unc << endl;

	infile.close();
}

size_t Genome::getSVTypeSingleLine(string &line){
	size_t sv_type = VAR_UNC;
	vector<string> str_vec;
	string str_tmp;

	str_vec = split(line, "\t");
	str_tmp = str_vec.at(4);
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
		if(str_vec.size()>MAX_BED_COLS_NUM and (str_vec.at(7).compare(VAR_TRA_STR)==0 or str_vec.at(7).compare(VAR_TRA_STR1)==0 or str_vec.at(7).compare(VAR_BND_STR)==0)){
			if(str_vec.at(7).compare(VAR_TRA_STR)==0 or str_vec.at(7).compare(VAR_TRA_STR1)==0) sv_type = VAR_TRA;
			else sv_type = VAR_BND;
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

	subsets = constructSubsetByChr(sv_vec);
	sortSVitem(subsets);
	//rmDupSVitem(subsets, paras->gt_size_ratio_match, paras->gt_min_identity_merge); // remove identical items, deleted on 2025-03-28
	rmDupSVitem(subsets, MIN_SIZE_RATIO_REDUNDANT_FILTER, MIN_IDENTITY_REDUNDANT_FILTER); // remove identical items
	outputResult(outfilename, subsets, filetype);

	remove(infilename.c_str());
	rename(outfilename.c_str(), infilename.c_str());

	destroyData(sv_vec);
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
