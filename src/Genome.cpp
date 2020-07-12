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

// Constructor with parameters
Genome::Genome(Paras *paras){
	this->paras = paras;
	init();
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
	string chrname_tmp, result_prefix;

	out_dir = paras->outDir;
	if(out_dir.size()>0){
		mkdir(out_dir.c_str(), S_IRWXU | S_IROTH);  // create the output directory

		out_dir_detect = out_dir + "/" + paras->out_dir_detect;
		out_dir_assemble = out_dir + "/" + paras->out_dir_assemble;
		out_dir_call = out_dir + "/" + paras->out_dir_call;
		out_dir_tra = out_dir + "/" + paras->out_dir_tra;
		out_dir_result = out_dir + "/" + paras->out_dir_result;
	}else{
		out_dir_detect = paras->out_dir_detect;
		out_dir_assemble = paras->out_dir_assemble;
		out_dir_call = paras->out_dir_call;
		out_dir_tra = paras->out_dir_tra;
		out_dir_result = paras->out_dir_result;
	}

	result_prefix = "";
	if(paras->outFilePrefix.size()) result_prefix = paras->outFilePrefix + "_";
	out_filename_detect_snv = out_dir_detect + "/" + result_prefix + "genome_SNV_candidates";
	out_filename_detect_indel = out_dir_detect + "/" + result_prefix + "genome_INDEL_candidate";
	out_filename_detect_clipReg = out_dir_detect + "/" + result_prefix + "genome_clipReg_candidate";
	out_filename_result_snv = out_dir_result + "/" + result_prefix + "genome_SNV";
	out_filename_result_indel = out_dir_result + "/" + result_prefix + "genome_INDEL.bed";
	out_filename_result_clipReg = out_dir_result + "/" + result_prefix + "genome_clipReg.bed";
	out_filename_result_tra = out_dir_result + "/" + result_prefix + "genome_TRA.bedpe";
	out_filename_result_vars = out_dir_result + "/" + result_prefix + "genome_variants.bed";
	blat_aln_info_filename_tra  = out_dir_tra + "/" + "blat_aln_info_tra";

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
		chr->setOutputDir(out_dir_detect, out_dir_assemble, out_dir_call);
		chromeVector.push_back(chr);
	}
	chromeVector.shrink_to_fit();
}

// allocate the Chrome node
Chrome* Genome::allocateChrome(string& chrname, int chrlen, faidx_t *fai){
	Chrome *chr_tmp;
	chr_tmp = new Chrome(chrname, chrlen, fai, paras);
	if(!chr_tmp){ cerr << "Genome: cannot allocate memory" << endl; exit(1); }
	return chr_tmp;
}

// get genome size
int Genome::getGenomeSize()
{
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
	for(i=0; i<chromeVector.size(); i++){
		chr = chromeVector.at(i);
		chr->chrFillDataEst(SIZE_EST_OP);
		if(paras->reg_sum_size_est>=paras->max_reg_sum_size_est) break;
	}
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

// detect variations for genome
int Genome::genomeDetect(){
	Time time;

	mkdir(out_dir_detect.c_str(), S_IRWXU | S_IROTH);  // create the directory for detect command
	Chrome *chr;
	for(size_t i=0; i<chromeVector.size(); i++){
		chr = chromeVector.at(i);
		//if(chr->chrname.compare("chr1")==0)
		{
			cout << "[" << time.getTime() << "]: processing Chr: " << chr->chrname << ", size: " << chr->chrlen << " bp" << endl;
			chr->chrDetect();
		}
	}

	// remove redundant translocations
	removeRedundantTra();

	// remove overlapped indels from mate clipping regions
	removeOverlappedIndelFromMateClipReg();

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
				if(clip_reg->leftClipReg) delete clip_reg->leftClipReg;
				if(clip_reg->leftClipReg2) delete clip_reg->leftClipReg2;
				if(clip_reg->rightClipReg) delete clip_reg->rightClipReg;
				if(clip_reg->rightClipReg2) delete clip_reg->rightClipReg2;
				if(clip_reg->var_cand) chr->removeVarCandNodeClipReg(clip_reg->var_cand);  // free item
				if(clip_reg->left_var_cand_tra) chr->removeVarCandNodeClipReg(clip_reg->left_var_cand_tra);  // free item
				if(clip_reg->right_var_cand_tra) chr->removeVarCandNodeClipReg(clip_reg->right_var_cand_tra);  // free item
				delete clip_reg;
				chr->mateClipRegVector.erase(chr->mateClipRegVector.begin()+j);
			}else j++;
		}
	}
}

void Genome::removeOverlappedIndelFromMateClipReg(){
	size_t i, j;
	Chrome *chr, *chr_tmp;
	for(i=0; i<chromeVector.size(); i++){
		chr = chromeVector.at(i);
		for(j=0; j<chromeVector.size(); j++){
			if(j!=i){
				chr_tmp = chromeVector.at(j);
				chr->removeFPIndelSnvInClipReg(chr_tmp->mateClipRegVector);
			}
		}
	}
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

	vector<Chrome*>::iterator chr;
	for(chr=chromeVector.begin(); chr!=chromeVector.end(); chr++){
		copySingleFile((*chr)->out_filename_detect_indel, out_file_indel); // indel
		copySingleFile((*chr)->out_filename_detect_snv, out_file_snv); // snv
		copySingleFile((*chr)->out_filename_detect_clipReg, out_file_clipReg); // clip regions
	}

	out_file_snv.close();
	out_file_indel.close();
	out_file_clipReg.close();
}

// local assembly for genome
int Genome::genomeLocalAssemble(){
	mkdir(out_dir_assemble.c_str(), S_IRWXU | S_IROTH);  // create directory for assemble command
	Chrome *chr;
	Time time;

	for(size_t i=0; i<chromeVector.size(); i++){
		chr = chromeVector.at(i);
		//if(chr->chrname.compare("1")==0)
		{
			chr->chrLoadDataAssemble();  // load the variation data
			chr->chrGenerateLocalAssembleWorkOpt();     // generate local assemble work
		}
	}

	cout << "Number of previously assembled regions: " << paras->assemble_reg_preDone_num << endl;
	cout << "Number of regions to be assembled: " << paras->assemble_reg_work_total << endl;

	//outputAssemWorkOptToFile(paras->assem_work_vec);

	// begin assemble
	cout << "[" << time.getTime() << "]: start local assemble..." << endl;
	processAssembleWork();

	computeVarNumStatAssemble(); // compute statistics for assemble command

	if(!paras->assem_work_vec.empty()) destroyAssembleWorkOptVec(paras->assem_work_vec);

	// reset assemble data
	for(size_t i=0; i<chromeVector.size(); i++){
		chr = chromeVector.at(i);
		chr->chrResetAssembleData();
	}
	return 0;
}


// process assemble work using thread pool
int Genome::processAssembleWork(){
	assembleWork_opt *assem_work_opt;
	assembleWork *assem_work;
	ofstream *var_cand_file;
	size_t num_work, num_work_per_ten_percent;

	hts_tpool *p = hts_tpool_init(paras->num_threads);
	hts_tpool_process *q = hts_tpool_process_init(p, paras->num_threads*2, 1);

	pthread_mutex_init(&paras->mtx_assemble_reg_workDone_num, NULL);

	paras->assemble_reg_workDone_num = 0;
	num_work = paras->assem_work_vec.size();
	num_work_per_ten_percent = num_work / paras->num_parts_progress;
	for(size_t i=0; i<num_work; i++){
		assem_work_opt = paras->assem_work_vec.at(i);
		var_cand_file = getVarcandFile(assem_work_opt->chrname, chromeVector, assem_work_opt->clip_reg_flag);
		if(var_cand_file==NULL){
			cerr << __func__ << ", line=" << __LINE__ << ": cannot get car_cand file for CHR: " << assem_work_opt->chrname << ", error!" << endl;
			exit(1);
		}

		assem_work = new assembleWork();
		assem_work->assem_work_opt = assem_work_opt;
		assem_work->work_id = i;
		assem_work->num_work = num_work;
		assem_work->num_work_per_ten_percent = num_work_per_ten_percent;
		assem_work->p_assemble_reg_workDone_num = &(paras->assemble_reg_workDone_num);
		assem_work->p_mtx_assemble_reg_workDone_num = &(paras->mtx_assemble_reg_workDone_num);
		assem_work->num_threads_per_assem_work = paras->num_threads_per_assem_work;
		assem_work->inBamFile = paras->inBamFile;
		assem_work->fai = fai;
		assem_work->var_cand_file = var_cand_file;
		assem_work->expected_cov_assemble = paras->expected_cov_assemble;
		assem_work->delete_reads_flag = paras->delete_reads_flag;

		hts_tpool_dispatch(p, q, processSingleAssembleWork, assem_work);
	}

    hts_tpool_process_flush(q);
    hts_tpool_process_destroy(q);
    hts_tpool_destroy(p);

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

// call variants for genome
int Genome::genomeCall(){
	size_t i;
	Chrome *chr;
	Time time;

	mkdir(out_dir_call.c_str(), S_IRWXU | S_IROTH);  // create directory for call command

	cout << "Begin calling variants ..." << endl;

	// call variants
	for(i=0; i<chromeVector.size(); i++){
		chr = chromeVector.at(i);
		//if(chr->chrname.compare("chr7")==0)
		{
			chr->chrLoadDataCall();
			chr->chrCall();
		}
	}

	// call TRA according to mate clip regions
	genomeCallTra();

	// fill variation sequence
	cout << "4444444444444" << endl;
	cout << "[" << time.getTime() << "]: fill variant sequences... " << endl;
	genomeFillVarseq();

	// save SV to file
	cout << "5555555555555" << endl;
	cout << "[" << time.getTime() << "]: save call result to file... " << endl;
	genomeSaveCallSV2File();

	// merge call results into single file
	cout << "6666666666666" << endl;
	cout << "[" << time.getTime() << "]: merge call result... " << endl;
	mergeCallResult();

	// compute statistics for call command
	cout << "[" << time.getTime() << "]: compute variant NUMBER statistics... " << endl;
	computeVarNumStatCall();

	return 0;
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
	reg_t *reg, *reg1, *reg2;
	string chrname_tmp;
	bool overlap_falg;

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
			reg1 = reg2 = NULL;
			if(clip_reg->valid_flag and clip_reg->call_success_flag==false and clip_reg->sv_type==VAR_TRA){ // valid TRA
				for(round_num=0; round_num<2; round_num++){
					if(round_num==0){
						var_cand = clip_reg->left_var_cand_tra;
						chrname_tmp = clip_reg->leftClipReg ? clip_reg->leftClipReg->chrname : clip_reg->leftClipReg2->chrname;
					}else{
						var_cand = clip_reg->right_var_cand_tra;
						chrname_tmp = clip_reg->rightClipReg ? clip_reg->rightClipReg->chrname : clip_reg->rightClipReg2->chrname;
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
									overlap_falg = isOverlappedPos(var_cand->leftClipRefPos, var_cand->rightClipRefPos, seg1->ref_end, seg2->ref_start);
									if(overlap_falg){ // indel
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
										reg->blat_aln_id = i;
										reg->call_success_status = true;
										reg->short_sv_flag = false;
										reg->zero_cov_flag = false;

										ref_dist = reg->endLocalRefPos - reg->startLocalRefPos + 1;
										query_dist = reg->endQueryPos - reg->startQueryPos + 1;
										reg->sv_len = query_dist - ref_dist + 1;
										if(reg->sv_len>0) reg->var_type = VAR_INS;
										else reg->var_type = VAR_DEL;

										if(round_num==0) reg1 = reg;
										else reg2 = reg;

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
				if(reg2){
					// right part
					clip_reg->right_var_cand_tra->varVec.push_back(reg2);

					if(clip_reg->rightClipReg) { delete clip_reg->rightClipReg; clip_reg->rightClipReg = NULL; clip_reg->rightMeanClipPos = 0; clip_reg->rightClipPosNum = 0; }
					if(clip_reg->rightClipReg2) { delete clip_reg->rightClipReg2; clip_reg->rightClipReg2 = NULL; clip_reg->rightMeanClipPos2 = 0; clip_reg->rightClipPosNum2 = 0; }
					clip_reg->rightClipRegNum = 0;
					clip_reg->right_var_cand_tra = NULL;
				}

				if(reg1 or reg2){
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

// call TRA according to mate clip regions
void Genome::genomeCallTra(){
	Time time;
	vector<blatAlnTra*> blat_aln_tra_vec;

	cout << "[" << time.getTime() << "]: call genomic translocations ..." << endl;

	mkdir(out_dir_tra.c_str(), S_IRWXU | S_IROTH);  // create the directory for TRA

	generateBlatAlnFilenameTra();
	blat_aln_tra_vec = loadBlatAlnDataTra();

	if(paras->num_threads<=1){
		cout << "Blat alignment using single thread ..." << endl;
		blatAlnTra_st(&blat_aln_tra_vec);
		cout << "Blat alignment using single thread finished." << endl;
	}else{
		// blat aln using multiple threads
		cout << "Blat alignment using " << paras->num_threads << " threads ..." << endl;
		blatAlnTra_mt(&blat_aln_tra_vec);
		cout << "Blat alignment using multiple threads finished." << endl;
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
	aln_file_tra.close();

	//cout << "File generated." << endl;
}

// call TRA according to mate clip regions
void Genome::genomeCallTraOp(){
	size_t i, j, k, t, round_num, success_num, total;
	Chrome *chr;
	vector<mateClipReg_t*> mate_clipReg_vec;
	mateClipReg_t *clip_reg;
	varCand *var_cand, *var_cand_tmp;
	blat_aln_t *item_blat;
	vector<int32_t> tra_loc_vec;

	total = 0;
	for(i=0; i<chromeVector.size(); i++){
		success_num = 0;
		chr = chromeVector.at(i);
		mate_clipReg_vec = chr->mateClipRegVector;
		for(j=0; j<mate_clipReg_vec.size(); j++){
//			if(i>=23 and j>=31)
//				cout << "gggggggggggggggggggg, i=" << i << ", j=" << j << ", " << chr->chrname << endl;

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
						// load align data
						var_cand->loadBlatAlnData();
						var_cand_tmp->loadBlatAlnData();

						// compute TRA locations
						tra_loc_vec = computeTraLoc(var_cand, var_cand_tmp, clip_reg);

						if(tra_loc_vec.size()>0){
							saveTraLoc2ClipReg(clip_reg, tra_loc_vec, var_cand, var_cand_tmp, round_num); // save TRA location to clip region
							clip_reg->call_success_flag = true;
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
			blat_aln_tra = new blatAlnTra(str_vec.at(0), str_vec.at(1), str_vec.at(2));
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
	MultiThread mt[paras->num_threads];
	for(size_t i=0; i<paras->num_threads; i++){
		mt[i].setNumThreads(paras->num_threads);
		mt[i].setBlatAlnTraVec(blat_aln_tra_vec);
		mt[i].setUserThreadID(i);
		if(!mt[i].startBlatAlnTra()){
			cerr << __func__ << ", line=" << __LINE__ << ": unable to create thread, error!" << endl;
			exit(1);
		}
	}
	for(size_t i=0; i<paras->num_threads; i++){
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

		var_cand_new->refseqfilename = var_cand_tmp->refseqfilename;  // refseq file name
		var_cand_new->ctgfilename = var_cand->ctgfilename;  // contig file name
		var_cand_new->readsfilename = var_cand->readsfilename;  // reads file name
		var_cand_new->ref_left_shift_size = var_cand_tmp->ref_left_shift_size;  // ref_left_shift_size
		var_cand_new->ref_right_shift_size = var_cand_tmp->ref_right_shift_size;  // ref_right_shift_size

		var_cand_new->blat_aligned_info_vec = NULL;
		var_cand_new->blat_var_cand_file = NULL;

		var_cand_new->assem_success = var_cand->assem_success;
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
	}

	return var_cand_new;
}

// compute variant location
vector<int32_t> Genome::computeTraLoc(varCand *var_cand, varCand *var_cand_tmp, mateClipReg_t *clip_reg){
	size_t i, j, aln_orient, missing_size_head, missing_size_inner, missing_size_tail;
	int32_t start_query_pos, end_query_pos, start_mate_query_pos, end_mate_query_pos;
	int32_t maxValue, maxIdx, secValue, query_dist, query_part_flag_maxValue, query_part_flag_secValue;  // query_part_flag: 0-head, 1-inner, 2-tail
	int32_t tra_blat_aln_idx, start_tra_aln_seg_idx, end_tra_aln_seg_idx, value_tmp, query_part_flag;
	blat_aln_t *blat_aln, *blat_aln_mate;
	aln_seg_t *aln_seg1, *aln_seg2;
	vector<int32_t> aln_idx_vec, mate_aln_idx_vec, tra_loc_vec, tra_loc_vec2;
	bool call_success_flag;

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
						aln_seg1 = blat_aln_mate->aln_segs.at(mate_aln_idx_vec.at(1));
						aln_seg2 = blat_aln_mate->aln_segs.at(mate_aln_idx_vec.at(2));
						if(blat_aln_mate->aln_orient==ALN_PLUS_ORIENT){ // plus orient
							start_mate_query_pos = aln_seg1->query_start;
							end_mate_query_pos = aln_seg2->query_end;
						}else{ // minus orient
							start_mate_query_pos = blat_aln_mate->query_len - aln_seg1->query_end + 1;
							end_mate_query_pos = blat_aln_mate->query_len - aln_seg2->query_start + 1;
						}

						//cout << "mate TRA: " << "start_mate_query_pos=" << start_mate_query_pos << ", end_mate_query_pos=" << end_mate_query_pos << endl;

						// compute TRA clipping locations
						tra_loc_vec = computeTraClippingLoc(query_part_flag, var_cand, aln_idx_vec, var_cand_tmp, mate_aln_idx_vec);

						//cout << "TRA location: " << var_cand->chrname << ":" << tra_loc_vec.at(0) << "-" << tra_loc_vec.at(1) << ", " << var_cand_tmp->chrname << ":" << tra_loc_vec.at(2) << "-" << tra_loc_vec.at(3) << endl;

						call_success_flag = computeTraCallSuccessFlag(tra_loc_vec, tra_loc_vec2, var_cand, var_cand_tmp, clip_reg);

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

	return tra_loc_vec;
}

// get the missing region in the mated clipping region: [0]-blat_aln_id, [1]-start_aln_seg_id, [2]-end_aln_seg_id
vector<int32_t> Genome::getMateTraReg(size_t query_id, size_t start_query_pos, size_t end_query_pos, varCand *var_cand_tmp){
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

	left_tra_refPos1 = right_tra_refPos1 = left_tra_refPos2 = right_tra_refPos2 = -1;
	blat_aln = var_cand->blat_aln_vec.at(aln_idx_vec.at(0));
	blat_aln_tmp = var_cand_tmp->blat_aln_vec.at(mate_aln_idx_vec.at(0));

	if(query_clip_part_flag==0){ // query head
		aln_seg1 = blat_aln->aln_segs.at(aln_idx_vec.at(2));
		left_tra_refPos1 = -1;
		if(blat_aln->aln_orient==ALN_PLUS_ORIENT)
			right_tra_refPos1 = aln_seg1->ref_start;
		else
			right_tra_refPos1 = aln_seg1->ref_end;
	}else if(query_clip_part_flag==1){ // query inner
		aln_seg1 = blat_aln->aln_segs.at(aln_idx_vec.at(1));
		aln_seg2 = blat_aln->aln_segs.at(aln_idx_vec.at(2));
		if(blat_aln->aln_orient==ALN_PLUS_ORIENT){
			left_tra_refPos1 = aln_seg1->ref_end;
			right_tra_refPos1 = aln_seg2->ref_start;
		}else{
			left_tra_refPos1 = aln_seg1->ref_start;
			right_tra_refPos1 = aln_seg2->ref_end;
		}
	}else if(query_clip_part_flag==2){ // query tail
		aln_seg2 = blat_aln->aln_segs.at(aln_idx_vec.at(1));
		if(blat_aln->aln_orient==ALN_PLUS_ORIENT)
			left_tra_refPos1 = aln_seg2->ref_end;
		else
			left_tra_refPos1 = aln_seg2->ref_start;
		right_tra_refPos1 = -1;
	}
	//cout << "left_tra_refPos1=" << left_tra_refPos1 << ", right_tra_refPos1=" << right_tra_refPos1 << endl;

	if(query_clip_part_flag==0){ // query head
		aln_seg2 = blat_aln_tmp->aln_segs.at(mate_aln_idx_vec.at(2));
		left_tra_refPos2 = -1;
		if(blat_aln_tmp->aln_orient==ALN_PLUS_ORIENT){ // plus orient
			right_tra_refPos2 = aln_seg2->ref_end;
		}else{ // minus orient
			right_tra_refPos2 = aln_seg2->ref_start;
		}
	}else if(query_clip_part_flag==1){ // query inner
		aln_seg1 = blat_aln_tmp->aln_segs.at(mate_aln_idx_vec.at(1));
		aln_seg2 = blat_aln_tmp->aln_segs.at(mate_aln_idx_vec.at(2));
		if(blat_aln_tmp->aln_orient==ALN_PLUS_ORIENT){ // plus orient
			left_tra_refPos2 = aln_seg1->ref_start;
			right_tra_refPos2 = aln_seg2->ref_end;
		}else{ // minus orient
			left_tra_refPos2 = aln_seg1->ref_end;
			right_tra_refPos2 = aln_seg2->ref_start;
		}
	}else if(query_clip_part_flag==2){ // query tail
		aln_seg1 = blat_aln_tmp->aln_segs.at(mate_aln_idx_vec.at(1));
		if(blat_aln_tmp->aln_orient==ALN_PLUS_ORIENT){ // plus orient
			left_tra_refPos2 = aln_seg1->ref_start;
		}else{ // minus orient
			left_tra_refPos2 = aln_seg1->ref_end;
		}
		right_tra_refPos2 = -1;
	}
	//cout << "left_tra_refPos2=" << left_tra_refPos2 << ", right_tra_refPos2=" << right_tra_refPos2 << endl;

	tra_loc_vec.push_back(left_tra_refPos1);
	tra_loc_vec.push_back(right_tra_refPos1);
	tra_loc_vec.push_back(left_tra_refPos2);
	tra_loc_vec.push_back(right_tra_refPos2);

	return tra_loc_vec;
}


// determine whether the TRA is called successfully according to TRA location vector and mate cliping information
bool Genome::computeTraCallSuccessFlag(vector<int32_t> &tra_loc_vec, vector<int32_t> &tra_loc_vec2, varCand *var_cand, varCand *var_cand_tmp, mateClipReg_t *clip_reg){
	bool flag = false, overlap_flag1, overlap_flag2, overlap_flag3, overlap_flag4;
	size_t i, tmp;
	reg_t *reg1_clipReg, *reg2_clipReg;
	int32_t ref_dist1, ref_dist2, tra_loc1, tra_loc2, tra_loc3, tra_loc4, tra_loc1_clipReg, tra_loc2_clipReg, tra_loc3_clipReg, tra_loc4_clipReg;
	double ratio;

	if(tra_loc_vec.size()>0){
		if(clip_reg->leftClipRegNum==1 and clip_reg->rightClipRegNum==1){
			if((tra_loc_vec.at(0)!=-1 or tra_loc_vec.at(1)!=-1) and (tra_loc_vec.at(2)!=-1 or tra_loc_vec.at(3)!=-1)){ // in correct region
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
					tra_loc2 = tra_loc_vec.at(0);
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
	if(tra_loc_vec.at(0)!=-1 and tra_loc_vec.at(1)!=-1 and tra_loc_vec.at(0)>tra_loc_vec.at(1)){  // INV_TRA
		if(round_num==0){
			clip_reg->chrname_leftTra1 = var_cand->chrname;
			clip_reg->leftClipPosTra1 = tra_loc_vec.at(1);
			clip_reg->chrname_rightTra1 = var_cand->chrname;
			clip_reg->rightClipPosTra1 = tra_loc_vec.at(0);
		}else{
			clip_reg->chrname_leftTra2 = var_cand->chrname;
			clip_reg->leftClipPosTra2 = tra_loc_vec.at(1);
			clip_reg->chrname_rightTra2 = var_cand->chrname;
			clip_reg->rightClipPosTra2 = tra_loc_vec.at(0);
		}
	}else{  // TRA
		if(round_num==0){
			if(tra_loc_vec.at(0)!=-1) clip_reg->chrname_leftTra1 = var_cand->chrname;
			clip_reg->leftClipPosTra1 = tra_loc_vec.at(0);
			if(tra_loc_vec.at(1)!=-1) clip_reg->chrname_rightTra1 = var_cand->chrname;
			clip_reg->rightClipPosTra1 = tra_loc_vec.at(1);
		}else{
			if(tra_loc_vec.at(0)!=-1) clip_reg->chrname_leftTra2 = var_cand->chrname;
			clip_reg->leftClipPosTra2 = tra_loc_vec.at(0);
			if(tra_loc_vec.at(1)!=-1) clip_reg->chrname_rightTra2 = var_cand->chrname;
			clip_reg->rightClipPosTra2 = tra_loc_vec.at(1);
		}
	}
	if(tra_loc_vec.at(2)!=-1 and tra_loc_vec.at(3)!=-1 and tra_loc_vec.at(2)>tra_loc_vec.at(3)){  // INV_TRA
		if(round_num==0){
			clip_reg->chrname_leftTra2 = var_cand_tmp->chrname;
			clip_reg->leftClipPosTra2 = tra_loc_vec.at(3);
			clip_reg->chrname_rightTra2 = var_cand_tmp->chrname;
			clip_reg->rightClipPosTra2 = tra_loc_vec.at(2);
		}else{
			clip_reg->chrname_leftTra1 = var_cand_tmp->chrname;
			clip_reg->leftClipPosTra1 = tra_loc_vec.at(3);
			clip_reg->chrname_rightTra1 = var_cand_tmp->chrname;
			clip_reg->rightClipPosTra1 = tra_loc_vec.at(2);
		}
	}else{  // TRA
		if(round_num==0){
			if(tra_loc_vec.at(2)!=-1) clip_reg->chrname_leftTra2 = var_cand_tmp->chrname;
			clip_reg->leftClipPosTra2 = tra_loc_vec.at(2);
			if(tra_loc_vec.at(3)!=-1) clip_reg->chrname_rightTra2 = var_cand_tmp->chrname;
			clip_reg->rightClipPosTra2 = tra_loc_vec.at(3);
		}else{
			if(tra_loc_vec.at(2)!=-1) clip_reg->chrname_leftTra1 = var_cand_tmp->chrname;
			clip_reg->leftClipPosTra1 = tra_loc_vec.at(2);
			if(tra_loc_vec.at(3)!=-1) clip_reg->chrname_rightTra1 = var_cand_tmp->chrname;
			clip_reg->rightClipPosTra1 = tra_loc_vec.at(3);
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

				if(min_dist1<(int32_t)paras->maxClipRegSize and min_dist2<(int32_t)paras->maxClipRegSize){
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
	Chrome *chr;
	for(size_t i=0; i<chromeVector.size(); i++){
		chr = chromeVector.at(i);
		//if(chr->chrname.compare("12")==0)
			chr->chrFillVarseq();
	}

	// fill sequences for TRA
	genomeFillVarseqTra();
}

void Genome::genomeFillVarseqTra(){
	size_t i, j;
	Chrome *chr;
	mateClipReg_t *clip_reg;
	ofstream assembly_info_file;
	string assembly_info_filename_tra;

	assembly_info_filename_tra = out_dir_tra + "/" + "tra_assembly_info";
	assembly_info_file.open(assembly_info_filename_tra);
	if(!assembly_info_file.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << assembly_info_filename_tra << endl;
		exit(1);
	}

	for(i=0; i<chromeVector.size(); i++){
		chr = chromeVector.at(i);
		for(j=0; j<chr->mateClipRegVector.size(); j++){
			clip_reg = chr->mateClipRegVector.at(j);
			if(clip_reg->valid_flag and clip_reg->call_success_flag and clip_reg->sv_type==VAR_TRA)
				fillVarseqSingleMateClipReg(clip_reg, assembly_info_file);
		}
	}

	assembly_info_file.close();
}

void Genome::fillVarseqSingleMateClipReg(mateClipReg_t *clip_reg, ofstream &assembly_info_file){
	varCand *var_cand_tmp;
	string tmpdir;
	reg_t *reg;
	vector<int32_t> ref_shift_size_vec;
	vector<size_t> query_loc_vec; // [0]-blat_aln_id, [1]-query_id, [2]-start_query_pos, [3]-end_query_pos, [4]-start_local_ref_pos, [5]-end_loal_ref_pos, [6]-start_ref_pos, [7]-end_ref_pos
	blat_aln_t *blat_aln;
	size_t i, assembly_extend_size;

	if(clip_reg->valid_flag){
		if(clip_reg->sv_type==VAR_TRA){
			if(clip_reg->leftClipPosTra1!=-1 and clip_reg->rightClipPosTra1!=-1 and clip_reg->chrname_leftTra1.compare(clip_reg->chrname_rightTra1)==0){

				// construct var_cand
				var_cand_tmp = new varCand();

				var_cand_tmp->chrname = clip_reg->chrname_leftTra1;
				var_cand_tmp->var_cand_filename = "";
				var_cand_tmp->out_dir_call = out_dir_tra;
				var_cand_tmp->misAln_filename = "";
				var_cand_tmp->inBamFile = paras->inBamFile;
				var_cand_tmp->fai = fai;
				var_cand_tmp->call_success = false;

				var_cand_tmp->blat_aligned_info_vec = NULL;
				var_cand_tmp->blat_var_cand_file = NULL;

				var_cand_tmp->leftClipRefPos = clip_reg->leftClipPosTra1;
				var_cand_tmp->rightClipRefPos = clip_reg->rightClipPosTra1;

				reg = new reg_t();
				reg->chrname = clip_reg->chrname_leftTra1;
				reg->startRefPos = clip_reg->leftClipPosTra1;
				reg->endRefPos = clip_reg->rightClipPosTra1;
				reg->var_type = clip_reg->sv_type;
				reg->call_success_status = false;
				reg->short_sv_flag = false;
				reg->zero_cov_flag = false;

				var_cand_tmp->varVec.push_back(reg);

				// construct the variant region
				var_cand_tmp->readsfilename = out_dir_tra  + "/" + "tra_reads_" + reg->chrname + "_" + to_string(reg->startRefPos) + "-" + to_string(reg->endRefPos) + ".fq";
				var_cand_tmp->ctgfilename = out_dir_tra  + "/" + "tra_contig_" + reg->chrname + "_" + to_string(reg->startRefPos) + "-" + to_string(reg->endRefPos) + ".fa";
				var_cand_tmp->refseqfilename = out_dir_tra  + "/" + "tra_refseq_" + reg->chrname + "_" + to_string(reg->startRefPos) + "-" + to_string(reg->endRefPos) + ".fa";
				var_cand_tmp->alnfilename = out_dir_tra  + "/" + "tra_blat_" + reg->chrname + "_" + to_string(reg->startRefPos) + "-" + to_string(reg->endRefPos) + ".sim4";
				tmpdir = out_dir_tra + "/" + "tmp_tra_" + reg->chrname + "_" + to_string(reg->startRefPos) + "-" + to_string(reg->endRefPos);

				for(i=0; i<3; i++){
					assembly_extend_size = ASSEMBLY_SIDE_LEN * i;
					// local assembly
					performLocalAssemblyTra(var_cand_tmp->readsfilename, var_cand_tmp->ctgfilename, var_cand_tmp->refseqfilename, tmpdir, paras->num_threads_per_assem_work, var_cand_tmp->varVec, reg->chrname, paras->inBamFile, fai, assembly_extend_size, assembly_info_file);

					ref_shift_size_vec = getRefShiftSize(var_cand_tmp->refseqfilename);
					var_cand_tmp->ref_left_shift_size = ref_shift_size_vec.at(0);
					var_cand_tmp->ref_right_shift_size = ref_shift_size_vec.at(1);

					var_cand_tmp->sv_type = clip_reg->sv_type;
					var_cand_tmp->dup_num = clip_reg->dup_num;
					var_cand_tmp->margin_adjusted_flag = false;

					// assembly status
					struct stat fileStat;
					if (stat(var_cand_tmp->ctgfilename.c_str(), &fileStat) == 0)
						var_cand_tmp->assem_success = true;
					else
						var_cand_tmp->assem_success = false;
					if(var_cand_tmp->assem_success){

						var_cand_tmp->loadBlatAlnData(); // BLAT align
						if(var_cand_tmp->align_success){
							// compute the aligned query locations for TRA
							query_loc_vec = computeQueryLocTra(var_cand_tmp, clip_reg, LEFT_END);
							if(query_loc_vec.size()>0){
								clip_reg->leftClipPosTra1 = query_loc_vec.at(6);
								clip_reg->rightClipPosTra1 = query_loc_vec.at(7);

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

			if(clip_reg->leftClipPosTra2!=-1 and clip_reg->rightClipPosTra2!=-1 and clip_reg->chrname_leftTra2.compare(clip_reg->chrname_rightTra2)==0){
				// construct var_cand
				var_cand_tmp = new varCand();

				var_cand_tmp->chrname = clip_reg->chrname_leftTra2;
				var_cand_tmp->var_cand_filename = "";
				var_cand_tmp->out_dir_call = out_dir_tra;
				var_cand_tmp->misAln_filename = "";
				var_cand_tmp->inBamFile = paras->inBamFile;
				var_cand_tmp->fai = fai;
				var_cand_tmp->align_success = false;

				var_cand_tmp->blat_aligned_info_vec = NULL;
				var_cand_tmp->blat_var_cand_file = NULL;

				var_cand_tmp->leftClipRefPos = clip_reg->leftClipPosTra2;
				var_cand_tmp->rightClipRefPos = clip_reg->rightClipPosTra2;

				reg = new reg_t();
				reg->chrname = clip_reg->chrname_leftTra2;
				reg->startRefPos = clip_reg->leftClipPosTra2;
				reg->endRefPos = clip_reg->rightClipPosTra2;
				reg->var_type = clip_reg->sv_type;
				reg->call_success_status = false;
				reg->short_sv_flag = false;
				reg->zero_cov_flag = false;

				var_cand_tmp->varVec.push_back(reg);

				// construct the variant region
				var_cand_tmp->readsfilename = out_dir_tra  + "/" + "tra_reads_" + reg->chrname + "_" + to_string(reg->startRefPos) + "-" + to_string(reg->endRefPos) + ".fq";
				var_cand_tmp->ctgfilename = out_dir_tra  + "/" + "tra_contig_" + reg->chrname + "_" + to_string(reg->startRefPos) + "-" + to_string(reg->endRefPos) + ".fa";
				var_cand_tmp->refseqfilename = out_dir_tra  + "/" + "tra_refseq_" + reg->chrname + "_" + to_string(reg->startRefPos) + "-" + to_string(reg->endRefPos) + ".fa";
				var_cand_tmp->alnfilename = out_dir_tra  + "/" + "tra_blat_" + reg->chrname + "_" + to_string(reg->startRefPos) + "-" + to_string(reg->endRefPos) + ".sim4";
				tmpdir = out_dir_tra + "/" + "tmp_tra_" + reg->chrname + "_" + to_string(reg->startRefPos) + "-" + to_string(reg->endRefPos);

				for(i=0; i<3; i++){
					assembly_extend_size = ASSEMBLY_SIDE_LEN * i;
					// local assembly
					performLocalAssemblyTra(var_cand_tmp->readsfilename, var_cand_tmp->ctgfilename, var_cand_tmp->refseqfilename, tmpdir, paras->num_threads_per_assem_work, var_cand_tmp->varVec, reg->chrname, paras->inBamFile, fai, assembly_extend_size, assembly_info_file);

					ref_shift_size_vec = getRefShiftSize(var_cand_tmp->refseqfilename);
					var_cand_tmp->ref_left_shift_size = ref_shift_size_vec.at(0);
					var_cand_tmp->ref_right_shift_size = ref_shift_size_vec.at(1);

					var_cand_tmp->sv_type = clip_reg->sv_type;
					var_cand_tmp->dup_num = clip_reg->dup_num;
					var_cand_tmp->margin_adjusted_flag = false;

					// assembly status
					struct stat fileStat;
					if (stat(var_cand_tmp->ctgfilename.c_str(), &fileStat) == 0)
						var_cand_tmp->assem_success = true;
					else
						var_cand_tmp->assem_success = false;

					if(var_cand_tmp->assem_success){

						var_cand_tmp->loadBlatAlnData(); // BLAT align
						if(var_cand_tmp->align_success){
							// compute the aligned query locations for TRA
							query_loc_vec = computeQueryLocTra(var_cand_tmp, clip_reg, RIGHT_END);
							if(query_loc_vec.size()>0){
								clip_reg->leftClipPosTra2 = query_loc_vec.at(6);
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

// perform local assembly
void Genome::performLocalAssemblyTra(string &readsfilename, string &contigfilename, string &refseqfilename, string &tmpdir, size_t num_threads_per_assem_work, vector<reg_t*> &varVec, string &chrname, string &inBamFile, faidx_t *fai, size_t assembly_extend_size, ofstream &assembly_info_file){

	LocalAssembly local_assembly(readsfilename, contigfilename, refseqfilename, tmpdir, num_threads_per_assem_work, varVec, chrname, inBamFile, fai, assembly_extend_size, paras->expected_cov_assemble, paras->delete_reads_flag);

	// extract the corresponding refseq from reference
	local_assembly.extractRefseq();

	// extract the reads data from BAM file
	local_assembly.extractReadsDataFromBAM();

	// local assembly using Canu
	local_assembly.localAssembleCanu();

	// record assembly information
	local_assembly.recordAssemblyInfo(assembly_info_file);

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

	if(end_flag==LEFT_END) { start_ref_pos_tra = clip_reg->leftClipPosTra1; end_ref_pos_tra = clip_reg->rightClipPosTra1; }
	else { start_ref_pos_tra = clip_reg->leftClipPosTra2; end_ref_pos_tra = clip_reg->rightClipPosTra2; }

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
	mkdir(out_dir_result.c_str(), S_IRWXU | S_IROTH);  // create the directory for final results

	Chrome *chr;
	for(size_t i=0; i<chromeVector.size(); i++){
		chr = chromeVector.at(i);
		chr->saveCallSV2File();
	}

	// save TRA
	saveTraCall2File();
}

// save TRA
void Genome::saveTraCall2File(){
	size_t i, j;
	Chrome *chr;
	ofstream outfile_tra;
	vector<mateClipReg_t*> mate_clipReg_vec;
	mateClipReg_t *clip_reg;
	string line, header_line_bedpe;;
	int32_t sv_len;

	outfile_tra.open(out_filename_result_tra);
	if(!outfile_tra.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << out_filename_result_tra << endl;
		exit(1);
	}

	header_line_bedpe = getCallFileHeaderBedpe();
	outfile_tra << header_line_bedpe << endl;

	for(i=0; i<chromeVector.size(); i++){
		chr = chromeVector.at(i);
		mate_clipReg_vec = chr->mateClipRegVector;
		for(j=0; j<mate_clipReg_vec.size(); j++){
			clip_reg = mate_clipReg_vec.at(j);
			if(clip_reg->valid_flag and clip_reg->call_success_flag and clip_reg->sv_type==VAR_TRA and (clip_reg->leftClipPosTra1>0 or clip_reg->rightClipPosTra1>0) and (clip_reg->leftClipPosTra2>0 or clip_reg->rightClipPosTra2>0)){
				line = clip_reg->left_var_cand_tra->chrname;
				if(clip_reg->leftClipPosTra1>0) line += "\t" + to_string(clip_reg->leftClipPosTra1);
				else line += "\t-";
				if(clip_reg->rightClipPosTra1>0) line += "\t" + to_string(clip_reg->rightClipPosTra1);
				else line += "\t-";
				line += "\t" + clip_reg->right_var_cand_tra->chrname;
				if(clip_reg->leftClipPosTra2>0) line += "\t" + to_string(clip_reg->leftClipPosTra2);
				else line += "\t-";
				if(clip_reg->rightClipPosTra2>0) line += "\t" + to_string(clip_reg->rightClipPosTra2);
				else line += "\t-";
				line += "\tTRA";

				if(clip_reg->leftClipPosTra1>0 and clip_reg->rightClipPosTra1>0){
					sv_len = clip_reg->rightClipPosTra1 - clip_reg->leftClipPosTra1 + 1;
					line += "\t" + to_string(sv_len);
				}else line += "\t-";
				if(clip_reg->leftClipPosTra2>0 and clip_reg->rightClipPosTra2>0){
					sv_len = clip_reg->rightClipPosTra2 - clip_reg->leftClipPosTra2 + 1;
					line += "\t" + to_string(sv_len);
				}else line += "\t-";

				if(clip_reg->refseq_tra.size()>0) line += "\t" + clip_reg->refseq_tra;
				else line += "\t-";
				if(clip_reg->altseq_tra.size()>0) line += "\t" + clip_reg->altseq_tra;
				else line += "\t-";
				if(clip_reg->refseq_tra2.size()>0) line += "\t" + clip_reg->refseq_tra2;
				else line += "\t-";
				if(clip_reg->altseq_tra2.size()>0) line += "\t" + clip_reg->altseq_tra2;
				else line += "\t-";

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
	header_line_bed = getCallFileHeaderBed();
	out_file_indel << header_line_bed << endl;
	out_file_clipReg << header_line_bed << endl;

	for(size_t i=0; i<chromeVector.size(); i++){
		copySingleFile(chromeVector.at(i)->out_filename_call_indel, out_file_indel); // indel
		copySingleFile(chromeVector.at(i)->out_filename_call_clipReg, out_file_clipReg); // clip_reg
	}
	out_file_indel.close();
	out_file_clipReg.close();

	// merge indels, clipping variants, translocations to single file
	header_line_bedpe = getCallFileHeaderBedpe();
	out_file_vars << header_line_bed << endl;
	out_file_vars << header_line_bedpe << endl;
	copySingleFile(out_filename_result_indel, out_file_vars); // indels
	copySingleFile(out_filename_result_clipReg, out_file_vars); // INV, DUP
	copySingleFile(out_filename_result_tra, out_file_vars); // TRA
	out_file_vars.close();
}

// save results in VCF file format for detect command
void Genome::saveResultVCF(){
	string outIndelFile, outSnvFile;

	cout << "Output results to file ..." << endl;

	// initialize the output filenames
	if(paras->outFilePrefix.size()>0) {
		outIndelFile = out_dir_call + "/" + paras->outFilePrefix + "_INDEL.vcf";
		outSnvFile = out_dir_call + "/" + paras->outFilePrefix + "_SNV.vcf";
	}else {
		outIndelFile = out_dir_call + "/" + "INDEL.vcf";
		outSnvFile = out_dir_call + "/" + "SNV.vcf";
	}

	// save results
	//saveIndelVCFDetect(out_filename_detect_indel, outIndelFile);
	//saveSnvVCFDetect(out_filename_detect_snv, outSnvFile);
}

// save results in VCF file format for detect command
void Genome::saveResultVCFDetect(){
	string outIndelFile, outSnvFile;

	cout << "Output results to file ..." << endl;

	// initialize the output filenames
	if(paras->outFilePrefix.size()>0) {
		outIndelFile = out_dir_detect + "/" + paras->outFilePrefix + "_INDEL.vcf";
		outSnvFile = out_dir_detect + "/" + paras->outFilePrefix + "_SNV.vcf";
	}else {
		outIndelFile = out_dir_detect + "/" + "INDEL.vcf";
		outSnvFile = out_dir_detect + "/" + "SNV.vcf";
	}

	// save results
	saveIndelVCFDetect(out_filename_detect_indel, outIndelFile);
	saveSnvVCFDetect(out_filename_detect_snv, outSnvFile);
}

// save indel in VCF file format for detect command
void Genome::saveIndelVCFDetect(string &in, string &out_vcf){
	string line, data_vcf, ref, reg;
	ifstream infile;
	ofstream outfile;
	vector<string> str_vec;
	char *seq;
	int seq_len;

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
	saveVCFHeader(outfile);

	// save results
	while(getline(infile, line))
		if(line.size()>0){
			str_vec = split(line, "\t");
			data_vcf = str_vec[0] + "\t" + str_vec[1] + "\t."; // CHROM, POS, ID

			reg = str_vec[0] + ":" + str_vec[1] + "-" + str_vec[2];
			seq = fai_fetch(fai, reg.c_str(), &seq_len);
			data_vcf += "\t" + (string)seq + "\t."; // REF, ALT

			// QUAL, FILTER, INFO0, FORMAT ........


			outfile << data_vcf << endl;

			free(seq);
		}
	infile.close();
	outfile.close();
}

// save SNV in VCF file format for detect command
void Genome::saveSnvVCFDetect(string &in, string &out_vcf){
	string line, data_vcf, ref, reg;
	ifstream infile;
	ofstream outfile;
	vector<string> str_vec;
	char *seq;
	int seq_len;

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
	saveVCFHeader(outfile);

	// save results
	while(getline(infile, line))
		if(line.size()>0){
			str_vec = split(line, "\t");
			data_vcf = str_vec[0] + "\t" + str_vec[1] + "\t."; // CHROM, POS, ID

			reg = str_vec[0] + ":" + str_vec[1] + "-" + str_vec[1];
			seq = fai_fetch(fai, reg.c_str(), &seq_len);
			data_vcf += "\t" + (string)seq + "\t."; // REF, ALT

			// QUAL, FILTER, INFO0, FORMAT ........


			outfile << data_vcf << endl;

			free(seq);
		}
	infile.close();
	outfile.close();
}

// save VCF header
void Genome::saveVCFHeader(ofstream &fp){
	time_t rawtime;
	struct tm *timeinfo;
	char buffer [128];

	// get the local time
	time(&rawtime);
	timeinfo = localtime (&rawtime);
	strftime(buffer, sizeof(buffer), "%Y%m%d", timeinfo);

	fp << "##fileformat=VCFv4.2" << endl;
	fp << "##fileDate=" << buffer << endl;
	fp << "##source=" << PROG_NAME << "V" << PROG_VERSION << endl;
	fp << "##reference=" << paras->refFile << endl;
	fp << "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">" << endl;
	fp << "##FILTER=<ID=q10,Description=\"Quality below 10\">" << endl;
	fp << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl;
	fp << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT" << endl;

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

// compute statistics for 'assemble' command
void Genome::computeVarNumStatAssemble(){
	int32_t total, total_indel, total_clipReg, total_succ_indel, total_fail_indel, total_succ_clipReg, total_fail_clipReg, succ_tmp, fail_tmp;
	vector<int32_t> num_vec;
	Chrome *chr;
	string filename;

	total_indel = total_clipReg = total_succ_indel = total_fail_indel = total_succ_clipReg = total_fail_clipReg = 0;
	for(size_t i=0; i<chromeVector.size(); i++){
		chr = chromeVector.at(i);

		filename = chr->getVarcandIndelFilename();
		num_vec = getSuccFailNumAssemble(filename);
		succ_tmp = num_vec.at(0);
		fail_tmp = num_vec.at(1);
		total_succ_indel += succ_tmp;
		total_fail_indel += fail_tmp;
		total_indel += succ_tmp + fail_tmp;

		filename = chr->getVarcandClipregFilename();
		num_vec = getSuccFailNumAssemble(filename);
		succ_tmp = num_vec.at(0);
		fail_tmp = num_vec.at(1);
		total_succ_clipReg += succ_tmp;
		total_fail_clipReg += fail_tmp;
		total_clipReg += succ_tmp + fail_tmp;
	}
	total = total_indel + total_clipReg;

	cout << "############## Brief statistics for local assembly ##############" << endl;
	cout << "There are " << total << " local assembled regions in total:" << endl;
	cout << "\t" << total_indel << " indel regions: " << total_succ_indel << " successful, " << total_fail_indel << " failed" << endl;
	cout << "\t" << total_clipReg << " clipping regions: " << total_succ_clipReg << " successful, " << total_fail_clipReg << " failed" << endl;
}

vector<int32_t> Genome::getSuccFailNumAssemble(string &filename){
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
				if(str_vec.at(5).compare(ASSEMBLY_SUCCESS)==0) num_succ ++;
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

	total = num_ins + num_del + num_inv + num_dup + num_tra;

	cout << "############## Brief statistics for variants call ##############" << endl;
	cout << "There are " << total << " variants in total:" << endl;
	cout << "\t" << "insertions: " << num_ins << endl;
	cout << "\t" << "deletions: " << num_del << endl;
	cout << "\t" << "tandem duplications: " << num_dup << endl;
	cout << "\t" << "inversions: " << num_inv << endl;
	cout << "\t" << "translocations: " << num_tra << endl;
	//cout << "\t" << "unresolved: " << num_unc << endl;

	infile.close();
}

size_t Genome::getSVTypeSingleLine(string &line){
	size_t sv_type = VAR_UNC;
	vector<string> str_vec;
	string str_tmp;

	str_vec = split(line, "\t");
	str_tmp = str_vec.at(3);
	if(str_tmp.compare("INS")==0 or str_tmp.compare("insertion")==0){
		sv_type = VAR_INS;
	}else if(str_tmp.compare("DEL")==0 or str_tmp.compare("deletion")==0){
		sv_type = VAR_DEL;
	}else if(str_tmp.compare("DUP")==0 or str_tmp.compare("duplication")==0){
		sv_type = VAR_DUP;
	}else if(str_tmp.compare("INV")==0 or str_tmp.compare("inversion")==0){
		sv_type = VAR_INV;
	}else{
		if(str_vec.size()>=8 and (str_vec.at(6).compare("TRA")==0 or str_vec.at(6).compare("translocation")==0 or str_vec.at(6).compare("BND")==0)){
			if(str_vec.at(6).compare("TRA")==0 or str_vec.at(6).compare("translocation")==0) sv_type = VAR_TRA;
			else sv_type = VAR_BND;
		}else{
			sv_type = VAR_MIX;
		}
	}

	return sv_type;
}

