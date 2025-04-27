#include <sys/stat.h>
#include <pthread.h>
#include "Chrome.h"
#include "Thread.h"
#include "util.h"
#include <pthread.h>
#include <htslib/thread_pool.h>

pthread_mutex_t mutex_mate_clip_reg = PTHREAD_MUTEX_INITIALIZER;
extern pthread_mutex_t mutex_fai;

// Constructor with parameters
Chrome::Chrome(string& chrname, int chrlen, faidx_t *fai, Paras *paras, vector<Chrome*> *chr_vec){
	this->paras = paras;
	this->chrname = chrname;
	this->chrlen = chrlen;
	this->fai = fai;
	this->chr_vec = chr_vec;
	init();

	blocks_out_file = chrname + "_blocks.bed";
	out_dir_detect = "";
	out_dir_cns = "";
	out_dir_call = "";

	out_filename_detect_snv = "";
	out_filename_detect_indel = "";
	out_filename_detect_clipReg = "";
	out_filename_call_snv = "";
	out_filename_call_indel = "";
	out_filename_call_clipReg = "";
}

//Destructor
Chrome::~Chrome(){
	if(!blockVector.empty()) destroyBlockVector();
	if(!var_cand_vec.empty()) destroyVarCandVector(var_cand_vec);
	if(!mis_aln_vec.empty()) destroyMisAlnVector();
	if(!clipRegVector.empty()) destroyClipRegVector();
	if(!mateClipRegVector.empty()) destroyMateClipRegVector();
	if(!var_cand_clipReg_vec.empty()) destroyVarCandVector(var_cand_clipReg_vec);
}

// initialization
void Chrome::init(){
	blockNum = 0;
	process_block_num = 0;
	print_flag = true;
	decoy_flag = isDecoyChr(chrname);
	alt_flag = isAltChr(chrname);
}

// set the output directory
void Chrome::setOutputDir(string& out_dir_detect_prefix, string& out_dir_cns_prefix, string& out_dir_call_prefix){
	out_dir_detect = out_dir_detect_prefix + "/" + chrname;
	out_dir_cns = out_dir_cns_prefix + "/" + chrname;
	out_dir_call = out_dir_call_prefix + "/" + chrname;

	out_dir_detect = preprocessPipeChar(out_dir_detect);
	out_dir_cns = preprocessPipeChar(out_dir_cns);
	out_dir_call = preprocessPipeChar(out_dir_call);

	out_filename_detect_snv = out_dir_detect + "_SNV_candidate";
	out_filename_detect_indel = out_dir_detect + "_INDEL_candidate";
	out_filename_detect_clipReg = out_dir_detect + "_clipReg_candidate";
	out_filename_call_snv = out_dir_call + "_SNV";
	out_filename_call_indel = out_dir_call + "_INDEL";
	misAln_reg_filename = out_dir_detect + "_misaln_reg";
	out_filename_call_clipReg = out_dir_call + "_clipReg";

	var_cand_indel_filename = out_dir_cns + "_var_cand_indel";
	var_cand_clipReg_filename = out_dir_cns + "_var_cand_clipReg";

	blat_var_cand_indel_filename = out_dir_call_prefix + "/" + "blat_aln_info_" + chrname + "_indel";
	blat_var_cand_clipReg_filename = out_dir_call_prefix + "/" + "blat_aln_info_" + chrname + "_clipReg";
}

string Chrome::getVarcandIndelFilename(){
	return var_cand_indel_filename;
}
string Chrome::getVarcandClipregFilename(){
	return var_cand_clipReg_filename;
}

// generate the chromosome blocks
int Chrome::generateChrBlocks(){
	int64_t pos, begPos, endPos;
	Block *block_tmp;
	bool headIgnFlag, tailIgnFlag, block_process_flag;
	vector<simpleReg_t*> sub_limit_reg_vec;

	process_block_num = 0;
	blockNum = 0;
	pos = 1;
	while(pos<=chrlen){
		begPos = pos;
		endPos = pos + paras->blockSize - 1;

//		if(begPos<14781633 and endPos>14781633){
//			cout<<begPos<<endl;
//		}

		if(chrlen-endPos<=paras->blockSize*0.5){
			endPos = chrlen;
			pos = endPos + 1;
		}else
			pos += paras->blockSize - 2 * paras->slideSize;
		blockNum ++;

		if(begPos==1) headIgnFlag = false;
		else headIgnFlag = true;
		if(endPos==chrlen) tailIgnFlag = false;
		else tailIgnFlag = true;

		// allocate block according to 'process_flag'
		block_process_flag = true;
		if(paras->limit_reg_process_flag){
			//sub_limit_reg_vec = getSimpleRegs(chrname, begPos, endPos, paras->limit_reg_vec);
			sub_limit_reg_vec = getOverlappedSimpleRegs(chrname, begPos, endPos, paras->limit_reg_vec);
			if(sub_limit_reg_vec.size()==0) block_process_flag = false;
		}

		block_tmp = allocateBlock(chrname, begPos, endPos, fai, headIgnFlag, tailIgnFlag, block_process_flag);
		block_tmp->setOutputDir(out_dir_detect, out_dir_cns, out_dir_call);
		if(sub_limit_reg_vec.size()) block_tmp->setLimitRegs(sub_limit_reg_vec);
		blockVector.push_back(block_tmp);

		if(block_process_flag) process_block_num ++;
	}
	blockVector.shrink_to_fit();

	if(process_block_num>0) print_flag = true;
	else  print_flag = false;

	return 0;
}

// allocate the memory for one block
Block *Chrome::allocateBlock(string& chrname, int64_t startPos, int64_t endPos, faidx_t *fai, bool headIgnFlag, bool tailIgnFlag, bool block_process_flag){
	Block *block_tmp = new Block(chrname, startPos, endPos, fai, paras);
	if(!block_tmp){
		cerr << "Chrome: cannot allocate memory" << endl;
		exit(1);
	}
	block_tmp->setRegIngFlag(headIgnFlag, tailIgnFlag);
	block_tmp->setProcessFlag(block_process_flag);
	return block_tmp;
}

// free the memory of blockVector
void Chrome::destroyBlockVector(){
	vector<Block*>::iterator bloc;
	for(bloc=blockVector.begin(); bloc!=blockVector.end(); bloc++)
		delete (*bloc);   // free each block
	vector<Block*>().swap(blockVector);
}

// free the memory of var_cand_vec
void Chrome::destroyVarCandVector(vector<varCand*> &var_cand_vec){
	varCand *var_cand;
	size_t i, j, k;
	for(i=0; i<var_cand_vec.size(); i++){
		var_cand = var_cand_vec.at(i);
		for(j=0; j<var_cand->varVec.size(); j++)  // free varVec
			delete var_cand->varVec.at(j);
		for(j=0; j<var_cand->newVarVec.size(); j++)  // free newVarVec
			delete var_cand->newVarVec.at(j);

		blat_aln_t *item_blat;
		for(j=0; j<var_cand->blat_aln_vec.size(); j++){ // free blat_aln_vec
			item_blat = var_cand->blat_aln_vec.at(j);
			for(k=0; k<item_blat->aln_segs.size(); k++) // free aln_segs
				delete item_blat->aln_segs.at(k);
			delete item_blat;
		}
		delete var_cand;   // free each item
	}
	vector<varCand*>().swap(var_cand_vec);
}

// free the memory of mis_aln_vec
void Chrome::destroyMisAlnVector(){
	vector<reg_t*>::iterator it;
	for(it=mis_aln_vec.begin(); it!=mis_aln_vec.end(); it++)  // free mis_aln_vec
		delete *it;
	vector<reg_t*>().swap(mis_aln_vec);
}

// free the memory of clipRegVector
void Chrome::destroyClipRegVector(){
	vector<reg_t*>::iterator it;
	for(it=clipRegVector.begin(); it!=clipRegVector.end(); it++)
		delete *it;
	vector<reg_t*>().swap(clipRegVector);
}

// free the memory of mateClipRegVector
void Chrome::destroyMateClipRegVector(){
	vector<mateClipReg_t*>::iterator it;
	for(it=mateClipRegVector.begin(); it!=mateClipRegVector.end(); it++){
		if((*it)->leftClipReg) delete (*it)->leftClipReg;
		if((*it)->leftClipReg2) delete (*it)->leftClipReg2;
		if((*it)->rightClipReg) delete (*it)->rightClipReg;
		if((*it)->rightClipReg2) delete (*it)->rightClipReg2;
		if((*it)->largeIndelClipReg) delete (*it)->largeIndelClipReg;
		delete *it;
	}
	vector<mateClipReg_t*>().swap(mateClipRegVector);
}

// save the blocks to file
void Chrome::saveChrBlocksToFile(){
	ofstream outfile(blocks_out_file);
	if(!outfile.is_open()){
		cerr << "Chrome: Cannot open file " << blocks_out_file << endl;
		exit(1);
	}

	Block* bloc;
	for(size_t i=0; i<blockVector.size(); i++){
		bloc = blockVector.at(i);
		outfile << bloc->chrname << "\t" << bloc->startPos << "\t" << bloc->endPos << endl;
	}
	outfile.close();
}

// fill the estimation data
void Chrome::chrFillDataEst(size_t op_est){
	int i, k, pos, begPos, endPos, seq_len;
	Block *block_tmp;
	string reg;
	char *seq;
	bool flag;
	double factor;

	if(chrlen>=MIN_CHR_SIZE_EST){ // minimal chromosome size 120kb
		for(k=0; k<2; k++){
			if(k==0) factor = 0.25;
			else factor = 0.75;

			//pos = chrlen / 2;
			//while(pos<chrlen-MIN_CHR_SIZE_EST/2){
			pos = chrlen * factor;
			while(pos<chrlen-MIN_CHR_SIZE_EST/4){
				flag = true;

				// construct a block with a 20kb size
				begPos = pos - BLOCK_SIZE_EST/2;
				endPos = begPos + BLOCK_SIZE_EST - 1;
				reg = chrname + ":" + to_string(begPos) + "-" + to_string(endPos);
				pthread_mutex_lock(&mutex_fai);
				seq = fai_fetch(fai, reg.c_str(), &seq_len);
				pthread_mutex_unlock(&mutex_fai);
				for(i=0; i<seq_len; i++){
					if(seq[i]=='N' or seq[i]=='n'){
						flag = false;
						break;
					}
				}
				free(seq);

				if(!flag){ // invalid region, slide to next region
					pos = endPos + 1;
				}else break;
			}

			if(flag){ // valid region
				//cout << "Est region: " << chrname << ":" << begPos << "-" << endPos << endl;
				block_tmp = allocateBlock(chrname, begPos, endPos, fai, false, false, true);
				// fill the data
				block_tmp->blockFillDataEst(op_est);
				delete block_tmp;

				paras->reg_sum_size_est += endPos - begPos + 1;
			}
		}
	}
}

// detect indels for chrome
int Chrome::chrDetect(){
	Time time;

	if(print_flag) cout << "[" << time.getTime() << "]: processing Chr: " << chrname << ", size: " << chrlen << " bp" << endl;

	//mkdir(out_dir_detect.c_str(), S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);  // create the directory for detect command // delete on 2024-04-06

	chrSetMisAlnRegFile();

	if(paras->num_threads<=1) chrDetect_st();  // single thread
	else chrDetect_mt();  // multiple threads

	// detect mated clip regions
	//cout << "[" << time.getTime() << "]: compute mate clip region on chromosome " << chrname << " ..." << endl;
	chrComputeMateClipReg();

	// remove FP indels and Snvs in clipping regions
	//cout << "[" << time.getTime() << "]: remove FP indels and SNVs in mate clip region on chromosome " << chrname << " ..." << endl;
	removeFPIndelSnvInClipReg(mateClipRegVector);

	// remove redundant Indels for 'detect' command
	//cout << "[" << time.getTime() << "]: remove redundant indels on chromosome " << chrname << " ..." << endl;
	removeRedundantIndelDetect();

	// merge the results to single file
	//chrMergeDetectResultToFile();

	chrResetMisAlnRegFile();

	return 0;
}

// single thread
int Chrome::chrDetect_st(){
	Time time;
	Block* bloc;
	for(size_t i=0; i<blockVector.size(); i++){
		bloc = blockVector.at(i);
		if(bloc->process_flag)
		{
			//cout << "[" << time.getTime() << "]: [" << i << "]: detect files:" << bloc->indelFilenameDetect << endl;
			bloc->blockDetect();
		}
	}

	return 0;
}

// multiple threads
int Chrome::chrDetect_mt(){
	int32_t i;
	MultiThread mt[paras->num_threads];
	for(i=0; i<paras->num_threads; i++){
		mt[i].setNumThreads(paras->num_threads);
		mt[i].setBlockVec(&blockVector);
		mt[i].setUserThreadID(i);
		if(!mt[i].startDetect()){
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

	return 0;
}

// remove repeatedly detected indels
void Chrome::removeRedundantIndelDetect(){
	size_t i, j;
	Block *bloc;
	reg_t *reg;

	for(i=0; i<blockVector.size(); i++){
		bloc = blockVector.at(i);
		for(j=0; j<bloc->indelVector.size(); j++){
			reg = bloc->indelVector.at(j);
			removeRedundantIndelItemDetect(reg, i, j);
		}
	}
}

void Chrome::removeRedundantIndelItemDetect(reg_t *reg, int32_t bloc_idx, int32_t indel_vec_idx){
	int64_t i, j, start_indel_vec_idx;
	Block *bloc;
	reg_t *reg_tmp;
	bool flag;
	int64_t startPos1, endPos1, startPos2, endPos2;

	// find the left side
	for(i=bloc_idx; i>=0; i--){
		bloc = blockVector.at(i);
		startPos1 = reg->startRefPos - MAX_VAR_REG_SIZE;
		endPos1 = reg->endRefPos + MAX_VAR_REG_SIZE;
		startPos2 = bloc->startPos - MAX_VAR_REG_SIZE;
		endPos2 = bloc->endPos + MAX_VAR_REG_SIZE;
		if(startPos1<1) startPos1 = 1;
		if(startPos2<1) startPos2 = 1;
		flag = isOverlappedPos(startPos1, endPos1, startPos2, endPos2);
		if(flag){ // overlapped
			if(i==bloc_idx) start_indel_vec_idx = indel_vec_idx + 1;
			else start_indel_vec_idx = 0;
			for(j=start_indel_vec_idx; j<(int64_t)bloc->indelVector.size(); ){
				reg_tmp = bloc->indelVector.at(j);
				if(reg_tmp->chrname.compare(reg->chrname)==0 and reg_tmp->startRefPos==reg->startRefPos and reg_tmp->endRefPos==reg->endRefPos){
					// invalid item, and delete it
					delete reg_tmp;
					bloc->indelVector.erase(bloc->indelVector.begin()+j);
				}else j++;
			}
		}else break;
	}

	// find the right side
	for(i=bloc_idx+1; i<(int64_t)blockVector.size(); i++){
		bloc = blockVector.at(i);
		if(bloc->indelVector.size()>0){
			startPos1 = reg->startRefPos - MAX_VAR_REG_SIZE;
			endPos1 = reg->endRefPos + MAX_VAR_REG_SIZE;
			startPos2 = bloc->startPos - MAX_VAR_REG_SIZE;
			endPos2 = bloc->endPos + MAX_VAR_REG_SIZE;
			if(startPos1<1) startPos1 = 1;
			if(startPos2<1) startPos2 = 1;
			flag = isOverlappedPos(startPos1, endPos1, startPos2, endPos2);
			if(flag){ // overlapped
				for(j=0; j<(int64_t)bloc->indelVector.size(); ){
					reg_tmp = bloc->indelVector.at(j);
					if(reg_tmp->chrname.compare(reg->chrname)==0 and reg_tmp->startRefPos==reg->startRefPos and reg_tmp->endRefPos==reg->endRefPos){
						// invalid item, and delete it
						delete reg_tmp;
						bloc->indelVector.erase(bloc->indelVector.begin()+j);
					}else j++;
				}
			}else break;
		}
	}
}

// detect mated clip regions for duplications and inversions
void Chrome::chrComputeMateClipReg(){
	//Time time;
	size_t i;
	Block *bloc;
	//reg_t *reg;
	//vector<bool> clip_processed_flag_vec;

	for(i=0; i<blockVector.size(); i++){
		bloc = blockVector.at(i);
		for(size_t j=0; j<bloc->clipRegVector.size(); j++) clipRegVector.push_back(bloc->clipRegVector.at(j));
	}
	clipRegVector.shrink_to_fit();

	sortRegVec(clipRegVector);
	if(isRegSorted(clipRegVector)==false){
		cout << "line=" << __LINE__ << ", region is not sorted in order." << endl;
		exit(1);
	}

//	for(i=0; i<clipRegVector.size(); i++) clip_processed_flag_vec.push_back(false);
//	clip_processed_flag_vec.shrink_to_fit();

	// remove redundant items
	//removeRedundantItemsClipReg(clipRegVector, clip_processed_flag_vec);
	removeRedundantItemsClipReg(clipRegVector);

	if(isRegSorted(clipRegVector)==false){
		cout << "line=" << __LINE__ << ", region is not sorted in order." << endl;
		exit(1);
	}

	//processMateClipRegDetectWork(clip_processed_flag_vec);
	processMateClipRegDetectWork();

	// sort
	sortMateClipRegsByWorkId();
	bubbleSortMateClipRegs();

	// check sort
	checkSortMateClipRegs();

	// merge same TRA region
	//mergeSameRegTRA();

	// remove FPs by detecting false overlapped clip regions for DUPs and INVs
	removeFPClipRegsDupInv();
}

// detect mate clipping regions using thread pool
//void Chrome::processMateClipRegDetectWork(vector<bool> &clip_processed_flag_vec){
void Chrome::processMateClipRegDetectWork(){
	int32_t num_work;
	reg_t *reg;
	mateClipRegDetectWork_opt *mate_clip_reg_work_opt;
	//Time time;

	if(clipRegVector.empty()) return;  // no work, then return directly

	//cout << "[" << time.getTime() << "], Chr " << chrname << ", detect mate clipping regions using " << paras->num_threads << " threads ..." << endl;

	hts_tpool *p = hts_tpool_init(paras->num_threads);
	hts_tpool_process *q = hts_tpool_process_init(p, paras->num_threads*2, 1);

	pthread_mutex_init(&mutex_mate_clip_reg, NULL);

	num_work = clipRegVector.size();
	for(int32_t i=0; i<num_work; i++){
		//if(clip_processed_flag_vec.at(i)==false){
			reg = clipRegVector.at(i);

			// cout << "[" << i << "]: " << reg->chrname << ":" << reg->startRefPos << "-" << reg->endRefPos << endl;

			// if(i>=36 and i<=36) // 36, 60
			// 	cout << reg->chrname << ":" << reg->startRefPos << "-" << reg->endRefPos << endl;
			// else continue;

			mate_clip_reg_work_opt = new mateClipRegDetectWork_opt();
			mate_clip_reg_work_opt->reg = reg;
			mate_clip_reg_work_opt->work_id = i;
			mate_clip_reg_work_opt->num_work = num_work;
			mate_clip_reg_work_opt->fai = fai;
			mate_clip_reg_work_opt->paras = paras;
			//mate_clip_reg_work_opt->clip_processed_flag_vec = &clip_processed_flag_vec;
			mate_clip_reg_work_opt->mateClipRegVector = &mateClipRegVector;
			mate_clip_reg_work_opt->clipRegVector = &clipRegVector;
			mate_clip_reg_work_opt->var_cand_clipReg_vec = &var_cand_clipReg_vec;
			mate_clip_reg_work_opt->blockVector = &blockVector;
			mate_clip_reg_work_opt->p_mutex_mate_clip_reg = &mutex_mate_clip_reg;
			//mate_clip_reg_work_opt->p_mate_clip_reg_fail_num = &mate_clip_reg_fail_num;

			hts_tpool_dispatch(p, q, processSingleMateClipRegDetectWork, mate_clip_reg_work_opt);
		//}
	}

    hts_tpool_process_flush(q);
    hts_tpool_process_destroy(q);
    hts_tpool_destroy(p);

    //cout << "\tChr " << chrname << ", num_work=" << num_work << ", mate_clip_reg_fail_num=" << mate_clip_reg_fail_num << endl;
}

// sort mate clipping region items by work id
void Chrome::sortMateClipRegsByWorkId(){
	size_t i;
	mateClipReg_t *mate_clip_1, *mate_clip_2, *tmp;

	if(mateClipRegVector.size()<2) return;

	for(i=0; i<mateClipRegVector.size()-1; i++){
		mate_clip_1 = mateClipRegVector.at(i);
		mate_clip_2 = mateClipRegVector.at(i+1);
		if(mate_clip_1->work_id>mate_clip_2->work_id){
			tmp = mateClipRegVector.at(i);
			mateClipRegVector.at(i) = mateClipRegVector.at(i+1);
			mateClipRegVector.at(i+1) = tmp;
		}
	}
}

// sort mate clipping region items by bubble sort
void Chrome::bubbleSortMateClipRegs(){
	size_t i, j;
	mateClipReg_t *mate_clip_1, *mate_clip_2, *tmp;
	int64_t left_pos1, left_pos2, other_left_pos1, other_left_pos2, left_pos, other_left_pos;
	int64_t right_pos1, right_pos2, other_right_pos1, other_right_pos2, right_pos, other_right_pos;
	int64_t left_end_pos1, left_end_pos2, other_left_end_pos1, other_left_end_pos2, left_end_pos, other_left_end_pos;
	int64_t right_end_pos1, right_end_pos2, other_right_end_pos1, other_right_end_pos2, right_end_pos, other_right_end_pos;
	int64_t clipPosNum, other_clipPosNum;
	int64_t min_start_pos, min_end_pos, min_other_start_pos, min_other_end_pos;
	string chrname_tmp;
	bool swap_flag, equal_flag;

	if(mateClipRegVector.size()<2) return;

	for(i=0; i<mateClipRegVector.size(); i++){
		for(j=0; j<mateClipRegVector.size()-i-1; j++){
			mate_clip_1 = mateClipRegVector.at(j);
			mate_clip_2 = mateClipRegVector.at(j+1);

			if(mate_clip_1->large_indel_flag==false and mate_clip_2->large_indel_flag==false){ // both clipping regions
				left_pos1 = left_pos2 = other_left_pos1 = other_left_pos2 = left_pos = other_left_pos = -1;
				right_pos1 = right_pos2 = other_right_pos1 = other_right_pos2 = right_pos = other_right_pos = -1;

				if(mate_clip_1->leftClipReg) {
					left_pos1 = mate_clip_1->leftClipReg->startRefPos;
					left_end_pos1 = mate_clip_1->leftClipReg->endRefPos;
				}
				if(mate_clip_1->leftClipReg2) {
					left_pos2 = mate_clip_1->leftClipReg2->startRefPos;
					left_end_pos2 = mate_clip_1->leftClipReg2->endRefPos;
				}
				if(mate_clip_1->rightClipReg){
					right_pos1 = mate_clip_1->rightClipReg->startRefPos;
					right_end_pos1 = mate_clip_1->rightClipReg->endRefPos;
				}
				if(mate_clip_1->rightClipReg2){
					right_pos2 = mate_clip_1->rightClipReg2->startRefPos;
					right_end_pos2 = mate_clip_1->rightClipReg2->endRefPos;
				}

				if(mate_clip_2->leftClipReg){
					other_left_pos1 = mate_clip_2->leftClipReg->startRefPos;
					other_left_end_pos1 = mate_clip_2->leftClipReg->endRefPos;
				}
				if(mate_clip_2->leftClipReg2){
					other_left_pos2 = mate_clip_2->leftClipReg2->startRefPos;
					other_left_end_pos2 = mate_clip_2->leftClipReg2->endRefPos;
				}
				if(mate_clip_2->rightClipReg){
					other_right_pos1 = mate_clip_2->rightClipReg->startRefPos;
					other_right_end_pos1 = mate_clip_2->rightClipReg->endRefPos;
				}
				if(mate_clip_2->rightClipReg2){
					other_right_pos2 = mate_clip_2->rightClipReg2->startRefPos;
					other_right_end_pos2 = mate_clip_2->rightClipReg2->endRefPos;
				}

				// ignore the second region if there are two regions
				if(left_pos1!=-1) { left_pos = left_pos1; left_end_pos = left_end_pos1; }
				else{ left_pos = left_pos2; left_end_pos = left_end_pos2; }
				if(right_pos1!=-1){ right_pos = right_pos1; right_end_pos = right_end_pos1; }
				else{ right_pos = right_pos2; right_end_pos = right_end_pos2; }

				if(other_left_pos1!=-1){ other_left_pos = other_left_pos1; other_left_end_pos = other_left_end_pos1; }
				else{ other_left_pos = other_left_pos2; other_left_end_pos = other_left_end_pos2; }
				if(other_right_pos1!=-1){ other_right_pos = other_right_pos1; other_right_end_pos = other_right_end_pos1; }
				else{ other_right_pos = other_right_pos2; other_right_end_pos = other_right_end_pos2; }

				clipPosNum = mate_clip_1->leftClipPosNum + mate_clip_1->leftClipPosNum2 + mate_clip_1->rightClipPosNum + mate_clip_1->rightClipPosNum2;
				other_clipPosNum = mate_clip_2->leftClipPosNum + mate_clip_2->leftClipPosNum2 + mate_clip_2->rightClipPosNum + mate_clip_2->rightClipPosNum2;

				swap_flag = equal_flag = false;
				if(left_pos>other_left_pos) swap_flag = true;
				else if(left_pos==other_left_pos){ // equal left loci, then compare right loci
					if(right_pos>other_right_pos) swap_flag = true;
					else if(right_pos==other_right_pos){ // equal right loci, then compare left end loci
						if(left_end_pos>other_left_end_pos) swap_flag = true;
						else if(left_end_pos==other_left_end_pos){ // equal left end loci, then compare right end loci
							if(right_end_pos>other_right_end_pos) swap_flag = true;
							else if(right_end_pos==other_right_end_pos){ // check clipPosNum
								if(clipPosNum<other_clipPosNum) swap_flag = true;
								else if(clipPosNum==other_clipPosNum) equal_flag = true;
							}
						}
					}
				}

				// sort by the second region
				if(equal_flag and mate_clip_1->leftClipRegNum==2 and mate_clip_2->leftClipRegNum==2){
					if(left_pos2>other_left_pos2) swap_flag = true;
					else if(left_pos2==other_left_pos2){
						if(left_end_pos2>other_left_end_pos2) swap_flag = true;
						else if(left_end_pos2==other_left_end_pos2) equal_flag = true;
					}

				}
				if(equal_flag and mate_clip_1->rightClipRegNum==2 and mate_clip_2->rightClipRegNum==2){
					if(right_pos2>other_right_pos2) swap_flag = true;
					else if(right_pos2==other_right_pos2){
						if(right_end_pos2>other_right_end_pos2) swap_flag = true;
						else if(right_end_pos2==other_right_end_pos2) equal_flag = true;
					}
				}

				if(swap_flag){ // swap
					tmp = mateClipRegVector.at(j);
					mateClipRegVector.at(j) = mateClipRegVector.at(j+1);
					mateClipRegVector.at(j+1) = tmp;
				}
			}else{ // at least one large indel region
				min_start_pos = min_end_pos = LONG_MAX;
				if(mate_clip_1->large_indel_flag==false){ // clipping region
					chrname_tmp = mate_clip_2->largeIndelClipReg->chrname;
					if(mate_clip_1->leftClipReg and mate_clip_1->leftClipReg->chrname.compare(chrname_tmp)==0) {
						left_pos1 = mate_clip_1->leftClipReg->startRefPos;
						left_end_pos1 = mate_clip_1->leftClipReg->endRefPos;
						if(left_pos1<min_start_pos) min_start_pos = left_pos1;
						if(left_end_pos1<min_end_pos) min_end_pos = left_end_pos1;
					}
					if(mate_clip_1->leftClipReg2 and mate_clip_1->leftClipReg2->chrname.compare(chrname_tmp)==0) {
						left_pos2 = mate_clip_1->leftClipReg2->startRefPos;
						left_end_pos2 = mate_clip_1->leftClipReg2->endRefPos;
						if(left_pos2<min_start_pos) min_start_pos = left_pos2;
						if(left_end_pos2<min_end_pos) min_end_pos = left_end_pos2;
					}
					if(mate_clip_1->rightClipReg and mate_clip_1->rightClipReg->chrname.compare(chrname_tmp)==0){
						right_pos1 = mate_clip_1->rightClipReg->startRefPos;
						right_end_pos1 = mate_clip_1->rightClipReg->endRefPos;
						if(right_pos1<min_start_pos) min_start_pos = right_pos1;
						if(right_end_pos1<min_end_pos) min_end_pos = right_end_pos1;
					}
					if(mate_clip_1->rightClipReg2 and mate_clip_1->rightClipReg2->chrname.compare(chrname_tmp)==0){
						right_pos2 = mate_clip_1->rightClipReg2->startRefPos;
						right_end_pos2 = mate_clip_1->rightClipReg2->endRefPos;
						if(right_pos2<min_start_pos) min_start_pos = right_pos2;
						if(right_end_pos2<min_end_pos) min_end_pos = right_end_pos2;
					}
				}else{ // large indel
					min_start_pos = mate_clip_1->largeIndelClipReg->startRefPos;
					min_end_pos = mate_clip_1->largeIndelClipReg->endRefPos;
				}

				min_other_start_pos = min_other_end_pos = LONG_MAX;
				if(mate_clip_2->large_indel_flag==false){ // clipping region
					chrname_tmp = mate_clip_1->largeIndelClipReg->chrname;
					if(mate_clip_2->leftClipReg and mate_clip_2->leftClipReg->chrname.compare(chrname_tmp)==0) {
						other_left_pos1 = mate_clip_2->leftClipReg->startRefPos;
						other_left_end_pos1 = mate_clip_2->leftClipReg->endRefPos;
						if(other_left_pos1<min_other_start_pos) min_other_start_pos = other_left_pos1;
						if(other_left_end_pos1<min_other_end_pos) min_other_end_pos = other_left_end_pos1;
					}
					if(mate_clip_2->leftClipReg2 and mate_clip_2->leftClipReg2->chrname.compare(chrname_tmp)==0) {
						other_left_pos2 = mate_clip_2->leftClipReg2->startRefPos;
						other_left_end_pos2 = mate_clip_2->leftClipReg2->endRefPos;
						if(other_left_pos2<min_start_pos) min_other_start_pos = other_left_pos2;
						if(other_left_end_pos2<min_end_pos) min_other_end_pos = other_left_end_pos2;
					}
					if(mate_clip_2->rightClipReg and mate_clip_2->rightClipReg->chrname.compare(chrname_tmp)==0){
						other_right_pos1 = mate_clip_2->rightClipReg->startRefPos;
						other_right_end_pos1 = mate_clip_2->rightClipReg->endRefPos;
						if(other_right_pos1<min_other_start_pos) min_other_start_pos = other_right_pos1;
						if(other_right_end_pos1<min_other_end_pos) min_other_end_pos = other_right_end_pos1;
					}
					if(mate_clip_2->rightClipReg2 and mate_clip_2->rightClipReg2->chrname.compare(chrname_tmp)==0){
						other_right_pos2 = mate_clip_2->rightClipReg2->startRefPos;
						other_right_end_pos2 = mate_clip_2->rightClipReg2->endRefPos;
						if(other_right_pos2<min_other_start_pos) min_other_start_pos = other_right_pos2;
						if(other_right_end_pos2<min_other_end_pos) min_other_end_pos = other_right_end_pos2;
					}
				}else{ // large indel
					min_other_start_pos = mate_clip_2->largeIndelClipReg->startRefPos;
					min_other_end_pos = mate_clip_2->largeIndelClipReg->endRefPos;
				}

				swap_flag = false;
				if(min_start_pos>min_other_start_pos) swap_flag = true;
				else if(min_start_pos==min_other_start_pos){
					if(min_end_pos>min_other_end_pos) swap_flag = true;
				}

				if(swap_flag){ // swap
					tmp = mateClipRegVector.at(j);
					mateClipRegVector.at(j) = mateClipRegVector.at(j+1);
					mateClipRegVector.at(j+1) = tmp;
				}
			}
		}
	}
}

// check sort mate clipping region items
void Chrome::checkSortMateClipRegs(){
	size_t i;
	mateClipReg_t *mate_clip_1, *mate_clip_2;
	int64_t left_pos1, left_pos2, other_left_pos1, other_left_pos2, left_pos, other_left_pos, err_pos1, err_pos2;
	int64_t right_pos1, right_pos2, other_right_pos1, other_right_pos2, right_pos, other_right_pos;
	int64_t left_end_pos1, left_end_pos2, other_left_end_pos1, other_left_end_pos2, left_end_pos, other_left_end_pos;
	int64_t right_end_pos1, right_end_pos2, other_right_end_pos1, other_right_end_pos2, right_end_pos, other_right_end_pos;
	int64_t clipPosNum, other_clipPosNum;
	int64_t min_start_pos, min_end_pos, min_other_start_pos, min_other_end_pos;
	string chrname_tmp;
	bool swap_flag, equal_flag;

	if(mateClipRegVector.size()<2) return;

	for(i=0; i<mateClipRegVector.size()-1; i++){
		mate_clip_1 = mateClipRegVector.at(i);
		mate_clip_2 = mateClipRegVector.at(i+1);

		if(mate_clip_1->large_indel_flag==false and mate_clip_2->large_indel_flag==false){ // both clipping regions
			left_pos1 = left_pos2 = other_left_pos1 = other_left_pos2 = left_pos = other_left_pos = -1;
			right_pos1 = right_pos2 = other_right_pos1 = other_right_pos2 = right_pos = other_right_pos = -1;

			if(mate_clip_1->leftClipReg) {
				left_pos1 = mate_clip_1->leftClipReg->startRefPos;
				left_end_pos1 = mate_clip_1->leftClipReg->endRefPos;
			}
			if(mate_clip_1->leftClipReg2) {
				left_pos2 = mate_clip_1->leftClipReg2->startRefPos;
				left_end_pos2 = mate_clip_1->leftClipReg2->endRefPos;
			}
			if(mate_clip_1->rightClipReg){
				right_pos1 = mate_clip_1->rightClipReg->startRefPos;
				right_end_pos1 = mate_clip_1->rightClipReg->endRefPos;
			}
			if(mate_clip_1->rightClipReg2){
				right_pos2 = mate_clip_1->rightClipReg2->startRefPos;
				right_end_pos2 = mate_clip_1->rightClipReg2->endRefPos;
			}

			if(mate_clip_2->leftClipReg){
				other_left_pos1 = mate_clip_2->leftClipReg->startRefPos;
				other_left_end_pos1 = mate_clip_2->leftClipReg->endRefPos;
			}
			if(mate_clip_2->leftClipReg2){
				other_left_pos2 = mate_clip_2->leftClipReg2->startRefPos;
				other_left_end_pos2 = mate_clip_2->leftClipReg2->endRefPos;
			}
			if(mate_clip_2->rightClipReg){
				other_right_pos1 = mate_clip_2->rightClipReg->startRefPos;
				other_right_end_pos1 = mate_clip_2->rightClipReg->endRefPos;
			}
			if(mate_clip_2->rightClipReg2){
				other_right_pos2 = mate_clip_2->rightClipReg2->startRefPos;
				other_right_end_pos2 = mate_clip_2->rightClipReg2->endRefPos;
			}

			// ignore the second region if there are two regions
			if(left_pos1!=-1) { left_pos = left_pos1; left_end_pos = left_end_pos1; }
			else{ left_pos = left_pos2; left_end_pos = left_end_pos2; }
			if(right_pos1!=-1){ right_pos = right_pos1; right_end_pos = right_end_pos1; }
			else{ right_pos = right_pos2; right_end_pos = right_end_pos2; }

			if(other_left_pos1!=-1){ other_left_pos = other_left_pos1; other_left_end_pos = other_left_end_pos1; }
			else{ other_left_pos = other_left_pos2; other_left_end_pos = other_left_end_pos2; }
			if(other_right_pos1!=-1){ other_right_pos = other_right_pos1; other_right_end_pos = other_right_end_pos1; }
			else{ other_right_pos = other_right_pos2; other_right_end_pos = other_right_end_pos2; }

			clipPosNum = mate_clip_1->leftClipPosNum + mate_clip_1->leftClipPosNum2 + mate_clip_1->rightClipPosNum + mate_clip_1->rightClipPosNum2;
			other_clipPosNum = mate_clip_2->leftClipPosNum + mate_clip_2->leftClipPosNum2 + mate_clip_2->rightClipPosNum + mate_clip_2->rightClipPosNum2;

			swap_flag = equal_flag = false;
			err_pos1 = err_pos2 = -1;
			if(left_pos>other_left_pos){
				swap_flag = true;
				err_pos1 = left_pos;
				err_pos2 = other_left_pos;
			}else if(left_pos==other_left_pos){ // equal left loci, then compare right loci
				if(right_pos>other_right_pos) {
					swap_flag = true;
					err_pos1 = right_pos;
					err_pos2 = other_right_pos;
				}else if(right_pos==other_right_pos){ // equal right loci, then compare left end loci
					if(left_end_pos>other_left_end_pos){
						swap_flag = true;
						err_pos1 = left_end_pos;
						err_pos2 = other_left_end_pos;
					}else if(left_end_pos==other_left_end_pos){ // equal left end loci, then compare right end loci
						if(right_end_pos>other_right_end_pos){
							swap_flag = true;
							err_pos1 = right_end_pos;
							err_pos2 = other_right_end_pos;
						}else if(right_end_pos==other_right_end_pos){ // check clipPosNum
							if(clipPosNum<other_clipPosNum) {
								swap_flag = true;
								err_pos1 = right_end_pos;
								err_pos2 = other_right_end_pos;
							}else if(clipPosNum==other_clipPosNum) equal_flag = true;
						}
					}
				}
			}

			// check the second region
			if(equal_flag and mate_clip_1->leftClipRegNum==2 and mate_clip_2->leftClipRegNum==2){
				if(left_pos2>other_left_pos2){
					swap_flag = true;
					err_pos1 = left_pos2;
					err_pos2 = other_left_pos2;
				}else if(left_pos2==other_left_pos2) {
					if(left_end_pos2>other_left_end_pos2){
						swap_flag = true;
						err_pos1 = left_end_pos2;
						err_pos2 = other_left_end_pos2;
					}else if(left_end_pos2==other_left_end_pos2) equal_flag = true;
				}
			}
			if(equal_flag and mate_clip_1->rightClipRegNum==2 and mate_clip_2->rightClipRegNum==2){
				if(right_pos2>other_right_pos2){
					swap_flag = true;
					err_pos1 = right_pos2;
					err_pos2 = other_right_pos2;
				}else if(right_pos2==other_right_pos2){
					if(right_end_pos2>other_right_end_pos2){
						swap_flag = true;
						err_pos1 = right_end_pos2;
						err_pos2 = other_right_end_pos2;
					}else if(right_end_pos2==other_right_end_pos2) equal_flag = true;
				}
			}

			if(swap_flag){
				cerr << "line=" << __LINE__ << ", " << __func__ << ": i=" << i << ", " << err_pos1 << " > " << err_pos2 << ", order error!" << endl;
				exit(1);
			}
		}else{ // at least one large indel region
			min_start_pos = min_end_pos = LONG_MAX;
			if(mate_clip_1->large_indel_flag==false){ // clipping region
				chrname_tmp = mate_clip_2->largeIndelClipReg->chrname;
				if(mate_clip_1->leftClipReg and mate_clip_1->leftClipReg->chrname.compare(chrname_tmp)==0) {
					left_pos1 = mate_clip_1->leftClipReg->startRefPos;
					left_end_pos1 = mate_clip_1->leftClipReg->endRefPos;
					if(left_pos1<min_start_pos) min_start_pos = left_pos1;
					if(left_end_pos1<min_end_pos) min_end_pos = left_end_pos1;
				}
				if(mate_clip_1->leftClipReg2 and mate_clip_1->leftClipReg2->chrname.compare(chrname_tmp)==0) {
					left_pos2 = mate_clip_1->leftClipReg2->startRefPos;
					left_end_pos2 = mate_clip_1->leftClipReg2->endRefPos;
					if(left_pos2<min_start_pos) min_start_pos = left_pos2;
					if(left_end_pos2<min_end_pos) min_end_pos = left_end_pos2;
				}
				if(mate_clip_1->rightClipReg and mate_clip_1->rightClipReg->chrname.compare(chrname_tmp)==0){
					right_pos1 = mate_clip_1->rightClipReg->startRefPos;
					right_end_pos1 = mate_clip_1->rightClipReg->endRefPos;
					if(right_pos1<min_start_pos) min_start_pos = right_pos1;
					if(right_end_pos1<min_end_pos) min_end_pos = right_end_pos1;
				}
				if(mate_clip_1->rightClipReg2 and mate_clip_1->rightClipReg2->chrname.compare(chrname_tmp)==0){
					right_pos2 = mate_clip_1->rightClipReg2->startRefPos;
					right_end_pos2 = mate_clip_1->rightClipReg2->endRefPos;
					if(right_pos2<min_start_pos) min_start_pos = right_pos2;
					if(right_end_pos2<min_end_pos) min_end_pos = right_end_pos2;
				}
			}else{ // large indel
				min_start_pos = mate_clip_1->largeIndelClipReg->startRefPos;
				min_end_pos = mate_clip_1->largeIndelClipReg->endRefPos;
			}

			min_other_start_pos = min_other_end_pos = LONG_MAX;
			if(mate_clip_2->large_indel_flag==false){ // clipping region
				chrname_tmp = mate_clip_1->largeIndelClipReg->chrname;
				if(mate_clip_2->leftClipReg and mate_clip_2->leftClipReg->chrname.compare(chrname_tmp)==0) {
					other_left_pos1 = mate_clip_2->leftClipReg->startRefPos;
					other_left_end_pos1 = mate_clip_2->leftClipReg->endRefPos;
					if(other_left_pos1<min_other_start_pos) min_other_start_pos = other_left_pos1;
					if(other_left_end_pos1<min_other_end_pos) min_other_end_pos = other_left_end_pos1;
				}
				if(mate_clip_2->leftClipReg2 and mate_clip_2->leftClipReg2->chrname.compare(chrname_tmp)==0) {
					other_left_pos2 = mate_clip_2->leftClipReg2->startRefPos;
					other_left_end_pos2 = mate_clip_2->leftClipReg2->endRefPos;
					if(other_left_pos2<min_start_pos) min_other_start_pos = other_left_pos2;
					if(other_left_end_pos2<min_end_pos) min_other_end_pos = other_left_end_pos2;
				}
				if(mate_clip_2->rightClipReg and mate_clip_2->rightClipReg->chrname.compare(chrname_tmp)==0){
					other_right_pos1 = mate_clip_2->rightClipReg->startRefPos;
					other_right_end_pos1 = mate_clip_2->rightClipReg->endRefPos;
					if(other_right_pos1<min_other_start_pos) min_other_start_pos = other_right_pos1;
					if(other_right_end_pos1<min_other_end_pos) min_other_end_pos = other_right_end_pos1;
				}
				if(mate_clip_2->rightClipReg2 and mate_clip_2->rightClipReg2->chrname.compare(chrname_tmp)==0){
					other_right_pos2 = mate_clip_2->rightClipReg2->startRefPos;
					other_right_end_pos2 = mate_clip_2->rightClipReg2->endRefPos;
					if(other_right_pos2<min_other_start_pos) min_other_start_pos = other_right_pos2;
					if(other_right_end_pos2<min_other_end_pos) min_other_end_pos = other_right_end_pos2;
				}
			}else{ // large indel
				min_other_start_pos = mate_clip_2->largeIndelClipReg->startRefPos;
				min_other_end_pos = mate_clip_2->largeIndelClipReg->endRefPos;
			}

			swap_flag = false;
			if(min_start_pos>min_other_start_pos) {
				swap_flag = true;
				err_pos1 = min_start_pos;
				err_pos2 = min_other_start_pos;
			}else if(min_start_pos==min_other_start_pos){
				if(min_end_pos>min_other_end_pos){
					swap_flag = true;
					err_pos1 = min_end_pos;
					err_pos2 = min_other_end_pos;
				}
			}

			if(swap_flag){ // swap
				cerr << "line=" << __LINE__ << ", " << __func__ << ": i=" << i << ", " << err_pos1 << " > " << err_pos2 << ", order error!" << endl;
				exit(1);
			}
		}
	}

	//cout << "mate clipping regions are sorted in order." << endl;
}

//void Chrome::removeRedundantItemsClipReg(vector<reg_t*> &clipReg_vec, vector<bool> &processed_flag_vec){
void Chrome::removeRedundantItemsClipReg(vector<reg_t*> &clipReg_vec){
	size_t i, j;
	reg_t *reg, *reg_tmp;

	for(i=0; i<clipReg_vec.size(); i++){
		reg = clipReg_vec.at(i);
		for(j=i+1; j<clipReg_vec.size(); ){
			reg_tmp = clipReg_vec.at(j);
			if(reg_tmp->chrname.compare(reg->chrname)==0 and reg_tmp->startRefPos==reg->startRefPos and reg_tmp->endRefPos==reg->endRefPos){ // redundant item
				//cout << "redundant: " << reg_tmp->chrname << ":" << reg_tmp->startRefPos << "-" << reg_tmp->endRefPos << endl;
				delete reg_tmp;
				clipReg_vec.erase(clipReg_vec.begin()+j);
				//processed_flag_vec.erase(processed_flag_vec.begin()+j);
			}else j++;
		}
	}
}

// merge same TRA region
void Chrome::mergeSameRegTRA(){
	size_t i, end_flag1, pos_num1, pos_num2, mean_pos1, mean_pos2, tmp;
	mateClipReg_t *clip_reg, *clip_reg_mate_end, *clip_reg_retained, *clip_reg_redudant, *clip_reg_mate_end_tmp;
	reg_t *reg1, *reg2, *reg_tmp;

	for(i=0; i<mateClipRegVector.size(); i++){
		clip_reg = mateClipRegVector.at(i);
		if(clip_reg->valid_flag and clip_reg->reg_mated_flag and clip_reg->sv_type==VAR_TRA and (clip_reg->leftClipRegNum==2 or clip_reg->rightClipRegNum==2)){
			// get other end region according to the same clipping region
			while(1){
				clip_reg_mate_end = getMateRegEndSameClipReg(clip_reg, mateClipRegVector);
				if(clip_reg_mate_end and clip_reg_mate_end!=clip_reg_mate_end_tmp){
					//cout << "line=" << __LINE__ << ", i=" << i << endl;

					clip_reg_mate_end_tmp = clip_reg_mate_end;

					// merge nodes
					end_flag1 = 0;
					reg1 = NULL; pos_num1 = mean_pos1 = 0;
					if(clip_reg->leftClipRegNum==1){
						if(clip_reg->leftClipReg){
							reg1 = clip_reg->leftClipReg;
							pos_num1 = clip_reg->leftClipPosNum;
							mean_pos1 = clip_reg->leftMeanClipPos;
						}else{
							reg1 = clip_reg->leftClipReg2;
							pos_num1 = clip_reg->leftClipPosNum2;
							mean_pos1 = clip_reg->leftMeanClipPos2;
						}
						end_flag1 = 1;
					}else if(clip_reg->rightClipRegNum==1){
						if(clip_reg->rightClipReg){
							reg1 = clip_reg->rightClipReg;
							pos_num1 = clip_reg->rightClipPosNum;
							mean_pos1 = clip_reg->rightMeanClipPos;
						}else{
							reg1 = clip_reg->rightClipReg2;
							pos_num1 = clip_reg->rightClipPosNum2;
							mean_pos1 = clip_reg->rightMeanClipPos2;
						}
						end_flag1 = 2;
					}

					reg2 = NULL; pos_num2 = mean_pos2 = 0;
					if(clip_reg_mate_end->leftClipRegNum==1){
						if(clip_reg_mate_end->leftClipReg){
							reg2 = clip_reg_mate_end->leftClipReg;
							pos_num2 = clip_reg_mate_end->leftClipPosNum;
							mean_pos2 = clip_reg_mate_end->leftMeanClipPos;
						}else{
							reg2 = clip_reg_mate_end->leftClipReg2;
							pos_num2 = clip_reg_mate_end->leftClipPosNum2;
							mean_pos2 = clip_reg_mate_end->leftMeanClipPos2;
						}
					}else if(clip_reg_mate_end->rightClipRegNum==1){
						if(clip_reg_mate_end->rightClipReg){
							reg2 = clip_reg_mate_end->rightClipReg;
							pos_num2 = clip_reg_mate_end->rightClipPosNum;
							mean_pos2 = clip_reg_mate_end->rightMeanClipPos;
						}else{
							reg2 = clip_reg_mate_end->rightClipReg2;
							pos_num2 = clip_reg_mate_end->rightClipPosNum2;
							mean_pos2 = clip_reg_mate_end->rightMeanClipPos2;
						}
					}

					if(reg1 and reg2){ // merge
						clip_reg_retained = clip_reg;
						clip_reg_redudant = clip_reg_mate_end;
						if(end_flag1==1){
							// process double-region end
							if(clip_reg_retained->rightClipRegNum==2 and clip_reg_retained->rightClipReg->startRefPos>clip_reg_retained->rightClipReg2->startRefPos){
								reg_tmp = clip_reg_retained->rightClipReg;
								clip_reg_retained->rightClipReg = clip_reg_retained->rightClipReg2;
								clip_reg_retained->rightClipReg2 = reg_tmp;
								tmp = clip_reg_retained->rightClipPosNum;
								clip_reg_retained->rightClipPosNum = clip_reg_retained->rightClipPosNum2;
								clip_reg_retained->rightClipPosNum2 = tmp;
								tmp = clip_reg_retained->rightMeanClipPos;
								clip_reg_retained->rightMeanClipPos = clip_reg_retained->rightMeanClipPos2;
								clip_reg_retained->rightMeanClipPos2 = tmp;

								cerr << "************** line=" << __LINE__ << ", invalid swap, error!" << endl;
							}

							// process single-region end
							if(reg1->startRefPos<=reg2->startRefPos){
								clip_reg_retained->leftClipReg = reg1;
								clip_reg_retained->leftClipPosNum = pos_num1;
								clip_reg_retained->leftMeanClipPos = mean_pos1;
								clip_reg_retained->leftClipReg2 = reg2;
								clip_reg_retained->leftClipPosNum2 = pos_num2;
								clip_reg_retained->leftMeanClipPos2 = mean_pos2;
							}else{
								clip_reg_retained->leftClipReg = reg2;
								clip_reg_retained->leftClipPosNum = pos_num2;
								clip_reg_retained->leftMeanClipPos = mean_pos2;
								clip_reg_retained->leftClipReg2 = reg1;
								clip_reg_retained->leftClipPosNum2 = pos_num1;
								clip_reg_retained->leftMeanClipPos2 = mean_pos1;
							}
							clip_reg_retained->leftClipRegNum = 2;

						}else{
							// process double-region end
							if(clip_reg_retained->leftClipReg->startRefPos>clip_reg_retained->leftClipReg2->startRefPos){
								reg_tmp = clip_reg_retained->leftClipReg;
								clip_reg_retained->leftClipReg = clip_reg_retained->leftClipReg2;
								clip_reg_retained->leftClipReg2 = reg_tmp;
								tmp = clip_reg_retained->leftClipPosNum;
								clip_reg_retained->leftClipPosNum = clip_reg_retained->leftClipPosNum2;
								clip_reg_retained->leftClipPosNum2 = tmp;
								tmp = clip_reg_retained->leftMeanClipPos;
								clip_reg_retained->leftMeanClipPos = clip_reg_retained->leftMeanClipPos2;
								clip_reg_retained->leftMeanClipPos2 = tmp;
							}

							// process single-region end
							if(reg1->startRefPos<=reg2->startRefPos){
								clip_reg_retained->rightClipReg = reg1;
								clip_reg_retained->rightClipPosNum = pos_num1;
								clip_reg_retained->rightMeanClipPos = mean_pos1;
								clip_reg_retained->rightClipReg2 = reg2;
								clip_reg_retained->rightClipPosNum2 = pos_num2;
								clip_reg_retained->rightMeanClipPos2 = mean_pos2;
							}else{
								clip_reg_retained->rightClipReg = reg2;
								clip_reg_retained->rightClipPosNum = pos_num2;
								clip_reg_retained->rightMeanClipPos = mean_pos2;
								clip_reg_retained->rightClipReg2 = reg1;
								clip_reg_retained->rightClipPosNum2 = pos_num1;
								clip_reg_retained->rightMeanClipPos2 = mean_pos1;
							}
							clip_reg_retained->rightClipRegNum = 2;
						}

						// update redundant region information
						if(clip_reg_redudant->leftClipReg==reg2) { clip_reg_redudant->leftClipReg = NULL; clip_reg_redudant->leftClipRegNum--; }
						else if(clip_reg_redudant->leftClipReg2==reg2) { clip_reg_redudant->leftClipReg2 = NULL; clip_reg_redudant->leftClipRegNum--; }
						else if(clip_reg_redudant->rightClipReg==reg2) { clip_reg_redudant->rightClipReg = NULL; clip_reg_redudant->rightClipRegNum--; }
						else if(clip_reg_redudant->rightClipReg2==reg2) { clip_reg_redudant->rightClipReg2 = NULL; clip_reg_redudant->rightClipRegNum--; }
						else{
							cerr << "line=" << __LINE__ << ", invalid region, error!" << endl;
							exit(1);
						}
						clip_reg_redudant->valid_flag = false;
					}
				}else break;
			}
		}
	}

	// remove redundant node
	for(i=0; i<mateClipRegVector.size(); ){
		clip_reg = mateClipRegVector.at(i);
		if(clip_reg->valid_flag==false){
			if(clip_reg->leftClipReg) delete clip_reg->leftClipReg;
			if(clip_reg->rightClipReg) delete clip_reg->rightClipReg;
			if(clip_reg->leftClipReg2) delete clip_reg->leftClipReg2;
			if(clip_reg->rightClipReg2) delete clip_reg->rightClipReg2;
			delete clip_reg;
			mateClipRegVector.erase(mateClipRegVector.begin()+i);
		}else i++;
	}
}

// get other end region according to the same clipping region
mateClipReg_t* Chrome::getMateRegEndSameClipReg(mateClipReg_t *clip_reg, vector<mateClipReg_t*> &mate_clipReg_vec){
	mateClipReg_t *clip_reg_mate_end, *clip_reg_tmp;
	reg_t *reg1, *reg2, *reg1_tmp, *reg2_tmp;

	reg1 = reg2 = NULL;
	if(clip_reg->leftClipRegNum==2) {
		reg1 = clip_reg->leftClipReg;
		reg2 = clip_reg->leftClipReg2;
	}else if(clip_reg->rightClipRegNum==2){
		reg1 = clip_reg->rightClipReg;
		reg2 = clip_reg->rightClipReg2;
	}

	clip_reg_mate_end = NULL;
	for(size_t i=0; i<mate_clipReg_vec.size(); i++){
		clip_reg_tmp = mate_clipReg_vec.at(i);
		if(clip_reg_tmp!=clip_reg){
			if(clip_reg_tmp->valid_flag and clip_reg_tmp->reg_mated_flag and clip_reg_tmp->sv_type==clip_reg->sv_type){
				reg1_tmp = reg2_tmp = NULL;
				if(clip_reg_tmp->leftClipRegNum==2){
					reg1_tmp = clip_reg_tmp->leftClipReg;
					reg2_tmp = clip_reg_tmp->leftClipReg2;
				}else if(clip_reg->rightClipRegNum==2){
					reg1_tmp = clip_reg_tmp->rightClipReg;
					reg2_tmp = clip_reg_tmp->rightClipReg2;
				}

				if(reg1_tmp and reg2_tmp){
					if((isOverlappedReg(reg1, reg1_tmp) or isOverlappedReg(reg1, reg2_tmp)) and (isOverlappedReg(reg2, reg1_tmp) or isOverlappedReg(reg2, reg2_tmp))){
						clip_reg_mate_end = clip_reg_tmp;
						break;
					}
				}
			}
		}
	}

	return clip_reg_mate_end;
}

// remove FP clip regions
void Chrome::removeFPClipRegsDupInv(){
	mateClipReg_t *mate_clip_reg, *mate_clip_reg_overlapped;
	reg_t *reg, *reg2;
	size_t i, clipPosNum, clipPosNum_overlapped, left_pos_num, right_pos_num, left_pos_num_overlapped, right_pos_num_overlapped;
	int32_t dist;
	vector<size_t> overlap_result;

	// remove FP clip regions by detecting too long distance
	for(i=0; i<mateClipRegVector.size(); i++){
		mate_clip_reg = mateClipRegVector.at(i);
		if(mate_clip_reg->valid_flag and mate_clip_reg->reg_mated_flag){
			if(mate_clip_reg->sv_type==VAR_DUP or mate_clip_reg->sv_type==VAR_INV){
//				if(mate_clip_reg->leftClipRegNum!=1 or mate_clip_reg->rightClipRegNum!=1){
//					cerr << "line=" << __LINE__ << ", invalid clipping region!" << endl;
//					exit(1);
//				}
				reg = mate_clip_reg->leftClipReg ? mate_clip_reg->leftClipReg : mate_clip_reg->leftClipReg2;
				reg2 = mate_clip_reg->rightClipReg ? mate_clip_reg->rightClipReg : mate_clip_reg->rightClipReg2;
				if(reg->chrname.compare(reg2->chrname)==0){  // same chrome
					dist = reg2->startRefPos - reg->startRefPos;
					if(dist<0) dist = -dist;
					if(dist>paras->maxVarRegSize) mate_clip_reg->valid_flag = false;
					if(reg->startRefPos>reg2->endRefPos) mate_clip_reg->valid_flag = false;
				} // different chrome
			}
		}else if(mate_clip_reg->large_indel_flag==false) mate_clip_reg->valid_flag = false;
	}

	// remove FP clip regions by detecting false overlapped regions
	for(i=0; i<mateClipRegVector.size(); ){
		mate_clip_reg = mateClipRegVector.at(i);
		if(mate_clip_reg->valid_flag and mate_clip_reg->reg_mated_flag and (mate_clip_reg->sv_type==VAR_DUP or mate_clip_reg->sv_type==VAR_INV)){
			mate_clip_reg_overlapped = getOverlappedMateClipReg(mate_clip_reg, mateClipRegVector);
			if(mate_clip_reg_overlapped){
				left_pos_num = mate_clip_reg->leftClipReg ? mate_clip_reg->leftClipPosNum : mate_clip_reg->leftClipPosNum2;
				right_pos_num = mate_clip_reg->rightClipReg ? mate_clip_reg->rightClipPosNum : mate_clip_reg->rightClipPosNum2;
				clipPosNum = left_pos_num + right_pos_num;

				left_pos_num_overlapped = mate_clip_reg_overlapped->leftClipReg ? mate_clip_reg_overlapped->leftClipPosNum : mate_clip_reg_overlapped->leftClipPosNum2;
				right_pos_num_overlapped = mate_clip_reg_overlapped->rightClipReg ? mate_clip_reg_overlapped->rightClipPosNum : mate_clip_reg_overlapped->rightClipPosNum2;
				clipPosNum_overlapped = left_pos_num_overlapped + right_pos_num_overlapped;
				if(clipPosNum>=clipPosNum_overlapped)
					mate_clip_reg_overlapped->valid_flag = false;
				else{
					mate_clip_reg->valid_flag = false;
					i++;
				}
			}else
				i++;
		}else i++;
	}

//	for(i=0; i<mateClipRegVector.size(); i++){
//		mate_clip_reg = mateClipRegVector.at(i);
//
//		for(j=0; j<2; j++){
//			if(j==0){ // left region
//				reg = mate_clip_reg->leftClipReg;
//				clipPosNum = mate_clip_reg->leftClipPosNum;
//			}else{  // right region
//				reg = mate_clip_reg->rightClipReg;
//				clipPosNum = mate_clip_reg->rightClipPosNum;
//			}
//
//			if(reg){
//				overlap_result = getOverlapClipReg(reg);
//				if(overlap_result.size()>0){
//					idx = overlap_result.at(0);
//					end_flag = overlap_result.at(1);
//					if(end_flag==LEFT_END) clipPosNum_overlap = mateClipRegVector.at(idx)->leftClipPosNum;
//					else clipPosNum_overlap = mateClipRegVector.at(idx)->rightClipPosNum;
//
//					if(clipPosNum>=clipPosNum_overlap) mateClipRegVector.at(idx)->valid_flag = false;
//					else mateClipRegVector.at(i)->valid_flag = false;
//				}
//			}
//		}
//	}

	// remove invalid items
	for(i=0; i<mateClipRegVector.size(); ){
		mate_clip_reg = mateClipRegVector.at(i);
		if(mate_clip_reg->valid_flag==false){
			if(mate_clip_reg->leftClipReg) delete mate_clip_reg->leftClipReg;
			if(mate_clip_reg->leftClipReg2) delete mate_clip_reg->leftClipReg2;
			if(mate_clip_reg->rightClipReg) delete mate_clip_reg->rightClipReg;
			if(mate_clip_reg->rightClipReg2) delete mate_clip_reg->rightClipReg2;
			delete mate_clip_reg;
			mateClipRegVector.erase(mateClipRegVector.begin()+i);
		}else i++;
	}
}

// get the overlapped clip region
vector<size_t> Chrome::getOverlapClipReg(reg_t *given_reg){
	vector<size_t> overlap_result;
	reg_t *reg;
	for(size_t i=0; i<mateClipRegVector.size(); i++){
		reg = mateClipRegVector.at(i)->leftClipReg;
		if(reg!=given_reg and isOverlappedReg(reg, given_reg)){
			overlap_result.push_back(i);
			overlap_result.push_back(LEFT_END);
			break;
		}

		reg = mateClipRegVector.at(i)->rightClipReg;
		if(reg!=given_reg and isOverlappedReg(reg, given_reg)){
			overlap_result.push_back(i);
			overlap_result.push_back(RIGHT_END);
			break;
		}
	}
	return overlap_result;
}

// remove FP indels in clipping regions
void Chrome::removeFPIndelSnvInClipReg(vector<mateClipReg_t*> &mate_clipReg_vec){
	size_t i, j, pos;
	Block *block;
	reg_t *reg;
	bool flag;

	//size_t num = 0;
	for(i=0; i<blockVector.size(); i++){
		block = blockVector.at(i);
		for(j=0; j<block->indelVector.size(); ){
			reg = block->indelVector.at(j);
			flag = isIndelInClipReg(reg, mate_clipReg_vec);
			if(flag){
				//num ++;
				//printRegVec(block->indelVector, "indelVector");
				delete reg;
				block->indelVector.erase(block->indelVector.begin()+j);
			}else j++;
		}
		for(j=0; j<block->snvVector.size(); ){
			pos = block->snvVector.at(j);
			flag = isSnvInClipReg(chrname, pos, mate_clipReg_vec);
			if(flag) block->snvVector.erase(block->snvVector.begin()+j);
			else j++;
		}
	}

	//cout << "+++++++++++++++++ num=" << num << endl;
}

// merge detect result to single file
void Chrome::chrMergeDetectResultToFile(){
	size_t i, j, pos;
	ofstream out_file_snv, out_file_indel, out_file_clipReg;
	string filename, line, var_type, mate_str, tmp_str;
	vector<reg_t*> indel_vec;
	vector<size_t> snv_vec;
	mateClipReg_t *mate_clip_reg;
	reg_t *reg;
	Block* bloc;

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

	for(i=0; i<blockVector.size(); i++){
		bloc = blockVector.at(i);
		indel_vec = bloc->indelVector;
		for(j=0; j<indel_vec.size(); j++){
			reg = indel_vec.at(j);
			out_file_indel << reg->chrname << "\t" << reg->startRefPos << "\t" <<  reg->endRefPos << "\t" << reg->sv_len << endl;
		}
		snv_vec = bloc->snvVector;
		for(j=0; j<snv_vec.size(); j++){
			pos = snv_vec.at(j);
			out_file_snv << chrname << "\t" << pos << endl;
		}
	}
	out_file_snv.close();
	out_file_indel.close();

	// clipping regions
	out_file_clipReg.open(out_filename_detect_clipReg);
	if(!out_file_clipReg.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << out_filename_detect_clipReg << endl;
		exit(1);
	}
	for(i=0; i<mateClipRegVector.size(); i++){
		mate_clip_reg = mateClipRegVector.at(i);
		if(mate_clip_reg->valid_flag){
			if(mate_clip_reg->large_indel_flag==false){
				reg = mate_clip_reg->leftClipReg;
				line = "";
				if(reg) line += reg->chrname + "\t" + to_string(reg->startRefPos) + "\t" + to_string(reg->endRefPos);
				else line += "-\t-\t-";

				reg = mate_clip_reg->leftClipReg2;
				if(reg) line += "\t" + reg->chrname + "\t" + to_string(reg->startRefPos) + "\t" + to_string(reg->endRefPos);
				else line += "\t-\t-\t-";

				reg = mate_clip_reg->rightClipReg;
				if(reg) line += "\t" + reg->chrname + "\t" + to_string(reg->startRefPos) + "\t" + to_string(reg->endRefPos);
				else line += "\t-\t-\t-";

				reg = mate_clip_reg->rightClipReg2;
				if(reg) line += "\t" + reg->chrname + "\t" + to_string(reg->startRefPos) + "\t" + to_string(reg->endRefPos);
				else line += "\t-\t-\t-";

				if(mateClipRegVector.at(i)->reg_mated_flag==true) line += "\t1";
				else line += "\t0";

				switch(mate_clip_reg->sv_type){
					case VAR_UNC: var_type = "UNC"; break;
					case VAR_DUP: var_type = "DUP"; break;
					case VAR_INV: var_type = "INV"; break;
					case VAR_TRA: var_type = "TRA"; break;
					case VAR_MIX: var_type = "MIX"; break;
					default: cerr << __func__ << ", line=" << __LINE__ << ": invalid var_type: " << mate_clip_reg->sv_type << endl; exit(1);
				}

				line += "\t####\t" + to_string(mate_clip_reg->leftMeanClipPos) + "\t" + to_string(mate_clip_reg->leftMeanClipPos2) + "\t" + to_string(mate_clip_reg->rightMeanClipPos) + "\t" + to_string(mate_clip_reg->rightMeanClipPos2);
				line += "\t" + var_type;

				if(mate_clip_reg->sv_type==VAR_DUP)
					line += "\t" + to_string(mate_clip_reg->dup_num);
				else
					line += "\t-";

				mate_str = mate_clip_reg->bnd_mate_reg_strs[0];
				for(j=1; j<4; j++) mate_str += ";" + mate_clip_reg->bnd_mate_reg_strs[j];
				line += "\t" + mate_str;

				line += "\t" + to_string(mate_clip_reg->leftClipPosNum) + "\t" + to_string(mate_clip_reg->leftClipPosNum2) + "\t" + to_string(mate_clip_reg->rightClipPosNum) + "\t" + to_string(mate_clip_reg->rightClipPosNum2);
			}else{ // large indel
				reg = mate_clip_reg->largeIndelClipReg;
				line = "";
				if(reg) line += reg->chrname + "\t" + to_string(reg->startRefPos) + "\t" + to_string(reg->endRefPos);// left_reg
				else line += "-\t-\t-";
				line += "\t-\t-\t-"; // left_reg2
				line += "\t-\t-\t-"; // right_reg
				line += "\t-\t-\t-"; // right_reg2
				line += "\t0"; // reg_mated_flag

				switch(mate_clip_reg->sv_type){
					case VAR_UNC: var_type = "UNC"; break;
					case VAR_INS: var_type = "INS"; break; //
					case VAR_DEL: var_type = "DEL"; break; //
					case VAR_DUP: var_type = "DUP"; break;
					case VAR_INV: var_type = "INV"; break;
					case VAR_TRA: var_type = "TRA"; break;
					case VAR_MIX: var_type = "MIX"; break;
					default: cerr << __func__ << ", line=" << __LINE__ << ": invalid var_type: " << mate_clip_reg->sv_type << endl; exit(1);
				}

				line += "\t####\t0\t0\t0\t0"; // mean clippos
				line += "\t" + var_type;

				if(mate_clip_reg->sv_type==VAR_DUP)
					line += "\t" + to_string(mate_clip_reg->dup_num);
				else
					line += "\t-";

				mate_str = mate_clip_reg->bnd_mate_reg_strs[0];
				for(j=1; j<4; j++) mate_str += ";" + mate_clip_reg->bnd_mate_reg_strs[j];
				line += "\t" + mate_str;

				line += "\t" + to_string(reg->sv_len) + "\t" + to_string(mate_clip_reg->supp_num_largeIndel) + "\t" + to_string(mate_clip_reg->depth_largeIndel) + "\t0";
			}

			out_file_clipReg << line << endl;
			//cout << line << endl;
		}
	}
	out_file_clipReg.close();
}

void Chrome::chrMergeDetectResultToFile_debug(){
	size_t i, j, pos;
	ofstream out_file_snv, out_file_indel, out_file_clipReg;
	string filename, line, var_type, mate_str, tmp_str;
	string out_filename_detect_snv_tmp, out_filename_detect_indel_tmp, out_filename_detect_clipReg_tmp;
	vector<reg_t*> indel_vec;
	vector<size_t> snv_vec;
	mateClipReg_t *mate_clip_reg;
	reg_t *reg;
	Block* bloc;

	out_filename_detect_snv_tmp = out_filename_detect_snv + "_tmp";
	out_filename_detect_indel_tmp = out_filename_detect_indel + "_tmp";
	out_filename_detect_clipReg_tmp = out_filename_detect_clipReg + "_tmp";

	out_file_snv.open(out_filename_detect_snv_tmp);
	if(!out_file_snv.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << out_filename_detect_snv << endl;
		exit(1);
	}
	out_file_indel.open(out_filename_detect_indel_tmp);
	if(!out_file_indel.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << out_filename_detect_indel << endl;
		exit(1);
	}

	for(i=0; i<blockVector.size(); i++){
		bloc = blockVector.at(i);
		indel_vec = bloc->indelVector;
		for(j=0; j<indel_vec.size(); j++){
			reg = indel_vec.at(j);
			out_file_indel << reg->chrname << "\t" << reg->startRefPos << "\t" <<  reg->endRefPos << endl;
		}
		snv_vec = bloc->snvVector;
		for(j=0; j<snv_vec.size(); j++){
			pos = snv_vec.at(j);
			out_file_snv << chrname << "\t" << pos << endl;
		}
	}
	out_file_snv.close();
	out_file_indel.close();

	// clipping regions
	out_file_clipReg.open(out_filename_detect_clipReg_tmp);
	if(!out_file_clipReg.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << out_filename_detect_clipReg << endl;
		exit(1);
	}
	for(i=0; i<mateClipRegVector.size(); i++){
		mate_clip_reg = mateClipRegVector.at(i);
		reg = mate_clip_reg->leftClipReg;
		line = "";
		if(reg) line += reg->chrname + "\t" + to_string(reg->startRefPos) + "\t" + to_string(reg->endRefPos);
		else line += "-\t-\t-";

		reg = mate_clip_reg->leftClipReg2;
		if(reg) line += "\t" + reg->chrname + "\t" + to_string(reg->startRefPos) + "\t" + to_string(reg->endRefPos);
		else line += "\t-\t-\t-";

		reg = mate_clip_reg->rightClipReg;
		if(reg) line += "\t" + reg->chrname + "\t" + to_string(reg->startRefPos) + "\t" + to_string(reg->endRefPos);
		else line += "\t-\t-\t-";

		reg = mate_clip_reg->rightClipReg2;
		if(reg) line += "\t" + reg->chrname + "\t" + to_string(reg->startRefPos) + "\t" + to_string(reg->endRefPos);
		else line += "\t-\t-\t-";

		if(mateClipRegVector.at(i)->reg_mated_flag==true) line += "\t1";
		else line += "\t0";

		switch(mate_clip_reg->sv_type){
			case VAR_UNC: var_type = "UNC"; break;
			case VAR_DUP: var_type = "DUP"; break;
			case VAR_INV: var_type = "INV"; break;
			case VAR_TRA: var_type = "TRA"; break;
			case VAR_MIX: var_type = "MIX"; break;
			default: cerr << __func__ << ", line=" << __LINE__ << ": invalid var_type: " << mate_clip_reg->sv_type << endl; exit(1);
		}

		line += "\t####\t" + to_string(mate_clip_reg->leftMeanClipPos) + "\t" + to_string(mate_clip_reg->leftMeanClipPos2) + "\t" + to_string(mate_clip_reg->rightMeanClipPos) + "\t" + to_string(mate_clip_reg->rightMeanClipPos2);
		line += "\t" + var_type;

		if(mate_clip_reg->sv_type==VAR_DUP)
			line += "\t" + to_string(mate_clip_reg->dup_num);
		else
			line += "\t-";

		mate_str = mate_clip_reg->bnd_mate_reg_strs[0];
		for(j=1; j<4; j++) mate_str += ";" + mate_clip_reg->bnd_mate_reg_strs[j];
		line += "\t" + mate_str;

		line += "\t" + to_string(mate_clip_reg->leftClipPosNum) + "\t" + to_string(mate_clip_reg->leftClipPosNum2) + "\t" + to_string(mate_clip_reg->rightClipPosNum) + "\t" + to_string(mate_clip_reg->rightClipPosNum2);

		out_file_clipReg << line << endl;
		//cout << line << endl;
	}
	out_file_clipReg.close();
}

// set file information  for local consensus
void Chrome::chrSetVarCandFiles(){
	string line, new_line, tmp_filename, header_line, old_out_dir, refseqfilename, contigfilename, readsfilename, pattern_str, limit_reg_str;
	vector<string> str_vec, var_str, var_str1, var_str2;
	ifstream infile;
	size_t i;
	simpleReg_t *simple_reg, *prev_simple_reg;
	vector<simpleReg_t*> sub_limit_reg_vec, prev_limit_reg_vec, prev_limit_reg_vec_tmp, pos_limit_reg_vec;
	bool flag, pos_contained_flag;

	// indel
	if(isFileExist(var_cand_indel_filename)){

		if(paras->limit_reg_process_flag){
			pattern_str = REFSEQ_PATTERN;
			simple_reg = new simpleReg_t();
		}

		tmp_filename = var_cand_indel_filename + "_tmp";
		rename(var_cand_indel_filename.c_str(), tmp_filename.c_str());

		infile.open(tmp_filename);
		if(!infile.is_open()){
			cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << tmp_filename << endl;
			exit(1);
		}

		var_cand_indel_file.open(var_cand_indel_filename);
		if(!var_cand_indel_file.is_open()){
			cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << var_cand_indel_filename << endl;
			exit(1);
		}

		header_line = getCnsFileHeaderLine();
		var_cand_indel_file << header_line << endl;

		old_out_dir = "";
		while(getline(infile, line)){
			if(line.size()>0 and line.at(0)!='#'){
				str_vec = split(line, "\t");

				if(str_vec.at(str_vec.size()-1).compare(DONE_STR)==0){
					if(old_out_dir.size()==0) old_out_dir = getOldOutDirname(str_vec.at(0), paras->out_dir_cns);

					refseqfilename = getUpdatedItemFilename(str_vec.at(0), paras->outDir, old_out_dir);
					contigfilename = getUpdatedItemFilename(str_vec.at(1), paras->outDir, old_out_dir);
					readsfilename = getUpdatedItemFilename(str_vec.at(2), paras->outDir, old_out_dir);

					flag = true;
					pos_contained_flag = false;
					if(paras->limit_reg_process_flag) {
						getRegByFilename(simple_reg, refseqfilename, pattern_str);
						//sub_limit_reg_vec = getSimpleRegs(simple_reg->chrname, simple_reg->startPos, simple_reg->endPos, paras->limit_reg_vec);
						sub_limit_reg_vec = getOverlappedSimpleRegs(simple_reg->chrname, simple_reg->startPos, simple_reg->endPos, paras->limit_reg_vec);
						if(sub_limit_reg_vec.size()==0){
							flag = false;
							prev_limit_reg_vec = extractSimpleRegsByStr(str_vec.at(8));
							for(i=0; i<prev_limit_reg_vec.size(); i++){
								prev_simple_reg = prev_limit_reg_vec.at(i);
								//prev_limit_reg_vec_tmp = getSimpleRegs(prev_simple_reg->chrname, prev_simple_reg->startPos, prev_simple_reg->endPos, paras->limit_reg_vec);
								prev_limit_reg_vec_tmp = getFullyContainedSimpleRegs(prev_simple_reg->chrname, prev_simple_reg->startPos, prev_simple_reg->endPos, paras->limit_reg_vec);
								if(prev_limit_reg_vec_tmp.size()){
									flag = true;
									break;
								}
							}

							if(flag==false){ // position contained
								pos_limit_reg_vec = getPosContainedSimpleRegs(simple_reg->chrname, simple_reg->startPos, simple_reg->endPos, paras->limit_reg_vec);
								if(pos_limit_reg_vec.size()){
									flag = true;
									pos_contained_flag = true;
								}
							}
						}
					}

					if(flag){
						new_line = refseqfilename + "\t" + contigfilename + "\t" + readsfilename;
						for(i=3; i<=7; i++) new_line += "\t" + str_vec.at(i);

						// limit regions
						if(paras->limit_reg_process_flag){
							if(sub_limit_reg_vec.size()) limit_reg_str = getLimitRegStr(sub_limit_reg_vec);
							else {
								if(pos_contained_flag==false) limit_reg_str = getLimitRegStr(prev_limit_reg_vec);
								else limit_reg_str = getLimitRegStr(pos_limit_reg_vec);
							}
						}else limit_reg_str = LIMIT_REG_ALL_STR;
						new_line += "\t" + limit_reg_str;

						// other fields
						for(i=9; i<str_vec.size(); i++) new_line += "\t" + str_vec.at(i);

						var_cand_indel_file << new_line << endl;
					}

					if(!prev_limit_reg_vec.empty()) destroyLimitRegVector(prev_limit_reg_vec);
				}else
					cout << "line=" << __LINE__ << "," << var_cand_indel_filename << ": line does not end with 'DONE', skipped!" << endl;
			}
		}

		infile.close();
		remove(tmp_filename.c_str());
		if(paras->limit_reg_process_flag) delete simple_reg;
	}else{
		var_cand_indel_file.open(var_cand_indel_filename);
		if(!var_cand_indel_file.is_open()){
			cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << var_cand_indel_filename << endl;
			exit(1);
		}

		header_line = getCnsFileHeaderLine();
		var_cand_indel_file << header_line << endl;
	}

	// clipReg
	if(isFileExist(var_cand_clipReg_filename)){

		if(paras->limit_reg_process_flag){
			pattern_str = CLIPREG_PATTERN;
			simple_reg = new simpleReg_t();
		}

		tmp_filename = var_cand_clipReg_filename + "_tmp";
		rename(var_cand_clipReg_filename.c_str(), tmp_filename.c_str());

		infile.open(tmp_filename);
		if(!infile.is_open()){
			cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << tmp_filename << endl;
			exit(1);
		}

		var_cand_clipReg_file.open(var_cand_clipReg_filename);
		if(!var_cand_clipReg_file.is_open()){
			cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << var_cand_clipReg_filename << endl;
			exit(1);
		}

		header_line = getCnsFileHeaderLine();
		var_cand_clipReg_file << header_line << endl;

		old_out_dir = "";
		while(getline(infile, line)){
			if(line.size()>0 and line.at(0)!='#'){
				str_vec = split(line, "\t");

				if(str_vec.at(str_vec.size()-1).compare(DONE_STR)==0){
					if(old_out_dir.size()==0) old_out_dir = getOldOutDirname(str_vec.at(0), paras->out_dir_cns);

					refseqfilename = getUpdatedItemFilename(str_vec.at(0), paras->outDir, old_out_dir);
					contigfilename = getUpdatedItemFilename(str_vec.at(1), paras->outDir, old_out_dir);
					readsfilename = getUpdatedItemFilename(str_vec.at(2), paras->outDir, old_out_dir);

					flag = true;
					pos_contained_flag = false;
					if(paras->limit_reg_process_flag) {
						getRegByFilename(simple_reg, refseqfilename, pattern_str);
						//sub_limit_reg_vec = getSimpleRegs(simple_reg->chrname, simple_reg->startPos, simple_reg->endPos, paras->limit_reg_vec);
						sub_limit_reg_vec = getOverlappedSimpleRegs(simple_reg->chrname, simple_reg->startPos, simple_reg->endPos, paras->limit_reg_vec);
						if(sub_limit_reg_vec.size()==0) {
							flag = false;
							prev_limit_reg_vec = extractSimpleRegsByStr(str_vec.at(8));
							for(i=0; i<prev_limit_reg_vec.size(); i++){
								prev_simple_reg = prev_limit_reg_vec.at(i);
								//prev_limit_reg_vec_tmp = getSimpleRegs(prev_simple_reg->chrname, prev_simple_reg->startPos, prev_simple_reg->endPos, paras->limit_reg_vec);
								prev_limit_reg_vec_tmp = getFullyContainedSimpleRegs(prev_simple_reg->chrname, prev_simple_reg->startPos, prev_simple_reg->endPos, paras->limit_reg_vec);
								if(prev_limit_reg_vec_tmp.size()){
									flag = true;
									break;
								}
							}
							if(flag==false){ // position contained
								pos_limit_reg_vec = getPosContainedSimpleRegs(simple_reg->chrname, simple_reg->startPos, simple_reg->endPos, paras->limit_reg_vec);
								if(pos_limit_reg_vec.size()){
									flag = true;
									pos_contained_flag = true;
								}
							}
						}
					}

					if(flag){
						new_line = refseqfilename + "\t" + contigfilename + "\t" + readsfilename;
						for(i=3; i<=7; i++) new_line += "\t" + str_vec.at(i);

						// limit regions
						if(paras->limit_reg_process_flag){
							if(sub_limit_reg_vec.size()) limit_reg_str = getLimitRegStr(sub_limit_reg_vec);
							else {
								if(pos_contained_flag==false) limit_reg_str = getLimitRegStr(prev_limit_reg_vec);
								else limit_reg_str = getLimitRegStr(pos_limit_reg_vec);
							}
						}else limit_reg_str = LIMIT_REG_ALL_STR;
						new_line += "\t" + limit_reg_str;

						// other fields
						for(i=9; i<str_vec.size(); i++) new_line += "\t" + str_vec.at(i);

						var_cand_clipReg_file << new_line << endl;
					}

					if(!prev_limit_reg_vec.empty()) destroyLimitRegVector(prev_limit_reg_vec);
				}else
					cout << "line=" << __LINE__ << "," << var_cand_indel_filename << ": line does not end with 'DONE'!" << endl;
			}
		}

		infile.close();
		remove(tmp_filename.c_str());
		if(paras->limit_reg_process_flag) delete simple_reg;
	}else{
		var_cand_clipReg_file.open(var_cand_clipReg_filename);
		if(!var_cand_clipReg_file.is_open()){
			cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << var_cand_clipReg_filename << endl;
			exit(1);
		}

		header_line = getCnsFileHeaderLine();
		var_cand_clipReg_file << header_line << endl;
	}

	vector<Block*>::iterator bloc;
	for(bloc=blockVector.begin(); bloc!=blockVector.end(); bloc++)
		(*bloc)->setVarCandFiles(&var_cand_indel_file, &var_cand_clipReg_file);
}

// reset file information for local consensus
void Chrome::chrResetVarCandFiles(){
	var_cand_indel_file.close();
	var_cand_clipReg_file.close();

	vector<Block*>::iterator bloc;
	for(bloc=blockVector.begin(); bloc!=blockVector.end(); bloc++)
		(*bloc)->resetVarCandFiles();
}

// set misAln region file
void Chrome::chrSetMisAlnRegFile(){
	misAln_reg_file.open(misAln_reg_filename);
	if(!misAln_reg_file.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << misAln_reg_filename << endl;
		exit(1);
	}

	vector<Block*>::iterator bloc;
	for(bloc=blockVector.begin(); bloc!=blockVector.end(); bloc++)
		(*bloc)->setMisAlnRegFile(&misAln_reg_file);
}

// reset misAln region file
void Chrome::chrResetMisAlnRegFile(){
	misAln_reg_file.close();

	vector<Block*>::iterator bloc;
	for(bloc=blockVector.begin(); bloc!=blockVector.end(); bloc++)
		(*bloc)->resetMisAlnRegFile();
}

// load detected variant data for local consensus
void Chrome::chrLoadDataCons(){
	size_t i;
	mateClipReg_t* mate_clip_reg;
	Block* tmp_bloc;

	if(paras->load_from_file_flag){ // load data from file
		chrLoadIndelDataCons();
		chrLoadClipRegDataCons();
	}

	// initialize block information
	for(i=0; i<mateClipRegVector.size(); i++){
		mate_clip_reg = mateClipRegVector.at(i);
		tmp_bloc = NULL;
		if(mate_clip_reg->valid_flag){
			if(mate_clip_reg->large_indel_flag==false){ // clip region
				if(mate_clip_reg->leftClipReg and mate_clip_reg->leftClipReg->chrname.compare(chrname)==0) tmp_bloc = computeBlocByPos(mate_clip_reg->leftClipReg->startRefPos, blockVector);  // get the block
				else if(mate_clip_reg->rightClipReg and mate_clip_reg->rightClipReg->chrname.compare(chrname)==0) tmp_bloc = computeBlocByPos(mate_clip_reg->rightClipReg->startRefPos, blockVector);  // get the block
				else if(mate_clip_reg->leftClipReg2 and mate_clip_reg->leftClipReg2->chrname.compare(chrname)==0) tmp_bloc = computeBlocByPos(mate_clip_reg->leftClipReg2->startRefPos, blockVector);  // get the block
				else if(mate_clip_reg->rightClipReg2 and mate_clip_reg->rightClipReg2->chrname.compare(chrname)==0) tmp_bloc = computeBlocByPos(mate_clip_reg->rightClipReg2->startRefPos, blockVector);  // get the block
			}else{ // large indel region
				if(mate_clip_reg->largeIndelClipReg and mate_clip_reg->largeIndelClipReg->chrname.compare(chrname)==0) tmp_bloc = computeBlocByPos(mate_clip_reg->largeIndelClipReg->startRefPos, blockVector);  // get the block
			}
			if(tmp_bloc) tmp_bloc->mateClipRegVector.push_back(mate_clip_reg);
		}
	}
}

// load detected indel data for local consensus
void Chrome::chrLoadIndelDataCons(){
	vector<simpleReg_t*> limit_reg_vec;
	//if(paras->limit_reg_process_flag) limit_reg_vec = getSimpleRegs(chrname, -1, -1, paras->limit_reg_vec);
	if(paras->limit_reg_process_flag) limit_reg_vec = getOverlappedSimpleRegs(chrname, -1, -1, paras->limit_reg_vec);
	chrLoadIndelData(paras->limit_reg_process_flag, limit_reg_vec);
}

// load detected indel data
void Chrome::chrLoadIndelData(bool limit_reg_process_flag, vector<simpleReg_t*> &limit_reg_vec){
	string line;
	vector<string> str_vec;
	ifstream infile;
	int64_t begPos, endPos, sv_len;
	Block* tmp_bloc;
	reg_t *reg;
	simpleReg_t *simple_reg;
	vector<simpleReg_t*> sub_limit_reg_vec;
	bool flag;

	if(isFileExist(out_filename_detect_indel)==false) return;

	infile.open(out_filename_detect_indel);
	if(!infile.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << out_filename_detect_indel << endl;
		exit(1);
	}

	if(limit_reg_process_flag) simple_reg = new simpleReg_t();

	while(getline(infile, line)){
		if(line.size()>0 and line.at(0)!='#'){
			str_vec = split(line, "\t");
			begPos = stoi(str_vec.at(1));
			endPos = stoi(str_vec.at(2));
			sv_len = stoi(str_vec.at(3));
			tmp_bloc = computeBlocByPos(begPos, blockVector);  // get the block

			// deal with limit regions
			flag = true;
			if(limit_reg_process_flag) {
				//sub_limit_reg_vec = getSimpleRegs(chrname, begPos, endPos, limit_reg_vec);
				sub_limit_reg_vec = getOverlappedSimpleRegsExt(chrname, begPos, endPos, limit_reg_vec, CNS_SIDE_EXT_SIZE);
				if(sub_limit_reg_vec.size()==0) flag = false;
			}

			if(flag){
				reg = new reg_t();
				reg->chrname = chrname;
				reg->startRefPos = begPos;
				reg->endRefPos = endPos;
				reg->startLocalRefPos = reg->endLocalRefPos = 0;
				reg->startQueryPos = reg->endQueryPos = 0;
				reg->var_type = VAR_UNC;
				reg->sv_len = sv_len;
				reg->query_id = -1;
				reg->blat_aln_id = -1;
				reg->minimap2_aln_id = -1;
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
				//cout << "blocID=" << computeBlocID(begPos, blockVector) << ", reg:" << begPos << "-" << endPos << endl;

				tmp_bloc->indelVector.push_back(reg);  // add the variant
			}
		}
	}
	infile.close();
	if(limit_reg_process_flag) delete simple_reg;
}

// load detected indel data for local consensus
void Chrome::chrLoadClipRegDataCons(){
	chrLoadMateClipRegDataOp(paras->limit_reg_process_flag, paras->limit_reg_vec);
}

void Chrome::chrLoadMateClipRegDataOp(bool limit_reg_process_flag, vector<simpleReg_t*> &limit_reg_vec){
	string line, chrname1, chrname2, chrname3, chrname4;
	vector<string> str_vec, bnd_mate_reg_str_vec;
	ifstream infile;
	reg_t *reg1, *reg2, *reg3, *reg4;
	bool mate_flag, flag, flag1, flag2, flag3, flag4;
	mateClipReg_t* mate_clip_reg;
	size_t var_type, left_size, left_size2, right_size, right_size2, dup_num_tmp, rd_num_left, rd_num_left2, rd_num_right, rd_num_right2;
	//simpleReg_t *simple_reg;
	vector<simpleReg_t*> sub_limit_reg_vec;

	if(isFileExist(out_filename_detect_clipReg)==false) return;

	infile.open(out_filename_detect_clipReg);
	if(!infile.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << out_filename_detect_clipReg << endl;
		exit(1);
	}

	//if(limit_reg_process_flag) simple_reg = new simpleReg_t();

	while(getline(infile, line)){
		if(line.size()>0 and line.at(0)!='#'){
			reg1 = reg2 = reg3 = reg4 = NULL;
			str_vec = split(line, "\t");

			chrname1 = str_vec.at(0);
			if(chrname1.compare("-")!=0) { reg1 = new reg_t(); reg1->chrname = chrname1; reg1->startRefPos = stoi(str_vec.at(1)); reg1->endRefPos = stoi(str_vec.at(2)); reg1->zero_cov_flag = false; reg1->aln_seg_end_flag = false; reg1->query_pos_invalid_flag = false; reg1->large_indel_flag = false; reg1->gt_type = -1; reg1->gt_seq = ""; reg1->AF = 0; reg1->supp_num = reg1->DP = 0; reg1->discover_level = VAR_DISCOV_L_UNUSED; }
			chrname2 = str_vec.at(3);
			if(chrname2.compare("-")!=0) { reg2 = new reg_t(); reg2->chrname = chrname2; reg2->startRefPos = stoi(str_vec.at(4)); reg2->endRefPos = stoi(str_vec.at(5)); reg2->zero_cov_flag = false; reg2->aln_seg_end_flag = false; reg2->query_pos_invalid_flag = false; reg2->large_indel_flag = false; reg2->gt_type = -1; reg2->gt_seq = ""; reg2->AF = 0; reg2->supp_num = reg2->DP = 0; reg2->discover_level = VAR_DISCOV_L_UNUSED; }
			chrname3 = str_vec.at(6);
			if(chrname3.compare("-")!=0) { reg3 = new reg_t(); reg3->chrname = chrname3; reg3->startRefPos = stoi(str_vec.at(7)); reg3->endRefPos = stoi(str_vec.at(8)); reg3->zero_cov_flag = false; reg3->aln_seg_end_flag = false; reg3->query_pos_invalid_flag = false; reg3->large_indel_flag = false; reg3->gt_type = -1; reg3->gt_seq = ""; reg3->AF = 0; reg3->supp_num = reg3->DP = 0; reg3->discover_level = VAR_DISCOV_L_UNUSED; }
			chrname4 = str_vec.at(9);
			if(chrname4.compare("-")!=0) { reg4 = new reg_t(); reg4->chrname = chrname4; reg4->startRefPos = stoi(str_vec.at(10)); reg4->endRefPos = stoi(str_vec.at(11)); reg4->zero_cov_flag = false; reg4->aln_seg_end_flag = false; reg4->query_pos_invalid_flag = false; reg4->large_indel_flag = false; reg4->gt_type = -1; reg4->gt_seq = ""; reg4->AF = 0; reg4->supp_num = reg4->DP = 0; reg4->discover_level = VAR_DISCOV_L_UNUSED; }

			if(str_vec.at(12).compare("0")==0) mate_flag = false;
			else mate_flag = true;

			// deal with limit regions
			flag1 = flag2 = flag3 = flag4 = true;
			if(reg1){
				if(limit_reg_process_flag) {
					//sub_limit_reg_vec = getSimpleRegs(reg1->chrname, reg1->startRefPos, reg1->endRefPos, limit_reg_vec);
					sub_limit_reg_vec = getOverlappedSimpleRegsExt(reg1->chrname, reg1->startRefPos, reg1->endRefPos, limit_reg_vec, CNS_SIDE_EXT_SIZE);
					if(sub_limit_reg_vec.size()==0) flag1 = false;
				}
			}else flag1 = false;
			if(reg2){
				if(limit_reg_process_flag) {
					//sub_limit_reg_vec = getSimpleRegs(reg2->chrname, reg2->startRefPos, reg2->endRefPos, limit_reg_vec);
					sub_limit_reg_vec = getOverlappedSimpleRegsExt(reg2->chrname, reg2->startRefPos, reg2->endRefPos, limit_reg_vec, CNS_SIDE_EXT_SIZE);
					if(sub_limit_reg_vec.size()==0) flag2 = false;
				}
			}else flag2 = false;
			if(reg3){
				if(limit_reg_process_flag) {
					//sub_limit_reg_vec = getSimpleRegs(reg3->chrname, reg3->startRefPos, reg3->endRefPos, limit_reg_vec);
					sub_limit_reg_vec = getOverlappedSimpleRegsExt(reg3->chrname, reg3->startRefPos, reg3->endRefPos, limit_reg_vec, CNS_SIDE_EXT_SIZE);
					if(sub_limit_reg_vec.size()==0) flag3 = false;
				}
			}else flag3 = false;
			if(reg4){
				if(limit_reg_process_flag) {
					//sub_limit_reg_vec = getSimpleRegs(reg4->chrname, reg4->startRefPos, reg4->endRefPos, limit_reg_vec);
					sub_limit_reg_vec = getOverlappedSimpleRegsExt(reg4->chrname, reg4->startRefPos, reg4->endRefPos, limit_reg_vec, CNS_SIDE_EXT_SIZE);
					if(sub_limit_reg_vec.size()==0) flag4 = false;
				}
			}else flag4 = false;

			flag = false;
			if(flag1 or flag2 or flag3 or flag4) flag = true;

			if(flag){
				left_size = left_size2 = right_size = right_size2 = 0; var_type = VAR_UNC; dup_num_tmp = 0;
				rd_num_left = rd_num_left2 = rd_num_right = rd_num_right2 = 0;
				if(str_vec.at(13).compare("####")==0){
					if(str_vec.at(14).compare("-")!=0) { left_size = stoi(str_vec.at(14)); }
					if(str_vec.at(15).compare("-")!=0) { left_size2 = stoi(str_vec.at(15)); }
					if(str_vec.at(16).compare("-")!=0) { right_size = stoi(str_vec.at(16)); }
					if(str_vec.at(17).compare("-")!=0) { right_size2 = stoi(str_vec.at(17)); }

					if(str_vec.at(18).compare("UNC")==0) var_type = VAR_UNC;
					else if(str_vec.at(18).compare("INS")==0) var_type = VAR_INS;
					else if(str_vec.at(18).compare("DEL")==0) var_type = VAR_DEL;
					else if(str_vec.at(18).compare("DUP")==0) var_type = VAR_DUP;
					else if(str_vec.at(18).compare("INV")==0) var_type = VAR_INV;
					else if(str_vec.at(18).compare("TRA")==0) var_type = VAR_TRA;
					else if(str_vec.at(18).compare("MIX")==0) var_type = VAR_MIX;

					if(str_vec.at(19).compare("-")!=0) { dup_num_tmp = stoi(str_vec.at(19)); }

					// BND mate information
					bnd_mate_reg_str_vec = split(str_vec.at(20), ";");
					if(bnd_mate_reg_str_vec.size()!=4){
						cerr << __func__ << ", line=" << __LINE__ << ": invalid number of bnd_mate_str_vec: " << bnd_mate_reg_str_vec.size() << ", which should be 4." << endl;
						exit(1);
					}

					// support reads number information
					rd_num_left = stoi(str_vec.at(21));
					rd_num_left2 = stoi(str_vec.at(22));
					rd_num_right = stoi(str_vec.at(23));
					rd_num_right2 = stoi(str_vec.at(24));
				}

				mate_clip_reg = new mateClipReg_t();
				if(var_type==VAR_INS or var_type==VAR_DEL){ // large indel
					reg1->var_type = var_type;
					reg1->sv_len = rd_num_left;
					mate_clip_reg->leftClipReg = NULL;
					mate_clip_reg->large_indel_flag = true;
					mate_clip_reg->largeIndelClipReg = reg1;
					mate_clip_reg->supp_num_largeIndel = rd_num_left2;
					mate_clip_reg->depth_largeIndel = rd_num_right;
					mate_clip_reg->leftClipPosNum = 0;
					mate_clip_reg->leftClipPosNum2 = 0;
				}else{ // clips
					mate_clip_reg->leftClipReg = reg1;
					mate_clip_reg->large_indel_flag = false;
					mate_clip_reg->largeIndelClipReg = NULL;
					mate_clip_reg->supp_num_largeIndel = 0;
					mate_clip_reg->depth_largeIndel = 0;
					mate_clip_reg->leftClipPosNum = rd_num_left;
					mate_clip_reg->leftClipPosNum2 = rd_num_left2;
				}
				mate_clip_reg->leftClipReg2 = reg2;
				mate_clip_reg->rightClipReg = reg3;
				mate_clip_reg->rightClipReg2 = reg4;
				mate_clip_reg->reg_mated_flag = mate_flag;
				mate_clip_reg->valid_flag = true;
				mate_clip_reg->call_success_flag = false;
				mate_clip_reg->supp_num_valid_flag = true;
				mate_clip_reg->sv_type = var_type;
				mate_clip_reg->leftMeanClipPos = left_size;
				mate_clip_reg->leftMeanClipPos2 = left_size2;
				mate_clip_reg->rightMeanClipPos = right_size;
				mate_clip_reg->rightMeanClipPos2 = right_size2;

				mate_clip_reg->leftClipRegNum = 0;
				if(mate_clip_reg->leftClipReg) mate_clip_reg->leftClipRegNum ++;
				if(mate_clip_reg->leftClipReg2) mate_clip_reg->leftClipRegNum ++;
				mate_clip_reg->rightClipRegNum = 0;
				if(mate_clip_reg->rightClipReg) mate_clip_reg->rightClipRegNum ++;
				if(mate_clip_reg->rightClipReg2) mate_clip_reg->rightClipRegNum ++;

				// support reads number information
//				mate_clip_reg->leftClipPosNum = rd_num_left;
//				mate_clip_reg->leftClipPosNum2 = rd_num_left2;
				mate_clip_reg->rightClipPosNum = rd_num_right;
				mate_clip_reg->rightClipPosNum2 = rd_num_right2;

				mate_clip_reg->tra_rescue_success_flag = false;
				mate_clip_reg->chrname_leftTra1 = mate_clip_reg->chrname_rightTra1 = mate_clip_reg->chrname_leftTra2 = mate_clip_reg->chrname_rightTra2 = "";
				mate_clip_reg->leftClipPosTra1 = mate_clip_reg->rightClipPosTra1 = mate_clip_reg->leftClipPosTra2 = mate_clip_reg->rightClipPosTra2 = -1;
				mate_clip_reg->dup_num = dup_num_tmp;
				for(int i=0; i<4; i++) mate_clip_reg->bnd_mate_reg_strs[i] = bnd_mate_reg_str_vec.at(i);

				mateClipRegVector.push_back(mate_clip_reg);
			}else{
				if(reg1) delete reg1;
				if(reg2) delete reg2;
				if(reg3) delete reg3;
				if(reg4) delete reg4;
			}
		}
	}
	infile.close();

	//if(limit_reg_process_flag) delete simple_reg;
}

// compute block by reference position
Block* Chrome::computeBlocByPos(int64_t begPos, vector<Block*> &block_vec){
	int32_t bloc_ID = computeBlocID(begPos, block_vec);
	return block_vec.at(bloc_ID);
}

// compute block ID
int32_t Chrome::computeBlocID(int64_t begPos, vector<Block*> &block_vec){
	int32_t bloc_ID, bloc_ID_tmp;
	size_t i;
	Block *bloc, *mid_bloc;

	bloc_ID_tmp = begPos / (paras->blockSize - 2 * paras->slideSize);
	if(bloc_ID_tmp>=blockNum) bloc_ID_tmp = blockNum - 1;

	bloc_ID = bloc_ID_tmp;

	mid_bloc = block_vec.at(bloc_ID_tmp);
	if(begPos<=mid_bloc->startPos){
		for(i=bloc_ID_tmp; i>=0; i--){
			bloc = block_vec.at(i);
			if(bloc->startPos<=begPos){
				bloc_ID = i;
				break;
			}
		}
	}else{
		for(i=bloc_ID_tmp+1; i<(size_t)blockNum; i++){
			bloc = block_vec.at(i);
			if(bloc->startPos>begPos){
				bloc_ID = i - 1;
				break;
			}
		}
	}

	return bloc_ID;
}

// load previously consensus information
void Chrome::loadPrevConsInfo(){
	vector<simpleReg_t*> limit_reg_vec;
	//if(paras->limit_reg_process_flag) limit_reg_vec = getSimpleRegs(chrname, -1, -1, paras->limit_reg_vec);
	if(paras->limit_reg_process_flag) limit_reg_vec = getOverlappedSimpleRegs(chrname, -1, -1, paras->limit_reg_vec);
	loadPrevConsInfo2(false, paras->limit_reg_process_flag, limit_reg_vec, paras->cns_work_vec);
	loadPrevConsInfo2(true, paras->limit_reg_process_flag, paras->limit_reg_vec, paras->cns_work_vec);

	// clean previously consensued temporary folders
	string dir_prefix = "tmp_";
	cleanPrevConsTmpDir(out_dir_cns, dir_prefix);

	// open the consensus information file
	//chrSetVarCandFiles();
}

// load previously consensued information
void Chrome::loadPrevConsInfo2(bool clipReg_flag, bool limit_reg_process_flag, vector<simpleReg_t*> &limit_reg_vec, vector<cnsWork_opt*> &cns_work_vec){
	ifstream infile;
	ofstream *out_file;
	string infilename, line, done_str, old_out_dir, refseqfilename, contigfilename, readsfilename, pattern_str, limit_reg_str;
	string tmp_filename, header_line, new_line;
	vector<string> var_str, var_str1, var_str2;
	vector<string> str_vec, str_vec2, str_vec3;
	size_t i;
	int64_t item_id;
	simpleReg_t *simple_reg, *prev_simple_reg;
	vector<simpleReg_t*> sub_limit_reg_vec, prev_limit_reg_vec, prev_limit_reg_vec_tmp, pos_limit_reg_vec;
	bool flag, pos_contained_flag, prev_delete_flag, done_flag, recns_flag;

	if(clipReg_flag){
		infilename = var_cand_clipReg_filename;
		out_file = &var_cand_clipReg_file;
	}else{
		infilename = var_cand_indel_filename;
		out_file = &var_cand_indel_file;
	}

	if(isFileExist(infilename)){
		if(limit_reg_process_flag){
			if(clipReg_flag) pattern_str = CLIPREG_PATTERN;
			else pattern_str = REFSEQ_PATTERN;
			simple_reg = new simpleReg_t();
		}

		tmp_filename = infilename + "_tmp";
		rename(infilename.c_str(), tmp_filename.c_str());

		infile.open(tmp_filename);
		if(!infile.is_open()){
			cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << tmp_filename << endl;
			exit(1);
		}

		out_file->open(infilename);
		if(!out_file->is_open()){
			cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << infilename << endl;
			exit(1);
		}

		header_line = getCnsFileHeaderLine();
		*out_file << header_line << endl;

		old_out_dir = "";
		while(getline(infile, line)){
			if(line.size()>0 and line.at(0)!='#'){
				str_vec = split(line, "\t");

				// update item file name
				if(old_out_dir.size()==0) old_out_dir = getOldOutDirname(str_vec.at(0), paras->out_dir_cns);
				refseqfilename = getUpdatedItemFilename(str_vec.at(0), paras->outDir, old_out_dir);
				contigfilename = getUpdatedItemFilename(str_vec.at(1), paras->outDir, old_out_dir);
				readsfilename = getUpdatedItemFilename(str_vec.at(2), paras->outDir, old_out_dir);

	//			if(refseqfilename.compare("output_test_limit_reg_20200807/output_chr2/2_cns/chr2/refseq_chr2_149997126-150006879.fa")==0){
	//				cout << __LINE__ << ": " << refseqfilename << endl;
	//			}

				flag = true;
				pos_contained_flag = false;
				prev_delete_flag = true;
				if(limit_reg_process_flag) {
					getRegByFilename(simple_reg, refseqfilename, pattern_str);
					//sub_limit_reg_vec = getSimpleRegs(simple_reg->chrname, simple_reg->startPos, simple_reg->endPos, limit_reg_vec);
					sub_limit_reg_vec = getOverlappedSimpleRegs(simple_reg->chrname, simple_reg->startPos, simple_reg->endPos, limit_reg_vec);
					if(sub_limit_reg_vec.size()==0){ // further to check previously recorded limited process regions
						flag = false;
						prev_limit_reg_vec = extractSimpleRegsByStr(str_vec.at(8));
						for(i=0; i<prev_limit_reg_vec.size(); i++){
							prev_simple_reg = prev_limit_reg_vec.at(i);
							//prev_limit_reg_vec_tmp = getSimpleRegs(prev_simple_reg->chrname, prev_simple_reg->startPos, prev_simple_reg->endPos, limit_reg_vec);
							prev_limit_reg_vec_tmp = getFullyContainedSimpleRegs(prev_simple_reg->chrname, prev_simple_reg->startPos, prev_simple_reg->endPos, limit_reg_vec);
							if(prev_limit_reg_vec_tmp.size()){ // region fully contained
								flag = true;
								break;
							}
						}

						if(flag==false){ // position contained
							pos_limit_reg_vec = getPosContainedSimpleRegs(simple_reg->chrname, simple_reg->startPos, simple_reg->endPos, limit_reg_vec);
							if(pos_limit_reg_vec.size()){
								flag = true;
								pos_contained_flag = true;
							}
						}
					}
				}

				if(flag){
					// remove the item from consensus work vector
					done_flag = recns_flag = false;
					item_id = getItemIDFromCnsWorkVec(contigfilename, cns_work_vec);
					if(item_id!=-1){ // already consensused
						done_str = str_vec.at(str_vec.size()-1);
						if(done_str.compare(DONE_STR)==0 and isFileExist(refseqfilename)) {
							if(isReadableFile(contigfilename) and isReadableFile(refseqfilename))
								done_flag = true;
						}
					}
					if(paras->recns_failed_work_flag and str_vec.at(5).compare(CNS_FAILURE)==0) recns_flag = true;
					if(done_flag and recns_flag==false){ // consensus done and NOT re-consensus
						deleteItemFromCnsWorkVec(item_id, cns_work_vec); // already consensused, then remove it from the work vector
						paras->cns_reg_preDone_num ++;
						paras->cns_reg_work_total --;

						// save to file
						new_line = refseqfilename + "\t" + contigfilename + "\t" + readsfilename;
						for(i=3; i<=7; i++) new_line += "\t" + str_vec.at(i);

						// limit regions
						if(limit_reg_process_flag){
							if(sub_limit_reg_vec.size()) limit_reg_str = getLimitRegStr(sub_limit_reg_vec);
							else {
								if(pos_contained_flag==false) limit_reg_str = getLimitRegStr(prev_limit_reg_vec);
								else limit_reg_str = getLimitRegStr(pos_limit_reg_vec);
							}
						}else limit_reg_str = LIMIT_REG_ALL_STR;
						new_line += "\t" + limit_reg_str;

						for(i=9; i<str_vec.size(); i++) new_line += "\t" + str_vec.at(i); // other fields
						*out_file << new_line << endl;
					}
				}

				if(!prev_limit_reg_vec.empty() and prev_delete_flag) destroyLimitRegVector(prev_limit_reg_vec);
			}
		}

		infile.close();
		remove(tmp_filename.c_str());
		if(limit_reg_process_flag) delete simple_reg;
	}else{
		out_file->open(infilename);
		if(!out_file->is_open()){
			cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << infilename << endl;
			exit(1);
		}

		header_line = getCnsFileHeaderLine();
		*out_file << header_line << endl;
	}
}

// get the consensus file header line which starts with '#'
string Chrome::getCnsFileHeaderLine(){
	string header_line = "#RefFile\tContigFile\tReadsFile\tLeftShiftSize\tRightShiftSize\tStatus\tVarRegs\tCovSamplingInfo\tLimitRegs\tDoneFlag";
	return header_line;
}

// generate local consensus work for chrome
int Chrome::chrGenerateLocalConsWorkOpt(){
	Time time;

	mkdir(out_dir_cns.c_str(), S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);  // create the directory for consensus command

	//cout << "[" << time.getTime() << "]: processing Chr: " << chrname << ", size: " << chrlen << " bp" << endl;

	// generate local consensus work
	if(paras->num_threads<=1) chrGenerateLocalConsWorkOpt_st();  // single thread
	else chrGenerateLocalConsWorkOpt_mt();  // multiple threads

//	// close and reset the consensus information file
//	chrResetVarCandFiles();

	return 0;
}

// generate local consensus work for chrome using single thread
int Chrome::chrGenerateLocalConsWorkOpt_st(){
	Block *bloc;
	for(size_t i=0; i<blockVector.size(); i++){
		bloc = blockVector.at(i);
		if(bloc->process_flag)
			bloc->blockGenerateLocalConsWorkOpt();
	}
	return 0;
}

// generate local consensus work for chrome using multiple threads
int Chrome::chrGenerateLocalConsWorkOpt_mt(){
	int32_t i;
	MultiThread mt[paras->num_threads];
	for(i=0; i<paras->num_threads; i++){
		mt[i].setNumThreads(paras->num_threads);
		mt[i].setBlockVec(&blockVector);
		mt[i].setUserThreadID(i);
		if(!mt[i].startGenConsWorkOpt()){
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
	return 0;
}

// reset chrome consensus data
void Chrome::chrResetConsData(){
	// close and reset the file information
	chrResetVarCandFiles();
}

// output consensus data to file
void Chrome::outputCnsDataToFile(string &filename){
	string line, cns_status, header, left_shift_size_str, right_shift_size_str;
	varCand *item;
	vector<reg_t*> varVec;
	reg_t *reg;

	ofstream outfile(filename);
	if(!outfile.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file" << filename << endl;
		exit(1);
	}

	for(size_t i=0; i<var_cand_vec.size(); i++){
		item = var_cand_vec[i];
		if(item->cns_success) cns_status = CNS_SUCCESS;
		else cns_status = CNS_FAILURE;
		line = item->refseqfilename + "\t" + item->ctgfilename + "\t" + item->readsfilename + "\t" + to_string(item->ref_left_shift_size) + "\t" + to_string(item->ref_right_shift_size) + "\t" + cns_status;

		varVec = item->varVec;
		for(size_t j=0; j<varVec.size(); j++){
			reg = varVec[j];
			line += "\t" + reg->chrname + ":" + to_string(reg->startRefPos) + "-" + to_string(reg->endRefPos);
		}

		outfile << line << endl;
	}
	outfile.close();
}

// collect call works
void Chrome::chrCollectCallWork(){
	int32_t total_chr;
	// set minimap2 var_cand files
	//setBlatVarcandFiles();

	total_chr = var_cand_vec.size() + var_cand_clipReg_vec.size();
	if(total_chr>0){
		paras->call_work_num += total_chr;
		for(size_t i=0; i<var_cand_vec.size(); i++) paras->call_work_vec.push_back(var_cand_vec.at(i));
		for(size_t i=0; i<var_cand_clipReg_vec.size(); i++) paras->call_work_vec.push_back(var_cand_clipReg_vec.at(i));
		paras->call_work_vec.shrink_to_fit();
		//cout << "\tChr " << chrname << ": " << var_cand_vec.size() << " indel regions, " << var_cand_clipReg_vec.size() << " clipping regions, total regions: " << total_chr << endl;
	}
}


int Chrome::chrCall(){
	Time time;

//	mkdir(out_dir_call.c_str(), S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);  // create the directory for call command

	if(print_flag) cout << "[" << time.getTime() << "]: processing Chr: " << chrname << ", size: " << chrlen << " bp" << endl;

	// set blat var_cand files
	setBlatVarcandFiles();

	// call variants
	if(paras->num_threads<=1) chrCall_st();  // single thread
	else chrCall_mt();  // multiple threads

	// reset blat var_cand files
	resetBlatVarcandFiles();

	// release blat align information
	if(!blat_aligned_chr_varCand_vec.empty()) destroyVarCandVector(blat_aligned_chr_varCand_vec);
	if(!blat_aligned_chr_clipReg_varCand_vec.empty()) destroyVarCandVector(blat_aligned_chr_clipReg_varCand_vec);

//	// remove redundant new called variants
//	cout << "--[" << time.getTime() << "]: aaaaaaaaaaaaa: removeRedundantVar..." << endl;
//	removeRedundantVar();
//
//	// remove FPs of new called variants
//	cout << "--[" << time.getTime() << "]: bbbbbbbbbbbbb: removeFPNewVarVec..." << endl;
//	removeFPNewVarVec();

	return 0;
}

// call variants for chrome using single thread
void Chrome::chrCall_st(){
	chrCallVariants(var_cand_vec);
	chrCallVariants(var_cand_clipReg_vec);
}

// call variants according to variant candidates
void Chrome::chrCallVariants(vector<varCand*> &var_cand_vec){
	varCand *var_cand;
	for(size_t i=0; i<var_cand_vec.size(); i++){
		var_cand = var_cand_vec.at(i);
		//if(var_cand->alnfilename.compare("output_v0.10.1_20210511/3_call/hs37d5/blat_hs37d5_34358416-34368733.sim4")==0)
		//if(i>=148)
		{
			//cout << ">>>>>>>>> " << i << "/" << var_cand_vec.size() << ", " << var_cand->alnfilename << ", " << var_cand->ctgfilename << endl;
			var_cand->callVariants();
		}
	}
}

// call variants for chrome using multiple threads
void Chrome::chrCall_mt(){
	int32_t i;
	MultiThread mt[paras->num_threads];
	for(i=0; i<paras->num_threads; i++){
		mt[i].setNumThreads(paras->num_threads);
		mt[i].setVarCandVec(&var_cand_vec, &var_cand_clipReg_vec);
		mt[i].setUserThreadID(i);
		if(!mt[i].startCall()){
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


void Chrome::chrLoadDataCall(){
	// load misAln regions
	if(paras->maskMisAlnRegFlag and mis_aln_vec.size()==0){
		loadMisAlnRegData();
		sortMisAlnRegData();
	}

	loadVarCandData();
	if(!isVarCandDataSorted(var_cand_vec)) { sortVarCandData(var_cand_vec); /*var_cand.outputCnsDataToFile(outfilename);*/ }
	if(!isVarCandDataSorted(var_cand_clipReg_vec)) { sortVarCandData(var_cand_clipReg_vec); /*var_cand.outputCnsDataToFile(outfilename);*/ }
}

void Chrome::loadVarCandData(){
	vector<simpleReg_t*> limit_reg_vec;
	//if(paras->limit_reg_process_flag) limit_reg_vec = getSimpleRegs(chrname, -1, -1, paras->limit_reg_vec);
	if(paras->limit_reg_process_flag) limit_reg_vec = getOverlappedSimpleRegs(chrname, -1, -1, paras->limit_reg_vec);

	//if(mateClipRegVector.size()==0) chrLoadMateClipRegDataOp(paras->limit_reg_process_flag, limit_reg_vec);
	loadVarCandDataFromFile(var_cand_vec, var_cand_indel_filename, false, paras->limit_reg_process_flag, limit_reg_vec);
	loadVarCandDataFromFile(var_cand_clipReg_vec, var_cand_clipReg_filename, true, paras->limit_reg_process_flag, paras->limit_reg_vec);

	// load previously minimap2 aligned information
	//loadPrevMinimapAlnItems(false, paras->limit_reg_process_flag, limit_reg_vec);
	//loadPrevMinimapAlnItems(true, paras->limit_reg_process_flag, paras->limit_reg_vec);
}

void Chrome::chrLoadMateClipRegData(){
	mkdir(out_dir_call.c_str(), S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);  // create the directory for call command

	vector<simpleReg_t*> limit_reg_vec;
	if(mateClipRegVector.size()==0){
		//if(paras->limit_reg_process_flag) limit_reg_vec = getSimpleRegs(chrname, -1, -1, paras->limit_reg_vec);
		if(paras->limit_reg_process_flag) limit_reg_vec = getOverlappedSimpleRegs(chrname, -1, -1, paras->limit_reg_vec);
		chrLoadMateClipRegDataOp(paras->limit_reg_process_flag, limit_reg_vec);
	}
}

void Chrome::loadVarCandDataFromFile(vector<varCand*> &var_cand_vec, string &var_cand_filename, bool clipReg_flag, bool limit_reg_process_flag, vector<simpleReg_t*> &limit_reg_vec){
	string line, alnfilename, str_tmp, chrname_str, old_out_dir, refseqfilename, contigfilename, readsfilename, clusterfilename, pattern_str;
	string chrname_mate_clip_reg, dirname_call_mate_clip_reg;
	vector<string> line_vec, var_str, var_str1, var_str2;
	vector<string> str_vec, str_vec2, str_vec3;
	ifstream infile;
	size_t i, lineNum;
	reg_t *reg, *reg1, *reg2;
	varCand *var_cand_tmp;
	mateClipReg_t *mate_clip_reg;
	int32_t clip_reg_idx_tra;
	simpleReg_t *simple_reg, *prev_simple_reg;
	vector<simpleReg_t*> sub_limit_reg_vec, prev_limit_reg_vec, prev_limit_reg_vec_tmp, pos_limit_reg_vec;
	bool flag, pos_contained_flag, prev_delete_flag;

	infile.open(var_cand_filename);
	if(!infile.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file: " << var_cand_filename << endl;
		exit(1);
	}

	if(limit_reg_process_flag){
		if(clipReg_flag) pattern_str = CLIPREG_PATTERN;
		else pattern_str = REFSEQ_PATTERN;
		simple_reg = new simpleReg_t();
	}

	lineNum = 0;
	while(getline(infile, line)){
		if(line.size()>0 and line.at(0)!='#'){
			line_vec = split(line, "\t");

			// update item file name
			if(old_out_dir.size()==0) old_out_dir = getOldOutDirname(line_vec.at(0), paras->out_dir_cns);
			refseqfilename = getUpdatedItemFilename(line_vec.at(0), paras->outDir, old_out_dir);
			contigfilename = getUpdatedItemFilename(line_vec.at(1), paras->outDir, old_out_dir);
			readsfilename = getUpdatedItemFilename(line_vec.at(2), paras->outDir, old_out_dir);
			clusterfilename = getUpdatedItemFilename(line_vec.at(3), paras->outDir, old_out_dir);
			chrname_mate_clip_reg = getChrnameByFilename(contigfilename);

			//if(refseqfilename.compare("output_20210418/2_cns/chr1/clipReg_refseq_chr1_1127942-1132735.fa")==0){
			//	cout << line << endl;
			//}

			flag = true;
			pos_contained_flag = false;
			prev_delete_flag = true;
			if(limit_reg_process_flag) {
				getRegByFilename(simple_reg, refseqfilename, pattern_str);
				//sub_limit_reg_vec = getSimpleRegs(simple_reg->chrname, simple_reg->startPos, simple_reg->endPos, limit_reg_vec);
				sub_limit_reg_vec = getOverlappedSimpleRegs(simple_reg->chrname, simple_reg->startPos, simple_reg->endPos, limit_reg_vec);
				if(sub_limit_reg_vec.size()==0) {
					flag = false;
					prev_limit_reg_vec = extractSimpleRegsByStr(line_vec.at(9));
					for(i=0; i<prev_limit_reg_vec.size(); i++){
						prev_simple_reg = prev_limit_reg_vec.at(i);
						//prev_limit_reg_vec_tmp = getSimpleRegs(prev_simple_reg->chrname, prev_simple_reg->startPos, prev_simple_reg->endPos, limit_reg_vec);
						prev_limit_reg_vec_tmp = getFullyContainedSimpleRegs(prev_simple_reg->chrname, prev_simple_reg->startPos, prev_simple_reg->endPos, limit_reg_vec);
						if(prev_limit_reg_vec_tmp.size()){ // region fully contained
							flag = true;
							break;
						}
					}

					if(flag==false){ // position contained
						pos_limit_reg_vec = getPosContainedSimpleRegs(simple_reg->chrname, simple_reg->startPos, simple_reg->endPos, limit_reg_vec);
						if(pos_limit_reg_vec.size()){
							flag = true;
							pos_contained_flag = true;
						}
					}
				}
			}

			if(flag){
				// allocate memory
				var_cand_tmp = new varCand();
				var_cand_tmp->chrname = chrname;
				var_cand_tmp->var_cand_filename = var_cand_filename;
				var_cand_tmp->out_dir_call = out_dir_call;
				var_cand_tmp->misAln_filename = misAln_reg_filename;
				var_cand_tmp->inBamFile = paras->inBamFile;
				var_cand_tmp->fai = fai;
				var_cand_tmp->technology = paras->technology;

				var_cand_tmp->refseqfilename = refseqfilename;	// refseq file name
				var_cand_tmp->ctgfilename = contigfilename;	// contig file name
				var_cand_tmp->readsfilename = readsfilename;	// reads file name
				var_cand_tmp->clusterfilename = clusterfilename; // cluster file name
				var_cand_tmp->ref_left_shift_size = stoi(line_vec.at(4));	// ref_left_shift_size
				var_cand_tmp->ref_right_shift_size = stoi(line_vec.at(5));	// ref_right_shift_size

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

				var_cand_tmp->blat_aligned_info_vec = NULL;
				var_cand_tmp->blat_var_cand_file = NULL;

				var_cand_tmp->minimap2_aligned_info_vec = NULL;
				var_cand_tmp->minimap2_var_cand_file = NULL;

				if(line_vec.at(6).compare(CNS_SUCCESS)==0) var_cand_tmp->cns_success = true;
				else var_cand_tmp->cns_success = false;

				// get the number of contigs
				var_cand_tmp->ctg_num = getCtgCount(var_cand_tmp->ctgfilename);

				// load variants
				if(line_vec.at(7).compare("-")!=0){
					var_str = split(line_vec.at(7), ";");
					for(i=0; i<var_str.size(); i++){
						var_str1 = split(var_str.at(i), ":");
						var_str2 = split(var_str1.at(1), "-");
						reg = new reg_t();
						reg->chrname = var_str1.at(0);
						reg->startRefPos = stoi(var_str2.at(0));
						reg->endRefPos = stoi(var_str2.at(1));
						reg->startLocalRefPos = reg->endLocalRefPos = 0;
						reg->startQueryPos = reg->endQueryPos = 0;
						reg->sv_len = 0;
						reg->dup_num = 0;
						reg->var_type = VAR_UNC;
						reg->query_id = -1;
						reg->blat_aln_id = -1;
						reg->minimap2_aln_id = -1;
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
						var_cand_tmp->varVec.push_back(reg);  // variant vector
					}
					var_cand_tmp->varVec.shrink_to_fit();
				}

				// limit regions
				var_cand_tmp->limit_reg_process_flag = limit_reg_process_flag;
				if(sub_limit_reg_vec.size()) for(i=0; i<sub_limit_reg_vec.size(); i++) var_cand_tmp->sub_limit_reg_vec.push_back(sub_limit_reg_vec.at(i));
				else{
					if(pos_contained_flag==false){
						for(i=0; i<prev_limit_reg_vec.size(); i++) var_cand_tmp->sub_limit_reg_vec.push_back(prev_limit_reg_vec.at(i));
						prev_limit_reg_vec.clear();
						var_cand_tmp->limit_reg_delete_flag = true;
						prev_delete_flag = false;
					}else for(i=0; i<pos_limit_reg_vec.size(); i++) var_cand_tmp->sub_limit_reg_vec.push_back(pos_limit_reg_vec.at(i));
				}

				// generate alignment file names
				str_vec = split(line_vec.at(1), "/");
				str_tmp = str_vec[str_vec.size()-1];  // contig file
				str_vec2 = split(str_tmp, "_");

				dirname_call_mate_clip_reg = getDirnameCall(chrname_mate_clip_reg);
//				if(clipReg_flag) alnfilename = dirname_call_mate_clip_reg + "/blat";
//				else
					alnfilename = dirname_call_mate_clip_reg + "/minimap2";
				for(i=1; i<str_vec2.size()-1; i++)
					alnfilename += "_" + str_vec2.at(i);

				str_tmp = str_vec2.at(str_vec2.size()-1);
				str_vec3 = split(str_tmp, ".");

				alnfilename += "_" + str_vec3.at(0);
				for(i=1; i<str_vec3.size()-1; i++)
					alnfilename += "." + str_vec3.at(i);
//				if(clipReg_flag) alnfilename += ".sim4";
//				else
					alnfilename += ".paf";

				var_cand_tmp->alnfilename = alnfilename;
				var_cand_tmp->align_success = false;
				var_cand_tmp->clip_reg_flag = clipReg_flag;

				// assign clipping information
				mate_clip_reg = NULL;
				if(clipReg_flag) {
					reg1 = var_cand_tmp->varVec.at(0);
					reg2 = (var_cand_tmp->varVec.size()>=2) ? var_cand_tmp->varVec.at(var_cand_tmp->varVec.size()-1) : NULL;
					mate_clip_reg = getMateClipReg(reg1, reg2, &clip_reg_idx_tra, chrname_mate_clip_reg);
				}
				if(mate_clip_reg){
					var_cand_tmp->large_indel_flag = mate_clip_reg->large_indel_flag;
					var_cand_tmp->depth_largeIndel = mate_clip_reg->depth_largeIndel;
					if(mate_clip_reg->large_indel_flag==false){ // clip region
						if(mate_clip_reg->sv_type==VAR_DUP or mate_clip_reg->sv_type==VAR_INV){ // DUP or INV
							var_cand_tmp->leftClipRefPos = mate_clip_reg->leftMeanClipPos>0 ? mate_clip_reg->leftMeanClipPos : mate_clip_reg->leftMeanClipPos2;
							var_cand_tmp->rightClipRefPos = mate_clip_reg->rightMeanClipPos>0 ? mate_clip_reg->rightMeanClipPos : mate_clip_reg->rightMeanClipPos2;
							var_cand_tmp->sv_type = mate_clip_reg->sv_type;
							var_cand_tmp->dup_num = mate_clip_reg->dup_num;
							mate_clip_reg->var_cand = var_cand_tmp;
						}else if(mate_clip_reg->sv_type==VAR_TRA){ // TRA
							if(clip_reg_idx_tra==0){
								var_cand_tmp->leftClipRefPos = mate_clip_reg->leftMeanClipPos;
								var_cand_tmp->rightClipRefPos = mate_clip_reg->rightMeanClipPos;
								mate_clip_reg->var_cand = var_cand_tmp;
							}else if(clip_reg_idx_tra==1){
								var_cand_tmp->chrname = mate_clip_reg->leftClipReg->chrname;
								var_cand_tmp->leftClipRefPos = mate_clip_reg->leftClipReg->startRefPos;
								var_cand_tmp->rightClipRefPos = mate_clip_reg->leftClipReg->endRefPos;
								mate_clip_reg->left_var_cand_tra = var_cand_tmp;
							}else if(clip_reg_idx_tra==2){
								var_cand_tmp->chrname = mate_clip_reg->rightClipReg->chrname;
								var_cand_tmp->leftClipRefPos = mate_clip_reg->rightClipReg->startRefPos;
								var_cand_tmp->rightClipRefPos = mate_clip_reg->rightClipReg->endRefPos;
								mate_clip_reg->right_var_cand_tra = var_cand_tmp;
							}else if(clip_reg_idx_tra==3){
								var_cand_tmp->chrname = mate_clip_reg->leftClipReg->chrname;
								var_cand_tmp->leftClipRefPos = mate_clip_reg->leftMeanClipPos;
								var_cand_tmp->rightClipRefPos = mate_clip_reg->leftMeanClipPos2;
								mate_clip_reg->left_var_cand_tra = var_cand_tmp;
							}else if(clip_reg_idx_tra==4){
								var_cand_tmp->chrname = mate_clip_reg->rightClipReg->chrname;
								var_cand_tmp->leftClipRefPos = mate_clip_reg->rightMeanClipPos;
								var_cand_tmp->rightClipRefPos = mate_clip_reg->rightMeanClipPos2;
								mate_clip_reg->right_var_cand_tra = var_cand_tmp;
							}
							var_cand_tmp->sv_type = mate_clip_reg->sv_type;
							var_cand_tmp->dup_num = 0;
						}else{
							cerr << __func__ << ", line=" << __LINE__ << ", invalid variant type=" << mate_clip_reg->sv_type << ", error!" << endl;
							exit(1);
						}
					}else{ // large indel
						if(mate_clip_reg->sv_type==VAR_INS or mate_clip_reg->sv_type==VAR_DEL){ // large INS or DEL
							var_cand_tmp->leftClipRefPos = mate_clip_reg->largeIndelClipReg->startRefPos;
							var_cand_tmp->rightClipRefPos = mate_clip_reg->largeIndelClipReg->endRefPos;
							var_cand_tmp->sv_type = mate_clip_reg->sv_type;
							var_cand_tmp->dup_num = mate_clip_reg->dup_num;
							mate_clip_reg->var_cand = var_cand_tmp;
						}
					}
					for(i=0; i<4; i++) var_cand_tmp->bnd_mate_reg_strs[i] = mate_clip_reg->bnd_mate_reg_strs[i];
				}else{
					var_cand_tmp->leftClipRefPos = var_cand_tmp->rightClipRefPos = 0;
					var_cand_tmp->sv_type = VAR_UNC;
					var_cand_tmp->dup_num = 0;
					for(i=0; i<4; i++) var_cand_tmp->bnd_mate_reg_strs[i] = "-";
				}

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

				var_cand_vec.push_back(var_cand_tmp);
				lineNum ++;
			}

			if(!prev_limit_reg_vec.empty() and prev_delete_flag) destroyLimitRegVector(prev_limit_reg_vec);
		}
	}
	var_cand_vec.shrink_to_fit();
	infile.close();
	if(limit_reg_process_flag) delete simple_reg;

	//cout << var_cand_filename << "\tlineNum=" << lineNum << endl;
}

string Chrome::getDirnameCall(string &chrname_given){
	string dirname_call_ret;
	Chrome *chr;
	for(size_t i=0; i<chr_vec->size(); i++){
		chr = chr_vec->at(i);
		if(chr->chrname.compare(chrname_given)==0){
			dirname_call_ret = chr->out_dir_call;
			break;
		}
	}
	return dirname_call_ret;
}

// get the mate clip region operation
mateClipReg_t* Chrome::getMateClipReg(reg_t *reg1, reg_t *reg2, int32_t *clip_reg_idx_tra, string &chrname_mate_clip_reg){
	Chrome *chr, *chr_mate_clip_reg = NULL;
	mateClipReg_t *mate_clip_reg = NULL;

	if(chrname_mate_clip_reg.compare(chrname)==0)
		mate_clip_reg = getMateClipRegOp(reg1, reg2, clip_reg_idx_tra, mateClipRegVector);
	else{ // search other chromosomes
		for(size_t i=0; i<chr_vec->size(); i++){
			chr = chr_vec->at(i);
			if(chr->chrname.compare(chrname_mate_clip_reg)==0){
				chr_mate_clip_reg = chr;
				break;
			}
		}
		if(chr_mate_clip_reg)
			mate_clip_reg = getMateClipRegOp(reg1, reg2, clip_reg_idx_tra, chr_mate_clip_reg->mateClipRegVector);
	}

	return mate_clip_reg;
}

// get the mate clip region operation
mateClipReg_t* Chrome::getMateClipRegOp(reg_t *reg1, reg_t *reg2, int32_t *clip_reg_idx_tra, vector<mateClipReg_t*> &mate_clipReg_vec){
	mateClipReg_t *mate_clip_reg;
	size_t i;
	int32_t idx;
	reg_t *reg, *reg1_tmp, *reg2_tmp;

	idx = -1; *clip_reg_idx_tra = -1;
	for(i=0; i<mate_clipReg_vec.size(); i++){
		mate_clip_reg = mate_clipReg_vec.at(i);
		if(mate_clip_reg->valid_flag){
			if(mate_clip_reg->large_indel_flag==false){ // clip region
				if(mate_clip_reg->sv_type==VAR_INV or mate_clip_reg->sv_type==VAR_DUP){ // DUP or INV
					//if(((reg1==NULL and mate_clip_reg->leftClipReg==NULL) or (reg1 and mate_clip_reg->leftClipReg and reg1->chrname.compare(mate_clip_reg->leftClipReg->chrname)==0 and reg1->startRefPos==mate_clip_reg->leftClipReg->startRefPos and reg1->endRefPos==mate_clip_reg->leftClipReg->endRefPos))
					//	and ((reg2==NULL and mate_clip_reg->rightClipReg==NULL) or (reg2 and mate_clip_reg->rightClipReg and reg2->chrname.compare(mate_clip_reg->rightClipReg->chrname)==0 and reg2->startRefPos==mate_clip_reg->rightClipReg->startRefPos and reg2->endRefPos==mate_clip_reg->rightClipReg->endRefPos)))
					if(reg1 and reg2){
						reg1_tmp = mate_clip_reg->leftClipReg ? mate_clip_reg->leftClipReg : mate_clip_reg->leftClipReg2;
						reg2_tmp = mate_clip_reg->rightClipReg ? mate_clip_reg->rightClipReg : mate_clip_reg->rightClipReg2;
						if((reg1->chrname.compare(reg1_tmp->chrname)==0 and reg1->startRefPos==reg1_tmp->startRefPos and reg1->endRefPos==reg1_tmp->endRefPos)
							and (reg2->chrname.compare(reg2_tmp->chrname)==0 and reg2->startRefPos==reg2_tmp->startRefPos and reg2->endRefPos==reg2_tmp->endRefPos))
							idx = i;
					}
				}else if(mate_clip_reg->sv_type==VAR_TRA){ // TRA
					if(reg1 and reg2){
						if((mate_clip_reg->leftClipReg and mate_clip_reg->rightClipReg)
							and (reg1->chrname.compare(mate_clip_reg->leftClipReg->chrname)==0 and reg1->startRefPos==mate_clip_reg->leftClipReg->startRefPos and reg1->endRefPos==mate_clip_reg->leftClipReg->endRefPos)
							and (reg2->chrname.compare(mate_clip_reg->rightClipReg->chrname)==0 and reg2->startRefPos==mate_clip_reg->rightClipReg->startRefPos and reg2->endRefPos==mate_clip_reg->rightClipReg->endRefPos)){
							idx = i; *clip_reg_idx_tra = 0;
						}else if((mate_clip_reg->leftClipReg and mate_clip_reg->leftClipReg2)
							and (reg1->chrname.compare(mate_clip_reg->leftClipReg->chrname)==0 and reg1->startRefPos==mate_clip_reg->leftClipReg->startRefPos and reg1->endRefPos==mate_clip_reg->leftClipReg->endRefPos)
							and (reg2->chrname.compare(mate_clip_reg->leftClipReg2->chrname)==0 and reg2->startRefPos==mate_clip_reg->leftClipReg2->startRefPos and reg2->endRefPos==mate_clip_reg->leftClipReg2->endRefPos)){
							idx = i; *clip_reg_idx_tra = 3;
						}else if((mate_clip_reg->rightClipReg and mate_clip_reg->rightClipReg2)
							and (reg1->chrname.compare(mate_clip_reg->rightClipReg->chrname)==0 and reg1->startRefPos==mate_clip_reg->rightClipReg->startRefPos and reg1->endRefPos==mate_clip_reg->rightClipReg->endRefPos)
							and (reg2->chrname.compare(mate_clip_reg->rightClipReg2->chrname)==0 and reg2->startRefPos==mate_clip_reg->rightClipReg2->startRefPos and reg2->endRefPos==mate_clip_reg->rightClipReg2->endRefPos)){
								idx = i; *clip_reg_idx_tra = 4;
							}
					}else if(reg1 or reg2){
						reg = reg1 ? reg1 : reg2;
						if(mate_clip_reg->leftClipReg and reg->chrname.compare(mate_clip_reg->leftClipReg->chrname)==0 and reg->startRefPos==mate_clip_reg->leftClipReg->startRefPos and reg->endRefPos==mate_clip_reg->leftClipReg->endRefPos){
							idx = i; *clip_reg_idx_tra = 1;
						}else if(mate_clip_reg->rightClipReg and reg->chrname.compare(mate_clip_reg->rightClipReg->chrname)==0 and reg->startRefPos==mate_clip_reg->rightClipReg->startRefPos and reg->endRefPos==mate_clip_reg->rightClipReg->endRefPos){
							idx = i; *clip_reg_idx_tra = 2;
						}
					}
				}else{
					cerr << __func__ << ", line=" << __LINE__ << ", invalid sv_type=" << mate_clip_reg->sv_type << ", error!" << endl;
					exit(1);
				}
				if(idx!=-1) break;
			}else{ // large indel
				if(mate_clip_reg->sv_type==VAR_INS or mate_clip_reg->sv_type==VAR_DEL){ // large INS or DEL
					if((reg1 or reg2) and mate_clip_reg->largeIndelClipReg){
						reg = reg1 ? reg1 : reg2;
						if(mate_clip_reg->largeIndelClipReg->chrname.compare(reg->chrname)==0 and reg->startRefPos==mate_clip_reg->largeIndelClipReg->startRefPos and reg->endRefPos==mate_clip_reg->largeIndelClipReg->endRefPos){
							idx = i;
						}
					}
				}
			}
		}
	}

	if(idx!=-1) mate_clip_reg = mate_clipReg_vec.at(idx);
	else mate_clip_reg = NULL;

	return mate_clip_reg;
}

// sort the items in ascending order
void Chrome::sortVarCandData(vector<varCand*> &var_cand_vec){
	varCand *item;
	size_t minIdx;
	for(size_t i=0; i<var_cand_vec.size(); i++){
		minIdx = i;
		for(size_t j=i+1; j<var_cand_vec.size(); j++)
			if(var_cand_vec[j]->varVec[0]->startRefPos < var_cand_vec[minIdx]->varVec[0]->startRefPos)
				minIdx = j;

		if(minIdx!=i){
			item = var_cand_vec[i];
			var_cand_vec[i] = var_cand_vec[minIdx];
			var_cand_vec[minIdx] = item;
		}
	}
}

// determine whether the var_cand_vec were sorted
bool Chrome::isVarCandDataSorted(vector<varCand*> &var_cand_vec){
	bool flag = true;
	vector<reg_t*> varVec1, varVec2;
	for(size_t i=1; i<var_cand_vec.size(); i++){
		varVec1 = var_cand_vec[i-1]->varVec;
		varVec2 = var_cand_vec[i]->varVec;
		if(varVec1[varVec1.size()-1]->endRefPos >= varVec2[0]->startRefPos){
			flag = false;
			break;
		}
	}
	return flag;
}

// load misAln region data
void Chrome::loadMisAlnRegData(){
	ifstream infile;
	string line;
	vector<string> line_vec;
	reg_t *reg;

	if(isFileExist(misAln_reg_filename)==false) return;

	infile.open(misAln_reg_filename);
	if(!infile.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file " << misAln_reg_filename << endl;
		exit(1);
	}

	while(getline(infile, line)){
		if(line.size()){
			reg = new reg_t();

			line_vec = split(line, "\t");
			reg->chrname = line_vec[0];
			reg->startRefPos = stoi(line_vec[1]);
			reg->endRefPos = stoi(line_vec[2]);
			reg->var_type = VAR_UNC;
			reg->query_id = -1;
			reg->sv_len = 0;
			reg->blat_aln_id = -1;
			reg->minimap2_aln_id = -1;
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
			mis_aln_vec.push_back(reg);
		}
	}

	infile.close();
}

// sort misAln region data
void Chrome::sortMisAlnRegData(){
	reg_t *item;
	size_t minIdx;
	for(size_t i=0; i<mis_aln_vec.size(); i++){
		minIdx = i;
		for(size_t j=i+1; j<mis_aln_vec.size(); j++)
			if(mis_aln_vec[j]->startRefPos < mis_aln_vec[minIdx]->startRefPos)
				minIdx = j;

		if(minIdx!=i){
			item = mis_aln_vec[i];
			mis_aln_vec[i] = mis_aln_vec[minIdx];
			mis_aln_vec[minIdx] = item;
		}
	}
}

// load previously minimap2 aligned information
void Chrome::loadPrevMinimapAlnItems(bool clipReg_flag, bool limit_reg_process_flag, vector<simpleReg_t*> &limit_reg_vec){
	ifstream infile;
	string infilename, line, done_str, old_out_dir, alnfilename, refseqfilename, contigfilename, pattern_str;
	vector<string> line_vec;
	varCand *var_cand_tmp;
	simpleReg_t *simple_reg, *prev_simple_reg;
	vector<simpleReg_t*> sub_limit_reg_vec, prev_limit_reg_vec, prev_limit_reg_vec_tmp, pos_limit_reg_vec;
	bool flag, pos_contained_flag, prev_delete_flag, aln_ctg_match_flag; //, aln_completed_flag;
	size_t i;

	if(clipReg_flag) infilename = blat_var_cand_clipReg_filename;
	else infilename = blat_var_cand_indel_filename;

	if(isFileExist(infilename)==false) return;

	infile.open(infilename);
	if(!infile.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << infilename << endl;
		exit(1);
	}

	if(limit_reg_process_flag){
		if(clipReg_flag) pattern_str = CLIPREG_PATTERN;
		else pattern_str = REFSEQ_PATTERN;
		simple_reg = new simpleReg_t();
	}

	old_out_dir = "";
	while(getline(infile, line)){
		if(line.size()>0 and line.at(0)!='#'){
			line_vec = split(line, "\t");

			// update item file name
			if(old_out_dir.size()==0) old_out_dir = getOldOutDirname(line_vec.at(0), paras->out_dir_call);
			alnfilename = getUpdatedItemFilename(line_vec.at(0), paras->outDir, old_out_dir);
			contigfilename = getUpdatedItemFilename(line_vec.at(1), paras->outDir, old_out_dir);
			refseqfilename = getUpdatedItemFilename(line_vec.at(2), paras->outDir, old_out_dir);

			flag = true;
			pos_contained_flag = false;
			prev_delete_flag = true;
			if(limit_reg_process_flag) {
				getRegByFilename(simple_reg, refseqfilename, pattern_str);
				//sub_limit_reg_vec = getSimpleRegs(simple_reg->chrname, simple_reg->startPos, simple_reg->endPos, limit_reg_vec);
				sub_limit_reg_vec = getOverlappedSimpleRegs(simple_reg->chrname, simple_reg->startPos, simple_reg->endPos, limit_reg_vec);
				if(sub_limit_reg_vec.size()==0) {
					flag = false;
					prev_limit_reg_vec = extractSimpleRegsByStr(line_vec.at(4));
					for(i=0; i<prev_limit_reg_vec.size(); i++){
						prev_simple_reg = prev_limit_reg_vec.at(i);
						//prev_limit_reg_vec_tmp = getSimpleRegs(prev_simple_reg->chrname, prev_simple_reg->startPos, prev_simple_reg->endPos, limit_reg_vec);
						prev_limit_reg_vec_tmp = getFullyContainedSimpleRegs(prev_simple_reg->chrname, prev_simple_reg->startPos, prev_simple_reg->endPos, limit_reg_vec);
						if(prev_limit_reg_vec_tmp.size()){  // region fully contained
							flag = true;
							break;
						}
					}

					if(flag==false){  // position contained
						pos_limit_reg_vec = getPosContainedSimpleRegs(simple_reg->chrname, simple_reg->startPos, simple_reg->endPos, limit_reg_vec);
						if(pos_limit_reg_vec.size()){
							flag = true;
							pos_contained_flag = true;
						}
					}
				}
			}

			if(flag){
				// allocate memory
				var_cand_tmp = new varCand();

				var_cand_tmp->chrname = "";
				var_cand_tmp->var_cand_filename = "";
				var_cand_tmp->out_dir_call = "";
				var_cand_tmp->misAln_filename = "";
				var_cand_tmp->inBamFile = "";
				var_cand_tmp->fai = NULL;
				var_cand_tmp->technology = paras->technology;

				var_cand_tmp->alnfilename = alnfilename;  // align file name
				var_cand_tmp->ctgfilename = contigfilename;  // contig file name
				var_cand_tmp->refseqfilename = refseqfilename;  // refseq file name
				var_cand_tmp->ref_left_shift_size = 0;  // ref_left_shift_size
				var_cand_tmp->ref_right_shift_size = 0;  // ref_right_shift_size

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

				if(line_vec[3].compare(ALIGN_SUCCESS)==0) var_cand_tmp->align_success = true;
				else var_cand_tmp->align_success = false;

				var_cand_tmp->ctg_num = 0;

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

				// limit regions
				var_cand_tmp->limit_reg_process_flag = limit_reg_process_flag;
				if(sub_limit_reg_vec.size()>0) for(i=0; i<sub_limit_reg_vec.size(); i++) var_cand_tmp->sub_limit_reg_vec.push_back(sub_limit_reg_vec.at(i));
				else{
					if(pos_contained_flag==false){
						for(i=0; i<prev_limit_reg_vec.size(); i++) var_cand_tmp->sub_limit_reg_vec.push_back(prev_limit_reg_vec.at(i));
						var_cand_tmp->limit_reg_delete_flag = true;
						prev_delete_flag = false;
					}else for(i=0; i<pos_limit_reg_vec.size(); i++) var_cand_tmp->sub_limit_reg_vec.push_back(pos_limit_reg_vec.at(i));
				}

				// deal with 'DONE' string
				aln_ctg_match_flag = false;
				//aln_completed_flag = false;
				done_str = line_vec.at(line_vec.size()-1);
				if(done_str.compare(DONE_STR)==0) {
					aln_ctg_match_flag = isBlatAlnResultMatch(contigfilename, alnfilename);
					//aln_completed_flag = isBlatAlnCompleted(alnfilename);
				}
				//if(aln_ctg_match_flag and aln_completed_flag){
				if(aln_ctg_match_flag){
					if(clipReg_flag) blat_aligned_chr_clipReg_varCand_vec.push_back(var_cand_tmp);
					else blat_aligned_chr_varCand_vec.push_back(var_cand_tmp);
				}else delete var_cand_tmp;
			}

			if(!prev_limit_reg_vec.empty() and prev_delete_flag) destroyLimitRegVector(prev_limit_reg_vec);
		}
	}

	if(clipReg_flag) blat_aligned_chr_clipReg_varCand_vec.shrink_to_fit();
	else blat_aligned_chr_varCand_vec.shrink_to_fit();

	infile.close();
	if(limit_reg_process_flag) delete simple_reg;
}

// set blat var_cand files
void Chrome::setBlatVarcandFiles(){
	ifstream infile;
	string line, new_line, tmp_filename, header_line, old_out_dir, alnfilename, contigfilename, refseqfilename, pattern_str, limit_reg_str;
	vector<string> str_vec;
	varCand *var_cand;
	size_t i;
	simpleReg_t *simple_reg, *prev_simple_reg;
	vector<simpleReg_t*> sub_limit_reg_vec, prev_limit_reg_vec, prev_limit_reg_vec_tmp, pos_limit_reg_vec;
	bool flag, pos_contained_flag;

	// indel
	if(isFileExist(blat_var_cand_indel_filename)){

		if(paras->limit_reg_process_flag){
			pattern_str = REFSEQ_PATTERN;
			simple_reg = new simpleReg_t();
		}

		tmp_filename = blat_var_cand_indel_filename + "_tmp";
		rename(blat_var_cand_indel_filename.c_str(), tmp_filename.c_str());

		infile.open(tmp_filename);
		if(!infile.is_open()){
			cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << tmp_filename << endl;
			exit(1);
		}

		blat_var_cand_indel_file.open(blat_var_cand_indel_filename);
		if(!blat_var_cand_indel_file.is_open()){
			cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << blat_var_cand_indel_filename << endl;
			exit(1);
		}

		header_line = getBlatVarcandFileHeaderLine();
		blat_var_cand_indel_file << header_line << endl;

		old_out_dir = "";
		while(getline(infile, line)){
			if(line.size()>0 and line.at(0)!='#'){
				str_vec = split(line, "\t");
				if(str_vec.at(str_vec.size()-1).compare(DONE_STR)==0){
					if(old_out_dir.size()==0) old_out_dir = getOldOutDirname(str_vec.at(0), paras->out_dir_call);

					alnfilename = getUpdatedItemFilename(str_vec.at(0), paras->outDir, old_out_dir);
					contigfilename = getUpdatedItemFilename(str_vec.at(1), paras->outDir, old_out_dir);
					refseqfilename = getUpdatedItemFilename(str_vec.at(2), paras->outDir, old_out_dir);

					flag = true;
					pos_contained_flag = false;
					if(paras->limit_reg_process_flag) {
						getRegByFilename(simple_reg, refseqfilename, pattern_str);
						//sub_limit_reg_vec = getSimpleRegs(simple_reg->chrname, simple_reg->startPos, simple_reg->endPos, paras->limit_reg_vec);
						sub_limit_reg_vec = getOverlappedSimpleRegs(simple_reg->chrname, simple_reg->startPos, simple_reg->endPos, paras->limit_reg_vec);
						if(sub_limit_reg_vec.size()==0) {
							flag = false;
							prev_limit_reg_vec = extractSimpleRegsByStr(str_vec.at(4));
							for(i=0; i<prev_limit_reg_vec.size(); i++){
								prev_simple_reg = prev_limit_reg_vec.at(i);
								//prev_limit_reg_vec_tmp = getSimpleRegs(prev_simple_reg->chrname, prev_simple_reg->startPos, prev_simple_reg->endPos, paras->limit_reg_vec);
								prev_limit_reg_vec_tmp = getFullyContainedSimpleRegs(prev_simple_reg->chrname, prev_simple_reg->startPos, prev_simple_reg->endPos, paras->limit_reg_vec);
								if(prev_limit_reg_vec_tmp.size()){  // region fully contained
									flag = true;
									break;
								}
							}

							if(flag==false){  // position contained
								pos_limit_reg_vec = getPosContainedSimpleRegs(simple_reg->chrname, simple_reg->startPos, simple_reg->endPos, paras->limit_reg_vec);
								if(pos_limit_reg_vec.size()){
									flag = true;
									pos_contained_flag = true;
								}
							}
						}
					}

					if(flag){
						new_line = alnfilename + "\t" + contigfilename + "\t" + refseqfilename;
						new_line += "\t" + str_vec.at(3);

						// limit regions
						if(paras->limit_reg_process_flag){
							if(sub_limit_reg_vec.size()) limit_reg_str = getLimitRegStr(sub_limit_reg_vec);
							else {
								if(pos_contained_flag==false) limit_reg_str = getLimitRegStr(prev_limit_reg_vec);
								else limit_reg_str = getLimitRegStr(pos_limit_reg_vec);
							}
						}else limit_reg_str = LIMIT_REG_ALL_STR;
						new_line += "\t" + limit_reg_str;

						// other fields
						for(i=5; i<str_vec.size(); i++) new_line += "\t" + str_vec.at(i);

						blat_var_cand_indel_file << new_line << endl;
					}

					if(!prev_limit_reg_vec.empty()) destroyLimitRegVector(prev_limit_reg_vec);
				}else
					cout << "line=" << __LINE__ << "," << var_cand_indel_filename << ": line does not end with 'DONE'!" << endl;
			}
		}

		infile.close();
		remove(tmp_filename.c_str());
		if(paras->limit_reg_process_flag) delete simple_reg;
	}else{
		blat_var_cand_indel_file.open(blat_var_cand_indel_filename);
		if(!blat_var_cand_indel_file.is_open()){
			cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << blat_var_cand_indel_filename << endl;
			exit(1);
		}
		header_line = getBlatVarcandFileHeaderLine();
		blat_var_cand_indel_file << header_line << endl;
	}

	for(i=0; i<var_cand_vec.size(); i++){
		var_cand = var_cand_vec.at(i);
		var_cand->setBlatVarcandFile(&blat_var_cand_indel_file, &blat_aligned_chr_varCand_vec);
	}

	// clipReg
	if(isFileExist(blat_var_cand_clipReg_filename)){

		if(paras->limit_reg_process_flag){
			pattern_str = CLIPREG_PATTERN;
			simple_reg = new simpleReg_t();
		}

		tmp_filename = blat_var_cand_clipReg_filename + "_tmp";
		rename(blat_var_cand_clipReg_filename.c_str(), tmp_filename.c_str());

		infile.open(tmp_filename);
		if(!infile.is_open()){
			cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << tmp_filename << endl;
			exit(1);
		}

		blat_var_cand_clipReg_file.open(blat_var_cand_clipReg_filename);
		if(!blat_var_cand_clipReg_file.is_open()){
			cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << blat_var_cand_clipReg_filename << endl;
			exit(1);
		}

		header_line = getBlatVarcandFileHeaderLine();
		blat_var_cand_clipReg_file << header_line << endl;

		old_out_dir = "";
		while(getline(infile, line)){
			if(line.size()>0 and line.at(0)!='#'){
				str_vec = split(line, "\t");
				if(str_vec.at(str_vec.size()-1).compare(DONE_STR)==0){
					if(old_out_dir.size()==0) old_out_dir = getOldOutDirname(str_vec.at(0), paras->out_dir_call);

					alnfilename = getUpdatedItemFilename(str_vec.at(0), paras->outDir, old_out_dir);
					contigfilename = getUpdatedItemFilename(str_vec.at(1), paras->outDir, old_out_dir);
					refseqfilename = getUpdatedItemFilename(str_vec.at(2), paras->outDir, old_out_dir);

					flag = true;
					pos_contained_flag = false;
					if(paras->limit_reg_process_flag) {
						getRegByFilename(simple_reg, refseqfilename, pattern_str);
						//sub_limit_reg_vec = getSimpleRegs(simple_reg->chrname, simple_reg->startPos, simple_reg->endPos, paras->limit_reg_vec);
						sub_limit_reg_vec = getOverlappedSimpleRegs(simple_reg->chrname, simple_reg->startPos, simple_reg->endPos, paras->limit_reg_vec);
						if(sub_limit_reg_vec.size()==0){
							flag = false;
							prev_limit_reg_vec = extractSimpleRegsByStr(str_vec.at(4));
							for(i=0; i<prev_limit_reg_vec.size(); i++){
								prev_simple_reg = prev_limit_reg_vec.at(i);
								//prev_limit_reg_vec_tmp = getSimpleRegs(prev_simple_reg->chrname, prev_simple_reg->startPos, prev_simple_reg->endPos, paras->limit_reg_vec);
								prev_limit_reg_vec_tmp = getFullyContainedSimpleRegs(prev_simple_reg->chrname, prev_simple_reg->startPos, prev_simple_reg->endPos, paras->limit_reg_vec);
								if(prev_limit_reg_vec_tmp.size()){ // region fully contained
									flag = true;
									break;
								}
							}

							if(flag==false){ // position contained
								pos_limit_reg_vec = getPosContainedSimpleRegs(simple_reg->chrname, simple_reg->startPos, simple_reg->endPos, paras->limit_reg_vec);
								if(pos_limit_reg_vec.size()){
									flag = true;
									pos_contained_flag = true;
								}
							}
						}
					}

					if(flag){
						new_line = alnfilename + "\t" + contigfilename + "\t" + refseqfilename;
						new_line += "\t" + str_vec.at(3);

						// limit regions
						if(paras->limit_reg_process_flag){
							if(sub_limit_reg_vec.size()) limit_reg_str = getLimitRegStr(sub_limit_reg_vec);
							else{
								if(pos_contained_flag==false) limit_reg_str = getLimitRegStr(prev_limit_reg_vec);
								else limit_reg_str = getLimitRegStr(pos_limit_reg_vec);
							}
						}else limit_reg_str = LIMIT_REG_ALL_STR;
						new_line += "\t" + limit_reg_str;

						// other fields
						for(i=5; i<str_vec.size(); i++) new_line += "\t" + str_vec.at(i);

						blat_var_cand_clipReg_file << new_line << endl;
					}

					if(!prev_limit_reg_vec.empty()) destroyLimitRegVector(prev_limit_reg_vec);
				}else
					cout << "line=" << __LINE__ << "," << var_cand_indel_filename << ": line does not end with 'DONE'!" << endl;
			}
		}

		infile.close();
		remove(tmp_filename.c_str());
		if(paras->limit_reg_process_flag) delete simple_reg;
	}else{
		blat_var_cand_clipReg_file.open(blat_var_cand_clipReg_filename);
		if(!blat_var_cand_clipReg_file.is_open()){
			cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << blat_var_cand_clipReg_filename << endl;
			exit(1);
		}
		header_line = getBlatVarcandFileHeaderLine();
		blat_var_cand_clipReg_file << header_line << endl;
	}

	for(i=0; i<var_cand_clipReg_vec.size(); i++){
		var_cand = var_cand_clipReg_vec.at(i);
		var_cand->setBlatVarcandFile(&blat_var_cand_clipReg_file, &blat_aligned_chr_clipReg_varCand_vec);
	}
}

// reset blat var_cand files
void Chrome::resetBlatVarcandFiles(){
	size_t i;
	varCand *var_cand;
	for(i=0; i<var_cand_vec.size(); i++){
		var_cand = var_cand_vec.at(i);
		var_cand->resetBlatVarcandFile();
	}
	for(i=0; i<var_cand_clipReg_vec.size(); i++){
		var_cand = var_cand_clipReg_vec.at(i);
		var_cand->resetBlatVarcandFile();
	}
}

// get blat align file header line which starts with '#'
string Chrome::getBlatVarcandFileHeaderLine(){
	string header_line = "#BlatAlnFile\tContigFile\tRefFile\tStatus\tLimitRegs\tDoneFlag";
	return header_line;
}

void Chrome::removeVarCandNodeIndel(varCand *var_cand){
	removeVarCandNode(var_cand, var_cand_vec);
}

void Chrome::removeVarCandNodeClipReg(varCand *var_cand){
	removeVarCandNode(var_cand, var_cand_clipReg_vec);
}

// fill the sequence, including reference sequence and contig sequence
void Chrome::chrFillVarseq(){
	//chrFillVarseqSingleVec(var_cand_vec);
	chrFillVarseqSingleVec(var_cand_clipReg_vec);
}

// fill the sequence, including reference sequence and contig sequence
void Chrome::chrFillVarseqSingleVec(vector<varCand*> &var_cand_vec){
	varCand *var_cand;
	for(size_t i=0; i<var_cand_vec.size(); i++){
		var_cand = var_cand_vec[i];
		//if(var_cand->ctgfilename.compare("output_hg19_v1.7_2.0_M0_20200909/2_cns/chr9_gl000199_random/contig_chr9_gl000199_random_64601-72429.fa")==0)
		{
			//cout << ">>>>>>>>> " << i << ", " << var_cand->alnfilename << ", " << var_cand->ctgfilename << endl;
			var_cand->fillVarseq();
		}
	}
}

// remove redundant called variants
void Chrome::removeRedundantVar(){
	removeRedundantIndel(var_cand_vec);
	//removeRedundantClipReg(var_cand_clipReg_vec);
}

// remove redundant called indels
void Chrome::removeRedundantIndel(vector<varCand*> &var_cand_vec){
	size_t i, j;
	varCand *var_cand;
	reg_t *reg;

	for(i=0; i<var_cand_vec.size(); i++){
		var_cand = var_cand_vec[i];
		j = 0;
		while(j<var_cand->newVarVec.size()){
			reg = var_cand->newVarVec[j];
			if(isRedundantVarItemBinSearch(reg, var_cand_vec)){ // remove item
				var_cand->newVarVec.erase(var_cand->newVarVec.begin()+j);
				delete reg;
			}else j++;
		}
	}

	// check newVarVector
	for(i=0; i<var_cand_vec.size(); i++){
		var_cand = var_cand_vec[i];
		for(j=0; j<var_cand->newVarVec.size(); j++){
			reg = var_cand->newVarVec[j];
			removeRedundantVarItemsInNewCalledVarvec(reg, i, var_cand_vec);
		}
	}
}

// determine whether the variant is redundant by checking varVector using binary search
bool Chrome::isRedundantVarItemBinSearch(reg_t *reg, vector<varCand*> &var_cand_vec){
	reg_t *foundReg;
	int32_t low, high, mid;

	low = 0;
	high = var_cand_vec.size() - 1;
	while(low<=high){
		mid = (low + high) / 2;
		foundReg = findVarvecItem(reg->startRefPos, reg->endRefPos, var_cand_vec[mid]->varVec);
		if(foundReg) return true;
		if(reg->startRefPos>var_cand_vec[mid]->varVec[var_cand_vec[mid]->varVec.size()-1]->endRefPos){
			low = mid + 1;
		}else if(reg->endRefPos<var_cand_vec[mid]->varVec[0]->startRefPos){
			high = mid - 1;
		}else break;
	}

	return false;
}

// determine whether the variant is redundant by checking newVarVector
void Chrome::removeRedundantVarItemsInNewCalledVarvec(reg_t *reg_given, int32_t idx_given, vector<varCand*> &var_cand_vec){
	varCand *var_cand;
	reg_t *foundReg;
	int32_t reg_idx;

	for(size_t i=0; i<var_cand_vec.size(); i++){
		if(i!=(size_t)idx_given){
			var_cand = var_cand_vec[i];
			foundReg = findVarvecItem(reg_given->startRefPos, reg_given->endRefPos, var_cand_vec[i]->newVarVec);
			if(foundReg){
				reg_idx = getVectorIdx(foundReg, var_cand_vec[i]->newVarVec);
				if(reg_idx!=-1){
					var_cand->newVarVec.erase(var_cand->newVarVec.begin()+reg_idx);
					delete foundReg;
				}else{
					cerr << "line=" << __LINE__ << ", error reg_idx=-1!" << endl;
					exit(1);
				}
			}
		}
	}
}

void Chrome::removeFPNewVarVec(){
	removeFPNewVarVecIndel(var_cand_vec);
}

// remove FPs of new called variants
void Chrome::removeFPNewVarVecIndel(vector<varCand*> &var_cand_vec){
	int32_t i, j, startRefPos, endRefPos;
	varCand *var_cand;
	reg_t *reg, *foundReg;
	bool FP_flag;
	vector<int32_t> numVec;

	//chrlen = faidx_seq_len(fai, chrname.c_str()); // get the reference length

	for(i=0; i<(int32_t)var_cand_vec.size(); i++){
		var_cand = var_cand_vec[i];
//		cout << i << ": " << var_cand->alnfilename << endl;
		j = 0;
		while(j<(int32_t)var_cand->newVarVec.size()){
			FP_flag = false;
			reg = var_cand->newVarVec[j];

			// debug
			//if(reg->startRefPos==41770036){
			//	cout << reg->startRefPos << endl;
			//}

			foundReg = findVarvecItem(reg->startRefPos, reg->endRefPos, mis_aln_vec);
			if(foundReg){
				FP_flag = true;
			}else{
				foundReg = findVarvecItemExtSize(reg->startRefPos, reg->endRefPos, var_cand->varVec, 5, 5);
				if(foundReg){
					FP_flag = true;
					//cout << reg->startRefPos << "-" << reg->endRefPos << endl;
				}else{
					startRefPos = reg->startRefPos - 10;
					endRefPos = reg->endRefPos + 10;
					if(startRefPos<1) startRefPos = 1;
					if(endRefPos>chrlen) endRefPos = chrlen;

					numVec = var_cand->computeDisagreeNumAndHighIndelBaseNum(reg->chrname, startRefPos, endRefPos, paras->inBamFile, fai);

					if(numVec[0]==0 and numVec[1]==0)
						FP_flag = true;
				}
			}

			// remove FP item
			if(FP_flag){
				var_cand->newVarVec.erase(var_cand->newVarVec.begin()+j);
				delete reg;
			}else j++;
		}
	}
}

void Chrome::saveCallSV2File(){
	saveCallIndelClipReg2File(out_filename_call_indel, out_filename_call_clipReg);
}

void Chrome::saveCallSV2File02(){
	saveCallIndelClipReg2File02(out_filename_call_indel, out_filename_call_clipReg);
}

void Chrome::saveCallIndelClipReg2File(string &outfilename_indel, string &outfilename_clipReg){
	size_t i, j, file_id;
	ofstream outfile_indel, outfile_clipReg;
	varCand *var_cand;
	reg_t *reg;
	string line, sv_type, header_line_bed;
	int32_t ref_dist, query_dist;
	bool size_satisfied;

	outfile_indel.open(outfilename_indel);
	if(!outfile_indel.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << outfilename_indel << endl;
		exit(1);
	}
	outfile_clipReg.open(outfilename_clipReg);
	if(!outfile_clipReg.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << outfilename_clipReg << endl;
		exit(1);
	}

	// header line
	header_line_bed = getCallFileHeaderBed(paras->sample);
	outfile_indel << header_line_bed << endl;
	outfile_clipReg << header_line_bed << endl;

	for(i=0; i<var_cand_vec.size(); i++){
		var_cand = var_cand_vec[i];
		for(j=0; j<var_cand->varVec.size(); j++){
			reg = var_cand->varVec[j];

			// choose the size-selected variants
			ref_dist = reg->endRefPos - reg->startRefPos + 1;
			query_dist = reg->endQueryPos - reg->startQueryPos + 1;
			size_satisfied = isSizeSatisfied(ref_dist, query_dist, paras->min_sv_size_usr, paras->max_sv_size_usr);

			if(reg->var_type!=VAR_UNC and reg->call_success_status and size_satisfied){
				file_id = 0;
				switch(reg->var_type){
					case VAR_UNC: sv_type = VAR_UNC_STR; break;
					case VAR_INS: sv_type = VAR_INS_STR; break;
					case VAR_DEL: sv_type = VAR_DEL_STR; break;
					case VAR_DUP: sv_type = VAR_DUP_STR; file_id = 1; break;
					case VAR_INV: sv_type = VAR_INV_STR; file_id = 1; break;
					case VAR_TRA: sv_type = VAR_TRA_STR; file_id = 1; break;
					default: sv_type = VAR_MIX_STR; break;
				}
				line = reg->chrname + "\t" + to_string(reg->startRefPos) + "\t" + to_string(reg->endRefPos) + "\t" + sv_type;
				if(reg->var_type!=VAR_TRA)
					line += "\t" + to_string(reg->sv_len);
				else
					line += "\t-";
				if(reg->var_type==VAR_DUP)
					line += "\t" + to_string(reg->dup_num);
				else
					line += "\t-";

				if(reg->refseq.size()) line += "\t" + reg->refseq;
				else line += "\t-";
				if(reg->altseq.size()) line += "\t" + reg->altseq;
				else line += "\t-";

				if(reg->gt_seq.size()==0) reg->gt_seq = GT_STR_DEFAULT;
				line += "\t" + reg->gt_seq;


				if(reg->short_sv_flag) line += "\tShortSV";

				//cout << "line=" << __LINE__ << ": " << line << endl;

				if(file_id==0) outfile_indel << line << endl;
				else outfile_clipReg << line << endl;
			}
		}
	}

	for(i=0; i<var_cand_clipReg_vec.size(); i++){
		var_cand = var_cand_clipReg_vec[i];

//		if(var_cand->alnfilename.compare("output_ccs_v1.0.1_20210528/3_call/1/blat_contig_1_1016017-1016409.sim4")==0){
//			cout << "line=" << __LINE__ << ", alnfilename=" << var_cand->alnfilename << endl;
//		}

		if(var_cand->clip_reg_flag==false){ // indel
			for(j=0; j<var_cand->varVec.size(); j++){
				reg = var_cand->varVec[j];

				// choose the size-selected variants
				ref_dist = reg->endRefPos - reg->startRefPos + 1;
				query_dist = reg->endQueryPos - reg->startQueryPos + 1;
				size_satisfied = isSizeSatisfied(ref_dist, query_dist, paras->min_sv_size_usr, paras->max_sv_size_usr);

				if(reg->var_type!=VAR_UNC and reg->call_success_status and size_satisfied){
					file_id = 0;
					switch(reg->var_type){
						case VAR_UNC: sv_type = VAR_UNC_STR; break;
						case VAR_INS: sv_type = VAR_INS_STR; break;
						case VAR_DEL: sv_type = VAR_DEL_STR; break;
						case VAR_DUP: sv_type = VAR_DUP_STR; file_id = 1; break;
						case VAR_INV: sv_type = VAR_INV_STR; file_id = 1; break;
						case VAR_TRA: sv_type = VAR_TRA_STR; file_id = 1; break;
						default: sv_type = VAR_MIX_STR; break;
					}
					line = reg->chrname + "\t" + to_string(reg->startRefPos) + "\t" + to_string(reg->endRefPos) + "\t" + sv_type;
					if(reg->var_type!=VAR_TRA)
						line += "\t" + to_string(reg->sv_len);
					else
						line += "\t-";
					if(reg->var_type==VAR_DUP)
						line += "\t" + to_string(reg->dup_num);
					else
						line += "\t-";
					if(reg->refseq.size()) line += "\t" + reg->refseq;
					else line += "\t-";
					if(reg->altseq.size()) line += "\t" + reg->altseq;
					else line += "\t-";

					if(reg->gt_seq.size()==0) reg->gt_seq = GT_STR_DEFAULT;
					line += "\t" + reg->gt_seq;

					if(reg->short_sv_flag) line += "\tShortSV";

					if(reg->var_type==VAR_UNC){
						cout << "line=" << __LINE__ << ": " << line << endl << endl << endl;
					}

					//cout << "line=" << __LINE__ << ": " << line << endl;

					if(file_id==0) outfile_indel << line << endl;
					else outfile_clipReg << line << endl;
				}
			}
		}else{ // cliping region
			reg = var_cand->clip_reg;

			// choose the size-selected variants
			if(reg){
				ref_dist = reg->endRefPos - reg->startRefPos + 1;
				query_dist = reg->endQueryPos - reg->startQueryPos + 1;
				size_satisfied = isSizeSatisfied(ref_dist, query_dist, paras->min_sv_size_usr, paras->max_sv_size_usr);
			}

			if(reg and reg->var_type!=VAR_UNC and var_cand->call_success and size_satisfied){
				file_id = 1;
				switch(reg->var_type){
					case VAR_UNC: sv_type = VAR_UNC_STR; break;
					case VAR_INS: sv_type = VAR_INS_STR; file_id = 0; break;
					case VAR_DEL: sv_type = VAR_DEL_STR; file_id = 0; break;
					case VAR_DUP: sv_type = VAR_DUP_STR; break;
					case VAR_INV: sv_type = VAR_INV_STR; break;
					case VAR_TRA: sv_type = VAR_TRA_STR; break;
					default: sv_type = VAR_MIX_STR; break;
				}

				line = reg->chrname + "\t" + to_string(reg->startRefPos) + "\t" + to_string(reg->endRefPos) + "\t" + sv_type;
				if(reg->var_type!=VAR_TRA)
					line += "\t" + to_string(reg->sv_len);
				else
					line += "\t-";
				if(reg->var_type==VAR_DUP)
					line += "\t" + to_string(reg->dup_num);
				else
					line += "\t-";
				if(reg->refseq.size()) line += "\t" + reg->refseq;
				else line += "\t-";
				if(reg->altseq.size()) line += "\t" + reg->altseq;
				else line += "\t-";
				//line += "\t" + reg->refseq + "\t" + reg->altseq;

				if(reg->gt_seq.size()==0) reg->gt_seq = GT_STR_DEFAULT;
				line += "\t" + reg->gt_seq;

				if(reg->short_sv_flag) line += "\tShortSV";

				if(reg->var_type==VAR_UNC){
					cout << "line=" << __LINE__ << ": " << line << endl << endl << endl;
				}

				//cout << "line=" << __LINE__ << ": " << line << endl;

				if(file_id==1) outfile_clipReg << line << endl;
				else outfile_indel << line << endl;
			}else{ // not confirmed SV
//				if(var_cand->sv_type!=VAR_UNC){
//					file_id = 1;
//					switch(var_cand->sv_type){
//						case VAR_UNC: sv_type = "UNCERTAIN"; break;
//						case VAR_INS: sv_type = "INS"; file_id = 0; break;
//						case VAR_DEL: sv_type = "DEL"; file_id = 0; break;
//						case VAR_DUP: sv_type = "DUP"; break;
//						case VAR_INV: sv_type = "INV"; break;
//						case VAR_TRA: sv_type = "TRA"; break;
//						default: sv_type = "MIX"; break;
//					}
//					line = var_cand->chrname;
//					if(var_cand->leftClipRefPos>0 and var_cand->rightClipRefPos>0) line += "\t" + to_string(var_cand->leftClipRefPos) + "\t" + to_string(var_cand->rightClipRefPos);
//					else line += "\t-\t-";
//					line += "\t.\t.\t" + sv_type + "\t.";
//					if(file_id==1) outfile_clipReg << line << endl;
//					else outfile_indel << line << endl;
//				}
			}
		}
	}

	outfile_indel.close();
	outfile_clipReg.close();
}

void Chrome::saveCallIndelClipReg2File02(string &outfilename_indel, string &outfilename_clipReg){
	size_t i, j;
	ofstream outfile;
	varCand *var_cand; //, *var_cand_pre, *var_cand_pre_pre;
	reg_t *reg;
	string line, sv_type, header_line_bed, discov_level_str;
	int32_t ref_dist, query_dist;
	static int32_t ins_num = 0, del_num = 0, dup_num= 0, inv_num = 0, tra_num = 0, sv_num; 
	bool size_satisfied; //, no_existed;

	outfile.open(outfilename_indel);
	if(!outfile.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << outfilename_indel << endl;
		exit(1);
	}

	// header line
	header_line_bed = getCallFileHeaderBed(paras->sample);
	outfile << header_line_bed << endl;

	for(i=0; i<var_cand_vec.size(); i++){ // indel region
		var_cand = var_cand_vec.at(i);
//		if(i>0) var_cand_pre = var_cand_vec.at(i-1);
//		if(i>1) var_cand_pre_pre = var_cand_vec.at(i-2);
		for(j=0; j<var_cand->newVarVec.size(); j++){
			//no_existed = true;
			reg = var_cand->newVarVec.at(j);
			// choose the size-selected variants
			//size_satisfied = isSizeSatisfied2(reg->sv_len, paras->min_sv_size_usr, paras->max_sv_size_usr); // deleted on 2025-03-20
			size_satisfied = isSizeSatisfied2(reg->sv_len, paras->min_sv_size_usr_final, paras->max_sv_size_usr);

//			if(i>0) no_existed = isNotAlreadyExists(var_cand_pre->newVarVec, reg);
//			if(i>1 and no_existed) no_existed = isNotAlreadyExists(var_cand_pre_pre->newVarVec, reg);
//			if(no_existed) no_existed = isNotAlreadyExists(var_cand->newVarVec, reg, j);

			//if(reg->var_type!=VAR_UNC and reg->call_success_status and size_satisfied and no_existed){ // deleted on 2024-08-02
			sv_num = 0;
			if(reg->var_type!=VAR_UNC and reg->call_success_status and size_satisfied){
				switch(reg->var_type){
					case VAR_UNC: sv_type = VAR_UNC_STR; break;
					case VAR_INS: sv_type = VAR_INS_STR; ins_num++; sv_num = ins_num; break;
					case VAR_DEL: sv_type = VAR_DEL_STR; del_num++; sv_num = del_num; break;
					case VAR_DUP: sv_type = VAR_DUP_STR; dup_num++; sv_num = dup_num; break;
					case VAR_INV: sv_type = VAR_INV_STR; inv_num++; sv_num = inv_num; break;
					case VAR_TRA: sv_type = VAR_TRA_STR; tra_num++; sv_num = tra_num; break;
					default: sv_type = VAR_MIX_STR; break;
				}
				// line = reg->chrname + "\t" + to_string(reg->startRefPos) + "\t" + to_string(reg->endRefPos) + "\t" + sv_type;
				line = reg->chrname + "\t" + to_string(reg->startRefPos) + "\t" + to_string(reg->endRefPos) + "\t" + "ASVCLR." + sv_type + "." + to_string(sv_num) + "\t" + sv_type;
				if(reg->var_type!=VAR_TRA)
					line += "\t" + to_string(reg->sv_len);
				else
					line += "\t-";
				if(reg->var_type==VAR_DUP)
					line += "\t" + to_string(reg->dup_num);
				else
					line += "\t-";

				//upperSeq
				if(reg->refseq.size()) {
					upperSeq(reg->refseq);
					line += "\t" + reg->refseq;
				}
				else line += "\t-";
				if(reg->altseq.size()){
					upperSeq(reg->altseq);
					line += "\t" + reg->altseq;
				}
				else line += "\t-";

				// extra information, split with ";"
				stringstream ss;
				ss << setprecision(3) << reg->AF;
				line += "\tDP=" + to_string(reg->DP) + ";AF=" + ss.str();

				discov_level_str = "";
				switch(reg->discover_level){
					case VAR_DISCOV_L_UNUSED:
						break;
					case VAR_DISCOV_L_CNS_ALN:
						discov_level_str = VAR_DISCOV_L_CNS_ALN_STR;
						break;
					case VAR_DISCOV_L_RESCUE_CNS_ALN:
						discov_level_str = VAR_DISCOV_L_RESCUE_CNS_ALN_STR;
						break;
					case VAR_DISCOV_L_READS:
						discov_level_str = VAR_DISCOV_L_READS_STR;
						break;
				}
				if(discov_level_str.size()>0) line = line + ";" + VAR_DISCOV_L_TITLE_STR + "=" + discov_level_str;

				if(reg->gt_seq.size()==0) reg->gt_seq = GT_STR_DEFAULT;
				line += "\t" + reg->gt_seq;


				if(reg->short_sv_flag) line += "\tShortSV";

				//cout << "line=" << __LINE__ << ": " << line << endl;

				outfile << line << endl;
			}
		}
	}
	outfile.close();

	outfile.open(outfilename_clipReg);
	if(!outfile.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << outfilename_clipReg << endl;
		exit(1);
	}
	outfile << header_line_bed << endl;

	for(i=0; i<var_cand_clipReg_vec.size(); i++){ // cliping region
		var_cand = var_cand_clipReg_vec.at(i);

//		if(var_cand->alnfilename.compare("output_ccs_v1.0.1_20210528/3_call/1/blat_contig_1_1016017-1016409.sim4")==0){
//			cout << "line=" << __LINE__ << ", alnfilename=" << var_cand->alnfilename << endl;
//		}

		for(j=0; j<2; j++){
			if(j==0) reg = var_cand->clip_reg;
			else reg = var_cand->clip_reg_allele;

			// choose the size-selected variants
			if(reg){
				ref_dist = reg->endRefPos - reg->startRefPos + 1;
				query_dist = reg->endQueryPos - reg->startQueryPos + 1;
				//size_satisfied = isSizeSatisfied(ref_dist, query_dist, paras->min_sv_size_usr, paras->max_sv_size_usr); // deleted on 2025-03-20
				size_satisfied = isSizeSatisfied(ref_dist, query_dist, paras->min_sv_size_usr_final, paras->max_sv_size_usr);
			}

			if(reg and reg->var_type!=VAR_UNC and var_cand->call_success and size_satisfied){
				switch(reg->var_type){
					case VAR_UNC: sv_type = VAR_UNC_STR; break;
					case VAR_INS: sv_type = VAR_INS_STR; ins_num++; sv_num = ins_num; break;
					case VAR_DEL: sv_type = VAR_DEL_STR; del_num++; sv_num = del_num; break;
					case VAR_DUP: sv_type = VAR_DUP_STR; dup_num++; sv_num = dup_num; break;
					case VAR_INV: sv_type = VAR_INV_STR; inv_num++; sv_num = inv_num; break;
					case VAR_TRA: sv_type = VAR_TRA_STR; tra_num++; sv_num = tra_num; break;
					default: sv_type = VAR_MIX_STR; break;
				}

				// line = reg->chrname + "\t" + to_string(reg->startRefPos) + "\t" + to_string(reg->endRefPos) + "\t" + sv_type;
				line = reg->chrname + "\t" + to_string(reg->startRefPos) + "\t" + to_string(reg->endRefPos) + "\t" + "ASVCLR." + sv_type + "." + to_string(sv_num) + "\t" + sv_type;
				if(reg->var_type!=VAR_TRA)
					line += "\t" + to_string(reg->sv_len);
				else
					line += "\t-";
				if(reg->var_type==VAR_DUP)
					line += "\t" + to_string(reg->dup_num);
				else
					line += "\t-";
				if(reg->refseq.size()) line += "\t" + reg->refseq;
				else line += "\t-";
				if(reg->altseq.size()) line += "\t" + reg->altseq;
				else line += "\t-";

				// extra information, split with ";"
				stringstream ss;
				ss << setprecision(3) << reg->AF;
				line += "\tDP=" + to_string(reg->DP) + ";AF=" + ss.str();

				discov_level_str = "";
				switch(reg->discover_level){
					case VAR_DISCOV_L_UNUSED:
						break;
					case VAR_DISCOV_L_CNS_ALN:
						discov_level_str = VAR_DISCOV_L_CNS_ALN_STR;
						break;
					case VAR_DISCOV_L_RESCUE_CNS_ALN:
						discov_level_str = VAR_DISCOV_L_RESCUE_CNS_ALN_STR;
						break;
					case VAR_DISCOV_L_READS:
						discov_level_str = VAR_DISCOV_L_READS_STR;
						break;
				}
				if(discov_level_str.size()>0) line = line + ";" + VAR_DISCOV_L_TITLE_STR + "=" + discov_level_str;

				if(reg->gt_seq.size()==0) reg->gt_seq = GT_STR_DEFAULT;
				line += "\t" + reg->gt_seq;

				if(reg->short_sv_flag) line += "\tShortSV";

				if(reg->var_type==VAR_UNC){
					cout << "line=" << __LINE__ << ": " << line << endl << endl << endl;
				}

				//cout << "line=" << __LINE__ << ": " << line << endl;

				outfile << line << endl;
			}else{ // not confirmed SV
//				if(var_cand->sv_type!=VAR_UNC){
//					file_id = 1;
//					switch(var_cand->sv_type){
//						case VAR_UNC: sv_type = "UNCERTAIN"; break;
//						case VAR_INS: sv_type = "INS"; file_id = 0; break;
//						case VAR_DEL: sv_type = "DEL"; file_id = 0; break;
//						case VAR_DUP: sv_type = "DUP"; break;
//						case VAR_INV: sv_type = "INV"; break;
//						case VAR_TRA: sv_type = "TRA"; break;
//						default: sv_type = "MIX"; break;
//					}
//					line = var_cand->chrname;
//					if(var_cand->leftClipRefPos>0 and var_cand->rightClipRefPos>0) line += "\t" + to_string(var_cand->leftClipRefPos) + "\t" + to_string(var_cand->rightClipRefPos);
//					else line += "\t-\t-";
//					line += "\t.\t.\t" + sv_type + "\t.";
//					if(file_id==1) outfile_clipReg << line << endl;
//					else outfile_indel << line << endl;
//				}
			}
		}
	}
	outfile.close();
}


