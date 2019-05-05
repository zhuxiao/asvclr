#include <sys/stat.h>
#include <pthread.h>
#include "Chrome.h"
#include "Thread.h"
#include "varCand.h"
#include "util.h"

// Constructor with parameters
Chrome::Chrome(string& chrname, int chrlen, faidx_t *fai, Paras *paras){
	this->paras = paras;
	this->chrname = chrname;
	this->chrlen = chrlen;
	this->fai = fai;
	init();

	blocks_out_file = chrname + "_blocks.bed";
	out_dir_detect = "";
	out_dir_assemble = "";
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
}

// set the output directory
void Chrome::setOutputDir(string& out_dir_detect_prefix, string& out_dir_assemble_prefix, string& out_dir_call_prefix){
	out_dir_detect = out_dir_detect_prefix + "/" + chrname;
	out_dir_assemble = out_dir_assemble_prefix + "/" + chrname;
	out_dir_call = out_dir_call_prefix + "/" + chrname;

	out_dir_detect = preprocessPipeChar(out_dir_detect);
	out_dir_assemble = preprocessPipeChar(out_dir_assemble);
	out_dir_call = preprocessPipeChar(out_dir_call);

	out_filename_detect_snv = out_dir_detect + "_SNV_candidate";
	out_filename_detect_indel = out_dir_detect + "_INDEL_candidate";
	out_filename_detect_clipReg = out_dir_detect + "_clipReg_candidate";
	out_filename_call_snv = out_dir_call + "_SNV";
	out_filename_call_indel = out_dir_call + "_INDEL";
	misAln_reg_filename = out_dir_detect + "_misaln_reg";
	out_filename_call_clipReg = out_dir_call + "_clipReg";

	var_cand_indel_filename = out_dir_assemble + "_var_cand_indel";
	var_cand_clipReg_filename = out_dir_assemble + "_var_cand_clipReg";
}

// generate the chromosome blocks
int Chrome::generateChrBlocks(){
	int pos, begPos, endPos;
	Block *block_tmp;
	bool headIgnFlag, tailIgnFlag;

	blockNum = 0;
	pos = 1;
	while(pos<=chrlen){
		begPos = pos;
		endPos = pos + paras->blockSize - 1;
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

		block_tmp = allocateBlock(chrname, chrlen, begPos, endPos, fai, headIgnFlag, tailIgnFlag);
		block_tmp->setOutputDir(out_dir_detect, out_dir_assemble, out_dir_call);
		blockVector.push_back(block_tmp);
	}
	blockVector.shrink_to_fit();

	return 0;
}

// allocate the memory for one block
Block *Chrome::allocateBlock(string& chrname, size_t chrlen, int startPos, int endPos, faidx_t *fai, bool headIgnFlag, bool tailIgnFlag){
	Block *block_tmp = new Block(chrname, chrlen, startPos, endPos, fai, paras);
	if(!block_tmp){
		cerr << "Chrome: cannot allocate memory" << endl;
		exit(1);
	}
	block_tmp->setRegIngFlag(headIgnFlag, tailIgnFlag);
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
	int i, pos, begPos, endPos, seq_len;
	Block *block_tmp;
	string reg;
	char *seq;
	bool flag;

	if(chrlen>=MIN_CHR_SIZE_EST){ // minimal chromosome size 50kb
		pos = chrlen / 2;
		while(pos<chrlen-(MIN_CHR_SIZE_EST-MIN_CHR_SIZE_EST)/2){
			flag = true;

			// construct a block with a 10kb size
			begPos = pos - BLOCK_SIZE_EST/2;
			endPos = begPos + BLOCK_SIZE_EST - 1;

			reg = chrname + ":" + to_string(begPos) + "-" + to_string(endPos);
			seq = fai_fetch(fai, reg.c_str(), &seq_len);
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
			cout << "Est region: " << chrname << ":" << begPos << "-" << endPos << endl;
			block_tmp = allocateBlock(chrname, chrlen, begPos, endPos, fai, false, false);
			// fill the data
			block_tmp->blockFillDataEst(op_est);
			delete block_tmp;

			paras->reg_sum_size_est += endPos - begPos + 1;
		}
	}
}

// detect indels for chrome
int Chrome::chrDetect(){
	mkdir(out_dir_detect.c_str(), S_IRWXU | S_IROTH);  // create the directory for detect command

	chrSetMisAlnRegFile();

	if(paras->num_threads<=1) chrDetect_st();  // single thread
	else chrDetect_mt();  // multiple threads

	// detect mated clip regions
	chrComputeMateClipReg();

	// remove FP indels and Snvs in clipping regions
	removeFPIndelSnvInClipReg(mateClipRegVector);

	// merge the results to single file
	//chrMergeDetectResultToFile();

	chrResetMisAlnRegFile();

	return 0;
}

// single thread
int Chrome::chrDetect_st(){
	Block* bloc;
	for(size_t i=0; i<blockVector.size(); i++){
		//if(i>=142)
		{
			bloc = blockVector.at(i);
			//cout << "detect files:" << bloc->out_dir_detect << "/" << bloc->snvFilenameDetect << ", " << bloc->indelFilenameDetect << endl;
			bloc->blockDetect();
		}
	}

	return 0;
}

// multiple threads
int Chrome::chrDetect_mt(){
	MultiThread mt[paras->num_threads];
	for(size_t i=0; i<paras->num_threads; i++){
		mt[i].setNumThreads(paras->num_threads);
		mt[i].setBlockVec(&blockVector);
		mt[i].setUserThreadID(i);
		if(!mt[i].startDetect()){
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

	return 0;
}

// detect mated clip regions for duplications and inversions
void Chrome::chrComputeMateClipReg(){
	size_t i;
	Block *bloc;
	reg_t *reg;
	vector<bool> clip_processed_flag_vec;

	for(i=0; i<blockVector.size(); i++){
		bloc = blockVector.at(i);
		for(size_t j=0; j<bloc->clipRegVector.size(); j++){
			clipRegVector.push_back(bloc->clipRegVector.at(j));
			clip_processed_flag_vec.push_back(false);
		}
	}
	clipRegVector.shrink_to_fit();

	// compute the mate flag for duplications and inversions
	for(i=0; i<clipRegVector.size(); i++){
		//if(i>=63)
		{
			reg = clipRegVector.at(i);
			if(clip_processed_flag_vec.at(i)==false){
				cout << i << ": " << reg->chrname << ":" << reg->startRefPos << "-" << reg->endRefPos << endl;

				clipReg clip_reg(reg->chrname, reg->startRefPos, reg->endRefPos, chrlen, paras->inBamFile, fai, paras);
				clip_reg.computeMateClipReg();

				processClipRegs(i, clip_processed_flag_vec, clip_reg.mate_clip_reg, reg);
			}
		}
	}

	// merge same TRA region
	mergeSameRegTRA();

	// remove FPs by detecting false overlapped clip regions for DUPs and INVs
	removeFPClipRegsDupInv();
}

// process clip regions and mate clip regions
void Chrome::processClipRegs(size_t idx, vector<bool> &clip_processed_flag_vec, mateClipReg_t &mate_clip_reg, reg_t *clip_reg){
	size_t i;
	reg_t *reg, *reg_tmp;
	mateClipReg_t *clip_reg_new;
	Block *bloc;
	int32_t idx_tmp;;

	if(mate_clip_reg.valid_flag){
		// add mate clip region
		if(mate_clip_reg.leftClipReg or mate_clip_reg.leftClipReg2 or mate_clip_reg.rightClipReg or mate_clip_reg.rightClipReg2){
			clip_reg_new = new mateClipReg_t();
			clip_reg_new->leftClipReg = mate_clip_reg.leftClipReg;
			clip_reg_new->leftClipPosNum = mate_clip_reg.leftClipPosNum;
			clip_reg_new->rightClipReg = mate_clip_reg.rightClipReg;
			clip_reg_new->rightClipPosNum = mate_clip_reg.rightClipPosNum;
			clip_reg_new->leftMeanClipPos = mate_clip_reg.leftMeanClipPos;
			clip_reg_new->rightMeanClipPos = mate_clip_reg.rightMeanClipPos;

			clip_reg_new->leftClipReg2 = mate_clip_reg.leftClipReg2;
			clip_reg_new->rightClipReg2 = mate_clip_reg.rightClipReg2;
			clip_reg_new->leftClipPosNum2 = mate_clip_reg.leftClipPosNum2;
			clip_reg_new->rightClipPosNum2 = mate_clip_reg.rightClipPosNum2;
			clip_reg_new->leftClipRegNum = mate_clip_reg.leftClipRegNum;
			clip_reg_new->rightClipRegNum = mate_clip_reg.rightClipRegNum;
			clip_reg_new->leftMeanClipPos2 = mate_clip_reg.leftMeanClipPos2;
			clip_reg_new->rightMeanClipPos2 = mate_clip_reg.rightMeanClipPos2;

			clip_reg_new->reg_mated_flag = mate_clip_reg.reg_mated_flag;
			clip_reg_new->chrname_leftTra1 = clip_reg_new->chrname_rightTra1 = clip_reg_new->chrname_leftTra2 = clip_reg_new->chrname_rightTra2 = "";
			clip_reg_new->leftClipPosTra1 = clip_reg_new->rightClipPosTra1 = clip_reg_new->leftClipPosTra2 = clip_reg_new->rightClipPosTra2 = -1;
			clip_reg_new->sv_type = mate_clip_reg.sv_type;
			clip_reg_new->dup_num = mate_clip_reg.dup_num;
			clip_reg_new->valid_flag = mate_clip_reg.valid_flag;
			clip_reg_new->call_success_flag = false;
			mateClipRegVector.push_back(clip_reg_new);
		}

		// detect overlapped clip regions
		if(mate_clip_reg.leftClipReg){ // left region 1
			for(i=0; i<clipRegVector.size(); i++){
				reg = clipRegVector.at(i);
				if(isOverlappedReg(mate_clip_reg.leftClipReg, reg)) // overlapped
					clip_processed_flag_vec.at(i) = true;
			}
		}
		if(mate_clip_reg.leftClipReg2){ // left region 2
			for(i=0; i<clipRegVector.size(); i++){
				reg = clipRegVector.at(i);
				if(isOverlappedReg(mate_clip_reg.leftClipReg2, reg)) // overlapped
					clip_processed_flag_vec.at(i) = true;
			}
		}
		if(mate_clip_reg.rightClipReg){ // right region 1
			for(i=0; i<clipRegVector.size(); i++){
				reg = clipRegVector.at(i);
				if(isOverlappedReg(mate_clip_reg.rightClipReg, reg)) // overlapped
					clip_processed_flag_vec.at(i) = true;
			}
		}
		if(mate_clip_reg.rightClipReg2){ // right region 2
			for(i=0; i<clipRegVector.size(); i++){
				reg = clipRegVector.at(i);
				if(isOverlappedReg(mate_clip_reg.rightClipReg2, reg)) // overlapped
					clip_processed_flag_vec.at(i) = true;
			}
		}
	}else{
		reg = new reg_t();
		reg->chrname = chrname;
		reg->startRefPos = clip_reg->startRefPos;
		reg->endRefPos = clip_reg->endRefPos;
		reg->startLocalRefPos = reg->endLocalRefPos = 0;
		reg->startQueryPos = reg->endQueryPos = 0;
		reg->var_type = VAR_UNC;
		reg->sv_len = 0;
		reg->query_id = -1;
		reg->blat_aln_id = -1;
		reg->call_success_status = false;
		reg->short_sv_flag = false;

		// get position
		idx_tmp = -1;
		bloc = blockVector.at(computeBlocID(clip_reg->startRefPos));
		for(i=0; i<bloc->indelVector.size(); i++){
			reg_tmp = bloc->indelVector.at(i);
			if(reg->startRefPos<reg_tmp->startRefPos){
				idx_tmp = i;
				break;
			}
		}
		// add item
		if(idx_tmp!=-1) bloc->indelVector.insert(bloc->indelVector.begin()+idx_tmp, reg);
		else bloc->indelVector.push_back(reg);
	}
}

// merge same TRA region
void Chrome::mergeSameRegTRA(){
	size_t i, end_flag1, pos_num1, pos_num2, mean_pos1, mean_pos2, tmp;
	mateClipReg_t *clip_reg, *clip_reg_mate_end, *clip_reg_retained, *clip_reg_redudant;
	reg_t *reg1, *reg2, *reg_tmp;

	for(i=0; i<mateClipRegVector.size(); i++){
		clip_reg = mateClipRegVector.at(i);
		if(clip_reg->valid_flag and clip_reg->reg_mated_flag and clip_reg->sv_type==VAR_TRA and (clip_reg->leftClipRegNum==2 or clip_reg->rightClipRegNum==2)){
			// get other end region according to the same clipping region
			clip_reg_mate_end = getMateRegEndSameClipReg(clip_reg, mateClipRegVector);
			if(clip_reg_mate_end){
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
					dist = (int32_t)reg2->startRefPos - (int32_t)reg->startRefPos;
					if(dist<0) dist = -dist;
					if(dist>(int32_t)paras->maxClipRegSize) mate_clip_reg->valid_flag = false;
					if(reg->startRefPos>reg2->endRefPos) mate_clip_reg->valid_flag = false;
				} // different chrome
			}
		}else mate_clip_reg->valid_flag = false;
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

	for(i=0; i<blockVector.size(); i++){

		//if(i>=246)
		{
			//cout << "<><><><><><><> " << i << endl;

			block = blockVector.at(i);
			for(j=0; j<block->indelVector.size(); ){
				reg = block->indelVector.at(j);
				flag = isIndelInClipReg(reg, mate_clipReg_vec);
				if(flag){
					delete reg;
					block->indelVector.erase(block->indelVector.begin()+j);
				}else j++;
			}
			for(j=0; j<block->snvVector.size(); ){
				pos = block->snvVector.at(j);
				flag = isSnvInClipReg(pos, mate_clipReg_vec);
				if(flag) block->snvVector.erase(block->snvVector.begin()+j);
				else j++;
			}
		}
	}
}

// determine whether the indel region in a clipping region
bool Chrome::isIndelInClipReg(reg_t *reg, vector<mateClipReg_t*> &mate_clipReg_vec){
	bool flag = false;
	mateClipReg_t *clip_reg;
	reg_t *reg_tmp;
	string chrname_tmp, chrname_tmp2;
	size_t start_pos, end_pos;

	for(size_t i=0; i<mate_clipReg_vec.size(); i++){
//		if(i>=0)
//			cout << ",,,,,,,,,,,, " << i << endl;

		clip_reg = mate_clipReg_vec.at(i);
		if(clip_reg->valid_flag and clip_reg->reg_mated_flag){
			if(clip_reg->sv_type==VAR_DUP or clip_reg->sv_type==VAR_INV){ // DUP or INV
				chrname_tmp = "";
				start_pos = end_pos = 0;
				if(clip_reg->leftClipReg){
					chrname_tmp = clip_reg->leftClipReg->chrname;
					start_pos = clip_reg->leftClipReg->startRefPos;
				}else if(clip_reg->leftClipReg2){
					chrname_tmp = clip_reg->leftClipReg2->chrname;
					start_pos = clip_reg->leftClipReg2->startRefPos;
				}
				if(clip_reg->rightClipReg){
					chrname_tmp2 = clip_reg->rightClipReg->chrname;
					end_pos = clip_reg->rightClipReg->endRefPos;
				}else if(clip_reg->rightClipReg2){
					chrname_tmp2 = clip_reg->rightClipReg2->chrname;
					end_pos = clip_reg->rightClipReg2->endRefPos;
				}

				if(reg->chrname.compare(chrname_tmp)==0 and chrname_tmp.compare(chrname_tmp2)==0 and start_pos>0 and end_pos>0){
					if(isOverlappedPos(reg->startRefPos, reg->endRefPos, start_pos, end_pos)){
						flag = true;
						break;
					}
				}
			}else if(clip_reg->sv_type==VAR_TRA){  // TRA
				// check left part
				chrname_tmp = "";
				start_pos = end_pos = 0;
				if(clip_reg->leftClipRegNum==2){
					chrname_tmp = clip_reg->leftClipReg->chrname;
					start_pos = clip_reg->leftClipReg->startRefPos;
					end_pos = clip_reg->leftClipReg2->endRefPos;
				}else if(clip_reg->leftClipRegNum==1){
					reg_tmp = clip_reg->leftClipReg ? clip_reg->leftClipReg : clip_reg->leftClipReg2;
					chrname_tmp = reg_tmp->chrname;
					start_pos = reg_tmp->startRefPos;
					end_pos = reg_tmp->endRefPos;
				}
				if(reg->chrname.compare(chrname_tmp)==0 and start_pos>0 and end_pos>0){
					if(isOverlappedPos(reg->startRefPos, reg->endRefPos, start_pos, end_pos)){
						flag = true;
						break;
					}
				}

				// check right part
				if(flag==false){
					chrname_tmp = "";
					start_pos = end_pos = 0;
					if(clip_reg->rightClipRegNum==2){
						chrname_tmp = clip_reg->rightClipReg->chrname;
						start_pos = clip_reg->rightClipReg->startRefPos;
						end_pos = clip_reg->rightClipReg2->endRefPos;
					}else if(clip_reg->rightClipRegNum==1){
						reg_tmp = clip_reg->rightClipReg ? clip_reg->rightClipReg : clip_reg->rightClipReg2;
						chrname_tmp = reg_tmp->chrname;
						start_pos = reg_tmp->startRefPos;
						end_pos = reg_tmp->endRefPos;
					}
					if(reg->chrname.compare(chrname_tmp)==0 and start_pos>0 and end_pos>0){
						if(isOverlappedPos(reg->startRefPos, reg->endRefPos, start_pos, end_pos)){
							flag = true;
							break;
						}
					}
				}
			}
		}
	}

	return flag;
}

// determine whether the SNV position in a clipping region
bool Chrome::isSnvInClipReg(size_t pos, vector<mateClipReg_t*> &mate_clipReg_vec){
	bool flag = false;
	mateClipReg_t *clip_reg;
	reg_t *reg_tmp;
	string chrname_tmp, chrname_tmp2;
	size_t start_pos, end_pos;

	for(size_t i=0; i<mate_clipReg_vec.size(); i++){
//		if(i>=0)
//			cout << ".............. " << i << endl;

		clip_reg = mate_clipReg_vec.at(i);
		if(clip_reg->valid_flag and clip_reg->reg_mated_flag){
			if(clip_reg->sv_type==VAR_DUP or clip_reg->sv_type==VAR_INV){ // DUP or INV
				chrname_tmp = "";
				start_pos = end_pos = 0;
				if(clip_reg->leftClipReg){
					chrname_tmp = clip_reg->leftClipReg->chrname;
					start_pos = clip_reg->leftClipReg->startRefPos;
				}else if(clip_reg->leftClipReg2){
					chrname_tmp = clip_reg->leftClipReg2->chrname;
					start_pos = clip_reg->leftClipReg2->startRefPos;
				}
				if(clip_reg->rightClipReg){
					chrname_tmp2 = clip_reg->rightClipReg->chrname;
					end_pos = clip_reg->rightClipReg->endRefPos;
				}else if(clip_reg->rightClipReg2){
					chrname_tmp2 = clip_reg->rightClipReg2->chrname;
					end_pos = clip_reg->rightClipReg2->endRefPos;
				}

				if(chrname.compare(chrname_tmp)==0 and chrname_tmp.compare(chrname_tmp2)==0 and start_pos>0 and end_pos>0){
					if(pos>=start_pos and pos<=end_pos){
						flag = true;
						break;
					}
				}
			}else if(clip_reg->sv_type==VAR_TRA){  // TRA
				// check left part
				chrname_tmp = "";
				start_pos = end_pos = 0;
				if(clip_reg->leftClipRegNum==2){
					chrname_tmp = clip_reg->leftClipReg->chrname;
					start_pos = clip_reg->leftClipReg->startRefPos;
					end_pos = clip_reg->leftClipReg2->endRefPos;
				}else if(clip_reg->leftClipRegNum==1){
					reg_tmp = clip_reg->leftClipReg ? clip_reg->leftClipReg : clip_reg->leftClipReg2;
					chrname_tmp = reg_tmp->chrname;
					start_pos = reg_tmp->startRefPos;
					end_pos = reg_tmp->endRefPos;
				}

				if(chrname.compare(chrname_tmp)==0 and start_pos>0 and end_pos>0){
					if(pos>=start_pos and pos<=end_pos){
						flag = true;
						break;
					}
				}

				// check right part
				if(flag==false){
					chrname_tmp = "";
					start_pos = end_pos = 0;
					if(clip_reg->rightClipRegNum==2){
						chrname_tmp = clip_reg->rightClipReg->chrname;
						start_pos = clip_reg->rightClipReg->startRefPos;
						end_pos = clip_reg->rightClipReg2->endRefPos;
					}else if(clip_reg->rightClipRegNum==1){
						reg_tmp = clip_reg->rightClipReg ? clip_reg->rightClipReg : clip_reg->rightClipReg2;
						chrname_tmp = reg_tmp->chrname;
						start_pos = reg_tmp->startRefPos;
						end_pos = reg_tmp->endRefPos;
					}

					if(chrname.compare(chrname_tmp)==0 and start_pos>0 and end_pos>0){
						if(pos>=start_pos and pos<=end_pos){
							flag = true;
							break;
						}
					}
				}
			}
		}
	}

	return flag;
}

// merge detect result to single file
void Chrome::chrMergeDetectResultToFile(){
	size_t i, j, pos;
	ofstream out_file_snv, out_file_indel, out_file_clipReg;
	string filename, line, var_type;
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
	out_file_clipReg.open(out_filename_detect_clipReg);
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

		line += "\t" + to_string(mate_clip_reg->leftClipPosNum) + "\t" + to_string(mate_clip_reg->leftClipPosNum2) + "\t" + to_string(mate_clip_reg->rightClipPosNum) + "\t" + to_string(mate_clip_reg->rightClipPosNum2);

		out_file_clipReg << line << endl;
		cout << line << endl;
	}
	out_file_clipReg.close();
}

// set assembly information file for local assembly
void Chrome::chrSetVarCandFiles(){
	var_cand_indel_file.open(var_cand_indel_filename);
	if(!var_cand_indel_file.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << var_cand_indel_filename << endl;
		exit(1);
	}
	var_cand_clipReg_file.open(var_cand_clipReg_filename);
	if(!var_cand_clipReg_file.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << var_cand_clipReg_filename << endl;
		exit(1);
	}

	vector<Block*>::iterator bloc;
	for(bloc=blockVector.begin(); bloc!=blockVector.end(); bloc++){
		(*bloc)->setVarCandFiles(&var_cand_indel_file, &var_cand_clipReg_file);
	}
}
// reset assembly information file for local assembly
void Chrome::chrResetVarCandFiles(){
	var_cand_indel_file.close();

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

// load detected variation data for local assembly
void Chrome::chrLoadDataAssemble(){
	size_t i;
	mateClipReg_t* mate_clip_reg;
	Block* tmp_bloc;

	if(paras->load_from_file_flag){ // load data from file
		chrLoadIndelDataAssemble();
		chrLoadClipRegDataAssemble();
	}

	// initialize block information
	for(i=0; i<mateClipRegVector.size(); i++){
		mate_clip_reg = mateClipRegVector.at(i);
		if(mate_clip_reg->leftClipReg and mate_clip_reg->leftClipReg->chrname.compare(chrname)==0) tmp_bloc = blockVector.at(computeBlocID(mate_clip_reg->leftClipReg->startRefPos));  // get the block
		else if(mate_clip_reg->rightClipReg and mate_clip_reg->rightClipReg->chrname.compare(chrname)==0) tmp_bloc = blockVector.at(computeBlocID(mate_clip_reg->rightClipReg->startRefPos));  // get the block
		else if(mate_clip_reg->leftClipReg2 and mate_clip_reg->leftClipReg2->chrname.compare(chrname)==0) tmp_bloc = blockVector.at(computeBlocID(mate_clip_reg->leftClipReg2->startRefPos));  // get the block
		else if(mate_clip_reg->rightClipReg2 and mate_clip_reg->rightClipReg2->chrname.compare(chrname)==0) tmp_bloc = blockVector.at(computeBlocID(mate_clip_reg->rightClipReg2->startRefPos));  // get the block
		tmp_bloc->mateClipRegVector.push_back(mate_clip_reg);
	}
}

// load detected indel data for local assembly
void Chrome::chrLoadIndelDataAssemble(){
	string line;
	vector<string> str_vec;
	ifstream infile;
	size_t begPos, endPos;
	Block* tmp_bloc;
	reg_t *reg;

	infile.open(out_filename_detect_indel);
	if(!infile.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << out_filename_detect_indel << endl;
		exit(1);
	}

	while(getline(infile, line))
	{
		if(line.size()){
			str_vec = split(line, "\t");
			begPos = stoi(str_vec[1]);
			endPos = stoi(str_vec[2]);
			tmp_bloc = blockVector.at(computeBlocID(begPos));  // get the block

			reg = new reg_t();
			reg->chrname = chrname;
			reg->startRefPos = begPos;
			reg->endRefPos = endPos;
			reg->startLocalRefPos = reg->endLocalRefPos = 0;
			reg->startQueryPos = reg->endQueryPos = 0;
			reg->var_type = VAR_UNC;
			reg->sv_len = 0;
			reg->query_id = -1;
			reg->blat_aln_id = -1;
			reg->call_success_status = false;
			reg->short_sv_flag = false;
			//cout << "blocID=" << computeBlocID(begPos) << ", reg:" << begPos << "-" << endPos << endl;

			tmp_bloc->indelVector.push_back(reg);  // add the variation
		}
	}
	infile.close();
}

// load detected indel data for local assembly
void Chrome::chrLoadClipRegDataAssemble(){
	chrLoadMateClipRegData();
}

void Chrome::chrLoadMateClipRegData(){
	string line, chrname1, chrname2, chrname3, chrname4;
	vector<string> str_vec;
	ifstream infile;
	//size_t startRefPos1, endRefPos1, startRefPos2, endRefPos2;
	reg_t *reg1, *reg2, *reg3, *reg4;
	bool mate_flag;
	mateClipReg_t* mate_clip_reg;
	size_t var_type, left_size, left_size2, right_size, right_size2, dup_num_tmp;

	infile.open(out_filename_detect_clipReg);
	if(!infile.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << out_filename_detect_clipReg << endl;
		exit(1);
	}

	while(getline(infile, line))
	{
		if(line.size()){
			reg1 = reg2 = reg3 = reg4 = NULL;
			str_vec = split(line, "\t");
			chrname1 = str_vec.at(0);
			if(chrname1.compare("-")!=0) { reg1 = new reg_t(); reg1->chrname = chrname1; reg1->startRefPos = stoi(str_vec.at(1)); reg1->endRefPos = stoi(str_vec.at(2)); }
			chrname2 = str_vec.at(3);
			if(chrname2.compare("-")!=0) { reg2 = new reg_t(); reg2->chrname = chrname2; reg2->startRefPos = stoi(str_vec.at(4)); reg2->endRefPos = stoi(str_vec.at(5)); }
			chrname3 = str_vec.at(6);
			if(chrname3.compare("-")!=0) { reg3 = new reg_t(); reg3->chrname = chrname3; reg3->startRefPos = stoi(str_vec.at(7)); reg3->endRefPos = stoi(str_vec.at(8)); }
			chrname4 = str_vec.at(9);
			if(chrname4.compare("-")!=0) { reg4 = new reg_t(); reg4->chrname = chrname4; reg4->startRefPos = stoi(str_vec.at(10)); reg4->endRefPos = stoi(str_vec.at(11)); }

			if(str_vec.at(12).compare("0")==0) mate_flag = false;
			else mate_flag = true;

			left_size = left_size2 = right_size = right_size2 = 0; var_type = VAR_UNC; dup_num_tmp = 0;
			if(str_vec.at(13).compare("####")==0){
				if(str_vec.at(14).compare("-")!=0) { left_size = stoi(str_vec.at(14)); }
				if(str_vec.at(15).compare("-")!=0) { left_size2 = stoi(str_vec.at(15)); }
				if(str_vec.at(16).compare("-")!=0) { right_size = stoi(str_vec.at(16)); }
				if(str_vec.at(17).compare("-")!=0) { right_size2 = stoi(str_vec.at(17)); }

				if(str_vec.at(18).compare("UNC")==0) var_type = VAR_UNC;
				else if(str_vec.at(18).compare("DUP")==0) var_type = VAR_DUP;
				else if(str_vec.at(18).compare("INV")==0) var_type = VAR_INV;
				else if(str_vec.at(18).compare("TRA")==0) var_type = VAR_TRA;
				else if(str_vec.at(18).compare("MIX")==0) var_type = VAR_MIX;

				if(str_vec.at(19).compare("-")!=0) { dup_num_tmp = stoi(str_vec.at(19)); }
			}

			mate_clip_reg = new mateClipReg_t();
			mate_clip_reg->leftClipReg = reg1;
			mate_clip_reg->leftClipReg2 = reg2;
			mate_clip_reg->rightClipReg = reg3;
			mate_clip_reg->rightClipReg2 = reg4;
			mate_clip_reg->reg_mated_flag = mate_flag;
			mate_clip_reg->valid_flag = true;
			mate_clip_reg->call_success_flag = false;
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

			//
			mate_clip_reg->chrname_leftTra1 = mate_clip_reg->chrname_rightTra1 = mate_clip_reg->chrname_leftTra2 = mate_clip_reg->chrname_rightTra2 = "";
			mate_clip_reg->leftClipPosTra1 = mate_clip_reg->rightClipPosTra1 = mate_clip_reg->leftClipPosTra2 = mate_clip_reg->rightClipPosTra2 = -1;
			mate_clip_reg->dup_num = dup_num_tmp;
			mateClipRegVector.push_back(mate_clip_reg);
		}
	}
	infile.close();
}

// compute block ID
int32_t Chrome::computeBlocID(size_t begPos){
	int32_t bloc_ID;

	bloc_ID = begPos / (paras->blockSize - 2 * paras->slideSize);
	if(bloc_ID>=blockNum) bloc_ID = blockNum - 1;

	return bloc_ID;
}

// local assembly for chrome
int Chrome::chrLocalAssemble(){
	mkdir(out_dir_assemble.c_str(), S_IRWXU | S_IROTH);  // create the directory for assemble command

	// open the assembly information file
	chrSetVarCandFiles();

	// local assembly
	if(paras->num_threads<=1) chrLocalAssemble_st();  // single thread
	else chrLocalAssemble_mt();  // multiple threads

	// close and reset the assembly information file
	chrResetVarCandFiles();

	return 0;
}

// local assembly for chrome using single thread
int Chrome::chrLocalAssemble_st(){
	Block *bloc;
	for(size_t i=0; i<blockVector.size(); i++){
//		if(i==1)
		{
			bloc = blockVector.at(i);
			bloc->blockLocalAssemble();
		}
	}
	return 0;
}

// local assembly for chrome using multiple threads
int Chrome::chrLocalAssemble_mt(){
	MultiThread mt[paras->num_threads];
	for(size_t i=0; i<paras->num_threads; i++){
		mt[i].setNumThreads(paras->num_threads);
		mt[i].setBlockVec(&blockVector);
		mt[i].setUserThreadID(i);
		if(!mt[i].startAssemble()){
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
	return 0;
}

// output assem data to file
void Chrome::outputAssemDataToFile(string &filename){
	string line, assembly_status, header, left_shift_size_str, right_shift_size_str;
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
		if(item->assem_success) assembly_status = ASSEMBLY_SUCCESS;
		else assembly_status = ASSEMBLY_FAILURE;
		line = item->refseqfilename + "\t" + item->ctgfilename + "\t" + item->readsfilename + "\t" + to_string(item->ref_left_shift_size) + "\t" + to_string(item->ref_right_shift_size) + "\t" + assembly_status;

		varVec = item->varVec;
		for(size_t j=0; j<varVec.size(); j++){
			reg = varVec[j];
			line += "\t" + reg->chrname + ":" + to_string(reg->startRefPos) + "-" + to_string(reg->endRefPos);
		}

		outfile << line << endl;
	}
	outfile.close();
}


int Chrome::chrCall(){
	mkdir(out_dir_call.c_str(), S_IRWXU | S_IROTH);  // create the directory for call command

	cout << "Begin calling variants ..." << endl;

	// call variants
	if(paras->num_threads<=1) chrCall_st();  // single thread
	else chrCall_mt();  // multiple threads

	// remove redundant new called variants
	cout << "11111111111" << endl;
	removeRedundantVar();

	// remove FPs of new called variants
	cout << "2222222222222" << endl;
	removeFPNewVarVec();

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
		//if(i>=350)
		{
			//cout << ">>>>>>>>> " << i << ", " << var_cand_vec.at(i)->alnfilename << endl;
			var_cand = var_cand_vec.at(i);
			var_cand->callVariants();
		}
	}
}

// call variants for chrome using multiple threads
void Chrome::chrCall_mt(){
	MultiThread mt[paras->num_threads];
	for(size_t i=0; i<paras->num_threads; i++){
		mt[i].setNumThreads(paras->num_threads);
		mt[i].setVarCandVec(&var_cand_vec, &var_cand_clipReg_vec);
		mt[i].setUserThreadID(i);
		if(!mt[i].startCall()){
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


void Chrome::chrLoadDataCall(){
	// load misAln regions
	if(mis_aln_vec.size()==0){
		loadMisAlnRegData();
		sortMisAlnRegData();
	}

	loadVarCandData();
	if(!isVarCandDataSorted(var_cand_vec)) { sortVarCandData(var_cand_vec); /*var_cand.outputAssemDataToFile(outfilename);*/ }
	if(!isVarCandDataSorted(var_cand_clipReg_vec)) { sortVarCandData(var_cand_clipReg_vec); /*var_cand.outputAssemDataToFile(outfilename);*/ }
}

void Chrome::loadVarCandData(){
	if(mateClipRegVector.size()==0) chrLoadMateClipRegData();
	loadVarCandDataFromFile(var_cand_vec, var_cand_indel_filename, false);
	loadVarCandDataFromFile(var_cand_clipReg_vec, var_cand_clipReg_filename, true);
}

void Chrome::loadVarCandDataFromFile(vector<varCand*> &var_cand_vec, string &var_cand_filename, bool clipReg_flag){
	string line, ctg_assembly_str, alnfilename, str_tmp, chrname_str;
	vector<string> line_vec, var_str, var_str2;
	vector<string> str_vec, str_vec2, str_vec3;
	ifstream infile;
	size_t i, lineNum;
	reg_t *reg, *reg1, *reg2;
	varCand *var_cand_tmp;
	mateClipReg_t *mate_clip_reg;
	int32_t clip_reg_idx_tra;

	infile.open(var_cand_filename);
	if(!infile.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << var_cand_filename << endl;
		exit(1);
	}

	lineNum = 0;
	while(getline(infile, line)){
		if(line.size()){
			// allocate memory
			var_cand_tmp = new varCand();

			var_cand_tmp->chrname = chrname;
			var_cand_tmp->var_cand_filename = var_cand_filename;
			var_cand_tmp->out_dir_call = out_dir_call;
			var_cand_tmp->misAln_filename = misAln_reg_filename;
			var_cand_tmp->inBamFile = paras->inBamFile;
			var_cand_tmp->fai = fai;

			line_vec = split(line, "\t");
			var_cand_tmp->refseqfilename = line_vec[0];  // refseq file name
			var_cand_tmp->ctgfilename = line_vec[1];  // contig file name
			var_cand_tmp->readsfilename = line_vec[2];  // reads file name
			var_cand_tmp->ref_left_shift_size = stoi(line_vec[3]);  // ref_left_shift_size
			var_cand_tmp->ref_right_shift_size = stoi(line_vec[4]);  // ref_right_shift_size

			if(line_vec[5].compare(ASSEMBLY_SUCCESS)==0) var_cand_tmp->assem_success = true;
			else var_cand_tmp->assem_success = false;

			// get the number of contigs
			var_cand_tmp->ctg_num = getCtgCount(var_cand_tmp->ctgfilename);

			// load variations
			for(i=6; i<line_vec.size(); i++){
				var_str = split(line_vec[i], ":");
				var_str2 = split(var_str[1], "-");
				reg = new reg_t();
				reg->chrname = var_str[0];
				reg->startRefPos = stoi(var_str2[0]);
				reg->endRefPos = stoi(var_str2[1]);
				reg->startLocalRefPos = reg->endLocalRefPos = 0;
				reg->startQueryPos = reg->endQueryPos = 0;
				reg->sv_len = 0;
				reg->dup_num = 0;
				reg->var_type = VAR_UNC;
				reg->query_id = -1;
				reg->blat_aln_id = -1;
				reg->call_success_status = false;
				reg->short_sv_flag = false;
				var_cand_tmp->varVec.push_back(reg);  // variation vector
			}
			var_cand_tmp->varVec.shrink_to_fit();

			// generate alignment file names
			str_vec = split(line_vec[1], "/");
			str_tmp = str_vec[str_vec.size()-1];  // contig file
			str_vec2 = split(str_tmp, "_");

			alnfilename = out_dir_call + "/blat";
			for(i=1; i<str_vec2.size()-1; i++)
				alnfilename += "_" + str_vec2[i];

			str_tmp = str_vec2[str_vec2.size()-1];
			str_vec3 = split(str_tmp, ".");

			alnfilename += "_" + str_vec3[0];
			for(i=1; i<str_vec3.size()-1; i++)
				alnfilename += "." + str_vec3[i];
			alnfilename += ".sim4";

			var_cand_tmp->alnfilename = alnfilename;
			var_cand_tmp->align_success = false;
			var_cand_tmp->clip_reg_flag = clipReg_flag;

			// assign clipping information
			mate_clip_reg = NULL;
			if(clipReg_flag) {
				reg1 = var_cand_tmp->varVec.at(0);
				reg2 = (var_cand_tmp->varVec.size()>=2) ? var_cand_tmp->varVec.at(var_cand_tmp->varVec.size()-1) : NULL;
				mate_clip_reg = getMateClipReg(reg1, reg2, &clip_reg_idx_tra);
			}
			if(mate_clip_reg){
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
					cerr << __func__ << ", line=" << __LINE__ << ", invalid variation type=" << mate_clip_reg->sv_type << ", error!" << endl;
					exit(1);
				}
			}else{
				var_cand_tmp->leftClipRefPos = var_cand_tmp->rightClipRefPos = 0;
				var_cand_tmp->sv_type = VAR_UNC;
				var_cand_tmp->dup_num = 0;
			}

			var_cand_vec.push_back(var_cand_tmp);
			lineNum ++;
		}
	}
	infile.close();

	cout << var_cand_filename << "\tlineNum=" << lineNum << endl;
}

// get the mate clip region
mateClipReg_t* Chrome::getMateClipReg(reg_t *reg1, reg_t *reg2, int32_t *clip_reg_idx_tra){
	mateClipReg_t *mate_clip_reg;
	size_t i;
	int32_t idx;
	reg_t *reg, *reg1_tmp, *reg2_tmp;

	idx = -1; *clip_reg_idx_tra = -1;
	for(i=0; i<mateClipRegVector.size(); i++){
		mate_clip_reg = mateClipRegVector.at(i);
		if(mate_clip_reg->valid_flag){
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
				cerr << __func__ << ", line=" << __LINE__ << ", invalid variation type=" << mate_clip_reg->sv_type << ", error!" << endl;
				exit(1);
			}
			if(idx!=-1) break;
		}
	}

	if(idx!=-1) mate_clip_reg = mateClipRegVector.at(idx);
	else mate_clip_reg = NULL;

	return mate_clip_reg;
}

// sort the assem in ascending order
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

	infile.open(misAln_reg_filename);
	if(!infile.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file" << misAln_reg_filename << endl;
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
			reg->call_success_status = false;
			reg->short_sv_flag = false;
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

void Chrome::removeVarCandNodeIndel(varCand *var_cand){
	removeVarCandNode(var_cand, var_cand_vec);
}

void Chrome::removeVarCandNodeClipReg(varCand *var_cand){
	removeVarCandNode(var_cand, var_cand_clipReg_vec);
}

void Chrome::removeVarCandNode(varCand *var_cand, vector<varCand*> &var_cand_vec){
	varCand *var_cand_tmp;
	for(size_t i=0; i<var_cand_vec.size(); i++){
		var_cand_tmp = var_cand_vec.at(i);
		if(var_cand_tmp==var_cand){
			var_cand_tmp->destroyVarCand();
			delete var_cand_tmp;
			var_cand_vec.erase(var_cand_vec.begin()+i);
			break;
		}
	}
}

// fill the sequence, including reference sequence and contig sequence
void Chrome::chrFillVarseq(){
	chrFillVarseqSingleVec(var_cand_vec);
	chrFillVarseqSingleVec(var_cand_clipReg_vec);
}

// fill the sequence, including reference sequence and contig sequence
void Chrome::chrFillVarseqSingleVec(vector<varCand*> &var_cand_vec){
	varCand *var_cand;
	for(size_t i=0; i<var_cand_vec.size(); i++){
		var_cand = var_cand_vec[i];
		var_cand->fillVarseq();
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

void Chrome::saveCallIndelClipReg2File(string &outfilename_indel, string &outfilename_clipReg){
	size_t i, j, file_id;
	ofstream outfile_indel, outfile_clipReg;
	varCand *var_cand;
	reg_t *reg;
	string line, sv_type;

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

	for(i=0; i<var_cand_vec.size(); i++){
		var_cand = var_cand_vec[i];
		for(j=0; j<var_cand->varVec.size(); j++){
			reg = var_cand->varVec[j];
			if(reg->call_success_status){
				file_id = 0;
				switch(reg->var_type){
					case VAR_UNC: sv_type = "UNCERTAIN"; break;
					case VAR_INS: sv_type = "INS"; break;
					case VAR_DEL: sv_type = "DEL"; break;
					case VAR_DUP: sv_type = "DUP"; file_id = 1; break;
					case VAR_INV: sv_type = "INV"; file_id = 1; break;
					case VAR_TRA: sv_type = "TRA"; file_id = 1; break;
					default: sv_type = "MIX"; break;
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
				line += "\t" + reg->refseq + "\t" + reg->altseq;

				if(reg->var_type==VAR_UNC){
					cout << "line=" << __LINE__ << ": " << line << ", short_sv_flag=" << reg->short_sv_flag << endl << endl << endl;
				}

//				line = reg->chrname + "\t" + to_string(reg->startRefPos) + "\t" + to_string(reg->endRefPos) + "\t" + sv_type + "\t" + reg->refseq + "\t" + reg->altseq;
//				if(reg->var_type==VAR_INS or reg->var_type==VAR_DEL)
//					line += "\t" + to_string(reg->sv_len);
//				else if(reg->var_type==VAR_DUP)
//					line += "\t" + to_string(reg->sv_len) + "\t" + to_string(reg->dup_num);
//				else
//					line += "\t.";
				if(file_id==0) outfile_indel << line << endl;
				else outfile_clipReg << line << endl;
			}
		}
	}

	for(i=0; i<var_cand_clipReg_vec.size(); i++){
		var_cand = var_cand_clipReg_vec[i];
		if(var_cand->clip_reg_flag==false){ // indel
			for(j=0; j<var_cand->varVec.size(); j++){
				reg = var_cand->varVec[j];
				if(reg->call_success_status){
					file_id = 0;
					switch(reg->var_type){
						case VAR_UNC: sv_type = "UNCERTAIN"; break;
						case VAR_INS: sv_type = "INS"; break;
						case VAR_DEL: sv_type = "DEL"; break;
						case VAR_DUP: sv_type = "DUP"; file_id = 1; break;
						case VAR_INV: sv_type = "INV"; file_id = 1; break;
						case VAR_TRA: sv_type = "TRA"; file_id = 1; break;
						default: sv_type = "MIX"; break;
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
					line += "\t" + reg->refseq + "\t" + reg->altseq;

					if(reg->var_type==VAR_UNC){
						cout << "line=" << __LINE__ << ": " << line << endl << endl << endl;
					}

//					line = reg->chrname + "\t" + to_string(reg->startRefPos) + "\t" + to_string(reg->endRefPos) + "\t" + sv_type + "\t" + reg->refseq + "\t" + reg->altseq;
//					if(reg->var_type==VAR_INS or reg->var_type==VAR_DEL)
//						line += "\t" + to_string(reg->sv_len);
//					else if(reg->var_type==VAR_DUP)
//						line += "\t" + to_string(reg->sv_len) + "\t" + to_string(reg->dup_num);
//					else
//						line += "\t.";

					if(file_id==0) outfile_indel << line << endl;
					else outfile_clipReg << line << endl;
				}
			}
		}else{ // cliping region
			reg = var_cand->clip_reg;
			if(var_cand->call_success){
				file_id = 1;
				switch(reg->var_type){
					case VAR_UNC: sv_type = "UNCERTAIN"; break;
					case VAR_INS: sv_type = "INS"; file_id = 0; break;
					case VAR_DEL: sv_type = "DEL"; file_id = 0; break;
					case VAR_DUP: sv_type = "DUP"; break;
					case VAR_INV: sv_type = "INV"; break;
					case VAR_TRA: sv_type = "TRA"; break;
					default: sv_type = "MIX"; break;
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
				line += "\t" + reg->refseq + "\t" + reg->altseq;

				if(reg->var_type==VAR_UNC){
					cout << "line=" << __LINE__ << ": " << line << endl << endl << endl;
				}

//				line = reg->chrname + "\t" + to_string(reg->startRefPos) + "\t" + to_string(reg->endRefPos) + "\t" + sv_type + "\t" + reg->refseq + "\t" + reg->altseq;
//				if(reg->var_type==VAR_DUP)
//					line += "\t" + to_string(reg->sv_len) + "\t" + to_string(reg->dup_num);
//				else
//					line += "\t.";
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

