#include <math.h>
#include "clipReg.h"
#include "clipAlnDataLoader.h"
#include "util.h"


clipReg::clipReg(string &chrname, size_t startRefPos, size_t endRefPos, size_t chrlen, string &inBamFile, faidx_t *fai, Paras *paras){
	this->chrname = chrname;
	this->startRefPos = startRefPos - CLIP_END_EXTEND_SIZE;
	this->endRefPos = endRefPos + CLIP_END_EXTEND_SIZE;
	this->chrlen = chrlen;
	this->inBamFile = inBamFile;
	this->fai = fai;
	this->paras = paras;

	if(this->startRefPos<1) this->startRefPos = 1;
	if(this->endRefPos>chrlen) this->endRefPos = chrlen;

	mate_clip_reg.leftClipReg = mate_clip_reg.leftClipReg2 = mate_clip_reg.rightClipReg = mate_clip_reg.rightClipReg2 = NULL;
	mate_clip_reg.leftClipRegNum = mate_clip_reg.rightClipRegNum = 0;
	mate_clip_reg.reg_mated_flag = false;
	mate_clip_reg.valid_flag = false;
	mate_clip_reg.sv_type = VAR_UNC;
	mate_clip_reg.dup_num = 0;
}

clipReg::~clipReg() {
	if(!clipAlnDataVector.empty()) destroyClipAlnDataVector();
}

void clipReg::destroyClipAlnDataVector(){
	for(size_t i=0; i<clipAlnDataVector.size(); i++){
		bam_destroy1(clipAlnDataVector.at(i)->bam);
		delete clipAlnDataVector.at(i);
	}
	vector<clipAlnData_t*>().swap(clipAlnDataVector);
}

// compute the mate clipping region
void clipReg::computeMateClipReg(){

	// fill clip align data vector
	fillClipAlnDataVectorWithSATag();

	// remove query having no clippings
	removeNonclipItems();

	// determine whether the clip region is valid
	if(isValidClipReg())
		computeMateAlnClipReg();
}

void clipReg::fillClipAlnDataVectorWithSATag(){
	clipAlnDataLoader clip_aln_data_loader(chrname, startRefPos, endRefPos, inBamFile);
	clip_aln_data_loader.loadClipAlnDataWithSATag(clipAlnDataVector, paras->max_ultra_high_cov);
}

// remove query having no clippings
void clipReg::removeNonclipItems(){
	clipAlnData_t *clip_aln;
	for(size_t i=0; i<clipAlnDataVector.size(); ){
		clip_aln = clipAlnDataVector.at(i);
		if(clip_aln->leftClipSize<MIN_CLIP_END_SIZE and clip_aln->rightClipSize<MIN_CLIP_END_SIZE){
			bam_destroy1(clip_aln->bam);
			delete clip_aln;
			clipAlnDataVector.erase(clipAlnDataVector.begin()+i);
		}else i++;
	}
	clipAlnDataVector.shrink_to_fit();
}

// determine whether the clip region is valid
bool clipReg::isValidClipReg(){
	bool valid_flag = false;
	clipAlnData_t *clip_aln;
	size_t i, aln_query_size, seg_num_non_SA, longer_seg_num_non_SA;

	seg_num_non_SA = 0; longer_seg_num_non_SA = 0;
	for(i=0; i<clipAlnDataVector.size(); i++){
		clip_aln = clipAlnDataVector.at(i);
		if(clip_aln->SA_tag_flag==false){
			seg_num_non_SA ++;
			if(clip_aln->aln_orient==ALN_PLUS_ORIENT) aln_query_size = clip_aln->endQueryPos - clip_aln->startQueryPos + 1;
			else aln_query_size = clip_aln->startQueryPos - clip_aln->endQueryPos + 1;
			if((double)aln_query_size/clip_aln->querylen>=0.5) {
				longer_seg_num_non_SA ++;
				//cout << i << ":" << clip_aln->queryname << ":" << clip_aln->chrname << ":" << clip_aln->startRefPos << "-" << clip_aln->endRefPos << "; " << "startQueryPos=" << clip_aln->startQueryPos << ", endQueryPos=" << clip_aln->endQueryPos << endl;
			}
		}
	}

	if(seg_num_non_SA>0 and longer_seg_num_non_SA>=LONGER_NON_SA_SEG_NUM_THRES and (double)longer_seg_num_non_SA/seg_num_non_SA>=LONGER_NON_SA_SEG_RATIO_THRES) valid_flag = true;

	//double ratio = (double)longer_seg_num_non_SA / seg_num_non_SA;
	//cout << "longer_non_seg_ratio=" << ratio << ", valid_flag=" << valid_flag << endl;

	return valid_flag;
}

// compute the mate clip region
void clipReg::computeMateAlnClipReg(){

	extractClipPosVec();  // get the clip pos vector

	splitClipPosVec(); // split vector

	//printClipVecs();  // print clips

	sortClipPos(); // sort clips

	//printClipVecs();  // print clips

	removeFakeClips();  // remove fake clips

	//printClipVecs();  // print clips

	computeClipRegs();  // get the mediate and the around region

	//removeFalseOverlappedMateClipReg();

	//printResultClipRegs();
}


// compute the mate clip region
void clipReg::extractClipPosVec(){
	size_t i, j;
	vector<clipAlnData_t*> query_aln_segs;
	clipAlnData_t *clip_aln_seg, *mate_clip_aln_seg;
	string queryname, clip_pos_str;
	vector<int32_t> adjClipAlnSegInfo;
	clipPos_t clip_pos_item, mate_clip_pos_item;
	int32_t arr_idx, clip_end_flag, mate_clip_end_flag, clip_vec_idx, mate_clip_vec_idx;
	bool valid_query_end_flag, valid_mate_query_end_flag, same_orient_flag, self_overlap_flag;

	for(i=0; i<clipAlnDataVector.size(); i++){
		if(clipAlnDataVector.at(i)->query_checked_flag==false){
			queryname = clipAlnDataVector.at(i)->queryname;
			query_aln_segs = getQueryClipAlnSegs(queryname, clipAlnDataVector);  // get query clip align segments

			for(j=0; j<query_aln_segs.size(); j++){
				clip_aln_seg = query_aln_segs.at(j);
				valid_query_end_flag = false;
				clip_end_flag = -1;
				if(clip_aln_seg->left_clip_checked_flag==false and clip_aln_seg->leftClipSize>=MIN_CLIP_END_SIZE and clip_aln_seg->startRefPos>=startRefPos and clip_aln_seg->startRefPos<=endRefPos){ // left end
					valid_query_end_flag = true;
					clip_end_flag = LEFT_END;
				}else if(clip_aln_seg->right_clip_checked_flag==false and clip_aln_seg->rightClipSize>=MIN_CLIP_END_SIZE and clip_aln_seg->endRefPos>=startRefPos and clip_aln_seg->endRefPos<=endRefPos){ // right end
					valid_query_end_flag = true;
					clip_end_flag = RIGHT_END;
				}

				valid_mate_query_end_flag = false;
				clip_vec_idx = mate_clip_vec_idx = -1;
				same_orient_flag = true;
				if(valid_query_end_flag){
					//cout << "\tseg_len:" << clip_aln_seg->endRefPos - clip_aln_seg->startRefPos << endl;
					// deal with clip end itself
					clip_pos_item.chrname = clip_aln_seg->chrname;
					if(clip_end_flag==LEFT_END){ // left clip end
						clip_pos_item.clipRefPos = clip_aln_seg->startRefPos;
						clip_aln_seg->left_clip_checked_flag = true;
						clip_vec_idx = 0;
					}else{ // right clip end
						clip_pos_item.clipRefPos = clip_aln_seg->endRefPos;
						clip_aln_seg->right_clip_checked_flag = true;
						clip_vec_idx = 1;
					}

					// deal with the mate clip end
					adjClipAlnSegInfo = getAdjacentClipAlnSeg(j, clip_end_flag, query_aln_segs);
					arr_idx = adjClipAlnSegInfo.at(0);
					mate_clip_end_flag = adjClipAlnSegInfo.at(1);
					if(arr_idx!=-1){ // mated
						mate_clip_aln_seg = query_aln_segs.at(arr_idx);
						if((mate_clip_end_flag==LEFT_END and mate_clip_aln_seg->left_clip_checked_flag==false) or (mate_clip_end_flag==RIGHT_END and mate_clip_aln_seg->right_clip_checked_flag==false)){
							//cout << "\tmate_seg_len:" << mate_clip_aln_seg->endRefPos - mate_clip_aln_seg->startRefPos << endl;

							mate_clip_pos_item.chrname = mate_clip_aln_seg->chrname;
							if(mate_clip_end_flag==LEFT_END){
								mate_clip_pos_item.clipRefPos = mate_clip_aln_seg->startRefPos;
								mate_clip_aln_seg->left_clip_checked_flag = true;
								mate_clip_vec_idx = 0;
							}else{
								mate_clip_pos_item.clipRefPos = mate_clip_aln_seg->endRefPos;
								mate_clip_aln_seg->right_clip_checked_flag = true;
								mate_clip_vec_idx = 1;
							}
							valid_mate_query_end_flag = true;

							if(clip_aln_seg->aln_orient!=mate_clip_aln_seg->aln_orient) same_orient_flag = false;

							// determine which is the left and which is the right clip region for INV and TRA
							self_overlap_flag = isSegSelfOverlap(clip_aln_seg, mate_clip_aln_seg);
							if(self_overlap_flag==false){
								if(clip_pos_item.clipRefPos<=mate_clip_pos_item.clipRefPos){
									clip_vec_idx = 0;
									mate_clip_vec_idx = 1;
								}else{
									clip_vec_idx = 1;
									mate_clip_vec_idx = 0;
								}
							}
						}
					}

					// save to vector
					if(valid_query_end_flag){
						clip_pos_item.same_orient_flag = same_orient_flag;
						if(clip_vec_idx==0)
							leftClipPosVector.push_back(clip_pos_item);
						else if(clip_vec_idx==1)
							rightClipPosVector.push_back(clip_pos_item);
					}
					if(valid_mate_query_end_flag){
						mate_clip_pos_item.same_orient_flag = same_orient_flag;
						if(mate_clip_vec_idx==0)
							leftClipPosVector.push_back(mate_clip_pos_item);
						else if(mate_clip_vec_idx==1)
							rightClipPosVector.push_back(mate_clip_pos_item);
					}
				}
			}
			for(j=0; j<query_aln_segs.size(); j++) query_aln_segs.at(j)->query_checked_flag = true;
		}
	}
}

// get query clip align segments
vector<clipAlnData_t*> clipReg::getQueryClipAlnSegs(string &queryname, vector<clipAlnData_t*> &clipAlnDataVector){
	vector<clipAlnData_t*> query_aln_segs;
	for(size_t i=0; i<clipAlnDataVector.size(); i++)
		if(clipAlnDataVector.at(i)->query_checked_flag==false and clipAlnDataVector.at(i)->queryname==queryname)
			query_aln_segs.push_back(clipAlnDataVector.at(i));
	return query_aln_segs;
}

// get adjacent clip segment according to query positions
vector<int32_t> clipReg::getAdjacentClipAlnSeg(int32_t arr_idx, int32_t clip_end_flag, vector<clipAlnData_t*> &query_aln_segs){
	clipAlnData_t *clip_aln;
	size_t i, clip_pos_based;
	int32_t dist, min_dist, idx_min_dist, end_flag;
	vector<int32_t> adjClipAlnSegInfo; // [0]: array index of minimal distance; [1]: segment end flag

	if(clip_end_flag==LEFT_END) clip_pos_based = query_aln_segs.at(arr_idx)->startQueryPos;
	else clip_pos_based = query_aln_segs.at(arr_idx)->endQueryPos;

	min_dist = INT_MAX;
	idx_min_dist = -1;
	end_flag = -1;
	for(i=0; i<query_aln_segs.size(); i++){
		if(i!=(size_t)arr_idx){
			clip_aln = query_aln_segs.at(i);
			// left end
			dist = clip_aln->startQueryPos - clip_pos_based;
			if(dist<0) dist = -dist;
			if(dist<min_dist) {
				min_dist = dist;
				idx_min_dist = i;
				end_flag = LEFT_END;
			}
			// right end
			dist = clip_aln->endQueryPos - clip_pos_based;
			if(dist<0) dist = -dist;
			if(dist<min_dist) {
				min_dist = dist;
				idx_min_dist = i;
				end_flag = RIGHT_END;
			}
		}
	}

	adjClipAlnSegInfo.push_back(idx_min_dist);
	adjClipAlnSegInfo.push_back(end_flag);

	return adjClipAlnSegInfo;
}

void clipReg::splitClipPosVec(){
	size_t i;
	clipPos_t clip_pos;

	for(i=0; i<leftClipPosVector.size(); ){
		clip_pos = leftClipPosVector.at(i);
		if(clip_pos.chrname.compare(chrname)==0 and (clip_pos.clipRefPos>=startRefPos and clip_pos.clipRefPos<=endRefPos)) i++;
		else if(clip_pos.same_orient_flag==false){
			leftClipPosVector2.push_back(clip_pos);
			leftClipPosVector.erase(leftClipPosVector.begin()+i);
		}else i++;
	}

	for(i=0; i<rightClipPosVector.size(); ){
		clip_pos = rightClipPosVector.at(i);
		if(clip_pos.chrname.compare(chrname)==0 and (clip_pos.clipRefPos>=startRefPos and clip_pos.clipRefPos<=endRefPos)) i++;
		else if(clip_pos.same_orient_flag==false){
			rightClipPosVector2.push_back(clip_pos);
			rightClipPosVector.erase(rightClipPosVector.begin()+i);
		}else i++;
	}
}

void clipReg::sortClipPos(){
	sortClipPosSingleVec(leftClipPosVector);
	sortClipPosSingleVec(leftClipPosVector2);
	sortClipPosSingleVec(rightClipPosVector);
	sortClipPosSingleVec(rightClipPosVector2);
}

void clipReg::sortClipPosSingleVec(vector<clipPos_t> &clipPosVec){
	size_t i, j, minPos, minIdx, num, tmp;
	string tmp_str;

	// sort vector ascending
	for(i=0; i<clipPosVec.size(); i++){
		minPos = clipPosVec.at(i).clipRefPos;
		minIdx = i;
		for(j=i+1; j<clipPosVec.size(); j++){
			num = clipPosVec.at(j).clipRefPos;
			if(num<minPos){
				minPos = num;
				minIdx = j;
			}
		}

		tmp_str = clipPosVec.at(i).chrname;
		tmp = clipPosVec.at(i).clipRefPos;
		clipPosVec.at(i).chrname = clipPosVec.at(minIdx).chrname;
		clipPosVec.at(i).clipRefPos = clipPosVec.at(minIdx).clipRefPos;
		clipPosVec.at(minIdx).chrname = tmp_str;
		clipPosVec.at(minIdx).clipRefPos = tmp;
	}
}

// remove fake clips
void clipReg::removeFakeClips(){
	removeFakeClipsDifferentChrSingleVec(leftClipPosVector);
	removeFakeClipsDifferentChrSingleVec(leftClipPosVector2);
	removeFakeClipsDifferentChrSingleVec(rightClipPosVector);
	removeFakeClipsDifferentChrSingleVec(rightClipPosVector2);

	removeFakeClipsLongDistSameOrientSingleVec(leftClipPosVector);
	removeFakeClipsLongDistSameOrientSingleVec(leftClipPosVector2);
	removeFakeClipsLongDistSameOrientSingleVec(rightClipPosVector);
	removeFakeClipsLongDistSameOrientSingleVec(rightClipPosVector2);

	removeFakeClipsLowCov(leftClipPosVector, paras->minClipReadsNumSupportSV);
	removeFakeClipsLowCov(leftClipPosVector2, paras->minClipReadsNumSupportSV);
	removeFakeClipsLowCov(rightClipPosVector, paras->minClipReadsNumSupportSV);
	removeFakeClipsLowCov(rightClipPosVector2, paras->minClipReadsNumSupportSV);
}

// remove fake clips with alignment onto different chromes based on single vector
void clipReg::removeFakeClipsDifferentChrSingleVec(vector<clipPos_t> &clipPosVector){
	size_t i, j, sub_sum, maxValue;
	int32_t maxIdx;
	vector<string> chrname_vec;
	string chrname_tmp;
	vector<size_t> num_vec;

	// get chrnames
	for(i=0; i<clipPosVector.size(); i++)
		if(isExistStr(clipPosVector.at(i).chrname, chrname_vec)==false) // not exist in vector, add it to vector
			chrname_vec.push_back(clipPosVector.at(i).chrname);

	// compute the number of clips for each chrome
	for(i=0; i<chrname_vec.size(); i++){
		chrname_tmp = chrname_vec.at(i);
		sub_sum = 0;
		for(j=0; j<clipPosVector.size(); j++) if(chrname_tmp.compare(clipPosVector.at(j).chrname)==0) sub_sum ++;
		num_vec.push_back(sub_sum);
	}

	// get the chrome with the most clips
	maxIdx = -1;
	maxValue = 0;
	for(i=0; i<num_vec.size(); i++){
		if(num_vec.at(i)>maxValue){
			maxValue = num_vec.at(i);
			maxIdx = i;
		}
	}

	// remove fake clips
	if(maxIdx>=0){
		chrname_tmp = chrname_vec.at(maxIdx);
		for(i=0; i<clipPosVector.size(); ){
			if(chrname_tmp.compare(clipPosVector.at(i).chrname)!=0){
				clipPosVector.erase(clipPosVector.begin()+i);
			}else i++;
		}
	}
}

// remove fake clips with long dist based on single vector
void clipReg::removeFakeClipsLongDistSameOrientSingleVec(vector<clipPos_t> &clipPosVector){
	size_t i, j;
	int32_t idx, dist;
	vector<size_t> clip_pos_num_vec;
	vector<clipPos_t> clip_pos_vec;
	clipPos_t clip_pos_item, clip_pos_item_tmp;

	size_t minPos, minIdx, maxValue, num, tmp;
	string tmp_str;

	for(i=0; i<clipPosVector.size(); i++){
		clip_pos_item = clipPosVector.at(i);
		idx = getItemIdxClipPosVec(clip_pos_item, clip_pos_vec);
		if(idx==-1){ // new item
			clip_pos_vec.push_back(clip_pos_item);
			clip_pos_num_vec.push_back(1);
		}else
			clip_pos_num_vec.at(idx) ++;
	}

	// sort vector ascending
	for(i=0; i<clip_pos_vec.size(); i++){
		minPos = clip_pos_vec.at(i).clipRefPos;
		minIdx = i;
		for(j=i+1; j<clip_pos_vec.size(); j++){
			num = clip_pos_vec.at(j).clipRefPos;
			if(num<minPos){
				minPos = num;
				minIdx = j;
			}
		}

		tmp = clip_pos_num_vec.at(i);
		clip_pos_num_vec.at(i) = clip_pos_num_vec.at(minIdx);
		clip_pos_num_vec.at(minIdx) = tmp;

		tmp_str = clip_pos_vec.at(i).chrname;
		tmp = clip_pos_vec.at(i).clipRefPos;
		clip_pos_vec.at(i).chrname = clip_pos_vec.at(minIdx).chrname;
		clip_pos_vec.at(i).clipRefPos = clip_pos_vec.at(minIdx).clipRefPos;
		clip_pos_vec.at(minIdx).chrname = tmp_str;
		clip_pos_vec.at(minIdx).clipRefPos = tmp;
	}

	// get maxIdx and maxValue
	idx = -1;
	maxValue = 0;
	for(i=0; i<clip_pos_num_vec.size(); i++){
		num = clip_pos_num_vec.at(i);
		if(maxValue<num){
			maxValue = num;
			idx = i;
		}
	}

	// remove false clips
	if(idx>=0){
		clip_pos_item = clip_pos_vec.at(idx);
		for(i=0; i<clipPosVector.size(); ){
			clip_pos_item_tmp = clipPosVector.at(i);
			dist = clip_pos_item_tmp.clipRefPos - clip_pos_item.clipRefPos;
			if(dist<0) dist = -dist;
			if(dist>MIN_CLIP_DIST_THRES) // invalid
				clipPosVector.erase(clipPosVector.begin()+i);
			else i++;
		}
	}
}

// remove fake clips with long dist based on single vector
void clipReg::removeFakeClipsLowCov(vector<clipPos_t> &clipPosVector, size_t min_clip_reads_num){
	// remove fake clips
	if(clipPosVector.size()<min_clip_reads_num) clipPosVector.clear();
}

// get the index
int32_t clipReg::getItemIdxClipPosVec(clipPos_t &item, vector<clipPos_t> &vec){
	int32_t idx = -1;
	for(size_t i=0; i<vec.size(); i++)
		if(item.chrname.compare(vec.at(i).chrname)==0 and item.clipRefPos==vec.at(i).clipRefPos){
			idx = i;
			break;
		}
	return idx;
}

// compute the clip region
void clipReg::computeClipRegs(){
	mate_clip_reg.leftClipReg = computeClipRegSingleVec(leftClipPosVector);
	mate_clip_reg.rightClipReg = computeClipRegSingleVec(rightClipPosVector);
	mate_clip_reg.leftClipPosNum = leftClipPosVector.size();
	mate_clip_reg.rightClipPosNum = rightClipPosVector.size();
	mate_clip_reg.leftMeanClipPos = computeMeanClipPos(leftClipPosVector);
	mate_clip_reg.rightMeanClipPos = computeMeanClipPos(rightClipPosVector);

	mate_clip_reg.leftClipReg2 = computeClipRegSingleVec(leftClipPosVector2);
	mate_clip_reg.rightClipReg2 = computeClipRegSingleVec(rightClipPosVector2);
	mate_clip_reg.leftClipPosNum2 = leftClipPosVector2.size();
	mate_clip_reg.rightClipPosNum2 = rightClipPosVector2.size();
	mate_clip_reg.leftMeanClipPos2 = computeMeanClipPos(leftClipPosVector2);
	mate_clip_reg.rightMeanClipPos2 = computeMeanClipPos(rightClipPosVector2);

	mate_clip_reg.leftClipRegNum = 0;
	if(mate_clip_reg.leftClipReg) mate_clip_reg.leftClipRegNum ++;
	if(mate_clip_reg.leftClipReg2) mate_clip_reg.leftClipRegNum ++;
	mate_clip_reg.rightClipRegNum = 0;
	if(mate_clip_reg.rightClipReg) mate_clip_reg.rightClipRegNum ++;
	if(mate_clip_reg.rightClipReg2) mate_clip_reg.rightClipRegNum ++;

	// remove FP region for single clip end
	removeFPClipSingleEnd(mate_clip_reg);

	// check location order
	checkLocOrder(mate_clip_reg);

	mate_clip_reg.valid_flag = true;
	if((mate_clip_reg.leftClipReg or mate_clip_reg.leftClipReg2) and (mate_clip_reg.rightClipReg or mate_clip_reg.rightClipReg2)) mate_clip_reg.reg_mated_flag = true;

	// sv_type and dup_num
	mate_clip_reg.sv_type = VAR_UNC;
	mate_clip_reg.dup_num = 0;
	if(mate_clip_reg.reg_mated_flag) // sv_type
		computeVarTypeClipReg(mate_clip_reg, inBamFile);

	if(mate_clip_reg.sv_type==VAR_UNC) mate_clip_reg.valid_flag = false;
}

// compute the clip region based on single vector
reg_t* clipReg::computeClipRegSingleVec(vector<clipPos_t> &clipPosVector){
	reg_t *reg = NULL;
	size_t i, minValue, maxValue;
	int32_t minIdx, maxIdx;

	minIdx = -1, minValue = INT_MAX;
	maxIdx = -1, maxValue = 0;
	for(i=0; i<clipPosVector.size(); i++){
		if(minValue>clipPosVector.at(i).clipRefPos){
			minValue = clipPosVector.at(i).clipRefPos;
			minIdx = i;
		}
		if(maxValue<clipPosVector.at(i).clipRefPos){
			maxValue = clipPosVector.at(i).clipRefPos;
			maxIdx = i;
		}
	}

	if(minIdx!=-1 and maxIdx!=-1){
		reg = new reg_t();
		reg->chrname = clipPosVector.at(minIdx).chrname;
		reg->startRefPos = clipPosVector.at(minIdx).clipRefPos;
		reg->endRefPos = clipPosVector.at(maxIdx).clipRefPos;
		reg->zero_cov_flag = false;
	}

	return reg;
}

// remove false overlapped clip region
void clipReg::removeFalseOverlappedMateClipReg(){
	if(mate_clip_reg.reg_mated_flag){
		if(isOverlappedReg(mate_clip_reg.leftClipReg, mate_clip_reg.rightClipReg)){
			delete mate_clip_reg.leftClipReg;
			delete mate_clip_reg.rightClipReg;
			mate_clip_reg.leftClipReg = NULL;
			mate_clip_reg.rightClipReg = NULL;
		}
	}
}

void clipReg::removeFPClipSingleEnd(mateClipReg_t &mate_clip_reg){
	size_t dist;
	bool valid_flag;

	if(mate_clip_reg.leftClipRegNum==2){
		valid_flag = true;
		if(mate_clip_reg.leftClipReg->chrname.compare(mate_clip_reg.leftClipReg2->chrname)==0){
			if(mate_clip_reg.leftClipReg->startRefPos<mate_clip_reg.leftClipReg2->startRefPos) dist = mate_clip_reg.leftClipReg2->startRefPos - mate_clip_reg.leftClipReg->startRefPos;
			else dist = mate_clip_reg.leftClipReg->startRefPos - mate_clip_reg.leftClipReg2->startRefPos;
			if(dist>MAX_DIST_SAME_CLIP_END) valid_flag = false;
		}else valid_flag = false;
		if(valid_flag==false){
			if(mate_clip_reg.leftClipPosNum>mate_clip_reg.leftClipPosNum2){
				delete mate_clip_reg.leftClipReg2;
				mate_clip_reg.leftClipReg2 = NULL;
				mate_clip_reg.leftClipPosNum2 = 0;
				mate_clip_reg.leftMeanClipPos2 = 0;
			}else{
				delete mate_clip_reg.leftClipReg;
				mate_clip_reg.leftClipReg = NULL;
				mate_clip_reg.leftClipPosNum = 0;
				mate_clip_reg.leftMeanClipPos = 0;
			}
			mate_clip_reg.leftClipRegNum --;
		}
	}
	if(mate_clip_reg.rightClipRegNum==2){
		valid_flag = true;
		if(mate_clip_reg.rightClipReg->chrname.compare(mate_clip_reg.rightClipReg2->chrname)==0){
			if(mate_clip_reg.rightClipReg->startRefPos<mate_clip_reg.rightClipReg2->startRefPos) dist = mate_clip_reg.rightClipReg2->startRefPos - mate_clip_reg.rightClipReg->startRefPos;
			else dist = mate_clip_reg.rightClipReg->startRefPos - mate_clip_reg.rightClipReg2->startRefPos;
			if(dist>MAX_DIST_SAME_CLIP_END) valid_flag = false;
		}else valid_flag = false;
		if(valid_flag==false){
			if(mate_clip_reg.rightClipPosNum>mate_clip_reg.rightClipPosNum2){
				delete mate_clip_reg.rightClipReg2;
				mate_clip_reg.rightClipReg2 = NULL;
				mate_clip_reg.rightClipPosNum2 = 0;
				mate_clip_reg.rightMeanClipPos2 = 0;
			}else{
				delete mate_clip_reg.rightClipReg;
				mate_clip_reg.rightClipReg = NULL;
				mate_clip_reg.rightClipPosNum = 0;
				mate_clip_reg.rightMeanClipPos = 0;
			}
			mate_clip_reg.rightClipRegNum --;
		}
	}
}

// compute the mean size of the clippings
size_t clipReg::computeMeanClipPos(vector<clipPos_t> &clipPosVector){
	size_t i, total, mean_pos;
	total = 0;
	for(i=0; i<clipPosVector.size(); i++) total += clipPosVector.at(i).clipRefPos;
	if(total>0) mean_pos = round((double)total/clipPosVector.size());
	else mean_pos = 0;
	return mean_pos;
}

// check SV location order
void clipReg::checkLocOrder(mateClipReg_t &mate_clip_reg){
	size_t left_loc, right_loc, tmp;
	reg_t *reg_tmp;
	string chrname_left, chrname_right;

	if(mate_clip_reg.leftClipRegNum==2){
		chrname_left = chrname_right = "";
		left_loc = right_loc = 0;
		if(mate_clip_reg.leftMeanClipPos>mate_clip_reg.leftMeanClipPos2){
			left_loc = mate_clip_reg.leftMeanClipPos;
			chrname_left = mate_clip_reg.leftClipReg->chrname;
			right_loc = mate_clip_reg.leftMeanClipPos2;
			chrname_right = mate_clip_reg.leftClipReg2->chrname;
			if(chrname_left.size()>0 and chrname_left.compare(chrname_right)==0 and left_loc>right_loc){
				// SV region
				reg_tmp = mate_clip_reg.leftClipReg;
				mate_clip_reg.leftClipReg = mate_clip_reg.leftClipReg2;
				mate_clip_reg.leftClipReg2 = reg_tmp;
				// leftClipPosNum
				tmp = mate_clip_reg.leftClipPosNum;
				mate_clip_reg.leftClipPosNum = mate_clip_reg.leftClipPosNum2;
				mate_clip_reg.leftClipPosNum2 = tmp;
				// leftMeanClipPos
				tmp = mate_clip_reg.leftMeanClipPos;
				mate_clip_reg.leftMeanClipPos = mate_clip_reg.leftMeanClipPos2;
				mate_clip_reg.leftMeanClipPos2 = tmp;
			}
		}
	}
	if(mate_clip_reg.rightClipRegNum==2){
		chrname_left = chrname_right = "";
		left_loc = right_loc = 0;
		if(mate_clip_reg.rightMeanClipPos>mate_clip_reg.rightMeanClipPos2){
			left_loc = mate_clip_reg.rightMeanClipPos;
			chrname_left = mate_clip_reg.rightClipReg->chrname;
			right_loc = mate_clip_reg.rightMeanClipPos2;
			chrname_right = mate_clip_reg.rightClipReg2->chrname;
			if(chrname_left.size()>0 and chrname_left.compare(chrname_right)==0 and left_loc>right_loc){
				// SV region
				reg_tmp = mate_clip_reg.rightClipReg;
				mate_clip_reg.rightClipReg = mate_clip_reg.rightClipReg2;
				mate_clip_reg.rightClipReg2 = reg_tmp;
				// rightClipPosNum
				tmp = mate_clip_reg.rightClipPosNum;
				mate_clip_reg.rightClipPosNum = mate_clip_reg.rightClipPosNum2;
				mate_clip_reg.rightClipPosNum2 = tmp;
				// rightMeanClipPos
				tmp = mate_clip_reg.rightMeanClipPos;
				mate_clip_reg.rightMeanClipPos = mate_clip_reg.rightMeanClipPos2;
				mate_clip_reg.rightMeanClipPos2 = tmp;
			}
		}
	}

	chrname_left = chrname_right = "";
	left_loc = right_loc = 0;
	if(mate_clip_reg.leftClipRegNum==2 and mate_clip_reg.leftMeanClipPos<mate_clip_reg.leftMeanClipPos2) { left_loc = mate_clip_reg.leftMeanClipPos2; chrname_left = mate_clip_reg.leftClipReg2->chrname; }
	else if(mate_clip_reg.leftClipRegNum==1){
		if(mate_clip_reg.leftMeanClipPos>0) { left_loc = mate_clip_reg.leftMeanClipPos; chrname_left = mate_clip_reg.leftClipReg->chrname; }
		else if(mate_clip_reg.leftMeanClipPos2>0) { left_loc = mate_clip_reg.leftMeanClipPos2; chrname_left = mate_clip_reg.leftClipReg2->chrname; }
	}
	if(mate_clip_reg.rightClipRegNum==2 and mate_clip_reg.rightMeanClipPos<mate_clip_reg.rightMeanClipPos2) { right_loc = mate_clip_reg.rightMeanClipPos; chrname_right = mate_clip_reg.rightClipReg->chrname; }
	else if(mate_clip_reg.rightClipRegNum==1){
		if(mate_clip_reg.rightMeanClipPos>0) { right_loc = mate_clip_reg.rightMeanClipPos; chrname_right = mate_clip_reg.rightClipReg->chrname; }
		else if(mate_clip_reg.rightMeanClipPos2>0) { right_loc = mate_clip_reg.rightMeanClipPos2; chrname_right = mate_clip_reg.rightClipReg2->chrname; }
	}
	if(chrname_left.size()>0 and chrname_left.compare(chrname_right)==0 and left_loc>right_loc){ // exchange
		reg_tmp = mate_clip_reg.leftClipReg;
		mate_clip_reg.leftClipReg = mate_clip_reg.rightClipReg;
		mate_clip_reg.rightClipReg = reg_tmp;
		reg_tmp = mate_clip_reg.leftClipReg2;
		mate_clip_reg.leftClipReg2 = mate_clip_reg.rightClipReg2;
		mate_clip_reg.rightClipReg2 = reg_tmp;

		tmp = mate_clip_reg.leftClipPosNum;
		mate_clip_reg.leftClipPosNum = mate_clip_reg.rightClipPosNum;
		mate_clip_reg.rightClipPosNum = tmp;
		tmp = mate_clip_reg.leftClipPosNum2;
		mate_clip_reg.leftClipPosNum2 = mate_clip_reg.rightClipPosNum2;
		mate_clip_reg.rightClipPosNum2 = tmp;

		tmp = mate_clip_reg.leftMeanClipPos;
		mate_clip_reg.leftMeanClipPos = mate_clip_reg.rightMeanClipPos;
		mate_clip_reg.rightMeanClipPos = tmp;
		tmp = mate_clip_reg.leftMeanClipPos2;
		mate_clip_reg.leftMeanClipPos2 = mate_clip_reg.rightMeanClipPos2;
		mate_clip_reg.rightMeanClipPos2 = tmp;

		tmp = mate_clip_reg.leftClipRegNum;
		mate_clip_reg.leftClipRegNum = mate_clip_reg.rightClipRegNum;
		mate_clip_reg.rightClipRegNum = tmp;
	}
}


void clipReg::computeVarTypeClipReg(mateClipReg_t &mate_clip_reg, string &inBamFile){
	size_t i, j, var_type, var_type_tmp, dup_type_num, inv_type_num, tra_type_num, dup_num_tmp;
	vector<size_t> dup_num_vec;
	clipAlnData_t *clip_aln_seg, *mate_clip_aln_seg;
	vector<clipAlnData_t*> query_aln_segs;
	string queryname;
	vector<int32_t> adjClipAlnSegInfo;
	clipPos_t clip_pos_item, mate_clip_pos_item;
	int32_t arr_idx, clip_end_flag, mate_clip_end_flag;
	bool same_chr_flag, same_orient_flag, self_overlap_flag, same_aln_reg_flag; // four features
	bool valid_query_end_flag;

	resetClipCheckFlag(clipAlnDataVector); // reset clipping check flag

	dup_type_num = inv_type_num = tra_type_num = 0;
	var_type = VAR_UNC;
	if(mate_clip_reg.reg_mated_flag){

		for(i=0; i<clipAlnDataVector.size(); i++){
			clip_aln_seg = clipAlnDataVector.at(i);
			if(clip_aln_seg->query_checked_flag==false){
				//cout << i << "\t" << clip_aln_seg->chrname << "\t" << clip_aln_seg->startRefPos << "\t" << clip_aln_seg->endRefPos << endl;
				queryname = clipAlnDataVector.at(i)->queryname;
				query_aln_segs = getQueryClipAlnSegs(queryname, clipAlnDataVector);  // get query clip align segments
				if(query_aln_segs.size()==0) continue;

				// sort according to queryPos
				sortQueryAlnSegs(query_aln_segs);

				same_chr_flag = isSameChrome(query_aln_segs);
				same_orient_flag = isSameOrient(query_aln_segs);
				self_overlap_flag = isQuerySelfOverlap(query_aln_segs);
				same_aln_reg_flag = isSameAlnReg(query_aln_segs);

				// predict the variant type
				var_type_tmp = VAR_UNC;
				if(same_chr_flag and same_orient_flag and self_overlap_flag and same_aln_reg_flag){ // DUP
					var_type_tmp = VAR_DUP;
					dup_type_num ++;
				}else if(same_chr_flag and same_orient_flag==false and self_overlap_flag==false and same_aln_reg_flag){ // INV
					var_type_tmp = VAR_INV;
					inv_type_num ++;
				}else if(self_overlap_flag==false and same_aln_reg_flag==false)	{ // TRA
					var_type_tmp = VAR_TRA;
					tra_type_num ++;
				}

				// compute dup_num
				if(var_type_tmp==VAR_DUP){ // DUP
					dup_num_tmp = 0;
					for(j=0; j<query_aln_segs.size(); j++){
						clip_aln_seg = query_aln_segs.at(j);
						valid_query_end_flag = false;
						clip_end_flag = -1;
						if(clip_aln_seg->left_clip_checked_flag==false and clip_aln_seg->leftClipSize>=MIN_CLIP_END_SIZE and clip_aln_seg->startRefPos>=startRefPos and clip_aln_seg->startRefPos<=endRefPos){ // left end
							valid_query_end_flag = true;
							clip_end_flag = LEFT_END;
						}else if(clip_aln_seg->right_clip_checked_flag==false and clip_aln_seg->rightClipSize>=MIN_CLIP_END_SIZE and clip_aln_seg->endRefPos>=startRefPos and clip_aln_seg->endRefPos<=endRefPos){ // right end
							valid_query_end_flag = true;
							clip_end_flag = RIGHT_END;
						}

						if(valid_query_end_flag){
							//cout << "\tseg_len:" << clip_aln_seg->endRefPos - clip_aln_seg->startRefPos << endl;
							// deal with clip end itself
							clip_pos_item.chrname = clip_aln_seg->chrname;
							if(clip_end_flag==LEFT_END){ // left clip end
								clip_pos_item.clipRefPos = clip_aln_seg->startRefPos;
								clip_aln_seg->left_clip_checked_flag = true;
							}else{ // right clip end
								clip_pos_item.clipRefPos = clip_aln_seg->endRefPos;
								clip_aln_seg->right_clip_checked_flag = true;
							}

							// deal with the mate clip end
							adjClipAlnSegInfo = getAdjacentClipAlnSeg(j, clip_end_flag, query_aln_segs);
							arr_idx = adjClipAlnSegInfo.at(0);
							mate_clip_end_flag = adjClipAlnSegInfo.at(1);
							if(arr_idx!=-1){ // mated
								mate_clip_aln_seg = query_aln_segs.at(arr_idx);
								if((mate_clip_end_flag==LEFT_END and mate_clip_aln_seg->left_clip_checked_flag==false) or (mate_clip_end_flag==RIGHT_END and mate_clip_aln_seg->right_clip_checked_flag==false)){
									//cout << "\tmate_seg_len:" << mate_clip_aln_seg->endRefPos - mate_clip_aln_seg->startRefPos << endl;

									mate_clip_pos_item.chrname = mate_clip_aln_seg->chrname;
									if(mate_clip_end_flag==LEFT_END){
										mate_clip_pos_item.clipRefPos = mate_clip_aln_seg->startRefPos;
										mate_clip_aln_seg->left_clip_checked_flag = true;
									}else{
										mate_clip_pos_item.clipRefPos = mate_clip_aln_seg->endRefPos;
										mate_clip_aln_seg->right_clip_checked_flag = true;
									}

									if((clip_end_flag==LEFT_END and mate_clip_end_flag==RIGHT_END and clip_aln_seg->startRefPos<mate_clip_aln_seg->endRefPos)
										or (clip_end_flag==RIGHT_END and mate_clip_end_flag==LEFT_END and clip_aln_seg->endRefPos>mate_clip_aln_seg->startRefPos)){
										dup_num_tmp ++;
										//cout << "\t >>>>>>>>>>>> DUP <<<<<<<<<<<<<" << endl;
										//cout << "\t" << j << ": " << clip_pos_item.clipRefPos << "\t" << mate_clip_pos_item.clipRefPos << endl;
									}
								}
							}
						}
						clip_aln_seg->query_checked_flag = true;
					}

					if(dup_num_tmp>0 and containCompleteDup(query_aln_segs, mate_clip_reg))
						dup_num_vec.push_back(dup_num_tmp);
				}
				for(j=0; j<query_aln_segs.size(); j++) query_aln_segs.at(j)->query_checked_flag = true;
			}
		}
		// extract var_type
		var_type = extractVarType(dup_type_num, inv_type_num, tra_type_num, paras->minClipReadsNumSupportSV);
	}

	mate_clip_reg.sv_type = var_type;
	if(var_type==VAR_DUP) mate_clip_reg.dup_num = computeDupNumClipReg(dup_num_vec); //compute the dup_num
}

void clipReg::resetClipCheckFlag(vector<clipAlnData_t*> &clipAlnDataVector){
	for(size_t i=0; i<clipAlnDataVector.size(); i++){
		clipAlnDataVector.at(i)->left_clip_checked_flag = false;
		clipAlnDataVector.at(i)->right_clip_checked_flag = false;
		clipAlnDataVector.at(i)->query_checked_flag = false;
	}
}

bool clipReg::isSameChrome(vector<clipAlnData_t*> &query_aln_segs){
	bool flag = true;
	string chrname_tmp;
	clipAlnData_t *clip_aln;

	// same chrome
	chrname_tmp = query_aln_segs.at(0)->chrname;
	for(size_t i=1; i<query_aln_segs.size(); i++){
		clip_aln = query_aln_segs.at(i);
		if(clip_aln->chrname.compare(chrname_tmp)!=0){
			flag = false;
			break;
		}
	}

	return flag;
}

bool clipReg::isSameOrient(vector<clipAlnData_t*> &query_aln_segs){
	bool flag = true;
	for(size_t i=1; i<query_aln_segs.size(); i++){
		if(query_aln_segs.at(i)->aln_orient!=query_aln_segs.at(i-1)->aln_orient){
			flag = false;
			break;
		}
	}
	return flag;
}

bool clipReg::isQuerySelfOverlap(vector<clipAlnData_t*> &query_aln_segs){
	bool flag;
	clipAlnData_t *clip_aln1, *clip_aln2;

	flag = false;
	for(size_t i=0; i<query_aln_segs.size()-1; i++){
		clip_aln1 = query_aln_segs.at(i);
		for(size_t j=i+1; j<query_aln_segs.size(); j++){
			clip_aln2 = query_aln_segs.at(j);
			flag = isSegSelfOverlap(clip_aln1, clip_aln2);
			if(flag) break;
		}
		if(flag) break;
	}

	return flag;
}

bool clipReg::isSegSelfOverlap(clipAlnData_t *clip_aln1, clipAlnData_t *clip_aln2){
	bool flag = false;
	int32_t overlap_size;

	if(clip_aln1->chrname.compare(clip_aln2->chrname)==0){
		if(isOverlappedPos(clip_aln1->startRefPos, clip_aln1->endRefPos, clip_aln2->startRefPos, clip_aln2->endRefPos)){
			overlap_size = getOverlapSize(clip_aln1->startRefPos, clip_aln1->endRefPos, clip_aln2->startRefPos, clip_aln2->endRefPos);
			if(overlap_size>=MIN_QUERY_SELF_OVERLAP_SIZE)
				flag = true;
		}else if(clip_aln2->startRefPos<clip_aln1->endRefPos){
			overlap_size = getOverlapSize(clip_aln1->startRefPos, clip_aln1->endRefPos, clip_aln2->startRefPos, clip_aln2->endRefPos);
			if(overlap_size>=-(int32_t)paras->maxClipRegSize)
				flag = true;
		}
	}

	return flag;
}

bool clipReg::isSameAlnReg(vector<clipAlnData_t*> &query_aln_segs){
	size_t i;
	string chrname_tmp;
	bool same_aln_reg_flag, same_chr_flag, adjacent_flag;

	same_aln_reg_flag = true;
	if(query_aln_segs.size()>=2){
		// same chrome
		same_chr_flag = isSameChrome(query_aln_segs);
		if(same_chr_flag==false) same_aln_reg_flag = false;

		// same region in the same chrome
		if(same_aln_reg_flag){
			//same_orient_flag = isSameOrient(query_aln_segs);

			//if(same_orient_flag){ // same orient
				for(i=1; i<query_aln_segs.size(); i++){
					adjacent_flag = isAdjacentClipAlnSeg(query_aln_segs.at(i-1), query_aln_segs.at(i), paras->maxClipRegSize);
					if(adjacent_flag==false){
						same_aln_reg_flag = false;
						break;
					}
				}
//			}else{ // different orient
//				for(i=1; i<query_aln_segs.size(); i++){
//					adjacent_flag = isAdjacentClipAlnSeg(query_aln_segs.at(i-1), query_aln_segs.at(i), paras->maxClipRegSize);
//					if(adjacent_flag==false){
//						same_aln_reg_flag = false;
//						break;
//					}
//				}
//			}
		}
	}

	// check jump segments
	if(same_aln_reg_flag and query_aln_segs.size()>=3){
		for(i=1; i<query_aln_segs.size()-1; i++){
			if(query_aln_segs.at(i)->endRefPos<query_aln_segs.at(i-1)->startRefPos or query_aln_segs.at(i)->startRefPos>query_aln_segs.at(i+1)->endRefPos){
				same_aln_reg_flag = false;
				break;
			}
		}
	}

	return same_aln_reg_flag;
}

void clipReg::sortQueryAlnSegs(vector<clipAlnData_t*> &query_aln_segs){
	size_t i, j, minIdx, minPos, query_pos1, query_pos2;
	clipAlnData_t *clip_aln1, *clip_aln2,  *tmp;

	for(i=0; i<query_aln_segs.size(); i++){
		clip_aln1 = query_aln_segs.at(i);
		if(clip_aln1->aln_orient==ALN_PLUS_ORIENT)
			query_pos1 = clip_aln1->startQueryPos;
		else
			query_pos1 = clip_aln1->querylen - clip_aln1->startQueryPos + 1;
		minIdx = i;
		minPos = query_pos1;
		for(j=i+1; j<query_aln_segs.size(); j++){
			clip_aln2 = query_aln_segs.at(j);
			if(clip_aln2->aln_orient==ALN_PLUS_ORIENT)
				query_pos2 = clip_aln2->startQueryPos;
			else
				query_pos2 = clip_aln2->querylen - clip_aln2->startQueryPos + 1;

			if(query_pos2<minPos){
				minPos = query_pos2;
				minIdx = j;
			}
		}
		if(minIdx!=i){
			tmp = query_aln_segs.at(i);
			query_aln_segs.at(i) = query_aln_segs.at(minIdx);
			query_aln_segs.at(minIdx) = tmp;
		}
	}
}

bool clipReg::isAdjacentClipAlnSeg(clipAlnData_t *clip_aln1, clipAlnData_t *clip_aln2, size_t dist_thres){
	bool flag = false;
	if(clip_aln1->chrname.compare(clip_aln2->chrname)==0)
		flag = isAdjacent(clip_aln1->startRefPos, clip_aln1->endRefPos, clip_aln2->startRefPos, clip_aln2->endRefPos, dist_thres);
	return flag;
}

bool clipReg::containCompleteDup(vector<clipAlnData_t*> &query_aln_segs, mateClipReg_t &mate_clip_reg){
	bool flag, left_end_valid_flag, right_end_valid_flag;
	size_t i;
	clipAlnData_t *clip_aln;

	flag = false;
	left_end_valid_flag = right_end_valid_flag = false;
	for(i=0; i<query_aln_segs.size(); i++){
		clip_aln = query_aln_segs.at(i);
		if(clip_aln->startRefPos<mate_clip_reg.leftMeanClipPos) left_end_valid_flag = true;
		if(clip_aln->endRefPos>mate_clip_reg.rightMeanClipPos) right_end_valid_flag = true;
		if(left_end_valid_flag and right_end_valid_flag){
			flag = true;
			break;
		}
	}

	return flag;
}

size_t clipReg::extractVarType(size_t dup_type_num, size_t inv_type_num, size_t tra_type_num, size_t min_reads_thres){
	size_t var_type, maxValue, maxIdx;

	maxIdx = 0;
	maxValue = dup_type_num;
	if(maxValue<inv_type_num){
		maxValue = inv_type_num;
		maxIdx = 1;
	}

	if(maxValue<tra_type_num){
		maxValue = tra_type_num;
		maxIdx = 2;
	}

	if(maxValue>=min_reads_thres){
		switch(maxIdx){
			case 0: var_type = VAR_DUP; break;
			case 1: var_type = VAR_INV; break;
			case 2: var_type = VAR_TRA; break;
			default: cerr << __func__ << ", line=" << __LINE__ << ": invalid index=" << maxIdx << ", error!" << endl; exit(1);
		}
	}else
		var_type = VAR_UNC;

	//cout << "In " << __func__ << "(): dup_type_num=" << dup_type_num << ", inv_type_num=" << inv_type_num << ", tra_type_num=" << tra_type_num << endl;

	return var_type;
}

size_t clipReg::computeDupNumClipReg(vector<size_t> &dup_num_vec){
	size_t i, dup_num_int, maxValue;

	maxValue = 0;
	for(i=0; i<dup_num_vec.size(); i++){
		if(maxValue<dup_num_vec.at(i)){
			maxValue = dup_num_vec.at(i);
		}
	}
	dup_num_int = maxValue;

	//cout << "leftMeanClipPos=" << mate_clip_reg.leftMeanClipPos << ", rightMeanClipPos=" << mate_clip_reg.rightMeanClipPos << endl;
	//cout << "dup_num_int=" << dup_num_int << endl;

	return dup_num_int;
}

// print clips
void clipReg::printClipVecs(){
	clipPos_t clip_pos;
	double ratio1, ratio2, total;

	cout << "left_mate_clip_pos_vec1: " << leftClipPosVector.size() << endl;
	printSingleClipVec(leftClipPosVector);

	cout << "left_mate_clip_pos_vec2: " << leftClipPosVector2.size() << endl;
	printSingleClipVec(leftClipPosVector2);

	cout << "right_mate_clip_pos_vec1: " << rightClipPosVector.size() << endl;
	printSingleClipVec(rightClipPosVector);

	cout << "right_mate_clip_pos_vec2: " << rightClipPosVector2.size() << endl;
	printSingleClipVec(rightClipPosVector2);

	ratio1 = -1, ratio2 = -1, total = leftClipPosVector.size() + rightClipPosVector.size();
	if(total>0) ratio1 = (double)leftClipPosVector.size() / total;
	if(total>0) ratio2 = (double)rightClipPosVector.size() / total;
	cout << "ratio1=" << ratio1 << "; " << "ratio2=" << ratio2 << endl;
}

void clipReg::printSingleClipVec(vector<clipPos_t> &clip_pos_vec){
	clipPos_t clip_pos;
	for(size_t i=0; i<clip_pos_vec.size(); i++){
		clip_pos = clip_pos_vec.at(i);
		cout << "\t" << clip_pos.chrname << ": " << clip_pos.clipRefPos << ", same_orient=" << clip_pos.same_orient_flag << endl;
	}
}

// print clip regions
void clipReg::printResultClipRegs(){
	reg_t *reg;
	cout << "The retult clip regions:" << endl;
	reg = mate_clip_reg.leftClipReg;
	if(reg){
		cout << "The left clip region 1:" << endl;
		cout << "\t" << reg->chrname << ": " << reg->startRefPos << "-" << reg->endRefPos << endl;
	}
	reg = mate_clip_reg.leftClipReg2;
	if(reg){
		cout << "The left clip region 2:" << endl;
		cout << "\t" << reg->chrname << ": " << reg->startRefPos << "-" << reg->endRefPos << endl;
	}
	reg = mate_clip_reg.rightClipReg;
	if(reg){
		cout << "The right clip region 1:" << endl;
		cout << "\t" << reg->chrname << ": " << reg->startRefPos << "-" << reg->endRefPos << endl;
	}
	reg = mate_clip_reg.rightClipReg2;
	if(reg){
		cout << "The right clip region 2:" << endl;
		cout << "\t" << reg->chrname << ": " << reg->startRefPos << "-" << reg->endRefPos << endl;
	}
}
