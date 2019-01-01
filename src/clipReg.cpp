#include <math.h>
#include "clipReg.h"
#include "util.h"


clipReg::clipReg(string &chrname, size_t startRefPos, size_t endRefPos, size_t chrlen, string &inBamFile, faidx_t *fai){
	this->chrname = chrname;
	this->startRefPos = startRefPos - CLIP_END_EXTEND_SIZE;
	this->endRefPos = endRefPos + CLIP_END_EXTEND_SIZE;
	this->chrlen = chrlen;
	this->inBamFile = inBamFile;
	this->fai = fai;

	if(this->startRefPos<1) this->startRefPos = 1;
	if(this->endRefPos>chrlen) this->endRefPos = chrlen;

	mate_clip_reg.leftClipReg = mate_clip_reg.rightClipReg = NULL;
	mate_clip_reg.reg_mated_flag = false;
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
	clip_aln_data_loader.loadClipAlnDataWithSATag(clipAlnDataVector);
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

	extractClipPosVec();  // get the clip pos vectorascending

	sortClipPos(); // sort clips

	//printClipVec();  // print clips

	removeFakeClips();  // remove fake clips

	//printClipVec();  // print clips

	computeClipRegs();  // get the mediate and the around region

	//removeFalseOverlappedMateClipReg();

	//printResultClipRegs();
}


// compute the mate clip region
void clipReg::extractClipPosVec(){
	size_t i, j;
	vector<clipAlnData_t*> query_aln_segs;
	clipAlnData_t *clip_aln_seg;
	string queryname, clip_pos_str;
	vector<int32_t> adjClipAlnSegInfo;
	clipPos_t clip_pos_item, mate_clip_pos_item;
	int32_t arr_idx, clip_end_flag, mate_clip_end_flag, clip_vec_idx, mate_clip_vec_idx;
	bool valid_query_end_flag, valid_mate_query_end_flag, same_orient_flag;

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
						if((mate_clip_end_flag==LEFT_END and query_aln_segs.at(arr_idx)->left_clip_checked_flag==false) or (mate_clip_end_flag==RIGHT_END and query_aln_segs.at(arr_idx)->right_clip_checked_flag==false)){
							//cout << "\tmate_seg_len:" << query_aln_segs.at(arr_idx)->endRefPos - query_aln_segs.at(arr_idx)->startRefPos << endl;

							mate_clip_pos_item.chrname = query_aln_segs.at(arr_idx)->chrname;
							if(mate_clip_end_flag==LEFT_END){
								mate_clip_pos_item.clipRefPos = query_aln_segs.at(arr_idx)->startRefPos;
								query_aln_segs.at(arr_idx)->left_clip_checked_flag = true;
								mate_clip_vec_idx = 0;
							}else{
								mate_clip_pos_item.clipRefPos = query_aln_segs.at(arr_idx)->endRefPos;
								query_aln_segs.at(arr_idx)->right_clip_checked_flag = true;
								mate_clip_vec_idx = 1;
							}
							valid_mate_query_end_flag = true;

							// determine which is the left and which is the right clip region for inversions
							if(clip_aln_seg->aln_orient!=query_aln_segs.at(arr_idx)->aln_orient){
								same_orient_flag = false;
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

void clipReg::sortClipPos(){
	sortClipPosSingleVec(leftClipPosVector);
	sortClipPosSingleVec(rightClipPosVector);
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
	removeFakeClipsDifferentChrSingleVec(rightClipPosVector);

	removeFakeClipsLongDistSingleVec(leftClipPosVector);
	removeFakeClipsLongDistSingleVec(rightClipPosVector);
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
void clipReg::removeFakeClipsLongDistSingleVec(vector<clipPos_t> &clipPosVector){
	size_t i, j;
	int32_t idx, dist;
	vector<size_t> clip_pos_num_vec;
	vector<clipPos_t> clip_pos_vec;
	clipPos_t clip_pos_item, clip_pos_item_tmp;

	size_t minPos, minIdx, maxValue, num, tmp, num0;
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

	// print
//	for(i=0; i<clip_pos_vec.size(); i++){
//		clip_pos_item = clip_pos_vec.at(i);
//		cout << "\t" << clip_pos_item.chrname << ":" << clip_pos_item.clipRefPos << "\t" << clip_pos_num_vec.at(i) << endl;
//	}

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
	num0 = clipPosVector.size();
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

//	double ratio = (double)clipPosVector.size() / num0;
//	if(ratio<MIN_VALID_CLIP_POS_RATIO)
//		clipPosVector.clear();
//	cout << ratio << endl;
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
	mate_clip_reg.valid_flag = true;
	if(mate_clip_reg.leftClipReg and mate_clip_reg.rightClipReg) mate_clip_reg.reg_mated_flag = true;

	// sv_type and dup_num
	mate_clip_reg.sv_type = VAR_UNC;
	mate_clip_reg.dup_num = 0;
	if(mate_clip_reg.reg_mated_flag) // sv_type
		computeVarTypeClipReg(mate_clip_reg, inBamFile);
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

// compute the mean size of the clippings
size_t clipReg::computeMeanClipPos(vector<clipPos_t> &clipPosVector){
	size_t i, total;
	total = 0;
	for(i=0; i<clipPosVector.size(); i++) total += clipPosVector.at(i).clipRefPos;
	return round((double)total/clipPosVector.size());
}


void clipReg::computeVarTypeClipReg(mateClipReg_t &mate_clip_reg, string &inBamFile){
	size_t i, j, var_type, dup_type_num, inv_type_num, tra_type_num, dup_num_tmp;
	vector<size_t> dup_num_vec;
	clipAlnData_t *clip_aln_seg, *mate_clip_aln_seg;
	vector<clipAlnData_t*> query_aln_segs;
	string queryname;
	vector<int32_t> adjClipAlnSegInfo;
	clipPos_t clip_pos_item, mate_clip_pos_item;
	int32_t arr_idx, clip_end_flag, mate_clip_end_flag, mate_clip_vec_idx;
	bool valid_query_end_flag;

	resetClipCheckFlag(clipAlnDataVector); // reset clipping check flag

	dup_type_num = inv_type_num = tra_type_num = 0;
	var_type = VAR_UNC;
	if(mate_clip_reg.reg_mated_flag){
		if(mate_clip_reg.leftClipReg->chrname.compare(mate_clip_reg.rightClipReg->chrname)==0){ // same chrome: DUP, INV or TRA
			for(i=0; i<clipAlnDataVector.size(); i++){
				clip_aln_seg = clipAlnDataVector.at(i);
				//cout << "\t" << clip_aln_seg->chrname << "\t" << clip_aln_seg->startRefPos << "\t" << clip_aln_seg->endRefPos << endl;

				queryname = clipAlnDataVector.at(i)->queryname;
				query_aln_segs = getQueryClipAlnSegs(queryname, clipAlnDataVector);  // get query clip align segments

				if(isSameOrient(query_aln_segs)){ // DUP or TRA


					dup_type_num ++;
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

						mate_clip_vec_idx = -1;
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
									//cout << "\tmate_seg_len:" << query_aln_segs.at(arr_idx)->endRefPos - query_aln_segs.at(arr_idx)->startRefPos << endl;

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
				}else{ // INV
					inv_type_num ++;
				}
			}
		}else{ // different chrome, TRA
			tra_type_num ++;
		}
		// extract var_type
		var_type = extractVarType(dup_type_num, inv_type_num, tra_type_num);
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

size_t clipReg::extractVarType(size_t dup_type_num, size_t inv_type_num, size_t tra_type_num){
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

	switch(maxIdx){
		case 0: var_type = VAR_DUP; break;
		case 1: var_type = VAR_INV; break;
		case 2: var_type = VAR_TRA; break;
		default: cerr << __func__ << ", line=" << __LINE__ << ": invalid index=" << maxIdx << ", error!" << endl; exit(1);
	}

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
void clipReg::printClipVec(){
	size_t i;
	double ratio1, ratio2, total;

	// print the mate clip position
	cout << "left_mate_clip_pos_vec: " << leftClipPosVector.size() << endl;
	for(i=0; i<leftClipPosVector.size(); i++)
		cout << "\t" << leftClipPosVector.at(i).chrname << ": " << leftClipPosVector.at(i).clipRefPos << ", same_orient=" << leftClipPosVector.at(i).same_orient_flag << endl;

	cout << "right_mate_clip_pos_vec: " << rightClipPosVector.size() << endl;
	for(i=0; i<rightClipPosVector.size(); i++)
		cout << "\t" << rightClipPosVector.at(i).chrname << ": " << rightClipPosVector.at(i).clipRefPos << ", same_orient=" << rightClipPosVector.at(i).same_orient_flag << endl;

	ratio1 = -1, ratio2 = -1, total = leftClipPosVector.size() + rightClipPosVector.size();
	if(total>0) ratio1 = (double)leftClipPosVector.size() / total;
	if(total>0) ratio2 = (double)rightClipPosVector.size() / total;
	cout << "ratio1=" << ratio1 << "; " << "ratio2=" << ratio2 << endl;
}

// print clip regions
void clipReg::printResultClipRegs(){
	reg_t *reg;
	cout << "The retult clip regions:" << endl;
	reg = mate_clip_reg.leftClipReg;
	if(reg){
		cout << "The left clip region:" << endl;
		cout << "\t" << reg->chrname << ": " << reg->startRefPos << "-" << reg->endRefPos << endl;
	}
	reg = mate_clip_reg.rightClipReg;
	if(reg){
		cout << "The right clip region:" << endl;
		cout << "\t" << reg->chrname << ": " << reg->startRefPos << "-" << reg->endRefPos << endl;
	}
}
