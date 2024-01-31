#include <math.h>
#include "clipReg.h"
#include "clipAlnDataLoader.h"
#include "util.h"

extern pthread_mutex_t mutex_fai;

clipReg::clipReg(string &chrname, int64_t startRefPos, int64_t endRefPos, string &inBamFile, faidx_t *fai, Paras *paras){
	this->chrname = chrname;
	this->startRefPos = startRefPos - CLIP_END_EXTEND_SIZE;
	this->endRefPos = endRefPos + CLIP_END_EXTEND_SIZE;
	this->chrlen = faidx_seq_len(fai, chrname.c_str()); // get the reference length;
	this->inBamFile = inBamFile;
	this->fai = fai;
	this->paras = paras;
	this->minClipEndSize = paras->minClipEndSize;
	maxVarRegSize = paras->maxVarRegSize;
	minAlnSize_same_orient = MIN_ALN_SIZE_SAME_ORIENT;
	minAlnSize_diff_orient = MIN_ALN_SIZE_DIFF_ORIENT;

	if(this->startRefPos<1) this->startRefPos = 1;
	if(this->endRefPos>chrlen) this->endRefPos = chrlen;
	left_part_changed = right_part_changed = false;

	largeIndelClipReg = NULL;
	large_indel_flag = false;
	supp_num_largeIndel = depth_largeIndel = 0;

	mate_clip_reg.leftClipReg = mate_clip_reg.leftClipReg2 = mate_clip_reg.rightClipReg = mate_clip_reg.rightClipReg2 = NULL;
	mate_clip_reg.leftClipRegNum = mate_clip_reg.rightClipRegNum = 0;
	mate_clip_reg.leftMeanClipPos = mate_clip_reg.rightMeanClipPos = mate_clip_reg.leftMeanClipPos2 = mate_clip_reg.rightMeanClipPos2 = 0;
	mate_clip_reg.leftClipPosNum = mate_clip_reg.rightClipPosNum = mate_clip_reg.leftClipPosNum2 = mate_clip_reg.rightClipPosNum2 = 0;
	mate_clip_reg.reg_mated_flag = false;
	mate_clip_reg.valid_flag = false;
	mate_clip_reg.sv_type = VAR_UNC;
	mate_clip_reg.dup_num = 0;
	mate_clip_reg.supp_num_valid_flag = false;
	mate_clip_reg.call_success_flag = false;
	mate_clip_reg.large_indel_flag = false;
	mate_clip_reg.largeIndelClipReg = NULL;

	mate_clip_reg.var_cand = mate_clip_reg.left_var_cand_tra = mate_clip_reg.right_var_cand_tra = NULL;
	mate_clip_reg.leftClipPosTra1 = mate_clip_reg.rightClipPosTra1 = mate_clip_reg.leftClipPosTra2 = mate_clip_reg.rightClipPosTra2 = 0;
}

clipReg::~clipReg() {
	if(!clipAlnDataVector.empty()) destroyClipAlnDataVector(clipAlnDataVector);
	if(!clipAlnDataVector2.empty()) destroyClipAlnDataVector(clipAlnDataVector2);
	if(!rightClipAlnDataVector.empty()) destroyClipAlnDataVector(rightClipAlnDataVector);
	if(!rightClipAlnDataVector2.empty()) destroyClipAlnDataVector(rightClipAlnDataVector2);
	if(!leftClipPosVector.empty()) destroyClipPosVec(leftClipPosVector);
	if(!leftClipPosVector2.empty()) destroyClipPosVec(leftClipPosVector2);
	if(!rightClipPosVector.empty()) destroyClipPosVec(rightClipPosVector);
	if(!rightClipPosVector2.empty()) destroyClipPosVec(rightClipPosVector2);
	if(!largeIndelClipPosVector.empty()) destroyClipPosVec(largeIndelClipPosVector);
	if(large_indel_flag) delete largeIndelClipReg;
}

void clipReg::destroyClipAlnDataVector(vector<clipAlnData_t*> &clipAlnDataVec){
	for(size_t i=0; i<clipAlnDataVec.size(); i++){
		bam_destroy1(clipAlnDataVec.at(i)->bam);
		delete clipAlnDataVec.at(i);
	}
	vector<clipAlnData_t*>().swap(clipAlnDataVec);
}

// destroy clipping position vector
void clipReg::destroyClipPosVec(vector<clipPos_t*> &clipPosVec){
	for(size_t i=0; i<clipPosVec.size(); i++) delete clipPosVec.at(i);
	vector<clipPos_t*>().swap(clipPosVec);
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
	clipAlnDataLoader clip_aln_data_loader(chrname, startRefPos, endRefPos, inBamFile, minClipEndSize);
	clip_aln_data_loader.loadClipAlnDataWithSATag(clipAlnDataVector, paras->max_ultra_high_cov); // removed 2023-12-08
//	clip_aln_data_loader.loadClipAlnDataWithSATagWithSegSize(clipAlnDataVector, paras->max_ultra_high_cov, paras->max_seg_size_ratio_usr); // modified 2023-12-08
}

void clipReg::removeNonclipItems(){
	removeNonclipItemsOp(clipAlnDataVector);
}

// remove query having no clippings
void clipReg::removeNonclipItemsOp(vector<clipAlnData_t*> &clipAlnDataVector){
	clipAlnData_t *clip_aln;
	for(size_t i=0; i<clipAlnDataVector.size(); ){
		clip_aln = clipAlnDataVector.at(i);
		if(clip_aln->leftClipSize<minClipEndSize and clip_aln->rightClipSize<minClipEndSize){
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

	//if(seg_num_non_SA>0 and longer_seg_num_non_SA>=LONGER_NON_SA_SEG_NUM_THRES and (double)longer_seg_num_non_SA/seg_num_non_SA>=LONGER_NON_SA_SEG_RATIO_THRES) valid_flag = true;
	if(longer_seg_num_non_SA>=paras->minReadsNumSupportSV*READS_NUM_SUPPORT_FACTOR and (double)longer_seg_num_non_SA/seg_num_non_SA>=LONGER_NON_SA_SEG_RATIO_THRES) valid_flag = true;
	if(clipAlnDataVector.size()<=10 and seg_num_non_SA>0 and (double)longer_seg_num_non_SA/seg_num_non_SA>=LONGER_NON_SA_SEG_RATIO_THRES) valid_flag = true;//The number of CliPVs is small
	//double ratio = (double)longer_seg_num_non_SA / seg_num_non_SA;
	//cout << "longer_non_seg_ratio=" << ratio << ", valid_flag=" << valid_flag << endl;

	return valid_flag;
}

// compute the mate clip region
void clipReg::computeMateAlnClipReg(){

	extractClipPosVec();  // get the clip pos vector

	//printClipVecs("After extract");  // print clips

	splitClipPosVec(); // split vector

	//printClipVecs("After split");  // print clips

	appendClipPos();  // append clipping positions

	sortClipPos(); // sort clips

	//printClipVecs("After sort");  // print clips

	removeFakeClips();  // remove fake clips

	//printClipVecs("After remove fake clips");  // print clips

	computeClipRegs();  // get the mediate and the around region

	//printResultClipRegs(); // print clip regions

	//cout << " ---- " << endl;
}

// compute the mate clip region
void clipReg::extractClipPosVec(){
	size_t i, j;
	vector<clipAlnData_t*> query_aln_segs;
	clipAlnData_t *clip_aln_seg, *mate_clip_aln_seg;
	string queryname, clip_pos_str;
	vector<int32_t> adjClipAlnSegInfo;
	clipPos_t *clip_pos_item, *mate_clip_pos_item, *clip_pos_tmp;
	int32_t mate_arr_idx, clip_end_flag, mate_clip_end_flag, clip_vec_idx, mate_clip_vec_idx, mate_clip_vec_idx_tmp, dist;
	int32_t mate_arr_idx_same_chr, mate_clip_end_flag_same_chr, min_dist_same_chr, min_ref_dist_same_chr, seg_skip_num, ref_skip_num;
	bool valid_query_end_flag, valid_mate_query_end_flag, same_orient_flag, self_overlap_flag, vec_id_changed_flag, skip_flag;

	resetClipCheckFlag(clipAlnDataVector); // reset clipping check flag

	seg_skip_num = ref_skip_num = 0;
	for(i=0; i<clipAlnDataVector.size(); i++){
		if(clipAlnDataVector.at(i)->query_checked_flag==false){
			queryname = clipAlnDataVector.at(i)->queryname;
			query_aln_segs = getQueryClipAlnSegsAll(queryname, clipAlnDataVector);  // get query clip align segments

			if(query_aln_segs.size()>MAX_ALN_SEG_NUM_PER_READ_TRA) { // ignore reads of too many align segments
				//cout << "clipReg: " << chrname << ":" << startRefPos << "-" << endRefPos << ", qname=" << queryname << ", align segment number=" << query_aln_segs.size() << endl;
				continue;
			}

//			if(queryname.compare("SRR8858449.1.125122")==0){
//				cout << queryname << endl;
//			}

			for(j=0; j<query_aln_segs.size(); j++){
				clip_aln_seg = query_aln_segs.at(j);
				valid_query_end_flag = false;
				clip_end_flag = -1;
				if(clip_aln_seg->left_clip_checked_flag==false and clip_aln_seg->leftClipSize>=minClipEndSize and clip_aln_seg->startRefPos>=startRefPos and clip_aln_seg->startRefPos<=endRefPos and clip_aln_seg->chrname.compare(chrname)==0){ // left end
					valid_query_end_flag = true;
					clip_end_flag = LEFT_END;
				}else if(clip_aln_seg->right_clip_checked_flag==false and clip_aln_seg->rightClipSize>=minClipEndSize and clip_aln_seg->endRefPos>=startRefPos and clip_aln_seg->endRefPos<=endRefPos and clip_aln_seg->chrname.compare(chrname)==0){ // right end
					valid_query_end_flag = true;
					clip_end_flag = RIGHT_END;
				}

				clip_pos_item = mate_clip_pos_item = NULL;
				valid_mate_query_end_flag = false;
				clip_vec_idx = mate_clip_vec_idx = -1;
				same_orient_flag = true;
				if(valid_query_end_flag){
					//cout << "\tseg_len:" << clip_aln_seg->endRefPos - clip_aln_seg->startRefPos << endl;
					// deal with clip end itself
					clip_pos_item = new clipPos_t();
					clip_pos_item->chrname = clip_aln_seg->chrname;
					clip_pos_item->clip_end = clip_end_flag;
					clip_pos_item->clip_aln = clip_aln_seg;
					if(clip_end_flag==LEFT_END){ // left clip end
						clip_pos_item->clipRefPos = clip_aln_seg->startRefPos;
						clip_pos_item->clipQueryPos = clip_aln_seg->startQueryPos;
						clip_pos_item->aln_orient = clip_aln_seg->aln_orient;
						clip_aln_seg->left_clip_checked_flag = true;
						clip_vec_idx = 0;
					}else{ // right clip end
						clip_pos_item->clipRefPos = clip_aln_seg->endRefPos;
						clip_pos_item->clipQueryPos = clip_aln_seg->endQueryPos;
						clip_pos_item->aln_orient = clip_aln_seg->aln_orient;
						clip_aln_seg->right_clip_checked_flag = true;
						clip_vec_idx = 1;
					}

					// deal with the mate clip end
					skip_flag = false;
					adjClipAlnSegInfo = getAdjacentClipAlnSeg(j, clip_end_flag, query_aln_segs, minClipEndSize, MAX_VAR_REG_SIZE);
					mate_arr_idx = adjClipAlnSegInfo.at(0);
					mate_clip_end_flag = adjClipAlnSegInfo.at(1);
					mate_arr_idx_same_chr = adjClipAlnSegInfo.at(4);
					mate_clip_end_flag_same_chr = adjClipAlnSegInfo.at(5);
					min_dist_same_chr = adjClipAlnSegInfo.at(6);
					min_ref_dist_same_chr = adjClipAlnSegInfo.at(7);
					if(mate_arr_idx!=-1){ // mated
						if(mate_arr_idx_same_chr!=-1){
							mate_clip_aln_seg = query_aln_segs.at(mate_arr_idx_same_chr);
							self_overlap_flag = isSegSelfOverlap(clip_aln_seg, mate_clip_aln_seg, maxVarRegSize);
							// prefer the segments on the same chromosome
							if(mate_arr_idx!=mate_arr_idx_same_chr and (abs(min_ref_dist_same_chr)<=MAX_REF_DIST_SAME_CHR or self_overlap_flag)){ // different chromosomes with near reference distance
							//if(abs(min_ref_dist_same_chr)<=MAX_REF_DIST_SAME_CHR){ // different chromosomes with near reference distance
								mate_arr_idx = mate_arr_idx_same_chr;
								mate_clip_end_flag = mate_clip_end_flag_same_chr;
								seg_skip_num ++;
								skip_flag = true;
							}else if(mate_arr_idx==mate_arr_idx_same_chr and mate_clip_end_flag==mate_clip_end_flag_same_chr and self_overlap_flag==false and abs(min_dist_same_chr)<=MAX_REF_DIST_SAME_CHR){
								ref_skip_num ++;
								skip_flag = true;
							}
						}

						mate_clip_aln_seg = query_aln_segs.at(mate_arr_idx);
						if((mate_clip_end_flag==LEFT_END and mate_clip_aln_seg->left_clip_checked_flag==false and mate_clip_aln_seg->chrname.compare(chrname)==0) or (mate_clip_end_flag==RIGHT_END and mate_clip_aln_seg->right_clip_checked_flag==false and mate_clip_aln_seg->chrname.compare(chrname)==0)){//
							//cout << "\tmate_seg_len:" << mate_clip_aln_seg->endRefPos - mate_clip_aln_seg->startRefPos << endl;
							mate_clip_pos_item = new clipPos_t();
							mate_clip_pos_item->chrname = mate_clip_aln_seg->chrname;
							mate_clip_pos_item->clip_aln = mate_clip_aln_seg;
							mate_clip_pos_item->clip_end = mate_clip_end_flag;
							if(mate_clip_end_flag==LEFT_END){
								mate_clip_pos_item->clipRefPos = mate_clip_aln_seg->startRefPos;
								mate_clip_pos_item->clipQueryPos = mate_clip_aln_seg->startQueryPos;
								mate_clip_pos_item->aln_orient = mate_clip_aln_seg->aln_orient;
								mate_clip_aln_seg->left_clip_checked_flag = true;
								mate_clip_vec_idx = 0;
							}else{
								mate_clip_pos_item->clipRefPos = mate_clip_aln_seg->endRefPos;
								mate_clip_pos_item->clipQueryPos = mate_clip_aln_seg->endQueryPos;
								mate_clip_pos_item->aln_orient = mate_clip_aln_seg->aln_orient;
								mate_clip_aln_seg->right_clip_checked_flag = true;
								mate_clip_vec_idx = 1;
							}
							valid_mate_query_end_flag = true;

							if(clip_aln_seg->aln_orient!=mate_clip_aln_seg->aln_orient) same_orient_flag = false;

							// determine which is the left and which is the right clip region for INV and TRA
							self_overlap_flag = isSegSelfOverlap(clip_aln_seg, mate_clip_aln_seg, maxVarRegSize);
							if(self_overlap_flag==false){
								if(clip_vec_idx==mate_clip_vec_idx){ // deal with same part
									mate_clip_vec_idx_tmp = getMinDistVecId(mate_clip_pos_item);
									if(mate_clip_vec_idx_tmp!=-1){
										mate_clip_vec_idx = mate_clip_vec_idx_tmp;
										if(mate_clip_vec_idx==2 or mate_clip_vec_idx==3) clip_vec_idx = 0;
									}else if(same_orient_flag==false){ // different orientation
										clip_vec_idx = 0;
										mate_clip_vec_idx = 1;
									}else
										valid_mate_query_end_flag = false;
								}else if(clip_pos_item->chrname.compare(mate_clip_pos_item->chrname)!=0 or clip_pos_item->clipRefPos<=mate_clip_pos_item->clipRefPos){
									clip_vec_idx = 0;
									mate_clip_vec_idx = 1;
								}else{
									clip_vec_idx = 1;
									mate_clip_vec_idx = 0;
								}
							}
							clip_pos_item->self_overlap_flag = mate_clip_pos_item->self_overlap_flag = self_overlap_flag;
						}
					}

					// adjust the vector id
					if(valid_query_end_flag and valid_mate_query_end_flag){
						vec_id_changed_flag = false;
						if(clip_pos_item->chrname.compare(mate_clip_pos_item->chrname)==0){
							dist = clip_pos_item->clipRefPos - mate_clip_pos_item->clipRefPos;
							if(dist<0) dist = -dist;
							if(dist<MAX_CLIP_REG_MERGE_DIST or (skip_flag and (dist<MAX_REF_DIST_SAME_CHR or self_overlap_flag))){
								clip_vec_idx = mate_clip_vec_idx = 0;
								vec_id_changed_flag = true;
							}
						}
						if(vec_id_changed_flag==false){
							if(clip_pos_item->chrname.compare(mate_clip_pos_item->chrname)!=0 or (mate_clip_pos_item->clip_end==RIGHT_END and (mate_clip_pos_item->clip_aln->left_aln==NULL or mate_clip_pos_item->clip_aln->left_aln->chrname.compare(clip_pos_item->chrname)!=0 or isOverlappedPos(mate_clip_pos_item->clip_aln->left_aln->startRefPos, mate_clip_pos_item->clip_aln->left_aln->endRefPos, startRefPos, startRefPos)==false))
							  or (mate_clip_pos_item->clip_end==LEFT_END and (mate_clip_pos_item->clip_aln->right_aln==NULL or mate_clip_pos_item->clip_aln->right_aln->chrname.compare(clip_pos_item->chrname)!=0 or isOverlappedPos(mate_clip_pos_item->clip_aln->right_aln->startRefPos, mate_clip_pos_item->clip_aln->right_aln->endRefPos, startRefPos, startRefPos)==false))){
								mate_clip_vec_idx_tmp = getMinDistVecId(mate_clip_pos_item);
								if(mate_clip_vec_idx_tmp!=-1)
									mate_clip_vec_idx = mate_clip_vec_idx_tmp;
							}
						}
//						if((same_orient_flag==false and dist>10000)or((clip_pos_item->chrname.compare(mate_clip_pos_item->chrname)!=0)and same_orient_flag==false)) valid_mate_query_end_flag = false;
					}

					// save to vector
					if(valid_query_end_flag){
						clip_pos_item->same_orient_flag = same_orient_flag;
						if(clip_vec_idx==0){
							// two clipping ends of the same align segments should be in different vector
							clip_pos_tmp = getClipPosSameAlnSeg(clip_pos_item, leftClipPosVector);
							if(clip_pos_tmp)
								leftClipPosVector2.push_back(clip_pos_item);
							else
								leftClipPosVector.push_back(clip_pos_item);
						}else if(clip_vec_idx==1){
							if(same_orient_flag){
								clip_pos_tmp = getClipPosSameAlnSeg(clip_pos_item, rightClipPosVector);
								if(clip_pos_tmp)
									rightClipPosVector2.push_back(clip_pos_item);
								else
									rightClipPosVector.push_back(clip_pos_item);
							}else{
								clip_pos_tmp = getClipPosSameAlnSeg(clip_pos_item, rightClipPosVector2);
								if(clip_pos_tmp)
									rightClipPosVector.push_back(clip_pos_item);
								else
									rightClipPosVector2.push_back(clip_pos_item);
							}
						}
					}else{
						delete clip_pos_item;
						clip_pos_item = NULL;
					}
					if(valid_mate_query_end_flag){
						mate_clip_pos_item->same_orient_flag = same_orient_flag;
						if(mate_clip_vec_idx==0){
							// two clipping ends of the same align segments should be in different vector
							clip_pos_tmp = getClipPosSameAlnSeg(mate_clip_pos_item, leftClipPosVector);
							if(clip_pos_tmp)
								leftClipPosVector2.push_back(mate_clip_pos_item);
							else
								leftClipPosVector.push_back(mate_clip_pos_item);
						}else if(mate_clip_vec_idx==1 or mate_clip_vec_idx==2){
							if(same_orient_flag) {
								clip_pos_tmp = getClipPosSameAlnSeg(mate_clip_pos_item, rightClipPosVector);
								if(clip_pos_tmp)
									rightClipPosVector2.push_back(mate_clip_pos_item);
								else
									rightClipPosVector.push_back(mate_clip_pos_item);
							}else{
								clip_pos_tmp = getClipPosSameAlnSeg(mate_clip_pos_item, rightClipPosVector2);
								if(clip_pos_tmp)
									rightClipPosVector.push_back(mate_clip_pos_item);
								else
									rightClipPosVector2.push_back(mate_clip_pos_item);
							}
							//rightClipPosVector.push_back(mate_clip_pos_item); // for mate_clip_vec_idx==2
						}else if(mate_clip_vec_idx==3)
							rightClipPosVector2.push_back(mate_clip_pos_item);
					}else{
						delete mate_clip_pos_item;
						mate_clip_pos_item = NULL;
					}
				}
			}
			for(j=0; j<query_aln_segs.size(); j++) query_aln_segs.at(j)->query_checked_flag = true;
		}
	}

	resetClipCheckFlag(clipAlnDataVector); // reset clipping check flag

	//cout << "seg_skip_num=" << seg_skip_num << ", ref_skip_num=" << ref_skip_num << endl;

	// adjust the single breakpoint
	if(seg_skip_num>=paras->minReadsNumSupportSV*READS_NUM_SUPPORT_FACTOR){
		if(rightClipPosVector2.size()>0) AdjustClipPosVecByBPConsistency(rightClipPosVector2, 3);
		if(rightClipPosVector.size()>0) AdjustClipPosVecByBPConsistency(rightClipPosVector, 2);
		if(leftClipPosVector2.size()>0) AdjustClipPosVecByBPConsistency(leftClipPosVector2, 1);
		if(leftClipPosVector.size()>0) AdjustClipPosVecByBPConsistency(leftClipPosVector, 0);
	}

	if(seg_skip_num>=paras->minReadsNumSupportSV or ref_skip_num>=paras->minReadsNumSupportSV){
		removeFakeClipsLowCov(leftClipPosVector, paras->minReadsNumSupportSV);
		removeFakeClipsLowCov(leftClipPosVector2, paras->minReadsNumSupportSV);
		removeFakeClipsLowCov(rightClipPosVector, paras->minReadsNumSupportSV);
		removeFakeClipsLowCov(rightClipPosVector2, paras->minReadsNumSupportSV);
	}
}

// get the clip position vector id of minimal distance
int32_t clipReg::getMinDistVecId(clipPos_t *clip_pos_item){
	int64_t i, dist_arr[4], minIdx, minValue;

	dist_arr[0] = getMeanDistSingleVec(clip_pos_item, leftClipPosVector);
	dist_arr[1] = getMeanDistSingleVec(clip_pos_item, leftClipPosVector2);
	dist_arr[2] = getMeanDistSingleVec(clip_pos_item, rightClipPosVector);
	dist_arr[3] = getMeanDistSingleVec(clip_pos_item, rightClipPosVector2);

	minIdx = -1;
	minValue = INT_MAX;
	for(i=0; i<4; i++){
		if(dist_arr[i]!=-1){
			if(minValue>dist_arr[i]){
				minValue = dist_arr[i];
				minIdx = i;
			}
		}
	}

	//if(minValue>MIN_CLIP_DIST_THRES){
	if(minValue>MAX_REF_DIST_SAME_CHR){
		//cout << __func__ << ": " << "minIdx=" << minIdx << ", minValue=" << minValue << ", invalid and is skipped!" << endl;
		minIdx = -1;
	}

	return minIdx;
}

// get the mean distance
int32_t clipReg::getMeanDistSingleVec(clipPos_t *clip_pos_item, vector<clipPos_t*> &clip_pos_vec){
	int64_t dist, num, sum;
	clipPos_t *clip_pos_tmp;

	sum = num = 0;
	for(size_t i=0; i<clip_pos_vec.size(); i++){
		clip_pos_tmp = clip_pos_vec.at(i);
		if(clip_pos_tmp->chrname.compare(clip_pos_item->chrname)==0){
			dist = clip_pos_item->clipRefPos - clip_pos_tmp->clipRefPos;
			if(dist<0) dist = -dist;
			sum += dist;
			num ++;
		}
	}

	if(num==0) return -1;
	else return sum / num;
}


// get the clipping position item with the same align segment from a vector
clipPos_t* clipReg::getClipPosSameAlnSeg(clipPos_t *clip_pos_item, vector<clipPos_t*> &clip_pos_vec){
	clipPos_t *clip_pos_ret = NULL, *clip_pos;

	if(clip_pos_item){
		for(size_t i=0; i<clip_pos_vec.size(); i++){
			clip_pos = clip_pos_vec.at(i);
			if(clip_pos->clip_aln==clip_pos_item->clip_aln){
				clip_pos_ret = clip_pos;
				break;
			}
		}
	}

	return clip_pos_ret;
}

void clipReg::splitClipPosVec(){
	size_t i;
	clipPos_t *clip_pos, *clip_pos2, *clip_pos3, *min_dist_clip_pos, *clip_pos_right_clip_end, *clip_pos_left_clip_end;
	clipAlnData_t *left_aln, *left_aln2, *right_aln, *right_aln2;
	bool seg_skip_flag, mid_inner_flag, exist_flag, clip_aln_mate_flag, inner_missing_valid_flag, valid_clip_pos_flag, size_consistent_flag;
	int32_t consistency_flag, vec_id, min_seg_size, end_flag, query_clip_pos, query_clip_pos2, min_dist_vec_id;
	int64_t query_dist, ref_dist;
	double ratio;
	vector<clipPos_t*> clip_pos_items, clip_pos_items_tmp;
	vector<clipAlnData_t*> query_aln_segs;
	vector<int32_t> adjClipAlnSegInfo; // [0]: array index of minimal distance; [1]: segment end flag; [2]: minimal query distance; [3]: minimal reference distance; [4]: array index of minimal distance in the same chromosome; [5]: segment end flag in the same chromosome; [6]: minimal query distance in the same chromosome; [7]: minimal reference distance in the same chromosome
	int32_t idx_vec, mate_arr_idx, mate_clip_end;
	int32_t mate_arr_idx_same_chr, mate_clip_end_same_chr, min_ref_dist_same_chr;

	for(i=0; i<leftClipPosVector.size(); ){
		clip_pos = leftClipPosVector.at(i);

//		if(clip_pos->clip_aln->queryname.compare("SRR8858452.1.86899")==0 or clip_pos->clip_aln->queryname.compare("S1_26622")==0){
//			cout << "line=" << __LINE__ << ", qname=" << clip_pos->clip_aln->queryname << endl;
//		}

		//query_aln_segs = getQueryClipAlnSegs(queryname, clipAlnDataVector);  // get query clip align segments
		query_aln_segs = getQueryClipAlnSegsAll(clip_pos->clip_aln->queryname, clipAlnDataVector);  // get query clip align segments

		if(clip_pos->chrname.compare(chrname)==0){
			if(clip_pos->clip_end==RIGHT_END){ // right end
				if(clip_pos->clip_aln and ((clip_pos->same_orient_flag and clip_pos->clip_aln->ref_dist>=minAlnSize_same_orient) or (clip_pos->same_orient_flag==false and clip_pos->clip_aln->ref_dist>=minAlnSize_diff_orient))){

//					right_aln = getAdjacentClipAlnSegPreferSameChr(clip_pos->clip_aln, clip_pos->clip_end, query_aln_segs, minClipEndSize, MAX_VAR_REG_SIZE);
//					if(right_aln==NULL) right_aln = clip_pos->clip_aln->right_aln;

					seg_skip_flag = false;
					idx_vec = getVecIdxClipAlnData(clip_pos->clip_aln, query_aln_segs);
					adjClipAlnSegInfo = getAdjacentClipAlnSeg(idx_vec, clip_pos->clip_end, query_aln_segs, minClipEndSize, MAX_VAR_REG_SIZE);
					mate_arr_idx = adjClipAlnSegInfo.at(0);
					mate_clip_end = adjClipAlnSegInfo.at(1);
					mate_arr_idx_same_chr = adjClipAlnSegInfo.at(4);
					mate_clip_end_same_chr = adjClipAlnSegInfo.at(5);
					min_ref_dist_same_chr = adjClipAlnSegInfo.at(7);
					if(mate_arr_idx!=-1){ // mated
						// prefer the segments on the same chromosome with near reference distance
						if(mate_arr_idx_same_chr!=-1 and mate_arr_idx!=mate_arr_idx_same_chr and abs(min_ref_dist_same_chr)<=MAX_REF_DIST_SAME_CHR){ // same chromosome with near reference distance
							mate_arr_idx = mate_arr_idx_same_chr;
							mate_clip_end = mate_clip_end_same_chr;
							seg_skip_flag = true;
						}
						right_aln = query_aln_segs.at(mate_arr_idx);
					}else
						right_aln = clip_pos->clip_aln->right_aln;

					// check breakpoint consistency, if true, do nothing
					if(seg_skip_flag or right_aln==NULL) { i++; continue; }
					else{
						consistency_flag = computeBPConsistencyFlag(clip_pos->clip_aln, clip_pos->clip_end, right_aln, mate_clip_end, leftClipPosVector);
						if(consistency_flag==1) { i++; continue; }
					}

					// compute query distance
					inner_missing_valid_flag = true;
					clip_pos_items_tmp = getClipPosItemsByQueryname(clip_pos->clip_aln->queryname, leftClipPosVector);
					if(clip_pos_items_tmp.size()>=2){ // more than two items
						//cout << clip_pos_items_tmp.size() << endl;

						min_dist_clip_pos = getMinDistClipPosItem(clip_pos, clip_pos_items_tmp);
						if(min_dist_clip_pos){
							clip_pos_right_clip_end = clip_pos;
							clip_pos_left_clip_end = min_dist_clip_pos;
							if(clip_pos_right_clip_end->clip_aln->right_aln!=clip_pos_left_clip_end->clip_aln or clip_pos_left_clip_end->clip_aln->left_aln!=clip_pos_right_clip_end->clip_aln){
								inner_missing_valid_flag = false;
							}else if(clip_pos_left_clip_end->clip_aln==right_aln){
								// compute query distance
								if(clip_pos_right_clip_end->clip_end==RIGHT_END) query_clip_pos = clip_pos_right_clip_end->clip_aln->endQueryPos;
								else query_clip_pos = clip_pos_right_clip_end->clip_aln->startQueryPos;
								if(clip_pos_left_clip_end->clip_end==RIGHT_END) query_clip_pos2 = clip_pos_left_clip_end->clip_aln->endQueryPos;
								else query_clip_pos2 = clip_pos_left_clip_end->clip_aln->startQueryPos;
								query_dist = query_clip_pos2 - query_clip_pos;
								if(query_dist<0) query_dist = -query_dist;
								if(query_dist>MAX_INNER_MISSING_IGNORE_SIZE)
									inner_missing_valid_flag = false;
							}
						}
					}

					if(inner_missing_valid_flag and right_aln and clip_pos->clip_aln->rightClipSize>=minClipEndSize and ((right_aln->aln_orient==clip_pos->clip_aln->aln_orient and right_aln->ref_dist>=minAlnSize_same_orient)
							or (right_aln->aln_orient!=clip_pos->clip_aln->aln_orient and right_aln->ref_dist>=minAlnSize_diff_orient))){
						if(clip_pos->clip_aln->aln_orient==right_aln->aln_orient){ // same orient

							right_aln2 = getAdjacentClipAlnSegPreferSameChr(right_aln, RIGHT_END, query_aln_segs, minClipEndSize, MAX_VAR_REG_SIZE);
							if(right_aln2==NULL) right_aln2 = right_aln->right_aln;

							if(right_aln2 and right_aln->rightClipSize>=minClipEndSize and ((right_aln2->aln_orient==right_aln->aln_orient and right_aln2->ref_dist>=minAlnSize_same_orient)
									or (right_aln2->aln_orient!=right_aln->aln_orient and right_aln2->ref_dist>=minAlnSize_diff_orient))){
								if(right_aln2->chrname.compare(clip_pos->chrname)==0){ // add new clip position item

									mid_inner_flag = exist_flag = false;
									if(right_aln->startRefPos>=clip_pos->clip_aln->startRefPos and right_aln->endRefPos<=right_aln2->endRefPos)
										mid_inner_flag = true;
									else if(isSegSelfOverlap(clip_pos->clip_aln, right_aln, maxVarRegSize) or isSegSelfOverlap(right_aln, right_aln2, maxVarRegSize))
										mid_inner_flag = true;
									vec_id = getClipPosVecId(right_aln2, LEFT_END);
									if(vec_id!=-1)
										exist_flag = true;

									size_consistent_flag = true;
									ref_dist = right_aln2->startRefPos - clip_pos->clipRefPos;
									if(ref_dist<0) ref_dist = -ref_dist;
									ratio = (double)ref_dist / right_aln->ref_dist;
									if(ratio<1.0-CLIP_DIFF_LEN_RATIO_SV or ratio>1.0+CLIP_DIFF_LEN_RATIO_SV)
										size_consistent_flag = false;

									if(mid_inner_flag==false and exist_flag==false and size_consistent_flag){
										clip_pos2 = new clipPos_t();
										clip_pos2->chrname = right_aln2->chrname;
										clip_pos2->clip_end = LEFT_END;
										clip_pos2->clip_aln = right_aln2;
										clip_pos2->clipRefPos = right_aln2->startRefPos;
										clip_pos2->clipQueryPos = right_aln2->startQueryPos;
										clip_pos2->aln_orient = right_aln2->aln_orient;
										if(right_aln->aln_orient==right_aln2->aln_orient) clip_pos2->same_orient_flag = true;
										else clip_pos2->same_orient_flag = false;
										clip_pos2->self_overlap_flag = isSegSelfOverlap(right_aln, right_aln2, maxVarRegSize);

	//									if(clip_pos->clipRefPos>=clip_pos2->clipRefPos){
	//										cout << "*********** line=" << __LINE__ << ", i=" << i << ", leftClipPos1=" << clip_pos->clipRefPos << ", leftClipPos2=" << clip_pos2->clipRefPos << endl;
	//									}

										// compute reference distance
										clip_pos_right_clip_end = clip_pos;
										clip_pos_left_clip_end = clip_pos2;
										ref_dist = clip_pos_left_clip_end->clipRefPos - clip_pos_right_clip_end->clipRefPos;
										if(ref_dist<0) ref_dist = -ref_dist;
										if(ref_dist<=MAX_INNER_MISSING_IGNORE_SIZE) // small distance
											leftClipPosVector.push_back(clip_pos2);
										else
											leftClipPosVector2.push_back(clip_pos2);
									}
								}
							}
						}else{ // different orient

							left_aln2 = getAdjacentClipAlnSegPreferSameChr(right_aln, LEFT_END, query_aln_segs, minClipEndSize, MAX_VAR_REG_SIZE);
							if(left_aln2==NULL) left_aln2 = right_aln->left_aln;

							if(left_aln2 and right_aln->leftClipSize>=minClipEndSize and ((left_aln2->aln_orient==right_aln->aln_orient and left_aln2->ref_dist>=minAlnSize_same_orient)
									or (left_aln2->aln_orient!=right_aln->aln_orient and left_aln2->ref_dist>=minAlnSize_diff_orient))){
								if(left_aln2->chrname.compare(clip_pos->chrname)==0){ // add new clip position item

									mid_inner_flag = exist_flag = false;
									if(right_aln->startRefPos>=clip_pos->clip_aln->startRefPos and right_aln->endRefPos<=left_aln2->endRefPos)
										mid_inner_flag = true;
									else if(isSegSelfOverlap(clip_pos->clip_aln, right_aln, maxVarRegSize) or isSegSelfOverlap(right_aln, left_aln2, maxVarRegSize))
										mid_inner_flag = true;
									vec_id = getClipPosVecId(left_aln2, LEFT_END);
									if(vec_id!=-1)
										exist_flag = true;

									size_consistent_flag = true;
									ref_dist = left_aln2->startRefPos - clip_pos->clipRefPos;
									if(ref_dist<0) ref_dist = -ref_dist;
									ratio = (double)ref_dist / right_aln->ref_dist;
									if(ratio<1.0-CLIP_DIFF_LEN_RATIO_SV or ratio>1.0+CLIP_DIFF_LEN_RATIO_SV)
										size_consistent_flag = false;

									if(mid_inner_flag==false and exist_flag==false and size_consistent_flag){
										clip_pos2 = new clipPos_t();
										clip_pos2->chrname = left_aln2->chrname;
										clip_pos2->clip_end = LEFT_END;
										clip_pos2->clip_aln = left_aln2;
										clip_pos2->clipRefPos = left_aln2->startRefPos;
										clip_pos2->clipQueryPos = left_aln2->startQueryPos;
										clip_pos2->aln_orient = left_aln2->aln_orient;
										if(right_aln->aln_orient==left_aln2->aln_orient) clip_pos2->same_orient_flag = true;
										else clip_pos2->same_orient_flag = false;
										clip_pos2->self_overlap_flag = isSegSelfOverlap(right_aln, left_aln2, maxVarRegSize);

	//									if(clip_pos->clipRefPos>=clip_pos2->clipRefPos){
	//										cout << "*********** line=" << __LINE__ << ", i=" << i << ", leftClipPos1=" << clip_pos->clipRefPos << ", leftClipPos2=" << clip_pos2->clipRefPos << endl;
	//									}

										// compute reference distance
										clip_pos_right_clip_end = clip_pos;
										clip_pos_left_clip_end = clip_pos2;
										ref_dist = clip_pos_left_clip_end->clipRefPos - clip_pos_right_clip_end->clipRefPos;
										if(ref_dist<0) ref_dist = -ref_dist;
										if(ref_dist<=MAX_INNER_MISSING_IGNORE_SIZE) // small distance
											leftClipPosVector.push_back(clip_pos2);
										else
											leftClipPosVector2.push_back(clip_pos2);
									}
								}
							}
						}

						min_seg_size = (clip_pos->clip_aln->aln_orient==right_aln->aln_orient) ? minAlnSize_same_orient : minAlnSize_diff_orient;
						// process the left end of right_aln, typically it should be in the third vector
						if(right_aln->leftClipSize>=min_seg_size){
							vec_id = getClipPosVecId(right_aln, LEFT_END);
							if(vec_id==-1){ // not exist, then add it to the third vector
								valid_clip_pos_flag = true;
								if(rightClipPosVector.size()>0 and rightClipPosVector.at(0)->chrname.compare(right_aln->chrname)!=0)
									valid_clip_pos_flag = false;
								else if(rightClipPosVector2.size()>0 and rightClipPosVector2.at(0)->chrname.compare(right_aln->chrname)!=0)
									valid_clip_pos_flag = false;
								if(valid_clip_pos_flag){
									clip_pos2 = new clipPos_t();
									clip_pos2->chrname = right_aln->chrname;
									clip_pos2->clip_end = LEFT_END;
									clip_pos2->clip_aln = right_aln;
									clip_pos2->clipRefPos = right_aln->startRefPos;
									clip_pos2->clipQueryPos = right_aln->startQueryPos;
									clip_pos2->aln_orient = right_aln->aln_orient;
									if(clip_pos->clip_aln->aln_orient==right_aln->aln_orient) clip_pos2->same_orient_flag = true;
									else clip_pos2->same_orient_flag = false;
									clip_pos2->self_overlap_flag = isSegSelfOverlap(clip_pos->clip_aln, right_aln, maxVarRegSize);
									rightClipPosVector.push_back(clip_pos2);
								}
							}
						}

						// process the right end of right_aln, typically should be in the fourth vector
						if(right_aln->rightClipSize>=min_seg_size){
							vec_id = getClipPosVecId(right_aln, RIGHT_END);
							if(vec_id==-1){ // not exist, then add it to the fourth vector
								if(rightClipPosVector2.size()>0 and rightClipPosVector2.at(0)->chrname.compare(right_aln->chrname)!=0)
									valid_clip_pos_flag = false;
								else if(rightClipPosVector.size()>0 and rightClipPosVector.at(0)->chrname.compare(right_aln->chrname)!=0)
									valid_clip_pos_flag = false;
								if(valid_clip_pos_flag){
									clip_pos2 = new clipPos_t();
									clip_pos2->chrname = right_aln->chrname;
									clip_pos2->clip_end = RIGHT_END;
									clip_pos2->clip_aln = right_aln;
									clip_pos2->clipRefPos = right_aln->endRefPos;
									clip_pos2->clipQueryPos = right_aln->endQueryPos;
									clip_pos2->aln_orient = right_aln->aln_orient;
									if(clip_pos->clip_aln->aln_orient==right_aln->aln_orient) clip_pos2->same_orient_flag = true;
									else clip_pos2->same_orient_flag = false;
									clip_pos2->self_overlap_flag = isSegSelfOverlap(clip_pos->clip_aln, right_aln, maxVarRegSize);
									rightClipPosVector2.push_back(clip_pos2);
								}
							}
						}
					}
				}
			}else{ // left end
				if(clip_pos->clip_aln and ((clip_pos->same_orient_flag and clip_pos->clip_aln->ref_dist>=minAlnSize_same_orient) or (clip_pos->same_orient_flag==false and clip_pos->clip_aln->ref_dist>=minAlnSize_diff_orient))){
//					left_aln = getAdjacentClipAlnSegPreferSameChr(clip_pos->clip_aln, clip_pos->clip_end, query_aln_segs, minClipEndSize, MAX_VAR_REG_SIZE);
//					if(left_aln==NULL) left_aln = clip_pos->clip_aln->left_aln;

					seg_skip_flag = false;
					idx_vec = getVecIdxClipAlnData(clip_pos->clip_aln, query_aln_segs);
					adjClipAlnSegInfo = getAdjacentClipAlnSeg(idx_vec, clip_pos->clip_end, query_aln_segs, minClipEndSize, MAX_VAR_REG_SIZE);
					mate_arr_idx = adjClipAlnSegInfo.at(0);
					mate_clip_end = adjClipAlnSegInfo.at(1);
					mate_arr_idx_same_chr = adjClipAlnSegInfo.at(4);
					mate_clip_end_same_chr = adjClipAlnSegInfo.at(5);
					min_ref_dist_same_chr = adjClipAlnSegInfo.at(7);
					if(mate_arr_idx!=-1){ // mated
						// prefer the segments on the same chromosome with near reference distance
						if(mate_arr_idx_same_chr!=-1 and mate_arr_idx!=mate_arr_idx_same_chr and abs(min_ref_dist_same_chr)<=MAX_REF_DIST_SAME_CHR){ // same chromosome with near reference distance
							mate_arr_idx = mate_arr_idx_same_chr;
							mate_clip_end = mate_clip_end_same_chr;
							seg_skip_flag = true;
						}
						left_aln = query_aln_segs.at(mate_arr_idx);
					}else
						left_aln = clip_pos->clip_aln->left_aln;

					// check breakpoint consistency, if consist with other breakpoints with skipped segment, do nothing
					if(seg_skip_flag or left_aln==NULL) { i++; continue; }
					else{
						consistency_flag = computeBPConsistencyFlag(clip_pos->clip_aln, clip_pos->clip_end, left_aln, mate_clip_end, leftClipPosVector);
						if(consistency_flag==1) { i++; continue; }
					}

					// compute query distance
					inner_missing_valid_flag = true;
					clip_pos_items_tmp = getClipPosItemsByQueryname(clip_pos->clip_aln->queryname, leftClipPosVector);
					if(clip_pos_items_tmp.size()>=2){ // more than two items
						//cout << clip_pos_items_tmp.size() << endl;

						min_dist_clip_pos = getMinDistClipPosItem(clip_pos, clip_pos_items_tmp);
						if(min_dist_clip_pos){
							clip_pos_right_clip_end = min_dist_clip_pos;
							clip_pos_left_clip_end = clip_pos;
							if(clip_pos_right_clip_end->clip_aln->right_aln!=clip_pos_left_clip_end->clip_aln or clip_pos_left_clip_end->clip_aln->left_aln!=clip_pos_right_clip_end->clip_aln){
								inner_missing_valid_flag = false;
							}else if(clip_pos_right_clip_end->clip_aln==left_aln){
								// compute query distance
								if(clip_pos_right_clip_end->clip_end==RIGHT_END) query_clip_pos = clip_pos_right_clip_end->clip_aln->endQueryPos;
								else query_clip_pos = clip_pos_right_clip_end->clip_aln->startQueryPos;
								if(clip_pos_left_clip_end->clip_end==RIGHT_END) query_clip_pos2 = clip_pos_left_clip_end->clip_aln->endQueryPos;
								else query_clip_pos2 = clip_pos_left_clip_end->clip_aln->startQueryPos;
								query_dist = query_clip_pos2 - query_clip_pos;
								if(query_dist<0) query_dist = -query_dist;
								if(query_dist>MAX_INNER_MISSING_IGNORE_SIZE)
									inner_missing_valid_flag = false;
							}
						}
					}

					if(inner_missing_valid_flag and left_aln and clip_pos->clip_aln->leftClipSize>=minClipEndSize and ((left_aln->aln_orient==clip_pos->clip_aln->aln_orient and left_aln->ref_dist>=minAlnSize_same_orient)
							or (left_aln->aln_orient!=clip_pos->clip_aln->aln_orient and left_aln->ref_dist>=minAlnSize_diff_orient))){
						if(clip_pos->clip_aln->aln_orient==left_aln->aln_orient){ // same orient
							left_aln2 = getAdjacentClipAlnSegPreferSameChr(left_aln, LEFT_END, query_aln_segs, minClipEndSize, MAX_VAR_REG_SIZE);
							if(left_aln2==NULL) left_aln2 = left_aln->left_aln;

							if(left_aln2 and left_aln->leftClipSize>=minClipEndSize and ((left_aln2->aln_orient==left_aln->aln_orient and left_aln2->ref_dist>=minAlnSize_same_orient)
									or (left_aln2->aln_orient!=left_aln->aln_orient and left_aln2->ref_dist>=minAlnSize_diff_orient))){
								if(left_aln2->chrname.compare(clip_pos->chrname)==0){ // add new clip position item

									mid_inner_flag = exist_flag = false;
									if(left_aln->startRefPos>=clip_pos->clip_aln->startRefPos and left_aln->endRefPos<=left_aln2->endRefPos)
										mid_inner_flag = true;
									else if(isSegSelfOverlap(clip_pos->clip_aln, left_aln, maxVarRegSize) or isSegSelfOverlap(left_aln, left_aln2, maxVarRegSize))
										mid_inner_flag = true;
									vec_id = getClipPosVecId(left_aln2, RIGHT_END);
									if(vec_id!=-1)
										exist_flag = true;

									size_consistent_flag = true;
									ref_dist = left_aln2->endRefPos - clip_pos->clipRefPos;
									if(ref_dist<0) ref_dist = -ref_dist;
									ratio = (double)ref_dist / left_aln->ref_dist;
									if(ratio<1.0-CLIP_DIFF_LEN_RATIO_SV or ratio>1.0+CLIP_DIFF_LEN_RATIO_SV)
										size_consistent_flag = false;

									if(mid_inner_flag==false and exist_flag==false and size_consistent_flag){
										clip_pos2 = new clipPos_t();
										clip_pos2->chrname = left_aln2->chrname;
										clip_pos2->clip_end = RIGHT_END;
										clip_pos2->clip_aln = left_aln2;
										clip_pos2->clipRefPos = left_aln2->endRefPos;
										clip_pos2->clipQueryPos = left_aln2->endQueryPos;
										clip_pos2->aln_orient = left_aln2->aln_orient;
										if(left_aln->aln_orient==left_aln2->aln_orient) clip_pos2->same_orient_flag = true;
										else clip_pos2->same_orient_flag = false;
										clip_pos2->self_overlap_flag = isSegSelfOverlap(left_aln, left_aln2, maxVarRegSize);

	//									if(clip_pos->clipRefPos>=clip_pos2->clipRefPos){
	//										cout << "*********** line=" << __LINE__ << ", i=" << i << ", leftClipPos1=" << clip_pos->clipRefPos << ", leftClipPos2=" << clip_pos2->clipRefPos << endl;
	//									}

										// compute reference distance
										clip_pos_right_clip_end = clip_pos2;
										clip_pos_left_clip_end = clip_pos;
										ref_dist = clip_pos_left_clip_end->clipRefPos - clip_pos_right_clip_end->clipRefPos;
										if(ref_dist<0) ref_dist = -ref_dist;
										if(ref_dist<=MAX_INNER_MISSING_IGNORE_SIZE) // small distance
											leftClipPosVector.push_back(clip_pos2);
										else
											leftClipPosVector2.push_back(clip_pos2);
									}
								}
							}
						}else{ // different orient
							right_aln2 = getAdjacentClipAlnSegPreferSameChr(left_aln, RIGHT_END, query_aln_segs, minClipEndSize, MAX_VAR_REG_SIZE);
							if(right_aln2==NULL) right_aln2 = left_aln->right_aln;

							if(right_aln2 and left_aln->rightClipSize>=minClipEndSize and ((right_aln2->aln_orient==left_aln->aln_orient and right_aln2->ref_dist>=minAlnSize_same_orient)
									or (right_aln2->aln_orient!=left_aln->aln_orient and right_aln2->ref_dist>=minAlnSize_diff_orient))){
								if(right_aln2->chrname.compare(clip_pos->chrname)==0){ // add new clip position item

									mid_inner_flag = exist_flag = false;
									if(left_aln->startRefPos>=clip_pos->clip_aln->startRefPos and left_aln->endRefPos<=right_aln2->endRefPos)
										mid_inner_flag = true;
									else if(isSegSelfOverlap(clip_pos->clip_aln, left_aln, maxVarRegSize) or isSegSelfOverlap(left_aln, right_aln2, maxVarRegSize))
										mid_inner_flag = true;
									vec_id = getClipPosVecId(right_aln2, RIGHT_END);
									if(vec_id!=-1)
										exist_flag = true;

									size_consistent_flag = true;
									ref_dist = right_aln2->endRefPos - clip_pos->clipRefPos;
									if(ref_dist<0) ref_dist = -ref_dist;
									ratio = (double)ref_dist / left_aln->ref_dist;
									if(ratio<1.0-CLIP_DIFF_LEN_RATIO_SV or ratio>1.0+CLIP_DIFF_LEN_RATIO_SV)
										size_consistent_flag = false;

									if(mid_inner_flag==false and exist_flag==false and size_consistent_flag){
										clip_pos2 = new clipPos_t();
										clip_pos2->chrname = right_aln2->chrname;
										clip_pos2->clip_end = RIGHT_END;
										clip_pos2->clip_aln = right_aln2;
										clip_pos2->clipRefPos = right_aln2->endRefPos;
										clip_pos2->clipQueryPos = right_aln2->endQueryPos;
										clip_pos2->aln_orient = right_aln2->aln_orient;
										if(left_aln->aln_orient==right_aln2->aln_orient) clip_pos2->same_orient_flag = true;
										else clip_pos2->same_orient_flag = false;
										clip_pos2->self_overlap_flag = isSegSelfOverlap(left_aln, right_aln2, maxVarRegSize);

	//									if(clip_pos->clipRefPos>=clip_pos2->clipRefPos){
	//										cout << "*********** line=" << __LINE__ << ", i=" << i << ", leftClipPos1=" << clip_pos->clipRefPos << ", leftClipPos2=" << clip_pos2->clipRefPos << endl;
	//									}

										// compute reference distance
										clip_pos_right_clip_end = clip_pos2;
										clip_pos_left_clip_end = clip_pos;
										ref_dist = clip_pos_left_clip_end->clipRefPos - clip_pos_right_clip_end->clipRefPos;
										if(ref_dist<0) ref_dist = -ref_dist;
										if(ref_dist<=MAX_INNER_MISSING_IGNORE_SIZE) // small distance
											leftClipPosVector.push_back(clip_pos2);
										else
											leftClipPosVector2.push_back(clip_pos2);
									}
								}
							}
						}

						// get clipping position items with same align segments by queryname
						clip_aln_mate_flag = true;
						clip_pos_items = getClipPosItemsByQueryname(clip_pos->clip_aln->queryname, leftClipPosVector);
						if(clip_pos_items.size()>=2){ // more than two items
							//cout << clip_pos_items.size() << endl;

							min_dist_clip_pos = getMinDistClipPosItem(clip_pos, clip_pos_items);
							if(min_dist_clip_pos){
								clip_pos_right_clip_end = min_dist_clip_pos;
								clip_pos_left_clip_end = clip_pos;
								if(clip_pos_right_clip_end->clip_aln->right_aln!=clip_pos_left_clip_end->clip_aln or clip_pos_left_clip_end->clip_aln->left_aln!=clip_pos_right_clip_end->clip_aln){
									clip_aln_mate_flag = false;
								}
							}
						}

						// process the right end of clip_aln, typically should be in the second vector
						min_seg_size = (clip_pos->same_orient_flag) ? minAlnSize_same_orient : minAlnSize_diff_orient;
						if(clip_aln_mate_flag and clip_pos->clip_aln->rightClipSize>=min_seg_size){
							min_dist_vec_id = -1;
							valid_clip_pos_flag = true;

							right_aln = getAdjacentClipAlnSegPreferSameChr(clip_pos->clip_aln, RIGHT_END, query_aln_segs, minClipEndSize, MAX_VAR_REG_SIZE);
							if(right_aln==NULL) right_aln = clip_pos->clip_aln->right_aln;

							vec_id = getClipPosVecId(clip_pos->clip_aln, RIGHT_END);
							if(vec_id==-1){ // not exist, then add it to the second vector
								clip_pos2 = new clipPos_t();
								clip_pos2->chrname = clip_pos->clip_aln->chrname;
								clip_pos2->clip_end = RIGHT_END;
								clip_pos2->clip_aln = clip_pos->clip_aln;
								clip_pos2->clipRefPos = clip_pos->clip_aln->endRefPos;
								clip_pos2->clipQueryPos = clip_pos->clip_aln->endQueryPos;
								clip_pos2->aln_orient = clip_pos->clip_aln->aln_orient;
								if(right_aln){
									if(clip_pos->clip_aln->aln_orient==right_aln->aln_orient) clip_pos2->same_orient_flag = true;
									else clip_pos2->same_orient_flag = false;
									clip_pos2->self_overlap_flag = isSegSelfOverlap(clip_pos->clip_aln, right_aln, maxVarRegSize);
								}else {
									clip_pos2->same_orient_flag = clip_pos->same_orient_flag;
									clip_pos2->self_overlap_flag = false;
								}

								// deal with left_aln2==NULL or right_aln2==NULL
								//if((clip_pos->clip_aln->aln_orient==left_aln->aln_orient and left_aln2==NULL) or (clip_pos->clip_aln->aln_orient!=left_aln->aln_orient and right_aln2==NULL)){ // deleted on 2024-01-25
									min_dist_vec_id = getMinDistVecId(clip_pos2);
									if(min_dist_vec_id==-1) valid_clip_pos_flag = false;
								//}
								if(valid_clip_pos_flag){
									if(min_dist_vec_id==-1) leftClipPosVector2.push_back(clip_pos2);
									else{
										if(min_dist_vec_id==0) leftClipPosVector.push_back(clip_pos2);
										else if(min_dist_vec_id==1) leftClipPosVector2.push_back(clip_pos2);
										else if(min_dist_vec_id==2) rightClipPosVector.push_back(clip_pos2);
										else if(min_dist_vec_id==3) rightClipPosVector2.push_back(clip_pos2);
									}
								}else { delete clip_pos2; clip_pos2 = NULL; }
							}

							// add the mate item of clip_pos2
							if(right_aln and (valid_clip_pos_flag and min_dist_vec_id==-1)){
								if(clip_pos->clip_aln->aln_orient==right_aln->aln_orient) end_flag = LEFT_END;
								else end_flag = RIGHT_END;
								vec_id = getClipPosVecId(right_aln, end_flag);
								if(vec_id==-1){ // not exist, then add it to the third or fourth vector
									clip_pos3 = new clipPos_t();
									clip_pos3->chrname = right_aln->chrname;
									clip_pos3->clip_end = end_flag;
									clip_pos3->clip_aln = right_aln;
									clip_pos3->self_overlap_flag = isSegSelfOverlap(clip_pos->clip_aln, right_aln, maxVarRegSize);;
									if(clip_pos->clip_aln->aln_orient==right_aln->aln_orient) {
										clip_pos3->clipRefPos = right_aln->startRefPos;
										clip_pos3->clipQueryPos = right_aln->startQueryPos;
										clip_pos3->aln_orient = right_aln->aln_orient;
										clip_pos3->same_orient_flag = true;
										//rightClipPosVector2.push_back(clip_pos3);
									}else {
										clip_pos3->clipRefPos = right_aln->endRefPos;
										clip_pos3->same_orient_flag = false;
										//rightClipPosVector.push_back(clip_pos3);
									}

									min_dist_vec_id = getMinDistVecId(clip_pos3);
									if(min_dist_vec_id==-1){
										if(clip_pos->clip_aln->aln_orient==right_aln->aln_orient)
											rightClipPosVector2.push_back(clip_pos3);
										else
											rightClipPosVector.push_back(clip_pos3);
									}else{
										if(min_dist_vec_id==0) leftClipPosVector.push_back(clip_pos3);
										else if(min_dist_vec_id==1) leftClipPosVector2.push_back(clip_pos3);
										else if(min_dist_vec_id==2) rightClipPosVector.push_back(clip_pos3);
										else if(min_dist_vec_id==3) rightClipPosVector2.push_back(clip_pos3);
									}
								}
							}
						}
					}
				}
			}
			i++;
		}
//		else if(clip_pos->same_orient_flag==false){
//
//			cout << "=========== line=" << __LINE__ << ", clip_pos=" << clip_pos->clipRefPos << ", same_orient_flag=false, error!" << endl;
//
//			leftClipPosVector2.push_back(clip_pos);
//			leftClipPosVector.erase(leftClipPosVector.begin()+i);
//		}
		else i++;
	}

	for(i=0; i<rightClipPosVector.size(); ){
		clip_pos = rightClipPosVector.at(i);

//		if(clip_pos->clip_aln->queryname.compare("SRR8858452.1.86899")==0 or clip_pos->clip_aln->queryname.compare("S1_343")==0){
//			cout << "line=" << __LINE__ << ", qname=" << clip_pos->clip_aln->queryname << endl;
//		}

		//query_aln_segs = getQueryClipAlnSegs(queryname, clipAlnDataVector);  // get query clip align segments
		query_aln_segs = getQueryClipAlnSegsAll(clip_pos->clip_aln->queryname, clipAlnDataVector);  // get query clip align segments

		//if(clip_pos->chrname.compare(chrname)==0){
			if(clip_pos->clip_end==RIGHT_END){ // right end
				if(clip_pos->clip_aln and ((clip_pos->same_orient_flag and clip_pos->clip_aln->ref_dist>=minAlnSize_same_orient) or (clip_pos->same_orient_flag==false and clip_pos->clip_aln->ref_dist>=minAlnSize_diff_orient))){
					//right_aln = getAdjacentClipAlnSegPreferSameChr(clip_pos->clip_aln, clip_pos->clip_end, query_aln_segs, minClipEndSize, MAX_VAR_REG_SIZE);

					seg_skip_flag = false;
					idx_vec = getVecIdxClipAlnData(clip_pos->clip_aln, query_aln_segs);
					adjClipAlnSegInfo = getAdjacentClipAlnSeg(idx_vec, clip_pos->clip_end, query_aln_segs, minClipEndSize, MAX_VAR_REG_SIZE);
					mate_arr_idx = adjClipAlnSegInfo.at(0);
					mate_clip_end = adjClipAlnSegInfo.at(1);
					mate_arr_idx_same_chr = adjClipAlnSegInfo.at(4);
					mate_clip_end_same_chr = adjClipAlnSegInfo.at(5);
					min_ref_dist_same_chr = adjClipAlnSegInfo.at(7);
					if(mate_arr_idx!=-1){ // mated
						// prefer the segments on the same chromosome with near reference distance
						if(mate_arr_idx_same_chr!=-1 and mate_arr_idx!=mate_arr_idx_same_chr and abs(min_ref_dist_same_chr)<=MAX_REF_DIST_SAME_CHR){ // same chromosome with near reference distance
							mate_arr_idx = mate_arr_idx_same_chr;
							mate_clip_end = mate_clip_end_same_chr;
							seg_skip_flag = true;
						}
						right_aln = query_aln_segs.at(mate_arr_idx);
					}else
						right_aln = clip_pos->clip_aln->right_aln;

					// check breakpoint consistency, if true, do nothing
					if(seg_skip_flag or right_aln==NULL) { i++; continue; }
					else{
						consistency_flag = computeBPConsistencyFlag(clip_pos->clip_aln, clip_pos->clip_end, right_aln, mate_clip_end, rightClipPosVector);
						if(consistency_flag==1) { i++; continue; }
					}

					// compute query distance
					inner_missing_valid_flag = true;
					clip_pos_items_tmp = getClipPosItemsByQueryname(clip_pos->clip_aln->queryname, rightClipPosVector);
					if(clip_pos_items_tmp.size()>=2){ // more than two items
						//cout << clip_pos_items_tmp.size() << endl;

						min_dist_clip_pos = getMinDistClipPosItem(clip_pos, clip_pos_items_tmp);
						if(min_dist_clip_pos){
							clip_pos_right_clip_end = clip_pos;
							clip_pos_left_clip_end = min_dist_clip_pos;
							if(clip_pos_right_clip_end->clip_aln->right_aln!=clip_pos_left_clip_end->clip_aln or clip_pos_left_clip_end->clip_aln->left_aln!=clip_pos_right_clip_end->clip_aln){
								inner_missing_valid_flag = false;
							}else if(clip_pos_left_clip_end->clip_aln==right_aln){
								// compute query distance
								if(clip_pos_right_clip_end->clip_end==RIGHT_END) query_clip_pos = clip_pos_right_clip_end->clip_aln->endQueryPos;
								else query_clip_pos = clip_pos_right_clip_end->clip_aln->startQueryPos;
								if(clip_pos_left_clip_end->clip_end==RIGHT_END) query_clip_pos2 = clip_pos_left_clip_end->clip_aln->endQueryPos;
								else query_clip_pos2 = clip_pos_left_clip_end->clip_aln->startQueryPos;
								query_dist = query_clip_pos2 - query_clip_pos;
								if(query_dist<0) query_dist = -query_dist;
								if(query_dist>MAX_INNER_MISSING_IGNORE_SIZE)
									inner_missing_valid_flag = false;
							}
						}
					}

					if(inner_missing_valid_flag and right_aln and clip_pos->clip_aln->rightClipSize>=minClipEndSize and ((right_aln->aln_orient==clip_pos->clip_aln->aln_orient and right_aln->ref_dist>=minAlnSize_same_orient)
							or (right_aln->aln_orient!=clip_pos->clip_aln->aln_orient and right_aln->ref_dist>=minAlnSize_diff_orient))){
						if(clip_pos->clip_aln->aln_orient==right_aln->aln_orient){ // same orient
							right_aln2 = getAdjacentClipAlnSegPreferSameChr(right_aln, RIGHT_END, query_aln_segs, minClipEndSize, MAX_VAR_REG_SIZE);
							if(right_aln2==NULL) right_aln2 = right_aln->right_aln;

							if(right_aln2 and right_aln->rightClipSize>=minClipEndSize and ((right_aln2->aln_orient==right_aln->aln_orient and right_aln2->ref_dist>=minAlnSize_same_orient)
									or (right_aln2->aln_orient!=right_aln->aln_orient and right_aln2->ref_dist>=minAlnSize_diff_orient))){
								if(right_aln2->chrname.compare(clip_pos->chrname)==0){ // add new clip position item

									mid_inner_flag = exist_flag = false;
									if(right_aln->startRefPos>=clip_pos->clip_aln->startRefPos and right_aln->endRefPos<=right_aln2->endRefPos)
										mid_inner_flag = true;
									else if(isSegSelfOverlap(clip_pos->clip_aln, right_aln, maxVarRegSize) or isSegSelfOverlap(right_aln, right_aln2, maxVarRegSize))
										mid_inner_flag = true;
									vec_id = getClipPosVecId(right_aln2, LEFT_END);
									if(vec_id!=-1)
										exist_flag = true;

									size_consistent_flag = true;
									ref_dist = right_aln2->startRefPos - clip_pos->clipRefPos;
									if(ref_dist<0) ref_dist = -ref_dist;
									ratio = (double)ref_dist / right_aln->ref_dist;
									if(ratio<1.0-CLIP_DIFF_LEN_RATIO_SV or ratio>1.0+CLIP_DIFF_LEN_RATIO_SV)
										size_consistent_flag = false;

									if(mid_inner_flag==false and exist_flag==false and size_consistent_flag){
										clip_pos2 = new clipPos_t();
										clip_pos2->chrname = right_aln2->chrname;
										clip_pos2->clip_end = LEFT_END;
										clip_pos2->clip_aln = right_aln2;
										clip_pos2->clipRefPos = right_aln2->startRefPos;
										clip_pos2->clipQueryPos = right_aln2->startQueryPos;
										clip_pos2->aln_orient = right_aln2->aln_orient;
										if(right_aln->aln_orient==right_aln2->aln_orient) clip_pos2->same_orient_flag = true;
										else clip_pos2->same_orient_flag = false;
										clip_pos2->self_overlap_flag = isSegSelfOverlap(right_aln, right_aln2, maxVarRegSize);

	//									if(clip_pos->clipRefPos>=clip_pos2->clipRefPos){
	//										cout << "----------- line=" << __LINE__ << ", i=" << i << ", leftClipPos1=" << clip_pos->clipRefPos << ", leftClipPos2=" << clip_pos2->clipRefPos << endl;
	//									}

										// compute reference distance
										clip_pos_right_clip_end = clip_pos;
										clip_pos_left_clip_end = clip_pos2;
										ref_dist = clip_pos_left_clip_end->clipRefPos - clip_pos_right_clip_end->clipRefPos;
										if(ref_dist<0) ref_dist = -ref_dist;
										if(ref_dist<=MAX_INNER_MISSING_IGNORE_SIZE) // small distance
											rightClipPosVector.push_back(clip_pos2);
										else
											rightClipPosVector2.push_back(clip_pos2);
									}
								}
							}
						}else{ // different orient
							left_aln2 = getAdjacentClipAlnSegPreferSameChr(right_aln, LEFT_END, query_aln_segs, minClipEndSize, MAX_VAR_REG_SIZE);
							if(left_aln2==NULL) left_aln2 = right_aln->left_aln;

							if(left_aln2 and right_aln->leftClipSize>=minClipEndSize and ((left_aln2->aln_orient==right_aln->aln_orient and left_aln2->ref_dist>=minAlnSize_same_orient)
									or (left_aln2->aln_orient!=right_aln->aln_orient and left_aln2->ref_dist>=minAlnSize_diff_orient))){
								if(left_aln2->chrname.compare(clip_pos->chrname)==0){ // add new clip position item

									mid_inner_flag = exist_flag = false;
									if(right_aln->startRefPos>=clip_pos->clip_aln->startRefPos and right_aln->endRefPos<=left_aln2->endRefPos)
										mid_inner_flag = true;
									else if(isSegSelfOverlap(clip_pos->clip_aln, right_aln, maxVarRegSize) or isSegSelfOverlap(right_aln, left_aln2, maxVarRegSize))
										mid_inner_flag = true;
									vec_id = getClipPosVecId(left_aln2, LEFT_END);
									if(vec_id!=-1)
										exist_flag = true;

									size_consistent_flag = true;
									ref_dist = left_aln2->startRefPos - clip_pos->clipRefPos;
									if(ref_dist<0) ref_dist = -ref_dist;
									ratio = (double)ref_dist / right_aln->ref_dist;
									if(ratio<1.0-CLIP_DIFF_LEN_RATIO_SV or ratio>1.0+CLIP_DIFF_LEN_RATIO_SV)
										size_consistent_flag = false;

									if(mid_inner_flag==false and exist_flag==false and size_consistent_flag){
										clip_pos2 = new clipPos_t();
										clip_pos2->chrname = left_aln2->chrname;
										clip_pos2->clip_end = LEFT_END;
										clip_pos2->clip_aln = left_aln2;
										clip_pos2->clipRefPos = left_aln2->startRefPos;
										clip_pos2->clipQueryPos = left_aln2->startQueryPos;
										clip_pos2->aln_orient = left_aln2->aln_orient;
										if(right_aln->aln_orient==left_aln2->aln_orient) clip_pos2->same_orient_flag = true;
										else clip_pos2->same_orient_flag = false;
										clip_pos2->self_overlap_flag = isSegSelfOverlap(right_aln, left_aln2, maxVarRegSize);

	//									if(clip_pos->clipRefPos>=clip_pos2->clipRefPos){
	//										cout << "----------- line=" << __LINE__ << ", i=" << i << ", leftClipPos1=" << clip_pos->clipRefPos << ", leftClipPos2=" << clip_pos2->clipRefPos << endl;
	//									}

										// compute reference distance
										clip_pos_right_clip_end = clip_pos;
										clip_pos_left_clip_end = clip_pos2;
										ref_dist = clip_pos_left_clip_end->clipRefPos - clip_pos_right_clip_end->clipRefPos;
										if(ref_dist<0) ref_dist = -ref_dist;
										if(ref_dist<=MAX_INNER_MISSING_IGNORE_SIZE) // small distance
											rightClipPosVector.push_back(clip_pos2);
										else
											rightClipPosVector2.push_back(clip_pos2);
									}
								}
							}
						}

						min_seg_size = (clip_pos->clip_aln->aln_orient==right_aln->aln_orient) ? minAlnSize_same_orient : minAlnSize_diff_orient;
						// process the left end of right_aln, typically it should be in the first vector
						if(right_aln->leftClipSize>=min_seg_size){
							vec_id = getClipPosVecId(right_aln, LEFT_END);
							if(vec_id==-1){ // not exist, then add it to the first vector
								valid_clip_pos_flag = true;
								if(leftClipPosVector.size()>0 and leftClipPosVector.at(0)->chrname.compare(right_aln->chrname)!=0)
									valid_clip_pos_flag = false;
								else if(leftClipPosVector2.size()>0 and leftClipPosVector2.at(0)->chrname.compare(right_aln->chrname)!=0)
									valid_clip_pos_flag = false;
								if(valid_clip_pos_flag){
									clip_pos2 = new clipPos_t();
									clip_pos2->chrname = right_aln->chrname;
									clip_pos2->clip_end = LEFT_END;
									clip_pos2->clip_aln = right_aln;
									clip_pos2->clipRefPos = right_aln->startRefPos;
									clip_pos2->clipQueryPos = right_aln->startQueryPos;
									clip_pos2->aln_orient = right_aln->aln_orient;
									if(clip_pos->clip_aln->aln_orient==right_aln->aln_orient) clip_pos2->same_orient_flag = true;
									else clip_pos2->same_orient_flag = false;
									clip_pos2->self_overlap_flag = isSegSelfOverlap(clip_pos->clip_aln, right_aln, maxVarRegSize);
									leftClipPosVector.push_back(clip_pos2);
								}
							}
						}

						// process the right end of right_aln, typically should be in the second vector
						if(right_aln->rightClipSize>=min_seg_size){
							vec_id = getClipPosVecId(right_aln, RIGHT_END);
							if(vec_id==-1){ // not exist, then add it to the second vector
								valid_clip_pos_flag = true;
								if(leftClipPosVector2.size()>0 and leftClipPosVector2.at(0)->chrname.compare(right_aln->chrname)!=0)
									valid_clip_pos_flag = false;
								else if(leftClipPosVector.size()>0 and leftClipPosVector.at(0)->chrname.compare(right_aln->chrname)!=0)
									valid_clip_pos_flag = false;
								if(valid_clip_pos_flag){
									clip_pos2 = new clipPos_t();
									clip_pos2->chrname = right_aln->chrname;
									clip_pos2->clip_end = RIGHT_END;
									clip_pos2->clip_aln = right_aln;
									clip_pos2->clipRefPos = right_aln->endRefPos;
									clip_pos2->clipQueryPos = right_aln->endQueryPos;
									clip_pos2->aln_orient = right_aln->aln_orient;
									if(clip_pos->clip_aln->aln_orient==right_aln->aln_orient) clip_pos2->same_orient_flag = true;
									else clip_pos2->same_orient_flag = false;
									clip_pos2->self_overlap_flag = isSegSelfOverlap(clip_pos->clip_aln, right_aln, maxVarRegSize);
									leftClipPosVector2.push_back(clip_pos2);
								}
							}
						}
					}
				}
			}else{ // left end
				if(clip_pos->clip_aln and ((clip_pos->same_orient_flag and clip_pos->clip_aln->ref_dist>=minAlnSize_same_orient) or (clip_pos->same_orient_flag==false and clip_pos->clip_aln->ref_dist>=minAlnSize_diff_orient))){
//					left_aln = getAdjacentClipAlnSegPreferSameChr(clip_pos->clip_aln, clip_pos->clip_end, query_aln_segs, minClipEndSize, MAX_VAR_REG_SIZE);
//					if(left_aln==NULL) left_aln = clip_pos->clip_aln->left_aln;

					seg_skip_flag = false;
					idx_vec = getVecIdxClipAlnData(clip_pos->clip_aln, query_aln_segs);
					adjClipAlnSegInfo = getAdjacentClipAlnSeg(idx_vec, clip_pos->clip_end, query_aln_segs, minClipEndSize, MAX_VAR_REG_SIZE);
					mate_arr_idx = adjClipAlnSegInfo.at(0);
					mate_clip_end = adjClipAlnSegInfo.at(1);
					mate_arr_idx_same_chr = adjClipAlnSegInfo.at(4);
					mate_clip_end_same_chr = adjClipAlnSegInfo.at(5);
					min_ref_dist_same_chr = adjClipAlnSegInfo.at(7);
					if(mate_arr_idx!=-1){ // mated
						// prefer the segments on the same chromosome with near reference distance
						if(mate_arr_idx_same_chr!=-1 and mate_arr_idx!=mate_arr_idx_same_chr and abs(min_ref_dist_same_chr)<=MAX_REF_DIST_SAME_CHR){ // same chromosome with near reference distance
							mate_arr_idx = mate_arr_idx_same_chr;
							mate_clip_end = mate_clip_end_same_chr;
							seg_skip_flag = true;
						}
						left_aln = query_aln_segs.at(mate_arr_idx);
					}else
						left_aln = clip_pos->clip_aln->left_aln;

					// check breakpoint consistency, if true, do nothing
					if(seg_skip_flag or left_aln==NULL) { i++; continue; }
					else{
						consistency_flag = computeBPConsistencyFlag(clip_pos->clip_aln, clip_pos->clip_end, left_aln, mate_clip_end, rightClipPosVector);
						if(consistency_flag==1) { i++; continue; }
					}

					// compute query distance
					inner_missing_valid_flag = true;
					clip_pos_items_tmp = getClipPosItemsByQueryname(clip_pos->clip_aln->queryname, rightClipPosVector);
					if(clip_pos_items_tmp.size()>=2){ // more than two items
						//cout << clip_pos_items_tmp.size() << endl;

						min_dist_clip_pos = getMinDistClipPosItem(clip_pos, clip_pos_items_tmp);
						if(min_dist_clip_pos){
							clip_pos_right_clip_end = min_dist_clip_pos;
							clip_pos_left_clip_end = clip_pos;
							if(clip_pos_right_clip_end->clip_aln->right_aln!=clip_pos_left_clip_end->clip_aln or clip_pos_left_clip_end->clip_aln->left_aln!=clip_pos_right_clip_end->clip_aln){
								inner_missing_valid_flag = false;
							}else if(clip_pos_right_clip_end->clip_aln==left_aln){
								// compute query distance
								if(clip_pos_right_clip_end->clip_end==RIGHT_END) query_clip_pos = clip_pos_right_clip_end->clip_aln->endQueryPos;
								else query_clip_pos = clip_pos_right_clip_end->clip_aln->startQueryPos;
								if(clip_pos_left_clip_end->clip_end==RIGHT_END) query_clip_pos2 = clip_pos_left_clip_end->clip_aln->endQueryPos;
								else query_clip_pos2 = clip_pos_left_clip_end->clip_aln->startQueryPos;
								query_dist = query_clip_pos2 - query_clip_pos;
								if(query_dist<0) query_dist = -query_dist;
								if(query_dist>MAX_INNER_MISSING_IGNORE_SIZE)
									inner_missing_valid_flag = false;
							}
						}
					}

					if(inner_missing_valid_flag and left_aln and clip_pos->clip_aln->leftClipSize>=minClipEndSize and ((left_aln->aln_orient==clip_pos->clip_aln->aln_orient and left_aln->ref_dist>=minAlnSize_same_orient)
							or (left_aln->aln_orient!=clip_pos->clip_aln->aln_orient and left_aln->ref_dist>=minAlnSize_diff_orient))){
						if(clip_pos->clip_aln->aln_orient==left_aln->aln_orient){ // same orient
							left_aln2 = getAdjacentClipAlnSegPreferSameChr(left_aln, LEFT_END, query_aln_segs, minClipEndSize, MAX_VAR_REG_SIZE);
							if(left_aln2==NULL) left_aln2 = left_aln->left_aln;

							if(left_aln2 and left_aln->leftClipSize>=minClipEndSize and ((left_aln2->aln_orient==left_aln->aln_orient and left_aln2->ref_dist>=minAlnSize_same_orient)
									or (left_aln2->aln_orient!=left_aln->aln_orient and left_aln2->ref_dist>=minAlnSize_diff_orient))){
								if(left_aln2->chrname.compare(clip_pos->chrname)==0){ // add new clip position item

									mid_inner_flag = exist_flag = false;
									if(left_aln->startRefPos>=clip_pos->clip_aln->startRefPos and left_aln->endRefPos<=left_aln2->endRefPos)
										mid_inner_flag = true;
									else if(isSegSelfOverlap(clip_pos->clip_aln, left_aln, maxVarRegSize) or isSegSelfOverlap(left_aln, left_aln2, maxVarRegSize))
										mid_inner_flag = true;
									vec_id = getClipPosVecId(left_aln2, RIGHT_END);
									if(vec_id!=-1)
										exist_flag = true;

									size_consistent_flag = true;
									ref_dist = left_aln2->endRefPos - clip_pos->clipRefPos;
									if(ref_dist<0) ref_dist = -ref_dist;
									ratio = (double)ref_dist / left_aln->ref_dist;
									if(ratio<1.0-CLIP_DIFF_LEN_RATIO_SV or ratio>1.0+CLIP_DIFF_LEN_RATIO_SV)
										size_consistent_flag = false;

									if(mid_inner_flag==false and exist_flag==false and size_consistent_flag){
										clip_pos2 = new clipPos_t();
										clip_pos2->chrname = left_aln2->chrname;
										clip_pos2->clip_end = RIGHT_END;
										clip_pos2->clip_aln = left_aln2;
										clip_pos2->clipRefPos = left_aln2->endRefPos;
										clip_pos2->clipQueryPos = left_aln2->endQueryPos;
										clip_pos2->aln_orient = left_aln2->aln_orient;
										if(left_aln->aln_orient==left_aln2->aln_orient) clip_pos2->same_orient_flag = true;
										else clip_pos2->same_orient_flag = false;
										clip_pos2->self_overlap_flag = isSegSelfOverlap(left_aln, left_aln2, maxVarRegSize);

	//									if(clip_pos->clipRefPos>=clip_pos2->clipRefPos){
	//										cout << "----------- line=" << __LINE__ << ", i=" << i << ", leftClipPos1=" << clip_pos->clipRefPos << ", leftClipPos2=" << clip_pos2->clipRefPos << endl;
	//									}

										// compute reference distance
										clip_pos_right_clip_end = clip_pos2;
										clip_pos_left_clip_end = clip_pos;
										ref_dist = clip_pos_left_clip_end->clipRefPos - clip_pos_right_clip_end->clipRefPos;
										if(ref_dist<0) ref_dist = -ref_dist;
										if(ref_dist<=MAX_INNER_MISSING_IGNORE_SIZE) // small distance
											rightClipPosVector.push_back(clip_pos2);
										else
											rightClipPosVector2.push_back(clip_pos2);
									}
								}
							}
						}else{ // different orient
							right_aln2 = getAdjacentClipAlnSegPreferSameChr(left_aln, RIGHT_END, query_aln_segs, minClipEndSize, MAX_VAR_REG_SIZE);
							if(right_aln2==NULL) right_aln2 = left_aln->right_aln;

							if(right_aln2 and left_aln->rightClipSize>=minClipEndSize and ((right_aln2->aln_orient==left_aln->aln_orient and right_aln2->ref_dist>=minAlnSize_same_orient)
									or (right_aln2->aln_orient!=left_aln->aln_orient and right_aln2->ref_dist>=minAlnSize_diff_orient))){
								if(right_aln2->chrname.compare(clip_pos->chrname)==0){ // add new clip position item

									mid_inner_flag = exist_flag = false;
									if(left_aln->startRefPos>=clip_pos->clip_aln->startRefPos and left_aln->endRefPos<=right_aln2->endRefPos)
										mid_inner_flag = true;
									else if(isSegSelfOverlap(clip_pos->clip_aln, left_aln, maxVarRegSize) or isSegSelfOverlap(left_aln, right_aln2, maxVarRegSize))
										mid_inner_flag = true;
									vec_id = getClipPosVecId(right_aln2, RIGHT_END);
									if(vec_id!=-1)
										exist_flag = true;

									size_consistent_flag = true;
									ref_dist = right_aln2->endRefPos - clip_pos->clipRefPos;
									if(ref_dist<0) ref_dist = -ref_dist;
									ratio = (double)ref_dist / left_aln->ref_dist;
									if(ratio<1.0-CLIP_DIFF_LEN_RATIO_SV or ratio>1.0+CLIP_DIFF_LEN_RATIO_SV)
										size_consistent_flag = false;

									if(mid_inner_flag==false and exist_flag==false and size_consistent_flag){
										clip_pos2 = new clipPos_t();
										clip_pos2->chrname = right_aln2->chrname;
										clip_pos2->clip_end = RIGHT_END;
										clip_pos2->clip_aln = right_aln2;
										clip_pos2->clipRefPos = right_aln2->endRefPos;
										clip_pos2->clipQueryPos = right_aln2->endQueryPos;
										clip_pos2->aln_orient = right_aln2->aln_orient;
										if(left_aln->aln_orient==right_aln2->aln_orient) clip_pos2->same_orient_flag = true;
										else clip_pos2->same_orient_flag = false;
										clip_pos2->self_overlap_flag = isSegSelfOverlap(left_aln, right_aln2, maxVarRegSize);

	//									if(clip_pos->clipRefPos>=clip_pos2->clipRefPos){
	//										cout << "----------- line=" << __LINE__ << ", i=" << i << ", leftClipPos1=" << clip_pos->clipRefPos << ", leftClipPos2=" << clip_pos2->clipRefPos << endl;
	//									}

										// compute reference distance
										clip_pos_right_clip_end = clip_pos2;
										clip_pos_left_clip_end = clip_pos;
										ref_dist = clip_pos_left_clip_end->clipRefPos - clip_pos_right_clip_end->clipRefPos;
										if(ref_dist<0) ref_dist = -ref_dist;
										if(ref_dist<=MAX_INNER_MISSING_IGNORE_SIZE) // small distance
											rightClipPosVector.push_back(clip_pos2);
										else
											rightClipPosVector2.push_back(clip_pos2);
									}
								}
							}
						}

						// get clipping position items with same align segments by queryname
						clip_aln_mate_flag = true;
						clip_pos_items = getClipPosItemsByQueryname(clip_pos->clip_aln->queryname, rightClipPosVector);
						if(clip_pos_items.size()>=2){ // more than two items
							//cout << clip_pos_items.size() << endl;

							min_dist_clip_pos = getMinDistClipPosItem(clip_pos, clip_pos_items);
							if(min_dist_clip_pos){
								clip_pos_right_clip_end = min_dist_clip_pos;
								clip_pos_left_clip_end = clip_pos;
								if(clip_pos_right_clip_end->clip_aln->right_aln!=clip_pos_left_clip_end->clip_aln or clip_pos_left_clip_end->clip_aln->left_aln!=clip_pos_right_clip_end->clip_aln){
									clip_aln_mate_flag = false;
								}
							}
						}

						min_seg_size = (clip_pos->same_orient_flag) ? minAlnSize_same_orient : minAlnSize_diff_orient;
						// process the left end of right_aln, typically it should be in the fourth vector
						if(clip_aln_mate_flag and clip_pos->clip_aln->rightClipSize>=min_seg_size){
							min_dist_vec_id = -1;
							valid_clip_pos_flag = true;

							right_aln = getAdjacentClipAlnSegPreferSameChr(clip_pos->clip_aln, RIGHT_END, query_aln_segs, minClipEndSize, MAX_VAR_REG_SIZE);
							if(right_aln==NULL) right_aln = clip_pos->clip_aln->right_aln;

							vec_id = getClipPosVecId(clip_pos->clip_aln, RIGHT_END);
							if(vec_id==-1){ // not exist, then add it to the second vector
								clip_pos2 = new clipPos_t();
								clip_pos2->chrname = clip_pos->clip_aln->chrname;
								clip_pos2->clip_end = RIGHT_END;
								clip_pos2->clip_aln = clip_pos->clip_aln;
								clip_pos2->clipRefPos = clip_pos->clip_aln->endRefPos;
								clip_pos2->clipQueryPos = clip_pos->clip_aln->endQueryPos;
								clip_pos2->aln_orient = clip_pos->clip_aln->aln_orient;
								if(right_aln){
									if(clip_pos->clip_aln->aln_orient==right_aln->aln_orient) clip_pos2->same_orient_flag = true;
									else clip_pos2->same_orient_flag = false;
									clip_pos2->self_overlap_flag = isSegSelfOverlap(clip_pos->clip_aln, right_aln, maxVarRegSize);
								}else {
									clip_pos2->same_orient_flag = clip_pos->same_orient_flag;
									clip_pos2->self_overlap_flag = false;
								}

								// deal with left_aln2==NULL or right_aln2==NULL
								//if((clip_pos->clip_aln->aln_orient==left_aln->aln_orient and left_aln2==NULL) or (clip_pos->clip_aln->aln_orient!=left_aln->aln_orient and right_aln2==NULL)){ // deleted on 2024-01-25
									min_dist_vec_id = getMinDistVecId(clip_pos2);
									if(min_dist_vec_id==-1) valid_clip_pos_flag = false;
								//}
								if(valid_clip_pos_flag){
									if(min_dist_vec_id==-1) rightClipPosVector2.push_back(clip_pos2);
									else{
										if(min_dist_vec_id==0) leftClipPosVector.push_back(clip_pos2);
										else if(min_dist_vec_id==1) leftClipPosVector2.push_back(clip_pos2);
										else if(min_dist_vec_id==2) rightClipPosVector.push_back(clip_pos2);
										else if(min_dist_vec_id==3) rightClipPosVector2.push_back(clip_pos2);
									}
								}else { delete clip_pos2; clip_pos2 = NULL; }

								//rightClipPosVector2.push_back(clip_pos2);
							}

							// add the mate item of clip_pos2
							if(right_aln and (valid_clip_pos_flag and min_dist_vec_id==-1)){
								if(clip_pos->clip_aln->aln_orient==right_aln->aln_orient) end_flag = LEFT_END;
								else end_flag = RIGHT_END;
								vec_id = getClipPosVecId(right_aln, end_flag);
								if(vec_id==-1){ // not exist, then add it to the first or second vector
									clip_pos3 = new clipPos_t();
									clip_pos3->chrname = right_aln->chrname;
									clip_pos3->clip_end = end_flag;
									clip_pos3->clip_aln = right_aln;
									clip_pos3->self_overlap_flag = isSegSelfOverlap(clip_pos->clip_aln, right_aln, maxVarRegSize);
									if(clip_pos->clip_aln->aln_orient==right_aln->aln_orient){
										clip_pos3->clipRefPos = right_aln->startRefPos;
										clip_pos3->clipQueryPos = right_aln->startQueryPos;
										clip_pos3->aln_orient = right_aln->aln_orient;
										clip_pos3->same_orient_flag = true;
										//leftClipPosVector2.push_back(clip_pos3);
									}else{
										clip_pos3->clipRefPos = right_aln->endRefPos;
										clip_pos3->clipQueryPos = right_aln->endQueryPos;
										clip_pos3->aln_orient = right_aln->aln_orient;
										clip_pos3->same_orient_flag = false;
										//leftClipPosVector.push_back(clip_pos3);
									}

									min_dist_vec_id = getMinDistVecId(clip_pos3);
									if(min_dist_vec_id==-1){
										if(clip_pos->clip_aln->aln_orient==right_aln->aln_orient)
											leftClipPosVector2.push_back(clip_pos3);
										else
											leftClipPosVector.push_back(clip_pos3);
									}else{
										if(min_dist_vec_id==0) leftClipPosVector.push_back(clip_pos3);
										else if(min_dist_vec_id==1) leftClipPosVector2.push_back(clip_pos3);
										else if(min_dist_vec_id==2) rightClipPosVector.push_back(clip_pos3);
										else if(min_dist_vec_id==3) rightClipPosVector2.push_back(clip_pos3);
									}
								}
							}
						}
					}
				}
			}
			i++;
//		}else if(clip_pos->same_orient_flag==false){
//			rightClipPosVector2.push_back(clip_pos);
//			rightClipPosVector.erase(rightClipPosVector.begin()+i);
//		}else i++;
	}
}

// compute the consistency flag for two clip align segments
int32_t clipReg::computeBPConsistencyFlag(clipAlnData_t *clip_aln, int32_t clip_end, clipAlnData_t *adj_aln, int32_t adj_clip_end, vector<clipPos_t*> &clip_pos_vec){
	int32_t consist_flag;
	int64_t clip_loc, adj_clip_loc, mate_clip_loc;
	size_t i;
	clipPos_t *clip_pos_tmp;
	clipAlnData_t *mate_clip_aln_tmp;
	vector<clipAlnData_t*> query_aln_segs;
	vector<int32_t> adjClipAlnSegInfo; // [0]: array index of minimal distance; [1]: segment end flag; [2]: minimal query distance; [3]: minimal reference distance; [4]: array index of minimal distance in the same chromosome; [5]: segment end flag in the same chromosome; [6]: minimal query distance in the same chromosome; [7]: minimal reference distance in the same chromosome
	int32_t idx_vec, mate_arr_idx, mate_clip_end, consist_num, reverse_consist_num;
	int32_t mate_arr_idx_same_chr, min_ref_dist_same_chr;

	if(clip_aln){
		if(clip_end==LEFT_END){
			clip_loc = clip_aln->startRefPos;
		}else{
			clip_loc = clip_aln->endRefPos;
		}
	}
	if(adj_aln){
		if(adj_clip_end==LEFT_END){
			adj_clip_loc = adj_aln->startRefPos;
		}else{
			adj_clip_loc = adj_aln->endRefPos;
		}
	}else{
		adj_clip_loc = -1;
	}

	consist_num = reverse_consist_num = 0;
	if(adj_aln){ // multiple segments
		for(i=0; i<clip_pos_vec.size(); i++){
			clip_pos_tmp = clip_pos_vec.at(i);
			if(clip_pos_tmp->clip_aln!=clip_aln and clip_pos_tmp->chrname.compare(clip_aln->chrname)==0 and ((clip_pos_tmp->clip_end==clip_end and abs(clip_pos_tmp->clipRefPos-clip_loc)<=GROUP_DIST_THRES) or (clip_pos_tmp->clip_end==adj_clip_end and abs(clip_pos_tmp->clipRefPos-adj_clip_loc)<=GROUP_DIST_THRES))){

				query_aln_segs = getQueryClipAlnSegsAll(clip_pos_tmp->clip_aln->queryname, clipAlnDataVector);  // get query clip align segments

				idx_vec = getVecIdxClipAlnData(clip_pos_tmp->clip_aln, query_aln_segs);
				adjClipAlnSegInfo = getAdjacentClipAlnSeg(idx_vec, clip_pos_tmp->clip_end, query_aln_segs, minClipEndSize, MAX_VAR_REG_SIZE);
				mate_arr_idx = adjClipAlnSegInfo.at(0);
				mate_clip_end = adjClipAlnSegInfo.at(1);
				mate_arr_idx_same_chr = adjClipAlnSegInfo.at(4);
	//			mate_clip_end_same_chr = adjClipAlnSegInfo.at(5);
				min_ref_dist_same_chr = adjClipAlnSegInfo.at(7);
				if(mate_arr_idx!=-1){ // mated
					mate_clip_aln_tmp = query_aln_segs.at(mate_arr_idx);
					if(mate_clip_end==LEFT_END){
						mate_clip_loc = mate_clip_aln_tmp->startRefPos;
					}else{
						mate_clip_loc = mate_clip_aln_tmp->endRefPos;
					}

					// only consider the clipPos item having skip segment
					if(mate_arr_idx_same_chr!=-1 and mate_arr_idx!=mate_arr_idx_same_chr){
						// prefer the segments on the same chromosome with near reference distance
						if(abs(min_ref_dist_same_chr)<=MAX_REF_DIST_SAME_CHR){ // same chromosome with near reference distance
		//					mate_arr_idx = mate_arr_idx_same_chr;
		//					mate_clip_end = mate_clip_end_same_chr;

							if(mate_clip_aln_tmp->chrname.compare(adj_aln->chrname)==0 and clip_pos_tmp->clip_end==clip_end and abs(clip_pos_tmp->clipRefPos-clip_loc)<=GROUP_DIST_THRES and mate_clip_end==adj_clip_end and abs(mate_clip_loc-adj_clip_loc)<=GROUP_DIST_THRES){ // consist
								consist_num ++;
							}
						}else { // same chromosome with large reference distance
							if(mate_clip_aln_tmp->chrname.compare(clip_aln->chrname)==0 and clip_pos_tmp->clip_end==adj_clip_end and abs(clip_pos_tmp->clipRefPos-adj_clip_loc)<=GROUP_DIST_THRES and mate_clip_end==clip_end and abs(mate_clip_loc-clip_loc)<=GROUP_DIST_THRES){
								reverse_consist_num ++;
							}
						}
					}
				}
			}
		}
	}else{ // only one segment
		for(i=0; i<clip_pos_vec.size(); i++){
			clip_pos_tmp = clip_pos_vec.at(i);
			if(clip_pos_tmp->clip_aln!=clip_aln and clip_pos_tmp->chrname.compare(clip_aln->chrname)==0 and clip_pos_tmp->clip_end==clip_end and abs(clip_pos_tmp->clipRefPos-clip_loc)<=GROUP_DIST_THRES){
				consist_num ++;
			}
		}
	}

	//cout << "consist_num=" << consist_num << ", reverse_consist_num=" << reverse_consist_num << endl;

	consist_flag = 0;
	//if(consist_num>=paras->minReadsNumSupportSV*READS_NUM_SUPPORT_FACTOR)
	if(consist_num>=paras->minReadsNumSupportSV)
		consist_flag = 1;
	//else if(reverse_consist_num>=paras->minReadsNumSupportSV*READS_NUM_SUPPORT_FACTOR)
	else if(reverse_consist_num>=paras->minReadsNumSupportSV)
		consist_flag = 2;

	return consist_flag;
}

// compute the consistency flag for multiple vectors
void clipReg::AdjustClipPosVecByBPConsistency(vector<clipPos_t*> &clip_pos_vec, int32_t vec_id){
	if(clip_pos_vec.size()>0){
		switch(vec_id){
			case 0:
				AdjustClipPosVecByBPConsistencySingleVec(clip_pos_vec, leftClipPosVector2);
				AdjustClipPosVecByBPConsistencySingleVec(clip_pos_vec, rightClipPosVector);
				AdjustClipPosVecByBPConsistencySingleVec(clip_pos_vec, rightClipPosVector2);
				break;
			case 1:
				AdjustClipPosVecByBPConsistencySingleVec(clip_pos_vec, leftClipPosVector);
				AdjustClipPosVecByBPConsistencySingleVec(clip_pos_vec, rightClipPosVector);
				AdjustClipPosVecByBPConsistencySingleVec(clip_pos_vec, rightClipPosVector2);
				break;
			case 2:
				AdjustClipPosVecByBPConsistencySingleVec(clip_pos_vec, rightClipPosVector2);
				AdjustClipPosVecByBPConsistencySingleVec(clip_pos_vec, leftClipPosVector);
				AdjustClipPosVecByBPConsistencySingleVec(clip_pos_vec, leftClipPosVector2);
				break;
			case 3:
				AdjustClipPosVecByBPConsistencySingleVec(clip_pos_vec, rightClipPosVector);
				AdjustClipPosVecByBPConsistencySingleVec(clip_pos_vec, leftClipPosVector);
				AdjustClipPosVecByBPConsistencySingleVec(clip_pos_vec, leftClipPosVector2);
				break;
			default: cerr << "invalid vec_id=" << vec_id << endl; exit(1);
		}
	}
}

// compute the consistency flag from vector
void clipReg::AdjustClipPosVecByBPConsistencySingleVec(vector<clipPos_t*> &clip_pos_vec, vector<clipPos_t*> &clip_pos_vec_main){
	//bool self_overlap_flag, skip_flag;
	clipPos_t *clip_pos_item;
	vector<clipAlnData_t*> query_aln_segs;
	//int32_t consist_flag, idx_vec, mate_arr_idx, mate_clip_end, mate_arr_idx_same_chr, mate_clip_end_same_chr, min_dist_same_chr, min_ref_dist_same_chr;
	int32_t consist_flag, idx_vec, mate_arr_idx, mate_clip_end;
	vector<int32_t> adjClipAlnSegInfo; // [0]: array index of minimal distance; [1]: segment end flag; [2]: minimal query distance; [3]: minimal reference distance; [4]: array index of minimal distance in the same chromosome; [5]: segment end flag in the same chromosome; [6]: minimal query distance in the same chromosome; [7]: minimal reference distance in the same chromosome
	clipAlnData_t *mate_clip_aln_item;

	if(clip_pos_vec_main.size()>0){
		for(size_t i=0; i<clip_pos_vec.size(); ){
			clip_pos_item = clip_pos_vec.at(i);

			query_aln_segs = getQueryClipAlnSegsAll(clip_pos_item->clip_aln->queryname, clipAlnDataVector);  // get query clip align segments

			//skip_flag = false;
			idx_vec = getVecIdxClipAlnData(clip_pos_item->clip_aln, query_aln_segs);
			adjClipAlnSegInfo = getAdjacentClipAlnSeg(idx_vec, clip_pos_item->clip_end, query_aln_segs, minClipEndSize, MAX_VAR_REG_SIZE);
			mate_arr_idx = adjClipAlnSegInfo.at(0);
			mate_clip_end = adjClipAlnSegInfo.at(1);
			//mate_arr_idx_same_chr = adjClipAlnSegInfo.at(4);
			//mate_clip_end_same_chr = adjClipAlnSegInfo.at(5);
			//min_dist_same_chr = adjClipAlnSegInfo.at(6);
			//min_ref_dist_same_chr = adjClipAlnSegInfo.at(7);
			if(mate_arr_idx!=-1){ // mated
//				if(mate_arr_idx_same_chr!=-1){
//					mate_clip_aln_item = query_aln_segs.at(mate_arr_idx_same_chr);
//					self_overlap_flag = isSegSelfOverlap(clip_pos_item->clip_aln, mate_clip_aln_item, maxVarRegSize);
//					// prefer the segments on the same chromosome
//					if(mate_arr_idx!=mate_arr_idx_same_chr and (abs(min_ref_dist_same_chr)<=MAX_REF_DIST_SAME_CHR or self_overlap_flag)){ // different chromosomes with near reference distance
//					//if(abs(min_ref_dist_same_chr)<=MAX_REF_DIST_SAME_CHR){ // different chromosomes with near reference distance
//						//mate_arr_idx = mate_arr_idx_same_chr;
//						//mate_clip_end_flag = mate_clip_end_flag_same_chr;
//						skip_flag = true;
//					}else if(mate_arr_idx==mate_arr_idx_same_chr and mate_clip_end==mate_clip_end_same_chr and abs(min_dist_same_chr)<=MAX_REF_DIST_SAME_CHR){
//						skip_flag = true;
//					}
//				}
				mate_clip_aln_item = query_aln_segs.at(mate_arr_idx);
			}else // only one segment
				mate_clip_aln_item = NULL;

			consist_flag = computeBPConsistencyFlag(clip_pos_item->clip_aln, clip_pos_item->clip_end, mate_clip_aln_item, mate_clip_end, clip_pos_vec_main);
			if(consist_flag==1) {
				clip_pos_vec_main.push_back(clip_pos_item);
				clip_pos_vec.erase(clip_pos_vec.begin()+i);
			}else if(consist_flag==2){
				clip_pos_vec.erase(clip_pos_vec.begin()+i);
			}else i++;
		}
	}
}

// get clipping position items with same align segments by queryname
vector<clipPos_t*> clipReg::getClipPosItemsByQueryname(string &queryname, vector<clipPos_t*> &clip_pos_vec){
	vector<clipPos_t*> clip_items_ret;
	clipPos_t *clip_pos_item;

	for(size_t i=0; i<clip_pos_vec.size(); i++){
		clip_pos_item = clip_pos_vec.at(i);
		if(clip_pos_item->clip_aln->queryname.compare(queryname)==0) clip_items_ret.push_back(clip_pos_item);
	}

	return clip_items_ret;
}

clipPos_t* clipReg::getMinDistClipPosItem(clipPos_t *clip_pos, vector<clipPos_t*> &clip_pos_vec){
	clipPos_t *clip_pos_ret = NULL, *clip_pos_item;
	int64_t dist, minIdx, minValue;

	minIdx = -1;
	minValue = INT_MAX;
	for(size_t i=0; i<clip_pos_vec.size(); i++){
		clip_pos_item = clip_pos_vec.at(i);
		if(clip_pos_item->clip_aln!=clip_pos->clip_aln and clip_pos_item->clip_end!=clip_pos->clip_end){
			dist = clip_pos_item->clipRefPos - clip_pos->clipRefPos;
			if(dist<0) dist = -dist;
			if(dist<minValue){
				minValue = dist;
				minIdx = i;
			}
		}
	}
	if(minIdx!=-1) clip_pos_ret = clip_pos_vec.at(minIdx);

	return clip_pos_ret;
}

// get the vector id of given clip position item
int32_t clipReg::getClipPosVecId(clipAlnData_t *clip_aln, int32_t clip_end){
	int32_t vec_id = -1;
	clipPos_t *clip_pos;

	clip_pos = getClipPosItemFromSingleVec(clip_aln, clip_end, leftClipPosVector);
	if(clip_pos) vec_id = 0;

	if(vec_id==-1){
		clip_pos = getClipPosItemFromSingleVec(clip_aln, clip_end, rightClipPosVector);
		if(clip_pos) vec_id = 2;
	}
	if(vec_id==-1){
		clip_pos = getClipPosItemFromSingleVec(clip_aln, clip_end, leftClipPosVector2);
		if(clip_pos) vec_id = 1;
	}
	if(vec_id==-1){
		clip_pos = getClipPosItemFromSingleVec(clip_aln, clip_end, rightClipPosVector2);
		if(clip_pos) vec_id = 3;
	}

	return vec_id;
}

// get clipping position item from single vector
clipPos_t* clipReg::getClipPosItemFromSingleVec(clipAlnData_t *clip_aln, int32_t clip_end, vector<clipPos_t*> &clip_pos_vec){
	clipPos_t *clip_pos_ret = NULL, *clip_pos;
	for(size_t i=0; i<clip_pos_vec.size(); i++){
		clip_pos = clip_pos_vec.at(i);
		if(clip_pos->clip_aln==clip_aln and clip_pos->clip_end==clip_end){ // compare the pointer
			clip_pos_ret = clip_pos;
			break;
		}else if(clip_pos->clip_aln->chrname.compare(clip_aln->chrname)==0 and clip_pos->clip_aln->queryname.compare(clip_aln->queryname)==0
				and clip_pos->clip_aln->startRefPos==clip_aln->startRefPos and clip_pos->clip_aln->endRefPos==clip_aln->endRefPos and clip_pos->clip_end==clip_end){
			// compare the content
			clip_pos_ret = clip_pos;
			break;
		}
	}
	return clip_pos_ret;
}

// append clipping positions
void clipReg::appendClipPos(){
	//appendClipPosSingleVec(leftClipPosVector);
	appendClipPosSingleVec(leftClipPosVector2, clipAlnDataVector2);
	appendClipPosSingleVec(rightClipPosVector, rightClipAlnDataVector);
	appendClipPosSingleVec(rightClipPosVector2, rightClipAlnDataVector2);
}

// append clipping positions for single vector
void clipReg::appendClipPosSingleVec(vector<clipPos_t*> &clipPosVec, vector<clipAlnData_t*> &clipAlnDataVec){
	//vector<clipAlnData_t*> clip_aln_data_vec;
	int64_t meanClipPos, start_pos, end_pos, chr_len;
	string chrname_max;
	vector<string> clip_pos_str_vec;

	if(clipPosVec.size()){
		// construct the region
		clip_pos_str_vec = getClipPosMaxOcc(clipPosVec);
		if(clip_pos_str_vec.size()){
			chrname_max = clip_pos_str_vec.at(0);
			meanClipPos = stoi(clip_pos_str_vec.at(1));
			start_pos = meanClipPos - CLIP_END_EXTEND_SIZE;
			end_pos = meanClipPos + CLIP_END_EXTEND_SIZE;

			chr_len = faidx_seq_len(fai, chrname_max.c_str()); // get the reference length
			if(start_pos<1) start_pos = 1;
			if(end_pos>chr_len) end_pos = chr_len;

			// load the clipping data
			clipAlnDataLoader clip_aln_data_loader(chrname_max, start_pos, end_pos, inBamFile, minClipEndSize);
			clip_aln_data_loader.loadClipAlnDataWithSATag(clipAlnDataVec, paras->max_ultra_high_cov);
			removeNonclipItemsOp(clipAlnDataVec);

			// add clipping end at the specified region
			addClipPosAtClipEnd(chrname_max, start_pos, end_pos, RIGHT_END, clipPosVec, clipAlnDataVec, clipAlnDataVector);
			addClipPosAtClipEnd(chrname_max, start_pos, end_pos, LEFT_END, clipPosVec, clipAlnDataVec, clipAlnDataVector);
		}
	}
}

// get the clipping position of the maximum occurrence
vector<string> clipReg::getClipPosMaxOcc(vector<clipPos_t*> &clipPosVec){
	vector<string> clipPos_ret;

	struct groupNode{
		int32_t group_id, num;
		double central, sum;
		string chrname;
	};

	size_t i, j;
	double dist, min_dist;
	int32_t min_group_id;
	vector<struct groupNode*> group_vec;
	struct groupNode *group_item;
	clipPos_t *clip_pos_item;

	int64_t maxValue, secValue;
	int32_t maxIdx, secIdx;

	// extract groups
	if(clipPosVec.size()>0){
		for(i=0; i<clipPosVec.size(); i++){
			clip_pos_item = clipPosVec.at(i);

			// get the minimal distance
			min_dist = INT_MAX;
			min_group_id = -1;
			for(j=0; j<group_vec.size(); j++){
				group_item = group_vec.at(j);
				if(group_item->chrname.compare(clip_pos_item->chrname)==0){
					dist = clip_pos_item->clipRefPos - group_item->central;
					if(dist<0) dist = -dist;
					if(min_dist>dist) {
						min_dist = dist;
						min_group_id = j;
					}
				}
			}

			if(min_group_id==-1 or min_dist>GROUP_DIST_THRES){ // add a group
				group_item = new struct groupNode();
				group_item->group_id = group_vec.size();
				group_item->central = group_item->sum = clip_pos_item->clipRefPos;
				group_item->num = 1;
				group_item->chrname = clip_pos_item->chrname;
				group_vec.push_back(group_item);
			}else{ // found the minimum
				group_item = group_vec.at(min_group_id);
				group_item->sum += clip_pos_item->clipRefPos;
				group_item->num ++;
				group_item->central = group_item->sum / group_item->num;
			}
		}

		// get maximum item and second maximum item
		maxIdx = secIdx = -1; maxValue = secValue = 0;
		for(i=0; i<group_vec.size(); i++){
			group_item = group_vec.at(i);
			if(maxValue<group_item->num) { secIdx = maxIdx; secValue = maxValue; maxIdx = i; maxValue = group_item->num; }
			else if(secValue<group_item->num) { secIdx = i; secValue = group_item->num; }
		}

		// return values
		if(maxIdx>=0){
			group_item = group_vec.at(maxIdx);
			clipPos_ret.push_back(group_item->chrname);
			maxValue = group_item->central;
			clipPos_ret.push_back(to_string(maxValue));
		}

		for(i=0; i<group_vec.size(); i++) delete group_vec.at(i);
		vector<struct groupNode*>().swap(group_vec);
	}

	return clipPos_ret;
}

// add clipping position at clipping end
void clipReg::addClipPosAtClipEnd(string chrname_clip, int64_t start_pos, int64_t end_pos, int32_t clip_end, vector<clipPos_t*> &clipPosVec, vector<clipAlnData_t*> &clip_aln_data_vec, vector<clipAlnData_t*> &mainClipAlnDataVec){
	size_t i;
	int32_t vec_id;
	clipAlnData_t *clip_aln, *clip_aln_tmp;
	clipPos_t *clip_pos;

	for(i=0; i<clip_aln_data_vec.size(); i++){
		clip_aln = clip_aln_data_vec.at(i);

//		if(clip_aln->queryname.compare("S1_11995")==0){
//			cout << clip_aln->queryname << endl;
//		}

		if(clip_aln->chrname.compare(chrname_clip)==0){
			if(clip_end==RIGHT_END){ // right end
				if(clip_aln->rightClipSize>=minClipEndSize and clip_aln->endRefPos>=start_pos and clip_aln->endRefPos<=end_pos){
					vec_id = getClipPosVecId(clip_aln, clip_end);

					// further to search the main clipping data vector
					if(vec_id==-1){ // not exist, then add new item
						clip_aln_tmp = getClipAlnItemFromMainVec(clip_aln, mainClipAlnDataVec);
						if(clip_aln_tmp==NULL) {
							clip_aln_tmp = clip_aln;
							//cout << "xxxxxxxxxxx exist in main clipping data vector: " << clip_aln->chrname << ":" << clip_aln->startRefPos << "-" << clip_aln->endRefPos << endl;
						}

						clip_pos = new clipPos_t();
						clip_pos->chrname = clip_aln_tmp->chrname;
						clip_pos->clip_end = RIGHT_END;
						clip_pos->clip_aln = clip_aln_tmp;
						clip_pos->clipRefPos = clip_aln_tmp->endRefPos;
						clip_pos->clipQueryPos = clip_aln_tmp->endQueryPos;
						clip_pos->aln_orient = clip_aln_tmp->aln_orient;
						if(clip_aln_tmp->right_aln and clip_aln_tmp->aln_orient==clip_aln_tmp->right_aln->aln_orient) clip_pos->same_orient_flag = true;
						else clip_pos->same_orient_flag = false;
						clip_pos->self_overlap_flag = isSegSelfOverlap(clip_aln_tmp, clip_aln_tmp->right_aln, maxVarRegSize);
						clipPosVec.push_back(clip_pos);

						//cout << "------- Right clipping end: " << clip_pos->chrname << ":" << clip_pos->clipRefPos << endl;

					}else{
						//cout << "******** already exist, right clipping end: " << clip_aln->chrname << ":" << clip_aln->endRefPos << endl;
					}
				}
			}else{ // left end
				if(clip_aln->leftClipSize>=minClipEndSize and clip_aln->startRefPos>=start_pos and clip_aln->startRefPos<=end_pos){
					vec_id = getClipPosVecId(clip_aln, clip_end);

					// further to search the main clipping data vector
					if(vec_id==-1){ // not exist, then add new item
						clip_aln_tmp = getClipAlnItemFromMainVec(clip_aln, mainClipAlnDataVec);
						if(clip_aln_tmp==NULL) {
							clip_aln_tmp = clip_aln;
							//cout << "xxxxxxxxxxx exist in main clipping data vector: " << clip_aln->chrname << ":" << clip_aln->startRefPos << "-" << clip_aln->endRefPos << endl;
						}

						clip_pos = new clipPos_t();
						clip_pos->chrname = clip_aln->chrname;
						clip_pos->clip_end = LEFT_END;
						clip_pos->clip_aln = clip_aln;
						clip_pos->clipRefPos = clip_aln->startRefPos;
						clip_pos->clipQueryPos = clip_aln->startQueryPos;
						clip_pos->aln_orient = clip_aln->aln_orient;
						if(clip_aln->left_aln and clip_aln->aln_orient==clip_aln->left_aln->aln_orient) clip_pos->same_orient_flag = true;
						else clip_pos->same_orient_flag = false;
						clip_pos->self_overlap_flag = isSegSelfOverlap(clip_aln, clip_aln->left_aln, maxVarRegSize);
						clipPosVec.push_back(clip_pos);

						//cout << "------- Left clipping end: " << clip_pos->chrname << ":" << clip_pos->clipRefPos << endl;
					}else{
						//cout << "******** already exist, left clipping end: " << clip_aln->chrname << ":" << clip_aln->startRefPos << endl;
					}
				}
			}
		}
	}
}

// search the align item from vector
clipAlnData_t* clipReg::getClipAlnItemFromMainVec(clipAlnData_t *clip_aln, vector<clipAlnData_t*> &mainClipAlnDataVec){
	clipAlnData_t *clip_aln_ret = NULL, *clip_aln_tmp;

	for(size_t i=0; i<mainClipAlnDataVec.size(); i++){
		clip_aln_tmp = mainClipAlnDataVec.at(i);
		if(clip_aln_tmp->chrname.compare(clip_aln->chrname)==0 and clip_aln_tmp->queryname.compare(clip_aln->queryname)==0 and clip_aln_tmp->startRefPos==clip_aln->startRefPos and clip_aln_tmp->endRefPos==clip_aln->endRefPos){
			clip_aln_ret = clip_aln_tmp;
			break;
		}
	}

	return clip_aln_ret;
}

// sort clipping position items ascending
void clipReg::sortClipPos(){
	sortClipPosSingleVec(leftClipPosVector);
	sortClipPosSingleVec(leftClipPosVector2);
	sortClipPosSingleVec(rightClipPosVector);
	sortClipPosSingleVec(rightClipPosVector2);
}

// sort clipping position items ascending
void clipReg::sortClipPosSingleVec(vector<clipPos_t*> &clipPosVec){
	size_t i, j, minPos, minIdx, num;
	string tmp_str;
	clipPos_t *clip_pos_tmp;

	// sort vector ascending
	for(i=0; i<clipPosVec.size(); i++){
		minPos = clipPosVec.at(i)->clipRefPos;
		minIdx = i;
		for(j=i+1; j<clipPosVec.size(); j++){
			num = clipPosVec.at(j)->clipRefPos;
			if(num<minPos){
				minPos = num;
				minIdx = j;
			}
		}

		if(minIdx!=i){
			clip_pos_tmp = clipPosVec.at(i);
			clipPosVec.at(i) = clipPosVec.at(minIdx);
			clipPosVec.at(minIdx) = clip_pos_tmp;
		}
	}
}

// remove fake clips
void clipReg::removeFakeClips(){
	string vec_name;

	removeFakeClipsDifferentChrSingleVec(leftClipPosVector);
	removeFakeClipsDifferentChrSingleVec(leftClipPosVector2);
	removeFakeClipsDifferentChrSingleVec(rightClipPosVector);
	removeFakeClipsDifferentChrSingleVec(rightClipPosVector2);

	removeFakeClipsLongDistSameOrientSingleVec(leftClipPosVector, vec_name="Left");
	removeFakeClipsLongDistSameOrientSingleVec(leftClipPosVector2, vec_name="Left2");
	removeFakeClipsLongDistSameOrientSingleVec(rightClipPosVector, vec_name="Right");
	removeFakeClipsLongDistSameOrientSingleVec(rightClipPosVector2, vec_name="Right2");

	removeFakeClipsLowCov(leftClipPosVector, int32_t(paras->minReadsNumSupportSV*READS_NUM_SUPPORT_FACTOR));
	removeFakeClipsLowCov(leftClipPosVector2, int32_t(paras->minReadsNumSupportSV*READS_NUM_SUPPORT_FACTOR));
	removeFakeClipsLowCov(rightClipPosVector, int32_t(paras->minReadsNumSupportSV*READS_NUM_SUPPORT_FACTOR));
	removeFakeClipsLowCov(rightClipPosVector2, int32_t(paras->minReadsNumSupportSV*READS_NUM_SUPPORT_FACTOR));
}

// remove fake clips with alignment onto different chromes based on single vector
void clipReg::removeFakeClipsDifferentChrSingleVec(vector<clipPos_t*> &clipPosVector){
	size_t i, j, sub_sum, maxValue;
	int32_t maxIdx;
	vector<string> chrname_vec;
	string chrname_tmp;
	vector<size_t> num_vec;

	// get chrnames
	for(i=0; i<clipPosVector.size(); i++)
		if(isExistStr(clipPosVector.at(i)->chrname, chrname_vec)==false) // not exist in vector, add it to vector
			chrname_vec.push_back(clipPosVector.at(i)->chrname);

	// compute the number of clips for each chrome
	for(i=0; i<chrname_vec.size(); i++){
		chrname_tmp = chrname_vec.at(i);
		sub_sum = 0;
		for(j=0; j<clipPosVector.size(); j++) if(chrname_tmp.compare(clipPosVector.at(j)->chrname)==0) sub_sum ++;
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
			if(chrname_tmp.compare(clipPosVector.at(i)->chrname)!=0){
				delete clipPosVector.at(i);
				clipPosVector.erase(clipPosVector.begin()+i);
			}else i++;
		}
	}
}

// remove fake clips with long dist based on single vector
void clipReg::removeFakeClipsLongDistSameOrientSingleVec(vector<clipPos_t*> &clipPosVector, string &vec_name){
	size_t i, group_id, group_id_max;
	int32_t dist;
	vector<size_t> clip_pos_group_id_vec, group_central_vec, group_num_vec;
	clipPos_t *clip_pos_item;

	size_t num, maxValue, secValue;
	int32_t maxIdx, secIdx;
	double sum, central_value, pre_central_value;

	// extract groups
	if(clipPosVector.size()>0){
		sum = central_value = pre_central_value = 0;
		group_id = 0;
		for(i=0; i<clipPosVector.size(); i++){
			clip_pos_item = clipPosVector.at(i);
			if(central_value==0){
				sum = central_value = clip_pos_item->clipRefPos;
				num = 1;
				clip_pos_group_id_vec.push_back(group_id);
			}else{
				dist = clip_pos_item->clipRefPos - central_value;
				if(dist<0){
					cerr << "line=" << __LINE__ << ", i=" << i << ", dist=" << dist << ", the vector was not sorted in ascending order, error!" << endl;
					exit(1);
				}

				if(dist<=GROUP_DIST_THRES){
					sum += clip_pos_item->clipRefPos;
					num ++;
					central_value = sum / num;
					clip_pos_group_id_vec.push_back(group_id);
				}else{
					group_central_vec.push_back(central_value);
					group_num_vec.push_back(num);
					pre_central_value = central_value;
					sum = central_value = clip_pos_item->clipRefPos;
					num = 1;
					clip_pos_group_id_vec.push_back(++group_id);
				}
			}
		}
		if(num>0){ // final item
			group_central_vec.push_back(central_value);
			group_num_vec.push_back(num);
		}

		// get maximum item and second maximum item
		maxIdx = secIdx = -1; maxValue = secValue = 0;
		for(i=0; i<group_num_vec.size(); i++){
			num = group_num_vec.at(i);
			if(maxValue<num) { secIdx = maxIdx; secValue = maxValue; maxIdx = i; maxValue = num; }
			else if(secValue<num) { secIdx = i; secValue = num; }
		}

		// remove false clips
		if(maxIdx>=0){
			group_id_max = maxIdx;
			for(i=0; i<clipPosVector.size(); ){
				group_id = clip_pos_group_id_vec.at(i);
				if(group_id!=group_id_max){ 	// remove non-maximum items
					delete clipPosVector.at(i);
					clipPosVector.erase(clipPosVector.begin()+i);
					clip_pos_group_id_vec.erase(clip_pos_group_id_vec.begin()+i);
				}else i++;
			}
		}
	}
}

// remove fake clips with long dist based on single vector
void clipReg::removeFakeClipsLowCov(vector<clipPos_t*> &clipPosVector, int32_t min_clip_reads_num){
	if(clipPosVector.size()<(size_t)min_clip_reads_num) destroyClipPosVec(clipPosVector); // remove fake clips Only this function modifies the minimum number of supports
}

// get the index
int32_t clipReg::getItemIdxClipPosVec(clipPos_t *item, vector<clipPos_t*> &vec){
	int32_t idx = -1;
	for(size_t i=0; i<vec.size(); i++)
		if(item->chrname.compare(vec.at(i)->chrname)==0 and item->clipRefPos==vec.at(i)->clipRefPos){
			idx = i;
			break;
		}
	return idx;
}

// compute single mate region for BND format
string clipReg::computeMateSingleRegStrForBND(vector<clipPos_t*> &clip_pos_vec, int32_t vec_id, string &chrname_clip, int64_t meanClipPos){
	string mate_reg_str, left_mate_str, right_mate_str, vec_id_str, orient_str, BND_str, mate_clip_pos_str;
	size_t i;
	clipPos_t *clip_pos;
	clipAlnData_t *clip_aln, *clip_aln2;
	int32_t j, clip_end2, vec_id2, num_arr[4], orient_num_arr[2][4];
	int32_t maxValue, maxIdx, support_num, cov_num;

	if(clip_pos_vec.size()==0) return "-";

	// coverage of right end clipping position
	cov_num = computeCovNumClipPos(chrname_clip, meanClipPos, RIGHT_END, fai, paras);

	// initialization
	support_num = 0;
	for(j=0; j<4; j++) num_arr[j] = orient_num_arr[0][j] = orient_num_arr[1][j] = 0;

	// right end
	for(i=0; i<clip_pos_vec.size(); i++){
		clip_pos = clip_pos_vec.at(i);
		clip_aln = clip_pos->clip_aln;

		// get the neighboring clipping segments 2+,3-;2-,3+;1-,0+;1+,0-
		clip_aln2 = NULL;
		clip_end2 = -1;
		if(clip_pos->clip_end==RIGHT_END){ // right end
			support_num ++;
			clip_aln2 = clip_aln->right_aln;
			if(clip_aln2){
				if(clip_aln2->aln_orient==clip_aln->aln_orient) // same orient
					clip_end2 = LEFT_END;
				else // different orient
					clip_end2 = RIGHT_END;
			}
		}

		// get vector id
		if(clip_aln2){
			vec_id2 = getClipPosVecId(clip_aln2, clip_end2);
			if(vec_id2!=-1){ // item exists
				if(vec_id2!=vec_id){
					num_arr[vec_id2] ++;
					if(clip_end2==LEFT_END){ // same orient
						orient_num_arr[0][vec_id2] ++;
						BND_str = "N[" + clip_aln2->chrname + ":" + to_string(clip_aln2->startRefPos) + "[";
						//cout << "Vec[" << vec_id << "]: qname=" << clip_aln->queryname << ", " << clip_aln->chrname << ":" << clip_aln->endRefPos << ", mate item: " << BND_str << endl;
					}else { // different orient
						orient_num_arr[1][vec_id2] ++;
						BND_str = "N]" + clip_aln2->chrname + ":" + to_string(clip_aln2->endRefPos) + "]";
						//cout << "Vec[" << vec_id << "]: qname=" << clip_aln->queryname << ", " << clip_aln->chrname << ":" << clip_aln->endRefPos << ", mate item: " << BND_str << endl;
					}
				}else{
					//cerr << __func__ << ", line=" << __LINE__ <<  ": invalid same vector id, vec_id=" << vec_id << ", vec_id2=" << vec_id2 << ", clip_vec_size=" << clip_pos_vec.size() << ", error!" << endl;
					//exit(1);
				}
			}
		}
	}

//	cout << "Right end: id_nums[" << vec_id << "]: ";
//	for(j=0; j<3; j++) cout << num_arr[j] << ", ";
//	cout << num_arr[3] << endl;

	// get the final vector id for right clipping end
	maxValue = maxIdx = -1;
	for(j=0; j<4; j++){
		if(num_arr[j]>0 and num_arr[j]>maxValue) {
			maxValue = num_arr[j];
			maxIdx = j;
		}
	}
	//cout << "maxIdx=" << maxIdx << ", maxValue=" << maxValue << ", support_num=" << support_num << endl;

	//if(maxIdx!=-1){
	if(maxIdx!=-1 and maxValue>=MIN_CLIP_REG_MATED_RATIO*support_num){
		//cout << "maxIdx=" << maxIdx << ", maxValue=" << maxValue << ", orient_nums=[" << orient_num_arr[0][maxIdx] << "," << orient_num_arr[1][maxIdx] << "]" << endl;

		vec_id_str = to_string(maxIdx);
		if(orient_num_arr[0][maxIdx]>orient_num_arr[1][maxIdx]){ // mate at left end
			orient_str = "+";
			mate_clip_pos_str = to_string(meanClipPos);
		}else{ // mate at right end
			orient_str = "-";
			mate_clip_pos_str = to_string(meanClipPos-1);
		}
		if(support_num>cov_num) cov_num = support_num;
		left_mate_str = vec_id_str + orient_str + "|" + mate_clip_pos_str + "|" + to_string(support_num) + "|" + to_string(cov_num);
	}else
		left_mate_str = "-";

	// coverage of left end clipping position
	cov_num = computeCovNumClipPos(chrname_clip, meanClipPos, LEFT_END, fai, paras);

	// initialization
	support_num = 0;
	for(j=0; j<4; j++) num_arr[j] = orient_num_arr[0][j] = orient_num_arr[1][j] = 0;

	// left end
	for(i=0; i<clip_pos_vec.size(); i++){
		clip_pos = clip_pos_vec.at(i);
		clip_aln = clip_pos->clip_aln;

		// get the neighboring clipping segments 2+|.|.|.,3-.;2-,3+;1-,0+;1+,0-
		clip_aln2 = NULL;
		clip_end2 = -1;
		if(clip_pos->clip_end==LEFT_END){ // left end
			support_num ++;
			clip_aln2 = clip_aln->left_aln;
			if(clip_aln2){
				if(clip_aln2->aln_orient==clip_aln->aln_orient) // same orient
					clip_end2 = RIGHT_END;
				else // different orient
					clip_end2 = LEFT_END;
			}
		}

		// get vector id
		if(clip_aln2){
			vec_id2 = getClipPosVecId(clip_aln2, clip_end2);
			if(vec_id2!=-1){ // item exists
				if(vec_id2!=vec_id){
					num_arr[vec_id2] ++;
					if(clip_end2==RIGHT_END) { // same orient
						orient_num_arr[0][vec_id2] ++;
						BND_str = "]" + clip_aln2->chrname + ":" + to_string(clip_aln2->endRefPos) + "]N";
						//cout << "Vec[" << vec_id << "]: qname=" << clip_aln->queryname << ", " << clip_aln->chrname << ":" << clip_aln->startRefPos << ", mate item: " << BND_str << endl;
					}else{ // different orient
						orient_num_arr[1][vec_id2] ++;
						BND_str = "[" + clip_aln2->chrname + ":" + to_string(clip_aln2->startRefPos) + "[N";
						//cout << "Vec[" << vec_id << "]: qname=" << clip_aln->queryname << ", " << clip_aln->chrname << ":" << clip_aln->startRefPos << ", mate item: " << BND_str << endl;
					}
				}else{
					//cerr << __func__ << ", line=" << __LINE__ <<  ": invalid same vector id, vec_id=" << vec_id << ", vec_id2=" << vec_id2 << ", clip_vec_size=" << clip_pos_vec.size() << ", error!" << endl;
					//exit(1);
				}
			}
		}
	}

//	cout << "Left end: id_nums[" << vec_id << "]: ";
//	for(j=0; j<3; j++) cout << num_arr[j] << ", ";
//	cout << num_arr[3] << endl;

	// get the final vector id for left clipping end
	maxValue = maxIdx = -1;
	for(j=0; j<4; j++){
		if(num_arr[j]>0 and num_arr[j]>maxValue) {
			maxValue = num_arr[j];
			maxIdx = j;
		}
	}

	//cout << "maxIdx=" << maxIdx << ", maxValue=" << maxValue << ", support_num=" << support_num << endl;

	//if(maxIdx!=-1){
	if(maxIdx!=-1 and maxValue>=MIN_CLIP_REG_MATED_RATIO*support_num){
		//cout << "maxIdx=" << maxIdx << ", maxValue=" << maxValue << ", orient_nums=[" << orient_num_arr[0][maxIdx] << "," << orient_num_arr[1][maxIdx] << "]" << endl;

		vec_id_str = to_string(maxIdx);
		if(orient_num_arr[0][maxIdx]>orient_num_arr[1][maxIdx]){ // mate at right end
			orient_str =  "+";
			mate_clip_pos_str = to_string(meanClipPos-1);
		}else{ // mate at left end
			orient_str =  "-";
			mate_clip_pos_str = to_string(meanClipPos);
		}
		if(support_num>cov_num) cov_num = support_num;
		right_mate_str = vec_id_str + orient_str + "|" + mate_clip_pos_str + "|" + to_string(support_num) + "|" + to_string(cov_num);
	}else
		right_mate_str = "-";

	if(left_mate_str.compare("-")==0 and right_mate_str.compare("-")==0) mate_reg_str = "-";
	else mate_reg_str = left_mate_str + "," + right_mate_str;

	//cout << mate_reg_str << endl;

	return mate_reg_str;
}

// compute the coverage of a clipping position
int32_t clipReg::computeCovNumClipPos(string &chrname, int64_t meanClipPos, int32_t clip_end, faidx_t *fai, Paras *paras){
	int64_t start_pos, end_pos, chr_len, pos, maxValue, num;
	Base *baseArray;
	vector<bam1_t*> alnDataVector;
	vector<bam1_t*>::iterator aln;

	if(clip_end==RIGHT_END) {
		start_pos = meanClipPos - CLIP_END_EXTEND_SIZE / 2;
		end_pos = meanClipPos - 1;
	}else{
		start_pos = meanClipPos;
		end_pos = meanClipPos + CLIP_END_EXTEND_SIZE / 2;
	}
	chr_len = faidx_seq_len(fai, chrname.c_str()); // get the reference length
	if(start_pos<1) start_pos = 1;
	if(end_pos>chr_len) end_pos = chr_len;

	covLoader cov_loader(chrname, start_pos, end_pos, fai, paras->min_ins_size_filt, paras->min_del_size_filt);
	baseArray = cov_loader.initBaseArray();
	alnDataLoader data_loader(chrname, start_pos, end_pos, paras->inBamFile);
	data_loader.loadAlnData(alnDataVector);

	// generate the base coverage array
	cov_loader.generateBaseCoverage(baseArray, alnDataVector);

	//totalReadBeseNum = totalRefBaseNum = 0;
	maxValue = 0;
	for(pos=start_pos; pos<=end_pos; pos++){
		// compute the meanCov excluding the gap regions
		if(baseArray[pos-start_pos].coverage.idx_RefBase!=4){ // excluding 'N'
//			totalReadBeseNum += baseArray[pos-start_pos].coverage.num_bases[5];
//			totalRefBaseNum ++;
			num = baseArray[pos-start_pos].coverage.num_bases[5];
			if(maxValue<num) maxValue = num;
		}
	}
//	if(totalRefBaseNum) mean_cov_num = (double)totalReadBeseNum/totalRefBaseNum;
//	else mean_cov_num = 0;

	// release memory
	if(baseArray) { delete[] baseArray; baseArray = NULL; }
	if(!alnDataVector.empty()){
		for(aln=alnDataVector.begin(); aln!=alnDataVector.end(); aln++)
			bam_destroy1(*aln);
		vector<bam1_t*>().swap(alnDataVector);
	}

	return maxValue;
}

// compute the clip region
void clipReg::computeClipRegs(){
	reg_t *reg_tmp;
	int64_t val_tmp, dist;
	size_t i;
	bool merge_flag;

	for(i=0; i<4; i++) mate_clip_reg.bnd_mate_reg_strs[i] = "-";

	// left part
	mate_clip_reg.leftMeanClipPos = computeMeanClipPos(leftClipPosVector);
	mate_clip_reg.leftMeanClipPos2 = computeMeanClipPos(leftClipPosVector2);
	mate_clip_reg.leftClipPosNum = leftClipPosVector.size();
	mate_clip_reg.leftClipPosNum2 = leftClipPosVector2.size();
	merge_flag = false;
	if(mate_clip_reg.leftClipPosNum>0 and mate_clip_reg.leftClipPosNum2>0){
		dist = mate_clip_reg.leftMeanClipPos2 - mate_clip_reg.leftMeanClipPos;
		if(dist<0) dist = -dist;
		if(dist<MAX_CLIP_REG_MERGE_DIST_ADJUST){ // short distance, then merge
			for(i=0; i<leftClipPosVector2.size(); i++) leftClipPosVector.push_back(leftClipPosVector2.at(i));
			vector<clipPos_t*>().swap(leftClipPosVector2);
			merge_flag = true;
		}
	}
	if(mate_clip_reg.leftClipPosNum==0 and mate_clip_reg.leftClipPosNum2>0){ // move items of vector2 to vector1
		for(i=0; i<leftClipPosVector2.size(); i++) leftClipPosVector.push_back(leftClipPosVector2.at(i));
		vector<clipPos_t*>().swap(leftClipPosVector2);
		mate_clip_reg.leftMeanClipPos = mate_clip_reg.leftMeanClipPos2;
		mate_clip_reg.leftMeanClipPos2 = 0;
		mate_clip_reg.leftClipPosNum = leftClipPosVector.size();
		mate_clip_reg.leftClipPosNum2 = leftClipPosVector2.size();
	}

	mate_clip_reg.leftClipReg = computeClipRegSingleVec(leftClipPosVector);
	mate_clip_reg.leftClipReg2 = computeClipRegSingleVec(leftClipPosVector2);
	if(merge_flag){
		mate_clip_reg.leftMeanClipPos = computeMeanClipPos(leftClipPosVector);
		mate_clip_reg.leftMeanClipPos2 = computeMeanClipPos(leftClipPosVector2);
	}
	mate_clip_reg.leftClipPosNum = leftClipPosVector.size();
	mate_clip_reg.leftClipPosNum2 = leftClipPosVector2.size();
	if(mate_clip_reg.leftClipReg and mate_clip_reg.leftClipReg2 and mate_clip_reg.leftMeanClipPos>mate_clip_reg.leftMeanClipPos2){
		reg_tmp = mate_clip_reg.leftClipReg; mate_clip_reg.leftClipReg = mate_clip_reg.leftClipReg2; mate_clip_reg.leftClipReg2 = reg_tmp;
		val_tmp = mate_clip_reg.leftMeanClipPos; mate_clip_reg.leftMeanClipPos = mate_clip_reg.leftMeanClipPos2; mate_clip_reg.leftMeanClipPos2 = val_tmp;
		val_tmp = mate_clip_reg.leftClipPosNum; mate_clip_reg.leftClipPosNum = mate_clip_reg.leftClipPosNum2; mate_clip_reg.leftClipPosNum2 = val_tmp;
		left_part_changed = true;
	}

	// right part
	mate_clip_reg.rightMeanClipPos = computeMeanClipPos(rightClipPosVector);
	mate_clip_reg.rightMeanClipPos2 = computeMeanClipPos(rightClipPosVector2);
	mate_clip_reg.rightClipPosNum = rightClipPosVector.size();
	mate_clip_reg.rightClipPosNum2 = rightClipPosVector2.size();
	merge_flag = false;
	if(mate_clip_reg.rightClipPosNum>0 and mate_clip_reg.rightClipPosNum2>0){
		dist = mate_clip_reg.rightMeanClipPos2 - mate_clip_reg.rightMeanClipPos;
		if(dist<0) dist = -dist;
		if(dist<MAX_CLIP_REG_MERGE_DIST_ADJUST){ // short distance, then merge
			for(i=0; i<rightClipPosVector2.size(); i++) rightClipPosVector.push_back(rightClipPosVector2.at(i));
			vector<clipPos_t*>().swap(rightClipPosVector2);
			merge_flag = true;
		}
	}
	if(mate_clip_reg.rightClipPosNum==0 and mate_clip_reg.rightClipPosNum2>0){ // move items of vector2 to vector1
		for(i=0; i<rightClipPosVector2.size(); i++) rightClipPosVector.push_back(rightClipPosVector2.at(i));
		vector<clipPos_t*>().swap(rightClipPosVector2);
		mate_clip_reg.rightMeanClipPos = mate_clip_reg.rightMeanClipPos2;
		mate_clip_reg.rightMeanClipPos2 = 0;
		mate_clip_reg.rightClipPosNum = rightClipPosVector.size();
		mate_clip_reg.rightClipPosNum2 = rightClipPosVector2.size();
	}

	mate_clip_reg.rightClipReg = computeClipRegSingleVec(rightClipPosVector);
	mate_clip_reg.rightClipReg2 = computeClipRegSingleVec(rightClipPosVector2);
	if(merge_flag){
		mate_clip_reg.rightMeanClipPos = computeMeanClipPos(rightClipPosVector);
		mate_clip_reg.rightMeanClipPos2 = computeMeanClipPos(rightClipPosVector2);
	}
	mate_clip_reg.rightClipPosNum = rightClipPosVector.size();
	mate_clip_reg.rightClipPosNum2 = rightClipPosVector2.size();
	if(mate_clip_reg.rightClipReg and mate_clip_reg.rightClipReg2 and mate_clip_reg.rightMeanClipPos>mate_clip_reg.rightMeanClipPos2){
		reg_tmp = mate_clip_reg.rightClipReg; mate_clip_reg.rightClipReg = mate_clip_reg.rightClipReg2; mate_clip_reg.rightClipReg2 = reg_tmp;
		val_tmp = mate_clip_reg.rightMeanClipPos; mate_clip_reg.rightMeanClipPos = mate_clip_reg.rightMeanClipPos2; mate_clip_reg.rightMeanClipPos2 = val_tmp;
		val_tmp = mate_clip_reg.rightClipPosNum; mate_clip_reg.rightClipPosNum = mate_clip_reg.rightClipPosNum2; mate_clip_reg.rightClipPosNum2 = val_tmp;
		right_part_changed = true;
	}

	// compute mate regions for BND format
	if(left_part_changed==false){
		if(mate_clip_reg.leftClipReg) mate_clip_reg.bnd_mate_reg_strs[0] = computeMateSingleRegStrForBND(leftClipPosVector, 0, mate_clip_reg.leftClipReg->chrname, mate_clip_reg.leftMeanClipPos);
		if(mate_clip_reg.leftClipReg2) mate_clip_reg.bnd_mate_reg_strs[1] = computeMateSingleRegStrForBND(leftClipPosVector2, 1, mate_clip_reg.leftClipReg2->chrname, mate_clip_reg.leftMeanClipPos2);
	}else{
		if(mate_clip_reg.leftClipReg) mate_clip_reg.bnd_mate_reg_strs[0] = computeMateSingleRegStrForBND(leftClipPosVector2, 0, mate_clip_reg.leftClipReg->chrname, mate_clip_reg.leftMeanClipPos);
		if(mate_clip_reg.leftClipReg2) mate_clip_reg.bnd_mate_reg_strs[1] = computeMateSingleRegStrForBND(leftClipPosVector, 1, mate_clip_reg.leftClipReg2->chrname, mate_clip_reg.leftMeanClipPos2);
	}
	if(right_part_changed==false){
		if(mate_clip_reg.rightClipReg) mate_clip_reg.bnd_mate_reg_strs[2] = computeMateSingleRegStrForBND(rightClipPosVector, 2, mate_clip_reg.rightClipReg->chrname, mate_clip_reg.rightMeanClipPos);
		if(mate_clip_reg.rightClipReg2) mate_clip_reg.bnd_mate_reg_strs[3] = computeMateSingleRegStrForBND(rightClipPosVector2, 3, mate_clip_reg.rightClipReg2->chrname, mate_clip_reg.rightMeanClipPos2);
	}else{
		if(mate_clip_reg.rightClipReg) mate_clip_reg.bnd_mate_reg_strs[2] = computeMateSingleRegStrForBND(rightClipPosVector2, 2, mate_clip_reg.rightClipReg->chrname, mate_clip_reg.rightMeanClipPos);
		if(mate_clip_reg.rightClipReg2) mate_clip_reg.bnd_mate_reg_strs[3] = computeMateSingleRegStrForBND(rightClipPosVector, 3, mate_clip_reg.rightClipReg2->chrname, mate_clip_reg.rightMeanClipPos2);
	}

	mate_clip_reg.leftClipRegNum = 0;
	if(mate_clip_reg.leftClipReg) mate_clip_reg.leftClipRegNum ++;
	if(mate_clip_reg.leftClipReg2) mate_clip_reg.leftClipRegNum ++;
	mate_clip_reg.rightClipRegNum = 0;
	if(mate_clip_reg.rightClipReg) mate_clip_reg.rightClipRegNum ++;
	if(mate_clip_reg.rightClipReg2) mate_clip_reg.rightClipRegNum ++;

	// remove unmated regions
	removeBNDUnmatedClipRegs();

	// determine large indel region
	determineLargeIndelReg();

	mate_clip_reg.var_cand = NULL;
	mate_clip_reg.left_var_cand_tra = mate_clip_reg.right_var_cand_tra = NULL;
	//for(int i=0; i<4; i++) mate_clip_reg.bnd_mate_reg_strs[i] = bnd_mate_reg_strs[i];

	// remove FP region for single clip end
	removeFPClipSingleEnd(mate_clip_reg);

	// check location order
	checkLocOrder(mate_clip_reg);

	mate_clip_reg.valid_flag = true;
	if((mate_clip_reg.leftClipReg or mate_clip_reg.leftClipReg2) and (mate_clip_reg.rightClipReg or mate_clip_reg.rightClipReg2)) mate_clip_reg.reg_mated_flag = true;

	// sv_type and dup_num
	mate_clip_reg.sv_type = VAR_UNC;
	mate_clip_reg.dup_num = 0;
	//if(mate_clip_reg.reg_mated_flag) // sv_type
		computeVarTypeClipReg(mate_clip_reg);

	if(mate_clip_reg.sv_type==VAR_UNC) mate_clip_reg.valid_flag = false;

	// check BND items
	checkBNDStrVec(mate_clip_reg);

	// check support num
	checkSuppNum(mate_clip_reg, int32_t(paras->minReadsNumSupportSV*READS_NUM_SUPPORT_FACTOR));
}

// determine large deletion region
void clipReg::determineLargeIndelReg(){
	reg_t *reg, *reg1, *reg2;
	vector<clipPos_t*> *clip_pos_vec1, *clip_pos_vec2, clip_pos_items_tmp1, clip_pos_items_tmp2;
	size_t i;
	vector<clipAlnData_t*> query_aln_segs;
	clipAlnData_t *mate_clip_aln_seg;
	string queryname, chrname_tmp;
	vector<int32_t> adjClipAlnSegInfo;
	clipPos_t *clip_pos;
	int32_t mate_arr_idx, mate_clip_end_flag, idx_vec, reg_num;
	int32_t mate_arr_idx_same_chr, mate_clip_end_flag_same_chr, min_dist_same_chr, valid_del_num;
	int32_t start_pos_tmp, end_pos_tmp;
	bool valid_mate_flag, same_orient_flag, self_overlap_flag, ref_skip_flag;

	// save the INS region if there are only one clipPos set
	reg_num = mate_clip_reg.leftClipRegNum + mate_clip_reg.rightClipRegNum;
	if(reg_num==1){ // only one clipPos set
		reg = NULL;
		if(mate_clip_reg.leftClipRegNum==1){ // left part
			if(mate_clip_reg.leftClipReg){
				reg = mate_clip_reg.leftClipReg;
				copyClipPosVec(leftClipPosVector, largeIndelClipPosVector);
			}else{
				reg = mate_clip_reg.leftClipReg2;
				copyClipPosVec(leftClipPosVector2, largeIndelClipPosVector);
			}
		}else{ // right part
			if(mate_clip_reg.rightClipReg){
				reg = mate_clip_reg.rightClipReg;
				copyClipPosVec(rightClipPosVector, largeIndelClipPosVector);
			}else{
				reg = mate_clip_reg.rightClipReg2;
				copyClipPosVec(rightClipPosVector2, largeIndelClipPosVector);
			}
		}

		// add the region into indel vector
		largeIndelClipReg = new reg_t();
		largeIndelClipReg->chrname = reg->chrname;
		largeIndelClipReg->startRefPos = reg->startRefPos;
		largeIndelClipReg->endRefPos = reg->endRefPos;
		largeIndelClipReg->startLocalRefPos = largeIndelClipReg->endLocalRefPos = 0;
		largeIndelClipReg->startQueryPos = largeIndelClipReg->endQueryPos = 0;
		largeIndelClipReg->var_type = VAR_UNC;
		largeIndelClipReg->sv_len = 0;
		largeIndelClipReg->query_id = -1;
		largeIndelClipReg->blat_aln_id = -1;
		largeIndelClipReg->minimap2_aln_id = -1;
		largeIndelClipReg->call_success_status = false;
		largeIndelClipReg->short_sv_flag = false;
		largeIndelClipReg->zero_cov_flag = false;
		largeIndelClipReg->aln_seg_end_flag = false;
		largeIndelClipReg->query_pos_invalid_flag = false;
		largeIndelClipReg->large_indel_flag = true;
		largeIndelClipReg->gt_type = -1;
		largeIndelClipReg->gt_seq = "";
		largeIndelClipReg->AF = 0;

		large_indel_flag = true;
		mate_clip_reg.large_indel_flag = true;
		mate_clip_reg.largeIndelClipReg = dupVarReg(largeIndelClipReg);
		mate_clip_reg.supp_num_largeIndel = supp_num_largeIndel;
	}
	else if(mate_clip_reg.leftClipRegNum==1 and mate_clip_reg.rightClipRegNum==1){
		if(leftClipPosVector.size()>0)clip_pos_vec1 = &leftClipPosVector;
		else clip_pos_vec1 = &leftClipPosVector2;

		resetClipCheckFlag(clipAlnDataVector);

		valid_del_num = 0;
		for(i=0; i<clip_pos_vec1->size(); i++){
			clip_pos = clip_pos_vec1->at(i);

			queryname = clip_pos->clip_aln->queryname;
			query_aln_segs = getQueryClipAlnSegsAll(queryname, clipAlnDataVector);  // get query clip align segments

			// deal with the mate clip end
			ref_skip_flag = false;
			idx_vec = getVecIdxClipAlnData(clip_pos->clip_aln, query_aln_segs);
			adjClipAlnSegInfo = getAdjacentClipAlnSeg(idx_vec, clip_pos->clip_end, query_aln_segs, minClipEndSize, MAX_VAR_REG_SIZE);
			mate_arr_idx = adjClipAlnSegInfo.at(0);
			mate_clip_end_flag = adjClipAlnSegInfo.at(1);
			mate_arr_idx_same_chr = adjClipAlnSegInfo.at(4);
			mate_clip_end_flag_same_chr = adjClipAlnSegInfo.at(5);
			min_dist_same_chr = adjClipAlnSegInfo.at(6);
			//min_ref_dist_same_chr = adjClipAlnSegInfo.at(7);
			if(mate_arr_idx!=-1){ // mated
				// prefer the segments on the same chromosome
				if(mate_arr_idx_same_chr!=-1 and mate_arr_idx==mate_arr_idx_same_chr and mate_clip_end_flag==mate_clip_end_flag_same_chr and abs(min_dist_same_chr)<=MAX_REF_DIST_SAME_CHR){
					ref_skip_flag = true;
				}

				if(ref_skip_flag){
					mate_clip_aln_seg = query_aln_segs.at(mate_arr_idx);

					if(clip_pos->clip_aln->aln_orient==mate_clip_aln_seg->aln_orient) same_orient_flag = true;
					else same_orient_flag = false;
					self_overlap_flag = isSegSelfOverlap(clip_pos->clip_aln, mate_clip_aln_seg, maxVarRegSize);
					clip_pos_items_tmp1 = getClipPosItemsByQueryname(clip_pos->clip_aln->queryname, *clip_pos_vec1);

					valid_mate_flag = false;
					if(clip_pos_items_tmp1.size()==1){ // only one segment in both vector, respectively
						clip_pos_vec2 = NULL;
						idx_vec = getClipPosVecId(mate_clip_aln_seg, mate_clip_end_flag);
						switch(idx_vec){
							case 0: clip_pos_vec2 = &leftClipPosVector; break;
							case 1: clip_pos_vec2 = &leftClipPosVector2; break;
							case 2: clip_pos_vec2 = &rightClipPosVector; break;
							case 3: clip_pos_vec2 = &rightClipPosVector2; break;
							//default:
								//printClipVecs("hahahah");
								//cerr << "invalid idx_vec=" << idx_vec << ", warning!" << endl;
								//exit(1);
						}
						if(clip_pos_vec2 and clip_pos_vec2!=clip_pos_vec1){
							clip_pos_items_tmp2 = getClipPosItemsByQueryname(clip_pos->clip_aln->queryname, *clip_pos_vec2);
							if(clip_pos_items_tmp2.size()==1) valid_mate_flag = true;
						}
					}

					if(same_orient_flag and self_overlap_flag==false and valid_mate_flag){
						valid_del_num ++;
					}
				}
			}
		}
		resetClipCheckFlag(clipAlnDataVector);

		// construct DEL region
		if(valid_del_num>paras->minReadsNumSupportSV*READS_NUM_SUPPORT_FACTOR){
			if(leftClipPosVector.size()>0){
				clip_pos_vec1 = &leftClipPosVector;
				reg1 = mate_clip_reg.leftClipReg;
				copyClipPosVec(leftClipPosVector, largeIndelClipPosVector);
			}else{
				clip_pos_vec1 = &leftClipPosVector2;
				reg1 = mate_clip_reg.leftClipReg2;
				copyClipPosVec(leftClipPosVector2, largeIndelClipPosVector);
			}

			if(rightClipPosVector.size()>0){
				clip_pos_vec2 = &rightClipPosVector;
				reg2 = mate_clip_reg.rightClipReg;
				copyClipPosVec(rightClipPosVector, largeIndelClipPosVector);
			}else{
				clip_pos_vec2 = &rightClipPosVector2;
				reg2 = mate_clip_reg.rightClipReg2;
				copyClipPosVec(rightClipPosVector2, largeIndelClipPosVector);
			}

			chrname_tmp = reg1->chrname;
			if(reg1->startRefPos<reg2->startRefPos) start_pos_tmp = reg1->startRefPos;
			else start_pos_tmp = reg2->startRefPos;
			if(reg1->endRefPos<reg2->endRefPos) end_pos_tmp = reg2->endRefPos;
			else end_pos_tmp = reg1->endRefPos;

			// add the DEL region into indel vector
			largeIndelClipReg = new reg_t();
			largeIndelClipReg->chrname = chrname_tmp;
			largeIndelClipReg->startRefPos = start_pos_tmp;
			largeIndelClipReg->endRefPos = end_pos_tmp;
			largeIndelClipReg->startLocalRefPos = largeIndelClipReg->endLocalRefPos = 0;
			largeIndelClipReg->startQueryPos = largeIndelClipReg->endQueryPos = 0;
			largeIndelClipReg->var_type = VAR_UNC;
			largeIndelClipReg->sv_len = 0;
			largeIndelClipReg->query_id = -1;
			largeIndelClipReg->blat_aln_id = -1;
			largeIndelClipReg->minimap2_aln_id = -1;
			largeIndelClipReg->call_success_status = false;
			largeIndelClipReg->short_sv_flag = false;
			largeIndelClipReg->zero_cov_flag = false;
			largeIndelClipReg->aln_seg_end_flag = false;
			largeIndelClipReg->query_pos_invalid_flag = false;
			largeIndelClipReg->large_indel_flag = true;
			largeIndelClipReg->gt_type = -1;
			largeIndelClipReg->gt_seq = "";
			largeIndelClipReg->AF = 0;

			large_indel_flag = true;
			mate_clip_reg.large_indel_flag = true;
			mate_clip_reg.largeIndelClipReg = dupVarReg(largeIndelClipReg);
			mate_clip_reg.supp_num_largeIndel = supp_num_largeIndel;
		}
	}
}

// compute the clip region based on single vector
reg_t* clipReg::computeClipRegSingleVec(vector<clipPos_t*> &clipPosVector){
	reg_t *reg = NULL;
	int64_t minValue, maxValue;
	int32_t minIdx, maxIdx;

	minIdx = -1, minValue = INT_MAX;
	maxIdx = -1, maxValue = 0;
	for(size_t i=0; i<clipPosVector.size(); i++){
		if(minValue>clipPosVector.at(i)->clipRefPos){
			minValue = clipPosVector.at(i)->clipRefPos;
			minIdx = i;
		}
		if(maxValue<clipPosVector.at(i)->clipRefPos){
			maxValue = clipPosVector.at(i)->clipRefPos;
			maxIdx = i;
		}
	}

	if(minIdx!=-1 and maxIdx!=-1){
		reg = new reg_t();
		reg->chrname = clipPosVector.at(minIdx)->chrname;
		reg->startRefPos = clipPosVector.at(minIdx)->clipRefPos;
		reg->endRefPos = clipPosVector.at(maxIdx)->clipRefPos;
		reg->zero_cov_flag = false;
		reg->aln_seg_end_flag = false;
		reg->query_pos_invalid_flag = false;
		reg->large_indel_flag = false;
		reg->gt_type = -1;
		reg->gt_seq = "";
		reg->AF = 0;
	}

	return reg;
}

// remove unmated clipping regions
void clipReg::removeBNDUnmatedClipRegs(){
	int32_t reg_id, clip_end;
	size_t i;
	reg_t *reg;
	vector<int32_t> mate_mate_reg_id_vec;
	string mate_reg_id_str, mate_reg_id_str2;

	if(mate_clip_reg.leftClipRegNum>=1 and mate_clip_reg.rightClipRegNum>=1){
		for(i=0; i<4; i++){
			reg_id = i;
			switch(reg_id){
				case 0: reg = mate_clip_reg.leftClipReg; break;
				case 1: reg = mate_clip_reg.leftClipReg2; break;
				case 2: reg = mate_clip_reg.rightClipReg; break;
				case 3: reg = mate_clip_reg.rightClipReg2; break;
				default: cerr << "line=" << __LINE__ << ", invalid reg_id=" << reg_id << ", error!" << endl; exit(1);
			}

			if(reg){
				if(mate_clip_reg.bnd_mate_reg_strs[reg_id].compare("-")!=0){
					// right end
					clip_end = RIGHT_END;
					mate_mate_reg_id_vec = getMateMateRegID(reg_id, clip_end, mate_clip_reg.bnd_mate_reg_strs);
					if(mate_mate_reg_id_vec.size()>0 and (mate_mate_reg_id_vec.at(0)!=reg_id or mate_mate_reg_id_vec.at(1)!=clip_end)){ // incorrectly mated
						//cout << "\tMate incorrectly: reg_id=" << reg_id << ", clip_end=" << clip_end << ", mate_mate_reg_id=" << mate_mate_reg_id_vec.at(0) << ", mate_mate_clip_end=" << mate_mate_reg_id_vec.at(1) << endl;
					}

					// left end
					clip_end = LEFT_END;
					mate_mate_reg_id_vec = getMateMateRegID(reg_id, clip_end, mate_clip_reg.bnd_mate_reg_strs);
					if(mate_mate_reg_id_vec.size()>0 and (mate_mate_reg_id_vec.at(0)!=reg_id or mate_mate_reg_id_vec.at(1)!=clip_end)){ // incorrectly mated
						//cout << "\tMate incorrectly: reg_id=" << reg_id << ", clip_end=" << clip_end << ", mate_mate_reg_id=" << mate_mate_reg_id_vec.at(0) << ", mate_mate_clip_end=" << mate_mate_reg_id_vec.at(1) << endl;
					}
				}else{ // no bnd information, then delete the region and its vector
					//cout << "\t===== There are no mate region: " << reg->chrname << ":" << reg->startRefPos << "-" << reg->endRefPos << endl;

					switch(reg_id){
						case 0:
							destroyClipPosVec(leftClipPosVector);
							mate_clip_reg.leftMeanClipPos = mate_clip_reg.leftClipPosNum = 0;
							mate_clip_reg.leftClipRegNum --;
							delete mate_clip_reg.leftClipReg;
							mate_clip_reg.leftClipReg = NULL;
							break;
						case 1:
							destroyClipPosVec(leftClipPosVector2);
							mate_clip_reg.leftMeanClipPos2 = mate_clip_reg.leftClipPosNum2 = 0;
							mate_clip_reg.leftClipRegNum --;
							delete mate_clip_reg.leftClipReg2;
							mate_clip_reg.leftClipReg2 = NULL;
							break;
						case 2:
							destroyClipPosVec(rightClipPosVector);
							mate_clip_reg.rightMeanClipPos = mate_clip_reg.rightClipPosNum = 0;
							mate_clip_reg.rightClipRegNum --;
							delete mate_clip_reg.rightClipReg;
							mate_clip_reg.rightClipReg = NULL;
							break;
						case 3:
							destroyClipPosVec(rightClipPosVector2);
							mate_clip_reg.rightMeanClipPos2 = mate_clip_reg.rightClipPosNum2 = 0;
							mate_clip_reg.rightClipRegNum --;
							delete mate_clip_reg.rightClipReg2;
							mate_clip_reg.rightClipReg2 = NULL;
							break;
						default: cerr << "line=" << __LINE__ << ", invalid reg_id=" << reg_id << ", error!" << endl; exit(1);
					}
				}
			}
		}
	}
}

// get mate of mate region id
vector<int32_t> clipReg::getMateMateRegID(int32_t reg_idx, int32_t end_flag, string bnd_mate_reg_strs[4]){
	vector<int32_t> mate_mate_reg_id_vec;
	int32_t mate_reg_id, mate_reg_id2, mate_end_flag, mate_end_flag2;
	char ch_orient;
	vector<string> str_vec, str_vec2;
	string mate_reg_id_str, mate_reg_id_str2;

	if(bnd_mate_reg_strs[reg_idx].compare("-")!=0){
		// right end
		mate_reg_id_str = getBNDMateRegStr(reg_idx, end_flag, bnd_mate_reg_strs);
		if(mate_reg_id_str.compare("-")!=0){
			mate_reg_id = stoi(mate_reg_id_str.substr(0, 1));
			ch_orient = mate_reg_id_str.at(1);
			if(ch_orient=='+'){ // same orient
				if(end_flag==RIGHT_END)
					mate_end_flag = LEFT_END;
				else
					mate_end_flag = RIGHT_END;
			}else{ // different orient
				if(end_flag==RIGHT_END)
					mate_end_flag = RIGHT_END;
				else
					mate_end_flag = LEFT_END;
			}
			mate_reg_id_str2 = getBNDMateRegStr(mate_reg_id, mate_end_flag, bnd_mate_reg_strs);
			if(mate_reg_id_str2.compare("-")!=0){
				mate_reg_id2 = stoi(mate_reg_id_str2.substr(0, 1));
				ch_orient = mate_reg_id_str2.at(1);
				if(ch_orient=='+'){ // same orient
					if(mate_end_flag==RIGHT_END)
						mate_end_flag2 = LEFT_END;
					else
						mate_end_flag2 = RIGHT_END;
				}else{ // different orient
					if(mate_end_flag==RIGHT_END)
						mate_end_flag2 = RIGHT_END;
					else
						mate_end_flag2 = LEFT_END;
				}

				mate_mate_reg_id_vec.push_back(mate_reg_id2);
				mate_mate_reg_id_vec.push_back(mate_end_flag2);
			}
		}
	}

	return mate_mate_reg_id_vec;
}

// get mate region id string
string clipReg::getBNDMateRegStr(int32_t reg_idx, int32_t end_flag, string bnd_mate_reg_strs[4]){
	vector<string> str_vec, str_vec2;
	string bnd_reg_str, mate_reg_id_str = "-";

	if(reg_idx<0 or reg_idx>=4){
		if(bnd_mate_reg_strs[reg_idx].compare("-")!=0){
			bnd_reg_str = "";
			str_vec = split(bnd_mate_reg_strs[reg_idx], ",");
			if(end_flag==RIGHT_END) // right end
				bnd_reg_str = str_vec.at(0);
			else // left end
				bnd_reg_str = str_vec.at(1);

			if(bnd_reg_str.compare("-")!=0){ // 2+|.|.|.
				str_vec2 = split(bnd_reg_str, "|");  // 2+
				mate_reg_id_str = str_vec2.at(0);
			}
		}
	}
	return mate_reg_id_str;
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
	int64_t dist;
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
size_t clipReg::computeMeanClipPos(vector<clipPos_t*> &clipPosVector){
	size_t i, total, mean_pos;
	total = 0;
	for(i=0; i<clipPosVector.size(); i++) total += clipPosVector.at(i)->clipRefPos;
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

// compute variant type for clipping region
void clipReg::computeVarTypeClipReg(mateClipReg_t &mate_clip_reg){
	size_t i, j, var_type, var_type_tmp, dup_type_num, inv_type_num, tra_type_num, dup_num_tmp;
	vector<size_t> dup_num_vec;
	clipAlnData_t *clip_aln_seg, *mate_clip_aln_seg;
	vector<clipAlnData_t*> query_aln_segs;
	string queryname;
	vector<int32_t> adjClipAlnSegInfo;
	clipPos_t clip_pos_item, mate_clip_pos_item, *clip_pos, *clip_pos_mate, *clip_pos_tmp;
	int32_t arr_idx, mate_arr_idx, clip_end_flag, mate_clip_end_flag;
	int32_t mate_arr_idx_same_chr, mate_clip_end_flag_same_chr, min_ref_dist_same_chr;
	bool same_chr_flag, same_orient_flag, self_overlap_flag, same_aln_reg_flag; // four features
	bool valid_query_end_flag, valid_query_flag;
	vector<clipPos_t*> query_pos_vec;
	int32_t ins_num, del_num, dist, orient_factor;
	int64_t sum_ins, sum_del, mean_size_ins, mean_size_del;
	vector<struct querySeqInfoNode*> query_seq_info_all;

	resetClipCheckFlag(clipAlnDataVector); // reset clipping check flag

	if(mate_clip_reg.large_indel_flag==false){
		dup_type_num = inv_type_num = tra_type_num = 0;
		var_type = VAR_UNC;
		if(mate_clip_reg.reg_mated_flag){ // clippings

			for(i=0; i<clipAlnDataVector.size(); i++){
				clip_aln_seg = clipAlnDataVector.at(i);
				if(clip_aln_seg->query_checked_flag==false){
	//				cout << i << "\t" << clip_aln_seg->chrname << "\t" << clip_aln_seg->startRefPos << "\t" << clip_aln_seg->endRefPos << endl;
					queryname = clipAlnDataVector.at(i)->queryname;
					//query_aln_segs = getQueryClipAlnSegs(queryname, clipAlnDataVector);  // get query clip align segments
					query_aln_segs = getQueryClipAlnSegsAll(queryname, clipAlnDataVector);  // get query clip align segments
					if(query_aln_segs.size()==0) continue;
					valid_query_flag = Filteredbychrname(query_aln_segs);
					if(valid_query_flag==false){
						for(j=0; j<query_aln_segs.size(); j++) query_aln_segs.at(j)->query_checked_flag = true;
						continue;
					}
					//sort according to queryPos in reference position
					sortQueryAlnSegs(query_aln_segs);
					FilteredbyAlignmentSegment(query_aln_segs);

					same_chr_flag = isSameChrome(query_aln_segs);
					same_orient_flag = isSameOrient(query_aln_segs);
					self_overlap_flag = isQuerySelfOverlap(query_aln_segs, maxVarRegSize);
					same_aln_reg_flag = isSameAlnReg(query_aln_segs);

					// predict the variant type
					var_type_tmp = VAR_UNC;
					if(same_chr_flag and same_orient_flag and self_overlap_flag and same_aln_reg_flag){ // DUP
						var_type_tmp = VAR_DUP;
						dup_type_num ++;
					}else if(same_chr_flag and same_orient_flag==false and self_overlap_flag==false){// and same_aln_reg_flag){ // INV
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
							if(clip_aln_seg->left_clip_checked_flag==false and clip_aln_seg->leftClipSize>=minClipEndSize and clip_aln_seg->startRefPos>=startRefPos and clip_aln_seg->startRefPos<=endRefPos){ // left end
								valid_query_end_flag = true;
								clip_end_flag = LEFT_END;
							}else if(clip_aln_seg->right_clip_checked_flag==false and clip_aln_seg->rightClipSize>=minClipEndSize and clip_aln_seg->endRefPos>=startRefPos and clip_aln_seg->endRefPos<=endRefPos){ // right end
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
								adjClipAlnSegInfo = getAdjacentClipAlnSeg(j, clip_end_flag, query_aln_segs, minClipEndSize, MAX_VAR_REG_SIZE);
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
	//			for(j=0; j<query_aln_segs.size(); j++) query_aln_segs.at(j)->query_checked_flag = true;
			}

			// extract var_type
			var_type = extractVarType(mate_clip_reg, dup_type_num, inv_type_num, tra_type_num, int32_t(paras->minReadsNumSupportSV*READS_NUM_SUPPORT_FACTOR));
		}

		mate_clip_reg.sv_type = var_type;
		if(var_type==VAR_DUP) mate_clip_reg.dup_num = computeDupNumClipReg(dup_num_vec); //compute the dup_num

	}else{ // large indels
		sum_ins = sum_del = 0;
		ins_num = del_num = 0;
		for(i=0; i<largeIndelClipPosVector.size(); i++){
			clip_pos = largeIndelClipPosVector.at(i);
			if(clip_pos->clip_aln->query_checked_flag==false){
				queryname = clip_pos->clip_aln->queryname;
				query_pos_vec = getClipPosItemsByQueryname(queryname, largeIndelClipPosVector);  // get query clip pos
				query_aln_segs = getQueryClipAlnSegsAll(queryname, clipAlnDataVector);  // get query clip align segments

				adjClipAlnSegInfo = getAdjacentClipAlnSeg(clip_pos->clip_aln, clip_pos->clip_end, query_aln_segs, minClipEndSize, MAX_VAR_REG_SIZE);
				mate_arr_idx = adjClipAlnSegInfo.at(0);
				mate_clip_end_flag = adjClipAlnSegInfo.at(1);
				mate_arr_idx_same_chr = adjClipAlnSegInfo.at(4);
				mate_clip_end_flag_same_chr = adjClipAlnSegInfo.at(5);
				//min_dist_same_chr = adjClipAlnSegInfo.at(6);
				min_ref_dist_same_chr = adjClipAlnSegInfo.at(7);
				if(mate_arr_idx!=-1){ // mated
					if(mate_arr_idx_same_chr!=-1){
						mate_clip_aln_seg = query_aln_segs.at(mate_arr_idx_same_chr);
						self_overlap_flag = isSegSelfOverlap(clip_pos->clip_aln, mate_clip_aln_seg, maxVarRegSize);
						// prefer the segments on the same chromosome
						if(mate_arr_idx!=mate_arr_idx_same_chr and (abs(min_ref_dist_same_chr)<=MAX_REF_DIST_SAME_CHR or self_overlap_flag)){ // different chromosomes with near reference distance
							mate_arr_idx = mate_arr_idx_same_chr;
							mate_clip_end_flag = mate_clip_end_flag_same_chr;
						}
					}

					mate_clip_aln_seg = query_aln_segs.at(mate_arr_idx);
					clip_pos_mate = NULL;
					for(j=0; j<query_pos_vec.size(); j++){
						clip_pos_tmp = query_pos_vec.at(j);
						if(clip_pos_tmp->clip_aln==mate_clip_aln_seg and clip_pos_tmp->clip_end==mate_clip_end_flag){
							clip_pos_mate = clip_pos_tmp;
							break;
						}
					}
					if(clip_pos_mate){
						if(clip_pos->aln_orient==clip_pos_mate->aln_orient){ // same orient
							if(clip_pos->aln_orient==ALN_PLUS_ORIENT and clip_pos_mate->aln_orient==ALN_PLUS_ORIENT) orient_factor = 1;
							else if(clip_pos->aln_orient==ALN_MINUS_ORIENT and clip_pos_mate->aln_orient==ALN_MINUS_ORIENT) orient_factor = -1;

							if(clip_pos->clip_end==RIGHT_END and clip_pos_mate->clip_end==LEFT_END)
								dist = (clip_pos_mate->clipQueryPos - clip_pos->clipQueryPos) * orient_factor - (clip_pos_mate->clipRefPos - clip_pos->clipRefPos);
							else
								dist = (clip_pos->clipQueryPos - clip_pos_mate->clipQueryPos) * orient_factor - (clip_pos->clipRefPos - clip_pos_mate->clipRefPos);

							if(dist>0) { ins_num ++; sum_ins += dist;}
							else { del_num ++; sum_del += dist; }
						}
					}
				}

				for(j=0; j<query_pos_vec.size(); j++) query_pos_vec.at(j)->clip_aln->query_checked_flag = true;
			}
		}

		if(ins_num>0) mean_size_ins = sum_ins / ins_num;
		else mean_size_ins = 0;
		if(del_num>0) mean_size_del = sum_del / del_num;
		else mean_size_del = 0;

		if(ins_num>=int32_t(paras->minReadsNumSupportSV*READS_NUM_SUPPORT_FACTOR) or del_num>=int32_t(paras->minReadsNumSupportSV*READS_NUM_SUPPORT_FACTOR)){
			if(ins_num>=del_num) { mate_clip_reg.sv_type = largeIndelClipReg->var_type = VAR_INS; mate_clip_reg.largeIndelClipReg->sv_len = largeIndelClipReg->sv_len = mean_size_ins; mate_clip_reg.supp_num_largeIndel = supp_num_largeIndel = ins_num; }
			else { mate_clip_reg.sv_type = largeIndelClipReg->var_type = VAR_DEL; mate_clip_reg.largeIndelClipReg->sv_len = largeIndelClipReg->sv_len = mean_size_del; mate_clip_reg.supp_num_largeIndel = supp_num_largeIndel = del_num; }

			// compute depth for large indel region
			mate_clip_reg.depth_largeIndel = computeCovNumReg(largeIndelClipReg->chrname, largeIndelClipReg->startRefPos, largeIndelClipReg->endRefPos, fai, inBamFile);
			if(mate_clip_reg.depth_largeIndel<mate_clip_reg.supp_num_largeIndel) mate_clip_reg.depth_largeIndel = mate_clip_reg.supp_num_largeIndel; // tolerate DEL region
			mate_clip_reg.supp_num_valid_flag = true;

			for(i=0; i<4; i++) mate_clip_reg.bnd_mate_reg_strs[i] = "-";

		}else{ // insufficient supporting reads
			delete largeIndelClipReg;
			largeIndelClipReg = NULL;
			large_indel_flag = false;

			delete mate_clip_reg.largeIndelClipReg;
			mate_clip_reg.largeIndelClipReg = NULL;
			mate_clip_reg.large_indel_flag = false;
			mate_clip_reg.valid_flag = false;
			mate_clip_reg.supp_num_valid_flag = false;

			if(!largeIndelClipPosVector.empty()) destroyClipPosVec(largeIndelClipPosVector);
		}
	}
	resetClipCheckFlag(clipAlnDataVector); // reset clipping check flag
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
bool clipReg::Filteredbychrname(vector<clipAlnData_t*> &query_aln_segs){
	bool flag =true;
	vector<string> chrnameV;
	size_t i, j;

	if(query_aln_segs.size()>=3){
		for(i=0; i<query_aln_segs.size(); i++){
//			if(chrnameV.size()==0){
//				chrnameV.push_back(query_aln_segs.at(i)->chrname);
//			}
//			else{
				for(j=0; j<chrnameV.size(); j++){
					if(query_aln_segs.at(i)->chrname.compare(chrnameV.at(j))==0)
						break;
				}
				if(j == chrnameV.size())
					chrnameV.push_back(query_aln_segs.at(i)->chrname);
//			}
		}
	}
	if(chrnameV.size()>=3)
		flag=false;
	return flag;
}

void clipReg::FilteredbyAlignmentSegment(vector<clipAlnData_t*> &query_aln_segs){
	size_t j;
	int32_t samechr_num;

	if(query_aln_segs.size()>3){
		samechr_num = 0;
		for(j=0; j<query_aln_segs.size(); j++){
			if(query_aln_segs.at(j)->chrname.compare(chrname)==0)
				samechr_num++;
		}
		bool orient_flag = ischeckSameOrient(query_aln_segs);
		if(samechr_num>=3 and orient_flag==false){		//Filter unreliable read
			for(j=0; j<query_aln_segs.size(); j++){
				if(query_aln_segs.at(j)->chrname.compare(chrname)!=0)
					query_aln_segs.erase(query_aln_segs.begin()+j);
			}
			query_aln_segs.shrink_to_fit();
		}
	}
}

bool clipReg::ischeckSameOrient(vector<clipAlnData_t*> &query_aln_segs){
	bool flag = true;
	int32_t Differentdirections_num = 0;
	for(size_t i=1; i<query_aln_segs.size(); i++){
		if(query_aln_segs.at(i)->aln_orient!=query_aln_segs.at(i-1)->aln_orient){
			Differentdirections_num++;
		}
	}
	if(Differentdirections_num >= 2) flag = false;
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
					adjacent_flag = isAdjacentClipAlnSeg(query_aln_segs.at(i-1), query_aln_segs.at(i), maxVarRegSize);
					if(adjacent_flag==false){
						same_aln_reg_flag = false;
						break;
					}
				}
//			}else{ // different orient
//				for(i=1; i<query_aln_segs.size(); i++){
//					adjacent_flag = isAdjacentClipAlnSeg(query_aln_segs.at(i-1), query_aln_segs.at(i), paras->maxVarRegSize);
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

// sort query align segments according to queryPos
void clipReg::sortQueryAlnSegs(vector<clipAlnData_t*> &query_aln_segs){
	size_t i, j, minIdx, minPos, query_pos1, query_pos2;
	clipAlnData_t *clip_aln1, *clip_aln2,  *tmp;

	for(i=0; i<query_aln_segs.size(); i++){
		clip_aln1 = query_aln_segs.at(i);
		if(clip_aln1->aln_orient==ALN_PLUS_ORIENT)
			query_pos1 = clip_aln1->startRefPos;
//			query_pos1 = clip_aln1->startQueryPos;
		else
			query_pos1 = clip_aln1->endRefPos;
//			query_pos1 = clip_aln1->querylen - clip_aln1->startQueryPos + 1;

		minIdx = i;
		minPos = query_pos1;
		for(j=i+1; j<query_aln_segs.size(); j++){
			clip_aln2 = query_aln_segs.at(j);
			if(clip_aln2->aln_orient==ALN_PLUS_ORIENT)
				query_pos2 = clip_aln2->startRefPos;
//				query_pos2 = clip_aln2->startQueryPos;
			else
				query_pos2 = clip_aln2->endRefPos;
//				query_pos2 = clip_aln2->querylen - clip_aln2->startQueryPos + 1;

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

size_t clipReg::extractVarType(mateClipReg_t &mate_clip_reg, size_t dup_type_num, size_t inv_type_num, size_t tra_type_num, size_t min_reads_thres){
	size_t var_type, maxValue, maxIdx;
	string chrname1, chrname2;

	if(mate_clip_reg.leftClipReg) chrname1 = mate_clip_reg.leftClipReg->chrname;
	else if(mate_clip_reg.leftClipReg2) chrname1 = mate_clip_reg.leftClipReg2->chrname;
	if(mate_clip_reg.rightClipReg) chrname2 = mate_clip_reg.rightClipReg->chrname;
	else if(mate_clip_reg.rightClipReg2) chrname2 = mate_clip_reg.rightClipReg2->chrname;

	if(chrname1.compare(chrname2)!=0){
		maxValue = tra_type_num;
		maxIdx = 2;
	}else{
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
void clipReg::printClipVecs(string head_info){
	clipPos_t clip_pos;
	double ratio1, ratio2, total;

	cout << "\n############## " << head_info << " ##############" << endl;

	cout << "left_mate_clip_pos_vec1: " << leftClipPosVector.size() << endl;
	printSingleClipVec(leftClipPosVector);

	cout << "left_mate_clip_pos_vec2: " << leftClipPosVector2.size() << endl;
	printSingleClipVec(leftClipPosVector2);

	cout << "right_mate_clip_pos_vec1: " << rightClipPosVector.size() << endl;
	printSingleClipVec(rightClipPosVector);

	cout << "right_mate_clip_pos_vec2: " << rightClipPosVector2.size() << endl;
	printSingleClipVec(rightClipPosVector2);

	ratio1 = -1, ratio2 = -1, total = leftClipPosVector.size() + leftClipPosVector2.size() + rightClipPosVector.size() + rightClipPosVector2.size();
	if(total>0) ratio1 = (double)(leftClipPosVector.size() + leftClipPosVector2.size()) / total;
	if(total>0) ratio2 = (double)(rightClipPosVector.size() + rightClipPosVector2.size()) / total;
	cout << "ratio1=" << ratio1 << "; " << "ratio2=" << ratio2 << endl;
}

void clipReg::printSingleClipVec(vector<clipPos_t*> &clip_pos_vec){
	clipPos_t *clip_pos;
	clipAlnData_t *left_clip_aln, *right_clip_aln;
	string line, line2, clip_end_str, flag_str, flag_str2, qname;
	int32_t left_clip_size, right_clip_size;

	for(size_t i=0; i<clip_pos_vec.size(); i++){
		clip_pos = clip_pos_vec.at(i);
		qname = clip_pos->clip_aln->queryname;
		if(clip_pos->clip_end==LEFT_END) clip_end_str = "Left";
		else clip_end_str = "Right";
		if(clip_pos->same_orient_flag) flag_str = "1";
		else flag_str = "0";
		if(clip_pos->self_overlap_flag) flag_str2 = "1";
		else flag_str2 = "0";
		line = "\t[" + to_string(i) + "]: " + clip_pos->chrname + ":" + to_string(clip_pos->clipRefPos) + ", qname=" + qname + ", clipping end=" + clip_end_str + ", same_orient=" + flag_str + ", self overlap flag=" + flag_str2;
		if(clip_pos->clip_aln) {
			left_clip_size = clip_pos->clip_aln->leftClipSize;
			right_clip_size = clip_pos->clip_aln->rightClipSize;
			if(clip_pos->clip_aln->aln_orient==ALN_PLUS_ORIENT) flag_str = "+";
			else flag_str = "-";
			line += ", left clipping size=" + to_string(left_clip_size) + ", right clipping size=" + to_string(right_clip_size) + ", strand=" + flag_str;
		}
		cout << line << endl;

		if(clip_pos->clip_aln) {
			// print left segment and right segment
			left_clip_aln = clip_pos->clip_aln->left_aln;
			if(left_clip_aln){
				line2 = "\t\t[left]: " + left_clip_aln->chrname + ":" + to_string(left_clip_aln->startRefPos) + "-" + to_string(left_clip_aln->endRefPos);
				left_clip_size = left_clip_aln->leftClipSize;
				right_clip_size = left_clip_aln->rightClipSize;
				if(left_clip_aln->aln_orient==ALN_PLUS_ORIENT) flag_str = "+";
				else flag_str = "-";
				line2 += ", left clipping size=" + to_string(left_clip_size) + ", right clipping size=" + to_string(right_clip_size) + ", strand=" + flag_str;
				cout << line2 << endl;
			}

			right_clip_aln = clip_pos->clip_aln->right_aln;
			if(right_clip_aln){
				line2 = "\t\t[right]: " + right_clip_aln->chrname + ":" + to_string(right_clip_aln->startRefPos) + "-" + to_string(right_clip_aln->endRefPos);
				left_clip_size = right_clip_aln->leftClipSize;
				right_clip_size = right_clip_aln->rightClipSize;
				if(right_clip_aln->aln_orient==ALN_PLUS_ORIENT) flag_str = "+";
				else flag_str = "-";
				line2 += ", left clipping size=" + to_string(left_clip_size) + ", right clipping size=" + to_string(right_clip_size) + ", strand=" + flag_str;
				cout << line2 << endl;
			}
		}
	}
}

// print clip regions
void clipReg::printResultClipRegs(){
	reg_t *reg;
	string sv_type_str;

	switch(mate_clip_reg.sv_type){
		case VAR_DUP: // DUP
			sv_type_str = VAR_DUP_STR;
			break;
		case VAR_INV: // INV
			sv_type_str = VAR_INV_STR;
			break;
		case VAR_INS: // large INS or DEL
			sv_type_str = VAR_INS_STR;
			break;
		case VAR_DEL:
			sv_type_str = VAR_DEL_STR;
			break;
		case VAR_TRA:
			sv_type_str = VAR_TRA_STR;
			break;
		default: sv_type_str = VAR_UNC_STR;
	}

	cout << "The retult clip regions: " << sv_type_str << endl;

	if(mate_clip_reg.large_indel_flag==false){ // clipping
		reg = mate_clip_reg.leftClipReg;
		if(reg){
			cout << "The left clip region 1:" << endl;
			cout << "\t" << reg->chrname << ":" << reg->startRefPos << "-" << reg->endRefPos << endl;
		}
		reg = mate_clip_reg.leftClipReg2;
		if(reg){
			cout << "The left clip region 2:" << endl;
			cout << "\t" << reg->chrname << ":" << reg->startRefPos << "-" << reg->endRefPos << endl;
		}
		reg = mate_clip_reg.rightClipReg;
		if(reg){
			cout << "The right clip region 1:" << endl;
			cout << "\t" << reg->chrname << ":" << reg->startRefPos << "-" << reg->endRefPos << endl;
		}
		reg = mate_clip_reg.rightClipReg2;
		if(reg){
			cout << "The right clip region 2:" << endl;
			cout << "\t" << reg->chrname << ":" << reg->startRefPos << "-" << reg->endRefPos << endl;
		}
	}else{ // large indel
		if(mate_clip_reg.large_indel_flag){
			cout << "The large indel region:" << endl;
			cout << "\t" << mate_clip_reg.largeIndelClipReg->chrname << ":" << mate_clip_reg.largeIndelClipReg->startRefPos << "-" << mate_clip_reg.largeIndelClipReg->endRefPos << ", sv_len=" << mate_clip_reg.largeIndelClipReg->sv_len << endl;
		}
	}
}
