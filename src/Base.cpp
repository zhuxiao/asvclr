#include "Base.h"
#include "util.h"

// Constructor
Base::Base(){
	init();
}

// Destructor
Base::~Base(){
	destroyBase();
}

// initialization
void Base::init(){
	int i;
	for(i=0; i<6; i++) coverage.num_bases[i] = 0;
	coverage.idx_RefBase = -1;
	coverage.num_max = -1;
	coverage.idx_max = -1;

	num_shortIns = num_shortdel = num_shortClip = 0;
}

// destroy base, including insVector, delVector and clipVector
void Base::destroyBase(){
	if(!insVector.empty()) destroyInsVector();
	if(!delVector.empty()) destroyDelVector();
	if(!clipVector.empty()) destroyClipVector();
}

// destroy the insVector
void Base::destroyInsVector(){
	vector<insEvent_t*>::iterator ins;
	for(ins=insVector.begin(); ins!=insVector.end(); ins++)
		delete (*ins);
	vector<insEvent_t*>().swap(insVector);
}

// destroy the delVector
void Base::destroyDelVector(){
	vector<delEvent_t*>::iterator del;
	for(del=delVector.begin(); del!=delVector.end(); del++)
		delete (*del);
	vector<delEvent_t*>().swap(delVector);
}

// destroy the clipVector
void Base::destroyClipVector(){
	vector<clipEvent_t*>::iterator clip;
	for(clip=clipVector.begin(); clip!=clipVector.end(); clip++)
		delete (*clip);
	vector<clipEvent_t*>().swap(clipVector);
}

// add an insertion event
void Base::addInsEvent(insEvent_t* insE){
	insVector.push_back(insE);
}

// add an deletion event
void Base::addDelEvent(delEvent_t* delE){
	delVector.push_back(delE);
}

// add a clip event
void Base::addClipEvent(clipEvent_t* clipE){
	clipVector.push_back(clipE);
}

// compute the maximal base and its count
void Base::updateCovInfo(){
	uint16_t i, maxIdx, maxNum, *cov_nums = coverage.num_bases;
	coverage.num_bases[5] = cov_nums[0] + cov_nums[1] + cov_nums[2] + cov_nums[3] + cov_nums[4];  // sum
	// get max and count
	maxIdx = 0; maxNum = cov_nums[0];
	for(i=1; i<5; i++) if(maxNum<cov_nums[i]){ maxIdx = i; maxNum = cov_nums[i]; }
	coverage.idx_max = maxIdx;
	coverage.num_max = maxNum;
}

// determine whether the base is a disagreements
bool Base::isDisagreeBase(){
	bool flag = false;
	int32_t baseNum_A, baseNum_C, baseNum_G, baseNum_T, total_cov;
	double total_tmp;

	if((double)coverage.num_max/coverage.num_bases[5]<=DISAGREE_THRES)
		flag = true;
	else if(coverage.idx_max!=coverage.idx_RefBase){
		if(coverage.idx_RefBase==5){
			baseNum_A = coverage.num_bases[0];
			baseNum_C = coverage.num_bases[1];
			baseNum_G = coverage.num_bases[2];
			baseNum_T = coverage.num_bases[3];
			total_cov = coverage.num_bases[5];
			switch(coverage.refBase){
				case 'M':
				case 'm':
					total_tmp = baseNum_A + baseNum_C;
					if((coverage.idx_max!=0 and coverage.idx_max!=1) or total_tmp/total_cov<=DISAGREE_THRES) flag = true;
					break;
				case 'R':
				case 'r':
					total_tmp = baseNum_A + baseNum_G;
					if((coverage.idx_max!=0 and coverage.idx_max!=2) or total_tmp/total_cov<=DISAGREE_THRES) flag = true;
					break;
				case 'S':
				case 's':
					total_tmp = baseNum_C + baseNum_G;
					if((coverage.idx_max!=1 and coverage.idx_max!=2) or total_tmp/total_cov<=DISAGREE_THRES) flag = true;
					break;
				case 'V':
				case 'v':
					total_tmp = baseNum_A + baseNum_C + baseNum_G;
					if((coverage.idx_max!=0 and coverage.idx_max!=1 and coverage.idx_max!=2) or total_tmp/total_cov<=DISAGREE_THRES) flag = true;
					break;
				case 'W':
				case 'w':
					total_tmp = baseNum_A + baseNum_T;
					if((coverage.idx_max!=0 and coverage.idx_max!=3) or total_tmp/total_cov<=DISAGREE_THRES) flag = true;
					break;
				case 'Y':
				case 'y':
					total_tmp = baseNum_C + baseNum_T;
					if((coverage.idx_max!=1 and coverage.idx_max!=3) or total_tmp/total_cov<=DISAGREE_THRES) flag = true;
					break;
				case 'H':
				case 'h':
					total_tmp = baseNum_A + baseNum_C + baseNum_T;
					if((coverage.idx_max!=0 and coverage.idx_max!=1 and coverage.idx_max!=3) or total_tmp/total_cov<=DISAGREE_THRES) flag = true;
					break;
				case 'K':
				case 'k':
					total_tmp = baseNum_G + baseNum_T;
					if((coverage.idx_max!=2 and coverage.idx_max!=3) or total_tmp/total_cov<=DISAGREE_THRES) flag = true;
					break;
				case 'D':
				case 'd':
					total_tmp = baseNum_A + baseNum_G + baseNum_T;
					if((coverage.idx_max!=0 and coverage.idx_max!=2 and coverage.idx_max!=3) or total_tmp/total_cov<=DISAGREE_THRES) flag = true;
					break;
				case 'B':
				case 'b':
					total_tmp = baseNum_C + baseNum_G + baseNum_T;
					if((coverage.idx_max!=1 and coverage.idx_max!=2 and coverage.idx_max!=3) or total_tmp/total_cov<=DISAGREE_THRES) flag = true;
					break;
				default: cerr << __func__ << ": unknown base: " << coverage.refBase << endl; exit(1);
			}
		}else
			flag = true;
	}
	return flag;
}

// determine whether the base is a zero coverage base
bool Base::isZeroCovBase(){
	bool flag = false;
	if(coverage.num_bases[5]<3) flag = true; // zero coverage
	return flag;
}

// determine whether the base is a high indel base
bool Base::isHighIndelBase(float threshold){
	bool flag = false;
	int32_t indelNum;

	indelNum = insVector.size() + delVector.size();
	if(coverage.num_bases[5]>0 and (double)indelNum/coverage.num_bases[5]>=threshold)
		flag = true;
	return flag;
}

bool Base::isMatchToRef(){
	bool flag = false, matchFlag;
	char ch_maxBase, ref_base;

	if(coverage.num_bases[5]>0){
		if(coverage.idx_max==coverage.idx_RefBase){
			flag = true;
		}else{
			switch(coverage.idx_max){
				case 0: ch_maxBase = 'A'; break;
				case 1: ch_maxBase = 'C'; break;
				case 2: ch_maxBase = 'G'; break;
				case 3: ch_maxBase = 'T'; break;
				case 4: ch_maxBase = 'N'; break;
				default: cerr << __func__ << ": unknown base idx_max: " << coverage.idx_max << endl; exit(1);
			}

			ref_base = coverage.refBase;
			matchFlag = isBaseMatch(ch_maxBase, ref_base);
			if(matchFlag)
				flag = true;
		}
	}
	return flag;
}

size_t Base::getLargeIndelNum(size_t thres){
	size_t i, large_indel_num;

	large_indel_num = 0;
	for(i=0; i<insVector.size(); i++){
		if(insVector.at(i)->seq.size()>=thres)
			large_indel_num ++;
	}
	for(i=0; i<delVector.size(); i++){
		if(delVector.at(i)->seq.size()>=thres)
			large_indel_num ++;
	}
	return large_indel_num;
}

size_t Base::getTotalIndelNum(){
	return insVector.size() + delVector.size() + num_shortIns + num_shortdel;
}

size_t Base::getTotalClipNum(){
	return clipVector.size() + num_shortClip;
}

size_t Base::getTotalCovNum(){
	return coverage.num_bases[5];
}
