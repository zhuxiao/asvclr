#include "Base.h"

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
	if(coverage.idx_max!=coverage.idx_RefBase or (double)coverage.num_max/coverage.num_bases[5]<=DISAGREE_THRES)
		flag = true;
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
