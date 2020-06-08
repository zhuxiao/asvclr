#include "misAlnReg.h"

misAlnReg::misAlnReg(size_t startPos, size_t endPos, size_t chrlen, Base *misAlnRegBaseArr) {
	this->startPos = startPos;
	this->endPos = endPos;
	this->chrlen = chrlen;
	this->misAlnRegBaseArr = misAlnRegBaseArr;
	this->disagrNum = 0;
	this->misAlnSubregNum = 0;
	this->subRegNum = 0;
	this->disagrRegRatio = 0;
	this->highClipBaseNum = 0;
	this->zeroCovBaseNum = 0;
	this->misAlnFlag = false;
	computeDisagrSubreg();
	computeHighClipBaseNum();
	computeZeroCovBaseNum();
}

// compute the disagreements of sub regions in single misAln region
void misAlnReg::computeDisagrSubreg(){
	int64_t pos, posIdx, regIdx;

	// divide into sub-regions
	subRegNum = (endPos - startPos) / SUB_MIS_ALN_REG_SIZE + 1;
	size_t disNumArray[subRegNum];
	for(pos=0; pos<subRegNum; pos++) disNumArray[pos] = 0;

	// compute disagreements for each sub-region
	disagrNum = 0;
	for(pos=startPos; pos<endPos; pos++){
		posIdx = pos - startPos;
		if(misAlnRegBaseArr[posIdx].coverage.idx_RefBase!=4){ // A, C, G, T, Mixed, but N
			if(misAlnRegBaseArr[posIdx].isDisagreeBase()){
				disagrNum ++;
				regIdx = (pos - startPos) / SUB_MIS_ALN_REG_SIZE + 1;
				disNumArray[regIdx] ++;
			}
		}
	}

	// compute the misAln sub-regions
	misAlnSubregNum = 0;
	for(regIdx=0; regIdx<subRegNum; regIdx++) if(disNumArray[regIdx]>0 and disNumArray[regIdx]<=SUB_REG_MAX_MIS_ALN_NUM_THRES) misAlnSubregNum ++;

	disagrRegRatio = (float)misAlnSubregNum / subRegNum;
}

// compute the number of bases which having much clipped events
void misAlnReg::computeHighClipBaseNum(){
	int64_t pos, posIdx;
	highClipBaseNum = 0;
	for(pos=startPos; pos<=endPos; pos++){
		posIdx = pos - startPos;
		if(misAlnRegBaseArr[posIdx].clipVector.size()>=CLIP_SUPPORT_READS_NUM_THRES and (double)misAlnRegBaseArr[posIdx].clipVector.size()/misAlnRegBaseArr[posIdx].coverage.num_bases[5]>=HIGH_CLIP_RATIO_THRES) highClipBaseNum ++;
	}
}

// compute the number of bases which having much clipped events
void misAlnReg::computeZeroCovBaseNum(){
	int64_t pos, posIdx;
	if(startPos>=REF_END_SKIP_SIZE and endPos+REF_END_SKIP_SIZE<=chrlen){
		zeroCovBaseNum = 0;
		for(pos=startPos; pos<=endPos; pos++){
			posIdx = pos - startPos;
			if(misAlnRegBaseArr[posIdx].isZeroCovBase()) zeroCovBaseNum ++;
		}
	}
}
