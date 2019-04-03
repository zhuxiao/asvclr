#include <algorithm>
#include "Region.h"
#include "clipReg.h"
#include "util.h"

//Constructor
Region::Region(string& chrname, size_t startRpos, size_t endRPos, size_t chrlen, Base *regBaseArr, size_t regFlag, Paras *paras) {
	this->paras = paras;
	this->chrname = chrname;
	this->startRPos = startRpos;
	this->endRPos = endRPos;
	this->chrlen = chrlen;
	this->regBaseArr = regBaseArr;
	this->regFlag = regFlag;
	fai = fai_load(paras->refFile.c_str());
	if(!fai){
		cerr << __func__ << ": could not load fai index of " << paras->refFile << endl;;
		exit(1);
	}

	switch(regFlag){
		case HEAD_REGION:
			startMidPartPos = startRpos;
			endMidPartPos = startMidPartPos + paras->slideSize - 1;
			if(endMidPartPos>endRPos) endMidPartPos = endRPos;
			break;
		case INNER_REGION:
			startMidPartPos = startRpos + paras->slideSize;
			if(startMidPartPos>endRPos) {
				cerr << __func__ << ", line=" << __LINE__ << ": error!" << endl; exit(1);
			}
			endMidPartPos = startMidPartPos + paras->slideSize - 1;
			if(endMidPartPos>endRPos) {
				cerr << __func__ << ", line=" << __LINE__ << ": error!" << endl; exit(1);
			}
			break;
		case TAIL_REGION:
			startMidPartPos = startRpos + paras->slideSize;
			if(startMidPartPos>endRPos) {
				cerr << __func__ << ", line=" << __LINE__ << ": error!" << endl; exit(1);
			}
			endMidPartPos = startMidPartPos + paras->slideSize - 1;
			if(endMidPartPos>endRPos) endMidPartPos = endRPos;
			break;
	}

	subRegSize = SUB_REG_SIZE;
	meanBlockCov = 0;

	// determine if all the base in the region are 'N' bases in reference
	wholeRefGapFlag = IsWholeRefGap();

	// compute the local region coverage
	localRegCov = computeMeanCovReg(startMidPartPos, endMidPartPos);
	highIndelSubRegNum = 0;

	indelCandFlag = false;
	difType = REG_DIF_NONE;
}

//Destructor
Region::~Region(){
	fai_destroy(fai);
	if(!disagrePosVector.empty()) destroyDisagrePosVector();
	if(!zeroCovPosVector.empty()) destroyZeroCovPosVector();
	if(!abCovRegVector.empty()) destroyAbCovRegVector();
	if(!snvVector.empty()) destroySnvVector();
	if(!indelVector.empty()) destroyIndelVector();
	if(!clipRegVector.empty()) destroyClipRegVector();
}

// set the output directory
void Region::setOutputDir(string& out_dir_assemble_prefix){
	out_dir_assemble = out_dir_assemble_prefix + "/" + chrname;
}

// determine if all the base in the region are 'N' bases in reference
bool Region::IsWholeRefGap(){
	bool flag = true;
	for(size_t i=startMidPartPos; i<=endMidPartPos; i++)
		if(regBaseArr[i-startRPos].coverage.idx_RefBase!=4){ // excluding 'N'
			flag = false;
			break;
		}
	return flag;
}

// set the block mean coverage
void Region::setMeanBlockCov(double meanBlockCov){
	this->meanBlockCov = meanBlockCov;
}

// compute abnormal signatures in region
int Region::computeRegAbSigs(){
	// compute disagreements
	computeDisagreements();

	// compute abnormal coverage, including high/low coverage
	computeAbCovReg();

	// compute the number sub-regions of high indel events in the region
	computeHighIndelEventRegNum();

	return 0;
}

// compute disagreements, only check the middle part sub-region
int Region::computeDisagreements(){
	size_t pos, regIdx;
	for(pos=startMidPartPos; pos<endMidPartPos; pos++){
		regIdx = pos - startRPos;
		if(regBaseArr[regIdx].coverage.idx_RefBase!=4){ // A, C, G, T, but N
			if(regBaseArr[regIdx].isDisagreeBase()) addDisagrePos(pos);
			if(regBaseArr[regIdx].isZeroCovBase()) addZeroCovPos(pos);
		}
	}
	disagrePosVector.shrink_to_fit();
	zeroCovPosVector.shrink_to_fit();
	return 0;
}

// add base position having disagreement to vector
void Region::addDisagrePos(size_t pos){
	disagrePosVector.push_back(pos);
}

// add base position having zero coverage to vector
void Region::addZeroCovPos(size_t pos){
	zeroCovPosVector.push_back(pos);
}

// destroy the disagreements vector
void Region::destroyDisagrePosVector(){
	vector<size_t>().swap(disagrePosVector);
}

// destroy the zero coverage vector
void Region::destroyZeroCovPosVector(){
	vector<size_t>().swap(zeroCovPosVector);
}

// destroy the abnormal coverage vector
void Region::destroyAbCovRegVector(){
	vector<reg_t*>::iterator it;
	for(it=abCovRegVector.begin(); it!=abCovRegVector.end(); it++)
		delete (*it);
	vector<reg_t*>().swap(abCovRegVector);
}

// destroy the SNV vector
void Region::destroySnvVector(){
	vector<size_t>().swap(snvVector);
}

// destroy the indel vector
void Region::destroyIndelVector(){
	vector<reg_t*>().swap(indelVector);
}

// destroy the clip region vector
void Region::destroyClipRegVector(){
	vector<reg_t*>().swap(clipRegVector);
}

// output abnormal signatures in region
void Region::printRegAbSigs(){
	vector<size_t>::iterator it;

	// output SNVs
	cout << "region: " << startMidPartPos << "-" << endMidPartPos << endl;
	cout << "    Number of SNVs: " << getSNVNum() << endl;
	for(it=snvVector.begin(); it!=snvVector.end(); it++)
		cout << "        " << (*it) << endl;

	// output Indels
	cout << "    Number of Indels: " << getIndelNum() << endl;
	vector<reg_t*>::iterator reg;
	for(reg=indelVector.begin(); reg!=indelVector.end(); reg++)
		cout << "        " << (*reg)->startRefPos << "-" << (*reg)->endRefPos << endl;

	// output position of disagreements
	cout << "    Number of disagreements: " << getDisagreeNum() << endl;
	for(it=disagrePosVector.begin(); it!=disagrePosVector.end(); it++)
		cout << "        " << (*it) << endl;

	// output position of zero coverage
	cout << "    Number of zero coverage: " << zeroCovPosVector.size() << endl;
	for(it=zeroCovPosVector.begin(); it!=zeroCovPosVector.end(); it++)
		cout << "        " << (*it) << endl;

	// output sub-regions of abnormal coverage
	cout << "    Number of sub-regions with abnormal coverage: " << abCovRegVector.size() << endl;
	for(reg=abCovRegVector.begin(); reg!=abCovRegVector.end(); reg++)
		cout << "        " << (*reg)->startRefPos << "-" << (*reg)->endRefPos << endl;

	// output the number of high indel sub-regions
	cout << "    Number of high indel sub-regions: " << highIndelSubRegNum << endl;
}

// compute the abnormal coverage region
int Region::computeAbCovReg(){
	size_t pos, startPosSubReg, endPosSubReg;
	double tmp_cov, covRatio2block, covRatio2local;
	reg_t *reg;

	for(pos=startMidPartPos; pos<=endMidPartPos; pos+=subRegSize){
		startPosSubReg = pos;
		endPosSubReg = startPosSubReg + subRegSize - 1;
		if(endPosSubReg>endMidPartPos) endPosSubReg = endMidPartPos;
		tmp_cov = computeMeanCovReg(startPosSubReg, endPosSubReg);
		covRatio2block = tmp_cov / meanBlockCov;
		covRatio2local = tmp_cov / localRegCov;
		if((covRatio2block<0.6 and covRatio2local<0.6) or (covRatio2block>2 and covRatio2local>2)) {  // abnormal coverage condition
			reg = allocateReg(startPosSubReg, endPosSubReg);
			abCovRegVector.push_back(reg);
		}
	}
	abCovRegVector.shrink_to_fit();

	return 0;
}

// allocate reg
reg_t* Region::allocateReg(size_t startPosReg, size_t endPosReg){
	reg_t *reg = new reg_t();
	if(!reg){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot allocate memory" << endl;
		exit(1);
	}
	reg->chrname = chrname;
	reg->startRefPos = startPosReg;
	reg->endRefPos = endPosReg;
	reg->startLocalRefPos = reg->endLocalRefPos = 0;
	reg->startQueryPos = reg->endQueryPos = 0;
	reg->var_type = VAR_UNC;
	reg->sv_len = 0;
	reg->query_id = -1;
	reg->blat_aln_id = -1;
	reg->call_success_status = false;
	reg->short_sv_flag = false;
	return reg;
}

// compute the mean coverage of the sub-region, excluding the gap region
double Region::computeMeanCovReg(size_t startPosReg, size_t endPosReg){
	uint64_t i, totalReadBeseNum = 0, totalRefBaseNum = 0;
	double mean_cov;
	for(i=startPosReg; i<=endPosReg; i++)
		if(regBaseArr[i-startRPos].coverage.idx_RefBase!=4){ // excluding 'N'
			totalReadBeseNum += regBaseArr[i-startRPos].coverage.num_bases[5];
			totalRefBaseNum ++;
		}
	if(totalRefBaseNum) mean_cov = (double)totalReadBeseNum/totalRefBaseNum;
	else mean_cov = 0;
	return mean_cov;
}

// compute the number of high indel events sub-regions
int Region::computeHighIndelEventRegNum(){
	size_t pos, startPosSubReg, endPosSubReg;
	highIndelSubRegNum = 0;
	for(pos=startMidPartPos; pos<=endMidPartPos; pos+=subRegSize){
		startPosSubReg = pos;
		endPosSubReg = startPosSubReg + subRegSize - 1;
		if(endPosSubReg>endMidPartPos) endPosSubReg = endMidPartPos;
		if(computeReadIndelEventNumReg(startPosSubReg, endPosSubReg)>=5) highIndelSubRegNum ++;
	}
	return 0;
}

// compute the total indel number in reads, including insertions, deletions and clippings
// return: the total number of the above indel events
size_t Region::computeReadIndelEventNumReg(size_t startPosReg, size_t endPosReg){
	size_t i, total = 0;
	Base *base;
	for(i=startPosReg; i<=endPosReg; i++){
		base = regBaseArr + i - startRPos;
		total += base->insVector.size() + base->delVector.size() + base->clipVector.size();
	}
	return total;
}

// detect SNVs for the region
void Region::detectSNV(){
	size_t startCheckPos, endCheckPos;
	bool SNV_flag;
	vector<size_t>::iterator disagr;
	for(disagr=disagrePosVector.begin(); disagr!=disagrePosVector.end();){
		startCheckPos = (*disagr) - subRegSize;
		endCheckPos = (*disagr) + subRegSize;
		if(startCheckPos<startRPos) startCheckPos = startRPos;
		if(endCheckPos>endRPos) endCheckPos = endRPos;

//		if(*disagr==18486270)
//			cout << *disagr << endl;

		SNV_flag = computeSNVFlag(*disagr, startCheckPos, endCheckPos);
		if(SNV_flag){
			if(haveMuchShortIndelsAround(startCheckPos, endCheckPos))  // further check around
				indelVector.push_back(allocateReg(*disagr, *disagr));
			else{
				snvVector.push_back(*disagr); // add the SNV
				disagr = disagrePosVector.erase(disagr);
			}
		}else disagr ++;
	}
	snvVector.shrink_to_fit();
}

// compute the SNV flag for a base position
bool Region::computeSNVFlag(size_t pos, size_t startCheckPos, size_t endCheckPos){
	bool SNV_flag = true;
	double cov_tmp;

	if(isInReg(pos, indelVector)) SNV_flag = false;

	if(SNV_flag) // check the maxBase and idx_ref
		if(regBaseArr[pos-startRPos].coverage.idx_max==regBaseArr[pos-startRPos].coverage.idx_RefBase
			or (double)regBaseArr[pos-startRPos].coverage.num_max/regBaseArr[pos-startRPos].coverage.num_bases[5]<MIN_RATIO_SNV)
			SNV_flag = false;

	// check the coverage
	if(SNV_flag){
		cov_tmp = computeMeanCovReg(startCheckPos, endCheckPos);
		if((cov_tmp/meanBlockCov<0.6 and cov_tmp/localRegCov<0.6) or (cov_tmp/meanBlockCov>2 and cov_tmp/localRegCov>2))
			SNV_flag = false;
	}
	return SNV_flag;
}

// determine whether there are much short indel events around
bool Region::haveMuchShortIndelsAround(size_t startCheckPos, size_t endCheckPos){
	bool flag = false;
	Base *base;

	for(size_t i=startCheckPos; i<=endCheckPos; i++){
		base = regBaseArr + i - startRPos;
		if(base->insVector.size()>=paras->min_ins_num_filt or base->delVector.size()>=paras->min_del_num_filt or base->clipVector.size()>=paras->min_clip_num_filt
			or base->num_shortIns>=paras->min_ins_num_filt or base->num_shortdel>=paras->min_del_num_filt or base->num_shortClip>=paras->min_clip_num_filt
			/*or base->getLargerInsNum(paras->min_ins_size_filt)>0 or base->getLargerDelNum(paras->min_del_size_filt)>0 or base->getLargerClipNum(paras->min_clip_size_filt)>0*/){
			flag = true;
			break;
		}
	}

	return flag;
}

// get the number of disagreements excluding the SNV items
size_t Region::getDisagreeNum(){
	return disagrePosVector.size();
}

// get the number of SNVs
size_t Region::getSNVNum(){
	return snvVector.size();
}

// get the number of Indels
size_t Region::getIndelNum(){
	return indelVector.size();
}

// detect the indel regions
void Region::detectIndelReg(){
	reg_t *reg = NULL;
	int32_t i = startMidPartPos - subRegSize;
	if(i<1) i = 1;
	while(i<(int32_t)endMidPartPos){
		//if(i>24935500)  //109500, 5851000, 11812500, 12601500, 14319500, 14868000, 18343500, <18710000>
		//	cout << i << endl;

		reg = getIndelReg(i);
		if(reg){
			indelVector.push_back(reg);
			i = reg->endRefPos + subRegSize;
		} else break;
	}
	indelVector.shrink_to_fit();
}

// get the indel region given the start checking chromosome position
reg_t* Region::getIndelReg(size_t startCheckPos){
	reg_t *reg = NULL;
	size_t i, checkPos, reg_size1, reg_size2, num1, num3, num4, extendSize;
	vector<double> num_vec;
	double high_indel_clip_ratio;
	int32_t startPos1, endPos1, startPos2;
	int32_t startPos_indel = -1, endPos_indel = -1;
	bool indel_beg_flag, indel_end_flag;

	checkPos = startCheckPos;
	while(checkPos<=endMidPartPos){
		startPos_indel = -1;
		endPos_indel = -1;
		indel_beg_flag = false;
		indel_end_flag = false;

		// get normal region1
		startPos1 = -1;
		reg_size1 = 0;
		for(i=checkPos; i<=endMidPartPos; i++){
			if(regBaseArr[i-startRPos].coverage.idx_RefBase==4){ // skip the Ns region
				startPos1 = -1;
				reg_size1 = 0;
				continue;
			}

			if(haveNoAbSigs(regBaseArr+i-startRPos, i)){
				if(startPos1==-1) startPos1 = i;
				reg_size1 ++;
			}else{
				if(reg_size1>=subRegSize) {
					indel_beg_flag = true;
					break;
				}else{
					startPos1 = -1;
					reg_size1 = 0;
				}
			}
		}
		if(indel_beg_flag){
			endPos1 = startPos1 + reg_size1 - 1;

			// get normal region2
			startPos2 = -1;
			reg_size2 = 0;
			extendSize = 0;
			for(i=endPos1+1; i<=endRPos+extendSize; i++){
				if(i>chrlen) break;

				if(regBaseArr[i-startRPos].coverage.idx_RefBase==4){ // skip the Ns region
					startPos2 = -1;
					reg_size2 = 0;
					continue;
				}

				if(haveNoAbSigs(regBaseArr+i-startRPos, i)){
					if(startPos2==-1)
						startPos2 = i;
					++reg_size2;
					if(reg_size2>=2*subRegSize){
						num3 = getLargeIndelNum(startPos2, i);
						if(num3<3){
							indel_end_flag = true;
							break;
						}else{
							startPos2 = -1;
							reg_size2 = 0;
						}
					}
				}else{
					startPos2 = -1;
					reg_size2 = 0;
				}

				if(i==endRPos and indel_end_flag==false)
					extendSize = paras->slideSize * 2;
			}

			if(indel_end_flag){
				startPos_indel = endPos1 + 1;
				endPos_indel = startPos2 - 1;
			}
		}
		// allocate the indel region
		if(indel_beg_flag and indel_end_flag)
			if(endPos_indel-startPos_indel+1>=(int32_t)paras->min_sv_size_usr){
				num1 = getDisZeroCovNum(startPos_indel, endPos_indel);
				//num2 = getMismatchBasesAround(startPos_indel, endPos_indel);
				num3 = getLargeIndelBaseNum(startPos_indel, endPos_indel);
				num_vec = getTotalHighIndelClipRatioBaseNum(startPos_indel, endPos_indel);
				num4 = num_vec.at(0);
				high_indel_clip_ratio = num_vec.at(1);
				//if(num1>0 or num2>=DISAGREE_NUM_THRES_REG or num3>0) {
				if(num1>0 or num3>0 or num4>0 or high_indel_clip_ratio>=HIGH_INDEL_CLIP_BASE_RATIO_THRES) {
					reg = allocateReg(startPos_indel, endPos_indel);
					break;
				}else checkPos = endPos_indel + 1;
			}else checkPos = endPos_indel + 1;
		else break;
	}

	return reg;
}

// determine whether the base have no abnormal signatures
bool Region::haveNoAbSigs(Base *base, size_t pos){
	if(base->isDisagreeBase())
		if(find(snvVector.begin(), snvVector.end(), pos)==snvVector.end()) return false;
	if(base->isZeroCovBase() or base->insVector.size()>=paras->min_ins_num_filt or base->delVector.size()>=paras->min_del_num_filt or base->clipVector.size()>=paras->min_clip_num_filt
		//or base->num_shortIns>=paras->min_ins_num_filt or base->num_shortdel>=paras->min_del_num_filt or base->num_shortClip>=paras->min_clip_num_filt
		/*or base->getLargerInsNum(paras->min_ins_size_filt)>0 or base->getLargerDelNum(paras->min_del_size_filt)>0 or base->getLargerClipNum(paras->min_clip_size_filt)>0*/)
		return false;
//	else if(getMismatchBasesAround(pos-DISAGREE_CHK_REG, pos+DISAGREE_CHK_REG)>=DISAGREE_NUM_THRES_REG)
//		return false;
	else if(base->getLargeIndelNum(paras->large_indel_size_thres)>=3 or (double)(base->getTotalIndelNum()+base->getTotalClipNum())/base->getTotalCovNum()>=HIGH_INDEL_CLIP_RATIO_THRES)
		return false;
	return true;
}

// check [-2, 2] region around
size_t Region::getMismatchBasesAround(size_t pos1, size_t pos2){
	int32_t i, num, startPos, endPos;
	Base *base;

	startPos = pos1;
	if(startPos<(int32_t)startMidPartPos) startPos = startMidPartPos;
	endPos = pos2;
	if(endPos>(int32_t)endMidPartPos) endPos = endMidPartPos;

	for(num=0, i=startPos; i<=endPos; i++){
		base = regBaseArr + i - startRPos;
		if(base->coverage.idx_max!=base->coverage.idx_RefBase or (double)base->coverage.num_max/base->coverage.num_bases[5]<=DISAGREE_THRES_REG)
			num ++;
	}
	return num;
}

size_t Region::getDisZeroCovNum(size_t startPos, size_t endPos){
	size_t i, total = 0;
	for(i=startPos; i<=endPos; i++)
		if(regBaseArr[i-startRPos].isDisagreeBase() or regBaseArr[i-startRPos].isZeroCovBase())
			total ++;
	return total;
}

// get the number of bases with long indels
size_t Region::getLargeIndelBaseNum(size_t startPos, size_t endPos){
	size_t i, large_indel_num, total;
	double ratio;

	total = 0;
	for(i=startPos; i<=endPos; i++){
		large_indel_num = regBaseArr[i-startRPos].getLargeIndelNum(paras->large_indel_size_thres);
		ratio = (double)large_indel_num / regBaseArr[i-startRPos].coverage.num_bases[5];
		if(ratio>=LARGE_INDEL_RATIO_THRES)
			total ++;
	}

	return total;
}

// get the number of bases with long indels
size_t Region::getLargeIndelNum(size_t startPos, size_t endPos){
	size_t i, large_indel_num;
	large_indel_num = 0;
	for(i=startPos; i<=endPos; i++){
		large_indel_num += regBaseArr[i-startRPos].getLargeIndelNum(paras->large_indel_size_thres);
	}
	return large_indel_num;
}

// get the number of high ratio indel bases
vector<double> Region::getTotalHighIndelClipRatioBaseNum(size_t startPos, size_t endPos){
	size_t i, indel_num, clip_num, total_cov, len;
	double ratio, total, total2;
	Base *base;
	vector<double> base_num_vec;

	total = total2 = 0;
	for(i=startPos; i<=endPos; i++){
		base  = regBaseArr + i - startRPos;
		indel_num = base->getTotalIndelNum();
		clip_num = base->getTotalClipNum();
		total_cov = base->getTotalCovNum();
		ratio = (double)(indel_num + clip_num) / total_cov;
		if(ratio>=HIGH_INDEL_CLIP_RATIO_THRES) total ++;
		if(ratio>=SECOND_INDEL_CLIP_RATIO_THRES) total2 ++;
	}

	len = endPos - startPos + 1;
	ratio = (double)total2 / len;

	base_num_vec.push_back(total);
	base_num_vec.push_back(ratio);

	return base_num_vec;
}

// determine the dif type for indel candidate
void Region::determineDifType(){
	if(getDisagreeNum()>0 or zeroCovPosVector.size()>0 or abCovRegVector.size()>0 or highIndelSubRegNum>0)
		indelCandFlag = true;
}

// get SNV vector
vector<size_t> Region::getSnvVector(){
	return snvVector;
}

// get indel vector
vector<reg_t*> Region::getIndelVector(){
	return indelVector;
}

// get indel vector
vector<reg_t*> Region::getClipRegVector(){
	return clipRegVector;
}

// get indel vector
vector<size_t> Region::getZeroCovPosVector(){
	return zeroCovPosVector;
}

//
void Region::detectHighClipReg(){
	int32_t i = startMidPartPos - SUB_CLIP_REG_SIZE;
	if(i<1) i = 1;
	while(i<(int32_t)endMidPartPos){
//		if(i>253000)
//			cout << i << endl;
		reg_t *reg = getClipReg(i);
		if(reg) {
			clipRegVector.push_back(reg);
			i = reg->endRefPos + SUB_CLIP_REG_SIZE;
		} else break;
	}
	clipRegVector.shrink_to_fit();
}

// get the clip region
reg_t* Region::getClipReg(size_t startCheckPos){
	reg_t *reg = NULL;
	size_t checkPos, subclipreg_size, startPos_tmp, endPos_tmp;
	int32_t startPos_clip, endPos_clip;
	bool normal_reg_flag, clip_reg_flag;

	//cout << chrname << ":" << to_string(startRPos) << "-" << to_string(endRPos) << ", localRegCov=" << localRegCov << endl;

	subclipreg_size = SUB_CLIP_REG_SIZE;
	clip_reg_flag = false;

	startPos_clip = endPos_clip = -1;
	checkPos = startCheckPos;
	while(checkPos<=endMidPartPos){
		startPos_tmp = checkPos;
		endPos_tmp = checkPos + subclipreg_size - 1;
		if(endPos_tmp>endMidPartPos) endPos_tmp = endMidPartPos;
		normal_reg_flag = haveNoClipSig(startPos_tmp, endPos_tmp, HIGH_CLIP_RATIO_THRES);
		if(normal_reg_flag==false){ // clip region
			if(clip_reg_flag==false){ // first clip region
				clip_reg_flag = true;
				startPos_clip = startPos_tmp;
				endPos_clip = endPos_tmp;
			}else{
				endPos_clip = endPos_tmp;
			}
		}else { // normal region
			if(clip_reg_flag==true){
				break;
			}
		}
		checkPos = endPos_tmp + 1;
	}

//	// compute clip positions by align data
//	startPos_clip = endPos_clip = -1;
//	if(startPos_clip==-1){
//		subreg_num = (endMidPartPos - startRPos) / SUB_CLIP_REG_SIZE;
//		for(i=0; i<subreg_num; i++){
//			startPos_tmp = startRPos + i * SUB_CLIP_REG_SIZE;
//			endPos_tmp = startPos_tmp + SUB_CLIP_REG_SIZE - 1;
//			if(endPos_tmp>endMidPartPos) endPos_tmp = endMidPartPos;
//			clip_num = 0;
//			for(j=startPos_tmp; j<=endPos_tmp; j++){
//				clip_vec = regBaseArr[j-startRPos].clipVector;
//				for(k=0; k<clip_vec.size(); k++){
//					if(stoi(clip_vec.at(k)->seq)>=MIN_CLIP_END_SIZE)
//						clip_num ++;
//				}
//			}
//			if(localRegCov>0){
//				ratio = clip_num / localRegCov;
//				if(ratio>=HIGH_CLIP_RATIO_THRES){
//					startPos_clip = startPos_tmp;
//					//cout << "####################### startPos_clip=" << startPos_clip << endl;
//					break;
//				}
//			}
//		}
//	}
//	if(startPos_clip==-1){
//		for(i=startRPos; i<=endMidPartPos; i++)
//			if((double)regBaseArr[i-startRPos].clipVector.size()/regBaseArr[i-startRPos].coverage.num_bases[5]>=HIGH_CLIP_RATIO_THRES){
//				startPos_clip = i;
//				break;
//			}
//	}
//
//	if(endPos_clip==-1){
//		subreg_num = (endRPos - startMidPartPos) / SUB_CLIP_REG_SIZE;
//		for(i=0; i<subreg_num; i++){
//			endPos_tmp = endRPos - i * SUB_CLIP_REG_SIZE;
//			startPos_tmp = endPos_tmp - SUB_CLIP_REG_SIZE + 1;
//			if(startPos_tmp<startMidPartPos) startPos_tmp = startMidPartPos;
//			clip_num = 0;
//			for(j=startPos_tmp; j<=endPos_tmp; j++){
//				clip_vec = regBaseArr[j-startRPos].clipVector;
//				for(k=0; k<clip_vec.size(); k++){
//					if(stoi(clip_vec.at(k)->seq)>=MIN_CLIP_END_SIZE)
//						clip_num ++;
//				}
//			}
//			if(localRegCov>0){
//				ratio = clip_num / localRegCov;
//				if(ratio>=HIGH_CLIP_RATIO_THRES){
//					endPos_clip = endPos_tmp;
//					//cout << "=#=#=#=#=#=#=#=#=#=#=#= endPos_clip=" << endPos_clip << endl;
//					break;
//				}
//			}
//		}
//	}
//
//
//	if(endPos_clip==-1){
//		for(i=endRPos; i>=startMidPartPos; i--)
//			if((double)regBaseArr[i-startRPos].clipVector.size()/regBaseArr[i-startRPos].coverage.num_bases[5]>=HIGH_CLIP_RATIO_THRES){
//				endPos_clip = i;
//				break;
//			}
//	}

//	// compute region by align data
//	if(startPos_clip==-1 or endPos_clip==-1){
//		// load the aligned reads data
//		clipAlnDataLoader data_loader(chrname, startRPos, endRPos, paras->inBamFile);
//		data_loader.loadClipAlnData(clipAlnDataVector);
//
//		if(clipAlnDataVector.size()>0){
//			// compute the start position
//			if(startPos_clip==-1){
//				subreg_num = (endMidPartPos - startRPos) / SUB_CLIP_REG_SIZE;
//				for(i=0; i<subreg_num; i++){
//					startPos_tmp = startRPos + i * SUB_CLIP_REG_SIZE;
//					endPos_tmp = startPos_tmp + SUB_CLIP_REG_SIZE - 1;
//					if(endPos_tmp>endMidPartPos) endPos_tmp = endMidPartPos;
//					total_reads_num = clip_num = 0;
//					for(j=0; j<clipAlnDataVector.size(); j++){
//						clip_aln = clipAlnDataVector.at(j);
//						if((clip_aln->startRefPos>=startPos_tmp and clip_aln->startRefPos<=endPos_tmp)
//							or (clip_aln->endRefPos>=startPos_tmp and clip_aln->endRefPos<=endPos_tmp)
//							or (startPos_tmp>=clip_aln->startRefPos and startPos_tmp<=clip_aln->endRefPos)
//							or (endPos_tmp>=clip_aln->startRefPos and endPos_tmp<=clip_aln->endRefPos)){ // overlap
//							total_reads_num ++;
//							if((clip_aln->leftClipSize>=MIN_CLIP_END_SIZE and clip_aln->startRefPos>=startPos_tmp and clip_aln->startRefPos<=endPos_tmp)
//								or (clip_aln->rightClipSize>=MIN_CLIP_END_SIZE and clip_aln->endRefPos>=startPos_tmp and clip_aln->endRefPos<=endPos_tmp))
//								clip_num ++;
//						}
//					}
//					if(total_reads_num>0){
//						ratio = (double)clip_num / total_reads_num;
//						if(ratio>=HIGH_CLIP_RATIO_THRES){
//							startPos_clip = startPos_tmp;
//							cout << "####################### startPos_clip=" << startPos_clip << endl;
//							break;
//						}
//					}
//				}
//			}
//
//			// compute the start position
//			if(endPos_clip==-1){
//				subreg_num = (endRPos - startMidPartPos) / SUB_CLIP_REG_SIZE;
//				for(i=0; i<subreg_num; i++){
//					endPos_tmp = endRPos - i * SUB_CLIP_REG_SIZE;
//					startPos_tmp = endPos_tmp - SUB_CLIP_REG_SIZE + 1;
//					if(startPos_tmp<startMidPartPos) startPos_tmp = startMidPartPos;
//					total_reads_num = clip_num = 0;
//					for(k=(int32_t)clipAlnDataVector.size()-1; k>=0; k--){
//						clip_aln = clipAlnDataVector.at(k);
//						if((clip_aln->startRefPos>=startPos_tmp and clip_aln->startRefPos<=endPos_tmp)
//							or (clip_aln->endRefPos>=startPos_tmp and clip_aln->endRefPos<=endPos_tmp)
//							or (startPos_tmp>=clip_aln->startRefPos and startPos_tmp<=clip_aln->endRefPos)
//							or (endPos_tmp>=clip_aln->startRefPos and endPos_tmp<=clip_aln->endRefPos)){ // overlap
//							total_reads_num ++;
//							if((clip_aln->leftClipSize>=MIN_CLIP_END_SIZE and clip_aln->startRefPos>=startPos_tmp and clip_aln->startRefPos<=endPos_tmp)
//								or (clip_aln->rightClipSize>=MIN_CLIP_END_SIZE and clip_aln->endRefPos>=startPos_tmp and clip_aln->endRefPos<=endPos_tmp))
//								clip_num ++;
//						}
//					}
//					if(total_reads_num>0){
//						ratio = (double)clip_num / total_reads_num;
//						if(ratio>=HIGH_CLIP_RATIO_THRES){
//							endPos_clip = endPos_tmp;
//							cout << "=#=#=#=#=#=#=#=#=#=#=#= endPos_clip=" << endPos_clip << endl;
//							break;
//						}
//					}
//				}
//			}
//		}
//	}

	// check disagreements
	if(startPos_clip!=-1 and endPos_clip!=-1){
		// compute the number of disagreements
		int32_t disagreeNum = computeDisagreeNum(regBaseArr+startPos_clip-startRPos, endPos_clip-startPos_clip+1);
		normal_reg_flag = haveNoClipSig(startPos_clip, endPos_clip, HIGH_CLIP_RATIO_THRES*3);

		if(disagreeNum>=1 or normal_reg_flag==false) reg = allocateReg(startPos_clip, endPos_clip);
	}

	return reg;
}


bool Region::haveNoClipSig(size_t startPos, size_t endPos, double clip_ratio_thres){
	bool flag = true;
	size_t i, j, clip_num;
	vector<clipEvent_t*> clip_vec;
	double ratio;

	clip_num = 0;
	for(i=startPos; i<=endPos; i++){
		clip_vec = regBaseArr[i-startRPos].clipVector;
		for(j=0; j<clip_vec.size(); j++){
			if(stoi(clip_vec.at(j)->seq)>=MIN_CLIP_END_SIZE)
				clip_num ++;
		}
	}

	if(localRegCov>0){
		ratio = clip_num / localRegCov;
		if(ratio>=clip_ratio_thres) flag = false;
	}

	if(flag){
		for(i=startPos; i<=endPos; i++)
			if((double)regBaseArr[i-startRPos].clipVector.size()/regBaseArr[i-startRPos].coverage.num_bases[5]>=clip_ratio_thres){
				flag = false;
				break;
			}
	}

	return flag;
}

