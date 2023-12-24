#include "covLoader.h"
#include "util.h"

covLoader::covLoader(string &chrname, int64_t startPos, int64_t endPos, faidx_t *fai) {
	this->chrname = chrname;
	this->startPos = startPos;
	this->endPos = endPos;
	this->fai = fai;
	min_ins_size_filt = min_del_size_filt = 0;
	left_ref_base = right_ref_base = '-';
	bam_type = BAM_INVALID;
}

covLoader::covLoader(string &chrname, int64_t startPos, int64_t endPos, faidx_t *fai, int32_t min_ins_size_filt, int32_t min_del_size_filt) {
	this->chrname = chrname;
	this->startPos = startPos;
	this->endPos = endPos;
	this->fai = fai;
	this->min_ins_size_filt = min_ins_size_filt;
	this->min_del_size_filt = min_del_size_filt;
	left_ref_base = right_ref_base = '-';
	bam_type = BAM_INVALID;
}

covLoader::~covLoader() {
}

// initialize the base array of the block
Base *covLoader::initBaseArray(){
	Base *baseArray = new Base[endPos-startPos+1]();
	if(!baseArray){
		cerr << __func__ << ": cannot allocate memory" << endl;
		exit(1);
	}
	// assign the ref base index in base coverage
	assignRefBase(baseArray, fai);
	assignPolymerFlag(baseArray);
	return baseArray;
}

void covLoader::freeBaseArray(Base *baseArray){
	if(baseArray){
		delete[] baseArray;
		baseArray = NULL;
	}
}

int covLoader::assignRefBase(Base *baseArray, faidx_t *fai){
	int64_t pos, start_pos, end_pos, chrlen, left_ext_size, right_ext_size;

	chrlen = faidx_seq_len(fai, chrname.c_str()); // get the reference length
	if(startPos==1) { start_pos = startPos;	left_ext_size = 0; }	// first base is the left base
	else { start_pos = startPos - 1; left_ext_size = 1; }
	if(endPos==chrlen) { end_pos = endPos; right_ext_size = 0; }	// last base is the right base
	else { end_pos = endPos + 1; right_ext_size = 1; }

	// get the region sequence
	string reg = chrname + ":" + to_string(start_pos) + "-" + to_string(end_pos);

	RefSeqLoader refseq_loader(reg, fai);
	refseq_loader.getRefSeq();

	// assign the left and right base around the region
	if(left_ext_size==1) left_ref_base = refseq_loader.refseq[0];
	else left_ref_base = '-';
	if(right_ext_size==1) right_ref_base = refseq_loader.refseq[refseq_loader.refseq_len-1];
	else right_ref_base = '-';

	// assign reference bases
	pos = startPos;
	for(int i=0; i<refseq_loader.refseq_len-left_ext_size-right_ext_size; i++){
		switch(refseq_loader.refseq[i+left_ext_size]){
			case 'A':
			case 'a': baseArray[i].coverage.idx_RefBase = 0; break;
			case 'C':
			case 'c': baseArray[i].coverage.idx_RefBase = 1; break;
			case 'G':
			case 'g': baseArray[i].coverage.idx_RefBase = 2; break;
			case 'T':
			case 't': baseArray[i].coverage.idx_RefBase = 3; break;
			case 'N':
			case 'n': baseArray[i].coverage.idx_RefBase = 4; break;
			case 'M':
			case 'm':
			case 'R':
			case 'r':
			case 'S':
			case 's':
			case 'V':
			case 'v':
			case 'W':
			case 'w':
			case 'Y':
			case 'y':
			case 'H':
			case 'h':
			case 'K':
			case 'k':
			case 'D':
			case 'd':
			case 'B':
			case 'b': baseArray[i].coverage.idx_RefBase = 5; break;
			default: cerr << __func__ << ": unknown base: " << refseq_loader.refseq[i+left_ext_size] << " at location: " << chrname << ":" << pos << endl; exit(1);
		}
		baseArray[i].coverage.refBase = refseq_loader.refseq[i+left_ext_size];
		pos ++;
	}
	return 0;
}

// assign polymer flag of bases
void covLoader::assignPolymerFlag(Base *baseArray){
	Base *base, *left_base, *right_base;
	int64_t i, len = endPos - startPos + 1;
	for(i=0; i<len; i++){
		base = baseArray + i;
		if(i>0 and i<len-1){ // inner items
			left_base = baseArray + i - 1;
			right_base = baseArray + i + 1;
			if(left_base->coverage.refBase==base->coverage.refBase or right_base->coverage.refBase==base->coverage.refBase)
				base->coverage.polymer_flag = true;
		}else if(i==0 and i<len-1){ // first item
			right_base = baseArray + i + 1;
			if(right_base->coverage.refBase==base->coverage.refBase or base->coverage.refBase==left_ref_base)
				base->coverage.polymer_flag = true;
		}else if(i>0 and i==len-1){ // last item
			left_base = baseArray + i - 1;
			if(left_base->coverage.refBase==base->coverage.refBase or base->coverage.refBase==right_ref_base)
				base->coverage.polymer_flag = true;
		}
	}
}

// assign base coverage
void covLoader::generateBaseCoverage(Base *baseArr, vector<bam1_t*> &alnDataVector){
	vector<struct alnSeg*> alnSegs;
	bam1_t *b;
	string qname,qname_com;

	if(alnDataVector.empty()) return; // tolerate zero coverage regions

	bam_type = getBamType(alnDataVector.at(0));
	if(bam_type==BAM_INVALID){
		cerr << __func__ << ": unknown bam type, error!" << endl;
		exit(1);
	}
	for(size_t i=0; i<alnDataVector.size(); i++){
		b = alnDataVector.at(i);
		if(!(b->core.flag & BAM_FUNMAP)){ // aligned
//			qname = bam_get_qname(b);
//			cin>>qname_com;
//			qname_com = "b219a8ba_98667_6199";
//			if(qname.compare(qname_com)==0)
//				cout << qname << endl;
//			else continue;

			switch(bam_type){
				case BAM_CIGAR_NO_DIFF_MD:
					//alnSegs = generateAlnSegs(b);
					alnSegs = generateAlnSegs2(b, startPos, endPos);
					break;
				case BAM_CIGAR_NO_DIFF_NO_MD:
				case BAM_CIGAR_DIFF_MD:
				case BAM_CIGAR_DIFF_NO_MD:
					//alnSegs = generateAlnSegs_no_MD(b, baseArr, startPos, endPos);
					alnSegs = generateAlnSegs_no_MD2(b, baseArr, startPos, endPos);
					break;
				default:
					cerr << __func__ << ": unknown bam type, error!" << endl;
					exit(1);
			}// generate align segments

			updateBaseInfo(baseArr, alnSegs); // update base information
			destroyAlnSegs(alnSegs);
		}
	}

	updateBaseCovInfo(baseArr);

	// compute number of deletions
	computeDelNumFromDelVec(baseArr);

	// compute consensus indel events
	computeConIndelEventRatio(baseArr);
}

// update the block base array information according to read alignments
int covLoader::updateBaseInfo(Base *baseArr, vector<struct alnSeg*> &alnSegs){
	baseCoverage_t *cover;
	int64_t pos, epos, tmp_endPos, misbase;
	int32_t idx, endflag, position, len;
	vector<struct alnSeg*>::iterator seg;
	string str;

	for(seg=alnSegs.begin(); seg!=alnSegs.end(); seg++){
		switch((*seg)->opflag){
//			case BAM_CMATCH:
//				//if((*seg)->startRpos>=startPos and (*seg)->startRpos<=endPos){
//				if(isOverlappedPos((*seg)->startRpos, (*seg)->startRpos+(*seg)->seglen-1, startPos, endPos)){
//					tmp_endPos = (*seg)->startRpos + (*seg)->seglen - 1;
//					for(pos=(*seg)->startRpos; pos<=tmp_endPos; pos++){
//						if(pos>=startPos and pos<=endPos){
//							cover = &(baseArr[pos-startPos].coverage);
//							if(cover->idx_RefBase>=0 and cover->idx_RefBase<=4)
//								cover->num_bases[cover->idx_RefBase] ++;
//							else{
//								cerr << __func__ << ", line=" << __LINE__ << ": invalid idx_RefBase " << cover->idx_RefBase << endl;
//								outputAlnSegs(alnSegs);
//								exit(1);
//							}
//						}
//					}
//				}
//				break;
//			//case MD_MISMATCH:  // single base mismatch
//			case BAM_CDIFF:  // single base mismatch
//					if((*seg)->startRpos>=startPos and (*seg)->startRpos<=endPos){
//						cover = &(baseArr[(*seg)->startRpos-startPos].coverage);
//						misbase = (*seg)->seg_MD.at(0);
//						switch(misbase){
//							case 'a':
//							case 'A': idx = 0; break;
//							case 'c':
//							case 'C': idx = 1; break;
//							case 'g':
//							case 'G': idx = 2; break;
//							case 't':
//							case 'T': idx = 3; break;
//							case 'n':
//							case 'N': idx = 4; break;
//							default:
//								cerr << __func__ << ", line=" << __LINE__ << ": invalid base " << int(misbase) << endl;
//								exit(1);
//						}
//						if(idx!=cover->idx_RefBase) cover->num_bases[idx] ++;
//						else if(idx==4) cover->num_bases[idx] ++;   // tolerate the mismatched base 'N'
//						else{
//							cerr << __func__ << ", line=" << __LINE__ << ": invalid array idx=" << idx << endl;
//							exit(1);
//						}
//					}
//				break;
			case BAM_CINS:  // insertion in query
				if((*seg)->startRpos>=startPos and (*seg)->startRpos<=endPos){
					if((*seg)->seglen>=min_ins_size_filt)
						baseArr[(*seg)->startRpos-startPos].addInsEvent(allocateInsEvent((*seg)->startRpos, (*seg)->seg_MD));
					else
						baseArr[(*seg)->startRpos-startPos].num_shortIns ++;
				}
				break;
			case BAM_CDEL:  // deletion in query
				//if((*seg)->startRpos>=startPos and (*seg)->startRpos<=endPos){
				if(isOverlappedPos((*seg)->startRpos, (*seg)->startRpos+(*seg)->seglen-1, startPos, endPos)){
					if((*seg)->startRpos>=startPos) pos = (*seg)->startRpos;
					else pos = startPos;
					if((*seg)->seglen>=min_del_size_filt){
						if((*seg)->startRpos<startPos or (*seg)->startRpos+(*seg)->seglen-1 > endPos){
							position = pos - (*seg)->startRpos;
							if((*seg)->startRpos+(*seg)->seglen-1>=endPos) epos = endPos;
							else epos = (*seg)->startRpos+(*seg)->seglen-1;
							len = epos - pos + 1;
							str = (*seg)->seg_MD.substr(position, len);
							baseArr[pos-startPos].addDelEvent(allocateDelEvent(pos, str));
						}
						else
							baseArr[pos-startPos].addDelEvent(allocateDelEvent((*seg)->startRpos, (*seg)->seg_MD));
					}else
						baseArr[pos-startPos].num_shortdel ++;

//					tmp_endPos = (*seg)->startRpos + (*seg)->seglen - 1;
//					for(pos=(*seg)->startRpos; pos<=tmp_endPos; pos++){
//						if(pos>=startPos and pos<=endPos){
//							if((*seg)->seglen>=min_del_size_filt){
//
//								baseArr[pos-startPos].addDelEvent(allocateDelEvent((*seg)->startRpos, (*seg)->seg_MD));
//							}else
//								baseArr[pos-startPos].num_shortdel ++;
//							break;
//						}
//					}
				}
				break;
			case BAM_CSOFT_CLIP:  // soft clipping in query
			case BAM_CHARD_CLIP:  // hard clipping in query
				if((*seg)->startRpos>=startPos and (*seg)->startRpos<=endPos){
					if((*seg)->startQpos==1) endflag = 0;
					else endflag = 1;
					baseArr[(*seg)->startRpos-startPos].addClipEvent(allocateClipEvent((*seg)->startRpos, (*seg)->opflag, endflag, (*seg)->seg_MD));
				}
				break;
			case BAM_CEQUAL:
				//if((*seg)->startRpos>=startPos and (*seg)->startRpos<=endPos){
				if(isOverlappedPos((*seg)->startRpos, (*seg)->startRpos+(*seg)->seglen-1, startPos, endPos)){
					tmp_endPos = (*seg)->startRpos + (*seg)->seglen - 1;
					for(pos=(*seg)->startRpos; pos<=tmp_endPos; pos++){
						if(pos>=startPos and pos<=endPos){
							cover = &(baseArr[pos-startPos].coverage);
							if(cover->idx_RefBase>=0 and cover->idx_RefBase<=4)
								cover->num_bases[cover->idx_RefBase] ++;
							else{
								cerr << __func__ << ", line=" << __LINE__ << ": invalid idx_RefBase " << cover->idx_RefBase << endl;
								outputAlnSegs(alnSegs);
								exit(1);
							}
						}
					}
				}
				break;
			case BAM_CDIFF:
				if((*seg)->startRpos>=startPos and (*seg)->startRpos<=endPos){
					tmp_endPos = (*seg)->startRpos + (*seg)->seglen - 1;
					for(pos=(*seg)->startRpos; pos<=tmp_endPos; pos++){
						if(pos>=startPos and pos<=endPos){
							cover = &(baseArr[pos-startPos].coverage);
							misbase = (*seg)->seg_MD.at(pos-(*seg)->startRpos);
							switch(misbase){
								case 'a':
								case 'A': idx = 0; break;
								case 'c':
								case 'C': idx = 1; break;
								case 'g':
								case 'G': idx = 2; break;
								case 't':
								case 'T': idx = 3; break;
								case 'n':
								case 'N': idx = 4; break;
								default:
									cerr << __func__ << ", line=" << __LINE__ << ": invalid base " << int(misbase) << endl;
									exit(1);
							}
							if(idx!=cover->idx_RefBase) cover->num_bases[idx] ++;
							else if(idx==4) cover->num_bases[idx] ++;   // tolerate the mismatched base 'N'
							else{
								cerr << __func__ << ", line=" << __LINE__ << ": invalid array idx=" << idx << endl;
								exit(1);
							}
						}
					}
				}
				break;
			case BAM_CREF_SKIP:  // unexpected events
			case BAM_CPAD:
			case BAM_CBACK:
				cerr << __func__ << ", line=" << __LINE__ << ": invalid opflag " << (*seg)->opflag << endl;
				exit(1);
				break;
			default:
				break;
		}
	}

	return 0;
}

// update the coverage information for each base in this block
void covLoader::updateBaseCovInfo(Base *baseArr){
	int64_t pos;
	for(pos=startPos; pos<=endPos; pos++){
		baseArr[pos-startPos].insVector.shrink_to_fit();
		baseArr[pos-startPos].delVector.shrink_to_fit();
		baseArr[pos-startPos].clipVector.shrink_to_fit();
		baseArr[pos-startPos].updateCovInfo();
	}
}

// compute number of deletions
void covLoader::computeDelNumFromDelVec(Base *baseArr){
	Base *base;
	delEvent_t *del_event;

	int64_t j, pos, pos_tmp, len;
	for(pos=startPos; pos<=endPos; pos++){
		base = baseArr + pos - startPos;
		for(size_t i=0; i<base->delVector.size(); i++){
			del_event = base->delVector.at(i);
			len = del_event->seq.size();
			for(j=0; j<len; j++) {
				pos_tmp = pos + j;
				if(pos_tmp<=endPos) baseArr[pos_tmp-startPos].del_num_from_del_vec ++;
				else break;
			}
			for(j=1; j<len; j++) { // 20220531
				pos_tmp = pos + j;
				if(pos_tmp<=endPos) baseArr[pos_tmp-startPos].extendDelVector.push_back(del_event);
				else break;
			}
			baseArr[pos_tmp-startPos].extendDelVector.shrink_to_fit();
		}
	}
}

// compute consensus indel event ratio
void covLoader::computeConIndelEventRatio(Base *baseArr){
	struct indelCountNode{
		indelEvent_t *indel_item;
		uint32_t count;
	};

	size_t i, j;
	int64_t pos, total_cov, maxValue, maxValue_ins, num_del, max_con_type;
	Base *base;
	insEvent_t *ins_event;
	vector<struct indelCountNode*> insCount_vec;
	vector<delEvent_t*> del_vec;
	bool polymer_flag1, polymer_flag2;
	struct indelCountNode *ins_count, *indel_count_item;

	for(pos=startPos; pos<=endPos; pos++){
		// insertions
		base = baseArr + pos - startPos;
		for(i=0; i<base->insVector.size(); i++){
			ins_event = base->insVector.at(i);
			indel_count_item = NULL;
			for(j=0; j<insCount_vec.size(); j++){
				ins_count = insCount_vec.at(j);
				if(ins_event->seq.compare(ins_count->indel_item->seq)==0){ // identical
					indel_count_item = ins_count;
					break;
				}else{ // check polymer
					polymer_flag1 = isPolymerSeq(ins_event->seq);
					polymer_flag2 = isPolymerSeq(ins_count->indel_item->seq);
					if(polymer_flag1 and polymer_flag2 and ins_event->seq.at(0)==ins_count->indel_item->seq.at(0)){
						indel_count_item = ins_count;
						break;
					}
				}
			}
			if(indel_count_item){ // found
				indel_count_item->count ++;
			}else{ // not found, then create it
				indel_count_item = new struct indelCountNode();
				indel_count_item->indel_item = ins_event;
				indel_count_item->count = 1;
				insCount_vec.push_back(indel_count_item);
			}
		}

		// compute maximum and total count
		maxValue_ins = 0;
		for(i=0; i<insCount_vec.size(); i++){
			ins_count = insCount_vec.at(i);
			if(ins_count->count>maxValue_ins) maxValue_ins = ins_count->count;
		}
		num_del = base->del_num_from_del_vec + base->num_shortdel;

		if(maxValue_ins==0 and num_del==0) { maxValue = 0; max_con_type = BASE_INDEL_CON_UNUSED; } // unused flag
		else if(maxValue_ins>num_del) { maxValue = maxValue_ins; max_con_type = BAM_CINS; }  // insertion majority
		else { maxValue = num_del; max_con_type = BAM_CDEL; }  // deletion majority
		base->maxConIndelEventNum = maxValue;
		base->max_con_type = max_con_type;
		total_cov = base->getTotalCovNum() + num_del;  // coverage and shadow coverage
		if(total_cov>0) base->maxConIndelEventRatio = (float)maxValue / total_cov;
		else base->maxConIndelEventRatio = 0;

		// release items
		for(j=0; j<insCount_vec.size(); j++) delete insCount_vec.at(j);
		insCount_vec.clear();
	}
}

