#include "genotyping.h"

genotyping::genotyping(reg_t *reg, faidx_t *fai, string &inBamFile, int32_t sig_size_thres, double size_ratio_match_thres, double min_dip_ratio_thres, double max_dip_ratio_thres, int32_t min_sup_num_recover){
	this->reg = reg;
	this->fai = fai;
	this->inBamFile = inBamFile;
	this->sig_size_thres = sig_size_thres;
	clip_size_thres = 50;
	clip_extend_match_thres = 200;
	this->size_ratio_match_thres = size_ratio_match_thres;
	valid_summed_size_ratio_read = 0.5;
	noProfileQueryNum = 0;
	this->min_alle_ratio_thres = min_dip_ratio_thres;
	this->max_alle_ratio_thres = max_dip_ratio_thres;
	this->min_sup_num_recover = min_sup_num_recover;
	seed_gtQuery = NULL;
}

genotyping::~genotyping(){
	destroyQueryGtSigVec(queryGtSig_vec);
	destroyMatchProfilePatVec(match_profile_pat_vec);
}

// free match profile pattern memory
void genotyping::destroyMatchProfilePatVec(vector<profile_pat_t*> &match_profile_pat_vec){
	for(size_t i=0; i<match_profile_pat_vec.size(); i++)
		delete match_profile_pat_vec.at(i);
	vector<profile_pat_t*>().swap(match_profile_pat_vec);
}

// compute genotype of variants
void genotyping::computeGenotype(){
	alnDataLoader data_loader(reg->chrname, reg->startRefPos, reg->endRefPos, inBamFile);
	data_loader.loadAlnData(alnDataVector);

	filterInvalidAlnData(alnDataVector, valid_summed_size_ratio_read);

	// extract genotyping signatures
	queryGtSig_vec = extractGtSigVec();

	// matching
	gtSigMatch(queryGtSig_vec);

	// compute valid gtSig columns
	validGtSigFlagVec = computeValidGtSigFlags(seed_gtQuery, reg);

	// compute genotype string
	gt_str_vec = computeGenotypeString(queryGtSig_vec, validGtSigFlagVec);
	reg->gt_seq = gt_str_vec.at(0) + "\t" + gt_str_vec.at(1);
	//cout << "gt_str_vec: " << gt_str_vec.at(0) << "\t" << gt_str_vec.at(1) << endl;

//	if(match_profile_pat_vec.size()>0 and match_profile_pat_vec.size()!=seed_gtQuery->gtSig_vec.size()){
//		cout << "match_profile_pat_vec.size=" << match_profile_pat_vec.size() << ", seed_gtQuery->gtSig_vec.size=" << seed_gtQuery->gtSig_vec.size() << endl;
//	}

	// recover variants
	if(queryGtSig_vec.size()>min_sup_num_recover)
		recoverVariants(reg, match_profile_pat_vec, seed_gtQuery, queryGtSig_vec, validGtSigFlagVec);

	data_loader.freeAlnData(alnDataVector);
}

// filter invalid short reads
void genotyping::filterInvalidAlnData(vector<bam1_t*> &alnDataVector, double valid_summed_size_ratio_read){
	bam1_t *b, *b2;
	string queryname, queryname2;
	int64_t start_pos, end_pos, aln_size, querylen, leftClipSize, rightClipSize;
	double len_ratio;
	vector<bam1_t*> query_aln_segs;

	for(size_t i=0; i<alnDataVector.size(); ){
		b = alnDataVector.at(i);
		queryname = bam_get_qname(b);
		start_pos = b->core.pos;
		end_pos = bam_endpos(b);

		query_aln_segs = getQueryAlnSegs(alnDataVector, queryname);
		aln_size = getAlnSizeSingleQuery(query_aln_segs);
		querylen = getOriginalQueryLen(b);
		len_ratio = (double)aln_size / querylen;

		if(len_ratio<valid_summed_size_ratio_read){ // invalid segments, then delete them

			//cout << "--------- deleted [" << i << "], queryname=" << queryname << ", querylen=" << querylen << ", align_size=" << aln_size << ", len_ratio=" << len_ratio << ", start_pos=" << start_pos << ", end_pos=" << end_pos << endl;

			for(size_t j=i; j<alnDataVector.size(); ){
				b2 = alnDataVector.at(j);
				queryname2 = bam_get_qname(b2);
				if(queryname.compare(queryname2)==0){
					bam_destroy1(b2);
					alnDataVector.erase(alnDataVector.begin()+j);
				}else j++;
			}
		}else{
			//cout << "[" << i << "], queryname=" << queryname << ", querylen=" << querylen << ", align_size=" << aln_size << ", len_ratio=" << len_ratio << ", start_pos=" << start_pos << ", end_pos=" << end_pos << endl;
			i++;
		}
	}
}

vector<bam1_t*> genotyping::getQueryAlnSegs(vector<bam1_t*> &alnDataVector, string &queryname){
	vector<bam1_t*> query_aln_segs;
	bam1_t *b;
	string qname;

	for(size_t i=0; i<alnDataVector.size(); i++){
		b = alnDataVector.at(i);
		qname = bam_get_qname(b);
		if(qname.compare(queryname)==0) query_aln_segs.push_back(b);
	}

	return query_aln_segs;
}

int32_t genotyping::getAlnSizeSingleQuery(vector<bam1_t*> &query_aln_segs){
	bam1_t *b;
	int32_t total = 0;
	vector<int32_t> aln_size_vec;

	for(size_t i=0; i<query_aln_segs.size(); i++){
		b = query_aln_segs.at(i);
		aln_size_vec = getAlnSizeSingleSeg(b);
		total += aln_size_vec.at(0);
	}

	return total;
}

// get align size vevtor: [0]-align size, [1]-querylen, [2]-left clip size, [3]-right clip size
vector<int32_t> genotyping::getAlnSizeSingleSeg(bam1_t *b){
	int32_t querylen, aln_size, left_clip_size, right_clip_size;
	uint32_t *c, op1, op2;
	vector<int32_t> aln_size_vec;

	c = bam_get_cigar(b);  // CIGAR
	// left clip
	op1 = bam_cigar_op(c[0]);
	if(op1==BAM_CSOFT_CLIP or op1==BAM_CHARD_CLIP) left_clip_size = bam_cigar_oplen(c[0]);
	else left_clip_size = 0;
	// left clip
	op2 = bam_cigar_op(c[b->core.n_cigar-1]);
	if(op2==BAM_CSOFT_CLIP or op2==BAM_CHARD_CLIP) right_clip_size = bam_cigar_oplen(c[b->core.n_cigar-1]);
	else right_clip_size = 0;

	// querylen
	querylen = b->core.l_qseq;
	if(op1==BAM_CHARD_CLIP)
		querylen += left_clip_size;
	if(op2==BAM_CHARD_CLIP)
		querylen += right_clip_size;

	// align size
	aln_size = querylen - left_clip_size - right_clip_size;

	aln_size_vec.push_back(aln_size);
	aln_size_vec.push_back(querylen);
	aln_size_vec.push_back(left_clip_size);
	aln_size_vec.push_back(right_clip_size);

	return aln_size_vec;
}

// extract CIGAR signatures
vector<queryGtSig_t*> genotyping::extractGtSigVec(){
	bam1_t *b;
	queryGtSig_t *queryGtSig;
	vector<gtSig_t*> sigVec_singleQuery;
	vector<struct alnSeg*> alnSegs;
	int32_t bam_type, seq_len;
	int64_t startPos, endPos, chrlen_tmp;
	string refseq, reg_str;
	char *seq;

	if(alnDataVector.size()>0){
		chrlen_tmp = faidx_seq_len(fai, reg->chrname.c_str()); // get the reference length

	//	startPos = reg->startRefPos - VAR_ALN_EXTEND_SIZE;
	//	endPos = reg->endRefPos + VAR_ALN_EXTEND_SIZE;
		startPos = reg->startRefPos - clip_extend_match_thres;
		endPos = reg->endRefPos + clip_extend_match_thres;
		if(startPos<1) startPos = 1;
		if(endPos>chrlen_tmp) endPos = chrlen_tmp;

		reg_str = reg->chrname + ":" + to_string(startPos) + "-" + to_string(endPos);
		seq = fai_fetch(fai, reg_str.c_str(), &seq_len);
		refseq = seq;
		free(seq);

		bam_type = getBamType(alnDataVector.at(0));
		if(bam_type==BAM_INVALID){
			cerr << __func__ << ": unknown bam type, error!" << endl;
			exit(1);
		}

		for(size_t i=0; i<alnDataVector.size(); i++){
			b = alnDataVector.at(i);

			if(!(b->core.flag & BAM_FUNMAP)){ // aligned
				switch(bam_type){
					case BAM_CIGAR_NO_DIFF_MD:
						//alnSegs = generateAlnSegs(b);
						alnSegs = generateAlnSegs2(b, startPos, endPos);
						break;
					case BAM_CIGAR_NO_DIFF_NO_MD:
					case BAM_CIGAR_DIFF_MD:
					case BAM_CIGAR_DIFF_NO_MD:
						//alnSegs = generateAlnSegs_no_MD(b, refseq, startPos, endPos);
						alnSegs = generateAlnSegs_no_MD2(b, refseq, startPos, endPos);
						break;
					default:
						cerr << __func__ << ": unknown bam type, error!" << endl;
						exit(1);
				}// generate align segments

				queryGtSig = new queryGtSig_t();
				queryGtSig->group_id = -1;
				//queryGtSig->score = -1;
				queryGtSig->seed_flag = false;
				queryGtSig->queryname = bam_get_qname(b);
				queryGtSig->gtSig_vec = extractGtSigsFromAlnSegsSingleQuery(alnSegs, startPos, endPos, refseq, sig_size_thres, clip_size_thres);
				queryGtSig_vec.push_back(queryGtSig);

				destroyAlnSegs(alnSegs);
			}
		}
	}

	// assign the reg_contain_flag
	//assignRegContainFlag(queryGtSig_vec, reg);

	return queryGtSig_vec;
}

// extract CIGAR signatures for single query
vector<gtSig_t*> genotyping::extractGtSigsFromAlnSegsSingleQuery(vector<struct alnSeg*> &alnSegs, int64_t startPos, int64_t endPos, string &refseq, int32_t sig_size_thres, int32_t clip_size_thres){
	vector<gtSig_t*> sig_vec;
	gtSig_t *gt_sig;
	int64_t pos, epos, tmp_endPos, misbase;
	int32_t idx, endflag, position, len;
	vector<struct alnSeg*>::iterator seg;
	string str;

	for(seg=alnSegs.begin(); seg!=alnSegs.end(); seg++){
		switch((*seg)->opflag){
//			case BAM_CMATCH:
//				if(isOverlappedPos((*seg)->startRpos, (*seg)->startRpos+(*seg)->seglen-1, startPos, endPos)){
//					if((*seg)->seglen>=sig_size_thres){
//						gt_sig = allocateGtSigNode(*seg);
//						sig_vec.push_back(gt_sig);
//					}
//				}
//				break;
//			//case MD_MISMATCH:  // single base mismatch
//			case BAM_CDIFF:  // single base mismatch
//				if((*seg)->startRpos>=startPos and (*seg)->startRpos<=endPos){
//					if((*seg)->seglen>=sig_size_thres){
//						gt_sig = allocateGtSigNode(*seg);
//						sig_vec.push_back(gt_sig);
//					}
//				}
//				break;
			case BAM_CINS:  // insertion in query
				if((*seg)->startRpos>=startPos and (*seg)->startRpos<=endPos){
					if((*seg)->seglen>=sig_size_thres){
						gt_sig = allocateGtSigNode(*seg);
						gt_sig->chrname = reg->chrname;
						sig_vec.push_back(gt_sig);
					}
				}
				break;
			case BAM_CDEL:  // deletion in query
				if(isOverlappedPos((*seg)->startRpos, (*seg)->startRpos+(*seg)->seglen-1, startPos, endPos)){
					if((*seg)->seglen>=sig_size_thres){
						gt_sig = allocateGtSigNode(*seg);
						gt_sig->chrname = reg->chrname;
						sig_vec.push_back(gt_sig);
					}
				}
				break;
			case BAM_CSOFT_CLIP:  // soft clipping in query
			case BAM_CHARD_CLIP:  // hard clipping in query
				if((*seg)->startRpos>=startPos and (*seg)->startRpos<=endPos){
					if((*seg)->seglen>=clip_size_thres){
						gt_sig = allocateGtSigNode(*seg);
						gt_sig->chrname = reg->chrname;
						sig_vec.push_back(gt_sig);
					}
				}
				break;
			case BAM_CEQUAL:
//				if(isOverlappedPos((*seg)->startRpos, (*seg)->startRpos+(*seg)->seglen-1, startPos, endPos)){
//					if((*seg)->seglen>=sig_size_thres){
//						gt_sig = allocateGtSigNode(*seg);
//						sig_vec.push_back(gt_sig);
//					}
//				}
				break;
			case BAM_CDIFF:
				if((*seg)->startRpos>=startPos and (*seg)->startRpos<=endPos){
					if((*seg)->seglen>=sig_size_thres){
						gt_sig = allocateGtSigNode(*seg);
						gt_sig->chrname = reg->chrname;
						sig_vec.push_back(gt_sig);
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
				cerr << __func__ << ", line=" << __LINE__ << ": invalid opflag " << (*seg)->opflag << endl;
				exit(1);
				break;
		}
	}

	return sig_vec;
}

// allocate genotype signature
gtSig_t *genotyping::allocateGtSigNode(struct alnSeg *aln_seg){
	gtSig_t *gt_sig = new gtSig_t();
	gt_sig->cigar_op = aln_seg->opflag;
	gt_sig->cigar_op_len = aln_seg->seglen;
	gt_sig->reg_contain_flag = false;
	gt_sig->ref_pos = aln_seg->startRpos;
	gt_sig->sig_id = -1;
	gt_sig->chrname = "";
	gt_sig->mate_gtSig = NULL;
	return gt_sig;
}

// destroy genotype signatures for queries
void genotyping::destroyQueryGtSigVec(vector<queryGtSig_t*> &queryGtSig_vec){
	queryGtSig_t *queryGtSig;
	vector<gtSig_t*> sig_vec_singleQuery;
	gtSig_t *gt_sig;
	for(size_t i=0; i<queryGtSig_vec.size(); i++){
		queryGtSig = queryGtSig_vec.at(i);
		for(size_t j=0; j<queryGtSig->gtSig_vec.size(); j++){
			gt_sig = queryGtSig->gtSig_vec.at(j);
			delete gt_sig;
		}
		vector<gtSig_t*>().swap(queryGtSig->gtSig_vec);
		delete queryGtSig;
	}
	vector<queryGtSig_t*>().swap(queryGtSig_vec);
}

// assign the reg_contain_flag
void genotyping::assignRegContainFlag(vector<queryGtSig_t*> &queryGtSig_vec, reg_t *reg){
	queryGtSig_t *queryGtSig;
	vector<gtSig_t*> gtSig_vec;
	gtSig_t *gt_sig;
	int64_t left_most_pos, right_most_pos;

	left_most_pos = reg->startRefPos - clip_extend_match_thres;
	right_most_pos = reg->endRefPos + clip_extend_match_thres;
	if(left_most_pos<1) left_most_pos = 1;

	for(size_t i=0; i<queryGtSig_vec.size(); i++){
		queryGtSig = queryGtSig_vec.at(i);
		gtSig_vec = queryGtSig->gtSig_vec;

		for(size_t j=0; j<gtSig_vec.size(); j++){
			gt_sig = gtSig_vec.at(j);
			if(gt_sig->chrname.compare(reg->chrname)==0){
				switch(gt_sig->cigar_op){
					case BAM_CINS:  // insertion in query
					case BAM_CSOFT_CLIP:  // soft clipping in query
					case BAM_CHARD_CLIP:  // hard clipping in query
						if(gt_sig->ref_pos>=left_most_pos and gt_sig->ref_pos<=right_most_pos)
							gt_sig->reg_contain_flag = true;
						break;
					case BAM_CEQUAL:
						break;
					case BAM_CDIFF:
					case BAM_CDEL:  // deletion in query
						if(gt_sig->ref_pos>=left_most_pos and gt_sig->ref_pos+gt_sig->cigar_op_len-1<=right_most_pos)
							gt_sig->reg_contain_flag = true;
						break;
					case BAM_CREF_SKIP:  // unexpected events
					case BAM_CPAD:
					case BAM_CBACK:
						cerr << __func__ << ", line=" << __LINE__ << ": invalid opflag " << gt_sig->cigar_op << endl;
						exit(1);
						break;
					default:
						cerr << __func__ << ", line=" << __LINE__ << ": invalid opflag " << gt_sig->cigar_op << endl;
						exit(1);
						break;
				}
			}
		}
	}
}

// genotype signature matching
void genotyping::gtSigMatch(vector<queryGtSig_t*> &queryGtSig_vec){
	// choose seed query
	seed_gtQuery = chooseSeedGtQuery(queryGtSig_vec);

	//cout << "seed query: " << seed_gtQuery->queryname << ", queryGtSig_vec.size=" << queryGtSig_vec.size() << endl;

	// compute the match profile of signatures of each query
	computeGtMatchProfile(queryGtSig_vec, seed_gtQuery);
}

// choose seed query: select the most number of signatures
queryGtSig_t *genotyping::chooseSeedGtQuery(vector<queryGtSig_t*> &queryGtSig_vec){
	queryGtSig_t *seed_gtQuery = NULL, *gt_query;
	int32_t maxValue;

	maxValue = -1;
	for(size_t i=0; i<queryGtSig_vec.size(); i++){
		gt_query = queryGtSig_vec.at(i);
		if(maxValue<(int32_t)gt_query->gtSig_vec.size()){
			maxValue = gt_query->gtSig_vec.size();
			seed_gtQuery = gt_query;
		}
	}
	if(seed_gtQuery) seed_gtQuery->seed_flag = true;

	return seed_gtQuery;
}

// compute the match profile of signatures of each query
void genotyping::computeGtMatchProfile(vector<queryGtSig_t*> &queryGtSig_vec, queryGtSig_t *seed_gtQuery){
	queryGtSig_t *queryGtSig;
	if(seed_gtQuery){
		for(size_t i=0; i<queryGtSig_vec.size(); i++){
			queryGtSig = queryGtSig_vec.at(i);
			if(queryGtSig->gtSig_vec.size()==0){
				noProfileQueryNum++;
			}
		}
	}
	if(seed_gtQuery){
		for(size_t i=0; i<queryGtSig_vec.size(); i++){
			queryGtSig = queryGtSig_vec.at(i);
			if(queryGtSig!=seed_gtQuery and queryGtSig->gtSig_vec.size()>0){
				//cout << "[" << i << "]: queryname=" << queryGtSig->queryname << endl;
				queryGtSig->match_profile_vec = computeGtMatchProfileSingleQuery(queryGtSig, seed_gtQuery);
			}else if(queryGtSig->gtSig_vec.size()==0){
				//cout << "gtSig_vec.size=0, [" << i << "]: queryname=" << queryGtSig->queryname << endl;
			}
		}
	}
}

// compute the match score of signatures of each query
vector<int8_t> genotyping::computeGtMatchProfileSingleQuery(queryGtSig_t *queryGtSig, queryGtSig_t *seed_gtQuery){
	vector<int8_t> match_profile_vec;
	int32_t rowsNum, colsNum, matchScore, mismatchScore, gapScore, gapOpenScore;
	int32_t scoreIJ, tmp_gapScore1, tmp_gapScore2, maxValue, path_val, op_len_tmp;
	struct alnScoreNode *scoreArr;
	int64_t i, j, arrSize;
	bool matchFlag;

	matchScore = GT_SIG_MATCH_SCORE;
	mismatchScore = GT_SIG_MISMATCH_SCORE;
	gapScore = GT_SIG_GAP_SCORE;
	gapOpenScore = GT_SIG_GAP_OPEN_SCORE;

	rowsNum = queryGtSig->gtSig_vec.size() + 1;
	colsNum = seed_gtQuery->gtSig_vec.size() + 1;

	arrSize = rowsNum * colsNum;
	scoreArr = (struct alnScoreNode*) calloc (arrSize, sizeof(struct alnScoreNode));
	if(scoreArr==NULL){
		cerr << "line=" << __LINE__ << ", rowsNum=" << rowsNum << ", colsNum=" << colsNum << ", cannot allocate memory, error!" << endl;
		exit(1);
	}

	// set the elements of the first row and the first column to be zero
	scoreArr[0].path_val = 0;
	for(j=1; j<colsNum; j++) scoreArr[j].path_val = 2;
	for(i=1; i<rowsNum; i++) scoreArr[i*colsNum].path_val = 1;

	// compute the scores of each element
	for(i=1; i<rowsNum; i++){
		for(j=1; j<colsNum; j++){
			matchFlag = isGtSigMatch(queryGtSig->gtSig_vec.at(i-1), seed_gtQuery->gtSig_vec.at(j-1), size_ratio_match_thres);
			if(matchFlag) scoreIJ = matchScore;
			else scoreIJ = mismatchScore;

			if(scoreArr[(i-1)*colsNum+j].path_val!=1)
				tmp_gapScore1 = gapOpenScore;
			else
				tmp_gapScore1 = gapScore;

			if(scoreArr[i*colsNum+j-1].path_val!=2)
				tmp_gapScore2 = gapOpenScore;
			else
				tmp_gapScore2 = gapScore;

			maxValue = INT_MIN;
			path_val = -1;
			// compute the maximal score
			if(scoreArr[(i-1)*colsNum+j-1].score+scoreIJ>maxValue) {// from (i-1, j-1)
				maxValue = scoreArr[(i-1)*colsNum+j-1].score + scoreIJ;
				path_val = 0;
			}
			if(scoreArr[(i-1)*colsNum+j].score+tmp_gapScore1>maxValue) {// from (i-1, j)
				maxValue = scoreArr[(i-1)*colsNum+j].score + tmp_gapScore1;
				path_val = 1;
			}
			if(scoreArr[i*colsNum+j-1].score+tmp_gapScore2>maxValue) {// from (i, j-1)
				maxValue = scoreArr[i*colsNum+j-1].score + tmp_gapScore2;
				path_val = 2;
			}

			scoreArr[i*colsNum+j].score = maxValue;
			scoreArr[i*colsNum+j].path_val = path_val;
		}
	}
/*
	// print score array
	cout << "\t";
	for(j=1; j<colsNum; j++) {
		if(seed_gtQuery->gtSig_vec.at(j-1)->cigar_op==BAM_CSOFT_CLIP or seed_gtQuery->gtSig_vec.at(j-1)->cigar_op==BAM_CHARD_CLIP)
			op_len_tmp = seed_gtQuery->gtSig_vec.at(j-1)->ref_pos;
		else op_len_tmp = seed_gtQuery->gtSig_vec.at(j-1)->cigar_op_len;
		cout << "\t" + to_string(seed_gtQuery->gtSig_vec.at(j-1)->cigar_op) + "_" + to_string(op_len_tmp);
	}
	cout << endl;
	for(j=0; j<colsNum; j++) cout << "\t" << scoreArr[j].score;
	cout << endl;

	for(i=1; i<rowsNum; i++){
		if(queryGtSig->gtSig_vec.at(i-1)->cigar_op==BAM_CSOFT_CLIP or queryGtSig->gtSig_vec.at(i-1)->cigar_op==BAM_CHARD_CLIP)
			op_len_tmp = queryGtSig->gtSig_vec.at(i-1)->ref_pos;
		else op_len_tmp = queryGtSig->gtSig_vec.at(i-1)->cigar_op_len;
		cout << to_string(queryGtSig->gtSig_vec.at(i-1)->cigar_op) + "_" + to_string(op_len_tmp);
		for(j=0; j<colsNum; j++) cout << "\t" << scoreArr[i*colsNum+j].score;
		cout << endl;
	}
	cout << "score finished." << endl;
*/

	// compute signature match profile
	match_profile_vec = computeSigMatchProfile(scoreArr, rowsNum, colsNum, queryGtSig, seed_gtQuery);

/*
	//print profile
	cout << "profile start from here: " << endl;
	if((noProfileQueryNum>=queryGtSig_vec.size()*min_alle_ratio_thres) and (noProfileQueryNum<=queryGtSig_vec.size()*max_alle_ratio_thres)){
		for(i=0; i<(int64_t)match_profile_vec.size(); i++) cout << to_string(match_profile_vec.at(i)) << "\t";
		cout << endl;
	}
*/
	return match_profile_vec;
}

// determine whether two genotype signatures are matched
bool genotyping::isGtSigMatch(gtSig_t *gt_sig, gtSig_t *seed_gt_sig, double size_ratio_match_thres){
	bool match_flag = false, flag;
	double size_ratio;
	int64_t start_pos1, end_pos1, start_pos2, end_pos2;

	if(gt_sig and seed_gt_sig){
		if(gt_sig->chrname.compare(seed_gt_sig->chrname)==0){
			if((gt_sig->cigar_op==BAM_CSOFT_CLIP or gt_sig->cigar_op==BAM_CHARD_CLIP) and (seed_gt_sig->cigar_op==BAM_CSOFT_CLIP or seed_gt_sig->cigar_op==BAM_CHARD_CLIP)){
				// clipping, and consider the breakpoint location
				//cout << "gt_sig.cigar_op=" << gt_sig->cigar_op << ", pos=" << gt_sig->ref_pos << ", cigar_op_len=" << gt_sig->cigar_op_len << "; seed_gt_sig.cigar_op=" << seed_gt_sig->cigar_op << ", pos=" << seed_gt_sig->ref_pos << ", cigar_op_len=" << seed_gt_sig->cigar_op_len << endl;

				start_pos1 = gt_sig->ref_pos - clip_extend_match_thres;
				end_pos1 = gt_sig->ref_pos + clip_extend_match_thres;
				if(start_pos1<1) start_pos1 = 1;

				start_pos2 = seed_gt_sig->ref_pos - clip_extend_match_thres;
				end_pos2 = seed_gt_sig->ref_pos + clip_extend_match_thres;
				if(start_pos2<1) start_pos2 = 1;

				flag = isOverlappedPos(start_pos1, end_pos1, start_pos2, end_pos2);
				if(flag) match_flag = true;

			}else{ // same cigar_op, ratio of signature size less than 0.7
				if(gt_sig->cigar_op==seed_gt_sig->cigar_op){
					if(gt_sig->cigar_op_len>0 and seed_gt_sig->cigar_op_len>0) {
						if(gt_sig->cigar_op_len <= seed_gt_sig->cigar_op_len) size_ratio = (double)gt_sig->cigar_op_len / seed_gt_sig->cigar_op_len;
						else size_ratio = (double)seed_gt_sig->cigar_op_len / gt_sig->cigar_op_len;
						if(size_ratio>=size_ratio_match_thres)
							match_flag = true;
					}else match_flag = false;
				}
			}
		}
	}

	return match_flag;
}

// compute signature match profile
vector<int8_t> genotyping::computeSigMatchProfile(struct alnScoreNode *scoreArr, int32_t rowsNum, int32_t colsNum, queryGtSig_t *queryGtSig, queryGtSig_t *seed_gtQuery){
	vector<int8_t> match_profile_vec;
	int32_t i = rowsNum - 1, j = colsNum - 1, value;

	while(i>=0 and j>=0){
		value = scoreArr[i*colsNum+j].path_val;
		if(value==0){ //from (i-1, j-1)
			match_profile_vec.push_back(1);
			queryGtSig->gtSig_vec.at(i-1)->mate_gtSig = seed_gtQuery->gtSig_vec.at(j-1);
			i--;
			j--;
		}
		else if(value==1){ //from (i-1, j)
			i--;
		}
		else{ //from (i, j-1)
			match_profile_vec.push_back(0);
			j--;
		}
		if(i==0 and j==0)
			break;
	}

	reverse(match_profile_vec.begin(), match_profile_vec.end());

	return match_profile_vec;
}

// compute genotype string
vector<string> genotyping::computeGenotypeString(vector<queryGtSig_t*> &queryGtSig_vec, vector<bool> &validGtSigFlagVec){
	vector<string> gt_dp_ad_str;
	string gt_str;
	size_t i, j, k;
	double noProfileQueryRatio;
	profile_pat_t *profile_pat, *profile_pat_tmp;
	queryGtSig_t *queryGtSig;
	bool exist_flag, match_flag;
	int64_t maxNum, secNum, items_have_profile, valid_colNum;
	string dp_str, ad_str1, ad_str2;
	double maxRatio, secRatio, ratio;

	gt_str = "./.";
	dp_str = ".";
	ad_str1 = ".";
	ad_str2 = ".";

	valid_colNum = 0;
	for(i=0; i<validGtSigFlagVec.size(); i++) if(validGtSigFlagVec.at(i)==true) valid_colNum ++;

	if(queryGtSig_vec.size()>0 and valid_colNum>0){
		noProfileQueryRatio = (double)noProfileQueryNum / queryGtSig_vec.size();
		if(noProfileQueryRatio>0.75) gt_str = "./.";
		else{
			items_have_profile = 0;
			for(i=0; i<queryGtSig_vec.size(); i++){
				exist_flag = false;
				queryGtSig = queryGtSig_vec.at(i);
				if(queryGtSig->match_profile_vec.size()){
					items_have_profile ++;
					for(j=0; j<match_profile_pat_vec.size(); j++){
						match_flag = true;
						profile_pat = match_profile_pat_vec.at(j);
						for(k=0; k<profile_pat->match_profile_vec.size(); k++){
							if(validGtSigFlagVec.at(k)==true and queryGtSig->match_profile_vec.at(k)!=profile_pat->match_profile_vec.at(k)){
								match_flag = false;
								break;
							}
						}
						if(match_flag) {
							profile_pat->num ++;
							exist_flag = true;
							break;
						}
					}
					if(exist_flag==false){ // create new node
						profile_pat_tmp = new struct profile_pat_node();
						profile_pat_tmp->num = 1;
						for(k=0; k<queryGtSig->match_profile_vec.size(); k++){
							if(validGtSigFlagVec.at(k)==true) profile_pat_tmp->match_profile_vec.push_back(queryGtSig->match_profile_vec.at(k));
						}
						match_profile_pat_vec.push_back(profile_pat_tmp);
					}
				}
			}

			// bubble sort
			if(match_profile_pat_vec.size()>=2){
				for(i=0; i<match_profile_pat_vec.size(); i++){
					for(j=0; j<match_profile_pat_vec.size()-i-1; j++){
						if(match_profile_pat_vec.at(j)->num<match_profile_pat_vec.at(j+1)->num){ // swap
							profile_pat_tmp = match_profile_pat_vec.at(j);
							match_profile_pat_vec.at(j) = match_profile_pat_vec.at(j+1);
							match_profile_pat_vec.at(j+1) = profile_pat_tmp;
						}
					}
				}
			}
/*
			// print profile information
			cout << "match_profile_pat_vec.size=" << match_profile_pat_vec.size() << endl;
			for(i=0; i<match_profile_pat_vec.size(); i++){
				profile_pat_tmp = match_profile_pat_vec.at(i);
				cout << "num=" << profile_pat_tmp->num;
				for(j=0; j<profile_pat_tmp->match_profile_vec.size(); j++){
					cout << "\t" << to_string(profile_pat_tmp->match_profile_vec.at(j));
				}
				cout << endl;
			}
*/
			// get the maximum and second maximum
			maxNum = secNum = 0;
			maxRatio = secRatio = 0;
			if(items_have_profile>0){
				if(match_profile_pat_vec.size()==1) maxNum = match_profile_pat_vec.at(0)->num;
				else if(match_profile_pat_vec.size()>1) {
					maxNum = match_profile_pat_vec.at(0)->num;
					secNum = match_profile_pat_vec.at(1)->num;
				}
				maxRatio = (double)maxNum / items_have_profile;
				secRatio = (double)secNum / items_have_profile;
			}

			if(noProfileQueryRatio>=0.25 and noProfileQueryRatio<=0.75){
				if(maxRatio>0.7) gt_str = "0/1";
				else{
					cout << "line= " << __LINE__ << ", noProfileQueryRatio=" << noProfileQueryRatio << ", maxNum=" << maxNum << ", secNum=" << secNum << ", maxRatio=" << maxRatio << ", secRatio=" << secRatio << endl;
				}
			}else if(noProfileQueryRatio<0.25){
				ratio = (double)secNum/maxNum;
				if(ratio<min_alle_ratio_thres) gt_str = "1/1";
				else if(ratio>=min_alle_ratio_thres and ratio<=max_alle_ratio_thres){
					gt_str = "0/1";
				}else{
					cout << "line= " << __LINE__ << ", noProfileQueryRatio=" << noProfileQueryRatio << ", ratio=" << ratio << ", maxNum=" << maxNum << ", secNum=" << secNum << ", maxRatio=" << maxRatio << ", secRatio=" << secRatio << endl;
				}
			}else{
				cout << "line= " << __LINE__ << ", noProfileQueryRatio=" << noProfileQueryRatio << ", maxNum=" << maxNum << ", secNum=" << secNum << ", maxRatio=" << maxRatio << ", secRatio=" << secRatio << endl;
			}
		}
	}

	if(gt_str.compare("./.")!=0){
		dp_str = to_string(queryGtSig_vec.size());
		ad_str1 = to_string(noProfileQueryNum);
		ad_str2 = to_string(queryGtSig_vec.size() - noProfileQueryNum);
	}else{
		dp_str = ".";
		ad_str1 = ".";
		ad_str2 = ".";
	}

	gt_dp_ad_str.push_back("GT:DP:AD");
	gt_dp_ad_str.push_back(gt_str + ":" + dp_str + ":" + ad_str1 + "," + ad_str2);

	return gt_dp_ad_str;
}

// compute valid gtSig columns
vector<bool> genotyping::computeValidGtSigFlags(queryGtSig_t *seed_gtQuery, reg_t *reg){
	struct loc_dif_node{
		int32_t col_id, loc_dif;
	};

	vector<bool> valid_gtSig_flag_vec;
	size_t i, j;
	gtSig_t *gtSig, *seed_gt_sig;
	int64_t end_pos, num, reg_size, gtSig_size, loc_dif;
	bool overlap_flag, exist_overlap_flag = false, type_match_flag, size_match_flag;
	double ratio;
	struct loc_dif_node *loc_dif_item;
	vector<struct loc_dif_node*> loc_dif_vec;

	if(seed_gtQuery){
		//cout << "seed_gtQuery->gtSig_vec.size() = " << seed_gtQuery->gtSig_vec.size() << endl;
		for(i=0; i<seed_gtQuery->gtSig_vec.size(); i++){
			gtSig = seed_gtQuery->gtSig_vec.at(i);
			end_pos = getSigEndPos(gtSig);
			overlap_flag = isOverlappedPos(gtSig->ref_pos, end_pos,reg->startRefPos, reg->endRefPos);
			if(overlap_flag){
				valid_gtSig_flag_vec.push_back(true);
				exist_overlap_flag = true;
			}else{
				valid_gtSig_flag_vec.push_back(false);
			}
			//cout << "valid_gtSig_flag_vec[i] = " << valid_gtSig_flag_vec.at(i) << endl;
		}

		if(exist_overlap_flag == false){

			for(i=0; i<seed_gtQuery->gtSig_vec.size(); i++){
				gtSig = seed_gtQuery->gtSig_vec.at(i);
				end_pos = getSigEndPos(gtSig);
				if(gtSig->ref_pos<reg->startRefPos) // gtSig is on the left side
					loc_dif = reg->startRefPos - end_pos;
				else  // gtSig is on the right side
					loc_dif = gtSig->ref_pos - reg->endRefPos;

				loc_dif_item = new struct loc_dif_node();
				loc_dif_item->col_id = i;
				loc_dif_item->loc_dif = loc_dif;
				loc_dif_vec.push_back(loc_dif_item);
			}

			// bubble sort according to 'loc_dif' field
			if(loc_dif_vec.size()>=2){
				for(i=0; i<loc_dif_vec.size(); i++){
					for(j=0; j<loc_dif_vec.size()-i-1; j++){
						if(loc_dif_vec.at(j)->loc_dif>loc_dif_vec.at(j+1)->loc_dif){ // swap
							loc_dif_item = loc_dif_vec.at(j);
							loc_dif_vec.at(j) = loc_dif_vec.at(j+1);
							loc_dif_vec.at(j+1) = loc_dif_item;
						}
					}
				}
			}

			// check each gtSig one by one according to distance
			for(i=0; i<loc_dif_vec.size(); i++){
				loc_dif_item = loc_dif_vec.at(i);
				gtSig = seed_gtQuery->gtSig_vec.at(loc_dif_item->col_id);
				end_pos = getSigEndPos(gtSig);

				// type condition
				type_match_flag = true;
				if(reg->var_type==VAR_INS){ // insertion
					if(gtSig->cigar_op!=BAM_CINS) type_match_flag = false;
				}else if(reg->var_type==VAR_DEL){ // deletion
					if(gtSig->cigar_op!=BAM_CDEL) type_match_flag = false;
				}

				// size condition
				size_match_flag = false;
				reg_size = reg->endRefPos - reg->startRefPos + 1;
				gtSig_size = end_pos - gtSig->ref_pos + 1;
				if(reg_size<gtSig_size) ratio = (double)reg_size / gtSig_size;
				else ratio = (double)gtSig_size / reg_size;
				if(ratio>=size_ratio_match_thres) size_match_flag = true;

				// choose the one with the nearest distance
				if(type_match_flag and size_match_flag){ // found and choose it
					valid_gtSig_flag_vec.at(loc_dif_item->col_id) = true;
					break;
				}
			}

			// free memory
			for(i=0; i<loc_dif_vec.size(); i++){
				loc_dif_item = loc_dif_vec.at(i);
				delete loc_dif_item;
			}
			vector<struct loc_dif_node*>().swap(loc_dif_vec);
		}
	}

	return valid_gtSig_flag_vec;
}

// recover variants
void genotyping::recoverVariants(reg_t *reg, vector<profile_pat_t*> &match_profile_pat_vec, queryGtSig_t *seed_gtQuery, vector<queryGtSig_t*> &queryGtSig_vec, vector<bool> &validGtSigFlagVec){
	int32_t pat_sig_idx, sig_col_id;
	gtSig_t *target_sig;
	vector<int64_t> aver_pos_vec;

	//get the target signature
	pat_sig_idx = getTargetPatSigIDFromProfile(match_profile_pat_vec);
	if(pat_sig_idx>=0){
		//compute average variant region
		sig_col_id = computeSigColID(pat_sig_idx, validGtSigFlagVec);
		target_sig = seed_gtQuery->gtSig_vec.at(sig_col_id);
		aver_pos_vec = computeAverVarReg(target_sig, seed_gtQuery, queryGtSig_vec);
		//compare and make decision
		compareAndUpdateVarReg(reg, target_sig, aver_pos_vec);
	}
}

// get target pattern signature
int32_t genotyping::getTargetPatSigIDFromProfile(vector<profile_pat_t*> &match_profile_pat_vec){
	profile_pat_t *profile_pat;
	size_t i, j, arr_size;
	int32_t sig_id = -1, max_idx = -1, maxValue;

	if(match_profile_pat_vec.size()>0){
		arr_size = match_profile_pat_vec.at(0)->match_profile_vec.size();
		int32_t num_arr[arr_size];
		for(i=0; i<arr_size; i++) num_arr[i] = 0;

		for(i=0; i<match_profile_pat_vec.size(); i++){
			profile_pat = match_profile_pat_vec.at(i);
			for(j=0; j<profile_pat->match_profile_vec.size(); j++){
				if(profile_pat->match_profile_vec.at(j)==1){
					num_arr[j] += profile_pat->num;
				}
			}
		}

//		cout << "num_arr: " << num_arr[0];
//		for(i=1; i<arr_size; i++) cout << "\t" << num_arr[i];
//		cout << endl;

		// select the maximum
		max_idx = -1;
		maxValue = 0;
		for(i=0; i<arr_size; i++){
			if(maxValue<num_arr[i]){
				max_idx = i;
				maxValue = num_arr[i];
			}
		}
		sig_id = max_idx;
	}

	//cout << "sig_id=" << sig_id << ", max_idx=" << max_idx << endl;

	return sig_id;
}

// compute gtSig column index
int32_t genotyping::computeSigColID(int32_t pat_sig_idx, vector<bool> &validGtSigFlagVec){
	int32_t sig_col_id = -1, pat_sig_id_tmp;

	pat_sig_id_tmp = 0;
	for(size_t i=0; i<validGtSigFlagVec.size(); i++){
		if(validGtSigFlagVec.at(i)==true){
			if(pat_sig_id_tmp==pat_sig_idx){
				sig_col_id = i;
				break;
			}
			pat_sig_id_tmp ++;
		}
	}

	return sig_col_id;
}

//compute average variant region
vector<int64_t> genotyping::computeAverVarReg(gtSig_t *target_sig, queryGtSig_t *seed_gtQuery, vector<queryGtSig_t*> &queryGtSig_vec){
	vector<int64_t> reg_ret;
	size_t i, j;
	int64_t sum_startpos = 0, sum_endpos = 0, num = 0, sig_len = 0, sum_sig_len = 0, end_pos, aver_startpos, aver_endpos, aver_siglen;
	queryGtSig_t *queryGtSig;
	gtSig_t *gt_sig;

	for(i=0; i<queryGtSig_vec.size(); i++){
		queryGtSig = queryGtSig_vec.at(i);
		if(queryGtSig!=seed_gtQuery){ // not the seed query
			for(j=0; j<queryGtSig->gtSig_vec.size(); j++){
				gt_sig = queryGtSig->gtSig_vec.at(j);
				if(gt_sig->mate_gtSig==target_sig){ // found the mated signature
					sum_startpos += gt_sig->ref_pos;
					end_pos = getSigEndPos(gt_sig);
					if(end_pos>0) sum_endpos += end_pos;
					else{
						cerr << __func__ << ", line=" << __LINE__ << ": invalid end_pos=" << end_pos << endl;
						exit(1);
					}
					sig_len = getSigLen(gt_sig);
					sum_sig_len += sig_len;
					num ++;
				}
			}
		}else{ // seed query
			sum_startpos += target_sig->ref_pos;
			end_pos = getSigEndPos(target_sig);
			if(end_pos>0) sum_endpos += end_pos;
			else{
				cerr << __func__ << ", line=" << __LINE__ << ": invalid end_pos=" << end_pos << endl;
				exit(1);
			}
			sig_len = getSigLen(target_sig);
			sum_sig_len += sig_len;
			num ++;
		}
	}

	// compute the average value
	if(num>0){
		aver_startpos = sum_startpos / num;
		aver_endpos = sum_endpos / num;
		aver_siglen = sum_sig_len / num;
	}else{
		aver_startpos = aver_endpos = aver_siglen = 0;
	}
	reg_ret.push_back(aver_startpos);
	reg_ret.push_back(aver_endpos);
	reg_ret.push_back(aver_siglen);

	//cout << "cigar_op=" << target_sig->cigar_op << ", num=" << num << ", average region: [" << aver_startpos << "," << aver_endpos << "], aver_siglen=" << aver_siglen << endl;

	return reg_ret;
}

// get the endpos
int64_t genotyping::getSigEndPos(gtSig_t *gt_sig){
	int64_t end_pos = 0;

	switch(gt_sig->cigar_op){
		case BAM_CINS:  // insertion in query
			end_pos = gt_sig->ref_pos;
			break;
		case BAM_CDEL:  // deletion in query
			end_pos += gt_sig->ref_pos + gt_sig->cigar_op_len - 1;
			break;
		case BAM_CEQUAL:
			break;
		case BAM_CDIFF:
			end_pos = gt_sig->ref_pos;
			break;
		case BAM_CSOFT_CLIP:  // soft clipping in query
		case BAM_CHARD_CLIP:  // hard clipping in query
			end_pos = gt_sig->ref_pos;
			break;
		case BAM_CREF_SKIP:  // unexpected events
		case BAM_CPAD:
		case BAM_CBACK:
			cerr << __func__ << ", line=" << __LINE__ << ": invalid opflag " << gt_sig->cigar_op << endl;
			exit(1);
			break;
		default:
			cerr << __func__ << ", line=" << __LINE__ << ": invalid opflag " << gt_sig->cigar_op << endl;
			exit(1);
			break;
	}

	return end_pos;
}

// get siglen
int32_t genotyping::getSigLen(gtSig_t *gt_sig){
	int32_t sig_len;

	switch(gt_sig->cigar_op){
		case BAM_CINS:  // insertion in query
		case BAM_CDEL:  // deletion in query
		case BAM_CDIFF:
			sig_len = gt_sig->cigar_op_len;
			break;
		case BAM_CEQUAL:
			break;
		case BAM_CSOFT_CLIP:  // soft clipping in query
		case BAM_CHARD_CLIP:  // hard clipping in query
			sig_len = 0;
			break;
		case BAM_CREF_SKIP:  // unexpected events
		case BAM_CPAD:
		case BAM_CBACK:
			cerr << __func__ << ", line=" << __LINE__ << ": invalid opflag " << gt_sig->cigar_op << endl;
			exit(1);
			break;
		default:
			cerr << __func__ << ", line=" << __LINE__ << ": invalid opflag " << gt_sig->cigar_op << endl;
			exit(1);
			break;
	}

	return sig_len;
}

//compare and make decision
void genotyping::compareAndUpdateVarReg(reg_t *reg, gtSig_t *target_sig, vector<int64_t> &aver_pos_vec){
	int64_t start_pos_aver, end_pos_aver, dist_aver, sig_len;
	int64_t start_pos_reg_tmp, end_pos_reg_tmp, start_pos_aver_tmp, end_pos_aver_tmp, reg_size, aver_reg_size;
	bool valid_flag = false, overlap_flag, var_type_flag = false;
	double ratio;

	start_pos_aver = aver_pos_vec.at(0);
	end_pos_aver = aver_pos_vec.at(1);
	sig_len = aver_pos_vec.at(2);
	dist_aver = end_pos_aver - start_pos_aver + 1;

	start_pos_reg_tmp = reg->startRefPos - clip_extend_match_thres;
	end_pos_reg_tmp = reg->endRefPos + clip_extend_match_thres;
	if(start_pos_reg_tmp<1) start_pos_reg_tmp = 1;

	start_pos_aver_tmp = start_pos_aver - clip_extend_match_thres;
	end_pos_aver_tmp = end_pos_aver + clip_extend_match_thres;
	if(start_pos_aver_tmp<1) start_pos_aver_tmp = 1;

	if(reg->call_success_status==false){ // recovery

		// check the overlap
		overlap_flag = isOverlappedPos(start_pos_reg_tmp, end_pos_reg_tmp, start_pos_aver_tmp, end_pos_aver_tmp);

		// check the region size
		if(target_sig->cigar_op!=BAM_CSOFT_CLIP and target_sig->cigar_op!=BAM_CHARD_CLIP){ // non-clipping
			reg_size = reg->endRefPos - reg->startRefPos + 1;
			aver_reg_size = end_pos_aver - start_pos_aver + 1;
			if(reg_size<aver_reg_size) ratio = (double)reg_size / aver_reg_size;
			else ratio = (double)aver_reg_size / reg_size;

			if(overlap_flag or ratio>=size_ratio_match_thres) valid_flag = true;
		}else{ // clipping
			if(overlap_flag) valid_flag = true;
		}

		if(valid_flag){
			// var_type
			switch(target_sig->cigar_op){
				case BAM_CINS:  // insertion in query
					reg->var_type = VAR_INS;
					//reg->sv_len = sig_len;
					break;
				case BAM_CDEL:  // deletion in query
					reg->var_type = VAR_DEL;
					//reg->sv_len = -dist_aver;
					break;
				case BAM_CSOFT_CLIP:  // soft clipping in query
				case BAM_CHARD_CLIP:  // hard clipping in query
					reg->var_type = VAR_BND;
					break;
				case BAM_CDIFF:
					cerr << __func__ << ", line=" << __LINE__ << ": invalid opflag: " << target_sig->cigar_op << endl;
					exit(1);
					break;
				default:
					cerr << __func__ << ", line=" << __LINE__ << ": invalid opflag: " << target_sig->cigar_op << endl;
					exit(1);
					break;
			}
			reg->call_success_status = true;
			//reg->short_sv_flag = false;
			reg->query_pos_invalid_flag = true;
		}
	}else{ // call success, then update

//		// check the overlap
//		overlap_flag = isOverlappedPos(start_pos_reg_tmp, end_pos_reg_tmp, start_pos_aver_tmp, end_pos_aver_tmp);

		// check the region size
		if (target_sig->cigar_op != BAM_CSOFT_CLIP and target_sig->cigar_op != BAM_CHARD_CLIP) { // non-clipping
			reg_size = reg->endRefPos - reg->startRefPos + 1;
			aver_reg_size = end_pos_aver - start_pos_aver + 1;
			if (reg_size < aver_reg_size)
				ratio = (double) reg_size / aver_reg_size;
			else
				ratio = (double) aver_reg_size / reg_size;

			if (ratio < size_ratio_match_thres) { //small ratio, then update pos and sv_len
				//reg->startRefPos = start_pos_aver;
				//reg->endRefPos = end_pos_aver;
				switch (target_sig->cigar_op) { // update sv_len
					case BAM_CINS:  // insertion in query
						//reg->sv_len = sig_len;
						break;
					case BAM_CDEL:  // deletion in query
						//reg->sv_len = -dist_aver;
						break;
					case BAM_CSOFT_CLIP:  // soft clipping in query
					case BAM_CHARD_CLIP:  // hard clipping in query
						break;
					case BAM_CDIFF:
						cerr << __func__ << ", line=" << __LINE__ << ": invalid opflag: " << target_sig->cigar_op << endl;
						exit(1);
						break;
					default:
						cerr << __func__ << ", line=" << __LINE__ << ": invalid opflag: " << target_sig->cigar_op << endl;
						exit(1);
						break;
				}
				reg->query_pos_invalid_flag = true;
			}
		}
	}
}
