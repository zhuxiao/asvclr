#include "alnDataLoader.h"
#include "clipAlnDataLoader.h"
#include "util.h"

//pthread_mutex_t mutex_down_sample = PTHREAD_MUTEX_INITIALIZER;

clipAlnDataLoader::clipAlnDataLoader(string &chrname, int64_t startRefPos, int64_t endRefPos, string &inBamFile, int32_t minClipEndSize, int32_t minMapQ, int32_t minHighMapQ) {
	this->chrname = chrname;
	this->startRefPos = startRefPos;
	this->endRefPos = endRefPos;
	this->inBamFile = inBamFile;
	this->minClipEndSize = minClipEndSize;
	this->minMapQ = minMapQ;
	this->minHighMapQ = minHighMapQ;
}

clipAlnDataLoader::~clipAlnDataLoader() {
}

// load clipping align data without down-sampling
void clipAlnDataLoader::loadClipAlnData(vector<clipAlnData_t*> &clipAlnDataVector){
	loadClipAlnData(clipAlnDataVector, 0);
}

void clipAlnDataLoader::loadClipAlnData(vector<clipAlnData_t*> &clipAlnDataVector, double max_ultra_high_cov){
	size_t i;
	vector<bam1_t*> alnDataVector;
	clipAlnData_t *clip_aln;
	bam_hdr_t *header;

	// load the align data
	alnDataLoader data_loader(chrname, startRefPos, endRefPos, inBamFile, minMapQ, minHighMapQ);
	data_loader.loadAlnData(alnDataVector, max_ultra_high_cov);

//	if(max_ultra_high_cov>0){
//		samplingAlnData(alnDataVector, data_loader.mean_read_len, max_ultra_high_cov);
//	}

	// load the sam/bam header
	header = loadSamHeader(inBamFile);

	// compute the aligned region
	for(i=0; i<alnDataVector.size(); i++){
		clip_aln = generateClipAlnData(alnDataVector.at(i), header);
		if(clip_aln){
			clipAlnDataVector.push_back(clip_aln);
		}else{
			cerr << __func__ << ", line=" << __LINE__ << ": cannot generate clip align item, error!" << endl;
			exit(1);
		}
	}

	bam_hdr_destroy(header);
}

void clipAlnDataLoader::loadClipAlnData(vector<clipAlnData_t*> &clipAlnDataVector, double max_ultra_high_cov, vector<string> &qname_vec){
	size_t i;
	vector<bam1_t*> alnDataVector;
	clipAlnData_t *clip_aln;
	bam_hdr_t *header;

	// load the align data
//	cout << __func__ << ", line=" << __LINE__ << "minMapQ :  " << minMapQ << endl;

	alnDataLoader data_loader(chrname, startRefPos, endRefPos, inBamFile, minMapQ, minHighMapQ);
//	data_loader.loadAlnData(alnDataVector);
	data_loader.loadAlnData(alnDataVector, max_ultra_high_cov, qname_vec);

//	if(max_ultra_high_cov>0){
//		samplingAlnData(alnDataVector, data_loader.mean_read_len, max_ultra_high_cov);
//	}

	// load the sam/bam header
	header = loadSamHeader(inBamFile);

	// compute the aligned region
	for(i=0; i<alnDataVector.size(); i++){
		clip_aln = generateClipAlnData(alnDataVector.at(i), header);
		if(clip_aln){
			clipAlnDataVector.push_back(clip_aln);
		}else{
			cerr << __func__ << ", line=" << __LINE__ << ": cannot generate clip align item, error!" << endl;
			exit(1);
		}
	}

	bam_hdr_destroy(header);
}

void clipAlnDataLoader::loadClipAlnData(vector<clipAlnData_t*> &clipAlnDataVector, vector<string> &qname_vec){
	size_t i;
	vector<bam1_t*> alnDataVector;
	clipAlnData_t *clip_aln;
	bam_hdr_t *header;

	// load the align data
	alnDataLoader data_loader(chrname, startRefPos, endRefPos, inBamFile, minMapQ, minHighMapQ);
	data_loader.loadAlnData(alnDataVector, qname_vec);

	// load the sam/bam header
	header = loadSamHeader(inBamFile);

	// compute the aligned region
	for(i=0; i<alnDataVector.size(); i++){
		clip_aln = generateClipAlnData(alnDataVector.at(i), header);
		if(clip_aln){
			clipAlnDataVector.push_back(clip_aln);
		}else{
			cerr << __func__ << ", line=" << __LINE__ << ": cannot generate clip align item, error!" << endl;
			exit(1);
		}
	}

	bam_hdr_destroy(header);
}

void clipAlnDataLoader::loadClipAlnDataWithSATag(vector<clipAlnData_t*> &clipAlnDataVector){
	loadClipAlnData(clipAlnDataVector);
	fillClipAlnDataBySATag(clipAlnDataVector);
	addAdjacentInfo(clipAlnDataVector); // order clipping segments
}

void clipAlnDataLoader::loadClipAlnDataWithSATag(vector<clipAlnData_t*> &clipAlnDataVector, double max_ultra_high_cov){
	loadClipAlnData(clipAlnDataVector, max_ultra_high_cov);
	fillClipAlnDataBySATag(clipAlnDataVector);
	addAdjacentInfo(clipAlnDataVector); // order clipping segments
}

void clipAlnDataLoader::loadClipAlnDataWithSATagWithSegSize(vector<clipAlnData_t*> &clipAlnDataVector, double max_ultra_high_cov, double primary_seg_size_ratio, double primary_seg_nm_ratio){
	loadClipAlnDataWithSATag(clipAlnDataVector, max_ultra_high_cov);
	removeClipAlnDataWithLowPrimarySegSizeRatio(clipAlnDataVector, primary_seg_size_ratio, primary_seg_nm_ratio);
}

void clipAlnDataLoader::loadClipAlnDataWithSATag(vector<clipAlnData_t*> &clipAlnDataVector, double max_ultra_high_cov, vector<string> &qname_vec){
	loadClipAlnData(clipAlnDataVector, max_ultra_high_cov, qname_vec);
	fillClipAlnDataBySATag(clipAlnDataVector);
	addAdjacentInfo(clipAlnDataVector);
}

void clipAlnDataLoader::loadClipAlnDataWithSATag(vector<clipAlnData_t*> &clipAlnDataVector, vector<string> &qname_vec){
	loadClipAlnData(clipAlnDataVector, qname_vec);
	fillClipAlnDataBySATag(clipAlnDataVector);
	addAdjacentInfo(clipAlnDataVector);
}

clipAlnData_t* clipAlnDataLoader::generateClipAlnData(bam1_t* b, bam_hdr_t *header){
	clipAlnData_t *clip_aln = NULL;
	uint32_t *c, op1, op2;
	uint8_t *nm_ptr;

	clip_aln = new clipAlnData_t();
	clip_aln->bam = b;
	clip_aln->queryname = bam_get_qname(b);
	clip_aln->chrname = header->target_name[b->core.tid];
	clip_aln->startRefPos = b->core.pos + 1;
	clip_aln->endRefPos = bam_endpos(b);
	clip_aln->left_aln = clip_aln->right_aln = NULL;
	clip_aln->leftmost_flag = clip_aln->rightmost_flag = false;
	clip_aln->NM = -1;

	c = bam_get_cigar(b);  // CIGAR
	// left clip
	clip_aln->leftHardClippedFlag = false;
	op1 = bam_cigar_op(c[0]);
	if(op1==BAM_CSOFT_CLIP or op1==BAM_CHARD_CLIP){
		clip_aln->leftClipSize = bam_cigar_oplen(c[0]);
		if(op1==BAM_CHARD_CLIP) clip_aln->leftHardClippedFlag = true;
	}else
		clip_aln->leftClipSize = 0;
	// right clip
	clip_aln->rightHardClippedFlag = false;
	op2 = bam_cigar_op(c[b->core.n_cigar-1]);
	if(op2==BAM_CSOFT_CLIP or op2==BAM_CHARD_CLIP){
		clip_aln->rightClipSize = bam_cigar_oplen(c[b->core.n_cigar-1]);
		if(op2==BAM_CHARD_CLIP) clip_aln->rightHardClippedFlag = true;
	}else
		clip_aln->rightClipSize = 0;

	// NM
	nm_ptr = bam_aux_get(b, "NM");
	if(nm_ptr){
		clip_aln->NM = bam_aux2i(nm_ptr);
		// cout << "clip_aln->startRefPos=" << clip_aln->startRefPos << ",clip_aln->queryname=" << clip_aln->queryname << ",clip_aln->NM=" << clip_aln->NM << endl;
	}

	// querylen
	clip_aln->querylen = b->core.l_qseq;
	if(op1==BAM_CHARD_CLIP)
		clip_aln->querylen += clip_aln->leftClipSize;
	if(op2==BAM_CHARD_CLIP)
		clip_aln->querylen += clip_aln->rightClipSize;

	// cout << "clip_aln->querylen=" << clip_aln->querylen << endl;
	//clip_aln->alnsize = clip_aln->querylen - clip_aln->leftClipSize - clip_aln->rightClipSize;

	// orientation
	if(bam_is_rev(b)) clip_aln->aln_orient = ALN_MINUS_ORIENT;
	else clip_aln->aln_orient = ALN_PLUS_ORIENT;

	// query positions
	if(clip_aln->aln_orient==ALN_PLUS_ORIENT){ // plus orient
		clip_aln->startQueryPos = clip_aln->leftClipSize + 1;
		clip_aln->endQueryPos = clip_aln->querylen - clip_aln->rightClipSize;
	}else{ // minus orient
		clip_aln->startQueryPos = clip_aln->querylen - clip_aln->leftClipSize;
		clip_aln->endQueryPos = clip_aln->rightClipSize + 1;
	}

	clip_aln->ref_dist = clip_aln->endRefPos - clip_aln->startRefPos + 1;
	if(clip_aln->aln_orient==ALN_PLUS_ORIENT) clip_aln->query_dist = clip_aln->endQueryPos - clip_aln->startQueryPos + 1;
	else clip_aln->query_dist = clip_aln->startQueryPos - clip_aln->endQueryPos + 1;

	clip_aln->query_checked_flag = false;
	clip_aln->left_clip_checked_flag = false;
	clip_aln->right_clip_checked_flag = false;
	clip_aln->SA_tag_flag = false;

	return clip_aln;
}

// align data sampling
void clipAlnDataLoader::samplingAlnData(vector<bam1_t*> &alnDataVector, double mean_read_len, double max_ultra_high_cov){
	double compensation_coefficient, local_cov_original;

	compensation_coefficient = computeCompensationCoefficient(startRefPos, endRefPos, mean_read_len);
	local_cov_original = computeLocalCov(alnDataVector, mean_read_len, compensation_coefficient);
	if(local_cov_original>max_ultra_high_cov)  // sampling
		samplingAlnDataOp(alnDataVector, mean_read_len, max_ultra_high_cov);
}

// align data sampling operation
double clipAlnDataLoader::samplingAlnDataOp(vector<bam1_t*> &alnDataVector, double mean_read_len, double expect_cov_val){
	double expected_total_bases, total_bases, sampled_cov;
	size_t index, max_reads_num, reg_size, reads_count;
	bam1_t *b;
	int8_t *selected_flag_array;
	int64_t i;

	sampled_cov = 0;
	selected_flag_array = (int8_t*) calloc(alnDataVector.size(), sizeof(int8_t));
	if(selected_flag_array==NULL){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot allocate memory, error!" << endl;
		exit(1);
	}

	reg_size = endRefPos - startRefPos + 1 + mean_read_len;
	expected_total_bases = reg_size * expect_cov_val;
	max_reads_num = alnDataVector.size();

	// make sure each down-sampling is equivalent
	//pthread_mutex_lock(&mutex_down_sample);
	srand(1);
	reads_count = 0;
	total_bases = 0;
	while(total_bases<=expected_total_bases and reads_count<=max_reads_num){
		index = rand() % max_reads_num;
		if(selected_flag_array[index]==0){
			selected_flag_array[index] = 1;

			b = alnDataVector.at(index);
			total_bases += b->core.l_qseq;
			reads_count ++;
		}
	}
	//pthread_mutex_unlock(&mutex_down_sample);
	sampled_cov = total_bases / reg_size;

	// remove unselected items
	for(i=alnDataVector.size()-1; i>=0; i--){
		if(selected_flag_array[i]==0){
			b = alnDataVector.at(i);
			bam_destroy1(b);
			alnDataVector.erase(alnDataVector.begin()+i);
		}
	}

	free(selected_flag_array);

	//cout << "After sampling, reads count: " << reads_count << ", total bases: " << (int64_t)total_bases << ", alnDataVector.size: " << alnDataVector.size() << endl;

	return sampled_cov;
}

double clipAlnDataLoader::computeLocalCov(vector<bam1_t*> &alnDataVector, double mean_read_len, double compensation_coefficient){
	double cov = 0, total;
	bam1_t *b;
	size_t refined_reg_size;

	refined_reg_size = endRefPos - startRefPos + 1 + 2 * mean_read_len;
	if(refined_reg_size>0){
		total = 0;
		for(size_t i=0; i<alnDataVector.size(); i++){
			b = alnDataVector.at(i);
			total += b->core.l_qseq;
		}
		cov = total / refined_reg_size * compensation_coefficient;
		//cout << "total bases: " << total << " bp, local coverage: " << cov << endl;
	}else{
		cov = -1;
		//cerr << "ERR: ref_seq_size=" << ref_seq_size << endl;
	}
	return cov;
}

double clipAlnDataLoader::computeCompensationCoefficient(size_t startRefPos, size_t endRefPos, double mean_read_len){
	size_t total_reg_size, refined_reg_size, reg_size_cns;
	double comp_coefficient;

	reg_size_cns = endRefPos - startRefPos + 1;
	total_reg_size = reg_size_cns + 2 * mean_read_len;
	refined_reg_size = reg_size_cns + mean_read_len;  // flanking_area = (2 * mean_read_len) / 2
	comp_coefficient = (double)total_reg_size / refined_reg_size;

	return comp_coefficient;
}

// fill data according to 'SA' tag
void clipAlnDataLoader::fillClipAlnDataBySATag(vector<clipAlnData_t*> &clipAlnDataVector){
	size_t i, j;
	clipAlnData_t *clip_aln;
	uint8_t *cigar_int;
	string cigar_str, aln_seg_info_str;
	vector<string> aln_seg_vec;

	for(i=0; i<clipAlnDataVector.size(); i++){
		clip_aln = clipAlnDataVector.at(i);
		if(clip_aln->SA_tag_flag==false){
			cigar_int = bam_aux_get(clip_aln->bam, "SA"); // SA
			if(cigar_int) {
				cigar_str = bam_aux2Z(cigar_int);
				aln_seg_vec = split(cigar_str, ";");
				for(j=0; j<aln_seg_vec.size(); j++){
					aln_seg_info_str = aln_seg_vec.at(j);
					//cout << aln_seg_vec.at(j) << endl;
					addNewSAItemToClipAlnDataVec(clip_aln->queryname,  aln_seg_info_str, clipAlnDataVector);
				}

				// cout << clip_aln->queryname << ":" << clip_aln->startRefPos << "-" << clip_aln->endRefPos << endl;
				// cout << cigar_str << endl;
			}
		}
	}
}

// add new SA item to clipAlnDataVector
clipAlnData_t* clipAlnDataLoader::addNewSAItemToClipAlnDataVec(string &queryname, string &aln_seg_info_str, vector<clipAlnData_t*> &clipAlnDataVector){
	clipAlnData_t clip_aln_tmp, *clip_aln = NULL, *clip_aln_new = NULL;
	size_t i;
	bool new_flag;

	// parse alignments in the 'SA' tag
	clip_aln_tmp.queryname = queryname;
	parseSingleAlnStrSA(clip_aln_tmp, aln_seg_info_str);

	new_flag = true;
	for(i=0; i<clipAlnDataVector.size(); i++){
		clip_aln = clipAlnDataVector.at(i);
		if(isSameClipAlnSeg(clip_aln, &clip_aln_tmp)) { new_flag = false; break;}
	}

	// add new item
	if(new_flag){
		clip_aln_new = new clipAlnData_t();
		clip_aln_new->bam = NULL;
		clip_aln_new->queryname = queryname;
		clip_aln_new->chrname = clip_aln_tmp.chrname;
		clip_aln_new->querylen = clip_aln_tmp.querylen;
		clip_aln_new->aln_orient = clip_aln_tmp.aln_orient;
		clip_aln_new->startRefPos = clip_aln_tmp.startRefPos;
		clip_aln_new->endRefPos = clip_aln_tmp.endRefPos;
		clip_aln_new->startQueryPos = clip_aln_tmp.startQueryPos;
		clip_aln_new->endQueryPos = clip_aln_tmp.endQueryPos;
		clip_aln_new->leftClipSize = clip_aln_tmp.leftClipSize;
		clip_aln_new->rightClipSize = clip_aln_tmp.rightClipSize;
		clip_aln_new->ref_dist = clip_aln_tmp.endRefPos - clip_aln_tmp.startRefPos + 1;
		if(clip_aln_new->aln_orient==ALN_PLUS_ORIENT) clip_aln_new->query_dist = clip_aln_tmp.endQueryPos - clip_aln_tmp.startQueryPos + 1;
		else clip_aln_new->query_dist = clip_aln_tmp.startQueryPos - clip_aln_tmp.endQueryPos + 1;
		clip_aln_new->left_clip_checked_flag = false;
		clip_aln_new->right_clip_checked_flag = false;
		clip_aln_new->query_checked_flag = false;
		clip_aln_new->SA_tag_flag = true;
		clip_aln_new->left_aln = clip_aln_new->right_aln = NULL;
		clipAlnDataVector.push_back(clip_aln_new);
	}

	return clip_aln_new;
}

// parse cigar string
void clipAlnDataLoader::parseSingleAlnStrSA(clipAlnData_t &clip_aln_ret, string &aln_seg_info_str){
	size_t i, ref_aln_size, query_aln_size, op_len;
	char ch, op_ch;
	string cigar_str, str_tmp;
	vector<string> aln_info_vec;
	vector<size_t> op_len_vec;
	vector<char> op_vec;

	aln_info_vec = split(aln_seg_info_str, ",");
	clip_aln_ret.bam = NULL;
	clip_aln_ret.chrname = aln_info_vec.at(0);   // chrname
	clip_aln_ret.startRefPos = stoi(aln_info_vec.at(1));  // startRefPos
	if(aln_info_vec.at(2).compare("+")==0) clip_aln_ret.aln_orient = ALN_PLUS_ORIENT;  // aln_orient
	else clip_aln_ret.aln_orient = ALN_MINUS_ORIENT;

	// cigar
	cigar_str = aln_info_vec.at(3);
	str_tmp = "";
	for(i=0; i<cigar_str.size(); i++){
		ch = cigar_str.at(i);
		str_tmp += ch;
		if(ch>='A' and ch<='Z') {
			op_len = stoi(str_tmp.substr(0, str_tmp.size()-1));
			op_ch = str_tmp.at(str_tmp.size()-1);

			op_len_vec.push_back(op_len);
			op_vec.push_back(op_ch);

			str_tmp = "";
		}
	}

	// NM
	// clip_aln_ret.NM = stoi(aln_info_vec.at(aln_info_vec.size()-1));

	// leftClipSize and leftClipSize
	if(op_vec.at(0)=='S' or op_vec.at(0)=='H') clip_aln_ret.leftClipSize = op_len_vec.at(0);
	else clip_aln_ret.leftClipSize = 0;
	if(op_vec.at(op_vec.size()-1)=='S' or op_vec.at(op_vec.size()-1)=='H') clip_aln_ret.rightClipSize = op_len_vec.at(op_vec.size()-1);
	else clip_aln_ret.rightClipSize = 0;

	// compute the align positions
	ref_aln_size = query_aln_size = 0;
	for(i=0; i<op_vec.size(); i++){
		op_len = op_len_vec.at(i);
		switch(op_vec.at(i)){
			case 'S':
			case 'H':
				break;
			case '=':
			case 'X':
			case 'M':
				ref_aln_size += op_len;
				query_aln_size += op_len;
				break;
			case 'I':
				query_aln_size += op_len;
				break;
			case 'D':
				ref_aln_size += op_len;
				break;
			default: cerr << __func__ << ": invalid op flag: " << op_vec.at(i) << endl; exit(1);
		}
	}

	clip_aln_ret.querylen = clip_aln_ret.leftClipSize + query_aln_size + clip_aln_ret.rightClipSize;
	clip_aln_ret.endRefPos = clip_aln_ret.startRefPos + ref_aln_size - 1;
	if(clip_aln_ret.aln_orient==ALN_PLUS_ORIENT){
		clip_aln_ret.startQueryPos = clip_aln_ret.leftClipSize + 1;
		clip_aln_ret.endQueryPos = clip_aln_ret.querylen - clip_aln_ret.rightClipSize;
	}else{
		clip_aln_ret.startQueryPos = clip_aln_ret.querylen - clip_aln_ret.leftClipSize;
		clip_aln_ret.endQueryPos = clip_aln_ret.rightClipSize + 1;
	}

	clip_aln_ret.left_clip_checked_flag = clip_aln_ret.right_clip_checked_flag = clip_aln_ret.query_checked_flag = false;
}

// determine whether the two clip align segments are the same
bool clipAlnDataLoader::isSameClipAlnSeg(clipAlnData_t *clip_aln1, clipAlnData_t *clip_aln2){
	bool flag = false;
	if(clip_aln1->queryname.compare(clip_aln2->queryname)==0 and clip_aln1->chrname.compare(clip_aln2->chrname)==0
		and clip_aln1->startRefPos==clip_aln2->startRefPos and clip_aln1->endRefPos==clip_aln2->endRefPos
		and clip_aln1->startQueryPos==clip_aln2->startQueryPos and clip_aln1->endQueryPos==clip_aln2->endQueryPos
		and clip_aln1->aln_orient==clip_aln2->aln_orient){
		flag = true;
	}
	return flag;
}

// release the memory
void clipAlnDataLoader::freeClipAlnData(vector<clipAlnData_t*> &clipAlnDataVector){
	for(size_t i=0; i<clipAlnDataVector.size(); i++){
		bam_destroy1(clipAlnDataVector.at(i)->bam);
		delete clipAlnDataVector.at(i);
	}
	vector<clipAlnData_t*>().swap(clipAlnDataVector);
}

// remove single item in clipAlnDataVector
void clipAlnDataLoader::removeSingleItemClipAlnData(vector<clipAlnData_t*> &clipAlnDataVector, int32_t idx){
	bam_destroy1(clipAlnDataVector.at(idx)->bam);
	delete clipAlnDataVector.at(idx);
	clipAlnDataVector.erase(clipAlnDataVector.begin() + idx);
}

// add adjacent info of align segments
void clipAlnDataLoader::addAdjacentInfo(vector<clipAlnData_t*> &clipAlnDataVector){
	size_t i;
	string queryname;
	vector<clipAlnData_t*> query_aln_segs;
	clipAlnData_t *clip_aln_seg;

	for(i=0; i<clipAlnDataVector.size(); i++){
		queryname = clipAlnDataVector.at(i)->queryname;

		//query_aln_segs = getQueryClipAlnSegs(queryname, clipAlnDataVector);  // get query clip align segments
		query_aln_segs = getQueryClipAlnSegsAll(queryname, clipAlnDataVector);  // get query clip align segments

		// order clipping segments
		orderClipAlnSegsSingleQuery(query_aln_segs);

		// assign side most flag
		assignSideMostFlag(query_aln_segs);
	}

	// reset flags
	for(i=0; i<clipAlnDataVector.size(); i++){
		clip_aln_seg = clipAlnDataVector.at(i);
		clip_aln_seg->left_clip_checked_flag = clip_aln_seg->right_clip_checked_flag = false;
	}
}

// order clipping segments
void clipAlnDataLoader::orderClipAlnSegsSingleQuery(vector<clipAlnData_t*> &query_aln_vec){
	size_t i;
	int32_t clip_end_flag, mate_arr_idx, mate_clip_end_flag;
	bool valid_query_end_flag;
	clipAlnData_t *clip_aln_seg, *mate_clip_aln_seg;
	vector<int32_t> adjClipAlnSegInfo;

	for(i=0; i<query_aln_vec.size(); i++){
		clip_aln_seg = query_aln_vec.at(i);

		// clipping at right end
		if(clip_aln_seg->right_clip_checked_flag==false){
			valid_query_end_flag = false;
			clip_end_flag = -1;
			if(clip_aln_seg->rightClipSize>=minClipEndSize){
				valid_query_end_flag = true;
				clip_end_flag = RIGHT_END;
			}

			if(valid_query_end_flag){
				// deal with the mate clip end
				adjClipAlnSegInfo = getAdjacentClipAlnSeg(i, clip_end_flag, query_aln_vec, minClipEndSize, MAX_VAR_REG_SIZE);
				mate_arr_idx = adjClipAlnSegInfo.at(0);
				mate_clip_end_flag = adjClipAlnSegInfo.at(1);
				if(mate_arr_idx!=-1){ // mated
					mate_clip_aln_seg = query_aln_vec.at(mate_arr_idx);
					clip_aln_seg->right_aln = mate_clip_aln_seg;
					if(mate_clip_end_flag==LEFT_END and mate_clip_aln_seg->left_clip_checked_flag==false){ // mate at left end
						mate_clip_aln_seg->left_aln = clip_aln_seg;
						mate_clip_aln_seg->left_clip_checked_flag = true;
					}else if(mate_clip_end_flag==RIGHT_END and mate_clip_aln_seg->right_clip_checked_flag==false){  // mate at right end
						mate_clip_aln_seg->right_aln = clip_aln_seg;
						mate_clip_aln_seg->right_clip_checked_flag = true;
					}else {
						cerr << __func__ << ": invalid mate_clip_end_flag=" << mate_clip_end_flag << endl;
						exit(1);
					}
				}
			}
			clip_aln_seg->right_clip_checked_flag = true;
		}

		// clipping at left end
		if(clip_aln_seg->left_clip_checked_flag==false){
			valid_query_end_flag = false;
			clip_end_flag = -1;
			if(clip_aln_seg->leftClipSize>=minClipEndSize){
				valid_query_end_flag = true;
				clip_end_flag = LEFT_END;
			}

			if(valid_query_end_flag){
				// deal with the mate clip end
				adjClipAlnSegInfo = getAdjacentClipAlnSeg(i, clip_end_flag, query_aln_vec, minClipEndSize, MAX_VAR_REG_SIZE);
				mate_arr_idx = adjClipAlnSegInfo.at(0);
				mate_clip_end_flag = adjClipAlnSegInfo.at(1);
				if(mate_arr_idx!=-1){ // mated
					mate_clip_aln_seg = query_aln_vec.at(mate_arr_idx);
					clip_aln_seg->left_aln = mate_clip_aln_seg;
					if(mate_clip_end_flag==LEFT_END and mate_clip_aln_seg->left_clip_checked_flag==false){ // mate at left end
						mate_clip_aln_seg->left_aln = clip_aln_seg;
						mate_clip_aln_seg->left_clip_checked_flag = true;
					}else if(mate_clip_end_flag==RIGHT_END and mate_clip_aln_seg->right_clip_checked_flag==false){  // mate at right end
						mate_clip_aln_seg->right_aln = clip_aln_seg;
						mate_clip_aln_seg->right_clip_checked_flag = true;
					}else {
						cerr << __func__ << ": invalid mate_clip_end_flag=" << mate_clip_end_flag << endl;
						exit(1);
					}
				}
			}
			clip_aln_seg->left_clip_checked_flag = true;
		}
	}
}

// assign the side most flag
void clipAlnDataLoader::assignSideMostFlag(vector<clipAlnData_t*> &query_aln_vec){
	size_t i;
	clipAlnData_t *clip_aln_seg;
	bool flag;
	int64_t min_id, max_id, minValue, maxValue, tmp_min, tmp_max;

	flag = false;
	for(i=0; i<query_aln_vec.size(); i++){
		clip_aln_seg = query_aln_vec.at(i);
		if(clip_aln_seg->leftmost_flag==true or clip_aln_seg->rightmost_flag==true){
			flag = true;
			break;
		}
	}

	if(flag==false){
		min_id = max_id = -1;
		minValue = INT_MAX, maxValue = INT_MIN;
		for(i=0; i<query_aln_vec.size(); i++){
			clip_aln_seg = query_aln_vec.at(i);

			if(clip_aln_seg->startQueryPos<clip_aln_seg->endQueryPos){
				tmp_min = clip_aln_seg->startQueryPos;
				tmp_max = clip_aln_seg->endQueryPos;
			}else{
				tmp_min = clip_aln_seg->endQueryPos;
				tmp_max = clip_aln_seg->startQueryPos;
			}

			if(tmp_min<minValue){
				min_id = i;
				minValue = tmp_min;
			}
			if(tmp_max>maxValue){
				max_id = i;
				maxValue = tmp_max;
			}
		}
		if(min_id!=-1 and max_id!=-1){
			if(query_aln_vec.at(min_id)->startRefPos>query_aln_vec.at(max_id)->startRefPos){ // swap
				query_aln_vec.at(max_id)->leftmost_flag = true;
				query_aln_vec.at(min_id)->rightmost_flag = true;
			}else{
				query_aln_vec.at(min_id)->leftmost_flag = true;
				query_aln_vec.at(max_id)->rightmost_flag = true;
			}
		}else{
			cout << "line=" << __LINE__ << ", can not find the left and right most align segments, error!" << endl;
			exit(1);
		}
	}
}

// remove clip alignment data with low primary segment size ratio
void clipAlnDataLoader::removeClipAlnDataWithLowPrimarySegSizeRatio(vector<clipAlnData_t*> &clipAlnDataVector, double primary_seg_size_ratio, double primary_seg_nm_ratio){
	size_t i, j, n;
	string queryname;
	vector<clipAlnData_t*> query_aln_segs;//query_aln_segs_no_SA;
	int64_t max_id, max_len;
	vector<string> remove_qname_vec;
	double max_primary_size_ratio, max_primary_seg_nm_ratio;
	vector<string>::iterator q;

	if(primary_seg_size_ratio>0){

		for(i=0; i<clipAlnDataVector.size(); i++){
			queryname = clipAlnDataVector.at(i)->queryname;
			
			q = find(remove_qname_vec.begin(), remove_qname_vec.end(), queryname);
			if(q == remove_qname_vec.end()){
				// get query_aln_segs and query_aln_segs_no_SA
				query_aln_segs = getQueryClipAlnSegsAll(queryname, clipAlnDataVector);  // get all query clip align segments
				//query_aln_segs_no_SA = getQueryClipAlnSegs(queryname, clipAlnDataVector);  // get all query clip align segments without SA tag

				if(query_aln_segs.size()>1){
					// if the max_primary_size_ratio < primary_seg_size_ratio, then remove these segs from clipAlnDataVector
					max_len = INT_MIN;
					max_id = -1;
					for(j=0; j<query_aln_segs.size(); j++){
						if(query_aln_segs.at(j)->SA_tag_flag==false and query_aln_segs.at(j)->query_dist>max_len){
							max_len = query_aln_segs.at(j)->query_dist;
							max_id = j;
						}
					}
					max_primary_size_ratio = (double)max_len/clipAlnDataVector.at(i)->querylen;

					if(query_aln_segs.at(max_id)->SA_tag_flag==false){
						max_primary_seg_nm_ratio = (double)query_aln_segs.at(max_id)->NM/query_aln_segs.at(max_id)->query_dist;
						// cout << clipAlnDataVector.at(i)->chrname << ":" << queryname << " query_aln_segs.at(max_id)->NM=" << query_aln_segs.at(max_id)->NM << ",query_aln_segs.at(max_id)->query_dist=" << query_aln_segs.at(max_id)->query_dist <<",primary_seg_nm_ratio=" << max_primary_seg_nm_ratio << endl;					
					}

					if(max_primary_size_ratio<primary_seg_size_ratio and max_primary_seg_nm_ratio>primary_seg_nm_ratio){ // primary_seg_nm_ratio ccs 0.1, clr/rs/ont 0.2
						remove_qname_vec.push_back(query_aln_segs.at(max_id)->queryname);
					}
				}
			}
		}

		//remove unreliable alignment segments
		for(n=0; n<clipAlnDataVector.size();){
			q = find(remove_qname_vec.begin(), remove_qname_vec.end(), clipAlnDataVector.at(n)->queryname);
			if(q != remove_qname_vec.end()) removeSingleItemClipAlnData(clipAlnDataVector, n);
			else n++;
		}
	}
}
