#include "alnDataLoader.h"
#include "clipAlnDataLoader.h"
#include "util.h"

clipAlnDataLoader::clipAlnDataLoader(string &chrname, int32_t startRefPos, int32_t endRefPos, string &inBamFile) {
	this->chrname = chrname;
	this->startRefPos = startRefPos;
	this->endRefPos = endRefPos;
	this->inBamFile = inBamFile;
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
	alnDataLoader data_loader(chrname, startRefPos, endRefPos, inBamFile);
	data_loader.loadAlnData(alnDataVector);

	if(max_ultra_high_cov>0){
		samplingAlnData(alnDataVector, data_loader.mean_read_len, max_ultra_high_cov);
	}

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
}

void clipAlnDataLoader::loadClipAlnDataWithSATag(vector<clipAlnData_t*> &clipAlnDataVector, double max_ultra_high_cov){
	loadClipAlnData(clipAlnDataVector, max_ultra_high_cov);
	fillClipAlnDataBySATag(clipAlnDataVector);
}

clipAlnData_t* clipAlnDataLoader::generateClipAlnData(bam1_t* b, bam_hdr_t *header){
	clipAlnData_t *clip_aln = NULL;
	uint32_t *c, op1, op2;

	clip_aln = new clipAlnData_t();
	clip_aln->bam = b;
	clip_aln->queryname = bam_get_qname(b);
	clip_aln->chrname = header->target_name[b->core.tid];
	clip_aln->startRefPos = b->core.pos + 1;
	clip_aln->endRefPos = bam_endpos(b);

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

	// querylen
	clip_aln->querylen = b->core.l_qseq;
	if(op1==BAM_CHARD_CLIP)
		clip_aln->querylen += clip_aln->leftClipSize;
	if(op2==BAM_CHARD_CLIP)
		clip_aln->querylen += clip_aln->rightClipSize;

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
	size_t total_reg_size, refined_reg_size, reg_size_assemble;
	double comp_coefficient;

	reg_size_assemble = endRefPos - startRefPos + 1;
	total_reg_size = reg_size_assemble + 2 * mean_read_len;
	refined_reg_size = reg_size_assemble + mean_read_len;  // flanking_area = (2 * mean_read_len) / 2
	comp_coefficient = (double)total_reg_size / refined_reg_size;

	return comp_coefficient;
}

// fill data according to 'SA' tag
void clipAlnDataLoader::fillClipAlnDataBySATag(vector<clipAlnData_t*> &clipAlnDataVector){
	size_t i, j;
	clipAlnData_t *clip_aln;
	uint8_t *ciger_int;
	string cigar_str, aln_seg_info_str;
	vector<string> aln_seg_vec;

	for(i=0; i<clipAlnDataVector.size(); i++){
		clip_aln = clipAlnDataVector.at(i);
		if(clip_aln->SA_tag_flag==false){
			ciger_int = bam_aux_get(clip_aln->bam, "SA"); // SA
			if(ciger_int) {
				cigar_str = bam_aux2Z(ciger_int);

				aln_seg_vec = split(cigar_str, ";");
				for(j=0; j<aln_seg_vec.size(); j++){
					aln_seg_info_str = aln_seg_vec.at(j);
					//cout << aln_seg_vec.at(j) << endl;
					addNewSAItemToClipAlnDataVec(clip_aln->queryname,  aln_seg_info_str, clipAlnDataVector);
				}

				//cout << clip_aln->queryname << ":" << clip_aln->startRefPos << "-" << clip_aln->endRefPos << endl;
				//cout << cigar_str << endl;
			}
		}
	}
}

// add new SA item to clipAlnDataVector
clipAlnData_t* clipAlnDataLoader::addNewSAItemToClipAlnDataVec(string &queryname, string &aln_seg_info_str, vector<clipAlnData_t*> &clipAlnDataVector){
	clipAlnData_t clip_aln_tmp, *clip_aln = NULL, *clip_aln_new;
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
		clip_aln_new->left_clip_checked_flag = false;
		clip_aln_new->right_clip_checked_flag = false;
		clip_aln_new->query_checked_flag = false;
		clip_aln_new->SA_tag_flag = true;
		clipAlnDataVector.push_back(clip_aln_new);
	}

	return clip_aln;
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
