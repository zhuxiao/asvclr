#include <math.h>
#include "varCand.h"
#include "covLoader.h"
#include "util.h"

//pthread_mutex_t mutex_print_var_cand = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mutex_write_blat_aln = PTHREAD_MUTEX_INITIALIZER;

varCand::varCand(){
	this->chrname = "";
	this->var_cand_filename = "";
	this->out_dir_call = "";
	inBamFile = "";
	fai = NULL;

	init();
}

varCand::~varCand() {
	if(clip_reg) delete clip_reg;
	if(limit_reg_delete_flag){
		for(size_t i=0; i<sub_limit_reg_vec.size(); i++)
			delete sub_limit_reg_vec.at(i);
		vector<simpleReg_t*>().swap(sub_limit_reg_vec);
		limit_reg_delete_flag = false;
	}
}

void varCand::init(){
	refseqfilename = "";
	ctgfilename = "";
	readsfilename = "";
	alnfilename = "";
	ref_left_shift_size = ref_right_shift_size = ctg_num = 0;
	assem_success = align_success = call_success = clip_reg_flag = false;
	clip_reg = NULL;
	leftClipRefPos = rightClipRefPos = 0;
	sv_type = VAR_UNC;
	dup_num = 0;
	margin_adjusted_flag = false;
	blat_var_cand_file = NULL;
	limit_reg_delete_flag = false;
}

void varCand::destroyVarCand(){
	size_t i, j, k;
	blat_aln_t *item_blat;

	for(i=0; i<varVec.size(); i++)  // free varVec
		delete varVec.at(i);
	for(i=0; i<newVarVec.size(); i++)  // free newVarVec
		delete newVarVec.at(i);

	for(j=0; j<blat_aln_vec.size(); j++){ // free blat_aln_vec
		item_blat = blat_aln_vec.at(j);
		for(k=0; k<item_blat->aln_segs.size(); k++) // free aln_segs
			delete item_blat->aln_segs.at(k);
		delete item_blat;
	}
}

void varCand::setBlatVarcandFile(ofstream *blat_var_cand_file, vector<varCand*> *blat_aligned_info_vec){
	this->blat_var_cand_file = blat_var_cand_file;
	this->blat_aligned_info_vec = blat_aligned_info_vec;
}

void varCand::resetBlatVarcandFile(){
	blat_var_cand_file = NULL;
	blat_aligned_info_vec = NULL;
}

// align contig to refseq
void varCand::alnCtg2Refseq(){
	//cout << "BLAT align: " << alnfilename << endl;

	bool blat_aln_done_flag = getBlatAlnDoneFlag();

	//if(assem_success and !isFileExist(alnfilename)){
	if(assem_success and blat_aln_done_flag==false){
		if(!isFileExist(alnfilename) or !isBlatAlnResultMatch(ctgfilename, alnfilename) or !isBlatAlnCompleted(alnfilename)) // file not exist, or query names not match, or blat align uncompleted
			blatAln(alnfilename, ctgfilename, refseqfilename); // BLAT alignment

		// record blat aligned information
		recordBlatAlnInfo();
	}
}

// record blat aligned information
void varCand::recordBlatAlnInfo(){
	string line, aln_flag_str, limit_reg_str;
	bool aln_flag;

	if(blat_var_cand_file){
		aln_flag = isFileExist(alnfilename);
		if(aln_flag) aln_flag_str = ALIGN_SUCCESS;
		else aln_flag_str = ALIGN_FAILURE;

		// limit regions
		if(limit_reg_process_flag) limit_reg_str = getLimitRegStr(sub_limit_reg_vec);
		else limit_reg_str = LIMIT_REG_ALL_STR;

		line = alnfilename + "\t" + ctgfilename + "\t" + refseqfilename + "\t" + aln_flag_str + "\t" + limit_reg_str + "\t" + DONE_STR;

		pthread_mutex_lock(&mutex_write_blat_aln);
		*blat_var_cand_file << line << endl;
		pthread_mutex_unlock(&mutex_write_blat_aln);
	}
}

// get blat align done flag
bool varCand::getBlatAlnDoneFlag(){
	bool blat_aln_done_flag;
	varCand *blat_aligned_info_node;

	blat_aln_done_flag = false;
	if(blat_aligned_info_vec){
		for(size_t i=0; i<blat_aligned_info_vec->size(); i++){
			blat_aligned_info_node = blat_aligned_info_vec->at(i);
			if(blat_aligned_info_node->alnfilename.compare(alnfilename)==0){
				blat_aln_done_flag = true;
				break;
			}
		}
	}

	return blat_aln_done_flag;
}


// call variants according to BLAT alignment
void varCand::callVariants(){
	if(sv_type!=VAR_TRA){
//		pthread_mutex_lock(&mutex_print_var_cand);
//		cout << "process: " << alnfilename << endl;
//		pthread_mutex_unlock(&mutex_print_var_cand);

		loadBlatAlnData(); // load align data

		if(clip_reg_flag==false)	// call indel variants
			callIndelVariants();
		else // call clipping variants
			callClipRegVariants();

		// call variants at align segment end
		callVariantsAlnSegEnd();

		// destroy local alignment items
		destroyLocalAlnVec(local_aln_vec);
	}
}

// load align data for TRA
void varCand::loadBlatAlnData(){
	// BLAT align
#if BLAT_ALN
	alnCtg2Refseq();
#endif
	assignBlatAlnStatus();
	if(align_success){
		// parse alignment
		blatParse();

		// filter the blat align data
		blatFilter();
	}else{
//		pthread_mutex_lock(&mutex_print_var_cand);
//		cout << "Align Failure: " << ctgfilename << endl;
//		pthread_mutex_unlock(&mutex_print_var_cand);
	}
}

// call indel variants
void varCand::callIndelVariants(){
	if(align_success){
		// call indel type
		determineIndelType();

		// call short variants whose alignment cannot be split at the variant location
		callShortVariants();

		// distinguish short duplications and inversions from indels
		distinguishShortDupInvFromIndels();

		// print SV
		//printSV();
	}
}

// parse and call clipReg variants
void varCand::callClipRegVariants(){
	if(align_success){
		// call variants according to variant types
		determineClipRegVarType();

		// print SV
		//printSV();
	}
}

// update align status
void varCand::assignBlatAlnStatus(){
	align_success = isFileExist(alnfilename);
}

// parse BLAT alignments
void varCand::blatParse(){
	string line, q_reg_str, s_reg_str, ident_percent_str, orient_str, query_name;
	vector<string> line_vec, line_vec1, len_vec, str_vec, query_name_vec;
	ifstream infile;
	int32_t ref_start_all, query_loc, blat_aln_id;
	blat_aln_t *blat_aln_item;
	aln_seg_t *aln_seg;

	if(blat_aln_vec.size()>0) return;

	// load query names
	FastaSeqLoader fa_loader(ctgfilename);
	query_name_vec = fa_loader.getFastaSeqNames();

	infile.open(alnfilename);
	if(!infile.is_open()){
		cerr << __func__ << "(), line=" << __LINE__ << ": cannot open file:" << alnfilename << endl;
		exit(1);
	}

	// reference start position
	ref_start_all = varVec[0]->startRefPos - ref_left_shift_size;

	// fill the data
	blat_aln_id = 0;
	while(getline(infile, line)){
		if(line.size()){
			if(line.substr(0, 3).compare("seq")==0){  // query and subject titles
				line_vec = split(line, "=");
				line_vec1 = split(line_vec[1].substr(1), ",");

				len_vec = split(line_vec1[line_vec1.size()-1].substr(1), " ");
				if(line_vec[0][3]=='1') {
					// get the query name
					query_loc = getQueryNameLoc(line_vec1.at(0), query_name_vec);
					if(query_loc==-1){
						//query_name_vec.push_back(line_vec1[0]);
						//query_loc = query_name_vec.size() - 1;
						cerr << __func__ << "(), line=" << __LINE__ << ": cannot get queryloc by query_name: " << line_vec1.at(0) << endl;
						exit(1);
					}

					blat_aln_item = new blat_aln_t;
					blat_aln_item->blat_aln_id = blat_aln_id++;
					blat_aln_item->query_id = query_loc;
					blat_aln_item->query_len = stoi(len_vec[0]); // query length
					blat_aln_item->head_hanging = blat_aln_item->tail_hanging = false;
					blat_aln_item->best_aln = false;
					blat_aln_item->valid_aln = false;
					blat_aln_item->aln_orient = ALN_PLUS_ORIENT;
					blat_aln_vec.push_back(blat_aln_item);
				}
				else blat_aln_item->subject_len = stoi(len_vec[0]); // subject length
			}else if(line[0]=='('){ // complement
				if(line.substr(1, 10).compare("complement")==0)
					blat_aln_item->aln_orient = ALN_MINUS_ORIENT;
			}else{    // segments
				line_vec = split(line, " ");

				if(line_vec.size()==4){ // tolerate killed BLAT results
					q_reg_str = line_vec[0];   // query
					s_reg_str = line_vec[1];   // subject
					ident_percent_str = line_vec[2];  // identity percent
					orient_str = line_vec[3];  // orientation

					aln_seg = new aln_seg_t;   // allocate memory

					str_vec = split(q_reg_str, "-");
					aln_seg->query_start = stoi(str_vec[0]);  // query start
					aln_seg->query_end = stoi(str_vec[1]);  // query end

					str_vec = split(s_reg_str, "-");
					aln_seg->subject_start = stoi(str_vec[0].substr(1));  // subject start
					aln_seg->subject_end = stoi(str_vec[1].substr(0, str_vec[1].size()-1));  // subject end

					try{
						aln_seg->ident_percent = stof(ident_percent_str) / 100.0;  // identity percent
					}catch(const std::invalid_argument& ia){
						cerr << "Invalid argument: " << ia.what() << endl;
					}

					if(orient_str.compare("->")==0 or orient_str.at(0)=='-') aln_seg->aln_orient = ALN_PLUS_ORIENT;  // orientation
					else aln_seg->aln_orient = ALN_MINUS_ORIENT;

					// compute the reference position
					aln_seg->ref_start = ref_start_all + aln_seg->subject_start - 1;
					aln_seg->ref_end = ref_start_all + aln_seg->subject_end - 1;

					blat_aln_item->aln_segs.push_back(aln_seg);
				}
			}
		}
	}
	infile.close();

	// assign best_aln flag
	assignBestAlnFlag();

	// assign align orientation flag
	assignAlnOrientFlag();

	// resize the aln_segs
	for(size_t i=0; i<blat_aln_vec.size(); i++) {
		blat_aln_item = blat_aln_vec.at(i);
		blat_aln_item->aln_segs.shrink_to_fit();
	}
}

// assign best_aln flag
void varCand::assignBestAlnFlag(){
	int32_t best_aln_id, query_id, query_num, sum, max_sum, max_sum_id, secMax_sum, secMax_sum_id, valid_seg_num, max_valid_seg_num, secMax_valid_seg_num, tmp_len;
	vector<aln_seg_t*> aln_segs;
	blat_aln_t *blat_aln;
	double ratio, aver_len_max, aver_len_secMax;

	if(blat_aln_vec.size()>0){
		blat_aln = blat_aln_vec.at(blat_aln_vec.size()-1);
		query_num = blat_aln->query_id + 1;
		for(query_id=0; query_id<query_num; query_id++){
			max_sum = 0; max_sum_id = 0;
			secMax_sum = 0; secMax_sum_id = 0;
			max_valid_seg_num = 0; secMax_valid_seg_num = 0;
			for(size_t i=0; i<blat_aln_vec.size(); i++){
				blat_aln = blat_aln_vec.at(i);
				if(blat_aln->query_id==query_id){
					aln_segs = blat_aln->aln_segs;
					sum = 0;
					valid_seg_num = 0;
					for(size_t j=0; j<aln_segs.size(); j++){
						tmp_len = aln_segs.at(j)->subject_end - aln_segs.at(j)->subject_start + 1;
						if(tmp_len>=MIN_VALID_BLAT_SEG_SIZE){
							valid_seg_num ++;
							sum += tmp_len;
						}
					}
					if(sum>max_sum){
						secMax_sum = max_sum; secMax_sum_id = max_sum_id; secMax_valid_seg_num = max_valid_seg_num;
						max_sum = sum; max_sum_id = i; max_valid_seg_num = valid_seg_num;
					}else if(sum>secMax_sum){secMax_sum = sum; secMax_sum_id = i; secMax_valid_seg_num = valid_seg_num; }
				}
			}
			best_aln_id = max_sum_id;
			if(secMax_sum!=max_sum) {
				ratio = (double)secMax_sum / max_sum;
				if(ratio>0.8){
					aver_len_max = (double)max_sum / max_valid_seg_num;
					aver_len_secMax = (double)secMax_sum / secMax_valid_seg_num;
					if(aver_len_max<aver_len_secMax)
						best_aln_id = secMax_sum_id;
				}
			}

			blat_aln = blat_aln_vec.at(best_aln_id);
			blat_aln->best_aln = true;
		}
	}
}

// assign align orientation flag
void varCand::assignAlnOrientFlag(){
	int32_t aln_orient = -1, maxQueryLen;
	blat_aln_t *blat_aln;

	maxQueryLen = 0;
	for(size_t i=0; i<blat_aln_vec.size(); i++){
		blat_aln = blat_aln_vec.at(i);
		if(blat_aln->best_aln==true and blat_aln->query_len>maxQueryLen){
			aln_orient = blat_aln->aln_orient;
			maxQueryLen = blat_aln->query_len;
		}
	}
	if(aln_orient!=-1){
		for(size_t i=0; i<varVec.size(); i++)
			varVec[i]->aln_orient = aln_orient;
		for(size_t i=0; i<newVarVec.size(); i++)
			newVarVec[i]->aln_orient = aln_orient;
	}else{
		cerr << "aln_orient=" << aln_orient << endl;
		exit(1);
	}
}

// filter the blat alignments
void varCand::blatFilter(){
	// filter by short align segments
	blatFilterByShortAlnSegs();

	// filter by query coverage
	blatFilterByQueryCoverage();

	// filter by genome coverage
	blatFilterByGenomeCoverage();
}

// filter by short align segments
void varCand::blatFilterByShortAlnSegs(){
	blat_aln_t *blat_aln;
	aln_seg_t *aln_seg;
	int32_t valid_aln_num, valid_seg_num, query_len, subject_len, sum;

	// filter by short align segments
	valid_aln_num = 0;
	for(size_t i=0; i<blat_aln_vec.size(); i++){
		blat_aln = blat_aln_vec[i];
		if(blat_aln->valid_aln==false){
			valid_seg_num = 0; sum = 0;
			for(size_t j=0; j<blat_aln->aln_segs.size(); j++){
				aln_seg = blat_aln->aln_segs[j];
				query_len = aln_seg->query_end - aln_seg->query_start + 1;
				subject_len = aln_seg->subject_end - aln_seg->subject_start + 1;
				sum += subject_len;
				if(query_len>=MIN_VALID_BLAT_SEG_SIZE or subject_len>=MIN_VALID_BLAT_SEG_SIZE) valid_seg_num ++;
			}
			if((double)valid_seg_num/blat_aln->aln_segs.size()>=MIN_VALID_BLAT_SEG_FRATCION or (double)sum/blat_aln->aln_segs.size()>MIN_AVER_SIZE_ALN_SEG) {
				blat_aln->valid_aln = true;
				valid_aln_num ++;
			}
		}
	}

	// process the most reasonable segments
	if(valid_aln_num==0){
		for(size_t i=0; i<blat_aln_vec.size(); i++){
			blat_aln = blat_aln_vec[i];
			if(blat_aln->best_aln==true and blat_aln->valid_aln==false) blat_aln->valid_aln = true;
		}
	}
}

// filter by query coverage
void varCand::blatFilterByQueryCoverage(){
	size_t i, j;
	int32_t start_query_pos, end_query_pos;
	vector<string> query_seq_vec;
	int8_t **query_cov_array;
	blat_aln_t *blat_aln;
	aln_seg_t *aln_seg;
	float covRatio;

	if(ctg_num==0) {
		cout << "ctgfilename=" << ctgfilename << ", ctg_num=" << ctg_num << endl;
		return;
	}

	// allocate array
	query_cov_array = new int8_t *[ctg_num]();
	FastaSeqLoader fa_loader(ctgfilename);
	for(i=0; i<(size_t)ctg_num; i++) query_cov_array[i] = new int8_t[fa_loader.getFastaSeq(i, ALN_PLUS_ORIENT).size()+1]();

//	ctg_id_tmp = -1;
//	for(i=0; i<blat_aln_vec.size(); i++){
//		blat_aln = blat_aln_vec.at(i);
//		//if((int32_t)blat_aln->query_id!=ctg_id_tmp) query_cov_array[++ctg_id_tmp] = new int8_t[blat_aln->query_len+1]();
//		if(ctg_id_tmp>=ctg_num){
//			cerr << "line=" << __LINE__ << ", ctg_id_tmp=" << ctg_id_tmp << ", ctg_num=" << ctg_num << ", error!" << endl;
//			exit(1);
//		}
//	}

	// fill coverage array according to the best blat alignments
	for(i=0; i<blat_aln_vec.size(); i++){
		blat_aln = blat_aln_vec.at(i);
		if(blat_aln->best_aln){
			// compute the query coverage based on the best blat alignment block for each query
			for(j=0; j<blat_aln->aln_segs.size(); j++){
				aln_seg = blat_aln->aln_segs.at(j);
				if(aln_seg->aln_orient==ALN_PLUS_ORIENT){
					start_query_pos = aln_seg->query_start;
					end_query_pos = aln_seg->query_end;
				}else{
					start_query_pos = blat_aln->query_len - aln_seg->query_end + 1;
					end_query_pos = blat_aln->query_len - aln_seg->query_start + 1;
				}
				memset(query_cov_array[blat_aln->query_id]+start_query_pos-1, 1, sizeof(int8_t)*(end_query_pos-start_query_pos+1));
			}
		}
	}

	// compute the coverage for other blat alignment blocks
	for(i=0; i<blat_aln_vec.size(); i++){
		blat_aln = blat_aln_vec.at(i);
		if(blat_aln->valid_aln and blat_aln->best_aln==false){
			covRatio = computeRepeatCovRatio(blat_aln, query_cov_array[blat_aln->query_id], true); // compute the coverage
			if(covRatio>=INVALID_COV_RATIO_THRES) blat_aln->valid_aln = false;
			else updateCovArray(blat_aln, query_cov_array[blat_aln->query_id], true); // update the coverage
		}
	}

	for(i=0; i<(size_t)ctg_num; i++)
		delete[] query_cov_array[i];
	delete[] query_cov_array;

//	size_t i, j, validNum, best_aln_num, start_query_pos, end_query_pos;
//	int8_t **query_cov_array;
//	blat_aln_t *blat_aln;
//	aln_seg_t *aln_seg;
//	float covRatio;
//
//	validNum = 0;
//	for(i=0; i<blat_aln_vec.size(); i++) if(blat_aln_vec[i]->valid_aln) validNum ++;
//	if(validNum>(size_t)ctg_num and ctg_num>0){
//		query_cov_array = new int8_t *[validNum]();
//		// get the number of best blat alignments
//		best_aln_num = 0;
//		for(i=0; i<blat_aln_vec.size(); i++){
//			blat_aln = blat_aln_vec.at(i);
//			if(blat_aln->best_aln){
//				query_cov_array[best_aln_num] = new int8_t[blat_aln->query_len+1]();
//				// compute the query coverage based on the first blat alignment block for each query
//				for(j=0; j<blat_aln->aln_segs.size(); j++){
//					aln_seg = blat_aln->aln_segs.at(j);
//					if(aln_seg->aln_orient==ALN_PLUS_ORIENT){
//						start_query_pos = aln_seg->query_start;
//						end_query_pos = aln_seg->query_end;
//					}else{
//						start_query_pos = blat_aln->query_len - aln_seg->query_end + 1;
//						end_query_pos = blat_aln->query_len - aln_seg->query_start + 1;
//					}
//					memset(query_cov_array[best_aln_num]+start_query_pos-1, 1, sizeof(int8_t)*(end_query_pos-start_query_pos+1));
//				}
//				best_aln_num ++;
//				if(best_aln_num>validNum){
//					cerr << "line=" << __LINE__ << ", best_aln_num=" << best_aln_num << ", validNum=" << validNum << ", error!" << endl;
//					exit(1);
//				}
//			}
//		}
//
//		// compute the coverage of other blat alignment blocks
//		for(i=0; i<blat_aln_vec.size(); i++){
//			blat_aln = blat_aln_vec.at(i);
//			if(blat_aln->valid_aln and blat_aln->best_aln==false){
//				covRatio = computeRepeatCovRatio(blat_aln, query_cov_array[blat_aln->query_id], true); // compute the coverage
//				if(covRatio>=INVALID_COV_RATIO_THRES) blat_aln->valid_aln = false;
//				else updateCovArray(blat_aln, query_cov_array[blat_aln->query_id], true); // update the coverage
//			}
//		}
//
//		//for(i=0; i<(size_t)ctg_num; i++) delete[] query_cov_array[i];
//		for(i=0; i<best_aln_num; i++)
//			delete[] query_cov_array[i];
//		delete[] query_cov_array;
//	}
}

// filter by genome coverage
void varCand::blatFilterByGenomeCoverage(){
	size_t i, j, validNum, subject_len, maxIdx, maxValue, sum;
	int8_t *subject_cov_array;
	blat_aln_t *blat_aln;
	aln_seg_t *aln_seg;
	float covRatio;

	validNum = 0;
	for(i=0; i<blat_aln_vec.size(); i++) if(blat_aln_vec.at(i)->valid_aln) validNum ++;
	if(validNum>=2){
		// get the largest blat alignment block
		maxIdx = 0;
		maxValue = 0;
		for(i=0; i<blat_aln_vec.size(); i++){
			blat_aln = blat_aln_vec.at(i);
			if(blat_aln->valid_aln){
				sum = 0;
				for(j=0; j<blat_aln->aln_segs.size(); j++){
					aln_seg = blat_aln->aln_segs.at(j);
					sum += aln_seg->subject_end - aln_seg->subject_start + 1;
				}
				if(maxValue<sum){
					maxValue = sum;
					maxIdx = i;
				}
			}
		}

		// initialize the genome coverage
		blat_aln = blat_aln_vec.at(maxIdx);
		subject_len = blat_aln->subject_len;
		subject_cov_array = new int8_t[subject_len+1]();
		memset(subject_cov_array, 0, sizeof(int8_t)*subject_len); // reset to 0: not covered

		// compute the genome coverage based on the largest blat alignment block
		for(i=0; i<blat_aln->aln_segs.size(); i++){
			aln_seg = blat_aln->aln_segs.at(i);
			memset(subject_cov_array+aln_seg->subject_start-1, 1, sizeof(int8_t)*(aln_seg->subject_end-aln_seg->subject_start+1));
		}

		// compute the coverage of other blat alignment blocks
		for(i=0; i<blat_aln_vec.size(); i++){
			blat_aln = blat_aln_vec.at(i);
			if(i!=maxIdx and blat_aln->valid_aln){
				covRatio = computeRepeatCovRatio(blat_aln, subject_cov_array, false); // compute the coverage
				if(covRatio>=INVALID_COV_RATIO_THRES) blat_aln->valid_aln = false;
				else updateCovArray(blat_aln, subject_cov_array, false); // update the coverage
			}
		}

		delete[] subject_cov_array;
	}
}

// compute the coverage percent of blat alignment blocks
float varCand::computeRepeatCovRatio(blat_aln_t *blat_aln, int8_t *cov_array, bool query_flag){
	size_t total, repeatSum, startPos, endPos;
	aln_seg_t *aln_seg;

	total = 0, repeatSum = 0;
	for(size_t i=0; i<blat_aln->aln_segs.size(); i++){
		aln_seg = blat_aln->aln_segs[i];
		if(query_flag){
			if(aln_seg->aln_orient==ALN_PLUS_ORIENT){
				startPos = aln_seg->query_start;
				endPos = aln_seg->query_end;
			}else{
				startPos = blat_aln->query_len - aln_seg->query_end + 1;
				endPos = blat_aln->query_len - aln_seg->query_start + 1;
			}
		}else{
			startPos = aln_seg->subject_start;
			endPos = aln_seg->subject_end;
		}
		for(size_t j=startPos-1; j<endPos; j++){
			total ++;
			if(cov_array[j]==1)	repeatSum ++;
		}
	}
	return (float)repeatSum/total;
}

// update the coverage
void varCand::updateCovArray(blat_aln_t *blat_aln, int8_t *cov_array, bool query_flag){
	size_t startPos, endPos;
	aln_seg_t *aln_seg;
	for(size_t i=0; i<blat_aln->aln_segs.size(); i++){
		aln_seg = blat_aln->aln_segs.at(i);
		if(query_flag){
			if(aln_seg->aln_orient==ALN_PLUS_ORIENT){
				startPos = aln_seg->query_start;
				endPos = aln_seg->query_end;
			}else{
				startPos = blat_aln->query_len - aln_seg->query_end + 1;
				endPos = blat_aln->query_len - aln_seg->query_start + 1;
			}
		}
		else{ startPos = aln_seg->subject_start; endPos = aln_seg->subject_end; }
		for(size_t j=startPos-1; j<endPos; j++) if(cov_array[j]==0)	cov_array[j] = 1;
	}
}

// call indel types
void varCand::determineIndelType(){
	int32_t  i, j, k, start_seg_idx, end_seg_idx, query_dist, ref_dist, decrease_size, var_type;
	int64_t startRefPos, endRefPos, startLocalRefPos, endLocalRefPos, startQueryPos, endQueryPos;
	blat_aln_t *blat_aln;
	aln_seg_t *seg1, *seg2;
	reg_t *reg, *reg_tmp;
	string reg_str;
	bool newFlag, updateFlag;
	vector<reg_t*> foundRegVec, candRegVec, reg_vec_tmp;
	//vector<int64_t> ref_query_len_vec; // ref_len, local_ref_len, query_len

	for(i=0; i<(int32_t)blat_aln_vec.size(); i++){
		blat_aln = blat_aln_vec[i];
		if(blat_aln->valid_aln==true and blat_aln->aln_segs.size()>1){
			for(j=0; j<(int32_t)blat_aln->aln_segs.size()-1; j++){
				// get the start and end segment idx
				start_seg_idx = end_seg_idx = -1;
				for(k=j; k<(int32_t)blat_aln->aln_segs.size()-1; k++){
					seg1 = blat_aln->aln_segs[k];
					query_dist = seg1->query_end - seg1->query_start + 1;
					ref_dist = seg1->ref_end - seg1->ref_start + 1;
					if(query_dist<0) query_dist = -query_dist;
					if(ref_dist<0) ref_dist = -ref_dist;
					if(query_dist>=MIN_INNER_BLAT_SEG_SIZE and ref_dist>=MIN_INNER_BLAT_SEG_SIZE){
						start_seg_idx = k;
						break;
					}
				}
				if(start_seg_idx!=-1){
					for(k=start_seg_idx+1; k<(int32_t)blat_aln->aln_segs.size(); k++){
						seg2 = blat_aln->aln_segs[k];
						query_dist = seg2->query_end - seg2->query_start + 1;
						ref_dist = seg2->ref_end - seg2->ref_start + 1;
						if(query_dist<0) query_dist = -query_dist;
						if(ref_dist<0) ref_dist = -ref_dist;
						if(query_dist>=MIN_INNER_BLAT_SEG_SIZE and ref_dist>=MIN_INNER_BLAT_SEG_SIZE){
							end_seg_idx = k;
							break;
						}
					}
				}
				if(start_seg_idx==-1 or end_seg_idx==-1) continue;

				var_type = VAR_UNC;
				seg1 = blat_aln->aln_segs[start_seg_idx];
				seg2 = blat_aln->aln_segs[end_seg_idx];
				if(seg1->query_end>CTG_END_SKIP_SIZE and seg2->query_start>CTG_END_SKIP_SIZE){  // skip contig ends
					query_dist = seg2->query_start - seg1->query_end;
					ref_dist = seg2->ref_start - seg1->ref_end;
					if(ref_dist<0 and query_dist>=0)
						var_type = VAR_DUP;
					else {
						if(query_dist>=ref_dist) var_type = VAR_INS;  // insertion
						else var_type = VAR_DEL;  // deletion
					}
				}

				if(var_type!=VAR_UNC){
					// find the variation and update its information, including variation type, positions, sequences
					if(seg1->ref_end+1<seg2->ref_start){
						startRefPos = seg1->ref_end + 1;
						endRefPos = seg2->ref_start - 1;
						startLocalRefPos = seg1->subject_end + 1;
						endLocalRefPos = seg2->subject_start - 1;
						startQueryPos = seg1->query_end + 1;
						endQueryPos = seg2->query_start - 1;
					}else if(seg1->ref_end+1==seg2->ref_start){
						startRefPos = seg1->ref_end;
						endRefPos = seg2->ref_start - 1;
						startLocalRefPos = seg1->subject_end;
						endLocalRefPos = seg2->subject_start - 1;
						startQueryPos = seg1->query_end;
						endQueryPos = seg2->query_start - 1;
					}else{
						startRefPos = seg1->ref_end;
						endRefPos = seg2->ref_start;
						startLocalRefPos = seg1->subject_end;
						endLocalRefPos = seg2->subject_start;
						startQueryPos = seg1->query_end;
						endQueryPos = seg2->query_start;
					}

					decrease_size = -1;
					ref_dist = endRefPos - startRefPos;
					query_dist = endQueryPos - startQueryPos;
					if(ref_dist<0 or query_dist<0){
						if(ref_dist<query_dist)
							decrease_size = -ref_dist;
						else
							decrease_size = -query_dist;
						startRefPos -= decrease_size;
						startLocalRefPos -= decrease_size;
						startQueryPos -= decrease_size;
					}

					newFlag = false;
					reg_vec_tmp = findVarvecItemAll(startRefPos, endRefPos, varVec);
					if(reg_vec_tmp.size()==0) reg_vec_tmp = findVarvecItemAllExtSize(startRefPos, endRefPos, varVec, EXT_SIZE_CHK_VAR_LOC, EXT_SIZE_CHK_VAR_LOC); // try extend size

					if(reg_vec_tmp.size()>0){ // found, true positive, and then update its information
						for(k=0; k<(int32_t)reg_vec_tmp.size(); k++){
							reg = reg_vec_tmp.at(k);
							if(reg->blat_aln_id==-1)
								reg->blat_aln_id = i;
							if(reg->blat_aln_id==i){
								foundRegVec.push_back(reg);

								reg_tmp = new reg_t();
								reg_tmp->var_type = var_type;
								reg_tmp->chrname = chrname;
								reg_tmp->startRefPos = startRefPos;
								reg_tmp->endRefPos = endRefPos;
								reg_tmp->startLocalRefPos = startLocalRefPos;
								reg_tmp->endLocalRefPos = endLocalRefPos;
								reg_tmp->startQueryPos = startQueryPos;
								reg_tmp->endQueryPos = endQueryPos;
								reg_tmp->query_id = blat_aln->query_id;
								reg_tmp->sv_len = (endQueryPos-startQueryPos)-(endRefPos-startRefPos);
								reg_tmp->blat_aln_id = i;
								reg_tmp->call_success_status = true;
								reg_tmp->short_sv_flag = false;
								reg_tmp->zero_cov_flag = false;
								reg_tmp->aln_seg_end_flag = false;
								reg_tmp->aln_orient = blat_aln->aln_orient;

								candRegVec.push_back(reg_tmp);
							}
						}
					}else { // if it is a false negative discovery, search the newVarVec
						reg = findVarvecItem(startRefPos, endRefPos, newVarVec);
						if(reg==NULL) {  // if not exist, add it to newVarVec
							reg = new reg_t();
							reg->var_type = var_type;
							reg->chrname = chrname;
							reg->startRefPos = startRefPos;
							reg->endRefPos = endRefPos;
							reg->startLocalRefPos = startLocalRefPos;
							reg->endLocalRefPos = endLocalRefPos;
							reg->startQueryPos = startQueryPos;
							reg->endQueryPos = endQueryPos;
							reg->query_id = blat_aln->query_id;
							reg->sv_len = (endQueryPos-startQueryPos)-(endRefPos-startRefPos);
							reg->blat_aln_id = i;
							reg->call_success_status = true;
							reg->short_sv_flag = false;
							reg->zero_cov_flag = false;
							reg->aln_seg_end_flag = false;
							reg->aln_orient = blat_aln->aln_orient;
							newVarVec.push_back(reg);
							newFlag = true;
						}
					}

					// confirm the variant position

				}
			}
		}else{
			// check the small differences
			//cout << var_cand->alnfilename << endl;
		}
	}

	// deal with the new called variants in exist larger variant regions
	processInnerContainedRegs(foundRegVec, candRegVec);

	mergeOverlappedReg(varVec);

	// release candRegVec
	vector<reg_t*>::iterator it;
	for(it=candRegVec.begin(); it!=candRegVec.end(); it++)
		delete (*it);
	vector<reg_t*>().swap(candRegVec);
}

// call short variants that cannot be split at the variant location
void varCand::callShortVariants(){
	reg_t *reg, *reg_overlapped;
	aln_seg_t *aln_seg;
	localAln_t *local_aln;
	int64_t blat_aln_id;

	for(size_t i=0; i<varVec.size(); i++){
		reg = varVec[i];
		if(reg->call_success_status==false){
			reg_overlapped = getOverlappedRegByCallFlag(reg, varVec);
			if(reg_overlapped==NULL and isUnsplitAln(reg)){
				reg->short_sv_flag = true;
				if(determineQueryidReg(reg)){ // determine query_id
					// get the aln_seg
					aln_seg = getAlnSeg(reg, &blat_aln_id);
					if(aln_seg){
						// compute local locations
						local_aln = new localAln_t();
						local_aln->reg = reg;
						local_aln->aln_seg = aln_seg;
						local_aln->blat_aln_id = blat_aln_id;
						local_aln->start_seg_extend = local_aln->end_seg_extend = NULL;
						local_aln->startRefPos = local_aln->endRefPos = -1;
						local_aln->startLocalRefPos = local_aln->startQueryPos = local_aln->endLocalRefPos = local_aln->endQueryPos = -1;
						local_aln->queryLeftShiftLen = local_aln->queryRightShiftLen = local_aln->localRefLeftShiftLen = local_aln->localRefRightShiftLen = -1;
						local_aln->chrlen = faidx_seq_len(fai, local_aln->reg->chrname.c_str()); // get the reference length

						// extend the align segments if the given segment is short (<1kb)
						if(aln_seg->ref_end-aln_seg->ref_start<MIN_AVER_SIZE_ALN_SEG)
							computeExtendAlnSegs(local_aln);

						computeLocalLocsAlnShortVar(local_aln);

						// get sub local reference sequence
						FastaSeqLoader refseqloader(refseqfilename);
						local_aln->refseq = refseqloader.getFastaSeqByPos(0, local_aln->startLocalRefPos, local_aln->endLocalRefPos, ALN_PLUS_ORIENT);

						//cout << local_aln->refseq << endl;

						FastaSeqLoader ctgseqloader(ctgfilename);
						local_aln->ctgseq = ctgseqloader.getFastaSeqByPos(reg->query_id, local_aln->startQueryPos, local_aln->endQueryPos, reg->aln_orient);
						//local_aln->ctgseq = ctgseqloader.getFastaSeqByPos(reg->query_id, local_aln->startQueryPos, local_aln->endQueryPos, ALN_PLUS_ORIENT);

						//cout << local_aln->ctgseq << endl;

						// pairwise alignment
						upperSeq(local_aln->ctgseq);
						upperSeq(local_aln->refseq);
						computeSeqAlignment(local_aln);

						// compute variant location
						computeVarLoc(local_aln);

						// adjust the short variant locations according to alignment
						adjustVarLocShortVar(local_aln);

						//confirmVarLoc(local_aln);

						// confirm short variation
						confirmShortVar(local_aln);

						delete local_aln;
					}
				}
			}
		}
	}
}

// check whether the variant region is unsplit
bool varCand::isUnsplitAln(reg_t *reg){
	blat_aln_t *blat_aln;
	aln_seg_t *aln_seg;
	bool flag = false;

	for(size_t i=0; i<blat_aln_vec.size(); i++){
		blat_aln = blat_aln_vec[i];
		if(blat_aln->valid_aln==true){
			for(size_t j=0; j<blat_aln->aln_segs.size(); j++){
				aln_seg = blat_aln->aln_segs[j];
				if(reg->startRefPos>aln_seg->ref_start and reg->endRefPos<aln_seg->ref_end){
					flag = true;
					break;
				}
			}
			if(flag) break;
		}
	}

	return flag;
}

// determine query_id
bool varCand::determineQueryidReg(reg_t *reg){
	blat_aln_t *blat_aln;
	aln_seg_t *aln_seg;
	bool flag = false;

	reg->query_id = -1;
	if(ctg_num==1){
		flag = true;
		reg->query_id = 0;
	}else if(ctg_num>=2){
		for(size_t i=0; i<varVec.size(); i++){
			if(reg==varVec[i]){
				if((i>0 and i+1<varVec.size()) and varVec[i-1]->query_id==varVec[i+1]->query_id){
					flag = true;
					reg->query_id = varVec[i-1]->query_id;
					break;
				}
			}
		}
		if(flag==false){
			for(size_t i=0; i<blat_aln_vec.size(); i++){
				blat_aln = blat_aln_vec[i];
				if(blat_aln->valid_aln==true){
					for(size_t j=0; j<blat_aln->aln_segs.size(); j++){
						aln_seg = blat_aln->aln_segs[j];
						if(reg->startRefPos>aln_seg->ref_start and reg->endRefPos<aln_seg->ref_end){
							reg->query_id = blat_aln->query_id;
							flag = true;
							break;
						}
					}
					if(flag) break;
				}
			}
		}
	}
	return flag;
}

// get the aln_seg
aln_seg_t* varCand::getAlnSeg(reg_t *reg, int64_t *blat_aln_id){
	blat_aln_t *blat_aln;
	aln_seg_t *aln_seg, *target;

	*blat_aln_id = -1;
	target = NULL;
	if(reg->query_id!=-1){
		for(size_t i=0; i<blat_aln_vec.size(); i++){
			blat_aln = blat_aln_vec.at(i);
			if(blat_aln->valid_aln==true and blat_aln->query_id==reg->query_id){
				for(size_t j=0; j<blat_aln->aln_segs.size(); j++){
					aln_seg = blat_aln->aln_segs[j];
					if(reg->startRefPos>aln_seg->ref_start and reg->endRefPos<aln_seg->ref_end){
						target = aln_seg;
						*blat_aln_id = i;
						break;
					}
				}
			}
			if(target) break;
		}
	}
	return target;
}

// extend the align segments if the seed segment is short (<1kb)
void varCand::computeExtendAlnSegs(localAln_t *local_aln){
	size_t i, seglen1, seglen2;
	blat_aln_t *blat_aln;
	aln_seg_t *aln_seg, *aln_seg1, *aln_seg2, *start_seg_extend, *end_seg_extend;
	int64_t startRefPos_extend, endRefPos_extend, chrlen_tmp, dist;

	start_seg_extend = end_seg_extend = NULL;
	if(local_aln->blat_aln_id!=-1 and local_aln->aln_seg and local_aln->reg){
		blat_aln = blat_aln_vec.at(local_aln->blat_aln_id);
		aln_seg = local_aln->aln_seg;
		chrlen_tmp = local_aln->chrlen;

		startRefPos_extend = aln_seg->ref_start - VAR_ALN_EXTEND_SIZE;
		endRefPos_extend = aln_seg->ref_end + VAR_ALN_EXTEND_SIZE;
		if(startRefPos_extend<1) startRefPos_extend = 1;
		if(endRefPos_extend>chrlen_tmp) endRefPos_extend = chrlen_tmp;

		if(aln_seg->ref_end-aln_seg->ref_start<MIN_AVER_SIZE_ALN_SEG){
			// extend on both sides
			if(blat_aln->aln_segs.at(0)->ref_start+VAR_ALN_EXTEND_SIZE>=aln_seg->ref_start)
				start_seg_extend = blat_aln->aln_segs.at(0);
			else{
				for(i=0; i<blat_aln->aln_segs.size(); i++){
					aln_seg1 = blat_aln->aln_segs.at(i);
					aln_seg2 = blat_aln->aln_segs.at(i+1);
					seglen1 = aln_seg1->ref_end - aln_seg1->ref_start + 1;
					seglen2 = aln_seg2->ref_end - aln_seg2->ref_start + 1;
					dist = aln_seg2->ref_start - aln_seg1->ref_end;
					if(seglen1>=MIN_AVER_SIZE_ALN_SEG and seglen2>=MIN_AVER_SIZE_ALN_SEG and dist<=MIN_AVER_SIZE_ALN_SEG){
						if((startRefPos_extend>=aln_seg1->ref_start and startRefPos_extend<=aln_seg1->ref_end) or (startRefPos_extend>=aln_seg1->ref_end and startRefPos_extend<=aln_seg2->ref_start)){
							start_seg_extend = aln_seg1;
							break;
						}else if(startRefPos_extend>=aln_seg2->ref_start and startRefPos_extend<=aln_seg2->ref_end){
							start_seg_extend = aln_seg2;
							break;
						}
					}
					if(aln_seg2==aln_seg) break;
				}
			}
			if(aln_seg->ref_end+VAR_ALN_EXTEND_SIZE>=blat_aln->aln_segs.at(blat_aln->aln_segs.size()-1)->ref_end)
				end_seg_extend = blat_aln->aln_segs.at(blat_aln->aln_segs.size()-1);
			else{
				for(i=blat_aln->aln_segs.size()-1; i>0; i--){
					aln_seg1 = blat_aln->aln_segs.at(i-1);
					aln_seg2 = blat_aln->aln_segs.at(i);
					seglen1 = aln_seg1->ref_end - aln_seg1->ref_start + 1;
					seglen2 = aln_seg2->ref_end - aln_seg2->ref_start + 1;
					dist = aln_seg2->ref_start - aln_seg1->ref_end;
					if(seglen1>=MIN_AVER_SIZE_ALN_SEG and seglen2>=MIN_AVER_SIZE_ALN_SEG and dist<=MIN_AVER_SIZE_ALN_SEG){
						if(endRefPos_extend>=aln_seg1->ref_start and endRefPos_extend<=aln_seg1->ref_end){
							end_seg_extend = aln_seg1;
							break;
						}else if((endRefPos_extend>=aln_seg1->ref_end and endRefPos_extend<=aln_seg2->ref_start) or (endRefPos_extend>=aln_seg2->ref_start and endRefPos_extend<=aln_seg2->ref_end)){
							end_seg_extend = aln_seg2;
							break;
						}
					}
					if(aln_seg1==aln_seg) break;
				}
			}
		}
	}

	local_aln->start_seg_extend = start_seg_extend;
	local_aln->end_seg_extend = end_seg_extend;
}

// compute local locations
void varCand::computeLocalLocsAlnShortVar(localAln_t *local_aln){
	reg_t *reg = local_aln->reg;
	aln_seg_t *aln_seg, *start_seg_extend, *end_seg_extend;;
	int32_t left_dist, right_dist; //, querylen; //, tmp_pos;
	blat_aln_t *blat_aln;

	aln_seg = local_aln->aln_seg;
	start_seg_extend = local_aln->start_seg_extend;
	end_seg_extend = local_aln->end_seg_extend;
	blat_aln = blat_aln_vec.at(local_aln->blat_aln_id);

	if(start_seg_extend==NULL or end_seg_extend==NULL){ // not two extensions on both sides
		if(aln_seg->ref_start+VAR_ALN_EXTEND_SIZE<reg->startRefPos)
			left_dist = reg->startRefPos - aln_seg->ref_start - VAR_ALN_EXTEND_SIZE;
		else
			left_dist = 0;
		local_aln->startRefPos = aln_seg->ref_start + left_dist;
		local_aln->startLocalRefPos = aln_seg->subject_start + left_dist;
		local_aln->startQueryPos = aln_seg->query_start + left_dist;
		if(local_aln->startLocalRefPos<1) local_aln->startLocalRefPos = 1;
		if(local_aln->startQueryPos<1) local_aln->startQueryPos = 1;

		if(reg->endRefPos+VAR_ALN_EXTEND_SIZE<aln_seg->ref_end)
			right_dist = aln_seg->ref_end - reg->endRefPos - VAR_ALN_EXTEND_SIZE;
		else
			right_dist = 0;
		local_aln->endRefPos = aln_seg->ref_end - right_dist;
		local_aln->endLocalRefPos = aln_seg->subject_end - right_dist;
		local_aln->endQueryPos = aln_seg->query_end - right_dist;
		if(local_aln->endLocalRefPos>blat_aln->subject_len) local_aln->endLocalRefPos = blat_aln->subject_len;
		if(local_aln->endQueryPos>blat_aln->query_len) local_aln->endQueryPos = blat_aln->query_len;

//		if(reg->aln_orient==ALN_MINUS_ORIENT) {
//			FastaSeqLoader ctgseqloader(ctgfilename);
//			querylen = ctgseqloader.getFastaSeqLen(reg->query_id);
//			tmp_pos = local_aln->startQueryPos;
//			local_aln->startQueryPos = querylen - local_aln->endQueryPos + 1;
//			local_aln->endQueryPos = querylen - tmp_pos + 1;
//		}
	}else{ // two extensions on both sides
		if(start_seg_extend->ref_start+VAR_ALN_EXTEND_SIZE<aln_seg->ref_start)
			left_dist = aln_seg->ref_start - start_seg_extend->ref_start - VAR_ALN_EXTEND_SIZE;
		else
			left_dist = 0;
		local_aln->startRefPos = start_seg_extend->ref_start + left_dist;
		local_aln->startLocalRefPos = start_seg_extend->subject_start + left_dist;
		local_aln->startQueryPos = start_seg_extend->query_start + left_dist;
		if(local_aln->startLocalRefPos<1) local_aln->startLocalRefPos = 1;
		if(local_aln->startQueryPos<1) local_aln->startQueryPos = 1;

		if(aln_seg->ref_end+VAR_ALN_EXTEND_SIZE<end_seg_extend->ref_end)
			right_dist = end_seg_extend->ref_end - aln_seg->ref_end - VAR_ALN_EXTEND_SIZE;
		else
			right_dist = 0;
		local_aln->endRefPos = end_seg_extend->ref_end - right_dist;
		local_aln->endLocalRefPos = end_seg_extend->subject_end - right_dist;
		local_aln->endQueryPos = end_seg_extend->query_end - right_dist;
		if(local_aln->endLocalRefPos>blat_aln->subject_len) local_aln->endLocalRefPos = blat_aln->subject_len;
		if(local_aln->endQueryPos>blat_aln->query_len) local_aln->endQueryPos = blat_aln->query_len;
	}
}

// compute sequence alignment (LCS), all characters are in upper case
void varCand::computeSeqAlignment(localAln_t *local_aln){
	int64_t rowsNum, colsNum, arrSize, mem_cost; // memory cost is measured in kB (1024 bytes)
	localAln_t *local_aln_tmp = NULL;
	bool exist_flag, aln_flag;

	// check local alignment vector
	exist_flag = aln_flag = false;
	if(isLocalAlnInfoComplete(local_aln)){
		local_aln_tmp = getIdenticalLocalAlnItem(local_aln, local_aln_vec);
		if(local_aln_tmp)
			exist_flag = true;
	}

	if(exist_flag){ // already exist, then copy the alignment information
		if(local_aln_tmp){
			copyLocalAlnInInfo(local_aln, local_aln_tmp);
		}else{
			cerr << __func__ << ": local_aln_tmp=" << local_aln_tmp << ", invalid." << endl;
			exit(1);
		}
	}else{ // not exist, then prepare to compute the alignment information
		// check memory consumption
		rowsNum = local_aln->ctgseq.size() + 1;
		colsNum = local_aln->refseq.size() + 1;
		arrSize = rowsNum * colsNum;
		mem_cost = (arrSize * sizeof(struct alnScoreNode)) >> 10; // divide by 1024

		while(1){
			pthread_mutex_lock(&mutex_mem);
			if(mem_seqAln+mem_cost<=mem_total*mem_use_block_factor+extend_total*extend_use_block_factor){ // prepare for alignment computation
				mem_seqAln += mem_cost;
				work_num ++;
				aln_flag = true;
			}
			pthread_mutex_unlock(&mutex_mem);

			if(aln_flag==false){ // block the alignment computation
				//cout << "\t" << __func__ << ": wait " << mem_wait_seconds << " seconds, rowsNum=" << rowsNum << ", colsNum=" << colsNum << ", " << alnfilename << endl;
				sleep(mem_wait_seconds);
			}else break;
		}
	}

	if(aln_flag){
		computeSeqAlignmentOp(local_aln);

		// update memory consumption
		pthread_mutex_lock(&mutex_mem);
		mem_seqAln -= mem_cost;
		work_num --;
		if(mem_seqAln<0 or work_num<0){
			cerr << "line=" << __LINE__ << ", mem_seqAln=" << mem_seqAln << ", work_num=" << work_num << ", error." << endl;
			exit(1);
		}
		pthread_mutex_unlock(&mutex_mem);

		// save local alignment information to vector
		local_aln_tmp = generateNewLocalAlnItem_OnlyAlnInfo(local_aln);
		addLocalAlnItemToVec(local_aln_tmp, local_aln_vec);
	}
}

// worker of compute sequence alignment (LCS), all characters are in upper case
void varCand::computeSeqAlignmentOp(localAln_t *local_aln){
	int64_t i, j, matchScore, mismatchScore, gapScore, gapOpenScore, tmp_gapScore1, tmp_gapScore2, maxValue, scoreIJ;
	int64_t rowsNum = local_aln->ctgseq.size() + 1, colsNum = local_aln->refseq.size() + 1, arrSize;
	struct alnScoreNode *scoreArr;
	int8_t path_val;
	int32_t maxValueLastRow, maxValueLastCol, maxCol, maxRow;
	int32_t mismatchNum, itemNum;
	string queryAlnResult, midAlnResult, subjectAlnResult;
	bool baseMatchFlag;

	matchScore = MATCH_SCORE;
	mismatchScore = MISMATCH_SCORE;
	gapScore = GAP_SCORE;
	gapOpenScore = GAP_OPEN_SCORE;

//	if(rowsNum>=30000 or colsNum>=30000){
		//cout << "\t" << __func__ << ": work_num=" << work_num << ", rowsNum=" << rowsNum << ", colsNum=" << colsNum << ", " << alnfilename << endl;
//	}

	arrSize = rowsNum * colsNum;
	scoreArr = (struct alnScoreNode*) calloc (arrSize, sizeof(struct alnScoreNode));
	if(scoreArr==NULL){
		cerr << "line=" << __LINE__ << ", rowsNum=" << rowsNum << ", colsNum=" << colsNum << ", cannot allocate memory, error!" << endl;
		exit(1);
	}

	for(i=0; i<rowsNum; i++) for(j=0; j<colsNum; j++) scoreArr[i*colsNum+j].path_val = -1;

	// compute the scores of each element
	for(i=1; i<rowsNum; i++){
		for(j=1; j<colsNum; j++){
			baseMatchFlag = isBaseMatch(local_aln->ctgseq[i-1], local_aln->refseq[j-1]);
			if(baseMatchFlag) scoreIJ = matchScore;
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

	// get the row and col of the maximal element in last row and column
	maxValueLastRow = INT_MIN; maxRow = maxCol = -1;
	for(j=0; j<colsNum; j++)
		if(scoreArr[(rowsNum-1)*colsNum+j].score>maxValueLastRow){
			maxValueLastRow = scoreArr[(rowsNum-1)*colsNum+j].score;
			maxCol = j;
		}
	maxValueLastCol = INT_MIN;
	for(i=0; i<rowsNum; i++)
		if(scoreArr[i*colsNum+colsNum-1].score>maxValueLastCol){
			maxValueLastCol = scoreArr[i*colsNum+colsNum-1].score;
			maxRow = i;
		}
	if(maxValueLastRow>=maxValueLastCol){
		maxValue = maxValueLastRow;
		maxRow = rowsNum - 1;
	}else{
		maxValue = maxValueLastCol;
		maxCol = colsNum - 1;
	}

	queryAlnResult = midAlnResult = subjectAlnResult = "";
	mismatchNum = 0;
	itemNum = 0;
	i = maxRow;
	j = maxCol;
	while(i>0 && j>0){
		if(scoreArr[i*colsNum+j].path_val==0){ // from (i-1, j-1)
			queryAlnResult += local_aln->ctgseq[i-1];
			baseMatchFlag = isBaseMatch(local_aln->ctgseq[i-1], local_aln->refseq[j-1]);
			if(baseMatchFlag) midAlnResult += '|';
			else{
				mismatchNum ++;
				midAlnResult += ' ';
			}

			subjectAlnResult += local_aln->refseq[j-1];
			i --;
			j --;
		}else if(scoreArr[i*colsNum+j].path_val==1){ // from (i-1, j)
			queryAlnResult += local_aln->ctgseq[i-1];
			midAlnResult += ' ';
			subjectAlnResult += '-';
			mismatchNum ++;
			i --;
		}else { // from (i, j-1)
			queryAlnResult += '-';
			midAlnResult += ' ';
			subjectAlnResult += local_aln->refseq[j-1];
			mismatchNum ++;
			j --;
		}
		itemNum ++;
	}

	local_aln->queryLeftShiftLen = i;
	local_aln->localRefLeftShiftLen = j;
	local_aln->queryRightShiftLen = rowsNum - 1 - maxRow;
	local_aln->localRefRightShiftLen = colsNum - 1 - maxCol;
	local_aln->overlapLen = itemNum;

	// recover the alignment result
	reverseSeq(queryAlnResult);
	reverseSeq(midAlnResult);
	reverseSeq(subjectAlnResult);

	queryAlnResult.shrink_to_fit();
	midAlnResult.shrink_to_fit();
	subjectAlnResult.shrink_to_fit();

	local_aln->alignResultVec.push_back(queryAlnResult);
	local_aln->alignResultVec.push_back(midAlnResult);
	local_aln->alignResultVec.push_back(subjectAlnResult);

#if ALIGN_DEBUG
	cout << "Query: " << local_aln->alignResultVec[0] << endl;
	cout << "       " << local_aln->alignResultVec[1] << endl;
	cout << "Sbjct: " << local_aln->alignResultVec[2] << endl;
#endif

	// adjust alignment
	adjustAlnResult(local_aln);

#if ALIGN_DEBUG
	cout << "Query: " << local_aln->alignResultVec[0] << endl;
	cout << "       " << local_aln->alignResultVec[1] << endl;
	cout << "Sbjct: " << local_aln->alignResultVec[2] << endl;
#endif

	free(scoreArr);
}

// adjust alignment
void varCand::adjustAlnResult(localAln_t *local_aln){
	int32_t i, j, k, gap_size, sub_seqlen, sub_seqlen2, pos, dif;
	size_t match_pos;
	string queryseq, refseq, sub_seq, sub_seq2;
	bool baseMatchFlag, baseMatchFlag2, found_flag;
	char ch_tmp;
	vector<mismatchReg_t*> misReg_vec;
	mismatchReg_t *mis_reg, *mis_reg2;

	string &query_aln_seq = local_aln->alignResultVec.at(0);
	string &mid_aln_seq = local_aln->alignResultVec.at(1);
	string &subject_aln_seq = local_aln->alignResultVec.at(2);

	// align to left side
	gap_size = 0;
	for(i=0; i<(int32_t)mid_aln_seq.size(); i++){
		if(mid_aln_seq.at(i)!='|'){ // gap or mismatch
			if(gap_size==0)
				gap_size = 1;

			if(subject_aln_seq.at(i)=='-'){ // gap in subject
				for(j=i+gap_size; j<(int32_t)subject_aln_seq.size(); j++){
					if(subject_aln_seq.at(j)!='-'){
						baseMatchFlag = isBaseMatch(query_aln_seq.at(i), subject_aln_seq.at(j));
						if(baseMatchFlag){ // exchange subject_aln_seq[i] and subject_aln_seq[j]
							ch_tmp = subject_aln_seq.at(i);
							subject_aln_seq.at(i) = subject_aln_seq.at(j);
							subject_aln_seq.at(j) = ch_tmp;

							mid_aln_seq.at(i) = '|';
							baseMatchFlag2 = isBaseMatch(query_aln_seq.at(j), subject_aln_seq.at(j));
							if(baseMatchFlag2) mid_aln_seq.at(j) = '|';
							else mid_aln_seq.at(j) = ' ';
						}else{ // gap end
							i += gap_size -1;
							gap_size = 0;
						}
						break;
					}else
						gap_size ++;
				}
			}else if(query_aln_seq.at(i)=='-'){  // gap in query
				for(j=i+gap_size; j<(int32_t)query_aln_seq.size(); j++){
					if(query_aln_seq.at(j)!='-'){
						baseMatchFlag = isBaseMatch(query_aln_seq.at(j), subject_aln_seq.at(i));
						if(baseMatchFlag){ // exchange query_aln_seq[i] and query_aln_seq[j]
							ch_tmp = query_aln_seq.at(i);
							query_aln_seq.at(i) = query_aln_seq.at(j);
							query_aln_seq.at(j) = ch_tmp;

							mid_aln_seq.at(i) = '|';
							baseMatchFlag2 = isBaseMatch(query_aln_seq.at(j), subject_aln_seq.at(j));
							if(baseMatchFlag2) mid_aln_seq.at(j) = '|';
							else mid_aln_seq.at(j) = ' ';
						}else{ // gap end
							i += gap_size -1;
							gap_size = 0;
						}
						break;
					}else
						gap_size ++;
				}
			}else{ // mismatch
				gap_size = 0;
			}
		}
	}

	// align to right side
	gap_size = 0;
	for(i=(int32_t)mid_aln_seq.size()-1; i>=0; i--){
		if(mid_aln_seq.at(i)!='|'){ // gap or mismatch
			if(gap_size==0)
				gap_size = 1;

			if(subject_aln_seq.at(i)=='-'){ // gap in subject
				for(j=i-gap_size; j>=0; j--){
					if(subject_aln_seq.at(j)!='-'){
						baseMatchFlag = isBaseMatch(query_aln_seq.at(i), subject_aln_seq.at(j));
						if(baseMatchFlag){ // exchange subject_aln_seq[i] and subject_aln_seq[j]
							ch_tmp = subject_aln_seq.at(i);
							subject_aln_seq.at(i) = subject_aln_seq.at(j);
							subject_aln_seq.at(j) = ch_tmp;

							mid_aln_seq.at(i) = '|';
							baseMatchFlag2 = isBaseMatch(query_aln_seq.at(j), subject_aln_seq.at(j));
							if(baseMatchFlag2) mid_aln_seq.at(j) = '|';
							else mid_aln_seq.at(j) = ' ';
						}else{ // gap end
							i -= gap_size -1;
							gap_size = 0;
						}
						break;
					}else
						gap_size ++;
				}
			}else if(query_aln_seq.at(i)=='-'){  // gap in query
				for(j=i-gap_size; j>=0; j--){
					if(query_aln_seq.at(j)!='-'){
						baseMatchFlag = isBaseMatch(query_aln_seq.at(j), subject_aln_seq.at(i));
						if(baseMatchFlag){ // exchange query_aln_seq[i] and query_aln_seq[j]
							ch_tmp = query_aln_seq.at(i);
							query_aln_seq.at(i) = query_aln_seq.at(j);
							query_aln_seq.at(j) = ch_tmp;

							mid_aln_seq.at(i) = '|';
							baseMatchFlag2 = isBaseMatch(query_aln_seq.at(j), subject_aln_seq.at(j));
							if(baseMatchFlag2) mid_aln_seq.at(j) = '|';
							else mid_aln_seq.at(j) = ' ';
						}else{ // gap end
							i -= gap_size - 1;
							gap_size = 0;
						}
						break;
					}else
						gap_size ++;
				}
			}else{ // mismatch
				gap_size = 0;
			}
		}
	}

	// adjust the inner matched sequence
	misReg_vec = getMismatchRegVecWithoutPos(local_aln);
	for(i=1; i<(int32_t)misReg_vec.size(); i++){
		mis_reg = misReg_vec.at(i-1);
		mis_reg2 = misReg_vec.at(i);

		if(mis_reg->gap_flag and mis_reg2->gap_flag and mis_reg2->start_aln_idx-mis_reg->end_aln_idx<=20 and (mis_reg->end_aln_idx-mis_reg->start_aln_idx>=10 or mis_reg2->end_aln_idx-mis_reg2->start_aln_idx>=10)){
//			cout << i << ": before adjust ..." << endl;
//
//			cout << "Query: " << local_aln->alignResultVec[0] << endl;
//			cout << "       " << local_aln->alignResultVec[1] << endl;
//			cout << "Sbjct: " << local_aln->alignResultVec[2] << endl;

			// match to left side
			j = mis_reg->start_aln_idx;
			k = mis_reg->end_aln_idx + 1;
			while(j<=mis_reg->end_aln_idx and k<=mis_reg2->start_aln_idx-1){
				if(query_aln_seq.at(j)=='-'){ // gap in query
					if(subject_aln_seq.at(j)==query_aln_seq.at(k)){ // swap and update misReg
						ch_tmp = query_aln_seq.at(j); query_aln_seq.at(j) = query_aln_seq.at(k); query_aln_seq.at(k) = ch_tmp;
						ch_tmp = mid_aln_seq.at(j); mid_aln_seq.at(j) = mid_aln_seq.at(k); mid_aln_seq.at(k) = ch_tmp;

						mis_reg->start_aln_idx ++; mis_reg->end_aln_idx ++;
					}else break;
				}else if(subject_aln_seq.at(j)=='-'){ // gap in reference
					if(query_aln_seq.at(j)==subject_aln_seq.at(k)){ // swap and update misReg
						ch_tmp = subject_aln_seq.at(j); subject_aln_seq.at(j) = subject_aln_seq.at(k); subject_aln_seq.at(k) = ch_tmp;
						ch_tmp = mid_aln_seq.at(j); mid_aln_seq.at(j) = mid_aln_seq.at(k); mid_aln_seq.at(k) = ch_tmp;

						mis_reg->start_aln_idx ++; mis_reg->end_aln_idx ++;
					}else break;
				}else break;
				j++; k++;
			}

			// search the right most match position
			found_flag = false;
			sub_seqlen = mis_reg2->start_aln_idx - mis_reg->end_aln_idx - 1;
			sub_seq = query_aln_seq.substr(mis_reg->end_aln_idx+1, sub_seqlen);
			sub_seqlen2 = mis_reg2->reg_size;
			sub_seq2 = subject_aln_seq.substr(mis_reg2->start_aln_idx, sub_seqlen2);
			match_pos = sub_seq2.rfind(sub_seq);
			if(match_pos!=string::npos){ // found in reference
				if(match_pos>=(size_t)sub_seqlen and mis_reg->reg_size>=sub_seqlen2-(int32_t)match_pos-sub_seqlen){
					found_flag = true;
					match_pos += mis_reg2->start_aln_idx;
					for(j=0; j<sub_seqlen; j++){
						pos = mis_reg->end_aln_idx + 1 + j;
						ch_tmp = query_aln_seq.at(pos); query_aln_seq.at(pos) = query_aln_seq.at(match_pos+j); query_aln_seq.at(match_pos+j) = ch_tmp;
						ch_tmp = mid_aln_seq.at(pos); mid_aln_seq.at(pos) = mid_aln_seq.at(match_pos+j); mid_aln_seq.at(match_pos+j) = ch_tmp;
					}
					// update misReg
					dif = match_pos - mis_reg->end_aln_idx - 1;
					mis_reg->end_aln_idx += dif;
					mis_reg->reg_size += dif;
					mis_reg2->start_aln_idx += dif;
					mis_reg2->reg_size -= dif;
				}
			}else{ // search query
				sub_seq2 = query_aln_seq.substr(mis_reg2->start_aln_idx, sub_seqlen2);
				match_pos = sub_seq2.rfind(sub_seq);
				if(match_pos!=string::npos){ // found in query
					if(match_pos>=(size_t)sub_seqlen and mis_reg->reg_size>=sub_seqlen2-(int32_t)match_pos-sub_seqlen){
						found_flag = true;
						match_pos += mis_reg2->start_aln_idx;
						for(j=0; j<sub_seqlen; j++){
							pos = mis_reg->end_aln_idx + 1 + j;
							ch_tmp = subject_aln_seq.at(pos); subject_aln_seq.at(pos) = subject_aln_seq.at(match_pos+j); subject_aln_seq.at(match_pos+j) = ch_tmp;
							ch_tmp = mid_aln_seq.at(pos); mid_aln_seq.at(pos) = mid_aln_seq.at(match_pos+j); mid_aln_seq.at(match_pos+j) = ch_tmp;
						}
						// update misReg
						dif = match_pos - mis_reg->end_aln_idx - 1;
						mis_reg->end_aln_idx += dif;
						mis_reg->reg_size += dif;
						mis_reg2->start_aln_idx += dif;
						mis_reg2->reg_size -= dif;
					}
				}
			}

			// search the left most match position
			if(found_flag==false){
				sub_seqlen2 = mis_reg->reg_size;
				sub_seq2 = subject_aln_seq.substr(mis_reg->start_aln_idx, sub_seqlen2);
				match_pos = sub_seq2.find(sub_seq);
				if(match_pos!=string::npos){ // found in reference
					if(match_pos>=(size_t)sub_seqlen and match_pos<(size_t)mis_reg2->reg_size){
						found_flag = true;
						match_pos += mis_reg->start_aln_idx;
						for(j=0; j<sub_seqlen; j++){
							pos = mis_reg->end_aln_idx + 1 + j;
							ch_tmp = query_aln_seq.at(pos); query_aln_seq.at(pos) = query_aln_seq.at(match_pos+j); query_aln_seq.at(match_pos+j) = ch_tmp;
							ch_tmp = mid_aln_seq.at(pos); mid_aln_seq.at(pos) = mid_aln_seq.at(match_pos+j); mid_aln_seq.at(match_pos+j) = ch_tmp;
						}
						// update misReg
						dif = mis_reg->end_aln_idx + 1 - match_pos;
						mis_reg->end_aln_idx -= dif;
						mis_reg->reg_size -= dif;
						mis_reg2->start_aln_idx -= dif;
						mis_reg2->reg_size += dif;
					}
				}else{ // search query
					sub_seq2 = query_aln_seq.substr(mis_reg->start_aln_idx, sub_seqlen2);
					match_pos = sub_seq2.find(sub_seq);
					if(match_pos!=string::npos){ // found in query
						if(match_pos>=(size_t)sub_seqlen and match_pos<(size_t)mis_reg2->reg_size){
							found_flag = true;
							match_pos += mis_reg->start_aln_idx;
							for(j=0; j<sub_seqlen; j++){
								pos = mis_reg->end_aln_idx + 1 + j;
								ch_tmp = subject_aln_seq.at(pos); subject_aln_seq.at(pos) = subject_aln_seq.at(match_pos+j); subject_aln_seq.at(match_pos+j) = ch_tmp;
								ch_tmp = mid_aln_seq.at(pos); mid_aln_seq.at(pos) = mid_aln_seq.at(match_pos+j); mid_aln_seq.at(match_pos+j) = ch_tmp;
							}
							// update misReg
							dif = mis_reg->end_aln_idx + 1 - match_pos;
							mis_reg->end_aln_idx -= dif;
							mis_reg->reg_size -= dif;
							mis_reg2->start_aln_idx -= dif;
							mis_reg2->reg_size += dif;
						}
					}
				}
			}

//			cout << ": after adjust ..." << endl;
//
//			cout << "Query: " << local_aln->alignResultVec[0] << endl;
//			cout << "       " << local_aln->alignResultVec[1] << endl;
//			cout << "Sbjct: " << local_aln->alignResultVec[2] << endl;
		}
	}

	// release memory
	releaseMismatchRegVec(misReg_vec);
}

// compute variant location
void varCand::computeVarLoc(localAln_t *local_aln){
	vector<string> var_loc_vec;
	string ctgseq_aln, midseq_aln, refseq_aln;
	int32_t queryPos, refPos, localRefPos, start_aln_idx_var, end_aln_idx_var, start_aln_idx_var_new, end_aln_idx_var_new;
	int32_t i, start_var_pos, end_var_pos, tmp_len;
	int32_t startRefPos_aln, endRefPos_aln, startLocalRefPos_aln, endLocalRefPos_aln, startQueryPos_aln, endQueryPos_aln;
	int32_t startRefPos_aln_new, endRefPos_aln_new, startLocalRefPos_aln_new, endLocalRefPos_aln_new, startQueryPos_aln_new, endQueryPos_aln_new;

	ctgseq_aln = local_aln->alignResultVec[0];
	midseq_aln = local_aln->alignResultVec[1];
	refseq_aln = local_aln->alignResultVec[2];

	start_var_pos = local_aln->reg->startRefPos;
	end_var_pos = local_aln->reg->endRefPos;
	refPos = local_aln->startRefPos + local_aln->localRefLeftShiftLen;
	localRefPos = local_aln->startLocalRefPos + local_aln->localRefLeftShiftLen;
	queryPos = local_aln->startQueryPos + local_aln->queryLeftShiftLen;

	start_aln_idx_var = end_aln_idx_var = -1;
	startRefPos_aln = endRefPos_aln = -1;
	startLocalRefPos_aln = endLocalRefPos_aln = -1;
	startQueryPos_aln = endQueryPos_aln = -1;

	for(i=0; i<(int32_t)refseq_aln.size(); i++){
		if(refPos==start_var_pos){
			start_aln_idx_var = i;
			startRefPos_aln = refPos;
			startLocalRefPos_aln = localRefPos;
			startQueryPos_aln = queryPos;
			break;
		}
		if(refseq_aln[i]!='-') { refPos ++; localRefPos ++; }
		if(ctgseq_aln[i]!='-') queryPos ++;
	}
	if(start_aln_idx_var!=-1){
		for(; i<(int32_t)refseq_aln.size(); i++){
			if(refPos==end_var_pos){
				tmp_len = 0;
				for(size_t j=i; j<refseq_aln.size(); j++){
					if(refseq_aln[j]=='-') tmp_len ++;
					else break;
				}
				end_aln_idx_var = i + tmp_len;
				endRefPos_aln = refPos;
				endLocalRefPos_aln = localRefPos;
				endQueryPos_aln = queryPos + tmp_len;
				break;
			}
			if(refseq_aln[i]!='-') { refPos ++; localRefPos ++; }
			if(ctgseq_aln[i]!='-') queryPos ++;
		}
	}

	// print the alignment around variations
	if(start_aln_idx_var==-1 or end_aln_idx_var==-1){
		//cerr << "line=" << __LINE__ << ", start_aln_idx_var=" << start_aln_idx_var << ", end_aln_idx_var=" << end_aln_idx_var << endl;
		//exit(1);
	}else {
		// compute the extended gap size on both sides of the region
		// left side
		start_aln_idx_var_new = -1;
		startRefPos_aln_new = startRefPos_aln;
		startLocalRefPos_aln_new = startLocalRefPos_aln;
		startQueryPos_aln_new = startQueryPos_aln;
		for(i=start_aln_idx_var-1; i>=0; i--){
			if(refseq_aln[i]=='-'){ // gap in reference
				start_aln_idx_var_new = i;
				startQueryPos_aln_new --;
			}else if(ctgseq_aln[i]=='-'){ // gap in query
				start_aln_idx_var_new = i;
				startRefPos_aln_new --;
				startLocalRefPos_aln_new --;
			}else
				break;
		}
		// right side
		end_aln_idx_var_new = -1;
		endRefPos_aln_new = endRefPos_aln;
		endLocalRefPos_aln_new = endLocalRefPos_aln;
		endQueryPos_aln_new = endQueryPos_aln;
		for(i=end_aln_idx_var+1; i<(int32_t)midseq_aln.size(); i++){
			if(refseq_aln[i]=='-'){ // gap in reference
				end_aln_idx_var_new = i;
				endQueryPos_aln_new ++;
			}else if(ctgseq_aln[i]=='-'){ // gap in query
				end_aln_idx_var_new = i;
				endRefPos_aln_new ++;
				endLocalRefPos_aln_new ++;
			}else
				break;
		}

		// update information
		if(start_aln_idx_var_new!=-1){
			//cout << "'start_aln_idx_var' is changed from " << start_aln_idx_var << " to " << start_aln_idx_var_new << endl;
			start_aln_idx_var = start_aln_idx_var_new;
			startRefPos_aln = startRefPos_aln_new;
			startLocalRefPos_aln = startLocalRefPos_aln_new;
			startQueryPos_aln = startQueryPos_aln_new;
		}
		if(end_aln_idx_var_new!=-1)
		{
			//cout << "'end_aln_idx_var' is changed from " << end_aln_idx_var << " to " << end_aln_idx_var_new << endl;
			end_aln_idx_var = end_aln_idx_var_new;
			endRefPos_aln = endRefPos_aln_new;
			endLocalRefPos_aln = endLocalRefPos_aln_new;
			endQueryPos_aln = endQueryPos_aln_new;
		}

#if ALIGN_DEBUG
		int32_t len = end_aln_idx_var - start_aln_idx_var + 1;
		if(len>0){
			int32_t start_aln_idx_new, end_aln_idx_new;
			start_aln_idx_new = start_aln_idx_var - 100;
			if(start_aln_idx_new<0) start_aln_idx_new = 0;
			end_aln_idx_new = len + 300;
			if(end_aln_idx_new>(int32_t)midseq_aln.size()) end_aln_idx_new = midseq_aln.size();
			cout << startRefPos_aln << "-" << endRefPos_aln << ", " << start_aln_idx_var << "-" << end_aln_idx_var  << endl;
			cout << "\t" << ctgseq_aln.substr(start_aln_idx_new, end_aln_idx_new) << endl;
			cout << "\t" << midseq_aln.substr(start_aln_idx_new, end_aln_idx_new) << endl;
			cout << "\t" << refseq_aln.substr(start_aln_idx_new, end_aln_idx_new) << endl;
		}
#endif
	}

	local_aln->start_aln_idx_var = start_aln_idx_var;
	local_aln->end_aln_idx_var = end_aln_idx_var;

	if(local_aln->cand_reg) local_aln->reg->query_id = local_aln->cand_reg->query_id;
	if(start_aln_idx_var!=-1 and end_aln_idx_var!=-1){
		local_aln->reg->startRefPos = startRefPos_aln;
		local_aln->reg->endRefPos = endRefPos_aln;
		local_aln->reg->startLocalRefPos = startLocalRefPos_aln;
		local_aln->reg->endLocalRefPos = endLocalRefPos_aln;
		local_aln->reg->startQueryPos = startQueryPos_aln;
		local_aln->reg->endQueryPos = endQueryPos_aln;
	}else{
		local_aln->reg->startLocalRefPos = 0;
		local_aln->reg->endLocalRefPos = 0;
		local_aln->reg->startQueryPos = 0;
		local_aln->reg->endQueryPos = 0;
		local_aln->start_aln_idx_var = local_aln->end_aln_idx_var = -1;
	}
}

// confirm the variant locations, i.e. there are no mismatches on both sides of the region
//void varCand::confirmVarLoc(localAln_t *local_aln){
//	int32_t i, beg_check_idx, end_check_idx, mismatchNum1, mismatchNum2;
//	int32_t refPos_tmp, startRefPos1, endRefPos1, startRefPos2, endRefPos2;
//	string midseq_aln, refseq_aln;
//	reg_t *reg_tmp1, *reg_tmp2;
//	bool valid_flag;
//
//	midseq_aln = local_aln->alignResultVec.at(1);
//	refseq_aln = local_aln->alignResultVec.at(2);
//
//	beg_check_idx = local_aln->start_aln_idx_var - EXT_SIZE_CHK_VAR_LOC;
//	if(beg_check_idx<0) beg_check_idx = 0;
//
//	if(local_aln->alignResultVec[1].size()>0){
//		end_check_idx = local_aln->end_aln_idx_var + EXT_SIZE_CHK_VAR_LOC;
//		if(end_check_idx>local_aln->alignResultVec[1].size()-1) end_check_idx = local_aln->alignResultVec[1].size() - 1;
//	}else
//		end_check_idx = -1;
//
//	valid_flag = true;
//	if(beg_check_idx!=-1 and end_check_idx!=-1){
//		// left side
//		mismatchNum1 = 0;
//		refPos_tmp = local_aln->reg->startRefPos - 1;
//		for(i=local_aln->start_aln_idx_var-1; i>=beg_check_idx; i--){
//			if(midseq_aln.at(i)=='|'){
//				refPos_tmp --;
//			}else {
//				mismatchNum1 ++;
//				if(refseq_aln.at(i)!='-') refPos_tmp --;
//			}
//		}
//		startRefPos1 = refPos_tmp;
//		endRefPos1 = local_aln->reg->startRefPos - 1;
//		reg_tmp1 = findVarvecItem(startRefPos1, endRefPos1, varVec);
//
//		// right side
//		mismatchNum2 = 0;
//		refPos_tmp = local_aln->reg->endRefPos + 1;
//		for(i=local_aln->start_aln_idx_var+1; i<=end_check_idx; i++){
//			if(midseq_aln.at(i)=='|'){
//				refPos_tmp ++;
//			}else{
//				mismatchNum2 ++;
//				if(refseq_aln.at(i)!='-') refPos_tmp ++;
//			}
//		}
//		startRefPos2 = local_aln->reg->endRefPos + 1;
//		endRefPos2 = refPos_tmp;
//		reg_tmp2 = findVarvecItem(startRefPos2, endRefPos2, varVec);
//
//		if((mismatchNum1>0 and reg_tmp1==NULL) or (mismatchNum2>0 and reg_tmp2==NULL))
//			valid_flag = false;
//	}else valid_flag = false;
//
//	if(valid_flag==false){
//		local_aln->start_aln_idx_var = -1;
//		local_aln->end_aln_idx_var = -1;
//
//		local_aln->reg->startLocalRefPos = 0;
//		local_aln->reg->endLocalRefPos = 0;
//		local_aln->reg->startQueryPos = 0;
//		local_aln->reg->endQueryPos = 0;
//		local_aln->reg->query_id = -1;
//	}
//}

// confirm short variation
void varCand::confirmShortVar(localAln_t *local_aln){
	int32_t i, chrlen_tmp, check_extend_size, start_check_idx, end_check_idx, startRefPos, endRefPos, mismatchNum_aln, size_match_misreg_num;
	string refseq_aln;
	vector<int32_t> numVec, misNum_vec;
	vector<mismatchReg_t*> misReg_vec;

	if(local_aln->start_aln_idx_var!=-1 and local_aln->end_aln_idx_var!=-1){
		chrlen_tmp = local_aln->chrlen;

		check_extend_size = SHORT_VAR_ALN_CHECK_EXTEND_SIZE;
		start_check_idx = local_aln->start_aln_idx_var - check_extend_size;
		end_check_idx = local_aln->end_aln_idx_var + check_extend_size;
		if(start_check_idx<0)
			start_check_idx = 0;
		if(end_check_idx>=local_aln->overlapLen)
			end_check_idx = local_aln->overlapLen;

		// compute the start and end reference positions
		refseq_aln = local_aln->alignResultVec[2];
		startRefPos = local_aln->reg->startRefPos;
		for(i=local_aln->start_aln_idx_var; i>=start_check_idx; i--)
			if(refseq_aln[i]!='-')
				startRefPos --;
		if(startRefPos<1) startRefPos = 1;
		endRefPos = local_aln->reg->endRefPos;
		for(i=local_aln->end_aln_idx_var; i<=end_check_idx; i++)
			if(refseq_aln[i]!='-')
				endRefPos ++;
		if(endRefPos>chrlen_tmp) endRefPos = chrlen_tmp;

		// extract mismatch regions and remove short items and short polymer items
		misReg_vec = getMismatchRegVec(local_aln);
		removeShortPolymerMismatchRegItems(local_aln, misReg_vec, inBamFile, fai);

		// compute the number of mismatched bases excluding polymers
		misNum_vec = getMismatchNumAln(local_aln->alignResultVec.at(1), start_check_idx, end_check_idx, misReg_vec, MIN_VALID_POLYMER_SIZE);
		mismatchNum_aln = misNum_vec.at(0);
		size_match_misreg_num = misNum_vec.at(1);

		if(mismatchNum_aln>=MIN_MISMATCH_NUM_SHORT_VAR_CALL){
			numVec = computeDisagreeNumAndHighIndelBaseNum(local_aln->reg->chrname, startRefPos, endRefPos, inBamFile, fai);

			// confirm and update the information
			if(numVec[0]>=MIN_DISAGR_NUM_SHORT_VAR_CALL or numVec[1]>0 or size_match_misreg_num>0) // confirmed
				computeVarType(local_aln->reg);  // compute variant type
		}

		// release memory
		releaseMismatchRegVec(misReg_vec);
	}
}

// get the number of mismatched bases
vector<int32_t> varCand::getMismatchNumAln(string &mid_seq, int32_t start_check_idx, int32_t end_check_idx, vector<mismatchReg_t*> &misReg_vec, int32_t min_match_misReg_size){
	int32_t i, j, mismatchNum, size_match_misreg_num;
	mismatchReg_t *mismatch_reg, *mismatch_reg2;
	vector<mismatchReg_t*> misReg_vec_tmp;
	bool exist_flag;
	vector<int32_t> misNum_vec;

	mismatchNum = size_match_misreg_num = 0;
	for(i=start_check_idx; i<=end_check_idx; i++){
		if(mid_seq[i]!='|'){
			mismatch_reg = getMismatchReg(i, misReg_vec);
			if(mismatch_reg){
				mismatchNum ++;
				if(mismatch_reg->reg_size>=min_match_misReg_size){
					exist_flag = false;
					for(j=0; j<(int32_t)misReg_vec_tmp.size(); j++){
						mismatch_reg2 = misReg_vec_tmp.at(j);
						if(mismatch_reg2==mismatch_reg){
							exist_flag = true;
							break;
						}
					}
					if(exist_flag==false) {
						misReg_vec_tmp.push_back(mismatch_reg);
						size_match_misreg_num ++;
					}
				}
			}
		}
	}

	misNum_vec.push_back(mismatchNum);
	misNum_vec.push_back(size_match_misreg_num);

	return misNum_vec;
}

// compute the number of high indel bases
int32_t varCand::computeHighIndelBaseNum(Base *baseArray, int32_t arr_size, float threshold, float polymer_ignore_ratio_thres){
	int32_t num = 0;
	for(int32_t i=0; i<arr_size; i++)
		if(baseArray[i].isHighConIndelBase(threshold, polymer_ignore_ratio_thres))
			num ++;
	return num;
}

// adjust the short variant locations according to alignment
void varCand::adjustVarLocShortVar(localAln_t *local_aln){
	int32_t start_aln_idx_var_new, end_aln_idx_var_new, tmp;

	if(local_aln->start_aln_idx_var!=-1 and local_aln->end_aln_idx_var!=-1){
		// compute the adjusted locations
		start_aln_idx_var_new = getAdjustedStartAlnIdxVar(local_aln);
		end_aln_idx_var_new = getAdjustedEndAlnIdxVar(local_aln);

		if(start_aln_idx_var_new!=-1 and end_aln_idx_var_new!=-1)
			if(start_aln_idx_var_new>end_aln_idx_var_new){
				tmp = start_aln_idx_var_new;
				start_aln_idx_var_new = end_aln_idx_var_new;
				end_aln_idx_var_new = tmp;

//				pthread_mutex_lock(&mutex_print_var_cand);
//				cout << "adjusted variant locations have been exchanged!" << endl;
//				pthread_mutex_unlock(&mutex_print_var_cand);
			}

		// adjust the locations
		if(start_aln_idx_var_new!=-1 and start_aln_idx_var_new!=local_aln->start_aln_idx_var) // left location
			adjustLeftLocShortVar(local_aln, start_aln_idx_var_new);
		if(end_aln_idx_var_new!=-1 and end_aln_idx_var_new!=local_aln->end_aln_idx_var) // right location
			adjustRightLocShortVar(local_aln, end_aln_idx_var_new);
	}
}

// compute the left mismatch location
int32_t varCand::getAdjustedStartAlnIdxVar(localAln_t *local_aln){
	int32_t i, minCheckIdx, maxCheckIdx, minMismatchIdx;
	string midseq_aln;

	minMismatchIdx = -1;
	if(local_aln->start_aln_idx_var!=-1 and local_aln->end_aln_idx_var!=-1){
		midseq_aln = local_aln->alignResultVec[1];
		minCheckIdx = local_aln->start_aln_idx_var - SHORT_VAR_ALN_CHECK_EXTEND_SIZE;
		maxCheckIdx = local_aln->end_aln_idx_var + SHORT_VAR_ALN_CHECK_EXTEND_SIZE * 3;
		if(minCheckIdx<0) minCheckIdx = 0;
		if(maxCheckIdx>=(int32_t)midseq_aln.size()) maxCheckIdx = midseq_aln.size() - 1;

		// slide to match position
		if(minCheckIdx>0) while(midseq_aln.at(minCheckIdx)!='|') { minCheckIdx--; if(minCheckIdx<=0) break; }
		if(maxCheckIdx<(int32_t)midseq_aln.size()-1) while(midseq_aln.at(maxCheckIdx)!='|') { maxCheckIdx++; if(maxCheckIdx>=(int32_t)midseq_aln.size()-1) break; }

		for(i=local_aln->start_aln_idx_var; i>=minCheckIdx; i--)
			if(midseq_aln[i]!='|')
				minMismatchIdx = i;
		if(minMismatchIdx==-1){
			for(i=local_aln->start_aln_idx_var+1; i<=maxCheckIdx; i++)
				if(midseq_aln[i]!='|'){
					minMismatchIdx = i;
					break;
				}
		}
	}
	if(minMismatchIdx>0) minMismatchIdx --;
	if(minMismatchIdx<0) minMismatchIdx = 0;

	return minMismatchIdx;
}

// compute the right mismatch location
int32_t varCand::getAdjustedEndAlnIdxVar(localAln_t *local_aln){
	int32_t i, minCheckIdx, maxCheckIdx, maxMismatchIdx;
	string midseq_aln;

	midseq_aln = local_aln->alignResultVec[1];
	minCheckIdx = local_aln->start_aln_idx_var - SHORT_VAR_ALN_CHECK_EXTEND_SIZE;
	maxCheckIdx = local_aln->end_aln_idx_var + SHORT_VAR_ALN_CHECK_EXTEND_SIZE * 3;
	if(minCheckIdx<0) minCheckIdx = 0;
	if(maxCheckIdx>=(int32_t)midseq_aln.size()) maxCheckIdx = midseq_aln.size() - 1;

	// slide to match position
	if(minCheckIdx>0) while(midseq_aln.at(minCheckIdx)!='|') { minCheckIdx--; if(minCheckIdx<=0) break; }
	if(maxCheckIdx<(int32_t)midseq_aln.size()-1) while(midseq_aln.at(maxCheckIdx)!='|') { maxCheckIdx++; if(maxCheckIdx>=(int32_t)midseq_aln.size()-1) break; }

	maxMismatchIdx = -1;
	for(i=local_aln->end_aln_idx_var; i<=maxCheckIdx; i++)
		if(midseq_aln[i]!='|')
			maxMismatchIdx = i;
	if(maxMismatchIdx==-1){
		for(i=local_aln->end_aln_idx_var-1; i>=minCheckIdx; i--)
			if(midseq_aln[i]!='|'){
				maxMismatchIdx = i;
				break;
			}
	}
	if(maxMismatchIdx!=-1) maxMismatchIdx ++;
	if(maxMismatchIdx>=(int32_t)midseq_aln.size()) maxMismatchIdx = midseq_aln.size() - 1;

	return maxMismatchIdx;
}

// adjust left variant location
void varCand::adjustLeftLocShortVar(localAln_t *local_aln, int32_t newStartAlnIdx){
	int32_t i, refPos, localRefPos, queryPos;
	string ctgseq_aln, refseq_aln;

	ctgseq_aln = local_aln->alignResultVec[0];
	refseq_aln = local_aln->alignResultVec[2];

	queryPos = local_aln->reg->startQueryPos;
	refPos = local_aln->reg->startRefPos;
	localRefPos = local_aln->reg->startLocalRefPos;

	if(newStartAlnIdx<local_aln->start_aln_idx_var){
		// adjust the start locations
		if(ctgseq_aln[local_aln->start_aln_idx_var]=='-') queryPos--;
		else if(refseq_aln[local_aln->start_aln_idx_var]=='-') { refPos--; localRefPos--; }

		for(i=local_aln->start_aln_idx_var; i>=newStartAlnIdx; i--){
			if(ctgseq_aln[i]!='-') queryPos --;
			if(refseq_aln[i]!='-') { refPos --; localRefPos --; }
		}
		queryPos ++; refPos ++; localRefPos ++;  // return to correct location
	}else{
		for(i=local_aln->start_aln_idx_var; i<=newStartAlnIdx; i++){
			if(ctgseq_aln[i]!='-') queryPos ++;
			if(refseq_aln[i]!='-') { refPos ++; localRefPos ++; }
		}
		queryPos --; refPos --; localRefPos --; // return to correct location
	}

	local_aln->reg->startQueryPos = queryPos;
	local_aln->reg->startRefPos = refPos;
	local_aln->reg->startLocalRefPos = localRefPos;
	local_aln->start_aln_idx_var = newStartAlnIdx;
}

// adjust right variant location
void varCand::adjustRightLocShortVar(localAln_t *local_aln, int32_t newEndAlnIdx){
	int32_t i, refPos, localRefPos, queryPos;
	string ctgseq_aln, refseq_aln;

	ctgseq_aln = local_aln->alignResultVec[0];
	refseq_aln = local_aln->alignResultVec[2];

	queryPos = local_aln->reg->endQueryPos;
	refPos = local_aln->reg->endRefPos;
	localRefPos = local_aln->reg->endLocalRefPos;

	if(newEndAlnIdx<local_aln->end_aln_idx_var){
		// adjust the start locations
		if(ctgseq_aln[local_aln->end_aln_idx_var]=='-') queryPos--;
		else if(refseq_aln[local_aln->end_aln_idx_var]=='-') { refPos--; localRefPos--; }

		for(i=local_aln->end_aln_idx_var; i>=newEndAlnIdx; i--){
			if(ctgseq_aln[i]!='-') queryPos --;
			if(refseq_aln[i]!='-') { refPos --; localRefPos --; }
		}
		queryPos ++; refPos ++; localRefPos ++;  // return to correct location
	}else{
		for(i=local_aln->end_aln_idx_var; i<=newEndAlnIdx; i++){
			if(ctgseq_aln[i]!='-') queryPos ++;
			if(refseq_aln[i]!='-') { refPos ++; localRefPos ++; }
		}
		queryPos --; refPos --; localRefPos --; // return to correct location
	}

	// backward to the last mismatch location
	if(local_aln->alignResultVec[1].at(newEndAlnIdx)=='|'){
		queryPos --; refPos --; localRefPos --; // return to the last mismatch location
		newEndAlnIdx --;
	}

	local_aln->reg->endQueryPos = queryPos;
	local_aln->reg->endRefPos = refPos;
	local_aln->reg->endLocalRefPos = localRefPos;
	local_aln->end_aln_idx_var = newEndAlnIdx;
}

// process the inner contained regions
void varCand::processInnerContainedRegs(vector<reg_t*> &foundRegVec, vector<reg_t*> &candRegVec){
	vector< vector<reg_t*> > regVec;

	// remove variants due to duplicated assembled contigs
	removeVariantsDuplicatedContigs(foundRegVec, candRegVec);

	// merge neighboring variants
	mergeNeighboringVariants(foundRegVec, candRegVec);

	// deal with the new called variants in exist larger variant regions
	regVec = dealWithTwoVariantSets(foundRegVec, candRegVec);

	// update varVec vector
	updateVarVec(regVec, foundRegVec, varVec);
}

// remove variants due to duplicated assembled contigs
void varCand::removeVariantsDuplicatedContigs(vector<reg_t*> &foundRegVec, vector<reg_t*> &candRegVec){
	int32_t i = 1;
	while(i<(int32_t)foundRegVec.size()){
		if(foundRegVec[i]==foundRegVec[i-1] and candRegVec[i]->query_id!=candRegVec[i-1]->query_id){
			delete candRegVec[i];
			candRegVec.erase(candRegVec.begin()+i);
			foundRegVec.erase(foundRegVec.begin()+i);
		}else i++;
	}
}

// merge neighboring variants
void varCand::mergeNeighboringVariants(vector<reg_t*> &foundRegVec, vector<reg_t*> &candRegVec){
	int32_t i, query_dist, ref_dist, var_type, startRefPos, endRefPos, tmp_len;
	bool mergeFlag;
	vector<int32_t> numVec;

	if(foundRegVec.size()>=2){
		i = 1;
		while(i<(int32_t)foundRegVec.size()){
			mergeFlag = false;
			if(foundRegVec[i]==foundRegVec[i-1] and candRegVec[i]->query_id==candRegVec[i-1]->query_id and candRegVec[i]->blat_aln_id==candRegVec[i-1]->blat_aln_id){
				startRefPos = candRegVec[i-1]->endRefPos + 1;
				endRefPos = candRegVec[i]->startRefPos - 1;
				tmp_len = endRefPos - startRefPos + 2;
				if(tmp_len>=SV_MIN_DIST){
					numVec = computeDisagreeNumAndHighIndelBaseNum(candRegVec[i-1]->chrname, startRefPos, endRefPos, inBamFile, fai);
					if(numVec[0]>=MIN_DISAGR_NUM_SHORT_VAR_CALL or numVec[1]>1) // found variant features, merge; otherwise, keep unchanged
						mergeFlag = true;
				}else mergeFlag = true;
			}

			if(mergeFlag){ // merge
				query_dist = candRegVec[i]->endQueryPos - candRegVec[i-1]->startQueryPos;
				ref_dist = candRegVec[i]->endRefPos - candRegVec[i-1]->startRefPos;

				if(query_dist>=ref_dist) var_type = VAR_INS;  // insertion
				else var_type = VAR_DEL;  // deletion

				candRegVec[i-1]->var_type = var_type;
				candRegVec[i-1]->endRefPos = candRegVec[i]->endRefPos;
				candRegVec[i-1]->endLocalRefPos = candRegVec[i]->endLocalRefPos;
				candRegVec[i-1]->endQueryPos = candRegVec[i]->endQueryPos;
				candRegVec[i-1]->sv_len = query_dist - ref_dist;
				candRegVec[i-1]->call_success_status = true;
				candRegVec[i-1]->short_sv_flag = false;

				delete candRegVec[i];
				candRegVec.erase(candRegVec.begin()+i);
				foundRegVec.erase(foundRegVec.begin()+i);
			}else
				i ++;
		}
	}
}

// deal with the two variant sets
vector< vector<reg_t*> > varCand::dealWithTwoVariantSets(vector<reg_t*> &foundRegVec, vector<reg_t*> &candRegVec){
	reg_t *reg, *reg_tmp, *reg_new_1, *reg_new_2, *reg_new_3, *reg_new_4;
	int32_t i, j, idx, same_count, opID1, opID2, startShiftLen, endShiftLen, maxShiftLen, remainShiftLen, ref_dist, query_dist;
	int64_t startRefPos1, endRefPos1, startQueryPos1, endQueryPos1, startRefPos2, endRefPos2, startQueryPos2, endQueryPos2, tmp_pos, query_len;
	int64_t new_startRefPos, new_endRefPos, new_startLocalRefPos, new_endLocalRefPos, new_startQueryPos, new_endQueryPos;
	vector<int32_t> numVec1, numVec2;
	vector<reg_t*> updatedRegVec[foundRegVec.size()], regVec_tmp, processed_reg_vec;
	vector< vector<reg_t*> > regVec;
	bool flag, ref_exchange_flag1, ref_exchange_flag2, query_exchange_flag1, query_exchange_flag2, overlap_flag;

	// deal with the two sets
	i = 0;
	while(i<(int32_t)foundRegVec.size()){
		same_count = 1;
		for(j=i+1; j<(int32_t)foundRegVec.size(); j++)
			if(foundRegVec[j]==foundRegVec[i] and candRegVec[j]->query_id==candRegVec[i]->query_id and candRegVec[j]->blat_aln_id==candRegVec[i]->blat_aln_id) same_count ++;
			else break;

		if(same_count>=2){ // split the variant region in foundRegVec
			for(j=0; j<same_count; j++){
				reg_tmp = candRegVec[i+j];
				reg = new reg_t();
				reg->chrname = reg_tmp->chrname;
				reg->startRefPos = reg_tmp->startRefPos;
				reg->endRefPos = reg_tmp->endRefPos;
				reg->startLocalRefPos = reg_tmp->startLocalRefPos;
				reg->endLocalRefPos = reg_tmp->endLocalRefPos;
				reg->startQueryPos = reg_tmp->startQueryPos;
				reg->endQueryPos = reg_tmp->endQueryPos;
//				reg->startRefPos = reg_tmp->startRefPos - 1;
//				reg->endRefPos = reg_tmp->endRefPos + 1;
//				reg->startLocalRefPos = reg_tmp->startLocalRefPos - 1;
//				reg->endLocalRefPos = reg_tmp->endLocalRefPos + 1;
//				reg->startQueryPos = reg_tmp->startQueryPos - 1;
//				reg->endQueryPos = reg_tmp->endQueryPos + 1;
				reg->query_id = reg_tmp->query_id;
				reg->blat_aln_id = reg_tmp->blat_aln_id;
				reg->aln_orient = reg_tmp->aln_orient;
				reg->zero_cov_flag = false;
				reg->aln_seg_end_flag = false;

				updatedRegVec[i].push_back(reg);
			}
			for(j=1; j<same_count; j++){
				delete candRegVec[i+1];
				candRegVec.erase(candRegVec.begin()+i+1);
				foundRegVec.erase(foundRegVec.begin()+i+1);
			}
		}else{
			idx = getVectorIdx(foundRegVec[i], varVec);
			reg = varVec[idx];
			reg_tmp = candRegVec[i];

			flag = ref_exchange_flag1 = ref_exchange_flag2 = query_exchange_flag1 = query_exchange_flag2 = false;
			for(j=0; j<(int32_t)processed_reg_vec.size(); j++) if(processed_reg_vec.at(j)==reg){ flag = true; break; }

			if(flag==false){
				overlap_flag = isOverlappedRegExtSize(reg, reg_tmp, SHORT_VAR_ALN_CHECK_EXTEND_SIZE, SHORT_VAR_ALN_CHECK_EXTEND_SIZE);
				if(overlap_flag and (reg_tmp->endRefPos-reg_tmp->startRefPos+1<=3 or reg_tmp->endQueryPos-reg_tmp->startQueryPos+1<=3)){
					reg->startRefPos = reg_tmp->startRefPos;
					reg->endRefPos = reg_tmp->endRefPos;
					reg->startLocalRefPos = reg_tmp->startLocalRefPos;
					reg->endLocalRefPos = reg_tmp->endLocalRefPos;
					reg->startQueryPos = reg_tmp->startQueryPos;
					reg->endQueryPos = reg_tmp->endQueryPos;
					reg->query_id = reg_tmp->query_id;
					reg->blat_aln_id = reg_tmp->blat_aln_id;
					reg->aln_orient = reg_tmp->aln_orient;
					updatedRegVec[i].push_back(reg);
					processed_reg_vec.push_back(reg);
					continue;
				}

				// compute the variant locations
				computeVarRegLoc(reg, reg_tmp);

				// set the 'maxShiftLen' to the minimum of 'ref_dist' and 'query_dist'
				ref_dist = reg->endRefPos - reg->startRefPos;
				query_dist = reg->endQueryPos - reg->startQueryPos;
				if(ref_dist<query_dist) maxShiftLen = ref_dist;
				else maxShiftLen = query_dist;
				remainShiftLen = maxShiftLen;

				//cout << "maxShiftLen=" << maxShiftLen << endl;

				FastaSeqLoader ctgseqloader(ctgfilename);
				query_len = ctgseqloader.getFastaSeqLen(reg_tmp->query_id);

				if(reg->startLocalRefPos>0 and reg->endLocalRefPos>0 and reg->startQueryPos>0 and reg->endQueryPos>0){
					if(isOverlappedReg(reg, reg_tmp)==false){ // no overlap and add the variant region
						updatedRegVec[i].push_back(reg);
						processed_reg_vec.push_back(reg);
					}else{ // overlap and try to split the region
						// check left marginal regions
						opID1 = -1;
						if(reg->startRefPos==reg_tmp->startRefPos){
							opID1 = 0;
							numVec1.insert(numVec1.begin(), 6, 0);
							startRefPos1 = endRefPos1 = startQueryPos1 = endQueryPos1 = 0; // will not be used
						}else if(reg->startRefPos<reg_tmp->startRefPos){
							startRefPos1 = reg->startRefPos;
							endRefPos1 = reg_tmp->startRefPos - 1;
							startQueryPos1 = reg->startQueryPos;
							endQueryPos1 = reg_tmp->startQueryPos - 1;
						}else{
							startRefPos1 = reg_tmp->startRefPos;
							endRefPos1 = reg->startRefPos - 1;
							startQueryPos1 = reg_tmp->startQueryPos;
							endQueryPos1 = reg->startQueryPos - 1;
						}
						if(startRefPos1>endRefPos1){ // exchange
							//cout << "aaaaaaaaaaa line=" << __LINE__ << ", ref1 exchanged: " << reg->chrname << ":" << startRefPos1 << "-" << endRefPos1 << endl;
							tmp_pos = startRefPos1;
							startRefPos1 = endRefPos1;
							endRefPos1 = tmp_pos;
							ref_exchange_flag1 = true;
						}
						if(startQueryPos1>endQueryPos1){ // exchange
							//cout << "bbbbbbbbbbb line=" << __LINE__ << ", query1 exchanged: " << reg_tmp->query_id << ":" << startQueryPos1 << "-" << endQueryPos1 << ", ref_reg: " << reg->chrname << ":" << startRefPos1 << "-" << endRefPos1 << endl;
							tmp_pos = startQueryPos1;
							startQueryPos1 = endQueryPos1;
							endQueryPos1 = tmp_pos;
							query_exchange_flag1 = true;
						}

						// check subject length and query length
						if(startQueryPos1<1) startQueryPos1 = 1;
						if(endQueryPos1>query_len) endQueryPos1 = query_len;

						//if(opID1!=0 and startRefPos1<endRefPos1 and startQueryPos1<endQueryPos1){
						if(opID1!=0){
							numVec1 = computeDisagreeNumAndHighIndelBaseNumAndMarginDist(reg->chrname, startRefPos1, endRefPos1, reg_tmp->query_id, startQueryPos1, endQueryPos1, reg->aln_orient, inBamFile, fai);
							if(numVec1[0]>=MIN_DISAGR_NUM_SHORT_VAR_CALL or numVec1[1]>0){ // found variant features
								//if(numVec1[3]>=SV_MIN_DIST and numVec1[3]<endRefPos1-startRefPos1+1)
								if(numVec1[3]>=SV_MIN_DIST and numVec1[4]==0)
									opID1 = 2;
								else
									opID1 = 1;
							}else
								opID1 = 1;
						}

						// check right marginal regions
						opID2 = -1;
						if(reg_tmp->endRefPos==reg->endRefPos){
							opID2 = 0;
							numVec2.insert(numVec2.begin(), 6, 0);
							startRefPos2 = endRefPos2 = startQueryPos2 = endQueryPos2 = 0; // will be not used
						}else if(reg_tmp->endRefPos<reg->endRefPos){
							startRefPos2 = reg_tmp->endRefPos + 1;
							endRefPos2 = reg->endRefPos;
							startQueryPos2 = reg_tmp->endQueryPos + 1;
							endQueryPos2 = reg->endQueryPos;
						}else{
							startRefPos2 = reg->endRefPos + 1;
							endRefPos2 = reg_tmp->endRefPos;
							startQueryPos2 = reg->endQueryPos + 1;
							endQueryPos2 = reg_tmp->endQueryPos;
						}

						if(startRefPos2>endRefPos2){ // exchange to tolerate errors
							//cout << "ccccccccccc line=" << __LINE__ << ", ref2 exchanged: " << reg->chrname << ":" << startRefPos2 << "-" << endRefPos2 << endl;
							tmp_pos = startRefPos2;
							startRefPos2 = endRefPos2;
							endRefPos2 = tmp_pos;
							ref_exchange_flag2 = true;
						}
						if(startQueryPos2>endQueryPos2){ // exchange
							//cout << "ddddddddddd line=" << __LINE__ << ", query2 exchanged: " << reg_tmp->query_id << ":" << startQueryPos2 << "-" << endQueryPos2 << ", ref_reg: " << reg->chrname << ":" << startRefPos2 << "-" << endRefPos2 << endl;
							tmp_pos = startQueryPos2;
							startQueryPos2 = endQueryPos2;
							endQueryPos2 = tmp_pos;
							query_exchange_flag2 = true;
						}

//						if(ref_exchange_flag1 or ref_exchange_flag2 or query_exchange_flag1 or query_exchange_flag2){
//							cout << reg_tmp->chrname << ":" << reg_tmp->startRefPos << "-" << reg_tmp->endRefPos << ", " << reg_tmp->query_id << ":" << reg_tmp->startQueryPos << "-" << reg_tmp->endQueryPos << endl;
//						}

						// check subject length and query length
						if(startQueryPos2<1) startQueryPos2 = 1;
						if(endQueryPos2>query_len) endQueryPos2 = query_len;

						//if(opID2!=0 and startRefPos2<endRefPos2 and startQueryPos2<endQueryPos2){
						if(opID2!=0){
							numVec2 = computeDisagreeNumAndHighIndelBaseNumAndMarginDist(reg->chrname, startRefPos2, endRefPos2, reg_tmp->query_id, startQueryPos2, endQueryPos2, reg->aln_orient, inBamFile, fai);
							if(numVec2[0]>=MIN_DISAGR_NUM_SHORT_VAR_CALL or numVec2[1]>0){ // found variant features
								//if(numVec2[2]>=SV_MIN_DIST and numVec2[2]<endRefPos2-startRefPos2+1)
								if(numVec2[2]>=SV_MIN_DIST and numVec2[4]==0)
									opID2 = 2;
								else
									opID2 = 1;
							}else
								opID2 = 1;
						}

						// update their information at the same time
						if(opID1!=2 and opID2!=2){
							new_startRefPos = reg->startRefPos;
							new_startLocalRefPos = reg->startLocalRefPos;
							new_startQueryPos = reg->startQueryPos;
							new_endRefPos = reg->startRefPos;
							new_endLocalRefPos = reg->startLocalRefPos;
							new_endQueryPos = reg->startQueryPos;

							if(opID1==0){
								new_startRefPos = reg_tmp->startRefPos;
								new_startLocalRefPos = reg_tmp->startLocalRefPos;
								new_startQueryPos = reg_tmp->startQueryPos;
								//reg->startRefPos = reg_tmp->startRefPos;
								//reg->startLocalRefPos = reg_tmp->startLocalRefPos;
								//reg->startQueryPos = reg_tmp->startQueryPos;
							}else if(opID1==1){ // update information
								startShiftLen = getEndShiftLenFromNumVec(numVec1, 1);
								if(startShiftLen>remainShiftLen) startShiftLen = remainShiftLen;
								if(reg->startRefPos+startShiftLen>=reg->endRefPos or reg->startQueryPos+startShiftLen>=reg->endQueryPos
										or reg_tmp->startRefPos+startShiftLen>=reg_tmp->endRefPos or reg_tmp->startQueryPos+startShiftLen>=reg_tmp->endQueryPos)
									startShiftLen = 0;
								if(startShiftLen>0){
									remainShiftLen -= startShiftLen;
									if(reg->startRefPos<reg_tmp->startRefPos){
										new_startRefPos = reg->startRefPos + startShiftLen;
										new_startLocalRefPos = reg->startLocalRefPos + startShiftLen;
										new_startQueryPos = reg->startQueryPos + startShiftLen;
										//reg->startRefPos += numVec1[2];
										//reg->startLocalRefPos += numVec1[2];
										//reg->startQueryPos += numVec1[2];
									}else{
										new_startRefPos = reg_tmp->startRefPos + startShiftLen;
										new_startLocalRefPos = reg_tmp->startLocalRefPos + startShiftLen;
										new_startQueryPos = reg_tmp->startQueryPos + startShiftLen;
										//reg->startRefPos = reg_tmp->startRefPos + numVec1[2];
										//reg->startLocalRefPos = reg_tmp->startLocalRefPos + numVec1[2];
										//reg->startQueryPos = reg_tmp->startQueryPos + numVec1[2];
									}
								}else{
									new_startRefPos = reg->startRefPos;
									new_startLocalRefPos = reg->startLocalRefPos;
									new_startQueryPos = reg->startQueryPos;
								}
							}

							if(opID2==0){
								new_endRefPos = reg_tmp->endRefPos;
								new_endLocalRefPos = reg_tmp->endLocalRefPos;
								new_endQueryPos = reg_tmp->endQueryPos;
								//reg->endRefPos = reg_tmp->endRefPos;
								//reg->endLocalRefPos = reg_tmp->endLocalRefPos;
								//reg->endQueryPos = reg_tmp->endQueryPos;
							}else if(opID2==1){
								endShiftLen = getEndShiftLenFromNumVec(numVec2, 2);
								if(endShiftLen>remainShiftLen) endShiftLen = remainShiftLen;
								if(reg->endRefPos-endShiftLen<=reg->startRefPos or reg->endQueryPos-endShiftLen<=reg->startQueryPos
										or reg_tmp->endRefPos-endShiftLen<=reg_tmp->startRefPos or reg_tmp->endQueryPos-endShiftLen<=reg_tmp->startQueryPos)
									endShiftLen = 0;
								if(endShiftLen>0){
									remainShiftLen -= endShiftLen;
									if(reg_tmp->endRefPos<reg->endRefPos){
										new_endRefPos = reg->endRefPos - endShiftLen;
										new_endLocalRefPos = reg->endLocalRefPos - endShiftLen;
										new_endQueryPos = reg->endQueryPos - endShiftLen;
										//reg->endRefPos -= numVec2[3];
										//reg->endLocalRefPos -= numVec2[3];
										//reg->endQueryPos -= numVec2[3];
									}else{
										new_endRefPos = reg_tmp->endRefPos - endShiftLen;
										new_endLocalRefPos = reg_tmp->endLocalRefPos - endShiftLen;
										new_endQueryPos = reg_tmp->endQueryPos - endShiftLen;
										//reg->endRefPos = reg_tmp->endRefPos - numVec2[3];
										//reg->endLocalRefPos = reg_tmp->endLocalRefPos - numVec2[3];
										//reg->endQueryPos = reg_tmp->endQueryPos - numVec2[3];
									}
								}else{
									new_endRefPos = reg->endRefPos;
									new_endLocalRefPos = reg->endLocalRefPos;
									new_endQueryPos = reg->endQueryPos;
								}
							}

							if(new_startRefPos<=new_endRefPos and new_startQueryPos<=new_endQueryPos){
								reg->startRefPos = new_startRefPos;
								reg->startLocalRefPos = new_startLocalRefPos;
								reg->startQueryPos = new_startQueryPos;
								reg->endRefPos = new_endRefPos;
								reg->endLocalRefPos = new_endLocalRefPos;
								reg->endQueryPos = new_endQueryPos;

								if(isRegValid(reg, 0)){
									exchangeRegLoc(reg);
									updatedRegVec[i].push_back(reg);
									processed_reg_vec.push_back(reg);
								}
							}
						}else if((opID1==1 or opID1==0) and opID2==2){ // check and split the region
							reg_new_3 = new reg_t();
							reg_new_3->startRefPos = reg_tmp->startRefPos;
							reg_new_3->startLocalRefPos = reg_tmp->startLocalRefPos;
							reg_new_3->startQueryPos = reg_tmp->startQueryPos;
							if(reg_tmp->endRefPos<reg->endRefPos){
								reg_new_3->endRefPos = reg_tmp->endRefPos;
								reg_new_3->endLocalRefPos = reg_tmp->endLocalRefPos;
								reg_new_3->endQueryPos = reg_tmp->endQueryPos;
							}else{
								reg_new_3->endRefPos = reg->endRefPos;
								reg_new_3->endLocalRefPos = reg->endLocalRefPos;
								reg_new_3->endQueryPos = reg->endQueryPos;
							}
							reg_new_3->chrname = reg_tmp->chrname;
							reg_new_3->query_id = reg_tmp->query_id;
							reg_new_3->blat_aln_id = reg_tmp->blat_aln_id;
							reg_new_3->aln_orient = reg_tmp->aln_orient;
							reg_new_3->zero_cov_flag = false;
							reg_new_3->aln_seg_end_flag = false;

							reg_new_4 = new reg_t();
							if(reg_tmp->endRefPos<reg->endRefPos){
								reg_new_4->startRefPos = reg_tmp->endRefPos + numVec2[2];
								reg_new_4->startLocalRefPos = reg_tmp->endLocalRefPos + numVec2[2];
								reg_new_4->startQueryPos = reg_tmp->endQueryPos + numVec2[2];

								reg_new_4->endRefPos = reg->endRefPos - numVec2[3];
								reg_new_4->endLocalRefPos = reg->endLocalRefPos - numVec2[3];
								reg_new_4->endQueryPos = reg->endQueryPos - numVec2[3];
							}else{
								reg_new_4->startRefPos = reg->endRefPos + numVec2[2];
								reg_new_4->startLocalRefPos = reg->endLocalRefPos + numVec2[2];
								reg_new_4->startQueryPos = reg->endQueryPos + numVec2[2];

								reg_new_4->endRefPos = reg_tmp->endRefPos - numVec2[3];
								reg_new_4->endLocalRefPos = reg_tmp->endLocalRefPos - numVec2[3];
								reg_new_4->endQueryPos = reg_tmp->endQueryPos - numVec2[3];
							}
							reg_new_4->chrname = reg_tmp->chrname;
							reg_new_4->query_id = reg_tmp->query_id;
							reg_new_4->blat_aln_id = reg_tmp->blat_aln_id;
							reg_new_4->aln_orient = reg_tmp->aln_orient;
							reg_new_4->zero_cov_flag = false;
							reg_new_4->aln_seg_end_flag = false;

							// update information according to opID1==1
							startShiftLen = getEndShiftLenFromNumVec(numVec1, 1);
							if(reg_tmp->startRefPos<reg->startRefPos){
								reg_new_3->startRefPos = reg_tmp->startRefPos + numVec1[2];
								reg_new_3->startLocalRefPos = reg_tmp->startLocalRefPos + numVec1[2];
								reg_new_3->startQueryPos = reg_tmp->startQueryPos + numVec1[2];
							}else{
								reg_new_3->startRefPos = reg->startRefPos + numVec1[2];
								reg_new_3->startLocalRefPos = reg->startLocalRefPos + numVec1[2];
								reg_new_3->startQueryPos = reg->startQueryPos + numVec1[2];
							}

							flag = false;
							if(isRegValid(reg_new_3, SV_MIN_DIST)){
								exchangeRegLoc(reg_new_3);
								updatedRegVec[i].push_back(reg_new_3);
								flag = true;
							}else delete reg_new_3;
							if(isRegValid(reg_new_4, SV_MIN_DIST)){
								exchangeRegLoc(reg_new_4);
								updatedRegVec[i].push_back(reg_new_4);
								flag = true;
							}else delete reg_new_4;

							if(flag) processed_reg_vec.push_back(reg);

						}else if(opID1==2 and (opID2==1 or opID2==0)){ // check and split the region
							reg_new_1 = new reg_t();
							if(reg->startRefPos<reg_tmp->startRefPos){
								reg_new_1->startRefPos = reg->startRefPos + numVec1[2];
								reg_new_1->startLocalRefPos = reg->startLocalRefPos + numVec1[2];
								reg_new_1->startQueryPos = reg->startQueryPos + numVec1[2];

								reg_new_1->endRefPos = reg_tmp->startRefPos - 1 - numVec1[3];
								reg_new_1->endLocalRefPos = reg_tmp->startLocalRefPos - 1 - numVec1[3];
								reg_new_1->endQueryPos = reg_tmp->startQueryPos - 1 - numVec1[3];
							}else{
								reg_new_1->startRefPos = reg_tmp->startRefPos + numVec1[2];
								reg_new_1->startLocalRefPos = reg_tmp->startLocalRefPos + numVec1[2];
								reg_new_1->startQueryPos = reg_tmp->startQueryPos + numVec1[2];

								reg_new_1->endRefPos = reg->startRefPos - 1 - numVec1[3];
								reg_new_1->endLocalRefPos = reg->startLocalRefPos - 1 - numVec1[3];
								reg_new_1->endQueryPos = reg->startQueryPos - 1 - numVec1[3];
							}
							reg_new_1->chrname = reg->chrname;
							reg_new_1->query_id = reg_tmp->query_id;
							reg_new_1->blat_aln_id = reg->blat_aln_id;
							reg_new_1->aln_orient = reg->aln_orient;
							reg_new_1->zero_cov_flag = false;
							reg_new_1->aln_seg_end_flag = false;

							reg_new_2 = new reg_t();
							if(reg->startRefPos<reg_tmp->startRefPos){
								reg_new_2->startRefPos = reg_tmp->startRefPos;
								reg_new_2->startLocalRefPos = reg_tmp->startLocalRefPos;
								reg_new_2->startQueryPos = reg_tmp->startQueryPos;
							}else{
								reg_new_2->startRefPos = reg->startRefPos;
								reg_new_2->startLocalRefPos = reg->startLocalRefPos;
								reg_new_2->startQueryPos = reg->startQueryPos;
							}
							reg_new_2->chrname = reg_tmp->chrname;
							reg_new_2->endRefPos = reg_tmp->endRefPos;
							reg_new_2->endLocalRefPos = reg_tmp->endLocalRefPos;
							reg_new_2->endQueryPos = reg_tmp->endQueryPos;
							reg_new_2->query_id = reg_tmp->query_id;
							reg_new_2->blat_aln_id = reg_tmp->blat_aln_id;
							reg_new_2->aln_orient = reg_tmp->aln_orient;
							reg_new_2->zero_cov_flag = false;
							reg_new_2->aln_seg_end_flag = false;

							// update information according to opID2==1
							endShiftLen = getEndShiftLenFromNumVec(numVec2, 2);
							if(reg_tmp->endRefPos<reg->endRefPos){
								reg_new_2->endRefPos = reg->endRefPos - numVec2[3];
								reg_new_2->endLocalRefPos = reg->endLocalRefPos - numVec2[3];
								reg_new_2->endQueryPos = reg->endQueryPos - numVec2[3];
							}else{
								reg_new_2->endRefPos = reg_tmp->endRefPos - numVec2[3];
								reg_new_2->endLocalRefPos = reg_tmp->endLocalRefPos - numVec2[3];
								reg_new_2->endQueryPos = reg_tmp->endQueryPos - numVec2[3];
							}

							flag = false;
							if(isRegValid(reg_new_1, SV_MIN_DIST)){
								exchangeRegLoc(reg_new_1);
								updatedRegVec[i].push_back(reg_new_1);
								flag = true;
							}else delete reg_new_1;
							if(isRegValid(reg_new_2, SV_MIN_DIST)){
								exchangeRegLoc(reg_new_2);
								updatedRegVec[i].push_back(reg_new_2);
								flag = true;
							}else delete reg_new_2;

							if(flag) processed_reg_vec.push_back(reg);

						}else if(opID1==2 and opID2==2){ // check and split the region
							reg_new_1 = new reg_t();
							if(reg->startRefPos<reg_tmp->startRefPos){
								reg_new_1->startRefPos = reg->startRefPos + numVec1[2];
								reg_new_1->startLocalRefPos = reg->startLocalRefPos + numVec1[2];
								reg_new_1->startQueryPos = reg->startQueryPos + numVec1[2];

								reg_new_1->endRefPos = reg_tmp->startRefPos - 1 - numVec1[3];
								reg_new_1->endLocalRefPos = reg_tmp->startLocalRefPos - 1 - numVec1[3];
								reg_new_1->endQueryPos = reg_tmp->startQueryPos - 1 - numVec1[3];
							}else{
								reg_new_1->startRefPos = reg_tmp->startRefPos + numVec1[2];
								reg_new_1->startLocalRefPos = reg_tmp->startLocalRefPos + numVec1[2];
								reg_new_1->startQueryPos = reg_tmp->startQueryPos + numVec1[2];

								reg_new_1->endRefPos = reg->startRefPos - 1 - numVec1[3];
								reg_new_1->endLocalRefPos = reg->startLocalRefPos - 1 - numVec1[3];
								reg_new_1->endQueryPos = reg->startQueryPos - 1 - numVec1[3];
							}
							reg_new_1->chrname = reg->chrname;
							reg_new_1->query_id = reg_tmp->query_id;
							reg_new_1->blat_aln_id = reg->blat_aln_id;
							reg_new_1->aln_orient = reg->aln_orient;
							reg_new_1->zero_cov_flag = false;
							reg_new_1->aln_seg_end_flag = false;

							reg_new_2 = new reg_t();
							if(reg->startRefPos<reg_tmp->startRefPos){
								reg_new_2->startRefPos = reg_tmp->startRefPos;
								reg_new_2->startLocalRefPos = reg_tmp->startLocalRefPos;
								reg_new_2->startQueryPos = reg_tmp->startQueryPos;
							}else{
								reg_new_2->startRefPos = reg->startRefPos;
								reg_new_2->startLocalRefPos = reg->startLocalRefPos;
								reg_new_2->startQueryPos = reg->startQueryPos;
							}
							if(reg->endRefPos<reg_tmp->endRefPos){
								reg_new_2->endRefPos = reg->endRefPos;
								reg_new_2->endLocalRefPos = reg->endLocalRefPos;
								reg_new_2->endQueryPos = reg->endQueryPos;
							}else{
								reg_new_2->endRefPos = reg_tmp->endRefPos;
								reg_new_2->endLocalRefPos = reg_tmp->endLocalRefPos;
								reg_new_2->endQueryPos = reg_tmp->endQueryPos;
							}
							reg_new_2->chrname = reg_tmp->chrname;
							reg_new_2->query_id = reg_tmp->query_id;
							reg_new_2->blat_aln_id = reg_tmp->blat_aln_id;
							reg_new_2->aln_orient = reg_tmp->aln_orient;
							reg_new_2->zero_cov_flag = false;
							reg_new_2->aln_seg_end_flag = false;

							reg_new_4 = new reg_t();
							if(reg_tmp->endRefPos<reg->endRefPos){
								reg_new_4->startRefPos = reg_tmp->endRefPos + numVec2[2];
								reg_new_4->startLocalRefPos = reg_tmp->endLocalRefPos + numVec2[2];
								reg_new_4->startQueryPos = reg_tmp->endQueryPos + numVec2[2];

								reg_new_4->endRefPos = reg->endRefPos - numVec2[3];
								reg_new_4->endLocalRefPos = reg->endLocalRefPos - numVec2[3];
								reg_new_4->endQueryPos = reg->endQueryPos - numVec2[3];
							}else{
								reg_new_4->startRefPos = reg->endRefPos + numVec2[2];
								reg_new_4->startLocalRefPos = reg->endLocalRefPos + numVec2[2];
								reg_new_4->startQueryPos = reg->endQueryPos + numVec2[2];

								reg_new_4->endRefPos = reg_tmp->endRefPos - numVec2[3];
								reg_new_4->endLocalRefPos = reg_tmp->endLocalRefPos - numVec2[3];
								reg_new_4->endQueryPos = reg_tmp->endQueryPos - numVec2[3];
							}
							reg_new_4->chrname = reg_tmp->chrname;
							reg_new_4->query_id = reg_tmp->query_id;
							reg_new_4->blat_aln_id = reg_tmp->blat_aln_id;
							reg_new_4->aln_orient = reg_tmp->aln_orient;
							reg_new_4->zero_cov_flag = false;
							reg_new_4->aln_seg_end_flag = false;

							flag = false;
							if(isRegValid(reg_new_1, SV_MIN_DIST)){
								exchangeRegLoc(reg_new_1);
								updatedRegVec[i].push_back(reg_new_1);
								flag = true;
							}else delete reg_new_1;
							if(isRegValid(reg_new_2, SV_MIN_DIST)){
								exchangeRegLoc(reg_new_2);
								updatedRegVec[i].push_back(reg_new_2);
								flag = true;
							}else delete reg_new_2;
							if(isRegValid(reg_new_4, SV_MIN_DIST)){
								exchangeRegLoc(reg_new_4);
								updatedRegVec[i].push_back(reg_new_4);
								flag = true;
							}else delete reg_new_4;

							if(flag) processed_reg_vec.push_back(reg);
						}
					}
				}else{
//					reg_new_tmp = new reg_t();
//					reg_new_tmp->chrname = reg_tmp->chrname;
//					reg_new_tmp->startRefPos = reg_tmp->startRefPos - 1;
//					reg_new_tmp->endRefPos = reg_tmp->endRefPos + 1;
//					reg_new_tmp->startLocalRefPos = reg_tmp->startLocalRefPos - 1;
//					reg_new_tmp->endLocalRefPos = reg_tmp->endLocalRefPos + 1;
//					reg_new_tmp->startQueryPos = reg_tmp->startQueryPos - 1;
//					reg_new_tmp->endQueryPos = reg_tmp->endQueryPos + 1;
//					reg_new_tmp->query_id = reg_tmp->query_id;
//					reg_new_tmp->blat_aln_id = reg_tmp->blat_aln_id;
//					reg_new_tmp->aln_orient = reg_tmp->aln_orient;
//					reg_new_tmp->zero_cov_flag = false;
//					reg_new_tmp->aln_seg_end_flag = false;
//
//					updatedRegVec[i].push_back(reg_new_tmp);

					//cout << "******** line=" << __LINE__ << ", new variant region: " << reg_tmp->chrname << ":" << reg_tmp->startRefPos << "-" << reg_tmp->endRefPos << endl;
				}
			}
		}
		i++;
	}

	// push back to regVec
	for(i=0; i<(int32_t)foundRegVec.size(); i++)
		regVec.push_back(updatedRegVec[i]);

	// determine variant type
	for(i=0; i<(int32_t)regVec.size(); i++){
		regVec_tmp = regVec[i];
		for(j=0; j<(int32_t)regVec_tmp.size(); j++){
			reg = regVec_tmp[j];
			if(reg->startLocalRefPos>1 and reg->startQueryPos>1){
				reg->startRefPos --;
				reg->startLocalRefPos --;
				reg->startQueryPos --;
			}
			computeVarType(reg);  // compute variant type
			reg->short_sv_flag = false;
		}
	}

	return regVec;
}

// compute variant type
void varCand::computeVarType(reg_t *reg){
	int64_t query_dist, ref_dist, var_type, decrease_size;
	query_dist = reg->endQueryPos - reg->startQueryPos;
	ref_dist = reg->endRefPos - reg->startRefPos;

	if(query_dist>=ref_dist) var_type = VAR_INS;  // insertion
	else var_type = VAR_DEL;  // deletion

	decrease_size = -1;
	if(ref_dist<0 or query_dist<0){
		if(ref_dist<query_dist)
			decrease_size = -ref_dist;
		else
			decrease_size = -query_dist;
		reg->startRefPos -= decrease_size;
		reg->startLocalRefPos -= decrease_size;
		reg->startQueryPos -= decrease_size;
	}

	reg->var_type = var_type;
	if(query_dist>=ref_dist) reg->sv_len = query_dist - ref_dist; //reg->sv_len = query_dist - ref_dist + 1;
	else reg->sv_len = query_dist - ref_dist; //reg->sv_len = query_dist - ref_dist - 1;
	reg->call_success_status = true;
}

// update varVec vector
void varCand::updateVarVec(vector<vector<reg_t*>> &regVec, vector<reg_t*> &foundRegVec, vector<reg_t*> &varVec){
	size_t i, j;
	int32_t reg_idx;
	reg_t *reg;
	vector<reg_t*> reg_vec_item;
	vector<bool> flag_vec;

	if(regVec.size()>0){
		for(i=0; i<varVec.size(); i++) flag_vec.push_back(true);
		for(i=0; i<foundRegVec.size(); i++){
			reg = foundRegVec[i];
			reg_vec_item = regVec[i];
			if(reg_vec_item.size()>1 or (reg_vec_item.size()>0 and reg_vec_item.at(0)!=reg)){
				reg_idx = getVectorIdx(reg, varVec);
				if(reg_idx!=-1){
					for(j=0; j<reg_vec_item.size(); j++)
					{
						varVec.insert(varVec.begin()+reg_idx+j, reg_vec_item[j]);
						flag_vec.insert(flag_vec.begin()+reg_idx+j, true);
					}
					flag_vec.at(reg_idx+reg_vec_item.size()) = false;
				}else{
					cerr << "line=" << __LINE__ << ", error reg_idx=-1!" << endl;
					exit(1);
				}
			}
		}

		for(i=0; i<varVec.size(); ){
			if(flag_vec.at(i)==false){
				delete varVec.at(i);
				varVec.erase(varVec.begin()+i);
				flag_vec.erase(flag_vec.begin()+i);
			}else i++;
		}
		varVec.shrink_to_fit();
	}
}

// compute the number of disagreements, high indel bases and margin dist for a reference region
vector<int32_t> varCand::computeDisagreeNumAndHighIndelBaseNumAndMarginDist(string &chrname, size_t startRefPos, size_t endRefPos, int32_t query_id, size_t startQueryPos, size_t endQueryPos, size_t aln_orient, string &inBamFile, faidx_t *fai){
	vector<int32_t> numVec, distVec;
	numVec = computeDisagreeNumAndHighIndelBaseNum(chrname, startRefPos, endRefPos, inBamFile, fai);

	// compute margins of variations
	distVec = confirmVarMargins(numVec[2], numVec[3], chrname, startRefPos, endRefPos, query_id, startQueryPos, endQueryPos, aln_orient, fai);
	numVec[2] = distVec[0];
	numVec[3] = distVec[1];
	numVec.push_back(distVec[2]);
	numVec.push_back(distVec[3]);

	return numVec;
}

// compute the number of disagreements and high indel bases in a reference region
vector<int32_t> varCand::computeDisagreeNumAndHighIndelBaseNum(string &chrname, size_t startRefPos, size_t endRefPos, string &inBamFile, faidx_t *fai){
	vector<int32_t> numVec, distVec;
	vector<bam1_t*> alnDataVector;
	int32_t disagreeNum, highIndelBaseNum;

	alnDataLoader data_loader(chrname, startRefPos, endRefPos, inBamFile);
	data_loader.loadAlnData(alnDataVector);

	// load coverage
	covLoader cov_loader(chrname, startRefPos, endRefPos, fai);
	Base *baseArray = cov_loader.initBaseArray();
	cov_loader.generateBaseCoverage(baseArray, alnDataVector);

	// compute the number of disagreements
	disagreeNum = computeDisagreeNum(baseArray, endRefPos-startRefPos+1);

	// compute the number of high indel bases
	highIndelBaseNum = computeHighIndelBaseNum(baseArray, endRefPos-startRefPos+1, MIN_HIGH_INDEL_BASE_RATIO, IGNORE_POLYMER_RATIO_THRES);

	// compute margins of variations
	distVec = computeVarMargins(baseArray, endRefPos-startRefPos+1, MIN_HIGH_INDEL_BASE_RATIO, IGNORE_POLYMER_RATIO_THRES);

	// release the memory
	data_loader.freeAlnData(alnDataVector);
	cov_loader.freeBaseArray(baseArray);

	numVec.push_back(disagreeNum);
	numVec.push_back(highIndelBaseNum);
	numVec.push_back(distVec[0]);
	numVec.push_back(distVec[1]);
	return numVec;
}

// compute the variant locations in VarVec vector
void varCand::computeVarRegLoc(reg_t *reg, reg_t *cand_reg){
	localAln_t *local_aln;

	// prepare align information
	local_aln = new localAln_t();
	local_aln->reg = reg;
	local_aln->cand_reg = cand_reg;
	local_aln->blat_aln_id = -1;
	local_aln->aln_seg = local_aln->start_seg_extend = local_aln->end_seg_extend = NULL;
	local_aln->startRefPos = local_aln->endRefPos = -1;
	local_aln->startLocalRefPos = local_aln->startQueryPos = local_aln->endLocalRefPos = local_aln->endQueryPos = -1;
	local_aln->queryLeftShiftLen = local_aln->queryRightShiftLen = local_aln->localRefLeftShiftLen = local_aln->localRefRightShiftLen = -1;
	local_aln->start_aln_idx_var = local_aln->end_aln_idx_var = -1;
	local_aln->chrlen = faidx_seq_len(fai, local_aln->reg->chrname.c_str()); // get the reference length
	computeLocalLocsAln(local_aln);

	// get align sequences
	FastaSeqLoader refseqloader(refseqfilename);
	local_aln->refseq = refseqloader.getFastaSeqByPos(0, local_aln->startLocalRefPos, local_aln->endLocalRefPos, ALN_PLUS_ORIENT);

	FastaSeqLoader ctgseqloader(ctgfilename);
	local_aln->ctgseq = ctgseqloader.getFastaSeqByPos(cand_reg->query_id, local_aln->startQueryPos, local_aln->endQueryPos, reg->aln_orient);

	// pairwise alignment
	upperSeq(local_aln->ctgseq);
	upperSeq(local_aln->refseq);
	computeSeqAlignment(local_aln);

	// compute locations
	computeVarLoc(local_aln);

	// adjust the variant locations by computing the around gap location
	adjustVarLoc(local_aln);

	// confirm the variant locations, i.e. there are no mismatches on both sides of the region
	//confirmVarLoc(local_aln);

	delete local_aln;
}

// adjust variant locations according to alignment
void varCand::adjustVarLoc(localAln_t *local_aln){
	vector<mismatchReg_t*> misReg_vec;	// all the mismatch regions including gap regions

	// get mismatch regions including gaps according to alignments
	misReg_vec = getMismatchRegVec(local_aln);

	// filter out short (e.g. <= 1 bp) mismatched polymer items
	removeShortPolymerMismatchRegItems(local_aln, misReg_vec, inBamFile, fai);

	// adjust variant locations according to mismatch regions
	adjustVarLocByMismatchRegs(local_aln->reg, misReg_vec, local_aln->start_aln_idx_var, local_aln->end_aln_idx_var);

	// release memory
	releaseMismatchRegVec(misReg_vec);
}

// compute margins of variations
vector<int32_t> varCand::computeVarMargins(Base *baseArray, int32_t arr_size, float threshold, float polymer_ignore_ratio_thres){
	vector<int32_t> distVec;
	int32_t i, j, leftDist, rightDist, discorNum;

	// compute left distance
	leftDist = 0;
	discorNum = 0;
	for(i=0; i<arr_size; i++){
		//if(!baseArray[i].isDisagreeBase() and !baseArray[i].isHighIndelBase(threshold))
		if((!baseArray[i].isDisagreeBase() and !baseArray[i].isHighIndelBase(threshold, polymer_ignore_ratio_thres)) or baseArray[i].isMatchToRef())
			leftDist ++;
		else{
			discorNum ++;
			if(discorNum<MARGIN_MISMATCH_NUM_IGNORE_THRES and leftDist>MIN_MARGIN_LEN_MISMATCH_IGNORE) leftDist ++;
			else break;
		}
	}

	// compute right distance
	rightDist = 0;
	discorNum = 0;
	for(j=arr_size-1; j>i; j--){
		if((!baseArray[j].isDisagreeBase() and !baseArray[j].isHighIndelBase(threshold, polymer_ignore_ratio_thres)) or baseArray[j].isMatchToRef())
			rightDist ++;
		else{
			discorNum ++;
			if(discorNum<MARGIN_MISMATCH_NUM_IGNORE_THRES and rightDist>MIN_MARGIN_LEN_MISMATCH_IGNORE) rightDist ++;
			else break;
		}
	}

	distVec.push_back(leftDist);
	distVec.push_back(rightDist);
	return distVec;
}

// confirm variant margins
vector<int32_t> varCand::confirmVarMargins(int32_t left_dist, int32_t right_dist, string &chrname, size_t startRefPos, size_t endRefPos, int32_t query_id, size_t startQueryPos, size_t endQueryPos, size_t aln_orient, faidx_t *fai){
	vector<int32_t> distVec;
	int32_t i, j, leftDist, rightDist, tmp_len, ref_len, ctg_len, misNum, all_match_flag, all_match_end;
	int32_t left_max, right_max, leftDist_new, rightDist_new;
	string refseq, ctgseq, reg_str;
	bool baseMatchFlag;

	// get the region sequence
	reg_str = chrname + ":" + to_string(startRefPos) + "-" + to_string(endRefPos);
	RefSeqLoader refseq_loader(reg_str, fai);
	refseq_loader.getRefSeq();
	refseq = refseq_loader.refseq;

	FastaSeqLoader ctgseqloader(ctgfilename);
	ctgseq = ctgseqloader.getFastaSeqByPos(query_id, startQueryPos, endQueryPos, aln_orient);

	// pairwise alignment
	upperSeq(ctgseq);
	upperSeq(refseq);

	ref_len = refseq.size();
	ctg_len = ctgseq.size();
	if(ctg_len<ref_len)
		tmp_len = ctg_len;
	else
		tmp_len = ref_len;
	leftDist = 0;
	misNum = 0;
	for(i=0; i<tmp_len; i++){
		baseMatchFlag = isBaseMatch(ctgseq[i], refseq[i]);
		if(baseMatchFlag)
			leftDist ++;
		else{
			misNum ++;
			if(misNum<MARGIN_MISMATCH_NUM_IGNORE_THRES and leftDist>MIN_MARGIN_LEN_MISMATCH_IGNORE) leftDist ++;
			else break;
		}
	}
	rightDist = 0;
	misNum = 0;
	for(j=0; j<tmp_len-i; j++){
		baseMatchFlag = isBaseMatch(ctgseq[ctg_len-1-j], refseq[ref_len-1-j]);
		if(baseMatchFlag)
			rightDist ++;
		else{
			misNum ++;
			if(misNum<MARGIN_MISMATCH_NUM_IGNORE_THRES and rightDist>MIN_MARGIN_LEN_MISMATCH_IGNORE) rightDist ++;
			else break;
		}
	}

	// set the all_base_match flag
	if(leftDist==tmp_len or rightDist==tmp_len){
		all_match_flag = 1;
		if(leftDist==tmp_len)
			all_match_end = 1;
		else if(rightDist==tmp_len)
			all_match_end = 2;
		else{
			cerr << "line=" << __LINE__ << ", tmp_len=" << tmp_len << ", leftDist=" << leftDist << ", rightDist=" << rightDist << ", error!" << endl;
			exit(1);
		}
	}else{
		all_match_flag = 0;
		all_match_end = 0;
	}

	if(leftDist<left_dist) left_max = left_dist;
	else left_max = leftDist;
	if(rightDist<right_dist) right_max = right_dist;
	else right_max = rightDist;

	if(tmp_len<left_max+right_max){
		leftDist_new = leftDist;
		rightDist_new = rightDist;
		if(leftDist>left_dist)
			rightDist_new = tmp_len - leftDist;
		if(rightDist>right_dist)
			leftDist_new = tmp_len - rightDist;

		leftDist = leftDist_new;
		rightDist = rightDist_new;
	}

////	if(leftDist>left_dist)
////		leftDist = left_dist;
////	if(rightDist>right_dist)
////		rightDist = right_dist;
//	if(leftDist<left_dist)
//		leftDist = left_dist;
//	if(rightDist<right_dist)
//		rightDist = right_dist;

	leftDist--; rightDist--;
	if(leftDist<0) leftDist = 0;
	if(rightDist<0) rightDist = 0;

	distVec.push_back(leftDist);
	distVec.push_back(rightDist);
	distVec.push_back(all_match_flag);
	distVec.push_back(all_match_end);

	return distVec;
}

// get the end shift length
int32_t varCand::getEndShiftLenFromNumVec(vector<int32_t> &numVec, size_t end_flag){
	int32_t endShiftLen = 0;

	if(numVec[4]==1){
		if(numVec[5]==1 and end_flag==1)
			endShiftLen = numVec[2];
		else if(numVec[5]==2 and end_flag==2)
			endShiftLen = numVec[3];
	}else{
		if(end_flag==1)
			endShiftLen = numVec[2];
		else if(end_flag==2){
			endShiftLen = numVec[3];
		}else{
			cerr << "line=" << __LINE__ << ", end_flag=" << end_flag << ", error!" << endl;
			exit(1);
		}
	}
	if(endShiftLen>0) endShiftLen --;

	return endShiftLen;
}

// compute local locations
void varCand::computeLocalLocsAln(localAln_t *local_aln){
	reg_t *reg = local_aln->reg, *cand_reg = local_aln->cand_reg;
	int64_t dist0, dist1, dist2, dif_cand_reg;
	int64_t querylen, subject_len;
	blat_aln_t *blat_aln;

	// get the maximal difference size
	if(cand_reg->endRefPos-cand_reg->startRefPos>cand_reg->endQueryPos-cand_reg->startQueryPos)
		dif_cand_reg = cand_reg->endRefPos - cand_reg->startRefPos;
	else
		dif_cand_reg = cand_reg->endQueryPos - cand_reg->startQueryPos;

	if(reg->endRefPos-reg->startRefPos>cand_reg->endRefPos-cand_reg->startRefPos)
		dist0 = ((reg->endRefPos - reg->startRefPos) / VAR_ALN_EXTEND_SIZE + 1) * VAR_ALN_EXTEND_SIZE;
	else
		dist0 = ((cand_reg->endRefPos - cand_reg->startRefPos) / VAR_ALN_EXTEND_SIZE + 1) * VAR_ALN_EXTEND_SIZE;

	if(reg->startRefPos<cand_reg->startRefPos)
		dist1 = dist0 + cand_reg->startRefPos - reg->startRefPos;
	else
		dist1 = dist0;

	dist1 += dif_cand_reg;
	if(cand_reg->startLocalRefPos<cand_reg->startQueryPos){
		if(cand_reg->startLocalRefPos<=dist1) dist1 = cand_reg->startLocalRefPos - 1;
	}else{
		if(cand_reg->startQueryPos<=dist1) dist1 = cand_reg->startQueryPos - 1;
	}

	if(cand_reg->endRefPos<reg->endRefPos)
		dist2 = dist0;
	else
		dist2 = dist0 + cand_reg->endRefPos - reg->endRefPos;

	dist2 += dif_cand_reg;
	blat_aln = blat_aln_vec.at(cand_reg->blat_aln_id);
	querylen = blat_aln->query_len;
	if(cand_reg->endQueryPos+dist2>=querylen) dist2 = querylen - cand_reg->endQueryPos;
	subject_len = blat_aln->subject_len;
	if(cand_reg->endLocalRefPos+dist2>=subject_len) dist2 = subject_len - cand_reg->endLocalRefPos;

	local_aln->startRefPos = cand_reg->startRefPos - dist1;
	local_aln->startLocalRefPos = cand_reg->startLocalRefPos - dist1;
	local_aln->startQueryPos = cand_reg->startQueryPos - dist1;

	local_aln->endRefPos = cand_reg->endRefPos + dist2;
	local_aln->endLocalRefPos = cand_reg->endLocalRefPos + dist2;
	local_aln->endQueryPos = cand_reg->endQueryPos + dist2;
}

// distinguish short duplications and inversions from indels
void varCand::distinguishShortDupInvFromIndels(){
	size_t i, j;
	reg_t *reg;
	string refseq, queryseq;
	vector<size_t> left_shift_size_vec, right_shift_size_vec;
	size_t leftRefShiftSize, leftQueryShiftSize, rightRefShiftSize, rightQueryShiftSize;
	size_t leftVarRefPos, leftVarLocalRefPos, leftVarQueryPos, rightVarRefPos, rightVarLocalRefPos, rightVarQueryPos;
	int32_t dup_num_int, ref_dist, query_dist, sv_len_tmp;
	vector<int32_t> blat_aln_idx_inv_vec;
	double dup_num_tmp, sim_ratio_inv;
	localAln_t *local_aln;

	for(i=0; i<varVec.size(); i++){
		reg = varVec.at(i);

		//cout << "startRefPos=" << reg->startRefPos << ", endRefPos=" << reg->endRefPos << ", startLocalRefPos=" << reg->startLocalRefPos << ", endLocalRefPos=" << reg->endLocalRefPos << ", startQueryPos=" << reg->startQueryPos << ", endQueryPos=" << reg->endQueryPos << ", aln_orient=" << reg->aln_orient << ", var_type=" << reg->var_type << ", sv_len=" << reg->sv_len << ", short_sv_flag=" << reg->short_sv_flag << endl;

		if(reg->call_success_status and reg->short_sv_flag==false){
			if(reg->var_type==VAR_INS and reg->sv_len>=MIN_SHORT_DUP_SIZE){ // distinguish DUP from indels
				// get sequences
				FastaSeqLoader refseqloader(refseqfilename);
				refseq = refseqloader.getFastaSeq(0, ALN_PLUS_ORIENT);
				FastaSeqLoader ctgseqloader(ctgfilename);
				queryseq = ctgseqloader.getFastaSeq(reg->query_id, reg->aln_orient);

				//cout << refseq << endl << endl;
				//cout << queryseq << endl << endl;

				// compute left shift size and right shift size
				left_shift_size_vec = computeLeftShiftSizeDup(reg, NULL, NULL, refseq, queryseq);
				right_shift_size_vec = computeRightShiftSizeDup(reg, NULL, NULL, refseq, queryseq);

				leftRefShiftSize = left_shift_size_vec.at(0);
				leftQueryShiftSize = left_shift_size_vec.at(1);
				rightRefShiftSize = right_shift_size_vec.at(0);
				rightQueryShiftSize = right_shift_size_vec.at(1);

				// compute dup_num
				leftVarRefPos = reg->startRefPos - leftRefShiftSize;
				leftVarLocalRefPos = reg->startLocalRefPos - leftRefShiftSize;
				leftVarQueryPos = reg->startQueryPos - leftQueryShiftSize;
				rightVarRefPos = reg->endRefPos + rightRefShiftSize;
				rightVarLocalRefPos = reg->endLocalRefPos + rightRefShiftSize;
				rightVarQueryPos = reg->endQueryPos + rightQueryShiftSize;

				ref_dist = rightVarRefPos - leftVarRefPos + 1;
				query_dist = rightVarQueryPos - leftVarQueryPos + 1;
				dup_num_tmp = (double)query_dist / ref_dist;
				dup_num_int = round(dup_num_tmp) - 1;
				sv_len_tmp = query_dist - ref_dist;

				//cout << "leftVarRefPos=" << leftVarRefPos << ", rightVarRefPos=" << rightVarRefPos << ", leftVarQueryPos=" << leftVarQueryPos << ", rightVarQueryPos=" << rightVarQueryPos << endl;
				//cout << "leftVarLocalRefPos=" << leftVarLocalRefPos << ", rightVarLocalRefPos=" << rightVarLocalRefPos << endl;
				//cout << "dup_num_int=" << dup_num_int << ", dup_num_tmp=" << dup_num_tmp << ", ref_dist=" << ref_dist << ", query_dist=" << query_dist << ", sv_len_tmp=" << sv_len_tmp << endl;

				if(ref_dist>=MIN_SHORT_DUP_SIZE and sv_len_tmp>=ref_dist){
					reg->var_type = VAR_DUP;
					if(leftVarRefPos>1 and leftVarLocalRefPos>1 and leftVarQueryPos>1){
						reg->startRefPos = leftVarRefPos - 1;
						reg->startLocalRefPos = leftVarLocalRefPos - 1;
						reg->startQueryPos = leftVarQueryPos - 1;
					}else{
						reg->startRefPos = leftVarRefPos;
						reg->startLocalRefPos = leftVarLocalRefPos;
						reg->startQueryPos = leftVarQueryPos;
					}
					reg->endRefPos = rightVarRefPos;
					reg->endLocalRefPos = rightVarLocalRefPos;
					reg->endQueryPos = rightVarQueryPos;
					reg->call_success_status = true;
					reg->short_sv_flag = true;
					reg->dup_num = dup_num_int;
					reg->sv_len = query_dist - ref_dist;

//					pthread_mutex_lock(&mutex_print_var_cand);
//					cout << "----------> Short DUP: " << reg->chrname << ":" << reg->startRefPos << "-" << reg->endRefPos << ", startLocalRefPos=" << reg->startLocalRefPos << ", endLocalRefPos=" << reg->endLocalRefPos << ", startQueryPos=" << reg->startQueryPos << ", endQueryPos=" << reg->endQueryPos << ", sv_len=" << reg->sv_len << endl;
//					pthread_mutex_unlock(&mutex_print_var_cand);
				}

			}else if(reg->sv_len<0.2*(reg->endRefPos-reg->startRefPos+1)){ // distinguish INV from indels
				ref_dist = reg->endRefPos - reg->startRefPos + 1;
				query_dist = reg->endQueryPos - reg->startQueryPos + 1;
				if(ref_dist>=MIN_SHORT_INV_SIZE and query_dist>=MIN_SHORT_INV_SIZE){ // size limitation
					// compute the corresponding blat_aln item
					blat_aln_idx_inv_vec = getInvBlatAlnItemIdx(reg, blat_aln_vec, reg->blat_aln_id, NULL, NULL);
					if(blat_aln_idx_inv_vec.size()>0){
						local_aln = new localAln_t();
						local_aln->reg = NULL;
						local_aln->blat_aln_id = -1;
						local_aln->aln_seg = local_aln->start_seg_extend = local_aln->end_seg_extend = NULL;
						//local_aln->startRefPos = local_aln->endRefPos = -1;
						//local_aln->startLocalRefPos = local_aln->startQueryPos = local_aln->endLocalRefPos = local_aln->endQueryPos = -1;
						local_aln->queryLeftShiftLen = local_aln->queryRightShiftLen = local_aln->localRefLeftShiftLen = local_aln->localRefRightShiftLen = -1;
						//local_aln->chrlen = 0;

						local_aln->startRefPos = reg->startRefPos;
						local_aln->endRefPos = reg->endRefPos;
						local_aln->startLocalRefPos = reg->startLocalRefPos;
						local_aln->endLocalRefPos = reg->endLocalRefPos;
						local_aln->startQueryPos = reg->startQueryPos;
						local_aln->endQueryPos = reg->endQueryPos;
						local_aln->chrlen = faidx_seq_len(fai, reg->chrname.c_str()); // get the reference length

						// get local sequences
						FastaSeqLoader refseqloader(refseqfilename);
						local_aln->refseq = refseqloader.getFastaSeqByPos(0, reg->startLocalRefPos, reg->endLocalRefPos, ALN_PLUS_ORIENT);
						FastaSeqLoader ctgseqloader(ctgfilename);
						local_aln->ctgseq = ctgseqloader.getFastaSeqByPos(reg->query_id, reg->startQueryPos, reg->endQueryPos, reg->aln_orient);
						reverseComplement(local_aln->ctgseq);

						// compute alignment
						upperSeq(local_aln->ctgseq);
						upperSeq(local_aln->refseq);
						computeSeqAlignment(local_aln);

						// compute sequence alignment similarity
						sim_ratio_inv = computeSimilaritySeqAln(local_aln);

						//cout << "sim_ratio_inv=" << sim_ratio_inv << endl;

						if(sim_ratio_inv>=SIMILARITY_THRES_INV){
							if(reg->startRefPos>1 and reg->startLocalRefPos>1 and reg->startQueryPos>1){
								reg->startRefPos --;
								reg->startLocalRefPos --;
								reg->startQueryPos --;
							}
							reg->var_type = VAR_INV;
							reg->call_success_status = true;
							reg->short_sv_flag = true;

//							pthread_mutex_lock(&mutex_print_var_cand);
//							cout << "==========> Short INV: " << reg->chrname << ":" << reg->startRefPos << "-" << reg->endRefPos << ", startLocalRefPos=" << reg->startLocalRefPos << ", endLocalRefPos=" << reg->endLocalRefPos << ", startQueryPos=" << reg->startQueryPos << ", endQueryPos=" << reg->endQueryPos << ", sv_len=" << reg->sv_len << ", sim_ratio_inv=" << sim_ratio_inv << endl;
//							pthread_mutex_unlock(&mutex_print_var_cand);
						}else{
							// compute inversion region
							determineClipRegInvType();
							if(clip_reg){
								for(j=0; j<varVec.size(); ){
									reg = varVec.at(j);
									if(isOverlappedReg(clip_reg, reg)){
										delete reg;
										varVec.erase(varVec.begin()+j);
									}else j++;
								}
								varVec.push_back(clip_reg);
								clip_reg = NULL;
							}
						}

						delete local_aln;
					}
				}
			}
		}
	}
}

// compute sequence alignment similarity
double varCand::computeSimilaritySeqAln(localAln_t *local_aln){
	double sim_ratio;
	size_t leftShiftSize, rightShiftSize, misNum, misNum_total, aln_size_total;

	if(local_aln->localRefLeftShiftLen<local_aln->queryLeftShiftLen) leftShiftSize = local_aln->queryLeftShiftLen;
	else leftShiftSize = local_aln->localRefLeftShiftLen;
	if(local_aln->localRefRightShiftLen<local_aln->queryRightShiftLen) rightShiftSize = local_aln->queryRightShiftLen;
	else rightShiftSize = local_aln->localRefRightShiftLen;

	misNum = computeMismatchNumLocalAln(local_aln);

	misNum_total = misNum + leftShiftSize + rightShiftSize;
	aln_size_total = local_aln->alignResultVec.at(1).size() + leftShiftSize + rightShiftSize;

	sim_ratio = 1 - (double)misNum_total / aln_size_total;

	return sim_ratio;
}

// call variants at align segment ends
void varCand::callVariantsAlnSegEnd(){
	size_t i, j, k;
	string chrname_tmp;
	int64_t startCheckPos_var, endCheckPos_var, chrlen_tmp;
	reg_t *reg, *reg_tmp;
	vector<string> query_name_vec;
	vector<int32_t> query_id_vec_tail_overlap, query_id_vec_head_overlap;
	vector<double> num_vec;
	blat_aln_t *blat_aln;
	aln_seg_t *aln_seg1, *aln_seg2;
	bool fully_contained_flag, exist_flag;
	int32_t disagreeNum, highIndelBaseNum, highIndelClipBaseNum;
	double high_indel_clip_ratio;

	if(assem_success==false or align_success==false) return;	// skip failed items

	// load query names
	FastaSeqLoader fa_loader(ctgfilename);
	query_name_vec = fa_loader.getFastaSeqNames();

	for(i=0; i<varVec.size(); i++){
		reg = varVec.at(i);
		if(reg->call_success_status==false and reg->short_sv_flag==false){
			// the region does not overlap other successfully called regions
			reg_tmp = getOverlappedRegByCallFlag(reg, varVec);
			if(reg_tmp) continue;

			chrlen_tmp = faidx_seq_len(fai, reg->chrname.c_str()); // get the reference length
			//startCheckPos_var = reg->startRefPos - EXT_SIZE_CHK_VAR_LOC * 2;
			//endCheckPos_var = reg->endRefPos + EXT_SIZE_CHK_VAR_LOC * 2;
			startCheckPos_var = reg->startRefPos - EXT_SIZE_CHK_VAR_LOC;
			endCheckPos_var = reg->endRefPos + EXT_SIZE_CHK_VAR_LOC;
			if(startCheckPos_var<1) startCheckPos_var = 1;
			if(endCheckPos_var>chrlen_tmp) endCheckPos_var = chrlen_tmp;

			fully_contained_flag = false;
			for(j=0; j<blat_aln_vec.size(); j++){
				blat_aln = blat_aln_vec.at(j);
				chrname_tmp = fa_loader.getFastaSeqNameByID(blat_aln->query_id);
				if(chrname_tmp.compare(reg->chrname)){
					// check whether no contig contains the variant check region as its sub-region
					aln_seg1 = blat_aln->aln_segs.at(0);
					aln_seg2 = blat_aln->aln_segs.at(blat_aln->aln_segs.size()-1);
					fully_contained_flag = isFullyContainedReg(reg->chrname, startCheckPos_var, endCheckPos_var, chrname_tmp, aln_seg1->ref_start, aln_seg1->ref_end);
					if(fully_contained_flag) break;
				}
			}

			if(fully_contained_flag==false){
				// get contigs whose right ends overlap the variant check region
				for(j=0; j<blat_aln_vec.size(); j++){
					blat_aln = blat_aln_vec.at(j);
					chrname_tmp = fa_loader.getFastaSeqNameByID(blat_aln->query_id);
					if(chrname_tmp.compare(reg->chrname)){
						aln_seg2 = blat_aln->aln_segs.at(blat_aln->aln_segs.size()-1);
						if(aln_seg2->ref_end>=startCheckPos_var and aln_seg2->ref_end<=endCheckPos_var){
							exist_flag = false;
							for(k=0; k<query_id_vec_tail_overlap.size(); k++){
								if(query_id_vec_tail_overlap.at(k)==blat_aln->query_id){
									exist_flag = true;
									break;
								}
							}
							if(exist_flag==false)
								query_id_vec_tail_overlap.push_back(blat_aln->query_id);
						}
					}
				}

				// get contigs whose left ends overlap the variant check region
				for(j=0; j<blat_aln_vec.size(); j++){
					blat_aln = blat_aln_vec.at(j);
					chrname_tmp = fa_loader.getFastaSeqNameByID(blat_aln->query_id);
					if(chrname_tmp.compare(reg->chrname)){
						aln_seg1 = blat_aln->aln_segs.at(0);
						if(aln_seg1->ref_start>=startCheckPos_var and aln_seg1->ref_start<=endCheckPos_var){
							exist_flag = false;
							for(k=0; k<query_id_vec_head_overlap.size(); k++){
								if(query_id_vec_head_overlap.at(k)==blat_aln->query_id){
									exist_flag = true;
									break;
								}
							}
							if(exist_flag==false)
								query_id_vec_head_overlap.push_back(blat_aln->query_id);
						}
					}
				}

				// check variant signatures
				if(!query_id_vec_tail_overlap.empty() and !query_id_vec_head_overlap.empty()){
					// get variant signature
					num_vec = computeDisagreeNumAndHighIndelBaseNumAndClipNum(reg->chrname, startCheckPos_var, endCheckPos_var, inBamFile, fai);
					disagreeNum = num_vec.at(0);
					highIndelBaseNum = num_vec.at(1);
					highIndelClipBaseNum = num_vec.at(2);
					high_indel_clip_ratio = num_vec.at(3);

					if(disagreeNum>0 or highIndelBaseNum>=1 or highIndelClipBaseNum>0 or high_indel_clip_ratio>=HIGH_INDEL_CLIP_BASE_RATIO_THRES){
//						if(reg->startLocalRefPos>1 and reg->startQueryPos>1){
//							reg->startRefPos --;
//							reg->startLocalRefPos --;
//							reg->startQueryPos --;
//						}

						reg->call_success_status = true;
						reg->aln_seg_end_flag = true;
					}
				}
			}
		}
	}
}

// compute the number of disagreements and high indel bases in a reference region
vector<double> varCand::computeDisagreeNumAndHighIndelBaseNumAndClipNum(string &chrname, size_t startRefPos, size_t endRefPos, string &inBamFile, faidx_t *fai){
	vector<int32_t> distVec;
	vector<bam1_t*> alnDataVector;
	int32_t disagreeNum, highIndelBaseNum;
	vector<double> numVec, clipNum_vec;

	alnDataLoader data_loader(chrname, startRefPos, endRefPos, inBamFile);
	data_loader.loadAlnData(alnDataVector);

	// load coverage
	covLoader cov_loader(chrname, startRefPos, endRefPos, fai);
	Base *baseArray = cov_loader.initBaseArray();
	cov_loader.generateBaseCoverage(baseArray, alnDataVector);

	// compute the number of disagreements
	disagreeNum = computeDisagreeNum(baseArray, endRefPos-startRefPos+1);

	// compute the number of high indel bases
	highIndelBaseNum = computeHighIndelBaseNum(baseArray, endRefPos-startRefPos+1, MIN_HIGH_INDEL_BASE_RATIO, IGNORE_POLYMER_RATIO_THRES);

	// compute the number of high ratio indel bases
	clipNum_vec = getTotalHighIndelClipRatioBaseNum(baseArray, endRefPos-startRefPos+1);

	// compute margins of variations
	//distVec = computeVarMargins(baseArray, endRefPos-startRefPos+1, MIN_HIGH_INDEL_BASE_RATIO);

	// release the memory
	data_loader.freeAlnData(alnDataVector);
	cov_loader.freeBaseArray(baseArray);

	numVec.push_back(disagreeNum);
	numVec.push_back(highIndelBaseNum);
	numVec.push_back(clipNum_vec.at(0));
	numVec.push_back(clipNum_vec.at(1));
//	numVec.push_back(distVec[0]);
//	numVec.push_back(distVec[1]);
	return numVec;
}

// adjust variant locations slightly at left end
//void varCand::adjustVarLocSlightly(){
//	size_t i;
//	int64_t j, query_len, start_idx, start_ref_pos, end_ref_pos, start_query_pos, end_query_pos, left_ext_size = EXT_SIZE_CHK_VAR_LOC, ext_size_ref, ext_size_query, ext_size;
//	reg_t *reg;
//	string refseq, altseq;
//
//	if(align_success){
//		if(clip_reg_flag==false){ // indels
//
//			// get the reference sequence
//			FastaSeqLoader refseqloader(refseqfilename);
//			FastaSeqLoader ctgseqloader(ctgfilename);
//			for(i=0; i<varVec.size(); i++){
//				reg = varVec[i];
//				if(reg->call_success_status and reg->aln_seg_end_flag==false){
//					query_len = ctgseqloader.getFastaSeqLen(reg->query_id);
//					start_ref_pos = reg->startLocalRefPos - left_ext_size;
//					end_ref_pos = reg->endLocalRefPos;
//					if(start_ref_pos<1) start_ref_pos = 1;
//					ext_size_ref = reg->startLocalRefPos - start_ref_pos;
//
//					if(reg->aln_orient==ALN_PLUS_ORIENT){
//						start_query_pos = reg->startQueryPos - left_ext_size;
//						end_query_pos = reg->endQueryPos;
//						if(start_query_pos<1) start_query_pos = 1;
//						ext_size_query = reg->startQueryPos - start_query_pos;
//					}else{
//						start_query_pos = reg->startQueryPos;
//						end_query_pos = reg->endQueryPos + left_ext_size;
//						if(end_query_pos>query_len) end_query_pos = query_len;
//						ext_size_query = end_query_pos - reg->endQueryPos;
//					}
//
//					if(ext_size_ref<ext_size_query) ext_size = ext_size_ref;
//					else ext_size = ext_size_query;
//					start_ref_pos = reg->startLocalRefPos - ext_size;
//					if(reg->aln_orient==ALN_PLUS_ORIENT) start_query_pos = reg->startQueryPos - ext_size;
//					else end_query_pos = reg->endQueryPos + ext_size;
//
//					refseq = refseqloader.getFastaSeqByPos(0, start_ref_pos, end_ref_pos, ALN_PLUS_ORIENT);
//					altseq = ctgseqloader.getFastaSeqByPos(reg->query_id, start_query_pos, end_query_pos, reg->aln_orient);
//
//					// compute variant location
//					start_idx = -1;
//					for(j=ext_size-1; j>=0; j--){
//						if(refseq.at(j)==altseq.at(j)){
//							start_idx = j;
//							break;
//						}
//					}
//					if(start_idx!=-1){
//						reg->startRefPos -= ext_size - start_idx;
//						reg->startLocalRefPos -= ext_size - start_idx;
//						if(reg->aln_orient==ALN_PLUS_ORIENT) reg->startQueryPos -= ext_size - start_idx;
//						else reg->endQueryPos += ext_size - start_idx;
//
//						reg->refseq = refseq.substr(start_idx);
//						reg->altseq = altseq.substr(start_idx);
//					}
//				}
//			}
////			for(i=0; i<newVarVec.size(); i++){
////				reg = newVarVec[i];
////				if(reg->call_success_status and reg->aln_seg_end_flag==false)
////					reg->refseq = refseqloader.getFastaSeqByPos(0, reg->startLocalRefPos, reg->endLocalRefPos, ALN_PLUS_ORIENT);
////			}
//
//			// get the query sequence
////			FastaSeqLoader ctgseqloader(ctgfilename);
////			for(i=0; i<varVec.size(); i++){
////				reg = varVec[i];
////				if(reg->call_success_status and reg->aln_seg_end_flag==false)
////					reg->altseq = ctgseqloader.getFastaSeqByPos(reg->query_id, reg->startQueryPos, reg->endQueryPos, reg->aln_orient);
////			}
////			for(i=0; i<newVarVec.size(); i++){
////				reg = newVarVec[i];
////				if(reg->call_success_status and reg->aln_seg_end_flag==false)
////					reg->altseq = ctgseqloader.getFastaSeqByPos(reg->query_id, reg->startQueryPos, reg->endQueryPos, reg->aln_orient);
////			}
//		}else{ // SVs
//			// get the reference and query sequence
//			FastaSeqLoader refseqloader(refseqfilename);
//			FastaSeqLoader ctgseqloader(ctgfilename);
//			if(call_success){
//				if((clip_reg->var_type==VAR_DUP or clip_reg->var_type==VAR_INV) and clip_reg->aln_seg_end_flag==false){
//					query_len = ctgseqloader.getFastaSeqLen(clip_reg->query_id);
//					start_ref_pos = clip_reg->startLocalRefPos - left_ext_size;
//					end_ref_pos = clip_reg->endLocalRefPos;
//					if(start_ref_pos<1) start_ref_pos = 1;
//					ext_size_ref = clip_reg->startLocalRefPos - start_ref_pos;
//
//					if(clip_reg->aln_orient==ALN_PLUS_ORIENT){
//						start_query_pos = clip_reg->startQueryPos - left_ext_size;
//						end_query_pos = clip_reg->endQueryPos;
//						if(start_query_pos<1) start_query_pos = 1;
//						ext_size_query = clip_reg->startQueryPos - start_query_pos;
//					}else{
//						start_query_pos = clip_reg->startQueryPos;
//						end_query_pos = clip_reg->endQueryPos + left_ext_size;
//						if(end_query_pos>query_len) end_query_pos = query_len;
//						ext_size_query = end_query_pos - clip_reg->endQueryPos;
//					}
//
//					if(ext_size_ref<ext_size_query) ext_size = ext_size_ref;
//					else ext_size = ext_size_query;
//					start_ref_pos = clip_reg->startLocalRefPos - ext_size;
//					if(clip_reg->aln_orient==ALN_PLUS_ORIENT) start_query_pos = clip_reg->startQueryPos - ext_size;
//					else end_query_pos = clip_reg->endQueryPos + ext_size;
//
//					refseq = refseqloader.getFastaSeqByPos(0, start_ref_pos, end_ref_pos, ALN_PLUS_ORIENT);
//					altseq = ctgseqloader.getFastaSeqByPos(clip_reg->query_id, start_query_pos, end_query_pos, clip_reg->aln_orient);
//
//					// compute variant location
//					start_idx = -1;
//					for(j=ext_size-1; j>=0; j--){
//						if(refseq.at(j)==altseq.at(j)){
//							start_idx = j;
//							break;
//						}
//					}
//					if(start_idx!=-1){
//						clip_reg->startRefPos -= ext_size - start_idx;
//						clip_reg->startLocalRefPos -= ext_size - start_idx;
//						if(clip_reg->aln_orient==ALN_PLUS_ORIENT) clip_reg->startQueryPos -= ext_size - start_idx;
//						else clip_reg->endQueryPos += ext_size - start_idx;
//
//						clip_reg->refseq = refseq.substr(start_idx);
//						clip_reg->altseq = altseq.substr(start_idx);
//					}
//				}
//			}
//		}
//	}
//
///*
//	// print the seqs
//	if(align_success){
//		for(i=0; i<varVec.size(); i++){
//			reg = varVec[i];
//			if(reg->call_success_status){
//				cout << reg->chrname << ":" << reg->startRefPos << "-" << reg->endRefPos << endl;
//				cout << "\tseqs: " << reg->refseq << "\t" << reg->altseq << endl;
//			}else
//				cout << "cannot be called" << endl;
//		}
//		for(i=0; i<newVarVec.size(); i++){
//			reg = newVarVec[i];
//			if(reg->call_success_status){
//				cout << "New SV: " << reg->chrname << ":" << reg->startRefPos << "-" << reg->endRefPos << endl;
//				cout << "\tseqs: " << reg->refseq << "\t" << reg->altseq << endl;
//			}else
//				cout << "cannot be called" << endl;
//		}
//	}
//*/
//}

// fill the sequence, including reference sequence and contig sequence
void varCand::fillVarseq(){
	size_t i;
	reg_t *reg;

	if(align_success){
		if(clip_reg_flag==false){ // indels
			// get the reference sequence
			FastaSeqLoader refseqloader(refseqfilename);
			for(i=0; i<varVec.size(); i++){
				reg = varVec[i];
				if(reg->call_success_status and reg->aln_seg_end_flag==false)
					reg->refseq = refseqloader.getFastaSeqByPos(0, reg->startLocalRefPos, reg->endLocalRefPos, ALN_PLUS_ORIENT);
			}
			for(i=0; i<newVarVec.size(); i++){
				reg = newVarVec[i];
				if(reg->call_success_status and reg->aln_seg_end_flag==false)
					reg->refseq = refseqloader.getFastaSeqByPos(0, reg->startLocalRefPos, reg->endLocalRefPos, ALN_PLUS_ORIENT);
			}

			// get the query sequence
			FastaSeqLoader ctgseqloader(ctgfilename);
			for(i=0; i<varVec.size(); i++){
				reg = varVec[i];
				if(reg->call_success_status and reg->aln_seg_end_flag==false)
					reg->altseq = ctgseqloader.getFastaSeqByPos(reg->query_id, reg->startQueryPos, reg->endQueryPos, reg->aln_orient);
			}
			for(i=0; i<newVarVec.size(); i++){
				reg = newVarVec[i];
				if(reg->call_success_status and reg->aln_seg_end_flag==false)
					reg->altseq = ctgseqloader.getFastaSeqByPos(reg->query_id, reg->startQueryPos, reg->endQueryPos, reg->aln_orient);
			}
		}else{ // SVs
			// get the reference and query sequence
			FastaSeqLoader refseqloader(refseqfilename);
			FastaSeqLoader ctgseqloader(ctgfilename);
			if(call_success){
				if((clip_reg->var_type==VAR_DUP or clip_reg->var_type==VAR_INV or clip_reg->var_type==VAR_INS or clip_reg->var_type==VAR_DEL) and clip_reg->aln_seg_end_flag==false){
					clip_reg->refseq = refseqloader.getFastaSeqByPos(0, clip_reg->startLocalRefPos, clip_reg->endLocalRefPos, ALN_PLUS_ORIENT);
					clip_reg->altseq = ctgseqloader.getFastaSeqByPos(clip_reg->query_id, clip_reg->startQueryPos, clip_reg->endQueryPos, clip_reg->aln_orient);
				}
			}
		}
	}
}

// determine clipReg variant type
void varCand::determineClipRegVarType(){
	switch(sv_type){
		case VAR_DUP: determineClipRegDupType(); break;
		case VAR_INV: determineClipRegInvType(); break;
		//case VAR_TRA: determineClipRegTraType(); break;
	}
}

// determine clipReg DUP type
void varCand::determineClipRegDupType(){
	size_t i, j;
	blat_aln_t *blat_aln;
	aln_seg_t *seg1, *seg2;
	int32_t segIdx, startClipPos, endClipPos, ref_dist, query_dist, sv_len_tmp;
	string refseq, queryseq;

	startClipPos = endClipPos = -1;
	if(varVec.size()>0){
		startClipPos = varVec.at(0)->startRefPos - 1;
		endClipPos = varVec.at(varVec.size()-1)->endRefPos + 1;
		if(startClipPos<1) startClipPos = 1;
	}else{
		cerr << "varVec.size=" << varVec.size() << ", error!" << endl;
		exit(1);
	}

	for(i=0; i<blat_aln_vec.size(); i++){
		// compute the start and end clipping positions, respectively
		segIdx = -1;
		blat_aln = blat_aln_vec.at(i);
		if(blat_aln->valid_aln==true){
			for(j=1; j<blat_aln->aln_segs.size(); j++){
				seg1 = blat_aln->aln_segs[j-1];
				seg2 = blat_aln->aln_segs[j];

				ref_dist = seg2->ref_start - seg1->ref_end;
				query_dist = seg2->query_start - seg1->query_end;
				sv_len_tmp = query_dist - ref_dist;
				if((ref_dist<MIN_SHORT_DUP_SIZE and sv_len_tmp>=MIN_SHORT_DUP_SIZE) and (seg1->ref_end>=startClipPos and seg1->ref_end<=endClipPos) and (seg2->ref_start>=startClipPos and seg2->ref_start<=endClipPos)){
					segIdx = j - 1;
					break;
				}
			}
			if(segIdx!=-1){
				seg1 = blat_aln->aln_segs.at(segIdx);
				seg2 = blat_aln->aln_segs.at(segIdx+1);

				FastaSeqLoader refseqloader(refseqfilename);
				refseq = refseqloader.getFastaSeq(0, ALN_PLUS_ORIENT);
				FastaSeqLoader ctgseqloader(ctgfilename);
				queryseq = ctgseqloader.getFastaSeq(blat_aln->query_id, blat_aln->aln_orient);

				//cout << refseq << endl << endl;
				//cout << queryseq << endl << endl;

				// compute left and right clipping positions
				clip_reg = computeClipPos(blat_aln, seg1, seg2, refseq, queryseq, sv_type);
			}
		}
	}
}

// determine clipReg INV type
void varCand::determineClipRegInvType(){
	size_t i, j, k;
	blat_aln_t *blat_aln;
	aln_seg_t *seg1, *seg2;
	int32_t segIdx_start, segIdx_end, startClipPos, endClipPos, ref_dist, query_dist;
	vector<int32_t> blat_aln_idx_inv_vec;
	string refseq, queryseq;
	localAln_t *local_aln;
	double sim_ratio_inv;

	startClipPos = endClipPos = -1;
	if(varVec.size()>0){
		startClipPos = varVec.at(0)->startRefPos;
		endClipPos = varVec.at(varVec.size()-1)->endRefPos;
	}else{
		cerr << "varVec.size=" << varVec.size() << ", error!" << endl;
		exit(1);
	}

	for(i=0; i<blat_aln_vec.size(); i++){
		// compute the start and end clipping positions, respectively
		segIdx_start = -1; segIdx_end = -1;
		blat_aln = blat_aln_vec.at(i);
		if(blat_aln->valid_aln==true){
			for(j=1; j<blat_aln->aln_segs.size(); j++){
				seg1 = blat_aln->aln_segs[j-1];
				seg2 = blat_aln->aln_segs[j];
				if((seg1->ref_end>=startClipPos-CLIP_END_EXTEND_SIZE and seg1->ref_end<=endClipPos+CLIP_END_EXTEND_SIZE) and (seg2->ref_start>=startClipPos-CLIP_END_EXTEND_SIZE and seg2->ref_start<=endClipPos+CLIP_END_EXTEND_SIZE)){
					segIdx_start = j - 1;

					// compute segIdx_end
					segIdx_end = -1;
					for(k=segIdx_start+1; k<blat_aln->aln_segs.size(); k++){
						seg2 = blat_aln->aln_segs.at(k);
						if(seg2->ref_end-seg2->ref_start+1>=MIN_VALID_BLAT_SEG_SIZE){ // ignore short align segments
							segIdx_end = k;
							break;
						}
					}
					if(segIdx_end!=-1) break;
				}
			}

			if(segIdx_start!=-1 and segIdx_end!=-1){
				seg1 = blat_aln->aln_segs.at(segIdx_start);
				seg2 = blat_aln->aln_segs.at(segIdx_end);

				// compute the corresponding blat_aln item
				blat_aln_idx_inv_vec = getInvBlatAlnItemIdx(NULL, blat_aln_vec, i, seg1, seg2);
				if(blat_aln_idx_inv_vec.size()>0){
					clip_reg = new reg_t();
					clip_reg->chrname = chrname;
					clip_reg->startRefPos =  seg1->ref_end + 1;
					clip_reg->endRefPos =  seg2->ref_start - 1;
					clip_reg->startLocalRefPos =  seg1->subject_end + 1;
					clip_reg->endLocalRefPos =  seg2->subject_start - 1;
					clip_reg->startQueryPos =  seg1->query_end + 1;
					clip_reg->endQueryPos =  seg2->query_start - 1;
					clip_reg->aln_orient = blat_aln->aln_orient;
					clip_reg->var_type = sv_type;
					clip_reg->query_id = blat_aln->query_id;
					clip_reg->blat_aln_id = i;
					clip_reg->call_success_status = true;
					clip_reg->short_sv_flag = false;
					clip_reg->zero_cov_flag = false;
					clip_reg->aln_seg_end_flag = false;

					ref_dist = clip_reg->endLocalRefPos - clip_reg->startLocalRefPos + 1;
					query_dist = clip_reg->endQueryPos - clip_reg->startQueryPos + 1;
					clip_reg->sv_len = query_dist - ref_dist + 1;

					// check the inverted sequences
					local_aln = new localAln_t();
					local_aln->reg = NULL;
					local_aln->cand_reg = NULL;
					local_aln->blat_aln_id = -1;
					local_aln->aln_seg = local_aln->start_seg_extend = local_aln->end_seg_extend = NULL;
					//local_aln->startRefPos = local_aln->endRefPos = -1;
					//local_aln->startLocalRefPos = local_aln->startQueryPos = local_aln->endLocalRefPos = local_aln->endQueryPos = -1;
					local_aln->queryLeftShiftLen = local_aln->queryRightShiftLen = local_aln->localRefLeftShiftLen = local_aln->localRefRightShiftLen = -1;
					//local_aln->chrlen = 0;

					local_aln->startRefPos = clip_reg->startRefPos;
					local_aln->endRefPos = clip_reg->endRefPos;
					local_aln->startLocalRefPos = clip_reg->startLocalRefPos;
					local_aln->endLocalRefPos = clip_reg->endLocalRefPos;
					local_aln->startQueryPos = clip_reg->startQueryPos;
					local_aln->endQueryPos = clip_reg->endQueryPos;
					local_aln->chrlen = faidx_seq_len(fai, clip_reg->chrname.c_str()); // get the reference length

					// get local sequences
					FastaSeqLoader refseqloader(refseqfilename);
					local_aln->refseq = refseqloader.getFastaSeqByPos(0, clip_reg->startLocalRefPos, clip_reg->endLocalRefPos, ALN_PLUS_ORIENT);
					FastaSeqLoader ctgseqloader(ctgfilename);
					local_aln->ctgseq = ctgseqloader.getFastaSeqByPos(clip_reg->query_id, clip_reg->startQueryPos, clip_reg->endQueryPos, clip_reg->aln_orient);
					reverseComplement(local_aln->ctgseq);

					// compute alignment
					upperSeq(local_aln->ctgseq);
					upperSeq(local_aln->refseq);
					computeSeqAlignment(local_aln);

					// compute sequence alignment similarity
					sim_ratio_inv = computeSimilaritySeqAln(local_aln);
					if(sim_ratio_inv>=SIMILARITY_THRES_INV){
						if(clip_reg->startLocalRefPos>1 and clip_reg->startQueryPos>1){
							clip_reg->startRefPos --;
							clip_reg->startLocalRefPos --;
							clip_reg->startQueryPos --;
						}

						clip_reg->var_type = VAR_INV;
						clip_reg->call_success_status = true;
						clip_reg->short_sv_flag = false;
						call_success = true;

//						pthread_mutex_lock(&mutex_print_var_cand);
//						cout << "----------> INV: " << clip_reg->chrname << ":" << clip_reg->startRefPos << "-" << clip_reg->endRefPos << ", startLocalRefPos=" << clip_reg->startLocalRefPos << ", endLocalRefPos=" << clip_reg->endLocalRefPos << ", startQueryPos=" << clip_reg->startQueryPos << ", endQueryPos=" << clip_reg->endQueryPos << ", sv_len=" << clip_reg->sv_len << ", sim_ratio_inv=" << sim_ratio_inv << endl;
//						pthread_mutex_unlock(&mutex_print_var_cand);
					}else{
						delete clip_reg;
						clip_reg = NULL;
					}

					delete local_aln;

					break;
				}
			}
		}
	}
}

vector<int32_t> varCand::getInvBlatAlnItemIdx(reg_t *reg, vector<blat_aln_t*> &blat_aln_vec, size_t given_idx, aln_seg_t *seg1, aln_seg_t *seg2){
	int32_t blat_aln_idx = -1, aln_seg_idx = -1;
	vector<int32_t> blat_aln_idx_vec;
	size_t i, j, aln_orient1, aln_orient2, startRefPos_Gap, endRefPos_Gap, startRefPos_Inv, endRefPos_Inv;
	blat_aln_t *blat_aln, *blat_aln2;
	vector<aln_seg_t*> aln_segs;
	aln_seg_t *aln_seg;
	bool flag;

	blat_aln = blat_aln_vec.at(given_idx);
	aln_orient1 = blat_aln->aln_orient;

	if(reg==NULL){ // for inverted clipping regions
		startRefPos_Gap = seg1->ref_end + 1;
		endRefPos_Gap = seg2->ref_start - 1;
	}else{ // for short inverted segments
		startRefPos_Gap = reg->startRefPos;
		endRefPos_Gap = reg->endRefPos;
	}

	for(i=0; i<blat_aln_vec.size(); i++){
		if(i!=given_idx){
			blat_aln2 = blat_aln_vec.at(i);
			aln_orient2 = blat_aln2->aln_orient;
			if(blat_aln2->query_id==blat_aln->query_id and aln_orient1!=aln_orient2){
				// get the start and end locations of the inversion sequence
				aln_segs = blat_aln2->aln_segs;
				for(j=0; j<aln_segs.size(); j++){
					aln_seg = aln_segs.at(j);
					startRefPos_Inv = aln_seg->ref_start;
					endRefPos_Inv = aln_seg->ref_end;
					flag = isOverlappedPos(startRefPos_Gap, endRefPos_Gap, startRefPos_Inv, endRefPos_Inv);
					if(flag){ // overlapped
						blat_aln_idx = i;
						aln_seg_idx = j;
						break;
					}
				}
			}
			if(blat_aln_idx!=-1){
				blat_aln_idx_vec.push_back(blat_aln_idx);
				blat_aln_idx_vec.push_back(aln_seg_idx);
				break;
			}
		}
	}

	return blat_aln_idx_vec;
}

// compute left and right clipping positions
reg_t* varCand::computeClipPos(blat_aln_t *blat_aln, aln_seg_t *seg1, aln_seg_t *seg2, string &refseq, string &queryseq, size_t var_type){
	int64_t leftRefShiftSize, leftQueryShiftSize, rightRefShiftSize, rightQueryShiftSize;
	int64_t leftClipRefPos, leftClipLocalRefPos, leftClipQueryPos, rightClipRefPos, rightClipLocalRefPos, rightClipQueryPos;
	int64_t dup_num_int, ref_dist, query_dist, minLeftClipPos, maxLeftClipPos, minRightClipPos, maxRightClipPos;
	double dup_num_tmp;
	vector<size_t> left_shift_size_vec, right_shift_size_vec;
	vector<clipPos_t> clip_pos_vec;
	clipPos_t clip_pos;
	reg_t *clip_reg_ret = NULL, *reg_new;
	bool call_success_flag = true;
	vector<size_t> pos_vec;
	//localAln_t *local_aln;
	//int64_t localRefPos_start, queryPos_start, subseq_len, left_ext_size, right_ext_size, dist_max, chrlen_tmp;

	if(var_type==VAR_DUP){ // duplication
		left_shift_size_vec = computeLeftShiftSizeDup(varVec.at(0), seg1, seg2, refseq, queryseq);
		right_shift_size_vec = computeRightShiftSizeDup(varVec.at(varVec.size()-1), seg1, seg2, refseq, queryseq);

		if((left_shift_size_vec.at(0)>=5 or left_shift_size_vec.at(1)>=5) and (right_shift_size_vec.at(0)>=5 or right_shift_size_vec.at(1)>=5)){
			leftRefShiftSize = left_shift_size_vec.at(0);
			leftQueryShiftSize = left_shift_size_vec.at(1);
			rightRefShiftSize = right_shift_size_vec.at(0);
			rightQueryShiftSize = right_shift_size_vec.at(1);

			leftClipRefPos = seg1->ref_end - leftRefShiftSize + 1;
			leftClipLocalRefPos = seg1->subject_end - leftRefShiftSize + 1;
			leftClipQueryPos = seg1->query_end - leftQueryShiftSize + 1;
			rightClipRefPos = seg2->ref_start + rightRefShiftSize - 1;
			rightClipLocalRefPos = seg2->subject_start + rightRefShiftSize - 1;
			rightClipQueryPos = seg2->query_start + rightQueryShiftSize - 1;

			ref_dist = rightClipRefPos - leftClipRefPos + 1;
			query_dist = rightClipQueryPos - leftClipQueryPos + 1;
			dup_num_tmp = (double)query_dist / ref_dist;
			dup_num_int = round(dup_num_tmp) - 1;

			//cout << "leftClipRefPos=" << leftClipRefPos << ", rightClipRefPos=" << rightClipRefPos << ", leftClipQueryPos=" << leftClipQueryPos << ", rightClipQueryPos=" << rightClipQueryPos << endl;
			//cout << "leftClipLocalRefPos=" << leftClipLocalRefPos << ", rightClipLocalRefPos=" << rightClipLocalRefPos << endl;
			//cout << "dup_num_int=" << dup_num_int << ", dup_num_tmp=" << dup_num_tmp << ", ref_dist=" << ref_dist << ", query_dist=" << query_dist << endl;

	//		if(dup_num_int!=(int32_t)dup_num){
	//			cout << "************* dup_num_int=" << dup_num_int << ", dup_num=" << dup_num << " *************" << endl;
	//		}

			minLeftClipPos = varVec.at(0)->startRefPos - VAR_ALN_EXTEND_SIZE;
			maxLeftClipPos = varVec.at(0)->endRefPos + VAR_ALN_EXTEND_SIZE;
			minRightClipPos = varVec.at(varVec.size()-1)->startRefPos - VAR_ALN_EXTEND_SIZE;
			maxRightClipPos = varVec.at(varVec.size()-1)->endRefPos + VAR_ALN_EXTEND_SIZE;
			if(minLeftClipPos<1) minLeftClipPos = 1;
			if(minRightClipPos<1) minRightClipPos = 1;

			margin_adjusted_flag = false;
			if(leftClipRefPos>minLeftClipPos and (int32_t)leftClipRefPos<maxLeftClipPos){
				clip_pos.chrname = varVec.at(0)->chrname;
				clip_pos.clipRefPos = leftClipRefPos;
				clip_pos.clipLocalRefPos = leftClipLocalRefPos;
				clip_pos.clipQueryPos = leftClipQueryPos;
				clip_pos.aln_orient = seg1->aln_orient;
			}else{ // use the detected clip position instead
				pos_vec = computeQueryClipPosDup(blat_aln, this->leftClipRefPos, refseq, queryseq); // compute query position
				if(pos_vec.size()>0){
					clip_pos.chrname = varVec.at(0)->chrname;
					clip_pos.clipRefPos = pos_vec.at(0);
					clip_pos.clipLocalRefPos = pos_vec.at(1);
					clip_pos.clipQueryPos = pos_vec.at(2);
					clip_pos.aln_orient = blat_aln->aln_orient;
				}else{
					clip_pos.chrname = "";
					clip_pos.clipRefPos = 0;
					clip_pos.clipLocalRefPos = 0;
					clip_pos.clipQueryPos = 0;
					clip_pos.aln_orient = 0;
					call_success_flag = false;
				}
				margin_adjusted_flag = true;
			}
			clip_pos_vec.push_back(clip_pos);

			if(rightClipRefPos>minRightClipPos and rightClipRefPos<maxRightClipPos){
				clip_pos.chrname = varVec.at(varVec.size()-1)->chrname;
				clip_pos.clipRefPos = rightClipRefPos;
				clip_pos.clipLocalRefPos = rightClipLocalRefPos;
				clip_pos.clipQueryPos = rightClipQueryPos;
				clip_pos.aln_orient = blat_aln->aln_orient;
			}else{ // use the detected clip position instead
				pos_vec = computeQueryClipPosDup(blat_aln, this->rightClipRefPos, refseq, queryseq); // compute query position
				if(pos_vec.size()>0){
					clip_pos.chrname = varVec.at(varVec.size()-1)->chrname;
					clip_pos.clipRefPos = pos_vec.at(0);
					clip_pos.clipLocalRefPos = pos_vec.at(1);
					clip_pos.clipQueryPos = pos_vec.at(2);
					clip_pos.aln_orient = blat_aln->aln_orient;
				}else{
					clip_pos.chrname = "";
					clip_pos.clipRefPos = 0;
					clip_pos.clipLocalRefPos = 0;
					clip_pos.clipQueryPos = 0;
					clip_pos.aln_orient = 0;
					call_success_flag = false;
				}
				margin_adjusted_flag = true;
			}
			clip_pos_vec.push_back(clip_pos);

			if(call_success_flag){
				call_success = true;
				clip_reg_ret = new reg_t();
				clip_reg_ret->chrname = clip_pos_vec.at(0).chrname;
				clip_reg_ret->startRefPos =  clip_pos_vec.at(0).clipRefPos - 1;
				clip_reg_ret->endRefPos =  clip_pos_vec.at(1).clipRefPos + 1;
				clip_reg_ret->startLocalRefPos =  clip_pos_vec.at(0).clipLocalRefPos - 1;
				clip_reg_ret->endLocalRefPos =  clip_pos_vec.at(1).clipLocalRefPos + 1;
				clip_reg_ret->startQueryPos =  clip_pos_vec.at(0).clipQueryPos - 1;
				clip_reg_ret->endQueryPos =  clip_pos_vec.at(1).clipQueryPos + 1;
				clip_reg_ret->aln_orient = clip_pos_vec.at(0).aln_orient;
				clip_reg_ret->var_type = var_type;
				clip_reg_ret->query_id = blat_aln->query_id;
				clip_reg_ret->blat_aln_id = blat_aln->blat_aln_id;
				clip_reg_ret->call_success_status = true;
				clip_reg_ret->short_sv_flag = false;
				clip_reg_ret->zero_cov_flag = false;
				clip_reg_ret->aln_seg_end_flag = false;

				ref_dist = clip_reg_ret->endLocalRefPos - clip_reg_ret->startLocalRefPos + 1;
				query_dist = clip_reg_ret->endQueryPos - clip_reg_ret->startQueryPos + 1;
				clip_reg_ret->sv_len = query_dist - ref_dist + 1;
				dup_num_tmp = (double)query_dist / ref_dist;
				dup_num_int = round(dup_num_tmp) - 1;
				//cout << "->->->->->-> dup_num_int=" << dup_num_int << ", dup_num_tmp=" << dup_num_tmp << ", ref_dist=" << ref_dist << ", query_dist=" << query_dist << ", margin adjusted=" << margin_adjusted_flag << endl;

				if(margin_adjusted_flag==false){
					clip_reg_ret->dup_num = dup_num_int;
				}else{
					clip_reg_ret->dup_num = dup_num;
				}

				//cout << "Final Clip Pos: chrname=" << clip_reg_ret->chrname << ", leftClipRefPos=" << clip_reg_ret->startRefPos << ", rightClipRefPos=" << clip_reg_ret->endRefPos << ", dup_num=" << clip_reg_ret->dup_num << ", sv_len=" << clip_reg_ret->sv_len << endl;
			}
		}else{ // it may be insertion
/*
			left_ext_size = 4 * EXT_SIZE_CHK_VAR_LOC;
			if(seg1->ref_end-1<left_ext_size) left_ext_size = seg1->ref_end - 1;
			if(seg1->subject_end-1<left_ext_size) left_ext_size = seg1->subject_end - 1;
			if(seg1->query_end-1<left_ext_size) left_ext_size = seg1->query_end - 1;

			right_ext_size = 4 * EXT_SIZE_CHK_VAR_LOC;
			chrlen_tmp = faidx_seq_len(fai, chrname.c_str()); // get the reference length
			if(chrlen_tmp-seg2->subject_start<right_ext_size) right_ext_size = chrlen_tmp - seg2->subject_start;
			if((int64_t)refseq.size()-seg2->subject_start<right_ext_size) right_ext_size = refseq.size() - seg2->subject_start;
			if((int64_t)queryseq.size()-seg2->query_start<right_ext_size) right_ext_size = queryseq.size() - seg2->query_start;

			dist_max = seg2->query_start - seg1->query_end;
			if(seg2->subject_start-seg1->subject_end>dist_max) dist_max = seg2->subject_start - seg1->subject_end;

			localRefPos_start = seg1->subject_end - left_ext_size;
			queryPos_start = seg1->query_end - left_ext_size;

			subseq_len = left_ext_size + right_ext_size + dist_max;
*/
			// try new indel
			reg_new = new reg_t();
			reg_new->chrname = chrname;
			reg_new->startRefPos = seg1->ref_end;
			reg_new->startLocalRefPos = seg1->subject_end;
			reg_new->startQueryPos = seg1->query_end;
			reg_new->endRefPos = seg2->ref_start;
			reg_new->endLocalRefPos = seg2->subject_start;
			reg_new->endQueryPos = seg2->query_start;
			reg_new->blat_aln_id = blat_aln->blat_aln_id;
			reg_new->query_id = blat_aln->query_id;
			reg_new->aln_orient = blat_aln->aln_orient;
			reg_new->short_sv_flag = false;
			reg_new->zero_cov_flag = false;
			reg_new->aln_seg_end_flag = false;
/*
			// compute local locations
			local_aln = new localAln_t();
			local_aln->reg = reg_new;
			local_aln->blat_aln_id = blat_aln->blat_aln_id;
			local_aln->aln_seg = local_aln->start_seg_extend = local_aln->end_seg_extend = NULL;
			local_aln->startRefPos = seg1->ref_end - left_ext_size;
			local_aln->startLocalRefPos = seg1->subject_end - left_ext_size;
			local_aln->startQueryPos = seg1->query_end - left_ext_size;
			local_aln->endRefPos = seg2->ref_start + right_ext_size;
			local_aln->endLocalRefPos = seg2->subject_start + right_ext_size;
			local_aln->endQueryPos = seg2->query_start + right_ext_size;
			local_aln->queryLeftShiftLen = local_aln->queryRightShiftLen = local_aln->localRefLeftShiftLen = local_aln->localRefRightShiftLen = -1;
			local_aln->chrlen = chrlen_tmp; // reference length

			//computeLocalLocsAln(local_aln);

			// get sub local sequence
			local_aln->refseq = refseq.substr(localRefPos_start-1, subseq_len);
			local_aln->ctgseq = queryseq.substr(queryPos_start-1, subseq_len);

			// pairwise alignment
			upperSeq(local_aln->ctgseq);
			upperSeq(local_aln->refseq);
			computeSeqAlignment(local_aln);

			// compute locations
			computeVarLoc(local_aln);

			// adjust the variant locations by computing the around gap location
			adjustVarLoc(local_aln);
*/
//			if(local_aln->start_aln_idx_var!=-1 and local_aln->end_aln_idx_var!=-1){
//				if(reg_new->startLocalRefPos>1 and reg_new->startQueryPos>1){
//					reg_new->startRefPos --;
//					reg_new->startLocalRefPos --;
//					reg_new->startQueryPos --;
//				}
				computeVarType(reg_new);  // compute variant type
				reg_new->short_sv_flag = false;

				if(reg_new->endLocalRefPos-reg_new->startLocalRefPos>=MIN_DUP_SIZE or reg_new->endQueryPos-reg_new->startQueryPos>=MIN_DUP_SIZE){
					call_success = true;
					clip_reg_ret = reg_new;
					// replace items
//					for(size_t i=0; i<varVec.size(); i++) delete varVec.at(i);
//					vector<reg_t*>().swap(varVec);
//					varVec.push_back(reg_new);
				}else delete reg_new;
//			}else delete reg_new;

			//delete local_aln;
		}

	}else if(var_type==VAR_INV){ // inversion

	}else if(var_type==VAR_TRA){ // translocation

	}

	return clip_reg_ret;
}

vector<size_t> varCand::computeLeftShiftSizeDup(reg_t *reg, aln_seg_t *seg1, aln_seg_t *seg2, string &refseq, string &queryseq){
	int64_t startCheckLocalRefPos, startCheckRefPos, startCheckQueryPos, localRefPos, refPos, queryPos, misNum, refShiftSize, queryShiftSize, subseq_len;
	int64_t i, refPos_start, localRefPos_start, queryPos_start, tmp_len, minCheckRefPos, maxShiftRefSize, maxShiftQuerySize, shiftRefSize, shiftQuerySize, query_decrease_size, ref_decrease_size;
	int64_t ref_dist, query_dist, max_dist;
	localAln_t *local_aln;
	string refseq_aln, midseq_aln, queryseq_aln;
	vector<size_t> shift_size_vec;
	bool baseMatchFlag;

	minCheckRefPos = reg->startRefPos - 200;
	if(minCheckRefPos<1) minCheckRefPos = 1;

	ref_dist = reg->endRefPos - reg->startRefPos + 1;
	query_dist = reg->endQueryPos - reg->endQueryPos + 1;
	if(ref_dist<query_dist) max_dist = query_dist;
	else max_dist = ref_dist;

	if(max_dist<EXT_SIZE_CHK_VAR_LOC)
		subseq_len = EXT_SIZE_CHK_VAR_LOC * 2;  // 50
	else
		subseq_len = 2 * max_dist;

	if(seg1 and seg2){
		startCheckLocalRefPos = seg1->subject_end;
		startCheckRefPos = seg1->ref_end;
		startCheckQueryPos = seg2->query_start - 1;

		maxShiftRefSize = seg1->subject_end - 1;
		maxShiftQuerySize = seg1->query_end - 1;
	}else{
		startCheckLocalRefPos = reg->startLocalRefPos;
		startCheckRefPos = reg->startRefPos;
		startCheckQueryPos = reg->endQueryPos - 1;

		maxShiftRefSize = reg->startLocalRefPos - 1;
		maxShiftQuerySize = reg->startQueryPos - 1;
	}

	misNum = 0;
	localRefPos = startCheckLocalRefPos;
	refPos = startCheckRefPos;
	queryPos = startCheckQueryPos;
	shiftRefSize = shiftQuerySize = 0;

	if(maxShiftRefSize>0 and maxShiftQuerySize>0){
		while(queryPos>=1 and localRefPos>=1){
			baseMatchFlag = isBaseMatch(queryseq.at(queryPos-1), refseq.at(localRefPos-1));
			if(baseMatchFlag){
				query_decrease_size = ref_decrease_size = 1;
			}else{
				misNum = 0;

				if(localRefPos<subseq_len) subseq_len = localRefPos;
				if(queryPos<subseq_len) subseq_len = queryPos;

				localRefPos_start = localRefPos - subseq_len + 1;
				queryPos_start = queryPos - subseq_len + 1;
				refPos_start = refPos - subseq_len + 1;

				// compute local locations
				local_aln = new localAln_t();
				local_aln->reg = NULL;
				local_aln->blat_aln_id = -1;
				local_aln->aln_seg = local_aln->start_seg_extend = local_aln->end_seg_extend = NULL;
				//local_aln->startRefPos = local_aln->endRefPos = -1;
				//local_aln->startLocalRefPos = local_aln->startQueryPos = local_aln->endLocalRefPos = local_aln->endQueryPos = -1;
				local_aln->queryLeftShiftLen = local_aln->queryRightShiftLen = local_aln->localRefLeftShiftLen = local_aln->localRefRightShiftLen = -1;
				//local_aln->chrlen = 0;

				local_aln->startRefPos = refPos_start;
				local_aln->endRefPos = refPos;
				local_aln->startLocalRefPos = localRefPos_start;
				local_aln->endLocalRefPos = localRefPos;
				local_aln->startQueryPos = queryPos_start;
				local_aln->endQueryPos = queryPos;
				local_aln->chrlen = faidx_seq_len(fai, reg->chrname.c_str()); // get the reference length

				// get sub local sequence
				local_aln->refseq = refseq.substr(localRefPos_start-1, subseq_len);
				local_aln->ctgseq = queryseq.substr(queryPos_start-1, subseq_len);

				// pairwise alignment
				upperSeq(local_aln->ctgseq);
				upperSeq(local_aln->refseq);
				computeSeqAlignment(local_aln);

				// compute the number of mismatches
				misNum = computeMismatchNumLocalAln(local_aln);
				if(local_aln->localRefRightShiftLen<local_aln->queryRightShiftLen)
					misNum += local_aln->queryRightShiftLen;
				else
					misNum += local_aln->localRefRightShiftLen;

				if(misNum<=5){
					query_decrease_size = subseq_len - local_aln->queryLeftShiftLen;
					ref_decrease_size = subseq_len - local_aln->localRefLeftShiftLen;
					delete local_aln;
				}else {
					if(local_aln->localRefRightShiftLen<=2 and local_aln->queryRightShiftLen<=2){
						// compute the number of shift bases
						tmp_len = 0;
						midseq_aln = local_aln->alignResultVec.at(1);
						for(i=midseq_aln.size()-1; i>=0; i--){
							if(midseq_aln.at(i)=='|') tmp_len ++; // match
							else break; // mismatch
						}

						query_decrease_size = tmp_len + local_aln->queryRightShiftLen;
						ref_decrease_size = tmp_len + local_aln->localRefRightShiftLen;
					}
					delete local_aln;
					break;
				}
			}

			if(refPos>minCheckRefPos+ref_decrease_size and shiftRefSize+ref_decrease_size<maxShiftRefSize and shiftQuerySize+query_decrease_size<maxShiftQuerySize){
				queryPos -= query_decrease_size;
				localRefPos -= ref_decrease_size;
				refPos -= ref_decrease_size;
				shiftRefSize += ref_decrease_size;
				shiftQuerySize += query_decrease_size;
			}else break;
		}
	}

	refShiftSize = startCheckLocalRefPos - localRefPos;
	queryShiftSize = startCheckQueryPos - queryPos;
	shift_size_vec.push_back(refShiftSize);
	shift_size_vec.push_back(queryShiftSize);

	//cout << "refShiftSize=" << refShiftSize << ", queryShiftSize=" << queryShiftSize << ", left localRefPos=" << localRefPos << ", left refPos=" << refPos << ", left queryPos=" << queryPos << endl;

	return shift_size_vec;
}

vector<size_t> varCand::computeRightShiftSizeDup(reg_t *reg, aln_seg_t *seg1, aln_seg_t *seg2, string &refseq, string &queryseq){
	int64_t startCheckLocalRefPos, startCheckRefPos, startCheckQueryPos, localRefPos, refPos, queryPos, misNum, refShiftSize, queryShiftSize, subseq_len;
	int64_t i, refPos_start, localRefPos_start, queryPos_start, tmp_len, maxCheckRefPos, chrlen_tmp, maxShiftRefSize, maxShiftQuerySize, shiftRefSize, shiftQuerySize, query_increase_size, ref_increase_size;
	int64_t ref_dist, query_dist, max_dist;
	localAln_t *local_aln;
	string midseq_aln;
	vector<size_t> shift_size_vec;
	bool baseMatchFlag;

	chrlen_tmp = faidx_seq_len(fai, reg->chrname.c_str()); // get the reference length
	maxCheckRefPos = reg->endRefPos + 200;
	if(maxCheckRefPos>chrlen_tmp) maxCheckRefPos = chrlen_tmp;

	ref_dist = reg->endRefPos - reg->startRefPos + 1;
	query_dist = reg->endQueryPos - reg->endQueryPos + 1;
	if(ref_dist<query_dist) max_dist = query_dist;
	else max_dist = ref_dist;

	if(max_dist<EXT_SIZE_CHK_VAR_LOC)
		subseq_len = EXT_SIZE_CHK_VAR_LOC * 2;  // 50
	else
		subseq_len = 2 * max_dist;

	if(seg1 and seg2){
		startCheckLocalRefPos = seg2->subject_start;
		startCheckRefPos = seg2->ref_start;
		startCheckQueryPos = seg1->query_end + 1;

		maxShiftRefSize = refseq.size() - seg2->subject_start;
		maxShiftQuerySize = queryseq.size() - seg2->query_start;
	}else{
		startCheckLocalRefPos = reg->endLocalRefPos;
		startCheckRefPos = reg->endRefPos;
		startCheckQueryPos = reg->startQueryPos + 1;

		maxShiftRefSize = refseq.size() - reg->endLocalRefPos;
		maxShiftQuerySize = queryseq.size() - reg->endQueryPos;
	}

	misNum = 0;
	localRefPos = startCheckLocalRefPos;
	refPos = startCheckRefPos;
	queryPos = startCheckQueryPos;
	shiftRefSize = shiftQuerySize = 0;

	if(maxShiftRefSize>0 and maxShiftQuerySize>0){
		while(queryPos<=(int64_t)queryseq.size() and localRefPos<=(int64_t)refseq.size()){
			baseMatchFlag = isBaseMatch(queryseq.at(queryPos-1), refseq.at(localRefPos-1));
			if(baseMatchFlag){
				query_increase_size = ref_increase_size = 1;
			}else{
				misNum = 0;

				localRefPos_start = localRefPos;
				queryPos_start = queryPos;
				refPos_start = refPos;

				tmp_len = refseq.size() - localRefPos;
				if(tmp_len<subseq_len) subseq_len = tmp_len;
				tmp_len = queryseq.size() - queryPos;
				if(tmp_len<subseq_len) subseq_len = tmp_len;

				// compute local locations
				local_aln = new localAln_t();
				local_aln->reg = NULL;
				local_aln->blat_aln_id = -1;
				local_aln->aln_seg = local_aln->start_seg_extend = local_aln->end_seg_extend = NULL;
				//local_aln->startRefPos = local_aln->endRefPos = -1;
				//local_aln->startLocalRefPos = local_aln->startQueryPos = local_aln->endLocalRefPos = local_aln->endQueryPos = -1;
				local_aln->queryLeftShiftLen = local_aln->queryRightShiftLen = local_aln->localRefLeftShiftLen = local_aln->localRefRightShiftLen = -1;
				//local_aln->chrlen = 0;

				local_aln->startRefPos = refPos_start;
				local_aln->endRefPos = refPos_start + subseq_len - 1;
				local_aln->startLocalRefPos = localRefPos_start;
				local_aln->endLocalRefPos = localRefPos_start + subseq_len - 1;
				local_aln->startQueryPos = queryPos_start;
				local_aln->endQueryPos = queryPos_start + subseq_len - 1;
				local_aln->chrlen = chrlen_tmp;

				local_aln->refseq = refseq.substr(localRefPos_start-1, subseq_len);
				local_aln->ctgseq = queryseq.substr(queryPos_start-1, subseq_len);

				// pairwise alignment
				upperSeq(local_aln->ctgseq);
				upperSeq(local_aln->refseq);
				computeSeqAlignment(local_aln);

				// compute the number of mismatches
				misNum = computeMismatchNumLocalAln(local_aln);
				if(local_aln->localRefLeftShiftLen<local_aln->queryLeftShiftLen)
					misNum += local_aln->queryLeftShiftLen;
				else
					misNum += local_aln->localRefLeftShiftLen;

				if(misNum<=5){
					query_increase_size = subseq_len - local_aln->queryRightShiftLen;
					ref_increase_size = subseq_len - local_aln->localRefRightShiftLen;
					delete local_aln;
				}else{
					if(local_aln->localRefLeftShiftLen<=2 and local_aln->queryLeftShiftLen<=2){
						// compute the number of shift bases
						tmp_len = 0;
						midseq_aln = local_aln->alignResultVec.at(1);
						for(i=0; i<(int64_t)midseq_aln.size(); i++){
							if(midseq_aln.at(i)=='|') tmp_len ++; // match
							else break; // mismatch
						}

						query_increase_size = tmp_len + local_aln->queryLeftShiftLen;
						ref_increase_size = tmp_len + local_aln->localRefLeftShiftLen;
					}
					delete local_aln;
					break;
				}
			}

			if(refPos+ref_increase_size<maxCheckRefPos and shiftRefSize+ref_increase_size<maxShiftRefSize and shiftQuerySize+query_increase_size<maxShiftQuerySize) {
				queryPos += query_increase_size;
				localRefPos +=  ref_increase_size;
				refPos += ref_increase_size;
				shiftRefSize += ref_increase_size;
				shiftQuerySize += query_increase_size;
			}else break;
		}
	}

	refShiftSize = localRefPos - startCheckLocalRefPos;
	queryShiftSize = queryPos - startCheckQueryPos;
	shift_size_vec.push_back(refShiftSize);
	shift_size_vec.push_back(queryShiftSize);

	//cout << "refShiftSize=" << refShiftSize << ", queryShiftSize=" << queryShiftSize << ", right localRefPos=" << localRefPos << ", right refPos=" << refPos << ", right queryPos=" << queryPos << endl;

	return shift_size_vec;
}


size_t varCand::computeMismatchNumLocalAln(localAln_t *local_aln){
	size_t i, misNum = 0;
	string midseq;
	midseq = local_aln->alignResultVec.at(1);
	for(i=0; i<midseq.size(); i++) if(midseq.at(i)!='|') misNum ++;
	return misNum;
}

// compute query position
vector<size_t> varCand::computeQueryClipPosDup(blat_aln_t *blat_aln, int32_t clipRefPos, string &refseq, string &queryseq){
	size_t i;
	int32_t localRefPos_start, refPos_start, queryPos_start;
	int32_t dist, queryPos, refPos, localRefPos, sub_refseq_len, sub_queryseq_len;
	localAln_t *local_aln;
	aln_seg_t *aln_seg;
	string query_aln_seq, mid_aln_seq, ref_aln_seq;
	vector<size_t> pos_vec;

	// get aln_seg
	for(i=0; i<blat_aln->aln_segs.size(); i++){
		aln_seg = blat_aln->aln_segs[i];
		if(aln_seg->ref_start<=clipRefPos and aln_seg->ref_end>=clipRefPos) break;
	}

	if(i<blat_aln->aln_segs.size()){
		dist = clipRefPos - aln_seg->ref_start;
		refPos_start = clipRefPos - 200;
		queryPos_start = aln_seg->query_start + dist - 200;
		localRefPos_start = aln_seg->subject_start + dist - 200;
		if(queryPos_start<1) queryPos_start = 1;
		if(localRefPos_start<1) localRefPos_start = 1;

		if(queryPos_start<1 or queryPos_start>(int32_t)queryseq.size()){
			cerr << "queryPos_start=" << queryPos_start << ", queryseq_len=" << queryseq.size() << ", error!" << endl;
			exit(1);
		}
		if(localRefPos_start<1 or localRefPos_start>(int32_t)refseq.size()){
			cerr << "localRefPos_start=" << localRefPos_start << ", refseq_len=" << refseq.size() << ", error!" << endl;
			exit(1);
		}

		// compute local locations
		local_aln = new localAln_t();
		local_aln->reg = NULL;
		local_aln->blat_aln_id = -1;
		local_aln->aln_seg = local_aln->start_seg_extend = local_aln->end_seg_extend = NULL;
		local_aln->startRefPos = local_aln->endRefPos = -1;
		local_aln->startLocalRefPos = local_aln->startQueryPos = local_aln->endLocalRefPos = local_aln->endQueryPos = -1;
		local_aln->queryLeftShiftLen = local_aln->queryRightShiftLen = local_aln->localRefLeftShiftLen = local_aln->localRefRightShiftLen = -1;
		local_aln->chrlen = 0;

		//sub_refseq_len = sub_queryseq_len = 400; // MIN_AVER_SIZE_ALN_SEG;
		sub_refseq_len = sub_queryseq_len = EXT_SIZE_CHK_VAR_LOC; // MIN_AVER_SIZE_ALN_SEG;
		if(localRefPos_start-1+sub_refseq_len>(int32_t)refseq.size()) sub_refseq_len = refseq.size() - (localRefPos_start - 1);
		if(queryPos_start-1+sub_queryseq_len>(int32_t)queryseq.size()) sub_queryseq_len = queryseq.size() - (queryPos_start - 1);

		local_aln->refseq = refseq.substr(localRefPos_start-1, sub_refseq_len);
		local_aln->ctgseq = queryseq.substr(queryPos_start-1, sub_queryseq_len);

		// pairwise alignment
		upperSeq(local_aln->ctgseq);
		upperSeq(local_aln->refseq);
		computeSeqAlignment(local_aln);

		query_aln_seq = local_aln->alignResultVec.at(0);
		mid_aln_seq = local_aln->alignResultVec.at(1);
		ref_aln_seq = local_aln->alignResultVec.at(2);

		// compute query position
		queryPos = queryPos_start + local_aln->queryLeftShiftLen;
		refPos = refPos_start + local_aln->localRefLeftShiftLen;
		localRefPos = localRefPos_start + local_aln->localRefLeftShiftLen;
		for(i=0; i<mid_aln_seq.size(); i++){
			if(mid_aln_seq.at(i)=='|'){ // match
				queryPos ++; refPos ++; localRefPos ++;
			}else{
				if(query_aln_seq.at(i)=='-'){
					refPos ++; localRefPos ++;
				}else if(ref_aln_seq.at(i)=='-'){
					queryPos ++;
				}else{
					queryPos ++; refPos ++; localRefPos ++;
				}
			}
			if(refPos==clipRefPos) break;
		}
		refPos--; localRefPos--; queryPos--; // back 1 base
		pos_vec.push_back(refPos);
		pos_vec.push_back(localRefPos);
		pos_vec.push_back(queryPos);

		delete local_aln;
	}

	return pos_vec;
}

// generate new local align item with only alignment information
localAln_t* varCand::generateNewLocalAlnItem_OnlyAlnInfo(localAln_t *local_aln){
	localAln_t *local_aln_new = NULL;

	if(local_aln){
		local_aln_new = new localAln_t();
		local_aln_new->reg = local_aln_new->cand_reg = NULL;
		local_aln_new->blat_aln_id = -1;
		local_aln_new->aln_seg = local_aln_new->start_seg_extend = local_aln_new->end_seg_extend = NULL;

		local_aln_new->startRefPos = local_aln->startRefPos;
		local_aln_new->endRefPos = local_aln->endRefPos;
		local_aln_new->startLocalRefPos = local_aln->startLocalRefPos;
		local_aln_new->endLocalRefPos = local_aln->endLocalRefPos;
		local_aln_new->startQueryPos = local_aln->startQueryPos;
		local_aln_new->endQueryPos = local_aln->endQueryPos;
		local_aln_new->queryLeftShiftLen = local_aln->queryLeftShiftLen;
		local_aln_new->queryRightShiftLen = local_aln->queryRightShiftLen;
		local_aln_new->localRefLeftShiftLen = local_aln->localRefLeftShiftLen;
		local_aln_new->localRefRightShiftLen = local_aln->localRefRightShiftLen;
		local_aln_new->chrlen = local_aln->chrlen;
		local_aln_new->ctgseq = local_aln->ctgseq;
		local_aln_new->refseq = local_aln->refseq;
		local_aln_new->alignResultVec.push_back(local_aln->alignResultVec.at(0));
		local_aln_new->alignResultVec.push_back(local_aln->alignResultVec.at(1));
		local_aln_new->alignResultVec.push_back(local_aln->alignResultVec.at(2));
		local_aln_new->overlapLen = local_aln->overlapLen;
		local_aln_new->start_aln_idx_var = local_aln->start_aln_idx_var;
		local_aln_new->end_aln_idx_var = local_aln->end_aln_idx_var;
	}

	return local_aln_new;
}

// copy local alignment information
void varCand::copyLocalAlnInInfo(localAln_t *local_aln_dest, localAln_t *local_aln_src){
	if(local_aln_dest and local_aln_src){
//		local_aln_dest->startRefPos = local_aln_src->startRefPos;
//		local_aln_dest->endRefPos = local_aln_src->endRefPos;
//		local_aln_dest->startLocalRefPos = local_aln_src->startLocalRefPos;
//		local_aln_dest->endLocalRefPos = local_aln_src->endLocalRefPos;
//		local_aln_dest->startQueryPos = local_aln_src->startQueryPos;
//		local_aln_dest->endQueryPos = local_aln_src->endQueryPos;
		local_aln_dest->queryLeftShiftLen = local_aln_src->queryLeftShiftLen;
		local_aln_dest->queryRightShiftLen = local_aln_src->queryRightShiftLen;
		local_aln_dest->localRefLeftShiftLen = local_aln_src->localRefLeftShiftLen;
		local_aln_dest->localRefRightShiftLen = local_aln_src->localRefRightShiftLen;
		//local_aln_dest->chrlen = local_aln_src->chrlen;
		//local_aln_dest->ctgseq = local_aln_src->ctgseq;
		//local_aln_dest->refseq = local_aln_src->refseq;
		local_aln_dest->alignResultVec.push_back(local_aln_src->alignResultVec.at(0));
		local_aln_dest->alignResultVec.push_back(local_aln_src->alignResultVec.at(1));
		local_aln_dest->alignResultVec.push_back(local_aln_src->alignResultVec.at(2));
		local_aln_dest->overlapLen = local_aln_src->overlapLen;
		local_aln_dest->start_aln_idx_var = local_aln_src->start_aln_idx_var;
		local_aln_dest->end_aln_idx_var = local_aln_src->end_aln_idx_var;
	}else{
		cerr << __func__ << ": local_aln_dest=" << local_aln_dest << ", local_aln_src=" << local_aln_src << ", invalid." << endl;
		exit(1);
	}
}

// determine whether the information of the given local alignment item is complete
bool varCand::isLocalAlnInfoComplete(localAln_t *local_aln){
	bool flag = true;

	if(local_aln){
		if(local_aln->startRefPos<=0 or local_aln->endRefPos<=0 or local_aln->startLocalRefPos<=0 or local_aln->endLocalRefPos<=0
				or local_aln->startQueryPos<=0 or local_aln->endQueryPos<=0 or local_aln->chrlen<=0
				or local_aln->ctgseq.size()==0 or local_aln->refseq.size()==0)
			flag = false;
	}else flag = false;

	return flag;
}

// determine whether the two align items are identical
bool varCand::isIdenticalLoclaAlnItems(localAln_t *local_aln1, localAln_t *local_aln2){
	bool flag = true;

	if(local_aln1 and local_aln2){
		if(local_aln1->startRefPos!=local_aln2->startRefPos or local_aln1->endRefPos!=local_aln2->endRefPos
				or local_aln1->startLocalRefPos!=local_aln2->startLocalRefPos or local_aln1->endLocalRefPos!=local_aln2->endLocalRefPos
				or local_aln1->startQueryPos!=local_aln2->startQueryPos or local_aln1->endQueryPos!=local_aln2->endQueryPos
				or local_aln1->chrlen!=local_aln2->chrlen
				or local_aln1->ctgseq.compare(local_aln2->ctgseq)!=0 or local_aln1->refseq.compare(local_aln2->refseq)!=0){
			flag = false;
		}
	}else if((local_aln1==NULL and local_aln2) or (local_aln1 and local_aln2==NULL)) flag = false;

	return flag;
}

// get identical local align item form vector
localAln_t* varCand::getIdenticalLocalAlnItem(localAln_t *local_aln, vector<localAln_t*> &local_aln_vec){
	localAln_t *local_aln_ret = NULL, *local_aln_tmp;

	if(local_aln){
		for(size_t i=0; i<local_aln_vec.size(); i++){
			local_aln_tmp = local_aln_vec.at(i);
			if(isIdenticalLoclaAlnItems(local_aln, local_aln_tmp)){
				local_aln_ret = local_aln_tmp;
				break;
			}
		}
	}

	return local_aln_ret;
}

// add local alignment item to vector
void varCand::addLocalAlnItemToVec(localAln_t *local_aln, vector<localAln_t*> &local_aln_vec){
	if(local_aln) local_aln_vec.push_back(local_aln);
	else{
		cerr << __func__ << ": local_aln=" << local_aln << ", invalid." << endl;
		exit(1);
	}
}

// destroy local alignment items
void varCand::destroyLocalAlnVec(vector<localAln_t*> &local_aln_vec){
	localAln_t *local_aln;
	for(size_t i=0; i<local_aln_vec.size(); i++){
		local_aln = local_aln_vec.at(i);
		delete local_aln;
	}
	vector<localAln_t*>().swap(local_aln_vec);
}


void varCand::printSV(){
	size_t j;
	reg_t *reg;
	string line, sv_type;

	for(j=0; j<varVec.size(); j++){
		reg = varVec[j];
		switch(reg->var_type){
			case VAR_UNC: sv_type = VAR_UNC_STR; break;
			case VAR_INS: sv_type = VAR_INS_STR; break;
			case VAR_DEL: sv_type = VAR_DEL_STR; break;
			case VAR_DUP: sv_type = VAR_DUP_STR; break;
			case VAR_INV: sv_type = VAR_INV_STR; break;
			case VAR_TRA: sv_type = VAR_TRA_STR; break;
			default: sv_type = VAR_MIX_STR; break;
		}
		line = reg->chrname + "\t" + to_string(reg->startRefPos) + "\t" + to_string(reg->endRefPos) + "\t" + reg->refseq + "\t" + reg->altseq + "\t" + sv_type;
		if(reg->var_type==VAR_INS or reg->var_type==VAR_DEL or reg->var_type==VAR_DUP)
			line += "\t" + to_string(reg->sv_len);
		else
			line += "\t.";
		cout << line << endl;
	}
}

// output newVarVector
void varCand::outputNewVarcandVec(){
	reg_t *reg;

	if(newVarVec.size()>0){
		cout << alnfilename << endl;
		for(size_t i=0; i<newVarVec.size(); i++){
			reg = newVarVec[i];
			cout << "\t" << reg->chrname << ":" << reg->startRefPos << "-" << reg->endRefPos << endl;
		}
	}
}

// output multiple blat align blocks
void varCand::outputMultiBlatAlnBlocks(){
	if(blat_aln_vec.size()>=2)
		cout << alnfilename << ": blat align blocks: " << blat_aln_vec.size() << endl;
	else if(blat_aln_vec.size()==0)
		cout << alnfilename << ": no blat align blocks" << endl;
}
