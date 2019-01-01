#include <math.h>
#include "varCand.h"
#include "util.h"

pthread_mutex_t mutex_print_var_cand = PTHREAD_MUTEX_INITIALIZER;

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
}

void varCand::init(){
	refseqfilename = "";
	ctgfilename = "";
	readsfilename = "";
	alnfilename = "";
	ref_left_shift_size = ref_right_shift_size = ctg_num = 0;
	assem_success = align_success = call_success = clip_reg_flag = false;
	clip_reg = NULL;
	leftClipRefPos = rightClipRefPos = 9;
	sv_type = VAR_UNC;
	dup_num = 0;
	margin_adjusted_flag = false;
}

// align contig to refseq
void varCand::alnCtg2Refseq(){
	//cout << "BLAT align: " << ctgfilename << endl;
	if(assem_success)
		blatAln(alnfilename, ctgfilename, refseqfilename); // BLAT alignment
}

// BLAT alignment, and the output is in sim4 format
void varCand::blatAln(string &alnfilename, string &contigfilename, string &refseqfilename){
	string blat_cmd, out_opt;
	int ret;

	out_opt = "-out=sim4 " + alnfilename;
	blat_cmd = "blat " + refseqfilename + " " + contigfilename + " " + out_opt + " > /dev/null 2>&1";

	//cout << "blat_cmd: " + blat_cmd << endl;

	ret = system(blat_cmd.c_str());
	if(ret!=0){
		cout << "Please run the correct blat command or check whether blat was correctly installed." << endl;
		exit(1);
	}
}

// call variants according to BLAT alignment
void varCand::callVariants(){
	if(sv_type!=VAR_TRA){
		pthread_mutex_lock(&mutex_print_var_cand);
		cout << "process: " << ctgfilename << endl;
		pthread_mutex_unlock(&mutex_print_var_cand);

		loadBlatAlnData(); // load align data

		if(clip_reg_flag==false)	// call indel variants
			callIndelVariants();
		else // call clipping variants
			callClipRegVariants();
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
		cout << "Align Failure: " << ctgfilename << endl;
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
	bool flag = false;
	struct stat fileStat;
	if (stat(alnfilename.c_str(), &fileStat) == 0)
		if(fileStat.st_size>0)
			flag = true;
	align_success = flag;
}

// parse BLAT alignments
void varCand::blatParse(){
	string line, q_reg_str, s_reg_str, ident_percent_str, orient_str, query_name;
	vector<string> line_vec, line_vec1, len_vec, str_vec, query_name_vec;
	ifstream infile;
	int32_t ref_start_all, query_loc, blat_aln_id;
	blat_aln_t *blat_aln_item;
	aln_seg_t *aln_seg;

	infile.open(alnfilename);
	if(!infile.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << alnfilename << endl;
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

				//line_vec = split(line, ",");
				len_vec = split(line_vec1[line_vec1.size()-1].substr(1), " ");
				if(line_vec[0][3]=='1') {
					// get the query name
					//best_aln = false;
					query_loc = getQueryNameLoc(line_vec1[0], query_name_vec);
					if(query_loc==-1){
						query_name_vec.push_back(line_vec1[0]);
						query_loc = query_name_vec.size() - 1;
						//best_aln = true;
					}

					blat_aln_item = new blat_aln_t;
					blat_aln_item->blat_aln_id = blat_aln_id++;
					blat_aln_item->query_id = query_loc;
					blat_aln_item->query_len = stoi(len_vec[0]); // query length
					blat_aln_item->head_hanging = blat_aln_item->tail_hanging = false;
					//blat_aln_item->best_aln = best_aln;
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
				q_reg_str = line_vec[0];   // query
				s_reg_str = line_vec[2];   // subject
				ident_percent_str = line_vec[5];  // identity percent
				orient_str = line_vec[6];  // orientation

				aln_seg = new aln_seg_t;   // allocate memory

				str_vec = split(q_reg_str, "-");
				aln_seg->query_start = stoi(str_vec[0]);  // query start
				aln_seg->query_end = stoi(str_vec[1]);  // query end

				str_vec = split(s_reg_str, "-");
				aln_seg->subject_start = stoi(str_vec[0].substr(1));  // subject start
				aln_seg->subject_end = stoi(str_vec[1].substr(0, str_vec[1].size()-1));  // subject end

				aln_seg->ident_percent = stof(ident_percent_str) / 100.0;  // identity percent

				if(orient_str.compare("->")==0) aln_seg->aln_orient = ALN_PLUS_ORIENT;  // orientation
				else aln_seg->aln_orient = ALN_MINUS_ORIENT;

				// compute the reference position
				aln_seg->ref_start = ref_start_all + aln_seg->subject_start - 1;
				aln_seg->ref_end = ref_start_all + aln_seg->subject_end - 1;

				blat_aln_item->aln_segs.push_back(aln_seg);
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

// get query name location
int32_t varCand::getQueryNameLoc(string query_name, vector<string> query_name_vec){
	int32_t i, loc = -1;
	for(i=query_name_vec.size()-1; i>=0; i--)
		if(query_name_vec[i].compare(query_name)==0) { loc = i; break; }
	return loc;
}

// assign best_aln flag
void varCand::assignBestAlnFlag(){
	size_t best_aln_id, query_id, query_num, sum, max_sum, max_sum_id, secMax_sum, secMax_sum_id, valid_seg_num, max_valid_seg_num, secMax_valid_seg_num, tmp_len;
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
	int32_t aln_orient = -1;
	blat_aln_t *blat_aln;
	for(size_t i=0; i<blat_aln_vec.size(); i++){
		blat_aln = blat_aln_vec.at(i);
		if(blat_aln->best_aln==true)
			aln_orient = blat_aln->aln_orient;
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

	// remove invalid blat align segments
	//removeInvalidAlnSegs();
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
	size_t i, j, validNum, best_aln_num, start_query_pos, end_query_pos;
	int8_t **query_cov_array;
	blat_aln_t *blat_aln;
	aln_seg_t *aln_seg;
	float covRatio;

	validNum = 0;
	for(i=0; i<blat_aln_vec.size(); i++) if(blat_aln_vec[i]->valid_aln) validNum ++;
	if(validNum>(size_t)ctg_num and ctg_num>0){
		query_cov_array = new int8_t *[ctg_num]();
		// get the number of best blat alignments
		best_aln_num = 0;
		for(i=0; i<blat_aln_vec.size(); i++){
			blat_aln = blat_aln_vec[i];
			if(blat_aln->best_aln){
				query_cov_array[best_aln_num] = new int8_t[blat_aln->query_len+1]();
				// compute the query coverage based on the first blat alignment block for each query
				for(j=0; j<blat_aln->aln_segs.size(); j++){
					aln_seg = blat_aln->aln_segs[j];
					if(aln_seg->aln_orient==ALN_PLUS_ORIENT){
						start_query_pos = aln_seg->query_start;
						end_query_pos = aln_seg->query_end;
					}else{
						start_query_pos = blat_aln->query_len - aln_seg->query_end + 1;
						end_query_pos = blat_aln->query_len - aln_seg->query_start + 1;
					}
					memset(query_cov_array[best_aln_num]+start_query_pos-1, 1, sizeof(int8_t)*((int32_t)end_query_pos-start_query_pos+1));
				}
				best_aln_num ++;
			}
		}

		// compute the coverage of other blat alignment blocks
		for(i=0; i<blat_aln_vec.size(); i++){
			blat_aln = blat_aln_vec[i];
			if(blat_aln->valid_aln and !blat_aln->best_aln){
				covRatio = computeRepeatCovRatio(blat_aln, query_cov_array[blat_aln->query_id], true); // compute the coverage
				if(covRatio>=INVALID_COV_RATIO_THRES) blat_aln->valid_aln = false;
				else updateCovArray(blat_aln, query_cov_array[blat_aln->query_id], true); // update the coverage
			}
		}

		for(i=0; i<(size_t)ctg_num; i++) delete[] query_cov_array[i];
		delete[] query_cov_array;
	}
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
					aln_seg = blat_aln->aln_segs[j];
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
			aln_seg = blat_aln->aln_segs[i];
			memset(subject_cov_array+aln_seg->subject_start-1, 1, sizeof(int8_t)*((int32_t)aln_seg->subject_end-aln_seg->subject_start+1));
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
		aln_seg = blat_aln->aln_segs[i];
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

// remove invalid blat align segments
//void varCand::removeInvalidAlnSegs(){
//	blat_aln_t* blat_item;
//	for(size_t i=0; i<blat_aln_vec.size(); ){
//		blat_item = blat_aln_vec.at(i);
//		if(!blat_item->valid_aln) {  // remove blat align segments
//			vector<aln_seg_t*>::iterator seg;
//			for(seg=blat_item->aln_segs.begin(); seg!=blat_item->aln_segs.end(); seg++) // free aln_segs
//				delete *seg;
//			delete blat_item;
//			blat_aln_vec.erase(blat_aln_vec.begin()+i);
//		}else i ++;
//	}
//	blat_aln_vec.shrink_to_fit();
//}

// call indel types
void varCand::determineIndelType(){
	int32_t  i, j, query_dist, ref_dist, decrease_size, var_type, startRefPos, endRefPos, startLocalRefPos, endLocalRefPos, startQueryPos, endQueryPos;
	blat_aln_t *blat_aln;
	aln_seg_t *seg1, *seg2;
	reg_t *reg, *reg_tmp;
	string reg_str;
	bool newFlag, updateFlag;
	vector<reg_t*> foundRegVec, candRegVec;

	for(i=0; i<(int32_t)blat_aln_vec.size(); i++){
		var_type = VAR_UNC;
		blat_aln = blat_aln_vec[i];
		if(blat_aln->valid_aln==true and blat_aln->aln_segs.size()>1){
			for(j=0; j<(int32_t)blat_aln->aln_segs.size()-1; j++){
				seg1 = blat_aln->aln_segs[j];
				seg2 = blat_aln->aln_segs[j+1];
				if(seg1->query_end>CTG_END_SKIP_SIZE and seg2->query_start>CTG_END_SKIP_SIZE){  // skip contig ends
					query_dist = (int32_t)seg2->query_start - (int32_t)seg1->query_end;
					ref_dist = (int32_t)seg2->ref_start - (int32_t)seg1->ref_end;
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
					reg = findVarvecItem(startRefPos, endRefPos, varVec);
					if(reg){ // found, true positive, and then update its information
						if(reg->blat_aln_id==-1)
							reg->blat_aln_id = i;
						if(reg->blat_aln_id==i)
						{
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
							reg_tmp->aln_orient = blat_aln->aln_orient;
							candRegVec.push_back(reg_tmp);
						}
					}else { // if it is a false negative discovery, search the newVarVec
						reg = findVarvecItem(startRefPos, endRefPos, newVarVec);
						if(reg==NULL) {  // if not exist, add it to newVarVec
							reg = new reg_t;
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

	// release candRegVec
	vector<reg_t*>::iterator it;
	for(it=candRegVec.begin(); it!=candRegVec.end(); it++)
		delete (*it);
	vector<reg_t*>().swap(candRegVec);
}

// call short variants that cannot be split at the variant location
void varCand::callShortVariants(){
	reg_t *reg;
	aln_seg_t *aln_seg;
	//string refseq_aln, ctgseq_aln;
	localAln_t *local_aln;

	for(size_t i=0; i<varVec.size(); i++){
		reg = varVec[i];
		if(reg->call_success_status==false){
			if(isUnsplitAln(reg)){
				reg->short_sv_flag = true;
				if(determineQueryidReg(reg)){ // determine query_id
					// get the aln_seg
					aln_seg = getAlnSeg(reg);
					if(aln_seg){
						// compute local locations
						local_aln = new localAln_t;
						local_aln->reg = reg;
						local_aln->aln_seg = aln_seg;
						local_aln->startRefPos = local_aln->endRefPos = -1;
						local_aln->startLocalRefPos = local_aln->startQueryPos = local_aln->endLocalRefPos = local_aln->endQueryPos = -1;
						local_aln->queryLeftShiftLen = local_aln->queryRightShiftLen = local_aln->localRefLeftShiftLen = local_aln->localRefRightShiftLen = -1;
						computeLocalLocsAlnShortVar(local_aln);

						// get sub local reference sequence
						FastaSeqLoader refseqloader(refseqfilename);
						local_aln->refseq = refseqloader.getFastaSeqByPos(0, local_aln->startLocalRefPos, local_aln->endLocalRefPos, ALN_PLUS_ORIENT);
						if(reg->aln_orient==ALN_MINUS_ORIENT) reverseComplement(local_aln->refseq);
						FastaSeqLoader ctgseqloader(ctgfilename);
						local_aln->ctgseq = ctgseqloader.getFastaSeqByPos(reg->query_id, local_aln->startQueryPos, local_aln->endQueryPos, reg->aln_orient);

						// pairwise alignment
						upperSeq(local_aln->ctgseq);
						upperSeq(local_aln->refseq);
						computeSeqAlignment(local_aln);

						// compute variant location
						computeVarLoc(local_aln);

						// adjust the variant locations according to alignment
						adjustVarLoc(local_aln);

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
aln_seg_t* varCand::getAlnSeg(reg_t *reg){
	blat_aln_t *blat_aln;
	aln_seg_t *aln_seg, *target;

	target = NULL;
	if(reg->query_id!=-1){
		for(size_t i=0; i<blat_aln_vec.size(); i++){
			blat_aln = blat_aln_vec.at(i);
			if(blat_aln->valid_aln==true and blat_aln->query_id==(size_t)reg->query_id){
				for(size_t j=0; j<blat_aln->aln_segs.size(); j++){
					aln_seg = blat_aln->aln_segs[j];
					if(reg->startRefPos>aln_seg->ref_start and reg->endRefPos<aln_seg->ref_end){
						target = aln_seg;
						break;
					}
				}
			}
		}
	}
	return target;
}

// compute local locations
void varCand::computeLocalLocsAlnShortVar(localAln_t *local_aln){
	reg_t *reg = local_aln->reg;
	aln_seg_t *aln_seg = local_aln->aln_seg;
	int32_t left_dist, right_dist;

	if(aln_seg->ref_start+VAR_ALN_EXTEND_SIZE<reg->startRefPos)
		left_dist = reg->startRefPos - aln_seg->ref_start - VAR_ALN_EXTEND_SIZE;
	else
		left_dist = 0;
	local_aln->startRefPos = aln_seg->ref_start + left_dist;
	local_aln->startLocalRefPos = aln_seg->subject_start + left_dist;
	local_aln->startQueryPos = aln_seg->query_start + left_dist;

	if(reg->endRefPos+VAR_ALN_EXTEND_SIZE<aln_seg->ref_end)
		right_dist = aln_seg->ref_end - reg->endRefPos - VAR_ALN_EXTEND_SIZE;
	else
		right_dist = 0;
	local_aln->endRefPos = aln_seg->ref_end - right_dist;
	local_aln->endLocalRefPos = aln_seg->subject_end - right_dist;
	local_aln->endQueryPos = aln_seg->query_end - right_dist;
}

// compute sequence alignment (LCS)
void varCand::computeSeqAlignment(localAln_t *local_aln){
	int32_t i, j, matchScore, mismatchScore, gapScore, gapOpenScore, tmp_gapScore1, tmp_gapScore2, maxValue, scoreIJ;
	int32_t rowsNum = local_aln->ctgseq.size() + 1, colsNum = local_aln->refseq.size() + 1;
	int32_t *scoreArr, *pathArr, path_val, arrSize;
	int32_t maxValueLastRow, maxValueLastCol, maxCol, maxRow;
	int32_t mismatchNum, itemNum;
	string queryAlnResult, midAlnResult, subjectAlnResult;

	matchScore = MATCH_SCORE;
	mismatchScore = MISMATCH_SCORE;
	gapScore = GAP_SCORE;
	gapOpenScore = GAP_OPEN_SCORE;

	arrSize = rowsNum * colsNum;
	scoreArr = (int32_t*) calloc (arrSize, sizeof(int32_t));
	pathArr = (int32_t*) calloc (arrSize, sizeof(int32_t));

	for(i=0; i<rowsNum; i++) scoreArr[i*colsNum] = 0;
	for(j=1; j<colsNum; j++) scoreArr[j] = 0;

	for(i=0; i<rowsNum; i++) for(j=0; j<colsNum; j++) pathArr[i*colsNum+j] = -1;

	// compute the scores of each element
	for(i=1; i<rowsNum; i++){
		for(j=1; j<colsNum; j++){
			if(local_aln->ctgseq[i-1]==local_aln->refseq[j-1])
				scoreIJ = matchScore;
			else
				scoreIJ = mismatchScore;

			if(pathArr[(i-1)*colsNum+j]!=1)
				tmp_gapScore1 = gapOpenScore;
			else
				tmp_gapScore1 = gapScore;

			if(pathArr[i*colsNum+j-1]!=2)
				tmp_gapScore2 = gapOpenScore;
			else
				tmp_gapScore2 = gapScore;

			maxValue = INT_MIN;
			path_val = -1;
			// compute the maximal score
			if(scoreArr[(i-1)*colsNum+j-1]+scoreIJ>maxValue) {// from (i-1, j-1)
				maxValue = scoreArr[(i-1)*colsNum+j-1] + scoreIJ;
				path_val = 0;
			}
			if(scoreArr[(i-1)*colsNum+j]+tmp_gapScore1>maxValue) {// from (i-1, j)
				maxValue = scoreArr[(i-1)*colsNum+j] + tmp_gapScore1;
				path_val = 1;
			}
			if(scoreArr[i*colsNum+j-1]+tmp_gapScore2>maxValue) {// from (i, j-1)
				maxValue = scoreArr[i*colsNum+j-1] + tmp_gapScore2;
				path_val = 2;
			}

			scoreArr[i*colsNum+j] = maxValue;
			pathArr[i*colsNum+j] = path_val;
		}
	}

	// get the row and col of the maximal element in last row and column
	maxValueLastRow = INT_MIN;
	for(j=0; j<colsNum; j++)
		if(scoreArr[(rowsNum-1)*colsNum+j]>maxValueLastRow){
			maxValueLastRow = scoreArr[(rowsNum-1)*colsNum+j];
			maxCol = j;
		}
	maxValueLastCol = INT_MIN;
	for(i=0; i<rowsNum; i++)
		if(scoreArr[i*colsNum+colsNum-1]>maxValueLastCol){
			maxValueLastCol = scoreArr[i*colsNum+colsNum-1];
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
		if(pathArr[i*colsNum+j]==0){ // from (i-1, j-1)
			queryAlnResult += local_aln->ctgseq[i-1];
			if(local_aln->ctgseq[i-1]!=local_aln->refseq[j-1]) {
				mismatchNum ++;
				midAlnResult += ' ';
			}else
				midAlnResult += '|';
			subjectAlnResult += local_aln->refseq[j-1];
			i --;
			j --;
		}else if(pathArr[i*colsNum+j]==1){ // from (i-1, j)
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

#if ALIGN_DEBUG
	cout << queryAlnResult << endl;
	cout << midAlnResult << endl;
	cout << subjectAlnResult << endl;
#endif

	local_aln->alignResultVec.push_back(queryAlnResult);
	local_aln->alignResultVec.push_back(midAlnResult);
	local_aln->alignResultVec.push_back(subjectAlnResult);

	free(scoreArr);
	free(pathArr);
}

// compute variant location
void varCand::computeVarLoc(localAln_t *local_aln){
	vector<string> var_loc_vec;
	string ctgseq_aln, midseq_aln, refseq_aln;
	int32_t queryPos, refPos, localRefPos, start_aln_idx_var, end_aln_idx_var;
	int32_t i, start_var_pos, end_var_pos, len, tmp_len;
	int32_t startRefPos_aln, endRefPos_aln, startLocalRefPos_aln, endLocalRefPos_aln, startQueryPos_aln, endQueryPos_aln;

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
		cerr << "start_aln_idx_var=" << start_aln_idx_var << ", end_aln_idx_var=" << end_aln_idx_var << endl;
		exit(1);
	}else {
#if ALIGN_DEBUG
			len = end_aln_idx_var - start_aln_idx_var + 1;
			cout << startRefPos_aln << "-" << endRefPos_aln << ", " << start_aln_idx_var << "-" << end_aln_idx_var  << endl;
			cout << "\t" << ctgseq_aln.substr(start_aln_idx_var-10, len+30) << endl;
			cout << "\t" << midseq_aln.substr(start_aln_idx_var-10, len+30) << endl;
			cout << "\t" << refseq_aln.substr(start_aln_idx_var-10, len+30) << endl;

#endif
	}

	local_aln->start_aln_idx_var = start_aln_idx_var;
	local_aln->end_aln_idx_var = end_aln_idx_var;

	local_aln->reg->startLocalRefPos = startLocalRefPos_aln;
	local_aln->reg->endLocalRefPos = endLocalRefPos_aln;
	local_aln->reg->startQueryPos = startQueryPos_aln;
	local_aln->reg->endQueryPos = endQueryPos_aln;
}

// confirm short variation
void varCand::confirmShortVar(localAln_t *local_aln){
	int32_t i, check_extend_size, start_check_idx, end_check_idx, startRefPos, endRefPos, mismatchNum_aln;
	string refseq_aln;
	vector<int32_t> numVec;

	check_extend_size = SHORT_VAR_ALN_CHECK_EXTEND_SIZE;
	start_check_idx = local_aln->start_aln_idx_var - check_extend_size;
	end_check_idx = local_aln->end_aln_idx_var + check_extend_size;
	if(start_check_idx<0)
		start_check_idx = 0;
	if(end_check_idx>=local_aln->overlapLen)
		end_check_idx = local_aln->overlapLen;

	// compute the refPos
	refseq_aln = local_aln->alignResultVec[2];
	startRefPos = local_aln->reg->startRefPos;
	for(i=local_aln->start_aln_idx_var; i>=start_check_idx; i--)
		if(refseq_aln[i]!='-')
			startRefPos --;
	endRefPos = local_aln->reg->endRefPos;
	for(i=local_aln->end_aln_idx_var; i<=end_check_idx; i++)
		if(refseq_aln[i]!='-')
			endRefPos ++;

	// compute the number of mismatched bases
	mismatchNum_aln = getMismatchNumAln(local_aln->alignResultVec, start_check_idx, end_check_idx);

	if(mismatchNum_aln>=MIN_MISMATCH_NUM_SHORT_VAR_CALL){
		numVec = computeDisagreeNumAndHighIndelBaseNum(local_aln->reg->chrname, startRefPos, endRefPos, inBamFile, fai);

		// confirm and update the information
		if(numVec[0]>=MIN_DISAGR_NUM_SHORT_VAR_CALL or numVec[1]>0) // confirmed
			computeVarType(local_aln->reg);  // compute variant type
	}
}

// get the number of mismatched bases
int32_t varCand::getMismatchNumAln(vector<string> &alignResultVec, int32_t start_check_idx, int32_t end_check_idx){
	int32_t i, mismatchNum, len;
	string check_str;

	mismatchNum = 0;
	check_str = alignResultVec[1].substr(start_check_idx, end_check_idx-start_check_idx+1);
	len = check_str.size();
	for(i=0; i<len; i++)
		if(check_str[i]!='|')
			mismatchNum ++;

	return mismatchNum;
}

// compute the number of disagreements
int32_t varCand::getDisagreeNum(Base *baseArray, int32_t arr_size){
	int32_t i, disagreeNum = 0;

	for(i=0; i<arr_size; i++)
		if(baseArray[i].isZeroCovBase() or baseArray[i].isDisagreeBase())
			disagreeNum ++;

	return disagreeNum;
}

// compute the number of high indel bases
int32_t varCand::computeHighIndelBaseNum(Base *baseArray, int32_t arr_size, float threshold){
	int32_t num = 0;
	for(int32_t i=0; i<arr_size; i++)
		if(baseArray[i].isHighIndelBase(threshold))
			num ++;
	return num;
}

// adjust the variant locations according to alignment
void varCand::adjustVarLoc(localAln_t *local_aln){
	int32_t start_aln_idx_var_new, end_aln_idx_var_new, tmp;

	// compute the adjusted locations
	start_aln_idx_var_new = getAdjustedStartAlnIdxVar(local_aln);
	end_aln_idx_var_new = getAdjustedEndAlnIdxVar(local_aln);

	if(start_aln_idx_var_new!=-1 and end_aln_idx_var_new!=-1)
		if(start_aln_idx_var_new>end_aln_idx_var_new){
			tmp = start_aln_idx_var_new;
			start_aln_idx_var_new = end_aln_idx_var_new;
			end_aln_idx_var_new = tmp;
			cout << "adjusted variant locations have been exchanged!" << endl;
		}

	// adjust the locations
	if(start_aln_idx_var_new!=-1 and start_aln_idx_var_new!=local_aln->start_aln_idx_var) // left location
		adjustLeftLocShortVar(local_aln, start_aln_idx_var_new);
	if(end_aln_idx_var_new!=-1 and end_aln_idx_var_new!=local_aln->end_aln_idx_var) // right location
		adjustRightLocShortVar(local_aln, end_aln_idx_var_new);
}

// compute the left mismatch location
int32_t varCand::getAdjustedStartAlnIdxVar(localAln_t *local_aln){
	int32_t i, minCheckIdx, maxCheckIdx, minMismatchIdx;
	string midseq_aln;

	midseq_aln = local_aln->alignResultVec[1];
	minCheckIdx = local_aln->start_aln_idx_var - SHORT_VAR_ALN_CHECK_EXTEND_SIZE;
	maxCheckIdx = local_aln->end_aln_idx_var + SHORT_VAR_ALN_CHECK_EXTEND_SIZE;
	if(minCheckIdx<0) minCheckIdx = 0;
	if(maxCheckIdx>=(int32_t)midseq_aln.size()) maxCheckIdx = midseq_aln.size() - 1;

	minMismatchIdx = -1;
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

	return minMismatchIdx;
}

// compute the right mismatch location
int32_t varCand::getAdjustedEndAlnIdxVar(localAln_t *local_aln){
	int32_t i, minCheckIdx, maxCheckIdx, maxMismatchIdx;
	string midseq_aln;

	midseq_aln = local_aln->alignResultVec[1];
	minCheckIdx = local_aln->start_aln_idx_var - SHORT_VAR_ALN_CHECK_EXTEND_SIZE;
	maxCheckIdx = local_aln->end_aln_idx_var + SHORT_VAR_ALN_CHECK_EXTEND_SIZE;
	if(minCheckIdx<0) minCheckIdx = 0;
	if(maxCheckIdx>=(int32_t)midseq_aln.size()) maxCheckIdx = midseq_aln.size() - 1;

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
		for(i=local_aln->start_aln_idx_var; i>=newStartAlnIdx; i--){
			if(ctgseq_aln[i]!='-') queryPos --;
			if(refseq_aln[i]!='-') { refPos --; localRefPos --; }
		}
		queryPos ++; refPos ++; localRefPos ++; // return to correct location
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

	if(newEndAlnIdx<local_aln->start_aln_idx_var){
		for(i=local_aln->end_aln_idx_var; i>=newEndAlnIdx; i--){
			if(ctgseq_aln[i]!='-') queryPos --;
			if(refseq_aln[i]!='-') { refPos --; localRefPos --; }
		}
		queryPos ++; refPos ++; localRefPos ++; // return to correct location
	}else{
		for(i=local_aln->end_aln_idx_var; i<=newEndAlnIdx; i++){
			if(ctgseq_aln[i]!='-') queryPos ++;
			if(refseq_aln[i]!='-') { refPos ++; localRefPos ++; }
		}
		queryPos --; refPos --; localRefPos --; // return to correct location
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
				if(tmp_len>=8){
					numVec = computeDisagreeNumAndHighIndelBaseNum(candRegVec[i-1]->chrname, startRefPos, endRefPos, inBamFile, fai);
					if(numVec[0]>=MIN_DISAGR_NUM_SHORT_VAR_CALL or numVec[1]>0) // found variant features, merge; otherwise, keep unchanged
						mergeFlag = true;
				}else mergeFlag = true;
			}

			if(mergeFlag){ // merge
				query_dist = (int32_t)candRegVec[i]->endQueryPos - (int32_t)candRegVec[i-1]->startQueryPos;
				ref_dist = (int32_t)candRegVec[i]->endRefPos - (int32_t)candRegVec[i-1]->startRefPos;

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
	int32_t i, j, idx, same_count, startRefPos, endRefPos, startQueryPos, endQueryPos, opID1, opID2;
	vector<int32_t> numVec1, numVec2;
	vector<reg_t*> updatedRegVec[foundRegVec.size()], regVec_tmp;
	vector< vector<reg_t*> > regVec;

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
				reg->startRefPos = reg_tmp->startRefPos - 1;
				reg->endRefPos = reg_tmp->endRefPos + 1;
				reg->startLocalRefPos = reg_tmp->startLocalRefPos - 1;
				reg->endLocalRefPos = reg_tmp->endLocalRefPos + 1;
				reg->startQueryPos = reg_tmp->startQueryPos - 1;
				reg->endQueryPos = reg_tmp->endQueryPos + 1;
				reg->query_id = reg_tmp->query_id;
				reg->blat_aln_id = reg_tmp->blat_aln_id;
				reg->aln_orient = reg_tmp->aln_orient;

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

			// compute the variant locations
			computeVarRegLoc(reg, reg_tmp);

			// check left marginal regions
			opID1 = -1;
			if(reg->startRefPos==reg_tmp->startRefPos){
				opID1 = 0;
				numVec1.insert(numVec1.begin(), 4, 0);
			}else if(reg->startRefPos<reg_tmp->startRefPos){
				startRefPos = reg->startRefPos;
				endRefPos = reg_tmp->startRefPos - 1;
				startQueryPos = reg->startQueryPos;
				endQueryPos = reg_tmp->startQueryPos - 1;
			}else{
				startRefPos = reg_tmp->startRefPos;
				endRefPos = reg->startRefPos - 1;
				startQueryPos = reg_tmp->startQueryPos;
				endQueryPos = reg->startQueryPos - 1;
			}
			if(opID1!=0){
				numVec1 = computeDisagreeNumAndHighIndelBaseNumAndMarginDist(reg->chrname, startRefPos, endRefPos, reg_tmp->query_id, startQueryPos, endQueryPos, reg->aln_orient, inBamFile, fai);
				if(numVec1[0]>=MIN_DISAGR_NUM_SHORT_VAR_CALL or numVec1[1]>0){ // found variant features
					if(numVec1[3]>=8 and numVec1[3]<endRefPos-startRefPos+1)
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
				numVec2.insert(numVec2.begin(), 4, 0);
			}else if(reg_tmp->endRefPos<reg->endRefPos){
				startRefPos = reg_tmp->endRefPos + 1;
				endRefPos = reg->endRefPos;
				startQueryPos = reg_tmp->endQueryPos + 1;
				endQueryPos = reg->endQueryPos;
			}else{
				startRefPos = reg->endRefPos + 1;
				endRefPos = reg_tmp->endRefPos;
				startQueryPos = reg->endQueryPos + 1;
				endQueryPos = reg_tmp->endQueryPos;
			}

			if(opID2!=0){
				numVec2 = computeDisagreeNumAndHighIndelBaseNumAndMarginDist(reg->chrname, startRefPos, endRefPos, reg_tmp->query_id, startQueryPos, endQueryPos, reg->aln_orient, inBamFile, fai);
				if(numVec2[0]>=MIN_DISAGR_NUM_SHORT_VAR_CALL or numVec2[1]>0){ // found variant features
					if(numVec2[2]>=8 and numVec2[2]<endRefPos-startRefPos+1)
						opID2 = 2;
					else
						opID2 = 1;
				}else
					opID2 = 1;
			}

			// update their information at the same time
			if(opID1!=2 and opID2!=2){
				if(opID1==0){
					reg->startRefPos = reg_tmp->startRefPos;
					reg->startLocalRefPos = reg_tmp->startLocalRefPos;
					reg->startQueryPos = reg_tmp->startQueryPos;
				}else if(opID1==1){ // update information
					if(reg->startRefPos<reg_tmp->startRefPos){
						reg->startRefPos += numVec1[2];
						reg->startLocalRefPos += numVec1[2];
						reg->startQueryPos += numVec1[2];
					}else{
						reg->startRefPos = reg_tmp->startRefPos + numVec1[2];
						reg->startLocalRefPos = reg_tmp->startLocalRefPos + numVec1[2];
						reg->startQueryPos = reg_tmp->startQueryPos + numVec1[2];
					}
				}

				if(opID2==0){
					reg->endRefPos = reg_tmp->endRefPos;
					reg->endLocalRefPos = reg_tmp->endLocalRefPos;
					reg->endQueryPos = reg_tmp->endQueryPos;
				}else if(opID2==1){
					if(reg_tmp->endRefPos<reg->endRefPos){
						reg->endRefPos -= numVec2[3];
						reg->endLocalRefPos -= numVec2[3];
						reg->endQueryPos -= numVec2[3];
					}else{
						reg->endRefPos = reg_tmp->endRefPos - numVec2[3];
						reg->endLocalRefPos = reg_tmp->endLocalRefPos - numVec2[3];
						reg->endQueryPos = reg_tmp->endQueryPos - numVec2[3];
					}
				}
				reg->query_id = reg_tmp->query_id;
				updatedRegVec[i].push_back(reg);

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

				// update information according to opID1==1
				if(reg_tmp->startRefPos<reg->startRefPos){
					reg_new_3->startRefPos = reg_tmp->startRefPos + numVec1[2];
					reg_new_3->startLocalRefPos = reg_tmp->startLocalRefPos + numVec1[2];
					reg_new_3->startQueryPos = reg_tmp->startQueryPos + numVec1[2];
				}else{
					reg_new_3->startRefPos = reg->startRefPos + numVec1[2];
					reg_new_3->startLocalRefPos = reg->startLocalRefPos + numVec1[2];
					reg_new_3->startQueryPos = reg->startQueryPos + numVec1[2];
				}

				updatedRegVec[i].push_back(reg_new_3);
				updatedRegVec[i].push_back(reg_new_4);

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

				// update information according to opID2==1
				if(reg_tmp->endRefPos<reg->endRefPos){
					reg_new_2->endRefPos = reg->endRefPos - numVec2[3];
					reg_new_2->endLocalRefPos = reg->endLocalRefPos - numVec2[3];
					reg_new_2->endQueryPos = reg->endQueryPos - numVec2[3];
				}else{
					reg_new_2->endRefPos = reg_tmp->endRefPos - numVec2[3];
					reg_new_2->endLocalRefPos = reg_tmp->endLocalRefPos - numVec2[3];
					reg_new_2->endQueryPos = reg_tmp->endQueryPos - numVec2[3];
				}

				updatedRegVec[i].push_back(reg_new_1);
				updatedRegVec[i].push_back(reg_new_2);

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

				updatedRegVec[i].push_back(reg_new_1);
				updatedRegVec[i].push_back(reg_new_2);
				updatedRegVec[i].push_back(reg_new_4);
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
			computeVarType(reg);  // compute variant type
			reg->short_sv_flag = false;
		}
	}

	return regVec;
}

// compute variant type
void varCand::computeVarType(reg_t *reg){
	int32_t query_dist, ref_dist, var_type, decrease_size;
	query_dist = (int32_t)reg->endQueryPos - (int32_t)reg->startQueryPos;
	ref_dist = (int32_t)reg->endRefPos - (int32_t)reg->startRefPos;

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
	reg->sv_len = query_dist - ref_dist;
	reg->call_success_status = true;
}

// update varVec vector
void varCand::updateVarVec(vector<vector<reg_t*>> &regVec, vector<reg_t*> &foundRegVec, vector<reg_t*> &varVec){
	int32_t i, j, reg_idx;
	reg_t *reg;
	vector<reg_t*> reg_vec_item;

	if(regVec.size()>0){
		for(i=0; i<(int32_t)foundRegVec.size(); i++){
			reg = foundRegVec[i];
			reg_vec_item = regVec[i];
			if(reg_vec_item.size()>1 or reg_vec_item[0]!=reg){
				reg_idx = getVectorIdx(reg, varVec);
				if(reg_idx!=-1){
					for(j=0; j<(int32_t)reg_vec_item.size(); j++)
						varVec.insert(varVec.begin()+reg_idx+j, reg_vec_item[j]);
					varVec.erase(varVec.begin()+reg_idx+reg_vec_item.size());
					varVec.shrink_to_fit();
					delete reg;
				}else{
					cerr << "line=" << __LINE__ << ", error reg_idx=-1!" << endl;
					exit(1);
				}
			}
		}
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
	disagreeNum = getDisagreeNum(baseArray, endRefPos-startRefPos+1);

	// further to do: compute the number of high indel bases
	highIndelBaseNum = computeHighIndelBaseNum(baseArray, endRefPos-startRefPos+1, MIN_HIGH_INDEL_BASE_RATIO);

	// compute margins of variations
	distVec = computeVarMargins(baseArray, endRefPos-startRefPos+1, MIN_HIGH_INDEL_BASE_RATIO);

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
	local_aln = new localAln_t;
	local_aln->reg = reg;
	local_aln->cand_reg = cand_reg;
	local_aln->startRefPos = local_aln->endRefPos = -1;
	local_aln->startLocalRefPos = local_aln->startQueryPos = local_aln->endLocalRefPos = local_aln->endQueryPos = -1;
	local_aln->queryLeftShiftLen = local_aln->queryRightShiftLen = local_aln->localRefLeftShiftLen = local_aln->localRefRightShiftLen = -1;
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

	delete local_aln;
}

// compute margins of variations
vector<int32_t> varCand::computeVarMargins(Base *baseArray, int32_t arr_size, float threshold){
	vector<int32_t> distVec;
	int32_t i, leftDist, rightDist, discorNum;

	// compute left distance
	leftDist = 0;
	discorNum = 0;
	for(i=0; i<arr_size; i++){
		if(!baseArray[i].isDisagreeBase() and !baseArray[i].isHighIndelBase(threshold))
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
	for(i=arr_size-1; i>=0; i--){
		if(!baseArray[i].isDisagreeBase() and !baseArray[i].isHighIndelBase(threshold))
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
	int32_t i, leftDist, rightDist, tmp_len, ref_len, ctg_len, misNum;
	string refseq, ctgseq, reg_str;

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
		if(refseq[i]==ctgseq[i])
			leftDist ++;
		else{
			misNum ++;
			if(misNum<MARGIN_MISMATCH_NUM_IGNORE_THRES and leftDist>MIN_MARGIN_LEN_MISMATCH_IGNORE) leftDist ++;
			else break;
		}
	}
	rightDist = 0;
	misNum = 0;
	for(i=0; i<tmp_len; i++){
		if(refseq[ref_len-1-i]==ctgseq[ctg_len-1-i])
			rightDist ++;
		else{
			misNum ++;
			if(misNum<MARGIN_MISMATCH_NUM_IGNORE_THRES and rightDist>MIN_MARGIN_LEN_MISMATCH_IGNORE) rightDist ++;
			else break;
		}
	}

//	if(leftDist>left_dist)
//		leftDist = left_dist;
//	if(rightDist>right_dist)
//		rightDist = right_dist;
	if(leftDist<left_dist)
		leftDist = left_dist;
	if(rightDist<right_dist)
		rightDist = right_dist;
	distVec.push_back(leftDist);
	distVec.push_back(rightDist);
	return distVec;
}

// compute local locations
void varCand::computeLocalLocsAln(localAln_t *local_aln){
	reg_t *reg = local_aln->reg, *cand_reg = local_aln->cand_reg;
	int32_t dist0, dist1, dist2;
	size_t querylen;
	blat_aln_t *blat_aln;

	if(reg->endRefPos-reg->startRefPos>cand_reg->endRefPos-cand_reg->startRefPos)
		dist0 = ((reg->endRefPos - reg->startRefPos) / VAR_ALN_EXTEND_SIZE + 1) * VAR_ALN_EXTEND_SIZE;
	else
		dist0 = ((cand_reg->endRefPos - cand_reg->startRefPos) / VAR_ALN_EXTEND_SIZE + 1) * VAR_ALN_EXTEND_SIZE;

	if(reg->startRefPos<cand_reg->startRefPos)
		dist1 = dist0 + cand_reg->startRefPos - reg->startRefPos;
	else
		dist1 = dist0;

	if(cand_reg->startLocalRefPos<cand_reg->startQueryPos){
		if((int32_t)cand_reg->startLocalRefPos<dist1) dist1 = cand_reg->startLocalRefPos - 1;
	}else{
		if((int32_t)cand_reg->startQueryPos<dist1) dist1 = cand_reg->startQueryPos - 1;
	}

	if(cand_reg->endRefPos<reg->endRefPos)
		dist2 = dist0;
	else
		dist2 = dist0 + cand_reg->endRefPos - reg->endRefPos;
	blat_aln = blat_aln_vec.at(cand_reg->blat_aln_id);
	querylen = blat_aln->query_len;
	if(cand_reg->endQueryPos+dist2>querylen) dist2 = querylen - cand_reg->endRefPos;

	local_aln->startRefPos = cand_reg->startRefPos - dist1;
	local_aln->startLocalRefPos = cand_reg->startLocalRefPos - dist1;
	local_aln->startQueryPos = cand_reg->startQueryPos - dist1;

	local_aln->endRefPos = cand_reg->endRefPos + dist2;
	local_aln->endLocalRefPos = cand_reg->endLocalRefPos + dist2;
	local_aln->endQueryPos = cand_reg->endQueryPos + dist2;
}

// distinguish short duplications and inversions from indels
void varCand::distinguishShortDupInvFromIndels(){
	size_t i;
	reg_t *reg;
	string refseq, queryseq;
	vector<size_t> left_shift_size_vec, right_shift_size_vec;
	size_t leftRefShiftSize, leftQueryShiftSize, rightRefShiftSize, rightQueryShiftSize;
	size_t leftVarRefPos, leftVarLocalRefPos, leftVarQueryPos, rightVarRefPos, rightVarLocalRefPos, rightVarQueryPos;
	int32_t dup_num_int, ref_dist, query_dist, sv_len_tmp, blat_aln_idx;
	double dup_num_tmp, sim_ratio_inv;
	localAln_t *local_aln;

	for(i=0; i<varVec.size(); i++){
		reg = varVec.at(i);

		//cout << "startRefPos=" << reg->startRefPos << ", endRefPos=" << reg->endRefPos << ", startLocalRefPos=" << reg->startLocalRefPos << ", endLocalRefPos=" << reg->endLocalRefPos << ", startQueryPos=" << reg->startQueryPos << ", endQueryPos=" << reg->endQueryPos << ", aln_orient=" << reg->aln_orient << ", var_type=" << reg->var_type << ", sv_len=" << reg->sv_len << ", short_sv_flag=" << reg->short_sv_flag << endl;

		if(reg->short_sv_flag==false){
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
					reg->startRefPos = leftVarRefPos;
					reg->endRefPos = rightVarRefPos;
					reg->startLocalRefPos = leftVarLocalRefPos;
					reg->endLocalRefPos = rightVarLocalRefPos;
					reg->startQueryPos = leftVarQueryPos;
					reg->endQueryPos = rightVarQueryPos;
					reg->call_success_status = true;
					reg->short_sv_flag = true;
					reg->dup_num = dup_num_int;
					reg->sv_len = query_dist - ref_dist;

					cout << "----------> Short DUP: " << reg->chrname << ":" << reg->startRefPos << "-" << reg->endRefPos << ", startLocalRefPos=" << reg->startLocalRefPos << ", endLocalRefPos=" << reg->endLocalRefPos << ", startQueryPos=" << reg->startQueryPos << ", endQueryPos=" << reg->endQueryPos << ", sv_len=" << reg->sv_len << endl;
				}

			}else if(reg->sv_len<0.2*(reg->endRefPos-reg->startRefPos+1)){ // inversions? translocations?
				ref_dist = reg->endRefPos - reg->startRefPos + 1;
				query_dist = reg->endQueryPos - reg->startQueryPos + 1;
				if(ref_dist>=MIN_SHORT_INV_SIZE and query_dist>=MIN_SHORT_INV_SIZE){ // size limitation
					// compute the corresponding blat_aln item
					blat_aln_idx = getInvBlatAlnItemIdx(reg, blat_aln_vec, reg->blat_aln_id, NULL, NULL);
					if(blat_aln_idx!=-1){
						local_aln = new localAln_t();
						local_aln->reg = NULL;
						local_aln->aln_seg = NULL;
						local_aln->startRefPos = local_aln->endRefPos = -1;
						local_aln->startLocalRefPos = local_aln->startQueryPos = local_aln->endLocalRefPos = local_aln->endQueryPos = -1;
						local_aln->queryLeftShiftLen = local_aln->queryRightShiftLen = local_aln->localRefLeftShiftLen = local_aln->localRefRightShiftLen = -1;

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
						if(sim_ratio_inv>=SIMILARITY_THRES_INV){
							reg->var_type = VAR_INV;
							reg->call_success_status = true;
							reg->short_sv_flag = true;

							cout << "==========> Short INV: " << reg->chrname << ":" << reg->startRefPos << "-" << reg->endRefPos << ", startLocalRefPos=" << reg->startLocalRefPos << ", endLocalRefPos=" << reg->endLocalRefPos << ", startQueryPos=" << reg->startQueryPos << ", endQueryPos=" << reg->endQueryPos << ", sv_len=" << reg->sv_len << ", sim_ratio_inv=" << sim_ratio_inv << endl;
						}

						cout << "sim_ratio_inv=" << sim_ratio_inv << endl;

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
	size_t i, leftShiftSize, rightShiftSize, misNum, misNum_total, aln_size_total;

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
				if(reg->call_success_status)
					reg->refseq = refseqloader.getFastaSeqByPos(0, reg->startLocalRefPos, reg->endLocalRefPos, ALN_PLUS_ORIENT);
			}
			for(i=0; i<newVarVec.size(); i++){
				reg = newVarVec[i];
				if(reg->call_success_status)
					reg->refseq = refseqloader.getFastaSeqByPos(0, reg->startLocalRefPos, reg->endLocalRefPos, ALN_PLUS_ORIENT);
			}

			// get the query sequence
			FastaSeqLoader ctgseqloader(ctgfilename);
			for(i=0; i<varVec.size(); i++){
				reg = varVec[i];
				if(reg->call_success_status)
					reg->altseq = ctgseqloader.getFastaSeqByPos(reg->query_id, reg->startQueryPos, reg->endQueryPos, reg->aln_orient);
			}
			for(i=0; i<newVarVec.size(); i++){
				reg = newVarVec[i];
				if(reg->call_success_status)
					reg->altseq = ctgseqloader.getFastaSeqByPos(reg->query_id, reg->startQueryPos, reg->endQueryPos, reg->aln_orient);
			}
		}else{ // SVs
			// get the reference and query sequence
			FastaSeqLoader refseqloader(refseqfilename);
			FastaSeqLoader ctgseqloader(ctgfilename);
			if(call_success){
				if(clip_reg->var_type==VAR_DUP or clip_reg->var_type==VAR_INV){
					clip_reg->refseq = refseqloader.getFastaSeqByPos(0, clip_reg->startLocalRefPos, clip_reg->endLocalRefPos, ALN_PLUS_ORIENT);
					clip_reg->altseq = ctgseqloader.getFastaSeqByPos(clip_reg->query_id, clip_reg->startQueryPos, clip_reg->endQueryPos, clip_reg->aln_orient);
				}
			}
		}
	}

/*
	// print the seqs
	if(align_success){
		for(i=0; i<varVec.size(); i++){
			reg = varVec[i];
			if(reg->call_success_status){
				cout << reg->chrname << ":" << reg->startRefPos << "-" << reg->endRefPos << endl;
				cout << "\tseqs: " << reg->refseq << "\t" << reg->altseq << endl;
			}else
				cout << "cannot be called" << endl;
		}
		for(i=0; i<newVarVec.size(); i++){
			reg = newVarVec[i];
			if(reg->call_success_status){
				cout << "New SV: " << reg->chrname << ":" << reg->startRefPos << "-" << reg->endRefPos << endl;
				cout << "\tseqs: " << reg->refseq << "\t" << reg->altseq << endl;
			}else
				cout << "cannot be called" << endl;
		}
	}
*/
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
		startClipPos = varVec.at(0)->startRefPos;
		endClipPos = varVec.at(varVec.size()-1)->endRefPos;
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
				if((ref_dist<MIN_SHORT_DUP_SIZE and sv_len_tmp>=MIN_SHORT_DUP_SIZE) and (seg1->ref_end>=(size_t)startClipPos and seg1->ref_end<=(size_t)endClipPos) and (seg2->ref_start>=(size_t)startClipPos and seg2->ref_start<=(size_t)endClipPos)){
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
	int32_t segIdx_start, segIdx_end, startClipPos, endClipPos, ref_dist, query_dist, blat_aln_idx;
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
		segIdx_start = -1;
		blat_aln = blat_aln_vec.at(i);
		if(blat_aln->valid_aln==true){
			for(j=1; j<blat_aln->aln_segs.size(); j++){
				seg1 = blat_aln->aln_segs[j-1];
				seg2 = blat_aln->aln_segs[j];
				if((seg1->ref_end>=(size_t)startClipPos and seg1->ref_end<=(size_t)endClipPos) and (seg2->ref_start>=(size_t)startClipPos and seg2->ref_start<=(size_t)endClipPos)){
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
				blat_aln_idx = getInvBlatAlnItemIdx(NULL, blat_aln_vec, i, seg1, seg2);
				if(blat_aln_idx!=-1){
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

					ref_dist = clip_reg->endLocalRefPos - clip_reg->startLocalRefPos + 1;
					query_dist = clip_reg->endQueryPos - clip_reg->startQueryPos + 1;
					clip_reg->sv_len = query_dist - ref_dist + 1;

					// check the inverted sequences
					local_aln = new localAln_t();
					local_aln->reg = NULL;
					local_aln->cand_reg = NULL;
					local_aln->aln_seg = NULL;
					local_aln->startRefPos = local_aln->endRefPos = -1;
					local_aln->startLocalRefPos = local_aln->startQueryPos = local_aln->endLocalRefPos = local_aln->endQueryPos = -1;
					local_aln->queryLeftShiftLen = local_aln->queryRightShiftLen = local_aln->localRefLeftShiftLen = local_aln->localRefRightShiftLen = -1;

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
						clip_reg->var_type = VAR_INV;
						clip_reg->call_success_status = true;
						clip_reg->short_sv_flag = false;
						call_success = true;

						cout << "----------> INV: " << clip_reg->chrname << ":" << clip_reg->startRefPos << "-" << clip_reg->endRefPos << ", startLocalRefPos=" << clip_reg->startLocalRefPos << ", endLocalRefPos=" << clip_reg->endLocalRefPos << ", startQueryPos=" << clip_reg->startQueryPos << ", endQueryPos=" << clip_reg->endQueryPos << ", sv_len=" << clip_reg->sv_len << ", sim_ratio_inv=" << sim_ratio_inv << endl;
					}

					delete local_aln;

					break;
				}
			}
		}
	}
}

int32_t varCand::getInvBlatAlnItemIdx(reg_t *reg, vector<blat_aln_t*> &blat_aln_vec, size_t given_idx, aln_seg_t *seg1, aln_seg_t *seg2){
	int32_t blat_aln_idx = -1;
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
						break;
					}
				}
			}
			if(blat_aln_idx!=-1) break;
		}
	}

	return blat_aln_idx;
}

// compute left and right clipping positions
reg_t* varCand::computeClipPos(blat_aln_t *blat_aln, aln_seg_t *seg1, aln_seg_t *seg2, string &refseq, string &queryseq, size_t var_type){
	size_t leftRefShiftSize, leftQueryShiftSize, rightRefShiftSize, rightQueryShiftSize;
	size_t leftClipRefPos, leftClipLocalRefPos, leftClipQueryPos, rightClipRefPos, rightClipLocalRefPos, rightClipQueryPos;
	int32_t dup_num_int, ref_dist, query_dist, minLeftClipPos, maxLeftClipPos, minRightClipPos, maxRightClipPos;
	double dup_num_tmp;
	vector<size_t> left_shift_size_vec, right_shift_size_vec;
	vector<clipPos_t> clip_pos_vec;
	clipPos_t clip_pos;
	reg_t *clip_reg_ret = NULL;
	bool call_success_flag = true;
	vector<size_t> pos_vec;

	if(var_type==VAR_DUP){ // duplication
		left_shift_size_vec = computeLeftShiftSizeDup(varVec.at(0), seg1, seg2, refseq, queryseq);
		right_shift_size_vec = computeRightShiftSizeDup(varVec.at(varVec.size()-1), seg1, seg2, refseq, queryseq);

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
		if((int32_t)leftClipRefPos>minLeftClipPos and (int32_t)leftClipRefPos<maxLeftClipPos){
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

		if((int32_t)rightClipRefPos>minRightClipPos and (int32_t)rightClipRefPos<maxRightClipPos){
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
			clip_reg_ret->startRefPos =  clip_pos_vec.at(0).clipRefPos;
			clip_reg_ret->endRefPos =  clip_pos_vec.at(1).clipRefPos;
			clip_reg_ret->startLocalRefPos =  clip_pos_vec.at(0).clipLocalRefPos;
			clip_reg_ret->endLocalRefPos =  clip_pos_vec.at(1).clipLocalRefPos;
			clip_reg_ret->startQueryPos =  clip_pos_vec.at(0).clipQueryPos;
			clip_reg_ret->endQueryPos =  clip_pos_vec.at(1).clipQueryPos;
			clip_reg_ret->aln_orient = clip_pos_vec.at(0).aln_orient;
			clip_reg_ret->var_type = var_type;
			clip_reg_ret->query_id = blat_aln->query_id;
			clip_reg_ret->blat_aln_id = blat_aln->blat_aln_id;
			clip_reg_ret->call_success_status = true;
			clip_reg_ret->short_sv_flag = false;

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
	}else if(var_type==VAR_INV){ // inversion

	}else if(var_type==VAR_TRA){ // translocation

	}

	return clip_reg_ret;
}

vector<size_t> varCand::computeLeftShiftSizeDup(reg_t *reg, aln_seg_t *seg1, aln_seg_t *seg2, string &refseq, string &queryseq){
	size_t startCheckLocalRefPos, startCheckRefPos, startCheckQueryPos, localRefPos, refPos, queryPos, misNum, refShiftSize, queryShiftSize, subseq_len;
	size_t i, localRefPos_start, queryPos_start, tmp_len;
	int32_t minCheckRefPos;
	localAln_t *local_aln;
	string refseq_aln, midseq_aln, queryseq_aln;
	vector<size_t> shift_size_vec;

	minCheckRefPos = reg->startRefPos - 200;
	if(minCheckRefPos<1) minCheckRefPos = 1;

	subseq_len = 50;
	if(seg1 and seg2){
		startCheckLocalRefPos = seg1->subject_end;
		startCheckRefPos = seg1->ref_end;
		startCheckQueryPos = seg2->query_start - 1;
	}else{
		startCheckLocalRefPos = reg->startLocalRefPos;
		startCheckRefPos = reg->startRefPos;
		startCheckQueryPos = reg->endQueryPos - 1;
	}

	misNum = 0;
	localRefPos = startCheckLocalRefPos;
	refPos = startCheckRefPos;
	queryPos = startCheckQueryPos;
	while(queryPos>=1 and localRefPos>=1){
		if(queryseq.at(queryPos-1)==refseq.at(localRefPos-1)){
			queryPos --; localRefPos --; refPos --;
		}else{
			misNum = 0;


			if(localRefPos<subseq_len) subseq_len = localRefPos;
			if(queryPos<subseq_len) subseq_len = queryPos;

			localRefPos_start = localRefPos - subseq_len + 1;
			queryPos_start = queryPos - subseq_len + 1;



			// compute local locations
			local_aln = new localAln_t;
			local_aln->reg = NULL;
			local_aln->aln_seg = NULL;
			local_aln->startRefPos = local_aln->endRefPos = -1;
			local_aln->startLocalRefPos = local_aln->startQueryPos = local_aln->endLocalRefPos = local_aln->endQueryPos = -1;
			local_aln->queryLeftShiftLen = local_aln->queryRightShiftLen = local_aln->localRefLeftShiftLen = local_aln->localRefRightShiftLen = -1;

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
				queryPos -=  subseq_len - local_aln->queryLeftShiftLen;
				localRefPos -=  subseq_len - local_aln->localRefLeftShiftLen;
				refPos -= subseq_len - local_aln->localRefLeftShiftLen;
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
					queryPos -= tmp_len + local_aln->queryRightShiftLen;
					localRefPos -= tmp_len + local_aln->localRefRightShiftLen;
					refPos -= tmp_len + local_aln->localRefRightShiftLen;
				}

				delete local_aln;
				break;
			}
		}

		if((int32_t)refPos<minCheckRefPos) break;
	}

	refShiftSize = startCheckLocalRefPos - localRefPos;
	queryShiftSize = startCheckQueryPos - queryPos;
	shift_size_vec.push_back(refShiftSize);
	shift_size_vec.push_back(queryShiftSize);

	//cout << "refShiftSize=" << refShiftSize << ", queryShiftSize=" << queryShiftSize << ", left localRefPos=" << localRefPos << ", left refPos=" << refPos << ", left queryPos=" << queryPos << endl;

	return shift_size_vec;
}

vector<size_t> varCand::computeRightShiftSizeDup(reg_t *reg, aln_seg_t *seg1, aln_seg_t *seg2, string &refseq, string &queryseq){
	size_t startCheckLocalRefPos, startCheckRefPos, startCheckQueryPos, localRefPos, refPos, queryPos, misNum, refShiftSize, queryShiftSize, subseq_len;
	size_t i, localRefPos_start, queryPos_start, tmp_len;
	int32_t maxCheckRefPos;
	localAln_t *local_aln;
	string midseq_aln;
	vector<size_t> shift_size_vec;

	maxCheckRefPos = reg->endRefPos + 200;

	subseq_len = 50;
	if(seg1 and seg2){
		startCheckLocalRefPos = seg2->subject_start;
		startCheckRefPos = seg2->ref_start;
		startCheckQueryPos = seg1->query_end + 1;
	}else{
		startCheckLocalRefPos = reg->endLocalRefPos;
		startCheckRefPos = reg->endRefPos;
		startCheckQueryPos = reg->startQueryPos + 1;
	}

	misNum = 0;
	localRefPos = startCheckLocalRefPos;
	refPos = startCheckRefPos;
	queryPos = startCheckQueryPos;
	while(queryPos<=queryseq.size() and localRefPos<=refseq.size()){
		if(queryseq.at(queryPos-1)==refseq.at(localRefPos-1)){
			queryPos ++; localRefPos ++; refPos ++;
		}else{
			misNum = 0;

			localRefPos_start = localRefPos;
			queryPos_start = queryPos;

			tmp_len = refseq.size() - localRefPos;
			if(tmp_len<subseq_len) subseq_len = tmp_len;
			tmp_len = queryseq.size() - queryPos;
			if(tmp_len<subseq_len) subseq_len = tmp_len;

			// compute local locations
			local_aln = new localAln_t;
			local_aln->reg = NULL;
			local_aln->aln_seg = NULL;
			local_aln->startRefPos = local_aln->endRefPos = -1;
			local_aln->startLocalRefPos = local_aln->startQueryPos = local_aln->endLocalRefPos = local_aln->endQueryPos = -1;
			local_aln->queryLeftShiftLen = local_aln->queryRightShiftLen = local_aln->localRefLeftShiftLen = local_aln->localRefRightShiftLen = -1;

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
				queryPos +=  subseq_len - local_aln->queryRightShiftLen;
				localRefPos +=  subseq_len - local_aln->localRefRightShiftLen;
				refPos += subseq_len - local_aln->localRefRightShiftLen;
				delete local_aln;
			}else{
				if(local_aln->localRefLeftShiftLen<=2 and local_aln->queryLeftShiftLen<=2){
					// compute the number of shift bases
					tmp_len = 0;
					midseq_aln = local_aln->alignResultVec.at(1);
					for(i=0; i<midseq_aln.size(); i++){
						if(midseq_aln.at(i)=='|') tmp_len ++; // match
						else break; // mismatch
					}
					queryPos += tmp_len + local_aln->queryLeftShiftLen;
					localRefPos += tmp_len + local_aln->localRefLeftShiftLen;
					refPos += tmp_len + local_aln->localRefLeftShiftLen;
				}

				delete local_aln;
				break;
			}
		}

		if((int32_t)refPos>maxCheckRefPos) break;
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
vector<size_t> varCand::computeQueryClipPosDup(blat_aln_t *blat_aln, size_t clipRefPos, string &refseq, string &queryseq){
	int32_t localRefPos_start, refPos_start, queryPos_start;
	size_t i, subseq_len, dist, queryPos, refPos, localRefPos;
	localAln_t *local_aln;
	aln_seg_t *aln_seg;
	string query_aln_seq, mid_aln_seq, ref_aln_seq;
	vector<size_t> pos_vec;

	// get aln_seg
	aln_seg = NULL;
	for(i=0; i<blat_aln->aln_segs.size(); i++){
		aln_seg = blat_aln->aln_segs[i];
		if(aln_seg->ref_start<=clipRefPos and aln_seg->ref_end>=clipRefPos) break;
	}

	if(aln_seg){
		dist = clipRefPos - aln_seg->ref_start;
		refPos_start = clipRefPos - 200;
		queryPos_start = aln_seg->query_start + dist - 200;
		localRefPos_start = aln_seg->subject_start + dist - 200;
		if(queryPos_start<1) queryPos_start = 1;
		if(localRefPos_start<1) localRefPos_start = 1;

		// compute local locations
		local_aln = new localAln_t;
		local_aln->reg = NULL;
		local_aln->aln_seg = NULL;
		local_aln->startRefPos = local_aln->endRefPos = -1;
		local_aln->startLocalRefPos = local_aln->startQueryPos = local_aln->endLocalRefPos = local_aln->endQueryPos = -1;
		local_aln->queryLeftShiftLen = local_aln->queryRightShiftLen = local_aln->localRefLeftShiftLen = local_aln->localRefRightShiftLen = -1;

		subseq_len = 400;
		local_aln->refseq = refseq.substr(localRefPos_start-1, subseq_len);
		local_aln->ctgseq = queryseq.substr(queryPos_start-1, subseq_len);

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
				if(query_aln_seq.at(i)!='-'){
					refPos ++; localRefPos ++;
				}else{
					queryPos ++;
				}
			}
			if(refPos==clipRefPos) break;
		}
		pos_vec.push_back(refPos);
		pos_vec.push_back(localRefPos);
		pos_vec.push_back(queryPos);

		delete local_aln;
	}

	return pos_vec;
}

void varCand::printSV(){
	size_t j;
	reg_t *reg;
	string line, sv_type;

	for(j=0; j<varVec.size(); j++){
		reg = varVec[j];
		switch(reg->var_type){
			case VAR_UNC: sv_type = "UNCERTAIN"; break;
			case VAR_INS: sv_type = "INS"; break;
			case VAR_DEL: sv_type = "DEL"; break;
			case VAR_DUP: sv_type = "DUP"; break;
			case VAR_INV: sv_type = "INV"; break;
			case VAR_TRA: sv_type = "TRA"; break;
			default: sv_type = "MIX"; break;
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
