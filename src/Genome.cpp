#include <sys/stat.h>
#include "Genome.h"

// Constructor with parameters
Genome::Genome(Paras *paras){
	this->paras = paras;
	init();
}

//Destructor
Genome::~Genome(){
	destroyChromeVector();
	fai_destroy(fai);
	bam_hdr_destroy(header);
}

// initialization
void Genome::init(){
	Chrome *chr;
	string chrname_tmp;

	out_filename_detect_snv = out_dir_detect + "/" + "genome_SNV_candidate";
	out_filename_detect_indel = out_dir_detect + "/" + "genome_INDEL_candidate";
	out_filename_detect_clipReg = out_dir_detect + "/" + "genome_clipReg_candidate";
	out_filename_call_snv = out_dir_call + "/" + "genome_SNV";
	out_filename_call_indel = out_dir_call + "/" + "genome_INDEL";
	out_filename_call_clipReg = out_dir_call + "/" + "genome_clipReg";
	out_filename_call_tra = out_dir_call + "/" + "genome_TRA";

	// load the fai
	fai = fai_load(paras->refFile.c_str());
	if ( !fai ) {
		cerr << __func__ << ": could not load fai index of " << paras->refFile << endl;;
		exit(1);
	}

	// load the sam/bam header
	header = loadSamHeader(paras->inBamFile);

	// allocate each genome
	for(int i=0; i<header->n_targets; i++){
		chrname_tmp = header->target_name[i];
		chr = allocateChrome(chrname_tmp, header->target_len[i], fai);
		chr->setOutputDir(out_dir_detect, out_dir_assemble, out_dir_call);
		chromeVector.push_back(chr);
	}
	chromeVector.shrink_to_fit();
}

// allocate the Chrome node
Chrome* Genome::allocateChrome(string& chrname, int chrlen, faidx_t *fai){
	Chrome *chr_tmp;
	chr_tmp = new Chrome(chrname, chrlen, fai, paras);
	if(!chr_tmp){ cerr << "Genome: cannot allocate memory" << endl; exit(1); }
	return chr_tmp;
}

// get genome size
int Genome::getGenomeSize()
{
    genomeSize = 0;
    for(int i=0; i<header->n_targets; i++) genomeSize += header->target_len[i];
	return 0;
}

// release each chrome in chromeVector
void Genome::destroyChromeVector(){
	vector<Chrome*>::iterator chr;
	for(chr=chromeVector.begin(); chr!=chromeVector.end(); chr++)
		delete (*chr);   // release each chrome
	vector<Chrome*>().swap(chromeVector);
}

// generate the genome blocks
int Genome::generateGenomeBlocks(){
	vector<Chrome*>::iterator chr;
	for(chr=chromeVector.begin(); chr!=chromeVector.end(); chr++)
		(*chr)->generateChrBlocks();
	return 0;
}

// save the genome blocks to file
void Genome::saveGenomeBlocksToFile(){
	vector<Chrome*>::iterator chr;
	for(chr=chromeVector.begin(); chr!=chromeVector.end(); chr++)
		(*chr)->saveChrBlocksToFile();
}

// estimate the insertion/deletion/clipping parameters
void Genome::estimateSVSizeNum(){
	// initialize the data
	paras->initEst();

	// fill the size data using each chromosome
	vector<Chrome*>::iterator chr;
	for(chr=chromeVector.begin(); chr!=chromeVector.end(); chr++)
		(*chr)->chrFillDataEst(SIZE_EST_OP);
	// size estimate
	paras->estimate(SIZE_EST_OP);

	// fill the num data using each chromosome
	for(chr=chromeVector.begin(); chr!=chromeVector.end(); chr++)
		(*chr)->chrFillDataEst(NUM_EST_OP);
	// size estimate
	paras->estimate(NUM_EST_OP);
}

// detect variations for genome
int Genome::genomeDetect(){
	mkdir(out_dir_detect.c_str(), S_IRWXU | S_IROTH);  // create the directory for detect command
	Chrome *chr;
	for(size_t i=0; i<chromeVector.size(); i++){
//		if(i==1)
		{
			chr = chromeVector.at(i);
			chr->chrDetect();
		}
	}

	// remove redundant translocations
	removeRedundantTra();

	// save detect result to file for each chrome
	saveDetectResultToFile();

	mergeDetectResult();
	return 0;
}

// remove repeatedly detected translocations
void Genome::removeRedundantTra(){
	size_t i, j, clipPosNum, clipPosNum_overlapped;;
	Chrome *chr;
	mateClipReg_t *clip_reg, *clip_reg_overlapped;
	for(i=0; i<chromeVector.size(); i++){
		chr = chromeVector.at(i);
		for(j=0; j<chr->mateClipRegVector.size(); j++){
			clip_reg = chr->mateClipRegVector.at(j);
			if(clip_reg->valid_flag and clip_reg->reg_mated_flag){
				clip_reg_overlapped = genomeGetOverlappedMateClipReg(clip_reg, chromeVector);
				if(clip_reg_overlapped){
					clipPosNum = clip_reg->leftClipPosNum + clip_reg->rightClipPosNum;
					clipPosNum_overlapped = clip_reg_overlapped->leftClipPosNum + clip_reg_overlapped->rightClipPosNum;
					if(clipPosNum>=clipPosNum_overlapped)
						clip_reg_overlapped->valid_flag = false;
					else
						clip_reg->valid_flag = false;
				}
			}
		}
	}

	// remove invalid elements
	for(i=0; i<chromeVector.size(); i++){
		chr = chromeVector.at(i);
		for(j=0; j<chr->mateClipRegVector.size(); ){
			clip_reg = chr->mateClipRegVector.at(j);
			if(clip_reg->valid_flag==false){
				delete clip_reg->leftClipReg;
				delete clip_reg->rightClipReg;
				delete clip_reg;
				chr->mateClipRegVector.erase(chr->mateClipRegVector.begin()+j);
			}else j++;
		}
	}
}

// get overlapped mate clip regions
mateClipReg_t* Genome::genomeGetOverlappedMateClipReg(mateClipReg_t *clip_reg_given, vector<Chrome*> &chrome_vec){
	mateClipReg_t *clip_reg_overlapped = NULL;
	Chrome *chr;
	for(size_t i=0; i<chrome_vec.size(); i++){
		chr = chrome_vec.at(i);
		clip_reg_overlapped = getOverlappedMateClipReg(clip_reg_given, chr->mateClipRegVector);
		if(clip_reg_overlapped)
			break;
	}
	return clip_reg_overlapped;
}

// save detect result to file for each chrome
void Genome::saveDetectResultToFile(){
	Chrome *chr;
	for(size_t i=0; i<chromeVector.size(); i++){
		chr = chromeVector.at(i);
		chr->chrMergeDetectResultToFile();
	}
}

// merge detect result to single file
void Genome::mergeDetectResult(){
	ofstream out_file_snv, out_file_indel, out_file_clipReg;

	out_file_snv.open(out_filename_detect_snv);
	if(!out_file_snv.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << out_filename_detect_snv << endl;
		exit(1);
	}
	out_file_indel.open(out_filename_detect_indel);
	if(!out_file_indel.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << out_filename_detect_indel << endl;
		exit(1);
	}
	out_file_clipReg.open(out_filename_detect_clipReg);
	if(!out_file_clipReg.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << out_filename_detect_clipReg << endl;
		exit(1);
	}

	vector<Chrome*>::iterator chr;
	for(chr=chromeVector.begin(); chr!=chromeVector.end(); chr++){
		copySingleFile((*chr)->out_filename_detect_indel, out_file_indel); // indel
		copySingleFile((*chr)->out_filename_detect_snv, out_file_snv); // snv
		copySingleFile((*chr)->out_filename_detect_clipReg, out_file_clipReg); // clip regions
	}

	out_file_snv.close();
	out_file_indel.close();
	out_file_clipReg.close();
}

// local assembly for genome
int Genome::genomeLocalAssemble(){
	mkdir(out_dir_assemble.c_str(), S_IRWXU | S_IROTH);  // create directory for assemble command
	for(size_t i=0; i<chromeVector.size(); i++){
		//if(i==2)
		{
			chromeVector.at(i)->chrLoadDataAssemble();  // load the variation data
			chromeVector.at(i)->chrLocalAssemble();     // local assembly
		}
	}
	return 0;
}

// call variants for genome
int Genome::genomeCall(){
	size_t i;
	Chrome *chr;

	mkdir(out_dir_call.c_str(), S_IRWXU | S_IROTH);  // create directory for call command

	// load data for call command
	loadCallData();

	// call variants
	for(i=0; i<chromeVector.size(); i++){
		chr = chromeVector.at(i);
		//if(i==1)
			chr->chrCall();
	}

	// call TRA according to mate clip regions
	genomeCallTra();

	// fill variation sequence
	genomeFillVarseq();

	// save SV to file
	genomeSaveCallSV2File();

	// merge call results into single file
	mergeCallResult();

	return 0;
}

// load data for call command
void Genome::loadCallData(){
	Chrome *chr;
	for(size_t i=0; i<chromeVector.size(); i++){
		chr = chromeVector.at(i);
		chr->chrLoadDataCall();
	}
}

// call TRA according to mate clip regions
void Genome::genomeCallTra(){
	size_t i, j, k, t, round_num;
	Chrome *chr;
	vector<mateClipReg_t*> mate_clipReg_vec;
	mateClipReg_t *clip_reg;
	varCand *var_cand, *var_cand_tmp;
	blat_aln_t *item_blat;
	bool call_success_flag;
	vector<int32_t> tra_loc_vec;

	for(i=0; i<chromeVector.size(); i++){
		chr = chromeVector.at(i);
		mate_clipReg_vec = chr->mateClipRegVector;
		for(j=0; j<mate_clipReg_vec.size(); j++){
			clip_reg = chr->mateClipRegVector.at(j);
			call_success_flag = false;
			if(clip_reg->valid_flag and clip_reg->sv_type==VAR_TRA){ // valid TRA
				for(round_num=0; round_num<2; round_num++){
					// construct var_cand
					if(round_num==0){
						var_cand = clip_reg->left_var_cand_tra;
						var_cand_tmp = constructNewVarCand(clip_reg->left_var_cand_tra, clip_reg->right_var_cand_tra);
					}else if(round_num==1){
						var_cand = clip_reg->right_var_cand_tra;
						var_cand_tmp = constructNewVarCand(clip_reg->right_var_cand_tra, clip_reg->left_var_cand_tra);
					}

					// load align data
					var_cand->loadBlatAlnData();
					var_cand_tmp->loadBlatAlnData();

					// compute TRA locations
					tra_loc_vec = computeTraLoc(var_cand, var_cand_tmp);
					if(tra_loc_vec.size()>0){
						call_success_flag = true;
						clip_reg->leftClipPosTra1 = tra_loc_vec.at(0);
						clip_reg->rightClipPosTra1 = tra_loc_vec.at(1);
						clip_reg->leftClipPosTra2 = tra_loc_vec.at(2);
						clip_reg->rightClipPosTra2 = tra_loc_vec.at(3);
					}

					// release memory
					for(k=0; k<var_cand_tmp->blat_aln_vec.size(); k++){ // free blat_aln_vec
						item_blat = var_cand_tmp->blat_aln_vec.at(k);
						for(t=0; t<item_blat->aln_segs.size(); t++) // free aln_segs
							delete item_blat->aln_segs.at(t);
						delete item_blat;
					}
					delete var_cand_tmp;

					if(call_success_flag) break;
				}
			}
		}
	}
}

varCand* Genome::constructNewVarCand(varCand *var_cand, varCand *var_cand_tmp){
	varCand *var_cand_new = NULL;
	if(var_cand and var_cand_tmp){
		var_cand_new = new varCand();
		var_cand_new->chrname = var_cand_tmp->chrname;
		var_cand_new->var_cand_filename = "";
		var_cand_new->out_dir_call = out_dir_call;
		var_cand_new->misAln_filename = "";
		var_cand_new->inBamFile = paras->inBamFile;
		var_cand_new->fai = fai;

		var_cand_new->refseqfilename = var_cand_tmp->refseqfilename;  // refseq file name
		var_cand_new->ctgfilename = var_cand->ctgfilename;  // contig file name
		var_cand_new->readsfilename = var_cand->readsfilename;  // reads file name
		var_cand_new->ref_left_shift_size = var_cand_tmp->ref_left_shift_size;  // ref_left_shift_size
		var_cand_new->ref_right_shift_size = var_cand_tmp->ref_right_shift_size;  // ref_right_shift_size

		var_cand_new->assem_success = var_cand->assem_success;
		var_cand_new->ctg_num = var_cand->ctg_num;

		for(size_t k=0; k<var_cand_tmp->varVec.size(); k++) var_cand_new->varVec.push_back(var_cand_tmp->varVec.at(k));
		var_cand_new->varVec.shrink_to_fit();

		var_cand_new->alnfilename = var_cand->alnfilename.substr(0, var_cand->alnfilename.size()-5) + "_" + var_cand_tmp->chrname + "_" + to_string(var_cand_tmp->leftClipRefPos) + "-" + to_string(var_cand_tmp->rightClipRefPos) + ".sim4";
		var_cand_new->align_success = false;
		var_cand_new->clip_reg_flag = true;

		var_cand_new->leftClipRefPos = var_cand_tmp->leftClipRefPos;
		var_cand_new->rightClipRefPos = var_cand_tmp->rightClipRefPos;
		var_cand_new->sv_type = var_cand_tmp->sv_type;
		var_cand_new->dup_num = var_cand_tmp->dup_num;
	}

	return var_cand_new;
}

// compute variant type
vector<int32_t> Genome::computeTraLoc(varCand *var_cand, varCand *var_cand_tmp){
	size_t i, j, aln_orient, missing_size_head, missing_size_inner, missing_size_tail;
	size_t start_query_pos, end_query_pos, start_mate_query_pos, end_mate_query_pos;
	int32_t maxValue, maxIdx, query_dist, query_part_flag;  // query_part_flag: 0-head, 1-inner, 2-tail
	int32_t left_tra_pos, right_tra_pos, tra_blat_aln_idx, start_tra_aln_seg_idx, end_tra_aln_seg_idx;
	blat_aln_t *blat_aln, *blat_aln_mate;
	aln_seg_t *aln_seg1, *aln_seg2;
	vector<int32_t> aln_idx_vec, mate_aln_idx_vec, tra_loc_vec;
	bool call_success_flag;

	call_success_flag = false;
	left_tra_pos = right_tra_pos = -1;
	for(i=0; i<var_cand->blat_aln_vec.size(); i++){
		blat_aln = var_cand->blat_aln_vec.at(i);
		if(blat_aln->valid_aln){
			aln_orient = blat_aln->aln_orient;

			// query head
			if(aln_orient==ALN_PLUS_ORIENT)
				missing_size_head = blat_aln->aln_segs.at(0)->query_start - 1;
			else
				missing_size_head = blat_aln->query_len - blat_aln->aln_segs.at(blat_aln->aln_segs.size()-1)->query_end;

			// query inner part
			maxValue = 0; maxIdx = -1;
			for(j=0; j<blat_aln->aln_segs.size()-1; j++){
				aln_seg1 = blat_aln->aln_segs.at(j);
				aln_seg2 = blat_aln->aln_segs.at(j+1);
				query_dist = aln_seg2->query_start - aln_seg1->query_end;
				if(query_dist>maxValue) {
					maxValue = query_dist;
					maxIdx = j;
				}
			}
			missing_size_inner = maxValue;

			// query tail
			if(aln_orient==ALN_PLUS_ORIENT)
				missing_size_tail = blat_aln->query_len - blat_aln->aln_segs.at(blat_aln->aln_segs.size()-1)->query_end;
			else
				missing_size_tail = blat_aln->aln_segs.at(0)->query_start - 1;

			// choose the maximum
			query_part_flag = 0; maxValue = missing_size_head;
			if(maxValue<(int32_t)missing_size_inner) { query_part_flag = 1; maxValue = missing_size_inner; }
			if(maxValue<(int32_t)missing_size_tail) { query_part_flag = 2; maxValue = missing_size_tail; }

			if(maxValue>MIN_CLIP_DIST_THRES){
				tra_blat_aln_idx = i;
				if(query_part_flag==0){ // query head
					start_query_pos = 1;
					end_query_pos = missing_size_head;
					start_tra_aln_seg_idx = -1;
					if(aln_orient==ALN_PLUS_ORIENT)  end_tra_aln_seg_idx = 0;
					else end_tra_aln_seg_idx = blat_aln->aln_segs.size() - 1;

					cout << "head: " << missing_size_head << ", " << var_cand->chrname << ":" << start_query_pos << "-" << end_query_pos << endl;
				}else if(query_part_flag==1){ // query inner
					aln_seg1 = blat_aln->aln_segs.at(maxIdx);
					aln_seg2 = blat_aln->aln_segs.at(maxIdx+1);
					if(aln_orient==ALN_PLUS_ORIENT){ // plus orient
						start_query_pos = aln_seg1->query_end + 1;
						end_query_pos = aln_seg2->query_start - 1;
					}else{ // minus orient
						start_query_pos = blat_aln->query_len - (aln_seg2->query_start - 1) + 1;
						end_query_pos = blat_aln->query_len - (aln_seg1->query_end + 1) + 1;
					}
					if(aln_orient==ALN_PLUS_ORIENT){
						start_tra_aln_seg_idx = maxIdx;
						end_tra_aln_seg_idx = maxIdx + 1;
					}else{
						start_tra_aln_seg_idx = maxIdx + 1;
						end_tra_aln_seg_idx = maxIdx;
					}

					cout << "inner: " << missing_size_inner << ", " << var_cand->chrname << ":" << start_query_pos << "-" << end_query_pos << endl;
				}else if(query_part_flag==2){ // query tail
					start_query_pos = blat_aln->query_len - missing_size_tail + 1;
					end_query_pos = blat_aln->query_len;

					if(aln_orient==ALN_PLUS_ORIENT) start_tra_aln_seg_idx = blat_aln->aln_segs.size() - 1;
					end_tra_aln_seg_idx = -1;

					cout << "tail: " << missing_size_tail << ", " << var_cand->chrname << ":" << start_query_pos << "-" << end_query_pos << endl;
				}
				aln_idx_vec.push_back(tra_blat_aln_idx);
				aln_idx_vec.push_back(start_tra_aln_seg_idx);
				aln_idx_vec.push_back(end_tra_aln_seg_idx);

				// get the missing region in the mated clipping region
				mate_aln_idx_vec = getMateTraReg(blat_aln->query_id, start_query_pos, end_query_pos, var_cand_tmp);
				if(mate_aln_idx_vec.at(0)!=-1){
					blat_aln_mate = var_cand_tmp->blat_aln_vec.at(mate_aln_idx_vec.at(0));
					aln_seg1 = blat_aln_mate->aln_segs.at(mate_aln_idx_vec.at(1));
					aln_seg2 = blat_aln_mate->aln_segs.at(mate_aln_idx_vec.at(2));
					if(blat_aln_mate->aln_orient==ALN_PLUS_ORIENT){ // plus orient
						start_mate_query_pos = aln_seg1->query_start;
						end_mate_query_pos = aln_seg2->query_end;
					}else{ // minus orient
						start_mate_query_pos = blat_aln_mate->query_len - aln_seg1->query_end + 1;
						end_mate_query_pos = blat_aln_mate->query_len - aln_seg2->query_start + 1;
					}

					cout << "mate TRA: " << "start_mate_query_pos=" << start_mate_query_pos << ", end_mate_query_pos=" << end_mate_query_pos << endl;

					// compute TRA clipping locations
					tra_loc_vec = computeTraClippingLoc(query_part_flag, var_cand, aln_idx_vec, var_cand_tmp, mate_aln_idx_vec);

					cout << "TRA location: " << var_cand->chrname << ":" << tra_loc_vec.at(0) << "-" << tra_loc_vec.at(1) << ", " << var_cand_tmp->chrname << ":" << tra_loc_vec.at(2) << "-" << tra_loc_vec.at(3) << endl;

					call_success_flag = computeTraCallSuccessFlag(tra_loc_vec);
					if(call_success_flag) break;
				}
			}
		}
	}

	return tra_loc_vec;
}

// get the missing region in the mated clipping region: [0]-blat_aln_id, [1]-start_aln_seg_id, [2]-end_aln_seg_id
vector<int32_t> Genome::getMateTraReg(size_t query_id, size_t start_query_pos, size_t end_query_pos, varCand *var_cand_tmp){
	vector<int32_t> mate_idx_vec;
	size_t i, aln_orient, start_pos_tmp, end_pos_tmp, start_mate_query_pos, end_mate_query_pos;
	int32_t j, blat_aln_idx, start_seg_idx, end_seg_idx;
	blat_aln_t *blat_aln;
	aln_seg_t *aln_seg;
	bool valid_tra_flag;

	for(i=0; i<var_cand_tmp->blat_aln_vec.size(); i++){
		blat_aln_idx = start_seg_idx = end_seg_idx = -1;
		start_mate_query_pos = end_mate_query_pos = 0;

		blat_aln = var_cand_tmp->blat_aln_vec.at(i);
		if(blat_aln->valid_aln and blat_aln->query_id==query_id){
			aln_orient = blat_aln->aln_orient;
			if(aln_orient==ALN_PLUS_ORIENT){ // plus orient
				// compute the start pos
				for(j=0; j<(int32_t)blat_aln->aln_segs.size(); j++){
					aln_seg = blat_aln->aln_segs.at(j);
					start_pos_tmp = aln_seg->query_start;
					end_pos_tmp = aln_seg->query_end;
					if(isOverlappedPos(start_pos_tmp, end_pos_tmp, start_query_pos, end_query_pos)){
						start_mate_query_pos = start_pos_tmp;
						start_seg_idx = j;
						break;
					}
				}

				// compute the end pos
				if(start_seg_idx!=-1){
					for(j=(int32_t)blat_aln->aln_segs.size()-1; j>=start_seg_idx; j--){
						aln_seg = blat_aln->aln_segs.at(j);
						start_pos_tmp = aln_seg->query_start;
						end_pos_tmp = aln_seg->query_end;
						if(isOverlappedPos(start_pos_tmp, end_pos_tmp, start_query_pos, end_query_pos)){
							end_mate_query_pos = end_pos_tmp;
							end_seg_idx = j;
							break;
						}
					}
				}
			}else{ // minus orient
				// compute the start pos
				for(j=(int32_t)blat_aln->aln_segs.size()-1; j>=0; j--){
					aln_seg = blat_aln->aln_segs.at(j);
					start_pos_tmp = blat_aln->query_len - aln_seg->query_end + 1;
					end_pos_tmp = blat_aln->query_len - aln_seg->query_start + 1;
					if(isOverlappedPos(start_pos_tmp, end_pos_tmp, start_query_pos, end_query_pos)){
						start_mate_query_pos = start_pos_tmp;
						start_seg_idx = j;
						break;
					}
				}

				// compute the end pos
				if(start_seg_idx!=-1){
					for(j=0; j<(int32_t)blat_aln->aln_segs.size(); j++){
						aln_seg = blat_aln->aln_segs.at(j);
						start_pos_tmp = blat_aln->query_len - aln_seg->query_end + 1;
						end_pos_tmp = blat_aln->query_len - aln_seg->query_start + 1;
						if(isOverlappedPos(start_pos_tmp, end_pos_tmp, start_query_pos, end_query_pos)){
							end_mate_query_pos = end_pos_tmp;
							end_seg_idx = j;
							break;
						}
					}
				}
			}

			if(start_seg_idx!=-1 and end_seg_idx!=-1){
				valid_tra_flag = isValidTraReg(start_query_pos, end_query_pos, start_mate_query_pos, end_mate_query_pos);
				if(valid_tra_flag){
					blat_aln_idx = i;
					break;
				}
			}else
				blat_aln_idx = start_seg_idx = end_seg_idx = -1;
		}
	}

	mate_idx_vec.push_back(blat_aln_idx);
	mate_idx_vec.push_back(start_seg_idx);
	mate_idx_vec.push_back(end_seg_idx);

	return mate_idx_vec;
}

// determine whether the TRA region is valid
bool Genome::isValidTraReg(size_t start_query_pos, size_t end_query_pos, size_t start_mate_query_pos, size_t end_mate_query_pos){
	bool valid_flag;
	double shared_ratio, shared_len, total_len;
	size_t minPos, maxPos, min_shared_pos, max_shared_pos;

	// compute the ratio of shared region
	if(start_query_pos<start_mate_query_pos){
		minPos = start_query_pos;
		min_shared_pos = start_mate_query_pos;
	}else{
		minPos = start_mate_query_pos;
		min_shared_pos = start_query_pos;
	}
	if(end_query_pos<end_mate_query_pos){
		maxPos = end_mate_query_pos;
		max_shared_pos = end_query_pos;
	}else{
		maxPos = end_query_pos;
		max_shared_pos = end_mate_query_pos;
	}
	shared_len = max_shared_pos - min_shared_pos + 1;
	total_len = maxPos - minPos + 1;
	shared_ratio = shared_len / total_len;

	cout << "shared_len=" << shared_len << ", total_len=" << total_len << ", shared_ratio=" << shared_ratio << endl;

	if(shared_ratio>0.9) valid_flag = true;
	else valid_flag = false;

	return valid_flag;
}

// compute TRA clipping locations
vector<int32_t> Genome::computeTraClippingLoc(size_t query_clip_part_flag, varCand *var_cand, vector<int32_t> &aln_idx_vec, varCand *var_cand_tmp, vector<int32_t> &mate_aln_idx_vec){
	vector<int32_t> tra_loc_vec;
	int32_t tra_refPos, left_tra_refPos1, right_tra_refPos1, left_tra_refPos2, right_tra_refPos2;
	blat_aln_t *blat_aln, *blat_aln_tmp;
	aln_seg_t *aln_seg1, *aln_seg2;

	tra_refPos = left_tra_refPos1 = right_tra_refPos1 = left_tra_refPos2 = right_tra_refPos2 = -1;
	blat_aln = var_cand->blat_aln_vec.at(aln_idx_vec.at(0));
	blat_aln_tmp = var_cand_tmp->blat_aln_vec.at(mate_aln_idx_vec.at(0));

	if(query_clip_part_flag==0){ // query head
		aln_seg1 = blat_aln->aln_segs.at(aln_idx_vec.at(2));
		left_tra_refPos1 = -1;
		if(blat_aln->aln_orient==ALN_PLUS_ORIENT)
			right_tra_refPos1 = aln_seg1->ref_start;
		else
			right_tra_refPos1 = aln_seg1->ref_end;
	}else if(query_clip_part_flag==1){ // query inner
		aln_seg1 = blat_aln->aln_segs.at(aln_idx_vec.at(1));
		aln_seg2 = blat_aln->aln_segs.at(aln_idx_vec.at(2));
		if(blat_aln->aln_orient==ALN_PLUS_ORIENT){
			left_tra_refPos1 = aln_seg1->ref_end;
			right_tra_refPos1 = aln_seg2->ref_start;
		}else{
			left_tra_refPos1 = aln_seg1->ref_start;
			right_tra_refPos1 = aln_seg2->ref_end;
		}
	}else if(query_clip_part_flag==2){ // query tail
		aln_seg2 = blat_aln->aln_segs.at(aln_idx_vec.at(1));
		if(blat_aln->aln_orient==ALN_PLUS_ORIENT)
			left_tra_refPos1 = aln_seg2->ref_end;
		else
			left_tra_refPos1 = aln_seg2->ref_start;
		right_tra_refPos1 = -1;
	}
	cout << "left_tra_refPos1=" << left_tra_refPos1 << ", right_tra_refPos1=" << right_tra_refPos1 << endl;

	if(query_clip_part_flag==0){ // query head
		aln_seg2 = blat_aln_tmp->aln_segs.at(mate_aln_idx_vec.at(2));
		left_tra_refPos2 = -1;
		if(blat_aln_tmp->aln_orient==ALN_PLUS_ORIENT){ // plus orient
			right_tra_refPos2 = aln_seg2->ref_end;
		}else{ // minus orient
			right_tra_refPos2 = aln_seg2->ref_start;
		}
	}else if(query_clip_part_flag==1){ // query inner
		aln_seg1 = blat_aln_tmp->aln_segs.at(mate_aln_idx_vec.at(1));
		aln_seg2 = blat_aln_tmp->aln_segs.at(mate_aln_idx_vec.at(2));
		if(blat_aln_tmp->aln_orient==ALN_PLUS_ORIENT){ // plus orient
			left_tra_refPos2 = aln_seg1->ref_start;
			right_tra_refPos2 = aln_seg2->ref_end;
		}else{ // minus orient
			left_tra_refPos2 = aln_seg1->ref_end;
			right_tra_refPos2 = aln_seg2->ref_start;
		}
	}else if(query_clip_part_flag==2){ // query tail
		aln_seg1 = blat_aln_tmp->aln_segs.at(mate_aln_idx_vec.at(1));
		if(blat_aln_tmp->aln_orient==ALN_PLUS_ORIENT){ // plus orient
			left_tra_refPos2 = aln_seg1->ref_start;
		}else{ // minus orient
			left_tra_refPos2 = aln_seg1->ref_end;
		}
		right_tra_refPos2 = -1;
	}
	cout << "left_tra_refPos2=" << left_tra_refPos2 << ", right_tra_refPos2=" << right_tra_refPos2 << endl;

	tra_loc_vec.push_back(left_tra_refPos1);
	tra_loc_vec.push_back(right_tra_refPos1);
	tra_loc_vec.push_back(left_tra_refPos2);
	tra_loc_vec.push_back(right_tra_refPos2);

	return tra_loc_vec;
}


// determine whether the TRA is called successfully
bool Genome::computeTraCallSuccessFlag(vector<int32_t> &tra_loc_vec){
	bool flag = false;
	if((tra_loc_vec.at(0)!=-1 or tra_loc_vec.at(1)!=-1) and (tra_loc_vec.at(2)!=-1 or tra_loc_vec.at(3)!=-1)){
		flag = true;
	}
	return flag;
}

// fill variant sequences
void Genome::genomeFillVarseq(){
	Chrome *chr;
	for(size_t i=0; i<chromeVector.size(); i++){
		chr = chromeVector.at(i);
		chr->chrFillVarseq();
	}
}

// save variants to file
void Genome::genomeSaveCallSV2File(){
	Chrome *chr;
	for(size_t i=0; i<chromeVector.size(); i++){
		chr = chromeVector.at(i);
		chr->saveCallSV2File();
	}

	// save TRA
	saveTraCall2File();
}

// save TRA
void Genome::saveTraCall2File(){
	size_t i, j;
	Chrome *chr;
	ofstream outfile_tra;
	vector<mateClipReg_t*> mate_clipReg_vec;
	mateClipReg_t *clip_reg;
	string line;
	int32_t sv_len;

	outfile_tra.open(out_filename_call_tra);
	if(!outfile_tra.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << out_filename_call_tra << endl;
		exit(1);
	}

	for(i=0; i<chromeVector.size(); i++){
		chr = chromeVector.at(i);
		mate_clipReg_vec = chr->mateClipRegVector;
		for(j=0; j<mate_clipReg_vec.size(); j++){
			clip_reg = mate_clipReg_vec.at(j);
			if(clip_reg->valid_flag and clip_reg->sv_type==VAR_TRA and (clip_reg->leftClipPosTra1>0 or clip_reg->rightClipPosTra1>0) and (clip_reg->leftClipPosTra2>0 or clip_reg->rightClipPosTra2>0)){
				line = clip_reg->left_var_cand_tra->chrname;
				if(clip_reg->leftClipPosTra1>0) line += "\t" + to_string(clip_reg->leftClipPosTra1);
				else line += "\t-";
				if(clip_reg->rightClipPosTra1>0) line += "\t" + to_string(clip_reg->rightClipPosTra1);
				else line += "\t-";
				line += "\t" + clip_reg->right_var_cand_tra->chrname;
				if(clip_reg->leftClipPosTra2>0) line += "\t" + to_string(clip_reg->leftClipPosTra2);
				else line += "\t-";
				if(clip_reg->rightClipPosTra2>0) line += "\t" + to_string(clip_reg->rightClipPosTra2);
				else line += "\t-";
				if(clip_reg->refseq_tra.size()>0) line += "\t" + clip_reg->refseq_tra;
				else line += "\t-";
				if(clip_reg->altseq_tra.size()>0) line += "\t" + clip_reg->altseq_tra;
				else line += "\t-";
				line += "\tTRA";

				if(clip_reg->leftClipPosTra1>0 and clip_reg->rightClipPosTra1>0){
					sv_len = clip_reg->rightClipPosTra1 - clip_reg->leftClipPosTra1 + 1;
					line += "\t" + to_string(sv_len);
				}else line += "\t-";
				if(clip_reg->leftClipPosTra2>0 and clip_reg->rightClipPosTra2>0){
					sv_len = clip_reg->rightClipPosTra2 - clip_reg->leftClipPosTra2 + 1;
					line += "\t" + to_string(sv_len);
				}else line += "\t-";

				outfile_tra << line << endl;
				cout << line << endl;
			}
		}
	}

	outfile_tra.close();
}

// merge call results into single file
void Genome::mergeCallResult(){
	ofstream out_file_indel;

	out_file_indel.open(out_filename_call_indel);
	if(!out_file_indel.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << out_filename_call_indel << endl;
		exit(1);
	}

	for(size_t i=0; i<chromeVector.size(); i++)
		copySingleFile(chromeVector.at(i)->out_filename_call_indel, out_file_indel); // indel

	out_file_indel.close();
}

// save results in VCF file format for detect command
void Genome::saveResultVCF(){
	string outIndelFile, outSnvFile;

	cout << "Output results to file ..." << endl;

	// initialize the output filenames
	if(paras->outFilePrefix.size()>0) {
		outIndelFile = out_dir_call + "/" + paras->outFilePrefix + "_INDEL.vcf";
		outSnvFile = out_dir_call + "/" + paras->outFilePrefix + "_SNV.vcf";
	}else {
		outIndelFile = out_dir_call + "/" + "INDEL.vcf";
		outSnvFile = out_dir_call + "/" + "SNV.vcf";
	}

	// save results
	//saveIndelVCFDetect(out_filename_detect_indel, outIndelFile);
	//saveSnvVCFDetect(out_filename_detect_snv, outSnvFile);
}

// save results in VCF file format for detect command
void Genome::saveResultVCFDetect(){
	string outIndelFile, outSnvFile;

	cout << "Output results to file ..." << endl;

	// initialize the output filenames
	if(paras->outFilePrefix.size()>0) {
		outIndelFile = out_dir_detect + "/" + paras->outFilePrefix + "_INDEL.vcf";
		outSnvFile = out_dir_detect + "/" + paras->outFilePrefix + "_SNV.vcf";
	}else {
		outIndelFile = out_dir_detect + "/" + "INDEL.vcf";
		outSnvFile = out_dir_detect + "/" + "SNV.vcf";
	}

	// save results
	saveIndelVCFDetect(out_filename_detect_indel, outIndelFile);
	saveSnvVCFDetect(out_filename_detect_snv, outSnvFile);
}

// save indel in VCF file format for detect command
void Genome::saveIndelVCFDetect(string &in, string &out_vcf){
	string line, data_vcf, ref, reg;
	ifstream infile;
	ofstream outfile;
	vector<string> str_vec;
	char *seq;
	int seq_len;

	// open files
	infile.open(in);
	if(!infile.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << in << endl;
		exit(1);
	}

	outfile.open(out_vcf);
	if(!outfile.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << out_vcf << endl;
		exit(1);
	}

	// save VCF header
	saveVCFHeader(outfile);

	// save results
	while(getline(infile, line))
		if(line.size()>0){
			str_vec = split(line, "\t");
			data_vcf = str_vec[0] + "\t" + str_vec[1] + "\t."; // CHROM, POS, ID

			reg = str_vec[0] + ":" + str_vec[1] + "-" + str_vec[2];
			seq = fai_fetch(fai, reg.c_str(), &seq_len);
			data_vcf += "\t" + (string)seq + "\t."; // REF, ALT

			// QUAL, FILTER, INFO0, FORMAT ........


			outfile << data_vcf << endl;

			free(seq);
		}
	infile.close();
	outfile.close();
}

// save SNV in VCF file format for detect command
void Genome::saveSnvVCFDetect(string &in, string &out_vcf){
	string line, data_vcf, ref, reg;
	ifstream infile;
	ofstream outfile;
	vector<string> str_vec;
	char *seq;
	int seq_len;

	// open files
	infile.open(in);
	if(!infile.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << in << endl;
		exit(1);
	}

	outfile.open(out_vcf);
	if(!outfile.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << out_vcf << endl;
		exit(1);
	}

	// save VCF header
	saveVCFHeader(outfile);

	// save results
	while(getline(infile, line))
		if(line.size()>0){
			str_vec = split(line, "\t");
			data_vcf = str_vec[0] + "\t" + str_vec[1] + "\t."; // CHROM, POS, ID

			reg = str_vec[0] + ":" + str_vec[1] + "-" + str_vec[1];
			seq = fai_fetch(fai, reg.c_str(), &seq_len);
			data_vcf += "\t" + (string)seq + "\t."; // REF, ALT

			// QUAL, FILTER, INFO0, FORMAT ........


			outfile << data_vcf << endl;

			free(seq);
		}
	infile.close();
	outfile.close();
}

// save VCF header
void Genome::saveVCFHeader(ofstream &fp){
	time_t rawtime;
	struct tm *timeinfo;
	char buffer [128];

	// get the local time
	time(&rawtime);
	timeinfo = localtime (&rawtime);
	strftime(buffer, sizeof(buffer), "%Y%m%d", timeinfo);

	fp << "##fileformat=VCFv4.2" << endl;
	fp << "##fileDate=" << buffer << endl;
	fp << "##source=" << PROG_NAME << "V" << PROG_VERSION << endl;
	fp << "##reference=" << paras->refFile << endl;
	fp << "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">" << endl;
	fp << "##FILTER=<ID=q10,Description=\"Quality below 10\">" << endl;
	fp << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl;
	fp << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT" << endl;

}

