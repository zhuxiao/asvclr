#include <math.h>
#include "varCand.h"
#include "covLoader.h"
#include "util.h"


pthread_mutex_t mutex_print_var_cand = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mutex_write_minimap2_aln = PTHREAD_MUTEX_INITIALIZER;
extern pthread_mutex_t mutex_fai;
extern pthread_mutex_t mutex_write;

varCand::varCand(){
	this->chrname = "";
	this->var_cand_filename = "";
	this->out_dir_call = "";
	this->out_dir_cns = "";
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

	if(!minimap2_aln_vec.empty()) destroyMinimap2AlnVec(minimap2_aln_vec);
}

void varCand::init(){
	refseqfilename = "";
	ctgfilename = "";
	readsfilename = "";
	alnfilename = "";
	clusterfilename = "";
	rescue_refseqfilename = rescue_cnsfilename = rescue_readsfilename = rescue_alnfilename = "";
	technology = SEQUENCING_TECH_DEFAULT;
	ref_left_shift_size = ref_right_shift_size = ctg_num = 0;
	cns_success = align_success = call_success = clip_reg_flag = false;
	clip_reg = clip_reg_allele = NULL;
	leftClipRefPos = rightClipRefPos = 0;
	sv_type = VAR_UNC;
	dup_num = 0;
	margin_adjusted_flag = false;
	minimap2_var_cand_file = NULL;
	limit_reg_delete_flag = false;
	killed_flag = false;
	minClipEndSize = MIN_CLIP_END_SIZE;
	max_ultra_high_cov = MAX_ULTRA_HIGH_COV_THRES;
	minMapQ = MIN_MAPQ_THRES;
	minHighMapQ = MIN_HIGH_MAPQ_THRES;

	max_seg_size_ratio = max_seg_nm_ratio = max_absig_density = -1;

	gt_min_sig_size = GT_SIG_SIZE_THRES;
	gt_size_ratio_match = GT_SIZE_RATIO_MATCH_THRES;
	gt_min_seqsim_merge = GT_MIN_SEQSIM_MERGE_THRES_CCS;
	gt_homo_ratio_thres = GT_HOMO_RATIO_THRES;
	gt_hete_ratio_thres = GT_HETE_RATIO_THRES;
}

void varCand::destroyVarCand(){
	size_t i;

	for(i=0; i<varVec.size(); i++)  // free varVec
		delete varVec.at(i);
	for(i=0; i<newVarVec.size(); i++)  // free newVarVec
		delete newVarVec.at(i);

	//delete minimap2_aln_vec
	destroyMinimap2AlnVec(minimap2_aln_vec);
}

void varCand::destroyMinimap2AlnVec(vector<minimap2_aln_t*> &minimap2_aln_vec){
	minimap2_aln_t *minimap2_aln_node;
	for(size_t i=0; i<minimap2_aln_vec.size(); i++){
		minimap2_aln_node = minimap2_aln_vec.at(i);
		struct pafalnSeg *pafalnSeg_node;
		for(size_t j=0; j<minimap2_aln_node->pafalnsegs.size();j++){
			pafalnSeg_node = minimap2_aln_node->pafalnsegs.at(j);
			delete pafalnSeg_node;
		}
		vector<struct pafalnSeg*>().swap(minimap2_aln_node->pafalnsegs);
		delete minimap2_aln_node;
	}
	vector<minimap2_aln_t*>().swap(minimap2_aln_vec);
}

// destroy the alignment data of the block
void varCand::destoryClipAlnData(vector<clipAlnData_t*> &clipAlnDataVector){
	for(size_t i=0; i<clipAlnDataVector.size(); i++){
		bam_destroy1(clipAlnDataVector.at(i)->bam);
		delete clipAlnDataVector.at(i);
	}
	vector<clipAlnData_t*>().swap(clipAlnDataVector);
}

void varCand::setGtParas(int32_t gt_min_sig_size, double gt_size_ratio_match, double gt_min_seqsim_merge, double gt_homo_ratio_thres, double gt_hete_ratio_thres, int32_t gt_min_sup_num_recover){
	this->gt_min_sig_size = gt_min_sig_size;
	this->gt_size_ratio_match = gt_size_ratio_match;
	this->gt_min_seqsim_merge = gt_min_seqsim_merge;
	this->gt_homo_ratio_thres = gt_homo_ratio_thres;
	this->gt_hete_ratio_thres = gt_hete_ratio_thres;
	this->gt_min_sup_num_recover = gt_min_sup_num_recover;
}

// align contig to refseq
void varCand::alnCtg2Refseq(){
	//cout << "Minimap2 align: " << alnfilename << endl;

	bool minimap2_aln_done_flag = getMinimap2AlnDoneFlag();
	if(cns_success and minimap2_aln_done_flag==false){
		if(!isFileExist(alnfilename) or !isMinimap2AlnResultMatch(ctgfilename, alnfilename)){ // file not exist, or query names not match
			minimap2Aln(alnfilename, ctgfilename, refseqfilename, max_proc_running_minutes); // MINIMAP2 alignment
		}

		// record minimap2 aligned information
		recordMinimap2AlnInfo();
	}
}

// record minimap2 aligned information
void varCand::recordMinimap2AlnInfo(){
	string line, aln_flag_str, limit_reg_str;
	bool aln_flag;

	if(minimap2_var_cand_file){
		aln_flag = isFileExist(alnfilename);
		if(aln_flag) aln_flag_str = ALIGN_SUCCESS;
		else aln_flag_str = ALIGN_FAILURE;

		// limit regions
		if(limit_reg_process_flag) limit_reg_str = getLimitRegStr(sub_limit_reg_vec);
		else limit_reg_str = LIMIT_REG_ALL_STR;

		line = alnfilename + "\t" + ctgfilename + "\t" + refseqfilename + "\t" + aln_flag_str + "\t" + limit_reg_str + "\t" + DONE_STR;

		pthread_mutex_lock(&mutex_write_minimap2_aln);
		*minimap2_var_cand_file << line << endl;
		pthread_mutex_unlock(&mutex_write_minimap2_aln);
	}
}

// get minimap2 align done flag
bool varCand::getMinimap2AlnDoneFlag(){
	bool minimap2_aln_done_flag;
	varCand *minimap2_aligned_info_node;

	minimap2_aln_done_flag = false;
	if(minimap2_aligned_info_vec){
		for(size_t i=0; i<minimap2_aligned_info_vec->size(); i++){
			minimap2_aligned_info_node = minimap2_aligned_info_vec->at(i);
			if(minimap2_aligned_info_node->alnfilename.compare(alnfilename)==0){
				minimap2_aln_done_flag = true;
				break;
			}
		}
	}

	return minimap2_aln_done_flag;
}

// call variants according to MINIMAP alignment
void varCand::callVariants(){
	if(sv_type!=VAR_TRA){
//		pthread_mutex_lock(&mutex_print_var_cand);
//		cout << "process: " << alnfilename << endl;
//		pthread_mutex_unlock(&mutex_print_var_cand);

		loadMinimap2AlnData(); // load align data

		if(clip_reg_flag==false){	// call indel variants
			callIndelVariants();
		}else{ // call clipping variants
			callClipRegVariants();
		}

//		if(call_success==false){
//			cout << "\tFail: " << alnfilename << ", ctg_num=" << ctg_num << ", clip_reg_flag=" << clip_reg_flag << endl;
//		}

		if(!minimap2_aln_vec.empty()) destroyMinimap2AlnVec(minimap2_aln_vec);
	}
}

// load align data
void varCand::loadMinimap2AlnData(){
	// minimap2 align
#if MINIMAP2_ALN
	alnCtg2Refseq();
#endif
	assignMinimap2AlnStatus();
	if(align_success){
		// parse alignment
		minimap2Parse();
	}else{
//		pthread_mutex_lock(&mutex_print_var_cand);
//		cout << "Align Failure: " << ctgfilename << endl;
//		pthread_mutex_unlock(&mutex_print_var_cand);
	}
}

//call indel variants from paf
void varCand::callIndelVariants(){
	vector<reg_t*> var_vec, var_vec_rescue;
	vector<int32_t> clusterId_incomplete;

	if(align_success){
		var_vec = computeIndelVarLoc(minimap2_aln_vec, ctgfilename, refseqfilename, false);

		// determine whether the allelic variants are called completed
		clusterId_incomplete = getClusterIdCalledIncomplete(var_vec, clusterfilename, minimap2_aln_vec);

		// rescue the variant by wtdbg2
		//if(var_vec.size()==0)
		if(clusterId_incomplete.size()>0){
			var_vec_rescue = rescueIndelVarLoc(clusterId_incomplete);

			// update cluster id incomplete
			updateClusterIdIncomplete(var_vec, var_vec_rescue, clusterId_incomplete);
		}
	}else
		clusterId_incomplete = getClusterIdCalledIncomplete(var_vec, clusterfilename, minimap2_aln_vec);

	 // call from reads
	//if(var_vec.size()==0){
	if(clusterId_incomplete.size()>0){
		var_vec_rescue = callIndelFromReadsIndelReg(clusterId_incomplete);

		// update cluster id incomplete
		updateClusterIdIncomplete(var_vec, var_vec_rescue, clusterId_incomplete);
	}

	// genotyping
	if(var_vec.size()>0){
		call_success = true;
		computeGenotypeIndelReg(var_vec);

		// refine allele frequency
		refineAlleleFreq(newVarVec);

		// sort variants
		sortRegVec(newVarVec);
	}
}

// parse and call clipReg variants
void varCand::callClipRegVariants(){
	vector<reg_t*> var_vec, var_vec_rescue, new_var_vec;
	vector<int32_t> clusterId_incomplete;

	if(align_success){
		// call variants according to variant types
		var_vec = computeClipRegVarLoc(alnfilename, ctgfilename, refseqfilename, clusterfilename, minimap2_aln_vec, NULL, false);

		// determine whether the allelic variants are called completed
		clusterId_incomplete = getClusterIdCalledIncomplete(var_vec, clusterfilename, minimap2_aln_vec);

		//if(var_vec.size()==0){ // rescue the calling
		if(clusterId_incomplete.size()>0){
			if(large_indel_flag==false){ // clipping
				//cout << "\tRescue DUPs and INVs for '" << ctgfilename << "'..." << endl;
				var_vec_rescue = rescueDupInvClipReg(clusterId_incomplete);

			}else{ // large indel
				//cout << "\tRescue large indels for '" << ctgfilename << "'..." << endl;
				var_vec_rescue = rescueLargeIndelClipReg(clusterId_incomplete);
			}
			// update cluster id incomplete
			updateClusterIdIncomplete(var_vec, var_vec_rescue, clusterId_incomplete);
		}
	}else
		clusterId_incomplete = getClusterIdCalledIncomplete(var_vec, clusterfilename, minimap2_aln_vec);

	 // call from reads
	//if(var_vec.size()==0){
	if(clusterId_incomplete.size()>0){
		var_vec_rescue = callVarFromReadsClipReg(clusterId_incomplete);

		// update cluster id incomplete
		updateClusterIdIncomplete(var_vec, var_vec_rescue, clusterId_incomplete);
	}

	if(var_vec.size()>0){
		call_success = true;
		computeGenotypeClipReg(var_vec);

		refineAlleleFreq(newVarVec);

		sortRegVec(newVarVec);
	}
}

// update align status
void varCand::assignMinimap2AlnStatus(){
	align_success = isFileExist(alnfilename);
}

void varCand::minimap2Parse(){
	if(minimap2_aln_vec.empty())
		minimap2_aln_vec = minimap2Parse(alnfilename, ctgfilename, refseqfilename);
}

// parse Minimap2 alignments (paf format)
vector<minimap2_aln_t*> varCand::minimap2Parse(string &alnfilename, string &ctgfilename, string &refseqfilename){
	vector<minimap2_aln_t*> minimap2_aln_vec;
	string line, ref_name, ref_name_split, ref_region, relative_srtand, cigar_tag, cigar, chrname, cons_seq, ref_seq, type, type_tag;
	vector<string> line_vec, line_vec1, qname_vec, ref_name_vec1, ref_name_vec2, ref_name_vec3, cigar_tag_split, query_name_vec, type_tag_split;
	ifstream infile;
	int32_t i, ref_start_all, ref_end_all, query_loc, minimap2_aln_id;//r_dis, q_dis;//subject_dis
	minimap2_aln_t *minimap2_aln_item;
	vector<struct pafalnSeg*> pafalnsegs;

	// load query names
	FastaSeqLoader fa_loader(ctgfilename);
	query_name_vec = fa_loader.getFastaSeqNames();

	// load query names
	FastaSeqLoader ref_loader(refseqfilename);
	ref_seq = ref_loader.getFastaSeq(0);

	infile.open(alnfilename);
	if(!infile.is_open()){
		cerr << __func__ << "(), line=" << __LINE__ << ": cannot open file:" << alnfilename << endl;
		exit(1);
	}

	// fill the data
	minimap2_aln_id = 0;
	i = 0;
	while(getline(infile, line)){
		if(line.size()){
			line_vec = split(line, "\t");
			// get the query name
			query_loc = getQueryNameLoc(line_vec.at(0), query_name_vec);
			if(query_loc==-1){
				//cerr << __func__ << "(), line=" << __LINE__ << ": cannot get queryloc by query_name: " << line_vec.at(0) << endl;
				//exit(1);
				continue;
			}
			i += 1;

			// get ref_start_all
			ref_name = line_vec.at(5);
			ref_name_vec1 = split(ref_name, ":");
			chrname = ref_name_vec1.at(0);
			ref_name_split = ref_name_vec1.at(1);
			ref_name_vec2 = split(ref_name_split, "___");
			ref_region = ref_name_vec2.at(0);
			ref_name_vec3 = split(ref_region, "-");
			ref_start_all = stoi(ref_name_vec3.at(0)) - 1;
			ref_end_all = stoi(ref_name_vec3.at(1)) - 1;

			//get qname vec
			line_vec1 = split(line_vec.at(0), "___");
			qname_vec = split(line_vec1.at(1), ";");
			//for(j=0; j<(int32_t).size(); j++) qname_vec.push_back(line_vec1.at(j));

			//query_len = stoi(line_vec.at(1));

			minimap2_aln_item = new minimap2_aln_t;
			minimap2_aln_item->minimap2_aln_id = minimap2_aln_id++;
			minimap2_aln_item->query_len = stoi(line_vec.at(1));
			minimap2_aln_item->query_id = query_loc;
			minimap2_aln_item->chrname = chrname;
			//minimap2_aln_item->query_start = stoi(line_vec.at(2)) + 1; //start from 1
			minimap2_aln_item->query_start = stoi(line_vec.at(2)); //start from 0
			minimap2_aln_item->query_end = stoi(line_vec.at(3));
			minimap2_aln_item->region_startRefPos = ref_start_all;
			minimap2_aln_item->region_endRefPos = ref_end_all;
			minimap2_aln_item->qname = qname_vec;

			//get relative strand
			relative_srtand = line_vec.at(4);
			if(relative_srtand.compare("+")==0){
				minimap2_aln_item->relative_strand = ALN_PLUS_ORIENT;
			}else{// if(relative_srtand.compare("-")==0)
				minimap2_aln_item->relative_strand = ALN_MINUS_ORIENT;
			}
			minimap2_aln_item->subject_len = stoi(line_vec.at(6));
			//minimap2_aln_item->subject_start = stoi(line_vec.at(7)) + 1; //start from 1
			minimap2_aln_item->subject_start = stoi(line_vec.at(7)); //start from 0
			minimap2_aln_item->subject_end = stoi(line_vec.at(8));
			//match_base_num = stoi(line_vec[9]);
			//match_ref_len = stoi(line_vec[10]);
			//match_base_num = minimap2_aln_item->query_end - minimap2_aln_item->query_start + 1;
			//match_ref_len = minimap2_aln_item->subject_end - minimap2_aln_item->subject_start + 1;
			//minimap2_aln_item->ref_start = ref_start_all + minimap2_aln_item->subject_start - 1;
			//minimap2_aln_item->ref_end = ref_end_all + minimap2_aln_item->subject_end -1;

			//get type
			type = "", type_tag = "";
			type_tag = line_vec.at(16);
			type_tag_split = split(type_tag, ":");
			if(type_tag_split.at(0).compare("tp")==0){
				type = type_tag_split.at(2);
				minimap2_aln_item->type = type;
			}else{
				cerr << "In " << __func__ << "(), cannot find 'tp' tag in PAF alignment file, error!" << endl;
				exit(1);
			}

			//get cigar
			cigar_tag = line_vec.at(line_vec.size()-1);
			cigar_tag_split = split(cigar_tag, ":");
			if(cigar_tag_split.at(0).compare("cg")==0){
				cigar = cigar_tag_split.at(2);
				minimap2_aln_item->cigar = cigar;
			}else{
				cerr << "In " << __func__ << "(), cannot find 'cg' tag in PAF alignment file, error!" << endl;
				exit(1);
			}

			cons_seq = fa_loader.getFastaSeq(query_loc, minimap2_aln_item->relative_strand);

			//parse cigar and generate pafalnsegs
			pafalnsegs = generatePafAlnSegs(minimap2_aln_item, cons_seq, ref_seq);
			minimap2_aln_item->pafalnsegs = pafalnsegs;

//			//2022.9.21 miss DEL
//			if(pafalnsegs.size()>0)
//				getMissingPafInDelsAtAlnSegEnd(minimap2_aln_item->pafalnsegs, pafalnsegs, minimap2_aln_item->region_startRefPos, minimap2_aln_item->region_endRefPos, minimap2_aln_item->query_start, minimap2_aln_item->query_end, minimap2_aln_item->query_len, minimap2_aln_item->subject_len, minimap2_aln_item->subject_start, minimap2_aln_item->subject_end, match_base_num, match_ref_len, cons_seq, ref_seq, minimap2_aln_item->relative_strand);
			minimap2_aln_vec.push_back(minimap2_aln_item);
		}
	}
	infile.close();

	return minimap2_aln_vec;
}
//check the missing large InDels at alignment ends
void varCand::getMissingPafInDelsAtAlnSegEnd(vector<struct pafalnSeg*> &all_pafalnsegs, vector<struct pafalnSeg*> &pafalnsegs, int32_t region_refstart, int32_t region_refend, int32_t querystart, int32_t queryend, int32_t query_len, int32_t subjectlen, int32_t subjectstart, int32_t subjectend, int32_t match_base_num, int32_t match_ref_len, string cons_seq, string ref_seq, int32_t aln_orient){
	//struct pafalnSeg* seg;
	int32_t query_misslen, ref_misslen, miss_del_len, miss_ins_len, left_qdis, right_qdis, miss_refstart, miss_startsubject, left_refdis, right_refdis, miss_querystart;
	uint32_t opflag;

	query_misslen = query_len - match_base_num;
	ref_misslen = subjectlen - match_ref_len;
	if(ref_misslen > query_misslen){//DEL
		opflag = BAM_CDEL;
		miss_del_len = ref_misslen - query_misslen;
		if(ref_misslen > miss_del_len){
			left_qdis = querystart - 1;
			right_qdis = query_len - queryend;
			//if(left_qdis>200 or right_qdis>200){
			left_refdis = subjectstart - 1;
			right_refdis = subjectlen - subjectend;
			if(left_qdis>10 or right_qdis>10){
				if(right_refdis > left_refdis){
					miss_refstart = pafalnsegs.at(pafalnsegs.size() - 1)->startRpos + pafalnsegs.at(pafalnsegs.size() - 1)->seglen;
					miss_startsubject = subjectend;
					//subject_dis = minimap2_aln_item->subject_len - miss_startsubject;
					if(right_refdis>=miss_del_len){
						miss_querystart = queryend;
						all_pafalnsegs.push_back(allocatePafAlnSeg(miss_refstart, miss_querystart, miss_startsubject, miss_del_len, opflag, cons_seq, ref_seq, aln_orient));
					}
				}else{
					miss_refstart = pafalnsegs.at(0)->startRpos - miss_del_len;
					//miss_startsubject = subjectstart - miss_del_len;
					miss_startsubject = left_refdis - miss_del_len;
					if(miss_refstart>region_refstart and miss_refstart<region_refend and miss_startsubject>1){
						miss_querystart = querystart;
						all_pafalnsegs.push_back(allocatePafAlnSeg(miss_refstart, miss_querystart, miss_startsubject, miss_del_len, opflag, cons_seq, ref_seq, aln_orient));
					}
				}
			}
		}
	}else{
		if(ref_misslen < query_misslen){//ins
			if(ref_misslen>0){
				opflag = BAM_CINS;
				miss_ins_len = query_misslen - ref_misslen;
				//2022.10.16 left_qdis = querystart - 1;
				if(query_misslen>miss_ins_len){
					left_refdis = subjectstart - 1;
					right_refdis = subjectlen - subjectend;
					left_qdis = querystart - 1;
					right_qdis = query_len - queryend;
					if(left_refdis>10 or right_refdis>10){
						//if(left_qdis>200 or right_qdis>200){
						if(left_qdis>right_qdis){
							//miss_querystart = querystart - miss_ins_len;
							miss_querystart = left_qdis - miss_ins_len;
							if(miss_querystart>1){
								miss_refstart = pafalnsegs.at(0)->startRpos;
								miss_startsubject = subjectstart;
								all_pafalnsegs.push_back(allocatePafAlnSeg(miss_refstart, miss_querystart, miss_startsubject, miss_ins_len, opflag, cons_seq, ref_seq, aln_orient));
							}
						}else{
							miss_querystart = queryend;
							if(right_qdis>=miss_ins_len){
								miss_refstart = pafalnsegs.at(pafalnsegs.size() - 1)->startRpos + pafalnsegs.at(pafalnsegs.size() - 1)->seglen;
								miss_startsubject = subjectend;
								all_pafalnsegs.push_back(allocatePafAlnSeg(miss_refstart, miss_querystart, miss_startsubject, miss_ins_len, opflag, cons_seq, ref_seq, aln_orient));
							}
						}
					}
				}
			}
		}
	}
}

vector<reg_t*> varCand::computeIndelVarLoc(vector<minimap2_aln_t*> &minimap2_aln_vec, string &ctgfilename_para, string &refseqfilename_para, bool rescue_flag){
	vector<reg_t*> var_vec, var_vec_sub;
	size_t  i, j;
	int32_t k, overlap_size, sv_size_sum, dist, seq_len;
	uint32_t op;
	char *p_seq;
	vector<clipAlnData_t*> clipAlnDataVector;
	string chrname_tmp, gt_header, gt_str, dp_str, ad_str1, ad_str2, reg_str, refseq;
	minimap2_aln_t *minimap2_aln;
	vector<int32_t> supp_num_vec;
	struct pafalnSeg *paf_alnseg, *paf_alnseg2;
	bool flag, pre_large_match_flag;
	int64_t startRefPos_cns, endRefPos_cns, end_ref_pos, end_ref_pos2, start_var_pos, end_var_pos, start_var_pos_pure, end_var_pos_pure, chrlen_tmp;
	double ratio;
	reg_t *reg;
	vector<string> qname_vec;
	vector<struct pafalnSeg*> cand_paf_alnseg_vec, cand_paf_alnseg_vec_sub;
	vector<int32_t> cand_mm2_aln_idx_vec, cand_mm2_aln_idx_vec_sub, recovered_qid_vec;

	chrlen_tmp = faidx_seq_len64(fai, chrname.c_str());
	start_var_pos_pure = varVec.at(0)->startRefPos - SHORT_VAR_ALN_CHECK_EXTEND_SIZE;
	end_var_pos_pure = varVec.at(varVec.size()-1)->endRefPos + SHORT_VAR_ALN_CHECK_EXTEND_SIZE;
	if(start_var_pos_pure<1) start_var_pos_pure = 1;
	if(end_var_pos_pure>chrlen_tmp) end_var_pos_pure = chrlen_tmp;
//	start_var_pos = start_var_pos_pure - EXT_SIZE_CHK_VAR_LOC_SMALL;
//	end_var_pos = end_var_pos_pure + EXT_SIZE_CHK_VAR_LOC_SMALL;  // deleted on 2025-04-03
	start_var_pos = start_var_pos_pure - EXT_SIZE_CHK_VAR_LOC;
	end_var_pos = end_var_pos_pure + EXT_SIZE_CHK_VAR_LOC;
	if(start_var_pos<1) start_var_pos = 1;
	if(end_var_pos>chrlen_tmp) end_var_pos = chrlen_tmp;

	reg_str = chrname + ":" + to_string(start_var_pos) + "-" + to_string(end_var_pos);
	pthread_mutex_lock(&mutex_fai);
	p_seq = fai_fetch(fai, reg_str.c_str(), &seq_len);
	pthread_mutex_unlock(&mutex_fai);
	refseq = p_seq;
	free(p_seq);

	for(i=0; i<minimap2_aln_vec.size(); i++){
		minimap2_aln = minimap2_aln_vec.at(i);
		if(minimap2_aln->type.compare("S")==0) continue;
		if(find(recovered_qid_vec.begin(), recovered_qid_vec.end(), minimap2_aln->query_id)!=recovered_qid_vec.end()) continue; // skip successfully recovered query

		startRefPos_cns = minimap2_aln->region_startRefPos + 1;
		endRefPos_cns = minimap2_aln->region_endRefPos;

		var_vec_sub.clear();
		cand_paf_alnseg_vec_sub.clear();
		cand_mm2_aln_idx_vec_sub.clear();
		pre_large_match_flag = false;
		chrname_tmp = minimap2_aln->chrname;
		qname_vec = minimap2_aln->qname;

		// load the clipping data
		clipAlnDataLoader data_loader(chrname_tmp, startRefPos_cns, endRefPos_cns, inBamFile, minClipEndSize, maxVarRegSize, minMapQ, minHighMapQ);
		data_loader.loadClipAlnDataWithSATag(clipAlnDataVector, max_ultra_high_cov, qname_vec);

		if(minimap2_aln->pafalnsegs.size()>0){
			for(j=0; j<minimap2_aln->pafalnsegs.size(); j++){
				paf_alnseg = minimap2_aln->pafalnsegs.at(j);
				op = paf_alnseg->opflag;

				if(op==BAM_CMATCH or op==BAM_CEQUAL) {
					if(paf_alnseg->seglen>2*MIN_VALID_MINIMAP2_SEG_SIZE) pre_large_match_flag = true;
					continue;
				}

				if(pre_large_match_flag==false and paf_alnseg->seglen<MIN_VALID_MINIMAP2_SEG_SIZE) continue; // skip segments having no previous large match segments

				flag = false;
				if(paf_alnseg->seglen>=min_sv_size){
					end_ref_pos = getEndRefPosAlnSeg(paf_alnseg->startRpos, paf_alnseg->opflag, paf_alnseg->seglen);
					if(isOverlappedPos(paf_alnseg->startRpos, end_ref_pos, start_var_pos, end_var_pos)){
						if(isOverlappedPos(paf_alnseg->startRpos, end_ref_pos, start_var_pos_pure, end_var_pos_pure)==false){

							if(paf_alnseg->opflag==BAM_CINS) flag = true; // INS
							else{ // DEL
								overlap_size = getOverlapSize(paf_alnseg->startRpos, end_ref_pos, start_var_pos, end_var_pos);
								ratio = (double)overlap_size / paf_alnseg->seglen;
								if(ratio>=0.1) flag = true;
								else{
									if(end_ref_pos<start_var_pos_pure){ // left side
										sv_size_sum = 0; dist = -1;
										for(k=j+1; k<(int32_t)minimap2_aln->pafalnsegs.size(); k++){
											paf_alnseg2 = minimap2_aln->pafalnsegs.at(k);
											if(paf_alnseg2->opflag==paf_alnseg->opflag and paf_alnseg2->startRpos<=end_var_pos and paf_alnseg2->seglen>=min_sv_size){
												sv_size_sum += paf_alnseg2->seglen;
												if(dist==-1) dist = paf_alnseg2->startRpos - end_ref_pos;
											}
											if(paf_alnseg2->startRpos>end_var_pos) break;
										}
										if(dist!=-1){
											ratio = (double)dist / (paf_alnseg->seglen + dist + sv_size_sum);
											if(ratio<0.2) flag = true;
										}
									}else if(paf_alnseg->startRpos>start_var_pos_pure){ // right side
										sv_size_sum = 0; dist = -1;
										for(k=(int32_t)j-1; k>=0; k--){
											paf_alnseg2 = minimap2_aln->pafalnsegs.at(k);
											end_ref_pos2 = getEndRefPosAlnSeg(paf_alnseg2->startRpos, paf_alnseg2->opflag, paf_alnseg2->seglen);
											if(paf_alnseg2->opflag==paf_alnseg->opflag and end_ref_pos2>=start_var_pos and paf_alnseg2->seglen>=min_sv_size){
												sv_size_sum += paf_alnseg2->seglen;
												if(dist==-1) dist = paf_alnseg->startRpos - end_ref_pos2;
											}
											if(end_ref_pos2<start_var_pos) break;
										}
										if(dist!=-1){
											ratio = (double)dist / (paf_alnseg->seglen + dist + sv_size_sum);
											if(ratio<0.2) flag = true;
										}
									}
									//cout << "overlap_size=" << overlap_size << ", ratio=" << ratio << ", opflag=" << paf_alnseg->opflag << ", seglen=" << paf_alnseg->seglen << ", " << alnfilename << endl;
								}
							}
						}else flag = true;
					}else if(var_vec_sub.size()>0){ // try neighbouring items
						reg = var_vec_sub.at(var_vec_sub.size()-1);
						if(minimap2_aln->minimap2_aln_id==reg->minimap2_aln_id){
							dist = abs(paf_alnseg->startRpos - reg->startRefPos);
							if(dist<=min_distance_merge) flag = true;
						}
					}
				}
				if(flag){
					switch (op) {
						case BAM_CMATCH:
							break;
						case BAM_CINS:
							//search
//								cout << __func__ << ", line=" << __LINE__ << "minMapQ :  " << minMapQ << endl;
							supp_num_vec = computeSuppNumFromRegionAlnSegs(qname_vec, paf_alnseg, clipAlnDataVector, chrname_tmp, startRefPos_cns, endRefPos_cns, QC_SIZE_RATIO_MATCH_THRES_INDEL, QC_SIZE_RATIO_MATCH_THRES_INDEL_MERGE);
							//if(supp_num_vec.at(0)>=minReadsNumSupportSV){ // deleted on 2025-05-16
							if(supp_num_vec.at(0)>=minReadsNumSupportSV or (supp_num_vec.at(1)>0 and supp_num_vec.at(0)>=supp_num_vec.at(1)*READS_NUM_SUPPORT_FACTOR_LOW_COV) or
									(paf_alnseg->seglen>=2*min_sv_size and ((supp_num_vec.at(0)>=MIN_SUPPORT_READS_NUM_LEVEL_1 and supp_num_vec.at(1)<=COV_DEPTH_LEVEL_1) or (supp_num_vec.at(0)>=MIN_SUPPORT_READS_NUM_LEVEL_2 and supp_num_vec.at(1)<=COV_DEPTH_LEVEL_2)))){
								//write to newVarVec
								reg = new reg_t();
								reg->var_type = VAR_INS;
								reg->chrname = chrname_tmp;
								reg->startRefPos = paf_alnseg->startRpos + 1;
								reg->endRefPos = paf_alnseg->startRpos + 1;
								reg->startLocalRefPos = paf_alnseg->startSubjectPos + 1;
								reg->endLocalRefPos = paf_alnseg->startSubjectPos + 1;
								reg->startQueryPos = paf_alnseg->startQpos + 1;
								reg->endQueryPos = paf_alnseg->startQpos + paf_alnseg->seglen;
								reg->query_id = minimap2_aln->query_id;
								reg->sv_len = paf_alnseg->seglen;
								reg->minimap2_aln_id = minimap2_aln->minimap2_aln_id;
								reg->call_success_status = true;
								reg->short_sv_flag = false;
								reg->zero_cov_flag = false;
								reg->aln_seg_end_flag = false;
								reg->query_pos_invalid_flag = false;
								reg->large_indel_flag = false;
								reg->merge_flag = false;
								reg->rescue_flag = rescue_flag;
								reg->aln_orient = minimap2_aln->relative_strand;
								reg->gt_type = -1;
								reg->gt_seq = "";
								reg->refseq =  paf_alnseg->ref_seq;
								reg->altseq = paf_alnseg->alt_seq;
								reg->AF = 0;
								reg->supp_num = supp_num_vec.at(0);
								reg->DP = supp_num_vec.at(1);
								reg->discover_level = VAR_DISCOV_L_CNS_ALN;
								//svPosCorrection(reg);
								var_vec_sub.push_back(reg);
							}else{
								cand_mm2_aln_idx_vec_sub.push_back(minimap2_aln->minimap2_aln_id);
								cand_paf_alnseg_vec_sub.push_back(paf_alnseg);
//									cout << "Nsupp=" << supp_num_vec.at(0) << ", " << chrname << ";" << startRefPos_cns << "-" << endRefPos_cns << ", op=INS, oplen=" << paf_alnseg->seglen << endl;
							}
							break;
						case BAM_CDEL:
//							cout << __func__ << ", line=" << __LINE__ << "minMapQ :  " << minMapQ << endl;
							supp_num_vec = computeSuppNumFromRegionAlnSegs(qname_vec, paf_alnseg, clipAlnDataVector, chrname_tmp, startRefPos_cns, endRefPos_cns, QC_SIZE_RATIO_MATCH_THRES_INDEL, QC_SIZE_RATIO_MATCH_THRES_INDEL_MERGE);
							//if(supp_num_vec.at(0)>=minReadsNumSupportSV){ // deleted on 2025-05-16
							if(supp_num_vec.at(0)>=minReadsNumSupportSV or (supp_num_vec.at(1)>0 and supp_num_vec.at(0)>=supp_num_vec.at(1)*READS_NUM_SUPPORT_FACTOR_LOW_COV) or
									(paf_alnseg->seglen>=2*min_sv_size and ((supp_num_vec.at(0)>=MIN_SUPPORT_READS_NUM_LEVEL_1 and supp_num_vec.at(1)<=COV_DEPTH_LEVEL_1) or (supp_num_vec.at(0)>=MIN_SUPPORT_READS_NUM_LEVEL_2 and supp_num_vec.at(1)<=COV_DEPTH_LEVEL_2)))){
								//write to newVarVec
								reg = new reg_t();
								reg->var_type = VAR_DEL;
								reg->chrname = chrname_tmp;
								reg->startRefPos = paf_alnseg->startRpos + 1;
								reg->endRefPos = paf_alnseg->startRpos + paf_alnseg->seglen;
								reg->startLocalRefPos = paf_alnseg->startSubjectPos + 1;
								reg->endLocalRefPos = paf_alnseg->startSubjectPos + paf_alnseg->seglen;
								reg->startQueryPos = paf_alnseg->startQpos + 1;
								reg->endQueryPos = paf_alnseg->startQpos + 1;
								reg->query_id = minimap2_aln->query_id;
								reg->sv_len = -paf_alnseg->seglen;
								reg->minimap2_aln_id = minimap2_aln->minimap2_aln_id;
								reg->call_success_status = true;
								reg->short_sv_flag = false;
								reg->zero_cov_flag = false;
								reg->aln_seg_end_flag = false;
								reg->query_pos_invalid_flag = false;
								reg->large_indel_flag = false;
								reg->merge_flag = false;
								reg->rescue_flag = rescue_flag;
								reg->aln_orient = minimap2_aln->relative_strand;
								reg->gt_type = -1;
								reg->gt_seq = "";
								reg->refseq =  paf_alnseg->ref_seq;
								reg->altseq = paf_alnseg->alt_seq;
								reg->AF = 0;
								reg->supp_num = supp_num_vec.at(0);
								reg->DP = supp_num_vec.at(1);
								reg->discover_level = VAR_DISCOV_L_CNS_ALN;
								//svPosCorrection(reg);
								var_vec_sub.push_back(reg);
							}else{
								cand_mm2_aln_idx_vec_sub.push_back(minimap2_aln->minimap2_aln_id);
								cand_paf_alnseg_vec_sub.push_back(paf_alnseg);
//									cout << "Nsupp=" << supp_num_vec.at(0) << ", " << chrname << ";" << startRefPos_cns << "-" << endRefPos_cns << ", op=DEL, oplen=" << paf_alnseg->seglen << endl;
							}
							break;
						case BAM_CSOFT_CLIP:
						case BAM_CHARD_CLIP:
						case BAM_CEQUAL:
						case BAM_CDIFF:
							break;
						case BAM_CREF_SKIP:
							// unexpected events
						case BAM_CPAD:
						case BAM_CBACK:
						default:
							cerr << __func__ << ", line=" << __LINE__ << ": invalid opflag " << op << endl;
							exit(1);
					}
				}
			}

			// recover candidate items for noisy regions
			//if(cand_mm2_aln_idx_vec_sub.size()>0 or var_vec_sub.empty())
			if(find(recovered_qid_vec.begin(), recovered_qid_vec.end(), minimap2_aln->query_id)==recovered_qid_vec.end()) // each query only recover once
				recoverCandVarItems(var_vec_sub, varVec, cand_mm2_aln_idx_vec_sub, cand_paf_alnseg_vec_sub, recovered_qid_vec, minimap2_aln->minimap2_aln_id, minimap2_aln_vec, refseq, start_var_pos, end_var_pos, start_var_pos_pure, end_var_pos_pure, qname_vec, clipAlnDataVector, ctgfilename_para, rescue_flag);
			for(auto const& item : var_vec_sub) var_vec.push_back(item);
			for(auto const& item : cand_mm2_aln_idx_vec_sub) cand_mm2_aln_idx_vec.push_back(item);
			for(auto const& item : cand_paf_alnseg_vec_sub) cand_paf_alnseg_vec.push_back(item);
		}
		if(!clipAlnDataVector.empty()) destoryClipAlnData(clipAlnDataVector);
	}

	// process missing intra-aligned segments
	computeIntraMisPafAlnSeg(var_vec, varVec, cand_mm2_aln_idx_vec, cand_paf_alnseg_vec, minimap2_aln_vec, ctgfilename_para, rescue_flag);

	// process the split inter-aligned segments
	computeIndelFromSplitSegs(var_vec, minimap2_aln_vec, ctgfilename_para, startRefPos_cns, endRefPos_cns, QC_SIZE_RATIO_MATCH_THRES_INDEL, QC_SIZE_RATIO_MATCH_THRES_INDEL_MERGE, rescue_flag);

	if(var_vec.size()>0)
		mergeNeighbouringVars(var_vec, MAX_DIST_MERGE_ARBITARY, min_distance_merge, min_seqsim_merge, MIN_VALID_SIG_SIZE_RATIO_THRES, fai, ctgfilename_para, refseqfilename_para);

	var_vec.shrink_to_fit();
	
	if(!clipAlnDataVector.empty()) destoryClipAlnData(clipAlnDataVector);

	// delete minimap2_aln_vec
	if(rescue_flag) destroyMinimap2AlnVec(minimap2_aln_vec);

	return var_vec;
}

// recover candidate items for noisy regions
void varCand::recoverCandVarItems(vector<reg_t*> &var_vec_sub, vector<reg_t*> &varVec, vector<int32_t> &cand_mm2_aln_idx_vec, vector<struct pafalnSeg*> &cand_paf_alnseg_vec, vector<int32_t> &recovered_qid_vec, int32_t mm2_aln_idx, vector<minimap2_aln_t*> &minimap2_aln_vec, string &refseq, int64_t start_var_pos, int64_t end_var_pos, int64_t start_var_pos_pure, int64_t end_var_pos_pure, vector<string> &qname_vec, vector<clipAlnData_t*> &clipAlnDataVector, string &ctgfilename_para, bool rescue_flag){
	reg_t *reg_new;
	size_t i;
	int32_t max_merge_span, sdev, item_num, ins_sum_mean, del_sum_mean, ins_sum_num, del_sum_num, min_idx, target_var_op;
	int64_t aln_ins_sum, aln_del_sum, cand_ins_sum, cand_del_sum, aln_cand_ins_sum, aln_cand_del_sum, end_ref_pos, dist, min_dist;
	vector<clipAlnData_t*> clipAlnDataVector_sub;
	vector<struct querySeqInfoNode*> query_seq_info_all;
	struct querySeqInfoNode *qnode;
	vector<struct alnSeg*> query_alnSegs;
	struct alnSeg *seg;
	vector<int32_t> cand_mm2_aln_idx_vec_sub;
	vector<struct pafalnSeg*> cand_paf_alnseg_vec_sub;
	double size_ratio, dif, mean_size;
	bool valid_flag, success_recover_flag;
#if RECOVER_OUTPUT_DEBUG
	bool recover_flag;
#endif
	if(qname_vec.size()<(size_t)minReadsNumSupportSV) return;

	max_merge_span = min_distance_merge;
	success_recover_flag = false;

	aln_ins_sum = aln_del_sum = 0;
	for(auto const& reg : var_vec_sub){
		if(reg->var_type==VAR_INS) aln_ins_sum += reg->sv_len;
		else if(reg->var_type==VAR_DEL) aln_del_sum += abs(reg->sv_len);
	}
	cand_ins_sum = cand_del_sum = 0;
	for(auto const& pafseg : cand_paf_alnseg_vec){
		if(pafseg->opflag==BAM_CINS) cand_ins_sum += pafseg->seglen;
		else if(pafseg->opflag==BAM_CDEL) cand_del_sum += pafseg->seglen;
	}
	aln_cand_ins_sum = aln_ins_sum + cand_ins_sum;
	aln_cand_del_sum = aln_del_sum + cand_del_sum;

#if READS_CALL_DEBUG
	cout << "aln_ins_sum=" << aln_ins_sum << ", aln_del_sum=" << aln_del_sum << ", cand_ins_sum=" << cand_ins_sum << ", cand_del_sum=" << cand_del_sum << ", aln_cand_ins_sum=" << aln_cand_ins_sum << ", aln_cand_del_sum=" << aln_cand_del_sum << endl;
#endif

	// get the target group of align data
	for(auto const& item : clipAlnDataVector) if(find(qname_vec.begin(), qname_vec.end(), item->queryname)!=qname_vec.end()) clipAlnDataVector_sub.push_back(item);

	// extract queries from clip align data vector
	query_seq_info_all = extractQueriesFromClipAlnDataVec(clipAlnDataVector_sub, refseq, chrname, start_var_pos, end_var_pos, start_var_pos_pure, end_var_pos_pure, fai, minConReadLen, clip_reg_flag, &mutex_fai);

	for(auto const& qnode : query_seq_info_all){
//		if(qnode->qname.compare("SRR8858457.1.68744")==0){
//			cout << "qname=" << qnode->qname << endl;
//		}

		query_alnSegs = qnode->query_alnSegs;
		if(query_alnSegs.size()>0){
			seg = query_alnSegs.at(query_alnSegs.size()-1);
			end_ref_pos = getEndRefPosAlnSeg(seg->startRpos, seg->opflag, seg->seglen);
			if(query_alnSegs.at(0)->startRpos<=start_var_pos_pure+SHORT_VAR_ALN_CHECK_EXTEND_SIZE and end_ref_pos>=end_var_pos_pure-SHORT_VAR_ALN_CHECK_EXTEND_SIZE) qnode->entire_flanking_flag = true;
			else qnode->entire_flanking_flag = false;
			qnode->qcSig_vec = extractQcSigsFromAlnSegsSingleQuery(qnode, qnode->clip_aln->chrname, start_var_pos, end_var_pos, min_sv_size*2);
			// remove neighboring false signatures
			rmNeighborFalseSig(qnode->qcSig_vec, MIN_AVER_SIZE_ALN_SEG, MIN_SIZE_RATIO_MATCH_CLIP_POS);
			qnode->qcSig_merge_flag = mergeNeighbouringSigsFlag(qnode, qnode->qcSig_vec, MAX_DIST_MERGE_ARBITARY, max_merge_span, CALL_MIN_SEQSIM_MERGE_THRES2, MIN_VALID_SIG_SIZE_RATIO_THRES, fai);
			filterSmallSigs(qnode->qcSig_vec, MIN_SIG_SIZE_RATIO_FILTER_THRES, MAX_INNER_MISSING_IGNORE_SIZE);

			qnode->ins_sum = qnode->del_sum = 0;
			for(auto const& item : qnode->qcSig_vec){
				if(item->cigar_op==BAM_CINS) qnode->ins_sum += item->cigar_op_len;
				else qnode->del_sum += item->cigar_op_len;
			}
		}
	}

#if READS_CALL_DEBUG
	printQcSigIndel(query_seq_info_all);
#endif

	// compute mean size
	ins_sum_mean = del_sum_mean = ins_sum_num = del_sum_num = 0;
	for(auto const& qnode : query_seq_info_all){
		if(qnode->entire_flanking_flag){
			if(qnode->ins_sum>0){
				ins_sum_mean += qnode->ins_sum;
				ins_sum_num ++;
			}
			if(qnode->del_sum>0){
				del_sum_mean += qnode->del_sum;
				del_sum_num ++;
			}
		}
	}
	if(ins_sum_num>0) ins_sum_mean /= ins_sum_num;
	if(del_sum_num>0) del_sum_mean /= del_sum_num;

	if(ins_sum_mean>del_sum_mean) { target_var_op = BAM_CINS; mean_size = ins_sum_mean; }
	else { target_var_op = BAM_CDEL; mean_size = del_sum_mean; }

	if(mean_size>0){
		// compute sdev
		item_num = sdev = 0;
		for(auto const& qnode : query_seq_info_all){
			if(qnode->entire_flanking_flag){
				if(target_var_op==BAM_CINS and qnode->ins_sum>0){
					sdev += (qnode->ins_sum - mean_size) * (qnode->ins_sum - mean_size);
					item_num ++;
				}else if(target_var_op==BAM_CDEL and qnode->del_sum>0){
					sdev += (qnode->del_sum - del_sum_mean) * (qnode->del_sum - del_sum_mean);
					item_num ++;
				}
			}
		}
		sdev = sqrt(sdev/item_num);

		// remove the outliers
		if(item_num>0){
			for(i=0; i<query_seq_info_all.size(); ){
				qnode = query_seq_info_all.at(i);
				valid_flag = true;
				if(qnode->entire_flanking_flag){
					if(target_var_op==BAM_CINS and qnode->ins_sum>0){
						dif = abs(qnode->ins_sum - mean_size);
#if READS_CALL_DEBUG
						cout << "\tqname=" << qnode->qname << ", dif=" << dif << endl;
#endif
						if(dif>2*sdev and dif>min_sv_size){
#if READS_CALL_DEBUG
							cout << "\txxxxxxxxx dif=" << dif << ", qname=" << qnode->qname << endl;
#endif
							valid_flag = false;
						}
					}else if(target_var_op==BAM_CDEL and qnode->del_sum>0){
						dif = abs(qnode->del_sum - mean_size);
#if READS_CALL_DEBUG
						cout << "\tqname=" << qnode->qname << ", dif=" << dif << endl;
#endif
						if(dif>2*sdev and dif>min_sv_size){
#if READS_CALL_DEBUG
							cout << "\txxxxxxxxx dif=" << dif << ", qname=" << qnode->qname << endl;
#endif
							valid_flag = false;
						}
					}
				}
				if(valid_flag==false){
					destroyAlnSegs(qnode->query_alnSegs);
					destroyQueryQcSig(qnode->qcSig_vec);
					delete qnode;
					query_seq_info_all.erase(query_seq_info_all.begin()+i);
				}else i++;
			}
		}
	}

	// recompute mean size
	ins_sum_mean = del_sum_mean = ins_sum_num = del_sum_num = 0;
	for(auto const& qnode : query_seq_info_all){
		if(qnode->entire_flanking_flag){
			if(qnode->ins_sum>0){
				ins_sum_mean += qnode->ins_sum;
				ins_sum_num ++;
			}
			if(qnode->del_sum>0){
				del_sum_mean += qnode->del_sum;
				del_sum_num ++;
			}
		}
	}
	if(ins_sum_num>0) ins_sum_mean /= ins_sum_num;
	if(del_sum_num>0) del_sum_mean /= del_sum_num;

	if(ins_sum_mean>del_sum_mean and ins_sum_mean>=EXT_SIZE_CHK_VAR_LOC_SMALL and ins_sum_num>=minReadsNumSupportSV){ // INS
		if(aln_ins_sum<=ins_sum_mean) size_ratio = (double)aln_ins_sum / ins_sum_mean;
		else size_ratio = (double)ins_sum_mean / aln_ins_sum;
		if(size_ratio>=MIN_ALN_SIZE_RATIO_RECOVER){ // discard candidate items
			cand_paf_alnseg_vec.clear();
			cand_mm2_aln_idx_vec.clear();
			success_recover_flag = true;
		}else{
			if(aln_cand_ins_sum<=ins_sum_mean) size_ratio = (double)aln_cand_ins_sum / ins_sum_mean;
			else size_ratio = (double)ins_sum_mean / aln_cand_ins_sum;
			if(size_ratio<MIN_SIZE_RATIO_RECOVER){ // incomplete

#if READS_CALL_DEBUG
				cout << "size_ratio=" << size_ratio << ", incomplete." << endl;
#endif

				// get representative idx
				min_idx = -1;
				min_dist = INT_MAX;
				for(i=0; i<query_seq_info_all.size(); i++){
					qnode = query_seq_info_all.at(i);
					if(qnode->entire_flanking_flag and qnode->ins_sum>0 and qnode->qcSig_vec.size()==1){
						dist = abs(qnode->ins_sum-ins_sum_mean);
						if(dist<min_dist){
							min_idx = i;
							min_dist = dist;
						}
					}
				}
				if(min_idx==-1){
					for(i=0; i<query_seq_info_all.size(); i++){
						qnode = query_seq_info_all.at(i);
						if(qnode->entire_flanking_flag and qnode->ins_sum>0){
							dist = abs(qnode->ins_sum-ins_sum_mean);
							if(dist<min_dist){
								min_idx = i;
								min_dist = dist;
							}
						}
					}
				}

#if READS_CALL_DEBUG
				cout << "min_idx=" << min_idx << ", min_dist=" << min_dist << endl;
#endif

				if(min_idx!=-1 and query_seq_info_all.size()>=(size_t)minReadsNumSupportSV){
					// remove original items
					for(auto const& item : var_vec_sub) delete item;
					var_vec_sub.clear();

					// add new items
#if RECOVER_OUTPUT_DEBUG
					recover_flag = false;
#endif
					qnode = query_seq_info_all.at(min_idx);
					for(auto const& sig : qnode->qcSig_vec){
						if(sig->cigar_op==BAM_CINS){
							//write to newVarVec
							reg_new = new reg_t();
							reg_new->var_type = VAR_INS;
							reg_new->chrname = chrname;
							reg_new->startRefPos = sig->ref_pos + 1;
							reg_new->endRefPos = getEndRefPosAlnSeg(sig->ref_pos, sig->cigar_op, sig->cigar_op_len) + 1;
							reg_new->startLocalRefPos = -1;
							reg_new->endLocalRefPos = -1;
							reg_new->startQueryPos = -1;
							reg_new->endQueryPos = -1;
							reg_new->query_id = minimap2_aln_vec.at(mm2_aln_idx)->query_id;
							reg_new->sv_len = sig->cigar_op_len;
							reg_new->minimap2_aln_id = mm2_aln_idx;
							reg_new->call_success_status = true;
							reg_new->short_sv_flag = false;
							reg_new->zero_cov_flag = false;
							reg_new->aln_seg_end_flag = false;
							reg_new->query_pos_invalid_flag = false;
							reg_new->large_indel_flag = false;
							reg_new->merge_flag = false;
							reg_new->rescue_flag = rescue_flag;
							reg_new->aln_orient = minimap2_aln_vec.at(mm2_aln_idx)->relative_strand;
							reg_new->gt_type = -1;
							reg_new->gt_seq = "";
							reg_new->refseq =  sig->refseq;
							reg_new->altseq = sig->altseq;
							reg_new->AF = 0;
							reg_new->supp_num = qname_vec.size();
							reg_new->DP = computeCovNumReg(reg_new->chrname, reg_new->startRefPos, reg_new->endRefPos, fai, inBamFile, minMapQ, minHighMapQ, max_ultra_high_cov);
							reg_new->discover_level = VAR_DISCOV_L_READS;

							var_vec_sub.push_back(reg_new);
							success_recover_flag = true;
#if RECOVER_OUTPUT_DEBUG
							recover_flag = true;
#endif
						}
					}
					cand_paf_alnseg_vec.clear();
					cand_mm2_aln_idx_vec.clear();

#if RECOVER_OUTPUT_DEBUG
					if(recover_flag)
						cout << "################ line=" << __LINE__ << ", INS recover: " << alnfilename << endl;
#endif
				}
			}else{ // valid
				if(query_seq_info_all.size()>=(size_t)minReadsNumSupportSV){
#if RECOVER_OUTPUT_DEBUG
					recover_flag = false;
#endif
					for(auto const& pafseg : cand_paf_alnseg_vec){
						if(pafseg->opflag==BAM_CINS){
							//write to newVarVec
							reg_new = new reg_t();
							reg_new->var_type = VAR_INS;
							reg_new->chrname = chrname;
							reg_new->startRefPos = pafseg->startRpos + 1;
							reg_new->endRefPos = getEndRefPosAlnSeg(pafseg->startRpos, pafseg->opflag, pafseg->seglen) + 1;
							reg_new->startLocalRefPos = pafseg->startSubjectPos + 1;
							reg_new->endLocalRefPos = getEndSubjectPosAlnSeg(pafseg->startSubjectPos, pafseg->opflag, pafseg->seglen) + 1;
							reg_new->startQueryPos = pafseg->startQpos + 1;
							reg_new->endQueryPos = getEndQueryPosAlnSeg(pafseg->startQpos, pafseg->opflag, pafseg->seglen) + 1;
							reg_new->query_id = minimap2_aln_vec.at(mm2_aln_idx)->query_id;
							reg_new->sv_len = pafseg->seglen;
							reg_new->minimap2_aln_id = mm2_aln_idx;
							reg_new->call_success_status = true;
							reg_new->short_sv_flag = false;
							reg_new->zero_cov_flag = false;
							reg_new->aln_seg_end_flag = false;
							reg_new->query_pos_invalid_flag = false;
							reg_new->large_indel_flag = false;
							reg_new->merge_flag = false;
							reg_new->rescue_flag = rescue_flag;
							reg_new->aln_orient = minimap2_aln_vec.at(mm2_aln_idx)->relative_strand;
							reg_new->gt_type = -1;
							reg_new->gt_seq = "";
							reg_new->refseq =  pafseg->ref_seq;
							reg_new->altseq = pafseg->alt_seq;
							reg_new->AF = 0;
							reg_new->supp_num = qname_vec.size();
							reg_new->DP = computeCovNumReg(reg_new->chrname, reg_new->startRefPos, reg_new->endRefPos, fai, inBamFile, minMapQ, minHighMapQ, max_ultra_high_cov);
							reg_new->discover_level = VAR_DISCOV_L_CNS_ALN;

							var_vec_sub.push_back(reg_new);
							success_recover_flag = true;

#if RECOVER_OUTPUT_DEBUG
							recover_flag = true;
#endif
						}
					}
					cand_paf_alnseg_vec.clear();
					cand_mm2_aln_idx_vec.clear();

					sortRegVec(var_vec_sub);

#if RECOVER_OUTPUT_DEBUG
					if(recover_flag)
						cout << "################ line=" << __LINE__ << ", INS recover: " << alnfilename << endl;
#endif
				}
			}
		}
	}else if(ins_sum_mean<del_sum_mean and del_sum_mean>=EXT_SIZE_CHK_VAR_LOC_SMALL and del_sum_num>=minReadsNumSupportSV){ // DEL
		if(aln_del_sum<=del_sum_mean) size_ratio = (double)aln_del_sum / del_sum_mean;
		else size_ratio = (double)del_sum_mean / aln_del_sum;
		if(size_ratio>=MIN_ALN_SIZE_RATIO_RECOVER){ // discard candidate items
			cand_paf_alnseg_vec.clear();
			cand_mm2_aln_idx_vec.clear();
			success_recover_flag = true;
		}else{
			if(aln_cand_del_sum<=del_sum_mean) size_ratio = (double)aln_cand_del_sum / del_sum_mean;
			else size_ratio = (double)del_sum_mean / aln_cand_del_sum;
			if(size_ratio<MIN_SIZE_RATIO_RECOVER){ // incomplete

#if READS_CALL_DEBUG
				cout << "size_ratio=" << size_ratio << ", incomplete." << endl;
#endif

				// get representative idx
				min_idx = -1;
				min_dist = INT_MAX;
				for(i=0; i<query_seq_info_all.size(); i++){
					qnode = query_seq_info_all.at(i);
					if(qnode->entire_flanking_flag and qnode->del_sum>0 and qnode->qcSig_vec.size()==1){
						dist = abs(qnode->del_sum-del_sum_mean);
						if(dist<min_dist){
							min_idx = i;
							min_dist = dist;
						}
					}
				}
				if(min_idx==-1){
					for(i=0; i<query_seq_info_all.size(); i++){
						qnode = query_seq_info_all.at(i);
						if(qnode->entire_flanking_flag and qnode->del_sum>0){
							dist = abs(qnode->del_sum-del_sum_mean);
							if(dist<min_dist){
								min_idx = i;
								min_dist = dist;
							}
						}
					}
				}

#if READS_CALL_DEBUG
				cout << "min_idx=" << min_idx << ", min_dist=" << min_dist << endl;
#endif

				if(min_idx!=-1 and query_seq_info_all.size()>=(size_t)minReadsNumSupportSV){
					// remove original items
					for(auto const& item : var_vec_sub) delete item;
					var_vec_sub.clear();

					// add new items
#if RECOVER_OUTPUT_DEBUG
					recover_flag = false;
#endif
					qnode = query_seq_info_all.at(min_idx);
					for(auto const& sig : qnode->qcSig_vec){
						if(sig->cigar_op==BAM_CDEL){
							//write to newVarVec
							reg_new = new reg_t();
							reg_new->var_type = VAR_DEL;
							reg_new->chrname = chrname;
							reg_new->startRefPos = sig->ref_pos + 1;
							reg_new->endRefPos = getEndRefPosAlnSeg(sig->ref_pos, sig->cigar_op, sig->cigar_op_len) + 1;
							reg_new->startLocalRefPos = -1;
							reg_new->endLocalRefPos = -1;
							reg_new->startQueryPos = -1;
							reg_new->endQueryPos = -1;
							reg_new->query_id = minimap2_aln_vec.at(mm2_aln_idx)->query_id;
							reg_new->sv_len = sig->cigar_op_len;
							reg_new->minimap2_aln_id = mm2_aln_idx;
							reg_new->call_success_status = true;
							reg_new->short_sv_flag = false;
							reg_new->zero_cov_flag = false;
							reg_new->aln_seg_end_flag = false;
							reg_new->query_pos_invalid_flag = false;
							reg_new->large_indel_flag = false;
							reg_new->merge_flag = false;
							reg_new->rescue_flag = rescue_flag;
							reg_new->aln_orient = -1;
							reg_new->gt_type = minimap2_aln_vec.at(mm2_aln_idx)->relative_strand;
							reg_new->gt_seq = "";
							reg_new->refseq =  sig->refseq;
							reg_new->altseq = sig->altseq;
							reg_new->AF = 0;
							reg_new->supp_num = qname_vec.size();
							reg_new->DP = computeCovNumReg(reg_new->chrname, reg_new->startRefPos, reg_new->endRefPos, fai, inBamFile, minMapQ, minHighMapQ, max_ultra_high_cov);
							reg_new->discover_level = VAR_DISCOV_L_READS;

							var_vec_sub.push_back(reg_new);
							success_recover_flag = true;

#if RECOVER_OUTPUT_DEBUG
							recover_flag = true;
#endif
						}
					}
					cand_paf_alnseg_vec.clear();
					cand_mm2_aln_idx_vec.clear();

#if RECOVER_OUTPUT_DEBUG
					if(recover_flag)
						cout << "################ line=" << __LINE__ << ", DEL recover: " << alnfilename << endl;
#endif
				}
			}else{ //valid
				if(query_seq_info_all.size()>=(size_t)minReadsNumSupportSV){
#if RECOVER_OUTPUT_DEBUG
					recover_flag = false;
#endif
					for(auto const& pafseg : cand_paf_alnseg_vec){
						if(pafseg->opflag==BAM_CDEL){
							//write to newVarVec
							reg_new = new reg_t();
							reg_new->var_type = VAR_DEL;
							reg_new->chrname = chrname;
							reg_new->startRefPos = pafseg->startRpos + 1;
							reg_new->endRefPos = getEndRefPosAlnSeg(pafseg->startRpos, pafseg->opflag, pafseg->seglen) + 1;
							reg_new->startLocalRefPos = pafseg->startSubjectPos + 1;
							reg_new->endLocalRefPos = getEndSubjectPosAlnSeg(pafseg->startSubjectPos, pafseg->opflag, pafseg->seglen) + 1;
							reg_new->startQueryPos = pafseg->startQpos + 1;
							reg_new->endQueryPos = getEndQueryPosAlnSeg(pafseg->startQpos, pafseg->opflag, pafseg->seglen) + 1;
							reg_new->query_id = minimap2_aln_vec.at(mm2_aln_idx)->query_id;
							reg_new->sv_len = pafseg->seglen;
							reg_new->minimap2_aln_id = mm2_aln_idx;
							reg_new->call_success_status = true;
							reg_new->short_sv_flag = false;
							reg_new->zero_cov_flag = false;
							reg_new->aln_seg_end_flag = false;
							reg_new->query_pos_invalid_flag = false;
							reg_new->large_indel_flag = false;
							reg_new->merge_flag = false;
							reg_new->rescue_flag = rescue_flag;
							reg_new->aln_orient = minimap2_aln_vec.at(mm2_aln_idx)->relative_strand;
							reg_new->gt_type = -1;
							reg_new->gt_seq = "";
							reg_new->refseq =  pafseg->ref_seq;
							reg_new->altseq = pafseg->alt_seq;
							reg_new->AF = 0;
							reg_new->supp_num = qname_vec.size();
							reg_new->DP = computeCovNumReg(reg_new->chrname, reg_new->startRefPos, reg_new->endRefPos, fai, inBamFile, minMapQ, minHighMapQ, max_ultra_high_cov);
							reg_new->discover_level = VAR_DISCOV_L_CNS_ALN;

							var_vec_sub.push_back(reg_new);
							success_recover_flag = true;

#if RECOVER_OUTPUT_DEBUG
							recover_flag = true;
#endif
						}
					}
					cand_paf_alnseg_vec.clear();
					cand_mm2_aln_idx_vec.clear();

					sortRegVec(var_vec_sub);

#if RECOVER_OUTPUT_DEBUG
					if(recover_flag)
						cout << "################ line=" << __LINE__ << ", DEL recover: " << alnfilename << endl;
#endif
				}
			}
		}
	}

	if(success_recover_flag) recovered_qid_vec.push_back(minimap2_aln_vec.at(mm2_aln_idx)->query_id);

	if(!query_seq_info_all.empty()) destoryQuerySeqInfoAll(query_seq_info_all);
}

// process missing intra-aligned segments
void varCand::computeIntraMisPafAlnSeg(vector<reg_t*> &var_vec, vector<reg_t*> &varVec, vector<int32_t> &cand_mm2_aln_idx_vec, vector<struct pafalnSeg*> &cand_paf_alnseg_vec, vector<minimap2_aln_t*> &minimap2_aln_vec, string &ctgfilename_para, bool rescue_flag){
	size_t i, j;
	reg_t *reg, *reg_new, *target_reg;
	struct pafalnSeg* paf_alnseg, *paf_alnseg_tmp, *paf_alnseg_tmp2;
	int64_t end_ref_pos, dist, start_pos_paf_ext, end_pos_paf_ext, start_pos_paf, seq_len;
	int64_t start_ref_pos_new, end_ref_pos_new, start_localref_pos_new, end_localref_pos_new, start_query_pos_new, end_query_pos_new;
	int32_t k, minimap2_aln_id, start_var_idx, end_var_idx, start_paf_idx, end_paf_idx, var_len, normal_len;
	double size_ratio;
	bool change_flag, contain_flag;
	minimap2_aln_t *mm2_aln;
	string reg_str, refseq_tmp, queryseq_tmp;
	char *p_seq;
	vector<reg_t*> var_vec_contain;
	vector<struct pafalnSeg*> pafaln_cand_vec;
	int64_t start_ref_pos_tmp, end_ref_pos_tmp, sum_len, sum_len2, sum_total, dif;
	int32_t max_size;

	change_flag = false;

	for(i=0; i<cand_paf_alnseg_vec.size(); i++){
		paf_alnseg = cand_paf_alnseg_vec.at(i);
		minimap2_aln_id = cand_mm2_aln_idx_vec.at(i);

		if(var_vec.size()>0){
			for(j=0; j<var_vec.size(); j++){
				dist = INT_MAX;
				reg = var_vec.at(j);
				if(reg->minimap2_aln_id==minimap2_aln_id){
					if(reg->var_type==VAR_INS and paf_alnseg->opflag==BAM_CINS){ // INS
						dist = reg->startRefPos - (paf_alnseg->startRpos + 1);
						if(dist<0) dist = -dist;
					}else if(reg->var_type==VAR_DEL and paf_alnseg->opflag==BAM_CDEL){ // DEL
						end_ref_pos = getEndRefPosAlnSeg(paf_alnseg->startRpos, paf_alnseg->opflag, paf_alnseg->seglen);
						if(paf_alnseg->startRpos<reg->startRefPos){ // candidate on the left
							dist = reg->startRefPos - end_ref_pos;
						}else{ // candidate on the right
							dist = paf_alnseg->startRpos - reg->endRefPos;
						}
					}

					if(dist<VAR_ALN_EXTEND_SIZE){
						size_ratio = (double)dist / (paf_alnseg->seglen + abs(reg->sv_len) + dist);
						//if(size_ratio<MIN_VALID_SIG_SIZE_RATIO_THRES){
						if(size_ratio<MIN_VALID_SIG_SIZE_RATIO_THRES or dist<=MAX_DIST_MERGE_ARBITARY){
							//write to newVarVec
							reg_new = new reg_t();
							reg_new->var_type = reg->var_type;
							reg_new->chrname = reg->chrname;
							reg_new->startRefPos = paf_alnseg->startRpos + 1;
							reg_new->endRefPos = getEndRefPosAlnSeg(paf_alnseg->startRpos, paf_alnseg->opflag, paf_alnseg->seglen) + 1;
							reg_new->startLocalRefPos = paf_alnseg->startSubjectPos + 1;
							reg_new->endLocalRefPos = getEndSubjectPosAlnSeg(paf_alnseg->startSubjectPos, paf_alnseg->opflag, paf_alnseg->seglen) + 1;
							reg_new->startQueryPos = paf_alnseg->startQpos + 1;
							reg_new->endQueryPos = getEndQueryPosAlnSeg(paf_alnseg->startQpos, paf_alnseg->opflag, paf_alnseg->seglen) + 1;
							reg_new->query_id = reg->query_id;
							reg_new->sv_len = paf_alnseg->seglen;
							if(reg->var_type==VAR_DEL) reg_new->sv_len = -reg_new->sv_len;
							reg_new->minimap2_aln_id = reg->minimap2_aln_id;
							reg_new->call_success_status = true;
							reg_new->short_sv_flag = false;
							reg_new->zero_cov_flag = false;
							reg_new->aln_seg_end_flag = false;
							reg_new->query_pos_invalid_flag = false;
							reg_new->large_indel_flag = false;
							reg_new->merge_flag = false;
							reg_new->rescue_flag = rescue_flag;
							reg_new->aln_orient = reg->aln_orient;
							reg_new->gt_type = -1;
							reg_new->gt_seq = "";
							reg_new->refseq =  paf_alnseg->ref_seq;
							reg_new->altseq = paf_alnseg->alt_seq;
							reg_new->AF = 0;
							reg_new->supp_num = reg->supp_num;
							reg_new->DP = reg->DP;
							reg_new->discover_level = VAR_DISCOV_L_CNS_ALN;

							var_vec.push_back(reg_new);

							change_flag = true;
							break;
						}
					}
				}
			}
		}else{ // added on 2026-02-18
			mm2_aln = minimap2_aln_vec.at(minimap2_aln_id);
			start_pos_paf_ext = paf_alnseg->startRpos - EXT_SIZE_CHK_VAR_LOC_SMALL;
			end_ref_pos = getEndRefPosAlnSeg(paf_alnseg->startRpos, paf_alnseg->opflag, paf_alnseg->seglen);
			end_pos_paf_ext = end_ref_pos + EXT_SIZE_CHK_VAR_LOC_SMALL;
			if(start_pos_paf_ext<1) start_pos_paf_ext = 1;

			start_var_idx = end_var_idx = -1;
			for(j=0; j<varVec.size(); j++){
				reg = varVec.at(j);
				if((reg->var_type==VAR_INS and paf_alnseg->opflag==BAM_CINS) or (reg->var_type==VAR_DEL and paf_alnseg->opflag==BAM_CDEL)){
					if(start_var_idx==-1 and reg->startRefPos>=start_pos_paf_ext){
						start_var_idx = j;
						break;
					}
				}
			}
			for(k=(int32_t)varVec.size()-1; k>=0; k--){
				reg = varVec.at(k);
				if((reg->var_type==VAR_INS and paf_alnseg->opflag==BAM_CINS) or (reg->var_type==VAR_DEL and paf_alnseg->opflag==BAM_CDEL)){
					if(end_var_idx==-1 and reg->endRefPos<=end_pos_paf_ext){
						end_var_idx = k;
						break;
					}
				}
			}
			if(start_var_idx!=-1 and end_var_idx!=-1 and start_var_idx==end_var_idx){
				reg = varVec.at(start_var_idx);
				if(abs(reg->sv_len)>paf_alnseg->seglen){
					start_paf_idx = end_paf_idx = -1;
					for(j=0; j<mm2_aln->pafalnsegs.size(); j++){
						paf_alnseg_tmp = mm2_aln->pafalnsegs.at(j);
						if((reg->var_type==VAR_INS and paf_alnseg_tmp->opflag==BAM_CINS) or (reg->var_type==VAR_DEL and paf_alnseg_tmp->opflag==BAM_CDEL)){
							if(start_paf_idx==-1 and paf_alnseg_tmp->startRpos>=start_pos_paf_ext){
								start_paf_idx = j;
								break;
							}
						}
					}
					for(k=(int32_t)mm2_aln->pafalnsegs.size()-1; k>=0; k--){
						paf_alnseg_tmp = mm2_aln->pafalnsegs.at(k);
						if((reg->var_type==VAR_INS and paf_alnseg_tmp->opflag==BAM_CINS) or (reg->var_type==VAR_DEL and paf_alnseg_tmp->opflag==BAM_CDEL)){
							end_ref_pos = getEndRefPosAlnSeg(paf_alnseg_tmp->startRpos, paf_alnseg_tmp->opflag, paf_alnseg_tmp->seglen);
							if(end_paf_idx==-1 and end_ref_pos<=end_pos_paf_ext){
								end_paf_idx = k;
								break;
							}
						}
					}
					if(start_paf_idx!=-1 and end_paf_idx!=-1 and start_paf_idx<=end_paf_idx){
						contain_flag = false;
						for(k=start_paf_idx; k<=end_paf_idx; k++){
							paf_alnseg_tmp = mm2_aln->pafalnsegs.at(k);
							if(paf_alnseg_tmp==paf_alnseg){
								contain_flag = true;
								break;
							}
						}
						if(contain_flag){
							var_len = normal_len = 0;
							for(k=start_paf_idx; k<=end_paf_idx; k++){
								paf_alnseg_tmp = mm2_aln->pafalnsegs.at(k);
								if((reg->var_type==VAR_INS and paf_alnseg_tmp->opflag==BAM_CINS) or (reg->var_type==VAR_DEL and paf_alnseg_tmp->opflag==BAM_CDEL))
									var_len += paf_alnseg_tmp->seglen;
								else
									normal_len += paf_alnseg_tmp->seglen;
							}

							if(var_len<abs(reg->sv_len)) size_ratio = (double)var_len / abs(reg->sv_len);
							else size_ratio = (double)abs(reg->sv_len) / var_len;
							if(size_ratio>=CALL_MIN_SEQSIM_MERGE_THRES2){
								//adjust the location
								paf_alnseg_tmp = mm2_aln->pafalnsegs.at(start_paf_idx);
								paf_alnseg_tmp2 = mm2_aln->pafalnsegs.at(end_paf_idx);
								start_pos_paf = paf_alnseg_tmp->startRpos + 1;
								start_ref_pos_new = reg->startRefPos;
								end_ref_pos_new = reg->endRefPos;
								if(start_pos_paf<reg->startRefPos){
									dist = reg->startRefPos - start_pos_paf;
									start_localref_pos_new = paf_alnseg_tmp->startSubjectPos + 1 + dist;
									start_query_pos_new = paf_alnseg_tmp->startQpos + 1 + dist;
									end_localref_pos_new = getEndSubjectPosAlnSeg(paf_alnseg_tmp2->startSubjectPos, paf_alnseg_tmp2->opflag, paf_alnseg_tmp2->seglen) + 1 + dist - normal_len;
									end_query_pos_new = getEndQueryPosAlnSeg(paf_alnseg_tmp2->startQpos, paf_alnseg_tmp2->opflag, paf_alnseg_tmp2->seglen) + 1 + dist - normal_len;
								}else{
									dist = start_pos_paf - reg->startRefPos;
									start_localref_pos_new = paf_alnseg_tmp->startSubjectPos + 1 - dist;
									start_query_pos_new = paf_alnseg_tmp->startQpos + 1 - dist;
									end_localref_pos_new = getEndSubjectPosAlnSeg(paf_alnseg_tmp2->startSubjectPos, paf_alnseg_tmp2->opflag, paf_alnseg_tmp2->seglen) + 1 - dist - normal_len;
									end_query_pos_new = getEndQueryPosAlnSeg(paf_alnseg_tmp2->startQpos, paf_alnseg_tmp2->opflag, paf_alnseg_tmp2->seglen) + 1 - dist - normal_len;
								}
								if(start_query_pos_new<1) start_query_pos_new = 1;
								if(end_query_pos_new>mm2_aln->query_len) end_query_pos_new = mm2_aln->query_len;

								// get refseq and altseq
								reg_str = reg->chrname + ":" + to_string(start_ref_pos_new) + "-" + to_string(end_ref_pos_new);
								pthread_mutex_lock(&mutex_fai);
								p_seq = fai_fetch64(fai, reg_str.c_str(), &seq_len);
								pthread_mutex_unlock(&mutex_fai);
								refseq_tmp = p_seq;
								free(p_seq);

								FastaSeqLoader fa_loader(ctgfilename_para);
								queryseq_tmp = fa_loader.getFastaSeqByPos(mm2_aln->query_id, start_query_pos_new, end_query_pos_new, mm2_aln->relative_strand);

								//write to newVarVec
								reg_new = new reg_t();
								reg_new->var_type = reg->var_type;
								reg_new->chrname = reg->chrname;
								reg_new->startRefPos = start_ref_pos_new;
								reg_new->endRefPos = end_ref_pos_new;
								reg_new->startLocalRefPos = start_localref_pos_new;
								reg_new->endLocalRefPos = end_localref_pos_new;
								reg_new->startQueryPos = start_query_pos_new;
								reg_new->endQueryPos = end_query_pos_new;
								reg_new->query_id = mm2_aln->query_id;
								reg_new->sv_len = reg->sv_len;
								reg_new->minimap2_aln_id = mm2_aln->minimap2_aln_id;
								reg_new->call_success_status = true;
								reg_new->short_sv_flag = false;
								reg_new->zero_cov_flag = false;
								reg_new->aln_seg_end_flag = false;
								reg_new->query_pos_invalid_flag = false;
								reg_new->large_indel_flag = false;
								reg_new->merge_flag = false;
								reg_new->rescue_flag = rescue_flag;
								reg_new->aln_orient = reg->aln_orient;
								reg_new->gt_type = -1;
								reg_new->gt_seq = "";
								reg_new->refseq =  refseq_tmp;
								reg_new->altseq = queryseq_tmp;
								reg_new->AF = 0;
								reg_new->supp_num = reg->supp_num;
								reg_new->DP = reg->DP;
								reg_new->discover_level = VAR_DISCOV_L_CNS_ALN;

								var_vec.push_back(reg_new);

								change_flag = true;
								break;
							}
						}
					}
				}
			}
		}
	}

	// recover the split variant item around
	if(var_vec.size()>0){
		for(auto const& var: varVec){
			if(abs(var->sv_len)<MAX_INNER_MISSING_IGNORE_SIZE) continue;

			start_ref_pos_tmp = var->startRefPos - EXT_SIZE_CHK_VAR_LOC_SMALL;
			end_ref_pos_tmp = var->endRefPos + EXT_SIZE_CHK_VAR_LOC_SMALL;
			if(start_ref_pos_tmp<1) start_ref_pos_tmp = 1;

			var_vec_contain.clear();
			for(auto const& var2: var_vec){
				if(var2->var_type==var->var_type and isOverlappedPos(var2->startRefPos, var2->endRefPos, start_ref_pos_tmp, end_ref_pos_tmp)){
					var_vec_contain.push_back(var2);
				}
			}

			target_reg = NULL;
			max_size = 0;
			for(auto const& var2: var_vec_contain){
				if(abs(var2->sv_len)>max_size){
					max_size = abs(var2->sv_len);
					target_reg = var2;
				}
			}

			if(target_reg){
				sum_len = 0;
				for(i=0; i<var_vec_contain.size(); ){
					reg = var_vec_contain.at(i);
					if(reg->query_id==target_reg->query_id) { sum_len += abs(reg->sv_len); i++; }
					else var_vec_contain.erase(var_vec_contain.begin()+i);
				}

				pafaln_cand_vec.clear();
				if(sum_len<abs(var->sv_len)){
					size_ratio = (double)sum_len / abs(var->sv_len);
					dif = abs(var->sv_len) - sum_len;

					//cout << "sum_len=" << sum_len << ", sv_len=" << var->sv_len << ", size_ratio=" << size_ratio << ", dif=" << dif << endl;

					if(size_ratio<CALL_MIN_SEQSIM_MERGE_THRES or dif>MAX_INNER_MISSING_IGNORE_SIZE){ // try recovering
						mm2_aln = minimap2_aln_vec.at(target_reg->minimap2_aln_id);
						start_paf_idx = -1;
						for(i=0; i<mm2_aln->pafalnsegs.size(); i++){
							paf_alnseg = mm2_aln->pafalnsegs.at(i);
							if(paf_alnseg->seglen==abs(target_reg->sv_len) and ((paf_alnseg->opflag==BAM_CINS and target_reg->var_type==VAR_INS) or (paf_alnseg->opflag==BAM_CDEL and target_reg->var_type==VAR_DEL))){
								start_paf_idx = i;
								break;
							}
						}
						if(start_paf_idx!=-1){
							// left part
							for(k=start_paf_idx-1; k>=0; k--){
								paf_alnseg = mm2_aln->pafalnsegs.at(k);
								if(((paf_alnseg->opflag==BAM_CINS and target_reg->var_type==VAR_INS) or (paf_alnseg->opflag==BAM_CDEL and target_reg->var_type==VAR_DEL)) and paf_alnseg->seglen>=min_sv_size){
									if(paf_alnseg->startRpos>=start_ref_pos_tmp){
										contain_flag = false;
										for(auto const& var_tmp : var_vec_contain){
											if(var_tmp->var_type==target_reg->var_type and var_tmp->sv_len==paf_alnseg->seglen){
												contain_flag = true;
												break;
											}
										}
										if(contain_flag==false)
											pafaln_cand_vec.push_back(paf_alnseg);
									}else break;
								}
							}
							// right part
							for(i=start_paf_idx+1; i<mm2_aln->pafalnsegs.size(); i++){
								paf_alnseg = mm2_aln->pafalnsegs.at(i);
								if(((paf_alnseg->opflag==BAM_CINS and target_reg->var_type==VAR_INS) or (paf_alnseg->opflag==BAM_CDEL and target_reg->var_type==VAR_DEL)) and paf_alnseg->seglen>=min_sv_size){
									if(paf_alnseg->startRpos<=start_ref_pos_tmp){
										contain_flag = false;
										for(auto const& var_tmp : var_vec_contain){
											if(var_tmp->var_type==target_reg->var_type and var_tmp->sv_len==paf_alnseg->seglen){
												contain_flag = true;
												break;
											}
										}
										if(contain_flag==false)
											pafaln_cand_vec.push_back(paf_alnseg);
									}else break;
								}
							}

							sum_len2 = 0;
							for(auto const& pafaln_tmp : pafaln_cand_vec) sum_len2 += pafaln_tmp->seglen;

							sum_total = sum_len + sum_len2;
							if(sum_total<abs(var->sv_len)) size_ratio = (double)sum_total / abs(var->sv_len);
							else size_ratio = (double)abs(var->sv_len) / sum_total;

							//cout << "sum_total=" << sum_total << ", sv_len=" << var->sv_len << ", sum_len2=" << sum_len2 << ", size_ratio=" << size_ratio << endl;

							if(size_ratio>=CALL_MIN_SEQSIM_MERGE_THRES){
								for(auto const& pafaln_tmp : pafaln_cand_vec){
									// get refseq and altseq
									start_ref_pos_new = pafaln_tmp->startRpos + 1;
									end_ref_pos_new = getEndRefPosAlnSeg(pafaln_tmp->startRpos, pafaln_tmp->opflag, pafaln_tmp->seglen) + 1;

									reg_str = var->chrname + ":" + to_string(start_ref_pos_new) + "-" + to_string(start_ref_pos_new);
									pthread_mutex_lock(&mutex_fai);
									p_seq = fai_fetch64(fai, reg_str.c_str(), &seq_len);
									pthread_mutex_unlock(&mutex_fai);
									refseq_tmp = p_seq;
									free(p_seq);

									start_localref_pos_new = pafaln_tmp->startSubjectPos + 1;
									end_localref_pos_new = getEndSubjectPosAlnSeg(pafaln_tmp->startSubjectPos, pafaln_tmp->opflag, pafaln_tmp->seglen) + 1;

									start_query_pos_new = pafaln_tmp->startQpos + 1;
									end_query_pos_new = getEndQueryPosAlnSeg(pafaln_tmp->startQpos, pafaln_tmp->opflag, pafaln_tmp->seglen) + 1;
									FastaSeqLoader fa_loader(ctgfilename_para);
									queryseq_tmp = fa_loader.getFastaSeqByPos(mm2_aln->query_id, start_query_pos_new, end_query_pos_new, mm2_aln->relative_strand);

									//write to newVarVec
									reg_new = new reg_t();
									reg_new->var_type = var->var_type;
									reg_new->chrname = var->chrname;
									reg_new->startRefPos = start_ref_pos_new;
									reg_new->endRefPos = end_ref_pos_new;
									reg_new->startLocalRefPos = start_localref_pos_new;
									reg_new->endLocalRefPos = end_localref_pos_new;
									reg_new->startQueryPos = start_query_pos_new;
									reg_new->endQueryPos = end_query_pos_new;
									reg_new->query_id = mm2_aln->query_id;
									if(pafaln_tmp->opflag==BAM_CINS) reg_new->sv_len = pafaln_tmp->seglen;
									else reg_new->sv_len = -pafaln_tmp->seglen;
									reg_new->minimap2_aln_id = mm2_aln->minimap2_aln_id;
									reg_new->call_success_status = true;
									reg_new->short_sv_flag = false;
									reg_new->zero_cov_flag = false;
									reg_new->aln_seg_end_flag = false;
									reg_new->query_pos_invalid_flag = false;
									reg_new->large_indel_flag = false;
									reg_new->merge_flag = false;
									reg_new->rescue_flag = rescue_flag;
									reg_new->aln_orient = mm2_aln->relative_strand;
									reg_new->gt_type = -1;
									reg_new->gt_seq = "";
									reg_new->refseq =  refseq_tmp;
									reg_new->altseq = queryseq_tmp;
									reg_new->AF = 0;
									reg_new->supp_num = target_reg->supp_num;
									reg_new->DP = target_reg->DP;
									reg_new->discover_level = VAR_DISCOV_L_CNS_ALN;

									var_vec.push_back(reg_new);

									change_flag = true;
								}
							}
						}
					}
				}
			}
		}
	}

	if(change_flag) sortRegVec(var_vec); // sort
}

// process the split inter-aligned segments
void varCand::computeIndelFromSplitSegs(vector<reg_t*> &var_vec, vector<minimap2_aln_t*> &minimap2_aln_vec, string &ctgfilename_para, int64_t startRefPos_cns, int64_t endRefPos_cns, double size_ratio_match_thres, double size_ratio_match_thres_merge, bool rescue_flag){
	vector<reg_t*> var_vec2;
	size_t i, j, k;
	minimap2_aln_t *minimap2_aln;
	bool flag;
	vector<int32_t> query_id_vec;
	int32_t query_id;
	reg_t *reg, *reg2;
	double size_ratio;

	// sort according to query_id
	for(i=0; i<minimap2_aln_vec.size(); i++){
		minimap2_aln = minimap2_aln_vec.at(i);
		if(find(query_id_vec.begin(), query_id_vec.end(), minimap2_aln->query_id)==query_id_vec.end()) query_id_vec.push_back(minimap2_aln->query_id);
	}

	for(i=0; i<query_id_vec.size(); i++){
		query_id = query_id_vec.at(i);
		var_vec2 = computeIndelFromSplitSegsSingleQuery(minimap2_aln_vec, ctgfilename_para, query_id, startRefPos_cns, endRefPos_cns, size_ratio_match_thres, size_ratio_match_thres_merge, rescue_flag);
		for(j=0; j<var_vec2.size(); j++) {
			reg = var_vec2.at(j);
			// validate
			flag = false;
			for(k=0; k<var_vec.size(); k++){
				reg2 = var_vec.at(k);
				if(reg->var_type==reg2->var_type){
					if(abs(reg->sv_len)<abs(reg2->sv_len)) size_ratio = (double)abs(reg->sv_len) / abs(reg2->sv_len);
					else size_ratio = (double)abs(reg2->sv_len) / abs(reg->sv_len);
					//cout << "------ size_ratio=" << size_ratio << ", " << reg2->chrname << ":" << reg2->startRefPos << "-" << reg2->endRefPos << ", sv_len=" << reg2->sv_len << ", " << reg->chrname << ":" << reg->startRefPos << "-" << reg->endRefPos << ", sv_len=" << reg->sv_len << endl;
					//if(size_ratio>0.95){ // deleted on 2026-05-19
					if(size_ratio>MIN_ALN_SIZE_RATIO_HIGHSIM){
						flag = true;
						break;
					}
				}
			}
			if(flag==false){
				var_vec.push_back(reg);
			}else delete reg;
		}
	}
}

// process the split inter-aligned segments
vector<reg_t*> varCand::computeIndelFromSplitSegsSingleQuery(vector<minimap2_aln_t*> &minimap2_aln_vec, string &ctgfilename_para, int32_t query_id_para, int64_t startRefPos_cns, int64_t endRefPos_cns, double size_ratio_match_thres, double size_ratio_match_thres_merge, bool rescue_flag){
	vector<reg_t*> var_vec;
	size_t  i, j;
	uint32_t op;
	string chrname_tmp, reg_str, refseq_tmp, queryseq_tmp;
	minimap2_aln_t *minimap2_aln, *minimap2_aln1, *minimap2_aln2;
	struct pafalnSeg* paf_alnseg, *paf_alnseg2;
	int64_t end_ref_pos, max_val, min_val, left_aln_idx, right_aln_idx, left_cg_idx, right_cg_idx;
	int64_t start_ref_pos_tmp, end_ref_pos_tmp, start_query_pos_tmp, end_query_pos_tmp, start_sub_pos_tmp, end_sub_pos_tmp, DP_tmp, suppNum_tmp, var_type_tmp, sv_len_tmp;
	int64_t chrlen_tmp, start_var_pos, end_var_pos;
	char *p_seq;
	bool overlap_flag;
	int32_t seq_len;
	reg_t *reg;
	vector<int32_t> supp_num_vec;

	chrlen_tmp = faidx_seq_len64(fai, chrname.c_str());
	start_var_pos = varVec.at(0)->startRefPos - EXT_SIZE_CHK_VAR_LOC;
	end_var_pos = varVec.at(varVec.size()-1)->endRefPos + EXT_SIZE_CHK_VAR_LOC;
	if(start_var_pos<1) start_var_pos = 1;
	if(end_var_pos>chrlen_tmp) end_var_pos = chrlen_tmp;

	if(minimap2_aln_vec.size()>0){
		// process the split segments
		DP_tmp = 0;
		for(i=0; i<minimap2_aln_vec.size(); i++) DP_tmp += minimap2_aln_vec.at(i)->qname.size();
		// initialize
		max_val = INT_MIN;
		min_val = INT_MAX;
		left_aln_idx = right_aln_idx = left_cg_idx = right_cg_idx = -1;
		for(i=0; i<minimap2_aln_vec.size(); i++){
			minimap2_aln = minimap2_aln_vec.at(i);
			if(minimap2_aln->query_id==query_id_para){ // target query
				if(minimap2_aln->query_end-minimap2_aln->query_start>=VAR_ALN_EXTEND_SIZE or minimap2_aln->subject_end-minimap2_aln->subject_start>=VAR_ALN_EXTEND_SIZE){
					if(minimap2_aln->pafalnsegs.size()>0){
						for(j=0; j<minimap2_aln->pafalnsegs.size(); j++){
							paf_alnseg = minimap2_aln->pafalnsegs.at(j);
							op = paf_alnseg->opflag;

							if((op==BAM_CMATCH or op==BAM_CEQUAL) and paf_alnseg->seglen>MIN_AVER_SIZE_ALN_SEG){
								// get the maximum on the left side, and the minimum on the right side
								end_ref_pos = getEndRefPosAlnSeg(paf_alnseg->startRpos, paf_alnseg->opflag, paf_alnseg->seglen);
								if(paf_alnseg->startRpos<varVec.at(0)->startRefPos and end_ref_pos>=start_var_pos){
									if(max_val<=paf_alnseg->startRpos){
										max_val = paf_alnseg->startRpos;
										left_aln_idx = i;
										left_cg_idx = j;
									}
								}
								if(end_ref_pos>varVec.at(varVec.size()-1)->endRefPos and paf_alnseg->startRpos<=end_var_pos){
									if(min_val>=end_ref_pos and left_aln_idx!=(int64_t)i){
										min_val = end_ref_pos;
										right_aln_idx = i;
										right_cg_idx = j;
									}
								}
							}
						}
					}
				}
			}
		}

		// process the final query
		if(left_aln_idx!=-1 and right_aln_idx!=-1){
			minimap2_aln1 = minimap2_aln_vec.at(left_aln_idx);
			minimap2_aln2 = minimap2_aln_vec.at(right_aln_idx);
			overlap_flag = isOverlappedPos(minimap2_aln1->query_start, minimap2_aln1->query_end, minimap2_aln2->query_start, minimap2_aln2->query_end);
			if(minimap2_aln1->relative_strand==minimap2_aln2->relative_strand and overlap_flag==false){ // same strand and no overlap in query
				paf_alnseg = minimap2_aln1->pafalnsegs.at(left_cg_idx);
				paf_alnseg2 = minimap2_aln2->pafalnsegs.at(right_cg_idx);
				start_ref_pos_tmp = getEndRefPosAlnSeg(paf_alnseg->startRpos, paf_alnseg->opflag, paf_alnseg->seglen) + 1;
				start_sub_pos_tmp = getEndSubjectPosAlnSeg(paf_alnseg->startSubjectPos, paf_alnseg->opflag, paf_alnseg->seglen) + 1;
				start_query_pos_tmp = getEndQueryPosAlnSeg(paf_alnseg->startQpos, paf_alnseg->opflag, paf_alnseg->seglen) + 1;
				end_ref_pos_tmp = paf_alnseg2->startRpos;
				end_query_pos_tmp = paf_alnseg2->startQpos;
				end_sub_pos_tmp = paf_alnseg2->startSubjectPos;
				chrname_tmp = minimap2_aln1->chrname;

				if(start_ref_pos_tmp<=end_ref_pos_tmp and start_query_pos_tmp<=end_query_pos_tmp){

					sv_len_tmp = (end_query_pos_tmp - start_query_pos_tmp) - (end_ref_pos_tmp - start_ref_pos_tmp);
					if(sv_len_tmp>=0) { var_type_tmp = VAR_INS; op = BAM_CINS; }
					else { var_type_tmp = VAR_DEL; op = BAM_CDEL; sv_len_tmp = -sv_len_tmp; }

					supp_num_vec = computeSuppNumFromRegionAlnSegs(minimap2_aln1->qname, op, sv_len_tmp, chrname_tmp, start_ref_pos_tmp, end_ref_pos_tmp, startRefPos_cns, endRefPos_cns, size_ratio_match_thres, size_ratio_match_thres_merge);
					suppNum_tmp = supp_num_vec.at(0);
					DP_tmp = supp_num_vec.at(1);
					if(suppNum_tmp>=minReadsNumSupportSV){
						reg_str = chrname_tmp + ":" + to_string(start_ref_pos_tmp) + "-" + to_string(end_ref_pos_tmp);
						pthread_mutex_lock(&mutex_fai);
						p_seq = fai_fetch(fai, reg_str.c_str(), &seq_len);
						pthread_mutex_unlock(&mutex_fai);
						refseq_tmp = p_seq;
						free(p_seq);

						FastaSeqLoader fa_loader(ctgfilename_para);
						queryseq_tmp = fa_loader.getFastaSeqByPos(minimap2_aln1->query_id, start_query_pos_tmp, end_query_pos_tmp, minimap2_aln1->relative_strand);

						//write to newVarVec
						reg = new reg_t();
						reg->var_type = var_type_tmp;
						reg->chrname = chrname_tmp;
						reg->startRefPos = start_ref_pos_tmp;
						reg->endRefPos = end_ref_pos_tmp;
						reg->startLocalRefPos = start_sub_pos_tmp;
						reg->endLocalRefPos = end_sub_pos_tmp;
						reg->startQueryPos = start_query_pos_tmp;
						reg->endQueryPos = end_query_pos_tmp;
						reg->query_id = minimap2_aln1->query_id;
						reg->sv_len = sv_len_tmp;
						reg->minimap2_aln_id = minimap2_aln1->minimap2_aln_id;
						reg->call_success_status = true;
						reg->short_sv_flag = false;
						reg->zero_cov_flag = false;
						reg->aln_seg_end_flag = false;
						reg->query_pos_invalid_flag = false;
						reg->large_indel_flag = false;
						reg->merge_flag = false;
						reg->rescue_flag = rescue_flag;
						reg->aln_orient = minimap2_aln1->relative_strand;
						reg->gt_type = -1;
						reg->gt_seq = "";
						reg->refseq =  refseq_tmp;
						reg->altseq = queryseq_tmp;
						reg->AF = 0;
						reg->supp_num = minimap2_aln1->qname.size();
						reg->DP = DP_tmp;
						reg->discover_level = VAR_DISCOV_L_CNS_ALN;

						var_vec.push_back(reg);
					}
				}
			}
		}
	}

	return var_vec;
}

// merge neighbouring indels
void varCand::mergeNeighbouringVars(vector<reg_t*> &regVector, int32_t max_ref_dist_arbitary_thres, int32_t max_ref_dist_thres, double min_merge_seqsim_thres, double min_valid_sig_size_ratio_thres, faidx_t *fai, string &contigfilename, string &reffilename){
	size_t i, j;
	int32_t query_id, seq_len, dist;
	int64_t merge_ref_dist_thres, start_querypos_comp, start_refpos_comp, end_refpos_comp, pre_refpos_comp, overlap_len, max_size;
	reg_t *reg_1, *reg_2;
	string comp_queryseq, comp_refseq, reg_str, new_queryseq;
    double sim_val, merge_seqsim_thres, overlap_ratio;
	//int32_t sv_supp_num2 = 0, sv_size = 0, total_sv_num = 0, sv_supp_num1 = 0;
	bool flag;
	char *seq;
	vector<int32_t> query_id_vec;
	vector<reg_t*> reg_vec;
	vector<reg_t*> sorted_reg_vec; 

	dist = -1;
	pre_refpos_comp = -1;

	merge_seqsim_thres = min_merge_seqsim_thres;
	merge_ref_dist_thres = max_ref_dist_thres;

//	complex_reg_flag = false;
//	if(regVector.size()>=2) complex_reg_flag = ComplexRegionFlag(regVector, 10, 5, 800);

	// sort according to query_id
	for(i=0; i<regVector.size(); i++){
		reg_1 = regVector.at(i);
		if(find(query_id_vec.begin(), query_id_vec.end(), reg_1->query_id)==query_id_vec.end()) query_id_vec.push_back(reg_1->query_id);
	}
	for(i=0; i<query_id_vec.size(); i++){
		query_id = query_id_vec.at(i);
		vector<reg_t*> group;
		for(j = 0; j < regVector.size(); j++){
			reg_1 = regVector.at(j);
			if(reg_1->query_id == query_id) {
				group.push_back(reg_1);
			}
		}
		sort(group.begin(), group.end(), [](reg_t* a, reg_t* b){
			return a->startRefPos < b->startRefPos;
		});

		sorted_reg_vec.insert(sorted_reg_vec.end(), group.begin(), group.end());
	}
	regVector = sorted_reg_vec;
	
	// merge
	for(i=0; i<regVector.size(); i++){
		reg_1 = regVector.at(i);
		query_id = reg_1->query_id;

		if(i+1<regVector.size() and reg_1->startQueryPos!=-1 and reg_1->endQueryPos!=-1){
			reg_2 = regVector.at(i+1);

			if(query_id == reg_2->query_id and reg_1->var_type == reg_2->var_type and reg_2->startQueryPos!=-1 and reg_2->endQueryPos!=-1){
				if(pre_refpos_comp<0) dist = reg_2->startRefPos - reg_1->endRefPos - 1;
				else dist = reg_2->startRefPos - pre_refpos_comp - 1;

				overlap_len = 0;
				overlap_ratio = 0 ;
				max_size = max(abs(reg_1->sv_len), abs(reg_2->sv_len));
				
				overlap_len = reg_1->endRefPos - reg_2->startRefPos + 1;
				if(overlap_len > 0){
					dist = - overlap_len;
				}else if(pre_refpos_comp<0) dist = reg_2->startRefPos - reg_1->endRefPos - 1;
				else dist = reg_2->startRefPos - pre_refpos_comp - 1;
				

				//cout << "svtype=" << to_string(reg_1->var_type) << ", reg_1: sv_len=" << reg_1->sv_len << ", reg_2: sv_len=" << reg_2->sv_len << ", dist=" << dist << endl;

				//update merge_distance_thres, deleted on 2026-05-12
//				if(dist > merge_ref_dist_thres){
//					sv_size = abs(reg_2->sv_len);
//					total_sv_num = regVector.size();
//
//					sv_supp_num1 = reg_1->supp_num;
//					sv_supp_num2 = reg_2->supp_num;
//					sv_supp_diff = 1 - (double)min(sv_supp_num1, sv_supp_num2)/max(sv_supp_num1, sv_supp_num2);
//					merge_ref_dist_thres = calculate_merge_distance_threshold(complex_reg_flag, sv_supp_num2, sv_size, total_sv_num, sv_supp_diff);
//					//cout << "sv_supp_diff=" << sv_supp_diff << ", dist_thres=" << merge_ref_dist_thres << endl;
//				}

				// overlaped deletion and insertion
				if(dist<0 and reg_1->var_type == BAM_CDEL){ // overlaped deletion
					// cout << "svtype:" << to_string(reg_1->var_type) << endl;
					// cout << "[" << complex_reg_flag << "] (" << reg_1->startRefPos << "," << reg_1->sv_len << "),(" << reg_2->startRefPos << "," << reg_2->sv_len << ") " << dist << "," << pre_refpos_comp << endl;
					flag = false;
					overlap_ratio = (double)overlap_len/abs(reg_1->sv_len);
					if((double)overlap_len/abs(reg_2->sv_len) < overlap_ratio){
						flag = true;
						overlap_ratio = (double)overlap_len/abs(reg_2->sv_len);
					}
					if(flag){
						delete reg_1;
						regVector.erase(regVector.begin()+i);
						i--;
					}else{
						delete reg_2;
						regVector.erase(regVector.begin()+i+1);
						i--;
					}
					//pre_refpos_comp = reg_2->endRefPos;
					continue;
				}else if(reg_1->endQueryPos > reg_2->startQueryPos and reg_1->var_type == BAM_CINS){ // overlaped insertion
					// cout << "svtype:" << reg_1->var_type << endl;
					// cout << "[" << complex_reg_flag << "] (" << reg_1->startRefPos << "," << reg_1->sv_len << "),(" << reg_2->startRefPos << "," << reg_2->sv_len << ") " << dist << "," << pre_refpos_comp << endl;			
					if(i+2<regVector.size()){
						reg_t *reg_3 = regVector.at(i+2);
						if(reg_3->startQueryPos!=-1 and reg_3->endQueryPos!=-1 and reg_1->endQueryPos > reg_3->startQueryPos){
							if(reg_2->endQueryPos > reg_3->startQueryPos){
								delete reg_2;
								regVector.erase(regVector.begin()+i+1);
								delete reg_1;
								regVector.erase(regVector.begin()+i);
								i++;
								//pre_refpos_comp = reg_3->endRefPos;
								continue;
							}else{
								delete reg_1;
								regVector.erase(regVector.begin()+i);
								// reg_1 = reg_2;
								// reg_2 = reg_3;
								// dist = reg_2->startRefPos - reg_1->startRefPos;
								i--;
								pre_refpos_comp = reg_2->endRefPos;
								continue;
							}
						}else{
							//pre_refpos_comp = reg_2->endRefPos;
							delete reg_2;
							regVector.erase(regVector.begin()+i+1);
							// reg_2 = reg_3;
							// dist = reg_2->startRefPos - reg_1->startRefPos;
							i--;
							continue;
						}
					}else{
						if(abs(reg_1->sv_len)>=abs(reg_2->sv_len)){
							//pre_refpos_comp = reg_2->endRefPos;
							delete reg_2;
							regVector.erase(regVector.begin()+i+1);
							break;
						}
					}
				}
				// cout << dist << endl;
				if(dist <= merge_ref_dist_thres and dist>=0){
					// min_seqsim_merge = 0.85? 0.9?
					// update merge_seqsim_thres
//					if(dist >= reg_2->sv_len){ // deleted on 2025-02-11
//						if(pre_refpos_comp<0){
//							common_len = dist - reg_2->sv_len;
//							merge_seqsim_thres = calculate_merge_seqsim_threshold(common_len, dist, min_merge_seqsim_thres);
//						}else{
//							pos1_tmp = reg_1->endQueryPos + dist;
//							pos2_tmp = reg_2->startQueryPos - (dist - reg_2->sv_len);
//							if(pos1_tmp > pos2_tmp){
//								common_len = pos1_tmp - pos2_tmp;
//								merge_seqsim_thres = calculate_merge_seqsim_threshold(common_len, dist, min_merge_seqsim_thres);
//								// merge_seqsim_thres = calculate_merge_seqsim_threshold2(common_len, distance, query_sv_type);
//							}
//						}
//					}
					if(dist<max_ref_dist_arbitary_thres) merge_seqsim_thres = QC_SEQSIM_RATIO_THRES_FACTOR * min_merge_seqsim_thres;

					// get compared queryseq
					start_querypos_comp = reg_2->endQueryPos - dist;
					FastaSeqLoader cns_fa_loader(contigfilename);
					comp_queryseq = cns_fa_loader.getFastaSeq(query_id).substr(start_querypos_comp, dist + 1);

					// get ref_comp seq
					start_refpos_comp = reg_2->endRefPos - dist;
					end_refpos_comp = reg_2->endRefPos;
					reg_str = reg_2->chrname + ":" + to_string(start_refpos_comp) + "-" + to_string(end_refpos_comp);
					pthread_mutex_lock(&mutex_fai);
					seq = fai_fetch(fai, reg_str.c_str(), &seq_len);
					pthread_mutex_unlock(&mutex_fai);
					comp_refseq = seq;
					free(seq);

					sim_val = computeVarseqSim(comp_queryseq, comp_refseq);
					//cout << "sim_val=" << sim_val << ", comp_queryseq.size=" << comp_queryseq.size() << ", comp_refseq.size=" << comp_refseq.size() << ", comp_queryseq=" << comp_queryseq << ", comp_refseq=" << comp_refseq << endl;

					if(reg_1->var_type==BAM_CDEL){ // DEL
						//cout << "DEL:sim_val=" << sim_val << ", merge_seqsim_thres=" << merge_seqsim_thres << endl;

						//if(((double)dist/(abs(reg_1->sv_len)+abs(reg_2->sv_len)+dist)<3*min_valid_sig_size_ratio_thres and dist<max_ref_dist_arbitary_thres) or sim_val >= merge_seqsim_thres) {
						if(((double)dist/(abs(reg_1->sv_len)+abs(reg_2->sv_len)+dist)<3*min_valid_sig_size_ratio_thres and dist<max_ref_dist_arbitary_thres) or (sim_val >= merge_seqsim_thres and (1-sim_val)*seq_len<max_size)) {
							reg_1->endRefPos += abs(reg_2->sv_len) - 1;
							reg_1->endLocalRefPos += abs(reg_2->sv_len) - 1;
							reg_1->sv_len += reg_2->sv_len;
							reg_1->refseq += reg_2->refseq;
							reg_1->merge_flag = true;
							pre_refpos_comp = reg_2->endRefPos;
							delete reg_2;
							regVector.erase(regVector.begin()+i+1);
							i--;
						}else pre_refpos_comp = -1;
					}else if(reg_1->var_type==BAM_CINS){ // INS
						//cout << "INS:sim_val=" << sim_val << ", merge_seqsim_thres=" << merge_seqsim_thres << endl;

						//if(((double)dist/(abs(reg_1->sv_len)+abs(reg_2->sv_len)+dist)<3*min_valid_sig_size_ratio_thres and dist<max_ref_dist_arbitary_thres) or sim_val >= merge_seqsim_thres){
						if(((double)dist/(abs(reg_1->sv_len)+abs(reg_2->sv_len)+dist)<3*min_valid_sig_size_ratio_thres and dist<max_ref_dist_arbitary_thres) or (sim_val >= merge_seqsim_thres and (1-sim_val)*seq_len<max_size)){
							// if(reg_1->endQueryPos > reg_2->startQueryPos){
							// 	cout << "1:" << reg_1->startRefPos << "," << reg_1->endRefPos << " 2:" << reg_2->startRefPos << "," << reg_2->endRefPos << endl;
							// 	cout << "1:" << reg_1->startQueryPos << "," << reg_1->endQueryPos << " 2:" << reg_2->startQueryPos << "," << reg_2->endQueryPos << endl;
							// 	cout << "1:" << reg_1->sv_len << " 2:" << reg_2->sv_len << endl;
							// 	cout << "-------------------" << endl;
							// }
							new_queryseq = "";
							reg_1->sv_len += reg_2->sv_len;
							pre_refpos_comp = reg_2->endRefPos;

							start_querypos_comp = reg_1->endQueryPos;
							FastaSeqLoader cns_fa_loader(contigfilename);
							new_queryseq = cns_fa_loader.getFastaSeq(query_id).substr(start_querypos_comp, abs(reg_2->sv_len) - 1);

							reg_1->altseq += new_queryseq;
							reg_1->endQueryPos = reg_1->endQueryPos + abs(reg_2->sv_len);
							reg_1->merge_flag = true;
							delete reg_2;
							regVector.erase(regVector.begin()+i+1);
							i--;
						}else pre_refpos_comp = -1;
					}
				}else pre_refpos_comp = -1;
			}else pre_refpos_comp = -1;
		}
	}
	regVector.shrink_to_fit();
}

// rescue indel variants
vector<reg_t*> varCand::rescueIndelVarLoc(vector<int32_t> &clusterId_incomplete){
	vector<reg_t*> var_vec;
	vector<vector<string>> qnames_vec;
	size_t i, j, cluster_id, id, id2, n_seqs;
	vector<string> qname_vec;
	vector<clipAlnData_t*> clipAlnDataVector;
	bool flag;
	int64_t start_var_pos, end_var_pos, startRefPos_cns, endRefPos_cns, chrlen_tmp, mem_avail;
	ofstream outfile_rescue_reads, outfile_rescue_refseq, outfile_rescue_cns;
	string tmp_cns_filename, refseq, reg_str, cmd, cmd2, cons_header, minimap2_cmd, output_prefix, tmp_reg_str;
	int32_t seq_len, ret_status, status, serial_number, left_shift_size, right_shift_size, len;
	char *p_seq;
	vector<struct querySeqInfoNode*> query_seq_info_all;
	struct seqsVec *smoothed_seqs;

	// load cluster reads
	qnames_vec = getClusterInfo(clusterfilename);
	if(qnames_vec.size()==0) return var_vec;

	start_var_pos = varVec.at(0)->startRefPos;
	end_var_pos = varVec.at(varVec.size()-1)->endRefPos;

	chrlen_tmp = faidx_seq_len64(fai, chrname.c_str());  // get reference size
	startRefPos_cns = start_var_pos - ref_left_shift_size;
	endRefPos_cns = end_var_pos + ref_right_shift_size;
	if(startRefPos_cns<1) startRefPos_cns = 1;
	if(endRefPos_cns>chrlen_tmp) endRefPos_cns = chrlen_tmp;
	left_shift_size = start_var_pos - startRefPos_cns;  // left shift size
	right_shift_size = endRefPos_cns - end_var_pos;  // right shift size

	rescue_refseqfilename = out_dir_call + "/rescue_refseq_" + chrname + "_" + to_string(start_var_pos) + "-" + to_string(end_var_pos) + ".fa";
	reg_str = chrname + ":" + to_string(startRefPos_cns) + "-" + to_string(endRefPos_cns);
	pthread_mutex_lock(&mutex_fai);
	p_seq = fai_fetch(fai, reg_str.c_str(), &seq_len);
	pthread_mutex_unlock(&mutex_fai);
	refseq = p_seq;
	free(p_seq);

	outfile_rescue_refseq.open(rescue_refseqfilename);
	if(outfile_rescue_refseq.is_open()==false){
		cerr << "In " << __func__ << "(), cannot open file '" << rescue_refseqfilename << "', error!" << endl;
		exit(1);
	}

	// refseq file
	cons_header = ">" + reg_str + "___" + to_string(left_shift_size) + "___" + to_string(right_shift_size) + "___" + rescue_refseqfilename;
	outfile_rescue_refseq << cons_header << endl;
	outfile_rescue_refseq << refseq << endl;
	outfile_rescue_refseq.close();

	// cnsfile
	rescue_cnsfilename = out_dir_call + "/rescue_cns_" + chrname + "_" + to_string(start_var_pos) + "-" + to_string(end_var_pos) + ".fa";
	outfile_rescue_cns.open(rescue_cnsfilename);
	if(outfile_rescue_cns.is_open()==false){
		cerr << "In " << __func__ << "(), cannot open file '" << rescue_cnsfilename << "', error!" << endl;
		exit(1);
	}

	rescue_alnfilename = out_dir_call + "/rescue_minimap2_" + chrname + "_" + to_string(start_var_pos) + "-" + to_string(end_var_pos) + ".paf";

	// load data according to clusters
	serial_number = 1;
	for(cluster_id=0; cluster_id<qnames_vec.size(); cluster_id++){
		if(find(clusterId_incomplete.begin(), clusterId_incomplete.end(), cluster_id) == clusterId_incomplete.end()){  // call complete, then skip
			continue;
		}

		qname_vec = qnames_vec.at(cluster_id);
		clipAlnDataLoader data_loader(chrname, startRefPos_cns, endRefPos_cns, inBamFile, minClipEndSize, maxVarRegSize, minMapQ, minHighMapQ);
		data_loader.loadClipAlnDataWithSATag(clipAlnDataVector, qname_vec);

		// extract queries from clip align data vector
		query_seq_info_all = extractQueriesFromClipAlnDataVec(clipAlnDataVector, refseq, chrname, startRefPos_cns, endRefPos_cns, start_var_pos, end_var_pos, fai, minConReadLen, clip_reg_flag, &mutex_fai);

		if(query_seq_info_all.size()>0){
			// construct the reads file
			rescue_readsfilename = out_dir_call + "/rescue_reads_" + chrname + "_" + to_string(start_var_pos) + "-" + to_string(end_var_pos) + "_" + to_string(cluster_id) + ".fa";
//			cout << rescue_readsfilename << endl;

			outfile_rescue_reads.open(rescue_readsfilename);
			if(outfile_rescue_reads.is_open()==false){
				cerr << "In " << __func__ << "(), cannot open file '" << rescue_readsfilename << "', error!" << endl;
				exit(1);
			}

			//smooth read segments
			smoothed_seqs = smoothQuerySeqData(refseq, query_seq_info_all, startRefPos_cns, min_sv_size);

			n_seqs = smoothed_seqs->seqs.size();
			if(n_seqs>=(size_t)minReadsNumSupportSV*READS_NUM_SUPPORT_FACTOR){ // only use sufficient reads data to generate consensus sequence
				// reads file
				for(i=0; i<n_seqs; i++) outfile_rescue_reads << ">" << smoothed_seqs->qname.at(i) << endl << smoothed_seqs->seqs.at(i) << endl;
				outfile_rescue_reads.close();

				flag = false;
				while(1){
					pthread_mutex_lock(&mutex_mem);
					mem_avail = getMemInfo("MemAvailable", 2);
					// prepare for alignment computation
					if(mem_total-mem_avail<=mem_total*mem_use_block_factor+extend_total*extend_use_block_factor){
						work_num ++;
						flag = true;
					}
					pthread_mutex_unlock(&mutex_mem);

					if(flag==false){ // block the alignment computation
						//cout << "\t" << __func__ << ": wait " << mem_wait_seconds << " seconds, mem_avail=" << (mem_avail >> 10) << " MB, work_num=" << work_num << ", " << contigfilename << endl;
						sleep(mem_wait_seconds);
					}else{
						id = rescue_readsfilename.find_last_of("/");
						id2 = rescue_readsfilename.find_last_of(".");
						len = id2 - id;
						if(id!=string::npos and id2!=string::npos) tmp_reg_str = rescue_readsfilename.substr(id+1, len-1);
						else{
							cerr << "line=" << __LINE__ << ", cannot find the / in " << rescue_readsfilename << ", error." << endl;
							exit(1);
						}

						//output_prefix = "tmp_rescue_wtdbg2_" + tmp_reg_str; // deleted on 2026-01-22
						output_prefix = pg_runid_str + "_tmp_rescue_wtdbg2_" + tmp_reg_str;
						cmd = "wtdbg2.pl -t 8 -x " + technology + " -o " + output_prefix + " -a -q " + rescue_readsfilename + " > /dev/null 2>&1";
						//cout << cmd << endl;
						system(cmd.c_str());

						pthread_mutex_lock(&mutex_mem);
						work_num --;
						if(work_num<0){
							cerr << "line=" << __LINE__ << ", work_num=" << work_num << ", error." << endl;
							exit(1);
						}
						pthread_mutex_unlock(&mutex_mem);

						break;
					}
				}

				tmp_cns_filename = output_prefix + ".cns.fa";
				flag = isFileExist(tmp_cns_filename);
				if(flag){
					FastaSeqLoader fa_loader(tmp_cns_filename);
					for(i=0; i<fa_loader.getFastaSeqCount(); i++){
						if(fa_loader.getFastaSeqLen(i)>0){
							cons_header = ">wtdbg2_rescueCns_";
							cons_header += to_string(serial_number) + "_"; //serial number
							cons_header += to_string(fa_loader.getFastaSeqLen(i)) + "_"; //consensus sequence length
							cons_header += to_string(smoothed_seqs->seqs.size()) + "___";  //read count

							for(j=0; j<smoothed_seqs->qname.size()-1; j++) cons_header += smoothed_seqs->qname.at(j) + ";";
							cons_header += smoothed_seqs->qname.at(smoothed_seqs->qname.size()-1);

							pthread_mutex_lock(&mutex_write);
							outfile_rescue_cns << cons_header << endl; // header
							outfile_rescue_cns << fa_loader.getFastaSeq(i) << endl;  // seq

//							cout << cons_header << endl; // header
//							cout << fa_loader.getFastaSeq(i) << endl;  // seq
							//cons_file.flush();
							pthread_mutex_unlock(&mutex_write);

							serial_number ++;
						}
					}
				}

				cmd2 = "rm -rf " + output_prefix + "*";
				system(cmd2.c_str());  // remove temporary files
				cmd2 = "rm -rf " + rescue_readsfilename;
				system(cmd2.c_str());  // remove reads file
			}else{
				outfile_rescue_reads.close();
			}
			delete smoothed_seqs;
		}

		if(!clipAlnDataVector.empty()) destoryClipAlnData(clipAlnDataVector);
		if(!query_seq_info_all.empty()) destoryQuerySeqInfoAll(query_seq_info_all);
	}
	outfile_rescue_cns.close();

	flag = isFileExist(rescue_cnsfilename);
	if(flag){
		minimap2_cmd = "minimap2 -c -x asm20 -o " + rescue_alnfilename + " " + rescue_refseqfilename + " " + rescue_cnsfilename + " > /dev/null 2>&1";
		status = system(minimap2_cmd.c_str());
		ret_status = getSuccessStatusSystemCmd(status);
		if(ret_status==0){ // command executed successfully
			// parse alignment information
			minimap2_aln_vec = minimap2Parse(rescue_alnfilename, rescue_cnsfilename, rescue_refseqfilename);

			// get reference sequence and query sequence
			var_vec = computeIndelVarLoc(minimap2_aln_vec, rescue_cnsfilename, rescue_refseqfilename, true);
		}
	}

	return var_vec;
}

vector<int32_t> varCand::computeSuppNumFromRegionAlnSegs(vector<string> &clu_qname_vec, struct pafalnSeg* paf_alnseg, vector<clipAlnData_t*> &clipAlnDataVector, string &chrname, int64_t startRefPos_cns, int64_t endRefPos_cns, double size_ratio_match_thres, double size_ratio_match_thres_merge){
	vector<int32_t> supp_num_vec, left_check_idx_vec, right_check_idx_vec;
	int32_t var_num, query_num, num, idx;
	size_t i, j;
	string qname, reg_str, refseq, clu_qname;
	int32_t seq_len, bam_type;
	char *p_seq;
	vector<clipAlnData_t*> query_aln_segs;
	bool no_otherchrname_flag, flag, overlap_flag, left_flag, right_flag;
	double size_ratio;
	int64_t noHardClipIdx, ref_dist, min_ref_dist, end_ref_pos, end_paf_pos, end_seg_pos;
	int64_t start_ref_pos_check, end_ref_pos_check, start_idx_check, sum_len;
	vector<struct alnSeg*> query_alnSegs;
	struct alnSeg *seg;
	vector<string>::iterator q;

	var_num = query_num = 0;
	if (clipAlnDataVector.size() > 0){
		for (i = 0; i < clipAlnDataVector.size(); i++) clipAlnDataVector.at(i)->query_checked_flag = false;

		reg_str = chrname + ":" + to_string(startRefPos_cns) + "-" + to_string(endRefPos_cns);
		pthread_mutex_lock(&mutex_fai);
		p_seq = fai_fetch(fai, reg_str.c_str(), &seq_len);
		pthread_mutex_unlock(&mutex_fai);
		refseq = p_seq;
		free(p_seq);

		if(seq_len>0){
			bam_type = getBamType(clipAlnDataVector);
			if(bam_type!=BAM_INVALID){
				//average_varlen = 0;
				for (i = 0; i < clipAlnDataVector.size(); i++) {//for each read query
					qname = clipAlnDataVector.at(i)->queryname;
					if (clipAlnDataVector.at(i)->query_checked_flag == false) {
						query_aln_segs = getQueryClipAlnSegs(qname, clipAlnDataVector); // get query clip align segments
						if(query_aln_segs.size()==1 and query_aln_segs.at(0)->query_dist<MIN_SPAN_SINGLE_QUERY and query_aln_segs.at(0)->query_dist<0.5*query_aln_segs.at(0)->querylen) continue; // skip low quality segments

						no_otherchrname_flag = true;
						end_paf_pos = getEndRefPosAlnSeg(paf_alnseg->startRpos, paf_alnseg->opflag, paf_alnseg->seglen);
						noHardClipIdx = getNoHardClipAlnItem(query_aln_segs, chrname, paf_alnseg->startRpos, end_paf_pos);
						if (noHardClipIdx != -1) {
							for (int32_t k = 0; k < (int32_t)query_aln_segs.size(); k++) {
								if (k != noHardClipIdx) {
									if (query_aln_segs.at(k)->chrname.compare(varVec[0]->chrname) != 0) {
										no_otherchrname_flag = false;
									}
								}
							}
							if (no_otherchrname_flag) {
								query_aln_segs.at(noHardClipIdx)->query_checked_flag = true;
								if (!(query_aln_segs.at(noHardClipIdx)->bam->core.flag & BAM_FUNMAP)) { // aligned
									switch (bam_type) {
										case BAM_CIGAR_NO_DIFF_MD:
											//alnSegs = generateAlnSegs(b);
											query_alnSegs = generateAlnSegs2(query_aln_segs.at(noHardClipIdx)->bam, startRefPos_cns, endRefPos_cns);
											//queryseq = "=ACMGRSVTWYHKDBN"[bam_get_seq(query_aln_segs.at(noHardClipIdx)->bam)];
											break;
										case BAM_CIGAR_NO_DIFF_NO_MD:
										case BAM_CIGAR_DIFF_MD:
										case BAM_CIGAR_DIFF_NO_MD:
											//alnSegs = generateAlnSegs_no_MD(b, refseq, startPos, endPos);
											query_alnSegs = generateAlnSegs_no_MD2(query_aln_segs.at(noHardClipIdx)->bam, refseq, startRefPos_cns, endRefPos_cns);
											break;
										default:
											cerr << __func__ << ": unknown bam type, error!" << endl;
											exit(1);
									}
								}

								//compute query_num
								//end_paf_pos = getEndRefPosAlnSeg(paf_alnseg->startRpos, paf_alnseg->opflag, paf_alnseg->seglen);
								end_seg_pos = getEndRefPosAlnSeg(query_alnSegs.at(query_alnSegs.size()-1)->startRpos, query_alnSegs.at(query_alnSegs.size()-1)->opflag, query_alnSegs.at(query_alnSegs.size()-1)->seglen);
								if(paf_alnseg->startRpos>query_alnSegs.at(0)->startRpos and end_paf_pos<end_seg_pos){
								//if(paf_alnseg->startRpos>query_alnSegs.at(0)->startRpos and (paf_alnseg->startRpos+paf_alnseg->seglen)<(query_alnSegs.at(query_alnSegs.size()-1)->startRpos + query_alnSegs.at(query_alnSegs.size()-1)->seglen)){
									//cout << "query_num=" << query_num << ", qname=" << qname << ", i=" << i << endl;
									query_num ++;
								}

								//find qname
								q = find(clu_qname_vec.begin(), clu_qname_vec.end(), qname);
								if(q==clu_qname_vec.end()){ // not found
									clipAlnDataVector.at(i)->query_checked_flag = true;
								}else{ // found
									//search
									min_ref_dist = INT_MAX;
									for(j=0; j<query_alnSegs.size(); j++){
										seg = query_alnSegs.at(j);
										// recompute the overlap
										end_ref_pos = getEndRefPosAlnSeg(seg->startRpos, seg->opflag, seg->seglen);
										flag = false;
										if(seg->startRpos-1>=startRefPos_cns and end_ref_pos-1<=endRefPos_cns) flag = true;

										if(flag and paf_alnseg->opflag==seg->opflag){
											if(paf_alnseg->seglen>=min_sv_size and seg->seglen>=min_sv_size){
												if(paf_alnseg->seglen <= seg->seglen) size_ratio = (double)paf_alnseg->seglen / seg->seglen;
												else size_ratio = (double)seg->seglen / paf_alnseg->seglen;
												if(size_ratio>=size_ratio_match_thres){
													if(paf_alnseg->opflag==BAM_CINS)  // INS
														ref_dist = abs(paf_alnseg->startRpos - seg->startRpos);
													else{ // DEL
														overlap_flag = isOverlappedPos(paf_alnseg->startRpos, end_paf_pos, seg->startRpos, end_ref_pos);
														if(overlap_flag) ref_dist = 0;
														else{
															if(paf_alnseg->startRpos<seg->startRpos) ref_dist = seg->startRpos - end_paf_pos;
															else ref_dist = paf_alnseg->startRpos - end_ref_pos;
														}
													}
													//if(ref_distance<200){
													//if(ref_dist<=MAX_REF_DIST_SAME_CHR){
													if(ref_dist<=EXT_SIZE_CHK_VAR_LOC_SMALL){
														//cout << "\txxxxxxx var_num=" << var_num << ", ref_dist=" << ref_dist << ", size_ratio=" << size_ratio << ", seg.pos=" << seg->startRpos << ", seg.len=" << seg->seglen << ", qname=" << qname << ", startRpos=" << paf_alnseg->startRpos << ", seglen=" << paf_alnseg->seglen << ", i=" << i << ", j=" << j << endl;
														if(ref_dist<min_ref_dist) min_ref_dist = ref_dist;
													}
												}
											}
										}
									}
									if(min_ref_dist!=INT_MAX){
										var_num ++;
										//cout << "\tyyyyyyy var_num=" << var_num << ", min_ref_dist=" << min_ref_dist << ", qname=" << qname << ", i=" << i << ", j=" << j << endl;
									}else{ // added on 2024-12-31, try to deal with the partial signatures
										//cout << qname << ": query_alnSegs.size=" << query_alnSegs.size() << ", paf_alnseg.start=" << paf_alnseg->startRpos << ", paf.seglen=" << paf_alnseg->seglen << endl;

										// compute the start_idx, end_idx, and then search partial variants around 500 bp regions
										start_ref_pos_check = paf_alnseg->startRpos - EXT_SIZE_CHK_VAR_LOC;
										end_ref_pos_check = end_paf_pos + EXT_SIZE_CHK_VAR_LOC;
										if(start_ref_pos_check<1) start_ref_pos_check = 1;

										start_idx_check = 0;
										for(j=0; j<query_alnSegs.size(); j++){
											seg = query_alnSegs.at(j);
											if(seg->startRpos>=paf_alnseg->startRpos) {
												start_idx_check = j;
												break;
											}
										}

										left_check_idx_vec.clear();
										right_check_idx_vec.clear();
										left_flag = right_flag = true;
										num = sum_len = 0;
										while(start_idx_check-num>=0 or start_idx_check+num<(int64_t)query_alnSegs.size()){
											flag = false;
											// check left side
											idx = start_idx_check - num;
											if(idx>=0){
												seg = query_alnSegs.at(idx);
												if(left_flag and seg->opflag==paf_alnseg->opflag and seg->seglen>=min_sv_size){
													if(seg->startRpos>start_ref_pos_check){
														left_check_idx_vec.insert(left_check_idx_vec.begin(), idx);
														sum_len += seg->seglen;
														//cout << "^^^^^^^ sum_len=" << sum_len << ", seg.startRpos=" << seg->startRpos << ", seg.seglen=" << seg->seglen << ", seg.opflag=" << seg->opflag << endl;
														flag = true;
													}else left_flag = false;
												}
											}

											// check right side
											if(num>0){
												idx = start_idx_check + num;
												if(idx<(int64_t)query_alnSegs.size()){
													seg = query_alnSegs.at(idx);
													if(right_flag and seg->opflag==paf_alnseg->opflag and seg->seglen>=min_sv_size){
														if(seg->startRpos<end_ref_pos_check){
															right_check_idx_vec.insert(right_check_idx_vec.begin(), idx);
															sum_len += seg->seglen;
															//cout << "&&&&&&&&&&&& sum_len=" << sum_len << ", seg.startRpos=" << seg->startRpos << ", seg.seglen=" << seg->seglen << ", seg.opflag=" << seg->opflag << endl;
															flag = true;
														}else right_flag = false;
													}
												}
											}

											// validate the sum_len
											if(flag){
												if(paf_alnseg->seglen <= sum_len) size_ratio = (double)paf_alnseg->seglen / sum_len;
												else size_ratio = (double)sum_len / paf_alnseg->seglen;
												if(size_ratio>=size_ratio_match_thres_merge){
													var_num ++;
													//cout << "\t-------- var_num=" << var_num << ", sum_len=" << sum_len << ", size_ratio=" << size_ratio << ", i=" << i << ", j=" << j << endl;
													break;
												}
											}

											num ++;
											if(left_flag==false and right_flag==false) break;
										}
									}
								}
							}
						}else{//all hard clip
							for(j=0; j<query_aln_segs.size(); j++) query_aln_segs.at(j)->query_checked_flag = true;
						}
						destroyAlnSegs(query_alnSegs);
					}
				}
			}
		}
	}

	supp_num_vec.push_back(var_num);
	supp_num_vec.push_back(query_num);

	return supp_num_vec;
}

vector<int32_t> varCand::computeSuppNumFromRegionAlnSegs(vector<string> &clu_qname_vec, int32_t opflag_para, int32_t oplen_para, string &chrname_para, int64_t start_var_pos_para, int64_t end_var_pos_para, int64_t startRefPos_cns, int64_t endRefPos_cns, double size_ratio_match_thres, double size_ratio_match_thres_merge){
	vector<int32_t> supp_num_vec, left_check_idx_vec, right_check_idx_vec;
	int32_t var_num, query_num, num, idx;
	size_t i, j;
	string qname, reg_str, refseq, clu_qname;
	int32_t seq_len, bam_type;
	char *p_seq;
	vector<clipAlnData_t*> query_aln_segs, clipAlnDataVector;
	bool no_otherchrname_flag, flag, left_flag, right_flag;
	double size_ratio;
	int64_t noHardClipIdx, ref_dist, min_ref_dist, end_ref_pos, end_seg_pos;
	int64_t start_ref_pos_check, end_ref_pos_check, start_idx_check, sum_len;
	vector<struct alnSeg*> query_alnSegs;
	struct alnSeg *seg;
	vector<string>::iterator q;

	// load the clipping data
	clipAlnDataLoader data_loader(chrname, startRefPos_cns, endRefPos_cns, inBamFile, minClipEndSize, maxVarRegSize, minMapQ, minHighMapQ);
	data_loader.loadClipAlnDataWithSATag(clipAlnDataVector, max_ultra_high_cov, clu_qname_vec);

	var_num = query_num = 0;
	if (clipAlnDataVector.size() > 0){
		for (i = 0; i < clipAlnDataVector.size(); i++) clipAlnDataVector.at(i)->query_checked_flag = false;

		reg_str = chrname_para + ":" + to_string(startRefPos_cns) + "-" + to_string(endRefPos_cns);
		pthread_mutex_lock(&mutex_fai);
		p_seq = fai_fetch(fai, reg_str.c_str(), &seq_len);
		pthread_mutex_unlock(&mutex_fai);
		refseq = p_seq;
		free(p_seq);

		if(seq_len>0){
			bam_type = getBamType(clipAlnDataVector);
			if(bam_type!= BAM_INVALID){
				for (i = 0; i < clipAlnDataVector.size(); i++) {//for each read query
					qname = clipAlnDataVector.at(i)->queryname;
					//find qname
		//			q=find(clu_qname_vec.begin(), clu_qname_vec.end(), qname);
		//			if(q==clu_qname_vec.end()){
		//				clipAlnDataVector.at(i)->query_checked_flag = true;
		//			}

					if (clipAlnDataVector.at(i)->query_checked_flag == false) {
						query_aln_segs = getQueryClipAlnSegs(qname, clipAlnDataVector); // get query clip align segments
						no_otherchrname_flag = true;
						noHardClipIdx = getNoHardClipAlnItem(query_aln_segs, chrname_para, start_var_pos_para, end_var_pos_para);
						if (noHardClipIdx != -1) {
							for (int32_t k = 0; k < (int32_t)query_aln_segs.size(); k++) {
								if (k != noHardClipIdx) {
									if (query_aln_segs.at(k)->chrname.compare(varVec[0]->chrname) != 0) {
										no_otherchrname_flag = false;
									}
								}
							}
							if (no_otherchrname_flag) {
								query_aln_segs.at(noHardClipIdx)->query_checked_flag = true;
								if (!(query_aln_segs.at(noHardClipIdx)->bam->core.flag & BAM_FUNMAP)) { // aligned
									switch (bam_type) {
										case BAM_CIGAR_NO_DIFF_MD:
											//alnSegs = generateAlnSegs(b);
											query_alnSegs = generateAlnSegs2(query_aln_segs.at(noHardClipIdx)->bam, startRefPos_cns, endRefPos_cns);
											//queryseq = "=ACMGRSVTWYHKDBN"[bam_get_seq(query_aln_segs.at(noHardClipIdx)->bam)];
											break;
										case BAM_CIGAR_NO_DIFF_NO_MD:
										case BAM_CIGAR_DIFF_MD:
										case BAM_CIGAR_DIFF_NO_MD:
											//alnSegs = generateAlnSegs_no_MD(b, refseq, startPos, endPos);
											query_alnSegs = generateAlnSegs_no_MD2(query_aln_segs.at(noHardClipIdx)->bam, refseq, startRefPos_cns, endRefPos_cns);
											break;
										default:
											cerr << __func__ << ": unknown bam type, error!" << endl;
											exit(1);
									}
								}

	//							if(qname.compare("SRR8858450.1.35388")==0){
	//								cout << "-- qname=" << qname << endl;
	//							}

								//compute query_num
								//end_paf_pos = getEndRefPosAlnSeg(paf_alnseg->startRpos, paf_alnseg->opflag, paf_alnseg->seglen);
								end_seg_pos = getEndRefPosAlnSeg(query_alnSegs.at(query_alnSegs.size()-1)->startRpos, query_alnSegs.at(query_alnSegs.size()-1)->opflag, query_alnSegs.at(query_alnSegs.size()-1)->seglen);
								if(start_var_pos_para>query_alnSegs.at(0)->startRpos and end_var_pos_para<end_seg_pos){
								//if(paf_alnseg->startRpos>query_alnSegs.at(0)->startRpos and (paf_alnseg->startRpos+paf_alnseg->seglen)<(query_alnSegs.at(query_alnSegs.size()-1)->startRpos + query_alnSegs.at(query_alnSegs.size()-1)->seglen)){
									//cout << "query_num=" << query_num << ", qname=" << qname << ", i=" << i << endl;
									query_num ++;
								}

								//find qname
								q = find(clu_qname_vec.begin(), clu_qname_vec.end(), qname);
								if(q==clu_qname_vec.end()){ // not found
									clipAlnDataVector.at(i)->query_checked_flag = true;
								}else{ // found
									//search
									min_ref_dist = INT_MAX;
									for(j=0; j<query_alnSegs.size(); j++){
										seg = query_alnSegs.at(j);
										// recompute the overlap
										end_ref_pos = getEndRefPosAlnSeg(seg->startRpos, seg->opflag, seg->seglen);
										flag = false;
										if(seg->startRpos-1>=startRefPos_cns and end_ref_pos-1<=endRefPos_cns) flag = true;

										if(flag and opflag_para==seg->opflag){
											if(oplen_para>=min_sv_size and seg->seglen>=min_sv_size){
												if(oplen_para <= seg->seglen) size_ratio = (double)oplen_para / seg->seglen;
												else size_ratio = (double)seg->seglen / oplen_para;
												if(size_ratio>=size_ratio_match_thres){
													ref_dist = abs(start_var_pos_para - seg->startRpos);
													//if(ref_distance<200){
													if(ref_dist<=MAX_REF_DIST_SAME_CHR){
														//cout << "\tvar_num=" << var_num << ", ref_dist=" << ref_dist << ", seg.pos=" << seg->startRpos << ", seg.len=" << seg->seglen << ", qname=" << qname << ", startRpos=" << paf_alnseg->startRpos << ", seglen=" << paf_alnseg->seglen << ", i=" << i << ", j=" << j << endl;
														//var_num ++;
														if(ref_dist<min_ref_dist) min_ref_dist = ref_dist;
														//len_sum+=(*seg)->seglen;
													}
												}
												else{
													// cout << "opflag:" << opflag_para << "," << start_var_pos_para << "," << oplen_para << " seg:" << seg->startRpos << "," <<  seg->seglen << endl;
												}
											}
										}
									}
									if(min_ref_dist!=INT_MAX){
										var_num ++;
									}else{ // added on 2024-12-31, try to deal with the partial signatures
										//cout << qname << ": query_alnSegs.size=" << query_alnSegs.size() << ", paf_alnseg.start=" << paf_alnseg->startRpos << ", paf.seglen=" << paf_alnseg->seglen << endl;

										// compute the start_idx, end_idx, and then search partial variants around 1-kb regions
										start_ref_pos_check = start_var_pos_para - VAR_ALN_EXTEND_SIZE;
										end_ref_pos_check = end_var_pos_para + VAR_ALN_EXTEND_SIZE;
										if(start_ref_pos_check<1) start_ref_pos_check = 1;

										start_idx_check = 0;
										for(j=0; j<query_alnSegs.size(); j++){
											seg = query_alnSegs.at(j);
											if(seg->startRpos>=start_var_pos_para) {
												start_idx_check = j;
												break;
											}
										}

										left_check_idx_vec.clear();
										right_check_idx_vec.clear();
										left_flag = right_flag = true;
										num = sum_len = 0;
										while(start_idx_check-num>=0 or start_idx_check+num<(int64_t)query_alnSegs.size()){
											flag = false;
											// check left side
											idx = start_idx_check - num;
											if(idx>=0){
												seg = query_alnSegs.at(idx);
												if(left_flag and seg->opflag==opflag_para and seg->seglen>=min_sv_size){
													if(seg->startRpos>start_ref_pos_check){
														left_check_idx_vec.insert(left_check_idx_vec.begin(), idx);
														sum_len += seg->seglen;
														flag = true;
													}else left_flag = false;
												}
											}

											// check right side
											if(num>0){
												idx = start_idx_check + num;
												if(idx<(int64_t)query_alnSegs.size()){
													seg = query_alnSegs.at(idx);
													if(right_flag and seg->opflag==opflag_para and seg->seglen>=min_sv_size){
														if(seg->startRpos<end_ref_pos_check){
															right_check_idx_vec.insert(right_check_idx_vec.begin(), idx);
															sum_len += seg->seglen;
															flag = true;
														}else right_flag = false;
													}
												}
											}

											// validate the sum_len
											if(flag){
												if(oplen_para <= sum_len) size_ratio = (double)oplen_para / sum_len;
												else size_ratio = (double)sum_len / oplen_para;
												if(size_ratio>=size_ratio_match_thres_merge){
													var_num ++;
													break;
												}
											}

											num ++;
											if(left_flag==false and right_flag==false) break;
										}
									}

									//query_num+=1;
									//destroyAlnSegs(query_alnSegs);
								}
							}
						}else{//all hard clip
							for(j=0; j<query_aln_segs.size(); j++) query_aln_segs.at(j)->query_checked_flag = true;
						}

						destroyAlnSegs(query_alnSegs);
					}
				}
			}
		}
	}

	if(!clipAlnDataVector.empty()) destoryClipAlnData(clipAlnDataVector);

	supp_num_vec.push_back(var_num);
	supp_num_vec.push_back(query_num);

	return supp_num_vec;
}

// refine allele frequency
void varCand::refineAlleleFreq(vector<reg_t*> &var_vec){
	bool freq_ok_flag;

	freq_ok_flag = isAlleleFreqOk(var_vec);

	if(freq_ok_flag==false) refineAlleleFreqOp(var_vec);
}

bool varCand::isAlleleFreqOk(vector<reg_t*> &var_vec){
	bool freq_ok_flag;
	size_t i;
	reg_t *reg;
	set<int32_t> ctg_id_vec;

	freq_ok_flag = false;
	// OK cases: alleles in different contigs
	if(var_vec.size()>=2){
		// deal with items of rescue_flag==false
		for(i=0; i<var_vec.size(); i++){
			reg = var_vec.at(i);
			if(reg->rescue_flag==false and ctg_id_vec.find(reg->query_id)==ctg_id_vec.end()){
				ctg_id_vec.insert(reg->query_id);
			}
		}

		if(ctg_id_vec.size()>=2){
			freq_ok_flag = true;
		}
	}

	return freq_ok_flag;
}

// allele frequency refine operation
void varCand::refineAlleleFreqOp(vector<reg_t*> &var_vec){
	size_t i, j, k, ki, start_idx;
	reg_t *reg, *reg_tmp;
	int32_t target_query_id, seq_len, max_merge_span, DP_num, sum, m, idx, dist1;
	vector<clipAlnData_t*> clipAlnDataVector;
	int64_t end_ref_pos, end_query_pos, start_var_pos, end_var_pos, start_var_pos_pure0, end_var_pos_pure0, start_var_pos_pure, end_var_pos_pure, start_var_pos_check, end_var_pos_check, chrlen_tmp;
	int64_t start_qpos, end_qpos, start_qpos_reg, end_qpos_reg;
	uint8_t *seq_int;
	vector<struct querySeqInfoNode*> query_seq_info_all;
	struct querySeqInfoNode *qnode;
	vector<struct alnSeg*> query_alnSegs;
	struct alnSeg *seg, *clip_seg;
	vector<qcSig_t*> overlap_qcSig_vec;
	qcSig_t *qcSig;
	string reg_str, refseq, cnsseq, refbase, queryseq_comp, varseq_comp;
	char *p_seq;
	bool flag, seq_valid_flag;
	map<int32_t, vector<int32_t> > occ_map;
	set<string> flanking_qname_set, success_qname_set;
	string gt_header, gt_str, dp_str, ad_str1, ad_str2;

	int32_t item_num, minVal;
	double min_sv_size_factor, len_ratio, sim_val;

	chrlen_tmp = faidx_seq_len64(fai, chrname.c_str());

	start_var_pos_pure0 = varVec.at(0)->startRefPos;
	end_var_pos_pure0 = varVec.at(varVec.size()-1)->endRefPos;
	start_var_pos_pure = start_var_pos_pure0 - SHORT_VAR_ALN_CHECK_EXTEND_SIZE;
	end_var_pos_pure = end_var_pos_pure0 + SHORT_VAR_ALN_CHECK_EXTEND_SIZE;
	if(start_var_pos_pure<1) start_var_pos_pure = 1;
	if(end_var_pos_pure>chrlen_tmp) end_var_pos_pure = chrlen_tmp;

	start_var_pos = start_var_pos_pure - VAR_ALN_EXTEND_SIZE; // EXT_SIZE_CHK_VAR_LOC;
	end_var_pos = end_var_pos_pure + VAR_ALN_EXTEND_SIZE; // EXT_SIZE_CHK_VAR_LOC;
	if(start_var_pos<1) start_var_pos = 1;
	if(end_var_pos>chrlen_tmp) end_var_pos = chrlen_tmp;

	max_merge_span = min_distance_merge; // CNS_MIN_DIST_MERGE_THRES;
//	if(sv_len_est>CNS_MIN_DIST_MERGE_THRES) max_merge_span = 1.5 * CNS_MIN_DIST_MERGE_THRES;

	seq_valid_flag = false;
	minVal = INT_MAX;
	for(i=0; i<var_vec.size(); i++){
		reg = var_vec.at(i);
		if(reg->call_success_status and abs(reg->sv_len)<minVal) minVal = abs(reg->sv_len);
		if(seq_valid_flag==false and reg->altseq.size()>0 and reg->refseq.size()>0) seq_valid_flag = true;
	}
	if(seq_valid_flag==false) return; // skip

	min_sv_size_factor = 0.98;
	if(minVal<SHORT_VAR_ALN_CHECK_EXTEND_SIZE) min_sv_size_factor = ULTRA_SHORT_SV_SIZE_FACTOR;

#if GT_REFINE_DEBUG
	cout << "min_sv_size_factor=" << min_sv_size_factor << endl;
#endif

	// load the clipping data
	clipAlnDataLoader data_loader(chrname, start_var_pos, end_var_pos, inBamFile, minClipEndSize, maxVarRegSize, minMapQ, minHighMapQ);
	if(clip_reg_flag) data_loader.loadClipAlnDataWithSATag(clipAlnDataVector, max_ultra_high_cov);
	else data_loader.loadClipAlnDataWithSATagWithSegSize(clipAlnDataVector, max_ultra_high_cov, max_seg_size_ratio, max_seg_nm_ratio, fai, max_absig_density);

	reg_str = chrname + ":" + to_string(start_var_pos) + "-" + to_string(end_var_pos);
	pthread_mutex_lock(&mutex_fai);
	p_seq = fai_fetch(fai, reg_str.c_str(), &seq_len);
	pthread_mutex_unlock(&mutex_fai);
	refseq = p_seq;
	free(p_seq);

	// extract queries from clip align data vector
	query_seq_info_all = extractQueriesFromClipAlnDataVec(clipAlnDataVector, refseq, chrname, start_var_pos, end_var_pos, start_var_pos_pure, end_var_pos_pure, fai, minConReadLen, clip_reg_flag, &mutex_fai);

	for(i=0; i<query_seq_info_all.size(); i++){
		qnode = query_seq_info_all.at(i);
		query_alnSegs = qnode->query_alnSegs;
		seg = query_alnSegs.at(query_alnSegs.size()-1);
		end_ref_pos = getEndRefPosAlnSeg(seg->startRpos, seg->opflag, seg->seglen);
		if(query_alnSegs.at(0)->startRpos<=start_var_pos_pure0+SHORT_VAR_ALN_CHECK_EXTEND_SIZE and end_ref_pos>=end_var_pos_pure0-SHORT_VAR_ALN_CHECK_EXTEND_SIZE) { qnode->entire_flanking_flag = true; flanking_qname_set.insert(qnode->qname); }
		else qnode->entire_flanking_flag = false;

		qnode->qcSig_vec = extractQcSigsFromAlnSegsSingleQuery(qnode, qnode->clip_aln->chrname, start_var_pos, end_var_pos, min_sv_size*min_sv_size_factor);
		// remove neighboring false signatures
		rmNeighborFalseSig(qnode->qcSig_vec, MIN_AVER_SIZE_ALN_SEG, MIN_SIZE_RATIO_MATCH_CLIP_POS);
		qnode->qcSig_merge_flag = mergeNeighbouringSigsFlag(qnode, qnode->qcSig_vec, MAX_DIST_MERGE_ARBITARY, max_merge_span, CALL_MIN_SEQSIM_MERGE_THRES2, MIN_VALID_SIG_SIZE_RATIO_THRES, fai);
		filterSmallSigs(qnode->qcSig_vec, MIN_SIG_SIZE_RATIO_FILTER_THRES, MAX_INNER_MISSING_IGNORE_SIZE);
	}
	DP_num = flanking_qname_set.size();

#if GT_REFINE_DEBUG
	cout << "DP_num=" << DP_num << endl;
	printQcSigIndel(query_seq_info_all);
#endif

	target_query_id = -1;
	cnsseq = "";
	for(i=0; i<var_vec.size(); i++){
		reg = var_vec.at(i);
		if(reg->call_success_status==true and reg->altseq.size()>0 and reg->refseq.size()>0){
			success_qname_set.clear();

			start_var_pos_check = reg->startRefPos - EXT_SIZE_CHK_VAR_LOC;
			end_var_pos_check = reg->endRefPos + EXT_SIZE_CHK_VAR_LOC;
			if(start_var_pos_check<1) start_var_pos_check = 1;
			if(end_var_pos_check>chrlen_tmp) end_var_pos_check = chrlen_tmp;

			// get consensus sequence
			if(reg->query_id!=-1 and reg->query_id!=target_query_id){
				if(reg->rescue_flag){
					FastaSeqLoader fa_loader(rescue_cnsfilename);
					cnsseq = fa_loader.getFastaSeq(reg->query_id, reg->aln_orient);
				}else{
					FastaSeqLoader fa_loader(ctgfilename);
					cnsseq = fa_loader.getFastaSeq(reg->query_id, reg->aln_orient);
				}
				target_query_id = reg->query_id;
			}else if(reg->query_id==-1 and cnsseq.size()>0){
				cnsseq = "";
				target_query_id = -1;
			}

			// try size sum
			item_num = 0;
			for(j=0; j<query_seq_info_all.size(); j++){
				qnode = query_seq_info_all.at(j);
				if(success_qname_set.find(qnode->qname)!=success_qname_set.end()) continue;

//				if(qnode->qname.compare("m84039_230414_235240_s2/96536133/ccs")==0){
//					cout << "j=" << j << ", qname=" << qnode->qname << endl;
//				}

				seq_int = bam_get_seq(qnode->clip_aln->bam);
				flag = false;
				for(start_idx=0; start_idx<qnode->qcSig_vec.size(); start_idx++){
					sum = 0;
					flag = false;
					for(k=start_idx; k<qnode->qcSig_vec.size(); k++){
						qcSig = qnode->qcSig_vec.at(k);
						end_ref_pos = getEndRefPosAlnSeg(qcSig->ref_pos, qcSig->cigar_op, qcSig->cigar_op_len);
						if(qcSig->ref_pos>=start_var_pos_check and end_ref_pos<=end_var_pos_check){
							if(qcSig->cigar_op==reg->var_type or (qcSig->cigar_op==BAM_CINS and reg->var_type==VAR_DUP)){
								sum += qcSig->cigar_op_len;
								if(abs(reg->sv_len)<sum) len_ratio = (double)abs(reg->sv_len) / sum;
								else len_ratio = (double)sum / abs(reg->sv_len);
								if(len_ratio>=QC_SIZE_RATIO_MATCH_THRES_INDEL_MERGE*min_sv_size_factor){
									item_num ++;
									flag = true;
									success_qname_set.insert(qnode->qname);
									break;
								}else if(sum>abs(reg->sv_len)) break;
							}
						}
					}
					if(flag) break;
				}

				// consider the larger one in the reads, try to check the merge counterparts
				if(flag==false){
					for(k=0; k<qnode->qcSig_vec.size(); k++){
						qcSig = qnode->qcSig_vec.at(k);
						end_ref_pos = getEndRefPosAlnSeg(qcSig->ref_pos, qcSig->cigar_op, qcSig->cigar_op_len);
						if(qcSig->ref_pos>=start_var_pos_check and end_ref_pos<=end_var_pos_check){
							if(qcSig->cigar_op==reg->var_type or (qcSig->cigar_op==BAM_CINS and reg->var_type==VAR_DUP)){
								if(qcSig->merge_flag and qcSig->cigar_op_len>abs(reg->sv_len)){
									for(start_idx=0; start_idx<var_vec.size(); start_idx++){
										sum = 0;
										for(ki=start_idx; ki<var_vec.size(); ki++){
											reg_tmp = var_vec.at(ki);
											if(reg_tmp->call_success_status==true){
												if(reg->var_type==reg_tmp->var_type or (reg->var_type==VAR_INS and reg_tmp->var_type==VAR_DUP)){
													if(reg_tmp->startRefPos>=start_var_pos_check and reg_tmp->endRefPos<=end_var_pos_check){
														sum += abs(reg_tmp->sv_len);
														if(qcSig->cigar_op_len<sum) len_ratio = (double)qcSig->cigar_op_len / sum;
														else len_ratio = (double)sum / qcSig->cigar_op_len;
														if(len_ratio>=QC_SIZE_RATIO_MATCH_THRES_INDEL_MERGE*min_sv_size_factor){
															if(i>=start_idx and i<=ki){
																item_num ++;
																flag = true;
																success_qname_set.insert(qnode->qname);
																break;
															}
														}else if(sum>qcSig->cigar_op_len) break;
													}
												}
											}
										}
										if(flag) break;
									}
									if(flag) break;
								}
							}
						}
					}
				}

				// deal with clipping ends
				if(flag==false and (reg->var_type==VAR_INS or reg->var_type==VAR_DUP) and reg->query_id!=-1 and abs(reg->sv_len)>=minClipEndSize){ // only consider INS, what about DEL?
					// try left end overlap
					clip_seg = qnode->query_alnSegs.at(0);
					if((clip_seg->opflag==BAM_CSOFT_CLIP or clip_seg->opflag==BAM_CHARD_CLIP) and clip_seg->seglen>=minClipEndSize and clip_seg->startRpos>=start_var_pos_pure and clip_seg->startRpos<=end_var_pos_pure){ // maybe there are errors

						// get start query pos
						if(clip_seg->opflag==BAM_CSOFT_CLIP) start_qpos = 1;
						else start_qpos = clip_seg->startQpos;

						idx = -1;
						for(m=0; m<(int32_t)qnode->query_alnSegs.size(); m++){
							seg = qnode->query_alnSegs.at(m);
							if(seg->opflag==BAM_CMATCH or seg->opflag==BAM_CEQUAL){
								end_ref_pos = getEndRefPosAlnSeg(seg->startRpos, seg->opflag, seg->seglen);
								if(reg->endRefPos+SHORT_VAR_ALN_CHECK_EXTEND_SIZE>=seg->startRpos and reg->endRefPos<=end_ref_pos+SHORT_VAR_ALN_CHECK_EXTEND_SIZE){
									idx = m;
									break;
								}else if(seg->startRpos>reg->endRefPos+EXT_SIZE_CHK_VAR_LOC) break;
							}
						}

						if(idx!=-1){
							seg = qnode->query_alnSegs.at(idx);
							end_ref_pos = getEndRefPosAlnSeg(seg->startRpos, seg->opflag, seg->seglen);
							end_query_pos = getEndQueryPosAlnSeg(seg->startQpos, seg->opflag, seg->seglen);
							dist1 = end_ref_pos - reg->endRefPos;
							end_qpos = end_query_pos - dist1;
							seq_len = end_qpos - start_qpos + 1;
							if(seq_len>0){
								//end_qpos_reg = reg->endQueryPos;
								end_qpos_reg = reg->endQueryPos;
								if(seq_len>=abs(reg->sv_len)) seq_len = abs(reg->sv_len);
								if(end_qpos_reg<seq_len) seq_len = end_qpos_reg;
								start_qpos = end_qpos - seq_len + 1;
								start_qpos_reg = end_qpos_reg - seq_len + 1;

								// compare to variant sequence
								queryseq_comp = varseq_comp = "";
								if(cnsseq.size()>0){
									for(m=start_qpos-1; m<end_qpos; m++) queryseq_comp += "=ACMGRSVTWYHKDBN"[bam_seqi(seq_int, m)];  // seq
									//varseq_comp = reg->altseq.substr(start_qpos_reg-reg->startQueryPos, seq_len);
									varseq_comp = cnsseq.substr(start_qpos_reg, seq_len);
								}else if(start_qpos_reg>=reg->startQueryPos){
									for(m=start_qpos-1; m<end_qpos; m++) queryseq_comp += "=ACMGRSVTWYHKDBN"[bam_seqi(seq_int, m)];  // seq
									varseq_comp = reg->altseq.substr(start_qpos_reg-reg->startQueryPos, seq_len);
									//varseq_comp = cnsseq.substr(start_qpos_reg, seq_len);
								}

								if(queryseq_comp.size()>0 and varseq_comp.size()>0){
									sim_val = computeVarseqSim(queryseq_comp, varseq_comp);

#if GT_REFINE_DEBUG
									cout << "LEFT END: queryseq_comp=" << queryseq_comp << ", size=" << queryseq_comp.size() << endl;
									cout << "LEFT END:   varseq_comp=" << varseq_comp << ", size=" << varseq_comp.size() << endl;
									cout << "LEFT END:       sim_val=" << sim_val << endl;
#endif

									// confirm
									if(sim_val>=min_seqsim_match) {
										item_num ++;
										flag = true;
										success_qname_set.insert(qnode->qname);
									}
								}
							}
						}
					}

					// try right end overlap
					if(flag==false){
						clip_seg = qnode->query_alnSegs.at(qnode->query_alnSegs.size()-1);
						if((clip_seg->opflag==BAM_CSOFT_CLIP or clip_seg->opflag==BAM_CHARD_CLIP) and clip_seg->seglen>=minClipEndSize and clip_seg->startRpos>=start_var_pos_pure and clip_seg->startRpos<=end_var_pos_pure){ // maybe there are errors
							// get query sequence
							// get start query pos
							end_qpos = getEndQueryPosAlnSeg(clip_seg->startQpos, clip_seg->opflag, clip_seg->seglen);
							idx = -1;
							for(m=qnode->query_alnSegs.size()-1; m>=0; m--){
								seg = qnode->query_alnSegs.at(m);
								if(seg->opflag==BAM_CMATCH or seg->opflag==BAM_CEQUAL){
									end_ref_pos = getEndRefPosAlnSeg(seg->startRpos, seg->opflag, seg->seglen);
									if(reg->startRefPos+SHORT_VAR_ALN_CHECK_EXTEND_SIZE>=seg->startRpos and reg->startRefPos<=end_ref_pos+SHORT_VAR_ALN_CHECK_EXTEND_SIZE){
										idx = m;
										break;
									}else if(end_ref_pos+EXT_SIZE_CHK_VAR_LOC<reg->startRefPos) break;
								}
							}
							if(idx!=-1){
								seg = qnode->query_alnSegs.at(idx);
								end_ref_pos = getEndRefPosAlnSeg(seg->startRpos, seg->opflag, seg->seglen);
								dist1 = reg->startRefPos - seg->startRpos;
								start_qpos = seg->startQpos + dist1;
								seq_len = end_qpos - start_qpos + 1;
								if(seq_len>0){
									start_qpos_reg = reg->startQueryPos;
									if(seq_len>=abs(reg->sv_len)) seq_len = abs(reg->sv_len);
									if(start_qpos_reg<seq_len) seq_len = start_qpos_reg;
									end_qpos = start_qpos + seq_len - 1;
									end_qpos_reg = start_qpos_reg + seq_len - 1;

									// compare to variant sequence
									queryseq_comp = varseq_comp = "";
									if(cnsseq.size()>0){
										for(m=start_qpos-1; m<end_qpos; m++) queryseq_comp += "=ACMGRSVTWYHKDBN"[bam_seqi(seq_int, m)];  // seq
										//varseq_comp = reg->altseq.substr(start_qpos_reg-reg->startQueryPos, seq_len);
										varseq_comp = cnsseq.substr(start_qpos_reg, seq_len);
									}else if(start_qpos_reg>=reg->startQueryPos){
										for(m=start_qpos-1; m<end_qpos; m++) queryseq_comp += "=ACMGRSVTWYHKDBN"[bam_seqi(seq_int, m)];  // seq
										varseq_comp = reg->altseq.substr(start_qpos_reg-reg->startQueryPos, seq_len);
									}

									if(queryseq_comp.size()>0 and varseq_comp.size()>0){
										sim_val = computeVarseqSim(queryseq_comp, varseq_comp);

#if GT_REFINE_DEBUG
										cout << "RIGHT END: queryseq_comp=" << queryseq_comp << ", size=" << queryseq_comp.size() << endl;
										cout << "RIGHT END:   varseq_comp=" << varseq_comp << ", size=" << varseq_comp.size() << endl;
										cout << "RIGHT END:       sim_val=" << sim_val << endl;
#endif

										// confirm
										if(sim_val>=min_seqsim_match) {
											item_num ++;
											flag = true;
											success_qname_set.insert(qnode->qname);
										}
									}
								}
							}
						}
					}
				}

#if GT_REFINE_DEBUG
				cout << "j=" << j << ", qname=" << qnode->qname << ", start_idx=" << start_idx << ", item_num=" << item_num << endl;
#endif

			}

#if GT_REFINE_DEBUG
			cout << "item_num=" << item_num << ", DP_num=" << DP_num << "; AD=" << reg->supp_num << ", DP=" << reg->DP << ", AF=" << reg->AF << endl;
#endif

			if(item_num>DP_num) DP_num = item_num;
			if(DP_num>0 and (reg->DP==0 or reg->supp_num==0 or (double)item_num/DP_num>(double)reg->supp_num/reg->DP)){ // update
				reg->supp_num = item_num;
				reg->DP = DP_num;
				reg->AF = (double)reg->supp_num / reg->DP;
				if(reg->AF>=gt_homo_ratio_thres){
					reg->gt_type = GT_HOMOZYGOUS;
				}else if(reg->AF>=gt_hete_ratio_thres){
					reg->gt_type = GT_HETEROZYGOUS;
				}else{
					reg->gt_type = GT_NOZYGOUS;
				}

				gt_header = "GT:AD:DP";
				gt_str = "./.";
				ad_str1 = ".";
				ad_str2 = ".";
				dp_str = ".";

				if(reg->gt_type==GT_HOMOZYGOUS){
					gt_str = GT_HOMOZYGOUS_STR;
				}else if(reg->gt_type==GT_HETEROZYGOUS){
					gt_str = GT_HETEROZYGOUS_STR;
				}else if(reg->gt_type==GT_NOZYGOUS){
					gt_str = GT_NOZYGOUS_STR;
				}else {
					cout << "line=" << __LINE__ << ", unknown genotype=" << reg->gt_type << ", error!" << endl;
					exit(1);
				}
				ad_str1 = to_string(reg->DP-reg->supp_num);
				ad_str2 = to_string(reg->supp_num);
				dp_str = to_string(reg->DP);

				reg->gt_seq = gt_header + "\t" + gt_str + ":" + ad_str1 + "," + ad_str2 + ":" + dp_str;
			}
		}
	}

	// release memory
	if(!clipAlnDataVector.empty()) destoryClipAlnData(clipAlnDataVector);
	if(!query_seq_info_all.empty()) destoryQuerySeqInfoAll(query_seq_info_all);

}

// compute the genotypes of indel regions
void varCand::computeGenotypeIndelReg(vector<reg_t*> &var_vec){

	struct alleleNode{
		reg_t *reg;
		int32_t supp_num, depth : 27, gt_type: 5;
		double AF_val;
		struct alleleNode *mate_allele;
		bool valid_flag;
	};

	size_t  i, j;
	string gt_header, gt_str, dp_str, ad_str1, ad_str2;
	bool flag, overlap_flag;

	vector<struct alleleNode*> allele_vec;
	struct alleleNode *allele_node, *allele_node2;
	double val;
	string comp_seq1, comp_seq2, side_seq, reg_str, chrname_tmp;
	int64_t min_pos, max_pos;
	char *p_seq;
	int32_t seq_len;

	if(var_vec.size()>=0){
		for(i=0; i<var_vec.size(); i++){
			allele_node = new struct alleleNode();
			allele_node->reg = var_vec.at(i);
			allele_node->mate_allele = NULL;
			allele_node->supp_num = var_vec.at(i)->supp_num;
			allele_node->depth = var_vec.at(i)->DP;
			allele_node->valid_flag = true;
			allele_vec.push_back(allele_node);
		}
	}

//	reg_t *reg;
//	cout << "\talleles: ";
//	for(i=0; i<allele_vec.size(); i++){
//		reg = allele_vec.at(i)->reg;
//		cout << "[" << reg->chrname << ":" << reg->startRefPos << "-" << reg->endRefPos << ", " << to_string(reg->var_type) << ", " << reg->sv_len << "], ";
//	}
//	cout << endl;
//
//	cout << "\tvarVec: ";
//	for(i=0; i<varVec.size(); i++){
//		reg = varVec.at(i);
//		cout << "[" << reg->chrname << ":" << reg->startRefPos << "-" << reg->endRefPos << ", " << to_string(reg->var_type) << ", " << reg->sv_len << "], ";
//	}
//	cout << endl;

	// genotyping
	// get the mate allele
	for(i=0; i<allele_vec.size(); i++){
		allele_node = allele_vec.at(i);
		flag = false;
		if(allele_node->mate_allele==NULL){
			for(j=i+1; j<allele_vec.size(); j++){
				allele_node2 = allele_vec.at(j);
				if(allele_node->valid_flag and allele_node2->valid_flag){
					if(allele_node->reg->query_id!=allele_node2->reg->query_id){ // on different queries
						//overlap_flag = isOverlappedRegExtSize(allele_node->reg, allele_node2->reg, 5, 5); // deleted on 2024-10-01
						overlap_flag = isOverlappedRegExtSize(allele_node->reg, allele_node2->reg, 100, 100);
						if(overlap_flag){
							if(allele_node->reg->var_type==allele_node2->reg->var_type){ // same variant type
								val = 0;
								comp_seq1 = comp_seq2 = "";
								if(allele_node->reg->var_type==VAR_INS){ // INS
									// left side part
									min_pos = max_pos = -1;
									if(allele_node->reg->startRefPos!=allele_node2->reg->startRefPos){
										if(allele_node->reg->startRefPos<allele_node2->reg->startRefPos){
											min_pos = allele_node->reg->startRefPos;
											max_pos = allele_node2->reg->startRefPos - 1;
											chrname_tmp = allele_node2->reg->chrname;
										}else{
											min_pos = allele_node2->reg->startRefPos;
											max_pos = allele_node->reg->startRefPos - 1;
											chrname_tmp = allele_node->reg->chrname;
										}
										reg_str = chrname_tmp + ":" + to_string(min_pos) + "-" + to_string(max_pos);
										pthread_mutex_lock(&mutex_fai);
										p_seq = fai_fetch(fai, reg_str.c_str(), &seq_len);
										pthread_mutex_unlock(&mutex_fai);
										side_seq = p_seq;
										free(p_seq);

										if(allele_node->reg->startRefPos<allele_node2->reg->startRefPos)
											comp_seq2 += side_seq;
										else if(allele_node->reg->startRefPos>allele_node2->reg->startRefPos)
											comp_seq1 += side_seq;
									}

									// variant part
									comp_seq1 += allele_node->reg->altseq;
									comp_seq2 += allele_node2->reg->altseq;

									// right side part
									if(allele_node->reg->endRefPos!=allele_node2->reg->endRefPos){
										if(allele_node->reg->endRefPos<allele_node2->reg->endRefPos){
											min_pos = allele_node->reg->endRefPos + 1;
											max_pos = allele_node2->reg->endRefPos;
											chrname_tmp = allele_node->reg->chrname;
										}else{
											min_pos = allele_node2->reg->endRefPos + 1;
											max_pos = allele_node->reg->endRefPos;
											chrname_tmp = allele_node2->reg->chrname;
										}
										reg_str = chrname_tmp + ":" + to_string(min_pos) + "-" + to_string(max_pos);
										pthread_mutex_lock(&mutex_fai);
										p_seq = fai_fetch(fai, reg_str.c_str(), &seq_len);
										pthread_mutex_unlock(&mutex_fai);
										side_seq = p_seq;
										free(p_seq);

										if(allele_node->reg->endRefPos<allele_node2->reg->endRefPos)
											comp_seq1 += side_seq;
										else
											comp_seq2 += side_seq;
									}

									// compute similarity
									//val = computeVarseqSim(allele_node->reg->altseq, allele_node2->reg->altseq); // ?
									val = computeVarseqSim(comp_seq1, comp_seq2);
									//cout << "val=" << val << ", comp_seq1.size=" << comp_seq1.size() << ", comp_seq2.size=" << comp_seq2.size() << ", altseq1.size=" << allele_node->reg->altseq.size() << ", altseq2.size=" << allele_node2->reg->altseq.size() << endl;
								}else{ // DEL
									if(allele_node->reg->refseq.size()<=allele_node2->reg->refseq.size()) val = (double)allele_node->reg->refseq.size() / allele_node2->reg->refseq.size();
									else val = (double)allele_node2->reg->refseq.size() / allele_node->reg->refseq.size();
								}
								if(val<gt_min_seqsim_merge){ // heterozygous
									allele_node->mate_allele = allele_node2;
									allele_node2->mate_allele = allele_node;
									allele_node->reg->gt_type = GT_HETEROZYGOUS;
									allele_node2->reg->gt_type = GT_HETEROZYGOUS;
									allele_node->gt_type = GT_HETEROZYGOUS;
									allele_node2->gt_type = GT_HETEROZYGOUS;
									allele_node->AF_val = (double)allele_node->supp_num / allele_node->depth;
									allele_node2->AF_val = (double)allele_node2->supp_num / allele_node2->depth;
									allele_node->reg->DP = allele_node->depth;
									allele_node2->reg->DP = allele_node2->depth;
									allele_node->reg->AF = allele_node->AF_val;
									allele_node2->reg->AF = allele_node2->AF_val;
								}else{ // homozygous, merge the same variant, prefer the larger supp_num
									if(allele_node->supp_num>=allele_node2->supp_num){
										allele_node->reg->gt_type = GT_HOMOZYGOUS;
										allele_node->gt_type = GT_HOMOZYGOUS;
										allele_node->supp_num += allele_node2->supp_num;
										allele_node->depth = max(allele_node->supp_num, allele_node->depth);
										allele_node->AF_val = (double)allele_node->supp_num / allele_node->depth;
										allele_node->reg->supp_num = allele_node->supp_num;
										allele_node->reg->DP = allele_node->depth;
										allele_node->reg->AF = allele_node->AF_val;
										allele_node2->reg->call_success_status = false;
										//delete allele_node2;
										//allele_vec.erase(allele_vec.begin()+j);
										allele_node2->valid_flag = false;
									}else{
										allele_node2->reg->gt_type = GT_HOMOZYGOUS;
										allele_node2->gt_type = GT_HOMOZYGOUS;
										allele_node2->supp_num += allele_node->supp_num;
										allele_node2->depth = max(allele_node2->supp_num, allele_node2->depth);
										allele_node2->AF_val = (double)allele_node2->supp_num / allele_node2->depth;
										allele_node2->reg->supp_num = allele_node2->supp_num;
										allele_node2->reg->DP = allele_node2->depth;
										allele_node2->reg->AF = allele_node2->AF_val;
										allele_node->reg->call_success_status = false;
										//delete allele_node;
										//allele_vec.erase(allele_vec.begin()+i);
										allele_node->valid_flag = false;
									}
									flag = true;
								}
							}else{ // different variant type, heterozygous
								allele_node->mate_allele = allele_node2;
								allele_node2->mate_allele = allele_node;
								allele_node->reg->gt_type = GT_HETEROZYGOUS;
								allele_node2->reg->gt_type = GT_HETEROZYGOUS;
								allele_node->gt_type = GT_HETEROZYGOUS;
								allele_node2->gt_type = GT_HETEROZYGOUS;
								allele_node->AF_val = (double)allele_node->supp_num / allele_node->depth;
								allele_node2->AF_val = (double)allele_node2->supp_num / allele_node2->depth;
								allele_node->reg->DP = allele_node->depth;
								allele_node2->reg->DP = allele_node2->depth;
								allele_node->reg->AF = allele_node->AF_val;
								allele_node2->reg->AF = allele_node2->AF_val;
							}
							if(flag==true) break;
							//break;
						}
					}
				}
			}

			// no mate-allele, needs to genotyping according to Nsupp and Depth
			if(flag==false){ //
				val = (double)allele_node->supp_num / allele_node->depth;
				if(val>=gt_homo_ratio_thres){
					allele_node->reg->gt_type = GT_HOMOZYGOUS;
					allele_node->gt_type = GT_HOMOZYGOUS;
				}else if(val>=gt_hete_ratio_thres){
					allele_node->reg->gt_type = GT_HETEROZYGOUS;
					allele_node->gt_type = GT_HETEROZYGOUS;
				}else{
					allele_node->reg->gt_type = GT_NOZYGOUS;
					allele_node->gt_type = GT_NOZYGOUS;
				}
				allele_node->reg->DP = allele_node->depth;
				allele_node->AF_val = allele_node->reg->AF = val;
			}
		}
	}

	// compute the genotype
	for(i=0; i<allele_vec.size(); i++){
		allele_node = allele_vec.at(i);
		if(allele_node->valid_flag){
			val = (double)allele_node->supp_num / allele_node->depth;
			if(val>=gt_homo_ratio_thres){
				allele_node->reg->gt_type = GT_HOMOZYGOUS;
				allele_node->gt_type = GT_HOMOZYGOUS;
			}else if(val>=gt_hete_ratio_thres){
				allele_node->reg->gt_type = GT_HETEROZYGOUS;
				allele_node->gt_type = GT_HETEROZYGOUS;
			}else{
				allele_node->reg->gt_type = GT_NOZYGOUS;
				allele_node->gt_type = GT_NOZYGOUS;
			}
			allele_node->reg->DP = allele_node->depth;
			allele_node->AF_val = allele_node->reg->AF = val;

			gt_header = "GT:AD:DP";
			gt_str = "./.";
			ad_str1 = ".";
			ad_str2 = ".";
			dp_str = ".";

			if(allele_node->gt_type==GT_HOMOZYGOUS){
				gt_str = GT_HOMOZYGOUS_STR;
			}else if(allele_node->gt_type==GT_HETEROZYGOUS){
				gt_str = GT_HETEROZYGOUS_STR;
			}else if(allele_node->gt_type==GT_NOZYGOUS){
				gt_str = GT_NOZYGOUS_STR;
			}else {
				cout << "line=" << __LINE__ << ", unknown genotype=" << allele_node->gt_type << ", error!" << endl;
				exit(1);
			}
			ad_str1 = to_string(allele_node->depth-allele_node->supp_num);
			ad_str2 = to_string(allele_node->supp_num);
			dp_str = to_string(allele_node->depth);

			allele_node->reg->gt_seq = gt_header + "\t" + gt_str + ":" + ad_str1 + "," + ad_str2 + ":" + dp_str;

			newVarVec.push_back(allele_node->reg);
		}
	}

	// release memory
	for(i=0; i<allele_vec.size(); i++) delete allele_vec.at(i);
	vector<struct alleleNode*>().swap(allele_vec);
}

vector<reg_t*> varCand::computeClipRegVarLoc(string &alnfilename, string &cnsfilename, string &refseqfilename, string &clusterfilename, vector<minimap2_aln_t*> &minimap2_aln_vec, vector<int32_t> *clusterId_incomplete, bool rescue_flag){
	vector<reg_t*> var_vec;
	vector<vector<string>> cluster_querynames_vec;
	vector<string> queryname_vec;
	size_t  i, cluster_id;
	string ref_seq, query_seq, tmp_refseq, tmp_queryseq;
	minimap2_aln_t *minimap2_aln, *minimap2_aln_left, *minimap2_aln_right;
	struct pafalnSeg *paf_alnseg, *paf_alnseg2;
	int64_t j, start_aln_seg_pos, end_aln_seg_pos, end_ref_pos, begin_var_idx, end_var_idx, begin_var_idx_new, end_var_idx_new, sv_len, dist_begin, dist_end, dist_query, dist_subject, dist_subject_tmp;
	int64_t start_var_ref_pos, end_var_ref_pos, start_var_subject_pos, end_var_subject_pos, start_var_query_pos, end_var_query_pos, end_qpos_tmp;
	reg_t *clip_reg_tmp;
	vector<int32_t> minimap2_item_id_vec, sidemost_item_id_vec, supp_cov_vec;
	bool tmp_call_success_flag, entire_flanking_flag, single_seg_flag, valid_sig_flag;
	int32_t query_id, idx_alnSegs_left, idx_alnSegs_right, margin_dist1, margin_dist2, new_sv_type, overlap_size, ref_dist, aln_orient, maxIdx, var_ext_size;
	int64_t left_rpos, left_subjectpos, left_qpos, right_rpos, right_subjectpos, right_qpos;
	double val, val2, maxVal;

	if(sv_type!=VAR_TRA){
		//var_ext_size = EXT_SIZE_CHK_VAR_LOC;
		var_ext_size = VAR_ALN_EXTEND_SIZE;

		// load query
		FastaSeqLoader fa_loader(cnsfilename);

		// load reference
		FastaSeqLoader ref_loader(refseqfilename);
		ref_seq = ref_loader.getFastaSeq(0);

		cluster_querynames_vec = getClusterInfo(clusterfilename);

		for(cluster_id=0; cluster_id<cluster_querynames_vec.size(); cluster_id++){
			if(clusterId_incomplete!=NULL and find((*clusterId_incomplete).begin(), (*clusterId_incomplete).end(), cluster_id) == (*clusterId_incomplete).end()){  // call complete, then skip
				continue;
			}

			queryname_vec = cluster_querynames_vec.at(cluster_id);
			minimap2_item_id_vec = getMinimapItemIdVec(queryname_vec, minimap2_aln_vec);

			// search for inverted align type
			tmp_call_success_flag = false;
			if(sv_type==VAR_INV){
				for(i=0; i<minimap2_item_id_vec.size(); i++){
					minimap2_aln = minimap2_aln_vec.at(minimap2_item_id_vec.at(i));
					if(minimap2_aln->type.compare("i")==0){
						start_aln_seg_pos = minimap2_aln->pafalnsegs.at(0)->startRpos;
						paf_alnseg = minimap2_aln->pafalnsegs.at(minimap2_aln->pafalnsegs.size()-1);
						end_aln_seg_pos = getEndRefPosAlnSeg(paf_alnseg->startRpos, paf_alnseg->opflag, paf_alnseg->seglen);
						overlap_size = getOverlapSize(leftClipRefPos, rightClipRefPos, start_aln_seg_pos, end_aln_seg_pos);
						if(overlap_size>0){
							if(abs(end_aln_seg_pos - start_aln_seg_pos + 1) < (rightClipRefPos - leftClipRefPos + 1))
								val = (double)abs(end_aln_seg_pos - start_aln_seg_pos + 1) / (rightClipRefPos - leftClipRefPos + 1);
							else
								val = (double)(rightClipRefPos - leftClipRefPos + 1) / abs(end_aln_seg_pos - start_aln_seg_pos + 1);
							if(val>=VALID_DUP_INV_SIZE_RATIO_CLIP_REG){
								query_seq = fa_loader.getFastaSeq(minimap2_aln->query_id, minimap2_aln->relative_strand);

								// allocate new node
								clip_reg_tmp = new reg_t();
								clip_reg_tmp->var_type = sv_type;
								clip_reg_tmp->chrname = minimap2_aln->chrname;
								clip_reg_tmp->startRefPos = start_aln_seg_pos + 1;
								clip_reg_tmp->endRefPos = end_aln_seg_pos;
								clip_reg_tmp->startLocalRefPos = minimap2_aln->pafalnsegs.at(0)->startSubjectPos + 1;
								clip_reg_tmp->endLocalRefPos = getEndSubjectPosAlnSeg(paf_alnseg->startSubjectPos, paf_alnseg->opflag, paf_alnseg->seglen);
								clip_reg_tmp->startQueryPos = minimap2_aln->pafalnsegs.at(0)->startQpos + 1;
								clip_reg_tmp->endQueryPos = getEndQueryPosAlnSeg(paf_alnseg->startQpos, paf_alnseg->opflag, paf_alnseg->seglen);
								clip_reg_tmp->query_id = minimap2_aln->query_id;
								clip_reg_tmp->sv_len = end_aln_seg_pos - start_aln_seg_pos + 1;
								clip_reg_tmp->dup_num = 0;
								clip_reg_tmp->minimap2_aln_id = -1;
								clip_reg_tmp->call_success_status = true;
								clip_reg_tmp->short_sv_flag = false;
								clip_reg_tmp->zero_cov_flag = false;
								clip_reg_tmp->aln_seg_end_flag = false;
								clip_reg_tmp->query_pos_invalid_flag = false;
								clip_reg_tmp->large_indel_flag = false;
								clip_reg_tmp->merge_flag = false;
								clip_reg_tmp->rescue_flag = rescue_flag;
								clip_reg_tmp->aln_orient = minimap2_aln->relative_strand;
								clip_reg_tmp->gt_type = -1;
								clip_reg_tmp->gt_seq = "";
								clip_reg_tmp->AF = 0;
								//clip_reg_tmp->supp_num = queryname_vec.size();
								//clip_reg_tmp->supp_num = computeSuppNum(clip_reg_tmp->chrname, clip_reg_tmp->startRefPos, clip_reg_tmp->endRefPos, inBamFile, minMapQ, minHighMapQ, max_ultra_high_cov, queryname_vec);
								//clip_reg_tmp->DP = 0;
								supp_cov_vec = getSuppNumCovClipReg(bnd_mate_reg_strs, supp_num_largeIndel, depth_largeIndel, queryname_vec, large_indel_flag);
								clip_reg_tmp->supp_num = supp_cov_vec.at(0);
								clip_reg_tmp->DP = supp_cov_vec.at(1);
								if(rescue_flag==false) clip_reg_tmp->discover_level = VAR_DISCOV_L_CNS_ALN;
								else clip_reg_tmp->discover_level = VAR_DISCOV_L_RESCUE_CNS_ALN;
								clip_reg_tmp->refseq = ref_seq.substr(clip_reg_tmp->startLocalRefPos-1, clip_reg_tmp->endLocalRefPos-clip_reg_tmp->startLocalRefPos+1);
								clip_reg_tmp->altseq = query_seq.substr(clip_reg_tmp->startQueryPos-1, clip_reg_tmp->endQueryPos-clip_reg_tmp->startQueryPos+1);

								clip_reg_tmp->dup_num = 0;
								clip_reg_tmp->large_indel_flag = false;

								var_vec.push_back(clip_reg_tmp);
								tmp_call_success_flag = true;
								break;
							}
						}
					}
				}
			}

			// analyze the intra-signatures
			if(tmp_call_success_flag==false){
				for(i=0; i<minimap2_item_id_vec.size(); i++){
					minimap2_aln = minimap2_aln_vec.at(minimap2_item_id_vec.at(i));
					if(minimap2_aln->pafalnsegs.size()>0){
						start_aln_seg_pos = minimap2_aln->pafalnsegs.at(0)->startRpos;
						paf_alnseg = minimap2_aln->pafalnsegs.at(minimap2_aln->pafalnsegs.size()-1);
						end_aln_seg_pos = getEndRefPosAlnSeg(paf_alnseg->startRpos, paf_alnseg->opflag, paf_alnseg->seglen);
						if(start_aln_seg_pos+1<=varVec.at(0)->startRefPos and end_aln_seg_pos+1>=varVec.at(varVec.size()-1)->endRefPos and start_aln_seg_pos+var_ext_size<=leftClipRefPos and rightClipRefPos+var_ext_size<=end_aln_seg_pos){
							// search for match align segments for begin_var_idx
							begin_var_idx = end_var_idx = -1;
							for(j=0; j<(int64_t)minimap2_aln->pafalnsegs.size(); j++){
								paf_alnseg = minimap2_aln->pafalnsegs.at(j);
								if(paf_alnseg->opflag==BAM_CMATCH){
									end_ref_pos = getEndRefPosAlnSeg(paf_alnseg->startRpos, paf_alnseg->opflag, paf_alnseg->seglen);
									//if(paf_alnseg->startRpos+1-100<=leftClipRefPos and end_ref_pos+1+100>=leftClipRefPos){
									if(paf_alnseg->startRpos+1-var_ext_size<=leftClipRefPos and end_ref_pos+1+var_ext_size>=leftClipRefPos){ // updated on 2025-02-10
										begin_var_idx = j;
										break;
									}
								}
							}
							if(begin_var_idx!=-1){
								paf_alnseg = minimap2_aln->pafalnsegs.at(begin_var_idx);
								end_ref_pos = getEndRefPosAlnSeg(paf_alnseg->startRpos, paf_alnseg->opflag, paf_alnseg->seglen);
								if(end_ref_pos+1-var_ext_size>rightClipRefPos){
									for(j=begin_var_idx-1; j>=0; j--){
										paf_alnseg = minimap2_aln->pafalnsegs.at(j);
										//end_ref_pos = getEndRefPosAlnSeg(paf_alnseg->startRpos, paf_alnseg->opflag, paf_alnseg->seglen);
										if(paf_alnseg->opflag==BAM_CMATCH and end_ref_pos+1+var_ext_size<rightClipRefPos){
											//if(paf_alnseg->startRpos+1-100<=leftClipRefPos and end_ref_pos+1+100>=leftClipRefPos){ // added on 2025-02-10
											if(paf_alnseg->startRpos+1-var_ext_size<=leftClipRefPos and end_ref_pos+1+var_ext_size>=leftClipRefPos){
												begin_var_idx = j;
												break;
											}else if(end_ref_pos<leftClipRefPos) // added on 2025-02-10
												break;
										}
									}
									if(begin_var_idx<0) begin_var_idx = 0;
								}
							}

							// search for match align segments for end_var_idx
							if(begin_var_idx!=-1){
								for(j=(int64_t)minimap2_aln->pafalnsegs.size()-1; j>begin_var_idx; j--){
									paf_alnseg = minimap2_aln->pafalnsegs.at(j);
									if(paf_alnseg->opflag==BAM_CMATCH){
										end_ref_pos = getEndRefPosAlnSeg(paf_alnseg->startRpos, paf_alnseg->opflag, paf_alnseg->seglen);
										//if(paf_alnseg->startRpos+1-100<=rightClipRefPos and end_ref_pos+1+100>=rightClipRefPos){
										if(paf_alnseg->startRpos+1-var_ext_size<=rightClipRefPos and end_ref_pos+1+var_ext_size>=rightClipRefPos){ // updated on 2025-02-10
											end_var_idx = j;
											break;
										}
									}
								}
								if(end_var_idx!=-1){
									paf_alnseg = minimap2_aln->pafalnsegs.at(end_var_idx);
									if(paf_alnseg->startRpos+1+var_ext_size<leftClipRefPos){
										for(j=end_var_idx+1; j<(int64_t)minimap2_aln->pafalnsegs.size(); j++){
											paf_alnseg = minimap2_aln->pafalnsegs.at(j);
											if(paf_alnseg->opflag==BAM_CMATCH and paf_alnseg->startRpos+1-100>=leftClipRefPos){
												//if(paf_alnseg->startRpos+1-100<=rightClipRefPos and end_ref_pos+1+100>=rightClipRefPos){ // added on 2025-02-10
												if(paf_alnseg->startRpos+1-var_ext_size<=rightClipRefPos and end_ref_pos+1+var_ext_size>=rightClipRefPos){ // updated on 2025-02-10
													end_var_idx = j;
													break;
												}else if(paf_alnseg->startRpos>rightClipRefPos) // added on 2025-02-10
													break;
											}
										}
										if(j>=(int64_t)minimap2_aln->pafalnsegs.size()) end_var_idx = minimap2_aln->pafalnsegs.size() - 1;
										if(end_var_idx>=(int64_t)minimap2_aln->pafalnsegs.size()) end_var_idx = minimap2_aln->pafalnsegs.size() - 1;
									}
								}
							}

							// compute the locations
							if(begin_var_idx!=-1 and end_var_idx!=-1){
								// perfer one segment
								valid_sig_flag = single_seg_flag = false;
								for(j=begin_var_idx; j<=end_var_idx; j++){
									paf_alnseg = minimap2_aln->pafalnsegs.at(j);
									if((paf_alnseg->opflag==BAM_CINS and sv_type==VAR_INS) or (paf_alnseg->opflag==BAM_CDEL and sv_type==VAR_DEL)){
										if(varVec.size()==1 and paf_alnseg->seglen>=min_sv_size){
											if(abs(varVec.at(0)->sv_len)>=paf_alnseg->seglen) val = (double)paf_alnseg->seglen / abs(varVec.at(0)->sv_len);
											else val = (double)abs(varVec.at(0)->sv_len) / paf_alnseg->seglen;
											if(val>=MIN_ALN_SIZE_RATIO_HIGHSIM){
												begin_var_idx = j - 1;
												end_var_idx = j + 1;
												valid_sig_flag = single_seg_flag = true;
												break;
											}
										}
									}
								}

								if(single_seg_flag==false){
									// shrink to fit for indels
									if(sv_type==VAR_INS or sv_type==VAR_DUP){
										// start location adjust
										begin_var_idx_new = begin_var_idx;
										for(j=begin_var_idx+1; j<end_var_idx; j++){
											paf_alnseg = minimap2_aln->pafalnsegs.at(j-1);
											paf_alnseg2 = minimap2_aln->pafalnsegs.at(j);
											if(paf_alnseg->opflag==BAM_CMATCH and (paf_alnseg2->opflag==BAM_CINS) and paf_alnseg2->seglen>5){
												if(paf_alnseg2->seglen>=min_sv_size) valid_sig_flag = true;
												begin_var_idx_new = j - 1;
												break;
											}
										}
										begin_var_idx = begin_var_idx_new;

										if(valid_sig_flag){
											end_var_idx_new = end_var_idx;
											for(j=end_var_idx-1; j>begin_var_idx; j--){
												paf_alnseg = minimap2_aln->pafalnsegs.at(j);
												paf_alnseg2 = minimap2_aln->pafalnsegs.at(j+1);
												if(paf_alnseg2->opflag==BAM_CMATCH and (paf_alnseg->opflag==BAM_CINS) and paf_alnseg->seglen>5){
													if(paf_alnseg2->seglen>=min_sv_size) valid_sig_flag = true;
													end_var_idx_new = j + 1;
													break;
												}
											}
											end_var_idx = end_var_idx_new;
										}
									}else if(sv_type==VAR_DEL){
										// start location adjust
										begin_var_idx_new = begin_var_idx;
										for(j=begin_var_idx+1; j<end_var_idx; j++){
											paf_alnseg = minimap2_aln->pafalnsegs.at(j-1);
											paf_alnseg2 = minimap2_aln->pafalnsegs.at(j);
											if(paf_alnseg->opflag==BAM_CMATCH and (paf_alnseg2->opflag==BAM_CDEL) and paf_alnseg2->seglen>5){
												if(paf_alnseg2->seglen>=min_sv_size) valid_sig_flag = true;
												begin_var_idx_new = j - 1;
												break;
											}
										}
										begin_var_idx = begin_var_idx_new;

										if(valid_sig_flag){
											end_var_idx_new = end_var_idx;
											for(j=end_var_idx-1; j>begin_var_idx; j--){
												paf_alnseg = minimap2_aln->pafalnsegs.at(j);
												paf_alnseg2 = minimap2_aln->pafalnsegs.at(j+1);
												if(paf_alnseg2->opflag==BAM_CMATCH and (paf_alnseg->opflag==BAM_CDEL) and paf_alnseg->seglen>5){
													if(paf_alnseg2->seglen>=min_sv_size) valid_sig_flag = true;
													end_var_idx_new = j + 1;
													break;
												}
											}
											end_var_idx = end_var_idx_new;
										}
									}
								}

								if(valid_sig_flag){
									new_sv_type = -1;
									if(end_var_idx-begin_var_idx==2){ // only one signature among them, added on 2024-08-05
		//								paf_alnseg = minimap2_aln->pafalnsegs.at(begin_var_idx);
		//								start_var_ref_pos = getEndRefPosAlnSeg(paf_alnseg->startRpos, paf_alnseg->opflag, paf_alnseg->seglen);
		//								start_var_subject_pos = getEndSubjectPosAlnSeg(paf_alnseg->startSubjectPos, paf_alnseg->opflag, paf_alnseg->seglen);
		//								start_var_query_pos = getEndQueryPosAlnSeg(paf_alnseg->startQpos, paf_alnseg->opflag, paf_alnseg->seglen);

										paf_alnseg = minimap2_aln->pafalnsegs.at(begin_var_idx+1);
										start_var_ref_pos = paf_alnseg->startRpos + 1;
										start_var_subject_pos = paf_alnseg->startSubjectPos + 1;
										start_var_query_pos = paf_alnseg->startQpos + 1;

										//paf_alnseg = minimap2_aln->pafalnsegs.at(begin_var_idx+1);
										end_var_ref_pos = getEndRefPosAlnSeg(paf_alnseg->startRpos, paf_alnseg->opflag, paf_alnseg->seglen) + 1;
										end_var_subject_pos = getEndSubjectPosAlnSeg(paf_alnseg->startSubjectPos, paf_alnseg->opflag, paf_alnseg->seglen) + 1;
										end_var_query_pos = getEndQueryPosAlnSeg(paf_alnseg->startQpos, paf_alnseg->opflag, paf_alnseg->seglen) + 1;

										if(paf_alnseg->opflag==BAM_CINS) new_sv_type = VAR_INS;
										else if(paf_alnseg->opflag==BAM_CDEL) new_sv_type = VAR_DEL;
									}else{
										// begin_pos
										paf_alnseg = minimap2_aln->pafalnsegs.at(begin_var_idx);
										end_ref_pos = getEndRefPosAlnSeg(paf_alnseg->startRpos, paf_alnseg->opflag, paf_alnseg->seglen);
										start_var_ref_pos = start_var_subject_pos = start_var_query_pos = -1;
										switch(paf_alnseg->opflag){
											case BAM_CMATCH:
												if(end_ref_pos+1>=leftClipRefPos){
													dist_begin = leftClipRefPos - paf_alnseg->startRpos - 1;
													if(dist_begin>=0){ // beg_paf_alnseg is in the clipping region
														start_var_ref_pos = leftClipRefPos;
														start_var_subject_pos = paf_alnseg->startSubjectPos + 1 + dist_begin;
														start_var_query_pos = paf_alnseg->startQpos + 1 + dist_begin;
													}else{ // beg_paf_alnseg is not in the clipping region, update the locations according to alignment
														start_var_ref_pos = getEndRefPosAlnSeg(paf_alnseg->startRpos, paf_alnseg->opflag, paf_alnseg->seglen) + 1;
														start_var_subject_pos = getEndSubjectPosAlnSeg(paf_alnseg->startSubjectPos, paf_alnseg->opflag, paf_alnseg->seglen) + 1;
														start_var_query_pos = getEndQueryPosAlnSeg(paf_alnseg->startQpos, paf_alnseg->opflag, paf_alnseg->seglen) + 1;
													}
												}else{
													start_var_ref_pos = getEndRefPosAlnSeg(paf_alnseg->startRpos, paf_alnseg->opflag, paf_alnseg->seglen) + 1;
													start_var_subject_pos = getEndSubjectPosAlnSeg(paf_alnseg->startSubjectPos, paf_alnseg->opflag, paf_alnseg->seglen) + 1;
													start_var_query_pos = getEndQueryPosAlnSeg(paf_alnseg->startQpos, paf_alnseg->opflag, paf_alnseg->seglen) + 1;
												}
												break;
											case BAM_CINS:
											case BAM_CDEL:
											case BAM_CSOFT_CLIP:
											case BAM_CHARD_CLIP:
												break;
											case BAM_CEQUAL:
											case BAM_CDIFF:
												cout << "line=" << __LINE__ << ", opflag=diff" << ", startRpos=" << start_var_ref_pos << endl;
												break;
											case BAM_CREF_SKIP:
												// unexpected events
											case BAM_CPAD:
											case BAM_CBACK:
											default:
												cerr << __func__ << ", line=" << __LINE__ << ": invalid opflag " << paf_alnseg->opflag << endl;
												exit(1);
										}
										if(start_var_ref_pos<0 or start_var_subject_pos<0 or start_var_query_pos<0){
											cerr << __func__ << ", line=" << __LINE__ << ": invalid start locations: start_var_ref_pos=" << start_var_ref_pos << ", start_var_subject_pos=" << start_var_subject_pos << ", start_var_query_pos=" << start_var_query_pos << endl;
											exit(1);
										}

										// end_pos
										paf_alnseg = minimap2_aln->pafalnsegs.at(end_var_idx);
										end_ref_pos = getEndRefPosAlnSeg(paf_alnseg->startRpos, paf_alnseg->opflag, paf_alnseg->seglen) + 1;
										end_var_ref_pos = end_var_subject_pos = end_var_query_pos = -1;
										switch(paf_alnseg->opflag){
											case BAM_CMATCH:
												if(end_ref_pos>=rightClipRefPos){
													dist_end = rightClipRefPos - paf_alnseg->startRpos - 1;
													if(dist_end>=0){ // end_paf_alnseg is in the clipping region
														end_var_ref_pos = rightClipRefPos;
														end_var_subject_pos = paf_alnseg->startSubjectPos + 1 + dist_end;
														end_var_query_pos = paf_alnseg->startQpos + 1 + dist_end;
													}else{ // end_paf_alnseg is not in the clipping region, update the locations according to alignment
														end_var_ref_pos = paf_alnseg->startRpos + 1;
														end_var_subject_pos = paf_alnseg->startSubjectPos + 1;
														end_var_query_pos = paf_alnseg->startQpos + 1;
													}
												}else{
													end_var_ref_pos = getEndRefPosAlnSeg(paf_alnseg->startRpos, paf_alnseg->opflag, paf_alnseg->seglen) + 1;
													end_var_subject_pos = getEndSubjectPosAlnSeg(paf_alnseg->startSubjectPos, paf_alnseg->opflag, paf_alnseg->seglen) + 1;
													end_var_query_pos = getEndQueryPosAlnSeg(paf_alnseg->startQpos, paf_alnseg->opflag, paf_alnseg->seglen) + 1;
												}
												break;
											case BAM_CINS:
											case BAM_CDEL:
											case BAM_CSOFT_CLIP:
											case BAM_CHARD_CLIP:
												break;
											case BAM_CEQUAL:
											case BAM_CDIFF:
												cout << "line=" << __LINE__ << ", opflag=diff" << ", startRpos=" << start_var_ref_pos << endl;
												break;
											case BAM_CREF_SKIP:
												// unexpected events
											case BAM_CPAD:
											case BAM_CBACK:
											default:
												cerr << __func__ << ", line=" << __LINE__ << ": invalid opflag " << paf_alnseg->opflag << endl;
												exit(1);
										}

										if(end_var_ref_pos<0 or end_var_subject_pos<0 or end_var_query_pos<0){
											cerr << __func__ << ", line=" << __LINE__ << ": invalid start locations: end_var_ref_pos=" << end_var_ref_pos << ", end_var_subject_pos=" << end_var_subject_pos << ", end_var_query_pos=" << end_var_query_pos << endl;
											exit(1);
										}
									}

									// confirm and save
									if((end_var_query_pos - start_var_query_pos) > (end_var_ref_pos - start_var_ref_pos))
										sv_len = (end_var_query_pos - start_var_query_pos) - (end_var_ref_pos - start_var_ref_pos) + 1; // INS
									else{
										sv_len = (end_var_query_pos - start_var_query_pos) - (end_var_ref_pos - start_var_ref_pos) - 1; // DEL
									}
									if(sv_len>=0) new_sv_type = VAR_INS;
									else new_sv_type = VAR_DEL;

									//if(sv_type==VAR_DUP){ // added on 2025-08-20
									if(new_sv_type==sv_type or (abs(sv_len)>=MAX_INNER_MISSING_IGNORE_SIZE and new_sv_type!=sv_type)){ // ignore short signatures for different types
										val = (double)abs(sv_len) / (rightClipRefPos - leftClipRefPos + 1);

										//if(val>=VALID_DUP_INV_SIZE_RATIO_CLIP_REG){ // deleted on 2026-01-29
										if(val>=VALID_DUP_INV_SIZE_RATIO_CLIP_REG or ((new_sv_type==VAR_INS or new_sv_type==VAR_DEL) and cluster_querynames_vec.size()>=2)){ // added on 2026-01-29
											query_seq = fa_loader.getFastaSeq(minimap2_aln->query_id, minimap2_aln->relative_strand);

											// allocate new node
											clip_reg_tmp = new reg_t();
											if(new_sv_type!=-1) clip_reg_tmp->var_type = new_sv_type;
											else clip_reg_tmp->var_type = sv_type;
											clip_reg_tmp->chrname = minimap2_aln->chrname;
											clip_reg_tmp->startRefPos = start_var_ref_pos;
											clip_reg_tmp->endRefPos = end_var_ref_pos;
											clip_reg_tmp->startLocalRefPos = start_var_subject_pos;
											clip_reg_tmp->endLocalRefPos = end_var_subject_pos;
											clip_reg_tmp->startQueryPos = start_var_query_pos;
											clip_reg_tmp->endQueryPos = end_var_query_pos;
											clip_reg_tmp->query_id = minimap2_aln->query_id;
											if(clip_reg_tmp->var_type==VAR_INV) clip_reg_tmp->sv_len = clip_reg_tmp->endRefPos - clip_reg_tmp->startRefPos + 1;
											else clip_reg_tmp->sv_len = sv_len;
											clip_reg_tmp->minimap2_aln_id = minimap2_aln->minimap2_aln_id;
											clip_reg_tmp->call_success_status = true;
											clip_reg_tmp->short_sv_flag = false;
											clip_reg_tmp->zero_cov_flag = false;
											clip_reg_tmp->aln_seg_end_flag = false;
											clip_reg_tmp->query_pos_invalid_flag = false;
											clip_reg_tmp->large_indel_flag = false;
											clip_reg_tmp->merge_flag = false;
											clip_reg_tmp->rescue_flag = rescue_flag;
											clip_reg_tmp->aln_orient = minimap2_aln->relative_strand;
											clip_reg_tmp->gt_type = -1;
											clip_reg_tmp->gt_seq = "";
											clip_reg_tmp->AF = 0;
											//clip_reg_tmp->supp_num = minimap2_aln->qname.size();
											//clip_reg_tmp->supp_num = computeSuppNum(clip_reg_tmp->chrname, clip_reg_tmp->startRefPos, clip_reg_tmp->endRefPos, inBamFile, minMapQ, minHighMapQ, max_ultra_high_cov, minimap2_aln->qname);
											//clip_reg_tmp->DP = 0;
											supp_cov_vec = getSuppNumCovClipReg(bnd_mate_reg_strs, supp_num_largeIndel, depth_largeIndel, queryname_vec, large_indel_flag);
											clip_reg_tmp->supp_num = supp_cov_vec.at(0);
											clip_reg_tmp->DP = supp_cov_vec.at(1);
											if(rescue_flag==false) clip_reg_tmp->discover_level = VAR_DISCOV_L_CNS_ALN;
											else clip_reg_tmp->discover_level = VAR_DISCOV_L_RESCUE_CNS_ALN;
											clip_reg_tmp->refseq =  ref_seq.substr(clip_reg_tmp->startLocalRefPos-1, clip_reg_tmp->endLocalRefPos-clip_reg_tmp->startLocalRefPos+1);
											clip_reg_tmp->altseq = query_seq.substr(clip_reg_tmp->startQueryPos-1, clip_reg_tmp->endQueryPos-clip_reg_tmp->startQueryPos+1);

											if(new_sv_type==-1){
												switch(sv_type){
													case VAR_DUP: // DUP
														clip_reg_tmp->dup_num = dup_num;
														clip_reg_tmp->large_indel_flag = false;
														break;
													case VAR_INV: // INV
														clip_reg_tmp->dup_num = 0;
														clip_reg_tmp->large_indel_flag = false;
														break;
													case VAR_INS: // large INS or DEL
													case VAR_DEL:
														clip_reg_tmp->dup_num = 0;
														clip_reg_tmp->large_indel_flag = true;
														break;
												}
											}

											var_vec.push_back(clip_reg_tmp);
											tmp_call_success_flag = true;
										}
									}//
								}
							}
						}
					}
				}
			}

			// analyze multiple minimap2 items
			if(tmp_call_success_flag==false and minimap2_item_id_vec.size()>1){
				sidemost_item_id_vec = getLeftRightMinimapItemId(leftClipRefPos, rightClipRefPos, minimap2_item_id_vec, minimap2_aln_vec);
				if(sidemost_item_id_vec.size()>0){
					query_id = sidemost_item_id_vec.at(0);
					minimap2_aln_left = minimap2_aln_vec.at(sidemost_item_id_vec.at(1));
					minimap2_aln_right = minimap2_aln_vec.at(sidemost_item_id_vec.at(2));
					if(sv_type==VAR_DUP){ // DUP
						if(minimap2_aln_left->relative_strand==minimap2_aln_right->relative_strand){
							paf_alnseg = minimap2_aln_left->pafalnsegs.at(0);
							paf_alnseg2 = minimap2_aln_right->pafalnsegs.at(minimap2_aln_right->pafalnsegs.size()-1);
							end_qpos_tmp = getEndQueryPosAlnSeg(paf_alnseg2->startQpos, paf_alnseg2->opflag, paf_alnseg2->seglen);
							end_ref_pos = getEndRefPosAlnSeg(paf_alnseg2->startRpos, paf_alnseg2->opflag, paf_alnseg2->seglen);
							dist_query = paf_alnseg->startQpos - end_qpos_tmp;

							entire_flanking_flag = false;
							if(minimap2_aln_left->relative_strand==ALN_PLUS_ORIENT and paf_alnseg->startQpos<end_qpos_tmp){
								entire_flanking_flag = true;
							}else if(minimap2_aln_left->relative_strand==ALN_MINUS_ORIENT and paf_alnseg->startQpos>end_qpos_tmp){
								entire_flanking_flag = true;
							}

							//val = (double)abs(dist_query - dist_ref + 1) / (rightClipRefPos - leftClipRefPos + 1);
							//if(abs(dist_query)<MIN_ALN_SIZE_SAME_ORIENT and val>=MIN_SIZE_RATIO_MATCH_CLIP_POS){ // prefer the original paf information
							if(abs(dist_query)<MIN_ALN_SIZE_SAME_ORIENT and entire_flanking_flag==false){ // prefer the original paf information
								left_rpos = paf_alnseg->startRpos;
								left_subjectpos = paf_alnseg->startSubjectPos;
								left_qpos = paf_alnseg->startQpos;

								right_rpos = end_ref_pos;
								right_subjectpos = getEndSubjectPosAlnSeg(paf_alnseg2->startSubjectPos, paf_alnseg2->opflag, paf_alnseg2->seglen);;
								right_qpos = end_qpos_tmp;

								// save
								sv_len = (right_qpos - left_qpos) - (right_rpos - left_rpos) + 1;
								if(sv_len<0) sv_len = -sv_len; // deleted on 2026-02-26

								// allocate new node
								clip_reg_tmp = new reg_t();
								//if(sv_len<0 and abs(sv_len)>MIN_AVER_SIZE_ALN_SEG) clip_reg_tmp->var_type = VAR_DEL; // deleted on 2026-04-12
								//else
									clip_reg_tmp->var_type = sv_type;
								clip_reg_tmp->chrname = minimap2_aln_left->chrname;
								clip_reg_tmp->startRefPos = left_rpos + 1;
								clip_reg_tmp->endRefPos = right_rpos;
								clip_reg_tmp->startLocalRefPos = left_subjectpos + 1;
								clip_reg_tmp->endLocalRefPos = right_subjectpos;
								clip_reg_tmp->startQueryPos = left_qpos + 1;
								clip_reg_tmp->endQueryPos = right_qpos;
								clip_reg_tmp->query_id = minimap2_aln_left->query_id;
								clip_reg_tmp->sv_len = sv_len;
								clip_reg_tmp->dup_num = dup_num;
								clip_reg_tmp->minimap2_aln_id = -1;
								clip_reg_tmp->call_success_status = true;
								clip_reg_tmp->short_sv_flag = false;
								clip_reg_tmp->zero_cov_flag = false;
								clip_reg_tmp->aln_seg_end_flag = false;
								clip_reg_tmp->query_pos_invalid_flag = false;
								clip_reg_tmp->large_indel_flag = false;
								clip_reg_tmp->merge_flag = false;
								clip_reg_tmp->rescue_flag = rescue_flag;
								clip_reg_tmp->aln_orient = minimap2_aln_left->relative_strand;
								clip_reg_tmp->gt_type = -1;
								clip_reg_tmp->gt_seq = "";
								clip_reg_tmp->AF = 0;
								//clip_reg_tmp->supp_num = queryname_vec.size();
								//clip_reg_tmp->supp_num = computeSuppNum(clip_reg_tmp->chrname, clip_reg_tmp->startRefPos, clip_reg_tmp->endRefPos, inBamFile, minMapQ, minHighMapQ, max_ultra_high_cov, queryname_vec);
								//clip_reg_tmp->DP = 0;
								supp_cov_vec = getSuppNumCovClipReg(bnd_mate_reg_strs, supp_num_largeIndel, depth_largeIndel, queryname_vec, large_indel_flag);
								clip_reg_tmp->supp_num = supp_cov_vec.at(0);
								clip_reg_tmp->DP = supp_cov_vec.at(1);
								if(rescue_flag==false) clip_reg_tmp->discover_level = VAR_DISCOV_L_CNS_ALN;
								else clip_reg_tmp->discover_level = VAR_DISCOV_L_RESCUE_CNS_ALN;

								clip_reg_tmp->refseq = ref_seq.at(clip_reg_tmp->startLocalRefPos-1);
								//query_seq = fa_loader.getFastaSeq(minimap2_aln_left->query_id, minimap2_aln_left->relative_strand);
								//clip_reg_tmp->refseq = ref_seq.substr(clip_reg_tmp->startLocalRefPos-1, clip_reg_tmp->endLocalRefPos-clip_reg_tmp->startLocalRefPos+1);// ????
								//clip_reg_tmp->altseq = query_seq.substr(clip_reg_tmp->startQueryPos-1, clip_reg_tmp->endQueryPos-clip_reg_tmp->startQueryPos+1);

								clip_reg_tmp->dup_num = dup_num;
								clip_reg_tmp->large_indel_flag = false;

								var_vec.push_back(clip_reg_tmp);
								tmp_call_success_flag = true;

							}else{ // compute the clipping information
								// left side
								idx_alnSegs_left = -1;
								for(i=0; i<minimap2_aln_left->pafalnsegs.size(); i++){
									paf_alnseg = minimap2_aln_left->pafalnsegs.at(i);
									if((paf_alnseg->opflag==BAM_CMATCH or paf_alnseg->opflag==BAM_CEQUAL)){
										if(paf_alnseg->startRpos<=leftClipRefPos-1){
											idx_alnSegs_left = i;
										}else{
											if(idx_alnSegs_left==-1 and paf_alnseg->startRpos-leftClipRefPos<=MAX_REF_DIST_SAME_CHR)
												idx_alnSegs_left = i;
											break;
										}
									}
								}

								// right side
								idx_alnSegs_right = -1;
								for(i=0; i<minimap2_aln_right->pafalnsegs.size(); i++){
									paf_alnseg = minimap2_aln_right->pafalnsegs.at(i);
									if((paf_alnseg->opflag==BAM_CMATCH or paf_alnseg->opflag==BAM_CEQUAL)){
										if(paf_alnseg->startRpos<=rightClipRefPos){
											idx_alnSegs_right = i;
										}else{
											if(idx_alnSegs_right==-1 and paf_alnseg->startRpos-rightClipRefPos<=MAX_REF_DIST_SAME_CHR)
												idx_alnSegs_right = i;
											break;
										}
									}
								}

								if(idx_alnSegs_left!=-1 and idx_alnSegs_right!=-1){ // found both margins
									new_sv_type = -1;
									margin_dist1 = leftClipRefPos - 1 - minimap2_aln_left->pafalnsegs.at(idx_alnSegs_left)->startRpos;
									if(margin_dist1<0) margin_dist1 = 0;
									left_rpos = minimap2_aln_left->pafalnsegs.at(idx_alnSegs_left)->startRpos + margin_dist1;
									left_subjectpos = minimap2_aln_left->pafalnsegs.at(idx_alnSegs_left)->startSubjectPos + margin_dist1;
									left_qpos = minimap2_aln_left->pafalnsegs.at(idx_alnSegs_left)->startQpos + margin_dist1;

									margin_dist2 = rightClipRefPos - minimap2_aln_right->pafalnsegs.at(idx_alnSegs_right)->startRpos;
									if(margin_dist2<0) margin_dist2 = 0;
									right_rpos = minimap2_aln_right->pafalnsegs.at(idx_alnSegs_right)->startRpos + margin_dist2;
									right_subjectpos = minimap2_aln_right->pafalnsegs.at(idx_alnSegs_right)->startSubjectPos + margin_dist2;
									right_qpos = minimap2_aln_right->pafalnsegs.at(idx_alnSegs_right)->startQpos + margin_dist2;

									if(left_qpos<right_qpos and left_subjectpos<right_subjectpos){
										// save
										sv_len = (right_qpos - left_qpos) - (right_rpos - left_rpos) + 1;
										//if(sv_len<0) sv_len = -sv_len;
										query_seq = fa_loader.getFastaSeq(minimap2_aln_left->query_id, minimap2_aln_left->relative_strand);

										if(sv_len<0 and abs(sv_len)>MIN_AVER_SIZE_ALN_SEG)
											new_sv_type = VAR_DEL;

										val = (double)abs(sv_len) / (rightClipRefPos - leftClipRefPos + 1); // added on 2025-08-20
										if(val>=VALID_DUP_INV_SIZE_RATIO_CLIP_REG){
											// allocate new node
											clip_reg_tmp = new reg_t();
											if(new_sv_type!=-1) clip_reg_tmp->var_type = new_sv_type;
											else clip_reg_tmp->var_type = sv_type;
											clip_reg_tmp->chrname = minimap2_aln_left->chrname;
											clip_reg_tmp->startRefPos = left_rpos + 1;
											clip_reg_tmp->endRefPos = right_rpos;
											clip_reg_tmp->startLocalRefPos = left_subjectpos + 1;
											clip_reg_tmp->endLocalRefPos = right_subjectpos;
											clip_reg_tmp->startQueryPos = left_qpos + 1;
											clip_reg_tmp->endQueryPos = right_qpos;
											clip_reg_tmp->query_id = minimap2_aln_left->query_id;
											clip_reg_tmp->sv_len = sv_len;
											clip_reg_tmp->dup_num = dup_num;
											clip_reg_tmp->minimap2_aln_id = -1;
											clip_reg_tmp->call_success_status = true;
											clip_reg_tmp->short_sv_flag = false;
											clip_reg_tmp->zero_cov_flag = false;
											clip_reg_tmp->aln_seg_end_flag = false;
											clip_reg_tmp->query_pos_invalid_flag = false;
											clip_reg_tmp->large_indel_flag = false;
											clip_reg_tmp->merge_flag = false;
											clip_reg_tmp->rescue_flag = rescue_flag;
											clip_reg_tmp->aln_orient = minimap2_aln_left->relative_strand;
											clip_reg_tmp->gt_type = -1;
											clip_reg_tmp->gt_seq = "";
											clip_reg_tmp->AF = 0;
											//clip_reg_tmp->supp_num = queryname_vec.size();
											//clip_reg_tmp->supp_num = computeSuppNum(clip_reg_tmp->chrname, clip_reg_tmp->startRefPos, clip_reg_tmp->endRefPos, inBamFile, minMapQ, minHighMapQ, max_ultra_high_cov, queryname_vec);
											//clip_reg_tmp->DP = 0;
											supp_cov_vec = getSuppNumCovClipReg(bnd_mate_reg_strs, supp_num_largeIndel, depth_largeIndel, queryname_vec, large_indel_flag);
											clip_reg_tmp->supp_num = supp_cov_vec.at(0);
											clip_reg_tmp->DP = supp_cov_vec.at(1);
											if(rescue_flag==false) clip_reg_tmp->discover_level = VAR_DISCOV_L_CNS_ALN;
											else clip_reg_tmp->discover_level = VAR_DISCOV_L_RESCUE_CNS_ALN;
											clip_reg_tmp->refseq = ref_seq.substr(clip_reg_tmp->startLocalRefPos-1, clip_reg_tmp->endLocalRefPos-clip_reg_tmp->startLocalRefPos+1);
											clip_reg_tmp->altseq = query_seq.substr(clip_reg_tmp->startQueryPos-1, clip_reg_tmp->endQueryPos-clip_reg_tmp->startQueryPos+1);

											if(new_sv_type==-1){
												switch(sv_type){
													case VAR_DUP: // DUP
														clip_reg_tmp->dup_num = dup_num;
														clip_reg_tmp->large_indel_flag = false;
														break;
													case VAR_INV: // INV
														clip_reg_tmp->dup_num = 0;
														clip_reg_tmp->large_indel_flag = false;
														break;
													case VAR_INS: // large INS or DEL
													case VAR_DEL:
														clip_reg_tmp->dup_num = 0;
														clip_reg_tmp->large_indel_flag = true;
														break;
												}
											}

											var_vec.push_back(clip_reg_tmp);
											tmp_call_success_flag = true;
										}
									}
								}
							}
						}
					}else if(sv_type==VAR_INV){ // INV
						entire_flanking_flag = false;
						minimap2_aln = NULL;

//						if(minimap2_aln_left->subject_start>minimap2_aln_right->subject_start){
//							cout << "********** Exchanged: minimap2_aln_left.subject_start=" << minimap2_aln_left->subject_start << ", minimap2_aln_right.subject_start=" << minimap2_aln_right->subject_start << endl;
//						}

						if(minimap2_aln_left->relative_strand!=minimap2_aln_right->relative_strand){ // different orient
							// get the inverted item
							for(i=0; i<2; i++){
								if(i==0) minimap2_aln = minimap2_aln_left;
								else minimap2_aln = minimap2_aln_right;
								paf_alnseg = minimap2_aln->pafalnsegs.at(minimap2_aln->pafalnsegs.size()-1);
								end_ref_pos = getEndRefPosAlnSeg(paf_alnseg->startRpos, paf_alnseg->opflag, paf_alnseg->seglen);
								overlap_size = getOverlapSize(minimap2_aln->pafalnsegs.at(0)->startRpos, end_ref_pos, leftClipRefPos, rightClipRefPos);
								if(overlap_size<rightClipRefPos-leftClipRefPos+1)
									val = (double)overlap_size / (rightClipRefPos - leftClipRefPos + 1);
								else
									val = (double)(rightClipRefPos - leftClipRefPos + 1) / overlap_size;
								// added on 2026-04-14
								ref_dist = end_ref_pos - minimap2_aln->pafalnsegs.at(0)->startRpos + 1;
								if(ref_dist<rightClipRefPos-leftClipRefPos+1)
									val2 = (double)ref_dist / (rightClipRefPos - leftClipRefPos + 1);
								else
									val2 = (double)(rightClipRefPos - leftClipRefPos + 1) / ref_dist;
								//if(val>=max_seg_size_ratio){ // deleted on 2026-04-14
								if(val>=max_seg_size_ratio and val2>=VALID_DUP_INV_SIZE_RATIO_CLIP_REG){
									entire_flanking_flag = true;
									//cout << "INV sim_val=" << val << endl;
									break;
								}
							}
						}else{ // same orient, then try to find the entire flanking item
							start_var_query_pos = minimap2_aln_left->query_end + 1;
							end_var_query_pos = minimap2_aln_right->query_start - 1;
							start_var_subject_pos = minimap2_aln_left->subject_end + 1;
							end_var_subject_pos = minimap2_aln_right->subject_start - 1;
							dist_query = end_var_query_pos - start_var_query_pos + 1;
							dist_subject = end_var_subject_pos - start_var_subject_pos + 1;
							if(dist_query<rightClipRefPos-leftClipRefPos+1)
								val = (double)dist_query / (rightClipRefPos - leftClipRefPos + 1);
							else
								val = (double)(rightClipRefPos - leftClipRefPos + 1) / dist_query;
							if(val>=max_seg_size_ratio){
								maxVal = 0;
								maxIdx = -1;
								for(i=0; i<minimap2_aln_vec.size(); i++){
									minimap2_aln = minimap2_aln_vec.at(i);
									if(minimap2_aln->query_id==query_id and minimap2_aln->relative_strand!=minimap2_aln_left->relative_strand){
										overlap_size = getOverlapSize(minimap2_aln->query_start, minimap2_aln->query_end, start_var_query_pos, end_var_query_pos);
										if(overlap_size<dist_query)
											val = (double)overlap_size / dist_query;
										else
											val = (double)dist_query / overlap_size;
										dist_subject_tmp = minimap2_aln->subject_end - minimap2_aln->subject_start + 1;
										if(dist_subject_tmp<dist_subject)
											val2 = (double)dist_subject_tmp / dist_subject;
										else
											val2 = (double)dist_subject / dist_subject_tmp;
										if(val>=max_seg_size_ratio and val2>=VALID_DUP_INV_SIZE_RATIO_CLIP_REG){
											if(maxVal<val2) {
												maxVal = val2;
												maxIdx = i;
											}
										}
									}
								}
								if(maxVal>=max_seg_size_ratio){
									//cout << "INV sim_val=" << maxVal << endl;
									entire_flanking_flag = true;
									minimap2_aln = minimap2_aln_vec.at(maxIdx);
								}
							}
						}

						if(entire_flanking_flag and minimap2_aln!=NULL){ // found
							if(minimap2_aln->relative_strand==ALN_PLUS_ORIENT) aln_orient = ALN_MINUS_ORIENT;
							else aln_orient = ALN_PLUS_ORIENT;
							query_seq = fa_loader.getFastaSeq(minimap2_aln->query_id, aln_orient);

							paf_alnseg = minimap2_aln->pafalnsegs.at(0);
							paf_alnseg2 = minimap2_aln->pafalnsegs.at(minimap2_aln->pafalnsegs.size()-1);
							end_ref_pos = getEndRefPosAlnSeg(paf_alnseg2->startRpos, paf_alnseg2->opflag, paf_alnseg2->seglen);

							// allocate new node
							clip_reg_tmp = new reg_t();
							clip_reg_tmp->var_type = sv_type;
							clip_reg_tmp->chrname = minimap2_aln->chrname;
							clip_reg_tmp->startRefPos = paf_alnseg->startRpos + 1;
							clip_reg_tmp->endRefPos = end_ref_pos;
							clip_reg_tmp->startLocalRefPos = minimap2_aln->subject_start + 1;
							clip_reg_tmp->endLocalRefPos = minimap2_aln->subject_end;
							clip_reg_tmp->startQueryPos = minimap2_aln->query_start + 1;
							clip_reg_tmp->endQueryPos = minimap2_aln->query_end;
							clip_reg_tmp->query_id = minimap2_aln->query_id;
							clip_reg_tmp->sv_len = minimap2_aln->subject_end - minimap2_aln->subject_start + 1;
							//if(sv_len<0) clip_reg_tmp->sv_len = -sv_len;
							//clip_reg_tmp->dup_num = dup_num; //--------------
							clip_reg_tmp->minimap2_aln_id = minimap2_aln->minimap2_aln_id;
							clip_reg_tmp->call_success_status = true;
							clip_reg_tmp->short_sv_flag = false;
							clip_reg_tmp->zero_cov_flag = false;
							clip_reg_tmp->aln_seg_end_flag = false;
							clip_reg_tmp->query_pos_invalid_flag = false;
							clip_reg_tmp->large_indel_flag = false;
							clip_reg_tmp->merge_flag = false;
							clip_reg_tmp->rescue_flag = rescue_flag;
							clip_reg_tmp->aln_orient = minimap2_aln->relative_strand;
							clip_reg_tmp->gt_type = -1;
							clip_reg_tmp->gt_seq = "";
							clip_reg_tmp->AF = 0;
							//clip_reg_tmp->supp_num = queryname_vec.size();
							//clip_reg_tmp->supp_num = computeSuppNum(clip_reg_tmp->chrname, clip_reg_tmp->startRefPos, clip_reg_tmp->endRefPos, inBamFile, minMapQ, minHighMapQ, max_ultra_high_cov, queryname_vec);
							//clip_reg_tmp->DP = 0;
							supp_cov_vec = getSuppNumCovClipReg(bnd_mate_reg_strs, supp_num_largeIndel, depth_largeIndel, queryname_vec, large_indel_flag);
							clip_reg_tmp->supp_num = supp_cov_vec.at(0);
							clip_reg_tmp->DP = supp_cov_vec.at(1);
							if(rescue_flag==false) clip_reg_tmp->discover_level = VAR_DISCOV_L_CNS_ALN;
							else clip_reg_tmp->discover_level = VAR_DISCOV_L_RESCUE_CNS_ALN;
							clip_reg_tmp->refseq = ref_seq.substr(clip_reg_tmp->startLocalRefPos-1, clip_reg_tmp->endLocalRefPos-clip_reg_tmp->startLocalRefPos);
							clip_reg_tmp->altseq = query_seq.substr(clip_reg_tmp->startQueryPos-1, clip_reg_tmp->endQueryPos-clip_reg_tmp->startQueryPos);

							switch(sv_type){
								case VAR_DUP: // DUP
									clip_reg_tmp->dup_num = dup_num;
									clip_reg_tmp->large_indel_flag = false;
									break;
								case VAR_INV: // INV
									clip_reg_tmp->dup_num = 0;
									clip_reg_tmp->large_indel_flag = false;
									break;
								case VAR_INS: // large INS or DEL
								case VAR_DEL:
									clip_reg_tmp->dup_num = 0;
									clip_reg_tmp->large_indel_flag = true;
									break;
							}

							var_vec.push_back(clip_reg_tmp);
							tmp_call_success_flag = true;
						}else{
							// left side
							idx_alnSegs_left = -1;
							for(i=0; i<minimap2_aln_left->pafalnsegs.size(); i++){
								paf_alnseg = minimap2_aln_left->pafalnsegs.at(i);
								if(paf_alnseg->opflag==BAM_CMATCH or paf_alnseg->opflag==BAM_CEQUAL){
									if(paf_alnseg->startRpos<=leftClipRefPos-1){
										idx_alnSegs_left = i;
									}else{
										if(idx_alnSegs_left==-1 and paf_alnseg->startRpos-leftClipRefPos<=MAX_REF_DIST_SAME_CHR)
											idx_alnSegs_left = i;
										break;
									}
								}
							}

							// right side
							idx_alnSegs_right = -1;
							for(i=0; i<minimap2_aln_right->pafalnsegs.size(); i++){
								paf_alnseg = minimap2_aln_right->pafalnsegs.at(i);
								if(paf_alnseg->opflag==BAM_CMATCH or paf_alnseg->opflag==BAM_CEQUAL){
									if(paf_alnseg->startRpos<=rightClipRefPos){
										idx_alnSegs_right = i;
									}else{
										if(idx_alnSegs_right==-1 and paf_alnseg->startRpos-rightClipRefPos<=MAX_REF_DIST_SAME_CHR)
											idx_alnSegs_right = i;
										break;
									}
								}
							}

							if(idx_alnSegs_left!=-1 and idx_alnSegs_right!=-1){ // found both margins
								margin_dist1 = leftClipRefPos - 1 - minimap2_aln_left->pafalnsegs.at(idx_alnSegs_left)->startRpos;
								if(margin_dist1<0) margin_dist1 = 0;
								left_rpos = minimap2_aln_left->pafalnsegs.at(idx_alnSegs_left)->startRpos + margin_dist1 + 1;
								left_subjectpos = minimap2_aln_left->pafalnsegs.at(idx_alnSegs_left)->startSubjectPos + margin_dist1 + 1;
								left_qpos = minimap2_aln_left->pafalnsegs.at(idx_alnSegs_left)->startQpos + margin_dist1 + 1;

								margin_dist2 = rightClipRefPos - minimap2_aln_right->pafalnsegs.at(idx_alnSegs_right)->startRpos;
								if(margin_dist2<0) margin_dist2 = 0;
								right_rpos = minimap2_aln_right->pafalnsegs.at(idx_alnSegs_right)->startRpos + margin_dist2 + 1;
								right_subjectpos = minimap2_aln_right->pafalnsegs.at(idx_alnSegs_right)->startSubjectPos + margin_dist2 + 1;
								right_qpos = minimap2_aln_right->pafalnsegs.at(idx_alnSegs_right)->startQpos + margin_dist2 + 1;

								if(left_qpos<right_qpos and left_subjectpos<right_subjectpos){
									// save
									//sv_len = (right_qpos - left_qpos) - (right_rpos - left_rpos) + 1;
									query_seq = fa_loader.getFastaSeq(minimap2_aln_left->query_id, minimap2_aln_left->relative_strand);

									tmp_refseq = ref_seq.substr(left_subjectpos-1, right_subjectpos-left_subjectpos+1);
									tmp_queryseq = query_seq.substr(left_qpos-1, right_qpos-left_qpos+1);

									if(sidemost_item_id_vec.at(1)==sidemost_item_id_vec.at(2)) // same segment
										val = computeVarseqSim(tmp_queryseq, tmp_refseq);
									else{ // different segments
										reverseComplement(tmp_queryseq);
										val = computeVarseqSim(tmp_queryseq, tmp_refseq);
										//cout << "seqsim=" << val << endl;
									}
									//if(val>=QC_SEQSIM_RATIO_MATCH_THRES){ // similarity confirm
									if(val>=min_seqsim_match){ // similarity confirm
										// allocate new node
										clip_reg_tmp = new reg_t();
										clip_reg_tmp->var_type = sv_type;
										clip_reg_tmp->chrname = minimap2_aln_left->chrname;
										clip_reg_tmp->startRefPos = left_rpos;
										clip_reg_tmp->endRefPos = right_rpos;
										clip_reg_tmp->startLocalRefPos = left_subjectpos;
										clip_reg_tmp->endLocalRefPos = right_subjectpos;
										clip_reg_tmp->startQueryPos = left_qpos;
										clip_reg_tmp->endQueryPos = right_qpos;
										clip_reg_tmp->query_id = minimap2_aln_left->query_id;
										clip_reg_tmp->sv_len = clip_reg_tmp->endRefPos - clip_reg_tmp->startRefPos + 1;
										//if(sv_len<0) clip_reg_tmp->sv_len = -sv_len;
										//clip_reg_tmp->dup_num = dup_num; //--------------
										clip_reg_tmp->minimap2_aln_id = -1;
										clip_reg_tmp->call_success_status = true;
										clip_reg_tmp->short_sv_flag = false;
										clip_reg_tmp->zero_cov_flag = false;
										clip_reg_tmp->aln_seg_end_flag = false;
										clip_reg_tmp->query_pos_invalid_flag = false;
										clip_reg_tmp->large_indel_flag = false;
										clip_reg_tmp->merge_flag = false;
										clip_reg_tmp->rescue_flag = rescue_flag;
										clip_reg_tmp->aln_orient = minimap2_aln_left->relative_strand;
										clip_reg_tmp->gt_type = -1;
										clip_reg_tmp->gt_seq = "";
										clip_reg_tmp->AF = 0;
										//clip_reg_tmp->supp_num = queryname_vec.size();
										//clip_reg_tmp->supp_num = computeSuppNum(clip_reg_tmp->chrname, clip_reg_tmp->startRefPos, clip_reg_tmp->endRefPos, inBamFile, minMapQ, minHighMapQ, max_ultra_high_cov, queryname_vec);;
										//clip_reg_tmp->DP = 0;
										supp_cov_vec = getSuppNumCovClipReg(bnd_mate_reg_strs, supp_num_largeIndel, depth_largeIndel, queryname_vec, large_indel_flag);
										clip_reg_tmp->supp_num = supp_cov_vec.at(0);
										clip_reg_tmp->DP = supp_cov_vec.at(1);
										if(rescue_flag==false) clip_reg_tmp->discover_level = VAR_DISCOV_L_CNS_ALN;
										else clip_reg_tmp->discover_level = VAR_DISCOV_L_RESCUE_CNS_ALN;
//										clip_reg_tmp->refseq = ref_seq.substr(clip_reg_tmp->startLocalRefPos-1, clip_reg_tmp->endLocalRefPos-clip_reg_tmp->startLocalRefPos);
//										clip_reg_tmp->altseq = query_seq.substr(clip_reg_tmp->startQueryPos-1, clip_reg_tmp->endQueryPos-clip_reg_tmp->startQueryPos);
										clip_reg_tmp->refseq = tmp_refseq;
										clip_reg_tmp->altseq = tmp_queryseq;

										// deal with identical ALT and REF alleles, change to symbolic type, added on 2025-11-16
										if(val==1){
											upperSeq(tmp_refseq);
											upperSeq(tmp_queryseq);
											if(tmp_refseq.compare(tmp_queryseq)==0){
												clip_reg_tmp->refseq = clip_reg_tmp->refseq.at(0);
												clip_reg_tmp->altseq = "";
											}
										}

										switch(sv_type){
											case VAR_DUP: // DUP
												clip_reg_tmp->dup_num = dup_num;
												clip_reg_tmp->large_indel_flag = false;
												break;
											case VAR_INV: // INV
												clip_reg_tmp->dup_num = 0;
												clip_reg_tmp->large_indel_flag = false;
												break;
											case VAR_INS: // large INS or DEL
											case VAR_DEL:
												clip_reg_tmp->dup_num = 0;
												clip_reg_tmp->large_indel_flag = true;
												break;
										}

										var_vec.push_back(clip_reg_tmp);
										tmp_call_success_flag = true;
									}
								}
							}
						}
					}else if(sv_type==VAR_DEL or sv_type==VAR_INS){ // DEL, INS
//						if(sv_type==VAR_INS) cout << "line=" << __LINE__ << ", sv_type=" << to_string(sv_type) << ", " << alnfilename << ", large INS." << endl;
//						else cout << "line=" << __LINE__ << ", sv_type=" << to_string(sv_type) << ", " << alnfilename << ", large DEL." << endl;

						// left side
						idx_alnSegs_left = -1;
						for(i=0; i<minimap2_aln_left->pafalnsegs.size(); i++){
							paf_alnseg = minimap2_aln_left->pafalnsegs.at(i);
							if((paf_alnseg->opflag==BAM_CMATCH or paf_alnseg->opflag==BAM_CEQUAL)){
								if(paf_alnseg->startRpos<=leftClipRefPos-1){
									idx_alnSegs_left = i;
								}else{
									if(idx_alnSegs_left==-1 and paf_alnseg->startRpos-leftClipRefPos<=MAX_REF_DIST_SAME_CHR)
										idx_alnSegs_left = i;
									break;
								}
							}
						}

						// right side
						idx_alnSegs_right = -1;
						for(i=0; i<minimap2_aln_right->pafalnsegs.size(); i++){
							paf_alnseg = minimap2_aln_right->pafalnsegs.at(i);
							if((paf_alnseg->opflag==BAM_CMATCH or paf_alnseg->opflag==BAM_CEQUAL)){
								if(paf_alnseg->startRpos<=rightClipRefPos){
									idx_alnSegs_right = i;
								}else{
									if(idx_alnSegs_right==-1 and paf_alnseg->startRpos-rightClipRefPos<=MAX_REF_DIST_SAME_CHR)
										idx_alnSegs_right = i;
									break;
								}
							}
						}

						if(idx_alnSegs_left!=-1 and idx_alnSegs_right!=-1){ // found both margins
							//cout << "idx_alnSegs_left=" << idx_alnSegs_left << ", idx_alnSegs_right=" << idx_alnSegs_right << endl;

							margin_dist1 = leftClipRefPos - 1 - minimap2_aln_left->pafalnsegs.at(idx_alnSegs_left)->startRpos;
							if(margin_dist1<0) margin_dist1 = 0;
							left_rpos = minimap2_aln_left->pafalnsegs.at(idx_alnSegs_left)->startRpos + margin_dist1 + 1;
							left_subjectpos = minimap2_aln_left->pafalnsegs.at(idx_alnSegs_left)->startSubjectPos + margin_dist1 + 1;
							left_qpos = minimap2_aln_left->pafalnsegs.at(idx_alnSegs_left)->startQpos + margin_dist1 + 1;

							margin_dist2 = rightClipRefPos - minimap2_aln_right->pafalnsegs.at(idx_alnSegs_right)->startRpos;
							if(margin_dist2<0) margin_dist2 = 0;
							right_rpos = minimap2_aln_right->pafalnsegs.at(idx_alnSegs_right)->startRpos + margin_dist2 + 1;
							right_subjectpos = minimap2_aln_right->pafalnsegs.at(idx_alnSegs_right)->startSubjectPos + margin_dist2 + 1;
							right_qpos = minimap2_aln_right->pafalnsegs.at(idx_alnSegs_right)->startQpos + margin_dist2 + 1;

							if(left_qpos<right_qpos and left_subjectpos<right_subjectpos){
								sv_len = (right_qpos - left_qpos) - (right_rpos - left_rpos) + 1;
								if((sv_len>=0 and sv_type==VAR_INS) or (sv_len<=0 and sv_type==VAR_DEL)){ // consistent type
									if(large_indel_flag){
										if(abs(sv_len)<=abs(varVec.at(0)->sv_len)) val = (double)abs(sv_len) / abs(varVec.at(0)->sv_len);
										else val = (double)abs(varVec.at(0)->sv_len) / abs(sv_len);
									}else{
										if(abs(sv_len)<=(rightClipRefPos - leftClipRefPos + 1)) val = (double)abs(sv_len) / (rightClipRefPos - leftClipRefPos + 1);
										else val = (double)(rightClipRefPos - leftClipRefPos + 1) / abs(sv_len);
									}

									// save
									if(val>=VALID_DUP_INV_SIZE_RATIO_CLIP_REG){
										query_seq = fa_loader.getFastaSeq(minimap2_aln_left->query_id, minimap2_aln_left->relative_strand);

										// allocate new node
										clip_reg_tmp = new reg_t();
										if(sv_len>0) clip_reg_tmp->var_type = VAR_INS;
										else clip_reg_tmp->var_type = VAR_DEL;
										//clip_reg_tmp->var_type = sv_type;
										clip_reg_tmp->chrname = minimap2_aln_left->chrname;
										clip_reg_tmp->startRefPos = left_rpos;
										clip_reg_tmp->endRefPos = right_rpos;
										clip_reg_tmp->startLocalRefPos = left_subjectpos;
										clip_reg_tmp->endLocalRefPos = right_subjectpos;
										clip_reg_tmp->startQueryPos = left_qpos;
										clip_reg_tmp->endQueryPos = right_qpos;
										clip_reg_tmp->query_id = minimap2_aln_left->query_id;
										clip_reg_tmp->sv_len = sv_len;
										//if(sv_len<0) clip_reg_tmp->sv_len = -sv_len;
										//clip_reg_tmp->dup_num = dup_num; //--------------
										clip_reg_tmp->minimap2_aln_id = -1;
										clip_reg_tmp->call_success_status = true;
										clip_reg_tmp->short_sv_flag = false;
										clip_reg_tmp->zero_cov_flag = false;
										clip_reg_tmp->aln_seg_end_flag = false;
										clip_reg_tmp->query_pos_invalid_flag = false;
										clip_reg_tmp->large_indel_flag = false;
										clip_reg_tmp->merge_flag = false;
										clip_reg_tmp->rescue_flag = rescue_flag;
										clip_reg_tmp->aln_orient = minimap2_aln_left->relative_strand;
										clip_reg_tmp->gt_type = -1;
										clip_reg_tmp->gt_seq = "";
										clip_reg_tmp->AF = 0;
										//clip_reg_tmp->supp_num = queryname_vec.size();
										//clip_reg_tmp->supp_num = computeSuppNum(clip_reg_tmp->chrname, clip_reg_tmp->startRefPos, clip_reg_tmp->endRefPos, inBamFile, minMapQ, minHighMapQ, max_ultra_high_cov, queryname_vec);;
										//clip_reg_tmp->DP = 0;
										supp_cov_vec = getSuppNumCovClipReg(bnd_mate_reg_strs, supp_num_largeIndel, depth_largeIndel, queryname_vec, large_indel_flag);
										clip_reg_tmp->supp_num = supp_cov_vec.at(0);
										clip_reg_tmp->DP = supp_cov_vec.at(1);
										if(rescue_flag==false) clip_reg_tmp->discover_level = VAR_DISCOV_L_CNS_ALN;
										else clip_reg_tmp->discover_level = VAR_DISCOV_L_RESCUE_CNS_ALN;
										clip_reg_tmp->refseq = ref_seq.substr(clip_reg_tmp->startLocalRefPos-1, clip_reg_tmp->endLocalRefPos-clip_reg_tmp->startLocalRefPos);
										clip_reg_tmp->altseq = query_seq.substr(clip_reg_tmp->startQueryPos-1, clip_reg_tmp->endQueryPos-clip_reg_tmp->startQueryPos);

										switch(sv_type){
											case VAR_DUP: // DUP
												clip_reg_tmp->dup_num = dup_num;
												clip_reg_tmp->large_indel_flag = false;
												break;
											case VAR_INV: // INV
												clip_reg_tmp->dup_num = 0;
												clip_reg_tmp->large_indel_flag = false;
												break;
											case VAR_INS: // large INS or DEL
											case VAR_DEL:
												clip_reg_tmp->dup_num = 0;
												clip_reg_tmp->large_indel_flag = true;
												break;
										}

										var_vec.push_back(clip_reg_tmp);
										tmp_call_success_flag = true;
									}
								}
							}
						}
					}
				}
			}
		}

		if(rescue_flag){
			//delete minimap2_aln_vec
			destroyMinimap2AlnVec(minimap2_aln_vec);
		}
	}

	return var_vec;
}

vector<int32_t> varCand::getClusterIdCalledIncomplete(vector<reg_t*> &var_vec, string &clusterfilename, vector<minimap2_aln_t*> &minimap2_aln_vec){
	bool complete_flag;
	size_t i, cluster_id;
	vector<string> queryname_vec;
	minimap2_aln_t *mm2_aln_item;
	vector<vector<string>> cluster_querynames_vec;
	vector<int32_t> minimap2_item_id_vec, clusterId_incomplete;
	int32_t query_id;
	reg_t *reg;

	cluster_querynames_vec = getClusterInfo(clusterfilename);

	for(cluster_id=0; cluster_id<cluster_querynames_vec.size(); cluster_id++){
		queryname_vec = cluster_querynames_vec.at(cluster_id);
		minimap2_item_id_vec = getMinimapItemIdVec(queryname_vec, minimap2_aln_vec);
		complete_flag = false;
		if(minimap2_item_id_vec.size()>0){
			mm2_aln_item = minimap2_aln_vec.at(minimap2_item_id_vec.at(0));
			query_id = mm2_aln_item->query_id;
			for(i=0; i<var_vec.size(); i++){
				reg = var_vec.at(i);
				if(reg->query_id==query_id){
					complete_flag = true;
					break;
				}
			}
		}
		if(complete_flag==false){
			clusterId_incomplete.push_back(cluster_id);
		}
	}

	return clusterId_incomplete;
}

// get the support number and coverage for clipping regions
vector<int32_t> varCand::getSuppNumCovClipReg(string bnd_mate_reg_strs[4], int32_t supp_num_largeIndel, int32_t depth_largeIndel, vector<string> &target_qname_vec, bool large_indel_flag){
	vector<int32_t> supp_cov_vec;
	int32_t i, sum, supp_val, cov_val, min_supp_val, max_cov_val;
	size_t j;
	vector<string> str_vec, str_vec2, str_vec2_tmp;

	supp_val = cov_val = 0;
	if(large_indel_flag==false){
		min_supp_val = INT_MAX;
		max_cov_val = INT_MIN;
		for(i=0; i<4; i++){
			if(bnd_mate_reg_strs[i].compare("-")!=0){
				str_vec = split(bnd_mate_reg_strs[i], ",");
				for(j=0; j<str_vec.size(); j++){
					if(str_vec.at(j).compare("-")!=0){
						if(str_vec.at(j).find('_')==string::npos){ // single mate region
							str_vec2 = split(str_vec.at(j), "|");
							supp_val = stoi(str_vec2.at(2));
							cov_val = stoi(str_vec2.at(3));
							if(supp_val<min_supp_val) min_supp_val = supp_val;
							if(cov_val>max_cov_val) max_cov_val = cov_val;
						}else{ // two mate regions
							sum = 0;
							str_vec2_tmp = split(str_vec.at(j), "_");
							for(const auto& item : str_vec2_tmp){
								str_vec2 = split(item, "|");
								supp_val = stoi(str_vec2.at(2));
								cov_val = stoi(str_vec2.at(3));
								sum += supp_val;
								//if(supp_val<min_supp_val) min_supp_val = supp_val;
								if(cov_val>max_cov_val) max_cov_val = cov_val;
							}
							if(sum<min_supp_val) min_supp_val = sum;
						}
					}
				}
			}
		}

		supp_val = max(min_supp_val, (int32_t)target_qname_vec.size());
		cov_val = max(max_cov_val, min_supp_val);
	}else{ // large indel
		supp_val = supp_num_largeIndel;
		cov_val = depth_largeIndel;
	}

	supp_cov_vec.push_back(supp_val);
	supp_cov_vec.push_back(cov_val);

	return supp_cov_vec;
}

// compute genotype
void varCand::computeGenotypeClipReg(vector<reg_t*> &var_vec){

	struct alleleNode{
		reg_t *reg;
		int32_t supp_num, depth : 27, gt_type: 5;
		struct alleleNode *mate_allele;
	};

	size_t  i, j;
	reg_t *reg_tmp;
	int32_t depth;
	double val;
	string gt_header, gt_str, dp_str, ad_str1, ad_str2;

	vector<struct alleleNode*> allele_vec;
	struct alleleNode *allele_node, *allele_node2;

	if(var_vec.size()>0){
		for(i=0; i<var_vec.size(); i++){
			allele_node = new struct alleleNode();
			allele_node->reg = var_vec.at(i);
			allele_node->mate_allele = NULL;
			allele_node->supp_num = var_vec.at(i)->supp_num;
			allele_node->depth = var_vec.at(i)->DP;
			allele_vec.push_back(allele_node);
		}
	}

	for(i=0; i<allele_vec.size(); i++){
		if(large_indel_flag==false){ // clipping, removed on 2025-07-31
			reg_tmp = allele_vec.at(i)->reg;
			allele_vec.at(i)->depth = computeCovNumReg(reg_tmp->chrname, reg_tmp->startRefPos, reg_tmp->endRefPos, fai, inBamFile, minMapQ, minHighMapQ, max_ultra_high_cov);
			if(allele_vec.at(i)->depth<allele_vec.at(i)->supp_num) allele_vec.at(i)->depth = allele_vec.at(i)->supp_num; // tolerate DEL region
		}else // large indel
			allele_vec.at(i)->depth = depth_largeIndel;

		if(allele_vec.at(i)->depth<allele_vec.at(i)->supp_num) // refine the supp_num
			allele_vec.at(i)->supp_num = allele_vec.at(i)->reg->supp_num = allele_vec.at(i)->depth;
	}

	// sort
	for(i=0; i<allele_vec.size(); i++){
		for(j=i+1; j<allele_vec.size(); j++){
			if(allele_vec.at(i)->supp_num<allele_vec.at(j)->supp_num){
				allele_node = allele_vec.at(i);
				allele_vec.at(i) = allele_vec.at(j);
				allele_vec.at(j) = allele_node;
			}
		}
	}

	if(allele_vec.size()>0){
		for(i=0; i<allele_vec.size(); i++){
			allele_node = allele_vec.at(i);
			for(j=i+1; j<allele_vec.size(); ){
				allele_node2 = allele_vec.at(j);
				if(allele_node->reg->var_type==VAR_DEL) val = computeVarseqSim(allele_node->reg->refseq, allele_node2->reg->refseq);
				else val = computeVarseqSim(allele_node->reg->altseq, allele_node2->reg->altseq);
				//cout << "allele seqsim=" << val << endl;
				if(val>=gt_min_seqsim_merge){ // merge
					allele_node->supp_num += allele_node2->supp_num;
					allele_node->depth = allele_node->supp_num;
					delete allele_node2->reg;
					delete allele_node2;
					allele_vec.erase(allele_vec.begin()+j);
				}else{
					j++;
				}
			}
		}
		if(allele_vec.size()>1){
			depth = 0;
			for(i=0; i<allele_vec.size(); i++) depth += allele_vec.at(i)->supp_num;
			for(i=0; i<allele_vec.size(); i++) allele_vec.at(i)->depth = depth;
		}

		if(allele_vec.size()==1){
			// compute the genotype
			gt_header = "GT:AD:DP";
			gt_str = "./.";
			ad_str1 = ".";
			ad_str2 = ".";
			dp_str = ".";

			allele_node = allele_vec.at(0);
			val = (double)allele_node->supp_num / allele_node->depth;
			if(val>=gt_homo_ratio_thres){
				allele_node->reg->gt_type = GT_HOMOZYGOUS;
				allele_node->gt_type = GT_HOMOZYGOUS;
				gt_str = GT_HOMOZYGOUS_STR;
			}else if(val>=gt_hete_ratio_thres){
				allele_node->reg->gt_type = GT_HETEROZYGOUS;
				allele_node->gt_type = GT_HETEROZYGOUS;
				gt_str = GT_HETEROZYGOUS_STR;
			}else{
				allele_node->reg->gt_type = GT_NOZYGOUS;
				allele_node->gt_type = GT_NOZYGOUS;
				gt_str = GT_NOZYGOUS_STR;
			}
			allele_node->reg->DP = allele_node->depth;
			allele_node->reg->AF = val;

			ad_str1 = to_string(allele_node->depth-allele_node->supp_num);
			ad_str2 = to_string(allele_node->supp_num);
			dp_str = to_string(allele_node->depth);

			allele_node->reg->gt_seq = gt_header + "\t" + gt_str + ":" + ad_str1 + "," + ad_str2 + ":" + dp_str;

		}else if(allele_vec.size()>1){
			for(i=0; i<allele_vec.size(); i++){
				allele_node = allele_vec.at(i);
				reg_tmp = allele_node->reg;
				reg_tmp->DP = allele_node->depth;
				reg_tmp->AF = (double)allele_node->supp_num / allele_node->depth;

				if(reg_tmp->call_success_status){
					// compute the genotype
					gt_header = "GT:AD:DP";
					gt_str = "./.";
					ad_str1 = ".";
					ad_str2 = ".";
					dp_str = ".";

					if(reg_tmp->AF>=gt_homo_ratio_thres){
						reg_tmp->gt_type = GT_HOMOZYGOUS;
						allele_node->gt_type = GT_HOMOZYGOUS;
						gt_str = GT_HOMOZYGOUS_STR;
					}else if(reg_tmp->AF>=gt_hete_ratio_thres){
						reg_tmp->gt_type = GT_HETEROZYGOUS;
						allele_node->gt_type = GT_HETEROZYGOUS;
						gt_str = GT_HETEROZYGOUS_STR;
					}else{
						reg_tmp->gt_type = GT_NOZYGOUS;
						allele_node->gt_type = GT_NOZYGOUS;
						gt_str = GT_NOZYGOUS_STR;
					}

					ad_str1 = to_string(allele_node->depth-allele_node->supp_num);
					ad_str2 = to_string(allele_node->supp_num);
					dp_str = to_string(allele_node->depth);

					reg_tmp->gt_seq = gt_header + "\t" + gt_str + ":" + ad_str1 + "," + ad_str2 + ":" + dp_str;
				}
			}
		}
	}

	for(i=0; i<allele_vec.size(); i++){
		allele_node = allele_vec.at(i);
		if(allele_node->reg) newVarVec.push_back(allele_node->reg);
	}

	for(i=0; i<allele_vec.size(); i++) delete allele_vec.at(i);
	vector<struct alleleNode*>().swap(allele_vec);
}

vector<reg_t*> varCand::rescueDupInvClipReg(vector<int32_t> &clusterId_incomplete){
	struct rescueSeqNode{
		string chrname, refseq;
		int64_t startRefPos, endRefPos;
		int32_t startQueryPos, endQueryPos: 28, aln_orient: 4;
		string qname, qseq, qual;
		vector<struct alnSeg*> alnSegs_left, alnSegs_right;

		int64_t startRefPos_clip, endRefPos_clip;
		int32_t startQueryPos_clip, endQueryPos_clip;
		string refseq_clip, qseq_clip;
	};

	vector<reg_t*> rescue_var_vec;
	int64_t startRefPos_cns, endRefPos_cns;
	size_t i, j, cluster_id;
	string qname, reg_str, reg_str_tmp, refseq, queryseq;
	char *p_seq;
	uint8_t *seq_int_whole;
	int32_t seq_len, bam_type;
	vector<vector<string>> querynames_vec;
	vector<string> qname_vec;
	vector<clipAlnData_t*> clipAlnDataVector, query_aln_segs;
	clipAlnData_t *left_aln_seg, *right_aln_seg, *clip_aln;
	int32_t margin_dist1, margin_dist2;
	int64_t left_qpos, right_qpos, left_rpos, right_rpos, sv_len, start_qpos_tmp, end_qpos_tmp;
	vector<vector<struct rescueSeqNode*>> rescue_seqs_vec;
	vector<struct rescueSeqNode*> re_seq_vec;
	struct rescueSeqNode *rescue_seq_node;
	reg_t *reg;
	vector<struct alnSeg*> alnSegs_left, alnSegs_right;
	struct alnSeg *aln_seg;
	int32_t idx_alnSegs_left, idx_alnSegs_right;
	vector<int32_t> sidemost_item_id_vec;

	ofstream outfile_rescue_reads, outfile_rescue_refseq, outfile_rescue_cns;
	string tmp_rescue_cnsfilename, cons_header, abpoa_cmd, minimap2_cmd, cmd_str;
	int32_t ret_status, status, serial_number, left_shift_size, right_shift_size, chrlen_tmp, query_orient, query_len_whole;
	bool flag, valid_flag;
	vector<minimap2_aln_t*> minimap2_aln_vec_rescue;
	double val;

	chrlen_tmp = faidx_seq_len64(fai, chrname.c_str()); // get reference size
	startRefPos_cns = leftClipRefPos - EXT_SIZE_CHK_VAR_LOC * 2;
	endRefPos_cns = rightClipRefPos + EXT_SIZE_CHK_VAR_LOC * 2;
	if(startRefPos_cns<1) startRefPos_cns = 1;
	if(endRefPos_cns>chrlen_tmp) endRefPos_cns = chrlen_tmp;
	left_shift_size = leftClipRefPos - startRefPos_cns;  // left shift size
	right_shift_size = endRefPos_cns - rightClipRefPos;  // right shift size

	if(sv_type==VAR_DUP or sv_type==VAR_INV){
		// load the clipping data
		clipAlnDataLoader data_loader(chrname, startRefPos_cns, endRefPos_cns, inBamFile, minClipEndSize, maxVarRegSize, minMapQ, minHighMapQ);
		data_loader.loadClipAlnDataWithSATag(clipAlnDataVector, max_ultra_high_cov);
		//data_loader.loadClipAlnDataWithSATag(clipAlnDataVector, 0); //deleted on 2024-03-22

		if(clipAlnDataVector.size() > 0){
			reg_str = chrname + ":" + to_string(startRefPos_cns) + "-" + to_string(endRefPos_cns);
			pthread_mutex_lock(&mutex_fai);
			p_seq = fai_fetch(fai, reg_str.c_str(), &seq_len);
			pthread_mutex_unlock(&mutex_fai);
			refseq = p_seq;
			free(p_seq);

			if(seq_len>0){
				bam_type = getBamType(clipAlnDataVector);
				if(bam_type!=BAM_INVALID){
					// get cluster information
					querynames_vec = getClusterInfo(clusterfilename);

					for(cluster_id=0; cluster_id<querynames_vec.size(); cluster_id++){
						if(find(clusterId_incomplete.begin(), clusterId_incomplete.end(), cluster_id) == clusterId_incomplete.end()){  // call complete, then skip
							continue;
						}

						qname_vec = querynames_vec.at(cluster_id);
						re_seq_vec.clear();

						//num = 0;
						for (i = 0; i < clipAlnDataVector.size(); i++) {//for each read query
							qname = clipAlnDataVector.at(i)->queryname;
							if(find(qname_vec.begin(), qname_vec.end(), qname) != qname_vec.end()){

								query_aln_segs = getQueryClipAlnSegsAll(qname, clipAlnDataVector);  // get query clip align segments
								//if(query_aln_segs.size()>MAX_ALN_SEG_NUM_PER_READ_TRA) { // ignore reads of too many align segments
								if(query_aln_segs.size()>(size_t)max_seg_num_per_read) { // ignore reads of too many align segments
									//cout << "clipReg: " << chrname << ":" << startRefPos << "-" << endRefPos << ", qname=" << queryname << ", align segment number=" << query_aln_segs.size() << endl;
									continue;
								}

//								if(qname.compare("SRR8858435.1.12303")==0){
//									cout << "qname=" << qname << endl;
//								}

								if (clipAlnDataVector.at(i)->query_checked_flag == false) {

									query_len_whole = query_orient = -1;
									for(j=0; j<query_aln_segs.size(); j++){
										clip_aln = query_aln_segs.at(j);
										if(clip_aln->leftHardClippedFlag==false and clip_aln->rightHardClippedFlag==false and clip_aln->bam){
											seq_int_whole = bam_get_seq(clip_aln->bam);
											query_len_whole = clip_aln->bam->core.l_qseq;
											query_orient = clip_aln->aln_orient;
											break;
										}
									}
									if(query_len_whole==-1) continue; // ignore all hard clip segments

									left_aln_seg = right_aln_seg = NULL;
									sidemost_item_id_vec = getLeftRightClipAlnId(startRefPos_cns, endRefPos_cns, query_aln_segs);
									if(sidemost_item_id_vec.size()>0){
										left_aln_seg = query_aln_segs.at(sidemost_item_id_vec.at(0));
										right_aln_seg = query_aln_segs.at(sidemost_item_id_vec.at(1));

										// alnSegs_left
										switch(bam_type){
											case BAM_CIGAR_NO_DIFF_MD:
												//alnSegs = generateAlnSegs(b);
												alnSegs_left = generateAlnSegs2(left_aln_seg->bam, startRefPos_cns, endRefPos_cns);
												//alnSegs_left_clip = generateAlnSegs2(left_aln_seg->bam, leftClipRefPos-1, rightClipRefPos);
												break;
											case BAM_CIGAR_NO_DIFF_NO_MD:
											case BAM_CIGAR_DIFF_MD:
											case BAM_CIGAR_DIFF_NO_MD:
												//alnSegs = generateAlnSegs_no_MD(b, refseq, startPos, endPos);
												alnSegs_left = generateAlnSegs_no_MD2(left_aln_seg->bam, refseq, startRefPos_cns, endRefPos_cns);
												//alnSegs_left_clip = generateAlnSegs_no_MD2(left_aln_seg->bam, refseq, leftClipRefPos-1, rightClipRefPos);
												break;
											default:
												cerr << __func__ << ": unknown bam type, error!" << endl;
												exit(1);
										}// generate align segments
										//if(clip_reg_flag==false)
											removeHeadTailAlnSegs(alnSegs_left, true);
										// only consider the query wholly , contain the variants
										idx_alnSegs_left = -1;
										for(j=0; j<alnSegs_left.size(); j++){
											aln_seg = alnSegs_left.at(j);
											if((aln_seg->opflag==BAM_CMATCH or aln_seg->opflag==BAM_CEQUAL)){
												if(aln_seg->startRpos<=leftClipRefPos-1){
													idx_alnSegs_left = j;
												}else{
													if(idx_alnSegs_left==-1 and aln_seg->startRpos-leftClipRefPos<=MAX_REF_DIST_SAME_CHR)
														idx_alnSegs_left = j;
													break;
												}
											}
										}

										// alnSegs_right
										switch(bam_type){
											case BAM_CIGAR_NO_DIFF_MD:
												//alnSegs = generateAlnSegs(b);
												alnSegs_right = generateAlnSegs2(right_aln_seg->bam, startRefPos_cns, endRefPos_cns);
												//alnSegs_right_clip = generateAlnSegs2(right_aln_seg->bam, leftClipRefPos-1, rightClipRefPos);
												break;
											case BAM_CIGAR_NO_DIFF_NO_MD:
											case BAM_CIGAR_DIFF_MD:
											case BAM_CIGAR_DIFF_NO_MD:
												//alnSegs = generateAlnSegs_no_MD(b, refseq, startPos, endPos);
												alnSegs_right = generateAlnSegs_no_MD2(right_aln_seg->bam, refseq, startRefPos_cns, endRefPos_cns);
												//alnSegs_right_clip = generateAlnSegs_no_MD2(right_aln_seg->bam, refseq, leftClipRefPos-1, rightClipRefPos);
												break;
											default:
												cerr << __func__ << ": unknown bam type, error!" << endl;
												exit(1);
										}// generate align segments
										//if(clip_reg_flag==false)
											removeHeadTailAlnSegs(alnSegs_right, true);
										idx_alnSegs_right = -1;
										for(j=0; j<alnSegs_right.size(); j++){
											aln_seg = alnSegs_right.at(j);
											if((aln_seg->opflag==BAM_CMATCH or aln_seg->opflag==BAM_CEQUAL)){
												if(aln_seg->startRpos<=rightClipRefPos){
													idx_alnSegs_right = j;
												}else{
													if(idx_alnSegs_right==-1 and aln_seg->startRpos-rightClipRefPos<=MAX_REF_DIST_SAME_CHR)
														idx_alnSegs_right = j;
													break;
												}
											}
										}

										if(alnSegs_left.size()>0 and alnSegs_right.size()>0){

											left_qpos = alnSegs_left.at(0)->startQpos;
											right_qpos = getEndQueryPosAlnSeg(alnSegs_right.at(alnSegs_right.size()-1)->startQpos, alnSegs_right.at(alnSegs_right.size()-1)->opflag, alnSegs_right.at(alnSegs_right.size()-1)->seglen);
											if(right_qpos>query_len_whole) right_qpos = query_len_whole;

											valid_flag = true;
											queryseq = "";
											if(left_aln_seg->aln_orient==query_orient){
												if(left_qpos<=right_qpos) for(j=left_qpos-1; j<(size_t)right_qpos; j++) queryseq += "=ACMGRSVTWYHKDBN"[bam_seqi(seq_int_whole, j)];  // seq
												else valid_flag = false;
											}else{
												start_qpos_tmp = query_len_whole - right_qpos;
												end_qpos_tmp = query_len_whole - left_qpos;
												if(start_qpos_tmp<=end_qpos_tmp) {
													for(j=start_qpos_tmp-1; j<(size_t)end_qpos_tmp; j++) queryseq += "=ACMGRSVTWYHKDBN"[bam_seqi(seq_int_whole, j)];  // seq
													reverseComplement(queryseq);
												}else valid_flag = false;
											}

											if(valid_flag){
												rescue_seq_node = new struct rescueSeqNode();
												rescue_seq_node->chrname = chrname;
												rescue_seq_node->startRefPos = alnSegs_left.at(0)->startRpos;
												rescue_seq_node->endRefPos = getEndRefPosAlnSeg(alnSegs_right.at(alnSegs_right.size()-1)->startRpos, alnSegs_right.at(alnSegs_right.size()-1)->opflag, alnSegs_right.at(alnSegs_right.size()-1)->seglen);
												rescue_seq_node->startQueryPos = left_qpos;
												rescue_seq_node->endQueryPos = right_qpos;
												rescue_seq_node->aln_orient = left_aln_seg->aln_orient;
												rescue_seq_node->qname = qname;
												rescue_seq_node->qseq = queryseq;
												rescue_seq_node->refseq = refseq;
												rescue_seq_node->alnSegs_left = alnSegs_left;
												rescue_seq_node->alnSegs_right = alnSegs_right;


												if(idx_alnSegs_left!=-1 and idx_alnSegs_right!=-1){
													margin_dist1 = leftClipRefPos - 1 - alnSegs_left.at(idx_alnSegs_left)->startRpos;
													if(margin_dist1<0) margin_dist1 = 0;
													left_rpos = alnSegs_left.at(idx_alnSegs_left)->startRpos + margin_dist1;
													left_qpos = alnSegs_left.at(idx_alnSegs_left)->startQpos + margin_dist1;
													if(left_qpos>query_len_whole) left_qpos = query_len_whole;

													margin_dist2 = rightClipRefPos - alnSegs_right.at(idx_alnSegs_right)->startRpos;
													if(margin_dist2<0) margin_dist2 = 0;
													right_rpos = alnSegs_right.at(idx_alnSegs_right)->startRpos + margin_dist2;
													right_qpos = alnSegs_right.at(idx_alnSegs_right)->startQpos + margin_dist2;
													if(right_qpos>query_len_whole) right_qpos = query_len_whole;

													queryseq = "";
													if(left_aln_seg->aln_orient==query_orient){
														if(left_qpos<=right_qpos) for(j=left_qpos-1; j<(size_t)right_qpos; j++) queryseq += "=ACMGRSVTWYHKDBN"[bam_seqi(seq_int_whole, j)];  // seq
														else valid_flag = false;
													}else{
														start_qpos_tmp = query_len_whole - right_qpos;
														end_qpos_tmp = query_len_whole - left_qpos;
														if(start_qpos_tmp<1) start_qpos_tmp = 1;
														if(end_qpos_tmp>query_len_whole) end_qpos_tmp = query_len_whole;
														if(start_qpos_tmp<=end_qpos_tmp) {
															for(j=start_qpos_tmp-1; j<(size_t)end_qpos_tmp; j++) queryseq += "=ACMGRSVTWYHKDBN"[bam_seqi(seq_int_whole, j)];  // seq
															reverseComplement(queryseq);
														}else valid_flag = false;
													}

													if(valid_flag){
														rescue_seq_node->startRefPos_clip = left_rpos;
														rescue_seq_node->endRefPos_clip = right_rpos;
														rescue_seq_node->startQueryPos_clip = left_qpos;
														rescue_seq_node->endQueryPos_clip = right_qpos;
														rescue_seq_node->qseq_clip = queryseq;

														reg_str_tmp = chrname + ":" + to_string(rescue_seq_node->startRefPos_clip) + "-" + to_string(rescue_seq_node->endRefPos_clip);
														pthread_mutex_lock(&mutex_fai);
														p_seq = fai_fetch(fai, reg_str_tmp.c_str(), &seq_len);
														pthread_mutex_unlock(&mutex_fai);
														rescue_seq_node->refseq_clip = p_seq;
														free(p_seq);
													}else {
														delete rescue_seq_node;
													}
												}
												if(valid_flag) re_seq_vec.push_back(rescue_seq_node);
											}
										}
									}
									for(j=0; j<query_aln_segs.size(); j++) query_aln_segs.at(j)->query_checked_flag = true;
								}
							}
						}
						if(re_seq_vec.size()>0) rescue_seqs_vec.push_back(re_seq_vec);
					}
				}
			}
		}

		if(rightClipRefPos-leftClipRefPos<=MAX_RESCUE_VAR_SIZE) {
			if(rescue_seqs_vec.size()>0){
				//rescue_readsfilename = out_dir_call + "/rescue_reads_" + chrname + "_" + to_string(leftClipRefPos) + "-" + to_string(rightClipRefPos) + ".fa";
				rescue_refseqfilename = out_dir_call + "/rescue_refseq_" + chrname + "_" + to_string(leftClipRefPos) + "-" + to_string(rightClipRefPos) + ".fa";
				rescue_cnsfilename = out_dir_call + "/rescue_cns_" + chrname + "_" + to_string(leftClipRefPos) + "-" + to_string(rightClipRefPos) + ".fa";
				//tmp_rescue_cnsfilename = out_dir_call + "/tmp_rescue_cns_" + chrname + "_" + to_string(leftClipRefPos) + "-" + to_string(rightClipRefPos) + ".fa";
				rescue_alnfilename = out_dir_call + "/rescue_minimap2_" + chrname + "_" + to_string(leftClipRefPos) + "-" + to_string(rightClipRefPos) + ".paf";

				outfile_rescue_refseq.open(rescue_refseqfilename);
				if(outfile_rescue_refseq.is_open()==false){
					cerr << "In " << __func__ << "(), cannot open file '" << rescue_refseqfilename << "', error!" << endl;
					exit(1);
				}
				outfile_rescue_cns.open(rescue_cnsfilename);
				if(outfile_rescue_cns.is_open()==false){
					cerr << "In " << __func__ << "(), cannot open file '" << rescue_cnsfilename << "', error!" << endl;
					exit(1);
				}

				// refseq file
				cons_header = ">" + reg_str + "___" + to_string(left_shift_size) + "___" + to_string(right_shift_size) + "___" + rescue_refseqfilename;
				outfile_rescue_refseq << cons_header << endl;
				outfile_rescue_refseq << refseq << endl;
				outfile_rescue_refseq.close();
			}

			serial_number = 1;
			for(cluster_id=0; cluster_id<rescue_seqs_vec.size(); cluster_id++){
				rescue_readsfilename = out_dir_call + "/rescue_reads_" + chrname + "_" + to_string(leftClipRefPos) + "-" + to_string(rightClipRefPos) + "_" + to_string(cluster_id) + ".fa";
				tmp_rescue_cnsfilename = out_dir_call + "/tmp_rescue_cns_" + chrname + "_" + to_string(leftClipRefPos) + "-" + to_string(rightClipRefPos) + "_" + to_string(cluster_id) + ".fa";

				re_seq_vec = rescue_seqs_vec.at(cluster_id);
				if(re_seq_vec.size()>0){
					outfile_rescue_reads.open(rescue_readsfilename);
					if(outfile_rescue_reads.is_open()==false){
						cerr << "In " << __func__ << "(), cannot open file '" << rescue_readsfilename << "', error!" << endl;
						exit(1);
					}

					// reads file
					for(i=0; i<re_seq_vec.size(); i++){
						rescue_seq_node = re_seq_vec.at(i);
			//			cout << "[" << i << "]: " << rescue_seq_node->chrname << ":" << rescue_seq_node->startRefPos << "-" << rescue_seq_node->endRefPos << endl;
			//			cout << rescue_seq_node->refseq << endl;
						outfile_rescue_reads << ">" << rescue_seq_node->qname << ":" << rescue_seq_node->startQueryPos << "-" << rescue_seq_node->endQueryPos << ":" << rescue_seq_node->aln_orient << endl;
						outfile_rescue_reads << rescue_seq_node->qseq << endl;
					}
					outfile_rescue_reads.close();

					// perform abpoa consensus
					abpoa_cmd = "abpoa -S -o " + tmp_rescue_cnsfilename + " " + rescue_readsfilename + " > /dev/null 2>&1";
					system(abpoa_cmd.c_str());

					flag = isFileExist(tmp_rescue_cnsfilename);
					if(flag){ // cons generated successfully
						// reorganize the consensus sequence
						FastaSeqLoader fa_loader(tmp_rescue_cnsfilename);
						for(i=0; i<fa_loader.getFastaSeqCount(); i++){
							cons_header = ">abpoa_rescueCns_";
							cons_header += to_string(serial_number) + "_"; //serial number
							cons_header += to_string(fa_loader.getFastaSeqLen(i)) + "_"; //consensus sequence length
							cons_header += to_string(re_seq_vec.size()) + "___";  //reads count

							for(j=0; j<re_seq_vec.size()-1; j++) cons_header += re_seq_vec.at(j)->qname + ";";
							cons_header += re_seq_vec.at(re_seq_vec.size()-1)->qname;

							outfile_rescue_cns << cons_header << endl; // header
							outfile_rescue_cns << fa_loader.getFastaSeq(i) << endl;  // seq

							serial_number ++;
						}
					}

					cmd_str = "rm -rf " + rescue_readsfilename;
					system(cmd_str.c_str());  // remove reads file
				}
			}
			outfile_rescue_cns.close();

			flag = isFileExist(rescue_cnsfilename);
			if(flag){
				minimap2_cmd = "minimap2 -c -x asm20 -o " + rescue_alnfilename + " " + rescue_refseqfilename + " " + rescue_cnsfilename + " > /dev/null 2>&1";
				status = system(minimap2_cmd.c_str());
				ret_status = getSuccessStatusSystemCmd(status);
				if(ret_status==0){ // command executed successfully
					// parse alignment information
					minimap2_aln_vec_rescue = minimap2Parse(rescue_alnfilename, rescue_cnsfilename, rescue_refseqfilename);

					// get reference sequence and query sequence
					rescue_var_vec = computeClipRegVarLoc(rescue_alnfilename, rescue_cnsfilename, rescue_refseqfilename, clusterfilename, minimap2_aln_vec_rescue, &clusterId_incomplete, true);
				}
			}
		}

		if(rescue_var_vec.size()==0){ // failed, then pick the first query sequence
			for(cluster_id=0; cluster_id<rescue_seqs_vec.size(); cluster_id++){
				re_seq_vec = rescue_seqs_vec.at(cluster_id);
				for(i=0; i<re_seq_vec.size(); i++){
					rescue_seq_node = re_seq_vec.at(i);
					//if(rescue_seq_node->startRefPos==leftClipRefPos and rescue_seq_node->endRefPos==rightClipRefPos){
					if(rescue_seq_node->qseq_clip.size()>0 and rescue_seq_node->refseq_clip.size()>0){
						sv_len = rescue_seq_node->qseq_clip.size() - rescue_seq_node->refseq_clip.size();
						if(sv_type==VAR_DUP){ // added on 2025-08-20
							val = (double)abs(sv_len) / (rightClipRefPos - leftClipRefPos + 1);
							if(val>=VALID_DUP_INV_SIZE_RATIO_CLIP_REG){
								// allocate new node
								reg = new reg_t();
								reg->var_type = sv_type;
								reg->chrname = rescue_seq_node->chrname;
								reg->startRefPos = rescue_seq_node->startRefPos_clip;
								reg->endRefPos = rescue_seq_node->endRefPos_clip;
								reg->startLocalRefPos = -1;
								reg->endLocalRefPos = -1;
								reg->startQueryPos = rescue_seq_node->startQueryPos_clip;
								reg->endQueryPos = rescue_seq_node->endQueryPos_clip;
								reg->query_id = -1;
								reg->sv_len = sv_len;
								//if(sv_len<0) reg->sv_len = -reg->sv_len;
								reg->minimap2_aln_id = -1;
								reg->call_success_status = true;
								reg->short_sv_flag = false;
								reg->zero_cov_flag = false;
								reg->aln_seg_end_flag = false;
								reg->query_pos_invalid_flag = false;
								reg->merge_flag = false;
								reg->rescue_flag = true;
								reg->aln_orient = rescue_seq_node->aln_orient;
								reg->gt_type = -1;
								reg->gt_seq = "";
								reg->supp_num = re_seq_vec.size();
								reg->DP = reg->AF = 0;
								reg->discover_level = VAR_DISCOV_L_READS;
								reg->refseq =  rescue_seq_node->refseq_clip;
								reg->altseq = rescue_seq_node->qseq_clip;

								switch(sv_type){
									case VAR_DUP: // DUP
										reg->dup_num = dup_num;
										reg->large_indel_flag = false;
										break;
									case VAR_INV: // INV
										reg->dup_num = 0;
										reg->large_indel_flag = false;
										break;
									case VAR_INS: // large INS or DEL
									case VAR_DEL:
										reg->dup_num = 0;
										reg->large_indel_flag = true;
										break;
								}

								rescue_var_vec.push_back(reg); // pick the first one
								break;
							}
						}
					}
				}
			}
		}

		// release memory
		for(cluster_id=0; cluster_id<rescue_seqs_vec.size(); cluster_id++){
			for(i=0; i<re_seq_vec.size(); i++){
				rescue_seq_node = re_seq_vec.at(i);
				destroyAlnSegs(rescue_seq_node->alnSegs_left);
				destroyAlnSegs(rescue_seq_node->alnSegs_right);
				//destroyAlnSegs(rescue_seq_node->alnSegs_left_clip);
				//destroyAlnSegs(rescue_seq_node->alnSegs_right_clip);
				delete rescue_seq_node;
			}
			vector<struct rescueSeqNode*>().swap(re_seq_vec);
		}
		vector<vector<struct rescueSeqNode*>>().swap(rescue_seqs_vec);

		if(!clipAlnDataVector.empty()) destoryClipAlnData(clipAlnDataVector);
	}

	return rescue_var_vec;
}

vector<reg_t*> varCand::rescueLargeIndelClipReg(vector<int32_t> &clusterId_incomplete){
	struct rescueSeqNode{
		string chrname, refseq;
		int64_t startRefPos, endRefPos;
		int32_t startQueryPos, endQueryPos: 28, aln_orient: 4;
		string qname, qseq, qual;
		vector<struct alnSeg*> alnSegs_left, alnSegs_right;

		int64_t startRefPos_clip, endRefPos_clip;
		int32_t startQueryPos_clip, endQueryPos_clip;
		string refseq_clip, qseq_clip;
		//vector<struct alnSeg*> alnSegs_left_clip, alnSegs_right_clip;
	};

	vector<reg_t*> rescue_var_vec;
	int64_t startRefPos_cns, endRefPos_cns;
	size_t i, j, cluster_id;
	string qname, reg_str, reg_str_tmp, refseq, queryseq;
	char *p_seq;
	uint8_t *seq_int_whole;
	int32_t k, seq_len, bam_type, query_len_whole, query_orient, new_sv_type;
	vector<vector<string>> querynames_vec;
	vector<string> qname_vec;
	vector<clipAlnData_t*> clipAlnDataVector, query_aln_segs;
	clipAlnData_t *clip_aln_seg, *mate_clip_aln_seg, *left_aln_seg, *right_aln_seg, *clip_aln;
	vector<int32_t> adjClipAlnSegInfo;
	int32_t mate_arr_idx, clip_end_flag, mate_clip_end_flag, dist, margin_dist1, margin_dist2, increase_direction;
	int32_t mate_arr_idx_same_chr, mate_clip_end_flag_same_chr, min_dist, min_ref_dist, min_dist_same_chr, min_ref_dist_same_chr, seg_skip_num, ref_skip_num, orient_factor;
	int64_t clip_refpos, clip_querypos, mate_clip_refpos, mate_clip_querypos;
	int64_t left_qpos, right_qpos, left_rpos, right_rpos, sv_len, start_qpos_tmp, end_qpos_tmp;
	bool valid_query_end_flag, self_overlap_flag, seg_skip_flag, ref_skip_flag;
	vector<vector<struct rescueSeqNode*>> rescue_seqs_vec;
	vector<struct rescueSeqNode*> re_seq_vec;
	struct rescueSeqNode *rescue_seq_node;
	reg_t *reg;
	vector<struct alnSeg*> alnSegs_left, alnSegs_right;
	struct alnSeg *aln_seg;
	int32_t idx_alnSegs_left, idx_alnSegs_right;

	ofstream outfile_rescue_reads, outfile_rescue_refseq, outfile_rescue_cns;
	string tmp_rescue_cnsfilename, cons_header, abpoa_cmd, minimap2_cmd, cmd_str;
	int32_t ret_status, status, serial_number, left_shift_size, right_shift_size, chrlen_tmp;
	bool flag, valid_flag;
	vector<minimap2_aln_t*> minimap2_aln_vec_rescue;
	double val;

	chrlen_tmp = faidx_seq_len64(fai, chrname.c_str()); // get reference size
	startRefPos_cns = leftClipRefPos - VAR_ALN_EXTEND_SIZE * 2;
	endRefPos_cns = rightClipRefPos + VAR_ALN_EXTEND_SIZE * 2;
	if(startRefPos_cns<1) startRefPos_cns = 1;
	if(endRefPos_cns>chrlen_tmp) endRefPos_cns = chrlen_tmp;
	left_shift_size = leftClipRefPos - startRefPos_cns;  // left shift size
	right_shift_size = endRefPos_cns - rightClipRefPos;  // right shift size

	if(sv_type==VAR_INS or sv_type==VAR_DEL){
		// load the clipping data

		clipAlnDataLoader data_loader(chrname, startRefPos_cns, endRefPos_cns, inBamFile, minClipEndSize, maxVarRegSize, minMapQ, minHighMapQ);
		data_loader.loadClipAlnDataWithSATag(clipAlnDataVector, max_ultra_high_cov);
		//data_loader.loadClipAlnDataWithSATag(clipAlnDataVector, 0); //deleted on 2024-03-22

		if(clipAlnDataVector.size() > 0){
			reg_str = chrname + ":" + to_string(startRefPos_cns) + "-" + to_string(endRefPos_cns);
			pthread_mutex_lock(&mutex_fai);
			p_seq = fai_fetch(fai, reg_str.c_str(), &seq_len);
			pthread_mutex_unlock(&mutex_fai);
			refseq = p_seq;
			free(p_seq);

			if(seq_len>0){
				bam_type = getBamType(clipAlnDataVector);
				if(bam_type!=BAM_INVALID) {
					// get cluster information
					querynames_vec = getClusterInfo(clusterfilename);

					for(cluster_id=0; cluster_id<querynames_vec.size(); cluster_id++){
						if(find(clusterId_incomplete.begin(), clusterId_incomplete.end(), cluster_id) == clusterId_incomplete.end()){  // call complete, then skip
							continue;
						}

						qname_vec = querynames_vec.at(cluster_id);
						re_seq_vec.clear();

						//num = 0;
						for (i = 0; i < clipAlnDataVector.size(); i++) {//for each read query
							qname = clipAlnDataVector.at(i)->queryname;
							if(find(qname_vec.begin(), qname_vec.end(), qname) != qname_vec.end()){

								query_aln_segs = getQueryClipAlnSegsAll(qname, clipAlnDataVector);  // get query clip align segments
								//if(query_aln_segs.size()>MAX_ALN_SEG_NUM_PER_READ_TRA) { // ignore reads of too many align segments
								if(query_aln_segs.size()>(size_t)max_seg_num_per_read) { // ignore reads of too many align segments
									//cout << "clipReg: " << chrname << ":" << startRefPos << "-" << endRefPos << ", qname=" << queryname << ", align segment number=" << query_aln_segs.size() << endl;
									continue;
								}

	//							if(qname.compare("SRR8955953.1886784")==0){
	//								cout << "qname=" << qname << endl;
	//							}

								if (clipAlnDataVector.at(i)->query_checked_flag == false) {

									query_len_whole = query_orient = -1;
									for(j=0; j<query_aln_segs.size(); j++){
										clip_aln = query_aln_segs.at(j);
										if(clip_aln->leftHardClippedFlag==false and clip_aln->rightHardClippedFlag==false and clip_aln->bam){
											seq_int_whole = bam_get_seq(clip_aln->bam);
											query_len_whole = clip_aln->bam->core.l_qseq;
											query_orient = clip_aln->aln_orient;
											break;
										}
									}
									if(query_len_whole==-1) continue; // ignore all hard clip segments

									for(j=0; j<query_aln_segs.size(); j++){
										clip_aln_seg = query_aln_segs.at(j);
										valid_query_end_flag = false;
										clip_end_flag = -1;
										if(clip_aln_seg->left_clip_checked_flag==false and clip_aln_seg->leftClipSize>=minClipEndSize and clip_aln_seg->startRefPos>=startRefPos_cns and clip_aln_seg->startRefPos<=endRefPos_cns and clip_aln_seg->chrname.compare(chrname)==0){ // left end
											valid_query_end_flag = true;
											clip_refpos = clip_aln_seg->startRefPos;
											clip_querypos = clip_aln_seg->startQueryPos;
											clip_aln_seg->left_clip_checked_flag = true;
											clip_end_flag = LEFT_END;
										}else if(clip_aln_seg->right_clip_checked_flag==false and clip_aln_seg->rightClipSize>=minClipEndSize and clip_aln_seg->endRefPos>=startRefPos_cns and clip_aln_seg->endRefPos<=endRefPos_cns and clip_aln_seg->chrname.compare(chrname)==0){ // right end
											valid_query_end_flag = true;
											clip_refpos = clip_aln_seg->endRefPos;
											clip_querypos = clip_aln_seg->endQueryPos;
											clip_aln_seg->right_clip_checked_flag = true;
											clip_end_flag = RIGHT_END;
										}

										if(valid_query_end_flag){

											// deal with the mate clip end
											seg_skip_flag = ref_skip_flag = false;
											adjClipAlnSegInfo = getAdjacentClipAlnSeg(j, clip_end_flag, query_aln_segs, minClipEndSize, maxVarRegSize, minMapQ);
											mate_arr_idx = adjClipAlnSegInfo.at(0);
											mate_clip_end_flag = adjClipAlnSegInfo.at(1);
											min_dist = adjClipAlnSegInfo.at(2);
											min_ref_dist = adjClipAlnSegInfo.at(3);
											mate_arr_idx_same_chr = adjClipAlnSegInfo.at(4);
											mate_clip_end_flag_same_chr = adjClipAlnSegInfo.at(5);
											min_dist_same_chr = adjClipAlnSegInfo.at(6);
											min_ref_dist_same_chr = adjClipAlnSegInfo.at(7);
											increase_direction = adjClipAlnSegInfo.at(8);
											if(mate_arr_idx!=-1){ // mated
												if(mate_arr_idx_same_chr!=-1){
													mate_clip_aln_seg = query_aln_segs.at(mate_arr_idx_same_chr);
													if((mate_clip_end_flag_same_chr==LEFT_END and mate_clip_aln_seg->left_aln) or (mate_clip_end_flag_same_chr==RIGHT_END and mate_clip_aln_seg->right_aln)){
														self_overlap_flag = isSegSelfOverlap(clip_aln_seg, mate_clip_aln_seg, MAX_SV_SIZE_USR);
														dist = min_dist_same_chr - min_ref_dist_same_chr;
														if(increase_direction==1) dist = -dist;

														// prefer the segments on the same chromosome
														//if((mate_arr_idx!=mate_arr_idx_same_chr and (abs(min_ref_dist_same_chr)<=MAX_REF_DIST_SAME_CHR or self_overlap_flag)) or (mate_arr_idx==mate_arr_idx_same_chr and abs(min_dist_same_chr)>MAX_REF_DIST_SAME_CHR and abs(dist)>MAX_REF_DIST_SAME_CHR)){ // different chromosomes with near reference distance
														if((mate_arr_idx!=mate_arr_idx_same_chr and (abs(min_ref_dist_same_chr)<=MAX_REF_DIST_SAME_CHR or self_overlap_flag)) or (mate_arr_idx==mate_arr_idx_same_chr and abs(min_dist_same_chr)>MAX_REF_DIST_SAME_CHR and dist>MAX_REF_DIST_SAME_CHR)){ // different chromosomes with near reference distance
														//if(abs(min_ref_dist_same_chr)<=MAX_REF_DIST_SAME_CHR){ // different chromosomes with near reference distance
															mate_arr_idx = mate_arr_idx_same_chr;
															mate_clip_end_flag = mate_clip_end_flag_same_chr;
															seg_skip_num ++;
															seg_skip_flag = true;
														}else if(mate_arr_idx==mate_arr_idx_same_chr and mate_clip_end_flag==mate_clip_end_flag_same_chr and self_overlap_flag==false and abs(min_dist_same_chr)<=MAX_REF_DIST_SAME_CHR and abs(min_ref_dist_same_chr)>MAX_REF_DIST_SAME_CHR){
															mate_arr_idx = mate_arr_idx_same_chr;
															mate_clip_end_flag = mate_clip_end_flag_same_chr;
															ref_skip_num ++;
															ref_skip_flag = true;
														}else{
															if(abs(min_dist_same_chr)>=MAX_REF_DIST_SAME_CHR and abs(min_dist_same_chr)<=maxVarRegSize and abs(min_ref_dist_same_chr)>MAX_REF_DIST_SAME_CHR and abs(min_ref_dist_same_chr)<=maxVarRegSize){
																if(clip_aln_seg->aln_orient==mate_clip_aln_seg->aln_orient and abs(min_dist)>MAX_REF_DIST_SAME_CHR and abs(min_ref_dist)>MAX_REF_DIST_SAME_CHR){
																	if(dist>0){
																		mate_arr_idx = mate_arr_idx_same_chr;
																		mate_clip_end_flag = mate_clip_end_flag_same_chr;
																		seg_skip_num ++;
																		seg_skip_flag = true;
																	}else{
																		mate_arr_idx = mate_arr_idx_same_chr;
																		mate_clip_end_flag = mate_clip_end_flag_same_chr;
																		ref_skip_num ++;
																		ref_skip_flag = true;
																	}
																}
															}
														}
													}
												}

												if(seg_skip_flag or ref_skip_flag){
													mate_clip_aln_seg = query_aln_segs.at(mate_arr_idx_same_chr);
													if(mate_clip_end_flag==LEFT_END){
														mate_clip_refpos = mate_clip_aln_seg->startRefPos;
														mate_clip_querypos = mate_clip_aln_seg->startQueryPos;
														mate_clip_aln_seg->left_clip_checked_flag = true;
													}else{
														mate_clip_refpos = mate_clip_aln_seg->startRefPos;
														mate_clip_querypos = mate_clip_aln_seg->startQueryPos;
														mate_clip_aln_seg->right_clip_checked_flag = true;
													}

													orient_factor = 0;
													if(clip_aln_seg->aln_orient==ALN_PLUS_ORIENT and mate_clip_aln_seg->aln_orient==ALN_PLUS_ORIENT) orient_factor = 1;
													else if(clip_aln_seg->aln_orient==ALN_MINUS_ORIENT and mate_clip_aln_seg->aln_orient==ALN_MINUS_ORIENT) orient_factor = -1;

													left_aln_seg = right_aln_seg = NULL;
													if(clip_end_flag==RIGHT_END and mate_clip_end_flag==LEFT_END and clip_aln_seg->right_aln and mate_clip_aln_seg->left_aln){
														dist = (mate_clip_querypos - clip_querypos) * orient_factor - (mate_clip_refpos - clip_refpos);
														left_aln_seg = clip_aln_seg;
														right_aln_seg = mate_clip_aln_seg;
													}else if(clip_end_flag==LEFT_END and mate_clip_end_flag==RIGHT_END and clip_aln_seg->left_aln and mate_clip_aln_seg->right_aln){
														dist = (clip_querypos - mate_clip_querypos) * orient_factor - (clip_refpos - mate_clip_refpos);
	//													left_aln_seg = clip_aln_seg;
	//													right_aln_seg = mate_clip_aln_seg;
														left_aln_seg = mate_clip_aln_seg;
														right_aln_seg = clip_aln_seg;
													}

													if(left_aln_seg and right_aln_seg and left_aln_seg->bam and right_aln_seg->bam and ((left_aln_seg->leftHardClippedFlag==false and left_aln_seg->rightHardClippedFlag==false) or (right_aln_seg->leftHardClippedFlag==false and right_aln_seg->rightHardClippedFlag==false))){
														// alnSegs_left
														switch(bam_type){
															case BAM_CIGAR_NO_DIFF_MD:
																//alnSegs = generateAlnSegs(b);
																alnSegs_left = generateAlnSegs2(left_aln_seg->bam, startRefPos_cns, endRefPos_cns);
																//alnSegs_left_clip = generateAlnSegs2(left_aln_seg->bam, leftClipRefPos-1, rightClipRefPos);
																break;
															case BAM_CIGAR_NO_DIFF_NO_MD:
															case BAM_CIGAR_DIFF_MD:
															case BAM_CIGAR_DIFF_NO_MD:
																//alnSegs = generateAlnSegs_no_MD(b, refseq, startPos, endPos);
																alnSegs_left = generateAlnSegs_no_MD2(left_aln_seg->bam, refseq, startRefPos_cns, endRefPos_cns);
																//alnSegs_left_clip = generateAlnSegs_no_MD2(left_aln_seg->bam, refseq, leftClipRefPos-1, rightClipRefPos);
																break;
															default:
																cerr << __func__ << ": unknown bam type, error!" << endl;
																exit(1);
														}// generate align segments
														removeHeadTailAlnSegs(alnSegs_left, true);
														if(alnSegs_left.size()==0){
															cout << "alnSegs_left.size=" << alnSegs_left.size() << ", error!" << endl;
															exit(1);
														}

														// only consider the query wholly contain the variants
														idx_alnSegs_left = -1;
														for(k=0; k<(int64_t)alnSegs_left.size(); k++){
															aln_seg = alnSegs_left.at(k);
															if((aln_seg->opflag==BAM_CMATCH or aln_seg->opflag==BAM_CEQUAL)){
																if(aln_seg->startRpos<=leftClipRefPos-1){
																	idx_alnSegs_left = k;
																}else{
																	if(idx_alnSegs_left==-1 and aln_seg->startRpos-leftClipRefPos<=MAX_REF_DIST_SAME_CHR)
																		idx_alnSegs_left = k;
																	break;
																}
															}
														}

														// alnSegs_right
														switch(bam_type){
															case BAM_CIGAR_NO_DIFF_MD:
																alnSegs_right = generateAlnSegs2(right_aln_seg->bam, startRefPos_cns, endRefPos_cns);
																break;
															case BAM_CIGAR_NO_DIFF_NO_MD:
															case BAM_CIGAR_DIFF_MD:
															case BAM_CIGAR_DIFF_NO_MD:
																alnSegs_right = generateAlnSegs_no_MD2(right_aln_seg->bam, refseq, startRefPos_cns, endRefPos_cns);
																break;
															default:
																cerr << __func__ << ": unknown bam type, error!" << endl;
																exit(1);
														}// generate align segments
														removeHeadTailAlnSegs(alnSegs_right, true);
														if(alnSegs_right.size()==0){
															cout << "alnSegs_right.size=" << alnSegs_right.size() << ", error!" << endl;
															exit(1);
														}

														idx_alnSegs_right = -1;
														for(k=0; k<(int64_t)alnSegs_right.size(); k++){
															aln_seg = alnSegs_right.at(k);
															if((aln_seg->opflag==BAM_CMATCH or aln_seg->opflag==BAM_CEQUAL)){
																if(aln_seg->startRpos<=rightClipRefPos){
																	idx_alnSegs_right = k;
																}else{
																	if(idx_alnSegs_right==-1 and aln_seg->startRpos-rightClipRefPos<=MAX_REF_DIST_SAME_CHR)
																		idx_alnSegs_right = k;
																	break;
																}
															}
														}

														left_rpos = alnSegs_left.at(0)->startRpos;
														right_rpos = getEndRefPosAlnSeg(alnSegs_right.at(alnSegs_right.size()-1)->startRpos, alnSegs_right.at(alnSegs_right.size()-1)->opflag, alnSegs_right.at(alnSegs_right.size()-1)->seglen);
														left_qpos = alnSegs_left.at(0)->startQpos;
														right_qpos = getEndQueryPosAlnSeg(alnSegs_right.at(alnSegs_right.size()-1)->startQpos, alnSegs_right.at(alnSegs_right.size()-1)->opflag, alnSegs_right.at(alnSegs_right.size()-1)->seglen);

	//													queryseq = "";
	//													if(left_aln_seg->bam and left_aln_seg->leftHardClippedFlag==false and left_aln_seg->rightHardClippedFlag==false)
	//														seq_int = bam_get_seq(left_aln_seg->bam);
	//													else if(right_aln_seg->bam and right_aln_seg->leftHardClippedFlag==false and right_aln_seg->rightHardClippedFlag==false)
	//														seq_int = bam_get_seq(right_aln_seg->bam);
	//													for(k=left_qpos-1; k<right_qpos; k++) queryseq += "=ACMGRSVTWYHKDBN"[bam_seqi(seq_int, k)];  // seq

														valid_flag = true;
														queryseq = "";
														if(left_aln_seg->aln_orient==query_orient){
															if(left_qpos<=right_qpos) for(j=left_qpos-1; j<(size_t)right_qpos; j++) queryseq += "=ACMGRSVTWYHKDBN"[bam_seqi(seq_int_whole, j)];  // seq
															else valid_flag = false;
														}else{
															start_qpos_tmp = query_len_whole - right_qpos;
															end_qpos_tmp = query_len_whole - left_qpos;
															if(start_qpos_tmp<=end_qpos_tmp){
																for(j=start_qpos_tmp-1; j<(size_t)end_qpos_tmp; j++) queryseq += "=ACMGRSVTWYHKDBN"[bam_seqi(seq_int_whole, j)];  // seq
																reverseComplement(queryseq);
															}else valid_flag = false;
														}

														if(valid_flag){
															rescue_seq_node = new struct rescueSeqNode();
															rescue_seq_node->chrname = chrname;
															rescue_seq_node->startRefPos = alnSegs_left.at(0)->startRpos;
															rescue_seq_node->endRefPos = getEndRefPosAlnSeg(alnSegs_right.at(alnSegs_right.size()-1)->startRpos, alnSegs_right.at(alnSegs_right.size()-1)->opflag, alnSegs_right.at(alnSegs_right.size()-1)->seglen);
															rescue_seq_node->startQueryPos = left_qpos;
															rescue_seq_node->endQueryPos = right_qpos;
															rescue_seq_node->aln_orient = clip_aln_seg->aln_orient;
															rescue_seq_node->qname = qname;
															rescue_seq_node->qseq = queryseq;
															rescue_seq_node->refseq = refseq;
															rescue_seq_node->alnSegs_left = alnSegs_left;
															rescue_seq_node->alnSegs_right = alnSegs_right;

															if(idx_alnSegs_left!=-1 and idx_alnSegs_right!=-1){
																margin_dist1 = leftClipRefPos - 1 - alnSegs_left.at(idx_alnSegs_left)->startRpos;
																if(margin_dist1<0) margin_dist1 = 0;
																left_rpos = alnSegs_left.at(idx_alnSegs_left)->startRpos + margin_dist1;
																left_qpos = alnSegs_left.at(idx_alnSegs_left)->startQpos + margin_dist1;

																margin_dist2 = rightClipRefPos - alnSegs_right.at(idx_alnSegs_right)->startRpos;
																if(margin_dist2<0) margin_dist2 = 0;
																right_rpos = alnSegs_right.at(idx_alnSegs_right)->startRpos + margin_dist2;
																right_qpos = alnSegs_right.at(idx_alnSegs_right)->startQpos + margin_dist2;
																end_qpos_tmp = getEndQueryPosAlnSeg(alnSegs_right.at(idx_alnSegs_right)->startQpos, alnSegs_right.at(idx_alnSegs_right)->opflag, alnSegs_right.at(idx_alnSegs_right)->seglen);

																if(left_rpos<=right_rpos and alnSegs_left.at(idx_alnSegs_left)->startRpos<=leftClipRefPos and end_qpos_tmp>=right_qpos){
		//															queryseq = "";
		//															if(left_aln_seg->bam and left_aln_seg->leftHardClippedFlag==false and left_aln_seg->rightHardClippedFlag==false)
		//																seq_int = bam_get_seq(left_aln_seg->bam);
		//															else if(right_aln_seg->bam and right_aln_seg->leftHardClippedFlag==false and right_aln_seg->rightHardClippedFlag==false)
		//																seq_int = bam_get_seq(right_aln_seg->bam);
		//															for(k=left_qpos-1; k<right_qpos; k++) queryseq += "=ACMGRSVTWYHKDBN"[bam_seqi(seq_int, k)];  // seq

																	queryseq = "";
																	if(left_aln_seg->aln_orient==query_orient){
																		if(left_qpos<=right_qpos) for(j=left_qpos-1; j<(size_t)right_qpos; j++) queryseq += "=ACMGRSVTWYHKDBN"[bam_seqi(seq_int_whole, j)];  // seq
																		else valid_flag = false;
																	}else{
																		start_qpos_tmp = query_len_whole - right_qpos;
																		end_qpos_tmp = query_len_whole - left_qpos;
																		if(start_qpos_tmp<=end_qpos_tmp){
																			for(j=start_qpos_tmp-1; j<(size_t)end_qpos_tmp; j++) queryseq += "=ACMGRSVTWYHKDBN"[bam_seqi(seq_int_whole, j)];  // seq
																			reverseComplement(queryseq);
																		}else valid_flag = false;
																	}

																	if(valid_flag){
																		rescue_seq_node->startRefPos_clip = left_rpos;
																		rescue_seq_node->endRefPos_clip = right_rpos;
																		rescue_seq_node->startQueryPos_clip = left_qpos;
																		rescue_seq_node->endQueryPos_clip = right_qpos;
																		rescue_seq_node->qseq_clip = queryseq;

																		reg_str_tmp = chrname + ":" + to_string(rescue_seq_node->startRefPos_clip) + "-" + to_string(rescue_seq_node->endRefPos_clip);
																		pthread_mutex_lock(&mutex_fai);
																		p_seq = fai_fetch(fai, reg_str_tmp.c_str(), &seq_len);
																		pthread_mutex_unlock(&mutex_fai);
																		rescue_seq_node->refseq_clip = p_seq;
																		free(p_seq);
																	}else {
																		delete rescue_seq_node;
																	}
																}
															}
															if(valid_flag) re_seq_vec.push_back(rescue_seq_node);
														}
													}
												}
											}

										}
									}
									for(j=0; j<query_aln_segs.size(); j++) query_aln_segs.at(j)->query_checked_flag = true;
								}
							}
						}
						if(re_seq_vec.size()>0) rescue_seqs_vec.push_back(re_seq_vec);
					}
				}
			}
		}

		if(rightClipRefPos-leftClipRefPos<=MAX_RESCUE_VAR_SIZE) {
			if(rescue_seqs_vec.size()>0){
				//rescue_readsfilename = out_dir_call + "/rescue_reads_" + chrname + "_" + to_string(leftClipRefPos) + "-" + to_string(rightClipRefPos) + ".fa";
				rescue_refseqfilename = out_dir_call + "/rescue_refseq_" + chrname + "_" + to_string(leftClipRefPos) + "-" + to_string(rightClipRefPos) + ".fa";
				rescue_cnsfilename = out_dir_call + "/rescue_cns_" + chrname + "_" + to_string(leftClipRefPos) + "-" + to_string(rightClipRefPos) + ".fa";
				//tmp_rescue_cnsfilename = out_dir_call + "/tmp_rescue_cns_" + chrname + "_" + to_string(leftClipRefPos) + "-" + to_string(rightClipRefPos) + ".fa";
				rescue_alnfilename = out_dir_call + "/rescue_minimap2_" + chrname + "_" + to_string(leftClipRefPos) + "-" + to_string(rightClipRefPos) + ".paf";

				outfile_rescue_refseq.open(rescue_refseqfilename);
				if(outfile_rescue_refseq.is_open()==false){
					cerr << "In " << __func__ << "(), cannot open file '" << rescue_refseqfilename << "', error!" << endl;
					exit(1);
				}
				outfile_rescue_cns.open(rescue_cnsfilename);
				if(outfile_rescue_cns.is_open()==false){
					cerr << "In " << __func__ << "(), cannot open file '" << rescue_cnsfilename << "', error!" << endl;
					exit(1);
				}

				// refseq file
				cons_header = ">" + reg_str + "___" + to_string(left_shift_size) + "___" + to_string(right_shift_size) + "___" + rescue_refseqfilename;
				outfile_rescue_refseq << cons_header << endl;
				outfile_rescue_refseq << refseq << endl;
				outfile_rescue_refseq.close();
			}

			serial_number = 1;
			for(cluster_id=0; cluster_id<rescue_seqs_vec.size(); cluster_id++){
				rescue_readsfilename = out_dir_call + "/rescue_reads_" + chrname + "_" + to_string(leftClipRefPos) + "-" + to_string(rightClipRefPos) + "_" + to_string(cluster_id) + ".fa";
				tmp_rescue_cnsfilename = out_dir_call + "/tmp_rescue_cns_" + chrname + "_" + to_string(leftClipRefPos) + "-" + to_string(rightClipRefPos) + "_" + to_string(cluster_id) + ".fa";

				re_seq_vec = rescue_seqs_vec.at(cluster_id);
				if(re_seq_vec.size()>0){
					outfile_rescue_reads.open(rescue_readsfilename);
					if(outfile_rescue_reads.is_open()==false){
						cerr << "In " << __func__ << "(), cannot open file '" << rescue_readsfilename << "', error!" << endl;
						exit(1);
					}

					// reads file
					for(i=0; i<re_seq_vec.size(); i++){
						rescue_seq_node = re_seq_vec.at(i);
			//			cout << "[" << i << "]: " << rescue_seq_node->chrname << ":" << rescue_seq_node->startRefPos << "-" << rescue_seq_node->endRefPos << endl;
			//			cout << rescue_seq_node->refseq << endl;
						outfile_rescue_reads << ">" << rescue_seq_node->qname << ":" << rescue_seq_node->startQueryPos << "-" << rescue_seq_node->endQueryPos << ":" << rescue_seq_node->aln_orient << endl;
						outfile_rescue_reads << rescue_seq_node->qseq << endl;
					}
					outfile_rescue_reads.close();

					// perform abpoa consensus
					abpoa_cmd = "abpoa -S -o " + tmp_rescue_cnsfilename + " " + rescue_readsfilename + " > /dev/null 2>&1";
					system(abpoa_cmd.c_str());

					flag = isFileExist(tmp_rescue_cnsfilename);
					if(flag){ // cons generated successfully
						// reorganize the consensus sequence
						FastaSeqLoader fa_loader(tmp_rescue_cnsfilename);
						for(i=0; i<fa_loader.getFastaSeqCount(); i++){
							cons_header = ">abpoa_rescueCns_";
							cons_header += to_string(serial_number) + "_"; //serial number
							cons_header += to_string(fa_loader.getFastaSeqLen(i)) + "_"; //consensus sequence length
							cons_header += to_string(re_seq_vec.size()) + "___";  //reads count

							for(j=0; j<re_seq_vec.size()-1; j++) cons_header += re_seq_vec.at(j)->qname + ";";
							cons_header += re_seq_vec.at(re_seq_vec.size()-1)->qname;

							outfile_rescue_cns << cons_header << endl; // header
							outfile_rescue_cns << fa_loader.getFastaSeq(i) << endl;  // seq

							serial_number ++;
						}
					}

					cmd_str = "rm -rf " + rescue_readsfilename;
					system(cmd_str.c_str());  // remove reads file
				}
			}
			outfile_rescue_cns.close();

			flag = isFileExist(rescue_cnsfilename);
			if(flag){
				minimap2_cmd = "minimap2 -c -x asm20 -o " + rescue_alnfilename + " " + rescue_refseqfilename + " " + rescue_cnsfilename + " > /dev/null 2>&1";
				status = system(minimap2_cmd.c_str());
				ret_status = getSuccessStatusSystemCmd(status);
				if(ret_status==0){ // command executed successfully
					// parse alignment information
					minimap2_aln_vec_rescue = minimap2Parse(rescue_alnfilename, rescue_cnsfilename, rescue_refseqfilename);

					// get reference sequence and query sequence
					rescue_var_vec = computeClipRegVarLoc(rescue_alnfilename, rescue_cnsfilename, rescue_refseqfilename, clusterfilename, minimap2_aln_vec_rescue, &clusterId_incomplete, true);
				}
			}
		}

		if(rescue_var_vec.size()==0){ // failed, then pick the first query sequence
			for(cluster_id=0; cluster_id<rescue_seqs_vec.size(); cluster_id++){
				re_seq_vec = rescue_seqs_vec.at(cluster_id);
				new_sv_type = -1;
				for(i=0; i<re_seq_vec.size(); i++){
					rescue_seq_node = re_seq_vec.at(i);
					//if(rescue_seq_node->startRefPos==leftClipRefPos and rescue_seq_node->endRefPos==rightClipRefPos){
					if(rescue_seq_node->qseq_clip.size()>0 and rescue_seq_node->refseq_clip.size()>0){

						sv_len = rescue_seq_node->qseq_clip.size() - rescue_seq_node->refseq_clip.size();
						//if(sv_len<0 and sv_type==VAR_DUP) new_sv_type = VAR_DEL;

						if(sv_type==VAR_DUP){ // added on 2025-08-20
							val = (double)abs(sv_len) / (rightClipRefPos - leftClipRefPos + 1);
							if(val>=VALID_DUP_INV_SIZE_RATIO_CLIP_REG){
								// allocate new node
								reg = new reg_t();
								if(new_sv_type!=-1) reg->var_type = new_sv_type;
								else reg->var_type = sv_type;
								reg->chrname = rescue_seq_node->chrname;
								reg->startRefPos = rescue_seq_node->startRefPos_clip;
								reg->endRefPos = rescue_seq_node->endRefPos_clip;
								reg->startLocalRefPos = -1;
								reg->endLocalRefPos = -1;
								reg->startQueryPos = rescue_seq_node->startQueryPos_clip;
								reg->endQueryPos = rescue_seq_node->endQueryPos_clip;
								reg->query_id = -1;
								reg->sv_len = sv_len;
								//if(sv_len<0) reg->sv_len = -reg->sv_len;
								reg->minimap2_aln_id = -1;
								reg->call_success_status = true;
								reg->short_sv_flag = false;
								reg->zero_cov_flag = false;
								reg->aln_seg_end_flag = false;
								reg->query_pos_invalid_flag = false;
								reg->merge_flag = false;
								reg->rescue_flag = true;
								reg->aln_orient = rescue_seq_node->aln_orient;
								reg->gt_type = -1;
								reg->gt_seq = "";
								reg->supp_num = re_seq_vec.size();
								reg->DP = reg->AF = 0;
								reg->discover_level = VAR_DISCOV_L_READS;
								reg->refseq =  rescue_seq_node->refseq_clip;
								reg->altseq = rescue_seq_node->qseq_clip;

								if(new_sv_type==-1){
									switch(sv_type){
										case VAR_DUP: // DUP
											reg->dup_num = dup_num;
											reg->large_indel_flag = false;
											break;
										case VAR_INV: // INV
											reg->dup_num = 0;
											reg->large_indel_flag = false;
											break;
										case VAR_INS: // large INS or DEL
										case VAR_DEL:
											reg->dup_num = 0;
											reg->large_indel_flag = true;
											break;
									}
								}

								rescue_var_vec.push_back(reg); // pick the first one
								break;
							}
						}
					}
				}
			}
		}

		// release memory
		for(cluster_id=0; cluster_id<rescue_seqs_vec.size(); cluster_id++){
			for(i=0; i<re_seq_vec.size(); i++){
				rescue_seq_node = re_seq_vec.at(i);
				destroyAlnSegs(rescue_seq_node->alnSegs_left);
				destroyAlnSegs(rescue_seq_node->alnSegs_right);
				delete rescue_seq_node;
			}
			vector<struct rescueSeqNode*>().swap(re_seq_vec);
		}
		vector<vector<struct rescueSeqNode*>>().swap(rescue_seqs_vec);

		if(!clipAlnDataVector.empty()) destoryClipAlnData(clipAlnDataVector);
	}

	return rescue_var_vec;
}

vector<int32_t> varCand::getMinimapItemIdVec(vector<string> &queryname_vec, vector<minimap2_aln_t*> &minimap2_aln_vec){
	vector<int32_t> minimap2_item_id_vec;
	minimap2_aln_t *minimap2_aln_item;
	size_t i, j, num;
	string queryname;

	num = 0;
	for(i=0; i<minimap2_aln_vec.size(); i++){
		minimap2_aln_item = minimap2_aln_vec.at(i);
		num = 0;
		for(j=0; j<minimap2_aln_item->qname.size(); j++){
			queryname = minimap2_aln_item->qname.at(j);
			if(find(queryname_vec.begin(), queryname_vec.end(), queryname)!=queryname_vec.end()){ // found
				num ++;
			}
		}
		if(num>=queryname_vec.size()*0.6) minimap2_item_id_vec.push_back(i);
	}

	return minimap2_item_id_vec;
}

vector<int32_t> varCand::getLeftRightMinimapItemId(int64_t leftClipRefPos, int64_t rightClipRefPos, vector<int32_t> &minimap2_item_id_vec, vector<minimap2_aln_t*> &minimap2_aln_vec){
	vector<int32_t> side_id_vec, sub_id_vec, query_id_vec;
	size_t i, k;
	minimap2_aln_t *minimap2_aln;
	struct pafalnSeg* paf_alnseg;
	int64_t start_aln_seg_pos, end_aln_seg_pos;
	int32_t id_vec, query_id, left_id, right_id; // left_orient, right_orient;
	int64_t max_left_query_pos, min_right_query_pos, leftClipRefPos_tmp, rightClipRefPos_tmp;

	leftClipRefPos_tmp = leftClipRefPos + 50;
	rightClipRefPos_tmp = rightClipRefPos - 50;
	if(rightClipRefPos_tmp<1) rightClipRefPos_tmp = 1;

	for(i=0; i<minimap2_item_id_vec.size(); i++){
		minimap2_aln = minimap2_aln_vec.at(minimap2_item_id_vec.at(i));
		if(find(query_id_vec.begin(), query_id_vec.end(), minimap2_aln->query_id)==query_id_vec.end()){ // not found, then add it
			query_id_vec.push_back(i);
		}
	}

	for(k=0; k<query_id_vec.size(); k++){
		query_id = query_id_vec.at(k);
		sub_id_vec.clear();
		for(i=0; i<minimap2_item_id_vec.size(); i++){
			minimap2_aln = minimap2_aln_vec.at(minimap2_item_id_vec.at(i));
			if(minimap2_aln->query_id==query_id) sub_id_vec.push_back(minimap2_item_id_vec.at(i));
		}

		left_id = right_id = -1;
		max_left_query_pos = INT_MIN;
		min_right_query_pos = INT_MAX;
		//left_orient = right_orient = -1;
		for(i=0; i<sub_id_vec.size(); i++){
			id_vec = sub_id_vec.at(i);
			minimap2_aln = minimap2_aln_vec.at(id_vec);
			if(minimap2_aln->pafalnsegs.size()>0 and minimap2_aln->type.compare("S")!=0){
				start_aln_seg_pos = minimap2_aln->pafalnsegs.at(0)->startRpos;
				paf_alnseg = minimap2_aln->pafalnsegs.at(minimap2_aln->pafalnsegs.size()-1);
				end_aln_seg_pos = getEndRefPosAlnSeg(paf_alnseg->startRpos, paf_alnseg->opflag, paf_alnseg->seglen);
				if(start_aln_seg_pos+1<=leftClipRefPos_tmp-VAR_ALN_EXTEND_SIZE and end_aln_seg_pos>=rightClipRefPos_tmp+VAR_ALN_EXTEND_SIZE){ // entirely flanking
					left_id = right_id = id_vec;
					//left_orient = right_orient = minimap2_aln->relative_strand;
					break;
				}else{
					if(start_aln_seg_pos+1<=leftClipRefPos_tmp){
						if(max_left_query_pos<minimap2_aln->query_start){ // pick the larger
//							if(right_orient==-1 or minimap2_aln->relative_strand==right_orient){
								max_left_query_pos = minimap2_aln->query_start;
								left_id = id_vec;
								//left_orient = minimap2_aln->relative_strand;
//							}
						}
					}

					if(end_aln_seg_pos+1>=rightClipRefPos_tmp){
						if(min_right_query_pos>minimap2_aln->query_end){ // pick the smaller
//							if(left_orient==-1 or minimap2_aln->relative_strand==left_orient){
								min_right_query_pos = minimap2_aln->query_end;
								right_id = id_vec;
								//right_orient = minimap2_aln->relative_strand;
//							}
						}
					}
				}
			}
		}
		if(left_id!=-1 and right_id!=-1){
			side_id_vec.push_back(query_id);
			side_id_vec.push_back(left_id);
			side_id_vec.push_back(right_id);
			break;
		}
	}

	return side_id_vec;
}

vector<int32_t> varCand::getLeftRightClipAlnId(int64_t leftClipRefPos, int64_t rightClipRefPos, vector<clipAlnData_t*> &query_aln_segs){
	vector<int32_t> side_id_vec;
	size_t i;
	clipAlnData_t *clip_aln;
	int64_t start_aln_seg_pos, end_aln_seg_pos;
	int32_t left_id, right_id, query_start, query_end, extend_size; // left_orient, right_orient;
	int64_t max_left_query_pos, min_right_query_pos, leftClipRefPos_tmp, rightClipRefPos_tmp;

	//extend_size = min<int64_t>(EXT_SIZE_CHK_VAR_LOC, (rightClipRefPos-leftClipRefPos+1)/4); // deleted on 2026-02-01
	extend_size = min<int64_t>(EXT_SIZE_CHK_VAR_LOC*2, (rightClipRefPos-leftClipRefPos+1)/4);

	leftClipRefPos_tmp = leftClipRefPos + extend_size;
	rightClipRefPos_tmp = rightClipRefPos - extend_size;
	if(rightClipRefPos_tmp<1) rightClipRefPos_tmp = 1;

	left_id = right_id = -1;
	max_left_query_pos = INT_MIN;
	min_right_query_pos = INT_MAX;
	//left_orient = right_orient = -1;
	for(i=0; i<query_aln_segs.size(); i++){
		clip_aln = query_aln_segs.at(i);

		if(clip_aln->bam and (clip_aln->ref_dist>=EXT_SIZE_CHK_VAR_LOC_SMALL or clip_aln->query_dist>=EXT_SIZE_CHK_VAR_LOC_SMALL)){ // added on 2025-09-16
			start_aln_seg_pos = clip_aln->startRefPos;
			end_aln_seg_pos = clip_aln->endRefPos;
			if(start_aln_seg_pos<=leftClipRefPos_tmp-VAR_ALN_EXTEND_SIZE and end_aln_seg_pos>=rightClipRefPos_tmp+VAR_ALN_EXTEND_SIZE){
				left_id = right_id = i;
				//left_orient = right_orient = clip_aln->aln_orient;
				break;
			}else{
				if(clip_aln->aln_orient==ALN_PLUS_ORIENT){
					query_start = clip_aln->startQueryPos;
					query_end = clip_aln->endQueryPos;
				}else{
//					query_start = clip_aln->querylen - clip_aln->startQueryPos + 1;
//					query_end = clip_aln->querylen - clip_aln->endQueryPos + 1;
					query_start = clip_aln->endQueryPos;
					query_end = clip_aln->startQueryPos;
				}

				if(start_aln_seg_pos<=leftClipRefPos_tmp){
					if(max_left_query_pos<query_start){ // pick the larger
//						if(right_orient==-1 or clip_aln->aln_orient==right_orient){
							if(left_id!=-1 and right_id==-1){
								min_right_query_pos = max_left_query_pos;
								right_id = left_id;
//
							}
							max_left_query_pos = query_start;
							left_id = i;
							//left_orient = clip_aln->aln_orient;
//						}
					}
				}else if(end_aln_seg_pos>=rightClipRefPos_tmp){
					if(min_right_query_pos>query_end){ // pick the smaller
//						if(left_orient==-1 or clip_aln->aln_orient==left_orient){
							if(right_id!=-1 and left_id==-1){
								max_left_query_pos = min_right_query_pos;
								left_id = right_id;
							}
							min_right_query_pos = query_end;
							right_id = i;
							//right_orient = clip_aln->aln_orient;

//						}
					}
				}
			}
		}
	}

	if(left_id!=-1 and right_id!=-1){
		if(query_aln_segs.at(left_id)->startRefPos<=query_aln_segs.at(right_id)->startRefPos){
			side_id_vec.push_back(left_id);
			side_id_vec.push_back(right_id);
		}else{
			side_id_vec.push_back(right_id);
			side_id_vec.push_back(left_id);
		}
	}

	return side_id_vec;
}

vector<int32_t> varCand::getLeftRightClipAlnIdINV(int64_t leftClipRefPos, int64_t rightClipRefPos, vector<clipAlnData_t*> &query_aln_segs){
	vector<int32_t> side_id_vec;
	size_t i;
	clipAlnData_t *clip_aln;
	int64_t start_aln_seg_pos, end_aln_seg_pos;
	int32_t left_id, right_id, query_start, query_end, extend_size; // left_orient, right_orient;
	int64_t max_left_query_pos, min_right_query_pos, leftClipRefPos_tmp, rightClipRefPos_tmp;
	bool overlap_flag;

	//extend_size = min<int64_t>(EXT_SIZE_CHK_VAR_LOC, (rightClipRefPos-leftClipRefPos+1)/4); // deleted on 2026-02-01
	extend_size = min<int64_t>(EXT_SIZE_CHK_VAR_LOC*2, (rightClipRefPos-leftClipRefPos+1)/4);

	leftClipRefPos_tmp = leftClipRefPos + extend_size;
	rightClipRefPos_tmp = rightClipRefPos - extend_size;
	if(rightClipRefPos_tmp<1) rightClipRefPos_tmp = 1;

	left_id = right_id = -1;
	max_left_query_pos = INT_MIN;
	min_right_query_pos = INT_MAX;
	//left_orient = right_orient = -1;
	for(i=0; i<query_aln_segs.size(); i++){
		clip_aln = query_aln_segs.at(i);
		overlap_flag = isOverlappedPos(clip_aln->startRefPos, clip_aln->endRefPos, leftClipRefPos-extend_size, rightClipRefPos+extend_size);

		//if(clip_aln->bam and (clip_aln->ref_dist>=EXT_SIZE_CHK_VAR_LOC_SMALL or clip_aln->query_dist>=EXT_SIZE_CHK_VAR_LOC_SMALL)){ // added on 2025-09-16
		if(overlap_flag){
			start_aln_seg_pos = clip_aln->startRefPos;
			end_aln_seg_pos = clip_aln->endRefPos;
//			if(start_aln_seg_pos<=leftClipRefPos_tmp-VAR_ALN_EXTEND_SIZE and end_aln_seg_pos>=rightClipRefPos_tmp+VAR_ALN_EXTEND_SIZE){ // deleted on 2026-04-17
//				left_id = right_id = i;
//				//left_orient = right_orient = clip_aln->aln_orient;
//				break;
//			}else{
				if(clip_aln->aln_orient==ALN_PLUS_ORIENT){
					query_start = clip_aln->startQueryPos;
					query_end = clip_aln->endQueryPos;
				}else{
//					query_start = clip_aln->querylen - clip_aln->startQueryPos + 1;
//					query_end = clip_aln->querylen - clip_aln->endQueryPos + 1;
					query_start = clip_aln->endQueryPos;
					query_end = clip_aln->startQueryPos;
				}

				if(start_aln_seg_pos<=leftClipRefPos_tmp){
					if(max_left_query_pos<query_start){ // pick the larger
//						if(right_orient==-1 or clip_aln->aln_orient==right_orient){
							if(left_id!=-1 and right_id==-1){
								min_right_query_pos = max_left_query_pos;
								right_id = left_id;
							}
							max_left_query_pos = query_start;
							left_id = i;
							//left_orient = clip_aln->aln_orient;
//						}
					}
				}else if(end_aln_seg_pos>=rightClipRefPos_tmp){
					if(min_right_query_pos>query_end){ // pick the smaller
//						if(left_orient==-1 or clip_aln->aln_orient==left_orient){
							if(right_id!=-1 and left_id==-1){
								max_left_query_pos = min_right_query_pos;
								left_id = right_id;
							}
							min_right_query_pos = query_end;
							right_id = i;
							//right_orient = clip_aln->aln_orient;
//						}
					}
				}
//			}
		}
	}
	if(left_id!=-1 and right_id!=-1){
		if(query_aln_segs.at(left_id)->startRefPos<=query_aln_segs.at(right_id)->startRefPos){
			side_id_vec.push_back(left_id);
			side_id_vec.push_back(right_id);
		}else{
			side_id_vec.push_back(right_id);
			side_id_vec.push_back(left_id);
		}
	}

	return side_id_vec;
}

vector<reg_t*> varCand::callIndelFromReadsIndelReg(vector<int32_t> &clusterId_incomplete){
	vector<reg_t*> var_vec;
	size_t i, j, k, cluster_id;
	reg_t *reg, *reg_new;
	int32_t seq_len, max_merge_span, DP_num;
	vector<vector<string>> querynames_vec;
	vector<clipAlnData_t*> clipAlnDataVector, clipAlnDataVector_sub;
	int64_t end_ref_pos, start_var_pos, end_var_pos, start_var_pos_pure0, end_var_pos_pure0, start_var_pos_pure, end_var_pos_pure, chrlen_tmp;
	vector<struct querySeqInfoNode*> query_seq_info_all;
	struct querySeqInfoNode *qnode;
	vector<struct alnSeg*> query_alnSegs;
	struct alnSeg *seg;
	vector<qcSig_t*> overlap_qcSig_vec;
	qcSig_t *qcSig;
	string queryname, reg_str, refseq, refbase;
	vector<string> qname_vec;
	char *p_seq;
	bool overlap_flag;
	map<int32_t, vector<int32_t> > occ_map;

	int32_t ins_sum_num, del_sum_num, target_var_op, item_num, maxVal, maxIdx;
	double len_ratio, ins_sum_mean0, del_sum_mean0, ins_sum_mean, del_sum_mean, mean_size, sdev, dif;

	chrlen_tmp = faidx_seq_len64(fai, chrname.c_str());
	start_var_pos_pure0 = varVec.at(0)->startRefPos;
	end_var_pos_pure0 = varVec.at(varVec.size()-1)->endRefPos;
	start_var_pos_pure = start_var_pos_pure0 - SHORT_VAR_ALN_CHECK_EXTEND_SIZE;
	end_var_pos_pure = end_var_pos_pure0 + SHORT_VAR_ALN_CHECK_EXTEND_SIZE;
	if(start_var_pos_pure<1) start_var_pos_pure = 1;
	if(end_var_pos_pure>chrlen_tmp) end_var_pos_pure = chrlen_tmp;

	start_var_pos = start_var_pos_pure - EXT_SIZE_CHK_VAR_LOC;
	end_var_pos = end_var_pos_pure + EXT_SIZE_CHK_VAR_LOC;
	if(start_var_pos<1) start_var_pos = 1;
	if(end_var_pos>chrlen_tmp) end_var_pos = chrlen_tmp;

	max_merge_span = min_distance_merge; // CNS_MIN_DIST_MERGE_THRES;
//	if(sv_len_est>CNS_MIN_DIST_MERGE_THRES) max_merge_span = 1.5 * CNS_MIN_DIST_MERGE_THRES;

	reg_str = chrname + ":" + to_string(start_var_pos) + "-" + to_string(end_var_pos);
	pthread_mutex_lock(&mutex_fai);
	p_seq = fai_fetch(fai, reg_str.c_str(), &seq_len);
	pthread_mutex_unlock(&mutex_fai);
	refseq = p_seq;
	free(p_seq);

	// get cluster information
	querynames_vec = getClusterInfo(clusterfilename);

	// load the clipping data
	clipAlnDataLoader data_loader(chrname, start_var_pos, end_var_pos, inBamFile, minClipEndSize, maxVarRegSize, minMapQ, minHighMapQ);
	if(clip_reg_flag) data_loader.loadClipAlnDataWithSATag(clipAlnDataVector, max_ultra_high_cov);
	else data_loader.loadClipAlnDataWithSATagWithSegSize(clipAlnDataVector, max_ultra_high_cov, max_seg_size_ratio, max_seg_nm_ratio, fai, max_absig_density);

	for(cluster_id=0; cluster_id<querynames_vec.size(); cluster_id++){
		if(find(clusterId_incomplete.begin(), clusterId_incomplete.end(), cluster_id) == clusterId_incomplete.end()){  // call complete, then skip
			continue;
		}

		qname_vec = querynames_vec.at(cluster_id);

		// get the subset
		clipAlnDataVector_sub.clear();
		for(i=0; i<clipAlnDataVector.size(); i++){
			queryname = clipAlnDataVector.at(i)->queryname;
			if(find(qname_vec.begin(), qname_vec.end(), queryname) != qname_vec.end()) clipAlnDataVector_sub.push_back(clipAlnDataVector.at(i));
		}

		// extract queries from clip align data vector
		query_seq_info_all = extractQueriesFromClipAlnDataVec(clipAlnDataVector_sub, refseq, chrname, start_var_pos, end_var_pos, start_var_pos_pure, end_var_pos_pure, fai, minConReadLen, clip_reg_flag, &mutex_fai);

		for(i=0; i<query_seq_info_all.size(); i++){
			qnode = query_seq_info_all.at(i);
			query_alnSegs = qnode->query_alnSegs;
			if(query_alnSegs.size()>0){
				seg = query_alnSegs.at(query_alnSegs.size()-1);
				end_ref_pos = getEndRefPosAlnSeg(seg->startRpos, seg->opflag, seg->seglen);
				if(query_alnSegs.at(0)->startRpos<=start_var_pos_pure0+SHORT_VAR_ALN_CHECK_EXTEND_SIZE and end_ref_pos>=end_var_pos_pure0-SHORT_VAR_ALN_CHECK_EXTEND_SIZE) qnode->entire_flanking_flag = true;
				else qnode->entire_flanking_flag = false;

				qnode->qcSig_vec = extractQcSigsFromAlnSegsSingleQuery(qnode, qnode->clip_aln->chrname, start_var_pos, end_var_pos, min_sv_size*2);
				// remove neighboring false signatures
				rmNeighborFalseSig(qnode->qcSig_vec, MIN_AVER_SIZE_ALN_SEG, MIN_SIZE_RATIO_MATCH_CLIP_POS);
				qnode->qcSig_merge_flag = mergeNeighbouringSigsFlag(qnode, qnode->qcSig_vec, MAX_DIST_MERGE_ARBITARY, max_merge_span, CALL_MIN_SEQSIM_MERGE_THRES2, MIN_VALID_SIG_SIZE_RATIO_THRES, fai);
				filterSmallSigs(qnode->qcSig_vec, MIN_SIG_SIZE_RATIO_FILTER_THRES, MAX_INNER_MISSING_IGNORE_SIZE);
			}
		}

#if READS_CALL_DEBUG
		printQcSigIndel(query_seq_info_all);
#endif

		for(i=0; i<varVec.size(); i++){
			reg = varVec.at(i);
			if(reg->sv_len!=0){
				overlap_qcSig_vec.clear();
				occ_map.clear();
				DP_num = 0;
				for(j=0; j<query_seq_info_all.size(); j++){
					qnode = query_seq_info_all.at(j);
					for(k=0; k<qnode->qcSig_vec.size(); k++){
						qcSig = qnode->qcSig_vec.at(k);
						if(qcSig->cigar_op==reg->var_type){ // same variant type
							if(abs(reg->sv_len)>qcSig->cigar_op_len) len_ratio = (double)qcSig->cigar_op_len / abs(reg->sv_len);
							else len_ratio = (double)abs(reg->sv_len) / qcSig->cigar_op_len;
							if(len_ratio>=QC_SIZE_RATIO_MATCH_THRES_INDEL){
								end_ref_pos = getEndRefPosAlnSeg(qcSig->ref_pos, qcSig->cigar_op, qcSig->cigar_op_len);
								overlap_flag = isOverlappedPos(reg->startRefPos-MAX_DIST_MERGE_ARBITARY, reg->endRefPos+MAX_DIST_MERGE_ARBITARY, qcSig->ref_pos-MAX_DIST_MERGE_ARBITARY, end_ref_pos+MAX_DIST_MERGE_ARBITARY);
								if(overlap_flag){
									overlap_qcSig_vec.push_back(qcSig);
								}
							}
						}
					}
					if(qnode->entire_flanking_flag==true) DP_num ++;
				}
#if READS_CALL_DEBUG
				cout << "overlap_qcSig_vec.size=" << overlap_qcSig_vec.size() << endl;
				for(j=0; j<overlap_qcSig_vec.size(); j++){
					qcSig = overlap_qcSig_vec.at(j);
					cout << qcSig->cigar_op << "," << qcSig->chrname << "," << qcSig->ref_pos << "," << qcSig->cigar_op_len << ", altseq=" << qcSig->altseq << ", refseq=" << qcSig->refseq << endl;
				}
#endif

				// compute
				if(overlap_qcSig_vec.size()>=1.5*minReadsNumSupportSV){
					// compute mean size
					ins_sum_mean0 = del_sum_mean0 = ins_sum_num = del_sum_num = 0;
					for(j=0; j<overlap_qcSig_vec.size(); j++){
						qcSig = overlap_qcSig_vec.at(j);
						if(qcSig->cigar_op==BAM_CINS){
							ins_sum_mean0 += qcSig->cigar_op_len;
							ins_sum_num ++;
						}else if(qcSig->cigar_op==BAM_CDEL){
							del_sum_mean0 += qcSig->cigar_op_len;
							del_sum_num ++;
						}
					}
					if(ins_sum_num>0) ins_sum_mean0 /= ins_sum_num;
					if(del_sum_num>0) del_sum_mean0 /= del_sum_num;

					// skip outliers
					ins_sum_mean = del_sum_mean = ins_sum_num = del_sum_num = 0;
					for(j=0; j<overlap_qcSig_vec.size(); j++){
						qcSig = overlap_qcSig_vec.at(j);
						if(qcSig->cigar_op==BAM_CINS){
							if(qcSig->cigar_op_len<ins_sum_mean0) len_ratio = (double)qcSig->cigar_op_len / ins_sum_mean0;
							else len_ratio = (double)ins_sum_mean0 / qcSig->cigar_op_len;
							if(len_ratio>=QC_SIZE_RATIO_MATCH_THRES_INDEL){
								ins_sum_mean += qcSig->cigar_op_len;
								ins_sum_num ++;
							}
						}else if(qcSig->cigar_op==BAM_CDEL){
							if(qcSig->cigar_op_len<del_sum_mean0) len_ratio = (double)qcSig->cigar_op_len / del_sum_mean0;
							else len_ratio = (double)del_sum_mean0 / qcSig->cigar_op_len;
							if(len_ratio>=QC_SIZE_RATIO_MATCH_THRES_INDEL){
								del_sum_mean += qcSig->cigar_op_len;
								del_sum_num ++;
							}
						}
					}
					if(ins_sum_num>0) ins_sum_mean /= ins_sum_num;
					if(del_sum_num>0) del_sum_mean /= del_sum_num;

					if(ins_sum_mean0>del_sum_mean) { target_var_op = BAM_CINS; mean_size = ins_sum_mean; item_num = ins_sum_num; }
					else { target_var_op = BAM_CDEL; mean_size = del_sum_mean; item_num = del_sum_num; }
					if(item_num>DP_num) DP_num = item_num;

					// compute sdev
					if(mean_size>0){
						item_num = sdev = 0;
						for(j=0; j<overlap_qcSig_vec.size(); j++){
							qcSig = overlap_qcSig_vec.at(j);
							if(target_var_op==BAM_CINS and qcSig->cigar_op==BAM_CINS){
								if(qcSig->cigar_op_len<ins_sum_mean0) len_ratio = (double)qcSig->cigar_op_len / ins_sum_mean0;
								else len_ratio = (double)ins_sum_mean0 / qcSig->cigar_op_len;
								if(len_ratio>=QC_SIZE_RATIO_MATCH_THRES_INDEL){
									sdev += (qcSig->cigar_op_len - mean_size) * (qcSig->cigar_op_len - mean_size);
									item_num ++;
								}
							}else if(target_var_op==BAM_CDEL and qcSig->cigar_op==BAM_CDEL){
								if(qcSig->cigar_op_len<del_sum_mean0) len_ratio = (double)qcSig->cigar_op_len / del_sum_mean0;
								else len_ratio = (double)del_sum_mean0 / qcSig->cigar_op_len;
								if(len_ratio>=QC_SIZE_RATIO_MATCH_THRES_INDEL){
									sdev += (qcSig->cigar_op_len - mean_size) * (qcSig->cigar_op_len - mean_size);
									item_num ++;
								}
							}
						}
						sdev = sqrt(sdev/item_num);

						// select the minimum difference
						if(item_num>0){
							item_num = 0;
							for(j=0; j<overlap_qcSig_vec.size(); j++){
								qcSig = overlap_qcSig_vec.at(j);
								if(target_var_op==qcSig->cigar_op) {
									dif = abs(qcSig->cigar_op_len - mean_size);
									if(dif<2*sdev){
										item_num ++;
										occ_map[qcSig->cigar_op_len].push_back(j);
									}
								}
							}

#if READS_CALL_DEBUG
							cout << "target_var_op=" << target_var_op << ", mean_size=" << mean_size << ", sdev=" << sdev << ", item_num=" << item_num << ", DP=" << DP_num << endl;
#endif

							maxVal = maxIdx = -1;
							for(auto const& pair : occ_map){
								const vector<int32_t>& idx_vec = pair.second;
								if(idx_vec.size()>0 and (int32_t)idx_vec.size()>maxVal){
									maxVal = idx_vec.size();
									maxIdx = idx_vec.at(0);
								}
							}
							if(maxIdx!=-1 and maxVal>=2 and (double)item_num/DP_num>=MIN_RATIO_ALLELE_CALL){ // found maximum occurrence
								qcSig = overlap_qcSig_vec.at(maxIdx);
								if(qcSig->ref_pos>start_var_pos) refbase = refseq.at(qcSig->ref_pos-start_var_pos-1);
								else refbase = "";

								//write to newVarVec
								reg_new = new reg_t();
								reg_new->var_type = target_var_op;
								reg_new->chrname = qcSig->chrname;
								reg_new->startRefPos = qcSig->ref_pos - 1;
								if(reg_new->var_type==VAR_INS) reg_new->endRefPos = reg_new->startRefPos;
								else reg_new->endRefPos = getEndRefPosAlnSeg(qcSig->ref_pos, qcSig->cigar_op, qcSig->cigar_op_len);
								reg_new->startLocalRefPos = -1;
								reg_new->endLocalRefPos = -1;
								reg_new->startQueryPos = -1;
								reg_new->endQueryPos = -1;
								reg_new->query_id = -1;
								reg_new->sv_len = qcSig->cigar_op_len;
								if(reg_new->var_type==VAR_DEL) reg_new->sv_len = -reg_new->sv_len;
								reg_new->minimap2_aln_id = -1;
								reg_new->call_success_status = true;
								reg_new->short_sv_flag = false;
								reg_new->zero_cov_flag = false;
								reg_new->aln_seg_end_flag = false;
								reg_new->query_pos_invalid_flag = false;
								reg_new->large_indel_flag = false;
								reg_new->merge_flag = false;
								reg_new->rescue_flag = false;
								reg_new->aln_orient = ALN_PLUS_ORIENT;
								reg_new->gt_type = -1;
								reg_new->gt_seq = "";
								reg_new->refseq =  refbase + qcSig->refseq;
								reg_new->altseq = refbase + qcSig->altseq;
								reg_new->AF = 0;
								reg_new->supp_num = item_num;
								reg_new->DP = DP_num;
								reg_new->discover_level = VAR_DISCOV_L_READS2;

								var_vec.push_back(reg_new);
							}
						}
					}
				}
			}
		}
		if(!query_seq_info_all.empty()) destoryQuerySeqInfoAll(query_seq_info_all);
	}
	// release memory
	if(!clipAlnDataVector.empty()) destoryClipAlnData(clipAlnDataVector);

#if READS_CALL_DEBUG
	cout << __func__ << ": " << alnfilename << ", var_vec.size=" << var_vec.size() << endl;
#endif

	return var_vec;
}

vector<reg_t*> varCand::callVarFromReadsClipReg(vector<int32_t> &clusterId_incomplete){
	vector<reg_t*> var_vec, var_vec_inner, var_vec_sub;
	size_t i, j, cluster_id;
	int32_t sv_len, seq_len, DP_num, sum, supp_num, supp_num_tmp, DP_num_tmp, ins_type_num, del_type_num, dup_type_num, inv_type_num, invdup_type_num, tra_type_num, maxIdx, overlap_size;
	vector<clipAlnData_t*> clipAlnDataVector, clipAlnDataVector_sub, query_aln_segs;
	int64_t start_var_pos, end_var_pos, start_var_pos_pure0, end_var_pos_pure0, start_var_pos_pure, end_var_pos_pure, chrlen_tmp, dist_ref, dist_query;
	string queryname, reg_str, refseq, queryseq, tmp_refseq, tmp_queryseq;
	vector<vector<string>> querynames_vec;
	vector<string> str_vec, str_vec2, str_vec2_tmp, qname_vec, qname_dup_vec, qname_inv_vec, qname_invdup_vec, qname_tra_vec;
	char *p_seq;
	clipAlnData_t *clip_aln, *clip_aln_left, *clip_aln_right;
	bool same_chr_flag, same_orient_flag, self_overlap_flag, same_aln_reg_flag, both_ends_overlap, turn_round_flag, valid_query_flag, entire_flanking_flag; // six features
	int64_t start_var_ref_pos, end_var_ref_pos, start_var_query_pos, end_var_query_pos, start_qpos_tmp, end_qpos_tmp;
	reg_t *reg_new, *reg_tmp;
	vector<struct querySeqInfoNode*> query_seq_info_all;
	struct querySeqInfoNode *qnode_left, *qnode_right;
	vector<qcSigList_t*> qcSigList_vec;
	qcSigList_t *qcSigList_node;
	double val, val2, maxVal, mean_left_clip_rpos, mean_right_clip_rpos;
	vector<int32_t> id_vec, left_clip_pos_vec, right_clip_pos_vec;
	uint8_t *seq_int;
	int32_t idx_alnSegs_left, idx_alnSegs_right, margin_dist1, margin_dist2, left_clip_qpos, right_clip_qpos, clip_end1, clip_end2, noHardClipIdx;
	int64_t left_rpos, left_qpos, right_rpos, right_qpos;
	struct alnSeg *aln_seg;
	int64_t left_clip_rpos, right_clip_rpos, left_clip_rpos_new, right_clip_rpos_new, pos_tmp, ins_sum_mean, del_sum_mean;
	vector<qcSig_t*> qc_sig_vec;
	qcSig_t *qc_sig;
	int8_t var_flag_arr[varVec.size()];
	int32_t min_dif_idx;
	double sum_len, mean_size, dif, min_dif;

	chrlen_tmp = faidx_seq_len64(fai, chrname.c_str());
	start_var_pos_pure0 = leftClipRefPos;
	end_var_pos_pure0 = rightClipRefPos;
	start_var_pos_pure = start_var_pos_pure0 - SHORT_VAR_ALN_CHECK_EXTEND_SIZE;
	end_var_pos_pure = end_var_pos_pure0 + SHORT_VAR_ALN_CHECK_EXTEND_SIZE;
	if(start_var_pos_pure<1) start_var_pos_pure = 1;
	if(end_var_pos_pure>chrlen_tmp) end_var_pos_pure = chrlen_tmp;

	start_var_pos = start_var_pos_pure - EXT_SIZE_CHK_VAR_LOC;
	end_var_pos = end_var_pos_pure + EXT_SIZE_CHK_VAR_LOC;
	if(start_var_pos<1) start_var_pos = 1;
	if(end_var_pos>chrlen_tmp) end_var_pos = chrlen_tmp;

	// load the clipping data
	clipAlnDataLoader data_loader(chrname, start_var_pos, end_var_pos, inBamFile, minClipEndSize, maxVarRegSize, minMapQ, minHighMapQ);
	if(clip_reg_flag) data_loader.loadClipAlnDataWithSATag(clipAlnDataVector, max_ultra_high_cov);
	else data_loader.loadClipAlnDataWithSATagWithSegSize(clipAlnDataVector, max_ultra_high_cov, max_seg_size_ratio, max_seg_nm_ratio, fai, max_absig_density);

	// get cluster information
	querynames_vec = getClusterInfo(clusterfilename);

	reg_str = chrname + ":" + to_string(start_var_pos) + "-" + to_string(end_var_pos);
	pthread_mutex_lock(&mutex_fai);
	p_seq = fai_fetch(fai, reg_str.c_str(), &seq_len);
	pthread_mutex_unlock(&mutex_fai);
	refseq = p_seq;
	free(p_seq);

	for(i=0; i<varVec.size(); i++) var_flag_arr[i] = 0;
	for(i=0; i<clipAlnDataVector.size(); i++) clipAlnDataVector.at(i)->query_checked_flag = false;

	for(cluster_id=0; cluster_id<querynames_vec.size(); cluster_id++){
		if(find(clusterId_incomplete.begin(), clusterId_incomplete.end(), cluster_id) == clusterId_incomplete.end()){  // call complete, then skip
			continue;
		}

		qname_vec = querynames_vec.at(cluster_id);

		// get the target group of align data
		clipAlnDataVector_sub.clear();
		for(auto const& item : clipAlnDataVector) if(find(qname_vec.begin(), qname_vec.end(), item->queryname)!=qname_vec.end()) clipAlnDataVector_sub.push_back(item);

		if(sv_type==VAR_INS or sv_type==VAR_DUP or sv_type==VAR_INV or sv_type==VAR_TRA or sv_type==VAR_BND){ // INS, DUP, INV, TRA/BND
			dup_type_num = inv_type_num = invdup_type_num = tra_type_num = 0;
			for(auto const& clip_aln : clipAlnDataVector_sub){
				if(clip_aln->query_checked_flag==false){
					queryname = clip_aln->queryname;

//					if(queryname.compare("m84039_230414_235240_s2/173539883/ccs")==0){
//						cout << "line=" << __LINE__ << ", queryname=" << queryname << endl;
//					}

					if(find(qname_vec.begin(), qname_vec.end(), queryname) != qname_vec.end()){ // found

						query_aln_segs = getQueryClipAlnSegsAll(queryname, clipAlnDataVector_sub);
						if(query_aln_segs.size()==0) continue;
						valid_query_flag = Filteredbychrname(query_aln_segs);
						if(valid_query_flag==false){
							for(j=0; j<query_aln_segs.size(); j++) query_aln_segs.at(j)->query_checked_flag = true;
							continue;
						}

						//sort according to queryPos in reference position
						sortQueryAlnSegs(query_aln_segs);
						FilteredbyAlignmentSegment(query_aln_segs, chrname, minMapQ, MIN_ALN_SEG_SIZE_CLIPREG);
						if(query_aln_segs.size()==0) continue;

						same_chr_flag = isSameChrome(query_aln_segs);
						//same_orient_flag = isSameOrient(query_aln_segs, MIN_ALN_SIZE_DIFF_ORIENT);
						same_orient_flag = isSameOrient(query_aln_segs, MIN_ALN_SEG_SIZE_CLIPREG);
						self_overlap_flag = isQuerySelfOverlap(query_aln_segs, maxVarRegSize);
						same_aln_reg_flag = isSameAlnReg(query_aln_segs, maxVarRegSize);
						both_ends_overlap = isBothEndsOverlap(query_aln_segs, leftClipRefPos, rightClipRefPos);
						turn_round_flag = isTurnRound(query_aln_segs, leftClipRefPos, rightClipRefPos);

						// predict the variant type
						//if(same_chr_flag and same_orient_flag and self_overlap_flag and same_aln_reg_flag and both_ends_overlap){ // DUP, deleted on 2026-04-09
						if(same_chr_flag and same_orient_flag and both_ends_overlap and turn_round_flag){ // DUP
#if READS_CALL_DEBUG
							cout << "DUP: " << queryname << ", aln_size=" << query_aln_segs.size() << endl;
#endif
							qname_dup_vec.push_back(queryname);
							dup_type_num ++;
						}else if(same_chr_flag and same_orient_flag==false and both_ends_overlap){// and same_aln_reg_flag){ // INV
//							if(self_overlap_flag==false){
#if READS_CALL_DEBUG
								cout << "INV: " << queryname << ", aln_size=" << query_aln_segs.size() << endl;
#endif
								qname_inv_vec.push_back(queryname);
								inv_type_num ++;
//							}else{
//#if READS_CALL_DEBUG
//								cout << "INVDUP: " << queryname << ", aln_size=" << query_aln_segs.size() << endl;
//#endif
//								qname_invdup_vec.push_back(queryname);
//								invdup_type_num ++;
//							}
						//}else if(self_overlap_flag==false and same_aln_reg_flag==false)	{ // TRA
						}else if(self_overlap_flag==false and same_aln_reg_flag==false and both_ends_overlap)	{ // TRA
#if READS_CALL_DEBUG
							cout << "TRA: " << queryname << ", aln_size=" << query_aln_segs.size() << endl;
#endif
							qname_tra_vec.push_back(queryname);
							tra_type_num ++;
						}
						for(j=0; j<query_aln_segs.size(); j++) query_aln_segs.at(j)->query_checked_flag = true;
					}
				}
			}
			for(auto const& clip_aln : clipAlnDataVector_sub) clip_aln->query_checked_flag = false;

#if READS_CALL_DEBUG
			cout << "sv_type=" << sv_type << ", dup_type_num=" << dup_type_num << ", inv_type_num=" << inv_type_num << ", invdup_type_num=" << invdup_type_num << ", tra_type_num=" << tra_type_num << endl;
#endif

			supp_num = DP_num = 0;
			for(i=0; i<4; i++){
				if(bnd_mate_reg_strs[i].compare("-")!=0){
					str_vec = split(bnd_mate_reg_strs[i], ",");
					for(j=0; j<str_vec.size(); j++){
						if(str_vec.at(j).compare("-")!=0){
							str_vec2_tmp = split(str_vec.at(j), "_");
							for(const auto& item : str_vec2_tmp){
								str_vec2 = split(item, "|");
								DP_num_tmp = stoi(str_vec2.at(3));
								if(DP_num_tmp>DP_num) DP_num = DP_num_tmp;
							}
						}
					}
				}
			}

			switch(sv_type){
				case VAR_INS:
				case VAR_DUP:
					supp_num = dup_type_num;
					if(supp_num>DP_num) DP_num = supp_num;
					//if(supp_num>=minReadsNumSupportSV){ // deleted on 2026-04-13
					if(supp_num>=minReadsNumSupportSV*READS_NUM_SUPPORT_FACTOR){
						reg_new = new reg_t();
						//reg_new->var_type = sv_type;
						reg_new->var_type = VAR_DUP;
						reg_new->chrname = chrname;
						reg_new->startRefPos = leftClipRefPos - 1;
						reg_new->endRefPos = rightClipRefPos;
						reg_new->startLocalRefPos = -1;
						reg_new->endLocalRefPos = -1;
						reg_new->startQueryPos = -1;
						reg_new->endQueryPos = -1;
						reg_new->query_id = -1;
						reg_new->sv_len = rightClipRefPos - leftClipRefPos + 1;
						reg_new->minimap2_aln_id = -1;
						reg_new->call_success_status = true;
						reg_new->short_sv_flag = false;
						reg_new->zero_cov_flag = false;
						reg_new->aln_seg_end_flag = false;
						reg_new->query_pos_invalid_flag = false;
						reg_new->large_indel_flag = false;
						reg_new->merge_flag = false;
						reg_new->rescue_flag = true;
						reg_new->aln_orient = ALN_PLUS_ORIENT;
						reg_new->gt_type = -1;
						reg_new->gt_seq = "";
						reg_new->refseq =  refseq.at(leftClipRefPos-start_var_pos);
						reg_new->altseq = "";
						reg_new->AF = 0;
						reg_new->supp_num = supp_num;
						reg_new->DP = DP_num;
						reg_new->dup_num = dup_num;
						reg_new->discover_level = VAR_DISCOV_L_READS2;

						var_vec.push_back(reg_new);
					}
					break;
				case VAR_INV:
	//				cout << "line=" << __LINE__ << ", sv_type=" << sv_type << ", chrname=" << chrname << ", leftClipRefPos=" << leftClipRefPos << ", rightClipRefPos=" << rightClipRefPos << ", large_indel_flag=" << large_indel_flag << endl;
	//				cout << "sv_type=" << sv_type << ", dup_type_num=" << dup_type_num << ", inv_type_num=" << inv_type_num << ", tra_type_num=" << tra_type_num << endl;
					supp_num = inv_type_num;
					if(supp_num>DP_num) DP_num = supp_num;
					if(supp_num>=minReadsNumSupportSV){
						// get element of both margins matched
						for(i=0; i<qname_inv_vec.size(); i++){
							queryname = qname_inv_vec.at(i);
							query_aln_segs = getQueryClipAlnSegsAll(queryname, clipAlnDataVector_sub);
							//query_aln_segs = getQueryClipAlnSegs(queryname, clipAlnDataVector_sub); // deleted on 2026-04-16
							if(query_aln_segs.size()>=3){

	//							if(queryname.compare("SRR8955953.10734928")==0){
	//								cout << "line=" << __LINE__ << ", i=" << i << ", queryname=" << queryname << ", aln_size=" << query_aln_segs.size() << endl;
	//							}else continue;

								sortQueryAlnSegs(query_aln_segs);
								FilteredbyAlignmentSegment(query_aln_segs, chrname, minMapQ, MIN_ALN_SEG_SIZE_CLIPREG);
								if(query_aln_segs.size()==0) continue;

								//id_vec = getLeftRightClipAlnId(leftClipRefPos, rightClipRefPos, query_aln_segs); // deleted on 2026-04-16
								id_vec = getLeftRightClipAlnIdINV(leftClipRefPos, rightClipRefPos, query_aln_segs);
								if(id_vec.size()>0){
									clip_aln_left = query_aln_segs.at(id_vec.at(0));
									clip_aln_right = query_aln_segs.at(id_vec.at(1));

									entire_flanking_flag = false;
									clip_aln = NULL;
									if(clip_aln_left->aln_orient!=clip_aln_right->aln_orient){ // different orient
										// get the inverted item
										for(j=0; j<2; j++){
											if(j==0) clip_aln = clip_aln_left;
											else clip_aln = clip_aln_right;
											if(clip_aln->left_aln and clip_aln->right_aln){
												overlap_size = getOverlapSize(clip_aln->startRefPos, clip_aln->endRefPos, leftClipRefPos, rightClipRefPos);
												if(overlap_size<rightClipRefPos-leftClipRefPos+1)
													val = (double)overlap_size / (rightClipRefPos - leftClipRefPos + 1);
												else
													val = (double)(rightClipRefPos - leftClipRefPos + 1) / overlap_size;
												if(clip_aln->ref_dist<rightClipRefPos-leftClipRefPos+1)
													val2 = (double)clip_aln->ref_dist / (rightClipRefPos - leftClipRefPos + 1);
												else
													val2 = (double)(rightClipRefPos - leftClipRefPos + 1) / clip_aln->ref_dist;

												if(val>=max_seg_size_ratio and val2>=VALID_DUP_INV_SIZE_RATIO_CLIP_REG){
													entire_flanking_flag = true;
													//cout << "INV sim_val=" << val << endl;
													break;
												}
											}
										}
									}else{ // same orient, then try to find the entire flanking item
										start_var_query_pos = clip_aln_left->endQueryPos;
										end_var_query_pos = clip_aln_right->startQueryPos;
										start_var_ref_pos = clip_aln_left->endRefPos;
										end_var_ref_pos = clip_aln_right->startRefPos;
										if(clip_aln_left->aln_orient==ALN_PLUS_ORIENT) dist_query = end_var_query_pos - start_var_query_pos + 1;
										else dist_query = start_var_query_pos - end_var_query_pos + 1;
										dist_ref = end_var_ref_pos - start_var_ref_pos + 1;
										if(dist_query<rightClipRefPos-leftClipRefPos+1)
											val = (double)dist_query / (rightClipRefPos - leftClipRefPos + 1);
										else
											val = (double)(rightClipRefPos - leftClipRefPos + 1) / dist_query;
										if(val>=max_seg_size_ratio){
											maxVal = 0;
											maxIdx = -1;
											for(j=0; j<query_aln_segs.size(); j++){
												clip_aln = query_aln_segs.at(j);
												if(clip_aln->aln_orient!=clip_aln_left->aln_orient){
													overlap_size = getOverlapSize(clip_aln->startQueryPos, clip_aln->endQueryPos, start_var_query_pos, end_var_query_pos);
													if(overlap_size<dist_query)
														val = (double)overlap_size / dist_query;
													else
														val = (double)dist_query / overlap_size;
													if(clip_aln->ref_dist<dist_ref)
														val2 = (double)clip_aln->ref_dist / dist_ref;
													else
														val2 = (double)dist_ref / clip_aln->ref_dist;
													if(val>=max_seg_size_ratio and val2>=VALID_DUP_INV_SIZE_RATIO_CLIP_REG){
														if(maxVal<val2) {
															maxVal = val2;
															maxIdx = j;
														}
													}
												}
											}
											if(maxVal>=max_seg_size_ratio){
	//											cout << "INV sim_val=" << maxVal << endl;
												entire_flanking_flag = true;
												clip_aln = query_aln_segs.at(maxIdx);
											}
										}
									}

									if(entire_flanking_flag and clip_aln!=NULL){ // found
										if(clip_aln->startRefPos>=start_var_pos){
											queryseq = "";
											if(clip_aln->bam){
												seq_int = bam_get_seq(clip_aln->bam);
												start_qpos_tmp = clip_aln->startQueryPos;
												end_qpos_tmp = clip_aln->endQueryPos;
												if(clip_aln->leftHardClippedFlag){ // if hard-clipping, get the correct loci
													if(clip_aln->aln_orient==ALN_PLUS_ORIENT){
														start_qpos_tmp -= clip_aln->leftClipSize;
														end_qpos_tmp -= clip_aln->leftClipSize;
													}else{
														start_qpos_tmp = clip_aln->querylen - clip_aln->endQueryPos - clip_aln->leftClipSize + 1;
														end_qpos_tmp = clip_aln->querylen - clip_aln->startQueryPos - clip_aln->leftClipSize + 1;
													}
												}
												//if(clip_aln->aln_orient==ALN_PLUS_ORIENT) for(j=clip_aln->startQueryPos - 1; j<(size_t)clip_aln->endQueryPos; j++) queryseq += "=ACMGRSVTWYHKDBN"[bam_seqi(seq_int, j)];  // seq
												//else for(j=clip_aln->endQueryPos - 1; j<(size_t)clip_aln->startQueryPos; j++) queryseq += "=ACMGRSVTWYHKDBN"[bam_seqi(seq_int, j)];  // seq
												if(clip_aln->aln_orient==ALN_PLUS_ORIENT) for(j=start_qpos_tmp - 1; j<(size_t)end_qpos_tmp; j++) queryseq += "=ACMGRSVTWYHKDBN"[bam_seqi(seq_int, j)];  // seq
												else for(j=end_qpos_tmp - 1; j<(size_t)start_qpos_tmp; j++) queryseq += "=ACMGRSVTWYHKDBN"[bam_seqi(seq_int, j)];  // seq

											}

											// allocate new node
											reg_new = new reg_t();
											reg_new->var_type = sv_type;
											reg_new->chrname = clip_aln->chrname;
											reg_new->startRefPos = clip_aln->startRefPos;
											reg_new->endRefPos = clip_aln->endRefPos;
											reg_new->startLocalRefPos = -1;
											reg_new->endLocalRefPos = -1;
											reg_new->startQueryPos = -1;
											reg_new->endQueryPos = -1;
											reg_new->query_id = -1;
											reg_new->sv_len = clip_aln->ref_dist;
											reg_new->minimap2_aln_id = -1;
											reg_new->call_success_status = true;
											reg_new->short_sv_flag = false;
											reg_new->zero_cov_flag = false;
											reg_new->aln_seg_end_flag = false;
											reg_new->query_pos_invalid_flag = false;
											reg_new->large_indel_flag = false;
											reg_new->merge_flag = false;
											reg_new->rescue_flag = true;
											reg_new->aln_orient = clip_aln->aln_orient;
											reg_new->gt_type = -1;
											reg_new->gt_seq = "";
											reg_new->AF = 0;
											//clip_reg_tmp->supp_num = queryname_vec.size();
											reg_new->supp_num = supp_num;
											reg_new->DP = 0;
											reg_new->discover_level = VAR_DISCOV_L_READS2;
											reg_new->refseq = refseq.substr(reg_new->startRefPos-start_var_pos, reg_new->endRefPos-reg_new->startRefPos + 1);
											reg_new->altseq = queryseq;

											switch(sv_type){
												case VAR_DUP: // DUP
													reg_new->dup_num = dup_num;
													reg_new->large_indel_flag = false;
													break;
												case VAR_INV: // INV
													reg_new->dup_num = 0;
													reg_new->large_indel_flag = false;
													break;
												case VAR_INS: // large INS or DEL
												case VAR_DEL:
													reg_new->dup_num = 0;
													reg_new->large_indel_flag = true;
													break;
											}

											var_vec.push_back(reg_new);
										}
										break;
									}else{
										// extract queries from clip align data vector
										query_seq_info_all = extractQueriesFromClipAlnDataVec(query_aln_segs, refseq, chrname, start_var_pos, end_var_pos, start_var_pos_pure, end_var_pos_pure, fai, minConReadLen, true, &mutex_fai);
	//									cout << "query_seq_info_all.size=" << query_seq_info_all.size() << endl;

										if(query_seq_info_all.size()==query_aln_segs.size()){
											qnode_left = query_seq_info_all.at(id_vec.at(0));
											qnode_right = query_seq_info_all.at(id_vec.at(1));

											// left side
											idx_alnSegs_left = -1;
											for(j=0; j<qnode_left->query_alnSegs.size(); j++){
												aln_seg = qnode_left->query_alnSegs.at(j);
												if(aln_seg->opflag==BAM_CMATCH or aln_seg->opflag==BAM_CEQUAL){
													if(aln_seg->startRpos<=leftClipRefPos-1){
														idx_alnSegs_left = j;
													}else{
														if(idx_alnSegs_left==-1 and aln_seg->startRpos-leftClipRefPos<=MAX_REF_DIST_SAME_CHR)
															idx_alnSegs_left = j;
														break;
													}
												}
											}

											// right side
											idx_alnSegs_right = -1;
											for(j=0; j<qnode_right->query_alnSegs.size(); j++){
												aln_seg = qnode_right->query_alnSegs.at(j);
												if(aln_seg->opflag==BAM_CMATCH or aln_seg->opflag==BAM_CEQUAL){
													if(aln_seg->startRpos<=rightClipRefPos){
														idx_alnSegs_right = j;
													}else{
														if(idx_alnSegs_right==-1 and aln_seg->startRpos-rightClipRefPos<=MAX_REF_DIST_SAME_CHR)
															idx_alnSegs_right = j;
														break;
													}
												}
											}

											if(idx_alnSegs_left!=-1 and idx_alnSegs_right!=-1){ // found both margins
												margin_dist1 = leftClipRefPos - 1 - qnode_left->query_alnSegs.at(idx_alnSegs_left)->startRpos;
												if(margin_dist1<0) margin_dist1 = 0;
												left_rpos = qnode_left->query_alnSegs.at(idx_alnSegs_left)->startRpos + margin_dist1;
												//left_subjectpos = qnode_left->query_alnSegs.at(idx_alnSegs_left)->startSubjectPos + margin_dist1;
												left_qpos = qnode_left->query_alnSegs.at(idx_alnSegs_left)->startQpos + margin_dist1;

												margin_dist2 = rightClipRefPos - qnode_right->query_alnSegs.at(idx_alnSegs_right)->startRpos;
												if(margin_dist2<0) margin_dist2 = 0;
												right_rpos = qnode_right->query_alnSegs.at(idx_alnSegs_right)->startRpos + margin_dist2;
		//										right_subjectpos = minimap2_aln_right->pafalnsegs.at(idx_alnSegs_right)->startSubjectPos + margin_dist2;
												right_qpos = qnode_right->query_alnSegs.at(idx_alnSegs_right)->startQpos + margin_dist2;

												if(left_qpos<right_qpos and left_rpos<right_rpos){
													//sv_len = (right_qpos - left_qpos) - (right_rpos - left_rpos) + 1;
													noHardClipIdx = getNoHardClipAlnItem(query_aln_segs, chrname, start_var_pos, end_var_pos);
													//queryseq = qnode_left->seq;
													if(noHardClipIdx!=-1){
														queryseq = query_seq_info_all.at(noHardClipIdx)->seq;
														if(query_aln_segs.at(noHardClipIdx)->aln_orient!=qnode_left->clip_aln->aln_orient)
															reverseComplement(queryseq);
													}else
														queryseq = qnode_left->seq;

													if(!query_seq_info_all.empty()) destoryQuerySeqInfoAll(query_seq_info_all);

													if(left_qpos-1+right_qpos-left_qpos+1<(int64_t)queryseq.size()){
														tmp_refseq = refseq.substr(left_rpos-start_var_pos, right_rpos-left_rpos+1);
														tmp_queryseq = queryseq.substr(left_qpos-1, right_qpos-left_qpos+1);

														if(id_vec.at(0)==id_vec.at(1)) // same segment
															val = computeVarseqSim(tmp_queryseq, tmp_refseq);
														else{ // different segments
															reverseComplement(tmp_queryseq);
															val = computeVarseqSim(tmp_queryseq, tmp_refseq);
															//cout << "seqsim=" << val << endl;
														}

														if(val>=min_seqsim_match){ // similarity confirm
															// allocate new node
															reg_new = new reg_t();
															reg_new->var_type = sv_type;
															reg_new->chrname = clip_aln_left->chrname;
															reg_new->startRefPos = left_rpos;
															reg_new->endRefPos = right_rpos;
															reg_new->startLocalRefPos = -1;
															reg_new->endLocalRefPos = -1;
															reg_new->startQueryPos = left_qpos;
															reg_new->endQueryPos = right_qpos;
															reg_new->query_id = -1;
															reg_new->sv_len = right_rpos - left_rpos + 1;
															//if(sv_len<0) clip_reg_tmp->sv_len = -sv_len;
															//clip_reg_tmp->dup_num = dup_num; //--------------
															reg_new->minimap2_aln_id = -1;
															reg_new->call_success_status = true;
															reg_new->short_sv_flag = false;
															reg_new->zero_cov_flag = false;
															reg_new->aln_seg_end_flag = false;
															reg_new->query_pos_invalid_flag = false;
															reg_new->large_indel_flag = false;
															reg_new->merge_flag = false;
															reg_new->rescue_flag = true;
															reg_new->aln_orient = clip_aln_left->aln_orient;
															reg_new->gt_type = -1;
															reg_new->gt_seq = "";
															reg_new->AF = 0;
															//clip_reg_tmp->supp_num = queryname_vec.size();
															reg_new->supp_num = supp_num;
															reg_new->DP = 0;
															reg_new->discover_level = VAR_DISCOV_L_CNS_ALN;
															reg_new->refseq = tmp_refseq;
															reg_new->altseq = tmp_queryseq;

															// deal with identical ALT and REF alleles, change to symbolic type, added on 2025-11-16
															if(val==1){
																upperSeq(tmp_refseq);
																upperSeq(tmp_queryseq);
																if(tmp_refseq.compare(tmp_queryseq)==0){
																	reg_new->refseq = reg_new->refseq.at(0);
																	reg_new->altseq = "";
																}
															}

															switch(sv_type){
																case VAR_DUP: // DUP
																	reg_new->dup_num = dup_num;
																	reg_new->large_indel_flag = false;
																	break;
																case VAR_INV: // INV
																	reg_new->dup_num = 0;
																	reg_new->large_indel_flag = false;
																	break;
																case VAR_INS: // large INS or DEL
																case VAR_DEL:
																	reg_new->dup_num = 0;
																	reg_new->large_indel_flag = true;
																	break;
															}

															var_vec.push_back(reg_new);
															break;
														}
													}
												}
											}
											if(!query_seq_info_all.empty()) destoryQuerySeqInfoAll(query_seq_info_all);
										}
									}
								}
							}
						}

						// otherwise, get both margins
						if(var_vec.size()==0){
							for(i=0; i<qname_inv_vec.size(); i++){
								queryname = qname_inv_vec.at(i);
								//query_aln_segs = getQueryClipAlnSegsAll(queryname, clipAlnDataVector_sub);
								query_aln_segs = getQueryClipAlnSegs(queryname, clipAlnDataVector_sub);
								if(query_aln_segs.size()>=2){

//									if(queryname.compare("m84039_230414_235240_s2/160567057/ccs")==0){
//										cout << "queryname=" << queryname << ", aln_size=" << query_aln_segs.size() << endl;
//									}else continue;

									sortQueryAlnSegs(query_aln_segs);
									FilteredbyAlignmentSegment(query_aln_segs, chrname, minMapQ, MIN_ALN_SEG_SIZE_CLIPREG);
									if(query_aln_segs.size()==0) continue;

									//id_vec = getLeftRightClipAlnId(leftClipRefPos, rightClipRefPos, query_aln_segs);
									id_vec = getLeftRightClipAlnIdINV(leftClipRefPos, rightClipRefPos, query_aln_segs);
									//cout << "id_vec.size=" << id_vec.size() << endl;
									if(id_vec.size()>0){
										left_clip_rpos = right_clip_rpos = left_clip_qpos = right_clip_qpos = clip_end1 = clip_end2 = -1;
										clip_aln_left = query_aln_segs.at(id_vec.at(0));
										clip_aln_right = query_aln_segs.at(id_vec.at(1));

										if(clip_aln_left->startRefPos>clip_aln_right->startRefPos){
	//										cout << "********** Exchanged: queryname=" << queryname << ", clip_aln_left.startRefPos=" << clip_aln_left->startRefPos << ", clip_aln_right.startRefPos=" << clip_aln_right->startRefPos << endl;
											clip_aln_left = query_aln_segs.at(id_vec.at(1));
											clip_aln_right = query_aln_segs.at(id_vec.at(0));
										}

										margin_dist1 = clip_aln_left->startRefPos - leftClipRefPos;
										margin_dist2 = clip_aln_left->endRefPos - leftClipRefPos;
										if(abs(margin_dist1)<MAX_INV_ALN_MARGIN_DIST or abs(margin_dist2)<MAX_INV_ALN_MARGIN_DIST){ // otherwise, there maybe INVDUP
											if(clip_aln_left->right_aln and (abs(margin_dist1)>abs(margin_dist2) or clip_aln_left->left_aln==NULL)){
												left_clip_rpos = clip_aln_left->endRefPos;
												left_clip_qpos = clip_aln_left->endQueryPos;
												clip_end1 = RIGHT_END;
											}else if(clip_aln_left->left_aln and (abs(margin_dist1)<abs(margin_dist2) or clip_aln_left->right_aln==NULL)){
												left_clip_rpos = clip_aln_left->startRefPos;
												left_clip_qpos = clip_aln_left->startQueryPos;
												clip_end1 = LEFT_END;
											}
										}

										margin_dist1 = clip_aln_right->startRefPos - rightClipRefPos;
										margin_dist2 = clip_aln_right->endRefPos - rightClipRefPos;
										if(abs(margin_dist1)<MAX_INV_ALN_MARGIN_DIST or abs(margin_dist2)<MAX_INV_ALN_MARGIN_DIST){ // otherwise, there maybe INVDUP
											if(clip_aln_right->right_aln and (abs(margin_dist1)>abs(margin_dist2) or clip_aln_right->left_aln==NULL)){
												right_clip_rpos = clip_aln_right->endRefPos;
												right_clip_qpos = clip_aln_right->endQueryPos;
												clip_end2 = RIGHT_END;
											}else if(clip_aln_right->left_aln and (abs(margin_dist1)<abs(margin_dist2) or clip_aln_right->right_aln==NULL)){
												right_clip_rpos = clip_aln_right->startRefPos;
												right_clip_qpos = clip_aln_right->startQueryPos;
												clip_end2 = LEFT_END;
											}
										}

										if(clip_end1!=-1 and clip_end2!=-1){
											overlap_size = getOverlapSize(clip_aln_left->startQueryPos, clip_aln_left->endQueryPos, clip_aln_right->startQueryPos, clip_aln_right->endQueryPos);
											//cout << "overlap_size=" << overlap_size << ", clip_end1=" << clip_end1 << ", clip_end2=" << clip_end2 << endl;

											if(overlap_size>EXT_SIZE_CHK_VAR_LOC){
												// extract queries from clip align data vector
												query_seq_info_all = extractQueriesFromClipAlnDataVec(query_aln_segs, refseq, chrname, start_var_pos, end_var_pos, start_var_pos_pure, end_var_pos_pure, fai, minConReadLen, true, &mutex_fai);
		//										cout << "query_seq_info_all.size=" << query_seq_info_all.size() << endl;

												if(query_seq_info_all.at(0)->clip_aln==clip_aln_left){
													qnode_left = query_seq_info_all.at(0);
													qnode_right = query_seq_info_all.at(1);
												}else{
													qnode_left = query_seq_info_all.at(1);
													qnode_right = query_seq_info_all.at(0);
												}

												// choose the segment to refine
												if(clip_end1==RIGHT_END and clip_end2==RIGHT_END){ // choose the second segment
		//											cout << "111111 overlap_size=" << overlap_size << ", clip_end1=" << clip_end1 << ", clip_end2=" << clip_end2 << endl;

													// right side
													right_clip_rpos_new = -1;
													for(j=0; j<qnode_right->query_alnSegs.size(); j++){
														aln_seg = qnode_right->query_alnSegs.at(j);
														if(clip_aln_right->aln_orient==ALN_PLUS_ORIENT) {
															start_qpos_tmp = aln_seg->startQpos;
															end_qpos_tmp = getEndQueryPosAlnSeg(aln_seg->startQpos, aln_seg->opflag, aln_seg->seglen);
														}else{
															start_qpos_tmp = clip_aln_right->querylen - getEndQueryPosAlnSeg(aln_seg->startQpos, aln_seg->opflag, aln_seg->seglen) + 1;
															end_qpos_tmp = clip_aln_right->querylen - aln_seg->startQpos + 1;
														}
														//cout << "start_qpos_tmp=" << start_qpos_tmp << ", end_qpos_tmp=" << end_qpos_tmp << endl;
														if(start_qpos_tmp<=left_clip_qpos and end_qpos_tmp>=left_clip_qpos){
															margin_dist1 = left_clip_qpos - start_qpos_tmp;
															if(clip_aln_right->aln_orient==ALN_PLUS_ORIENT)
																right_clip_rpos_new = aln_seg->startRpos + margin_dist1;
															else
																right_clip_rpos_new = aln_seg->startRpos + (aln_seg->seglen - margin_dist1);
															break;
														}
													}
		//											cout << "right_clip_rpos=" << right_clip_rpos << ", right_clip_rpos_new=" << right_clip_rpos_new << endl;
													if(right_clip_rpos_new!=-1)
														right_clip_rpos = right_clip_rpos_new;
												}else if(clip_end1==LEFT_END and clip_end2==LEFT_END){ // choose the first segment
		//											cout << "222222 overlap_size=" << overlap_size << ", clip_end1=" << clip_end1 << ", clip_end2=" << clip_end2 << endl;

													// left side
													left_clip_rpos_new = -1;
													for(j=0; j<qnode_left->query_alnSegs.size(); j++){
														aln_seg = qnode_left->query_alnSegs.at(j);
														if(clip_aln_left->aln_orient==ALN_PLUS_ORIENT) {
															start_qpos_tmp = aln_seg->startQpos;
															end_qpos_tmp = getEndQueryPosAlnSeg(aln_seg->startQpos, aln_seg->opflag, aln_seg->seglen);
														}else{
															start_qpos_tmp = clip_aln_left->querylen - getEndQueryPosAlnSeg(aln_seg->startQpos, aln_seg->opflag, aln_seg->seglen) + 1;
															end_qpos_tmp = clip_aln_left->querylen - aln_seg->startQpos + 1;
														}
														//cout << "start_qpos_tmp=" << start_qpos_tmp << ", end_qpos_tmp=" << end_qpos_tmp << endl;
														if(start_qpos_tmp<=right_clip_qpos and end_qpos_tmp>=right_clip_qpos){
															margin_dist1 = right_clip_qpos - start_qpos_tmp;
															if(clip_aln_left->aln_orient==ALN_PLUS_ORIENT)
																left_clip_rpos_new = aln_seg->startRpos + margin_dist1;
															else
																left_clip_rpos_new = aln_seg->startRpos + (aln_seg->seglen - margin_dist1);
															break;
														}
													}
		//											cout << "left_clip_rpos=" << left_clip_rpos << ", left_clip_rpos_new=" << left_clip_rpos_new << endl;
													if(left_clip_rpos_new!=-1)
														left_clip_rpos = left_clip_rpos_new;
												}else{
		//											cout << "333333 overlap_size=" << overlap_size << ", clip_end1=" << clip_end1 << ", clip_end2=" << clip_end2 << endl;

												}
												if(!query_seq_info_all.empty()) destoryQuerySeqInfoAll(query_seq_info_all);
											}
											//cout << "-------- left_clip_rpos=" << left_clip_rpos << ", right_clip_rpos=" << right_clip_rpos << endl;
											if(left_clip_rpos<=right_clip_rpos){
												left_clip_pos_vec.push_back(left_clip_rpos);
												right_clip_pos_vec.push_back(right_clip_rpos);
											}else{
												left_clip_pos_vec.push_back(right_clip_rpos);
												right_clip_pos_vec.push_back(left_clip_rpos);
											}
										}
									}
								}
							}
	//						cout << "left_clip_pos_vec.size=" << left_clip_pos_vec.size() << ", right_clip_pos_vec.size=" << right_clip_pos_vec.size() << endl;
							if(left_clip_pos_vec.size()>0 and right_clip_pos_vec.size()>0){
								mean_left_clip_rpos = accumulate(left_clip_pos_vec.begin(), left_clip_pos_vec.end(), 0.0) / left_clip_pos_vec.size();
								mean_right_clip_rpos = accumulate(right_clip_pos_vec.begin(), right_clip_pos_vec.end(), 0.0) / right_clip_pos_vec.size();
	//							cout << "mean_left_clip_rpos=" << mean_left_clip_rpos << ", mean_right_clip_rpos=" << mean_right_clip_rpos << endl;

								double squared_diff_sum = 0.0;
								for (int64_t value : left_clip_pos_vec) squared_diff_sum += (value - mean_left_clip_rpos) * (value - mean_left_clip_rpos);
								double sdev1 = sqrt(squared_diff_sum / left_clip_pos_vec.size());

								squared_diff_sum = 0.0;
								for (int64_t value : right_clip_pos_vec) squared_diff_sum += (value - mean_right_clip_rpos) * (value - mean_right_clip_rpos);
								double sdev2 = sqrt(squared_diff_sum / right_clip_pos_vec.size());

	//							cout << "sdev1=" << sdev1 << ", sdev2=" << sdev2 << endl;

								double dif;
								for (i=0; i<left_clip_pos_vec.size(); ){
									left_clip_rpos = left_clip_pos_vec.at(i);
									dif = abs(left_clip_rpos - mean_left_clip_rpos);
									if((dif>2*sdev1 and dif>CLIP_END_EXTEND_SIZE) or dif>VAR_ALN_EXTEND_SIZE){
										left_clip_pos_vec.erase(left_clip_pos_vec.begin()+i);
									}else i++;
								}
								for (i=0; i<right_clip_pos_vec.size(); ){
									right_clip_rpos = right_clip_pos_vec.at(i);
									dif = abs(right_clip_rpos - mean_right_clip_rpos);
									if((dif>2*sdev2 and dif>CLIP_END_EXTEND_SIZE) or dif>VAR_ALN_EXTEND_SIZE){
										right_clip_pos_vec.erase(right_clip_pos_vec.begin()+i);
									}else i++;
								}

								if(left_clip_pos_vec.size()>=(size_t)minReadsNumSupportSV and right_clip_pos_vec.size()>=(size_t)minReadsNumSupportSV){
									mean_left_clip_rpos = accumulate(left_clip_pos_vec.begin(), left_clip_pos_vec.end(), 0.0) / left_clip_pos_vec.size();
									mean_right_clip_rpos = accumulate(right_clip_pos_vec.begin(), right_clip_pos_vec.end(), 0.0) / right_clip_pos_vec.size();
		//							cout << "mean_left_clip_rpos=" << mean_left_clip_rpos << ", mean_right_clip_rpos=" << mean_right_clip_rpos << endl;
									if(mean_left_clip_rpos>mean_right_clip_rpos){ pos_tmp = mean_left_clip_rpos; mean_left_clip_rpos = mean_right_clip_rpos; mean_right_clip_rpos = pos_tmp; }

									if((int64_t)mean_left_clip_rpos>=start_var_pos){
										// allocate new node
										reg_new = new reg_t();
										reg_new->var_type = sv_type;
										reg_new->chrname = chrname;
										reg_new->startRefPos = mean_left_clip_rpos;
										reg_new->endRefPos = mean_right_clip_rpos;
										reg_new->startLocalRefPos = -1;
										reg_new->endLocalRefPos = -1;
										reg_new->startQueryPos = -1;
										reg_new->endQueryPos = -1;
										reg_new->query_id = -1;
										reg_new->sv_len = mean_right_clip_rpos - mean_left_clip_rpos + 1;
										reg_new->minimap2_aln_id = -1;
										reg_new->call_success_status = true;
										reg_new->short_sv_flag = false;
										reg_new->zero_cov_flag = false;
										reg_new->aln_seg_end_flag = false;
										reg_new->query_pos_invalid_flag = false;
										reg_new->large_indel_flag = false;
										reg_new->merge_flag = false;
										reg_new->rescue_flag = true;
										reg_new->aln_orient = -1;
										reg_new->gt_type = -1;
										reg_new->gt_seq = "";
										reg_new->AF = 0;
										//clip_reg_tmp->supp_num = queryname_vec.size();
										reg_new->supp_num = supp_num;
										reg_new->DP = 0;
										reg_new->discover_level = VAR_DISCOV_L_READS2;
										if(mean_left_clip_rpos-start_var_pos<refseq.size()) reg_new->refseq = refseq.at(reg_new->startRefPos-start_var_pos);
										else{
											reg_str = chrname + ":" + to_string(reg_new->startRefPos) + "-" + to_string(reg_new->startRefPos);
											pthread_mutex_lock(&mutex_fai);
											p_seq = fai_fetch(fai, reg_str.c_str(), &seq_len);
											pthread_mutex_unlock(&mutex_fai);
											reg_new->refseq = p_seq;
											free(p_seq);
										}
										reg_new->altseq = "";

										reg_new->dup_num = 0;
										reg_new->large_indel_flag = false;

										var_vec.push_back(reg_new);
									}
								}
							}
						}
					}
					break;
				case VAR_TRA:
				case VAR_BND:
					cout << "line=" << __LINE__ << ", sv_type=" << sv_type << ", chrname=" << chrname << ", leftClipRefPos=" << leftClipRefPos << ", rightClipRefPos=" << rightClipRefPos << ", large_indel_flag=" << large_indel_flag << endl;

					break;
			}
		}

		if(var_vec.size()==0 and (sv_type==VAR_INS or sv_type==VAR_DEL)){ // INS, DEL
			//if(large_indel_flag and abs(varVec.at(0)->sv_len)>=EXT_SIZE_CHK_VAR_LOC_SMALL){
			if(large_indel_flag){

#if READS_CALL_DEBUG
				cout << "line=" << __LINE__ << ", sv_type=" << sv_type << ", sv_len=" << varVec.at(0)->sv_len << ", chrname=" << chrname << ", leftClipRefPos=" << leftClipRefPos << ", rightClipRefPos=" << rightClipRefPos << ", large_indel_flag=" << large_indel_flag << endl;
#endif
				// extract queries from clip align data vector
				query_seq_info_all = extractQueriesFromClipAlnDataVec(clipAlnDataVector_sub, refseq, chrname, start_var_pos, end_var_pos, start_var_pos_pure, end_var_pos_pure, fai, minConReadLen, clip_reg_flag, &mutex_fai);

				// extract the features
				clipRegCluster clip_reg_cluster(chrname, leftClipRefPos, rightClipRefPos, minClipEndSize, maxVarRegSize, min_sv_size, minReadsNumSupportSV, minMapQ, min_seqsim_match, technology, fai);
				qcSigList_vec = clip_reg_cluster.extractQcSigsClipReg(query_seq_info_all);
				clip_reg_cluster.prepareQcSigListInfoClipRegForCluster(qcSigList_vec);

#if READS_CALL_DEBUG
				clip_reg_cluster.printQcSigListVec(qcSigList_vec);
#endif

				var_vec_sub.clear();
				for(i=0; i<varVec.size(); i++){
					if(var_flag_arr[i]==0){
						reg_tmp = varVec.at(i);
						sv_len = abs(reg_tmp->sv_len);
						ins_type_num = del_type_num = 0;
						qc_sig_vec.clear();
						for(j=0; j<qcSigList_vec.size(); j++){
							qcSigList_node = qcSigList_vec.at(j);
							if(qcSigList_node->qcSig_vec.size()>0){
								val = -1;
								if(sv_type==VAR_INS and qcSigList_node->ins_sum>0){ // INS
									for(const auto& sig : qcSigList_node->qcSig_vec){
										if(sig->cigar_op==BAM_CINS){
											if(sig->cigar_op_len<sv_len) val = (double)sig->cigar_op_len / sv_len;
											else val = (double)sv_len / sig->cigar_op_len;
											if(val>=min_seqsim_match){
												ins_type_num ++;
												qc_sig_vec.push_back(sig);
											}
										}
									}
								}else if(sv_type==VAR_DEL and qcSigList_node->del_sum>0){ // DEL
									for(const auto& sig : qcSigList_node->qcSig_vec){
										if(sig->cigar_op==BAM_CDEL){
											if(sig->cigar_op_len<sv_len) val = (double)sig->cigar_op_len / sv_len;
											else val = (double)sv_len / sig->cigar_op_len;
											if(val>=min_seqsim_match){
												del_type_num ++;
												qc_sig_vec.push_back(sig);
											}
										}
									}
								}
							}
						}

#if READS_CALL_DEBUG
						cout << "sv_type=" << sv_type << ", ins_type_num=" << ins_type_num << ", del_type_num=" << del_type_num << endl;
#endif

						if(qc_sig_vec.size()>0){
							// calculate the mean size
							sum_len = 0;
							for(auto const& sig : qc_sig_vec) sum_len += sig->cigar_op_len;
							mean_size = (double)sum_len / qc_sig_vec.size();

							min_dif_idx = -1;
							min_dif = LONG_MAX;
							for(j=0; j<qc_sig_vec.size(); j++){
								dif = abs(qc_sig_vec.at(j)->cigar_op_len-mean_size);
								//cout << "dif=" << dif << ", mean_size=" << mean_size << ", cigar_op_len=" << qc_sig_vec.at(j)->cigar_op_len << endl;
								if(dif<min_dif){
									min_dif = dif;
									min_dif_idx = j;
								}
							}

							if(min_dif_idx!=-1){
								supp_num = sv_len = 0;
								DP_num = depth_largeIndel;
								qc_sig = qc_sig_vec.at(min_dif_idx);
								switch(sv_type){
									case VAR_INS:
										supp_num = ins_type_num;
										if(supp_num>DP_num) DP_num = supp_num;
										sv_len = qc_sig->cigar_op_len;
										break;
									case VAR_DEL:
										supp_num = del_type_num;
										if(supp_num>DP_num) DP_num = supp_num;
										sv_len = -qc_sig->cigar_op_len;
										break;
								}

								//if(supp_num>=minReadsNumSupportSV){ // deleted on 2026-04-13
								if(supp_num>=minReadsNumSupportSV*READS_NUM_SUPPORT_FACTOR){
									reg_new = new reg_t();
									reg_new->var_type = sv_type;
									reg_new->chrname = chrname;
									reg_new->startRefPos = qc_sig->ref_pos;
									reg_new->endRefPos = getEndRefPosAlnSeg(qc_sig->ref_pos, qc_sig->cigar_op, qc_sig->cigar_op_len);
									reg_new->startLocalRefPos = -1;
									reg_new->endLocalRefPos = -1;
									reg_new->startQueryPos = -1;
									reg_new->endQueryPos = -1;
									reg_new->query_id = -1;
									reg_new->sv_len = sv_len;
									reg_new->minimap2_aln_id = -1;
									reg_new->call_success_status = true;
									reg_new->short_sv_flag = false;
									reg_new->zero_cov_flag = false;
									reg_new->aln_seg_end_flag = false;
									reg_new->query_pos_invalid_flag = false;
									reg_new->large_indel_flag = false;
									reg_new->merge_flag = false;
									reg_new->rescue_flag = true;
									reg_new->aln_orient = ALN_PLUS_ORIENT;
									reg_new->gt_type = -1;
									reg_new->gt_seq = "";
									reg_new->refseq = qc_sig->refseq;
									reg_new->altseq = qc_sig->altseq;
									reg_new->AF = 0;
									reg_new->supp_num = supp_num;
									reg_new->DP = DP_num;
									reg_new->dup_num = dup_num;
									reg_new->discover_level = VAR_DISCOV_L_READS2;

									var_vec_sub.push_back(reg_new);
									var_flag_arr[i] = 1;
								}
							}
						}
					}
				}

				if(var_vec_sub.size()>0){
					for(auto const& item : var_vec_sub) var_vec.push_back(item);
				}else{ // local inner variants, ignore sv_type
					ins_sum_mean = del_sum_mean = ins_type_num = del_type_num = 0;
					for(j=0; j<qcSigList_vec.size(); j++){
						qcSigList_node = qcSigList_vec.at(j);
						if(qcSigList_node->ins_sum>0){ // INS
							ins_sum_mean += qcSigList_node->ins_sum;
							ins_type_num ++;
						}
						if(qcSigList_node->del_sum>0){ // DEL
							del_sum_mean += qcSigList_node->del_sum;
							del_type_num ++;
						}
					}

					// mean size
					if(ins_type_num>0) ins_sum_mean /= ins_type_num;
					if(del_type_num>0) del_sum_mean /= del_type_num;

#if READS_CALL_DEBUG
					cout << "sv_type=" << sv_type << ", ins_sum_mean=" << ins_sum_mean << ", ins_type_num=" << ins_type_num << ", del_sum_mean=" << del_sum_mean << ", del_type_num=" << del_type_num << endl;
#endif

					if(ins_sum_mean>del_sum_mean and ins_sum_mean>=2*min_sv_size and ins_type_num>=minReadsNumSupportSV){ // INS
						// get representative idx
						min_dif_idx = -1;
						min_dif = INT_MAX;
						for(j=0; j<qcSigList_vec.size(); j++){
							qcSigList_node = qcSigList_vec.at(j);
							if(qcSigList_node->ins_sum>0 and qcSigList_node->qcSig_vec.size()==1){
								dif = abs(qcSigList_node->ins_sum-ins_sum_mean);
								if(dif<min_dif){
									min_dif_idx = j;
									min_dif = dif;
								}
							}
						}
						if(min_dif_idx==-1){
							for(j=0; j<qcSigList_vec.size(); j++){
								qcSigList_node = qcSigList_vec.at(j);
								if(qcSigList_node->ins_sum>0){
									dif = abs(qcSigList_node->ins_sum-ins_sum_mean);
									if(dif<min_dif){
										min_dif_idx = j;
										min_dif = dif;
									}
								}
							}
						}

#if READS_CALL_DEBUG
						cout << "min_idx=" << min_dif_idx << ", min_dif=" << min_dif << endl;
#endif

						if(min_dif_idx!=-1){
							qcSigList_node = qcSigList_vec.at(min_dif_idx);
							for(auto const& sig : qcSigList_node->qcSig_vec){
								if(sig->cigar_op==BAM_CINS){
									//save
									reg_new = new reg_t();
									reg_new->var_type = VAR_INS;
									reg_new->chrname = chrname;
									reg_new->startRefPos = sig->ref_pos + 1;
									reg_new->endRefPos = getEndRefPosAlnSeg(sig->ref_pos, sig->cigar_op, sig->cigar_op_len) + 1;
									reg_new->startLocalRefPos = -1;
									reg_new->endLocalRefPos = -1;
									reg_new->startQueryPos = -1;
									reg_new->endQueryPos = -1;
									reg_new->query_id = -1;
									reg_new->sv_len = sig->cigar_op_len;
									reg_new->minimap2_aln_id = -1;
									reg_new->call_success_status = true;
									reg_new->short_sv_flag = false;
									reg_new->zero_cov_flag = false;
									reg_new->aln_seg_end_flag = false;
									reg_new->query_pos_invalid_flag = false;
									reg_new->large_indel_flag = false;
									reg_new->merge_flag = false;
									reg_new->rescue_flag = true;
									reg_new->aln_orient = -1;
									reg_new->gt_type = -1;
									reg_new->gt_seq = "";
									reg_new->refseq =  sig->refseq;
									reg_new->altseq = sig->altseq;
									reg_new->AF = 0;
									reg_new->supp_num = ins_type_num;
									reg_new->DP = computeCovNumReg(reg_new->chrname, reg_new->startRefPos, reg_new->endRefPos, fai, inBamFile, minMapQ, minHighMapQ, max_ultra_high_cov);
									reg_new->discover_level = VAR_DISCOV_L_READS;

									var_vec_inner.push_back(reg_new);
								}
							}
						}
					}else if(ins_sum_mean<del_sum_mean and del_sum_mean>=2*min_sv_size and del_type_num>=minReadsNumSupportSV){ // DEL
						// get representative idx
						min_dif_idx = -1;
						min_dif = INT_MAX;
						for(j=0; j<qcSigList_vec.size(); j++){
							qcSigList_node = qcSigList_vec.at(j);
							if(qcSigList_node->del_sum>0 and qcSigList_node->qcSig_vec.size()==1){
								dif = abs(qcSigList_node->del_sum-del_sum_mean);
								if(dif<min_dif){
									min_dif_idx = j;
									min_dif = dif;
								}
							}
						}
						if(min_dif_idx==-1){
							for(j=0; j<qcSigList_vec.size(); j++){
								qcSigList_node = qcSigList_vec.at(j);
								if(qcSigList_node->del_sum>0){
									dif = abs(qcSigList_node->del_sum-del_sum_mean);
									if(dif<min_dif){
										min_dif_idx = j;
										min_dif = dif;
									}
								}
							}
						}

#if READS_CALL_DEBUG
						cout << "min_idx=" << min_dif_idx << ", min_dif=" << min_dif << endl;
#endif

						if(min_dif_idx!=-1){
							qcSigList_node = qcSigList_vec.at(min_dif_idx);
							for(auto const& sig : qcSigList_node->qcSig_vec){
								if(sig->cigar_op==BAM_CDEL){
									//save
									reg_new = new reg_t();
									reg_new->var_type = VAR_DEL;
									reg_new->chrname = chrname;
									reg_new->startRefPos = sig->ref_pos + 1;
									reg_new->endRefPos = getEndRefPosAlnSeg(sig->ref_pos, sig->cigar_op, sig->cigar_op_len) + 1;
									reg_new->startLocalRefPos = -1;
									reg_new->endLocalRefPos = -1;
									reg_new->startQueryPos = -1;
									reg_new->endQueryPos = -1;
									reg_new->query_id = -1;
									reg_new->sv_len = sig->cigar_op_len;
									reg_new->minimap2_aln_id = -1;
									reg_new->call_success_status = true;
									reg_new->short_sv_flag = false;
									reg_new->zero_cov_flag = false;
									reg_new->aln_seg_end_flag = false;
									reg_new->query_pos_invalid_flag = false;
									reg_new->large_indel_flag = false;
									reg_new->merge_flag = false;
									reg_new->rescue_flag = true;
									reg_new->aln_orient = -1;
									reg_new->gt_type = -1;
									reg_new->gt_seq = "";
									reg_new->refseq =  sig->refseq;
									reg_new->altseq = sig->altseq;
									reg_new->AF = 0;
									reg_new->supp_num = del_type_num;
									reg_new->DP = computeCovNumReg(reg_new->chrname, reg_new->startRefPos, reg_new->endRefPos, fai, inBamFile, minMapQ, minHighMapQ, max_ultra_high_cov);
									reg_new->discover_level = VAR_DISCOV_L_READS;

									var_vec_inner.push_back(reg_new);
								}
							}
						}
					}
				}

				// release memory
				if(!query_seq_info_all.empty()) destoryQuerySeqInfoAll(query_seq_info_all);
				if(!qcSigList_vec.empty()) clip_reg_cluster.destoryQcSigList(qcSigList_vec);
			}
		}
	}

	// extract from mateClipReg item
	if(clip_reg_flag and var_vec.size()==0 and clusterId_incomplete.size()==querynames_vec.size()){
		if(sv_type==VAR_INV){
			supp_num = DP_num = INT_MAX;
			for(i=0; i<4; i++){
				if(bnd_mate_reg_strs[i].compare("-")!=0){
					str_vec = split(bnd_mate_reg_strs[i], ",");
					for(j=0; j<str_vec.size(); j++){
						if(str_vec.at(j).compare("-")!=0){
							sum = 0;
							str_vec2_tmp = split(str_vec.at(j), "_");
							for(const auto& item : str_vec2_tmp){
								str_vec2 = split(item, "|");
								supp_num_tmp = stoi(str_vec2.at(2));
								DP_num_tmp = stoi(str_vec2.at(3));
								sum += supp_num_tmp;
								//if(supp_num_tmp<supp_num) supp_num = supp_num_tmp;
								if(DP_num_tmp<DP_num) DP_num = DP_num_tmp;
							}
							if(sum<supp_num) supp_num = sum;
						}
					}
				}
			}
			//cout << "supp_num=" << supp_num << ", DP_num=" << DP_num << endl;

			if(supp_num!=INT_MAX and supp_num>=minReadsNumSupportSV){
				// allocate new node
				reg_new = new reg_t();
				reg_new->var_type = sv_type;
				reg_new->chrname = chrname;
				reg_new->startRefPos = leftClipRefPos;
				reg_new->endRefPos = rightClipRefPos;
				reg_new->startLocalRefPos = -1;
				reg_new->endLocalRefPos = -1;
				reg_new->startQueryPos = -1;
				reg_new->endQueryPos = -1;
				reg_new->query_id = -1;
				reg_new->sv_len = rightClipRefPos - leftClipRefPos + 1;
				reg_new->minimap2_aln_id = -1;
				reg_new->call_success_status = true;
				reg_new->short_sv_flag = false;
				reg_new->zero_cov_flag = false;
				reg_new->aln_seg_end_flag = false;
				reg_new->query_pos_invalid_flag = false;
				reg_new->large_indel_flag = false;
				reg_new->merge_flag = false;
				reg_new->rescue_flag = true;
				reg_new->aln_orient = ALN_PLUS_ORIENT;
				reg_new->gt_type = -1;
				reg_new->gt_seq = "";
				reg_new->AF = 0;
				//clip_reg_tmp->supp_num = queryname_vec.size();
				reg_new->supp_num = supp_num;
				reg_new->DP = DP_num;
				reg_new->discover_level = VAR_DISCOV_L_READS2;
				reg_new->refseq = refseq.at(reg_new->startRefPos-start_var_pos);
				reg_new->altseq = "";

				switch(sv_type){
					case VAR_DUP: // DUP
						reg_new->dup_num = dup_num;
						reg_new->large_indel_flag = false;
						break;
					case VAR_INV: // INV
						reg_new->dup_num = 0;
						reg_new->large_indel_flag = false;
						break;
					case VAR_INS: // large INS or DEL
					case VAR_DEL:
						reg_new->dup_num = 0;
						reg_new->large_indel_flag = true;
						break;
				}

				var_vec.push_back(reg_new);
			}
			for(auto const& item : var_vec_inner) var_vec.push_back(item);
		}else if(large_indel_flag and (sv_type==VAR_INS or sv_type==VAR_DEL)){
			if(varVec.size()==1 and supp_num_largeIndel>=minReadsNumSupportSV){

				reg_tmp = varVec.at(0);
				sum_len = 0;
				for(auto const& item : var_vec_inner) if(reg_tmp->var_type==item->var_type) sum_len += abs(item->sv_len);

				val = (double)sum_len / abs(reg_tmp->sv_len);
				if(val<=MIN_SIZE_RATIO_RECOVER){
					reg_new = dupVarReg(varVec.at(0));
					reg_new->supp_num = supp_num_largeIndel;
					reg_new->DP = depth_largeIndel;
					reg_new->discover_level = VAR_DISCOV_L_READS2;
					reg_new->gt_type = -1;
					reg_new->gt_seq = "";
					reg_new->AF = 0;
					reg_new->refseq = refseq.at(reg_new->startRefPos-start_var_pos);
					reg_new->altseq = "";
					reg_new->dup_num = 0;
					reg_new->large_indel_flag = true;
					reg_new->call_success_status = true;

					var_vec.push_back(reg_new);
				}
			}
			for(auto const& item : var_vec_inner) var_vec.push_back(item);
		}
	}

	// release memory
	if(!clipAlnDataVector.empty()) destoryClipAlnData(clipAlnDataVector);
	if(!query_seq_info_all.empty()) destoryQuerySeqInfoAll(query_seq_info_all);

	return var_vec;
}

void varCand::updateClusterIdIncomplete(vector<reg_t*> &var_vec, vector<reg_t*> &var_vec_rescue, vector<int32_t> &clusterId_incomplete){
	size_t i, j;
	reg_t *reg;
	vector< vector<string> > qnames_vec_cluster;
	int32_t group_id;

	qnames_vec_cluster = getClusterInfo(clusterfilename);

	for(i=0; i<var_vec_rescue.size(); i++) {
		reg = var_vec_rescue.at(i);
		var_vec.push_back(reg);
		group_id = getReadGroupId(reg, qnames_vec_cluster, this);
		for(j=0; j<clusterId_incomplete.size(); j++){
			if(clusterId_incomplete.at(j)==group_id){
				clusterId_incomplete.erase(clusterId_incomplete.begin()+j);
				break;
			}
		}
	}
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
