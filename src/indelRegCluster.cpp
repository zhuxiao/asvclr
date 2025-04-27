#include "indelRegCluster.h"
#include <numeric>

extern pthread_mutex_t mutex_fai;

indelRegCluster::indelRegCluster(string &chrname, int64_t var_startRefPos, int64_t var_endRefPos, int32_t sv_len_est, int64_t startRefPos_cns, int64_t endRefPos_cns, int32_t min_sv_size, int32_t min_supp_num, double min_identity_match, faidx_t *fai) {
	this->chrname = chrname;
	this->var_startRefPos = var_startRefPos;
	this->var_endRefPos = var_endRefPos;
	this->sv_len_est = sv_len_est;
	this->startRefPos_cns = startRefPos_cns;
	this->endRefPos_cns = endRefPos_cns;
	this->min_sv_size = min_sv_size;
	this->min_supp_num = min_supp_num;
	this->chrlen = faidx_seq_len(fai, chrname.c_str()); // get reference size
	this->fai = fai;
	this->min_identity_match = min_identity_match;

	max_merge_span = CNS_MIN_DIST_MERGE_THRES;
	if(sv_len_est>CNS_MIN_DIST_MERGE_THRES) max_merge_span = 1.5 * CNS_MIN_DIST_MERGE_THRES;
}

indelRegCluster::~indelRegCluster() {
}

// query clustering for indel
vector<struct querySeqInfoVec*> indelRegCluster::queryCluster(vector<struct querySeqInfoNode*> &query_seq_info_all){
	vector<struct querySeqInfoVec*> query_clu_vec;
	size_t i;
	//int32_t max_group_num;

	for(i=0; i<query_seq_info_all.size(); i++){
		query_seq_info_all.at(i)->cluster_finished_flag = false;
		extractQcSigsFromAlnSegs(query_seq_info_all.at(i));
	}

	prepareQueryInfoForCluster(query_seq_info_all); // prepare query information for cluster

	// estimate number of the group
	//max_group_num = estimateClusterGroupNum(query_seq_info_all);

#if POA_ALIGN_DEBUG
	// print extracted signatures
	printQcSigIndel(query_seq_info_all);
	//cout << "Estimated max_group_num=" << max_group_num << endl;
#endif

//	if(max_group_num==1)
//		query_clu_vec = queryClusterSingleGroup(query_seq_info_all);
//	else
		query_clu_vec = queryClusterDoubleGroup(query_seq_info_all);

#if POA_ALIGN_DEBUG
	printQcSigIndelCluster(query_clu_vec);
#endif

	return query_clu_vec;
}

// extract query cluster signatures from alnSegs
void indelRegCluster::extractQcSigsFromAlnSegs(struct querySeqInfoNode* &query_seq_info_node){
	int64_t startSpanPos_extend, endSpanPos_extend, start_var_pos, end_var_pos;

//	start_var_pos = var_startRefPos - EXT_SIZE_CHK_VAR_LOC; // deleted on 2025-02-28
//	end_var_pos = var_endRefPos + EXT_SIZE_CHK_VAR_LOC;
	start_var_pos = var_startRefPos - max_merge_span;
	end_var_pos = var_endRefPos + max_merge_span;
	if(start_var_pos<1) start_var_pos = 1;
	if(end_var_pos>chrlen) end_var_pos = chrlen;

	// extract
//	startSpanPos_extend = startRefPos_cns;
//	endSpanPos_extend = endRefPos_cns;
	startSpanPos_extend = start_var_pos;
	endSpanPos_extend = end_var_pos;

	query_seq_info_node->qcSig_vec = extractQcSigsFromAlnSegsSingleQuery(query_seq_info_node, query_seq_info_node->clip_aln->chrname, startSpanPos_extend, endSpanPos_extend, min_sv_size);
	query_seq_info_node->qcSig_merge_flag = mergeNeighbouringSigsFlag(query_seq_info_node, query_seq_info_node->qcSig_vec, MAX_DIST_MERGE_ARBITARY, max_merge_span, CNS_MIN_IDENTITY_MERGE_THRES2, MIN_VALID_SIG_SIZE_RATIO_THRES, fai);
	filterSmallSigs(query_seq_info_node->qcSig_vec, MIN_SIG_SIZE_RATIO_FILTER_THRES);
	computeVarSums(query_seq_info_node, start_var_pos, end_var_pos);
}

// merge neighbouring signatures for indels
bool indelRegCluster::mergeNeighbouringSigsFlag(struct querySeqInfoNode* query_seq_info_node, vector<qcSig_t*> &sig_vec, int32_t min_ref_dist_thres, int32_t max_ref_dist_thres, double min_merge_identity_thres, double min_valid_sig_size_ratio_thres, faidx_t *fai){
	size_t i, j;
	int32_t dist, query_opflag, seq_len;
	int64_t l, query1_startRpos, query1_endRpos, query1_svlen, query1_startQpos, query1_endQpos, query2_startRpos, query2_svlen, query2_startQpos, query2_endQpos, start_querypos_comp, end_querypos_comp, start_refpos_comp, end_refpos_comp, start_querypos_merged, end_querypos_merged, query2_endRpos;
	int64_t pre_distance;
	uint8_t *seq_int;
	string comp_refseq, comp_queryseq, reg_str, merging_seq;
	double iden_val;
	qcSig_t *sig1, *sig2;
	char *seq;
	bool flag;

	dist = -1;
	flag = false;
	// if(query_seq_info_node->qname.compare("SRR11008518.1.11068758")==0)
		// cout << "--------------qname=" << query_seq_info_node->qname << "--------------" << endl;
	if(query_seq_info_node->clip_aln and query_seq_info_node->clip_aln->bam){
		seq_int = bam_get_seq(query_seq_info_node->clip_aln->bam);
		for(i=0; i<sig_vec.size(); i++){
			sig1 = sig_vec.at(i);
			// skip normal signature
			if(sig1->cigar_op!=BAM_CINS and sig1->cigar_op!=BAM_CDEL) continue;

			query_opflag = sig1->cigar_op;
			// pre_refpos_comp = -1;
			pre_distance = -1;

			query1_startRpos = sig1->ref_pos;
			query1_svlen = sig1->cigar_op_len;
			query1_endRpos = getEndRefPosAlnSeg(query1_startRpos, query_opflag, query1_svlen);
			query1_startQpos = sig1->query_pos;
			query1_endQpos = getEndQueryPosAlnSeg(query1_startQpos, query_opflag, query1_svlen);

			for(j=i+1; j<sig_vec.size(); j++){
				sig2 = sig_vec.at(j);
				// skip normal signature
				if(sig2->cigar_op!=BAM_CINS and sig2->cigar_op!=BAM_CDEL) continue;
				merging_seq = "";
				// compute distance
				query2_startRpos = sig2->ref_pos;
				// if(pre_refpos_comp < 0) distance = query2_startRpos - query1_endRpos; // -1;
				// else distance = query2_startRpos - pre_refpos_comp; // -1;
				dist = query2_startRpos - query1_endRpos;

				if(dist > max_ref_dist_thres) break;
				if(query_opflag == sig2->cigar_op and dist <= max_ref_dist_thres){
					query2_svlen = sig2->cigar_op_len;
					query2_startQpos = sig2->query_pos;
					query2_endQpos = getEndQueryPosAlnSeg(query2_startQpos, sig2->cigar_op, query2_svlen);
					query2_endRpos = getEndRefPosAlnSeg(query2_startRpos, sig2->cigar_op, query2_svlen);
					if(query_opflag==BAM_CDEL){ // DEL
						// get compared query pos
						// start_querypos_comp = query1_endQpos;	// deleted on 2025-01-03
						// end_querypos_comp = query2_startQpos - 1;	// deleted on 2025-01-03
						start_querypos_comp = query2_startQpos - dist;
						end_querypos_comp = query2_startQpos - 1;
						// get compared query seq
						comp_queryseq = "";
						for(l=start_querypos_comp-1; l<end_querypos_comp; l++) comp_queryseq += "=ACMGRSVTWYHKDBN"[bam_seqi(seq_int, l)];  // seq

						// get compared ref pos
						// start_refpos_comp = query1_endRpos + query2_svlen;	// deleted on 2025-01-03
						// end_refpos_comp = start_refpos_comp + distance;	// deleted on 2025-01-03
						start_refpos_comp = query2_endRpos - dist;
						end_refpos_comp = query2_endRpos - 1;

						// get compared ref seq
						reg_str = query_seq_info_node->clip_aln->chrname + ":" + to_string(start_refpos_comp) + "-" + to_string(end_refpos_comp);
						pthread_mutex_lock(&mutex_fai);
						seq = fai_fetch(fai, reg_str.c_str(), &seq_len);
						pthread_mutex_unlock(&mutex_fai);
						comp_refseq = seq;
						free(seq);

						// cout << "sig1_op_len=" << sig1->cigar_op_len << "_" << sig1->cigar_op << ", sig2_op_len=" << sig2->cigar_op_len << "_" << sig2->cigar_op << ", distance=" << distance << endl;
						// cout << "comp_queryseq=" << comp_queryseq << ", comp_refseq=" << comp_refseq << endl;
						iden_val = computeVarseqIdentity(comp_queryseq, comp_refseq);
						// cout << "iden_val=" << iden_val << endl;

						if(pre_distance>0) min_ref_dist_thres += pre_distance;
						// cout << "min_ref_dist_thres=" << min_ref_dist_thres << endl;

						if(((double)dist/(sig1->cigar_op_len+sig2->cigar_op_len+dist)<min_valid_sig_size_ratio_thres and dist<min_ref_dist_thres) or iden_val>=min_merge_identity_thres){
							flag = true;
							sig1->cigar_op_len += query2_svlen;

							query1_startRpos = sig1->ref_pos;
							query1_svlen = sig1->cigar_op_len;
							query1_endRpos = getEndRefPosAlnSeg(query1_startRpos, query_opflag, query1_svlen);

							sig1->altseq += sig2->altseq;
							// pre_refpos_comp = query2_endRpos;
							pre_distance = dist;
							delete sig2;
							sig_vec.erase(sig_vec.begin()+j);
							j--;
						}else
							break;
					}else if(query_opflag==BAM_CINS){ // INS
						// get compared query pos
						start_querypos_comp = query2_endQpos - dist;
						end_querypos_comp = query2_endQpos;
						// get compared query seq
						comp_queryseq = "";
						for(l=start_querypos_comp; l<end_querypos_comp; l++) comp_queryseq += "=ACMGRSVTWYHKDBN"[bam_seqi(seq_int, l)];  // seq

						// get compared ref pos
						start_refpos_comp = query2_startRpos - dist;
						end_refpos_comp = query2_startRpos - 1;
						// get compared ref seq
						reg_str = query_seq_info_node->clip_aln->chrname + ":" + to_string(start_refpos_comp) + "-" + to_string(end_refpos_comp);
						pthread_mutex_lock(&mutex_fai);
						seq = fai_fetch(fai, reg_str.c_str(), &seq_len);
						pthread_mutex_unlock(&mutex_fai);
						comp_refseq = seq;
						free(seq);

						// cout << "sig1_op_len=" << sig1->cigar_op_len << "_" << sig1->cigar_op << ", sig2_op_len=" << sig2->cigar_op_len << "_" << sig2->cigar_op << ", distance=" << distance << endl;
						// cout << "comp_queryseq=" << comp_queryseq << ", comp_refseq=" << comp_refseq << endl;
						iden_val = computeVarseqIdentity(comp_queryseq, comp_refseq);
						// cout << "iden_val=" << iden_val << endl;

						if(pre_distance>0)	min_ref_dist_thres += pre_distance;
						// cout << "min_ref_dist_thres=" << min_ref_dist_thres << endl;

						if(((double)dist/(sig1->cigar_op_len+sig2->cigar_op_len+dist)<min_valid_sig_size_ratio_thres and dist < min_ref_dist_thres) or iden_val >= min_merge_identity_thres){

							query1_startQpos = sig1->query_pos;
							query1_endQpos = getEndQueryPosAlnSeg(query1_startQpos, query_opflag, query1_svlen);

							// get altseq after merge
							start_querypos_merged = query1_endQpos;
							end_querypos_merged = query1_endQpos + query2_svlen;
							for(l=start_querypos_merged-1; l<end_querypos_merged; l++) merging_seq += "=ACMGRSVTWYHKDBN"[bam_seqi(seq_int, l)];
							sig1->altseq += merging_seq;

							sig1->cigar_op_len += query2_svlen;
							query1_svlen = sig1->cigar_op_len;

							pre_distance = dist;
							delete sig2;
							sig_vec.erase(sig_vec.begin()+j);
							j--;
							flag = true;
						}else
							break;
					}
				}else
					break;
			}
		}
		// if(sig_vec.size()>0 and flag){
		// 	cout << "~~~~~~~~~~~~" << endl;
		// 	for (i=0; i<sig_vec.size(); i++){
		// 		// get reads seq based on position
		// 		start_querypos_comp = sig_vec.at(i)->query_pos;
		// 		end_querypos_comp = getEndQueryPosAlnSeg(start_querypos_comp, sig_vec.at(i)->cigar_op, sig_vec.at(i)->cigar_op_len);
		// 		// get compared query seq
		// 		comp_queryseq = "";
		// 		for(l=start_querypos_comp-1; l<end_querypos_comp; l++) comp_queryseq += "=ACMGRSVTWYHKDBN"[bam_seqi(seq_int, l)];  // seq
		// 		iden_val = computeVarseqIdentity(comp_queryseq, sig_vec.at(i)->altseq);
		// 		cout << "iden=" << iden_val << ", comp_queryseq=" << comp_queryseq << ", sig_vec.at("<< i << ")->altseq=" << sig_vec.at(i)->altseq << endl;
		// 	}
		// }
		sig_vec.shrink_to_fit();
	}
	// if(flag)
	// 	cout << "merge_flag=" << flag << " sig_vec.size()=" << sig_vec.size() << endl;
	return flag;
}

// filter small signatures (l/max_size<0.2)
void indelRegCluster::filterSmallSigs(vector<qcSig_t *> &qcSig_vec, double valid_size_ratio_thres){
	size_t i;
	int32_t max_size_ins, max_size_del;
	qcSig_t *sig;
	bool valid_flag;

	max_size_ins = max_size_del = 0;
	for(i=0; i<qcSig_vec.size(); i++){
		sig = qcSig_vec.at(i);
		if(sig->cigar_op==BAM_CINS and max_size_ins<sig->cigar_op_len) max_size_ins = sig->cigar_op_len;
		else if(sig->cigar_op==BAM_CDEL and max_size_del<sig->cigar_op_len) max_size_del = sig->cigar_op_len;
	}

	if(max_size_ins>0 or max_size_del>0){
		for(i=0; i<qcSig_vec.size(); ){
			sig = qcSig_vec.at(i);
			valid_flag = true;
			if(sig->cigar_op==BAM_CINS and (double)sig->cigar_op_len/max_size_ins<valid_size_ratio_thres) valid_flag = false;
			else if(sig->cigar_op==BAM_CDEL and (double)sig->cigar_op_len/max_size_del<valid_size_ratio_thres) valid_flag = false;
			if(valid_flag==false){
				delete sig;
				qcSig_vec.erase(qcSig_vec.begin()+i);
			}else i++;
		}
	}
}

// compute the sum of insertions and deletions in the
void indelRegCluster::computeVarSums(struct querySeqInfoNode* &query_seq_info_node, int64_t start_var_pos, int64_t end_var_pos){
	size_t i;
	int64_t ins_sum, del_sum, end_ref_pos;
	qcSig_t *sig;
	bool flag;

	ins_sum = del_sum = 0;
	for(i=0; i<query_seq_info_node->qcSig_vec.size(); i++){
		sig = query_seq_info_node->qcSig_vec.at(i);
		end_ref_pos = getEndRefPosAlnSeg(sig->ref_pos, sig->cigar_op, sig->cigar_op_len);
		flag = isOverlappedPos(sig->ref_pos, end_ref_pos, start_var_pos, end_var_pos);
		if(flag){
			if(sig->cigar_op==BAM_CINS) ins_sum += sig->cigar_op_len;
			else if(sig->cigar_op==BAM_CDEL) del_sum += sig->cigar_op_len;
		}
	}
	query_seq_info_node->ins_sum = ins_sum;
	query_seq_info_node->del_sum = del_sum;
}

// prepare query information for cluster
void indelRegCluster::prepareQueryInfoForCluster(vector<struct querySeqInfoNode*> &query_seq_info_vec){
	int32_t span, span_next;
	struct querySeqInfoNode *q_seq_node;
	vector<struct alnSeg*> query_alnSegs;
	struct alnSeg *seg;
	size_t i, j;
	int64_t end_ref_pos, start_var_pos, end_var_pos, sum;
	bool flag;
	qcSig_t *sig;

	//sort in descending order by |query_endRpos - query_startRpos|
	//query_seq_info_item = new struct querySeqInfoNode();
	for(i=0; i<query_seq_info_vec.size(); i++){
		for(j=0; j<query_seq_info_vec.size()-i-1; j++){
			span = query_seq_info_vec.at(j)->query_alnSegs.at(query_seq_info_vec.at(j)->query_alnSegs.size() - 1)->startRpos +  query_seq_info_vec.at(j)->query_alnSegs.at(query_seq_info_vec.at(j)->query_alnSegs.size() - 1)->seglen - query_seq_info_vec.at(j)->query_alnSegs.at(0)->startRpos;
			span_next = query_seq_info_vec.at(j+1)->query_alnSegs.at(query_seq_info_vec.at(j+1)->query_alnSegs.size() - 1)->startRpos +  query_seq_info_vec.at(j+1)->query_alnSegs.at(query_seq_info_vec.at(j+1)->query_alnSegs.size() - 1)->seglen - query_seq_info_vec.at(j+1)->query_alnSegs.at(0)->startRpos;
			if(span<span_next){ // swap
				q_seq_node = query_seq_info_vec.at(j);
				query_seq_info_vec.at(j) = query_seq_info_vec.at(j+1);
				query_seq_info_vec.at(j+1) = q_seq_node;
			}
		}
	}

	// count the number of overlapped signatures
//	start_var_pos = var_startRefPos - EXT_SIZE_CHK_VAR_LOC; // deleted on 2025-02-28
//	end_var_pos = var_endRefPos + EXT_SIZE_CHK_VAR_LOC;
	start_var_pos = var_startRefPos - max_merge_span;
	end_var_pos = var_endRefPos + max_merge_span;
	for(i=0; i<query_seq_info_vec.size(); i++){
		q_seq_node = query_seq_info_vec.at(i);
		q_seq_node->overlap_sig_num = 0;
		sum = 0; // added on 2025-01-31
		for(j=0; j<q_seq_node->qcSig_vec.size(); j++){
			sig = q_seq_node->qcSig_vec.at(j);
			if(sig->cigar_op_len>=min_sv_size and (sig->cigar_op==BAM_CINS or sig->cigar_op==BAM_CDEL)){
				end_ref_pos = getEndRefPosAlnSeg(sig->ref_pos, sig->cigar_op, sig->cigar_op_len);
				flag = isOverlappedPos(sig->ref_pos, end_ref_pos, start_var_pos, end_var_pos);
//				for(k=0; k<varVec.size(); k++){
//					reg = varVec.at(k);
//					if(isOverlappedPos(sig->ref_pos, end_ref_pos, reg->startRefPos, reg->endRefPos)){
//						flag = true;
//						break;
//					}
//				}
				if(flag){
					//cout << "i=" << i << ", j=" << j << ", startpos=" << seg->startRpos << ", endpos=" << end_ref_pos << ", op=" << seg->opflag << endl;
					q_seq_node->overlap_sig_num ++;
					sum += sig->cigar_op_len;  // added on 2025-01-05
				}
			}
		}
		if(q_seq_node->overlap_sig_num>0) q_seq_node->aver_sig_size = (double)sum / q_seq_node->overlap_sig_num;
		else q_seq_node->aver_sig_size = 0;
	}

	// sort descendingly according to the number of overlapped signatures
	for(i=0; i<query_seq_info_vec.size(); i++){
		for(j=i+1; j<query_seq_info_vec.size(); j++){
			if(query_seq_info_vec.at(i)->overlap_sig_num<query_seq_info_vec.at(j)->overlap_sig_num){
				q_seq_node = query_seq_info_vec.at(i);
				query_seq_info_vec.at(i) = query_seq_info_vec.at(j);
				query_seq_info_vec.at(j) = q_seq_node;
			}
		}
	}

	// sort according average length of overlapped signatures
	sortQueryInfoByAverSigSize(query_seq_info_vec);

	// sort according to the category of the number of signatures
	sortQueryInfoByNumCategory(query_seq_info_vec);

	// calculate the entire flanking flag
//	start_var_pos = var_startRefPos - EXT_SIZE_ENTIRE_FLANKING;
//	end_var_pos = var_endRefPos + EXT_SIZE_ENTIRE_FLANKING;  // deleted on 2025-03-02
	start_var_pos = var_startRefPos - EXT_SIZE_CHK_VAR_LOC;
	end_var_pos = var_endRefPos + EXT_SIZE_CHK_VAR_LOC;
	if(start_var_pos<1) start_var_pos = 1;
	if(end_var_pos>chrlen) end_var_pos = chrlen;
	for(i=0; i<query_seq_info_vec.size(); i++){
		query_alnSegs = query_seq_info_vec.at(i)->query_alnSegs;
		seg = query_alnSegs.at(query_alnSegs.size()-1);
		end_ref_pos = getEndRefPosAlnSeg(seg->startRpos, seg->opflag, seg->seglen);
		if(query_alnSegs.at(0)->startRpos<=start_var_pos and end_ref_pos>=end_var_pos) query_seq_info_vec.at(i)->entire_flanking_flag = true;
		else query_seq_info_vec.at(i)->entire_flanking_flag = false;
	}

	// prefer the alignments flanking the entire variant region
	for(i=0; i<query_seq_info_vec.size(); i++){
		if(query_seq_info_vec.at(i)->entire_flanking_flag==false){
			flag = false;
			for(j=i+1; j<query_seq_info_vec.size(); j++){
				if(query_seq_info_vec.at(j)->entire_flanking_flag){
					flag = true;
					q_seq_node = query_seq_info_vec.at(i);
					query_seq_info_vec.at(i) = query_seq_info_vec.at(j);
					query_seq_info_vec.at(j) = q_seq_node;
					break;
				}
			}
			if(flag==false) break;
		}
	}

}

// sort according to average length of overlapped signatures
void indelRegCluster::sortQueryInfoByAverSigSize(vector<struct querySeqInfoNode*> &query_seq_info_vec){
	size_t i, j;
	struct querySeqInfoNode *node_tmp;

	// sort descendingly according to the average size of overlapped signatures
	for(i=0; i<query_seq_info_vec.size(); i++){
		for(j=i+1; j<query_seq_info_vec.size(); j++){
			if(query_seq_info_vec.at(i)->aver_sig_size<query_seq_info_vec.at(j)->aver_sig_size){
				node_tmp = query_seq_info_vec.at(i);
				query_seq_info_vec.at(i) = query_seq_info_vec.at(j);
				query_seq_info_vec.at(j) = node_tmp;
			}
		}
	}
}

// sort according to number category of overlapped signatures
void indelRegCluster::sortQueryInfoByNumCategory(vector<struct querySeqInfoNode*> &query_seq_info_vec){
	size_t i, j, sum;
	vector<struct querySeqInfoNode*> query_seq_info_vec_tmp;
	int32_t idx;
	vector<mrmin_match_t*> num_id_vec;
	mrmin_match_t *num_id_node;
	bool flag;

	for(i=0; i<query_seq_info_vec.size(); i++){
		flag = false;
		for(j=0; j<num_id_vec.size(); j++){
			if(num_id_vec.at(j)->num==(int32_t)query_seq_info_vec.at(i)->overlap_sig_num){
				idx = j;
				flag = true;
				break;
			}
		}
		if(flag){ // exist, update item
			num_id_node  = num_id_vec.at(idx);
			num_id_node->pos_id.push_back(i);
		}else{ // add item
			num_id_node = new mrmin_match_t();
			num_id_node->num = query_seq_info_vec.at(i)->overlap_sig_num;
			num_id_node->pos_id.push_back(i);
			num_id_node->MR_recluster = -1;
			num_id_node->candidate_flag = false;
			num_id_vec.push_back(num_id_node);
		}
	}

	// sort descendingly according to the number category of overlapped signatures
	for(i=0; i<num_id_vec.size(); i++){
		if(num_id_vec.at(i)->num>0){
			for(j=i+1; j<num_id_vec.size(); j++){
				if(num_id_vec.at(j)->num>0){
					if(num_id_vec.at(i)->pos_id.size()<num_id_vec.at(j)->pos_id.size()){
						num_id_node = num_id_vec.at(i);
						num_id_vec.at(i) = num_id_vec.at(j);
						num_id_vec.at(j) = num_id_node;
					}
				}
			}
		}
	}

	// update the query vector
	sum = 0;
	for(i=0; i<num_id_vec.size(); i++) sum += num_id_vec.at(i)->pos_id.size();
	if(sum==query_seq_info_vec.size()){
		for(i=0; i<num_id_vec.size(); i++){
			num_id_node = num_id_vec.at(i);
			for(j=0; j<(size_t)num_id_node->pos_id.size(); j++) query_seq_info_vec_tmp.push_back(query_seq_info_vec.at(num_id_node->pos_id.at(j)));
		}
		for(i=0; i<query_seq_info_vec_tmp.size(); i++) query_seq_info_vec.at(i) = query_seq_info_vec_tmp.at(i);
	}else{
		cerr << "sum=" << sum << ", query_seq_info_vec.size=" << query_seq_info_vec.size() << ", error!" << endl;
		exit(1);

	}

	// release memory
	for(i=0; i<num_id_vec.size(); i++) delete num_id_vec.at(i);
	vector<mrmin_match_t*>().swap(num_id_vec);
}

// estimate number of the group
int32_t indelRegCluster::estimateClusterGroupNum(vector<struct querySeqInfoNode*> &query_seq_info_vec){
	int32_t max_group_num, empty_num, item_num;
	size_t i;
	double empty_ratio;
	struct querySeqInfoNode *q_node;

	empty_num = item_num = 0;
	// for(i=0; i<query_seq_info_vec.size(); i++) if(query_seq_info_vec.at(i)->qcSig_vec.empty()) empty_num ++; // deleted on 2025-02-24
	//for(i=0; i<query_seq_info_vec.size(); i++) if(query_seq_info_vec.at(i)->ins_sum == 0 and query_seq_info_vec.at(i)->del_sum == 0) empty_num ++; // deleted on 2025-03-02
	for(i=0; i<query_seq_info_vec.size(); i++){
		q_node = query_seq_info_vec.at(i);
		if(q_node->entire_flanking_flag==true){
			if(q_node->ins_sum==0 and q_node->del_sum==0) empty_num ++;
			item_num ++;
		}
	}

	//empty_ratio = (double)empty_num / query_seq_info_vec.size();
	if(item_num>0) empty_ratio = (double)empty_num / item_num;
	else empty_ratio = 0;
	if(empty_ratio>=MIN_SINGLE_CLUSTER_GROUP_RATIO_THRES) max_group_num = 1;
	else max_group_num = 2;

#if POA_ALIGN_DEBUG
	cout << "max_group_num=" << max_group_num << ", empty_num=" << empty_num << ", empty_ratio=" << empty_ratio << ", valid_item_num=" << item_num << ", vec_size=" << query_seq_info_vec.size() << endl;
#endif

	return max_group_num;
}

vector<struct querySeqInfoVec*> indelRegCluster::queryClusterSingleGroup(vector<struct querySeqInfoNode*> &query_seq_info_all){
	vector<struct querySeqInfoVec*> query_clu_vec;
	size_t i;
	struct querySeqInfoVec *q_cluster_node;
	struct querySeqInfoNode *qnode;
	vector<struct querySeqInfoNode*> q_cluster_a;

	int32_t ins_sum_num, del_sum_num, target_var_op, item_num;
	double ins_sum_mean0, del_sum_mean0, ins_sum_mean, del_sum_mean, mean_size, sdev, dif;

	// compute mean size
	ins_sum_mean0 = del_sum_mean0 = ins_sum_num = del_sum_num = 0;
	for(i=0; i<query_seq_info_all.size(); i++){
		qnode = query_seq_info_all.at(i);
		if(qnode->overlap_sig_num==0) continue; // skip items of no signatures

		if(qnode->ins_sum>0){
			ins_sum_mean0 += qnode->ins_sum;
			ins_sum_num ++;
		}
		if(qnode->del_sum>0){
			del_sum_mean0 += qnode->del_sum;
			del_sum_num ++;
		}
	}
	if(ins_sum_num>0) ins_sum_mean0 /= ins_sum_num;
	if(del_sum_num>0) del_sum_mean0 /= del_sum_num;

	// skip outliers
	ins_sum_mean = del_sum_mean = ins_sum_num = del_sum_num = 0;
	for(i=0; i<query_seq_info_all.size(); i++){
		qnode = query_seq_info_all.at(i);
		if(qnode->overlap_sig_num==0) continue; // skip items of no signatures

		if(qnode->ins_sum>0 and qnode->ins_sum<3*ins_sum_mean0){
			ins_sum_mean += qnode->ins_sum;
			ins_sum_num ++;
		}
		if(qnode->del_sum>0 and qnode->del_sum<3*del_sum_mean0){
			del_sum_mean += qnode->del_sum;
			del_sum_num ++;
		}
	}
	if(ins_sum_num>0) ins_sum_mean /= ins_sum_num;
	if(del_sum_num>0) del_sum_mean /= del_sum_num;

	if(ins_sum_mean>del_sum_mean) { target_var_op = BAM_CINS; mean_size = ins_sum_mean; }
	else { target_var_op = BAM_CDEL; mean_size = del_sum_mean; }

#if POA_ALIGN_DEBUG
	cout << "target_var_op=" << target_var_op << ", mean_size=" << mean_size << ", ins_sum_mean=" << ins_sum_mean << ", ins_sum_num=" << ins_sum_num << ", del_sum_mean=" << del_sum_mean << ", del_sum_num=" << del_sum_num << endl;
#endif

	// compute sdev
	if(mean_size>0){
		item_num = sdev = 0;
		for(i=0; i<query_seq_info_all.size(); i++){
			qnode = query_seq_info_all.at(i);
			if(qnode->overlap_sig_num==0) continue; // skip items of no signatures

			if(target_var_op==BAM_CINS and qnode->ins_sum>0){
				sdev += (qnode->ins_sum - mean_size) * (qnode->ins_sum - mean_size);
				item_num ++;
			}else if(target_var_op==BAM_CDEL and qnode->del_sum>0){
				sdev += (qnode->del_sum - mean_size) * (qnode->del_sum - mean_size);
				item_num ++;
			}
		}
		sdev = sqrt(sdev/item_num);

#if POA_ALIGN_DEBUG
		cout << "target_var_op=" << target_var_op << ", mean_size=" << mean_size << ", sdev=" << sdev << ", item_num=" << item_num << endl;
#endif

		// remove the outliers
		if(item_num>0){
			for(i=0; i<query_seq_info_all.size(); i++){
				qnode = query_seq_info_all.at(i);
				if(qnode->overlap_sig_num==0) continue; // skip items of no signatures

				if(target_var_op==BAM_CINS and qnode->ins_sum>0){
					dif =  abs(qnode->ins_sum - mean_size);
#if POA_ALIGN_DEBUG
					cout << "\tdif=" << dif << endl;
#endif
					if(dif<=2*sdev){
						q_cluster_a.push_back(qnode);
					}
#if POA_ALIGN_DEBUG
					else{
						cout << "\txxxxxxxxx dif=" << dif << ", qname=" << qnode->qname << endl;
					}
#endif
				}else if(target_var_op==BAM_CDEL and qnode->del_sum>0){
					dif =  abs(qnode->del_sum - mean_size);
#if POA_ALIGN_DEBUG
					cout << "\tdif=" << dif << endl;
#endif
					if(dif<2*sdev){
						q_cluster_a.push_back(qnode);
					}
#if POA_ALIGN_DEBUG
					else{
						cout << "\txxxxxxxxx dif=" << dif << ", qname=" << qnode->qname << endl;
					}
#endif
				}
			}
		}

		if(q_cluster_a.size()>0){
			q_cluster_node = new struct querySeqInfoVec();
			q_cluster_node->query_seq_info_all = q_cluster_a;
			q_cluster_node->rescue_flag = false;
			query_clu_vec.push_back(q_cluster_node);
		}
	}

	return query_clu_vec;
}

void indelRegCluster::queryClusterAppendSingleGroup(vector<struct querySeqInfoNode*> &q_cluster_a, vector<struct querySeqInfoNode*> &query_seq_info_all){
	size_t i;
	struct querySeqInfoNode *qnode;
	int32_t ins_sum_num, del_sum_num, target_var_op, item_num;
	double ins_sum_mean0, del_sum_mean0, ins_sum_mean, del_sum_mean, mean_size, sdev, dif;

	// compute mean size
	ins_sum_mean0 = del_sum_mean0 = ins_sum_num = del_sum_num = 0;
	for(i=0; i<q_cluster_a.size(); i++){
		qnode = q_cluster_a.at(i);
		if(qnode->overlap_sig_num==0) continue; // skip items of no signatures

		if(qnode->ins_sum>0){
			ins_sum_mean0 += qnode->ins_sum;
			ins_sum_num ++;
		}
		if(qnode->del_sum>0){
			del_sum_mean0 += qnode->del_sum;
			del_sum_num ++;
		}
	}
	if(ins_sum_num>0) ins_sum_mean0 /= ins_sum_num;
	if(del_sum_num>0) del_sum_mean0 /= del_sum_num;

	// skip outliers
	ins_sum_mean = del_sum_mean = ins_sum_num = del_sum_num = 0;
	for(i=0; i<q_cluster_a.size(); i++){
		qnode = q_cluster_a.at(i);
		if(qnode->overlap_sig_num==0) continue; // skip items of no signatures

		if(qnode->ins_sum>0 and qnode->ins_sum<3*ins_sum_mean0){
			ins_sum_mean += qnode->ins_sum;
			ins_sum_num ++;
		}
		if(qnode->del_sum>0 and qnode->del_sum<3*del_sum_mean0){
			del_sum_mean += qnode->del_sum;
			del_sum_num ++;
		}
	}
	if(ins_sum_num>0) ins_sum_mean /= ins_sum_num;
	if(del_sum_num>0) del_sum_mean /= del_sum_num;

	if(ins_sum_mean>del_sum_mean) { target_var_op = BAM_CINS; mean_size = ins_sum_mean; }
	else { target_var_op = BAM_CDEL; mean_size = del_sum_mean; }

#if POA_ALIGN_DEBUG
	cout << "target_var_op=" << target_var_op << ", mean_size=" << mean_size << ", ins_sum_mean=" << ins_sum_mean << ", ins_sum_num=" << ins_sum_num << ", del_sum_mean=" << del_sum_mean << ", del_sum_num=" << del_sum_num << endl;
#endif

	// compute sdev
	if(mean_size>0){
		item_num = sdev = 0;
		for(i=0; i<q_cluster_a.size(); i++){
			qnode = q_cluster_a.at(i);
			if(qnode->overlap_sig_num==0) continue; // skip items of no signatures

			if(target_var_op==BAM_CINS and qnode->ins_sum>0){
				sdev += (qnode->ins_sum - mean_size) * (qnode->ins_sum - mean_size);
				item_num ++;
			}else if(target_var_op==BAM_CDEL and qnode->del_sum>0){
				sdev += (qnode->del_sum - mean_size) * (qnode->del_sum - mean_size);
				item_num ++;
			}
		}
		sdev = sqrt(sdev/item_num);

#if POA_ALIGN_DEBUG
		cout << "target_var_op=" << target_var_op << ", mean_size=" << mean_size << ", sdev=" << sdev << ", item_num=" << item_num << endl;
#endif

		// remove the outliers
		if(item_num>0){
			for(i=0; i<query_seq_info_all.size(); i++){
				qnode = query_seq_info_all.at(i);
				if(qnode->cluster_finished_flag or qnode->overlap_sig_num==0 or qnode->qcSig_vec.size()==0) continue; // skip items of no signatures

				if(target_var_op==BAM_CINS and qnode->ins_sum>0){
					dif =  abs(qnode->ins_sum - mean_size);
#if POA_ALIGN_DEBUG
					cout << "\tdif=" << dif << endl;
#endif
					if(dif<=2*sdev){
						q_cluster_a.push_back(qnode);
					}
#if POA_ALIGN_DEBUG
					else{
						cout << "\txxxxxxxxx dif=" << dif << ", qname=" << qnode->qname << endl;
					}
#endif
				}else if(target_var_op==BAM_CDEL and qnode->del_sum>0){
					dif =  abs(qnode->del_sum - mean_size);
#if POA_ALIGN_DEBUG
					cout << "\tdif=" << dif << endl;
#endif
					if(dif<2*sdev){
						q_cluster_a.push_back(qnode);
					}
#if POA_ALIGN_DEBUG
					else{
						cout << "\txxxxxxxxx dif=" << dif << ", qname=" << qnode->qname << endl;
					}
#endif
				}
			}
		}
	}
}

vector<struct querySeqInfoVec*> indelRegCluster::queryClusterDoubleGroup(vector<struct querySeqInfoNode*> &query_seq_info_all){
	vector<struct querySeqInfoVec*> query_clu_vec;
	struct querySeqInfoVec *q_cluster_node;
	vector<struct querySeqInfoNode*> q_cluster_a, q_cluster_b, q_cluster_rescue, q_cluster_a_final, q_cluster_b_final, q_recluster_a;
	// vector<indelCluster_t*> q_cluster_a;
	struct seedQueryInfo *seed_info_a, *seed_info_b;
	size_t i, j, k;
	int32_t smin_id, mr_id_a, mr_id_b, qseq_id_a, qseq_id_b;
	double match_ratio_a, match_ratio_b, sco_min;
	vector<int32_t> min_id_vec, min_id_vec_recluster;
	vector<double> match_ratio_vec;
	vector<mrmin_match_t*> mrmin_match_vec, mrmin_match_vec_recluster;
	mrmin_match_t *mrmin_match_node;
	bool flag, q_cluster_b_null_valid_flag;
	bool recluster_flag, mrmin_match_vec_begin_flag;
	vector<int32_t> init_idx_vec;

//	skip_num = round(query_seq_info_all.size()*0.1);
//	if(query_seq_info_all.at(skip_num)->qcSig_vec.size()==0) skip_num = 0;

	// select initial clusters
	for(k=0; k<query_seq_info_all.size(); k++){
		// validate the signature
		if(query_seq_info_all.at(k)->overlap_sig_num==0) continue; // skip items of no signatures

		//initializing a_cluster
		query_seq_info_all.at(k)->cluster_finished_flag = true;
		q_cluster_a.push_back(query_seq_info_all.at(k));

		//initializing b_cluster: choose one that is most different query from a_cluster add in b_cluster
		sco_min = INT_MAX; //0.99
		match_ratio_vec.clear();
		min_id_vec.clear();
		q_cluster_b_null_valid_flag = false;

		for(i=0; i<query_seq_info_all.size(); i++){

			// if(query_seq_info_all.at(k)->qname.compare("SRR11008518.1.822753")==0){
			// 	cout << query_seq_info_all.at(k)->qname << endl;
			// }

			if(query_seq_info_all.at(i)->overlap_sig_num==0 or query_seq_info_all.at(i)->entire_flanking_flag==false) continue; // skip items of no signatures

//			if(query_seq_info_all.at(i)->cluster_finished_flag == false){
				seed_info_a = chooseSeedClusterQuery(query_seq_info_all.at(i), q_cluster_a);
				if(seed_info_a->span_intersection>0){
					if(i!=k){
						match_ratio_a = computeMatchRatio(query_seq_info_all.at(i), q_cluster_a.at(seed_info_a->id), seed_info_a->startSpanPos, seed_info_a->endSpanPos, min_sv_size, min_identity_match);
						match_ratio_vec.push_back(match_ratio_a);
					}else{
						match_ratio_a = 1;
						match_ratio_vec.push_back(match_ratio_a);
					}
					min_id_vec.push_back(i);
					if(match_ratio_a<sco_min) sco_min = match_ratio_a;

//					cout << "query.qname=[" << query_seq_info_all.at(i)->qname << "], q_cluster_a.qname=[" << q_cluster_a.at(seed_info_a->id)->qname << "], score_ratio_a=" << score_ratio_a << ", sco_min=" << sco_min << endl;

				}else{
					match_ratio_vec.push_back(2);
				}
				delete seed_info_a;
//			}else{
//				score_ratio_vec.push_back(1);
//			}
		}
		query_seq_info_all.at(k)->cluster_finished_flag = false;
		q_cluster_a.pop_back();

		//exist_flag = false;
		if(sco_min < 1){
			min_id_vec.clear();
			for(i=0; i<match_ratio_vec.size(); i++){
				if(match_ratio_vec.at(i)<=1)
					min_id_vec.push_back(i);
			}
			for(i=0; i<min_id_vec.size(); i++){
				flag = false;
				for(j=0; j<mrmin_match_vec.size(); j++){
					if(match_ratio_vec.at(min_id_vec.at(i))==mrmin_match_vec.at(j)->MR){
						smin_id = j;
						flag = true;
						break;
					}
				}
				if(flag==false){
					mrmin_match_node = new mrmin_match_t();
					mrmin_match_node->pos_id.push_back(min_id_vec.at(i));
					mrmin_match_node->num = 1;
					mrmin_match_node->MR = match_ratio_vec.at(min_id_vec.at(i));
					mrmin_match_node->MR_recluster = -1;
					mrmin_match_node->candidate_flag = false;
					mrmin_match_vec.push_back(mrmin_match_node);
				}else{
					mrmin_match_vec.at(smin_id)->num += 1;
					mrmin_match_vec.at(smin_id)->pos_id.push_back(min_id_vec.at(i));
				}
			}

			// sort descendingly and initialize clusters adaptively
			sortMRVec(mrmin_match_vec);

#if POA_ALIGN_DEBUG
			printMRVec(mrmin_match_vec);
#endif

			flag = true;
			if(mrmin_match_vec.size()==2 and mrmin_match_vec.at(0)->MR==0 and mrmin_match_vec.at(1)->MR==1 and mrmin_match_vec.at(1)->num==1) flag = false; // skip, and try next one
			else if(mrmin_match_vec.size()>2 and mrmin_match_vec.at(0)->MR==0 and mrmin_match_vec.at(0)->num>MIN_INVLAID_SECOND_GROUP_RATIO_THRES*min_id_vec.size()) flag = false; // skip, and try next one

			if(flag){ // recluster for items of MR=0
				recluster_flag = getReClusterFlag(mrmin_match_vec);
				if(recluster_flag){
					min_id_vec_recluster.clear();
					mrmin_match_vec_begin_flag = initReCluster(mrmin_match_vec, min_id_vec_recluster); //[0] valid [1] begin
					mrmin_match_vec_recluster = reClusterByMR(query_seq_info_all, min_id_vec_recluster, mrmin_match_vec_begin_flag);
					for(i=0; i<mrmin_match_vec_recluster.size(); i++){
						mrmin_match_vec_recluster.at(i)->MR_recluster = mrmin_match_vec_recluster.at(i)->MR;
						mrmin_match_vec_recluster.at(i)->MR = 0;
						mrmin_match_vec.push_back(mrmin_match_vec_recluster.at(i));
					}

					//sort
					sortMRVec(mrmin_match_vec);

#if POA_ALIGN_DEBUG
					printMRVec(mrmin_match_vec);
#endif
				}

				// select the target sub-set
				init_idx_vec = getInitGroupIdx(mrmin_match_vec, query_seq_info_all);

				mr_id_a = mr_id_b = qseq_id_a = qseq_id_b = -1;
				if(init_idx_vec.size()>=2) {
					mr_id_a = init_idx_vec.at(0);
					mr_id_b = init_idx_vec.at(1);
					seed_info_a = getInitClusterItem(mrmin_match_vec.at(mr_id_a)->pos_id, query_seq_info_all);
					seed_info_b = getInitClusterItem(mrmin_match_vec.at(mr_id_b)->pos_id, query_seq_info_all);
					qseq_id_a = seed_info_a->id;
					qseq_id_b = seed_info_b->id;

					match_ratio_a = computeMatchRatio(query_seq_info_all.at(qseq_id_a), query_seq_info_all.at(qseq_id_b), seed_info_a->startSpanPos, seed_info_a->endSpanPos, min_sv_size, CNS_MIN_IDENTITY_CLUSTER_THRES);

					//cout << "match_ratio_a=" << match_ratio_a << endl;

					if(match_ratio_a==1){ // same group, then only select one group
						qseq_id_b = -1;
						if((mrmin_match_vec.at(mr_id_b)->MR!=0 and recluster_flag==false) or (mrmin_match_vec.at(mr_id_b)->MR_recluster!=0 and recluster_flag==true))
							q_cluster_b_null_valid_flag = true;
					}

					delete seed_info_a;
					delete seed_info_b;
				}else if(init_idx_vec.size()==1){
					mr_id_a = init_idx_vec.at(0);
					seed_info_a = getInitClusterItem(mrmin_match_vec.at(mr_id_a)->pos_id, query_seq_info_all);
					qseq_id_a = seed_info_a->id;
					delete seed_info_a;
					q_cluster_b_null_valid_flag = true;
				}

//				query_seq_info_all.at(qseq_id_a)->cluster_finished_flag = true;
//				q_cluster_a.push_back(query_seq_info_all.at(qseq_id_a));
				if(qseq_id_b!=-1){
					if(init_idx_vec.size()>1 and mrmin_match_vec.at(mr_id_b)->num>=0.6*min_supp_num){
						if((double)mrmin_match_vec.at(mr_id_b)->num/mrmin_match_vec.at(mr_id_a)->num<0.2)
							q_cluster_b_null_valid_flag = true;
//						else{ // deleted on 2025-03-02
//							query_seq_info_all.at(qseq_id_b)->cluster_finished_flag = true;
//							q_cluster_b.push_back(query_seq_info_all.at(qseq_id_b));
//						}
					}else{
						if(init_idx_vec.size()==1){
							if(mrmin_match_vec.at(mr_id_a)->num>=0.6*min_supp_num){
								q_cluster_b_null_valid_flag = true;
							}
						}else if((double)mrmin_match_vec.at(mr_id_b)->num/mrmin_match_vec.at(mr_id_a)->num<0.2){
							q_cluster_b_null_valid_flag = true;
						}
					}
				}

				// copy items to initial clusters
				if(mr_id_a!=-1){
					mrmin_match_node = mrmin_match_vec.at(mr_id_a);
					if(mrmin_match_node->MR==0 and mrmin_match_node->MR_recluster==0){
						// choose seed item with recaluated MR=1 for (MR=0 and MR_recluster=0)
						q_cluster_a = getInitClusterItemsMR0(mrmin_match_vec.at(mr_id_a)->pos_id, query_seq_info_all);
					}else{
						for(i=0; i<(size_t)mrmin_match_node->num; i++){
							qseq_id_a = mrmin_match_node->pos_id.at(i);
							query_seq_info_all.at(qseq_id_a)->cluster_finished_flag = true;
							q_cluster_a.push_back(query_seq_info_all.at(qseq_id_a));
						}
					}
				}
				if(mr_id_b!=-1){
					mrmin_match_node = mrmin_match_vec.at(mr_id_b);
					if(mrmin_match_node->MR==0 and mrmin_match_node->MR_recluster==0){
						// choose seed item with recaluated MR=1 for (MR=0 and MR_recluster=0)
						q_cluster_b = getInitClusterItemsMR0(mrmin_match_vec.at(mr_id_b)->pos_id, query_seq_info_all);
					}else{
						for(i=0; i<(size_t)mrmin_match_node->num; i++){
							qseq_id_b = mrmin_match_node->pos_id.at(i);
							query_seq_info_all.at(qseq_id_b)->cluster_finished_flag = true;
							q_cluster_b.push_back(query_seq_info_all.at(qseq_id_b));
						}
					}
				}
			}
		}else{ // single cluster
			// deleted on 2025-04-13
//			query_seq_info_all.at(k)->cluster_finished_flag = true;
//			q_cluster_a.push_back(query_seq_info_all.at(k));
			// fill single cluster
			for(i=0; i<min_id_vec.size(); i++){
				query_seq_info_all.at(min_id_vec.at(i))->cluster_finished_flag = true;
				q_cluster_a.push_back(query_seq_info_all.at(min_id_vec.at(i)));
			}
			q_cluster_b_null_valid_flag = true;
		}

		if(mrmin_match_vec.size()>0){
			vector<mrmin_match_t*>::iterator svp;
			for(svp=mrmin_match_vec.begin(); svp!=mrmin_match_vec.end(); svp++) delete *svp;
			vector<mrmin_match_t*>().swap(mrmin_match_vec);
		}

		if(q_cluster_b.size()==0 and !q_cluster_b_null_valid_flag){
			if(q_cluster_a.size()>0) q_cluster_a.clear(); // update on 2025-01-12
			query_seq_info_all.at(k)->cluster_finished_flag = false;
		}else{ // two initial sets or one valid set
			break;
		}
	}

#if POA_ALIGN_DEBUG
	// print extracted signatures
	if(q_cluster_a.size()>0){
		cout << "initial cluster A:" << endl;
		printQcSigIndel(q_cluster_a);
	}
	if(q_cluster_b.size()>0){
		cout << "initial cluster B:" << endl;
		printQcSigIndel(q_cluster_b);
	}
#endif

	if(q_cluster_b.size()==0){
		if(q_cluster_a.size()>0)
			queryClusterAppendSingleGroup(q_cluster_a, query_seq_info_all);
		else{
			for(i=0; i<query_seq_info_all.size(); i++){
				//if(query_seq_info_all.at(i)->cluster_finished_flag == false and query_seq_info_all.at(i)->overlap_sig_num>0){ // deleted on 2025-01-30
				if(query_seq_info_all.at(i)->cluster_finished_flag == false and query_seq_info_all.at(i)->overlap_sig_num>0 and query_seq_info_all.at(i)->qcSig_vec.size()>0){
					query_seq_info_all.at(i)->cluster_finished_flag = true;
					q_cluster_a.push_back(query_seq_info_all.at(i));
				}
			}
		}
	}else{
		for(i=0; i<query_seq_info_all.size(); i++){
			if(query_seq_info_all.at(i)->cluster_finished_flag==false and query_seq_info_all.at(i)->overlap_sig_num>0){
				//choose seed query from a_cluster
				seed_info_a = chooseSeedClusterQuery(query_seq_info_all.at(i), q_cluster_a);
				seed_info_b = chooseSeedClusterQuery(query_seq_info_all.at(i), q_cluster_b);
				if(seed_info_a->span_intersection>0){
					match_ratio_a = computeMatchRatio(query_seq_info_all.at(i), q_cluster_a.at(seed_info_a->id), seed_info_a->startSpanPos, seed_info_a->endSpanPos, min_sv_size, min_identity_match);
					if(seed_info_b->span_intersection>0){
						match_ratio_b = computeMatchRatio(query_seq_info_all.at(i), q_cluster_b.at(seed_info_b->id), seed_info_b->startSpanPos, seed_info_b->endSpanPos, min_sv_size, min_identity_match);
						if(match_ratio_a>match_ratio_b){
							query_seq_info_all.at(i)->cluster_finished_flag = true;
							q_cluster_a.push_back(query_seq_info_all.at(i));
						}else{
							if(match_ratio_b==match_ratio_a){
								query_seq_info_all.at(i)->cluster_finished_flag = true;
								q_cluster_rescue.push_back(query_seq_info_all.at(i));
							}else{
								query_seq_info_all.at(i)->cluster_finished_flag = true;
								q_cluster_b.push_back(query_seq_info_all.at(i));
							}
						}
					}else{
						if(match_ratio_a<1){
							query_seq_info_all.at(i)->cluster_finished_flag = true;
							q_cluster_b.push_back(query_seq_info_all.at(i));
						}else{
							query_seq_info_all.at(i)->cluster_finished_flag = true;
							q_cluster_a.push_back(query_seq_info_all.at(i));
						}
					}
				}else{
					if(seed_info_b->span_intersection>0){
						match_ratio_b = computeMatchRatio(query_seq_info_all.at(i), q_cluster_b.at(seed_info_b->id), seed_info_b->startSpanPos, seed_info_b->endSpanPos, min_sv_size, min_identity_match);
						if(match_ratio_b<1){
							query_seq_info_all.at(i)->cluster_finished_flag = true;
							q_cluster_a.push_back(query_seq_info_all.at(i));
						}else{
							query_seq_info_all.at(i)->cluster_finished_flag = true;
							q_cluster_b.push_back(query_seq_info_all.at(i));
						}
					}
				}
				delete seed_info_a;
				delete seed_info_b; // added on 2023-11-23
			}
		}
	}

	if(q_cluster_rescue.size()>0){
		for(i=0; i<q_cluster_rescue.size(); ){
			seed_info_a = chooseSeedClusterQuery(q_cluster_rescue.at(i),q_cluster_a);
			seed_info_b = chooseSeedClusterQuery(q_cluster_rescue.at(i), q_cluster_b);
			match_ratio_a = computeMatchRatio(q_cluster_rescue.at(i), q_cluster_a.at(seed_info_a->id), seed_info_a->startSpanPos, seed_info_a->endSpanPos, min_sv_size, min_identity_match);
			match_ratio_b = computeMatchRatio(q_cluster_rescue.at(i), q_cluster_b.at(seed_info_b->id), seed_info_b->startSpanPos, seed_info_b->endSpanPos, min_sv_size, min_identity_match);
			if(match_ratio_a>match_ratio_b){
				q_cluster_a.push_back(q_cluster_rescue.at(i));
				q_cluster_rescue.erase(q_cluster_rescue.begin()+i);
			}else{
				if(match_ratio_a!=match_ratio_b){
					q_cluster_b.push_back(q_cluster_rescue.at(i));
					q_cluster_rescue.erase(q_cluster_rescue.begin()+i);
				}else{
					query_seq_info_all.at(i)->cluster_finished_flag = false;
					i++;
				}
			}
			// added on 2023-11-23
			delete seed_info_a;
			delete seed_info_b;
		}
	}

	for(i=0; i<q_cluster_a.size(); i++) if(q_cluster_a.at(i)->overlap_sig_num>0) q_cluster_a_final.push_back(q_cluster_a.at(i));
	for(i=0; i<q_cluster_b.size(); i++) if(q_cluster_b.at(i)->overlap_sig_num>0) q_cluster_b_final.push_back(q_cluster_b.at(i));

	if(q_cluster_a_final.size()>0){
		q_cluster_node = new struct querySeqInfoVec();
		q_cluster_node->query_seq_info_all = q_cluster_a_final;
		q_cluster_node->rescue_flag = false;
		query_clu_vec.push_back(q_cluster_node);
	}
	if(q_cluster_b_final.size()>0){
		q_cluster_node = new struct querySeqInfoVec();
		q_cluster_node->query_seq_info_all = q_cluster_b_final;
		q_cluster_node->rescue_flag = false;
		query_clu_vec.push_back(q_cluster_node);
	}
	if(q_cluster_rescue.size()>0){
		q_cluster_node = new struct querySeqInfoVec();
		q_cluster_node->query_seq_info_all = q_cluster_rescue;
		q_cluster_node->rescue_flag = true;
		query_clu_vec.push_back(q_cluster_node);
	}

	return query_clu_vec;
}

// choose the seed query for indel clustering
struct seedQueryInfo* indelRegCluster::chooseSeedClusterQuery(struct querySeqInfoNode* query_seq_info_node, vector<struct querySeqInfoNode*> &q_cluster){
	size_t seed_id, i;
	int64_t startSpanPos_tmp, endSpanPos_tmp, startSpanPos, endSpanPos, span, end_ref_pos1, end_ref_pos2;
	struct seedQueryInfo* seed_info;
	struct querySeqInfoNode *q_node;
	int32_t min_idx, min_dif;
	vector<double> dif_vec;
	struct alnSeg *aln_seg1, *aln_seg2;

	int32_t ins_sum_num, del_sum_num, target_var_op, mode_num[3], indel_mode, mode_tmp, max_idx, max_val;
	double ins_sum_mean0, del_sum_mean0, ins_sum_mean, del_sum_mean, mean_size, dif;

	// compute mean size
	ins_sum_mean0 = del_sum_mean0 = ins_sum_num = del_sum_num = 0;
	for(i=0; i<3; i++) mode_num[i] = 0;
	for(i=0; i<q_cluster.size(); i++){
		q_node = q_cluster.at(i);
		if(q_node->overlap_sig_num==0) continue; // skip items of no signatures

		if(q_node->ins_sum>0){
			ins_sum_mean0 += q_node->ins_sum;
			ins_sum_num ++;
		}
		if(q_node->del_sum>0){
			del_sum_mean0 += q_node->del_sum;
			del_sum_num ++;
		}
	}
	if(ins_sum_num>0) ins_sum_mean0 /= ins_sum_num;
	if(del_sum_num>0) del_sum_mean0 /= del_sum_num;

	// skip outliers
	ins_sum_mean = del_sum_mean = ins_sum_num = del_sum_num = 0;
	for(i=0; i<q_cluster.size(); i++){
		q_node = q_cluster.at(i);
		if(q_node->overlap_sig_num==0) continue; // skip items of no signatures

		if(q_node->ins_sum>0 and q_node->ins_sum<3*ins_sum_mean0){
			ins_sum_mean += q_node->ins_sum;
			ins_sum_num ++;
		}
		if(q_node->del_sum>0 and q_node->del_sum<3*del_sum_mean0){
			del_sum_mean += q_node->del_sum;
			del_sum_num ++;
		}
		if(q_node->ins_sum>0 and q_node->ins_sum<3*ins_sum_mean0 and q_node->del_sum>0 and q_node->del_sum<3*del_sum_mean0) mode_num[0] ++;
		else if(q_node->ins_sum>0 and q_node->ins_sum<3*ins_sum_mean0) mode_num[1] ++;
		else if(q_node->del_sum>0 and q_node->del_sum<3*del_sum_mean0) mode_num[2] ++;
	}
	if(ins_sum_num>0) ins_sum_mean /= ins_sum_num;
	if(del_sum_num>0) del_sum_mean /= del_sum_num;

	if(ins_sum_mean>del_sum_mean) { target_var_op = BAM_CINS; mean_size = ins_sum_mean; }
	else { target_var_op = BAM_CDEL; mean_size = del_sum_mean; }

	// compute mode
	max_idx = -1;
	max_val = INT_MIN;
	for(i=0; i<3; i++){
		if(mode_num[i]>=min_supp_num and mode_num[i]>max_val) {
			max_val = mode_num[i];
			max_idx = i;
		}
	}
	if(max_idx==0) indel_mode = MIX_MODE;
	else if(max_idx==1) indel_mode = INS_ONLY_MODE;
	else if(max_idx==2) indel_mode = DEL_ONLY_MODE;

	// compute difference
	for(i=0; i<q_cluster.size(); i++){
		q_node = q_cluster.at(i);
		if(target_var_op==BAM_CINS) dif = q_node->ins_sum - mean_size;
		else dif = q_node->del_sum - mean_size;
		dif_vec.push_back(dif);
	}

	// consider the number of signatures
	min_dif = INT_MAX;
	min_idx = -1;
	for(i=0; i<q_cluster.size(); i++){
		q_node = q_cluster.at(i);
		if(q_node->qcSig_vec.size()==query_seq_info_node->qcSig_vec.size() and q_node->entire_flanking_flag and min_dif>abs(dif_vec.at(i))){
			if(q_node->ins_sum>0 and q_node->del_sum>0) mode_tmp = MIX_MODE;
			else if(q_node->ins_sum>0) mode_tmp = INS_ONLY_MODE;
			else if(q_node->del_sum>0) mode_tmp = DEL_ONLY_MODE;
			else mode_tmp = UNUSED_MODE;
			if(mode_tmp==indel_mode){
				min_dif = abs(dif_vec.at(i));
				min_idx = i;
			}
		}
	}
	if(min_idx==-1){
		min_dif = INT_MAX;
		min_idx = -1;
		for(i=0; i<q_cluster.size(); i++){
			q_node = q_cluster.at(i);
			if(q_node->entire_flanking_flag and min_dif>abs(dif_vec.at(i))){
				if(q_node->ins_sum>0 and q_node->del_sum>0) mode_tmp = MIX_MODE;
				else if(q_node->ins_sum>0) mode_tmp = INS_ONLY_MODE;
				else if(q_node->del_sum>0) mode_tmp = DEL_ONLY_MODE;
				else mode_tmp = UNUSED_MODE;
				if(mode_tmp==indel_mode){
					min_dif = abs(dif_vec.at(i));
					min_idx = i;
				}
			}
		}
		if(min_idx==-1) min_idx = q_cluster.size() / 2;
	}

	seed_id = -1;
	span = 0;
	startSpanPos = -1;
	endSpanPos = -1;
	if(min_idx!=-1){
		q_node = q_cluster.at(min_idx);
		aln_seg1 = q_node->query_alnSegs.at(q_node->query_alnSegs.size() - 1);
		aln_seg2 = query_seq_info_node->query_alnSegs.at(query_seq_info_node->query_alnSegs.size() - 1);
		end_ref_pos1 = getEndRefPosAlnSeg(aln_seg1->startRpos, aln_seg1->opflag, aln_seg1->seglen);
		end_ref_pos2 = getEndRefPosAlnSeg(aln_seg2->startRpos, aln_seg2->opflag, aln_seg2->seglen);
		startSpanPos_tmp = max(q_node->query_alnSegs.at(0)->startRpos, query_seq_info_node->query_alnSegs.at(0)->startRpos);
		endSpanPos_tmp = min(end_ref_pos1, end_ref_pos2);

		span = endSpanPos_tmp - startSpanPos_tmp;
		startSpanPos = startSpanPos_tmp;
		endSpanPos = endSpanPos_tmp;
		seed_id = min_idx;
	}

	seed_info = new struct seedQueryInfo();
	seed_info->startSpanPos = startSpanPos;
	seed_info->endSpanPos = endSpanPos;
	seed_info->span_intersection = span;
	seed_info->id = seed_id;
	seed_info->indel_mode = indel_mode;

	return seed_info;
}

// choose the seed query for indel clustering
//struct seedQueryInfo* indelRegCluster::chooseSeedClusterQuery2(struct querySeqInfoNode* query_seq_info_node, vector<struct querySeqInfoNode*> &q_cluster){
//	size_t seed_id, j;
//	int64_t startSpanPos_tmp, endSpanPos_tmp, startSpanPos, endSpanPos, span, span_max;
//	struct seedQueryInfo* seed_info;
//
//	span_max = 0, seed_id = 0;
//	startSpanPos = 0;
//	endSpanPos = 0;
//	for(j=0;j<q_cluster.size();j++){
//		startSpanPos_tmp = max(q_cluster.at(j)->query_alnSegs.at(0)->startRpos, query_seq_info_node->query_alnSegs.at(0)->startRpos);
//		endSpanPos_tmp = min(q_cluster.at(j)->query_alnSegs.at(q_cluster.at(j)->query_alnSegs.size() - 1)->startRpos + q_cluster.at(j)->query_alnSegs.at(q_cluster.at(j)->query_alnSegs.size() - 1)->seglen -1, query_seq_info_node->query_alnSegs.at(query_seq_info_node->query_alnSegs.size() - 1)->startRpos + query_seq_info_node->query_alnSegs.at(query_seq_info_node->query_alnSegs.size() - 1)->seglen - 1);
//		if(startSpanPos_tmp < endSpanPos_tmp){
//			span = endSpanPos_tmp - startSpanPos_tmp;
//			if(span>span_max){
//				span_max = span;
//				startSpanPos = startSpanPos_tmp;
//				endSpanPos = endSpanPos_tmp;
//				seed_id = j;
//			}
//		}
//	}
//
//	seed_info = new struct seedQueryInfo();
//	seed_info->startSpanPos = startSpanPos;
//	seed_info->endSpanPos = endSpanPos;
//	seed_info->span_intersection = span_max;
//	seed_info->id = seed_id;
//	return seed_info;
//}

// choose the seed query for indel reclustering
struct seedQueryInfo* indelRegCluster::chooseSeedClusterQueryRecluster(struct seedQueryInfo* seedinfo, vector<struct querySeqInfoNode*> &q_cluster){
	size_t seed_id, j;
	int64_t startSpanPos_tmp, endSpanPos_tmp, startSpanPos, endSpanPos, span, span_max, query_start, query_end;
	struct seedQueryInfo* seed_info;

	span_max = 0, seed_id = 0;
	startSpanPos = 0;
	endSpanPos = 0;
	for(j=0;j<q_cluster.size();j++){
		query_start = q_cluster.at(j)->query_alnSegs.at(0)->startRpos;
		query_end = q_cluster.at(j)->query_alnSegs.at(q_cluster.at(j)->query_alnSegs.size() - 1)->startRpos + q_cluster.at(j)->query_alnSegs.at(q_cluster.at(j)->query_alnSegs.size() - 1)->seglen -1;
		if(query_start<seedinfo->startSpanPos){
			if(query_end>seedinfo->endSpanPos){
				startSpanPos_tmp = seedinfo->startSpanPos;
				endSpanPos_tmp = seedinfo->endSpanPos;
				span = endSpanPos_tmp - startSpanPos_tmp;
			}
			if(query_end>seedinfo->startSpanPos and query_end<=seedinfo->endSpanPos){
				startSpanPos_tmp = seedinfo->startSpanPos;
				endSpanPos_tmp = query_end;
				span = endSpanPos_tmp - startSpanPos_tmp;
			}
			if(query_end<seedinfo->startSpanPos){
				startSpanPos_tmp = 0;
				endSpanPos_tmp = 0;
				span = 0;
			}
		}else{
			if(query_start<=seedinfo->endSpanPos){
				if(query_end>seedinfo->endSpanPos){
					startSpanPos_tmp = query_start;
					endSpanPos_tmp = seedinfo->endSpanPos;
					span = endSpanPos_tmp - startSpanPos_tmp;
				}else{
					startSpanPos_tmp = query_start;
					endSpanPos_tmp = query_end;
					span = endSpanPos_tmp - startSpanPos_tmp;
				}
			}else{
				startSpanPos_tmp = 0;
				endSpanPos_tmp = 0;
				span = 0;
			}
		}

		if(span>span_max){
			span_max = span;
			startSpanPos = startSpanPos_tmp;
			endSpanPos = endSpanPos_tmp;
			seed_id = j;
		}
	}

	seed_info = new struct seedQueryInfo();
	seed_info->startSpanPos = startSpanPos;
	seed_info->endSpanPos = endSpanPos;
	seed_info->span_intersection = span_max;
	seed_info->id = seed_id;
	return seed_info;
}

vector<int32_t> indelRegCluster::getInitGroupIdx(vector<mrmin_match_t *> &mrmin_match_vec, vector<struct querySeqInfoNode*> &qseq_info_all){
	vector<int32_t> init_idx_vec;
	mrmin_match_t *mrmin_match_node;
	size_t i, j, num, max_num_MR0;
	double num_ratio;
	struct querySeqInfoNode *q_node;

	int32_t ins_sum_num, del_sum_num, mode_num[3], max_idx, sec_idx, max_val, sec_val;
	double ins_sum_mean0, del_sum_mean0, ins_sum_mean, del_sum_mean;
	bool cons_flag;

	// select the target sub-set
	num = 0;
	for(i=0; i<mrmin_match_vec.size(); i++){
		mrmin_match_node = mrmin_match_vec.at(i);
		if(mrmin_match_node->MR>0 and mrmin_match_node->num>=min_supp_num) num ++;
	}

	max_num_MR0 = 1;
	if(num==0) max_num_MR0 = 2;
	num = 0;
	for(i=0; i<mrmin_match_vec.size(); i++){
		mrmin_match_node = mrmin_match_vec.at(i);
		if(mrmin_match_node->MR>0 and mrmin_match_node->num>=min_supp_num)
			init_idx_vec.push_back(i);
		else{
			if(num<max_num_MR0){
				// prefer the item of reclustered with MR_recluster>0
				if(mrmin_match_node->MR==0 and mrmin_match_node->MR_recluster==0){
					for(j=i+1; j<mrmin_match_vec.size(); j++){
						if(mrmin_match_vec.at(j)->MR==0 and mrmin_match_vec.at(j)->MR_recluster>0 and mrmin_match_vec.at(j)->num>=min_supp_num){
							num_ratio = 0;
							if(mrmin_match_vec.at(j)->num<mrmin_match_node->num) num_ratio = (double)mrmin_match_vec.at(j)->num / mrmin_match_node->num;
							else num_ratio = (double)mrmin_match_node->num / mrmin_match_vec.at(j)->num;
							if(num_ratio>=0.6){
								mrmin_match_node = mrmin_match_vec.at(i);
								mrmin_match_vec.at(i) = mrmin_match_vec.at(j);
								mrmin_match_vec.at(j) = mrmin_match_node;
							}
							break;
						}
					}

					// check the consistency, including the mode and divergence
					// compute mean size
					ins_sum_mean0 = del_sum_mean0 = ins_sum_num = del_sum_num = 0;
					for(j=0; j<3; j++) mode_num[j] = 0;
					for(j=0; j<(size_t)mrmin_match_node->num; j++){
						q_node = qseq_info_all.at(mrmin_match_node->pos_id.at(j));
						if(q_node->overlap_sig_num==0) continue; // skip items of no signatures

						if(q_node->ins_sum>0){
							ins_sum_mean0 += q_node->ins_sum;
							ins_sum_num ++;
						}
						if(q_node->del_sum>0){
							del_sum_mean0 += q_node->del_sum;
							del_sum_num ++;
						}
					}
					if(ins_sum_num>0) ins_sum_mean0 /= ins_sum_num;
					if(del_sum_num>0) del_sum_mean0 /= del_sum_num;

					// skip outliers
					ins_sum_mean = del_sum_mean = ins_sum_num = del_sum_num = 0;
					for(j=0; j<(size_t)mrmin_match_node->num; j++){
						q_node = qseq_info_all.at(mrmin_match_node->pos_id.at(j));
						if(q_node->overlap_sig_num==0) continue; // skip items of no signatures

						if(q_node->ins_sum>0 and q_node->ins_sum<2*ins_sum_mean0){
							ins_sum_mean += q_node->ins_sum;
							ins_sum_num ++;
						}
						if(q_node->del_sum>0 and q_node->del_sum<2*del_sum_mean0){
							del_sum_mean += q_node->del_sum;
							del_sum_num ++;
						}
						if(q_node->ins_sum>0 and q_node->ins_sum<3*ins_sum_mean0 and q_node->del_sum>0 and q_node->del_sum<3*del_sum_mean0) mode_num[0] ++;
						else if(q_node->ins_sum>0 and q_node->ins_sum<3*ins_sum_mean0) mode_num[1] ++;
						else if(q_node->del_sum>0 and q_node->del_sum<3*del_sum_mean0) mode_num[2] ++;
					}
					if(ins_sum_num>0) ins_sum_mean /= ins_sum_num;
					if(del_sum_num>0) del_sum_mean /= del_sum_num;

					// compute mode
					max_idx = sec_idx = -1;
					max_val = sec_val = 0;
					for(j=0; j<3; j++){
						if(mode_num[j]>max_val) {
							sec_val = max_val;
							sec_idx = max_idx;
							max_val = mode_num[j];
							max_idx = j;
						}else if(mode_num[j]>sec_val){
							sec_val = mode_num[j];
							sec_idx = j;
						}
					}

					cons_flag = true;
					num_ratio = (double)sec_val / max_val;
					//max_pct = (double)max_val / mrmin_match_node->num;
					if(sec_val>=min_supp_num and num_ratio>0.6){
						cons_flag = false;
					}

#if POA_ALIGN_DEBUG
					cout << "ins_sum_mean=" << ins_sum_mean << ", ins_sum_num=" << ins_sum_num << endl;
					cout << "del_sum_mean=" << del_sum_mean << ", del_sum_num=" << del_sum_num << endl;
					cout << "cons_flag=" << cons_flag << ", sec/max=" << num_ratio << ", mode_arr=[" << mode_num[0] << "," << mode_num[1] << "," << mode_num[2] << "]" << endl;
#endif

					//if(i==0 and mrmin_match_node->num>=min_supp_num and mrmin_match_node->num>mrmin_match_vec.at(i+1)->num*2){
					if(cons_flag and mrmin_match_node->num>=min_supp_num and (i+1<mrmin_match_vec.size() and mrmin_match_node->num>=mrmin_match_vec.at(i+1)->num*1.5)){
						init_idx_vec.push_back(i);
						num ++;
					}
				}else if(mrmin_match_node->num>=0.6*min_supp_num){
					init_idx_vec.push_back(i);
					num ++;
				}
			}
		}
		if(init_idx_vec.size()>=2) break;
	}

#if POA_ALIGN_DEBUG
	if(init_idx_vec.size()>0){
		cout << "init_idx_vec: [" << init_idx_vec.at(0);
		for(i=1; i<init_idx_vec.size(); i++) cout << "," << init_idx_vec.at(i);
		cout << "]" << endl;
	}
#endif

	return init_idx_vec;
}

// choose the seed query for indel clustering
struct seedQueryInfo* indelRegCluster::getInitClusterItem(vector<int32_t> &qseq_id_vec_sub_group, vector<struct querySeqInfoNode*> &qseq_info_all){
	size_t seed_id, i;
	int64_t startSpanPos, endSpanPos, span, end_ref_pos;
	struct seedQueryInfo* seed_info;
	struct querySeqInfoNode *q_node;
	int32_t dif, min_idx, min_dif;
	vector<double> dif_vec;
	struct alnSeg *aln_seg;

	int32_t ins_sum_num, del_sum_num, target_var_op;
	double ins_sum_mean, del_sum_mean, mean_size;

	// compute mean size by two rounds processing
	ins_sum_mean = del_sum_mean = ins_sum_num = del_sum_num = 0;
	for(i=0; i<qseq_id_vec_sub_group.size(); i++){
		q_node = qseq_info_all.at(qseq_id_vec_sub_group.at(i));
		if(q_node->overlap_sig_num==0) continue; // skip items of no signatures

		if(q_node->ins_sum>MIN_SV_SIZE_USR*2){
			ins_sum_mean += q_node->ins_sum;
			ins_sum_num ++;
		}
		if(q_node->del_sum>MIN_SV_SIZE_USR*2){
			del_sum_mean += q_node->del_sum;
			del_sum_num ++;
		}
	}
	if(ins_sum_num>=min_supp_num) ins_sum_mean /= ins_sum_num;
	else ins_sum_mean = 0;
	if(del_sum_num>=min_supp_num) del_sum_mean /= del_sum_num;
	else del_sum_mean = 0;

#if POA_ALIGN_DEBUG
	cout << "ins_sum_mean=" << ins_sum_mean << ", ins_sum_num=" << ins_sum_num << ", del_sum_mean=" << del_sum_mean << ", del_sum_num=" << del_sum_num << endl;
#endif

	if(ins_sum_mean>0 or del_sum_mean>0){ // found
		if(ins_sum_mean>del_sum_mean) { target_var_op = BAM_CINS; mean_size = ins_sum_mean; }
		else { target_var_op = BAM_CDEL; mean_size = del_sum_mean; }

#if POA_ALIGN_DEBUG
	cout << "round=0" << ", target_var_op=" << target_var_op << ", mean_size=" << mean_size << endl;
#endif
	}else{ // try next round
		target_var_op = -1; mean_size = 0;
		ins_sum_mean = del_sum_mean = ins_sum_num = del_sum_num = 0;
		for(i=0; i<qseq_id_vec_sub_group.size(); i++){
			q_node = qseq_info_all.at(qseq_id_vec_sub_group.at(i));
			if(q_node->overlap_sig_num==0) continue; // skip items of no signatures

			if(q_node->ins_sum>0){
				ins_sum_mean += q_node->ins_sum;
				ins_sum_num ++;
			}
			if(q_node->del_sum>0){
				del_sum_mean += q_node->del_sum;
				del_sum_num ++;
			}
		}
		if(ins_sum_num>0) ins_sum_mean /= ins_sum_num;
		if(del_sum_num>0) del_sum_mean /= del_sum_num;

		if(ins_sum_mean>del_sum_mean) { target_var_op = BAM_CINS; mean_size = ins_sum_mean; }
		else { target_var_op = BAM_CDEL; mean_size = del_sum_mean; }

#if POA_ALIGN_DEBUG
	cout << "ins_sum_mean=" << ins_sum_mean << ", ins_sum_num=" << ins_sum_num << ", del_sum_mean=" << del_sum_mean << ", del_sum_num=" << del_sum_num << endl;
	cout << "round=1" << ", target_var_op=" << target_var_op << ", mean_size=" << mean_size << endl;
#endif
	}

	// compute difference
	for(i=0; i<qseq_id_vec_sub_group.size(); i++){
		q_node = qseq_info_all.at(qseq_id_vec_sub_group.at(i));
		if(target_var_op==BAM_CINS) dif = q_node->ins_sum - mean_size;
		else dif = q_node->del_sum - mean_size;
		dif_vec.push_back(dif);

#if POA_ALIGN_DEBUG
		cout << "\tqname=" << q_node->qname << ", target_var_op=" << target_var_op << ", dif=" << dif << ", mean_size=" << mean_size << endl;
#endif
	}

//	ins_sum = del_sum = 0;
//	for(i=0; i<qseq_id_vec_sub_group.size(); i++){
//		q_node = qseq_info_all.at(qseq_id_vec_sub_group.at(i));
//		if(q_node->ins_sum>=q_node->del_sum) ins_sum += q_node->ins_sum;
//		else del_sum += q_node->del_sum;
//	}
//
//	if(ins_sum>=del_sum) target_sv_type = BAM_CINS;
//	else target_sv_type = BAM_CDEL;
//
//	// compute difference
//	for(i=0; i<qseq_id_vec_sub_group.size(); i++){
//		q_node = qseq_info_all.at(qseq_id_vec_sub_group.at(i));
//		if(target_var_op==BAM_CINS) dif = q_node->ins_sum - mean_size;
//		else dif = q_node->del_sum - mean_size;
//		dif_vec.push_back(dif);
//	}

//	aver_size = 0;
//	for(i=0; i<qseq_id_vec_sub_group.size(); i++){
//		q_node = qseq_info_all.at(qseq_id_vec_sub_group.at(i));
//		if(target_sv_type==BAM_CINS) aver_size += q_node->ins_sum;
//		else aver_size += q_node->del_sum;
//	}
//	if(aver_size>0) aver_size /= qseq_id_vec_sub_group.size();
//
//	for(i=0; i<qseq_id_vec_sub_group.size(); i++){
//		q_node = qseq_info_all.at(qseq_id_vec_sub_group.at(i));
//		if(target_sv_type==BAM_CINS) dif = q_node->ins_sum - aver_size;
//		else dif = q_node->del_sum - aver_size;
//		dif_vec.push_back(dif);
//	}

	min_dif = INT_MAX;
	min_idx = -1;
	for(i=0; i<qseq_id_vec_sub_group.size(); i++){
		q_node = qseq_info_all.at(qseq_id_vec_sub_group.at(i));
		if(q_node->entire_flanking_flag and min_dif>abs(dif_vec.at(i))){
			min_dif = abs(dif_vec.at(i));
			min_idx = i;
		}
	}
	if(min_idx==-1) min_idx = qseq_id_vec_sub_group.size() / 2;

	seed_id = -1;
	span = 0;
	startSpanPos = -1;
	endSpanPos = -1;
	if(min_idx!=-1){
		q_node = qseq_info_all.at(qseq_id_vec_sub_group.at(min_idx));
		aln_seg = q_node->query_alnSegs.at(q_node->query_alnSegs.size() - 1);
		end_ref_pos = getEndRefPosAlnSeg(aln_seg->startRpos, aln_seg->opflag, aln_seg->seglen);
		startSpanPos = q_node->query_alnSegs.at(0)->startRpos;
		endSpanPos = end_ref_pos;
		span = endSpanPos - startSpanPos;
		seed_id = qseq_id_vec_sub_group.at(min_idx);
	}

	seed_info = new struct seedQueryInfo();
	seed_info->startSpanPos = startSpanPos;
	seed_info->endSpanPos = endSpanPos;
	seed_info->span_intersection = span;
	seed_info->id = seed_id;
	if(target_var_op==BAM_CINS) seed_info->indel_mode = INS_ONLY_MODE;
	else if(target_var_op==BAM_CDEL) seed_info->indel_mode = DEL_ONLY_MODE;
	else seed_info->indel_mode = UNUSED_MODE;

#if POA_ALIGN_DEBUG
	cout << "min_idx=" << min_idx << ", id=" << seed_info->id << ", startSpanPos=" << seed_info->startSpanPos << ", endSpanPos=" << seed_info->endSpanPos << ", span_intersection=" << seed_info->span_intersection << endl;
#endif

	return seed_info;
}

// choose the seed query for indel clustering
vector<struct querySeqInfoNode*> indelRegCluster::getInitClusterItemsMR0(vector<int32_t> &qseq_id_vec, vector<struct querySeqInfoNode*> &qseq_info_all){
	vector<struct querySeqInfoNode*> qseq_info_vec;
	struct querySeqInfoNode *qnode1, *qnode2;
	size_t i;
	struct seedQueryInfo *seed_info;
	int32_t var_size, target_op;
	double match_ratio, size_ratio;

	seed_info = getInitClusterItem(qseq_id_vec, qseq_info_all);
	qnode1 = qseq_info_all.at(seed_info->id);

	if(seed_info->indel_mode==INS_ONLY_MODE) { target_op = BAM_CINS; var_size = qnode1->ins_sum; }
	else if(seed_info->indel_mode==DEL_ONLY_MODE) { target_op = BAM_CDEL; var_size = qnode1->del_sum; }
	else{
		cout << "indel_mode=" << seed_info->indel_mode << ", error!" << endl;
		exit(1);
	}

	for(i=0; i<qseq_id_vec.size(); i++){
		qnode2 = qseq_info_all.at(qseq_id_vec.at(i));
		match_ratio = computeMatchRatio(qnode1, qnode2, seed_info->startSpanPos, seed_info->endSpanPos, min_sv_size, CNS_MIN_IDENTITY_CLUSTER_THRES);
		if(match_ratio>0){
			if(target_op==BAM_CINS){ // INS
				if(var_size<qnode2->ins_sum) size_ratio = (double)var_size / qnode2->ins_sum;
				else size_ratio = (double)qnode2->ins_sum / var_size;
			}else{ // DEL
				if(var_size<qnode2->del_sum) size_ratio = (double)var_size / qnode2->del_sum;
				else size_ratio = (double)qnode2->del_sum / var_size;
			}

			if(size_ratio>=CNS_MIN_SIZE_RATIO_MR0_THRES){
				qnode2->cluster_finished_flag = true;
				qseq_info_vec.push_back(qnode2);

#if POA_ALIGN_DEBUG
				cout << "\tqname=" << qnode2->qname << ", added, size_ratio=" << size_ratio << ", qseq_id=" << qseq_id_vec.at(i) << endl;
#endif
			}

#if POA_ALIGN_DEBUG
		cout << "size_ratio=" << size_ratio << ", qseq_id=" << qseq_id_vec.at(i) << ", match_ratio=" << match_ratio << endl;
#endif
		}
	}
	delete seed_info;

	return qseq_info_vec;
}


double indelRegCluster::computeMatchRatio(struct querySeqInfoNode* query_seq_info_node, struct querySeqInfoNode* q_cluster_node, int64_t startSpanPos, int64_t endSpanPos, int32_t min_sv_size, double min_identity_match){
	double score_ratio;
	int32_t score_sum;
	size_t i;
	queryCluSig_t queryCluSig, seed_qcQuery;
	bool mergeFlag, mergeFlag_tmp1, mergeFlag_tmp2;

#if POA_ALIGN_DEBUG
	cout << "qname1=" << query_seq_info_node->qname << ", qname2=" << q_cluster_node->qname << endl;
#endif

	score_ratio = 0;
	mergeFlag_tmp1 = mergeFlag_tmp2 = mergeFlag = false;
	// matching prepare
	// queryCluSig = new queryCluSig_t();
	// queryCluSig->qcSig_vec = extractQcSigsFromAlnSegsSingleQuery(query_seq_info_node, query_seq_info_node->clip_aln->chrname, startSpanPos, endSpanPos, min_sv_size); // deleted on 2025-01-11
	queryCluSig.qcSig_vec = query_seq_info_node->qcSig_vec;
	// seed_qcQuery = new queryCluSig_t();
	// seed_qcQuery->qcSig_vec = extractQcSigsFromAlnSegsSingleQuery(q_cluster_node, q_cluster_node->clip_aln->chrname, startSpanPos, endSpanPos, min_sv_size); // deleted on 2025-01-11
	seed_qcQuery.qcSig_vec = q_cluster_node->qcSig_vec;

	//merge queryCluSig->qcSig_vec and seed_qcQuery->qcSig_vec
	// mergeFlag_tmp1 = mergeNeighbouringSigsFlag(query_seq_info_node, queryCluSig->qcSig_vec, 800, 0.7, fai);
	// mergeFlag_tmp2 = mergeNeighbouringSigsFlag(q_cluster_node, seed_qcQuery->qcSig_vec, 800, 0.7, fai);
	mergeFlag_tmp1 = query_seq_info_node->qcSig_merge_flag;
	mergeFlag_tmp2 = q_cluster_node->qcSig_merge_flag;
	if(mergeFlag_tmp1 or mergeFlag_tmp2) mergeFlag = true;

	//matching
	// queryCluSig->match_profile_vec = computeQcMatchProfileSingleQuery(queryCluSig, seed_qcQuery);
	queryCluSig.match_profile_vec = computeQcMatchProfileSingleQuery(&queryCluSig, &seed_qcQuery, mergeFlag, min_identity_match);
	if(queryCluSig.match_profile_vec.size()==0){
		//score_ratio = 2;//have no varsig
		score_ratio = 0; //0.99
	}else{
		score_sum = 0;
		for(i=0; i<queryCluSig.match_profile_vec.size(); i++) score_sum += queryCluSig.match_profile_vec.at(i);
		score_ratio = (double)score_sum / queryCluSig.match_profile_vec.size();
	}
#if POA_ALIGN_DEBUG
	cout << "score_ratio=" << score_ratio << endl;
#endif

	// destroyQueryQcSig(queryCluSig.qcSig_vec);
	// destroyQueryQcSig(seed_qcQuery.qcSig_vec);

	return score_ratio;
}

vector<int8_t> indelRegCluster::computeQcMatchProfileSingleQuery(queryCluSig_t *queryCluSig, queryCluSig_t *seed_qcQuery, bool mergeFlag, double min_identity_match){
	vector<int8_t> match_profile_vec;
	int32_t rowsNum, colsNum, matchScore, mismatchScore, gapScore, gapOpenScore;
	int32_t scoreIJ, tmp_gapScore1, tmp_gapScore2, maxValue, path_val, maxValue_ul, maxValue_l, maxValue_u;
	struct alnScoreNode *scoreArr;
	int64_t i, j, arrSize;
	bool matchFlag;

	matchScore = GT_SIG_MATCH_SCORE;
	mismatchScore = GT_SIG_MISMATCH_SCORE;
	gapScore = GT_SIG_GAP_SCORE;
	gapOpenScore = GT_SIG_GAP_OPEN_SCORE;


	if(queryCluSig->qcSig_vec.size()>0 or seed_qcQuery->qcSig_vec.size()>0){
		rowsNum = queryCluSig->qcSig_vec.size() + 1;
		colsNum = seed_qcQuery->qcSig_vec.size() + 1;

		arrSize = rowsNum * colsNum;
		scoreArr = (struct alnScoreNode*) calloc (arrSize, sizeof(struct alnScoreNode));
		if(scoreArr==NULL){
			cerr << "line=" << __LINE__ << ", rowsNum=" << rowsNum << ", colsNum=" << colsNum << ", cannot allocate memory, error!" << endl;
			exit(1);
		}

		// set the elements of the first row and the first column to be zero
		scoreArr[0].path_val = 0;
		for(j=1; j<colsNum; j++) { scoreArr[j].path_val = 2; scoreArr[j].score = -j; }
		for(i=1; i<rowsNum; i++) { scoreArr[i*colsNum].path_val = 1; scoreArr[i*colsNum].score = -i; }

		// compute the scores of each element
		for(i=1; i<rowsNum; i++){
			for(j=1; j<colsNum; j++){
				if(mergeFlag)
					matchFlag = isQcSigMatch(queryCluSig->qcSig_vec.at(i-1), seed_qcQuery->qcSig_vec.at(j-1), MAX_DIST_MATCH_INDEL_MERGE, QC_SIZE_RATIO_MATCH_THRES_INDEL, min_identity_match, fai);
				else
					// matchFlag = isQcSigMatch(queryCluSig->qcSig_vec.at(i-1), seed_qcQuery->qcSig_vec.at(j-1), MAX_DIST_MATCH_INDEL, QC_SIZE_RATIO_MATCH_THRES_INDEL, QC_IDENTITY_RATIO_MATCH_THRES, fai);
					matchFlag = isQcSigMatch(queryCluSig->qcSig_vec.at(i-1), seed_qcQuery->qcSig_vec.at(j-1), MAX_DIST_MATCH_INDEL, QC_SIZE_RATIO_MATCH_THRES_INDEL, min_identity_match, fai);
				if(matchFlag) scoreIJ = matchScore;
				else scoreIJ = mismatchScore;

				if(scoreArr[(i-1)*colsNum+j].path_val!=1)//up
					tmp_gapScore1 = gapOpenScore;
				else
					tmp_gapScore1 = gapScore;

				if(scoreArr[i*colsNum+j-1].path_val!=2)//left
					tmp_gapScore2 = gapOpenScore;//-4
				else
					tmp_gapScore2 = gapScore;//-2

				maxValue = maxValue_ul = maxValue_u = maxValue_l = INT_MIN;
				path_val = -1;
				// compute the maximal score
				if(scoreArr[(i-1)*colsNum+j-1].score+scoreIJ>maxValue) {// from (i-1, j-1)
					maxValue = scoreArr[(i-1)*colsNum+j-1].score + scoreIJ;
					path_val = 0;
					maxValue_ul = maxValue;
				}
				if(scoreArr[(i-1)*colsNum+j].score+tmp_gapScore1>maxValue) {// from (i-1, j) up
					maxValue = scoreArr[(i-1)*colsNum+j].score + tmp_gapScore1;
					path_val = 1;
					maxValue_u = maxValue;
				}
				if(scoreArr[i*colsNum+j-1].score+tmp_gapScore2>maxValue) {// from (i, j-1) left
					maxValue = scoreArr[i*colsNum+j-1].score + tmp_gapScore2;
					path_val = 2;
					maxValue_l = maxValue;
				}

				if(maxValue_ul>=maxValue_u and maxValue_ul>=maxValue_l){
					if(matchFlag==false){
						scoreArr[i*colsNum+j].ismismatch = 1;
					}else{
						scoreArr[i*colsNum+j].ismismatch = 0;
					}
				}
				scoreArr[i*colsNum+j].score = maxValue;
				scoreArr[i*colsNum+j].path_val = path_val;
			}
		}

		// cout << scoreArr[rowsNum * colsNum - 1].score << endl;
		// for (i = 0; i < rowsNum; ++i) {
		// 	for (j = 0; j < colsNum; ++j) {
		// 		cout << scoreArr[i * colsNum + j].score << "_" << scoreArr[i * colsNum + j].path_val << "\t";
		// 		// cout << scoreArr[i * colsNum + j].path_val << "\t";
		// 	}
		// 	cout << endl;
		// }

		// compute signature match profile
		match_profile_vec = qComputeSigMatchProfile(scoreArr, rowsNum, colsNum, queryCluSig, seed_qcQuery);
		free(scoreArr);
	}

	return match_profile_vec;
}

vector<int8_t> indelRegCluster::qComputeSigMatchProfile(struct alnScoreNode *scoreArr, int32_t rowsNum, int32_t colsNum, queryCluSig_t *queryCluSig, queryCluSig_t *seed_qcQuery){
	vector<int8_t> match_profile_vec;
	int32_t i = rowsNum - 1, j = colsNum - 1, value;
	if(i==0 or j==0){
		match_profile_vec.push_back(0);
	}else{
		while(i>0 or j>0){
			value = scoreArr[i*colsNum+j].path_val;
			if(value==0){ //from (i-1, j-1)
				if(scoreArr[i*colsNum+j].ismismatch==1){
					match_profile_vec.push_back(0);
				}else{
					match_profile_vec.push_back(1);
				}
				queryCluSig->qcSig_vec.at(i-1)->mate_qcSig = seed_qcQuery->qcSig_vec.at(j-1);
				i--;
				j--;
			}
			else if(value==1){ //from (i-1, j)
				//?
				match_profile_vec.push_back(0);
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
	}

//============================================
#if POA_ALIGN_DEBUG
	cout << "------------------- Profile -------------------" << endl;
	// output the profile
	vector<string> sigseq_vec1, sigseq_vec2;
	string sig_str1, sig_str2;
	i = rowsNum - 1; j = colsNum - 1;
	while(i>0 or j>0){
		value = scoreArr[i*colsNum+j].path_val;
		if(value==0){ //from (i-1, j-1)
			sig_str1 = to_string(queryCluSig->qcSig_vec.at(i-1)->cigar_op_len) + "_" + to_string(queryCluSig->qcSig_vec.at(i-1)->cigar_op);
			sig_str2 = to_string(seed_qcQuery->qcSig_vec.at(j-1)->cigar_op_len) + "_" + to_string(seed_qcQuery->qcSig_vec.at(j-1)->cigar_op);

			sigseq_vec1.push_back(sig_str1);
			sigseq_vec2.push_back(sig_str2);

			i--;
			j--;
		}
		else if(value==1){ //from (i-1, j)
			//?
			//match_profile_vec.push_back(0);
			sig_str1 = to_string(queryCluSig->qcSig_vec.at(i-1)->cigar_op_len) + "_" + to_string(queryCluSig->qcSig_vec.at(i-1)->cigar_op);
			sigseq_vec1.push_back(sig_str1);
			sigseq_vec2.push_back("-");

			i--;
		}
		else{ //from (i, j-1)
			//match_profile_vec.push_back(0);
			sig_str2 = to_string(seed_qcQuery->qcSig_vec.at(j-1)->cigar_op_len) + "_" + to_string(seed_qcQuery->qcSig_vec.at(j-1)->cigar_op);
			sigseq_vec1.push_back("-");
			sigseq_vec2.push_back(sig_str2);

			j--;
		}
		if(i==0 and j==0)
			break;
	}

	reverse(sigseq_vec1.begin(), sigseq_vec1.end());
	reverse(sigseq_vec2.begin(), sigseq_vec2.end());

	size_t k;
	for(k=0; k<sigseq_vec1.size(); k++) cout << sigseq_vec1.at(k) << "\t";
	cout << endl;
	for(k=0; k<sigseq_vec2.size(); k++) cout << sigseq_vec2.at(k) << "\t";
	cout << endl;
	for(k=0; k<match_profile_vec.size(); k++) cout << match_profile_vec.at(k) << "\t";
	cout << endl;
#endif
//=====================================

	return match_profile_vec;
}

bool indelRegCluster::getReClusterFlag(vector<mrmin_match_t *> &mr_min_match_vec){
	bool recluster_flag;
	size_t i;

	recluster_flag = false;
	if(mr_min_match_vec.size()>1){
		for(i=0;i<2;i++){
			if(mr_min_match_vec.at(i)->MR==0 and mr_min_match_vec.at(i)->num>=min_supp_num){
				recluster_flag = true;
				break;
			}
		}
	}else
		recluster_flag = false;

	return recluster_flag;
}

bool indelRegCluster::initReCluster(vector<mrmin_match_t *> &mr_vec_recluster, vector<int32_t> &min_id_for_recluster_vec){
	bool first_flag;
	size_t i, j;

	first_flag = false;
	if(mr_vec_recluster.size()>1){
		for(i=0;i<2;i++){
			if(mr_vec_recluster.at(i)->MR==0){
				for(j=0; j<mr_vec_recluster.at(i)->pos_id.size(); j++)
					min_id_for_recluster_vec.push_back(mr_vec_recluster.at(i)->pos_id.at(j));
				delete mr_vec_recluster.at(i);
				mr_vec_recluster.erase(mr_vec_recluster.begin()+i);
				if(i==0) first_flag = true;
				break;
			}
		}
	}else{
		first_flag = true;
	}

	return first_flag;
}

vector<mrmin_match_t *> indelRegCluster::reClusterByMR(vector<struct querySeqInfoNode*> &query_seq_info_all, vector<int32_t> &min_id_vec_recluster, bool mrmin_match_vec_begin_flag){
	vector<struct querySeqInfoNode*> q_recluster_a;
	vector<mrmin_match_t *> mrmin_match_vec, mrmin_match_vec_tmp;
	seedQueryInfo *seed_info_recluster;
	size_t i, j, k, round, mmin_id, total_num;
	double mr_recluster_a, sco_min;
	vector<double> mr_vec_recluster;
	bool flag, recluster_valid_flag, next_round_flag, in_break_flag, copy_flag;
	mrmin_match_t *mrmin_match_node;
	vector<int32_t> min_id_vec, *p_min_id_vec;

	recluster_valid_flag = false;
	total_num = min_id_vec_recluster.size();

	for(round=0; round<2; round++){
		if(round==0){
			p_min_id_vec = &min_id_vec_recluster;
		}else{
			p_min_id_vec = &min_id_vec;
		}
		if(p_min_id_vec->size()<(size_t)min_supp_num) continue;

		for(k=0; k<p_min_id_vec->size(); k++){
			mrmin_match_vec_tmp.clear();
			mr_vec_recluster.clear();
			sco_min = INT_MAX;
			in_break_flag = false;
			next_round_flag = false;

			//recluster based on MR
			query_seq_info_all.at(p_min_id_vec->at(k))->cluster_finished_flag = true;
			q_recluster_a.push_back(query_seq_info_all.at(p_min_id_vec->at(k)));
			//cout << "----------------k=" << p_min_id_vec->at(k) << " qname=" << query_seq_info_all.at(p_min_id_vec->at(k))->qname << "----------------" << endl;
			//compute score_ratio_vec_recluster and sco_min
			for(i=0;i<p_min_id_vec->size();i++){
				seed_info_recluster = chooseSeedClusterQuery(query_seq_info_all.at(p_min_id_vec->at(i)), q_recluster_a);
				if(seed_info_recluster->span_intersection>0){
					if(i!=k) {
						mr_recluster_a = computeMatchRatio(query_seq_info_all.at(p_min_id_vec->at(i)), q_recluster_a.at(seed_info_recluster->id), seed_info_recluster->startSpanPos, seed_info_recluster->endSpanPos, min_sv_size, min_identity_match);
						mr_vec_recluster.push_back(mr_recluster_a);
					}else{
						mr_recluster_a = 1;
						mr_vec_recluster.push_back(mr_recluster_a);
					}
					if(mr_recluster_a<sco_min){
						sco_min = mr_recluster_a;
					}
				}else{
					mr_vec_recluster.push_back(2);
				}
				delete seed_info_recluster;
			}
			query_seq_info_all.at(p_min_id_vec->at(k))->cluster_finished_flag = false;
			q_recluster_a.pop_back();

			if(sco_min<1){
				for(i=0; i<p_min_id_vec->size(); i++){
					flag = false;
					for(j=0; j<mrmin_match_vec_tmp.size(); j++){
						if(mr_vec_recluster.at(i)==mrmin_match_vec_tmp.at(j)->MR){
							mmin_id = j;
							flag = true;
							break;
						}
					}
					if(flag==false){
						mrmin_match_node = new mrmin_match_t();
						mrmin_match_node->pos_id.push_back(p_min_id_vec->at(i));
						mrmin_match_node->num = 1;
						mrmin_match_node->MR = mr_vec_recluster.at(i);
						mrmin_match_node->MR_recluster = -1;
						mrmin_match_node->candidate_flag = false;
						mrmin_match_vec_tmp.push_back(mrmin_match_node);
					}else{
						mrmin_match_vec_tmp.at(mmin_id)->num += 1;
						mrmin_match_vec_tmp.at(mmin_id)->pos_id.push_back(p_min_id_vec->at(i));
					}
				}

				//sort
				sortMRVec(mrmin_match_vec_tmp);

#if POA_ALIGN_DEBUG
				printMRVec(mrmin_match_vec_tmp);
#endif

				flag = true;
				if(mrmin_match_vec_tmp.size()==2 and mrmin_match_vec_tmp.at(0)->MR==0 and mrmin_match_vec_tmp.at(1)->MR==1 and mrmin_match_vec_tmp.at(1)->num==1) flag = false; // skip, and try next one
				else if(mrmin_match_vec_tmp.size()>2 and mrmin_match_vec_tmp.at(0)->MR==0 and mrmin_match_vec_tmp.at(0)->num>MIN_INVLAID_SECOND_GROUP_RATIO_THRES*p_min_id_vec->size()) flag = false; // skip, and try next one

				copy_flag = false;
				if(flag){
					// try next round
					if(round==0 and mrmin_match_vec_tmp.size()>=2 and ((mrmin_match_vec_tmp.at(0)->MR==0 and mrmin_match_vec_tmp.at(0)->num>=MIN_SINGLE_CLUSTER_GROUP_RATIO_THRES*total_num and mrmin_match_vec_tmp.at(0)->num>=min_supp_num) or (mrmin_match_vec_tmp.at(1)->MR==0 and mrmin_match_vec_tmp.at(1)->num>=MIN_SINGLE_CLUSTER_GROUP_RATIO_THRES*total_num and mrmin_match_vec_tmp.at(0)->num>=min_supp_num))){
						copy_flag = true;
					}

//					recluster_valid_flag = matchVecValidFlag(mrmin_match_vec_tmp, mrmin_match_vec_begin_flag);
//					if(recluster_valid_flag) break;
				}

				if(copy_flag){ // copy part of items, but need next round clustering for MR=0 items
					// copy items of MR>0
					for(i=0; i<mrmin_match_vec_tmp.size(); i++){
						mrmin_match_node = mrmin_match_vec_tmp.at(i);
						if(mrmin_match_node->MR!=0) mrmin_match_vec.push_back(mrmin_match_node);
						else{ // get query ids for MR=0
							for(j=0; j<mrmin_match_node->pos_id.size(); j++) min_id_vec.push_back(mrmin_match_node->pos_id.at(j));
							delete mrmin_match_node;
						}
					}
					mrmin_match_vec_tmp.clear();

					in_break_flag = true;
					next_round_flag = true;
				}else if(flag){ // validate and copy all items if necessary
					// validate
					recluster_valid_flag = matchVecValidFlag(mrmin_match_vec_tmp, mrmin_match_vec_begin_flag);
					if(recluster_valid_flag){ // valid, then save items
						in_break_flag = true;
						next_round_flag = false;
						for(i=0; i<mrmin_match_vec_tmp.size(); i++) mrmin_match_vec.push_back(mrmin_match_vec_tmp.at(i));
						mrmin_match_vec_tmp.clear();
					}else{ // invalid, then try next query, and release memory
						in_break_flag = false;
						for(i=0; i<mrmin_match_vec_tmp.size(); i++) delete mrmin_match_vec_tmp.at(i);
						mrmin_match_vec_tmp.clear();
					}
				}else{
					in_break_flag = false;
					recluster_valid_flag = false;
					// release memory
					for(i=0; i<mrmin_match_vec_tmp.size(); i++) delete mrmin_match_vec_tmp.at(i);
					mrmin_match_vec_tmp.clear();
				}
				if(in_break_flag) break;
			}else{ // valid single group with MR=1
				recluster_valid_flag = true;
				next_round_flag = false;
				in_break_flag = true;

				mrmin_match_node = new mrmin_match_t();
				for(i=0; i<p_min_id_vec->size(); i++) mrmin_match_node->pos_id.push_back(p_min_id_vec->at(i));
				mrmin_match_node->num = p_min_id_vec->size();
				mrmin_match_node->MR = 1;
				mrmin_match_node->MR_recluster = -1;
				mrmin_match_node->candidate_flag = false;
				mrmin_match_vec.push_back(mrmin_match_node);

				//break;
			}
			if(in_break_flag) break;
		}

		if(next_round_flag==false) break;
	}

	// cout << "" << endl;
	if(recluster_valid_flag==false){
		if(round==0){ // round one
			for(i=0; i<mrmin_match_vec.size(); i++) delete mrmin_match_vec.at(i);
			mrmin_match_vec.clear();

			mrmin_match_node = new mrmin_match_t();
			for(i=0; i<min_id_vec_recluster.size(); i++) mrmin_match_node->pos_id.push_back(min_id_vec_recluster.at(i));
			mrmin_match_node->num = min_id_vec_recluster.size();
			mrmin_match_node->MR = 0;
			mrmin_match_node->MR_recluster = -1;
			mrmin_match_node->candidate_flag = false;
			mrmin_match_vec.push_back(mrmin_match_node);
		}else{ // round two
			mrmin_match_node = new mrmin_match_t();
			for(i=0; i<min_id_vec.size(); i++) mrmin_match_node->pos_id.push_back(min_id_vec.at(i));
			mrmin_match_node->num = min_id_vec.size();
			mrmin_match_node->MR = 0;
			mrmin_match_node->MR_recluster = -1;
			mrmin_match_node->candidate_flag = false;
			mrmin_match_vec.push_back(mrmin_match_node);
		}
	}

#if POA_ALIGN_DEBUG
		cout << "recluster_valid_flag=" << recluster_valid_flag << endl;
		printMRVec(mrmin_match_vec);
#endif

	return mrmin_match_vec;
}

// sort match vector descendingly
void indelRegCluster::sortMRVec(vector<mrmin_match_t *> &match_vec){
	bool flag;
	size_t i, j;
	mrmin_match_t *mrmin_match_node;

	//sort descendingly
	for(i=0; i<match_vec.size(); i++){
		for(j=i+1; j<match_vec.size(); j++){
			flag = false;
			if(match_vec.at(i)->num<match_vec.at(j)->num) flag = true;
			else if(match_vec.at(i)->num==match_vec.at(j)->num and match_vec.at(i)->MR<match_vec.at(j)->MR) flag = true;
			if(flag){
				mrmin_match_node = match_vec.at(i);
				match_vec.at(i) = match_vec.at(j);
				match_vec.at(j) = mrmin_match_node;
			}
		}
	}
}

// sort match vector descendingly
void indelRegCluster::printMRVec(vector<mrmin_match_t *> &match_vec){
	size_t i, j;
	mrmin_match_t *node;

	cout << "match_vec.size=" << match_vec.size() << endl;
	for(i=0; i<match_vec.size(); i++){
		node = match_vec.at(i);
		cout << "[" << i << "]: num=" << node->num << ", MR=" << node->MR << ", MR_re=" << node->MR_recluster << ", pos_id={" << node->pos_id.at(0);
		for(j=1; j<(size_t)node->num; j++) cout << "," << node->pos_id.at(j);
		cout << "}" << endl;
	}
}

bool indelRegCluster::matchVecValidFlag(vector<mrmin_match_t *> &match_vec_tmp, bool match_vec_begin_flag){
	bool flag;
	int32_t equally_reads_num;
	double max_id_vec_num_mr;
	size_t i;

	flag = false;
	//max_id_vec_num = match_vec_tmp.at(0)->num;
	max_id_vec_num_mr = match_vec_tmp.at(0)->MR;
	equally_reads_num = 0;

	for(i=0;i<match_vec_tmp.size();i++){
//		if(match_vec_tmp.at(i)->num > max_id_vec_num) {
//			max_id_vec_num = match_vec_tmp.at(i)->num;
//			max_id_vec_num_mr = match_vec_tmp.at(i)->MR;
//		}
		if(match_vec_tmp.at(i)->MR==1) { equally_reads_num = match_vec_tmp.at(i)->num; break; }
	}

	if(max_id_vec_num_mr==0 and equally_reads_num < min_supp_num and match_vec_begin_flag) flag = false;
	else flag = true;

	return flag;
}

//================================= debug operation =================================
// debug operation
void indelRegCluster::printQcSigIndelCluster(vector<struct querySeqInfoVec*> &query_clu_vec){
	vector<struct querySeqInfoNode*> query_seq_info_all;
	size_t i;

	for(i=0; i<query_clu_vec.size(); i++){
		query_seq_info_all = query_clu_vec.at(i)->query_seq_info_all;
		cout << "Information for cluster " << i << ", rescue_flag=" << query_clu_vec.at(i)->rescue_flag << endl;
		printQcSigIndel(query_seq_info_all);
	}
}

// debug operation
void indelRegCluster::printQcSigIndel(vector<struct querySeqInfoNode*> &query_seq_info_all){
	size_t i, aver_ins_sum, aver_del_sum, num_ins_sum, num_del_sum;
	struct querySeqInfoNode *query_seq_info_node;

	aver_ins_sum = aver_del_sum = num_ins_sum = num_del_sum = 0;
	cout << "Query_size=" << query_seq_info_all.size() << ": " << endl;
	for(i=0; i<query_seq_info_all.size(); i++){
		query_seq_info_node = query_seq_info_all.at(i);
		printQcSigIndelSingleQuery(query_seq_info_node, i);
		if(query_seq_info_node->ins_sum>0) {
			aver_ins_sum += query_seq_info_node->ins_sum;
			num_ins_sum ++;
		}
		if(query_seq_info_node->del_sum>0){
			aver_del_sum += query_seq_info_node->del_sum;
			num_del_sum ++;
		}
	}
	if(num_ins_sum>0) aver_ins_sum = (double)aver_ins_sum / num_ins_sum;
	if(num_del_sum>0) aver_del_sum = (double)aver_del_sum / num_del_sum;
	cout << "aver_ins_sum=" << aver_ins_sum << ", num_ins_sum=" << num_ins_sum << ", aver_del_sum=" << aver_del_sum << ", num_del_sum=" << num_del_sum << endl;
}

// debug operation
void indelRegCluster::printQcSigIndelSingleQuery(struct querySeqInfoNode *query_seq_info_node, int32_t idx){
	size_t i;
	qcSig_t *sig;

	cout << "[" << idx << "]: " << query_seq_info_node->qname << ", overlap_sig_num=" << query_seq_info_node->overlap_sig_num << ", entire_flanking=" << query_seq_info_node->entire_flanking_flag << ", ins_sum=" << query_seq_info_node->ins_sum << ", del_sum=" << query_seq_info_node->del_sum << ", aver_sig_size=" << query_seq_info_node->aver_sig_size << ": ";
	for(i=0; i<query_seq_info_node->qcSig_vec.size(); i++){
		sig = query_seq_info_node->qcSig_vec.at(i);
		cout << "\t(" << sig->cigar_op << "," << sig->chrname << "," << sig->ref_pos << "," << sig->cigar_op_len << ")";
	}
	cout << endl;
}
