#include "clipRegCluster.h"

extern pthread_mutex_t mutex_fai;

clipRegCluster::clipRegCluster(string &chrname, int64_t var_startRefPos, int64_t var_endRefPos, int32_t minClipEndSize, int32_t min_sv_size, int32_t min_supp_num, faidx_t *fai) {
	this->chrname = chrname;
	this->var_startRefPos = var_startRefPos;
	this->var_endRefPos = var_endRefPos;
	this->minClipEndSize = minClipEndSize;
	this->min_sv_size = min_sv_size;
	this->min_supp_num = min_supp_num;
	this->chrlen = faidx_seq_len(fai, chrname.c_str()); // get reference size
	this->fai = fai;
}

clipRegCluster::~clipRegCluster() {
}

void clipRegCluster::destoryQcSigList(vector<qcSigList_t*> &qcSigList_vec){
	size_t i, j;
	qcSigList_t *qcSig_list;
	for(i=0; i<qcSigList_vec.size(); i++){
		qcSig_list = qcSigList_vec.at(i);
		for(j=0; j<qcSig_list->qcSig_vec.size(); j++)
			delete qcSig_list->qcSig_vec.at(j);
		vector<qcSig_t*>().swap(qcSig_list->qcSig_vec);
		delete qcSig_list;
	}
	vector<qcSigList_t*>().swap(qcSigList_vec);
}

vector<qcSigListVec_t*> clipRegCluster::queryCluster(vector<struct querySeqInfoNode*> &query_seq_info_vec){
	vector<qcSigListVec_t*> query_clu_vec;
	vector<qcSigList_t*> qcSigList_vec;
	qcSigListVec_t *q_cluster_node_a, *q_cluster_node_b, *tmp;
	size_t i, j, k, smin_id;
	qcSigList_t *qcSigList_node, *qcSigList_node2, *tmp_qcSigList_node;
	vector<struct querySeqInfoNode*> query_seqs;
	vector<qcSigList_t*> q_cluster_a, q_cluster_b, q_cluster_rescue, q_cluster_a_final, q_cluster_b_final;
	double match_ratio_a, match_ratio_b, match_min;
	vector<size_t> min_id_vec;// span_vec2_id;
	vector<double> match_ratio_vec;
	struct seedQueryInfo *seed_info_a, *seed_info_b;
	bool exist_flag;
	vector<srmin_match_t*> mrmin_match_vec;
	srmin_match_t *mrmin_match_node;

	// extract the features
	qcSigList_vec = extractQcSigsClipReg(query_seq_info_vec);

	for(i=0; i<qcSigList_vec.size(); i++) qcSigList_vec.at(i)->cluster_finished_flag = false;

	// sort in descending order according to the number of features
	for(i=0; i<qcSigList_vec.size(); i++){
		for(j=0; j<qcSigList_vec.size()-i-1; j++){
			qcSigList_node = qcSigList_vec.at(j);
			qcSigList_node2 = qcSigList_vec.at(j+1);
			if(qcSigList_node->qcSig_vec.size()<qcSigList_node2->qcSig_vec.size()){
				tmp_qcSigList_node = qcSigList_vec.at(j);
				qcSigList_vec.at(j) = qcSigList_vec.at(j+1);
				qcSigList_vec.at(j+1) = tmp_qcSigList_node;
			}
		}
	}

#if CLIP_REG_CLUSTER_DEBUG
	cout << "*********** Before cluster, size=" << qcSigList_vec.size() << ":" << endl;
	printQcSigListVec(qcSigList_vec);
#endif

	// initialize a_cluster and b_cluster
	for(k=0; k<qcSigList_vec.size(); k++){

		if(qcSigList_vec.at(k)->qcSig_vec.size()==0) continue;

		qcSigList_vec.at(k)->cluster_finished_flag = true;
		q_cluster_a.push_back(qcSigList_vec.at(k));

		//initializing b_cluster: choose one that is most different query from a_cluster add in b_cluster
		match_min = INT_MAX; //0.99
		match_ratio_vec.clear();
		min_id_vec.clear();

		match_ratio_vec.push_back(1);
		for(i=1; i<qcSigList_vec.size(); i++){
			if(qcSigList_vec.at(i)->qcSig_vec.size()>0){
				seed_info_a = chooseSeedClusterQueryClipReg(qcSigList_vec.at(i), q_cluster_a);
				if(seed_info_a->span_intersection>0){
					match_ratio_a = computeMatchRatioClipReg(qcSigList_vec.at(i), q_cluster_a.at(seed_info_a->id), seed_info_a->startSpanPos, seed_info_a->endSpanPos);
					match_ratio_vec.push_back(match_ratio_a);
					if(match_ratio_a<match_min){
						match_min = match_ratio_a;
						//smin_id = i;
						//sam_min_id.push_back(i);
					}
				}else{
					match_ratio_vec.push_back(2);
				}
				delete seed_info_a;
			}
		}
		qcSigList_vec.at(k)->cluster_finished_flag = false;
		q_cluster_a.pop_back();

//		exist_flag = false;
		if(match_min < 1){
			for(i=0; i<match_ratio_vec.size(); i++){
				if(match_ratio_vec.at(i)<=1)
					min_id_vec.push_back(i);
			}
			for(i=0; i<min_id_vec.size(); i++){
				exist_flag = false;
				for(j=0; j<mrmin_match_vec.size(); j++){
					if(match_ratio_vec.at(min_id_vec.at(i))==mrmin_match_vec.at(j)->SR){
						smin_id = j;
						exist_flag = true;
						break;
					}
				}
				if(exist_flag==false){
					mrmin_match_node = new srmin_match_t();
					mrmin_match_node->pos_id.push_back(min_id_vec.at(i));
					mrmin_match_node->num = 1;
					mrmin_match_node->SR = match_ratio_vec.at(min_id_vec.at(i));
					mrmin_match_vec.push_back(mrmin_match_node);
				}else{
					mrmin_match_vec.at(smin_id)->num += 1;
					mrmin_match_vec.at(smin_id)->pos_id.push_back(min_id_vec.at(i));
				}
			}

			// sort descendingly and initialize clusters adaptively
			for(i=0; i<mrmin_match_vec.size(); i++){
				for(j=i+1; j<mrmin_match_vec.size(); j++){
					if(mrmin_match_vec.at(i)->num<mrmin_match_vec.at(j)->num){
						mrmin_match_node = mrmin_match_vec.at(i);
						mrmin_match_vec.at(i) = mrmin_match_vec.at(j);
						mrmin_match_vec.at(j) = mrmin_match_node;
					}
				}
			}

			// only two initial clusters, deleted on 2024-07-04
//			mrsame_maxnum = mrmin_match_vec.at(0)->num;
//			smin_id = 0;
//			for(i=0; i<mrmin_match_vec.size(); i++){
//				if(mrsame_maxnum<mrmin_match_vec.at(i)->num){
//					smin_id = i;
//					mrsame_maxnum = mrmin_match_vec.at(i)->num;
//				}
//			}

			qcSigList_vec.at(mrmin_match_vec.at(0)->pos_id.at(0))->cluster_finished_flag = true;
			q_cluster_a.push_back(qcSigList_vec.at(mrmin_match_vec.at(0)->pos_id.at(0)));
			if(mrmin_match_vec.size()>1 and mrmin_match_vec.at(1)->num>=0.6*min_supp_num){
				qcSigList_vec.at(mrmin_match_vec.at(1)->pos_id.at(0))->cluster_finished_flag = true;
				q_cluster_b.push_back(qcSigList_vec.at(mrmin_match_vec.at(1)->pos_id.at(0)));
			}
		}else{
			qcSigList_vec.at(k)->cluster_finished_flag = true;
			q_cluster_a.push_back(qcSigList_vec.at(k));
		}

		if(mrmin_match_vec.size()>0){
			vector<srmin_match_t*>::iterator svp;
			for(svp=mrmin_match_vec.begin(); svp!=mrmin_match_vec.end(); svp++) delete *svp;
			vector<srmin_match_t*>().swap(mrmin_match_vec);
		}
		if(q_cluster_b.size()==0){
			q_cluster_a.pop_back();
			qcSigList_vec.at(k)->cluster_finished_flag = false;
		}else{
			break;
		}
	}

	if(q_cluster_b.size()==0){
		for(i=0; i<qcSigList_vec.size(); i++){
			if(qcSigList_vec.at(i)->qcSig_vec.size()>0 and qcSigList_vec.at(i)->cluster_finished_flag == false){
				qcSigList_vec.at(i)->cluster_finished_flag = true;
				q_cluster_a.push_back(qcSigList_vec.at(i));
			}
		}
	}else{
		for(i=0; i<qcSigList_vec.size(); i++){
			if(qcSigList_vec.at(i)->qcSig_vec.size()>0 and qcSigList_vec.at(i)->cluster_finished_flag == false){
				//choose seed query from a_cluster
				seed_info_a = chooseSeedClusterQueryClipReg(qcSigList_vec.at(i), q_cluster_a);
				seed_info_b = chooseSeedClusterQueryClipReg(qcSigList_vec.at(i), q_cluster_b);
				if(seed_info_a->span_intersection>0){
					match_ratio_a = computeMatchRatioClipReg(qcSigList_vec.at(i), q_cluster_a.at(seed_info_a->id), seed_info_a->startSpanPos, seed_info_a->endSpanPos);
					if(seed_info_b->span_intersection>0){
						match_ratio_b = computeMatchRatioClipReg(qcSigList_vec.at(i), q_cluster_b.at(seed_info_b->id), seed_info_b->startSpanPos, seed_info_b->endSpanPos);
						if(match_ratio_a>match_ratio_b){
							qcSigList_vec.at(i)->cluster_finished_flag = true;
							q_cluster_a.push_back(qcSigList_vec.at(i));
						}else{
							if(match_ratio_b==match_ratio_a){
								qcSigList_vec.at(i)->cluster_finished_flag = true;
								q_cluster_rescue.push_back(qcSigList_vec.at(i));
							}else{
								qcSigList_vec.at(i)->cluster_finished_flag = true;
								q_cluster_b.push_back(qcSigList_vec.at(i));
							}
						}
					}else{
						if(match_ratio_a<1){
							qcSigList_vec.at(i)->cluster_finished_flag = true;
							q_cluster_b.push_back(qcSigList_vec.at(i));
						}else{
							qcSigList_vec.at(i)->cluster_finished_flag = true;
							q_cluster_a.push_back(qcSigList_vec.at(i));
						}
					}
				}else{
					if(seed_info_b->span_intersection>0){
						match_ratio_b = computeMatchRatioClipReg(qcSigList_vec.at(i), q_cluster_b.at(seed_info_b->id), seed_info_b->startSpanPos, seed_info_b->endSpanPos);
						if(match_ratio_b<1){
							qcSigList_vec.at(i)->cluster_finished_flag = true;
							q_cluster_a.push_back(qcSigList_vec.at(i));
						}else{
							qcSigList_vec.at(i)->cluster_finished_flag = true;
							q_cluster_b.push_back(qcSigList_vec.at(i));
						}
					}
				}
				delete seed_info_a;
				delete seed_info_b;
			}
		}
	}

	if(q_cluster_rescue.size()>0){
		for(i=0;i<q_cluster_rescue.size();i++){
			seed_info_a = chooseSeedClusterQueryClipReg(q_cluster_rescue.at(i),q_cluster_a);
			seed_info_b = chooseSeedClusterQueryClipReg(q_cluster_rescue.at(i), q_cluster_b);
			match_ratio_a = computeMatchRatioClipReg(q_cluster_rescue.at(i), q_cluster_a.at(seed_info_a->id), seed_info_a->startSpanPos, seed_info_a->endSpanPos);
			match_ratio_b = computeMatchRatioClipReg(q_cluster_rescue.at(i), q_cluster_b.at(seed_info_b->id), seed_info_b->startSpanPos, seed_info_b->endSpanPos);
			if(match_ratio_a>match_ratio_b){
				q_cluster_a.push_back(q_cluster_rescue.at(i));
			}else{
				if(match_ratio_a!=match_ratio_b){
					q_cluster_b.push_back(q_cluster_rescue.at(i));
				}else{
					qcSigList_vec.at(i)->cluster_finished_flag = false;
				}
			}
			delete seed_info_a;
			delete seed_info_b;
		}
	}

	for(i=0; i<q_cluster_a.size(); i++) if(q_cluster_a.at(i)->qcSig_vec.size()>0) q_cluster_a_final.push_back(q_cluster_a.at(i));
	for(i=0; i<q_cluster_b.size(); i++) if(q_cluster_b.at(i)->qcSig_vec.size()>0) q_cluster_b_final.push_back(q_cluster_b.at(i));

	if(q_cluster_a_final.size()>0){
		q_cluster_node_a = new qcSigListVec_t();
		q_cluster_node_a->qcSigList = q_cluster_a_final;
		query_clu_vec.push_back(q_cluster_node_a);
	}

	if(q_cluster_b_final.size()>0){
		q_cluster_node_b = new qcSigListVec_t();
		q_cluster_node_b->qcSigList = q_cluster_b_final;
		query_clu_vec.push_back(q_cluster_node_b);
	}

	// sort
	for(i=0; i<query_clu_vec.size(); i++){
		for(j=i+1; j<query_clu_vec.size(); j++){
			if(query_clu_vec.at(i)->qcSigList.size()<query_clu_vec.at(j)->qcSigList.size()){
				tmp = query_clu_vec.at(i);
				query_clu_vec.at(i) = query_clu_vec.at(j);
				query_clu_vec.at(j) = tmp;
			}
		}
	}

#if CLIP_REG_CLUSTER_DEBUG
	for(i=0; i<query_clu_vec.size(); i++){
		cout << "*********** Queries of cluster " << i << ", size=" << query_clu_vec.at(i)->qcSigList.size() << ":" << endl;
		printQcSigListVec(query_clu_vec.at(i)->qcSigList);
	}
#endif

	// remove unclustered queries
	int32_t total = qcSigList_vec.size();
	int32_t num_clustered = 0;
	for(i=0; i<query_clu_vec.size(); i++) num_clustered += query_clu_vec.at(i)->qcSigList.size();

	int32_t num_unclustered = removeUnclusteredQueries(qcSigList_vec, query_clu_vec);
	if(num_clustered+num_unclustered != total){
		cout << "line=" << __LINE__ << ", total=" << total << ", num_clustered=" << num_clustered << ", num_unclustered=" << num_unclustered << ", error!" << endl;
		exit(1);
	}

	return query_clu_vec;
}


vector<struct querySeqInfoNode*> clipRegCluster::getQuerySeqs(string &qname, vector<struct querySeqInfoNode*> &query_seq_info_vec){
	vector<struct querySeqInfoNode*> queryseq_vec;
	struct querySeqInfoNode *queryseq_node;
	size_t i;

	for(i=0; i<query_seq_info_vec.size(); i++){
		queryseq_node = query_seq_info_vec.at(i);
		if(queryseq_node->qname.compare(qname)==0) queryseq_vec.push_back(queryseq_node);
	}

	return queryseq_vec;
}

// remove unclustered queries
int32_t clipRegCluster::removeUnclusteredQueries(vector<qcSigList_t*> &qcSigList_vec, vector<qcSigListVec_t*> &query_clu_vec){
	int32_t num_unclustered, clu_id;
	size_t i, j;
	qcSigList_t *qcSigList_node;

	num_unclustered = 0;
	for(i=0; i<qcSigList_vec.size(); ){
		qcSigList_node = qcSigList_vec.at(i);
		clu_id = getCluIDByQuery(qcSigList_node, query_clu_vec);
		if(clu_id==-1){
			num_unclustered ++;
			for(j=0; j<qcSigList_node->qcSig_vec.size(); j++)
				delete qcSigList_node->qcSig_vec.at(j);
			vector<qcSig_t*>().swap(qcSigList_node->qcSig_vec);
			delete qcSigList_node;
			qcSigList_vec.erase(qcSigList_vec.begin()+i);
		}else i++;
	}

	return num_unclustered;
}

// get the cluster id by query
int32_t clipRegCluster::getCluIDByQuery(qcSigList_t *qcSigList_node, vector<qcSigListVec_t*> &query_clu_vec){
	int32_t clu_id, clu_num;
	size_t i;
	vector<qcSigList_t*> cluster_feature_vec;

	clu_id = -1;
	clu_num = 0;
	for(i=0; i<query_clu_vec.size(); i++){
		cluster_feature_vec = query_clu_vec.at(i)->qcSigList;
		if(find(cluster_feature_vec.begin(), cluster_feature_vec.end(), qcSigList_node)!=cluster_feature_vec.end()){
			clu_id = i;
			clu_num ++;
			//break;
		}
	}

	if(clu_num>1){
		cout << "line=" << __LINE__ << ", clu_num=" << clu_num << ", error!" << endl;
		exit(1);
	}

	return clu_id;
}

// print the feature
void clipRegCluster::printQcSigListVec(vector<qcSigList_t*> &qcSigList_vec){
	size_t i, j;
	qcSigList_t *qcSigList_node;
	qcSig_t *qc_sig;

	// print the feature
	for(i=0; i<qcSigList_vec.size(); i++){
		qcSigList_node = qcSigList_vec.at(i);
		cout << "[" << i << "]: " << qcSigList_node->query_seqs_vec.at(0)->qname << ", cluster_finished=" << qcSigList_node->cluster_finished_flag << ", qcSig_vec.size=" << qcSigList_node->qcSig_vec.size() << ":";
		for(j=0; j<qcSigList_node->qcSig_vec.size(); j++){
			qc_sig = qcSigList_node->qcSig_vec.at(j);
			if(qc_sig->cigar_op==BAM_CSOFT_CLIP or qc_sig->cigar_op==BAM_CHARD_CLIP) // clippings
				cout << "\t(" << qc_sig->cigar_op << "," << qc_sig->chrname << ":" << qc_sig->ref_pos << "," << qc_sig->dist_next_clip << ")";
			else // indels
				cout << "\t(" << qc_sig->cigar_op << "," << qc_sig->chrname << ":" << qc_sig->ref_pos << "," << qc_sig->cigar_op_len << ")";
		}
		cout << endl;
	}
}

vector<qcSigList_t*> clipRegCluster::extractQcSigsClipReg(vector<struct querySeqInfoNode*> &query_seq_info_vec){
	vector<qcSigList_t*> qcSigList_vec;
	vector<qcSig_t*> qcSig_vec;
	qcSigList_t *qcSigList_node;
	struct querySeqInfoNode *queryseq_node;
	vector<struct querySeqInfoNode*> query_seqs;
	size_t i, j;

	for(i=0; i<query_seq_info_vec.size(); i++) { query_seq_info_vec.at(i)->cluster_finished_flag = false; query_seq_info_vec.at(i)->selected_flag = false; }

	// extract the features
	for(i=0; i<query_seq_info_vec.size(); i++){
		queryseq_node = query_seq_info_vec.at(i);
		if(queryseq_node->selected_flag==false){

//			if(queryseq_node->qname.compare("SRR8858445.1.146248")==0){
//				cout << "i=" << i << ", " << queryseq_node->qname << endl;
//			}

			query_seqs = getQuerySeqs(queryseq_node->qname, query_seq_info_vec);
			qcSigList_node = extractQcSigsSingleQueryClipReg(query_seqs);
			if(qcSigList_node) qcSigList_vec.push_back(qcSigList_node);

			for(j=0; j<query_seqs.size(); j++) query_seqs.at(j)->selected_flag = true;
		}
	}

	for(i=0; i<query_seq_info_vec.size(); i++) query_seq_info_vec.at(i)->selected_flag = false;

	return qcSigList_vec;
}

qcSigList_t* clipRegCluster::extractQcSigsSingleQueryClipReg(vector<struct querySeqInfoNode*> &query_seqs){
	vector<qcSig_t*> qcSig_vec, qcSig_vec_indel;
	qcSigList_t *qcSigList_node;
	vector<int32_t> sidemost_idx_vec, adjClipAlnSegInfo;
	struct querySeqInfoNode *queryseq_node, *mate_queryseq_node;
	qcSig_t *qc_sig, *qc_sig_next, *last_clip_sig;
	int32_t id, leftmost_id, increase_direction;
	int32_t left_clip_end_flag, right_clip_end_flag, end_flag;
	bool overlap_flag, valid_query_left_end_flag, valid_query_right_end_flag, self_overlap_flag, ref_skip_flag, query_skip_flag;
	size_t i, j;
	int32_t mate_arr_idx, mate_clip_end_flag, mate_arr_idx_same_chr, mate_clip_end_flag_same_chr, min_dist, min_ref_dist, min_dist_same_chr, min_ref_dist_same_chr, seq_len, ignore_end_flag;
	int64_t var_startRefPos_tmp, var_endRefPos_tmp, end_ref_pos, dist;
	string reg_str, refseq;
	char *p_seq;
	int32_t seq_id_arr[query_seqs.size()];

	var_startRefPos_tmp = var_startRefPos - GROUP_DIST_THRES;
	var_endRefPos_tmp = var_endRefPos + GROUP_DIST_THRES;
	if(var_startRefPos_tmp<1) var_startRefPos_tmp = 1;
	if(var_endRefPos_tmp>chrlen) var_endRefPos_tmp = chrlen;

	for(i=0; i<query_seqs.size(); i++) seq_id_arr[i] = 0;

	sidemost_idx_vec = getLeftRightMostAlnSegIdx(query_seqs);
	//cout << "leftmost_id=" << sidemost_idx_vec.at(0) << ", rightmost_id=" << sidemost_idx_vec.at(1) << endl;

	leftmost_id = sidemost_idx_vec.at(0);
	if(leftmost_id==-1) leftmost_id = 0;
	id = leftmost_id;
	end_flag = ignore_end_flag = -1;
	mate_arr_idx_same_chr = mate_clip_end_flag_same_chr = min_dist_same_chr = min_ref_dist_same_chr = -1;
	while(1){
		queryseq_node = query_seqs.at(id);
		seq_id_arr[id] = 1;

		// clipping at left end
		last_clip_sig = NULL;
		valid_query_left_end_flag = false;
		left_clip_end_flag = -1;
		if(end_flag!=RIGHT_END){
			if(ignore_end_flag!=LEFT_END){
				end_ref_pos = getEndRefPosAlnSeg(queryseq_node->query_alnSegs.at(0)->startRpos, queryseq_node->query_alnSegs.at(0)->opflag, queryseq_node->query_alnSegs.at(0)->seglen);
				overlap_flag = isOverlappedPos(queryseq_node->query_alnSegs.at(0)->startRpos, end_ref_pos, var_startRefPos_tmp, var_endRefPos_tmp);
				if(overlap_flag and queryseq_node->clip_aln->left_aln and queryseq_node->clip_aln->leftClipSize>=minClipEndSize){
					valid_query_left_end_flag = true;
					left_clip_end_flag = LEFT_END;
				}
				if(valid_query_left_end_flag){ // left end
					qc_sig = new qcSig_t();
					qc_sig->chrname = queryseq_node->clip_aln->chrname;
					qc_sig->cigar_op = queryseq_node->query_alnSegs.at(0)->opflag;
					qc_sig->cigar_op_len = queryseq_node->query_alnSegs.at(0)->seglen;
					qc_sig->aln_orient = queryseq_node->clip_aln->aln_orient;
					qc_sig->reg_contain_flag = false;
					qc_sig->ref_pos = queryseq_node->query_alnSegs.at(0)->startRpos;
					qc_sig->query_pos = queryseq_node->query_alnSegs.at(0)->startQpos;
					qc_sig->chrname_next_clip = "";
					qc_sig->dist_next_clip = 0;
					qc_sig->have_next_clip = false;
					qc_sig->mate_qcSig = NULL;
					qcSig_vec.push_back(qc_sig);
				}
			}

			// indel
			qcSig_vec_indel = extractQcSigsFromAlnSegsSingleQuery(queryseq_node, queryseq_node->clip_aln->chrname, var_startRefPos_tmp, var_endRefPos_tmp, min_sv_size);
			for(i=0; i<qcSig_vec_indel.size(); i++) {
				qc_sig = qcSig_vec_indel.at(i);
				if(qc_sig->cigar_op==BAM_CINS){
					reg_str = qc_sig->chrname + ":" + to_string(qc_sig->ref_pos) + "-" + to_string(qc_sig->ref_pos);
					pthread_mutex_lock(&mutex_fai);
					p_seq = fai_fetch(fai, reg_str.c_str(), &seq_len);
					pthread_mutex_unlock(&mutex_fai);
					qc_sig->refseq = p_seq;
					free(p_seq);
				}else if(qc_sig->cigar_op==BAM_CDEL){
					qc_sig->altseq = queryseq_node->seq.substr(qc_sig->query_pos, 1);
				}
				qcSig_vec.push_back(qc_sig);
			}

			// clipping at right end
			if(ignore_end_flag!=RIGHT_END){
				valid_query_right_end_flag = false;
				right_clip_end_flag = -1;
				end_ref_pos = getEndRefPosAlnSeg(queryseq_node->query_alnSegs.at(queryseq_node->query_alnSegs.size()-1)->startRpos, queryseq_node->query_alnSegs.at(queryseq_node->query_alnSegs.size()-1)->opflag, queryseq_node->query_alnSegs.at(queryseq_node->query_alnSegs.size()-1)->seglen);
				overlap_flag = isOverlappedPos(queryseq_node->query_alnSegs.at(queryseq_node->query_alnSegs.size()-1)->startRpos, end_ref_pos, var_startRefPos_tmp, var_endRefPos_tmp);
				if(overlap_flag and queryseq_node->clip_aln->right_aln and queryseq_node->clip_aln->rightClipSize>=minClipEndSize){
					valid_query_right_end_flag = true;
					right_clip_end_flag = RIGHT_END;
				}
				if(valid_query_right_end_flag){ // right end
					qc_sig = new qcSig_t();
					qc_sig->chrname = queryseq_node->clip_aln->chrname;
					qc_sig->cigar_op = queryseq_node->query_alnSegs.at(queryseq_node->query_alnSegs.size()-1)->opflag;
					qc_sig->cigar_op_len = queryseq_node->query_alnSegs.at(queryseq_node->query_alnSegs.size()-1)->seglen;
					qc_sig->aln_orient = queryseq_node->clip_aln->aln_orient;
					qc_sig->reg_contain_flag = false;
					qc_sig->ref_pos = queryseq_node->query_alnSegs.at(queryseq_node->query_alnSegs.size()-1)->startRpos;
					qc_sig->query_pos = queryseq_node->query_alnSegs.at(queryseq_node->query_alnSegs.size()-1)->startQpos;
					qc_sig->chrname_next_clip = "";
					qc_sig->dist_next_clip = 0;
					qc_sig->have_next_clip = false;
					qc_sig->mate_qcSig = NULL;
					qcSig_vec.push_back(qc_sig);
					last_clip_sig = qc_sig;
				}else break;
			}
		}else{
			// clipping at right end
			valid_query_right_end_flag = false;
			right_clip_end_flag = -1;
			if(end_flag!=LEFT_END){
				if(ignore_end_flag!=RIGHT_END){
					end_ref_pos = getEndRefPosAlnSeg(queryseq_node->query_alnSegs.at(queryseq_node->query_alnSegs.size()-1)->startRpos, queryseq_node->query_alnSegs.at(queryseq_node->query_alnSegs.size()-1)->opflag, queryseq_node->query_alnSegs.at(queryseq_node->query_alnSegs.size()-1)->seglen);
					overlap_flag = isOverlappedPos(queryseq_node->query_alnSegs.at(queryseq_node->query_alnSegs.size()-1)->startRpos, end_ref_pos, var_startRefPos_tmp, var_endRefPos_tmp);
					if(overlap_flag and queryseq_node->clip_aln->right_aln and queryseq_node->clip_aln->rightClipSize>=minClipEndSize){
						valid_query_right_end_flag = true;
						right_clip_end_flag = RIGHT_END;
					}
					if(valid_query_right_end_flag){ // right end
						qc_sig = new qcSig_t();
						qc_sig->chrname = queryseq_node->clip_aln->chrname;
						qc_sig->cigar_op = queryseq_node->query_alnSegs.at(queryseq_node->query_alnSegs.size()-1)->opflag;
						qc_sig->cigar_op_len = queryseq_node->query_alnSegs.at(queryseq_node->query_alnSegs.size()-1)->seglen;
						qc_sig->aln_orient = queryseq_node->clip_aln->aln_orient;
						qc_sig->reg_contain_flag = false;
						qc_sig->ref_pos = queryseq_node->query_alnSegs.at(queryseq_node->query_alnSegs.size()-1)->startRpos;
						qc_sig->query_pos = queryseq_node->query_alnSegs.at(queryseq_node->query_alnSegs.size()-1)->startQpos;
						qc_sig->chrname_next_clip = "";
						qc_sig->dist_next_clip = 0;
						qc_sig->have_next_clip = false;
						qc_sig->mate_qcSig = NULL;
						qcSig_vec.push_back(qc_sig);
					}
				}

				// indel
				qcSig_vec_indel = extractQcSigsFromAlnSegsSingleQuery(queryseq_node, queryseq_node->clip_aln->chrname, var_startRefPos_tmp, var_endRefPos_tmp, min_sv_size);
				for(i=0; i<qcSig_vec_indel.size(); i++) {
					qc_sig = qcSig_vec_indel.at(i);
					if(qc_sig->cigar_op==BAM_CINS){
						reg_str = qc_sig->chrname + ":" + to_string(qc_sig->ref_pos) + "-" + to_string(qc_sig->ref_pos);
						pthread_mutex_lock(&mutex_fai);
						p_seq = fai_fetch(fai, reg_str.c_str(), &seq_len);
						pthread_mutex_unlock(&mutex_fai);
						qc_sig->refseq = p_seq;
						free(p_seq);
					}else if(qc_sig->cigar_op==BAM_CDEL){
						qc_sig->altseq = queryseq_node->seq.substr(qc_sig->query_pos, 1);
					}
					qcSig_vec.push_back(qc_sig);
				}

				// clipping at left end
				if(ignore_end_flag!=LEFT_END){
					valid_query_left_end_flag = false;
					left_clip_end_flag = -1;
					end_ref_pos = getEndRefPosAlnSeg(queryseq_node->query_alnSegs.at(0)->startRpos, queryseq_node->query_alnSegs.at(0)->opflag, queryseq_node->query_alnSegs.at(0)->seglen);
					overlap_flag = isOverlappedPos(queryseq_node->query_alnSegs.at(0)->startRpos, end_ref_pos, var_startRefPos_tmp, var_endRefPos_tmp);
					if(overlap_flag and queryseq_node->clip_aln->left_aln and queryseq_node->clip_aln->leftClipSize>=minClipEndSize){
						valid_query_left_end_flag = true;
						left_clip_end_flag = LEFT_END;
					}
					if(valid_query_left_end_flag){ // left end
						qc_sig = new qcSig_t();
						qc_sig->chrname = queryseq_node->clip_aln->chrname;
						qc_sig->cigar_op = queryseq_node->query_alnSegs.at(0)->opflag;
						qc_sig->cigar_op_len = queryseq_node->query_alnSegs.at(0)->seglen;
						qc_sig->aln_orient = queryseq_node->clip_aln->aln_orient;
						qc_sig->reg_contain_flag = false;
						qc_sig->ref_pos = queryseq_node->query_alnSegs.at(0)->startRpos;
						qc_sig->query_pos = queryseq_node->query_alnSegs.at(0)->startQpos;
						qc_sig->chrname_next_clip = "";
						qc_sig->dist_next_clip = 0;
						qc_sig->have_next_clip = false;
						qc_sig->mate_qcSig = NULL;
						qcSig_vec.push_back(qc_sig);
						last_clip_sig = qc_sig;
					}else break;
				}
			}
		}

		if(queryseq_node->clip_aln->rightmost_flag) break;
		else{
			ignore_end_flag = -1;
			adjClipAlnSegInfo.clear();
			mate_arr_idx = -1;
			mate_clip_end_flag = -1;
			if(end_flag==-1){ // the first segment
				if(valid_query_right_end_flag){
					adjClipAlnSegInfo = getAdjacentAlnSegClipReg(id, right_clip_end_flag, query_seqs, minClipEndSize);
					mate_arr_idx = adjClipAlnSegInfo.at(0);
					mate_clip_end_flag = adjClipAlnSegInfo.at(1);
				}else if(valid_query_left_end_flag){
					adjClipAlnSegInfo = getAdjacentAlnSegClipReg(id, left_clip_end_flag, query_seqs, minClipEndSize);
					mate_arr_idx = adjClipAlnSegInfo.at(0);
					mate_clip_end_flag = adjClipAlnSegInfo.at(1);
				}
			}else{ // not the first segment
				if(end_flag==LEFT_END){
					adjClipAlnSegInfo = getAdjacentAlnSegClipReg(id, RIGHT_END, query_seqs, minClipEndSize);
					mate_arr_idx = adjClipAlnSegInfo.at(0);
					mate_clip_end_flag = adjClipAlnSegInfo.at(1);
				}else if(end_flag==RIGHT_END){
					adjClipAlnSegInfo = getAdjacentAlnSegClipReg(id, LEFT_END, query_seqs, minClipEndSize);
					mate_arr_idx = adjClipAlnSegInfo.at(0);
					mate_clip_end_flag = adjClipAlnSegInfo.at(1);
				}
			}
			if(mate_arr_idx!=-1 and mate_clip_end_flag!=-1){
				ref_skip_flag = query_skip_flag = false;
				min_dist = adjClipAlnSegInfo.at(2);
				min_ref_dist = adjClipAlnSegInfo.at(3);
				mate_arr_idx_same_chr = adjClipAlnSegInfo.at(4);
				mate_clip_end_flag_same_chr = adjClipAlnSegInfo.at(5);
				min_dist_same_chr = adjClipAlnSegInfo.at(6);
				min_ref_dist_same_chr = adjClipAlnSegInfo.at(7);
				increase_direction = adjClipAlnSegInfo.at(8);
				if(mate_arr_idx_same_chr!=-1){
					mate_queryseq_node = query_seqs.at(mate_arr_idx_same_chr);
					if((mate_clip_end_flag_same_chr==LEFT_END and mate_queryseq_node->clip_aln->left_aln) or (mate_clip_end_flag_same_chr==RIGHT_END and mate_queryseq_node->clip_aln->right_aln)){
						self_overlap_flag = isSegSelfOverlap(queryseq_node->clip_aln, mate_queryseq_node->clip_aln, MAX_VAR_REG_SIZE);
						dist = min_dist_same_chr - min_ref_dist_same_chr;
						if(increase_direction==1) dist = -dist;
						// prefer the segments on the same chromosome
						//if((mate_arr_idx!=mate_arr_idx_same_chr and (abs(min_ref_dist_same_chr)<=MAX_REF_DIST_SAME_CHR or self_overlap_flag)) or (mate_arr_idx==mate_arr_idx_same_chr and abs(min_dist_same_chr)>MAX_REF_DIST_SAME_CHR and abs(dist)>MAX_REF_DIST_SAME_CHR)){ // different chromosomes with near reference distance, deleted on 2024-02-25
						if((abs(min_ref_dist_same_chr)<=MAX_REF_DIST_SAME_CHR or self_overlap_flag) or (mate_arr_idx==mate_arr_idx_same_chr and abs(min_dist_same_chr)>MAX_REF_DIST_SAME_CHR and abs(dist)>MAX_REF_DIST_SAME_CHR)){ // different chromosomes with near reference distance
							mate_arr_idx = mate_arr_idx_same_chr;
							mate_clip_end_flag = mate_clip_end_flag_same_chr;
							query_skip_flag = true;
						}else if(mate_arr_idx==mate_arr_idx_same_chr and mate_clip_end_flag==mate_clip_end_flag_same_chr and self_overlap_flag==false and abs(min_dist_same_chr)<=MAX_REF_DIST_SAME_CHR and abs(min_ref_dist_same_chr)>MAX_REF_DIST_SAME_CHR){
							mate_arr_idx = mate_arr_idx_same_chr;
							mate_clip_end_flag = mate_clip_end_flag_same_chr;
							ref_skip_flag = true;
						}else{
							if(abs(min_dist_same_chr)>=MAX_REF_DIST_SAME_CHR and abs(min_dist_same_chr)<=MAX_VAR_REG_SIZE and abs(min_ref_dist_same_chr)>MAX_REF_DIST_SAME_CHR and abs(min_ref_dist_same_chr)<=MAX_VAR_REG_SIZE){
								if(queryseq_node->clip_aln->aln_orient==mate_queryseq_node->clip_aln->aln_orient and abs(min_dist)>MAX_REF_DIST_SAME_CHR and abs(min_ref_dist)>MAX_REF_DIST_SAME_CHR){
									if(dist>0){
										mate_arr_idx = mate_arr_idx_same_chr;
										mate_clip_end_flag = mate_clip_end_flag_same_chr;
										query_skip_flag = true;
									}else{
										mate_arr_idx = mate_arr_idx_same_chr;
										mate_clip_end_flag = mate_clip_end_flag_same_chr;
										ref_skip_flag = true;
									}
								}
							}
						}
					}
				}

				if(last_clip_sig){
					if(query_skip_flag){
						last_clip_sig->dist_next_clip = last_clip_sig->cigar_op_len = abs(dist) + 1;
						last_clip_sig->cigar_op = BAM_CINS;
						//end_ref_pos = getEndRefPosAlnSeg(queryseq_node->query_alnSegs.at(0)->startRpos, queryseq_node->query_alnSegs.at(0)->opflag, queryseq_node->query_alnSegs.at(0)->seglen);
						reg_str = last_clip_sig->chrname + ":" + to_string(last_clip_sig->ref_pos) + "-" + to_string(last_clip_sig->ref_pos);
						pthread_mutex_lock(&mutex_fai);
						p_seq = fai_fetch(fai, reg_str.c_str(), &seq_len);
						pthread_mutex_unlock(&mutex_fai);
						last_clip_sig->refseq = p_seq;
						free(p_seq);

						last_clip_sig->altseq = queryseq_node->seq.substr(last_clip_sig->query_pos-1, last_clip_sig->cigar_op_len);

					}else if(ref_skip_flag){
						last_clip_sig->dist_next_clip = last_clip_sig->cigar_op_len = abs(dist) + 1;
						last_clip_sig->cigar_op = BAM_CDEL;

						reg_str = last_clip_sig->chrname + ":" + to_string(last_clip_sig->ref_pos) + "-" + to_string(last_clip_sig->ref_pos+last_clip_sig->cigar_op_len-1);
						pthread_mutex_lock(&mutex_fai);
						p_seq = fai_fetch(fai, reg_str.c_str(), &seq_len);
						pthread_mutex_unlock(&mutex_fai);
						last_clip_sig->refseq = p_seq;
						free(p_seq);

						last_clip_sig->altseq = queryseq_node->seq.substr(last_clip_sig->query_pos-1, 1);
					}
				}else{
					query_skip_flag = ref_skip_flag = false;
				}
				//for(i=0; i<adjClipAlnSegInfo.size(); i++) last_clip_sig->adjClipAlnSegInfo.push_back(adjClipAlnSegInfo.at(i));

				id = mate_arr_idx;
				end_flag = mate_clip_end_flag;
				if(query_skip_flag or ref_skip_flag) ignore_end_flag = mate_clip_end_flag;

				if(seq_id_arr[id]==1) break;
			}else
				break;
		}
	}

	// assign neighbour information
	for(i=1; i<qcSig_vec.size(); i++){
		qc_sig = qcSig_vec.at(i-1);
		if(qc_sig->cigar_op==BAM_CSOFT_CLIP or qc_sig->cigar_op==BAM_CHARD_CLIP){ // clipping
			// find the next clip
			for(j=i; j<qcSig_vec.size(); j++){
				qc_sig_next = qcSig_vec.at(j);
				if(qc_sig_next->cigar_op==BAM_CSOFT_CLIP or qc_sig_next->cigar_op==BAM_CHARD_CLIP){ // clipping
					qc_sig->chrname_next_clip = qc_sig_next->chrname;
					qc_sig->dist_next_clip = qc_sig_next->ref_pos - qc_sig->ref_pos;
					qc_sig->have_next_clip = true;
					i = j;
					break;
				}
			}
		}
	}

	qcSigList_node = new qcSigList_t();
	qcSigList_node->startSpanRefPos = sidemost_idx_vec.at(2);
	qcSigList_node->endSpanRefPos = sidemost_idx_vec.at(3);
	qcSigList_node->query_seqs_vec = query_seqs;
	qcSigList_node->qcSig_vec = qcSig_vec;

	return qcSigList_node;
}

// get the left and right most align segments
vector<int32_t> clipRegCluster::getLeftRightMostAlnSegIdx(vector<struct querySeqInfoNode*> &query_seqs){
	vector<int32_t> ret_idx_vec;
	size_t i;
	int32_t leftmost_id, rightmost_id;
	int64_t startSpanRefPos, endSpanRefPos, maxValue, minValue;
	struct querySeqInfoNode *query_seq_node;

	leftmost_id = rightmost_id = -1;
	startSpanRefPos = endSpanRefPos = -1;
	minValue = LONG_MAX;
	maxValue = LONG_MIN;
	for(i=0; i<query_seqs.size(); i++){
		query_seq_node = query_seqs.at(i);
		if(query_seq_node->clip_aln->leftmost_flag) leftmost_id = i;
		if(query_seq_node->clip_aln->rightmost_flag) rightmost_id = i;
		if(query_seq_node->clip_aln->startRefPos<minValue) minValue = query_seq_node->clip_aln->startRefPos;
		if(query_seq_node->clip_aln->endRefPos>maxValue) maxValue = query_seq_node->clip_aln->endRefPos;
	}

	if(minValue!=LONG_MAX and maxValue!=LONG_MIN){
		startSpanRefPos = minValue;
		endSpanRefPos = maxValue;
	}else{
		cout << "line=" << __LINE__ << ", qname=" << query_seq_node->qname << ", minValue=" << minValue << ", maxValue=" << maxValue << ", error!" << endl;
		exit(1);
	}

//	if(leftmost_id==-1 or rightmost_id==-1){
//		cout << "line=" << __LINE__ << ", qname=" << query_seq_node->qname << ", leftmost_id=" << leftmost_id << ", rightmost_id=" << rightmost_id << ", can not find the left and right most align segments, error!" << endl;
//		//exit(1);
//	}

	ret_idx_vec.push_back(leftmost_id);
	ret_idx_vec.push_back(rightmost_id);
	ret_idx_vec.push_back(startSpanRefPos);
	ret_idx_vec.push_back(endSpanRefPos);

	return ret_idx_vec;
}

// get adjacent clip segment according to query positions
vector<int32_t> clipRegCluster::getAdjacentAlnSegClipReg(int32_t arr_idx, int32_t clip_end_flag, vector<struct querySeqInfoNode*> &query_seqs, int32_t minClipEndSize){
	vector<int32_t> adjClipAlnSegInfo; // [0]: array index of minimal distance; [1]: segment end flag
	vector<clipAlnData_t*> query_aln_segs;
	size_t i;

	for(i=0; i<query_seqs.size(); i++) if(query_seqs.at(i)->clip_aln) query_aln_segs.push_back(query_seqs.at(i)->clip_aln);
	adjClipAlnSegInfo = getAdjacentClipAlnSeg(arr_idx, clip_end_flag, query_aln_segs, minClipEndSize, MAX_VAR_REG_SIZE);

	return adjClipAlnSegInfo;
}

struct seedQueryInfo* clipRegCluster::chooseSeedClusterQueryClipReg(qcSigList_t* qcSigList_node, vector<qcSigList_t*> &q_cluster){
	size_t seed_id, j;
	int64_t startSpanPos_tmp, endSpanPos_tmp, startSpanPos, endSpanPos, span, span_max;

	span_max = 0, seed_id = 0;
	startSpanPos = 0;
	endSpanPos = 0;
	struct seedQueryInfo* seed_info;
	//if(query_seq_info_node->query_alnSegs.size()>0){
	for(j=0;j<q_cluster.size();j++){
		startSpanPos_tmp = q_cluster.at(j)->startSpanRefPos;
		endSpanPos_tmp = q_cluster.at(j)->endSpanRefPos;
		if(startSpanPos_tmp < endSpanPos_tmp){
			span = endSpanPos_tmp - startSpanPos_tmp;
			if(span>span_max){
				span_max = span;
				startSpanPos = startSpanPos_tmp;
				endSpanPos = endSpanPos_tmp;
				seed_id = j;
			}
		}
		//}
	}

	seed_info = new struct seedQueryInfo();
	seed_info->startSpanPos = startSpanPos;
	seed_info->endSpanPos = endSpanPos;
	seed_info->span_intersection = span_max;
	seed_info->id = seed_id;
	return seed_info;
}

double clipRegCluster::computeMatchRatioClipReg(qcSigList_t* qcSigList_node, qcSigList_t* q_cluster_node, int64_t startSpanPos, int64_t endSpanPos){
	double match_ratio = 0;
	size_t i;
	int32_t score_sum;

#if CLIP_REG_CLUSTER_DEBUG
	qcSig_t *qc_sig;
	cout << "----------" << qcSigList_node->query_seqs_vec.at(0)->qname << ", " << q_cluster_node->query_seqs_vec.at(0)->qname << endl;
	for(size_t j=0; j<qcSigList_node->qcSig_vec.size(); j++){
		qc_sig = qcSigList_node->qcSig_vec.at(j);
		if(qc_sig->cigar_op==BAM_CSOFT_CLIP or qc_sig->cigar_op==BAM_CHARD_CLIP) // clippings
			cout << "\t(" << qc_sig->cigar_op << "," << qc_sig->chrname << ":" << qc_sig->ref_pos << "," << qc_sig->dist_next_clip << ")";
		else // indels
			cout << "\t(" << qc_sig->cigar_op << "," << qc_sig->chrname << ":" << qc_sig->ref_pos << "," << qc_sig->cigar_op_len << ")";
	}
	cout << endl;
	for(size_t j=0; j<q_cluster_node->qcSig_vec.size(); j++){
		qc_sig = q_cluster_node->qcSig_vec.at(j);
		if(qc_sig->cigar_op==BAM_CSOFT_CLIP or qc_sig->cigar_op==BAM_CHARD_CLIP) // clippings
			cout << "\t(" << qc_sig->cigar_op << "," << qc_sig->chrname << ":" << qc_sig->ref_pos << "," << qc_sig->dist_next_clip << ")";
		else // indels
			cout << "\t(" << qc_sig->cigar_op << "," << qc_sig->chrname << ":" << qc_sig->ref_pos << "," << qc_sig->cigar_op_len << ")";
	}
	cout << endl;
#endif

	//matching
	qcSigList_node->match_profile_vec = computeQcMatchProfileSingleQueryClipReg(qcSigList_node, q_cluster_node);
	if(qcSigList_node->match_profile_vec.size()==0){
		//score_ratio = 2;//have no varsig
		match_ratio = 0; //0.99
	}else{
		score_sum = 0;
		for(i=0; i<qcSigList_node->match_profile_vec.size(); i++){
			if(qcSigList_node->match_profile_vec.at(i)==1)
				score_sum += qcSigList_node->match_profile_vec.at(i);
		}
		match_ratio = (double)score_sum / (qcSigList_node->match_profile_vec.size());

		// skip head and tail gaps
//		skip_head_num = 0;
//		for(i=0; i<qcSigList_node->match_profile_vec.size(); i++){
//			if(qcSigList_node->match_profile_vec.at(i)==2) skip_head_num ++;
//			else break;
//		}
//		skip_tail_num = 0;
//		for(j=(int32_t)qcSigList_node->match_profile_vec.size()-1; j>=0; j--){
//			if(qcSigList_node->match_profile_vec.at(j)==2) skip_tail_num ++;
//			else break;
//		}
//
//		if(skip_head_num+skip_tail_num<(int32_t)qcSigList_node->match_profile_vec.size()){
//			score_sum = 0;
//			for(i=skip_head_num; i<qcSigList_node->match_profile_vec.size()-skip_tail_num; i++){
//				if(qcSigList_node->match_profile_vec.at(i)==1)
//					score_sum += qcSigList_node->match_profile_vec.at(i);
//			}
//			match_ratio = (double)score_sum / (qcSigList_node->match_profile_vec.size() - skip_head_num - skip_tail_num);
//		}
	}

#if CLIP_REG_CLUSTER_DEBUG
	cout << "match_ratio=" << match_ratio << ", total=" << qcSigList_node->match_profile_vec.size() << endl;
	//cout << "match_ratio=" << match_ratio << ", skip_head_num=" << skip_head_num << ", skip_tail_num=" << skip_tail_num << ", total=" << qcSigList_node->match_profile_vec.size() << endl;
#endif

	return match_ratio;
}

// return values: 0-mismatch, 1-match, 2-gap
vector<int8_t> clipRegCluster::computeQcMatchProfileSingleQueryClipReg(qcSigList_t *queryCluSig, qcSigList_t *seed_qcQuery){
	vector<int8_t> match_profile_vec;
	int32_t rowsNum, colsNum, matchScore, mismatchScore, gapScore, gapOpenScore;
	int32_t scoreIJ, tmp_gapScore1, tmp_gapScore2, maxValue, path_val, maxValue_ul, maxValue_l, maxValue_u;
	struct alnScoreNode *scoreArr;
	int64_t i, j, arrSize;
	qcSig_t *qc_sig1, *qc_sig2;
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
		for(j=1; j<colsNum; j++) scoreArr[j].path_val = 2;
		for(i=1; i<rowsNum; i++) scoreArr[i*colsNum].path_val = 1;

		// compute the scores of each element
		for(i=1; i<rowsNum; i++){
			for(j=1; j<colsNum; j++){
				qc_sig1 = queryCluSig->qcSig_vec.at(i-1);
				qc_sig2 = seed_qcQuery->qcSig_vec.at(j-1);
				matchFlag = isQcSigMatchClipReg(qc_sig1, qc_sig2, MAX_DIST_MATCH_CLIP_POS, MIN_SIZE_RATIO_MATCH_CLIP_POS, QC_SIZE_RATIO_MATCH_THRES_INDEL, QC_CONSIST_RATIO_MATCH_THRES);
				if(matchFlag) scoreIJ = matchScore;
				else scoreIJ = mismatchScore;//

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

		// compute signature match profile
		match_profile_vec = qComputeSigMatchProfileClipReg(scoreArr, rowsNum, colsNum, queryCluSig, seed_qcQuery);
		free(scoreArr);
	}

	return match_profile_vec;
}

bool clipRegCluster::isQcSigMatchClipReg(qcSig_t *qc_sig, qcSig_t *seed_qc_sig, int64_t max_dist_match_clip_pos, double min_size_ratio_match_thres_clip, double size_ratio_match_thres, double consist_ratio_match_thres){
	bool match_flag = false;
	double size_ratio;
	int64_t ref_dist;

	if(qc_sig and seed_qc_sig){
		if((qc_sig->cigar_op==BAM_CINS and seed_qc_sig->cigar_op==BAM_CINS) or (qc_sig->cigar_op==BAM_CDEL and seed_qc_sig->cigar_op==BAM_CDEL)){ // indel
			match_flag = isQcSigMatch(qc_sig, seed_qc_sig, max_dist_match_clip_pos, size_ratio_match_thres, consist_ratio_match_thres, fai);
		}else if((qc_sig->cigar_op==BAM_CSOFT_CLIP or qc_sig->cigar_op==BAM_CHARD_CLIP) and (seed_qc_sig->cigar_op==BAM_CSOFT_CLIP or seed_qc_sig->cigar_op==BAM_CHARD_CLIP)){ // clipping
			if(qc_sig->chrname.compare(seed_qc_sig->chrname)==0){
				ref_dist = abs(qc_sig->ref_pos - seed_qc_sig->ref_pos);
				//cout << "ref_dist=" << ref_dist << endl;
				if(qc_sig->have_next_clip and seed_qc_sig->have_next_clip){ // pos and distance
					size_ratio = 0;
					if(qc_sig->dist_next_clip>0 and seed_qc_sig->dist_next_clip>0){
						if(qc_sig->dist_next_clip<seed_qc_sig->dist_next_clip) size_ratio = (double)qc_sig->dist_next_clip / seed_qc_sig->dist_next_clip;
						else size_ratio = (double)seed_qc_sig->dist_next_clip / qc_sig->dist_next_clip;
					}else if(qc_sig->dist_next_clip<0 and seed_qc_sig->dist_next_clip<0){
						if(qc_sig->dist_next_clip>seed_qc_sig->dist_next_clip) size_ratio = (double)qc_sig->dist_next_clip / seed_qc_sig->dist_next_clip;
						else size_ratio = (double)seed_qc_sig->dist_next_clip / qc_sig->dist_next_clip;
					}
					//cout << "qc_sig->dist_next_clip=" << qc_sig->dist_next_clip << ", seed_qc_sig->dist_next_clip=" << seed_qc_sig->dist_next_clip << ", size_ratio=" << size_ratio << endl;
					if(ref_dist<=max_dist_match_clip_pos and size_ratio>=min_size_ratio_match_thres_clip){
						match_flag = true;
					}
				}else{ // pos
					if(ref_dist<=max_dist_match_clip_pos){
						match_flag = true;
					}
				}
			}
		}
	}

	return match_flag;
}

// return values: 0-mismatch, 1-match, 2-gap
vector<int8_t> clipRegCluster::qComputeSigMatchProfileClipReg(struct alnScoreNode *scoreArr, int32_t rowsNum, int32_t colsNum, qcSigList_t *queryCluSig, qcSigList_t *seed_qcQuery){
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
				//match_profile_vec.push_back(0);
				match_profile_vec.push_back(2);
				i--;
			}
			else{ //from (i, j-1)
				//match_profile_vec.push_back(0);
				match_profile_vec.push_back(2);
				j--;
			}
			if(i==0 and j==0)
				break;
		}

		reverse(match_profile_vec.begin(), match_profile_vec.end());
	}

//============================================
#if CLIP_REG_CLUSTER_DEBUG
	cout << "------------------- Profile -------------------" << endl;
	// output the profile
	vector<string> sigseq_vec1, sigseq_vec2;
	string sig_str1, sig_str2;
	i = rowsNum - 1; j = colsNum - 1;
	while(i>0 or j>0){
		value = scoreArr[i*colsNum+j].path_val;
		if(value==0){ //from (i-1, j-1)
			sig_str1 = to_string(queryCluSig->qcSig_vec.at(i-1)->ref_pos) + "_" + to_string(queryCluSig->qcSig_vec.at(i-1)->dist_next_clip);
			sig_str2 = to_string(seed_qcQuery->qcSig_vec.at(j-1)->ref_pos) + "_" + to_string(seed_qcQuery->qcSig_vec.at(j-1)->dist_next_clip);

			sigseq_vec1.push_back(sig_str1);
			sigseq_vec2.push_back(sig_str2);

			i--;
			j--;
		}
		else if(value==1){ //from (i-1, j)
			//?
			//match_profile_vec.push_back(0);
			sig_str1 = to_string(queryCluSig->qcSig_vec.at(i-1)->ref_pos) + "_" + to_string(queryCluSig->qcSig_vec.at(i-1)->dist_next_clip);
			sigseq_vec1.push_back(sig_str1);
			sigseq_vec2.push_back("-");

			i--;
		}
		else{ //from (i, j-1)
			//match_profile_vec.push_back(0);
			sig_str2 = to_string(seed_qcQuery->qcSig_vec.at(j-1)->ref_pos) + "_" + to_string(seed_qcQuery->qcSig_vec.at(j-1)->dist_next_clip);
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

