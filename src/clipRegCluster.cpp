#include "clipRegCluster.h"

clipRegCluster::clipRegCluster(int32_t minClipEndSize, int32_t min_sv_size) {
	this->minClipEndSize = minClipEndSize;
	this->min_sv_size = min_sv_size;
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
	qcSigListVec_t *q_cluster_node_a, *q_cluster_node_b;
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
	int32_t mrsame_maxnum;

	// extract the features
	qcSigList_vec = extractQcSigsClipReg(query_seq_info_vec);

	for(i=0; i<qcSigList_vec.size(); i++) qcSigList_vec.at(i)->cluster_finished_flag = false;

	// sort in ascending order according to the number of features
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
		match_min = 0.99;
		match_ratio_vec.clear();
		min_id_vec.clear();
		for(i=0; i<qcSigList_vec.size(); i++){
			if(qcSigList_vec.at(i)->qcSig_vec.size()>0 and qcSigList_vec.at(i)->cluster_finished_flag==false){
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
			}else{
				//span_vec1.push_back(0);
				match_ratio_vec.push_back(1);
			}
		}

		exist_flag = false;
		if(match_min < 0.99){
			for(i=0; i<match_ratio_vec.size(); i++){
				if(match_ratio_vec.at(i)<0.99)
					min_id_vec.push_back(i);
			}
			for(i=0; i<min_id_vec.size(); i++){
				for(j=0; j<mrmin_match_vec.size(); j++){
					if(match_ratio_vec.at(min_id_vec.at(i))==mrmin_match_vec.at(j)->SR){
						smin_id = j;
						exist_flag = true;
					}else{
						exist_flag = false;
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

			mrsame_maxnum = mrmin_match_vec.at(0)->num;
			smin_id = 0;
			for(i=0; i<mrmin_match_vec.size(); i++){
				if(mrsame_maxnum<mrmin_match_vec.at(i)->num){
					smin_id = i;
					mrsame_maxnum = mrmin_match_vec.at(i)->num;
				}
			}

			qcSigList_vec.at(mrmin_match_vec.at(smin_id)->pos_id.at(0))->cluster_finished_flag = true;
			q_cluster_b.push_back(qcSigList_vec.at(mrmin_match_vec.at(smin_id)->pos_id.at(0)));
		}

		if(mrmin_match_vec.size()>0){
			//for(i=0; i<mrmin_match_vec.size(); i++){
				vector<srmin_match_t*>::iterator svp;
				for(svp=mrmin_match_vec.begin(); svp!=mrmin_match_vec.end(); svp++) delete *svp;
				vector<srmin_match_t*>().swap(mrmin_match_vec);
			//}
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
						if(match_ratio_a<0.99){
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
						if(match_ratio_b<0.99){
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

//	cout << "*********** Total queries, size=" << qcSigList_vec.size() << ":" << endl;
//	printQcSigListVec(qcSigList_vec);

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
	struct querySeqInfoNode *queryseq_node;
	qcSig_t *qc_sig, *qc_sig_next;
	int32_t id, neighbor_id, leftmost_id, neighbor_clip_end_flag;
	int32_t left_clip_end_flag, right_clip_end_flag, end_flag;
	bool valid_query_left_end_flag, valid_query_right_end_flag;
	size_t i, j;

	sidemost_idx_vec = getLeftRightMostAlnSegIdx(query_seqs);
	//cout << "leftmost_id=" << sidemost_idx_vec.at(0) << ", rightmost_id=" << sidemost_idx_vec.at(1) << endl;

	leftmost_id = sidemost_idx_vec.at(0);
	if(leftmost_id==-1) leftmost_id = 0;
	id = leftmost_id;
	end_flag = -1;
	while(1){
		queryseq_node = query_seqs.at(id);

		// clipping at left end
		valid_query_left_end_flag = false;
		left_clip_end_flag = -1;
		if(end_flag!=RIGHT_END){
			if(queryseq_node->clip_aln->left_aln and queryseq_node->clip_aln->leftClipSize>=minClipEndSize){
				valid_query_left_end_flag = true;
				left_clip_end_flag = LEFT_END;
			}
			if(valid_query_left_end_flag){ // left end
				qc_sig = new qcSig_t();
				qc_sig->chrname = queryseq_node->clip_aln->chrname;
				qc_sig->cigar_op = queryseq_node->query_alnSegs.at(0)->opflag;
				qc_sig->cigar_op_len = queryseq_node->query_alnSegs.at(0)->seglen;
				qc_sig->reg_contain_flag = false;
				qc_sig->ref_pos = queryseq_node->query_alnSegs.at(0)->startRpos;
				qc_sig->chrname_next_clip = "";
				qc_sig->dist_next_clip = 0;
				qc_sig->have_next_clip = false;
				qc_sig->mate_qcSig = NULL;
				qcSig_vec.push_back(qc_sig);
			}

			// indel
			qcSig_vec_indel = extractQcSigsFromAlnSegsSingleQuery(queryseq_node, queryseq_node->clip_aln->chrname, queryseq_node->clip_aln->startRefPos, queryseq_node->clip_aln->endRefPos, min_sv_size);
			for(i=0; i<qcSig_vec_indel.size(); i++) qcSig_vec.push_back(qcSig_vec_indel.at(i));

			// clipping at right end
			valid_query_right_end_flag = false;
			right_clip_end_flag = -1;
			if(queryseq_node->clip_aln->right_aln and queryseq_node->clip_aln->rightClipSize>=minClipEndSize){
				valid_query_right_end_flag = true;
				right_clip_end_flag = RIGHT_END;
			}
			if(valid_query_right_end_flag){ // right end
				qc_sig = new qcSig_t();
				qc_sig->chrname = queryseq_node->clip_aln->chrname;
				qc_sig->cigar_op = queryseq_node->query_alnSegs.at(queryseq_node->query_alnSegs.size()-1)->opflag;
				qc_sig->cigar_op_len = queryseq_node->query_alnSegs.at(queryseq_node->query_alnSegs.size()-1)->seglen;
				qc_sig->reg_contain_flag = false;
				qc_sig->ref_pos = queryseq_node->query_alnSegs.at(queryseq_node->query_alnSegs.size()-1)->startRpos;
				qc_sig->chrname_next_clip = "";
				qc_sig->dist_next_clip = 0;
				qc_sig->have_next_clip = false;
				qc_sig->mate_qcSig = NULL;
				qcSig_vec.push_back(qc_sig);
			}
		}else{
			// clipping at right end
			valid_query_right_end_flag = false;
			right_clip_end_flag = -1;
			if(end_flag!=LEFT_END){
				if(queryseq_node->clip_aln->right_aln and queryseq_node->clip_aln->rightClipSize>=minClipEndSize){
					valid_query_right_end_flag = true;
					right_clip_end_flag = RIGHT_END;
				}
				if(valid_query_right_end_flag){ // right end
					qc_sig = new qcSig_t();
					qc_sig->chrname = queryseq_node->clip_aln->chrname;
					qc_sig->cigar_op = queryseq_node->query_alnSegs.at(queryseq_node->query_alnSegs.size()-1)->opflag;
					qc_sig->cigar_op_len = queryseq_node->query_alnSegs.at(queryseq_node->query_alnSegs.size()-1)->seglen;
					qc_sig->reg_contain_flag = false;
					qc_sig->ref_pos = queryseq_node->query_alnSegs.at(queryseq_node->query_alnSegs.size()-1)->startRpos;
					qc_sig->chrname_next_clip = "";
					qc_sig->dist_next_clip = 0;
					qc_sig->have_next_clip = false;
					qc_sig->mate_qcSig = NULL;
					qcSig_vec.push_back(qc_sig);
				}

				// indel
				qcSig_vec_indel = extractQcSigsFromAlnSegsSingleQuery(queryseq_node, queryseq_node->clip_aln->chrname, queryseq_node->clip_aln->startRefPos, queryseq_node->clip_aln->endRefPos, min_sv_size);
				for(i=0; i<qcSig_vec_indel.size(); i++) qcSig_vec.push_back(qcSig_vec_indel.at(i));

				// clipping at left end
				valid_query_left_end_flag = false;
				left_clip_end_flag = -1;
				if(queryseq_node->clip_aln->left_aln and queryseq_node->clip_aln->leftClipSize>=minClipEndSize){
					valid_query_left_end_flag = true;
					left_clip_end_flag = LEFT_END;
				}
				if(valid_query_left_end_flag){ // left end
					qc_sig = new qcSig_t();
					qc_sig->chrname = queryseq_node->clip_aln->chrname;
					qc_sig->cigar_op = queryseq_node->query_alnSegs.at(0)->opflag;
					qc_sig->cigar_op_len = queryseq_node->query_alnSegs.at(0)->seglen;
					qc_sig->reg_contain_flag = false;
					qc_sig->ref_pos = queryseq_node->query_alnSegs.at(0)->startRpos;
					qc_sig->chrname_next_clip = "";
					qc_sig->dist_next_clip = 0;
					qc_sig->have_next_clip = false;
					qc_sig->mate_qcSig = NULL;
					qcSig_vec.push_back(qc_sig);
				}
			}
		}

		if(queryseq_node->clip_aln->rightmost_flag) break;
		else{
			neighbor_id = -1;
			neighbor_clip_end_flag = -1;
			if(end_flag==-1){ // the first segment
				if(valid_query_right_end_flag){
					adjClipAlnSegInfo = getAdjacentAlnSegClipReg(id, right_clip_end_flag, query_seqs, minClipEndSize);
					neighbor_id = adjClipAlnSegInfo.at(0);
					neighbor_clip_end_flag = adjClipAlnSegInfo.at(1);
				}else if(valid_query_left_end_flag){
					adjClipAlnSegInfo = getAdjacentAlnSegClipReg(id, left_clip_end_flag, query_seqs, minClipEndSize);
					neighbor_id = adjClipAlnSegInfo.at(0);
					neighbor_clip_end_flag = adjClipAlnSegInfo.at(1);
				}
			}else{ // not the first segment
				if(end_flag==LEFT_END){
					adjClipAlnSegInfo = getAdjacentAlnSegClipReg(id, RIGHT_END, query_seqs, minClipEndSize);
					neighbor_id = adjClipAlnSegInfo.at(0);
					neighbor_clip_end_flag = adjClipAlnSegInfo.at(1);
				}else if(end_flag==RIGHT_END){
					adjClipAlnSegInfo = getAdjacentAlnSegClipReg(id, LEFT_END, query_seqs, minClipEndSize);
					neighbor_id = adjClipAlnSegInfo.at(0);
					neighbor_clip_end_flag = adjClipAlnSegInfo.at(1);
				}
			}
			if(neighbor_id!=-1 and neighbor_clip_end_flag!=-1){
				id = neighbor_id;
				end_flag = neighbor_clip_end_flag;
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
	int32_t j, score_sum, skip_head_num, skip_tail_num;

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
		match_ratio = 0.99;
	}else{
		// skip head and tail gaps
		skip_head_num = 0;
		for(i=0; i<qcSigList_node->match_profile_vec.size(); i++){
			if(qcSigList_node->match_profile_vec.at(i)==2) skip_head_num ++;
			else break;
		}
		skip_tail_num = 0;
		for(j=(int32_t)qcSigList_node->match_profile_vec.size()-1; j>=0; j--){
			if(qcSigList_node->match_profile_vec.at(j)==2) skip_tail_num ++;
			else break;
		}

		if(skip_head_num+skip_tail_num<(int32_t)qcSigList_node->match_profile_vec.size()){
			score_sum = 0;
			for(i=skip_head_num; i<qcSigList_node->match_profile_vec.size()-skip_tail_num; i++){
				if(qcSigList_node->match_profile_vec.at(i)==1)
					score_sum += qcSigList_node->match_profile_vec.at(i);
			}
			match_ratio = (double)score_sum / (qcSigList_node->match_profile_vec.size() - skip_head_num - skip_tail_num);
		}
	}

#if CLIP_REG_CLUSTER_DEBUG
	cout << "match_ratio=" << match_ratio << endl;
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
				matchFlag = isQcSigMatchClipReg(qc_sig1, qc_sig2, MAX_DIST_MATCH_CLIP_POS, MIN_SIZE_RATIO_MATCH_CLIP_POS, QC_SIZE_RATIO_MATCH_THRES_INDEL);
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

bool clipRegCluster::isQcSigMatchClipReg(qcSig_t *qc_sig, qcSig_t *seed_qc_sig, int64_t max_dist_match_clip_pos, double min_size_ratio_match_thres_clip, double size_ratio_match_thres){
	bool match_flag = false;
	double size_ratio;
	int64_t ref_dist;

	if(qc_sig and seed_qc_sig){
		if((qc_sig->cigar_op==BAM_CINS and seed_qc_sig->cigar_op==BAM_CINS) or (qc_sig->cigar_op==BAM_CDEL and seed_qc_sig->cigar_op==BAM_CDEL)){ // indel
			match_flag = isQcSigMatch(qc_sig, seed_qc_sig, QC_SIZE_RATIO_MATCH_THRES_INDEL);
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

