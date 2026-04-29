#include "clipRegCluster.h"

extern pthread_mutex_t mutex_fai;

clipRegCluster::clipRegCluster(string &chrname, int64_t var_startRefPos, int64_t var_endRefPos, int32_t minClipEndSize, int32_t maxVarRegSize, int32_t min_sv_size, int32_t min_supp_num, int32_t minMapQ, double min_seqsim_match, string &technology, faidx_t *fai) {
	this->chrname = chrname;
	this->var_startRefPos = var_startRefPos;
	this->var_endRefPos = var_endRefPos;
	this->minClipEndSize = minClipEndSize;
	this->maxVarRegSize = maxVarRegSize;
	this->min_sv_size = min_sv_size;
	this->min_supp_num = min_supp_num;
	this->minMapQ = minMapQ;
	this->chrlen = faidx_seq_len64(fai, chrname.c_str()); // get reference size
	this->fai = fai;
	this->min_seqsim_match = min_seqsim_match;
	this->technology = technology;

	max_merge_span = CNS_MIN_DIST_MERGE_THRES;
	if(var_endRefPos-var_startRefPos+1>CNS_MIN_DIST_MERGE_THRES) max_merge_span = 2 * CNS_MIN_DIST_MERGE_THRES;
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
	size_t i;

	// extract the features
	qcSigList_vec = extractQcSigsClipReg(query_seq_info_vec);

	for(i=0; i<qcSigList_vec.size(); i++) qcSigList_vec.at(i)->cluster_finished_flag = false;

	// prepare qcSigList information of clipping region for clustering
	prepareQcSigListInfoClipRegForCluster(qcSigList_vec);

#if CLIP_REG_CLUSTER_DEBUG
	cout << "*********** Before cluster, size=" << qcSigList_vec.size() << ":" << endl;
	printQcSigListVec(qcSigList_vec);
#endif

	query_clu_vec = queryClusterOp(qcSigList_vec);

	// append overlap reads and hanging reads
//	appendOverlapHangingReadsClipReg(qcSigList_vec, query_clu_vec);
//
//#if CLIP_REG_CLUSTER_DEBUG
//	for(i=0; i<query_clu_vec.size(); i++){
//		cout << "*********** Queries of cluster " << i << ", size=" << query_clu_vec.at(i)->qcSigList.size() << ":" << endl;
//		printQcSigListVec(query_clu_vec.at(i)->qcSigList);
//	}
//#endif

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

vector<qcSigListVec_t*> clipRegCluster::queryClusterOp(vector<qcSigList_t*> &qcSigList_vec){
	struct MR_pair{
		string qname1, qname2;
		double MR;
	};

	vector<qcSigListVec_t*> query_clu_vec;
	qcSigList_t *qcSigList_node;
	qcSigListVec_t *q_cluster_node_a, *q_cluster_node_b, *tmp;
	size_t i, j, k, smin_id, match_id, num;
	vector<qcSigList_t*> q_cluster_a, q_cluster_b, q_cluster_rescue, q_cluster_a_final, q_cluster_b_final;
	double match_ratio_a, match_ratio_b, match_min, collapse_ratio;
	vector<size_t> min_id_vec, match_id_vec;// span_vec2_id;
	vector<double> match_ratio_vec;
	struct seedQueryInfo *seed_info_a, *seed_info_b;
	bool exist_flag, single_cluster_flag;
	vector<mrmin_match_t*> mrmin_match_vec;
	mrmin_match_t *mrmin_match_node;
	int32_t id_MR1, id_MR_tmp, item_num_collapsed;
	vector<MR_pair> mr_pair_vec;
	MR_pair mr_pair_item;
	int32_t mr_pair_idx;
	string qname1, qname2;

	if(qcSigList_vec.size()>0){
		item_num_collapsed = 0;
		for(i=0; i<qcSigList_vec.size(); i++){
			qcSigList_node = qcSigList_vec.at(i);
			if(qcSigList_node->collapse_to_single_reg_flag==true) item_num_collapsed ++;
		}
		collapse_ratio = (double)item_num_collapsed / qcSigList_vec.size();

		if(collapse_ratio>MIN_RATIO_COLLAPSE_READS){ // sufficient collapsed reads
			//cout << chrname << ":" << var_startRefPos << "-" << var_endRefPos << ", collapse_ratio=" << collapse_ratio << ", single cluster." << endl;

			q_cluster_a_final.clear();
			q_cluster_b_final.clear();
			// copy reads to single cluster
			for(i=0; i<qcSigList_vec.size(); i++){
				qcSigList_node = qcSigList_vec.at(i);
				q_cluster_a_final.push_back(qcSigList_node);
			}
		}else{ // normal clustering
			// initialize a_cluster and b_cluster
			for(k=0; k<qcSigList_vec.size(); k++){

				if(qcSigList_vec.at(k)->qcSig_vec.size()==0) continue;

				qcSigList_vec.at(k)->cluster_finished_flag = true;
				q_cluster_a.push_back(qcSigList_vec.at(k));

				//initialize b_cluster: choose one that is most different query from a_cluster add in b_cluster
				match_min = INT_MAX; //0.99
				match_ratio_vec.clear();
				match_id_vec.clear();
				num = 0;

				match_ratio_vec.push_back(1);
				match_id_vec.push_back(k);
				for(i=1; i<qcSigList_vec.size(); i++){
//					if(qcSigList_vec.at(i)->query_seqs_vec.at(0)->qname.compare("m84039_230415_002321_s3/222761402/ccs")==0){
//						cout << i << ", " << qcSigList_vec.at(i)->query_seqs_vec.at(0)->qname << endl;
//					}
					if(qcSigList_vec.at(i)->qcSig_vec.size()>0 and (qcSigList_vec.at(i)->entire_flanking_flag or qcSigList_vec.at(i)->single_aln_seg_flag)){ // prefer the items flanking the entire region, or the items with single align segment
						seed_info_a = chooseSeedClusterQueryClipReg(qcSigList_vec.at(i), q_cluster_a);
						if(seed_info_a->span_intersection>0){
							qname1 = qcSigList_vec.at(i)->query_seqs_vec.at(0)->qname;
							qname2 = q_cluster_a.at(seed_info_a->id)->query_seqs_vec.at(0)->qname;
							mr_pair_idx = -1;
							for(j=0; j<mr_pair_vec.size(); j++){
								if((mr_pair_vec.at(j).qname1.compare(qname1)==0 and mr_pair_vec.at(j).qname2.compare(qname2)==0) or (mr_pair_vec.at(j).qname1.compare(qname2)==0 and mr_pair_vec.at(j).qname2.compare(qname1)==0)){
									mr_pair_idx = j;
									break;
								}
							}
							if(mr_pair_idx!=-1) // found
								match_ratio_a = mr_pair_vec.at(mr_pair_idx).MR;
							else{
								match_ratio_a = computeMatchRatioClipReg(qcSigList_vec.at(i), q_cluster_a.at(seed_info_a->id), seed_info_a->startSpanPos, seed_info_a->endSpanPos);
								mr_pair_item = {qname1, qname2, match_ratio_a};
								mr_pair_vec.push_back(mr_pair_item);
							}
							match_ratio_vec.push_back(match_ratio_a);
							match_id_vec.push_back(i);
							if(match_ratio_a<match_min){
								match_min = match_ratio_a;
								//smin_id = i;
								//sam_min_id.push_back(i);
							}
							if(match_ratio_a==1) num++;
						}else{
							match_ratio_vec.push_back(2);
						}
						delete seed_info_a;
					}
				}

				if(num<(size_t)min_supp_num){ // insufficient, then try other conditions
					match_min = INT_MAX; //0.99
					match_ratio_vec.clear();
					match_id_vec.clear();

					match_ratio_vec.push_back(1);
					match_id_vec.push_back(k);
					for(i=1; i<qcSigList_vec.size(); i++){
						if(qcSigList_vec.at(i)->qcSig_vec.size()>0){
							seed_info_a = chooseSeedClusterQueryClipReg(qcSigList_vec.at(i), q_cluster_a);
							if(seed_info_a->span_intersection>0){

								qname1 = qcSigList_vec.at(i)->query_seqs_vec.at(0)->qname;
								qname2 = q_cluster_a.at(seed_info_a->id)->query_seqs_vec.at(0)->qname;
								mr_pair_idx = -1;
								for(j=0; j<mr_pair_vec.size(); j++){
									if((mr_pair_vec.at(j).qname1.compare(qname1)==0 and mr_pair_vec.at(j).qname2.compare(qname2)==0) or (mr_pair_vec.at(j).qname1.compare(qname2)==0 and mr_pair_vec.at(j).qname2.compare(qname1)==0)){
										mr_pair_idx = j;
										break;
									}
								}
								if(mr_pair_idx!=-1) // found
									match_ratio_a = mr_pair_vec.at(mr_pair_idx).MR;
								else{
									match_ratio_a = computeMatchRatioClipReg(qcSigList_vec.at(i), q_cluster_a.at(seed_info_a->id), seed_info_a->startSpanPos, seed_info_a->endSpanPos);
									mr_pair_item = {qname1, qname2, match_ratio_a};
									mr_pair_vec.push_back(mr_pair_item);
								}
								match_ratio_vec.push_back(match_ratio_a);
								match_id_vec.push_back(i);
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
				}

				qcSigList_vec.at(k)->cluster_finished_flag = false;
				q_cluster_a.pop_back();

				id_MR1 = -1;
				if(match_min < 1){
					min_id_vec.clear();
					for(i=0; i<match_ratio_vec.size(); i++) if(match_ratio_vec.at(i)<=1) min_id_vec.push_back(i);
					for(i=0; i<min_id_vec.size(); i++){
						exist_flag = false;
						for(j=0; j<mrmin_match_vec.size(); j++){
							if(match_ratio_vec.at(min_id_vec.at(i))==mrmin_match_vec.at(j)->MR){
								smin_id = j;
								exist_flag = true;
								break;
							}
						}
						if(exist_flag==false){
							mrmin_match_node = new mrmin_match_t();
							mrmin_match_node->pos_id.push_back(match_id_vec.at(min_id_vec.at(i)));
							mrmin_match_node->num = 1;
							mrmin_match_node->MR = match_ratio_vec.at(min_id_vec.at(i));
							mrmin_match_node->MR_recluster = -1;
							mrmin_match_node->candidate_flag = false;
							mrmin_match_vec.push_back(mrmin_match_node);
						}else{
							mrmin_match_vec.at(smin_id)->num += 1;
							mrmin_match_vec.at(smin_id)->pos_id.push_back(match_id_vec.at(min_id_vec.at(i)));
						}
					}

					// sort descendingly and initialize clusters adaptively
					sortMRVec(mrmin_match_vec);

#if CLIP_REG_CLUSTER_DEBUG
					printMRVec(mrmin_match_vec);
#endif

					id_MR1 = -1;
					if(mrmin_match_vec.at(0)->MR==1) id_MR1 = 0;
					else if(mrmin_match_vec.size()>=2 and mrmin_match_vec.at(1)->MR==1) id_MR1 = 1;
					if(id_MR1!=-1){ // have MR=1
						for(i=0; i<mrmin_match_vec.at(id_MR1)->pos_id.size(); i++){
							qcSigList_vec.at(mrmin_match_vec.at(id_MR1)->pos_id.at(i))->cluster_finished_flag = true;
							q_cluster_a.push_back(qcSigList_vec.at(mrmin_match_vec.at(id_MR1)->pos_id.at(i)));
						}
						if(mrmin_match_vec.size()>=2){
							if(id_MR1==0) id_MR_tmp = 1;
							else id_MR_tmp = 0;
							//if(mrmin_match_vec.at(id_MR_tmp)->num>=0.6*min_supp_num){ // deleted on 2025-08-15
							if(mrmin_match_vec.at(id_MR_tmp)->num>=READS_NUM_SUPPORT_FACTOR*min_supp_num){
								qcSigList_vec.at(mrmin_match_vec.at(id_MR_tmp)->pos_id.at(0))->cluster_finished_flag = true;
								q_cluster_b.push_back(qcSigList_vec.at(mrmin_match_vec.at(id_MR_tmp)->pos_id.at(0)));
							}
						}
					}else{ // no MR=1
						qcSigList_vec.at(mrmin_match_vec.at(0)->pos_id.at(0))->cluster_finished_flag = true;
						q_cluster_a.push_back(qcSigList_vec.at(mrmin_match_vec.at(0)->pos_id.at(0)));
						//if(mrmin_match_vec.size()>1 and mrmin_match_vec.at(1)->num>=0.6*min_supp_num){ // deleted on 2025-08-15
						if(mrmin_match_vec.size()>1 and mrmin_match_vec.at(1)->num>=READS_NUM_SUPPORT_FACTOR*min_supp_num){
							qcSigList_vec.at(mrmin_match_vec.at(1)->pos_id.at(0))->cluster_finished_flag = true;
							q_cluster_b.push_back(qcSigList_vec.at(mrmin_match_vec.at(1)->pos_id.at(0)));
						}
					}
				}else{ // all MR=1, single cluster
					for(i=0; i<match_id_vec.size(); i++){
						qcSigList_vec.at(match_id_vec.at(i))->cluster_finished_flag = true;
						q_cluster_a.push_back(qcSigList_vec.at(match_id_vec.at(i)));
					}
				}

				single_cluster_flag = false;
				//if(match_min==1 or (mrmin_match_vec.size()>=2 and ((mrmin_match_vec.at(1)->num<0.6*min_supp_num) or (mrmin_match_vec.at(0)->MR==1 and (double)mrmin_match_vec.at(0)->num/(mrmin_match_vec.at(0)->num+mrmin_match_vec.at(1)->num)>=MIN_NUM_RATIO_SINGLE_CLUSTER)))) // deleted on 2025-08-15
				if(match_min==1 or (mrmin_match_vec.size()>=2 and ((mrmin_match_vec.at(1)->num<READS_NUM_SUPPORT_FACTOR*min_supp_num) or (id_MR1!=-1 and mrmin_match_vec.at(id_MR1)->MR==1 and (double)mrmin_match_vec.at(id_MR1)->num/(mrmin_match_vec.at(0)->num+mrmin_match_vec.at(1)->num)>=MIN_NUM_RATIO_SINGLE_CLUSTER))))
					single_cluster_flag = true;

				if(mrmin_match_vec.size()>0){
					vector<mrmin_match_t*>::iterator svp;
					for(svp=mrmin_match_vec.begin(); svp!=mrmin_match_vec.end(); svp++) delete *svp;
					vector<mrmin_match_t*>().swap(mrmin_match_vec);
				}
				// if(q_cluster_b.size()==0){ // deleted on 2024-12-24
				if(q_cluster_b.size()==0 and single_cluster_flag == false){ // try next round
					if(q_cluster_a.size()==1){
						q_cluster_a.pop_back();
						qcSigList_vec.at(k)->cluster_finished_flag = false;
					}else{
						q_cluster_a.clear();
						for(i=0; i<qcSigList_vec.size(); i++) qcSigList_vec.at(i)->cluster_finished_flag = false;
					}
				}else{
					break;
				}
			}

			if(q_cluster_b.size()==0){
				for(i=0; i<match_id_vec.size(); i++){
					match_id = match_id_vec.at(i);
					if(qcSigList_vec.at(match_id)->qcSig_vec.size()>0 and qcSigList_vec.at(match_id)->cluster_finished_flag == false){
						qcSigList_vec.at(match_id)->cluster_finished_flag = true;
						q_cluster_a.push_back(qcSigList_vec.at(match_id));
					}
				}
			}else{
				for(i=0; i<qcSigList_vec.size(); i++){
					if(qcSigList_vec.at(i)->qcSig_vec.size()>0 and qcSigList_vec.at(i)->cluster_finished_flag == false){
						//choose seed query from a_cluster
						seed_info_a = chooseSeedClusterQueryClipReg(qcSigList_vec.at(i), q_cluster_a);
						seed_info_b = chooseSeedClusterQueryClipReg(qcSigList_vec.at(i), q_cluster_b);
						if(seed_info_a->span_intersection>0){
							qname1 = qcSigList_vec.at(i)->query_seqs_vec.at(0)->qname;
							qname2 = q_cluster_a.at(seed_info_a->id)->query_seqs_vec.at(0)->qname;
							mr_pair_idx = -1;
							for(j=0; j<mr_pair_vec.size(); j++){
								if((mr_pair_vec.at(j).qname1.compare(qname1)==0 and mr_pair_vec.at(j).qname2.compare(qname2)==0) or (mr_pair_vec.at(j).qname1.compare(qname2)==0 and mr_pair_vec.at(j).qname2.compare(qname1)==0)){
									mr_pair_idx = j;
									break;
								}
							}
							if(mr_pair_idx!=-1) // found
								match_ratio_a = mr_pair_vec.at(mr_pair_idx).MR;
							else{
								match_ratio_a = computeMatchRatioClipReg(qcSigList_vec.at(i), q_cluster_a.at(seed_info_a->id), seed_info_a->startSpanPos, seed_info_a->endSpanPos);
								mr_pair_item = {qname1, qname2, match_ratio_a};
								mr_pair_vec.push_back(mr_pair_item);
							}
							if(seed_info_b->span_intersection>0){
								qname1 = qcSigList_vec.at(i)->query_seqs_vec.at(0)->qname;
								qname2 = q_cluster_b.at(seed_info_b->id)->query_seqs_vec.at(0)->qname;
								mr_pair_idx = -1;
								for(j=0; j<mr_pair_vec.size(); j++){
									if((mr_pair_vec.at(j).qname1.compare(qname1)==0 and mr_pair_vec.at(j).qname2.compare(qname2)==0) or (mr_pair_vec.at(j).qname1.compare(qname2)==0 and mr_pair_vec.at(j).qname2.compare(qname1)==0)){
										mr_pair_idx = j;
										break;
									}
								}
								if(mr_pair_idx!=-1) // found
									match_ratio_b = mr_pair_vec.at(mr_pair_idx).MR;
								else{
									match_ratio_b = computeMatchRatioClipReg(qcSigList_vec.at(i), q_cluster_b.at(seed_info_b->id), seed_info_b->startSpanPos, seed_info_b->endSpanPos);
									mr_pair_item = {qname1, qname2, match_ratio_b};
									mr_pair_vec.push_back(mr_pair_item);
								}
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
								qname1 = qcSigList_vec.at(i)->query_seqs_vec.at(0)->qname;
								qname2 = q_cluster_b.at(seed_info_b->id)->query_seqs_vec.at(0)->qname;
								mr_pair_idx = -1;
								for(j=0; j<mr_pair_vec.size(); j++){
									if((mr_pair_vec.at(j).qname1.compare(qname1)==0 and mr_pair_vec.at(j).qname2.compare(qname2)==0) or (mr_pair_vec.at(j).qname1.compare(qname2)==0 and mr_pair_vec.at(j).qname2.compare(qname1)==0)){
										mr_pair_idx = j;
										break;
									}
								}
								if(mr_pair_idx!=-1) // found
									match_ratio_b = mr_pair_vec.at(mr_pair_idx).MR;
								else{
									match_ratio_b = computeMatchRatioClipReg(qcSigList_vec.at(i), q_cluster_b.at(seed_info_b->id), seed_info_b->startSpanPos, seed_info_b->endSpanPos);
									mr_pair_item = {qname1, qname2, match_ratio_b};
									mr_pair_vec.push_back(mr_pair_item);
								}
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

					qname1 = q_cluster_rescue.at(i)->query_seqs_vec.at(0)->qname;
					qname2 = q_cluster_a.at(seed_info_a->id)->query_seqs_vec.at(0)->qname;
					mr_pair_idx = -1;
					for(j=0; j<mr_pair_vec.size(); j++){
						if((mr_pair_vec.at(j).qname1.compare(qname1)==0 and mr_pair_vec.at(j).qname2.compare(qname2)==0) or (mr_pair_vec.at(j).qname1.compare(qname2)==0 and mr_pair_vec.at(j).qname2.compare(qname1)==0)){
							mr_pair_idx = j;
							break;
						}
					}
					if(mr_pair_idx!=-1) // found
						match_ratio_a = mr_pair_vec.at(mr_pair_idx).MR;
					else{
						match_ratio_a = computeMatchRatioClipReg(q_cluster_rescue.at(i), q_cluster_a.at(seed_info_a->id), seed_info_a->startSpanPos, seed_info_a->endSpanPos);
						mr_pair_item = {qname1, qname2, match_ratio_a};
						mr_pair_vec.push_back(mr_pair_item);
					}

					qname1 = q_cluster_rescue.at(i)->query_seqs_vec.at(0)->qname;
					qname2 = q_cluster_b.at(seed_info_b->id)->query_seqs_vec.at(0)->qname;
					mr_pair_idx = -1;
					for(j=0; j<mr_pair_vec.size(); j++){
						if((mr_pair_vec.at(j).qname1.compare(qname1)==0 and mr_pair_vec.at(j).qname2.compare(qname2)==0) or (mr_pair_vec.at(j).qname1.compare(qname2)==0 and mr_pair_vec.at(j).qname2.compare(qname1)==0)){
							mr_pair_idx = j;
							break;
						}
					}
					if(mr_pair_idx!=-1) // found
						match_ratio_b = mr_pair_vec.at(mr_pair_idx).MR;
					else{
						match_ratio_b = computeMatchRatioClipReg(q_cluster_rescue.at(i), q_cluster_b.at(seed_info_b->id), seed_info_b->startSpanPos, seed_info_b->endSpanPos);
						mr_pair_item = {qname1, qname2, match_ratio_b};
						mr_pair_vec.push_back(mr_pair_item);
					}
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
		}

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
	}

#if CLIP_REG_CLUSTER_DEBUG
	for(i=0; i<query_clu_vec.size(); i++){
		cout << "*********** Queries of cluster " << i << ", size=" << query_clu_vec.at(i)->qcSigList.size() << ":" << endl;
		printQcSigListVec(query_clu_vec.at(i)->qcSigList);
	}
#endif

	return query_clu_vec;
}

//  append overlap reads and hanging reads
void clipRegCluster::appendOverlapHangingReadsClipReg(vector<qcSigList_t*> &qcSigList_vec, vector<qcSigListVec_t*> &query_clu_vec){
	size_t i;
	clipAlnData_t *clip_aln;
	qcSigList_t *qcSigList_node;
	int64_t target_idx, var_startRefPos_tmp, var_endRefPos_tmp, start_checkPos_left, end_checkPos_left, start_checkPos_right, end_checkPos_right;
	bool valid_flag;

	var_startRefPos_tmp = var_startRefPos - GROUP_DIST_THRES;
	var_endRefPos_tmp = var_endRefPos + GROUP_DIST_THRES;
	if(var_startRefPos_tmp<1) var_startRefPos_tmp = 1;
	if(var_endRefPos_tmp>chrlen) var_endRefPos_tmp = chrlen;

	start_checkPos_left = var_startRefPos_tmp;
	end_checkPos_left = var_startRefPos + GROUP_DIST_THRES;
	if(end_checkPos_left>chrlen) end_checkPos_left = chrlen;

	start_checkPos_right = var_endRefPos - GROUP_DIST_THRES;
	if(start_checkPos_right<1) start_checkPos_right = 1;
	end_checkPos_right = var_endRefPos_tmp;

	target_idx = -1;
	if(query_clu_vec.size()>=1) target_idx = 0;

	if(target_idx!=-1){
		for(i=0; i<qcSigList_vec.size(); i++){
			qcSigList_node = qcSigList_vec.at(i);
			if(qcSigList_node->cluster_finished_flag==false and qcSigList_node->query_seqs_vec.size()==1){
				clip_aln = qcSigList_node->query_seqs_vec.at(0)->clip_aln;
				valid_flag = false;
				if(clip_aln->leftClipSize>=minClipEndSize and clip_aln->startRefPos>=start_checkPos_left and clip_aln->startRefPos<=end_checkPos_left)  // left side
					valid_flag = true;
				else if(clip_aln->rightClipSize>=minClipEndSize and clip_aln->endRefPos>=start_checkPos_right and clip_aln->endRefPos<=end_checkPos_right) // right side
					valid_flag = true;

				//if((qcSigList_node->startSpanRefPos<var_startRefPos and qcSigList_node->endSpanRefPos>var_startRefPos_tmp and qcSigList_node->endSpanRefPos<var_endRefPos_tmp) or (qcSigList_node->endSpanRefPos>var_endRefPos and qcSigList_node->startSpanRefPos>var_startRefPos_tmp and qcSigList_node->startSpanRefPos<var_endRefPos_tmp)){
				if(valid_flag){
					cout << i << ", qname=" << qcSigList_node->query_seqs_vec.at(0)->qname << endl;
					qcSigList_node->cluster_finished_flag = true;
					query_clu_vec.at(target_idx)->qcSigList.push_back(qcSigList_node);
				}
			}
		}
	}
}

// prepare qcSigList information of clipping region for clustering
void clipRegCluster::prepareQcSigListInfoClipRegForCluster(vector<qcSigList_t*> &qcSigList_vec){
	size_t i, j;
	qcSigList_t *qcSigList_node, *qcSigList_node2, *tmp_qcSigList_node;
	struct querySeqInfoNode *q_info;
	qcSig_t *qc_sig;
	bool exist_flag;
	double ratio;
	int64_t var_startRefPos_tmp, var_endRefPos_tmp;

	var_startRefPos_tmp = var_startRefPos - VAR_ALN_EXTEND_SIZE;
	if(var_startRefPos_tmp<1) var_startRefPos_tmp = 1;
	var_endRefPos_tmp = var_endRefPos + VAR_ALN_EXTEND_SIZE;
	if(var_endRefPos_tmp>chrlen) var_endRefPos_tmp = chrlen;

	for(i=0; i<qcSigList_vec.size(); i++){
		qcSigList_node = qcSigList_vec.at(i);
		qcSigList_node->collapse_to_single_reg_flag = false;
		qcSigList_node->total_valid_seg_num = qcSigList_node->seg_num_aln_to_single_reg = 0;
		for(j=0; j<qcSigList_node->query_seqs_vec.size(); j++){
			q_info = qcSigList_node->query_seqs_vec.at(j);
			if(q_info->clip_aln and q_info->clip_aln->query_dist>=VAR_ALN_EXTEND_SIZE and q_info->clip_aln->ref_dist>=VAR_ALN_EXTEND_SIZE){
				qcSigList_node->total_valid_seg_num ++;
				if(q_info->clip_aln->startRefPos>var_startRefPos_tmp and q_info->clip_aln->endRefPos<var_endRefPos_tmp)
					qcSigList_node->seg_num_aln_to_single_reg ++;
			}
		}
		ratio = (double)qcSigList_node->seg_num_aln_to_single_reg / qcSigList_node->total_valid_seg_num;
		if(qcSigList_node->query_seqs_vec.size()>=MIN_SEG_NUM_READ_COLLAPSE and ratio>=MIN_SEG_RATIO_READ_COLLAPSE) qcSigList_node->collapse_to_single_reg_flag = true;

		//cout << "qname=" << qcSigList_node->query_seqs_vec.at(0)->qname << ", q_vec.size=" << qcSigList_node->query_seqs_vec.size() << ", seg_num_aln_to_single_reg=" << qcSigList_node->seg_num_aln_to_single_reg << ", ratio=" << ratio << ", collapse_to_single_reg_flag=" << qcSigList_node->collapse_to_single_reg_flag << endl;
	}
	//cout << "item_num=" << item_num << ", ratio=" << (double)item_num/qcSigList_vec.size() << endl;

	// sort in descending order according to the number of features
	for(i=0; i<qcSigList_vec.size(); i++){
		qcSigList_node = qcSigList_vec.at(i);
		for(j=i+1; j<qcSigList_vec.size(); j++){
			qcSigList_node2 = qcSigList_vec.at(j);
			if(qcSigList_node->qcSig_vec.size()<qcSigList_node2->qcSig_vec.size()){
				tmp_qcSigList_node = qcSigList_vec.at(i);
				qcSigList_vec.at(i) = qcSigList_vec.at(j);
				qcSigList_vec.at(j) = tmp_qcSigList_node;
			}
		}
	}

	// fill ins_sum, del_sum and inv_sum
	for(i=0; i<qcSigList_vec.size(); i++){
		qcSigList_node = qcSigList_vec.at(i);
		qcSigList_node->ins_sum = qcSigList_node->del_sum = qcSigList_node->inv_sum = 0;
		for(j=0; j<qcSigList_node->qcSig_vec.size(); j++){
			qc_sig = qcSigList_node->qcSig_vec.at(j);
			if(qc_sig->cigar_op==BAM_CINS) qcSigList_node->ins_sum += qc_sig->cigar_op_len;
			else if(qc_sig->cigar_op==BAM_CDEL) qcSigList_node->del_sum += qc_sig->cigar_op_len;
			else if(qc_sig->cigar_op==CIGAR_OP_INV_CLIP) qcSigList_node->inv_sum += qc_sig->cigar_op_len;
		}
	}

	// sort according to the category of the number of signatures
	sortQueryInfoByNumCategoryClipReg(qcSigList_vec);

	// prefer the alignments flanking the entire variant region
	for(i=0; i<qcSigList_vec.size(); i++){
		qcSigList_node = qcSigList_vec.at(i);
		if(qcSigList_node->qcSig_vec.size()>0 and qcSigList_node->entire_flanking_flag==false){
			exist_flag = false;
			for(j=i+1; j<qcSigList_vec.size(); j++){
				qcSigList_node2 = qcSigList_vec.at(j);
				if(qcSigList_node2->qcSig_vec.size()>0 and qcSigList_node2->entire_flanking_flag){
					exist_flag = true;
					tmp_qcSigList_node = qcSigList_vec.at(i);
					qcSigList_vec.at(i) = qcSigList_vec.at(j);
					qcSigList_vec.at(j) = tmp_qcSigList_node;
					break;
				}
			}
			if(exist_flag==false) break;
		}
	}

	// sort according to the category of the size of signatures
	sortQueryInfoBySizeCategoryClipReg(qcSigList_vec);
}

void clipRegCluster::sortQueryInfoByNumCategoryClipReg(vector<qcSigList_t*> &qcSigList_vec){
	size_t i, j, sum;
	vector<qcSigList_t*> qcSigList_vec_tmp;
	int32_t idx;
	vector<mrmin_match_t*> num_id_vec;
	mrmin_match_t *num_id_node;
	bool flag;

	for(i=0; i<qcSigList_vec.size(); i++){
		flag = false;
		for(j=0; j<num_id_vec.size(); j++){
			if(num_id_vec.at(j)->num==(int32_t)qcSigList_vec.at(i)->qcSig_vec.size()){
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
			num_id_node->num = qcSigList_vec.at(i)->qcSig_vec.size();
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
	if(sum==qcSigList_vec.size()){
		for(i=0; i<num_id_vec.size(); i++){
			num_id_node = num_id_vec.at(i);
			for(j=0; j<(size_t)num_id_node->pos_id.size(); j++) qcSigList_vec_tmp.push_back(qcSigList_vec.at(num_id_node->pos_id.at(j)));
		}
		qcSigList_vec.assign(qcSigList_vec_tmp.begin(), qcSigList_vec_tmp.end());
	}else{
		cerr << "sum=" << sum << ", qcSigList_vec.size=" << qcSigList_vec.size() << ", error!" << endl;
		exit(1);
	}

	// release memory
	for(i=0; i<num_id_vec.size(); i++) delete num_id_vec.at(i);
	vector<mrmin_match_t*>().swap(num_id_vec);
}

void clipRegCluster::sortQueryInfoBySizeCategoryClipReg(vector<qcSigList_t*> &qcSigList_vec){
	size_t i;
	qcSigList_t *qcSigList_node;
	vector<qcSigList_t*> qcSigList_vec_tmp;
	double size_ratio;
	vector<pair<int32_t, vector<int32_t>>> size_vec;
	pair<int32_t, vector<int32_t>> pair_size_node;
	vector<int32_t> size_vec_tmp;
	bool flag;

	// compute mean size
	for(i=0; i<qcSigList_vec.size(); i++){
		qcSigList_node = qcSigList_vec.at(i);
		if(qcSigList_node->qcSig_vec.size()==0 or (qcSigList_node->ins_sum==0 and qcSigList_node->del_sum==0 and qcSigList_node->inv_sum==0)) continue; // skip items of no signatures

		if(qcSigList_node->ins_sum>0 and qcSigList_node->ins_sum>=qcSigList_node->del_sum and qcSigList_node->ins_sum>qcSigList_node->inv_sum){ // max: INS
			flag = false;
			for(auto& pair1 : size_vec){
				const int32_t& size_key = pair1.first;
				vector<int32_t>& id_vec = pair1.second;
				if(qcSigList_node->ins_sum>size_key) size_ratio = (double)size_key / qcSigList_node->ins_sum;
				else size_ratio = (double)qcSigList_node->ins_sum / size_key;
				if(size_ratio>=QC_SIZE_RATIO_MATCH_THRES_INDEL){
					id_vec.push_back(i);
					flag = true;
					break;
				}
			}
			if(flag==false){
				size_vec_tmp.clear();
				size_vec_tmp.push_back(i);
				pair_size_node.first = qcSigList_node->ins_sum;
				pair_size_node.second = size_vec_tmp;
				size_vec.push_back(pair_size_node);
			}
		}else if(qcSigList_node->del_sum>0 and qcSigList_node->del_sum>=qcSigList_node->ins_sum and qcSigList_node->del_sum>qcSigList_node->inv_sum){ // max: DEL
			flag = false;
			for(auto& pair1 : size_vec){
				const int32_t& size_key = pair1.first;
				vector<int32_t>& id_vec = pair1.second;
				if(qcSigList_node->del_sum>abs(size_key)) size_ratio = (double)abs(size_key) / qcSigList_node->del_sum;
				else size_ratio = (double)qcSigList_node->del_sum / abs(size_key);
				if(size_ratio>=QC_SIZE_RATIO_MATCH_THRES_INDEL){
					id_vec.push_back(i);
					flag = true;
					break;
				}
			}
			if(flag==false){
				size_vec_tmp.clear();
				size_vec_tmp.push_back(i);
				pair_size_node.first = -qcSigList_node->del_sum;
				pair_size_node.second = size_vec_tmp;
				size_vec.push_back(pair_size_node);
			}
		}else if(qcSigList_node->inv_sum>0 and qcSigList_node->inv_sum>=qcSigList_node->ins_sum and qcSigList_node->inv_sum>=qcSigList_node->del_sum){ // max: INV
			flag = false;
			for(auto& pair1 : size_vec){
				const int32_t& size_key = pair1.first;
				vector<int32_t>& id_vec = pair1.second;
				if(qcSigList_node->inv_sum>abs(size_key)) size_ratio = (double)abs(size_key) / qcSigList_node->inv_sum;
				else size_ratio = (double)qcSigList_node->inv_sum / abs(size_key);
				if(size_ratio>=QC_SIZE_RATIO_MATCH_THRES_INDEL){
					id_vec.push_back(i);
					flag = true;
					break;
				}
			}
			if(flag==false){
				size_vec_tmp.clear();
				size_vec_tmp.push_back(i);
				pair_size_node.first = qcSigList_node->inv_sum;
				pair_size_node.second = size_vec_tmp;
				size_vec.push_back(pair_size_node);
			}
		}
	}

	// sort according to size count category
	sort(size_vec.begin(), size_vec.end(), desc_cmp_size_vec);

	// copy items
	for(auto& pair1 : size_vec){
		vector<int32_t>& id_vec = pair1.second;
		for(auto& id: id_vec){
			qcSigList_node = qcSigList_vec.at(id);
			qcSigList_vec_tmp.push_back(qcSigList_node);
		}
	}
	for(i=0; i<qcSigList_vec.size(); i++){
		qcSigList_node = qcSigList_vec.at(i);
		if(qcSigList_node->qcSig_vec.size()==0 or (qcSigList_node->ins_sum==0 and qcSigList_node->del_sum==0 and qcSigList_node->inv_sum==0)){
			qcSigList_vec_tmp.push_back(qcSigList_node);
		}
	}
	if(qcSigList_vec_tmp.size()!=qcSigList_vec.size()){
		cerr << "line=" << __LINE__ << ", qcSigList_vec_tmp.size=" << qcSigList_vec_tmp.size() << ", qcSigList_vec.size=" << qcSigList_vec.size() << ", error!" << endl;
		exit(1);
	}

	qcSigList_vec.assign(qcSigList_vec_tmp.begin(), qcSigList_vec_tmp.end());
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
	int64_t item_num_collapse;

	if(qcSigList_vec.size()>0){
		// print the feature
		item_num_collapse = 0;
		for(i=0; i<qcSigList_vec.size(); i++){
			qcSigList_node = qcSigList_vec.at(i);
			cout << "[" << i << "]: " << qcSigList_node->query_seqs_vec.at(0)->qname << ", cluster_finished=" << qcSigList_node->cluster_finished_flag << ", entire_flanking=" << qcSigList_node->entire_flanking_flag << ", single_aln_seg=" << qcSigList_node->single_aln_seg_flag << ", ins_sum=" << qcSigList_node->ins_sum << ", del_sum=" << qcSigList_node->del_sum << ", qcSig_vec.size=" << qcSigList_node->qcSig_vec.size() << ", total_valid_seg_num=" << qcSigList_node->total_valid_seg_num << ", seg_num_aln_to_single_reg=" << qcSigList_node->seg_num_aln_to_single_reg << ", collapse_to_single_reg_flag=" << qcSigList_node->collapse_to_single_reg_flag << ":";
			for(j=0; j<qcSigList_node->qcSig_vec.size(); j++){
				qc_sig = qcSigList_node->qcSig_vec.at(j);
				if(qc_sig->cigar_op==BAM_CSOFT_CLIP or qc_sig->cigar_op==BAM_CHARD_CLIP) // clippings
					cout << "\t(" << qc_sig->cigar_op << "," << qc_sig->chrname << ":" << qc_sig->ref_pos << "," << qc_sig->dist_next_clip << "," << qc_sig->inv_flag << ")";
				else // indels
					cout << "\t(" << qc_sig->cigar_op << "," << qc_sig->chrname << ":" << qc_sig->ref_pos << "," << qc_sig->cigar_op_len << "," << qc_sig->inv_flag << ")";
			}
			cout << endl;

			if(qcSigList_node->collapse_to_single_reg_flag) item_num_collapse ++;
		}
		cout << "item_num_collapse=" << item_num_collapse << ", ratio=" << (double)item_num_collapse/qcSigList_vec.size() << endl;
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

//			if(queryseq_node->qname.compare("m84039_230414_235240_s2/162599892/ccs")==0){
//				cout << "i=" << i << ", " << queryseq_node->qname << endl;
//			}else continue;

			query_seqs = getQuerySeqs(queryseq_node->qname, query_seq_info_vec);
			qcSigList_node = extractQcSigsSingleQueryClipReg(query_seqs);
			if(qcSigList_node){
				if(qcSigList_node->qcSig_vec.size()>=2){
					rmNeighborFalseSig(qcSigList_node->qcSig_vec, MIN_AVER_SIZE_ALN_SEG, MIN_SIZE_RATIO_MATCH_CLIP_POS);
					//qcSigList_node->qcSig_merge_flag = mergeNeighbouringSigsFlagClipReg(qcSigList_node, qcSigList_node->qcSig_vec, MAX_DIST_MERGE_ARBITARY, var_endRefPos-var_startRefPos, CNS_MIN_SEQSIM_MERGE_THRES2, MIN_VALID_SIG_SIZE_RATIO_THRES, fai);
					qcSigList_node->qcSig_merge_flag = mergeNeighbouringSigsFlagClipReg(qcSigList_node, qcSigList_node->qcSig_vec, MAX_DIST_MERGE_ARBITARY, max_merge_span, CNS_MIN_SEQSIM_MERGE_THRES2, MIN_VALID_SIG_SIZE_RATIO_THRES, fai);
				}
				filterSmallSigs(qcSigList_node->qcSig_vec, MIN_SIG_SIZE_RATIO_FILTER_THRES);

				qcSigList_vec.push_back(qcSigList_node);
			}

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
	int32_t id, leftmost_id, rightmost_id, increase_direction, id3;
	int32_t left_clip_end_flag, right_clip_end_flag, end_flag;
	bool overlap_flag, same_orient_flag, valid_query_left_end_flag, valid_query_right_end_flag, self_overlap_flag, ref_skip_flag, query_skip_flag, valid_flag, inv_flag;
	bool overlap_flag_left_end, overlap_flag_right_end;
	size_t i, j;
	int32_t mate_arr_idx, mate_clip_end_flag, mate_arr_idx_same_chr, mate_clip_end_flag_same_chr, min_dist, min_ref_dist, min_dist_same_chr, min_ref_dist_same_chr, seq_len, ignore_end_flag, clip_end_flag3;
	int64_t var_startRefPos_tmp, var_endRefPos_tmp, end_ref_pos, dist;
	int64_t start_ref_pos_var, end_ref_pos_var, start_query_pos_var, end_query_pos_var, tmp_pos, min_idx, min_val;
	string reg_str, refseq, queryseq_whole, reverse_seq;
	char *p_seq;
	int32_t seq_id_arr[query_seqs.size()];
	clipAlnData_t *clip_aln3;

	var_startRefPos_tmp = var_startRefPos - GROUP_DIST_THRES * 2;
	var_endRefPos_tmp = var_endRefPos + GROUP_DIST_THRES * 2;
	if(var_startRefPos_tmp<1) var_startRefPos_tmp = 1;
	if(var_endRefPos_tmp>chrlen) var_endRefPos_tmp = chrlen;

	queryseq_whole = "";
	//query_orient = -1;
	for(i=0; i<query_seqs.size(); i++){
		queryseq_node = query_seqs.at(i);
		if(queryseq_node->clip_aln and queryseq_node->clip_aln->leftHardClippedFlag==false and queryseq_node->clip_aln->rightHardClippedFlag==false){
			queryseq_whole = queryseq_node->seq;
			//query_orient = queryseq_node->clip_aln->aln_orient;
			break;
		}
	}
	//cout << "queryseq_whole=" << queryseq_whole << endl;
	if(queryseq_whole.size()==0) return NULL; // ignore all hard clip segments

	for(i=0; i<query_seqs.size(); i++) seq_id_arr[i] = 0;

	sidemost_idx_vec = getLeftRightMostAlnSegIdx(query_seqs);
	//cout << "qname=" << query_seqs.at(0)->qname << ", query_orient=" << to_string(query_orient) << ", minQposIdx=" << sidemost_idx_vec.at(0) << ", maxQposIdx=" << sidemost_idx_vec.at(1) << ", seq_num=" << query_seqs.size() << endl;

	valid_flag = true;
//	if(query_orient==ALN_PLUS_ORIENT) { // plus orient
		leftmost_id = sidemost_idx_vec.at(0);
		rightmost_id = sidemost_idx_vec.at(1);
//	}else{ // minus orient
//		leftmost_id = sidemost_idx_vec.at(1);
//		rightmost_id = sidemost_idx_vec.at(0);
//	}
	if(leftmost_id==-1 or rightmost_id==-1) {
		cerr << "line=" << __LINE__ << ", leftmost_id=" << leftmost_id << ", rightmost_id=" << rightmost_id << ", error!" << endl;
		leftmost_id = 0;
		exit(1);
	}
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
				overlap_flag_left_end = overlap_flag = isOverlappedPos(queryseq_node->query_alnSegs.at(0)->startRpos, end_ref_pos, var_startRefPos_tmp, var_endRefPos_tmp);
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
					qc_sig->reg_contain_flag = qc_sig->have_next_clip = qc_sig->merge_flag = qc_sig->inv_flag = false;
					qc_sig->altseq_complete_flag = true;
					//qc_sig->ref_pos = queryseq_node->query_alnSegs.at(0)->startRpos;
					qc_sig->ref_pos = queryseq_node->clip_aln->startRefPos - 1;
					if(queryseq_node->clip_aln->aln_orient==ALN_PLUS_ORIENT) qc_sig->query_pos = queryseq_node->clip_aln->startQueryPos;
					else qc_sig->query_pos = queryseq_node->clip_aln->querylen - queryseq_node->clip_aln->startQueryPos + 1;
//					end_query_pos = getEndQueryPosAlnSeg(queryseq_node->query_alnSegs.at(0)->startQpos, queryseq_node->query_alnSegs.at(0)->opflag, queryseq_node->query_alnSegs.at(0)->seglen);
//					if(queryseq_node->query_alnSegs.at(0)->opflag==BAM_CSOFT_CLIP) qc_sig->query_pos = end_query_pos;
//					else if(queryseq_node->query_alnSegs.at(0)->opflag==BAM_CHARD_CLIP) qc_sig->query_pos = queryseq_node->clip_aln->leftClipSize + end_query_pos;
//					else qc_sig->query_pos = queryseq_node->query_alnSegs.at(0)->startQpos;
					qc_sig->chrname_next_clip = "";
					qc_sig->dist_next_clip = 0;
					qc_sig->have_next_clip = false;
					qc_sig->mate_qcSig = NULL;
					qcSig_vec.push_back(qc_sig);
					last_clip_sig = qc_sig;

					if(qc_sig->cigar_op==BAM_CSOFT_CLIP or qc_sig->cigar_op==BAM_CHARD_CLIP){
						adjClipAlnSegInfo = getAdjacentAlnSegClipReg(id, left_clip_end_flag, query_seqs, minClipEndSize);
						mate_arr_idx = adjClipAlnSegInfo.at(0);
						mate_clip_end_flag = adjClipAlnSegInfo.at(1);
						mate_arr_idx_same_chr = adjClipAlnSegInfo.at(4);
						mate_clip_end_flag_same_chr = adjClipAlnSegInfo.at(5);
						if(mate_arr_idx_same_chr!=-1){
							mate_queryseq_node = query_seqs.at(mate_arr_idx_same_chr);
							// get another segment
							id3 = clip_end_flag3 = -1;
							clip_aln3 = NULL;
							if(mate_clip_end_flag_same_chr==RIGHT_END){
								clip_end_flag3 = LEFT_END;
								clip_aln3 = mate_queryseq_node->clip_aln->left_aln;
							}else{
								clip_end_flag3 = RIGHT_END;
								clip_aln3 = mate_queryseq_node->clip_aln->right_aln;
							}
							if(clip_aln3){
								for(i=0; i<query_seqs.size(); i++){
									if(clip_aln3==query_seqs.at(i)->clip_aln and seq_id_arr[i]==0){
										id3 = i;
										break;
									}
								}
							}

							if(queryseq_node->clip_aln->aln_orient!=mate_queryseq_node->clip_aln->aln_orient){ // different orientation
								if(clip_aln3==NULL){
									end_ref_pos = getEndRefPosAlnSeg(queryseq_node->query_alnSegs.at(queryseq_node->query_alnSegs.size()-1)->startRpos, queryseq_node->query_alnSegs.at(queryseq_node->query_alnSegs.size()-1)->opflag, queryseq_node->query_alnSegs.at(queryseq_node->query_alnSegs.size()-1)->seglen);
									overlap_flag_right_end = isOverlappedPos(queryseq_node->query_alnSegs.at(queryseq_node->query_alnSegs.size()-1)->startRpos, end_ref_pos, var_startRefPos_tmp, var_endRefPos_tmp);
									if(overlap_flag_right_end){ // both ends overlap, then queryseq_node is the main inverted node
										//start_ref_pos_var = queryseq_node->query_alnSegs.at(0)->startRpos;
										start_ref_pos_var = queryseq_node->clip_aln->startRefPos;
										if(queryseq_node->clip_aln->aln_orient==ALN_PLUS_ORIENT) qc_sig->query_pos = queryseq_node->clip_aln->startQueryPos;
										else qc_sig->query_pos = queryseq_node->clip_aln->querylen - queryseq_node->clip_aln->startQueryPos + 1;

//										end_query_pos = getEndQueryPosAlnSeg(queryseq_node->query_alnSegs.at(0)->startQpos, queryseq_node->query_alnSegs.at(0)->opflag, queryseq_node->query_alnSegs.at(0)->seglen);
//										if(queryseq_node->query_alnSegs.at(0)->opflag==BAM_CSOFT_CLIP) start_query_pos_var = end_query_pos;
//										else if(queryseq_node->query_alnSegs.at(0)->opflag==BAM_CHARD_CLIP) start_query_pos_var = queryseq_node->clip_aln->leftClipSize + end_query_pos;
//										else start_query_pos_var = queryseq_node->query_alnSegs.at(0)->startQpos;

										if(mate_clip_end_flag_same_chr==LEFT_END){
											end_ref_pos_var = mate_queryseq_node->clip_aln->startRefPos;
											//end_query_pos_var = mate_queryseq_node->clip_aln->startQueryPos;
											if(mate_queryseq_node->clip_aln->aln_orient==ALN_PLUS_ORIENT) end_query_pos_var = mate_queryseq_node->clip_aln->startQueryPos;
											else end_query_pos_var = mate_queryseq_node->clip_aln->querylen - mate_queryseq_node->clip_aln->startQueryPos + 1;
										}else{
											end_ref_pos_var = mate_queryseq_node->clip_aln->endRefPos;
											//end_query_pos_var = mate_queryseq_node->clip_aln->endQueryPos;
											if(mate_queryseq_node->clip_aln->aln_orient==ALN_PLUS_ORIENT) end_query_pos_var = mate_queryseq_node->clip_aln->endQueryPos;
											else end_query_pos_var = mate_queryseq_node->clip_aln->querylen - mate_queryseq_node->clip_aln->endQueryPos + 1;
										}
//										start_ref_pos_var = queryseq_node->clip_aln->startRefPos;
//										start_query_pos_var = queryseq_node->clip_aln->startQueryPos;
//										if(mate_clip_end_flag_same_chr==LEFT_END){
//											end_ref_pos_var = mate_queryseq_node->clip_aln->startRefPos;
//											end_query_pos_var = mate_queryseq_node->clip_aln->startQueryPos;
//										}else{
//											end_ref_pos_var = mate_queryseq_node->clip_aln->endRefPos;
//											end_query_pos_var = mate_queryseq_node->clip_aln->endQueryPos;
//										}
										qc_sig->chrname = queryseq_node->clip_aln->chrname;
										qc_sig->aln_orient = queryseq_node->clip_aln->aln_orient;
									}else{ // else, mate_queryseq_node is the main inverted node
										if(mate_clip_end_flag_same_chr==LEFT_END){
											start_ref_pos_var = mate_queryseq_node->clip_aln->startRefPos;
											//start_query_pos_var = mate_queryseq_node->clip_aln->startQueryPos;
											if(mate_queryseq_node->clip_aln->aln_orient==ALN_PLUS_ORIENT) start_query_pos_var = mate_queryseq_node->clip_aln->startQueryPos;
											else start_query_pos_var = mate_queryseq_node->clip_aln->querylen - mate_queryseq_node->clip_aln->startQueryPos + 1;
										}else{
											start_ref_pos_var = mate_queryseq_node->clip_aln->endRefPos;
											//start_query_pos_var = mate_queryseq_node->clip_aln->endQueryPos;
											if(mate_queryseq_node->clip_aln->aln_orient==ALN_PLUS_ORIENT) start_query_pos_var = mate_queryseq_node->clip_aln->endQueryPos;
											else start_query_pos_var = mate_queryseq_node->clip_aln->querylen - mate_queryseq_node->clip_aln->endQueryPos + 1;
										}
										end_ref_pos_var = queryseq_node->clip_aln->startRefPos;
										//end_query_pos_var = queryseq_node->clip_aln->startQueryPos;
										if(queryseq_node->clip_aln->aln_orient==ALN_PLUS_ORIENT) end_query_pos_var = queryseq_node->clip_aln->startQueryPos;
										else end_query_pos_var = queryseq_node->clip_aln->querylen - queryseq_node->clip_aln->startQueryPos + 1;

										qc_sig->chrname = mate_queryseq_node->clip_aln->chrname;
										qc_sig->aln_orient = mate_queryseq_node->clip_aln->aln_orient;
									}

									if(start_ref_pos_var>end_ref_pos_var){
										tmp_pos = start_ref_pos_var;
										start_ref_pos_var = end_ref_pos_var;
										end_ref_pos_var = tmp_pos;
										tmp_pos = start_query_pos_var;
										start_query_pos_var = end_query_pos_var;
										end_query_pos_var = tmp_pos;
									}

									qc_sig->cigar_op = CIGAR_OP_INV_CLIP;
									if(start_ref_pos_var<end_ref_pos_var) qc_sig->cigar_op_len = end_ref_pos_var - start_ref_pos_var + 1;
									else qc_sig->cigar_op_len = start_ref_pos_var - end_ref_pos_var + 1;
									qc_sig->inv_flag = true;
									//last_clip_sig->reg_contain_flag = qc_sig->have_next_clip = qc_sig->merge_flag = false;
									qc_sig->ref_pos = start_ref_pos_var;
									qc_sig->end_ref_pos = end_ref_pos_var;
									qc_sig->query_pos = start_query_pos_var;
									qc_sig->end_query_pos = end_query_pos_var;

									seq_id_arr[mate_arr_idx_same_chr] = 1;
//									mate_arr_idx = id3;
//									mate_clip_end_flag = clip_end_flag3;
								}
							}

						}
					}
				}
			}

			// indel
			qcSig_vec_indel = extractQcSigsFromAlnSegsSingleQuery(queryseq_node, queryseq_node->clip_aln->chrname, var_startRefPos_tmp, var_endRefPos_tmp, MIN_SV_SIZE_CLIP_REG); // min_sv_size
			rmNeighborFalseSig(qcSig_vec_indel, MIN_AVER_SIZE_ALN_SEG, MIN_SIZE_RATIO_MATCH_CLIP_POS);
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
					if(queryseq_node->clip_aln->leftHardClippedFlag==false) qc_sig->altseq = queryseq_node->seq.substr(qc_sig->query_pos, 1);
					else qc_sig->altseq = queryseq_node->seq.substr(qc_sig->query_pos-queryseq_node->clip_aln->leftClipSize, 1);
				}
				qcSig_vec.push_back(qc_sig);
			}

			// clipping at right end
			if(ignore_end_flag!=RIGHT_END){
				valid_query_right_end_flag = false;
				right_clip_end_flag = -1;
				end_ref_pos = getEndRefPosAlnSeg(queryseq_node->query_alnSegs.at(queryseq_node->query_alnSegs.size()-1)->startRpos, queryseq_node->query_alnSegs.at(queryseq_node->query_alnSegs.size()-1)->opflag, queryseq_node->query_alnSegs.at(queryseq_node->query_alnSegs.size()-1)->seglen);
				overlap_flag_right_end = overlap_flag = isOverlappedPos(queryseq_node->query_alnSegs.at(queryseq_node->query_alnSegs.size()-1)->startRpos, end_ref_pos, var_startRefPos_tmp, var_endRefPos_tmp);
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
					qc_sig->reg_contain_flag = qc_sig->have_next_clip = qc_sig->merge_flag = qc_sig->inv_flag = false;
					qc_sig->altseq_complete_flag = true;
					//qc_sig->ref_pos = queryseq_node->query_alnSegs.at(queryseq_node->query_alnSegs.size()-1)->startRpos;
					//qc_sig->query_pos = queryseq_node->query_alnSegs.at(queryseq_node->query_alnSegs.size()-1)->startQpos;
//					if(queryseq_node->query_alnSegs.at(queryseq_node->query_alnSegs.size()-1)->opflag==BAM_CHARD_CLIP) qc_sig->query_pos = queryseq_node->clip_aln->leftClipSize + queryseq_node->query_alnSegs.at(queryseq_node->query_alnSegs.size()-1)->startQpos;
//					else qc_sig->query_pos = queryseq_node->query_alnSegs.at(queryseq_node->query_alnSegs.size()-1)->startQpos;
					qc_sig->ref_pos = queryseq_node->clip_aln->endRefPos;
					if(queryseq_node->clip_aln->aln_orient==ALN_PLUS_ORIENT) qc_sig->query_pos = queryseq_node->clip_aln->endQueryPos;
					else qc_sig->query_pos = queryseq_node->clip_aln->querylen - queryseq_node->clip_aln->endQueryPos + 1;
					qc_sig->chrname_next_clip = "";
					qc_sig->dist_next_clip = 0;
					qc_sig->have_next_clip = false;
					qc_sig->mate_qcSig = NULL;
					qcSig_vec.push_back(qc_sig);
					last_clip_sig = qc_sig;

					if(qc_sig->cigar_op==BAM_CSOFT_CLIP or qc_sig->cigar_op==BAM_CHARD_CLIP){
						adjClipAlnSegInfo = getAdjacentAlnSegClipReg(id, right_clip_end_flag, query_seqs, minClipEndSize);
						mate_arr_idx = adjClipAlnSegInfo.at(0);
						mate_clip_end_flag = adjClipAlnSegInfo.at(1);
						mate_arr_idx_same_chr = adjClipAlnSegInfo.at(4);
						mate_clip_end_flag_same_chr = adjClipAlnSegInfo.at(5);
						if(mate_arr_idx_same_chr!=-1){
							mate_queryseq_node = query_seqs.at(mate_arr_idx_same_chr);
							// get another segment
							id3 = clip_end_flag3 = -1;
							clip_aln3 = NULL;
							if(mate_clip_end_flag_same_chr==RIGHT_END){
								clip_end_flag3 = LEFT_END;
								clip_aln3 = mate_queryseq_node->clip_aln->left_aln;
							}else{
								clip_end_flag3 = RIGHT_END;
								clip_aln3 = mate_queryseq_node->clip_aln->right_aln;
							}
							if(clip_aln3){
								for(i=0; i<query_seqs.size(); i++){
									if(clip_aln3==query_seqs.at(i)->clip_aln and seq_id_arr[i]==0){
										id3 = i;
										break;
									}
								}
							}

							if(queryseq_node->clip_aln->aln_orient!=mate_queryseq_node->clip_aln->aln_orient){ // different orientation
								if(clip_aln3==NULL){
									end_ref_pos = getEndRefPosAlnSeg(queryseq_node->query_alnSegs.at(0)->startRpos, queryseq_node->query_alnSegs.at(0)->opflag, queryseq_node->query_alnSegs.at(0)->seglen);
									overlap_flag_left_end = isOverlappedPos(queryseq_node->query_alnSegs.at(0)->startRpos, end_ref_pos, var_startRefPos_tmp, var_endRefPos_tmp);
									if(overlap_flag_left_end){ // both ends overlap, then queryseq_node is the main inverted node
										if(mate_clip_end_flag_same_chr==LEFT_END){
											start_ref_pos_var = mate_queryseq_node->clip_aln->startRefPos;
											//start_query_pos_var = mate_queryseq_node->clip_aln->startQueryPos;
											if(mate_queryseq_node->clip_aln->aln_orient==ALN_PLUS_ORIENT) start_query_pos_var = mate_queryseq_node->clip_aln->startQueryPos;
											else start_query_pos_var = mate_queryseq_node->clip_aln->querylen - mate_queryseq_node->clip_aln->startQueryPos + 1;
										}else{
											start_ref_pos_var = mate_queryseq_node->clip_aln->endRefPos;
											//start_query_pos_var = mate_queryseq_node->clip_aln->endQueryPos;
											if(mate_queryseq_node->clip_aln->aln_orient==ALN_PLUS_ORIENT) start_query_pos_var = mate_queryseq_node->clip_aln->endQueryPos;
											else start_query_pos_var = mate_queryseq_node->clip_aln->querylen - mate_queryseq_node->clip_aln->endQueryPos + 1;
										}
										end_ref_pos_var = queryseq_node->clip_aln->endRefPos;
										//end_query_pos_var = queryseq_node->clip_aln->endQueryPos;
										if(queryseq_node->clip_aln->aln_orient==ALN_PLUS_ORIENT) end_query_pos_var = queryseq_node->clip_aln->endQueryPos;
										else end_query_pos_var = queryseq_node->clip_aln->querylen - queryseq_node->clip_aln->endQueryPos + 1;

										qc_sig->chrname = queryseq_node->clip_aln->chrname;
										qc_sig->aln_orient = queryseq_node->clip_aln->aln_orient;
									}else{ // else, mate_queryseq_node is the main inverted node
										start_ref_pos_var = queryseq_node->clip_aln->endRefPos;
										//start_query_pos_var = queryseq_node->clip_aln->endQueryPos;
										if(queryseq_node->clip_aln->aln_orient==ALN_PLUS_ORIENT) start_query_pos_var = queryseq_node->clip_aln->endQueryPos;
										else start_query_pos_var = queryseq_node->clip_aln->querylen - queryseq_node->clip_aln->endQueryPos + 1;
										if(mate_clip_end_flag_same_chr==RIGHT_END){
											end_ref_pos_var = mate_queryseq_node->clip_aln->endRefPos;
											//end_query_pos_var = mate_queryseq_node->clip_aln->endQueryPos;
											if(mate_queryseq_node->clip_aln->aln_orient==ALN_PLUS_ORIENT) end_query_pos_var = mate_queryseq_node->clip_aln->endQueryPos;
											else end_query_pos_var = mate_queryseq_node->clip_aln->querylen - mate_queryseq_node->clip_aln->endQueryPos + 1;
										}else{
											end_ref_pos_var = mate_queryseq_node->clip_aln->startRefPos;
											//end_query_pos_var = mate_queryseq_node->clip_aln->startQueryPos;
											if(mate_queryseq_node->clip_aln->aln_orient==ALN_PLUS_ORIENT) end_query_pos_var = mate_queryseq_node->clip_aln->startQueryPos;
											else end_query_pos_var = mate_queryseq_node->clip_aln->querylen - mate_queryseq_node->clip_aln->startQueryPos + 1;
										}

										qc_sig->chrname = mate_queryseq_node->clip_aln->chrname;
										qc_sig->aln_orient = mate_queryseq_node->clip_aln->aln_orient;
									}

									if(start_ref_pos_var>end_ref_pos_var){
										tmp_pos = start_ref_pos_var;
										start_ref_pos_var = end_ref_pos_var;
										end_ref_pos_var = tmp_pos;
										tmp_pos = start_query_pos_var;
										start_query_pos_var = end_query_pos_var;
										end_query_pos_var = tmp_pos;
									}

									qc_sig->cigar_op = CIGAR_OP_INV_CLIP;
									if(start_ref_pos_var<end_ref_pos_var) qc_sig->cigar_op_len = end_ref_pos_var - start_ref_pos_var + 1;
									else qc_sig->cigar_op_len = start_ref_pos_var - end_ref_pos_var + 1;
									//qc_sig->aln_orient = queryseq_node->clip_aln->aln_orient;
									qc_sig->inv_flag = true;
									//last_clip_sig->reg_contain_flag = qc_sig->have_next_clip = qc_sig->merge_flag = false;
									qc_sig->ref_pos = start_ref_pos_var;
									qc_sig->end_ref_pos = end_ref_pos_var;
									qc_sig->query_pos = start_query_pos_var;
									qc_sig->end_query_pos = end_query_pos_var;

//									if(qc_sig->aln_orient==ALN_PLUS_ORIENT){
//										qc_sig->query_pos = queryseq_node->clip_aln->startQueryPos;
//										qc_sig->end_query_pos = queryseq_node->clip_aln->endQueryPos;
//									}else{
//										qc_sig->query_pos = queryseq_node->clip_aln->endQueryPos;
//										qc_sig->end_query_pos = queryseq_node->clip_aln->startQueryPos;
//									}

									seq_id_arr[mate_arr_idx_same_chr] = 1;
//									mate_arr_idx = mate_arr_idx_same_chr;
//									mate_clip_end_flag = mate_clip_end_flag_same_chr;
								}
							}
						}
					}

				//}else break; // deleted on 2026-03-18
				}
			}
		}else{
//			cout << "queryname=" << queryseq_node->qname << endl;

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
						qc_sig->reg_contain_flag = qc_sig->have_next_clip = qc_sig->merge_flag = qc_sig->inv_flag = false;
						qc_sig->altseq_complete_flag = true;
						//qc_sig->ref_pos = queryseq_node->query_alnSegs.at(queryseq_node->query_alnSegs.size()-1)->startRpos;
//						if(queryseq_node->query_alnSegs.at(queryseq_node->query_alnSegs.size()-1)->opflag==BAM_CHARD_CLIP) qc_sig->query_pos = queryseq_node->clip_aln->leftClipSize + queryseq_node->query_alnSegs.at(queryseq_node->query_alnSegs.size()-1)->startQpos;
//						else qc_sig->query_pos = queryseq_node->query_alnSegs.at(queryseq_node->query_alnSegs.size()-1)->startQpos;
						qc_sig->ref_pos = queryseq_node->clip_aln->endRefPos;
						if(queryseq_node->clip_aln->aln_orient==ALN_PLUS_ORIENT) qc_sig->query_pos = queryseq_node->clip_aln->endQueryPos;
						else qc_sig->query_pos = queryseq_node->clip_aln->querylen - queryseq_node->clip_aln->endQueryPos + 1;
						qc_sig->chrname_next_clip = "";
						qc_sig->dist_next_clip = 0;
						qc_sig->have_next_clip = false;
						qc_sig->mate_qcSig = NULL;
						qcSig_vec.push_back(qc_sig);
					}
				}

				// indel
				qcSig_vec_indel = extractQcSigsFromAlnSegsSingleQuery(queryseq_node, queryseq_node->clip_aln->chrname, var_startRefPos_tmp, var_endRefPos_tmp, MIN_SV_SIZE_CLIP_REG); // min_sv_size
				rmNeighborFalseSig(qcSig_vec_indel, MIN_AVER_SIZE_ALN_SEG, MIN_SIZE_RATIO_MATCH_CLIP_POS);
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
						if(queryseq_node->clip_aln->leftHardClippedFlag==false) qc_sig->altseq = queryseq_node->seq.substr(qc_sig->query_pos, 1);
						else qc_sig->altseq = queryseq_node->seq.substr(qc_sig->query_pos-queryseq_node->clip_aln->leftClipSize, 1);
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
						qc_sig->reg_contain_flag = qc_sig->have_next_clip = qc_sig->merge_flag = qc_sig->inv_flag = false;
						qc_sig->altseq_complete_flag = true;
						//qc_sig->ref_pos = queryseq_node->query_alnSegs.at(0)->startRpos;
						qc_sig->ref_pos = queryseq_node->clip_aln->startRefPos;
//						end_query_pos = getEndQueryPosAlnSeg(queryseq_node->query_alnSegs.at(0)->startQpos, queryseq_node->query_alnSegs.at(0)->opflag, queryseq_node->query_alnSegs.at(0)->seglen);
//						if(queryseq_node->query_alnSegs.at(0)->opflag==BAM_CSOFT_CLIP) qc_sig->query_pos = end_query_pos;
//						else if(queryseq_node->query_alnSegs.at(0)->opflag==BAM_CHARD_CLIP) qc_sig->query_pos = queryseq_node->clip_aln->leftClipSize + end_query_pos;
//						else qc_sig->query_pos = queryseq_node->query_alnSegs.at(0)->startQpos;
						//qc_sig->query_pos = queryseq_node->query_alnSegs.at(0)->startQpos;
						if(queryseq_node->clip_aln->aln_orient==ALN_PLUS_ORIENT) qc_sig->query_pos = queryseq_node->clip_aln->startQueryPos;
						else qc_sig->query_pos = queryseq_node->clip_aln->querylen - queryseq_node->clip_aln->startQueryPos + 1;
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

		// unify the TRA and indels
		//if(queryseq_node->clip_aln->rightmost_flag) break; // deleted on 2025-08-12
		if(id==rightmost_id) break;
		else{
			ignore_end_flag = -1;
			adjClipAlnSegInfo.clear();
			mate_arr_idx = -1;
			mate_clip_end_flag = -1;
			id3 = -1;
			clip_aln3 = NULL;
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
				ref_skip_flag = query_skip_flag = inv_flag = false;
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
						self_overlap_flag = isSegSelfOverlap(queryseq_node->clip_aln, mate_queryseq_node->clip_aln, maxVarRegSize);
						dist = min_dist_same_chr - min_ref_dist_same_chr;
						if(increase_direction==1) dist = -dist;
						if(queryseq_node->clip_aln->aln_orient==mate_queryseq_node->clip_aln->aln_orient) same_orient_flag = true;
						else same_orient_flag = false;

						overlap_flag = isOverlappedPos(mate_queryseq_node->clip_aln->startRefPos, mate_queryseq_node->clip_aln->endRefPos, var_startRefPos_tmp, var_endRefPos_tmp);

						// prefer the segments on the same chromosome
						//if((mate_arr_idx!=mate_arr_idx_same_chr and (abs(min_ref_dist_same_chr)<=MAX_REF_DIST_SAME_CHR or self_overlap_flag)) or (mate_arr_idx==mate_arr_idx_same_chr and abs(min_dist_same_chr)>MAX_REF_DIST_SAME_CHR and abs(dist)>MAX_REF_DIST_SAME_CHR)){ // different chromosomes with near reference distance, deleted on 2024-02-25
						//if((abs(min_ref_dist_same_chr)<=MAX_REF_DIST_SAME_CHR or self_overlap_flag) or (mate_arr_idx==mate_arr_idx_same_chr and abs(min_dist_same_chr)>MAX_REF_DIST_SAME_CHR and abs(dist)>MAX_REF_DIST_SAME_CHR)){ // different chromosomes with near reference distance, deleted on 2025-08-13
						//if((abs(min_ref_dist_same_chr)<=MAX_REF_DIST_SAME_CHR or self_overlap_flag) or (mate_arr_idx==mate_arr_idx_same_chr and ((abs(min_dist_same_chr)>MAX_REF_DIST_SAME_CHR and abs(dist)>MAX_REF_DIST_SAME_CHR) or (overlap_flag and min_ref_dist_same_chr<0 and mate_queryseq_node->clip_aln->querylen<abs(min_ref_dist_same_chr))))){ // DUP, different chromosomes with near reference distance
						if(same_orient_flag and (abs(dist)>MAX_REF_DIST_SAME_CHR and ((abs(min_ref_dist_same_chr)<=MAX_REF_DIST_SAME_CHR or self_overlap_flag) or (mate_arr_idx==mate_arr_idx_same_chr and abs(min_dist_same_chr)>MAX_REF_DIST_SAME_CHR) or (overlap_flag and ((min_ref_dist_same_chr<0 and increase_direction==0) or (min_ref_dist_same_chr>0 and increase_direction==1)) and (self_overlap_flag==false or mate_queryseq_node->clip_aln->querylen<abs(min_ref_dist_same_chr)))))){ // DUP, different chromosomes with near reference distance
							mate_arr_idx = mate_arr_idx_same_chr;
							mate_clip_end_flag = mate_clip_end_flag_same_chr;
							query_skip_flag = true;
						}else if(same_orient_flag and (mate_arr_idx==mate_arr_idx_same_chr and mate_clip_end_flag==mate_clip_end_flag_same_chr and self_overlap_flag==false and abs(min_dist_same_chr)<=MAX_REF_DIST_SAME_CHR and abs(min_ref_dist_same_chr)>MAX_REF_DIST_SAME_CHR)){
							mate_arr_idx = mate_arr_idx_same_chr;
							mate_clip_end_flag = mate_clip_end_flag_same_chr;
							ref_skip_flag = true;
						}else{
							//if(abs(dist)>MAX_REF_DIST_SAME_CHR and abs(min_dist_same_chr)>=MAX_REF_DIST_SAME_CHR and abs(min_dist_same_chr)<=maxVarRegSize and abs(min_ref_dist_same_chr)>MAX_REF_DIST_SAME_CHR and abs(min_ref_dist_same_chr)<=maxVarRegSize){ // deleted on 2026-03-19
							if((abs(dist)>=MAX_REF_DIST_SAME_CHR or abs(min_dist_same_chr)>=MAX_REF_DIST_SAME_CHR) and abs(min_dist_same_chr)<=maxVarRegSize and abs(min_ref_dist_same_chr)>MAX_REF_DIST_SAME_CHR and abs(min_ref_dist_same_chr)<=maxVarRegSize){
								if(same_orient_flag){
									if(abs(min_dist)>MAX_REF_DIST_SAME_CHR and abs(min_ref_dist)>MAX_REF_DIST_SAME_CHR){
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
								}else{ // INV
									// get another segment
									if(mate_clip_end_flag==RIGHT_END){
										clip_end_flag3 = LEFT_END;
										clip_aln3 = mate_queryseq_node->clip_aln->left_aln;
									}else{
										clip_end_flag3 = RIGHT_END;
										clip_aln3 = mate_queryseq_node->clip_aln->right_aln;
									}
									if(clip_aln3){
										for(i=0; i<query_seqs.size(); i++){
											if(clip_aln3==query_seqs.at(i)->clip_aln and seq_id_arr[i]==0){
												id3 = i;
												break;
											}
										}

										// if found, then update the variant signature
										if(id3!=-1 and queryseq_node->clip_aln->aln_orient==clip_aln3->aln_orient){
											inv_flag = true;
										}
									}
								}
							}
						}
					}
				}

				if(last_clip_sig){
					if(query_skip_flag){ // INS
						if(increase_direction==1){ // decrease direction, then update the clip position
							if(mate_clip_end_flag_same_chr==RIGHT_END){ // right end
								last_clip_sig->ref_pos = mate_queryseq_node->clip_aln->endRefPos;
								if(mate_queryseq_node->clip_aln->aln_orient==ALN_PLUS_ORIENT) last_clip_sig->query_pos = mate_queryseq_node->clip_aln->endQueryPos;
								else last_clip_sig->query_pos = mate_queryseq_node->clip_aln->querylen - mate_queryseq_node->clip_aln->endQueryPos + 1;
							}else{ // left end
								last_clip_sig->ref_pos = mate_queryseq_node->clip_aln->startRefPos;
								if(mate_queryseq_node->clip_aln->aln_orient==ALN_PLUS_ORIENT) last_clip_sig->query_pos = mate_queryseq_node->clip_aln->startQueryPos;
								else last_clip_sig->query_pos = mate_queryseq_node->clip_aln->querylen - mate_queryseq_node->clip_aln->startQueryPos + 1;
							}
						}

						last_clip_sig->dist_next_clip = last_clip_sig->cigar_op_len = abs(dist) + 1;
						last_clip_sig->cigar_op = BAM_CINS;
						//end_ref_pos = getEndRefPosAlnSeg(queryseq_node->query_alnSegs.at(0)->startRpos, queryseq_node->query_alnSegs.at(0)->opflag, queryseq_node->query_alnSegs.at(0)->seglen);
						reg_str = last_clip_sig->chrname + ":" + to_string(last_clip_sig->ref_pos) + "-" + to_string(last_clip_sig->ref_pos);
						pthread_mutex_lock(&mutex_fai);
						p_seq = fai_fetch(fai, reg_str.c_str(), &seq_len);
						pthread_mutex_unlock(&mutex_fai);
						last_clip_sig->refseq = p_seq;
						free(p_seq);

						//last_clip_sig->altseq = queryseq_node->seq.substr(last_clip_sig->query_pos-1, last_clip_sig->cigar_op_len);
//						if(last_clip_sig->aln_orient==query_orient){
							if(last_clip_sig->query_pos+last_clip_sig->cigar_op_len-1<=(int64_t)queryseq_whole.size())
								last_clip_sig->altseq = queryseq_whole.substr(last_clip_sig->query_pos-1, last_clip_sig->cigar_op_len);
							else{
								last_clip_sig->altseq = queryseq_whole.substr(last_clip_sig->query_pos-1);
								last_clip_sig->altseq_complete_flag = false;
							}
//						}else{
//							//cout << "line=" << __LINE__ << ", " << chrname << ":" << var_startRefPos << "-" << var_endRefPos << ": " << queryseq_node->qname << ", aln_orient is not the same!" << endl;
//
//							end_query_pos = getEndQueryPosAlnSeg(last_clip_sig->query_pos, last_clip_sig->cigar_op, last_clip_sig->cigar_op_len);
//							if(end_query_pos-1<(int64_t)queryseq_whole.size()){
//								start_qpos_tmp = queryseq_whole.size() - end_query_pos;
//								reverse_seq = queryseq_whole.substr(start_qpos_tmp-1, last_clip_sig->cigar_op_len);
//								reverseComplement(reverse_seq);
//								last_clip_sig->altseq = reverse_seq;
//							}else valid_flag = false;
//						}
//						if(valid_flag==false) break;

					}else if(ref_skip_flag){ // DEL
						if(increase_direction==1){ // decrease direction, then update the clip position
							if(mate_clip_end_flag_same_chr==RIGHT_END){ // right end
								last_clip_sig->ref_pos = mate_queryseq_node->clip_aln->endRefPos;
								if(mate_queryseq_node->clip_aln->aln_orient==ALN_PLUS_ORIENT) last_clip_sig->query_pos = mate_queryseq_node->clip_aln->endQueryPos;
								else last_clip_sig->query_pos = mate_queryseq_node->clip_aln->querylen - mate_queryseq_node->clip_aln->endQueryPos + 1;
							}else{ // left end
								last_clip_sig->ref_pos = mate_queryseq_node->clip_aln->startRefPos;
								if(mate_queryseq_node->clip_aln->aln_orient==ALN_PLUS_ORIENT) last_clip_sig->query_pos = mate_queryseq_node->clip_aln->startQueryPos;
								else last_clip_sig->query_pos = mate_queryseq_node->clip_aln->querylen - mate_queryseq_node->clip_aln->startQueryPos + 1;
							}
						}

						last_clip_sig->dist_next_clip = last_clip_sig->cigar_op_len = abs(dist) + 1;
						last_clip_sig->cigar_op = BAM_CDEL;

						reg_str = last_clip_sig->chrname + ":" + to_string(last_clip_sig->ref_pos) + "-" + to_string(last_clip_sig->ref_pos+last_clip_sig->cigar_op_len-1);
						pthread_mutex_lock(&mutex_fai);
						p_seq = fai_fetch(fai, reg_str.c_str(), &seq_len);
						pthread_mutex_unlock(&mutex_fai);
						last_clip_sig->refseq = p_seq;
						free(p_seq);

						//last_clip_sig->altseq = queryseq_node->seq.substr(last_clip_sig->query_pos-1, 1);
//						if(last_clip_sig->aln_orient==query_orient){
							if(last_clip_sig->query_pos<=(int64_t)queryseq_whole.size())
								last_clip_sig->altseq = queryseq_whole.substr(last_clip_sig->query_pos-1, 1);
							else valid_flag = false; // deleted on 2025-08-13
//						}else{
//							//cout << "line=" << __LINE__ << ", " << chrname << ":" << var_startRefPos << "-" << var_endRefPos << ": " << queryseq_node->qname << ", aln_orient is not the same!" << endl;
//
//							end_query_pos = getEndQueryPosAlnSeg(last_clip_sig->query_pos, last_clip_sig->cigar_op, last_clip_sig->cigar_op_len);
//							if(end_query_pos-1<(int64_t)queryseq_whole.size()){
//								start_qpos_tmp = queryseq_whole.size() - end_query_pos;
//								reverse_seq = queryseq_whole.substr(start_qpos_tmp-1, 1);
//								reverseComplement(reverse_seq);
//								last_clip_sig->altseq = reverse_seq;
//							}else valid_flag = false;
//						}
						if(valid_flag==false) break;
					}else if(inv_flag){ // INV
						// update the variant signature
						if(last_clip_sig->cigar_op==BAM_CSOFT_CLIP or last_clip_sig->cigar_op==BAM_CHARD_CLIP){
							last_clip_sig->cigar_op = CIGAR_OP_INV_CLIP;
							last_clip_sig->chrname = mate_queryseq_node->clip_aln->chrname;
							last_clip_sig->cigar_op_len = mate_queryseq_node->clip_aln->ref_dist;
							last_clip_sig->aln_orient = mate_queryseq_node->clip_aln->aln_orient;
							last_clip_sig->inv_flag = true;
							//last_clip_sig->reg_contain_flag = qc_sig->have_next_clip = qc_sig->merge_flag = false;
							last_clip_sig->ref_pos = mate_queryseq_node->clip_aln->startRefPos;
							last_clip_sig->end_ref_pos = mate_queryseq_node->clip_aln->endRefPos;
							if(last_clip_sig->aln_orient==ALN_PLUS_ORIENT){
								last_clip_sig->query_pos = mate_queryseq_node->clip_aln->startQueryPos;
								last_clip_sig->end_query_pos = mate_queryseq_node->clip_aln->endQueryPos;
							}else{
								//last_clip_sig->query_pos = mate_queryseq_node->clip_aln->endQueryPos;
								//last_clip_sig->end_query_pos = mate_queryseq_node->clip_aln->startQueryPos;
								last_clip_sig->query_pos = mate_queryseq_node->clip_aln->querylen - mate_queryseq_node->clip_aln->endQueryPos + 1;
								last_clip_sig->end_query_pos = mate_queryseq_node->clip_aln->querylen - mate_queryseq_node->clip_aln->startQueryPos + 1;
							}
							last_clip_sig->chrname_next_clip = "";
							last_clip_sig->dist_next_clip = 0;
							last_clip_sig->have_next_clip = false;
							last_clip_sig->mate_qcSig = NULL;

							seq_id_arr[mate_arr_idx] = 1;
							//seq_id_arr[id3] = 1;
							mate_arr_idx = id3;
							mate_clip_end_flag = clip_end_flag3;
						}
					}
				}else{
					query_skip_flag = ref_skip_flag = inv_flag = false;
				}

				id = mate_arr_idx;
				end_flag = mate_clip_end_flag;
				if(query_skip_flag or ref_skip_flag or inv_flag) ignore_end_flag = mate_clip_end_flag;

				if(seq_id_arr[id]==1) break;
			}else
				break;
		}
	}

	if(valid_flag){
		// merge nearest signatures for non-ccs data, deleted on 2026-03-20
//		if(technology.compare(PACBIO_CCS_TECH_STR)!=0){
//			for(i=1; i<qcSig_vec.size(); ){
//				qc_sig = qcSig_vec.at(i-1);
//				qc_sig_next = qcSig_vec.at(i);
//				merge_flag = false;
//				if(qc_sig->cigar_op==qc_sig_next->cigar_op and qc_sig->aln_orient==qc_sig_next->aln_orient and qc_sig->cigar_op!=BAM_CSOFT_CLIP and qc_sig_next->cigar_op!=BAM_CHARD_CLIP){
//					end_ref_pos = getEndRefPosAlnSeg(qc_sig->ref_pos, qc_sig->cigar_op, qc_sig->cigar_op_len);
//					end_ref_pos2 = getEndRefPosAlnSeg(qc_sig_next->ref_pos, qc_sig_next->cigar_op, qc_sig_next->cigar_op_len);
//					//end_query_pos = getEndQueryPosAlnSeg(qc_sig->query_pos, qc_sig->cigar_op, qc_sig->cigar_op_len);
//					end_query_pos2 = getEndQueryPosAlnSeg(qc_sig_next->query_pos, qc_sig_next->cigar_op, qc_sig_next->cigar_op_len);
//					overlap_flag = isOverlappedPos(qc_sig->ref_pos-MAX_DIST_MATCH_INDEL, end_ref_pos+MAX_DIST_MATCH_INDEL, qc_sig_next->ref_pos-MAX_DIST_MATCH_INDEL, end_ref_pos2+MAX_DIST_MATCH_INDEL);
//					if(overlap_flag){
//						dist = abs(end_query_pos2 - qc_sig->query_pos - (end_ref_pos2 - qc_sig->ref_pos));
//						qc_sig->cigar_op_len = dist + 1;
//
//						if(end_ref_pos2<qc_sig->ref_pos) end_ref_pos2 = qc_sig->ref_pos;
//						reg_str = qc_sig->chrname + ":" + to_string(qc_sig->ref_pos) + "-" + to_string(end_ref_pos2);
//						pthread_mutex_lock(&mutex_fai);
//						p_seq = fai_fetch(fai, reg_str.c_str(), &seq_len);
//						pthread_mutex_unlock(&mutex_fai);
//						qc_sig->refseq = p_seq;
//						free(p_seq);
//
//						dist = abs(end_query_pos2 - qc_sig->query_pos);
//						if(qc_sig->aln_orient==query_orient)
//							qc_sig->altseq = queryseq_whole.substr(qc_sig->query_pos-1, dist+1);
//						else{
//							//cout << "line=" << __LINE__ << ", " << chrname << ":" << var_startRefPos << "-" << var_endRefPos << ": " << queryseq_node->qname << ", aln_orient is not the same!" << endl;
//
//							end_query_pos = getEndQueryPosAlnSeg(qc_sig->query_pos, qc_sig->cigar_op, qc_sig->cigar_op_len);
//							start_qpos_tmp = queryseq_whole.size() - end_query_pos;
//							reverse_seq = queryseq_whole.substr(start_qpos_tmp-1, dist+1);
//							reverseComplement(reverse_seq);
//							qc_sig->altseq = reverse_seq;
//						}
//
//						delete qc_sig_next;
//						qcSig_vec.erase(qcSig_vec.begin()+i);
//
//						merge_flag = true;
//					}
//				}
//				if(merge_flag==false) i++;
//			}
//		}

		// sort according to queryPos
		if(qcSig_vec.size()>=2){
			for(i=0; i<qcSig_vec.size(); i++){
				qc_sig = qcSig_vec.at(i);
				min_idx = -1;
				min_val = LONG_MAX;
				for(j=i; j<qcSig_vec.size(); j++){
					qc_sig_next = qcSig_vec.at(j);
					if(qc_sig_next->query_pos<min_val){
						min_idx = j;
						min_val = qc_sig_next->query_pos;
					}
				}
				if(min_idx!=(int64_t)i){ // swap
					last_clip_sig = qcSig_vec.at(i);
					qcSig_vec.at(i) = qcSig_vec.at(min_idx);
					qcSig_vec.at(min_idx) = last_clip_sig;
				}
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
		qcSigList_node->startSpanRefPos = sidemost_idx_vec.at(4);
		qcSigList_node->endSpanRefPos = sidemost_idx_vec.at(5);
		qcSigList_node->query_seqs_vec = query_seqs;
		qcSigList_node->qcSig_vec = qcSig_vec;
		qcSigList_node->ins_sum = qcSigList_node->del_sum = 0;
		qcSigList_node->total_valid_seg_num = qcSigList_node->seg_num_aln_to_single_reg = 0;
		qcSigList_node->cluster_finished_flag = qcSigList_node->collapse_to_single_reg_flag = false;
		qcSigList_node->seg_num_aln_to_single_reg = 0;

		var_startRefPos_tmp = var_startRefPos + EXT_SIZE_CHK_VAR_LOC;
		if(var_startRefPos_tmp<1) var_startRefPos_tmp = 1;
		var_endRefPos_tmp = var_endRefPos - EXT_SIZE_CHK_VAR_LOC;
		if(var_endRefPos_tmp>chrlen) var_endRefPos_tmp = chrlen;
		if(qcSigList_node->startSpanRefPos<=var_startRefPos_tmp and qcSigList_node->endSpanRefPos>=var_endRefPos_tmp) qcSigList_node->entire_flanking_flag = true;
		else qcSigList_node->entire_flanking_flag = false;

		if(query_seqs.size()==1) qcSigList_node->single_aln_seg_flag = true;
		else qcSigList_node->single_aln_seg_flag = false;

		return qcSigList_node;
	}else{
		for(i=1; i<qcSig_vec.size(); i++) delete qcSig_vec.at(i);
		vector<qcSig_t*>().swap(qcSig_vec);
		return NULL;
	}
}

// get the left and right most align segments: [0]-minQposIdx, [1]-maxQposIdx, [2]-minQpos, [3]-maxQpos, [4]-minRpos, [5]-maxRpos.
// Updated on 2025-08-12
vector<int32_t> clipRegCluster::getLeftRightMostAlnSegIdx(vector<struct querySeqInfoNode*> &query_seqs){
	vector<int32_t> ret_idx_vec;
	size_t i;
	int32_t minQpos_id, maxQpos_id;
	int64_t minQpos, maxQpos, minRpos, maxRpos;
	struct querySeqInfoNode *query_seq_node;

	minQpos_id = maxQpos_id = -1;
	minQpos = minRpos = LONG_MAX;
	maxQpos = maxRpos = LONG_MIN;
	for(i=0; i<query_seqs.size(); i++){
		query_seq_node = query_seqs.at(i);
//		if(query_seq_node->clip_aln->leftmost_flag) leftmost_id = i;
//		if(query_seq_node->clip_aln->rightmost_flag) rightmost_id = i;
		if(query_seq_node->clip_aln->aln_orient==ALN_PLUS_ORIENT){
			if(minQpos>query_seq_node->clip_aln->startQueryPos) {
				minQpos = query_seq_node->clip_aln->startQueryPos;
				minQpos_id = i;
			}
			if(maxQpos<query_seq_node->clip_aln->endQueryPos) {
				maxQpos = query_seq_node->clip_aln->endQueryPos;
				maxQpos_id = i;
			}
		}else if(query_seq_node->clip_aln->aln_orient==ALN_MINUS_ORIENT){
			if(minQpos>query_seq_node->clip_aln->endQueryPos) {
				minQpos = query_seq_node->clip_aln->endQueryPos;
				minQpos_id = i;
			}
			if(maxQpos<query_seq_node->clip_aln->startQueryPos) {
				maxQpos = query_seq_node->clip_aln->startQueryPos;
				maxQpos_id = i;
			}
		}
		if(query_seq_node->clip_aln->startRefPos<minRpos) minRpos = query_seq_node->clip_aln->startRefPos;
		if(query_seq_node->clip_aln->endRefPos>maxRpos) maxRpos = query_seq_node->clip_aln->endRefPos;
	}

	ret_idx_vec.push_back(minQpos_id);
	ret_idx_vec.push_back(maxQpos_id);
	ret_idx_vec.push_back(minQpos);
	ret_idx_vec.push_back(maxQpos);
	ret_idx_vec.push_back(minRpos);
	ret_idx_vec.push_back(maxRpos);

	return ret_idx_vec;
}

/**
 * get adjacent clip segment according to query positions
 * @return: the adjacent clip segment information
 *  [0]: array index of minimal distance;
 *  [1]: segment end flag;
 *  [2]: minimal query distance;
 *  [3]: minimal reference distance;
 *  [4]: array index of minimal distance in the same chromosome;
 *  [5]: segment end flag in the same chromosome;
 *  [6]: minimal query distance in the same chromosome;
 *  [7]: minimal reference distance in the same chromosome;
 *  [8]: rpos increase direction: -1 for unused, 0 for increase, 1 for decrease
 */
vector<int32_t> clipRegCluster::getAdjacentAlnSegClipReg(int32_t arr_idx, int32_t clip_end_flag, vector<struct querySeqInfoNode*> &query_seqs, int32_t minClipEndSize){
	vector<int32_t> adjClipAlnSegInfo; // [0]: array index of minimal distance; [1]: segment end flag
	vector<clipAlnData_t*> query_aln_segs;
	size_t i;

	for(i=0; i<query_seqs.size(); i++) if(query_seqs.at(i)->clip_aln) query_aln_segs.push_back(query_seqs.at(i)->clip_aln);
	adjClipAlnSegInfo = getAdjacentClipAlnSeg(arr_idx, clip_end_flag, query_aln_segs, minClipEndSize, maxVarRegSize, minMapQ);

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
				//matchFlag = isQcSigMatchClipReg(qc_sig1, qc_sig2, MAX_DIST_MATCH_CLIP_POS, MIN_SIZE_RATIO_MATCH_CLIP_POS, QC_SIZE_RATIO_MATCH_THRES_INDEL, QC_SEQSIM_RATIO_MATCH_THRES);
				matchFlag = isQcSigMatchClipReg(qc_sig1, qc_sig2, var_startRefPos, var_endRefPos, MAX_DIST_MATCH_CLIP_POS, MIN_SIZE_IGNORE_SEQSIM, MIN_SIZE_RATIO_MATCH_CLIP_POS, QC_SIZE_RATIO_MATCH_THRES_INDEL, min_seqsim_match);
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

bool clipRegCluster::isQcSigMatchClipReg(qcSig_t *qc_sig, qcSig_t *seed_qc_sig, int64_t var_startRefPos, int64_t var_endRefPos, int64_t max_dist_match_clip_pos, int32_t min_size_ignore_seqsim, double min_size_ratio_match_thres_clip, double size_ratio_match_thres, double seqsim_ratio_match_thres){
	bool match_flag = false;
	double size_ratio, overlap_ratio1, overlap_ratio2;
	int64_t ref_dist, ref_dist2, overlap_size;

	if(qc_sig and seed_qc_sig){
		if((qc_sig->cigar_op==BAM_CINS and seed_qc_sig->cigar_op==BAM_CINS) or (qc_sig->cigar_op==BAM_CDEL and seed_qc_sig->cigar_op==BAM_CDEL)){ // indel
			match_flag = isQcSigMatch(qc_sig, seed_qc_sig, var_startRefPos, var_endRefPos, max_dist_match_clip_pos, min_size_ignore_seqsim, size_ratio_match_thres, seqsim_ratio_match_thres, fai);
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
		}else if(qc_sig->inv_flag and seed_qc_sig->inv_flag and qc_sig->cigar_op==CIGAR_OP_INV_CLIP and seed_qc_sig->cigar_op==CIGAR_OP_INV_CLIP){
			if(qc_sig->chrname.compare(seed_qc_sig->chrname)==0){
				ref_dist = qc_sig->ref_pos - seed_qc_sig->ref_pos;
				ref_dist2 = qc_sig->end_ref_pos - seed_qc_sig->end_ref_pos;
				overlap_size = getOverlapSize(qc_sig->ref_pos, qc_sig->end_ref_pos, seed_qc_sig->ref_pos, seed_qc_sig->end_ref_pos);
				overlap_ratio1 = (double)overlap_size / (qc_sig->end_ref_pos - qc_sig->ref_pos + 1);
				overlap_ratio2 = (double)overlap_size / (seed_qc_sig->end_ref_pos - seed_qc_sig->ref_pos + 1);
				if(qc_sig->cigar_op_len<seed_qc_sig->cigar_op_len) size_ratio = (double)qc_sig->cigar_op_len / seed_qc_sig->cigar_op_len;
				else size_ratio = (double)seed_qc_sig->cigar_op_len / qc_sig->cigar_op_len;
				if((abs(ref_dist)<=max_dist_match_clip_pos or abs(ref_dist2)<=max_dist_match_clip_pos or (overlap_ratio1>=min_size_ratio_match_thres_clip and overlap_ratio2>=min_size_ratio_match_thres_clip)) and size_ratio>=min_size_ratio_match_thres_clip)
					match_flag = true;
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
			if(queryCluSig->qcSig_vec.at(i-1)->cigar_op==BAM_CINS or queryCluSig->qcSig_vec.at(i-1)->cigar_op==BAM_CDEL or queryCluSig->qcSig_vec.at(i-1)->cigar_op==CIGAR_OP_INV_CLIP)
				sig_str1 = to_string(queryCluSig->qcSig_vec.at(i-1)->ref_pos) + "_" + to_string(queryCluSig->qcSig_vec.at(i-1)->cigar_op_len) + "_" + to_string(queryCluSig->qcSig_vec.at(i-1)->cigar_op);
			else
				sig_str1 = to_string(queryCluSig->qcSig_vec.at(i-1)->ref_pos) + "_" + to_string(queryCluSig->qcSig_vec.at(i-1)->dist_next_clip) + "_" + to_string(queryCluSig->qcSig_vec.at(i-1)->cigar_op);
			if(seed_qcQuery->qcSig_vec.at(j-1)->cigar_op==BAM_CINS or seed_qcQuery->qcSig_vec.at(j-1)->cigar_op==BAM_CDEL or seed_qcQuery->qcSig_vec.at(j-1)->cigar_op==CIGAR_OP_INV_CLIP)
				sig_str2 = to_string(seed_qcQuery->qcSig_vec.at(j-1)->ref_pos) + "_" + to_string(seed_qcQuery->qcSig_vec.at(j-1)->cigar_op_len) + "_" + to_string(seed_qcQuery->qcSig_vec.at(j-1)->cigar_op);
			else
				sig_str2 = to_string(seed_qcQuery->qcSig_vec.at(j-1)->ref_pos) + "_" + to_string(seed_qcQuery->qcSig_vec.at(j-1)->dist_next_clip) + "_" + to_string(seed_qcQuery->qcSig_vec.at(j-1)->cigar_op);

			sigseq_vec1.push_back(sig_str1);
			sigseq_vec2.push_back(sig_str2);

			i--;
			j--;
		}
		else if(value==1){ //from (i-1, j)
			//match_profile_vec.push_back(0);
			if(queryCluSig->qcSig_vec.at(i-1)->cigar_op==BAM_CINS or queryCluSig->qcSig_vec.at(i-1)->cigar_op==BAM_CDEL or queryCluSig->qcSig_vec.at(i-1)->cigar_op==CIGAR_OP_INV_CLIP)
				sig_str1 = to_string(queryCluSig->qcSig_vec.at(i-1)->ref_pos) + "_" + to_string(queryCluSig->qcSig_vec.at(i-1)->cigar_op_len) + "_" + to_string(queryCluSig->qcSig_vec.at(i-1)->cigar_op);
			else
				sig_str1 = to_string(queryCluSig->qcSig_vec.at(i-1)->ref_pos) + "_" + to_string(queryCluSig->qcSig_vec.at(i-1)->dist_next_clip) + "_" + to_string(queryCluSig->qcSig_vec.at(i-1)->cigar_op);
			sigseq_vec1.push_back(sig_str1);
			sigseq_vec2.push_back("-");

			i--;
		}
		else{ //from (i, j-1)
			//match_profile_vec.push_back(0);
			if(seed_qcQuery->qcSig_vec.at(j-1)->cigar_op==BAM_CINS or seed_qcQuery->qcSig_vec.at(j-1)->cigar_op==BAM_CDEL or seed_qcQuery->qcSig_vec.at(j-1)->cigar_op==CIGAR_OP_INV_CLIP)
				sig_str2 = to_string(seed_qcQuery->qcSig_vec.at(j-1)->ref_pos) + "_" + to_string(seed_qcQuery->qcSig_vec.at(j-1)->cigar_op_len) + "_" + to_string(seed_qcQuery->qcSig_vec.at(j-1)->cigar_op);
			else
				sig_str2 = to_string(seed_qcQuery->qcSig_vec.at(j-1)->ref_pos) + "_" + to_string(seed_qcQuery->qcSig_vec.at(j-1)->dist_next_clip) + "_" + to_string(seed_qcQuery->qcSig_vec.at(j-1)->cigar_op);
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

