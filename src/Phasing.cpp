#include "Phasing.h"

Phasing::Phasing(int32_t min_sv_size, int32_t max_sv_size, double size_ratio_thres, double seqsim_thres, vector<varCand*> &var_cand_vec, vector<varCand*> &var_cand_clipReg_vec) {
	this->min_sv_size = min_sv_size;
	this->max_sv_size = max_sv_size;
	this->size_ratio_thres = size_ratio_thres;
	this->seqsim_thres = seqsim_thres;
	this->var_cand_vec = var_cand_vec;
	this->var_cand_clipReg_vec = var_cand_clipReg_vec;
}

Phasing::~Phasing() {
}

void Phasing::releaseHapNodes(vector< vector<phasing_reg_t*> > &hapnodes_vec){
	size_t i, j;
	vector<phasing_reg_t*> hapnode_vec;
	for(i=0; i<hapnodes_vec.size(); i++) {
		hapnode_vec = hapnodes_vec.at(i);
		for(j=0; j<hapnode_vec.size(); j++) delete hapnode_vec.at(j);
		vector<phasing_reg_t*>().swap(hapnode_vec);
	}
	vector< vector<phasing_reg_t*> >().swap(hapnodes_vec);
}

// operation for phasing
void Phasing::performPhasing(){
	bool success_flag, size_satisfied, overlap_flag;
	varCand *var_cand, *var_cand2;
	size_t i, j, k, m;
	phasing_reg_t *phasing_reg;
	vector<phasing_reg_t*> hapnode_vec1, hapnode_vec2;
	vector< vector<phasing_reg_t*> > hapnodes_vec, hapnodes_vec2;
	int32_t end_idx, max_len, secmax_len;
	reg_t *reg, *reg2;
	map<int64_t, int32_t> ps_id_map;
	double size_ratio, val, min_len;

	// get sub-groups
	i = 0;
	while(i<var_cand_vec.size()){
		end_idx = i;
		var_cand = var_cand_vec.at(i);
		if(var_cand->call_success){

#if PHASING_DEBUG
			if(var_cand->newVarVec.at(0)->startRefPos==1108719){
				cout << var_cand->newVarVec.at(0)->startRefPos << endl;
			}
			cout << var_cand->chrname << ":" << var_cand->newVarVec.at(0)->startRefPos << "-" << var_cand->newVarVec.at(var_cand->newVarVec.size()-1)->endRefPos << ", ctg_num=" << var_cand->ctg_num << ", var_num=" << var_cand->newVarVec.size() << endl;
#endif

			// initialize the hap nodes
			if(hapnodes_vec.empty()){
				generateHapNodeSingleVarCand(hapnodes_vec, var_cand);
				if(isValidStartLinkHapNode(hapnodes_vec)==false) { i++; continue; }
			}

			// link neighbors
			for(j=i+1; j<var_cand_vec.size(); j++){
				var_cand2 = var_cand_vec.at(j);
				if(var_cand2->call_success){

#if PHASING_DEBUG
					cout << var_cand2->chrname << ":" << var_cand2->newVarVec.at(0)->startRefPos << "-" << var_cand2->newVarVec.at(var_cand2->newVarVec.size()-1)->endRefPos << ", ctg_num=" << var_cand2->ctg_num << ", var_num=" << var_cand2->newVarVec.size() << endl;
					if(var_cand2->newVarVec.at(0)->startRefPos==1108719){
						cout << var_cand2->newVarVec.at(0)->startRefPos << endl;
					}
#endif

					generateHapNodeSingleVarCand(hapnodes_vec2, var_cand2);
					if(isValidStartLinkHapNode(hapnodes_vec2)==false) continue;

					success_flag = linkHapNodes(hapnodes_vec, hapnodes_vec2);

#if PHASING_DEBUG
					cout << "success_flag=" << success_flag << endl;
					printHapNodes(hapnodes_vec);
#endif

					if(success_flag) {
						end_idx = j + 1;
						releaseHapNodes(hapnodes_vec2); // release memory
					}else {
						end_idx = j;
						break;
					}
				}
			}

			ps_id_map.clear();
			for(j=0; j<hapnodes_vec.size(); j++){
				hapnode_vec1 = hapnodes_vec.at(j);
				if(hapnode_vec1.size()>=2){ // ignore single node
					// remove redundant items
					for(k=1; k<hapnode_vec1.size(); k++){
						reg = hapnode_vec1.at(k-1)->reg;
						reg2 = hapnode_vec1.at(k)->reg;
						if(reg->var_type==reg2->var_type and reg->chrname.compare(reg2->chrname)==0){
							overlap_flag = false;
							if(reg->startRefPos==reg2->startRefPos and reg->endRefPos==reg2->endRefPos) overlap_flag = true;
							else overlap_flag = isOverlappedPos(reg->startRefPos-SHORT_VAR_ALN_CHECK_EXTEND_SIZE, reg->endRefPos+SHORT_VAR_ALN_CHECK_EXTEND_SIZE, reg2->startRefPos-SHORT_VAR_ALN_CHECK_EXTEND_SIZE, reg2->endRefPos+SHORT_VAR_ALN_CHECK_EXTEND_SIZE);

							if(overlap_flag){
								if(reg->sv_len==reg2->sv_len) size_ratio = 1;
								else{
									if(abs(reg->sv_len)>abs(reg2->sv_len)){ max_len = abs(reg->sv_len); secmax_len = abs(reg2->sv_len); }
									else{ max_len = abs(reg2->sv_len); secmax_len = abs(reg->sv_len); }
									size_ratio = secmax_len / max_len;
								}

								if(size_ratio>=size_ratio_thres){
									if(reg->var_type!=VAR_DEL){
										val = computeVarseqSim(reg->altseq, reg2->altseq);
										min_len = min(abs(reg->sv_len), abs(reg2->sv_len));
										if(val>=seqsim_thres or (val>=QC_SEQSIM_RATIO_THRES_FACTOR*seqsim_thres and min_len<MAX_DIST_MERGE_ARBITARY)) hapnode_vec1.at(k)->redundant_flag = true;
									}else{
										hapnode_vec1.at(k)->redundant_flag = true;
									}
								}
							}
						}
					}

					// count
					for(k=0; k<hapnode_vec1.size(); k++){
						phasing_reg = hapnode_vec1.at(k);
						if(phasing_reg->redundant_flag==false){
							//size_satisfied = isSizeSatisfied2(phasing_reg->reg->sv_len, min_sv_size, max_sv_size); // deleted on 2025-09-17
							size_satisfied = isSizeSatisfied2(phasing_reg->reg->sv_len, min_sv_size);
							if(size_satisfied) ps_id_map[phasing_reg->ps_id] ++;
						}
					}
				}
			}

			// process the haplotype information
			for(j=0; j<hapnodes_vec.size(); j++){
				hapnode_vec1 = hapnodes_vec.at(j);
				if(hapnode_vec1.size()>=2){ // ignore single node
					for(k=0; k<hapnode_vec1.size(); k++){
						phasing_reg = hapnode_vec1.at(k);
						if(phasing_reg->redundant_flag==false and ps_id_map[phasing_reg->ps_id]>=2){
							reg = phasing_reg->reg;
							//size_satisfied = isSizeSatisfied2(reg->sv_len, min_sv_size, max_sv_size); // deleted on 2025-09-17
							size_satisfied = isSizeSatisfied2(reg->sv_len, min_sv_size);
							if(size_satisfied){
								//cout << phasing_reg->ps_id << ", " << reg->chrname << ":" << reg->startRefPos << "-" << reg->endRefPos << ", type=" << to_string(reg->var_type) << ", svlen=" << reg->sv_len << ", redundant=" << phasing_reg->redundant_flag << endl;
								if(reg->PS.size()==0) reg->PS = reg->chrname + ":" + to_string(phasing_reg->ps_id);
								else reg->PS += "," + reg->chrname + ":" + to_string(phasing_reg->ps_id);
								for(m=0; m<reg->gt_seq.size(); m++) if(reg->gt_seq.at(m)=='/') reg->gt_seq.at(m) = '|';
							}
						}
					}
				}
			}

#if PHASING_DEBUG
			//cout << "ps_num_tmp=" << ps_num_tmp << ", ps_num=" << ps_num << ", phased_var_num_tmp=" << phased_var_num_tmp << ", phased_var_num=" << phased_var_num << endl;

			cout << "After local extension:" << endl;
			printHapNodes(hapnodes_vec);
#endif

			if(i!=(size_t)end_idx) i = end_idx;
			else i = end_idx + 1;
			releaseHapNodes(hapnodes_vec); // release memory
			hapnodes_vec = hapnodes_vec2;
			hapnodes_vec2.clear();
		}else{
			releaseHapNodes(hapnodes_vec); // release memory
			i++;
		}
	}
	if(!hapnodes_vec.empty()) releaseHapNodes(hapnodes_vec);
	if(!hapnodes_vec2.empty()) releaseHapNodes(hapnodes_vec2);

	// search overlapped mate clip regions


}

// initialize haplotype node
void Phasing::generateHapNodeSingleVarCand(vector< vector<phasing_reg_t*> > &hapnodes_vec, varCand *var_cand){
	int32_t query_id, group_id, id_vec;
	int64_t min_val;
	size_t i, j;
	reg_t *reg;
	phasing_reg_t *phasing_reg;
	vector< vector<string> > qnames_vec_cluster;
	vector<int32_t> query_id_vec, group_id_vec, ps_id_vec;
	vector<phasing_reg_t*> hapnode_vec1, hapnode_vec2;

//	if(hapnodes_vec.empty()){ // initialize the hap node
	if(hapnodes_vec.empty() or (hapnodes_vec.at(0).empty() and hapnodes_vec.at(1).empty())){ // initialize the hap node
		if(hapnodes_vec.empty()){
			hapnodes_vec.push_back(hapnode_vec1);
			hapnodes_vec.push_back(hapnode_vec2);
		}

		qnames_vec_cluster = getClusterInfo(var_cand->clusterfilename);

		// compute the group id for each query
		for(i=0; i<var_cand->newVarVec.size(); i++){
			reg = var_cand->newVarVec.at(i);
			if(reg->var_type!=VAR_UNC and reg->call_success_status and reg->query_id!=-1)
				if(find(query_id_vec.begin(), query_id_vec.end(), reg->query_id) == query_id_vec.end()) query_id_vec.push_back(reg->query_id); // not found
		}
		for(i=0; i<query_id_vec.size(); i++){
			query_id = query_id_vec.at(i);
			group_id = getReadGroupId(query_id, qnames_vec_cluster, var_cand);
			group_id_vec.push_back(group_id);
		}

		if(group_id_vec.size()==1){
			ps_id_vec.push_back(var_cand->newVarVec.at(0)->startRefPos); // ??????
		//}else if(group_id_vec.size()==2){
		}else if(group_id_vec.size()>=2){
			// get minimal value
			for(i=0; i<query_id_vec.size(); i++){
				//group_id = group_id_vec.at(i);
				query_id = query_id_vec.at(i);
				min_val = LONG_MAX;
				for(j=0; j<var_cand->newVarVec.size(); j++){
					reg = var_cand->newVarVec.at(j);
					if(reg->var_type!=VAR_UNC and reg->call_success_status and reg->query_id==query_id){
						if(reg->startRefPos<min_val) min_val = reg->startRefPos;
					}
				}
				ps_id_vec.push_back(min_val);
			}
			if(ps_id_vec.size()>=2 and ps_id_vec.at(0)==ps_id_vec.at(1) and ps_id_vec.at(0)!=LONG_MAX) ps_id_vec.at(1)++;  // ??????????
		}
//		else{
//			cerr << "group count: " << group_id_vec.size() << ", error!" << endl;
//			exit(1);
//		}

		for(i=0; i<var_cand->newVarVec.size(); i++){
			reg = var_cand->newVarVec.at(i);
			if(reg->var_type!=VAR_UNC and reg->call_success_status and reg->query_id!=-1){
				group_id = id_vec = -1;
				for(j=0; j<query_id_vec.size(); j++){
					if(query_id_vec.at(j)==reg->query_id){
						group_id = group_id_vec.at(j);
						id_vec = j;
						break;
					}
				}

				if(group_id!=-1){
					phasing_reg = new phasing_reg_t();
					phasing_reg->reg = reg;
					phasing_reg->inner_group_id = group_id;
					phasing_reg->next_group_num = 0;
					phasing_reg->next_group_id = phasing_reg->next_group_id2 = -1;
					phasing_reg->ps_id = ps_id_vec.at(id_vec);
					phasing_reg->chain_link_flag = false;
					phasing_reg->redundant_flag = false;
					phasing_reg->read_name_vec = qnames_vec_cluster.at(group_id);
					hapnodes_vec.at(group_id).push_back(phasing_reg);
				}
			}
		}
	}

#if PHASING_DEBUG
	printHapNodes(hapnodes_vec);
#endif

}

// link haplotype node
bool Phasing::linkHapNodes(vector< vector<phasing_reg_t*> > &hapnodes_vec, vector< vector<phasing_reg_t*> > &hapnodes_vec2){
	bool success_flag;
	int32_t link_num, id_hapnode_vec, id_hapnode_vec_nolink, num1, num2, hap_num, max_val, sec_val, max_idx, sec_idx, other_idx;
	size_t i, j, k;
	phasing_reg_t *phasing_reg, *phasing_reg2;
	vector<phasing_reg_t*> hapnode_vec, hapnode_vec2;
	vector<string> qname_vec, qname_vec2;
	vector< vector<int32_t> > links_num_vec;
	vector<int32_t> link_nums, stop_linkid_vec, success_nextlinkid_vec;

	// process links
	num1 = num2 = 0;
	for(i=0; i<hapnodes_vec.size(); i++) if(hapnodes_vec.at(i).size()>0) num1 ++;
	for(i=0; i<hapnodes_vec2.size(); i++) if(hapnodes_vec2.at(i).size()>0) num2 ++;

#if PHASING_DEBUG
		cout << __LINE__ << ", num1=" << num1 << ", num2=" << num2 << endl;
#endif

	// different cases
	success_flag = false;
	if(num1==1 and num2==1){ // 1--1
		hapnode_vec = hapnodes_vec.at(0).size()>0 ? hapnodes_vec.at(0) : hapnodes_vec.at(1);
		hapnode_vec2 = hapnodes_vec2.at(0).size()>0 ? hapnodes_vec2.at(0) : hapnodes_vec2.at(1);
		qname_vec = hapnode_vec.at(hapnode_vec.size()-1)->read_name_vec;
		qname_vec2 = hapnode_vec2.at(0)->read_name_vec;
		link_num = computeLinkNum(qname_vec, qname_vec2);

#if PHASING_DEBUG
		cout << __LINE__ << ", link_num=" << link_num << endl;
#endif

		if(link_num>=MIN_LINK_READS_NUM) {
			phasing_reg = hapnode_vec.at(hapnode_vec.size()-1);
			phasing_reg->next_group_id = 0;
			phasing_reg->next_group_num = 1;
			phasing_reg->chain_link_flag = true;
			for(j=0; j<hapnode_vec2.size(); j++) {
				phasing_reg2 = dupPhasingRegItem(hapnode_vec2.at(j));
				phasing_reg2->ps_id = phasing_reg->ps_id;
				hapnodes_vec.at(0).push_back(phasing_reg2);
			}
			success_flag = true;
		}else{ // link stop
			phasing_reg = hapnode_vec.at(hapnode_vec.size()-1);
			phasing_reg->next_group_id = -1;
			phasing_reg->next_group_num = 0;
			phasing_reg->chain_link_flag = false;
//			for(j=0; j<hapnode_vec2.size(); j++) {
//				phasing_reg2 = dupPhasingRegItem(hapnode_vec2.at(j));
//				hapnodes_vec.at(0).push_back(phasing_reg2);
//			}
		}
	}else{
		// link nodes
		id_hapnode_vec = -1;
		for(i=0; i<hapnodes_vec.size(); i++){
			hapnode_vec = hapnodes_vec.at(i);
			link_nums.clear();
			if(hapnode_vec.size()>0){
				id_hapnode_vec = i;
				qname_vec = hapnode_vec.at(hapnode_vec.size()-1)->read_name_vec;
				for(j=0; j<hapnodes_vec2.size(); j++){
					hapnode_vec2 = hapnodes_vec2.at(j);
					if(hapnode_vec2.size()>0){
						qname_vec2 = hapnode_vec2.at(0)->read_name_vec;
						link_num = computeLinkNum(qname_vec, qname_vec2);

#if PHASING_DEBUG
						cout << __LINE__ << ", link_num=" << link_num << endl;
#endif

						link_nums.push_back(link_num);
					}
				}
			}
			links_num_vec.push_back(link_nums);
		}

		if(num1==1 and num2==2){ // 1--2
			// sort
			hapnode_vec = hapnodes_vec.at(id_hapnode_vec);
			if(hapnode_vec.size()>0){
				link_nums = links_num_vec.at(id_hapnode_vec);
				hap_num = 0;
				max_val = sec_val = max_idx = sec_idx = -1;
				for(j=0; j<link_nums.size(); j++){
					if(link_nums.at(j)>0){
						hap_num ++;
						if(max_val<link_nums.at(j)){
							sec_val = max_val;
							sec_idx = max_idx;
							max_val = link_nums.at(j);
							max_idx = j;
						}else if(sec_val<link_nums.at(j)){
							sec_val = link_nums.at(j);
							sec_idx = j;
						}
					}
				}
			}

			link_nums = links_num_vec.at(id_hapnode_vec);
			if(max_val>=MIN_LINK_READS_NUM){
				hapnode_vec = hapnodes_vec.at(id_hapnode_vec);
				hapnode_vec2 = hapnodes_vec2.at(max_idx);
				phasing_reg = hapnode_vec.at(hapnode_vec.size()-1);
				phasing_reg->next_group_id = max_idx;
				phasing_reg->next_group_id2 = sec_idx;
				phasing_reg->next_group_num = hap_num;
				phasing_reg->chain_link_flag = true;
				for(j=0; j<hapnode_vec2.size(); j++){
					phasing_reg2 = dupPhasingRegItem(hapnode_vec2.at(j));
					phasing_reg2->ps_id = phasing_reg->ps_id;
					hapnodes_vec.at(id_hapnode_vec).push_back(phasing_reg2);
				}
				success_flag = true;
			}else{
				hapnode_vec = hapnodes_vec.at(id_hapnode_vec);
				phasing_reg = hapnode_vec.at(hapnode_vec.size()-1);
				phasing_reg->next_group_id = -1;
				phasing_reg->next_group_id2 = -1;
				phasing_reg->next_group_num = 0;
				phasing_reg->chain_link_flag = false;
			}

			if(success_flag){
				// link stop on other hap strand
				id_hapnode_vec_nolink = id_hapnode_vec==0 ? 1 : 0;
				other_idx = max_idx==0 ? 1 : 0;
				hapnode_vec = hapnodes_vec.at(id_hapnode_vec_nolink);
				hapnode_vec2 = hapnodes_vec2.at(other_idx);
				if(hapnode_vec.size()>0){ // update
					phasing_reg = hapnode_vec.at(hapnode_vec.size()-1);
					phasing_reg->next_group_id = -1;
					phasing_reg->next_group_id2 = -1;
					phasing_reg->next_group_num = 0;
					phasing_reg->chain_link_flag = false;
					for(j=0; j<hapnode_vec2.size(); j++){
						phasing_reg2 = dupPhasingRegItem(hapnode_vec2.at(j));
						hapnodes_vec.at(id_hapnode_vec_nolink).push_back(phasing_reg2);
					}
				}else{ // empty list, then add new item
					for(j=0; j<hapnode_vec2.size(); j++){
						phasing_reg2 = dupPhasingRegItem(hapnode_vec2.at(j));

						hapnodes_vec.at(id_hapnode_vec_nolink).push_back(phasing_reg2);
					}
				}
			}
		}else if(num1==2){
			// sort
			stop_linkid_vec.clear();
			success_nextlinkid_vec.clear();
			for(i=0; i<hapnodes_vec.size(); i++){
				hapnode_vec = hapnodes_vec.at(i);
				if(hapnode_vec.size()>0){
					link_nums = links_num_vec.at(i);
					max_val = sec_val = max_idx = sec_idx = -1;
					for(j=0; j<link_nums.size(); j++){
						if(link_nums.at(j)>0){
							if(max_val<link_nums.at(j)){
								sec_val = max_val;
								sec_idx = max_idx;
								max_val = link_nums.at(j);
								max_idx = j;
							}else if(sec_val<link_nums.at(j)){
								sec_val = link_nums.at(j);
								sec_idx = j;
							}
						}
					}

					if(max_val>=MIN_LINK_READS_NUM){
//					if(num2==1){ // 2--1
						phasing_reg = hapnode_vec.at(hapnode_vec.size()-1);
						phasing_reg->next_group_id = max_idx;
						phasing_reg->next_group_id2 = -1;
						phasing_reg->next_group_num = 1;
						phasing_reg->chain_link_flag = true;
						hapnode_vec2 = hapnodes_vec2.at(max_idx);
						for(j=0; j<hapnode_vec2.size(); j++){
							phasing_reg2 = dupPhasingRegItem(hapnode_vec2.at(j));
							phasing_reg2->ps_id = phasing_reg->ps_id;
							hapnodes_vec.at(i).push_back(phasing_reg2);
						}
						success_flag = true;
						success_nextlinkid_vec.push_back(max_idx);
//					}else if(num2==2){ // 2--2
//
//					}
					}else{ // link stop
						stop_linkid_vec.push_back(i);
					}
				}
			}

			// next round, record stop links
			for(i=0; i<stop_linkid_vec.size(); i++){
				hapnode_vec = hapnodes_vec.at(stop_linkid_vec.at(i));
				if(hapnode_vec.size()>0){
					phasing_reg = hapnode_vec.at(hapnode_vec.size()-1);
					phasing_reg->next_group_id = -1;
					phasing_reg->next_group_id2 = -1;
					phasing_reg->next_group_num = 0;
					phasing_reg->chain_link_flag = false;
					if(success_flag){
						for(j=0; j<hapnodes_vec2.size(); j++){
							if(find(success_nextlinkid_vec.begin(), success_nextlinkid_vec.end(), j) == success_nextlinkid_vec.end()){  // ignore success items
								hapnode_vec2 = hapnodes_vec2.at(j);
								for(k=0; k<hapnode_vec2.size(); k++){
									phasing_reg2 = dupPhasingRegItem(hapnode_vec2.at(k));
									//phasing_reg2->ps_id = phasing_reg->ps_id;
									hapnodes_vec.at(stop_linkid_vec.at(i)).push_back(phasing_reg2);
								}
							}
						}
					}
				}
			}
		}else{ // do nothing

		}
	}

	return success_flag;
}

// compute number of links
int32_t Phasing::computeLinkNum(vector<string> &qname_vec, vector<string> &qname_vec2){
	int32_t link_num;
	size_t i;
	string qname;

	link_num = 0;
	for(i=0; i<qname_vec.size(); i++){
		qname = qname_vec.at(i);
		if(find(qname_vec2.begin(), qname_vec2.end(), qname) != qname_vec2.end()) link_num ++; // found
	}

	return link_num;
}

// get group id
int32_t Phasing::getReadGroupId(int32_t query_id, vector< vector<string> > &qnames_vec_cluster, varCand *var_cand){
	int32_t group_id, check_num;
	vector< vector<string> > qnames_vec;
	vector<string> qname_vec, qname_vec_cluster;
	string qname, tmp_cnsfilename;
	size_t i, j;

	if(var_cand->rescue_flag) tmp_cnsfilename = var_cand->rescue_cnsfilename;
	else tmp_cnsfilename = var_cand->ctgfilename;

	FastaSeqLoader fa_loader(tmp_cnsfilename);
	qnames_vec = fa_loader.getReadNamesFromCnsHeaders();
//	if(query_id>=(int32_t)qnames_vec_cluster.size()){
//		cout << "query_id=" << query_id << ", qnames_vec_cluster.size=" << qnames_vec_cluster.size() << endl;
//	}
	qname_vec = qnames_vec.at(query_id);

	group_id = -1;
	for(i=0; i<qnames_vec_cluster.size(); i++){
		qname_vec_cluster = qnames_vec_cluster.at(i);
		check_num = 0;
		for(j=0; j<qname_vec.size(); j++){
			qname = qname_vec.at(j);
			if(find(qname_vec_cluster.begin(), qname_vec_cluster.end(), qname) != qname_vec_cluster.end()){  // found
				check_num ++;
				if(check_num>=MAX_CHECK_READS_NUM){
					group_id = i;
					break;
				}
			}
		}

#if PHASING_DEBUG
		cout << __LINE__ << ", In " << __func__ << "(), check_num=" << check_num << endl;
#endif

		if(group_id!=-1) break;
	}

	return group_id;
}

phasing_reg_t* Phasing::dupPhasingRegItem(phasing_reg_t *phasing_reg){
	phasing_reg_t *phasing_reg_new;

	phasing_reg_new = new phasing_reg_t();
	phasing_reg_new->reg = phasing_reg->reg;
	phasing_reg_new->inner_group_id = phasing_reg->inner_group_id;
	phasing_reg_new->next_group_num = phasing_reg->next_group_num;
	phasing_reg_new->next_group_id = phasing_reg->next_group_id;
	phasing_reg_new->next_group_id2 = phasing_reg->next_group_id2;
	phasing_reg_new->ps_id = phasing_reg->ps_id;
	phasing_reg_new->chain_link_flag = phasing_reg->chain_link_flag;
	phasing_reg_new->redundant_flag = phasing_reg->redundant_flag;
	phasing_reg_new->read_name_vec = phasing_reg->read_name_vec;

	return phasing_reg_new;
}

// check whether the start link hap node is valid
bool Phasing::isValidStartLinkHapNode(vector< vector<phasing_reg_t*> > &hapnodes_vec){
	bool flag = false;
	for(size_t i=0; i<hapnodes_vec.size(); i++){
		vector<phasing_reg_t*> &hapnode_vec = hapnodes_vec.at(i);
		if(hapnode_vec.size()>0){
			flag = true;
			break;
		}
	}
	return flag;
}

// print hap nodes
void Phasing::printHapNodes(vector< vector<phasing_reg_t*> > &hapnodes_vec){
	size_t i, j;
	phasing_reg_t *phasing_reg;
	vector<phasing_reg_t*> hapnode_vec;
	int32_t hap_id;

	hap_id = 0;
	for(i=0; i<hapnodes_vec.size(); i++){
		hapnode_vec = hapnodes_vec.at(i);
		if(hapnode_vec.size()>0){
			cout << "Hap " << hap_id << ", size=" << hapnode_vec.size() << " :";
			hap_id ++;
			for(j=0; j<hapnode_vec.size(); j++){
				phasing_reg = hapnode_vec.at(j);
				cout << "\t(" << phasing_reg->ps_id << "," << phasing_reg->reg->chrname << ":" << phasing_reg->reg->startRefPos << "-" << phasing_reg->reg->endRefPos
					<< ",type=" << to_string(phasing_reg->reg->var_type) << ",len=" << phasing_reg->reg->sv_len << ",PS=" << phasing_reg->reg->PS << ",gt_seq=" << phasing_reg->reg->gt_seq << "), chain_link=" << phasing_reg->chain_link_flag;
			}
			cout << endl;
		}
	}
}
