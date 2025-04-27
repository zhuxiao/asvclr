#include "localCns.h"

#include "clipAlnDataLoader.h"

pthread_mutex_t mutex_write = PTHREAD_MUTEX_INITIALIZER;
//extern pthread_mutex_t mutex_down_sample;
extern pthread_mutex_t mutex_fai;

localCns::localCns(string &readsfilename, string &contigfilename, string &refseqfilename, string &clusterfilename, string &tmpdir, string &technology, double min_identity_match, int32_t sv_len_est, size_t num_threads_per_cns_work, vector<reg_t*> &varVec, string &chrname, string &inBamFile, faidx_t *fai, size_t cns_extend_size, double expected_cov, double min_input_cov, double max_ultra_high_cov, int32_t minMapQ, int32_t minHighMapQ, bool delete_reads_flag, bool keep_failed_reads_flag, bool clip_reg_flag, int32_t minClipEndSize, int32_t minConReadLen, int32_t min_sv_size, int32_t min_supp_num, double max_seg_size_ratio, double max_seg_nm_ratio){

	this->chrname = chrname;
	this->chrlen = faidx_seq_len(fai, chrname.c_str()); // get reference size
	this->readsfilename = preprocessPipeChar(readsfilename);
	this->contigfilename = preprocessPipeChar(contigfilename);
	this->refseqfilename = preprocessPipeChar(refseqfilename);
	this->clusterfilename = preprocessPipeChar(clusterfilename);
	this->tmpdir = preprocessPipeChar(tmpdir);
	this->technology = technology;
	this->sv_len_est = sv_len_est;
	this->num_threads_per_cns_work = num_threads_per_cns_work;
	this->minClipEndSize = minClipEndSize;
	this->minMapQ = minMapQ;
	this->minHighMapQ = minHighMapQ;
	this->minConReadLen = minConReadLen;
	this->min_sv_size = min_sv_size;
	this->min_supp_num = min_supp_num;
	this->max_seg_size_ratio = max_seg_size_ratio;
	this->max_seg_nm_ratio = max_seg_nm_ratio;
	this->max_ultra_high_cov = max_ultra_high_cov;
	this->min_identity_match = min_identity_match;
	this->varVec = varVec;
	this->fai = fai;
	this->inBamFile = inBamFile;
	//this->cns_extend_size = CNS_SIDE_EXT_SIZE + cns_extend_size;
	this->cns_extend_size = cns_extend_size;
	startRefPos_cns = endRefPos_cns = 0;
	mean_read_len = 0;

	if(sv_len_est<=MAX_SV_LEN_USING_POA and clip_reg_flag==false) use_poa_flag = true;
//	if(sv_len_est<=MAX_SV_LEN_USING_POA and varVec.at(varVec.size()-1)->endRefPos-varVec.at(0)->startRefPos<=MIN_SEQ_LEN_USING_MINIMIZER and clip_reg_flag==false) use_poa_flag = true;
	else use_poa_flag = false;
	//cout << "sv_len_est=" << sv_len_est << ", use_poa_flag=" << use_poa_flag << endl;

	readsfilename_prefix = this->readsfilename.substr(0, this->readsfilename.size()-3);
	readsfilename_suffix = this->readsfilename.substr(this->readsfilename.size()-3, 3);

	ref_seq_size = reads_count_original = total_bases_original = reads_count = total_bases = 0;
	local_cov_original = sampled_cov = 0;
	this->expected_cov = expected_cov;
	this->compensation_coefficient = 1;
	sampling_flag = false;
	this->delete_reads_flag = delete_reads_flag;
	this->keep_failed_reads_flag = keep_failed_reads_flag;
	cns_success_flag = false;
	this->clip_reg_flag = clip_reg_flag;

	if(isFileExist(contigfilename) and isFileExist(refseqfilename)) cns_success_preDone_flag = true;
	else cns_success_preDone_flag = false;

	limit_reg_process_flag = false;
	this->min_input_cov_canu = min_input_cov;

	time(&start_time);
	end_time = 0;
}

localCns::~localCns() {
	if(delete_reads_flag){
		if(keep_failed_reads_flag==false or (keep_failed_reads_flag and cns_success_flag)){
//			if(clip_reg_flag==false){ // indel region
				for(size_t i=0; i<readsfilename_vec.size(); i++) remove(readsfilename_vec.at(i).c_str());
//			}else{
//				remove(readsfilename.c_str());	// delete the reads file to save disk space
//			}
		}
	}
}

// destroy the alignment data of the block
void localCns::destoryClipAlnData(){
	for(size_t i=0; i<clipAlnDataVector.size(); i++){
		bam_destroy1(clipAlnDataVector.at(i)->bam);
		delete clipAlnDataVector.at(i);
	}
	vector<clipAlnData_t*>().swap(clipAlnDataVector);
}

void localCns::destorySeqsVec(vector<struct seqsVec*> &seqs_vec){
	struct seqsVec* seqs_vec_node;
	for(size_t i=0; i<seqs_vec.size(); i++){
		seqs_vec_node = seqs_vec.at(i);
		delete seqs_vec_node;
	}
	vector<struct seqsVec*>().swap(seqs_vec);
}
void localCns::destoryQueryCluVec(vector<struct querySeqInfoVec*> &query_clu_vec){
	struct querySeqInfoVec *qc_vec_node;
//	struct querySeqInfoNode *q_info_node;

	for(size_t i=0; i<query_clu_vec.size(); i++){
		qc_vec_node = query_clu_vec.at(i);
//		for(size_t j; j<query_clu_vec.at(i)->query_seq_info_all.size();j++){
//			q_info_node = query_clu_vec.at(i)->query_seq_info_all.at(j);
//			delete q_info_node;
//		}
//		vector<struct querySeqInfoNode*>().swap(query_clu_vec.at(i)->query_seq_info_all);
		delete qc_vec_node;
	}
	vector<struct querySeqInfoVec*>().swap(query_clu_vec);
}
void localCns::destoryQueryCluVecClipReg(vector<qcSigListVec_t*> &query_clu_vec_clipReg){
	qcSigListVec_t *qc_vec_node;
	size_t i, j, k;
	qcSigList_t *qcSig_list;

	for(i=0; i<query_clu_vec_clipReg.size(); i++){
		qc_vec_node = query_clu_vec_clipReg.at(i);
		for(j=0; j<qc_vec_node->qcSigList.size(); j++){
			qcSig_list = qc_vec_node->qcSigList.at(j);
			for(k=0; k<qcSig_list->qcSig_vec.size(); k++)
				delete qcSig_list->qcSig_vec.at(k);
			vector<qcSig_t*>().swap(qcSig_list->qcSig_vec);
			delete qcSig_list;
		}
		vector<qcSigList_t*>().swap(qc_vec_node->qcSigList);
		delete qc_vec_node;
	}
	vector<qcSigListVec_t*>().swap(query_clu_vec_clipReg);
}

// set limit process regions
void localCns::setLimitRegs(bool limit_reg_process_flag, vector<simpleReg_t*> limit_reg_vec){
	this->limit_reg_process_flag = limit_reg_process_flag;
	for(size_t i=0; i<limit_reg_vec.size(); i++) this->limit_reg_vec.push_back(limit_reg_vec.at(i));
}

// extract the corresponding refseq from reference
void localCns::extractRefseq(){
	int32_t startRefPos, endRefPos, left_shift_size, right_shift_size;
	string reg, header;
	ofstream outfile;

	// generate the refseq region
	startRefPos = varVec[0]->startRefPos - cns_extend_size;
	if(startRefPos<1) startRefPos = 1;
	left_shift_size = varVec[0]->startRefPos - startRefPos;  // left shift size
	endRefPos = varVec[varVec.size()-1]->endRefPos + cns_extend_size;
	if(endRefPos>chrlen) endRefPos = chrlen;
	right_shift_size = endRefPos - varVec[varVec.size()-1]->endRefPos;  // right shift size

	// get the region sequence
	reg = chrname + ":" + to_string(startRefPos) + "-" + to_string(endRefPos);
	RefSeqLoader refseq_loader(reg, fai);
	refseq_loader.getRefSeq();
	ref_seq_size = refseq_loader.refseq_len;

	// save the refseq to file
	outfile.open(refseqfilename);
	if(!outfile.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file " << refseqfilename << endl;
		exit(1);
	}

	outfile << ">" + reg + "___" + to_string(left_shift_size) + "___" + to_string(right_shift_size) + "___" + contigfilename << endl;  // header
	outfile << refseq_loader.refseq << endl;  // seq

	outfile.close();
}

// extract the reads data from BAM file
void localCns::extractReadsDataFromBAM(){
	size_t i;
	string qname, seq, qual, reg_str, refseq;
	vector<struct querySeqInfoVec*> query_clu_vec;
	vector<qcSigListVec_t*> query_clu_vec_clipReg;
	vector<struct querySeqInfoNode*> query_seq_info_all;
	int32_t seq_len, maxValue;
	char *p_seq;
	struct seqsVec *smoothed_seqs;
	struct seqsVec *seqs_vec_node;
	struct querySeqInfoVec *q_cluster_node;
	vector<struct alnSeg*> alnSegs;
	double expected_cov_cons;
	qcSigListVec_t *qc_vec_node;
	bool insufficient_flag;

	startRefPos_cns = varVec[0]->startRefPos - cns_extend_size;
	if(startRefPos_cns<1) startRefPos_cns = 1;
	endRefPos_cns = varVec[varVec.size()-1]->endRefPos + cns_extend_size;
	if(endRefPos_cns>chrlen) endRefPos_cns = chrlen;

	// load the clipping data ---20220705
	clipAlnDataLoader data_loader(varVec[0]->chrname, startRefPos_cns, endRefPos_cns, inBamFile, minClipEndSize, minMapQ, minHighMapQ);
	if(clip_reg_flag) data_loader.loadClipAlnDataWithSATag(clipAlnDataVector, max_ultra_high_cov);
	else data_loader.loadClipAlnDataWithSATagWithSegSize(clipAlnDataVector, max_ultra_high_cov, max_seg_size_ratio, max_seg_nm_ratio);

	// update the use_poa_flag
//	if(use_poa_flag) {
//		use_poa_flag = updateUsepoaFlag(clipAlnDataVector, MAX_VAR_REG_SIZE, 0.8);
////		if(use_poa_flag==false){
////			cout << "line=" << __LINE__ << ", clipAlnDataVector.size=" << clipAlnDataVector.size() << ", use_poa_flag=" << use_poa_flag << endl;
////		}
//	}

	if(clipAlnDataVector.size()>0){

		for(i=0; i<clipAlnDataVector.size(); i++) clipAlnDataVector.at(i)->query_checked_flag = false;

		reg_str = varVec[0]->chrname + ":" + to_string(startRefPos_cns) + "-" + to_string(endRefPos_cns);
		pthread_mutex_lock(&mutex_fai);
		p_seq = fai_fetch(fai, reg_str.c_str(), &seq_len);
		pthread_mutex_unlock(&mutex_fai);
		refseq = p_seq;
		free(p_seq);

		// extract queries from clip align data vector
		query_seq_info_all = extractQueriesFromClipAlnDataVec(clipAlnDataVector, refseq, chrname, startRefPos_cns, endRefPos_cns, fai, minConReadLen, clip_reg_flag, &mutex_fai);

		// sampling for ultra-high coverage regions
		expected_cov_cons = expected_cov;
		reads_count_original = query_seq_info_all.size();
		total_bases_original = 0;
		for(i=0; i<query_seq_info_all.size(); i++) total_bases_original += query_seq_info_all.at(i)->seq.size();
		mean_read_len = (double)total_bases_original / query_seq_info_all.size();

		if(expected_cov_cons!=0){
			compensation_coefficient = computeCompensationCoefficient(startRefPos_cns, endRefPos_cns, mean_read_len);
			//cout << compensation_coefficient << endl;
			samplingReads(query_seq_info_all, expected_cov_cons, compensation_coefficient);
		}

		if(clip_reg_flag){ // clipping region
			if(query_seq_info_all.size()>0){
				insufficient_flag = false;
				clipRegCluster clip_reg_cluster(varVec[0]->chrname, varVec[0]->startRefPos, varVec[varVec.size()-1]->endRefPos, minClipEndSize, min_sv_size, min_supp_num, min_identity_match, technology, fai);
				//if(query_seq_info_all.size()>1){
				if(query_seq_info_all.size()>=(size_t)min_supp_num){ // 2024-01-28
					query_clu_vec_clipReg = clip_reg_cluster.queryCluster(query_seq_info_all);

					maxValue = 0;
					for(i=0; i<query_clu_vec_clipReg.size(); i++)
						if(maxValue<(int32_t)query_clu_vec_clipReg.at(i)->qcSigList.size()) maxValue = query_clu_vec_clipReg.at(i)->qcSigList.size();
					if(maxValue<min_supp_num) { // remove insufficient items
						insufficient_flag = true;
						destoryQueryCluVecClipReg(query_clu_vec_clipReg);
					}
				}else
					insufficient_flag = true;

				if(insufficient_flag){
					qc_vec_node = new qcSigListVec_t();
					qc_vec_node->qcSigList = clip_reg_cluster.extractQcSigsClipReg(query_seq_info_all);
					query_clu_vec_clipReg.push_back(qc_vec_node);
				}
				//cout << "~~~~~~~~~~ query_clu_vec_clipReg.size=" << query_clu_vec_clipReg.size() << endl;
				//save read segments
				for(i=0; i<query_clu_vec_clipReg.size(); i++){
					seqs_vec_node = collectQuerySeqDataClipReg(query_clu_vec_clipReg.at(i)->qcSigList);
					seqs_vec.push_back(seqs_vec_node);
				}
				destoryQueryCluVecClipReg(query_clu_vec_clipReg);
			}
		}else{ // indel region
			if(query_seq_info_all.size()>0){
				insufficient_flag = false;
				if(query_seq_info_all.size()>(size_t)min_supp_num){

					indelRegCluster indel_reg_cluster(varVec[0]->chrname, varVec[0]->startRefPos, varVec[varVec.size()-1]->endRefPos, sv_len_est, startRefPos_cns, endRefPos_cns, min_sv_size, min_supp_num, min_identity_match, fai);

					query_clu_vec = indel_reg_cluster.queryCluster(query_seq_info_all);

					maxValue = 0;
					for(i=0; i<query_clu_vec.size(); i++){
						q_cluster_node = query_clu_vec.at(i);
						if(q_cluster_node->rescue_flag==false and maxValue<(int32_t)q_cluster_node->query_seq_info_all.size()) maxValue = q_cluster_node->query_seq_info_all.size();
					}
					if(maxValue<min_supp_num) { // remove insufficient items
						insufficient_flag = true;
						destoryQueryCluVec(query_clu_vec);
					}
				}else
					insufficient_flag = true;

				if(insufficient_flag){
					q_cluster_node = new struct querySeqInfoVec();
					q_cluster_node->query_seq_info_all = query_seq_info_all;
					query_clu_vec.push_back(q_cluster_node);
				}
				//smooth read segments
				for(i=0; i<query_clu_vec.size(); i++){
					q_cluster_node = query_clu_vec.at(i);
					if(q_cluster_node->rescue_flag==false){
						smoothed_seqs = smoothQuerySeqData(refseq, query_clu_vec.at(i)->query_seq_info_all, startRefPos_cns, min_sv_size);
						seqs_vec.push_back(smoothed_seqs);
					}
				}
			}
		}

		// save cluster information
		if(seqs_vec.size()>0) saveClusterInfo(clusterfilename, seqs_vec);
	}
	if(!clipAlnDataVector.empty()) destoryClipAlnData();

	// release memory
	if(!query_seq_info_all.empty()) destoryQuerySeqInfoAll(query_seq_info_all);
	if(!query_clu_vec.empty()) destoryQueryCluVec(query_clu_vec);

}

void localCns::saveClusterInfo(string &clusterfilename, vector<struct seqsVec*> &seqs_vec){
	ofstream outfile;
	size_t i, j;
	struct seqsVec *seq_vec;
	string line, qname_str;

	outfile.open(clusterfilename);
	if(outfile.is_open()==false){
		cerr << "line=" << __LINE__ << ", In " << __func__ << "(), cannot open file '" << clusterfilename << "', error!" << endl;
		exit(1);
	}

	outfile << "#clusterID\tqueryNum\treadsInfo" << endl;
	for(i=0; i<seqs_vec.size(); i++){
		seq_vec = seqs_vec.at(i);
		if(seq_vec->qname.size()>=(size_t)min_supp_num*READS_NUM_SUPPORT_FACTOR){ // only use sufficient reads data to generate consensus sequence
			line = to_string(i+1) + "\t" + to_string(seq_vec->qname.size());
			qname_str = seq_vec->qname.at(0);
			for(j=1; j<seq_vec->qname.size(); j++) qname_str += ";" + seq_vec->qname.at(j);
			line += "\t" + qname_str;
			outfile << line << endl;
		}
	}
	outfile.close();
}

bool localCns::updateUsepoaFlag(vector<clipAlnData_t*> &clipAlnDataVector, int32_t maxVarRegSize, double repeat_reads_thres){
	bool flag = true;
	size_t i, j, total_reads_num;
	string query_name_tmp;
	double repeat_ratio;
	int32_t repeat_reads_num = 0;
	vector<clipAlnData_t*> query_aln_segs;

	for(i=0; i<clipAlnDataVector.size(); i++) clipAlnDataVector.at(i)->query_checked_flag = false;

	total_reads_num = clipAlnDataVector.size();
	for(i=0; i<total_reads_num; i++){
		if(clipAlnDataVector.at(i)->query_checked_flag==false){
			query_name_tmp = clipAlnDataVector.at(i)->queryname;
			query_aln_segs = getQueryClipAlnSegsAll(query_name_tmp, clipAlnDataVector);
			if(query_aln_segs.size()>=2) {
				if(isQuerySelfOverlap(query_aln_segs, MAX_VAR_REG_SIZE)){
					repeat_reads_num += query_aln_segs.size();
				}
			}
			for(j=0; j<query_aln_segs.size(); j++) query_aln_segs.at(j)->query_checked_flag = true;
		}
	}
	for(i=0; i<clipAlnDataVector.size(); i++) clipAlnDataVector.at(i)->query_checked_flag = false;

	if(repeat_reads_num>0 and total_reads_num>0){
		repeat_ratio = (double)repeat_reads_num/total_reads_num;
		//cout << "repeat_reads_num=" << repeat_reads_num << ", repeat_ratio=" << repeat_ratio << ", clipAlnDataVector.size=" << clipAlnDataVector.size() << endl;
	}else
		repeat_ratio = 0;
	if(repeat_ratio>repeat_reads_thres) {
		flag = false;
		cout << contigfilename << ", repeat_reads_num=" << repeat_reads_num << ", repeat_ratio=" << repeat_ratio << ", clipAlnDataVector.size=" << clipAlnDataVector.size() << endl;
	}
	return flag;
}

int32_t localCns::calculateSignalDistance(alnSeg *query_alnSeg1, alnSeg *query_alnSeg2){
	int32_t alnSeg1_end, alnSeg2_start, alnSegs_distance;

	alnSeg1_end = getEndRefPosAlnSeg(query_alnSeg1->startRpos, query_alnSeg1->opflag, query_alnSeg1->seglen);
	alnSeg2_start = query_alnSeg2->startRpos;
	alnSegs_distance = abs(alnSeg1_end - alnSeg2_start);
    return alnSegs_distance;
}


// collect query sequence data
struct seqsVec* localCns::collectQuerySeqDataClipReg(vector<qcSigList_t*> &qcSigList){
	struct seqsVec *seqs_info;
	size_t i, j;
	qcSigList_t *q_node;
	string seq, qname;
	int32_t plus_orient_id;

	seqs_info = new struct seqsVec();

	for(i=0; i<qcSigList.size(); i++){
		q_node = qcSigList.at(i);
		plus_orient_id = -1;
		for(j=0; j<q_node->query_seqs_vec.size(); j++){
			if(q_node->query_seqs_vec.at(j)->clip_aln->aln_orient==ALN_PLUS_ORIENT){
				plus_orient_id = j;
				break;
			}
		}
		if(plus_orient_id==-1) plus_orient_id = 0;
		seq = q_node->query_seqs_vec.at(plus_orient_id)->seq;
		qname = q_node->query_seqs_vec.at(plus_orient_id)->qname;

		seqs_info->qname.push_back(qname);
		seqs_info->seqs.push_back(seq);
	}

	return seqs_info;
}

double localCns::computeCompensationCoefficient(size_t startRefPos_cns, size_t endRefPos_cns, double mean_read_len){
	size_t total_reg_size, refined_reg_size, reg_size_cns;
	double comp_coefficient = 1.0;

	if(clip_reg_flag){
		reg_size_cns = endRefPos_cns - startRefPos_cns + 1;
		total_reg_size = reg_size_cns + 2 * mean_read_len;
		refined_reg_size = reg_size_cns + mean_read_len;  // flanking_area = (2 * mean_read_len) / 2
		comp_coefficient = (double)total_reg_size / refined_reg_size;
	}
	return comp_coefficient;
}

double localCns::computeLocalCovClipReg(vector<struct fqSeqNode*> &fq_seq_full_vec, double compensation_coefficient){
	double cov = 0, total;
	struct fqSeqNode* fq_node;
	size_t refined_reg_size;

	refined_reg_size = endRefPos_cns - startRefPos_cns + 1 + 2 * mean_read_len;
	if(refined_reg_size>0){
		total = 0;
		for(size_t i=0; i<fq_seq_full_vec.size(); i++){
			fq_node = fq_seq_full_vec.at(i);
			total += fq_node->seq.size();
		}
		cov = total / refined_reg_size * compensation_coefficient;
		//cout << "total bases: " << total << " bp, local coverage: " << cov << endl;
	}else{
		cov = -1;
		//cerr << "ERR: ref_seq_size=" << ref_seq_size << endl;
	}
	return cov;
}

double localCns::computeLocalCovIndelReg(vector<struct querySeqInfoNode*> &query_seq_info_all, double compensation_coefficient){
	double cov = 0, total;
	struct querySeqInfoNode* qseq_node;
	size_t refined_reg_size;

	refined_reg_size = endRefPos_cns - startRefPos_cns + 1 + 2 * mean_read_len;
	if(refined_reg_size>0){
		total = 0;
		for(size_t i=0; i<query_seq_info_all.size(); i++){
			qseq_node = query_seq_info_all.at(i);
			total += qseq_node->seq.size();
		}
		cov = total / refined_reg_size * compensation_coefficient;
		//cout << "total bases: " << total << " bp, local coverage: " << cov << endl;
	}else{
		cov = -1;
		//cerr << "ERR: ref_seq_size=" << ref_seq_size << endl;
	}
	return cov;
}

void localCns::samplingReadsClipReg(vector<struct fqSeqNode*> &fq_seq_vec, double expect_cov_val, double compensation_coefficient){
	local_cov_original = computeLocalCovClipReg(fq_seq_vec, compensation_coefficient);

	if(local_cov_original>expect_cov_val){ // sampling
		//cout << "sampling for " << readsfilename << ", original coverage: " << local_cov_original << ", expected coverage: " << expect_cov_val << ", compensation_coefficient: " << compensation_coefficient << endl;
		samplingReadsClipRegOp(fq_seq_vec, expect_cov_val);
	}
}

void localCns::samplingReadsClipRegOp(vector<struct fqSeqNode*> &fq_seq_vec, double expect_cov_val){
	double expected_total_bases, total_bases;
	size_t index, max_reads_num, reg_size;
	struct fqSeqNode* fq_node;

	// reverse the select flag
	for(size_t i=0; i<fq_seq_vec.size(); i++){
		fq_node = fq_seq_vec.at(i);
		fq_node->selected_flag = false;
	}

	reg_size = endRefPos_cns - startRefPos_cns + 1 + mean_read_len;
	expected_total_bases = reg_size * expect_cov_val;
	max_reads_num = fq_seq_vec.size();

	//pthread_mutex_lock(&mutex_down_sample);
	srand(1);
	reads_count = 0;
	total_bases = 0;
	while(total_bases<=expected_total_bases and reads_count<=max_reads_num){
		index = rand() % max_reads_num;
		fq_node = fq_seq_vec.at(index);
		if(fq_node->selected_flag==false){
			fq_node->selected_flag = true;

			total_bases += fq_node->seq.size();
			reads_count ++;
		}
	}
	//pthread_mutex_unlock(&mutex_down_sample);
	sampled_cov = total_bases / reg_size;
	sampling_flag = true;

//	/cout << "After sampling, reads count: " << reads_count << ", total bases: " << total_bases << endl;
}

void localCns::samplingReads(vector<struct querySeqInfoNode*> &query_seq_info_all, double expect_cov_val, double compensation_coefficient){
	local_cov_original = computeLocalCovIndelReg(query_seq_info_all, compensation_coefficient);

	if(local_cov_original>expect_cov_val){ // sampling
		//cout << "sampling for " << readsfilename << ", original coverage: " << local_cov_original << ", expected coverage: " << expect_cov_val << ", compensation_coefficient: " << compensation_coefficient << endl;
		samplingReadsOp(query_seq_info_all, expect_cov_val);
	}
}

void localCns::samplingReadsOp(vector<struct querySeqInfoNode*> &query_seq_info_all, double expect_cov_val){
	double expected_total_bases, total_bases;
	size_t i, index, num, count, max_reads_num, reg_size;
	struct querySeqInfoNode* qseq_node;
	set<string> selected_qname_vec;
	//int32_t k, min_id, max_id;
	//clipAlnData_t *clip_aln;
	//vector<int32_t> high_qual_id_vec;
	//string qname;

	// reset the select flag
	for(i=0; i<query_seq_info_all.size(); i++) query_seq_info_all.at(i)->selected_flag = false;

	reg_size = endRefPos_cns - startRefPos_cns + 1 + mean_read_len;
	expected_total_bases = reg_size * expect_cov_val;
	max_reads_num = query_seq_info_all.size();

//	for(i=0; i<query_seq_info_all.size(); i++) {
//		clip_aln = query_seq_info_all.at(i)->clip_aln;
//		if(clip_aln and clip_aln->bam and clip_aln->bam->core.qual>=minHighMapQ and clip_aln->bam->core.qual!=255)
//			high_qual_id_vec.push_back(i);
//	}
//
//	num = count = total_bases = 0;
//
//	if(high_qual_id_vec.size()>0){
//		// [1/4, 1/2]
//		min_id = high_qual_id_vec.size() / 4;
//		max_id = high_qual_id_vec.size() / 2;
//		for(k=max_id; k>=min_id; k--){
//			index = high_qual_id_vec.at(k);
//			qseq_node = query_seq_info_all.at(index);
//			if(qseq_node->selected_flag==false){
//				qname = qseq_node->qname;
//				qseq_node->selected_flag = true;
//				count ++;
//				if(selected_qname_vec.find(qname)==selected_qname_vec.end()){ // new item
//					selected_qname_vec.insert(qname);
//					total_bases += qseq_node->clip_aln->querylen;
//					num ++;
//					if(total_bases >= expected_total_bases) break;
//				}
//			}
//		}
//		// [1/2, 3/4]
//		if(total_bases < expected_total_bases){
//			min_id = high_qual_id_vec.size() / 2 + 1;
//			max_id = high_qual_id_vec.size() / 4 * 3;
//			if(max_id+1>(int32_t)high_qual_id_vec.size()) max_id = (int32_t)high_qual_id_vec.size() - 1;
//			for(k=min_id; k<=max_id; k++){
//				index = high_qual_id_vec.at(k);
//				qseq_node = query_seq_info_all.at(index);
//				if(qseq_node->selected_flag==false){
//					qname = qseq_node->qname;
//					qseq_node->selected_flag = true;
//					count ++;
//					if(selected_qname_vec.find(qname)==selected_qname_vec.end()){ // new item
//						selected_qname_vec.insert(qname);
//						total_bases += qseq_node->clip_aln->querylen;
//						num ++;
//						if(total_bases >= expected_total_bases) break;
//					}
//				}
//			}
//		}
//		// [0, 1/4]
//		if(total_bases < expected_total_bases){
//			min_id = 0;
//			max_id = high_qual_id_vec.size() / 4;
//			for(k=max_id; k>=min_id; k--){
//				index = high_qual_id_vec.at(k);
//				if(qseq_node->selected_flag==false){
//					qname = qseq_node->qname;
//					qseq_node->selected_flag = true;
//					count ++;
//					if(selected_qname_vec.find(qname)==selected_qname_vec.end()){ // new item
//						selected_qname_vec.insert(qname);
//						total_bases += qseq_node->clip_aln->querylen;
//						num ++;
//						if(total_bases >= expected_total_bases) break;
//					}
//				}
//			}
//		}
//		// [3/4, 1]
//		if(total_bases < expected_total_bases){
//			min_id = high_qual_id_vec.size() / 4 * 3 + 1;
//			max_id = high_qual_id_vec.size() - 1;
//			if(max_id+1>(int32_t)high_qual_id_vec.size()) max_id = (int32_t)high_qual_id_vec.size() - 1;
//			for(k=min_id; k<=max_id; k++){
//				index = high_qual_id_vec.at(k);
//				if(qseq_node->selected_flag==false){
//					qname = qseq_node->qname;
//					qseq_node->selected_flag = true;
//					count ++;
//					if(selected_qname_vec.find(qname)==selected_qname_vec.end()){ // new item
//						selected_qname_vec.insert(qname);
//						total_bases += qseq_node->clip_aln->querylen;
//						num ++;
//						if(total_bases >= expected_total_bases) break;
//					}
//				}
//			}
//		}
//	}
//	else{
//		cout << "high_qual_id_vec.size=" << high_qual_id_vec.size() << ", " << chrname << ":" << startRefPos_cns << "-" << endRefPos_cns << endl;
//	}

	//pthread_mutex_lock(&mutex_down_sample);
//	count = qual_num = 0;
//	while(total_bases<=expected_total_bases and qual_num<total_high_qual_num){
//		index = rand() % max_reads_num;
//		qseq_node = query_seq_info_all.at(index);
//		if(qseq_node->selected_flag==false){
//			clip_aln = qseq_node->clip_aln;
//			if(clip_aln and clip_aln->bam and clip_aln->bam->core.qual>=minHighMapQ and clip_aln->bam->core.qual!=255){
//				qseq_node->selected_flag = true;
//				if(selected_qname_vec.find(qseq_node->qname)==selected_qname_vec.end()){ // new item
//					selected_qname_vec.insert(qseq_node->qname);
//					total_bases += qseq_node->seq.size();
//					num ++;
//				}
//				count ++;
//				qual_num ++;
//			}
//		}
//	}

	srand(1);
	num = count = total_bases = 0;
	while(total_bases<=expected_total_bases and count<max_reads_num){
		index = rand() % max_reads_num;
		qseq_node = query_seq_info_all.at(index);
		if(qseq_node->selected_flag==false){
			qseq_node->selected_flag = true;
			if(selected_qname_vec.find(qseq_node->qname)==selected_qname_vec.end()){ // new item
				selected_qname_vec.insert(qseq_node->qname);
				total_bases += qseq_node->seq.size();
				num ++;
			}
			count ++;
		}
	}
	//pthread_mutex_unlock(&mutex_down_sample);

	// append remaining align items of the selected reads
	for(i=0; i<query_seq_info_all.size(); i++){
		qseq_node = query_seq_info_all.at(i);
		if(qseq_node->selected_flag==false){
			if(selected_qname_vec.find(qseq_node->qname)!=selected_qname_vec.end()) // found, then add item
				qseq_node->selected_flag = true;
		}
	}

	reads_count = num;
	sampled_cov = total_bases / reg_size;
	sampling_flag = true;

//	/cout << "After sampling, reads count: " << reads_count << ", total bases: " << total_bases << endl;

	// remove unselected items
	for(i=0; i<query_seq_info_all.size(); ){
		qseq_node = query_seq_info_all.at(i);
		if(qseq_node->selected_flag==false){
			destroyAlnSegs(qseq_node->query_alnSegs);
			delete qseq_node;
			query_seq_info_all.erase(query_seq_info_all.begin()+i);
		}else i++;
	}
}

void localCns::saveSampledReads(string &readsfilename, vector<struct fqSeqNode*> &fq_seq_vec){
	struct fqSeqNode *fq_node;
	ofstream outfile;
	size_t selected_num;

	outfile.open(readsfilename);
	if(!outfile.is_open()) {
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file " << readsfilename << endl;
		exit(1);
	}

	selected_num = 0;
	for(size_t i=0; i<fq_seq_vec.size(); i++){
		fq_node = fq_seq_vec.at(i);
		if(fq_node->selected_flag){
			outfile << "@" << fq_node->seq_name << endl;
			outfile << fq_node->seq << endl;
			outfile << "+" << endl;
			outfile << fq_node->qual << endl;
			selected_num ++;
		}
	}
	outfile.close();

	//cout <<"\t" << readsfilename << ": clip_aln_data_size=" << clipAlnDataVector.size() << ", reads_num=" << fq_seq_vec.size() << ", selected_num=" << selected_num << "; ref_size=" << endRefPos_cns-startRefPos_cns+1 << ", total_bases_original=" << total_bases_original << ", local_cov_original=" << local_cov_original << ", sampled_cov=" << sampled_cov << endl;
}

// mark Hard clipped segments from vector
void localCns::markHardClipSegs(size_t idx, vector<clipAlnData_t*> &query_aln_segs){
	for(size_t i=0; i<query_aln_segs.size(); ){
		if(i!=idx) {
			query_aln_segs.at(i)->query_checked_flag = true;
			query_aln_segs.erase(query_aln_segs.begin()+i);
		}else i++;
	}
}

// get query without soft clipped sequence
//vector<string> LocalCns::getSeqWithoutSoftClipSeqs(clipAlnData_t* clip_aln){
//	vector<string> query_info_vec;
//	uint8_t *seq_int, *qual_int;
//	string qseq, qual;
//	size_t i, query_aln_size, baseNum, startQPos, endQPos;
//
//	seq_int = bam_get_seq(clip_aln->bam);
//	qual_int = bam_get_qual(clip_aln->bam);
//
//	if(clip_aln->aln_orient==ALN_PLUS_ORIENT) query_aln_size = clip_aln->endQueryPos - clip_aln->startQueryPos + 1;
//	else query_aln_size = clip_aln->startQueryPos - clip_aln->endQueryPos + 1;
//
//	if(query_aln_size<clip_aln->bam->core.l_qseq) baseNum = query_aln_size;
//	else baseNum = clip_aln->bam->core.l_qseq;
//
//	startQPos = 0;
//	if(clip_aln->leftClipSize)
//
//	qseq = qual = "";
//
//
//	query_info_vec.push_back(qseq);
//	query_info_vec.push_back(qual);
//
//	return query_info_vec;
//}

// join query align segments
vector<string> localCns::joinQueryAlnSegs(vector<clipAlnData_t*> &query_aln_segs){
	string qseq, qual, qseq_result, qual_result, gap_seq, gap_qual;
	vector<string> query_seq_qual_vec;
	uint8_t *seq_int, *qual_int;
	int32_t j, join_orient, seq_len, left_most_idx, overlap_size, gap_size, startQueryPos, endQueryPos;
	clipAlnData_t *clip_aln, *pre_clip_aln;

	pre_clip_aln = NULL; join_orient = ALN_PLUS_ORIENT;
	startQueryPos = endQueryPos = -1;
	qseq_result = qual_result = "";
	for(size_t i=0; i<query_aln_segs.size(); i++){
		// get left most segment
		left_most_idx = getLeftMostAlnSeg(query_aln_segs);
		if(left_most_idx!=-1){
			clip_aln = query_aln_segs.at(left_most_idx);
			seq_int = bam_get_seq(clip_aln->bam);
			qual_int = bam_get_qual(clip_aln->bam);

			qseq = qual = "";
			for(j=0; j<clip_aln->bam->core.l_qseq; j++) qseq += "=ACMGRSVTWYHKDBN"[bam_seqi(seq_int, j)];  // seq
			for(j=0; j<clip_aln->bam->core.l_qseq; j++) qual += qual_int[j] + 33;  // qual

			cout << "qseq=" << qseq << endl;
			cout << "qual=" << qual << endl;

			clip_aln->query_checked_flag = true;

			if(pre_clip_aln){
				if(clip_aln->aln_orient!=join_orient){ // different orient
					reverseComplement(qseq);
					reverseSeq(qual);
				}

				if(join_orient==ALN_PLUS_ORIENT){ // plus orient
					if(clip_aln->aln_orient==ALN_PLUS_ORIENT) // same orient
						overlap_size = endQueryPos - clip_aln->startQueryPos + 1;
					else
						overlap_size = endQueryPos - clip_aln->endQueryPos + 1;

					if(overlap_size>0){
						qseq_result += qseq.substr(overlap_size);
						qual_result += qual.substr(overlap_size);
						endQueryPos += qseq.size() - overlap_size;
					}else{
						gap_size = -overlap_size;
						gap_seq = gap_qual = "";
						for(j=0; j<gap_size; j++) { gap_seq += 'N'; gap_qual += '!'; }
						qseq_result += gap_seq + qseq;
						qual_result += gap_qual + qual;
						endQueryPos += gap_size + qseq.size();
					}

				}else{ // minus orient
					if(clip_aln->aln_orient==ALN_MINUS_ORIENT) //  same orient
						overlap_size = startQueryPos - clip_aln->endQueryPos + 1;
					else
						overlap_size = startQueryPos - clip_aln->startQueryPos  + 1;

					if(overlap_size>0){
						seq_len = qseq.size() - overlap_size;
						qseq_result = qseq.substr(0, seq_len) + qseq_result;
						qual_result = qual.substr(0, seq_len) + qual_result;
						startQueryPos += qseq.size() - overlap_size;
					}else{
						gap_size = -overlap_size;
						gap_seq = gap_qual = "";
						for(j=0; j<gap_size; j++) { gap_seq += 'N'; gap_qual += '!'; }
						qseq_result = qseq + gap_seq + qseq_result;
						qual_result = qual + gap_seq + qual_result;
						startQueryPos += gap_size + qseq.size();
					}
				}
			}else{
				qseq_result = qseq;
				qual_result = qual;
				join_orient = clip_aln->aln_orient;  // join orient
				startQueryPos = clip_aln->startQueryPos;
				endQueryPos = clip_aln->endQueryPos;
			}
			pre_clip_aln = clip_aln;

		}else break;
	}

	cout << "qseq_result=" << qseq_result << endl;
	cout << "qual_result=" << qual_result << endl;

	return query_seq_qual_vec;
}

int32_t localCns::getLeftMostAlnSeg(vector<clipAlnData_t*> &query_aln_segs){
	int32_t left_most_idx;
	int64_t minPos;
	clipAlnData_t *clip_aln;

	left_most_idx = -1;
	minPos = INT_MAX;
	for(size_t i=0; i<query_aln_segs.size(); i++){
		clip_aln = query_aln_segs.at(i);
		if(clip_aln->query_checked_flag==false){
			if(clip_aln->startQueryPos<minPos){
				minPos = clip_aln->startQueryPos;
				left_most_idx = i;
			}
		}
	}

	return left_most_idx;
}

bool localCns::localConsensus(){
	bool flag;

	if(use_poa_flag){
		flag = cnsByPoa(); // local consensus using abPOA
		if(flag==false){
			flag = localCnsWtdbg2();
			//cout << "+++++++++ " << flag << ": rescued by wtdbg2, " << contigfilename << ", sv_len_est=" << sv_len_est << endl;
		}
	}else{
		flag = localCnsWtdbg2(); // local consensus using wtdbg2
		// if(flag==false){
		// 	flag = cnsByPoa();
		// 	//cout << "--------- " << flag << ": rescued by abpoa, " << contigfilename  << ", sv_len_est=" << sv_len_est << endl;
		// }
	}

	destorySeqsVec(seqs_vec);

	return flag;
}

// get heter local consensus using abPOA
bool localCns::cnsByPoa(){
	bool flag, abpoa_flag;
	string tmp_cons_filename, tmp_reads_filename, cons_header, cmd1, cmd2, cmd4;
	ofstream cons_file, reads_file;
	size_t k, i, j, n_seqs, serial_number;
	int64_t mem_avail;
	// check the file
	flag = isFileExist(contigfilename);
	if(flag) return true; // cons was generated successfully previously

	cmd1 = "mkdir -p " + tmpdir;
	cmd4 = "rm -rf " + tmpdir;  // delete files

	system(cmd1.c_str());

	// save the cons to file
	cons_file.open(contigfilename);
	if(!cons_file.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file " << contigfilename << endl;
		exit(1);
	}

	serial_number = 1;
	for(k=0; k<seqs_vec.size(); k++){
		tmp_reads_filename = readsfilename_prefix + "_" + to_string(k) + readsfilename_suffix;
		readsfilename_vec.push_back(tmp_reads_filename);  // save read file names
		reads_file.open(tmp_reads_filename);
		if(!reads_file.is_open()){
			cerr << __func__ << ", line=" << __LINE__ << ": cannot open file " << tmp_reads_filename << endl;
			exit(1);
		}

		n_seqs = seqs_vec.at(k)->seqs.size();
		if(n_seqs>=(size_t)min_supp_num*READS_NUM_SUPPORT_FACTOR){ // only use sufficient reads data to generate consensus sequence
			for (i = 0; i < n_seqs; ++i) reads_file << ">" << seqs_vec.at(k)->qname[i] << endl << seqs_vec.at(k)->seqs[i] << endl;
			reads_file.close();

			abpoa_flag = false;
			while(1){
				pthread_mutex_lock(&mutex_mem);
				mem_avail = getMemInfo("MemAvailable", 2);
				// prepare for alignment computation
				if(mem_avail>=min_mem_avail and mem_total-mem_avail<=mem_total*mem_use_block_factor+extend_total*extend_use_block_factor){
//				if(mem_seqAln+mem_cost<=mem_total*mem_use_block_factor+extend_total*extend_use_block_factor){
					work_num ++;
					abpoa_flag = true;
				}
				pthread_mutex_unlock(&mutex_mem);

				if(abpoa_flag==false){ // block the alignment computation
					//cout << "\t" << __func__ << ": wait " << mem_wait_seconds << " seconds, mem_avail=" << (mem_avail >> 10) << " MB, work_num=" << work_num << ", " << contigfilename << endl;
					sleep(mem_wait_seconds);
				}else{
					tmp_cons_filename = tmpdir + "/consensus_" + to_string(k) + ".fa";
					//cmd2 = "abpoa -o " + tmp_cons_filename + " " + tmp_reads_filename + " > /dev/null 2>&1";
					cmd2 = "abpoa";
					//if(mean_read_len>=MIN_SEQ_LEN_USING_MINIMIZER) cmd2 += " -S"; // deleted on 2025-03-16
					if(sv_len_est>MIN_SEQ_LEN_USING_MINIMIZER) cmd2 += " -S";
					cmd2 += " -o " + tmp_cons_filename + " " + tmp_reads_filename + " > /dev/null 2>&1";

					//cout << cmd2 << endl;
					system(cmd2.c_str());

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

			flag = isFileExist(tmp_cons_filename);
			if(flag){
				FastaSeqLoader fa_loader(tmp_cons_filename);
				for(i=0; i<fa_loader.getFastaSeqCount(); i++){
					if(fa_loader.getFastaSeqLen(i)>0){
						cons_header = ">abpoa_cns_";
						cons_header += to_string(serial_number) + "_"; //serial number
						cons_header += to_string(fa_loader.getFastaSeqLen(i)) + "_"; //consensus sequence length
						cons_header += to_string(seqs_vec.at(k)->seqs.size()) + "-";  // reads count

						for(j=0; j<seqs_vec.at(k)->qname.size()-1; j++) cons_header += seqs_vec.at(k)->qname.at(j) + "-";
						cons_header += seqs_vec.at(k)->qname.at(seqs_vec.at(k)->qname.size()-1);

						pthread_mutex_lock(&mutex_write);
						cons_file << cons_header << endl; // header
						cons_file << fa_loader.getFastaSeq(i) << endl;  // seq
						//cons_file.flush();
						pthread_mutex_unlock(&mutex_write);

						serial_number ++;
					}
				}
			}
		}else{
			reads_file.close();
		}
	}

	cons_file.close();
	//destorySeqsVec(seqs_vec);

	flag = isFileExist(contigfilename);
	if(flag){ // cons generated successfully
		cns_success_flag = true;
		//rename(tmp_cons_filename.c_str(), contigfilename.c_str());
	}
	system(cmd4.c_str());  // remove temporary files

	return flag;
}

// local consensus using wtdbg2
bool localCns::localCnsWtdbg2(){
	bool flag, wtdbg2_flag;
	string output_prefix, tmp_reg_str, tmp_cns_filename, cmd2, cmd3;
	string tmp_reads_filename, cons_header;
	size_t id, k, i, j, n_seqs, serial_number;
	ofstream cons_file, reads_file;
	int64_t mem_avail;

	// check the file
	flag = isFileExist(contigfilename);
	if(flag) return true; // cons was generated successfully previously

	id = tmpdir.find_last_of("/");
	if(id!=string::npos) tmp_reg_str = tmpdir.substr(id+1);
	else{
		cerr << "line=" << __LINE__ << ", cannot find the / in " << tmpdir << ", error." << endl;
		exit(1);
	}

	// save the cons to file
	cons_file.open(contigfilename);
	if(!cons_file.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file " << contigfilename << endl;
		exit(1);
	}

	serial_number = 1;
	for(k=0; k<seqs_vec.size(); k++){
		tmp_reads_filename = readsfilename_prefix + "_" + to_string(k) + readsfilename_suffix;
		readsfilename_vec.push_back(tmp_reads_filename);  // save read file names
		reads_file.open(tmp_reads_filename);
		if(!reads_file.is_open()){
			cerr << __func__ << ", line=" << __LINE__ << ": cannot open file " << tmp_reads_filename << endl;
			exit(1);
		}

		n_seqs = seqs_vec.at(k)->seqs.size();
		if(n_seqs>=(size_t)min_supp_num*READS_NUM_SUPPORT_FACTOR){ // only use sufficient reads data to generate consensus sequence
			for (i = 0; i < n_seqs; ++i) reads_file << ">" << seqs_vec.at(k)->qname[i] << endl << seqs_vec.at(k)->seqs[i] << endl;
			reads_file.close();

			wtdbg2_flag = false;
			while(1){
				pthread_mutex_lock(&mutex_mem);
				mem_avail = getMemInfo("MemAvailable", 2);
				// prepare for alignment computation
				if(mem_avail>=min_mem_avail and mem_total-mem_avail<=mem_total*mem_use_block_factor+extend_total*extend_use_block_factor){
					work_num ++;
					wtdbg2_flag = true;
				}
				pthread_mutex_unlock(&mutex_mem);

				if(wtdbg2_flag==false){ // block the alignment computation
					//cout << "\t" << __func__ << ": wait " << mem_wait_seconds << " seconds, mem_avail=" << (mem_avail >> 10) << " MB, work_num=" << work_num << ", " << contigfilename << endl;
					sleep(mem_wait_seconds);
				}else{
					output_prefix = "tmp_wtdbg2_" + tmp_reg_str;
					cmd2 = "wtdbg2.pl -t 8 -x " + technology + " -o " + output_prefix + " -a -q " + tmp_reads_filename + " > /dev/null 2>&1";
					//cout << cmd2 << endl;
					system(cmd2.c_str());

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
						cons_header = ">wtdbg2_cns_";
						cons_header += to_string(serial_number) + "_"; //serial number
						cons_header += to_string(fa_loader.getFastaSeqLen(i)) + "_"; //consensus sequence length
						cons_header += to_string(seqs_vec.at(k)->seqs.size()) + "-";  //read count

						for(j=0; j<seqs_vec.at(k)->qname.size()-1; j++) cons_header += seqs_vec.at(k)->qname.at(j) + "-";
						cons_header += seqs_vec.at(k)->qname.at(seqs_vec.at(k)->qname.size()-1);

						pthread_mutex_lock(&mutex_write);
						cons_file << cons_header << endl; // header
						cons_file << fa_loader.getFastaSeq(i) << endl;  // seq
						//cons_file.flush();
						pthread_mutex_unlock(&mutex_write);

						serial_number ++;
					}
				}
			}

			cmd3 = "rm -rf " + output_prefix + "*";
			system(cmd3.c_str());  // remove temporary files
		}else{
			reads_file.close();
		}
	}

	cons_file.close();
	//destorySeqsVec(seqs_vec);

	flag = isFileExist(contigfilename);
	if(flag){ // consensus generated successfully
		cns_success_flag = true;
		//rename(tmp_cns_filename.c_str(), contigfilename.c_str());
	}

	return flag;
}


// local consensus using Canu with more than one time
bool localCns::localCnsCanu(){
	bool flag = localCnsCanu_IncreaseGenomeSize();
	if(!flag) flag = localCnsCanu_DecreaseGenomeSize();
	//if(!flag) cout << "ASS_FAILURE: " << contigfilename << endl;
	return flag;
}

// local consensus using Canu with increasing genome size parameter
bool localCns::localCnsCanu_IncreaseGenomeSize(){
	string canu_cmd, cns_prefix, cmd2, cmd3, cmd4, tmp_ctg_filename, gnuplotTested_str, min_inout_cov_str, technology_str; // fast_option;
	int i, genomeSize_Canu, step_size;
	bool flag;
	string cmd_limited_threads_str, limited_threads_str;
	double cost_min;

	// gnuplotTested_str
	gnuplotTested_str = "";
	if(canu_version.compare(MIN_CANU_VERSION_NO_GNPPLOT)<0)
		gnuplotTested_str = " gnuplotTested=true ";

	min_inout_cov_str = " minInputCoverage=" + to_string(min_input_cov_canu) + " stopOnLowCoverage=" + to_string(min_input_cov_canu) + " ";

//	fast_option = "";
//	if(canu_version.compare("1.8")==0)
//		fast_option = " -fast";

	// check the file
	flag = isFileExist(contigfilename);
	if(flag) return true; // contig was generated successfully previously

	// generate command string
	cns_prefix = "cns";
	tmp_ctg_filename = tmpdir + "/" + cns_prefix + ".contigs.fasta";
	cmd2 = "mv " + tmp_ctg_filename + " " + contigfilename;   // move and rename file
	cmd4 = "rm -rf " + tmpdir;  // delete files

	// limited number threads
	cmd_limited_threads_str = "";
	if(num_threads_per_cns_work>0){
		limited_threads_str = to_string(num_threads_per_cns_work);
		cmd_limited_threads_str = " maxThreads=" + limited_threads_str + " obtovlThreads=" + limited_threads_str + " corThreads=" + limited_threads_str + " utgovlThreads=" + limited_threads_str + " redThreads=" + limited_threads_str + " batThreads=" + limited_threads_str + " ";
	}

	// technology
	technology_str = "";
	if(technology.compare(PACBIO_RS_TECH_STR)==0 or technology.compare(PACBIO_SQ_TECH_STR)==0){
		if(canu_version.compare(MIN_CANU_VERSION_NO_RAW)>=0)
			technology_str = " -pacbio ";
		else
			technology_str = " -pacbio-raw ";
	}else if(technology.compare(NANOPORE_TECH_STR)==0){
		if(canu_version.compare(MIN_CANU_VERSION_NO_RAW)>=0)
			technology_str = " -nanopore ";
		else
			technology_str = " -nanopore-raw ";
	}else if(technology.compare(PACBIO_CCS_TECH_STR)==0){
		if(canu_version.compare(MIN_CANU_VERSION_HIFI)>=0)
			technology_str = " -pacbio-hifi ";
		else{
			cerr << __func__ << ", line=" << __LINE__ << ": installed canu version is " << canu_version << ", however, it is required to use the canu version " << MIN_CANU_VERSION_HIFI << " or higher to consensus Pacbio CCS sequencing data." << endl;
			exit(1);
		}
	}else{
		cerr << __func__ << ", line=" << __LINE__ << ": invalid sequencing technology: " << technology << endl;
		exit(1);
	}

	// increase the genome size
	genomeSize_Canu = CNS_GENOME_SIZE_INITIAL;
	step_size = CNS_STEP_SIZE;
	flag = false;
	for(i=1; i<=3; i++){
		// try canu1.7
		////canu_cmd = "canu -p " + cns_prefix + " -d " + tmpdir + " genomeSize=" + to_string(genomeSize_Canu) + " gnuplotTested=true -pacbio-raw " + readsfilename + " > /dev/null 2>&1";
		//canu_cmd = "canu -p " + cns_prefix + " -d " + tmpdir + " genomeSize=" + to_string(genomeSize_Canu) + gnuplotTested_str + " -pacbio-raw " + readsfilename + " > /dev/null 2>&1";
		////canu_cmd = "canu -p " + cns_prefix + " -d " + tmpdir + " genomeSize=" + to_string(genomeSize_Canu) + fast_option + gnuplotTested_str + " -pacbio-raw " + readsfilename + " > /dev/null 2>&1";
//		canu_cmd = "canu1.7 -p " + cns_prefix + " -d " + tmpdir + cmd_limited_threads_str + " genomeSize=" + to_string(genomeSize_Canu) + " gnuplotTested=true -pacbio-raw " + readsfilename + " > /dev/null 2>&1";
//		system(canu_cmd.c_str());  // local consensus, invoke Canu command
//
//		// save consensus result and remove temporary files if successfully consensued
//		flag = isFileExist(tmp_ctg_filename);
//		if(flag){ // contig generated successfully
//			rename(tmp_ctg_filename.c_str(), contigfilename.c_str());
//			system(cmd4.c_str());  // remove temporary files
//			break;
//		}else { // contig generated failed
//			system(cmd4.c_str());  // remove temporary files

			time(&end_time);
			cost_min = difftime(end_time, start_time) / 60.0;
			if(cost_min<=MAX_CNS_MINUTES){
				// try canu1.8, or 2.0
				canu_cmd = "canu -p " + cns_prefix + " -d " + tmpdir + cmd_limited_threads_str + " genomeSize=" + to_string(genomeSize_Canu) + gnuplotTested_str + min_inout_cov_str + technology_str + readsfilename + " > /dev/null 2>&1";
				//cout << canu_cmd << endl;
				system(canu_cmd.c_str());  // local consensus, invoke Canu command

				// save consensus result and remove temporary files if successfully consensused
				flag = isFileExist(tmp_ctg_filename);
				if(flag){ // contig generated successfully
					cns_success_flag = true;
					rename(tmp_ctg_filename.c_str(), contigfilename.c_str());
					system(cmd4.c_str());  // remove temporary files
					break;
				}else { // contig generated failed
					system(cmd4.c_str());  // remove temporary files
					genomeSize_Canu += i * step_size;
				}
			}else break;
//		}
	}
	return flag;
}

// local consensus using Canu with decreasing genome size parameter
bool localCns::localCnsCanu_DecreaseGenomeSize(){
	string canu_cmd, cns_prefix, cmd2, cmd3, cmd4, tmp_ctg_filename, gnuplotTested_str, min_inout_cov_str, technology_str; //, fast_option;
	int i, genomeSize_Canu, step_size;
	bool flag;
	string cmd_limited_threads_str, limited_threads_str;
	double cost_min;

	// gnuplotTested_str
	gnuplotTested_str = "";
	if(canu_version.compare(MIN_CANU_VERSION_NO_GNPPLOT)<0)
		gnuplotTested_str = " gnuplotTested=true ";

	min_inout_cov_str = " minInputCoverage=" + to_string(min_input_cov_canu) + " stopOnLowCoverage=" + to_string(min_input_cov_canu) + " ";

//	fast_option = "";
//	if(canu_version.compare("1.8")==0 or canu_version.compare("1.9")==0)
//		fast_option = " -fast";

	// check the file
	flag = isFileExist(contigfilename);
	if(flag) return true;  // contig was generated successfully previously

	// generate command string
	cns_prefix = "cns";
	tmp_ctg_filename = tmpdir + "/" + cns_prefix + ".contigs.fasta";
	cmd2 = "mv " + tmp_ctg_filename + " " + contigfilename;   // move and rename file
	cmd4 = "rm -rf " + tmpdir;  // delete files

	// limited number threads
	cmd_limited_threads_str = "";
	if(num_threads_per_cns_work>0){
		limited_threads_str = to_string(num_threads_per_cns_work);
		cmd_limited_threads_str = " maxThreads=" + limited_threads_str + " obtovlThreads=" + limited_threads_str + " corThreads=" + limited_threads_str + " utgovlThreads=" + limited_threads_str + " redThreads=" + limited_threads_str + " batThreads=" + limited_threads_str + " ";
	}

	// technology
	technology_str = "";
	if(technology.compare(PACBIO_RS_TECH_STR)==0 or technology.compare(PACBIO_SQ_TECH_STR)==0){
		if(canu_version.compare(MIN_CANU_VERSION_NO_RAW)>=0)
			technology_str = " -pacbio ";
		else
			technology_str = " -pacbio-raw ";
	}else if(technology.compare(NANOPORE_TECH_STR)==0){
		if(canu_version.compare(MIN_CANU_VERSION_NO_RAW)>=0)
			technology_str = " -nanopore ";
		else
			technology_str = " -nanopore-raw ";
	}else if(technology.compare(PACBIO_CCS_TECH_STR)==0){
		if(canu_version.compare(MIN_CANU_VERSION_HIFI)>=0)
			technology_str = " -pacbio-hifi ";
		else{
			cerr << __func__ << ", line=" << __LINE__ << ": installed canu version is " << canu_version << ", however, it is required to use the canu version " << MIN_CANU_VERSION_HIFI << " or higher to consensues Pacbio CCS sequencing data." << endl;
			exit(1);
		}
	}else{
		cerr << __func__ << ", line=" << __LINE__ << ": invalid sequencing technology: " << technology << endl;
		exit(1);
	}

	// decrease the genome size
	step_size = CNS_STEP_SIZE;
	genomeSize_Canu = CNS_GENOME_SIZE_INITIAL - step_size;
	flag = false;
	for(i=1; i<=3 and genomeSize_Canu>=step_size; i++){
		// try canu1.7
		////canu_cmd = "canu -p " + cns_prefix + " -d " + tmpdir + " genomeSize=" + to_string(genomeSize_Canu) + " gnuplotTested=true -pacbio-raw " + readsfilename + " > /dev/null 2>&1";
		//canu_cmd = "canu -p " + cns_prefix + " -d " + tmpdir + " genomeSize=" + to_string(genomeSize_Canu) + gnuplotTested_str + " -pacbio-raw " + readsfilename + " > /dev/null 2>&1";
		////canu_cmd = "canu -p " + cns_prefix + " -d " + tmpdir + " genomeSize=" + to_string(genomeSize_Canu) + fast_option + gnuplotTested_str + " -pacbio-raw " + readsfilename + " > /dev/null 2>&1";
//		canu_cmd = "canu1.7 -p " + cns_prefix + " -d " + tmpdir + cmd_limited_threads_str + " genomeSize=" + to_string(genomeSize_Canu) + " gnuplotTested=true -pacbio-raw " + readsfilename + " > /dev/null 2>&1";
//		system(canu_cmd.c_str());  // local consensus, invoke Canu command
//
//		// save consensus result and remove temporary files if successfully consensued
//		flag = isFileExist(tmp_ctg_filename);
//		if(flag){ // contig generated successfully
//			rename(tmp_ctg_filename.c_str(), contigfilename.c_str());
//			system(cmd4.c_str());  // remove temporary files
//			break;
//		}else { // contig generated failed
//			system(cmd4.c_str());  // remove temporary files

			time(&end_time);
			cost_min = difftime(end_time, start_time) / 60.0;
			if(cost_min<=MAX_CNS_MINUTES){
				// try canu1.8, or 2.0
				canu_cmd = "canu -p " + cns_prefix + " -d " + tmpdir + cmd_limited_threads_str + " genomeSize=" + to_string(genomeSize_Canu) + gnuplotTested_str + min_inout_cov_str + technology_str + readsfilename + " > /dev/null 2>&1";
				//cout << canu_cmd << endl;
				system(canu_cmd.c_str());  // local consensus, invoke Canu command

				// save consensus result and remove temporary files if successfully consensused
				flag = isFileExist(tmp_ctg_filename);
				if(flag){ // contig generated successfully
					cns_success_flag = true;
					rename(tmp_ctg_filename.c_str(), contigfilename.c_str());
					system(cmd4.c_str());  // remove temporary files
					break;
				}else { // contig generated failed
					system(cmd4.c_str());  // remove temporary files
					genomeSize_Canu -= i * step_size;
				}
			}else break;
//		}
	}
	return flag;
}

// record consensus information
void localCns::recordCnsInfo(ofstream &cns_info_file){
	string line, cns_status, header, left_shift_size_str, right_shift_size_str, reg_str, sampling_str, limit_reg_str, limit_reg_str2;
	reg_t *reg;
	ifstream infile;
	vector<string> str_vec;
	simpleReg_t *simple_reg;

	// ref shift size
	infile.open(refseqfilename);
	if(!infile.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file " << refseqfilename << endl;
		exit(1);
	}

	getline(infile, header);
	str_vec = split(header, "___");
	left_shift_size_str = str_vec[1];
	right_shift_size_str = str_vec[2];

	infile.close();

	// consensus status
	if (isFileExist(contigfilename)) cns_status = CNS_SUCCESS;
	else cns_status = CNS_FAILURE;

	line = refseqfilename + "\t" + contigfilename + "\t" + readsfilename + "\t" + clusterfilename + "\t" + left_shift_size_str + "\t" + right_shift_size_str + "\t" + cns_status;

	reg_str = "";
	if(varVec.size()>0){
		reg = varVec.at(0);
		reg_str = reg->chrname + ":" + to_string(reg->startRefPos) + "-" + to_string(reg->endRefPos);
		for(size_t i=1; i<varVec.size(); i++){
			reg = varVec.at(i);
			reg_str += ";" + reg->chrname + ":" + to_string(reg->startRefPos) + "-" + to_string(reg->endRefPos);
		}
		line += "\t" + reg_str;
	}else{
		line += "\t-";
	}

	// coverage sampling information
	if(sampling_flag)
		sampling_str = to_string(local_cov_original) + ";" + to_string(sampled_cov) + ";" + to_string(compensation_coefficient) + ";" + SAMPLED_STR;
	else{
		sampling_str = "-;-;-;-;";
		sampling_str = sampling_str + UNSAMPLED_STR;
	}
	line += "\t" + sampling_str;

	// limit process regions
	if(limit_reg_process_flag){
		if(limit_reg_vec.size()){
			simple_reg = limit_reg_vec.at(0);
			limit_reg_str = simple_reg->chrname;
			if(simple_reg->startPos!=-1 and simple_reg->endPos!=-1) limit_reg_str += ":" + to_string(simple_reg->startPos) + "-" + to_string(simple_reg->endPos);
			for(size_t i=1; i<limit_reg_vec.size(); i++){
				simple_reg = limit_reg_vec.at(i);
				limit_reg_str2 = simple_reg->chrname;
				if(simple_reg->startPos!=-1 and simple_reg->endPos!=-1) limit_reg_str2 += ":" + to_string(simple_reg->startPos) + "-" + to_string(simple_reg->endPos);
				limit_reg_str += ";" + limit_reg_str2;
			}
		}else limit_reg_str = "-";
	}else{
		limit_reg_str = LIMIT_REG_ALL_STR;
	}
	line += "\t" + limit_reg_str;

	// done string
	line = line + "\t" + DONE_STR;

	pthread_mutex_lock(&mutex_write);
	cns_info_file << line << endl;
	pthread_mutex_unlock(&mutex_write);
}
