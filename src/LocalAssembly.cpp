#include "LocalAssembly.h"
#include "clipAlnDataLoader.h"

pthread_mutex_t mutex_write = PTHREAD_MUTEX_INITIALIZER;
extern pthread_mutex_t mutex_down_sample;
//static int32_t MIN_SVLEN = 20;

LocalAssembly::LocalAssembly(string &readsfilename, string &contigfilename, string &refseqfilename, string &tmpdir, string &technology, string &canu_version, size_t num_threads_per_assem_work, vector<reg_t*> &varVec, string &chrname, string &inBamFile, faidx_t *fai, size_t assembly_extend_size, double expected_cov, double min_input_cov, bool delete_reads_flag, bool keep_failed_reads_flag, bool clip_reg_flag, int32_t minClipEndSize, int32_t minConReadLen, int32_t min_sv_size, double max_seg_size_ratio){
	this->chrname = chrname;
	this->chrlen = faidx_seq_len(fai, chrname.c_str()); // get reference size
	this->readsfilename = preprocessPipeChar(readsfilename);
	this->contigfilename = preprocessPipeChar(contigfilename);
	this->refseqfilename = preprocessPipeChar(refseqfilename);
	this->tmpdir = preprocessPipeChar(tmpdir);
	this->technology = technology;
	this->canu_version = canu_version;
	this->num_threads_per_assem_work = num_threads_per_assem_work;
	this->minClipEndSize = minClipEndSize;
	this->minConReadLen = minConReadLen;
	this->min_sv_size = min_sv_size;
	this->max_seg_size_ratio = max_seg_size_ratio;
	this->varVec = varVec;
	this->fai = fai;
	this->inBamFile = inBamFile;
	//this->assembly_extend_size = ASSEMBLE_SIDE_EXT_SIZE + assembly_extend_size;
	this->assembly_extend_size = assembly_extend_size;
	startRefPos_assembly = endRefPos_assembly = 0;
	mean_read_len = 0;
	use_poa_flag = true;

	ref_seq_size = reads_count_original = total_bases_original = reads_count = total_bases = 0;
	local_cov_original = sampled_cov = 0;
	this->expected_cov = expected_cov;
	this->compensation_coefficient = 1;
	sampling_flag = false;
	this->delete_reads_flag = delete_reads_flag;
	this->keep_failed_reads_flag = keep_failed_reads_flag;
	assem_success_flag = false;
	this->clip_reg_flag = clip_reg_flag;

	if(isFileExist(contigfilename) and isFileExist(refseqfilename)) assem_success_preDone_flag = true;
	else assem_success_preDone_flag = false;

	limit_reg_process_flag = false;
	this->min_input_cov_canu = min_input_cov;

	time(&start_time);
	end_time = 0;
}

LocalAssembly::~LocalAssembly() {
	if(delete_reads_flag){
		if(keep_failed_reads_flag==false or (keep_failed_reads_flag and assem_success_flag))
			remove(readsfilename.c_str());	// delete the reads file to save disk space
	}
}

void LocalAssembly::destoryAlnData(){
	for(size_t i=0; i<alnDataVector.size(); i++)
		bam_destroy1(alnDataVector.at(i));
	vector<bam1_t*>().swap(alnDataVector);
}

// destroy the alignment data of the block
void LocalAssembly::destoryClipAlnData(){
	for(size_t i=0; i<clipAlnDataVector.size(); i++){
		bam_destroy1(clipAlnDataVector.at(i)->bam);
		delete clipAlnDataVector.at(i);
	}
	vector<clipAlnData_t*>().swap(clipAlnDataVector);
}

void LocalAssembly::destoryFqSeqs(vector<struct fqSeqNode*> &fq_seq_vec){
	struct fqSeqNode *fq_node;
	for(size_t i=0; i<fq_seq_vec.size(); i++){
		fq_node = fq_seq_vec.at(i);
		delete fq_node;
	}
	vector<struct fqSeqNode*>().swap(fq_seq_vec);
}
void LocalAssembly::destorySeqsVec(vector<struct seqsVec*> &seqs_vec){
	struct seqsVec* seqs_vec_node;
	for(size_t i=0; i<seqs_vec.size(); i++){
		seqs_vec_node = seqs_vec.at(i);
		delete seqs_vec_node;
	}
	vector<struct seqsVec*>().swap(seqs_vec);
}
void LocalAssembly::destoryQueryCluVec(vector<struct querySeqInfoVec*> &query_clu_vec){
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
void LocalAssembly::destoryQuerySeqInfoAll(vector<struct querySeqInfoNode*> &query_seq_info_all){
	struct querySeqInfoNode *q_node;
	//vector<struct alnSeg*> query_alnSegs;
	struct alnSeg *q_alnSeg;
	for(size_t i=0; i<query_seq_info_all.size(); i++){
		q_node = query_seq_info_all.at(i);
		//destroyAlnSegs(query_seq_info_all.at(i)->query_alnSegs);
		for(size_t j=0; j<query_seq_info_all.at(i)->query_alnSegs.size(); j++){
			q_alnSeg = query_seq_info_all.at(i)->query_alnSegs.at(j);
			delete q_alnSeg;
		}
		vector<struct alnSeg*>().swap(query_seq_info_all.at(i)->query_alnSegs);
		delete q_node;
	}
	vector<struct querySeqInfoNode*>().swap(query_seq_info_all);
}

// set limit process regions
void LocalAssembly::setLimitRegs(bool limit_reg_process_flag, vector<simpleReg_t*> limit_reg_vec){
	this->limit_reg_process_flag = limit_reg_process_flag;
	for(size_t i=0; i<limit_reg_vec.size(); i++) this->limit_reg_vec.push_back(limit_reg_vec.at(i));
}

// extract the corresponding refseq from reference
void LocalAssembly::extractRefseq(){
	int32_t startRefPos, endRefPos, left_shift_size, right_shift_size;
	string reg, header;
	ofstream outfile;

	// generate the refseq region
	//startRefPos = varVec[0]->startRefPos - REFSEQ_SIDE_LEN;
	startRefPos = varVec[0]->startRefPos - assembly_extend_size;
	if(startRefPos<1) startRefPos = 1;
	left_shift_size = varVec[0]->startRefPos - startRefPos;  // left shift size
	//endRefPos = varVec[varVec.size()-1]->endRefPos + REFSEQ_SIDE_LEN;
	endRefPos = varVec[varVec.size()-1]->endRefPos + assembly_extend_size;
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
void LocalAssembly::extractReadsDataFromBAM(){
	size_t i, j, k, seq_id;
	string qname, seq, qual, reg_str, refseq;
	vector<struct querySeqInfoVec*> query_clu_vec;
	vector<clipAlnData_t*> query_aln_segs;
	vector<string> query_seq_qual_vec;
	vector<struct fqSeqNode*> fq_seq_vec; // [0]: query name; [1]: sequence; [2]: query name (optional); [3]: quality
	struct fqSeqNode *fq_node;
	vector<querySeqInfoNode*> query_seq_info_all;
	struct querySeqInfoNode *q_node;
	vector<struct query_seq_info_cluster*> query_seq_info_all_cluster;
	int64_t noHardClipIdx;
	int32_t seq_len, bam_type;
	char *p_seq;
	//bool no_otherchrname_flag;
	struct seqsVec * smoothed_seqs;
	struct seqsVec *seqs_vec_node;
	struct alnSeg *q_alnSeg;
	struct querySeqInfoVec *q_cluster_node;
	vector<struct alnSeg*> alnSegs;

	startRefPos_assembly = varVec[0]->startRefPos - assembly_extend_size;
	if(startRefPos_assembly<1) startRefPos_assembly = 1;
	endRefPos_assembly = varVec[varVec.size()-1]->endRefPos + assembly_extend_size;
	if(endRefPos_assembly>chrlen) endRefPos_assembly = chrlen;

	// load the aligned reads data
//	clipAlnDataLoader data_loader(varVec[0]->chrname, startRefPos_assembly, endRefPos_assembly, inBamFile, minClipEndSize);
//	data_loader.loadClipAlnData(clipAlnDataVector);

	// load the clipping data ---20220705
	clipAlnDataLoader data_loader(varVec[0]->chrname, startRefPos_assembly, endRefPos_assembly, inBamFile, minClipEndSize);
	//data_loader.loadClipAlnDataWithSATag(clipAlnDataVector, 0);
	if(clip_reg_flag) data_loader.loadClipAlnDataWithSATag(clipAlnDataVector);
	else data_loader.loadClipAlnDataWithSATagWithSegSize(clipAlnDataVector, 0, max_seg_size_ratio);

	if(clipAlnDataVector.size()>0){// and startRefPos_assembly==12377901){

		//cout << "start_pos_assembly=" << startRefPos_assembly << ", end_pos_assembly=" << endRefPos_assembly << ", clipAlnDataVector.size=" << clipAlnDataVector.size() << endl;

		for(i=0; i<clipAlnDataVector.size(); i++) clipAlnDataVector.at(i)->query_checked_flag = false;

		reg_str = varVec[0]->chrname + ":" + to_string(startRefPos_assembly) + "-" + to_string(endRefPos_assembly);
		p_seq = fai_fetch(fai, reg_str.c_str(), &seq_len);
		refseq = p_seq;
		free(p_seq);

		bam_type = getBamType(clipAlnDataVector.at(0)->bam);
		if(bam_type==BAM_INVALID){
			cerr << __func__ << ": unknown bam type, error!" << endl;
			exit(1);
		}

		// join query clip align segments
		seq_id = 0;
		total_bases_original = 0;
		for(i=0; i<clipAlnDataVector.size(); i++){
			if(clipAlnDataVector.at(i)->query_checked_flag==false){
				qname = clipAlnDataVector.at(i)->queryname;
				query_aln_segs = getQueryClipAlnSegs(qname, clipAlnDataVector);  // get query clip align segments

//				if(qname.compare("S1_168094")==0){
//					cout << "qname=" << qname << endl;
//				}else continue;

			//	no_otherchrname_flag = true;
				noHardClipIdx = getNoHardClipAlnItem(query_aln_segs);
				if(noHardClipIdx!=-1){
//					for(size_t t=0; t<query_aln_segs.size(); t++){
//						if(t!=noHardClipIdx) {
//							if(query_aln_segs.at(t)->chrname.compare(varVec[0]->chrname)!=0){
//								no_otherchrname_flag = false;
//							}
//						}
//					}
					//if(no_otherchrname_flag){
						//for(j=0; j<query_aln_segs.size(); j++){
					query_seq_qual_vec = getQuerySeqWithSoftClipSeqs(query_aln_segs.at(noHardClipIdx), refseq, bam_type);
					//markHardClipSegs(noHardClipIdx, query_aln_segs);
					query_aln_segs.at(noHardClipIdx)->query_checked_flag = true;

					if(query_seq_qual_vec.size()>0 and query_seq_qual_vec.at(0).length()>minConReadLen){
						seq = query_seq_qual_vec.at(0);
						if(clip_reg_flag){
							qual = query_seq_qual_vec.at(1);
							fq_node = new struct fqSeqNode();
							fq_node->seq_id = seq_id++;
							fq_node->seq_name = qname;
							fq_node->seq = seq;
							fq_node->qual = qual;
							fq_node->selected_flag = true;
							fq_seq_vec.push_back(fq_node);

							total_bases_original += seq.size();
						}else{
							//seqs.push_back(seq);
							//prepare for smoothing
							seq = query_seq_qual_vec.at(0);
							q_node = new querySeqInfoNode();
							q_node->qname = qname;
							q_node->seq = seq;
							//q_node->seq_id = fq_node->seq_id;
							switch(bam_type){
								case BAM_CIGAR_NO_DIFF_MD:
									//alnSegs = generateAlnSegs(b);
									alnSegs = generateAlnSegs2(query_aln_segs.at(noHardClipIdx)->bam, startRefPos_assembly, endRefPos_assembly);
									break;
								case BAM_CIGAR_NO_DIFF_NO_MD:
								case BAM_CIGAR_DIFF_MD:
								case BAM_CIGAR_DIFF_NO_MD:
									//alnSegs = generateAlnSegs_no_MD(b, refseq, startPos, endPos);
									alnSegs = generateAlnSegs_no_MD2(query_aln_segs.at(noHardClipIdx)->bam, refseq, startRefPos_assembly, endRefPos_assembly);
									break;
								default:
									cerr << __func__ << ": unknown bam type, error!" << endl;
									exit(1);
							}// generate align segments
							//if(clip_reg_flag==false)
							removeHeadTailAlnSegs(alnSegs);
							q_node->query_alnSegs = alnSegs;

							query_seq_info_all.push_back(q_node);
						}
					}
						//}
				//	}
				}else{ // all hard clipping
					//cout << "qname=" << qname << ", querylen=" << clipAlnDataVector.at(i)->querylen << endl;
					// join query align segments without 'SA' tags
//					query_seq_qual_vec = joinQueryAlnSegs(query_aln_segs);
					for(j=0; j<query_aln_segs.size(); j++){
						query_seq_qual_vec = getQuerySeqWithSoftClipSeqs(query_aln_segs.at(j), refseq, bam_type);
						if(query_seq_qual_vec.size()>0 and query_seq_qual_vec.at(0).length()>minConReadLen){
							seq = query_seq_qual_vec.at(0);
							if(clip_reg_flag){
								qual = query_seq_qual_vec.at(1);
								fq_node = new struct fqSeqNode();
								fq_node->seq_id = seq_id++;
								fq_node->seq_name = qname;
								fq_node->seq = seq;
								fq_node->qual = qual;
								fq_node->selected_flag = true;
								fq_seq_vec.push_back(fq_node);

								total_bases_original += seq.size();
							}else{
								//seqs.push_back(seq);
								//prepare for smoothing
								q_node = new struct querySeqInfoNode();
								q_node->qname = qname;
								q_node->seq = seq;
								//q_node->seq_id = fq_node->seq_id;

								switch(bam_type){
									case BAM_CIGAR_NO_DIFF_MD:
										//alnSegs = generateAlnSegs(b);
										alnSegs = generateAlnSegs2(query_aln_segs.at(j)->bam, startRefPos_assembly, endRefPos_assembly);
										break;
									case BAM_CIGAR_NO_DIFF_NO_MD:
									case BAM_CIGAR_DIFF_MD:
									case BAM_CIGAR_DIFF_NO_MD:
										//alnSegs = generateAlnSegs_no_MD(b, refseq, startPos, endPos);
										alnSegs = generateAlnSegs_no_MD2(query_aln_segs.at(j)->bam, refseq, startRefPos_assembly, endRefPos_assembly);
										break;
									default:
										cerr << __func__ << ": unknown bam type, error!" << endl;
										exit(1);
								}// generate align segments
								//q_node->query_alnSegs = generateAlnSegs2(query_aln_segs.at(j)->bam,startRefPos_assembly, endRefPos_assembly,);
								//if(clip_reg_flag==false)
								removeHeadTailAlnSegs(alnSegs);
								q_node->query_alnSegs = alnSegs;
								query_seq_info_all.push_back(q_node);
							}
							//destroyAlnSegs(query_alnSegs);
						}
					}
					for(j=0; j<query_aln_segs.size(); j++) query_aln_segs.at(j)->query_checked_flag = true;
				}
			}
		}

		if(clip_reg_flag){ // clipping region
			reads_count_original = fq_seq_vec.size();
			mean_read_len = (double)total_bases_original / fq_seq_vec.size();

			//cout << "mean_read_len=" << mean_read_len << endl;

			// sampling to expected coverage
			if(expected_cov!=0){
				compensation_coefficient = computeCompensationCoefficient(startRefPos_assembly, endRefPos_assembly, mean_read_len);
				//cout << compensation_coefficient << endl;
				samplingReads(fq_seq_vec, expected_cov, compensation_coefficient);
			}

			// save sampled reads to file
			saveSampledReads(readsfilename, fq_seq_vec);

		}else{ // indel region
			if(query_seq_info_all.size()){
				if(query_seq_info_all.size()>1){
					query_clu_vec = queryCluster(query_seq_info_all);
				}else{
					q_cluster_node = new struct querySeqInfoVec();
					q_cluster_node->query_seq_info_all = query_seq_info_all;
					query_clu_vec.push_back(q_cluster_node);
				}
				//smooth read segments
				for(j=0;j<query_clu_vec.size();j++){
					smoothed_seqs = smoothQuerySeqData(refseq, query_clu_vec.at(j)->query_seq_info_all);
					//seqs_vec smoothed_seqs;
					//seqs_vec_node = new struct seqsVec();
					seqs_vec_node = smoothed_seqs;
					seqs_vec.push_back(seqs_vec_node);
				}
			}
		}


		//not smooth read segments
	//	for(j=0;j<query_clu_vec.size();j++){
	//		for(k=0;k<query_clu_vec.at(j)->query_seq_info_all.size();k++){
	//			smoothed_seqs.push_back(query_clu_vec.at(j)->query_seq_info_all.at(k)->seq);
	//		}
	//		seqs_vec_node = new struct seqsVec();
	//		seqs_vec_node->seqs = smoothed_seqs;
	//		seqs_vec.push_back(seqs_vec_node);
	//	}
	}
	if(!clipAlnDataVector.empty()) destoryClipAlnData();


	// release memory
	if(!fq_seq_vec.empty()) destoryFqSeqs(fq_seq_vec);
	//if(!query_alnSegs.empty()) destroyAlnSegs(query_alnSegs);
	if(!query_seq_info_all.empty()) destoryQuerySeqInfoAll(query_seq_info_all);
	if(!query_clu_vec.empty()) destoryQueryCluVec(query_clu_vec);

}
void LocalAssembly::destroyQueryQcSig(queryCluSig_t *qc_Sig){
	qcSig_t *qcluSig_t;
	for(size_t i=0; i<qc_Sig->qcSig_vec.size(); i++){
		qcluSig_t = qc_Sig->qcSig_vec.at(i);
		delete qcluSig_t;
	}
	vector<qcSig_t*>().swap(qc_Sig->qcSig_vec);
	delete qc_Sig;
}

vector<struct querySeqInfoVec*> LocalAssembly::queryCluster(vector<struct querySeqInfoNode*> &query_seq_info_all){
	vector<struct querySeqInfoVec*> query_clu_vec;
	struct querySeqInfoVec *q_cluster_node_a;
	struct querySeqInfoVec *q_cluster_node_b;
	struct querySeqInfoVec *q_cluster_node_rescue;
	vector<struct querySeqInfoNode*> q_cluster_a, q_cluster_b, q_cluster_rescue, q_cluster_a_final, q_cluster_b_final;
	//vector<queryCluSig_t*> queryCluSig_vec;
	int32_t query_endRpos, query_startRpos, span, span_max, span_next, span_startpos, span_endpos;
	struct querySeqInfoNode *query_seq_info_item;
	struct seedQueryInfo *seed_info_a;
	struct seedQueryInfo *seed_info_b;
	size_t i, j, k, seed_a_id,seed_b_id, smin_id, smin_id2, span_max_id;
	double score_ratio_a, score_ratio_b, sco_min, same_min_num;
	queryCluSig_t *queryCluSig, *queryCluSig1;
	int32_t max_csig_num, max_csig_id, srsame_maxnum;
	vector<size_t> min_id_vec;// span_vec2_id;
	vector<double> score_ratio_vec;
	//vector<int32_t> span_vec1, span_vec2;
	vector<srmin_match_t*> srmin_match_vec;
	srmin_match_t *srmin_match_node;
	bool exist_flag;


	for(i=0; i<query_seq_info_all.size(); i++){
		query_seq_info_all.at(i)->cluster_finished_flag = false;
//		if(query_seq_info_all.at(i)->qname.compare("S1_68264")==0){
//			cout<<"S1_68264: "<<i<<endl;
//		}
	}
	//sort in ascending order by |query_endRpos - query_startRpos|
	//query_seq_info_item = new struct querySeqInfoNode();
	for(i=0; i<query_seq_info_all.size(); i++){
		for(j=0; j<query_seq_info_all.size()-i-1; j++){
		span = query_seq_info_all.at(j)->query_alnSegs.at(query_seq_info_all.at(j)->query_alnSegs.size() - 1)->startRpos +  query_seq_info_all.at(j)->query_alnSegs.at(query_seq_info_all.at(j)->query_alnSegs.size() - 1)->seglen - query_seq_info_all.at(j)->query_alnSegs.at(0)->startRpos;
		span_next = query_seq_info_all.at(j+1)->query_alnSegs.at(query_seq_info_all.at(j+1)->query_alnSegs.size() - 1)->startRpos +  query_seq_info_all.at(j+1)->query_alnSegs.at(query_seq_info_all.at(j+1)->query_alnSegs.size() - 1)->seglen - query_seq_info_all.at(j+1)->query_alnSegs.at(0)->startRpos;
		if(span<span_next){ // swap
			query_seq_info_item = query_seq_info_all.at(j);
			query_seq_info_all.at(j) = query_seq_info_all.at(j+1);
			query_seq_info_all.at(j+1) = query_seq_info_item;
			}
		}
	}
	//delete query_seq_info_item;
//	for(i=0; i<query_seq_info_all.size(); i++){
//		if(query_seq_info_all.at(i)->qname.compare("S1_168094")==0){
//			cout<<"S1_168094: "<<i<<endl;
//		}
//	}
	for(k=0;k<query_seq_info_all.size();k++){
		//initializing a_cluster
		query_seq_info_all.at(k)->cluster_finished_flag = true;
		q_cluster_a.push_back(query_seq_info_all.at(k));

		//initializing b_cluster: choose one that is most different query from a_cluster add in b_cluster
		sco_min = 0.99;
		score_ratio_vec.clear();
		min_id_vec.clear();
		for(i=0; i<query_seq_info_all.size(); i++){
			if(query_seq_info_all.at(i)->cluster_finished_flag == false){
				seed_info_a = chooseSeedClusterQuery(query_seq_info_all.at(i), q_cluster_a);
				//span_vec1.push_back(seed_info_a->span_intersection);
				if(seed_info_a->span_intersection>0){
					score_ratio_a = computeScoreRatio(query_seq_info_all.at(i), q_cluster_a.at(seed_info_a->id), seed_info_a->startSpanPos, seed_info_a->endSpanPos);
					score_ratio_vec.push_back(score_ratio_a);
					if(score_ratio_a<sco_min){
						sco_min = score_ratio_a;
						//smin_id = i;
						//sam_min_id.push_back(i);
					}
				}else{
					score_ratio_vec.push_back(2);
				}
				delete seed_info_a;
			}else{
				//span_vec1.push_back(0);
				score_ratio_vec.push_back(1);
			}
		}

		exist_flag = false;
		if(sco_min < 0.99){
			for(i=0;i<score_ratio_vec.size();i++){
				if(score_ratio_vec.at(i)<0.99)
					min_id_vec.push_back(i);
			}
			for(i=0; i<min_id_vec.size(); i++){
				for(j=0; j<srmin_match_vec.size(); j++){
					if(score_ratio_vec.at(min_id_vec.at(i))==srmin_match_vec.at(j)->SR){
						smin_id = j;
						exist_flag = true;
					}else{
						exist_flag = false;
					}
				}
				if(exist_flag==false){
					srmin_match_node = new srmin_match_t();
					srmin_match_node->pos_id.push_back(min_id_vec.at(i));
					srmin_match_node->num = 1;
					srmin_match_node->SR = score_ratio_vec.at(min_id_vec.at(i));
					srmin_match_vec.push_back(srmin_match_node);
				}else{
					srmin_match_vec.at(smin_id)->num += 1;
					srmin_match_vec.at(smin_id)->pos_id.push_back(min_id_vec.at(i));
				}
			}

			srsame_maxnum = srmin_match_vec.at(0)->num;
			smin_id = 0;
			for(i=0; i<srmin_match_vec.size(); i++){
				if(srsame_maxnum<srmin_match_vec.at(i)->num){
					smin_id = i;
					srsame_maxnum = srmin_match_vec.at(i)->num;
				}
			}

			query_seq_info_all.at(srmin_match_vec.at(smin_id)->pos_id.at(0))->cluster_finished_flag = true;
			q_cluster_b.push_back(query_seq_info_all.at(srmin_match_vec.at(smin_id)->pos_id.at(0)));
		}

		if(srmin_match_vec.size()>0){
			for(i=0; i<srmin_match_vec.size(); i++){
				vector<srmin_match_t*>::iterator svp;
				for(svp=srmin_match_vec.begin(); svp!=srmin_match_vec.end(); svp++) delete *svp;
				vector<srmin_match_t*>().swap(srmin_match_vec);
			}
		}

		if(q_cluster_b.size()==0){
			query_seq_info_all.at(k)->cluster_finished_flag = false;
		}else{
			break;
		}
	}

	if(q_cluster_b.size()==0){
		for(i=0; i<query_seq_info_all.size(); i++){
			if(query_seq_info_all.at(i)->cluster_finished_flag == false){
				queryCluSig = new queryCluSig_t();
				queryCluSig->qcSig_vec = extractQcSigsFromAlnSegsSingleQuery(query_seq_info_all.at(i), query_seq_info_all.at(i)->query_alnSegs.at(0)->startRpos, query_seq_info_all.at(i)->query_alnSegs.at(query_seq_info_all.at(i)->query_alnSegs.size() - 1)->startRpos + query_seq_info_all.at(i)->query_alnSegs.at(query_seq_info_all.at(i)->query_alnSegs.size() - 1)->seglen);
				if(queryCluSig->qcSig_vec.size()>0){
					query_seq_info_all.at(i)->cluster_finished_flag = true;
					q_cluster_a.push_back(query_seq_info_all.at(i));
				}
				//delete queryCluSig;
				destroyQueryQcSig(queryCluSig);
			}
		}
	}else{
		for(i=0; i<query_seq_info_all.size(); i++){
			if(query_seq_info_all.at(i)->cluster_finished_flag == false){
				//choose seed query from a_cluster
				seed_info_a = chooseSeedClusterQuery(query_seq_info_all.at(i),q_cluster_a);
				seed_info_b = chooseSeedClusterQuery(query_seq_info_all.at(i), q_cluster_b);
				if(seed_info_a->span_intersection>0){
					score_ratio_a = computeScoreRatio(query_seq_info_all.at(i), q_cluster_a.at(seed_info_a->id), seed_info_a->startSpanPos, seed_info_a->endSpanPos);
					if(seed_info_b->span_intersection>0){
						score_ratio_b = computeScoreRatio(query_seq_info_all.at(i), q_cluster_b.at(seed_info_b->id), seed_info_b->startSpanPos, seed_info_b->endSpanPos);
						if(score_ratio_a>score_ratio_b){
							query_seq_info_all.at(i)->cluster_finished_flag = true;
							q_cluster_a.push_back(query_seq_info_all.at(i));
						}else{
							if(score_ratio_b==score_ratio_a){
								query_seq_info_all.at(i)->cluster_finished_flag = true;
								q_cluster_rescue.push_back(query_seq_info_all.at(i));
							}else{
								query_seq_info_all.at(i)->cluster_finished_flag = true;
								q_cluster_b.push_back(query_seq_info_all.at(i));
							}
						}
					}else{
						if(score_ratio_a<0.99){
							query_seq_info_all.at(i)->cluster_finished_flag = true;
							q_cluster_b.push_back(query_seq_info_all.at(i));
						}else{
							query_seq_info_all.at(i)->cluster_finished_flag = true;
							q_cluster_a.push_back(query_seq_info_all.at(i));
						}
					}
					delete seed_info_b;
				}else{
					if(seed_info_b->span_intersection>0){
						score_ratio_b = computeScoreRatio(query_seq_info_all.at(i), q_cluster_b.at(seed_info_b->id), seed_info_b->startSpanPos, seed_info_b->endSpanPos);
						if(score_ratio_b<0.99){
							query_seq_info_all.at(i)->cluster_finished_flag = true;
							q_cluster_a.push_back(query_seq_info_all.at(i));
						}else{
							query_seq_info_all.at(i)->cluster_finished_flag = true;
							q_cluster_b.push_back(query_seq_info_all.at(i));
						}
					}
//					else{
//						//no need to deal with it
//						query_seq_info_all.at(i)->cluster_finished_flag = false;
//					}
				}
				delete seed_info_a;
			}
		}
	}

	if(q_cluster_rescue.size()>0){
		for(i=0;i<q_cluster_rescue.size();i++){
			seed_info_a = chooseSeedClusterQuery(q_cluster_rescue.at(i),q_cluster_a);
			seed_info_b = chooseSeedClusterQuery(q_cluster_rescue.at(i), q_cluster_b);
			score_ratio_a = computeScoreRatio(q_cluster_rescue.at(i), q_cluster_a.at(seed_info_a->id), seed_info_a->startSpanPos, seed_info_a->endSpanPos);
			score_ratio_b = computeScoreRatio(q_cluster_rescue.at(i), q_cluster_b.at(seed_info_b->id), seed_info_b->startSpanPos, seed_info_b->endSpanPos);
			if(score_ratio_a>score_ratio_b){
				q_cluster_a.push_back(q_cluster_rescue.at(i));
			}else{
				if(score_ratio_a!=score_ratio_b){
					q_cluster_b.push_back(q_cluster_rescue.at(i));
				}else{
					query_seq_info_all.at(i)->cluster_finished_flag = false;
				}
			}
		}
	}

	for(i=0; i<q_cluster_a.size(); i++){
		queryCluSig = new queryCluSig_t();
		queryCluSig->qcSig_vec = extractQcSigsFromAlnSegsSingleQuery(q_cluster_a.at(i), q_cluster_a.at(i)->query_alnSegs.at(0)->startRpos, q_cluster_a.at(i)->query_alnSegs.at(q_cluster_a.at(i)->query_alnSegs.size() - 1)->startRpos + q_cluster_a.at(i)->query_alnSegs.at(q_cluster_a.at(i)->query_alnSegs.size() - 1)->seglen);
		if(queryCluSig->qcSig_vec.size()>0){
			q_cluster_a_final.push_back(q_cluster_a.at(i));
		}
		//delete queryCluSig;
		destroyQueryQcSig(queryCluSig);
	}

	for(i=0; i<q_cluster_b.size(); i++){
		queryCluSig = new queryCluSig_t();
		queryCluSig->qcSig_vec = extractQcSigsFromAlnSegsSingleQuery(q_cluster_b.at(i), q_cluster_b.at(i)->query_alnSegs.at(0)->startRpos, q_cluster_b.at(i)->query_alnSegs.at(q_cluster_b.at(i)->query_alnSegs.size() - 1)->startRpos + q_cluster_b.at(i)->query_alnSegs.at(q_cluster_b.at(i)->query_alnSegs.size() - 1)->seglen);
		if(queryCluSig->qcSig_vec.size()>0){
			q_cluster_b_final.push_back(q_cluster_b.at(i));
		}
		//delete queryCluSig;
		destroyQueryQcSig(queryCluSig);
	}

	if(q_cluster_a_final.size()>0){
		q_cluster_node_a = new struct querySeqInfoVec();
		q_cluster_node_a->query_seq_info_all = q_cluster_a_final;
		query_clu_vec.push_back(q_cluster_node_a);
	}

	if(q_cluster_b_final.size()>0){
		q_cluster_node_b = new struct querySeqInfoVec();
		q_cluster_node_b->query_seq_info_all = q_cluster_b_final;
		query_clu_vec.push_back(q_cluster_node_b);
	}

	return query_clu_vec;
}

void LocalAssembly::resucueCluster(vector<struct querySeqInfoNode*> &query_seq_info_all, vector<struct querySeqInfoNode*> &q_cluster_a, vector<struct querySeqInfoNode*> &q_cluster_b){
	//initializing b_cluster: choose one that is most different query from a_cluster add in b_cluster
	struct seedQueryInfo *seed_info_a;
	struct seedQueryInfo *seed_info_b;
	double score_ratio_a, score_ratio_b, sco_min;
	size_t i, j, seed_a_id,seed_b_id, smin_id;

	sco_min = 1;
	for(i=0; i<query_seq_info_all.size(); i++){
		if(query_seq_info_all.at(i)->cluster_finished_flag == false){
			seed_info_a = chooseSeedClusterQuery(query_seq_info_all.at(i), q_cluster_a);
			if(seed_info_a->span_intersection>0){
				score_ratio_a = computeScoreRatio(query_seq_info_all.at(i), q_cluster_a.at(seed_info_a->id), seed_info_a->startSpanPos, seed_info_a->endSpanPos);
				if(score_ratio_a<sco_min){
					sco_min = score_ratio_a;
					smin_id = i;
				}
			}
			delete seed_info_a;
		}
	}
	if(sco_min < 1){
		query_seq_info_all.at(smin_id)->cluster_finished_flag = true;
		q_cluster_b.push_back(query_seq_info_all.at(smin_id));
	}

	//classify
	for(i=0; i<query_seq_info_all.size(); i++){
		if(query_seq_info_all.at(i)->cluster_finished_flag == false){
			//choose seed query from a_cluster
			seed_info_a = chooseSeedClusterQuery(query_seq_info_all.at(i),q_cluster_a);
			if(seed_info_a->span_intersection>0){
				seed_info_b = chooseSeedClusterQuery02(seed_info_a, q_cluster_b);
				if(seed_info_b->span_intersection>0){
					score_ratio_a = computeScoreRatio(query_seq_info_all.at(i), q_cluster_a.at(seed_info_a->id), seed_info_b->startSpanPos, seed_info_b->endSpanPos);
					score_ratio_b = computeScoreRatio(query_seq_info_all.at(i), q_cluster_b.at(seed_info_b->id), seed_info_b->startSpanPos, seed_info_b->endSpanPos);
					if(score_ratio_a>score_ratio_b){
						query_seq_info_all.at(i)->cluster_finished_flag = true;
						q_cluster_a.push_back(query_seq_info_all.at(i));
					}
					if(score_ratio_a<score_ratio_b){
						query_seq_info_all.at(i)->cluster_finished_flag = true;
						q_cluster_b.push_back(query_seq_info_all.at(i));
					}
				}
				delete seed_info_b;
			}
			delete seed_info_a;
		}
	}
}

double LocalAssembly::computeScoreRatio(struct querySeqInfoNode* query_seq_info_node, struct querySeqInfoNode* q_cluster_node, int64_t startSpanPos, int64_t endSpanPos){
	double score_ratio;
	int32_t score_sum;
	size_t i;
	queryCluSig_t *queryCluSig, *seed_qcQuery;
	score_ratio = 0;
	//if(query_seq_info_node->query_alnSegs.size()>0){
		//matching prepare
	queryCluSig = new queryCluSig_t();
	queryCluSig->qcSig_vec = extractQcSigsFromAlnSegsSingleQuery(query_seq_info_node, startSpanPos, endSpanPos);
	seed_qcQuery = new queryCluSig_t();
	seed_qcQuery->qcSig_vec = extractQcSigsFromAlnSegsSingleQuery(q_cluster_node, startSpanPos, endSpanPos);
	//matching
	queryCluSig->match_profile_vec = computeQcMatchProfileSingleQuery(queryCluSig, seed_qcQuery);
	if(queryCluSig->match_profile_vec.size()==0){
		//score_ratio = 2;//have no varsig
		score_ratio = 0.99;
	}else{
		score_sum = 0;
		for(i=0;i<queryCluSig->match_profile_vec.size();i++){
			score_sum += queryCluSig->match_profile_vec.at(i);
		}
		score_ratio = (double)score_sum / queryCluSig->match_profile_vec.size();
	}

//	delete queryCluSig;
//	delete seed_qcQuery;
	destroyQueryQcSig(queryCluSig);
	destroyQueryQcSig(seed_qcQuery);
	//}

	return score_ratio;
}

vector<int8_t> LocalAssembly::computeQcMatchProfileSingleQuery(queryCluSig_t *queryCluSig, queryCluSig_t *seed_qcQuery){
	vector<int8_t> match_profile_vec;
	int32_t rowsNum, colsNum, matchScore, mismatchScore, gapScore, gapOpenScore;
	int32_t scoreIJ, tmp_gapScore1, tmp_gapScore2, maxValue, path_val, op_len_tmp, maxValue_ul, maxValue_l, maxValue_u;
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
		for(j=1; j<colsNum; j++) scoreArr[j].path_val = 2;
		for(i=1; i<rowsNum; i++) scoreArr[i*colsNum].path_val = 1;
		// compute the scores of each element
			for(i=1; i<rowsNum; i++){
			for(j=1; j<colsNum; j++){
				matchFlag = isQcSigMatch(queryCluSig->qcSig_vec.at(i-1), seed_qcQuery->qcSig_vec.at(j-1), QC_SIZE_RATIO_MATCH_THRES);
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

				maxValue = INT_MIN;
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
						scoreArr[i*colsNum+j].ismissmatch = true;
					}else{
						scoreArr[i*colsNum+j].ismissmatch = false;
					}
				}
				scoreArr[i*colsNum+j].score = maxValue;
				scoreArr[i*colsNum+j].path_val = path_val;
			}
		}

		// compute signature match profile
		match_profile_vec = qComputeSigMatchProfile(scoreArr, rowsNum, colsNum, queryCluSig, seed_qcQuery);
		free(scoreArr);
	}

	return match_profile_vec;
}

vector<int8_t> LocalAssembly::qComputeSigMatchProfile(struct alnScoreNode *scoreArr, int32_t rowsNum, int32_t colsNum, queryCluSig_t *queryCluSig, queryCluSig_t *seed_qcQuery){
	vector<int8_t> match_profile_vec;
	int32_t i = rowsNum - 1, j = colsNum - 1, value;

//	if(i==0 and j==0){
//		match_profile_vec.push_back(1);
//	}else if(i!=0 and j==0){
//		match_profile_vec.push_back(0);
//	}else if(i==0 and j!=0){
//		match_profile_vec.push_back(0);
//	}
	if(i==0 or j==0){
		match_profile_vec.push_back(0);
	}
	while(i>0 or j>0){
		value = scoreArr[i*colsNum+j].path_val;
		if(value==0){ //from (i-1, j-1)
			if(scoreArr[i*colsNum+j].ismissmatch){
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

	return match_profile_vec;
}

bool LocalAssembly::isQcSigMatch(qcSig_t *qc_sig, qcSig_t *seed_qc_sig, double size_ratio_match_thres){
	bool match_flag = false, flag;
	double size_ratio;
	int64_t start_pos1, end_pos1, start_pos2, end_pos2, ref_distance;

	if(qc_sig and seed_qc_sig){
		if(qc_sig->chrname.compare(seed_qc_sig->chrname)==0){
		 // same cigar_op, ratio of signature size less than 0.7
			if(qc_sig->cigar_op==seed_qc_sig->cigar_op){
				if(qc_sig->cigar_op_len>0 and seed_qc_sig->cigar_op_len>0) {
					if(qc_sig->cigar_op_len <= seed_qc_sig->cigar_op_len) size_ratio = (double)qc_sig->cigar_op_len / seed_qc_sig->cigar_op_len;
					else size_ratio = (double)seed_qc_sig->cigar_op_len / qc_sig->cigar_op_len;
					if(size_ratio>=size_ratio_match_thres){
						ref_distance = abs(qc_sig->ref_pos - seed_qc_sig->ref_pos);
						if(ref_distance<100){
							match_flag = true;
						}else match_flag = false;
					}else match_flag = false;
				}else match_flag = false;
			}
		}
	}
	return match_flag;
}

struct seedQueryInfo* LocalAssembly::chooseSeedClusterQuery(struct querySeqInfoNode* query_seq_info_node, vector<struct querySeqInfoNode*> &q_cluster){
	size_t seed_id, j;
	int64_t startSpanPos_tmp, endSpanPos_tmp, startSpanPos, endSpanPos, span, span_max;
	span_max = 0, seed_id = 0;
	startSpanPos = 0;
	endSpanPos = 0;
	struct seedQueryInfo* seed_info;
	//if(query_seq_info_node->query_alnSegs.size()>0){
	for(j=0;j<q_cluster.size();j++){
		startSpanPos_tmp = max(q_cluster.at(j)->query_alnSegs.at(0)->startRpos,query_seq_info_node->query_alnSegs.at(0)->startRpos);
		endSpanPos_tmp = min(q_cluster.at(j)->query_alnSegs.at(q_cluster.at(j)->query_alnSegs.size() - 1)->startRpos + q_cluster.at(j)->query_alnSegs.at(q_cluster.at(j)->query_alnSegs.size() - 1)->seglen -1, query_seq_info_node->query_alnSegs.at(query_seq_info_node->query_alnSegs.size() - 1)->startRpos + query_seq_info_node->query_alnSegs.at(query_seq_info_node->query_alnSegs.size() - 1)->seglen - 1);
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

struct seedQueryInfo* LocalAssembly::chooseSeedClusterQuery02(struct seedQueryInfo* seedinfo, vector<struct querySeqInfoNode*> &q_cluster){
	size_t seed_id, j;
	int64_t startSpanPos_tmp, endSpanPos_tmp, startSpanPos, endSpanPos, span, span_max, query_start, query_end;
	span_max = 0, seed_id = 0;
	startSpanPos = 0;
	endSpanPos = 0;
	struct seedQueryInfo* seed_info;
	//if(query_seq_info_node->query_alnSegs.size()>0){
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
		//}
	}

	seed_info = new struct seedQueryInfo();
	seed_info->startSpanPos = startSpanPos;
	seed_info->endSpanPos = endSpanPos;
	seed_info->span_intersection = span_max;
	seed_info->id = seed_id;
	return seed_info;
}

vector<qcSig_t*> LocalAssembly::extractQcSigsFromAlnSegsSingleQuery(struct querySeqInfoNode *query_seq_info_node, int64_t startSpanPos, int64_t endSpanPos){
	vector<qcSig_t*> sig_vec;
	qcSig_t *qc_sig;
	vector<struct alnSeg*>::iterator seg;
	size_t i;
	//int64_t startSpanPos_extend, endSpanPos_extend;

	//startSpanPos_extend = startSpanPos - 10;
	//endSpanPos_extend = endSpanPos + 10;

	//bool corrupt = false;
	for(seg=query_seq_info_node->query_alnSegs.begin(); seg!=query_seq_info_node->query_alnSegs.end(); seg++){
		switch((*seg)->opflag){
			case BAM_CINS:  // insertion in query
				if((*seg)->startRpos>=startSpanPos and (*seg)->startRpos<=endSpanPos){
					if((*seg)->seglen>=min_sv_size){
						qc_sig = allocateQcSigNode(*seg);
						qc_sig->chrname = varVec[0]->chrname;
						sig_vec.push_back(qc_sig);
					}
				}
				break;
			case BAM_CDEL:  // deletion in query
				if((*seg)->startRpos>=startSpanPos and ((*seg)->startRpos + (*seg)->seglen -1)<=endSpanPos){
					if((*seg)->seglen>=min_sv_size){
						qc_sig = allocateQcSigNode(*seg);
						qc_sig->chrname = varVec[0]->chrname;
						sig_vec.push_back(qc_sig);
					}
				}
				break;
			case BAM_CEQUAL:
				break;
			case BAM_CDIFF:
				break;
			case BAM_CSOFT_CLIP:// unexpected events
			case BAM_CHARD_CLIP:
			case BAM_CREF_SKIP:
			case BAM_CMATCH:
				break;
			case BAM_CPAD:
			case BAM_CBACK:
				cerr << __func__ << ", line=" << __LINE__ << ": invalid opflag " << (*seg)->opflag << endl;
				exit(1);
				break;
			default:
				cerr << __func__ << ", line=" << __LINE__ << ": invalid opflag " << (*seg)->opflag << endl;
//				cerr <<"op_len:" << (*seg)->seglen << endl;
//				cerr <<"startRpos:" << (*seg)->startRpos << endl;
//				cerr <<"qname:"<<query_seq_info_node->qname<<endl;
//				corrupt = true;
				exit(1);
				break;
		}
//		if(corrupt){
//			break;
//		}
	}

	return sig_vec;
}

// allocate query cluster signature
qcSig_t *LocalAssembly::allocateQcSigNode(struct alnSeg *aln_seg){
	qcSig_t *qc_sig = new qcSig_t();
	qc_sig->cigar_op = aln_seg->opflag;
	qc_sig->cigar_op_len = aln_seg->seglen;
	qc_sig->reg_contain_flag = false;
	qc_sig->ref_pos = aln_seg->startRpos;
	qc_sig->sig_id = -1;
	qc_sig->chrname = "";
	qc_sig->mate_qcSig = NULL;
	return qc_sig;
}

// destroy cluster signatures for queries
//void LocalAssembly::destroyQueryCluSigVec(vector<queryCluSig_t*> &queryCluSig_vec){
//	queryCluSig_t *queryCluSig;
//	vector<qcSig_t*> sig_vec_singleQuery;
//	qcSig_t *qc_sig;
//	for(size_t i=0; i<queryCluSig_vec.size(); i++){
//		queryCluSig = queryCluSig_vec.at(i);
//		for(size_t j=0; j<queryCluSig->qcSig_vec.size(); j++){
//			qc_sig = queryCluSig->qcSig_vec.at(j);
//			delete qc_sig;
//		}
//		vector<qcSig_t*>().swap(queryCluSig->qcSig_vec);
//		delete queryCluSig;
//	}
//	vector<queryCluSig_t*>().swap(queryCluSig_vec);
//}

struct seqsVec *LocalAssembly::smoothQuerySeqData(string &refseq, vector<struct querySeqInfoNode*> &query_seq_info_vec){
	size_t i, j;
	uint32_t op;
	struct seqsVec *smmothed_seqs_info;
	struct querySeqInfoNode *q_node;
	uint32_t query_startQpos;
	char refbase;
	string seq ,delete_seq, qname;

	smmothed_seqs_info = new struct seqsVec();

	for(i=0;i<query_seq_info_vec.size();i++){
		q_node = query_seq_info_vec.at(i);
		query_startQpos = 1;
		seq = q_node->seq;
		qname = q_node->qname;
		for(j=0;j<q_node->query_alnSegs.size();j++){
			op = q_node->query_alnSegs.at(j)->opflag;
			//query_endQpos += q_node->query_alnSegs.at(j)->seglen;
			switch (op) {
				case BAM_CMATCH:
					query_startQpos += q_node->query_alnSegs.at(j)->seglen;
					break;
				case BAM_CINS:
					if(q_node->query_alnSegs.at(j)->seglen<min_sv_size){
						seq.erase(query_startQpos - 1, q_node->query_alnSegs.at(j)->seglen);
					}else{
						query_startQpos += q_node->query_alnSegs.at(j)->seglen;
					}
					break;
				case BAM_CDEL:
					if(q_node->query_alnSegs.at(j)->seglen<min_sv_size){
						delete_seq = refseq.substr(q_node->query_alnSegs.at(j)->startRpos - startRefPos_assembly, q_node->query_alnSegs.at(j)->seglen);
						seq.insert(query_startQpos - 1, delete_seq);
						query_startQpos += q_node->query_alnSegs.at(j)->seglen;
					}
					break;
				case BAM_CSOFT_CLIP:
				case BAM_CHARD_CLIP:
					break;
				case BAM_CEQUAL:
					query_startQpos += q_node->query_alnSegs.at(j)->seglen;
					break;
				case BAM_CDIFF:
					refbase = refseq.at(q_node->query_alnSegs.at(j)->startRpos - startRefPos_assembly);
					seq[query_startQpos - 1] = refbase;
					query_startQpos += q_node->query_alnSegs.at(j)->seglen;
					break;
				case BAM_CREF_SKIP:
					// unexpected events
				case BAM_CPAD:
				case BAM_CBACK:
				default:
					cerr << __func__ << ", line=" << __LINE__ << ": invalid opflag " << op << endl;
					cerr << "startRpos" << q_node->query_alnSegs.at(j)->startRpos << endl;
					exit(1);
			}
		}
		smmothed_seqs_info->seqs.push_back(seq);
		smmothed_seqs_info->qname.push_back(qname);
	}
	return smmothed_seqs_info;
}

double LocalAssembly::computeCompensationCoefficient(size_t startRefPos_assembly, size_t endRefPos_assembly, double mean_read_len){
	size_t total_reg_size, refined_reg_size, reg_size_assemble;
	double comp_coefficient;

	reg_size_assemble = endRefPos_assembly - startRefPos_assembly + 1;
	total_reg_size = reg_size_assemble + 2 * mean_read_len;
	refined_reg_size = reg_size_assemble + mean_read_len;  // flanking_area = (2 * mean_read_len) / 2
	comp_coefficient = (double)total_reg_size / refined_reg_size;

	return comp_coefficient;
}

double LocalAssembly::computeLocalCov(vector<struct fqSeqNode*> &fq_seq_full_vec, double compensation_coefficient){
	double cov = 0, total;
	struct fqSeqNode* fq_node;
	size_t refined_reg_size;

	refined_reg_size = endRefPos_assembly - startRefPos_assembly + 1 + 2 * mean_read_len;
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

void LocalAssembly::samplingReads(vector<struct fqSeqNode*> &fq_seq_vec, double expect_cov_val, double compensation_coefficient){
	local_cov_original = computeLocalCov(fq_seq_vec, compensation_coefficient);

	if(local_cov_original>expect_cov_val){ // sampling
		//cout << "sampling for " << readsfilename << ", original coverage: " << local_cov_original << ", expected coverage: " << expect_cov_val << ", compensation_coefficient: " << compensation_coefficient << endl;
		samplingReadsOp(fq_seq_vec, expect_cov_val);
	}
}

void LocalAssembly::samplingReadsOp(vector<struct fqSeqNode*> &fq_seq_vec, double expect_cov_val){
	double expected_total_bases, total_bases;
	size_t index, max_reads_num, reg_size;
	struct fqSeqNode* fq_node;

	// reverse the select flag
	for(size_t i=0; i<fq_seq_vec.size(); i++){
		fq_node = fq_seq_vec.at(i);
		fq_node->selected_flag = false;
	}

	reg_size = endRefPos_assembly - startRefPos_assembly + 1 + mean_read_len;
	expected_total_bases = reg_size * expect_cov_val;
	max_reads_num = fq_seq_vec.size();

	pthread_mutex_lock(&mutex_down_sample);
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
	pthread_mutex_unlock(&mutex_down_sample);
	sampled_cov = total_bases / reg_size;
	sampling_flag = true;

//	/cout << "After sampling, reads count: " << reads_count << ", total bases: " << total_bases << endl;
}

void LocalAssembly::saveSampledReads(string &readsfilename, vector<struct fqSeqNode*> &fq_seq_vec){
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

	//cout <<"\t" << readsfilename << ": clip_aln_data_size=" << clipAlnDataVector.size() << ", reads_num=" << fq_seq_vec.size() << ", selected_num=" << selected_num << "; ref_size=" << endRefPos_assembly-startRefPos_assembly+1 << ", total_bases_original=" << total_bases_original << ", local_cov_original=" << local_cov_original << ", sampled_cov=" << sampled_cov << endl;
}

// get the non hard-clip align item
int32_t LocalAssembly::getNoHardClipAlnItem(vector<clipAlnData_t*> &query_aln_segs){
	int32_t idx;
	clipAlnData_t *clip_aln;

	idx = -1;
	for(size_t i=0; i<query_aln_segs.size(); i++){
		clip_aln = query_aln_segs.at(i);
		if(clip_aln->leftHardClippedFlag==false and clip_aln->rightHardClippedFlag==false){
			idx = i;
			break;
		}
	}

	return idx;
}

// mark Hard clipped segments from vector
void LocalAssembly::markHardClipSegs(size_t idx, vector<clipAlnData_t*> &query_aln_segs){
	for(size_t i=0; i<query_aln_segs.size(); ){
		if(i!=idx) {
			query_aln_segs.at(i)->query_checked_flag = true;
			query_aln_segs.erase(query_aln_segs.begin()+i);
		}else i++;
	}
}

// get query with soft clipped sequence
vector<string> LocalAssembly::getQuerySeqWithSoftClipSeqs(clipAlnData_t* clip_aln, string &refseq, int32_t bam_type){
	vector<string> query_info_vec;
	uint8_t *seq_int, *qual_int;
	string qseq, qual;
	int32_t i, query_start_loc, query_end_loc; // 0-based
	vector<int32_t> query_loc_vec;

	seq_int = bam_get_seq(clip_aln->bam);
	qual_int = bam_get_qual(clip_aln->bam);

	qseq = qual = "";
	if(clip_reg_flag==false){
		//if(seq_int[0]!=255){
			//queryStartEndLoc_vec = computeQueryStartEndLocByRefPos(clipAlnData_t *clip_aln, startRefPos_assembly, endRefPos_assembly);
		query_loc_vec = computeQueryStartEndLocByRefPos(clip_aln->bam, startRefPos_assembly, endRefPos_assembly, refseq, bam_type);
		if(query_loc_vec.size()>0){
			query_start_loc = query_loc_vec.at(0);
			query_end_loc = query_loc_vec.at(1);

	//			for(i=0; i<clip_aln->bam->core.l_qseq; i++) qseq += "=ACMGRSVTWYHKDBN"[bam_seqi(seq_int, i)];  // seq
	//
	//			if(qual_int[0]!=255) for(i=0; i<clip_aln->bam->core.l_qseq; i++) qual += qual_int[i] + 33;  // qual
	//			else for(i=0; i<clip_aln->bam->core.l_qseq; i++) qual += 34;  // qual

			for(i=query_start_loc - 1; i<query_end_loc; i++) qseq += "=ACMGRSVTWYHKDBN"[bam_seqi(seq_int, i)];  // seq

			if(qual_int[0]!=255) for(i=query_start_loc - 1; i<query_end_loc; i++) qual += qual_int[i] + 33;  // qual
			else for(i=query_start_loc - 1; i<query_end_loc; i++) qual += 34;  // qual

			query_info_vec.push_back(qseq);
			query_info_vec.push_back(qual);
		}
	}else{
//		for(i=0; i<sizeof(seq_int)/sizeof(uint8_t); i++) qseq += "=ACMGRSVTWYHKDBN"[bam_seqi(seq_int, i)];
//		if(qual_int[0]!=255) for(i=0; i<sizeof(qual_int)/sizeof(uint8_t); i++) qual += qual_int[i] + 33;  // qual
//		else for(i=0; i<sizeof(qual_int)/sizeof(uint8_t); i++) qual += 34;  // qual

		for(i=0; i<clip_aln->bam->core.l_qseq; i++) qseq += "=ACMGRSVTWYHKDBN"[bam_seqi(seq_int, i)];  // seq

		if(qual_int[0]!=255) for(i=0; i<clip_aln->bam->core.l_qseq; i++) qual += qual_int[i] + 33;  // qual
		else for(i=0; i<clip_aln->bam->core.l_qseq; i++) qual += 34;  // qual

		query_info_vec.push_back(qseq);
		query_info_vec.push_back(qual);
	}

	//}

	return query_info_vec;
}
void LocalAssembly::removeHeadTailAlnSegs(vector<struct alnSeg*> &alnSegs){
	uint32_t op;
	op = alnSegs.at(0)->opflag;
	switch (op) {
		case BAM_CMATCH:
			break;
		case BAM_CINS:
			delete *alnSegs.begin();
			alnSegs.erase(alnSegs.begin());
			break;
		case BAM_CDEL:
			delete *alnSegs.begin();
			alnSegs.erase(alnSegs.begin());
			break;
		case BAM_CSOFT_CLIP:
		case BAM_CHARD_CLIP:
			delete *alnSegs.begin();
			alnSegs.erase(alnSegs.begin());
			break;
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
	if(alnSegs.size()>0){
		op = alnSegs.at(alnSegs.size() - 1)->opflag;
		switch(op){
			case BAM_CMATCH:
				break;
			case BAM_CINS:
				delete *(alnSegs.end() - 1);
				alnSegs.erase(alnSegs.end() - 1);
				break;
			case BAM_CDEL:
				delete *(alnSegs.end() - 1);
				alnSegs.erase(alnSegs.end() - 1);
				break;
			case BAM_CSOFT_CLIP:
			case BAM_CHARD_CLIP:
				delete *(alnSegs.end() - 1);
				alnSegs.erase(alnSegs.end() - 1);
				break;
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

// compute query locations according to given reference positions
vector<int32_t> LocalAssembly::computeQueryStartEndLocByRefPos(bam1_t *b, int64_t startRefPos, int64_t endRefPos, string &refseq, int32_t bam_type){
	vector<int32_t> query_loc_vec;
	uint32_t op;
	int32_t target_startQpos = 0, target_endQpos = 0;
	//vector<struct alnSeg*> alnSegs;
	vector<struct alnSeg*> query_alnSegs;

	if(!(b->core.flag & BAM_FUNMAP)){ // aligned

//		string qname = bam_get_qname(b);
//		cout << "qname=" << qname << endl;

		switch(bam_type){
			case BAM_CIGAR_NO_DIFF_MD:
				//alnSegs = generateAlnSegs(b);
				query_alnSegs = generateAlnSegs2(b, startRefPos, endRefPos);
				break;
			case BAM_CIGAR_NO_DIFF_NO_MD:
			case BAM_CIGAR_DIFF_MD:
			case BAM_CIGAR_DIFF_NO_MD:
				//alnSegs = generateAlnSegs_no_MD(b, refseq, startPos, endPos);
				query_alnSegs = generateAlnSegs_no_MD2(b, refseq, startRefPos, endRefPos);
				break;
			default:
				cerr << __func__ << ": unknown bam type, error!" << endl;
				exit(1);
		}// generate align segments
		//outputAlnSegs(alnSegs);

		op = query_alnSegs.at(0)->opflag;
		if(query_alnSegs.size()>0){
			switch (op) {
				case BAM_CMATCH:
					target_startQpos = query_alnSegs.at(0)->startQpos;
					break;
				case BAM_CINS:
					target_startQpos = query_alnSegs.at(0)->startQpos + query_alnSegs.at(0)->seglen;
					delete *query_alnSegs.begin();
					query_alnSegs.erase(query_alnSegs.begin());
					break;
				case BAM_CDEL:
					target_startQpos = query_alnSegs.at(0)->startQpos;
					delete *query_alnSegs.begin();
					query_alnSegs.erase(query_alnSegs.begin());
					break;
				case BAM_CSOFT_CLIP:
				case BAM_CHARD_CLIP:
					target_startQpos = query_alnSegs.at(0)->startQpos + query_alnSegs.at(0)->seglen;
					delete *query_alnSegs.begin();
					query_alnSegs.erase(query_alnSegs.begin());
					break;
				case BAM_CEQUAL:
				case BAM_CDIFF:
					target_startQpos = query_alnSegs.at(0)->startQpos;
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

		if(query_alnSegs.size()>0){
			op = query_alnSegs.at(query_alnSegs.size() - 1)->opflag;
			switch(op){
				case BAM_CMATCH:
					target_endQpos = query_alnSegs.at(query_alnSegs.size() - 1)->startQpos +  query_alnSegs.at(query_alnSegs.size() - 1)->seglen - 1;
					break;
				case BAM_CINS:
					target_endQpos = query_alnSegs.at(query_alnSegs.size() - 1)->startQpos - 1;
					delete *(query_alnSegs.end() - 1);
					query_alnSegs.erase(query_alnSegs.end() - 1);
					break;
				case BAM_CDEL:
					target_endQpos = query_alnSegs.at(query_alnSegs.size() - 1)->startQpos - 1;
					delete *(query_alnSegs.end() - 1);
					query_alnSegs.erase(query_alnSegs.end() - 1);
					break;
				case BAM_CSOFT_CLIP:
				case BAM_CHARD_CLIP:
					target_endQpos = query_alnSegs.at(query_alnSegs.size() - 1)->startQpos - 1;
					delete *(query_alnSegs.end() - 1);
					query_alnSegs.erase(query_alnSegs.end() - 1);
					break;
				case BAM_CEQUAL:
				case BAM_CDIFF:
					target_endQpos = query_alnSegs.at(query_alnSegs.size() - 1)->startQpos +  query_alnSegs.at(query_alnSegs.size() - 1)->seglen - 1;
					break;
				case BAM_CREF_SKIP:
					// unexpected events
				case BAM_CPAD:
				case BAM_CBACK:
				default:
					cerr << __func__ << ", line=" << __LINE__ << ": invalid opflag " << op << endl;
					exit(1);
			}

			if(query_alnSegs.size()>0){
				destroyAlnSegs(query_alnSegs);
			}
			// save query locations
			query_loc_vec.push_back(target_startQpos);
			query_loc_vec.push_back(target_endQpos);
		}
	}

	return query_loc_vec;
}

// get query without soft clipped sequence
//vector<string> LocalAssembly::getSeqWithoutSoftClipSeqs(clipAlnData_t* clip_aln){
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
vector<string> LocalAssembly::joinQueryAlnSegs(vector<clipAlnData_t*> &query_aln_segs){
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

int32_t LocalAssembly::getLeftMostAlnSeg(vector<clipAlnData_t*> &query_aln_segs){
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

//// get heter local consensus using abPOA
//bool LocalAssembly::cnsByPoa(){
//	bool flag;
//	string poa_cmd, tmp_cons_filename, cmd1, cmd4;
//
//	// check the file
//	flag = isFileExist(contigfilename);
//	if(flag) return true; // cons was generated successfully previously
//
//	cmd1 = "mkdir -p " + tmpdir;
//	cmd4 = "rm -rf " + tmpdir;  // delete files
//
//	tmp_cons_filename = tmpdir + "/cons.fa";
//	system(cmd1.c_str());
//
//	poa_cmd = "abpoa -d " + readsfilename + " -o " + tmp_cons_filename + " > /dev/null 2>&1";
//	system(poa_cmd.c_str());
//
//	flag = isFileExist(tmp_cons_filename);
//	if(flag){ // consgenerated successfully
//		assem_success_flag = true;
//		rename(tmp_cons_filename.c_str(), contigfilename.c_str());
//		system(cmd4.c_str());
//	}else { // poa generated failed
//		system(cmd4.c_str());  // remove temporary files
//	}
//	return flag;
//}

// get heter local consensus using abPOA
bool LocalAssembly::cnsByPoa(){
	bool flag;
	string tmp_cons_filename, cmd1, cmd4;
	ofstream outfile;
	size_t k, serial_number;

	// check the file
	flag = isFileExist(contigfilename);
	if(flag) return true; // cons was generated successfully previously

	cmd1 = "mkdir -p " + tmpdir;
	cmd4 = "rm -rf " + tmpdir;  // delete files
	tmp_cons_filename = tmpdir + "/cons.fa";

	system(cmd1.c_str());

// save the cons to file
	outfile.open(tmp_cons_filename);
	if(!outfile.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file " << tmp_cons_filename << endl;
		exit(1);
	}

	for(k=0;k<seqs_vec.size();k++){
		 // --- POA
		uint n_seqs = seqs_vec.at(k)->seqs.size();
		abpoa_t *ab = abpoa_init();
		abpoa_para_t *abpt = abpoa_init_para();
		abpt->disable_seeding = 1;
		abpt->align_mode = 0; // global
		abpt->out_msa = 0;
		abpt->out_cons = 1;
		abpt->out_gfa = 0;
		// abpt->w = 6, abpt->k = 9;
		// abpt->min_w = 10; // minimizer-based seeding and partition
		// abpt->is_diploid = 1;
		// abpt->min_freq = 0.3;
		abpt->progressive_poa = 0;
		abpoa_post_set_para(abpt);

		int *seq_lens = (int *)malloc(sizeof(int) * n_seqs);
		uint8_t **bseqs = (uint8_t **)malloc(sizeof(uint8_t *) * n_seqs);
		for (uint i = 0; i < n_seqs; ++i) {
			seq_lens[i] = seqs_vec.at(k)->seqs[i].size();
			bseqs[i] = (uint8_t *)malloc(sizeof(uint8_t) * seq_lens[i]);
			for (int j = 0; j < seq_lens[i]; ++j) {
				bseqs[i][j] = nt4_table[(int)seqs_vec.at(k)->seqs[i][j]];
			}
		}

		abpoa_msa(ab, abpt, n_seqs, NULL, seq_lens, bseqs, NULL, NULL);

		abpoa_cons_t *abc = ab->abc;

		string cons = "";
		if (abc->n_cons > 0) {
			for (int j = 0; j < abc->cons_len[0]; ++j) {
				cons += "ACGTN"[abc->cons_base[0][j]];
			}
		}

		serial_number = k + 1;
		outfile << ">Consensus_sequence_" ;
		outfile << serial_number << "_"; //serial number
		outfile << seqs_vec.at(k)->seqs.size() << "-";  //cluster_number

		for(uint i=0; i<seqs_vec.at(k)->qname.size()-1; i++){
			outfile << seqs_vec.at(k)->qname.at(i) << "-";
		}
		outfile << seqs_vec.at(k)->qname.at(seqs_vec.at(k)->qname.size()-1) << endl;

		outfile << cons << endl;  // seq


		for (uint i = 0; i < n_seqs; ++i){
			free(bseqs[i]);
		}
		free(bseqs);
		free(seq_lens);
		abpoa_free(ab);
		abpoa_free_para(abpt);
	}

	outfile.close();
	destorySeqsVec(seqs_vec);

	flag = isFileExist(tmp_cons_filename);
	if(flag){ // consgenerated successfully
		assem_success_flag = true;
		rename(tmp_cons_filename.c_str(), contigfilename.c_str());
		system(cmd4.c_str());
	}else { // poa generated failed
		system(cmd4.c_str());  // remove temporary files
	}

	return flag;
}

// local assembly using Canu with more than one time
bool LocalAssembly::localAssembleCanu(){
	bool flag = localAssembleCanu_IncreaseGenomeSize();
	if(!flag) flag = localAssembleCanu_DecreaseGenomeSize();
	//if(!flag) cout << "ASS_FAILURE: " << contigfilename << endl;
	return flag;
}

// local assembly using Canu with increasing genome size parameter
bool LocalAssembly::localAssembleCanu_IncreaseGenomeSize(){
	string canu_cmd, assem_prefix, cmd2, cmd3, cmd4, tmp_ctg_filename, gnuplotTested_str, min_inout_cov_str, technology_str; // fast_option;
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
	assem_prefix = "assembly";
	tmp_ctg_filename = tmpdir + "/" + assem_prefix + ".contigs.fasta";
	cmd2 = "mv " + tmp_ctg_filename + " " + contigfilename;   // move and rename file
	cmd4 = "rm -rf " + tmpdir;  // delete files

	// limited number threads
	cmd_limited_threads_str = "";
	if(num_threads_per_assem_work>0){
		limited_threads_str = to_string(num_threads_per_assem_work);
		cmd_limited_threads_str = " maxThreads=" + limited_threads_str + " obtovlThreads=" + limited_threads_str + " corThreads=" + limited_threads_str + " utgovlThreads=" + limited_threads_str + " redThreads=" + limited_threads_str + " batThreads=" + limited_threads_str + " ";
	}

	// technology
	technology_str = "";
	if(technology.compare(PACBIO_CLR_TECH_STR)==0){
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
			cerr << __func__ << ", line=" << __LINE__ << ": installed canu version is " << canu_version << ", however, it is required to use the canu version " << MIN_CANU_VERSION_HIFI << " or higher to assemble Pacbio CCS sequencing data." << endl;
			exit(1);
		}
	}else{
		cerr << __func__ << ", line=" << __LINE__ << ": invalid sequencing technology: " << technology << endl;
		exit(1);
	}

	// increase the genome size
	genomeSize_Canu = ASSEMBLE_GENOME_SIZE_INITIAL;
	step_size = ASSEMBLE_STEP_SIZE;
	flag = false;
	for(i=1; i<=3; i++){
		// try canu1.7
		////canu_cmd = "canu -p " + assem_prefix + " -d " + tmpdir + " genomeSize=" + to_string(genomeSize_Canu) + " gnuplotTested=true -pacbio-raw " + readsfilename + " > /dev/null 2>&1";
		//canu_cmd = "canu -p " + assem_prefix + " -d " + tmpdir + " genomeSize=" + to_string(genomeSize_Canu) + gnuplotTested_str + " -pacbio-raw " + readsfilename + " > /dev/null 2>&1";
		////canu_cmd = "canu -p " + assem_prefix + " -d " + tmpdir + " genomeSize=" + to_string(genomeSize_Canu) + fast_option + gnuplotTested_str + " -pacbio-raw " + readsfilename + " > /dev/null 2>&1";
//		canu_cmd = "canu1.7 -p " + assem_prefix + " -d " + tmpdir + cmd_limited_threads_str + " genomeSize=" + to_string(genomeSize_Canu) + " gnuplotTested=true -pacbio-raw " + readsfilename + " > /dev/null 2>&1";
//		system(canu_cmd.c_str());  // local assembly, invoke Canu command
//
//		// save assembly result and remove temporary files if successfully assembled
//		flag = isFileExist(tmp_ctg_filename);
//		if(flag){ // contig generated successfully
//			rename(tmp_ctg_filename.c_str(), contigfilename.c_str());
//			system(cmd4.c_str());  // remove temporary files
//			break;
//		}else { // contig generated failed
//			system(cmd4.c_str());  // remove temporary files

			time(&end_time);
			cost_min = difftime(end_time, start_time) / 60.0;
			if(cost_min<=MAX_ASSEMBLE_MINUTES){
				// try canu1.8, or 2.0
				canu_cmd = "canu -p " + assem_prefix + " -d " + tmpdir + cmd_limited_threads_str + " genomeSize=" + to_string(genomeSize_Canu) + gnuplotTested_str + min_inout_cov_str + technology_str + readsfilename + " > /dev/null 2>&1";
				//cout << canu_cmd << endl;
				system(canu_cmd.c_str());  // local assembly, invoke Canu command

				// save assembly result and remove temporary files if successfully assembled
				flag = isFileExist(tmp_ctg_filename);
				if(flag){ // contig generated successfully
					assem_success_flag = true;
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

// local assembly using Canu with decreasing genome size parameter
bool LocalAssembly::localAssembleCanu_DecreaseGenomeSize(){
	string canu_cmd, assem_prefix, cmd2, cmd3, cmd4, tmp_ctg_filename, gnuplotTested_str, min_inout_cov_str, technology_str; //, fast_option;
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
	assem_prefix = "assembly";
	tmp_ctg_filename = tmpdir + "/" + assem_prefix + ".contigs.fasta";
	cmd2 = "mv " + tmp_ctg_filename + " " + contigfilename;   // move and rename file
	cmd4 = "rm -rf " + tmpdir;  // delete files

	// limited number threads
	cmd_limited_threads_str = "";
	if(num_threads_per_assem_work>0){
		limited_threads_str = to_string(num_threads_per_assem_work);
		cmd_limited_threads_str = " maxThreads=" + limited_threads_str + " obtovlThreads=" + limited_threads_str + " corThreads=" + limited_threads_str + " utgovlThreads=" + limited_threads_str + " redThreads=" + limited_threads_str + " batThreads=" + limited_threads_str + " ";
	}

	// technology
	technology_str = "";
	if(technology.compare(PACBIO_CLR_TECH_STR)==0){
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
			cerr << __func__ << ", line=" << __LINE__ << ": installed canu version is " << canu_version << ", however, it is required to use the canu version " << MIN_CANU_VERSION_HIFI << " or higher to assemble Pacbio CCS sequencing data." << endl;
			exit(1);
		}
	}else{
		cerr << __func__ << ", line=" << __LINE__ << ": invalid sequencing technology: " << technology << endl;
		exit(1);
	}

	// decrease the genome size
	step_size = ASSEMBLE_STEP_SIZE;
	genomeSize_Canu = ASSEMBLE_GENOME_SIZE_INITIAL - step_size;
	flag = false;
	for(i=1; i<=3 and genomeSize_Canu>=step_size; i++){
		// try canu1.7
		////canu_cmd = "canu -p " + assem_prefix + " -d " + tmpdir + " genomeSize=" + to_string(genomeSize_Canu) + " gnuplotTested=true -pacbio-raw " + readsfilename + " > /dev/null 2>&1";
		//canu_cmd = "canu -p " + assem_prefix + " -d " + tmpdir + " genomeSize=" + to_string(genomeSize_Canu) + gnuplotTested_str + " -pacbio-raw " + readsfilename + " > /dev/null 2>&1";
		////canu_cmd = "canu -p " + assem_prefix + " -d " + tmpdir + " genomeSize=" + to_string(genomeSize_Canu) + fast_option + gnuplotTested_str + " -pacbio-raw " + readsfilename + " > /dev/null 2>&1";
//		canu_cmd = "canu1.7 -p " + assem_prefix + " -d " + tmpdir + cmd_limited_threads_str + " genomeSize=" + to_string(genomeSize_Canu) + " gnuplotTested=true -pacbio-raw " + readsfilename + " > /dev/null 2>&1";
//		system(canu_cmd.c_str());  // local assembly, invoke Canu command
//
//		// save assembly result and remove temporary files if successfully assembled
//		flag = isFileExist(tmp_ctg_filename);
//		if(flag){ // contig generated successfully
//			rename(tmp_ctg_filename.c_str(), contigfilename.c_str());
//			system(cmd4.c_str());  // remove temporary files
//			break;
//		}else { // contig generated failed
//			system(cmd4.c_str());  // remove temporary files

			time(&end_time);
			cost_min = difftime(end_time, start_time) / 60.0;
			if(cost_min<=MAX_ASSEMBLE_MINUTES){
				// try canu1.8, or 2.0
				canu_cmd = "canu -p " + assem_prefix + " -d " + tmpdir + cmd_limited_threads_str + " genomeSize=" + to_string(genomeSize_Canu) + gnuplotTested_str + min_inout_cov_str + technology_str + readsfilename + " > /dev/null 2>&1";
				//cout << canu_cmd << endl;
				system(canu_cmd.c_str());  // local assembly, invoke Canu command

				// save assembly result and remove temporary files if successfully assembled
				flag = isFileExist(tmp_ctg_filename);
				if(flag){ // contig generated successfully
					assem_success_flag = true;
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

// record assembly information
void LocalAssembly::recordAssemblyInfo(ofstream &assembly_info_file){
	string line, assembly_status, header, left_shift_size_str, right_shift_size_str, reg_str, sampling_str, limit_reg_str, limit_reg_str2;
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

	// assembly status
	if (isFileExist(contigfilename)) assembly_status = ASSEMBLY_SUCCESS;
	else assembly_status = ASSEMBLY_FAILURE;

	line = refseqfilename + "\t" + contigfilename + "\t" + readsfilename + "\t" + left_shift_size_str + "\t" + right_shift_size_str + "\t" + assembly_status;

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
	assembly_info_file << line << endl;
	pthread_mutex_unlock(&mutex_write);
}
