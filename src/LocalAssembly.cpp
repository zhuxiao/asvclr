#include "LocalAssembly.h"
#include "clipAlnDataLoader.h"

pthread_mutex_t mutex_write = PTHREAD_MUTEX_INITIALIZER;

LocalAssembly::LocalAssembly(string &readsfilename, string &contigfilename, string &refseqfilename, string &tmpdir, size_t num_threads_per_assem_work, vector<reg_t*> &varVec, string &chrname, string &inBamFile, faidx_t *fai, size_t assembly_extend_size, double expected_cov, bool delete_reads_flag){
	this->chrname = chrname;
	this->chrlen = faidx_seq_len(fai, chrname.c_str()); // get reference size
	this->readsfilename = preprocessPipeChar(readsfilename);
	this->contigfilename = preprocessPipeChar(contigfilename);
	this->refseqfilename = preprocessPipeChar(refseqfilename);
	this->tmpdir = preprocessPipeChar(tmpdir);
	this->num_threads_per_assem_work = num_threads_per_assem_work;
	this->varVec = varVec;
	this->fai = fai;
	this->inBamFile = inBamFile;
	this->assembly_extend_size = ASSEMBLY_SIDE_EXT_SIZE + assembly_extend_size;
	startRefPos_assembly = endRefPos_assembly = 0;
	mean_read_len = 0;
	//this->canu_version = canu_version;

	ref_seq_size = reads_count_original = total_bases_original = reads_count = total_bases = 0;
	local_cov_original = sampled_cov = 0;
	this->expected_cov = expected_cov;
	this->compensation_coefficient = 1;
	sampling_flag = false;
	this->delete_reads_flag = delete_reads_flag;

	if(isFileExist(contigfilename)) assem_success_preDone_flag = true;
	else assem_success_preDone_flag = false;

	limit_reg_process_flag = false;
}

LocalAssembly::~LocalAssembly() {
	if(delete_reads_flag) remove(readsfilename.c_str());	// delete the reads file to save disk space
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
	startRefPos = varVec[0]->startRefPos - REFSEQ_SIDE_LEN;
	if(startRefPos<1) startRefPos = 1;
	left_shift_size = varVec[0]->startRefPos - startRefPos;  // left shift size
	endRefPos = varVec[varVec.size()-1]->endRefPos + REFSEQ_SIDE_LEN;
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
	size_t i, j, seq_id;
	string qname, seq, qual;
	vector<clipAlnData_t*> query_aln_segs;
	vector<string> query_seq_qual_vec;
	vector<struct fqSeqNode*> fq_seq_vec; // [0]: query name; [1]: sequence; [2]: query name (optional); [3]: quality
	struct fqSeqNode *fq_node;
	int64_t noHardClipIdx;

	startRefPos_assembly = varVec[0]->startRefPos - assembly_extend_size;
	if(startRefPos_assembly<1) startRefPos_assembly = 1;
	endRefPos_assembly = varVec[varVec.size()-1]->endRefPos + assembly_extend_size;
	if(endRefPos_assembly>chrlen) endRefPos_assembly = chrlen;

	// load the aligned reads data
	clipAlnDataLoader data_loader(varVec[0]->chrname, startRefPos_assembly, endRefPos_assembly, inBamFile);
	data_loader.loadClipAlnData(clipAlnDataVector);

	//cout << "start_pos_assembly=" << startRefPos_assembly << ", end_pos_assembly=" << endRefPos_assembly << ", clipAlnDataVector.size=" << clipAlnDataVector.size() << endl;

	for(i=0; i<clipAlnDataVector.size(); i++) clipAlnDataVector.at(i)->query_checked_flag = false;

	// join query clip align segments
	seq_id = 0;
	total_bases_original = 0;
	for(i=0; i<clipAlnDataVector.size(); i++){
		if(clipAlnDataVector.at(i)->query_checked_flag==false){
			qname = clipAlnDataVector.at(i)->queryname;
			query_aln_segs = getQueryClipAlnSegs(qname, clipAlnDataVector);  // get query clip align segments

			noHardClipIdx = getNoHardClipAlnItem(query_aln_segs);
			if(noHardClipIdx!=-1){
				query_seq_qual_vec = getQuerySeqWithSoftClipSeqs(query_aln_segs.at(noHardClipIdx));
				markHardClipSegs(noHardClipIdx, query_aln_segs);

				if(query_seq_qual_vec.size()>0){
					seq = query_seq_qual_vec.at(0);
					qual = query_seq_qual_vec.at(1);

					fq_node = new struct fqSeqNode();
					fq_node->seq_id = seq_id++;
					fq_node->seq_name = qname;
					fq_node->seq = seq;
					fq_node->qual = qual;
					fq_node->selected_flag = true;
					fq_seq_vec.push_back(fq_node);

					total_bases_original += seq.size();
				}
			}else{
				//cout << "qname=" << qname << ", querylen=" << clipAlnDataVector.at(i)->querylen << endl;
				// join query align segments without 'SA' tags
//				query_seq_qual_vec = joinQueryAlnSegs(query_aln_segs);

				for(j=0; j<query_aln_segs.size(); j++){
					query_seq_qual_vec = getQuerySeqWithSoftClipSeqs(query_aln_segs.at(j));
					if(query_seq_qual_vec.size()>0){
						seq = query_seq_qual_vec.at(0);
						qual = query_seq_qual_vec.at(1);

						fq_node = new struct fqSeqNode();
						fq_node->seq_id = seq_id++;
						fq_node->seq_name = qname;
						fq_node->seq = seq;
						fq_node->qual = qual;
						fq_node->selected_flag = true;
						fq_seq_vec.push_back(fq_node);

						total_bases_original += seq.size();
					}
				}
			}

			for(j=0; j<query_aln_segs.size(); j++) query_aln_segs.at(j)->query_checked_flag = true;
		}
	}
	reads_count_original = fq_seq_vec.size();
	mean_read_len = (double)total_bases_original / fq_seq_vec.size();

	//cout << "mean_read_len=" << mean_read_len << endl;

	if(!clipAlnDataVector.empty()) destoryClipAlnData();

	// sampling to expected coverage
	if(expected_cov!=0){
		compensation_coefficient = computeCompensationCoefficient(startRefPos_assembly, endRefPos_assembly, mean_read_len);
		//cout << compensation_coefficient << endl;
		samplingReads(fq_seq_vec, expected_cov, compensation_coefficient);
	}

	// save sampled reads to file
	saveSampledReads(readsfilename, fq_seq_vec);

	// release memory
	if(!fq_seq_vec.empty()) destoryFqSeqs(fq_seq_vec);
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

// get query clip align segments
vector<clipAlnData_t*> LocalAssembly::getQueryClipAlnSegs(string &queryname, vector<clipAlnData_t*> &clipAlnDataVector){
	vector<clipAlnData_t*> query_aln_segs;
	for(size_t i=0; i<clipAlnDataVector.size(); i++)
		if(clipAlnDataVector.at(i)->query_checked_flag==false and clipAlnDataVector.at(i)->queryname==queryname)
			query_aln_segs.push_back(clipAlnDataVector.at(i));
	return query_aln_segs;
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
vector<string> LocalAssembly::getQuerySeqWithSoftClipSeqs(clipAlnData_t* clip_aln){
	vector<string> query_info_vec;
	uint8_t *seq_int, *qual_int;
	string qseq, qual;
	int32_t i;

	seq_int = bam_get_seq(clip_aln->bam);
	qual_int = bam_get_qual(clip_aln->bam);

	qseq = qual = "";
	for(i=0; i<clip_aln->bam->core.l_qseq; i++) qseq += "=ACMGRSVTWYHKDBN"[bam_seqi(seq_int, i)];  // seq
	for(i=0; i<clip_aln->bam->core.l_qseq; i++) qual += qual_int[i] + 33;  // qual

	query_info_vec.push_back(qseq);
	query_info_vec.push_back(qual);

	return query_info_vec;
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
	size_t i, join_orient, seq_len;
	uint8_t *seq_int, *qual_int;
	int32_t j, left_most_idx, overlap_size, gap_size, startQueryPos, endQueryPos;
	clipAlnData_t *clip_aln, *pre_clip_aln;

	pre_clip_aln = NULL; join_orient = ALN_PLUS_ORIENT;
	startQueryPos = endQueryPos = -1;
	qseq_result = qual_result = "";
	for(i=0; i<query_aln_segs.size(); i++){
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
	size_t minPos;
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

// local assembly using Canu with more than one time
bool LocalAssembly::localAssembleCanu(){
	bool flag = localAssembleCanu_IncreaseGenomeSize();
	if(!flag) flag = localAssembleCanu_DecreaseGenomeSize();
	//if(!flag) cout << "ASS_FAILURE: " << contigfilename << endl;
	return flag;
}

// local assembly using Canu with increasing genome size parameter
bool LocalAssembly::localAssembleCanu_IncreaseGenomeSize(){
	string canu_cmd, assem_prefix, cmd2, cmd3, cmd4, tmp_ctg_filename, gnuplotTested_str, fast_option;
	int i, genomeSize_Canu, step_size;
	bool flag;
	string cmd_limited_threads_str, limited_threads_str;

//	gnuplotTested_str = "";
//	if(canu_version.compare("1.7.1")==0 or canu_version.compare("1.7")==0 or canu_version.compare("1.6")==0)
//		gnuplotTested_str = " gnuplotTested=true";

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

	// increase the genome size
	genomeSize_Canu = ASSEMBLY_GENOME_SIZE_INITIAL;
	step_size = ASSEMBLY_STEP_SIZE;
	flag = false;
	for(i=1; i<=3; i++){
		// try canu1.7
		////canu_cmd = "canu -p " + assem_prefix + " -d " + tmpdir + " genomeSize=" + to_string(genomeSize_Canu) + " gnuplotTested=true -pacbio-raw " + readsfilename + " > /dev/null 2>&1";
		//canu_cmd = "canu -p " + assem_prefix + " -d " + tmpdir + " genomeSize=" + to_string(genomeSize_Canu) + gnuplotTested_str + " -pacbio-raw " + readsfilename + " > /dev/null 2>&1";
		////canu_cmd = "canu -p " + assem_prefix + " -d " + tmpdir + " genomeSize=" + to_string(genomeSize_Canu) + fast_option + gnuplotTested_str + " -pacbio-raw " + readsfilename + " > /dev/null 2>&1";
		canu_cmd = "canu1.7 -p " + assem_prefix + " -d " + tmpdir + cmd_limited_threads_str + " genomeSize=" + to_string(genomeSize_Canu) + " gnuplotTested=true -pacbio-raw " + readsfilename + " > /dev/null 2>&1";
		system(canu_cmd.c_str());  // local assembly, invoke Canu command

		// save assembly result and remove temporary files if successfully assembled
		flag = isFileExist(tmp_ctg_filename);
		if(flag){ // contig generated successfully
			rename(tmp_ctg_filename.c_str(), contigfilename.c_str());
			system(cmd4.c_str());  // remove temporary files
			break;
		}else { // contig generated failed
			system(cmd4.c_str());  // remove temporary files

			// try canu1.8, or 2.0
			canu_cmd = "canu2.0 -p " + assem_prefix + " -d " + tmpdir + cmd_limited_threads_str + " genomeSize=" + to_string(genomeSize_Canu) + " -pacbio-raw " + readsfilename + " > /dev/null 2>&1";
			//cout << canu_cmd << endl;
			system(canu_cmd.c_str());  // local assembly, invoke Canu command

			// save assembly result and remove temporary files if successfully assembled
			flag = isFileExist(tmp_ctg_filename);
			if(flag){ // contig generated successfully
				rename(tmp_ctg_filename.c_str(), contigfilename.c_str());
				system(cmd4.c_str());  // remove temporary files
				break;
			}else { // contig generated failed
				system(cmd4.c_str());  // remove temporary files
				genomeSize_Canu += i * step_size;
			}
		}
	}
	return flag;
}

// local assembly using Canu with decreasing genome size parameter
bool LocalAssembly::localAssembleCanu_DecreaseGenomeSize(){
	string canu_cmd, assem_prefix, cmd2, cmd3, cmd4, tmp_ctg_filename, gnuplotTested_str, fast_option;
	int i, genomeSize_Canu, step_size;
	bool flag;
	string cmd_limited_threads_str, limited_threads_str;

//	gnuplotTested_str = "";
//	if(canu_version.compare("1.7.1")==0 or canu_version.compare("1.7")==0 or canu_version.compare("1.6")==0)
//		gnuplotTested_str = " gnuplotTested=true";

//	fast_option = "";
//	if(canu_version.compare("1.8")==0)
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

	// decrease the genome size
	step_size = ASSEMBLY_STEP_SIZE;
	genomeSize_Canu = ASSEMBLY_GENOME_SIZE_INITIAL - step_size;
	flag = false;
	for(i=1; i<=3 and genomeSize_Canu>=step_size; i++){
		// try canu1.7
		////canu_cmd = "canu -p " + assem_prefix + " -d " + tmpdir + " genomeSize=" + to_string(genomeSize_Canu) + " gnuplotTested=true -pacbio-raw " + readsfilename + " > /dev/null 2>&1";
		//canu_cmd = "canu -p " + assem_prefix + " -d " + tmpdir + " genomeSize=" + to_string(genomeSize_Canu) + gnuplotTested_str + " -pacbio-raw " + readsfilename + " > /dev/null 2>&1";
		////canu_cmd = "canu -p " + assem_prefix + " -d " + tmpdir + " genomeSize=" + to_string(genomeSize_Canu) + fast_option + gnuplotTested_str + " -pacbio-raw " + readsfilename + " > /dev/null 2>&1";
		canu_cmd = "canu1.7 -p " + assem_prefix + " -d " + tmpdir + cmd_limited_threads_str + " genomeSize=" + to_string(genomeSize_Canu) + " gnuplotTested=true -pacbio-raw " + readsfilename + " > /dev/null 2>&1";
		system(canu_cmd.c_str());  // local assembly, invoke Canu command

		// save assembly result and remove temporary files if successfully assembled
		flag = isFileExist(tmp_ctg_filename);
		if(flag){ // contig generated successfully
			rename(tmp_ctg_filename.c_str(), contigfilename.c_str());
			system(cmd4.c_str());  // remove temporary files
			break;
		}else { // contig generated failed
			system(cmd4.c_str());  // remove temporary files

			// try canu1.8, or 2.0
			canu_cmd = "canu2.0 -p " + assem_prefix + " -d " + tmpdir + cmd_limited_threads_str + " genomeSize=" + to_string(genomeSize_Canu) + " -pacbio-raw " + readsfilename + " > /dev/null 2>&1";
			//cout << canu_cmd << endl;
			system(canu_cmd.c_str());  // local assembly, invoke Canu command

			// save assembly result and remove temporary files if successfully assembled
			flag = isFileExist(tmp_ctg_filename);
			if(flag){ // contig generated successfully
				rename(tmp_ctg_filename.c_str(), contigfilename.c_str());
				system(cmd4.c_str());  // remove temporary files
				break;
			}else { // contig generated failed
				system(cmd4.c_str());  // remove temporary files
				genomeSize_Canu -= i * step_size;
			}
		}
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
