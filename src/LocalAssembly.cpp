#include <sys/stat.h>

#include "LocalAssembly.h"

pthread_mutex_t mutex_write = PTHREAD_MUTEX_INITIALIZER;

LocalAssembly::LocalAssembly(string &readsfilename, string &contigfilename, string &refseqfilename, string &tmpdir, vector<reg_t*> &varVec, string &chrname, string &inBamFile, faidx_t *fai) {
	this->chrname = chrname;
	this->chrlen = faidx_seq_len(fai, chrname.c_str()); // get reference size
	this->readsfilename = readsfilename;
	this->contigfilename = contigfilename;
	this->refseqfilename = refseqfilename;
	this->tmpdir = tmpdir;
	this->varVec = varVec;
	this->fai = fai;
	this->inBamFile = inBamFile;
}

LocalAssembly::~LocalAssembly() {
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
	if(endRefPos>(int32_t)chrlen) endRefPos = chrlen;
	right_shift_size = endRefPos - varVec[varVec.size()-1]->endRefPos;  // right shift size

	// get the region sequence
	reg = chrname + ":" + to_string(startRefPos) + "-" + to_string(endRefPos);
	RefSeqLoader refseq_loader(reg, fai);
	refseq_loader.getRefSeq();

	// save the refseq to file
	outfile.open(refseqfilename);
	if(!outfile.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file " << refseqfilename << endl;
		exit(1);
	}

	outfile << ">" + reg + "|" + to_string(left_shift_size) + "|" + to_string(right_shift_size) + "|" + contigfilename << endl;  // header
	outfile << refseq_loader.refseq << endl;  // seq

	outfile.close();
}

// extract the reads data from BAM file
void LocalAssembly::extractReadsDataFromBAM(){
	size_t i, j;
	string qname, qseq, qual;
	ofstream outfile;
	vector<clipAlnData_t*> query_aln_segs;
	vector<string> query_seq_qual_vec;
	int32_t noHardClipIdx;

	// load the aligned reads data
	clipAlnDataLoader data_loader(varVec[0]->chrname, varVec[0]->startRefPos, varVec[varVec.size()-1]->endRefPos, inBamFile);
	data_loader.loadClipAlnData(clipAlnDataVector);

	// save to file
	outfile.open(readsfilename);
	if(!outfile.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file " << readsfilename << endl;
		exit(1);
	}

	for(i=0; i<clipAlnDataVector.size(); i++) clipAlnDataVector.at(i)->query_checked_flag = false;

	// join query clip align segments
	for(i=0; i<clipAlnDataVector.size(); i++){
		if(clipAlnDataVector.at(i)->query_checked_flag==false){
			qname = clipAlnDataVector.at(i)->queryname;
			query_aln_segs = getQueryClipAlnSegs(qname, clipAlnDataVector);  // get query clip align segments

			noHardClipIdx = getNoHardClipAlnItem(query_aln_segs);
			if(noHardClipIdx!=-1){
				query_seq_qual_vec = getQuerySeqWithSoftClipSeqs(query_aln_segs.at(noHardClipIdx));
				markHardClipSegs(noHardClipIdx, query_aln_segs);

				if(query_seq_qual_vec.size()){
					qseq = query_seq_qual_vec.at(0);
					qual = query_seq_qual_vec.at(1);

					// save the read to file
					outfile << "@" + qname << endl;
					outfile << qseq << endl;
					outfile << "+" << endl;
					outfile << qual << endl;
				}
			}else{
				//cout << "qname=" << qname << ", querylen=" << clipAlnDataVector.at(i)->querylen << endl;
				// join query align segments without 'SA' tags
//				query_seq_qual_vec = joinQueryAlnSegs(query_aln_segs);

				for(j=0; j<query_aln_segs.size(); j++){
					query_seq_qual_vec = getQuerySeqWithSoftClipSeqs(query_aln_segs.at(j));
					if(query_seq_qual_vec.size()){
						qseq = query_seq_qual_vec.at(0);
						qual = query_seq_qual_vec.at(1);

						// save the read to file
						outfile << "@" + qname << endl;
						outfile << qseq << endl;
						outfile << "+" << endl;
						outfile << qual << endl;
					}
				}
			}

			for(j=0; j<query_aln_segs.size(); j++) query_aln_segs.at(j)->query_checked_flag = true;
		}
	}
	outfile.close();

	if(!clipAlnDataVector.empty()) destoryClipAlnData();


//	size_t i, j;
//	string qname, qseq, qual;
//	uint8_t *qual_int;
//	ofstream outfile;
//	bam1_t *b;
//
//	// load the aligned reads data
//	alnDataLoader data_loader(varVec[0]->chrname, varVec[0]->startRefPos, varVec[varVec.size()-1]->endRefPos, inBamFile);
//	data_loader.loadAlnData(alnDataVector);
//
//	// save to file
//	outfile.open(readsfilename);
//	if(!outfile.is_open()){
//		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file " << readsfilename << endl;
//		exit(1);
//	}
//
//	// convert to fastq and save to file
//	for(i=0; i<(int32_t)alnDataVector.size(); i++){
//		b = alnDataVector[i];
//		qname = bam_get_qname(b);   // qname
//		qual_int = bam_get_qual(b);
//
//		qseq = qual = "";
//		for(j=0; j<b->core.l_qseq; j++) qseq += "=ACMGRSVTWYHKDBN"[bam_seqi(bam_get_seq(b), j)];  // seq
//		for(j=0; j<b->core.l_qseq; j++) qual += qual_int[j] + 33;  // qual
//
//		// save the read to file
//		outfile << "@" + qname << endl;
//		outfile << qseq << endl;
//		outfile << "+" << endl;
//		outfile << qual << endl;
//	}
//	outfile.close();
//
//	if(!alnDataVector.empty()) destoryAlnData();

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

	pre_clip_aln = NULL;
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
	if(!flag) cout << "ASS_FAILURE: " << contigfilename << endl;
	return flag;
}

// local assembly using Canu with increasing genome size parameter
bool LocalAssembly::localAssembleCanu_IncreaseGenomeSize(){
	string canu_cmd, assem_prefix, cmd2, cmd3, cmd4, tmp_ctg_filename;
	int i, genomeSize_Canu, step_size;
	bool flag;

	// generate command string
	assem_prefix = "assembly";
	tmp_ctg_filename = tmpdir + "/" + assem_prefix + ".contigs.fasta";
	cmd2 = "mv " + tmp_ctg_filename + " " + contigfilename;   // move and rename file
	cmd4 = "rm -rf " + tmpdir;  // delete files

	// increase the genome size
	genomeSize_Canu = ASSEMBLY_GENOME_SIZE_INITIAL;
	step_size = ASSEMBLY_STEP_SIZE;
	flag = false;
	for(i=1; i<=5; i++){
		canu_cmd = "canu -p " + assem_prefix + " -d " + tmpdir + " genomeSize=" + to_string(genomeSize_Canu) + " gnuplotTested=true -pacbio-raw " + readsfilename + " > /dev/null 2>&1";
		system(canu_cmd.c_str());  // local assembly, invoke Canu command

		// save assembly result and remove temporary files if successfully assembled
		struct stat fileStat;
		if (stat(tmp_ctg_filename.c_str(), &fileStat) == 0)
			if(fileStat.st_size>0) flag = true; // contig generated successfully

		if(flag){ // contig generated successfully
			rename(tmp_ctg_filename.c_str(), contigfilename.c_str());
			system(cmd4.c_str());  // remove temporary files
			break;
		}else { // contig generated failed
			genomeSize_Canu += i * step_size;
			system(cmd4.c_str());  // remove temporary files
		}
	}
	return flag;
}

// local assembly using Canu with decreasing genome size parameter
bool LocalAssembly::localAssembleCanu_DecreaseGenomeSize(){
	string canu_cmd, assem_prefix, cmd2, cmd3, cmd4, tmp_ctg_filename;
	int i, genomeSize_Canu, step_size;
	bool flag;

	// generate command string
	assem_prefix = "assembly";
	tmp_ctg_filename = tmpdir + "/" + assem_prefix + ".contigs.fasta";
	cmd2 = "mv " + tmp_ctg_filename + " " + contigfilename;   // move and rename file
	cmd4 = "rm -rf " + tmpdir;  // delete files

	// increase the genome size
	step_size = ASSEMBLY_STEP_SIZE;
	genomeSize_Canu = ASSEMBLY_GENOME_SIZE_INITIAL - step_size;
	flag = false;
	for(i=1; i<=5 and genomeSize_Canu>=step_size; i++){
		canu_cmd = "canu -p " + assem_prefix + " -d " + tmpdir + " genomeSize=" + to_string(genomeSize_Canu) + " gnuplotTested=true -pacbio-raw " + readsfilename + " > /dev/null 2>&1";
		system(canu_cmd.c_str());  // local assembly, invoke Canu command

		// save assembly result and remove temporary files if successfully assembled
		struct stat fileStat;
		if (stat(tmp_ctg_filename.c_str(), &fileStat) == 0)
			if(fileStat.st_size>0) flag = true; // contig generated successfully

		if(flag){ // contig generated successfully
			rename(tmp_ctg_filename.c_str(), contigfilename.c_str());
			system(cmd4.c_str());  // remove temporary files
			break;
		}else { // contig generated failed
			genomeSize_Canu -= i * step_size;
			system(cmd4.c_str());  // remove temporary files
		}
	}
	return flag;
}

// record assembly information
void LocalAssembly::recordAssemblyInfo(ofstream &assembly_info_file){
	string line, assembly_status, header, left_shift_size_str, right_shift_size_str;
	reg_t *reg;
	ifstream infile;
	vector<string> str_vec;

	// ref shift size
	infile.open(refseqfilename);
	if(!infile.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file " << refseqfilename << endl;
		exit(1);
	}

	getline(infile, header);
	str_vec = split(header, "|");
	left_shift_size_str = str_vec[1];
	right_shift_size_str = str_vec[2];

	infile.close();

	// assembly status
	struct stat fileStat;
	if (stat(contigfilename.c_str(), &fileStat) == 0)
		assembly_status = ASSEMBLY_SUCCESS;
	else
		assembly_status = ASSEMBLY_FAILURE;

	line = refseqfilename + "\t" + contigfilename + "\t" + readsfilename + "\t" + left_shift_size_str + "\t" + right_shift_size_str + "\t" + assembly_status;

	for(int32_t i=0; i<(int32_t)varVec.size(); i++){
		reg = varVec[i];
		line += "\t" + reg->chrname + ":" + to_string(reg->startRefPos) + "-" + to_string(reg->endRefPos);
	}

	pthread_mutex_lock(&mutex_write);
	assembly_info_file << line << endl;
	pthread_mutex_unlock(&mutex_write);
}
