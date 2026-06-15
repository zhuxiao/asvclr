#include "alnDataLoader.h"

alnDataLoader::alnDataLoader(string &chrname, int32_t startRefPos, int32_t endRefPos, string &inBamFile, int32_t minMapQ, int32_t minHighMapQ) {
	this->reg_str = chrname + ":" + to_string(startRefPos) + "-" + to_string(endRefPos);
	this->inBamFile = inBamFile;
	mean_read_len = 0;
	localCov = -1;
	this->minMapQ = minMapQ;
	this->minHighMapQ = minHighMapQ;
	this->startRefPos = startRefPos;
	this->endRefPos = endRefPos;
}

alnDataLoader::~alnDataLoader() {
}

void alnDataLoader::loadAlnData(vector<bam1_t*> &alnDataVector, double max_ultra_high_cov){
	samFile *in = NULL;
	bam_hdr_t *header;
	vector<string> qname_vec;
	vector<int32_t> qlen_vec;
	vector<int32_t> qual_vec;
	size_t total_len = 0, total_num = 0;

	if ((in = sam_open(inBamFile.c_str(), "r")) == 0) {
		cerr << __func__ << ": failed to open " << inBamFile << " for reading" << endl;
		exit(1);
	}

	if ((header = sam_hdr_read(in)) == 0) {
		cerr << __func__ << ": fail to read the header from " << inBamFile << endl;
		exit(1);
	}

	hts_idx_t *idx = sam_index_load(in, inBamFile.c_str()); // load index
	if (idx == 0) { // index is unavailable
		cerr << __func__ << ": random alignment retrieval only works for indexed BAM files.\n" << endl;
		exit(1);
	}

	hts_itr_t *iter = sam_itr_querys(idx, header, reg_str.c_str()); // parse a region in the format like `chr2:100-200'
	if (iter == NULL) { // region invalid or reference name not found
		int beg1, end1;
		if (hts_parse_reg(reg_str.c_str(), &beg1, &end1))
			cerr <<  __func__ << ": region " << reg_str << " specifies an unknown reference name." << endl;
		else
			cerr <<  __func__ << ": region " << reg_str << " could not be parsed." << endl;
		exit(1);
	}
	computeAlnDataNumFromIter(in, header, iter, reg_str, qname_vec, qlen_vec, qual_vec, total_len, total_num); // total_num, total_len, q_len
	//cout << "total_len="<< total_len << " total_num=" << total_num << " q_len.size()=" << qlen_vec.size() << endl;

	hts_itr_destroy(iter);

	if(localCov<=MAX_ULTRA_ULTRA_HIGH_COV){
		iter = sam_itr_querys(idx, header, reg_str.c_str());
		if (iter == NULL) { // region invalid or reference name not found
			int beg2, end2;
			if (hts_parse_reg(reg_str.c_str(), &beg2, &end2))
				cerr <<  __func__ << ": region " << reg_str << " specifies an unknown reference name." << endl;
			else
				cerr <<  __func__ << ": region " << reg_str << " could not be parsed." << endl;
			exit(1);
		}

		// load align data from iteration
		loadAlnDataFromIter(alnDataVector, in, header, iter, reg_str, max_ultra_high_cov, qname_vec, qlen_vec, qual_vec, total_len, total_num);

		hts_itr_destroy(iter);
	}

	hts_idx_destroy(idx); // destroy the BAM index
	bam_hdr_destroy(header);
	if(in) sam_close(in);
}

void alnDataLoader::loadAlnData(vector<bam1_t*> &alnDataVector, double max_ultra_high_cov, faidx_t *fai, double max_absig_density, double primary_seg_size_ratio){
	samFile *in = NULL;
	bam_hdr_t *header;
	vector<string> qname_vec;
	vector<int32_t> qlen_vec, qual_vec;
	vector<double> sig_density_vec;
	int64_t startRpos, endRpos;
	size_t total_len = 0, total_num = 0;

	if ((in = sam_open(inBamFile.c_str(), "r")) == 0) {
		cerr << __func__ << ": failed to open " << inBamFile << " for reading" << endl;
		exit(1);
	}

	if ((header = sam_hdr_read(in)) == 0) {
		cerr << __func__ << ": fail to read the header from " << inBamFile << endl;
		exit(1);
	}

	hts_idx_t *idx = sam_index_load(in, inBamFile.c_str()); // load index
	if (idx == 0) { // index is unavailable
		cerr << __func__ << ": random alignment retrieval only works for indexed BAM files.\n" << endl;
		exit(1);
	}

	// compute refspan
	hts_itr_t *iter = sam_itr_querys(idx, header, reg_str.c_str()); // parse a region in the format like `chr2:100-200'
	if (iter == NULL) { // region invalid or reference name not found
		int beg1, end1;
		if (hts_parse_reg(reg_str.c_str(), &beg1, &end1))
			cerr <<  __func__ << ": region " << reg_str << " specifies an unknown reference name." << endl;
		else
			cerr <<  __func__ << ": region " << reg_str << " could not be parsed." << endl;
		exit(1);
	}
	computeRefSpanFromIter(in, header, iter, reg_str, startRpos, endRpos);
	// cout << startRpos << "-" << endRpos << endl;
	hts_itr_destroy(iter);

	// computeAlnDataNumFromIter
	iter = sam_itr_querys(idx, header, reg_str.c_str());
	if (iter == NULL) { // region invalid or reference name not found
		int beg2, end2;
		if (hts_parse_reg(reg_str.c_str(), &beg2, &end2))
			cerr <<  __func__ << ": region " << reg_str << " specifies an unknown reference name." << endl;
		else
			cerr <<  __func__ << ": region " << reg_str << " could not be parsed." << endl;
		exit(1);
	}
	// computeAlnDataNumFromIter(in, header, iter, reg_str, qname_vec, qlen_vec, qual_vec, total_len, total_num); // total_num, total_len, q_len // deleted on 2025-07-03
	computeAlnDataNumFromIter(in, header, iter, reg_str, startRpos, endRpos, fai, max_absig_density, qname_vec, qlen_vec, qual_vec, sig_density_vec, total_len, total_num, primary_seg_size_ratio);
	//computeAlnDataNumFromIter(in, header, iter, reg_str, qlen_vec, total_len, total_num); // total_num, total_len, q_len
	// cout << "total_len="<< total_len << " total_num=" << total_num << " q_len.size()=" << qlen_vec.size() << endl;
	hts_itr_destroy(iter);

	iter = sam_itr_querys(idx, header, reg_str.c_str());
	if (iter == NULL) { // region invalid or reference name not found
		int beg2, end2;
		if (hts_parse_reg(reg_str.c_str(), &beg2, &end2))
			cerr <<  __func__ << ": region " << reg_str << " specifies an unknown reference name." << endl;
		else
			cerr <<  __func__ << ": region " << reg_str << " could not be parsed." << endl;
		exit(1);
	}

	// load align data from iteration
	// loadAlnDataFromIter(alnDataVector, in, header, iter, reg_str, max_ultra_high_cov, qname_vec, qlen_vec, qual_vec, total_len, total_num);
	// loadAlnDataFromIter(alnDataVector, in, header, iter, reg_str, max_ultra_high_cov, qname_vec, qlen_vec, qual_vec, total_len, total_num, fai, max_absig_density, startRpos, endRpos); // deleted on 2025-07-10
	loadAlnDataFromIter(alnDataVector, in, header, iter, reg_str, max_ultra_high_cov, qname_vec, qlen_vec, qual_vec, total_len, total_num);
	// cout << "alnDataVector.size()=" << alnDataVector.size() << endl;

	hts_itr_destroy(iter);
	hts_idx_destroy(idx); // destroy the BAM index
	bam_hdr_destroy(header);
	if(in) sam_close(in);
}

void alnDataLoader::loadAlnData(vector<bam1_t*> &alnDataVector, double max_ultra_high_cov, vector<string> &target_qname_vec){
	samFile *in = NULL;
	bam_hdr_t *header;
	vector<string> qname_vec;
	vector<int32_t> qlen_vec;
	vector<int32_t> qual_vec;
	size_t total_len=0, total_num=0;

	if ((in = sam_open(inBamFile.c_str(), "r")) == 0) {
		cerr << __func__ << ": failed to open " << inBamFile << " for reading" << endl;
		exit(1);
	}

	if ((header = sam_hdr_read(in)) == 0) {
		cerr << __func__ << ": fail to read the header from " << inBamFile << endl;
		exit(1);
	}

	hts_idx_t *idx = sam_index_load(in, inBamFile.c_str()); // load index
	if (idx == 0) { // index is unavailable
		cerr << __func__ << ": random alignment retrieval only works for indexed BAM files.\n" << endl;
		exit(1);
	}

	hts_itr_t *iter = sam_itr_querys(idx, header, reg_str.c_str()); // parse a region in the format like `chr2:100-200'
	if (iter == NULL) { // region invalid or reference name not found
		int beg1, end1;
		if (hts_parse_reg(reg_str.c_str(), &beg1, &end1))
			cerr <<  __func__ << ": region " << reg_str << " specifies an unknown reference name." << endl;
		else
			cerr <<  __func__ << ": region " << reg_str << " could not be parsed." << endl;
		exit(1);
	}
	computeAlnDataNumFromIter(in, header, iter, reg_str, qname_vec, qlen_vec, qual_vec, total_len, total_num); // total_num, total_len, q_len
	//cout << "total_len="<< total_len << " total_num=" << total_num << " q_len.size()=" << qlen_vec.size() << endl;

	hts_itr_destroy(iter);

	if(localCov<=MAX_ULTRA_ULTRA_HIGH_COV){
		iter = sam_itr_querys(idx, header, reg_str.c_str());
		if (iter == NULL) { // region invalid or reference name not found
			int beg2, end2;
			if (hts_parse_reg(reg_str.c_str(), &beg2, &end2))
				cerr <<  __func__ << ": region " << reg_str << " specifies an unknown reference name." << endl;
			else
				cerr <<  __func__ << ": region " << reg_str << " could not be parsed." << endl;
			exit(1);
		}

		// load align data from iteration
		loadAlnDataFromIter(alnDataVector, in, header, iter, reg_str, max_ultra_high_cov, target_qname_vec, qname_vec, qlen_vec, qual_vec, total_len, total_num);

		hts_itr_destroy(iter);
	}

	hts_idx_destroy(idx); // destroy the BAM index
	bam_hdr_destroy(header);
	if(in) sam_close(in);
}

void alnDataLoader::loadAlnData(vector<bam1_t*> &alnDataVector, vector<string> &target_qname_vec){
	samFile *in = NULL;
	bam_hdr_t *header;

	if ((in = sam_open(inBamFile.c_str(), "r")) == 0) {
		cerr << __func__ << ": failed to open " << inBamFile << " for reading" << endl;
		exit(1);
	}

	if ((header = sam_hdr_read(in)) == 0) {
		cerr << __func__ << ": fail to read the header from " << inBamFile << endl;
		exit(1);
	}

	hts_idx_t *idx = sam_index_load(in, inBamFile.c_str()); // load index
	if (idx == 0) { // index is unavailable
		cerr << __func__ << ": random alignment retrieval only works for indexed BAM files.\n" << endl;
		exit(1);
	}

	hts_itr_t *iter = sam_itr_querys(idx, header, reg_str.c_str()); // parse a region in the format like `chr2:100-200'
	if (iter == NULL) { // region invalid or reference name not found
		int beg1, end1;
		if (hts_parse_reg(reg_str.c_str(), &beg1, &end1))
			cerr <<  __func__ << ": region " << reg_str << " specifies an unknown reference name." << endl;
		else
			cerr <<  __func__ << ": region " << reg_str << " could not be parsed." << endl;
		exit(1);
	}

	// load align data from iteration
	loadAlnDataFromIter(alnDataVector, in, header, iter, reg_str, target_qname_vec);

	hts_itr_destroy(iter);
	hts_idx_destroy(idx); // destroy the BAM index
	bam_hdr_destroy(header);
	if(in) sam_close(in);
}

// compute align data number from iteration
double alnDataLoader::computeAlnDataNumFromIter(samFile *in, bam_hdr_t *header, hts_itr_t *iter, string& reg, vector<string> &qname_vec, vector<int32_t> &qlen_vec, vector<int32_t> &qual_vec, size_t &total_len, size_t &total_num){
	int result;
	bam1_t *b;
	string qname;
	double total;
	int64_t num, min_pos, max_pos, begin_pos, end_pos, ref_span;

	total = 0;
	num = 0;
	min_pos = LONG_MAX;
	max_pos = LONG_MIN;
	b = bam_init1();
	while ((result = sam_itr_next(in, iter, b)) >= 0) {
		if(b->core.l_qseq>0 and b->core.qual>=minMapQ and b->core.qual!=255 and !(b->core.flag & BAM_FSECONDARY)){
			total_len += b->core.l_qseq;
			total_num++;
			qname = bam_get_qname(b);
			qname_vec.push_back(qname);
			qlen_vec.push_back(b->core.l_qseq);
			//qual_vec.push_back(b->core.qual);
		}
		if(b->core.l_qseq>0){
			begin_pos = b->core.pos;
			end_pos = bam_endpos(b);
			if(min_pos>begin_pos) min_pos = begin_pos;
			if(max_pos<end_pos) max_pos = end_pos;
			num++;
			total += b->core.l_qseq;
		}
		bam_destroy1(b);
		b = bam_init1();
	}
	bam_destroy1(b);

	mean_read_len = (double)total_len / total_num;
	ref_span = max_pos - min_pos + 1;
	if(num>0) localCov = 4 * (double)total / ref_span;

//	cout << reg << ", total_len=" << total_len << " ; total_num=" << total_num << " ; mean_read_len=" << mean_read_len << ", localCov=" << localCov << ", ref_span=" << ref_span << endl;

	if (result < -1) {
		cerr <<  __func__ << ": retrieval of region " << reg << " failed due to truncated file or corrupt BAM index file." << endl;
		exit(1);
	}

	return localCov;
}

void alnDataLoader::computeRefSpanFromIter(samFile *in, bam_hdr_t *header, hts_itr_t *iter, string& reg, int64_t &startRpos, int64_t &endRpos){
	int result;
	bam1_t *b;
	string qname;
	int64_t current_end_pos;

	b = bam_init1();
	startRpos = LONG_MAX;
	endRpos = LONG_MIN;

	while ((result = sam_itr_next(in, iter, b)) >= 0) {
		if(b->core.l_qseq>0 and b->core.qual>=minMapQ and b->core.qual!=255 and !(b->core.flag & BAM_FSECONDARY)){
			// compute position
			current_end_pos = bam_endpos(b);
			if(b->core.pos < startRpos) startRpos = b->core.pos;
			if(current_end_pos > endRpos) endRpos = current_end_pos;
		}
		bam_destroy1(b);
		b = bam_init1();
	}
	bam_destroy1(b);

	if (result < -1) {
		cerr <<  __func__ << ": retrieval of region " << reg << " failed due to truncated file or corrupt BAM index file." << endl;
		exit(1);
	}
}

void alnDataLoader::computeAlnDataNumFromIter(samFile *in, bam_hdr_t *header, hts_itr_t *iter, string& reg, int64_t startRpos, int64_t endRpos, faidx_t *fai, double max_absig_density, vector<string> &qname_vec, vector<int32_t> &qlen_vec, vector<int32_t> &qual_vec, vector<double> &sig_density_vec, size_t &total_len, size_t &total_num, double primary_seg_size_ratio){
	int result;
	bam1_t *b;
	uint8_t *cigar_int;
	string qname;
	double sig_density;
	map<string, double> qname_max_sig_density_map;
	set<string> mono_seg_set, qname_passed_set;
	vector<string> qname_vec_tmp;
	vector<int32_t> qlen_vec_tmp;
	bool mono_flag;

	qname_vec.clear();
	qlen_vec.clear();
	qual_vec.clear();
	sig_density_vec.clear();
	total_len = 0;
	total_num = 0;

	b = bam_init1();
	while((result = sam_itr_next(in, iter, b)) >= 0){
		if(b->core.l_qseq>0 and b->core.qual>=minMapQ and b->core.qual!=255 and !(b->core.flag & BAM_FSECONDARY)){
			qname = bam_get_qname(b);

			sig_density = computeSigDensityFromAlnData(b, fai, header, startRpos, endRpos);

			cigar_int = bam_aux_get(b, "SA"); // SA
			if(cigar_int==NULL) mono_seg_set.insert(qname);

			//if(sig_density > max_absig_density and bam_endpos(b)-b->core.pos < primary_seg_size_ratio*b->core.l_qseq and cigar_int) sig_density = -1; // exclusion
			//if(sig_density > max_absig_density and bam_endpos(b)-b->core.pos < primary_seg_size_ratio*b->core.l_qseq) sig_density = -1; // exclusion
			if(bam_endpos(b)-b->core.pos < primary_seg_size_ratio*b->core.l_qseq and ((cigar_int and sig_density > max_absig_density) or (cigar_int==NULL and sig_density > ABSIG_DENSITY_SINGLE_ALN_SEG_FACTOR*max_absig_density))) sig_density = -1; // exclusion
			if(qname_max_sig_density_map.find(qname) == qname_max_sig_density_map.end())
				qname_max_sig_density_map[qname] = sig_density;
			else if(sig_density > qname_max_sig_density_map[qname])
				qname_max_sig_density_map[qname] = sig_density;

			qname_vec_tmp.push_back(qname);
			qlen_vec_tmp.push_back(b->core.l_qseq);
			//cout << "qname=" << qname << ", sig_density=" << sig_density << ", max_absig_density=" << max_absig_density << endl;
		}
		bam_destroy1(b);
		b = bam_init1();
	}

	for(auto const& pair : qname_max_sig_density_map){
		const string& qname_key = pair.first;
		const double& max_density = pair.second;
		mono_flag = true;
		if(mono_seg_set.find(qname_key)==mono_seg_set.end()) mono_flag = false;
		if((max_density <= max_absig_density and max_density >= 0) or (mono_flag and max_density <= ABSIG_DENSITY_SINGLE_ALN_SEG_FACTOR*max_absig_density and max_density >= 0)){
//		if(max_density <= max_absig_density and max_density >= 0){
			qname_passed_set.insert(qname_key);
		}
	}

	for(size_t i = 0; i < qname_vec_tmp.size(); ++i){
		if(qname_passed_set.count(qname_vec_tmp.at(i))){
			qname_vec.push_back(qname_vec_tmp.at(i));
			qlen_vec.push_back(qlen_vec_tmp.at(i));

			total_len += qlen_vec_tmp.at(i);
			total_num ++;
		}
	}
	if(total_num > 0) mean_read_len = (double) total_len / total_num;
	else mean_read_len = 0;

//	cout << "total_len=" << total_len << " ;total_num=" << total_num << " ;mean_read_len=" << mean_read_len << endl;
	bam_destroy1(b);
	if (result < -1) {
		cerr <<  __func__ << ": retrieval of region " << reg << " failed due to truncated file or corrupt BAM index file." << endl;
		exit(1);
	}
}

void alnDataLoader::loadAlnDataFromIter(vector<bam1_t*> &alnDataVector, samFile *in, bam_hdr_t *header, hts_itr_t *iter, string& reg, double max_ultra_high_cov, vector<string> &qname_vec, vector<int32_t> &qlen_vec, vector<int32_t> &qual_vec, size_t total_len, size_t total_num){
	int result;
	bam1_t *b;
	double expected_total_bases;//, sampled_cov;
	double compensation_coefficient, local_cov_original;
	size_t i, index, reg_size, num, max_reads_num, total_bases;
	int8_t *selected_flag_array;
	set<string> selected_qname_vec;

	alnDataVector.clear();

	if(qname_vec.size()!=qlen_vec.size()){
		cerr << __func__ << ", line=" << __LINE__ << ": qname_vec.size=" << qname_vec.size() << ", qlen_vec.size=" << qlen_vec.size() << ", error!" << endl;
		exit(1);
	}

	selected_flag_array = (int8_t*) calloc(qlen_vec.size(), sizeof(int8_t));
	if(selected_flag_array==NULL){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot allocate memory, error!" << endl;
		exit(1);
	}

	compensation_coefficient = computeCompensationCoefficient(startRefPos, endRefPos);
	local_cov_original = computeLocalCov(total_len, compensation_coefficient);
	if(max_ultra_high_cov>0 and local_cov_original>max_ultra_high_cov){  // down sample
		reg_size = endRefPos - startRefPos + 1 + mean_read_len;
		expected_total_bases = reg_size * max_ultra_high_cov;
		max_reads_num = qlen_vec.size();

		// ensure each down-sampling is equivalent
		num = total_bases = 0;
		srand(1);
		while(total_bases <= expected_total_bases and num < max_reads_num){
			index = rand() % max_reads_num;
			if(selected_flag_array[index]==0){
				selected_flag_array[index] = 1;
				if(selected_qname_vec.find(qname_vec.at(index))==selected_qname_vec.end()){ // new item
					selected_qname_vec.insert(qname_vec.at(index));
					total_bases += qlen_vec.at(index);
				}
				num ++;
				// cout << "index=" << index << " selected_flag_array[index]=" << selected_flag_array[index] << " total_bases=" << total_bases << ", num=" << num << endl;
			}
		}

		// append remaining align items of the selected reads
		for(i=0; i<qname_vec.size(); i++){
			if(selected_flag_array[i]==0){
				if(selected_qname_vec.find(qname_vec.at(i))!=selected_qname_vec.end()) // found, then add item
					selected_flag_array[i] = 1;
			}
		}
		//sampled_cov = total_bases / reg_size;
		//cout << "sampled_cov=" << sampled_cov << endl;

		b = bam_init1();
		i = 0;
		while((result = sam_itr_next(in, iter, b)) >= 0) {
			if(b->core.l_qseq>0 and b->core.qual>=minMapQ and b->core.qual!=255 and !(b->core.flag & BAM_FSECONDARY)){
				if(selected_flag_array[i]==1) alnDataVector.push_back(b);
				else bam_destroy1(b);
				i++;
			}else bam_destroy1(b);
			b = bam_init1();
		}
	}else{ // load all data
		for(const string& qname_val : qname_vec) { selected_qname_vec.insert(qname_val); /*cout << "qname=" << qname_val << endl;*/ }
		b = bam_init1();
		while((result = sam_itr_next(in, iter, b)) >= 0) {
			if(b->core.l_qseq>0 and b->core.qual>=minMapQ and b->core.qual!=255 and !(b->core.flag & BAM_FSECONDARY)){
				if(selected_qname_vec.count(bam_get_qname(b))) alnDataVector.push_back(b);
				else bam_destroy1(b);
				b = bam_init1();
			}else bam_destroy1(b);
			b = bam_init1();
		}
	}

	bam_destroy1(b);
	if (result < -1) {
		cerr <<  __func__ << ": retrieval of region " << reg << " failed due to truncated file or corrupt BAM index file." << endl;
		exit(1);
	}
	alnDataVector.shrink_to_fit();

	free(selected_flag_array);
}

void alnDataLoader::loadAlnDataFromIter(vector<bam1_t*> &alnDataVector, samFile *in, bam_hdr_t *header, hts_itr_t *iter, string& reg, double max_ultra_high_cov, vector<string> &target_qname_vec, vector<string> &qname_vec, vector<int32_t> &qlen_vec, vector<int32_t> &qual_vec, size_t total_len, size_t total_num){
	int result;
	bam1_t *b;
	double expected_total_bases;
	double compensation_coefficient, local_cov_original;
	size_t i, index, reg_size, num, max_reads_num, total_bases;
	int8_t *selected_flag_array;
	set<string> selected_qname_vec;
	string qname;

	if(qname_vec.size()!=qlen_vec.size()){
		cerr << __func__ << ", line=" << __LINE__ << ": qname_vec.size=" << qname_vec.size() << ", qlen_vec.size=" << qlen_vec.size() << ", error!" << endl;
		exit(1);
	}

	selected_flag_array = (int8_t*) calloc(qlen_vec.size(), sizeof(int8_t));
	if(selected_flag_array==NULL){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot allocate memory, error!" << endl;
		exit(1);
	}

	compensation_coefficient = computeCompensationCoefficient(startRefPos, endRefPos);
	local_cov_original = computeLocalCov(total_len, compensation_coefficient);
	if(max_ultra_high_cov>0 and local_cov_original>max_ultra_high_cov){  // down sample
		reg_size = endRefPos - startRefPos + 1 + mean_read_len;
		expected_total_bases = reg_size * max_ultra_high_cov;
		max_reads_num = qlen_vec.size();

		// mark the target reads
		num = total_bases = 0;
		for(i=0; i<max_reads_num; i++){
			qname = qname_vec.at(i);
			if(find(target_qname_vec.begin(), target_qname_vec.end(), qname) != target_qname_vec.end()){  // found
				selected_flag_array[i] = 1;
				if(selected_qname_vec.find(qname)==selected_qname_vec.end()){ // new item
					selected_qname_vec.insert(qname);
					total_bases += qlen_vec.at(i);
				}
				num ++;
			}
		}

		// ensure each down-sampling is equivalent
		srand(1);
		while(total_bases <= expected_total_bases and num < max_reads_num){
			index = rand() % max_reads_num;
			if(selected_flag_array[index]==0){
				selected_flag_array[index] = 1;
				if(selected_qname_vec.find(qname_vec.at(index))==selected_qname_vec.end()){ // new item
					selected_qname_vec.insert(qname_vec.at(index));
					total_bases += qlen_vec.at(index);
				}
				num ++;
//				cout << "index=" << index << " selected_flag_array[index]=" << selected_flag_array[index] << " total_bases=" << total_bases << ", num=" << num << endl;
			}
		}

		// append remaining align items of the selected reads
		for(i=0; i<qname_vec.size(); i++){
			if(selected_flag_array[i]==0){
				if(selected_qname_vec.find(qname_vec.at(i))!=selected_qname_vec.end()) // found, then add item
					selected_flag_array[i] = 1;
			}
		}

		b = bam_init1();
		i = 0;
		while((result = sam_itr_next(in, iter, b)) >= 0) {
			if(b->core.l_qseq>0 and b->core.qual>=minMapQ and b->core.qual!=255 and !(b->core.flag & BAM_FSECONDARY)){
				if(selected_flag_array[i]==1) alnDataVector.push_back(b);
				else bam_destroy1(b);
				i++;
			}else bam_destroy1(b);
			b = bam_init1();
		}
	}else{ // load all data according to target_qname_vec
		b = bam_init1();
		while ((result = sam_itr_next(in, iter, b)) >= 0) {
			if(b->core.l_qseq>0 and b->core.qual>=minMapQ and b->core.qual!=255 and !(b->core.flag & BAM_FSECONDARY)) alnDataVector.push_back(b); // found
			else bam_destroy1(b);
			b = bam_init1();
		}
	}

	bam_destroy1(b);
	if (result < -1) {
		cerr <<  __func__ << ": retrieval of region " << reg << " failed due to truncated file or corrupt BAM index file." << endl;
		exit(1);
	}
	alnDataVector.shrink_to_fit();

	free(selected_flag_array);
}


void alnDataLoader::loadAlnDataFromIter(vector<bam1_t*> &alnDataVector, samFile *in, bam_hdr_t *header, hts_itr_t *iter, string& reg, vector<string> &qname_vec){
	int result;
	size_t sum, count;
	bam1_t *b;
	string qname;
	bool flag;

	// fetch alignments
	sum = count = 0;
	b = bam_init1();
	while ((result = sam_itr_next(in, iter, b)) >= 0) {
		flag = false;
		if(b->core.l_qseq>0 and b->core.qual>=minMapQ and b->core.qual!=255 and !(b->core.flag & BAM_FSECONDARY)){
			qname = bam_get_qname(b);
			if(find(qname_vec.begin(), qname_vec.end(), qname) != qname_vec.end()){  // found
				sum += b->core.l_qseq;
				count++;
				alnDataVector.push_back(b);
				flag = true;
			}
		}
		if(flag == false) bam_destroy1(b);
		b = bam_init1();
	}
	mean_read_len = (double) sum / count;

	alnDataVector.shrink_to_fit();
	bam_destroy1(b);
	if (result < -1) {
		cerr <<  __func__ << ": retrieval of region " << reg << " failed due to truncated file or corrupt BAM index file." << endl;
		exit(1);
	}
}


double alnDataLoader::computeLocalCov(size_t total_len, double compensation_coefficient){
	double cov = 0;
//	bam1_t *b;
	size_t refined_reg_size;

	refined_reg_size = endRefPos - startRefPos + 1 + 2 * mean_read_len;
	if(refined_reg_size>0){
		cov = total_len / refined_reg_size * compensation_coefficient;
		//cout << "total bases: " << total << " bp, local coverage: " << cov << endl;
	}else{
		cov = -1;
		//cerr << "ERR: ref_seq_size=" << ref_seq_size << endl;
	}
	return cov;
}

double alnDataLoader::computeCompensationCoefficient(size_t startRefPos, size_t endRefPos){
	size_t total_reg_size, refined_reg_size, reg_size_cns;
	double comp_coefficient;

	reg_size_cns = endRefPos - startRefPos + 1;
	total_reg_size = reg_size_cns + 2 * mean_read_len;
	refined_reg_size = reg_size_cns + mean_read_len;  // flanking_area = (2 * mean_read_len) / 2
	comp_coefficient = (double)total_reg_size / refined_reg_size;

	return comp_coefficient;
}

// release the memory
void alnDataLoader::freeAlnData(vector<bam1_t*> &alnDataVector){
	if(!alnDataVector.empty()){
		vector<bam1_t*>::iterator aln;
		for(aln=alnDataVector.begin(); aln!=alnDataVector.end(); aln++)
			bam_destroy1(*aln);
		vector<bam1_t*>().swap(alnDataVector);
		mean_read_len = 0;
	}
}
