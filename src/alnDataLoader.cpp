
#include "alnDataLoader.h"

//alnDataLoader::alnDataLoader(){
//}

alnDataLoader::alnDataLoader(string &chrname, int32_t startRefPos, int32_t endRefPos, string &inBamFile) {
	this->reg_str = chrname + ":" + to_string(startRefPos) + "-" + to_string(endRefPos);
	this->inBamFile = inBamFile;
	mean_read_len = 0;
}

alnDataLoader::~alnDataLoader() {
}

void alnDataLoader::loadAlnData(vector<bam1_t*> &alnDataVector){
	samFile *in = 0;
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
		int beg, end;
		if (hts_parse_reg(reg_str.c_str(), &beg, &end))
			cerr <<  __func__ << ": region " << reg_str << " specifies an unknown reference name." << endl;
		else
			cerr <<  __func__ << ": region " << reg_str << " could not be parsed." << endl;
		exit(1);
	}

	// load align data from iteration
	loadAlnDataFromIter(alnDataVector, in, header, iter, reg_str);

	hts_itr_destroy(iter);
	hts_idx_destroy(idx); // destroy the BAM index
	bam_hdr_destroy(header);
	if(in) sam_close(in);
}

// load align data from iteration
void alnDataLoader::loadAlnDataFromIter(vector<bam1_t*> &alnDataVector, samFile *in, bam_hdr_t *header, hts_itr_t *iter, string& reg){
	int result;
	size_t sum, count;
	bam1_t *b;

	// fetch alignments
	sum = count = 0;
	b = bam_init1();
	while ((result = sam_itr_next(in, iter, b)) >= 0) {
		sum += b->core.l_qseq;
		count++;
		alnDataVector.push_back(b);
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
