
#ifndef SRC_ALNDATALOADER_H_
#define SRC_ALNDATALOADER_H_

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/faidx.h>

//#include "Paras.h"

using namespace std;

class alnDataLoader {
	public:
		string reg_str, inBamFile;
		double mean_read_len;
		int32_t startRefPos, endRefPos, minMapQ;

	public:
		//alnDataLoader();
		alnDataLoader(string &chrname, int32_t startRefPos, int32_t endRefPos, string &inBamFile, int32_t minMapQ);
		virtual ~alnDataLoader();
//		void loadAlnData(vector<bam1_t*> &alnDataVector);
		void loadAlnData(vector<bam1_t*> &alnDataVector, double max_ultra_high_cov);
		void loadAlnData(vector<bam1_t*> &alnDataVector, vector<string> &qname_vec);
		void freeAlnData(vector<bam1_t*> &alnDataVector);

	private:
		void computeAlnDataNumFromIter(samFile *in, bam_hdr_t *header, hts_itr_t *iter, string& reg, vector<int32_t> &qlen_vec, size_t &total_len, size_t &total_num);
//		void loadAlnDataFromIter(vector<bam1_t*> &alnDataVector, samFile *in, bam_hdr_t *header, hts_itr_t *iter, string& reg);
		void loadAlnDataFromIter(vector<bam1_t*> &alnDataVector, samFile *in, bam_hdr_t *header, hts_itr_t *iter, string& reg, double max_ultra_high_cov, vector<int32_t> &qlen_vec, size_t total_len, size_t total_num);
		void loadAlnDataFromIter(vector<bam1_t*> &alnDataVector, samFile *in, bam_hdr_t *header, hts_itr_t *iter, string& reg, vector<string> &qname_vec);
		double computeLocalCov(size_t total_len, double compensation_coefficient);
		double computeCompensationCoefficient(size_t startRefPos, size_t endRefPos);
};

#endif /* SRC_ALNDATALOADER_H_ */
