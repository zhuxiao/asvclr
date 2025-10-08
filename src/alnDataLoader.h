
#ifndef SRC_ALNDATALOADER_H_
#define SRC_ALNDATALOADER_H_

#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <algorithm>
#include <map>

#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/faidx.h>

#include "util.h"

using namespace std;

#define ABSIG_DENSITY_SINGLE_ALN_SEG_FACTOR		(1.5f)

class alnDataLoader {
	public:
		string reg_str, inBamFile;
		double mean_read_len, localCov;
		int32_t startRefPos, endRefPos, minMapQ, minHighMapQ;

	public:
		//alnDataLoader();
		alnDataLoader(string &chrname, int32_t startRefPos, int32_t endRefPos, string &inBamFile, int32_t minMapQ, int32_t minHighMapQ);
		virtual ~alnDataLoader();
//		void loadAlnData(vector<bam1_t*> &alnDataVector);
		void loadAlnData(vector<bam1_t*> &alnDataVector, double max_ultra_high_cov);
		void loadAlnData(vector<bam1_t*> &alnDataVector, double max_ultra_high_cov, faidx_t *fai, double max_absig_density, double primary_seg_size_ratio);
		void loadAlnData(vector<bam1_t*> &alnDataVector, double max_ultra_high_cov, vector<string> &qname_vec);
		void loadAlnData(vector<bam1_t*> &alnDataVector, vector<string> &qname_vec);
		void freeAlnData(vector<bam1_t*> &alnDataVector);

	private:
		//void computeAlnDataNumFromIter(samFile *in, bam_hdr_t *header, hts_itr_t *iter, string& reg, vector<int32_t> &qlen_vec, size_t &total_len, size_t &total_num);
		double computeAlnDataNumFromIter(samFile *in, bam_hdr_t *header, hts_itr_t *iter, string& reg, vector<string> &qname_vec, vector<int32_t> &qlen_vec, vector<int32_t> &qual_vec, size_t &total_len, size_t &total_num);
		void computeAlnDataNumFromIter(samFile *in, bam_hdr_t *header, hts_itr_t *iter, string& reg, int64_t startRpos, int64_t endRpos, faidx_t *fai, double max_absig_density, vector<string> &qname_vec, vector<int32_t> &qlen_vec, vector<int32_t> &qual_vec, vector<double> &sig_density_vec, size_t &total_len, size_t &total_num, double primary_seg_size_ratio);
//		void loadAlnDataFromIter(vector<bam1_t*> &alnDataVector, samFile *in, bam_hdr_t *header, hts_itr_t *iter, string& reg);
		// void loadAlnDataFromIter(vector<bam1_t*> &alnDataVector, samFile *in, bam_hdr_t *header, hts_itr_t *iter, string& reg, double max_ultra_high_cov, vector<string> &qname_vec, vector<int32_t> &qlen_vec, vector<int32_t> &qual_vec, size_t total_len, size_t total_num);
		void loadAlnDataFromIter(vector<bam1_t*> &alnDataVector, samFile *in, bam_hdr_t *header, hts_itr_t *iter, string& reg, double max_ultra_high_cov, vector<string> &qname_vec, vector<int32_t> &qlen_vec, vector<int32_t> &qual_vec, size_t total_len, size_t total_num);
		void loadAlnDataFromIter(vector<bam1_t*> &alnDataVector, samFile *in, bam_hdr_t *header, hts_itr_t *iter, string& reg, double max_ultra_high_cov, vector<string> &target_qname_vec, vector<string> &qname_vec, vector<int32_t> &qlen_vec, vector<int32_t> &qual_vec, size_t total_len, size_t total_num);
		void loadAlnDataFromIter(vector<bam1_t*> &alnDataVector, samFile *in, bam_hdr_t *header, hts_itr_t *iter, string& reg, vector<string> &qname_vec);
		void computeRefSpanFromIter(samFile *in, bam_hdr_t *header, hts_itr_t *iter, string& reg, int64_t &startRpos, int64_t &endRpos);
		double computeLocalCov(size_t total_len, double compensation_coefficient);
		double computeCompensationCoefficient(size_t startRefPos, size_t endRefPos);
};

#endif /* SRC_ALNDATALOADER_H_ */
