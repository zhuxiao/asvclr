
#ifndef SRC_ALNDATALOADER_H_
#define SRC_ALNDATALOADER_H_

#include <iostream>
#include <string>
#include <vector>

#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/faidx.h>

#include "Paras.h"

using namespace std;

class alnDataLoader {
	public:
		string reg_str, inBamFile;
		double mean_read_len;

	public:
		//alnDataLoader();
		alnDataLoader(string &chrname, int32_t startRefPos, int32_t endRefPos, string &inBamFile);
		virtual ~alnDataLoader();
		void loadAlnData(vector<bam1_t*> &alnDataVector);
		void freeAlnData(vector<bam1_t*> &alnDataVector);

	private:
		void loadAlnDataFromIter(vector<bam1_t*> &alnDataVector, samFile *in, bam_hdr_t *header, hts_itr_t *iter, string& reg);
};

#endif /* SRC_ALNDATALOADER_H_ */
