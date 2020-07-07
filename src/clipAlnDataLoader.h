#ifndef SRC_CLIPALNDATALOADER_H_
#define SRC_CLIPALNDATALOADER_H_

#include "structures.h"

using namespace std;


class clipAlnDataLoader {
	public:
		string chrname, inBamFile;
		size_t startRefPos, endRefPos;

	public:
		clipAlnDataLoader(string &chrname, int32_t startRefPos, int32_t endRefPos, string &inBamFile);
		virtual ~clipAlnDataLoader();
		void loadClipAlnData(vector<clipAlnData_t*> &clipAlnDataVector);
		void loadClipAlnData(vector<clipAlnData_t*> &clipAlnDataVector, double max_ultra_high_cov);
		void loadClipAlnDataWithSATag(vector<clipAlnData_t*> &clipAlnDataVector);
		void loadClipAlnDataWithSATag(vector<clipAlnData_t*> &clipAlnDataVector, double max_ultra_high_cov);
		void freeClipAlnData(vector<clipAlnData_t*> &clipAlnDataVector);

	private:
		void samplingAlnData(vector<bam1_t*> &alnDataVector, double mean_read_len, double max_ultra_high_cov);
		double samplingAlnDataOp(vector<bam1_t*> &alnDataVector, double mean_read_len, double expect_cov_val);
		double computeLocalCov(vector<bam1_t*> &alnDataVector, double mean_read_len, double compensation_coefficient);
		double computeCompensationCoefficient(size_t startRefPos, size_t endRefPos, double mean_read_len);
		clipAlnData_t* generateClipAlnData(bam1_t* bam, bam_hdr_t *header);
		void fillClipAlnDataBySATag(vector<clipAlnData_t*> &clipAlnDataVector);
		clipAlnData_t* addNewSAItemToClipAlnDataVec(string &queryname, string &aln_seg_info_str, vector<clipAlnData_t*> &clipAlnDataVector);
		void parseSingleAlnStrSA(clipAlnData_t &clip_aln_ret, string &aln_seg_info_str);
		bool isSameClipAlnSeg(clipAlnData_t *clip_aln1, clipAlnData_t *clip_aln2);
};

#endif /* SRC_CLIPALNDATALOADER_H_ */
