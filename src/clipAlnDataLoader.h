#ifndef SRC_CLIPALNDATALOADER_H_
#define SRC_CLIPALNDATALOADER_H_

#include "structures.h"

using namespace std;


class clipAlnDataLoader {
	public:
		string chrname, inBamFile;
		int64_t startRefPos, endRefPos;
		int32_t minClipEndSize;
		int32_t minMapQ, minHighMapQ;
	public:
		clipAlnDataLoader(string &chrname, int64_t startRefPos, int64_t endRefPos, string &inBamFile, int32_t minClipEndSize, int32_t minMapQ, int32_t minHighMapQ);
		virtual ~clipAlnDataLoader();
		void loadClipAlnData(vector<clipAlnData_t*> &clipAlnDataVector);
		void loadClipAlnData(vector<clipAlnData_t*> &clipAlnDataVector, double max_ultra_high_cov);
		void loadClipAlnData(vector<clipAlnData_t*> &clipAlnDataVector, double max_ultra_high_cov, vector<string> &qname_vec);
		void loadClipAlnData(vector<clipAlnData_t*> &clipAlnDataVector, vector<string> &qname_vec);
		void loadClipAlnDataWithSATag(vector<clipAlnData_t*> &clipAlnDataVector);
		void loadClipAlnDataWithSATag(vector<clipAlnData_t*> &clipAlnDataVector, double max_ultra_high_cov);
		void loadClipAlnDataWithSATagWithSegSize(vector<clipAlnData_t*> &clipAlnDataVector, double max_ultra_high_cov, double primary_seg_size_ratio, double primary_seg_nm_ratio);
		void loadClipAlnDataWithSATag(vector<clipAlnData_t*> &clipAlnDataVector, double max_ultra_high_cov, vector<string> &qname_vec);
		void loadClipAlnDataWithSATag(vector<clipAlnData_t*> &clipAlnDataVector, vector<string> &qname_vec);

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
		void addAdjacentInfo(vector<clipAlnData_t*> &clipAlnDataVector);
		void orderClipAlnSegsSingleQuery(vector<clipAlnData_t*> &query_aln_vec);
		void assignSideMostFlag(vector<clipAlnData_t*> &query_aln_vec);
		void removeClipAlnDataWithLowPrimarySegSizeRatio(vector<clipAlnData_t*> &clipAlnDataVector, double primary_seg_size_ratio, double primary_seg_nm_ratio);
		void removeSingleItemClipAlnData(vector<clipAlnData_t*> &clipAlnDataVector, int32_t idx);
};

#endif /* SRC_CLIPALNDATALOADER_H_ */
