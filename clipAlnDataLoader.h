#ifndef SRC_CLIPALNDATALOADER_H_
#define SRC_CLIPALNDATALOADER_H_

#include <iostream>
#include <string>
#include <vector>
#include <htslib/sam.h>

using namespace std;


typedef struct{
	bam1_t *bam;
	string queryname, chrname;
	size_t querylen, aln_orient, startRefPos, endRefPos, startQueryPos, endQueryPos, leftClipSize, rightClipSize;
	bool leftHardClippedFlag, rightHardClippedFlag;
	bool left_clip_checked_flag, right_clip_checked_flag, query_checked_flag, SA_tag_flag;
}clipAlnData_t;


class clipAlnDataLoader {
	public:
		string chrname, inBamFile;
		size_t startRefPos, endRefPos;

	public:
		clipAlnDataLoader(string &chrname, int32_t startRefPos, int32_t endRefPos, string &inBamFile);
		virtual ~clipAlnDataLoader();
		void loadClipAlnData(vector<clipAlnData_t*> &clipAlnDataVector);
		void loadClipAlnDataWithSATag(vector<clipAlnData_t*> &clipAlnDataVector);
		void freeClipAlnData(vector<clipAlnData_t*> &clipAlnDataVector);

	private:
		clipAlnData_t* generateClipAlnData(bam1_t* bam, bam_hdr_t *header);
		void fillClipAlnDataBySATag(vector<clipAlnData_t*> &clipAlnDataVector);
		clipAlnData_t* addNewSAItemToClipAlnDataVec(string &queryname, string &aln_seg_info_str, vector<clipAlnData_t*> &clipAlnDataVector);
		void parseSingleAlnStrSA(clipAlnData_t &clip_aln_ret, string &aln_seg_info_str);
		bool isSameClipAlnSeg(clipAlnData_t *clip_aln1, clipAlnData_t *clip_aln2);
};

#endif /* SRC_CLIPALNDATALOADER_H_ */
