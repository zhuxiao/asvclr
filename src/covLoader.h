#ifndef SRC_COVLOADER_H_
#define SRC_COVLOADER_H_

#include <htslib/sam.h>

#include "Base.h"
#include "RefSeqLoader.h"

#include "structures.h"

using namespace std;

#define BAM_INVALID					0
#define BAM_CIGAR_NO_DIFF_MD		1
#define BAM_CIGAR_NO_DIFF_NO_MD		2
#define BAM_CIGAR_DIFF_NO_MD		3
#define BAM_CIGAR_DIFF_MD			4

class covLoader {
	public:
		string chrname;
		int64_t startPos, endPos;
		int64_t min_ins_size_filt:20, min_del_size_filt:20, bam_type:24;
		char left_ref_base, right_ref_base;
		faidx_t *fai;

	public:
		covLoader(string &chrname, int64_t startPos, int64_t endPos, faidx_t *fai);
		covLoader(string &chrname, int64_t startPos, int64_t endPos, faidx_t *fai, int32_t min_ins_size_filt, int32_t min_del_size_filt);
		virtual ~covLoader();
		Base *initBaseArray();
		void freeBaseArray(Base *baseArray);
		void generateBaseCoverage(Base *baseArr, vector<bam1_t*> &alnDataVector);

	private:
		int assignRefBase(Base *baseArray, faidx_t *fai);
		void assignPolymerFlag(Base *baseArray);
//		vector<struct alnSeg*> generateAlnSegs(bam1_t* b);
//		vector<struct alnSeg*> generateAlnSegs_no_MD(bam1_t* b, Base* BaseArr);
//		struct alnSeg* allocateAlnSeg(int64_t startRpos, int64_t startQpos, int32_t seglen, uint32_t opFlag, string &seg_MD);
//		void destroyAlnSegs(vector<struct alnSeg*> &alnSegs);
//		int getBamType(bam1_t *b);
//		vector<struct MD_seg*> extractMDSegs(bam1_t* b);
//		struct MD_seg* allocateMDSeg(string& seg, uint32_t opflag);
//		void destroyMDSeg(vector<struct MD_seg*> &segs_MD);
		int updateBaseInfo(Base *baseArr, vector<struct alnSeg*> &alnSegs);
		void updateBaseCovInfo(Base *baseArr);
		void computeDelNumFromDelVec(Base *baseArr);
		void computeConIndelEventRatio(Base *baseArr);
};

#endif /* SRC_COVLOADER_H_ */
