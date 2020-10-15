#ifndef SRC_COVLOADER_H_
#define SRC_COVLOADER_H_

#include <htslib/sam.h>

#include "Base.h"
#include "RefSeqLoader.h"

#include "structures.h"

using namespace std;

#define MD_MISMATCH   				10

class covLoader {
	public:
		string chrname;
		size_t startPos, endPos;
		size_t min_ins_size_filt, min_del_size_filt;
		faidx_t *fai;

	public:
		covLoader(string &chrname, size_t startPos, size_t endPos, faidx_t *fai);
		covLoader(string &chrname, size_t startPos, size_t endPos, faidx_t *fai, size_t min_ins_size_filt, size_t min_del_size_filt);
		virtual ~covLoader();
		Base *initBaseArray();
		void freeBaseArray(Base *baseArray);
		void generateBaseCoverage(Base *baseArr, vector<bam1_t*> alnDataVector);

	private:
		int assignRefBase(Base *baseArray, faidx_t *fai);
		vector<struct alnSeg*> generateAlnSegs(bam1_t* b);
		struct alnSeg* allocateAlnSeg(size_t startRpos, size_t startQpos, size_t seglen, size_t opFlag, string seg_MD);
		void destroyAlnSegs(vector<struct alnSeg*> alnSegs);
		vector<struct MD_seg*> extractMDSegs(bam1_t* b);
		struct MD_seg* allocateMDSeg(string& seg, size_t opflag);
		void destroyMDSeg(vector<struct MD_seg*> segs_MD);
		int updateBaseInfo(Base *baseArr, vector<struct alnSeg*> alnSegs);
		void updateBaseCovInfo(Base *baseArr);
		void computeDelNumFromDelVec(Base *baseArr);

		// testing functions
		void outputAlnSegs(vector<struct alnSeg*> alnSegs);
		void outputMDSeg(vector<struct MD_seg*> segs_MD);
		bool haveOpflagCigar(bam1_t* b, size_t opfalg);
};

#endif /* SRC_COVLOADER_H_ */
