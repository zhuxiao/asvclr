#ifndef SRC_REGION_H_
#define SRC_REGION_H_

#include <iostream>
#include <fstream>
#include <htslib/faidx.h>

#include "Paras.h"
#include "Base.h"
#include "misAlnReg.h"
#include "clipAlnDataLoader.h"

using namespace std;

#define HEAD_REGION		0
#define INNER_REGION	1
#define TAIL_REGION		2

#define REG_DIF_NONE		0
#define REG_DIF_SNV			1
#define REG_DIF_INDEL		2
#define REG_DIF_BOTH		3
#define REG_DIF_UNCERTAIN	4

#define SUB_REG_SIZE			10
#define SUB_CLIP_REG_SIZE		100

#define LARGE_INDEL_RATIO_THRES		(0.1f)

typedef struct{
	string chrname;
	size_t startRefPos, endRefPos, startLocalRefPos, endLocalRefPos, startQueryPos, endQueryPos;
	size_t var_type, aln_orient, dup_num;
	int32_t query_id, sv_len, blat_aln_id;
	string refseq, altseq;
	bool call_success_status, short_sv_flag;
}reg_t;


class Region {
	public:
		Paras *paras;
		faidx_t *fai;
		string chrname;
		size_t startRPos, endRPos, startMidPartPos, endMidPartPos, chrlen;
		Base *regBaseArr;  // the local region base array
		size_t regFlag, subRegSize;
		bool wholeRefGapFlag; // true -- if all the bases in the region of reference are 'N'; false -- otherwise

		double meanBlockCov, localRegCov;
		size_t highIndelSubRegNum;

		// candidate flag
		bool indelCandFlag;
		size_t difType;

		// output file
		string out_dir_assemble;

	private:
		vector<size_t> disagrePosVector;  // element: the reference position with disagreement
		vector<size_t> zeroCovPosVector;  // element: the reference position with zero coverage
		vector<reg_t*> abCovRegVector;    // element: the reference region with abnormal coverage, format: pos1-pos2

		// SNV and indel
		vector<size_t> snvVector;
		vector<reg_t*> indelVector;
		vector<reg_t*> clipRegVector;  // only for duplication and inversion

	public:
		Region(string& chrname, size_t startRpos, size_t endRPos, size_t chrlen, Base *regBaseArr, size_t regFlag, Paras *paras);
		virtual ~Region();
		void setOutputDir(string& out_dir_assemble_prefix);

	public:
		void setMeanBlockCov(double meanBlockCov);
		int computeRegAbSigs();
		int computeDisagreements();
		void printRegAbSigs();
		void detectSNV();
		void detectIndelReg();
		void detectHighClipReg();
		size_t getDisagreeNum();
		size_t getSNVNum();
		size_t getIndelNum();
		void determineDifType();
		vector<size_t> getSnvVector();
		vector<reg_t*> getIndelVector();
		vector<reg_t*> getClipRegVector();
		vector<size_t> getZeroCovPosVector();

	private:
		bool IsWholeRefGap();
		void addDisagrePos(size_t pos);
		void addZeroCovPos(size_t pos);
		void destroyDisagrePosVector();
		void destroyZeroCovPosVector();
		void destroyAbCovRegVector();
		void destroySnvVector();
		void destroyIndelVector();
		void destroyClipRegVector();
		int computeAbCovReg();
		reg_t* allocateReg(size_t startPosReg, size_t endPosReg);
		double computeMeanCovReg(size_t startPosReg, size_t endPosReg);
		int computeHighIndelEventRegNum();
		size_t computeReadIndelEventNumReg(size_t startPosReg, size_t endPosReg);
		bool computeSNVFlag(size_t pos, size_t startCheckPos, size_t endCheckPos);
		bool haveMuchShortIndelsAround(size_t startCheckPos, size_t endCheckPos);
		reg_t* getIndelReg(size_t startCheckPos);
		bool haveNoAbSigs(Base *base, size_t pos);
		size_t getMismatchBasesAround(size_t pos1, size_t pos2);
		size_t getDisZeroCovNum(size_t startPos, size_t endPos);
		size_t getLargeIndelBaseNum(size_t startPos, size_t endPos);
		size_t getLargeIndelNum(size_t startPos, size_t endPos);

		// duplication and inversion
		reg_t* getClipReg(size_t startCheckPos);
		bool haveNoClipSig(size_t startPos, size_t endPos, double clip_ratio_thres);
};

#endif /* SRC_REGION_H_ */
