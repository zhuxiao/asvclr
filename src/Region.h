#ifndef SRC_REGION_H_
#define SRC_REGION_H_

#include <iostream>
#include <fstream>
#include <htslib/faidx.h>

#include "structures.h"
#include "Paras.h"
#include "Base.h"
#include "misAlnReg.h"

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

#define MIN_REG_SIZE_EXTRACT_SIG			30 // added on 2025-03-06

#define LARGE_INDEL_RATIO_THRES				(0.1f)
#define HIGH_INDEL_CLIP_RATIO_THRES			(0.6f)
#define SECOND_INDEL_CLIP_RATIO_THRES		(0.3f)
#define HIGH_INDEL_CLIP_BASE_RATIO_THRES	(0.1f)
#define MIN_VALID_SIG_COV_RATIO				(0.05f)


class Region {
	public:
		Paras *paras;
		faidx_t *fai;
		string chrname;
		int64_t startRPos, endRPos, startMidPartPos, endMidPartPos, chrlen, minRPos, maxRPos;
		Base *regBaseArr;  // the local region base array
		size_t regFlag, subRegSize;
		bool wholeRefGapFlag; // true -- if all the bases in the region of reference are 'N'; false -- otherwise

		double meanBlockCov, localRegCov, refinedLocalRegCov;
		size_t highIndelSubRegNum;

		// candidate flag
		bool indelCandFlag;
		size_t difType;

		// output file
		string out_dir_cns;

	private:
		vector<int64_t> disagrePosVector;  // element: the reference position with disagreement
		vector<int64_t> zeroCovPosVector;  // element: the reference position with zero coverage
		vector<reg_t*> abCovRegVector;    // element: the reference region with abnormal coverage, format: pos1-pos2

		// SNV and indel
		vector<int64_t> snvVector;
		vector<reg_t*> indelVector;
		vector<reg_t*> clipRegVector;  // only for duplication and inversion

	public:
		Region(string& chrname, int64_t startRpos, int64_t endRPos, int64_t chrlen, int64_t minRPos, int64_t maxRPos, Base *regBaseArr, size_t regFlag, Paras *paras, faidx_t *fai);
		virtual ~Region();
		void setOutputDir(string& out_dir_cns_prefix);

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
		vector<int64_t> getSnvVector();
		vector<reg_t*> getIndelVector();
		vector<reg_t*> getClipRegVector();
		vector<int64_t> getZeroCovPosVector();

	private:
		bool IsWholeRefGap();
		void addDisagrePos(int64_t pos);
		void addZeroCovPos(int64_t pos);
		void destroyDisagrePosVector();
		void destroyZeroCovPosVector();
		void destroyAbCovRegVector();
		void destroySnvVector();
		void destroyIndelVector();
		void destroyClipRegVector();
		int computeAbCovReg();
		reg_t* allocateReg(string &chrname, int64_t startPosReg, int64_t endPosReg, int32_t sv_len);
		double computeMeanCovReg(int64_t startPosReg, int64_t endPosReg);
		double computeRefinedMeanCovReg(int64_t startPosReg, int64_t endPosReg);
		int computeHighIndelEventRegNum();
		int32_t computeReadIndelEventNumReg(int64_t startPosReg, int64_t endPosReg);
		bool computeSNVFlag(int64_t pos, int64_t startCheckPos, int64_t endCheckPos);
		bool haveMuchShortIndelsAround(int64_t startCheckPos, int64_t endCheckPos);
		int32_t computeValidSigNumReg(int64_t startPosReg, int64_t endPosReg, int32_t min_sig_size);
		reg_t* getIndelReg(int64_t startCheckPos);
		bool haveNoAbSigs(Base *base, int64_t pos);
		int32_t getMismatchBasesAround(int64_t pos1, int64_t pos2);
		int32_t getDisZeroCovNum(int64_t startPos, int64_t endPos);
		int32_t getLargeIndelBaseNum(int64_t startPos, int64_t endPos);
		int32_t getLargeIndelNum(int64_t startPos, int64_t endPos);
		int32_t rescueByLargeIndelNum(int64_t startPos, int64_t endPos, int32_t large_indel_size_thres);
		int32_t getHighConIndelNum(int64_t startPos, int64_t endPos, float threshold, float polymer_ignore_ratio_thres);
		int32_t computeEstSVLen(int64_t startPos, int64_t endPos);

		// duplication and inversion
		reg_t* getClipReg(int64_t startCheckPos);
		bool haveNoClipSig(int64_t startPos, int64_t endPos, double clip_ratio_thres);
};

#endif /* SRC_REGION_H_ */
