#ifndef SRC_BASE_H_
#define SRC_BASE_H_

#include <vector>
#include "events.h"

#define DISAGREE_THRES			(0.80)  // important, and need to be estimated
#define DISAGREE_THRES_REG		(0.85)
#define DISAGREE_NUM_THRES_REG	2
#define DISAGREE_CHK_REG		2
#define MIN_RATIO_SNV			(0.75)


// Base
class Base
{
	public:
		baseCoverage_t coverage;
		vector<insEvent_t*> insVector;  // insertion event vector
		vector<delEvent_t*> delVector;  // deletion event vector
		vector<delEvent_t*> extendDelVector;  // extend deletion event vector, 20220530     (10k)
		vector<clipEvent_t*> clipVector;  // clip event vector

		int32_t num_shortIns, num_shortdel, num_shortClip, del_num_from_del_vec;
		int32_t maxConIndelEventNum: 27, max_con_type: 5;		// max_con_type: BAM_CINS, BAM_CDEL, BAM_CSOFT_CLIP, BAM_CHARD_CLIP, BAM_CDIFF, etc.
		float maxConIndelEventRatio;  // the ratio of maximum number of consensus indel events to the sum of total coverage and deletions (i.e. shadow coverage)

	public:
		Base();
		virtual ~Base();
		void destroyBase();
		void addInsEvent(insEvent_t* insE);
		void addDelEvent(delEvent_t* delE);
		void addClipEvent(clipEvent_t* clipE);
		void updateCovInfo();
		bool isDisagreeBase();
		bool isZeroCovBase();
		bool isHighIndelBase(float threshold, float ignore_polymer_ratio_threshold);
		bool isHighConIndelBase(float threshold, float ignore_polymer_ratio_threshold);
		bool isMatchToRef();
		size_t getLargeIndelNum(size_t thres);
		size_t getRegLargeIndelNum(size_t thres);
		size_t getTotalIndelNum();
		size_t getTotalClipNum();
		size_t getTotalCovNum();

	private:
		void init();
		void destroyInsVector();
		void destroyDelVector();
		void destroyClipVector();
};

#endif /* SRC_BASE_H_ */
