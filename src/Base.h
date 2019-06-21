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
		vector<clipEvent_t*> clipVector;  // deletion event vector

		size_t num_shortIns, num_shortdel, num_shortClip;

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
		bool isHighIndelBase(float threshold);
		size_t getLargeIndelNum(size_t thres);
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
