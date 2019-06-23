#ifndef SRC_EVENTS_H_
#define SRC_EVENTS_H_

#include <iostream>
#include <string>

using namespace std;

#define ALN_PLUS_ORIENT				0   // plus orientation
#define ALN_MINUS_ORIENT			1   // minus orientation

#define VAR_UNC				0	// uncertain
#define VAR_INS				1	// insertion
#define VAR_DEL				2	// deletion
#define VAR_DUP				3	// duplication
#define VAR_INV				4	// inversion
#define VAR_TRA				5	// translocation
#define VAR_BND				6	// translocation
#define VAR_INV_TRA			7	// inverted translocation
#define VAR_MIX				10	// mixed variation

#define CTG_END_SKIP_SIZE		1000

// base coverage structure
typedef struct
{
	// [0..4]: A, C, G, T, N, sum(A+C+G+T+N)
	// idx_RefBase points to the element in the num_bases[] array
	uint16_t num_bases[6], idx_RefBase;
	uint16_t idx_max, num_max;  // the maximal base index and the corresponding base count
}baseCoverage_t;

typedef struct{
	uint32_t startPos;
	string seq;
}indelEvent_t;

typedef struct{
	uint32_t startPos;
	uint16_t opflag;
	uint16_t endFlag;  // 0: head; 1: tail
	string seq;
}clipEvent_t;

typedef indelEvent_t insEvent_t;
typedef indelEvent_t delEvent_t;

#define allocateInsEvent(pos, seq)  (allocateIndelEvent(pos, seq))
#define allocateDelEvent(pos, seq)  (allocateIndelEvent(pos, seq))

indelEvent_t* allocateIndelEvent(uint32_t startPos, string& seq);
clipEvent_t* allocateClipEvent(uint32_t startPos, uint16_t opflag, uint16_t endFlag, string& seq);

#endif /* SRC_EVENTS_H_ */
