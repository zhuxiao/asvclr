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
	uint32_t num_bases[6];
	int8_t idx_RefBase; // idx_RefBase: 5 for mixed base symbols
	int8_t idx_max;
	int32_t num_max;  // the maximal base index and the corresponding base count
	char refBase; // A, C, G, T, N(A+C+G+T), M(A+C), R(A+G), S(C+G), V(A+C+G), W(A+T), Y(C+T), H(A+C+T), K(G+T), D(A+G+T), B(C+G+T), ACMGRSVTWYHKDBN
	bool polymer_flag;
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
