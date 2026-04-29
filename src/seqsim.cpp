#include "seqsim.h"

// computeVarseqSim	upperSeq
double computeVarseqSim(const string& seqA, const string& seqB) {
	double similarity = 0;
	int edit_dist, max_len;
	string  seqA_tmp, seqB_tmp;

	seqA_tmp = seqA;
	seqB_tmp = seqB;
	upperSeq(seqA_tmp);
	upperSeq(seqB_tmp);

	EdlibAlignResult result = edlibAlign(seqA_tmp.c_str(), seqA_tmp.size(), seqB_tmp.c_str(), seqB_tmp.size(), edlibDefaultAlignConfig());
	if (result.status == EDLIB_STATUS_OK) {
		edit_dist = result.editDistance;
		max_len = max(seqA_tmp.size(), seqB_tmp.size());
		similarity = 1.0 - (double)edit_dist / max_len;

	}
	edlibFreeAlignResult(result);

	return similarity;
}

