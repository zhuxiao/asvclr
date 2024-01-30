#ifndef SRC_CLIPREGCLUSTER_H_
#define SRC_CLIPREGCLUSTER_H_

#include "structures.h"
#include "RefSeqLoader.h"
#include "util.h"

using namespace std;

#define CLIP_REG_CLUSTER_DEBUG				0

#define MAX_DIST_MATCH_CLIP_POS				500
#define MIN_SIZE_RATIO_MATCH_CLIP_POS		(0.8f)
#define QC_SIZE_RATIO_MATCH_THRES_INDEL		(0.7f)

class clipRegCluster {
private:
	int32_t minClipEndSize, min_sv_size;

public:
	clipRegCluster(int32_t minClipEndSize, int32_t min_sv_size);
	virtual ~clipRegCluster();
	void destoryQcSigList(vector<qcSigList_t*> &qcSigList_vec);
	vector<qcSigList_t*> extractQcSigsClipReg(vector<struct querySeqInfoNode*> &query_seq_info_vec);
	vector<qcSigListVec_t*> queryCluster(vector<struct querySeqInfoNode*> &query_seq_info_vec);

private:
	vector<struct querySeqInfoNode*> getQuerySeqs(string &qname, vector<struct querySeqInfoNode*> &query_seq_info_vec);
	void printQcSigListVec(vector<qcSigList_t*> &qcSigList_vec);
	int32_t removeUnclusteredQueries(vector<qcSigList_t*> &qcSigList_vec, vector<qcSigListVec_t*> &query_clu_vec);
	int32_t getCluIDByQuery(qcSigList_t *qcSigList_node, vector<qcSigListVec_t*> &query_clu_vec);
	qcSigList_t* extractQcSigsSingleQueryClipReg(vector<struct querySeqInfoNode*> &query_seqs);
	vector<int32_t> getLeftRightMostAlnSegIdx(vector<struct querySeqInfoNode*> &query_seqs);
	vector<int32_t> getAdjacentAlnSegClipReg(int32_t arr_idx, int32_t clip_end_flag, vector<struct querySeqInfoNode*> &query_seqs, int32_t minClipEndSize);
	struct seedQueryInfo* chooseSeedClusterQueryClipReg(qcSigList_t* qcSigList_node, vector<qcSigList_t*> &q_cluster);
	double computeMatchRatioClipReg(qcSigList_t* query_seq_info_node, qcSigList_t* q_cluster_node, int64_t startSpanPos, int64_t endSpanPos);
	vector<int8_t> computeQcMatchProfileSingleQueryClipReg(qcSigList_t *queryCluSig, qcSigList_t *seed_qcQuery);
	bool isQcSigMatchClipReg(qcSig_t *qc_sig, qcSig_t *seed_qc_sig, int64_t max_dist_match_clip_pos, double min_size_ratio_match_thres_clip, double size_ratio_match_thres);
	vector<int8_t> qComputeSigMatchProfileClipReg(struct alnScoreNode *scoreArr, int32_t rowsNum, int32_t colsNum, qcSigList_t *queryCluSig, qcSigList_t *seed_qcQuery);
};

#endif /* SRC_CLIPREGCLUSTER_H_ */
