#ifndef SRC_INDELREGCLUSTER_H_
#define SRC_INDELREGCLUSTER_H_

#include "structures.h"
#include "util.h"

using namespace std;

#define POA_ALIGN_DEBUG					0

#define EXT_SIZE_ENTIRE_FLANKING		50

#define MIN_SINGLE_CLUSTER_GROUP_RATIO_THRES		(0.4f)
#define MIN_INVLAID_SECOND_GROUP_RATIO_THRES		(0.8f)

#define CNS_MIN_IDENTITY_CLUSTER_THRES				(0.8f)
#define CNS_MIN_SIZE_RATIO_MR0_THRES				(0.7f)

#define UNUSED_MODE			0
#define MIX_MODE			1
#define INS_ONLY_MODE		2
#define DEL_ONLY_MODE		3


class indelRegCluster {
private:
	string chrname;
	int64_t var_startRefPos, var_endRefPos, chrlen, startRefPos_cns, endRefPos_cns;
	int32_t sv_len_est, min_sv_size, min_supp_num, max_merge_span;
	double min_identity_match;
	faidx_t *fai;

public:
	indelRegCluster(string &chrname, int64_t var_startRefPos, int64_t var_endRefPos, int32_t sv_len_est, int64_t startRefPos_cns, int64_t endRefPos_cns, int32_t min_sv_size, int32_t min_supp_num, double min_identity_match, faidx_t *fai);
	virtual ~indelRegCluster();
	vector<struct querySeqInfoVec*> queryCluster(vector<struct querySeqInfoNode*> &query_seq_info_all);

private:
	void extractQcSigsFromAlnSegs(struct querySeqInfoNode* &query_seq_info_node);

	// void mergeNeighbouringSigs(vector<struct querySeqInfoNode*> &query_seq_info_vec, int32_t max_ref_dist_thres, double min_merge_identity_thres, faidx_t *fai);
	bool mergeNeighbouringSigsFlag(struct querySeqInfoNode* query_seq_info_node, vector<qcSig_t*> &sig_vec, int32_t min_ref_dist_thres, int32_t max_ref_dist_thres, double min_merge_identity_thres, double min_valid_sig_size_ratio_thres, faidx_t *fai);
	void filterSmallSigs(vector<qcSig_t *> &qcSig_vec, double valid_size_ratio_thres);
	void computeVarSums(struct querySeqInfoNode* &query_seq_info_node, int64_t start_var_pos, int64_t end_var_pos);
	void prepareQueryInfoForCluster(vector<struct querySeqInfoNode*> &query_seq_info_vec);
	void sortQueryInfoByAverSigSize(vector<struct querySeqInfoNode*> &query_seq_info_vec);
	void sortQueryInfoByNumCategory(vector<struct querySeqInfoNode*> &query_seq_info_vec);
	
	int32_t estimateClusterGroupNum(vector<struct querySeqInfoNode*> &query_seq_info_vec);
	vector<struct querySeqInfoVec*> queryClusterSingleGroup(vector<struct querySeqInfoNode*> &query_seq_info_all);
	void queryClusterAppendSingleGroup(vector<struct querySeqInfoNode*> &q_cluster_a, vector<struct querySeqInfoNode*> &query_seq_info_all);
	vector<struct querySeqInfoVec*> queryClusterDoubleGroup(vector<struct querySeqInfoNode*> &query_seq_info_all);

	struct seedQueryInfo* chooseSeedClusterQuery(struct querySeqInfoNode* query_seq_info_node, vector<struct querySeqInfoNode*> &q_cluster);
	//struct seedQueryInfo* chooseSeedClusterQuery2(struct querySeqInfoNode* query_seq_info_node, vector<struct querySeqInfoNode*> &q_cluster);
	struct seedQueryInfo* chooseSeedClusterQueryRecluster(struct seedQueryInfo* seedinfo, vector<struct querySeqInfoNode*> &q_cluster);

	vector<int32_t> getInitGroupIdx(vector<mrmin_match_t *> &mr_min_match_vec, vector<struct querySeqInfoNode*> &qseq_info_all);

	struct seedQueryInfo* getInitClusterItem(vector<int32_t> &qseq_id_vec_sub_group, vector<struct querySeqInfoNode*> &qseq_info_all);
	vector<struct querySeqInfoNode*> getInitClusterItemsMR0(vector<int32_t> &qseq_id_vec_sub_group, vector<struct querySeqInfoNode*> &qseq_info_all);

	double computeMatchRatio(struct querySeqInfoNode *query_seq_info_node, struct querySeqInfoNode *q_cluster_node, int64_t startSpanPos, int64_t endSpanPos, int32_t min_sv_size, double min_identity_match);
	double computeMatchRatioBySVlen(struct querySeqInfoNode* query_seq_info_node, struct querySeqInfoNode* q_cluster_node, int64_t startSpanPos, int64_t endSpanPos, int32_t min_sv_size);
	vector<int8_t> computeQcMatchProfileSingleQuery(queryCluSig_t *queryCluSig, queryCluSig_t *seed_qcQuery, bool mergeFlag, double min_identity_match);
	vector<int8_t> qComputeSigMatchProfile(struct alnScoreNode *scoreArr, int32_t rowsNum, int32_t colsNum, queryCluSig_t *queryCluSig, queryCluSig_t *seed_qcQuery);
	bool getReClusterFlag(vector<mrmin_match_t *> &mr_min_match_vec);
	bool initReCluster(vector<mrmin_match_t *> &mr_vec_recluster, vector<int32_t> &min_id_for_recluster_vec);
	vector<mrmin_match_t *> reClusterByMR(vector<struct querySeqInfoNode*> &query_seq_info_all, vector<int32_t> &min_id_vec_recluster, bool mrmin_match_vec_begin_flag);
	void sortMRVec(vector<mrmin_match_t *> &match_vec);
	void printMRVec(vector<mrmin_match_t *> &match_vec);
	bool matchVecValidFlag(vector<mrmin_match_t *> &match_vec_tmp, bool match_vec_begin_flag);

	// debug operation
	void printQcSigIndelCluster(vector<struct querySeqInfoVec*> &query_clu_vec);
	void printQcSigIndel(vector<struct querySeqInfoNode*> &query_seq_info_all);
	void printQcSigIndelSingleQuery(struct querySeqInfoNode *query_seq_info_node, int32_t idx);

};


#endif
