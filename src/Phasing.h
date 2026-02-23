#ifndef SRC_PHASING_H_
#define SRC_PHASING_H_

#include <map>

#include "structures.h"
#include "varCand.h"

#define PHASING_DEBUG			0

#define MIN_LINK_READS_NUM		3

class Phasing {
public:
	int32_t min_sv_size, max_sv_size;
	double size_ratio_thres, seqsim_thres;

private:
	vector<varCand*> var_cand_vec, var_cand_clipReg_vec;

public:
	Phasing(int32_t min_sv_size, int32_t max_sv_size, double size_ratio_thres, double seqsim_thres, vector<varCand*> &var_cand_vec, vector<varCand*> &var_cand_clipReg_vec);
	virtual ~Phasing();
	void performPhasing();

private:
	void releaseHapNodes(vector< vector<phasing_reg_t*> > &hapnodes_vec);
	void generateHapNodeSingleVarCand(vector< vector<phasing_reg_t*> > &hapnodes_vec, varCand *var_cand);
	bool linkHapNodes(vector< vector<phasing_reg_t*> > &hapnodes_vec, vector< vector<phasing_reg_t*> > &hapnodes_vec2);
	int32_t computeLinkNum(vector<string> &qname_vec, vector<string> &qname_vec2);
	vector<string> getReadNamesFromCnsHeader(string &cns_header_name);
	phasing_reg_t* dupPhasingRegItem(phasing_reg_t *phasing_reg);
	bool isValidStartLinkHapNode(vector< vector<phasing_reg_t*> > &hapnodes_vec);
	void printHapNodes(vector< vector<phasing_reg_t*> > &hapnodes_vec);
};


#endif /* SRC_PHASING_H_ */
