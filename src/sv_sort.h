#ifndef SRC_SV_SORT_H_
#define SRC_SV_SORT_H_

#include <iostream>
#include <string.h>
#include <vector>
#include <set>
#include <fstream>
#include <algorithm>

using namespace std;

#define VAR_UNC				0	// uncertain
#define VAR_INS				1	// insertion
#define VAR_DEL				2	// deletion
#define VAR_DUP				3	// duplication
#define VAR_INV				4	// inversion
#define VAR_TRA				5	// translocation
#define VAR_BND				6	// BND
#define VAR_INV_TRA			7	// inverted translocation
#define VAR_MIX				10	// mixed variation
#define VAR_MNP				11	// MNP
#define VAR_SNV				12	// SNV
#define VAR_CNV				13	// CNV

typedef struct{
	string chrname, chrname2, info;
	int64_t startPos, endPos, startPos2, endPos2:56, sv_type:8;
}SV_item;


SV_item *constructSVItem(string &line);
vector<SV_item*> loadDataBED(string &filename);
bool isComma(string &seq);
bool isSeq(string &seq);
bool isBase(const char ch);
vector<string> getSVType(vector<string> &str_vec);
SV_item *allocateSVItem(string &chrname, size_t startPos, size_t endPos, string &chrname2, size_t startPos2, size_t endPos2, string &sv_type_str, string &line);
vector<SV_item*> loadDataVcf(string &filename);
void destroyData(vector<SV_item*> &sv_vec);
set<string> getChrnames(vector<SV_item*> &dataset);
vector<string> sortChrnames(set<string> &chrname_set);
bool IsSameChrname(string &chrname1, string &chrname2);
vector<SV_item*> getItemsByChr(string &chrname, vector<SV_item*> &dataset);
vector<vector<SV_item*>> constructSubsetByChrOp(vector<SV_item*> &sv_vec, vector<string> &chrname_vec);
vector<vector<SV_item*>> constructSubsetByChr(vector<SV_item*> &sv_vec);
bool sortFunSameChr(const SV_item *item1, const SV_item *item2);
void sortSVitem(vector<vector<SV_item*>> &subsets);
void sortSubset(vector<SV_item*> &sv_vec);


#endif /* SRC_SV_SORT_H_ */
