#ifndef SRC_FASTASEQLOADER_H_
#define SRC_FASTASEQLOADER_H_

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sys/stat.h>

using namespace std;

class FastaSeqLoader {
	public:
		string fastafilename;
		vector<string> fastaSeqNameVec;
		vector<string> fastaSeqVec;

	public:
		FastaSeqLoader();
		FastaSeqLoader(string &fastafilename);
		virtual ~FastaSeqLoader();
		string getFastaSeq(size_t fa_id);
		string getFastaSeq(size_t fa_id, size_t aln_orient);
		string getFastaSeqByPos(size_t fa_id, size_t startPos, size_t endPos, size_t aln_orient); // ctg_id starts from 0, pos starts from 1
		size_t getFastaSeqLen(size_t fa_id);
		size_t getFastaSeqCount();
		vector<string> getFastaSeqNames();
		string getFastaSeqNameByID(int32_t fa_id);
		bool isValidSeq(string &fastafilename);
		bool isValidSeq();
		vector< vector<string> > getReadNamesFromCnsHeaders();
		vector< vector<string> > getReadNamesFromCnsHeaders(vector<string> &cns_header_vec);
		vector<string> getReadNamesFromSingleCnsHeader(string &cns_header);

	private:
		void initFastaSeq();
};

#endif /* SRC_FASTASEQLOADER_H_ */
