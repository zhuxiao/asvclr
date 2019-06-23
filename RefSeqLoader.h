#ifndef SRC_REFSEQLOADER_H_
#define SRC_REFSEQLOADER_H_

#include <iostream>
#include <string>

#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/faidx.h>

using namespace std;

class RefSeqLoader {
	public:
		string reg;
		faidx_t *fai;
		char *refseq;
		int refseq_len;

	public:
		RefSeqLoader(string &reg, faidx_t *fai);
		virtual ~RefSeqLoader();
		void getRefSeq();

};

#endif /* SRC_REFSEQLOADER_H_ */
