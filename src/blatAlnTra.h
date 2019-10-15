#ifndef BLATALNTRA_H_
#define BLATALNTRA_H_

#include <iostream>
#include <string>

using namespace std;

extern void blatAln(string &alnfilename, string &contigfilename, string &refseqfilename);

class blatAlnTra {
	public:
		string alnfilename, ctgfilename, refseqfilename;

	public:
		//blatAlnTra();
		blatAlnTra(string &alnfilename, string &contigfilename, string &refseqfilename);
		virtual ~blatAlnTra();
		//void blatAlnTra_mt(vector<blatAlnTra*> *blat_aln_tra_vec);
		void generateBlatResult();
};

#endif /* BLATALNTRA_H_ */
