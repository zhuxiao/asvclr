#ifndef BLATALNTRA_H_
#define BLATALNTRA_H_

#include <iostream>
#include <string>

using namespace std;

class blatAlnTra {
	public:
		string alnfilename, ctgfilename, refseqfilename;

	public:
		blatAlnTra(string &alnfilename, string &contigfilename, string &refseqfilename);
		virtual ~blatAlnTra();
		void generateBlatResult();
};

#endif /* BLATALNTRA_H_ */
