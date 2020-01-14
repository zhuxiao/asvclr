#include "blatAlnTra.h"
#include "util.h"

blatAlnTra::blatAlnTra(string &alnfilename, string &contigfilename, string &refseqfilename) {
	this->alnfilename = alnfilename;
	this->ctgfilename = contigfilename;
	this->refseqfilename = refseqfilename;
}

blatAlnTra::~blatAlnTra() {
}

void blatAlnTra::generateBlatResult(){
	if(isFileExist(ctgfilename) and isFileExist(refseqfilename) and !isFileExist(alnfilename))
		blatAln(alnfilename, ctgfilename, refseqfilename);
}
