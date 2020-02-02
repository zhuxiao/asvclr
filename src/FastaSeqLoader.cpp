#include "FastaSeqLoader.h"
#include "util.h"

FastaSeqLoader::FastaSeqLoader(string &fastafilename) {
	this->fastafilename = fastafilename;
	initFastaSeq();
}

FastaSeqLoader::~FastaSeqLoader() {
	// TODO Auto-generated destructor stub
}

void FastaSeqLoader::initFastaSeq(){
	string fastaSeq, line;
	ifstream infile(fastafilename);
	if(!infile.is_open()) {
		cerr << "FastaSeqLoader: Cannot open file " << fastafilename << endl;
		exit(1);
	}

	fastaSeq = "";
	while(getline(infile, line))
		if(line.size()>0){
			if(line[0]=='>'){ // header line
				fastaSeqNameVec.push_back(line.substr(1)); // save header name

				if(fastaSeq.size()>0){
					fastaSeqVec.push_back(fastaSeq); // save sequence
					fastaSeq = "";
				}
			}else { // seq
				fastaSeq += line;
			}
		}
	if(fastaSeq.size()>0) fastaSeqVec.push_back(fastaSeq); // save sequence

	infile.close();
}

// get the query sequence
string FastaSeqLoader::getFastaSeq(size_t fa_id, size_t aln_orient){
	string ret;
	if(fa_id<0 or fa_id>=fastaSeqVec.size()) {
		cerr << "FastaSeqLoader: invalid ctg_idx=" << fa_id << endl;
		exit(1);
	}
	ret = fastaSeqVec[fa_id];
	if(aln_orient==ALN_MINUS_ORIENT) reverseComplement(ret);
	return ret;
}

// get the contig sequence by pos
string FastaSeqLoader::getFastaSeqByPos(size_t fa_id, size_t startPos, size_t endPos, size_t aln_orient){
	size_t tmp_startPos, tmp_endPos, tmp_startPos2, tmp_endPos2;
	string ret;

	if(fa_id<0 or fa_id>=fastaSeqVec.size()) {
		cerr << "FastaSeqLoader: invalid ctg_idx=" << fa_id << endl;
		exit(1);
	}
	if((startPos<1 or startPos>fastaSeqVec[fa_id].size()) or (endPos<1 or endPos>fastaSeqVec[fa_id].size())){
		cerr << "FastaSeqLoader: invalid region with startPos=" << startPos << ", endPos=" << endPos << ", filename=" << fastafilename << endl;
		exit(1);
	}

	if(startPos<=endPos){
		tmp_startPos = startPos;
		tmp_endPos = endPos;
	}else{
		tmp_startPos = endPos;
		tmp_endPos = startPos;
	}

	if(aln_orient==ALN_PLUS_ORIENT){
		ret = fastaSeqVec[fa_id].substr(tmp_startPos-1, tmp_endPos-tmp_startPos+1);
	}else{
		tmp_startPos2 = fastaSeqVec[fa_id].size() - tmp_endPos + 1;
		tmp_endPos2 = fastaSeqVec[fa_id].size() - tmp_startPos + 1;
		ret = fastaSeqVec[fa_id].substr(tmp_startPos2-1, tmp_endPos2-tmp_startPos2+1);
		reverseComplement(ret);
	}

	return ret;
}

size_t FastaSeqLoader::getFastaSeqLen(size_t fa_id){
	return fastaSeqVec[fa_id].size();
}

size_t FastaSeqLoader::getFastaSeqCount(){
	return fastaSeqNameVec.size();
}

vector<string> FastaSeqLoader::getFastaSeqNames(){
	return fastaSeqNameVec;
}

string FastaSeqLoader::getFastaSeqNameByID(int32_t fa_id){
	if(fa_id<0 or fa_id+1>(int32_t)fastaSeqNameVec.size()){
		cerr << "FastaSeqLoader: invalid fa_id=" << fa_id << endl;
		exit(1);
	}
	return fastaSeqNameVec.at(fa_id);
}
