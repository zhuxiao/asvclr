#include "RefSeqLoader.h"

pthread_mutex_t mutex_fai = PTHREAD_MUTEX_INITIALIZER;

RefSeqLoader::RefSeqLoader(string &reg, faidx_t *fai) {
	this->reg = reg;
	this->fai = fai;
	this->refseq = NULL;
	this->refseq_len = 0;
}

RefSeqLoader::~RefSeqLoader() {
	free(refseq);
	refseq_len = 0;
}

void RefSeqLoader::getRefSeq(){
	pthread_mutex_lock(&mutex_fai);
	refseq = fai_fetch(fai, reg.c_str(), &refseq_len);
	pthread_mutex_unlock(&mutex_fai);
	if ( refseq_len < 0 ) {
		cerr << __func__ << ": failed to fetch sequence in " << reg << endl;
		exit(1);
	}
}
