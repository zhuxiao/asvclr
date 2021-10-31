#include "blatAlnTra.h"
#include "util.h"

blatAlnTra::blatAlnTra(string &alnfilename, string &contigfilename, string &refseqfilename, int32_t max_proc_running_minutes) {
	this->alnfilename = alnfilename;
	this->ctgfilename = contigfilename;
	this->refseqfilename = refseqfilename;
	this->max_proc_running_minutes = max_proc_running_minutes;
	killed_blat_work_vec = NULL;
	killed_blat_work_file = NULL;
	mtx_killed_blat_work = NULL;
	killed_flag = false;
}

blatAlnTra::~blatAlnTra() {
}

void blatAlnTra::generateBlatResult(){
	int32_t ret;
	killedBlatWork_t *killed_blat_work;

	if(killed_flag==false) killed_flag = getBlatWorkKilledFlag();

	if(isFileExist(ctgfilename) and isFileExist(refseqfilename) and !isFileExist(alnfilename)){
		if(killed_flag==false){
			ret = blatAln(alnfilename, ctgfilename, refseqfilename, max_proc_running_minutes);
			if(ret==-2){ // killed work, then add it to vector and record it to file
				if(killed_blat_work_vec and killed_blat_work_file and mtx_killed_blat_work){ // pointer should not be NULL
					killed_flag = true;

					killed_blat_work = new killedBlatWork_t();
					killed_blat_work->alnfilename = alnfilename;
					killed_blat_work->ctgfilename = ctgfilename;
					killed_blat_work->refseqfilename = refseqfilename;

					pthread_mutex_lock(mtx_killed_blat_work);
					(*killed_blat_work_vec).push_back(killed_blat_work);
					*killed_blat_work_file << alnfilename << "\t" << ctgfilename << "\t" << refseqfilename << endl;
					pthread_mutex_unlock(mtx_killed_blat_work);
				}
			}
		}
	}
}

bool blatAlnTra::getBlatWorkKilledFlag(){
	killedBlatWork_t* killed_blat_work;
	bool killed_flag = false;

	if(killed_blat_work_vec and killed_blat_work_file and mtx_killed_blat_work){ // pointer should not be NULL
		for(size_t i=0; i<killed_blat_work_vec->size(); i++){
			killed_blat_work = killed_blat_work_vec->at(i);
			if(alnfilename.compare(killed_blat_work->alnfilename)==0 and ctgfilename.compare(killed_blat_work->ctgfilename)==0 and refseqfilename.compare(killed_blat_work->refseqfilename)==0){
				killed_flag = true;
				break;
			}
		}
	}
	return killed_flag;
}
