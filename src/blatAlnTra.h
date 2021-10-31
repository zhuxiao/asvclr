#ifndef BLATALNTRA_H_
#define BLATALNTRA_H_

#include <iostream>
#include <string>

#include "structures.h"

using namespace std;

class blatAlnTra {
	public:
		string alnfilename, ctgfilename, refseqfilename;
		int32_t max_proc_running_minutes;
		bool killed_flag;

		// process monitor killed blat work
		vector<killedBlatWork_t*> *killed_blat_work_vec;
		ofstream *killed_blat_work_file;
		pthread_mutex_t *mtx_killed_blat_work;

	public:
		blatAlnTra(string &alnfilename, string &contigfilename, string &refseqfilename, int32_t max_proc_running_minutes);
		virtual ~blatAlnTra();
		void generateBlatResult();
		bool getBlatWorkKilledFlag();
};

#endif /* BLATALNTRA_H_ */
