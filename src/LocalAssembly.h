#ifndef SRC_LOCALASSEMBLY_H_
#define SRC_LOCALASSEMBLY_H_

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <limits.h>
#include <sys/stat.h>

#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/faidx.h>

#include "Region.h"
//#include "alnDataLoader.h"
#include "clipAlnDataLoader.h"
#include "RefSeqLoader.h"
#include "util.h"

using namespace std;

#define REFSEQ_SIDE_LEN					10000
#define ASSEMBLY_GENOME_SIZE_INITIAL	30000
#define ASSEMBLY_STEP_SIZE				10000

class LocalAssembly {
	public:
		string chrname, readsfilename, contigfilename, refseqfilename, tmpdir, inBamFile;
		size_t chrlen;
		vector<reg_t*> varVec;
		faidx_t *fai;

	private:
		vector<bam1_t*> alnDataVector;
		vector<clipAlnData_t*> clipAlnDataVector;

	public:
		LocalAssembly(string &readsfilename, string &contigfilename, string &refseqfilename, string &tmpdir, vector<reg_t*> &varVec, string &chrname, string &inBamFile, faidx_t *fai);
		virtual ~LocalAssembly();
		void extractRefseq();
		void extractReadsDataFromBAM();
		bool localAssembleCanu();
		void recordAssemblyInfo(ofstream &assembly_info_file);

	private:
		void destoryAlnData();
		void destoryClipAlnData();
		vector<clipAlnData_t*> getQueryClipAlnSegs(string &queryname, vector<clipAlnData_t*> &clipAlnDataVector);
		int32_t getNoHardClipAlnItem(vector<clipAlnData_t*> &clipAlnDataVector);
		void markHardClipSegs(size_t idx, vector<clipAlnData_t*> &query_aln_segs);
		vector<string> getQuerySeqWithSoftClipSeqs(clipAlnData_t* clip_aln);
		vector<string> getQuerySeqWithoutSoftClipSeqs(clipAlnData_t* clip_aln);
		bool localAssembleCanu_IncreaseGenomeSize();
		bool localAssembleCanu_DecreaseGenomeSize();
		vector<string> joinQueryAlnSegs(vector<clipAlnData_t*> &query_aln_segs);
		int32_t getLeftMostAlnSeg(vector<clipAlnData_t*> &query_aln_segs);
};

#endif /* SRC_LOCALASSEMBLY_H_ */
