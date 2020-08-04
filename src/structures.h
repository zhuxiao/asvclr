#ifndef SRC_STRUCTURES_H_
#define SRC_STRUCTURES_H_

#include <iostream>
#include <string>
#include <vector>
#include <htslib/sam.h>
#include <htslib/faidx.h>

using namespace std;

// from clipAlnDataLoader.h
typedef struct{
	bam1_t *bam;
	string queryname, chrname;
	size_t querylen, aln_orient, startRefPos, endRefPos, startQueryPos, endQueryPos, leftClipSize, rightClipSize;
	bool leftHardClippedFlag, rightHardClippedFlag;
	bool left_clip_checked_flag, right_clip_checked_flag, query_checked_flag, SA_tag_flag;
}clipAlnData_t;

// from Region.h
typedef struct{
	string chrname;
	int64_t startRefPos, endRefPos, startLocalRefPos, endLocalRefPos, startQueryPos, endQueryPos;
	int32_t var_type, aln_orient, dup_num;
	int32_t query_id, sv_len, blat_aln_id;
	int32_t leftExtGapSize, rightExtGapSize;		// used for extended gap size on both sides of the region according to alignments
	string refseq, altseq;
	bool call_success_status, short_sv_flag, zero_cov_flag;
}reg_t;

// from clipReg.h
typedef struct{
	string chrname;
	size_t clipRefPos, clipLocalRefPos, clipQueryPos, aln_orient;
	bool same_orient_flag;  // true: ++, --; false: -+, +-
}clipPos_t;


//from covLoader.h
struct alnSeg{
	size_t startRpos, startQpos;
	size_t seglen: 26, opflag : 6;
	string seg_MD;
};
//from covLoader.h
struct MD_seg{
	string seg;
	size_t seglen: 26, opflag: 6;
};

// Genome.h
typedef struct {
	string refseqfilename, ctgfilename, alnfilename;
}blat_aln_filename_t;


// from Paras.h
typedef struct {
	string chrname;
	int64_t startPos, endPos;	// -1 for CHR format
}simpleReg_t;

// from Paras.h
typedef struct {
	string chrname, readsfilename, contigfilename, refseqfilename, tmpdir;
	reg_t **var_array;
	uint32_t arr_size;
	bool clip_reg_flag;
}assembleWork_opt;

// from Paras.h
typedef struct {
	assembleWork_opt *assem_work_opt;

	size_t work_id, num_work, num_work_per_ten_percent;  // 'work_id' starts from 1
	size_t *p_assemble_reg_workDone_num;   // pointer to the global variable which was declared in Paras.h
	pthread_mutex_t *p_mtx_assemble_reg_workDone_num; // pointer to the global variable which was declared in Paras.h
	size_t num_threads_per_assem_work;

	string inBamFile;
	faidx_t *fai;
	ofstream *var_cand_file;
	double expected_cov_assemble;
	bool delete_reads_flag;
}assembleWork;

// from varCand.h
typedef struct{
	int64_t query_start, query_end, subject_start, subject_end;
	float ident_percent;
	int32_t aln_orient;
	int64_t ref_start, ref_end;  // reference positions
}aln_seg_t;

// from varCand.h
typedef struct{
	//string query_name, subject_name;
	size_t blat_aln_id, query_len, subject_len, query_id, aln_orient;
	bool head_hanging, tail_hanging, best_aln, valid_aln;  // default: false
	vector<aln_seg_t*> aln_segs;
} blat_aln_t;

// from varCand.h
typedef struct{
	reg_t *reg, *cand_reg;
	aln_seg_t *aln_seg;
	int32_t startRefPos, endRefPos, startLocalRefPos, endLocalRefPos, startQueryPos, endQueryPos;
	string ctgseq, refseq;
	vector<string> alignResultVec;		// [0]: aligned ctgseq, [1]: match character, [2]: aligned refseq
	int32_t overlapLen, queryLeftShiftLen, queryRightShiftLen, localRefLeftShiftLen, localRefRightShiftLen;
	int32_t start_aln_idx_var, end_aln_idx_var;
}localAln_t;



#endif /* SRC_STRUCTURES_H_ */
