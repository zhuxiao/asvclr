#ifndef SRC_STRUCTURES_H_
#define SRC_STRUCTURES_H_

#include <iostream>
#include <string>
#include <vector>
#include <htslib/sam.h>
#include <htslib/faidx.h>

using namespace std;

// from clipAlnDataLoader.h
typedef struct clipAlnData_node{
	bam1_t *bam;
	string queryname, chrname;
	int64_t startRefPos, endRefPos:60, aln_orient:4;
	int32_t querylen, startQueryPos, endQueryPos, leftClipSize, rightClipSize;
	int32_t ref_dist, query_dist;
	bool leftHardClippedFlag, rightHardClippedFlag;
	bool left_clip_checked_flag, right_clip_checked_flag, query_checked_flag, SA_tag_flag;
	struct clipAlnData_node *left_aln, *right_aln;
}clipAlnData_t;

// from Region.h
typedef struct{
	string chrname;
	int64_t startRefPos, endRefPos;
	int32_t startLocalRefPos, endLocalRefPos, startQueryPos, endQueryPos, sv_len, dup_num; //, support_num, cov_num;
	//int64_t ref_len, local_ref_len, query_len;
	int8_t var_type, aln_orient;
	int16_t query_id, blat_aln_id;
	int32_t leftExtGapSize, rightExtGapSize;		// used for extended gap size on both sides of the region according to alignments
	string refseq, altseq;
	bool call_success_status, short_sv_flag, zero_cov_flag, aln_seg_end_flag;
}reg_t;

// from clipReg.h
typedef struct{
	string chrname;
	int64_t clipRefPos:60, clip_end:4;
	int32_t clipLocalRefPos, clipQueryPos:28, aln_orient:4;
	clipAlnData_t *clip_aln;
	bool same_orient_flag;  // true: ++, --; false: -+, +-
	bool self_overlap_flag;
}clipPos_t;


//from covLoader.h
struct alnSeg{
	int64_t startRpos, startQpos;
	int32_t seglen: 26, opflag : 6;
	string seg_MD;
};
//from covLoader.h
struct MD_seg{
	string seg;
	int32_t seglen: 26, opflag: 6;
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
	simpleReg_t **limit_reg_array;
	uint32_t arr_size, limit_reg_array_size;
	bool clip_reg_flag, limit_reg_process_flag;
}assembleWork_opt;

// from Paras.h
typedef struct {
	assembleWork_opt *assem_work_opt;

	int32_t work_id, num_work, num_work_per_ten_percent;  // 'work_id' starts from 1
	int32_t *p_assemble_reg_workDone_num;   // pointer to the global variable which was declared in Paras.h
	pthread_mutex_t *p_mtx_assemble_reg_workDone_num; // pointer to the global variable which was declared in Paras.h
	int32_t num_threads_per_assem_work, minClipEndSize;

	string inBamFile, technology, canu_version;
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
	int32_t blat_aln_id, query_len, subject_len;
	int32_t query_id:24, aln_orient:8;
	bool head_hanging, tail_hanging, best_aln, valid_aln;  // default: false
	vector<aln_seg_t*> aln_segs;
} blat_aln_t;

// from varCand.h
typedef struct{
	reg_t *reg, *cand_reg;
	aln_seg_t *aln_seg, *start_seg_extend, *end_seg_extend;  // for short variants
	int64_t chrlen, startRefPos, endRefPos;
	int32_t blat_aln_id, startLocalRefPos, endLocalRefPos, startQueryPos, endQueryPos;
	string ctgseq, refseq;
	vector<string> alignResultVec;		// [0]: aligned ctgseq, [1]: match character, [2]: aligned refseq
	int32_t overlapLen, queryLeftShiftLen, queryRightShiftLen, localRefLeftShiftLen, localRefRightShiftLen;
	int32_t start_aln_idx_var, end_aln_idx_var;
}localAln_t;

typedef struct{
	int32_t start_aln_idx, end_aln_idx, reg_size;
	int64_t startRefPos, endRefPos, startLocalRefPos, endLocalRefPos, startQueryPos, endQueryPos;
	bool gap_flag, valid_flag;
}mismatchReg_t;

// from Genome.h
typedef struct BND_Node{
	string bnd_id;
	int8_t reg_id, clip_end, mate_reg_id, mate_clip_end;
	char mate_orient_ch;
	string chrname, mate_chrname;
	int64_t bnd_pos, mate_bnd_pos;
	string seq, bnd_str, vcf_line;
	int32_t support_num, mate_support_num, cov_num, mate_cov_num;
	struct BND_Node *mate_node;
}BND_t;

#endif /* SRC_STRUCTURES_H_ */
