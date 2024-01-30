#include <iostream>
#include <cstring>
#include <unistd.h>
#include <pthread.h>
#include <cerrno>
#include <htslib/thread_pool.h>

#include "covLoader.h"
#include "util.h"
#include "Block.h"
#include "LocalCns.h"

// string split function
vector<string> split(const string& s, const string& delim)
{
	vector<string> elems;
	string tmp;
	size_t pos = 0;
	size_t len = s.length();
	size_t delim_len = delim.length();
	int find_pos;

	if (delim_len == 0) return elems;
	while (pos < len)
	{
		find_pos = s.find(delim, pos);
		if(find_pos < 0){
			elems.push_back(s.substr(pos, len - pos));
			break;
		}
		tmp = s.substr(pos, find_pos - pos);
		if(!tmp.empty()) elems.push_back(tmp);
		pos = find_pos + delim_len;
	}
	return elems;
}

// determine whether the str exist in string vector
bool isExistStr(string &str, vector<string> &vec){
	bool exist_flag = false;
	for(size_t i=0; i<vec.size(); i++)
		if(str.compare(vec.at(i))==0){
			exist_flag = true;
			break;
		}
	return exist_flag;
}

// copy single file
int copySingleFile(string &infilename, ofstream &outfile){
	ifstream infile;
	string line;

	if(isFileExist(infilename)==false) return 0;

	infile.open(infilename);
	if(!infile.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << infilename << endl;
		exit(1);
	}
	// read each line and save to the output file
	while(getline(infile, line))
		if(line.size()>0 and line.at(0)!='#'){
			outfile << line << endl;
		}
	infile.close();

	return 0;
}

char getBase(string &seq, size_t pos, size_t orient){
	char ch;

	if(pos<1 or pos>seq.size()){
		cerr << __func__ << ", line=" << __LINE__ << ": invalid pos: " << pos << endl;
		exit(1);
	}

	ch = seq[pos-1];
	if(orient==ALN_MINUS_ORIENT){
		switch(ch){
			case 'A': ch = 'T'; break;
			case 'C': ch = 'G'; break;
			case 'G': ch = 'C'; break;
			case 'T': ch = 'A'; break;
			case 'a': ch = 't'; break;
			case 'c': ch = 'g'; break;
			case 'g': ch = 'c'; break;
			case 't': ch = 'a'; break;
			case 'N':
			case 'n': break;
			default: ch = 'N'; break;
		}
	}

	return ch;
}

// reverse the sequence
void reverseSeq(string &seq){
	char tmp;
	int32_t len = seq.size();
	for(int32_t i=0; i<len/2; i++){
		tmp = seq[i];
		seq[i] = seq[len-1-i];
		seq[len-1-i] = tmp;
	}
}

// reverse complement the sequence
void reverseComplement(string &seq){
	reverseSeq(seq);

	int32_t len = seq.size();
	for(int32_t i=0; i<len; i++){
		switch(seq[i]){
			case 'a': seq[i] = 't'; break;
			case 'c': seq[i] = 'g'; break;
			case 'g': seq[i] = 'c'; break;
			case 't': seq[i] = 'a'; break;
			case 'A': seq[i] = 'T'; break;
			case 'C': seq[i] = 'G'; break;
			case 'G': seq[i] = 'C'; break;
			case 'T': seq[i] = 'A'; break;
			case 'N': seq[i] = 'N'; break;
			case 'n': seq[i] = 'n'; break;

			// mixed bases
			case 'M': seq[i] = 'K'; break;
			case 'm': seq[i] = 'k'; break;
			case 'R': seq[i] = 'Y'; break;
			case 'r': seq[i] = 'y'; break;
			case 'S': seq[i] = 'S'; break;
			case 's': seq[i] = 's'; break;
			case 'V': seq[i] = 'B'; break;
			case 'v': seq[i] = 'b'; break;
			case 'W': seq[i] = 'W'; break;
			case 'w': seq[i] = 'w'; break;
			case 'Y': seq[i] = 'R'; break;
			case 'y': seq[i] = 'r'; break;
			case 'H': seq[i] = 'D'; break;
			case 'h': seq[i] = 'd'; break;
			case 'K': seq[i] = 'M'; break;
			case 'k': seq[i] = 'm'; break;
			case 'D': seq[i] = 'H'; break;
			case 'd': seq[i] = 'h'; break;
			case 'B': seq[i] = 'V'; break;
			case 'b': seq[i] = 'v'; break;
			default:
				cout << "i=" << i << ", seq_len=" << len << endl;
				cout << seq << endl;
				cerr << __func__ << ": unknown base: " << seq[i] << endl;
				exit(1);
		}
	}
}

// upper the sequence
void upperSeq(string &seq){
	for(size_t i=0; i<seq.size(); i++)
		if(seq[i]>='a' and seq[i]<='z') seq[i] -= 32;
}

// get the number of contigs
size_t getCtgCount(string &contigfilename){
	size_t ctg_num;
	string line;
	ifstream infile;
	bool flag = false;

	flag = isFileExist(contigfilename);
	ctg_num = 0;
	if(flag){
		infile.open(contigfilename);
		if(!infile.is_open()){
			cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << contigfilename << endl;
			exit(1);
		}

		while(getline(infile, line))
			if(line.size())
				if(line[0]=='>')	ctg_num ++;
		infile.close();
	}

	return ctg_num;
}

// duplicate the variant region
reg_t* dupVarReg(reg_t *reg){
	reg_t *reg_ret = NULL;

	if(reg){
		reg_ret = new reg_t();
		reg_ret->chrname = reg->chrname;
		reg_ret->startRefPos = reg->startRefPos;
		reg_ret->endRefPos = reg->endRefPos;
		reg_ret->startLocalRefPos = reg->startLocalRefPos;
		reg_ret->endLocalRefPos = reg->endLocalRefPos;
		reg_ret->startQueryPos = reg->startQueryPos;
		reg_ret->endQueryPos = reg->endQueryPos;
		reg_ret->var_type = reg->var_type;
		reg_ret->sv_len = reg->sv_len;
		reg_ret->aln_orient = reg->aln_orient;
		reg_ret->dup_num = reg->dup_num;
		reg_ret->query_id = reg->query_id;
		reg_ret->blat_aln_id = reg->blat_aln_id;
		reg_ret->minimap2_aln_id = reg->minimap2_aln_id;
		reg_ret->leftExtGapSize = reg->leftExtGapSize;
		reg_ret->rightExtGapSize = reg->rightExtGapSize;
		reg_ret->call_success_status = reg->call_success_status;
		reg_ret->short_sv_flag = reg->short_sv_flag;
		reg_ret->zero_cov_flag = reg->zero_cov_flag;
		reg_ret->aln_seg_end_flag = reg->aln_seg_end_flag;
		reg_ret->query_pos_invalid_flag = reg->query_pos_invalid_flag;
		reg_ret->large_indel_flag = reg->large_indel_flag;
		reg_ret->gt_type = reg->gt_type;
		reg_ret->gt_seq = reg->gt_seq;
		reg_ret->refseq = reg->refseq;
		reg_ret->altseq = reg->altseq;
		reg_ret->AF = reg->AF;
	}

	return reg_ret;
}

// find the vector item
reg_t* findVarvecItem(int32_t startRefPos, int32_t endRefPos, vector<reg_t*> &varVec){
	reg_t *reg, *target_reg = NULL;
	for(size_t i=0; i<varVec.size(); i++){
		reg = varVec[i];
		if((startRefPos>=reg->startRefPos and startRefPos<=reg->endRefPos) or (endRefPos>=reg->startRefPos and endRefPos<=reg->endRefPos)
			or (reg->startRefPos>=startRefPos and reg->startRefPos<=endRefPos) or (reg->endRefPos>=startRefPos and reg->endRefPos<=endRefPos)){
			// overlap
			target_reg = reg;
			break;
		}
	}
	return target_reg;
}

// find all the vector item
vector<reg_t*> findVarvecItemAll(int32_t startRefPos, int32_t endRefPos, vector<reg_t*> &varVec){
	reg_t *reg;
	vector<reg_t*> reg_vec_ret;
	for(size_t i=0; i<varVec.size(); i++){
		reg = varVec[i];
		if((startRefPos>=reg->startRefPos and startRefPos<=reg->endRefPos) or (endRefPos>=reg->startRefPos and endRefPos<=reg->endRefPos)
			or (reg->startRefPos>=startRefPos and reg->startRefPos<=endRefPos) or (reg->endRefPos>=startRefPos and reg->endRefPos<=endRefPos)){
			// overlap
			reg_vec_ret.push_back(reg);
		}
	}
	return reg_vec_ret;
}

// find the vector item according to extended margin size
reg_t* findVarvecItemExtSize(int32_t startRefPos, int32_t endRefPos, vector<reg_t*> &varVec, int32_t leftExtSize, int32_t rightExtSize){
	reg_t *reg, *target_reg = NULL;
	for(size_t i=0; i<varVec.size(); i++){
		reg = varVec[i];
		if((startRefPos+leftExtSize>=reg->startRefPos and startRefPos<=reg->endRefPos+rightExtSize) or (endRefPos+leftExtSize>=reg->startRefPos and endRefPos<=reg->endRefPos+rightExtSize)
			or (reg->startRefPos+leftExtSize>=startRefPos and reg->startRefPos<=endRefPos+rightExtSize) or (reg->endRefPos+leftExtSize>=startRefPos and reg->endRefPos<=endRefPos+rightExtSize)){
			// overlap
			target_reg = reg;
			break;
		}
	}
	return target_reg;
}

// find all the vector item according to extended margin size
vector<reg_t*> findVarvecItemAllExtSize(int32_t startRefPos, int32_t endRefPos, vector<reg_t*> &varVec, int32_t leftExtSize, int32_t rightExtSize){
	reg_t *reg;
	vector<reg_t*> reg_vec_ret;

	for(size_t i=0; i<varVec.size(); i++){
		reg = varVec[i];
		if((startRefPos+leftExtSize>=reg->startRefPos and startRefPos<=reg->endRefPos+rightExtSize) or (endRefPos+leftExtSize>=reg->startRefPos and endRefPos<=reg->endRefPos+rightExtSize)
			or (reg->startRefPos+leftExtSize>=startRefPos and reg->startRefPos<=endRefPos+rightExtSize) or (reg->endRefPos+leftExtSize>=startRefPos and reg->endRefPos<=endRefPos+rightExtSize)){
			// overlap
			reg_vec_ret.push_back(reg);
		}
	}
	return reg_vec_ret;
}

// get the vector index
int32_t getVectorIdx(reg_t *reg, vector<reg_t*> &varVec){
	int32_t idx = -1;
	for(size_t i=0; i<varVec.size(); i++)
		if(varVec[i]==reg){
			idx = i;
			break;
		}
	return idx;
}

reg_t* getOverlappedRegByCallFlag(reg_t *reg, vector<reg_t*> &varVec){
	reg_t *reg_ret = NULL;
	int32_t idx = getOverlappedRegIdxByCallFlag(reg, varVec);
	if(idx!=-1) reg_ret = varVec.at(idx);
	return reg_ret;
}

reg_t* getOverlappedReg(reg_t *reg, vector<reg_t*> &varVec){
	reg_t *reg_ret = NULL;
	int32_t idx = getOverlappedRegIdx(reg, varVec);
	if(idx!=-1) reg_ret = varVec.at(idx);
	return reg_ret;
}


// get overlapped region
int32_t getOverlappedRegIdx(reg_t *reg, vector<reg_t*> &varVec){
	reg_t *reg_tmp;
	bool flag;
	int32_t idx_ret = -1;
	for(size_t i=0; i<varVec.size(); i++){
		reg_tmp = varVec.at(i);
		flag = isOverlappedReg(reg, reg_tmp);
		if(flag){
			idx_ret = i;
			break;
		}
	}
	return idx_ret;
}

// get overlapped region
int32_t getOverlappedRegIdxByCallFlag(reg_t *reg, vector<reg_t*> &varVec){
	reg_t *reg_tmp;
	bool flag;
	int32_t idx_ret = -1;
	for(size_t i=0; i<varVec.size(); i++){
		reg_tmp = varVec.at(i);
		if(reg_tmp->call_success_status){
			flag = isOverlappedReg(reg, reg_tmp);
			if(flag){
				idx_ret = i;
				break;
			}
		}
	}
	return idx_ret;
}

// determine whether two regions are overlapped
bool isOverlappedReg(reg_t* reg1, reg_t* reg2){
	bool flag = false;
	if(reg1->chrname.compare(reg2->chrname)==0){
		if(isOverlappedPos(reg1->startRefPos, reg1->endRefPos, reg2->startRefPos, reg2->endRefPos))
			flag = true;
	}
	return flag;
}

// determine whether two regions are overlapped considering extend positions
bool isOverlappedRegExtSize(reg_t* reg1, reg_t* reg2, int32_t leftExtSize, int32_t rightExtSize){
	bool flag = false;
	if(reg1->chrname.compare(reg2->chrname)==0){
		if(isOverlappedPos(reg1->startRefPos-leftExtSize, reg1->endRefPos+rightExtSize, reg2->startRefPos-leftExtSize, reg2->endRefPos+rightExtSize))
			flag = true;
	}
	return flag;
}

// determine whether the given positions of two regions are overlapped
bool isOverlappedPos(size_t startPos1, size_t endPos1, size_t startPos2, size_t endPos2){
	bool flag = false;
//	if((startPos1>=startPos2 and startPos1<=endPos2)
//		or (endPos2>=startPos1 and endPos2<=endPos1)
//		or (startPos2>=startPos1 and startPos2<=endPos1)
//		or (endPos1>=startPos2 and endPos1<=endPos2))

	if((endPos2>=startPos1 and endPos2<=endPos1) or (endPos1>=startPos2 and endPos1<=endPos2))
		flag = true;
	return flag;
}

int32_t getOverlapSize(size_t startPos1, size_t endPos1, size_t startPos2, size_t endPos2){
	int32_t overlap_size;

	if(startPos1<=startPos2) overlap_size = endPos1 - startPos2 + 1;
	else overlap_size = endPos2 - startPos1 + 1;

	return overlap_size;
}

bool isAdjacent(size_t startPos1, size_t endPos1, size_t startPos2, size_t endPos2, int32_t dist_thres){
	bool flag;
	int32_t overlap_size;

	if(dist_thres<=0){
		cerr << __func__ << ", line=" << __LINE__ << ", invalid distance threshold: " << dist_thres << endl;
		exit(1);
	}

	overlap_size = getOverlapSize(startPos1, endPos1, startPos2, endPos2);
	//if(overlap_size>=-dist_thres and overlap_size<=dist_thres) flag = true;
	if(overlap_size>=-dist_thres) flag = true;
	else flag = false;

	return flag;
}

// determine whether the two mate chip regions are overlapped
bool isOverlappedMateClipReg(mateClipReg_t *mate_clip_reg1, mateClipReg_t *mate_clip_reg2){
	bool overlap_flag, overlap_flag1, overlap_flag2;
	reg_t *left_reg1, *right_reg1, *left_reg2, *right_reg2;
	size_t mate_end_flag, start_pos1, end_pos1, start_pos2, end_pos2;
	string chrname1, chrname2;

	overlap_flag = overlap_flag1 = overlap_flag2 = false;
	if(mate_clip_reg1->large_indel_flag==false and mate_clip_reg2->large_indel_flag==false){
		if(mate_clip_reg1->leftClipRegNum==1 and mate_clip_reg1->rightClipRegNum==1 and mate_clip_reg2->leftClipRegNum==1 and mate_clip_reg2->rightClipRegNum==1){
			left_reg1 = mate_clip_reg1->leftClipReg ? mate_clip_reg1->leftClipReg : mate_clip_reg1->leftClipReg2;
			right_reg1 = mate_clip_reg1->rightClipReg ? mate_clip_reg1->rightClipReg : mate_clip_reg1->rightClipReg2;
			left_reg2 = mate_clip_reg2->leftClipReg ? mate_clip_reg2->leftClipReg : mate_clip_reg2->leftClipReg2;
			right_reg2 = mate_clip_reg2->rightClipReg ? mate_clip_reg2->rightClipReg : mate_clip_reg2->rightClipReg2;

			if((left_reg1->chrname.compare(left_reg2->chrname)==0 and right_reg1->chrname.compare(right_reg2->chrname)==0 and isOverlappedReg(left_reg1, left_reg2) and isOverlappedReg(right_reg1, right_reg2))
				or (left_reg1->chrname.compare(right_reg2->chrname)==0 and right_reg1->chrname.compare(left_reg2->chrname)==0 and isOverlappedReg(left_reg1, right_reg2) and isOverlappedReg(right_reg1, left_reg2)))
				overlap_flag = true;

		}else if(mate_clip_reg1->leftClipRegNum==2 and mate_clip_reg1->rightClipRegNum==2 and mate_clip_reg2->leftClipRegNum==2 and mate_clip_reg2->rightClipRegNum==2){
			chrname1 = chrname2 = "";
			start_pos1 = end_pos1 = start_pos2 = end_pos2 = 0;

			if(mate_clip_reg1->leftClipReg->chrname.compare(mate_clip_reg1->leftClipReg2->chrname)==0 and mate_clip_reg1->rightClipReg->chrname.compare(mate_clip_reg1->rightClipReg2->chrname)==0
				and mate_clip_reg2->leftClipReg->chrname.compare(mate_clip_reg2->leftClipReg2->chrname)==0 and mate_clip_reg2->rightClipReg->chrname.compare(mate_clip_reg2->rightClipReg2->chrname)==0){
				// check left region
				chrname1 = mate_clip_reg1->leftClipReg->chrname;
				if(mate_clip_reg1->leftMeanClipPos<mate_clip_reg1->leftMeanClipPos2){
					start_pos1 = mate_clip_reg1->leftMeanClipPos;
					end_pos1 = mate_clip_reg1->leftMeanClipPos2;
				}else{
					start_pos1 = mate_clip_reg1->leftMeanClipPos2;
					end_pos1 = mate_clip_reg1->leftMeanClipPos;
				}

				chrname2 = mate_clip_reg2->leftClipReg->chrname;
				if(mate_clip_reg2->leftMeanClipPos<mate_clip_reg2->leftMeanClipPos2){
					start_pos2 = mate_clip_reg2->leftMeanClipPos;
					end_pos2 = mate_clip_reg2->leftMeanClipPos2;
				}else{
					start_pos2 = mate_clip_reg2->leftMeanClipPos2;
					end_pos2 = mate_clip_reg2->leftMeanClipPos;
				}

				if(chrname1.size()>0 and chrname1.compare(chrname2)==0) overlap_flag1 = isOverlappedPos(start_pos1, end_pos1, start_pos2, end_pos2); // same chrome
				else overlap_flag1 = false;

				// check right region
				chrname1 = mate_clip_reg1->rightClipReg->chrname;
				if(mate_clip_reg1->rightMeanClipPos<mate_clip_reg1->rightMeanClipPos2){
					start_pos1 = mate_clip_reg1->rightMeanClipPos;
					end_pos1 = mate_clip_reg1->rightMeanClipPos2;
				}else{
					start_pos1 = mate_clip_reg1->rightMeanClipPos2;
					end_pos1 = mate_clip_reg1->rightMeanClipPos;
				}

				chrname2 = mate_clip_reg2->rightClipReg->chrname;
				if(mate_clip_reg2->rightMeanClipPos<mate_clip_reg2->rightMeanClipPos2){
					start_pos2 = mate_clip_reg2->rightMeanClipPos;
					end_pos2 = mate_clip_reg2->rightMeanClipPos2;
				}else{
					start_pos2 = mate_clip_reg2->rightMeanClipPos2;
					end_pos2 = mate_clip_reg2->rightMeanClipPos;
				}

				if(chrname1.size()>0 and chrname1.compare(chrname2)==0) overlap_flag2 = isOverlappedPos(start_pos1, end_pos1, start_pos2, end_pos2); // same chrome
				else overlap_flag2 = false;

				if(overlap_flag1 and overlap_flag2) overlap_flag = true;
			}
		}else{
			if(mate_clip_reg1->leftClipRegNum==1 and mate_clip_reg1->rightClipRegNum==2){
				// check left region
				chrname1 = chrname2 = "";
				start_pos1 = end_pos1 = start_pos2 = end_pos2 = 0;

				if(mate_clip_reg1->leftClipReg){
					chrname1 = mate_clip_reg1->leftClipReg->chrname;
					start_pos1 = mate_clip_reg1->leftClipReg->startRefPos;
					end_pos1 = mate_clip_reg1->leftClipReg->endRefPos;
				}else{
					chrname1 = mate_clip_reg1->leftClipReg2->chrname;
					start_pos1 = mate_clip_reg1->leftClipReg2->startRefPos;
					end_pos1 = mate_clip_reg1->leftClipReg2->endRefPos;
				}

				// get the mate end
				mate_end_flag = 0;
				if(mate_clip_reg2->leftClipRegNum==1 and mate_clip_reg2->rightClipRegNum==2){
					mate_end_flag = 1;
				}else if(mate_clip_reg2->leftClipRegNum==2 and mate_clip_reg2->rightClipRegNum==1){
					mate_end_flag = 2;
				}

				if(mate_end_flag==1){
					if(mate_clip_reg2->leftClipReg){
						chrname2 = mate_clip_reg2->leftClipReg->chrname;
						start_pos2 = mate_clip_reg2->leftClipReg->startRefPos;
						end_pos2 = mate_clip_reg2->leftClipReg->endRefPos;
					}else{
						chrname2 = mate_clip_reg2->leftClipReg2->chrname;
						start_pos2 = mate_clip_reg2->leftClipReg2->startRefPos;
						end_pos2 = mate_clip_reg2->leftClipReg2->endRefPos;
					}
				}else if(mate_end_flag==2){
					if(mate_clip_reg2->rightClipReg){
						chrname2 = mate_clip_reg2->rightClipReg->chrname;
						start_pos2 = mate_clip_reg2->rightClipReg->startRefPos;
						end_pos2 = mate_clip_reg2->rightClipReg->endRefPos;
					}else{
						chrname2 = mate_clip_reg2->rightClipReg2->chrname;
						start_pos2 = mate_clip_reg2->rightClipReg2->startRefPos;
						end_pos2 = mate_clip_reg2->rightClipReg2->endRefPos;
					}
				}

				if(chrname1.size()>0 and chrname1.compare(chrname2)==0) overlap_flag1 = isOverlappedPos(start_pos1, end_pos1, start_pos2, end_pos2); // same chrome
				else overlap_flag1 = false;

				// check right regions
				chrname1 = chrname2 = "";
				start_pos1 = end_pos1 = start_pos2 = end_pos2 = 0;

				if(mate_end_flag==1) mate_end_flag = 2;
				else if(mate_end_flag==2) mate_end_flag = 1;

				if(mate_clip_reg1->rightClipReg->chrname.compare(mate_clip_reg1->rightClipReg2->chrname)==0){
					chrname1 = mate_clip_reg1->rightClipReg->chrname;
					if(mate_clip_reg1->rightMeanClipPos<mate_clip_reg1->rightMeanClipPos2){
						start_pos1 = mate_clip_reg1->rightMeanClipPos;
						end_pos1 = mate_clip_reg1->rightMeanClipPos2;
					}else{
						start_pos1 = mate_clip_reg1->rightMeanClipPos2;
						end_pos1 = mate_clip_reg1->rightMeanClipPos;
					}

					if(mate_end_flag==1){
						if(mate_clip_reg2->leftClipRegNum==1){
							if(mate_clip_reg2->leftClipReg){
								chrname2 = mate_clip_reg2->leftClipReg->chrname;
								start_pos2 = mate_clip_reg2->leftClipReg->startRefPos;
								end_pos2 = mate_clip_reg2->leftClipReg->endRefPos;
							}else{
								chrname2 = mate_clip_reg2->leftClipReg2->chrname;
								start_pos2 = mate_clip_reg2->leftClipReg2->startRefPos;
								end_pos2 = mate_clip_reg2->leftClipReg2->endRefPos;
							}
						}else if(mate_clip_reg2->leftClipRegNum==2){
							if(mate_clip_reg2->leftClipReg->chrname.compare(mate_clip_reg2->leftClipReg2->chrname)==0){
								chrname2 = mate_clip_reg2->leftClipReg->chrname;
								if(mate_clip_reg2->leftMeanClipPos<mate_clip_reg2->leftMeanClipPos2){
									start_pos2 = mate_clip_reg2->leftMeanClipPos;
									end_pos2 = mate_clip_reg2->leftMeanClipPos2;
								}else{
									start_pos2 = mate_clip_reg2->leftMeanClipPos2;
									end_pos2 = mate_clip_reg2->leftMeanClipPos;
								}
							}
						}
					}else if(mate_end_flag==2){
						if(mate_clip_reg2->rightClipRegNum==1){
							if(mate_clip_reg2->rightClipReg){
								chrname2 = mate_clip_reg2->rightClipReg->chrname;
								start_pos2 = mate_clip_reg2->rightClipReg->startRefPos;
								end_pos2 = mate_clip_reg2->rightClipReg->endRefPos;
							}else{
								chrname2 = mate_clip_reg2->rightClipReg2->chrname;
								start_pos2 = mate_clip_reg2->rightClipReg2->startRefPos;
								end_pos2 = mate_clip_reg2->rightClipReg2->endRefPos;
							}
						}else if(mate_clip_reg2->rightClipRegNum==2){
							if(mate_clip_reg2->rightClipReg->chrname.compare(mate_clip_reg2->rightClipReg2->chrname)==0){
								chrname2 = mate_clip_reg2->rightClipReg->chrname;
								if(mate_clip_reg2->rightMeanClipPos<mate_clip_reg2->rightMeanClipPos2){
									start_pos2 = mate_clip_reg2->rightMeanClipPos;
									end_pos2 = mate_clip_reg2->rightMeanClipPos2;
								}else{
									start_pos2 = mate_clip_reg2->rightMeanClipPos2;
									end_pos2 = mate_clip_reg2->rightMeanClipPos;
								}
							}
						}
					}
					if(chrname1.size()>0 and chrname1.compare(chrname2)==0) overlap_flag2 = isOverlappedPos(start_pos1, end_pos1, start_pos2, end_pos2); // same chrome
					else overlap_flag2 = false;
				}

				if(overlap_flag1 and overlap_flag2) overlap_flag = true;
			}else if(mate_clip_reg1->leftClipRegNum==2 and mate_clip_reg1->rightClipRegNum==1){
				// check right region
				chrname1 = chrname2 = "";
				start_pos1 = end_pos1 = start_pos2 = end_pos2 = 0;

				if(mate_clip_reg1->rightClipReg){
					chrname1 = mate_clip_reg1->rightClipReg->chrname;
					start_pos1 = mate_clip_reg1->rightClipReg->startRefPos;
					end_pos1 = mate_clip_reg1->rightClipReg->endRefPos;
				}else{
					chrname1 = mate_clip_reg1->rightClipReg2->chrname;
					start_pos1 = mate_clip_reg1->rightClipReg2->startRefPos;
					end_pos1 = mate_clip_reg1->rightClipReg2->endRefPos;
				}

				// get the mate end
				mate_end_flag = 0;
				if(mate_clip_reg2->leftClipRegNum==1 and mate_clip_reg2->rightClipRegNum==2){
					mate_end_flag = 1;
				}else if(mate_clip_reg2->leftClipRegNum==2 and mate_clip_reg2->rightClipRegNum==1){
					mate_end_flag = 2;
				}

				if(mate_end_flag==1){
					if(mate_clip_reg2->leftClipReg){
						chrname2 = mate_clip_reg2->leftClipReg->chrname;
						start_pos2 = mate_clip_reg2->leftClipReg->startRefPos;
						end_pos2 = mate_clip_reg2->leftClipReg->endRefPos;
					}else{
						chrname2 = mate_clip_reg2->leftClipReg2->chrname;
						start_pos2 = mate_clip_reg2->leftClipReg2->startRefPos;
						end_pos2 = mate_clip_reg2->leftClipReg2->endRefPos;
					}
				}else if(mate_end_flag==2){
					if(mate_clip_reg2->rightClipReg){
						chrname2 = mate_clip_reg2->rightClipReg->chrname;
						start_pos2 = mate_clip_reg2->rightClipReg->startRefPos;
						end_pos2 = mate_clip_reg2->rightClipReg->endRefPos;
					}else{
						chrname2 = mate_clip_reg2->rightClipReg2->chrname;
						start_pos2 = mate_clip_reg2->rightClipReg2->startRefPos;
						end_pos2 = mate_clip_reg2->rightClipReg2->endRefPos;
					}
				}

				if(chrname1.size()>0 and chrname1.compare(chrname2)==0) overlap_flag1 = isOverlappedPos(start_pos1, end_pos1, start_pos2, end_pos2); // same chrome
				else overlap_flag1 = false;

				// check left regions
				chrname1 = chrname2 = "";
				start_pos1 = end_pos1 = start_pos2 = end_pos2 = 0;

				if(mate_end_flag==1) mate_end_flag = 2;
				else if(mate_end_flag==2) mate_end_flag = 1;

				if(mate_clip_reg1->leftClipReg->chrname.compare(mate_clip_reg1->leftClipReg2->chrname)==0){
					chrname1 = mate_clip_reg1->leftClipReg->chrname;
					if(mate_clip_reg1->leftMeanClipPos<mate_clip_reg1->leftMeanClipPos2){
						start_pos1 = mate_clip_reg1->leftMeanClipPos;
						end_pos1 = mate_clip_reg1->leftMeanClipPos2;
					}else{
						start_pos1 = mate_clip_reg1->leftMeanClipPos2;
						end_pos1 = mate_clip_reg1->leftMeanClipPos;
					}

					if(mate_end_flag==1){
						if(mate_clip_reg2->leftClipRegNum==1){
							if(mate_clip_reg2->leftClipReg){
								chrname2 = mate_clip_reg2->leftClipReg->chrname;
								start_pos2 = mate_clip_reg2->leftClipReg->startRefPos;
								end_pos2 = mate_clip_reg2->leftClipReg->endRefPos;
							}else{
								chrname2 = mate_clip_reg2->leftClipReg2->chrname;
								start_pos2 = mate_clip_reg2->leftClipReg2->startRefPos;
								end_pos2 = mate_clip_reg2->leftClipReg2->endRefPos;
							}
						}else if(mate_clip_reg2->leftClipRegNum==2){
							if(mate_clip_reg2->leftClipReg->chrname.compare(mate_clip_reg2->leftClipReg2->chrname)==0){
								chrname2 = mate_clip_reg2->leftClipReg->chrname;
								if(mate_clip_reg2->leftMeanClipPos<mate_clip_reg2->leftMeanClipPos2){
									start_pos2 = mate_clip_reg2->leftMeanClipPos;
									end_pos2 = mate_clip_reg2->leftMeanClipPos2;
								}else{
									start_pos2 = mate_clip_reg2->leftMeanClipPos2;
									end_pos2 = mate_clip_reg2->leftMeanClipPos;
								}
							}
						}
					}else if(mate_end_flag==2){
						if(mate_clip_reg2->rightClipRegNum==1){
							if(mate_clip_reg2->rightClipReg){
								chrname2 = mate_clip_reg2->rightClipReg->chrname;
								start_pos2 = mate_clip_reg2->rightClipReg->startRefPos;
								end_pos2 = mate_clip_reg2->rightClipReg->endRefPos;
							}else{
								chrname2 = mate_clip_reg2->rightClipReg2->chrname;
								start_pos2 = mate_clip_reg2->rightClipReg2->startRefPos;
								end_pos2 = mate_clip_reg2->rightClipReg2->endRefPos;
							}
						}else if(mate_clip_reg2->rightClipRegNum==2){
							if(mate_clip_reg2->rightClipReg->chrname.compare(mate_clip_reg2->rightClipReg2->chrname)==0){
								chrname2 = mate_clip_reg2->rightClipReg->chrname;
								if(mate_clip_reg2->rightMeanClipPos<mate_clip_reg2->rightMeanClipPos2){
									start_pos2 = mate_clip_reg2->rightMeanClipPos;
									end_pos2 = mate_clip_reg2->rightMeanClipPos2;
								}else{
									start_pos2 = mate_clip_reg2->rightMeanClipPos2;
									end_pos2 = mate_clip_reg2->rightMeanClipPos;
								}
							}
						}
					}
					if(chrname1.size()>0 and chrname1.compare(chrname2)==0) overlap_flag2 = isOverlappedPos(start_pos1, end_pos1, start_pos2, end_pos2); // same chrome
					else overlap_flag2 = false;
				}

				if(overlap_flag1 and overlap_flag2) overlap_flag = true;
			}
		}
	}else{ // both large indel regions
		overlap_flag = isOverlappedReg(mate_clip_reg1->largeIndelClipReg, mate_clip_reg2->largeIndelClipReg);
	}

	return overlap_flag;
}

// get the overlapped mate clip region
mateClipReg_t* getOverlappedMateClipReg(mateClipReg_t *mate_clip_reg_given, vector<mateClipReg_t*> &mateClipRegVec){
	mateClipReg_t *mate_clip_reg_overlapped = NULL, *clip_reg;
	bool valid_search_flag;

	for(size_t i=0; i<mateClipRegVec.size(); i++){
		clip_reg = mateClipRegVec.at(i);
		valid_search_flag = false;
		if(clip_reg!=mate_clip_reg_given and clip_reg->valid_flag and clip_reg->sv_type==mate_clip_reg_given->sv_type){
			if(clip_reg->large_indel_flag==false){ // clipping regions
				if(clip_reg->reg_mated_flag)
					valid_search_flag = true;;
			}else{ // large indel
				valid_search_flag = clip_reg;
			}
		}

		if(valid_search_flag){
			if(isOverlappedMateClipReg(mate_clip_reg_given, clip_reg)==true){ // same chrome
				mate_clip_reg_overlapped = clip_reg;
				break;
			}
		}
	}

	return mate_clip_reg_overlapped;
}

// determine whether the indel region in a clipping region
bool isIndelInClipReg(reg_t *reg, vector<mateClipReg_t*> &mate_clipReg_vec){
	bool flag = false;
	mateClipReg_t *clip_reg;

	for(size_t i=0; i<mate_clipReg_vec.size(); i++){
		clip_reg = mate_clipReg_vec.at(i);
		flag = isIndelInSingleClipReg(reg, clip_reg);
		if(flag) break;
	}

	return flag;
}

// determine whether the indel region in a clipping region
bool isIndelInSingleClipReg(reg_t *reg, mateClipReg_t *clip_reg){
	bool flag = false;
	reg_t *reg_tmp;
	string chrname_tmp, chrname_tmp2;
	size_t start_pos, end_pos;

	if(clip_reg->valid_flag and clip_reg->reg_mated_flag){
		if(clip_reg->sv_type==VAR_DUP or clip_reg->sv_type==VAR_INV){ // DUP or INV
			chrname_tmp = "";
			start_pos = end_pos = 0;
			if(clip_reg->leftClipReg){
				chrname_tmp = clip_reg->leftClipReg->chrname;
				start_pos = clip_reg->leftClipReg->startRefPos;
			}else if(clip_reg->leftClipReg2){
				chrname_tmp = clip_reg->leftClipReg2->chrname;
				start_pos = clip_reg->leftClipReg2->startRefPos;
			}
			if(clip_reg->rightClipReg){
				chrname_tmp2 = clip_reg->rightClipReg->chrname;
				end_pos = clip_reg->rightClipReg->endRefPos;
			}else if(clip_reg->rightClipReg2){
				chrname_tmp2 = clip_reg->rightClipReg2->chrname;
				end_pos = clip_reg->rightClipReg2->endRefPos;
			}

			if(reg->chrname.compare(chrname_tmp)==0 and chrname_tmp.compare(chrname_tmp2)==0 and start_pos>0 and end_pos>0){
				if(isOverlappedPos(reg->startRefPos, reg->endRefPos, start_pos, end_pos)){
					flag = true;
				}
			}
		}else if(clip_reg->sv_type==VAR_TRA){  // TRA
			// check left part
			chrname_tmp = "";
			start_pos = end_pos = 0;
			if(clip_reg->leftClipRegNum==2){
				chrname_tmp = clip_reg->leftClipReg->chrname;
				start_pos = clip_reg->leftClipReg->startRefPos;
				end_pos = clip_reg->leftClipReg2->endRefPos;
			}else if(clip_reg->leftClipRegNum==1){
				reg_tmp = clip_reg->leftClipReg ? clip_reg->leftClipReg : clip_reg->leftClipReg2;
				chrname_tmp = reg_tmp->chrname;
				start_pos = reg_tmp->startRefPos;
				end_pos = reg_tmp->endRefPos;
			}
			if(reg->chrname.compare(chrname_tmp)==0 and start_pos>0 and end_pos>0){
				if(isOverlappedPos(reg->startRefPos, reg->endRefPos, start_pos, end_pos)){
					flag = true;
				}
			}

			// check right part
			if(flag==false){
				chrname_tmp = "";
				start_pos = end_pos = 0;
				if(clip_reg->rightClipRegNum==2){
					chrname_tmp = clip_reg->rightClipReg->chrname;
					start_pos = clip_reg->rightClipReg->startRefPos;
					end_pos = clip_reg->rightClipReg2->endRefPos;
				}else if(clip_reg->rightClipRegNum==1){
					reg_tmp = clip_reg->rightClipReg ? clip_reg->rightClipReg : clip_reg->rightClipReg2;
					chrname_tmp = reg_tmp->chrname;
					start_pos = reg_tmp->startRefPos;
					end_pos = reg_tmp->endRefPos;
				}
				if(reg->chrname.compare(chrname_tmp)==0 and start_pos>0 and end_pos>0){
					if(isOverlappedPos(reg->startRefPos, reg->endRefPos, start_pos, end_pos)){
						flag = true;
					}
				}
			}
		}
	}

	return flag;
}

// determine whether the SNV position in a clipping region
bool isSnvInClipReg(string &chrname, size_t pos, vector<mateClipReg_t*> &mate_clipReg_vec){
	bool flag = false;
	mateClipReg_t *clip_reg;

	for(size_t i=0; i<mate_clipReg_vec.size(); i++){
		clip_reg = mate_clipReg_vec.at(i);
		flag = isSnvInSingleClipReg(chrname, pos, clip_reg);
		if(flag) break;
	}

	return flag;
}

// determine whether the SNV position in a clipping region
bool isSnvInSingleClipReg(string &chrname, size_t pos, mateClipReg_t *clip_reg){
	bool flag = false;
	reg_t *reg_tmp;
	string chrname_tmp, chrname_tmp2;
	size_t start_pos, end_pos;

	if(clip_reg->valid_flag and clip_reg->reg_mated_flag){
		if(clip_reg->sv_type==VAR_DUP or clip_reg->sv_type==VAR_INV){ // DUP or INV
			chrname_tmp = "";
			start_pos = end_pos = 0;
			if(clip_reg->leftClipReg){
				chrname_tmp = clip_reg->leftClipReg->chrname;
				start_pos = clip_reg->leftClipReg->startRefPos;
			}else if(clip_reg->leftClipReg2){
				chrname_tmp = clip_reg->leftClipReg2->chrname;
				start_pos = clip_reg->leftClipReg2->startRefPos;
			}
			if(clip_reg->rightClipReg){
				chrname_tmp2 = clip_reg->rightClipReg->chrname;
				end_pos = clip_reg->rightClipReg->endRefPos;
			}else if(clip_reg->rightClipReg2){
				chrname_tmp2 = clip_reg->rightClipReg2->chrname;
				end_pos = clip_reg->rightClipReg2->endRefPos;
			}

			if(chrname.compare(chrname_tmp)==0 and chrname_tmp.compare(chrname_tmp2)==0 and start_pos>0 and end_pos>0){
				if(pos>=start_pos and pos<=end_pos){
					flag = true;
				}
			}
		}else if(clip_reg->sv_type==VAR_TRA){  // TRA
			// check left part
			chrname_tmp = "";
			start_pos = end_pos = 0;
			if(clip_reg->leftClipRegNum==2){
				chrname_tmp = clip_reg->leftClipReg->chrname;
				start_pos = clip_reg->leftClipReg->startRefPos;
				end_pos = clip_reg->leftClipReg2->endRefPos;
			}else if(clip_reg->leftClipRegNum==1){
				reg_tmp = clip_reg->leftClipReg ? clip_reg->leftClipReg : clip_reg->leftClipReg2;
				chrname_tmp = reg_tmp->chrname;
				start_pos = reg_tmp->startRefPos;
				end_pos = reg_tmp->endRefPos;
			}

			if(chrname.compare(chrname_tmp)==0 and start_pos>0 and end_pos>0){
				if(pos>=start_pos and pos<=end_pos){
					flag = true;
				}
			}

			// check right part
			if(flag==false){
				chrname_tmp = "";
				start_pos = end_pos = 0;
				if(clip_reg->rightClipRegNum==2){
					chrname_tmp = clip_reg->rightClipReg->chrname;
					start_pos = clip_reg->rightClipReg->startRefPos;
					end_pos = clip_reg->rightClipReg2->endRefPos;
				}else if(clip_reg->rightClipRegNum==1){
					reg_tmp = clip_reg->rightClipReg ? clip_reg->rightClipReg : clip_reg->rightClipReg2;
					chrname_tmp = reg_tmp->chrname;
					start_pos = reg_tmp->startRefPos;
					end_pos = reg_tmp->endRefPos;
				}

				if(chrname.compare(chrname_tmp)==0 and start_pos>0 and end_pos>0){
					if(pos>=start_pos and pos<=end_pos){
						flag = true;
					}
				}
			}
		}
	}

	return flag;
}

// load sam/bam header
bam_hdr_t* loadSamHeader(string &inBamFile)
{
	bam_hdr_t *header = NULL;
	samFile *in = 0;

    if ((in = sam_open(inBamFile.c_str(), "r")) == 0) {
        cerr << __func__ << ": failed to open " << inBamFile << " for reading" << endl;
        exit(1);
    }

    if ((header = sam_hdr_read(in)) == 0) {
        cerr << "fail to read the header from " << inBamFile << endl;
        exit(1);
    }
    sam_close(in);

	return header;
}

// determine whether a position in an region or not.
bool isInReg(int32_t pos, vector<reg_t*> &vec){
	bool flag = false;
	reg_t* reg;
	for(size_t i=0; i<vec.size(); i++){
		reg = vec.at(i);
		if(pos>=reg->startRefPos and pos<=reg->endRefPos) { flag = true; break; }
	}
	return flag;
}

// compute the number of disagreements
int32_t computeDisagreeNum(Base *baseArray, int32_t arr_size){
	int32_t i, disagreeNum = 0;
	for(i=0; i<arr_size; i++)
		if(baseArray[i].isZeroCovBase() or baseArray[i].isDisagreeBase())
			disagreeNum ++;
	return disagreeNum;
}

// merge overlapped indels
void mergeOverlappedReg(vector<reg_t*> &regVector){
	size_t i, j;
	reg_t *reg1, *reg2;

	for(i=0; i<regVector.size(); i++){
		reg1 = regVector.at(i);
		if(reg1->call_success_status){
			for(j=i+1; j<regVector.size(); ){
				reg2 = regVector.at(j);
				if(reg2->call_success_status and isOverlappedReg(reg1, reg2) and reg1->query_id!=-1 and reg1->query_id==reg2->query_id){ // same reference region and query
					updateReg(reg1, reg2);
					delete reg2;
					regVector.erase(regVector.begin()+j);
				}else j++;
			}
		}
	}
	regVector.shrink_to_fit();
}

// merge two overlapped indels
void updateReg(reg_t* reg1, reg_t* reg2){
	size_t new_startRefPos, new_endRefPos, new_startLocalRefPos, new_endLocalRefPos, new_startQueryPos, new_endQueryPos;
	if(reg1->startRefPos<=reg2->startRefPos) {
		new_startRefPos = reg1->startRefPos;
		new_startLocalRefPos = reg1->startLocalRefPos;
		new_startQueryPos = reg1->startQueryPos;
	}else{
		new_startRefPos = reg2->startRefPos;
		new_startLocalRefPos = reg2->startLocalRefPos;
		new_startQueryPos = reg2->startQueryPos;
	}
	if(reg1->endRefPos>=reg2->endRefPos){
		new_endRefPos = reg1->endRefPos;
		new_endLocalRefPos = reg1->endLocalRefPos;
		new_endQueryPos = reg1->endQueryPos;
	}else{
		new_endRefPos = reg2->endRefPos;
		new_endLocalRefPos = reg2->endLocalRefPos;
		new_endQueryPos = reg2->endQueryPos;
	}
	reg1->startRefPos = new_startRefPos;
	reg1->endRefPos = new_endRefPos;
	reg1->startLocalRefPos = new_startLocalRefPos;
	reg1->endLocalRefPos = new_endLocalRefPos;
	reg1->startQueryPos = new_startQueryPos;
	reg1->endQueryPos = new_endQueryPos;
}

// merge adjacent regions
void mergeAdjacentReg(vector<reg_t*> &regVec, size_t dist_thres){
	size_t i;
	reg_t *reg_tmp, *reg_tmp2;
	bool flag;
	for(i=1; i<regVec.size(); ){
		reg_tmp = regVec.at(i-1);
		reg_tmp2 = regVec.at(i);
		flag = isAdjacent(reg_tmp->startRefPos, reg_tmp->endRefPos, reg_tmp2->startRefPos, reg_tmp2->endRefPos, dist_thres);
		if(flag){ // merge
			if(reg_tmp->startRefPos>reg_tmp2->startRefPos) reg_tmp->startRefPos = reg_tmp2->startRefPos;
			if(reg_tmp->endRefPos<reg_tmp2->endRefPos) reg_tmp->endRefPos = reg_tmp2->endRefPos;
			delete reg_tmp2;
			regVec.erase(regVec.begin()+i);
		}else i ++;
	}
}

void printRegVec(vector<reg_t*> &regVec, string header){
	reg_t *reg;
	cout << header << " region size: " << regVec.size() << endl;
	for(size_t i=0; i<regVec.size(); i++){
		reg = regVec.at(i);
		cout << "[" << i << "]: " << reg->chrname << ":" << reg->startRefPos << "-" << reg->endRefPos << ", localRef: " << reg->startLocalRefPos << "-" << reg->endLocalRefPos << ", Query :" << reg->startQueryPos << "-" << reg->endQueryPos << endl;
	}
}

void printMateClipReg(mateClipReg_t *mate_clip_reg){
	string chrname_left, chrname_right;
	size_t start_pos_left, end_pos_left, start_pos_right, end_pos_right;

	// left part
	if(mate_clip_reg->leftClipRegNum==2){
		chrname_left = mate_clip_reg->leftClipReg->chrname;
		start_pos_left = mate_clip_reg->leftMeanClipPos;
		end_pos_left = mate_clip_reg->leftMeanClipPos2;
	}else if(mate_clip_reg->leftClipRegNum==1){
		if(mate_clip_reg->leftClipReg){
			chrname_left = mate_clip_reg->leftClipReg->chrname;
			start_pos_left = mate_clip_reg->leftClipReg->startRefPos;
			end_pos_left = mate_clip_reg->leftClipReg->endRefPos;
		}else if(mate_clip_reg->leftClipReg2){
			chrname_left = mate_clip_reg->leftClipReg2->chrname;
			start_pos_left = mate_clip_reg->leftClipReg2->startRefPos;
			end_pos_left = mate_clip_reg->leftClipReg2->endRefPos;
		}
	}else
		cout << "No region in left part" << endl;

	// right part
	if(mate_clip_reg->rightClipRegNum==2){
		chrname_right = mate_clip_reg->rightClipReg->chrname;
		start_pos_right = mate_clip_reg->rightMeanClipPos;
		end_pos_right = mate_clip_reg->rightMeanClipPos2;
	}else if(mate_clip_reg->rightClipRegNum==1){
		if(mate_clip_reg->rightClipReg){
			chrname_right = mate_clip_reg->rightClipReg->chrname;
			start_pos_right = mate_clip_reg->rightClipReg->startRefPos;
			end_pos_right = mate_clip_reg->rightClipReg->endRefPos;
		}else if(mate_clip_reg->rightClipReg2){
			chrname_right = mate_clip_reg->rightClipReg2->chrname;
			start_pos_right = mate_clip_reg->rightClipReg2->startRefPos;
			end_pos_right = mate_clip_reg->rightClipReg2->endRefPos;
		}
	}else
		cout << "No region in right part" << endl;

	cout << chrname_left << ":" << start_pos_left << "-" << end_pos_left << ", " << chrname_right << ":" << start_pos_right << "-" << end_pos_right << " :" << endl;
	cout << "\tvalid_flag=" << mate_clip_reg->valid_flag << ", call_success_flag=" << mate_clip_reg->call_success_flag << ", leftClipRegNum=" << mate_clip_reg->leftClipRegNum << ", rightClipRegNum=" << mate_clip_reg->rightClipRegNum << ", sv_type=" << mate_clip_reg->sv_type << ", dup_num=" << mate_clip_reg->dup_num << endl;
}

vector<string> getLeftRightPartChrname(mateClipReg_t *mate_clip_reg){
	vector<string> chrname_vec;
	string chrname_left = "", chrname_right = "";

	// left part
	if(mate_clip_reg->leftClipRegNum==2){
		chrname_left = mate_clip_reg->leftClipReg->chrname;
	}else if(mate_clip_reg->leftClipRegNum==1){
		if(mate_clip_reg->leftClipReg)
			chrname_left = mate_clip_reg->leftClipReg->chrname;
		else if(mate_clip_reg->leftClipReg2)
			chrname_left = mate_clip_reg->leftClipReg2->chrname;
	}

	// right part
	if(mate_clip_reg->rightClipRegNum==2){
		chrname_right = mate_clip_reg->rightClipReg->chrname;
	}else if(mate_clip_reg->rightClipRegNum==1){
		if(mate_clip_reg->rightClipReg)
			chrname_right = mate_clip_reg->rightClipReg->chrname;
		else if(mate_clip_reg->rightClipReg2)
			chrname_right = mate_clip_reg->rightClipReg2->chrname;
	}

	chrname_vec.push_back(chrname_left);
	chrname_vec.push_back(chrname_right);

	return chrname_vec;
}

// preprocess pipe command characters
string preprocessPipeChar(string &cmd_str){
	char ch;
	string ret_str = "";

	for(size_t i=0; i<cmd_str.size(); i++){
		ch = cmd_str.at(i);
		if(ch=='|') // pipe character
			ret_str += "_";
		else
			ret_str += ch;
	}
	return ret_str;
}

bool isFileExist(const string &filename){
	bool flag = false;
	struct stat fileStat;
	if (stat(filename.c_str(), &fileStat) == 0)
		if(fileStat.st_size>0)
			flag = true;
	return flag;
}

void removeRedundantItems(vector<reg_t*> &reg_vec){
	size_t i, j;
	reg_t *reg, *reg_tmp;

	for(i=0; i<reg_vec.size(); i++){
		reg = reg_vec.at(i);
		for(j=i+1; j<reg_vec.size(); ){
			reg_tmp = reg_vec.at(j);
			if(reg_tmp->chrname.compare(reg->chrname)==0 and reg_tmp->startRefPos==reg->startRefPos and reg_tmp->endRefPos==reg->endRefPos){ // redundant item
				delete reg_tmp;
				reg_vec.erase(reg_vec.begin()+j);
			}else j++;
		}
	}
}

// get line count of file
int32_t getLineCount(string &filename){
	int32_t num;
	ifstream infile;
	string line;

	infile.open(filename);
	if(!infile.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << filename << endl;
		exit(1);
	}

	num = 0;
	while(getline(infile, line)) if(line.size()>0 and line.at(0)!='#') num ++;
	infile.close();

	return num;
}

bool isBaseMatch(char ctgBase, char refBase){
	bool match_flag = false;

	// Upper case
	if(ctgBase>='a' and ctgBase<='z') ctgBase -= 32;
	if(refBase>='a' and refBase<='z') refBase -= 32;

	if(ctgBase==refBase){
		match_flag = true;
	}else if(ctgBase=='-' or refBase=='-'){
		match_flag = false;
	}else{
		switch(refBase){
			case 'A':
			case 'a':
			case 'C':
			case 'c':
			case 'G':
			case 'g':
			case 'T':
			case 't':
			case 'N':
			case 'n':
				match_flag = false;
				break;
			case 'M':
			case 'm':
				if(ctgBase=='A' or ctgBase=='C') match_flag = true;
				break;
			case 'R':
			case 'r':
				if(ctgBase=='A' or ctgBase=='G') match_flag = true;
				break;
			case 'S':
			case 's':
				if(ctgBase=='C' or ctgBase=='G') match_flag = true;
				break;
			case 'V':
			case 'v':
				if(ctgBase=='A' or ctgBase=='C' or ctgBase=='G') match_flag = true;
				break;
			case 'W':
			case 'w':
				if(ctgBase=='A' or ctgBase=='T') match_flag = true;
				break;
			case 'Y':
			case 'y':
				if(ctgBase=='C' or ctgBase=='T') match_flag = true;
				break;
			case 'H':
			case 'h':
				if(ctgBase=='A' or ctgBase=='C' or ctgBase=='T') match_flag = true;
				break;
			case 'K':
			case 'k':
				if(ctgBase=='G' or ctgBase=='T') match_flag = true;
				break;
			case 'D':
			case 'd':
				if(ctgBase=='A' or ctgBase=='G' or ctgBase=='T') match_flag = true;
				break;
			case 'B':
			case 'b':
				if(ctgBase=='C' or ctgBase=='G' or ctgBase=='T') match_flag = true;
				break;
			default: cerr << __func__ << ": unknown base: " << refBase << endl; exit(1);
		}
	}

	return match_flag;
}

bool isRegValid(reg_t *reg, int32_t min_size){
	bool flag = false;
	if(reg->startRefPos<=reg->endRefPos and reg->startLocalRefPos<=reg->endLocalRefPos and reg->startQueryPos<=reg->endQueryPos){
		if(reg->endRefPos-reg->startRefPos+1>=min_size or reg->endQueryPos-reg->startQueryPos+1>=min_size)
			flag = true;
	}
//	if(flag==false){
//		cout << "========= Invalid region: " << reg->chrname << ":" << reg->startRefPos << "-" << reg->endRefPos << ", local_Loc: " << reg->startLocalRefPos << "-" << reg->endLocalRefPos << ", query_Loc: " << reg->startQueryPos << "-" << reg->endQueryPos << endl;
//	}
	return flag;
}

void exchangeRegLoc(reg_t *reg){
	int64_t tmp;
	if(reg->startRefPos>reg->endRefPos){
		//cout << "startRefPos<-->endRefPos and startLocalRefPos<-->endLocalRefPos has been exchanged." << endl;
		tmp = reg->startRefPos;
		reg->startRefPos = reg->endRefPos;
		reg->endRefPos = tmp;
		tmp = reg->startLocalRefPos;
		reg->startLocalRefPos = reg->endLocalRefPos;
		reg->endLocalRefPos = tmp;
	}

	if(reg->startQueryPos>reg->endQueryPos){
		tmp = reg->startQueryPos;
		reg->startQueryPos = reg->endQueryPos;
		reg->endQueryPos = tmp;
	}
}

//Minimap2 alignment, and the output is in paf format
int32_t minimap2Aln(string &alnfilename, string &contigfilename, string &refseqfilename, int32_t max_proc_running_minutes){
	string minimap2_cmd, out_opt;
	int32_t i, ret_status, status, sleep_sec;
	time_t start_time, end_time; // start time and end time
	double cost_min;

	//minimap2_cmd = "minimap2 -c --secondary=no -o " + alnfilename + " " + refseqfilename + " " + contigfilename + " > /dev/null 2>&1";
	minimap2_cmd = "minimap2 -c -x asm5 -o " + alnfilename + " " + refseqfilename + " " + contigfilename + " > /dev/null 2>&1";

	//cout << "blat_cmd: " + blat_cmd << endl;
	time(&start_time);
	end_time = 0;

	ret_status = 0;
	for(i=0; i<3; i++){
		time(&end_time);
		cost_min = difftime(end_time, start_time) / 60.0;
		if(cost_min<=MAX_ALN_MINUTES){
			status = system(minimap2_cmd.c_str());
			ret_status = getSuccessStatusSystemCmd(status);
			if(ret_status!=0){ // command executed failed
				//cout << __func__ << ", line=" << __LINE__ << ": ret_status=" << ret_status << ", blat_cmd=" << blat_cmd << endl;

				if(i<2){
					time(&end_time);
					cost_min = difftime(end_time, start_time) / 60.0;
					if(cost_min<=MAX_ALN_MINUTES){
						sleep_sec = (i + 1) * (i + 2);
						sleep(sleep_sec);
						cout << __func__ << ": retry aligning " << contigfilename << endl;
					}
				}else{
					cerr << "Please run the correct minimap2 command or check whether minimap2 was correctly installed when aligning " << contigfilename <<"." << endl;
					//exit(1);
				}
			}else break;
		}else break;
	}

//		if(ret_status!=0){
//			time(&end_time);
//			cost_min = difftime(end_time, start_time) / 60.0;
//			if(cost_min>max_proc_running_minutes) { // killed work
//				ret_status = -2;
//			}
//		}

	return ret_status;
}

// BLAT alignment, and the output is in sim4 format
int32_t blatAln(string &alnfilename, string &contigfilename, string &refseqfilename, int32_t max_proc_running_minutes){
	string blat_cmd, out_opt;
	int32_t i, ret_status, status, sleep_sec;
	time_t start_time, end_time; // start time and end time
	double cost_min;

	out_opt = "-out=sim4 " + alnfilename;
	blat_cmd = "blat " + refseqfilename + " " + contigfilename + " " + out_opt + " > /dev/null 2>&1";

	//cout << "blat_cmd: " + blat_cmd << endl;
	time(&start_time);
	end_time = 0;

	ret_status = 0;
	for(i=0; i<3; i++){
		time(&end_time);
		cost_min = difftime(end_time, start_time) / 60.0;
		if(cost_min<=MAX_ALN_MINUTES){
			status = system(blat_cmd.c_str());
			ret_status = getSuccessStatusSystemCmd(status);
			if(ret_status!=0){ // command executed failed
				//cout << __func__ << ", line=" << __LINE__ << ": ret_status=" << ret_status << ", blat_cmd=" << blat_cmd << endl;

				if(i<2){
					time(&end_time);
					cost_min = difftime(end_time, start_time) / 60.0;
					if(cost_min<=MAX_ALN_MINUTES){
						sleep_sec = (i + 1) * (i + 2);
						sleep(sleep_sec);
						cout << __func__ << ": retry aligning " << contigfilename << endl;
					}
				}else{
					cerr << "Please run the correct blat command or check whether blat was correctly installed when aligning " << contigfilename <<"." << endl;
					//exit(1);
				}
			}else break;
		}else break;
	}

	if(ret_status!=0){
		time(&end_time);
		cost_min = difftime(end_time, start_time) / 60.0;
		if(cost_min>max_proc_running_minutes) { // killed work
			ret_status = -2;
		}
	}

	/*
	 https://www.cnblogs.com/alantu2018/p/8554281.html
	 http://blog.chinaunix.net/uid-22263887-id-3023260.html
	 https://www.cnblogs.com/zhaoyl/p/4963551.html
	 https://blog.csdn.net/zybasjj/article/details/8231837
	 */

	return ret_status;
}

// get the success status to determine whether the command is executed successfully
// 		0: successful; -1: failed;
int32_t getSuccessStatusSystemCmd(int32_t status){
	int32_t ret_status = 0;

	if(status==-1) ret_status = -1;
	else{
		if(WIFEXITED(status)){
			if(WEXITSTATUS(status)!=0)
				ret_status = -1;
		}else ret_status = -1;
	}

	return ret_status;
}

// check whether the blat alignment is matched to the contig
bool isMinimap2AlnResultMatch(string &contigfilename, string &alnfilename){
	bool flag;
	vector<string> query_name_vec_ctg, query_name_vec_minimap2, line_vec;
	vector<size_t> query_len_vec_minimap2;
	ifstream infile;
	string line, query_name_minimap2, query_name_ctg, query_name;
	size_t i, query_len_minimap2, query_len_ctg ,query_len;
	int32_t query_loc;

	infile.open(alnfilename);
	if(!infile.is_open()){
		cerr << __func__ << "(), line=" << __LINE__ << ": cannot open file:" << alnfilename << endl;
		exit(1);
	}

	// load query names
	FastaSeqLoader fa_loader(contigfilename);
	query_name_vec_ctg = fa_loader.getFastaSeqNames();


	// get blat alignment information
	while(getline(infile, line)){
		if(line.size()){
			line_vec = split(line, "	");
			query_name = line_vec.at(0);
			query_len = stoi(line_vec.at(1));

			// get the query name and length
			query_name_vec_minimap2.push_back(query_name);
			query_len_vec_minimap2.push_back(query_len);

		}
	}
	infile.close();

	// check the match items
	flag = true;
	for(i=0; i<query_name_vec_minimap2.size(); i++){
		query_name_minimap2 = query_name_vec_minimap2.at(i);
		query_len_minimap2 = query_len_vec_minimap2.at(i);

		// search the query from contig vector
		query_loc = getQueryNameLoc(query_name_minimap2, query_name_vec_ctg);
		if(query_loc!=-1){ // query name matched
			query_len_ctg = fa_loader.getFastaSeqLen(query_loc);
			if(query_len_minimap2!=query_len_ctg) {
				flag = false;
				break;
			}
		}else{ // query name unmatched
			flag = false;
			break;
		}
	}

	if(flag==false){
		cerr << __func__ << ", line=" << __LINE__ << ", it is not matched between the minimap2 alignment and its corresponding contig, contigfile=" << contigfilename << ", alnfile=" << alnfilename << ", error!" << endl;
	}
	return flag;
}

// check whether the blat alignment is matched to the contig
bool isBlatAlnResultMatch(string &contigfilename, string &alnfilename){
	bool flag;
	vector<string> query_name_vec_ctg, query_name_vec_blat, line_vec, line_vec1, len_vec;
	vector<size_t> query_len_vec_blat;
	ifstream infile;
	string line, query_name_blat, query_name_ctg;
	size_t i, query_len_blat, query_len_ctg;
	int32_t query_loc;

	infile.open(alnfilename);
	if(!infile.is_open()){
		cerr << __func__ << "(), line=" << __LINE__ << ": cannot open file:" << alnfilename << endl;
		exit(1);
	}

	// load query names
	FastaSeqLoader fa_loader(contigfilename);
	query_name_vec_ctg = fa_loader.getFastaSeqNames();


	// get blat alignment information
	while(getline(infile, line)){
		if(line.size()){
			if(line.substr(0, 3).compare("seq")==0){  // query and subject titles
				line_vec = split(line, "=");
				line_vec1 = split(line_vec[1].substr(1), ",");

				if(line_vec1.size()==2){
					if(line_vec[0][3]=='1') {
						// get the query name and length
						query_name_vec_blat.push_back(line_vec1.at(0));
						len_vec = split(line_vec1[line_vec1.size()-1].substr(1), " ");
						query_len_vec_blat.push_back(stoi(len_vec[0]));
					}
				}
			}
		}
	}
	infile.close();

	// check the match items
	flag = true;
	for(i=0; i<query_name_vec_blat.size(); i++){
		query_name_blat = query_name_vec_blat.at(i);
		query_len_blat = query_len_vec_blat.at(i);

		// search the query from contig vector
		query_loc = getQueryNameLoc(query_name_blat, query_name_vec_ctg);
		if(query_loc!=-1){ // query name matched
			query_len_ctg = fa_loader.getFastaSeqLen(query_loc);
			if(query_len_blat!=query_len_ctg) {
				flag = false;
				break;
			}
		}else{ // query name unmatched
			flag = false;
			break;
		}
	}

	if(flag==false){
		cerr << __func__ << ", line=" << __LINE__ << ", it is not matched between the blat alignment and its corresponding contig, contigfile=" << contigfilename << ", alnfile=" << alnfilename << ", error!" << endl;
	}
	return flag;
}

// get query name location
int32_t getQueryNameLoc(string &query_name, vector<string> &query_name_vec){
	int32_t i, loc = -1, query_name_len;
	string query_name_tmp;
	query_name_len = query_name.size();
	if(query_name_vec.size()){
		for(i=query_name_vec.size()-1; i>=0; i--){
			query_name_tmp = query_name_vec.at(i).substr(0, query_name_len);
			if(query_name_tmp.compare(query_name)==0) { loc = i; break; }
		}
	}
	return loc;
}

// determine whether the blat align is completed
bool isBlatAlnCompleted(string &alnfilename){
	bool flag = true;
	string line;
	ifstream infile;
	vector<string> line_vec, line_vec1, len_vec;

	if(!isFileExist(alnfilename)) return false;

	infile.open(alnfilename);
	if(!infile.is_open()){
		cerr << __func__ << "(), line=" << __LINE__ << ": cannot open file:" << alnfilename << endl;
		exit(1);
	}

	// get data
	while(getline(infile, line)){
		if(line.size()){
			if(line.substr(0, 3).compare("seq")==0){  // query and subject titles
				line_vec = split(line, "=");
				line_vec1 = split(line_vec[1].substr(1), ",");

				if(line_vec1.size()<2){
//					cout << line_vec1.size() << endl;
					flag = false;
					break;
				}

				len_vec = split(line_vec1[line_vec1.size()-1].substr(1), " ");
				if(line_vec[0][3]=='1' or line_vec[0][3]=='2') {
					if(len_vec.at(1).compare("bp")!=0){
						flag = false;
						break;
					}
				}else{
					flag = false;
					break;
				}
			}else if(line[0]=='('){ // complement
				if(line.size()==12){
					if(line.substr(1, 10).compare("complement")!=0){
						flag = false;
						break;
					}
				}else{
					flag = false;
					break;
				}
			}else{    // segments
				line_vec = split(line, " ");
				if(line_vec.size()==4){
					if(line_vec.at(3).compare("->")!=0 and line_vec.at(3).compare("-")!=0 and line_vec.at(3).compare("<-")!=0 and line_vec.at(3).compare("<")!=0){
						flag = false;
						break;
					}
				}else{
					flag = false;
					break;
				}
			}
		}
	}
	infile.close();

	return flag;
}

// clean consensus temporary folders
void cleanPrevConsTmpDir(const string &cns_dir_str, const string &dir_prefix)
{
	DIR *dp;
	struct dirent *entry;
	struct stat statbuf;
	string cmd_str, dir_str;
	char *path;

	// get current work directory
	path = getcwd(NULL, 0);

	if(isFileExist(cns_dir_str)){
		if((dp=opendir(cns_dir_str.c_str()))==NULL){
			cerr << "cannot open directory: " << cns_dir_str << endl;
			exit(1);
		}
		chdir(cns_dir_str.c_str());
		while((entry = readdir(dp)) != NULL){
			lstat(entry->d_name, &statbuf);
			if(S_ISDIR(statbuf.st_mode)){
				if(strcmp(".", entry->d_name) == 0 || strcmp("..", entry->d_name) == 0) continue;

				if(strlen(entry->d_name)>4){
					dir_str = entry->d_name;
					if(dir_str.substr(0, 4).compare(dir_prefix)==0){
						cmd_str = "rm -rf ";
						cmd_str += entry->d_name;
						system(cmd_str.c_str());
					}
				}
			}
		}
		closedir(dp);
		chdir(path);
	}
	free(path);
}

// get call file header line for INDEL which starts with '#'
string getCallFileHeaderBed(string &sample){
	string header_line;
	header_line = "#Chr\tStart\tEnd\tSVType\tSVLen\tDupNum\tRef\tAlt\tExtra\tFORMAT\t" + sample;
	return header_line;
}

// get call file header line for INDEL which starts with '#'
string getCallFileHeaderBedpe(string &sample){
	string header_line;
	header_line = "#Chr1\tStart1\tEnd1\tChr2\tStart2\tEnd2\tSVType\tSVLen1\tSVLen2\tMateReg\tRef1\tAlt1\tRef2\tAlt2\tFORMAT\t" + sample;
	return header_line;
}

// extract queries from clip align data vector
vector<struct querySeqInfoNode*> extractQueriesFromClipAlnDataVec(vector<clipAlnData_t*> &clipAlnDataVector, const string &refseq_given, string &chrname, int64_t startRefPos, int64_t endRefPos, faidx_t *fai, int32_t minConReadLen, bool clip_reg_flag, pthread_mutex_t *p_mutex_fai){
	vector<struct querySeqInfoNode*> query_seq_info_all;
	size_t i, j;
	string qname, seq, reg_str, refseq;
	int64_t noHardClipIdx;
	int32_t seq_len, bam_type;
	char *p_seq;
	vector<clipAlnData_t*> query_aln_segs;
	vector<string> query_seq_qual_vec;
	struct querySeqInfoNode *q_node;
	vector<struct alnSeg*> alnSegs;

	if(clipAlnDataVector.size()>0){
		for(i=0; i<clipAlnDataVector.size(); i++) clipAlnDataVector.at(i)->query_checked_flag = false;

		reg_str = chrname + ":" + to_string(startRefPos) + "-" + to_string(endRefPos);
		pthread_mutex_lock(p_mutex_fai);
		p_seq = fai_fetch(fai, reg_str.c_str(), &seq_len);
		pthread_mutex_unlock(p_mutex_fai);
		refseq = p_seq;
		free(p_seq);

		// check reference identical
		if(refseq_given.size()>0){
			if(refseq.compare(refseq_given)!=0){
				cerr << __func__ << ": different reference sequences, error!" << endl;
				exit(1);
			}
		}

		bam_type = getBamType(clipAlnDataVector);
		if(bam_type==BAM_INVALID){
			cerr << __func__ << ": unknown bam type, error!" << endl;
			exit(1);
		}

		// join query clip align segments
		for(i=0; i<clipAlnDataVector.size(); i++){
			if(clipAlnDataVector.at(i)->query_checked_flag==false){
				qname = clipAlnDataVector.at(i)->queryname;
				query_aln_segs = getQueryClipAlnSegs(qname, clipAlnDataVector);  // get query clip align segments

//				if(qname.compare("SRR8858443.1.96142")==0){
//					cout << "i=" << i << ", qname=" << qname << endl;
//				}else continue;

			//	no_otherchrname_flag = true;
				noHardClipIdx = getNoHardClipAlnItem(query_aln_segs);
				if(noHardClipIdx!=-1){
//					for(size_t t=0; t<query_aln_segs.size(); t++){
//						if(t!=noHardClipIdx) {
//							if(query_aln_segs.at(t)->chrname.compare(varVec[0]->chrname)!=0){
//								no_otherchrname_flag = false;
//							}
//						}
//					}
					//if(no_otherchrname_flag){
						//for(j=0; j<query_aln_segs.size(); j++){
					query_seq_qual_vec = getQuerySeqWithSoftClipSeqs(query_aln_segs.at(noHardClipIdx), refseq, startRefPos, endRefPos, bam_type, clip_reg_flag);
					//markHardClipSegs(noHardClipIdx, query_aln_segs);
					query_aln_segs.at(noHardClipIdx)->query_checked_flag = true;

					if(query_seq_qual_vec.size()>0 and query_seq_qual_vec.at(0).length()>(size_t)minConReadLen){
						seq = query_seq_qual_vec.at(0);
//						if(clip_reg_flag){
//							qual = query_seq_qual_vec.at(1);
//							fq_node = new struct fqSeqNode();
//							fq_node->seq_id = seq_id++;
//							fq_node->seq_name = qname;
//							fq_node->seq = seq;
//							fq_node->qual = qual;
//							fq_node->selected_flag = true;
//							fq_seq_vec.push_back(fq_node);
//
//							total_bases_original += seq.size();
//						}else{
							//seqs.push_back(seq);
							//prepare for smoothing
							seq = query_seq_qual_vec.at(0);
							q_node = new struct querySeqInfoNode();
							q_node->clip_aln = query_aln_segs.at(noHardClipIdx);
							q_node->qname = qname;
							q_node->seq = seq;
							q_node->selected_flag = true;
							//q_node->seq_id = fq_node->seq_id;
							switch(bam_type){
								case BAM_CIGAR_NO_DIFF_MD:
									//alnSegs = generateAlnSegs(b);
									alnSegs = generateAlnSegs2(query_aln_segs.at(noHardClipIdx)->bam, startRefPos, endRefPos);
									break;
								case BAM_CIGAR_NO_DIFF_NO_MD:
								case BAM_CIGAR_DIFF_MD:
								case BAM_CIGAR_DIFF_NO_MD:
									//alnSegs = generateAlnSegs_no_MD(b, refseq, startPos, endPos);
									alnSegs = generateAlnSegs_no_MD2(query_aln_segs.at(noHardClipIdx)->bam, refseq, startRefPos, endRefPos);
									break;
								default:
									cerr << __func__ << ": unknown bam type, error!" << endl;
									exit(1);
							}// generate align segments
							if(clip_reg_flag==false)
								removeHeadTailAlnSegs(alnSegs);
							q_node->query_alnSegs = alnSegs;

							query_seq_info_all.push_back(q_node);
//						}
					}
						//}
				//	}
				}else{ // all hard clipping
					//cout << "qname=" << qname << ", querylen=" << clipAlnDataVector.at(i)->querylen << endl;
					// join query align segments without 'SA' tags
//					query_seq_qual_vec = joinQueryAlnSegs(query_aln_segs);
					for(j=0; j<query_aln_segs.size(); j++){
						query_seq_qual_vec = getQuerySeqWithSoftClipSeqs(query_aln_segs.at(j), refseq, startRefPos, endRefPos, bam_type, clip_reg_flag);
						if(query_seq_qual_vec.size()>0 and query_seq_qual_vec.at(0).length()>(size_t)minConReadLen){
							seq = query_seq_qual_vec.at(0);
//							if(clip_reg_flag){
//								qual = query_seq_qual_vec.at(1);
//								fq_node = new struct fqSeqNode();
//								fq_node->seq_id = seq_id++;
//								fq_node->seq_name = qname;
//								fq_node->seq = seq;
//								fq_node->qual = qual;
//								fq_node->selected_flag = true;
//								fq_seq_vec.push_back(fq_node);
//
//								total_bases_original += seq.size();
//							}else{
								//seqs.push_back(seq);
								//prepare for smoothing
								q_node = new struct querySeqInfoNode();
								q_node->clip_aln = query_aln_segs.at(j);
								q_node->qname = qname;
								q_node->seq = seq;
								q_node->selected_flag = true;
								//q_node->seq_id = fq_node->seq_id;

								// generate align segments
								switch(bam_type){
									case BAM_CIGAR_NO_DIFF_MD:
										//alnSegs = generateAlnSegs(b);
										alnSegs = generateAlnSegs2(query_aln_segs.at(j)->bam, startRefPos, endRefPos);
										break;
									case BAM_CIGAR_NO_DIFF_NO_MD:
									case BAM_CIGAR_DIFF_MD:
									case BAM_CIGAR_DIFF_NO_MD:
										//alnSegs = generateAlnSegs_no_MD(b, refseq, startPos, endPos);
										alnSegs = generateAlnSegs_no_MD2(query_aln_segs.at(j)->bam, refseq, startRefPos, endRefPos);
										break;
									default:
										cerr << __func__ << ": unknown bam type, error!" << endl;
										exit(1);
								}
								if(clip_reg_flag==false)
									removeHeadTailAlnSegs(alnSegs);
								q_node->query_alnSegs = alnSegs;
								query_seq_info_all.push_back(q_node);
//							}
							//destroyAlnSegs(query_alnSegs);
						}
					}
					for(j=0; j<query_aln_segs.size(); j++) query_aln_segs.at(j)->query_checked_flag = true;
				}
			}
		}
		for(i=0; i<clipAlnDataVector.size(); i++) clipAlnDataVector.at(i)->query_checked_flag = false;
	}

	return query_seq_info_all;
}

// get the non hard-clip align item
int32_t getNoHardClipAlnItem(vector<clipAlnData_t*> &query_aln_segs){
	int32_t idx;
	clipAlnData_t *clip_aln;

	idx = -1;
	for(size_t i=0; i<query_aln_segs.size(); i++){
		clip_aln = query_aln_segs.at(i);
		if(clip_aln->leftHardClippedFlag==false and clip_aln->rightHardClippedFlag==false){
			idx = i;
			break;
		}
	}

	return idx;
}

// get query with soft clipped sequence
vector<string> getQuerySeqWithSoftClipSeqs(clipAlnData_t* clip_aln, string &refseq, int64_t startRefPos, int64_t endRefPos, int32_t bam_type, bool clip_reg_flag){
	vector<string> query_info_vec;
	uint8_t *seq_int, *qual_int;
	string qseq, qual;
	int32_t i, query_start_loc, query_end_loc; // 0-based
	vector<int32_t> query_loc_vec;

	seq_int = bam_get_seq(clip_aln->bam);
	qual_int = bam_get_qual(clip_aln->bam);

	qseq = qual = "";
	if(clip_reg_flag==false){
		//if(seq_int[0]!=255){
			//queryStartEndLoc_vec = computeQueryStartEndLocByRefPos(clipAlnData_t *clip_aln, startRefPos, endRefPos);
		query_loc_vec = computeQueryStartEndLocByRefPos(clip_aln->bam, startRefPos, endRefPos, refseq, bam_type);
		if(query_loc_vec.size()>0){
			query_start_loc = query_loc_vec.at(0);
			query_end_loc = query_loc_vec.at(1);

	//			for(i=0; i<clip_aln->bam->core.l_qseq; i++) qseq += "=ACMGRSVTWYHKDBN"[bam_seqi(seq_int, i)];  // seq
	//
	//			if(qual_int[0]!=255) for(i=0; i<clip_aln->bam->core.l_qseq; i++) qual += qual_int[i] + 33;  // qual
	//			else for(i=0; i<clip_aln->bam->core.l_qseq; i++) qual += 34;  // qual

			for(i=query_start_loc - 1; i<query_end_loc; i++) qseq += "=ACMGRSVTWYHKDBN"[bam_seqi(seq_int, i)];  // seq

			if(qual_int[0]!=255) for(i=query_start_loc - 1; i<query_end_loc; i++) qual += qual_int[i] + 33;  // qual
			else for(i=query_start_loc - 1; i<query_end_loc; i++) qual += 34;  // qual

			query_info_vec.push_back(qseq);
			query_info_vec.push_back(qual);
		}
	}else{
//		for(i=0; i<sizeof(seq_int)/sizeof(uint8_t); i++) qseq += "=ACMGRSVTWYHKDBN"[bam_seqi(seq_int, i)];
//		if(qual_int[0]!=255) for(i=0; i<sizeof(qual_int)/sizeof(uint8_t); i++) qual += qual_int[i] + 33;  // qual
//		else for(i=0; i<sizeof(qual_int)/sizeof(uint8_t); i++) qual += 34;  // qual

		for(i=0; i<clip_aln->bam->core.l_qseq; i++) qseq += "=ACMGRSVTWYHKDBN"[bam_seqi(seq_int, i)];  // seq

		if(qual_int[0]!=255) for(i=0; i<clip_aln->bam->core.l_qseq; i++) qual += qual_int[i] + 33;  // qual
		else for(i=0; i<clip_aln->bam->core.l_qseq; i++) qual += 34;  // qual

		query_info_vec.push_back(qseq);
		query_info_vec.push_back(qual);
	}

	//}

	return query_info_vec;
}

// compute query locations according to given reference positions
vector<int32_t> computeQueryStartEndLocByRefPos(bam1_t *b, int64_t startRefPos, int64_t endRefPos, string &refseq, int32_t bam_type){
	vector<int32_t> query_loc_vec;
	uint32_t op;
	int32_t target_startQpos = 0, target_endQpos = 0;
	//vector<struct alnSeg*> alnSegs;
	vector<struct alnSeg*> query_alnSegs;

	if(!(b->core.flag & BAM_FUNMAP)){ // aligned

//		string qname = bam_get_qname(b);
//		cout << "qname=" << qname << endl;

		switch(bam_type){
			case BAM_CIGAR_NO_DIFF_MD:
				//alnSegs = generateAlnSegs(b);
				query_alnSegs = generateAlnSegs2(b, startRefPos, endRefPos);
				break;
			case BAM_CIGAR_NO_DIFF_NO_MD:
			case BAM_CIGAR_DIFF_MD:
			case BAM_CIGAR_DIFF_NO_MD:
				//alnSegs = generateAlnSegs_no_MD(b, refseq, startPos, endPos);
				query_alnSegs = generateAlnSegs_no_MD2(b, refseq, startRefPos, endRefPos);
				break;
			default:
				cerr << __func__ << ": unknown bam type, error!" << endl;
				exit(1);
		}// generate align segments
		//outputAlnSegs(alnSegs);

		op = query_alnSegs.at(0)->opflag;
		if(query_alnSegs.size()>0){
			switch (op) {
				case BAM_CMATCH:
					target_startQpos = query_alnSegs.at(0)->startQpos;
					break;
				case BAM_CINS:
					target_startQpos = query_alnSegs.at(0)->startQpos + query_alnSegs.at(0)->seglen;
					delete *query_alnSegs.begin();
					query_alnSegs.erase(query_alnSegs.begin());
					break;
				case BAM_CDEL:
					target_startQpos = query_alnSegs.at(0)->startQpos;
					delete *query_alnSegs.begin();
					query_alnSegs.erase(query_alnSegs.begin());
					break;
				case BAM_CSOFT_CLIP:
				case BAM_CHARD_CLIP:
					target_startQpos = query_alnSegs.at(0)->startQpos + query_alnSegs.at(0)->seglen;
					delete *query_alnSegs.begin();
					query_alnSegs.erase(query_alnSegs.begin());
					break;
				case BAM_CEQUAL:
				case BAM_CDIFF:
					target_startQpos = query_alnSegs.at(0)->startQpos;
					break;
				case BAM_CREF_SKIP:
					// unexpected events
				case BAM_CPAD:
				case BAM_CBACK:
				default:
					cerr << __func__ << ", line=" << __LINE__ << ": invalid opflag " << op << endl;
					exit(1);
			}
		}

		if(query_alnSegs.size()>0){
			op = query_alnSegs.at(query_alnSegs.size() - 1)->opflag;
			switch(op){
				case BAM_CMATCH:
					target_endQpos = query_alnSegs.at(query_alnSegs.size() - 1)->startQpos +  query_alnSegs.at(query_alnSegs.size() - 1)->seglen - 1;
					break;
				case BAM_CINS:
					target_endQpos = query_alnSegs.at(query_alnSegs.size() - 1)->startQpos - 1;
					delete *(query_alnSegs.end() - 1);
					query_alnSegs.erase(query_alnSegs.end() - 1);
					break;
				case BAM_CDEL:
					target_endQpos = query_alnSegs.at(query_alnSegs.size() - 1)->startQpos - 1;
					delete *(query_alnSegs.end() - 1);
					query_alnSegs.erase(query_alnSegs.end() - 1);
					break;
				case BAM_CSOFT_CLIP:
				case BAM_CHARD_CLIP:
					target_endQpos = query_alnSegs.at(query_alnSegs.size() - 1)->startQpos - 1;
					delete *(query_alnSegs.end() - 1);
					query_alnSegs.erase(query_alnSegs.end() - 1);
					break;
				case BAM_CEQUAL:
				case BAM_CDIFF:
					target_endQpos = query_alnSegs.at(query_alnSegs.size() - 1)->startQpos +  query_alnSegs.at(query_alnSegs.size() - 1)->seglen - 1;
					break;
				case BAM_CREF_SKIP:
					// unexpected events
				case BAM_CPAD:
				case BAM_CBACK:
				default:
					cerr << __func__ << ", line=" << __LINE__ << ": invalid opflag " << op << endl;
					exit(1);
			}

			if(query_alnSegs.size()>0){
				destroyAlnSegs(query_alnSegs);
			}
			// save query locations
			query_loc_vec.push_back(target_startQpos);
			query_loc_vec.push_back(target_endQpos);
		}
	}

	return query_loc_vec;
}

// remove head align segments
void removeHeadTailAlnSegs(vector<struct alnSeg*> &alnSegs){
	uint32_t op;
	op = alnSegs.at(0)->opflag;
	switch (op) {
		case BAM_CMATCH:
			break;
		case BAM_CINS:
			delete *alnSegs.begin();
			alnSegs.erase(alnSegs.begin());
			break;
		case BAM_CDEL:
			delete *alnSegs.begin();
			alnSegs.erase(alnSegs.begin());
			break;
		case BAM_CSOFT_CLIP:
		case BAM_CHARD_CLIP:
			delete *alnSegs.begin();
			alnSegs.erase(alnSegs.begin());
			break;
		case BAM_CEQUAL:
		case BAM_CDIFF:
			break;
		case BAM_CREF_SKIP:
			// unexpected events
		case BAM_CPAD:
		case BAM_CBACK:
		default:
			cerr << __func__ << ", line=" << __LINE__ << ": invalid opflag " << op << endl;
			exit(1);
	}
	if(alnSegs.size()>0){
		op = alnSegs.at(alnSegs.size() - 1)->opflag;
		switch(op){
			case BAM_CMATCH:
				break;
			case BAM_CINS:
				delete *(alnSegs.end() - 1);
				alnSegs.erase(alnSegs.end() - 1);
				break;
			case BAM_CDEL:
				delete *(alnSegs.end() - 1);
				alnSegs.erase(alnSegs.end() - 1);
				break;
			case BAM_CSOFT_CLIP:
			case BAM_CHARD_CLIP:
				delete *(alnSegs.end() - 1);
				alnSegs.erase(alnSegs.end() - 1);
				break;
			case BAM_CEQUAL:
			case BAM_CDIFF:
				break;
			case BAM_CREF_SKIP:
				// unexpected events
			case BAM_CPAD:
			case BAM_CBACK:
			default:
				cerr << __func__ << ", line=" << __LINE__ << ": invalid opflag " << op << endl;
				exit(1);
		}
	}

}

cnsWork_opt* allocateCnsWorkOpt(string &chrname, string &readsfilename, string &contigfilename, string &refseqfilename, string &tmpdir, vector<reg_t*> &varVec, bool clip_reg_flag, bool limit_reg_process_flag, vector<simpleReg_t*> &limit_reg_vec){
	cnsWork_opt *cns_work_opt;
	cns_work_opt = new cnsWork_opt();
	cns_work_opt->chrname = chrname;
	cns_work_opt->readsfilename = readsfilename;
	cns_work_opt->contigfilename = contigfilename;
	cns_work_opt->refseqfilename = refseqfilename;
	cns_work_opt->tmpdir = tmpdir;
	cns_work_opt->clip_reg_flag = clip_reg_flag;

	// sub-regions
	cns_work_opt->var_array = (reg_t**)malloc(varVec.size()*sizeof(reg_t*));
	cns_work_opt->arr_size = varVec.size();
	for(size_t i=0; i<varVec.size(); i++) cns_work_opt->var_array[i] = varVec.at(i);

	// limit process regions
	cns_work_opt->limit_reg_process_flag = limit_reg_process_flag;
	cns_work_opt->limit_reg_array = NULL;
	cns_work_opt->limit_reg_array_size = 0;
	if(limit_reg_process_flag and limit_reg_vec.size()){
		cns_work_opt->limit_reg_array = (simpleReg_t**)malloc(limit_reg_vec.size()*sizeof(simpleReg_t*));
		cns_work_opt->limit_reg_array_size = limit_reg_vec.size();
		for(size_t i=0; i<limit_reg_vec.size(); i++) cns_work_opt->limit_reg_array[i] = limit_reg_vec.at(i);
	}

	return cns_work_opt;
}

void releaseConsWorkOpt(cnsWork_opt *cns_work_opt){
	free(cns_work_opt->var_array); cns_work_opt->var_array = NULL;
	if(cns_work_opt->limit_reg_array_size>0) { free(cns_work_opt->limit_reg_array); cns_work_opt->limit_reg_array = NULL; cns_work_opt->limit_reg_array_size = 0; }
	delete cns_work_opt;
}

void destroyConsWorkOptVec(vector<cnsWork_opt*> &cns_work_vec){
	cnsWork_opt *cns_work_opt;
	for(size_t i=0; i<cns_work_vec.size(); i++){
		cns_work_opt = cns_work_vec.at(i);
		releaseConsWorkOpt(cns_work_opt);
	}
	vector<cnsWork_opt*>().swap(cns_work_vec);
}

void deleteItemFromCnsWorkVec(int32_t item_id, vector<cnsWork_opt*> &cns_work_vec){
	if(item_id!=-1) {
		releaseConsWorkOpt(cns_work_vec.at(item_id));
		cns_work_vec.erase(cns_work_vec.begin()+item_id);
	}
}

// get consensus work id from vector
int32_t getItemIDFromCnsWorkVec(string &contigfilename, vector<cnsWork_opt*> &cns_work_vec){
	cnsWork_opt *cns_work_item;
	int32_t vec_idx = -1;

	for(size_t i=0; i<cns_work_vec.size(); i++){
		cns_work_item = cns_work_vec.at(i);
		if(cns_work_item->contigfilename.compare(contigfilename)==0){
			vec_idx = i;
			break;
		}
	}

	return vec_idx;
}

// process single consensus work
void* processSingleConsWork(void *arg){
	cnsWork *cns_work = (cnsWork *)arg;
	cnsWork_opt *cns_work_opt = cns_work->cns_work_opt;
	vector<reg_t*> varVec;
	vector<simpleReg_t*> sub_limit_reg_vec;
	size_t i, num_done, num_work, num_work_percent;
	double percentage;
	Time time;

//	pthread_mutex_lock(cns_work->p_mtx_cns_reg_workDone_num);
//	cout << "consensus region [" << cns_work->work_id << "]: " << cns_work_opt->readsfilename << endl;
//	pthread_mutex_unlock(cns_work->p_mtx_cns_reg_workDone_num);

	for(i=0; i<cns_work_opt->arr_size; i++) varVec.push_back(cns_work_opt->var_array[i]);
	for(i=0; i<cns_work_opt->limit_reg_array_size; i++) sub_limit_reg_vec.push_back(cns_work_opt->limit_reg_array[i]);

	performLocalCons(cns_work_opt->readsfilename, cns_work_opt->contigfilename, cns_work_opt->refseqfilename, cns_work_opt->tmpdir, cns_work->technology, cns_work->canu_version, cns_work->num_threads_per_cns_work, varVec, cns_work_opt->chrname, cns_work->inBamFile, cns_work->fai, cns_work->cnsSideExtSize, *(cns_work->var_cand_file), cns_work->expected_cov_cns, cns_work->min_input_cov_canu, cns_work->delete_reads_flag, cns_work->keep_failed_reads_flag, cns_work_opt->clip_reg_flag, cns_work->minClipEndSize, cns_work->minConReadLen, cns_work->min_sv_size, cns_work->min_supp_num, cns_work->max_seg_size_ratio, cns_work_opt->limit_reg_process_flag, sub_limit_reg_vec);

	// release memory
	sub_limit_reg_vec.clear();
	if(cns_work_opt->arr_size>0) { free(cns_work_opt->var_array); cns_work_opt->var_array = NULL; cns_work_opt->arr_size = 0; }
	if(cns_work_opt->limit_reg_array_size>0) { free(cns_work_opt->limit_reg_array); cns_work_opt->limit_reg_array = NULL; cns_work_opt->limit_reg_array_size = 0; }

	// output progress information
	pthread_mutex_lock(cns_work->p_mtx_cns_reg_workDone_num);
	(*cns_work->p_cns_reg_workDone_num) ++;
	num_done = *cns_work->p_cns_reg_workDone_num;

	num_work = cns_work->num_work;
	num_work_percent = cns_work->num_work_percent;
	if(num_done==1 and num_done!=num_work){
		num_done = 0;
		percentage = (double)num_done / num_work * 100;
		cout << "[" << time.getTime() << "]: processed regions: " << num_done << "/" << num_work << " (" << percentage << "%)" << endl;
	}else if(num_done%num_work_percent==0 or num_done==num_work){
		percentage = (double)num_done / num_work * 100;
		cout << "[" << time.getTime() << "]: processed regions: " << num_done << "/" << num_work << " (" << percentage << "%)" << endl;
	}
	pthread_mutex_unlock(cns_work->p_mtx_cns_reg_workDone_num);

	delete (cnsWork *)arg;

	return NULL;
}


void performLocalCons(string &readsfilename, string &contigfilename, string &refseqfilename, string &tmpdir, string &technology, string &canu_version, size_t num_threads_per_cns_work, vector<reg_t*> &varVec, string &chrname, string &inBamFile, faidx_t *fai, int32_t cns_extend_size, ofstream &cns_info_file, double expected_cov_cns, double min_input_cov_canu, bool delete_reads_flag, bool keep_failed_reads_flag, bool clip_reg_flag, int32_t minClipEndSize, int32_t minConReadLen, int32_t min_sv_size, int32_t min_supp_num, double max_seg_size_ratio, bool limit_reg_process_flag, vector<simpleReg_t*> &limit_reg_vec){

	LocalCns local_cns(readsfilename, contigfilename, refseqfilename, tmpdir, technology, canu_version, num_threads_per_cns_work, varVec, chrname, inBamFile, fai, cns_extend_size, expected_cov_cns, min_input_cov_canu, delete_reads_flag, keep_failed_reads_flag, clip_reg_flag, minClipEndSize, minConReadLen, min_sv_size, min_supp_num, max_seg_size_ratio);

	local_cns.setLimitRegs(limit_reg_process_flag, limit_reg_vec);
	if(local_cns.cns_success_preDone_flag==false){
		// extract the corresponding refseq from reference
		local_cns.extractRefseq();

		// extract the reads data from BAM file
		local_cns.extractReadsDataFromBAM();

		if(local_cns.clip_reg_flag==false){
			local_cns.cnsByPoa(); // local consensus using abPOA
		}else
			local_cns.localCnsWtdbg2(); // local consensus using wtdbg2
	}

	// record consensus information
	local_cns.recordCnsInfo(cns_info_file);

	// empty the varVec
	varVec.clear();
}

// determine whether a file is readable
bool isReadableFile(string &filename){
	bool flag = true;;
	ifstream infile;
	string line;

	if(!isFileExist(filename)) flag = false;

	infile.open(filename);
	if(!infile.is_open()) flag = false;

	// get data
	getline(infile, line);
	if(line.size()==0) flag = false;

	return flag;
}

// process single minimap2 work
void* processSingleMinimap2AlnWork(void *arg){
	callWork_opt *call_work_opt = (callWork_opt *)arg;
	varCand *var_cand = call_work_opt->var_cand;
	size_t num_done, num_work, num_work_percent;
	double percentage;
	Time time;

	//cout << var_cand->alnfilename << endl;

	var_cand->alnCtg2Refseq02();

	// output progress information
	pthread_mutex_lock(call_work_opt->p_mtx_call_workDone_num);
	(*call_work_opt->p_call_workDone_num) ++;
	num_done = *call_work_opt->p_call_workDone_num;

	num_work = call_work_opt->num_work;
	num_work_percent = call_work_opt->num_work_percent;
	if(num_done==1 and num_done!=num_work){
		num_done = 0;
		percentage = (double)num_done / num_work * 100;
		cout << "[" << time.getTime() << "]: processed regions: " << num_done << "/" << num_work << " (" << percentage << "%)" << endl;
	}else if(num_done%num_work_percent==0 or num_done==num_work){
		percentage = (double)num_done / num_work * 100;
		cout << "[" << time.getTime() << "]: processed regions: " << num_done << "/" << num_work << " (" << percentage << "%)" << endl;
	}
	//cout << "\t[" << call_work_opt->work_id << "]: " << var_cand->alnfilename << ", sv_type=" << var_cand->sv_type << endl;
	pthread_mutex_unlock(call_work_opt->p_mtx_call_workDone_num);

	delete (callWork_opt *)arg;

	return NULL;
}

// process single blat align work
void* processSingleBlatAlnWork(void *arg){
	callWork_opt *call_work_opt = (callWork_opt *)arg;
	varCand *var_cand = call_work_opt->var_cand;
	size_t num_done, num_work, num_work_percent;
	double percentage;
	Time time;

	//cout << var_cand->alnfilename << endl;

	var_cand->alnCtg2Refseq();

	// output progress information
	pthread_mutex_lock(call_work_opt->p_mtx_call_workDone_num);
	(*call_work_opt->p_call_workDone_num) ++;
	num_done = *call_work_opt->p_call_workDone_num;

	num_work = call_work_opt->num_work;
	num_work_percent = call_work_opt->num_work_percent;
	if(num_done==1 and num_done!=num_work){
		num_done = 0;
		percentage = (double)num_done / num_work * 100;
		cout << "[" << time.getTime() << "]: processed regions: " << num_done << "/" << num_work << " (" << percentage << "%)" << endl;
	}else if(num_done%num_work_percent==0 or num_done==num_work){
		percentage = (double)num_done / num_work * 100;
		cout << "[" << time.getTime() << "]: processed regions: " << num_done << "/" << num_work << " (" << percentage << "%)" << endl;
	}
	//cout << "\t[" << call_work_opt->work_id << "]: " << var_cand->alnfilename << ", sv_type=" << var_cand->sv_type << endl;
	pthread_mutex_unlock(call_work_opt->p_mtx_call_workDone_num);

	delete (callWork_opt *)arg;

	return NULL;
}

// process single call work
void* processSingleCallWork(void *arg){
	callWork_opt *call_work_opt = (callWork_opt *)arg;
	varCand *var_cand = call_work_opt->var_cand;
	size_t num_done, num_work, num_work_percent;
	double percentage;
	Time time;

//	cout << var_cand->alnfilename << endl;

	//var_cand->callVariants();
	var_cand->callVariants02();

	// output progress information
	pthread_mutex_lock(call_work_opt->p_mtx_call_workDone_num);
	(*call_work_opt->p_call_workDone_num) ++;
	num_done = *call_work_opt->p_call_workDone_num;

	num_work = call_work_opt->num_work;
	num_work_percent = call_work_opt->num_work_percent;
	if(num_done==1 and num_done!=num_work){
		num_done = 0;
		percentage = (double)num_done / num_work * 100;
		cout << "[" << time.getTime() << "]: processed regions: " << num_done << "/" << num_work << " (" << percentage << "%)" << endl;
	}else if(num_done%num_work_percent==0 or num_done==num_work){
		percentage = (double)num_done / num_work * 100;
		cout << "[" << time.getTime() << "]: processed regions: " << num_done << "/" << num_work << " (" << percentage << "%)" << endl;
	}
	//cout << "\t[" << call_work_opt->work_id << "]: " << var_cand->alnfilename << ", sv_type=" << var_cand->sv_type << endl;
	pthread_mutex_unlock(call_work_opt->p_mtx_call_workDone_num);

	delete (callWork_opt *)arg;

	return NULL;
}

void outputCnsWorkOptToFile_debug(vector<cnsWork_opt*> &cns_work_opt_vec){
	cnsWork_opt *cns_work_opt;
//	reg_t *reg;
	ofstream outfile;

	outfile.open("cns_work_list");
	if(!outfile.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file" << endl;
		exit(1);
	}

	for(size_t i=0; i<cns_work_opt_vec.size(); i++){
		cns_work_opt = cns_work_opt_vec.at(i);
		outfile << cns_work_opt->refseqfilename << "\t" << cns_work_opt->contigfilename << "\t" << cns_work_opt->readsfilename << endl;
	}

	outfile.close();
}

// get old output directory name
string getOldOutDirname(string &filename, string &sub_work_dir){
	string old_dir = "";
	size_t pos = filename.find(sub_work_dir);
	if(pos!=filename.npos){
		if(pos==0) old_dir = "";
		else{
			old_dir = filename.substr(0, pos);
			old_dir = deleteTailPathChar(old_dir);
		}
	}else{
		cerr << __func__ << ", line=" << __LINE__ << ", filename=" << filename << ", cannot find old output directory for " << sub_work_dir << ", error!" << endl;
		exit(1);
	}
	return old_dir;
}

// get updated item file name
string getUpdatedItemFilename(string &filename, string &out_dir, string &old_out_dir){
	string new_filename;

	if(old_out_dir.size()==0)
		new_filename = out_dir + "/" + filename;
	else if(filename.at(old_out_dir.size())=='/')
		new_filename = out_dir + filename.substr(old_out_dir.size());
	else
		new_filename = out_dir + "/" + filename.substr(old_out_dir.size());

	return new_filename;
}

// get chromosome name by given file name
string getChrnameByFilename(string &filename){
	string chrname;
	int32_t start_idx, end_idx;

	if(filename.size()>0){
		start_idx = end_idx = -1;
		end_idx = filename.find_last_of("/") - 1;
		for(int32_t i=end_idx; i>=0; i--){
			if(filename.at(i)=='/'){
				start_idx = i + 1;
				break;
			}
		}
		if(start_idx!=-1 and end_idx!=-1)
			chrname = filename.substr(start_idx, end_idx-start_idx+1);
	}

	return chrname;
}

// delete the tail '/' path character
string deleteTailPathChar(string &dirname){
	if(dirname.size()>0 and dirname.at(dirname.size()-1)=='/'){
		return dirname.substr(0, dirname.size()-1);
	}else{
		return dirname;
	}
}

// get overlapped simple regions
vector<simpleReg_t*> getOverlappedSimpleRegsExt(string &chrname, int64_t begPos, int64_t endPos, vector<simpleReg_t*> &limit_reg_vec, int32_t ext_size){
	simpleReg_t *simple_reg;
	vector<simpleReg_t*> sub_simple_reg_vec;
	bool flag;
	int64_t new_startRegPos, new_endRegPos;

	if(chrname.size()){
		for(size_t i=0; i<limit_reg_vec.size(); i++){
			simple_reg = limit_reg_vec.at(i);
			flag = false;
			if(simple_reg->chrname.compare(chrname)==0){
				if((begPos==-1 and endPos==-1) or (simple_reg->startPos==-1 and simple_reg->endPos==-1)) flag = true;
				else{
					new_startRegPos = simple_reg->startPos - ext_size;
					if(new_startRegPos<1) new_startRegPos = 1;
					new_endRegPos = simple_reg->endPos + ext_size;
					if(isOverlappedPos(begPos, endPos, new_startRegPos, new_endRegPos)) flag = true;
				}

				if(flag) sub_simple_reg_vec.push_back(simple_reg);
			}
		}
	}

	return sub_simple_reg_vec;
}

// get overlapped simple regions
vector<simpleReg_t*> getOverlappedSimpleRegs(string &chrname, int64_t begPos, int64_t endPos, vector<simpleReg_t*> &limit_reg_vec){
	simpleReg_t *simple_reg;
	vector<simpleReg_t*> sub_simple_reg_vec;
	bool flag;

	if(chrname.size()){
		for(size_t i=0; i<limit_reg_vec.size(); i++){
			simple_reg = limit_reg_vec.at(i);
			flag = false;
			if(simple_reg->chrname.compare(chrname)==0){
				if((begPos==-1 and endPos==-1) or (simple_reg->startPos==-1 and simple_reg->endPos==-1)) flag = true;
				else if(isOverlappedPos(begPos, endPos, simple_reg->startPos, simple_reg->endPos)) flag = true;

				if(flag) sub_simple_reg_vec.push_back(simple_reg);
			}
		}
	}

	return sub_simple_reg_vec;
}

// get fully contained simple regions
vector<simpleReg_t*> getFullyContainedSimpleRegs(string &chrname, int64_t begPos, int64_t endPos, vector<simpleReg_t*> &limit_reg_vec){
	simpleReg_t *simple_reg;
	vector<simpleReg_t*> sub_simple_reg_vec;
	bool flag;

	if(chrname.size()){
		for(size_t i=0; i<limit_reg_vec.size(); i++){
			simple_reg = limit_reg_vec.at(i);
			flag = false;
			if(simple_reg->chrname.compare(chrname)==0){
				if((begPos==-1 and endPos==-1) and (simple_reg->startPos==-1 and simple_reg->endPos==-1)){ // chr, chr --> true
					flag = true;
				}else if(simple_reg->startPos==-1 and simple_reg->endPos==-1){ // chr:start-end, chr --> true
					flag = true;
				}else{ // CHR:START-END fully contained in 'simple_reg'
					if(isFullyContainedReg(chrname, begPos, endPos, simple_reg->chrname, simple_reg->startPos, simple_reg->endPos))
						flag = true;
				}
			}
			if(flag) sub_simple_reg_vec.push_back(simple_reg);
		}
	}
	return sub_simple_reg_vec;
}

// determine whether region1 is fully contained in region2
bool isFullyContainedReg(string &chrname1, int64_t begPos1,  int64_t endPos1, string &chrname2, int64_t startPos2, int64_t endPos2){
	bool flag = false;
	if(chrname1.compare(chrname2)==0 and begPos1>=startPos2 and endPos1<=endPos2) flag = true;
	return flag;
}

// get simple regions whose region contain the given start or end position
vector<simpleReg_t*> getPosContainedSimpleRegs(string &chrname, int64_t begPos, int64_t endPos, vector<simpleReg_t*> &limit_reg_vec){
	simpleReg_t *simple_reg;
	vector<simpleReg_t*> sub_simple_reg_vec;
	bool flag;

	if(chrname.size()){
		for(size_t i=0; i<limit_reg_vec.size(); i++){
			simple_reg = limit_reg_vec.at(i);
			flag = false;
			if(simple_reg->chrname.compare(chrname)==0){
				if(simple_reg->startPos==-1 and simple_reg->endPos==-1){ // pos, chr --> true
					flag = true;
				}else if((begPos>=simple_reg->startPos and begPos<=simple_reg->endPos) or (endPos>=simple_reg->startPos and endPos<=simple_reg->endPos)){
					flag = true;
				}
			}
			if(flag) sub_simple_reg_vec.push_back(simple_reg);
		}
	}

	return sub_simple_reg_vec;
}

// determine whether the start position is contained in region2
bool isPosContained(string &chrname1, int64_t begPos1, string &chrname2, int64_t startPos2, int64_t endPos2){
	bool flag = false;
	if(chrname1.compare(chrname2)==0 and begPos1>=startPos2 and begPos1<=endPos2) flag = true;
	return flag;
}

// allocate simple region node
simpleReg_t* allocateSimpleReg(string &simple_reg_str){
	simpleReg_t *simple_reg = NULL;
	vector<string> str_vec, pos_vec;
	string chrname_tmp;
	int64_t pos1, pos2;

	str_vec = split(simple_reg_str, ":");
	if(str_vec.size()){
		chrname_tmp = str_vec.at(0);
		if(str_vec.size()==1){
			pos1 = pos2 = -1;
		}else if(str_vec.size()==2){
			pos_vec = split(str_vec.at(1), "-");
			if(pos_vec.size()==2){
				pos1 = stoi(pos_vec.at(0));
				pos2 = stoi(pos_vec.at(1));
				if(pos1>pos2) goto fail;
			}else goto fail;
		}else goto fail;
	}else goto fail;

	simple_reg = new simpleReg_t();
	simple_reg->chrname = chrname_tmp;
	simple_reg->startPos = pos1;
	simple_reg->endPos = pos2;

	return simple_reg;

fail:
	cout << "Error: Invalid region: " << simple_reg_str << endl;
	return NULL;
}

// free the memory of limited regions
void destroyLimitRegVector(vector<simpleReg_t*> &limit_reg_vec){
	vector<simpleReg_t*>::iterator s_reg;
	for(s_reg=limit_reg_vec.begin(); s_reg!=limit_reg_vec.end(); s_reg++)
		delete (*s_reg);   // free each node
	vector<simpleReg_t*>().swap(limit_reg_vec);
}

// print limit regions
void printLimitRegs(vector<simpleReg_t*> &limit_reg_vec, string &description){
	simpleReg_t *limit_reg;
	string limit_reg_str;

	if(limit_reg_vec.size()){
		cout << endl << description << endl;
		for(size_t i=0; i<limit_reg_vec.size(); i++){
			limit_reg = limit_reg_vec.at(i);
			limit_reg_str = limit_reg->chrname;
			if(limit_reg->startPos!=-1 and limit_reg->endPos!=-1) limit_reg_str += ":" + to_string(limit_reg->startPos) + "-" + to_string(limit_reg->endPos);
			cout << "region [" << i << "]: " << limit_reg_str << endl;
		}
		cout << endl;
	}
}

string getLimitRegStr(vector<simpleReg_t*> &limit_reg_vec){
	simpleReg_t *limit_reg;
	string limit_reg_str;

	if(limit_reg_vec.size()){
		limit_reg = limit_reg_vec.at(0);
		limit_reg_str = limit_reg->chrname;
		if(limit_reg->startPos!=-1 and limit_reg->endPos!=-1) limit_reg_str += ":" + to_string(limit_reg->startPos) + "-" + to_string(limit_reg->endPos);
		for(size_t i=1; i<limit_reg_vec.size(); i++){
			limit_reg = limit_reg_vec.at(i);
			limit_reg_str = limit_reg->chrname;
			if(limit_reg->startPos!=-1 and limit_reg->endPos!=-1) limit_reg_str += ":" + to_string(limit_reg->startPos) + "-" + to_string(limit_reg->endPos);
			limit_reg_str += ":" + limit_reg_str;
		}
	}

	return limit_reg_str;
}

// get region by given file name
void getRegByFilename(simpleReg_t *reg, string &filename, string &pattern_str){
	string simple_filename, chr_pos_str, chr_str, pos_str;
	vector<string> pos_vec;
	size_t pos, pattern_size;

	reg->chrname = "";
	reg->startPos = reg->endPos = -1;
	pattern_size = pattern_str.size();

	pos = filename.find_last_of("/");			// path/pattern_chr_start-end
	simple_filename = filename.substr(pos+1);	// pattern_chr_start-end
	pos = simple_filename.find(pattern_str);
	chr_pos_str = simple_filename.substr(pos+pattern_size+1); // chr_start-end
	pos = chr_pos_str.find_last_of("_");
	chr_str = chr_pos_str.substr(0, pos);
	pos_str = chr_pos_str.substr(pos+1);
	pos_vec = split(pos_str, "-");

	reg->chrname = chr_str;
	reg->startPos = stoi(pos_vec.at(0));
	reg->endPos = stoi(pos_vec.at(1));
}

// extract simple regions by string formated regions
vector<simpleReg_t*> extractSimpleRegsByStr(string &regs_str){
	vector<simpleReg_t*> limit_reg_vec;
	simpleReg_t *simple_reg;
	vector<string> reg_vec, chr_pos_vec, pos_vec;
	string reg_str, chrname_reg;
	int64_t start_pos, end_pos;

	reg_vec = split(regs_str, ";");
	for(size_t i=0; i<reg_vec.size(); i++){
		reg_str = reg_vec.at(i);	// CHR or CHR:START-END

		chr_pos_vec = split(reg_str, ":");
		chrname_reg = chr_pos_vec.at(0);
		if(chr_pos_vec.size()==1)
			start_pos = end_pos = -1;
		else if(chr_pos_vec.size()==2){
			pos_vec = split(chr_pos_vec.at(1), "-");
			start_pos = stoi(pos_vec.at(0));
			end_pos = stoi(pos_vec.at(1));
		}else{
			cerr << __func__ << ", line=" << __LINE__ << ", invalid simple region: " << reg_str << ", error!" << endl;
			exit(1);
		}

		simple_reg = new simpleReg_t();
		simple_reg->chrname = chrname_reg;
		simple_reg->startPos = start_pos;
		simple_reg->endPos = end_pos;

		limit_reg_vec.push_back(simple_reg);
	}

	return limit_reg_vec;
}

// create directory
void createDir(string &dirname){
	string dirname_tmp;
	vector<string> str_vec;
	int32_t ret;

	str_vec = split(dirname, "/");
	dirname_tmp = str_vec.at(0);
	for(size_t i=0; i<str_vec.size(); i++){
		ret = mkdir(dirname_tmp.c_str(), S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);  // create the output directory
		if(ret==-1 and errno!=EEXIST){
			cerr << __func__ << ", line=" << __LINE__ << ": cannot create directory: " << dirname_tmp << endl;
			exit(1);
		}
		if(i+1<str_vec.size()) dirname_tmp = dirname_tmp + "/" + str_vec.at(i+1);
	}
}

// get the number of high ratio indel bases
vector<double> getTotalHighIndelClipRatioBaseNum(Base *regBaseArr, int64_t arr_size){
	int32_t i, indel_num, clip_num, total_cov;
	double ratio, total, total2;
	Base *base;
	vector<double> base_num_vec;

	total = total2 = 0;
	for(i=0; i<arr_size; i++){
		base  = regBaseArr + i;
		indel_num = base->getTotalIndelNum();
		clip_num = base->getTotalClipNum();
		total_cov = base->getTotalCovNum();
		ratio = (double)(indel_num + clip_num) / total_cov;
		if(ratio>=HIGH_INDEL_CLIP_RATIO_THRES) total ++;
		//if(ratio>=SECOND_INDEL_CLIP_RATIO_THRES) total2 ++;
		if(ratio>=0.1f) total2 ++;
	}

	ratio = (double)total2 / arr_size;

	base_num_vec.push_back(total);
	base_num_vec.push_back(ratio);

	return base_num_vec;
}

// get mismatch regions
vector<mismatchReg_t*> getMismatchRegVec(localAln_t *local_aln){
	vector<mismatchReg_t*> misReg_vec;	// all the mismatch regions including gap regions
	int64_t ref_pos, local_ref_pos, query_pos, tmp, startMisRefPos, endMisRefPos, startMisLocalRefPos, endMisLocalRefPos, startMisQueryPos, endMisQueryPos;
	int32_t i, j, start_aln_idx, start_check_idx, end_check_idx, mis_reg_size, begin_mismatch_aln_idx, end_mismatch_aln_idx;
	string ctgseq_aln, midseq_aln, refseq_aln;
	mismatchReg_t *mis_reg, *mis_reg2;
	reg_t *reg;
	bool gap_flag, extend_flag, flag;
	int16_t mis_loc_flag;	// -1 for unused, 0 for gap in query, 1 for gap in reference, 2 for mismatch

	if(local_aln->reg==NULL or local_aln->start_aln_idx_var==-1 or local_aln->end_aln_idx_var==-1) return misReg_vec;

	ctgseq_aln = local_aln->alignResultVec[0];
	midseq_aln = local_aln->alignResultVec[1];
	refseq_aln = local_aln->alignResultVec[2];

	// compute refPos and queryPos at start_aln_idx
	start_check_idx = local_aln->start_aln_idx_var - 2* EXT_SIZE_CHK_VAR_LOC;
	end_check_idx = local_aln->end_aln_idx_var + 2 * EXT_SIZE_CHK_VAR_LOC;
	if(start_check_idx<0) start_check_idx = 0;
	if(end_check_idx>(int32_t)midseq_aln.size()-1) end_check_idx = midseq_aln.size() - 1;
	extend_flag = false;

	reg = local_aln->reg;
	ref_pos = reg->startRefPos;
	local_ref_pos = reg->startLocalRefPos;
	query_pos = reg->startQueryPos;

	// initialize the start positions
	start_aln_idx = local_aln->start_aln_idx_var;
	if(ctgseq_aln[start_aln_idx]=='-'){ // gap in query
		query_pos --;
		while(ctgseq_aln[start_aln_idx]=='-' and start_aln_idx>0){
			ref_pos --;
			local_ref_pos --;
			start_aln_idx --;
		}
	}else if(refseq_aln[start_aln_idx]=='-'){ // gap in reference
		while(refseq_aln[start_aln_idx]=='-' and start_aln_idx>0){
			query_pos --;
			start_aln_idx --;
		}
		ref_pos --;
		local_ref_pos --;
	}

	for(i=start_aln_idx; i>=start_check_idx; i--){
		if(midseq_aln[i]==' '){ // mismatch including gap
			if(ctgseq_aln[i]=='-'){ // gap in query
				ref_pos --;
				local_ref_pos --;
			}else if(refseq_aln[i]=='-'){ // gap in reference
				query_pos --;
			}else{ // mismatch
				query_pos --;
				ref_pos --;
				local_ref_pos --;
			}

			if(i==end_check_idx and extend_flag==false){
				end_check_idx -= EXT_SIZE_CHK_VAR_LOC;
				if(end_check_idx<0) end_check_idx = 0;
				extend_flag = true;
				//cout << "line=" << __LINE__ << ", end_check_idx=" << end_check_idx << endl;
			}
		}else{ // match
			if(extend_flag) break;

			query_pos --;
			ref_pos --;
			local_ref_pos --;
		}
	}
	if(i>=start_check_idx) {
		start_check_idx = i;
	}else{ // adjust to correct value
		start_check_idx = i + 1;
		query_pos ++;
		ref_pos ++;
		local_ref_pos ++;
	}

	startMisRefPos = endMisRefPos = startMisLocalRefPos = endMisLocalRefPos = startMisQueryPos = endMisQueryPos = -1;
	begin_mismatch_aln_idx = end_mismatch_aln_idx = -1;
	mis_reg_size = 0;
	gap_flag = extend_flag = false;
	mis_loc_flag = -1;
	for(i=start_check_idx; i<=end_check_idx; i++){
		if(midseq_aln[i]==' '){ // mismatch including gap
			if(begin_mismatch_aln_idx==-1){
				begin_mismatch_aln_idx = i;
				startMisRefPos = ref_pos;
				startMisLocalRefPos = local_ref_pos;
				startMisQueryPos = query_pos;

				if(ctgseq_aln[i]=='-'){ // gap in query
					mis_loc_flag = 0;
				}else if(refseq_aln[i]=='-'){ // gap in reference
					mis_loc_flag = 1;
				}else{ // mismatch
					mis_loc_flag = 2;
				}
				mis_reg_size = 1;
			}else{
				mis_reg_size ++;
			}

			if(ctgseq_aln[i]=='-'){ // gap in query
				ref_pos ++;
				local_ref_pos ++;
				gap_flag = true;
			}else if(refseq_aln[i]=='-'){ // gap in reference
				query_pos ++;
				gap_flag = true;
			}else{ // mismatch
				query_pos ++;
				ref_pos ++;
				local_ref_pos ++;
			}

			if(i==end_check_idx){
				if(extend_flag==false){
					end_check_idx += EXT_SIZE_CHK_VAR_LOC;
					if(end_check_idx>(int32_t)midseq_aln.size()-1) end_check_idx = midseq_aln.size() - 1;
					extend_flag = true;
					//cout << "line=" << __LINE__ << ", end_check_idx=" << end_check_idx << endl;
				}else if(begin_mismatch_aln_idx!=-1){ // mismatch region not closed
					end_check_idx += EXT_SIZE_CHK_VAR_LOC;
					if(end_check_idx>(int32_t)midseq_aln.size()-1) end_check_idx = midseq_aln.size() - 1;
					//cout << "line=" << __LINE__ << ", end_check_idx=" << end_check_idx << endl;
				}
			}
		}else{ // match
			if(begin_mismatch_aln_idx!=-1){
				end_mismatch_aln_idx = i - 1;
				switch(mis_loc_flag){
					case 0:  // gap in query
						endMisRefPos = ref_pos - 1;
						endMisLocalRefPos = local_ref_pos - 1;
						endMisQueryPos = query_pos;
						break;
					case 1:  // gap in reference
						endMisRefPos = ref_pos;
						endMisLocalRefPos = local_ref_pos;
						endMisQueryPos = query_pos - 1;
						break;
					case 2:  // mismatch
						endMisRefPos = ref_pos - 1;
						endMisLocalRefPos = local_ref_pos - 1;
						endMisQueryPos = query_pos - 1;
						break;
					default: // error
						cerr << __func__ << ", line=" << __LINE__ << ", invalid mis_loc_flag: " << mis_loc_flag << ", error!" << endl;
						exit(1);
				}

				mis_reg = new mismatchReg_t();
				mis_reg->start_aln_idx = begin_mismatch_aln_idx;
				mis_reg->end_aln_idx = end_mismatch_aln_idx;
				mis_reg->reg_size = mis_reg_size;
				mis_reg->startRefPos = startMisRefPos;
				mis_reg->endRefPos = endMisRefPos;
				mis_reg->startLocalRefPos = startMisLocalRefPos;
				mis_reg->endLocalRefPos = endMisLocalRefPos;
				mis_reg->startQueryPos = startMisQueryPos;
				mis_reg->endQueryPos = endMisQueryPos;
				mis_reg->valid_flag = true;
				mis_reg->gap_flag = gap_flag;
				misReg_vec.push_back(mis_reg);

				begin_mismatch_aln_idx = end_mismatch_aln_idx = -1;
				startMisRefPos = endMisRefPos = startMisLocalRefPos = endMisLocalRefPos = startMisQueryPos = endMisQueryPos = -1;
				mis_reg_size = 0;
				mis_loc_flag = -1;

				if(extend_flag) break;
			}
			if(extend_flag) break;

			query_pos ++;
			ref_pos ++;
			local_ref_pos ++;
			gap_flag = false;
		}
	}

	// sort
	for(i=0; i<(int32_t)misReg_vec.size(); i++){
		mis_reg = misReg_vec.at(i);
		for(j=i+1; j<(int32_t)misReg_vec.size(); j++){
			mis_reg2 = misReg_vec.at(j);
			if(mis_reg->start_aln_idx>mis_reg2->start_aln_idx){ // swap
				tmp = mis_reg->start_aln_idx; mis_reg->start_aln_idx = mis_reg2->start_aln_idx; mis_reg2->start_aln_idx = tmp;
				tmp = mis_reg->end_aln_idx; mis_reg->end_aln_idx = mis_reg2->end_aln_idx; mis_reg2->end_aln_idx = tmp;
				tmp = mis_reg->reg_size; mis_reg->reg_size = mis_reg2->reg_size; mis_reg2->reg_size = tmp;
				tmp = mis_reg->startRefPos; mis_reg->startRefPos = mis_reg2->startRefPos; mis_reg2->startRefPos = tmp;
				tmp = mis_reg->endRefPos; mis_reg->endRefPos = mis_reg2->endRefPos; mis_reg2->endRefPos = tmp;
				tmp = mis_reg->startLocalRefPos; mis_reg->startLocalRefPos = mis_reg2->startLocalRefPos; mis_reg2->startLocalRefPos = tmp;
				tmp = mis_reg->endLocalRefPos; mis_reg->endLocalRefPos = mis_reg2->endLocalRefPos; mis_reg2->endLocalRefPos = tmp;
				tmp = mis_reg->startQueryPos; mis_reg->startQueryPos = mis_reg2->startQueryPos; mis_reg2->startQueryPos = tmp;
				tmp = mis_reg->endQueryPos; mis_reg->endQueryPos = mis_reg2->endQueryPos; mis_reg2->endQueryPos = tmp;
				flag = mis_reg->valid_flag; mis_reg->valid_flag = mis_reg2->valid_flag; mis_reg2->valid_flag = flag;
				flag = mis_reg->gap_flag; mis_reg->gap_flag = mis_reg2->gap_flag; mis_reg2->gap_flag = flag;
			}
		}
	}

	return misReg_vec;
}

// get mismatch regions
vector<mismatchReg_t*> getMismatchRegVecWithoutPos(localAln_t *local_aln){
	vector<mismatchReg_t*> misReg_vec;	// all the mismatch regions including gap regions
	string ctgseq_aln, midseq_aln, refseq_aln;
	int32_t i, j, tmp, start_check_idx, end_check_idx, mis_reg_size, begin_mismatch_aln_idx, end_mismatch_aln_idx;
	mismatchReg_t *mis_reg, *mis_reg2;
	bool gap_flag, flag;

	ctgseq_aln = local_aln->alignResultVec[0];
	midseq_aln = local_aln->alignResultVec[1];
	refseq_aln = local_aln->alignResultVec[2];

	// compute refPos and queryPos at start_aln_idx
	start_check_idx = 0;
	end_check_idx = midseq_aln.size() - 1;

	begin_mismatch_aln_idx = end_mismatch_aln_idx = -1;
	mis_reg_size = 0;
	gap_flag = false;
	for(i=start_check_idx; i<=end_check_idx; i++){
		if(midseq_aln[i]==' '){ // mismatch including gap
			if(begin_mismatch_aln_idx==-1){
				begin_mismatch_aln_idx = i;
				mis_reg_size = 1;
			}else{
				mis_reg_size ++;
			}

			if(ctgseq_aln[i]=='-' or refseq_aln[i]=='-'){ // gap in query or reference
				gap_flag = true;
			}
		}else{ // match
			if(begin_mismatch_aln_idx!=-1){
				end_mismatch_aln_idx = i - 1;

				mis_reg = new mismatchReg_t();
				mis_reg->start_aln_idx = begin_mismatch_aln_idx;
				mis_reg->end_aln_idx = end_mismatch_aln_idx;
				mis_reg->reg_size = mis_reg_size;
				mis_reg->startRefPos = mis_reg->endRefPos = mis_reg->startLocalRefPos = mis_reg->endLocalRefPos = mis_reg->startQueryPos = mis_reg->endQueryPos = -1;
				mis_reg->valid_flag = true;
				mis_reg->gap_flag = gap_flag;
				misReg_vec.push_back(mis_reg);

				begin_mismatch_aln_idx = end_mismatch_aln_idx = -1;
				mis_reg_size = 0;
			}
			gap_flag = false;
		}
	}

	// sort
	for(i=0; i<(int32_t)misReg_vec.size(); i++){
		mis_reg = misReg_vec.at(i);
		for(j=i+1; j<(int32_t)misReg_vec.size(); j++){
			mis_reg2 = misReg_vec.at(j);
			if(mis_reg->start_aln_idx>mis_reg2->start_aln_idx){ // swap
				tmp = mis_reg->start_aln_idx; mis_reg->start_aln_idx = mis_reg2->start_aln_idx; mis_reg2->start_aln_idx = tmp;
				tmp = mis_reg->end_aln_idx; mis_reg->end_aln_idx = mis_reg2->end_aln_idx; mis_reg2->end_aln_idx = tmp;
				tmp = mis_reg->reg_size; mis_reg->reg_size = mis_reg2->reg_size; mis_reg2->reg_size = tmp;
				flag = mis_reg->valid_flag; mis_reg->valid_flag = mis_reg2->valid_flag; mis_reg2->valid_flag = flag;
				flag = mis_reg->gap_flag; mis_reg->gap_flag = mis_reg2->gap_flag; mis_reg2->gap_flag = flag;
			}
		}
	}

	return misReg_vec;
}

void removeShortPolymerMismatchRegItems(localAln_t *local_aln, vector<mismatchReg_t*> &misReg_vec, string &inBamFile, faidx_t *fai){
	size_t m;
	int64_t start_ref_pos, end_ref_pos, chrlen_tmp;
	int32_t i, j, start_check_idx, end_check_idx;
	string chrname_tmp, ctgseq_aln, refseq_aln, query_substr, subject_substr, substr, substr_gap;
	mismatchReg_t *mis_reg, *mis_reg0, *mis_reg2, *tmp_mis_reg;
	char ch_left, ch_right;
	bool flag, large_neighbor_dist_flag, high_con_ratio_flag, disagree_flag;
	vector<bam1_t*> alnDataVector;
	Base *start_base, *end_base;

	if(misReg_vec.empty()) return;

	if(local_aln->reg) chrname_tmp = local_aln->reg->chrname;
	else if(local_aln->cand_reg) chrname_tmp = local_aln->cand_reg->chrname;
	else return;
	chrlen_tmp = faidx_seq_len(fai, chrname_tmp.c_str()); // get the reference length

	start_ref_pos = misReg_vec.at(0)->startRefPos - 2 * EXT_SIZE_CHK_VAR_LOC;
	end_ref_pos = misReg_vec.at(misReg_vec.size()-1)->endRefPos + 2* EXT_SIZE_CHK_VAR_LOC;
	if(start_ref_pos<1) start_ref_pos = 1;
	if(end_ref_pos>chrlen_tmp) end_ref_pos = chrlen_tmp;

	// load align data
	alnDataLoader data_loader(chrname_tmp, start_ref_pos, end_ref_pos, inBamFile);
	data_loader.loadAlnData(alnDataVector);

	// load coverage
	covLoader cov_loader(chrname_tmp, start_ref_pos, end_ref_pos, fai);
	Base *baseArray = cov_loader.initBaseArray();
	cov_loader.generateBaseCoverage(baseArray, alnDataVector);

	ctgseq_aln = local_aln->alignResultVec[0];
	refseq_aln = local_aln->alignResultVec[2];

	// filter out short (e.g. <= 1 bp) mismatched polymer items
	for(i=0; i<(int32_t)misReg_vec.size(); i++){
		mis_reg = misReg_vec.at(i);
		start_check_idx = mis_reg->start_aln_idx;
		end_check_idx = mis_reg->end_aln_idx;

		// determine whether the base contains many insertions or deletions
		start_base = baseArray + mis_reg->startRefPos - start_ref_pos;
		end_base = baseArray + mis_reg->endRefPos - start_ref_pos;
		high_con_ratio_flag = disagree_flag = false;

		if(((start_base->max_con_type==BAM_CINS and start_base->maxConIndelEventRatio>=MIN_HIGH_CONSENSUS_INS_RATIO) or (start_base->max_con_type==BAM_CDEL and start_base->maxConIndelEventRatio>=MIN_HIGH_CONSENSUS_DEL_RATIO)) or ((end_base->max_con_type==BAM_CINS and end_base->maxConIndelEventRatio>=MIN_HIGH_CONSENSUS_INS_RATIO) or (end_base->max_con_type==BAM_CDEL and end_base->maxConIndelEventRatio>=MIN_HIGH_CONSENSUS_DEL_RATIO)))
			high_con_ratio_flag = true;
		if(start_base->isDisagreeBase() or end_base->isDisagreeBase())
			disagree_flag = true;

		if(high_con_ratio_flag or disagree_flag) continue; // skip below operations

		if(mis_reg->gap_flag){
			query_substr = ctgseq_aln.substr(start_check_idx, end_check_idx-start_check_idx+1);
			subject_substr = refseq_aln.substr(start_check_idx, end_check_idx-start_check_idx+1);

			// check polymer
			if(query_substr.find("-")!=string::npos){ // gap in query
				substr_gap = query_substr;
				substr = subject_substr;
				ch_left = refseq_aln.at(start_check_idx-1);
				ch_right = refseq_aln.at(end_check_idx+1);
			}else{ // gap in subject
				substr = query_substr;
				substr_gap = subject_substr;
				ch_left = ctgseq_aln.at(start_check_idx-1);
				ch_right = ctgseq_aln.at(end_check_idx+1);
			}

			flag = isPolymerSeq(substr);
			if(flag){ // check
				if(high_con_ratio_flag==false and disagree_flag==false and (substr.at(0)==ch_left or substr.at(substr.size()-1)==ch_right) and substr.size()<MIN_VALID_POLYMER_SIZE)
					mis_reg->valid_flag = false;
			}
		}

		// further to check the distance of its neighbors
		large_neighbor_dist_flag = false;
		if(mis_reg->valid_flag and mis_reg->reg_size<=MAX_SHORT_MIS_REG_SIZE){
			if(i==0){ // first item
				mis_reg2 = NULL;
				for(m=i+1; m<misReg_vec.size(); m++){ // search right side
					tmp_mis_reg = misReg_vec.at(m);
					if(tmp_mis_reg->valid_flag){
						mis_reg2 = tmp_mis_reg;
						break;
					}
				}
				if(mis_reg2==NULL or (mis_reg2->start_aln_idx-mis_reg->end_aln_idx>MIN_ADJACENT_REG_DIST))
					large_neighbor_dist_flag = true;
			}else if(i==(int32_t)misReg_vec.size()-1){ // last item
				mis_reg0 = NULL;
				for(j=i-1; j>=0; j--){ // search left side
					tmp_mis_reg = misReg_vec.at(j);
					if(tmp_mis_reg->valid_flag){
						mis_reg0 = tmp_mis_reg;
						break;
					}
				}
				if(mis_reg0==NULL or (mis_reg->start_aln_idx-mis_reg0->end_aln_idx>MIN_ADJACENT_REG_DIST))
					large_neighbor_dist_flag = true;
			}else{ // middle items
				mis_reg0 = NULL;
				for(j=i-1; j>=0; j--){ // search left side
					tmp_mis_reg = misReg_vec.at(j);
					if(tmp_mis_reg->valid_flag){
						mis_reg0 = tmp_mis_reg;
						break;
					}
				}
				mis_reg2 = NULL;
				for(m=i+1; m<misReg_vec.size(); m++){ // search right side
					tmp_mis_reg = misReg_vec.at(m);
					if(tmp_mis_reg->valid_flag){
						mis_reg2 = tmp_mis_reg;
						break;
					}
				}
				if((mis_reg0==NULL or (mis_reg->start_aln_idx-mis_reg0->end_aln_idx>MIN_ADJACENT_REG_DIST)) and (mis_reg2==NULL or (mis_reg2->start_aln_idx-mis_reg->end_aln_idx>MIN_ADJACENT_REG_DIST)))
					large_neighbor_dist_flag = true;
			}

			if(large_neighbor_dist_flag==true and high_con_ratio_flag==false and disagree_flag==false)
				mis_reg->valid_flag = false;
		}
	}

	// remove invalid items
	for(m=0; m<misReg_vec.size(); ){
		mis_reg = misReg_vec.at(m);
		if(mis_reg->valid_flag==false){
			delete mis_reg;
			misReg_vec.erase(misReg_vec.begin()+m);
		}else m++;
	}

	// release memory
	data_loader.freeAlnData(alnDataVector);
	cov_loader.freeBaseArray(baseArray);
}

void adjustVarLocByMismatchRegs(reg_t *reg, vector<mismatchReg_t*> &misReg_vec, int32_t start_aln_idx_var, int32_t end_aln_idx_var){
	int32_t i, max_reg_idx, maxValue, secMax_reg_idx, secMaxValue, left_reg_idx, right_reg_idx;
	mismatchReg_t *mis_reg, *mis_reg2;

	if(start_aln_idx_var!=-1 and end_aln_idx_var!=-1){
		// get maximum and second maximum mismatch region
		max_reg_idx = secMax_reg_idx = -1;
		maxValue = secMaxValue = 0;
		for(i=0; i<(int32_t)misReg_vec.size(); i++){
			mis_reg = misReg_vec.at(i);
			if(mis_reg->reg_size>maxValue){
				max_reg_idx = i;
				maxValue = mis_reg->reg_size;
			}
			if(mis_reg->reg_size>=secMaxValue){
				secMax_reg_idx = i;
				secMaxValue = mis_reg->reg_size;
			}
		}

		// check its neighbors
		if(max_reg_idx!=-1){
			left_reg_idx = right_reg_idx = -1;

			// left side extend
			mis_reg = misReg_vec.at(max_reg_idx);
			for(i=max_reg_idx-1; i>=0; i--){
				mis_reg2 = misReg_vec.at(i);
				//if((mis_reg->start_aln_idx-mis_reg2->end_aln_idx<=MIN_ADJACENT_REG_DIST) or (mis_reg->start_aln_idx>=end_aln_idx_var and mis_reg2->end_aln_idx<=start_aln_idx_var)){ // distance less than threshold
				if((mis_reg->start_aln_idx-mis_reg2->end_aln_idx<=MIN_ADJACENT_REG_DIST) and (mis_reg->start_aln_idx>=end_aln_idx_var and mis_reg2->end_aln_idx<=start_aln_idx_var)){ // distance less than threshold
					left_reg_idx = i;
					mis_reg = mis_reg2;
					mis_reg2 = NULL;
				}else break;
			}
			if(left_reg_idx==-1) left_reg_idx = max_reg_idx;

			// right side extend
			mis_reg = misReg_vec.at(max_reg_idx);
			for(i=max_reg_idx+1; i<(int32_t)misReg_vec.size(); i++){
				mis_reg2 = misReg_vec.at(i);
				//if((mis_reg2->start_aln_idx-mis_reg->end_aln_idx<=MIN_ADJACENT_REG_DIST) or (mis_reg->end_aln_idx<=start_aln_idx_var and mis_reg2->start_aln_idx>=end_aln_idx_var)){ // distance less than threshold
				if((mis_reg2->start_aln_idx-mis_reg->end_aln_idx<=MIN_ADJACENT_REG_DIST) and (mis_reg->end_aln_idx<=start_aln_idx_var and mis_reg2->start_aln_idx>=end_aln_idx_var)){ // distance less than threshold
					right_reg_idx = i;
					mis_reg = mis_reg2;
					mis_reg2 = NULL;
				}else break;
			}
			if(right_reg_idx==-1) right_reg_idx = max_reg_idx;

			mis_reg = misReg_vec.at(left_reg_idx);
			mis_reg2 = misReg_vec.at(right_reg_idx);
			reg->startRefPos = mis_reg->startRefPos - 1;
			reg->endRefPos = mis_reg2->endRefPos;
			reg->startLocalRefPos = mis_reg->startLocalRefPos - 1;
			reg->endLocalRefPos = mis_reg2->endLocalRefPos;
			reg->startQueryPos = mis_reg->startQueryPos - 1;
			reg->endQueryPos = mis_reg2->endQueryPos;
		}
	}
}

void releaseMismatchRegVec(vector<mismatchReg_t*> &misReg_vec){
	// release memory
	for(size_t i=0; i<misReg_vec.size(); i++) delete misReg_vec.at(i);
	vector<mismatchReg_t*>().swap(misReg_vec);
}

mismatchReg_t *getMismatchReg(int32_t aln_idx, vector<mismatchReg_t*> &misReg_vec){
	mismatchReg_t *mis_item_ret = NULL, *mis_item;
	for(size_t i=0; i<misReg_vec.size(); i++){
		mis_item = misReg_vec.at(i);
		if(mis_item->valid_flag){
			if(aln_idx>=mis_item->start_aln_idx and aln_idx<=mis_item->end_aln_idx){
				mis_item_ret = mis_item;
				break;
			}
		}
	}
	return mis_item_ret;
}

bool isPolymerSeq(string &seq){
	bool flag = true;
	for(size_t m=1; m<seq.size(); m++){
		if(seq.at(m-1)!=seq.at(m)){
			flag = false;
			break;
		}
	}
	return flag;
}

// get query clip align segments
vector<clipAlnData_t*> getQueryClipAlnSegs(string &queryname, vector<clipAlnData_t*> &clipAlnDataVector){
	vector<clipAlnData_t*> query_aln_segs;
	for(size_t i=0; i<clipAlnDataVector.size(); i++)
		if(clipAlnDataVector.at(i)->query_checked_flag==false and clipAlnDataVector.at(i)->queryname.compare(queryname)==0 and clipAlnDataVector.at(i)->bam!=NULL){
			query_aln_segs.push_back(clipAlnDataVector.at(i));
		}

	return query_aln_segs;
}

// get query clip align segments
vector<clipAlnData_t*> getQueryClipAlnSegsAll(string &queryname, vector<clipAlnData_t*> &clipAlnDataVector){
	vector<clipAlnData_t*> query_aln_segs;
	for(size_t i=0; i<clipAlnDataVector.size(); i++)
		if(clipAlnDataVector.at(i)->query_checked_flag==false and clipAlnDataVector.at(i)->queryname.compare(queryname)==0){
			query_aln_segs.push_back(clipAlnDataVector.at(i));
		}

	return query_aln_segs;
}

bool isQuerySelfOverlap(vector<clipAlnData_t*> &query_aln_segs, int32_t maxVarRegSize){
	bool flag;
	clipAlnData_t *clip_aln1, *clip_aln2;

	flag = false;
	for(size_t i=0; i<query_aln_segs.size()-1; i++){
		clip_aln1 = query_aln_segs.at(i);
		for(size_t j=i+1; j<query_aln_segs.size(); j++){
			clip_aln2 = query_aln_segs.at(j);
			flag = isSegSelfOverlap(clip_aln1, clip_aln2, maxVarRegSize);
			if(flag) break;
		}
		if(flag) break;
	}

	return flag;
}

bool isSegSelfOverlap(clipAlnData_t *clip_aln1, clipAlnData_t *clip_aln2, int32_t maxVarRegSize){
	bool flag = false;
	int32_t overlap_size;

	if(clip_aln1 and clip_aln2){
		if(clip_aln1->chrname.compare(clip_aln2->chrname)==0){
			if(isOverlappedPos(clip_aln1->startRefPos, clip_aln1->endRefPos, clip_aln2->startRefPos, clip_aln2->endRefPos)){
				overlap_size = getOverlapSize(clip_aln1->startRefPos, clip_aln1->endRefPos, clip_aln2->startRefPos, clip_aln2->endRefPos);
				if(overlap_size>=MIN_QUERY_SELF_OVERLAP_SIZE)
					flag = true;
			}
//			else if(clip_aln2->startRefPos<clip_aln1->endRefPos){
//				overlap_size = getOverlapSize(clip_aln1->startRefPos, clip_aln1->endRefPos, clip_aln2->startRefPos, clip_aln2->endRefPos);
//				if(overlap_size>=-maxVarRegSize)
//					flag = true;
//			}
		}
	}

	return flag;
}

int32_t getVecIdxClipAlnData(clipAlnData_t *clip_aln, vector<clipAlnData_t*> &clipAlnDataVector){
	int32_t idx = -1;
	clipAlnData_t *clip_aln_tmp;

	for(size_t i=0; i<clipAlnDataVector.size(); i++){
		clip_aln_tmp = clipAlnDataVector.at(i);
		if(clip_aln_tmp==clip_aln){
			idx = i;
			break;
		}
	}

	return idx;
}

// get adjacent clip segment according to query positions
vector<int32_t> getAdjacentClipAlnSeg(int32_t arr_idx, int32_t clip_end_flag, vector<clipAlnData_t*> &query_aln_segs, int32_t minClipEndSize, int32_t max_reg_size){
	clipAlnData_t *clip_aln, *clip_aln_based, *clip_aln_other_end;
	int32_t clip_pos_based, clip_ref_pos_based, mate_clip_end, dist, ref_dist, min_dist, min_ref_dist, idx_min_dist, end_flag, min_dist_same_chr, min_ref_dist_same_chr, idx_min_dist_same_chr, end_flag_same_chr;
	string pre_chrname;
	int32_t qpos_increase_direction; // -1 for unused, 0 for increase, 1 for decrease
	vector<int32_t> adjClipAlnSegInfo; // [0]: array index of minimal distance; [1]: segment end flag; [2]: minimal query distance; [3]: minimal reference distance; [4]: array index of minimal distance in the same chromosome; [5]: segment end flag in the same chromosome; [6]: minimal query distance in the same chromosome; [7]: minimal reference distance in the same chromosome

	if(arr_idx!=-1){ // valid based segment
		clip_aln_based = query_aln_segs.at(arr_idx);
		qpos_increase_direction = -1;
		if(clip_end_flag==LEFT_END) {
			clip_pos_based = clip_aln_based->startQueryPos;
			clip_ref_pos_based = clip_aln_based->startRefPos;
			clip_aln_other_end = clip_aln_based->right_aln;
			if(clip_aln_based->aln_orient==ALN_PLUS_ORIENT)
				qpos_increase_direction = 1; // decrease direction
			else
				qpos_increase_direction = 0; // increase direction
		}else {
			clip_pos_based = clip_aln_based->endQueryPos;
			clip_ref_pos_based = clip_aln_based->endRefPos;
			clip_aln_other_end = clip_aln_based->left_aln;
			if(clip_aln_based->aln_orient==ALN_PLUS_ORIENT)
				qpos_increase_direction = 0; // increase direction
			else
				qpos_increase_direction = 1; // decrease direction
		}

		pre_chrname = clip_aln_based->chrname;
		min_dist = min_ref_dist = min_dist_same_chr = min_ref_dist_same_chr = INT_MAX;
		idx_min_dist = idx_min_dist_same_chr = -1;
		end_flag = end_flag_same_chr = -1;
		for(size_t i=0; i<query_aln_segs.size(); i++){
			if(i!=(size_t)arr_idx){
				clip_aln = query_aln_segs.at(i);

				// determine the valid mate end
				if(clip_aln->aln_orient==clip_aln_based->aln_orient){
					if(clip_end_flag==LEFT_END) mate_clip_end = RIGHT_END;
					else mate_clip_end = LEFT_END;
				}else {
					if(clip_end_flag==LEFT_END) mate_clip_end = LEFT_END;
					else mate_clip_end = RIGHT_END;
				}

				// left end
				if(mate_clip_end==LEFT_END and clip_aln->left_clip_checked_flag==false and clip_aln->leftClipSize>=minClipEndSize){
					dist = clip_aln->startQueryPos - clip_pos_based;
					ref_dist = clip_aln->startRefPos - clip_ref_pos_based;
					if(clip_aln->aln_orient==ALN_MINUS_ORIENT) dist = -dist;
					if(abs(dist)<abs(min_dist)) {
						min_dist = dist;
						min_ref_dist = ref_dist;
						idx_min_dist = i;
						end_flag = LEFT_END;
					}
					if(abs(ref_dist)<=max_reg_size and abs(dist)<abs(min_dist_same_chr) and clip_aln->chrname.compare(pre_chrname)==0 and clip_aln->aln_orient==clip_aln_based->aln_orient){
						min_dist_same_chr = dist;
						min_ref_dist_same_chr = ref_dist;
						idx_min_dist_same_chr = i;
						end_flag_same_chr = LEFT_END;
					}
				}
				// right end
				if(mate_clip_end==RIGHT_END and clip_aln->right_clip_checked_flag==false and clip_aln->rightClipSize>=minClipEndSize){
					dist = clip_aln->endQueryPos - clip_pos_based;
					ref_dist = clip_aln->endRefPos - clip_ref_pos_based;
					if(clip_aln->aln_orient==ALN_MINUS_ORIENT) dist = -dist;
					if(abs(dist)<abs(min_dist)) {
						min_dist = dist;
						min_ref_dist = ref_dist;
						idx_min_dist = i;
						end_flag = RIGHT_END;
					}
					if(abs(ref_dist)<=max_reg_size and abs(dist)<abs(min_dist_same_chr) and clip_aln->chrname.compare(pre_chrname)==0 and clip_aln->aln_orient==clip_aln_based->aln_orient){
						min_dist_same_chr = dist;
						min_ref_dist_same_chr = ref_dist;
						idx_min_dist_same_chr = i;
						end_flag_same_chr = RIGHT_END;
					}
				}
			}
		}

		// confirm
		if(idx_min_dist!=-1){
			clip_aln = query_aln_segs.at(idx_min_dist);
			if(clip_aln==clip_aln_other_end){ // invalid if two ends of a item point to the same other item
				idx_min_dist = min_dist = min_ref_dist = end_flag = -1;
				idx_min_dist_same_chr = min_dist_same_chr = min_ref_dist_same_chr = end_flag_same_chr = -1;
			}else{
				// make sure the increase_direction of the position is consist with the clipping segments
				if(qpos_increase_direction==0){ // increase direction
					if((clip_aln->aln_orient==ALN_PLUS_ORIENT and clip_aln->endQueryPos<clip_pos_based) or (clip_aln->aln_orient==ALN_MINUS_ORIENT and clip_aln->startQueryPos<clip_pos_based)){
						idx_min_dist = min_dist = min_ref_dist = end_flag = -1;
						idx_min_dist_same_chr = min_dist_same_chr = min_ref_dist_same_chr = end_flag_same_chr = -1;
					}
				}else{ // decrease direction
					if((clip_aln->aln_orient==ALN_PLUS_ORIENT and clip_aln->startQueryPos>clip_pos_based) or (clip_aln->aln_orient==ALN_MINUS_ORIENT and clip_aln->endQueryPos>clip_pos_based)){
						idx_min_dist = min_dist = min_ref_dist = end_flag = -1;
						idx_min_dist_same_chr = min_dist_same_chr = min_ref_dist_same_chr = end_flag_same_chr = -1;
					}
				}
			}
			//if(idx_min_dist!=idx_min_dist_same_chr) min_ref_dist = -1;
		}else{
			idx_min_dist = min_dist = min_ref_dist = end_flag = -1;
			idx_min_dist_same_chr = min_dist_same_chr = min_ref_dist_same_chr = end_flag_same_chr = -1;
		}
	}else{ // invalid based segment
		idx_min_dist = min_dist = min_ref_dist = end_flag = -1;
		idx_min_dist_same_chr = min_dist_same_chr = min_ref_dist_same_chr = end_flag_same_chr = -1;
	}

	// save to vector
	adjClipAlnSegInfo.push_back(idx_min_dist);
	adjClipAlnSegInfo.push_back(end_flag);
	adjClipAlnSegInfo.push_back(min_dist);
	adjClipAlnSegInfo.push_back(min_ref_dist);
	adjClipAlnSegInfo.push_back(idx_min_dist_same_chr);
	adjClipAlnSegInfo.push_back(end_flag_same_chr);
	adjClipAlnSegInfo.push_back(min_dist_same_chr);
	adjClipAlnSegInfo.push_back(min_ref_dist_same_chr);

	return adjClipAlnSegInfo;
}

// get adjacent clip segment according to query positions
vector<int32_t> getAdjacentClipAlnSeg(clipAlnData_t *clip_aln, int32_t clip_end_flag, vector<clipAlnData_t*> &query_aln_segs, int32_t minClipEndSize, int32_t max_reg_size){
	vector<int32_t> adjClipAlnSegInfo; // [0]: array index of minimal distance; [1]: segment end flag; [2]: minimal query distance; [3]: minimal reference distance; [4]: array index of minimal distance in the same chromosome; [5]: segment end flag in the same chromosome; [6]: minimal query distance in the same chromosome; [7]: minimal reference distance in the same chromosome
	int32_t idx_vec;

	idx_vec = getVecIdxClipAlnData(clip_aln, query_aln_segs);
	adjClipAlnSegInfo = getAdjacentClipAlnSeg(idx_vec, clip_end_flag, query_aln_segs, minClipEndSize, max_reg_size);

	return adjClipAlnSegInfo;
}

// get adjacent clip segment with prefered same chromosome according to query positions
clipAlnData_t* getAdjacentClipAlnSegPreferSameChr(clipAlnData_t *clip_aln, int32_t clip_end_flag, vector<clipAlnData_t*> &query_aln_segs, int32_t minClipEndSize, int32_t max_reg_size){
	clipAlnData_t *adj_clip_aln = NULL;
	vector<int32_t> adjClipAlnSegInfo; // [0]: array index of minimal distance; [1]: segment end flag; [2]: minimal query distance; [3]: minimal reference distance; [4]: array index of minimal distance in the same chromosome; [5]: segment end flag in the same chromosome; [6]: minimal query distance in the same chromosome; [7]: minimal reference distance in the same chromosome
	int32_t idx_vec, mate_arr_idx;
	int32_t mate_arr_idx_same_chr, min_ref_dist_same_chr;

	idx_vec = getVecIdxClipAlnData(clip_aln, query_aln_segs);
	adjClipAlnSegInfo = getAdjacentClipAlnSeg(idx_vec, clip_end_flag, query_aln_segs, minClipEndSize, max_reg_size);

	mate_arr_idx = adjClipAlnSegInfo.at(0);
	mate_arr_idx_same_chr = adjClipAlnSegInfo.at(4);
	min_ref_dist_same_chr = adjClipAlnSegInfo.at(7);
	if(mate_arr_idx!=-1){ // mated
		// prefer the segments on the same chromosome with near reference distance
		if(mate_arr_idx_same_chr!=-1 and mate_arr_idx!=mate_arr_idx_same_chr and abs(min_ref_dist_same_chr)<=MAX_REF_DIST_SAME_CHR) // same chromosome with near reference distance
			mate_arr_idx = mate_arr_idx_same_chr;
		adj_clip_aln = query_aln_segs.at(mate_arr_idx);
	}

	return adj_clip_aln;
}

// generate BND items
vector<BND_t*> generateBNDItems(int32_t reg_id, int32_t clip_end, int32_t checked_arr[][2], string &chrname1, string &chrname2, int64_t tra_pos_arr[4], vector<string> &bnd_str_vec, faidx_t *fai){
	vector<BND_t*> bnd_vec;
	vector<string> bnd_str_vec2, bnd_str_vec2_tmp, mate_str_vec, mate_str_vec2;
	BND_t *bnd_item, *mate_bnd_item;
	int64_t bnd_pos, mate_bnd_pos, mate_dist;
	string bnd_id_str, mate_bnd_id_str, reg_str, mate_reg_str, mate_bnd_str, mate_bnd_str2, seq_str, seq2_str, support_num_str, cov_num_str;
	string mate_reg_str2;
	char *seq, *seq2, mate_orient_ch, mate_orient_ch2;
	int32_t mate_reg_id, mate_reg_id2, mate_clip_end, sub_vec_id, mate_sub_vec_id, seq_len, seq_len2, support_num, cov_num, support_num2, cov_num2;
	string id, ref, alt, qual, filter, info, format, format_val_str, sv_type, line_vcf;

	if(reg_id<0 or reg_id>=4){
		cerr << __func__ << ", line=" << __LINE__ << ": invalid region id: " << reg_id << ", which should be the intenger of [0,3], error." << endl;
		exit(1);
	}

	if(tra_pos_arr[reg_id]!=-1 or bnd_str_vec.at(reg_id).compare("-")!=0) {
		sub_vec_id = (clip_end + 1) % 2;
		if(checked_arr[reg_id][sub_vec_id]==0){ // unprocessed, then process it
			if(clip_end==RIGHT_END) bnd_pos = tra_pos_arr[reg_id] - 1;   // POS
			else bnd_pos = tra_pos_arr[reg_id];

			if(bnd_pos<=0)
				cout << "line=" << __LINE__ << ", bnd_pos=" << bnd_pos << ", error." << endl;

			reg_str = chrname1 + ":" + to_string(bnd_pos) + "-" + to_string(bnd_pos);
			seq = fai_fetch(fai, reg_str.c_str(), &seq_len);
			seq_str = seq;
			free(seq);

//			if(bnd_str_vec.at(reg_id).compare("-")!=0){
				bnd_str_vec2 = split(bnd_str_vec.at(reg_id), ",");
//			}else{
//				cerr << __func__ << ", line=" << __LINE__ << ": invalid BND information:" << bnd_str_vec.at(reg_id) << endl;
//				exit(1);
//			}

			mate_reg_str = bnd_str_vec2.at(sub_vec_id); // "2+|.|.|."
			if(mate_reg_str.compare("-")!=0){
				bnd_item = new struct BND_Node();
				mate_bnd_item = new struct BND_Node();

				mate_str_vec = split(mate_reg_str, "|");
				mate_reg_id = stoi(mate_str_vec.at(0).substr(0, 1));
				mate_orient_ch = mate_str_vec.at(0).at(1);
				support_num_str = mate_str_vec.at(2);
				cov_num_str = mate_str_vec.at(3);
				support_num = stoi(support_num_str);
				cov_num = stoi(cov_num_str);

				mate_bnd_pos = tra_pos_arr[mate_reg_id];
				if(mate_bnd_pos==-1){
					cout << chrname1 << ":" << bnd_pos << " <---> " << chrname2 << ":" << mate_bnd_pos << endl;
				}

				checked_arr[reg_id][sub_vec_id] = 1;
				if(clip_end==RIGHT_END){ // right end
					if(mate_orient_ch=='+'){ // same orient
						mate_clip_end = LEFT_END;
						reg_str = chrname2 + ":" + to_string(mate_bnd_pos) + "-" + to_string(mate_bnd_pos);
						seq2 = fai_fetch(fai, reg_str.c_str(), &seq_len2);

						seq2_str = seq2;
						mate_bnd_str = seq_str + "[" + chrname2 + ":" + to_string(mate_bnd_pos) + "[";
						mate_bnd_str2 = "]" + chrname1 + ":" + to_string(bnd_pos) + "]" + seq2_str;
					}else{ // different orient
						mate_clip_end = RIGHT_END;
						mate_bnd_pos --;
						reg_str = chrname2 + ":" + to_string(mate_bnd_pos) + "-" + to_string(mate_bnd_pos);
						seq2 = fai_fetch(fai, reg_str.c_str(), &seq_len2);

						seq2_str = seq2;
						mate_bnd_str = seq_str + "]" + chrname2 + ":" + to_string(mate_bnd_pos) + "]";
						mate_bnd_str2 =  "[" + chrname1 + ":" + to_string(bnd_pos) + "[" + seq2_str;
					}
				}else{ // left end
					if(mate_orient_ch=='+'){ // same orient
						mate_clip_end = RIGHT_END;
						mate_bnd_pos --;
						reg_str = chrname2 + ":" + to_string(mate_bnd_pos) + "-" + to_string(mate_bnd_pos);
						seq2 = fai_fetch(fai, reg_str.c_str(), &seq_len2);

						seq2_str = seq2;
						mate_bnd_str = "]" + chrname2 + ":" + to_string(mate_bnd_pos) + "]" + seq_str;
						mate_bnd_str2 = seq2_str + "[" + chrname1 + ":" + to_string(bnd_pos) + "[";
						//checked_arr[mate_reg_id][0] = 1;
					}else{ // different orient
						mate_clip_end = LEFT_END;
						reg_str = chrname2 + ":" + to_string(mate_bnd_pos) + "-" + to_string(mate_bnd_pos);
						seq2 = fai_fetch(fai, reg_str.c_str(), &seq_len2);

						seq2_str = seq2;
						mate_bnd_str =  "[" + chrname2 + ":" + to_string(mate_bnd_pos) + "[" + seq_str;
						mate_bnd_str2 = seq2_str + "[" + chrname1 + ":" + to_string(bnd_pos) + "[";
						//checked_arr[mate_reg_id][1] = 1;
					}
				}
				free(seq2);
				mate_sub_vec_id = (mate_clip_end + 1) % 2;
				if(checked_arr[mate_reg_id][mate_sub_vec_id]==0)
					checked_arr[mate_reg_id][mate_sub_vec_id] = 1;
				else{
					cerr << __func__ << ", line=" << __LINE__ << ": the BND information " << mate_bnd_str << " should not be processed twice, error." << endl;
					exit(1);
				}

				// mate coverage information
				if(bnd_str_vec.at(mate_reg_id).compare("-")!=0){
					bnd_str_vec2_tmp = split(bnd_str_vec.at(mate_reg_id), ",");
				}else{
					cerr << __func__ << ", line=" << __LINE__ << ": invalid BND information:" << bnd_str_vec.at(reg_id) << endl;
					exit(1);
				}

				mate_reg_str2 = bnd_str_vec2_tmp.at(mate_sub_vec_id); // "0+|..|..|.."
				if(mate_reg_str2.compare("-")!=0){
					mate_str_vec2 = split(mate_reg_str2, "|");
					mate_reg_id2 = stoi(mate_str_vec2.at(0).substr(0, 1));
					mate_orient_ch2 = mate_str_vec2.at(0).at(1);
					support_num2 = stoi(mate_str_vec.at(2));
					cov_num2 = stoi(mate_str_vec.at(3));

					if(mate_reg_id2!=reg_id){
						cerr << __func__ << ", line=" << __LINE__ << ": two regions do not match, reg_id=" << reg_id << ", mate_reg_id2=" << mate_reg_id2 << ", error!" << endl;
						exit(1);
					}else if(mate_orient_ch2!=mate_orient_ch){
						cerr << __func__ << ", line=" << __LINE__ << ": orientation of two regions do not match, mate_orient=" << mate_orient_ch << ", mate_orient2=" << mate_orient_ch2 << ", error!" << endl;
						exit(1);
					}
				}else{
					cerr << __func__ << ", line=" << __LINE__ << ": invalid BND information, error." << endl;
					exit(1);
				}

				bnd_id_str = "BND." + chrname1 + ":" + to_string(bnd_pos) + "-" + chrname2 + ":" + to_string(mate_bnd_pos);
				mate_bnd_id_str = "BND." + chrname2 + ":" + to_string(mate_bnd_pos) + "-" + chrname1 + ":" + to_string(bnd_pos);
				mate_dist = -1;
				if(chrname1.compare(chrname2)==0){
					mate_dist = mate_bnd_pos - bnd_pos;
					if(mate_dist<0) mate_dist = -mate_dist;
				}

				// bnd item
				id = bnd_id_str;			// ID
				ref = seq_str;					// REF
				alt = mate_bnd_str;			// ALT
				qual = ".";					// QUAL
				filter = "PASS";			// FILTER

				sv_type = VAR_BND_STR;
				format = "GT:AD:DP";		// FORMAT and the values
				format_val_str = "./.";
				format_val_str += ":" + support_num_str + "," + to_string(cov_num-support_num);
				format_val_str += ":" + cov_num_str;

				info = "SVTYPE=" + sv_type + ";" + "MATEID=" + mate_bnd_id_str;	// INFO: SVTYPE=sv_type;MATEID=mate_bnd_id
				if(chrname1.compare(chrname2)==0) info += ";MATEDIST=" + to_string(mate_dist);  // INFO: MATEDIST=mate_dist
				line_vcf = chrname1 + "\t" + to_string(bnd_pos) + "\t" + id + "\t" + ref + "\t" + alt + "\t" + qual + "\t" + filter + "\t" + info + "\t" + format + "\t" + format_val_str;

				bnd_item->reg_id = reg_id;
				bnd_item->clip_end = clip_end;
				bnd_item->mate_reg_id = mate_reg_id;
				bnd_item->mate_clip_end = mate_clip_end;
				bnd_item->mate_orient_ch = mate_orient_ch;
				bnd_item->chrname = chrname1;
				bnd_item->mate_chrname = chrname2;
				bnd_item->bnd_pos = bnd_pos;
				bnd_item->mate_bnd_pos = mate_bnd_pos;
				bnd_item->seq = seq_str;
				bnd_item->bnd_str = mate_bnd_str;
				bnd_item->vcf_line = line_vcf;
				bnd_item->support_num = support_num;
				bnd_item->mate_support_num = support_num2;
				bnd_item->cov_num = cov_num;
				bnd_item->mate_cov_num = cov_num2;
				bnd_item->mate_node = mate_bnd_item;

				// mate bnd item
				id = mate_bnd_id_str;		// ID
				ref = seq2_str;					// REF
				alt = mate_bnd_str2;		// ALT
				qual = ".";					// QUAL
				filter = "PASS";			// FILTER

				sv_type = VAR_BND_STR;
				format = "GT:AD:DP";		// FORMAT and the values
				format_val_str = "./.";
				format_val_str += ":" + to_string(support_num2) + "," + to_string(cov_num2-support_num2);
				format_val_str += ":" + to_string(cov_num2);

				info = "SVTYPE=" + sv_type + ";" + "MATEID=" + bnd_id_str;	// INFO: SVTYPE=sv_type;MATEID=mate_bnd_id
				if(chrname1.compare(chrname2)==0) info += ";MATEDIST=" + to_string(mate_dist);  // INFO: MATEDIST=mate_dist
				line_vcf = chrname2 + "\t" + to_string(mate_bnd_pos) + "\t" + id + "\t" + ref + "\t" + alt + "\t" + qual + "\t" + filter + "\t" + info + "\t" + format + "\t" + format_val_str;

				mate_bnd_item->reg_id = mate_reg_id;
				mate_bnd_item->clip_end = mate_clip_end;
				mate_bnd_item->mate_reg_id = reg_id;
				mate_bnd_item->mate_clip_end = clip_end;
				mate_bnd_item->mate_orient_ch = mate_orient_ch;
				mate_bnd_item->chrname = chrname2;
				mate_bnd_item->mate_chrname = chrname1;
				mate_bnd_item->bnd_pos = mate_bnd_pos;
				mate_bnd_item->mate_bnd_pos = bnd_pos;
				mate_bnd_item->seq = seq2_str;
				mate_bnd_item->bnd_str = mate_bnd_str2;
				mate_bnd_item->vcf_line = line_vcf;
				mate_bnd_item->support_num = support_num2;
				mate_bnd_item->mate_support_num = support_num;
				mate_bnd_item->cov_num = cov_num2;
				mate_bnd_item->mate_cov_num = cov_num;
				mate_bnd_item->mate_node = bnd_item;

				bnd_vec.push_back(bnd_item);
				bnd_vec.push_back(mate_bnd_item);
			}
		}
	}

	return bnd_vec;
}

// check BND String
void checkBNDStrVec(mateClipReg_t &mate_clip_reg){
	int64_t tra_pos_arr[4];
	int32_t i, checked_arr[4][2];
	string chrname1, chrname2;
	vector<string> bnd_str_vec;
	bool flag, sub_flag;

	flag = true;
	if(mate_clip_reg.valid_flag and mate_clip_reg.sv_type==VAR_TRA){ // TRA
		for(i=0; i<4; i++) { checked_arr[i][0] = checked_arr[i][1] = 0; } // initialize

		if(mate_clip_reg.leftClipReg) chrname1 = mate_clip_reg.leftClipReg->chrname;		// CHR1
		else if(mate_clip_reg.leftClipReg2) chrname1 = mate_clip_reg.leftClipReg2->chrname;		// CHR1
		else{
			//cout << __func__ << ", line=" << __LINE__ << ": null left part" << endl;
			flag = false;
			goto bnd_invalid_op;
		}
		if(mate_clip_reg.rightClipReg) chrname2 = mate_clip_reg.rightClipReg->chrname;		// CHR2
		else if(mate_clip_reg.rightClipReg2) chrname2 = mate_clip_reg.rightClipReg2->chrname;		// CHR2
		else{
			//cout << __func__ << ", line=" << __LINE__ << ": null right part" << endl;
			flag = false;
			goto bnd_invalid_op;
		}

		if((mate_clip_reg.leftClipReg==NULL and mate_clip_reg.bnd_mate_reg_strs[0].compare("-")!=0) or (mate_clip_reg.leftClipReg and mate_clip_reg.bnd_mate_reg_strs[0].compare("-")==0)
			or (mate_clip_reg.leftClipReg2==NULL and mate_clip_reg.bnd_mate_reg_strs[1].compare("-")!=0) or (mate_clip_reg.leftClipReg2 and mate_clip_reg.bnd_mate_reg_strs[1].compare("-")==0)
			or (mate_clip_reg.rightClipReg==NULL and mate_clip_reg.bnd_mate_reg_strs[2].compare("-")!=0) or (mate_clip_reg.rightClipReg and mate_clip_reg.bnd_mate_reg_strs[2].compare("-")==0)
			or (mate_clip_reg.rightClipReg2==NULL and mate_clip_reg.bnd_mate_reg_strs[3].compare("-")!=0) or (mate_clip_reg.rightClipReg2 and mate_clip_reg.bnd_mate_reg_strs[3].compare("-")==0))
		{
			flag = false;
			goto bnd_invalid_op;
		}

		if(mate_clip_reg.leftMeanClipPos!=0) tra_pos_arr[0] = mate_clip_reg.leftMeanClipPos;
		else tra_pos_arr[0] = -1;
		if(mate_clip_reg.leftMeanClipPos2!=0) tra_pos_arr[1] = mate_clip_reg.leftMeanClipPos2;
		else tra_pos_arr[1] = -1;
		if(mate_clip_reg.rightMeanClipPos!=0) tra_pos_arr[2] = mate_clip_reg.rightMeanClipPos;
		else tra_pos_arr[2] = -1;
		if(mate_clip_reg.rightMeanClipPos2!=0) tra_pos_arr[3] = mate_clip_reg.rightMeanClipPos2;
		else tra_pos_arr[3] = -1;

		for(i=0; i<4; i++) bnd_str_vec.push_back(mate_clip_reg.bnd_mate_reg_strs[i]);

		// four regions
		for(i=0; i<4; i++){
			if(mate_clip_reg.bnd_mate_reg_strs[i].compare("-")!=0){
				// right end
				sub_flag = isValidBNDStr(i, RIGHT_END, checked_arr, chrname1, chrname2, tra_pos_arr, bnd_str_vec);
				if(sub_flag==false){
					//cout << __func__ << ", line=" << __LINE__ << ": invalid BND item" << endl;
					flag = false;
					goto bnd_invalid_op;
				}

				// left end
				sub_flag = isValidBNDStr(i, LEFT_END, checked_arr, chrname1, chrname2, tra_pos_arr, bnd_str_vec);
				if(sub_flag==false){
					//cout << __func__ << ", line=" << __LINE__ << ": invalid BND item" << endl;
					flag = false;
					goto bnd_invalid_op;
				}
			}
		}
	}

	return;

bnd_invalid_op:
	if(flag==false)
		mate_clip_reg.valid_flag = false;
}

// check BND items
bool isValidBNDStr(int32_t reg_id, int32_t clip_end, int32_t checked_arr[][2], string &chrname1, string &chrname2, int64_t tra_pos_arr[4], vector<string> &bnd_str_vec){
	vector<string> bnd_str_vec2, bnd_str_vec2_tmp, mate_str_vec, mate_str_vec2;
	int64_t mate_bnd_pos;
	int32_t mate_reg_id, mate_reg_id2, mate_clip_end, sub_vec_id, mate_sub_vec_id;
	string mate_reg_str, mate_reg_str2, mate_bnd_str, seq_str;
	char mate_orient_ch, mate_orient_ch2;

	if(reg_id<0 or reg_id>=4) return false;

	if(tra_pos_arr[reg_id]!=-1 or bnd_str_vec.at(reg_id).compare("-")!=0) {
		sub_vec_id = (clip_end + 1) % 2;
		if(checked_arr[reg_id][sub_vec_id]==0){ // unprocessed, then process it

			bnd_str_vec2 = split(bnd_str_vec.at(reg_id), ",");

			mate_reg_str = bnd_str_vec2.at(sub_vec_id); // "2+|.|.|."
			if(mate_reg_str.compare("-")!=0){
				mate_str_vec = split(mate_reg_str, "|");
				mate_reg_id = stoi(mate_str_vec.at(0).substr(0, 1));
				mate_orient_ch = mate_str_vec.at(0).at(1);
//				support_num_str = mate_str_vec.at(2);
//				cov_num_str = mate_str_vec.at(3);
//				support_num = stoi(support_num_str);
//				cov_num = stoi(cov_num_str);

				mate_bnd_pos = tra_pos_arr[mate_reg_id];
				if(mate_bnd_pos==-1){
					//cout << chrname1 << ":" << bnd_pos << " <---> " << chrname2 << ":" << mate_bnd_pos << endl;
					return false;
				}

				checked_arr[reg_id][sub_vec_id] = 1;
				if(clip_end==RIGHT_END){ // right end
					if(mate_orient_ch=='+'){ // same orient
						mate_clip_end = LEFT_END;
						//reg_str = chrname2 + ":" + to_string(mate_bnd_pos) + "-" + to_string(mate_bnd_pos);
						//seq2 = fai_fetch(fai, reg_str.c_str(), &seq_len2);

						//seq2_str = seq2;
						mate_bnd_str = seq_str + "[" + chrname2 + ":" + to_string(mate_bnd_pos) + "[";
						//mate_bnd_str2 = "]" + chrname1 + ":" + to_string(bnd_pos) + "]" + seq2_str;
					}else{ // different orient
						mate_clip_end = RIGHT_END;
						mate_bnd_pos --;
						//reg_str = chrname2 + ":" + to_string(mate_bnd_pos) + "-" + to_string(mate_bnd_pos);
						//seq2 = fai_fetch(fai, reg_str.c_str(), &seq_len2);

						//seq2_str = seq2;
						mate_bnd_str = seq_str + "]" + chrname2 + ":" + to_string(mate_bnd_pos) + "]";
						//mate_bnd_str2 =  "[" + chrname1 + ":" + to_string(bnd_pos) + "[" + seq2_str;
					}
				}else{ // left end
					if(mate_orient_ch=='+'){ // same orient
						mate_clip_end = RIGHT_END;
						mate_bnd_pos --;
						//reg_str = chrname2 + ":" + to_string(mate_bnd_pos) + "-" + to_string(mate_bnd_pos);
						//seq2 = fai_fetch(fai, reg_str.c_str(), &seq_len2);

						//seq2_str = seq2;
						mate_bnd_str = "]" + chrname2 + ":" + to_string(mate_bnd_pos) + "]" + seq_str;
						//mate_bnd_str2 = seq2_str + "[" + chrname1 + ":" + to_string(bnd_pos) + "[";
					}else{ // different orient
						mate_clip_end = LEFT_END;
						//reg_str = chrname2 + ":" + to_string(mate_bnd_pos) + "-" + to_string(mate_bnd_pos);
						//seq2 = fai_fetch(fai, reg_str.c_str(), &seq_len2);

						//seq2_str = seq2;
						mate_bnd_str =  "[" + chrname2 + ":" + to_string(mate_bnd_pos) + "[" + seq_str;
						//mate_bnd_str2 = seq2_str + "[" + chrname1 + ":" + to_string(bnd_pos) + "[";
					}
				}
				//free(seq2);
				mate_sub_vec_id = (mate_clip_end + 1) % 2;
				if(checked_arr[mate_reg_id][mate_sub_vec_id]==0)
					checked_arr[mate_reg_id][mate_sub_vec_id] = 1;
				else{
					//cout << __func__ << ", line=" << __LINE__ << ": the BND information " << mate_bnd_str << " should not be processed twice, error." << endl;
					return false;
				}

				// mate coverage information
				if(bnd_str_vec.at(mate_reg_id).compare("-")!=0){
					bnd_str_vec2_tmp = split(bnd_str_vec.at(mate_reg_id), ",");
				}else{
					//cout << __func__ << ", line=" << __LINE__ << ": invalid BND information:" << bnd_str_vec.at(reg_id) << endl;
					return false;
				}

				mate_reg_str2 = bnd_str_vec2_tmp.at(mate_sub_vec_id); // "0+|..|..|.."
				if(mate_reg_str2.compare("-")!=0){
					mate_str_vec2 = split(mate_reg_str2, "|");
					mate_reg_id2 = stoi(mate_str_vec2.at(0).substr(0, 1));
					mate_orient_ch2 = mate_str_vec2.at(0).at(1);

					if(mate_reg_id2!=reg_id){
						//cout << __func__ << ", line=" << __LINE__ << ": two regions do not match, reg_id=" << reg_id << ", mate_reg_id2=" << mate_reg_id2 << ", error!" << endl;
						return false;
					}else if(mate_orient_ch2!=mate_orient_ch){
						//cout << __func__ << ", line=" << __LINE__ << ": orientation of two regions do not match, mate_orient=" << mate_orient_ch << ", mate_orient2=" << mate_orient_ch2 << ", error!" << endl;
						return false;
					}
				}else{
					//cout << __func__ << ", line=" << __LINE__ << ": invalid BND information, error." << endl;
					return false;
				}
			}
		}
	}

	return true;
}

void checkSuppNum(mateClipReg_t &mate_clip_reg, int32_t support_num_thres){
	vector<string> mate_reg_str_vec, mate_str_vec;
	size_t i, j;
	string support_num_str;
	int32_t support_num;

	if(mate_clip_reg.valid_flag){
		mate_clip_reg.supp_num_valid_flag = true;
		if(mate_clip_reg.large_indel_flag==false){
			for(i=0; i<4; i++){
				if(mate_clip_reg.bnd_mate_reg_strs[i].compare("-")!=0) {
					mate_reg_str_vec = split(mate_clip_reg.bnd_mate_reg_strs[i], ",");
					for(j=0; j<2; j++){
						if(mate_reg_str_vec.at(j).compare("-")!=0){
							mate_str_vec = split(mate_reg_str_vec.at(j), "|");
							support_num_str = mate_str_vec.at(2);
							//cov_num_str = mate_str_vec.at(3);
							support_num = stoi(support_num_str);
							//cov_num = stoi(cov_num_str);

							if(support_num<support_num_thres){
								mate_clip_reg.supp_num_valid_flag = false;
								return;
							}
						}
					}
				}
			}
		}else{
			if(mate_clip_reg.supp_num_largeIndel<support_num_thres)
				mate_clip_reg.supp_num_valid_flag = false;
		}
	}
}

// copy clipPos vector
void copyClipPosVec(vector<clipPos_t*> &sourceClipPosVector, vector<clipPos_t*> &destClipPosVector){
	size_t i;
	clipPos_t *clip_pos_item, *clip_pos_item_tmp;

	for(i=0; i<sourceClipPosVector.size(); i++){
		clip_pos_item = sourceClipPosVector.at(i);
		clip_pos_item_tmp = new clipPos_t();
		clip_pos_item_tmp->chrname = clip_pos_item->chrname;
		clip_pos_item_tmp->clipRefPos = clip_pos_item->clipRefPos;
		clip_pos_item_tmp->clipLocalRefPos = clip_pos_item->clipLocalRefPos;
		clip_pos_item_tmp->clipQueryPos = clip_pos_item->clipQueryPos;
		clip_pos_item_tmp->clip_aln = clip_pos_item->clip_aln;
		clip_pos_item_tmp->clip_end = clip_pos_item->clip_end;
		clip_pos_item_tmp->aln_orient = clip_pos_item->aln_orient;
		clip_pos_item_tmp->same_orient_flag = clip_pos_item->same_orient_flag;
		clip_pos_item_tmp->self_overlap_flag = clip_pos_item->self_overlap_flag;
		destClipPosVector.push_back(clip_pos_item_tmp);
	}
}

// compute the coverage of a clipping position
int32_t computeCovNumReg(string &chrname, int64_t startPos, int64_t endPos, faidx_t *fai, string &inBamFile){
	int64_t start_pos, end_pos, chr_len, pos, totalReadBeseNum, totalRefBaseNum;
	double mean_cov_num;
	Base *baseArray;
	vector<bam1_t*> alnDataVector;
	vector<bam1_t*>::iterator aln;

	start_pos = startPos - CLIP_END_EXTEND_SIZE / 2;
	end_pos = endPos + CLIP_END_EXTEND_SIZE / 2;

	chr_len = faidx_seq_len(fai, chrname.c_str()); // get the reference length
	if(start_pos<1) start_pos = 1;
	if(end_pos>chr_len) end_pos = chr_len;

	covLoader cov_loader(chrname, start_pos, end_pos, fai, 0, 0);
	baseArray = cov_loader.initBaseArray();
	alnDataLoader data_loader(chrname, start_pos, end_pos, inBamFile);
	data_loader.loadAlnData(alnDataVector);

	// generate the base coverage array
	cov_loader.generateBaseCoverage(baseArray, alnDataVector);

	totalReadBeseNum = totalRefBaseNum = 0;
	for(pos=start_pos; pos<=end_pos; pos++){
		// compute the meanCov excluding the gap regions
		if(baseArray[pos-start_pos].coverage.idx_RefBase!=4){ // excluding 'N'
			totalReadBeseNum += baseArray[pos-start_pos].coverage.num_bases[5];
			totalRefBaseNum ++;
//			num = baseArray[pos-start_pos].coverage.num_bases[5];
//			sum += num;
		}
	}
	if(totalRefBaseNum) mean_cov_num = (double)totalReadBeseNum/totalRefBaseNum;
	else mean_cov_num = 0;

	// release memory
	if(baseArray) { delete[] baseArray; baseArray = NULL; }
	if(!alnDataVector.empty()){
		for(aln=alnDataVector.begin(); aln!=alnDataVector.end(); aln++)
			bam_destroy1(*aln);
		vector<bam1_t*>().swap(alnDataVector);
	}

	return round(mean_cov_num);
}

// determine whether the variant is size satisfied
bool isSizeSatisfied(int64_t ref_dist, int64_t query_dist, int64_t min_sv_size_usr, int64_t max_sv_size_usr){
	int64_t dif;
	bool size_satisfy;

	dif = query_dist - ref_dist;
	if(dif<0) dif = -dif;

	size_satisfy = false;
	if((dif>=min_sv_size_usr and dif<=max_sv_size_usr) or (ref_dist>=min_sv_size_usr and ref_dist<=max_sv_size_usr) or (query_dist>=min_sv_size_usr and query_dist<=max_sv_size_usr))
		size_satisfy = true;

	return size_satisfy;
}

// determine whether the variant is size satisfied
bool isSizeSatisfied2(int64_t sv_len, int64_t min_sv_size_usr, int64_t max_sv_size_usr){
	int64_t size;
	bool size_satisfy;

	size = sv_len;

	size_satisfy = false;
	if(size<=max_sv_size_usr and size>=min_sv_size_usr)
		size_satisfy = true;

	return size_satisfy;
}

bool isNotAlreadyExists(vector<reg_t*> &varVec, reg_t *reg){
	bool not_existed = true;
	int32_t dis;
	size_t i;
	double dif_ratio;

//	for(i=0;i<varVec.size();i++){
//		if(reg->startRefPos==varVec[i]->startRefPos){
//			if(reg->var_type==varVec[i]->var_type){
//				if(reg->sv_len==varVec[i]->sv_len){
//					not_existed = false;
//				}
//			}
//		}
//	}

	for(i=0;i<varVec.size();i++){
		if(reg->sv_len<varVec[i]->sv_len){
			dif_ratio = (double)reg->sv_len/varVec[i]->sv_len;
		}else{
			dif_ratio = (double)varVec[i]->sv_len/reg->sv_len;
		}
		dis = abs(reg->startRefPos-varVec[i]->startRefPos);
		//dif_len = abs(reg->sv_len-varVec[i]->sv_len);
		if(reg->var_type==varVec[i]->var_type){
			if(dis<10 and dif_ratio>0.9){
				not_existed = false;
			}
		}
	}

	return not_existed;
}

bool isNotAlreadyExists(vector<reg_t*> &varVec, reg_t *reg, int32_t pos){
	bool not_existed = true;
	int32_t i, dis;
	double dif_ratio;

//	for(i=0;i<pos;i++){
//		if(reg->startRefPos==varVec[i]->startRefPos){
//			if(reg->var_type==varVec[i]->var_type){
//				if(reg->sv_len==varVec[i]->sv_len){
//					not_existed = false;
//				}
//			}
//		}
//	}
	//

	for(i=0; i<pos; i++){
		if(reg->sv_len<varVec[i]->sv_len){
			dif_ratio = (double)reg->sv_len/varVec[i]->sv_len;
		}else{
			dif_ratio = (double)varVec[i]->sv_len/reg->sv_len;
		}
		dis = abs(reg->startRefPos-varVec[i]->startRefPos);
		//dif_len = abs(reg->sv_len-varVec[i]->sv_len);
		if(reg->var_type==varVec[i]->var_type){
			if(dis<10 and dif_ratio>0.9){
				not_existed = false;
			}
		}
	}

	return not_existed;
}

// determine whether the chrome is decoy chrome
bool isDecoyChr(string &chrname){
	bool flag = false;
	size_t pos;

	pos = chrname.find(DECOY_PREFIX);
	if(pos==0) flag = true;
	else{
		pos = chrname.find(DECOY_PREFIX2);
		if(pos==0) flag = true;
	}

	return flag;
}

void removeVarCandNode(varCand *var_cand, vector<varCand*> &var_cand_vec){
	varCand *var_cand_tmp;
	for(size_t i=0; i<var_cand_vec.size(); i++){
		var_cand_tmp = var_cand_vec.at(i);
		if(var_cand_tmp==var_cand){
			var_cand_tmp->destroyVarCand();
			delete var_cand_tmp;
			var_cand_vec.erase(var_cand_vec.begin()+i);
			break;
		}
	}
}

// process single mate clipping region detection work
void *processSingleMateClipRegDetectWork(void *arg){
	mateClipRegDetectWork_opt *mate_clip_reg_work_opt = (mateClipRegDetectWork_opt *)arg;
	int32_t work_id = mate_clip_reg_work_opt->work_id;
	reg_t *reg = mate_clip_reg_work_opt->reg;
	faidx_t *fai = mate_clip_reg_work_opt->fai;
	Paras *paras = mate_clip_reg_work_opt->paras;
	//vector<bool> *clip_processed_flag_vec = mate_clip_reg_work_opt->clip_processed_flag_vec;
	vector<mateClipReg_t*> *mateClipRegVector = mate_clip_reg_work_opt->mateClipRegVector;
	vector<reg_t*> *clipRegVector = mate_clip_reg_work_opt->clipRegVector;
	vector<varCand*> *var_cand_clipReg_vec = mate_clip_reg_work_opt->var_cand_clipReg_vec;
	vector<Block*> *blockVector = mate_clip_reg_work_opt->blockVector;
	//int32_t *p_mate_clip_reg_fail_num = mate_clip_reg_work_opt->p_mate_clip_reg_fail_num;
	pthread_mutex_t *p_mutex_mate_clip_reg = mate_clip_reg_work_opt->p_mutex_mate_clip_reg;
	Time time;

	//cout << "\t[" << time.getTime() << "], [" << mate_clip_reg_work_opt->work_id << "]: " << reg->chrname << ":" << reg->startRefPos << "-" << reg->endRefPos << endl;

	clipReg clip_reg(reg->chrname, reg->startRefPos, reg->endRefPos, paras->inBamFile, fai, paras);
	clip_reg.computeMateClipReg();

	//cout << "\t[" << time.getTime() << "]: process clip regions " << reg->chrname << ":" << reg->startRefPos << "-" << reg->endRefPos << endl;
	//processClipRegs(work_id, *clip_processed_flag_vec, clip_reg.mate_clip_reg, reg, mateClipRegVector, clipRegVector, *var_cand_clipReg_vec, *blockVector, paras, p_mutex_mate_clip_reg, p_mate_clip_reg_fail_num);
	processClipRegs(work_id, clip_reg.mate_clip_reg, reg, mateClipRegVector, clipRegVector, *var_cand_clipReg_vec, blockVector, paras, p_mutex_mate_clip_reg);

	delete (mateClipRegDetectWork_opt *)arg;

	return NULL;
}

// process clip regions and mate clip regions
//void processClipRegs(int32_t work_id, vector<bool> &clip_processed_flag_vec, mateClipReg_t &mate_clip_reg, reg_t *clip_reg, vector<mateClipReg_t*> *mateClipRegVector, vector<reg_t*> *clipRegVector, vector<varCand*> &var_cand_clipReg_vec, vector<Block*> &blockVector, Paras *paras, pthread_mutex_t *p_mutex_mate_clip_reg, int32_t *mate_clip_reg_fail_num){
void processClipRegs(int32_t work_id, mateClipReg_t &mate_clip_reg, reg_t *clip_reg, vector<mateClipReg_t*> *mateClipRegVector, vector<reg_t*> *clipRegVector, vector<varCand*> &var_cand_clipReg_vec, vector<Block*> *blockVector, Paras *paras, pthread_mutex_t *p_mutex_mate_clip_reg){
	size_t i;
	reg_t *reg, *reg_tmp;
	mateClipReg_t *clip_reg_new;
	Block *bloc;
	int32_t idx_tmp;

	if(mate_clip_reg.valid_flag){
		// add mate clip region
		if(mate_clip_reg.leftClipReg or mate_clip_reg.leftClipReg2 or mate_clip_reg.rightClipReg or mate_clip_reg.rightClipReg2 or mate_clip_reg.large_indel_flag){
			clip_reg_new = new mateClipReg_t();
			clip_reg_new->work_id = work_id;
			clip_reg_new->leftClipReg = mate_clip_reg.leftClipReg;
			clip_reg_new->leftClipPosNum = mate_clip_reg.leftClipPosNum;
			clip_reg_new->rightClipReg = mate_clip_reg.rightClipReg;
			clip_reg_new->rightClipPosNum = mate_clip_reg.rightClipPosNum;
			clip_reg_new->leftMeanClipPos = mate_clip_reg.leftMeanClipPos;
			clip_reg_new->rightMeanClipPos = mate_clip_reg.rightMeanClipPos;

			clip_reg_new->leftClipReg2 = mate_clip_reg.leftClipReg2;
			clip_reg_new->rightClipReg2 = mate_clip_reg.rightClipReg2;
			clip_reg_new->leftClipPosNum2 = mate_clip_reg.leftClipPosNum2;
			clip_reg_new->rightClipPosNum2 = mate_clip_reg.rightClipPosNum2;
			clip_reg_new->leftClipRegNum = mate_clip_reg.leftClipRegNum;
			clip_reg_new->rightClipRegNum = mate_clip_reg.rightClipRegNum;
			clip_reg_new->leftMeanClipPos2 = mate_clip_reg.leftMeanClipPos2;
			clip_reg_new->rightMeanClipPos2 = mate_clip_reg.rightMeanClipPos2;

			clip_reg_new->reg_mated_flag = mate_clip_reg.reg_mated_flag;
			clip_reg_new->chrname_leftTra1 = clip_reg_new->chrname_rightTra1 = clip_reg_new->chrname_leftTra2 = clip_reg_new->chrname_rightTra2 = "";
			clip_reg_new->leftClipPosTra1 = clip_reg_new->rightClipPosTra1 = clip_reg_new->leftClipPosTra2 = clip_reg_new->rightClipPosTra2 = -1;
			clip_reg_new->sv_type = mate_clip_reg.sv_type;
			clip_reg_new->dup_num = mate_clip_reg.dup_num;
			clip_reg_new->valid_flag = mate_clip_reg.valid_flag;
			clip_reg_new->call_success_flag = false;
			clip_reg_new->tra_rescue_success_flag = false;
			clip_reg_new->supp_num_valid_flag = mate_clip_reg.supp_num_valid_flag;
			for(i=0; i<4; i++) clip_reg_new->bnd_mate_reg_strs[i] = mate_clip_reg.bnd_mate_reg_strs[i];

			clip_reg_new->large_indel_flag = mate_clip_reg.large_indel_flag;
			clip_reg_new->largeIndelClipReg = mate_clip_reg.largeIndelClipReg;
			clip_reg_new->supp_num_largeIndel = mate_clip_reg.supp_num_largeIndel;
			clip_reg_new->depth_largeIndel = mate_clip_reg.depth_largeIndel;

			pthread_mutex_lock(p_mutex_mate_clip_reg);
			mateClipRegVector->push_back(clip_reg_new);
			pthread_mutex_unlock(p_mutex_mate_clip_reg);
		}
	}else{
		// delete mate clip region
		if(mate_clip_reg.leftClipReg) { delete mate_clip_reg.leftClipReg; mate_clip_reg.leftClipReg = NULL; }
		if(mate_clip_reg.leftClipReg2) { delete mate_clip_reg.leftClipReg2; mate_clip_reg.leftClipReg2 = NULL; }
		if(mate_clip_reg.rightClipReg) { delete mate_clip_reg.rightClipReg; mate_clip_reg.rightClipReg = NULL; }
		if(mate_clip_reg.rightClipReg2) { delete mate_clip_reg.rightClipReg2; mate_clip_reg.rightClipReg2 = NULL; }

		pthread_mutex_lock(p_mutex_mate_clip_reg);
		if(mate_clip_reg.var_cand) { removeVarCandNode(mate_clip_reg.var_cand, var_cand_clipReg_vec); mate_clip_reg.var_cand = NULL; } // free item
		if(mate_clip_reg.left_var_cand_tra) { removeVarCandNode(mate_clip_reg.left_var_cand_tra, var_cand_clipReg_vec); mate_clip_reg.left_var_cand_tra = NULL; }  // free item
		if(mate_clip_reg.right_var_cand_tra) { removeVarCandNode(mate_clip_reg.right_var_cand_tra, var_cand_clipReg_vec); mate_clip_reg.right_var_cand_tra = NULL; }  // free item
		pthread_mutex_unlock(p_mutex_mate_clip_reg);

		if(mate_clip_reg.supp_num_valid_flag){
			// add the region into indel vector
			reg = new reg_t();
			reg->chrname = clip_reg->chrname;
			reg->startRefPos = clip_reg->startRefPos;
			reg->endRefPos = clip_reg->endRefPos;
			reg->startLocalRefPos = clip_reg->endLocalRefPos = 0;
			reg->startQueryPos = clip_reg->endQueryPos = 0;
			reg->var_type = VAR_UNC;
			reg->sv_len = 0;
			reg->query_id = -1;
			reg->blat_aln_id = -1;
			reg->minimap2_aln_id = -1;
			reg->call_success_status = false;
			reg->short_sv_flag = false;
			reg->zero_cov_flag = false;
			reg->aln_seg_end_flag = false;
			reg->query_pos_invalid_flag = false;
			reg->large_indel_flag = clip_reg->large_indel_flag;
			reg->gt_type = -1;
			reg->gt_seq = "";
			reg->AF = 0;

			// get position
			pthread_mutex_lock(p_mutex_mate_clip_reg);
			idx_tmp = -1;
			bloc = computeBlocByPos_util(clip_reg->startRefPos, *blockVector, paras);
			for(i=0; i<bloc->indelVector.size(); i++){
				reg_tmp = bloc->indelVector.at(i);
				if(reg->startRefPos<reg_tmp->startRefPos){
					idx_tmp = i;
					break;
				}
			}
			// add item
			if(idx_tmp!=-1) bloc->indelVector.insert(bloc->indelVector.begin()+idx_tmp, reg);
			else bloc->indelVector.push_back(reg);
			pthread_mutex_unlock(p_mutex_mate_clip_reg);
		}
	}
}

// compute block by reference position
Block* computeBlocByPos_util(int64_t begPos, vector<Block*> &block_vec, Paras *paras){
	int32_t bloc_ID = computeBlocID_util(begPos, block_vec, paras);
	return block_vec.at(bloc_ID);
}

// compute block ID
int32_t computeBlocID_util(int64_t begPos, vector<Block*> &block_vec, Paras *paras){
	int32_t bloc_ID, bloc_ID_tmp;
	size_t i;
	Block *bloc, *mid_bloc;

	bloc_ID_tmp = begPos / (paras->blockSize - 2 * paras->slideSize);
	if(bloc_ID_tmp>=(int32_t)block_vec.size()) bloc_ID_tmp = (int32_t)block_vec.size() - 1;

	bloc_ID = bloc_ID_tmp;

	mid_bloc = block_vec.at(bloc_ID_tmp);
	if(begPos<=mid_bloc->startPos){
		for(i=bloc_ID_tmp; i>=0; i--){
			bloc = block_vec.at(i);
			if(bloc->startPos<=begPos){
				bloc_ID = i;
				break;
			}
		}
	}else{
		for(i=bloc_ID_tmp+1; i<block_vec.size(); i++){
			bloc = block_vec.at(i);
			if(bloc->startPos>begPos){
				bloc_ID = i - 1;
				break;
			}
		}
	}

	return bloc_ID;
}

// sort the region vector
void sortRegVec(vector<reg_t*> &regVector){
	size_t min_idx;
	reg_t *tmp;

	if(regVector.size()>=2){
		for(size_t i=0; i<regVector.size(); i++){
			min_idx = i;
			for(size_t j=i+1; j<regVector.size(); j++)
				if(regVector[min_idx]->startRefPos>regVector[j]->startRefPos) min_idx = j;
			if(min_idx!=i){
				tmp = regVector[i];
				regVector[i] = regVector[min_idx];
				regVector[min_idx] = tmp;
				//cout << "i=" << i << ", min_idx=" << min_idx << ", swap!" << endl;
			}
		}
	}
}

bool isRegSorted(vector<reg_t*> &regVector){
	bool flag = true;
	if(regVector.size()>=2){
		for(size_t i=0; i<regVector.size()-1; i++){
			if(regVector[i]->startRefPos>regVector[i+1]->startRefPos){
				flag = false;
				break;
			}
		}
	}
	return flag;
}

// start work process monitor
void startWorkProcessMonitor(string &work_finish_filename, string &monitoring_proc_names, int32_t max_proc_running_minutes){
	if(isFileExist(work_finish_filename))
		remove(work_finish_filename.c_str());

	procMonitor_op *proc_monitor_op = new procMonitor_op();
	proc_monitor_op->work_finish_filename = work_finish_filename;
	proc_monitor_op->monitoring_proc_names = monitoring_proc_names;
	proc_monitor_op->max_proc_running_minutes = max_proc_running_minutes;

	// start new thread
	pthread_t tid;
	pthread_create(&tid, NULL, workProcessMonitor, proc_monitor_op);
	pthread_detach(tid);
}

// process monitor
void *workProcessMonitor(void *arg){
	procMonitor_op *proc_monitor_op = (procMonitor_op *)arg;
	FILE *fp;
	char buffer[256];
	string pscmd, line, tmp, work_finish_filename, monitoring_proc_names;
	int32_t uid, len, max_proc_running_minutes;
	vector<string> str_vec, tmp_vec, tmp_vec2;
	int32_t day, hour, minute, total;
	pid_t pid;

	work_finish_filename = proc_monitor_op->work_finish_filename;
	monitoring_proc_names = proc_monitor_op->monitoring_proc_names;
	max_proc_running_minutes = proc_monitor_op->max_proc_running_minutes;
	delete (procMonitor_op *)arg;

	while(1){
		if(isFileExist(work_finish_filename)){
			remove(work_finish_filename.c_str());
			//cout << "================ monitor file removed. ==============" << endl;
			break;
		}else{
			pscmd = "ps -C " + monitoring_proc_names;
			uid = getuid();
			pscmd += " -o comm,pid,time,uid | awk '$NF ~/" + to_string(uid) + "/'";
			fp = popen(pscmd.c_str(), "r");
			if(fp==NULL){
				//cout << __func__ << ": popen error, wait for " << MONITOR_WAIT_SECONDS << " seconds." << endl;
				sleep(MONITOR_WAIT_SECONDS);
				continue;
			}

			while(fgets(buffer, sizeof(buffer), fp)){
				//printf("%s", buffer);
				line = buffer;
				str_vec = split(line," ");
//				for(uint32_t i=0; i<str_vec.size(); i++){
//					cout << "["<< i << "]:" << str_vec.at(i) << endl;
//				}

				len = str_vec.size();
				if(len>=4){
					//compute cpu running time
					day = hour= minute = 0;
					if(str_vec.at(len-2).find('-')==string::npos){
						tmp_vec = split(str_vec.at(len-2), ":");
						hour = stoi(tmp_vec.at(0));
						minute = stoi(tmp_vec.at(1));
					}else{
						tmp_vec = split(str_vec.at(len-2), "-");
						day = stoi(tmp_vec.at(0));
						tmp_vec2 = split(tmp_vec.at(1), ":");
						hour = stoi(tmp_vec2.at(0));
						minute =  stoi(tmp_vec2.at(1));
					}
					total = day*24*60 + hour*60 + minute;
					//cout << total << endl;

					// kill process
					pid = stoi(str_vec.at(len-3));
					if(total>max_proc_running_minutes){
						kill(pid, SIGKILL);
						cout << __func__ << ": " << str_vec.at(0) << " (pid=" << pid << ") was killed because it exceeded the maximum running time (" << max_proc_running_minutes << " minutes). Please do NOT worry, this is programmed rather than error." << endl;
					}
				}
			}
			pclose(fp);
		}
		sleep(MONITOR_WAIT_SECONDS);
	}

	//cout << "-------------------- End of monitor. ----------------" << endl;

	return NULL;
}

// generate the alignment segments
vector<struct alnSeg*> generateAlnSegs(bam1_t* b){
	vector<struct alnSeg*> alnSegs;
	vector<struct MD_seg*> segs_MD;
	struct MD_seg* seg_MD;

	uint32_t *c, op, i = 0, j = 0, startRpos, startQpos;
	int32_t k, tmp_cigar_len, tmp_MD_len, common_len;  // common_len: used only for match flag
	string ins_str;

	startRpos = b->core.pos + 1;  // 1-based
	startQpos = 1;

	c = bam_get_cigar(b);  // CIGAR
	op = bam_cigar_op(c[i]);
	tmp_cigar_len = bam_cigar_oplen(c[i]);

	segs_MD = extractMDSegs(b);
	seg_MD = segs_MD[j];
	tmp_MD_len = seg_MD->seglen;
	while(i<b->core.n_cigar or j<segs_MD.size()){
		common_len = -1;
		switch(op){
			case BAM_CMATCH:
				if(tmp_cigar_len>=tmp_MD_len) // tmp_cigar_len >= seg_MD.len
					common_len = tmp_MD_len;
				else // tmp_cigar_len < seg_MD.len
					common_len = tmp_cigar_len;

				// add alnSeg item
				//if(common_len==1 and seg_MD->opflag==MD_MISMATCH){
				if(common_len==1 and seg_MD->opflag==BAM_CDIFF){
					// change the mismatched reference base to query base
					seg_MD->seg = "=ACMGRSVTWYHKDBN"[bam_seqi(bam_get_seq(b), startQpos-1)];
					//alnSegs.push_back(allocateAlnSeg(startRpos, startQpos, common_len, MD_MISMATCH, seg_MD->seg));
					alnSegs.push_back(allocateAlnSeg(startRpos, startQpos, common_len, BAM_CDIFF, seg_MD->seg));
				}else{
					//str_tmp = "";
					//alnSegs.push_back(allocateAlnSeg(startRpos, startQpos, common_len, BAM_CMATCH, str_tmp));
					alnSegs.push_back(allocateAlnSeg(startRpos, startQpos, common_len, BAM_CEQUAL, ""));
				}

				startRpos += common_len;
				startQpos += common_len;

				tmp_cigar_len -= common_len;
				tmp_MD_len -= common_len;

				if(tmp_cigar_len==0 and ++i<b->core.n_cigar){
					op = bam_cigar_op(c[i]);
					tmp_cigar_len = bam_cigar_oplen(c[i]);
				}
				if(tmp_MD_len==0 and ++j<segs_MD.size()){
					seg_MD = segs_MD[j];
					tmp_MD_len = seg_MD->seglen;
				}
				break;
			case BAM_CINS:
				// add alnSeg item
				ins_str = "";
				for(k=0; k<tmp_cigar_len; k++) ins_str += "=ACMGRSVTWYHKDBN"[bam_seqi(bam_get_seq(b), startQpos+k-1)];
				alnSegs.push_back(allocateAlnSeg(startRpos, startQpos, tmp_cigar_len, BAM_CINS, ins_str));
				startQpos += tmp_cigar_len;
				op = bam_cigar_op(c[++i]);
				tmp_cigar_len = bam_cigar_oplen(c[i]);

				break;
			case BAM_CDEL:
				// add alnSeg item
				if(tmp_cigar_len==tmp_MD_len){
					if(seg_MD->opflag==BAM_CDEL){
						alnSegs.push_back(allocateAlnSeg(startRpos, startQpos, tmp_cigar_len, BAM_CDEL, seg_MD->seg));
						startRpos += tmp_cigar_len;
						op = bam_cigar_op(c[++i]);
						tmp_cigar_len = bam_cigar_oplen(c[i]);
						seg_MD = segs_MD[++j];
						tmp_MD_len = seg_MD->seglen;
					}else{
						cerr << __func__ << ", line=" << __LINE__ << ": error" << endl;
						exit(1);
					}
				}else{
					cerr << __func__ << ", line=" << __LINE__ << ": error" << endl;
					exit(1);
				}
				break;
			case BAM_CSOFT_CLIP:
			case BAM_CHARD_CLIP:
				//str_tmp = to_string(tmp_cigar_len);
				alnSegs.push_back(allocateAlnSeg(startRpos, startQpos, tmp_cigar_len, op, to_string(tmp_cigar_len)));
				if(op==BAM_CSOFT_CLIP) startQpos += tmp_cigar_len;
				if(++i<b->core.n_cigar){
					op = bam_cigar_op(c[i]);
					tmp_cigar_len = bam_cigar_oplen(c[i]);
				}
				break;
			case BAM_CREF_SKIP:  // unexpected events
			case BAM_CPAD:
			case BAM_CEQUAL:
			case BAM_CDIFF:
			case BAM_CBACK:
			default:
				cerr << __func__ << ", line=" << __LINE__ << ": invalid opflag " << op << endl;
				exit(1);
		}
	}
	destroyMDSeg(segs_MD);

	return alnSegs;
}

// generate the alignment segments
vector<struct alnSeg*> generateAlnSegs2(bam1_t* b, int64_t startRefPos_paras, int64_t endRefPos_paras){
	vector<struct alnSeg*> alnSegs;
	vector<struct MD_seg*> segs_MD;
	struct MD_seg* seg_MD;

	uint32_t *c, op, i = 0, j = 0, startRpos, startQpos;
	int32_t k, tmp_cigar_len, tmp_MD_len, common_len, startmatch_idx, endmatch_idx;  // common_len: used only for match flag
	string ins_str;
	bool overlap_flag;
	int64_t end_ref_pos;


	startRpos = b->core.pos + 1;  // 1-based
	startQpos = 1;

	c = bam_get_cigar(b);  // CIGAR
	op = bam_cigar_op(c[i]);
	tmp_cigar_len = bam_cigar_oplen(c[i]);

	segs_MD = extractMDSegs(b);
	seg_MD = segs_MD[j];
	tmp_MD_len = seg_MD->seglen;
	while(i<b->core.n_cigar or j<segs_MD.size()){
		end_ref_pos = getEndRefPosAlnSeg(startRpos, op, tmp_cigar_len);
		overlap_flag = isOverlappedPos(startRpos, end_ref_pos, startRefPos_paras, endRefPos_paras);

		switch(op){
			case BAM_CMATCH:
				common_len = -1;
				if(tmp_cigar_len>=tmp_MD_len) // tmp_cigar_len >= seg_MD.len
					common_len = tmp_MD_len;
				else // tmp_cigar_len < seg_MD.len
					common_len = tmp_cigar_len;

				// recompute the overlap
				end_ref_pos = getEndRefPosAlnSeg(startRpos, op, common_len);
				overlap_flag = isOverlappedPos(startRpos, end_ref_pos, startRefPos_paras, endRefPos_paras);
				if(overlap_flag){
					// add alnSeg item
					//if(common_len==1 and seg_MD->opflag==MD_MISMATCH){
					if(common_len==1 and seg_MD->opflag==BAM_CDIFF){
						// change the mismatched reference base to query base
						seg_MD->seg = "=ACMGRSVTWYHKDBN"[bam_seqi(bam_get_seq(b), startQpos-1)];
						//alnSegs.push_back(allocateAlnSeg(startRpos, startQpos, common_len, MD_MISMATCH, seg_MD->seg));
						alnSegs.push_back(allocateAlnSeg(startRpos, startQpos, common_len, BAM_CDIFF, seg_MD->seg));
					}else{
						startmatch_idx = startRefPos_paras - startRpos;
						if(startmatch_idx<0) startmatch_idx = 0;
						endmatch_idx = common_len - 1;
						if(startRpos+common_len-1>endRefPos_paras)
							endmatch_idx -= (startRpos+common_len-1) - endRefPos_paras;
						//str_tmp = "";
						//alnSegs.push_back(allocateAlnSeg(startRpos, startQpos, common_len, BAM_CMATCH, str_tmp));
						alnSegs.push_back(allocateAlnSeg(startRpos+startmatch_idx, startQpos+startmatch_idx, endmatch_idx-startmatch_idx+1, BAM_CEQUAL, ""));
					}
				}

				startRpos += common_len;
				startQpos += common_len;

				tmp_cigar_len -= common_len;
				tmp_MD_len -= common_len;

				if(tmp_cigar_len==0 and ++i<b->core.n_cigar){
					op = bam_cigar_op(c[i]);
					tmp_cigar_len = bam_cigar_oplen(c[i]);
				}
				if(tmp_MD_len==0 and ++j<segs_MD.size()){
					seg_MD = segs_MD[j];
					tmp_MD_len = seg_MD->seglen;
				}
				break;
			case BAM_CINS:
				if(overlap_flag){
					// add alnSeg item
					ins_str = "";
					for(k=0; k<tmp_cigar_len; k++) ins_str += "=ACMGRSVTWYHKDBN"[bam_seqi(bam_get_seq(b), startQpos+k-1)];
					alnSegs.push_back(allocateAlnSeg(startRpos, startQpos, tmp_cigar_len, BAM_CINS, ins_str));
				}
				startQpos += tmp_cigar_len;
				op = bam_cigar_op(c[++i]);
				tmp_cigar_len = bam_cigar_oplen(c[i]);
				break;
			case BAM_CDEL:
				// add alnSeg item
				if(tmp_cigar_len==tmp_MD_len){
					if(seg_MD->opflag==BAM_CDEL){
						if(overlap_flag)
							alnSegs.push_back(allocateAlnSeg(startRpos, startQpos, tmp_cigar_len, BAM_CDEL, seg_MD->seg));
						startRpos += tmp_cigar_len;
						op = bam_cigar_op(c[++i]);
						tmp_cigar_len = bam_cigar_oplen(c[i]);
						seg_MD = segs_MD[++j];
						tmp_MD_len = seg_MD->seglen;
					}else{
						cerr << __func__ << ", line=" << __LINE__ << ": error" << endl;
						exit(1);
					}
				}else{
					cerr << __func__ << ", line=" << __LINE__ << ": error" << endl;
					exit(1);
				}
				break;
			case BAM_CSOFT_CLIP:
			case BAM_CHARD_CLIP:
				if(overlap_flag){
					//str_tmp = to_string(tmp_cigar_len);
					alnSegs.push_back(allocateAlnSeg(startRpos, startQpos, tmp_cigar_len, op, to_string(tmp_cigar_len)));
				}
				if(op==BAM_CSOFT_CLIP) startQpos += tmp_cigar_len;
				if(++i<b->core.n_cigar){
					op = bam_cigar_op(c[i]);
					tmp_cigar_len = bam_cigar_oplen(c[i]);
				}
				break;
			case BAM_CREF_SKIP:  // unexpected events
			case BAM_CPAD:
			case BAM_CEQUAL:
			case BAM_CDIFF:
			case BAM_CBACK:
			default:
				cerr << __func__ << ", line=" << __LINE__ << ": invalid opflag " << op << endl;
				exit(1);
		}
	}
	destroyMDSeg(segs_MD);

	return alnSegs;
}

vector<struct alnSeg*> generateAlnSegs_no_MD(bam1_t* b, Base* BaseArr, int64_t startRefPos_paras, int64_t endRefPos_paras){
	vector<struct alnSeg*> alnSegs;
	uint32_t *c, op, i = 0,startRpos, startQpos;
	int32_t k, tmp_cigar_len, mis_idx, startmatch_idx, endmatch_idx;
	char queBase, ref;
	string ins_str, diff, del;

	startRpos = b->core.pos + 1;  // 1-based
	startQpos = 1;

	c = bam_get_cigar(b);  // CIGAR
	op = bam_cigar_op(c[i]);
	tmp_cigar_len = bam_cigar_oplen(c[i]);
	//string qname = bam_get_qname(b);

	while(i<b->core.n_cigar){
		switch(op){
			case BAM_CMATCH:
				mis_idx = -1;
				for(k=0; k<tmp_cigar_len; k++){
					if(startRpos+k>=startRefPos_paras and startRpos+k<=endRefPos_paras){
						queBase = "=ACMGRSVTWYHKDBN"[bam_seqi(bam_get_seq(b), startQpos+k-1)];
						ref = BaseArr[startRpos+k-startRefPos_paras].coverage.refBase;
						if(isBaseMatch(queBase,ref)==false){
							if(mis_idx==-1) {
								startmatch_idx = 0;
								endmatch_idx = k;//
								//str_tmp = "";
								//alnSegs.push_back(allocateAlnSeg(startRpos, startQpos, endmatch_idx, BAM_CMATCH, str_tmp));
								alnSegs.push_back(allocateAlnSeg(startRpos, startQpos, endmatch_idx, BAM_CEQUAL, ""));
							}else {
								startmatch_idx = mis_idx + 1;
								endmatch_idx = k - 1;
								//str_tmp = "";
								//alnSegs.push_back(allocateAlnSeg(startRpos+startmatch_idx+1, startQpos+startmatch_idx, endmatch_idx-startmatch_idx, BAM_CMATCH, str_tmp));
								alnSegs.push_back(allocateAlnSeg(startRpos+startmatch_idx, startQpos+startmatch_idx, endmatch_idx-startmatch_idx+1, BAM_CEQUAL, ""));
							}
							//diff = "=ACMGRSVTWYHKDBN"[bam_seqi(bam_get_seq(b), startQpos+k-1)];
							diff = queBase;
							//alnSegs.push_back(allocateAlnSeg(startRpos+k, startQpos+k, 1, MD_MISMATCH, diff));
							alnSegs.push_back(allocateAlnSeg(startRpos+k, startQpos+k, 1, BAM_CDIFF, diff));
							mis_idx = k;
						}
					}
				}
				if(mis_idx==-1){
					//str_tmp = "";
					//alnSegs.push_back(allocateAlnSeg(startRpos, startQpos, tmp_cigar_len, BAM_CMATCH, str_tmp));
					alnSegs.push_back(allocateAlnSeg(startRpos, startQpos, tmp_cigar_len, BAM_CEQUAL, ""));
				}else {
					startmatch_idx = mis_idx + 1;
					//str_tmp = "";
					//alnSegs.push_back(allocateAlnSeg(startRpos+startmatch_idx+1, startQpos+startmatch_idx+1, tmp_cigar_len-startmatch_idx-1, BAM_CMATCH, str_tmp));
					alnSegs.push_back(allocateAlnSeg(startRpos+startmatch_idx, startQpos+startmatch_idx, tmp_cigar_len-startmatch_idx, BAM_CEQUAL, ""));
				}
				startQpos += tmp_cigar_len;
				startRpos += tmp_cigar_len;
				op = bam_cigar_op(c[++i]);
				tmp_cigar_len = bam_cigar_oplen(c[i]);
				break;
			case BAM_CINS:
				// add alnSeg item
				ins_str = "";
				for(k=0; k<tmp_cigar_len; k++) ins_str += "=ACMGRSVTWYHKDBN"[bam_seqi(bam_get_seq(b), startQpos+k-1)];
				alnSegs.push_back(allocateAlnSeg(startRpos, startQpos, tmp_cigar_len, BAM_CINS, ins_str));
				startQpos += tmp_cigar_len;
				op = bam_cigar_op(c[++i]);
				tmp_cigar_len = bam_cigar_oplen(c[i]);

				break;
			case BAM_CDEL:
				// add alnSeg item
				del = "";
				for(k=0; k<tmp_cigar_len; k++) del += "=ACMGRSVTWYHKDBN"[bam_seqi(bam_get_seq(b), startQpos+k-1)];//something wrong?
				alnSegs.push_back(allocateAlnSeg(startRpos, startQpos, tmp_cigar_len, BAM_CDEL, del));
				startRpos += tmp_cigar_len;
				op = bam_cigar_op(c[++i]);
				tmp_cigar_len = bam_cigar_oplen(c[i]);
				break;
			case BAM_CSOFT_CLIP:
			case BAM_CHARD_CLIP:
				//str_tmp = to_string(tmp_cigar_len);
				alnSegs.push_back(allocateAlnSeg(startRpos, startQpos, tmp_cigar_len, op, to_string(tmp_cigar_len)));
				if(op==BAM_CSOFT_CLIP) startQpos += tmp_cigar_len;
				if(++i<b->core.n_cigar){
					op = bam_cigar_op(c[i]);
					tmp_cigar_len = bam_cigar_oplen(c[i]);
				}
				break;
			case BAM_CEQUAL:
				//str_tmp = "";
				alnSegs.push_back(allocateAlnSeg(startRpos, startQpos, tmp_cigar_len, BAM_CEQUAL, ""));
				startQpos += tmp_cigar_len;
				startRpos += tmp_cigar_len;
				op = bam_cigar_op(c[++i]);
				tmp_cigar_len = bam_cigar_oplen(c[i]);
				break;
			case BAM_CDIFF:
				diff = "";
				for(k=0; k<tmp_cigar_len; k++) diff += "=ACMGRSVTWYHKDBN"[bam_seqi(bam_get_seq(b), startQpos+k-1)];
				alnSegs.push_back(allocateAlnSeg(startRpos, startQpos, tmp_cigar_len, BAM_CDIFF, diff));

				startQpos += tmp_cigar_len;
				startRpos += tmp_cigar_len;
				op = bam_cigar_op(c[++i]);
				tmp_cigar_len = bam_cigar_oplen(c[i]);
				break;
			case BAM_CREF_SKIP: // unexpected events
			case BAM_CPAD:
			case BAM_CBACK:
			default:
				cerr << __func__ << ", line=" << __LINE__ << ": invalid opflag " << op << endl;
				exit(1);
		}
	}
	return alnSegs;
}

vector<struct alnSeg*> generateAlnSegs_no_MD2(bam1_t* b, Base* BaseArr, int64_t startRefPos_paras, int64_t endRefPos_paras){
	vector<struct alnSeg*> alnSegs;
	uint32_t *c, op, i = 0,startRpos, startQpos;
	int32_t k, tmp_cigar_len, mis_idx, startmatch_idx, endmatch_idx;
	char queBase, ref;
	string ins_str, diff, del;
	bool overlap_flag;
	int64_t end_ref_pos;

	startRpos = b->core.pos + 1;  // 1-based
	startQpos = 1;

	c = bam_get_cigar(b);  // CIGAR
	op = bam_cigar_op(c[i]);
	tmp_cigar_len = bam_cigar_oplen(c[i]);
	//string qname = bam_get_qname(b);

	while(i<b->core.n_cigar){
		end_ref_pos = getEndRefPosAlnSeg(startRpos, op, tmp_cigar_len);
		overlap_flag = isOverlappedPos(startRpos, end_ref_pos, startRefPos_paras, endRefPos_paras);
		switch(op){
			case BAM_CMATCH:
				if(overlap_flag){
					mis_idx = -1;
					for(k=0; k<tmp_cigar_len; k++){
						if(startRpos+k>=startRefPos_paras and startRpos+k<=endRefPos_paras){
							queBase = "=ACMGRSVTWYHKDBN"[bam_seqi(bam_get_seq(b), startQpos+k-1)];
							ref = BaseArr[startRpos+k-startRefPos_paras].coverage.refBase;
							if(isBaseMatch(queBase,ref)==false){
								if(startRpos+k>startRefPos_paras){
									if(mis_idx==-1) {
										//startmatch_idx = 0;
										startmatch_idx = startRefPos_paras - startRpos;
										if(startmatch_idx<0) startmatch_idx = 0;
									}else
										startmatch_idx = mis_idx + 1;
									endmatch_idx = k - 1;
									//str_tmp = "";
									//alnSegs.push_back(allocateAlnSeg(startRpos+startmatch_idx+1, startQpos+startmatch_idx, endmatch_idx-startmatch_idx, BAM_CMATCH, str_tmp));
									alnSegs.push_back(allocateAlnSeg(startRpos+startmatch_idx, startQpos+startmatch_idx, endmatch_idx-startmatch_idx+1, BAM_CEQUAL, ""));
								}
								//diff = "=ACMGRSVTWYHKDBN"[bam_seqi(bam_get_seq(b), startQpos+k-1)];
								diff = queBase;
								//alnSegs.push_back(allocateAlnSeg(startRpos+k, startQpos+k, 1, MD_MISMATCH, diff));
								alnSegs.push_back(allocateAlnSeg(startRpos+k, startQpos+k, 1, BAM_CDIFF, diff));
								mis_idx = k;
							}
						}
					}
					if(mis_idx==-1){
						startmatch_idx = startRefPos_paras - startRpos;
						if(startmatch_idx<0) startmatch_idx = 0;
						endmatch_idx = tmp_cigar_len - 1;
						if(startRpos+tmp_cigar_len-1>endRefPos_paras)
							endmatch_idx -= (startRpos+tmp_cigar_len-1) - endRefPos_paras;
						//str_tmp = "";
						//alnSegs.push_back(allocateAlnSeg(startRpos, startQpos, tmp_cigar_len, BAM_CMATCH, str_tmp));
						alnSegs.push_back(allocateAlnSeg(startRpos+startmatch_idx, startQpos+startmatch_idx, endmatch_idx-startmatch_idx+1, BAM_CEQUAL, ""));
					}else {
						startmatch_idx = mis_idx + 1;
						endmatch_idx = tmp_cigar_len - 1;
						if(startRpos+tmp_cigar_len-1>endRefPos_paras)
							endmatch_idx -= (startRpos+tmp_cigar_len-1) - endRefPos_paras;
						//str_tmp = "";
						//alnSegs.push_back(allocateAlnSeg(startRpos+startmatch_idx+1, startQpos+startmatch_idx+1, tmp_cigar_len-startmatch_idx-1, BAM_CMATCH, str_tmp));
						//alnSegs.push_back(allocateAlnSeg(startRpos+startmatch_idx, startQpos+startmatch_idx, tmp_cigar_len-startmatch_idx, BAM_CEQUAL, ""));
						alnSegs.push_back(allocateAlnSeg(startRpos+startmatch_idx, startQpos+startmatch_idx, endmatch_idx-startmatch_idx+1, BAM_CEQUAL, ""));
					}
				}
				startQpos += tmp_cigar_len;
				startRpos += tmp_cigar_len;
				op = bam_cigar_op(c[++i]);
				tmp_cigar_len = bam_cigar_oplen(c[i]);
				break;
			case BAM_CINS:
				if(overlap_flag){
					// add alnSeg item
					ins_str = "";
					for(k=0; k<tmp_cigar_len; k++) ins_str += "=ACMGRSVTWYHKDBN"[bam_seqi(bam_get_seq(b), startQpos+k-1)];
					alnSegs.push_back(allocateAlnSeg(startRpos, startQpos, tmp_cigar_len, BAM_CINS, ins_str));
				}
				startQpos += tmp_cigar_len;
				op = bam_cigar_op(c[++i]);
				tmp_cigar_len = bam_cigar_oplen(c[i]);

				break;
			case BAM_CDEL:
				if(overlap_flag){
					// add alnSeg item
					del = "";
					for(k=0; k<tmp_cigar_len; k++) del += "=ACMGRSVTWYHKDBN"[bam_seqi(bam_get_seq(b), startQpos+k-1)];//something wrong?
					alnSegs.push_back(allocateAlnSeg(startRpos, startQpos, tmp_cigar_len, BAM_CDEL, del));
				}
				startRpos += tmp_cigar_len;
				op = bam_cigar_op(c[++i]);
				tmp_cigar_len = bam_cigar_oplen(c[i]);
				break;
			case BAM_CSOFT_CLIP:
			case BAM_CHARD_CLIP:
				if(overlap_flag){
					//str_tmp = to_string(tmp_cigar_len);
					alnSegs.push_back(allocateAlnSeg(startRpos, startQpos, tmp_cigar_len, op, to_string(tmp_cigar_len)));
				}
				if(op==BAM_CSOFT_CLIP) startQpos += tmp_cigar_len;
				if(++i<b->core.n_cigar){
					op = bam_cigar_op(c[i]);
					tmp_cigar_len = bam_cigar_oplen(c[i]);
				}
				break;
			case BAM_CEQUAL:
				if(overlap_flag){
					//str_tmp = "";
					alnSegs.push_back(allocateAlnSeg(startRpos, startQpos, tmp_cigar_len, BAM_CEQUAL, ""));
				}
				startQpos += tmp_cigar_len;
				startRpos += tmp_cigar_len;
				op = bam_cigar_op(c[++i]);
				tmp_cigar_len = bam_cigar_oplen(c[i]);
				break;
			case BAM_CDIFF:
				if(overlap_flag){
					diff = "";
					for(k=0; k<tmp_cigar_len; k++) diff += "=ACMGRSVTWYHKDBN"[bam_seqi(bam_get_seq(b), startQpos+k-1)];
					alnSegs.push_back(allocateAlnSeg(startRpos, startQpos, tmp_cigar_len, BAM_CDIFF, diff));
				}
				startQpos += tmp_cigar_len;
				startRpos += tmp_cigar_len;
				op = bam_cigar_op(c[++i]);
				tmp_cigar_len = bam_cigar_oplen(c[i]);
				break;
			case BAM_CREF_SKIP: // unexpected events
			case BAM_CPAD:
			case BAM_CBACK:
			default:
				cerr << __func__ << ", line=" << __LINE__ << ": invalid opflag " << op << endl;
				exit(1);
		}
	}
	return alnSegs;
}

vector<struct alnSeg*> generateAlnSegs_no_MD(bam1_t* b, string &refseq, int64_t startRefPos_paras, int64_t endRefPos_paras){
	vector<struct alnSeg*> alnSegs;
	uint32_t *c, op, i = 0,startRpos, startQpos;
	int32_t k, tmp_cigar_len, mis_idx, startmatch_idx, endmatch_idx;
	char queBase, ref;
	string ins_str, diff, del;

	startRpos = b->core.pos + 1;  // 1-based
	startQpos = 1;

	c = bam_get_cigar(b);  // CIGAR
	op = bam_cigar_op(c[i]);
	tmp_cigar_len = bam_cigar_oplen(c[i]);
	//string qname = bam_get_qname(b);

	while(i<b->core.n_cigar){
		switch(op){
			case BAM_CMATCH:
				mis_idx = -1;
				for(k=0; k<tmp_cigar_len; k++){
					if(startRpos+k>=startRefPos_paras and startRpos+k<=endRefPos_paras){
						queBase = "=ACMGRSVTWYHKDBN"[bam_seqi(bam_get_seq(b), startQpos+k-1)];
						ref = refseq.at(startRpos+k-startRefPos_paras);
						if(isBaseMatch(queBase,ref)==false){
							if(mis_idx==-1) {
								startmatch_idx = 0;
								endmatch_idx = k;//
								//str_tmp = "";
								//alnSegs.push_back(allocateAlnSeg(startRpos, startQpos, endmatch_idx, BAM_CMATCH, str_tmp));
								alnSegs.push_back(allocateAlnSeg(startRpos, startQpos, endmatch_idx, BAM_CEQUAL, ""));
							}else {
								startmatch_idx = mis_idx + 1;
								endmatch_idx = k - 1;
								//str_tmp = "";
								//alnSegs.push_back(allocateAlnSeg(startRpos+startmatch_idx+1, startQpos+startmatch_idx, endmatch_idx-startmatch_idx, BAM_CMATCH, str_tmp));
								alnSegs.push_back(allocateAlnSeg(startRpos+startmatch_idx, startQpos+startmatch_idx, endmatch_idx-startmatch_idx+1, BAM_CEQUAL, ""));
							}
							//diff = "=ACMGRSVTWYHKDBN"[bam_seqi(bam_get_seq(b), startQpos+k-1)];
							diff = queBase;
							//alnSegs.push_back(allocateAlnSeg(startRpos+k, startQpos+k, 1, MD_MISMATCH, diff));
							alnSegs.push_back(allocateAlnSeg(startRpos+k, startQpos+k, 1, BAM_CDIFF, diff));
							mis_idx = k;
						}
					}
				}
				if(mis_idx==-1){
					//str_tmp = "";
					//alnSegs.push_back(allocateAlnSeg(startRpos, startQpos, tmp_cigar_len, BAM_CMATCH, str_tmp));
					alnSegs.push_back(allocateAlnSeg(startRpos, startQpos, tmp_cigar_len, BAM_CEQUAL, ""));
				}else {
					startmatch_idx = mis_idx + 1;
					//str_tmp = "";
					//alnSegs.push_back(allocateAlnSeg(startRpos+startmatch_idx+1, startQpos+startmatch_idx+1, tmp_cigar_len-startmatch_idx-1, BAM_CMATCH, str_tmp));
					alnSegs.push_back(allocateAlnSeg(startRpos+startmatch_idx, startQpos+startmatch_idx, tmp_cigar_len-startmatch_idx, BAM_CEQUAL, ""));
				}
				startQpos += tmp_cigar_len;
				startRpos += tmp_cigar_len;
				op = bam_cigar_op(c[++i]);
				tmp_cigar_len = bam_cigar_oplen(c[i]);
				break;
			case BAM_CINS:
				// add alnSeg item
				ins_str = "";
				for(k=0; k<tmp_cigar_len; k++) ins_str += "=ACMGRSVTWYHKDBN"[bam_seqi(bam_get_seq(b), startQpos+k-1)];
				alnSegs.push_back(allocateAlnSeg(startRpos, startQpos, tmp_cigar_len, BAM_CINS, ins_str));
				startQpos += tmp_cigar_len;
				op = bam_cigar_op(c[++i]);
				tmp_cigar_len = bam_cigar_oplen(c[i]);

				break;
			case BAM_CDEL:
				// add alnSeg item
				del = "";
				for(k=0; k<tmp_cigar_len; k++) del += "=ACMGRSVTWYHKDBN"[bam_seqi(bam_get_seq(b), startQpos+k-1)];//something wrong?
				alnSegs.push_back(allocateAlnSeg(startRpos, startQpos, tmp_cigar_len, BAM_CDEL, del));
				startRpos += tmp_cigar_len;
				op = bam_cigar_op(c[++i]);
				tmp_cigar_len = bam_cigar_oplen(c[i]);
				break;
			case BAM_CSOFT_CLIP:
			case BAM_CHARD_CLIP:
				//str_tmp = to_string(tmp_cigar_len);
				alnSegs.push_back(allocateAlnSeg(startRpos, startQpos, tmp_cigar_len, op, to_string(tmp_cigar_len)));
				if(op==BAM_CSOFT_CLIP) startQpos += tmp_cigar_len;
				if(++i<b->core.n_cigar){
					op = bam_cigar_op(c[i]);
					tmp_cigar_len = bam_cigar_oplen(c[i]);
				}
				break;
			case BAM_CEQUAL:
				//str_tmp = "";
				alnSegs.push_back(allocateAlnSeg(startRpos, startQpos, tmp_cigar_len, BAM_CEQUAL, ""));
				startQpos += tmp_cigar_len;
				startRpos += tmp_cigar_len;
				op = bam_cigar_op(c[++i]);
				tmp_cigar_len = bam_cigar_oplen(c[i]);
				break;
			case BAM_CDIFF:
				diff = "";
				for(k=0; k<tmp_cigar_len; k++) diff += "=ACMGRSVTWYHKDBN"[bam_seqi(bam_get_seq(b), startQpos+k-1)];
				alnSegs.push_back(allocateAlnSeg(startRpos, startQpos, tmp_cigar_len, BAM_CDIFF, diff));

				startQpos += tmp_cigar_len;
				startRpos += tmp_cigar_len;
				op = bam_cigar_op(c[++i]);
				tmp_cigar_len = bam_cigar_oplen(c[i]);
				break;
			case BAM_CREF_SKIP: // unexpected events
			case BAM_CPAD:
			case BAM_CBACK:
			default:
				cerr << __func__ << ", line=" << __LINE__ << ": invalid opflag " << op << endl;
				exit(1);
		}
	}
	return alnSegs;
}

vector<struct alnSeg*> generateAlnSegs_no_MD2(bam1_t* b, string &refseq, int64_t startRefPos_paras, int64_t endRefPos_paras){
	vector<struct alnSeg*> alnSegs;
	uint32_t *c, op, i = 0,startRpos, startQpos;
	int32_t k, tmp_cigar_len, mis_idx, startmatch_idx, endmatch_idx;
	char queBase, ref;
	string ins_str, diff, del;
	bool overlap_flag;
	int64_t end_ref_pos;

	startRpos = b->core.pos + 1;  // 1-based
	startQpos = 1;

	c = bam_get_cigar(b);  // CIGAR
	op = bam_cigar_op(c[i]);
	tmp_cigar_len = bam_cigar_oplen(c[i]);
	//string qname = bam_get_qname(b);

	while(i<b->core.n_cigar){
		end_ref_pos = getEndRefPosAlnSeg(startRpos, op, tmp_cigar_len);
		overlap_flag = isOverlappedPos(startRpos, end_ref_pos, startRefPos_paras, endRefPos_paras);
		switch(op){
			case BAM_CMATCH:
				if(overlap_flag){
					mis_idx = -1;
					for(k=0; k<tmp_cigar_len; k++){
						if(startRpos+k>=startRefPos_paras and startRpos+k<=endRefPos_paras){
							queBase = "=ACMGRSVTWYHKDBN"[bam_seqi(bam_get_seq(b), startQpos+k-1)];
							ref = refseq.at(startRpos+k-startRefPos_paras);
							if(isBaseMatch(queBase,ref)==false){
								if(startRpos+k>startRefPos_paras){
									if(mis_idx==-1) {
										//startmatch_idx = 0;
										startmatch_idx = startRefPos_paras - startRpos;
										if(startmatch_idx<0) startmatch_idx = 0;
									}else
										startmatch_idx = mis_idx + 1;
									endmatch_idx = k - 1;
									//str_tmp = "";
									//alnSegs.push_back(allocateAlnSeg(startRpos+startmatch_idx+1, startQpos+startmatch_idx, endmatch_idx-startmatch_idx, BAM_CMATCH, str_tmp));
									alnSegs.push_back(allocateAlnSeg(startRpos+startmatch_idx, startQpos+startmatch_idx, endmatch_idx-startmatch_idx+1, BAM_CEQUAL, ""));
								}
								//diff = "=ACMGRSVTWYHKDBN"[bam_seqi(bam_get_seq(b), startQpos+k-1)];
								diff = queBase;
								//alnSegs.push_back(allocateAlnSeg(startRpos+k, startQpos+k, 1, MD_MISMATCH, diff));
								alnSegs.push_back(allocateAlnSeg(startRpos+k, startQpos+k, 1, BAM_CDIFF, diff));
								mis_idx = k;
							}
						}
					}
					if(mis_idx==-1){
						startmatch_idx = startRefPos_paras - startRpos;
						if(startmatch_idx<0) startmatch_idx = 0;
						endmatch_idx = tmp_cigar_len - 1;
						if(startRpos+tmp_cigar_len-1>endRefPos_paras)
							endmatch_idx -= (startRpos+tmp_cigar_len-1) - endRefPos_paras;
						//str_tmp = "";
						//alnSegs.push_back(allocateAlnSeg(startRpos, startQpos, tmp_cigar_len, BAM_CMATCH, str_tmp));
						//alnSegs.push_back(allocateAlnSeg(startRpos+startmatch_idx, startQpos+startmatch_idx, tmp_cigar_len-startmatch_idx, BAM_CEQUAL, ""));
						alnSegs.push_back(allocateAlnSeg(startRpos+startmatch_idx, startQpos+startmatch_idx, endmatch_idx-startmatch_idx+1, BAM_CEQUAL, ""));
					}else {
						startmatch_idx = mis_idx + 1;
						endmatch_idx = tmp_cigar_len - 1;
						if(startRpos+tmp_cigar_len-1>endRefPos_paras)
							endmatch_idx -= (startRpos+tmp_cigar_len-1) - endRefPos_paras;
						//str_tmp = "";
						//alnSegs.push_back(allocateAlnSeg(startRpos+startmatch_idx+1, startQpos+startmatch_idx+1, tmp_cigar_len-startmatch_idx-1, BAM_CMATCH, str_tmp));
						//alnSegs.push_back(allocateAlnSeg(startRpos+startmatch_idx, startQpos+startmatch_idx, tmp_cigar_len-startmatch_idx, BAM_CEQUAL, ""));
						alnSegs.push_back(allocateAlnSeg(startRpos+startmatch_idx, startQpos+startmatch_idx, endmatch_idx-startmatch_idx+1, BAM_CEQUAL, ""));
					}
				}

				startQpos += tmp_cigar_len;
				startRpos += tmp_cigar_len;
				op = bam_cigar_op(c[++i]);
				tmp_cigar_len = bam_cigar_oplen(c[i]);
				break;
			case BAM_CINS:
				if(overlap_flag){
					// add alnSeg item
					ins_str = "";
					for(k=0; k<tmp_cigar_len; k++) ins_str += "=ACMGRSVTWYHKDBN"[bam_seqi(bam_get_seq(b), startQpos+k-1)];
					alnSegs.push_back(allocateAlnSeg(startRpos, startQpos, tmp_cigar_len, BAM_CINS, ins_str));
				}

				startQpos += tmp_cigar_len;
				op = bam_cigar_op(c[++i]);
				tmp_cigar_len = bam_cigar_oplen(c[i]);
				break;
			case BAM_CDEL:
				if(overlap_flag){
					// add alnSeg item
					del = "";
					for(k=0; k<tmp_cigar_len; k++) del += "=ACMGRSVTWYHKDBN"[bam_seqi(bam_get_seq(b), startQpos+k-1)];//something wrong?
					alnSegs.push_back(allocateAlnSeg(startRpos, startQpos, tmp_cigar_len, BAM_CDEL, del));
				}

				startRpos += tmp_cigar_len;
				op = bam_cigar_op(c[++i]);
				tmp_cigar_len = bam_cigar_oplen(c[i]);
				break;
			case BAM_CSOFT_CLIP:
			case BAM_CHARD_CLIP:
				if(overlap_flag){
					//str_tmp = to_string(tmp_cigar_len);
					alnSegs.push_back(allocateAlnSeg(startRpos, startQpos, tmp_cigar_len, op, to_string(tmp_cigar_len)));
				}

				if(op==BAM_CSOFT_CLIP) startQpos += tmp_cigar_len;
				if(++i<b->core.n_cigar){
					op = bam_cigar_op(c[i]);
					tmp_cigar_len = bam_cigar_oplen(c[i]);
				}
				break;
			case BAM_CEQUAL:
				if(overlap_flag){
					//str_tmp = "";
					alnSegs.push_back(allocateAlnSeg(startRpos, startQpos, tmp_cigar_len, BAM_CEQUAL, ""));
				}

				startQpos += tmp_cigar_len;
				startRpos += tmp_cigar_len;
				op = bam_cigar_op(c[++i]);
				tmp_cigar_len = bam_cigar_oplen(c[i]);
				break;
			case BAM_CDIFF:
				if(overlap_flag){
					diff = "";
					for(k=0; k<tmp_cigar_len; k++) diff += "=ACMGRSVTWYHKDBN"[bam_seqi(bam_get_seq(b), startQpos+k-1)];
					alnSegs.push_back(allocateAlnSeg(startRpos, startQpos, tmp_cigar_len, BAM_CDIFF, diff));
				}

				startQpos += tmp_cigar_len;
				startRpos += tmp_cigar_len;
				op = bam_cigar_op(c[++i]);
				tmp_cigar_len = bam_cigar_oplen(c[i]);
				break;
			case BAM_CREF_SKIP: // unexpected events
			case BAM_CPAD:
			case BAM_CBACK:
			default:
				cerr << __func__ << ", line=" << __LINE__ << ": invalid opflag " << op << endl;
				exit(1);
		}
		if(startRpos>endRefPos_paras) break;
	}
	return alnSegs;
}

vector<struct pafalnSeg*> generatePafAlnSegs(minimap2_aln_t* minimap2_aln_item, string &cons_seq, string &ref_seq){
	vector<struct pafalnSeg*> pafalnSegs;
	string cigar;
	int32_t seglen, relative_strand;
	uint32_t opflag, i;
	bool cigar_int_flag = false;
	int64_t startRpos, startQpos, startSubjectPos;
	int32_t op_num=0, op_nums=0, *seglens;
	uint32_t *opflags;
	char *opflags_char;

	cigar = minimap2_aln_item->cigar;
	for(i = 0; i < cigar.length(); i++){
		if(cigar.at(i) >= 'A' and cigar.at(i) <= 'Z'){
			op_nums++;
		}
	}
	seglens = new int32_t[op_nums]();
	opflags = new uint32_t[op_nums]();
	opflags_char = new char[op_nums]();
	//op_num = 0;
	//cout<<"cigar length:"<<cigar.length()<<endl;
	for(i = 0; i < cigar.length(); i++){
		if(cigar.at(i) >= '0' and cigar.at(i) <= '9'){
			if(cigar_int_flag){//804M
				//seglens.push_back(seglens[op_num]*10 + cigar.at(i) - '0');
				seglens[op_num] = seglens[op_num]*10 + cigar.at(i) - '0';
			}else{//8M
				seglens[op_num] = cigar.at(i) - '0';
				cigar_int_flag = true;
			}
		}else if(cigar.at(i) >= 'A' and cigar.at(i) <= 'Z'){
			opflags_char[op_num] = cigar.at(i);
			op_num++;
			if(cigar_int_flag) cigar_int_flag = false;
		}
	}
//	cout<<"op_num:"<<op_num<<endl;
	for(i = 0; i < (uint32_t)op_num; i++){
		switch(opflags_char[i]){
			case 'M':
				opflags[i] = BAM_CMATCH;
				break;
			case 'I':
				opflags[i] = BAM_CINS;
				break;
			case 'D':
				opflags[i] = BAM_CDEL;
				break;
			case 'S':
				opflags[i] = BAM_CSOFT_CLIP;
				break;
			case 'H':
				opflags[i] = BAM_CHARD_CLIP;
				break;
			case BAM_CREF_SKIP: // unexpected events
			case BAM_CPAD:
			case BAM_CBACK:
			default:
				cerr << __func__ << ", line=" << __LINE__ << ": invalid opflag " << opflags[i] << endl;
				exit(1);
		}
	}

	if(minimap2_aln_item->relative_strand==ALN_PLUS_ORIENT){ // plus strand
		i = 0;
		while(i < (uint32_t)op_num){
			if(i == 0){
				//startRpos = minimap2_aln_item->region_startRefPos + minimap2_aln_item->subject_start - 1;
				startRpos = minimap2_aln_item->region_startRefPos + minimap2_aln_item->subject_start;
				startQpos = minimap2_aln_item->query_start;
				startSubjectPos = minimap2_aln_item->subject_start;
			}else{
				switch (opflags[i-1]) {
					case BAM_CMATCH:
						startRpos = startRpos + seglens[i-1];
						startQpos = startQpos + seglens[i-1];
						startSubjectPos = startSubjectPos + seglens[i-1];
//						startRpos = startRpos + seglens[i-1] + 1;
//						startQpos = startQpos + seglens[i-1] + 1;
//						startSubjectPos = startSubjectPos + seglens[i-1] + 1;
						break;
					case BAM_CINS:
						startQpos = startQpos + seglens[i-1];
//						startQpos = startQpos + seglens[i-1] - 1;
						//startRpos = startRpos - 1;
						//startSubjectPos = startSubjectPos - 1;
						break;
					case BAM_CDEL:
						startRpos = startRpos + seglens[i-1];
						startSubjectPos = startSubjectPos + seglens[i-1];
//						startRpos = startRpos + seglens[i-1] - 1;
//						startSubjectPos = startSubjectPos + seglens[i-1] - 1;
						//startQpos = startQpos - 1;
						break;
					case BAM_CSOFT_CLIP:
					case BAM_CHARD_CLIP:
					case BAM_CEQUAL:
					case BAM_CDIFF:
						break;
					case BAM_CREF_SKIP:
						// unexpected events
					case BAM_CPAD:
					case BAM_CBACK:
					default:
						cerr << __func__ << ", line=" << __LINE__ << ": invalid opflag " << opflags[i-1] << endl;
						exit(1);
				}
			}

			seglen = seglens[i];
			opflag = opflags[i];
			relative_strand = minimap2_aln_item->relative_strand;
			if(startRpos<=minimap2_aln_item->region_endRefPos and startQpos<=minimap2_aln_item->query_end and startSubjectPos<=minimap2_aln_item->subject_end)
				pafalnSegs.push_back(allocatePafAlnSeg(startRpos, startQpos, startSubjectPos, seglen, opflag, cons_seq, ref_seq, relative_strand));
			else{
				cout << "^^^^^^^^^^^^ line=" << __LINE__ << ", startRpos=" << startRpos << ", startQpos=" << startQpos << ", startSubjectPos=" << startSubjectPos << endl;
			}
			i++;
		}
	}else{ // minus strand
		i = 0;
		while(i < (uint32_t)op_num){
			if(i == 0){
				//startRpos = minimap2_aln_item->region_startRefPos + minimap2_aln_item->subject_start - 1;
				startRpos = minimap2_aln_item->region_startRefPos + minimap2_aln_item->subject_start;
				startQpos = minimap2_aln_item->query_len - minimap2_aln_item->query_end + 1;
				startSubjectPos = minimap2_aln_item->subject_start;
			}else{
				switch (opflags[i-1]) {
					case BAM_CMATCH:
						startRpos = startRpos + seglens[i-1];
						startQpos = startQpos + seglens[i-1];
						startSubjectPos = startSubjectPos + seglens[i-1];
//						startRpos = startRpos + seglens[i-1] + 1;
//						startQpos = startQpos + seglens[i-1] + 1;
//						startSubjectPos = startSubjectPos + seglens[i-1] + 1;
						break;
					case BAM_CINS:
						startQpos = startQpos + seglens[i-1];
//						startQpos = startQpos + seglens[i-1] - 1;
//						startRpos = startRpos - 1;
//						startSubjectPos = startSubjectPos - 1;
						break;
					case BAM_CDEL:
						startRpos = startRpos + seglens[i-1];
						startSubjectPos = startSubjectPos + seglens[i-1];
//						startRpos = startRpos + seglens[i-1] - 1;
//						startSubjectPos = startSubjectPos + seglens[i-1] - 1;
//						startQpos = startQpos - 1;
						break;
					case BAM_CSOFT_CLIP:
					case BAM_CHARD_CLIP:
					case BAM_CEQUAL:
					case BAM_CDIFF:
						break;
					case BAM_CREF_SKIP:
						// unexpected events
					case BAM_CPAD:
					case BAM_CBACK:
					default:
						cerr << __func__ << ", line=" << __LINE__ << ": invalid opflag " << opflags[i-1] << endl;
						exit(1);
				}
			}
			seglen = seglens[i];
			opflag = opflags[i];
			relative_strand = minimap2_aln_item->relative_strand;

			pafalnSegs.push_back(allocatePafAlnSeg(startRpos, startQpos, startSubjectPos, seglen, opflag, cons_seq, ref_seq, relative_strand));

			i++;
		}
	}

	delete []seglens;
	delete []opflags;
	delete []opflags_char;

	return pafalnSegs;
}

struct alnSeg* allocateAlnSeg(int64_t startRpos, int64_t startQpos, int32_t seglen, uint32_t opflag, const string &seg_MD){
	struct alnSeg* seg = new struct alnSeg;
	if(!seg){
		cerr << __func__ << ": cannot allocate memory" << endl;
		exit(1);
	}
	seg->startRpos = startRpos;
	seg->startQpos = startQpos;
	seg->seglen = seglen;
	seg->opflag = opflag;
	seg->seg_MD = seg_MD;
	return seg;
}

struct pafalnSeg* allocatePafAlnSeg(int64_t startRpos, int64_t startQpos, int64_t startSubjectPos, int32_t seglen, uint32_t opflag, string &cons_seq, string &ref_seq, int32_t aln_orient){
	struct pafalnSeg* seg = new struct pafalnSeg;
	string ret;
	int32_t left_shift_num;

	if(!seg){
		cerr << __func__ << ": cannot allocate memory" << endl;
		exit(1);
	}

	left_shift_num = 1;
	if(startRpos-1<0 or startQpos-1<0 or startSubjectPos-1<0){
		left_shift_num = 0;
	}

//	seg->startRpos = startRpos;
//	seg->startQpos = startQpos;
//	seg->startSubjectPos = startSubjectPos;
	seg->startRpos = startRpos - left_shift_num;
	seg->startQpos = startQpos - left_shift_num;
	seg->startSubjectPos = startSubjectPos - left_shift_num;
	seg->seglen = seglen;
	seg->opflag = opflag;

	switch (opflag) {
		case BAM_CMATCH:
			seg->alt_seq = "N";
			seg->ref_seq = "N";
			break;
		case BAM_CINS:
			//seg->ref_seq = ref_seq.substr(startSubjectPos - 1, 1);
			seg->ref_seq = ref_seq.substr(startSubjectPos - left_shift_num, 1);
//			if(aln_orient==ALN_PLUS_ORIENT){
			//seg->alt_seq = cons_seq.substr(startQpos - 1, seglen);
			seg->alt_seq = cons_seq.substr(startQpos - left_shift_num, seglen + 1);
//				}else{
//					tmp_endPos = startQpos + seglen - 1;
//					tmp_startPos = startQpos;
//					tmp_startPos2 = cons_seq.size() - tmp_endPos + 1;
//					tmp_endPos2 = cons_seq.size() - tmp_startPos + 1;
//					ret = cons_seq.substr(tmp_startPos2-1, tmp_endPos2-tmp_startPos2+1);
//					reverseComplement(ret);
//					seg->alt_seq = ret;
//				}
			break;
		case BAM_CDEL:
			//seg->alt_seq = cons_seq.substr(startQpos - 1, 1);
			seg->alt_seq = cons_seq.substr(startQpos - left_shift_num, 1);
//			if(aln_orient==ALN_PLUS_ORIENT){
			//seg->ref_seq = ref_seq.substr(startSubjectPos - 1, seglen);
			seg->ref_seq = ref_seq.substr(startSubjectPos - left_shift_num, seglen + 1);
//			}else{
////				tmp_endPos = startSubjectPos + seglen - 1;
////				tmp_startPos = startSubjectPos;
////				tmp_startPos2 = cons_seq.size() - tmp_endPos + 1;
////				tmp_endPos2 = cons_seq.size() - tmp_startPos + 1;
//				ret = cons_seq.substr(startQpos - 1, seglen);
//				reverseComplement(ret);
//				seg->ref_seq = ret;
//			}
			break;
		case BAM_CSOFT_CLIP:
		case BAM_CHARD_CLIP:
			break;
		case BAM_CEQUAL:
		case BAM_CDIFF:
			seg->alt_seq = "N";
			seg->ref_seq = "N";
			break;
		case BAM_CREF_SKIP:
			// unexpected events
		case BAM_CPAD:
		case BAM_CBACK:
		default:
			cerr << __func__ << ", line=" << __LINE__ << ": invalid opflag " << opflag << endl;
			exit(1);
	}
	return seg;
}

// release the alignment segments
void destroyAlnSegs(vector<struct alnSeg*> &alnSegs){
	vector<struct alnSeg*>::iterator seg;
	for(seg=alnSegs.begin(); seg!=alnSegs.end(); seg++) delete *seg;
	vector<struct alnSeg*>().swap(alnSegs);
}

// get bam type from align data
int32_t getBamType(vector<bam1_t*> &alnDataVector){
	int32_t bamtype;
	for(size_t i=0; i<alnDataVector.size(); i++){
		bamtype = getBamTypeSingleItem(alnDataVector.at(i));
		if(bamtype!=BAM_INVALID)
			return bamtype;
	}

	cout<<"invalid bam file"<<endl;
	return BAM_INVALID;
}

// get bam type from clip align data
int32_t getBamType(vector<clipAlnData_t*> &clipAlnDataVector){
	int32_t bamtype;
	for(size_t i=0; i<clipAlnDataVector.size(); i++){
		bamtype = getBamTypeSingleItem(clipAlnDataVector.at(i)->bam);
		if(bamtype!=BAM_INVALID)
			return bamtype;
	}

	cout<<"invalid bam file"<<endl;
	return BAM_INVALID;
}

//estimate the kind of bamfile
int32_t getBamTypeSingleItem(bam1_t *b){
	uint32_t *c, op,i;
	uint8_t *flag;
	flag = bam_aux_get(b,"MD");
	if ( flag ){
		for(i=0; i<b->core.n_cigar; i++){
			c = bam_get_cigar(b);  // CIGAR
			op = bam_cigar_op(c[i]);
			if(op == BAM_CMATCH){
				//cout<<"there are MD flag and CMATCH in the bam file"<<endl;
				return BAM_CIGAR_NO_DIFF_MD;
			}
		}
		for(i=0; i<b->core.n_cigar; i++){
			c = bam_get_cigar(b);  // CIGAR
			op = bam_cigar_op(c[i]);
			if(op == BAM_CEQUAL || op == BAM_CDIFF){
				//cout<<"there are MD flag and CDIFF in the bam file"<<endl;
				return BAM_CIGAR_DIFF_MD;
			}
		}
	}else {
		for(i=0; i<b->core.n_cigar; i++){
			c = bam_get_cigar(b);  // CIGAR
			op = bam_cigar_op(c[i]);
			if(op == BAM_CMATCH){
				//cout<<"there is no MD flag but having CMATCH in the bam file"<<endl;
				return BAM_CIGAR_NO_DIFF_NO_MD;
			}
		}
		for(i=0; i<b->core.n_cigar; i++){
			c = bam_get_cigar(b);  // CIGAR
			op = bam_cigar_op(c[i]);
			if(op == BAM_CEQUAL || op == BAM_CDIFF){
				//cout<<"there is no MD flag but having CDIFF in the bam file"<<endl;
				return BAM_CIGAR_DIFF_NO_MD;
			}
		}
	}

	//cout<<"invalid bam file"<<endl;
	return BAM_INVALID;
}

// get the MD segs
vector<struct MD_seg*> extractMDSegs(bam1_t* b){
	vector<struct MD_seg*> segs_MD;
	char *p;
	string seg;

	p = bam_aux2Z(bam_aux_get(b, "MD"));  // MD
	while(*p){
		seg = "";
		if(*p=='^'){  // deletion
			p++;  // omit the '^'
			while((*p>='A' and *p<='Z') or (*p>='a' and *p<='z')) seg += *p++;
			// allocate the MD seg node, and add it to vector
			segs_MD.push_back(allocateMDSeg(seg, BAM_CDEL));
		}else if(*p>='0' and *p<='9'){  // match (BAM_CEQUAL)
			seg += *p++;
			while(*p>='0' and *p<='9') seg += *p++;
			// allocate the MD seg node, and add it to vector
			//if(seg!="0") segs_MD.push_back(allocateMDSeg(seg, BAM_CMATCH));
			if(seg!="0") segs_MD.push_back(allocateMDSeg(seg, BAM_CEQUAL));
		}else if((*p>='A' and *p<='Z') or (*p>='a' and *p<='z')){ // mismatch (BAM_CDIFF)
			seg += *p++;
			if((*p>='A' and *p<='Z') or (*p>='a' and *p<='z')){
				cerr << __func__ << ": invalid seg" << endl;
				exit(1);
			}

			// allocate the MD seg node, and add it to vector
			//segs_MD.push_back(allocateMDSeg(seg, MD_MISMATCH));
			segs_MD.push_back(allocateMDSeg(seg, BAM_CDIFF));
		}else{
			cerr << __func__ << ": invalid seg" << endl;
			exit(1);
		}
	}

	return segs_MD;
}

// allocate MD seg node
struct MD_seg* allocateMDSeg(string& seg, uint32_t opflag){
	struct MD_seg* seg_MD = new struct MD_seg();
	if(!seg_MD){
		cerr << __func__ << ": cannot allocate memory" << endl;
		exit(1);
	}
	seg_MD->seg = seg;
	seg_MD->opflag = opflag;
	switch(opflag){
		//case BAM_CMATCH: seg_MD->seglen = stoi(seg_MD->seg); break;
		case BAM_CEQUAL: seg_MD->seglen = stoi(seg_MD->seg); break;
		case BAM_CDEL:
		//case MD_MISMATCH: seg_MD->seglen = seg_MD->seg.size(); break;
		case BAM_CDIFF: seg_MD->seglen = seg_MD->seg.size(); break;
		default:
			cerr << __func__ << ": invalid opflag:" << opflag << endl; exit(1);
	}
	return seg_MD;
}

// destroy MD seg vector
void destroyMDSeg(vector<struct MD_seg*> &segs_MD){
	vector<struct MD_seg*>::iterator seg;
	for(seg=segs_MD.begin(); seg!=segs_MD.end(); seg++) delete *seg;
	vector<struct MD_seg*>().swap(segs_MD);
}

// output the aln seg items
void outputAlnSegs(vector<struct alnSeg*> &alnSegs){
	size_t i;
	vector<struct alnSeg*>::iterator seg;
	cout << "count: " << alnSegs.size() << endl;
	for(i=0, seg=alnSegs.begin(); seg!=alnSegs.end(); seg++, i++)
		cout << "[" << i << "]: startRpos=" << (*seg)->startRpos << ", startQpos=" << (*seg)->startQpos <<
				", seglen=" << (*seg)->seglen << ", opflag=" << size_t((*seg)->opflag) <<
				", seg_MD=" << (*seg)->seg_MD << endl;
}

// output the MD seg items
void outputMDSeg(vector<struct MD_seg*> &segs_MD){
	size_t i;
	vector<struct MD_seg*>::iterator seg;
	cout << "count: " << segs_MD.size() << endl;
	for(i=0, seg=segs_MD.begin(); seg!=segs_MD.end(); seg++, i++)
		cout << "[" << i << "]: seg=" << (*seg)->seg << ", opflag=" << size_t((*seg)->opflag) << endl;
}

// determine whether the alignment have the given cigar flag
bool haveOpflagCigar(bam1_t* b, uint32_t opfalg){
	bool flag = false;
	uint32_t i, *c, op;

	for(i=0; i<b->core.n_cigar; i++){
		c = bam_get_cigar(b);  // CIGAR
		op = bam_cigar_op(c[i]);
		if(op==opfalg){ flag = true; break;}
	}
	return flag;
}


// test alnSeg vector for single chromosome
void testAlnSegVec(string &inBamFile, faidx_t *fai){
	samFile *in = 0;
	bam_hdr_t *header;
	int32_t i, result, bam_type, seq_len, startPos, endPos;
	size_t count;
	bam1_t *b;
	string reg_str, chr_name, refseq;
	vector<struct alnSeg*> alnSegs;
	char *seq;

	if ((in = sam_open(inBamFile.c_str(), "r")) == 0) {
		cerr << __func__ << ": failed to open " << inBamFile << " for reading" << endl;
		exit(1);
	}

	if ((header = sam_hdr_read(in)) == 0) {
		cerr << __func__ << ": fail to read the header from " << inBamFile << endl;
		exit(1);
	}

	hts_idx_t *idx = sam_index_load(in, inBamFile.c_str()); // load index
	if (idx == 0) { // index is unavailable
		cerr << __func__ << ": random alignment retrieval only works for indexed BAM files.\n" << endl;
		exit(1);
	}

	for(i=0; i<header->n_targets; i++){
		chr_name = header->target_name[i];
		startPos = 1;
		endPos = header->target_len[i];
		reg_str = chr_name + ":" + to_string(startPos) + "-" + to_string(endPos);
		seq = fai_fetch(fai, reg_str.c_str(), &seq_len);
		refseq = seq;
		free(seq);

		hts_itr_t *iter = sam_itr_querys(idx, header, reg_str.c_str()); // parse a region in the format like `chr2:100-200'
		if (iter == NULL) { // region invalid or reference name not found
			int beg, end;
			if (hts_parse_reg(reg_str.c_str(), &beg, &end))
				cerr <<  __func__ << ": region " << reg_str << " specifies an unknown reference name." << endl;
			else
				cerr <<  __func__ << ": region " << reg_str << " could not be parsed." << endl;
			exit(1);
		}

		// fetch alignments
		bam_type = -1;
		count = 0;
		b = bam_init1();
		while ((result = sam_itr_next(in, iter, b)) >= 0) {
			count++;

			if(bam_type==-1) bam_type = getBamTypeSingleItem(b);
			if(bam_type==BAM_INVALID){
				cerr << __func__ << ": unknown bam type, error!" << endl;
				//exit(1);
			}else{

				if(!(b->core.flag & BAM_FUNMAP)){ // aligned
	//				qname = bam_get_qname(b);
	//				cin>>qname_com;
	//				qname_com = "b219a8ba_98667_6199";
	//				if(qname.compare(qname_com)==0)
	//					cout << qname << endl;
	//				else continue;

					switch(bam_type){
						case BAM_CIGAR_NO_DIFF_MD:
							alnSegs = generateAlnSegs(b);
							break;
						case BAM_CIGAR_NO_DIFF_NO_MD:
						case BAM_CIGAR_DIFF_MD:
						case BAM_CIGAR_DIFF_NO_MD:
							alnSegs = generateAlnSegs_no_MD(b, refseq, startPos, endPos);
							break;
						default:
							cerr << __func__ << ": unknown bam type, error!" << endl;
							exit(1);
					}// generate align segments

					// test single alnSeg vector
					checkAlnSegVecSingleQuery(alnSegs);

					//updateBaseInfo(baseArr, alnSegs); // update base information
					destroyAlnSegs(alnSegs);
				}
			}

			b = bam_init1();
		}

		bam_destroy1(b);
		if (result < -1) {
			cerr <<  __func__ << ": retrieval of region " << reg_str << " failed due to truncated file or corrupt BAM index file." << endl;
			exit(1);
		}

		hts_itr_destroy(iter);
	}

	// load align data from iteration
	//loadAlnDataFromIter(alnDataVector, in, header, iter, reg_str);

	hts_idx_destroy(idx); // destroy the BAM index
	bam_hdr_destroy(header);
	if(in) sam_close(in);
}

void checkAlnSegVecSingleQuery(vector<struct alnSeg*> &alnSegs){
	struct alnSeg *seg, *seg_prev;

	seg = NULL;
	for(size_t i=0; i<alnSegs.size(); i++){
		seg = alnSegs.at(i);
		if(seg->seglen==0){
			cout << "line=" << __LINE__ << ", opflag=" << seg->opflag << ", seglen=" << seg->seglen << ", seg_MD=" << seg->seg_MD << ", startRpos=" << seg->startRpos << ", startQpos=" << seg->startQpos << endl;
		}
	}

	for(size_t i=1; i<alnSegs.size(); i++){
		seg_prev = alnSegs.at(i-1);
		seg = alnSegs.at(i);
		if((seg_prev->opflag==BAM_CDIFF or seg_prev->opflag==BAM_CDEL) and (seg->opflag==BAM_CDIFF or seg->opflag==BAM_CDEL)){
			cout << "line=" << __LINE__ << ", opflag=" << seg->opflag << ", seglen=" << seg->seglen << ", seg_MD=" << seg->seg_MD << ", startRpos=" << seg->startRpos << ", startQpos=" << seg->startQpos << endl;
		}
	}
}

// get original query length
int32_t getOriginalQueryLen(bam1_t *b){
	int32_t querylen, left_clip_size, right_clip_size;
	uint32_t *c, op1, op2;

	c = bam_get_cigar(b);  // CIGAR
	// left clip
	op1 = bam_cigar_op(c[0]);
	if(op1==BAM_CSOFT_CLIP or op1==BAM_CHARD_CLIP) left_clip_size = bam_cigar_oplen(c[0]);
	else left_clip_size = 0;
	// left clip
	op2 = bam_cigar_op(c[b->core.n_cigar-1]);
	if(op2==BAM_CSOFT_CLIP or op2==BAM_CHARD_CLIP) right_clip_size = bam_cigar_oplen(c[b->core.n_cigar-1]);
	else right_clip_size = 0;

	// querylen
	querylen = b->core.l_qseq;
	if(op1==BAM_CHARD_CLIP)
		querylen += left_clip_size;
	if(op2==BAM_CHARD_CLIP)
		querylen += right_clip_size;

	return querylen;
}

int64_t getEndRefPosAlnSeg(int64_t startRpos, int32_t opflag, int32_t op_len){
	int64_t end_ref_pos = 0;

	switch(opflag){
		case BAM_CMATCH:
			end_ref_pos = startRpos + op_len - 1;
			break;
		case BAM_CINS:
			end_ref_pos = startRpos;
			break;
		case BAM_CDEL:
			end_ref_pos = startRpos + op_len - 1;
			break;
		case BAM_CSOFT_CLIP:
		case BAM_CHARD_CLIP:
			end_ref_pos = startRpos;
			break;
		case BAM_CEQUAL:
		case BAM_CDIFF:
			end_ref_pos = startRpos + op_len - 1;
			break;
		case BAM_CREF_SKIP:
			// unexpected events
		case BAM_CPAD:
		case BAM_CBACK:
		default:
			cerr << __func__ << ", line=" << __LINE__ << ": invalid opflag " << opflag << endl;
			exit(1);
	}

	return end_ref_pos;
}

int64_t getEndSubjectPosAlnSeg(int64_t startSubpos, int32_t opflag, int32_t op_len){
	int64_t end_sub_pos = 0;

	switch(opflag){
		case BAM_CMATCH:
			end_sub_pos = startSubpos + op_len - 1;
			break;
		case BAM_CINS:
			end_sub_pos = startSubpos;
			break;
		case BAM_CDEL:
			end_sub_pos = startSubpos + op_len - 1;
			break;
		case BAM_CSOFT_CLIP:
		case BAM_CHARD_CLIP:
			end_sub_pos = startSubpos;
			break;
		case BAM_CEQUAL:
		case BAM_CDIFF:
			end_sub_pos = startSubpos + op_len - 1;
			break;
		case BAM_CREF_SKIP:
			// unexpected events
		case BAM_CPAD:
		case BAM_CBACK:
		default:
			cerr << __func__ << ", line=" << __LINE__ << ": invalid opflag " << opflag << endl;
			exit(1);
	}

	return end_sub_pos;
}

int64_t getEndQueryPosAlnSeg(int64_t startQpos, int32_t opflag, int32_t op_len){
	int64_t end_query_pos = 0;

	switch(opflag){
		case BAM_CMATCH:
			end_query_pos = startQpos + op_len - 1;
			break;
		case BAM_CINS:
			end_query_pos = startQpos + op_len - 1;
			break;
		case BAM_CDEL:
			end_query_pos = startQpos;
			break;
		case BAM_CSOFT_CLIP:
		case BAM_CHARD_CLIP:
			end_query_pos = startQpos;
			break;
		case BAM_CEQUAL:
		case BAM_CDIFF:
			end_query_pos = startQpos + op_len - 1;
			break;
		case BAM_CREF_SKIP:
			// unexpected events
		case BAM_CPAD:
		case BAM_CBACK:
		default:
			cerr << __func__ << ", line=" << __LINE__ << ": invalid opflag " << opflag << endl;
			exit(1);
	}

	return end_query_pos;
}

bool isQcSigMatch(qcSig_t *qc_sig, qcSig_t *seed_qc_sig, double size_ratio_match_thres){
	bool match_flag = false;
	double size_ratio;
	int64_t ref_distance;

	if(qc_sig and seed_qc_sig){
		if(qc_sig->chrname.compare(seed_qc_sig->chrname)==0){
		 // same cigar_op, ratio of signature size less than 0.7
			if(qc_sig->cigar_op==seed_qc_sig->cigar_op){
				if(qc_sig->cigar_op_len>0 and seed_qc_sig->cigar_op_len>0) {
					if(qc_sig->cigar_op_len <= seed_qc_sig->cigar_op_len) size_ratio = (double)qc_sig->cigar_op_len / seed_qc_sig->cigar_op_len;
					else size_ratio = (double)seed_qc_sig->cigar_op_len / qc_sig->cigar_op_len;
					if(size_ratio>=size_ratio_match_thres){
						ref_distance = abs(qc_sig->ref_pos - seed_qc_sig->ref_pos);
						if(ref_distance<100){
							match_flag = true;
						}
					}
				}
			}
		}
	}
	return match_flag;
}

vector<qcSig_t*> extractQcSigsFromAlnSegsSingleQuery(struct querySeqInfoNode *query_seq_info_node, string &chrname, int64_t startSpanPos, int64_t endSpanPos, int32_t min_sv_size){
	vector<qcSig_t*> sig_vec;
	qcSig_t *qc_sig;
	vector<struct alnSeg*>::iterator seg;

	//int64_t startSpanPos_extend, endSpanPos_extend;

	//startSpanPos_extend = startSpanPos - 10;
	//endSpanPos_extend = endSpanPos + 10;

	//bool corrupt = false;
	for(seg=query_seq_info_node->query_alnSegs.begin(); seg!=query_seq_info_node->query_alnSegs.end(); seg++){
		switch((*seg)->opflag){
			case BAM_CINS:  // insertion in query
				if((*seg)->startRpos>=startSpanPos and (*seg)->startRpos<=endSpanPos){
					if((*seg)->seglen>=min_sv_size){
						qc_sig = allocateQcSigNode(*seg);
						qc_sig->chrname = chrname;
						sig_vec.push_back(qc_sig);
					}
				}
				break;
			case BAM_CDEL:  // deletion in query
				if((*seg)->startRpos>=startSpanPos and ((*seg)->startRpos + (*seg)->seglen -1)<=endSpanPos){
					if((*seg)->seglen>=min_sv_size){
						qc_sig = allocateQcSigNode(*seg);
						qc_sig->chrname = chrname;
						sig_vec.push_back(qc_sig);
					}
				}
				break;
			case BAM_CEQUAL:
				break;
			case BAM_CDIFF:
				break;
			case BAM_CSOFT_CLIP:// unexpected events
			case BAM_CHARD_CLIP:
			case BAM_CREF_SKIP:
			case BAM_CMATCH:
				break;
			case BAM_CPAD:
			case BAM_CBACK:
				cerr << __func__ << ", line=" << __LINE__ << ": invalid opflag " << (*seg)->opflag << endl;
				exit(1);
				break;
			default:
				cerr << __func__ << ", line=" << __LINE__ << ": invalid opflag " << (*seg)->opflag << endl;
//				cerr <<"op_len:" << (*seg)->seglen << endl;
//				cerr <<"startRpos:" << (*seg)->startRpos << endl;
//				cerr <<"qname:"<<query_seq_info_node->qname<<endl;
//				corrupt = true;
				exit(1);
				break;
		}
//		if(corrupt){
//			break;
//		}
	}

	return sig_vec;
}

// allocate query cluster signature
qcSig_t *allocateQcSigNode(struct alnSeg *aln_seg){
	qcSig_t *qc_sig = new qcSig_t();
	qc_sig->cigar_op = aln_seg->opflag;
	qc_sig->cigar_op_len = aln_seg->seglen;
	qc_sig->reg_contain_flag = false;
	qc_sig->ref_pos = aln_seg->startRpos;
	qc_sig->sig_id = -1;
	qc_sig->chrname = "";
	qc_sig->mate_qcSig = NULL;
	return qc_sig;
}

void destoryQuerySeqInfoAll(vector<struct querySeqInfoNode*> &query_seq_info_all){
	struct querySeqInfoNode *q_node;
	for(size_t i=0; i<query_seq_info_all.size(); i++){
		q_node = query_seq_info_all.at(i);
		destroyAlnSegs(q_node->query_alnSegs);
		delete q_node;
	}
	vector<struct querySeqInfoNode*>().swap(query_seq_info_all);
}



Time::Time() {
	timestamp = start = end = start_all = end_all = 0;
	time_tm = NULL;
	time(&start_all);
}

Time::~Time() {
}

string Time::getTime(){
	string time_str;
	time(&timestamp);	// get the Unix time stamp
	time_tm = localtime(&timestamp); // convert the Unix time to Year-Mon-Day...
	time_str = ctime(&timestamp);
	time_str.erase(time_str.find("\n"), 1);
	return time_str;
}

void Time::printTime(){
	cout << getTime() << endl;
}

void Time::setStartTime(){
	time(&start);
}

void Time::printElapsedTime(){
	double cost_sec, cost_hour, cost_min, cost_day;

	if(start!=0){
		time(&end);
		cost_sec = difftime(end, start);
		cost_min = cost_sec / 60;
		cost_hour = cost_min / 60;
		cost_day = cost_hour / 24;

		start = end = 0; // reset time

		if(cost_day>1)
			cout << "Elapsed time: " << cost_day << " days = " << cost_hour << " hours = " << cost_min << " minutes = " << cost_sec << " seconds" << endl;
		else if(cost_hour>1)
			cout << "Elapsed time: " << cost_hour << " hours = " << cost_min << " minutes = " << cost_sec << " seconds" << endl;
		else if(cost_min>1)
			cout << "Elapsed time: " << cost_min << " minutes = " << cost_sec << " seconds" << endl;
		else
			cout << "Elapsed time: " << cost_sec << " seconds" << endl;
	}else{
		cerr << "Please set the start time before calculating the running time" << endl;
		exit(1);
	}
}

void Time::printSubCmdElapsedTime(){
	double cost_sec, cost_hour, cost_min, cost_day;

	if(start!=0){
		time(&end);
		cost_sec = difftime(end, start);
		cost_min = cost_sec / 60;
		cost_hour = cost_min / 60;
		cost_day = cost_hour / 24;

		start = end = 0; // reset time

		if(cost_day>1)
			cout << "Sub-command elapsed time: " << cost_day << " days = " << cost_hour << " hours = " << cost_min << " minutes = " << cost_sec << " seconds" << endl << endl;
		else if(cost_hour>1)
			cout << "Sub-command elapsed time: " << cost_hour << " hours = " << cost_min << " minutes = " << cost_sec << " seconds" << endl << endl;
		else if(cost_min>1)
			cout << "Sub-command elapsed time: " << cost_min << " minutes = " << cost_sec << " seconds" << endl << endl;
		else
			cout << "Sub-command elapsed time: " << cost_sec << " seconds" << endl << endl;
	}else{
		cerr << "Please set the start time before calculating the running time" << endl;
		exit(1);
	}
}

void Time::printOverallElapsedTime(){
	double cost_sec, cost_hour, cost_min, cost_day;

	if(start_all!=0){
		time(&end_all);
		cost_sec = difftime(end_all, start_all);
		cost_min = cost_sec / 60;
		cost_hour = cost_min / 60;
		cost_day = cost_hour / 24;

		if(cost_day>1)
			cout << "Overall elapsed time: " << cost_day << " days = " << cost_hour << " hours = " << cost_min << " minutes = " << cost_sec << " seconds" << endl << endl;
		else if(cost_hour>1)
			cout << "Overall elapsed time: " << cost_hour << " hours = " << cost_min << " minutes = " << cost_sec << " seconds" << endl << endl;
		else if(cost_min>1)
			cout << "Overall elapsed time: " << cost_min << " minutes = " << cost_sec << " seconds" << endl << endl;
		else
			cout << "Overall elapsed time: " << cost_sec << " seconds" << endl << endl;
	}else{
		cerr << "Please set the overall start time before calculating the overall running time" << endl;
		exit(1);
	}
}
