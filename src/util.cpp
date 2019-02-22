#include "util.h"

// string split function
vector<string> split(const string& s, const string& delim)
{
    vector<string> elems;
    size_t pos = 0;
    size_t len = s.length();
    size_t delim_len = delim.length();
    if (delim_len == 0) return elems;
    while (pos < len)
    {
        int find_pos = s.find(delim, pos);
        if (find_pos < 0)
        {
            elems.push_back(s.substr(pos, len - pos));
            break;
        }
        elems.push_back(s.substr(pos, find_pos - pos));
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

	infile.open(infilename);
	if(!infile.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << infilename << endl;
		exit(1);
	}
	// read each line and save to the output file
	while(getline(infile, line))
		if(line.size()>0){
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

	struct stat fileStat;
	if (stat(contigfilename.c_str(), &fileStat) == 0)
		if(fileStat.st_size>0)
			flag = true;

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

// find the vector item
reg_t* findVarvecItem(size_t startRefPos, size_t endRefPos, vector<reg_t*> &varVec){
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

// find the vector item
vector<reg_t*> findVarvecItemAll(size_t startRefPos, size_t endRefPos, vector<reg_t*> &varVec){
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
reg_t* findVarvecItemExtSize(size_t startRefPos, size_t endRefPos, vector<reg_t*> &varVec, int32_t leftExtSize, int32_t rightExtSize){
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

// determine whether two regions are overlapped
bool isOverlappedReg(reg_t* reg1, reg_t* reg2){
	bool flag = false;
	if(reg1->chrname.compare(reg2->chrname)==0){
		if(isOverlappedPos(reg1->startRefPos, reg1->endRefPos, reg2->startRefPos, reg2->endRefPos))
			flag = true;
	}
	return flag;
}

// determine whether the given positions of two regions are overlapped
bool isOverlappedPos(size_t startPos1, size_t endPos1, size_t startPos2, size_t endPos2){
	bool flag = false;
	if((startPos1>=startPos2 and startPos1<=endPos2)
		or (endPos2>=startPos1 and endPos2<=endPos1)
		or (startPos2>=startPos1 and startPos2<=endPos1)
		or (endPos1>=startPos2 and endPos1<=endPos2))
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

	return overlap_flag;
}

// get the overlapped mate clip region
mateClipReg_t* getOverlappedMateClipReg(mateClipReg_t *mate_clip_reg_given, vector<mateClipReg_t*> &mateClipRegVec){
	mateClipReg_t *mate_clip_reg_overlapped = NULL, *clip_reg;

	for(size_t i=0; i<mateClipRegVec.size(); i++){
		clip_reg = mateClipRegVec.at(i);
		if(clip_reg!=mate_clip_reg_given and clip_reg->valid_flag and clip_reg->reg_mated_flag and clip_reg->sv_type==mate_clip_reg_given->sv_type){
			if(isOverlappedMateClipReg(mate_clip_reg_given, clip_reg)==true){ // same chrome
				mate_clip_reg_overlapped = clip_reg;
				break;
			}
		}
	}

	return mate_clip_reg_overlapped;
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
bool isInReg(size_t pos, vector<reg_t*> &vec){
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
		for(j=i+1; j<regVector.size(); ){
			reg2 = regVector.at(j);
			if(isOverlappedReg(reg1, reg2)){
				updateReg(reg1, reg2);
				delete reg2;
				regVector.erase(regVector.begin()+j);
			}else j++;
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



Time::Time() {
	timestamp = 0;
	time_tm = NULL;
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


