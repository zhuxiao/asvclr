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

// determine whether the two mate chip regions are the same chrome
bool isSameChrome(mateClipReg_t *mate_clip_reg1, mateClipReg_t *mate_clip_reg2){
	bool flag = false;
	if(mate_clip_reg1->leftClipReg->chrname.compare(mate_clip_reg1->rightClipReg->chrname)==0 and mate_clip_reg2->leftClipReg->chrname.compare(mate_clip_reg2->rightClipReg->chrname)==0 and mate_clip_reg1->leftClipReg->chrname.compare(mate_clip_reg2->leftClipReg->chrname)==0)
		flag = true;
	return flag;
}

// determine whether the two mate chip regions are overlapped
bool isOverlappedMateClipReg(mateClipReg_t *mate_clip_reg1, mateClipReg_t *mate_clip_reg2){
	bool flag = false;
	reg_t *left_reg1, *right_reg1, *left_reg2, *right_reg2;

	left_reg1 = mate_clip_reg1->leftClipReg;
	right_reg1 = mate_clip_reg1->rightClipReg;
	left_reg2 = mate_clip_reg2->leftClipReg;
	right_reg2 = mate_clip_reg2->rightClipReg;

	if(isSameChrome(mate_clip_reg1, mate_clip_reg2)==true){ // same chrome
		if(isOverlappedPos(mate_clip_reg1->leftMeanClipPos, mate_clip_reg1->rightMeanClipPos, mate_clip_reg2->leftMeanClipPos, mate_clip_reg2->rightMeanClipPos))
			flag = true;
	}else{ // different chrome
		if((left_reg1->chrname.compare(left_reg2->chrname)==0 and right_reg1->chrname.compare(right_reg2->chrname)==0 and isOverlappedPos(left_reg1->startRefPos, left_reg1->endRefPos, left_reg2->startRefPos, left_reg2->endRefPos) and isOverlappedPos(right_reg1->startRefPos, right_reg1->endRefPos, right_reg2->startRefPos, right_reg2->endRefPos))
			or (left_reg1->chrname.compare(right_reg2->chrname)==0 and right_reg1->chrname.compare(left_reg2->chrname)==0 and isOverlappedPos(left_reg1->startRefPos, left_reg1->endRefPos, right_reg2->startRefPos, right_reg2->endRefPos) and isOverlappedPos(right_reg1->startRefPos, right_reg1->endRefPos, left_reg2->startRefPos, left_reg2->endRefPos)))
			flag = true;
	}

	return flag;
}

// get the overlapped mate clip region
mateClipReg_t* getOverlappedMateClipReg(mateClipReg_t *mate_clip_reg_given, vector<mateClipReg_t*> &mateClipRegVec){
	mateClipReg_t *mate_clip_reg_overlapped = NULL, *clip_reg;

	for(size_t i=0; i<mateClipRegVec.size(); i++){
		clip_reg = mateClipRegVec.at(i);
		if(clip_reg!=mate_clip_reg_given and clip_reg->valid_flag and clip_reg->reg_mated_flag){
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


