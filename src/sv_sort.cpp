#include "sv_sort.h"
#include "util.h"


SV_item *constructSVItem(string &line){
	SV_item *item = NULL;
	vector<string> str_vec;
	string str_tmp;

	if(line.size()){
		str_vec = split(line, "\t");
		item = new SV_item();
		item->chrname = str_vec.at(0);
		if(str_vec.at(1).compare("-")==0) item->startPos = 0;
		else item->startPos = stoi(str_vec.at(1));
		if(str_vec.at(2).compare("-")==0) item->endPos = 0;
		else item->endPos = stoi(str_vec.at(2));
		item->info = line;

		item->sv_type = VAR_UNC;
		item->chrname2 = "";
		item->startPos2 = item->endPos2 = 0;
		item->valid_flag = true;

		str_tmp = str_vec.at(4);
		if(str_tmp.compare("INS")==0 or str_tmp.compare("insertion")==0){
			item->sv_type = VAR_INS;
		}else if(str_tmp.compare("DEL")==0 or str_tmp.compare("deletion")==0){
			item->sv_type = VAR_DEL;
		}else if(str_tmp.compare("DUP")==0 or str_tmp.compare("duplication")==0){
			item->sv_type = VAR_DUP;
		}else if(str_tmp.compare("INV")==0 or str_tmp.compare("inversion")==0){
			item->sv_type = VAR_INV;
		}else if(str_tmp.compare("SNV")==0 or str_tmp.compare("snv")==0){
			item->sv_type = VAR_SNV;
		}else if(str_tmp.compare("CNV")==0 or str_tmp.compare("cnv")==0){
			item->sv_type = VAR_CNV;
		}else{
			if(str_vec.size()>=10 and (str_vec.at(7).compare("TRA")==0 or str_vec.at(7).compare("translocation")==0 or str_vec.at(7).compare("BND")==0)){
				item->chrname2 = str_vec.at(3);
				if(str_vec.at(4).compare("-")==0) item->startPos2 = 0;
				else item->startPos2 = stoi(str_vec.at(4));
				if(str_vec.at(5).compare("-")==0) item->endPos2 = 0;
				else item->endPos2 = stoi(str_vec.at(5));
				if(str_vec.at(7).compare("TRA")==0 or str_vec.at(7).compare("translocation")==0) item->sv_type = VAR_TRA;
				else item->sv_type = VAR_BND;
			}else
				item->sv_type = VAR_MIX;
		}

		if(item->sv_type!=VAR_TRA and item->sv_type!=VAR_BND and item->sv_type!=VAR_MIX){
			item->sv_len = stoi(str_vec.at(5));
			item->altseq = str_vec.at(8);
		}
	}
	return item;
}

vector<SV_item*> loadDataBED(string &filename){
	vector<SV_item*> sv_vec;
	string line;
	SV_item *item;
	ifstream infile;
	infile.open(filename);
	if(!infile.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << filename << endl;
		exit(1);
	}

	while(getline(infile, line)){
		if(line.size()){
			if(line.at(0)!='#'){
				item = constructSVItem(line);
				if(item) sv_vec.push_back(item);
			}
		}
	}

	infile.close();
	return sv_vec;
}

bool isComma(string &seq){
	bool flag = false;
	char ch;
	for(size_t i=0; i<seq.size(); i++){
		ch = seq.at(i);
		if(ch==','){
			flag = true;
			break;
		}
	}
	return flag;
}

bool isSeq(string &seq){
	bool flag = false;
	bool comma_flag = false;
	char ch;
	for(size_t i=0; i<seq.size(); i++){
		ch = seq.at(i);
		switch(ch){
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
					case 'M':
					case 'm':
					case 'R':
					case 'r':
					case 'S':
					case 's':
					case 'V':
					case 'v':
					case 'W':
					case 'w':
					case 'Y':
					case 'y':
					case 'H':
					case 'h':
					case 'K':
					case 'k':
					case 'D':
					case 'd':
					case 'B':
					case 'b':
						flag = true;
						break;
					case ',':
						flag = true;
						comma_flag = true;
						break;
					default: cerr << __func__ << ": unknown base: " << ch << endl; exit(1);
				}
	}
	if(comma_flag) flag = false;
	return flag;
}

// determine whether the character is a base
bool isBase(const char ch){
	bool flag = false;
	switch(ch){
		case 'a':
		case 'c':
		case 'g':
		case 't':
		case 'A':
		case 'C':
		case 'G':
		case 'T':
		case 'N':
		case 'n':

		// mixed bases
		case 'M':
		case 'm':
		case 'R':
		case 'r':
		case 'S':
		case 's':
		case 'V':
		case 'v':
		case 'W':
		case 'w':
		case 'Y':
		case 'y':
		case 'H':
		case 'h':
		case 'K':
		case 'k':
		case 'D':
		case 'd':
		case 'B':
		case 'b': flag = true; break;
		default: cerr << __func__ << ": unknown base: " << ch << endl; exit(1);
	}

	return flag;
}

// get sv type
vector<string> getSVType(vector<string> &str_vec){
	string sv_type_str, sv_type_str1, sv_type_str2, str_tmp;
	size_t i, j;
	bool is_seq_flag_ref, is_seq_flag_alt, is_comma_flag;
	vector<string> comma_str_vec;
	vector<string> info_vec, sv_type_vec;
	vector<int32_t> sv_len_vec;
	int32_t sv_type_pos, SVTYPE_pos, vec_idx, str_pos, start_pos, end_pos;

	sv_type_str = "";
	for(i=3; i<str_vec.size(); i++){
		str_tmp = str_vec.at(i);
		if(str_tmp.compare("INS")==0 or str_tmp.compare("DEL")==0 or str_tmp.compare("DUP")==0 or str_tmp.compare("INV")==0 or str_tmp.compare("TRA")==0 or str_tmp.compare("BND")==0 or str_tmp.compare("INVTRA")==0 or str_tmp.compare("MIX")==0 or str_tmp.compare("MNP")==0 or str_tmp.compare("SNV")==0){
			sv_type_str = str_tmp;
		}else if(str_tmp.compare("<INS>")==0 or str_tmp.compare("<DEL>")==0 or str_tmp.compare("<DUP>")==0 or str_tmp.compare("<INV>")==0 or str_tmp.compare("<TRA>")==0 or str_tmp.compare("<BND>")==0 or str_tmp.compare("<INVTRA>")==0 or str_tmp.compare("<MIX>")==0 or str_tmp.compare("<MNP>")==0 or str_tmp.compare("<SNV>")==0){
			sv_type_str = str_tmp.substr(1, str_tmp.size()-2);
		}else if(str_tmp.compare("insertion")==0 or str_tmp.compare("deletion")==0 or str_tmp.compare("duplication")==0 or str_tmp.compare("inversion")==0 or str_tmp.compare("translocation")==0 or str_tmp.compare("snv")==0){
			sv_type_str = str_tmp;
			sv_type_vec.push_back(sv_type_str);
		}
	}

	if(sv_type_str.size()==0){
		vec_idx = -1; str_pos = -1;
		for(i=3; i<str_vec.size(); i++){
			str_tmp = str_vec.at(i);
			sv_type_pos = str_tmp.find("sv_type=");
			if(sv_type_pos!=-1){
				vec_idx = i;
				str_pos = sv_type_pos + 8;
				break;
			}else{
				SVTYPE_pos = str_tmp.find("SVTYPE=");
				if(SVTYPE_pos!=-1){
					vec_idx = i;
					str_pos = SVTYPE_pos + 7;
					break;
				}
			}
		}

		if(vec_idx!=-1 and str_pos!=-1){ // found
			str_tmp = str_vec.at(vec_idx);
			start_pos = str_pos; end_pos = -1;
			for(j=start_pos; j<str_tmp.size(); j++){ // search ';'
				if(str_tmp.at(j)==';'){
					end_pos = j - 1;
					break;
				}
			}
			if(end_pos==-1) end_pos = str_tmp.size() - 1;
			sv_type_str = str_tmp.substr(start_pos, end_pos-start_pos+1);
		}
	}

	if(sv_type_str.size()==0){// not found
		// check sequence
		is_seq_flag_ref = isSeq(str_vec.at(3));
		if(is_seq_flag_ref){
			is_seq_flag_alt = isSeq(str_vec.at(4));
			if(is_seq_flag_alt){
				if(str_vec.at(3).size()==1 and str_vec.at(4).size()==1)sv_type_str = "SNV";
				if(str_vec.at(3).size()<str_vec.at(4).size())sv_type_str = "INS";
				if(str_vec.at(3).size()>str_vec.at(4).size())sv_type_str = "DEL";
				sv_type_vec.push_back(sv_type_str);
			}else{
				str_tmp = str_vec.at(4);
				is_comma_flag = isComma(str_vec.at(4));
				if(is_comma_flag){
					comma_str_vec = split(str_vec.at(4), ",");
					if(str_vec.at(3).size()<comma_str_vec.at(0).size())sv_type_str1 = "INS";
					if(str_vec.at(3).size()>comma_str_vec.at(0).size())sv_type_str1 = "DEL";
					if(str_vec.at(3).size()<comma_str_vec.at(1).size())sv_type_str2 = "INS";
					if(str_vec.at(3).size()>comma_str_vec.at(1).size())sv_type_str2 = "DEL";
					sv_type_vec.push_back(sv_type_str1);
					sv_type_vec.push_back(sv_type_str2);
				}
			}
		}
	}

	return sv_type_vec;
}

// allocate SV item
SV_item *allocateSVItem(string &chrname, size_t startPos, size_t endPos, string &chrname2, size_t startPos2, size_t endPos2, string &sv_type_str, size_t sv_len, string &altseq, string &line){
	size_t sv_type;

	if(sv_type_str.compare("UNC")==0){
		sv_type = VAR_UNC;
	}else if(sv_type_str.compare("INS")==0 or sv_type_str.compare("insertion")==0){
		sv_type = VAR_INS;
	}else if(sv_type_str.compare("DEL")==0 or sv_type_str.compare("deletion")==0){
		sv_type = VAR_DEL;
	}else if(sv_type_str.compare("DUP")==0 or sv_type_str.compare("duplication")==0){
		sv_type = VAR_DUP;
	}else if(sv_type_str.compare("INV")==0 or sv_type_str.compare("inversion")==0){
		sv_type = VAR_INV;
	}else if(sv_type_str.compare("TRA")==0 or sv_type_str.compare("translocation")==0){
		sv_type = VAR_TRA;
	}else if(sv_type_str.compare("BND")==0){
		sv_type = VAR_BND;
	}else if(sv_type_str.compare("INVTRA")==0){
		sv_type = VAR_INV_TRA;
	}else if(sv_type_str.compare("MIX")==0){
		sv_type = VAR_MIX;
	}else if(sv_type_str.compare("MNP")==0){
		sv_type = VAR_MNP;
	}else if(sv_type_str.compare("SNV")==0){
		sv_type = VAR_SNV;
	}else if(sv_type_str.compare("CNV")==0 or sv_type_str.compare("cnv")==0){
		sv_type = VAR_CNV;
	}else{
		sv_type = VAR_UNC;
	}

	SV_item *item = new SV_item();
	item->chrname = chrname;
	item->startPos = startPos;
	item->endPos = endPos;
	item->chrname2 = chrname2;
	item->startPos2 = startPos2;
	item->endPos2 = endPos2;
	item->sv_type = sv_type;
	item->sv_len = sv_len;
	item->altseq = altseq;
	item->info = line;
	item->valid_flag = true;
	return item;
}

vector<SV_item*> loadDataVcf(string &filename){
	vector<SV_item*> sv_vec;
	ifstream infile;
	string line, chrname, start_pos_str, endpos_str, chrname2, start_pos_str2, endpos_str2, sv_type_str, sv_type_str1, sv_type_str2, sv_len_str, altseq, bnd_str, bnd_pos_str;
	size_t i, start_pos, endpos, start_pos2, endpos2, sv_len;
	vector<string> str_vec, sv_type_vec, info_vec, sub_info_vec, bnd_pos_vec;
	SV_item *sv_item;

	infile.open(filename);
	if(!infile.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << filename << endl;
		exit(1);
	}

	// convert
	while(getline(infile, line)){
		if(line.size()){
			if(line.at(0)!='#'){

				str_vec = split(line, "\t");

				// get locations
				chrname = str_vec.at(0);
				start_pos_str = str_vec.at(1);
				altseq = str_vec.at(4);

				chrname2 = endpos_str = sv_type_str = "";
				info_vec = split(str_vec.at(7), ";");

				for(i=0; i<info_vec.size(); i++){
					sub_info_vec = split(info_vec.at(i), "=");
					if(sub_info_vec.at(0).compare("CHR2")==0)
						chrname2 = sub_info_vec.at(1);
					if(sub_info_vec.at(0).compare("END")==0)
						endpos_str = sub_info_vec.at(1);
					if(sub_info_vec.at(0).compare("SVTYPE")==0)
						sv_type_str = sub_info_vec.at(1);
					if(sub_info_vec.at(0).compare("SVLEN")==0)
						sv_len_str = sub_info_vec.at(1);
				}

				// get sv type
				if(sv_type_str.size()==0){
					sv_type_vec = getSVType(str_vec);
					if(sv_type_vec.size()==1){
						sv_type_str = sv_type_vec.at(0);
					}else{
						sv_type_str1 = sv_type_vec.at(0);
						sv_type_str2 = sv_type_vec.at(1);
					}
				}

				if(sv_type_str.size()>0){
					if(sv_type_str.compare("TRA")==0 or sv_type_str.compare("BND")==0){ // TRA or BND
						if(sv_type_str.compare("BND")==0){ // BND
							bnd_str = str_vec.at(4);
							if(bnd_str.at(0)==']' or bnd_str.at(0)=='[')
								bnd_pos_str = bnd_str.substr(1, bnd_str.size()-3);
							else if(isBase(bnd_str.at(0)))
								bnd_pos_str = bnd_str.substr(2, bnd_str.size()-3);
							bnd_pos_vec = split(bnd_pos_str, ":");
							chrname2 = bnd_pos_vec.at(0);
							endpos_str = bnd_pos_vec.at(1);
						}

						start_pos_str2 = endpos_str;
						endpos_str2 = start_pos_str2;
						endpos_str = start_pos_str;
					}else { // INS, DEL, INV, DUP, UNC
						start_pos_str2 = "0";
						endpos_str2 = "0";
					}

					start_pos = stoi(start_pos_str);
					endpos = stoi(endpos_str);
					start_pos2 = stoi(start_pos_str2);
					endpos2 = stoi(endpos_str2);
					if(sv_len_str.size()>0) sv_len = stoi(sv_len_str);
					else sv_len = 0;
					sv_item = allocateSVItem(chrname, start_pos, endpos, chrname2, start_pos2, endpos2, sv_type_str, sv_len, altseq, line);
					sv_vec.push_back(sv_item);

				}else{
					cout << line << endl;
					cout << "missing SVTYPE information" << endl;
				}
			}
		}
	}
	infile.close();

	return sv_vec;
}

void destroyData(vector<SV_item*> &sv_vec){
	vector<SV_item*>::iterator it;
	for(it=sv_vec.begin(); it!=sv_vec.end(); it++)
		delete *it;
	vector<SV_item*>().swap(sv_vec);
}

//get chrnames
set<string> getChrnames(vector<SV_item*> &dataset){
	set<string> chrname_set;
	size_t i;
	SV_item *item;
	set<string>::iterator iter;

	for(i=0;i<dataset.size();i++){
		item = dataset.at(i);
		if(chrname_set.find(item->chrname)==chrname_set.end()) chrname_set.insert(item->chrname);
	}

	return chrname_set;
}

vector<string> sortChrnames(set<string> &chrname_set){
	int8_t *selected_flag_array;
	size_t i, j;
	int32_t idx_chr;
	set<string>::iterator iter;
	string chr_str1, chr_str2, head_str, head_str_chr, chrname;
	vector<string> chrname_vec_sorted, chr_str_vec, chrname_vec;
	bool have_chr_prefix_flag;

	chr_str1 = "chr1_chr2_chr3_chr4_chr5_chr6_chr7_chr8_chr9_chr10_chr11_chr12_chr13_chr14_chr15_chr16_chr17_chr18_chr19_chr20_chr21_chr22_chrX_chrY_chrM";
	chr_str2 = "1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_X_Y_MT";

	for(iter=chrname_set.begin(); iter!=chrname_set.end(); ++iter) chrname_vec.push_back(*iter);
	selected_flag_array = (int8_t*) calloc(chrname_vec.size(), sizeof(int8_t));

	// determine the have_chr_prefix_flag
	have_chr_prefix_flag = false;
	for(i=0; i<chrname_vec.size(); i++){
		if(chrname_vec.at(i).compare("chr1")==0){
			have_chr_prefix_flag = true;
			break;
		}
	}

	if(have_chr_prefix_flag) chr_str_vec = split(chr_str1, "_");
	else chr_str_vec = split(chr_str2, "_");

	// select chromosomes by 'chr'
	for(i=0; i<chr_str_vec.size(); i++){
		idx_chr = -1;
		for(j=0; j<chrname_vec.size(); j++){
			chrname = chrname_vec.at(j);
			if(chrname.compare(chr_str_vec.at(i))==0){
				idx_chr = j;
				break;
			}
		}
		if(idx_chr!=-1){
			chrname = chrname_vec.at(idx_chr);
			chrname_vec_sorted.push_back(chrname);
			selected_flag_array[idx_chr] = 1;
		}
	}

	// select chromosomes by 'chrA_*'
	for(i=0; i<chr_str_vec.size(); i++){
		head_str = chr_str_vec.at(i) + "_";
		for(j=0; j<chrname_vec.size(); j++){
			chrname = chrname_vec.at(j);
			head_str_chr = chrname.substr(0, head_str.size());
			if(head_str_chr.compare(head_str)==0){
				chrname_vec_sorted.push_back(chrname);
				selected_flag_array[j] = 1;
			}
		}
	}

	// add unselected items
	for(i=0; i<chrname_vec.size(); i++){
		if(selected_flag_array[i]==0){
			chrname = chrname_vec.at(i);
			chrname_vec_sorted.push_back(chrname);
		}
	}

	vector<string>().swap(chrname_vec);	// free memory

	free(selected_flag_array);

	return chrname_vec_sorted;
}

//judge chrname is same
bool IsSameChrname(string &chrname1, string &chrname2){
	string chr = "chr";
	if(chrname1.compare(chrname2)==0 or (chr + chrname1).compare(chrname2)==0 or chrname1.compare(chr + chrname2)==0){
		return true;
	}
	return false;
}

vector<SV_item*> getItemsByChr(string &chrname, vector<SV_item*> &dataset){
	vector<SV_item*> result;
	size_t i;
	SV_item *item;

	for(i=0; i<dataset.size(); i++){
		item = dataset.at(i);
		if(IsSameChrname(item->chrname, chrname)) result.push_back(item);
	}

	return result;
}

vector<vector<SV_item*>> constructSubsetByChrOp(vector<SV_item*> &sv_vec, vector<string> &chrname_vec){
	vector< vector<SV_item*> > result;
	string chr;
	vector<SV_item*> subset;

	for(size_t i=0; i<chrname_vec.size(); i++){
		chr = chrname_vec.at(i);
		subset = getItemsByChr(chr, sv_vec);
		result.push_back(subset);
	}

	return result;
}

//generate datasets
vector<vector<SV_item*>> constructSubsetByChr(vector<SV_item*> &sv_vec){
	vector< vector<SV_item*> > result;
	set<string> chrname_set;
	vector<string> chrname_vec_sorted;

	chrname_set = getChrnames(sv_vec);
	chrname_vec_sorted = sortChrnames(chrname_set); // sort chromosome names
	result = constructSubsetByChrOp(sv_vec, chrname_vec_sorted);

	return result;
}

bool sortFunSameChr(const SV_item *item1, const SV_item *item2){
//	int32_t flag;

	if(item1->startPos!=item2->startPos){
		return item1->startPos < item2->startPos;
	}else { // equal

//		if(item1->endPos != item2->endPos){
			return item1->endPos < item2->endPos;
//		}else{ // equal
//			cout << item1->info << "---------" << item1->sv_type << endl;
//			cout << item2->info << "---------" << item2->sv_type << endl << endl;

//			if((item1->sv_type==VAR_TRA or item1->sv_type==VAR_BND) and (item2->sv_type==VAR_TRA or item2->sv_type==VAR_BND)){ // TRA or BND
//				if(item1->chrname2.size()==0 and item2->chrname2.size()==0) return true; // chrname2 is null
//				else{ // chrname2 is not null
//					flag = item1->chrname2.compare(item2->chrname2);
//					if(flag<0) return true;
//					else if(flag>0) return false;
//					else { // equal
//						if(item1->startPos2!=item2->startPos2) return item1->startPos2 < item2->startPos2;
//						else return item1->endPos2 < item2->endPos2;
//					}
//				}
////				return false;
//			}else{
//				if(item1->sv_type!=VAR_TRA and item1->sv_type!=VAR_BND) return true;
//				else return false;
//				return true;
//			}
//		}
	}
}

void sortSubset(vector<SV_item*> &sv_vec){
	sort(sv_vec.begin(), sv_vec.end(), sortFunSameChr);
}

void sortSVitem(vector<vector<SV_item*>> &subsets){
	for(size_t i=0;i<subsets.size();i++){
		sortSubset(subsets.at(i));
	}
}

// remove duplicated items
void rmDupSVitem(vector<vector<SV_item*>> &subsets, double size_ratio_thres, double identity_thres){
	cout << "remove duplicated items ..." << endl;
	for(size_t i=0;i<subsets.size();i++)
		rmDupSVitemSubset(subsets.at(i), size_ratio_thres, identity_thres);
}

void rmDupSVitemSubset(vector<SV_item*> &sv_vec, double size_ratio_thres, double identity_thres){
	SV_item *item1, *item2;
	double val, max_len, secmax_len, size_ratio;

	for(size_t i=1;i<sv_vec.size(); i++){
		item1 = sv_vec.at(i-1);
		item2 = sv_vec.at(i);
		if(item1->sv_type==item2->sv_type and item1->chrname.compare(item2->chrname)==0 and item1->startPos==item2->startPos and item1->endPos==item2->endPos){
			if(item1->sv_len==item2->sv_len) size_ratio = 1;
			else{
				max_len = item1->sv_len>item2->sv_len ? (max_len=item1->sv_len, secmax_len=item2->sv_len) : (max_len=item2->sv_len, secmax_len=item1->sv_len);
				size_ratio = secmax_len / max_len;
			}

			if(size_ratio>=size_ratio_thres){
				val = computeVarseqIdentity(item1->altseq, item2->altseq);
				if(val>=identity_thres) item2->valid_flag = false;
			}
		}
	}
}

