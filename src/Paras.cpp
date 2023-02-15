#include <string.h>
#include <vector>
#include <htslib/hts.h>
#include <cmath>
#include "Paras.h"
#include "util.h"
#include <algorithm>

// Constructor with parameters
Paras::Paras(){
	init();
}

// Constructor with parameters
Paras::Paras(int argc, char **argv){
	init();
	if(parseParas(argc, argv)!=0) exit(1);

	// check bam file
	if(checkBamFile()!=0) exit(1);
}

//Destructor
Paras::~Paras(){
	if(!limit_reg_vec.empty()) destroyLimitRegVector(limit_reg_vec);
}

// initialization
void Paras::init(){
	command = "";
	inBamFile = "";
	outFilePrefix = "";
	outDir = OUT_DIR;
	sample = SAMPLE_DEFAULT;
	pg_cmd_str = "";
	num_threads = 0;
	delete_reads_flag = true;
	keep_failed_reads_flag = reassemble_failed_work_flag = false;
	maskMisAlnRegFlag = false;
	assemChunkSize = ASM_CHUNK_SIZE_INDEL;
	assemSideExtSize = ASM_CHUNK_SIZE_EXT_INDEL;
	minConReadLen = MIN_CONS_READ_LEN;
	technology = SEQUENCING_TECH_DEFAULT;
	include_decoy = false;

	min_ins_size_filt = 0;
	min_del_size_filt = 0;
	min_clip_size_filt = 0;
	min_ins_num_filt = 0;
	min_del_num_filt = 0;
	min_clip_num_filt = 0;

	mean_read_len = total_read_num_est = 0;

	reg_sum_size_est = 0;
	max_reg_sum_size_est = MAX_REG_SUM_SIZE_EST;

	min_input_cov_canu = MIN_INPUT_COV_CANU;
	expected_cov_assemble = EXPECTED_COV_ASSEMBLE;
	max_ultra_high_cov = MAX_ULTRA_HIGH_COV_THRES;

	assemble_reg_preDone_num = assemble_reg_work_total = assemble_reg_workDone_num = 0;
	num_parts_progress = NUM_PARTS_PROGRESS;
	num_threads_per_assem_work = NUM_THREADS_PER_ASSEM_WORK;

	canu_version = getCanuVersion();
	if(canu_version.empty()){
		cerr << "Cannot find the 'Canu' assembler, please make sure it is correctly installed and the executable file 'canu' or its soft link is included in the '$PATH' directory." << endl;
		exit(1);
	}

	monitoring_proc_names_assemble = DEFAULT_MONITOR_PROC_NAMES_ASSEMBLE;
	max_proc_running_minutes_assemble = MAX_PROC_RUNNING_MINUTES_ASSEMBLE;
	monitoring_proc_names_call = DEFAULT_MONITOR_PROC_NAMES_CALL;
	max_proc_running_minutes_call = MAX_PROC_RUNNING_MINUTES_CALL;

	mem_total = getMemInfo("MemTotal", 2);
	swap_total = getMemInfo("SwapTotal", 2);
	if(mem_total<0 or swap_total<0){
		cerr << "line=" << __LINE__ << ", mem_total=" << mem_total << ", swap_total=" << swap_total << ", cannot get the supplied memory information, error." << endl;
		exit(1);
	}
	extend_total = mem_total<swap_total ? mem_total : swap_total;

	call_work_num = 0;

	gt_min_sig_size = GT_SIG_SIZE_THRES;
	gt_size_ratio_match = GT_SIZE_RATIO_MATCH_THRES;
	gt_min_alle_ratio = GT_MIN_ALLE_RATIO_THRES;
	gt_max_alle_ratio = GT_MAX_ALLE_RATIO_THRES;
}

// get the Canu program version
string Paras::getCanuVersion(){
	FILE *stream;
	char tmp[256], info[256] = {0};
	string canu_version_str = "";

	sprintf(tmp, "canu -version | awk '$1 ~/(Canu)|(canu)/' | awk '{print $2}'");
	stream = popen(tmp, "r");
	if(fread(info, 1, sizeof(info), stream)>0){
		canu_version_str = info;
		if(canu_version_str.at(canu_version_str.size()-1)=='\n') canu_version_str.at(canu_version_str.size()-1) = '\0';
	}
	pclose(stream);

	return canu_version_str;
}

// check whether the canu version is recommended
bool Paras::isRecommendCanuVersion(string &canu_version, const string &recommend_version){
	bool flag = true;
	int32_t result;

	result = canu_version.compare(recommend_version);
	if(result<0) flag = false;

	return flag;
}

// check Bam file, and generate the BAM index if it is unavailable
int Paras::checkBamFile(){
	samFile *in = 0;
	string idx_filename;
	int ret = 1;

	if ((in = sam_open(inBamFile.c_str(), "r")) == 0) {
		cerr << __func__ << ": failed to open " << inBamFile.c_str() << " for reading" << endl;
		exit(1);
	}

	hts_idx_t *idx = sam_index_load(in, inBamFile.c_str()); // load index
	if (idx == NULL) { // index is unavailable, then generate it
		if(in) sam_close(in);

		cout << __func__ << ": BAM index is unavailable, now generate it, please wait ...\n" << endl;

		// construct the index file name
		if(inBamFile.substr(inBamFile.size()-4).compare(".bam")==0)
			idx_filename = inBamFile.substr(0, inBamFile.size()-4) + ".bai";
		else
			idx_filename = inBamFile + ".bai";

		ret = sam_index_build3(inBamFile.c_str(), idx_filename.c_str(), 0, num_threads);
		switch (ret) {
			case 0:
				break;
			case -1:
				cerr << __func__ << ", error occurred while opening " << inBamFile << endl;
				exit(-1);
			case -2:
				cerr << __func__ << ", failed to open " << inBamFile << endl;
				exit(-2);
			case -3:
				cerr << __func__ << inBamFile << " is in a format that cannot be usefully indexed" << endl;
				exit(-3);
			case -4:
				cerr << __func__ << ", failed to create or write index" << endl;
				exit(-4);
			default:
				cerr << __func__ << ", failed to create index for " << inBamFile << endl;
				exit(1);
		}
	}else{
		ret = 0;
		hts_idx_destroy(idx); // destroy the BAM index
		if(in) sam_close(in);
	}

	return ret;
}

// parse the parameters
int Paras::parseParas(int argc, char **argv){
	if (argc < 2) { showUsage(); return 1; }

	if (strcmp(argv[1], "-h") == 0 or strcmp(argv[1], "help") == 0 or strcmp(argv[1], "--help") == 0) {
		if (argc == 2) { showUsage(); exit(0); }
		argv++;
		argc = 2;
	}

	// save program command line
	pg_cmd_str = getPgCmd(argc, argv);

	if (strcmp(argv[1], "detect")==0){
		if(argc==2){ showDetectUsage(); exit(0); }
		command = "detect";
		return parseDetectParas(argc-1, argv+1);
	}else if(strcmp(argv[1], "assemble")==0){
		if(argc==2){ showAssembleUsage(); exit(0); }
		command = "assemble";
		return parseAssembleParas(argc-1, argv+1);
	}else if(strcmp(argv[1], "call")==0){
		if(argc==2){ showCallUsage(); exit(0); }
		command = "call";
		return parseCallParas(argc-1, argv+1);
	}else if(strcmp(argv[1], "all")==0){
		if(argc==2){ showAllUsage("all"); exit(0); }
		command = "all";
		return parseAllParas(argc-1, argv+1, "all");
	}else if(strcmp(argv[1], "detect-assemble")==0){
		if(argc==2){ showDetectAssembleUsage(); exit(0); }
		command = "detect-assemble";
		return parseDetectAssembleParas(argc-1, argv+1);
	}else if (strcmp(argv[1], "--version") == 0 or strcmp(argv[1], "-v") == 0) {
		showVersion();
		exit(0);
	}else{
		cerr << "Error: invalid command " << argv[1] << endl << endl;
		showUsage(); return 1;
	}
}

// get program command string
string Paras::getPgCmd(int argc, char **argv){
	string pg_cmd_str, tmp_str;
	size_t found_idx;

	tmp_str = argv[0];
	if((found_idx=tmp_str.find_last_of("/"))!=string::npos)
		pg_cmd_str = tmp_str.substr(found_idx+1);
	else pg_cmd_str = tmp_str;
	for(int i=1; i<argc; i++) pg_cmd_str = pg_cmd_str + " " + argv[i];

	return pg_cmd_str;
}

// show version
void Paras::showVersion(){
	cout << PROG_NAME << " " << PROG_VERSION << endl;
	cout << "Using htslib " << hts_version() << endl;
}

// parse the parameters for detect command
int Paras::parseDetectParas(int argc, char **argv){
	int opt, threadNum_tmp = 0, option_index;
	blockSize = BLOCKSIZE;
	slideSize = SLIDESIZE;
	min_sv_size_usr = MIN_SV_SIZE_USR;
	max_sv_size_usr = MAX_SV_SIZE_USR;
	minReadsNumSupportSV = MIN_SUPPORT_READS_NUM_EST;
	//minReadsNumSupportSV = -1;
	//maxVarRegSize = MAX_VAR_REG_SIZE;
	minClipEndSize = MIN_CLIP_END_SIZE;
	outDir = OUT_DIR;
	//max_seg_size_ratio_usr = MAX_SEG_SIZE_RATIO;
	simpleReg_t *simple_reg;
	string simple_reg_str, opt_name_str;

	static struct option lopts[] = {
//		{ "daemon", no_argument, NULL, 'D' },
//		{ "dir", required_argument, NULL, 'd' },
//		{ "out", required_argument, NULL, 'o' },
//		{ "log", required_argument, NULL, 'l' },
		{ "sample", required_argument, NULL, 0 },
		{ "mask-noisy-region", no_argument, NULL, 0 },
		{ "include-decoy", no_argument, NULL, 0 },
		{ "version", no_argument, NULL, 'v' },
		{ "help", no_argument, NULL, 'h' },
		{ NULL, 0, NULL, 0 }
	};

	while( (opt = getopt_long(argc, argv, ":b:s:m:M:e:o:p:t:vh", lopts, &option_index)) != -1 ){
		switch(opt){
			case 'b': blockSize = stoi(optarg); break;
			case 's': slideSize = stoi(optarg); break;
			case 'm': min_sv_size_usr = stoi(optarg); break;
			case 'M': max_sv_size_usr = stoi(optarg); break;
			case 'n': minReadsNumSupportSV = stoi(optarg); break;
			case 'e': minClipEndSize = stoi(optarg); break;
			case 'o': outDir = optarg; break;
			case 'p': outFilePrefix = optarg; break;
			case 't': threadNum_tmp = stoi(optarg); break;
			case 'v': showVersion(); exit(0);
			case 'h': showDetectUsage(); exit(0);
			case '?':
				if(optopt) cout << "unknown option '-" << (char)optopt << "'" << endl;
				else{ // Bad long option
                    if (optind > 0 && strncmp(argv[optind-1], "--", 2) == 0)
                        cout << "unknown option '" << argv[optind-1] << "'" << endl;
				}
				exit(1);
			case ':':
				if(optopt) cout << "the option '-" << (char)optopt << "' needs a value" << endl;
				else{
                    if (optind > 0 && strncmp(argv[optind-1], "--", 2) == 0)
                        cout << "the option '" << argv[optind-1] << "' needs a value" << endl;
				}
				exit(1);
			case 0: // long options
				if(parse_long_opt(option_index, optarg, lopts)!=0){
					showDetectUsage();
					return 1;
				}
				break;
			default:
				cout << "Error: please specify the correct options for 'detect' command" << endl;
				showDetectUsage();
				return 1;
		}
	}

	load_from_file_flag = true;
	if(threadNum_tmp==0) num_threads = sysconf(_SC_NPROCESSORS_ONLN);
	else num_threads = (threadNum_tmp>=sysconf(_SC_NPROCESSORS_ONLN)) ? sysconf(_SC_NPROCESSORS_ONLN) : threadNum_tmp;
	maxVarRegSize = max_sv_size_usr;

	outDir = deleteTailPathChar(outDir);

	opt = argc - optind; // the reference and BAM file on the command line
	if(opt>=2) {
		refFile = argv[optind];
		inBamFile = argv[optind+1];

		for(int i=optind+2; i<argc; i++){
			simple_reg_str = argv[i];
			simple_reg = allocateSimpleReg(simple_reg_str);
			if(simple_reg) limit_reg_vec.push_back(simple_reg);
			else{
				cout << "Error: Please specify the correct genomic regions to be processed." << endl << endl;
				showDetectUsage();
				return 1;
			}
		}
		if(limit_reg_vec.size()) limit_reg_process_flag = true;
	}else{
		cout << "Error: Please specify the reference file and coordinate-sorted BAM file." << endl << endl;
		showDetectUsage();
		return 1;
	}

	return 0;
}

// parse the parameters for assemble command
int Paras::parseAssembleParas(int argc, char **argv){
	int opt, threadNum_tmp = 0, option_index;
	blockSize = BLOCKSIZE;
	slideSize = SLIDESIZE;
	//assemSlideSize = ASSEM_SLIDE_SIZE;
	min_sv_size_usr = MIN_SV_SIZE_USR;
	max_sv_size_usr = MAX_SV_SIZE_USR;
	minReadsNumSupportSV = MIN_SUPPORT_READS_NUM_EST;
	max_seg_size_ratio_usr = MAX_SEG_SIZE_RATIO;
	maxVarRegSize = MAX_VAR_REG_SIZE;
	assemChunkSize = ASM_CHUNK_SIZE_INDEL;
	assemSideExtSize = ASM_CHUNK_SIZE_EXT_INDEL;
	assemSideExtSizeClip = ASM_CHUNK_SIZE_EXT_CLIP;
	minClipEndSize = MIN_CLIP_END_SIZE;
	expected_cov_assemble = EXPECTED_COV_ASSEMBLE;
	num_threads_per_assem_work = NUM_THREADS_PER_ASSEM_WORK;
	outDir = OUT_DIR;

	static struct option lopts[] = {
//		{ "daemon", no_argument, NULL, 'D' },
//		{ "dir", required_argument, NULL, 'd' },
//		{ "out", required_argument, NULL, 'o' },
//		{ "log", required_argument, NULL, 'l' },
		{ "sample", required_argument, NULL, 0 },
		{ "threads-per-assem-work", required_argument, NULL, 0 },
		{ "assem-chunk-size", required_argument, NULL, 0 },
		{ "assem-side-ext-size", required_argument, NULL, 0 },
		{ "assem-side-ext-size-clip", required_argument, NULL, 0 },
		{ "min-cons-read-size", required_argument, NULL, 0 },
		{ "keep-assemble-reads", no_argument, NULL, 0 },
		{ "keep-failed-reads", no_argument, NULL, 0 },
		{ "reassemble-failed-work", no_argument, NULL, 0 },
		{ "min-input-cov-assemble", required_argument, NULL, 0 },
		{ "monitor_proc_names_assemble", required_argument, NULL, 0 },
		{ "max_proc_running_minutes_assemble", required_argument, NULL, 0 },
		{ "technology", required_argument, NULL, 0 },
		{ "include-decoy", no_argument, NULL, 0 },
		{ "version", no_argument, NULL, 'v' },
		{ "help", no_argument, NULL, 'h' },
		{ NULL, 0, NULL, 0 }
	};

	while( (opt = getopt_long(argc, argv, ":b:m:M:n:r:e:x:o:p:t:vh", lopts, &option_index)) != -1 ){
		switch(opt){
			case 'b': blockSize = stoi(optarg); break;
			//case 'S': assemSlideSize = stoi(optarg); break;
			case 'm': min_sv_size_usr = stoi(optarg); break;
			case 'M': max_sv_size_usr = stoi(optarg); break;
			case 'n': minReadsNumSupportSV = stoi(optarg); break;
			case 'r': max_seg_size_ratio_usr = stod(optarg); break;
			case 'e': minClipEndSize = stoi(optarg); break;
			case 'x': expected_cov_assemble = stod(optarg); break;
			case 'o': outDir = optarg; break;
			case 'p': outFilePrefix = optarg; break;
			case 't': threadNum_tmp = stoi(optarg); break;
			case 'v': showVersion(); exit(0);
			case 'h': showAssembleUsage(); exit(0);
			case '?':
				if(optopt) cout << "unknown option '-" << (char)optopt << "'" << endl;
				else{ // Bad long option
                    if (optind > 0 && strncmp(argv[optind-1], "--", 2) == 0)
                        cout << "unknown option '" << argv[optind-1] << "'" << endl;
				}
				exit(1);
			case ':':
				if(optopt) cout << "the option '-" << (char)optopt << "' needs a value" << endl;
				else{
                    if (optind > 0 && strncmp(argv[optind-1], "--", 2) == 0)
                        cout << "the option '" << argv[optind-1] << "' needs a value" << endl;
				}
				exit(1);
			case 0: // long options
				if(parse_long_opt(option_index, optarg, lopts)!=0){
					showAssembleUsage();
					return 1;
				}
				break;
			default:
				cout << "Error: please specify the correct options for 'assemble' command" << endl;
				showAssembleUsage();
				return 1;
		}
	}

	load_from_file_flag = true;
	if(threadNum_tmp==0) num_threads = sysconf(_SC_NPROCESSORS_ONLN);
	else num_threads = (threadNum_tmp>=sysconf(_SC_NPROCESSORS_ONLN)) ? sysconf(_SC_NPROCESSORS_ONLN) : threadNum_tmp;
	maxVarRegSize = max_sv_size_usr;

	if(num_threads*num_threads_per_assem_work>sysconf(_SC_NPROCESSORS_ONLN)){ // warning
		cout << "Warning: the user-specified total number of concurrent assemble work is " << num_threads << ", and the user-specified number of threads for each assemble work is " << num_threads_per_assem_work << ", which exceeds the total number of available processors on the machine (" << sysconf(_SC_NPROCESSORS_ONLN) << ")." << endl;
	}

	outDir = deleteTailPathChar(outDir);

	opt = argc - optind; // the reference file and BAM file on the command line
	if(opt>=2) {
		refFile = argv[optind];
		inBamFile = argv[optind+1];
	}else{
		cout << "Error: Please specify the reference file and coordinate-sorted BAM file." << endl << endl;
		showAssembleUsage();
		return 1;
	}

	return 0;
}

// parse the parameters for call command
int Paras::parseCallParas(int argc, char **argv){
	int opt, threadNum_tmp = 0, option_index;
	blockSize = BLOCKSIZE;
	slideSize = SLIDESIZE;
	//assemSlideSize = ASSEM_SLIDE_SIZE;
	min_sv_size_usr = MIN_SV_SIZE_USR;
	max_sv_size_usr = MAX_SV_SIZE_USR;
	minReadsNumSupportSV = MIN_SUPPORT_READS_NUM_EST;
	//maxVarRegSize = MAX_VAR_REG_SIZE;
	minClipEndSize = MIN_CLIP_END_SIZE;
	assemSideExtSize = ASM_CHUNK_SIZE_EXT_INDEL;
	assemSideExtSizeClip = ASM_CHUNK_SIZE_EXT_CLIP;
	outDir = OUT_DIR;

	static struct option lopts[] = {
//		{ "daemon", no_argument, NULL, 'D' },
//		{ "dir", required_argument, NULL, 'd' },
//		{ "out", required_argument, NULL, 'o' },
//		{ "log", required_argument, NULL, 'l' },
		{ "monitor_proc_names_call", required_argument, NULL, 0 },
		{ "max_proc_running_minutes_call", required_argument, NULL, 0 },
		{ "sample", required_argument, NULL, 0 },
		{ "include-decoy", no_argument, NULL, 0 },
		{ "gt_min_sig_size", required_argument, NULL, 0 },
		{ "gt_size_ratio_match", required_argument, NULL, 0 },
		{ "gt_min_alle_ratio", required_argument, NULL, 0 },
		{ "gt_max_alle_ratio", required_argument, NULL, 0 },
		{ "version", no_argument, NULL, 'v' },
		{ "help", no_argument, NULL, 'h' },
		{ NULL, 0, NULL, 0 }
	};

	while( (opt = getopt_long(argc, argv, ":b:m:M:n:e:o:p:t:vh", lopts, &option_index)) != -1 ){
		switch(opt){
			case 'b': blockSize = stoi(optarg); break;
			//case 'S': assemSlideSize = stoi(optarg); break;
			case 'm': min_sv_size_usr = stoi(optarg); break;
			case 'M': max_sv_size_usr = stoi(optarg); break;
			case 'n': minReadsNumSupportSV = stoi(optarg); break;
			case 'e': minClipEndSize = stoi(optarg); break;
			case 'o': outDir = optarg; break;
			case 'p': outFilePrefix = optarg; break;
			case 't': threadNum_tmp = stoi(optarg); break;
			case 'v': showVersion(); exit(0);
			case 'h': showCallUsage(); exit(0);
			case '?':
				if(optopt) cout << "unknown option '-" << (char)optopt << "'" << endl;
				else{ // Bad long option
                    if (optind > 0 && strncmp(argv[optind-1], "--", 2) == 0)
                        cout << "unknown option '" << argv[optind-1] << "'" << endl;
				}
				exit(1);
			case ':':
				if(optopt) cout << "the option -'" << (char)optopt << "' needs a value" << endl;
				else{
                    if (optind > 0 && strncmp(argv[optind-1], "--", 2) == 0)
                        cout << "the option '" << argv[optind-1] << "' needs a value" << endl;
				}
				exit(1);
			case 0: // long options
				if(parse_long_opt(option_index, optarg, lopts)!=0){
					showCallUsage();
					return 1;
				}
				break;
			default:
				cout << "Error: please specify the correct options for 'call' command" << endl;
				showCallUsage();
				return 1;
		}
	}

	load_from_file_flag = true;
	if(threadNum_tmp==0) num_threads = sysconf(_SC_NPROCESSORS_ONLN);
	else num_threads = (threadNum_tmp>=sysconf(_SC_NPROCESSORS_ONLN)) ? sysconf(_SC_NPROCESSORS_ONLN) : threadNum_tmp;
	maxVarRegSize = max_sv_size_usr;

	outDir = deleteTailPathChar(outDir);

	opt = argc - optind; // the reference file and BAM file on the command line
	if(opt>=2) {
		refFile = argv[optind];
		inBamFile = argv[optind+1];
	}else{
		cout << "Error: Please specify the reference file and coordinate-sorted BAM file." << endl << endl;
		showCallUsage();
		return 1;
	}

	return 0;
}

// parse the parameters for 'all' command
int Paras::parseAllParas(int argc, char **argv, const string &cmd_str){
	int opt, threadNum_tmp = 0, option_index;
	blockSize = BLOCKSIZE;
	slideSize = SLIDESIZE;
	//assemSlideSize = ASSEM_SLIDE_SIZE;
	min_sv_size_usr = MIN_SV_SIZE_USR;
	max_sv_size_usr = MAX_SV_SIZE_USR;
	minReadsNumSupportSV = MIN_SUPPORT_READS_NUM_EST;
	max_seg_size_ratio_usr = MAX_SEG_SIZE_RATIO;
	//maxVarRegSize = MAX_VAR_REG_SIZE;
	assemChunkSize = ASM_CHUNK_SIZE_INDEL;
	assemSideExtSize = ASM_CHUNK_SIZE_EXT_INDEL;
	assemSideExtSizeClip = ASM_CHUNK_SIZE_EXT_CLIP;
	minClipEndSize = MIN_CLIP_END_SIZE;
	expected_cov_assemble = EXPECTED_COV_ASSEMBLE;
	num_threads_per_assem_work = NUM_THREADS_PER_ASSEM_WORK;
	outDir = OUT_DIR;
	simpleReg_t *simple_reg;
	string simple_reg_str;

	static struct option lopts[] = {
//		{ "daemon", no_argument, NULL, 'D' },
//		{ "dir", required_argument, NULL, 'd' },
//		{ "out", required_argument, NULL, 'o' },
//		{ "log", required_argument, NULL, 'l' },
		{ "sample", required_argument, NULL, 0 },
		{ "threads-per-assem-work", required_argument, NULL, 0 },
		{ "assem-chunk-size", required_argument, NULL, 0 },
		{ "keep-assemble-reads", no_argument, NULL, 0 },
		{ "assem-side-ext-size", required_argument, NULL, 0 },
		{ "assem-side-ext-size-clip", required_argument, NULL, 0 },
		{ "keep-failed-reads", no_argument, NULL, 0 },
		{ "reassemble-failed-work", no_argument, NULL, 0 },
		{ "mask-noisy-region", required_argument, NULL, 0 },
		{ "monitor_proc_names_assemble", required_argument, NULL, 0 },
		{ "monitor_proc_names_call", required_argument, NULL, 0 },
		{ "max_proc_running_minutes_assemble", required_argument, NULL, 0 },
		{ "max_proc_running_minutes_call", required_argument, NULL, 0 },
		{ "technology", required_argument, NULL, 0 },
		{ "include-decoy", no_argument, NULL, 0 },
		{ "gt_min_sig_size", required_argument, NULL, 0 },
		{ "gt_size_ratio_match", required_argument, NULL, 0 },
		{ "gt_min_alle_ratio", required_argument, NULL, 0 },
		{ "gt_max_alle_ratio", required_argument, NULL, 0 },
		{ "version", no_argument, NULL, 'v' },
		{ "help", no_argument, NULL, 'h' },
		{ NULL, 0, NULL, 0 }
	};

	while( (opt = getopt_long(argc, argv, ":b:s:m:M:n:r:e:x:o:p:t:vh", lopts, &option_index)) != -1 ){
		switch(opt){
			case 'b': blockSize = stoi(optarg); break;
			case 's': slideSize = stoi(optarg); break;
			//case 'S': assemSlideSize = stoi(optarg); break;
			case 'm': min_sv_size_usr = stoi(optarg); break;
			case 'M': max_sv_size_usr = stoi(optarg); break;
			case 'n': minReadsNumSupportSV = stoi(optarg); break;
			case 'r': max_seg_size_ratio_usr = stod(optarg); break;
			case 'e': minClipEndSize = stoi(optarg); break;
			case 'x': expected_cov_assemble = stod(optarg); break;
			case 'o': outDir = optarg; break;
			case 'p': outFilePrefix = optarg; break;
			case 't': threadNum_tmp = stoi(optarg); break;
			case 'v': showVersion(); exit(0);
			case 'h': showAllUsage(cmd_str); exit(0);
			case '?':
				if(optopt) cout << "unknown option '-" << (char)optopt << "'" << endl;
				else{ // Bad long option
                    if (optind > 0 && strncmp(argv[optind-1], "--", 2) == 0) {
                        cout << "unknown option '" << argv[optind-1] << "'" << endl;
                    }
				}
				exit(1);
			case ':':
				if(optopt) cout << "the option '-" << (char)optopt << "' needs a value" << endl;
				else{
                    if (optind > 0 && strncmp(argv[optind-1], "--", 2) == 0)
                        cout << "the option '" << argv[optind-1] << "' needs a value" << endl;
				}
				exit(1);
			case 0: // long options
				if(parse_long_opt(option_index, optarg, lopts)!=0){
					showAllUsage(cmd_str);
					return 1;
				}
				break;
			default:
				cout << "Error: please specify the correct options for 'all' command" << endl;
				showAllUsage(cmd_str);
				return 1;
		}
	}

	load_from_file_flag = false;
	if(threadNum_tmp==0) num_threads = sysconf(_SC_NPROCESSORS_ONLN);
	else num_threads = (threadNum_tmp>=sysconf(_SC_NPROCESSORS_ONLN)) ? sysconf(_SC_NPROCESSORS_ONLN) : threadNum_tmp;
	maxVarRegSize = max_sv_size_usr;

	if(num_threads*num_threads_per_assem_work>sysconf(_SC_NPROCESSORS_ONLN)){ // warning
		cout << "Warning: the user-specified total number of concurrent assemble work is " << num_threads << ", and the user-specified number of threads for each assemble work is " << num_threads_per_assem_work << ", which exceeds the total number of available processors on the machine (" << sysconf(_SC_NPROCESSORS_ONLN) << ")." << endl;
	}

	outDir = deleteTailPathChar(outDir);

	opt = argc - optind; // the reference file and BAM file on the command line
	if(opt>=2) {
		refFile = argv[optind];
		inBamFile = argv[optind+1];

		for(int i=optind+2; i<argc; i++){
			simple_reg_str = argv[i];
			simple_reg = allocateSimpleReg(simple_reg_str);
			if(simple_reg) limit_reg_vec.push_back(simple_reg);
			else{
				cout << "Error: Please specify the correct genomic regions to be processed." << endl << endl;
				showAllUsage(cmd_str);
				return 1;
			}
		}
		if(limit_reg_vec.size()) limit_reg_process_flag = true;
	}else{
		cout << "Error: Please specify the reference file and coordinate-sorted BAM file." << endl << endl;
		showAllUsage(cmd_str);
		return 1;
	}

	if(refFile.size()==0){
		cout << "Error: Please specify the reference" << endl << endl;
		showAllUsage(cmd_str);
		return 1;
	}

	return 0;
}

// parse the parameters for 'detect-assemble' command
int Paras::parseDetectAssembleParas(int argc, char **argv){
	parseAllParas(argc, argv, "detect-assemble");
	return 0;
}

// show the usage
void Paras::showUsage(){
	cout << "Program: " << PROG_NAME << " (" << PROG_DESC << ")" << endl;
	cout << "Version: " << PROG_VERSION << " (using htslib " << hts_version() << ")" << endl << endl;
	cout << "Usage:  asvclr  <command> [options] <REF_FILE> <BAM_FILE> [Region ...]" << endl << endl;

	cout << "Description:" << endl;
	cout << "   REF_FILE          Reference file (required)" << endl;
	cout << "   BAM_FILE          Coordinate-sorted BAM file (required)" << endl;
	cout << "   Region            Reference regions to process: CHR|CHR:START-END." << endl;
	cout << "                     If unspecified, all reference regions will be " << endl;
	cout << "                     processed (optional)" << endl << endl;

	cout << "Commands:" << endl;
	cout << "   detect            detect indel signatures in aligned reads" << endl;
	cout << "   assemble          assemble candidate regions" << endl;
	cout << "   call              call indels by alignments of local genome assemblies" << endl;
	cout << "   all               run the above commands in turn" << endl;
	cout << "   detect-assemble   run 'detect' and 'assemble' commands in turn" << endl;
}

// show the usage for detect command
void Paras::showDetectUsage(){
	cout << "Program: " << PROG_NAME << " (" << PROG_DESC << ")" << endl;
	cout << "Version: " << PROG_VERSION << " (using htslib " << hts_version() << ")" << endl << endl;
	cout << "Usage: asvclr detect [options] <REF_FILE> <BAM_FILE> [Region ...]" << endl << endl;

	cout << "Description:" << endl;
	cout << "   REF_FILE      Reference file (required)" << endl;
	cout << "   BAM_FILE      Coordinate-sorted BAM file (required)" << endl;
	cout << "   Region        Reference regions to process: CHR|CHR:START-END." << endl;
	cout << "                 If unspecified, all reference regions will be " << endl;
	cout << "                 processed (optional)" << endl << endl;

	cout << "Options: " << endl;
	cout << "   -b INT        block size [1000000]" << endl;
	cout << "   -s INT        Slide window size [" << SLIDESIZE <<"]" << endl;
	cout << "   -m INT        minimal SV size to report [" << MIN_SV_SIZE_USR << "]" << endl;
	cout << "   -M INT        maximal SV size to report [" << MAX_SV_SIZE_USR << "]" << endl;
	cout << "                 Variants with size smaller than threshold will be ignored" << endl;
	cout << "   -n INT        minimal number of reads supporting a SV [" << MIN_SUPPORT_READS_NUM_EST << "]." << endl;
	cout << "                 When it is not specified, " << MIN_SUPPORT_READS_NUM_EST << " means the value will be" << endl;
	cout << "                 estimated to be "<< MIN_SUPPORT_READS_NUM_FACTOR << " times of the sequencing depth of" << endl;
	cout << "                 the data set (rounded)" << endl;
	cout << "   -e INT        minimal clipping end size [" << MIN_CLIP_END_SIZE << "]. Clipping events" << endl;
	cout << "                 with size smaller than threshold will be ignored" << endl;
	//cout << "   -r FILE       limit reference regions to process [null]: CHR|CHR:START-END" << endl;
	cout << "   -o DIR        output directory [output]" << endl;
	cout << "   -p STR        prefix of output result files [null]" << endl;
	cout << "   -t INT        number of concurrent work [0]. 0 for the maximal number" << endl;
	cout << "                 of threads in machine" << endl;
	cout << "   --include-decoy" << endl;
	cout << "                 include decoy items in result" << endl;
	cout << "   --sample STR  Sample name [\"" << SAMPLE_DEFAULT << "\"]" << endl;

	cout << "   -v,--version  show version information" << endl;
	cout << "   -h,--help     show this help message and exit" << endl;
}

// show the usage for assemble command
void Paras::showAssembleUsage(){
	cout << "Program: " << PROG_NAME << " (" << PROG_DESC << ")" << endl;
	cout << "Version: " << PROG_VERSION << " (using htslib " << hts_version() << ")" << endl << endl;
	cout << "Usage: asvclr assemble [options] <REF_FILE> <BAM_FILE>" << endl << endl;

	cout << "Description:" << endl;
	cout << "   REF_FILE      Reference file (required)" << endl;
	cout << "   BAM_FILE      Coordinate-sorted BAM file (required)" << endl << endl;

	cout << "Options: " << endl;
	//cout << "   -s INT        Slide window size for 'detect' command  [" << SLIDESIZE <<"]" << endl;
	//cout << "   -S INT        Slide window size for 'assemble' and 'call' command [" << ASSEM_SLIDE_SIZE << "]" << endl;
	cout << "   -e INT        minimal clipping end size [" << MIN_CLIP_END_SIZE << "]. Clipping events" << endl;
	cout << "                 with size smaller than threshold will be ignored" << endl;
	cout << "   -x FLOAT      expected sampling coverage for local assemble [" << EXPECTED_COV_ASSEMBLE << "], " << endl;
	cout << "                 0 for no coverage sampling" << endl;
	cout << "   -o DIR        output directory [output]" << endl;
	cout << "   -p STR        prefix of output result files [null]" << endl;
	cout << "   -t INT        number of concurrent work [0]. 0 for the maximal number" << endl;
	cout << "                 of threads in machine" << endl;
	cout << "   -r FLOAT      minimal ratio threshold of the largest split-alignment segment" << endl;
	cout << "                 of a read allowing for indel detection. [" << MAX_SEG_SIZE_RATIO << "]" << endl;
	cout << "   --threads-per-assem-work INT" << endl;
	cout << "                 Limited number of threads for each assemble work [0]:" << endl;
	cout << "                 0 for unlimited, and positive INT for the limited" << endl;
	cout << "                 number of threads for each assemble work" << endl;
	cout << "   --assem-chunk-size INT" << endl;
	cout << "                 maximal reference chunk size to collect reads data to perform" << endl;
	cout << "                 local assemble [" << ASM_CHUNK_SIZE_INDEL << "]. Reads of variants with reference" << endl;
	cout << "                 distance < INT will be collected to perform local assemble" << endl;
	cout << "   --assem-side-ext-size INT" << endl;
	cout << "                 region extend size on both sides while generating consensus sequences" << endl;
	cout << "                 [" << ASM_CHUNK_SIZE_EXT_INDEL << "] bp." << endl;
	cout << "   --assem-side-ext-size-clip INT" << endl;
	cout << "                 clip region extend size on both sides while generating assembly sequences" << endl;
	cout << "                 [" << ASM_CHUNK_SIZE_EXT_CLIP << "] bp." << endl;
	cout << "   --min-cons-read-size INT" << endl;
	cout << "                 minimal read segment size to generate consensus sequences " << "[" << MIN_CONS_READ_LEN << "] bp." << endl;
	cout << "   --keep-assemble-reads" << endl;
	cout << "                 Keep temporary reads from being deleted during local assemble." << endl;
	cout << "                 This may take some additional disk space" << endl;
	cout << "   --reassemble-failed-work" << endl;
	cout << "                 Reperform previously failed local assemble work." << endl;
	cout << "   --min-input-cov-assemble FLOAT" << endl;
	cout << "                 Minimum input coverage for local assemble [" << MIN_INPUT_COV_CANU << "]" << endl;
	cout << "   --monitor_proc_names_assemble STR" << endl;
	cout << "                 Process names to be monitored during Canu assemble. These processes may" << endl;
	cout << "                 have ultra-high CPU running time under some certain circumstances and" << endl;
	cout << "                 should be terminated in advance if they are computation intensive works." << endl;
	cout << "                 Note that the process names should be comma-delimited and without blanks:" << endl;
	cout << "                 [\"" << DEFAULT_MONITOR_PROC_NAMES_ASSEMBLE << "\"]" << endl;
	cout << "   --max_proc_running_minutes_assemble INT" << endl;
	cout << "                 Monitored processes for assemble will be terminated if their CPU running" << endl;
	cout << "                 time exceed INT minutes: [" << MAX_PROC_RUNNING_MINUTES_ASSEMBLE << "]" << endl;
	cout << "   --technology STR" << endl;
	cout << "                 Sequencing technology [pacbio]:" << endl;
	cout << "                   pacbio     : the PacBio CLR sequencing technology;" << endl;
	cout << "                   nanopore   : the Nanopore sequencing technology;" << endl;
	cout << "                   pacbio-hifi: the PacBio CCS sequencing technology." << endl;
	cout << "   --include-decoy" << endl;
	cout << "                 include decoy items in result" << endl;
	cout << "   --sample STR  Sample name [\"" << SAMPLE_DEFAULT << "\"]" << endl;
	cout << "   -v,--version  show version information" << endl;
	cout << "   -h,--help     show this help message and exit" << endl;
}

// show the usage for call command
void Paras::showCallUsage(){
	cout << "Program: " << PROG_NAME << " (" << PROG_DESC << ")" << endl;
	cout << "Version: " << PROG_VERSION << " (using htslib " << hts_version() << ")" << endl << endl;
	cout << "Usage: asvclr call [options] <REF_FILE> <BAM_FILE>" << endl << endl;

	cout << "Description:" << endl;
	cout << "   REF_FILE      Reference file (required)" << endl;
	cout << "   BAM_FILE      Coordinate-sorted BAM file (required)" << endl << endl;

	cout << "Options: " << endl;
	//cout << "   -s INT        Slide window size for 'detect' command  [" << SLIDESIZE <<"]" << endl;
	//cout << "   -S INT        Slide window size for 'assemble' and 'call' command [" << ASSEM_SLIDE_SIZE << "]" << endl;
	cout << "   -m INT        minimal SV size to report [" << MIN_SV_SIZE_USR << "]" << endl;
	cout << "   -M INT        maximal SV size to report [" << MAX_SV_SIZE_USR << "]" << endl;
	cout << "                 Variants with size smaller than threshold will be ignored" << endl;
	cout << "   -e INT        minimal clipping end size [" << MIN_CLIP_END_SIZE << "]. Clipping events" << endl;
	cout << "                 with size smaller than threshold will be ignored" << endl;
	cout << "   -o DIR        output directory [output]" << endl;
	cout << "   -p STR        prefix of output result files [null]" << endl;
	cout << "   -t INT        number of concurrent work [0]. 0 for the maximal number" << endl;
	cout << "                 of threads in machine" << endl;
	cout << "   --monitor_proc_names_call STR" << endl;
	cout << "                 Process names to be monitored during BLAT alignment. These processes may" << endl;
	cout << "                 have ultra-high CPU running time under some certain circumstances and" << endl;
	cout << "                 should be terminated in advance if they are computation intensive works." << endl;
	cout << "                 Note that the process names should be comma-delimited and without blanks:" << endl;
	cout << "                 [\"" << DEFAULT_MONITOR_PROC_NAMES_CALL << "\"]" << endl;
	cout << "   --max_proc_running_minutes_call INT" << endl;
	cout << "                 Monitored processes for call will be terminated if their CPU running time" << endl;
	cout << "                 exceed INT minutes: [" << MAX_PROC_RUNNING_MINUTES_CALL << "]" << endl;
	cout << "   --include-decoy" << endl;
	cout << "                 include decoy items in result" << endl;
	cout << "   --sample STR  Sample name [\"" << SAMPLE_DEFAULT << "\"]" << endl;

	cout << "   --gt_min_sig_size INT" << endl;
	cout << "                 minimal signature size threshold for genotyping [" << GT_SIG_SIZE_THRES << "]." << endl;
	cout << "                 Signatures with size larger than INT will be processed, otherwise, they will be ignored." << endl;
	cout << "   --gt_size_ratio_match FLOAT" << endl;
	cout << "                 signature size ratio threshold for genotyping [" << GT_SIZE_RATIO_MATCH_THRES << "]." << endl;
	cout << "                 Two signatures are match if the ratio of their sizes is larger than FLOAT." << endl;
	cout << "   --gt_min_alle_ratio FLOAT" << endl;
	cout << "                 minimal allele ratio threshold for genotyping [" << GT_MIN_ALLE_RATIO_THRES << "]." << endl;
	cout << "                 Variation is homozygous if the ratio of allele count is less than FLOAT. " << endl;
	cout << "                 Variation is heterozygous if the ratio of allele count is larger than FLOAT" << endl;
	cout << "                 and less than maximal allele ratio threshold (--gt_max_alle_ratio)." << endl;
	cout << "   --gt_max_alle_ratio FLOAT" << endl;
	cout << "                 maximal allele ratio threshold for genotyping [" << GT_MAX_ALLE_RATIO_THRES << "]." << endl;
	cout << "                 There is no genotype information if ratio of allele count is larger than FLOAT." << endl;
	cout << "                 Variation is heterozygous if the ratio of allele count is less than FLOAT" << endl;
	cout << "                 and larger than minimal allele ratio threshold (--gt_min_alle_ratio)." << endl;

	cout << "   -v,--version  show version information" << endl;
	cout << "   -h,--help     show this help message and exit" << endl;
}

// show the usage for all command
void Paras::showAllUsage(const string &cmd_str){
	cout << "Program: " << PROG_NAME << " (" << PROG_DESC << ")" << endl;
	cout << "Version: " << PROG_VERSION << " (using htslib " << hts_version() << ")" << endl << endl;
	cout << "Usage: asvclr " << cmd_str << " [options] <REF_FILE> <BAM_FILE> [Region ...]" << endl << endl;

	cout << "Description:" << endl;
	cout << "   REF_FILE      Reference file (required)" << endl;
	cout << "   BAM_FILE      Coordinate-sorted file (required)" << endl;
	cout << "   Region        Reference regions to process: CHR|CHR:START-END." << endl;
	cout << "                 If unspecified, all reference regions will be " << endl;
	cout << "                 processed (optional)" << endl << endl;

	cout << "Options: " << endl;
	cout << "   -b INT        block size [1000000]" << endl;
	cout << "   -s INT        Slide window size [" << SLIDESIZE <<"]" << endl;
	//cout << "   -S INT        Slide window size for 'assemble' and 'call' command [" << ASSEM_SLIDE_SIZE << "]" << endl;
	cout << "   -m INT        minimal SV size to report [" << MIN_SV_SIZE_USR << "]" << endl;
	cout << "   -M INT        maximal SV size to report [" << MAX_SV_SIZE_USR << "]" << endl;
	cout << "                 Variants with size smaller than threshold will be ignored" << endl;
	//cout << "   -n INT        minimal number of reads supporting a SV [" << MIN_SUPPORT_READS_NUM_EST << "]" << endl;
	cout << "   -n INT        minimal number of reads supporting a SV [" << MIN_SUPPORT_READS_NUM_EST << "]." << endl;
	cout << "                 When it is not specified, " << MIN_SUPPORT_READS_NUM_EST << " means the value will be" << endl;
	cout << "                 estimated to be "<< MIN_SUPPORT_READS_NUM_FACTOR << " times of the sequencing depth of" << endl;
	cout << "                 the data set (rounded)" << endl;
	cout << "   -r FLOAT      minimal ratio threshold of the largest split-alignment segment" << endl;
	cout << "                 of a read allowing for indel detection. [" << MAX_SEG_SIZE_RATIO << "]" << endl;
	cout << "   -e INT        minimal clipping end size [" << MIN_CLIP_END_SIZE << "]. Clipping events" << endl;
	cout << "                 with size smaller than threshold will be ignored" << endl;
	cout << "   -x FLOAT      expected sampling coverage for local assemble [" << EXPECTED_COV_ASSEMBLE << "], " << endl;
	cout << "                 0 for no coverage sampling" << endl;
	cout << "   -o DIR        output directory [output]" << endl;
	cout << "   -p STR        prefix of output result files [null]" << endl;
	cout << "   -t INT        number of concurrent work [0]. 0 for the maximal number" << endl;
	cout << "                 of threads in machine" << endl;
	cout << "   --threads-per-assem-work INT" << endl;
	cout << "                 Limited number of threads for each assemble work [0]:" << endl;
	cout << "                 0 for unlimited, and positive INT for the limited" << endl;
	cout << "                 number of threads for each assemble work" << endl;
	cout << "   --assem-chunk-size INT" << endl;
	cout << "                 maximal reference chunk size to collect reads data to perform" << endl;
	cout << "                 local assemble [" << ASM_CHUNK_SIZE_INDEL << "]. Reads of variants with reference" << endl;
	cout << "                 distance < INT will be collected to perform local assemble" << endl;
	cout << "   --assem-side-ext-size INT" << endl;
	cout << "                 region extend size on both sides while generating consensus sequences" << endl;
	cout << "                 [" << ASM_CHUNK_SIZE_EXT_INDEL << "] bp." << endl;
	cout << "   --assem-side-ext-size-clip INT" << endl;
	cout << "                 clip region extend size on both sides while generating assembly sequences" << endl;
	cout << "                 [" << ASM_CHUNK_SIZE_EXT_CLIP << "] bp." << endl;
	cout << "   --min-cons-read-size INT" << endl;
	cout << "                 minimal read segment size to generate consensus sequences " << "[" << MIN_CONS_READ_LEN << "] bp." << endl;
	cout << "   --keep-assemble-reads" << endl;
	cout << "                 Keep temporary reads from being deleted during local assemble." << endl;
	cout << "                 This may take some additional disk space" << endl;
	cout << "   --reassemble-failed-work" << endl;
	cout << "                 Reperform previously failed local assemble work." << endl;
	cout << "   --min-input-cov-assemble FLOAT" << endl;
	cout << "                 Minimum input coverage for local assemble [" << MIN_INPUT_COV_CANU << "]" << endl;
	cout << "   --monitor_proc_names_assemble STR" << endl;
	cout << "                 Process names to be monitored during Canu assemble. These processes may" << endl;
	cout << "                 have ultra-high CPU running time under some certain circumstances and" << endl;
	cout << "                 should be terminated in advance if they are computation intensive works." << endl;
	cout << "                 Note that the process names should be comma-delimited and without blanks:" << endl;
	cout << "                 [\"" << DEFAULT_MONITOR_PROC_NAMES_ASSEMBLE << "\"]" << endl;

	if(cmd_str.compare("all")==0){
		cout << "   --monitor_proc_names_call STR" << endl;
		cout << "                 Process names to be monitored during BLAT alignment. These processes may" << endl;
		cout << "                 have ultra-high CPU running time under some certain circumstances and" << endl;
		cout << "                 should be terminated in advance if they are computation intensive works." << endl;
		cout << "                 Note that the process names should be comma-delimited and without blanks:" << endl;
		cout << "                 [\"" << DEFAULT_MONITOR_PROC_NAMES_CALL << "\"]" << endl;
	}

	cout << "   --max_proc_running_minutes_assemble INT" << endl;
	cout << "                 Monitored processes for assemble will be terminated if their CPU running" << endl;
	cout << "                 time exceed INT minutes: [" << MAX_PROC_RUNNING_MINUTES_ASSEMBLE << "]" << endl;

	if(cmd_str.compare("all")==0){
		cout << "   --max_proc_running_minutes_call INT" << endl;
		cout << "                 Monitored processes for call will be terminated if their CPU running time" << endl;
		cout << "                 exceed INT minutes: [" << MAX_PROC_RUNNING_MINUTES_CALL << "]" << endl;
	}

	cout << "   --technology STR" << endl;
	cout << "                 Sequencing technology [pacbio]:" << endl;
	cout << "                   pacbio     : the PacBio CLR sequencing technology;" << endl;
	cout << "                   nanopore   : the Nanopore sequencing technology;" << endl;
	cout << "                   pacbio-hifi: the PacBio CCS sequencing technology." << endl;
	cout << "   --include-decoy" << endl;
	cout << "                 include decoy items in result" << endl;
	cout << "   --sample STR  Sample name [\"" << SAMPLE_DEFAULT << "\"]" << endl;

	cout << "   --gt_min_sig_size INT" << endl;
	cout << "                 minimal signature size threshold for genotyping [" << GT_SIG_SIZE_THRES << "]." << endl;
	cout << "                 Signatures with size larger than INT will be processed, otherwise, they will be ignored." << endl;
	cout << "   --gt_size_ratio_match FLOAT" << endl;
	cout << "                 signature size ratio threshold for genotyping [" << GT_SIZE_RATIO_MATCH_THRES << "]." << endl;
	cout << "                 Two signatures are match if the ratio of their sizes is larger than FLOAT." << endl;
	cout << "   --gt_min_alle_ratio FLOAT" << endl;
	cout << "                 minimal allele ratio threshold for genotyping [" << GT_MIN_ALLE_RATIO_THRES << "]." << endl;
	cout << "                 Variation is homozygous if the ratio of allele count is less than FLOAT. " << endl;
	cout << "                 Variation is heterozygous if the ratio of allele count is larger than FLOAT" << endl;
	cout << "                 and less than maximal allele ratio threshold (--gt_max_alle_ratio)." << endl;
	cout << "   --gt_max_alle_ratio FLOAT" << endl;
	cout << "                 maximal allele ratio threshold for genotyping [" << GT_MAX_ALLE_RATIO_THRES << "]." << endl;
	cout << "                 There is no genotype information if ratio of allele count is larger than FLOAT." << endl;
	cout << "                 Variation is heterozygous if the ratio of allele count is less than FLOAT" << endl;
	cout << "                 and larger than minimal allele ratio threshold (--gt_min_alle_ratio)." << endl;

	cout << "   -v,--version  show version information" << endl;
	cout << "   -h,--help     show this help message and exit" << endl;
}

// show the usage for detect-assemble command
void Paras::showDetectAssembleUsage(){
	showAllUsage("detect-assemble");
}

// output parameters
void Paras::outputParas(){
	char buf[256];
	string tmp;

	// check the paths: current path, outDir
	getcwd(buf, 256);
	tmp = buf;
	if(tmp.find(' ')!=string::npos){
		cout << "Error: The current work directory '" << tmp << "' should not contain space characters, please specify a correct one." << endl;
		exit(1);
	}
	if(outDir.find(' ')!=string::npos){
		cout << "Error: The output directory '" << outDir << "' should not contain space characters, please specify a correct one." << endl;
		exit(1);
	}

	cout << "Program: " << PROG_NAME << " (" << PROG_DESC << ")" << endl;
	cout << "Version: " << PROG_VERSION << " (using htslib " << hts_version() << ")" << endl << endl;

	if(refFile.size()) cout << "Reference file: " << refFile << endl;
	if(inBamFile.size()) cout << "Alignment file: " << inBamFile << endl;
	if(outDir.size()) cout << "Output directory: " << outDir << endl;
	if(outFilePrefix.size()) cout << "Output result file prefix: " << outFilePrefix << endl;

	cout << "Sample: " << sample << endl;
	if(minReadsNumSupportSV != -1){
		cout << "Number of reads supporting SV: " << minReadsNumSupportSV << endl;
	}
	cout << "Block size: " << blockSize << " bp" << endl;
	cout << "Slide size: " << slideSize << " bp" << endl;
	cout << "Minimal SV size to report: " << min_sv_size_usr << " bp" << endl;
	cout << "Maximal SV size to report: " << max_sv_size_usr << " bp" << endl;
	cout << "Minimal ratio of max-segment for indel calling: " << max_seg_size_ratio_usr << endl;
	cout << "Minimal clipping end size: " << minClipEndSize << " bp" << endl;
	cout << "Expected sampling coverage: " << expected_cov_assemble << endl;
	cout << "Local assemble chunk size : " << assemChunkSize << " bp" << endl;
	cout << "Local assemble chunk extend size for indels: " << assemSideExtSize << " bp" << endl;
	cout << "Local assemble chunk extend size for clippings: " << assemSideExtSizeClip << " bp" << endl;
	cout << "Minimal read segment size to generate consensus sequence: " << minConReadLen << " bp" << endl;
	cout << "Number of concurrent works: " << num_threads << endl;
	cout << "Limited number of threads for each assemble work: " << num_threads_per_assem_work << endl;
	if(maskMisAlnRegFlag) cout << "Mask noisy regions: yes" << endl;
	if(delete_reads_flag==false) cout << "Retain local temporary reads: yes" << endl;
	if(keep_failed_reads_flag) cout << "Retain failed local temporary reads: yes" << endl;
	if(reassemble_failed_work_flag) cout << "Reperform previously failed local assemble work: yes" << endl;
	cout << "Minimum input coverage for local assemble: " << min_input_cov_canu << endl;
	if(command.compare("assemble")==0 or command.compare("all")==0 or command.compare("detect-assemble")==0)
		cout << "Monitored process names for assemble: " << monitoring_proc_names_assemble << endl;
	if(command.compare("call")==0 or command.compare("all")==0)
		cout << "Monitored process names for call: " << monitoring_proc_names_call << endl;
	if(command.compare("assemble")==0 or command.compare("all")==0 or command.compare("detect-assemble")==0)
		cout << "Maximum monitored process running minutes for assemble: " << max_proc_running_minutes_assemble << endl;
	if(command.compare("call")==0 or command.compare("all")==0)
		cout << "Maximum monitored process running minutes for call: " << max_proc_running_minutes_call << endl;
	if(include_decoy) cout << "Include decoy items: yes" << endl;
	cout << "Sequencing technology: " << technology << endl;
	cout << "Canu version: " << canu_version << endl << endl;

	bool recommend_ver_flag = isRecommendCanuVersion(canu_version, CANU_RECOMMEND_VERSION);
	if(recommend_ver_flag==false)
		cout << "Warning: detected the installed canu version is " << canu_version << ", however, it is highly recommended to use the latest version or the version " << CANU_RECOMMEND_VERSION << " and higher." << endl;

	string desc = "Regions to process:";
	printLimitRegs(limit_reg_vec, desc);
}

// output the estimation parameters
void Paras::outputEstParas(string info){
	cout << info << endl;
	cout << "Mean read length: " << mean_read_len << " bp" << endl;
	if(minReadsNumSupportSV != -1){
		cout << "Number of reads supporting SV: " << minReadsNumSupportSV << endl;
	}

	cout << "min_ins_size_filt: " << min_ins_size_filt << " bp" << endl;
	cout << "min_del_size_filt: " << min_del_size_filt << " bp" << endl;
//	cout << "min_clip_size_filt: " << min_clip_size_filt << " bp" << endl;
	cout << "min_ins_num_filt: " << min_ins_num_filt << endl;
	cout << "min_del_num_filt: " << min_del_num_filt << endl;
	cout << "min_clip_num_filt: " << min_clip_num_filt << endl;
	cout << "large_indel_size_thres: " << large_indel_size_thres << " bp" << endl << endl;
}

// initialize the estimation auxiliary data
void Paras::initEst(){
	int32_t i;
	for(i=0; i<AUX_ARR_SIZE; i++) insSizeEstArr[i] = 0;
	for(i=0; i<AUX_ARR_SIZE; i++) delSizeEstArr[i] = 0;
	for(i=0; i<AUX_ARR_SIZE; i++) clipSizeEstArr[i] = 0;
	for(i=0; i<AUX_ARR_SIZE; i++) insNumEstArr[i] = 0;
	for(i=0; i<AUX_ARR_SIZE; i++) delNumEstArr[i] = 0;
	for(i=0; i<AUX_ARR_SIZE; i++) clipNumEstArr[i] = 0;
}

// estimate parameters
void Paras::estimate(size_t op_est){
	if(op_est==SIZE_EST_OP){
		// min_ins_size_filt
//		cout << "min_ins_size_filt:" << endl;
		min_ins_size_filt = estimateSinglePara(insSizeEstArr, AUX_ARR_SIZE, SIZE_PERCENTILE_EST, MIN_INDEL_EVENT_SIZE);

		// min_del_size_filt
//		cout << "min_del_size_filt:" << endl;
		min_del_size_filt = estimateSinglePara(delSizeEstArr, AUX_ARR_SIZE, SIZE_PERCENTILE_EST, MIN_INDEL_EVENT_SIZE);

		// min_clip_size_filt
//		min_clip_size_filt = estimateSinglePara(clipSizeEstArr, AUX_ARR_SIZE, SIZE_PERCENTILE_EST, MIN_INDEL_EVENT_SIZE);

		if(min_ins_size_filt>min_del_size_filt) large_indel_size_thres = min_ins_size_filt * LARGE_INDEL_SIZE_FACTOR;
		else large_indel_size_thres = min_del_size_filt * LARGE_INDEL_SIZE_FACTOR;

		// compute mean read length
		if(total_read_num_est>0) mean_read_len = round((double)mean_read_len / total_read_num_est);
		else mean_read_len = total_read_num_est = 0;

		//compute mean read depth
		//if(minReadsNumSupportSV==-1) minReadsNumSupportSV =  round(((double)total_depth/chrome_num)*0.1);
		if(minReadsNumSupportSV==-1) minReadsNumSupportSV = estimateMinReadsNumSupportSV(mean_depth_vec);


	}else if(op_est==NUM_EST_OP){
		// min_ins_num_filt
//		cout << "min_ins_num_filt:" << endl;
		min_ins_num_filt = estimateSinglePara(insNumEstArr, AUX_ARR_SIZE, NUM_PERCENTILE_EST, MIN_INDEL_EVENT_NUM);

		// min_del_num_filt
//		cout << "min_del_num_filt:" << endl;
		min_del_num_filt = estimateSinglePara(delNumEstArr, AUX_ARR_SIZE, NUM_PERCENTILE_EST, MIN_INDEL_EVENT_NUM);

		// min_clip_num_filt
//		cout << "min_clip_num_filt:" << endl;
		min_clip_num_filt = estimateSinglePara(clipNumEstArr, AUX_ARR_SIZE, NUM_PERCENTILE_EST, MIN_INDEL_EVENT_NUM);
	}else if(op_est==SNV_EST_OP){

	}else{

		cerr << __func__ << ", line=" << __LINE__ << ", invalid estimation op_flag: " << op_est << endl;
		exit(1);
	}
}

// estimate single parameter
int64_t Paras::estimateSinglePara(int64_t *arr, int32_t n, double threshold, int32_t min_val){
	int64_t i, total, val = 1;
	double total1;

	for(i=0, total=0; i<AUX_ARR_SIZE; i++) total += arr[i];
	if(total>0)
//		cout << "total=" << total << endl;
		for(i=0, total1=0; i<AUX_ARR_SIZE; i++){
			total1 += arr[i];
//			cout << "\t" << i << ": " << arr[i] << ", " << total1/total << endl;
			if(total1/total>=threshold){
				val = i + 1;
				break;
			}
		}
	if(val<min_val) val = min_val;

	return val;
}

int64_t Paras::estimateMinReadsNumSupportSV(vector<int64_t> &mean_depth_vec){
	struct support_match_NODE{
		int32_t num;
		//vector<int32_t> suport_id;
		vector<int64_t> support_num;
	};

	vector<struct support_match_NODE*> support_match_vec;
	struct support_match_NODE* support_match_node;
	size_t i, j;
	int64_t max_id, max_num, same_id, min_supp_reads_num, total_depth;
	bool exist_flag = false;
	if(mean_depth_vec.size()==0){
		cout<<"paras minReadsNumSupportSV wrong"<<endl;
	}
	for(i=0; i<mean_depth_vec.size(); i++){
		for(j=0; j<support_match_vec.size(); j++){
			if(abs(mean_depth_vec.at(i)-support_match_vec.at(j)->support_num.at(0))<10){
				same_id = j;
				exist_flag = true;
			}else{
				exist_flag = false;
			}
		}
		if(exist_flag==false){
			support_match_node = new struct support_match_NODE();
			//support_match_node->suport_id.push_back(i);
			support_match_node->num = 1;
			support_match_node->support_num.push_back(mean_depth_vec.at(i));
			support_match_vec.push_back(support_match_node);
		}else{
			support_match_vec.at(same_id)->num += 1;
			//support_match_vec.at(same_id)->suport_id.push_back(i);
			support_match_vec.at(same_id)->support_num.push_back(mean_depth_vec.at(i));
		}
	}
	max_num = support_match_vec.at(0)->num;
	max_id = 0;
	for(i=0; i<support_match_vec.size(); i++){
		if(max_num<support_match_vec.at(i)->num){
			max_id = i;
			max_num = support_match_vec.at(i)->num;
		}
	}
	total_depth = 0;
	for(i=0; i<support_match_vec.at(max_id)->support_num.size(); i++){
		total_depth += support_match_vec.at(max_id)->support_num.at(i);
	}
	min_supp_reads_num = 2 + round(((double)total_depth / support_match_vec.at(max_id)->support_num.size())*MIN_SUPPORT_READS_NUM_FACTOR);

	//destroy and support_match_vec
	for(i=0; i<support_match_vec.size(); i++){
		vector<struct support_match_NODE*>::iterator svp;
		for(svp=support_match_vec.begin(); svp!=support_match_vec.end(); svp++) delete *svp;
		vector<struct support_match_NODE*>().swap(support_match_vec);
	}

	return min_supp_reads_num;
}

// parse long options
int Paras::parse_long_opt(int32_t option_index, const char *optarg, const struct option *lopts){
	int ret = 0;
	size_t find_pos;
	string lower_str;

	string opt_name_str = lopts[option_index].name;
	if(opt_name_str.compare("sample")==0){ // "sample"
		if(optarg) sample = optarg;
		else{
			cout << "Error: Please specify the correct sample name using --sample option." << endl << endl;
			ret = 1;
		}
	}else if(opt_name_str.compare("threads-per-assem-work")==0){ // "threads-per-assem-work"
		num_threads_per_assem_work = true;
	}else if(opt_name_str.compare("keep-assemble-reads")==0){ // "keep-assemble-reads"
		delete_reads_flag = false;
	}else if(opt_name_str.compare("keep-failed-reads")==0){ // "keep-failed-reads"
		keep_failed_reads_flag = true;
	}else if(opt_name_str.compare("reassemble-failed-work")==0){ // "reassemble-failed-work"
		reassemble_failed_work_flag = true;
	}else if(opt_name_str.compare("min-input-cov-assemble")==0){ // "min-input-cov-assemble"
		min_input_cov_canu = stoi(optarg);
	}else if(opt_name_str.compare("mask-noisy-region")==0){ // "mask-noisy-region"
		maskMisAlnRegFlag = true;
	}else if(opt_name_str.compare("assem-chunk-size")==0){ // assem-chunk-size
		assemChunkSize = stoi(optarg);
	}else if(opt_name_str.compare("assem-side-ext-size")==0){ // assem-side-ext-size
		assemSideExtSize = stoi(optarg);
	}else if(opt_name_str.compare("assem-side-ext-size-clip")==0){ // assem-side-ext-size-clip
		assemSideExtSizeClip = stoi(optarg);
	}else if(opt_name_str.compare("min-cons-read-size")==0){ // min-cons-read-size
		minConReadLen = stoi(optarg);
	}else if(opt_name_str.compare("monitor_proc_names_assemble")==0){ // monitor_proc_names_assemble
		monitoring_proc_names_assemble = optarg;
		if(monitoring_proc_names_assemble.size()>0){
			find_pos = monitoring_proc_names_assemble.find(" ");
			if(find_pos!=string::npos){
				cout << "Error: Monitored process names for 'assemble' step should not include blank characters." << endl << endl;
				ret = 1;
			}
		}else{
			cout << "Error: Please specify the correct monitored process names using '--monitor_proc_names_assemble' option." << endl << endl;
			ret = 1;
		}
	}else if(opt_name_str.compare("monitor_proc_names_call")==0){ // monitor_proc_names_call
		monitoring_proc_names_call = optarg;
		if(monitoring_proc_names_call.size()>0){
			find_pos = monitoring_proc_names_call.find(" ");
			if(find_pos!=string::npos){
				cout << "Error: Monitored process names for 'call' step should not include blank characters." << endl << endl;
				ret = 1;
			}
		}else{
			cout << "Error: Please specify the correct monitored process names using '--monitor_proc_names_call' option." << endl << endl;
			ret = 1;
		}
	}else if(opt_name_str.compare("max_proc_running_minutes_assemble")==0){ // max_proc_running_minutes_assemble
		max_proc_running_minutes_assemble = stoi(optarg);
		if(max_proc_running_minutes_assemble<ULTRA_LOW_PROC_RUNNING_MINUTES){
			cout << "Error: The specified maximum process running minutes is too small '" << max_proc_running_minutes_assemble << "', please specify a larger one at least " << ULTRA_LOW_PROC_RUNNING_MINUTES << "." << endl << endl;
			ret = 1;
		}
	}else if(opt_name_str.compare("max_proc_running_minutes_call")==0){ // max_proc_running_minutes_call
		max_proc_running_minutes_call = stoi(optarg);
		if(max_proc_running_minutes_call<ULTRA_LOW_PROC_RUNNING_MINUTES){
			cout << "Error: The specified maximum process running minutes is too small '" << max_proc_running_minutes_call << "', please specify a larger one at least " << ULTRA_LOW_PROC_RUNNING_MINUTES << "." << endl << endl;
			ret = 1;
		}
	}else if(opt_name_str.compare("technology")==0){ // technology
		technology = optarg;
		lower_str.resize(technology.size());
		transform(technology.begin(), technology.end(), lower_str.begin(), ::tolower);
		if(lower_str.compare(PACBIO_CLR_TECH_STR)==0 or lower_str.compare(PACBIO_CCS_TECH_STR)==0 or lower_str.compare(NANOPORE_TECH_STR)==0){
			technology = lower_str;
		}else{
			cout << "Error: Please specify the correct sequencing technology using '--technology' option." << endl << endl;
			ret = 1;
		}
	}else if(opt_name_str.compare("include-decoy")==0){ // include-decoy
		include_decoy = true;
	}else if(opt_name_str.compare("gt_min_sig_size")==0){ // "gt_min_sig_size"
		gt_min_sig_size = stoi(optarg);
	}else if(opt_name_str.compare("gt_size_ratio_match")==0){ // "gt_size_ratio_match"
		gt_size_ratio_match = stof(optarg);
	}else if(opt_name_str.compare("gt_min_alle_ratio")==0){ // "gt_min_alle_ratio"
		gt_min_alle_ratio = stof(optarg);
	}else if(opt_name_str.compare("gt_max_alle_ratio")==0){ // "gt_max_alle_ratio"
		gt_max_alle_ratio = stof(optarg);
	}


	return ret;
}
