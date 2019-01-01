#include <string.h>
#include <vector>
#include <htslib/hts.h>
#include "Paras.h"

// Constructor with parameters
Paras::Paras(){
	init();
}

// Constructor with parameters
Paras::Paras(int argc, char **argv)
{
	init();
	if(parseParas(argc, argv)!=0) exit(1);
}

// initialization
void Paras::init()
{
	command = "";
	inBamFile = "";
	outFilePrefix = "";
	num_threads = 0;

	min_ins_size_filt = 0;
	min_del_size_filt = 0;
	min_clip_size_filt = 0;
	min_ins_num_filt = 0;
	min_del_num_filt = 0;
	min_clip_num_filt = 0;
}

// parse the parameters
int Paras::parseParas(int argc, char **argv)
{
	if (argc < 2) { showUsage(); return 1; }

    if (strcmp(argv[1], "-h") == 0 or strcmp(argv[1], "help") == 0 or strcmp(argv[1], "--help") == 0) {
        if (argc == 2) { showUsage(); return 0; }
        argv++;
        argc = 2;
    }

    if (strcmp(argv[1], "detect")==0){
    	command = "detect";
    	return parseDetectParas(argc-1, argv+1);
    }else if(strcmp(argv[1], "assemble")==0){
    	command = "assemble";
    	return parseAssembleParas(argc-1, argv+1);
    }else if(strcmp(argv[1], "call")==0){
    	command = "call";
    	return parseCallParas(argc-1, argv+1);
    }else if(strcmp(argv[1], "all")==0){
    	command = "all";
    	showUsage(); return 1;
    }else{
    	cerr << "invalid command " << argv[1] << endl;
    	showUsage(); return 1;
    }
}

// parse the parameters for detect command
int Paras::parseDetectParas(int argc, char **argv)
{
	int opt, threadNum_tmp = 0;
	blockSize = BLOCKSIZE;
	slideSize = SLIDESIZE;
	min_sv_size_usr = MIN_SV_SIZE_USR;
	maskMisAlnRegFlag = false;

	while( (opt = getopt(argc, argv, ":f:b:s:o:t:Mh")) != -1 ){
		switch(opt){
			case 'f': refFile = optarg; break;
			case 'b': blockSize = stoi(optarg); break;
			case 's': slideSize = stoi(optarg); break;
			case 'm': min_sv_size_usr = stoi(optarg); break;
			case 'o': outFilePrefix = optarg; break;
			case 't': threadNum_tmp = stoi(optarg); break;
			case 'M': maskMisAlnRegFlag = true; break;
			case 'h': showDetectUsage(); exit(0);
			case '?': cout << "unknown option -" << (char)optopt << endl; exit(1);
			case ':': cout << "the option -" << (char)optopt << " needs a value" << endl; exit(1);
		}
	}

	num_threads = (threadNum_tmp>=sysconf(_SC_NPROCESSORS_ONLN)) ? sysconf(_SC_NPROCESSORS_ONLN) : threadNum_tmp;

	opt = argc - optind; // the number of SAMs on the command line
	if(opt==1) inBamFile = argv[optind];
	else { showDetectUsage(); return 1; }

	if(refFile.size()==0){
		cout << "Error: Please specify the reference" << endl << endl;
		showDetectUsage();
		return 1;
	}

	return 0;
}

// parse the parameters for assemble command
int Paras::parseAssembleParas(int argc, char **argv)
{
	int opt, threadNum_tmp = 0;
	blockSize = BLOCKSIZE;
	slideSize = ASSEM_SLIDE_SIZE;
	min_sv_size_usr = MIN_SV_SIZE_USR;
	maskMisAlnRegFlag = false;

	while( (opt = getopt(argc, argv, ":f:b:S:o:t:Mh")) != -1 ){
		switch(opt){
			case 'f': refFile = optarg; break;
			case 'b': blockSize = stoi(optarg); break;
			case 'S': slideSize = stoi(optarg); break;
			case 'm': min_sv_size_usr = stoi(optarg); break;
			case 'o': outFilePrefix = optarg; break;
			case 't': threadNum_tmp = stoi(optarg); break;
			case 'M': maskMisAlnRegFlag = true; break;
			case 'h': showAssembleUsage(); exit(0);
			case '?': cout << "unknown option -" << (char)optopt << endl; exit(1);
			case ':': cout << "the option -" << (char)optopt << " needs a value" << endl; exit(1);
		}
	}

	num_threads = (threadNum_tmp>=sysconf(_SC_NPROCESSORS_ONLN)) ? sysconf(_SC_NPROCESSORS_ONLN) : threadNum_tmp;

	opt = argc - optind; // the number of SAMs on the command line
	if(opt==1) inBamFile = argv[optind];
	else { showAssembleUsage(); return 1; }

	if(refFile.size()==0){
		cout << "Error: Please specify the reference" << endl << endl;
		showAssembleUsage();
		return 1;
	}

	return 0;
}

// parse the parameters for call command
int Paras::parseCallParas(int argc, char **argv)
{
	int opt, threadNum_tmp = 0;
	blockSize = BLOCKSIZE;
	slideSize = ASSEM_SLIDE_SIZE;
	min_sv_size_usr = MIN_SV_SIZE_USR;
	maskMisAlnRegFlag = false;

	while( (opt = getopt(argc, argv, ":f:b:S:o:t:Mh")) != -1 ){
		switch(opt){
			case 'f': refFile = optarg; break;
			case 'b': blockSize = stoi(optarg); break;
			case 'S': slideSize = stoi(optarg); break;
			case 'm': min_sv_size_usr = stoi(optarg); break;
			case 'o': outFilePrefix = optarg; break;
			case 't': threadNum_tmp = stoi(optarg); break;
			case 'M': maskMisAlnRegFlag = true; break;
			case 'h': showCallUsage(); exit(0);
			case '?': cout << "unknown option -" << (char)optopt << endl; exit(1);
			case ':': cout << "the option -" << (char)optopt << " needs a value" << endl; exit(1);
		}
	}

	num_threads = (threadNum_tmp>=sysconf(_SC_NPROCESSORS_ONLN)) ? sysconf(_SC_NPROCESSORS_ONLN) : threadNum_tmp;

	opt = argc - optind; // the number of SAMs on the command line
	if(opt==1) inBamFile = argv[optind];
	else { showCallUsage(); return 1; }

	if(refFile.size()==0){
		cout << "Error: Please specify the reference" << endl << endl;
		showCallUsage();
		return 1;
	}

	return 0;
}

// show the usage
void Paras::showUsage()
{
	cout << "Program: cindel_lr (call the indels for genome alignments)" << endl;
	cout << "Version: 0.1.0 (using htslib " << hts_version() << ")" << endl << endl;
	cout << "Usage:  cindel_lr  <command> [options]" << endl << endl;

	cout << "Commands:" << endl;
	cout << "     detect       detect indel signatures in aligned reads" << endl;
	cout << "     assemble     assemble candidate regions and align assemblies" << endl;
	cout << "                  back to reference" << endl;
	cout << "     call         call indels by alignments of local genome assemblies" << endl;
	cout << "     all          run the above commands in turn" << endl;
}

// show the usage for detect command
void Paras::showDetectUsage()
{
	cout << "Program: cindel_lr (compute the coverage for chromosome)" << endl;
	cout << "Version: 0.1.0 (using htslib " << hts_version() << ")" << endl << endl;
	cout << "Usage: cindel_lr detect [options] <in.bam>|<in.sam> [region ...]?" << endl << endl;

	cout << "Options: " << endl;
	cout << "     -f FILE      reference file name (required)" << endl;
	cout << "     -b INT       block size [1000000]" << endl;
	cout << "     -s INT       detect slide size [500]" << endl;
	cout << "     -m INT       minimal SV size to detect [2]" << endl;
	cout << "     -o FILE      prefix of the output file [stdout]" << endl;
	cout << "     -t INT       number of threads [0]" << endl;
	cout << "     -M           Mask mis-aligned regions, default: false" << endl;
	cout << "     -h           show this help message and exit" << endl;
}

// show the usage for assemble command
void Paras::showAssembleUsage()
{
	cout << "Program: cindel_lr (compute the coverage for chromosome)" << endl;
	cout << "Version: 0.1.0 (using htslib " << hts_version() << ")" << endl << endl;
	cout << "Usage: cindel_lr assemble [options] <in.bam>|<in.sam> [region ...]?" << endl << endl;

	cout << "Options: " << endl;
	cout << "     -f FILE      reference file name (required)" << endl;
	cout << "     -b INT       block size [1000000]" << endl;
	cout << "     -S INT       assemble slide size [10000]" << endl;
	cout << "     -o FILE      prefix of the output file [stdout]" << endl;
	cout << "     -t INT       number of threads [0]" << endl;
	cout << "     -M           Mask mis-aligned regions, default: false" << endl;
	cout << "     -h           show this help message and exit" << endl;
}

// show the usage for call command
void Paras::showCallUsage()
{
	cout << "Program: cindel_lr (compute the coverage for chromosome)" << endl;
	cout << "Version: 0.1.0 (using htslib " << hts_version() << ")" << endl << endl;
	cout << "Usage: cindel_lr call [options] <in.bam>|<in.sam> [region ...]?" << endl << endl;

	cout << "Options: " << endl;
	cout << "     -f FILE      reference file name (required)" << endl;
	cout << "     -b INT       block size [1000000]" << endl;
	cout << "     -S INT       assemble slide size [10000]" << endl;
	cout << "     -o FILE      prefix of the output file [stdout]" << endl;
	cout << "     -t INT       number of threads [0]" << endl;
	cout << "     -M           Mask mis-aligned regions, default: false" << endl;
	cout << "     -h           show this help message and exit" << endl;
}

// output parameters
void Paras::outputParas(){
	if(refFile.size()) cout << "Ref file: " << refFile << endl;
	if(inBamFile.size()) cout << "In file: " << inBamFile << endl;
	if(outFilePrefix.size()) cout << "Out prefix: " << outFilePrefix << endl;

	cout << "Block size: " << blockSize << endl;
	cout << "Slide size: " << slideSize << endl;
	cout << "Num threads: " << num_threads << endl;
	if(maskMisAlnRegFlag) cout << "Mask mis-aligned regions: yes" << endl << endl;
	else cout << "Mask mis-aligned regions: no" << endl << endl;
}

// output the estimation parameters
void Paras::outputEstParas(string info){
	cout << info << endl;
	cout << "min_ins_size_filt: " << min_ins_size_filt << endl;
	cout << "min_del_size_filt: " << min_del_size_filt << endl;
//	cout << "min_clip_size_filt: " << min_clip_size_filt << endl;
	cout << "min_ins_num_filt: " << min_ins_num_filt << endl;
	cout << "min_del_num_filt: " << min_del_num_filt << endl;
	cout << "min_clip_num_filt: " << min_clip_num_filt << endl;
	cout << "large_indel_size_thres: " << large_indel_size_thres << endl << endl;
}

// initialize the estimation auxiliary data
void Paras::initEst(){
	size_t i;
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
size_t Paras::estimateSinglePara(size_t *arr, size_t n, double threshold, size_t min_val){
	size_t i, total, val = 1;
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

