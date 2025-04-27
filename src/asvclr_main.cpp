#include "Paras.h"
#include "Genome.h"
#include "util.h"

int main(int argc, char **argv) {
	Time time;

	Paras paras(argc, argv);

	if(paras.command.size()==0) return 1;

	// output parameters
	paras.outputParas();

	Genome genome(&paras);

	// estimate the parameters for noisy background
	genome.estimateSVSizeNum();
	paras.outputEstParas();

	genome.generateGenomeBlocks();

	// detect
	if(paras.command.compare(CMD_DET_STR)==0 or paras.command.compare(CMD_ALL_STR)==0 or paras.command.compare(CMD_DET_CNS_STR)==0){
		time.setStartTime();
		cout << "[" << time.getTime() << "]: detect structural variants ..." << endl;
		genome.genomeDetect();

		//cout << "Total misAln region size: " << paras.misAlnRegLenSum << " bp" << endl;
		cout << "[" << time.getTime() << "]: detect structural variants finished." << endl;
		time.printSubCmdElapsedTime();
	}

	// consensus
	if(paras.command.compare(CMD_CNS_STR)==0 or paras.command.compare(CMD_ALL_STR)==0 or paras.command.compare(CMD_DET_CNS_STR)==0){
		time.setStartTime();
		cout << "[" << time.getTime() << "]: local consensus ..." << endl;

		// local consensus
		genome.genomeLocalCons();

		cout << "[" << time.getTime() << "]: local consensus finished." << endl;
		time.printSubCmdElapsedTime();
	}

	// call
	if(paras.command.compare(CMD_CALL_STR)==0 or paras.command.compare(CMD_ALL_STR)==0){
		time.setStartTime();
		cout << "[" << time.getTime() << "]: call variants ..." << endl;

		// call variants
		genome.genomeCall();

		cout << "[" << time.getTime() << "]: call variants finished." << endl;
		time.printSubCmdElapsedTime();
	}

	time.printOverallElapsedTime();

	return 0;
}
