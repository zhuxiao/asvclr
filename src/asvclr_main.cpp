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
	paras.outputEstParas("After estimation:");

	genome.generateGenomeBlocks();

	// detect
	if(paras.command.compare("detect")==0 or paras.command.compare("all")==0){
		time.setStartTime();
		cout << "[" << time.getTime() << "]: detect structural variations ..." << endl;
		genome.genomeDetect();

		cout << "Total misAln region size: " << paras.misAlnRegLenSum << " bp" << endl;
		cout << "[" << time.getTime() << "]: detect structural variations finished." << endl;
		time.printSubCmdElapsedTime();
	}

	// assemble
	if(paras.command.compare("assemble")==0 or paras.command.compare("all")==0){
		time.setStartTime();
		cout << "[" << time.getTime() << "]: local assemble ..." << endl;

		// local assemble
		paras.slideSize = paras.assemSlideSize;
		genome.genomeLocalAssemble();

		cout << "[" << time.getTime() << "]: local assemble finished." << endl;
		time.printSubCmdElapsedTime();
	}

	// call
	if(paras.command.compare("call")==0 or paras.command.compare("all")==0){
		time.setStartTime();
		cout << "[" << time.getTime() << "]: call variants ..." << endl;

		// call variants
		paras.slideSize = paras.assemSlideSize;	// update slide size
		genome.genomeCall();

		cout << "[" << time.getTime() << "]: call variants finished." << endl;
		time.printSubCmdElapsedTime();
	}

	time.printOverallElapsedTime();

	return 0;
}
