#include "Paras.h"
#include "Chrome.h"
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

	// detect indels
	if(paras.command=="detect")
	{
		cout << "[" << time.getTime() << "]: detect structural variations ..." << endl;

		genome.genomeDetect();

		cout << "Total misAln region size: " << paras.misAlnRegLenSum << endl;
		cout << "[" << time.getTime() << "]: detect structural variations finished." << endl;
	}

	// assemble
	if(paras.command=="assemble")
	{
		cout << "[" << time.getTime() << "]: local assemble ..." << endl;

		// local assemble
		genome.genomeLocalAssemble();

		cout << "[" << time.getTime() << "]: local assemble finished." << endl;
	}

	// call
	if(paras.command=="call")
	{
		cout << "[" << time.getTime() << "]: call variants ..." << endl;

		// call variants
		genome.genomeCall();

		cout << "[" << time.getTime() << "]: call variants finished." << endl;
	}

	// save results in VCF file format
	//genome.saveResultVCF();

	return 0;
}
