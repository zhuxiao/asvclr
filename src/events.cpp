#include "events.h"

// allocate indel event
indelEvent_t* allocateIndelEvent(uint32_t startPos, string& seq){
	indelEvent_t* indelE = new indelEvent_t;
	if(!indelE){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot allocate memory" << endl;
		exit(1);
	}

	indelE->startPos = startPos;
	indelE->seq = seq;

	return indelE;
}

// allocate clipping event
clipEvent_t* allocateClipEvent(uint32_t startPos, uint16_t opflag, uint16_t endFlag, string& seq){
	clipEvent_t* clipE = new clipEvent_t;
	if(!clipE){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot allocate memory" << endl;
		exit(1);
	}

	clipE->startPos = startPos;
	clipE->opflag = opflag;
	clipE->endFlag = endFlag;
	clipE->seq = seq;

	return clipE;
}
