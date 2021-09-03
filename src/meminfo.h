#ifndef SRC_MEMINFO_H_
#define SRC_MEMINFO_H_

#include <iostream>

using namespace std;

// global variables
extern int64_t mem_total;		// in kB
//int64_t mem_free;		// in kB
//int64_t mem_available;	// in kB
extern int64_t swap_total;	// in kB

extern int64_t mem_seqAln;		// in kB
extern double mem_use_block_factor;
extern double swap_use_block_factor;
extern int32_t mem_block_seconds;	// block the computation by 10 second

extern int32_t work_num;

extern pthread_mutex_t mutex_mem;

int64_t getMemInfo(const char *name, int32_t index);


#endif /* SRC_MEMINFO_H_ */
