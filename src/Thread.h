#ifndef _THREAD_H
#define _THREAD_H

#include <iostream>
#include <string>
#include <vector>

#include <pthread.h>
#include <unistd.h>
#include "Block.h"
#include "blatAlnTra.h"

using namespace std;

class Thread
{
	private:
		pthread_t tid;
		size_t user_tid;   // 0-based: 0, 1, 2, 3, ...
		static void* runDetect0(void* pVoid);  // the pointer to executing function
		void* runDetect1();  // inner executing method
		static void* runGenConsWorkOpt0(void* pVoid);  // the pointer to executing function
		void* runGenConsWorkOpt1();  // inner executing method
		static void* runCall0(void* pVoid);  // the pointer to executing function
		void* runCall1();  // inner executing method
		static void* runBlatAlnTra0(void* pVoid);  // the pointer to executing function
		void* runBlatAlnTra1();  // inner executing method

	public:
		vector<Block*> *blockVector;
		vector<varCand*> *var_cand_vec;
		vector<varCand*> *var_cand_clipReg_vec;
		vector<blatAlnTra*> *blat_aln_tra_vec;
		Thread();
		virtual ~Thread();
		virtual void runDetect() = 0;  // thread running entity
		virtual void runGenConsWorkOpt() = 0;  // thread running entity
		virtual void runCall() = 0;  // thread running entity
		virtual void runBlatAlnTra() = 0;  // thread running entity
		bool startDetect();  // start the thread for detection
		bool startGenConsWorkOpt();  // start the thread for detection
		bool startCall();  // start the thread for detection
		bool startBlatAlnTra();  // start the thread for detection
		pthread_t getThreadID();
		void setUserThreadID(size_t user_tid); // set user tid
		void setBlockVec(vector<Block*> *vec);
		void setVarCandVec(vector<varCand*> *var_cand_vec, vector<varCand*> *var_cand_clipReg_vec);
		void setBlatAlnTraVec(vector<blatAlnTra*> *blat_aln_tra_vec);
		size_t getUserThreadID();
		bool join();
};


// multi-thread
class MultiThread: public Thread {
	public:
		size_t num_threads = 0;
		void runDetect();
		void runGenConsWorkOpt();
		void runCall();
		void runBlatAlnTra();
		//void runFillVarSeq();
		void setNumThreads(size_t n);
};


#endif /* _THREAD_H */
