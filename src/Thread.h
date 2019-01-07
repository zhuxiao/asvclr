#ifndef _THREAD_H
#define _THREAD_H

#include <iostream>
#include <string>
#include <vector>

#include <pthread.h>
#include <unistd.h>
#include "Block.h"

using namespace std;

class Thread
{
	private:
		pthread_t tid;
		size_t user_tid;   // 0-based: 0, 1, 2, 3, ...
		static void* runDetect0(void* pVoid);  // the pointer to executing function
		void* runDetect1();  // inner executing method
		static void* runAssemble0(void* pVoid);  // the pointer to executing function
		void* runAssemble1();  // inner executing method
		static void* runCall0(void* pVoid);  // the pointer to executing function
		void* runCall1();  // inner executing method

	public:
		vector<Block*> *blockVector;
		vector<varCand*> *var_cand_vec;
		vector<varCand*> *var_cand_clipReg_vec;
		Thread();
		virtual ~Thread();
		virtual void runDetect() = 0;  // thread running entity
		virtual void runAssemble() = 0;  // thread running entity
		virtual void runCall() = 0;  // thread running entity
		bool startDetect();  // start the thread for detection
		bool startAssemble();  // start the thread for detection
		bool startCall();  // start the thread for detection
		pthread_t getThreadID();
		void setUserThreadID(size_t user_tid); // set user tid
		void setBlockVec(vector<Block*> *vec);
		void setVarCandVec(vector<varCand*> *var_cand_vec, vector<varCand*> *var_cand_clipReg_vec);
		size_t getUserThreadID();
		bool join();
};

Thread::Thread()
{
	tid = 0;
	user_tid = 0;
	blockVector = NULL;
	var_cand_vec = NULL;
	var_cand_clipReg_vec = NULL;
}

Thread::~Thread(){
}

void* Thread::runDetect0(void* pVoid)
{
	Thread* p = (Thread*) pVoid;
	p->runDetect1();
	return p;
}

void* Thread::runDetect1()
{
	runDetect();
	pthread_exit(NULL);
}

bool Thread::startDetect()
{
	return pthread_create(&tid, NULL, runDetect0, this) == 0;
}

void* Thread::runAssemble0(void* pVoid)
{
	Thread* p = (Thread*) pVoid;
	p->runAssemble1();
	return p;
}

void* Thread::runAssemble1()
{
	runAssemble();
	pthread_exit(NULL);
}

bool Thread::startAssemble()
{
	return pthread_create(&tid, NULL, runAssemble0, this) == 0;
}

void* Thread::runCall0(void* pVoid)
{
	Thread* p = (Thread*) pVoid;
	p->runCall1();
	return p;
}

void* Thread::runCall1()
{
	runCall();
	pthread_exit(NULL);
}

bool Thread::startCall()
{
	return pthread_create(&tid, NULL, runCall0, this) == 0;
}

pthread_t Thread::getThreadID()
{
	return tid;
}

void Thread::setUserThreadID(size_t user_tid){
	this->user_tid = user_tid;
}

void Thread::setBlockVec(vector<Block*> *blockVector){
	this->blockVector = blockVector;
}
void Thread::setVarCandVec(vector<varCand*> *var_cand_vec, vector<varCand*> *var_cand_clipReg_vec){
	this->var_cand_vec = var_cand_vec;
	this->var_cand_clipReg_vec = var_cand_clipReg_vec;
}

size_t Thread::getUserThreadID(){
	return user_tid;
}

bool Thread::join()
{
	return pthread_join(tid, NULL) == 0;
}

// multi-thread
class MultiThread : public Thread {
	public:
		size_t num_threads = 0;
		void runDetect();
		void runAssemble();
		void runCall();
		void runFillVarSeq();
		void setNumThreads(size_t n);
};

void MultiThread::setNumThreads(size_t n){
	num_threads = n;
}

void MultiThread::runDetect(){
	size_t user_tid_tmp = getUserThreadID();
	Block *bloc;
	for (size_t i=0; i < blockVector->size(); i++){
		if(i%num_threads==user_tid_tmp){
			bloc = blockVector->at(i);
			bloc->blockDetect();
		}
	}
}

void MultiThread::runAssemble(){
	size_t user_tid_tmp = getUserThreadID();
	Block *bloc;
	for (size_t i=0; i < blockVector->size(); i++){
		if(i%num_threads==user_tid_tmp){
			bloc = blockVector->at(i);
			bloc->blockLocalAssemble();
		}
	}
}

void MultiThread::runCall(){
	size_t i, user_tid_tmp = getUserThreadID();
	varCand *var_cand;
	for (i=0; i < var_cand_vec->size(); i++){
		if(i%num_threads==user_tid_tmp){
			var_cand = var_cand_vec->at(i);
			var_cand->callVariants();
		}
	}
	for (i=0; i < var_cand_clipReg_vec->size(); i++){
		if(i%num_threads==user_tid_tmp){
			var_cand = var_cand_clipReg_vec->at(i);
			var_cand->callVariants();
		}
	}
}

#endif /* _THREAD_H */
