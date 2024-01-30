#include "Thread.h"


Thread::Thread()
{
	tid = 0;
	user_tid = 0;
	blockVector = NULL;
	var_cand_vec = NULL;
	var_cand_clipReg_vec = NULL;
	blat_aln_tra_vec = NULL;
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

void* Thread::runGenConsWorkOpt0(void* pVoid)
{
	Thread* p = (Thread*) pVoid;
	p->runGenConsWorkOpt1();
	return p;
}

void* Thread::runGenConsWorkOpt1()
{
	runGenConsWorkOpt();
	pthread_exit(NULL);
}

bool Thread::startGenConsWorkOpt()
{
	return pthread_create(&tid, NULL, runGenConsWorkOpt0, this) == 0;
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

void* Thread::runBlatAlnTra0(void* pVoid)
{
	Thread* p = (Thread*) pVoid;
	p->runBlatAlnTra1();
	return p;
}

void* Thread::runBlatAlnTra1()
{
	runBlatAlnTra();
	pthread_exit(NULL);
}

bool Thread::startBlatAlnTra()
{
	return pthread_create(&tid, NULL, runBlatAlnTra0, this) == 0;
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

void Thread::setBlatAlnTraVec(vector<blatAlnTra*> *blat_aln_tra_vec){
	this->blat_aln_tra_vec = blat_aln_tra_vec;
}

size_t Thread::getUserThreadID(){
	return user_tid;
}

bool Thread::join()
{
	return pthread_join(tid, NULL) == 0;
}

void MultiThread::setNumThreads(size_t n){
	num_threads = n;
}

void MultiThread::runDetect(){
	size_t user_tid_tmp = getUserThreadID();
	Block *bloc;
	for (size_t i=0; i < blockVector->size(); i++){
		if(i%num_threads==user_tid_tmp){
			bloc = blockVector->at(i);
			if(bloc->process_flag) bloc->blockDetect();
		}
	}
}

void MultiThread::runGenConsWorkOpt(){
	size_t user_tid_tmp = getUserThreadID();
	Block *bloc;
	for (size_t i=0; i < blockVector->size(); i++){
		if(i%num_threads==user_tid_tmp){
			bloc = blockVector->at(i);
			if(bloc->process_flag) bloc->blockGenerateLocalConsWorkOpt();
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

void MultiThread::runBlatAlnTra(){
	size_t i, user_tid_tmp = getUserThreadID();
	blatAlnTra *blat_aln_tra;
	for (i=0; i < blat_aln_tra_vec->size(); i++){
		if(i%num_threads==user_tid_tmp){
			blat_aln_tra = blat_aln_tra_vec->at(i);
			blat_aln_tra->generateBlatResult();
		}
	}
}
