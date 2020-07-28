/*
qLD - High performance computation of Linkage disequilibrium
Copyright (C) 2020  C. Theodoris, N. Alachiotis

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/
#define EXTERN_PTH

#include "../include/header.h"
#include "../include/pthreads.h"
#include "../include/correlate_function__.h"

/* Helper global vars for the sense reversal barrier */
static int count;
static bool sense;

pthread_mutex_t lock;
pthread_mutex_t fileAccess;

/* Sense reversal barrier functions */
static localsense_t *localsense_list = NULL;

static inline void syncThreadsBARRIER(threadData_t *threadData);
static inline void execFuncMaster(threadData_t *threadData, int op);

void sense_reversal_barrier_init(int num_threads)
{
	int i;
	sense=true;
	count=num_threads;

	if(localsense_list == NULL)
    {
		localsense_list=(localsense_t *)malloc(sizeof(localsense_t)*
                (unsigned long)num_threads);
		assert(localsense_list != NULL);
	}

	for(i=0; i < num_threads; i++)
    {
		localsense_list[i].lsense=true;
	}
}

void sense_reversal_barrier(int tid, int num_threads)
{
	int threadno=tid;
	localsense_list[threadno].lsense=!localsense_list[threadno].lsense;

	if (__sync_fetch_and_sub (&count, 1) == 1)
    {
		count=num_threads;
		sense=localsense_list[threadno].lsense;
	}
	else
    {
		while(sense != localsense_list[threadno].lsense) __sync_synchronize();
	}
}

void sense_reversal_barrier_destroy(void)
{
	if (localsense_list != NULL)
		free(localsense_list);
	localsense_list=NULL;
}

void initializeThreadData(threadData_t *cur, int i, int threads)
{
	cur->threadID=i;
	cur->threadTOTAL=threads;
	cur->threadOP=BUSYWAIT;
    cur->threadStats.job_counter=0;
    cur->threadStats.init_time=0.0;
    cur->threadStats.mlt_time=0.0;
    cur->threadStats.gemm_time=0.0;
    cur->threadStats.ld_time=0.0;
    cur->threadStats.write_time=0.0;
    cur->threadStats.total_lds=0;
    cur->threadLog=NULL;
    cur->r2limit=0.2;
    cur->ploidy=0;
    cur->gpu=0;
    cur->blis=0;
    cur->mdf=0;
    cur->compQ=0;
    cur->task_count=0;
}

void setThreadArgs(threadData_t * threadData,
                   char * logname,
                   int tid,
                   float r2limit,
                   int ploidy,
                   int gpu,
                   int blis,
                   int mdf,
                   int compQ,
                   int task_count)
{
#if defined(VERBOSE) || defined(BENCHMARK)
    threadData[tid].threadLog=fopen(logname, "w");
#else
    logname=logname;
#endif
    threadData[tid].r2limit=r2limit;
    threadData[tid].ploidy=ploidy;
    threadData[tid].gpu=gpu;
    threadData[tid].blis=blis;
    threadData[tid].mdf=mdf;
    threadData[tid].compQ=compQ;  
    threadData[tid].task_count=task_count;  
}

void updateThreadArgs(threadData_t * threadData,
                      float r2limit,
                      int ploidy,
                      int gpu,
                      int blis,
                      int mdf,
                      int compQ,
                      int task_count)
{
	int threadIndex=0;
    char logname[20];
    int gpu_l=0;
    int blis_l=0;
#if defined(VERBOSE) || defined(BENCHMARK)
    char threadNo[11];
#endif
	for(threadIndex=0;threadIndex < threadData->threadTOTAL;threadIndex++)
    {
#if defined(VERBOSE) || defined(BENCHMARK)
        if(!(sprintf(threadNo,"%d",threadIndex)))
            exit(1);
        strcpy(logname,"thread_");
        strcat(logname,threadNo);
        strcat(logname,".log");
#else
        strcpy(logname,"");
#endif

#ifdef GPU
        if(threadIndex == (threadData->threadTOTAL-1) && gpu)
            gpu_l=1;
        else
            gpu_l=0;
#else
        gpu=gpu;
#endif

#ifdef CBLAS_USE
        blis_l=blis;
#else
        blis=blis;
        blis_l=0;
#endif

        setThreadArgs(threadData,
                      logname,
                      threadIndex,
                      r2limit,
                      ploidy,
                      gpu_l,
                      blis_l,
                      mdf,
                      compQ,
                      task_count);
    }
}

static inline void syncThreadsBARRIER(threadData_t *threadData)
{
	int threads=threadData[0].threadTOTAL;
	threadData[0].threadOP=BUSYWAIT;
	sense_reversal_barrier(0, threads);
}

static inline void execFuncMaster(threadData_t *threadData, int op)
{
	if(op == CORRELATE)
    {
#ifdef GPU
        if(threadData[0].gpu)
		    correlate_function_gpu(&threadData[0]);
        else
		    correlate_function(&threadData[0]);
#else
		correlate_function(&threadData[0]);
#endif
    }
}

static inline void setThreadOP(threadData_t *threadData, int op){
	int i, threads=threadData[0].threadTOTAL;

	for(i=0; i < threads; i++)
		threadData[i].threadOP=op;
}

void startThreadOPS(threadData_t *threadData, int op){
	setThreadOP(threadData, op);
	execFuncMaster(threadData, op);
	syncThreadsBARRIER(threadData);
}

static void pinToCore(int tid){
	cpu_set_t cpuset;
	CPU_ZERO(&cpuset);
	CPU_SET(tid, &cpuset);

	if(pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset) != 0)
    {
		fprintf(stdout, "\n ERROR: Please specify a number of threads that is \
smaller or equal\n        to the number of available physical cores (%d).\n\n",\
				tid);
		exit(1);
	}
}

void *thread(void *x)
{
	threadData_t *currentThread=(threadData_t *)x;

	int tid=currentThread->threadID;
    pinToCore(tid);

	int threads=currentThread->threadTOTAL;

	while (1)
    {
		__sync_synchronize();

		if(currentThread->threadOP == EXIT)
			return NULL;

		if(currentThread->threadOP == CORRELATE)
        {
#ifdef GPU
            if(currentThread->gpu)
		        correlate_function_gpu(currentThread);
            else
		        correlate_function(currentThread);
#else
		    correlate_function(currentThread);
#endif
			currentThread->threadOP=BUSYWAIT;
			sense_reversal_barrier(tid, threads);
		}
	}
	return NULL;
}

void terminateWorkerThreads(pthread_t *workerThreadL, threadData_t *threadData)
{
	int i, threads=threadData[0].threadTOTAL;
#if defined(VERBOSE) || defined(BENCHMARK)
    double total_init=0.0, total_mlt=0.0, total_gemm=0.0, total_ld=0.0, total_write=0.0;
    double total_t;
#endif
	for(i=0; i < threads; i++)
    {
#if defined(VERBOSE) || defined(BENCHMARK)
        total_t=(threadData[i].threadStats.init_time+
                 threadData[i].threadStats.mlt_time+
                 threadData[i].threadStats.gemm_time+
                 threadData[i].threadStats.ld_time+
                 threadData[i].threadStats.write_time);
        printf("\
Closing Thread[%d] with %d jobs done, %lld lds calculated\n",
    i, threadData[i].threadStats.job_counter,threadData[i].threadStats.total_lds);

        fprintf(threadData[i].threadLog, "Closing Thread[%d] with %d jobs done\n\
------------------------\n\
Thread accumulated time:\n\
\tInit: \t%.3fs\n\
\tMLT:  \t%.3fs\n\
\tGEMM: \t%.3fs\n\
\tLD:   \t%.3fs\n\
\tWrite:\t%.3fs\n\
------------------------\n\
\tTotal:\t%.3fs\n\
------------------------\n",
                i,
                threadData[i].threadStats.job_counter,
                threadData[i].threadStats.init_time,
                threadData[i].threadStats.mlt_time,
                threadData[i].threadStats.gemm_time,
                threadData[i].threadStats.ld_time,
                threadData[i].threadStats.write_time,
                total_t);
        fclose(threadData[i].threadLog);

        total_init+=threadData[i].threadStats.init_time;
        total_mlt+=threadData[i].threadStats.mlt_time;
        total_gemm+=threadData[i].threadStats.gemm_time;
        total_ld+=threadData[i].threadStats.ld_time;
        total_write+=threadData[i].threadStats.write_time;
#endif
		threadData[i].threadOP=EXIT;
    }
	for(i=1; i < threads; i++)
		pthread_join(workerThreadL[i-1], NULL);

#if defined(VERBOSE) || defined(BENCHMARK)
    printf("\
------------------------\n\
Total accumulated time:\n\
\tInit: \t %.3fs\n\
\tMLT:  \t %.3fs\n\
\tGEMM: \t %.3fs\n\
\tLD:   \t %.3fs\n\
\tWrite:\t %.3fs\n\
------------------------\n\
\tTotal: \t %.3fs\n\
------------------------\n",
                total_init,
                total_mlt,
                total_gemm,
                total_ld,
                total_write,
                (total_init+
                 total_mlt+
                 total_gemm+
                 total_ld+
                 total_write));
#endif
}
