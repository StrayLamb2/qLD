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
#ifndef PTHREAD_H
#define PTHREAD_H

#include "header.h"
#include <pthread.h>
#include <unistd.h>

#include <string.h>
#include <limits.h>
#include <ctype.h>
#include <stdbool.h>

#define EXIT 127
#define BUSYWAIT 0
#define CORRELATE 1

EXTERN_PTH pthread_mutex_t lock;
EXTERN_PTH pthread_mutex_t fileAccess;

/* Sense of each thread */
typedef struct localsense_t {
	bool lsense;
}localsense_t;

void sense_reversal_barrier_init(int num_threads);
void sense_reversal_barrier(int tid, int num_threads);
void sense_reversal_barrier_destroy(void);

/*
 * Calculation data for each thread. They mostly use the [0] thread's values
 * 	as only max, min and avg are actually unique for each thread
 * Everything else is just pointers to-be-used in the threads
 *
 */

typedef struct{
    int job_counter;
    double init_time;
    double mlt_time;
    double gemm_time;
    double ld_time;
    double write_time;
    long long total_lds;
}threadStats_t;

/*
 * Thread data. These data are unique to each thread, and include variables for:
 * 		the id,
 *		the total number of threads (this doesn't need to be unique),
 * 		if they are on a barrier,
 * 		the current operation of the thread
 * 		and the above data (correlateData).
 *
 */
typedef struct{
	int threadID;
	int threadTOTAL;
	int threadOP;
    FILE *threadLog;
    threadStats_t threadStats;
	float r2limit;
    int ploidy;
    int gpu;
    int blis;
    int mdf;
}threadData_t;

/*
 * Thread functions. These include:
 *		the initialization,
 *		the barrier syncing,
 * 		the operation function that each thread runs for the omega calculation,
 * 		the master function that assigns operations in each thread (only one
 * 			operation here, could be removed, stayed during developement in
 *			case it is needed, not removed for future scalability)
 * 		the actual master function, which starts the threads
 * 		the thread function, which includes barriers (busy-wait mode)
 * 		the termination of the workerThreads (threads != 0)
 *
 */
void initializeThreadData(threadData_t *cur, int i, int threads);
void setThreadArgs(threadData_t * threadData,
                   char *logname,
                   int tid,
                   float r2limit,
                   int ploidy,
                   int gpu,
                   int blis,
                   int mdf);
void updateThreadArgs(threadData_t * threadData,
                      float r2limit,
                      int ploidy,
                      int gpu,
                      int blis,
                      int mdf);
void startThreadOPS(threadData_t *threadData, int op);
void *thread(void *x);
void terminateWorkerThreads(pthread_t *workerThreadL, threadData_t *threadData);

void sense_reversal_barrier_init(int num_threads);
void sense_reversal_barrier(int tid, int num_threads);
void sense_reversal_barrier_destroy(void);

#endif
