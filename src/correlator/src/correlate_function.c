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
#include "../include/header.h"
#include "../include/correlate_function__.h"
#include "../include/pthreads.h"
#include "../include/correlate_MEM.h"
#include "../include/input_list.h"
#include "../include/read_file.h"
#include "../include/correlate_IO.h"
#include "../include/correlate_IOx64.h"
#include "../include/gemm_cpu.h"

void correlate_function(threadData_t *threadData)
{
    int first_run=1; // Used to create SNPtables the first time the thread is called
    unsigned int tableSize1_s=BUFFERSIZE;
    unsigned int tableSize2_s=BUFFERSIZE;

#if defined(VERBOSE) || defined( BENCHMARK)
    double tet_s=gettime(), tet_e=0.0;          //total external time start/end
    double total_t=0.0;                         //total from internal timers
    double twc_s=0.0, twc_e=0.0, twc_na=0.0;    //time without correlate start/end
    double time1=0.0,time2=0.0,time3=0.0;
    int lds;

#endif

    char *writeBuffer=(char *)malloc(FILEBUFFERSIZE*sizeof(char));
    assert(writeBuffer);

    // Allocate memory
    task_t *nodeData=create_task_t(1);
    task_t *nodeData2=create_task_t(2);
    table_x64 *tableData=create_table_x64(tableSize1_s);
    table_x64 *tableData2=create_table_x64(tableSize2_s);
    helper_t *helperData=create_helper_t(threadData->ploidy);

    while(1)
    {
#if defined(VERBOSE) || defined(BENCHMARK)
        twc_s=gettime();
#endif
        // Competing queue
        pthread_mutex_lock(&lock);
        t_node *new_node=dequeue_task();
        pthread_mutex_unlock(&lock);

        if(!new_node) // If we run out of tasks
            break;

        threadData[0].threadStats.job_counter++;

        // Buffer output for bulk writes
        outFileType fpOut=FOPEN_OUT(new_node->output_file,"wb");
        assert(fpOut != NULL);
        setvbuf(fpOut, writeBuffer, _IOFBF, FILEBUFFERSIZE);

#if defined(VERBOSE) || defined( BENCHMARK)
        time1=gettime();
        fprintf(threadData[0].threadLog, "Thread[%d] Task: %d\n",
                threadData[0].threadID, new_node->id);
#endif
        // Check if inputs are the same (printing purposes)
        if(new_node->inPath2set == 1 && !strcmp(new_node->input_file,
                                                new_node->input_file2))
        {
#ifdef VERBOSE
            fprintf(threadData[0].threadLog,"Input 2 same as Input, using Input only.\n");
#endif
            new_node->inPath2set = 0;
        }
#ifdef VERBOSE
        fprintf(threadData[0].threadLog, "\nInput\n");
#endif

        // Sample Lists create validity map. 
        // If no Sample List is provided, the whole valid map is '1's
        valid_t *validData=NULL;
        valid_t *validData2=NULL;
        validData=create_valid_t(new_node->snipSize);
        makeValidList(new_node->sampleList, new_node->headerLine2, validData);
        validData2=create_valid_t(new_node->snipSize2);
        makeValidList(new_node->sampleList2, new_node->headerLine2_2, validData2);

        // Compressed size is calculated from size of valid samples
        unsigned int compSize, compSize2;
        compSize=(validData->valid_count/(sizeof(inputDataType_x64)*8)+
                (validData->valid_count%(sizeof(inputDataType_x64)*8)!=0?1:0));
        compSize2=(validData2->valid_count/(sizeof(inputDataType_x64)*8)+
                (validData2->valid_count%(sizeof(inputDataType_x64)*8)!=0?1:0));

        if(first_run)
        {
            first_run=0;
            // Reallocate memory of SNPtables according to the compressed size
            tableData->SNPtable=(inputDataType_x64*)realloc_buff(
                                    (void *)tableData->SNPtable,
                                    tableSize1_s*compSize,
                                    "inputDataType_x64");
            assert(tableData->SNPtable);
            tableData2->SNPtable=(inputDataType_x64*)realloc_buff(
                                    (void *)tableData2->SNPtable,
                                    tableSize2_s*compSize2,
                                    "inputDataType_x64");
            assert(tableData2->SNPtable);
        }

        // Fill nodeData struct for input1
        nodeData->tableSize=new_node->posWmax-new_node->posWmin;
        nodeData->filesList=new_node->filesList;
        nodeData->filesListNum=new_node->filesListNum;
        nodeData->snipSize=new_node->snipSize;
        nodeData->headerLine1=new_node->headerLine1;
        nodeData->headerLine2=new_node->headerLine2;
        nodeData->posWmin=new_node->posWmin;
        nodeData->posWmax=new_node->posWmax;
        nodeData->posWset=new_node->posWset1;
        nodeData->valid_count=validData->valid_count;
        nodeData->valid_mask=validData->validList;

        // (Re-)Initialize values on tableData & helperData
        tableData->tableIndex=0;
        tableData->compSize=compSize;
        tableData->compIndex=0;
        helperData->snipCharIndex=0;

        // Read tableA from input files to memory
        readTable_x64(threadData,
                      nodeData,
                      tableData,
                      helperData);

        // Read TableB from different input file, or with different pos window
        if(new_node->inPath2set == 1)
        {
            // Fill nodeData struct for input2
            nodeData2->tableSize=new_node->posWmax2-new_node->posWmin2;
            nodeData2->filesList=new_node->filesList2;
            nodeData2->filesListNum=new_node->filesListNum2;
            nodeData2->snipSize=new_node->snipSize2;
            nodeData2->headerLine1=new_node->headerLine1_2;
            nodeData2->headerLine2=new_node->headerLine2_2;
            nodeData2->posWmin=new_node->posWmin2;
            nodeData2->posWmax=new_node->posWmax2;
            nodeData2->posWset=new_node->posWset2;
            nodeData2->valid_count=validData2->valid_count;
            nodeData2->valid_mask=validData2->validList;

            // (Re-)Initialize values on tableData2
            tableData2->tableSize=tableSize2_s;
            tableData2->tableIndex=0;
            tableData2->compSize=compSize2;
            tableData2->compIndex=0;

            // Read tableB from input files to memory
            readTable_x64(threadData,
                          nodeData2,
                          tableData2,
                          helperData);
        }
        else if(new_node->posWset2 == 1)
        {
            // Fill nodeData struct for input2
            nodeData2->tableSize=new_node->posWmax-new_node->posWmin;
            nodeData2->filesList=new_node->filesList;
            nodeData2->filesListNum=new_node->filesListNum;
            nodeData2->snipSize=new_node->snipSize;
            nodeData2->headerLine1=new_node->headerLine1;
            nodeData2->headerLine2=new_node->headerLine2;
            nodeData2->posWmin=new_node->posWmin2;
            nodeData2->posWmax=new_node->posWmax2;
            nodeData2->posWset=new_node->posWset2;
            nodeData2->valid_count=validData2->valid_count;
            nodeData2->valid_mask=validData2->validList;
            
            // (Re-)Initialize values on tableData2
            tableData2->tableSize=tableSize1_s;
            tableData2->tableIndex=0;
            tableData2->compSize=compSize;
            tableData2->compIndex=0;
            
            // Read tableB from input files to memory
            readTable_x64(threadData,
                          nodeData2,
                          tableData2,
                          helperData);
        }
        else
        {
#ifdef VERBOSE
            fprintf(threadData[0].threadLog, "Table B same as Table A\n");
#else
            ;
#endif

        }

#if defined(VERBOSE) || defined(BENCHMARK)
        time2 = gettime();
        twc_na=gettime()-twc_s;
        twc_e+=twc_na;
        fprintf(threadData[0].threadLog,"\n\n[T] Disk Access times: %f\n", twc_na);
#endif

        // Correlate single or dual input modes
        if(new_node->posWset2 == 0 && new_node->inPath2set == 0)
        {
#ifndef COUNTSTATES
            correlate(threadData,
                      fpOut,
                      tableData->SNPtable,
                      tableData->POStable,
                      tableData->IDtable,
                      tableData->BCtable,
                      tableData->tableIndex,
                      tableData->SNPtable,
                      tableData->POStable,
                      tableData->IDtable,
                      tableData->BCtable,
                      tableData->tableIndex,
                      tableData->compSize,
                      new_node->snipSize,
                      new_node->posWset2);
#else
            countStates_mode(fpOut,
                             tableData,
                             tableData,
                             new_node->posWset2);

#endif
        }
        else
        {
#ifndef COUNTSTATES
            correlate(threadData,
                      fpOut,
                      tableData->SNPtable,
                      tableData->POStable,
                      tableData->IDtable,
                      tableData->BCtable,
                      tableData->tableIndex,
                      tableData2->SNPtable,
                      tableData2->POStable,
                      tableData2->IDtable,
                      tableData2->BCtable,
                      tableData2->tableIndex,
                      tableData2->compSize,
                      new_node->snipSize,
                      new_node->posWset2);
#else
            countStates_mode(fpOut,
                             tableData,
                             tableData2,
                             new_node->posWset2);

#endif
        }
#if defined(VERBOSE) || defined(BENCHMARK)
        time3=gettime();
        if(new_node->posWset2 == 0 && new_node->inPath2set == 0)
        {
            lds=tableData->tableIndex*tableData->tableIndex;
            threadData->threadStats.total_lds+=lds;
            fprintf(threadData->threadLog, "%d snips processed, results: %d\n",
                tableData->tableIndex*2, lds);
            fprintf(threadData->threadLog, "Create Tables: %fs, Correlations: %fs, "
                                             "Total Time: %fs \n",
                                             time2-time1, time3-time2, time3-time1);
            fprintf(threadData->threadLog, "----------------------------------------\n\n");

        }
        else
        {
            lds=tableData->tableIndex*tableData2->tableIndex;
            threadData->threadStats.total_lds+=tableData->tableIndex*tableData2->tableIndex;
            fprintf(threadData->threadLog, "%d snips processed, results: %d\n",
                tableData->tableIndex + tableData2->tableIndex, lds);
            fprintf(threadData->threadLog, "Create Tables: %fs, Correlations: %fs, "
                                             "Total Time: %fs \n",
                                             time2-time1, time3-time2, time3-time1);
            fprintf(threadData->threadLog, "----------------------------------------\n\n");

        }
#endif
        // Free and close task-centric allocated space and files
        fflush(fpOut);
        FCLOZE_OUT(fpOut);
        free_task(new_node);
        free_valid_t(validData);
        free_valid_t(validData2);
    }
    // Free thread-centric allocated space
    free(writeBuffer);
    free_task_t(nodeData);
    free_task_t(nodeData2);
    free_table_x64(tableData);
    free_table_x64(tableData2);
    free_helper_t(helperData);

#if defined(VERBOSE) || defined(BENCHMARK)
    tet_e=gettime()-tet_s;
    total_t=(threadData[0].threadStats.init_time+
             threadData[0].threadStats.mlt_time+
             threadData[0].threadStats.gemm_time+
             threadData[0].threadStats.ld_time+
             threadData[0].threadStats.write_time);
    printf("Thread[%d]:\n"
            "\tAccumulated Time within correlate:       %.3fs\n"
            "\tAccumulated Time outside of correlate:   %.3fs\n"
            "\tAccumulated Total time (external timer): %.3fs\n",
            threadData[0].threadID, total_t, twc_e, tet_e);
#endif
}
