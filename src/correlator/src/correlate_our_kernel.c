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
#include "../include/correlate_our_kernel.h"
#include "../include/correlate_MEM.h"
#include "../include/correlate_IO.h"
#include "../include/pthreads.h"
#include "../include/input_list.h"
#include "../include/load_balancer.h"
#include "../include/read_file.h"
#include "../include/gemm_gpu.h"

long get_corecount(void)
{
    FILE *cpuinfo=fopen("/proc/cpuinfo", "r");
    assert(cpuinfo);
    char line[STRINGLENGTH];
    long cores;
    while(fgets(line, STRINGLENGTH, cpuinfo))
    {
        if(strstr(line, "cpu cores"))
        {
            sscanf(line,"cpu cores\t\t: %ld", &cores);
            return cores;
        }
    }

    return 0;
}

int get_blis(void)
{
    DIR* dir=opendir("./blis/blis/lib");
    if(dir)
    {
        closedir(dir);
        return 1;
    }
    else if(ENOENT == errno)
    {   
        return 0;
    }
    return 0;
}

void printHelp(void)
{
    printf("qLD-compute manual\n"
           "------------------\n"
           "\t-input        first_input_Folder\n"
           "\t-input2       second_input_Folder\n"
           "\t-output       output_File\n"
           "\t-ploidy       {haploid/phased_diploid/unphased_diploid}\n"
           "\t-posWmin1     snip_pos\n"
           "\t-posWmax1     snip_pos\n"
           "\t-posWmin2     snip_pos\n"
           "\t-posWmax2     snip_pos\n"
           "\t-sampleList   input_File\n"
           "\t-sampleList2  input_File\n"
           "\t-inputList    input_File\n"
           "\t-r2limit      value\n"
           "\t-threads      value\n"
           "\t-sorted\n"
           "\t-blis\n"
#ifdef GPU
           "\t-gpu\n"
#endif
           "\nDescription\n"
           "\t-input       <STRING>  Specifies the path of the first input alignment parsed files\n"
           "\t-input2      <STRING>  Specifies the path of the second input alignment parsed files\n"
           "\t                       Optional, uses Position Window 2 if set, else\n"
           "\t                       Position Window 1 is used\n"
           "\t-output      <STRING>  Specifies the name of the output alignment file.\n"
           "\t-ploidy      <STRING>  Supported ploidy types:\n"
           "\t                       haploid:           single digit snip: ex. '0'\n"
           "\t                       phased_diploid:    double digit snip: ex. \"0|0\"\n"
           "\t                       unphased_diploid:  double digit snip:\tex. \"0/0\"\n"
           "\t-sampleList  <STRING>  txt file with Format:\n"
           "\t                       \"sample1\n"
           "\t                        sample2\n"
           "\t                        sample3\n"
           "\t                        ...\n"
           "\t                        sampleN\"\n"
           "\t                       Specifies the name of the file that includes\n"
           "\t                       a list of valid samples from first input\n"
           "\t                       that will be selected for processing\n"
           "\t-sampleList2 <STRING>  Specifies the name of the file that includes\n"
           "\t                       a list of valid samples from second input\n"
           "\t                       that will be selected for processing\n"
           "\t                       Optional, uses seperate files for each input if set,\n"
           "\t                       else uses the file from -sampleList for both inputs\n\n"
           "\tPosition Window 1 (mutually exclusive with input list)\n"
           "\t------------------------------------------------------\n"
           "\t-posWmin1    <INT>     pos of the minimum snip to be included to window 1\n"
           "\t-posWmax1    <INT>     pos of the maximum snip to be included to window 1\n\n"
           "\tPosition Window 2 (mutually exclusive with input list)\n"
           "\tOptional, uses Input 2 if set, else Input 1 is used\n"
           "\t------------------------------------------------------\n"
           "\t-posWmin2    <INT>     pos of the minimum snip to be included to window 2\n"
           "\t-posWmax2    <INT>     pos of the maximum snip to be included to window 2\n\n"
           "\tInput list (mutually exclusive with inputs and position windows)\n"
           "\t----------------------------------------------------------------\n"
           "\t-inputList   <STRING>  csv file with Format:\n"
           "\t                       \"input,posWmin1,posWmin2,input2,posWmin2,posWmax2\"\n"
           "\t                       If single input is needed, duplicate the first 3 args.\n"
           "\t-sorted                (requires input list) Sorts Input List\n\n"
           "\t-r2limit     <FLOAT>   the lowest r2 value to be included in the results (default 0.2)\n"
           "\t-threads     <INT>     Number of threads to run in parallel.\n"
           "\t                       Suggested to use physical core number at max.\n"
           "\t                       On your system this would be %ld.\n", get_corecount());
        printf(
           "\t-blis                  Use the blis framework for calculations\n");
    if(get_blis())
        printf(
           "\t                       (Your System is eligible for blis)\n");
    else
        printf(
           "\t                       (Your System is not eligible for blis)\n");
#ifdef GPU
    printf("\t-gpu                   use the gpu for calculations\n"
           "\t                       (there has to be a gpu in the system)\n");
#endif
    printf("\t-mdf                   use preprocessed files as input\n"
           "\t                       (needs mdf parsing your input before using this feature)\n"
           "\t                       Recommended, perfomance greatly enhanced\n\n");
}

void commandLineParser(int argc,
        char** argv,
        char * inPath,
        char * inPath2,
        int * inPath2set,
        char * outfile,
        char **sampleListName,
        char **sampleListName2,
        int *ploidy,
        int * posWmin1,
        int * posWmax1,
        int *posWset1,
        int * posWmin2,
        int * posWmax2,
        int *posWset2,
        char * inList,
        int * inListSet,
        float * r2limit,
        int *gpu,
        int *blis,
        int *table,
        int *threads,
        int *sorted,
        int *mdf)
{
    int i, pathSet=0, fileSet=0, sampleListSet=0, sampleList2Set=0;
    if(argc < 2)
    {
        fprintf(stderr, "qLD-compute needs some arguments to run. Run \"qLD-compute -help\" for more info\n");
        exit(1);
    }
    if(strcmp(argv[1], "-help") &&
           strcmp(argv[1], "-h") &&
           strcmp(argv[1], "--help"))
    {
        printf("qLD-compute running with arguments:\n");
    }
    for(i=1; i < argc; ++i)
    {
        if(!strcmp(argv[i], "-help") ||
           !strcmp(argv[i], "-h") ||
           !strcmp(argv[i], "--help"))
        {
            printHelp();
            exit(0);
        }
        if(!strcmp(argv[i], "-input"))
        {
            if (i != argc-1)
            {
                printf("\t%s %s\n",argv[i],argv[i+1]);
                fflush(stdout);
                getDir(inPath,argv[++i]);
                DIR* directory=opendir(inPath);
                if (directory == NULL) {
                    fprintf(stderr, "\n ERROR: Directory %s does not exist.\n\n",inPath);
                    exit(1);
                }
                else {
#ifdef VERBOSE
                    printf( "Input Directory %s opened.\n",inPath);
#endif
                    closedir(directory);
                }
                pathSet=1;
            }
            continue;
        }
        if(!strcmp(argv[i], "-input2"))
        {
            if (i != argc-1)
            {
                printf("\t%s %s\n",argv[i],argv[i+1]);
                fflush(stdout);
                getDir(inPath2,argv[++i]);
                DIR* directory=opendir(inPath2);

                if (directory == NULL)
                {
                    fprintf(stderr, "\n ERROR: Directory %s does not exist.\n\n",inPath2);
                    exit(1);
                }
                else
                {
#ifdef VERBOSE
                    printf( "Input 2 Directory %s opened.\n",inPath2);
#endif
                    closedir(directory);
                }
                *inPath2set=1;
            }
            continue;
        }
        if(!strcmp(argv[i], "-inputList"))
        {
            printf("\t%s %s\n",argv[i],argv[i+1]);
            fflush(stdout);
            *inListSet=1;
            strcpy(inList,argv[++i]);
            continue;
        }
        if(!strcmp(argv[i], "-output"))
        {
            if(i != argc-1)
            {
                printf("\t%s %s\n",argv[i],argv[i+1]);
                fflush(stdout);
                strcpy(outfile,argv[++i]);
                fileSet=1;
            }
            continue;
        }

        if(!strcmp(argv[i], "-sampleList"))
        {
            sampleListSet=1;
            char *temp=(char *)malloc(sizeof(char)*INFILENAMESIZE);
            assert(temp);
            (*sampleListName)=temp;
            assert((*sampleListName));
            if(i != argc-1)
            {
                printf("\t%s %s\n",argv[i],argv[i+1]);
                fflush(stdout);
                strcpy((*sampleListName), argv[++i]);
            }
            continue;
        }
        if(!strcmp(argv[i], "-sampleList2"))
        {
            sampleList2Set=1;
            char *temp=(char *)malloc(sizeof(char)*INFILENAMESIZE);
            assert(temp);
            (*sampleListName2)=temp;
            assert((*sampleListName2));
            if(i != argc-1)
            {
                printf("\t%s %s\n",argv[i],argv[i+1]);
                fflush(stdout);
                strcpy((*sampleListName2), argv[++i]);
            }
            continue;
        }
        if(!strcmp(argv[i], "-ploidy"))
        {
            printf("\t%s %s\n",argv[i],argv[i+1]);
            fflush(stdout);
            char p_val[20];
            strcpy(p_val,argv[++i]);
            if(!strcmp(p_val, "haploid"))
                *ploidy=1;
            else if(!strcmp(p_val, "phased_diploid"))
                *ploidy=2;
            else if(!strcmp(p_val, "unphased_diploid"))
                *ploidy=2;
            else
            {
                fprintf(stderr,"Please insert supported ploidy (see --help for more info)\n");
                exit(1);
            }
            continue;
        }
        if(!strcmp(argv[i], "-posWmin1"))
        {
            if (i != argc-1)
            {
                printf("\t%s %s\n",argv[i],argv[i+1]);
                fflush(stdout);
                *posWmin1=atoi(argv[++i]);
                *posWset1=1;
            }
            continue;
        }
        if(!strcmp(argv[i], "-posWmax1"))
        {
            if (i != argc-1)
            {
                printf("\t%s %s\n",argv[i],argv[i+1]);
                fflush(stdout);
                *posWmax1=atoi(argv[++i]);
                *posWset1=1;
            }
            continue;
        }
        if(!strcmp(argv[i], "-posWmin2"))
        {
            if (i != argc-1)
            {
                printf("\t%s %s\n",argv[i],argv[i+1]);
                fflush(stdout);
                *posWmin2=atoi(argv[++i]);
                *posWset2=1;
            }
            continue;
        }
        if(!strcmp(argv[i], "-posWmax2"))
        {
            if (i != argc-1)
            {
                printf("\t%s %s\n",argv[i],argv[i+1]);
                fflush(stdout);
                *posWmax2=atoi(argv[++i]);
                *posWset2=1;
            }
            continue;
        }
        if(!strcmp(argv[i], "-r2limit"))
        {
            if (i != argc-1)
            {
                printf("\t%s %s\n",argv[i],argv[i+1]);
                fflush(stdout);
                *r2limit=atof(argv[++i]);
            }
            continue;
        }
        if(!strcmp(argv[i], "-threads"))
        {
            if(i != argc-1)
            {
                printf("\t%s %s\n",argv[i],argv[i+1]);
                fflush(stdout);
                *threads=atoi(argv[++i]);
            }
            continue;
        }
        if(!strcmp(argv[i], "-sorted"))
        {
            printf("\t%s\n",argv[i]);
            fflush(stdout);
            *sorted=1;
            continue;
        }
        if(!strcmp(argv[i], "-rTable"))
        {
            printf("\t%s\n",argv[i]);
            fflush(stdout);
            *table = 1;
            continue;
        } 
#ifdef CBLAS_USE
        if(!strcmp(argv[i], "-blis"))
        {
            printf("\t%s\n",argv[i]);
            fflush(stdout);
            *blis=1;
            continue;
        }
#else
        *blis=0;
#endif
        //if flag GPU is given to the makefile and the code is run with -gpu
#ifdef GPU
        if(!strcmp(argv[i], "-gpu"))
        {
            printf("\t%s\n",argv[i]);
            fflush(stdout);
            *gpu=1;
            continue;
        }
#else
        *gpu=0;
#endif
        if(!strcmp(argv[i], "-mdf"))
        {
            printf("\t%s\n",argv[i]);
            fflush(stdout);
            *mdf=1;
            continue;
        }

        fprintf(stderr, "\nERROR: %s is not a valid command line parameter\n\n",argv[i]);
        fflush(stderr);
        exit(1);
    }
    if(!(*ploidy))
    {
        fprintf(stderr,"Please insert supported ploidy (see --help for more info)\n");
        exit(1);
    }
    if((*inListSet) == 1)
    {
        pathSet=1;
        *posWset1=1;
        *inPath2set=1;
        *posWset2=1;
    }
    if(pathSet == 0)
    {
        fprintf(stderr,"\n ERROR: Please specify a path for the input with -input\n\n");
        exit(1);
    }
    if(fileSet == 0)
    {
        fprintf(stderr, "\n ERROR: Please specify an alignment folder with -output\n\n");
        exit(1);
    }
    if(sampleListSet == 1 && sampleList2Set == 0)
    {
        char *temp=(char *)malloc(sizeof(char)*INFILENAMESIZE);
        assert(temp);
        (*sampleListName2)=temp;
        assert((*sampleListName2));
        strcpy((*sampleListName2),(*sampleListName));
    }
    if(sampleListSet == 0 && sampleList2Set == 1)
    {
        fprintf(stderr,"\n ERROR: Please set sampleList before setting sampleList2\n\n");
        exit(1);
    }
}

int main(int argc, char** argv)
{
    // Saving project path, in case something changes along the way
    t_head=NULL;
    t_tail=NULL;

    int inPath2set=0, posWmin1=0, posWmax1=0, posWset1=0;
    int posWmin2=0, posWmax2=0, posWset2=0;
    int inListSet=0, gpu=0, blis=0, rTable=0, threads=1, sorted=0, mdf=0;
    int ploidy=0;
    float r2limit = 0.2;
    char inputPathName[INFILENAMESIZE]="";
    char inputPath2Name[INFILENAMESIZE]="";
    char inputListName[INFILENAMESIZE]="";
    char outputFileName[INFILENAMESIZE]="";
    char *sampleListName=NULL;
    char *sampleListName2=NULL;

    commandLineParser(argc,
            argv,
            inputPathName,
            inputPath2Name,
            &inPath2set,
            outputFileName,
            &sampleListName,
            &sampleListName2,
            &ploidy,
            &posWmin1,
            &posWmax1,
            &posWset1,
            &posWmin2,
            &posWmax2,
            &posWset2,
            inputListName,
            &inListSet,
            &r2limit,
            &gpu,
            &blis,
            &rTable,
            &threads,
            &sorted,
            &mdf);

//#if defined(VERBOSE) || defined(BENCHMARK)
//    double prep_time_s=gettime();
    if(blis)
    {
        printf("Blis kernel active\n");
    }
    if(gpu)
    {
        printf("GPU kernel active\n");
    }
    if(mdf)
    {
        printf("MDF input active\n");
    }
//#endif
    sample_t *sampleList=create_sample_t(BUFFERSIZE,ploidy);
    sample_t *sampleList2=create_sample_t(BUFFERSIZE,ploidy);
    getFilter(sampleListName, sampleList);
    getFilter(sampleListName2, sampleList2);
    if(sampleListName!=NULL)
        free(sampleListName);
    if(sampleListName2!=NULL)
        free(sampleListName2);

    if(inListSet == 1)
    {
        create_task_queue(inputListName,
                          outputFileName,
                          sampleList,
                          sampleList2,
                          inPath2set,
                          posWset1,
                          posWset2,
                          mdf);
    }
    else
    {
        preprocess_data(inputPathName,
                        inputPath2Name,
                        outputFileName,
                        sampleList,
                        sampleList2,
                        inPath2set,
                        posWset1,
                        posWset2,
                        posWmin1,
                        posWmax1,
                        posWmin2,
                        posWmax2,
                        mdf);
    }
    fflush(stdout);

    if(sorted && inListSet)
        MergeSort(&t_head);

#ifdef GPU
    if(gpu)
        gpu_init();
#endif

//#if defined(VERBOSE) || defined(BENCHMARK)
//    double prep_time_e=gettime()-prep_time_s;
    double timeStart=gettime();
//#endif
    sense_reversal_barrier_init(threads);
    assert(threads > 0);
    static pthread_t *workerThread=NULL;
    if(threads > 1)
        workerThread=(pthread_t *)malloc(sizeof(pthread_t)*((unsigned long)(threads-1)));

    threadData_t *threadData=(threadData_t *)malloc(sizeof(threadData_t)*
            ((unsigned long)threads));
    assert(threadData!=NULL);

    for(int i=0; i < threads; i++)
    {
        initializeThreadData(&threadData[i], i, threads);
    }
    for(int i=1; i < threads; i++)
    {
        pthread_create(&workerThread[i-1],NULL,thread,(void *)(&threadData[i]));
    }

    updateThreadArgs(&threadData[0], r2limit, ploidy, gpu, blis, mdf);
    startThreadOPS(threadData, CORRELATE);
    terminateWorkerThreads(workerThread,threadData);
    sense_reversal_barrier_destroy();

#ifdef GPU
    if(gpu)
        gpu_release();
#endif

    free_sample_t(sampleList);
    free_sample_t(sampleList2);
    free(workerThread);
    if(threadData!=NULL)
        free(threadData);
    threadData=NULL;

//#if defined(VERBOSE) || defined(BENCHMARK)
    double timeTotal=gettime()-timeStart;
//    printf("------------------------\nPreprocess time: %.3fs\n", prep_time_e);
    printf("Finished in: %.3fs\n", timeTotal);
//#endif
    printf("Output saved in: %s\n", dirname(outputFileName));
    return 0;
}
