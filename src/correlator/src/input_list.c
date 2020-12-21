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
#define EXTERN_INP

#include "../include/header.h"
#include "../include/input_list.h"
#include "../include/correlate_IO.h"
#include "../include/load_balancer.h"

t_node *t_head, *t_tail;

// Snippet to strip extention of a file
void strip_ext(char *fname)
{
    char *end=NULL;
    end=fname + strlen(fname);

    while (end > fname && *end != '.' && *end != '\\' && *end != '/') {
        --end;
    }
    if ((end > fname && *end == '.') &&
        (*(end - 1) != '\\' && *(end - 1) != '/')) {
        *end = '\0';
    }
}

void report_input(FILE *fpInRep, 
                  int id, 
                  int posWmin, 
                  int posWmax, 
                  int posWmin2, 
                  int posWmax2,
                  int snipSize,
                  int snipSize2)
{
    int sites1=posWmax-posWmin;
    int sites2=posWmax2-posWmin2;
    
    fprintf(fpInRep, "%d)\n"
                     "\tSamples:  %10d\t%10d\n"
                     "\tSites:    %10d\t%10d\n",
                     id+1, snipSize, snipSize2, sites1, sites2);
    fflush(fpInRep);
}

// Function to create and enqueue tasks
void enqueue_task(FILE *fpInRep,
                  char* input,
                  int posWmin,
                  int posWmax,
                  char* input2,
                  int posWmin2,
                  int posWmax2,
                  char* output,
                  int inPath2set,
                  int posWset1,
                  int posWset2,
                  char **filesList,
                  int filesListNum,
                  char **filesList2,
                  int filesListNum2,
                  char *headerLine1,
                  char *headerLine2,
                  char *headerLine1_2,
                  char *headerLine2_2,
                  int snipSize,
                  int snipSize2,
                  sample_t *sampleList,
                  sample_t *sampleList2)
{
    static int id;
    if(!t_head)
    {
        id=0;
        t_head=(t_node *)malloc(sizeof(t_node));
        t_tail=t_head;
    }
    else
    {
        id++;
        t_tail->next=(t_node *)malloc(sizeof(t_node));
        t_tail=t_tail->next;
    }

    t_tail->id=id;
    strcpy(t_tail->input_file,input);
    t_tail->posWmin=posWmin;
    t_tail->posWmax=posWmax;
    strcpy(t_tail->input_file2,input2);
    t_tail->posWmin2=posWmin2;
    t_tail->posWmax2=posWmax2;
    strcpy(t_tail->output_file,output);
    t_tail->inPath2set=inPath2set;
    t_tail->posWset1=posWset1;
    t_tail->posWset2=posWset2;
    t_tail->headerLine1=headerLine1;
    t_tail->headerLine2=headerLine2;
    t_tail->headerLine1_2=headerLine1_2;
    t_tail->headerLine2_2=headerLine2_2;
    t_tail->snipSize=snipSize;
    t_tail->snipSize2=snipSize2;
    t_tail->sampleList=sampleList;
    t_tail->sampleList2=sampleList2;

    t_tail->filesList=malloc(filesListNum*sizeof(char*));
    for(int i=0; i < filesListNum; i++)
        t_tail->filesList[i]=malloc(INFILENAMESIZE*sizeof(char));
    for(int i=0; i < filesListNum; i++)
        strcpy(t_tail->filesList[i],filesList[i]);
    t_tail->filesListNum=filesListNum;

    t_tail->filesList2=malloc(filesListNum2*sizeof(char*));
    for(int i=0; i < filesListNum2; i++)
        t_tail->filesList2[i]=malloc(INFILENAMESIZE*sizeof(char));
    for(int i=0; i < filesListNum2; i++)
        strcpy(t_tail->filesList2[i],filesList2[i]);
    t_tail->filesListNum2=filesListNum2;

    t_tail->next=NULL;

    report_input(fpInRep, 
                 id, 
                 posWmin, 
                 posWmax, 
                 posWmin2, 
                 posWmax2,
                 snipSize,
                 snipSize2);
}

void free_task(t_node *task)
{
    int i;

    free(task->headerLine1);
    free(task->headerLine2);
    free(task->headerLine1_2);
    free(task->headerLine2_2);

    for(i=0; i < task->filesListNum; i++)
            free(task->filesList[i]);
    free(task->filesList);
    for(i=0; i < task->filesListNum2; i++)
            free(task->filesList2[i]);
    free(task->filesList2);

    free(task);
}

// Function to get a task from the queue
t_node* dequeue_task(void)
{
    if(t_head)
    {
        t_node* task=t_head;
        t_head=t_head->next;
        return task;
    }
    t_head=NULL;
    t_tail=NULL;

    return NULL;
}

t_node* get_task(int id)
{
    if(!t_head)
        return NULL;
    t_node* task=t_head;
    if(t_head->id == id)
    {
        task=t_head;
        t_head=t_head->next;
        return task;
    }
    if(!task->next)
        return NULL;
    while(task->next->id != id)
    {
        task=task->next;
        if(!task->next)
            return NULL;
    }
    t_node* node=task->next;
    task->next=task->next->next;
    return node;
}

// Function to get and validate task from the input list, to insert in the queue
int create_task_queue(FILE *fpInRep,
                      char *inp_list,
                      char *output_file,
                      sample_t *sampleList,
                      sample_t *sampleList2,
                      int inPath2set,
                      int posWset1,
                      int posWset2,
                      int mdf,
                      int *task_count)
{
    FILE *fp;
    int ret, lineNo=0, posWmin1, posWmax1, posWmin2, posWmax2, taskno=0;

    char input[INFILENAMESIZE];
    char posWmin1_s[INFILENAMESIZE];
    char posWmax1_s[INFILENAMESIZE];

    char input2[INFILENAMESIZE];
    char posWmin2_s[INFILENAMESIZE];
    char posWmax2_s[INFILENAMESIZE];

    char out_temp[INFILENAMESIZE];
    char output[INFILENAMESIZE];

    char line[STRINGLENGTH];

    fp=fopen(inp_list, "r");
    if(!fp)
    {
        fprintf(stderr, "ERROR: Input list file not existing!\n");
        exit(1);
    }
    // Strip extention of output, to inject the line number to the output name
    strip_ext(output_file);

    // Read input list line by line, where each line is a task
    while(++lineNo)
    {
        if(fgets(line, STRINGLENGTH, fp) == NULL)
            break;
        // Input list is a csv with input,posWmin1,posWmax1,input2,posWmin2,posWmax2
        // as values
        ret=sscanf(line, "%[^,],%[^,],%[^,],%[^,],%[^,],%[^\n]",
                   input, posWmin1_s, posWmax1_s, input2, posWmin2_s, posWmax2_s);

        // If scanning failed to get all the args
        if(ret!=6)
        {
            printf("ABORTING: Line %d ---Invalid structure---\n", lineNo);
            continue;
        }

        // Get input in correct form
        getDir(input, input);
        getDir(input2, input2);

        posWmin1=atoi(posWmin1_s);
        posWmax1=atoi(posWmax1_s);
        if(posWmin1>=posWmax1)
        {
            printf("ABORTING: Line %d ---posWmin greater than posWmax---\n", lineNo);
            continue;
        }

        posWmin2=atoi(posWmin2_s);
        posWmax2=atoi(posWmax2_s);
        if(posWmin2>=posWmax2)
        {
            printf("ABORTING: Line %d ---posWmin2 greater than posWmax2---\n", lineNo);
            continue;
        }

        // Manipulate output to include line number
        sprintf(out_temp,"%d",lineNo);
        strcpy(output, output_file);
        strcat(output, "_");
        strcat(output, out_temp);
        strcat(output, ".txt");

        // Check if input files exist, drop the task if not
        if((access(input, F_OK) != -1) && (access(input2, F_OK) != -1))
        {
            taskno++;
            preprocess_data(fpInRep,
                            input,
                            input2,
                            output,
                            sampleList,
                            sampleList2,
                            inPath2set,
                            posWset1,
                            posWset2,
                            posWmin1,
                            posWmax1,
                            posWmin2,
                            posWmax2,
                            mdf,
                            task_count);

        }
        else
        {
            printf("ABORTING: Line %d ---Input files do not exist---\n",lineNo);
            continue;
        }
    }
    fclose(fp);

    if(!t_head)
    {
        fprintf(stderr,"ERROR: All files in the list are invalid!\n");
        exit(1);
    }
#ifdef VERBOSE
    printf("\n");
#endif
    return 0;
}
