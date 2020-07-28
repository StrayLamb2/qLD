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
#ifndef IL_H
#define IL_H

#include "header.h"
#include "correlate_MEM.h"

typedef struct node_s{
    char input_file[INFILENAMESIZE];
    char input_file2[INFILENAMESIZE];
    char output_file[INFILENAMESIZE];
    int id, posWmin, posWmax, posWmin2, posWmax2;
    int inPath2set, posWset1, posWset2;
    char **filesList;
    int filesListNum;
    char **filesList2;
    int filesListNum2;
    char *headerLine1;
    char *headerLine2;
    char *headerLine1_2;
    char *headerLine2_2;
    sample_t *sampleList;
    sample_t *sampleList2;
    int snipSize;
    int snipSize2;
    struct node_s *next;
} t_node;

EXTERN_INP t_node *t_head, *t_tail;

void enqueue_task(char* input,
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
                  sample_t *sampleList2);

void free_task(t_node *task);

t_node* dequeue_task(void);
t_node* get_task(int id);

int create_task_queue(char *inp_list,
                      char *output_file,
                      sample_t * sampleList,
                      sample_t * sampleList2,
                      int inPath2set,
                      int posWset1,
                      int posWset2,
                      int mdf,
                      int *task_count);

void strip_ext(char *fname);
#endif
