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
#ifndef CMEM_H
#define CMEM_H

#include "header.h"

void *realloc_buff(void *buffer, unsigned int newSize, char type[32]);

void print_alloc_size(char table[32], unsigned int size);

typedef struct sample_t{
    int sampleListSize;
    int sampleListIndex;
    char **sampleList;
    int ploidy;
} sample_t;

typedef struct valid_t{
    int validListSize;
    int valid_count;
    int *validList;
} valid_t;

typedef struct rw_task_s{
    int id;
    unsigned int tableSize;
    char **filesList;
    int filesListNum;
    int snipSize;
    char *headerLine1;
    char *headerLine2;
    int posWmin;
    int posWmax;
    int posWset;
    int *valid_mask;
    int valid_count;
} task_t;

typedef struct table_x32{
    inputDataType_x32 *SNPtable;
    unsigned int *POStable;
    unsigned int *BCtable;
    char **IDtable;
    unsigned int tableSize;
    unsigned int tableIndex;
    unsigned int compSize;
    unsigned int compIndex;
}table_x32;

typedef struct table_x64{
    inputDataType_x64 *SNPtable;
    unsigned int *POStable;
    unsigned int *BCtable;
    char **IDtable;
    unsigned int tableSize;
    unsigned int tableIndex;
    unsigned int compSize;
    unsigned int compIndex;
}table_x64;

typedef struct helper_s{
    int ploidy;
    char *snipChar;
    int snipCharLength;
    int snipCharIndex;
    char *snipPart;
    char *line;
    int lineLength;
    char *word;
    int wordLength;
}helper_t;

sample_t* create_sample_t(int sampleListSize, int ploidy);
void free_sample_t(sample_t* sampleList);

valid_t *create_valid_t(int validSize);
void free_valid_t(valid_t *validdata);

task_t* create_task_t(int ploidy);
void free_task_t(task_t *nodeData);

table_x32* create_table_x32(unsigned int tableSize);
void free_table_x32(table_x32 *tableData);

table_x64* create_table_x64(unsigned int tableSize);
void free_table_x64(table_x64 *tableData);

helper_t* create_helper_t(int ploidy);
void free_helper_t(helper_t *helperData);

#endif
