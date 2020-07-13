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
#ifndef CIO_H
#define CIO_H

#include "header.h"
#include "structs.h"

int getValueFromFilename(const void* filepath);

int myCompare(const void* a, const void* b);

void sort_files(char** arr, int n);

header_t *init_header_struct(void);

void print_header_struct(header_t *headerData);

void free_header_struct(header_t *headerData);

char *getDir(char *output, char *input);

void getFilter(char *filterFile, char ***wordList, int *wordListSize);

int sample_isValid(char **list, int list_size, char *sample, int *valid_count);

void readHeaderFile(pre_t *preData, header_t *headerData);

void findFiles(pre_t *preData, header_t *headerData);

void writeResultsHeader(outFileType fpOut);

void writeResults(outFileType fpOut,
                  unsigned int *table_A_posIndex,
                  char **tableA_IDindex,
                  unsigned int *tableA_bitcount,
                  int tableASize,
                  unsigned int *table_B_posIndex,
                  char **tableB_IDindex,
                  unsigned int *tableB_bitcount,
                  int tableBSize,
                  int snp_size,
                  ResultDataType *results,
                  float r2limit,
                  int posWset2);
#endif
