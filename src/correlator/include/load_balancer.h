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
#ifndef LB_H
#define LB_H

#include "input_list.h"
#include "correlate_MEM.h"

t_node* SortedMerge(t_node* a, t_node* b);

void FrontBackSplit(t_node* source, t_node** frontRef,
        t_node** backRef);

void MergeSort(t_node** head);

void preprocess_data(FILE *fpInRep,
        char *inputPathName,
        char *inputPath2Name,
        char *outputFileName,
        sample_t *sampleList,
        sample_t *sampleList2,
        int inPath2set,
        int posWset1,
        int posWset2,
        int posWmin1,
        int posWmax1,
        int posWmin2,
        int posWmax2,
        int mdf,
        int *task_count);

#endif
