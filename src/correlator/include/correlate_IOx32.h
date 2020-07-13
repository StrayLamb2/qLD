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
#ifndef CIO32_H
#define CIO32_H

#include "header.h"
#include "correlate_MEM.h"
#include "pthreads.h"

int writeMDF_x32(table_x32 *tableData,
                 helper_t *helperData);

int writeToTable_x32(task_t *nodeData,
                     table_x32 *tableData,
                     helper_t *helperData,
                     int pos,
                     int *skippedLines,
                     int *counter);

void readTable_x32(threadData_t *threadData,
                   task_t *nodeData,
                   table_x32 *tableData,
                   helper_t *helperData);
#endif
