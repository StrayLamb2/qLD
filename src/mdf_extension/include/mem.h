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
#include "structs.h"

void *realloc_buff(void *buffer, unsigned int newSize, char type[32]);

void print_alloc_size(char table[32], unsigned int size);

table_x64* create_table_x64(unsigned int tableSize);
void free_table_x64(table_x64 *tableData);

helper_t* create_helper_t(int ploidy);
void free_helper_t(helper_t *helperData);

#endif
