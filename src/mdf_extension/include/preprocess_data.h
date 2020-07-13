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

#include "header.h"
#include "structs.h"

pre_t *init_preprocess_struct(void);

void print_preprocess_struct(pre_t *preData);

void free_preprocess_struct(pre_t *preData);

void preprocess_data(pre_t *preData, header_t *headerData);
#endif
