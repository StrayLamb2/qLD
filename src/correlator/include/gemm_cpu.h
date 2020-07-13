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
#ifndef GC_H
#define GC_H

#include "header.h"
#include "pthreads.h"

void correlate(threadData_t *threadData,
               outFileType fpOut,
               inputDataType_x64* tableA,
               unsigned int *tableA_posIndex,
               char **tableA_IDindex,
               unsigned int *tableA_bitcount,
               int tableAsize,
               inputDataType_x64* tableB,
               unsigned int *tableB_posIndex,
               char **tableB_IDindex,
               unsigned int *tableB_bitcount,
               int tableBsize,
               int compressed_snp_size,
               int snp_size,
               int posWset2);

void Pack_A(inputDataType_x64 *A,
            unsigned int lda,
            inputDataType_x64 *A_pack,
            unsigned int m,
            unsigned int k);

void Pack_B(inputDataType_x64 *B,
            unsigned int ldb,
            inputDataType_x64 *B_pack,
            unsigned int k,
            unsigned int n);

void dgemm_ref(int k,
               int mr_alg,
               int nr_alg,
               inputDataType_x64 alpha,
               inputDataType_x64* a,
               inputDataType_x64* b,
               inputDataType_x64* c,
               int rs_c,
               int cs_c);

void blis_gemm(unsigned int M,
               unsigned int N,
               unsigned int K,
               inputDataType_x64 alphap,
               inputDataType_x64 *a,
               inputDataType_x64 *b,
               inputDataType_x64 *c);

void gemm(unsigned int m,
          unsigned int n,
          unsigned int k,
		  inputDataType_x64 alphap,
          inputDataType_x64 * A,
          unsigned int lda,
          inputDataType_x64 * B,
          unsigned int ldb,
          inputDataType_x64 * C,
          unsigned int ldc,
          void * Ac_pack_v,
          void * Bc_pack_v);

void get_pairwise_ld_score(unsigned int * tableA_bitcount,
                           unsigned int * tableB_bitcount,
                           inputDataType_x64 * C,
                           int tableAsize,
                           int tableBsize,
                           int snp_size,
                           ResultDataType** results);

void mlt(unsigned int m,
         unsigned int k,
         inputDataType_x64* A,
         inputDataType_x64* tableA);
#endif
