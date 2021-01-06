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
#ifndef GG_H
#define GG_H

#ifdef GPU
#include "header.h"
#include "pthreads.h"

#define PROGRAM_FILE "bin/gpu_kernel/blislike.cl"
#define KERNEL_NAME "blis_like8x4"
#define RESULTS_PART_SIZE_GPU 4 //2*2

// macro define needed for deprecated warning
#define CL_TARGET_OPENCL_VERSION 120
#define CL_USE_DEPRECATED_OPENCL_1_2_APIS

#include <CL/cl.h>
#define DOUBLE unsigned int

EXTERN_GPU cl_platform_id *platforms;
EXTERN_GPU cl_device_id *devices;
EXTERN_GPU cl_context context;
EXTERN_GPU cl_program program;
EXTERN_GPU cl_mem    a_buffers[2];
EXTERN_GPU cl_mem    b_buffers[2];
EXTERN_GPU cl_mem    c_buffers[2];
EXTERN_GPU cl_kernel kernels[2];
EXTERN_GPU cl_command_queue io_queue;
EXTERN_GPU cl_command_queue compute_queue;
EXTERN_GPU cl_event events[4*2];
EXTERN_GPU unsigned int rs_c;
// NOTE: since we're writing into a smaller buffer and not the full output
// matrix, we mult by MC to get to the next row instead of the full "m".
EXTERN_GPU unsigned int cs_c;

void gpu_init(void);

void gpu_release(void);

void correlate_gpu(threadData_t *threadData,
                   outFileType fpOut,
                   inputDataType_x32* tableA,
                   unsigned int *tableA_posIndex,
                   char **tableA_IDindex,
                   unsigned int *tableA_bitcount,
                   int tableAsize,
                   inputDataType_x32* tableB,
                   unsigned int *tableB_posIndex,
                   char **tableB_IDindex,
                   unsigned int *tableB_bitcount,
                   int tableBsize,
                   int compressed_snp_size,
                   int snp_size,
                   int posWset2);

void printCLErr(cl_int err,int line, char* file);

void GPU_Pack_A(inputDataType_x32 *A,
                unsigned int lda,
                DOUBLE *A_pack,
                unsigned int m,
                unsigned int k);

void GPU_Pack_B(inputDataType_x32 *B,
                unsigned int ldb,
                DOUBLE *B_pack,
                unsigned int k,
                unsigned int n);

void gpu_gemm(unsigned int m,
              unsigned int n,
              unsigned int k,
              inputDataType_x32 *A,
              unsigned int lda,
              inputDataType_x32 *B,
              unsigned int ldb,
              inputDataType_x32 * C,
              unsigned int ldc,
              void * Ac_pack_v,
              void * Bc_pack_v);

void get_pairwise_ld_score_gpu(unsigned int * tableA_bitcount,
                           unsigned int * tableB_bitcount,
                           inputDataType_x32 * C,
                           int tableAsize,
                           int tableBsize,
                           int snp_size,
                           ResultDataType** results);

void mlt_gpu(unsigned int m,
         unsigned int k,
         inputDataType_x32 *A,
         inputDataType_x32* tableA);
#endif
#endif
