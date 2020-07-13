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
#ifdef GPU
#define EXTERN_GPU

#include "../include/header.h"
#include "../include/gemm_gpu.h"
#include "../include/correlate_IO.h"
#include "../include/correlate_IOx32.h"
#include "../include/snip_processing.h"
#include "../include/read_file.h"


static unsigned int round_up_mult(unsigned int x, unsigned int mult) {
    return ((x + mult - 1) / mult) * mult;
}

void printCLErr(cl_int err,int line, char* file)
{
    if(err == CL_SUCCESS) {
        return;
    }
    switch (err)
    {
        case CL_SUCCESS:
            printf("CL_SUCCESS\n");
            break;
        case CL_DEVICE_NOT_FOUND:
            printf("CL_DEVICE_NOT_FOUND\n");
            break;
        case CL_DEVICE_NOT_AVAILABLE:
            printf("CL_DEVICE_NOT_AVAILABLE\n");
            break;
        case CL_COMPILER_NOT_AVAILABLE:
            printf("CL_COMPILER_NOT_AVAILABLE\n");
            break;
        case CL_MEM_OBJECT_ALLOCATION_FAILURE:
            printf("CL_MEM_OBJECT_ALLOCATION_FAILURE\n");
            break;
        case CL_OUT_OF_RESOURCES:
            printf("CL_OUT_OF_RESOURCES\n");
            break;
        case CL_OUT_OF_HOST_MEMORY:
            printf("CL_OUT_OF_HOST_MEMORY\n");
            break;
        case CL_PROFILING_INFO_NOT_AVAILABLE:
            printf("CL_PROFILING_INFO_NOT_AVAILABLE\n");
            break;
        case CL_MEM_COPY_OVERLAP:
            printf("CL_MEM_COPY_OVERLAP\n");
            break;
        case CL_IMAGE_FORMAT_MISMATCH:
            printf("CL_IMAGE_FORMAT_MISMATCH\n");
            break;
        case CL_IMAGE_FORMAT_NOT_SUPPORTED:
            printf("CL_IMAGE_FORMAT_NOT_SUPPORTED\n");
            break;
        case CL_BUILD_PROGRAM_FAILURE:
            printf("CL_BUILD_PROGRAM_FAILURE\n");
            break;
        case CL_MAP_FAILURE:
            printf("CL_MAP_FAILURE\n");
            break;
        case CL_INVALID_VALUE:
            printf("CL_INVALID_VALUE\n");
            break;
        case CL_INVALID_DEVICE_TYPE:
            printf("CL_INVALID_DEVICE_TYPE\n");
            break;
        case CL_INVALID_PLATFORM:
            printf("CL_INVALID_PLATFORM\n");
            break;
        case CL_INVALID_DEVICE:
            printf("CL_INVALID_DEVICE\n");
            break;
        case CL_INVALID_CONTEXT:
            printf("CL_INVALID_CONTEXT\n");
            break;
        case CL_INVALID_QUEUE_PROPERTIES:
            printf("CL_INVALID_QUEUE_PROPERTIES\n");
            break;
        case CL_INVALID_COMMAND_QUEUE:
            printf("CL_INVALID_COMMAND_QUEUE\n");
            break;
        case CL_INVALID_HOST_PTR:
            printf("CL_INVALID_HOST_PTR\n");
            break;
        case CL_INVALID_MEM_OBJECT:
            printf("CL_INVALID_MEM_OBJECT\n");
            break;
        case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR:
            printf("CL_INVALID_IMAGE_FORMAT_DESCRIPTOR\n");
            break;
        case CL_INVALID_IMAGE_SIZE:
            printf("CL_INVALID_IMAGE_SIZE\n");
            break;
        case CL_INVALID_SAMPLER:
            printf("CL_INVALID_SAMPLER\n");
            break;
        case CL_INVALID_BINARY:
            printf("CL_INVALID_BINARY\n");
            break;
        case CL_INVALID_BUILD_OPTIONS:
            printf("CL_INVALID_BUILD_OPTIONS\n");
            break;
        case CL_INVALID_PROGRAM:
            printf("CL_INVALID_PROGRAM\n");
            break;
        case CL_INVALID_PROGRAM_EXECUTABLE:
            printf("CL_INVALID_PROGRAM_EXECUTABLE\n");
            break;
        case CL_INVALID_KERNEL_NAME:
            printf("CL_INVALID_KERNEL_NAME\n");
            break;
        case CL_INVALID_KERNEL_DEFINITION:
            printf("CL_INVALID_KERNEL_DEFINITION\n");
            break;
        case CL_INVALID_KERNEL:
            printf("CL_INVALID_KERNEL\n");
            break;
        case CL_INVALID_ARG_INDEX:
            printf("CL_INVALID_ARG_INDEX\n");
            break;
        case CL_INVALID_ARG_VALUE:
            printf("CL_INVALID_ARG_VALUE\n");
            break;
        case CL_INVALID_ARG_SIZE:
            printf("CL_INVALID_ARG_SIZE\n");
            break;
        case CL_INVALID_KERNEL_ARGS:
            printf("CL_INVALID_KERNEL_ARGS\n");
            break;
        case CL_INVALID_WORK_DIMENSION:
            printf("CL_INVALID_WORK_DIMENSION\n");
            break;
        case CL_INVALID_WORK_GROUP_SIZE:
            printf("CL_INVALID_WORK_GROUP_SIZE\n");
            break;
        case CL_INVALID_WORK_ITEM_SIZE:
            printf("CL_INVALID_WORK_ITEM_SIZE\n");
            break;
        case CL_INVALID_GLOBAL_OFFSET:
            printf("CL_INVALID_GLOBAL_OFFSET\n");
            break;
        case CL_INVALID_EVENT_WAIT_LIST:
            printf("CL_INVALID_EVENT_WAIT_LIST\n");
            break;
        case CL_INVALID_EVENT:
            printf("CL_INVALID_EVENT\n");
            break;
        case CL_INVALID_OPERATION:
            printf("CL_INVALID_OPERATION\n");
            break;
        case CL_INVALID_GL_OBJECT:
            printf("CL_INVALID_GL_OBJECT\n");
            break;
        case CL_INVALID_BUFFER_SIZE:
            printf("CL_INVALID_BUFFER_SIZE\n");
            break;
        case CL_INVALID_MIP_LEVEL:
            printf("CL_INVALID_MIP_LEVEL\n");
            break;
        case CL_INVALID_GLOBAL_WORK_SIZE:
            printf("CL_INVALID_GLOBAL_WORK_SIZE\n");
            break;
        default:
            printf("don't know what error that was\n");
            break;
    }
    printf ("Blah error at %s (%d)\n", file, line);
}

static void create_program_with_source(cl_program *program,
                                       cl_context *context,
                                       const char *program_file)
{
    // read in program source file and create `program`
    FILE *program_handle;
    size_t program_size;
    char *program_buffer;
    int err;

    program_handle=fopen(program_file, "r");
    //assert(program_handle);
    fseek(program_handle, 0, SEEK_END);
    program_size=ftell(program_handle);
    rewind(program_handle);

    program_buffer=(char*) malloc(program_size+1);
    assert(program_buffer);
    program_buffer[program_size]='\0';

    size_t ret_code=fread(program_buffer, sizeof(char), program_size, program_handle);
    if(ret_code != program_size)
    {
        printf("error reading file %s\n", program_file);
        if(feof(program_handle))
        {
            printf("Error reading file `%s`: unexpected end of file\n",program_file);
        }
        else if(ferror(program_handle))
        {
            perror("Error reading file `%s`");
        }
        exit(1);
    }
    fclose(program_handle);

    *program=clCreateProgramWithSource(*context, 1, (const char**) &program_buffer,
                                         &program_size, &err);

    free(program_buffer);

    printCLErr(err,__LINE__,__FILE__);
}

static cl_ulong get_event_timing(cl_event *event)
{
    cl_ulong p_start, p_end;
    clGetEventProfilingInfo(*event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong),
                            &p_start, NULL);
    clGetEventProfilingInfo(*event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong),
                            &p_end, NULL);
    return (p_end - p_start);
}

void GPU_Pack_A(inputDataType_x32 *A,
        unsigned int lda,
        DOUBLE *A_pack,
        unsigned int m,
        unsigned int k)
{
    DOUBLE *A_pack_local;
    unsigned int m_alg;
    //#pragma omp  parallel for num_threads(*n_threads) private(A_pack_local)
    //  #pragma omp parallel
    //  #pragma omp single
    ////    #pragma omp taskloop private(A_pack_local) num_tasks((*n_threads)) label(gemm_pack)
    //  #pragma omp taskloop private(A_pack_local) grainsize(16) label(gemm_pack)
    for(unsigned int ic=0;ic<m;ic+=GPU_BLOCK_MR)
    {
        A_pack_local=&A_pack[ic*k];
        m_alg=min(GPU_BLOCK_MR,m-ic);
        for(unsigned int pc=0;pc<k;pc++)
        {
            for(unsigned int ir=0;ir<m_alg;ir++)
            {
                A_pack_local[0]=A[(ic+ir)+pc*lda]; //auto xtypaei...
                A_pack_local++;
            }
        }
    }
}

void GPU_Pack_B(inputDataType_x32 *B,
        unsigned int ldb,
        DOUBLE *B_pack,
        unsigned int k,
        unsigned int n)
{
    DOUBLE *B_pack_local;
    unsigned int n_alg;
    //#pragma omp parallel for num_threads(*n_threads) private(B_pack_local)
    //  #pragma omp parallel
    //  #pragma omp single
    ////    #pragma omp taskloop private(B_pack_local) num_tasks((*n_threads)) label(gemm_pack)
    //  #pragma omp taskloop private(B_pack_local) grainsize(16) label(gemm_pack)
    for(unsigned int jc=0;jc<n;jc+=GPU_BLOCK_NR)
    {
        B_pack_local=&B_pack[jc*k];
        n_alg=min(GPU_BLOCK_NR,n-jc);
        for(unsigned int pc=0;pc<k;pc++)
        {
            for(unsigned int jr=0;jr<n_alg;jr++)
            {
                B_pack_local[0]=B[pc+jc*ldb+jr*ldb];
                B_pack_local++;
            }
        }
    }
}

void gpu_init(void)
{
    // ---- OpenCL stuff ---------------
    int err;

    // query number of platforms we have
    unsigned int num_platforms=0;
    err=clGetPlatformIDs(1, NULL, &num_platforms);
    printCLErr(err,__LINE__,__FILE__);

    // allocate array, an entry for each platform found
    platforms=(cl_platform_id *) malloc(sizeof(cl_platform_id) * num_platforms);
    assert(platforms);
    // place the `cl_platform_id` structures in the platforms array
    err=clGetPlatformIDs(num_platforms, platforms, NULL);
    printCLErr(err,__LINE__,__FILE__);

    // determine number of devices
    // NOTE: arbitrarily pick first platform
    unsigned int num_devices=0;
    err=clGetDeviceIDs(platforms[0], CL_DEVICE_TYPE_GPU, 1, NULL, &num_devices);
    printCLErr(err,__LINE__,__FILE__);
    // allocate memory for device array
    devices=(cl_device_id*) malloc(sizeof(cl_device_id) * num_devices);
    assert(devices);
    // populate device array
    // NOTE: arbitrarily pick first platform
    err=clGetDeviceIDs(platforms[0], CL_DEVICE_TYPE_GPU, num_devices, devices, NULL);
    printCLErr(err,__LINE__,__FILE__);

    context=clCreateContext(NULL, 1, &devices[0], NULL, NULL, &err);
    printCLErr(err,__LINE__,__FILE__);

    create_program_with_source(&program, &context, PROGRAM_FILE);

    // add any kernel compiler options to this string
    const char* options="-cl-mad-enable";
    err=clBuildProgram(program, 1, &devices[0], options, NULL, NULL);
    // print build errors
    if(err != CL_SUCCESS)
    {
        perror("error during build");
        size_t log_size=0;
        clGetProgramBuildInfo(program, devices[0], CL_PROGRAM_BUILD_LOG, 0, NULL,
                              &log_size);
        char *program_log=(char*)malloc(log_size+1);
        assert(program_log);
        program_log[log_size]='\0';
        clGetProgramBuildInfo(program, devices[0], CL_PROGRAM_BUILD_LOG,
                              log_size+1, program_log, NULL);
        printf("%s\n", program_log);
        free(program_log);
        exit(1);
    }

    io_queue=clCreateCommandQueue(context, devices[0], CL_QUEUE_PROFILING_ENABLE, &err);
    printCLErr(err,__LINE__,__FILE__);

    compute_queue=clCreateCommandQueue(context, devices[0], CL_QUEUE_PROFILING_ENABLE,
                                       &err);
    printCLErr(err,__LINE__,__FILE__);

    // NOTE: assume device has enough memory
    // TODO: does performance degrade if k is much less than KC?
    cl_ulong a_buffer_size=GPU_BLOCK_MC * GPU_BLOCK_KC * sizeof(inputDataType_x32);
    cl_ulong b_buffer_size=GPU_BLOCK_KC * GPU_BLOCK_NC * sizeof(inputDataType_x32);
    cl_ulong c_buffer_size=GPU_BLOCK_MC * GPU_BLOCK_NC * sizeof(inputDataType_x32);

    cl_ulong max_alloc=0;
    err=clGetDeviceInfo(devices[0], CL_DEVICE_MAX_MEM_ALLOC_SIZE,sizeof(max_alloc),
                        &max_alloc, NULL);
    printCLErr(err,__LINE__,__FILE__);

    cl_ulong global_mem=0;
    err=clGetDeviceInfo(devices[0], CL_DEVICE_GLOBAL_MEM_SIZE,sizeof(global_mem),
                        &global_mem, NULL);
    printCLErr(err,__LINE__,__FILE__);

    // printf("\n%lu %lu %lu\n", a_buffer_size, b_buffer_size, c_buffer_size);

    cl_ulong total=2*a_buffer_size+2*b_buffer_size+2*c_buffer_size;
    if(total > global_mem)
    {
        printf("not enough global storage!\n");
        exit(1);
    }
    if((a_buffer_size > max_alloc) ||
       (b_buffer_size > max_alloc) ||
       (c_buffer_size > max_alloc))
    {
        printf("some buffer is too big!\n");
        exit(1);
    }

    unsigned int i;
    // NOTE: this is what's being passed in for the other GEMM implementation
    rs_c=1;
    cs_c=GPU_BLOCK_MC;
    for(i=0; i < 2; i++)
    {
        // create double buffers for input A matrix
        // printf("a buffer %u\n", i);
        a_buffers[i]=clCreateBuffer(
                context, CL_MEM_READ_ONLY,
                a_buffer_size, NULL, &err
                );
        printCLErr(err,__LINE__,__FILE__);

        // create double buffers for input B matrix
        // printf("b buffer %u\n", i);
        b_buffers[i]=clCreateBuffer(
                context, CL_MEM_READ_ONLY,
                b_buffer_size, NULL, &err
                );
        printCLErr(err,__LINE__,__FILE__);

        // TODO: eventually change this to not just write only??
        // create output buffers for C matrix
        // printf("c buffer %u\n", i);
        c_buffers[i]=clCreateBuffer(
                context, CL_MEM_READ_WRITE,
                c_buffer_size, NULL, &err
                );
        printCLErr(err,__LINE__,__FILE__);

        /* TODO:
           What's better: using an intermediate buffer between C and c_buffer, or
           many calls to C clEnqueueReadBuffer/clEnqueueWriteBuffer (using offset in
           a for loop).
           when the writes are blocking, the intermediate looks better.  try using
           event queueing...
           */
        // c_sub_matrix[i]=calloc(c_buffer_size,1);

        // create kernels
        kernels[i]=clCreateKernel(program, kernel_name, &err);
        printCLErr(err,__LINE__,__FILE__);

        // TODO: might have to move this inside the loops depending on
        // what parameters we're using (ir, mr, jr, nr, pr ?) can this be
        // substituted by global ids if A is in col major, B is in row major?
        err |= clSetKernelArg(kernels[i], 3, sizeof(cl_mem), &a_buffers[i]);
        // NOTE: don't do B here since that alternates, note later.
        err |= clSetKernelArg(kernels[i], 5, sizeof(cl_mem), &c_buffers[i]);
        err |= clSetKernelArg(kernels[i], 6, sizeof(unsigned int), &rs_c);
        err |= clSetKernelArg(kernels[i], 7, sizeof(unsigned int), &cs_c);
        printCLErr(err,__LINE__,__FILE__);
    }
    // ---- end OpenCL stuff -----------
}

void gpu_release(void)
{
    for(int i=0; i < 2; i++)
    {
        clReleaseKernel(kernels[i]);
        clReleaseMemObject(c_buffers[i]);
        clReleaseMemObject(b_buffers[i]);
        clReleaseMemObject(a_buffers[i]);
        // free(c_sub_matrix[i]);
        //        for(int j=0; j < 4; j++)
        //            clReleaseEvent(events[(2*i)+j]);
    }
    // clReleaseCommandQueue(compute_queue);
    clReleaseCommandQueue(io_queue);
    clReleaseCommandQueue(compute_queue);
    clReleaseProgram(program);
    clReleaseContext(context);
    free(devices);
    free(platforms);
}

void gpu_gemm(unsigned int m,
        unsigned int n,
        unsigned int k,
        inputDataType_x32 *A,
        unsigned int lda,
        inputDataType_x32 *B,
        unsigned int ldb,
        inputDataType_x32 *C,
        unsigned int ldc,
        void *Ac_pack_v,
        void *Bc_pack_v)
{
    inputDataType_x32 *Ac, *Bc;
    inputDataType_x32 *Cc;
    //unsigned int *Ar, *Br;
    //unsigned int *Cr;
    //DOUBLE beta;
    unsigned int i=0;
    int err;
    const unsigned int work_dim=2;
    size_t local[work_dim];
    local[0]=LOCAL_0;
    local[1]=LOCAL_1;
    size_t global[work_dim];

    cl_ulong p_total=0;

    for(unsigned int jc=0; jc<n; jc+=GPU_BLOCK_NC)
    {
        unsigned int n_alg=min(GPU_BLOCK_NC,n-jc);

        for(i=0; i < 2; i++)
        {
            err=clSetKernelArg(kernels[i], 2, sizeof(inputDataType_x32), &n_alg);
            printCLErr(err,__LINE__,__FILE__);
        }

        unsigned int pc_iter=0;
        for(unsigned int pc=0; pc<k; pc+=GPU_BLOCK_KC)
        {
            unsigned int k_alg=min(GPU_BLOCK_KC,k-pc);

            for(i=0; i < 2; i++)
            {
                err=clSetKernelArg(kernels[i], 0, sizeof(inputDataType_x32), &k_alg);
                printCLErr(err,__LINE__,__FILE__);
            }

            //beta=*betap;

            Bc=&B[pc+jc*ldb];
            GPU_Pack_B(Bc, ldb, Bc_pack_v, k_alg, n_alg);

            // TODO: double buffer B, don't just use b_buffers[0]
            // printf("\nwriting b (k: %u, n: %u)\n", k_alg, n_alg);

            // NOTE: getting CL_MEM_OBJECT_ALLOCATION_FAILURE when
            // KC*NC*sizeof(uint) is too big... but not more than max_alloc??

            // NOTE: don't have to use an event here since
            //        io_queue is processed in-order
            err=clEnqueueWriteBuffer(
                    io_queue, b_buffers[pc_iter % 2], CL_FALSE, 0,
                    k_alg*n_alg*sizeof(inputDataType_x32), Bc_pack_v,
                    0, NULL, NULL
                    );
            printCLErr(err,__LINE__,__FILE__);

            // printf("write done\n");
            // set both kernels to use this iteration's b_buffer
            for(i=0; i < 2; i++)
            {
                err=clSetKernelArg(
                        kernels[i], 4, sizeof(cl_mem), &b_buffers[pc_iter % 2]
                        );
                printCLErr(err,__LINE__,__FILE__);
            }

            unsigned int ic_iter=0;
            // write A, write C, kernel, read C (*2 for double buffer)

            for(unsigned int ic=0; ic<m; ic+=GPU_BLOCK_MC)
            {

                unsigned int m_alg=min(GPU_BLOCK_MC,m-ic);

                // m_alg is MC until it is the tail
                // global is going to be derived from m_alg/n_alg.
                // at - least as big as NR/MR (1 block).
                // ceiling division
                // TODO: figure this out for larger MR/NR
                // global[0]=((max(m_alg, LOCAL_0) + LOCAL_0 - 1) / LOCAL_0) * LOCAL_0 / GPU_BLOCK_MR;
                global[0]=round_up_mult(
                        round_up_mult(m_alg, LOCAL_0), BLOCK_SIZE_X
                        ) / BLOCK_SIZE_X * LOCAL_0;
                // (((n_alg + LOCAL_0 - 1)/LOCAL_0) * LOCAL_0) * LOCAL_0 / BLOCK_SIZE_X;
                global[1]=round_up_mult(
                        round_up_mult(n_alg, LOCAL_1), BLOCK_SIZE_Y
                        ) / BLOCK_SIZE_Y * LOCAL_1;

                global[0]=max(global[0], local[0]);
                global[1]=max(global[1], local[1]);

                // global[0]=round_up_mult(
                //   round_up_mult(n_alg, LOCAL_0) * LOCAL_0, BLOCK_SIZE_X
                // ) / BLOCK_SIZE_X;
                // // (((n_alg + LOCAL_0 - 1)/LOCAL_0) * LOCAL_0) * LOCAL_0 / BLOCK_SIZE_X;
                // global[1]=round_up_mult(
                //   round_up_mult(m_alg, LOCAL_1) * LOCAL_1, BLOCK_SIZE_Y
                // ) / BLOCK_SIZE_Y;



                //(((m_alg + LOCAL_1 - 1)/LOCAL_1) * LOCAL_1) * LOCAL_1 / BLOCK_SIZE_Y;
                // global[1]=((max(n_alg, LOCAL_1) + LOCAL_1 - 1) / LOCAL_1) * LOCAL_1 / GPU_BLOCK_NR;

                err=clSetKernelArg(
                        kernels[ic_iter % 2], 1, sizeof(inputDataType_x32), &m_alg
                        );
                printCLErr(err,__LINE__,__FILE__);

                Ac=&A[ic+pc*lda];
                GPU_Pack_A(Ac,lda,Ac_pack_v,m_alg,k_alg);

                Cc=&C[ic+jc*ldc];

                for(i=0; i < n_alg; ++i)
                {
                    // for(i=0; i < GPU_BLOCK_NC; ++i) {
                    err=clEnqueueWriteBuffer(
                            io_queue, c_buffers[ic_iter % 2], CL_FALSE,
                            i*cs_c*sizeof(inputDataType_x32),
                            m_alg*sizeof(inputDataType_x32), &Cc[i*ldc],
                            //GPU_BLOCK_MC*sizeof(unsigned int), &Cc[i*ldc],
                            0, NULL, NULL
                            );
                    printCLErr(err,__LINE__,__FILE__);
                }

                // printf("writing a\n");

                err=clEnqueueWriteBuffer(
                        io_queue, a_buffers[ic_iter % 2], CL_FALSE, 0,
                        m_alg*k_alg*sizeof(inputDataType_x32), Ac_pack_v,
                        0, NULL, &events[(ic_iter % 2)*4]
                        );
                printCLErr(err,__LINE__,__FILE__);

                // assuming C starts cleared, but we want to += if iterating over k (pc)
                // for(i=0; i < n_alg; ++i) {
                //   memcpy(
                //     &(c_sub_matrix[ic_iter % 2][i*cs_c]), &Cc[i*ldc], m_alg*sizeof(unsigned int)
                //   );
                // }
                //
                // err=clEnqueueWriteBuffer(
                //   io_queue, c_buffers[ic_iter % 2], CL_FALSE, 0,
                //   GPU_BLOCK_MC*GPU_BLOCK_NC*sizeof(unsigned int),
                //                                             c_sub_matrix[ic_iter % 2],
                //   0, NULL, &events[(ic_iter % 2)*4 + 1]
                // );
                // printCLErr(err,__LINE__,__FILE__);


                // printf(
                //   "\nm: %u, n: %u, k: %u, ic: %u, jc: %u, loc: %lu,%lu, glob: %lu,%lu\n",
                //   m_alg, n_alg, k_alg, ic, jc, local[0], local[1], global[0], global[1]
                // );

                err=clEnqueueNDRangeKernel(
                        compute_queue, kernels[ic_iter % 2], work_dim, NULL, global, local,
                        1, &events[(ic_iter % 2)*4], &events[(ic_iter % 2)*4 + 2]
                        );
                printCLErr(err,__LINE__,__FILE__);

                // err=clEnqueueReadBuffer(
                //   io_queue, c_buffers[ic_iter % 2], CL_TRUE, 0,
                //   GPU_BLOCK_MC*GPU_BLOCK_NC*sizeof(unsigned int),
                //                                               c_sub_matrix[ic_iter % 2],
                //   1, &events[(ic_iter % 2)*4 + 2], &events[(ic_iter % 2)*4 + 3]
                // );
                // printCLErr(err,__LINE__,__FILE__);

                // printf("\n");
                // print_matrix(c_sub_matrix[0], GPU_BLOCK_MC, GPU_BLOCK_NC);
                // printf("\n");

                // add event timing (since kernel is done -- TODO: make this callback too)

                for(i=0; i < n_alg; ++i)
                {
                    // for(i=0; i < GPU_BLOCK_NC; ++i) {
                    err=clEnqueueReadBuffer(
                            io_queue, c_buffers[ic_iter % 2], CL_FALSE,
                            i*cs_c*sizeof(inputDataType_x32),
                            m_alg*sizeof(inputDataType_x32), &Cc[i*ldc],
                            // GPU_BLOCK_MC*sizeof(unsigned int), &Cc[i*ldc],
                            1, &events[(ic_iter % 2)*4 + 2], NULL
                            );
                    printCLErr(err,__LINE__,__FILE__);
                }

                //                getc(stdin); //used to pause program for monitoring - MPAMPIS -

                clFinish(io_queue);
                clFinish(compute_queue);
                p_total += get_event_timing(&events[(ic_iter % 2)*4 + 2]);

                // TODO: do this in a callback function and make read C asynchronous
                // n_alg rows, each row ldc or cs_c
                // for(i=0; i < n_alg; ++i) {
                //   memcpy(
                //     &Cc[i*ldc], &(c_sub_matrix[ic_iter % 2][i*cs_c]),
                //     m_alg*sizeof(unsigned int)
                //   );
                // }

                ic_iter += 1;
            }
            pc_iter += 1;
        }
    }

    // more timing stuff
#ifdef VERBOSE
    cl_double kernelExecTimeMs=(cl_double)p_total*(cl_double)(1e-09);
    printf("\nopencl kernel (%s): %lf seconds\n", kernel_name, kernelExecTimeMs);
#endif
}

void get_pairwise_ld_score_gpu(unsigned int * tableA_bitcount,
        unsigned int * tableB_bitcount,
        inputDataType_x32 * C,
        int tableAsize,
        int tableBsize,
        int snp_size,
        ResultDataType** results)
{
    int i,j;
    ResultDataType val_1, val_2, val_3;
    for(i=0;i<tableBsize;i++)
    {
        for(j=0;j<tableAsize;j++)
        {
            (*results)[i*tableAsize+j]=0.0f;
            if(tableB_bitcount[i] != 0 && tableA_bitcount[j] != 0)
            {
                val_1=((ResultDataType)tableA_bitcount[j])/snp_size;
                val_2=((ResultDataType)tableB_bitcount[i])/snp_size;
                val_3=((ResultDataType)C[i*tableAsize+j])/snp_size;
                (*results)[i*tableAsize+j]=((val_3-val_1*val_2)*(val_3-val_1*val_2));
                (*results)[i*tableAsize+j] /= (val_1*val_2*(1.0-val_1)*(1.0-val_2));

                if((*results)[i*tableAsize+j]>1.0001)
                {
                    fprintf(stderr, "qLD-compute: gemm_gpu.c:722: \
get_pairwise_ld_score: Entry i %d j %d greater than 1.\n\
%u %u %u, result = %f\n", i, j, tableA_bitcount[j], tableB_bitcount[i],
                            C[i*tableAsize+j], (*results)[i*tableAsize+j]);
                }
            }
        }
    }
    fflush(stderr);
}

void mlt_gpu(unsigned int m,
        unsigned int k,
        inputDataType_x32 *A,
        inputDataType_x32 *tableA)
{
    for(unsigned int i=0;i<m;i++)
    {
        for(unsigned int j=0;j<k;j++)
        {
            ((inputDataType_x32*)A)[j*m + i]=tableA[i*k + j];
        }
    }
}

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
        int posWset2)
{
    int m=tableAsize, n=tableBsize, k=compressed_snp_size,i;
    int pm;     //return value used for assert
    double GFLOPS_BLIS=0.0;
    struct timeval start,end;
    double time_ref=0.0;
    double ops=2.0*m*n*k;

    double time000, time001, time002, time003, time004, time005;
    double init_time=0.0, mlt_time=0.0, gemm_time=0.0, ld_time=0.0, write_time=0.0;

    time000=gettime(); //000 starting

    void *Ac_pack_v=NULL, *Bc_pack_v=NULL, *A=NULL, *C=NULL;
    pm=posix_memalign(&(Ac_pack_v), 4096,
            12*GPU_BLOCK_MC*GPU_BLOCK_KC*sizeof(inputDataType_x32));
    assert(!pm);
    pm=posix_memalign(&(Bc_pack_v), 4096,
            GPU_BLOCK_KC*GPU_BLOCK_NC*sizeof(inputDataType_x32));
    assert(!pm);

    long long int tableCsize=m*n;

    pm=posix_memalign(&A, 4096, m*k*sizeof(inputDataType_x32) +
            m*k*sizeof(inputDataType_x32)%4096);
    assert(!pm);
    pm=posix_memalign(&C, 4096, tableCsize*sizeof(inputDataType_x32) +
            tableCsize*sizeof(inputDataType_x32)%4096);
    assert(!pm);

    ResultDataType *results =(ResultDataType*)malloc(tableCsize*sizeof(ResultDataType));
    assert(results);
    gettimeofday( &start, NULL );

    writeResultsHeader(fpOut);

    init_time=gettime() - time000;

    time001=gettime(); //001 before mlt

    for(i=0;i<tableCsize;i++)
    {
        ((inputDataType_x32*)C)[i]=0;
        results[i]=0;
    }

    mlt_gpu(m, k, A, tableA);

    time002=gettime(); //002 after mlt
    mlt_time += time002-time001;

    //int jc=1, ic=1, jr=1, ir=1;
    //    DOUBLE alphap=1.0;
    //    DOUBLE betap=0.0;

    gpu_gemm(m,
            n,
            k,
            A,
            m,
            tableB,
            k,
            C,
            m,
            Ac_pack_v,
            Bc_pack_v);

    time003=gettime(); //004 after gemm
    gemm_time += time003 - time002;

    get_pairwise_ld_score_gpu(tableA_bitcount,
            tableB_bitcount,
            C,
            m,
            n,
            snp_size,
            &results);

    time004=gettime(); //004 after ldScore
    ld_time += time004 - time003;

    writeResults(fpOut,
            tableA_posIndex,
            tableA_IDindex,
            tableA_bitcount,
            m,
            tableB_posIndex,
            tableB_IDindex,
            tableB_bitcount,
            n,
            snp_size,
            results,
            threadData[0].r2limit,
            posWset2);

    time005=gettime(); //005 after writeRes
    write_time += time005 - time004;
    gettimeofday( &end, NULL );

#if defined(VERBOSE) || defined(BENCHMARK)
    threadData[0].threadStats.init_time+=init_time;
    threadData[0].threadStats.mlt_time+=mlt_time;
    threadData[0].threadStats.gemm_time+=gemm_time;
    threadData[0].threadStats.ld_time+=ld_time;
    threadData[0].threadStats.write_time+=write_time;

    fprintf(threadData[0].threadLog,"[T] Init times: %f\n",init_time);
    fprintf(threadData[0].threadLog,"[T] MLT times: %f\n", mlt_time);
    fprintf(threadData[0].threadLog,"[T] GEMM times: %f\n", gemm_time);
    fprintf(threadData[0].threadLog,"[T] LD times: %f\n", ld_time);
    fprintf(threadData[0].threadLog,"[T] Write times: %f\n\n", write_time);

    time_ref=(double)((end.tv_sec-start.tv_sec) * 1000000.0 +
            (end.tv_usec-start.tv_usec))/1000000.0;
    GFLOPS_BLIS=ops/(time_ref*1.0e9);
    fprintf(threadData[0].threadLog,"Table A size: %d, Table B size: %d,"
"Compressed snip size: %d, Time: %5.6fs, GFLOPS: %3.3f\n",m, n, k, time_ref, GFLOPS_BLIS);
#else
    threadData[0].threadID=threadData[0].threadID;
    init_time=init_time;
    mlt_time=mlt_time;
    gemm_time=gemm_time;
    ld_time=ld_time;
    write_time=write_time;
    ops=ops;
    time_ref=time_ref;
    GFLOPS_BLIS=GFLOPS_BLIS;
#endif
    if(Ac_pack_v)
        free(Ac_pack_v);
    if(Bc_pack_v)
        free(Bc_pack_v);
    if(A)
        free(A);
    if(C)
        free(C);
    if(results)
        free(results);
}
#endif
