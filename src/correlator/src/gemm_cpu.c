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
#include "../include/header.h"
#include "../include/gemm_cpu.h"
#include "../include/correlate_IO.h"
#include "../include/correlate_IOx64.h"
#include "../include/snip_processing.h"
#include "../include/read_file.h"

void Pack_A(inputDataType_x64 *A,
            unsigned int lda,
            inputDataType_x64 *A_pack,
            unsigned int m,
            unsigned int k)
{
	inputDataType_x64 *A_pack_local;
	for(unsigned int ic=0;ic<m;ic+=BLOCK_MR)
    {
		A_pack_local=&A_pack[ic*k];
		unsigned int m_alg=fmin(BLOCK_MR,m-ic);
		for(unsigned int pc=0;pc<k;pc++)
        {
			for(unsigned int ir=0;ir<m_alg;ir++)
            {
				A_pack_local[0]=A[(ic+ir)+pc*lda];
				A_pack_local++;
			}
		}
	}
}

void Pack_B(inputDataType_x64 *B,
            unsigned int ldb,
            inputDataType_x64 *B_pack,
            unsigned int k,
            unsigned int n)
{
    inputDataType_x64 *B_pack_local;
	for(unsigned int jc=0;jc<n;jc+=BLOCK_NR)
    {
	    B_pack_local=&B_pack[jc*k];
        unsigned int n_alg=fmin(BLOCK_NR,n-jc);
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

void dgemm_ref(int k,
               int mr_alg,
               int nr_alg,
               inputDataType_x64 alpha,
               inputDataType_x64* restrict a,
               inputDataType_x64* restrict b,
               inputDataType_x64* restrict c,
               int rs_c,
               int cs_c)
{
	inputDataType_x64 ab[mr_alg*nr_alg];
	inputDataType_x64 bj, ai;
	unsigned int l,i,j;
	for(i=0;i < (unsigned int)(mr_alg * nr_alg); ++i) //set 0s
        *(ab+i)=0;

	for(l=0;l < (unsigned int)k; ++l)
    {
		inputDataType_x64 *abij=ab;
        for(j=0;j < (unsigned int)nr_alg; ++j)
        {
            bj = *(b + j);
            for(i=0;i < (unsigned int)mr_alg; ++i)
            {
                ai = *(a + i);
                if(sizeof(inputDataType_x64) <= 4)
                    abij[0] += __builtin_popcount(ai & bj);
                else if(sizeof(inputDataType_x64) < 8)
                    abij[0] += __builtin_popcountl(ai & bj);
                else if(sizeof(inputDataType_x64) >= 8)
                    abij[0] += __builtin_popcountll(ai & bj);
                abij++;
            }
        }
        a+=mr_alg;
        b+=nr_alg;
    }

	for(i=0;i < (unsigned int)(mr_alg * nr_alg); ++i) //Scale by alpha
        ab[i]*=alpha;

	inputDataType_x64 *Cr_l;
	inputDataType_x64 *ab_l=ab;

	for(j=0;j < (unsigned int)nr_alg; ++j)
    {
        Cr_l=&c[j*cs_c];
        for(i=0;i < (unsigned int)mr_alg; ++i)
        {
            *Cr_l += *ab_l;
            ab_l++;
            Cr_l += rs_c;
        }
    }
}
#ifdef CBLAS_USE
void blis_gemm(unsigned int M,
               unsigned int N,
               unsigned int K,
               inputDataType_x64 alphap,
               inputDataType_x64 *a,
               inputDataType_x64 *b,
               inputDataType_x64 *c,
               int posWset)
{
    //inputDataType_x64 *Ac, *Bc;
	//inputDataType_x64 *Cc;
	//inputDataType_x64 *Ar, *Br;
	//inputDataType_x64 *Cr;

    //printf("table A:\n");
    //for(unsigned int i=0; i < M; i++)
    //{
    //    for(unsigned int tk=0; tk < K; tk++)
    //    {
    //        printf("%20lu ",a[tk+K*i]);
    //    }
    //    printf("\n");
    //}
    //printf("\n");
    //printf("table B:\n");
    //for(unsigned int i=0; i < N; i++)
    //{
    //    for(unsigned int tk=0; tk < K; tk++)
    //    {
    //        printf("%20lu ",b[tk+K*i]);
    //    }
    //    printf("\n");
    //}
    //printf("\n");

    //printf("table C before:\n");
    //for(unsigned int i=0; i < M; i++)
    //{
    //    for(unsigned int tk=0; tk < N; tk++)
    //    {
    //        printf("%4lu ",c[tk+N*i]);
    //    }
    //    printf("\n");
    //}
    //printf("\n");

    //for(unsigned int i=0; i < M; i++)
    //{
    //    for(unsigned int tk=0; tk < K; tk++)
    //    {
    //        if(a[tk+K*i] != b[tk+K*i])
    //        {
    //            printf("Mismatch in [%d,%d]\n", i, tk);
    //        }
 
    //        if(a[tk+K*i] != (double)a[tk+K*i])
    //        {
    //            printf("Uint != Double in A[%d,%d]\n", i, tk);
    //        }
    //        if(b[tk+K*i] != (double)b[tk+K*i])
    //        {
    //            printf("Uint != Double in B[%d,%d]\n", i, tk);
    //        }
    //    }
    //}

    //FILE* diag;
    //FILE* second_diag;
    //FILE* other_elements;

#ifdef SYRK
    if(posWset)
    {
#endif
        //diag=fopen("diagonal_blis.txt","w");
        //second_diag=fopen("second_diagonal_blis.txt","w");
        //other_elements=fopen("other_elements_blis.txt","w");
        
        cblas_dgemm(CblasColMajor,
                    CblasTrans,
                    CblasNoTrans,
                    M,
                    N,
                    K,
                    (double)alphap,
                    (double*)a,
                    K,
                    (double*)b,
                    K,
                    0.0,
                    (double*)c,
                    M);
#ifdef SYRK
    }
    else
    {
        //diag=fopen("diagonal_syrk.txt","w");
        //second_diag=fopen("second_diagonal_syrk.txt","w");
        //other_elements=fopen("other_elements_syrk.txt","w");

        printf("SYRK Enabled\n");
        cblas_dsyrk(CblasColMajor,
                    CblasLower,
                    CblasTrans,
                    M,
                    K,
                    (double)alphap,
                    (double*)a,
                    K,
                    0.0,
                    (double*)c,
                    M);
    }
#endif
    
    //printf("table C as double:\n");
    //for(unsigned int i=0; i < M; i++)
    //{
    //    for(unsigned int tk=0; tk < N; tk++)
    //    {
    //        printf("%6.2f ",(double)c[tk+N*i]);
    //    }
    //    printf("\n");
    //}
    //printf("\n");

    //printf("table C as long unsigned:\n");
    //for(unsigned int i=0; i < M; i++)
    //{
    //    for(unsigned int tk=0; tk < N; tk++)
    //    {
    //        printf("%4lu ",(long unsigned int)c[tk+N*i]);

    //        if(tk == i)
    //        {
    //            fprintf(diag, "%lu\n", c[tk+N*i]);
    //        }
    //        else if(tk == i+1)
    //        {
    //            fprintf(second_diag, "%lu\n",c[tk+N*i]);
    //        }
    //        else if(tk > i+1)
    //            fprintf(other_elements, "%lu\n",c[tk+N*i]);

    //    }
    //    printf("\n");
    //}
    //printf("\n");

    //fclose(diag);
    //fclose(second_diag);
    //fclose(other_elements);
}
#endif

void gemm(unsigned int m,
          unsigned int n,
          unsigned int k,
		  inputDataType_x64 alphap,
          inputDataType_x64 *A,
          unsigned int lda,
          inputDataType_x64 *B,
          unsigned int ldb,
          inputDataType_x64 *C,
          unsigned int ldc,
          void * Ac_pack_v,
          void * Bc_pack_v)
{
	inputDataType_x64 *Ac, *Bc;
	inputDataType_x64 *Cc;
    inputDataType_x64 *Ar, *Br;
	inputDataType_x64 *Cr;
    inputDataType_x64 *Ac_pack=(inputDataType_x64 *)Ac_pack_v;
	inputDataType_x64 *Bc_pack=(inputDataType_x64 *)Bc_pack_v;

    for (unsigned int jc=0; jc<n; jc+=BLOCK_NC)
	{
		unsigned int n_alg=fmin(BLOCK_NC,n-jc);
		for (unsigned int pc=0; pc<k; pc+=BLOCK_KC)
		{
			unsigned int k_alg=fmin(BLOCK_KC,k-pc);

			Bc=&B[pc+jc*ldb];
            Pack_B(Bc, ldb, Bc_pack, k_alg, n_alg);  //PACK B
			for (unsigned int ic=0; ic<m; ic+=BLOCK_MC)
			{
				unsigned int m_alg=fmin(BLOCK_MC,m-ic);
				Ac=&A[ic+pc*lda];
				inputDataType_x64 *Ac_pack_local=Ac_pack; // Ac pack pointer per Loop 3 thread
                Pack_A(Ac,lda,(inputDataType_x64*)Ac_pack_local,m_alg,k_alg); //PACK A
				Cc=&C[ic+jc*ldc];
                for(unsigned jr=0;jr<n_alg;jr+=BLOCK_NR)
				{
					unsigned int nr_alg=fmin(BLOCK_NR,n_alg-jr);
					for(unsigned int ir=0;ir<m_alg;ir+=BLOCK_MR)
					{
						unsigned int mr_alg=fmin(BLOCK_MR,m_alg-ir);
						Ar=&Ac_pack_local[ir*k_alg];
						Br=&Bc_pack[jr*k_alg];
						Cr=&Cc[ir+jr*ldc];

                        dgemm_ref(k_alg,mr_alg,nr_alg,alphap,Ar,Br,Cr,1,ldc);
					}
				}
			}
		}
	}
#ifdef CBLAS_USE
    Ar=Ar;
    Br=Br;
	Cr=Cr;
    Ac_pack=Ac_pack;
    Ac_pack_v=Ac_pack_v;
	Bc_pack=Bc_pack;
    Bc_pack_v=Bc_pack_v;
#endif
}

void get_pairwise_ld_score(unsigned int * tableA_bitcount,
                           unsigned int * tableB_bitcount,
                           inputDataType_x64 * C,
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
			(*results)[i*tableAsize+j] = 0.0f;
	        if(tableB_bitcount[i] != 0 && tableA_bitcount[j] != 0)
	        {
			    val_1=((ResultDataType)tableA_bitcount[j])/snp_size;
			    val_2=((ResultDataType)tableB_bitcount[i])/snp_size;
			    val_3=((ResultDataType)C[i*tableAsize+j])/snp_size;
			    (*results)[i*tableAsize+j] = ((val_3-val_1*val_2)*(val_3-val_1*val_2));
			    (*results)[i*tableAsize+j] /= (val_1*val_2*(1.0-val_1)*(1.0-val_2));
                
                //printf("\tA[%d]\tB[%d]\tC[%u]\n",
                //        j,i,i*tableAsize+j);
                //printf("%8u%8u%8lu\n\n",
                //        tableA_bitcount[j],
                //        tableB_bitcount[i],
                //        C[i*tableAsize+j]);
                
                //printf("A[%d]: %8u\t%f\n",j,tableA_bitcount[j],val_1);
                //printf("B[%d]: %8u\t%f\n",i,tableB_bitcount[i],val_2);
                //printf("C[%d]: %8lu\t%f\n\n",i*tableAsize+j,C[i*tableAsize+j],val_3);


                if((*results)[i*tableAsize+j]>1.0001)
		        {
//			        fprintf(stderr, "qLD-compute: gemm_cpu.c:247: "
//			                        "get_pairwise_ld_score: "
//			                        "Entry i %d j %d greater than 1.\n"
//                                  "%u %u %lu, result = %f\n", 
//                                  i, j, tableA_bitcount[j], tableB_bitcount[i],
//                                  C[i*tableAsize+j], (*results)[i*tableAsize+j]);
                    
                    //printf("Danger! High Voltage in %d:%d! %f %f %f\n",i,j,val_3, val_1, val_2);
                    (*results)[i*tableAsize+j]=123.456000000;
		        }
	        }
		}
	}

    //printf("Results:\n");
    //for(int i=0; i < tableAsize; i++)
    //{
    //    for(int tk=0; tk < tableBsize; tk++)
    //    {
    //        printf("%6.2f ",(*results)[tk+tableBsize*i]);
    //    }
    //    printf("\n");
    //}
    //printf("\n");
}

void mlt(unsigned int m,
         unsigned int k,
         inputDataType_x64* A,
         inputDataType_x64* tableA)
{
	for(unsigned int i=0;i<m;i++)
	{
		for(unsigned int j=0;j<k;j++)
		{
			((inputDataType_x64*)A)[j*m + i] = tableA[i*k + j];
		}
	}
}

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
               int posWset2)
{
	int m=tableAsize, n=tableBsize, k=compressed_snp_size;
    long long int i;
    long long int tableCsize = (long long int)m*(long long int)n;
    int pm;     //return value used for assert
	double GFLOPS_BLIS=0.0;
	struct timeval start,end;
	double time_ref=0.0;
	double ops=2.0*m*n*k;

    double time000, time001, time002, time003, time004, time005;
    double init_time=0.0, mlt_time=0.0, gemm_time=0.0, ld_time=0.0, write_time=0.0;

    time000 = gettime(); //000 starting

	void *Ac_pack_v=NULL, *Bc_pack_v=NULL, *A=NULL, *C=NULL;

    pm = posix_memalign(&(Ac_pack_v), 4096,
                        12*BLOCK_MC*BLOCK_KC*sizeof(inputDataType_x64));
    assert(!pm);
    pm = posix_memalign(&(Bc_pack_v), 4096,
                        BLOCK_KC*BLOCK_NC*sizeof(inputDataType_x64));
    assert(!pm);

    pm = posix_memalign(&A, 4096, m*k*sizeof(inputDataType_x64) +
                         m*k*sizeof(inputDataType_x64)%4096);
    assert(!pm);

    pm = posix_memalign(&C, 4096, tableCsize*sizeof(inputDataType_x64) +
                        tableCsize*sizeof(inputDataType_x64)%4096);
    assert(!pm);

	ResultDataType *results =(ResultDataType*)malloc(tableCsize*sizeof(ResultDataType));
    assert(results);

	inputDataType_x64  alphap=1.0;
	gettimeofday( &start, NULL );

    writeResultsHeader(fpOut);

    init_time = gettime() - time000;

    time001 = gettime(); //001 before mlt

	for(i=0;i<tableCsize;i++)
	{
		((inputDataType_x64*)C)[i] = 0;
		results[i]=0;
	}

    mlt(m, k, A, tableA);
    
    time002 = gettime(); //002 after mlt
    mlt_time += time002-time001;
    
#ifdef COUNTSTATES
    uint8_t *STtable=(uint8_t*)malloc(m*n*sizeof(uint8_t));
    assert(STtable);
    for(uint32_t indA=0; indA < (uint32_t)m; indA++)
    {
        for(uint32_t indB=0; indB < (uint32_t)n; indB++)
        {
            STtable[indA*n+indB]=0;
            for(uint32_t compD=0; compD < (uint32_t)k; compD++)
            {
                for(int bitP=0; bitP < 64; bitP++)
                {
                    if(STtable[indA*n+indB] == 15)
                        break;
                    // STtable OR 1 in position "bitA*2 + bitB" where
                    // states are: A=0,B=0 -> STtable[0]
                    //             A=0,B=1 -> STtable[1]
                    //             A=1,B=0 -> STtable[2]
                    //             A=1,B=1 -> STtable[3]
                    STtable[indA*n+indB]|=(1 << (((tableA[indA*k+compD] >> bitP) & 1)*2
                                                +((tableB[indB*k+compD] >> bitP) & 1)));
                }
                STtable[indA*n+indB]=__builtin_popcount(STtable[indA*n+indB]);
            }
        }
    }
#endif
    if(threadData[0].blis)
    {
#ifdef CBLAS_USE
        blis_gemm(m,
                  n,
                  k,
                  alphap,
                  tableA,
                  tableB,
                  C,
                  posWset2);
#else
        fprintf(stderr,"ERROR: Invalid state, cannot use blis without defining it\n");
        exit(1);
#endif
    }
    else
    {
        gemm(m,
             n,
             k,
	    	 alphap,
             A,
             m,
             tableB,
             k,
             C,
             m,
             Ac_pack_v,
             Bc_pack_v);
    }
    time003 = gettime(); //004 after gemm
    gemm_time += time003 - time002;

    //for(int i=0; i < m; i++)
    //{
    //    if((tableA_bitcount[i] == 0) || (tableA_bitcount[i] == (uint32_t) snp_size))
    //        printf("Ind: %d\t Pos: %d\tBC: %u\n",i,tableA_posIndex[i],tableA_bitcount[i]);
    //}

	get_pairwise_ld_score (tableA_bitcount,
                           tableB_bitcount,
                           C,
                           m,
                           n,
                           snp_size,
                           &results);

    time004 = gettime(); //004 after ldScore
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
	write_time+=time005 - time004;
    gettimeofday(&end, NULL);

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
	fprintf(threadData[0].threadLog,"Table A size: %d, Table B size: %d, \
Compressed snip size: %d, Time: %5.6fs, GFLOPS: %3.3f\n",m, n, k, time_ref, GFLOPS_BLIS);
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
	free(Ac_pack_v);
	free(Bc_pack_v);
    free(A);
    free(C);
    free(results);
#ifdef COUNTSTATES
    free(STtable);
#endif
}
