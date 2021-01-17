#include "blis.h"
#include "immintrin.h"

void bin_gemm_2x2(
                        dim_t              k,
                        double* restrict   alpha,
                        double* restrict   a,
                        double* restrict   b,
                        double* restrict   beta,
                        double* restrict   c, inc_t rs_c, inc_t cs_c,
                        auxinfo_t*         data
                      )
{
  	dim_t   k_iter = k / 4;
	dim_t   k_left = k % 4;

	register uint64_t *a_ptr, *b_ptr;
	register uint64_t p, q, r;
	register uint64_t t00, t10, t11, t01;

	a_ptr = (uint64_t*)a;
	b_ptr = (uint64_t*)b;

	/*	printf("%x %x == %llu\n%x %x == %llu\n", 
	       a_ptr[0], b_ptr[0],_mm_popcnt_u64(a_ptr[0] & b_ptr[0]), 
	       a_ptr[2], b_ptr[2],_mm_popcnt_u64(a_ptr[2] & b_ptr[2]) );
	*/
	t00 = 0; t01 = 0;
	t10 = 0; t11 = 0;

	for (int i = 0; i != k_iter; ++i){
	  p = (*a_ptr);
	  q = (*b_ptr);
	  r = (*(b_ptr+1));

	  t00 += _mm_popcnt_u64(p & q);
	  t01 += _mm_popcnt_u64(p & r);

	  p = (*(a_ptr+1));

	  t10 += _mm_popcnt_u64(p & q);
	  t11 += _mm_popcnt_u64(p & r);


	  p = (*(a_ptr+2));
	  q = (*(b_ptr+2));
	  r = (*(b_ptr+3));

	  t00 += _mm_popcnt_u64(p & q);
	  t01 += _mm_popcnt_u64(p & r);

	  p = (*(a_ptr+3));

	  t10 += _mm_popcnt_u64(p & q);
	  t11 += _mm_popcnt_u64(p & r);



	  p = (*(a_ptr+4));
	  q = (*(b_ptr+4));
	  r = (*(b_ptr+5));

	  t00 += _mm_popcnt_u64(p & q);
	  t01 += _mm_popcnt_u64(p & r);

	  p = (*(a_ptr+5));

	  t10 += _mm_popcnt_u64(p & q);
	  t11 += _mm_popcnt_u64(p & r);



	  p = (*(a_ptr+6));
	  q = (*(b_ptr+6));
	  r = (*(b_ptr+7));

	  t00 += _mm_popcnt_u64(p & q);
	  t01 += _mm_popcnt_u64(p & r);

	  p = (*(a_ptr+7));

	  t10 += _mm_popcnt_u64(p & q);
	  t11 += _mm_popcnt_u64(p & r);

	  a_ptr += 8;
	  b_ptr += 8;
	}


	for (int i = 0; i != k_left; ++i){
	  p = (*a_ptr);
	  q = (*b_ptr);
	  r = (*(b_ptr+1));

	  //	  printf("--%x--\n", p);

	  t00 += _mm_popcnt_u64(p & q);
	  t01 += _mm_popcnt_u64(p & r);

	  p = (*(a_ptr+1));

	  t10 += _mm_popcnt_u64(p & q);
	  t11 += _mm_popcnt_u64(p & r);

	  a_ptr += 2;
	  b_ptr += 2;


	  //	  printf("%llu %llu\n%llu %llu\n\n", t00, t01, t10, t11);
	}

	/*

	*(c)  += (*alpha) * t00;
	*((c + rs_c))  += (*alpha) * t01;
	*((c + cs_c))  += (*alpha) * t10;
	*((c + rs_c + cs_c))  += (*alpha) * t11;
	*/

    //printf("Before\tc:%lu \tt00:%lu\tt10:%lu\tt01:%lu\tt11:%lu\n",
    //      *((uint64_t*)c),t00,t10,t01,t11);

    // C = (a *) AB (+ bC)

    if((*beta) == 0.0)
    {
	    *((uint64_t*)c)  =  t00;
	    *((uint64_t*)(c + rs_c))  = t10;
	    *((uint64_t*)(c + cs_c))  = t01;
	    *((uint64_t*)(c + rs_c + cs_c))  = t11;
    }
    else
    {
	    *((uint64_t*)c)  +=  t00;
	    *((uint64_t*)(c + rs_c))  += t10;
	    *((uint64_t*)(c + cs_c))  += t01;
	    *((uint64_t*)(c + rs_c + cs_c))  += t11; 
    }
    //printf("After\tc:%lu\tt00:%lu\tt10:%lu\tt01:%lu\tt11:%lu\n",
    //      *((uint64_t*)c),t00,t10,t01,t11);
}
