/*

   BLIS    
   An object-based framework for developing high-performance BLAS-like
   libraries.

   Copyright (C) 2014, The University of Texas at Austin

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are
   met:
    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    - Neither the name of The University of Texas at Austin nor the names
      of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
   HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "blis.h"

void bli_gemm_blk_var3f( obj_t*  a,
                         obj_t*  b,
                         obj_t*  c,
                         gemm_t* cntl,
                         gemm_thrinfo_t* thread )
{
    obj_t  c_pack_s;
    obj_t  a1_pack_s, b1_pack_s;

    obj_t  a1, b1;
    obj_t* a1_pack = NULL;
    obj_t* b1_pack = NULL;
    obj_t* c_pack = NULL;

	dim_t  i;
	dim_t  b_alg;
	dim_t  k_trans;

    if( thread_am_ochief( thread ) ){
        // Initialize object for packing C
	    bli_obj_init_pack( &c_pack_s );
        bli_packm_init( c, &c_pack_s,
                        cntl_sub_packm_c( cntl ) );

        // Scale C by beta (if instructed).
        bli_scalm_int( &BLIS_ONE,
                       c,
                       cntl_sub_scalm( cntl ) );
    }
    c_pack = thread_obroadcast( thread, &c_pack_s );

    // Initialize pack objects for A and B that are passed into packm_init().
    if( thread_am_ichief( thread ) ){
        bli_obj_init_pack( &a1_pack_s );
        bli_obj_init_pack( &b1_pack_s );
    }
    a1_pack = thread_ibroadcast( thread, &a1_pack_s );
    b1_pack = thread_ibroadcast( thread, &b1_pack_s );

	// Pack C (if instructed).
	bli_packm_int( c, c_pack,
	               cntl_sub_packm_c( cntl ),
                   gemm_thread_sub_opackm( thread ) );

	// Query dimension in partitioning direction.
	k_trans = bli_obj_width_after_trans( *a );

	// Partition along the k dimension.
	for ( i = 0; i < k_trans; i += b_alg )
	{
		// Determine the current algorithmic blocksize.
		// NOTE: We call a gemm/hemm/symm-specific function to determine
		// the kc blocksize so that we can implement the "nudging" of kc
		// to be a multiple of mr or nr, as needed.
		b_alg = bli_gemm_determine_kc_f( i, k_trans, a, b,
		                                 cntl_blocksize( cntl ) );

		// Acquire partitions for A1 and B1.
		bli_acquire_mpart_l2r( BLIS_SUBPART1,
		                       i, b_alg, a, &a1 );
		bli_acquire_mpart_t2b( BLIS_SUBPART1,
		                       i, b_alg, b, &b1 );

		// Initialize objects for packing A1 and B1.
        if( thread_am_ichief( thread ) ) {
            bli_packm_init( &a1, a1_pack,
                            cntl_sub_packm_a( cntl ) );
            bli_packm_init( &b1, b1_pack,
                            cntl_sub_packm_b( cntl ) );
        }
        thread_ibarrier( thread );

		// Pack A1 (if instructed).
		bli_packm_int( &a1, a1_pack,
		               cntl_sub_packm_a( cntl ),
                       gemm_thread_sub_ipackm( thread ) );

		// Pack B1 (if instructed).
		bli_packm_int( &b1, b1_pack,
		               cntl_sub_packm_b( cntl ),
                       gemm_thread_sub_ipackm( thread ) );

		// Perform gemm subproblem.
		bli_gemm_int( &BLIS_ONE,
		              a1_pack,
		              b1_pack,
		              &BLIS_ONE,
		              c_pack,
		              cntl_sub_gemm( cntl ),
                      gemm_thread_sub_gemm( thread) );

		// This variant executes multiple rank-k updates. Therefore, if the
		// internal beta scalar on matrix C is non-zero, we must use it
		// only for the first iteration (and then BLIS_ONE for all others).
		// And since c_pack is a local obj_t, we can simply overwrite the
		// internal beta scalar with BLIS_ONE once it has been used in the
		// first iteration.
        thread_ibarrier( thread );
		if ( i == 0 && thread_am_ichief( thread ) ) bli_obj_scalar_reset( c_pack );

	}

    thread_obarrier( thread );

	// Unpack C (if C was packed).
    bli_unpackm_int( c_pack, c,
                     cntl_sub_unpackm_c( cntl ),
                     gemm_thread_sub_opackm( thread ) );

	// If any packing buffers were acquired within packm, release them back
	// to the memory manager.
    if( thread_am_ochief( thread ) )
        bli_packm_release( c_pack, cntl_sub_packm_c( cntl ) );
    if( thread_am_ichief( thread ) ){
        bli_packm_release( a1_pack, cntl_sub_packm_a( cntl ) );
        bli_packm_release( b1_pack, cntl_sub_packm_b( cntl ) );
    }
}

