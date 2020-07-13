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

extern scalm_t*   scalm_cntl;

blksz_t*          gemm4m1_mc;
blksz_t*          gemm4m1_nc;
blksz_t*          gemm4m1_kc;
blksz_t*          gemm4m1_mr;
blksz_t*          gemm4m1_nr;
blksz_t*          gemm4m1_kr;

func_t*           gemm4m1_ukrs;

packm_t*          gemm4m1_packa_cntl;
packm_t*          gemm4m1_packb_cntl;

gemm_t*           gemm4m1_cntl_bp_ke;
gemm_t*           gemm4m1_cntl_op_bp;
gemm_t*           gemm4m1_cntl_mm_op;
gemm_t*           gemm4m1_cntl_vl_mm;

gemm_t*           gemm4m1_cntl;


void bli_gemm4m1_cntl_init()
{
	// Create blocksize objects for each dimension.
	// NOTE: the complex blocksizes for 4m1 are generally equal to their
	// corresponding real domain counterparts. However, we want to promote
	// similar cache footprints for the micro-panels of A and B (when
	// compared to executing in the real domain), and since the complex
	// micro-panels are twice as "fat" (due to storing real and imaginary
	// parts), we reduce KC by a factor of 2 to compensate.
	gemm4m1_mc
	=
	bli_blksz_obj_create( 0,                 0,
	                      0,                 0,
	                      BLIS_DEFAULT_MC_S, BLIS_MAXIMUM_MC_S,
	                      BLIS_DEFAULT_MC_D, BLIS_MAXIMUM_MC_D );
	gemm4m1_nc
	=
	bli_blksz_obj_create( 0,                 0,
	                      0,                 0,
	                      BLIS_DEFAULT_NC_S, BLIS_MAXIMUM_NC_S,
	                      BLIS_DEFAULT_NC_D, BLIS_MAXIMUM_NC_D );
	gemm4m1_kc
	=
	bli_blksz_obj_create( 0,                   0,
	                      0,                   0,
	                      BLIS_DEFAULT_KC_S/2, BLIS_MAXIMUM_KC_S/2,
	                      BLIS_DEFAULT_KC_D/2, BLIS_MAXIMUM_KC_D/2 );
	gemm4m1_mr
	=
	bli_blksz_obj_create( 0,                 0,
	                      0,                 0,
	                      BLIS_DEFAULT_MR_S, BLIS_PACKDIM_MR_S,
	                      BLIS_DEFAULT_MR_D, BLIS_PACKDIM_MR_D );
	gemm4m1_nr
	=
	bli_blksz_obj_create( 0,                 0,
	                      0,                 0,
	                      BLIS_DEFAULT_NR_S, BLIS_PACKDIM_NR_S,
	                      BLIS_DEFAULT_NR_D, BLIS_PACKDIM_NR_D );
	gemm4m1_kr
	=
	bli_blksz_obj_create( 0,                 0,
	                      0,                 0,
	                      BLIS_DEFAULT_KR_S, BLIS_PACKDIM_KR_S,
	                      BLIS_DEFAULT_KR_D, BLIS_PACKDIM_KR_D );


	// Attach the register blksz_t objects as blocksize multiples to the cache
	// blksz_t objects.
	bli_blksz_obj_attach_mult_to( gemm4m1_mr, gemm4m1_mc );
	bli_blksz_obj_attach_mult_to( gemm4m1_nr, gemm4m1_nc );
	bli_blksz_obj_attach_mult_to( gemm4m1_kr, gemm4m1_kc );


	// Attach the mr and nr blksz_t objects to each cache blksz_t object.
	// The primary example of why this is needed relates to nudging kc.
	// In hemm, symm, trmm, or trmm3, we need to know both mr and nr,
	// since the multiple we target in nudging depends on whether the
	// structured matrix is on the left or the right.
	bli_blksz_obj_attach_mr_nr_to( gemm4m1_mr, gemm4m1_nr, gemm4m1_mc );
	bli_blksz_obj_attach_mr_nr_to( gemm4m1_mr, gemm4m1_nr, gemm4m1_nc );
	bli_blksz_obj_attach_mr_nr_to( gemm4m1_mr, gemm4m1_nr, gemm4m1_kc );


	// Create function pointer object for each datatype-specific gemm
	// micro-kernel.
	gemm4m1_ukrs
	=
	bli_func_obj_create(
	    NULL,                 FALSE,
	    NULL,                 FALSE,
	    BLIS_CGEMM4M1_UKERNEL, BLIS_CGEMM4M1_UKERNEL_PREFERS_CONTIG_ROWS,
	    BLIS_ZGEMM4M1_UKERNEL, BLIS_ZGEMM4M1_UKERNEL_PREFERS_CONTIG_ROWS );


	// Create control tree objects for packm operations.
	gemm4m1_packa_cntl
	=
	bli_packm_cntl_obj_create( BLIS_BLOCKED,
	                           BLIS_VARIANT1,
	                           gemm4m1_mr,
	                           gemm4m1_kr,
	                           FALSE, // do NOT invert diagonal
	                           FALSE, // reverse iteration if upper?
	                           FALSE, // reverse iteration if lower?
	                           BLIS_PACKED_ROW_PANELS_4MI,
	                           BLIS_BUFFER_FOR_A_BLOCK );

	gemm4m1_packb_cntl
	=
	bli_packm_cntl_obj_create( BLIS_BLOCKED,
	                           BLIS_VARIANT1,
	                           gemm4m1_kr,
	                           gemm4m1_nr,
	                           FALSE, // do NOT invert diagonal
	                           FALSE, // reverse iteration if upper?
	                           FALSE, // reverse iteration if lower?
	                           BLIS_PACKED_COL_PANELS_4MI,
	                           BLIS_BUFFER_FOR_B_PANEL );


	//
	// Create a control tree for packing A and B, and streaming C.
	//

	// Create control tree object for lowest-level block-panel kernel.
	gemm4m1_cntl_bp_ke
	=
	bli_gemm_cntl_obj_create( BLIS_UNB_OPT,
	                          BLIS_VARIANT2,
	                          NULL,
	                          gemm4m1_ukrs,
	                          NULL, NULL, NULL,
	                          NULL, NULL, NULL );

	// Create control tree object for outer panel (to block-panel)
	// problem.
	gemm4m1_cntl_op_bp
	=
	bli_gemm_cntl_obj_create( BLIS_BLOCKED,
	                          BLIS_VARIANT1,
	                          gemm4m1_mc,
	                          NULL,
	                          NULL,
	                          gemm4m1_packa_cntl,
	                          gemm4m1_packb_cntl,
	                          NULL,
	                          gemm4m1_cntl_bp_ke,
	                          NULL );

	// Create control tree object for general problem via multiple
	// rank-k (outer panel) updates.
	gemm4m1_cntl_mm_op
	=
	bli_gemm_cntl_obj_create( BLIS_BLOCKED,
	                          BLIS_VARIANT3,
	                          gemm4m1_kc,
	                          NULL,
	                          NULL,
	                          NULL,
	                          NULL,
	                          NULL,
	                          gemm4m1_cntl_op_bp,
	                          NULL );

	// Create control tree object for very large problem via multiple
	// general problems.
	gemm4m1_cntl_vl_mm
	=
	bli_gemm_cntl_obj_create( BLIS_BLOCKED,
	                          BLIS_VARIANT2,
	                          gemm4m1_nc,
	                          NULL,
	                          NULL,
	                          NULL,
	                          NULL,
	                          NULL,
	                          gemm4m1_cntl_mm_op,
	                          NULL );

	// Alias the "master" gemm control tree to a shorter name.
	gemm4m1_cntl = gemm4m1_cntl_vl_mm;

}

void bli_gemm4m1_cntl_finalize()
{
	bli_blksz_obj_free( gemm4m1_mc );
	bli_blksz_obj_free( gemm4m1_nc );
	bli_blksz_obj_free( gemm4m1_kc );
	bli_blksz_obj_free( gemm4m1_mr );
	bli_blksz_obj_free( gemm4m1_nr );
	bli_blksz_obj_free( gemm4m1_kr );

	bli_func_obj_free( gemm4m1_ukrs );

	bli_cntl_obj_free( gemm4m1_packa_cntl );
	bli_cntl_obj_free( gemm4m1_packb_cntl );

	bli_cntl_obj_free( gemm4m1_cntl_bp_ke );
	bli_cntl_obj_free( gemm4m1_cntl_op_bp );
	bli_cntl_obj_free( gemm4m1_cntl_mm_op );
	bli_cntl_obj_free( gemm4m1_cntl_vl_mm );
}

