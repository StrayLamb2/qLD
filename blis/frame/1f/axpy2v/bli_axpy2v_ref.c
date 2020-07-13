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

/*
#define FUNCPTR_T axpy2v_fp

typedef void (*FUNCPTR_T)(
                           conj_t conjx,
                           conj_t conjy,
                           dim_t  n,
                           void*  alpha1,
                           void*  alpha2,
                           void*  x, inc_t incx,
                           void*  y, inc_t incy,
                           void*  z, inc_t incz
                         );

// If some mixed datatype functions will not be compiled, we initialize
// the corresponding elements of the function array to NULL.
#ifdef BLIS_ENABLE_MIXED_PRECISION_SUPPORT
static FUNCPTR_T GENARRAY3_ALL(ftypes,axpy2v_ref);
#else
#ifdef BLIS_ENABLE_MIXED_DOMAIN_SUPPORT
static FUNCPTR_T GENARRAY3_EXT(ftypes,axpy2v_ref);
#else
static FUNCPTR_T GENARRAY3_MIN(ftypes,axpy2v_ref);
#endif
#endif


void bli_axpy2v_ref( obj_t*  alpha1,
                          obj_t*  alpha2,
                          obj_t*  x,
                          obj_t*  y,
                          obj_t*  z )
{
	num_t     dt_x      = bli_obj_datatype( *x );
	num_t     dt_y      = bli_obj_datatype( *y );

	conj_t    conjx     = bli_obj_conj_status( *x );
	conj_t    conjy     = bli_obj_conj_status( *y );
	dim_t     n         = bli_obj_vector_dim( *x );

	inc_t     inc_x     = bli_obj_vector_inc( *x );
	void*     buf_x     = bli_obj_buffer_at_off( *x );

	inc_t     inc_y     = bli_obj_vector_inc( *y );
	void*     buf_y     = bli_obj_buffer_at_off( *y );

	inc_t     inc_z     = bli_obj_vector_inc( *z );
	void*     buf_z     = bli_obj_buffer_at_off( *z );

	num_t     dt_alpha1;
	void*     buf_alpha1;

	num_t     dt_alpha2;
	void*     buf_alpha2;

	FUNCPTR_T f;

	// If alpha is a scalar constant, use dt_x to extract the address of the
	// corresponding constant value; otherwise, use the datatype encoded
	// within the alpha object and extract the buffer at the alpha offset.
	bli_set_scalar_dt_buffer( alpha1, dt_x, dt_alpha1, buf_alpha1 );
	bli_set_scalar_dt_buffer( alpha2, dt_x, dt_alpha2, buf_alpha2 );

	// Index into the type combination array to extract the correct
	// function pointer.
	f = ftypes[dt_alpha1][dt_x][dt_y];

	// Invoke the function.
	f( conjx,
	   conjy,
	   n,
	   buf_alpha1,
	   buf_alpha2,
	   buf_x, inc_x,
	   buf_y, inc_y,
	   buf_z, inc_z );
}
*/


#undef  GENTFUNC3U12
#define GENTFUNC3U12( ctype_x, ctype_y, ctype_z, ctype_xy, chx, chy, chz, chxy, varname, kername ) \
\
void PASTEMAC3(chx,chy,chz,varname) \
     ( \
       conj_t             conjx, \
       conj_t             conjy, \
       dim_t              n, \
       ctype_xy* restrict alpha1, \
       ctype_xy* restrict alpha2, \
       ctype_x*  restrict x, inc_t incx, \
       ctype_y*  restrict y, inc_t incy, \
       ctype_z*  restrict z, inc_t incz  \
     ) \
{ \
	ctype_xy* alpha1_cast = alpha1; \
	ctype_xy* alpha2_cast = alpha2; \
	ctype_x*  x_cast      = x; \
	ctype_y*  y_cast      = y; \
	ctype_z*  z_cast      = z; \
\
	PASTEMAC3(chxy,chx,chz,kername)( conjx, \
	                                 n, \
	                                 alpha1_cast, \
	                                 x_cast, incx, \
	                                 z_cast, incz ); \
	PASTEMAC3(chxy,chy,chz,kername)( conjy, \
	                                 n, \
	                                 alpha2_cast, \
	                                 y_cast, incy, \
	                                 z_cast, incz ); \
}

// Define the basic set of functions unconditionally, and then also some
// mixed datatype functions if requested.
INSERT_GENTFUNC3U12_BASIC( axpy2v_ref, AXPYV_KERNEL )

#ifdef BLIS_ENABLE_MIXED_DOMAIN_SUPPORT
INSERT_GENTFUNC3U12_MIX_D( axpy2v_ref, AXPYV_KERNEL )
#endif

#ifdef BLIS_ENABLE_MIXED_PRECISION_SUPPORT
INSERT_GENTFUNC3U12_MIX_P( axpy2v_ref, AXPYV_KERNEL )
#endif

