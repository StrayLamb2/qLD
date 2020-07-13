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

#define FUNCPTR_T normfv_fp

typedef void (*FUNCPTR_T)(
                           dim_t   m,
                           void*   x, inc_t incx,
                           void*   norm
                         );

static FUNCPTR_T GENARRAY(ftypes,normfv_unb_var1);


void bli_normfv_unb_var1( obj_t*  x,
                          obj_t*  norm )
{
	num_t     dt_x      = bli_obj_datatype( *x );

	dim_t     m         = bli_obj_vector_dim( *x );

	inc_t     incx      = bli_obj_vector_inc( *x );
	void*     buf_x     = bli_obj_buffer_at_off( *x );

	void*     buf_norm  = bli_obj_buffer_at_off( *norm );

	FUNCPTR_T f;

	// Index into the type combination array to extract the correct
	// function pointer.
	f = ftypes[dt_x];

	// Invoke the function.
	f( m,
	   buf_x, incx,
	   buf_norm );
}


#undef  GENTFUNCR
#define GENTFUNCR( ctype_x, ctype_xr, chx, chxr, varname, kername ) \
\
void PASTEMAC(chx,varname)( \
                            dim_t   m, \
                            void*   x, inc_t incx, \
                            void*   norm \
                          ) \
{ \
	ctype_x*  x_cast     = x; \
	ctype_xr* norm_cast  = norm; \
	ctype_xr* zero       = PASTEMAC(chxr,0); \
	ctype_xr* one        = PASTEMAC(chxr,1); \
	ctype_xr  scale; \
	ctype_xr  sumsq; \
	ctype_xr  sqrt_sumsq; \
\
	/* Return a norm of zero if either dimension is zero. */ \
	if ( bli_zero_dim1( m ) ) \
	{ \
		PASTEMAC(chxr,set0s)( *norm_cast ); \
		return; \
	} \
\
	/* Initialize scale and sumsq to begin the summation. */ \
	PASTEMAC2(chxr,chxr,copys)( *zero, scale ); \
	PASTEMAC2(chxr,chxr,copys)( *one,  sumsq ); \
\
	/* Compute the sum of the squares of the vector. */ \
	PASTEMAC(chx,kername)( m, \
	                       x_cast, incx, \
	                       &scale, \
	                       &sumsq ); \
\
	/* Compute: norm = scale * sqrt( sumsq ) */ \
	PASTEMAC2(chxr,chxr,sqrt2s)( sumsq, sqrt_sumsq ); \
	PASTEMAC2(chxr,chxr,scals)( scale, sqrt_sumsq ); \
\
	/* Store the final value to the output variable. */ \
	PASTEMAC2(chxr,chxr,copys)( sqrt_sumsq, *norm_cast ); \
}

INSERT_GENTFUNCR_BASIC( normfv_unb_var1, sumsqv_unb_var1 )

