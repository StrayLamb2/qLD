//MPAMPIS -not included before-
#include "bin/gpu_kernel/blislike.h"

__kernel void blis_like (
    unsigned int k, unsigned int m_alg, unsigned int n_alg,
    __global unsigned int *a,
    __global unsigned int *b,
    __global unsigned int *c,
    unsigned int rs_c, unsigned int cs_c
) {
  const size_t bx = get_group_id(0);
  const size_t by = get_group_id(1);
  const size_t tx = get_local_id(0);
  const size_t ty = get_local_id(1);

  unsigned int a_offset = bx * BLOCK_SIZE_X;
  unsigned int b_offset = by * BLOCK_SIZE_Y;

  unsigned int ir = bx * BLOCK_SIZE_X / GPU_BLOCK_MR * GPU_BLOCK_MR;
  unsigned int jr = by * BLOCK_SIZE_Y / GPU_BLOCK_NR * GPU_BLOCK_MR;

  unsigned int nr_alg = min((uint) GPU_BLOCK_NR, n_alg - jr);
  unsigned int mr_alg = min((uint) GPU_BLOCK_MR, m_alg - ir);

  unsigned int nr_block = min((uint) BLOCK_SIZE_Y, n_alg - b_offset);
  unsigned int mr_block = min((uint) BLOCK_SIZE_X, m_alg - a_offset);

	__global unsigned int *Ar = &a[ir*k + (a_offset % BLOCK_SIZE_X)];
	__global unsigned int *Br = &b[jr*k + (b_offset % BLOCK_SIZE_Y)];
	__global unsigned int *Cr = &c[a_offset+b_offset*cs_c];

  // can this be shared between threads in a warp??  does it even benefit us?
  __local unsigned int ab[BLOCK_SIZE_Y*BLOCK_SIZE_X];
  __local unsigned int *abij;

  unsigned int bj, ai;
  unsigned int l, i, j;
  for (j = ty; j < nr_block; j += LOCAL_1) {
    abij = &ab[j*mr_block];
    for (i = tx; i < mr_block; i += LOCAL_0) {
      abij[i] = 0;
    }
  }

  // __global unsigned int *Cr_l;
  for (l = 0; l < k; ++l) {
    abij=ab;
    for (j = ty; j < nr_block; j += LOCAL_1) {
      bj = Br[j];      
      abij = &ab[j*mr_block];
      for (i = tx; i < mr_block; i += LOCAL_0) {
        ai = Ar[i];
        abij[i] += popcount(ai & bj);        
      }
    }
    Ar += mr_alg;
    Br += nr_alg;
  }

  __global unsigned int *Cr_l;
  for (j = ty; j < nr_block; j += LOCAL_1) {
    Cr_l = &Cr[j*cs_c];
    abij = &ab[j*mr_block];
    for (i = tx; i < mr_block; i += LOCAL_0) {
      // Cr_l[i*rs_c] += abij[i];
      // since rs_c is assumed to be 1...
      Cr_l[i] += abij[i];
    }
  }
}

__kernel void blis_like_empty (
    unsigned int k, unsigned int mr_alg, unsigned int nr_alg,
    __global unsigned int *a,
    __global unsigned int *b,
    __global unsigned int *c,
    unsigned int rs_c, unsigned int cs_c
) {
    const size_t by = get_group_id(1);
    const size_t ty = get_local_id(1);  // should be zero?
    const size_t bx = get_group_id(0);
    const size_t tx = get_local_id(0);

    const size_t a_offset = by * GPU_BLOCK_MC + ty;
    const size_t b_offset = bx * GPU_BLOCK_NC + tx;

    // a += a_offset;
    // b += b_offset;
    // is this correct?
    // c += b_offset * cs_c;

    // unsigned int l, i;
    // unsigned int bj, ai;
    // for (l = 0; l < k; ++l) {
    //   bj = ~(*b);
    //   for (i = 0; i < mr_alg; ++i) {
    //     ai = *(a + i);
    //
    //     *(c + i) += popcount(ai & bj);
    //   }
    //   a += mr_alg;
    //   b += nr_alg;
    // }


}

// AxB
__kernel void blis_like8x4 (
    unsigned int k, unsigned int m_alg, unsigned int n_alg,
    __global unsigned int *a,
    __global unsigned int *b,
    __global unsigned int *c,
    unsigned int rs_c, unsigned int cs_c
) {
  const size_t bx = get_group_id(0);
  const size_t by = get_group_id(1);
  const size_t tx = get_local_id(0);
  const size_t ty = get_local_id(1);

  unsigned int ir = bx * BLOCK_SIZE_X / GPU_BLOCK_MR * GPU_BLOCK_MR;
  unsigned int jr = by * BLOCK_SIZE_Y / GPU_BLOCK_NR * GPU_BLOCK_NR;

  unsigned int mr_alg = min((uint) GPU_BLOCK_MR, m_alg - ir);
  unsigned int nr_alg = min((uint) GPU_BLOCK_NR, n_alg - jr);

  const size_t a_offset = bx * BLOCK_SIZE_X;
  const size_t b_offset = by * BLOCK_SIZE_Y;

  a += ir*k + (a_offset % GPU_BLOCK_MR) + tx;
  b += jr*k + (b_offset % GPU_BLOCK_NR) + ty;
  c += cs_c * (b_offset + ty) + (a_offset + tx);

  unsigned int a0, a1, a2, a3, a4, a5, a6, a7;
  unsigned int b0, b1, b2, b3;

  unsigned int t00, t01, t02, t03;
  unsigned int t10, t11, t12, t13;
  unsigned int t20, t21, t22, t23;
  unsigned int t30, t31, t32, t33;
  unsigned int t40, t41, t42, t43;
  unsigned int t50, t51, t52, t53;
  unsigned int t60, t61, t62, t63;
  unsigned int t70, t71, t72, t73;

  unsigned int c00, c01, c02, c03;
  unsigned int c10, c11, c12, c13;
  unsigned int c20, c21, c22, c23;
  unsigned int c30, c31, c32, c33;
  unsigned int c40, c41, c42, c43;
  unsigned int c50, c51, c52, c53;
  unsigned int c60, c61, c62, c63;
  unsigned int c70, c71, c72, c73;

  c00 = c01 = c02 = c03 = 0;
  c10 = c11 = c12 = c13 = 0;
  c20 = c21 = c22 = c23 = 0;
  c30 = c31 = c32 = c33 = 0;
  c40 = c41 = c42 = c43 = 0;
  c50 = c51 = c52 = c53 = 0;
  c60 = c61 = c62 = c63 = 0;
  c70 = c71 = c72 = c73 = 0;

  unsigned int l;
  for (l = 0; l < k; l++) {

    a0 = *(a + 0*BLOCK_SIZE_X/8);
    a1 = *(a + 1*BLOCK_SIZE_X/8);
    a2 = *(a + 2*BLOCK_SIZE_X/8);
    a3 = *(a + 3*BLOCK_SIZE_X/8);
    a4 = *(a + 4*BLOCK_SIZE_X/8);
    a5 = *(a + 5*BLOCK_SIZE_X/8);
    a6 = *(a + 6*BLOCK_SIZE_X/8);
    a7 = *(a + 7*BLOCK_SIZE_X/8);

    b0 = *(b + 0*BLOCK_SIZE_Y/4);
    b1 = *(b + 1*BLOCK_SIZE_Y/4);
    b2 = *(b + 2*BLOCK_SIZE_Y/4);
    b3 = *(b + 3*BLOCK_SIZE_Y/4);

    a += mr_alg;
    b += nr_alg;

    t00 = a0 & b0;
    t10 = a1 & b0;
    t20 = a2 & b0;
    t30 = a3 & b0;
    t40 = a4 & b0;
    t50 = a5 & b0;
    t60 = a6 & b0;
    t70 = a7 & b0;

    t01 = a0 & b1;
    t11 = a1 & b1;
    t21 = a2 & b1;
    t31 = a3 & b1;
    t41 = a4 & b1;
    t51 = a5 & b1;
    t61 = a6 & b1;
    t71 = a7 & b1;

    t02 = a0 & b2;
    t12 = a1 & b2;
    t22 = a2 & b2;
    t32 = a3 & b2;
    t42 = a4 & b2;
    t52 = a5 & b2;
    t62 = a6 & b2;
    t72 = a7 & b2;

    t03 = a0 & b3;
    t13 = a1 & b3;
    t23 = a2 & b3;
    t33 = a3 & b3;
    t43 = a4 & b3;
    t53 = a5 & b3;
    t63 = a6 & b3;
    t73 = a7 & b3;


    c00 += popcount(t00);
    c10 += popcount(t10);
    c20 += popcount(t20);
    c30 += popcount(t30);
    c40 += popcount(t40);
    c50 += popcount(t50);
    c60 += popcount(t60);
    c70 += popcount(t70);

    c01 += popcount(t01);
    c11 += popcount(t11);
    c21 += popcount(t21);
    c31 += popcount(t31);
    c41 += popcount(t41);
    c51 += popcount(t51);
    c61 += popcount(t61);
    c71 += popcount(t71);

    c02 += popcount(t02);
    c12 += popcount(t12);
    c22 += popcount(t22);
    c32 += popcount(t32);
    c42 += popcount(t42);
    c52 += popcount(t52);
    c62 += popcount(t62);
    c72 += popcount(t72);

    c03 += popcount(t03);
    c13 += popcount(t13);
    c23 += popcount(t23);
    c33 += popcount(t33);
    c43 += popcount(t43);
    c53 += popcount(t53);
    c63 += popcount(t63);
    c73 += popcount(t73);
  }

  *(c + 0*BLOCK_SIZE_X/8 + cs_c*0*BLOCK_SIZE_Y/4) += c00;
  *(c + 1*BLOCK_SIZE_X/8 + cs_c*0*BLOCK_SIZE_Y/4) += c10;
  *(c + 2*BLOCK_SIZE_X/8 + cs_c*0*BLOCK_SIZE_Y/4) += c20;
  *(c + 3*BLOCK_SIZE_X/8 + cs_c*0*BLOCK_SIZE_Y/4) += c30;
  *(c + 4*BLOCK_SIZE_X/8 + cs_c*0*BLOCK_SIZE_Y/4) += c40;
  *(c + 5*BLOCK_SIZE_X/8 + cs_c*0*BLOCK_SIZE_Y/4) += c50;
  *(c + 6*BLOCK_SIZE_X/8 + cs_c*0*BLOCK_SIZE_Y/4) += c60;
  *(c + 7*BLOCK_SIZE_X/8 + cs_c*0*BLOCK_SIZE_Y/4) += c70;

  *(c + 0*BLOCK_SIZE_X/8 + cs_c*1*BLOCK_SIZE_Y/4) += c01;
  *(c + 1*BLOCK_SIZE_X/8 + cs_c*1*BLOCK_SIZE_Y/4) += c11;
  *(c + 2*BLOCK_SIZE_X/8 + cs_c*1*BLOCK_SIZE_Y/4) += c21;
  *(c + 3*BLOCK_SIZE_X/8 + cs_c*1*BLOCK_SIZE_Y/4) += c31;
  *(c + 4*BLOCK_SIZE_X/8 + cs_c*1*BLOCK_SIZE_Y/4) += c41;
  *(c + 5*BLOCK_SIZE_X/8 + cs_c*1*BLOCK_SIZE_Y/4) += c51;
  *(c + 6*BLOCK_SIZE_X/8 + cs_c*1*BLOCK_SIZE_Y/4) += c61;
  *(c + 7*BLOCK_SIZE_X/8 + cs_c*1*BLOCK_SIZE_Y/4) += c71;

  *(c + 0*BLOCK_SIZE_X/8 + cs_c*2*BLOCK_SIZE_Y/4) += c02;
  *(c + 1*BLOCK_SIZE_X/8 + cs_c*2*BLOCK_SIZE_Y/4) += c12;
  *(c + 2*BLOCK_SIZE_X/8 + cs_c*2*BLOCK_SIZE_Y/4) += c22;
  *(c + 3*BLOCK_SIZE_X/8 + cs_c*2*BLOCK_SIZE_Y/4) += c32;
  *(c + 4*BLOCK_SIZE_X/8 + cs_c*2*BLOCK_SIZE_Y/4) += c42;
  *(c + 5*BLOCK_SIZE_X/8 + cs_c*2*BLOCK_SIZE_Y/4) += c52;
  *(c + 6*BLOCK_SIZE_X/8 + cs_c*2*BLOCK_SIZE_Y/4) += c62;
  *(c + 7*BLOCK_SIZE_X/8 + cs_c*2*BLOCK_SIZE_Y/4) += c72;

  *(c + 0*BLOCK_SIZE_X/8 + cs_c*3*BLOCK_SIZE_Y/4) += c03;
  *(c + 1*BLOCK_SIZE_X/8 + cs_c*3*BLOCK_SIZE_Y/4) += c13;
  *(c + 2*BLOCK_SIZE_X/8 + cs_c*3*BLOCK_SIZE_Y/4) += c23;
  *(c + 3*BLOCK_SIZE_X/8 + cs_c*3*BLOCK_SIZE_Y/4) += c33;
  *(c + 4*BLOCK_SIZE_X/8 + cs_c*3*BLOCK_SIZE_Y/4) += c43;
  *(c + 5*BLOCK_SIZE_X/8 + cs_c*3*BLOCK_SIZE_Y/4) += c53;
  *(c + 6*BLOCK_SIZE_X/8 + cs_c*3*BLOCK_SIZE_Y/4) += c63;
  *(c + 7*BLOCK_SIZE_X/8 + cs_c*3*BLOCK_SIZE_Y/4) += c73;
}

// now A is associated with Y, and B with X
__kernel void blis_like4x8 (
    unsigned int k, unsigned int m_alg, unsigned int n_alg,
    __global unsigned int *a,
    __global unsigned int *b,
    __global unsigned int *c,
    unsigned int rs_c, unsigned int cs_c
) {
  const size_t bx = get_group_id(0);
  const size_t by = get_group_id(1);
  const size_t tx = get_local_id(0);
  const size_t ty = get_local_id(1);

  unsigned int ir = by * BLOCK_SIZE_Y / GPU_BLOCK_MR * GPU_BLOCK_MR;
  unsigned int jr = bx * BLOCK_SIZE_X / GPU_BLOCK_NR * GPU_BLOCK_NR;

  unsigned int mr_alg = min((uint) GPU_BLOCK_MR, m_alg - ir);
  unsigned int nr_alg = min((uint) GPU_BLOCK_NR, n_alg - jr);

  const size_t a_offset = by * BLOCK_SIZE_Y;
  const size_t b_offset = bx * BLOCK_SIZE_X;

  a += ir*k + (a_offset % GPU_BLOCK_MR) + ty;
  b += jr*k + (b_offset % GPU_BLOCK_NR) + tx;
  c += cs_c * (b_offset + tx) + (a_offset + ty);

    unsigned int a0, a1, a2, a3;
    unsigned int b0, b1, b2, b3, b4, b5, b6, b7;

    unsigned int t00, t01, t02, t03, t04, t05, t06, t07;
    unsigned int t10, t11, t12, t13, t14, t15, t16, t17;
    unsigned int t20, t21, t22, t23, t24, t25, t26, t27;
    unsigned int t30, t31, t32, t33, t34, t35, t36, t37;

    unsigned int c00, c01, c02, c03, c04, c05, c06, c07;
    unsigned int c10, c11, c12, c13, c14, c15, c16, c17;
    unsigned int c20, c21, c22, c23, c24, c25, c26, c27;
    unsigned int c30, c31, c32, c33, c34, c35, c36, c37;

    c00 = c01 = c02 = c03 = c04 = c05 = c06 = c07 = 0;
    c10 = c11 = c12 = c13 = c14 = c15 = c16 = c17 = 0;
    c20 = c21 = c22 = c23 = c24 = c25 = c26 = c27 = 0;
    c30 = c31 = c32 = c33 = c34 = c35 = c36 = c37 = 0;

    unsigned int l;
    for (l = 0; l < k; l++) {
      a0 = *(a + 0*BLOCK_SIZE_Y/4);
      a1 = *(a + 1*BLOCK_SIZE_Y/4);
      a2 = *(a + 2*BLOCK_SIZE_Y/4);
      a3 = *(a + 3*BLOCK_SIZE_Y/4);

      b0 = *(b + 0*BLOCK_SIZE_X/8);
      b1 = *(b + 1*BLOCK_SIZE_X/8);
      b2 = *(b + 2*BLOCK_SIZE_X/8);
      b3 = *(b + 3*BLOCK_SIZE_X/8);
      b4 = *(b + 4*BLOCK_SIZE_X/8);
      b5 = *(b + 5*BLOCK_SIZE_X/8);
      b6 = *(b + 6*BLOCK_SIZE_X/8);
      b7 = *(b + 7*BLOCK_SIZE_X/8);

      a += mr_alg;
      b += nr_alg;

      t00 = a0 & b0;
      t01 = a0 & b1;
      t02 = a0 & b2;
      t03 = a0 & b3;
      t04 = a0 & b4;
      t05 = a0 & b5;
      t06 = a0 & b6;
      t07 = a0 & b7;

      t10 = a1 & b0;
      t11 = a1 & b1;
      t12 = a1 & b2;
      t13 = a1 & b3;
      t14 = a1 & b4;
      t15 = a1 & b5;
      t16 = a1 & b6;
      t17 = a1 & b7;

      t20 = a2 & b0;
      t21 = a2 & b1;
      t22 = a2 & b2;
      t23 = a2 & b3;
      t24 = a2 & b4;
      t25 = a2 & b5;
      t26 = a2 & b6;
      t27 = a2 & b7;

      t30 = a3 & b0;
      t31 = a3 & b1;
      t32 = a3 & b2;
      t33 = a3 & b3;
      t34 = a3 & b4;
      t35 = a3 & b5;
      t36 = a3 & b6;
      t37 = a3 & b7;

      c00 += popcount(t00);
      c01 += popcount(t01);
      c02 += popcount(t02);
      c03 += popcount(t03);
      c04 += popcount(t04);
      c05 += popcount(t05);
      c06 += popcount(t06);
      c07 += popcount(t07);

      c10 += popcount(t10);
      c11 += popcount(t11);
      c12 += popcount(t12);
      c13 += popcount(t13);
      c14 += popcount(t14);
      c15 += popcount(t15);
      c16 += popcount(t16);
      c17 += popcount(t17);

      c20 += popcount(t20);
      c21 += popcount(t21);
      c22 += popcount(t22);
      c23 += popcount(t23);
      c24 += popcount(t24);
      c25 += popcount(t25);
      c26 += popcount(t26);
      c27 += popcount(t27);

      c30 += popcount(t30);
      c31 += popcount(t31);
      c32 += popcount(t32);
      c33 += popcount(t33);
      c34 += popcount(t34);
      c35 += popcount(t35);
      c36 += popcount(t36);
      c37 += popcount(t37);
    }

  *(c + cs_c*0*BLOCK_SIZE_X/8 + 0*BLOCK_SIZE_Y/4) += c00;
  *(c + cs_c*1*BLOCK_SIZE_X/8 + 0*BLOCK_SIZE_Y/4) += c01;
  *(c + cs_c*2*BLOCK_SIZE_X/8 + 0*BLOCK_SIZE_Y/4) += c02;
  *(c + cs_c*3*BLOCK_SIZE_X/8 + 0*BLOCK_SIZE_Y/4) += c03;
  *(c + cs_c*4*BLOCK_SIZE_X/8 + 0*BLOCK_SIZE_Y/4) += c04;
  *(c + cs_c*5*BLOCK_SIZE_X/8 + 0*BLOCK_SIZE_Y/4) += c05;
  *(c + cs_c*6*BLOCK_SIZE_X/8 + 0*BLOCK_SIZE_Y/4) += c06;
  *(c + cs_c*7*BLOCK_SIZE_X/8 + 0*BLOCK_SIZE_Y/4) += c07;

  *(c + cs_c*0*BLOCK_SIZE_X/8 + 1*BLOCK_SIZE_Y/4) += c10;
  *(c + cs_c*1*BLOCK_SIZE_X/8 + 1*BLOCK_SIZE_Y/4) += c11;
  *(c + cs_c*2*BLOCK_SIZE_X/8 + 1*BLOCK_SIZE_Y/4) += c12;
  *(c + cs_c*3*BLOCK_SIZE_X/8 + 1*BLOCK_SIZE_Y/4) += c13;
  *(c + cs_c*4*BLOCK_SIZE_X/8 + 1*BLOCK_SIZE_Y/4) += c14;
  *(c + cs_c*5*BLOCK_SIZE_X/8 + 1*BLOCK_SIZE_Y/4) += c15;
  *(c + cs_c*6*BLOCK_SIZE_X/8 + 1*BLOCK_SIZE_Y/4) += c16;
  *(c + cs_c*7*BLOCK_SIZE_X/8 + 1*BLOCK_SIZE_Y/4) += c17;

  *(c + cs_c*0*BLOCK_SIZE_X/8 + 2*BLOCK_SIZE_Y/4) += c20;
  *(c + cs_c*1*BLOCK_SIZE_X/8 + 2*BLOCK_SIZE_Y/4) += c21;
  *(c + cs_c*2*BLOCK_SIZE_X/8 + 2*BLOCK_SIZE_Y/4) += c22;
  *(c + cs_c*3*BLOCK_SIZE_X/8 + 2*BLOCK_SIZE_Y/4) += c23;
  *(c + cs_c*4*BLOCK_SIZE_X/8 + 2*BLOCK_SIZE_Y/4) += c24;
  *(c + cs_c*5*BLOCK_SIZE_X/8 + 2*BLOCK_SIZE_Y/4) += c25;
  *(c + cs_c*6*BLOCK_SIZE_X/8 + 2*BLOCK_SIZE_Y/4) += c26;
  *(c + cs_c*7*BLOCK_SIZE_X/8 + 2*BLOCK_SIZE_Y/4) += c27;

  *(c + cs_c*0*BLOCK_SIZE_X/8 + 3*BLOCK_SIZE_Y/4) += c30;
  *(c + cs_c*1*BLOCK_SIZE_X/8 + 3*BLOCK_SIZE_Y/4) += c31;
  *(c + cs_c*2*BLOCK_SIZE_X/8 + 3*BLOCK_SIZE_Y/4) += c32;
  *(c + cs_c*3*BLOCK_SIZE_X/8 + 3*BLOCK_SIZE_Y/4) += c33;
  *(c + cs_c*4*BLOCK_SIZE_X/8 + 3*BLOCK_SIZE_Y/4) += c34;
  *(c + cs_c*5*BLOCK_SIZE_X/8 + 3*BLOCK_SIZE_Y/4) += c35;
  *(c + cs_c*6*BLOCK_SIZE_X/8 + 3*BLOCK_SIZE_Y/4) += c36;
  *(c + cs_c*7*BLOCK_SIZE_X/8 + 3*BLOCK_SIZE_Y/4) += c37;

}

// associate X with A, but use a 32x1 local size
__kernel void blis_like4x8v2 (
    unsigned int k, unsigned int m_alg, unsigned int n_alg,
    __global unsigned int *a,
    __global unsigned int *b,
    __global unsigned int *c,
    unsigned int rs_c, unsigned int cs_c
) {
  const size_t bx = get_group_id(0);
  const size_t by = get_group_id(1);
  const size_t tx = get_local_id(0);
  const size_t ty = get_local_id(1);

  unsigned int ir = bx * BLOCK_SIZE_X / GPU_BLOCK_MR * GPU_BLOCK_MR;
  unsigned int jr = by * BLOCK_SIZE_Y / GPU_BLOCK_NR * GPU_BLOCK_NR;

  unsigned int mr_alg = min((uint) GPU_BLOCK_MR, m_alg - ir);
  unsigned int nr_alg = min((uint) GPU_BLOCK_NR, n_alg - jr);

  const size_t a_offset = bx * BLOCK_SIZE_X;
  const size_t b_offset = by * BLOCK_SIZE_Y;

  a += ir*k + (a_offset % GPU_BLOCK_MR) + tx;
  b += jr*k + (b_offset % GPU_BLOCK_NR) + ty;
  c += cs_c * (b_offset + ty) + (a_offset + tx);

    unsigned int a0, a1, a2, a3;
    unsigned int b0, b1, b2, b3, b4, b5, b6, b7;

    unsigned int t00, t01, t02, t03, t04, t05, t06, t07;
    unsigned int t10, t11, t12, t13, t14, t15, t16, t17;
    unsigned int t20, t21, t22, t23, t24, t25, t26, t27;
    unsigned int t30, t31, t32, t33, t34, t35, t36, t37;

    unsigned int c00, c01, c02, c03, c04, c05, c06, c07;
    unsigned int c10, c11, c12, c13, c14, c15, c16, c17;
    unsigned int c20, c21, c22, c23, c24, c25, c26, c27;
    unsigned int c30, c31, c32, c33, c34, c35, c36, c37;

    c00 = c01 = c02 = c03 = c04 = c05 = c06 = c07 = 0;
    c10 = c11 = c12 = c13 = c14 = c15 = c16 = c17 = 0;
    c20 = c21 = c22 = c23 = c24 = c25 = c26 = c27 = 0;
    c30 = c31 = c32 = c33 = c34 = c35 = c36 = c37 = 0;

    unsigned int l;
    for (l = 0; l < k; l++) {
      a0 = *(a + 0*BLOCK_SIZE_X/4);
      a1 = *(a + 1*BLOCK_SIZE_X/4);
      a2 = *(a + 2*BLOCK_SIZE_X/4);
      a3 = *(a + 3*BLOCK_SIZE_X/4);

      b0 = *(b + 0*BLOCK_SIZE_Y/8);
      b1 = *(b + 1*BLOCK_SIZE_Y/8);
      b2 = *(b + 2*BLOCK_SIZE_Y/8);
      b3 = *(b + 3*BLOCK_SIZE_Y/8);
      b4 = *(b + 4*BLOCK_SIZE_Y/8);
      b5 = *(b + 5*BLOCK_SIZE_Y/8);
      b6 = *(b + 6*BLOCK_SIZE_Y/8);
      b7 = *(b + 7*BLOCK_SIZE_Y/8);

      a += mr_alg;
      b += nr_alg;

      t00 = a0 & b0;
      t01 = a0 & b1;
      t02 = a0 & b2;
      t03 = a0 & b3;
      t04 = a0 & b4;
      t05 = a0 & b5;
      t06 = a0 & b6;
      t07 = a0 & b7;

      t10 = a1 & b0;
      t11 = a1 & b1;
      t12 = a1 & b2;
      t13 = a1 & b3;
      t14 = a1 & b4;
      t15 = a1 & b5;
      t16 = a1 & b6;
      t17 = a1 & b7;

      t20 = a2 & b0;
      t21 = a2 & b1;
      t22 = a2 & b2;
      t23 = a2 & b3;
      t24 = a2 & b4;
      t25 = a2 & b5;
      t26 = a2 & b6;
      t27 = a2 & b7;

      t30 = a3 & b0;
      t31 = a3 & b1;
      t32 = a3 & b2;
      t33 = a3 & b3;
      t34 = a3 & b4;
      t35 = a3 & b5;
      t36 = a3 & b6;
      t37 = a3 & b7;

      c00 += popcount(t00);
      c01 += popcount(t01);
      c02 += popcount(t02);
      c03 += popcount(t03);
      c04 += popcount(t04);
      c05 += popcount(t05);
      c06 += popcount(t06);
      c07 += popcount(t07);

      c10 += popcount(t10);
      c11 += popcount(t11);
      c12 += popcount(t12);
      c13 += popcount(t13);
      c14 += popcount(t14);
      c15 += popcount(t15);
      c16 += popcount(t16);
      c17 += popcount(t17);

      c20 += popcount(t20);
      c21 += popcount(t21);
      c22 += popcount(t22);
      c23 += popcount(t23);
      c24 += popcount(t24);
      c25 += popcount(t25);
      c26 += popcount(t26);
      c27 += popcount(t27);

      c30 += popcount(t30);
      c31 += popcount(t31);
      c32 += popcount(t32);
      c33 += popcount(t33);
      c34 += popcount(t34);
      c35 += popcount(t35);
      c36 += popcount(t36);
      c37 += popcount(t37);
    }

  *(c + cs_c*0*BLOCK_SIZE_Y/8 + 0*BLOCK_SIZE_X/4) += c00;
  *(c + cs_c*1*BLOCK_SIZE_Y/8 + 0*BLOCK_SIZE_X/4) += c01;
  *(c + cs_c*2*BLOCK_SIZE_Y/8 + 0*BLOCK_SIZE_X/4) += c02;
  *(c + cs_c*3*BLOCK_SIZE_Y/8 + 0*BLOCK_SIZE_X/4) += c03;
  *(c + cs_c*4*BLOCK_SIZE_Y/8 + 0*BLOCK_SIZE_X/4) += c04;
  *(c + cs_c*5*BLOCK_SIZE_Y/8 + 0*BLOCK_SIZE_X/4) += c05;
  *(c + cs_c*6*BLOCK_SIZE_Y/8 + 0*BLOCK_SIZE_X/4) += c06;
  *(c + cs_c*7*BLOCK_SIZE_Y/8 + 0*BLOCK_SIZE_X/4) += c07;

  *(c + cs_c*0*BLOCK_SIZE_Y/8 + 1*BLOCK_SIZE_X/4) += c10;
  *(c + cs_c*1*BLOCK_SIZE_Y/8 + 1*BLOCK_SIZE_X/4) += c11;
  *(c + cs_c*2*BLOCK_SIZE_Y/8 + 1*BLOCK_SIZE_X/4) += c12;
  *(c + cs_c*3*BLOCK_SIZE_Y/8 + 1*BLOCK_SIZE_X/4) += c13;
  *(c + cs_c*4*BLOCK_SIZE_Y/8 + 1*BLOCK_SIZE_X/4) += c14;
  *(c + cs_c*5*BLOCK_SIZE_Y/8 + 1*BLOCK_SIZE_X/4) += c15;
  *(c + cs_c*6*BLOCK_SIZE_Y/8 + 1*BLOCK_SIZE_X/4) += c16;
  *(c + cs_c*7*BLOCK_SIZE_Y/8 + 1*BLOCK_SIZE_X/4) += c17;

  *(c + cs_c*0*BLOCK_SIZE_Y/8 + 2*BLOCK_SIZE_X/4) += c20;
  *(c + cs_c*1*BLOCK_SIZE_Y/8 + 2*BLOCK_SIZE_X/4) += c21;
  *(c + cs_c*2*BLOCK_SIZE_Y/8 + 2*BLOCK_SIZE_X/4) += c22;
  *(c + cs_c*3*BLOCK_SIZE_Y/8 + 2*BLOCK_SIZE_X/4) += c23;
  *(c + cs_c*4*BLOCK_SIZE_Y/8 + 2*BLOCK_SIZE_X/4) += c24;
  *(c + cs_c*5*BLOCK_SIZE_Y/8 + 2*BLOCK_SIZE_X/4) += c25;
  *(c + cs_c*6*BLOCK_SIZE_Y/8 + 2*BLOCK_SIZE_X/4) += c26;
  *(c + cs_c*7*BLOCK_SIZE_Y/8 + 2*BLOCK_SIZE_X/4) += c27;

  *(c + cs_c*0*BLOCK_SIZE_Y/8 + 3*BLOCK_SIZE_X/4) += c30;
  *(c + cs_c*1*BLOCK_SIZE_Y/8 + 3*BLOCK_SIZE_X/4) += c31;
  *(c + cs_c*2*BLOCK_SIZE_Y/8 + 3*BLOCK_SIZE_X/4) += c32;
  *(c + cs_c*3*BLOCK_SIZE_Y/8 + 3*BLOCK_SIZE_X/4) += c33;
  *(c + cs_c*4*BLOCK_SIZE_Y/8 + 3*BLOCK_SIZE_X/4) += c34;
  *(c + cs_c*5*BLOCK_SIZE_Y/8 + 3*BLOCK_SIZE_X/4) += c35;
  *(c + cs_c*6*BLOCK_SIZE_Y/8 + 3*BLOCK_SIZE_X/4) += c36;
  *(c + cs_c*7*BLOCK_SIZE_Y/8 + 3*BLOCK_SIZE_X/4) += c37;

}
