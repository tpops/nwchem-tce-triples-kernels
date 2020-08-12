/* These have been separated out from ccsd_t_singles_l.F and ccsd_t_doubles_l.F
 * and converted from Fortran 77 to C99 by hand. */

#include "ccsd_t_kernels.h"
typedef struct tensor2d{
  int nr;
  int nc;
  double * vals;
  void (*taco2tensor2d)(taco_tensor_t*,tensor2d*);
  taco_tensor_t* (*tensor2d2taco) (struct tensor2d*);
} tensor2d;


typedef struct tensor4d{
  int ni;
  int nj;
  int nk;
  int nl;
  double * vals;
  void (*taco2tensor4d)(taco_tensor_t*,tensor4d*);
  taco_tensor_t* (*tensor4d2taco) (struct tensor4d*);
} tensor4d;

typedef struct tensor6d{
  int ni;
  int nj;
  int nk;
  int nl;
  int nm;
  int nn;
  double * vals;
  void (*taco2tensor6d)(taco_tensor_t*,tensor6d*);
  taco_tensor_t* (*tensor6d2taco) (struct tensor6d*);
} tensor6d;


[[clang::syntax(taco)]] void c1d_sd_t_s1_1_taco(tensor6d* t3, tensor2d* t2,
	       tensor4d* v){
    t3(h3,h2,h1,p6,p5,p4)= t1(p4,h1) * v2(h3,h2,p6,p5)
}	

/* The possible downcast from fint to int is safe because array dimensions
 * are on the order of 30, i.e. nowhere near INT_MAX. */

void c1d_sd_t_s1_1_(fint * h3d, fint * h2d, fint * h1d,
                    fint * p6d, fint * p5d, fint * p4d,
                    double * restrict t3, const double * restrict t1, const double * restrict v2)
{
    const int h1u = (int)(*h1d);
    const int h2u = (int)(*h2d);
    const int h3u = (int)(*h3d);
    const int p4u = (int)(*p4d);
    const int p5u = (int)(*p5d);
    const int p6u = (int)(*p6d);
    //t3(h3,h2,h1,p6,p5,p4)+=t1(p4,h1)*v2(h3,h2,p6,p5);
    tensor6d t3_tensor;
    t3_tensor.ni = h3u;
    t3_tensor.nj = h2u; 
    t3_tensor.nk = h1u;
    t3_tensor.nl = p6u;
    t3_tensor.nm = p5u;
    t3_tensor.nn = p4u;
    t3_tensor.vals = t3



    
    return;
}

void c1d_sd_t_s1_2_(fint * h3d, fint * h2d, fint * h1d,
                    fint * p6d, fint * p5d, fint * p4d,
                    double * restrict t3, const double * restrict t1, const double * restrict v2)
{
    const int h1u = (int)(*h1d);
    const int h2u = (int)(*h2d);
    const int h3u = (int)(*h3d);
    const int p4u = (int)(*p4d);
    const int p5u = (int)(*p5d);
    const int p6u = (int)(*p6d);
    OMP_PARALLEL_FOR_COLLAPSE(OMP_COLLAPSE_LEVEL)
    for (int p4=0; p4<p4u; p4++)
    for (int p5=0; p5<p5u; p5++)
    for (int p6=0; p6<p6u; p6++)
    for (int h1=0; h1<h1u; h1++)
    PRAGMA_UNROLL
    for (int h2=0; h2<h2u; h2++)
    PRAGMA_UNROLL
    PRAGMA_SIMD
    for (int h3=0; h3<h3u; h3++)
    //t3(h3,h1,h2,p6,p5,p4)-=t1(p4,h1)*v2(h3,h2,p6,p5);
    t3[h3+h3u*(h1+h1u*(h2+h2u*(p6+p6u*(p5+p5u*p4))))] -= t1[p4+p4u*h1] * v2[h3+h3u*(h2+h2u*(p6+p6u*p5))];
    return;
}

void c1d_sd_t_s1_3_(fint * h3d, fint * h2d, fint * h1d,
                    fint * p6d, fint * p5d, fint * p4d,
                    double * restrict t3, const double * restrict t1, const double * restrict v2)
{
    const int h1u = (int)(*h1d);
    const int h2u = (int)(*h2d);
    const int h3u = (int)(*h3d);
    const int p4u = (int)(*p4d);
    const int p5u = (int)(*p5d);
    const int p6u = (int)(*p6d);
    OMP_PARALLEL_FOR_COLLAPSE(OMP_COLLAPSE_LEVEL)
    for (int p4=0; p4<p4u; p4++)
    for (int p5=0; p5<p5u; p5++)
    for (int p6=0; p6<p6u; p6++)
    for (int h2=0; h2<h2u; h2++)
    PRAGMA_UNROLL
    for (int h3=0; h3<h3u; h3++)
    PRAGMA_UNROLL
    PRAGMA_SIMD
    for (int h1=0; h1<h1u; h1++)
    //t3(h1,h3,h2,p6,p5,p4)+=t1(p4,h1)*v2(h3,h2,p6,p5);
    t3[h1+h1u*(h3+h3u*(h2+h2u*(p6+p6u*(p5+p5u*p4))))] += t1[p4+p4u*h1] * v2[h3+h3u*(h2+h2u*(p6+p6u*p5))];
    return;
}

void c1d_sd_t_s1_4_(fint * h3d, fint * h2d, fint * h1d,
                    fint * p6d, fint * p5d, fint * p4d,
                    double * restrict t3, const double * restrict t1, const double * restrict v2)
{
    const int h1u = (int)(*h1d);
    const int h2u = (int)(*h2d);
    const int h3u = (int)(*h3d);
    const int p4u = (int)(*p4d);
    const int p5u = (int)(*p5d);
    const int p6u = (int)(*p6d);
    OMP_PARALLEL_FOR_COLLAPSE(OMP_COLLAPSE_LEVEL)
    for (int p5=0; p5<p5u; p5++)
    for (int p4=0; p4<p4u; p4++)
    for (int p6=0; p6<p6u; p6++)
    for (int h1=0; h1<h1u; h1++)
    for (int h2=0; h2<h2u; h2++)
    PRAGMA_SIMD
    for (int h3=0; h3<h3u; h3++)
    //t3(h3,h2,h1,p6,p4,p5)-=t1(p4,h1)*v2(h3,h2,p6,p5);
    t3[h3+h3u*(h2+h2u*(h1+h1u*(p6+p6u*(p4+p4u*p5))))] -= t1[p4+p4u*h1] * v2[h3+h3u*(h2+h2u*(p6+p6u*p5))];
    return;
}

void c1d_sd_t_s1_5_(fint * h3d, fint * h2d, fint * h1d,
                    fint * p6d, fint * p5d, fint * p4d,
                    double * restrict t3, const double * restrict t1, const double * restrict v2)
{
    const int h1u = (int)(*h1d);
    const int h2u = (int)(*h2d);
    const int h3u = (int)(*h3d);
    const int p4u = (int)(*p4d);
    const int p5u = (int)(*p5d);
    const int p6u = (int)(*p6d);
    OMP_PARALLEL_FOR_COLLAPSE(OMP_COLLAPSE_LEVEL)
    for (int p5=0; p5<p5u; p5++)
    for (int p4=0; p4<p4u; p4++)
    for (int p6=0; p6<p6u; p6++)
    for (int h1=0; h1<h1u; h1++)
    PRAGMA_UNROLL
    for (int h2=0; h2<h2u; h2++)
    PRAGMA_UNROLL
    PRAGMA_SIMD
    for (int h3=0; h3<h3u; h3++)
    //t3(h3,h1,h2,p6,p4,p5)+=t1(p4,h1)*v2(h3,h2,p6,p5);
    t3[h3+h3u*(h1+h1u*(h2+h2u*(p6+p6u*(p4+p4u*p5))))] += t1[p4+p4u*h1] * v2[h3+h3u*(h2+h2u*(p6+p6u*p5))];
    return;
}

void c1d_sd_t_s1_6_(fint * h3d, fint * h2d, fint * h1d,
                    fint * p6d, fint * p5d, fint * p4d,
                    double * restrict t3, const double * restrict t1, const double * restrict v2)
{
    const int h1u = (int)(*h1d);
    const int h2u = (int)(*h2d);
    const int h3u = (int)(*h3d);
    const int p4u = (int)(*p4d);
    const int p5u = (int)(*p5d);
    const int p6u = (int)(*p6d);
    OMP_PARALLEL_FOR_COLLAPSE(OMP_COLLAPSE_LEVEL)
    for (int p5=0; p5<p5u; p5++)
    for (int p4=0; p4<p4u; p4++)
    for (int p6=0; p6<p6u; p6++)
    for (int h2=0; h2<h2u; h2++)
    PRAGMA_UNROLL
    for (int h3=0; h3<h3u; h3++)
    PRAGMA_UNROLL
    PRAGMA_SIMD
    for (int h1=0; h1<h1u; h1++)
    //t3(h1,h3,h2,p6,p4,p5)-=t1(p4,h1)*v2(h3,h2,p6,p5);
    t3[h1+h1u*(h3+h3u*(h2+h2u*(p6+p6u*(p4+p4u*p5))))] -= t1[p4+p4u*h1] * v2[h3+h3u*(h2+h2u*(p6+p6u*p5))];
    return;
}

void c1d_sd_t_s1_7_(fint * h3d, fint * h2d, fint * h1d,
                    fint * p6d, fint * p5d, fint * p4d,
                    double * restrict t3, const double * restrict t1, const double * restrict v2)
{
    const int h1u = (int)(*h1d);
    const int h2u = (int)(*h2d);
    const int h3u = (int)(*h3d);
    const int p4u = (int)(*p4d);
    const int p5u = (int)(*p5d);
    const int p6u = (int)(*p6d);
    OMP_PARALLEL_FOR_COLLAPSE(OMP_COLLAPSE_LEVEL)
    for (int p5=0; p5<p5u; p5++)
    for (int p6=0; p6<p6u; p6++)
    for (int p4=0; p4<p4u; p4++)
    for (int h1=0; h1<h1u; h1++)
    PRAGMA_UNROLL
    for (int h2=0; h2<h2u; h2++)
    PRAGMA_UNROLL
    PRAGMA_SIMD
    for (int h3=0; h3<h3u; h3++)
    //t3(h3,h2,h1,p4,p6,p5)+=t1(p4,h1)*v2(h3,h2,p6,p5);
    t3[h3+h3u*(h2+h2u*(h1+h1u*(p4+p4u*(p6+p6u*p5))))] += t1[p4+p4u*h1] * v2[h3+h3u*(h2+h2u*(p6+p6u*p5))];
    return;
}

void c1d_sd_t_s1_8_(fint * h3d, fint * h2d, fint * h1d,
                    fint * p6d, fint * p5d, fint * p4d,
                    double * restrict t3, const double * restrict t1, const double * restrict v2)
{
    const int h1u = (int)(*h1d);
    const int h2u = (int)(*h2d);
    const int h3u = (int)(*h3d);
    const int p4u = (int)(*p4d);
    const int p5u = (int)(*p5d);
    const int p6u = (int)(*p6d);
    OMP_PARALLEL_FOR_COLLAPSE(OMP_COLLAPSE_LEVEL)
    for (int p5=0; p5<p5u; p5++)
    for (int p6=0; p6<p6u; p6++)
    for (int p4=0; p4<p4u; p4++)
    for (int h2=0; h2<h2u; h2++)
    PRAGMA_UNROLL
    for (int h1=0; h1<h1u; h1++)
    PRAGMA_UNROLL
    PRAGMA_SIMD
    for (int h3=0; h3<h3u; h3++)
    //t3(h3,h1,h2,p4,p6,p5)-=t1(p4,h1)*v2(h3,h2,p6,p5);
    t3[h3+h3u*(h1+h1u*(h2+h2u*(p4+p4u*(p6+p6u*p5))))] -= t1[p4+p4u*h1] * v2[h3+h3u*(h2+h2u*(p6+p6u*p5))];
    return;
}

void c1d_sd_t_s1_9_(fint * h3d, fint * h2d, fint * h1d,
                    fint * p6d, fint * p5d, fint * p4d,
                    double * restrict t3, const double * restrict t1, const double * restrict v2)
{
    const int h1u = (int)(*h1d);
    const int h2u = (int)(*h2d);
    const int h3u = (int)(*h3d);
    const int p4u = (int)(*p4d);
    const int p5u = (int)(*p5d);
    const int p6u = (int)(*p6d);
    OMP_PARALLEL_FOR_COLLAPSE(OMP_COLLAPSE_LEVEL)
    for (int p5=0; p5<p5u; p5++)
    for (int p6=0; p6<p6u; p6++)
    for (int p4=0; p4<p4u; p4++)
    for (int h2=0; h2<h2u; h2++)
    PRAGMA_UNROLL
    for (int h3=0; h3<h3u; h3++)
    PRAGMA_UNROLL
    PRAGMA_SIMD
    for (int h1=0; h1<h1u; h1++)
    //t3(h1,h3,h2,p4,p6,p5)+=t1(p4,h1)*v2(h3,h2,p6,p5);
    t3[h1+h1u*(h3+h3u*(h2+h2u*(p4+p4u*(p6+p6u*p5))))] += t1[p4+p4u*h1] * v2[h3+h3u*(h2+h2u*(p6+p6u*p5))];
    return;
}

void c1d_sd_t_d1_1_(fint * h3d, fint * h2d, fint * h1d,
                    fint * p6d, fint * p5d, fint * p4d,
                    fint * h7d,
                    double * restrict t3, const double * restrict t2, const double * restrict v2)
{
    const int h1u = (int)(*h1d);
    const int h2u = (int)(*h2d);
    const int h3u = (int)(*h3d);
    const int p4u = (int)(*p4d);
    const int p5u = (int)(*p5d);
    const int p6u = (int)(*p6d);
    const int h7u = (int)(*h7d);
    OMP_PARALLEL_FOR_COLLAPSE(OMP_COLLAPSE_LEVEL)
    for (int p4=0; p4<p4u; p4++)
    for (int p5=0; p5<p5u; p5++)
    for (int p6=0; p6<p6u; p6++)
    for (int h1=0; h1<h1u; h1++)
    for (int h2=0; h2<h2u; h2++)
    for (int h3=0; h3<h3u; h3++)
    PRAGMA_SIMD
    for (int h7=0; h7<h7u; h7++)
    //t3(h3,h2,h1,p6,p5,p4)-=t2(h7,p4,p5,h1)*v2(h3,h2,p6,h7);
    t3[h3+h3u*(h2+h2u*(h1+h1u*(p6+p6u*(p5+p5u*p4))))] -= t2[h7+h7u*(p4+p4u*(p5+p5u*h1))] * v2[h3+h3u*(h2+h2u*(p6+p6u*h7))];
    return;
}

void c1d_sd_t_d1_2_(fint * h3d, fint * h2d, fint * h1d,
                    fint * p6d, fint * p5d, fint * p4d,
                    fint * h7d,
                    double * restrict t3, const double * restrict t2, const double * restrict v2)
{
    const int h1u = (int)(*h1d);
    const int h2u = (int)(*h2d);
    const int h3u = (int)(*h3d);
    const int p4u = (int)(*p4d);
    const int p5u = (int)(*p5d);
    const int p6u = (int)(*p6d);
    const int h7u = (int)(*h7d);
    OMP_PARALLEL_FOR_COLLAPSE(OMP_COLLAPSE_LEVEL)
    for (int p4=0; p4<p4u; p4++)
    for (int p5=0; p5<p5u; p5++)
    for (int p6=0; p6<p6u; p6++)
    for (int h2=0; h2<h2u; h2++)
    for (int h1=0; h1<h1u; h1++)
    for (int h3=0; h3<h3u; h3++)
    PRAGMA_SIMD
    for (int h7=0; h7<h7u; h7++)
    //t3(h3,h1,h2,p6,p5,p4)+=t2(h7,p4,p5,h1)*v2(h3,h2,p6,h7);
    t3[h3+h3u*(h1+h1u*(h2+h2u*(p6+p6u*(p5+p5u*p4))))] += t2[h7+h7u*(p4+p4u*(p5+p5u*h1))] * v2[h3+h3u*(h2+h2u*(p6+p6u*h7))];
    return;
}

void c1d_sd_t_d1_3_(fint * h3d, fint * h2d, fint * h1d,
                    fint * p6d, fint * p5d, fint * p4d,
                    fint * h7d,
                    double * restrict t3, const double * restrict t2, const double * restrict v2)
{
    const int h1u = (int)(*h1d);
    const int h2u = (int)(*h2d);
    const int h3u = (int)(*h3d);
    const int p4u = (int)(*p4d);
    const int p5u = (int)(*p5d);
    const int p6u = (int)(*p6d);
    const int h7u = (int)(*h7d);
    OMP_PARALLEL_FOR_COLLAPSE(OMP_COLLAPSE_LEVEL)
    for (int p4=0; p4<p4u; p4++)
    for (int p5=0; p5<p5u; p5++)
    for (int p6=0; p6<p6u; p6++)
    for (int h2=0; h2<h2u; h2++)
    for (int h3=0; h3<h3u; h3++)
    for (int h1=0; h1<h1u; h1++)
    PRAGMA_SIMD
    for (int h7=0; h7<h7u; h7++)
    //t3(h1,h3,h2,p6,p5,p4)-=t2(h7,p4,p5,h1)*v2(h3,h2,p6,h7);
    t3[h1+h1u*(h3+h3u*(h2+h2u*(p6+p6u*(p5+p5u*p4))))] -= t2[h7+h7u*(p4+p4u*(p5+p5u*h1))] * v2[h3+h3u*(h2+h2u*(p6+p6u*h7))];
    return;
}

void c1d_sd_t_d1_4_(fint * h3d, fint * h2d, fint * h1d,
                    fint * p6d, fint * p5d, fint * p4d,
                    fint * h7d,
                    double * restrict t3, const double * restrict t2, const double * restrict v2)
{
    const int h1u = (int)(*h1d);
    const int h2u = (int)(*h2d);
    const int h3u = (int)(*h3d);
    const int p4u = (int)(*p4d);
    const int p5u = (int)(*p5d);
    const int p6u = (int)(*p6d);
    const int h7u = (int)(*h7d);
    OMP_PARALLEL_FOR_COLLAPSE(OMP_COLLAPSE_LEVEL)
    for (int p6=0; p6<p6u; p6++)
    for (int p4=0; p4<p4u; p4++)
    for (int p5=0; p5<p5u; p5++)
    for (int h1=0; h1<h1u; h1++)
    for (int h2=0; h2<h2u; h2++)
    for (int h3=0; h3<h3u; h3++)
    PRAGMA_SIMD
    for (int h7=0; h7<h7u; h7++)
    //t3(h3,h2,h1,p5,p4,p6)-=t2(h7,p4,p5,h1)*v2(h3,h2,p6,h7);
    t3[h3+h3u*(h2+h2u*(h1+h1u*(p5+p5u*(p4+p4u*p6))))] -= t2[h7+h7u*(p4+p4u*(p5+p5u*h1))] * v2[h3+h3u*(h2+h2u*(p6+p6u*h7))];
    return;
}

void c1d_sd_t_d1_5_(fint * h3d, fint * h2d, fint * h1d,
                    fint * p6d, fint * p5d, fint * p4d,
                    fint * h7d,
                    double * restrict t3, const double * restrict t2, const double * restrict v2)
{
    const int h1u = (int)(*h1d);
    const int h2u = (int)(*h2d);
    const int h3u = (int)(*h3d);
    const int p4u = (int)(*p4d);
    const int p5u = (int)(*p5d);
    const int p6u = (int)(*p6d);
    const int h7u = (int)(*h7d);
    OMP_PARALLEL_FOR_COLLAPSE(OMP_COLLAPSE_LEVEL)
    for (int p6=0; p6<p6u; p6++)
    for (int p4=0; p4<p4u; p4++)
    for (int p5=0; p5<p5u; p5++)
    for (int h2=0; h2<h2u; h2++)
    for (int h1=0; h1<h1u; h1++)
    for (int h3=0; h3<h3u; h3++)
    PRAGMA_SIMD
    for (int h7=0; h7<h7u; h7++)
    //t3(h3,h1,h2,p5,p4,p6)+=t2(h7,p4,p5,h1)*v2(h3,h2,p6,h7);
    t3[h3+h3u*(h1+h1u*(h2+h2u*(p5+p5u*(p4+p4u*p6))))] += t2[h7+h7u*(p4+p4u*(p5+p5u*h1))] * v2[h3+h3u*(h2+h2u*(p6+p6u*h7))];
    return;
}

void c1d_sd_t_d1_6_(fint * h3d, fint * h2d, fint * h1d,
                    fint * p6d, fint * p5d, fint * p4d,
                    fint * h7d,
                    double * restrict t3, const double * restrict t2, const double * restrict v2)
{
    const int h1u = (int)(*h1d);
    const int h2u = (int)(*h2d);
    const int h3u = (int)(*h3d);
    const int p4u = (int)(*p4d);
    const int p5u = (int)(*p5d);
    const int p6u = (int)(*p6d);
    const int h7u = (int)(*h7d);
    OMP_PARALLEL_FOR_COLLAPSE(OMP_COLLAPSE_LEVEL)
    for (int p6=0; p6<p6u; p6++)
    for (int p4=0; p4<p4u; p4++)
    for (int p5=0; p5<p5u; p5++)
    for (int h2=0; h2<h2u; h2++)
    for (int h3=0; h3<h3u; h3++)
    for (int h1=0; h1<h1u; h1++)
    PRAGMA_SIMD
    for (int h7=0; h7<h7u; h7++)
    //t3(h1,h3,h2,p5,p4,p6)-=t2(h7,p4,p5,h1)*v2(h3,h2,p6,h7);
    t3[h1+h1u*(h3+h3u*(h2+h2u*(p5+p5u*(p4+p4u*p6))))] -= t2[h7+h7u*(p4+p4u*(p5+p5u*h1))] * v2[h3+h3u*(h2+h2u*(p6+p6u*h7))];
    return;
}

void c1d_sd_t_d1_7_(fint * h3d, fint * h2d, fint * h1d,
                    fint * p6d, fint * p5d, fint * p4d,
                    fint * h7d,
                    double * restrict t3, const double * restrict t2, const double * restrict v2)
{
    const int h1u = (int)(*h1d);
    const int h2u = (int)(*h2d);
    const int h3u = (int)(*h3d);
    const int p4u = (int)(*p4d);
    const int p5u = (int)(*p5d);
    const int p6u = (int)(*p6d);
    const int h7u = (int)(*h7d);
    OMP_PARALLEL_FOR_COLLAPSE(OMP_COLLAPSE_LEVEL)
    for (int p4=0; p4<p4u; p4++)
    for (int p6=0; p6<p6u; p6++)
    for (int p5=0; p5<p5u; p5++)
    for (int h1=0; h1<h1u; h1++)
    for (int h2=0; h2<h2u; h2++)
    for (int h3=0; h3<h3u; h3++)
    PRAGMA_SIMD
    for (int h7=0; h7<h7u; h7++)
    //t3(h3,h2,h1,p5,p6,p4)+=t2(h7,p4,p5,h1)*v2(h3,h2,p6,h7);
    t3[h3+h3u*(h2+h2u*(h1+h1u*(p5+p5u*(p6+p6u*p4))))] += t2[h7+h7u*(p4+p4u*(p5+p5u*h1))] * v2[h3+h3u*(h2+h2u*(p6+p6u*h7))];
    return;
}

void c1d_sd_t_d1_8_(fint * h3d, fint * h2d, fint * h1d,
                    fint * p6d, fint * p5d, fint * p4d,
                    fint * h7d,
                    double * restrict t3, const double * restrict t2, const double * restrict v2)
{
    const int h1u = (int)(*h1d);
    const int h2u = (int)(*h2d);
    const int h3u = (int)(*h3d);
    const int p4u = (int)(*p4d);
    const int p5u = (int)(*p5d);
    const int p6u = (int)(*p6d);
    const int h7u = (int)(*h7d);
    OMP_PARALLEL_FOR_COLLAPSE(OMP_COLLAPSE_LEVEL)
    for (int p4=0; p4<p4u; p4++)
    for (int p6=0; p6<p6u; p6++)
    for (int p5=0; p5<p5u; p5++)
    for (int h2=0; h2<h2u; h2++)
    for (int h1=0; h1<h1u; h1++)
    for (int h3=0; h3<h3u; h3++)
    PRAGMA_SIMD
    for (int h7=0; h7<h7u; h7++)
    //t3(h3,h1,h2,p5,p6,p4)-=t2(h7,p4,p5,h1)*v2(h3,h2,p6,h7);
    t3[h3+h3u*(h1+h1u*(h2+h2u*(p5+p5u*(p6+p6u*p4))))] -= t2[h7+h7u*(p4+p4u*(p5+p5u*h1))] * v2[h3+h3u*(h2+h2u*(p6+p6u*h7))];
    return;
}

void c1d_sd_t_d1_9_(fint * h3d, fint * h2d, fint * h1d,
                    fint * p6d, fint * p5d, fint * p4d,
                    fint * h7d,
                    double * restrict t3, const double * restrict t2, const double * restrict v2)
{
    const int h1u = (int)(*h1d);
    const int h2u = (int)(*h2d);
    const int h3u = (int)(*h3d);
    const int p4u = (int)(*p4d);
    const int p5u = (int)(*p5d);
    const int p6u = (int)(*p6d);
    const int h7u = (int)(*h7d);
    OMP_PARALLEL_FOR_COLLAPSE(OMP_COLLAPSE_LEVEL)
    for (int p4=0; p4<p4u; p4++)
    for (int p6=0; p6<p6u; p6++)
    for (int p5=0; p5<p5u; p5++)
    for (int h2=0; h2<h2u; h2++)
    for (int h3=0; h3<h3u; h3++)
    for (int h1=0; h1<h1u; h1++)
    PRAGMA_SIMD
    for (int h7=0; h7<h7u; h7++)
    //t3(h1,h3,h2,p5,p6,p4)+=t2(h7,p4,p5,h1)*v2(h3,h2,p6,h7);
    t3[h1+h1u*(h3+h3u*(h2+h2u*(p5+p5u*(p6+p6u*p4))))] += t2[h7+h7u*(p4+p4u*(p5+p5u*h1))] * v2[h3+h3u*(h2+h2u*(p6+p6u*h7))];
    return;
}

void c1d_sd_t_d2_1_(fint * h3d, fint * h2d, fint * h1d,
                    fint * p6d, fint * p5d, fint * p4d,
                    fint * h7d,
                    double * restrict t3, const double * restrict t2, const double * restrict v2)
{
    const int h1u = (int)(*h1d);
    const int h2u = (int)(*h2d);
    const int h3u = (int)(*h3d);
    const int p4u = (int)(*p4d);
    const int p5u = (int)(*p5d);
    const int p6u = (int)(*p6d);
    const int h7u = (int)(*h7d);
    OMP_PARALLEL_FOR_COLLAPSE(OMP_COLLAPSE_LEVEL)
    for (int p4=0; p4<p4u; p4++)
    for (int p5=0; p5<p5u; p5++)
    for (int p6=0; p6<p6u; p6++)
    for (int h1=0; h1<h1u; h1++)
    for (int h2=0; h2<h2u; h2++)
    for (int h3=0; h3<h3u; h3++)
    PRAGMA_SIMD
    for (int h7=0; h7<h7u; h7++)
    //t3(h3,h2,h1,p6,p5,p4)-=t2(h7,p4,h1,h2)*v2(h7,h3,p6,p5);
    t3[h3+h3u*(h2+h2u*(h1+h1u*(p6+p6u*(p5+p5u*p4))))] -= t2[h7+h7u*(p4+p4u*(h1+h1u*h2))] * v2[h7+h7u*(h3+h3u*(p6+p6u*p5))];
    return;
}

void c1d_sd_t_d2_2_(fint * h3d, fint * h2d, fint * h1d,
                    fint * p6d, fint * p5d, fint * p4d,
                    fint * h7d,
                    double * restrict t3, const double * restrict t2, const double * restrict v2)
{
    const int h1u = (int)(*h1d);
    const int h2u = (int)(*h2d);
    const int h3u = (int)(*h3d);
    const int p4u = (int)(*p4d);
    const int p5u = (int)(*p5d);
    const int p6u = (int)(*p6d);
    const int h7u = (int)(*h7d);
    OMP_PARALLEL_FOR_COLLAPSE(OMP_COLLAPSE_LEVEL)
    for (int p4=0; p4<p4u; p4++)
    for (int p5=0; p5<p5u; p5++)
    for (int p6=0; p6<p6u; p6++)
    for (int h3=0; h3<h3u; h3++)
    for (int h1=0; h1<h1u; h1++)
    for (int h2=0; h2<h2u; h2++)
    PRAGMA_SIMD
    for (int h7=0; h7<h7u; h7++)
    //t3(h2,h1,h3,p6,p5,p4)-=t2(h7,p4,h1,h2)*v2(h7,h3,p6,p5);
    t3[h2+h2u*(h1+h1u*(h3+h3u*(p6+p6u*(p5+p5u*p4))))] -= t2[h7+h7u*(p4+p4u*(h1+h1u*h2))] * v2[h7+h7u*(h3+h3u*(p6+p6u*p5))];
    return;
}

void c1d_sd_t_d2_3_(fint * h3d, fint * h2d, fint * h1d,
                    fint * p6d, fint * p5d, fint * p4d,
                    fint * h7d,
                    double * restrict t3, const double * restrict t2, const double * restrict v2)
{
    const int h1u = (int)(*h1d);
    const int h2u = (int)(*h2d);
    const int h3u = (int)(*h3d);
    const int p4u = (int)(*p4d);
    const int p5u = (int)(*p5d);
    const int p6u = (int)(*p6d);
    const int h7u = (int)(*h7d);
    OMP_PARALLEL_FOR_COLLAPSE(OMP_COLLAPSE_LEVEL)
    for (int p4=0; p4<p4u; p4++)
    for (int p5=0; p5<p5u; p5++)
    for (int p6=0; p6<p6u; p6++)
    for (int h1=0; h1<h1u; h1++)
    for (int h3=0; h3<h3u; h3++)
    for (int h2=0; h2<h2u; h2++)
    PRAGMA_SIMD
    for (int h7=0; h7<h7u; h7++)
    //t3(h2,h3,h1,p6,p5,p4)+=t2(h7,p4,h1,h2)*v2(h7,h3,p6,p5);
    t3[h2+h2u*(h3+h3u*(h1+h1u*(p6+p6u*(p5+p5u*p4))))] += t2[h7+h7u*(p4+p4u*(h1+h1u*h2))] * v2[h7+h7u*(h3+h3u*(p6+p6u*p5))];
    return;
}

void c1d_sd_t_d2_4_(fint * h3d, fint * h2d, fint * h1d,
                    fint * p6d, fint * p5d, fint * p4d,
                    fint * h7d,
                    double * restrict t3, const double * restrict t2, const double * restrict v2)
{
    const int h1u = (int)(*h1d);
    const int h2u = (int)(*h2d);
    const int h3u = (int)(*h3d);
    const int p4u = (int)(*p4d);
    const int p5u = (int)(*p5d);
    const int p6u = (int)(*p6d);
    const int h7u = (int)(*h7d);
    OMP_PARALLEL_FOR_COLLAPSE(OMP_COLLAPSE_LEVEL)
    for (int p5=0; p5<p5u; p5++)
    for (int p4=0; p4<p4u; p4++)
    for (int p6=0; p6<p6u; p6++)
    for (int h1=0; h1<h1u; h1++)
    for (int h2=0; h2<h2u; h2++)
    for (int h3=0; h3<h3u; h3++)
    PRAGMA_SIMD
    for (int h7=0; h7<h7u; h7++)
    //t3(h3,h2,h1,p6,p4,p5)+=t2(h7,p4,h1,h2)*v2(h7,h3,p6,p5);
    t3[h3+h3u*(h2+h2u*(h1+h1u*(p6+p6u*(p4+p4u*p5))))] += t2[h7+h7u*(p4+p4u*(h1+h1u*h2))] * v2[h7+h7u*(h3+h3u*(p6+p6u*p5))];
    return;
}

void c1d_sd_t_d2_5_(fint * h3d, fint * h2d, fint * h1d,
                    fint * p6d, fint * p5d, fint * p4d,
                    fint * h7d,
                    double * restrict t3, const double * restrict t2, const double * restrict v2)
{
    const int h1u = (int)(*h1d);
    const int h2u = (int)(*h2d);
    const int h3u = (int)(*h3d);
    const int p4u = (int)(*p4d);
    const int p5u = (int)(*p5d);
    const int p6u = (int)(*p6d);
    const int h7u = (int)(*h7d);
    OMP_PARALLEL_FOR_COLLAPSE(OMP_COLLAPSE_LEVEL)
    for (int p5=0; p5<p5u; p5++)
    for (int p4=0; p4<p4u; p4++)
    for (int p6=0; p6<p6u; p6++)
    for (int h3=0; h3<h3u; h3++)
    for (int h1=0; h1<h1u; h1++)
    for (int h2=0; h2<h2u; h2++)
    PRAGMA_SIMD
    for (int h7=0; h7<h7u; h7++)
    //t3(h2,h1,h3,p6,p4,p5)+=t2(h7,p4,h1,h2)*v2(h7,h3,p6,p5);
    t3[h2+h2u*(h1+h1u*(h3+h3u*(p6+p6u*(p4+p4u*p5))))] += t2[h7+h7u*(p4+p4u*(h1+h1u*h2))] * v2[h7+h7u*(h3+h3u*(p6+p6u*p5))];
    return;
}

void c1d_sd_t_d2_6_(fint * h3d, fint * h2d, fint * h1d,
                    fint * p6d, fint * p5d, fint * p4d,
                    fint * h7d,
                    double * restrict t3, const double * restrict t2, const double * restrict v2)
{
    const int h1u = (int)(*h1d);
    const int h2u = (int)(*h2d);
    const int h3u = (int)(*h3d);
    const int p4u = (int)(*p4d);
    const int p5u = (int)(*p5d);
    const int p6u = (int)(*p6d);
    const int h7u = (int)(*h7d);
    OMP_PARALLEL_FOR_COLLAPSE(OMP_COLLAPSE_LEVEL)
    for (int p5=0; p5<p5u; p5++)
    for (int p4=0; p4<p4u; p4++)
    for (int p6=0; p6<p6u; p6++)
    for (int h1=0; h1<h1u; h1++)
    for (int h3=0; h3<h3u; h3++)
    for (int h2=0; h2<h2u; h2++)
    PRAGMA_SIMD
    for (int h7=0; h7<h7u; h7++)
    //t3(h2,h3,h1,p6,p4,p5)-=t2(h7,p4,h1,h2)*v2(h7,h3,p6,p5);
    t3[h2+h2u*(h3+h3u*(h1+h1u*(p6+p6u*(p4+p4u*p5))))] -= t2[h7+h7u*(p4+p4u*(h1+h1u*h2))] * v2[h7+h7u*(h3+h3u*(p6+p6u*p5))];
    return;
}

void c1d_sd_t_d2_7_(fint * h3d, fint * h2d, fint * h1d,
                    fint * p6d, fint * p5d, fint * p4d,
                    fint * h7d,
                    double * restrict t3, const double * restrict t2, const double * restrict v2)
{
    const int h1u = (int)(*h1d);
    const int h2u = (int)(*h2d);
    const int h3u = (int)(*h3d);
    const int p4u = (int)(*p4d);
    const int p5u = (int)(*p5d);
    const int p6u = (int)(*p6d);
    const int h7u = (int)(*h7d);
    OMP_PARALLEL_FOR_COLLAPSE(OMP_COLLAPSE_LEVEL)
    for (int p5=0; p5<p5u; p5++)
    for (int p6=0; p6<p6u; p6++)
    for (int p4=0; p4<p4u; p4++)
    for (int h1=0; h1<h1u; h1++)
    for (int h2=0; h2<h2u; h2++)
    for (int h3=0; h3<h3u; h3++)
    PRAGMA_SIMD
    for (int h7=0; h7<h7u; h7++)
    //t3(h3,h2,h1,p4,p6,p5)-=t2(h7,p4,h1,h2)*v2(h7,h3,p6,p5);
    t3[h3+h3u*(h2+h2u*(h1+h1u*(p4+p4u*(p6+p6u*p5))))] -= t2[h7+h7u*(p4+p4u*(h1+h1u*h2))] * v2[h7+h7u*(h3+h3u*(p6+p6u*p5))];
    return;
}

void c1d_sd_t_d2_8_(fint * h3d, fint * h2d, fint * h1d,
                    fint * p6d, fint * p5d, fint * p4d,
                    fint * h7d,
                    double * restrict t3, const double * restrict t2, const double * restrict v2)
{
    const int h1u = (int)(*h1d);
    const int h2u = (int)(*h2d);
    const int h3u = (int)(*h3d);
    const int p4u = (int)(*p4d);
    const int p5u = (int)(*p5d);
    const int p6u = (int)(*p6d);
    const int h7u = (int)(*h7d);
    OMP_PARALLEL_FOR_COLLAPSE(OMP_COLLAPSE_LEVEL)
    for (int p5=0; p5<p5u; p5++)
    for (int p6=0; p6<p6u; p6++)
    for (int p4=0; p4<p4u; p4++)
    for (int h3=0; h3<h3u; h3++)
    for (int h1=0; h1<h1u; h1++)
    for (int h2=0; h2<h2u; h2++)
    PRAGMA_SIMD
    for (int h7=0; h7<h7u; h7++)
    //t3(h2,h1,h3,p4,p6,p5)-=t2(h7,p4,h1,h2)*v2(h7,h3,p6,p5);
    t3[h2+h2u*(h1+h1u*(h3+h3u*(p4+p4u*(p6+p6u*p5))))] -= t2[h7+h7u*(p4+p4u*(h1+h1u*h2))] * v2[h7+h7u*(h3+h3u*(p6+p6u*p5))];
    return;
}

void c1d_sd_t_d2_9_(fint * h3d, fint * h2d, fint * h1d,
                    fint * p6d, fint * p5d, fint * p4d,
                    fint * h7d,
                    double * restrict t3, const double * restrict t2, const double * restrict v2)
{
    const int h1u = (int)(*h1d);
    const int h2u = (int)(*h2d);
    const int h3u = (int)(*h3d);
    const int p4u = (int)(*p4d);
    const int p5u = (int)(*p5d);
    const int p6u = (int)(*p6d);
    const int h7u = (int)(*h7d);
    OMP_PARALLEL_FOR_COLLAPSE(OMP_COLLAPSE_LEVEL)
    for (int p5=0; p5<p5u; p5++)
    for (int p6=0; p6<p6u; p6++)
    for (int p4=0; p4<p4u; p4++)
    for (int h1=0; h1<h1u; h1++)
    for (int h3=0; h3<h3u; h3++)
    for (int h2=0; h2<h2u; h2++)
    PRAGMA_SIMD
    for (int h7=0; h7<h7u; h7++)
    //t3(h2,h3,h1,p4,p6,p5)+=t2(h7,p4,h1,h2)*v2(h7,h3,p6,p5);
    t3[h2+h2u*(h3+h3u*(h1+h1u*(p4+p4u*(p6+p6u*p5))))] += t2[h7+h7u*(p4+p4u*(h1+h1u*h2))] * v2[h7+h7u*(h3+h3u*(p6+p6u*p5))];
    return;
}
