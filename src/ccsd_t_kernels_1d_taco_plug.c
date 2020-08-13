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


typedef struct scalar{
  int val;
  void (*taco2scalar)(taco_tensor_t*,struct scalar*);
  taco_tensor_t* (*scalar2taco) (struct scalar*);
} scalar;


[[clang::syntax(taco)]] void c1d_sd_t_s1_1_taco(tensor6d* t3, tensor2d* t2,
	       tensor4d* v2,const char * format=""){
    t3(h3,h2,h1,p6,p5,p4)= t1(p4,h1) * v2(h3,h2,p6,p5)
}	

[[clang::syntax(taco)]] void c1d_sd_t_s1_2_taco(tensor6d* t3, tensor2d* t2,
	       tensor4d* v2, scalar * alpha,const char * format=""){
    t3(h3,h1,h2,p6,p5,p4) = alpha * t1(p4,h1) * v2(h3,h2,p6,p5)
}	

[[clang::syntax(taco)]] void c1d_sd_t_s1_3_taco(tensor6d* t3, tensor2d* t2,
	       tensor4d* v2,const char * format=""){

    t3(h1,h3,h2,p6,p5,p4) = t1(p4,h1) * v2(h3,h2,p6,p5)
}	


[[clang::syntax(taco)]] void c1d_sd_t_s1_4_taco(tensor6d* t3, tensor2d* t2,
	       tensor4d* v2, scalar * alpha,const char * format=""){
  t3(h3,h2,h1,p6,p4,p5)= alpha *  t1(p4,h1)*v2(h3,h2,p6,p5)
}

[[clang::syntax(taco)]] void c1d_sd_t_s1_5_taco(tensor6d* t3, tensor2d* t2,
	       tensor4d* v2,const char * format=""){
  t3(h3,h1,h2,p6,p4,p5)= t1(p4,h1) * v2(h3,h2,p6,p5)
}

[[clang::syntax(taco)]] void c1d_sd_t_s1_6_taco(tensor6d* t3, tensor2d* t2,
	       tensor4d* v2, scalar * alpha,const char * format=""){

  t3(h1,h3,h2,p6,p4,p5) = alpha * t1(p4,h1) * v2(h3,h2,p6,p5)
}


[[clang::syntax(taco)]] void c1d_sd_t_s1_7_taco(tensor6d* t3, tensor2d* t2,
	       tensor4d* v2,const char * format=""){

  t3(h3,h2,h1,p4,p6,p5) = t1(p4,h1) * v2(h3,h2,p6,p5)
}
[[clang::syntax(taco)]] void c1d_sd_t_s1_8_taco(tensor6d* t3, tensor2d* t2,
	       tensor4d* v2, scalar * alpha,const char * format=""){

  t3(h3,h1,h2,p4,p6,p5)= alpha * t1(p4,h1) * v2(h3,h2,p6,p5)
}

[[clang::syntax(taco)]] void c1d_sd_t_s1_9_taco(tensor6d* t3, tensor2d* t2,
	       tensor4d* v2,const char * format=""){
  t3(h1,h3,h2,p4,p6,p5) = t1(p4,h1) * v2(h3,h2,p6,p5)
}


/*Conversion routines.*/
static inline taco_tensor_t *tensor2d2taco(tensor2d *t2d) {
  taco_tensor_t *tensor = (taco_tensor_t *)malloc(sizeof(taco_tensor_t));
  tensor->order = 2;
  tensor->dimensions = (int32_t *)malloc(sizeof(int32_t) * 2);
  tensor->dimensions[0] = t2d->nr;
  tensor->dimensions[1] = t2d->nc;
  tensor->csize = sizeof(*t2d->vals);
  tensor->mode_ordering = (int32_t *)malloc(sizeof(int32_t) * 2);
  // storage data layout order
  if (t2d->order) {
    tensor->mode_ordering[0] = t2d->order[0];
    tensor->mode_ordering[1] = t2d->order[1];
  } else {
    tensor->mode_ordering[0] = 0;
    tensor->mode_ordering[1] = 1;
  }
  tensor->mode_types = (taco_mode_t *)malloc(2 * sizeof(taco_mode_t));
  tensor->indices = (uint8_t ***)malloc(2 * sizeof(uint8_t ***));

  // allocate memory for dense indices and store information
  // on dense indices for row
  tensor->indices[0] = (uint8_t **)malloc(1 * sizeof(uint8_t **));
  (tensor->indices[0])[0] = (uint8_t *)malloc(sizeof(uint8_t *));
  *(tensor->indices[0][0]) = t2d->nr;

  // allocate memory for compressed indices and store information
  // on compressed indices rptr, colidx
  tensor->indices[1] = (uint8_t **)malloc(1 * sizeof(uint8_t **));
  (tensor->indices[1])[0] = (uint8_t *)malloc(sizeof(uint8_t *));
  *(tensor->indices[1][0]) = t2d->nr;

  // m is dense sparse
  tensor->mode_types[0] = taco_mode_dense;
  tensor->mode_types[1] = taco_mode_dense;

  tensor->vals_size = t2d->nr * t2d->nc;
  tensor->vals = (uint8_t *)t2d->vals;
  return tensor;
}

static inline void taco2matrix(taco_tensor_t *t, tensor2d *t2d) {
  t2d->order[0] = t->mode_ordering[0];
  t2d->order[1] = t->mode_ordering[1];
  t2d->nr = t->dimensions[0];
  t2d->nc = t->dimensions[1];
  t2d->vals = (double *)t->vals;
}

static inline taco_tensor_t *tensor4d2taco(tensor4d *t4d) {
  taco_tensor_t *tensor = (taco_tensor_t *)malloc(sizeof(taco_tensor_t));
  tensor->order = 2;
  tensor->dimensions = (int32_t *)malloc(sizeof(int32_t) * 4);
  tensor->dimensions[0] = t4d->ni;
  tensor->dimensions[1] = t4d->nj;
  tensor->dimensions[2] = t4d->nk;
  tensor->dimensions[3] = t4d->nl;
  tensor->csize = sizeof(*t4d->vals);
  tensor->mode_ordering = (int32_t *)malloc(sizeof(int32_t) * 4);
  // storage data layout order
  if (t4d->order) {
    tensor->mode_ordering[0] = t4d->order[0];
    tensor->mode_ordering[1] = t4d->order[1];
    tensor->mode_ordering[2] = t4d->order[2];
    tensor->mode_ordering[3] = t4d->order[3];
  } else {
    tensor->mode_ordering[0] = 0;
    tensor->mode_ordering[1] = 1;
    tensor->mode_ordering[2] = 2;
    tensor->mode_ordering[3] = 3;
  }
  tensor->mode_types = (taco_mode_t *)malloc(4 * sizeof(taco_mode_t));
  tensor->indices = (uint8_t ***)malloc(4 * sizeof(uint8_t ***));

  // allocate memory for dense indices and store information
  // on dense indices for row
  tensor->indices[0] = (uint8_t **)malloc(1 * sizeof(uint8_t **));
  (tensor->indices[0])[0] = (uint8_t *)malloc(sizeof(uint8_t *));
  *(tensor->indices[0][0]) = t4d->ni;

  // allocate memory for compressed indices and store information
  // on compressed indices rptr, colidx
  tensor->indices[1] = (uint8_t **)malloc(1 * sizeof(uint8_t **));
  (tensor->indices[1])[0] = (uint8_t *)malloc(sizeof(uint8_t *));
  *(tensor->indices[1][0]) = t4d->nj;

  tensor->indices[2] = (uint8_t **)malloc(1 * sizeof(uint8_t **));
  (tensor->indices[2])[0] = (uint8_t *)malloc(sizeof(uint8_t *));
  *(tensor->indices[2][0]) = t4d->nk;


  tensor->indices[3] = (uint8_t **)malloc(1 * sizeof(uint8_t **));
  (tensor->indices[3])[0] = (uint8_t *)malloc(sizeof(uint8_t *));
  *(tensor->indices[3][0]) = t4d->nl;



  // m is dense sparse
  tensor->mode_types[0] = taco_mode_dense;
  tensor->mode_types[1] = taco_mode_dense;
  tensor->mode_types[2] = taco_mode_dense;
  tensor->mode_types[3] = taco_mode_dense;

  tensor->vals_size = t4d->ni * t4d->nj * t4d->nk * t4d->nl;
  tensor->vals = (uint8_t *)t4d->vals;
  return tensor;
}

static inline void taco2matrix(taco_tensor_t *t, tensor4d *t4d) {
  t4d->order[0] = t->mode_ordering[0];
  t4d->order[1] = t->mode_ordering[1];
  t4d->nr = t->dimensions[0];
  t4d->nc = t->dimensions[1];
  t4d->vals = (double *)t->vals;
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
    t3_tensor.vals = t3;

    tensor2d t1_tensor;
    t1_tensor.nr = p4u;
    t1_tensor.nc = h1u;
    t1_tensor.vals = t1;

    tensor4d v2_tensor;
    v2_tensor.ni = h3u;
    v2_tensor.nj = h2u;
    v2_tensor.nk = p6u;
    v2_tensor.nl = p5u;
    v2_tensor.vals = v2;

    c1d_sd_t_s1_1_taco(&t3_tensor,&t1_tensor,&v2_tensor);

    
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
    tensor6d t3_tensor;
    t3_tensor.ni = h3u;
    t3_tensor.nj = h1u; 
    t3_tensor.nk = h2u;
    t3_tensor.nl = p6u;
    t3_tensor.nm = p5u;
    t3_tensor.nn = p4u;
    t3_tensor.vals = t3;

    tensor2d t1_tensor;
    t1_tensor.nr = p4u;
    t1_tensor.nc = h1u;
    t1_tensor.vals = t1;

    tensor4d v2_tensor;
    v2_tensor.ni = h3u;
    v2_tensor.nj = h2u;
    v2_tensor.nk = p6u;
    v2_tensor.nl = p5u;
    v2_tensor.vals = v2;
    
    scalar alpha;
    alpha.val = -1;
    c1d_sd_t_s1_2_taco(&t3_tensor,&t1_tensor,&v2_tensor,&alpha);

    //t3(h3,h1,h2,p6,p5,p4)-=t1(p4,h1)*v2(h3,h2,p6,p5);
    
    
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
    
    tensor6d t3_tensor;
    t3_tensor.ni = h1u;
    t3_tensor.nj = h3u; 
    t3_tensor.nk = h2u;
    t3_tensor.nl = p6u;
    t3_tensor.nm = p5u;
    t3_tensor.nn = p4u;
    t3_tensor.vals = t3;

    tensor2d t1_tensor;
    t1_tensor.nr = p4u;
    t1_tensor.nc = h1u;
    t1_tensor.vals = t1;

    tensor4d v2_tensor;
    v2_tensor.ni = h3u;
    v2_tensor.nj = h2u;
    v2_tensor.nk = p6u;
    v2_tensor.nl = p5u;
    v2_tensor.vals = v2;
    
    c1d_sd_t_s1_3_taco(&t3_tensor,&t1_tensor,&v2_tensor);
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
    
    tensor6d t3_tensor;
    t3_tensor.ni = h3u;
    t3_tensor.nj = h2u; 
    t3_tensor.nk = h1u;
    t3_tensor.nl = p6u;
    t3_tensor.nm = p4u;
    t3_tensor.nn = p5u;
    t3_tensor.vals = t3;

    tensor2d t1_tensor;
    t1_tensor.nr = p4u;
    t1_tensor.nc = h1u;
    t1_tensor.vals = t1;

    tensor4d v2_tensor;
    v2_tensor.ni = h3u;
    v2_tensor.nj = h2u;
    v2_tensor.nk = p6u;
    v2_tensor.nl = p5u;
    v2_tensor.vals = v2;
    
    scalar alpha;
    alpha.val = -1;
    c1d_sd_t_s1_4_taco(&t3_tensor,&t1_tensor,&v2_tensor,&alpha);


    //t3(h3,h2,h1,p6,p4,p5)-=t1(p4,h1)*v2(h3,h2,p6,p5);

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
    

    tensor6d t3_tensor;
    t3_tensor.ni = h3u;
    t3_tensor.nj = h1u; 
    t3_tensor.nk = h2u;
    t3_tensor.nl = p6u;
    t3_tensor.nm = p4u;
    t3_tensor.nn = p5u;
    t3_tensor.vals = t3;

    tensor2d t1_tensor;
    t1_tensor.nr = p4u;
    t1_tensor.nc = h1u;
    t1_tensor.vals = t1;

    tensor4d v2_tensor;
    v2_tensor.ni = h3u;
    v2_tensor.nj = h2u;
    v2_tensor.nk = p6u;
    v2_tensor.nl = p5u;
    v2_tensor.vals = v2;
    
    //t3(h3,h1,h2,p6,p4,p5)+=t1(p4,h1)*v2(h3,h2,p6,p5);
    c1d_sd_t_s1_5_taco(&t3_tensor,&t1_tensor,&v2_tensor);
    
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
    
    tensor6d t3_tensor;
    t3_tensor.ni = h1u;
    t3_tensor.nj = h3u; 
    t3_tensor.nk = h2u;
    t3_tensor.nl = p6u;
    t3_tensor.nm = p4u;
    t3_tensor.nn = p5u;
    t3_tensor.vals = t3;

    tensor2d t1_tensor;
    t1_tensor.nr = p4u;
    t1_tensor.nc = h1u;
    t1_tensor.vals = t1;

    tensor4d v2_tensor;
    v2_tensor.ni = h3u;
    v2_tensor.nj = h2u;
    v2_tensor.nk = p6u;
    v2_tensor.nl = p5u;
    v2_tensor.vals = v2;
    
    scalar alpha;
    alpha.val = -1;
    //t3(h1,h3,h2,p6,p4,p5)-=t1(p4,h1)*v2(h3,h2,p6,p5);
    c1d_sd_t_s1_6_taco(&t3_tensor,&t1_tensor,&v2_tensor,&alpha);
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
    
    tensor6d t3_tensor;
    t3_tensor.ni = h3u;
    t3_tensor.nj = h2u; 
    t3_tensor.nk = h1u;
    t3_tensor.nl = p4u;
    t3_tensor.nm = p6u;
    t3_tensor.nn = p5u;
    t3_tensor.vals = t3;

    tensor2d t1_tensor;
    t1_tensor.nr = p4u;
    t1_tensor.nc = h1u;
    t1_tensor.vals = t1;

    tensor4d v2_tensor;
    v2_tensor.ni = h3u;
    v2_tensor.nj = h2u;
    v2_tensor.nk = p6u;
    v2_tensor.nl = p5u;
    v2_tensor.vals = v2;
    //t3(h3,h2,h1,p4,p6,p5)+=t1(p4,h1)*v2(h3,h2,p6,p5);
    c1d_sd_t_s1_7_taco(&t3_tensor,&t1_tensor,&v2_tensor);
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

    tensor6d t3_tensor;
    t3_tensor.ni = h1u;
    t3_tensor.nj = h3u; 
    t3_tensor.nk = h2u;
    t3_tensor.nl = p6u;
    t3_tensor.nm = p4u;
    t3_tensor.nn = p5u;
    t3_tensor.vals = t3;

    tensor2d t1_tensor;
    t1_tensor.nr = p4u;
    t1_tensor.nc = h1u;
    t1_tensor.vals = t1;

    tensor4d v2_tensor;
    v2_tensor.ni = h3u;
    v2_tensor.nj = h2u;
    v2_tensor.nk = p6u;
    v2_tensor.nl = p5u;
    v2_tensor.vals = v2;
    
    scalar alpha;
    alpha.val = -1;
   //t3(h3,h1,h2,p4,p6,p5)-=t1(p4,h1)*v2(h3,h2,p6,p5);
    c1d_sd_t_s1_8_taco(&t3_tensor,&t1_tensor,&v2_tensor,&alpha);
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

    tensor6d t3_tensor;
    t3_tensor.ni = h1u;
    t3_tensor.nj = h3u; 
    t3_tensor.nk = h2u;
    t3_tensor.nl = p4u;
    t3_tensor.nm = p6u;
    t3_tensor.nn = p5u;
    t3_tensor.vals = t3;

    tensor2d t1_tensor;
    t1_tensor.nr = p4u;
    t1_tensor.nc = h1u;
    t1_tensor.vals = t1;

    tensor4d v2_tensor;
    v2_tensor.ni = h3u;
    v2_tensor.nj = h2u;
    v2_tensor.nk = p6u;
    v2_tensor.nl = p5u;
    v2_tensor.vals = v2;
   //t3(h1,h3,h2,p4,p6,p5)+=t1(p4,h1)*v2(h3,h2,p6,p5);
    c1d_sd_t_s1_9_taco(&t3_tensor,&t1_tensor,&v2_tensor);
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
