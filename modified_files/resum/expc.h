#ifndef expc_h
#define expc_h

//#include "mellinint.h"

#include <complex>

using namespace std;

namespace expc
{
  extern complex <double> aexpH;
  extern complex <double> *aexpgg;

  extern complex <double> *aexpqq;
  extern complex <double> *aexpqg;
  extern complex <double> *aexpqqqed;
  extern complex <double> *aexpqgamma;
  extern complex <double> *aexpqqb;
  extern complex <double> *aexpqqp;
  extern complex <double> *aexpqqbp;
  extern complex <double> *aexpqq_1;
  extern complex <double> *aexpqg_1;
  extern complex <double> *aexpgq_1;
  extern complex <double> *aexpqqqed_1;
  extern complex <double> *aexpqgamma_1;
  extern complex <double> *aexpqqb_1;
  extern complex <double> *aexpqqp_1;
  extern complex <double> *aexpqqbp_1;
  extern complex <double> *aexpqq_2;
  extern complex <double> *aexpqg_2;
  extern complex <double> *aexpgq_2;
  extern complex <double> *aexpqqqed_2;
  extern complex <double> *aexpqgamma_2;
  extern complex <double> *aexpqqb_2;
  extern complex <double> *aexpqqp_2;
  extern complex <double> *aexpqqbp_2;
#pragma omp threadprivate(aexpqq,aexpqg,aexpqqqed,aexpqgamma,aexpqqb,aexpqqp,aexpqqbp,aexpqq_1,aexpqg_1,aexpqqqed_1,aexpqgamma_1,aexpqqb_1,aexpqqp_1,aexpqqbp_1,aexpqq_2,aexpqg_2,aexpqqqed_2,aexpqgamma_2,aexpqqb_2,aexpqqp_2,aexpqqbp_2)

  extern complex <double> *aexpqg_bb;
  extern complex <double> *aexpqg_1_bb;
  extern complex <double> *aexpqg_2_bb;
#pragma omp threadprivate(aexpqg_bb,aexpqg_1_bb,aexpqg_2_bb)
 
  extern complex <double> *qqb;
  extern complex <double> *qg;
  extern complex <double> *qqbqed;
  extern complex <double> *qgamma;
  extern complex <double> *qg_1;
  extern complex <double> *qg_2;
  extern complex <double> *qgamma_1;
  extern complex <double> *qgamma_2;
  extern complex <double> *gq_1;
  extern complex <double> *gq_2;
  extern complex <double> *qq;
  extern complex <double> *qq_1;
  extern complex <double> *qq_2;
  extern complex <double> *qqp;
  extern complex <double> *qqp_1;
  extern complex <double> *qqp_2;
  extern complex <double> *qqbp;
  extern complex <double> *qqbp_1;
  extern complex <double> *qqbp_2;
  extern complex <double> *gg;
  extern complex <double> *qbg;
  extern complex <double> *qbg_1;
  extern complex <double> *qbg_2;
  extern complex <double> *qpg;
  extern complex <double> *qpg_1;
  extern complex <double> *qpg_2;
  extern complex <double> *qbpg;
  extern complex <double> *qbpg_1;
  extern complex <double> *qbpg_2;
  extern complex <double> *qbq;
  extern complex <double> *qbqp;
  extern complex <double> *qbqp_1;
  extern complex <double> *qbqp_2;
  extern complex <double> *qbqbp;
  extern complex <double> *qbqbp_1;
  extern complex <double> *qbqbp_2;
  extern complex <double> *qpqp;
  extern complex <double> *qpqp_1;
  extern complex <double> *qpqp_2;
  extern complex <double> *qpqbp;
  extern complex <double> *qbpqp;
#pragma omp threadprivate(qqb,qg,qg_1,qg_2,qqbqed,qgamma,qgamma_1,qgamma_2,qq,qq_1,qq_2,qqp,qqp_1,qqp_2,qqbp,qqbp_1,qqbp_2,gg,qbg,qbg_1,qbg_2,qpg,qpg_1,qpg_2,qbpg,qbpg_1,qbpg_2,qbq,qbqp,qbqp_1,qbqp_2,qbqbp,qbqbp_1,qbqbp_2,qpqp,qpqp_1,qpqp_2,qpqbp,qbpqp)

  extern complex <double> *gg_bb;
#pragma omp threadprivate(gg_bb)
  
  //  extern complex <double> aexp;
  extern complex <double> xlambda;
  extern complex <double> lamB,aexpB;
  extern complex <double> xlambdaqed;
  extern complex <double> lamBqed,aexpBqed;
  extern complex <double> lamC,aexpC;
  extern complex <double> lamD,aexpD;
  extern complex <double> lamE,aexpE;
  extern complex <double> lamF,aexpF;
  extern complex <double> lamG,aexpG;
  extern complex <double> lamH;
#pragma omp threadprivate(xlambda,lamB,xlambdaqed,lamBqed,lamC,lamD,lamE,lamF,lamG,lamH)
#pragma omp threadprivate(aexpB,aexpBqed,aexpC,aexpD,aexpE,aexpF,aexpG,aexpH)

  extern complex <double> lamB_bb,aexpB_bb;
#pragma omp threadprivate(lamB_bb)
#pragma omp threadprivate(aexpB_bb)
  
  extern void allocate();
  extern void reset();
  extern void calc(complex <double> b);
  extern void free();

  extern void noexp();
  extern void delta();
  extern void diag();
  extern void full();
  extern void fullqed();
  extern void denlow();
  extern void denfull();
  extern void css(complex <double> b);
  extern void matrix();
  //extern void taylor();

}

#endif
