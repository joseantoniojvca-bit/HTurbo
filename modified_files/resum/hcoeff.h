#ifndef hcoeff_h
#define hcoeff_h

#include "mellinint.h"


#include <complex>

using namespace std;

double H1intHSchm(double costh, double s);

namespace hcoeff
{
//  extern complex <double> *Hqqb;
//  extern complex <double> *Hqg;
//  extern complex <double> *Hqq;
//  extern complex <double> *Hqqp;
//  extern complex <double> *Hqqbp;
//  extern complex <double> *Hgg;
//  extern complex <double> *Hqbg;
//  extern complex <double> *Hqpg;
//  extern complex <double> *Hqbpg;
//#pragma omp threadprivate(Hqqb,Hqg,Hqq,Hqqp,Hqqbp,Hgg,Hqbg,Hqpg,Hqbpg)

  extern complex <double> *Hqqb;
  extern complex <double> *Hqqb_qed;
  extern complex <double> *Hqg;
  extern complex <double> *Hqg_1;
  extern complex <double> *Hqg_2;
  extern complex <double> *Hgq_1;
  extern complex <double> *Hgq_2;
  extern complex <double> *Hqgamma;
  extern complex <double> *Hqgamma_1;
  extern complex <double> *Hqgamma_2;
  extern complex <double> *Hqq;
  extern complex <double> *Hqq_1;
  extern complex <double> *Hqq_2;
  extern complex <double> *Hqqp;
  extern complex <double> *Hqqp_1;
  extern complex <double> *Hqqp_2;
  extern complex <double> *Hqqbp;
  extern complex <double> *Hqqbp_1;
  extern complex <double> *Hqqbp_2;
  extern complex <double> *Hgg;
  extern complex <double> *Hqbg;
  extern complex <double> *Hqbg_1;
  extern complex <double> *Hqbg_2;
  extern complex <double> *Hqpg;
  extern complex <double> *Hqpg_1;
  extern complex <double> *Hqpg_2;
  extern complex <double> *Hqbpg;
  extern complex <double> *Hqbpg_1;
  extern complex <double> *Hqbpg_2;
  extern complex <double> Hqqb_nfz;
  extern complex <double> *Hqbq;
  extern complex <double> *Hqbqp;
  extern complex <double> *Hqbqp_1;
  extern complex <double> *Hqbqp_2;
  extern complex <double> *Hqbqbp;
  extern complex <double> *Hqbqbp_1;
  extern complex <double> *Hqbqbp_2;
  extern complex <double> *Hqpqp;
  extern complex <double> *Hqpqp_1;
  extern complex <double> *Hqpqp_2;
  extern complex <double> *Hqpqbp;
  extern complex <double> *Hqbpqp;
#pragma omp threadprivate(Hqqb,Hqqb_qed,Hqg,Hqg_1,Hqg_2,Hqgamma,Hqgamma_1,Hqgamma_2,Hqq,Hqq_1,Hqq_2,Hqqp,Hqqp_1,Hqqp_2,Hqqbp,Hqqbp_1,Hqqbp_2,Hgg,Hqbg,Hqbg_1,Hqbg_2,Hqpg,Hqpg_1,Hqpg_2,Hqbpg,Hqbpg_1,Hqbpg_2,Hqqb_nfz,Hqbq,Hqbqp,Hqbqp_1,Hqbqp_2,Hqbqbp,Hqbqbp_1,Hqbqbp_2,Hqpqp,Hqpqp_1,Hqpqp_2,Hqpqbp,Hqbpqp)
  
//  extern complex <double> *H1stqqb;
//  extern complex <double> *H1stqg;
//  extern complex <double> H1stgg;
//
//  extern complex <double> *H2stqq;
//  extern complex <double> *H2stqqp;
//  extern complex <double> *H2stqqb;
//  extern complex <double> *H2stqg;
//  extern complex <double> *H2stgg;
//
//  extern complex <double> *aexpqq;
//  extern complex <double> *aexpqg;
//#pragma omp threadprivate(H1stgg,H1stqg,H1stqqb,H2stqq,H2stqqp,H2stqqb,H2stqg,H2stgg,aexpqq,aexpqg)


  extern complex <double> Hqqbz;
  extern complex <double> Hqqb_qed_z;
  extern complex <double> Hqgz;
  extern complex <double> Hqgammaz;
  extern complex <double> Hqqz;
  extern complex <double> Hqqpz;
  extern complex <double> Hggz;
#pragma omp threadprivate(Hqqbz,Hqqb_qed_z,Hqgz,Hqgammaz,Hqqz,Hqqpz,Hggz)

  inline int index(int i1, int i2, int sign)
  {return i2 + mellinint::mdim*(i1 + mellinint::mdim*sign);}
  
  extern void allocate();

  extern void reset();
  
  // q2-dependent quantities
  extern void calc();
  extern void calc(double m2, double costh);

  // b-dependent quantities
  extern void calcb(double aass, complex <double> logmuf2q2, complex <double> loga, complex <double> alpq, complex <double> aexp, complex <double> aexpb);

  extern void truncate();
    
  extern void invert(double q2);
  
  extern void free();
}

#endif
