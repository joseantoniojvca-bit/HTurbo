#ifndef resconst_h
#define resconst_h

namespace resconst
{
  extern void init();
  extern void initcmw();

  //number of light flavours
  extern const int NF;
  //number of quarks up, down and leptons
  extern const int nups;
  extern const int ndowns;
  extern const int nleps;
  //square of the quark charge
  extern const double equ2;
  extern const double eqd2;
  // Number of flavors weighed with charges
  extern const double NFq2;
  extern const double NFq3;
  extern const double NFq4;
  extern const double NFq6;
  // Value of alpha_qed(0)
  extern const double a0qed;
  // Value of pcF coefficient (Hst) DY QED
  extern const double pcF_QED;


  //Mathematics constants
  extern const double Z2;
  extern const double Z3;
  extern const double Z4;
  extern const double Z5;
  extern const double b0;
  
  //Resummation coefficients
  extern double beta0, beta1, beta2, beta3, beta4, beta5, Kappa;
  extern double beta0p, beta1p, beta01_QCD, beta01_QED;
  extern double A1g, A2g, A3g, A4g, B1g, B2g, B3g, C1ggn;
  extern double A1q, A2q, A3q, A4q, A5q, A6q;
  extern double A1qp, A2qp;
  extern double B1q, B2q, B3q, B4q, B5q;
  extern double B1qp, B2qp, D1w_QED;
  extern double B1qHSCHM, B2qHSCHM, B3qHSCHM;
  extern double H1q, H2q, H3q, H4q, H1g, H2g, H3g;
  extern double H1qp;
  extern double K1, K2, K3, K4;
  extern double C1qqn, C2qqn;
  extern double C1qqnp;
  extern double s1, s2, s3, s4;
  extern double k1, k2, k3, k4, k3_nfv, k4_nfv;

  // Delta terms
  extern double C1ggdelta, C1qqdelta, Delta2qq, Delta2gg;

  // Coefficients of D0 and D1 in P*P (as/M_PI normalization)
  extern double D0qqqq, D1qqqq, D0gggg, D1gggg;

  // Coefficients of delta(1-z) in P*P
  extern double Deltaqqqq, Deltagggg;

  // H2qq contribution: coefficient of delta(1-z)
  extern double H2qqdelta;

  // H2qq contribution: coefficient of D0(z)
  extern double H2qqD0;
}

#endif
