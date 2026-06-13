#ifndef mellinint_h
#define mellinint_h
#include "fcomplex.h"

#include <map>
#include <complex>

//fortran interfaces
extern "C"
{
  void mellinint_pdf_mesq_expy_(int& i1, int& i2, int& sign);
  fcomplex mellinint_integrand_(int& i1, int& i2, int& sign);
}

using namespace std;
namespace mellinint
{
  extern int mdim;
  extern complex <double> CCp,CCm; //angles of the positive and negative branches in the unitary complex circle
  extern double cphi,sphi;
  
  extern complex <double> *wn; //weights of the gaussian quadrature rule
  extern complex <double> *Np; //nodes on the positive branch of the contour
  extern complex <double> *Nm; //nodes on the negative branch of the contour
#pragma omp threadprivate(wn,Np,Nm)

  //mellin 2d case
  extern complex <double> *wn_1, *wn_2; //weights of the gaussian quadrature rule
  extern complex <double> *Np_1, *Np_2; //nodes on the positive branch of the contour
  extern complex <double> *Nm_1, *Nm_2; //nodes on the negative branch of the contour
#pragma omp threadprivate(wn_1,Np_1,Nm_1,wn_2,Np_2,Nm_2)

  extern void initgauss();
  extern void update();
  extern void updategauss();
  extern void allocate();
  extern void free();
  extern void release();

  extern double Not;
  extern double Not_1;
  extern double Not_2;
  extern double sigma;
  extern double lambda;
  extern double nu;
  extern double alpha;
  extern double mu;
  //pragma omp threadprivate(Not,Not_1,Not_2)
  
  extern bool midpoint;
  extern bool weideman;
  
  extern void pdf_mesq_expy(int i1, int i2, int sign);
  extern double integrand2d(int i1, int i2, int sign);
  extern double integrand1d(int i);

  extern complex <double> calc1d();
  extern complex <double> calc1d(double eq2, std::string sector);
  // Andrea: trial for W boson production, no up-down separation 
  extern complex <double> calc1d(double eq2, double eqbarp2);
  extern complex <double> calc1d_muf();
  //
  extern complex <double> calc2d();
  extern complex <double> calc2d(double eq2, std::string sector);
  // Andrea for W
  extern complex <double> calc2d(double eq2, double eqbarp2);
  extern complex <double> calc2d_muf();
  //
  extern void reset();

  inline int index(int i1, int i2)
  {return i2 + mellinint::mdim*i1;}
  
  //luminosity times mesqij
  extern complex<double> GGN;
  extern complex<double> QGN_1;
  extern complex<double> QGN_2;
  extern complex<double> GQN_1;
  extern complex<double> GQN_2;
  extern complex<double> QGNUP_1;
  extern complex<double> QGNUP_2;
  extern complex<double> QGNDW_1;
  extern complex<double> QGNDW_2;
  extern complex<double> QANUP_1;
  extern complex<double> QANUP_2;
  extern complex<double> QANDW_1;
  extern complex<double> QANDW_2;
// ANDREA: photon channels for W 
  extern complex<double> QAN_1;
  extern complex<double> QAN_2;
  /// // 
  extern complex<double> QQBN;
  extern complex<double> QQBNUP;
  extern complex<double> QQBNDW;
  extern complex<double> QQN;
  extern complex<double> QQN_1;
  extern complex<double> QQN_2;
  extern complex<double> QQPN_1;
  extern complex<double> QQPN_2;
#pragma omp threadprivate(GGN,QGN_1,QGN_2,QGNUP_1,QGNUP_2,QGNDW_1,QGNDW_2,QANUP_1,QANUP_2,QANDW_1,QANDW_2,QAN_1,QAN_2,QQBN,QQBNUP,QQBNDW,QQN,QQN_1,QQN_2,QQPN_1,QQPN_2)

  extern complex<double> GGN_bb;
#pragma omp threadprivate(GGN_bb)
  
  extern complex <double> QQBN_nfz;
  extern complex <double> QQBPN_1;
  extern complex <double> QQBPN_2;
  extern complex <double> QBGN_1;
  extern complex <double> QBGN_2;
  extern complex <double> QPGN_1;
  extern complex <double> QPGN_2;
  extern complex <double> QBPGN_1;
  extern complex <double> QBPGN_2;
#pragma omp threadprivate(QQBN_nfz,QQBPN_1,QQBPN_2,QBGN_1,QBGN_2,QPGN_1,QPGN_2,QBPGN_1,QBPGN_2)

  extern complex <double> QBQN;
  extern complex <double> QBQPN_2;
  extern complex <double> QBQPN_1;
  extern complex <double> QBQBPN_2;
  extern complex <double> QBQBPN_1;
  extern complex <double> QPQPN_2;
  extern complex <double> QPQPN_1;
  extern complex <double> QPQBPN;
  extern complex <double> QBPQPN;
#pragma omp threadprivate(QBQN,QBQPN_2,QBQPN_1,QBQBPN_2,QBQBPN_1,QPQPN_2,QPQPN_1,QPQBPN,QBPQPN)

  //HF
  extern complex<double> BBRN,CCBN;
#pragma omp threadprivate(BBRN,CCBN)

  extern complex<double> bbfun, ccfun;
#pragma omp threadprivate(bbfun,ccfun)
  
}

#endif
