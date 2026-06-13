#include "ccoeff.h"
#include "dynnlo_interface.h"
#include "dyres_interface.h"
#include "anomalous.h"
#include "mellinint.h"
#include "mesq.h"
#include "settings.h"
#include "gaussrules.h"
#include "resconst.h"
#include "constants.h"
#include "icoeff.h"
#include "phasespace.h"
//#include "hs.h"
#include "psi.h"
#include "cmom.h"
#include "melfun.h"
#include "resint.h"
#include "accel.h"
#include <iostream>
#include <iomanip>


//fotran interface for dymellinh2.f
extern "C"
{
  fcomplex fun1_(fcomplex &xn)  {return fcx(fun1(cx(xn)));};
  fcomplex fun2_(fcomplex &xn)  {return fcx(fun2(cx(xn)));};
  fcomplex fun3_(fcomplex &xn)  {return fcx(fun3(cx(xn)));};
  fcomplex fun4_(fcomplex &xn)  {return fcx(fun4(cx(xn)));};
  fcomplex fun5_(fcomplex &xn)  {return fcx(fun5(cx(xn)));};
  fcomplex fun6_(fcomplex &xn)  {return fcx(fun6(cx(xn)));};
  fcomplex fun7_(fcomplex &xn)  {return fcx(fun7(cx(xn)));};
  fcomplex fun8_(fcomplex &xn)  {return fcx(fun8(cx(xn)));};
  fcomplex fun9_(fcomplex &xn)  {return fcx(fun9(cx(xn)));};
  fcomplex fun10_(fcomplex &xn) {return fcx(fun10(cx(xn)));};
  fcomplex fun11_(fcomplex &xn) {return fcx(fun11(cx(xn)));};
  fcomplex fun12_(fcomplex &xn) {return fcx(fun12(cx(xn)));};
  fcomplex fun13_(fcomplex &xn) {return fcx(fun13(cx(xn)));};
  fcomplex fun14_(fcomplex &xn) {return fcx(fun14(cx(xn)));};
  fcomplex fun15_(fcomplex &xn) {return fcx(fun15(cx(xn)));};
  fcomplex fun16_(fcomplex &xn) {return fcx(fun16(cx(xn)));};

  fcomplex dypsi0_(fcomplex &xn)  {return fcx(psi0(cx(xn)));};
  fcomplex dypsi1_(fcomplex &xn)  {return fcx(psi1(cx(xn)));};
  fcomplex dypsi2_(fcomplex &xn)  {return fcx(psi2(cx(xn)));};
  fcomplex dypsi3_(fcomplex &xn)  {return fcx(psi3(cx(xn)));};
  fcomplex dybe0_(fcomplex &xn)  {return fcx(be0(cx(xn)));};
  fcomplex dybe1_(fcomplex &xn)  {return fcx(be1(cx(xn)));};
  fcomplex dybe2_(fcomplex &xn)  {return fcx(be2(cx(xn)));};
  fcomplex dybe3_(fcomplex &xn)  {return fcx(be3(cx(xn)));};
}

complex <double> c1qgreg(complex <double> N)
{
   complex <double> Z1_,Z2_,Z3_,Z4_,Z5_,Z6_,Z7_,Z8_,Z9_,Z10_,Z11_,Z12_,Z13_,Z14_,Z15_,Z16_,Z17_,Z18_,Z19_,Z20_,Z21_,Z22_,Z23_,Z24_,Z25_,Z26_;

    Z1_=pow(N+1.,-1);
    Z2_=pow(N+2.,-1);
    Z3_=pow(N+3.,-1);
    Z4_=pow(N+4.,-1);
    Z5_=pow(N+5.,-1);
    Z6_=pow(N+6.,-1);
    Z7_=pow(N+7.,-1);
    Z8_=pow(N+8.,-1);
    Z9_=pow(N+9.,-1);
    Z10_=pow(N+10.,-1);
    Z11_=log(N+10.);
   Z12_= - 7.45739E+5 + 4.461336985917776E+5*Z10_;
   Z12_=Z12_*Z10_;
   Z12_=Z12_ - 5.172586451199571E+4;
   Z13_=2.784739626844890E-4*Z10_;
   Z12_=Z12_*Z13_;
   Z12_=Z12_ + 2.659558823529411E+1;
   Z14_=3.333333333333333E-1*Z10_;
   Z12_=Z12_*Z14_;
   Z12_=Z12_ + 7.013745687347923E-1;
   Z12_=Z12_*Z10_;
   Z12_=Z12_ - 1.444444444444444E+0;
   Z12_=Z12_*Z10_;
   Z12_=Z12_ - 1.331698776556260E-1;
   Z12_=Z12_*Z10_;
   Z12_=Z12_ + 3.093610093610093E-1;
   Z12_=Z12_*Z10_;
   Z12_=Z12_ + 3.426400318573169E-2;
   Z12_=Z12_*Z10_;
   Z12_=Z12_ - 9.090909090909090E-2;
   Z12_=Z12_*Z10_;
   Z12_=Z12_ - 1.267503187070597E-2;
   Z12_=Z12_*Z10_;
   Z12_=Z12_ + 3.888888888888888E-2;
   Z12_=Z12_*Z10_;
   Z12_=Z12_ + 7.387934195554950E-3;
   Z12_=Z12_*Z10_;
   Z12_=Z12_ - 2.645502645502645E-2;
   Z12_=Z12_*Z10_;
   Z12_=Z12_ - 7.959954125942146E-3;
   Z12_=Z12_*Z10_;
   Z12_=Z12_ + 3.333333333333333E-2;
   Z12_=Z12_*Z10_;
   Z12_=Z12_ + 2.208627403484887E-2;
   Z12_=Z12_*Z10_;
   Z12_=Z12_ - 1.111111111111111E-1;
   Z12_=Z12_*Z10_;
   Z12_=Z12_ - 4.616034810892295E-1;
   Z12_=Z12_*Z10_;
   Z12_=Z12_ - 2.102954219868710E+0;
   Z12_=Z12_*Z10_;
   Z12_=Z12_ + 2.637482654207926E+0;
   Z12_=2.5E-1*Z12_;
   Z15_=pow(Z10_,2);
   Z16_= - 3.664745196324143E+1 + 3.174745454545454E+2*Z15_;
   Z16_=Z16_*Z15_;
   Z16_=Z16_ + 5.319117647058823E+0;
   Z16_=Z16_*Z15_;
   Z16_=Z16_ - 1.E+0;
   Z16_=Z16_*Z15_;
   Z16_=Z16_ + 2.531135531135531E-1;
   Z16_=Z16_*Z15_;
   Z16_=Z16_ - 9.090909090909090E-2;
   Z16_=Z16_*Z15_;
   Z16_=Z16_ + 5.E-2;
   Z16_=Z16_*Z15_;
   Z16_=Z16_ - 4.761904761904761E-2;
   Z16_=Z16_*Z15_;
   Z16_=Z16_ + 1.E-1;
   Z15_=Z16_*Z15_;
   Z15_=Z15_ - 1.E+0;
   Z16_=1.666666666666666E-1*Z10_;
   Z15_=Z15_*Z16_;
   Z15_=Z15_ - 1.E+0;
   Z15_=Z15_*Z10_;
   Z15_=Z15_ + 1.154431329803065E+0;
   Z16_=Z11_ + Z15_;
   Z17_=3.333333333333333E-1*Z11_;
   Z16_=Z16_*Z17_;
   Z16_=Z12_ + Z16_;
   Z18_=5.E-1*Z16_;
   Z15_=5.E-1*Z15_;
   Z19_=Z4_ - Z11_;
   Z20_=Z15_ - Z19_;
   Z21_=3.333333333333333E-1*Z8_;
   Z21_=Z20_*Z21_;
   Z22_=Z15_ + Z11_;
   Z23_=3.333333333333333E-1*Z4_;
   Z22_=Z22_*Z23_;
   Z21_=Z21_ + Z22_;
   Z19_=Z3_ + Z19_ + Z8_;
   Z23_=Z15_ - Z19_;
   Z24_=3.333333333333333E-1*Z9_;
   Z23_=Z23_*Z24_;
   Z20_= - Z8_ + Z20_;
   Z24_=3.333333333333333E-1*Z3_;
   Z20_=Z20_*Z24_;
   Z23_=Z18_ - Z21_ - Z23_ - Z20_;
   Z24_= - Z5_*Z23_;
   Z19_=Z15_ - Z19_ - Z9_;
   Z25_=3.333333333333333E-1*Z5_;
   Z25_=Z19_*Z25_;
   Z23_=Z23_ - Z25_;
   Z25_= - Z6_*Z23_;
   Z19_= - Z5_ + Z19_;
   Z26_=3.333333333333333E-1*Z6_;
   Z26_=Z19_*Z26_;
   Z23_= - Z26_ + Z23_;
   Z26_= - Z7_*Z23_;
   Z19_= - Z6_ + Z19_;
   Z19_=Z7_*Z19_;
   Z19_=3.333333333333333E-1*Z19_ - Z23_;
   Z19_=Z2_*Z19_;
   Z15_=Z15_ + Z17_;
   Z15_=Z15_*Z17_;
   Z12_=Z12_ + Z15_;
   Z12_=Z11_*Z12_;
   Z15_= - Z4_*Z16_;
   Z16_= - 1.188595953600946E+49 - 1.152856814380635E+50*Z10_;
   Z16_=Z10_*Z16_;
   Z16_=1.010248376642855E+48 + 7.692307692307692E-2*Z16_;
   Z13_=Z16_*Z13_;
   Z13_=3.235822262175554E+43 + Z13_;
   Z13_=Z10_*Z13_;
   Z13_= - 4.013786172337223E+43 + Z13_;
   Z13_=Z10_*Z13_;
   Z13_= - 5.228637138997677E+42 + Z13_;
   Z13_=Z10_*Z13_;
   Z13_=7.373257772748786E+42 + Z13_;
   Z13_=Z10_*Z13_;
   Z13_=1.010068223206116E+41 + 9.090909090909090E-2*Z13_;
   Z13_=Z10_*Z13_;
   Z13_= - 1.641935539191188E+41 + Z13_;
   Z13_=Z10_*Z13_;
   Z13_= - 2.957007898981506E+40 + Z13_;
   Z13_=Z10_*Z13_;
   Z13_=8.020032667265607E+39 + 1.428571428571428E-1*Z13_;
   Z13_=Z10_*Z13_;
   Z13_=1.827199271599809E+39 + Z13_;
   Z13_=Z10_*Z13_;
   Z13_= - 4.072426281678678E+39 + Z13_;
   Z13_=Z10_*Z13_;
   Z13_= - 2.629598404293013E+38 + 2.E-1*Z13_;
   Z13_=Z10_*Z13_;
   Z13_=6.710715633826102E+38 + Z13_;
   Z13_=Z10_*Z13_;
   Z13_=3.963293989409197E+38 + Z13_;
   Z13_=Z10_*Z13_;
   Z13_= - 1.008174991051888E+39 + Z13_;
   Z13_=Z10_*Z13_;
   Z13_= - 2.154431329803065E+39 + Z13_;
   Z13_=Z10_*Z13_;
   Z13_= - 3.419517970721087E+39 + Z13_;
   Z13_=Z13_*Z14_;
   Z13_= - 1.253017328183604E+40 + Z13_;
   Z13_=Z10_*Z13_;
   Z13_=1.451966521729418E+40 + Z13_;
   Z12_=Z15_ + 4.166666666666666E-41*Z13_ + Z12_;
   Z13_=Z18_ - Z21_;
   Z14_= - Z3_*Z13_;
   Z13_=Z20_ - Z13_;
   Z13_=Z9_*Z13_;
   Z15_= - Z18_ + Z22_;
   Z15_=Z8_*Z15_;
   Z12_=Z19_ + Z26_ + Z25_ + Z24_ + Z13_ + Z14_ + 5.E-1*Z12_ + Z15_;
   Z12_=Z1_*Z12_;
   Z12_=4.292002547222587E-42 + Z12_;

    return 5.E-1*Z12_;

}

using namespace std;
using namespace resconst;
using namespace constants;

// overload int - complex operators to allow arithmetics
// std::complex<double> operator+( const int& lhs, const std::complex<double>& rhs)
// {
//   return double (lhs) + rhs;
// }
// std::complex<double> operator-( const int& lhs, const std::complex<double>& rhs)
// {
//   return double (lhs) - rhs;
// }
// std::complex<double> operator*( const int& lhs, const std::complex<double>& rhs)
// {
//   return double (lhs) * rhs;
// }
// std::complex<double> operator+( const std::complex<double>& lhs, const int& rhs)
// {
//   return lhs + double(rhs);
// }
// std::complex<double> operator-( const std::complex<double>& lhs, const int& rhs)
// {
//   return lhs - double(rhs);
// }
// std::complex<double> operator*( const std::complex<double>& lhs, const int& rhs )
// {
//   return lhs * double(rhs);
// }


double ccoeff::C1qq_delta;
double ccoeff::C2qq_delta;
double ccoeff::C3qq_delta;
double ccoeff::C3qq_delta_NFV;
double ccoeff::C4qq_delta;
double ccoeff::C4qq_delta_NFV;

double ccoeff::C1gg_delta;

complex <double> ccoeff::C1gg;
complex <double> *ccoeff::C1qg;
complex <double> *ccoeff::C1qq;
complex <double> *ccoeff::C1qgamma; // QED
complex <double> *ccoeff::C1qqqed; // QED
complex <double> *ccoeff::C1qqb;
complex <double> *ccoeff::C1qqp;
complex <double> *ccoeff::C1qqbp;
complex <double> *ccoeff::C1gq;

complex <double> *ccoeff::C2qg;
complex <double> *ccoeff::C2qq;
complex <double> *ccoeff::C2qqb;
complex <double> *ccoeff::C2qqp;
complex <double> *ccoeff::C2qqbp;
complex <double> *ccoeff::C2gq;

complex <double> *ccoeff::C3qg;
complex <double> *ccoeff::C3qq;
complex <double> *ccoeff::C3qqb;
complex <double> *ccoeff::C3qqp;
complex <double> *ccoeff::C3qqbp;

complex <double> *ccoeff::C4qg;
complex <double> *ccoeff::C4qq;
complex <double> *ccoeff::C4qqb;
complex <double> *ccoeff::C4qqp;
complex <double> *ccoeff::C4qqbp;

complex <double> *ccoeff::C1qg_1;
complex <double> *ccoeff::C1gq_1;
complex <double> *ccoeff::C1qgamma_1; // QED
complex <double> *ccoeff::C1qq_1;
complex <double> *ccoeff::C1qqqed_1; // QED
complex <double> *ccoeff::C1qqb_1;
complex <double> *ccoeff::C1qqp_1;
complex <double> *ccoeff::C1qqbp_1;
complex <double> *ccoeff::H2ggM_1;
complex <double> *ccoeff::H2gqM_1;
complex <double> *ccoeff::H2qqM_1;
complex <double> *ccoeff::G1N_1;
complex <double> *ccoeff::C2qg_1;
complex <double> *ccoeff::C2gq_1;
complex <double> *ccoeff::C2qq_1;
complex <double> *ccoeff::C2qqb_1;
complex <double> *ccoeff::C2qqp_1;
complex <double> *ccoeff::C2qqbp_1;
complex <double> *ccoeff::C3qg_1;
complex <double> *ccoeff::C3qq_1;
complex <double> *ccoeff::C3qqb_1;
complex <double> *ccoeff::C3qqp_1;
complex <double> *ccoeff::C3qqbp_1;
complex <double> *ccoeff::C4qg_1;
complex <double> *ccoeff::C4qq_1;
complex <double> *ccoeff::C4qqb_1;
complex <double> *ccoeff::C4qqp_1;
complex <double> *ccoeff::C4qqbp_1;

complex <double> *ccoeff::C1qg_2;
complex <double> *ccoeff::C1gq_2;
complex <double> *ccoeff::C1qgamma_2; // QED
complex <double> *ccoeff::C1qq_2;
complex <double> *ccoeff::C1qqqed_2; // QED
complex <double> *ccoeff::C1qqb_2;
complex <double> *ccoeff::C1qqp_2;
complex <double> *ccoeff::C1qqbp_2;
complex <double> *ccoeff::H2ggM_2;
complex <double> *ccoeff::H2gqM_2;
complex <double> *ccoeff::H2qqM_2;
complex <double> *ccoeff::G1N_2;
complex <double> *ccoeff::C2qg_2;
complex <double> *ccoeff::C2gq_2;
complex <double> *ccoeff::C2qq_2;
complex <double> *ccoeff::C2qqb_2;
complex <double> *ccoeff::C2qqp_2;
complex <double> *ccoeff::C2qqbp_2;
complex <double> *ccoeff::C3qg_2;
complex <double> *ccoeff::C3qq_2;
complex <double> *ccoeff::C3qqb_2;
complex <double> *ccoeff::C3qqp_2;
complex <double> *ccoeff::C3qqbp_2;
complex <double> *ccoeff::C4qg_2;
complex <double> *ccoeff::C4qq_2;
complex <double> *ccoeff::C4qqb_2;
complex <double> *ccoeff::C4qqp_2;
complex <double> *ccoeff::C4qqbp_2;

complex <double> *ccoeff::C1qg_reg;
complex <double> *ccoeff::C1qg_reg_1;
complex <double> *ccoeff::C1qg_reg_2;

//const int rule = 500;
//double t[rule];              //gauss nodes
//double fac[rule];            //overall factor
//complex <double> kern[rule]; //kernel of the Mellin transform

void ccoeff::allocate()
{
  if (opts.mellin1d)
    {
      C1qg = new complex <double>[mellinint::mdim*2];
      C1qgamma = new complex <double>[mellinint::mdim*2];
      C1qq = new complex <double>[mellinint::mdim*2];
      C1qqqed = new complex <double>[mellinint::mdim*2];
      C1qqb = new complex <double>[mellinint::mdim*2];
      C1qqp = new complex <double>[mellinint::mdim*2];
      C1qqbp = new complex <double>[mellinint::mdim*2];
      C1gq = new complex <double>[mellinint::mdim*2];

      C1qg_reg = new complex <double>[mellinint::mdim*2];
      
      C2qg = new complex <double>[mellinint::mdim*2];
      C2qq = new complex <double>[mellinint::mdim*2];
      C2qqb = new complex <double>[mellinint::mdim*2];
      C2qqp = new complex <double>[mellinint::mdim*2];
      C2qqbp = new complex <double>[mellinint::mdim*2];
      C2gq = new complex <double>[mellinint::mdim*2];

      C3qg = new complex <double>[mellinint::mdim*2];
      C3qq = new complex <double>[mellinint::mdim*2];
      C3qqb = new complex <double>[mellinint::mdim*2];
      C3qqp = new complex <double>[mellinint::mdim*2];
      C3qqbp = new complex <double>[mellinint::mdim*2];

      C4qg = new complex <double>[mellinint::mdim*2];
      C4qq = new complex <double>[mellinint::mdim*2];
      C4qqb = new complex <double>[mellinint::mdim*2];
      C4qqp = new complex <double>[mellinint::mdim*2];
      C4qqbp = new complex <double>[mellinint::mdim*2];
      
      //fill(C1qg,   C1qg+2*mellinint::mdim, 0.);
      //fill(C1qq,   C1qq+2*mellinint::mdim, 0.);
      //fill(C1qqb,  C1qqb+2*mellinint::mdim, 0.);
      //fill(C1qqp,  C1qqp+2*mellinint::mdim, 0.);
      //fill(C1qqbp, C1qqbp+2*mellinint::mdim, 0.);
      //
      //fill(C2qg,   C2qg+2*mellinint::mdim, 0.);
      //fill(C2qq,   C2qq+2*mellinint::mdim, 0.);
      //fill(C2qqb,  C2qqb+2*mellinint::mdim, 0.);
      //fill(C2qqp,  C2qqp+2*mellinint::mdim, 0.);
      //fill(C2qqbp, C2qqbp+2*mellinint::mdim, 0.);
      //
      //fill(C3qg,   C2qg+2*mellinint::mdim, 0.);
      //fill(C3qq,   C2qq+2*mellinint::mdim, 0.);
      //fill(C3qqb,  C2qqb+2*mellinint::mdim, 0.);
      //fill(C3qqp,  C2qqp+2*mellinint::mdim, 0.);
      //fill(C3qqbp, C2qqbp+2*mellinint::mdim, 0.);
    }
  else
    {
      C1qg_1 = new complex <double>[mellinint::mdim*2];
      C1gq_1 = new complex <double>[mellinint::mdim*2];
      C1qgamma_1 = new complex <double>[mellinint::mdim*2]; // QED
      C1qq_1 = new complex <double>[mellinint::mdim*2];
      C1qqqed_1 = new complex <double>[mellinint::mdim*2]; // QED
      C1qqb_1 = new complex <double>[mellinint::mdim*2];
      C1qqp_1 = new complex <double>[mellinint::mdim*2];
      C1qqbp_1 = new complex <double>[mellinint::mdim*2];
      H2ggM_1 = new complex <double>[mellinint::mdim*2];
      H2gqM_1 = new complex <double>[mellinint::mdim*2];
      H2qqM_1 = new complex <double>[mellinint::mdim*2];
      G1N_1 = new complex <double>[mellinint::mdim*2];
      C2qg_1 = new complex <double>[mellinint::mdim*2];
      C2gq_1 = new complex <double>[mellinint::mdim*2];
      C2qq_1 = new complex <double>[mellinint::mdim*2];
      C2qqb_1 = new complex <double>[mellinint::mdim*2];
      C2qqp_1 = new complex <double>[mellinint::mdim*2];
      C2qqbp_1 = new complex <double>[mellinint::mdim*2];
      C3qg_1 = new complex <double>[mellinint::mdim*2];
      C3qq_1 = new complex <double>[mellinint::mdim*2];
      C3qqb_1 = new complex <double>[mellinint::mdim*2];
      C3qqp_1 = new complex <double>[mellinint::mdim*2];
      C3qqbp_1 = new complex <double>[mellinint::mdim*2];
      C4qg_1 = new complex <double>[mellinint::mdim*2];
      C4qq_1 = new complex <double>[mellinint::mdim*2];
      C4qqb_1 = new complex <double>[mellinint::mdim*2];
      C4qqp_1 = new complex <double>[mellinint::mdim*2];
      C4qqbp_1 = new complex <double>[mellinint::mdim*2];
      C1qg_2 = new complex <double>[mellinint::mdim*2];
      C1gq_2 = new complex <double>[mellinint::mdim*2];
      C1qgamma_2 = new complex <double>[mellinint::mdim*2];
      C1qq_2 = new complex <double>[mellinint::mdim*2];
      C1qqqed_2 = new complex <double>[mellinint::mdim*2]; //QED
      C1qqb_2 = new complex <double>[mellinint::mdim*2];
      C1qqp_2 = new complex <double>[mellinint::mdim*2];
      C1qqbp_2 = new complex <double>[mellinint::mdim*2];
      H2ggM_2 = new complex <double>[mellinint::mdim*2];
      H2gqM_2 = new complex <double>[mellinint::mdim*2];
      H2qqM_2 = new complex <double>[mellinint::mdim*2];
      G1N_2 = new complex <double>[mellinint::mdim*2];
      C2qg_2 = new complex <double>[mellinint::mdim*2];
      C2gq_2 = new complex <double>[mellinint::mdim*2];
      C2qq_2 = new complex <double>[mellinint::mdim*2];
      C2qqb_2 = new complex <double>[mellinint::mdim*2];
      C2qqp_2 = new complex <double>[mellinint::mdim*2];
      C2qqbp_2 = new complex <double>[mellinint::mdim*2];
      C3qg_2 = new complex <double>[mellinint::mdim*2];
      C3qq_2 = new complex <double>[mellinint::mdim*2];
      C3qqb_2 = new complex <double>[mellinint::mdim*2];
      C3qqp_2 = new complex <double>[mellinint::mdim*2];
      C3qqbp_2 = new complex <double>[mellinint::mdim*2];
      C4qg_2 = new complex <double>[mellinint::mdim*2];
      C4qq_2 = new complex <double>[mellinint::mdim*2];
      C4qqb_2 = new complex <double>[mellinint::mdim*2];
      C4qqp_2 = new complex <double>[mellinint::mdim*2];
      C4qqbp_2 = new complex <double>[mellinint::mdim*2];
      C1qg_reg_1 = new complex <double>[mellinint::mdim*2];
      C1qg_reg_2 = new complex <double>[mellinint::mdim*2];
    }
}

void ccoeff::init()
{
  //Evaluate N-independent part
  delta();
  icoeff::allocate();
  icoeff::init();

  //  if (opts.melup <= 1)
  //    allocate();
  //if (opts.melup == 0)
  //  if (opts.mellin1d)
  //    calc1d();
  //  else
  //    calc2d();
}

void ccoeff::delta()
{
  //delta pieces -> N-independent part (proportional to delta(1-z))
  int NF2 = NF*NF; //Number of flavours from resconst

  //NFV is the charge weighted sum of the quark flavours (in the case of purely electromagnetic interactions we have NFgamma=1/eq*Sum_{i=,nf} ei).
  double NFV;
  NFV = 0.;

  //Soft function coefficients non-exponentiated (From Appendix 4 of https://arxiv.org/pdf/1912.05778.pdf)
  double s1 = - 1./2.*CF*zeta2;
  double s2 = - 41./162.*CF*NF + 7./36.*CF*NF*zeta3 + 5./24.*CF*NF*zeta2 + 607./324.*CA*CF + 5./8.*CA*CF*zeta4 - 77./72.*CA*CF*zeta3 - 67./48.*CA*CF*zeta2;
  double s3 = - 4./6561.*CF*pow(NF,2) - 11./432.*CF*pow(NF,2)*zeta4 - 35./972.*CF*pow(NF,2)*zeta3 - 17./216.*CF*pow(NF,2)*zeta2 - 42727./31104.*pow(CF,2)*NF + 7./18.*pow(CF,2)*NF*zeta5 + 19./72.*pow(CF,2)*NF*zeta4 + 109./162.*pow(CF,2)*NF*zeta3 + 275./576.*pow(CF,2)*NF*zeta2 - 5./12.*pow(CF,2)*NF*zeta2*zeta3 - 412765./419904.*CA*CF*NF - 23./24.*CA*CF*NF*zeta5 - 13./54.*CA*CF*NF*zeta4 + 1019./648.*CA*CF*NF*zeta3 + 37265./23328.*CA*CF*NF*zeta2 + 5./72.*CA*CF*NF*zeta2*zeta3 + 5211949./839808.*pow(CA,2)*CF + 451./144.*pow(CA,2)*CF*zeta5 + 3649./1728.*pow(CA,2)*CF*zeta4 - 37783./3888.*pow(CA,2)*CF*zeta3 + 29./18.*pow(CA,2)*CF*pow(zeta3,2) - 297481./46656.*pow(CA,2)*CF*zeta2 - 1543./1512.*pow(CA,2)*CF*zeta2*zeta4 + 275./144.*pow(CA,2)*CF*zeta2*zeta3;
  //double s4 = 0.;
  double s4 = (1.+opts.H4)*accel::levin(s1,s2,s3);
 
  double SS1 = s1;
  double SS2 = s2+s1*s1/2.;
  double SS3 = s3+s1*s2+pow(s1,3)/6.;
  double SS4 = s4+s1*s3+pow(s2,2)/2.+pow(s1,2)*s2/2.+pow(s1,4)/24.; //--> check this formula!
  
  double S1 = -CF*pi2/12.;
  double S2 = CF*(9.*CF*pi4 + CA*(4856. - 603.*pi2 + 18.*pi4 - 2772.*zeta3) + NF*(-656. + 90.*pi2 + 504.*zeta3))/2592.;
  double S3 = CF*(pow(CA,2)*(689661*pi4 - 55548*pi6 + 105*pi2*(-297481 + 89100*zeta3) + 35*(5211949 - 8161128*zeta3 + 1353024*pow(zeta3,2) + 2630232*zeta5))
		  - 7*(405*pow(CF,2)*pi6 + 4*NF2*(640 + 13770*pi2 + 297*pi4 + 37800*zeta3) - 27*CF*NF*(6*pi4 - 5*pi2*(-3131 + 2664*zeta3) + 5*(-42727 + 20928*zeta3 + 12096*zeta5)))
		  - 7*CA*(135*CF*pi2*(4856 - 603*pi2 + 18*pi4 - 2772*zeta3) + 2*NF*(5616*pi4 - 75*pi2*(7453 + 324*zeta3) + 5*(412765 - 660312*zeta3 + 402408*zeta5))))/29393280.;
  //double S4 = accel::levin(S1,S2,S3);
  double S4 = SS4;

  double K1 = CF*(-4 + 7*pi2/12.);
  double K2 = CF*(27.*CF*(7665. - 1660.*pi2 + 134.*pi4 - 3600.*zeta3) + 10.*NF*(4085. - 546.*pi2 + 72.*zeta3) + CA*(-255785. + 31830.*pi2 - 288.*pi4 + 112680.*zeta3))/25920.;
  double K3 = CF*(NF2*(-190931./419904. + (403*pi2)/3888. + (43*pi4)/38880. - (13*zeta3)/243.)
		  + CA*CF*(824281./20736. - (406507*pi2)/62208. + (92237*pi4)/155520. - (739*pi6)/54432. - (13141*zeta3)/432. + (845*pi2*zeta3)/288. + (37*pow(zeta3,2))/12. - (689*zeta5)/72.)
		  + pow(CA,2)*(-51082685./1679616. + (596513*pi2)/139968. - (4303*pi4)/311040. + (299*pi6)/102060. + (505087*zeta3)/15552. - (73*pi2*zeta3)/36. - (71*pow(zeta3,2))/18. - (217*zeta5)/144.)
		  + CF*NF*(-56963./31104. + (13705*pi2)/15552. - (1463*pi4)/15552. +  (815*zeta3)/162. - (37*pi2*zeta3)/144. - (13*zeta5)/9.)
		  + ((-4 + pow(NC,2))*NFV*(1./8. + (5*pi2)/96. - pi4/2880. + (7*zeta3)/48. - (5*zeta5)/6.))/NC
		  + CA*NF*(1700171./209952. - (201749*pi2)/139968. -  (35*pi4)/15552. - (134*zeta3)/27. + (37*pi2*zeta3)/144. -  zeta5/24.)
		  + pow(CF,2)*(-5599/384. + (4339*pi2)/2304. - (173*pi4)/480. +  (27403*pi6)/1088640. - (115*zeta3)/16. - (35*pi2*zeta3)/48. + pow(zeta3,2)/2. + (83*zeta5)/4.));
  double K3_NFV = ( (-4 + pow(NC,2)) *(1./8. + (5*pi2)/96. - pi4/2880. + (7*zeta3)/48. - (5*zeta5)/6.) ) /NC; 
  //double K4 = accel::levin(K1,K2,K3);

  double Hq4SCET = 321482899./26873856.+41601./(32768.*pow(NC,4))+21812995./(2654208.*pow(NC,2))+(583291805.*pow(NC,2))/71663616.-(6361535363.*pow(NC,4))/214990848.-(5528861.*NF)/(3981312.*pow(NC,3))-(677109931.*NF)/(214990848.*NC)-(3301359761.*NC*NF)/214990848.+(237612677.*pow(NC,3)*NF)/11943936.+(818632639.*pow(NF,2))/214990848.-(8087615.*pow(NF,2))/(26873856.*pow(NC,2))-(9307799.*pow(NC,2)*pow(NF,2))/2654208.-(1865531.*pow(NF,3))/(13436928.*NC)+(1865531.*NC*pow(NF,3))/13436928.+(73.*NFV)/(72.*pow(NC,3))+(617.*NFV)/(288.*NC)-(2309.*NC*NFV)/576.+(491.*pow(NC,3)*NFV)/576.+(245.*NF*NFV)/288.-(49.*NF*NFV)/(72.*pow(NC,2))-49./288.*pow(NC,2)*NF*NFV-(282994049.*pi2)/80621568.+(1411.*pi2)/(49152.*pow(NC,4))-(190811.*pi2)/(497664.*pow(NC,2))-(151186117.*pow(NC,2)*pi2)/161243136.+(1548736387.*pow(NC,4)*pi2)/322486272.+(2613235.*NF*pi2)/(5971968.*pow(NC,3))+(255090349.*NF*pi2)/(322486272.*NC)+(373605629.*NC*NF*pi2)/161243136.-(1143416297.*pow(NC,3)*NF*pi2)/322486272.-(112557385.*pow(NF,2)*pi2)/161243136.-(1064527.*pow(NF,2)*pi2)/(40310784.*pow(NC,2))+(116815493.*pow(NC,2)*pow(NF,2)*pi2)/161243136.+(19847.*pow(NF,3)*pi2)/(559872.*NC)-(19847.*NC*pow(NF,3)*pi2)/559872.+(25.*NFV*pi2)/(96.*pow(NC,3))+(2893.*NFV*pi2)/(3456.*NC)-(9595.*NC*NFV*pi2)/6912.+(2009.*pow(NC,3)*NFV*pi2)/6912.+(695.*NF*NFV*pi2)/3456.-(139.*NF*NFV*pi2)/(864.*pow(NC,2))-(139.*pow(NC,2)*NF*NFV*pi2)/3456.-(5021281.*pi4)/107495424.+(77933.*pi4)/(4423680.*pow(NC,4))+(656111.*pi4)/(4976640.*pow(NC,2))-(15813319.*pow(NC,2)*pi4)/26873856.+(104417575.*pow(NC,4)*pi4)/214990848.-(457717.*NF*pi4)/(9953280.*pow(NC,3))-(14829269.*NF*pi4)/(268738560.*NC)+(91615573.*NC*NF*pi4)/268738560.-(12885589.*pow(NC,3)*NF*pi4)/53747712.-(29591.*pow(NF,2)*pi4)/746496.+(4783.*pow(NF,2)*pi4)/(466560.*pow(NC,2))+(109691.*pow(NC,2)*pow(NF,2)*pi4)/3732480.+(949.*pow(NF,3)*pi4)/(933120.*NC)-(949.*NC*pow(NF,3)*pi4)/933120.-(41.*NFV*pi4)/(1920.*pow(NC,3))-(2.*NFV*pi4)/(405.*NC)+(28367.*NC*NFV*pi4)/829440.-(6559.*pow(NC,3)*NFV*pi4)/829440.-(95.*NF*NFV*pi4)/20736.+(19.*NF*NFV*pi4)/(5184.*pow(NC,2))+(19.*pow(NC,2)*NF*NFV*pi4)/20736.+(1099393.*pi6)/29859840.-(22679.*pi6)/(8709120.*pow(NC,4))-(722717.*pi6)/(59719680.*pow(NC,2))-(11045201.*pow(NC,2)*pi6)/418037760.+(25733.*pow(NC,4)*pi6)/5971968.+(90193.*NF*pi6)/(52254720.*pow(NC,3))-(75449.*NF*pi6)/(11612160.*NC)+(253727.*NC*NF*pi6)/52254720.-(419.*pow(NC,3)*NF*pi6)/4976640.-(1751.*pow(NF,2)*pi6)/26127360.+(1217.*pow(NF,2)*pi6)/(8709120.*pow(NC,2))-(95.*pow(NC,2)*pow(NF,2)*pi6)/1306368.+(943.*NFV*pi6)/(725760.*pow(NC,3))-(41449.*NFV*pi6)/(5225472.*NC)+(428999.*NC*NFV*pi6)/52254720.-(16481.*pow(NC,3)*NFV*pi6)/10450944.-(11.*NF*NFV*pi6)/11664.+(11.*NF*NFV*pi6)/(14580.*pow(NC,2))+(11.*pow(NC,2)*NF*NFV*pi6)/58320.+(6530779.*pi8)/6270566400.+(5926681.*pi8)/(20901888000.*pow(NC,4))-(743663.*pi8)/(870912000.*pow(NC,2))-(19452931.*pow(NC,2)*pi8)/31352832000.+(267479.*pow(NC,4)*pi8)/1791590400.+23./256.*zeta53-(69.*zeta53)/(320.*pow(NC,4))+(1277.*zeta53)/(1280.*pow(NC,2))-61./80.*pow(NC,2)*zeta53-7./64.*pow(NC,4)*zeta53-(211844861.*zeta3)/2985984.-(23965.*zeta3)/(6144.*pow(NC,4))+(65653.*zeta3)/(6144.*pow(NC,2))+(11207041.*pow(NC,2)*zeta3)/373248.+(33976055.*pow(NC,4)*zeta3)/995328.-(820819.*NF*zeta3)/(331776.*pow(NC,3))+(3797689.*NF*zeta3)/(186624.*NC)-(17374699.*NC*NF*zeta3)/2985984.-(2000053.*pow(NC,3)*NF*zeta3)/165888.+(35681.*pow(NF,2)*zeta3)/46656.-(1005107.*pow(NF,2)*zeta3)/(746496.*pow(NC,2))+(144737.*pow(NC,2)*pow(NF,2)*zeta3)/248832.-(13.*pow(NF,3)*zeta3)/(7776.*NC)+(13.*NC*pow(NF,3)*zeta3)/7776.+(19.*NFV*zeta3)/(16.*pow(NC,3))+(529.*NFV*zeta3)/(96.*NC)-6487./768.*NC*NFV*zeta3+1343./768.*pow(NC,3)*NFV*zeta3+35./64.*NF*NFV*zeta3-(7.*NF*NFV*zeta3)/(16.*pow(NC,2))-7./64.*pow(NC,2)*NF*NFV*zeta3+(1177831.*pi2*zeta3)/746496.-(3895.*pi2*zeta3)/(6144.*pow(NC,4))+(1058855.*pi2*zeta3)/(331776.*pow(NC,2))-(1343051.*pow(NC,2)*pi2*zeta3)/2985984.-(5502499.*pow(NC,4)*pi2*zeta3)/1492992.-(132493.*NF*pi2*zeta3)/(497664.*pow(NC,3))-(567893.*NF*pi2*zeta3)/(497664.*NC)-(51719.*NC*NF*pi2*zeta3)/497664.+(752105.*pow(NC,3)*NF*pi2*zeta3)/497664.+(19.*pow(NF,2)*pi2*zeta3)/1458.+(21509.*pow(NF,2)*pi2*zeta3)/(186624.*pow(NC,2))-(23941.*pow(NC,2)*pow(NF,2)*pi2*zeta3)/186624.-(pow(NF,3)*pi2*zeta3)/(3456.*NC)+(NC*pow(NF,3)*pi2*zeta3)/3456.+(19.*NFV*pi2*zeta3)/(288.*pow(NC,3))+(91.*NFV*pi2*zeta3)/(144.*NC)-(4039.*NC*NFV*pi2*zeta3)/4608.+(823.*pow(NC,3)*NFV*pi2*zeta3)/4608.+35./576.*NF*NFV*pi2*zeta3-(7.*NF*NFV*pi2*zeta3)/(144.*pow(NC,2))-7./576.*pow(NC,2)*NF*NFV*pi2*zeta3+(82177.*pi4*zeta3)/829440.+(863.*pi4*zeta3)/(46080.*pow(NC,4))-(61981.*pi4*zeta3)/(414720.*pow(NC,2))+71./576.*pow(NC,2)*pi4*zeta3-(75989.*pow(NC,4)*pi4*zeta3)/829440.+(6083.*NF*pi4*zeta3)/(207360.*pow(NC,3))-(12167.*NF*pi4*zeta3)/(414720.*NC)-(1403.*NC*NF*pi4*zeta3)/138240.+(421.*pow(NC,3)*NF*pi4*zeta3)/41472.-(11.*NFV*pi4*zeta3)/(1440.*pow(NC,3))+(5.*NFV*pi4*zeta3)/(576.*NC)-(19.*NC*NFV*pi4*zeta3)/23040.-(pow(NC,3)*NFV*pi4*zeta3)/4608.+(746603.*pow(zeta3,2))/82944.+(2697.*pow(zeta3,2))/(1024.*pow(NC,4))-(677.*pow(zeta3,2))/(1536.*pow(NC,2))-(70517.*pow(NC,2)*pow(zeta3,2))/41472.-(196867.*pow(NC,4)*pow(zeta3,2))/20736.+(113.*NF*pow(zeta3,2))/(128.*pow(NC,3))-(26489.*NF*pow(zeta3,2))/(10368.*NC)+(9517.*NC*NF*pow(zeta3,2))/41472.+(59827.*pow(NC,3)*NF*pow(zeta3,2))/41472.-(673.*pow(NF,2)*pow(zeta3,2))/5184.+(2279.*pow(NF,2)*pow(zeta3,2))/(20736.*pow(NC,2))+(413.*pow(NC,2)*pow(NF,2)*pow(zeta3,2))/20736.+(21.*NFV*pow(zeta3,2))/(32.*pow(NC,3))-(2755.*NFV*pow(zeta3,2))/(384.*NC)+3113./384.*NC*NFV*pow(zeta3,2)-305./192.*pow(NC,3)*NFV*pow(zeta3,2)-55./48.*NF*NFV*pow(zeta3,2)+(11.*NF*NFV*pow(zeta3,2))/(12.*pow(NC,2))+11./48.*pow(NC,2)*NF*NFV*pow(zeta3,2)+(1289.*pi2*pow(zeta3,2))/6912.+(131.*pi2*pow(zeta3,2))/(1536.*pow(NC,4))-(25.*pi2*pow(zeta3,2))/(288.*pow(NC,2))+(1387.*pow(NC,2)*pi2*pow(zeta3,2))/13824.-(493.*pow(NC,4)*pi2*pow(zeta3,2))/1728.+(1092277.*zeta5)/18432.+(29041.*zeta5)/(6144.*pow(NC,4))-(1110893.*zeta5)/(36864.*pow(NC,2))-(15801697.*pow(NC,2)*zeta5)/552960.-(729227.*pow(NC,4)*zeta5)/138240.+(3461.*NF*zeta5)/(18432.*pow(NC,3))-(235601.*NF*zeta5)/(55296.*NC)+(2457541.*NC*NF*zeta5)/276480.-(49313.*pow(NC,3)*NF*zeta5)/10240.-(20927.*pow(NF,2)*zeta5)/17280.+(949.*pow(NF,2)*zeta5)/(3456.*pow(NC,2))+899./960.*pow(NC,2)*pow(NF,2)*zeta5+(53.*pow(NF,3)*zeta5)/(17280.*NC)-(53.*NC*pow(NF,3)*zeta5)/17280.-(955.*NFV*zeta5)/(192.*pow(NC,3))+(19915.*NFV*zeta5)/(2304.*NC)-(19705.*NC*NFV*zeta5)/4608.+(2795.*pow(NC,3)*NFV*zeta5)/4608.-275./288.*NF*NFV*zeta5+(55.*NF*NFV*zeta5)/(72.*pow(NC,2))+55./288.*pow(NC,2)*NF*NFV*zeta5+35./432.*pi2*zeta5+(371.*pi2*zeta5)/(768.*pow(NC,4))-(37715.*pi2*zeta5)/(55296.*pow(NC,2))-(4821.*pow(NC,2)*pi2*zeta5)/2048.+(68345.*pow(NC,4)*pi2*zeta5)/27648.+(517.*NF*pi2*zeta5)/(13824.*pow(NC,3))-(NF*pi2*zeta5)/(9.*NC)+(857.*NC*NF*pi2*zeta5)/1728.-(5837.*pow(NC,3)*NF*pi2*zeta5)/13824.+(55.*NFV*pi2*zeta5)/(144.*pow(NC,3))-(35.*NFV*pi2*zeta5)/(24.*NC)+1015./768.*NC*NFV*pi2*zeta5-(565.*pow(NC,3)*NFV*pi2*zeta5)/2304.+149./256.*zeta3*zeta5-(7.*zeta3*zeta5)/(8.*pow(NC,4))-(7.*zeta3*zeta5)/(64.*pow(NC,2))-1429./256.*pow(NC,2)*zeta3*zeta5+383./64.*pow(NC,4)*zeta3*zeta5-(40251.*zeta7)/2048.-(29919.*zeta7)/(4096.*pow(NC,4))+(243851.*zeta7)/(24576.*pow(NC,2))+(1329235.*pow(NC,2)*zeta7)/73728.-(36605.*pow(NC,4)*zeta7)/36864.+(1707.*NF*zeta7)/(2048.*pow(NC,3))-(35.*NF*zeta7)/(512.*NC)-(19265.*NC*NF*zeta7)/9216.+(24427.*pow(NC,3)*NF*zeta7)/18432.-(721.*NFV*zeta7)/(192.*pow(NC,3))-(497.*NFV*zeta7)/(768.*NC)+(8813.*NC*NFV*zeta7)/1536.-(2051.*pow(NC,3)*NFV*zeta7)/1536.;
  double Hq3SCET=-((51082685.*pow(CA,2)*CF)/1679616.)+(824281.*CA*pow(CF,2))/20736.-(5599.*pow(CF,3))/384.+(1700171.*CA*CF*NF)/209952.-(56963.*pow(CF,2)*NF)/31104.-(190931.*CF*pow(NF,2))/419904.-(CF*NFV)/(2.*NC)+(CF*NFV*pow(NC,2))/(8.*NC)+(522155.*pow(CA,2)*CF*pi2)/139968.-(406507.*CA*pow(CF,2)*pi2)/62208.+(4339.*pow(CF,3)*pi2)/2304.-(179879.*CA*CF*NF*pi2)/139968.+(15163.*pow(CF,2)*NF*pi2)/15552.+(403.*CF*pow(NF,2)*pi2)/3888.+1./32.*CF*(17.*pow(CA,2)-5.*CA*NF-3.*CF*NF)*pi2-(5.*CF*NFV*pi2)/(24.*NC)+(5.*CF*NFV*pow(NC,2)*pi2)/(96.*NC)-(4303.*pow(CA,2)*CF*pi4)/311040.+(92237.*CA*pow(CF,2)*pi4)/155520.-(173.*pow(CF,3)*pi4)/480.-(35.*CA*CF*NF*pi4)/15552.-(1463.*pow(CF,2)*NF*pi4)/15552.+(43.*CF*pow(NF,2)*pi4)/38880.+(CF*NFV*pi4)/(720.*NC)-(CF*NFV*pow(NC,2)*pi4)/(2880.*NC)+(299.*pow(CA,2)*CF*pi6)/102060.-(739.*CA*pow(CF,2)*pi6)/54432.+(27403.*pow(CF,3)*pi6)/1088640.+11./48.*pow(CA,2)*CF*pi2*zeta3-7./24.*CA*CF*NF*pi2*zeta3+1./4.*pow(CF,2)*NF*pi2*zeta3+(505087.*pow(CA,2)*CF*zeta3)/15552.-13141./432.*CA*pow(CF,2)*zeta3-115./16.*pow(CF,3)*zeta3-134./27.*CA*CF*NF*zeta3+815./162.*pow(CF,2)*NF*zeta3-13./243.*CF*pow(NF,2)*zeta3-(7.*CF*NFV*zeta3)/(12.*NC)+(7.*CF*NFV*pow(NC,2)*zeta3)/(48.*NC)-325./144.*pow(CA,2)*CF*pi2*zeta3+845./288.*CA*pow(CF,2)*pi2*zeta3-35./48.*pow(CF,3)*pi2*zeta3+79./144.*CA*CF*NF*pi2*zeta3-73./144.*pow(CF,2)*NF*pi2*zeta3-71./18.*pow(CA,2)*CF*pow(zeta3,2)+37./12.*CA*pow(CF,2)*pow(zeta3,2)+1./2.*pow(CF,3)*pow(zeta3,2)-217./144.*pow(CA,2)*CF*zeta5-689./72.*CA*pow(CF,2)*zeta5+83./4.*pow(CF,3)*zeta5-1./24.*CA*CF*NF*zeta5-13./9.*pow(CF,2)*NF*zeta5+(10.*CF*NFV*zeta5)/(3.*NC)-(5.*CF*NFV*pow(NC,2)*zeta5)/(6.*NC);
  double Hq2SCET=-((51157.*CA*CF)/5184.)+(511.*pow(CF,2))/64.+(4085.*CF*NF)/2592.+1061./864.*CA*CF*pi2-(83.*pow(CF,2)*pi2)/48.-91./432.*CF*NF*pi2-1./90.*CA*CF*pi4+(67.*pow(CF,2)*pi4)/480.+313./72.*CA*CF*zeta3-15./4.*pow(CF,2)*zeta3+1./36.*CF*NF*zeta3;
  double Hq1SCET=-4.*CF+(7.*CF*pi2)/12.;

  double K4 = Hq4SCET;
  NFV = 1.;
  double K4_NFV = (73.*NFV)/(72.*pow(NC,3))+(617.*NFV)/(288.*NC)-(2309.*NC*NFV)/576.+(491.*pow(NC,3)*NFV)/576.+(245.*NF*NFV)/288.-(49.*NF*NFV)/(72.*pow(NC,2))-49./288.*pow(NC,2)*NF*NFV+(25.*NFV*pi2)/(96.*pow(NC,3))+(2893.*NFV*pi2)/(3456.*NC)-(9595.*NC*NFV*pi2)/6912.+(2009.*pow(NC,3)*NFV*pi2)/6912.+(695.*NF*NFV*pi2)/3456.-(139.*NF*NFV*pi2)/(864.*pow(NC,2))-(139.*pow(NC,2)*NF*NFV*pi2)/3456.-(41.*NFV*pi4)/(1920.*pow(NC,3))-(2.*NFV*pi4)/(405.*NC)+(28367.*NC*NFV*pi4)/829440.-(6559.*pow(NC,3)*NFV*pi4)/829440.-(95.*NF*NFV*pi4)/20736.+(19.*NF*NFV*pi4)/(5184.*pow(NC,2))+(19.*pow(NC,2)*NF*NFV*pi4)/20736.+(943.*NFV*pi6)/(725760.*pow(NC,3))-(41449.*NFV*pi6)/(5225472.*NC)+(428999.*NC*NFV*pi6)/52254720.-(16481.*pow(NC,3)*NFV*pi6)/10450944.-(11.*NF*NFV*pi6)/11664.+(11.*NF*NFV*pi6)/(14580.*pow(NC,2))+(11.*pow(NC,2)*NF*NFV*pi6)/58320.+(19.*NFV*zeta3)/(16.*pow(NC,3))+(529.*NFV*zeta3)/(96.*NC)-6487./768.*NC*NFV*zeta3+1343./768.*pow(NC,3)*NFV*zeta3+35./64.*NF*NFV*zeta3-(7.*NF*NFV*zeta3)/(16.*pow(NC,2))-7./64.*pow(NC,2)*NF*NFV*zeta3+(19.*NFV*pi2*zeta3)/(288.*pow(NC,3))+(91.*NFV*pi2*zeta3)/(144.*NC)-(4039.*NC*NFV*pi2*zeta3)/4608.+(823.*pow(NC,3)*NFV*pi2*zeta3)/4608.+35./576.*NF*NFV*pi2*zeta3-(7.*NF*NFV*pi2*zeta3)/(144.*pow(NC,2))-7./576.*pow(NC,2)*NF*NFV*pi2*zeta3-(11.*NFV*pi4*zeta3)/(1440.*pow(NC,3))+(5.*NFV*pi4*zeta3)/(576.*NC)-(19.*NC*NFV*pi4*zeta3)/23040.-(pow(NC,3)*NFV*pi4*zeta3)/4608.+(21.*NFV*pow(zeta3,2))/(32.*pow(NC,3))-(2755.*NFV*pow(zeta3,2))/(384.*NC)+3113./384.*NC*NFV*pow(zeta3,2)-305./192.*pow(NC,3)*NFV*pow(zeta3,2)-55./48.*NF*NFV*pow(zeta3,2)+(11.*NF*NFV*pow(zeta3,2))/(12.*pow(NC,2))+11./48.*pow(NC,2)*NF*NFV*pow(zeta3,2)-(955.*NFV*zeta5)/(192.*pow(NC,3))+(19915.*NFV*zeta5)/(2304.*NC)-(19705.*NC*NFV*zeta5)/4608.+(2795.*pow(NC,3)*NFV*zeta5)/4608.-275./288.*NF*NFV*zeta5+(55.*NF*NFV*zeta5)/(72.*pow(NC,2))+55./288.*pow(NC,2)*NF*NFV*zeta5+(55.*NFV*pi2*zeta5)/(144.*pow(NC,3))-(35.*NFV*pi2*zeta5)/(24.*NC)+1015./768.*NC*NFV*pi2*zeta5-(565.*pow(NC,3)*NFV*pi2*zeta5)/2304.-(721.*NFV*zeta7)/(192.*pow(NC,3))-(497.*NFV*zeta7)/(768.*NC)+(8813.*NC*NFV*zeta7)/1536.-(2051.*pow(NC,3)*NFV*zeta7)/1536.;

  C1qq_delta = (K1+S1)/2.;
  C2qq_delta = (K2+S2+K1*S1)/2. - pow(C1qq_delta,2)/2.;
  C3qq_delta = (K3 + S3 + K2*S1 + K1*S2)/2. - C1qq_delta*C2qq_delta;
  C4qq_delta = (K4 + S4 + K2*S2 + K1*S3 + K3*S1)/2. - C1qq_delta*C3qq_delta - pow(C2qq_delta,2)/2.;

  //C4qq_delta += opts.H4*(K4 + S4 + K2*S2 + K1*S3 + K3*S1)/2.;
  
  //cout << "Levin K4 " << accel::levin(K1,K2,K3) << endl;
  //cout << "K1 " << K1 << endl;
  //cout << "K2 " << K2 << endl;
  //cout << "K3 " << K3 << endl;
  //cout << "K4 " << K4 << endl;
  //cout << endl;
  //cout << "s1 " << s1 << endl;
  //cout << "s2 " << s2 << endl;
  //cout << "s3 " << s3 << endl;
  //cout << "s4 " << s4 << endl;
  //cout << endl;
  //cout << endl;
  //cout << "S1 " << S1 << endl;
  //cout << "S2 " << S2 << endl;
  //cout << "S3 " << S3 << endl;
  //cout << "S4 " << S4 << endl;
  //cout << endl;
  //cout << "SS1 " << SS1 << endl;
  //cout << "SS2 " << SS2 << endl;
  //cout << "SS3 " << SS3 << endl;
  //cout << "SS4 " << SS4 << endl;
  //cout << endl;
  //cout << "H1qqdelta " << (K1+S1)			    << endl;
  //cout << "H2qqdelta " << (K2+S2+K1*S1)                     << endl;
  //cout << "H3qqdelta " << (K3 + S3 + K2*S1 + K1*S2)         << endl;
  //cout << "H4qqdelta " << (K4 + S4 + K2*S2 + K1*S3 + K3*S1) << endl;
  //cout << endl;
  //cout << "H1qqdelta S1=0 " << (K1)			    << endl;
  //cout << "H2qqdelta S2=0 " << (K2+S2-s2+K1*S1)                    << endl;
  //cout << "H3qqdelta S3=0 " << (K3+S3-s3 + K2*S1 + K1*S2)          << endl;
  //cout << "H4qqdelta S4=0 " << (K4+S4-s4 + K2*S2 + K1*S3 + K3*S1)  << endl;
  //cout << endl;
  //cout << "C1qq_delta " << C1qq_delta << endl;
  //cout << "C2qq_delta " << C2qq_delta << endl;
  //cout << "C3qq_delta " << C3qq_delta << endl;
  //cout << "C4qq_delta " << C4qq_delta << endl;
  
  //Naive CA/CF scaling
  //double S1g = -CA*pi2/12.;
  //double K1g = CA*(-4 + 7*pi2/12.);
  //C1gg_delta = (K1g+S1g)/2.;

  C1gg_delta = (pi2/2.+11./2.+pi2)/2.; //from HRes
  
  if (opts.qbox)
    {
      C3qq_delta_NFV = K3_NFV/2.;
      C4qq_delta_NFV = K4_NFV/2.;
    }
  else
    {
      C3qq_delta_NFV = 0.;
      C4qq_delta_NFV = 0.;
    }
    
  icoeff::delta();
}

void ccoeff::calc()
{
  if (opts.mellin1d)
    calc1d();
  else
    calc2d();
}

void ccoeff::calc1d()
{
  //analytical C1 and C2 from h2calc
  if (false)
    for (int m = 0; m < mellinint::mdim; m++)
      {
	fcomplex fxn; //input of ancalc
	fxn.real = real(mellinint::Np[m]);
	fxn.imag = imag(mellinint::Np[m]);
	complex <double> cxn = mellinint::Np[m];
	fcomplex fC2qg,fC2NSqqb,fC2NSqq,fC2Sqqb;
	dyh2calc_(fC2qg,fC2NSqqb,fC2NSqq,fC2Sqqb,fxn);

	int idx = anomalous::index(m,mesq::positive);
	
	C1qg[idx] = 1./((cxn+1.)*(cxn+2.))                             /2.;
	C1qq[idx] = (2.*constants::pi2/3.-16./3.+4./3./(cxn*(cxn+1.))) /2.;
	C1qgamma[idx] = 1./((cxn+1.)*(cxn+2.))                             /2.*6.; // QED C1qgamma = C1qg*Nc/TR
	C1qqqed[idx] = (2.*constants::pi2/3.-16./3.+4./3./(cxn*(cxn+1.))) /2.*(3./4.);// QED C1qqqed = C1qq/CF

	C2qg[idx]   =  cx(fC2qg)                 /4.;
	C2qq[idx]   = (cx(fC2NSqq)+cx(fC2Sqqb))  /4.;
	C2qqb[idx]  = (cx(fC2NSqqb)+cx(fC2Sqqb)) /4.;
	C2qqp[idx]  =  cx(fC2Sqqb)               /4.;
	C2qqbp[idx] =  cx(fC2Sqqb)               /4.;
      }

  //Evaluate numerically C1, C2 and C3 with I functions
  if (false)
    {
      icoeff::calc1d();
      for (int m = 0; m < mellinint::mdim; m++)
	{
	  int idx = anomalous::index(m,mesq::positive);
	  C1qg[idx]   = icoeff::c1qg[idx];
	  C1qq[idx]   = icoeff::c1qq[idx];
	  C1qgamma[idx]   = icoeff::c1qg[idx]*6.; // QED
	  C1qqqed[idx]    = icoeff::c1qq[idx]*(3./4.); // QED
	  
	  C2qg[idx]   = icoeff::c2qg[idx];
	  C2qq[idx]   = icoeff::c2qq[idx];
	  C2qqb[idx]  = icoeff::c2qqb[idx];
	  C2qqp[idx]  = icoeff::c2qqp[idx];
	  C2qqbp[idx] = icoeff::c2qqbp[idx];
	  
	  C3qg[idx]   = icoeff::c3qg[idx];
	  C3qq[idx]   = icoeff::c3qq[idx];
	  C3qqb[idx]  = icoeff::c3qqb[idx];
	  C3qqp[idx]  = icoeff::c3qqp[idx];
	  C3qqbp[idx] = icoeff::c3qqbp[idx];
	}
    }

  //Analytic C1, C2 and C3 from systematic expansion of harmonic sums
  for (int m = 0; m < mellinint::mdim; m++)
    {
      int idx = anomalous::index(m,mesq::positive);
      complex <double> N = mellinint::Np[m]-1.;
      C1qg[idx]   = C1qgN(N);
      C1qq[idx]   = C1qqN(N);
      C1qgamma[idx]   = C1qgN(N)*6.;      // QED
      C1qqqed[idx]    = C1qqN(N)*(3./4.); // QED
      C1gq[idx]   = C1gqN(N);

      C2qg[idx]   = C2qgN(N);
      C2qq[idx]   = C2qqN(N);
      C2qqb[idx]  = C2qqbN(N);
      C2qqp[idx]  = C2qqpN(N);
      C2qqbp[idx] = C2qqp[idx];
      C2gq[idx]   = C2gqN(N);

      C3qg[idx]   = C3qgN(N);
      C3qq[idx]   = C3qqN(N);
      C3qqb[idx]  = C3qqbN(N);
      C3qqp[idx]  = C3qqpN(N);
      C3qqbp[idx] = C3qqbpN(N);

      //accelerate I4
      complex <double> I1qq   = C1qq[idx] - C1qq_delta;
      complex <double> I1qg   = C1qg[idx];
      complex <double> I2qq   = C2qq[idx] - C1qq_delta*I1qq - C2qq_delta;
      complex <double> I2qg   = C2qg[idx] - C1qq_delta*I1qg;
      complex <double> I2qqb  = C2qqb[idx];
      complex <double> I2qqp  = C2qqbp[idx];
      complex <double> I2qqbp = C2qqbp[idx];
      complex <double> I3qq   = C3qq[idx]   - C1qq_delta*I2qq - C2qq_delta*I1qq - C3qq_delta;
      complex <double> I3qg   = C3qg[idx]   - C1qq_delta*I2qg - C2qq_delta*I1qg;
      complex <double> I3qqb  = C3qqb[idx]  - C1qq_delta*I2qqb;
      complex <double> I3qqbp = C3qqbp[idx] - C1qq_delta*I2qqbp;
      complex <double> I3qqp  = C3qqp[idx]  - C1qq_delta*I2qqp;
      complex <double> I4qq   = (opts.C4qq  +1.)*accel::levin(I1qq,I2qq,I3qq);
      complex <double> I4qg   = (opts.C4qg  +1.)*accel::levin(I1qg,I2qg,I3qg);
      complex <double> I4qqb  = (opts.C4qqb +1.)*accel::levin(I2qqb,I3qqb);
      complex <double> I4qqp  = (opts.C4qqp +1.)*accel::levin(I2qqp,I3qqp);
      complex <double> I4qqbp = (opts.C4qqbp+1.)*accel::levin(I2qqbp,I3qqbp);
      C4qq[idx]   = I4qq  + C1qq_delta*I3qq  + C2qq_delta*I2qq + C3qq_delta*I1qq + C4qq_delta;
      C4qg[idx]   = I4qg  + C1qq_delta*I3qg  + C2qq_delta*I2qg + C3qq_delta*I1qg;
      C4qqb[idx]  = I4qqb + C1qq_delta*I3qqb + C2qq_delta*I2qqb;
      C4qqp[idx]  = I4qqp + C1qq_delta*I3qqp + C2qq_delta*I2qqp;
      C4qqbp[idx] = I4qqbp+ C1qq_delta*I3qqbp+ C2qq_delta*I2qqbp;

      //accelerate C4
      if (false)
	{
	  C4qq  [idx] = (opts.C4qq  +1.)*accel::levin(C1qq[idx],C2qq[idx],C3qq[idx]);
	  C4qg  [idx] = (opts.C4qg  +1.)*accel::levin(C1qg[idx],C2qg[idx],C3qg[idx]);
	  C4qqb [idx] = (opts.C4qqb +1.)*accel::levin(C2qqb[idx],C3qqb[idx]);
	  C4qqp [idx] = (opts.C4qqp +1.)*accel::levin(C2qqp[idx],C3qqp[idx]);
	  C4qqbp[idx] = (opts.C4qqbp+1.)*accel::levin(C2qqbp[idx],C3qqbp[idx]);
	}

      //cout << idx << endl;
      //cout << "C4qq "<< C4qq[idx] << endl;
      //cout << "C4qg "<< C4qg[idx] << endl;
      //cout << "C4qqb "<< C4qqb[idx] << endl;
      //cout << "C4qqbp "<< C4qqbp[idx] << endl;
      //cout << "C4qqp "<< C4qqp[idx] << endl;
      //cout << endl;
      
      C1qg_reg[idx] = c1qgreg(N); //Mellin transform of -CF/48*log(1-z)^3; --> This is the most divergent part of C2 at z -> 1
      //cout << idx << "  " << C1qg_reg[idx] << endl;
      //C1qg_reg[idx] = pow(log(N+1.),3)/(N+1.);
      
      //double regfac = 1.; //Regularisation factor
      //C1qg[idx] -= resint::aass*regfac*C1qg_reg[idx];
      //C2qg[idx] += regfac*C1qg_reg[idx];

      //small-x resummation Eq.(43,44) of https://arxiv.org/pdf/1511.06039.pdf
      bool smallx = false;
      if (smallx)
	{
	  C1qg[idx] += 1./12;

	  C2qq[idx] += CF/(N-1.)*(43./108.-pi2/36.);
	  C2qg[idx] += CA/(N-1.)*(43./108.-pi2/36.);
	}
    }

  if (opts.cmw)
    for (int m = 0; m < mellinint::mdim; m++)
      {
	int idx = anomalous::index(m,mesq::positive);
	complex <double> C2qgCMW = C2qg[idx] - K1*C1qg[idx];
	complex <double> C2qqCMW = C2qq[idx] - K1*C1qq[idx];
	complex <double> C3qgCMW = C3qg[idx] - 2.*K1*C2qg[idx] + 2.*pow(K1,2)*C1qg[idx] - K2*C1qg[idx];
	complex <double> C3qqCMW = C3qq[idx] - 2.*K1*C2qq[idx] + 2.*pow(K1,2)*C1qq[idx] - K2*C1qq[idx];
	
	C2qg[idx]  = C2qgCMW;
	C2qq[idx]  = C2qqCMW;
	C3qg[idx]  = C3qgCMW;
	C3qq[idx]  = C3qqCMW;
    }

  if (opts.resscheme > 0)
    for (int m = 0; m < mellinint::mdim; m++)
      {
	//see 3.73 of https://inspirehep.net/literature/1839718
	//These formulas are recursive, if the orders is changed they do not work correctly
	int idx = anomalous::index(m,mesq::positive);
	//NLL
	C1qq[idx]   += - H1q/2.;
	C1qqqed[idx]   += - H1qp/2.; // QED
	//NNLL
	C2qq[idx]   += - H1q/2.*C1qq[idx] - (H2q/2.-pow(H1q,2)/8.);
	C2qg[idx]   += - H1q/2.*C1qg[idx];
	C3qqb[idx]  += - H1q/2.*C2qqb[idx];
	C3qqbp[idx] += - H1q/2.*C2qqbp[idx];
	C3qqp[idx]  += - H1q/2.*C2qqp[idx];
	//N3LL
	C3qq[idx]   += - H1q/2.*C2qq[idx] - (H2q/2.-pow(H1q,2)/8.)*C1qq[idx] - (H3q/2.+pow(H1q,3)/16.-1./4.*H1q*H2q);
	C3qg[idx]   += - H1q/2.*C2qg[idx] - (H2q/2.-pow(H1q,2)/8.)*C1qg[idx];
	C4qqb[idx]  += - H1q/2.*C3qqb[idx]- (H2q/2.-pow(H1q,2)/8.)*C2qqb[idx];
	C4qqp[idx]  += - H1q/2.*C3qqp[idx]- (H2q/2.-pow(H1q,2)/8.)*C2qqp[idx];
	C4qqbp[idx] += - H1q/2.*C3qqbp[idx]- (H2q/2.-pow(H1q,2)/8.)*C2qqbp[idx];
	//N4LL
	C4qq[idx]   += - H1q/2.*C3qq[idx] - (H2q/2.-pow(H1q,2)/8.)*C2qq[idx] - (H3q/2.+pow(H1q,3)/16.-1./4.*H1q*H2q)*C1qq[idx] - (H4q/2. - 5./128.*pow(H1q,4)-H1q*H3q/4.-pow(H2q,2)/8. + H2q*pow(H1q,2)*3./16.);
	C4qg[idx]   += - H1q/2.*C3qg[idx] - (H2q/2.-pow(H1q,2)/8.)*C2qg[idx] - (H3q/2.+pow(H1q,3)/16.-1./4.*H1q*H2q)*C1qg[idx];
	//C5qqb  += - H1q/2.*C4qqb - (H2q/2.-pow(H1q,2)/8.)*C3qqb - (H3q/2.+pow(H1q,3)/16.-1./4.*H1q*H2q)*C2qqb;
	//C5qqp  += - H1q/2.*C4qqp - (H2q/2.-pow(H1q,2)/8.)*C3qqp - (H3q/2.+pow(H1q,3)/16.-1./4.*H1q*H2q)*C2qqp;
	//C5qqbp  += - H1q/2.*C4qqbp - (H2q/2.-pow(H1q,2)/8.)*C3qqbp - (H3q/2.+pow(H1q,3)/16.-1./4.*H1q*H2q)*C2qqbp;
      }
  
  //Compute negative branch
  for (int m = 0; m < mellinint::mdim; m++)
   {
     int idxp = anomalous::index(m,mesq::positive);
     int idxm = anomalous::index(m,mesq::negative);

     C1qg[idxm]   = conj(C1qg[idxp]);
     C1qq[idxm]   = conj(C1qq[idxp]);
     C1qgamma[idxm]   = conj(C1qgamma[idxp]);
     C1qqqed[idxm]    = conj(C1qqqed[idxp]);
     C1gq[idxm]   = conj(C1gq[idxp]);

     C2qg[idxm]   = conj(C2qg[idxp]);
     C2qq[idxm]   = conj(C2qq[idxp]);
     C2qqb[idxm]  = conj(C2qqb[idxp]);
     C2qqp[idxm]  = conj(C2qqp[idxp]);
     C2qqbp[idxm] = conj(C2qqbp[idxp]);
     C2gq[idxm]   = conj(C2gq[idxp]);

     C3qg[idxm]   = conj(C3qg[idxp]);
     C3qq[idxm]   = conj(C3qq[idxp]);
     C3qqb[idxm]  = conj(C3qqb[idxp]);
     C3qqp[idxm]  = conj(C3qqp[idxp]);
     C3qqbp[idxm] = conj(C3qqbp[idxp]);

     C4qg[idxm]   = conj(C4qg[idxp]);
     C4qq[idxm]   = conj(C4qq[idxp]);
     C4qqb[idxm]  = conj(C4qqb[idxp]);
     C4qqp[idxm]  = conj(C4qqp[idxp]);
     C4qqbp[idxm] = conj(C4qqbp[idxp]);
     
     C1qg_reg[idxm] = conj(C1qg_reg[idxp]);
   }

  //Checks
  /*
  for (int m = 0; m < mellinint::mdim; m++)
   {
     complex <double> cxn = mellinint::Np[m];// - 1.;
     fcomplex fxn; //input of ancalc
     fxn.real = real(cxn);
     fxn.imag = imag(cxn);
     fcomplex fC2qg,fC2NSqqb,fC2NSqq,fC2Sqqb;
     h2calc_(fC2qg,fC2NSqqb,fC2NSqq,fC2Sqqb,fxn);
     fcomplex dyfC2qg,dyfC2NSqqb,dyfC2NSqq,dyfC2Sqqb;
     dyh2calc_(dyfC2qg,dyfC2NSqqb,dyfC2NSqq,dyfC2Sqqb,fxn);
     
      int idx = anomalous::index(m,mesq::positive);
      complex <double> N = cxn-1.;
      cout << endl;
      cout << " anomalous    expanded " << endl;
      cout << cxn << " C1qg   " << 1./((cxn+1.)*(cxn+2.))                             /2.  - C1qgN(N)  << endl;
      cout << cxn << " C1qq   " << (2.*constants::pi2/3.-16./3.+4./3./(cxn*(cxn+1.))) /2.  - C1qqN(N)  << endl;


      cout << cxn << " C2qg   " << cx(fC2qg)/4.                 - C2qgN(N)   << "  " << cx(dyfC2qg)/4.                    - C2qgN(N)  << endl;
      cout << cxn << " C2qq   " << (cx(fC2NSqq)+cx(fC2Sqqb))/4. - C2qqN(N)   << "  " << (cx(dyfC2NSqq)+cx(dyfC2Sqqb))/4.  - C2qqN(N)  << endl;
      cout << cxn << " C2qqp  " << cx(fC2Sqqb)/4.               - C2qqpN(N)  << "  " << cx(dyfC2Sqqb)/4.                  - C2qqpN(N) << endl;
      cout << cxn << " C2qqb  " << (cx(fC2NSqqb)+cx(fC2Sqqb))/4.- C2qqbN(N)  << "  " << (cx(dyfC2NSqqb)+cx(dyfC2Sqqb))/4. - C2qqbN(N) << endl;
    }
  */
  
  /*
  for (int m = 0; m < mellinint::mdim; m++)
   {
      int idx = anomalous::index(m,mesq::positive);
      complex <double> N = mellinint::Np[m]-1.;
      cout << endl;
      cout << " anomalous    expanded    numerical" << endl;
      cout << mellinint::Np[m] << " C2qg   " << anomalous::C2qgM[idx]/4.                              << "  " << C2qgN(N)  << "  " << icoeff::c2qg[idx]  << endl;
      cout << mellinint::Np[m] << " C2qq   " << (anomalous::C2NSqqM[idx]+anomalous::C2SqqbM[idx])/4.  << "  " << C2qqN(N)  << "  " << icoeff::c2qq[idx]  << endl;
      cout << mellinint::Np[m] << " C2qqp  " << anomalous::C2SqqbM[idx]/4.                            << "  " << C2qqpN(N) << "  " << icoeff::c2qqp[idx] << endl;
      cout << mellinint::Np[m] << " C2qqb  " << (anomalous::C2NSqqbM[idx]+anomalous::C2SqqbM[idx])/4. << "  " << C2qqbN(N) << "  " << icoeff::c2qqb[idx] << endl;

      cout << endl;      
      cout << mellinint::Np[m] << " C3qg   " << icoeff::c3qg[idx] - C3qgN(N)	   << endl;
      cout << mellinint::Np[m] << " C3qq   " << icoeff::c3qq[idx] - C3qqN(N)	   << endl;
      cout << mellinint::Np[m] << " C3qqp  " << icoeff::c3qqp[idx] - C3qqpN(N)   << endl;
      cout << mellinint::Np[m] << " C3qqb  " << icoeff::c3qqb[idx] - C3qqbN(N)   << endl;
      cout << mellinint::Np[m] << " C3qqbp " << icoeff::c3qqbp[idx] - C3qqbpN(N) << endl;
    }
  */
  
/*
  int prec = 8;
  int dist = 2*prec+20;
  cout << setprecision(8) << endl;
  cout << setw(10) << "C2qqn = " << resconst::C2qqn  << endl;
  cout << setw(10) << "H2qqD0 "    << resconst::H2qqD0         << endl;
  cout << setw(10) << "H2qqdelta " << resconst::H2qqdelta      << endl;
  cout << setw(dist) << " N "            << setw(10) << "coeff "     << setw(dist) << " analytical " << " function " << setw(dist) << " numerical "       << setw(dist/2) <<" Rel Error " << endl;
  for (int m = 0; m < mellinint::mdim; m++)
   {
     fcomplex xn = fcx(mellinint::Np[m]);
     fcomplex ps0n;
     psi0_(xn,ps0n);
     complex <double> s1nm1 = cx(ps0n)+resconst::Euler;

     int idx = anomalous::index(m,mesq::positive);
     cout << setw(dist) << mellinint::Np[m] << setw(10) << "C2QG "      << setw(dist) << anomalous::C2qgM[idx]    << " 4*c2qg          " << setw(dist) << 4.*c2qg[idx]        << setw(dist/2) << fabs(anomalous::C2qgM[idx]/4./c2qg[idx]-1.) << endl;
     cout << setw(dist) << mellinint::Np[m] << setw(10) << "C2NSQQ "    << setw(dist) << anomalous::C2NSqqM[idx]  << " 4*(c2qqreg+C2qqn+H2qqD0/2-c2qqp)       " << setw(dist) << 4.*(c2NSqq[idx]+resconst::C2qqn+(-s1nm1*resconst::H2qqD0/2.))      << setw(dist/2) << fabs(anomalous::C2NSqqM[idx]/(4.*(c2NSqq[idx]+resconst::C2qqn+(-s1nm1*resconst::H2qqD0/2.)))-1.) << endl;
     cout << setw(dist) << mellinint::Np[m] << setw(10) << "C2SQQB "    << setw(dist) << anomalous::C2SqqbM[idx]  << " 4*c2qqp         " << setw(dist) << 4.*c2Sqqb[idx]      << setw(dist/2) << fabs(anomalous::C2SqqbM[idx]/4./c2Sqqb[idx]-1.) << endl;
     cout << setw(dist) << mellinint::Np[m] << setw(10) << "C2NSQQB "   << setw(dist) << anomalous::C2NSqqbM[idx] << " 4*(c2qqb-c2qqp)         " << setw(dist) << 4.*c2NSqqb[idx]     << setw(dist/2) << fabs(anomalous::C2NSqqbM[idx]/4./c2NSqqb[idx]-1.) << endl;
     cout << setw(dist) << mellinint::Np[m] << setw(10) << "ps0n "      << setw(dist) << cx(ps0n) << endl;
     cout << setw(dist) << mellinint::Np[m] << setw(10) << "s1nm1 "     << setw(dist) << s1nm1 << endl;
     cout << endl;
   }
*/

  /*
  int idx;
  for (int m = 0; m < mellinint::mdim; m++)
   {
     fcomplex xn = fcx(mellinint::Np[m]);
     fcomplex ps0n;
     psi0_(xn,ps0n);
     complex <double> s1nm1 = cx(ps0n)+resconst::Euler;

     int idxp = anomalous::index(m,mesq::positive);
     anomalous::C2qgM[idxp]    = 4.*c2qg[idxp];
     anomalous::C2NSqqM[idxp]  = 4.*(c2NSqq[idxp]+resconst::C2qqn+(-s1nm1*resconst::H2qqD0/2.));
     anomalous::C2SqqbM[idxp]  = 4.*c2Sqqb[idxp];
     anomalous::C2NSqqbM[idxp] = 4.*c2NSqqb[idxp];                                                
     int idxm = anomalous::index(m,mesq::negative);
     anomalous::C2qgM[idxm]    = conj(anomalous::C2qgM[idxp]);
     anomalous::C2NSqqM[idxm]  = conj(anomalous::C2NSqqM[idxp]);
     anomalous::C2SqqbM[idxm]  = conj(anomalous::C2SqqbM[idxp]);
     anomalous::C2NSqqbM[idxm] = conj(anomalous::C2NSqqbM[idxp]);

     //cout << conj(anomalous::C2qgM[idxp])    <<	 "  " << 4.*conj(c2qg[idxp]) << endl;                                                         
     //cout << conj(anomalous::C2NSqqM[idxp])  <<	 "  " << 4.*(conj(c2NSqq[idxp])+resconst::C2qqn+(-conj(s1nm1)*resconst::H2qqD0/2.)) << endl;  
     //cout << conj(anomalous::C2SqqbM[idxp])  <<	 "  " << 4.*conj(c2Sqqb[idxp]) << endl;						       
     //cout << conj(anomalous::C2NSqqbM[idxp]) <<	 "  " << 4.*conj(c2NSqqb[idxp]) << endl;						       
     //cout << endl;
   }
  */
}

void ccoeff::calc2d()
{
  //analytical C1 and C2
  if (false)
    {
      for (int m = 0; m < mellinint::mdim; m++)
	{
	  int idx = anomalous::index(m,mesq::positive);
      
	  C1qg_1[idx] = anomalous::C1QG_1[idx]/2.;
	  C1qq_1[idx] = anomalous::C1QQ_1[idx]/2.;
      C1gq_1[idx] = anomalous::C1GQ_1[idx]/2.;
	  C1qgamma_1[idx] = anomalous::C1QG_1[idx]/2.*6.; // QED C1qgamma = C1qg Nc/TR
	  C1qqqed_1[idx]  = anomalous::C1QQ_1[idx]/2.*(3./4.); // QED C1qqqed  = C1qq /CF
	  C2qg_1[idx]   = anomalous::C2qgM_1[idx]/4.;
	  C2qq_1[idx]   = (anomalous::C2NSqqM_1[idx]+anomalous::C2SqqbM_1[idx])/4.;
	  C2qqb_1[idx]  = (anomalous::C2NSqqbM_1[idx]+anomalous::C2SqqbM_1[idx])/4.;
	  C2qqp_1[idx]  = anomalous::C2SqqbM_1[idx]/4.;
	  C2qqbp_1[idx] = anomalous::C2SqqbM_1[idx]/4.;

	  C2qg_2[idx]   = anomalous::C2qgM_2[idx]/4.;
	  C2qq_2[idx]   = (anomalous::C2NSqqM_2[idx]+anomalous::C2SqqbM_2[idx])/4.;
      C1gq_2[idx] = anomalous::C1GQ_2[idx]/2.;
	  C2qqb_2[idx]  = (anomalous::C2NSqqbM_2[idx]+anomalous::C2SqqbM_2[idx])/4.;
	  C2qqp_2[idx]  = anomalous::C2SqqbM_2[idx]/4.;
	  C2qqbp_2[idx] = anomalous::C2SqqbM_2[idx]/4.;
	  C1qg_2[idx] = anomalous::C1QG_2[idx]/2.;
	  C1qq_2[idx] = anomalous::C1QQ_2[idx]/2.;
	  C1qgamma_2[idx] = anomalous::C1QG_2[idx]/2.*6.; // QED C1qgamma = C1qg Nc/TR
	  C1qqqed_2[idx] = anomalous::C1QQ_2[idx]/2.*(3./4.); // QED C1qqqed  = C1qq /CF 
	}

      //Evaluate C3 coefficients with I functions
      icoeff::calc2d();

      for (int m = 0; m < mellinint::mdim; m++)
	{
	  int idx = anomalous::index(m,mesq::positive);

	  C3qg_1[idx]   = icoeff::c3qg_1[idx];
	  C3qq_1[idx]   = icoeff::c3qq_1[idx];
	  C3qqb_1[idx]  = icoeff::c3qqb_1[idx];
	  C3qqp_1[idx]  = icoeff::c3qqp_1[idx];
	  C3qqbp_1[idx] = icoeff::c3qqbp_1[idx];

	  C3qg_2[idx]   = icoeff::c3qg_2[idx];
	  C3qq_2[idx]   = icoeff::c3qq_2[idx];
	  C3qqb_2[idx]  = icoeff::c3qqb_2[idx];
	  C3qqp_2[idx]  = icoeff::c3qqp_2[idx];
	  C3qqbp_2[idx] = icoeff::c3qqbp_2[idx];
	}
    }

  //Analytic C1, C2 and C3 from systematic expansion of harmonic sums
  //if (false)
    for (int m = 0; m < mellinint::mdim; m++)
      {
	int idx = anomalous::index(m,mesq::positive);

    fcomplex fxn_1; // Input of mellinh2
    complex <double> cxn_1; // Input of G1N
    fxn_1.real = real(mellinint::Np_1[m]);
    fxn_1.imag = imag(mellinint::Np_1[m]);
    cxn_1 = mellinint::Np_1[m];

    // NLL Higgs
    C1gg = C1ggN;

	complex <double> N_1 = mellinint::Np_1[m]-1.;
	C1qg_1[idx]   = C1qgN(N_1);
    C1gq_1[idx]   = C1gqN(N_1);
	C1qq_1[idx]   = C1qqN(N_1);
	C1qgamma_1[idx]   = C1qgN(N_1)*6.; // QED C1qgamma = C1qg Nc/TR  
	C1qqqed_1[idx]    = C1qqN(N_1)*(3./4.); // QED C1qqqed  = C1qq /CF
    H2ggM_1[idx]  = cx(mellinh2gg_(fxn_1));
    H2gqM_1[idx]  = cx(mellinh2gq_(fxn_1));
    H2qqM_1[idx]  = cx(mellinh2qq_(fxn_1));
    G1N_1[idx]    = G1N(N_1);
	C2qg_1[idx]   = C2qgN(N_1);
	C2qq_1[idx]   = C2qqN(N_1);
    C2gq_1[idx]   = H2gqM_1[idx]-C1gg*C1gq_1[idx];
	C2qqb_1[idx]  = C2qqbN(N_1);
	C2qqp_1[idx]  = C2qqpN(N_1);
	C2qqbp_1[idx] = C2qqp_1[idx];
	C3qg_1[idx]   = C3qgN(N_1);
	C3qq_1[idx]   = C3qqN(N_1);
	C3qqb_1[idx]  = C3qqbN(N_1);
	C3qqp_1[idx]  = C3qqpN(N_1);
	C3qqbp_1[idx] = C3qqbpN(N_1);
	complex <double> I1qq_1   = C1qq_1[idx] - C1qq_delta;
	complex <double> I1qg_1   = C1qg_1[idx];
	complex <double> I2qq_1   = C2qq_1[idx] - C1qq_delta*I1qq_1 - C2qq_delta;
	complex <double> I2qg_1   = C2qg_1[idx] - C1qq_delta*I1qg_1;
	complex <double> I2qqb_1  = C2qqb_1[idx];
	complex <double> I2qqp_1  = C2qqbp_1[idx];
	complex <double> I2qqbp_1 = C2qqbp_1[idx];
	complex <double> I3qq_1   = C3qq_1[idx]   - C1qq_delta*I2qq_1 - C2qq_delta*I1qq_1 - C3qq_delta;
	complex <double> I3qg_1   = C3qg_1[idx]   - C1qq_delta*I2qg_1 - C2qq_delta*I1qg_1;
	complex <double> I3qqb_1  = C3qqb_1[idx]  - C1qq_delta*I2qqb_1;
	complex <double> I3qqbp_1 = C3qqbp_1[idx] - C1qq_delta*I2qqbp_1;
	complex <double> I3qqp_1  = C3qqp_1[idx]  - C1qq_delta*I2qqp_1;
	complex <double> I4qq_1   = (opts.C4qq  +1.)*accel::levin(I1qq_1,I2qq_1,I3qq_1);
	complex <double> I4qg_1   = (opts.C4qg  +1.)*accel::levin(I1qg_1,I2qg_1,I3qg_1);
	complex <double> I4qqb_1  = (opts.C4qqb +1.)*accel::levin(I2qqb_1,I3qqb_1);
	complex <double> I4qqp_1  = (opts.C4qqp +1.)*accel::levin(I2qqp_1,I3qqp_1);
	complex <double> I4qqbp_1 = (opts.C4qqbp+1.)*accel::levin(I2qqbp_1,I3qqbp_1);
	C4qq_1[idx]   = I4qq_1  + C1qq_delta*I3qq_1  + C2qq_delta*I2qq_1 + C3qq_delta*I1qq_1 + C4qq_delta;
	C4qg_1[idx]   = I4qg_1  + C1qq_delta*I3qg_1  + C2qq_delta*I2qg_1 + C3qq_delta*I1qg_1;
	C4qqb_1[idx]  = I4qqb_1 + C1qq_delta*I3qqb_1 + C2qq_delta*I2qqb_1;
	C4qqp_1[idx]  = I4qqp_1 + C1qq_delta*I3qqp_1 + C2qq_delta*I2qqp_1;
	C4qqbp_1[idx] = I4qqbp_1+ C1qq_delta*I3qqbp_1+ C2qq_delta*I2qqbp_1;

	C1qg_reg_1[idx] = c1qgreg(N_1);

    fcomplex fxn_2; // Input of mellinh2
    complex <double> cxn_2; // Input of G1N
    fxn_2.real = real(mellinint::Np_2[m]);
    fxn_2.imag = imag(mellinint::Np_2[m]);
    cxn_2 = mellinint::Np_2[m];

	complex <double> N_2 = mellinint::Np_2[m]-1.;
	C1qg_2[idx]   = C1qgN(N_2);
	C1qq_2[idx]   = C1qqN(N_2);
    C1gq_2[idx]   = C1gqN(N_2);
	C1qgamma_2[idx]   = C1qgN(N_2)*6.; // QED
	C1qqqed_2[idx]    = C1qqN(N_2)*(3./4.); //QED
    H2ggM_2[idx]  = cx(mellinh2gg_(fxn_2));
    H2gqM_2[idx]  = cx(mellinh2gq_(fxn_2));
    H2qqM_2[idx]  = cx(mellinh2qq_(fxn_2));
    G1N_2[idx]    = G1N(N_2);
	C2qg_2[idx]   = C2qgN(N_2);
	C2qq_2[idx]   = C2qqN(N_2);
    C2gq_2[idx]   = H2gqM_2[idx]-C1gg*C1gq_2[idx];
	C2qqb_2[idx]  = C2qqbN(N_2);
	C2qqp_2[idx]  = C2qqpN(N_2);
	C2qqbp_2[idx] = C2qqp_2[idx];
	C3qg_2[idx]   = C3qgN(N_2);
	C3qq_2[idx]   = C3qqN(N_2);
	C3qqb_2[idx]  = C3qqbN(N_2);
	C3qqp_2[idx]  = C3qqpN(N_2);
	C3qqbp_2[idx] = C3qqbpN(N_2);
	complex <double> I1qq_2   = C1qq_2[idx] - C1qq_delta;
	complex <double> I1qg_2   = C1qg_2[idx];
	complex <double> I2qq_2   = C2qq_2[idx] - C1qq_delta*I1qq_2 - C2qq_delta;
	complex <double> I2qg_2   = C2qg_2[idx] - C1qq_delta*I1qg_2;
	complex <double> I2qqb_2  = C2qqb_2[idx];
	complex <double> I2qqp_2  = C2qqbp_2[idx];
	complex <double> I2qqbp_2 = C2qqbp_2[idx];
	complex <double> I3qq_2   = C3qq_2[idx]   - C1qq_delta*I2qq_2 - C2qq_delta*I1qq_2 - C3qq_delta;
	complex <double> I3qg_2   = C3qg_2[idx]   - C1qq_delta*I2qg_2 - C2qq_delta*I1qg_2;
	complex <double> I3qqb_2  = C3qqb_2[idx]  - C1qq_delta*I2qqb_2;
	complex <double> I3qqbp_2 = C3qqbp_2[idx] - C1qq_delta*I2qqbp_2;
	complex <double> I3qqp_2  = C3qqp_2[idx]  - C1qq_delta*I2qqp_2;
	complex <double> I4qq_2   = (opts.C4qq  +1.)*accel::levin(I1qq_2,I2qq_2,I3qq_2);
	complex <double> I4qg_2   = (opts.C4qg  +1.)*accel::levin(I1qg_2,I2qg_2,I3qg_2);
	complex <double> I4qqb_2  = (opts.C4qqb +1.)*accel::levin(I2qqb_2,I3qqb_2);
	complex <double> I4qqp_2  = (opts.C4qqp +1.)*accel::levin(I2qqp_2,I3qqp_2);
	complex <double> I4qqbp_2 = (opts.C4qqbp+1.)*accel::levin(I2qqbp_2,I3qqbp_2);
	C4qq_2[idx]   = I4qq_2  + C1qq_delta*I3qq_2  + C2qq_delta*I2qq_2 + C3qq_delta*I1qq_2 + C4qq_delta;
	C4qg_2[idx]   = I4qg_2  + C1qq_delta*I3qg_2  + C2qq_delta*I2qg_2 + C3qq_delta*I1qg_2;
	C4qqb_2[idx]  = I4qqb_2 + C1qq_delta*I3qqb_2 + C2qq_delta*I2qqb_2;
	C4qqp_2[idx]  = I4qqp_2 + C1qq_delta*I3qqp_2 + C2qq_delta*I2qqp_2;
	C4qqbp_2[idx] = I4qqbp_2+ C1qq_delta*I3qqbp_2+ C2qq_delta*I2qqbp_2;
	C1qg_reg_2[idx] = c1qgreg(N_2);
      }
  
  if (opts.cmw)
    for (int m = 0; m < mellinint::mdim; m++)
      {
	int idx = anomalous::index(m,mesq::positive);
	complex <double> C2qgCMW,C2qqCMW,C3qgCMW,C3qqCMW;

	C2qgCMW = C2qg_1[idx] - K1*C1qg_1[idx];
	C2qqCMW = C2qq_1[idx] - K1*C1qq_1[idx];
	C3qgCMW = C3qg_1[idx] - 2.*K1*C2qg_1[idx] + 2.*pow(K1,2)*C1qg_1[idx] - K2*C1qg_1[idx];
	C3qqCMW = C3qq_1[idx] - 2.*K1*C2qq_1[idx] + 2.*pow(K1,2)*C1qq_1[idx] - K2*C1qq_1[idx];
	C2qg_1[idx]  = C2qgCMW;
	C2qq_1[idx]  = C2qqCMW;
	C3qg_1[idx]  = C3qgCMW;
	C3qq_1[idx]  = C3qqCMW;

	C2qgCMW = C2qg_2[idx] - K1*C1qg_2[idx];
	C2qqCMW = C2qq_2[idx] - K1*C1qq_2[idx];
	C3qgCMW = C3qg_2[idx] - 2.*K1*C2qg_2[idx] + 2.*pow(K1,2)*C1qg_2[idx] - K2*C1qg_2[idx];
	C3qqCMW = C3qq_2[idx] - 2.*K1*C2qq_2[idx] + 2.*pow(K1,2)*C1qq_2[idx] - K2*C1qq_2[idx];
	C2qg_2[idx]  = C2qgCMW;
	C2qq_2[idx]  = C2qqCMW;
	C3qg_2[idx]  = C3qgCMW;
	C3qq_2[idx]  = C3qqCMW;
      }

  if (opts.resscheme == 5)
    for (int m = 0; m < mellinint::mdim; m++)
      {
    	int idx = anomalous::index(m,mesq::positive);
	//convert from C to I
	C2qq_1[idx]   += - C1qq_delta*(C1qq_1[idx]-C1qq_delta);
	C2qg_1[idx]   += - C1qq_delta*C1qg_1[idx];
	C2qq_2[idx]   += - C1qq_delta*(C1qq_2[idx]-C1qq_delta);
	C2qg_2[idx]   += - C1qq_delta*C1qg_2[idx];

	//delta part as in https://arxiv.org/pdf/1604.07869.pdf 7.2, 7.3
	double c1qd = (-CF*pi2/6.)/4.; // From Eq. (7.2) of https://arxiv.org/pdf/1604.07869.pdf
	double c2qd = (CF*CF*pi4/72. + CF*CA*(1214./81.-67.*pi2/36.-77./9.*zeta3+pi4/18.)+CF*TF*NF*(-328/81.+pi2*5./9.+28./9.*zeta3))/16.;
    	C1qq_1[idx]   += - H1q/2.;
    	C1qq_2[idx]   += - H1q/2.;
	C1qqqed_1[idx]   += - H1qp/2.; // QED
	C1qqqed_2[idx]   += - H1qp/2.; // QED
    	C2qq_1[idx]   += - H1q/2.*(C1qq_delta-H1q/2.) - (H2q/2.-pow(H1q,2)/8.);
    	C2qq_2[idx]   += - H1q/2.*(C1qq_delta-H1q/2.) - (H2q/2.-pow(H1q,2)/8.);
    	//C1qq_1[idx]   += c1qd - C1qq_delta;
    	//C1qq_2[idx]   += c1qd - C1qq_delta;
    	//C2qq_1[idx]   += c2qd - C2qq_delta;
    	//C2qq_2[idx]   += c2qd - C2qq_delta;
      }
  else if (opts.resscheme > 0)
    for (int m = 0; m < mellinint::mdim; m++)
      {
	//see 3.73 of https://inspirehep.net/literature/1839718
	//These formulas are recursive, if the orders is changed they do not work correctly
	int idx = anomalous::index(m,mesq::positive);
	//NLL
	C1qq_1[idx]   += - H1q/2.;
	//NNLL
	C2qq_1[idx]   += - H1q/2.*C1qq_1[idx] - (H2q/2.-pow(H1q,2)/8.);
	C2qg_1[idx]   += - H1q/2.*C1qg_1[idx];
	C3qqb_1[idx]  += - H1q/2.*C2qqb_1[idx];
	C3qqbp_1[idx] += - H1q/2.*C2qqbp_1[idx];
	C3qqp_1[idx]  += - H1q/2.*C2qqp_1[idx];
	//N3LL
	C3qq_1[idx]   += - H1q/2.*C2qq_1[idx] - (H2q/2.-pow(H1q,2)/8.)*C1qq_1[idx] - (H3q/2.+pow(H1q,3)/16.-1./4.*H1q*H2q);
	C3qg_1[idx]   += - H1q/2.*C2qg_1[idx] - (H2q/2.-pow(H1q,2)/8.)*C1qg_1[idx];
	C4qqb_1[idx]  += - H1q/2.*C3qqb_1[idx]- (H2q/2.-pow(H1q,2)/8.)*C2qqb_1[idx];
	C4qqp_1[idx]  += - H1q/2.*C3qqp_1[idx]- (H2q/2.-pow(H1q,2)/8.)*C2qqp_1[idx];
	C4qqbp_1[idx] += - H1q/2.*C3qqbp_1[idx]- (H2q/2.-pow(H1q,2)/8.)*C2qqbp_1[idx];
	//N4LL
	C4qq_1[idx]   += - H1q/2.*C3qq_1[idx] - (H2q/2.-pow(H1q,2)/8.)*C2qq_1[idx] - (H3q/2.+pow(H1q,3)/16.-1./4.*H1q*H2q)*C1qq_1[idx] - (H4q/2. - 5./128.*pow(H1q,4)-H1q*H3q/4.-pow(H2q,2)/8. + H2q*pow(H1q,2)*3./16.);
	C4qg_1[idx]   += - H1q/2.*C3qg_1[idx] - (H2q/2.-pow(H1q,2)/8.)*C2qg_1[idx] - (H3q/2.+pow(H1q,3)/16.-1./4.*H1q*H2q)*C1qg_1[idx];
	//C5qqb_1  += - H1q/2.*C4qqb_1 - (H2q/2.-pow(H1q,2)/8.)*C3qqb_1 - (H3q/2.+pow(H1q,3)/16.)*C2qqb_1;
	//NLL
	C1qq_2[idx]   += - H1q/2.;
	C1qqqed_2[idx]   += - H1qp/2.; // QED
	//NNLL
	C2qq_2[idx]   += - H1q/2.*C1qq_2[idx] - (H2q/2.-pow(H1q,2)/8.);
	C2qg_2[idx]   += - H1q/2.*C1qg_2[idx];
	C3qqb_2[idx]  += - H1q/2.*C2qqb_2[idx];
	C3qqbp_2[idx] += - H1q/2.*C2qqbp_2[idx];
	C3qqp_2[idx]  += - H1q/2.*C2qqp_2[idx];
	//N3LL
	C3qq_2[idx]   += - H1q/2.*C2qq_2[idx] - (H2q/2.-pow(H1q,2)/8.)*C1qq_2[idx] - (H3q/2.+pow(H1q,3)/16.-1./4.*H1q*H2q);
	C3qg_2[idx]   += - H1q/2.*C2qg_2[idx] - (H2q/2.-pow(H1q,2)/8.)*C1qg_2[idx];
	C4qqb_2[idx]  += - H1q/2.*C3qqb_2[idx]- (H2q/2.-pow(H1q,2)/8.)*C2qqb_2[idx];
	C4qqp_2[idx]  += - H1q/2.*C3qqp_2[idx]- (H2q/2.-pow(H1q,2)/8.)*C2qqp_2[idx];
	C4qqbp_2[idx]  += - H1q/2.*C3qqbp_2[idx]- (H2q/2.-pow(H1q,2)/8.)*C2qqbp_2[idx];
	//N4LL
	C4qq_2[idx]   += - H1q/2.*C3qq_2[idx] - (H2q/2.-pow(H1q,2)/8.)*C2qq_2[idx] - (H3q/2.+pow(H1q,3)/16.-1./4.*H1q*H2q)*C1qq_2[idx] - (H4q/2. - 5./128.*pow(H1q,4)-H1q*H3q/4.-pow(H2q,2)/8. + H2q*pow(H1q,2)*3./16.);
	C4qg_2[idx]   += - H1q/2.*C3qg_2[idx] - (H2q/2.-pow(H1q,2)/8.)*C2qg_2[idx] - (H3q/2.+pow(H1q,3)/16.-1./4.*H1q*H2q)*C1qg_2[idx];
	//C5qqb_2  += - H1q/2.*C4qqb_2 - (H2q/2.-pow(H1q,2)/8.)*C3qqb_2 - (H3q/2.+pow(H1q,3)/16.-1./4.*H1q*H2q)*C2qqb_2;
      }

  
  //Compute negative branch
  for (int m = 0; m < mellinint::mdim; m++)
   {
     int idxp = anomalous::index(m,mesq::positive);
     int idxm = anomalous::index(m,mesq::negative);

     C1qg_1[idxm]   = conj(C1qg_1[idxp]);
     C1qq_1[idxm]   = conj(C1qq_1[idxp]);
     C1gq_1[idxm]   = conj(C1gq_1[idxp]);
     C1qgamma_1[idxm]   = conj(C1qg_1[idxp])*6.;      // QED
     C1qqqed_1[idxm]    = conj(C1qq_1[idxp])*(3./4.); // QED
     H2ggM_1[idxm]  = conj(H2ggM_1[idxp]);
     H2gqM_1[idxm]  = conj(H2gqM_1[idxp]);
     H2qqM_1[idxm]  = conj(H2qqM_1[idxp]);
     G1N_1[idxm]    = conj(G1N_1[idxp]);
     C2qg_1[idxm]   = conj(C2qg_1[idxp]);
      C2gq_1[idxm]   = conj(C2gq_1[idxp]);
     C2qq_1[idxm]   = conj(C2qq_1[idxp]);
     C2qqb_1[idxm]  = conj(C2qqb_1[idxp]);
     C2qqp_1[idxm]  = conj(C2qqp_1[idxp]);
     C2qqbp_1[idxm] = conj(C2qqbp_1[idxp]);
     C3qg_1[idxm]   = conj(C3qg_1[idxp]);
     C3qq_1[idxm]   = conj(C3qq_1[idxp]);
     C3qqb_1[idxm]  = conj(C3qqb_1[idxp]);
     C3qqp_1[idxm]  = conj(C3qqp_1[idxp]);
     C3qqbp_1[idxm] = conj(C3qqbp_1[idxp]);
     C4qg_1[idxm]   = conj(C4qg_1[idxp]);
     C4qq_1[idxm]   = conj(C4qq_1[idxp]);
     C4qqb_1[idxm]  = conj(C4qqb_1[idxp]);
     C4qqp_1[idxm]  = conj(C4qqp_1[idxp]);
     C4qqbp_1[idxm] = conj(C4qqbp_1[idxp]);
     C1qg_reg_1[idxm]   = conj(C1qg_reg_1[idxp]);

     C1qg_2[idxm]   = conj(C1qg_2[idxp]);
     C1gq_2[idxm]   = conj(C1gq_2[idxp]);
     C1qq_2[idxm]   = conj(C1qq_2[idxp]);
     C1qgamma_2[idxm]   = conj(C1qg_2[idxp])*6.;      // QED
     C1qqqed_2[idxm ]   = conj(C1qq_2[idxp])*(3./4.); // QED
     H2ggM_2[idxm]  = conj(H2ggM_2[idxp]);
     H2gqM_2[idxm]  = conj(H2gqM_2[idxp]);
     H2qqM_2[idxm]  = conj(H2qqM_2[idxp]);
     G1N_2[idxm]    = conj(G1N_2[idxp]);
     C2qg_2[idxm]   = conj(C2qg_2[idxp]);
     C2gq_2[idxm]   = conj(C2gq_2[idxp]);
     C2qq_2[idxm]   = conj(C2qq_2[idxp]);
     C2qqb_2[idxm]  = conj(C2qqb_2[idxp]);
     C2qqp_2[idxm]  = conj(C2qqp_2[idxp]);
     C2qqbp_2[idxm] = conj(C2qqbp_2[idxp]);
     C3qg_2[idxm]   = conj(C3qg_2[idxp]);
     C3qq_2[idxm]   = conj(C3qq_2[idxp]);
     C3qqb_2[idxm]  = conj(C3qqb_2[idxp]);
     C3qqp_2[idxm]  = conj(C3qqp_2[idxp]);
     C3qqbp_2[idxm] = conj(C3qqbp_2[idxp]);
     C4qg_2[idxm]   = conj(C4qg_2[idxp]);
     C4qq_2[idxm]   = conj(C4qq_2[idxp]);
     C4qqb_2[idxm]  = conj(C4qqb_2[idxp]);
     C4qqp_2[idxm]  = conj(C4qqp_2[idxp]);
     C4qqbp_2[idxm] = conj(C4qqbp_2[idxp]);
     C1qg_reg_2[idxm]   = conj(C1qg_reg_2[idxp]);
   }
}

void ccoeff::num_calc()
{
  //Evaluate C1 and C2 coefficients numerically
  int rule = 500;
  double *tlog,*tlin;              //gauss nodes
  double *faclog,*faclin;            //overall factor
  complex <double> *kernlog,*kernlin; //kernel of the Mellin transform
  //memory leak here!!!
  tlog = new double [rule];
  faclog = new double [rule];
  kernlog = new complex <double> [rule*mellinint::mdim*2];
  tlin = new double [rule];
  faclin = new double [rule];
  kernlin = new complex <double> [rule*mellinint::mdim*2];

  // boundaries of integration
  double xmax = 1;
  //double xminlog = pow(bins.mbins.front()/opts.sroot,2);
  //double xminlin = pow(bins.mbins.front()/opts.sroot,2);
  double xminlog = 1e-14;
  double xminlin = 0;
  double ll = log(xmax/xminlog);
  double cc = 0.5;
  double mm = 0.5;
  for (int i = 0; i < rule; i++)
    {
      double x = cc+mm*gr::xxx[rule-1][i];
      tlog[i] = xminlog*exp(ll*x);
      tlin[i] = xminlin+(xmax-xminlin)*x;
      double jaclog = mm * tlog[i] * ll;
      double jaclin = mm*(xmax-xminlin);
      faclog[i] = jaclog * gr::www[rule-1][i];
      faclin[i] = jaclin * gr::www[rule-1][i];
      for (int m = 0; m < mellinint::mdim; m++)
	{
	  kernlog[i*mellinint::mdim+m] = pow(tlog[i], mellinint::Np[m]-1.);
	  kernlin[i*mellinint::mdim+m] = pow(tlin[i], mellinint::Np[m]-1.);
	}
    }

  fill(C1qg,   C1qg+2*mellinint::mdim, 0.);
  fill(C1qq,   C1qq+2*mellinint::mdim, 0.);
  fill(C1qgamma,   C1qgamma+2*mellinint::mdim, 0.); // QED
  fill(C1qqqed,   C1qqqed+2*mellinint::mdim, 0.);   // QED
  fill(C1qqb,  C1qqb+2*mellinint::mdim, 0.);
  fill(C1qqp,  C1qqp+2*mellinint::mdim, 0.);
  fill(C1qqbp, C1qqbp+2*mellinint::mdim, 0.);

  fill(C2qg,   C2qg+2*mellinint::mdim, 0.);
  fill(C2qq,   C2qq+2*mellinint::mdim, 0.);
  fill(C2qqb,  C2qqb+2*mellinint::mdim, 0.);
  fill(C2qqp,  C2qqp+2*mellinint::mdim, 0.);
  fill(C2qqbp, C2qqbp+2*mellinint::mdim, 0.);

  //Calculate Mellin moments as:
  //integral_0^1{ x^(N-1) fx dx}
  for (int i = 0; i < rule; i++)
    {
      double zlog = tlog[i];
      double zlin = tlin[i];

      double c1qgzlin = cqg_(zlin);
      double c1qqzlin = cqq_(zlin);
      double c1qgammazlin = cqg_(zlin)*6.;       // QED
      double c1qqqedzlin =  cqq_(zlin)*(3./4.);  // QED
      
      double c2qgzlin = c2qg_(zlin);
      double c2qqzlin = c2qqreg_(zlin);
      double c2qqbzlog = c2qqb_(zlog);
      double c2qqpzlog = c2qqp_(zlog);
      double c2qqbpzlog = c2qqpzlog; //=c2qqp_(zlog);
      for (int m = 0; m < mellinint::mdim; m++)
	{
	  int idx = anomalous::index(m,mesq::positive);
	  
	  C1qg[idx] += faclin[i]*c1qgzlin * kernlin[i*mellinint::mdim+m];
	  C1qq[idx] += faclin[i]*c1qqzlin * kernlin[i*mellinint::mdim+m];

	  C1qgamma[idx] += faclin[i]*c1qgammazlin * kernlin[i*mellinint::mdim+m];  // QED
	  C1qqqed[idx]  += faclin[i]*c1qqqedzlin * kernlin[i*mellinint::mdim+m];   // QED

	  C2qg[idx]   += faclin[i]*kernlin[i*mellinint::mdim+m] * c2qgzlin;
	  C2qq[idx]   += faclin[i]*kernlin[i*mellinint::mdim+m] * c2qqzlin;
	  C2qqb[idx]  += faclog[i]*kernlog[i*mellinint::mdim+m] * c2qqbzlog;
	  C2qqp[idx]  += faclog[i]*kernlog[i*mellinint::mdim+m] * c2qqpzlog;
	  C2qqbp[idx] += faclog[i]*kernlog[i*mellinint::mdim+m] * c2qqbpzlog;

	  //cout << i << "  " << kern[i*mellinint::mdim+m] << endl;
	}
    }

  //Add delta pieces
  for (int m = 0; m < mellinint::mdim; m++)
    {
      int idx = anomalous::index(m,mesq::positive);
      C1qq[idx] += resconst::C1qqn;
      C1qqqed[idx] += resconst::C1qqn*(3./4.); // QED
      C2qq[idx] += resconst::C2qqn;
    }

  //Add plus-distribution pieces: the Mellin transform of 1/(1-z)+ is -S1(N-1)
  for (int m = 0; m < mellinint::mdim; m++)
    {
      complex <double> s1nm1 = cpsi0(mellinint::Np[m])+constants::euler;
      int idx = anomalous::index(m,mesq::positive);
      C2qq[idx] += (-s1nm1)*resconst::H2qqD0/2.;
    }
}

void ccoeff::truncate()
{
  //start from original untruncated moments
  for (int m = 0; m < mellinint::mdim; m++)
    {
      int idxp = anomalous::index(m,mesq::positive);
	  
      C1qg[idxp] = anomalous::C1QG[idxp]/2.;
      C1qq[idxp] = anomalous::C1QQ[idxp]/2.;
      C1qgamma[idxp] = anomalous::C1QG[idxp]/2.*6.;       // QED
      C1qqqed[idxp]  = anomalous::C1QQ[idxp]/2.*(3./4.);  // QED

      int idxm = anomalous::index(m,mesq::negative);
      C1qg[idxm]   = conj(C1qg[idxp]);
      C1qq[idxm]   = conj(C1qq[idxp]);
      C1qgamma[idxm]   = conj(C1qgamma[idxp]);  // QED
      C1qqqed[idxm]    = conj(C1qqqed[idxp]);   // QED
    }
  
  //Calculate truncated moments
  complex <double> C1qg_tr[mellinint::mdim];
  complex <double> C1qq_tr[mellinint::mdim];
  complex <double> C1qgamma_tr[mellinint::mdim];
  complex <double> C1qqqed_tr[mellinint::mdim];

  double x1 = phasespace::m/opts.sroot*exp(phasespace::ymin);
  double x2 = phasespace::m/opts.sroot*exp(-phasespace::ymax);

  //double x1 = 1e-8;//pow(phasespace::m/opts.sroot,2);
  //double x2 = 1e-8;//pow(phasespace::m/opts.sroot,2);

  double lx1 = log(x1);
  double lx2 = log(x2);
  
  //cache x^(N) values
  complex <double> x1n[mellinint::mdim];
  complex <double> x2n[mellinint::mdim];
  for (int n = 0; n < mellinint::mdim; n++)
    {
      x1n[n] = pow(x1,mellinint::Np[n]);
      x2n[n] = pow(x2,mellinint::Np[n]);
    }

  //Normalisation times Jacobian
  complex <double> facp = mellinint::CCp/2./M_PI/complex <double>(0.,1);
  complex <double> facm = mellinint::CCm/2./M_PI/complex <double>(0.,1);
  
  //original moments times prefactor and weight
  complex <double> c1qq_p[mellinint::mdim];
  complex <double> c1qg_p[mellinint::mdim];
  complex <double> c1qq_m[mellinint::mdim];
  complex <double> c1qg_m[mellinint::mdim];

  complex <double> c1qqqed_p[mellinint::mdim];   // QED
  complex <double> c1qgamma_p[mellinint::mdim];  // QED 
  complex <double> c1qqqed_m[mellinint::mdim];   // QED 
  complex <double> c1qgamma_m[mellinint::mdim];  // QED 
  for (int m = 0; m < mellinint::mdim; m++)
    {
      c1qq_p[m] = facp * C1qq[m] * mellinint::wn[m];
      c1qg_p[m] = facp * C1qg[m] * mellinint::wn[m];
      c1qq_m[m] = facm * conj(C1qq[m]) * mellinint::wn[m];
      c1qg_m[m] = facm * conj(C1qg[m]) * mellinint::wn[m];

      c1qqqed_p[m] = facp * C1qqqed[m] * mellinint::wn[m];
      c1qgamma_p[m] = facp * C1qgamma[m] * mellinint::wn[m];
      c1qqqed_m[m] = facm * conj(C1qqqed[m]) * mellinint::wn[m];
      c1qgamma_m[m] = facm * conj(C1qgamma[m]) * mellinint::wn[m];
    }

  //cache factor (1-x^(N-M))/(N-M) which limit is ln(x) when N-M -> 0
  complex <double> llx1p[mellinint::mdim][mellinint::mdim];
  complex <double> llx2p[mellinint::mdim][mellinint::mdim];
  complex <double> llx1m[mellinint::mdim][mellinint::mdim];
  complex <double> llx2m[mellinint::mdim][mellinint::mdim];
  for (int n = 0; n < mellinint::mdim; n++)
    for (int m = 0; m < mellinint::mdim; m++)
      {
	llx1p[n][m] = (1.-x1n[n]/x1n[m])/(mellinint::Np[n]-mellinint::Np[m]);
	llx2p[n][m] = (1.-x2n[n]/x2n[m])/(mellinint::Np[n]-mellinint::Np[m]);
	llx1m[n][m] = (1.-x1n[n]/conj(x1n[m]))/(mellinint::Np[n]-mellinint::Nm[m]);
	llx2m[n][m] = (1.-x2n[n]/conj(x2n[m]))/(mellinint::Np[n]-mellinint::Nm[m]);
      }

  //overwrite divergent diagonal part
  for (int n = 0; n < mellinint::mdim; n++)
    {
      llx1p[n][n] = -lx1;
      llx2p[n][n] = -lx2;
    }
  
  for (int n = 0; n < mellinint::mdim; n++)
    {
      //positive branch
      for (int m = 0; m < mellinint::mdim; m++)
	{
	    /*
	    if (m == n)
	      {
		ffx1[n][f] += fm1p[m][f] * (-lx1);
		ffx2[n][f] += fm2p[m][f] * (-lx2);
	      }
	    else
	      {
		ffx1[n][f] += fm1p[m][f] * (1.-x1n[n]/x1n[m])/(mellinint::Np[n]-mellinint::Np[m]);
		ffx2[n][f] += fm2p[m][f] * (1.-x2n[n]/x2n[m])/(mellinint::Np[n]-mellinint::Np[m]);
	      }
	    */
	    //	    ffx1[n][f] += fm1p[m][f]*llx1p[n][m]; 
	    //	    ffx2[n][f] += fm2p[m][f]*llx2p[n][m];
	  C1qq_tr[n] += c1qq_p[m]*llx1p[n][m] - c1qq_m[m]*llx1m[n][m]; 
	  C1qg_tr[n] += c1qg_p[m]*llx2p[n][m] - c1qg_m[m]*llx2m[n][m];

	  C1qqqed_tr[n]  += c1qqqed_p[m]*llx1p[n][m] - c1qqqed_m[m]*llx1m[n][m];   // QED
	  C1qgamma_tr[n] += c1qgamma_p[m]*llx2p[n][m] - c1qgamma_m[m]*llx2m[n][m];   // QED
	}

      /*
      //negative branch
      for (int m = 0; m < mellinint::mdim; m++)
	for (int f = 0; f < 2*MAXNF+1; f++)
	  {
	    //	    ffx1[n][f] -= fm1m[m][f]*(1.-x1n[n]/conj(x1n[m]))/(mellinint::Np[n]-mellinint::Nm[m]);
	    //	    ffx2[n][f] -= fm2m[m][f]*(1.-x2n[n]/conj(x2n[m]))/(mellinint::Np[n]-mellinint::Nm[m]);
	    ffx1[n][f] -= fm1m[m][f]*llx1m[n][m];
	    ffx2[n][f] -= fm2m[m][f]*llx2m[n][m];
	  }
      */

  
      //cout << "truncated " << n << ffx1[5] << endl;
    }

  //replace moments
  for (int m = 0; m < mellinint::mdim; m++)
    {
      int idxp = anomalous::index(m,mesq::positive);
	  
      C1qg[idxp] = C1qg_tr[m];
      C1qq[idxp] = C1qq_tr[m];

      C1qgamma[idxp] = C1qgamma_tr[m];   // QED
      C1qqqed[idxp] = C1qqqed_tr[m];   // QED

      int idxm = anomalous::index(m,mesq::negative);
      C1qg[idxm]   = conj(C1qg[idxp]);
      C1qq[idxm]   = conj(C1qq[idxp]);

      C1qgamma[idxm]   = conj(C1qgamma[idxp]);
      C1qqqed[idxm]   = conj(C1qqqed[idxp]);
    }

}

void ccoeff::free()
{
  if (opts.mellin1d)
    {
      delete[] C1qg;
      delete[] C1qq;
      delete[] C1qgamma;
      delete[] C1qqqed;
      delete[] C1qqb;
      delete[] C1qqp;
      delete[] C1qqbp;
      delete[] C1gq;
      delete[] C2qg;
      delete[] C2qq;
      delete[] C2qqb;
      delete[] C2qqp;
      delete[] C2qqbp;
      delete[] C2gq;
      delete[] C3qg;
      delete[] C3qq;
      delete[] C3qqb;
      delete[] C3qqp;
      delete[] C3qqbp;
      delete[] C4qg;
      delete[] C4qq;
      delete[] C4qqb;
      delete[] C4qqp;
      delete[] C4qqbp;
      delete[] C1qg_reg;
    }
  else
    {
      delete[] C1qg_1;
      delete[] C1qq_1;
      delete[] C1qgamma_1;   // QED
      delete[] C1qqqed_1;    // QED
      delete[] C1qqb_1;
      delete[] C1qqp_1;
      delete[] C1qqbp_1;
      delete[] C2qg_1;
      delete[] C2qq_1;
      delete[] C2qqb_1;
      delete[] C2qqp_1;
      delete[] C2qqbp_1;
      delete[] C3qg_1;
      delete[] C3qq_1;
      delete[] C3qqb_1;
      delete[] C3qqp_1;
      delete[] C3qqbp_1;
      delete[] C4qg_1;
      delete[] C4qq_1;
      delete[] C4qqb_1;
      delete[] C4qqp_1;
      delete[] C4qqbp_1;

      delete[] C1qg_2;
      delete[] C1qq_2;
      delete[] C1qgamma_2;    // QED
      delete[] C1qqqed_2;     // QED
      delete[] C1qqb_2;
      delete[] C1qqp_2;
      delete[] C1qqbp_2;
      delete[] C2qg_2;
      delete[] C2qq_2;
      delete[] C2qqb_2;
      delete[] C2qqp_2;
      delete[] C2qqbp_2;
      delete[] C3qg_2;
      delete[] C3qq_2;
      delete[] C3qqb_2;
      delete[] C3qqp_2;
      delete[] C3qqbp_2;
      delete[] C4qg_2;
      delete[] C4qq_2;
      delete[] C4qqb_2;
      delete[] C4qqp_2;
      delete[] C4qqbp_2;
      delete[] C1qg_reg_1;
      delete[] C1qg_reg_2;
    }
}
void ccoeff::release()
{
  //  if (opts.melup <= 1)
  //    free();
  
  icoeff::release();
}
