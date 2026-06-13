#include "expc.h"
#include "resconst.h"
#include "mellinint.h"
#include "scales.h"
#include "blim.h"
#include "pdfevol.h"
#include "anomalous.h"
#include "hcoefficients.h"
#include "mesq.h"
#include "resint.h"
#include "settings.h"
#include "sudakovff.h"
#include "ccoeff.h"
#include "alphas.h"
#include "gint.h"
#include "pdf.h"
#include "Li2.hpp"
#include "isnan.h"
#include "psi.h"
#include "constants.h"
#include "phasespace.h"
#include "hcoeff.h"
#include "bstar.h"
#include "btilde.h"
#include <iostream>

#define cdefexp(c1,c2,c3) C1 = c1; C2 = c2; C3 = c3; rootC = sqrt(pow(C2/C1,2)-4.*C3/C1); n2llaexp = (1.+aass*C2/C1/(1.-lam))/(1.+aass*C2/C1); //n2llaexp = (1.-lam+aass*C2/C1)/(1.-lam)/(1.+aass*C2/C1);
#define n3llaexp() n2llaexp\
			*exp(log((aass*(rootC+C2/C1)+2.)*(aass*(rootC-C2/C1)-2.*(1.-lam))/(aass*(rootC+C2/C1)+2.*(1.-lam))/(aass*(rootC-C2/C1)-2.))*C2/C1/rootC) \
			*(pow(1.-lam,2)+aass*C2/C1*(1.-lam)+aass2*C3/C1)/pow(1.-lam,2)/(1.+aass*C2/C1+aass2*C3/C1) \
			* exp(-log(n2llaexp)*beta1/beta0*C1/C2) *pow(aexpB,beta1/beta0) \
			*exp(+2.*aass*beta1/beta0*(lam+log(1.-lam))/(lam-1.)-C1/C2*beta1/beta0 * log(1.-lam) *(log(1.-lam)+2.*log( aass*C2/C1/(1.-lam + aass*C2/C1)))+2.*C1/C2*beta1/beta0*(-Li2(-C1/C2/aass)+Li2(C1/C2*(lam-1.)/aass))) \
			*exp(log(n2llaexp)*2.*beta0*(C1/C2)*(LR-LQ))	\
			*pow(aexpB,-2.*beta0*(LR-LQ)) ;

//*exp(log(n2llaexp)*2.*beta0*(C1/C2)*(LQR)) \
//*pow(aexpB,-2.*beta0*(LQR)) ;

using namespace polylogarithm;

complex <double> *expc::aexpgg;
complex <double> *expc::aexpqq;
complex <double> *expc::aexpqg;
complex <double> *expc::aexpqqqed;
complex <double> *expc::aexpqgamma;
complex <double> *expc::aexpqqb;
complex <double> *expc::aexpqqp;
complex <double> *expc::aexpqqbp;
complex <double> *expc::aexpqq_1;
complex <double> *expc::aexpqg_1;
complex <double> *expc::aexpgq_1;
complex <double> *expc::aexpqqqed_1;
complex <double> *expc::aexpqgamma_1;
complex <double> *expc::aexpqqb_1;
complex <double> *expc::aexpqqp_1;
complex <double> *expc::aexpqqbp_1;
complex <double> *expc::aexpqq_2;
complex <double> *expc::aexpqg_2;
complex <double> *expc::aexpgq_2;
complex <double> *expc::aexpqqqed_2;
complex <double> *expc::aexpqgamma_2;
complex <double> *expc::aexpqqb_2;
complex <double> *expc::aexpqqp_2;
complex <double> *expc::aexpqqbp_2;

complex <double> *expc::aexpqg_bb;
complex <double> *expc::aexpqg_1_bb;
complex <double> *expc::aexpqg_2_bb;

complex <double> *expc::qqb;
complex <double> *expc::qg;
complex <double> *expc::qqbqed;
complex <double> *expc::qgamma;
complex <double> *expc::qg_1;
complex <double> *expc::qg_2;
complex <double> *expc::gq_1;
complex <double> *expc::gq_2;
complex <double> *expc::qgamma_1;
complex <double> *expc::qgamma_2;
complex <double> *expc::qq;
complex <double> *expc::qq_1;
complex <double> *expc::qq_2;
complex <double> *expc::qqp;
complex <double> *expc::qqp_1;
complex <double> *expc::qqp_2;
complex <double> *expc::qqbp;
complex <double> *expc::qqbp_1;
complex <double> *expc::qqbp_2;
complex <double> *expc::gg;
complex <double> *expc::qbg;
complex <double> *expc::qbg_1;
complex <double> *expc::qbg_2;
complex <double> *expc::qpg;
complex <double> *expc::qpg_1;
complex <double> *expc::qpg_2;
complex <double> *expc::qbpg;
complex <double> *expc::qbpg_1;
complex <double> *expc::qbpg_2;
complex <double> *expc::qbq;
complex <double> *expc::qbqp;
complex <double> *expc::qbqp_1;
complex <double> *expc::qbqp_2;
complex <double> *expc::qbqbp;
complex <double> *expc::qbqbp_1;
complex <double> *expc::qbqbp_2;
complex <double> *expc::qpqp;
complex <double> *expc::qpqp_1;
complex <double> *expc::qpqp_2;
complex <double> *expc::qpqbp;
complex <double> *expc::qbpqp;

complex <double> *expc::gg_bb;

//complex <double> expc::aexp;
complex <double> expc::xlambda;
complex <double> expc::lamB;
complex <double> expc::xlambdaqed;
complex <double> expc::lamBqed;
complex <double> expc::lamC;
complex <double> expc::lamD;
complex <double> expc::lamE;
complex <double> expc::lamF;
complex <double> expc::lamG;
complex <double> expc::lamH;
complex <double> expc::aexpB;
complex <double> expc::aexpBqed;
complex <double> expc::aexpC;
complex <double> expc::aexpD;
complex <double> expc::aexpE;
complex <double> expc::aexpF;
complex <double> expc::aexpG;
complex <double> expc::aexpH;

complex <double> expc::lamB_bb;
complex <double> expc::aexpB_bb;

using namespace anomalous;
using namespace resconst;
using namespace resint;
using namespace ccoeff;

void expc::allocate()
{
  //if (opts.order == 0)
  //return;

  //allocate memory
  if (opts.mellin1d)
    {
      aexpqq   = new complex <double> [mellinint::mdim*2];
      aexpqg   = new complex <double> [mellinint::mdim*2];
      aexpqqqed   = new complex <double> [mellinint::mdim*2];
      aexpqgamma   = new complex <double> [mellinint::mdim*2];
      aexpqqb  = new complex <double> [mellinint::mdim*2];
      aexpqqp  = new complex <double> [mellinint::mdim*2];
      aexpqqbp = new complex <double> [mellinint::mdim*2];

      aexpqg_bb   = new complex <double> [mellinint::mdim*2];
      
      qqb  = new complex <double> [mellinint::mdim*2];
      qg   = new complex <double> [mellinint::mdim*2];
      qqbqed  = new complex <double> [mellinint::mdim*2];
      qgamma   = new complex <double> [mellinint::mdim*2];
      qq   = new complex <double> [mellinint::mdim*2];
      qqp  = new complex <double> [mellinint::mdim*2];
      qqbp = new complex <double> [mellinint::mdim*2];
      gg   = new complex <double> [mellinint::mdim*2];
      qbg  = new complex <double> [mellinint::mdim*2];
      qpg  = new complex <double> [mellinint::mdim*2];
      qbpg = new complex <double> [mellinint::mdim*2];
      qbq   = new complex <double> [mellinint::mdim*2];
      qbqp  = new complex <double> [mellinint::mdim*2];
      qbqbp = new complex <double> [mellinint::mdim*2];
      qpqp  = new complex <double> [mellinint::mdim*2];
      qpqbp = new complex <double> [mellinint::mdim*2];
      qbpqp = new complex <double> [mellinint::mdim*2];

      gg_bb   = new complex <double> [mellinint::mdim*2];
    }
  else
    {
	  aexpgg   = new complex <double> [mellinint::mdim*2];
      aexpqq_1   = new complex <double> [mellinint::mdim*2];
      aexpqg_1   = new complex <double> [mellinint::mdim*2];
	  aexpgq_1   = new complex <double> [mellinint::mdim*2];
      aexpqqqed_1   = new complex <double> [mellinint::mdim*2];
      aexpqgamma_1   = new complex <double> [mellinint::mdim*2];
      aexpqqb_1  = new complex <double> [mellinint::mdim*2];
      aexpqqp_1  = new complex <double> [mellinint::mdim*2];
      aexpqqbp_1 = new complex <double> [mellinint::mdim*2];
      aexpqq_2   = new complex <double> [mellinint::mdim*2];
      aexpqg_2   = new complex <double> [mellinint::mdim*2];
      aexpgq_2   = new complex <double> [mellinint::mdim*2];
      aexpqqqed_2   = new complex <double> [mellinint::mdim*2];
      aexpqgamma_2   = new complex <double> [mellinint::mdim*2];
      aexpqqb_2  = new complex <double> [mellinint::mdim*2];
      aexpqqp_2  = new complex <double> [mellinint::mdim*2];
      aexpqqbp_2 = new complex <double> [mellinint::mdim*2];

      aexpqg_1_bb   = new complex <double> [mellinint::mdim*2];
      aexpqg_2_bb   = new complex <double> [mellinint::mdim*2];
      
      qqb    = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qg_1   = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qg_2   = new complex <double> [mellinint::mdim*mellinint::mdim*2];
	  gq_1   = new complex <double> [mellinint::mdim*mellinint::mdim*2];
	  gq_2   = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qqbqed    = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qgamma_1   = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qgamma_2   = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qq_1   = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qq_2   = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qqp_1  = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qqp_2  = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qqbp_1 = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qqbp_2 = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      gg     = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qbg_1  = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qbg_2  = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qpg_1  = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qpg_2  = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qbpg_1 = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qbpg_2 = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qbq     = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qbqp_1  = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qbqp_2  = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qbqbp_1 = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qbqbp_2 = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qpqp_1  = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qpqp_2  = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qpqbp   = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qbpqp   = new complex <double> [mellinint::mdim*mellinint::mdim*2];

      gg_bb     = new complex <double> [mellinint::mdim*mellinint::mdim*2];
    }
}

void expc::free()
{
  //if (opts.order == 0)
  //return;

  if (opts.mellin1d)
    {
      delete[] aexpqq;
      delete[] aexpqg;
      delete[] aexpqqqed;
      delete[] aexpqgamma;
      delete[] aexpqqb;
      delete[] aexpqqp;
      delete[] aexpqqbp;
      delete[] qqb ;
      delete[] qg  ;
      delete[] qqbqed ;
      delete[] qgamma  ;
      delete[] qq  ;
      delete[] qqp ;
      delete[] qqbp;
      delete[] gg  ;
      delete[] qbg ;
      delete[] qpg ;
      delete[] qbpg;
      delete[] qbq;
      delete[] qbqp;
      delete[] qbqbp;
      delete[] qpqp;
      delete[] qpqbp;
      delete[] qbpqp;
      delete[] aexpqg_bb;
      delete[] gg_bb  ;
    }
  else
    {
	  delete[] aexpgg;
	  delete[] aexpgq_1;
	  delete[] aexpgq_2;
      delete[] aexpqq_1;
      delete[] aexpqg_1;
      delete[] aexpqqqed_1;
      delete[] aexpqgamma_1;
      delete[] aexpqqb_1;
      delete[] aexpqqp_1;
      delete[] aexpqqbp_1;
      delete[] aexpqq_2;
      delete[] aexpqg_2;
      delete[] aexpqqqed_2;
      delete[] aexpqgamma_2;
      delete[] aexpqqb_2;
      delete[] aexpqqp_2;
      delete[] aexpqqbp_2;
      delete[] qqb   ;
      delete[] qg_1  ;
      delete[] qg_2  ;
	  delete[] gq_1  ;
	  delete[] gq_2  ;
      delete[] qqbqed   ;
      delete[] qgamma_1  ;
      delete[] qgamma_2  ;
      delete[] qq_1  ;
      delete[] qq_2  ;
      delete[] qqp_1 ;
      delete[] qqp_2 ;
      delete[] qqbp_1;
      delete[] qqbp_2;
      delete[] gg    ;
      delete[] qbg_1 ;
      delete[] qbg_2 ;
      delete[] qpg_1 ;
      delete[] qpg_2 ;
      delete[] qbpg_1;
      delete[] qbpg_2;
      delete[] qbq;
      delete[] qbqp_1;
      delete[] qbqp_2;
      delete[] qbqbp_1;
      delete[] qbqbp_2;
      delete[] qpqp_1;
      delete[] qpqp_2;
      delete[] qpqbp;
      delete[] qbpqp;
      delete[] aexpqg_1_bb;
      delete[] aexpqg_2_bb;
      delete[] gg_bb    ;
    }
}

void expc::reset()
{
  //if (opts.order_expc == 0)
  //return;

  if (opts.mellin1d)
    if (opts.sumlogs)
      {
	fill(qqb , qqb+mellinint::mdim*2, sudakov::S);
	fill(qg  , qg +mellinint::mdim*2, sudakov::S);
	fill(qqbqed , qqbqed+mellinint::mdim*2, sudakov::S);
	fill(qgamma  , qgamma +mellinint::mdim*2, sudakov::S);
	fill(qq  , qq +mellinint::mdim*2, sudakov::S);
	fill(qqp , qqp +mellinint::mdim*2, sudakov::S);
	fill(qqbp, qqbp +mellinint::mdim*2, sudakov::S);
	fill(gg  , gg +mellinint::mdim*2, sudakov::S);
	fill(qbg , qbg +mellinint::mdim*2, sudakov::S);
	fill(qpg , qpg +mellinint::mdim*2, sudakov::S);
	fill(qbpg, qbpg +mellinint::mdim*2, sudakov::S);
      }
    else
      {
	fill(aexpqq   ,aexpqq   +mellinint::mdim*2, 1.);
	fill(aexpqg   ,aexpqg   +mellinint::mdim*2, 1.);
	fill(aexpqqqed   ,aexpqqqed   +mellinint::mdim*2, 1.);
	fill(aexpqgamma   ,aexpqgamma   +mellinint::mdim*2, 1.);
	fill(aexpqq   ,aexpqq   +mellinint::mdim*2, 1.);
	fill(aexpqg   ,aexpqg   +mellinint::mdim*2, 1.);
	fill(aexpqqb  ,aexpqqb  +mellinint::mdim*2, 1.);
	fill(aexpqqp  ,aexpqqp  +mellinint::mdim*2, 1.);
	fill(aexpqqbp ,aexpqqbp +mellinint::mdim*2, 1.);
	fill(qqb , qqb+mellinint::mdim*2, 1.);
	fill(qg  , qg +mellinint::mdim*2, 1.);
	fill(qqbqed , qqbqed+mellinint::mdim*2, 1.);
	fill(qgamma  , qgamma +mellinint::mdim*2, 1.);
	fill(qq  , qq +mellinint::mdim*2, 1.);
	fill(qqp , qqp +mellinint::mdim*2, 1.);
	fill(qqbp, qqbp +mellinint::mdim*2, 1.);
	fill(gg  , gg +mellinint::mdim*2, 1.);
	fill(qbg , qbg +mellinint::mdim*2, 1.);
	fill(qpg , qpg +mellinint::mdim*2, 1.);
	fill(qbpg, qbpg +mellinint::mdim*2, 1.);
	fill(qbq  , qbq   + mellinint::mdim*2, 1.);
	fill(qbqp , qbqp  + mellinint::mdim*2, 1.);
	fill(qbqbp, qbqbp + mellinint::mdim*2, 1.);
	fill(qpqp , qpqp  + mellinint::mdim*2, 1.);
	fill(qpqbp, qpqbp + mellinint::mdim*2, 1.);
	fill(qbpqp, qbpqp + mellinint::mdim*2, 1.);
	fill(aexpqg_bb   ,aexpqg_bb   +mellinint::mdim*2, 1.);
	fill(gg_bb  , gg_bb +mellinint::mdim*2, 1.);
      }
  else
    if (opts.sumlogs)
      {
	fill(qqb,    qqb +mellinint::mdim*mellinint::mdim*2, sudakov::S);
	fill(qg_1,   qg_1+mellinint::mdim*mellinint::mdim*2, sudakov::S);
	fill(qg_2,   qg_2+mellinint::mdim*mellinint::mdim*2, sudakov::S);
	fill(qqbqed,    qqbqed +mellinint::mdim*mellinint::mdim*2, sudakov::S);
	fill(qgamma_1,   qgamma_1+mellinint::mdim*mellinint::mdim*2, sudakov::S);
	fill(qgamma_2,   qgamma_2+mellinint::mdim*mellinint::mdim*2, sudakov::S);
	fill(qq_1,   qq_1+mellinint::mdim*mellinint::mdim*2, sudakov::S);
	fill(qq_2,   qq_2+mellinint::mdim*mellinint::mdim*2, sudakov::S);
	fill(qqp_1,  qqp_1+mellinint::mdim*mellinint::mdim*2, sudakov::S);
	fill(qqp_2,  qqp_2+mellinint::mdim*mellinint::mdim*2, sudakov::S);
	fill(qqbp_1, qqbp_1+mellinint::mdim*mellinint::mdim*2, sudakov::S);
	fill(qqbp_2, qqbp_2+mellinint::mdim*mellinint::mdim*2, sudakov::S);
	fill(gg,     gg  +mellinint::mdim*mellinint::mdim*2, sudakov::S);

      }
    else
      {
	fill(aexpgg   ,aexpgg   +mellinint::mdim*2, 1.);
	fill(aexpqq_1   ,aexpqq_1   +mellinint::mdim*2, 1.);
	fill(aexpqg_1   ,aexpqg_1   +mellinint::mdim*2, 1.);
	fill(aexpgq_1   ,aexpgq_1   +mellinint::mdim*2, 1.);
	fill(aexpqqqed_1   ,aexpqqqed_1   +mellinint::mdim*2, 1.);
	fill(aexpqgamma_1   ,aexpqgamma_1   +mellinint::mdim*2, 1.);
	fill(aexpqqb_1  ,aexpqqb_1  +mellinint::mdim*2, 1.);
	fill(aexpqqp_1  ,aexpqqp_1  +mellinint::mdim*2, 1.);
	fill(aexpqqbp_1 ,aexpqqbp_1 +mellinint::mdim*2, 1.);
	fill(aexpqq_2   ,aexpqq_2   +mellinint::mdim*2, 1.);
	fill(aexpqg_2   ,aexpqg_2   +mellinint::mdim*2, 1.);
	fill(aexpgq_2   ,aexpgq_2   +mellinint::mdim*2, 1.);
	fill(aexpqqqed_2   ,aexpqqqed_2   +mellinint::mdim*2, 1.);
	fill(aexpqgamma_2   ,aexpqgamma_2   +mellinint::mdim*2, 1.);
	fill(aexpqqb_2  ,aexpqqb_2  +mellinint::mdim*2, 1.);
	fill(aexpqqp_2  ,aexpqqp_2  +mellinint::mdim*2, 1.);
	fill(aexpqqbp_2 ,aexpqqbp_2 +mellinint::mdim*2, 1.);
	fill(qqb,    qqb +mellinint::mdim*mellinint::mdim*2, 1.);
	fill(qg_1,   qg_1+mellinint::mdim*mellinint::mdim*2, 1.);
	fill(qg_2,   qg_2+mellinint::mdim*mellinint::mdim*2, 1.);
	fill(gq_1,   gq_1+mellinint::mdim*mellinint::mdim*2, 1.);
	fill(gq_2,   gq_2+mellinint::mdim*mellinint::mdim*2, 1.);
	fill(qqbqed,    qqbqed +mellinint::mdim*mellinint::mdim*2, 1.);
	fill(qgamma_1,   qgamma_1+mellinint::mdim*mellinint::mdim*2, 1.);
	fill(qgamma_2,   qgamma_2+mellinint::mdim*mellinint::mdim*2, 1.);
	fill(qq_1,   qq_1+mellinint::mdim*mellinint::mdim*2, 1.);
	fill(qq_2,   qq_2+mellinint::mdim*mellinint::mdim*2, 1.);
	fill(qqp_1,  qqp_1+mellinint::mdim*mellinint::mdim*2, 1.);
	fill(qqp_2,  qqp_2+mellinint::mdim*mellinint::mdim*2, 1.);
	fill(qqbp_1, qqbp_1+mellinint::mdim*mellinint::mdim*2, 1.);
	fill(qqbp_2, qqbp_2+mellinint::mdim*mellinint::mdim*2, 1.);
	fill(gg,     gg  +mellinint::mdim*mellinint::mdim*2, 1.);
	fill(qbq    , qbq     + mellinint::mdim*mellinint::mdim*2, 1.);
	fill(qbqp_1 , qbqp_1  + mellinint::mdim*mellinint::mdim*2, 1.);
	fill(qbqp_2 , qbqp_2  + mellinint::mdim*mellinint::mdim*2, 1.);
	fill(qbqbp_1, qbqbp_1 + mellinint::mdim*mellinint::mdim*2, 1.);
	fill(qbqbp_2, qbqbp_2 + mellinint::mdim*mellinint::mdim*2, 1.);
	fill(qpqp_1 , qpqp_1  + mellinint::mdim*mellinint::mdim*2, 1.);
	fill(qpqp_2 , qpqp_2  + mellinint::mdim*mellinint::mdim*2, 1.);
	fill(qpqbp  , qpqbp   + mellinint::mdim*mellinint::mdim*2, 1.);
	fill(qbpqp  , qbpqp   + mellinint::mdim*mellinint::mdim*2, 1.);
	fill(aexpqg_1_bb   ,aexpqg_1_bb   +mellinint::mdim*2, 1.);
	fill(aexpqg_2_bb   ,aexpqg_2_bb   +mellinint::mdim*2, 1.);
	fill(gg_bb,     gg_bb  +mellinint::mdim*mellinint::mdim*2, 1.);
      }
}

//b-dependent C exponentiation
void expc::calc(complex <double> b)
{

	if (opts.nproc == 4) {

		//Return 1 if b is close to the Landau pole
		double b_L;
		if (!opts.modlog)
			b_L = b0/scales::res * exp(1./(2.*aass*beta0));
		else if (opts.p == 1)
			b_L = b0/scales::res * sqrt(exp(1./(aass*beta0))-1.);
		else
			b_L = b0/scales::res * pow(sqrt(exp(opts.p/(aass*beta0))-1.),1./opts.p);

		if (fabs(b-b_L) < 1e-8)
			return;

		//  if (opts.expc == 0)
		//    return;

		//LL (no evolution)
		if (opts.order_expc == 0)
			return;

		double blim = blim::expc;  //better use blim::sudakov?

		//  the limit below implies xlambda<1/2 and then aa2<= 1
		//      blim=b0p*(1/q)*exp(1/(2*as*beta0))
		//      blim=b0p*(1/q)*exp(1/(4*as*beta0))
		// Set a limit to avoid very large values of b (= very small scales ~1/b)
		//       blim=b0p*(1/q)*exp(1/(2*aass*beta0)) ! avoid Landau pole
		//     write(*,*) "blim",blim
		//     without this additional blim some scale variations will fail (when mures > muren)
		//  blim=0.5; // --> allow this to a separate blim in the settings, or set it using muren instead of mures
		//      blim=1.1229190d0
		//      blim=3d0

		//Set b according to bstar or other prescriptions
		complex <double> bstar;
		if (opts.bprescription == 0 || opts.bprescription == 4 || opts.bstar_expc)
			bstar = real(b)/sqrt(1.+pow(real(b)/blim,2));
		else
			bstar = b;

		double Q = scales::res;
		double aass2 = pow(aass,2);

		complex <double> blog;
		if (!opts.modlog)
			blog = log(pow(Q*bstar/b0,2));   //normal sudakov
			else if (opts.p == 1)
				blog = log(pow(Q*bstar/b0,2) + 1.); //modified sudakov
				else
					blog = 1./opts.p*log(pow(Q*bstar/b0,2*opts.p) + 1.); //modified sudakov with exponent p

					// cout << "\nmodlog = " << opts.modlog << endl;
					// cout << "p = " << opts.p << endl;
					// cout << "bstar = " << bstar << endl;
					// cout << "blim = " << blim << endl;
					// cout << "b0 = " << b0 << endl;
					// cout << "blog = " << blog << endl;

					//if (opts.modlog)
					//  blog = log(pow(Q*bstar/b0,2) + 1.); //modified sudakov
					//else
					//  blog = log(pow(Q*bstar/b0,2));   //normal sudakov

					//lambda as defined in Eq. (25) of hep-ph/0508068.
					complex <double> xlambda = beta0*aass*blog;
					complex <double> log1xlambda = log(1.-xlambda);

					complex <double> logasl;

					//LL (no evolution)
					logasl = 0.;

					//NLL (LO evolution)
					if (opts.order_expc >= 1)
						logasl += log1xlambda;

		//NNLL (NLO evolution)
		if (opts.order_expc >= 2)
			logasl += aass* beta1/beta0*log1xlambda/(1.-xlambda);

		//NNNLL (NNLO evolution)
		if (opts.order_expc >= 3)
			logasl += aass2* ((pow(beta1/beta0,2)-beta2/beta0) *xlambda/pow(1.-xlambda,2)
			+ pow(beta1/beta0,2)             *log1xlambda/pow(1.-xlambda,2)
			- pow(beta1/beta0,2)             *pow(log1xlambda,2)/(2.*pow(1.-xlambda,2)));

		double LQR = real(logq2mur2-2.*loga);
		//QCD coupling scale dependence
		if (opts.order_expc >= 2)
			logasl += aass*LQR*beta0*xlambda/(1.-xlambda);

		//from c3new.txt
		if (opts.order_expc >= 3)
			logasl += aass2*(+LQR*beta1                   *(xlambda-log1xlambda)/pow(1.-xlambda,2)
			+LQR*beta1                   *xlambda/(1.-xlambda)                      //missing piece
			+0.5*pow(LQR,2)*pow(beta0,2) *xlambda*(xlambda-2.)/pow(1.-xlambda,2));  //missing piece

			//from g4 LQR part by replacing B1qbar->beta0 and B2qbar->beta1 (and all An pieces set to zero)
			//if (opts.order_expc >= 3)
			//  logasl += aass2*(LQR*beta1*(xlambda*(2.-xlambda)-log1xlambda)/pow(1.-xlambda,2)
			//		     +0.5*pow(LQR,2)*pow(beta0,2)*xlambda*(xlambda-2.)/pow(1.-xlambda,2));

			//cout << logasl << " logasl" << endl;
			//!!! gint::logasl is evaluated with blim::sudakov, not with blim:expc !!!
			if (opts.numexpc)
				logasl = gint::logasl_expc;
		//cout << gint::logasl_expc << " gint logasl" << endl;

		if (isnan_ofast(real(logasl)) || isnan_ofast(imag(logasl)))
		{
			cout << "Warning, logasl = " << logasl << ", b = "  << b << ", lambda = " << xlambda
			<< " log(1-xlambda) " << log1xlambda
			<< endl;
			logasl = 0.;
		}

		//For the aexp evolution operator use the same hard scale as for the PDF evolution, and the same blim
		//aexp should always reflect the PDF evolution, in evolmode 2 and 3 it should be log(ASF/ASI), plus eventually the LQR scale variation pieces
		//cout << "aexp " << logasl << " pdfs " << pdfevol::logasl << endl;
		//--> No aexp should be with Q and not muf

		//if (opts.mufevol)
		//complex <double> logasl = pdfevol::logasl;

		complex <double> aexp = exp(-logasl); //aexp is approximately alphas(Q^2)/alphas(mub^2)
		expc::aexpH = aexp; // To be used in spin correlations for process - Higgs

		//aexp = exp(-logasl); //aexp is approximately alphas(Q^2)/alphas(mub^2)

		if (isnan_ofast(real(aexp)) || isnan_ofast(imag(aexp)))
		{
			cout << "Warning, aexp = " << aexp << ", b = "  << b << ", lambda = " << xlambda
			<< " log(1-xlambda) " << log1xlambda
			<< endl;
			aexp = 1.;
		}

		/*
		 * //Set b according to bstar or other prescriptions
		 * //For the C exponentiation always use Q at the resummation scale, and blim::sudakov, to match what is done in the Sudakov
		 * blim = blim::expc; //better use blim::sudakov;
		 * complex <double> bstar;
		 * if (opts.bprescription == 0 || opts.bstar_expc) //better use opts.bstar_sudakov)
		 *	bstar = real(b)/sqrt(1.+pow(real(b)/blim,2));
		 * else
		 *	bstar = b;
		 *
		 * double Q = scales::res;
		 * complex <double> blog;
		 * if (opts.modlog)
		 *	blog = log(pow(Q*bstar/b0,2) + 1.); //modified sudakov
		 * else
		 *	blog = log(pow(Q*bstar/b0,2));   //normal sudakov
		 *
		 * //lambda as defined in Eq. (25) of hep-ph/0508068.
		 * complex <double> xlambda = beta0*aass*blog;
		 * complex <double> log1xlambda = log(1.-xlambda);
		 *
		 * Q = scales::res;
		 * if (opts.modlog)
		 *	blog = log(pow(Q*bstar/b0,2) + 1.); //modified sudakov
		 * else
		 *	blog = log(pow(Q*bstar/b0,2));   //normal sudakov
		 *
		 * xlambda = beta0*aass*blog;
		 * log1xlambda = log(1.-xlambda);
		 */

		//Multiply lamB, lamC, lamD by the corresponding power of alphas to avoid overflow/uunderflow problems in the exponential
		complex <double> lamB = aass*xlambda/(1.-xlambda);
		complex <double> aexpB = exp(lamB);

		// cout << "beta0 = " << beta0 << endl;
		// cout << "aass = " << aass << endl;
		// cout << "xlambda = " << xlambda << endl;
		// cout << "aexpB = " << aexpB << endl;

		complex <double> lamC = aass2*xlambda*(xlambda-2.)/pow(1.-xlambda,2);
		complex <double> aexpC = exp(lamC);

		complex <double> lamD = aass2*log1xlambda/pow(1.-xlambda,2);
		complex <double> aexpD = exp(lamD);


		//  if (opts.bprescription == 0)
		//    {
		//      //double lambdaqcd         = muren/(exp(1./(2.*aass*beta0))); //--> corresponds to a divergence in alphas
		//      double lambdaqcd     = mures/(exp(1./(2.*aass*beta0))); //--> correspond to a divergence in the Sudakov
		//      blim = b0/lambdaqcd;
		//      bstar = b/sqrt(1.+pow(b/blim,2));
		//      if (opts.modlog)
		//	blog = log(pow(Q*bstar/b0,2) + 1.); //modified sudakov
		//      else
		//	blog = log(pow(Q*bstar/b0,2));   //normal sudakov
		//
		//      xlambda = beta0*aass*blog;
		//      log1xlambda = log(1.-xlambda);
		//
		//      lamC = xlambda*(xlambda-2.)/pow(1.-xlambda,2);
		//      aexpC = exp(lamC);
		//    }
		//

		double aassh = aass/2.;
		double aasshsq = pow(aassh,2);

		complex <double> aexp2 = pow(aexp,2);

		// NLL Exponentiation of C coefficients .......................................................................................................................................................

		// CJU
		// cout << "order expc = " << opts.order_expc << endl;
		// cout << "expc = " << opts.expc << endl;

		if (opts.order_expc == 1) {

			if (opts.mellin1d) {

				if (opts.sumlogs) {

					fill(qg, qg+mellinint::mdim*2, exp(-logasl+sudakov::logS));
				}

				else {

					fill(aexpqg, aexpqg+mellinint::mdim*2, aexp);
					fill(qg, qg+mellinint::mdim*2, aexp);

				}
			}

			// Mellin 2d
			else {

				if (opts.sumlogs) {

					fill(qg_1, qg_1+mellinint::mdim*mellinint::mdim*2, exp(-logasl+sudakov::logS));
					fill(qg_2, qg_2+mellinint::mdim*mellinint::mdim*2, exp(-logasl+sudakov::logS));

				}

				else {

					fill(aexpqg_1, aexpqg_1+mellinint::mdim*2, aexp);
					fill(aexpqg_2, aexpqg_2+mellinint::mdim*2, aexp);
					fill(aexpgq_1, aexpgq_1+mellinint::mdim*2, aexp);
					fill(aexpgq_2, aexpgq_2+mellinint::mdim*2, aexp);
					fill(qg_1, qg_1+mellinint::mdim*mellinint::mdim*2, aexp);
					fill(qg_2, qg_2+mellinint::mdim*mellinint::mdim*2, aexp);
					fill(gq_1, gq_1+mellinint::mdim*mellinint::mdim*2, aexp);
					fill(gq_2, gq_2+mellinint::mdim*mellinint::mdim*2, aexp);

				}

			}

			//gint::calc(b);
			//cout << endl;
			//cout << " gint::logasl "   << gint::logasl      << " expc " << logasl      << endl;
			// NNLL
			//complex <double> c1delta = pow(aexpb,aassh*(-2.*C1qqn));

		}

		// NNLL Exponentiation of C coefficients ......................................................................................................................................................

		// CJU
		// cout << "order expc = " << opts.order_expc << endl;
		// cout << "expc = " << opts.expc << endl;

		if (opts.order_expc == 2) {

			if (opts.mellin1d) {

				//Use pointers to allocate memory on the heap
				complex <double> *alogqq   = new complex <double> [mellinint::mdim*2];
				complex <double> *alogqg   = new complex <double> [mellinint::mdim*2];
				complex <double> *alogqqb  = new complex <double> [mellinint::mdim*2];
				complex <double> *alogqqp  = new complex <double> [mellinint::mdim*2];
				complex <double> *alogqqbp = new complex <double> [mellinint::mdim*2];
				//complex <double> *aexpqq   = new complex <double> [mellinint::mdim*2];
				//complex <double> *aexpqg   = new complex <double> [mellinint::mdim*2];
				//complex <double> *aexpqqb  = new complex <double> [mellinint::mdim*2];
				//complex <double> *aexpqqp  = new complex <double> [mellinint::mdim*2];
				//complex <double> *aexpqqbp = new complex <double> [mellinint::mdim*2];
				//complex <double> alogqq[mellinint::mdim*2];
				//complex <double> alogqg[mellinint::mdim*2];
				//complex <double> alogqqb[mellinint::mdim*2];
				//complex <double> alogqqp[mellinint::mdim*2];
				//complex <double> alogqqbp[mellinint::mdim*2];
				//complex <double> aexpqq[mellinint::mdim*2];
				//complex <double> aexpqg[mellinint::mdim*2];
				//complex <double> aexpqqb[mellinint::mdim*2];
				//complex <double> aexpqqp[mellinint::mdim*2];
				//complex <double> aexpqqbp[mellinint::mdim*2];

				//No exponentiation of the C coefficients
				if (opts.expc == 0)
					if (opts.sumlogs)
					{
						fill(alogqq, alogqq+mellinint::mdim*2, 0.);
						fill(alogqg, alogqq+mellinint::mdim*2, 0.);
						fill(alogqqb, alogqqb+mellinint::mdim*2, 0.);
						fill(alogqqp, alogqqp+mellinint::mdim*2, 0.);
						fill(alogqqbp, alogqqbp+mellinint::mdim*2, 0.);
					}
					else
					{
						fill(aexpqq, aexpqq+mellinint::mdim*2, 1.);
						fill(aexpqg, aexpqg+mellinint::mdim*2, 1.);
						fill(aexpqqb, aexpqqb+mellinint::mdim*2, 1.);
						fill(aexpqqp, aexpqqp+mellinint::mdim*2, 1.);
						fill(aexpqqbp, aexpqqbp+mellinint::mdim*2, 1.);
					}
					//Only delta terms
					else if (opts.expc == 1)
						if (opts.sumlogs)
						{
							fill(alogqg, alogqq+mellinint::mdim*2, 0.);
							fill(alogqqb, alogqqb+mellinint::mdim*2, 0.);
							fill(alogqqp, alogqqp+mellinint::mdim*2, 0.);
							fill(alogqqbp, alogqqbp+mellinint::mdim*2, 0.);
							for (int sign = mesq::positive; sign <= mesq::negative; sign++)
								for (int i = 0; i < mellinint::mdim; i++)
								{
									int idx = anomalous::index(i,sign);
									alogqq[idx] = lamB*ccoeff::C1qq_delta;
								}
						}
						else
						{
							fill(aexpqg, aexpqg+mellinint::mdim*2, 1.);
							fill(aexpqqb, aexpqqb+mellinint::mdim*2, 1.);
							fill(aexpqqp, aexpqqp+mellinint::mdim*2, 1.);
							fill(aexpqqbp, aexpqqbp+mellinint::mdim*2, 1.);
							for (int sign = mesq::positive; sign <= mesq::negative; sign++)
								for (int i = 0; i < mellinint::mdim; i++)
								{
									int idx = anomalous::index(i,sign);
									aexpqq[idx] = pow(aexpB,ccoeff::C1qq_delta);
								}
						}
						//Taylor expansion
						else if (opts.expc == 2)
							for (int sign = mesq::positive; sign <= mesq::negative; sign++)
								for (int i = 0; i < mellinint::mdim; i++)
								{
									int idx = anomalous::index(i,sign);
									complex <double> qqtayl = 1.;
									complex <double> qgtayl = 1.;
									complex <double> qqbtayl = 1.;
									complex <double> qqptayl = 1.;
									complex <double> qqbptayl = 1.;
									int fact = 1;
									for (int t = 1; t <= opts.ntaylor; t++)
									{
										fact *= t;
										qqtayl += 1./double(fact) *pow(lamB*(ccoeff::C1qq[idx]-ccoeff::C1qq_delta),t);
										qgtayl += 1./double(fact) *pow(lamB*ccoeff::C2qg[idx]/ccoeff::C1qg[idx],t);
										qqbtayl += 1./double(fact) *pow(lamB*ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx],t);
										qqptayl += 1./double(fact) *pow(lamB*ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx],t);
										qqbptayl += 1./double(fact) *pow(lamB*ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx],t);
									}
									if (opts.sumlogs)
									{
										alogqq[idx] = lamB*ccoeff::C1qq_delta + log(qqtayl);
										alogqg[idx] = log(qgtayl);
										alogqqb[idx] = log(qqbtayl);
										alogqqp[idx] = log(qqptayl);
										alogqqbp[idx] = log(qqbptayl);
									}
									else
									{
										aexpqq[idx] = pow(aexpB,ccoeff::C1qq_delta)*qqtayl;
										aexpqg[idx] = qgtayl;
										aexpqqb[idx] = qqbtayl;
										aexpqqp[idx] = qqptayl;
										aexpqqbp[idx] = qqbptayl;
									}
								}
								//Standard formula
								else if (opts.expc == 3)
									for (int sign = mesq::positive; sign <= mesq::negative; sign++)
										for (int i = 0; i < mellinint::mdim; i++)
										{
											int idx = anomalous::index(i,sign);

											if (opts.sumlogs)
											{
												//alogqq = lamB*1./2.*(C1QQ[idx]-2.*C1qqn);
												//alogqg = lamB*1./2.*(C2qgM[idx]/C1QG[idx]-2.*C1qqn);
												alogqq[idx] = lamB*ccoeff::C1qq[idx];
												alogqg[idx] = lamB*ccoeff::C2qg[idx]/ccoeff::C1qg[idx];
												alogqqb[idx] = lamB*ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx];
												alogqqp[idx] = lamB*ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx];
												alogqqbp[idx] = lamB*ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx];
											}
											else
											{
												//aexpqq = pow(aexpB,1./2.*(C1QQ[idx]-2.*C1qqn));
												//aexpqg = pow(aexpB,1./2.*(C2qgM[idx]/C1QG[idx]-2.*C1qqn));
												aexpqq[idx] = pow(aexpB,ccoeff::C1qq[idx]);
												aexpqg[idx] = pow(aexpB,ccoeff::C2qg[idx]/ccoeff::C1qg[idx]);
												aexpqqb[idx] = pow(aexpB,ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx]);
												aexpqqp[idx] = pow(aexpB,ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx]);
												aexpqqbp[idx] = pow(aexpB,ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx]);
											}
										}
										//cout << endl;
										//cout << aexpqq[0]   << " expc::aexpqq " << endl;
										//cout << aexpqg[0]   << " expc::aexpqg " << endl;
										//cout << aexpqqb[0]  << " expc::aexpqqb " << endl;
										//cout << aexpqqp[0]  << " expc::aexpqqb " << endl;
										//cout << aexpqqbp[0] << " expc::aexpqqbp " << endl;
										//cout << endl;
										//At order n expand up to terms containing Cn (i.e. expand up to n for D, up to n-1 for OD, and up to n-2 for DOD)
										else if (opts.expc == 4)
										{
											if (opts.sumlogs)
											{
												fill(alogqg, alogqq+mellinint::mdim*2, 0.);
												fill(alogqqb, alogqqb+mellinint::mdim*2, 0.);
												fill(alogqqp, alogqqp+mellinint::mdim*2, 0.);
												fill(alogqqbp, alogqqbp+mellinint::mdim*2, 0.);
											}
											else
											{
												fill(aexpqg, aexpqg+mellinint::mdim*2, 1.);
												fill(aexpqqb, aexpqqb+mellinint::mdim*2, 1.);
												fill(aexpqqp, aexpqqp+mellinint::mdim*2, 1.);
												fill(aexpqqbp, aexpqqbp+mellinint::mdim*2, 1.);
											}
											for (int sign = mesq::positive; sign <= mesq::negative; sign++)
												for (int i = 0; i < mellinint::mdim; i++)
												{
													int idx = anomalous::index(i,sign);
													if (opts.sumlogs)
														alogqq[idx] = lamB*ccoeff::C1qq[idx];
													else
														aexpqq[idx] = pow(aexpB,ccoeff::C1qq[idx]);
												}
										}
										//Do not expand C at the denominator
										else if (opts.expc == 5)
											for (int sign = mesq::positive; sign <= mesq::negative; sign++)
												for (int i = 0; i < mellinint::mdim; i++)
												{
													int idx = anomalous::index(i,sign);
													if (opts.sumlogs)
													{
														alogqq[idx] = lamB*ccoeff::C1qq[idx];
														alogqg[idx]   = log(1.-xlambda/(1.+aass*ccoeff::C2qg[idx]/ccoeff::C1qg[idx])) - log(1.-xlambda);
														alogqqb[idx]  = log(1.-xlambda/(1.+aass*ccoeff::C3qqb[idx] /ccoeff::C2qqb[idx] )) - log(1.-xlambda);
														alogqqp[idx]  = log(1.-xlambda/(1.+aass*ccoeff::C3qqp[idx] /ccoeff::C2qqp[idx] )) - log(1.-xlambda);
														alogqqbp[idx] = log(1.-xlambda/(1.+aass*ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx])) - log(1.-xlambda);
														//alternative with as0 at the denominator
														//alogqg[idx]   = lamB*ccoeff::C2qg[idx]  /ccoeff::C1qg[idx]  /(1.+aass*ccoeff::C2qg[idx]  /ccoeff::C1qg[idx]  );
														//alogqqb[idx]  = lamB*ccoeff::C3qqb[idx] /ccoeff::C2qqb[idx] /(1.+aass*ccoeff::C3qqb[idx] /ccoeff::C2qqb[idx] );
														//alogqqp[idx]  = lamB*ccoeff::C3qqp[idx] /ccoeff::C2qqp[idx] /(1.+aass*ccoeff::C3qqp[idx] /ccoeff::C2qqp[idx] );
														//alogqqbp[idx] = lamB*ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx]/(1.+aass*ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx]);

													}
													else
													{
														aexpqq[idx] = pow(aexpB,ccoeff::C1qq[idx]);
														aexpqg[idx]   = (1.-xlambda/(1.+aass*ccoeff::C2qg[idx]/ccoeff::C1qg[idx]))/(1.-xlambda);
														aexpqqb[idx]  = (1.-xlambda/(1.+aass*ccoeff::C3qqb[idx] /ccoeff::C2qqb[idx] ))/(1.-xlambda);
														aexpqqp[idx]  = (1.-xlambda/(1.+aass*ccoeff::C3qqp[idx] /ccoeff::C2qqp[idx] ))/(1.-xlambda);
														aexpqqbp[idx] = (1.-xlambda/(1.+aass*ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx]))/(1.-xlambda);

													}
												}
												//gint::calc(b);
												//cout << endl;
												//cout << " gint::logasl "   << gint::logasl      << " expc " << logasl      << endl;
												//cout << " gint::alogqq "   << gint::alogqq[0]   << " expc " << alogqq[0]   << endl;
												//cout << " gint::alogqg "   << gint::alogqg[0]   << " expc " << alogqg[0]   << endl;
												//cout << " gint::alogqqb "  << gint::alogqqb[0]  << " expc " << alogqqb[0]  << endl;
												//cout << " gint::alogqqp "  << gint::alogqqp[0]  << " expc " << alogqqp[0]  << endl;
												//cout << " gint::alogqqbp " << gint::alogqqbp[0] << " expc " << alogqqbp[0] << endl;

												//cout << aexpqq[0]   << " gint::alogqq " << endl;
												//cout << aexpqg[0]   << " gint::alogqg " << endl;
												//cout << aexpqqb[0]  << " gint::alogqqb " << endl;
												//cout << aexpqqp[0]  << " gint::alogqqb " << endl;
												//cout << aexpqqbp[0] << " gint::alogqqbp " << endl;
												//cout << endl;

												//!!! gint::alog are evaluated with blim::sudakov, not with blim:expc !!!
												if (opts.numexpc)
													for (int sign = mesq::positive; sign <= mesq::negative; sign++)
														for (int i = 0; i < mellinint::mdim; i++)
														{
															int idx = anomalous::index(i,sign);

															if (opts.sumlogs)
															{
																alogqq[idx]   = gint::alogqq[idx];
																alogqg[idx]   = gint::alogqg[idx];
																alogqqb[idx]  = gint::alogqqb[idx];
																alogqqp[idx]  = gint::alogqqp[idx];
																alogqqbp[idx] = gint::alogqqbp[idx];
															}
															else
															{
																aexpqq[idx]   = exp(gint::alogqq[idx]);
																aexpqg[idx]   = exp(gint::alogqg[idx]);
																aexpqqb[idx]  = exp(gint::alogqqb[idx]);
																aexpqqp[idx]  = exp(gint::alogqqp[idx]);
																aexpqqbp[idx] = exp(gint::alogqqbp[idx]);
															}
														}

														for (int sign = mesq::positive; sign <= mesq::negative; sign++)
															for (int i = 0; i < mellinint::mdim; i++)
															{
																int idx = anomalous::index(i,sign);
																int ii = hcoefficients::index(i,sign);
																if (opts.sumlogs)
																{
																	qqb[ii] = exp(alogqq[idx]+alogqq[idx]+sudakov::logS);
																	qg[ii]  = exp(-logasl+alogqg[idx]+alogqq[idx]+sudakov::logS);
																	qq[ii]  = exp(alogqq[idx]+alogqqb[idx]-2.*logasl+sudakov::logS);
																	qqp[ii] = exp(alogqq[idx]+alogqqbp[idx]-2.*logasl+sudakov::logS);
																	qqbp[ii] = exp(alogqq[idx]+alogqqp[idx]-2.*logasl+sudakov::logS);
																	gg[ii]  = exp(-logasl+alogqg[idx]-logasl+alogqg[idx]+sudakov::logS);
																}
																else
																{
																	qqb[ii] = aexpqq[idx] * aexpqq[idx];
																	qg[ii]  = aexp * aexpqg[idx]*aexpqq[idx];
																	qq[ii]  = aexpqq[idx] * aexpqqb[idx] * aexp2; //*c1delta
																	qqp[ii] = aexpqq[idx] * aexpqqbp[idx] * aexp2; //*c1delta
																	qqbp[ii] = aexpqq[idx] * aexpqqp[idx] * aexp2; //*c1delta
																	gg[ii]  = aexp*aexpqg[idx] * aexp*aexpqg[idx];
																}
															}

															delete[] alogqq;
															delete[] alogqg;
															delete[] alogqqb;
															delete[] alogqqp;
															delete[] alogqqbp;

															//include aexp into aexpqg and aexpqq,qqp,qqb,qqbp
															for (int sign = mesq::positive; sign <= mesq::negative; sign++)
																for (int i = 0; i < mellinint::mdim; i++)
																{
																	int idx = anomalous::index(i,sign);
																	aexpqg[idx]   *= aexp;
																	aexpqqb[idx]  *= aexp2;
																	aexpqqp[idx]  *= aexp2;
																	aexpqqbp[idx] *= aexp2;
																}
																//delete[] aexpqq;
																//delete[] aexpqg;
																//delete[] aexpqqb;
																//delete[] aexpqqp;
																//delete[] aexpqqbp;

			}

			// Mellin2d - rapidity dependent

			else {

				//complex <double> aexpqq_1[mellinint::mdim*2];
				//complex <double> aexpqg_1[mellinint::mdim*2];
				//complex <double> aexpqqb_1[mellinint::mdim*2];
				//complex <double> aexpqqp_1[mellinint::mdim*2];
				//complex <double> aexpqqbp_1[mellinint::mdim*2];
				//complex <double> aexpqq_2[mellinint::mdim*2];
				//complex <double> aexpqg_2[mellinint::mdim*2];
				//complex <double> aexpqqb_2[mellinint::mdim*2];
				//complex <double> aexpqqp_2[mellinint::mdim*2];
				//complex <double> aexpqqbp_2[mellinint::mdim*2];

				//case opts.expc == 0
				fill(aexpqq_1, aexpqq_1+mellinint::mdim*2, 1.);
				fill(aexpqg_1, aexpqg_1+mellinint::mdim*2, 1.);
				fill(aexpqqb_1, aexpqqb_1+mellinint::mdim*2, 1.);
				fill(aexpqqp_1, aexpqqp_1+mellinint::mdim*2, 1.);
				fill(aexpqqbp_1, aexpqqbp_1+mellinint::mdim*2, 1.);
				fill(aexpqq_2, aexpqq_2+mellinint::mdim*2, 1.);
				fill(aexpqg_2, aexpqg_2+mellinint::mdim*2, 1.);
				fill(aexpqqb_2, aexpqqb_2+mellinint::mdim*2, 1.);
				fill(aexpqqp_2, aexpqqp_2+mellinint::mdim*2, 1.);
				fill(aexpqqbp_2, aexpqqbp_2+mellinint::mdim*2, 1.);

				//Only delta terms
				if (opts.expc == 1) {

					for (int sign = mesq::positive; sign <= mesq::negative; sign++)
						for (int i = 0; i < mellinint::mdim; i++)
						{
							int idx = anomalous::index(i,sign);
							aexpqq_1[idx] = pow(aexpB,ccoeff::C1qq_delta);
							aexpqq_2[idx] = pow(aexpB,ccoeff::C1qq_delta);
						}
				}

				//Taylor expansion
				else if (opts.expc == 2) {

					{
						for (int sign = mesq::positive; sign <= mesq::negative; sign++)
							for (int i = 0; i < mellinint::mdim; i++)
							{
								int idx = anomalous::index(i,sign);
								complex <double> qqtayl_1 = 1.;
								complex <double> qgtayl_1 = 1.;
								complex <double> qqbtayl_1 = 1.;
								complex <double> qqptayl_1 = 1.;
								complex <double> qqbptayl_1 = 1.;
								complex <double> qqtayl_2 = 1.;
								complex <double> qgtayl_2 = 1.;
								complex <double> qqbtayl_2 = 1.;
								complex <double> qqptayl_2 = 1.;
								complex <double> qqbptayl_2 = 1.;
								int fact = 1;
								for (int t = 1; t <= opts.ntaylor; t++)
								{
									fact *= t;
									qqtayl_1 += 1./double(fact) *pow(lamB*(ccoeff::C1qq_1[idx]-ccoeff::C1qq_delta),t);
									qgtayl_1 += 1./double(fact) *pow(lamB*ccoeff::C2qg_1[idx]/ccoeff::C1qg_1[idx],t);
									qqbtayl_1 += 1./double(fact) *pow(lamB*ccoeff::C3qqb_1[idx]/ccoeff::C2qqb_1[idx],t);
									qqptayl_1 += 1./double(fact) *pow(lamB*ccoeff::C3qqp_1[idx]/ccoeff::C2qqp_1[idx],t);
									qqbptayl_1 += 1./double(fact) *pow(lamB*ccoeff::C3qqbp_1[idx]/ccoeff::C2qqbp_1[idx],t);
									qqtayl_2 += 1./double(fact) *pow(lamB*(ccoeff::C1qq_2[idx]-ccoeff::C1qq_delta),t);
									qgtayl_2 += 1./double(fact) *pow(lamB*ccoeff::C2qg_2[idx]/ccoeff::C1qg_2[idx],t);
									qqbtayl_2 += 1./double(fact) *pow(lamB*ccoeff::C3qqb_2[idx]/ccoeff::C2qqb_2[idx],t);
									qqptayl_2 += 1./double(fact) *pow(lamB*ccoeff::C3qqp_2[idx]/ccoeff::C2qqp_2[idx],t);
									qqbptayl_2 += 1./double(fact) *pow(lamB*ccoeff::C3qqbp_2[idx]/ccoeff::C2qqbp_2[idx],t);
								}
								aexpqq_1[idx] = pow(aexpB,ccoeff::C1qq_delta)*qqtayl_1;
								aexpqg_1[idx] = qgtayl_1;
								aexpqqb_1[idx] = qqbtayl_1;
								aexpqqp_1[idx] = qqptayl_1;
								aexpqqbp_1[idx] = qqbptayl_1;
								aexpqq_2[idx] = pow(aexpB,ccoeff::C1qq_delta)*qqtayl_2;
								aexpqg_2[idx] = qgtayl_2;
								aexpqqb_2[idx] = qqbtayl_2;
								aexpqqp_2[idx] = qqptayl_2;
								aexpqqbp_2[idx] = qqbptayl_2;
							}
					}

				}

				//Standard expression
				else if (opts.expc == 3) {

					for (int sign = mesq::positive; sign <= mesq::negative; sign++) {

						for (int i = 0; i < mellinint::mdim; i++) {

							int idx = anomalous::index(i,sign);
							// aexpgq_1[idx] = pow(aexpB,(ccoeff::C2gq_1[idx]/ccoeff::C1gq_1[idx] - ccoeff::C1gg)/2.);
							aexpgq_1[idx] = pow(aexpB,(ccoeff::C2gq_1[idx]/ccoeff::C1gq_1[idx] - ccoeff::C1gg)/2.);
							// aexpgq_2[idx] = pow(aexpB,(ccoeff::C1gq_2[idx]/ccoeff::C1gq_2[idx] - ccoeff::C1gg)/2.);
							aexpgq_2[idx] = pow(aexpB,(ccoeff::C2gq_2[idx]/ccoeff::C1gq_2[idx] - ccoeff::C1gg)/2.);
							// aexpqq_1[idx] = pow(aexpB,(ccoeff::C2gq_1[idx]/ccoeff::C1gq_1[idx] - ccoeff::C1gg)/2.);
							aexpqq_1[idx] = pow(aexpB,(ccoeff::C2gq_2[idx]/ccoeff::C1gq_2[idx] - ccoeff::C1gg)/2.);
							// aexpqq_2[idx] = pow(aexpB,(ccoeff::C2gq_1[idx]/ccoeff::C1gq_1[idx] - ccoeff::C1gg)/2.);
							aexpqq_2[idx] = pow(aexpB,(ccoeff::C2gq_2[idx]/ccoeff::C1gq_2[idx] - ccoeff::C1gg)/2.);

						}

					}

				}

				//At order n expand up to terms containing Cn (i.e. expand up to n for D, up to n-1 for OD, and up to n-2 for DOD)
				else if (opts.expc == 4) {

					{
						for (int sign = mesq::positive; sign <= mesq::negative; sign++)
							for (int i = 0; i < mellinint::mdim; i++)
							{
								int idx = anomalous::index(i,sign);
								aexpqq_1[idx] = pow(aexpB,ccoeff::C1qq_1[idx]);
								aexpqq_2[idx] = pow(aexpB,ccoeff::C1qq_2[idx]);
							}
					}

				}

				//Do not expand C at the denominator in the off diagonal and double off diagonal channels
				else if (opts.expc == 5) {

					for (int sign = mesq::positive; sign <= mesq::negative; sign++)
						for (int i = 0; i < mellinint::mdim; i++)
						{
							int idx = anomalous::index(i,sign);
							aexpqq_1[idx] = pow(aexpB,ccoeff::C1qq_1[idx]);
							aexpqg_1[idx]   = (1. - xlambda/(1.+aass*ccoeff::C2qg_1[idx]/ccoeff::C1qg_1[idx]))/(1.-xlambda);
							aexpqqb_1[idx]  = (1. - xlambda/(1.+aass*ccoeff::C3qqb_1[idx] /ccoeff::C2qqb_1[idx] ))/(1.-xlambda);
							aexpqqp_1[idx]  = (1. - xlambda/(1.+aass*ccoeff::C3qqp_1[idx] /ccoeff::C2qqp_1[idx] ))/(1.-xlambda);
							aexpqqbp_1[idx] = (1. - xlambda/(1.+aass*ccoeff::C3qqbp_1[idx]/ccoeff::C2qqbp_1[idx]))/(1.-xlambda);
							aexpqq_2[idx] = pow(aexpB,ccoeff::C1qq_2[idx]);
							aexpqg_2[idx]   = (1. - xlambda/(1.+aass*ccoeff::C2qg_2[idx]/ccoeff::C1qg_2[idx]))/(1.-xlambda);
							aexpqqb_2[idx]  = (1. - xlambda/(1.+aass*ccoeff::C3qqb_2[idx] /ccoeff::C2qqb_2[idx] ))/(1.-xlambda);
							aexpqqp_2[idx]  = (1. - xlambda/(1.+aass*ccoeff::C3qqp_2[idx] /ccoeff::C2qqp_2[idx] ))/(1.-xlambda);
							aexpqqbp_2[idx] = (1. - xlambda/(1.+aass*ccoeff::C3qqbp_2[idx]/ccoeff::C2qqbp_2[idx]))/(1.-xlambda);
						}

				}

				// CJU
				// Compute the aexp, aexpgq and aexpqq coefficients
				for (int sign = mesq::positive; sign <= mesq::negative; sign++) {

					for (int i1 = 0; i1 < mellinint::mdim; i1++) {

						for (int i2 = 0; i2 < mellinint::mdim; i2++) {

							int idx1 = anomalous::index(i1,mesq::positive);
							int idx2 = anomalous::index(i2,sign);
							int ii12 = hcoefficients::index(i1,i2,sign);

							gg[ii12]   += pow(aass, 2)*pow(constants::CA, 2)*ccoeff::G1N_1[idx1]*ccoeff::G1N_2[idx2] * pow(aexp, 2);
							// gg[ii12]   = 0.;

							// gq_1[ii12] = aexp*aexpgq_1[idx1];
							gq_1[ii12] = aexp*aexpgq_1[idx1] + pow(aass, 2)*constants::CA*constants::CF*ccoeff::G1N_1[idx1]*ccoeff::G1N_2[idx2] * pow(aexp, 2);
							// gq_1[ii12]   = 0.;

							// gq_2[ii12] = aexp*aexpgq_2[idx2];
							gq_2[ii12] = aexp*aexpgq_2[idx2] + pow(aass, 2)*constants::CA*constants::CF*ccoeff::G1N_1[idx1]*ccoeff::G1N_2[idx2] * pow(aexp, 2);
							// gq_2[ii12]   = 0.;

							// qq_1[ii12] = aexp*aexpqq_1[idx1];
							qq_1[ii12] = aexp*aexpqq_1[idx1] + pow(aass, 2)*pow(constants::CF, 2)*ccoeff::G1N_1[idx1]*ccoeff::G1N_2[idx2] * pow(aexp, 2);
							// qq_1[ii12] += pow(aass, 2) * pow(constants::CF, 2) * ccoeff::G1N_1[i1] * ccoeff::G1N_2[i2] * pow(aexp, 2);
							// qq_1[ii12]   = 0.;

							// qq_2[ii12] = aexp*aexpqq_2[idx2];
							qq_2[ii12] = aexp*aexpqq_2[idx2] + pow(aass, 2)*pow(constants::CF, 2)*ccoeff::G1N_1[idx1]*ccoeff::G1N_2[idx2] * pow(aexp, 2);
							// qq_2[ii12] += pow(aass, 2) * pow(constants::CF, 2) * ccoeff::G1N_1[i1] * ccoeff::G1N_2[i2] * pow(aexp, 2);
							// qq_2[ii12]   = 0.;

							// Check numerical values
							// cout.precision(10);
							// cout << "aass = " << aass << endl;
							// cout << "aexp = " << aexp << endl;
							// cout << "\naexpB = " << aexpB << endl;
							// cout << "aexpgq_1 = " << aexpgq_1[idx1] << endl;
							// cout << "aexpgq_2 = " << aexpgq_2[idx2] << endl;
							// cout << "aexpqq_1 = " << aexpqq_1[idx1] << endl;
							// cout << "aexpqq_2 = " << aexpqq_2[idx2] << endl;
							// cout << "gq_1 = " << gq_1[ii12] << endl;
							// cout << "gq_2 = " << gq_2[ii12] << endl;
							// cout << "qq_1 = " << qq_1[ii12] << endl;
							// cout << "qq_2 = " << qq_2[ii12] << endl;

						}

					}

				}

				// Include aexp into aexpgq and aexpqq
				for (int sign = mesq::positive; sign <= mesq::negative; sign++) {

					for (int i = 0; i < mellinint::mdim; i++) {

						int idx = anomalous::index(i,sign);

						aexpgq_1[idx]   *= aexp;
						aexpgq_2[idx]   *= aexp;
						aexpqq_1[idx]   *= aexp;
						aexpqq_2[idx]   *= aexp;
						// aexpqg_1[idx]   *= aexp;
						// aexpqqp_1[idx]  *= aexp2;
						// aexpqqbp_1[idx] *= aexp2;
						// aexpqg_2[idx]   *= aexp;
						// aexpqqb_2[idx]  *= aexp2;

					}

				}

			}

		}

		// NNNLL

		//multiple parton case as in Eqs.(100-107) of https://arxiv.org/pdf/hep-ph/0508068.pdf
		//B3tilde additional pieces as in Eq. (50)
		//-2 beta1 C1 + 2 beta0 C1^2 - 4 beta0 C2

		//(B3qbar*(-2. + y)*y)/(2.*beta0*pow(1. - y,2))
		//1./beta0 * 1./2.  y*(y-2.)/pow(1.-y,2) * (-2 beta1 C1 + 2 beta0 C1^2 - 4 beta0 C2)
		//log(aexpC) * (C1^2 - beta1/beta0*C1 - 2*C2)

		//pieces from g4
		//(B3qbar*(-2. + y)*y)/(2.*beta0*pow(1. - y,2))
		//(B2qbar*beta1*((2. - y)*y + 2.*log1y))/(2.*pow(beta0,2)*pow(1. - y,2))
		//((rlogq2mur2-2.*rloga)*(3.*y*(-2.*B2qbar*pow(beta0,3)*(-2. + y)*(-1. + y)))/(6.*pow(beta0,3)*pow(-1. + y,3));

		//(1/2)*(1/beta0)*B3qbar       *y*(y-2)/(1-y)^2
		//-(1/2)*B2qbar*beta1/beta0^2  *y*(y-2)/(1-y)^2
		//B2qbar*beta1/beta0^2        *log1y  /(1-y)^2
		//-B2qbar((rlogq2mur2-2.*rloga)*y*(y-2)/(1-y)^2

		//(-beta1/beta0*C1+C1^2-2*C2)       *y*(y-2)/(1-y)^2
		//C1*beta1/beta0                      *y*(y-2)/(1-y)^2
		//-2*C1*beta1/beta0                    *log1y  /(1-y)^2
		//(2*beta0*C1)*((rlogq2mur2-2.*rloga) *y*(y-2)/(1-y)^2

		//Diagonal form factor
		//Gca_D =
		//  + (aS/Pi)*C1                           *(lam)/(1-lam)
		//
		//  + (aS^2/Pi^2)*(
		//		   + (C1^2/2 - C2)           *lam*(lam-2)/(1-lam)^2
		//		   - beta1/beta0*C1          *Log[1-lam]/(1-lam)^2
		//		   + beta0*C1*Log[Q2/muR2]   *lam*(lam-2)/(1-lam)^2
		//		   );

		//Off-diagonal form factor
		//Gca_OD =
		//  - Log[1 - lam]
		//
		//  + (aS/Pi)*(
		//	       + (C2/C1)                         *lam/(1-lam)
		//	       - (beta1/beta0                    *Log[1-lam])/(1-lam)
		//	       - (beta0*Log[Q2/muR2])            *lam/(1-lam)
		//	       )
		//
		//  + (aS^2/Pi^2)*(
		//		   - beta1^2/beta0^2                        *lam/(1-lam)^2
		//		   + beta2/beta0                            *lam/(1-lam)^2
		//		   - beta1^2/beta0^2                        *Log[1-lam]/(1-lam)^2
		//		   + 1/2*beta1^2/beta0^2                    *Log[1-lam]^2/(1-lam)^2
		//		   + (C2^2/C1^2/2 - C1*C3/C1^2)             *lam*(lam-2)/(1-lam)^2
		//		   - beta1/beta0*C2/C1                      *Log[1-lam]/(1-lam)^2
		//
		//		   - beta1                  *Log[Q2/muR2]   *lam/(1-lam)^2
		//		   + beta1                  *Log[Q2/muR2]   *Log[1-lam]/(1-lam)^2
		//		   + beta0*C2/C1            *Log[Q2/muR2]   *lam*(lam-2)/(1-lam)^2
		//		   );

		//Doubly off-diagonal form factor
		//Gca_DOD =
		//  - 2*Log[1 - lam]
		//  + (aS/Pi)* (
		//		+ (C3/C2)                      *lam/(1-lam)
		//		- (2*beta1/beta0               *Log[1-lam])/(1-lam)
		//		- (2*beta0*Log[Q2/muR2])       *lam/(1-lam)
		//		)
		//
		//  + (aS^2/Pi^2)* (
		//		    - 2*beta1^2/beta0^2                    *lam/(1-lam)^2
		//		    + 2*beta2/beta0                        *lam/(1-lam)^2
		//		    - 2*beta1^2/beta0^2                    *Log[1-lam]/(1-lam)^2
		//		    + beta1^2/beta0^2                      *Log[1-lam]^2/(1-lam)^2
		//		    + (C3^2/C2^2/2 - C4/C2)                *lam*(lam-2)/(1-lam)^2
		//		    + beta1/beta0*C3/C2                    *Log[1-lam]/(1-lam)^2
		//
		//		    - 2*beta1      *Log[Q2/muR2]          *lam/(1-lam)^2
		//		    + 2*beta1      *Log[Q2/muR2]          *Log[1-lam]/(1-lam)^2
		//		    + beta0*C3/C2  *Log[Q2/muR2]          *lam*(lam-2)/(1-lam)^2
		//		    );

		// CJU - N3LL
		if (opts.order_expc == 3) {

			if (opts.mellin1d) {

				cout << "Order 3 not implemented for Higgs" << endl;

			}

			// Mellin2d
			else {

				// cout << "Order 3 not implemented for Higgs" << endl;

				for (int sign = mesq::positive; sign <= mesq::negative; sign++) {

					for (int i = 0; i < mellinint::mdim; i++) {

						int idx = anomalous::index(i,sign);
						aexpgq_1[idx] += 0.;
						aexpgq_2[idx] += 0.;
						aexpqq_1[idx] += 0.;
						aexpqq_2[idx] += 0.;

					}

				}

			}

		}

		if (opts.order_expc >= 4) {

			cout << "Order not implemented" << endl;

		}


		//  return;
		if (opts.npff == 3)
			if (opts.mellin1d)
				for (int sign = mesq::positive; sign <= mesq::negative; sign++)
					for (int i = 0; i < mellinint::mdim; i++)
					{
						complex <double> N;
						if (sign == mesq::positive)
							N = mellinint::Np[i];
						else
							N = mellinint::Nm[i];
						complex <double> S1N = 0.;//cpsi0(N)+constants::euler;
						complex <double> dmw = exp(-pow(b,2)/2.*(opts.a2*(-2.*S1N-(4.*N+4.)/(N*(N+2.)))));
						int ii = hcoefficients::index(i,sign);
						qqb[ii]  *= dmw*dmw;
						qg[ii]   *= dmw*dmw;
						qq[ii]   *= dmw*dmw;
						qqp[ii]  *= dmw*dmw;
						qqbp[ii] *= dmw*dmw;
						gg[ii]   *= dmw*dmw;
						qbg[ii]  *= dmw*dmw;
						qpg[ii]  *= dmw*dmw;
						qbpg[ii] *= dmw*dmw;
					}
					else
						for (int sign = mesq::positive; sign <= mesq::negative; sign++)
							for (int i1 = 0; i1 < mellinint::mdim; i1++)
								for (int i2 = 0; i2 < mellinint::mdim; i2++)
								{
									complex <double> N_1,N_2;
									if (sign == mesq::positive)
									{
										N_1 = mellinint::Np_1[i1];
										N_2 = mellinint::Np_2[i2];
									}
									else
									{
										N_1 = mellinint::Nm_1[i1];
										N_2 = mellinint::Nm_2[i2];
									}
									complex <double> S1N_1 = cpsi0(N_1+1.)+constants::euler;
									complex <double> S1N_2 = cpsi0(N_2+1.)+constants::euler;
									complex <double> dmw_1 = exp(-pow(b,2)/2.*(opts.a2*(-2.*S1N_1-(4.*N_1+4.)/(N_1*(N_1+2.)))));
									complex <double> dmw_2 = exp(-pow(b,2)/2.*(opts.a2*(-2.*S1N_2-(4.*N_2+4.)/(N_2*(N_2+2.)))));
									int ii12 = hcoefficients::index(i1,i2,sign);
									qqb[ii12]    *= dmw_1*dmw_2;
									qg_1[ii12]   *= dmw_1;
									qg_2[ii12]   *= dmw_2;
									qq_1[ii12]   *= 1.;
									qq_2[ii12]   *= 1.;
									qqp_1[ii12]  *= 1.;
									qqp_2[ii12]  *= 1.;
									qqbp_1[ii12] *= 1.;
									qqbp_2[ii12] *= 1.;
									gg[ii12]     *= 1.;
									qbg_1[ii12]  *= 1.;
									qbg_2[ii12]  *= 1.;
									qpg_1[ii12]  *= 1.;
									qpg_2[ii12]  *= 1.;
									qbpg_1[ii12] *= 1.;
									qbpg_2[ii12] *= 1.;
								}

								//check nans
								if (opts.mellin1d) {

									for (int sign = mesq::positive; sign <= mesq::negative; sign++) {

										for (int i = 0; i < mellinint::mdim; i++) {

											int ii = hcoefficients::index(i,sign);
											if (isnan_ofast(real(qg[ii])) || isnan_ofast(imag(qg[ii]))) {

												int idx = anomalous::index(i,sign);
												cout << "Warning, expc::qg[" << i << "] = " << qg[ii] << ", b = "  << b << ", lambda = " << xlambda << " logasl " << logasl
												//<< " aexpqq " << aexpqq[idx] << " aexpqg " << aexpqg[idx]
												<< endl;
												if (opts.sumlogs)
													qg[ii] = sudakov::S;
												else
													qg[ii] = 1.;

											}

										}

									}

								}

	}
	// Drell-Yan
	else{
	//Return 1 if b is close to the Landau pole
	double b_L;
	if (!opts.modlog)
		b_L = b0/scales::res * exp(1./(2.*aass*beta0));
		else if (opts.p == 1)
			b_L = b0/scales::res * sqrt(exp(1./(aass*beta0))-1.);
		else
			b_L = b0/scales::res * pow(sqrt(exp(opts.p/(aass*beta0))-1.),1./opts.p);

	if (abs(b-b_L) < 1e-8)
		return;

	//  if (opts.expc == 0)
	//    return;

	//LL (no evolution)
	if (opts.order_expc == 0)
		return;

	double blim = blim::expc;  //better use blim::sudakov?

	//  the limit below implies xlambda<1/2 and then aa2<= 1
	//      blim=b0p*(1/q)*exp(1/(2*as*beta0))
	//      blim=b0p*(1/q)*exp(1/(4*as*beta0))
	// Set a limit to avoid very large values of b (= very small scales ~1/b)
	//       blim=b0p*(1/q)*exp(1/(2*aass*beta0)) ! avoid Landau pole
	//     write(*,*) "blim",blim
	//     without this additional blim some scale variations will fail (when mures > muren)
	//  blim=0.5; // --> allow this to a separate blim in the settings, or set it using muren instead of mures
	//      blim=1.1229190d0
	//      blim=3d0

	//Set b according to bstar or other prescriptions
	complex <double> bstar;
	if (opts.bprescription == 0 || opts.bprescription == 4 || opts.bstar_expc)
		//bstar = real(b)/sqrt(1.+pow(real(b)/blim,2));
		bstar = bstar::calc(b,blim);
	else
		bstar = b;

	double Q = scales::res;
	double aass2 = pow(aass,2);
	double aass3 = pow(aass,3);


	complex <double> blog;
	if (!opts.modlog)
		blog = log(pow(Q*bstar/b0,2));   //normal sudakov
	else if (opts.p == 1)
		blog = log(pow(Q*bstar/b0,2) + 1.); //modified sudakov
	else
		blog = 1./opts.p*log(pow(Q*bstar/b0,2*opts.p) + 1.); //modified sudakov with exponent p

	//if (opts.modlog)
	//  blog = log(pow(Q*bstar/b0,2) + 1.); //modified sudakov
	//else
	//  blog = log(pow(Q*bstar/b0,2));   //normal sudakov


	//lambda as defined in Eq. (25) of hep-ph/0508068.
	xlambda = beta0*aass*blog;
	complex <double> log1xlambda = log(1.-xlambda);
	xlambdaqed = beta0p*aaqed*blog;
	complex <double> log1xlambdaqed = log(1.-xlambdaqed);

	complex <double> logasl;
	complex <double> logaslqed;

	//LL (no evolution)
	logasl = 0.;
	logaslqed = 0.;

	//NLL (LO evolution)
	if (opts.order_expc >= 1)
		{
		logasl += log1xlambda;
		if(opts.qed) logaslqed += log1xlambdaqed;
		}
	//NNLL (NLO evolution)
	if (opts.order_expc >= 2)
		logasl += aass* beta1/beta0*log1xlambda/(1.-xlambda);

	//NNNLL (NNLO evolution)
	if (opts.order_expc >= 3)
		logasl += aass2* ((pow(beta1/beta0,2)-beta2/beta0) *xlambda/pow(1.-xlambda,2)
				+ pow(beta1/beta0,2)             *log1xlambda/pow(1.-xlambda,2)
				- pow(beta1/beta0,2)             *pow(log1xlambda,2)/(2.*pow(1.-xlambda,2)));

	//NNNNLL (NNNLO evolution)
	if (opts.order_expc >= 4)
		logasl += aass3*
		(3.*xlambda*(pow(beta0,2)*beta3*(-2. + xlambda)
			- 2.*beta0*beta1*beta2*(-1. + xlambda)
			+ pow(beta1,3)*xlambda)
		+ 6.*(-(pow(beta1,3)*xlambda)
			+ beta0*beta1*beta2*(1. + xlambda))*log1xlambda
		- 9.*pow(beta1,3)*pow(log1xlambda,2)
		+ 2.*pow(beta1,3)*pow(log1xlambda,3)
		)/(6.*pow(beta0,3)*pow(1. - xlambda,3));

	//double LQR = real(logq2mur2-2.*loga);
	double LQR = LR-LQ;
	//QCD coupling scale dependence
	if (opts.order_expc >= 2)
		logasl += aass*LQR*beta0*xlambda/(1.-xlambda);

	//from c3new.txt
	if (opts.order_expc >= 3)
		logasl += aass2*(+LQR*beta1                   *(xlambda-log1xlambda)/pow(1.-xlambda,2)
				+LQR*beta1                   *xlambda/(1.-xlambda)                      //missing piece
					+0.5*pow(LQR,2)*pow(beta0,2) *xlambda*(xlambda-2.)/pow(1.-xlambda,2));  //missing piece

	//from g4 LQR part by replacing B1qbar->beta0 and B2qbar->beta1 (and all An pieces set to zero)
	//if (opts.order_expc >= 3)
	//  logasl += aass2*(LQR*beta1*(xlambda*(2.-xlambda)-log1xlambda)/pow(1.-xlambda,2)
	//		     +0.5*pow(LQR,2)*pow(beta0,2)*xlambda*(xlambda-2.)/pow(1.-xlambda,2));

	if (opts.order_expc >= 4)
		logasl += aass3*
		(LQR*(6.*pow(beta1,2)*pow(log1xlambda,2)
			+ 6.*beta1*log1xlambda*(-3.*beta1 + pow(beta0,2)*LQR)
			+ xlambda*(-6.*pow(beta1,2) + 6.*beta0*beta2*(4. - 3.*xlambda + pow(xlambda,2))
				- 9.*pow(beta0,2)*beta1*(3. - 3.*xlambda + pow(xlambda,2))*LQR
				+ 2.*pow(beta0,4)*(3. - 3.*xlambda + pow(xlambda,2))*pow(LQR,2))))/(6.*beta0*pow(1. - xlambda,3));

	//cout << logasl << " logasl" << endl;
	//gint::logasl is evaluated with blim::sudakov, while gint::logasl_expc is evaluated with blim:expc
	if (opts.numexpc)
		logasl = gint::logasl_expc;
	//cout << gint::logasl_expc << " gint logasl" << endl;

	if (isnan_ofast(real(logasl)) || isnan_ofast(imag(logasl)))
		{
		cout << "Warning, logasl = " << logasl << ", b = "  << b << ", lambda = " << xlambda
		<< " log(1-xlambda) " << log1xlambda
		<< endl;
		logasl = 0.;
		}

	//For the aexp evolution operator use the same hard scale as for the PDF evolution, and the same blim
	//aexp should always reflect the PDF evolution, in evolmode 2 and 3 it should be log(ASF/ASI), plus eventually the LQR scale variation pieces
	//cout << "aexp " << logasl << " pdfs " << pdfevol::logasl << endl;
	//--> No aexp should be with Q and not muf

	//if (opts.mufevol)
	//complex <double> logasl = pdfevol::logasl;

	complex <double> aexp = exp(-logasl); //aexp is approximately alphas(Q^2)/alphas(mub^2)
	complex <double> aexpqed = exp(-logaslqed);
	//aexp = exp(-logasl); //aexp is approximately alphas(Q^2)/alphas(mub^2)

	if (isnan_ofast(real(aexp)) || isnan_ofast(imag(aexp)))
		{
		cout << "Warning, aexp = " << aexp << ", b = "  << b << ", lambda = " << xlambda
		<< " log(1-xlambda) " << log1xlambda
		<< endl;
		aexp = 1.;
		}


	/*
	//Set b according to bstar or other prescriptions
	//For the C exponentiation always use Q at the resummation scale, and blim::sudakov, to match what is done in the Sudakov
	blim = blim::expc; //better use blim::sudakov;
	complex <double> bstar;
	if (opts.bprescription == 0 || opts.bstar_expc) //better use opts.bstar_sudakov)
		bstar = real(b)/sqrt(1.+pow(real(b)/blim,2));
	else
		bstar = b;

	double Q = scales::res;
	complex <double> blog;
	if (opts.modlog)
		blog = log(pow(Q*bstar/b0,2) + 1.); //modified sudakov
	else
		blog = log(pow(Q*bstar/b0,2));   //normal sudakov

	//lambda as defined in Eq. (25) of hep-ph/0508068.
	complex <double> xlambda = beta0*aass*blog;
	complex <double> log1xlambda = log(1.-xlambda);

	Q = scales::res;
	if (opts.modlog)
		blog = log(pow(Q*bstar/b0,2) + 1.); //modified sudakov
	else
		blog = log(pow(Q*bstar/b0,2));   //normal sudakov

	xlambda = beta0*aass*blog;
	log1xlambda = log(1.-xlambda);
	*/

	//Multiply lamB, lamC, lamD by the corresponding power of alphas to avoid overflow/uunderflow problems in the exponential
	lamB = aass*xlambda/(1.-xlambda);
	lamBqed = aaqed*xlambdaqed/(1.-xlambdaqed);
	lamC = aass2*xlambda*(xlambda-2.)/pow(1.-xlambda,2);
	lamD = aass2*log1xlambda/pow(1.-xlambda,2);
	lamE = aass3*xlambda/pow(1.-xlambda,3);
	lamF = aass3*xlambda*(3.-3.*xlambda+pow(xlambda,2))/pow(1.-xlambda,3);
	lamG = aass3*log1xlambda/pow(1.-xlambda,3);
	lamH = aass3*pow(log1xlambda,2)/pow(1.-xlambda,3);

	aexpB = exp(lamB);
	aexpBqed = exp(lamBqed);
	aexpC = exp(lamC);
	aexpD = exp(lamD);
	aexpE = exp(lamE);
	aexpF = exp(lamF);
	aexpG = exp(lamG);
	aexpH = exp(lamH);

	complex <double> aexp_bb = aexp;
	aexpB_bb = aexpB;
	double mb = pdf::mb;
	if (opts.vfncqg)
		if (mb > b0/real(bstar))
		{
		//double lamR = beta0*aass*log(pow(scales::res/mb,2));
		//aexp_bb = exp(-log(1.-lamR));

		double blog_mb = log(pow(scales::res/mb,2));
		//double blog_mb = real(btilde::blog(b0/mb,scales::res));

		//double blog_mb = real(btilde::blog(bstar,Q));

		double xlambda_bb = beta0*aass*blog_mb;
		double log1xlambda_bb = log(1.-xlambda_bb);

		double logasl_bb;
		//LL (no evolution)
		logasl_bb = 0.;

		//NLL (LO evolution)
		if (opts.order_expc >= 1)
		logasl_bb += log1xlambda_bb;

		//NNLL (NLO evolution)
		if (opts.order_expc >= 2)
		logasl_bb += aass* beta1/beta0*log1xlambda_bb/(1.-xlambda_bb);

		aexp_bb = exp(-logasl_bb);

		lamB_bb = aass*xlambda_bb/(1.-xlambda_bb);
		aexpB_bb = exp(lamB_bb);
		}



	//  if (opts.bprescription == 0)
	//    {
	//      //double lambdaqcd         = muren/(exp(1./(2.*aass*beta0))); //--> corresponds to a divergence in alphas
	//      double lambdaqcd     = mures/(exp(1./(2.*aass*beta0))); //--> correspond to a divergence in the Sudakov
	//      blim = b0/lambdaqcd;
	//      bstar = b/sqrt(1.+pow(b/blim,2));
	//      if (opts.modlog)
	//	blog = log(pow(Q*bstar/b0,2) + 1.); //modified sudakov
	//      else
	//	blog = log(pow(Q*bstar/b0,2));   //normal sudakov
	//
	//      xlambda = beta0*aass*blog;
	//      log1xlambda = log(1.-xlambda);
	//
	//      lamC = xlambda*(xlambda-2.)/pow(1.-xlambda,2);
	//      aexpC = exp(lamC);
	//    }
	//

	//double aassh = aass/2.;
	//double aasshsq = pow(aassh,2);

	complex <double> aexp2 = pow(aexp,2);
	complex <double> aexp2qed = pow(aexpqed,2);

	if (!opts.numexpc)
		switch (opts.expc)
		{
		case 0: expc::noexp(); break; 	    //No exponentiation of the C coefficients
		case 1: expc::delta(); break;         //Exponentiate only delta part
		case 2: expc::diag(); break;          //Exponentiate only flavour diagonal part
		case 3: expc::full(); fullqed(); break;          //Full expansion (standard formalism)
		case 4: expc::denlow(); break;        //Do not expand the denominator, keep only terms with C up to order-1
		case 5: expc::denfull(); break;       //Do not expand the denominator
		case 6: expc::css(b); break;          //Take C(as) at mub
		case 7: expc::matrix(); break;       //Matrix logarithm
		//case 8: expc::taylor(); break;      //Taylor expansion
		default: break;
		}

	if (opts.expc == 7)
		{
		aexp = 1.;
		aexp2 = 1.;
		}

	if (opts.expc == 6)
		{
		aexp = 1.;
		aexp2 = 1.;

		double blim = blim::expc;
		complex <double> bstar;
		if (opts.bprescription == 0 || opts.bprescription == 4 || opts.bstar_expc)
		//bstar = real(b)/sqrt(1.+pow(real(b)/blim,2));
		bstar = bstar::calc(b,blim);
		else
		bstar = b;
		complex <double> mubstar = resconst::b0/bstar;

		/*
		complex <double> mubstartilde;
		if (!opts.modlog)
		mubstartilde = mubstar;
		else if (opts.p == 1)
		mubstartilde = mubstar * scales::res / sqrt((pow(mubstar,2) + pow(scales::res,2)));
		else
		mubstartilde = mubstar * scales::res / pow((pow(mubstar,2*opts.p) + pow(scales::res,2*opts.p)),0.5/opts.p);
		*/
		complex <double> mubstartilde = btilde::mubtl(bstar,scales::res);

		//iterative alphas
		//complex <double> asmub = alphas::calc(mubstartilde, opts.order_expc+1);
		//complex <double> asmub = alphas::calc(mubstartilde, opts.order_hcoef+1);
		complex <double> asmub = alphas::calc(mubstartilde, opts.order_expc+1);
		complex <double> asmubLO = alphas::asLO;
		complex <double> asmubNLO = alphas::asNLO;


		if (opts.mellin1d)
		for (int sign = mesq::positive; sign <= mesq::negative; sign++)
		for (int i = 0; i < mellinint::mdim; i++)
			{
			int idx = anomalous::index(i,sign);
			//NLLp and NNLL
			if (opts.order_hcoef == 1)
			{
			if (sign == mesq::positive)
				{
				qqb[idx]  = (1. + asmub*(2.*C1qq[idx])+aass*H1q)/hcoeff::Hqqb[i];
				qg[idx]   = asmub*C1qg[idx]*(1.+aass*H1q)       /hcoeff::Hqg[i];
				qqbqed[idx]  = (1. + aaqed*(2.*C1qqqed[idx])+aaqed*H1qp)/hcoeff::Hqqb_qed[i];
				qgamma[idx]   = aaqed*C1qgamma[idx]*(1.+aaqed*H1qp)       /hcoeff::Hqgamma[i];
				// cout << idx << "  " << aass*C1qg[idx] << "  " << hcoeff::Hqg[i] << endl;
				}
			else
				{
				qqb[idx]  = (1. + asmub*(2.*C1qq[idx])+aass*H1q)/conj(hcoeff::Hqqb[i]);
				qg[idx]   = asmub*C1qg[idx]*(1.+aass*H1q)       /conj(hcoeff::Hqg[i]);
				qqbqed[idx]  = (1. + aaqed*(2.*C1qqqed[idx])+aaqed*H1qp)/conj(hcoeff::Hqqb_qed[i]);
				qgamma[idx]   = aaqed*C1qgamma[idx]*(1.+aaqed*H1qp)       /conj(hcoeff::Hqgamma[i]);
				//cout << idx << "  " << aass*C1qg[idx] << "  " << conj(hcoeff::Hqg[i]) << endl;
				}
			}
			//NNLLp and N3LL
			if (opts.order_hcoef == 2)
			{
			double aass2 = aass*aass;
			complex <double> asmub2 = asmub*asmub;
			if (sign == mesq::positive)
				{
				qqb[idx]  = (1.+asmub*(C1qq[idx]+C1qq[idx])+asmub2*(C2qq[idx]+C2qq[idx]+C1qq[idx]*C1qq[idx])) *(1+aass*H1q+aass2*H2q)/hcoeff::Hqqb[i];
				qg[idx]   = (asmub*C1qg[idx]+asmub2*(C2qg[idx]+C1qg[idx]*C1qq[idx]))                          *(1+aass*H1q+aass2*H2q)/hcoeff::Hqg[i];
				qq[idx]   = asmub2*C2qqb[idx]                                                                 *(1+aass*H1q+aass2*H2q)/hcoeff::Hqq[i];
				qqp[idx]  = asmub2*C2qqp[idx]                                                                 *(1+aass*H1q+aass2*H2q)/hcoeff::Hqqp[i];
				qqbp[idx] = asmub2*C2qqbp[idx]                                                                *(1+aass*H1q+aass2*H2q)/hcoeff::Hqqbp[i];
				gg[idx]   = (asmub2*(C1qg[idx]*C1qg[idx]))                                                    *(1+aass*H1q+aass2*H2q)/hcoeff::Hgg[i];
				}
			else
				{
				qqb[idx]  = (1.+asmub*(C1qq[idx]+C1qq[idx])+asmub2*(C2qq[idx]+C2qq[idx]+C1qq[idx]*C1qq[idx])) *(1+aass*H1q+aass2*H2q)/conj(hcoeff::Hqqb[i]);
				qg[idx]   = (asmub*C1qg[idx]+asmub2*(C2qg[idx]+C1qg[idx]*C1qq[idx]))                          *(1+aass*H1q+aass2*H2q)/conj(hcoeff::Hqg[i]);
				qq[idx]   = asmub2*C2qqb[idx]                                                                 *(1+aass*H1q+aass2*H2q)/conj(hcoeff::Hqq[i]);
				qqp[idx]  = asmub2*C2qqp[idx]                                                                 *(1+aass*H1q+aass2*H2q)/conj(hcoeff::Hqqp[i]);
				qqbp[idx] = asmub2*C2qqbp[idx]                                                                *(1+aass*H1q+aass2*H2q)/conj(hcoeff::Hqqbp[i]);
				gg[idx]   = (asmub2*(C1qg[idx]*C1qg[idx]))                                                    *(1+aass*H1q+aass2*H2q)/conj(hcoeff::Hgg[i]);
				}
			}
			}
		else
		for (int sign = mesq::positive; sign <= mesq::negative; sign++)
		for (int i1 = 0; i1 < mellinint::mdim; i1++)
			for (int i2 = 0; i2 < mellinint::mdim; i2++)
			{
			int idx1 = anomalous::index(i1,mesq::positive);
			int idx2 = anomalous::index(i2,sign);
			int ii12 = hcoefficients::index(i1,i2,sign);
			//NLLp and NNLL
			if (opts.order_hcoef == 1)
			{
				//qqb[ii12]    = (1. + asmub*(ccoeff::C1qq_1[idx1]+ccoeff::C1qq_2[idx2])+aass*H1q)/hcoeff::Hqqb[ii12];
				//qqb[ii12]    = (1. + asmub*ccoeff::C1qq_1[idx1])/(1. + aass*ccoeff::C1qq_1[idx1])*(1. + asmub*ccoeff::C1qq_2[idx2])/(1. + aass*ccoeff::C1qq_2[idx2]);
				//qqb[ii12]    = ((1. + asmub*ccoeff::C1qq_1[idx1])*(1. + asmub*ccoeff::C1qq_2[idx2]))/hcoeff::Hqqb[ii12];
				//qqb[ii12]    = (1. + asmub*(ccoeff::C1qq_1[idx1]+ccoeff::C1qq_2[idx2]))/hcoeff::Hqqb[ii12];//(1. + aass*(ccoeff::C1qq_1[idx1]+ccoeff::C1qq_2[idx2]+H1q));


				//qg_1[ii12]   = asmub/aass;
				//qg_2[ii12]   = asmub/aass;

				//reproduce SCET: expand beam functions, factorise Hard function
				//qqb[ii12]    = (1. + asmub*(C1qq_1[idx1]+C1qq_2[idx2]))*(1.+aass*H1q)/hcoeff::Hqqb[ii12];
				//qg_1[ii12]   = asmub*C1qg_1[idx1]                      *(1.+aass*H1q)/hcoeff::Hqg_1[ii12];
				//qg_2[ii12]   = asmub*C1qg_2[idx2]                      *(1.+aass*H1q)/hcoeff::Hqg_2[ii12];

				//reproduce NangaParbat TMD: factorise each beam function and the Hard function (works better also for SCETlib)
				qqb[ii12]    = (1.+asmub*C1qq_1[idx1])*(1.+asmub*C1qq_2[idx2])  *(1.+aass*H1q)/hcoeff::Hqqb[ii12];
				qg_1[ii12]   = (asmub*C1qg_1[idx1])*(1.+asmub*C1qq_2[idx2])     *(1.+aass*H1q)/hcoeff::Hqg_1[ii12];
				qg_2[ii12]   = (asmub*C1qg_2[idx2])*(1.+asmub*C1qq_1[idx1])     *(1.+aass*H1q)/hcoeff::Hqg_2[ii12];
				qqbqed[ii12]    = (1.+aaqed*C1qq_1[idx1]*(3./4.))*(1.+aaqed*C1qq_2[idx2]*(3./4.))  *(1.+aaqed*H1qp)/hcoeff::Hqqb_qed[ii12];
				qgamma_1[ii12]   = (aaqed*C1qg_1[idx1]*6.)*(1.+aaqed*C1qq_2[idx2]*(3./4.))     *(1.+aaqed*H1qp)/hcoeff::Hqgamma_1[ii12];
				qgamma_2[ii12]   = (aaqed*C1qg_2[idx2]*6.)*(1.+aaqed*C1qq_1[idx1]*(3./4.))     *(1.+aaqed*H1qp)/hcoeff::Hqgamma_2[ii12];
			}
			//NNLLp and N3LL
			if (opts.order_hcoef == 2)
			{
				double aass2 = aass*aass;
				complex <double> asmub2 = asmub*asmub;
				//qqb[ii12] = (1. + asmub*(ccoeff::C1qq_1[idx1]+ccoeff::C1qq_2[idx2])+aass*H1q  + asmub2*(C2qq_1[idx1] + C2qq_2[idx2] + C1qq_1[idx1]*C1qq_2[idx2]) + aass2*H2q + aass*H1q*asmub*C1qq_1[idx1] + aass*H1q*asmub*C1qq_2[idx2])/(1. + aass*(ccoeff::C1qq_1[idx1]+ccoeff::C1qq_2[idx2])+aass*H1q  + aass2*(C2qq_1[idx1] + C2qq_2[idx2] + C1qq_1[idx1]*C1qq_2[idx2]) + aass2*H2q + aass*H1q*aass*C1qq_1[idx1] + aass*H1q*aass*C1qq_2[idx2]);
				////qqb[ii12] = (1. + asmub*ccoeff::C1qq_1[idx1] + pow(asmub,2)*ccoeff::C2qq_1[idx1])/(1. + aass*ccoeff::C1qq_1[idx1] + pow(aass,2)*ccoeff::C2qq_1[idx1]) * (1. + asmub*ccoeff::C1qq_2[idx2] + pow(asmub,2)*ccoeff::C2qq_2[idx2])/(1. + aass*ccoeff::C1qq_2[idx2] + pow(aass,2)*ccoeff::C2qq_2[idx2]);
				////qg_1[ii12] = (asmub*C1qg_1[idx1] + asmub2*(C2qg_1[idx1] + C1qq_2[idx2]*C1qg_1[idx1]) + aass*H1q*asmub*C1qg_1[idx1])/(aass*C1qg_1[idx1] + aass2*(C2qg_1[idx1] + C1qq_2[idx2]*C1qg_1[idx1]) + aass*H1q*aass*C1qg_1[idx1]);
				////qg_2[ii12] = (asmub*C1qg_2[idx2] + asmub2*(C2qg_2[idx2] + C1qq_1[idx1]*C1qg_1[idx2]) + aass*H1q*asmub*C1qg_2[idx2])/(aass*C1qg_2[idx2] + aass2*(C2qg_2[idx2] + C1qq_1[idx1]*C1qg_2[idx2]) + aass*H1q*aass*C1qg_2[idx2]);
				//qg_1[ii12]   = (1. + asmub*ccoeff::C2qg_1[idx1]/ccoeff::C1qg_1[idx1])/(1. + aass*ccoeff::C2qg_1[idx1]/ccoeff::C1qg_1[idx1]) * (1. + asmub*ccoeff::C1qq_2[idx2] + pow(asmub,2)*ccoeff::C2qq_2[idx2])/(1. + aass*ccoeff::C1qq_2[idx2] + pow(aass,2)*ccoeff::C2qq_2[idx2]);
				//qg_2[ii12]   = (1. + asmub*ccoeff::C2qg_2[idx2]/ccoeff::C1qg_2[idx2])/(1. + aass*ccoeff::C2qg_2[idx2]/ccoeff::C1qg_2[idx2]) * (1. + asmub*ccoeff::C1qq_1[idx1] + pow(asmub,2)*ccoeff::C2qq_1[idx1])/(1. + aass*ccoeff::C1qq_1[idx1] + pow(aass,2)*ccoeff::C2qq_1[idx1]);
				//qq_1[ii12]  = pow(asmub,2)/pow(aass,2);
				//qq_2[ii12]  = pow(asmub,2)/pow(aass,2);
				//qqp_1[ii12]  = pow(asmub,2)/pow(aass,2);
				//qqp_2[ii12]  = pow(asmub,2)/pow(aass,2);
				//qqbp_1[ii12] = pow(asmub,2)/pow(aass,2);
				//qqbp_2[ii12] = pow(asmub,2)/pow(aass,2);
				//gg[ii12]     = pow(asmub,2)/pow(aass,2);


				//qg_1[ii12]   = aass*(C1qg_1[idx1]+asmub*C2qg_1[idx1]+H1q*asmub*C1qg_1[idx1]) * (1.+asmub*C1qq_2[idx2]+pow(asmub,2)*C2qq_2[idx2])/hcoeff::Hqg_1[ii12];
				//qg_2[ii12]   = aass*(C1qg_2[idx2]+asmub*C2qg_2[idx2]+H1q*asmub*C1qg_2[idx2]) * (1.+asmub*C1qq_1[idx1]+pow(asmub,2)*C2qq_1[idx1])/hcoeff::Hqg_2[ii12];


				//reproduce SCET: expand beam functions, factorise Hard function
				//qqb[ii12]  = (1.+asmub*(C1qq_1[idx1]+C1qq_2[idx2])+asmub2*(C2qq_1[idx1]+C2qq_2[idx2]+C1qq_1[idx1]*C1qq_2[idx2])) *(1+aass*H1q+aass2*H2q)/hcoeff::Hqqb[ii12];
				//qg_1[ii12] = (asmub*C1qg_1[idx1]+asmub2*(C2qg_1[idx1]+C1qg_1[idx1]*C1qq_2[idx2]))                                *(1+aass*H1q+aass2*H2q)/hcoeff::Hqg_1[ii12];
				//qg_2[ii12] = (asmub*C1qg_2[idx2]+asmub2*(C2qg_2[idx2]+C1qg_2[idx2]*C1qq_1[idx1]))                                *(1+aass*H1q+aass2*H2q)/hcoeff::Hqg_2[ii12];
				//qq_1[ii12]   = asmub2*C2qqb_1[idx1]                                                                              *(1+aass*H1q+aass2*H2q)/hcoeff::Hqq_1[ii12];
				//qq_2[ii12]   = asmub2*C2qqb_2[idx2]                                                                              *(1+aass*H1q+aass2*H2q)/hcoeff::Hqq_2[ii12];
				//qqp_1[ii12]  = asmub2*C2qqp_1[idx1]                                                                              *(1+aass*H1q+aass2*H2q)/hcoeff::Hqqp_1[ii12];
				//qqp_2[ii12]  = asmub2*C2qqp_2[idx2]                                                                              *(1+aass*H1q+aass2*H2q)/hcoeff::Hqqp_2[ii12];
				//qqbp_1[ii12] = asmub2*C2qqbp_1[idx1]                                                                             *(1+aass*H1q+aass2*H2q)/hcoeff::Hqqbp_1[ii12];
				//qqbp_2[ii12] = asmub2*C2qqbp_2[idx2]                                                                             *(1+aass*H1q+aass2*H2q)/hcoeff::Hqqbp_2[ii12];
				//gg[ii12]     = (asmub2*(C1qg_1[idx1]+C1qg_2[idx2]))                                                              *(1+aass*H1q+aass2*H2q)/hcoeff::Hgg[ii12];

				//reproduce NangaParbat TMD: factorise each beam function and the Hard function (works better also for SCETlib)
				qqb[ii12]  = (1.+asmub*C1qq_1[idx1]+asmub2*C2qq_1[idx1])*(1.+asmub*C1qq_2[idx2]+asmub2*C2qq_2[idx2])*(1+aass*H1q+aass2*H2q)/hcoeff::Hqqb[ii12];
				qg_1[ii12] = (asmub*C1qg_1[idx1]+asmub2*C2qg_1[idx1])*(1.+asmub*C1qq_2[idx2]+asmub2*C2qq_2[idx2])   *(1+aass*H1q+aass2*H2q)/hcoeff::Hqg_1[ii12];
				qg_2[ii12] = (asmub*C1qg_2[idx2]+asmub2*C2qg_2[idx2])*(1.+asmub*C1qq_1[idx1]+asmub2*C2qq_1[idx1])   *(1+aass*H1q+aass2*H2q)/hcoeff::Hqg_2[ii12];
				qq_1[ii12]   = asmub2*C2qqb_1[idx1]*(1.+asmub*C1qq_2[idx2]+asmub2*C2qq_2[idx2])                     *(1+aass*H1q+aass2*H2q)/hcoeff::Hqq_1[ii12];
				qq_2[ii12]   = asmub2*C2qqb_2[idx2]*(1.+asmub*C1qq_1[idx1]+asmub2*C2qq_1[idx1])                     *(1+aass*H1q+aass2*H2q)/hcoeff::Hqq_2[ii12];
				qqp_1[ii12]  = asmub2*C2qqp_1[idx1]*(1.+asmub*C1qq_2[idx2]+asmub2*C2qq_2[idx2])                     *(1+aass*H1q+aass2*H2q)/hcoeff::Hqqp_1[ii12];
				qqp_2[ii12]  = asmub2*C2qqp_2[idx2]*(1.+asmub*C1qq_1[idx1]+asmub2*C2qq_1[idx1])                     *(1+aass*H1q+aass2*H2q)/hcoeff::Hqqp_2[ii12];
				qqbp_1[ii12] = asmub2*C2qqbp_1[idx1]*(1.+asmub*C1qq_2[idx2]+asmub2*C2qq_2[idx2])                    *(1+aass*H1q+aass2*H2q)/hcoeff::Hqqbp_1[ii12];
				qqbp_2[ii12] = asmub2*C2qqbp_2[idx2]*(1.+asmub*C1qq_1[idx1]+asmub2*C2qq_1[idx1])                    *(1+aass*H1q+aass2*H2q)/hcoeff::Hqqbp_2[ii12];
				gg[ii12]     = (asmub*C1qg_1[idx1]+asmub2*C2qg_1[idx1])*(asmub*C1qg_2[idx2]+asmub2*C2qg_2[idx2])    *(1+aass*H1q+aass2*H2q)/hcoeff::Hgg[ii12];
			}
			}
		return;
		}


	if (opts.numexpc)
		if (opts.mellin1d)
		for (int sign = mesq::positive; sign <= mesq::negative; sign++)
		for (int i = 0; i < mellinint::mdim; i++)
		{
			int idx = anomalous::index(i,sign);
			aexpqq[idx]   = exp(gint::alogqq[idx]);
			aexpqg[idx]   = exp(gint::alogqg[idx]);
			aexpqqqed[idx]   = exp(gint::alogqq[idx]); // QED
			aexpqgamma[idx]   = exp(gint::alogqg[idx]); //QED
			aexpqqb[idx]  = exp(gint::alogqqb[idx]);
			aexpqqp[idx]  = exp(gint::alogqqp[idx]);
			aexpqqbp[idx] = exp(gint::alogqqbp[idx]);
		}
		else
		for (int sign = mesq::positive; sign <= mesq::negative; sign++)
		for (int i = 0; i < mellinint::mdim; i++)
		{
			int idx = anomalous::index(i,sign);
			aexpqq_1[idx]   = exp(gint::alogqq_1[idx]);
			aexpqg_1[idx]   = exp(gint::alogqg_1[idx]);
			aexpqqqed_1[idx]   = exp(gint::alogqq_1[idx]); // QED not ready
			aexpqgamma_1[idx]   = exp(gint::alogqg_1[idx]);// QED not ready
			aexpqqb_1[idx]  = exp(gint::alogqqb_1[idx]);
			aexpqqp_1[idx]  = exp(gint::alogqqp_1[idx]);
			aexpqqbp_1[idx] = exp(gint::alogqqbp_1[idx]);
			aexpqq_2[idx]   = exp(gint::alogqq_2[idx]);
			aexpqg_2[idx]   = exp(gint::alogqg_2[idx]);
			aexpqqqed_2[idx]   = exp(gint::alogqq_2[idx]);// QED not ready
			aexpqgamma_2[idx]   = exp(gint::alogqg_2[idx]);// QED not ready
			aexpqqb_2[idx]  = exp(gint::alogqqb_2[idx]);
			aexpqqp_2[idx]  = exp(gint::alogqqp_2[idx]);
			aexpqqbp_2[idx] = exp(gint::alogqqbp_2[idx]);
		}

	//include aexp into aexpqg and aexpqq,qqp,qqb,qqbp
	//if (opts.expc > 0)
	if (opts.mellin1d)
		for (int sign = mesq::positive; sign <= mesq::negative; sign++)
		for (int i = 0; i < mellinint::mdim; i++)
		{
		int idx = anomalous::index(i,sign);
		aexpqg[idx]       *= aexp;
		aexpqgamma[idx]   *= aexpqed;
		aexpqqb[idx]      *= aexp2;
		aexpqqp[idx]      *= aexp2;
		aexpqqbp[idx]     *= aexp2;
		aexpqg_bb[idx]   *= aexp_bb;
		}
	else
		for (int sign = mesq::positive; sign <= mesq::negative; sign++)
		for (int i = 0; i < mellinint::mdim; i++)
		{
		int idx = anomalous::index(i,sign);
		aexpqg_1[idx]   *= aexp;
		aexpqgamma_1[idx]   *= aexpqed;
		aexpqqb_1[idx]  *= aexp2;
		aexpqqp_1[idx]  *= aexp2;
		aexpqqbp_1[idx] *= aexp2;
		aexpqg_2[idx]   *= aexp;
		aexpqgamma_2[idx]   *= aexpqed;
		aexpqqb_2[idx]  *= aexp2;
		aexpqqp_2[idx]  *= aexp2;
		aexpqqbp_2[idx] *= aexp2;

		aexpqg_1_bb[idx]   *= aexp_bb;
		aexpqg_2_bb[idx]   *= aexp_bb;
		}

	//include DMW
	if (opts.npff == 3)
		{
		//keep the N-independent part with b real
		complex <double> bb;
		bb = real(b); //b; //abs(b);
		//complex <double> dmw_fix = exp(-(opts.a2 * (2.*log(phasespace::m)+1.) - opts.a2p) * pow(bb,2)/2.);
		complex <double> dmw_fix = 1.;

		//Need a bstar with (large) blim otherwise values at large N will overflow
		double bmx = 30.;//blim::sudakov;
		complex <double> bst;
		if (opts.bprescription == 0 || opts.bprescription == 4 || opts.bstar_expc)
		bst = real(b)/sqrt(1.+pow(real(b)/bmx,2));
		//bstar = bstar::calc(b,bmx);
		else
		bst = b;

		if (opts.mellin1d)
		for (int sign = mesq::positive; sign <= mesq::negative; sign++)
		for (int i = 0; i < mellinint::mdim; i++)
			{
			complex <double> N;
			if (sign == mesq::positive)
			N = mellinint::Np[i];
			else
			N = mellinint::Nm[i];
			//complex <double> S1N = 0.;
			complex <double> S1N = cpsi0(N+1.)+constants::euler;
			//complex <double> S1N = cpsi0(N)+constants::euler;
			complex <double> dmw = exp(-pow(bst,2)/2.*(opts.a2*(-2.*S1N-(4.*N+4.)/(N*(N+2.)))));
			int idx = anomalous::index(i,sign);
			aexpqq[idx]   *= dmw*dmw_fix;
			//Include the N-independent part in the other channels
			aexpqg[idx]   *= dmw_fix;
			aexpqqb[idx]  *= dmw_fix;
			aexpqqp[idx]  *= dmw_fix;
			aexpqqbp[idx] *= dmw_fix;
			}
		else
		for (int sign = mesq::positive; sign <= mesq::negative; sign++)
		for (int i = 0; i < mellinint::mdim; i++)
			{
			complex <double> N_1,N_2;
			if (sign == mesq::positive)
			{
			N_1 = mellinint::Np_1[i];
			N_2 = mellinint::Np_2[i];
			}
			else
			{
			N_1 = mellinint::Nm_1[i];
			N_2 = mellinint::Nm_2[i];
			}
			complex <double> S1N_1 = cpsi0(N_1+1.)+constants::euler;
			complex <double> S1N_2 = cpsi0(N_2+1.)+constants::euler;
			//complex <double> S1N_1 = cpsi0(N_1)+constants::euler;
			//complex <double> S1N_2 = cpsi0(N_2)+constants::euler;
			complex <double> dmw_1 = exp(-pow(bst,2)/2.*(opts.a2*(-2.*S1N_1-(4.*N_1+4.)/(N_1*(N_1+2.)))));
			complex <double> dmw_2 = exp(-pow(bst,2)/2.*(opts.a2*(-2.*S1N_2-(4.*N_2+4.)/(N_2*(N_2+2.)))));
			int idx = anomalous::index(i,sign);
			aexpqq_1[idx]   *= dmw_1*dmw_fix;
			aexpqq_2[idx]   *= dmw_2*dmw_fix;
			//Include the N-independent part in the other channels
			aexpqg_1[idx]   *= dmw_fix;
			aexpqqb_1[idx]  *= dmw_fix;
			aexpqqp_1[idx]  *= dmw_fix;
			aexpqqbp_1[idx] *= dmw_fix;
			aexpqg_2[idx]   *= dmw_fix;
			aexpqqb_2[idx]  *= dmw_fix;
			aexpqqp_2[idx]  *= dmw_fix;
			aexpqqbp_2[idx] *= dmw_fix;
			}
		}

	if (opts.mellin1d)
		for (int sign = mesq::positive; sign <= mesq::negative; sign++)
		for (int i = 0; i < mellinint::mdim; i++)
		{
		int idx = anomalous::index(i,sign);
		int ii = hcoefficients::index(i,sign);
		//NLL
		qqb[ii] = aexpqq[idx] * aexpqq[idx];
		qg[ii]  = aexpqg[idx] * aexpqq[idx];
		qqbqed[ii] = aexpqqqed[idx] * aexpqqqed[idx];
		qgamma[ii]  = aexpqgamma[idx] * aexpqqqed[idx];
		//NNLL
		qq[ii]   = aexpqq[idx] * aexpqqb[idx];
		qqp[ii]  = aexpqq[idx] * aexpqqbp[idx];
		qqbp[ii] = aexpqq[idx] * aexpqqp[idx];
		gg[ii]   = aexpqg[idx] * aexpqg[idx];
		gg_bb[ii]= aexpqg_bb[idx] * aexpqg_bb[idx];
		//N3LL
		qbg[ii]  = aexpqg[idx] * aexpqqb[idx];
		qpg[ii]  = aexpqg[idx] * aexpqqp[idx];
		qbpg[ii] = aexpqg[idx] * aexpqqbp[idx];
		//N4LL
		qbq[ii]   = aexpqqb[idx] * aexpqqb[idx];
		qbqp[ii]  = aexpqqb[idx] * aexpqqbp[idx];
		qbqbp[ii] = aexpqqb[idx] * aexpqqp[idx];
		qpqp[ii]  = aexpqqp[idx] * aexpqqbp[idx];
		qpqbp[ii] = aexpqqp[idx] * aexpqqp[idx];
		qbpqp[ii] = aexpqqp[idx] * aexpqqp[idx];
		}
	else
		for (int sign = mesq::positive; sign <= mesq::negative; sign++)
		for (int i1 = 0; i1 < mellinint::mdim; i1++)
		for (int i2 = 0; i2 < mellinint::mdim; i2++)
		{
			int idx1 = anomalous::index(i1,mesq::positive);
			int idx2 = anomalous::index(i2,sign);
			int ii12 = hcoefficients::index(i1,i2,sign);
			//NLL
			qqb[ii12]    = aexpqq_1[idx1] * aexpqq_2[idx2];
			qg_1[ii12]   = aexpqg_1[idx1] * aexpqq_2[idx2];
			qg_2[ii12]   = aexpqg_2[idx2] * aexpqq_1[idx1];
			qqbqed[ii12]    = aexpqqqed_1[idx1] * aexpqqqed_2[idx2];
			qgamma_1[ii12]   = aexpqgamma_1[idx1] * aexpqqqed_2[idx2];
			qgamma_2[ii12]   = aexpqgamma_2[idx2] * aexpqqqed_1[idx1];
			//NNLL

			// //bug bug bug!!! _1 and _2 should be inverted
			// qq_1[ii12]   = aexpqq_1[idx1] * aexpqqb_2[idx2];
			// qq_2[ii12]   = aexpqq_2[idx2] * aexpqqb_1[idx1];
			// qqp_1[ii12]  = aexpqq_1[idx1] * aexpqqbp_2[idx2];
			// qqp_2[ii12]  = aexpqq_2[idx2] * aexpqqbp_1[idx1];
			// qqbp_1[ii12] = aexpqq_1[idx1] * aexpqqp_2[idx2];
			// qqbp_2[ii12] = aexpqq_2[idx2] * aexpqqp_1[idx1];

			//bug fix
			qq_1[ii12]   = aexpqqb_1[idx1]  * aexpqq_2[idx2];
			qq_2[ii12]   = aexpqqb_2[idx2]  * aexpqq_1[idx1];
			qqp_1[ii12]  = aexpqqbp_1[idx1] * aexpqq_2[idx2];
			qqp_2[ii12]  = aexpqqbp_2[idx2] * aexpqq_1[idx1];
			qqbp_1[ii12] = aexpqqp_1[idx1]  * aexpqq_2[idx2];
			qqbp_2[ii12] = aexpqqp_2[idx2]  * aexpqq_1[idx1];

			gg[ii12]     = aexpqg_1[idx1] * aexpqg_2[idx2];
			gg_bb[ii12]  = aexpqg_1_bb[idx1] * aexpqg_2_bb[idx2];
			//N3LL
			qbg_1[ii12]  = aexpqg_1[idx1] * aexpqqb_2[idx2];
			qbg_2[ii12]  = aexpqg_2[idx2] * aexpqqb_1[idx1];
			qpg_1[ii12]  = aexpqg_1[idx1] * aexpqqp_2[idx2];
			qpg_2[ii12]  = aexpqg_2[idx2] * aexpqqp_1[idx1];
			qbpg_1[ii12] = aexpqg_1[idx1] * aexpqqbp_2[idx2];
			qbpg_2[ii12] = aexpqg_2[idx2] * aexpqqbp_1[idx1];
		}

	//new code implemented for most options
	if (opts.expc == 0)
		return;
	if (opts.expc == 1)
		return;
	if (opts.expc == 2)
		return;
	if (opts.expc == 3)
		return;
	if (opts.expc == 4)
		return;
	if (opts.expc == 5)
		return;
	if (opts.expc == 6)
		return;
	if (opts.expc == 7)
		return;


	// NLL
	if (opts.order_expc == 1)
		if (opts.mellin1d)
		if (opts.sumlogs)
		fill(qg, qg+mellinint::mdim*2, exp(-logasl+sudakov::logS));
		else
		{
		fill(aexpqg, aexpqg+mellinint::mdim*2, aexp);
		fill(qg, qg+mellinint::mdim*2, aexp);
		}
		else
		if (opts.sumlogs)
		{
		fill(qg_1, qg_1+mellinint::mdim*mellinint::mdim*2, exp(-logasl+sudakov::logS));
		fill(qg_2, qg_2+mellinint::mdim*mellinint::mdim*2, exp(-logasl+sudakov::logS));
		}
		else
		{
		fill(aexpqg_1, aexpqg_1+mellinint::mdim*2, aexp);
		fill(aexpqg_2, aexpqg_2+mellinint::mdim*2, aexp);
		fill(qg_1, qg_1+mellinint::mdim*mellinint::mdim*2, aexp);
		fill(qg_2, qg_2+mellinint::mdim*mellinint::mdim*2, aexp);
		}

	//gint::calc(b);
	//cout << endl;
	//cout << " gint::logasl_expc "   << gint::logasl_expc      << " expc " << logasl      << endl;
	// NNLL
	//complex <double> c1delta = pow(aexpb,aassh*(-2.*C1qqn));
	if (opts.order_expc == 2)
		{
		if (opts.mellin1d)
		{
		//Use pointers to allocate memory on the heap
		complex <double> *alogqq   = new complex <double> [mellinint::mdim*2];
		complex <double> *alogqg   = new complex <double> [mellinint::mdim*2];
		complex <double> *alogqqb  = new complex <double> [mellinint::mdim*2];
		complex <double> *alogqqp  = new complex <double> [mellinint::mdim*2];
		complex <double> *alogqqbp = new complex <double> [mellinint::mdim*2];
		//complex <double> *aexpqq   = new complex <double> [mellinint::mdim*2];
		//complex <double> *aexpqg   = new complex <double> [mellinint::mdim*2];
		//complex <double> *aexpqqb  = new complex <double> [mellinint::mdim*2];
		//complex <double> *aexpqqp  = new complex <double> [mellinint::mdim*2];
		//complex <double> *aexpqqbp = new complex <double> [mellinint::mdim*2];
		//complex <double> alogqq[mellinint::mdim*2];
		//complex <double> alogqg[mellinint::mdim*2];
		//complex <double> alogqqb[mellinint::mdim*2];
		//complex <double> alogqqp[mellinint::mdim*2];
		//complex <double> alogqqbp[mellinint::mdim*2];
		//complex <double> aexpqq[mellinint::mdim*2];
		//complex <double> aexpqg[mellinint::mdim*2];
		//complex <double> aexpqqb[mellinint::mdim*2];
		//complex <double> aexpqqp[mellinint::mdim*2];
		//complex <double> aexpqqbp[mellinint::mdim*2];

		//No exponentiation of the C coefficients
		if (opts.expc == 0)
			if (opts.sumlogs)
			{
			fill(alogqq, alogqq+mellinint::mdim*2, 0.);
			fill(alogqg, alogqq+mellinint::mdim*2, 0.);
			fill(alogqqb, alogqqb+mellinint::mdim*2, 0.);
			fill(alogqqp, alogqqp+mellinint::mdim*2, 0.);
			fill(alogqqbp, alogqqbp+mellinint::mdim*2, 0.);
			}
			else
			{
			fill(aexpqq, aexpqq+mellinint::mdim*2, 1.);
			fill(aexpqg, aexpqg+mellinint::mdim*2, 1.);
			fill(aexpqqb, aexpqqb+mellinint::mdim*2, 1.);
			fill(aexpqqp, aexpqqp+mellinint::mdim*2, 1.);
			fill(aexpqqbp, aexpqqbp+mellinint::mdim*2, 1.);
			}
		//Only delta terms
		else if (opts.expc == 1)
			if (opts.sumlogs)
			{
			fill(alogqg, alogqq+mellinint::mdim*2, 0.);
			fill(alogqqb, alogqqb+mellinint::mdim*2, 0.);
			fill(alogqqp, alogqqp+mellinint::mdim*2, 0.);
			fill(alogqqbp, alogqqbp+mellinint::mdim*2, 0.);
			for (int sign = mesq::positive; sign <= mesq::negative; sign++)
			for (int i = 0; i < mellinint::mdim; i++)
				{
				int idx = anomalous::index(i,sign);
				alogqq[idx] = lamB*ccoeff::C1qq_delta;
				}
			}
			else
			{
			fill(aexpqg, aexpqg+mellinint::mdim*2, 1.);
			fill(aexpqqb, aexpqqb+mellinint::mdim*2, 1.);
			fill(aexpqqp, aexpqqp+mellinint::mdim*2, 1.);
			fill(aexpqqbp, aexpqqbp+mellinint::mdim*2, 1.);
			for (int sign = mesq::positive; sign <= mesq::negative; sign++)
			for (int i = 0; i < mellinint::mdim; i++)
				{
				int idx = anomalous::index(i,sign);
				aexpqq[idx] = pow(aexpB,ccoeff::C1qq_delta);
				}
			}
		//Taylor expansion
		else if (opts.expc == 2)
			for (int sign = mesq::positive; sign <= mesq::negative; sign++)
			for (int i = 0; i < mellinint::mdim; i++)
			{
			int idx = anomalous::index(i,sign);
			complex <double> qqtayl = 1.;
			complex <double> qgtayl = 1.;
			complex <double> qqbtayl = 1.;
			complex <double> qqptayl = 1.;
			complex <double> qqbptayl = 1.;
			int fact = 1;
			for (int t = 1; t <= opts.ntaylor; t++)
				{
				fact *= t;
				qqtayl += 1./double(fact) *pow(lamB*(ccoeff::C1qq[idx]-ccoeff::C1qq_delta),t);
				qgtayl += 1./double(fact) *pow(lamB*ccoeff::C2qg[idx]/ccoeff::C1qg[idx],t);
				qqbtayl += 1./double(fact) *pow(lamB*ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx],t);
				qqptayl += 1./double(fact) *pow(lamB*ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx],t);
				qqbptayl += 1./double(fact) *pow(lamB*ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx],t);
				}
			if (opts.sumlogs)
				{
				alogqq[idx] = lamB*ccoeff::C1qq_delta + log(qqtayl);
				alogqg[idx] = log(qgtayl);
				alogqqb[idx] = log(qqbtayl);
				alogqqp[idx] = log(qqptayl);
				alogqqbp[idx] = log(qqbptayl);
				}
			else
				{
				aexpqq[idx] = pow(aexpB,ccoeff::C1qq_delta)*qqtayl;
				aexpqg[idx] = qgtayl;
				aexpqqb[idx] = qqbtayl;
				aexpqqp[idx] = qqptayl;
				aexpqqbp[idx] = qqbptayl;
				}
			}
		//Standard formula
		else if (opts.expc == 3)
			for (int sign = mesq::positive; sign <= mesq::negative; sign++)
			for (int i = 0; i < mellinint::mdim; i++)
			{
			int idx = anomalous::index(i,sign);
			if (opts.sumlogs)
				{
				//alogqq = lamB*1./2.*(C1QQ[idx]-2.*C1qqn);
				//alogqg = lamB*1./2.*(C2qgM[idx]/C1QG[idx]-2.*C1qqn);
				alogqq[idx] = lamB*ccoeff::C1qq[idx];
				alogqg[idx] = lamB*ccoeff::C2qg[idx]/ccoeff::C1qg[idx];
				alogqqb[idx] = lamB*ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx];
				alogqqp[idx] = lamB*ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx];
				alogqqbp[idx] = lamB*ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx];
				}
			else
				{
				//aexpqq = pow(aexpB,1./2.*(C1QQ[idx]-2.*C1qqn));
				//aexpqg = pow(aexpB,1./2.*(C2qgM[idx]/C1QG[idx]-2.*C1qqn));
				aexpqq[idx] = pow(aexpB,ccoeff::C1qq[idx]);
				aexpqg[idx] = pow(aexpB,ccoeff::C2qg[idx]/ccoeff::C1qg[idx]);
				aexpqqb[idx] = pow(aexpB,ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx]);
				aexpqqp[idx] = pow(aexpB,ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx]);
				aexpqqbp[idx] = pow(aexpB,ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx]);
				}
			}
		//cout << endl;
		//cout << aexpqq[0]   << " expc::aexpqq " << endl;
		//cout << aexpqg[0]   << " expc::aexpqg " << endl;
		//cout << aexpqqb[0]  << " expc::aexpqqb " << endl;
		//cout << aexpqqp[0]  << " expc::aexpqqb " << endl;
		//cout << aexpqqbp[0] << " expc::aexpqqbp " << endl;
		//cout << endl;
		//At order n expand up to terms containing Cn (i.e. expand up to n for D, up to n-1 for OD, and up to n-2 for DOD)
		else if (opts.expc == 4)
			{
			if (opts.sumlogs)
			{
			fill(alogqg, alogqq+mellinint::mdim*2, 0.);
			fill(alogqqb, alogqqb+mellinint::mdim*2, 0.);
			fill(alogqqp, alogqqp+mellinint::mdim*2, 0.);
			fill(alogqqbp, alogqqbp+mellinint::mdim*2, 0.);
			}
			else
			{
			fill(aexpqg, aexpqg+mellinint::mdim*2, 1.);
			fill(aexpqqb, aexpqqb+mellinint::mdim*2, 1.);
			fill(aexpqqp, aexpqqp+mellinint::mdim*2, 1.);
			fill(aexpqqbp, aexpqqbp+mellinint::mdim*2, 1.);
			}
			for (int sign = mesq::positive; sign <= mesq::negative; sign++)
			for (int i = 0; i < mellinint::mdim; i++)
			{
				int idx = anomalous::index(i,sign);
				if (opts.sumlogs)
				alogqq[idx] = lamB*ccoeff::C1qq[idx];
				else
				aexpqq[idx] = pow(aexpB,ccoeff::C1qq[idx]);
			}
			}
		//Do not expand C at the denominator
		else if (opts.expc == 5)
			for (int sign = mesq::positive; sign <= mesq::negative; sign++)
			for (int i = 0; i < mellinint::mdim; i++)
			{
			int idx = anomalous::index(i,sign);
			if (opts.sumlogs)
				{
				alogqq[idx] = lamB*ccoeff::C1qq[idx];
				alogqg[idx]   = log(1.-xlambda/(1.+aass*ccoeff::C2qg[idx]/ccoeff::C1qg[idx])) - log(1.-xlambda);
				alogqqb[idx]  = log(1.-xlambda/(1.+aass*ccoeff::C3qqb[idx] /ccoeff::C2qqb[idx] )) - log(1.-xlambda);
				alogqqp[idx]  = log(1.-xlambda/(1.+aass*ccoeff::C3qqp[idx] /ccoeff::C2qqp[idx] )) - log(1.-xlambda);
				alogqqbp[idx] = log(1.-xlambda/(1.+aass*ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx])) - log(1.-xlambda);
				//alternative with as0 at the denominator
				//alogqg[idx]   = lamB*ccoeff::C2qg[idx]  /ccoeff::C1qg[idx]  /(1.+aass*ccoeff::C2qg[idx]  /ccoeff::C1qg[idx]  );
				//alogqqb[idx]  = lamB*ccoeff::C3qqb[idx] /ccoeff::C2qqb[idx] /(1.+aass*ccoeff::C3qqb[idx] /ccoeff::C2qqb[idx] );
				//alogqqp[idx]  = lamB*ccoeff::C3qqp[idx] /ccoeff::C2qqp[idx] /(1.+aass*ccoeff::C3qqp[idx] /ccoeff::C2qqp[idx] );
				//alogqqbp[idx] = lamB*ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx]/(1.+aass*ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx]);

				}
			else
				{
				aexpqq[idx] = pow(aexpB,ccoeff::C1qq[idx]);
				aexpqg[idx]   = (1.-xlambda/(1.+aass*ccoeff::C2qg[idx]/ccoeff::C1qg[idx]))/(1.-xlambda);
				aexpqqb[idx]  = (1.-xlambda/(1.+aass*ccoeff::C3qqb[idx] /ccoeff::C2qqb[idx] ))/(1.-xlambda);
				aexpqqp[idx]  = (1.-xlambda/(1.+aass*ccoeff::C3qqp[idx] /ccoeff::C2qqp[idx] ))/(1.-xlambda);
				aexpqqbp[idx] = (1.-xlambda/(1.+aass*ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx]))/(1.-xlambda);

				}
			}
		//gint::calc(b);
		//cout << endl;
		//cout << " gint::logasl "   << gint::logasl_expc      << " expc " << logasl      << endl;
		//cout << " gint::alogqq "   << gint::alogqq[0]   << " expc " << alogqq[0]   << endl;
		//cout << " gint::alogqg "   << gint::alogqg[0]   << " expc " << alogqg[0]   << endl;
		//cout << " gint::alogqqb "  << gint::alogqqb[0]  << " expc " << alogqqb[0]  << endl;
		//cout << " gint::alogqqp "  << gint::alogqqp[0]  << " expc " << alogqqp[0]  << endl;
		//cout << " gint::alogqqbp " << gint::alogqqbp[0] << " expc " << alogqqbp[0] << endl;

		//cout << aexpqq[0]   << " gint::alogqq " << endl;
		//cout << aexpqg[0]   << " gint::alogqg " << endl;
		//cout << aexpqqb[0]  << " gint::alogqqb " << endl;
		//cout << aexpqqp[0]  << " gint::alogqqb " << endl;
		//cout << aexpqqbp[0] << " gint::alogqqbp " << endl;
		//cout << endl;

		//!!! gint::alog are evaluated with blim::sudakov, not with blim:expc !!!
		if (opts.numexpc)
			for (int sign = mesq::positive; sign <= mesq::negative; sign++)
			for (int i = 0; i < mellinint::mdim; i++)
			{
			int idx = anomalous::index(i,sign);

			if (opts.sumlogs)
				{
				alogqq[idx]   = gint::alogqq[idx];
				alogqg[idx]   = gint::alogqg[idx];
				alogqqb[idx]  = gint::alogqqb[idx];
				alogqqp[idx]  = gint::alogqqp[idx];
				alogqqbp[idx] = gint::alogqqbp[idx];
				}
			else
				{
				aexpqq[idx]   = exp(gint::alogqq[idx]);
				aexpqg[idx]   = exp(gint::alogqg[idx]);
				aexpqqb[idx]  = exp(gint::alogqqb[idx]);
				aexpqqp[idx]  = exp(gint::alogqqp[idx]);
				aexpqqbp[idx] = exp(gint::alogqqbp[idx]);
				}
			}

		for (int sign = mesq::positive; sign <= mesq::negative; sign++)
			for (int i = 0; i < mellinint::mdim; i++)
			{
			int idx = anomalous::index(i,sign);
			int ii = hcoefficients::index(i,sign);
			if (opts.sumlogs)
			{
				qqb[ii] = exp(alogqq[idx]+alogqq[idx]+sudakov::logS);
				qg[ii]  = exp(-logasl+alogqg[idx]+alogqq[idx]+sudakov::logS);
				qq[ii]  = exp(alogqq[idx]+alogqqb[idx]-2.*logasl+sudakov::logS);
				qqp[ii] = exp(alogqq[idx]+alogqqbp[idx]-2.*logasl+sudakov::logS);
				qqbp[ii] = exp(alogqq[idx]+alogqqp[idx]-2.*logasl+sudakov::logS);
				gg[ii]  = exp(-logasl+alogqg[idx]-logasl+alogqg[idx]+sudakov::logS);
			}
			else
			{
				qqb[ii] = aexpqq[idx] * aexpqq[idx];
				qg[ii]  = aexp * aexpqg[idx]*aexpqq[idx];
				qq[ii]  = aexpqq[idx] * aexpqqb[idx] * aexp2; //*c1delta
				qqp[ii] = aexpqq[idx] * aexpqqbp[idx] * aexp2; //*c1delta
				qqbp[ii] = aexpqq[idx] * aexpqqp[idx] * aexp2; //*c1delta
				gg[ii]  = aexp*aexpqg[idx] * aexp*aexpqg[idx];
			}
			}

		delete[] alogqq;
		delete[] alogqg;
		delete[] alogqqb;
		delete[] alogqqp;
		delete[] alogqqbp;

		//include aexp into aexpqg and aexpqq,qqp,qqb,qqbp
		for (int sign = mesq::positive; sign <= mesq::negative; sign++)
			for (int i = 0; i < mellinint::mdim; i++)
			{
			int idx = anomalous::index(i,sign);
			aexpqg[idx]   *= aexp;
			aexpqqb[idx]  *= aexp2;
			aexpqqp[idx]  *= aexp2;
			aexpqqbp[idx] *= aexp2;
			}
		//delete[] aexpqq;
		//delete[] aexpqg;
		//delete[] aexpqqb;
		//delete[] aexpqqp;
		//delete[] aexpqqbp;
		}
		else //mellin2d
		{
		//complex <double> aexpqq_1[mellinint::mdim*2];
		//complex <double> aexpqg_1[mellinint::mdim*2];
		//complex <double> aexpqqb_1[mellinint::mdim*2];
		//complex <double> aexpqqp_1[mellinint::mdim*2];
		//complex <double> aexpqqbp_1[mellinint::mdim*2];
		//complex <double> aexpqq_2[mellinint::mdim*2];
		//complex <double> aexpqg_2[mellinint::mdim*2];
		//complex <double> aexpqqb_2[mellinint::mdim*2];
		//complex <double> aexpqqp_2[mellinint::mdim*2];
		//complex <double> aexpqqbp_2[mellinint::mdim*2];

		//case opts.expc == 0
		fill(aexpqq_1, aexpqq_1+mellinint::mdim*2, 1.);
		fill(aexpqg_1, aexpqg_1+mellinint::mdim*2, 1.);
		fill(aexpqqqed_1, aexpqqqed_1+mellinint::mdim*2, 1.);
		fill(aexpqgamma_1, aexpqgamma_1+mellinint::mdim*2, 1.);
		fill(aexpqqb_1, aexpqqb_1+mellinint::mdim*2, 1.);
		fill(aexpqqp_1, aexpqqp_1+mellinint::mdim*2, 1.);
		fill(aexpqqbp_1, aexpqqbp_1+mellinint::mdim*2, 1.);
		fill(aexpqq_2, aexpqq_2+mellinint::mdim*2, 1.);
		fill(aexpqg_2, aexpqg_2+mellinint::mdim*2, 1.);
		fill(aexpqqqed_2, aexpqqqed_2+mellinint::mdim*2, 1.);
		fill(aexpqgamma_2, aexpqgamma_2+mellinint::mdim*2, 1.);
		fill(aexpqqb_2, aexpqqb_2+mellinint::mdim*2, 1.);
		fill(aexpqqp_2, aexpqqp_2+mellinint::mdim*2, 1.);
		fill(aexpqqbp_2, aexpqqbp_2+mellinint::mdim*2, 1.);

		//Only delta terms
		if (opts.expc == 1)
			for (int sign = mesq::positive; sign <= mesq::negative; sign++)
			for (int i = 0; i < mellinint::mdim; i++)
			{
			int idx = anomalous::index(i,sign);
			aexpqq_1[idx] = pow(aexpB,ccoeff::C1qq_delta);
			aexpqq_2[idx] = pow(aexpB,ccoeff::C1qq_delta);
			}
		//Taylor expansion
		else if (opts.expc == 2)
			{
			for (int sign = mesq::positive; sign <= mesq::negative; sign++)
			for (int i = 0; i < mellinint::mdim; i++)
			{
				int idx = anomalous::index(i,sign);
				complex <double> qqtayl_1 = 1.;
				complex <double> qgtayl_1 = 1.;
				complex <double> qqbtayl_1 = 1.;
				complex <double> qqptayl_1 = 1.;
				complex <double> qqbptayl_1 = 1.;
				complex <double> qqtayl_2 = 1.;
				complex <double> qgtayl_2 = 1.;
				complex <double> qqbtayl_2 = 1.;
				complex <double> qqptayl_2 = 1.;
				complex <double> qqbptayl_2 = 1.;
				int fact = 1;
				for (int t = 1; t <= opts.ntaylor; t++)
				{
				fact *= t;
				qqtayl_1 += 1./double(fact) *pow(lamB*(ccoeff::C1qq_1[idx]-ccoeff::C1qq_delta),t);
				qgtayl_1 += 1./double(fact) *pow(lamB*ccoeff::C2qg_1[idx]/ccoeff::C1qg_1[idx],t);
				qqbtayl_1 += 1./double(fact) *pow(lamB*ccoeff::C3qqb_1[idx]/ccoeff::C2qqb_1[idx],t);
				qqptayl_1 += 1./double(fact) *pow(lamB*ccoeff::C3qqp_1[idx]/ccoeff::C2qqp_1[idx],t);
				qqbptayl_1 += 1./double(fact) *pow(lamB*ccoeff::C3qqbp_1[idx]/ccoeff::C2qqbp_1[idx],t);
				qqtayl_2 += 1./double(fact) *pow(lamB*(ccoeff::C1qq_2[idx]-ccoeff::C1qq_delta),t);
				qgtayl_2 += 1./double(fact) *pow(lamB*ccoeff::C2qg_2[idx]/ccoeff::C1qg_2[idx],t);
				qqbtayl_2 += 1./double(fact) *pow(lamB*ccoeff::C3qqb_2[idx]/ccoeff::C2qqb_2[idx],t);
				qqptayl_2 += 1./double(fact) *pow(lamB*ccoeff::C3qqp_2[idx]/ccoeff::C2qqp_2[idx],t);
				qqbptayl_2 += 1./double(fact) *pow(lamB*ccoeff::C3qqbp_2[idx]/ccoeff::C2qqbp_2[idx],t);
				}
				aexpqq_1[idx] = pow(aexpB,ccoeff::C1qq_delta)*qqtayl_1;
				aexpqg_1[idx] = qgtayl_1;
				aexpqqb_1[idx] = qqbtayl_1;
				aexpqqp_1[idx] = qqptayl_1;
				aexpqqbp_1[idx] = qqbptayl_1;
				aexpqq_2[idx] = pow(aexpB,ccoeff::C1qq_delta)*qqtayl_2;
				aexpqg_2[idx] = qgtayl_2;
				aexpqqb_2[idx] = qqbtayl_2;
				aexpqqp_2[idx] = qqptayl_2;
				aexpqqbp_2[idx] = qqbptayl_2;
			}
			}
		//Standard expression
		else if (opts.expc == 3)
			for (int sign = mesq::positive; sign <= mesq::negative; sign++)
			for (int i = 0; i < mellinint::mdim; i++)
			{
			int idx = anomalous::index(i,sign);
			aexpqq_1[idx] = pow(aexpB,ccoeff::C1qq_1[idx]);
			aexpqg_1[idx] = pow(aexpB,ccoeff::C2qg_1[idx]/ccoeff::C1qg_1[idx]);
			aexpqqb_1[idx] = pow(aexpB,ccoeff::C3qqb_1[idx]/ccoeff::C2qqb_1[idx]);
			aexpqqp_1[idx] = pow(aexpB,ccoeff::C3qqp_1[idx]/ccoeff::C2qqp_1[idx]);
			aexpqqbp_1[idx] = pow(aexpB,ccoeff::C3qqbp_1[idx]/ccoeff::C2qqbp_1[idx]);
			aexpqq_2[idx] = pow(aexpB,ccoeff::C1qq_2[idx]);
			aexpqg_2[idx] = pow(aexpB,ccoeff::C2qg_2[idx]/ccoeff::C1qg_2[idx]);
			aexpqqb_2[idx] = pow(aexpB,ccoeff::C3qqb_2[idx]/ccoeff::C2qqb_2[idx]);
			aexpqqp_2[idx] = pow(aexpB,ccoeff::C3qqp_2[idx]/ccoeff::C2qqp_2[idx]);
			aexpqqbp_2[idx] = pow(aexpB,ccoeff::C3qqbp_2[idx]/ccoeff::C2qqbp_2[idx]);
			}

		//At order n expand up to terms containing Cn (i.e. expand up to n for D, up to n-1 for OD, and up to n-2 for DOD)
		else if (opts.expc == 4)
			{
			for (int sign = mesq::positive; sign <= mesq::negative; sign++)
			for (int i = 0; i < mellinint::mdim; i++)
			{
				int idx = anomalous::index(i,sign);
				aexpqq_1[idx] = pow(aexpB,ccoeff::C1qq_1[idx]);
				aexpqq_2[idx] = pow(aexpB,ccoeff::C1qq_2[idx]);
			}
			}

		//Do not expand C at the denominator in the off diagonal and double off diagonal channels
		else if (opts.expc == 5)
			for (int sign = mesq::positive; sign <= mesq::negative; sign++)
			for (int i = 0; i < mellinint::mdim; i++)
			{
			int idx = anomalous::index(i,sign);
			aexpqq_1[idx] = pow(aexpB,ccoeff::C1qq_1[idx]);
			aexpqg_1[idx]   = (1. - xlambda/(1.+aass*ccoeff::C2qg_1[idx]/ccoeff::C1qg_1[idx]))/(1.-xlambda);
			aexpqqb_1[idx]  = (1. - xlambda/(1.+aass*ccoeff::C3qqb_1[idx] /ccoeff::C2qqb_1[idx] ))/(1.-xlambda);
			aexpqqp_1[idx]  = (1. - xlambda/(1.+aass*ccoeff::C3qqp_1[idx] /ccoeff::C2qqp_1[idx] ))/(1.-xlambda);
			aexpqqbp_1[idx] = (1. - xlambda/(1.+aass*ccoeff::C3qqbp_1[idx]/ccoeff::C2qqbp_1[idx]))/(1.-xlambda);
			aexpqq_2[idx] = pow(aexpB,ccoeff::C1qq_2[idx]);
			aexpqg_2[idx]   = (1. - xlambda/(1.+aass*ccoeff::C2qg_2[idx]/ccoeff::C1qg_2[idx]))/(1.-xlambda);
			aexpqqb_2[idx]  = (1. - xlambda/(1.+aass*ccoeff::C3qqb_2[idx] /ccoeff::C2qqb_2[idx] ))/(1.-xlambda);
			aexpqqp_2[idx]  = (1. - xlambda/(1.+aass*ccoeff::C3qqp_2[idx] /ccoeff::C2qqp_2[idx] ))/(1.-xlambda);
			aexpqqbp_2[idx] = (1. - xlambda/(1.+aass*ccoeff::C3qqbp_2[idx]/ccoeff::C2qqbp_2[idx]))/(1.-xlambda);
			}

		for (int sign = mesq::positive; sign <= mesq::negative; sign++)
			for (int i1 = 0; i1 < mellinint::mdim; i1++)
			for (int i2 = 0; i2 < mellinint::mdim; i2++)
			{
			int idx1 = anomalous::index(i1,mesq::positive);
			int idx2 = anomalous::index(i2,sign);
			int ii12 = hcoefficients::index(i1,i2,sign);
			qqb[ii12]  = aexpqq_1[idx1] * aexpqq_2[idx2];

			//qg_1[ii12] = aexp * aexpqg_1[idx1]*aexpqq_1[idx1]; //--> bug fix in dyres
			//qg_2[ii12] = aexp * aexpqg_2[idx2]*aexpqq_2[idx2]; //--> bug fix in dyres
			qg_1[ii12] = aexp * aexpqg_1[idx1]*aexpqq_2[idx2];
			qg_2[ii12] = aexp * aexpqg_2[idx2]*aexpqq_1[idx1];

			qq_1[ii12] = aexpqq_1[idx1] * aexpqqb_2[idx2] * aexp2; //*c1delta
			qq_2[ii12] = aexpqq_2[idx2] * aexpqqb_1[idx1] * aexp2; //*c1delta
			qqp_1[ii12] = aexpqq_1[idx1] * aexpqqbp_2[idx2] * aexp2; //*c1delta
			qqp_2[ii12] = aexpqq_2[idx2] * aexpqqbp_1[idx1] * aexp2; //*c1delta
			qqbp_1[ii12] = aexpqq_1[idx1] * aexpqqp_2[idx2] * aexp2; //*c1delta
			qqbp_2[ii12] = aexpqq_2[idx2] * aexpqqp_1[idx1] * aexp2; //*c1delta
			gg[ii12]   = aexp*aexpqg_1[idx1] * aexp*aexpqg_2[idx2];
			}

		//include aexp into aexpqg and aexpqq,qqp,qqb,qqbp
		for (int sign = mesq::positive; sign <= mesq::negative; sign++)
			for (int i = 0; i < mellinint::mdim; i++)
			{
			int idx = anomalous::index(i,sign);
			aexpqg_1[idx]   *= aexp;
			aexpqqb_1[idx]  *= aexp2;
			aexpqqp_1[idx]  *= aexp2;
			aexpqqbp_1[idx] *= aexp2;
			aexpqg_2[idx]   *= aexp;
			aexpqqb_2[idx]  *= aexp2;
			aexpqqp_2[idx]  *= aexp2;
			aexpqqbp_2[idx] *= aexp2;
			}

		}
		//cout << alogqq[0] << " expc::alogqq " << endl;
		//cout << alogqg[0] << " expc::alogqg " << endl;
		//cout << alogqqb[0] << " expc::alogqqb " << endl;
		//cout << endl;

		//cout << "cexp::qqb " << qqb[0] << endl;
		//cout << "cexp::qg " << qg[0]  << endl;
		//cout << "cexp::qq " << qq[0]  << endl;

		}


	// NNNLL
	//multiple parton case as in Eqs.(100-107) of https://arxiv.org/pdf/hep-ph/0508068.pdf
	//B3tilde additional pieces as in Eq. (50)
	//-2 beta1 C1 + 2 beta0 C1^2 - 4 beta0 C2

	//(B3qbar*(-2. + y)*y)/(2.*beta0*pow(1. - y,2))
	//1./beta0 * 1./2.  y*(y-2.)/pow(1.-y,2) * (-2 beta1 C1 + 2 beta0 C1^2 - 4 beta0 C2)
	//log(aexpC) * (C1^2 - beta1/beta0*C1 - 2*C2)

	//pieces from g4
	//(B3qbar*(-2. + y)*y)/(2.*beta0*pow(1. - y,2))
	//(B2qbar*beta1*((2. - y)*y + 2.*log1y))/(2.*pow(beta0,2)*pow(1. - y,2))
	//((rlogq2mur2-2.*rloga)*(3.*y*(-2.*B2qbar*pow(beta0,3)*(-2. + y)*(-1. + y)))/(6.*pow(beta0,3)*pow(-1. + y,3));

	//(1/2)*(1/beta0)*B3qbar       *y*(y-2)/(1-y)^2
	//-(1/2)*B2qbar*beta1/beta0^2  *y*(y-2)/(1-y)^2
	//B2qbar*beta1/beta0^2        *log1y  /(1-y)^2
	//-B2qbar((rlogq2mur2-2.*rloga)*y*(y-2)/(1-y)^2

	//(-beta1/beta0*C1+C1^2-2*C2)       *y*(y-2)/(1-y)^2
	//C1*beta1/beta0                      *y*(y-2)/(1-y)^2
	//-2*C1*beta1/beta0                    *log1y  /(1-y)^2
	//(2*beta0*C1)*((rlogq2mur2-2.*rloga) *y*(y-2)/(1-y)^2

	//Diagonal form factor
	//Gca_D =
	//  + (aS/Pi)*C1                           *(lam)/(1-lam)
	//
	//  + (aS^2/Pi^2)*(
	//		   + (C1^2/2 - C2)           *lam*(lam-2)/(1-lam)^2
	//		   - beta1/beta0*C1          *Log[1-lam]/(1-lam)^2
	//		   + beta0*C1*Log[Q2/muR2]   *lam*(lam-2)/(1-lam)^2
	//		   );

	//Off-diagonal form factor
	//Gca_OD =
	//  - Log[1 - lam]
	//
	//  + (aS/Pi)*(
	//	       + (C2/C1)                         *lam/(1-lam)
	//	       - (beta1/beta0                    *Log[1-lam])/(1-lam)
	//	       - (beta0*Log[Q2/muR2])            *lam/(1-lam)
	//	       )
	//
	//  + (aS^2/Pi^2)*(
	//		   - beta1^2/beta0^2                        *lam/(1-lam)^2
	//		   + beta2/beta0                            *lam/(1-lam)^2
	//		   - beta1^2/beta0^2                        *Log[1-lam]/(1-lam)^2
	//		   + 1/2*beta1^2/beta0^2                    *Log[1-lam]^2/(1-lam)^2
	//		   + (C2^2/C1^2/2 - C1*C3/C1^2)             *lam*(lam-2)/(1-lam)^2
	//		   - beta1/beta0*C2/C1                      *Log[1-lam]/(1-lam)^2
	//
	//		   - beta1                  *Log[Q2/muR2]   *lam/(1-lam)^2
	//		   + beta1                  *Log[Q2/muR2]   *Log[1-lam]/(1-lam)^2
	//		   + beta0*C2/C1            *Log[Q2/muR2]   *lam*(lam-2)/(1-lam)^2
	//		   );

	//Doubly off-diagonal form factor
	//Gca_DOD =
	//  - 2*Log[1 - lam]
	//  + (aS/Pi)* (
	//		+ (C3/C2)                      *lam/(1-lam)
	//		- (2*beta1/beta0               *Log[1-lam])/(1-lam)
	//		- (2*beta0*Log[Q2/muR2])       *lam/(1-lam)
	//		)
	//
	//  + (aS^2/Pi^2)* (
	//		    - 2*beta1^2/beta0^2                    *lam/(1-lam)^2
	//		    + 2*beta2/beta0                        *lam/(1-lam)^2
	//		    - 2*beta1^2/beta0^2                    *Log[1-lam]/(1-lam)^2
	//		    + beta1^2/beta0^2                      *Log[1-lam]^2/(1-lam)^2
	//		    + (C3^2/C2^2/2 - C4/C2)                *lam*(lam-2)/(1-lam)^2
	//		    + beta1/beta0*C3/C2                    *Log[1-lam]/(1-lam)^2
	//
	//		    - 2*beta1      *Log[Q2/muR2]          *lam/(1-lam)^2
	//		    + 2*beta1      *Log[Q2/muR2]          *Log[1-lam]/(1-lam)^2
	//		    + beta0*C3/C2  *Log[Q2/muR2]          *lam*(lam-2)/(1-lam)^2
	//		    );

	if (opts.order_expc == 3)
		{
		if (opts.mellin1d)
		{
		//Use pointers to allocate memory on the heap
		complex <double> *alogqq   = new complex <double> [mellinint::mdim*2];
		complex <double> *alogqg   = new complex <double> [mellinint::mdim*2];
		complex <double> *alogqqb  = new complex <double> [mellinint::mdim*2];
		complex <double> *alogqqp  = new complex <double> [mellinint::mdim*2];
		complex <double> *alogqqbp = new complex <double> [mellinint::mdim*2];
		//complex <double> *aexpqq   = new complex <double> [mellinint::mdim*2];
		//complex <double> *aexpqg   = new complex <double> [mellinint::mdim*2];
		//complex <double> *aexpqqb  = new complex <double> [mellinint::mdim*2];
		//complex <double> *aexpqqp  = new complex <double> [mellinint::mdim*2];
		//complex <double> *aexpqqbp = new complex <double> [mellinint::mdim*2];

		if (opts.expc == 0)
			if (opts.sumlogs)
			{
			fill(alogqq, alogqq+mellinint::mdim*2, 0.);
			fill(alogqg, alogqq+mellinint::mdim*2, 0.);
			fill(alogqqb, alogqqb+mellinint::mdim*2, 0.);
			fill(alogqqp, alogqqp+mellinint::mdim*2, 0.);
			fill(alogqqbp, alogqqbp+mellinint::mdim*2, 0.);
			}
			else
			{
			fill(aexpqq, aexpqq+mellinint::mdim*2, 1.);
			fill(aexpqg, aexpqg+mellinint::mdim*2, 1.);
			fill(aexpqqb, aexpqqb+mellinint::mdim*2, 1.);
			fill(aexpqqp, aexpqqp+mellinint::mdim*2, 1.);
			fill(aexpqqbp, aexpqqbp+mellinint::mdim*2, 1.);
			}
		else if (opts.expc == 1)
			{
			if (opts.sumlogs)
			{
			fill(alogqg, alogqq+mellinint::mdim*2, 0.);
			fill(alogqqb, alogqqb+mellinint::mdim*2, 0.);
			fill(alogqqp, alogqqp+mellinint::mdim*2, 0.);
			fill(alogqqbp, alogqqbp+mellinint::mdim*2, 0.);
			for (int sign = mesq::positive; sign <= mesq::negative; sign++)
				for (int i = 0; i < mellinint::mdim; i++)
				{
				int idx = anomalous::index(i,sign);
				alogqq[idx] = lamB*ccoeff::C1qq_delta
				+lamC*(0.5*pow(ccoeff::C1qq_delta,2) - ccoeff::C2qq_delta)
				-lamD*beta1/beta0*ccoeff::C1qq_delta
				+lamC*beta0*ccoeff::C1qq_delta*(resint::rlogq2mur2-2.*resint::rloga)
				;
				}
			}
			else
			{
			fill(aexpqg, aexpqg+mellinint::mdim*2, 1.);
			fill(aexpqqb, aexpqqb+mellinint::mdim*2, 1.);
			fill(aexpqqp, aexpqqp+mellinint::mdim*2, 1.);
			fill(aexpqqbp, aexpqqbp+mellinint::mdim*2, 1.);
			for (int sign = mesq::positive; sign <= mesq::negative; sign++)
				for (int i = 0; i < mellinint::mdim; i++)
				{
				int idx = anomalous::index(i,sign);
				aexpqq[idx] = pow(aexpB,ccoeff::C1qq_delta)
				*pow(aexpC,(0.5*pow(ccoeff::C1qq_delta,2) - ccoeff::C2qq_delta))
				*pow(aexpD,-beta1/beta0*ccoeff::C1qq_delta)
				*pow(aexpC,beta0*ccoeff::C1qq_delta*(resint::rlogq2mur2-2.*resint::rloga))
				;
				}
			}
			}
		else if (opts.expc == 2)
			{
			for (int sign = mesq::positive; sign <= mesq::negative; sign++)
			for (int i = 0; i < mellinint::mdim; i++)
			{
				int idx = anomalous::index(i,sign);
				complex <double> qqtayl = 1.;
				complex <double> qgtayl = 1.;
				complex <double> qqbtayl = 1.;
				complex <double> qqptayl = 1.;
				complex <double> qqbptayl = 1.;
				int fact = 1;
				for (int t = 1; t <= opts.ntaylor; t++)
				{
				fact *= t;
				qqtayl += 1./double(fact)* (pow(lamB*(ccoeff::C1qq[idx]-ccoeff::C1qq_delta),t)
								+pow(lamC*(0.5*pow((ccoeff::C1qq[idx]-ccoeff::C1qq_delta),2) - (ccoeff::C2qq[idx]-ccoeff::C2qq_delta)),t)
								+pow(-lamD*beta1/beta0*(ccoeff::C1qq[idx]-ccoeff::C1qq_delta),t)
								//+pow(lamC*beta0*(ccoeff::C1qq[idx]-ccoeff::C1qq_delta)*(resint::rlogq2mur2-2.*resint::rloga),t)
								);

				qgtayl += 1./double(fact)* (pow(lamB*ccoeff::C2qg[idx]/ccoeff::C1qg[idx],t)
								+pow(lamC*(pow(ccoeff::C2qg[idx]/ccoeff::C1qg[idx],2)/2. - ccoeff::C3qg[idx]/ccoeff::C1qg[idx]),t)
								+pow(-lamD*beta1/beta0*ccoeff::C2qg[idx]/ccoeff::C1qg[idx],t)
								//+pow(lamC*(beta0*ccoeff::C2qg[idx]/ccoeff::C1qg[idx])*(resint::rlogq2mur2-2.*resint::rloga),t)
								);

				qqbtayl += 1./double(fact) *(pow(lamB*ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx],t)
								+pow(lamC*(pow(ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx],2)/2.),t)
								-pow(lamD*beta1/beta0*ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx],t)
								//+pow(lamC*beta0*ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx]*(resint::rlogq2mur2-2.*resint::rloga),t)
								);

				qqptayl += 1./double(fact) *(pow(lamB*ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx],t)
								+pow(lamC*(pow(ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx],2)/2.),t)
								-pow(lamD*beta1/beta0*ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx],t)
								//+pow(lamC*beta0*ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx]*(resint::rlogq2mur2-2.*resint::rloga),t)
								);

				qqbptayl += 1./double(fact) *(pow(lamB*ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx],t)
								+pow(lamC*(pow(ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx],2)/2.),t)
								-pow(lamD*beta1/beta0*ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx],t)
								//+pow(lamC*beta0*ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx]*(resint::rlogq2mur2-2.*resint::rloga),t)
								);
				}
				if (opts.sumlogs)
				{
				alogqq[idx] = lamB*ccoeff::C1qq_delta
				+lamC*(0.5*pow(ccoeff::C1qq_delta,2) - ccoeff::C2qq_delta)
				-lamD*beta1/beta0*ccoeff::C1qq_delta
				+lamC*beta0*ccoeff::C1qq_delta*(resint::rlogq2mur2-2.*resint::rloga)
				+log(qqtayl);
				alogqg[idx] = log(qgtayl);
				alogqqb[idx] = log(qqbtayl);
				alogqqp[idx] = log(qqptayl);
				alogqqbp[idx] = log(qqbptayl);
				}
				else
				{
				aexpqq[idx] = pow(aexpB,ccoeff::C1qq_delta)
				*pow(aexpC,(0.5*pow(ccoeff::C1qq_delta,2) - ccoeff::C2qq_delta))
				*pow(aexpD,-beta1/beta0*ccoeff::C1qq_delta)
				*pow(aexpC,beta0*ccoeff::C1qq_delta*(resint::rlogq2mur2-2.*resint::rloga))
				*qqtayl;
				aexpqg[idx] = qgtayl;
				aexpqqb[idx] = qqbtayl;
				aexpqqp[idx] = qqptayl;
				aexpqqbp[idx] = qqbptayl;
				}
			}
			}
		else if (opts.expc == 3)
			for (int sign = mesq::positive; sign <= mesq::negative; sign++)
			for (int i = 0; i < mellinint::mdim; i++)
			{
			int idx = anomalous::index(i,sign);
			if (opts.sumlogs)
				{
	//		      alogqq[idx] = lamB*1./2.*C1QQ[idx]
	//			+lamC*1./4.*(0.5*pow(C1QQ[idx],2) - (C2NSqqM[idx] + C2SqqbM[idx]))
	//			-lamD*1./2.*beta1/beta0*C1QQ[idx]
	//			+lamC*1./2.*beta0*C1QQ[idx]*(resint::rlogq2mur2-2.*resint::rloga)
	//			;

				alogqq[idx] = lamB*ccoeff::C1qq[idx]
				+lamC*(0.5*pow(ccoeff::C1qq[idx],2) - ccoeff::C2qq[idx])
				-lamD*beta1/beta0*ccoeff::C1qq[idx]
				+lamC*beta0*ccoeff::C1qq[idx]*(LR-LQ)
				;

				alogqg[idx] = lamB*(ccoeff::C2qg[idx]/ccoeff::C1qg[idx])
				//+lamC*1./4.*(pow(C2qgM[idx]/C1QG[idx],2)/2. - 8.*ccoeff::C3qg[idx]/C1QG[idx])
				//-lamD*1./2.*beta1/beta0*C2qgM[idx]/C1QG[idx]
				//+lamC*1./2.*(beta0*C2qgM[idx]/C1QG[idx])*(resint::rlogq2mur2-2.*resint::rloga)
				+lamC*(pow(ccoeff::C2qg[idx]/ccoeff::C1qg[idx],2)/2. - ccoeff::C3qg[idx]/ccoeff::C1qg[idx])
				-lamD*beta1/beta0*ccoeff::C2qg[idx]/ccoeff::C1qg[idx]
				+lamC*(beta0*ccoeff::C2qg[idx]/ccoeff::C1qg[idx])*(LR-LQ)
				;

				alogqqb[idx] = lamB*ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx]
				+lamC*(pow(ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx],2)/2.)// - C4/C2)
				-lamD*beta1/beta0*ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx]
				+lamC*beta0*ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx]*(LR-LQ)
				;

				alogqqp[idx] = lamB*ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx]
				+lamC*(pow(ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx],2)/2.)// - C4/C2)
				-lamD*beta1/beta0*ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx]
				+lamC*beta0*ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx]*(LR-LQ)
				;

				alogqqbp[idx] = lamB*ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx]
				+lamC*(pow(ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx],2)/2.)// - C4/C2)
				-lamD*beta1/beta0*ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx]
				+lamC*beta0*ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx]*(LR-LQ)
				;

				}
			else
				{
				aexpqq[idx] = pow(aexpB,ccoeff::C1qq[idx])
				*pow(aexpC,(0.5*pow(ccoeff::C1qq[idx],2) - ccoeff::C2qq[idx]))
				*pow(aexpD,-beta1/beta0*ccoeff::C1qq[idx])
				*pow(aexpC,beta0*ccoeff::C1qq[idx]*(LR-LQ))
				;

				aexpqg[idx] = pow(aexpB,ccoeff::C2qg[idx]/ccoeff::C1qg[idx])
				*pow(aexpC,pow(ccoeff::C2qg[idx]/ccoeff::C1qg[idx],2)/2. - ccoeff::C3qg[idx]/ccoeff::C1qg[idx])
				*pow(aexpD,-beta1/beta0*ccoeff::C2qg[idx]/ccoeff::C1qg[idx])
				*pow(aexpC,beta0*ccoeff::C2qg[idx]/ccoeff::C1qg[idx]*(LR-LQ))
				;

				aexpqqb[idx] = pow(aexpB,ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx])
				*pow(aexpC,pow(ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx],2)/2.)// - C4/C2)
				*pow(aexpD,-beta1/beta0*ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx])
				*pow(aexpC,beta0*ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx]*(LR-LQ))
				;

				aexpqqp[idx] = pow(aexpB,ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx])
				*pow(aexpC,pow(ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx],2)/2.)// - C4/C2))
				*pow(aexpD,-beta1/beta0*ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx])
				*pow(aexpC,beta0*ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx]*(LR-LQ))
				;

				aexpqqbp[idx] = pow(aexpB,ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx])
				*pow(aexpC,pow(ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx],2)/2.)// - C4/C2))
				*pow(aexpD,-beta1/beta0*ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx])
				*pow(aexpC,beta0*ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx]*(LR-LQ))
				;

			//if (sign == mesq::positive)
			//	cout << mellinint::Np[i] << "  " << aexpqg[idx]
			//	     << " C2/C1 " << pow(ccoeff::C2qg[idx]/ccoeff::C1qg[idx],2)
			//	     << " C3/C1 " << ccoeff::C3qg[idx]/ccoeff::C2qg[idx]
			//	     << " C3/C2 qqb" << ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx]
			//	     << " C3/C2 qqp" << ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx]
			//	     << " C3/C2 qqbp" << ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx]
			//	     << endl;
			//else
			//	cout << mellinint::Nm[i] << "  " << aexpqg[idx]
			//	     << " C2/C1 " << pow(ccoeff::C2qg[idx]/ccoeff::C1qg[idx],2)
			//	     << " C3/C1 " << ccoeff::C3qg[idx]/ccoeff::C2qg[idx]
			//	     << " C3/C2 qqb" << ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx]
			//	     << " C3/C2 qqp" << ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx]
			//	     << " C3/C2 qqbp" << ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx]
			//	     << endl;

				}
			}

		//cout << endl;
		//cout << alogqq[0] << " expc::alogqq " << endl;
		//cout << alogqg[0] << " expc::alogqg " << endl;
		//cout << alogqqb[0] << " expc::alogqqb " << endl;
		//cout << alogqqp[0] << " expc::alogqqb " << endl;
		//cout << alogqqbp[0] << " expc::alogqqbp " << endl;
		//cout << endl;

		//{
		//  cout << endl;
		//  int i;
		//  i = 0;
		//  //i = opts.mellinrule - 1;
		//  cout << aexpqq[i]   << " expc::aexpqq " << endl;
		//  cout << aexpqg[i]   << " expc::aexpqg " << endl;
		//  cout << aexpqqb[i]  << " expc::aexpqqb " << endl;
		//  cout << aexpqqp[i]  << " expc::aexpqqb " << endl;
		//  cout << aexpqqbp[i] << " expc::aexpqqbp " << endl;
		//  cout << endl;
		//}
		else if (opts.expc == 4)
			{
			if (opts.sumlogs)
			{
			fill(alogqqb, alogqqb+mellinint::mdim*2, 0.);
			fill(alogqqp, alogqqp+mellinint::mdim*2, 0.);
			fill(alogqqbp, alogqqbp+mellinint::mdim*2, 0.);
			}
			else
			{
			fill(aexpqqb, aexpqqb+mellinint::mdim*2, 1.);
			fill(aexpqqp, aexpqqp+mellinint::mdim*2, 1.);
			fill(aexpqqbp, aexpqqbp+mellinint::mdim*2, 1.);
			}
			for (int sign = mesq::positive; sign <= mesq::negative; sign++)
			for (int i = 0; i < mellinint::mdim; i++)
			{
				int idx = anomalous::index(i,sign);
				if (opts.sumlogs)
				{
				alogqq[idx] = lamB*ccoeff::C1qq[idx]
				+lamC*(0.5*pow(ccoeff::C1qq[idx],2) - ccoeff::C2qq[idx])
				-lamD*beta1/beta0*ccoeff::C1qq[idx]
				+lamC*beta0*ccoeff::C1qq[idx]*(LR-LQ)
				;
				alogqg[idx]   = log(1.-xlambda/(1.+aass*ccoeff::C2qg[idx]/ccoeff::C1qg[idx])) - log(1.-xlambda);
				}
				else
				{
				aexpqq[idx] = pow(aexpB,ccoeff::C1qq[idx])
				*pow(aexpC,(0.5*pow(ccoeff::C1qq[idx],2) - ccoeff::C2qq[idx]))
				*pow(aexpD,-beta1/beta0*ccoeff::C1qq[idx])
				*pow(aexpC,beta0*ccoeff::C1qq[idx]*(LR-LQ))
				;
				aexpqg[idx]   = (1.-xlambda/(1.+aass*ccoeff::C2qg[idx]/ccoeff::C1qg[idx]))/(1.-xlambda);
				}
			}
			}
		//for (int i = 0; i < mellinint::mdim; i++)
		//{
		//  int i;
		//  i = 0;
		//  //i = opts.mellinrule - 1;
		//  cout << i << endl;
		//  cout << aexpqq[i]   << " gint::alogqq " << endl;
		//  cout << aexpqg[i]   << "  " << exp(gint::alogqg[i])  << " gint::alogqg " << endl;
		//  cout << aexpqqb[i]  << " gint::alogqqb " << endl;
		//  cout << aexpqqp[i]  << " gint::alogqqb " << endl;
		//  cout << aexpqqbp[i] << " gint::alogqqbp " << endl;
		//  cout << endl;
		//}
		else if (opts.expc == 5)
			{
			for (int sign = mesq::positive; sign <= mesq::negative; sign++)
			for (int i = 0; i < mellinint::mdim; i++)
			{
			int idx = anomalous::index(i,sign);
			if (opts.sumlogs)
				{
				if (!opts.numexpc) gint::calc(b);
				alogqq[idx] = lamB*ccoeff::C1qq[idx]
				+lamC*(0.5*pow(ccoeff::C1qq[idx],2) - ccoeff::C2qq[idx])
				-lamD*beta1/beta0*ccoeff::C1qq[idx]
				+lamC*beta0*ccoeff::C1qq[idx]*(LR-LQ)
				;

				complex <double> C1,C2,C3,rootC,n2llalog;
				C1 = ccoeff::C1qg[idx];
				C2 = ccoeff::C2qg[idx];
				C3 = ccoeff::C3qg[idx];
				rootC = sqrt(pow(C2/C1,2)-4.*C3/C1);
				complex <double> lam = xlambda;
				//lamB = aass*lam/(1.-lam);
				//lamC = aass2*lam*(lam-2.)/pow(1.-lam,2);
				//lamD = aass2*log1lam/pow(1.-lam,2);

				//original formulas
				complex <double> oldN2LL,oldN3LL,oldN3LLb,oldN3LLmur;
				oldN2LL    = lamB*(C2/C1);
				oldN3LL    = (pow(C2/C1,2)/2. - C3/C1) * lamC - beta1/beta0 * C2/C1 * lamC/2.;
				oldN3LLb   = - beta1/beta0 * C2/C1*lamD + beta1/beta0 * C2/C1 * lamC/2.;
				oldN3LLmur = beta0*C2/C1*lamC*(LR-LQ);

				//Formulas without expanding the denominator
				complex <double> newN2LL,newN3LL,newN3LLb,newN3LLmur;

				newN2LL = log(1.-lam/(1.+aass*C2/C1)) - log(1.-lam);
				//newN2LL = log((1.-lam+aass*C2/C1)/(1.-lam)/(1.+aass*C2/C1)); //--> Check

				/*
				newN3LL = 0.
				//These pieces are problematic because of as0 at the denominator -> does not converge for bprescription = 2
				//integral of 3*C3/C1*beta0*as3_3l/(1+C2/C1*as0^2)
				//			-3.*C3/C1/(1.+aass2*C3/C1) * lamC/2.
				//integral of 2*C2/C1*beta1*as3_3l/(1+C2/C1*as0) --> not used
				//-2.*beta1/beta0*C2/C1/(1.+aass*C2/C1) * lamC/2.

				//integral of  2.*C2/C1*beta1*as3_3l/(1+C2/C1*as1_1l)
				-2.*beta1/beta0*(C1/C2*log((1.+aass*C2/C1/(1.-lam))/(1.+aass*C2/C1))-lamB)

				//This piece is problematic because of as0 at the denominator -> does not converge for bprescription = 2
				//integral of 	3*C3/C1*beta0*as3_3l/(1+C2/C1*as0^2)+2*C2/C1*beta1*as3_3l/(1+C2/C1*as0)
				-(aass*C2/C1*C3/C1*(3.+2.*aass*beta1/beta0)+2.*beta1/beta0*C2/C1+3.*C3/C1) / (1.+aass*C2/C1)/(1.+aass2*C3/C1) * lamC/2.

				//integral of beta0*as1_1l/(1+C2/C1*as1_1l+C3/C1*as2_2l)
				+C2/C1/sqrt(-pow(C2/C1,2)+4.*C3/C1)*(+atan((aass*C2/C1 + 2.*(1.-lam))/(aass*sqrt(-pow(C2/C1,2) + 4.*C3/C1)))
									-atan((aass*C2/C1 + 2.)         /(aass*sqrt(-pow(C2/C1,2) + 4.*C3/C1))))
				+0.5*log((1.+aass*C2/C1+aass2*C3/C1)/(pow(lam-1.,2)+aass*C2/C1*(1.-lam)+aass2*C3/C1))

				//integral of beta1*as2_2l/(1+C2/C1*as1_1l)
				+beta1/beta0*C1/C2*log((1.+aass*C2/C1/(1.-lam))/(1.+aass*C2/C1))

				//integral of +2*C2/C1*beta0*as2_2l/(1+C2/C1*as1_1l)
				+2.*log((1.+aass*C2/C1/(1.-lam))/(1.+aass*C2/C1))

				//sum of the two pieces above
				//+(2.+beta1/beta0*C1/C2)*log((1.-lam+aass*C2/C1))/(1.-lam)/(1.+aass*C2/C1)) //--> Check

				//integral of -beta0*as1_1l
				+log(1.-lam)

				//integral of -beta1*as2_2l
				-beta1/beta0*lamB
				;
				*/

				//simplified alternative:
				newN3LL =

				//integral of 2*C3/C1*beta0*as3_3l/(1+C3/C1*as2_2l) --> this piece is not correct, should be 2*C3/C1*beta0*as3_3l/(1+C2/C1*as1_1l+C3/C1*as2_2l)
				//-log((1.+aass2*C3/C1)/(pow(lam-1.,2)+aass2*C3/C1))
				//-2.*log(1.-lam)
				// --> equivalent to:
				//+log((pow(1.-lam,2)+aass2*C3/C1)/(pow(1.-lam,2)*(1.+aass2*C3/C1)))

				//integral of 2*C3/C1*beta0*as3_3l /(1+C2/C1*as1_1l+C3/C1*as2_2l)
				//C2/C1/sqrt(pow(C2/C1,2)-4.*C3/C1)
				//*log((aass*(C2/C1+sqrt(pow(C2/C1,2)-4.*C3/C1))+2.)/(aass*(C2/C1+sqrt(pow(C2/C1,2)-4.*C3/C1))+2.*(1.-lam))
				//     *(aass*(-C2/C1+sqrt(pow(C2/C1,2)-4.*C3/C1))-2.*(1.-lam))/(aass*(-C2/C1+sqrt(pow(C2/C1,2)-4.*C3/C1))-2.))
				//+log((pow(1.-lam,2)+aass*C2/C1*(1.-lam) + aass2*C3/C1)/pow(1.-lam,2)/(1.+aass*C2/C1+aass2*C3/C1))

				+log(pow((aass*(rootC+C2/C1)+2.)*(aass*(rootC-C2/C1)-2.*(1.-lam))/(aass*(rootC+C2/C1)+2.*(1.-lam))/(aass*(rootC-C2/C1)-2.),C2/C1/rootC)
					*(pow(1.-lam,2)+aass*C2/C1*(1.-lam)+aass2*C3/C1)/pow(1.-lam,2)/(1.+aass*C2/C1+aass2*C3/C1))
					//*(1.+aass*C2/C1/(1.-lam)+aass2*C3/C1/pow(1.-lam,2))/(1.+aass*C2/C1+aass2*C3/C1))


				//integral of  C2/C1*beta1*as3_3l/(1+C2/C1*as1_1l)
				//-beta1/beta0*(C1/C2*log((1.+aass*C2/C1/(1.-lam))/(1.+aass*C2/C1))-lamB)
				;

			//newN3LLb = 1./(beta0*C2*(lam-1.))
			//	*beta1*(+2.*aass*C2*lam
			//		+log(1.-lam)*(2.*aass*C2-C1*(lam-1.)*(log(1.-lam)-2.*log(1.-(C1*(lam-1.))/(aass*C2))))
			//		+2.*C1*(lam-1.)*(-Li2(-((C1)/(aass*C2))) + Li2((C1*(lam-1.))/(aass*C2))));


				newN3LLb =

				+2.*aass*beta1/beta0*(lam+log(1.-lam))/(lam-1.)
				-C1/C2*beta1/beta0 * log(1.-lam) *(log(1.-lam)+2.*log( aass*C2/C1/(1.-lam + aass*C2/C1)))
				+2.*C1/C2*beta1/beta0*(-Li2(-C1/C2/aass) + Li2(C1/C2*(lam-1.)/aass));


				newN3LLmur =
				+2.*beta0*(C1/C2*log((1.+aass*C2/C1/(1.-lam))/(1.+aass*C2/C1))-lamB)*(LR-LQ);


				//alternative with as0 at the denominator as much as possible --> Not good, because asLO/as diverges at lam -> 1
				//newN2LL = lamB*C2/C1/(1.+aass*C2/C1);
				//newN3LL = 0.5*
				//	(-(aass*C2/C1*C3/C1*(3.+aass*beta1/beta0)+beta1/beta0*C2/C1+3.*C3/C1) / (1.+aass*C2/C1)/(1.+aass2*C3/C1) * lamC
				//
				//	 +2.*C2/C1/sqrt(-pow(C2/C1,2)+4.*C3/C1)*(+atan((aass*C2/C1 + 2.*(1.-lam))/(aass*sqrt(-pow(C2/C1,2) + 4.*C3/C1)))
				//						 -atan((aass*C2/C1 + 2.)         /(aass*sqrt(-pow(C2/C1,2) + 4.*C3/C1))))
				//	 +2.*log(1.-lam)
				//	 +4.*log((1.+aass*C2/C1/(1.-lam))/(1.+aass*C2/C1))
				//	 +log((1.+aass*C2/C1+aass2*C3/C1)/(pow(lam-1.,2)+aass*C2/C1*(1.-lam)+aass2*C3/C1)));
				//newN3LLb   = beta1/beta0*C2/C1*(lamC/2.-lamD)/(1.+aass*C2/C1);
				//newN3LLmur = beta0*C2/C1*lamC*(LR-LQ)/(1.+aass*C2/C1);


				/*** test formulas ***

				complex <double> lam = 0.001;

				complex <double> lB = aass*lam/(1.-lam);
				complex <double> lC = aass2*lam*(lam-2.)/pow(1.-lam,2);
				complex <double> lD = aass2*log(1.-lam)/pow(1.-lam,2);

				complex <double> C1 = 1.;
				complex <double> C2 = complex <double> (0.1, 0.1);
				complex <double> C3 = complex <double> (0.2, 0.2);

				complex <double> oldN3LLb = - beta1/beta0 * C2/C1*lD + beta1/beta0 * C2/C1 * lC/2.;
				complex <double> newN3LLb = 1./(beta0*C2*(lam-1.))
				*beta1*(2.*aass*C2*lam+log(1.-lam)*(2.*aass*C2-C1*(lam-1.)
											*(log(1.-lam)-2.*log(1.-(C1*(lam-1.))/(aass*C2))))
					+2.*C1*(lam-1.)*(- Li2(-((C1)/(aass*C2))) + Li2((C1*(lam-1.) )/(aass*C2))));
				*/

				/**********************************************************/
				//alogqg[idx] = newN2LL+newN3LL+newN3LLb+newN3LLmur;

				/*
				if (sign == mesq::positive)
				cout << mellinint::Np[i];
				else
				cout << mellinint::Nm[i];
				cout << " old " << oldN3LL//+oldN3LL+oldN3LLb+oldN3LLmur
				<< " new " << newN3LL//+newN3LL+newN3LLb+newN3LLmur
				<< " gint " << gint::alogqg[i] << endl;
				*/

	//#define rootC(C1,C2,C3) sqrt(pow(C2/C1,2)-4.*C3/C1)
	//#define logfinite(C1,C2) log((1.-lam+aass*C2/C1)/(1.-lam)/(1.+aass*C2/C1))
	//#define n3llexpc(C1,C2,C3) logfinite(C1,C2)+C2/C1/rootC(C1,C2,C3)*log((aass*(C2/C1+rootC(C1,C2,C3))+2.)/(aass*(C2/C1+rootC(C1,C2,C3))+2.*(1.-lam))*(aass*(-C2/C1+rootC(C1,C2,C3))-2.*(1.-lam))/(aass*(-C2/C1+rootC(C1,C2,C3))-2.))+log((pow(1.-lam,2)+aass*C2/C1*(1.-lam)+aass2*C3/C1)/pow(1.-lam,2)/(1.+aass*C2/C1+aass2*C3/C1))-beta1/beta0*(C1/C2*logfinite(C1,C2)-lamB)+1./(beta0*C2*(lam-1.))*beta1*(+2.*aass*C2*lam+log(1.-lam)*(2.*aass*C2-C1*(lam-1.)*(log(1.-lam)-2.*log(1.-(C1*(lam-1.))/(aass*C2))))+2.*C1*(lam-1.)*(-Li2(-((C1)/(aass*C2))) + Li2((C1*(lam-1.))/(aass*C2))))+2.*beta0*(C1/C2*logfinite(C1,C2)-lamB)*(LR-LQ);
	//		      alogqg[idx]   = n3llexpc(ccoeff::C1qg[idx],   ccoeff::C2qg[idx],   ccoeff::C3qg[idx]);
	//		      alogqqb[idx]  = n3llexpc(ccoeff::C2qqb[idx],  ccoeff::C3qqb[idx],  0.);
	//		      alogqqp[idx]  = n3llexpc(ccoeff::C2qqp[idx],  ccoeff::C3qqp[idx],  0.);
	//		      alogqqbp[idx] = n3llexpc(ccoeff::C2qqbp[idx], ccoeff::C3qqbp[idx], 0.);

	#define cdeflog(c1,c2,c3) C1 = c1; C2 = c2; C3 = c3; rootC = sqrt(pow(C2/C1,2)-4.*C3/C1); n2llalog = log((1.-lam+aass*C2/C1)/(1.-lam)/(1.+aass*C2/C1));
	#define n3llalog() n2llalog\
				+C2/C1/rootC*log((aass*(C2/C1+rootC)+2.)/(aass*(C2/C1+rootC)+2.*(1.-lam))*(aass*(-C2/C1+rootC)-2.*(1.-lam))/(aass*(-C2/C1+rootC)-2.)) \
				+log((pow(1.-lam,2)+aass*C2/C1*(1.-lam)+aass2*C3/C1)/pow(1.-lam,2)/(1.+aass*C2/C1+aass2*C3/C1)) \
				-beta1/beta0*(C1/C2*n2llalog-lamB)		\
				+1./(beta0*C2*(lam-1.))*beta1*(+2.*aass*C2*lam+log(1.-lam)*(2.*aass*C2-C1*(lam-1.)*(log(1.-lam)-2.*log(1.-(C1*(lam-1.))/(aass*C2))))+2.*C1*(lam-1.)*(-Li2(-((C1)/(aass*C2))) + Li2((C1*(lam-1.))/(aass*C2)))) \
				+2.*beta0*(C1/C2*n2llalog-lamB)*(LR-LQ);

				cdeflog(ccoeff::C1qg[idx],   ccoeff::C2qg[idx],   ccoeff::C3qg[idx]); alogqg[idx]   = n3llalog();
				//cdeflog(ccoeff::C1qg[idx],   ccoeff::C2qg[idx],   0.); alogqg[idx]   = n3llalog();
				cdeflog(ccoeff::C2qqb[idx],  ccoeff::C3qqb[idx],  0.               ); alogqqb[idx]  = n3llalog();
				cdeflog(ccoeff::C2qqp[idx],  ccoeff::C3qqp[idx],  0.               ); alogqqp[idx]  = n3llalog();
				cdeflog(ccoeff::C2qqbp[idx], ccoeff::C3qqbp[idx], 0.               ); alogqqbp[idx] = n3llalog();

	//		      C1 = ccoeff::C1qg[idx]; C2 = ccoeff::C2qg[idx]; C3 = ccoeff::C3qg[idx]; rootC = sqrt(pow(C2/C1,2)-4.*C3/C1); n2llalog = log((1.-lam+aass*C2/C1)/(1.-lam)/(1.+aass*C2/C1));
	//		      alogqg[idx]   = n2llalog+C2/C1/rootC*log((aass*(C2/C1+rootC)+2.)/(aass*(C2/C1+rootC)+2.*(1.-lam))*(aass*(-C2/C1+rootC)-2.*(1.-lam))/(aass*(-C2/C1+rootC)-2.))+log((pow(1.-lam,2)+aass*C2/C1*(1.-lam)+aass2*C3/C1)/pow(1.-lam,2)/(1.+aass*C2/C1+aass2*C3/C1))-beta1/beta0*(C1/C2*n2llalog-lamB)+1./(beta0*C2*(lam-1.))*beta1*(+2.*aass*C2*lam+log(1.-lam)*(2.*aass*C2-C1*(lam-1.)*(log(1.-lam)-2.*log(1.-(C1*(lam-1.))/(aass*C2))))+2.*C1*(lam-1.)*(-Li2(-((C1)/(aass*C2))) + Li2((C1*(lam-1.))/(aass*C2))))+2.*beta0*(C1/C2*n2llalog-lamB)*(LR-LQ);
	//		      C1 = ccoeff::C2qqb[idx];  C2 = ccoeff::C3qqb[idx];  C3 = 0.; rootC = sqrt(pow(C2/C1,2)-4.*C3/C1); n2llalog = log((1.-lam+aass*C2/C1)/(1.-lam)/(1.+aass*C2/C1));
	//		      alogqqb[idx]  = n2llalog+C2/C1/rootC*log((aass*(C2/C1+rootC)+2.)/(aass*(C2/C1+rootC)+2.*(1.-lam))*(aass*(-C2/C1+rootC)-2.*(1.-lam))/(aass*(-C2/C1+rootC)-2.))+log((pow(1.-lam,2)+aass*C2/C1*(1.-lam)+aass2*C3/C1)/pow(1.-lam,2)/(1.+aass*C2/C1+aass2*C3/C1))-beta1/beta0*(C1/C2*n2llalog-lamB)+1./(beta0*C2*(lam-1.))*beta1*(+2.*aass*C2*lam+log(1.-lam)*(2.*aass*C2-C1*(lam-1.)*(log(1.-lam)-2.*log(1.-(C1*(lam-1.))/(aass*C2))))+2.*C1*(lam-1.)*(-Li2(-((C1)/(aass*C2))) + Li2((C1*(lam-1.))/(aass*C2))))+2.*beta0*(C1/C2*n2llalog-lamB)*(LR-LQ);
	//		      C1 = ccoeff::C2qqp[idx];  C2 = ccoeff::C3qqp[idx];  C3 = 0.; rootC = sqrt(pow(C2/C1,2)-4.*C3/C1); n2llalog = log((1.-lam+aass*C2/C1)/(1.-lam)/(1.+aass*C2/C1));
	//		      alogqqp[idx]  = n2llalog+C2/C1/rootC*log((aass*(C2/C1+rootC)+2.)/(aass*(C2/C1+rootC)+2.*(1.-lam))*(aass*(-C2/C1+rootC)-2.*(1.-lam))/(aass*(-C2/C1+rootC)-2.))+log((pow(1.-lam,2)+aass*C2/C1*(1.-lam)+aass2*C3/C1)/pow(1.-lam,2)/(1.+aass*C2/C1+aass2*C3/C1))-beta1/beta0*(C1/C2*n2llalog-lamB)+1./(beta0*C2*(lam-1.))*beta1*(+2.*aass*C2*lam+log(1.-lam)*(2.*aass*C2-C1*(lam-1.)*(log(1.-lam)-2.*log(1.-(C1*(lam-1.))/(aass*C2))))+2.*C1*(lam-1.)*(-Li2(-((C1)/(aass*C2))) + Li2((C1*(lam-1.))/(aass*C2))))+2.*beta0*(C1/C2*n2llalog-lamB)*(LR-LQ);
	//		      C1 = ccoeff::C2qqbp[idx]; C2 = ccoeff::C3qqbp[idx]; C3 = 0.; rootC = sqrt(pow(C2/C1,2)-4.*C3/C1); n2llalog = log((1.-lam+aass*C2/C1)/(1.-lam)/(1.+aass*C2/C1));
	//		      alogqqbp[idx] = n2llalog+C2/C1/rootC*log((aass*(C2/C1+rootC)+2.)/(aass*(C2/C1+rootC)+2.*(1.-lam))*(aass*(-C2/C1+rootC)-2.*(1.-lam))/(aass*(-C2/C1+rootC)-2.))+log((pow(1.-lam,2)+aass*C2/C1*(1.-lam)+aass2*C3/C1)/pow(1.-lam,2)/(1.+aass*C2/C1+aass2*C3/C1))-beta1/beta0*(C1/C2*n2llalog-lamB)+1./(beta0*C2*(lam-1.))*beta1*(+2.*aass*C2*lam+log(1.-lam)*(2.*aass*C2-C1*(lam-1.)*(log(1.-lam)-2.*log(1.-(C1*(lam-1.))/(aass*C2))))+2.*C1*(lam-1.)*(-Li2(-((C1)/(aass*C2))) + Li2((C1*(lam-1.))/(aass*C2))))+2.*beta0*(C1/C2*n2llalog-lamB)*(LR-LQ);

				//alogqqb[idx]  = log(1.-xlambda/(1.+aass*ccoeff::C3qqb[idx] /ccoeff::C2qqb[idx] )) - log(1.-xlambda);
				//alogqqp[idx]  = log(1.-xlambda/(1.+aass*ccoeff::C3qqp[idx] /ccoeff::C2qqp[idx] )) - log(1.-xlambda);
				//alogqqbp[idx] = log(1.-xlambda/(1.+aass*ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx])) - log(1.-xlambda);
				}
			else
				{
	//		      if (opts.numexpc) gint::calc(b);
	//		      complex <double> C1,C2,C3,rootC,n2llaexp;
	//		      C1 = ccoeff::C1qg[idx];
	//		      C2 = ccoeff::C2qg[idx];
	//		      C3 = ccoeff::C3qg[idx];
	//		      rootC = sqrt(pow(C2/C1,2)-4.*C3/C1);
	//		      complex <double> lam = xlambda;
	//
	//		      //original formulas
	//		      complex <double> oldN2LL,oldN3LL,oldN3LLb,oldN3LLmur;
	//		      oldN2LL    = exp(lamB*(C2/C1));
	//		      oldN3LL    = exp((pow(C2/C1,2)/2. - C3/C1) * lamC - beta1/beta0 * C2/C1 * lamC/2.);
	//		      oldN3LLb   = exp(- beta1/beta0 * C2/C1*lamD + beta1/beta0 * C2/C1 * lamC/2.);
	//		      oldN3LLmur = exp(beta0*C2/C1*lamC*(LR-LQ));
	//
	//		      //Formulas without expanding the denominator
	//		      complex <double> newN2LL,newN3LL,newN3LLb,newN3LLmur;
	//
	//#define cdefexp(c1,c2,c3) C1 = c1; C2 = c2; C3 = c3; rootC = sqrt(pow(C2/C1,2)-4.*C3/C1); n2llaexp = (1.+aass*C2/C1/(1.-lam))/(1.+aass*C2/C1); //n2llaexp = (1.-lam+aass*C2/C1)/(1.-lam)/(1.+aass*C2/C1);
	//		      cdefexp(ccoeff::C1qg[idx],   ccoeff::C2qg[idx],   ccoeff::C3qg[idx]);
	//		      newN2LL = (1.+aass*C2/C1/(1.-lam))/(1.+aass*C2/C1);
	//		      newN3LL = exp(log((aass*(rootC+C2/C1)+2.)*(aass*(rootC-C2/C1)-2.*(1.-lam))/(aass*(rootC+C2/C1)+2.*(1.-lam))/(aass*(rootC-C2/C1)-2.))*C2/C1/rootC) \
	//			*(pow(1.-lam,2)+aass*C2/C1*(1.-lam)+aass2*C3/C1)/pow(1.-lam,2)/(1.+aass*C2/C1+aass2*C3/C1) \
	//			* exp(-log(n2llaexp)*beta1/beta0*C1/C2) *pow(aexpB,beta1/beta0);
	//		      newN3LLb = exp(+2.*aass*beta1/beta0*(lam+log(1.-lam))/(lam-1.)-C1/C2*beta1/beta0 * log(1.-lam) *(log(1.-lam)+2.*log( aass*C2/C1/(1.-lam + aass*C2/C1)))+2.*C1/C2*beta1/beta0*(-Li2(-C1/C2/aass)+Li2(C1/C2*(lam-1.)/aass)));
	//		      newN3LLmur = exp(log(n2llaexp)*2.*beta0*(C1/C2)*(LR-LQ))	\
	//			*pow(aexpB,-2.*beta0*(LR-LQ));
	//
	//		      if (sign == mesq::positive)
	//			cout << mellinint::Np[i];
	//		      else
	//			cout << mellinint::Nm[i];
	//		      cout << " old " << oldN3LLb//+oldN3LL+oldN3LLb+oldN3LLmur
	//			   << " new " << newN3LLb//+newN3LL+newN3LLb+newN3LLmur
	//			   << " gint " << exp(gint::alogqg[i]) << endl;


				aexpqq[idx] = pow(aexpB,ccoeff::C1qq[idx])
				*pow(aexpC,(0.5*pow(ccoeff::C1qq[idx],2) - ccoeff::C2qq[idx]))
				*pow(aexpD,-beta1/beta0*ccoeff::C1qq[idx])
				*pow(aexpC,beta0*ccoeff::C1qq[idx]*(LR-LQ))
				;

				if (isnan_ofast(real(aexpqq[idx])) || isnan_ofast(imag(aexpqq[idx])))
				{
				cout << "Warning, aexpqq[" << idx << "] = " << aexpqq[idx] << ", b = "  << b << ", lambda = " << xlambda
					<< " aexpB " << aexpB << " aexpC " << aexpC << " aexpD " << aexpD
					<< endl;
				}

				complex <double> C1,C2,C3,rootC,n2llaexp;
				complex <double> lam = xlambda;

				//cdefexp(ccoeff::C1qg[idx],   ccoeff::C2qg[idx],   ccoeff::C3qg[idx]); aexpqg[idx]   = n3llaexp();
				cdefexp(ccoeff::C1qg[idx],   ccoeff::C2qg[idx],   0.               ); aexpqg[idx]   = n3llaexp();
				cdefexp(ccoeff::C2qqb[idx],  ccoeff::C3qqb[idx],  0.               ); aexpqqb[idx]  = n3llaexp();
				cdefexp(ccoeff::C2qqp[idx],  ccoeff::C3qqp[idx],  0.               ); aexpqqp[idx]  = n3llaexp();
				cdefexp(ccoeff::C2qqbp[idx], ccoeff::C3qqbp[idx], 0.               ); aexpqqbp[idx] = n3llaexp();

				//aexpqg[idx]   = (1.-xlambda/(1.+aass*ccoeff::C2qg[idx]/ccoeff::C1qg[idx]))/(1.-xlambda);
				//aexpqqb[idx]  = (1.-xlambda/(1.+aass*ccoeff::C3qqb[idx] /ccoeff::C2qqb[idx] ))/(1.-xlambda);
				//aexpqqp[idx]  = (1.-xlambda/(1.+aass*ccoeff::C3qqp[idx] /ccoeff::C2qqp[idx] ))/(1.-xlambda);
				//aexpqqbp[idx] = (1.-xlambda/(1.+aass*ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx]))/(1.-xlambda);
				}
			}
			}
		/*
		//for (int i = 0; i < mellinint::mdim; i++)
		//{
		if (opts.numexpc) gint::calc(b);
		int i;
		i = 0;
		//  //i = opts.mellinrule - 1;
		//cout << i << endl;
		//  cout << aexpqq[i]   << " gint::alogqq " << endl;
		cout << alogqg[i]   << "  " << gint::alogqg[i]  << " gint::alogqg " << endl;
		//  cout << aexpqqb[i]  << " gint::alogqqb " << endl;
		//  cout << aexpqqp[i]  << " gint::alogqqb " << endl;
		//  cout << aexpqqbp[i] << " gint::alogqqbp " << endl;
		cout << endl;
		//}
		*/
		//cout << endl;
		//cout << " gint::logasl "   << gint::logasl_expc      << " expc " << logasl      << endl;
		//cout << " gint::alogqq "   << gint::alogqq[0]   << " expc " << log(aexpqq[0])   << endl;
		//cout << " gint::alogqg "   << gint::alogqg[0]   << " expc " << log(aexpqg[0])   << endl;
		//cout << " gint::alogqqb "  << gint::alogqqb[0]  << " expc " << log(aexpqqb[0])  << endl;
		//cout << " gint::alogqqp "  << gint::alogqqp[0]  << " expc " << log(aexpqqp[0])  << endl;
		//cout << " gint::alogqqbp " << gint::alogqqbp[0] << " expc " << log(aexpqqbp[0]) << endl;

		//!!! gint::alog are evaluated with blim::sudakov, not with blim:expc !!!
		if (opts.numexpc)
			for (int sign = mesq::positive; sign <= mesq::negative; sign++)
			for (int i = 0; i < mellinint::mdim; i++)
			{
			int idx = anomalous::index(i,sign);

			if (opts.sumlogs)
				{
				alogqq[idx]   = gint::alogqq[idx];
				alogqg[idx]   = gint::alogqg[idx];
				alogqqb[idx]  = gint::alogqqb[idx];
				alogqqp[idx]  = gint::alogqqp[idx];
				alogqqbp[idx] = gint::alogqqbp[idx];
				}
			else
				{
				aexpqq[idx]   = exp(gint::alogqq[idx]);
				aexpqg[idx]   = exp(gint::alogqg[idx]);
				aexpqqb[idx]  = exp(gint::alogqqb[idx]);
				aexpqqp[idx]  = exp(gint::alogqqp[idx]);
				aexpqqbp[idx] = exp(gint::alogqqbp[idx]);
				}
			}


		for (int sign = mesq::positive; sign <= mesq::negative; sign++)
			for (int i = 0; i < mellinint::mdim; i++)
			{
			int idx = anomalous::index(i,sign);
			int ii = hcoefficients::index(i,sign);
			if (opts.sumlogs)
			{
				qqb[ii] = exp(alogqq[idx]+alogqq[idx]+sudakov::logS);
				qg[ii]  = exp(-logasl+alogqg[idx]+alogqq[idx]+sudakov::logS);
				qq[ii]  = exp(alogqq[idx]+alogqqb[idx]-2.*logasl+sudakov::logS);
				qqp[ii] = exp(alogqq[idx]+alogqqbp[idx]-2.*logasl+sudakov::logS);
				qqbp[ii] = exp(alogqq[idx]+alogqqp[idx]-2.*logasl+sudakov::logS);
				gg[ii]   = exp(-logasl+alogqg[idx]-logasl+alogqg[idx]+sudakov::logS);
				qbg[ii]  = exp(-logasl+alogqg[idx]-2.*logasl+alogqqb[idx]+sudakov::logS);
				qpg[ii]  = exp(-logasl+alogqg[idx]-2.*logasl+alogqqp[idx]+sudakov::logS);
				qbpg[ii] = exp(-logasl+alogqg[idx]-2.*logasl+alogqqbp[idx]+sudakov::logS);
			}
			else
			{
				qqb[ii] = aexpqq[idx] * aexpqq[idx];
				qg[ii]  = aexp * aexpqg[idx]*aexpqq[idx];
				qq[ii]  = aexpqq[idx] * aexpqqb[idx] * aexp2; //*c1delta
				qqp[ii] = aexpqq[idx] * aexpqqbp[idx] * aexp2; //*c1delta
				qqbp[ii] = aexpqq[idx] * aexpqqp[idx] * aexp2; //*c1delta
				gg[ii]  = aexp*aexpqg[idx] * aexp*aexpqg[idx];
				qbg[ii] = aexp * aexpqg[idx] *aexp2*aexpqqb[idx];
				qpg[ii] = aexp * aexpqg[idx] *aexp2*aexpqqp[idx];
				qbpg[ii] = aexp * aexpqg[idx] *aexp2*aexpqqbp[idx];

				//		  cout << endl;
				//		  cout << bstar << "  " << i << "  " << aexp << "  " << aexpqq[idx] << "  " << aexpqg[idx] << endl;
				//		  cout << "  " << pow(aexpB,1./2.*(C2qgM[idx]/C1QG[idx]))
				//		       << "  " << pow(aexpC,1./4.*(pow(C2qgM[idx]/C1QG[idx],2)/2. - 8.*ccoeff::C3qg[idx]/C1QG[idx]))
				//		       << "  " << pow(aexpC,(pow(ccoeff::C2qg[idx]/ccoeff::C1qg[idx],2)/2. - ccoeff::C3qg[idx]/ccoeff::C1qg[idx]))
				//		       << "  " << pow(aexpD,-1./2.*beta1/beta0*C2qgM[idx]/C1QG[idx])
				//		       << "  " << pow(aexpC,1./2.*(beta0*C2qgM[idx]/C1QG[idx])*(resint::rlogq2mur2-2.*resint::rloga))
				//		       << endl;
			}
			if (isnan_ofast(real(qg[ii])) || isnan_ofast(imag(qg[ii])))
			{
				int idx = anomalous::index(i,sign);
				cout << "Warning, expc::qg[" << i << "] = " << qg[ii] << ", b = "  << b << ", lambda = " << xlambda << " logasl " << logasl
				<< " aexp " << aexp << " aexpqq " << aexpqq[idx] << " aexpqg " << aexpqg[idx]
				<< endl;
			}

			}
		delete[] alogqq;
		delete[] alogqg;
		delete[] alogqqb;
		delete[] alogqqp;
		delete[] alogqqbp;
		//include aexp into aexpqg and aexpqq,qqp,qqb,qqbp
		for (int sign = mesq::positive; sign <= mesq::negative; sign++)
			for (int i = 0; i < mellinint::mdim; i++)
			{
			int idx = anomalous::index(i,sign);
			aexpqg[idx]   *= aexp;
			aexpqqb[idx]  *= aexp2;
			aexpqqp[idx]  *= aexp2;
			aexpqqbp[idx] *= aexp2;
			}
		//delete[] aexpqq;
		//delete[] aexpqg;
		//delete[] aexpqqb;
		//delete[] aexpqqp;
		//delete[] aexpqqbp;

		//cout << "mellin1d" << endl;
		//for (int sign = mesq::positive; sign <= mesq::negative; sign++)
		//  for (int i = 0; i < mellinint::mdim; i++)
		//    {
		//	int idx = anomalous::index(i,sign);
		//	cout << i << " " << sign << " aexpqq   " << aexpqq[idx] << endl;
		//	cout << i << " " << sign << " aexpqg   " << aexpqg[idx] << endl;
		//	//cout << i << " " << sign << " aexpqqb  " << aexpqqb[idx] << endl;
		//	//cout << i << " " << sign << " aexpqqp  " << aexpqqp[idx] << endl;
		//	//cout << i << " " << sign << " aexpqqbp " << aexpqqbp[idx] << endl;
		//    }

		}
		else //mellin2d
		{
		//complex <double> aexpqq_1[mellinint::mdim*2];
		//complex <double> aexpqg_1[mellinint::mdim*2];
		//complex <double> aexpqqb_1[mellinint::mdim*2];
		//complex <double> aexpqqp_1[mellinint::mdim*2];
		//complex <double> aexpqqbp_1[mellinint::mdim*2];
		//complex <double> aexpqq_2[mellinint::mdim*2];
		//complex <double> aexpqg_2[mellinint::mdim*2];
		//complex <double> aexpqqb_2[mellinint::mdim*2];
		//complex <double> aexpqqp_2[mellinint::mdim*2];
		//complex <double> aexpqqbp_2[mellinint::mdim*2];

		//case opts.expc == 0
		fill(aexpqq_1, aexpqq_1+mellinint::mdim*2, 1.);
		fill(aexpqg_1, aexpqg_1+mellinint::mdim*2, 1.);
		fill(aexpqqqed_1, aexpqqqed_1+mellinint::mdim*2, 1.);
		fill(aexpqgamma_1, aexpqgamma_1+mellinint::mdim*2, 1.);
		fill(aexpqqb_1, aexpqqb_1+mellinint::mdim*2, 1.);
		fill(aexpqqp_1, aexpqqp_1+mellinint::mdim*2, 1.);
		fill(aexpqqbp_1, aexpqqbp_1+mellinint::mdim*2, 1.);
		fill(aexpqq_2, aexpqq_2+mellinint::mdim*2, 1.);
		fill(aexpqg_2, aexpqg_2+mellinint::mdim*2, 1.);
		fill(aexpqqqed_2, aexpqqqed_2+mellinint::mdim*2, 1.);
		fill(aexpqgamma_2, aexpqgamma_2+mellinint::mdim*2, 1.);
		fill(aexpqqb_2, aexpqqb_2+mellinint::mdim*2, 1.);
		fill(aexpqqp_2, aexpqqp_2+mellinint::mdim*2, 1.);
		fill(aexpqqbp_2, aexpqqbp_2+mellinint::mdim*2, 1.);


		if (opts.expc == 1)
			{
			for (int sign = mesq::positive; sign <= mesq::negative; sign++)
			for (int i = 0; i < mellinint::mdim; i++)
			{
				int idx = anomalous::index(i,sign);
				aexpqq_1[idx] = pow(aexpB,ccoeff::C1qq_delta)
				*pow(aexpC,(0.5*pow(ccoeff::C1qq_delta,2) - ccoeff::C2qq_delta))
				*pow(aexpD,-beta1/beta0*ccoeff::C1qq_delta)
				*pow(aexpC,beta0*ccoeff::C1qq_delta*(resint::rlogq2mur2-2.*resint::rloga)) //*(LR-LQ)
				;
				aexpqq_2[idx] = pow(aexpB,ccoeff::C1qq_delta)
				*pow(aexpC,(0.5*pow(ccoeff::C1qq_delta,2) - ccoeff::C2qq_delta))
				*pow(aexpD,-beta1/beta0*ccoeff::C1qq_delta)
				*pow(aexpC,beta0*ccoeff::C1qq_delta*(resint::rlogq2mur2-2.*resint::rloga)) //*(LR-LQ)
				;
			}
			}
		/*
		else if (opts.expc == 2)
			{
			for (int sign = mesq::positive; sign <= mesq::negative; sign++)
			for (int i = 0; i < mellinint::mdim; i++)
			{
				int idx = anomalous::index(i,sign);
				complex <double> qqtayl = 1.;
				complex <double> qgtayl = 1.;
				complex <double> qqbtayl = 1.;
				complex <double> qqptayl = 1.;
				complex <double> qqbptayl = 1.;
				int fact = 1;
				for (int t = 1; t <= opts.ntaylor; t++)
				{
				fact *= t;
				qqtayl += 1./double(fact)* (pow(lamB*(ccoeff::C1qq[idx]-ccoeff::C1qq_delta),t)
								+pow(lamC*(0.5*pow((ccoeff::C1qq[idx]-ccoeff::C1qq_delta),2) - (ccoeff::C2qq[idx]-ccoeff::C2qq_delta)),t)
								+pow(-lamD*beta1/beta0*(ccoeff::C1qq[idx]-ccoeff::C1qq_delta),t)
								//+pow(lamC*beta0*(ccoeff::C1qq[idx]-ccoeff::C1qq_delta)*(resint::rlogq2mur2-2.*resint::rloga),t)
								);

				qgtayl += 1./double(fact)* (pow(lamB*ccoeff::C2qg[idx]/ccoeff::C1qg[idx],t)
								+pow(lamC*(pow(ccoeff::C2qg[idx]/ccoeff::C1qg[idx],2)/2. - ccoeff::C3qg[idx]/ccoeff::C1qg[idx]),t)
								+pow(-lamD*beta1/beta0*ccoeff::C2qg[idx]/ccoeff::C1qg[idx],t)
								//+pow(lamC*(beta0*ccoeff::C2qg[idx]/ccoeff::C1qg[idx])*(resint::rlogq2mur2-2.*resint::rloga),t)
								);

				qqbtayl += 1./double(fact) *(pow(lamB*ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx],t)
								+pow(lamC*(pow(ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx],2)/2.),t)
								-pow(lamD*beta1/beta0*ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx],t)
								//+pow(lamC*beta0*ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx]*(resint::rlogq2mur2-2.*resint::rloga),t)
								);

				qqptayl += 1./double(fact) *(pow(lamB*ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx],t)
								+pow(lamC*(pow(ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx],2)/2.),t)
								-pow(lamD*beta1/beta0*ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx],t)
								//+pow(lamC*beta0*ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx]*(resint::rlogq2mur2-2.*resint::rloga),t)
								);

				qqbptayl += 1./double(fact) *(pow(lamB*ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx],t)
								+pow(lamC*(pow(ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx],2)/2.),t)
								-pow(lamD*beta1/beta0*ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx],t)
								//+pow(lamC*beta0*ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx]*(resint::rlogq2mur2-2.*resint::rloga),t)
								);
				}
				aexpqq[idx] = pow(aexpB,ccoeff::C1qq_delta)
				*pow(aexpC,(0.5*pow(ccoeff::C1qq_delta,2) - ccoeff::C2qq_delta))
				*pow(aexpD,-beta1/beta0*ccoeff::C1qq_delta)
				*pow(aexpC,beta0*ccoeff::C1qq_delta*(resint::rlogq2mur2-2.*resint::rloga))
				*qqtayl;
				aexpqg[idx] = qgtayl;
				aexpqqb[idx] = qqbtayl;
				aexpqqp[idx] = qqptayl;
				aexpqqbp[idx] = qqbptayl;
			}
			}
		*/
		else if (opts.expc == 3)
			for (int sign = mesq::positive; sign <= mesq::negative; sign++)
			for (int i = 0; i < mellinint::mdim; i++)
			{
			int idx = anomalous::index(i,sign);

			aexpqq_1[idx] = pow(aexpB,ccoeff::C1qq_1[idx])
				*pow(aexpC,(0.5*pow(ccoeff::C1qq_1[idx],2) - ccoeff::C2qq_1[idx]))
				*pow(aexpD,-beta1/beta0*ccoeff::C1qq_1[idx])
				*pow(aexpC,beta0*ccoeff::C1qq_1[idx]*(LR-LQ))
				;
			aexpqg_1[idx] = pow(aexpB,ccoeff::C2qg_1[idx]/ccoeff::C1qg_1[idx])
				*pow(aexpC,pow(ccoeff::C2qg_1[idx]/ccoeff::C1qg_1[idx],2)/2. - ccoeff::C3qg_1[idx]/ccoeff::C1qg_1[idx])
				*pow(aexpD,-beta1/beta0*ccoeff::C2qg_1[idx]/ccoeff::C1qg_1[idx])
				*pow(aexpC,beta0*ccoeff::C2qg_1[idx]/ccoeff::C1qg_1[idx]*(LR-LQ))
				;
			aexpqqb_1[idx] = pow(aexpB,ccoeff::C3qqb_1[idx]/ccoeff::C2qqb_1[idx])
				*pow(aexpC,pow(ccoeff::C3qqb_1[idx]/ccoeff::C2qqb_1[idx],2)/2.)// - C4/C2)
				*pow(aexpD,-beta1/beta0*ccoeff::C3qqb_1[idx]/ccoeff::C2qqb_1[idx])
				*pow(aexpC,beta0*ccoeff::C3qqb_1[idx]/ccoeff::C2qqb_1[idx]*(LR-LQ))
				;
			aexpqqp_1[idx] = pow(aexpB,ccoeff::C3qqp_1[idx]/ccoeff::C2qqp_1[idx])
				*pow(aexpC,pow(ccoeff::C3qqp_1[idx]/ccoeff::C2qqp_1[idx],2)/2.)// - C4/C2))
				*pow(aexpD,-beta1/beta0*ccoeff::C3qqp_1[idx]/ccoeff::C2qqp_1[idx])
				*pow(aexpC,beta0*ccoeff::C3qqp_1[idx]/ccoeff::C2qqp_1[idx]*(LR-LQ))
				;
			aexpqqbp_1[idx] = pow(aexpB,ccoeff::C3qqbp_1[idx]/ccoeff::C2qqbp_1[idx])
				*pow(aexpC,pow(ccoeff::C3qqbp_1[idx]/ccoeff::C2qqbp_1[idx],2)/2.)// - C4/C2))
				*pow(aexpD,-beta1/beta0*ccoeff::C3qqbp_1[idx]/ccoeff::C2qqbp_1[idx])
				*pow(aexpC,beta0*ccoeff::C3qqbp_1[idx]/ccoeff::C2qqbp_1[idx]*(LR-LQ))
				;
			aexpqq_2[idx] = pow(aexpB,ccoeff::C1qq_2[idx])
				*pow(aexpC,(0.5*pow(ccoeff::C1qq_2[idx],2) - ccoeff::C2qq_2[idx]))
				*pow(aexpD,-beta1/beta0*ccoeff::C1qq_2[idx])
				*pow(aexpC,beta0*ccoeff::C1qq_2[idx]*(LR-LQ))
				;
			aexpqg_2[idx] = pow(aexpB,ccoeff::C2qg_2[idx]/ccoeff::C1qg_2[idx])
				*pow(aexpC,pow(ccoeff::C2qg_2[idx]/ccoeff::C1qg_2[idx],2)/2. - ccoeff::C3qg_2[idx]/ccoeff::C1qg_2[idx])
				*pow(aexpD,-beta1/beta0*ccoeff::C2qg_2[idx]/ccoeff::C1qg_2[idx])
				*pow(aexpC,beta0*ccoeff::C2qg_2[idx]/ccoeff::C1qg_2[idx]*(LR-LQ))
				;
			aexpqqb_2[idx] = pow(aexpB,ccoeff::C3qqb_2[idx]/ccoeff::C2qqb_2[idx])
				*pow(aexpC,pow(ccoeff::C3qqb_2[idx]/ccoeff::C2qqb_2[idx],2)/2.)// - C4/C2)
				*pow(aexpD,-beta1/beta0*ccoeff::C3qqb_2[idx]/ccoeff::C2qqb_2[idx])
				*pow(aexpC,beta0*ccoeff::C3qqb_2[idx]/ccoeff::C2qqb_2[idx]*(LR-LQ))
				;
			aexpqqp_2[idx] = pow(aexpB,ccoeff::C3qqp_2[idx]/ccoeff::C2qqp_2[idx])
				*pow(aexpC,pow(ccoeff::C3qqp_2[idx]/ccoeff::C2qqp_2[idx],2)/2.)// - C4/C2))
				*pow(aexpD,-beta1/beta0*ccoeff::C3qqp_2[idx]/ccoeff::C2qqp_2[idx])
				*pow(aexpC,beta0*ccoeff::C3qqp_2[idx]/ccoeff::C2qqp_2[idx]*(LR-LQ))
				;
			aexpqqbp_2[idx] = pow(aexpB,ccoeff::C3qqbp_2[idx]/ccoeff::C2qqbp_2[idx])
				*pow(aexpC,pow(ccoeff::C3qqbp_2[idx]/ccoeff::C2qqbp_2[idx],2)/2.)// - C4/C2))
				*pow(aexpD,-beta1/beta0*ccoeff::C3qqbp_2[idx]/ccoeff::C2qqbp_2[idx])
				*pow(aexpC,beta0*ccoeff::C3qqbp_2[idx]/ccoeff::C2qqbp_2[idx]*(LR-LQ))
				;

			//cout << mellinint::Np_1[i] << "  " << aexpqg_1[idx] << "  " << aexpqg_2[idx]
			//     << " C1 " << fabs(ccoeff::C1qg_1[idx])
			//     << " C2 " << fabs(ccoeff::C2qg_1[idx])
			//     << " C3 " << fabs(ccoeff::C3qg_1[idx])
			//     << endl;
			}
		else if (opts.expc == 4)
			{
			fill(aexpqqb_1, aexpqqb_1+mellinint::mdim*2, 1.);
			fill(aexpqqp_1, aexpqqp_1+mellinint::mdim*2, 1.);
			fill(aexpqqbp_1, aexpqqbp_1+mellinint::mdim*2, 1.);
			fill(aexpqqb_2, aexpqqb_2+mellinint::mdim*2, 1.);
			fill(aexpqqp_2, aexpqqp_2+mellinint::mdim*2, 1.);
			fill(aexpqqbp_2, aexpqqbp_2+mellinint::mdim*2, 1.);
			for (int sign = mesq::positive; sign <= mesq::negative; sign++)
			for (int i = 0; i < mellinint::mdim; i++)
			{
				int idx = anomalous::index(i,sign);
				aexpqq_1[idx] = pow(aexpB,ccoeff::C1qq_1[idx])
				*pow(aexpC,(0.5*pow(ccoeff::C1qq_1[idx],2) - ccoeff::C2qq_1[idx]))
				*pow(aexpD,-beta1/beta0*ccoeff::C1qq_1[idx])
				*pow(aexpC,beta0*ccoeff::C1qq_1[idx]*(LR-LQ))
				;
				aexpqg_1[idx]   = (1.-xlambda/(1.+aass*ccoeff::C2qg_1[idx]/ccoeff::C1qg_1[idx]))/(1.-xlambda);
				aexpqq_2[idx] = pow(aexpB,ccoeff::C1qq_2[idx])
				*pow(aexpC,(0.5*pow(ccoeff::C1qq_2[idx],2) - ccoeff::C2qq_2[idx]))
				*pow(aexpD,-beta1/beta0*ccoeff::C1qq_2[idx])
				*pow(aexpC,beta0*ccoeff::C1qq_2[idx]*(LR-LQ))
				;
				aexpqg_2[idx]   = (1.-xlambda/(1.+aass*ccoeff::C2qg_2[idx]/ccoeff::C1qg_2[idx]))/(1.-xlambda);
			}
			}
		else if (opts.expc == 5)
			for (int sign = mesq::positive; sign <= mesq::negative; sign++)
			for (int i = 0; i < mellinint::mdim; i++)
			{
			int idx = anomalous::index(i,sign);
			aexpqq_1[idx] = pow(aexpB,ccoeff::C1qq_1[idx])
				*pow(aexpC,(0.5*pow(ccoeff::C1qq_1[idx],2) - ccoeff::C2qq_1[idx]))
				*pow(aexpD,-beta1/beta0*ccoeff::C1qq_1[idx])
				*pow(aexpC,beta0*ccoeff::C1qq_1[idx]*(LR-LQ))
				;
			aexpqq_2[idx] = pow(aexpB,ccoeff::C1qq_2[idx])
				*pow(aexpC,(0.5*pow(ccoeff::C1qq_2[idx],2) - ccoeff::C2qq_2[idx]))
				*pow(aexpD,-beta1/beta0*ccoeff::C1qq_2[idx])
				*pow(aexpC,beta0*ccoeff::C1qq_2[idx]*(LR-LQ))
				;

			complex <double> C1,C2,C3,rootC,n2llaexp;
			complex <double> lam = xlambda;

			//cdefexp(ccoeff::C1qg_1[idx],   ccoeff::C2qg_1[idx],   ccoeff::C3qg_1[idx]); aexpqg_1[idx]   = n3llaexp();
			cdefexp(ccoeff::C1qg_1[idx],   ccoeff::C2qg_1[idx],   0.               ); aexpqg_1[idx]   = n3llaexp();
			cdefexp(ccoeff::C2qqb_1[idx],  ccoeff::C3qqb_1[idx],  0.               ); aexpqqb_1[idx]  = n3llaexp();
			cdefexp(ccoeff::C2qqp_1[idx],  ccoeff::C3qqp_1[idx],  0.               ); aexpqqp_1[idx]  = n3llaexp();
			cdefexp(ccoeff::C2qqbp_1[idx], ccoeff::C3qqbp_1[idx], 0.               ); aexpqqbp_1[idx] = n3llaexp();

			//cdefexp(ccoeff::C1qg_2[idx],   ccoeff::C2qg_2[idx],   ccoeff::C3qg_2[idx]); aexpqg_2[idx]   = n3llaexp();
			cdefexp(ccoeff::C1qg_2[idx],   ccoeff::C2qg_2[idx],   0.               ); aexpqg_2[idx]   = n3llaexp();
			cdefexp(ccoeff::C2qqb_2[idx],  ccoeff::C3qqb_2[idx],  0.               ); aexpqqb_2[idx]  = n3llaexp();
			cdefexp(ccoeff::C2qqp_2[idx],  ccoeff::C3qqp_2[idx],  0.               ); aexpqqp_2[idx]  = n3llaexp();
			cdefexp(ccoeff::C2qqbp_2[idx], ccoeff::C3qqbp_2[idx], 0.               ); aexpqqbp_2[idx] = n3llaexp();
			}


		for (int sign = mesq::positive; sign <= mesq::negative; sign++)
			for (int i1 = 0; i1 < mellinint::mdim; i1++)
			for (int i2 = 0; i2 < mellinint::mdim; i2++)
			{
			int idx1 = anomalous::index(i1,mesq::positive);
			int idx2 = anomalous::index(i2,sign);
			int ii12 = hcoefficients::index(i1,i2,sign);
			qqb[ii12]  = aexpqq_1[idx1] * aexpqq_2[idx2];
			qg_1[ii12] = aexp * aexpqg_1[idx1]*aexpqq_2[idx2];
			qg_2[ii12] = aexp * aexpqg_2[idx2]*aexpqq_1[idx1];
			qq_1[ii12] = aexpqq_1[idx1] * aexpqqb_2[idx2] * aexp2; //*c1delta
			qq_2[ii12] = aexpqq_2[idx2] * aexpqqb_1[idx1] * aexp2; //*c1delta
			qqp_1[ii12] = aexpqq_1[idx1] * aexpqqbp_2[idx2] * aexp2; //*c1delta
			qqp_2[ii12] = aexpqq_2[idx2] * aexpqqbp_1[idx1] * aexp2; //*c1delta
			qqbp_1[ii12] = aexpqq_1[idx1] * aexpqqp_2[idx2] * aexp2; //*c1delta
			qqbp_2[ii12] = aexpqq_2[idx2] * aexpqqp_1[idx1] * aexp2; //*c1delta
			gg[ii12]   = aexp*aexpqg_1[idx1] * aexp*aexpqg_2[idx2];
			qbg_1[ii12] = aexp * aexpqg_1[idx1] *aexp2*aexpqqb_2[idx2];
			qbg_2[ii12] = aexp * aexpqg_2[idx2] *aexp2*aexpqqb_1[idx1];
			qpg_1[ii12] = aexp * aexpqg_1[idx1] *aexp2*aexpqqp_2[idx2];
			qpg_2[ii12] = aexp * aexpqg_2[idx2] *aexp2*aexpqqp_1[idx1];
			qbpg_1[ii12] = aexp * aexpqg_1[idx1] *aexp2*aexpqqbp_2[idx2];
			qbpg_2[ii12] = aexp * aexpqg_2[idx2] *aexp2*aexpqqbp_1[idx1];
			}
		//include aexp into aexpqg and aexpqq,qqp,qqb,qqbp
		for (int sign = mesq::positive; sign <= mesq::negative; sign++)
			for (int i = 0; i < mellinint::mdim; i++)
			{
			int idx = anomalous::index(i,sign);
			aexpqg_1[idx]   *= aexp;
			aexpqqb_1[idx]  *= aexp2;
			aexpqqp_1[idx]  *= aexp2;
			aexpqqbp_1[idx] *= aexp2;
			aexpqg_2[idx]   *= aexp;
			aexpqqb_2[idx]  *= aexp2;
			aexpqqp_2[idx]  *= aexp2;
			aexpqqbp_2[idx] *= aexp2;
			}
		//cout << "mellin2d" << endl;
		//for (int sign = mesq::positive; sign <= mesq::negative; sign++)
		//  for (int i = 0; i < mellinint::mdim; i++)
		//    {
		//	int idx = anomalous::index(i,sign);
		//	cout << i << " " << sign << " aexpqq   " << aexpqq_1[idx] << "  " << aexpqq_1[idx] - aexpqq_2[idx] << endl;
		//	cout << i << " " << sign << " aexpqg   " << aexpqg_1[idx] << "  " << aexpqg_1[idx] - aexpqg_1[idx] << endl;
		//	//cout << i << " " << sign << " aexpqqb_1  " << aexpqqb_1[idx] << endl;
		//	//cout << i << " " << sign << " aexpqqp_1  " << aexpqqp_1[idx] << endl;
		//	//cout << i << " " << sign << " aexpqqbp_1 " << aexpqqbp_1[idx] << endl;
		//	//cout << i << " " << sign << " aexpqqb_2  " << aexpqqb_2[idx] << endl;
		//	//cout << i << " " << sign << " aexpqqp_2  " << aexpqqp_2[idx] << endl;
		//	//cout << i << " " << sign << " aexpqqbp_2 " << aexpqqbp_2[idx] << endl;
		//    }


		}
		//cout << alogqq[0] << " expc::alogqq " << endl;
		//cout << alogqg[0] << " expc::alogqg " << endl;
		//cout << alogqqb[0] << " expc::alogqqb " << endl;
		//cout << alogqqp[0] << " expc::alogqqb " << endl;
		//cout << alogqqbp[0] << " expc::alogqqbp " << endl;
		//cout << endl;
		}


	if (opts.order_expc == 4)
		{
		if (opts.mellin1d)
		{
		if (opts.expc == 0)
			{
			fill(aexpqq, aexpqq+mellinint::mdim*2, 1.);
			fill(aexpqg, aexpqg+mellinint::mdim*2, 1.);
			fill(aexpqqqed, aexpqqqed+mellinint::mdim*2, 1.);
			fill(aexpqgamma, aexpqgamma+mellinint::mdim*2, 1.);
			fill(aexpqqb, aexpqqb+mellinint::mdim*2, 1.);
			fill(aexpqqp, aexpqqp+mellinint::mdim*2, 1.);
			fill(aexpqqbp, aexpqqbp+mellinint::mdim*2, 1.);
			}
		else if (opts.expc == 1)
			{
			fill(aexpqg, aexpqg+mellinint::mdim*2, 1.);
			fill(aexpqqb, aexpqqb+mellinint::mdim*2, 1.);
			fill(aexpqqp, aexpqqp+mellinint::mdim*2, 1.);
			fill(aexpqqbp, aexpqqbp+mellinint::mdim*2, 1.);
			for (int sign = mesq::positive; sign <= mesq::negative; sign++)
			for (int i = 0; i < mellinint::mdim; i++)
			{
				int idx = anomalous::index(i,sign);
				aexpqq[idx] = pow(aexpB,ccoeff::C1qq_delta)
				*pow(aexpC,(0.5*pow(ccoeff::C1qq_delta,2) - ccoeff::C2qq_delta))
				*pow(aexpD,-beta1/beta0*ccoeff::C1qq_delta)
				*pow(aexpC,beta0*ccoeff::C1qq_delta*LQR)
				*pow(aexpE,pow(beta1/beta0,2)*ccoeff::C1qq_delta)
				*pow(aexpE,-beta2/beta0*ccoeff::C1qq_delta)
				*pow(aexpF,-1./3.*(pow(ccoeff::C1qq_delta,3)-3.*ccoeff::C1qq_delta*ccoeff::C2qq_delta+3.*ccoeff::C3qq_delta))
				*pow(aexpG,beta1/pow(beta0,2)*(beta1*ccoeff::C1qq_delta - beta0*pow(ccoeff::C1qq_delta,2) + 2.*beta0*ccoeff::C2qq_delta))
				*pow(aexpH,-pow(beta1/beta0,2)*ccoeff::C1qq_delta)
				*pow(aexpG,2.*beta1*ccoeff::C1qq_delta*LQR)
				*pow(aexpF,(-beta1*ccoeff::C1qq_delta+beta0*pow(ccoeff::C1qq_delta,2)-2.*beta0*ccoeff::C2qq_delta)*LQR)
				*pow(aexpF,pow(beta0,2)*ccoeff::C1qq_delta*pow(LQR,2))
				;
			}
			}
		else if (opts.expc == 3)
			for (int sign = mesq::positive; sign <= mesq::negative; sign++)
			for (int i = 0; i < mellinint::mdim; i++)
			{
			int idx = anomalous::index(i,sign);
			aexpqq[idx] = pow(aexpB,ccoeff::C1qq[idx])
				*pow(aexpC,(0.5*pow(ccoeff::C1qq[idx],2) - ccoeff::C2qq[idx]))
				*pow(aexpD,-beta1/beta0*ccoeff::C1qq[idx])
				*pow(aexpC,beta0*ccoeff::C1qq[idx]*(LR-LQ))
				;

			aexpqg[idx] = pow(aexpB,ccoeff::C2qg[idx]/ccoeff::C1qg[idx])
				*pow(aexpC,pow(ccoeff::C2qg[idx]/ccoeff::C1qg[idx],2)/2. - ccoeff::C3qg[idx]/ccoeff::C1qg[idx])
				*pow(aexpD,-beta1/beta0*ccoeff::C2qg[idx]/ccoeff::C1qg[idx])
				*pow(aexpC,beta0*ccoeff::C2qg[idx]/ccoeff::C1qg[idx]*(LR-LQ))
				;

			aexpqqb[idx] = pow(aexpB,ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx])
				*pow(aexpC,pow(ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx],2)/2.)// - C4/C2)
				*pow(aexpD,-beta1/beta0*ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx])
				*pow(aexpC,beta0*ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx]*(LR-LQ))
				;

			aexpqqp[idx] = pow(aexpB,ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx])
				*pow(aexpC,pow(ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx],2)/2.)// - C4/C2))
				*pow(aexpD,-beta1/beta0*ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx])
				*pow(aexpC,beta0*ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx]*(LR-LQ))
				;

			aexpqqbp[idx] = pow(aexpB,ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx])
				*pow(aexpC,pow(ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx],2)/2.)// - C4/C2))
				*pow(aexpD,-beta1/beta0*ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx])
				*pow(aexpC,beta0*ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx]*(LR-LQ))
				;
			}
		else if (opts.expc == 4)
			{
			fill(aexpqqb, aexpqqb+mellinint::mdim*2, 1.);
			fill(aexpqqp, aexpqqp+mellinint::mdim*2, 1.);
			fill(aexpqqbp, aexpqqbp+mellinint::mdim*2, 1.);
			for (int sign = mesq::positive; sign <= mesq::negative; sign++)
			for (int i = 0; i < mellinint::mdim; i++)
			{
				int idx = anomalous::index(i,sign);

				//aexpqq[idx] = pow(aexpB,ccoeff::C1qq[idx])
				//  *pow(aexpC,(0.5*pow(ccoeff::C1qq[idx],2) - ccoeff::C2qq[idx]))
				//  *pow(aexpD,-beta1/beta0*ccoeff::C1qq[idx])
				//  *pow(aexpC,beta0*ccoeff::C1qq[idx]*(LQR))
				//  *pow(aexpE,-pow(beta1/beta0,2)*ccoeff::C1qq[idx])
				//  *pow(aexpE,beta2/beta0*ccoeff::C1qq[idx])
				//  *pow(aexpF,1./3.*(pow(ccoeff::C1qq[idx],3)-3.*ccoeff::C1qq[idx]*ccoeff::C2qq[idx]+3.*ccoeff::C3qq[idx]))
				//  *pow(aexpG,-beta1/pow(beta0,2)*(beta1*ccoeff::C1qq[idx] - beta0*pow(ccoeff::C1qq[idx],2) + 2.*beta0*ccoeff::C2qq[idx]))
				//  *pow(aexpH,pow(beta1/beta0,2)*ccoeff::C1qq[idx])
				//  *pow(aexpG,2.*beta1*ccoeff::C1qq[idx]*LQR)
				//  *pow(aexpF,(-beta1*ccoeff::C1qq[idx]+beta0*pow(ccoeff::C1qq[idx],2)-2.*beta0*ccoeff::C2qq[idx])*LQR)
				//  *pow(aexpF,pow(beta0,2)*ccoeff::C1qq[idx]*pow(LQR,2))
				//  ;

				aexpqq[idx] = exp(
						+lamB*(ccoeff::C1qq[idx])
						+lamC*((0.5*pow(ccoeff::C1qq[idx],2) - ccoeff::C2qq[idx]))
						+lamD*(-beta1/beta0*ccoeff::C1qq[idx])
						+lamC*(beta0*ccoeff::C1qq[idx]*(LQR))
						+lamE*(-pow(beta1/beta0,2)*ccoeff::C1qq[idx])
						+lamE*(beta2/beta0*ccoeff::C1qq[idx])
						+lamF*(1./3.*(pow(ccoeff::C1qq[idx],3)-3.*ccoeff::C1qq[idx]*ccoeff::C2qq[idx]+3.*ccoeff::C3qq[idx]))
						+lamG*(-beta1/pow(beta0,2)*(beta1*ccoeff::C1qq[idx] - beta0*pow(ccoeff::C1qq[idx],2) + 2.*beta0*ccoeff::C2qq[idx]))
						+lamH*(pow(beta1/beta0,2)*ccoeff::C1qq[idx])
						+lamG*(2.*beta1*ccoeff::C1qq[idx]*LQR)
						+lamF*((-beta1*ccoeff::C1qq[idx]+beta0*pow(ccoeff::C1qq[idx],2)-2.*beta0*ccoeff::C2qq[idx])*LQR)
						+lamF*(pow(beta0,2)*ccoeff::C1qq[idx]*pow(LQR,2))
						);

				complex <double> C1,C2,C3,rootC,n2llaexp;
				complex <double> lam = xlambda;

				//cdefexp(ccoeff::C1qg[idx],   ccoeff::C2qg[idx],   ccoeff::C3qg[idx]); aexpqg[idx]   = n3llaexp();
				cdefexp(ccoeff::C1qg[idx],   ccoeff::C2qg[idx],   0.); aexpqg[idx]   = n3llaexp();

				//aexpqg[idx]   = (1.-xlambda/(1.+aass*ccoeff::C2qg[idx]/ccoeff::C1qg[idx]))/(1.-xlambda);

				aexpqqb[idx]   = (1.-xlambda/(1.+aass*ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx]))/(1.-xlambda);
				aexpqqp[idx]   = (1.-xlambda/(1.+aass*ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx]))/(1.-xlambda);
				aexpqqbp[idx]   = (1.-xlambda/(1.+aass*ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx]))/(1.-xlambda);
			}
			}


		////gint::allocate();
		////gint::calc(b);
		//cout << endl;
		//cout << " gint::logasl "   << gint::logasl_expc      << " expc " << logasl      << endl;
		//cout << " gint::alogqq "   << gint::alogqq[0]   << " expc " << log(aexpqq[0])   << endl;
		//cout << " gint::alogqg "   << gint::alogqg[0]   << " expc " << log(aexpqg[0])   << endl;
		//cout << " gint::alogqqb "  << gint::alogqqb[0]  << " expc " << log(aexpqqb[0])  << endl;
		//cout << " gint::alogqqp "  << gint::alogqqp[0]  << " expc " << log(aexpqqp[0])  << endl;
		//cout << " gint::alogqqbp " << gint::alogqqbp[0] << " expc " << log(aexpqqbp[0]) << endl;
		////gint::free();


		for (int sign = mesq::positive; sign <= mesq::negative; sign++)
			for (int i = 0; i < mellinint::mdim; i++)
			{
			int idx = anomalous::index(i,sign);
			int ii = hcoefficients::index(i,sign);
			qqb[ii] = aexpqq[idx] * aexpqq[idx];
			qg[ii]  = aexp * aexpqg[idx]*aexpqq[idx];
			qq[ii]  = aexpqq[idx] * aexpqqb[idx] * aexp2;
			qqp[ii] = aexpqq[idx] * aexpqqbp[idx] * aexp2;
			qqbp[ii] = aexpqq[idx] * aexpqqp[idx] * aexp2;
			gg[ii]  = aexp*aexpqg[idx] * aexp*aexpqg[idx];
			qbg[ii] = aexp * aexpqg[idx] *aexp2*aexpqqb[idx];
			qpg[ii] = aexp * aexpqg[idx] *aexp2*aexpqqp[idx];
			qbpg[ii] = aexp * aexpqg[idx] *aexp2*aexpqqbp[idx];
			if (isnan_ofast(real(qg[ii])) || isnan_ofast(imag(qg[ii])))
			{
				int idx = anomalous::index(i,sign);
				cout << "Warning, expc::qg[" << i << "] = " << qg[ii] << ", b = "  << b << ", lambda = " << xlambda << " logasl " << logasl
				<< " aexp " << aexp << " aexpqq " << aexpqq[idx] << " aexpqg " << aexpqg[idx]
				<< endl;
			}
			}
		//include aexp into aexpqg and aexpqq,qqp,qqb,qqbp
		for (int sign = mesq::positive; sign <= mesq::negative; sign++)
			for (int i = 0; i < mellinint::mdim; i++)
			{
			int idx = anomalous::index(i,sign);
			aexpqg[idx]   *= aexp;
			aexpqqb[idx]  *= aexp2;
			aexpqqp[idx]  *= aexp2;
			aexpqqbp[idx] *= aexp2;
			}
		}
		}


	if (opts.order_expc >= 5 || opts.numexpc)
		{
		if (opts.mellin1d)
		{
		complex <double> alogqq[mellinint::mdim*2];
		complex <double> alogqg[mellinint::mdim*2];
		complex <double> alogqqb[mellinint::mdim*2];
		complex <double> alogqqp[mellinint::mdim*2];
		complex <double> alogqqbp[mellinint::mdim*2];
		//complex <double> aexpqq[mellinint::mdim*2];
		//complex <double> aexpqg[mellinint::mdim*2];
		//complex <double> aexpqqb[mellinint::mdim*2];
		//complex <double> aexpqqp[mellinint::mdim*2];
		//complex <double> aexpqqbp[mellinint::mdim*2];

		if (opts.expc == 0)
			if (opts.sumlogs)
			{
			fill(alogqq, alogqq+mellinint::mdim*2, 0.);
			fill(alogqg, alogqq+mellinint::mdim*2, 0.);
			fill(alogqqb, alogqqb+mellinint::mdim*2, 0.);
			fill(alogqqp, alogqqp+mellinint::mdim*2, 0.);
			fill(alogqqbp, alogqqbp+mellinint::mdim*2, 0.);
			}
			else
			{
			fill(aexpqq, aexpqq+mellinint::mdim*2, 1.);
			fill(aexpqg, aexpqg+mellinint::mdim*2, 1.);
			fill(aexpqqb, aexpqqb+mellinint::mdim*2, 1.);
			fill(aexpqqp, aexpqqp+mellinint::mdim*2, 1.);
			fill(aexpqqbp, aexpqqbp+mellinint::mdim*2, 1.);
			}

		//!!! gint::alog are evaluated with blim::sudakov, not with blim:expc !!!
		if (opts.numexpc)
			for (int sign = mesq::positive; sign <= mesq::negative; sign++)
			for (int i = 0; i < mellinint::mdim; i++)
			{
			int idx = anomalous::index(i,sign);

			if (opts.sumlogs)
				{
				alogqq[idx]   = gint::alogqq[idx];
				alogqg[idx]   = gint::alogqg[idx];
				alogqqb[idx]  = gint::alogqqb[idx];
				alogqqp[idx]  = gint::alogqqp[idx];
				alogqqbp[idx] = gint::alogqqbp[idx];
				}
			else
				{
				aexpqq[idx]   = exp(gint::alogqq[idx]);
				aexpqg[idx]   = exp(gint::alogqg[idx]);
				aexpqqb[idx]  = exp(gint::alogqqb[idx]);
				aexpqqp[idx]  = exp(gint::alogqqp[idx]);
				aexpqqbp[idx] = exp(gint::alogqqbp[idx]);
				}
			}

		for (int sign = mesq::positive; sign <= mesq::negative; sign++)
			for (int i = 0; i < mellinint::mdim; i++)
			{
			int idx = anomalous::index(i,sign);
			int ii = hcoefficients::index(i,sign);
			if (opts.sumlogs)
			{
				qqb[ii] = exp(alogqq[idx]+alogqq[idx]+sudakov::logS);
				qg[ii]  = exp(-logasl+alogqg[idx]+alogqq[idx]+sudakov::logS);
				qq[ii]  = exp(alogqq[idx]+alogqqb[idx]-2.*logasl+sudakov::logS);
				qqp[ii] = exp(alogqq[idx]+alogqqbp[idx]-2.*logasl+sudakov::logS);
				qqbp[ii] = exp(alogqq[idx]+alogqqp[idx]-2.*logasl+sudakov::logS);
				gg[ii]   = exp(-logasl+alogqg[idx]-logasl+alogqg[idx]+sudakov::logS);
				qbg[ii]  = exp(-logasl+alogqg[idx]-2.*logasl+alogqqb[idx]+sudakov::logS);
				qpg[ii]  = exp(-logasl+alogqg[idx]-2.*logasl+alogqqp[idx]+sudakov::logS);
				qbpg[ii] = exp(-logasl+alogqg[idx]-2.*logasl+alogqqbp[idx]+sudakov::logS);
			}
			else
			{
				qqb[ii] = aexpqq[idx] * aexpqq[idx];
				qg[ii]  = aexp * aexpqg[idx]*aexpqq[idx];
				qq[ii]  = aexpqq[idx] * aexpqqb[idx] * aexp2;
				qqp[ii] = aexpqq[idx] * aexpqqbp[idx] * aexp2;
				qqbp[ii] = aexpqq[idx] * aexpqqp[idx] * aexp2;
				gg[ii]  = aexp*aexpqg[idx] * aexp*aexpqg[idx];
				qbg[ii] = aexp * aexpqg[idx] *aexp2*aexpqqb[idx];
				qpg[ii] = aexp * aexpqg[idx] *aexp2*aexpqqp[idx];
				qbpg[ii] = aexp * aexpqg[idx] *aexp2*aexpqqbp[idx];
			}
			}
		//include aexp into aexpqg and aexpqq,qqp,qqb,qqbp
		for (int sign = mesq::positive; sign <= mesq::negative; sign++)
			for (int i = 0; i < mellinint::mdim; i++)
			{
			int idx = anomalous::index(i,sign);
			aexpqg[idx]   *= aexp;
			aexpqqb[idx]  *= aexp2;
			aexpqqp[idx]  *= aexp2;
			aexpqqbp[idx] *= aexp2;
			}
		}
		else //mellin2d
		{
		if (opts.expc == 0)
			{
			fill(aexpqq_1, aexpqq_1+mellinint::mdim*2, 1.);
			fill(aexpqg_1, aexpqg_1+mellinint::mdim*2, 1.);
			fill(aexpqqb_1, aexpqqb_1+mellinint::mdim*2, 1.);
			fill(aexpqqp_1, aexpqqp_1+mellinint::mdim*2, 1.);
			fill(aexpqqbp_1, aexpqqbp_1+mellinint::mdim*2, 1.);
			fill(aexpqq_2, aexpqq_2+mellinint::mdim*2, 1.);
			fill(aexpqg_2, aexpqg_2+mellinint::mdim*2, 1.);
			fill(aexpqqb_2, aexpqqb_2+mellinint::mdim*2, 1.);
			fill(aexpqqp_2, aexpqqp_2+mellinint::mdim*2, 1.);
			fill(aexpqqbp_2, aexpqqbp_2+mellinint::mdim*2, 1.);
			}

		if (opts.numexpc)
			for (int sign = mesq::positive; sign <= mesq::negative; sign++)
			for (int i = 0; i < mellinint::mdim; i++)
			{
			int idx = anomalous::index(i,sign);

			aexpqq_1[idx]   = exp(gint::alogqq_1[idx]);
			aexpqg_1[idx]   = exp(gint::alogqg_1[idx]);
			aexpqqb_1[idx]  = exp(gint::alogqqb_1[idx]);
			aexpqqp_1[idx]  = exp(gint::alogqqp_1[idx]);
			aexpqqbp_1[idx] = exp(gint::alogqqbp_1[idx]);
			aexpqq_2[idx]   = exp(gint::alogqq_2[idx]);
			aexpqg_2[idx]   = exp(gint::alogqg_2[idx]);
			aexpqqb_2[idx]  = exp(gint::alogqqb_2[idx]);
			aexpqqp_2[idx]  = exp(gint::alogqqp_2[idx]);
			aexpqqbp_2[idx] = exp(gint::alogqqbp_2[idx]);
			}
		}
		}


	//check nans
	if (opts.mellin1d)
		for (int sign = mesq::positive; sign <= mesq::negative; sign++)
		for (int i = 0; i < mellinint::mdim; i++)
		{
		int ii = hcoefficients::index(i,sign);
		if (isnan_ofast(real(qg[ii])) || isnan_ofast(imag(qg[ii])))
			{
			int idx = anomalous::index(i,sign);
			cout << "Warning, expc::qg[" << i << "] = " << qg[ii] << ", b = "  << b << ", lambda = " << xlambda << " logasl " << logasl
			//<< " aexpqq " << aexpqq[idx] << " aexpqg " << aexpqg[idx]
			<< endl;
			if (opts.sumlogs)
			qg[ii] = sudakov::S;
			else
			qg[ii] = 1.;
			}
		}
	}
}

void expc::noexp()
{
  if (opts.mellin1d)
    {
      fill(aexpqq, aexpqq+mellinint::mdim*2, 1.);
      fill(aexpqg, aexpqg+mellinint::mdim*2, 1.);
      fill(aexpqqb, aexpqqb+mellinint::mdim*2, 1.);
      fill(aexpqqp, aexpqqp+mellinint::mdim*2, 1.);
      fill(aexpqqbp, aexpqqbp+mellinint::mdim*2, 1.);
    }
  else
    {
      fill(aexpqq_1, aexpqq_1+mellinint::mdim*2, 1.);
      fill(aexpqg_1, aexpqg_1+mellinint::mdim*2, 1.);
      fill(aexpqqb_1, aexpqqb_1+mellinint::mdim*2, 1.);
      fill(aexpqqp_1, aexpqqp_1+mellinint::mdim*2, 1.);
      fill(aexpqqbp_1, aexpqqbp_1+mellinint::mdim*2, 1.);
      fill(aexpqq_2, aexpqq_2+mellinint::mdim*2, 1.);
      fill(aexpqg_2, aexpqg_2+mellinint::mdim*2, 1.);
      fill(aexpqqb_2, aexpqqb_2+mellinint::mdim*2, 1.);
      fill(aexpqqp_2, aexpqqp_2+mellinint::mdim*2, 1.);
      fill(aexpqqbp_2, aexpqqbp_2+mellinint::mdim*2, 1.);
    }
}

void expc::delta()
{
  if (opts.order_expc == 1)
    return;
 
  double LQR = LR-LQ;
  if (opts.mellin1d)
    {
      if (opts.order_expc == 2)
	for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	  for (int i = 0; i < mellinint::mdim; i++)
	    {
	      int idx = anomalous::index(i,sign);
	      aexpqq[idx] = pow(aexpB,ccoeff::C1qq_delta);
	    }
      else if (opts.order_expc == 3)
	for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	  for (int i = 0; i < mellinint::mdim; i++)
	    {
	      int idx = anomalous::index(i,sign);
	      aexpqq[idx] = pow(aexpB,ccoeff::C1qq_delta)
		*pow(aexpC,(0.5*pow(ccoeff::C1qq_delta,2) - ccoeff::C2qq_delta))
		*pow(aexpD,-beta1/beta0*ccoeff::C1qq_delta)
		*pow(aexpC,beta0*ccoeff::C1qq_delta*(LQR))
		;
	    }
      else if (opts.order_expc == 4)
	for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	  for (int i = 0; i < mellinint::mdim; i++)
	    {
	      int idx = anomalous::index(i,sign);
	      aexpqq[idx] = pow(aexpB,ccoeff::C1qq_delta)
		*pow(aexpC,(0.5*pow(ccoeff::C1qq_delta,2) - ccoeff::C2qq_delta))
		*pow(aexpD,-beta1/beta0*ccoeff::C1qq_delta)
		*pow(aexpC,beta0*ccoeff::C1qq_delta*LQR)
		*pow(aexpE,pow(beta1/beta0,2)*ccoeff::C1qq_delta)
		*pow(aexpE,-beta2/beta0*ccoeff::C1qq_delta)
		*pow(aexpF,-1./3.*(pow(ccoeff::C1qq_delta,3)-3.*ccoeff::C1qq_delta*ccoeff::C2qq_delta+3.*ccoeff::C3qq_delta))
		*pow(aexpG,beta1/pow(beta0,2)*(beta1*ccoeff::C1qq_delta - beta0*pow(ccoeff::C1qq_delta,2) + 2.*beta0*ccoeff::C2qq_delta))
		*pow(aexpH,-pow(beta1/beta0,2)*ccoeff::C1qq_delta)
		*pow(aexpG,2.*beta1*ccoeff::C1qq_delta*LQR)
		*pow(aexpF,(-beta1*ccoeff::C1qq_delta+beta0*pow(ccoeff::C1qq_delta,2)-2.*beta0*ccoeff::C2qq_delta)*LQR)
		*pow(aexpF,pow(beta0,2)*ccoeff::C1qq_delta*pow(LQR,2))
		;
	    }
    }
  else
    {
      if (opts.order_expc == 2)
	for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	  for (int i = 0; i < mellinint::mdim; i++)
	    {
	      int idx = anomalous::index(i,sign);
	      aexpqq_1[idx] = pow(aexpB,ccoeff::C1qq_delta);
	      aexpqq_2[idx] = pow(aexpB,ccoeff::C1qq_delta);
	    }
      else if (opts.order_expc == 3)
	for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	  for (int i = 0; i < mellinint::mdim; i++)
	    {
	      int idx = anomalous::index(i,sign);
	      aexpqq_1[idx] = pow(aexpB,ccoeff::C1qq_delta)
		*pow(aexpC,(0.5*pow(ccoeff::C1qq_delta,2) - ccoeff::C2qq_delta))
		*pow(aexpD,-beta1/beta0*ccoeff::C1qq_delta)
		*pow(aexpC,beta0*ccoeff::C1qq_delta*LQR)
		;
	      aexpqq_2[idx] = pow(aexpB,ccoeff::C1qq_delta)
		*pow(aexpC,(0.5*pow(ccoeff::C1qq_delta,2) - ccoeff::C2qq_delta))
		*pow(aexpD,-beta1/beta0*ccoeff::C1qq_delta)
		*pow(aexpC,beta0*ccoeff::C1qq_delta*LQR)
		;
	    }
      else if (opts.order_expc == 4)
	for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	  for (int i = 0; i < mellinint::mdim; i++)
	    {
	      int idx = anomalous::index(i,sign);
	      aexpqq_1[idx] = pow(aexpB,ccoeff::C1qq_delta)
		*pow(aexpC,(0.5*pow(ccoeff::C1qq_delta,2) - ccoeff::C2qq_delta))
		*pow(aexpD,-beta1/beta0*ccoeff::C1qq_delta)
		*pow(aexpC,beta0*ccoeff::C1qq_delta*LQR)
		*pow(aexpE,pow(beta1/beta0,2)*ccoeff::C1qq_delta)
		*pow(aexpE,-beta2/beta0*ccoeff::C1qq_delta)
		*pow(aexpF,-1./3.*(pow(ccoeff::C1qq_delta,3)-3.*ccoeff::C1qq_delta*ccoeff::C2qq_delta+3.*ccoeff::C3qq_delta))
		*pow(aexpG,beta1/pow(beta0,2)*(beta1*ccoeff::C1qq_delta - beta0*pow(ccoeff::C1qq_delta,2) + 2.*beta0*ccoeff::C2qq_delta))
		*pow(aexpH,-pow(beta1/beta0,2)*ccoeff::C1qq_delta)
		*pow(aexpG,2.*beta1*ccoeff::C1qq_delta*LQR)
		*pow(aexpF,(-beta1*ccoeff::C1qq_delta+beta0*pow(ccoeff::C1qq_delta,2)-2.*beta0*ccoeff::C2qq_delta)*LQR)
		*pow(aexpF,pow(beta0,2)*ccoeff::C1qq_delta*pow(LQR,2))
		;
	      aexpqq_2[idx] = pow(aexpB,ccoeff::C1qq_delta)
		*pow(aexpC,(0.5*pow(ccoeff::C1qq_delta,2) - ccoeff::C2qq_delta))
		*pow(aexpD,-beta1/beta0*ccoeff::C1qq_delta)
		*pow(aexpC,beta0*ccoeff::C1qq_delta*LQR)
		*pow(aexpE,pow(beta1/beta0,2)*ccoeff::C1qq_delta)
		*pow(aexpE,-beta2/beta0*ccoeff::C1qq_delta)
		*pow(aexpF,-1./3.*(pow(ccoeff::C1qq_delta,3)-3.*ccoeff::C1qq_delta*ccoeff::C2qq_delta+3.*ccoeff::C3qq_delta))
		*pow(aexpG,beta1/pow(beta0,2)*(beta1*ccoeff::C1qq_delta - beta0*pow(ccoeff::C1qq_delta,2) + 2.*beta0*ccoeff::C2qq_delta))
		*pow(aexpH,-pow(beta1/beta0,2)*ccoeff::C1qq_delta)
		*pow(aexpG,2.*beta1*ccoeff::C1qq_delta*LQR)
		*pow(aexpF,(-beta1*ccoeff::C1qq_delta+beta0*pow(ccoeff::C1qq_delta,2)-2.*beta0*ccoeff::C2qq_delta)*LQR)
		*pow(aexpF,pow(beta0,2)*ccoeff::C1qq_delta*pow(LQR,2))
		;
	    }
    }
}

void expc::diag()
{
  if (opts.order_expc == 1)
    return;
  
  double LQR = LR-LQ;
  if (opts.mellin1d)
    {
      if (opts.order_expc == 2)
	for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	  for (int i = 0; i < mellinint::mdim; i++)
	    {
	      int idx = anomalous::index(i,sign);
	      aexpqq[idx] = pow(aexpB,ccoeff::C1qq[idx]);
	    }
      else if (opts.order_expc == 3)
	for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	  for (int i = 0; i < mellinint::mdim; i++)
	    {
	      int idx = anomalous::index(i,sign);
	      aexpqq[idx] = exp(
				+lamB*(ccoeff::C1qq[idx])
				+lamC*((0.5*pow(ccoeff::C1qq[idx],2) - ccoeff::C2qq[idx]))
				+lamD*(-beta1/beta0*ccoeff::C1qq[idx])
				+lamC*(beta0*ccoeff::C1qq[idx]*(LQR))
				);
	    }
      else if (opts.order_expc == 4)
	for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	  for (int i = 0; i < mellinint::mdim; i++)
	    {
	      int idx = anomalous::index(i,sign);
	      complex <double> C1,C2,C3;	      
	      C1 = ccoeff::C1qq[idx]; C2 = ccoeff::C2qq[idx]; C3 = ccoeff::C3qq[idx];
	      aexpqq[idx] = exp(
				+lamB*(C1)
				+lamC*((0.5*pow(C1,2) - C2))
				+lamD*(-beta1/beta0*C1)
				+lamC*(beta0*C1*(LQR))
				+lamE*(-pow(beta1/beta0,2)*C1)
				+lamE*(beta2/beta0*C1)
				+lamF*(1./3.*(pow(C1,3)-3.*C1*C2+3.*C3))
				+lamG*(-beta1/pow(beta0,2)*(beta1*C1 - beta0*pow(C1,2) + 2.*beta0*C2))
				+lamH*(pow(beta1/beta0,2)*C1)
				+lamG*(2.*beta1*C1*LQR)
				+lamF*((-beta1*C1+beta0*pow(C1,2)-2.*beta0*C2)*LQR)
				+lamF*(pow(beta0,2)*C1*pow(LQR,2))
				);
	    }
    }
  //mellin 2d
  else
    {
      if (opts.order_expc == 2)
	for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	  for (int i = 0; i < mellinint::mdim; i++)
	    {
	      int idx = anomalous::index(i,sign);
	      aexpqq_1[idx] = pow(aexpB,ccoeff::C1qq_1[idx]);
	      aexpqq_2[idx] = pow(aexpB,ccoeff::C1qq_2[idx]);
	    }
      else if (opts.order_expc == 3)
	for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	  for (int i = 0; i < mellinint::mdim; i++)
	    {
	      int idx = anomalous::index(i,sign);
	      aexpqq_1[idx] = exp(
				+lamB*(ccoeff::C1qq_1[idx])
				+lamC*((0.5*pow(ccoeff::C1qq_1[idx],2) - ccoeff::C2qq_1[idx]))
				+lamD*(-beta1/beta0*ccoeff::C1qq_1[idx])
				+lamC*(beta0*ccoeff::C1qq_1[idx]*(LQR))
				);
	      aexpqq_2[idx] = exp(
				+lamB*(ccoeff::C1qq_2[idx])
				+lamC*((0.5*pow(ccoeff::C1qq_2[idx],2) - ccoeff::C2qq_2[idx]))
				+lamD*(-beta1/beta0*ccoeff::C1qq_2[idx])
				+lamC*(beta0*ccoeff::C1qq_2[idx]*(LQR))
				);
	    }
      else if (opts.order_expc == 4)
	for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	  for (int i = 0; i < mellinint::mdim; i++)
	    {
	      int idx = anomalous::index(i,sign);
	      complex <double> C1,C2,C3;	      
	      C1 = ccoeff::C1qq_1[idx]; C2 = ccoeff::C2qq_1[idx]; C3 = ccoeff::C3qq_1[idx];
	      aexpqq_1[idx] = exp(
				+lamB*(C1)
				+lamC*((0.5*pow(C1,2) - C2))
				+lamD*(-beta1/beta0*C1)
				+lamC*(beta0*C1*(LQR))
				+lamE*(-pow(beta1/beta0,2)*C1)
				+lamE*(beta2/beta0*C1)
				+lamF*(1./3.*(pow(C1,3)-3.*C1*C2+3.*C3))
				+lamG*(-beta1/pow(beta0,2)*(beta1*C1 - beta0*pow(C1,2) + 2.*beta0*C2))
				+lamH*(pow(beta1/beta0,2)*C1)
				+lamG*(2.*beta1*C1*LQR)
				+lamF*((-beta1*C1+beta0*pow(C1,2)-2.*beta0*C2)*LQR)
				+lamF*(pow(beta0,2)*C1*pow(LQR,2))
				);
	      C1 = ccoeff::C1qq_2[idx]; C2 = ccoeff::C2qq_2[idx]; C3 = ccoeff::C3qq_2[idx];
	      aexpqq_2[idx] = exp(
				+lamB*(C1)
				+lamC*((0.5*pow(C1,2) - C2))
				+lamD*(-beta1/beta0*C1)
				+lamC*(beta0*C1*(LQR))
				+lamE*(-pow(beta1/beta0,2)*C1)
				+lamE*(beta2/beta0*C1)
				+lamF*(1./3.*(pow(C1,3)-3.*C1*C2+3.*C3))
				+lamG*(-beta1/pow(beta0,2)*(beta1*C1 - beta0*pow(C1,2) + 2.*beta0*C2))
				+lamH*(pow(beta1/beta0,2)*C1)
				+lamG*(2.*beta1*C1*LQR)
				+lamF*((-beta1*C1+beta0*pow(C1,2)-2.*beta0*C2)*LQR)
				+lamF*(pow(beta0,2)*C1*pow(LQR,2))
				);
	    }
      
    }
}
 
void expc::css(complex <double> b)
{
  double blim = blim::expc;
  complex <double> bstar;
  if (opts.bprescription == 0 || opts.bprescription == 4 || opts.bstar_expc)
    //bstar = real(b)/sqrt(1.+pow(real(b)/blim,2));
    bstar = bstar::calc(b,blim);
  else
    bstar = b;
  complex <double> mubstar = resconst::b0/bstar;
  /*
  complex <double> mubstartilde;
  if (!opts.modlog)
    mubstartilde = mubstar;
  else if (opts.p == 1)
    mubstartilde = mubstar * scales::res / sqrt((pow(mubstar,2) + pow(scales::res,2)));
  else
    mubstartilde = mubstar * scales::res / pow((pow(mubstar,2*opts.p) + pow(scales::res,2*opts.p)),0.5/opts.p);
  */
  complex <double> mubstartilde = btilde::mubtl(bstar,scales::res);
    
  //iterative alphas
  complex <double> asmub = alphas::calc(mubstartilde, opts.order_expc+1);
  complex <double> asmubLO = alphas::asLO;
  complex <double> asmubNLO = alphas::asNLO;

  //numerical alphas
  //alphas::rgkt5(mubstartilde, scales::ren, pdf::alphas(scales::ren, opts.cmw));   
  //complex <double> asmub = alphas::asNLO/M_PI;

  //lhapdf alphas
  //complex <double> asmub = pdf::extalphas(real(mubstartilde))/M_PI;
  if (opts.mellin1d)
    {
      if (opts.order_expc == 1)
	for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	  for (int i = 0; i < mellinint::mdim; i++)
	    {
	      int idx = anomalous::index(i,sign);
	      aexpqq[idx] = (1. + asmub*(ccoeff::C1qq[idx]))/(1. + aass*ccoeff::C1qq[idx]);
	      aexpqg[idx] = asmub/aass;
	    }
      else if (opts.order_expc == 2)
	for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	  for (int i = 0; i < mellinint::mdim; i++)
	    {
	      int idx = anomalous::index(i,sign);
	      aexpqq[idx]   = (1. + asmub*ccoeff::C1qq[idx] + pow(asmub,2)*ccoeff::C2qq[idx])/(1. + aass*ccoeff::C1qq[idx] + pow(aass,2)*ccoeff::C2qq[idx]);
	      aexpqg[idx]   = (1. + asmub*ccoeff::C2qg[idx]/ccoeff::C1qg[idx])/(1. + aass*ccoeff::C2qg[idx]/ccoeff::C1qg[idx]);
	      aexpqqb[idx]  = 1.;//pow(asmub,2)/pow(aass,2);
	      aexpqqp[idx]  = 1.;//pow(asmub,2)/pow(aass,2);
	      aexpqqbp[idx] = 1.;//pow(asmub,2)/pow(aass,2);
	    }
      else if (opts.order_expc >= 3)
	for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	  for (int i = 0; i < mellinint::mdim; i++)
	    {
	      int idx = anomalous::index(i,sign);
		  
	      aexpqq[idx]   = (1. + asmub*ccoeff::C1qq[idx] + pow(asmub,2)*ccoeff::C2qq[idx] + pow(asmub,3)*ccoeff::C3qq[idx])/(1. + aass*ccoeff::C1qq[idx] + pow(aass,2)*ccoeff::C2qq[idx] + pow(aass,3)*ccoeff::C3qq[idx]);
	      aexpqg[idx]   = (1. + asmub*ccoeff::C2qg[idx]/ccoeff::C1qg[idx] + pow(asmub,2)*ccoeff::C3qg[idx]/ccoeff::C1qg[idx])/(1. + aass*ccoeff::C2qg[idx]/ccoeff::C1qg[idx] + pow(aass,2)*ccoeff::C3qg[idx]/ccoeff::C1qg[idx]);
	      aexpqqb[idx]  = (1. + asmub*ccoeff::C3qqb[idx] /ccoeff::C2qqb[idx]  )/(1. + aass*ccoeff::C3qqb[idx] /ccoeff::C2qqb[idx]  );
	      aexpqqp[idx]  = (1. + asmub*ccoeff::C3qqp[idx] /ccoeff::C2qqp[idx]  )/(1. + aass*ccoeff::C3qqp[idx] /ccoeff::C2qqp[idx]  );
	      aexpqqbp[idx] = (1. + asmub*ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx] )/(1. + aass*ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx] );
	    }
    }
  else
    {
      if (opts.order_expc == 1)
	for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	  for (int i = 0; i < mellinint::mdim; i++)
	    {
	      int idx = anomalous::index(i,sign);
	      aexpqq_1[idx] = (1. + asmub*ccoeff::C1qq_1[idx])/(1. + aass*ccoeff::C1qq_1[idx]);
	      aexpqq_2[idx] = (1. + asmub*ccoeff::C1qq_2[idx])/(1. + aass*ccoeff::C1qq_2[idx]);
	      aexpqg_1[idx] = asmub/aass;
	      aexpqg_2[idx] = asmub/aass;
	    }
      if (opts.order_expc == 2)
	for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	  for (int i = 0; i < mellinint::mdim; i++)
	    {
	      int idx = anomalous::index(i,sign);
	      aexpqq_1[idx]   = (1. + asmub*ccoeff::C1qq_1[idx] + pow(asmub,2)*ccoeff::C2qq_1[idx])/(1. + aass*ccoeff::C1qq_1[idx] + pow(aass,2)*ccoeff::C2qq_1[idx]);
	      aexpqg_1[idx]   = (1. + asmub*ccoeff::C2qg_1[idx]/ccoeff::C1qg_1[idx])/(1. + aass*ccoeff::C2qg_1[idx]/ccoeff::C1qg_1[idx]);
	      aexpqq_2[idx]   = (1. + asmub*ccoeff::C1qq_2[idx] + pow(asmub,2)*ccoeff::C2qq_2[idx])/(1. + aass*ccoeff::C1qq_2[idx] + pow(aass,2)*ccoeff::C2qq_2[idx]);
	      aexpqg_2[idx]   = (1. + asmub*ccoeff::C2qg_2[idx]/ccoeff::C1qg_2[idx])/(1. + aass*ccoeff::C2qg_2[idx]/ccoeff::C1qg_2[idx]);
	    }
      if (opts.order_expc >= 3)
	for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	  for (int i = 0; i < mellinint::mdim; i++)
	    {
	      int idx = anomalous::index(i,sign);
	      aexpqq_1[idx]   = (1. + asmub*ccoeff::C1qq_1[idx] + pow(asmub,2)*ccoeff::C2qq_1[idx] + pow(asmub,3)*ccoeff::C3qq_1[idx])/(1. + aass*ccoeff::C1qq_1[idx] + pow(aass,2)*ccoeff::C2qq_1[idx] + pow(aass,3)*ccoeff::C3qq_1[idx]);
	      aexpqg_1[idx]   = (1. + asmub*ccoeff::C2qg_1[idx]/ccoeff::C1qg_1[idx] + pow(asmub,2)*ccoeff::C3qg_1[idx]/ccoeff::C1qg_1[idx])/(1. + aass*ccoeff::C2qg_1[idx]/ccoeff::C1qg_1[idx] + pow(aass,2)*ccoeff::C3qg_1[idx]/ccoeff::C1qg_1[idx]);
	      aexpqqb_1[idx]  = (1. + asmub*ccoeff::C3qqb_1[idx] /ccoeff::C2qqb_1[idx]  )/(1. + aass*ccoeff::C3qqb_1[idx] /ccoeff::C2qqb_1[idx]  );
	      aexpqqp_1[idx]  = (1. + asmub*ccoeff::C3qqp_1[idx] /ccoeff::C2qqp_1[idx]  )/(1. + aass*ccoeff::C3qqp_1[idx] /ccoeff::C2qqp_1[idx]  );
	      aexpqqbp_1[idx] = (1. + asmub*ccoeff::C3qqbp_1[idx]/ccoeff::C2qqbp_1[idx] )/(1. + aass*ccoeff::C3qqbp_1[idx]/ccoeff::C2qqbp_1[idx] );
	      aexpqq_2[idx]   = (1. + asmub*ccoeff::C1qq_2[idx] + pow(asmub,2)*ccoeff::C2qq_2[idx] + pow(asmub,3)*ccoeff::C3qq_2[idx])/(1. + aass*ccoeff::C1qq_2[idx] + pow(aass,2)*ccoeff::C2qq_2[idx] + pow(aass,3)*ccoeff::C3qq_2[idx]);
	      aexpqg_2[idx]   = (1. + asmub*ccoeff::C2qg_2[idx]/ccoeff::C1qg_2[idx] + pow(asmub,2)*ccoeff::C3qg_2[idx]/ccoeff::C1qg_2[idx])/(1. + aass*ccoeff::C2qg_2[idx]/ccoeff::C1qg_2[idx] + pow(aass,2)*ccoeff::C3qg_2[idx]/ccoeff::C1qg_2[idx]);
	      aexpqqb_2[idx]  = (1. + asmub*ccoeff::C3qqb_2[idx] /ccoeff::C2qqb_2[idx]  )/(1. + aass*ccoeff::C3qqb_2[idx] /ccoeff::C2qqb_2[idx]  );
	      aexpqqp_2[idx]  = (1. + asmub*ccoeff::C3qqp_2[idx] /ccoeff::C2qqp_2[idx]  )/(1. + aass*ccoeff::C3qqp_2[idx] /ccoeff::C2qqp_2[idx]  );
	      aexpqqbp_2[idx] = (1. + asmub*ccoeff::C3qqbp_2[idx]/ccoeff::C2qqbp_2[idx] )/(1. + aass*ccoeff::C3qqbp_2[idx]/ccoeff::C2qqbp_2[idx] );
	    }
    }
}

//Standard formula for opts.expc == 3
void expc::full()
{
  if (opts.order_expc == 1)
    return;
  
  double regfac = opts.expcreg; //Regularisation factor

  double LQR = LR-LQ;
  if (opts.mellin1d)
    {
      if (opts.order_expc == 2)
	for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	  for (int i = 0; i < mellinint::mdim; i++)
	    {
	      int idx = anomalous::index(i,sign);
	      aexpqq[idx] = pow(aexpB,ccoeff::C1qq[idx]);
	      complex <double> C1 = (ccoeff::C2qg[idx]+regfac*ccoeff::C1qg_reg[idx])/(ccoeff::C1qg[idx]-aass*regfac*ccoeff::C1qg_reg[idx]);
	      aexpqg[idx] = pow(aexpB,C1);
	      aexpqqb[idx] = pow(aexpB,ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx]);
	      aexpqqp[idx] = pow(aexpB,ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx]);
	      aexpqqbp[idx] = pow(aexpB,ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx]);

	      aexpqg_bb[idx] = pow(aexpB_bb,C1);
	    }
      else if (opts.order_expc == 3)
	for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	  for (int i = 0; i < mellinint::mdim; i++)
	    {
	      int idx = anomalous::index(i,sign);
	      complex <double> C1,C2;
	      C1 = ccoeff::C1qq[idx]; C2 = ccoeff::C2qq[idx];
	      aexpqq[idx] = exp(
				+lamB*(C1)
				+lamC*((0.5*pow(C1,2) - C2))
				+lamD*(-beta1/beta0*C1)
				+lamC*(beta0*C1*(LQR))
				);

	      //C1 = ccoeff::C2qg[idx]/ccoeff::C1qg[idx]; C2 = ccoeff::C3qg[idx]/ccoeff::C1qg[idx];
	      C1 = (ccoeff::C2qg[idx]+regfac*ccoeff::C1qg_reg[idx])/(ccoeff::C1qg[idx]-aass*regfac*ccoeff::C1qg_reg[idx]);
	      C2 = ccoeff::C3qg[idx]/(ccoeff::C1qg[idx]-aass*regfac*ccoeff::C1qg_reg[idx]);

	//cout << idx
	//	   << " non reg " << ccoeff::C2qg[idx]/ccoeff::C1qg[idx]
	//	   << " reg " << C1 << endl;
	      
	      aexpqg[idx] = exp(+lamB*(C1)
				+lamC*((0.5*pow(C1,2) - C2))
				+lamD*(-beta1/beta0*C1)
				+lamC*(beta0*C1*(LQR)));
	      C1 = ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx]; C2 = - (H1q/2.*C1+(H2q/2.-pow(H1q,2)/8.)); //--> scheme invariance for C4/C2
	      aexpqqb[idx] = exp(+lamB*(C1)
				+lamC*((0.5*pow(C1,2) - C2))
				+lamD*(-beta1/beta0*C1)
				+lamC*(beta0*C1*(LQR)));
	      C1 = ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx]; C2 = - (H1q/2.*C1+(H2q/2.-pow(H1q,2)/8.)); //--> scheme invariance for C4/C2
	      aexpqqp[idx] = exp(+lamB*(C1)
				+lamC*((0.5*pow(C1,2) - C2))
				+lamD*(-beta1/beta0*C1)
				+lamC*(beta0*C1*(LQR)));
	      C1 = ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx]; C2 = - (H1q/2.*C1+(H2q/2.-pow(H1q,2)/8.)); //--> scheme invariance for C4/C2
	      aexpqqbp[idx] = exp(+lamB*(C1)
				+lamC*((0.5*pow(C1,2) - C2))
				+lamD*(-beta1/beta0*C1)
				+lamC*(beta0*C1*(LQR)));
	    }
      else if (opts.order_expc == 4)
	for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	  for (int i = 0; i < mellinint::mdim; i++)
	    { 
	      int idx = anomalous::index(i,sign);

	      complex <double> C1,C2,C3;	      
	      
	      C1 = ccoeff::C1qq[idx]; C2 = ccoeff::C2qq[idx]; C3 = ccoeff::C3qq[idx];
	     
	      aexpqq[idx] = exp(
				+lamB*(C1)
				+lamC*((0.5*pow(C1,2) - C2))
				+lamD*(-beta1/beta0*C1)
				+lamC*(beta0*C1*(LQR))
				+lamE*(-pow(beta1/beta0,2)*C1)
				+lamE*(beta2/beta0*C1)
				+lamF*(1./3.*(pow(C1,3)-3.*C1*C2+3.*C3))
				+lamG*(-beta1/pow(beta0,2)*(beta1*C1 - beta0*pow(C1,2) + 2.*beta0*C2))
				+lamH*(pow(beta1/beta0,2)*C1)
				+lamG*(2.*beta1*C1*LQR)
				+lamF*((-beta1*C1+beta0*pow(C1,2)-2.*beta0*C2)*LQR)
				+lamF*(pow(beta0,2)*C1*pow(LQR,2))
				);

	      C1 = (ccoeff::C2qg[idx]+regfac*ccoeff::C1qg_reg[idx])/(ccoeff::C1qg[idx]-aass*regfac*ccoeff::C1qg_reg[idx]);
	      C2 = ccoeff::C3qg[idx]/(ccoeff::C1qg[idx]-aass*regfac*ccoeff::C1qg_reg[idx]);
	      //C3 = - H1q/2.*C2  - (H2q/2.-pow(H1q,2)/8.)*C1  - (H3q/2.+pow(H1q,3)/16.);
	      C3 = ccoeff::C4qg[idx]/(ccoeff::C1qg[idx]-aass*regfac*ccoeff::C1qg_reg[idx]);
	      
	      //C1 = ccoeff::C2qg[idx]/ccoeff::C1qg[idx]; C2 = ccoeff::C3qg[idx]/ccoeff::C1qg[idx]; C3 = ccoeff::C4qg[idx]/ccoeff::C1qg[idx];

	      aexpqg[idx] = exp(
				+lamB*(C1)
				+lamC*((0.5*pow(C1,2) - C2))
				+lamD*(-beta1/beta0*C1)
				+lamC*(beta0*C1*(LQR))
				+lamE*(-pow(beta1/beta0,2)*C1)
				+lamE*(beta2/beta0*C1)
				+lamF*(1./3.*(pow(C1,3)-3.*C1*C2+3.*C3))
				+lamG*(-beta1/pow(beta0,2)*(beta1*C1 - beta0*pow(C1,2) + 2.*beta0*C2))
				+lamH*(pow(beta1/beta0,2)*C1)
				+lamG*(2.*beta1*C1*LQR)
				+lamF*((-beta1*C1+beta0*pow(C1,2)-2.*beta0*C2)*LQR)
				+lamF*(pow(beta0,2)*C1*pow(LQR,2))
				);
	      C1 = ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx]; C2 = - (H1q/2.*C1+(H2q/2.-pow(H1q,2)/8.)); C3 = - H1q/2.*C2  - (H2q/2.-pow(H1q,2)/8.)*C1  - (H3q/2.+pow(H1q,3)/16.-1./4.*H1q*H2q);
	      //C1 = ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx]; C2 = ccoeff::C4qqb[idx]/ccoeff::C2qqb[idx]; C3 = - H1q/2.*C2  - (H2q/2.-pow(H1q,2)/8.)*C1  - (H3q/2.+pow(H1q,3)/16.-1./4.*H1q*H2q);
	      aexpqqb[idx] = exp(
				 +lamB*(C1)
				 +lamC*((0.5*pow(C1,2) - C2))
				 +lamD*(-beta1/beta0*C1)
				 +lamC*(beta0*C1*(LQR))
				 +lamE*(-pow(beta1/beta0,2)*C1)
				 +lamE*(beta2/beta0*C1)
				 +lamF*(1./3.*(pow(C1,3)-3.*C1*C2+3.*C3))
				 +lamG*(-beta1/pow(beta0,2)*(beta1*C1 - beta0*pow(C1,2) + 2.*beta0*C2))
				 +lamH*(pow(beta1/beta0,2)*C1)
				 +lamG*(2.*beta1*C1*LQR)
				 +lamF*((-beta1*C1+beta0*pow(C1,2)-2.*beta0*C2)*LQR)
				 +lamF*(pow(beta0,2)*C1*pow(LQR,2))
				 );
		
	      C1 = ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx]; C2 = - (H1q/2.*C1+(H2q/2.-pow(H1q,2)/8.)); C3 = - H1q/2.*C2  - (H2q/2.-pow(H1q,2)/8.)*C1  - (H3q/2.+pow(H1q,3)/16.-1./4.*H1q*H2q);
	      //C1 = ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx]; C2 = ccoeff::C4qqp[idx]/ccoeff::C2qqp[idx]; C3 = - H1q/2.*C2  - (H2q/2.-pow(H1q,2)/8.)*C1  - (H3q/2.+pow(H1q,3)/16.-1./4.*H1q*H2q);
	      aexpqqp[idx] = exp(
				 +lamB*(C1)
				 +lamC*((0.5*pow(C1,2) - C2))
				 +lamD*(-beta1/beta0*C1)
				 +lamC*(beta0*C1*(LQR))
				 +lamE*(-pow(beta1/beta0,2)*C1)
				 +lamE*(beta2/beta0*C1)
				 +lamF*(1./3.*(pow(C1,3)-3.*C1*C2+3.*C3))
				 +lamG*(-beta1/pow(beta0,2)*(beta1*C1 - beta0*pow(C1,2) + 2.*beta0*C2))
				 +lamH*(pow(beta1/beta0,2)*C1)
				 +lamG*(2.*beta1*C1*LQR)
				 +lamF*((-beta1*C1+beta0*pow(C1,2)-2.*beta0*C2)*LQR)
				 +lamF*(pow(beta0,2)*C1*pow(LQR,2))
				 );
		
	      C1 = ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx]; C2 = - (H1q/2.*C1+(H2q/2.-pow(H1q,2)/8.)); C3 = - H1q/2.*C2  - (H2q/2.-pow(H1q,2)/8.)*C1  - (H3q/2.+pow(H1q,3)/16.-1./4.*H1q*H2q);
	      //C1 = ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx]; C2 = ccoeff::C4qqbp[idx]/ccoeff::C2qqbp[idx]; C3 = - H1q/2.*C2  - (H2q/2.-pow(H1q,2)/8.)*C1  - (H3q/2.+pow(H1q,3)/16.-1./4.*H1q*H2q);
	      aexpqqbp[idx] = exp(
				  +lamB*(C1)
				  +lamC*((0.5*pow(C1,2) - C2))
				  +lamD*(-beta1/beta0*C1)
				  +lamC*(beta0*C1*(LQR))
				  +lamE*(-pow(beta1/beta0,2)*C1)
				  +lamE*(beta2/beta0*C1)
				  +lamF*(1./3.*(pow(C1,3)-3.*C1*C2+3.*C3))
				  +lamG*(-beta1/pow(beta0,2)*(beta1*C1 - beta0*pow(C1,2) + 2.*beta0*C2))
				  +lamH*(pow(beta1/beta0,2)*C1)
				  +lamG*(2.*beta1*C1*LQR)
				  +lamF*((-beta1*C1+beta0*pow(C1,2)-2.*beta0*C2)*LQR)
				  +lamF*(pow(beta0,2)*C1*pow(LQR,2))
				 );
	    }
    }
  else
    {
      if (opts.order_expc == 2)
	for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	  for (int i = 0; i < mellinint::mdim; i++)
	    {
	      int idx = anomalous::index(i,sign);
	      complex <double> C1;
	      aexpqq_1[idx] = pow(aexpB,ccoeff::C1qq_1[idx]);
	      C1 = (ccoeff::C2qg_1[idx]+regfac*ccoeff::C1qg_reg_1[idx])/(ccoeff::C1qg_1[idx]-aass*regfac*ccoeff::C1qg_reg_1[idx]);
	      aexpqg_1[idx] = pow(aexpB,C1);
	      aexpqg_1_bb[idx] = pow(aexpB_bb,C1);
	      aexpqqb_1[idx] = pow(aexpB,ccoeff::C3qqb_1[idx]/ccoeff::C2qqb_1[idx]);
	      aexpqqp_1[idx] = pow(aexpB,ccoeff::C3qqp_1[idx]/ccoeff::C2qqp_1[idx]);
	      aexpqqbp_1[idx] = pow(aexpB,ccoeff::C3qqbp_1[idx]/ccoeff::C2qqbp_1[idx]);
	      aexpqq_2[idx] = pow(aexpB,ccoeff::C1qq_2[idx]);
	      C1 = (ccoeff::C2qg_2[idx]+regfac*ccoeff::C1qg_reg_2[idx])/(ccoeff::C1qg_2[idx]-aass*regfac*ccoeff::C1qg_reg_2[idx]);
	      aexpqg_2[idx] = pow(aexpB,C1);
	      aexpqg_2_bb[idx] = pow(aexpB_bb,C1);
	      aexpqqb_2[idx] = pow(aexpB,ccoeff::C3qqb_2[idx]/ccoeff::C2qqb_2[idx]);
	      aexpqqp_2[idx] = pow(aexpB,ccoeff::C3qqp_2[idx]/ccoeff::C2qqp_2[idx]);
	      aexpqqbp_2[idx] = pow(aexpB,ccoeff::C3qqbp_2[idx]/ccoeff::C2qqbp_2[idx]);
	    }
      else if (opts.order_expc == 3)
	for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	  for (int i = 0; i < mellinint::mdim; i++)
	    {
	      int idx = anomalous::index(i,sign);
	      complex <double> C1,C2;
	      C1 = ccoeff::C1qq_1[idx]; C2 = ccoeff::C2qq_1[idx];
	      aexpqq_1[idx] = exp(
				+lamB*(C1)
				+lamC*((0.5*pow(C1,2) - C2))
				+lamD*(-beta1/beta0*C1)
				+lamC*(beta0*C1*(LQR))
				);
	      C1 = ccoeff::C1qq_2[idx]; C2 = ccoeff::C2qq_2[idx];
	      aexpqq_2[idx] = exp(
				+lamB*(C1)
				+lamC*((0.5*pow(C1,2) - C2))
				+lamD*(-beta1/beta0*C1)
				+lamC*(beta0*C1*(LQR))
				);

	      C1 = (ccoeff::C2qg_1[idx]+regfac*ccoeff::C1qg_reg_1[idx])/(ccoeff::C1qg_1[idx]-aass*regfac*ccoeff::C1qg_reg_1[idx]);
	      C2 = ccoeff::C3qg_1[idx]/(ccoeff::C1qg_1[idx]-aass*regfac*ccoeff::C1qg_reg_1[idx]);
	      aexpqg_1[idx] = exp(
				+lamB*(C1)
				+lamC*((0.5*pow(C1,2) - C2))
				+lamD*(-beta1/beta0*C1)
				+lamC*(beta0*C1*(LQR))
				);
	      C1 = (ccoeff::C2qg_2[idx]+regfac*ccoeff::C1qg_reg_2[idx])/(ccoeff::C1qg_2[idx]-aass*regfac*ccoeff::C1qg_reg_2[idx]);
	      C2 = ccoeff::C3qg_2[idx]/(ccoeff::C1qg_2[idx]-aass*regfac*ccoeff::C1qg_reg_2[idx]);
	      aexpqg_2[idx] = exp(
				+lamB*(C1)
				+lamC*((0.5*pow(C1,2) - C2))
				+lamD*(-beta1/beta0*C1)
				+lamC*(beta0*C1*(LQR))
				);

	      aexpqqb_1[idx] = pow(aexpB,ccoeff::C3qqb_1[idx]/ccoeff::C2qqb_1[idx])
		*pow(aexpC,pow(ccoeff::C3qqb_1[idx]/ccoeff::C2qqb_1[idx],2)/2.)// - C4/C2)
		*pow(aexpD,-beta1/beta0*ccoeff::C3qqb_1[idx]/ccoeff::C2qqb_1[idx])
		*pow(aexpC,beta0*ccoeff::C3qqb_1[idx]/ccoeff::C2qqb_1[idx]*LQR)
		;
	      aexpqqp_1[idx] = pow(aexpB,ccoeff::C3qqp_1[idx]/ccoeff::C2qqp_1[idx])
		*pow(aexpC,pow(ccoeff::C3qqp_1[idx]/ccoeff::C2qqp_1[idx],2)/2.)// - C4/C2))
		*pow(aexpD,-beta1/beta0*ccoeff::C3qqp_1[idx]/ccoeff::C2qqp_1[idx])
		*pow(aexpC,beta0*ccoeff::C3qqp_1[idx]/ccoeff::C2qqp_1[idx]*LQR)
		;
	      aexpqqbp_1[idx] = pow(aexpB,ccoeff::C3qqbp_1[idx]/ccoeff::C2qqbp_1[idx])
		*pow(aexpC,pow(ccoeff::C3qqbp_1[idx]/ccoeff::C2qqbp_1[idx],2)/2.)// - C4/C2))
		*pow(aexpD,-beta1/beta0*ccoeff::C3qqbp_1[idx]/ccoeff::C2qqbp_1[idx])
		*pow(aexpC,beta0*ccoeff::C3qqbp_1[idx]/ccoeff::C2qqbp_1[idx]*LQR)
		;

	      aexpqqb_2[idx] = pow(aexpB,ccoeff::C3qqb_2[idx]/ccoeff::C2qqb_2[idx])
		*pow(aexpC,pow(ccoeff::C3qqb_2[idx]/ccoeff::C2qqb_2[idx],2)/2.)// - C4/C2)
		*pow(aexpD,-beta1/beta0*ccoeff::C3qqb_2[idx]/ccoeff::C2qqb_2[idx])
		*pow(aexpC,beta0*ccoeff::C3qqb_2[idx]/ccoeff::C2qqb_2[idx]*LQR)
		;
	      aexpqqp_2[idx] = pow(aexpB,ccoeff::C3qqp_2[idx]/ccoeff::C2qqp_2[idx])
		*pow(aexpC,pow(ccoeff::C3qqp_2[idx]/ccoeff::C2qqp_2[idx],2)/2.)// - C4/C2))
		*pow(aexpD,-beta1/beta0*ccoeff::C3qqp_2[idx]/ccoeff::C2qqp_2[idx])
		*pow(aexpC,beta0*ccoeff::C3qqp_2[idx]/ccoeff::C2qqp_2[idx]*LQR)
		;
	      aexpqqbp_2[idx] = pow(aexpB,ccoeff::C3qqbp_2[idx]/ccoeff::C2qqbp_2[idx])
		*pow(aexpC,pow(ccoeff::C3qqbp_2[idx]/ccoeff::C2qqbp_2[idx],2)/2.)// - C4/C2))
		*pow(aexpD,-beta1/beta0*ccoeff::C3qqbp_2[idx]/ccoeff::C2qqbp_2[idx])
		*pow(aexpC,beta0*ccoeff::C3qqbp_2[idx]/ccoeff::C2qqbp_2[idx]*LQR)
		;
	      
	    }
      else if (opts.order_expc == 4)
	for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	  for (int i = 0; i < mellinint::mdim; i++)
	    {
	      int idx = anomalous::index(i,sign);

	      complex <double> C1,C2,C3;	      
	      C1 = ccoeff::C1qq_1[idx]; C2 = ccoeff::C2qq_1[idx]; C3 = ccoeff::C3qq_1[idx];
	      aexpqq_1[idx] = exp(
				+lamB*(C1)
				+lamC*((0.5*pow(C1,2) - C2))
				+lamD*(-beta1/beta0*C1)
				+lamC*(beta0*C1*(LQR))
				+lamE*(-pow(beta1/beta0,2)*C1)
				+lamE*(beta2/beta0*C1)
				+lamF*(1./3.*(pow(C1,3)-3.*C1*C2+3.*C3))
				+lamG*(-beta1/pow(beta0,2)*(beta1*C1 - beta0*pow(C1,2) + 2.*beta0*C2))
				+lamH*(pow(beta1/beta0,2)*C1)
				+lamG*(2.*beta1*C1*LQR)
				+lamF*((-beta1*C1+beta0*pow(C1,2)-2.*beta0*C2)*LQR)
				+lamF*(pow(beta0,2)*C1*pow(LQR,2))
				);
	      C1 = ccoeff::C1qq_2[idx]; C2 = ccoeff::C2qq_2[idx]; C3 = ccoeff::C3qq_2[idx];
	      aexpqq_2[idx] = exp(
				+lamB*(C1)
				+lamC*((0.5*pow(C1,2) - C2))
				+lamD*(-beta1/beta0*C1)
				+lamC*(beta0*C1*(LQR))
				+lamE*(-pow(beta1/beta0,2)*C1)
				+lamE*(beta2/beta0*C1)
				+lamF*(1./3.*(pow(C1,3)-3.*C1*C2+3.*C3))
				+lamG*(-beta1/pow(beta0,2)*(beta1*C1 - beta0*pow(C1,2) + 2.*beta0*C2))
				+lamH*(pow(beta1/beta0,2)*C1)
				+lamG*(2.*beta1*C1*LQR)
				+lamF*((-beta1*C1+beta0*pow(C1,2)-2.*beta0*C2)*LQR)
				+lamF*(pow(beta0,2)*C1*pow(LQR,2))
				);

	      C1 = (ccoeff::C2qg_1[idx]+regfac*ccoeff::C1qg_reg_1[idx])/(ccoeff::C1qg_1[idx]-aass*regfac*ccoeff::C1qg_reg_1[idx]);
	      C2 = ccoeff::C3qg_1[idx]/(ccoeff::C1qg_1[idx]-aass*regfac*ccoeff::C1qg_reg_1[idx]);
	      C3 = ccoeff::C4qg_1[idx]/(ccoeff::C1qg_1[idx]-aass*regfac*ccoeff::C1qg_reg_1[idx]);
	      aexpqg_1[idx] = exp(
				+lamB*(C1)
				+lamC*((0.5*pow(C1,2) - C2))
				+lamD*(-beta1/beta0*C1)
				+lamC*(beta0*C1*(LQR))
				+lamE*(-pow(beta1/beta0,2)*C1)
				+lamE*(beta2/beta0*C1)
				+lamF*(1./3.*(pow(C1,3)-3.*C1*C2+3.*C3))
				+lamG*(-beta1/pow(beta0,2)*(beta1*C1 - beta0*pow(C1,2) + 2.*beta0*C2))
				+lamH*(pow(beta1/beta0,2)*C1)
				+lamG*(2.*beta1*C1*LQR)
				+lamF*((-beta1*C1+beta0*pow(C1,2)-2.*beta0*C2)*LQR)
				+lamF*(pow(beta0,2)*C1*pow(LQR,2))
				);
	      C1 = (ccoeff::C2qg_2[idx]+regfac*ccoeff::C1qg_reg_2[idx])/(ccoeff::C1qg_2[idx]-aass*regfac*ccoeff::C1qg_reg_2[idx]);
	      C2 = ccoeff::C3qg_2[idx]/(ccoeff::C1qg_2[idx]-aass*regfac*ccoeff::C1qg_reg_2[idx]);
	      C3 = ccoeff::C4qg_2[idx]/(ccoeff::C1qg_2[idx]-aass*regfac*ccoeff::C1qg_reg_2[idx]);
	      aexpqg_2[idx] = exp(
				+lamB*(C1)
				+lamC*((0.5*pow(C1,2) - C2))
				+lamD*(-beta1/beta0*C1)
				+lamC*(beta0*C1*(LQR))
				+lamE*(-pow(beta1/beta0,2)*C1)
				+lamE*(beta2/beta0*C1)
				+lamF*(1./3.*(pow(C1,3)-3.*C1*C2+3.*C3))
				+lamG*(-beta1/pow(beta0,2)*(beta1*C1 - beta0*pow(C1,2) + 2.*beta0*C2))
				+lamH*(pow(beta1/beta0,2)*C1)
				+lamG*(2.*beta1*C1*LQR)
				+lamF*((-beta1*C1+beta0*pow(C1,2)-2.*beta0*C2)*LQR)
				+lamF*(pow(beta0,2)*C1*pow(LQR,2))
				);
	      
	      //C1 = ccoeff::C3qqb_1[idx]/ccoeff::C2qqb_1[idx]; C2 = ccoeff::C4qqb_1[idx]/ccoeff::C2qqb_1[idx]; C3 = 0.;
	      C1 = ccoeff::C3qqb_1[idx]/ccoeff::C2qqb_1[idx]; C2 = 0.; C3 = 0.;
	      aexpqqb_1[idx] = exp(
				 +lamB*(C1)
				 +lamC*((0.5*pow(C1,2) - C2))
				 +lamD*(-beta1/beta0*C1)
				 +lamC*(beta0*C1*(LQR))
				 +lamE*(-pow(beta1/beta0,2)*C1)
				 +lamE*(beta2/beta0*C1)
				 +lamF*(1./3.*(pow(C1,3)-3.*C1*C2+3.*C3))
				 +lamG*(-beta1/pow(beta0,2)*(beta1*C1 - beta0*pow(C1,2) + 2.*beta0*C2))
				 +lamH*(pow(beta1/beta0,2)*C1)
				 +lamG*(2.*beta1*C1*LQR)
				 +lamF*((-beta1*C1+beta0*pow(C1,2)-2.*beta0*C2)*LQR)
				 +lamF*(pow(beta0,2)*C1*pow(LQR,2))
				 );
	      //C1 = ccoeff::C3qqp_1[idx]/ccoeff::C2qqp_1[idx]; C2 = ccoeff::C4qqp_1[idx]/ccoeff::C2qqp_1[idx]; C3 = 0.;
	      C1 = ccoeff::C3qqp_1[idx]/ccoeff::C2qqp_1[idx]; C2 = 0.; C3 = 0.;
	      aexpqqp_1[idx] = exp(
				 +lamB*(C1)
				 +lamC*((0.5*pow(C1,2) - C2))
				 +lamD*(-beta1/beta0*C1)
				 +lamC*(beta0*C1*(LQR))
				 +lamE*(-pow(beta1/beta0,2)*C1)
				 +lamE*(beta2/beta0*C1)
				 +lamF*(1./3.*(pow(C1,3)-3.*C1*C2+3.*C3))
				 +lamG*(-beta1/pow(beta0,2)*(beta1*C1 - beta0*pow(C1,2) + 2.*beta0*C2))
				 +lamH*(pow(beta1/beta0,2)*C1)
				 +lamG*(2.*beta1*C1*LQR)
				 +lamF*((-beta1*C1+beta0*pow(C1,2)-2.*beta0*C2)*LQR)
				 +lamF*(pow(beta0,2)*C1*pow(LQR,2))
				 );
	      //C1 = ccoeff::C3qqbp_1[idx]/ccoeff::C2qqbp_1[idx]; C2 = ccoeff::C4qqbp_1[idx]/ccoeff::C2qqbp_1[idx]; C3 = 0.;
	      C1 = ccoeff::C3qqbp_1[idx]/ccoeff::C2qqbp_1[idx]; C2 = 0.; C3 = 0.;
	      aexpqqbp_1[idx] = exp(
				  +lamB*(C1)
				  +lamC*((0.5*pow(C1,2) - C2))
				  +lamD*(-beta1/beta0*C1)
				  +lamC*(beta0*C1*(LQR))
				  +lamE*(-pow(beta1/beta0,2)*C1)
				  +lamE*(beta2/beta0*C1)
				  +lamF*(1./3.*(pow(C1,3)-3.*C1*C2+3.*C3))
				  +lamG*(-beta1/pow(beta0,2)*(beta1*C1 - beta0*pow(C1,2) + 2.*beta0*C2))
				  +lamH*(pow(beta1/beta0,2)*C1)
				  +lamG*(2.*beta1*C1*LQR)
				  +lamF*((-beta1*C1+beta0*pow(C1,2)-2.*beta0*C2)*LQR)
				  +lamF*(pow(beta0,2)*C1*pow(LQR,2))
				 );

	      //C1 = ccoeff::C3qqb_2[idx]/ccoeff::C2qqb_2[idx]; C2 = ccoeff::C4qqb_2[idx]/ccoeff::C2qqb_2[idx]; C3 = 0.;
	      C1 = ccoeff::C3qqb_2[idx]/ccoeff::C2qqb_2[idx]; C2 = 0.; C3 = 0.;
	      aexpqqb_2[idx] = exp(
				 +lamB*(C1)
				 +lamC*((0.5*pow(C1,2) - C2))
				 +lamD*(-beta1/beta0*C1)
				 +lamC*(beta0*C1*(LQR))
				 +lamE*(-pow(beta1/beta0,2)*C1)
				 +lamE*(beta2/beta0*C1)
				 +lamF*(1./3.*(pow(C1,3)-3.*C1*C2+3.*C3))
				 +lamG*(-beta1/pow(beta0,2)*(beta1*C1 - beta0*pow(C1,2) + 2.*beta0*C2))
				 +lamH*(pow(beta1/beta0,2)*C1)
				 +lamG*(2.*beta1*C1*LQR)
				 +lamF*((-beta1*C1+beta0*pow(C1,2)-2.*beta0*C2)*LQR)
				 +lamF*(pow(beta0,2)*C1*pow(LQR,2))
				 );
	      //C1 = ccoeff::C3qqp_2[idx]/ccoeff::C2qqp_2[idx]; C2 = ccoeff::C4qqp_2[idx]/ccoeff::C2qqp_2[idx]; C3 = 0.;
	      C1 = ccoeff::C3qqp_2[idx]/ccoeff::C2qqp_2[idx]; C2 = 0.; C3 = 0.;
	      aexpqqp_2[idx] = exp(
				 +lamB*(C1)
				 +lamC*((0.5*pow(C1,2) - C2))
				 +lamD*(-beta1/beta0*C1)
				 +lamC*(beta0*C1*(LQR))
				 +lamE*(-pow(beta1/beta0,2)*C1)
				 +lamE*(beta2/beta0*C1)
				 +lamF*(1./3.*(pow(C1,3)-3.*C1*C2+3.*C3))
				 +lamG*(-beta1/pow(beta0,2)*(beta1*C1 - beta0*pow(C1,2) + 2.*beta0*C2))
				 +lamH*(pow(beta1/beta0,2)*C1)
				 +lamG*(2.*beta1*C1*LQR)
				 +lamF*((-beta1*C1+beta0*pow(C1,2)-2.*beta0*C2)*LQR)
				 +lamF*(pow(beta0,2)*C1*pow(LQR,2))
				 );
	      //C1 = ccoeff::C3qqbp_2[idx]/ccoeff::C2qqbp_2[idx]; C2 = ccoeff::C4qqbp_2[idx]/ccoeff::C2qqbp_2[idx]; C3 = 0.;
	      C1 = ccoeff::C3qqbp_2[idx]/ccoeff::C2qqbp_2[idx]; C2 = 0.; C3 = 0.;
	      aexpqqbp_2[idx] = exp(
				  +lamB*(C1)
				  +lamC*((0.5*pow(C1,2) - C2))
				  +lamD*(-beta1/beta0*C1)
				  +lamC*(beta0*C1*(LQR))
				  +lamE*(-pow(beta1/beta0,2)*C1)
				  +lamE*(beta2/beta0*C1)
				  +lamF*(1./3.*(pow(C1,3)-3.*C1*C2+3.*C3))
				  +lamG*(-beta1/pow(beta0,2)*(beta1*C1 - beta0*pow(C1,2) + 2.*beta0*C2))
				  +lamH*(pow(beta1/beta0,2)*C1)
				  +lamG*(2.*beta1*C1*LQR)
				  +lamF*((-beta1*C1+beta0*pow(C1,2)-2.*beta0*C2)*LQR)
				  +lamF*(pow(beta0,2)*C1*pow(LQR,2))
				 );
	    }
    }
}

void expc::denlow()
{
  if (opts.order_expc == 1)
    return;
  
  double LQR = LR-LQ;
  if (opts.mellin1d)
    {
      if (opts.order_expc == 2)
	for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	  for (int i = 0; i < mellinint::mdim; i++)
	    {
	      int idx = anomalous::index(i,sign);
	      aexpqq[idx] = pow(aexpB,ccoeff::C1qq[idx]);
	    }
      else if (opts.order_expc == 3)
	for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	  for (int i = 0; i < mellinint::mdim; i++)
	    {
	      int idx = anomalous::index(i,sign);
	      aexpqq[idx] = pow(aexpB,ccoeff::C1qq[idx])
		*pow(aexpC,(0.5*pow(ccoeff::C1qq[idx],2) - ccoeff::C2qq[idx]))
		*pow(aexpD,-beta1/beta0*ccoeff::C1qq[idx])
		*pow(aexpC,beta0*ccoeff::C1qq[idx]*(LQR))
		;
	      aexpqg[idx]   = (1.-xlambda/(1.+aass*ccoeff::C2qg[idx]/ccoeff::C1qg[idx]))/(1.-xlambda);
	    }
      
    }
  else
    {
      if (opts.order_expc == 2)
	for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	  for (int i = 0; i < mellinint::mdim; i++)
	    {
	      int idx = anomalous::index(i,sign);
	      aexpqq_1[idx] = pow(aexpB,ccoeff::C1qq_1[idx]);
	      aexpqq_2[idx] = pow(aexpB,ccoeff::C1qq_2[idx]);
	    }
      else if (opts.order_expc == 3)
	for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	  for (int i = 0; i < mellinint::mdim; i++)
	    {
	      int idx = anomalous::index(i,sign);
	      aexpqq_1[idx] = pow(aexpB,ccoeff::C1qq_1[idx])
		*pow(aexpC,(0.5*pow(ccoeff::C1qq_1[idx],2) - ccoeff::C2qq_1[idx]))
		*pow(aexpD,-beta1/beta0*ccoeff::C1qq_1[idx])
		*pow(aexpC,beta0*ccoeff::C1qq_1[idx]*(LQR))
		;
	      aexpqq_2[idx] = pow(aexpB,ccoeff::C1qq_2[idx])
		*pow(aexpC,(0.5*pow(ccoeff::C1qq_2[idx],2) - ccoeff::C2qq_2[idx]))
		*pow(aexpD,-beta1/beta0*ccoeff::C1qq_2[idx])
		*pow(aexpC,beta0*ccoeff::C1qq_2[idx]*(LQR))
		;
	      aexpqg_1[idx]   = (1.-xlambda/(1.+aass*ccoeff::C2qg_1[idx]/ccoeff::C1qg_1[idx]))/(1.-xlambda);
	      aexpqg_2[idx]   = (1.-xlambda/(1.+aass*ccoeff::C2qg_2[idx]/ccoeff::C1qg_2[idx]))/(1.-xlambda);
	    }
    }
}

void expc::denfull()
{
  if (opts.order_expc == 1)
    return;
  
  double LQR = LR-LQ;
  double aass2 = pow(aass,2);
  if (opts.mellin1d)
    {
      if (opts.order_expc == 2)
	for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	  for (int i = 0; i < mellinint::mdim; i++)
	    {
	      int idx = anomalous::index(i,sign);
	      aexpqq[idx]   = pow(aexpB,ccoeff::C1qq[idx]);
	      aexpqg[idx]   = (1.-xlambda/(1.+aass*ccoeff::C2qg[idx]/ccoeff::C1qg[idx]))/(1.-xlambda);
	      aexpqqb[idx]  = (1.-xlambda/(1.+aass*ccoeff::C3qqb[idx] /ccoeff::C2qqb[idx] ))/(1.-xlambda);
	      aexpqqp[idx]  = (1.-xlambda/(1.+aass*ccoeff::C3qqp[idx] /ccoeff::C2qqp[idx] ))/(1.-xlambda);
	      aexpqqbp[idx] = (1.-xlambda/(1.+aass*ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx]))/(1.-xlambda);
	    }
      else if (opts.order_expc == 3)
	for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	  for (int i = 0; i < mellinint::mdim; i++)
	    {
	      int idx = anomalous::index(i,sign);
	      aexpqq[idx] = pow(aexpB,ccoeff::C1qq[idx])
		*pow(aexpC,(0.5*pow(ccoeff::C1qq[idx],2) - ccoeff::C2qq[idx]))
		*pow(aexpD,-beta1/beta0*ccoeff::C1qq[idx])
		*pow(aexpC,beta0*ccoeff::C1qq[idx]*(LQR))
		;

	      complex <double> C1,C2,C3,rootC,n2llaexp;
	      complex <double> lam = xlambda;

	      //cdefexp(ccoeff::C1qg[idx],   ccoeff::C2qg[idx],   ccoeff::C3qg[idx]); aexpqg[idx]   = n3llaexp();
	      cdefexp(ccoeff::C1qg[idx],   ccoeff::C2qg[idx],   0.               ); aexpqg[idx]   = n3llaexp();
	      cdefexp(ccoeff::C2qqb[idx],  ccoeff::C3qqb[idx],  0.               ); aexpqqb[idx]  = n3llaexp();
	      cdefexp(ccoeff::C2qqp[idx],  ccoeff::C3qqp[idx],  0.               ); aexpqqp[idx]  = n3llaexp();
	      cdefexp(ccoeff::C2qqbp[idx], ccoeff::C3qqbp[idx], 0.               ); aexpqqbp[idx] = n3llaexp();
	    }
      
    }
  else
    {
      if (opts.order_expc == 2)
	for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	  for (int i = 0; i < mellinint::mdim; i++)
	    {
	      int idx = anomalous::index(i,sign);
	      aexpqq_1[idx] = pow(aexpB,ccoeff::C1qq_1[idx]);
	      aexpqg_1[idx]   = (1. - xlambda/(1.+aass*ccoeff::C2qg_1[idx]/ccoeff::C1qg_1[idx]))/(1.-xlambda);
	      aexpqqb_1[idx]  = (1. - xlambda/(1.+aass*ccoeff::C3qqb_1[idx] /ccoeff::C2qqb_1[idx] ))/(1.-xlambda);
	      aexpqqp_1[idx]  = (1. - xlambda/(1.+aass*ccoeff::C3qqp_1[idx] /ccoeff::C2qqp_1[idx] ))/(1.-xlambda);
	      aexpqqbp_1[idx] = (1. - xlambda/(1.+aass*ccoeff::C3qqbp_1[idx]/ccoeff::C2qqbp_1[idx]))/(1.-xlambda);
	      aexpqq_2[idx] = pow(aexpB,ccoeff::C1qq_2[idx]);
	      aexpqg_2[idx]   = (1. - xlambda/(1.+aass*ccoeff::C2qg_2[idx]/ccoeff::C1qg_2[idx]))/(1.-xlambda);
	      aexpqqb_2[idx]  = (1. - xlambda/(1.+aass*ccoeff::C3qqb_2[idx] /ccoeff::C2qqb_2[idx] ))/(1.-xlambda);
	      aexpqqp_2[idx]  = (1. - xlambda/(1.+aass*ccoeff::C3qqp_2[idx] /ccoeff::C2qqp_2[idx] ))/(1.-xlambda);
	      aexpqqbp_2[idx] = (1. - xlambda/(1.+aass*ccoeff::C3qqbp_2[idx]/ccoeff::C2qqbp_2[idx]))/(1.-xlambda);
	    }
      else if (opts.order_expc == 3)
	for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	  for (int i = 0; i < mellinint::mdim; i++)
	    {
	      int idx = anomalous::index(i,sign);
	      aexpqq_1[idx] = pow(aexpB,ccoeff::C1qq_1[idx])
		*pow(aexpC,(0.5*pow(ccoeff::C1qq_1[idx],2) - ccoeff::C2qq_1[idx]))
		*pow(aexpD,-beta1/beta0*ccoeff::C1qq_1[idx])
		*pow(aexpC,beta0*ccoeff::C1qq_1[idx]*(LQR))
		;
	      aexpqq_2[idx] = pow(aexpB,ccoeff::C1qq_2[idx])
		*pow(aexpC,(0.5*pow(ccoeff::C1qq_2[idx],2) - ccoeff::C2qq_2[idx]))
		*pow(aexpD,-beta1/beta0*ccoeff::C1qq_2[idx])
		*pow(aexpC,beta0*ccoeff::C1qq_2[idx]*(LQR))
		;

	      complex <double> C1,C2,C3,rootC,n2llaexp;
	      complex <double> lam = xlambda;
	      
	      //cdefexp(ccoeff::C1qg_1[idx],   ccoeff::C2qg_1[idx],   ccoeff::C3qg_1[idx]); aexpqg_1[idx]   = n3llaexp();
	      cdefexp(ccoeff::C1qg_1[idx],   ccoeff::C2qg_1[idx],   0.               ); aexpqg_1[idx]   = n3llaexp();
	      cdefexp(ccoeff::C2qqb_1[idx],  ccoeff::C3qqb_1[idx],  0.               ); aexpqqb_1[idx]  = n3llaexp();
	      cdefexp(ccoeff::C2qqp_1[idx],  ccoeff::C3qqp_1[idx],  0.               ); aexpqqp_1[idx]  = n3llaexp();
	      cdefexp(ccoeff::C2qqbp_1[idx], ccoeff::C3qqbp_1[idx], 0.               ); aexpqqbp_1[idx] = n3llaexp();

	      //cdefexp(ccoeff::C1qg_2[idx],   ccoeff::C2qg_2[idx],   ccoeff::C3qg_2[idx]); aexpqg_2[idx]   = n3llaexp();
	      cdefexp(ccoeff::C1qg_2[idx],   ccoeff::C2qg_2[idx],   0.               ); aexpqg_2[idx]   = n3llaexp();
	      cdefexp(ccoeff::C2qqb_2[idx],  ccoeff::C3qqb_2[idx],  0.               ); aexpqqb_2[idx]  = n3llaexp();
	      cdefexp(ccoeff::C2qqp_2[idx],  ccoeff::C3qqp_2[idx],  0.               ); aexpqqp_2[idx]  = n3llaexp();
	      cdefexp(ccoeff::C2qqbp_2[idx], ccoeff::C3qqbp_2[idx], 0.               ); aexpqqbp_2[idx] = n3llaexp();
	    }
    }
}

void expc::matrix()
{
  if (opts.order_expc == 1)
    return;
  
  double LQR = LR-LQ;
  if (opts.mellin1d)
    {
      if (opts.order_expc == 2)
	for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	  for (int i = 0; i < mellinint::mdim; i++)
	    {
	      int idx = anomalous::index(i,sign);
	      aexpqq[idx] = pow(aexpB,ccoeff::C1qq[idx]);
	      aexpqg[idx] = pow(aexpB,ccoeff::C1qg[idx]);
	    }
      else if (opts.order_expc == 3)
	for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	  for (int i = 0; i < mellinint::mdim; i++)
	    {
	      int idx = anomalous::index(i,sign);
	      complex <double> C1,C2;
	      C1 = ccoeff::C1qq[idx]; C2 = ccoeff::C2qq[idx];
	      aexpqq[idx] = exp(
				+lamB*(C1)
				+lamC*((0.5*pow(C1,2) - C2))
				+lamD*(-beta1/beta0*C1)
				+lamC*(beta0*C1*(LQR))
				);
	      C1 = ccoeff::C1qg[idx]; C2 = ccoeff::C2qg[idx];
	      aexpqg[idx] = exp(+lamB*(C1)
				+lamC*((0.5*pow(C1,2) - C2))
				+lamD*(-beta1/beta0*C1)
				+lamC*(beta0*C1*(LQR)));
	      C1 = 0.; C2 = ccoeff::C2qqb[idx];
	      aexpqqb[idx] = exp(+lamB*(C1)
				+lamC*((0.5*pow(C1,2) - C2))
				+lamD*(-beta1/beta0*C1)
				+lamC*(beta0*C1*(LQR)));
	      C1 = 0.; C2 = ccoeff::C2qqp[idx];
	      aexpqqp[idx] = exp(+lamB*(C1)
				+lamC*((0.5*pow(C1,2) - C2))
				+lamD*(-beta1/beta0*C1)
				+lamC*(beta0*C1*(LQR)));
	      C1 = 0.; C2 = ccoeff::C2qqbp[idx];
	      aexpqqbp[idx] = exp(+lamB*(C1)
				+lamC*((0.5*pow(C1,2) - C2))
				+lamD*(-beta1/beta0*C1)
				+lamC*(beta0*C1*(LQR)));
	    }
      else if (opts.order_expc == 4)
	for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	  for (int i = 0; i < mellinint::mdim; i++)
	    {
	      int idx = anomalous::index(i,sign);

	      complex <double> C1,C2,C3;	      
	      
	      C1 = ccoeff::C1qq[idx]; C2 = ccoeff::C2qq[idx]; C3 = ccoeff::C3qq[idx];
	      aexpqq[idx] = exp(
				+lamB*(C1)
				+lamC*((0.5*pow(C1,2) - C2))
				+lamD*(-beta1/beta0*C1)
				+lamC*(beta0*C1*(LQR))
				+lamE*(-pow(beta1/beta0,2)*C1)
				+lamE*(beta2/beta0*C1)
				+lamF*(1./3.*(pow(C1,3)-3.*C1*C2+3.*C3))
				+lamG*(-beta1/pow(beta0,2)*(beta1*C1 - beta0*pow(C1,2) + 2.*beta0*C2))
				+lamH*(pow(beta1/beta0,2)*C1)
				+lamG*(2.*beta1*C1*LQR)
				+lamF*((-beta1*C1+beta0*pow(C1,2)-2.*beta0*C2)*LQR)
				+lamF*(pow(beta0,2)*C1*pow(LQR,2))
				);
	      C1 = ccoeff::C1qg[idx]; C2 = ccoeff::C2qg[idx]; C3 = ccoeff::C3qg[idx];
	      aexpqg[idx] = exp(
				+lamB*(C1)
				+lamC*((0.5*pow(C1,2) - C2))
				+lamD*(-beta1/beta0*C1)
				+lamC*(beta0*C1*(LQR))
				+lamE*(-pow(beta1/beta0,2)*C1)
				+lamE*(beta2/beta0*C1)
				+lamF*(1./3.*(pow(C1,3)-3.*C1*C2+3.*C3))
				+lamG*(-beta1/pow(beta0,2)*(beta1*C1 - beta0*pow(C1,2) + 2.*beta0*C2))
				+lamH*(pow(beta1/beta0,2)*C1)
				+lamG*(2.*beta1*C1*LQR)
				+lamF*((-beta1*C1+beta0*pow(C1,2)-2.*beta0*C2)*LQR)
				+lamF*(pow(beta0,2)*C1*pow(LQR,2))
				);
	      C1 = 0.; C2 = ccoeff::C2qqb[idx]; C3 = ccoeff::C3qqb[idx];
	      aexpqqb[idx] = exp(
				 +lamB*(C1)
				 +lamC*((0.5*pow(C1,2) - C2))
				 +lamD*(-beta1/beta0*C1)
				 +lamC*(beta0*C1*(LQR))
				 +lamE*(-pow(beta1/beta0,2)*C1)
				 +lamE*(beta2/beta0*C1)
				 +lamF*(1./3.*(pow(C1,3)-3.*C1*C2+3.*C3))
				 +lamG*(-beta1/pow(beta0,2)*(beta1*C1 - beta0*pow(C1,2) + 2.*beta0*C2))
				 +lamH*(pow(beta1/beta0,2)*C1)
				 +lamG*(2.*beta1*C1*LQR)
				 +lamF*((-beta1*C1+beta0*pow(C1,2)-2.*beta0*C2)*LQR)
				 +lamF*(pow(beta0,2)*C1*pow(LQR,2))
				 );
		
	      C1 = 0.; C2 = ccoeff::C2qqp[idx]; C3 = ccoeff::C3qqp[idx];
	      aexpqqp[idx] = exp(
				 +lamB*(C1)
				 +lamC*((0.5*pow(C1,2) - C2))
				 +lamD*(-beta1/beta0*C1)
				 +lamC*(beta0*C1*(LQR))
				 +lamE*(-pow(beta1/beta0,2)*C1)
				 +lamE*(beta2/beta0*C1)
				 +lamF*(1./3.*(pow(C1,3)-3.*C1*C2+3.*C3))
				 +lamG*(-beta1/pow(beta0,2)*(beta1*C1 - beta0*pow(C1,2) + 2.*beta0*C2))
				 +lamH*(pow(beta1/beta0,2)*C1)
				 +lamG*(2.*beta1*C1*LQR)
				 +lamF*((-beta1*C1+beta0*pow(C1,2)-2.*beta0*C2)*LQR)
				 +lamF*(pow(beta0,2)*C1*pow(LQR,2))
				 );
		
	      C1 = 0.; C2 = ccoeff::C2qqbp[idx]; C3 = ccoeff::C3qqbp[idx];
	      aexpqqbp[idx] = exp(
				  +lamB*(C1)
				  +lamC*((0.5*pow(C1,2) - C2))
				  +lamD*(-beta1/beta0*C1)
				  +lamC*(beta0*C1*(LQR))
				  +lamE*(-pow(beta1/beta0,2)*C1)
				  +lamE*(beta2/beta0*C1)
				  +lamF*(1./3.*(pow(C1,3)-3.*C1*C2+3.*C3))
				  +lamG*(-beta1/pow(beta0,2)*(beta1*C1 - beta0*pow(C1,2) + 2.*beta0*C2))
				  +lamH*(pow(beta1/beta0,2)*C1)
				  +lamG*(2.*beta1*C1*LQR)
				  +lamF*((-beta1*C1+beta0*pow(C1,2)-2.*beta0*C2)*LQR)
				  +lamF*(pow(beta0,2)*C1*pow(LQR,2))
				 );
	    }
    }
  else
    {
      if (opts.order_expc == 2)
	for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	  for (int i = 0; i < mellinint::mdim; i++)
	    {
	      int idx = anomalous::index(i,sign);
	      aexpqq_1[idx] = pow(aexpB,ccoeff::C1qq_1[idx]);
	      aexpqg_1[idx] = pow(aexpB,ccoeff::C1qg_1[idx]);
	      aexpqg_2[idx] = pow(aexpB,ccoeff::C1qg_2[idx]);
	    }
      else if (opts.order_expc == 3)
	for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	  for (int i = 0; i < mellinint::mdim; i++)
	    {
	      int idx = anomalous::index(i,sign);
	      complex <double> C1,C2;
	      C1 = ccoeff::C1qq_1[idx]; C2 = ccoeff::C2qq_1[idx];
	      aexpqq_1[idx] = exp(
				+lamB*(C1)
				+lamC*((0.5*pow(C1,2) - C2))
				+lamD*(-beta1/beta0*C1)
				+lamC*(beta0*C1*(LQR))
				);
	      C1 = ccoeff::C1qq_2[idx]; C2 = ccoeff::C2qq_2[idx];
	      aexpqq_2[idx] = exp(
				+lamB*(C1)
				+lamC*((0.5*pow(C1,2) - C2))
				+lamD*(-beta1/beta0*C1)
				+lamC*(beta0*C1*(LQR))
				);

	      C1 = ccoeff::C1qg_1[idx]; C2 = ccoeff::C2qg_1[idx];
	      aexpqg_1[idx] = exp(
				+lamB*(C1)
				+lamC*((0.5*pow(C1,2) - C2))
				+lamD*(-beta1/beta0*C1)
				+lamC*(beta0*C1*(LQR))
				);
	      C1 = ccoeff::C1qg_2[idx]; C2 = ccoeff::C2qg_2[idx];
	      aexpqg_2[idx] = exp(
				+lamB*(C1)
				+lamC*((0.5*pow(C1,2) - C2))
				+lamD*(-beta1/beta0*C1)
				+lamC*(beta0*C1*(LQR))
				);

	      aexpqqb_1[idx]  = exp(-lamC*ccoeff::C2qqb_1[idx]);
	      aexpqqp_1[idx]  = exp(-lamC*ccoeff::C2qqp_1[idx]);
	      aexpqqbp_1[idx] = exp(-lamC*ccoeff::C2qqbp_1[idx]);
	      aexpqqb_2[idx]  = exp(-lamC*ccoeff::C2qqb_2[idx]); 
	      aexpqqp_2[idx]  = exp(-lamC*ccoeff::C2qqp_2[idx]); 
	      aexpqqbp_2[idx] = exp(-lamC*ccoeff::C2qqbp_2[idx]);
	    }
      else if (opts.order_expc == 4)
	for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	  for (int i = 0; i < mellinint::mdim; i++)
	    {
	      int idx = anomalous::index(i,sign);

	      complex <double> C1,C2,C3;	      
	      C1 = ccoeff::C1qq_1[idx]; C2 = ccoeff::C2qq_1[idx]; C3 = ccoeff::C3qq_1[idx];
	      aexpqq_1[idx] = exp(
				+lamB*(C1)
				+lamC*((0.5*pow(C1,2) - C2))
				+lamD*(-beta1/beta0*C1)
				+lamC*(beta0*C1*(LQR))
				+lamE*(-pow(beta1/beta0,2)*C1)
				+lamE*(beta2/beta0*C1)
				+lamF*(1./3.*(pow(C1,3)-3.*C1*C2+3.*C3))
				+lamG*(-beta1/pow(beta0,2)*(beta1*C1 - beta0*pow(C1,2) + 2.*beta0*C2))
				+lamH*(pow(beta1/beta0,2)*C1)
				+lamG*(2.*beta1*C1*LQR)
				+lamF*((-beta1*C1+beta0*pow(C1,2)-2.*beta0*C2)*LQR)
				+lamF*(pow(beta0,2)*C1*pow(LQR,2))
				);
	      C1 = ccoeff::C1qq_2[idx]; C2 = ccoeff::C2qq_2[idx]; C3 = ccoeff::C3qq_2[idx];
	      aexpqq_2[idx] = exp(
				+lamB*(C1)
				+lamC*((0.5*pow(C1,2) - C2))
				+lamD*(-beta1/beta0*C1)
				+lamC*(beta0*C1*(LQR))
				+lamE*(-pow(beta1/beta0,2)*C1)
				+lamE*(beta2/beta0*C1)
				+lamF*(1./3.*(pow(C1,3)-3.*C1*C2+3.*C3))
				+lamG*(-beta1/pow(beta0,2)*(beta1*C1 - beta0*pow(C1,2) + 2.*beta0*C2))
				+lamH*(pow(beta1/beta0,2)*C1)
				+lamG*(2.*beta1*C1*LQR)
				+lamF*((-beta1*C1+beta0*pow(C1,2)-2.*beta0*C2)*LQR)
				+lamF*(pow(beta0,2)*C1*pow(LQR,2))
				);

	      C1 = ccoeff::C1qg_1[idx]; C2 = ccoeff::C2qg_1[idx]; C3 = ccoeff::C3qg_1[idx];
	      aexpqg_1[idx] = exp(
				+lamB*(C1)
				+lamC*((0.5*pow(C1,2) - C2))
				+lamD*(-beta1/beta0*C1)
				+lamC*(beta0*C1*(LQR))
				+lamE*(-pow(beta1/beta0,2)*C1)
				+lamE*(beta2/beta0*C1)
				+lamF*(1./3.*(pow(C1,3)-3.*C1*C2+3.*C3))
				+lamG*(-beta1/pow(beta0,2)*(beta1*C1 - beta0*pow(C1,2) + 2.*beta0*C2))
				+lamH*(pow(beta1/beta0,2)*C1)
				+lamG*(2.*beta1*C1*LQR)
				+lamF*((-beta1*C1+beta0*pow(C1,2)-2.*beta0*C2)*LQR)
				+lamF*(pow(beta0,2)*C1*pow(LQR,2))
				);
	      C1 = ccoeff::C1qg_2[idx]; C2 = ccoeff::C2qg_2[idx]; C3 = ccoeff::C3qg_2[idx];
	      aexpqg_2[idx] = exp(
				+lamB*(C1)
				+lamC*((0.5*pow(C1,2) - C2))
				+lamD*(-beta1/beta0*C1)
				+lamC*(beta0*C1*(LQR))
				+lamE*(-pow(beta1/beta0,2)*C1)
				+lamE*(beta2/beta0*C1)
				+lamF*(1./3.*(pow(C1,3)-3.*C1*C2+3.*C3))
				+lamG*(-beta1/pow(beta0,2)*(beta1*C1 - beta0*pow(C1,2) + 2.*beta0*C2))
				+lamH*(pow(beta1/beta0,2)*C1)
				+lamG*(2.*beta1*C1*LQR)
				+lamF*((-beta1*C1+beta0*pow(C1,2)-2.*beta0*C2)*LQR)
				+lamF*(pow(beta0,2)*C1*pow(LQR,2))
				);
	      
	      C1 = 0.; C2 = ccoeff::C2qqb_1[idx]; C3 = ccoeff::C3qqb_1[idx];
	      aexpqqb_1[idx] = exp(
				 +lamB*(C1)
				 +lamC*((0.5*pow(C1,2) - C2))
				 +lamD*(-beta1/beta0*C1)
				 +lamC*(beta0*C1*(LQR))
				 +lamE*(-pow(beta1/beta0,2)*C1)
				 +lamE*(beta2/beta0*C1)
				 +lamF*(1./3.*(pow(C1,3)-3.*C1*C2+3.*C3))
				 +lamG*(-beta1/pow(beta0,2)*(beta1*C1 - beta0*pow(C1,2) + 2.*beta0*C2))
				 +lamH*(pow(beta1/beta0,2)*C1)
				 +lamG*(2.*beta1*C1*LQR)
				 +lamF*((-beta1*C1+beta0*pow(C1,2)-2.*beta0*C2)*LQR)
				 +lamF*(pow(beta0,2)*C1*pow(LQR,2))
				 );
	      C1 = 0.; C2 = ccoeff::C2qqp_1[idx]; C3 = ccoeff::C3qqp_1[idx];
	      aexpqqp_1[idx] = exp(
				 +lamB*(C1)
				 +lamC*((0.5*pow(C1,2) - C2))
				 +lamD*(-beta1/beta0*C1)
				 +lamC*(beta0*C1*(LQR))
				 +lamE*(-pow(beta1/beta0,2)*C1)
				 +lamE*(beta2/beta0*C1)
				 +lamF*(1./3.*(pow(C1,3)-3.*C1*C2+3.*C3))
				 +lamG*(-beta1/pow(beta0,2)*(beta1*C1 - beta0*pow(C1,2) + 2.*beta0*C2))
				 +lamH*(pow(beta1/beta0,2)*C1)
				 +lamG*(2.*beta1*C1*LQR)
				 +lamF*((-beta1*C1+beta0*pow(C1,2)-2.*beta0*C2)*LQR)
				 +lamF*(pow(beta0,2)*C1*pow(LQR,2))
				 );
	      C1 = 0.; C2 = ccoeff::C2qqbp_1[idx]; C3 = ccoeff::C3qqbp_1[idx];
	      aexpqqbp_1[idx] = exp(
				  +lamB*(C1)
				  +lamC*((0.5*pow(C1,2) - C2))
				  +lamD*(-beta1/beta0*C1)
				  +lamC*(beta0*C1*(LQR))
				  +lamE*(-pow(beta1/beta0,2)*C1)
				  +lamE*(beta2/beta0*C1)
				  +lamF*(1./3.*(pow(C1,3)-3.*C1*C2+3.*C3))
				  +lamG*(-beta1/pow(beta0,2)*(beta1*C1 - beta0*pow(C1,2) + 2.*beta0*C2))
				  +lamH*(pow(beta1/beta0,2)*C1)
				  +lamG*(2.*beta1*C1*LQR)
				  +lamF*((-beta1*C1+beta0*pow(C1,2)-2.*beta0*C2)*LQR)
				  +lamF*(pow(beta0,2)*C1*pow(LQR,2))
				 );

	      C1 = 0.; C2 = ccoeff::C2qqb_2[idx]; C3 = ccoeff::C3qqb_2[idx];
	      aexpqqb_2[idx] = exp(
				 +lamB*(C1)
				 +lamC*((0.5*pow(C1,2) - C2))
				 +lamD*(-beta1/beta0*C1)
				 +lamC*(beta0*C1*(LQR))
				 +lamE*(-pow(beta1/beta0,2)*C1)
				 +lamE*(beta2/beta0*C1)
				 +lamF*(1./3.*(pow(C1,3)-3.*C1*C2+3.*C3))
				 +lamG*(-beta1/pow(beta0,2)*(beta1*C1 - beta0*pow(C1,2) + 2.*beta0*C2))
				 +lamH*(pow(beta1/beta0,2)*C1)
				 +lamG*(2.*beta1*C1*LQR)
				 +lamF*((-beta1*C1+beta0*pow(C1,2)-2.*beta0*C2)*LQR)
				 +lamF*(pow(beta0,2)*C1*pow(LQR,2))
				 );
	      C1 = 0.; C2 = ccoeff::C2qqp_2[idx]; C3 = ccoeff::C3qqp_2[idx];
	      aexpqqp_2[idx] = exp(
				 +lamB*(C1)
				 +lamC*((0.5*pow(C1,2) - C2))
				 +lamD*(-beta1/beta0*C1)
				 +lamC*(beta0*C1*(LQR))
				 +lamE*(-pow(beta1/beta0,2)*C1)
				 +lamE*(beta2/beta0*C1)
				 +lamF*(1./3.*(pow(C1,3)-3.*C1*C2+3.*C3))
				 +lamG*(-beta1/pow(beta0,2)*(beta1*C1 - beta0*pow(C1,2) + 2.*beta0*C2))
				 +lamH*(pow(beta1/beta0,2)*C1)
				 +lamG*(2.*beta1*C1*LQR)
				 +lamF*((-beta1*C1+beta0*pow(C1,2)-2.*beta0*C2)*LQR)
				 +lamF*(pow(beta0,2)*C1*pow(LQR,2))
				 );
	      C1 = 0.; C2 = ccoeff::C2qqbp_2[idx]; C3 = ccoeff::C3qqbp_2[idx];
	      aexpqqbp_2[idx] = exp(
				  +lamB*(C1)
				  +lamC*((0.5*pow(C1,2) - C2))
				  +lamD*(-beta1/beta0*C1)
				  +lamC*(beta0*C1*(LQR))
				  +lamE*(-pow(beta1/beta0,2)*C1)
				  +lamE*(beta2/beta0*C1)
				  +lamF*(1./3.*(pow(C1,3)-3.*C1*C2+3.*C3))
				  +lamG*(-beta1/pow(beta0,2)*(beta1*C1 - beta0*pow(C1,2) + 2.*beta0*C2))
				  +lamH*(pow(beta1/beta0,2)*C1)
				  +lamG*(2.*beta1*C1*LQR)
				  +lamF*((-beta1*C1+beta0*pow(C1,2)-2.*beta0*C2)*LQR)
				  +lamF*(pow(beta0,2)*C1*pow(LQR,2))
				 );
	    }
    }
}


//Standard formula for opts.expc == 3 QED
void expc::fullqed()
{
  if (opts.qedorder_expc == 1)
    return;
  
  double regfac = opts.expcreg; //Regularisation factor

  double LQR = LR-LQ;
  if (opts.mellin1d)
    {
      if (opts.qedorder_expc == 2)
	for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	  for (int i = 0; i < mellinint::mdim; i++)
	    {
	      int idx = anomalous::index(i,sign);
	      aexpqqqed[idx] = pow(aexpBqed,ccoeff::C1qqqed[idx]);
	      complex <double> C1 = (ccoeff::C2qg[idx]+regfac*ccoeff::C1qg_reg[idx])/(ccoeff::C1qg[idx]-aaqed*regfac*ccoeff::C1qg_reg[idx]);
	      aexpqg[idx] = pow(aexpBqed,C1);
	      aexpqqb[idx] = pow(aexpBqed,ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx]);
	      aexpqqp[idx] = pow(aexpBqed,ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx]);
	      aexpqqbp[idx] = pow(aexpBqed,ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx]);
	    }
  else
    {
      if (opts.qedorder_expc == 2)
	for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	  for (int i = 0; i < mellinint::mdim; i++)
	    {
	      int idx = anomalous::index(i,sign);
	      complex <double> C1;
	      aexpqqqed_1[idx] = pow(aexpBqed,ccoeff::C1qqqed_1[idx]);
	      C1 = (ccoeff::C2qg_1[idx]+regfac*ccoeff::C1qg_reg_1[idx])/(ccoeff::C1qg_1[idx]-aaqed*regfac*ccoeff::C1qg_reg_1[idx]);
	      aexpqgamma_1[idx] = pow(aexpBqed,C1);
	      aexpqqb_1[idx] = pow(aexpBqed,ccoeff::C3qqb_1[idx]/ccoeff::C2qqb_1[idx]);
	      aexpqqp_1[idx] = pow(aexpBqed,ccoeff::C3qqp_1[idx]/ccoeff::C2qqp_1[idx]);
	      aexpqqbp_1[idx] = pow(aexpBqed,ccoeff::C3qqbp_1[idx]/ccoeff::C2qqbp_1[idx]);
	      aexpqq_2[idx] = pow(aexpBqed,ccoeff::C1qq_2[idx]);
	      C1 = (ccoeff::C2qg_2[idx]+regfac*ccoeff::C1qg_reg_2[idx])/(ccoeff::C1qg_2[idx]-aass*regfac*ccoeff::C1qg_reg_2[idx]);
	      aexpqgamma_2[idx] = pow(aexpBqed,C1);
	      aexpqqb_2[idx] = pow(aexpBqed,ccoeff::C3qqb_2[idx]/ccoeff::C2qqb_2[idx]);
	      aexpqqp_2[idx] = pow(aexpBqed,ccoeff::C3qqp_2[idx]/ccoeff::C2qqp_2[idx]);
	      aexpqqbp_2[idx] = pow(aexpBqed,ccoeff::C3qqbp_2[idx]/ccoeff::C2qqbp_2[idx]);
	    }
    }
    }
}


























