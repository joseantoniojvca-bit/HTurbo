#include "muf.h"
#include "resint.h"
#include "pmom.h"
#include "anomalous.h"
#include "mesq.h"
#include "settings.h"
#include "scales.h"
#include "pdfevol.h"
#include "parton.h"
#include "ccoeff.h"
//#include "anomalous.h"
#include "pmom.h"
#include "resconst.h"
#include "expc.h"
#include "hcoeff.h"

#include <iostream>

double muf::LR  ;
double muf::LR2 ;
double muf::LR3 ;
double muf::LF  ;
double muf::LF2 ;
double muf::LF3 ;
double muf::LF4 ;
double muf::LQ  ;
double muf::LQ2 ;
double muf::LQ3 ;
double muf::LQ4 ;
double muf::LQ5 ;
double muf::LQ6 ;
double muf::LQ7 ;
double muf::LQF ;
double muf::LQF2;
double muf::LQF3;
double muf::LQF4;
double muf::A1  ;
double muf::A2  ;
double muf::A3  ;
double muf::A4  ;
double muf::B1  ;
double muf::B2  ;
double muf::B3  ;
double muf::B4  ;

complex <double> *muf::EC1qq   ;
complex <double> *muf::EC2qq   ;
complex <double> *muf::EC3qq   ;
complex <double> *muf::EC1qg   ;
complex <double> *muf::EC2qg   ;
complex <double> *muf::EC3qg   ;
complex <double> *muf::EC2qqp  ;
complex <double> *muf::EC2qqb  ;
complex <double> *muf::EC2qqbp ;
complex <double> *muf::EC3qqp  ;
complex <double> *muf::EC3qqb  ;
complex <double> *muf::EC3qqbp ;

//LL
complex <double> *muf::qqb;
complex <double> *muf::qqbqed;
//NLL
complex <double> *muf::qg;
complex <double> *muf::qg_1;
complex <double> *muf::qg_2;
complex <double> *muf::qgamma;
complex <double> *muf::qgamma_1;
complex <double> *muf::qgamma_2;
complex <double> *muf::gq_1;
complex <double> *muf::gq_2;
//NNLL
complex <double> *muf::qq;
complex <double> *muf::qq_1;
complex <double> *muf::qq_2;
complex <double> *muf::qqp;
complex <double> *muf::qqp_1;
complex <double> *muf::qqp_2;
complex <double> *muf::qqbp;
complex <double> *muf::qqbp_1;
complex <double> *muf::qqbp_2;
complex <double> *muf::gg;
//NNNLL
complex <double> *muf::qbg;
complex <double> *muf::qbg_1;
complex <double> *muf::qbg_2;
complex <double> *muf::qpg;
complex <double> *muf::qpg_1;
complex <double> *muf::qpg_2;
complex <double> *muf::qbpg;
complex <double> *muf::qbpg_1;
complex <double> *muf::qbpg_2;
//NNNNLL
complex <double> *muf::qbq;
complex <double> *muf::qpqbp;
complex <double> *muf::qbpqp;
complex <double> *muf::qbqp;
complex <double> *muf::qbqp_1;
complex <double> *muf::qbqp_2;
complex <double> *muf::qbqbp;
complex <double> *muf::qbqbp_1;
complex <double> *muf::qbqbp_2;
complex <double> *muf::qpqp;
complex <double> *muf::qpqp_1;
complex <double> *muf::qpqp_2;


//Split b-dependent C exponentiation
//mellin 1d
complex <double> *muf::qqb_eqqb;
complex <double> *muf::qg_eqqb;
complex <double> *muf::qqb_eqqbqed;
complex <double> *muf::qgamma_eqqb;
complex <double> *muf::qq_eqqb;
complex <double> *muf::qqp_eqqb;
complex <double> *muf::qqbp_eqqb;
complex <double> *muf::gg_eqqb;
complex <double> *muf::qbg_eqqb;
complex <double> *muf::qpg_eqqb;
complex <double> *muf::qbpg_eqqb;
complex <double> *muf::qbq_eqqb;
complex <double> *muf::qpqbp_eqqb;
complex <double> *muf::qbpqp_eqqb;
complex <double> *muf::qbqp_eqqb;
complex <double> *muf::qbqbp_eqqb;
complex <double> *muf::qpqp_eqqb;
complex <double> *muf::qqb_eqg;
complex <double> *muf::qg_eqg;
complex <double> *muf::qq_eqg;
complex <double> *muf::qqp_eqg;
complex <double> *muf::qqbp_eqg;
complex <double> *muf::gg_eqg;
complex <double> *muf::qbg_eqg;
complex <double> *muf::qpg_eqg;
complex <double> *muf::qbpg_eqg;
complex <double> *muf::qbq_eqg;
complex <double> *muf::qbqp_eqg;
complex <double> *muf::qbqbp_eqg;
complex <double> *muf::qpqp_eqg;
complex <double> *muf::qpqbp_eqg;
complex <double> *muf::qbpqp_eqg;
complex <double> *muf::qg_eqq;
complex <double> *muf::qq_eqq;
complex <double> *muf::qbg_eqq;
complex <double> *muf::gg_eqq;
complex <double> *muf::qbq_eqq;
complex <double> *muf::qbqp_eqq;
complex <double> *muf::qbqbp_eqq;
complex <double> *muf::qg_eqqp;
complex <double> *muf::qqp_eqqp;
complex <double> *muf::qbpg_eqqp;
complex <double> *muf::gg_eqqp;
complex <double> *muf::qbqp_eqqp;
complex <double> *muf::qpqp_eqqp;
complex <double> *muf::qbpqp_eqqp;
complex <double> *muf::qg_eqqbp;
complex <double> *muf::qqbp_eqqbp;
complex <double> *muf::qpg_eqqbp;
complex <double> *muf::gg_eqqbp;
complex <double> *muf::qbqbp_eqqbp;
complex <double> *muf::qpqp_eqqbp;
complex <double> *muf::qpqbp_eqqbp;
complex <double> *muf::qg_egg;
complex <double> *muf::gg_egg;
complex <double> *muf::qbg_egg;
complex <double> *muf::qpg_egg;
complex <double> *muf::qbpg_egg;
complex <double> *muf::qqb_egg;
complex <double> *muf::qq_egg;
complex <double> *muf::qqp_egg;
complex <double> *muf::qqbp_egg;
complex <double> *muf::qbq_egg;
complex <double> *muf::qbqp_egg;
complex <double> *muf::qbqbp_egg;
complex <double> *muf::qpqp_egg;
complex <double> *muf::qpqbp_egg;
complex <double> *muf::qbpqp_egg;
complex <double> *muf::qq_eqbg;
complex <double> *muf::gg_eqbg;
complex <double> *muf::qbg_eqbg;
complex <double> *muf::qbq_eqbg;
complex <double> *muf::qbqp_eqbg;
complex <double> *muf::qbqbp_eqbg;
complex <double> *muf::qqbp_eqpg;
complex <double> *muf::gg_eqpg;
complex <double> *muf::qpg_eqpg;
complex <double> *muf::qbqbp_eqpg;
complex <double> *muf::qpqp_eqpg;
complex <double> *muf::qpqbp_eqpg;
complex <double> *muf::qqp_eqbpg;
complex <double> *muf::gg_eqbpg;
complex <double> *muf::qbpg_eqbpg;
complex <double> *muf::qbqp_eqbpg;
complex <double> *muf::qpqp_eqbpg;
complex <double> *muf::qbpqp_eqbpg;

//mellin 2d
complex <double> *muf::qg_1_eqqb;
complex <double> *muf::qg_2_eqqb;
complex <double> *muf::qgamma_1_eqqb;
complex <double> *muf::qgamma_2_eqqb;
complex <double> *muf::qq_1_eqqb;
complex <double> *muf::qq_2_eqqb;
complex <double> *muf::qqp_1_eqqb;
complex <double> *muf::qqp_2_eqqb;
complex <double> *muf::qqbp_1_eqqb;
complex <double> *muf::qqbp_2_eqqb;
complex <double> *muf::qbg_1_eqqb;
complex <double> *muf::qbg_2_eqqb;
complex <double> *muf::qpg_1_eqqb;
complex <double> *muf::qpg_2_eqqb;
complex <double> *muf::qbpg_1_eqqb;
complex <double> *muf::qbpg_2_eqqb;
complex <double> *muf::qqb_eqg_1;
complex <double> *muf::qg_1_eqg_1;
complex <double> *muf::qg_2_eqg_1;
complex <double> *muf::qq_1_eqg_1;
complex <double> *muf::qqp_1_eqg_1;
complex <double> *muf::qqbp_1_eqg_1;
complex <double> *muf::gg_eqg_1;
complex <double> *muf::qbg_1_eqg_1;
complex <double> *muf::qbg_2_eqg_1;
complex <double> *muf::qpg_1_eqg_1;
complex <double> *muf::qpg_2_eqg_1;
complex <double> *muf::qbpg_1_eqg_1;
complex <double> *muf::qbpg_2_eqg_1;
complex <double> *muf::qqb_eqg_2;
complex <double> *muf::qg_1_eqg_2;
complex <double> *muf::qg_2_eqg_2;
complex <double> *muf::qq_2_eqg_2;
complex <double> *muf::qqp_2_eqg_2;
complex <double> *muf::qqbp_2_eqg_2;
complex <double> *muf::gg_eqg_2;
complex <double> *muf::qbg_1_eqg_2;
complex <double> *muf::qbg_2_eqg_2;
complex <double> *muf::qpg_1_eqg_2;
complex <double> *muf::qpg_2_eqg_2;
complex <double> *muf::qbpg_1_eqg_2;
complex <double> *muf::qbpg_2_eqg_2;
complex <double> *muf::qg_1_eqq_1;
complex <double> *muf::qq_1_eqq_1;
complex <double> *muf::qbg_2_eqq_1;
complex <double> *muf::qg_2_eqq_2;
complex <double> *muf::qq_2_eqq_2;
complex <double> *muf::qbg_1_eqq_2;
complex <double> *muf::qg_1_eqqp_1;
complex <double> *muf::qqp_1_eqqp_1;
complex <double> *muf::qbpg_2_eqqp_1;
complex <double> *muf::qg_2_eqqp_2;
complex <double> *muf::qqp_2_eqqp_2;
complex <double> *muf::qbpg_1_eqqp_2;
complex <double> *muf::qg_1_eqqbp_1;
complex <double> *muf::qqbp_1_eqqbp_1;
complex <double> *muf::qpg_2_eqqbp_1;
complex <double> *muf::qg_2_eqqbp_2;
complex <double> *muf::qqbp_2_eqqbp_2;
complex <double> *muf::qpg_1_eqqbp_2;
complex <double> *muf::qg_1_egg;
complex <double> *muf::qg_2_egg;
complex <double> *muf::qbg_1_egg;
complex <double> *muf::qbg_2_egg;
complex <double> *muf::qpg_1_egg;
complex <double> *muf::qpg_2_egg;
complex <double> *muf::qbpg_1_egg;
complex <double> *muf::qbpg_2_egg;

//using namespace ccoeff;
using namespace pmom;
using namespace resconst;
//using namespace resint;
using namespace ccoeff;

  //allocate memory
void muf::allocate()
{
  //if (opts.order == 0)
  //return;

  if (opts.mellin1d)
    {
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
      qpqbp = new complex <double> [mellinint::mdim*2];
      qbpqp = new complex <double> [mellinint::mdim*2];
      qbqp  = new complex <double> [mellinint::mdim*2];
      qbqbp = new complex <double> [mellinint::mdim*2];
      qpqp  = new complex <double> [mellinint::mdim*2];

      qqb_eqqb   = new complex <double> [mellinint::mdim*2];
      qg_eqqb    = new complex <double> [mellinint::mdim*2];
      qqb_eqqbqed   = new complex <double> [mellinint::mdim*2];
      qgamma_eqqb    = new complex <double> [mellinint::mdim*2];
      qq_eqqb    = new complex <double> [mellinint::mdim*2];
      qqp_eqqb   = new complex <double> [mellinint::mdim*2];
      qqbp_eqqb  = new complex <double> [mellinint::mdim*2];
      gg_eqqb    = new complex <double> [mellinint::mdim*2];
      qbg_eqqb   = new complex <double> [mellinint::mdim*2];
      qpg_eqqb   = new complex <double> [mellinint::mdim*2];
      qbpg_eqqb  = new complex <double> [mellinint::mdim*2];
      qbq_eqqb    = new complex <double> [mellinint::mdim*2];
      qpqbp_eqqb  = new complex <double> [mellinint::mdim*2];
      qbpqp_eqqb  = new complex <double> [mellinint::mdim*2];
      qbqp_eqqb   = new complex <double> [mellinint::mdim*2];
      qbqbp_eqqb  = new complex <double> [mellinint::mdim*2];
      qpqp_eqqb   = new complex <double> [mellinint::mdim*2];
      qqb_eqg    = new complex <double> [mellinint::mdim*2];
      qg_eqg     = new complex <double> [mellinint::mdim*2];
      qq_eqg     = new complex <double> [mellinint::mdim*2];
      qqp_eqg    = new complex <double> [mellinint::mdim*2];
      qqbp_eqg   = new complex <double> [mellinint::mdim*2];
      gg_eqg     = new complex <double> [mellinint::mdim*2];
      qbg_eqg    = new complex <double> [mellinint::mdim*2];
      qpg_eqg    = new complex <double> [mellinint::mdim*2];
      qbpg_eqg   = new complex <double> [mellinint::mdim*2];
      qbq_eqg   = new complex <double> [mellinint::mdim*2];
      qbqp_eqg  = new complex <double> [mellinint::mdim*2];
      qbqbp_eqg = new complex <double> [mellinint::mdim*2];
      qpqp_eqg  = new complex <double> [mellinint::mdim*2];
      qpqbp_eqg = new complex <double> [mellinint::mdim*2];
      qbpqp_eqg = new complex <double> [mellinint::mdim*2];
      qg_eqq     = new complex <double> [mellinint::mdim*2];
      qq_eqq     = new complex <double> [mellinint::mdim*2];
      qbg_eqq    = new complex <double> [mellinint::mdim*2];
      gg_eqq     = new complex <double> [mellinint::mdim*2];
      qbq_eqq    = new complex <double> [mellinint::mdim*2];
      qbqp_eqq   = new complex <double> [mellinint::mdim*2];
      qbqbp_eqq  = new complex <double> [mellinint::mdim*2];
      qg_eqqp    = new complex <double> [mellinint::mdim*2];
      qqp_eqqp  = new complex <double> [mellinint::mdim*2];
      qbpg_eqqp   = new complex <double> [mellinint::mdim*2];
      gg_eqqp     = new complex <double> [mellinint::mdim*2];
      qbqp_eqqp   = new complex <double> [mellinint::mdim*2];
      qpqp_eqqp   = new complex <double> [mellinint::mdim*2];
      qbpqp_eqqp  = new complex <double> [mellinint::mdim*2];
      qg_eqqbp   = new complex <double> [mellinint::mdim*2];
      qqbp_eqqbp  = new complex <double> [mellinint::mdim*2];
      qpg_eqqbp = new complex <double> [mellinint::mdim*2];
      gg_eqqbp     = new complex <double> [mellinint::mdim*2];
      qbqbp_eqqbp  = new complex <double> [mellinint::mdim*2];
      qpqp_eqqbp   = new complex <double> [mellinint::mdim*2];
      qpqbp_eqqbp  = new complex <double> [mellinint::mdim*2];
      qg_egg     = new complex <double> [mellinint::mdim*2];
      gg_egg     = new complex <double> [mellinint::mdim*2];
      qbg_egg    = new complex <double> [mellinint::mdim*2];
      qpg_egg    = new complex <double> [mellinint::mdim*2];
      qbpg_egg   = new complex <double> [mellinint::mdim*2];
      qqb_egg       = new complex <double> [mellinint::mdim*2];
      qq_egg	    = new complex <double> [mellinint::mdim*2];
      qqp_egg       = new complex <double> [mellinint::mdim*2];
      qqbp_egg      = new complex <double> [mellinint::mdim*2];
      qbq_egg       = new complex <double> [mellinint::mdim*2];
      qbqp_egg      = new complex <double> [mellinint::mdim*2];
      qbqbp_egg     = new complex <double> [mellinint::mdim*2];
      qpqp_egg      = new complex <double> [mellinint::mdim*2];
      qpqbp_egg     = new complex <double> [mellinint::mdim*2];
      qbpqp_egg     = new complex <double> [mellinint::mdim*2];
      qq_eqbg       = new complex <double> [mellinint::mdim*2];
      gg_eqbg       = new complex <double> [mellinint::mdim*2];
      qbg_eqbg      = new complex <double> [mellinint::mdim*2];
      qbq_eqbg      = new complex <double> [mellinint::mdim*2];
      qbqp_eqbg     = new complex <double> [mellinint::mdim*2];
      qbqbp_eqbg    = new complex <double> [mellinint::mdim*2];
      qqbp_eqpg     = new complex <double> [mellinint::mdim*2];
      gg_eqpg       = new complex <double> [mellinint::mdim*2];
      qpg_eqpg      = new complex <double> [mellinint::mdim*2];
      qbqbp_eqpg    = new complex <double> [mellinint::mdim*2];
      qpqp_eqpg     = new complex <double> [mellinint::mdim*2];
      qpqbp_eqpg    = new complex <double> [mellinint::mdim*2];
      qqp_eqbpg     = new complex <double> [mellinint::mdim*2];
      gg_eqbpg      = new complex <double> [mellinint::mdim*2];
      qbpg_eqbpg    = new complex <double> [mellinint::mdim*2];
      qbqp_eqbpg    = new complex <double> [mellinint::mdim*2];
      qpqp_eqbpg    = new complex <double> [mellinint::mdim*2];
      qbpqp_eqbpg   = new complex <double> [mellinint::mdim*2];

      EC1qq    = new complex <double> [mellinint::mdim*2];
      EC2qq    = new complex <double> [mellinint::mdim*2];
      EC3qq    = new complex <double> [mellinint::mdim*2];
      EC1qg    = new complex <double> [mellinint::mdim*2];
      EC2qg    = new complex <double> [mellinint::mdim*2];
      EC3qg    = new complex <double> [mellinint::mdim*2];
      EC2qqp   = new complex <double> [mellinint::mdim*2];
      EC2qqb   = new complex <double> [mellinint::mdim*2];
      EC2qqbp  = new complex <double> [mellinint::mdim*2];
      EC3qqp   = new complex <double> [mellinint::mdim*2];
      EC3qqb   = new complex <double> [mellinint::mdim*2];
      EC3qqbp  = new complex <double> [mellinint::mdim*2];
    }
  else
    {
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
      qpqbp   = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qbpqp   = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qbqp_1  = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qbqp_2  = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qbqbp_1 = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qbqbp_2 = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qpqp_1  = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qpqp_2  = new complex <double> [mellinint::mdim*mellinint::mdim*2];  

      qqb_eqqb       = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qqb_eqqbqed       = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      gg_eqqb	     = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      gg_egg	     = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qg_1_eqqb      = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qg_2_eqqb	     = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qgamma_1_eqqb      = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qgamma_2_eqqb	 = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qq_1_eqqb	     = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qq_2_eqqb	     = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qqp_1_eqqb     = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qqp_2_eqqb     = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qqbp_1_eqqb    = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qqbp_2_eqqb    = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qbg_1_eqqb     = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qbg_2_eqqb     = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qpg_1_eqqb     = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qpg_2_eqqb     = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qbpg_1_eqqb    = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qbpg_2_eqqb    = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qqb_eqg_1	     = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qg_1_eqg_1     = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qg_2_eqg_1     = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qq_1_eqg_1     = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qqp_1_eqg_1    = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qqbp_1_eqg_1   = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      gg_eqg_1	     = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qbg_1_eqg_1    = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qbg_2_eqg_1    = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qpg_1_eqg_1    = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qpg_2_eqg_1    = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qbpg_1_eqg_1   = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qbpg_2_eqg_1   = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qqb_eqg_2	     = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qg_1_eqg_2     = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qg_2_eqg_2     = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qq_2_eqg_2     = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qqp_2_eqg_2    = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qqbp_2_eqg_2   = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      gg_eqg_2	     = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qbg_1_eqg_2    = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qbg_2_eqg_2    = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qpg_1_eqg_2    = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qpg_2_eqg_2    = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qbpg_1_eqg_2   = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qbpg_2_eqg_2   = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qg_1_eqq_1     = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qq_1_eqq_1     = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qbg_2_eqq_1    = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qg_2_eqq_2     = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qq_2_eqq_2     = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qbg_1_eqq_2    = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qg_1_eqqp_1    = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qqp_1_eqqp_1   = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qbpg_2_eqqp_1  = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qg_2_eqqp_2    = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qqp_2_eqqp_2   = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qbpg_1_eqqp_2  = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qg_1_eqqbp_1   = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qqbp_1_eqqbp_1 = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qpg_2_eqqbp_1  = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qg_2_eqqbp_2   = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qqbp_2_eqqbp_2 = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qpg_1_eqqbp_2  = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qg_1_egg	     = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qg_2_egg	     = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qbg_1_egg	     = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qbg_2_egg	     = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qpg_1_egg	     = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qpg_2_egg	     = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qbpg_1_egg     = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qbpg_2_egg     = new complex <double> [mellinint::mdim*mellinint::mdim*2]; 
    }
}

void muf::free()
{
  //if (opts.order == 0)
  //return;

  if (opts.mellin1d)
    {
      delete[] qqb;
      delete[] qg;
      delete[] qqbqed;
      delete[] qgamma;
      delete[] qq;
      delete[] qqp;  
      delete[] qqbp;
      delete[] gg;
      delete[] qbg;
      delete[] qpg;
      delete[] qbpg;
      delete[] qbq;
      delete[] qpqbp;
      delete[] qbpqp;
      delete[] qbqp;
      delete[] qbqbp;
      delete[] qpqp;

      delete[] qqb_eqqb;
      delete[] qg_eqqb;
      delete[] qqb_eqqbqed;
      delete[] qgamma_eqqb;
      delete[] qq_eqqb;
      delete[] qqp_eqqb;
      delete[] qqbp_eqqb;
      delete[] gg_eqqb;
      delete[] qbg_eqqb;
      delete[] qpg_eqqb;
      delete[] qbpg_eqqb;
      delete[] qbq_eqqb;
      delete[] qpqbp_eqqb;
      delete[] qbpqp_eqqb;
      delete[] qbqp_eqqb;
      delete[] qbqbp_eqqb;
      delete[] qpqp_eqqb;
      delete[] qqb_eqg;
      delete[] qg_eqg;
      delete[] qq_eqg;
      delete[] qqp_eqg;
      delete[] qqbp_eqg;
      delete[] gg_eqg;
      delete[] qbg_eqg;
      delete[] qpg_eqg;
      delete[] qbpg_eqg;
      delete[] qbq_eqg;
      delete[] qbqp_eqg;
      delete[] qbqbp_eqg;
      delete[] qpqp_eqg;
      delete[] qpqbp_eqg;
      delete[] qbpqp_eqg;
      delete[] qg_eqq;
      delete[] qq_eqq;
      delete[] qbg_eqq;
      delete[] gg_eqq;
      delete[] qbq_eqq;
      delete[] qbqp_eqq;
      delete[] qbqbp_eqq;
      delete[] qg_eqqp;
      delete[] qqp_eqqp;
      delete[] qbpg_eqqp;
      delete[] gg_eqqp;
      delete[] qbqp_eqqp;
      delete[] qpqp_eqqp;
      delete[] qbpqp_eqqp;
      delete[] qg_eqqbp;
      delete[] qqbp_eqqbp;
      delete[] qpg_eqqbp;
      delete[] gg_eqqbp;
      delete[] qbqbp_eqqbp;
      delete[] qpqp_eqqbp;
      delete[] qpqbp_eqqbp;
      delete[] qg_egg;
      delete[] gg_egg;
      delete[] qbg_egg;
      delete[] qpg_egg;
      delete[] qbpg_egg;
      delete[] qqb_egg;
      delete[] qq_egg;
      delete[] qqp_egg;
      delete[] qqbp_egg;
      delete[] qbq_egg;
      delete[] qbqp_egg;
      delete[] qbqbp_egg;
      delete[] qpqp_egg;
      delete[] qpqbp_egg;
      delete[] qbpqp_egg;
      delete[] qq_eqbg;
      delete[] gg_eqbg;
      delete[] qbg_eqbg;
      delete[] qbq_eqbg;
      delete[] qbqp_eqbg;
      delete[] qbqbp_eqbg;
      delete[] qqbp_eqpg;
      delete[] gg_eqpg;
      delete[] qpg_eqpg;
      delete[] qbqbp_eqpg;
      delete[] qpqp_eqpg;
      delete[] qpqbp_eqpg;
      delete[] qqp_eqbpg;
      delete[] gg_eqbpg;
      delete[] qbpg_eqbpg;
      delete[] qbqp_eqbpg;
      delete[] qpqp_eqbpg;
      delete[] qbpqp_eqbpg;

      delete[] EC1qq  ;
      delete[] EC2qq  ;
      delete[] EC3qq  ;
      delete[] EC1qg  ;
      delete[] EC2qg  ;
      delete[] EC3qg  ;
      delete[] EC2qqp ;
      delete[] EC2qqb ;
      delete[] EC2qqbp;
      delete[] EC3qqp ;
      delete[] EC3qqb ;
      delete[] EC3qqbp;
    }
  else
    {
      delete[] qqb;
      delete[] qg_1;
      delete[] qg_2;
      delete[] gq_1;
      delete[] gq_2;
      delete[] qqbqed;
      delete[] qgamma_1;
      delete[] qgamma_2;
      delete[] qq_1;
      delete[] qq_2;
      delete[] qqp_1;
      delete[] qqp_2;
      delete[] qqbp_1;
      delete[] qqbp_2;
      delete[] gg;
      delete[] qbg_1;
      delete[] qbg_2;
      delete[] qpg_1;
      delete[] qpg_2;
      delete[] qbpg_1;
      delete[] qbpg_2;
      delete[] qbq;
      delete[] qpqbp;
      delete[] qbpqp;
      delete[] qbqp_1;
      delete[] qbqp_2;
      delete[] qbqbp_1;
      delete[] qbqbp_2;
      delete[] qpqp_1;
      delete[] qpqp_2;

      delete[] qqb_eqqb;
      delete[] qqb_eqqbqed;
      delete[] gg_eqqb;
      delete[] gg_egg;
      delete[] qg_1_eqqb;
      delete[] qg_2_eqqb;
      delete[] qgamma_1_eqqb;
      delete[] qgamma_2_eqqb;
      delete[] qq_1_eqqb;
      delete[] qq_2_eqqb;
      delete[] qqp_1_eqqb;
      delete[] qqp_2_eqqb;
      delete[] qqbp_1_eqqb;
      delete[] qqbp_2_eqqb;
      delete[] qbg_1_eqqb;
      delete[] qbg_2_eqqb;
      delete[] qpg_1_eqqb;
      delete[] qpg_2_eqqb;
      delete[] qbpg_1_eqqb;
      delete[] qbpg_2_eqqb;
      delete[] qqb_eqg_1;
      delete[] qg_1_eqg_1;
      delete[] qg_2_eqg_1;
      delete[] qq_1_eqg_1;
      delete[] qqp_1_eqg_1;
      delete[] qqbp_1_eqg_1;
      delete[] gg_eqg_1;
      delete[] qbg_1_eqg_1;
      delete[] qbg_2_eqg_1;
      delete[] qpg_1_eqg_1;
      delete[] qpg_2_eqg_1;
      delete[] qbpg_1_eqg_1;
      delete[] qbpg_2_eqg_1;
      delete[] qqb_eqg_2;
      delete[] qg_1_eqg_2;
      delete[] qg_2_eqg_2;
      delete[] qq_2_eqg_2;
      delete[] qqp_2_eqg_2;
      delete[] qqbp_2_eqg_2;
      delete[] gg_eqg_2;
      delete[] qbg_1_eqg_2;
      delete[] qbg_2_eqg_2;
      delete[] qpg_1_eqg_2;
      delete[] qpg_2_eqg_2;
      delete[] qbpg_1_eqg_2;
      delete[] qbpg_2_eqg_2;
      delete[] qg_1_eqq_1;
      delete[] qq_1_eqq_1;
      delete[] qbg_2_eqq_1;
      delete[] qg_2_eqq_2;
      delete[] qq_2_eqq_2;
      delete[] qbg_1_eqq_2;
      delete[] qg_1_eqqp_1;
      delete[] qqp_1_eqqp_1;
      delete[] qbpg_2_eqqp_1;
      delete[] qg_2_eqqp_2;
      delete[] qqp_2_eqqp_2;
      delete[] qbpg_1_eqqp_2;
      delete[] qg_1_eqqbp_1;
      delete[] qqbp_1_eqqbp_1;
      delete[] qpg_2_eqqbp_1;
      delete[] qg_2_eqqbp_2;
      delete[] qqbp_2_eqqbp_2;
      delete[] qpg_1_eqqbp_2;
      delete[] qg_1_egg;
      delete[] qg_2_egg;
      delete[] qbg_1_egg;
      delete[] qbg_2_egg;
      delete[] qpg_1_egg;
      delete[] qpg_2_egg;
      delete[] qbpg_1_egg;
      delete[] qbpg_2_egg;
    }
}

void muf::reset()
{
  //if (opts.order == 0)
  //return;

  if (opts.mellin1d)
    { 
      fill(qqb , qqb  + mellinint::mdim*2, 0.);
      fill(qg  , qg   + mellinint::mdim*2, 0.);
      fill(qqbqed , qqbqed  + mellinint::mdim*2, 0.);
      fill(qgamma  , qgamma   + mellinint::mdim*2, 0.);
      fill(qq  , qq   + mellinint::mdim*2, 0.);
      fill(qqp , qqp  + mellinint::mdim*2, 0.);
      fill(qqbp, qqbp + mellinint::mdim*2, 0.);
      fill(gg  , gg   + mellinint::mdim*2, 0.);
      fill(qbg , qbg  + mellinint::mdim*2, 0.);
      fill(qpg , qpg  + mellinint::mdim*2, 0.);
      fill(qbpg, qbpg + mellinint::mdim*2, 0.);
      fill(qbq   , qbq   + mellinint::mdim*2, 0.);
      fill(qpqbp , qpqbp + mellinint::mdim*2, 0.);
      fill(qbpqp , qbpqp + mellinint::mdim*2, 0.);
      fill(qbqp  , qbqp  + mellinint::mdim*2, 0.);
      fill(qbqbp , qbqbp + mellinint::mdim*2, 0.);
      fill(qpqp  , qpqp  + mellinint::mdim*2, 0.);

      fill(qqb_eqqb   , qqb_eqqb   + mellinint::mdim*2, 0.);
      fill(qg_eqqb    ,	qg_eqqb    + mellinint::mdim*2, 0.);
      fill(qqb_eqqbqed   , qqb_eqqbqed   + mellinint::mdim*2, 0.);
      fill(qgamma_eqqb    ,	qgamma_eqqb    + mellinint::mdim*2, 0.);
      fill(qq_eqqb    ,	qq_eqqb    + mellinint::mdim*2, 0.);
      fill(qqp_eqqb   ,	qqp_eqqb   + mellinint::mdim*2, 0.);
      fill(qqbp_eqqb  ,	qqbp_eqqb  + mellinint::mdim*2, 0.);
      fill(gg_eqqb    ,	gg_eqqb    + mellinint::mdim*2, 0.);
      fill(qbg_eqqb   ,	qbg_eqqb   + mellinint::mdim*2, 0.);
      fill(qpg_eqqb   ,	qpg_eqqb   + mellinint::mdim*2, 0.);
      fill(qbpg_eqqb  ,	qbpg_eqqb  + mellinint::mdim*2, 0.);
      fill(qbq_eqqb   , qbq_eqqb   + mellinint::mdim*2, 0.);
      fill(qpqbp_eqqb ,	qpqbp_eqqb + mellinint::mdim*2, 0.);
      fill(qbpqp_eqqb ,	qbpqp_eqqb + mellinint::mdim*2, 0.);
      fill(qbqp_eqqb  ,	qbqp_eqqb  + mellinint::mdim*2, 0.);
      fill(qbqbp_eqqb ,	qbqbp_eqqb + mellinint::mdim*2, 0.);
      fill(qpqp_eqqb  ,	qpqp_eqqb  + mellinint::mdim*2, 0.);
      fill(qqb_eqg    ,	qqb_eqg    + mellinint::mdim*2, 0.);
      fill(qg_eqg     ,	qg_eqg     + mellinint::mdim*2, 0.);
      fill(qq_eqg     ,	qq_eqg     + mellinint::mdim*2, 0.);
      fill(qqp_eqg    ,	qqp_eqg    + mellinint::mdim*2, 0.);
      fill(qqbp_eqg   ,	qqbp_eqg   + mellinint::mdim*2, 0.);
      fill(gg_eqg     ,	gg_eqg     + mellinint::mdim*2, 0.);
      fill(qbg_eqg    ,	qbg_eqg    + mellinint::mdim*2, 0.);
      fill(qpg_eqg    ,	qpg_eqg    + mellinint::mdim*2, 0.);
      fill(qbpg_eqg   ,	qbpg_eqg   + mellinint::mdim*2, 0.);
      fill(qbq_eqg    , qbq_eqg    + mellinint::mdim*2, 0.);
      fill(qbqp_eqg   ,	qbqp_eqg   + mellinint::mdim*2, 0.);
      fill(qbqbp_eqg  ,	qbqbp_eqg  + mellinint::mdim*2, 0.);
      fill(qpqp_eqg   ,	qpqp_eqg   + mellinint::mdim*2, 0.);
      fill(qpqbp_eqg  ,	qpqbp_eqg  + mellinint::mdim*2, 0.);
      fill(qbpqp_eqg  ,	qbpqp_eqg  + mellinint::mdim*2, 0.);
      fill(qg_eqq     ,	qg_eqq     + mellinint::mdim*2, 0.);
      fill(qq_eqq     ,	qq_eqq     + mellinint::mdim*2, 0.);
      fill(qbg_eqq    ,	qbg_eqq    + mellinint::mdim*2, 0.);
      fill(gg_eqq    , gg_eqq     + mellinint::mdim*2, 0.);
      fill(qbq_eqq   , qbq_eqq    + mellinint::mdim*2, 0.);
      fill(qbqp_eqq  , qbqp_eqq   + mellinint::mdim*2, 0.);
      fill(qbqbp_eqq , qbqbp_eqq  + mellinint::mdim*2, 0.);
      fill(qg_eqqp    ,	qg_eqqp    + mellinint::mdim*2, 0.);
      fill(qqp_eqqp   ,	qqp_eqqp  + mellinint::mdim*2, 0.);
      fill(qbpg_eqqp  ,	qbpg_eqqp   + mellinint::mdim*2, 0.);
      fill(gg_eqqp    ,gg_eqqp     + mellinint::mdim*2, 0.);
      fill(qbqp_eqqp  ,qbqp_eqqp   + mellinint::mdim*2, 0.);
      fill(qpqp_eqqp  ,qpqp_eqqp   + mellinint::mdim*2, 0.);
      fill(qbpqp_eqqp ,qbpqp_eqqp  + mellinint::mdim*2, 0.);
      fill(qg_eqqbp   ,	qg_eqqbp   + mellinint::mdim*2, 0.);
      fill(qqbp_eqqbp ,	qqbp_eqqbp  + mellinint::mdim*2, 0.);
      fill(qpg_eqqbp  ,	qpg_eqqbp + mellinint::mdim*2, 0.);
      fill(gg_eqqbp    , gg_eqqbp     + mellinint::mdim*2, 0.);
      fill(qbqbp_eqqbp , qbqbp_eqqbp  + mellinint::mdim*2, 0.);
      fill(qpqp_eqqbp  , qpqp_eqqbp   + mellinint::mdim*2, 0.);
      fill(qpqbp_eqqbp , qpqbp_eqqbp  + mellinint::mdim*2, 0.);
      fill(qg_egg     ,	qg_egg     + mellinint::mdim*2, 0.);
      fill(gg_egg     ,	gg_egg     + mellinint::mdim*2, 0.);
      fill(qbg_egg    ,	qbg_egg    + mellinint::mdim*2, 0.);
      fill(qpg_egg    ,	qpg_egg    + mellinint::mdim*2, 0.);
      fill(qbpg_egg   ,	qbpg_egg   + mellinint::mdim*2, 0.);
      fill(qqb_egg     , qqb_egg     + mellinint::mdim*2, 0.);
      fill(qq_egg      , qq_egg      + mellinint::mdim*2, 0.);
      fill(qqp_egg     , qqp_egg     + mellinint::mdim*2, 0.);
      fill(qqbp_egg    , qqbp_egg    + mellinint::mdim*2, 0.);
      fill(qbq_egg     , qbq_egg     + mellinint::mdim*2, 0.);
      fill(qbqp_egg    , qbqp_egg    + mellinint::mdim*2, 0.);
      fill(qbqbp_egg   , qbqbp_egg   + mellinint::mdim*2, 0.);
      fill(qpqp_egg    , qpqp_egg    + mellinint::mdim*2, 0.);
      fill(qpqbp_egg   , qpqbp_egg   + mellinint::mdim*2, 0.);
      fill(qbpqp_egg   , qbpqp_egg   + mellinint::mdim*2, 0.);
      fill(qq_eqbg     , qq_eqbg     + mellinint::mdim*2, 0.);
      fill(gg_eqbg     , gg_eqbg     + mellinint::mdim*2, 0.);
      fill(qbg_eqbg    , qbg_eqbg    + mellinint::mdim*2, 0.);
      fill(qbq_eqbg    , qbq_eqbg    + mellinint::mdim*2, 0.);
      fill(qbqp_eqbg   , qbqp_eqbg   + mellinint::mdim*2, 0.);
      fill(qbqbp_eqbg  , qbqbp_eqbg  + mellinint::mdim*2, 0.);
      fill(qqbp_eqpg   , qqbp_eqpg   + mellinint::mdim*2, 0.);
      fill(gg_eqpg     , gg_eqpg     + mellinint::mdim*2, 0.);
      fill(qpg_eqpg    , qpg_eqpg    + mellinint::mdim*2, 0.);
      fill(qbqbp_eqpg  , qbqbp_eqpg  + mellinint::mdim*2, 0.);
      fill(qpqp_eqpg   , qpqp_eqpg   + mellinint::mdim*2, 0.);
      fill(qpqbp_eqpg  , qpqbp_eqpg  + mellinint::mdim*2, 0.);
      fill(qqp_eqbpg   , qqp_eqbpg   + mellinint::mdim*2, 0.);
      fill(gg_eqbpg    , gg_eqbpg    + mellinint::mdim*2, 0.);
      fill(qbpg_eqbpg  , qbpg_eqbpg  + mellinint::mdim*2, 0.);
      fill(qbqp_eqbpg  , qbqp_eqbpg  + mellinint::mdim*2, 0.);
      fill(qpqp_eqbpg  , qpqp_eqbpg  + mellinint::mdim*2, 0.);
      fill(qbpqp_eqbpg , qbpqp_eqbpg + mellinint::mdim*2, 0.);
    }
  else
    {
      fill(qqb   ,qqb    + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qg_1  ,qg_1   + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qg_2  ,qg_2   + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(gq_1  ,gq_1   + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(gq_2  ,gq_2   + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qqbqed   ,qqbqed    + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qgamma_1  ,qgamma_1   + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qgamma_2  ,qgamma_2   + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qq_1  ,qq_1   + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qq_2  ,qq_2   + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qqp_1 ,qqp_1  + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qqp_2 ,qqp_2  + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qqbp_1,qqbp_1 + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qqbp_2,qqbp_2 + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(gg    ,gg     + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qbg_1 ,qbg_1  + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qbg_2 ,qbg_2  + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qpg_1 ,qpg_1  + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qpg_2 ,qpg_2  + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qbpg_1,qbpg_1 + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qbpg_2,qbpg_2 + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qbq     , qbq     + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qpqbp   , qpqbp   + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qbpqp   , qbpqp   + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qbqp_1  , qbqp_1  + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qbqp_2  , qbqp_2  + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qbqbp_1 , qbqbp_1 + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qbqbp_2 , qbqbp_2 + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qpqp_1  , qpqp_1  + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qpqp_2  , qpqp_2  + mellinint::mdim*mellinint::mdim*2, 0.);  

      fill(qqb_eqqb       ,qqb_eqqb       + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qqb_eqqbqed    ,qqb_eqqbqed       + mellinint::mdim*mellinint::mdim*2, 0.);

      fill(gg_eqqb	  ,gg_eqqb	  + mellinint::mdim*mellinint::mdim*2, 0.);  
      fill(gg_egg	  ,gg_egg	  + mellinint::mdim*mellinint::mdim*2, 0.);  
      fill(qg_1_eqqb      ,qg_1_eqqb      + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qg_2_eqqb	  ,qg_2_eqqb	  + mellinint::mdim*mellinint::mdim*2, 0.);  

      fill(qgamma_1_eqqb  ,qgamma_1_eqqb      + mellinint::mdim*mellinint::mdim*2, 0.);

      fill(qgamma_2_eqqb  ,qgamma_2_eqqb	  + mellinint::mdim*mellinint::mdim*2, 0.);  

      fill(qq_1_eqqb	  ,qq_1_eqqb	  + mellinint::mdim*mellinint::mdim*2, 0.);  
      fill(qq_2_eqqb	  ,qq_2_eqqb	  + mellinint::mdim*mellinint::mdim*2, 0.);  
      fill(qqp_1_eqqb     ,qqp_1_eqqb     + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qqp_2_eqqb     ,qqp_2_eqqb     + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qqbp_1_eqqb    ,qqbp_1_eqqb    + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qqbp_2_eqqb    ,qqbp_2_eqqb    + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qbg_1_eqqb     ,qbg_1_eqqb     + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qbg_2_eqqb     ,qbg_2_eqqb     + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qpg_1_eqqb     ,qpg_1_eqqb     + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qpg_2_eqqb     ,qpg_2_eqqb     + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qbpg_1_eqqb    ,qbpg_1_eqqb    + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qbpg_2_eqqb    ,qbpg_2_eqqb    + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qqb_eqg_1	  ,qqb_eqg_1	  + mellinint::mdim*mellinint::mdim*2, 0.);  
      fill(qg_1_eqg_1     ,qg_1_eqg_1     + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qg_2_eqg_1     ,qg_2_eqg_1     + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qq_1_eqg_1     ,qq_1_eqg_1     + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qqp_1_eqg_1    ,qqp_1_eqg_1    + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qqbp_1_eqg_1   ,qqbp_1_eqg_1   + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(gg_eqg_1	  ,gg_eqg_1	  + mellinint::mdim*mellinint::mdim*2, 0.);  
      fill(qbg_1_eqg_1    ,qbg_1_eqg_1    + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qbg_2_eqg_1    ,qbg_2_eqg_1    + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qpg_1_eqg_1    ,qpg_1_eqg_1    + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qpg_2_eqg_1    ,qpg_2_eqg_1    + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qbpg_1_eqg_1   ,qbpg_1_eqg_1   + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qbpg_2_eqg_1   ,qbpg_2_eqg_1   + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qqb_eqg_2	  ,qqb_eqg_2	  + mellinint::mdim*mellinint::mdim*2, 0.);  
      fill(qg_1_eqg_2     ,qg_1_eqg_2     + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qg_2_eqg_2     ,qg_2_eqg_2     + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qq_2_eqg_2     ,qq_2_eqg_2     + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qqp_2_eqg_2    ,qqp_2_eqg_2    + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qqbp_2_eqg_2   ,qqbp_2_eqg_2   + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(gg_eqg_2	  ,gg_eqg_2	  + mellinint::mdim*mellinint::mdim*2, 0.);  
      fill(qbg_1_eqg_2    ,qbg_1_eqg_2    + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qbg_2_eqg_2    ,qbg_2_eqg_2    + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qpg_1_eqg_2    ,qpg_1_eqg_2    + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qpg_2_eqg_2    ,qpg_2_eqg_2    + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qbpg_1_eqg_2   ,qbpg_1_eqg_2   + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qbpg_2_eqg_2   ,qbpg_2_eqg_2   + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qg_1_eqq_1     ,qg_1_eqq_1     + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qq_1_eqq_1     ,qq_1_eqq_1     + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qbg_2_eqq_1    ,qbg_2_eqq_1    + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qg_2_eqq_2     ,qg_2_eqq_2     + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qq_2_eqq_2     ,qq_2_eqq_2     + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qbg_1_eqq_2    ,qbg_1_eqq_2    + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qg_1_eqqp_1    ,qg_1_eqqp_1    + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qqp_1_eqqp_1   ,qqp_1_eqqp_1   + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qbpg_2_eqqp_1  ,qbpg_2_eqqp_1  + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qg_2_eqqp_2    ,qg_2_eqqp_2    + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qqp_2_eqqp_2   ,qqp_2_eqqp_2   + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qbpg_1_eqqp_2  ,qbpg_1_eqqp_2  + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qg_1_eqqbp_1   ,qg_1_eqqbp_1   + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qqbp_1_eqqbp_1 ,qqbp_1_eqqbp_1 + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qpg_2_eqqbp_1  ,qpg_2_eqqbp_1  + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qg_2_eqqbp_2   ,qg_2_eqqbp_2   + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qqbp_2_eqqbp_2 ,qqbp_2_eqqbp_2 + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qpg_1_eqqbp_2  ,qpg_1_eqqbp_2  + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qg_1_egg	  ,qg_1_egg	  + mellinint::mdim*mellinint::mdim*2, 0.);  
      fill(qg_2_egg	  ,qg_2_egg	  + mellinint::mdim*mellinint::mdim*2, 0.);  
      fill(qbg_1_egg	  ,qbg_1_egg	  + mellinint::mdim*mellinint::mdim*2, 0.);  
      fill(qbg_2_egg	  ,qbg_2_egg	  + mellinint::mdim*mellinint::mdim*2, 0.);  
      fill(qpg_1_egg	  ,qpg_1_egg	  + mellinint::mdim*mellinint::mdim*2, 0.);  
      fill(qpg_2_egg	  ,qpg_2_egg	  + mellinint::mdim*mellinint::mdim*2, 0.);  
      fill(qbpg_1_egg     ,qbpg_1_egg     + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qbpg_2_egg     ,qbpg_2_egg     + mellinint::mdim*mellinint::mdim*2, 0.);
      
    }
}

void muf::ccalc()
{
  //LO
  if (opts.order_hcoef == 0)
    return;

  if (opts.mellin1d) //rapidity integrated
    ccalc1d();
  else
    ccalc2d();
}

void muf::ccalcb(complex <double> b)
{//cout << " ccalcb " << endl;
  //LO
  if (opts.order_hcoef == 0)
    return;

  if (opts.mellin1d) //rapidity integrated
    for (int sign = mesq::positive; sign <= mesq::negative; sign++)
      for (int i = 0; i < mellinint::mdim; i++)
	{ //cout << " aqui " << endl;
	  int idx = pmom::index(i,sign);
	
	  qqb[idx]  = qqb_eqqb[idx]*expc::qqb[idx];
	  qg[idx]   = qg_eqqb[idx] *expc::qqb[idx];
	  qqbqed[idx]  = qqb_eqqbqed[idx]*expc::qqbqed[idx];
	  qgamma[idx]   = qgamma_eqqb[idx] *expc::qqbqed[idx];

	  if (max(opts.order_hcoef,opts.order_evol) == 1)
	    continue;
	  
	  qqb[idx]  += qqb_eqg[idx] *expc::qg[idx];
	  qg[idx]   += qg_eqg[idx]  *expc::qg[idx];
	  qq[idx]   = qq_eqqb[idx]  *expc::qqb[idx] + qq_eqg[idx]*expc::qg[idx];
	  qqp[idx]  = qqp_eqqb[idx] *expc::qqb[idx] + qqp_eqg[idx]*expc::qg[idx];
	  qqbp[idx] = qqbp_eqqb[idx]*expc::qqb[idx] + qqbp_eqg[idx]*expc::qg[idx];
	  gg[idx]   = gg_eqqb[idx]  *expc::qqb[idx] + gg_eqg[idx]*expc::qg[idx];

	  if (max(opts.order_hcoef,opts.order_evol) == 2)
	    continue;

	  qg[idx]   += qg_eqq[idx]  *expc::qq[idx] + qg_eqqp[idx]*expc::qqp[idx] + qg_eqqbp[idx]*expc::qqbp[idx] + qg_egg[idx]*expc::gg[idx];
	  qq[idx]   += qq_eqq[idx]  *expc::qq[idx];
	  qqp[idx]  += qqp_eqqp[idx]*expc::qqp[idx];
	  qqbp[idx] += qqbp_eqqbp[idx]*expc::qqbp[idx];
	  gg[idx]   += gg_egg[idx]  *expc::gg[idx];
	  qbg[idx]  = qbg_eqqb[idx] *expc::qqb[idx] + qbg_eqg[idx]*expc::qg[idx] + qbg_eqq[idx]*expc::qq[idx] + qbg_egg[idx]*expc::gg[idx];
	  qpg[idx]  = qpg_eqqb[idx] *expc::qqb[idx] + qpg_eqg[idx]*expc::qg[idx] + qpg_eqqbp[idx]*expc::qqbp[idx] + qpg_egg[idx]*expc::gg[idx];
	  qbpg[idx] = qbpg_eqqb[idx]*expc::qqb[idx] + qbpg_eqg[idx]*expc::qg[idx] + qbpg_eqqp[idx]*expc::qqp[idx] + qbpg_egg[idx]*expc::gg[idx];

	  if (max(opts.order_hcoef,opts.order_evol) == 3)
	    continue;

	  qqb[idx]  += qqb_egg[idx] *expc::gg[idx];
	  qq[idx]   += qq_egg[idx]  *expc::gg[idx] + qq_eqbg[idx]  *expc::qbg[idx];
	  qqp[idx]  += qqp_egg[idx] *expc::gg[idx] + qqp_eqbpg[idx]*expc::qbpg[idx];
	  qqbp[idx] += qqbp_egg[idx]*expc::gg[idx] + qqbp_eqpg[idx]*expc::qpg[idx];
	  gg[idx]   += gg_eqq[idx]  *expc::qq[idx] + gg_eqqp[idx]  *expc::qqp[idx] + gg_eqqbp[idx]  *expc::qqbp[idx] + gg_eqbg[idx]*expc::qbg[idx] + gg_eqpg[idx]*expc::qpg[idx] + gg_eqbpg[idx]*expc::qbpg[idx];
	  qbg[idx]  += qbg_eqbg[idx]  *expc::qbg[idx];
	  qpg[idx]  += qpg_eqpg[idx]  *expc::qpg[idx];
	  qbpg[idx] += qbpg_eqbpg[idx]*expc::qbpg[idx];
	  qbq[idx]   = qbq_eqqb[idx]   *expc::qqb[idx] + qbq_eqg[idx]   *expc::qg[idx] + qbq_eqq[idx]    *expc::qq[idx]   + qbq_egg[idx]   *expc::gg[idx] + qbq_eqbg[idx] *expc::qbg[idx];
	  qpqbp[idx] = qpqbp_eqqb[idx] *expc::qqb[idx] + qpqbp_eqg[idx] *expc::qg[idx] + qpqbp_eqqbp[idx]*expc::qqbp[idx] + qpqbp_egg[idx] *expc::gg[idx] + qpqbp_eqpg[idx]*expc::qpg[idx];
	  qbpqp[idx] = qbpqp_eqqb[idx] *expc::qqb[idx] + qbpqp_eqg[idx] *expc::qg[idx] + qbpqp_eqqp[idx] *expc::qqp[idx]  + qbpqp_egg[idx] *expc::gg[idx] + qbpqp_eqbpg[idx]*expc::qbpg[idx];
	  qbqp[idx]  = qbqp_eqqb[idx]  *expc::qqb[idx] + qbqp_eqg[idx]  *expc::qg[idx] + qbqp_eqq[idx]   *expc::qbqp[idx] + qbqp_eqqp[idx]  *expc::qqp[idx]  + qbqp_egg[idx]  *expc::gg[idx] + qbqp_eqbg[idx] *expc::qbg[idx] + qbqp_eqbpg[idx]*expc::qbpg[idx];
	  qbqbp[idx] = qbqbp_eqqb[idx] *expc::qqb[idx] + qbqbp_eqg[idx] *expc::qg[idx] + qbqbp_eqq[idx]  *expc::qbqbp[idx]+ qbqbp_eqqbp[idx]*expc::qqbp[idx] + qbqbp_egg[idx] *expc::gg[idx] + qbqbp_eqbg[idx]*expc::qbg[idx] + qbqbp_eqpg[idx]*expc::qpg[idx];
	  qpqp[idx]  = qpqp_eqqb[idx]  *expc::qqb[idx] + qpqp_eqg[idx]  *expc::qg[idx] + qpqp_eqqp[idx]  *expc::qqp[idx]  + qpqp_eqqbp[idx] *expc::qqbp[idx] + qpqp_egg[idx]  *expc::gg[idx] + qpqp_eqpg[idx] *expc::qpg[idx] + qpqp_eqbpg[idx]*expc::qbpg[idx];
	}  
  else //rapidity dependent
    for (int sign = mesq::positive; sign <= mesq::negative; sign++)
      for (int i1 = 0; i1 < mellinint::mdim; i1++)
	for (int i2 = 0; i2 < mellinint::mdim; i2++)
	  {
	    int idx = hcoeff::index(i1,i2,sign);

	    qqb[idx]   = qqb_eqqb[idx] *expc::qqb[idx];
	    qg_1[idx]  = qg_1_eqqb[idx]*expc::qqb[idx];
	    qg_2[idx]  = qg_2_eqqb[idx]*expc::qqb[idx];
	    qqbqed[idx]   = qqb_eqqbqed[idx] *expc::qqbqed[idx];
	    qgamma_1[idx]  = qgamma_1_eqqb[idx]*expc::qqbqed[idx];
	    qgamma_2[idx]  = qgamma_2_eqqb[idx]*expc::qqbqed[idx];
	    
	    if (max(opts.order_hcoef,opts.order_evol) == 1)
	      continue;

	    qqb[idx]    += qqb_eqg_1[idx] *expc::qg_1[idx] + qqb_eqg_2[idx]*expc::qg_2[idx];
	    qg_1[idx]   += qg_1_eqg_1[idx]*expc::qg_1[idx];
	    qg_2[idx]   += qg_2_eqg_2[idx]*expc::qg_2[idx];
	    qq_1[idx]   = qq_1_eqqb[idx]  *expc::qqb[idx] + qq_1_eqg_1[idx]*expc::qg_1[idx]  ;
	    qq_2[idx]   = qq_2_eqqb[idx]  *expc::qqb[idx] + qq_2_eqg_2[idx]*expc::qg_2[idx]  ;
	    qqp_1[idx]  = qqp_1_eqqb[idx] *expc::qqb[idx] + qqp_1_eqg_1[idx]*expc::qg_1[idx] ;
	    qqp_2[idx]  = qqp_2_eqqb[idx] *expc::qqb[idx] + qqp_2_eqg_2[idx]*expc::qg_2[idx] ;
	    qqbp_1[idx] = qqbp_1_eqqb[idx]*expc::qqb[idx] + qqbp_1_eqg_1[idx]*expc::qg_1[idx];
	    qqbp_2[idx] = qqbp_2_eqqb[idx]*expc::qqb[idx] + qqbp_2_eqg_2[idx]*expc::qg_2[idx];
	    gg[idx]     = gg_eqqb[idx]    *expc::qqb[idx] + gg_eqg_1[idx]*expc::qg_1[idx] + gg_eqg_2[idx]*expc::qg_2[idx];

	    if (max(opts.order_hcoef,opts.order_evol) == 2)
	      continue;

	    qg_1[idx]   += qg_1_eqg_2[idx]*expc::qg_2[idx] + qg_1_eqq_1[idx]*expc::qq_1[idx] + qg_1_eqqp_1[idx]*expc::qqp_1[idx] + qg_1_eqqbp_1[idx]*expc::qqbp_1[idx] + qg_1_egg[idx]*expc::gg[idx];
	    qg_2[idx]   += qg_2_eqg_1[idx]*expc::qg_1[idx] + qg_2_eqq_2[idx]*expc::qq_2[idx] + qg_2_eqqp_2[idx]*expc::qqp_2[idx] + qg_2_eqqbp_2[idx]*expc::qqbp_2[idx] + qg_2_egg[idx]*expc::gg[idx];
	    qq_1[idx]   += qq_1_eqq_1[idx]*expc::qq_1[idx];
	    qq_2[idx]   += qq_2_eqq_2[idx]*expc::qq_2[idx];
	    qqp_1[idx]  += qqp_1_eqqp_1[idx]*expc::qqp_1[idx];
	    qqp_2[idx]  += qqp_2_eqqp_2[idx]*expc::qqp_2[idx];
	    qqbp_1[idx] += qqbp_1_eqqbp_1[idx]*expc::qqbp_1[idx];
	    qqbp_2[idx] += qqbp_2_eqqbp_2[idx]*expc::qqbp_2[idx];
	    gg[idx]     += gg_egg[idx]*expc::gg[idx];
	    qbg_1[idx]  = qbg_1_eqqb[idx]  *expc::qqb[idx]+ qbg_1_eqg_1[idx]*expc::qg_1[idx] + qbg_1_eqg_2[idx]*expc::qg_2[idx] + qbg_1_eqq_2[idx]*expc::qq_2[idx] + qbg_1_egg[idx]*expc::gg[idx];
	    qbg_2[idx]  = qbg_2_eqqb[idx]  *expc::qqb[idx]+ qbg_2_eqg_1[idx]*expc::qg_1[idx] + qbg_2_eqg_2[idx]*expc::qg_2[idx] + qbg_2_eqq_1[idx]*expc::qq_1[idx] + qbg_2_egg[idx]*expc::gg[idx];
	    qpg_1[idx]  = qpg_1_eqqb[idx]  *expc::qqb[idx]+ qpg_1_eqg_1[idx]*expc::qg_1[idx] + qpg_1_eqg_2[idx]*expc::qg_2[idx] + qpg_1_eqqbp_2[idx]*expc::qqbp_2[idx] + qpg_1_egg[idx]*expc::gg[idx];
	    qpg_2[idx]  = qpg_2_eqqb[idx]  *expc::qqb[idx]+ qpg_2_eqg_1[idx]*expc::qg_1[idx] + qpg_2_eqg_2[idx]*expc::qg_2[idx] + qpg_2_eqqbp_1[idx]*expc::qqbp_1[idx] + qpg_2_egg[idx]*expc::gg[idx];
	    qbpg_1[idx] = qbpg_1_eqqb[idx] *expc::qqb[idx]+ qbpg_1_eqg_1[idx]*expc::qg_1[idx] + qbpg_1_eqg_2[idx]*expc::qg_2[idx] + qbpg_1_eqqp_2[idx]*expc::qqp_2[idx] + qbpg_1_egg[idx]*expc::gg[idx];
	    qbpg_2[idx] = qbpg_2_eqqb[idx] *expc::qqb[idx]+ qbpg_2_eqg_1[idx]*expc::qg_1[idx] + qbpg_2_eqg_2[idx]*expc::qg_2[idx] + qbpg_2_eqqp_1[idx]*expc::qqp_1[idx] + qbpg_2_egg[idx]*expc::gg[idx];

	    //NNNNLO to be implemented
	  }
}


void muf::calc(complex <double> b)
{
  //LO
  if (opts.order_hcoef == 0)
    return;

  //if (opts.evolmode == 1 && opts.mufevol)
  //return;

  double LQ2 = pow(LQ,2);
  double LQ3 = pow(LQ,3);
  double LQ4 = pow(LQ,4);
  double LQ5 = pow(LQ,5);
  double LQ6 = pow(LQ,6);

  double LR2 = pow(LR,2);

  double LF2 = pow(LF,2);
  double LF3 = pow(LF,3);
  double LQF  = LF - LQ;
  double LQF2 = pow(LQF,2);
  double LQF3 = pow(LQF,3);
  /*
   * double LQb = LQ;
   * double LFb = LF;
   * if (opts.mufevol && real(b) > 0)
   *   {
   *     double m = phasespace::m;
   *     double mub = real(pdfevol::mub);
   *     double muf = scales::fac;
   *     double Q = scales::res;
   *     double m2tilde   = pow(m*mub,2)/(pow(m,2) + pow(mub,2));
   *     double Q2tilde   = pow(Q*mub,2)/(pow(Q,2) + pow(mub,2));
   *     double muf2tilde = pow(muf*mub,2)/(pow(muf,2) + pow(mub,2));
   *     LQb = log(m2tilde/Q2tilde);
   *     LFb = log(m2tilde/muf2tilde);
   *
   *     //LQb = log( pow(phasespace::m,2)/(pow(scales::res*real(pdfevol::mub),2)/(pow(scales::res,2) + pow(real(pdfevol::mub),2))));
   *     //LFb = log( pow(phasespace::m,2)/(pow(scales::fac*real(pdfevol::mub),2)/(pow(scales::fac,2) + pow(real(pdfevol::mub),2))));
}
//LQb=0.;
//LFb=0.;
double LQb2 = pow(LQb,2);
double LQb3 = pow(LQb,3);
double LQb4 = pow(LQb,3);
double LQb5 = pow(LQb,3);

double LFb2 = pow(LFb,2);
double LFb3 = pow(LFb,3);

double LQFb = LFb - LQb;
*/

  double as = resint::aass;
  double as2 = pow(resint::aass,2);
  double as3 = pow(resint::aass,3);

  //cout << "mub " << real(pdfevol::mub) << "  " << LQb << "  " << LFb << "  " << LQFb << endl;

  if (opts.mellin1d) //rapidity integrated
    for (int sign = mesq::positive; sign <= mesq::negative; sign++)
      for (int i = 0; i < mellinint::mdim; i++)
      {
        int idx = pmom::index(i,sign);

        //Include exponentiation in the C coefficients
        complex <double> eC1qq   = ccoeff::C1qq[idx];
        complex <double> eC2qq   = ccoeff::C2qq[idx];
        complex <double> eC1qg   = ccoeff::C1qg[idx]*expc::aexpqg[idx]/expc::aexpqq[idx];
        complex <double> eC2qg   = ccoeff::C2qg[idx]*expc::aexpqg[idx]/expc::aexpqq[idx];
        complex <double> eC2qqp  = ccoeff::C2qqp[idx]*expc::aexpqqp[idx]/expc::aexpqq[idx];
        complex <double> eC2qqb  = ccoeff::C2qqb[idx]*expc::aexpqqb[idx]/expc::aexpqq[idx];
        complex <double> eC2qqbp = ccoeff::C2qqbp[idx]*expc::aexpqqbp[idx]/expc::aexpqq[idx];

        //NLO
        complex <double> H1st_qqb = 2.*gamma1qq[idx]*LQF;
        complex <double> H1st_qg  = gamma1qg[idx]*LQF;

        qqb[idx] += as*H1st_qqb;
        qg[idx]  += as*H1st_qg;

        if (opts.order_hcoef == 1)
          continue;

        //NNLO
        complex <double> H2st_qqb  = LQF*(2.*eC1qg*gamma1gq[idx] + 2.*gamma2qq[idx] + gamma1gq[idx]*gamma1qg[idx]*LQF + gamma1qq[idx]*(4.*eC1qq + 2.*gamma1qq[idx]*LF - 2.*(B1q + gamma1qq[idx])*LQ - A1q*LQ2 + beta0*(LF + LQ)));
        complex <double> H2st_qq   = ((2.*eC1qg*gamma1gq[idx] + 2.*gamma2qqb[idx] + gamma1gq[idx]*gamma1qg[idx]*LQF)*LQF)/2.;
        complex <double> H2st_qqp  = ((2.*eC1qg*gamma1gq[idx] + 2.*gamma2qqp[idx] + gamma1gq[idx]*gamma1qg[idx]*LQF)*LQF)/2.;
        complex <double> H2st_qqbp = ((2.*eC1qg*gamma1gq[idx] + 2.*gamma2qqp[idx] + gamma1gq[idx]*gamma1qg[idx]*LQF)*LQF)/2.;
        complex <double> H2st_qg   = (LQF*(4.*eC1qq*gamma1qg[idx] + 2.*eC1qg*(gamma1gg[idx] + gamma1qq[idx]) + 2.*gamma2qg[idx] + beta0*gamma1qg[idx]*LF + gamma1gg[idx]*gamma1qg[idx]*LF + 3.*gamma1qg[idx]*gamma1qq[idx]*LF - 2.*B1q*gamma1qg[idx]*LQ + beta0*gamma1qg[idx]*LQ - gamma1gg[idx]*gamma1qg[idx]*LQ - 3.*gamma1qg[idx]*gamma1qq[idx]*LQ - A1q*gamma1qg[idx]*LQ2))/2.;
        complex <double> H2st_gg   = gamma1qg[idx]*(2.*eC1qg + gamma1qg[idx]*LQF)*LQF;

        //Mixed LR*LQF variations
        H2st_qqb += -beta0*H1st_qqb*LR;
        H2st_qg  += -beta0*H1st_qg*LR;

        qqb[idx]  += as2*H2st_qqb;
        qg[idx]   += as2*H2st_qg;
        qq[idx]   += as2*H2st_qq;
        qqp[idx]  += as2*H2st_qqp;
        qqbp[idx] += as2*H2st_qqbp;
        gg[idx]   += as2*H2st_gg;

        if (opts.order_hcoef == 2)
          continue;

        //NNNLO
        double nf = 5.;
        complex <double> H3st_qqb  = (LQF*(24.*eC2qg*gamma1gq[idx] + 24.*pow(eC1qq,2)*gamma1qq[idx] + 48.*eC2qq*gamma1qq[idx] + 48.*eC1qq*gamma2qq[idx] + 24.*gamma3qq[idx] + 24.*eC1qq*gamma1gq[idx]*gamma1qg[idx]*LF + 12.*beta1*gamma1qq[idx]*LF + 24.*beta0*eC1qq*gamma1qq[idx]*LF + 48.*eC1qq*pow(gamma1qq[idx],2)*LF + 12.*gamma1qg[idx]*gamma2gq[idx]*LF + 12.*gamma1gq[idx]*gamma2qg[idx]*LF + 24.*beta0*gamma2qq[idx]*LF + 48.*gamma1qq[idx]*gamma2qq[idx]*LF + 12.*beta0*gamma1gq[idx]*gamma1qg[idx]*LF2 + 4.*gamma1gg[idx]*gamma1gq[idx]*gamma1qg[idx]*LF2 + 8.*pow(beta0,2)*gamma1qq[idx]*LF2 + 20.*gamma1gq[idx]*gamma1qg[idx]*gamma1qq[idx]*LF2 + 24.*beta0*pow(gamma1qq[idx],2)*LF2 + 16.*pow(gamma1qq[idx],3)*LF2 - 24.*eC1qq*gamma1gq[idx]*gamma1qg[idx]*LQ - 24.*B2q*gamma1qq[idx]*LQ + 12.*beta1*gamma1qq[idx]*LQ - 48.*B1q*eC1qq*gamma1qq[idx]*LQ + 72.*beta0*eC1qq*gamma1qq[idx]*LQ - 48.*eC1qq*pow(gamma1qq[idx],2)*LQ - 12.*gamma1qg[idx]*gamma2gq[idx]*LQ - 12.*gamma1gq[idx]*gamma2qg[idx]*LQ - 24.*B1q*gamma2qq[idx]*LQ + 24.*beta0*gamma2qq[idx]*LQ - 48.*gamma1qq[idx]*gamma2qq[idx]*LQ - 12.*B1q*gamma1gq[idx]*gamma1qg[idx]*LF*LQ - 8.*gamma1gg[idx]*gamma1gq[idx]*gamma1qg[idx]*LF*LQ - 12.*B1q*beta0*gamma1qq[idx]*LF*LQ + 8.*pow(beta0,2)*gamma1qq[idx]*LF*LQ - 40.*gamma1gq[idx]*gamma1qg[idx]*gamma1qq[idx]*LF*LQ - 24.*B1q*pow(gamma1qq[idx],2)*LF*LQ - 32.*pow(gamma1qq[idx],3)*LF*LQ + 12.*B1q*gamma1gq[idx]*gamma1qg[idx]*LQ2 - 12.*beta0*gamma1gq[idx]*gamma1qg[idx]*LQ2 + 4.*gamma1gg[idx]*gamma1gq[idx]*gamma1qg[idx]*LQ2 - 12.*A2q*gamma1qq[idx]*LQ2 + 12.*pow(B1q,2)*gamma1qq[idx]*LQ2 - 24.*B1q*beta0*gamma1qq[idx]*LQ2 + 8.*pow(beta0,2)*gamma1qq[idx]*LQ2 - 24.*A1q*eC1qq*gamma1qq[idx]*LQ2 + 20.*gamma1gq[idx]*gamma1qg[idx]*gamma1qq[idx]*LQ2 + 24.*B1q*pow(gamma1qq[idx],2)*LQ2 - 24.*beta0*pow(gamma1qq[idx],2)*LQ2 + 16.*pow(gamma1qq[idx],3)*LQ2 - 12.*A1q*gamma2qq[idx]*LQ2 - 6.*A1q*gamma1gq[idx]*gamma1qg[idx]*LF*LQ2 - 6.*A1q*beta0*gamma1qq[idx]*LF*LQ2 - 12.*A1q*pow(gamma1qq[idx],2)*LF*LQ2 + 6.*A1q*gamma1gq[idx]*gamma1qg[idx]*LQ3 + 12.*A1q*B1q*gamma1qq[idx]*LQ3 - 14.*A1q*beta0*gamma1qq[idx]*LQ3 + 12.*A1q*pow(gamma1qq[idx],2)*LQ3 + 3.*pow(A1q,2)*gamma1qq[idx]*LQ4 + 12.*eC1qg*(2.*eC1qq*gamma1gq[idx] + 2.*gamma2gq[idx] + gamma1gq[idx]*(beta0*LF + gamma1gg[idx]*LF + 3.*gamma1qq[idx]*LF - 2.*B1q*LQ + 3.*beta0*LQ - gamma1gg[idx]*LQ - 3.*gamma1qq[idx]*LQ - A1q*LQ2))))/12.;
        complex <double> H3st_qq   = (LQF*(12.*(eC2qg*gamma1gq[idx] + 2.*eC2qqb*gamma1qq[idx] + 2.*eC1qq*gamma2qqb[idx] + gamma3qqb[idx]) + 6.*(2.*eC1qq*gamma1gq[idx]*gamma1qg[idx] + gamma1qg[idx]*gamma2gq[idx] + gamma1gq[idx]*gamma2qg[idx] + 2.*beta0*gamma2qqb[idx] + 4.*gamma1qq[idx]*gamma2qqb[idx])*LF + 2.*gamma1gq[idx]*gamma1qg[idx]*(3.*beta0 + gamma1gg[idx] + 5.*gamma1qq[idx])*LF2 - 2.*(6.*eC1qq*gamma1gq[idx]*gamma1qg[idx] + 3.*gamma1qg[idx]*gamma2gq[idx] + 3.*gamma1gq[idx]*gamma2qg[idx] + 6.*(B1q - beta0 + 2.*gamma1qq[idx])*gamma2qqb[idx] + gamma1gq[idx]*gamma1qg[idx]*(3.*B1q + 2.*gamma1gg[idx] + 10.*gamma1qq[idx])*LF)*LQ + (-6.*A1q*gamma2qqb[idx] + gamma1gq[idx]*gamma1qg[idx]*(6.*B1q - 6.*beta0 + 2.*gamma1gg[idx] + 10.*gamma1qq[idx] - 3.*A1q*LF))*LQ2 + 3.*A1q*gamma1gq[idx]*gamma1qg[idx]*LQ3 + 6.*eC1qg*(2.*gamma2gq[idx] + gamma1gq[idx]*(2.*eC1qq + (beta0 + gamma1gg[idx] + 3.*gamma1qq[idx])*LF - (2.*B1q - 3.*beta0 + gamma1gg[idx] + 3.*gamma1qq[idx])*LQ - A1q*LQ2))))/12.;
        complex <double> H3st_qqp  = (LQF*(12.*eC2qg*gamma1gq[idx] + 24.*eC2qqbp*gamma1qq[idx] + 24.*eC1qq*gamma2qqbp[idx] + 12.*gamma3qqbp[idx] + 12.*eC1qq*gamma1gq[idx]*gamma1qg[idx]*LF + 6.*gamma1qg[idx]*gamma2gq[idx]*LF + 6.*gamma1gq[idx]*gamma2qg[idx]*LF + 12.*beta0*gamma2qqbp[idx]*LF + 24.*gamma1qq[idx]*gamma2qqbp[idx]*LF + 6.*beta0*gamma1gq[idx]*gamma1qg[idx]*LF2 + 2.*gamma1gg[idx]*gamma1gq[idx]*gamma1qg[idx]*LF2 + 10.*gamma1gq[idx]*gamma1qg[idx]*gamma1qq[idx]*LF2 - 12.*eC1qq*gamma1gq[idx]*gamma1qg[idx]*LQ - 6.*gamma1qg[idx]*gamma2gq[idx]*LQ - 6.*gamma1gq[idx]*gamma2qg[idx]*LQ - 12.*B1q*gamma2qqbp[idx]*LQ + 12.*beta0*gamma2qqbp[idx]*LQ - 24.*gamma1qq[idx]*gamma2qqbp[idx]*LQ - 6.*B1q*gamma1gq[idx]*gamma1qg[idx]*LF*LQ - 4.*gamma1gg[idx]*gamma1gq[idx]*gamma1qg[idx]*LF*LQ - 20.*gamma1gq[idx]*gamma1qg[idx]*gamma1qq[idx]*LF*LQ + 6.*B1q*gamma1gq[idx]*gamma1qg[idx]*LQ2 - 6.*beta0*gamma1gq[idx]*gamma1qg[idx]*LQ2 + 2.*gamma1gg[idx]*gamma1gq[idx]*gamma1qg[idx]*LQ2 + 10.*gamma1gq[idx]*gamma1qg[idx]*gamma1qq[idx]*LQ2 - 6.*A1q*gamma2qqbp[idx]*LQ2 - 3.*A1q*gamma1gq[idx]*gamma1qg[idx]*LF*LQ2 + 3.*A1q*gamma1gq[idx]*gamma1qg[idx]*LQ3 + 6.*eC1qg*(2.*eC1qq*gamma1gq[idx] + 2.*gamma2gq[idx] + gamma1gq[idx]*(beta0*LF + gamma1gg[idx]*LF + 3.*gamma1qq[idx]*LF - 2.*B1q*LQ + 3.*beta0*LQ - gamma1gg[idx]*LQ - 3.*gamma1qq[idx]*LQ - A1q*LQ2))))/12.;
        complex <double> H3st_qqbp = (LQF*(12.*(eC2qg*gamma1gq[idx] + 2.*eC2qqp*gamma1qq[idx] + 2.*eC1qq*gamma2qqp[idx] + gamma3qqp[idx]) + 6.*(2.*eC1qq*gamma1gq[idx]*gamma1qg[idx] + gamma1qg[idx]*gamma2gq[idx] + gamma1gq[idx]*gamma2qg[idx] + 2.*beta0*gamma2qqp[idx] + 4.*gamma1qq[idx]*gamma2qqp[idx])*LF + 2.*gamma1gq[idx]*gamma1qg[idx]*(3.*beta0 + gamma1gg[idx] + 5.*gamma1qq[idx])*LF2 - 2.*(6.*eC1qq*gamma1gq[idx]*gamma1qg[idx] + 3.*gamma1qg[idx]*gamma2gq[idx] + 3.*gamma1gq[idx]*gamma2qg[idx] + 6.*(B1q - beta0 + 2.*gamma1qq[idx])*gamma2qqp[idx] + gamma1gq[idx]*gamma1qg[idx]*(3.*B1q + 2.*gamma1gg[idx] + 10.*gamma1qq[idx])*LF)*LQ + (-6.*A1q*gamma2qqp[idx] + gamma1gq[idx]*gamma1qg[idx]*(6.*B1q - 6.*beta0 + 2.*gamma1gg[idx] + 10.*gamma1qq[idx] - 3.*A1q*LF))*LQ2 + 3.*A1q*gamma1gq[idx]*gamma1qg[idx]*LQ3 + 6.*eC1qg*(2.*gamma2gq[idx] + gamma1gq[idx]*(2.*eC1qq + (beta0 + gamma1gg[idx] + 3.*gamma1qq[idx])*LF - (2.*B1q - 3.*beta0 + gamma1gg[idx] + 3.*gamma1qq[idx])*LQ - A1q*LQ2))))/12.;
        complex <double> H3st_qg   = (LQF*(24.*pow(eC1qg,2)*gamma1gq[idx] + 24.*eC2qg*(gamma1gg[idx] + gamma1qq[idx]) + 24.*((pow(eC1qq,2) + 2.*eC2qq + eC2qqb - eC2qqbp - eC2qqp)*gamma1qg[idx] + 2.*eC1qq*gamma2qg[idx] + gamma3qg[idx]) + 12.*(beta1 + 2.*beta0*eC1qq)*gamma1qg[idx]*LF + 24.*eC1qq*gamma1gg[idx]*gamma1qg[idx]*LF + 72.*eC1qq*gamma1qg[idx]*gamma1qq[idx]*LF + 12.*gamma1qg[idx]*gamma2gg[idx]*LF + 24.*beta0*gamma2qg[idx]*LF + 12.*gamma1gg[idx]*gamma2qg[idx]*LF + 36.*gamma1qq[idx]*gamma2qg[idx]*LF + 36.*gamma1qg[idx]*gamma2qq[idx]*LF + 12.*gamma1qg[idx]*gamma2qqb[idx]*LF - 12.*gamma1qg[idx]*gamma2qqbp[idx]*LF - 12.*gamma1qg[idx]*gamma2qqp[idx]*LF + 8.*pow(beta0,2)*gamma1qg[idx]*LF2 + 12.*beta0*gamma1gg[idx]*gamma1qg[idx]*LF2 + 4.*pow(gamma1gg[idx],2)*gamma1qg[idx]*LF2 + 12.*gamma1gq[idx]*pow(gamma1qg[idx],2)*LF2 + 36.*beta0*gamma1qg[idx]*gamma1qq[idx]*LF2 + 16.*gamma1gg[idx]*gamma1qg[idx]*gamma1qq[idx]*LF2 + 28.*gamma1qg[idx]*pow(gamma1qq[idx],2)*LF2 - 24.*B2q*gamma1qg[idx]*LQ + 12.*beta1*gamma1qg[idx]*LQ - 48.*B1q*eC1qq*gamma1qg[idx]*LQ + 72.*beta0*eC1qq*gamma1qg[idx]*LQ - 24.*eC1qq*gamma1gg[idx]*gamma1qg[idx]*LQ - 72.*eC1qq*gamma1qg[idx]*gamma1qq[idx]*LQ - 12.*gamma1qg[idx]*gamma2gg[idx]*LQ - 24.*B1q*gamma2qg[idx]*LQ + 24.*beta0*gamma2qg[idx]*LQ - 12.*gamma1gg[idx]*gamma2qg[idx]*LQ - 36.*gamma1qq[idx]*gamma2qg[idx]*LQ - 36.*gamma1qg[idx]*gamma2qq[idx]*LQ - 12.*gamma1qg[idx]*gamma2qqb[idx]*LQ + 12.*gamma1qg[idx]*gamma2qqbp[idx]*LQ + 12.*gamma1qg[idx]*gamma2qqp[idx]*LQ - 12.*B1q*beta0*gamma1qg[idx]*LF*LQ + 8.*pow(beta0,2)*gamma1qg[idx]*LF*LQ - 12.*B1q*gamma1gg[idx]*gamma1qg[idx]*LF*LQ - 8.*pow(gamma1gg[idx],2)*gamma1qg[idx]*LF*LQ - 24.*gamma1gq[idx]*pow(gamma1qg[idx],2)*LF*LQ - 36.*B1q*gamma1qg[idx]*gamma1qq[idx]*LF*LQ - 32.*gamma1gg[idx]*gamma1qg[idx]*gamma1qq[idx]*LF*LQ - 56.*gamma1qg[idx]*pow(gamma1qq[idx],2)*LF*LQ - 12.*A2q*gamma1qg[idx]*LQ2 + 12.*pow(B1q,2)*gamma1qg[idx]*LQ2 - 24.*B1q*beta0*gamma1qg[idx]*LQ2 + 8.*pow(beta0,2)*gamma1qg[idx]*LQ2 - 24.*A1q*eC1qq*gamma1qg[idx]*LQ2 + 12.*B1q*gamma1gg[idx]*gamma1qg[idx]*LQ2 - 12.*beta0*gamma1gg[idx]*gamma1qg[idx]*LQ2 + 4.*pow(gamma1gg[idx],2)*gamma1qg[idx]*LQ2 + 12.*gamma1gq[idx]*pow(gamma1qg[idx],2)*LQ2 + 36.*B1q*gamma1qg[idx]*gamma1qq[idx]*LQ2 - 36.*beta0*gamma1qg[idx]*gamma1qq[idx]*LQ2 + 16.*gamma1gg[idx]*gamma1qg[idx]*gamma1qq[idx]*LQ2 + 28.*gamma1qg[idx]*pow(gamma1qq[idx],2)*LQ2 - 12.*A1q*gamma2qg[idx]*LQ2 - 6.*A1q*beta0*gamma1qg[idx]*LF*LQ2 - 6.*A1q*gamma1gg[idx]*gamma1qg[idx]*LF*LQ2 - 18.*A1q*gamma1qg[idx]*gamma1qq[idx]*LF*LQ2 + 12.*A1q*B1q*gamma1qg[idx]*LQ3 - 14.*A1q*beta0*gamma1qg[idx]*LQ3 + 6.*A1q*gamma1gg[idx]*gamma1qg[idx]*LQ3 + 18.*A1q*gamma1qg[idx]*gamma1qq[idx]*LQ3 + 3.*pow(A1q,2)*gamma1qg[idx]*LQ4 + 4.*gamma1qg[idx]*(6.*eC2qqbp + 6.*eC2qqp + (3.*gamma2qqbp[idx] + 3.*gamma2qqp[idx] + 2.*gamma1gq[idx]*gamma1qg[idx]*LQF)*LQF)*nf + 12.*eC1qg*(2.*eC1qq*(gamma1gg[idx] + gamma1qq[idx]) + 2.*(gamma2gg[idx] + gamma2qq[idx]) + pow(gamma1gg[idx],2)*LF + 3.*gamma1gq[idx]*gamma1qg[idx]*LF + 2.*gamma1gg[idx]*gamma1qq[idx]*LF + pow(gamma1qq[idx],2)*LF - 2.*B1q*gamma1gg[idx]*LQ - pow(gamma1gg[idx],2)*LQ - 3.*gamma1gq[idx]*gamma1qg[idx]*LQ - 2.*B1q*gamma1qq[idx]*LQ - 2.*gamma1gg[idx]*gamma1qq[idx]*LQ - pow(gamma1qq[idx],2)*LQ - A1q*gamma1gg[idx]*LQ2 - A1q*gamma1qq[idx]*LQ2 + beta0*(gamma1gg[idx] + gamma1qq[idx])*(LF + 3.*LQ) + 2.*gamma1gq[idx]*gamma1qg[idx]*LQF*nf)))/24.;
        complex <double> H3st_qbg  = ((2.*eC2qqb*gamma1qg[idx] + (eC1qg + gamma1qg[idx]*LQF)*(2.*eC1qg*gamma1gq[idx] + 2.*gamma2qqb[idx] + gamma1gq[idx]*gamma1qg[idx]*LQF))*LQF)/2.;
        complex <double> H3st_qpg  = ((2.*eC2qqp*gamma1qg[idx] + (eC1qg + gamma1qg[idx]*LQF)*(2.*eC1qg*gamma1gq[idx] + 2.*gamma2qqp[idx] + gamma1gq[idx]*gamma1qg[idx]*LQF))*LQF)/2.;
        complex <double> H3st_qbpg = ((2.*eC2qqbp*gamma1qg[idx] + (eC1qg + gamma1qg[idx]*LQF)*(2.*eC1qg*gamma1gq[idx] + 2.*gamma2qqbp[idx] + gamma1gq[idx]*gamma1qg[idx]*LQF))*LQF)/2.;
        complex <double> H3st_gg   = (LQF*(4.*pow(eC1qg,2)*gamma1gg[idx] + 4.*eC2qg*gamma1qg[idx] + 2.*eC1qg*(2.*eC1qq*gamma1qg[idx] + 2.*gamma2qg[idx] + gamma1qg[idx]*(beta0*LF + 3.*gamma1gg[idx]*LF + gamma1qq[idx]*LF - 2.*B1q*LQ + 3.*beta0*LQ - 3.*gamma1gg[idx]*LQ - gamma1qq[idx]*LQ - A1q*LQ2)) + gamma1qg[idx]*LQF*(4.*eC1qq*gamma1qg[idx] + 4.*gamma2qg[idx] + gamma1qg[idx]*(2.*gamma1gg[idx]*LF + 2.*gamma1qq[idx]*LF - 2.*B1q*LQ - 2.*gamma1gg[idx]*LQ - 2.*gamma1qq[idx]*LQ - A1q*LQ2 + 2.*beta0*(LF + LQ)))))/2.;

        //Mixed LR*LQF variations
        H3st_qqb  += -2*beta0*H2st_qqb*LR + (-beta1*LR - pow(beta0,2)*pow(LR,2))*H1st_qqb;
        H3st_qg   += -2*beta0*H2st_qg*LR  + (-beta1*LR - pow(beta0,2)*pow(LR,2))*H1st_qg;
        H3st_qq   += -2*beta0*H2st_qq*LR;
        H3st_qqp  += -2*beta0*H2st_qqp*LR;
        H3st_qqbp += -2*beta0*H2st_qqbp*LR;
        H3st_gg   += -2*beta0*H2st_gg*LR;

        qqb[idx]  += as3*H3st_qqb;
        qg[idx]   += as3*H3st_qg;
        qq[idx]   += as3*H3st_qq;
        qqp[idx]  += as3*H3st_qqp;
        qqbp[idx] += as3*H3st_qqbp;
        gg[idx]   += as3*H3st_gg;
        qbg[idx]  += as3*H3st_qbg;
        qpg[idx]  += as3*H3st_qpg;
        qbpg[idx] += as3*H3st_qbpg;
      }
      else //rapidity dependent
        for (int sign = mesq::positive; sign <= mesq::negative; sign++)
          for (int i1 = 0; i1 < mellinint::mdim; i1++)
            for (int i2 = 0; i2 < mellinint::mdim; i2++)
            {
              int ii1 = anomalous::index(i1,mesq::positive);
              int ii2 = anomalous::index(i2,sign);
              int idx = hcoeff::index(i1,i2,sign);

              //Include exponentiation in the C coefficients
              complex <double> eC1qq_1   = ccoeff::C1qq_1[ii1];
              complex <double> eC2qq_1   = ccoeff::C2qq_1[ii1];
              complex <double> eC1qg_1   = ccoeff::C1qg_1[ii1]*expc::aexpqg_1[ii1]/expc::aexpqq_1[ii1];
              complex <double> eC2qg_1   = ccoeff::C2qg_1[ii1]*expc::aexpqg_1[ii1]/expc::aexpqq_1[ii1];
              complex <double> eC2qqp_1  = ccoeff::C2qqp_1[ii1]*expc::aexpqqp_1[ii1]/expc::aexpqq_1[ii1];
              complex <double> eC2qqb_1  = ccoeff::C2qqb_1[ii1]*expc::aexpqqb_1[ii1]/expc::aexpqq_1[ii1];
              complex <double> eC2qqbp_1 = ccoeff::C2qqbp_1[ii1]*expc::aexpqqbp_1[ii1]/expc::aexpqq_1[ii1];
              complex <double> eC1qq_2   = ccoeff::C1qq_2[ii2];
              complex <double> eC2qq_2   = ccoeff::C2qq_2[ii2];
              complex <double> eC1qg_2   = ccoeff::C1qg_2[ii2]*expc::aexpqg_2[ii2]/expc::aexpqq_2[ii2];
              complex <double> eC2qg_2   = ccoeff::C2qg_2[ii2]*expc::aexpqg_2[ii2]/expc::aexpqq_2[ii2];
              complex <double> eC2qqp_2  = ccoeff::C2qqp_2[ii2]*expc::aexpqqp_2[ii2]/expc::aexpqq_2[ii2];
              complex <double> eC2qqb_2  = ccoeff::C2qqb_2[ii2]*expc::aexpqqb_2[ii2]/expc::aexpqq_2[ii2];
              complex <double> eC2qqbp_2 = ccoeff::C2qqbp_2[ii2]*expc::aexpqqbp_2[ii2]/expc::aexpqq_2[ii2];


              complex <double> H1st_qqb  = (gamma1qq_1[ii1] + gamma1qq_2[ii2])*LQF ;
              complex <double> H1st_qg_1 = gamma1qg_1[ii1]*LQF;
              complex <double> H1st_qg_2 = gamma1qg_2[ii2]*LQF;

              qqb[idx]   += as*H1st_qqb;
              qg_1[idx]  += as*H1st_qg_1;
              qg_2[idx]  += as*H1st_qg_2;

              if (opts.order_hcoef == 1)
                continue;

              complex <double> H2st_qqb    = LQF*(eC1qg_1*gamma1gq_1[ii1] + eC1qq_2*gamma1qq_1[ii1] + eC1qg_2*gamma1gq_2[ii2] + eC1qq_2*gamma1qq_2[ii2] + eC1qq_1*(gamma1qq_1[ii1] + gamma1qq_2[ii2]) + gamma2qq_1[ii1] + gamma2qq_2[ii2] + ((gamma1gq_1[ii1]*gamma1qg_1[ii1] + gamma1gq_2[ii2]*gamma1qg_2[ii2] + (gamma1qq_1[ii1] + gamma1qq_2[ii2])*(beta0 + gamma1qq_1[ii1] + gamma1qq_2[ii2]))*LQF)/ 2. - ((gamma1qq_1[ii1] + gamma1qq_2[ii2])*LQ*(2.*B1q - 2.*beta0 + A1q*LQ))/2.);
              complex <double> H2st_qg_1   = (LQF*(2.*eC1qg_1*(gamma1qq_2[ii2] + gamma1gg_1[ii1]) + 2.*gamma2qg_1[ii1] + gamma1qg_1[ii1]*(2.*eC1qq_2 + 2.*eC1qq_1 + (beta0 + 2.*gamma1qq_2[ii2] + gamma1gg_1[ii1] + gamma1qq_1[ii1])*LF - (2.*B1q - beta0 + 2.*gamma1qq_2[ii2] + gamma1gg_1[ii1] + gamma1qq_1[ii1])*LQ - A1q*LQ2)))/2.;
              complex <double> H2st_qg_2   = (LQF*(2.*eC1qg_2*(gamma1qq_1[ii1] + gamma1gg_2[ii2]) + 2.*gamma2qg_2[ii2] + gamma1qg_2[ii2]*(2.*eC1qq_1 + 2.*eC1qq_2 + (beta0 + 2.*gamma1qq_1[ii1] + gamma1gg_2[ii2] + gamma1qq_2[ii2])*LF - (2.*B1q - beta0 + 2.*gamma1qq_1[ii1] + gamma1gg_2[ii2] + gamma1qq_2[ii2])*LQ - A1q*LQ2)))/2.;
              complex <double> H2st_qq_1   = ((2.*eC1qg_1*gamma1gq_1[ii1] + 2.*gamma2qqb_1[ii1] + gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQF)* LQF)/2.;
              complex <double> H2st_qq_2   = ((2.*eC1qg_2*gamma1gq_2[ii2] + 2.*gamma2qqb_2[ii2] + gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQF)* LQF)/2.;
              complex <double> H2st_qqp_1  = ((2.*eC1qg_1*gamma1gq_1[ii1] + 2.*gamma2qqbp_1[ii1] + gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQF)* LQF)/2.;
              complex <double> H2st_qqp_2  = ((2.*eC1qg_2*gamma1gq_2[ii2] + 2.*gamma2qqbp_2[ii2] + gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQF)* LQF)/2.;
              complex <double> H2st_qqbp_1 = ((2.*eC1qg_1*gamma1gq_1[ii1] + 2.*gamma2qqp_1[ii1] + gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQF)* LQF)/2.;
              complex <double> H2st_qqbp_2 = ((2.*eC1qg_2*gamma1gq_2[ii2] + 2.*gamma2qqp_2[ii2] + gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQF)* LQF)/2.;
              complex <double> H2st_gg     = (eC1qg_2*gamma1qg_1[ii1] + gamma1qg_2[ii2]*(eC1qg_1 + gamma1qg_1[ii1]*LQF))*LQF;

              //Mixed LR*LQF variations
              H2st_qqb  += - beta0*H1st_qqb*LR;
              H2st_qg_1 += - beta0*H1st_qg_1*LR;
              H2st_qg_2 += - beta0*H1st_qg_2*LR;

              qqb[idx]    += as2*H2st_qqb;
              qg_1[idx]   += as2*H2st_qg_1;
              qg_2[idx]   += as2*H2st_qg_2;
              qq_1[idx]   += as2*H2st_qq_1;
              qq_2[idx]   += as2*H2st_qq_2;
              qqp_1[idx]  += as2*H2st_qqp_1;
              qqp_2[idx]  += as2*H2st_qqp_2;
              qqbp_1[idx] += as2*H2st_qqbp_1;
              qqbp_2[idx] += as2*H2st_qqbp_2;
              gg[idx]     += as2*H2st_gg;

              if (opts.order_hcoef == 2)
                continue;

              double nf = 5.;
              //bug
              //complex <double> H3st_qqb    = eC1qg_2*eC1qq_1*gamma1gq_2[ii2]*LF + eC2qg_2*gamma1gq_2[ii2]*LF + eC1qq_2*eC1qq_1*gamma1qq_2[ii2]*LF + eC2qq_2*gamma1qq_2[ii2]*LF + eC2qq_1*gamma1qq_2[ii2]*LF + eC1qq_2*eC1qg_1*gamma1gq_1[ii1]*LF + eC2qg_1*gamma1gq_1[ii1]*LF + eC1qq_2*eC1qq_1*gamma1qq_1[ii1]*LF + eC2qq_2*gamma1qq_1[ii1]*LF + eC2qq_1*gamma1qq_1[ii1]*LF + eC1qg_1*gamma2gq_1[ii1]*LF + eC1qg_2*gamma2gq_2[ii2]*LF + eC1qq_2*gamma2qq_2[ii2]*LF + eC1qq_1*gamma2qq_2[ii2]*LF + eC1qq_2*gamma2qq_1[ii1]*LF + eC1qq_1*gamma2qq_1[ii1]*LF + gamma3qq_2[ii2]*LF + gamma3qq_1[ii1]*LF + (beta0*eC1qg_2*gamma1gq_2[ii2]*LF2)/2. + (eC1qg_2*gamma1gg_2[ii2]*gamma1gq_2[ii2]*LF2)/2. + (eC1qq_2*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF2)/ 2. + (eC1qq_1*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF2)/2. + (beta1*gamma1qq_2[ii2]*LF2)/2. + (beta0*eC1qq_2*gamma1qq_2[ii2]*LF2)/2. + (beta0*eC1qq_1*gamma1qq_2[ii2]*LF2)/2. + (eC1qg_2*gamma1gq_2[ii2]*gamma1qq_2[ii2]*LF2)/2. + (eC1qq_2*pow(gamma1qq_2[ii2],2)*LF2)/2. + (eC1qq_1*pow(gamma1qq_2[ii2],2)*LF2)/2. + (beta0*eC1qg_1*gamma1gq_1[ii1]*LF2)/2. + eC1qg_1*gamma1qq_2[ii2]*gamma1gq_1[ii1]*LF2 + (eC1qq_2*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF2)/2. + (beta1*gamma1qq_1[ii1]*LF2)/2. + (beta0*eC1qq_2*gamma1qq_1[ii1]*LF2)/2. + (beta0*eC1qq_1*gamma1qq_1[ii1]*LF2)/2. + eC1qg_2*gamma1gq_2[ii2]*gamma1qq_1[ii1]*LF2 + eC1qq_2*gamma1qq_2[ii2]*gamma1qq_1[ii1]*LF2 + eC1qq_1*gamma1qq_2[ii2]*gamma1qq_1[ii1]*LF2 + (eC1qq_2*pow(gamma1qq_1[ii1],2)*LF2)/2. + gamma1qg_2[ii2]*gamma2gq_2[ii2]*LF2 + beta0*gamma2qq_2[ii2]*LF2 + gamma1qq_2[ii2]*gamma2qq_2[ii2]*LF2 + gamma1qq_1[ii1]*gamma2qq_2[ii2]*LF2 + gamma1qg_1[ii1]*gamma2gq_1[ii1]*LF2 + beta0*gamma2qq_1[ii1]*LF2 + gamma1qq_2[ii2]*gamma2qq_1[ii1]*LF2 + gamma1qq_1[ii1]*gamma2qq_1[ii1]*LF2 + (beta0*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF3)/2. + (pow(beta0,2)*gamma1qq_2[ii2]*LF3)/3. + (beta0*pow(gamma1qq_2[ii2],2)*LF3)/2. + (beta0*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF3)/2. + (gamma1qq_2[ii2]*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF3)/2. + (pow(beta0,2)*gamma1qq_1[ii1]*LF3)/3. + (gamma1gq_2[ii2]*gamma1qg_2[ii2]*gamma1qq_1[ii1]*LF3)/2. + beta0*gamma1qq_2[ii2]*gamma1qq_1[ii1]* LF3 + (pow(gamma1qq_2[ii2],2)*gamma1qq_1[ii1]*LF3)/2. + (beta0*pow(gamma1qq_1[ii1],2)*LF3)/2. + (gamma1qq_2[ii2]*pow(gamma1qq_1[ii1],2)*LF3)/2. + ((gamma1gg_1[ii1]*gamma1gq_1[ii1]*gamma1qg_1[ii1] + 2.*gamma1gq_1[ii1]*gamma1qg_1[ii1]*gamma1qq_1[ii1] + pow(gamma1qq_1[ii1],3))*LF3)/6. + ((gamma1gg_2[ii2]*gamma1gq_2[ii2]*gamma1qg_2[ii2] + 2.*gamma1gq_2[ii2]*gamma1qg_2[ii2]*gamma1qq_2[ii2] + pow(gamma1qq_2[ii2],3))*LQF3)/6. - eC1qg_2*eC1qq_1*gamma1gq_2[ii2]*LQ - eC2qg_2*gamma1gq_2[ii2]*LQ - eC1qq_2*eC1qq_1*gamma1qq_2[ii2]*LQ - eC2qq_2*gamma1qq_2[ii2]*LQ - eC2qq_1*gamma1qq_2[ii2]*LQ - eC1qq_2*eC1qg_1*gamma1gq_1[ii1]*LQ - eC2qg_1*gamma1gq_1[ii1]*LQ - eC1qq_2*eC1qq_1*gamma1qq_1[ii1]*LQ - eC2qq_2*gamma1qq_1[ii1]*LQ - eC2qq_1*gamma1qq_1[ii1]*LQ - eC1qg_1*gamma2gq_1[ii1]*LQ - eC1qg_2*gamma2gq_2[ii2]*LQ - eC1qq_2*gamma2qq_2[ii2]*LQ - eC1qq_1*gamma2qq_2[ii2]*LQ - eC1qq_2*gamma2qq_1[ii1]*LQ - eC1qq_1*gamma2qq_1[ii1]*LQ - gamma3qq_2[ii2]*LQ - gamma3qq_1[ii1]*LQ - B1q*eC1qg_2*gamma1gq_2[ii2]*LF*LQ + beta0*eC1qg_2*gamma1gq_2[ii2]*LF*LQ - eC1qg_2*gamma1gg_2[ii2]*gamma1gq_2[ii2]*LF*LQ - eC1qq_2*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF*LQ - eC1qq_1*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF*LQ - B2q*gamma1qq_2[ii2]*LF*LQ - B1q*eC1qq_2*gamma1qq_2[ii2]*LF*LQ + beta0*eC1qq_2*gamma1qq_2[ii2]*LF*LQ - B1q*eC1qq_1*gamma1qq_2[ii2]*LF*LQ + beta0*eC1qq_1*gamma1qq_2[ii2]*LF*LQ - eC1qg_2*gamma1gq_2[ii2]*gamma1qq_2[ii2]*LF*LQ - eC1qq_2*pow(gamma1qq_2[ii2],2)*LF*LQ - eC1qq_1*pow(gamma1qq_2[ii2],2)*LF*LQ - B1q*eC1qg_1*gamma1gq_1[ii1]*LF*LQ + beta0*eC1qg_1*gamma1gq_1[ii1]*LF*LQ - 2.*eC1qg_1*gamma1qq_2[ii2]*gamma1gq_1[ii1]*LF*LQ - eC1qg_1*gamma1gg_1[ii1]*gamma1gq_1[ii1]*LF*LQ - eC1qq_2*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF*LQ - eC1qq_1*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF*LQ - B2q*gamma1qq_1[ii1]*LF*LQ - B1q*eC1qq_2*gamma1qq_1[ii1]*LF*LQ + beta0*eC1qq_2*gamma1qq_1[ii1]*LF*LQ - B1q*eC1qq_1*gamma1qq_1[ii1]*LF*LQ + beta0*eC1qq_1*gamma1qq_1[ii1]*LF*LQ - 2.*eC1qg_2*gamma1gq_2[ii2]*gamma1qq_1[ii1]*LF*LQ - 2.*eC1qq_2*gamma1qq_2[ii2]*gamma1qq_1[ii1]*LF*LQ - 2.*eC1qq_1*gamma1qq_2[ii2]*gamma1qq_1[ii1]*LF* LQ - eC1qg_1*gamma1gq_1[ii1]*gamma1qq_1[ii1]*LF*LQ - eC1qq_2*pow(gamma1qq_1[ii1],2)*LF*LQ - eC1qq_1*pow(gamma1qq_1[ii1],2)*LF*LQ - gamma1qg_2[ii2]*gamma2gq_2[ii2]*LF*LQ - gamma1gq_2[ii2]*gamma2qg_2[ii2]*LF*LQ - B1q*gamma2qq_2[ii2]*LF*LQ - 2.*gamma1qq_2[ii2]*gamma2qq_2[ii2]*LF*LQ - 2.*gamma1qq_1[ii1]*gamma2qq_2[ii2]*LF*LQ - gamma1qg_1[ii1]*gamma2gq_1[ii1]*LF*LQ - gamma1gq_1[ii1]*gamma2qg_1[ii1]*LF*LQ - B1q*gamma2qq_1[ii1]*LF*LQ - 2.*gamma1qq_2[ii2]*gamma2qq_1[ii1]*LF*LQ - 2.*gamma1qq_1[ii1]*gamma2qq_1[ii1]*LF*LQ - (B1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF2*LQ)/2. - (beta0*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF2*LQ)/2. - (B1q*beta0*gamma1qq_2[ii2]*LF2*LQ)/2. - (B1q*pow(gamma1qq_2[ii2],2)*LF2*LQ)/2. - (beta0*pow(gamma1qq_2[ii2],2)*LF2*LQ)/2. - (B1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF2*LQ)/2. - (beta0*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF2*LQ)/ 2. - (3.*gamma1qq_2[ii2]*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF2*LQ)/2. - (B1q*beta0*gamma1qq_1[ii1]*LF2*LQ)/2. - (3.*gamma1gq_2[ii2]*gamma1qg_2[ii2]*gamma1qq_1[ii1]*LF2*LQ)/2. - B1q*gamma1qq_2[ii2]*gamma1qq_1[ii1]*LF2*LQ - beta0*gamma1qq_2[ii2]*gamma1qq_1[ii1]*LF2*LQ - (3.*pow(gamma1qq_2[ii2],2)*gamma1qq_1[ii1]*LF2*LQ)/2. - (B1q*pow(gamma1qq_1[ii1],2)*LF2*LQ)/2. - (beta0*pow(gamma1qq_1[ii1],2)*LF2*LQ)/2. - (3.*gamma1qq_2[ii2]*pow(gamma1qq_1[ii1],2)*LF2*LQ)/2. - ((gamma1gg_1[ii1]*gamma1gq_1[ii1]*gamma1qg_1[ii1] + 2.*gamma1gq_1[ii1]*gamma1qg_1[ii1]*gamma1qq_1[ii1] + pow(gamma1qq_1[ii1],3))*LF2*LQ)/2. + B1q*eC1qg_2*gamma1gq_2[ii2]*LQ2 - (3.*beta0*eC1qg_2*gamma1gq_2[ii2]*LQ2)/2. + (eC1qg_2*gamma1gg_2[ii2]*gamma1gq_2[ii2]*LQ2)/2. + (eC1qq_2*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQ2)/2. + (eC1qq_1*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQ2)/ 2. + B2q*gamma1qq_2[ii2]*LQ2 - (beta1*gamma1qq_2[ii2]*LQ2)/2. + B1q*eC1qq_2*gamma1qq_2[ii2]*LQ2 - (3.*beta0*eC1qq_2*gamma1qq_2[ii2]*LQ2)/2. + B1q*eC1qq_1*gamma1qq_2[ii2]*LQ2 - (3.*beta0*eC1qq_1*gamma1qq_2[ii2]*LQ2)/2. + (eC1qg_2*gamma1gq_2[ii2]*gamma1qq_2[ii2]*LQ2)/2. + (eC1qq_2*pow(gamma1qq_2[ii2],2)*LQ2)/2. + (eC1qq_1*pow(gamma1qq_2[ii2],2)*LQ2)/2. + B1q*eC1qg_1*gamma1gq_1[ii1]*LQ2 - (3.*beta0*eC1qg_1*gamma1gq_1[ii1]*LQ2)/2. + eC1qg_1*gamma1qq_2[ii2]*gamma1gq_1[ii1]*LQ2 + (eC1qq_2*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQ2)/2. + B2q*gamma1qq_1[ii1]*LQ2 - (beta1*gamma1qq_1[ii1]*LQ2)/2. + B1q*eC1qq_2*gamma1qq_1[ii1]*LQ2 - (3.*beta0*eC1qq_2*gamma1qq_1[ii1]*LQ2)/2. + B1q*eC1qq_1*gamma1qq_1[ii1]*LQ2 - (3.*beta0*eC1qq_1*gamma1qq_1[ii1]*LQ2)/2. + eC1qg_2*gamma1gq_2[ii2]*gamma1qq_1[ii1]*LQ2 + eC1qq_2*gamma1qq_2[ii2]*gamma1qq_1[ii1]*LQ2 + eC1qq_1*gamma1qq_2[ii2]*gamma1qq_1[ii1]*LQ2 + (eC1qq_2*pow(gamma1qq_1[ii1],2)*LQ2)/2. + gamma1qg_2[ii2]*gamma2gq_2[ii2]*LQ2 + B1q*gamma2qq_2[ii2]*LQ2 - beta0*gamma2qq_2[ii2]*LQ2 + gamma1qq_2[ii2]*gamma2qq_2[ii2]*LQ2 + gamma1qq_1[ii1]*gamma2qq_2[ii2]*LQ2 + gamma1qg_1[ii1]*gamma2gq_1[ii1]*LQ2 + B1q*gamma2qq_1[ii1]*LQ2 - beta0*gamma2qq_1[ii1]*LQ2 + gamma1qq_2[ii2]*gamma2qq_1[ii1]*LQ2 + gamma1qq_1[ii1]*gamma2qq_1[ii1]*LQ2 - (A1q*eC1qg_2*gamma1gq_2[ii2]*LF*LQ2)/2. + B1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF*LQ2 - (beta0*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF*LQ2)/2. - (A2q*gamma1qq_2[ii2]*LF*LQ2)/2. + (pow(B1q,2)*gamma1qq_2[ii2]*LF*LQ2)/2. - (B1q*beta0*gamma1qq_2[ii2]*LF*LQ2)/2. - (A1q*eC1qq_2*gamma1qq_2[ii2]*LF*LQ2)/2. - (A1q*eC1qq_1*gamma1qq_2[ii2]*LF*LQ2)/2. + B1q*pow(gamma1qq_2[ii2],2)*LF*LQ2 - (beta0*pow(gamma1qq_2[ii2],2)*LF*LQ2)/2. - (A1q*eC1qg_1*gamma1gq_1[ii1]*LF*LQ2)/2. + B1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF*LQ2 - (beta0*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF*LQ2)/2. + (3.*gamma1qq_2[ii2]*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF*LQ2)/2. - (A2q*gamma1qq_1[ii1]*LF*LQ2)/2. + (pow(B1q,2)*gamma1qq_1[ii1]*LF*LQ2)/2. - (B1q*beta0*gamma1qq_1[ii1]*LF*LQ2)/2. - (A1q*eC1qq_2*gamma1qq_1[ii1]*LF*LQ2)/2. - (A1q*eC1qq_1*gamma1qq_1[ii1]*LF*LQ2)/2. + (3.*gamma1gq_2[ii2]*gamma1qg_2[ii2]*gamma1qq_1[ii1]*LF*LQ2)/2. + 2.*B1q*gamma1qq_2[ii2]*gamma1qq_1[ii1]*LF*LQ2 - beta0*gamma1qq_2[ii2]*gamma1qq_1[ii1]*LF*LQ2 + (3.*pow(gamma1qq_2[ii2],2)*gamma1qq_1[ii1]*LF*LQ2)/2. + B1q*pow(gamma1qq_1[ii1],2)*LF*LQ2 - (beta0*pow(gamma1qq_1[ii1],2)*LF*LQ2)/2. + (3.*gamma1qq_2[ii2]*pow(gamma1qq_1[ii1],2)*LF*LQ2)/2. + ((gamma1gg_1[ii1]*gamma1gq_1[ii1]*gamma1qg_1[ii1] + 2.*gamma1gq_1[ii1]*gamma1qg_1[ii1]*gamma1qq_1[ii1] + pow(gamma1qq_1[ii1],3))*LF*LQ2)/2. - (A1q*gamma2qq_2[ii2]*LF*LQ2)/2. - (A1q*gamma2qq_1[ii1]*LF*LQ2)/2. - (A1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF2*LQ2)/4. - (A1q*beta0*gamma1qq_2[ii2]*LF2*LQ2)/4. - (A1q*pow(gamma1qq_2[ii2],2)*LF2*LQ2)/4. - (A1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF2*LQ2)/4. - (A1q*beta0*gamma1qq_1[ii1]*LF2*LQ2)/4. - (A1q*gamma1qq_2[ii2]*gamma1qq_1[ii1]*LF2*LQ2)/2. - (A1q*pow(gamma1qq_1[ii1],2)*LF2*LQ2)/4. + (A1q*eC1qg_2*gamma1gq_2[ii2]*LQ3)/2. - (B1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQ3)/2. + (beta0*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQ3)/2. + (A2q*gamma1qq_2[ii2]*LQ3)/2. - (pow(B1q,2)*gamma1qq_2[ii2]*LQ3)/2. + B1q*beta0*gamma1qq_2[ii2]*LQ3 - (pow(beta0,2)*gamma1qq_2[ii2]*LQ3)/3. + (A1q*eC1qq_2*gamma1qq_2[ii2]*LQ3)/2. + (A1q*eC1qq_1*gamma1qq_2[ii2]*LQ3)/2. - (B1q*pow(gamma1qq_2[ii2],2)*LQ3)/2. + (beta0*pow(gamma1qq_2[ii2],2)*LQ3)/2. + (A1q*eC1qg_1*gamma1gq_1[ii1]*LQ3)/2. - (B1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQ3)/2. + (beta0*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQ3)/2. - (gamma1qq_2[ii2]*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQ3)/2. + (A2q*gamma1qq_1[ii1]*LQ3)/2. - (pow(B1q,2)*gamma1qq_1[ii1]*LQ3)/2. + B1q*beta0*gamma1qq_1[ii1]*LQ3 - (pow(beta0,2)*gamma1qq_1[ii1]*LQ3)/3. + (A1q*eC1qq_2*gamma1qq_1[ii1]*LQ3)/2. + (A1q*eC1qq_1*gamma1qq_1[ii1]*LQ3)/2. - (gamma1gq_2[ii2]*gamma1qg_2[ii2]*gamma1qq_1[ii1]*LQ3)/2. - B1q*gamma1qq_2[ii2]*gamma1qq_1[ii1]*LQ3 + beta0*gamma1qq_2[ii2]*gamma1qq_1[ii1]*LQ3 - (pow(gamma1qq_2[ii2],2)*gamma1qq_1[ii1]*LQ3)/2. - (B1q*pow(gamma1qq_1[ii1],2)*LQ3)/2. + (beta0*pow(gamma1qq_1[ii1],2)*LQ3)/2. - (gamma1qq_2[ii2]*pow(gamma1qq_1[ii1],2)*LQ3)/2. - ((gamma1gg_1[ii1]*gamma1gq_1[ii1]*gamma1qg_1[ii1] + 2.*gamma1gq_1[ii1]*gamma1qg_1[ii1]*gamma1qq_1[ii1] + pow(gamma1qq_1[ii1],3))*LQ3)/6. + (A1q*gamma2qq_2[ii2]*LQ3)/2. + (A1q*gamma2qq_1[ii1]*LQ3)/2. + (A1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF*LQ3)/2. + (A1q*B1q*gamma1qq_2[ii2]*LF*LQ3)/2. - (A1q*beta0*gamma1qq_2[ii2]*LF*LQ3)/3. + (A1q*pow(gamma1qq_2[ii2],2)*LF*LQ3)/2. + (A1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF*LQ3)/2. + (A1q*B1q*gamma1qq_1[ii1]*LF*LQ3)/2. - (A1q*beta0*gamma1qq_1[ii1]*LF*LQ3)/3. + A1q*gamma1qq_2[ii2]*gamma1qq_1[ii1]*LF*LQ3 + (A1q*pow(gamma1qq_1[ii1],2)*LF*LQ3)/2. - (A1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQ4)/4. - (A1q*B1q*gamma1qq_2[ii2]*LQ4)/2. + (7.*A1q*beta0*gamma1qq_2[ii2]*LQ4)/12. - (A1q*pow(gamma1qq_2[ii2],2)*LQ4)/4. - (A1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQ4)/4. - (A1q*B1q*gamma1qq_1[ii1]*LQ4)/2. + (7.*A1q*beta0*gamma1qq_1[ii1]*LQ4)/12. - (A1q*gamma1qq_2[ii2]*gamma1qq_1[ii1]*LQ4)/2. - (A1q*pow(gamma1qq_1[ii1],2)*LQ4)/4. + (pow(A1q,2)*gamma1qq_2[ii2]*LF*LQ4)/8. + (pow(A1q,2)*gamma1qq_1[ii1]*LF*LQ4)/8. - (pow(A1q,2)*gamma1qq_2[ii2]*LQ5)/8. - (pow(A1q,2)*gamma1qq_1[ii1]*LQ5)/8. + ((eC1qg_1*gamma1gq_1[ii1]*(gamma1gg_1[ii1] + gamma1qq_1[ii1]) + eC1qq_1*(gamma1gq_1[ii1]*gamma1qg_1[ii1] + pow(gamma1qq_1[ii1],2)))*(LF2 + LQ2))/2.;
              //complex <double> H3st_qg_1   = eC2qg_1*gamma1qq_2[ii2]*LF + eC1qg_1*(eC1qg_2*gamma1gq_2[ii2] + eC1qq_2*gamma1qq_2[ii2])*LF + eC2qq_2*gamma1qg_1[ii1]*LF + eC1qq_2*(eC1qg_1*gamma1gg_1[ii1] + eC1qq_1*gamma1qg_1[ii1])*LF + eC1qg_1*gamma2qq_2[ii2]*LF + eC1qq_2*gamma2qg_1[ii1]*LF + (eC1qg_1*gamma2gg_1[ii1] + eC1qq_1*gamma2qg_1[ii1])*LF + gamma3qg_1[ii1]*LF + (beta0*eC1qg_1*gamma1qq_2[ii2]*LF2)/2. + (eC1qg_1*(gamma1gq_2[ii2]*gamma1qg_2[ii2] + pow(gamma1qq_2[ii2],2))*LF2)/2. + (beta1*gamma1qg_1[ii1]*LF2)/2. + (beta0*eC1qq_2*gamma1qg_1[ii1]*LF2)/2. + (eC1qg_2*gamma1gq_2[ii2] + eC1qq_2*gamma1qq_2[ii2])*gamma1qg_1[ii1]*LF2 + (beta0*(eC1qg_1*gamma1gg_1[ii1] + eC1qq_1*gamma1qg_1[ii1])*LF2)/2. + gamma1qq_2[ii2]*(eC1qg_1*gamma1gg_1[ii1] + eC1qq_1*gamma1qg_1[ii1])*LF2 + (eC1qq_2*gamma1qg_1[ii1]*(gamma1gg_1[ii1] + gamma1qq_1[ii1])*LF2)/2. + gamma1qg_1[ii1]*gamma2qq_2[ii2]*LF2 + beta0*gamma2qg_1[ii1]*LF2 + gamma1qq_2[ii2]*gamma2qg_1[ii1]*LF2 + (gamma1qg_1[ii1]*gamma2gg_1[ii1] + gamma1qq_1[ii1]*gamma2qg_1[ii1])* LF2 + (pow(beta0,2)*gamma1qg_1[ii1]*LF3)/3. + beta0*gamma1qq_2[ii2]*gamma1qg_1[ii1]*LF3 + ((gamma1gq_2[ii2]*gamma1qg_2[ii2] + pow(gamma1qq_2[ii2],2))*gamma1qg_1[ii1]*LF3)/2. + (beta0*gamma1qg_1[ii1]*(gamma1gg_1[ii1] + gamma1qq_1[ii1])*LF3)/2. + (gamma1qq_2[ii2]*gamma1qg_1[ii1]*(gamma1gg_1[ii1] + gamma1qq_1[ii1])*LF3)/2. - eC2qg_1*gamma1qq_2[ii2]*LQ - eC1qg_1*(eC1qg_2*gamma1gq_2[ii2] + eC1qq_2*gamma1qq_2[ii2])*LQ - eC2qq_2*gamma1qg_1[ii1]*LQ - eC1qq_2*(eC1qg_1*gamma1gg_1[ii1] + eC1qq_1*gamma1qg_1[ii1])*LQ - eC1qg_1*gamma2qq_2[ii2]*LQ - eC1qq_2*gamma2qg_1[ii1]*LQ - (eC1qg_1*gamma2gg_1[ii1] + eC1qq_1*gamma2qg_1[ii1])*LQ - gamma3qg_1[ii1]*LQ - B1q*eC1qg_1*gamma1qq_2[ii2]*LF*LQ + beta0*eC1qg_1*gamma1qq_2[ii2]*LF*LQ - eC1qg_1*(gamma1gq_2[ii2]*gamma1qg_2[ii2] + pow(gamma1qq_2[ii2],2))*LF*LQ - B2q*gamma1qg_1[ii1]*LF*LQ - B1q*eC1qq_2*gamma1qg_1[ii1]*LF*LQ + beta0*eC1qq_2*gamma1qg_1[ii1]*LF*LQ - 2.*(eC1qg_2*gamma1gq_2[ii2] + eC1qq_2*gamma1qq_2[ii2])*gamma1qg_1[ii1]*LF*LQ - B1q*(eC1qg_1*gamma1gg_1[ii1] + eC1qq_1*gamma1qg_1[ii1])*LF*LQ + beta0*(eC1qg_1*gamma1gg_1[ii1] + eC1qq_1*gamma1qg_1[ii1])*LF*LQ - 2.*gamma1qq_2[ii2]*(eC1qg_1*gamma1gg_1[ii1] + eC1qq_1*gamma1qg_1[ii1])*LF*LQ - eC1qq_2*gamma1qg_1[ii1]*(gamma1gg_1[ii1] + gamma1qq_1[ii1])*LF*LQ - 2.*gamma1qg_1[ii1]*gamma2qq_2[ii2]*LF*LQ - B1q*gamma2qg_1[ii1]*LF*LQ - 2.*gamma1qq_2[ii2]*gamma2qg_1[ii1]*LF*LQ - (gamma1qg_1[ii1]*gamma2gg_1[ii1] + gamma1qq_1[ii1]*gamma2qg_1[ii1])*LF*LQ - (B1q*beta0*gamma1qg_1[ii1]*LF2*LQ)/2. - B1q*gamma1qq_2[ii2]*gamma1qg_1[ii1]*LF2*LQ - beta0*gamma1qq_2[ii2]*gamma1qg_1[ii1]*LF2*LQ - (3.*(gamma1gq_2[ii2]*gamma1qg_2[ii2] + pow(gamma1qq_2[ii2],2))*gamma1qg_1[ii1]*LF2*LQ)/2. - (B1q*gamma1qg_1[ii1]*(gamma1gg_1[ii1] + gamma1qq_1[ii1])*LF2*LQ)/2. - (beta0*gamma1qg_1[ii1]*(gamma1gg_1[ii1] + gamma1qq_1[ii1])*LF2*LQ)/2. - (3.*gamma1qq_2[ii2]*gamma1qg_1[ii1]*(gamma1gg_1[ii1] + gamma1qq_1[ii1])*LF2*LQ)/2. + B1q*eC1qg_1*gamma1qq_2[ii2]*LQ2 - (3.*beta0*eC1qg_1*gamma1qq_2[ii2]*LQ2)/2. + (eC1qg_1*(gamma1gq_2[ii2]*gamma1qg_2[ii2] + pow(gamma1qq_2[ii2],2))*LQ2)/2. + B2q*gamma1qg_1[ii1]*LQ2 - (beta1*gamma1qg_1[ii1]*LQ2)/2. + B1q*eC1qq_2*gamma1qg_1[ii1]*LQ2 - (3.*beta0*eC1qq_2*gamma1qg_1[ii1]*LQ2)/2. + (eC1qg_2*gamma1gq_2[ii2] + eC1qq_2*gamma1qq_2[ii2])* gamma1qg_1[ii1]*LQ2 + B1q*(eC1qg_1*gamma1gg_1[ii1] + eC1qq_1*gamma1qg_1[ii1])*LQ2 - (3.*beta0*(eC1qg_1*gamma1gg_1[ii1] + eC1qq_1*gamma1qg_1[ii1])*LQ2)/2. + gamma1qq_2[ii2]*(eC1qg_1*gamma1gg_1[ii1] + eC1qq_1*gamma1qg_1[ii1])*LQ2 + (eC1qq_2*gamma1qg_1[ii1]*(gamma1gg_1[ii1] + gamma1qq_1[ii1])*LQ2)/2. + gamma1qg_1[ii1]*gamma2qq_2[ii2]*LQ2 + B1q*gamma2qg_1[ii1]*LQ2 - beta0*gamma2qg_1[ii1]*LQ2 + gamma1qq_2[ii2]*gamma2qg_1[ii1]*LQ2 + (gamma1qg_1[ii1]*gamma2gg_1[ii1] + gamma1qq_1[ii1]*gamma2qg_1[ii1])* LQ2 - (A1q*eC1qg_1*gamma1qq_2[ii2]*LF*LQ2)/2. - (A2q*gamma1qg_1[ii1]*LF*LQ2)/2. + (pow(B1q,2)*gamma1qg_1[ii1]*LF*LQ2)/2. - (B1q*beta0*gamma1qg_1[ii1]*LF*LQ2)/2. - (A1q*eC1qq_2*gamma1qg_1[ii1]*LF*LQ2)/2. + 2.*B1q*gamma1qq_2[ii2]*gamma1qg_1[ii1]*LF*LQ2 - beta0*gamma1qq_2[ii2]*gamma1qg_1[ii1]*LF*LQ2 + (3.*(gamma1gq_2[ii2]*gamma1qg_2[ii2] + pow(gamma1qq_2[ii2],2))*gamma1qg_1[ii1]*LF*LQ2)/2. - (A1q*(eC1qg_1*gamma1gg_1[ii1] + eC1qq_1*gamma1qg_1[ii1])*LF*LQ2)/2. + B1q*gamma1qg_1[ii1]*(gamma1gg_1[ii1] + gamma1qq_1[ii1])*LF*LQ2 - (beta0*gamma1qg_1[ii1]*(gamma1gg_1[ii1] + gamma1qq_1[ii1])*LF*LQ2)/2. + (3.*gamma1qq_2[ii2]*gamma1qg_1[ii1]*(gamma1gg_1[ii1] + gamma1qq_1[ii1])*LF*LQ2)/2. - (A1q*gamma2qg_1[ii1]*LF*LQ2)/2. - (A1q*beta0*gamma1qg_1[ii1]*LF2*LQ2)/4. - (A1q*gamma1qq_2[ii2]*gamma1qg_1[ii1]*LF2*LQ2)/2. - (A1q*gamma1qg_1[ii1]*(gamma1gg_1[ii1] + gamma1qq_1[ii1])*LF2*LQ2)/4. + (A1q*eC1qg_1*gamma1qq_2[ii2]*LQ3)/2. + (A2q*gamma1qg_1[ii1]*LQ3)/2. - (pow(B1q,2)*gamma1qg_1[ii1]*LQ3)/2. + B1q*beta0*gamma1qg_1[ii1]*LQ3 - (pow(beta0,2)*gamma1qg_1[ii1]*LQ3)/3. + (A1q*eC1qq_2*gamma1qg_1[ii1]*LQ3)/2. - B1q*gamma1qq_2[ii2]*gamma1qg_1[ii1]*LQ3 + beta0*gamma1qq_2[ii2]*gamma1qg_1[ii1]*LQ3 - ((gamma1gq_2[ii2]*gamma1qg_2[ii2] + pow(gamma1qq_2[ii2],2))*gamma1qg_1[ii1]*LQ3)/2. + (A1q*(eC1qg_1*gamma1gg_1[ii1] + eC1qq_1*gamma1qg_1[ii1])*LQ3)/2. - (B1q*gamma1qg_1[ii1]*(gamma1gg_1[ii1] + gamma1qq_1[ii1])*LQ3)/2. + (beta0*gamma1qg_1[ii1]*(gamma1gg_1[ii1] + gamma1qq_1[ii1])*LQ3)/2. - (gamma1qq_2[ii2]*gamma1qg_1[ii1]*(gamma1gg_1[ii1] + gamma1qq_1[ii1])*LQ3)/2. + (A1q*gamma2qg_1[ii1]*LQ3)/2. + (A1q*B1q*gamma1qg_1[ii1]*LF*LQ3)/2. - (A1q*beta0*gamma1qg_1[ii1]*LF*LQ3)/3. + A1q*gamma1qq_2[ii2]*gamma1qg_1[ii1]*LF*LQ3 + (A1q*gamma1qg_1[ii1]*(gamma1gg_1[ii1] + gamma1qq_1[ii1])*LF*LQ3)/2. - (A1q*B1q*gamma1qg_1[ii1]*LQ4)/2. + (7.*A1q*beta0*gamma1qg_1[ii1]*LQ4)/12. - (A1q*gamma1qq_2[ii2]*gamma1qg_1[ii1]*LQ4)/2. - (A1q*gamma1qg_1[ii1]*(gamma1gg_1[ii1] + gamma1qq_1[ii1])*LQ4)/4. + (pow(A1q,2)*gamma1qg_1[ii1]*LF*LQ4)/8. - (pow(A1q,2)*gamma1qg_1[ii1]*LQ5)/8. + LF*(eC2qg_1*gamma1gg_1[ii1] + gamma1qg_1[ii1]*(eC2qq_1 + eC2qqb_1 + (eC2qqbp_1 + eC2qqp_1)*(-1. + nf))) - LQ*(eC2qg_1*gamma1gg_1[ii1] + gamma1qg_1[ii1]*(eC2qq_1 + eC2qqb_1 + (eC2qqbp_1 + eC2qqp_1)*(-1. + nf))) - LF*LQ*(gamma1gg_1[ii1]*gamma2qg_1[ii1] + gamma1qg_1[ii1]*(gamma2qq_1[ii1] + gamma2qqb_1[ii1] + (gamma2qqbp_1[ii1] + gamma2qqp_1[ii1])*(-1. + nf))) + (gamma1qg_1[ii1]*LF3*(pow(gamma1gg_1[ii1],2) + gamma1gg_1[ii1]*gamma1qq_1[ii1] + pow(gamma1qq_1[ii1],2) + 2.*gamma1gq_1[ii1]*gamma1qg_1[ii1]*nf))/6. - (gamma1qg_1[ii1]*LF2*LQ*(pow(gamma1gg_1[ii1],2) + gamma1gg_1[ii1]*gamma1qq_1[ii1] + pow(gamma1qq_1[ii1],2) + 2.*gamma1gq_1[ii1]*gamma1qg_1[ii1]*nf))/2. + (gamma1qg_1[ii1]*LF*LQ2*(pow(gamma1gg_1[ii1],2) + gamma1gg_1[ii1]*gamma1qq_1[ii1] + pow(gamma1qq_1[ii1],2) + 2.*gamma1gq_1[ii1]*gamma1qg_1[ii1]*nf))/2. - (gamma1qg_1[ii1]*LQ3*(pow(gamma1gg_1[ii1],2) + gamma1gg_1[ii1]*gamma1qq_1[ii1] + pow(gamma1qq_1[ii1],2) + 2.*gamma1gq_1[ii1]*gamma1qg_1[ii1]*nf))/6. + (LF2*(eC1qq_1*gamma1qg_1[ii1]*(gamma1gg_1[ii1] + gamma1qq_1[ii1]) + eC1qg_1*(pow(gamma1gg_1[ii1],2) + 2.*gamma1gq_1[ii1]*gamma1qg_1[ii1]*nf)))/2. - LF*LQ*(eC1qq_1*gamma1qg_1[ii1]*(gamma1gg_1[ii1] + gamma1qq_1[ii1]) + eC1qg_1*(pow(gamma1gg_1[ii1],2) + 2.*gamma1gq_1[ii1]*gamma1qg_1[ii1]*nf)) + (LQ2*(eC1qq_1*gamma1qg_1[ii1]*(gamma1gg_1[ii1] + gamma1qq_1[ii1]) + eC1qg_1*(pow(gamma1gg_1[ii1],2) + 2.*gamma1gq_1[ii1]*gamma1qg_1[ii1]*nf)))/2.;
              //complex <double> H3st_qg_2   = eC2qg_2*gamma1qq_1[ii1]*LF + eC1qg_2*(eC1qg_1*gamma1gq_1[ii1] + eC1qq_1*gamma1qq_1[ii1])*LF + eC2qq_1*gamma1qg_2[ii2]*LF + eC1qq_1*(eC1qg_2*gamma1gg_2[ii2] + eC1qq_2*gamma1qg_2[ii2])*LF + eC1qg_2*gamma2qq_1[ii1]*LF + eC1qq_1*gamma2qg_2[ii2]*LF + (eC1qg_2*gamma2gg_2[ii2] + eC1qq_2*gamma2qg_2[ii2])*LF + gamma3qg_2[ii2]*LF + (beta0*eC1qg_2*gamma1qq_1[ii1]*LF2)/2. + (eC1qg_2*(gamma1gq_1[ii1]*gamma1qg_1[ii1] + pow(gamma1qq_1[ii1],2))*LF2)/2. + (beta1*gamma1qg_2[ii2]*LF2)/2. + (beta0*eC1qq_1*gamma1qg_2[ii2]*LF2)/2. + (eC1qg_1*gamma1gq_1[ii1] + eC1qq_1*gamma1qq_1[ii1])*gamma1qg_2[ii2]*LF2 + (beta0*(eC1qg_2*gamma1gg_2[ii2] + eC1qq_2*gamma1qg_2[ii2])*LF2)/2. + gamma1qq_1[ii1]*(eC1qg_2*gamma1gg_2[ii2] + eC1qq_2*gamma1qg_2[ii2])*LF2 + (eC1qq_1*gamma1qg_2[ii2]*(gamma1gg_2[ii2] + gamma1qq_2[ii2])*LF2)/2. + gamma1qg_2[ii2]*gamma2qq_1[ii1]*LF2 + beta0*gamma2qg_2[ii2]*LF2 + gamma1qq_1[ii1]*gamma2qg_2[ii2]*LF2 + (gamma1qg_2[ii2]*gamma2gg_2[ii2] + gamma1qq_2[ii2]*gamma2qg_2[ii2])* LF2 + (pow(beta0,2)*gamma1qg_2[ii2]*LF3)/3. + beta0*gamma1qq_1[ii1]*gamma1qg_2[ii2]*LF3 + ((gamma1gq_1[ii1]*gamma1qg_1[ii1] + pow(gamma1qq_1[ii1],2))*gamma1qg_2[ii2]*LF3)/2. + (beta0*gamma1qg_2[ii2]*(gamma1gg_2[ii2] + gamma1qq_2[ii2])*LF3)/2. + (gamma1qq_1[ii1]*gamma1qg_2[ii2]*(gamma1gg_2[ii2] + gamma1qq_2[ii2])*LF3)/2. - eC2qg_2*gamma1qq_1[ii1]*LQ - eC1qg_2*(eC1qg_1*gamma1gq_1[ii1] + eC1qq_1*gamma1qq_1[ii1])*LQ - eC2qq_1*gamma1qg_2[ii2]*LQ - eC1qq_1*(eC1qg_2*gamma1gg_2[ii2] + eC1qq_2*gamma1qg_2[ii2])*LQ - eC1qg_2*gamma2qq_1[ii1]*LQ - eC1qq_1*gamma2qg_2[ii2]*LQ - (eC1qg_2*gamma2gg_2[ii2] + eC1qq_2*gamma2qg_2[ii2])*LQ - gamma3qg_2[ii2]*LQ - B1q*eC1qg_2*gamma1qq_1[ii1]*LF*LQ + beta0*eC1qg_2*gamma1qq_1[ii1]*LF*LQ - eC1qg_2*(gamma1gq_1[ii1]*gamma1qg_1[ii1] + pow(gamma1qq_1[ii1],2))*LF*LQ - B2q*gamma1qg_2[ii2]*LF*LQ - B1q*eC1qq_1*gamma1qg_2[ii2]*LF*LQ + beta0*eC1qq_1*gamma1qg_2[ii2]*LF*LQ - 2.*(eC1qg_1*gamma1gq_1[ii1] + eC1qq_1*gamma1qq_1[ii1])*gamma1qg_2[ii2]*LF*LQ - B1q*(eC1qg_2*gamma1gg_2[ii2] + eC1qq_2*gamma1qg_2[ii2])*LF*LQ + beta0*(eC1qg_2*gamma1gg_2[ii2] + eC1qq_2*gamma1qg_2[ii2])*LF*LQ - 2.*gamma1qq_1[ii1]*(eC1qg_2*gamma1gg_2[ii2] + eC1qq_2*gamma1qg_2[ii2])*LF*LQ - eC1qq_1*gamma1qg_2[ii2]*(gamma1gg_2[ii2] + gamma1qq_2[ii2])*LF*LQ - 2.*gamma1qg_2[ii2]*gamma2qq_1[ii1]*LF*LQ - B1q*gamma2qg_2[ii2]*LF*LQ - 2.*gamma1qq_1[ii1]*gamma2qg_2[ii2]*LF*LQ - (gamma1qg_2[ii2]*gamma2gg_2[ii2] + gamma1qq_2[ii2]*gamma2qg_2[ii2])*LF*LQ - (B1q*beta0*gamma1qg_2[ii2]*LF2*LQ)/2. - B1q*gamma1qq_1[ii1]*gamma1qg_2[ii2]*LF2*LQ - beta0*gamma1qq_1[ii1]*gamma1qg_2[ii2]*LF2*LQ - (3.*(gamma1gq_1[ii1]*gamma1qg_1[ii1] + pow(gamma1qq_1[ii1],2))*gamma1qg_2[ii2]*LF2*LQ)/2. - (B1q*gamma1qg_2[ii2]*(gamma1gg_2[ii2] + gamma1qq_2[ii2])*LF2*LQ)/2. - (beta0*gamma1qg_2[ii2]*(gamma1gg_2[ii2] + gamma1qq_2[ii2])*LF2*LQ)/2. - (3.*gamma1qq_1[ii1]*gamma1qg_2[ii2]*(gamma1gg_2[ii2] + gamma1qq_2[ii2])*LF2*LQ)/2. + B1q*eC1qg_2*gamma1qq_1[ii1]*LQ2 - (3.*beta0*eC1qg_2*gamma1qq_1[ii1]*LQ2)/2. + (eC1qg_2*(gamma1gq_1[ii1]*gamma1qg_1[ii1] + pow(gamma1qq_1[ii1],2))*LQ2)/2. + B2q*gamma1qg_2[ii2]*LQ2 - (beta1*gamma1qg_2[ii2]*LQ2)/2. + B1q*eC1qq_1*gamma1qg_2[ii2]*LQ2 - (3.*beta0*eC1qq_1*gamma1qg_2[ii2]*LQ2)/2. + (eC1qg_1*gamma1gq_1[ii1] + eC1qq_1*gamma1qq_1[ii1])* gamma1qg_2[ii2]*LQ2 + B1q*(eC1qg_2*gamma1gg_2[ii2] + eC1qq_2*gamma1qg_2[ii2])*LQ2 - (3.*beta0*(eC1qg_2*gamma1gg_2[ii2] + eC1qq_2*gamma1qg_2[ii2])*LQ2)/2. + gamma1qq_1[ii1]*(eC1qg_2*gamma1gg_2[ii2] + eC1qq_2*gamma1qg_2[ii2])*LQ2 + (eC1qq_1*gamma1qg_2[ii2]*(gamma1gg_2[ii2] + gamma1qq_2[ii2])*LQ2)/2. + gamma1qg_2[ii2]*gamma2qq_1[ii1]*LQ2 + B1q*gamma2qg_2[ii2]*LQ2 - beta0*gamma2qg_2[ii2]*LQ2 + gamma1qq_1[ii1]*gamma2qg_2[ii2]*LQ2 + (gamma1qg_2[ii2]*gamma2gg_2[ii2] + gamma1qq_2[ii2]*gamma2qg_2[ii2])* LQ2 - (A1q*eC1qg_2*gamma1qq_1[ii1]*LF*LQ2)/2. - (A2q*gamma1qg_2[ii2]*LF*LQ2)/2. + (pow(B1q,2)*gamma1qg_2[ii2]*LF*LQ2)/2. - (B1q*beta0*gamma1qg_2[ii2]*LF*LQ2)/2. - (A1q*eC1qq_1*gamma1qg_2[ii2]*LF*LQ2)/2. + 2.*B1q*gamma1qq_1[ii1]*gamma1qg_2[ii2]*LF*LQ2 - beta0*gamma1qq_1[ii1]*gamma1qg_2[ii2]*LF*LQ2 + (3.*(gamma1gq_1[ii1]*gamma1qg_1[ii1] + pow(gamma1qq_1[ii1],2))*gamma1qg_2[ii2]*LF*LQ2)/2. - (A1q*(eC1qg_2*gamma1gg_2[ii2] + eC1qq_2*gamma1qg_2[ii2])*LF*LQ2)/2. + B1q*gamma1qg_2[ii2]*(gamma1gg_2[ii2] + gamma1qq_2[ii2])*LF*LQ2 - (beta0*gamma1qg_2[ii2]*(gamma1gg_2[ii2] + gamma1qq_2[ii2])*LF*LQ2)/2. + (3.*gamma1qq_1[ii1]*gamma1qg_2[ii2]*(gamma1gg_2[ii2] + gamma1qq_2[ii2])*LF*LQ2)/2. - (A1q*gamma2qg_2[ii2]*LF*LQ2)/2. - (A1q*beta0*gamma1qg_2[ii2]*LF2*LQ2)/4. - (A1q*gamma1qq_1[ii1]*gamma1qg_2[ii2]*LF2*LQ2)/2. - (A1q*gamma1qg_2[ii2]*(gamma1gg_2[ii2] + gamma1qq_2[ii2])*LF2*LQ2)/4. + (A1q*eC1qg_2*gamma1qq_1[ii1]*LQ3)/2. + (A2q*gamma1qg_2[ii2]*LQ3)/2. - (pow(B1q,2)*gamma1qg_2[ii2]*LQ3)/2. + B1q*beta0*gamma1qg_2[ii2]*LQ3 - (pow(beta0,2)*gamma1qg_2[ii2]*LQ3)/3. + (A1q*eC1qq_1*gamma1qg_2[ii2]*LQ3)/2. - B1q*gamma1qq_1[ii1]*gamma1qg_2[ii2]*LQ3 + beta0*gamma1qq_1[ii1]*gamma1qg_2[ii2]*LQ3 - ((gamma1gq_1[ii1]*gamma1qg_1[ii1] + pow(gamma1qq_1[ii1],2))*gamma1qg_2[ii2]*LQ3)/2. + (A1q*(eC1qg_2*gamma1gg_2[ii2] + eC1qq_2*gamma1qg_2[ii2])*LQ3)/2. - (B1q*gamma1qg_2[ii2]*(gamma1gg_2[ii2] + gamma1qq_2[ii2])*LQ3)/2. + (beta0*gamma1qg_2[ii2]*(gamma1gg_2[ii2] + gamma1qq_2[ii2])*LQ3)/2. - (gamma1qq_1[ii1]*gamma1qg_2[ii2]*(gamma1gg_2[ii2] + gamma1qq_2[ii2])*LQ3)/2. + (A1q*gamma2qg_2[ii2]*LQ3)/2. + (A1q*B1q*gamma1qg_2[ii2]*LF*LQ3)/2. - (A1q*beta0*gamma1qg_2[ii2]*LF*LQ3)/3. + A1q*gamma1qq_1[ii1]*gamma1qg_2[ii2]*LF*LQ3 + (A1q*gamma1qg_2[ii2]*(gamma1gg_2[ii2] + gamma1qq_2[ii2])*LF*LQ3)/2. - (A1q*B1q*gamma1qg_2[ii2]*LQ4)/2. + (7.*A1q*beta0*gamma1qg_2[ii2]*LQ4)/12. - (A1q*gamma1qq_1[ii1]*gamma1qg_2[ii2]*LQ4)/2. - (A1q*gamma1qg_2[ii2]*(gamma1gg_2[ii2] + gamma1qq_2[ii2])*LQ4)/4. + (pow(A1q,2)*gamma1qg_2[ii2]*LF*LQ4)/8. - (pow(A1q,2)*gamma1qg_2[ii2]*LQ5)/8. + LF*(eC2qg_2*gamma1gg_2[ii2] + gamma1qg_2[ii2]*(eC2qq_2 + eC2qqb_2 + (eC2qqbp_2 + eC2qqp_2)*(-1. + nf))) - LQ*(eC2qg_2*gamma1gg_2[ii2] + gamma1qg_2[ii2]*(eC2qq_2 + eC2qqb_2 + (eC2qqbp_2 + eC2qqp_2)*(-1. + nf))) - LF*LQ*(gamma1gg_2[ii2]*gamma2qg_2[ii2] + gamma1qg_2[ii2]*(gamma2qq_2[ii2] + gamma2qqb_2[ii2] + (gamma2qqbp_2[ii2] + gamma2qqp_2[ii2])*(-1. + nf))) + (gamma1qg_2[ii2]*LF3*(pow(gamma1gg_2[ii2],2) + gamma1gg_2[ii2]*gamma1qq_2[ii2] + pow(gamma1qq_2[ii2],2) + 2.*gamma1gq_2[ii2]*gamma1qg_2[ii2]*nf))/6. - (gamma1qg_2[ii2]*LF2*LQ*(pow(gamma1gg_2[ii2],2) + gamma1gg_2[ii2]*gamma1qq_2[ii2] + pow(gamma1qq_2[ii2],2) + 2.*gamma1gq_2[ii2]*gamma1qg_2[ii2]*nf))/2. + (gamma1qg_2[ii2]*LF*LQ2*(pow(gamma1gg_2[ii2],2) + gamma1gg_2[ii2]*gamma1qq_2[ii2] + pow(gamma1qq_2[ii2],2) + 2.*gamma1gq_2[ii2]*gamma1qg_2[ii2]*nf))/2. - (gamma1qg_2[ii2]*LQ3*(pow(gamma1gg_2[ii2],2) + gamma1gg_2[ii2]*gamma1qq_2[ii2] + pow(gamma1qq_2[ii2],2) + 2.*gamma1gq_2[ii2]*gamma1qg_2[ii2]*nf))/6. + (LF2*(eC1qq_2*gamma1qg_2[ii2]*(gamma1gg_2[ii2] + gamma1qq_2[ii2]) + eC1qg_2*(pow(gamma1gg_2[ii2],2) + 2.*gamma1gq_2[ii2]*gamma1qg_2[ii2]*nf)))/2. - LF*LQ*(eC1qq_2*gamma1qg_2[ii2]*(gamma1gg_2[ii2] + gamma1qq_2[ii2]) + eC1qg_2*(pow(gamma1gg_2[ii2],2) + 2.*gamma1gq_2[ii2]*gamma1qg_2[ii2]*nf)) + (LQ2*(eC1qq_2*gamma1qg_2[ii2]*(gamma1gg_2[ii2] + gamma1qq_2[ii2]) + eC1qg_2*(pow(gamma1gg_2[ii2],2) + 2.*gamma1gq_2[ii2]*gamma1qg_2[ii2]*nf)))/2.;
              //complex <double> H3st_qq_1   = eC2qqb_1*gamma1qq_2[ii2]*LF + eC1qq_2*eC1qg_1*gamma1gq_1[ii1]*LF + (eC2qg_1*gamma1gq_1[ii1] + eC2qqb_1*gamma1qq_1[ii1])*LF + eC1qq_2*gamma2qqb_1[ii1]*LF + (eC1qg_1*gamma2gq_1[ii1] + eC1qq_1*gamma2qqb_1[ii1])*LF + gamma3qqb_1[ii1]*LF + (beta0*eC1qg_1*gamma1gq_1[ii1]*LF2)/2. + eC1qg_1*gamma1qq_2[ii2]*gamma1gq_1[ii1]*LF2 + (eC1qq_2*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF2)/2. + (gamma1gq_1[ii1]*(eC1qq_1*gamma1qg_1[ii1] + eC1qg_1*(gamma1gg_1[ii1] + gamma1qq_1[ii1]))*LF2)/2. + beta0*gamma2qqb_1[ii1]*LF2 + gamma1qq_2[ii2]*gamma2qqb_1[ii1]*LF2 + (gamma1qg_1[ii1]*gamma2gq_1[ii1] + gamma1qq_1[ii1]*gamma2qqb_1[ii1])*LF2 + (beta0*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF3)/2. + (gamma1qq_2[ii2]*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF3)/2. + (gamma1gq_1[ii1]*gamma1qg_1[ii1]*(gamma1gg_1[ii1] + 2.*gamma1qq_1[ii1])*LF3)/6. - eC2qqb_1*gamma1qq_2[ii2]*LQ - eC1qq_2*eC1qg_1*gamma1gq_1[ii1]*LQ - (eC2qg_1*gamma1gq_1[ii1] + eC2qqb_1*gamma1qq_1[ii1])*LQ - eC1qq_2*gamma2qqb_1[ii1]*LQ - (eC1qg_1*gamma2gq_1[ii1] + eC1qq_1*gamma2qqb_1[ii1])*LQ - gamma3qqb_1[ii1]*LQ - B1q*eC1qg_1*gamma1gq_1[ii1]*LF*LQ + beta0*eC1qg_1*gamma1gq_1[ii1]*LF*LQ - 2.*eC1qg_1*gamma1qq_2[ii2]*gamma1gq_1[ii1]*LF*LQ - eC1qq_2*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF*LQ - gamma1gq_1[ii1]*(eC1qq_1*gamma1qg_1[ii1] + eC1qg_1*(gamma1gg_1[ii1] + gamma1qq_1[ii1]))*LF*LQ - B1q*gamma2qqb_1[ii1]*LF*LQ - 2.*gamma1qq_2[ii2]*gamma2qqb_1[ii1]*LF*LQ - (gamma1qg_1[ii1]*gamma2gq_1[ii1] + gamma1qq_1[ii1]*gamma2qqb_1[ii1])*LF*LQ - (gamma1gq_1[ii1]*gamma2qg_1[ii1] + gamma1qq_1[ii1]*gamma2qqb_1[ii1])*LF*LQ - (B1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF2*LQ)/2. - (beta0*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF2*LQ)/ 2. - (3.*gamma1qq_2[ii2]*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF2*LQ)/2. - (gamma1gq_1[ii1]*gamma1qg_1[ii1]*(gamma1gg_1[ii1] + 2.*gamma1qq_1[ii1])*LF2*LQ)/2. + B1q*eC1qg_1*gamma1gq_1[ii1]*LQ2 - (3.*beta0*eC1qg_1*gamma1gq_1[ii1]*LQ2)/2. + eC1qg_1*gamma1qq_2[ii2]*gamma1gq_1[ii1]*LQ2 + (eC1qq_2*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQ2)/2. + (gamma1gq_1[ii1]*(eC1qq_1*gamma1qg_1[ii1] + eC1qg_1*(gamma1gg_1[ii1] + gamma1qq_1[ii1]))*LQ2)/2. + B1q*gamma2qqb_1[ii1]*LQ2 - beta0*gamma2qqb_1[ii1]*LQ2 + gamma1qq_2[ii2]*gamma2qqb_1[ii1]*LQ2 + (gamma1qg_1[ii1]*gamma2gq_1[ii1] + gamma1qq_1[ii1]*gamma2qqb_1[ii1])*LQ2 - (A1q*eC1qg_1*gamma1gq_1[ii1]*LF*LQ2)/2. + B1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF*LQ2 - (beta0*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF*LQ2)/2. + (3.*gamma1qq_2[ii2]*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF*LQ2)/2. + (gamma1gq_1[ii1]*gamma1qg_1[ii1]*(gamma1gg_1[ii1] + 2.*gamma1qq_1[ii1])*LF*LQ2)/2. - (A1q*gamma2qqb_1[ii1]*LF*LQ2)/2. - (A1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF2*LQ2)/4. + (A1q*eC1qg_1*gamma1gq_1[ii1]*LQ3)/2. - (B1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQ3)/2. + (beta0*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQ3)/2. - (gamma1qq_2[ii2]*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQ3)/2. - (gamma1gq_1[ii1]*gamma1qg_1[ii1]*(gamma1gg_1[ii1] + 2.*gamma1qq_1[ii1])*LQ3)/6. + (A1q*gamma2qqb_1[ii1]*LQ3)/2. + (A1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF*LQ3)/2. - (A1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQ4)/4.;
              //complex <double> H3st_qq_2   = eC2qqb_2*gamma1qq_1[ii1]*LF + eC1qq_1*eC1qg_2*gamma1gq_2[ii2]*LF + (eC2qg_2*gamma1gq_2[ii2] + eC2qqb_2*gamma1qq_2[ii2])*LF + eC1qq_1*gamma2qqb_2[ii2]*LF + (eC1qg_2*gamma2gq_2[ii2] + eC1qq_2*gamma2qqb_2[ii2])*LF + gamma3qqb_2[ii2]*LF + (beta0*eC1qg_2*gamma1gq_2[ii2]*LF2)/2. + eC1qg_2*gamma1qq_1[ii1]*gamma1gq_2[ii2]*LF2 + (eC1qq_1*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF2)/2. + (gamma1gq_2[ii2]*(eC1qq_2*gamma1qg_2[ii2] + eC1qg_2*(gamma1gg_2[ii2] + gamma1qq_2[ii2]))*LF2)/2. + beta0*gamma2qqb_2[ii2]*LF2 + gamma1qq_1[ii1]*gamma2qqb_2[ii2]*LF2 + (gamma1qg_2[ii2]*gamma2gq_2[ii2] + gamma1qq_2[ii2]*gamma2qqb_2[ii2])*LF2 + (beta0*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF3)/2. + (gamma1qq_1[ii1]*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF3)/2. + (gamma1gq_2[ii2]*gamma1qg_2[ii2]*(gamma1gg_2[ii2] + 2.*gamma1qq_2[ii2])*LF3)/6. - eC2qqb_2*gamma1qq_1[ii1]*LQ - eC1qq_1*eC1qg_2*gamma1gq_2[ii2]*LQ - (eC2qg_2*gamma1gq_2[ii2] + eC2qqb_2*gamma1qq_2[ii2])*LQ - eC1qq_1*gamma2qqb_2[ii2]*LQ - (eC1qg_2*gamma2gq_2[ii2] + eC1qq_2*gamma2qqb_2[ii2])*LQ - gamma3qqb_2[ii2]*LQ - B1q*eC1qg_2*gamma1gq_2[ii2]*LF*LQ + beta0*eC1qg_2*gamma1gq_2[ii2]*LF*LQ - 2.*eC1qg_2*gamma1qq_1[ii1]*gamma1gq_2[ii2]*LF*LQ - eC1qq_1*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF*LQ - gamma1gq_2[ii2]*(eC1qq_2*gamma1qg_2[ii2] + eC1qg_2*(gamma1gg_2[ii2] + gamma1qq_2[ii2]))*LF*LQ - B1q*gamma2qqb_2[ii2]*LF*LQ - 2.*gamma1qq_1[ii1]*gamma2qqb_2[ii2]*LF*LQ - (gamma1qg_2[ii2]*gamma2gq_2[ii2] + gamma1qq_2[ii2]*gamma2qqb_2[ii2])*LF*LQ - (gamma1gq_2[ii2]*gamma2qg_2[ii2] + gamma1qq_2[ii2]*gamma2qqb_2[ii2])*LF*LQ - (B1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF2*LQ)/2. - (beta0*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF2*LQ)/ 2. - (3.*gamma1qq_1[ii1]*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF2*LQ)/2. - (gamma1gq_2[ii2]*gamma1qg_2[ii2]*(gamma1gg_2[ii2] + 2.*gamma1qq_2[ii2])*LF2*LQ)/2. + B1q*eC1qg_2*gamma1gq_2[ii2]*LQ2 - (3.*beta0*eC1qg_2*gamma1gq_2[ii2]*LQ2)/2. + eC1qg_2*gamma1qq_1[ii1]*gamma1gq_2[ii2]*LQ2 + (eC1qq_1*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQ2)/2. + (gamma1gq_2[ii2]*(eC1qq_2*gamma1qg_2[ii2] + eC1qg_2*(gamma1gg_2[ii2] + gamma1qq_2[ii2]))*LQ2)/2. + B1q*gamma2qqb_2[ii2]*LQ2 - beta0*gamma2qqb_2[ii2]*LQ2 + gamma1qq_1[ii1]*gamma2qqb_2[ii2]*LQ2 + (gamma1qg_2[ii2]*gamma2gq_2[ii2] + gamma1qq_2[ii2]*gamma2qqb_2[ii2])*LQ2 - (A1q*eC1qg_2*gamma1gq_2[ii2]*LF*LQ2)/2. + B1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF*LQ2 - (beta0*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF*LQ2)/2. + (3.*gamma1qq_1[ii1]*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF*LQ2)/2. + (gamma1gq_2[ii2]*gamma1qg_2[ii2]*(gamma1gg_2[ii2] + 2.*gamma1qq_2[ii2])*LF*LQ2)/2. - (A1q*gamma2qqb_2[ii2]*LF*LQ2)/2. - (A1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF2*LQ2)/4. + (A1q*eC1qg_2*gamma1gq_2[ii2]*LQ3)/2. - (B1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQ3)/2. + (beta0*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQ3)/2. - (gamma1qq_1[ii1]*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQ3)/2. - (gamma1gq_2[ii2]*gamma1qg_2[ii2]*(gamma1gg_2[ii2] + 2.*gamma1qq_2[ii2])*LQ3)/6. + (A1q*gamma2qqb_2[ii2]*LQ3)/2. + (A1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF*LQ3)/2. - (A1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQ4)/4.;
              //complex <double> H3st_qqp_1  = (2.*gamma1gq_1[ii1]*gamma1qg_1[ii1]*(gamma1gg_1[ii1] + 2.*gamma1qq_1[ii1])*LQF3 + 6.*gamma1gq_1[ii1]*(eC1qq_1*gamma1qg_1[ii1] + eC1qg_1*(gamma1gg_1[ii1] + gamma1qq_1[ii1]))* (LF2 + LQ2) - 3.*(-4.*eC2qg_1*gamma1gq_1[ii1]*LF - 4.*eC1qg_1*gamma2gq_1[ii1]*LF - 4.*eC1qq_1*gamma2qqbp_1[ii1]*LF - 4.*gamma3qqbp_1[ii1]*LF - 2.*beta0*eC1qg_1*gamma1gq_1[ii1]*LF2 - 4.*eC1qg_1*gamma1qq_2[ii2]*gamma1gq_1[ii1]*LF2 - 4.*gamma1qg_1[ii1]*gamma2gq_1[ii1]*LF2 - 4.*beta0*gamma2qqbp_1[ii1]*LF2 - 4.*gamma1qq_2[ii2]*gamma2qqbp_1[ii1]*LF2 - 4.*gamma1qq_1[ii1]*gamma2qqbp_1[ii1]*LF2 - 2.*beta0*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF3 - 2.*gamma1qq_2[ii2]*gamma1gq_1[ii1]*gamma1qg_1[ii1]* LF3 - 4.*eC2qqbp_1*(gamma1qq_2[ii2] + gamma1qq_1[ii1])*LQF - 2.*eC1qq_2*(2.*eC1qg_1*gamma1gq_1[ii1] + 2.*gamma2qqbp_1[ii1] + gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQF)*LQF + 4.*eC2qg_1*gamma1gq_1[ii1]*LQ + 4.*eC1qg_1*gamma2gq_1[ii1]*LQ + 4.*eC1qq_1*gamma2qqbp_1[ii1]*LQ + 4.*gamma3qqbp_1[ii1]*LQ + 4.*B1q*eC1qg_1*gamma1gq_1[ii1]*LF*LQ - 4.*beta0*eC1qg_1*gamma1gq_1[ii1]*LF*LQ + 8.*eC1qg_1*gamma1qq_2[ii2]*gamma1gq_1[ii1]*LF*LQ + 4.*eC1qg_1*gamma1gg_1[ii1]*gamma1gq_1[ii1]*LF* LQ + 4.*eC1qq_1*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF*LQ + 4.*eC1qg_1*gamma1gq_1[ii1]* gamma1qq_1[ii1]*LF*LQ + 4.*gamma1qg_1[ii1]*gamma2gq_1[ii1]*LF*LQ + 4.*gamma1gq_1[ii1]*gamma2qg_1[ii1]*LF*LQ + 4.*B1q*gamma2qqbp_1[ii1]*LF*LQ + 8.*gamma1qq_2[ii2]*gamma2qqbp_1[ii1]*LF*LQ + 8.*gamma1qq_1[ii1]*gamma2qqbp_1[ii1]*LF*LQ + 2.*B1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF2*LQ + 2.*beta0*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF2* LQ + 6.*gamma1qq_2[ii2]*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF2*LQ - 4.*B1q*eC1qg_1*gamma1gq_1[ii1]*LQ2 + 6.*beta0*eC1qg_1*gamma1gq_1[ii1]*LQ2 - 4.*eC1qg_1*gamma1qq_2[ii2]*gamma1gq_1[ii1]*LQ2 - 4.*gamma1qg_1[ii1]*gamma2gq_1[ii1]*LQ2 - 4.*B1q*gamma2qqbp_1[ii1]*LQ2 + 4.*beta0*gamma2qqbp_1[ii1]*LQ2 - 4.*gamma1qq_2[ii2]*gamma2qqbp_1[ii1]*LQ2 - 4.*gamma1qq_1[ii1]*gamma2qqbp_1[ii1]*LQ2 + 2.*A1q*eC1qg_1*gamma1gq_1[ii1]*LF*LQ2 - 4.*B1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF*LQ2 + 2.*beta0*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF*LQ2 - 6.*gamma1qq_2[ii2]*gamma1gq_1[ii1]* gamma1qg_1[ii1]*LF*LQ2 + 2.*A1q*gamma2qqbp_1[ii1]*LF*LQ2 + A1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF2*LQ2 - 2.*A1q*eC1qg_1*gamma1gq_1[ii1]*LQ3 + 2.*B1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQ3 - 2.*beta0*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQ3 + 2.*gamma1qq_2[ii2]*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQ3 - 2.*A1q*gamma2qqbp_1[ii1]*LQ3 - 2.*A1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF*LQ3 + A1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQ4))/12.;
              //complex <double> H3st_qqp_2  = (2.*gamma1gq_2[ii2]*gamma1qg_2[ii2]*(gamma1gg_2[ii2] + 2.*gamma1qq_2[ii2])*LQF3 + 6.*gamma1gq_2[ii2]*(eC1qq_2*gamma1qg_2[ii2] + eC1qg_2*(gamma1gg_2[ii2] + gamma1qq_2[ii2]))* (LF2 + LQ2) - 3.*(-4.*eC2qg_2*gamma1gq_2[ii2]*LF - 4.*eC1qg_2*gamma2gq_2[ii2]*LF - 4.*eC1qq_2*gamma2qqbp_2[ii2]*LF - 4.*gamma3qqbp_2[ii2]*LF - 2.*beta0*eC1qg_2*gamma1gq_2[ii2]*LF2 - 4.*eC1qg_2*gamma1qq_1[ii1]*gamma1gq_2[ii2]*LF2 - 4.*gamma1qg_2[ii2]*gamma2gq_2[ii2]*LF2 - 4.*beta0*gamma2qqbp_2[ii2]*LF2 - 4.*gamma1qq_1[ii1]*gamma2qqbp_2[ii2]*LF2 - 4.*gamma1qq_2[ii2]*gamma2qqbp_2[ii2]*LF2 - 2.*beta0*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF3 - 2.*gamma1qq_1[ii1]*gamma1gq_2[ii2]*gamma1qg_2[ii2]* LF3 - 4.*eC2qqbp_2*(gamma1qq_1[ii1] + gamma1qq_2[ii2])*LQF - 2.*eC1qq_1*(2.*eC1qg_2*gamma1gq_2[ii2] + 2.*gamma2qqbp_2[ii2] + gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQF)*LQF + 4.*eC2qg_2*gamma1gq_2[ii2]*LQ + 4.*eC1qg_2*gamma2gq_2[ii2]*LQ + 4.*eC1qq_2*gamma2qqbp_2[ii2]*LQ + 4.*gamma3qqbp_2[ii2]*LQ + 4.*B1q*eC1qg_2*gamma1gq_2[ii2]*LF*LQ - 4.*beta0*eC1qg_2*gamma1gq_2[ii2]*LF*LQ + 8.*eC1qg_2*gamma1qq_1[ii1]*gamma1gq_2[ii2]*LF*LQ + 4.*eC1qg_2*gamma1gg_2[ii2]*gamma1gq_2[ii2]*LF* LQ + 4.*eC1qq_2*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF*LQ + 4.*eC1qg_2*gamma1gq_2[ii2]* gamma1qq_2[ii2]*LF*LQ + 4.*gamma1qg_2[ii2]*gamma2gq_2[ii2]*LF*LQ + 4.*gamma1gq_2[ii2]*gamma2qg_2[ii2]*LF*LQ + 4.*B1q*gamma2qqbp_2[ii2]*LF*LQ + 8.*gamma1qq_1[ii1]*gamma2qqbp_2[ii2]*LF*LQ + 8.*gamma1qq_2[ii2]*gamma2qqbp_2[ii2]*LF*LQ + 2.*B1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF2*LQ + 2.*beta0*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF2* LQ + 6.*gamma1qq_1[ii1]*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF2*LQ - 4.*B1q*eC1qg_2*gamma1gq_2[ii2]*LQ2 + 6.*beta0*eC1qg_2*gamma1gq_2[ii2]*LQ2 - 4.*eC1qg_2*gamma1qq_1[ii1]*gamma1gq_2[ii2]*LQ2 - 4.*gamma1qg_2[ii2]*gamma2gq_2[ii2]*LQ2 - 4.*B1q*gamma2qqbp_2[ii2]*LQ2 + 4.*beta0*gamma2qqbp_2[ii2]*LQ2 - 4.*gamma1qq_1[ii1]*gamma2qqbp_2[ii2]*LQ2 - 4.*gamma1qq_2[ii2]*gamma2qqbp_2[ii2]*LQ2 + 2.*A1q*eC1qg_2*gamma1gq_2[ii2]*LF*LQ2 - 4.*B1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF*LQ2 + 2.*beta0*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF*LQ2 - 6.*gamma1qq_1[ii1]*gamma1gq_2[ii2]* gamma1qg_2[ii2]*LF*LQ2 + 2.*A1q*gamma2qqbp_2[ii2]*LF*LQ2 + A1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF2*LQ2 - 2.*A1q*eC1qg_2*gamma1gq_2[ii2]*LQ3 + 2.*B1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQ3 - 2.*beta0*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQ3 + 2.*gamma1qq_1[ii1]*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQ3 - 2.*A1q*gamma2qqbp_2[ii2]*LQ3 - 2.*A1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF*LQ3 + A1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQ4))/12.;
              //complex <double> H3st_qqbp_1 = (2.*gamma1gq_1[ii1]*gamma1qg_1[ii1]*(gamma1gg_1[ii1] + 2.*gamma1qq_1[ii1])*LQF3 + 6.*gamma1gq_1[ii1]*(eC1qq_1*gamma1qg_1[ii1] + eC1qg_1*(gamma1gg_1[ii1] + gamma1qq_1[ii1]))* (LF2 + LQ2) - 3.*(-4.*eC2qg_1*gamma1gq_1[ii1]*LF - 4.*eC1qg_1*gamma2gq_1[ii1]*LF - 4.*eC1qq_1*gamma2qqp_1[ii1]*LF - 4.*gamma3qqp_1[ii1]*LF - 2.*beta0*eC1qg_1*gamma1gq_1[ii1]* LF2 - 4.*eC1qg_1*gamma1qq_2[ii2]*gamma1gq_1[ii1]*LF2 - 4.*gamma1qg_1[ii1]*gamma2gq_1[ii1]* LF2 - 4.*beta0*gamma2qqp_1[ii1]*LF2 - 4.*gamma1qq_2[ii2]*gamma2qqp_1[ii1]*LF2 - 4.*gamma1qq_1[ii1]*gamma2qqp_1[ii1]*LF2 - 2.*beta0*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF3 - 2.*gamma1qq_2[ii2]*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF3 - 4.*eC2qqp_1*(gamma1qq_2[ii2] + gamma1qq_1[ii1])*LQF - 2.*eC1qq_2*(2.*eC1qg_1*gamma1gq_1[ii1] + 2.*gamma2qqp_1[ii1] + gamma1gq_1[ii1]*gamma1qg_1[ii1]* LQF)*LQF + 4.*eC2qg_1*gamma1gq_1[ii1]*LQ + 4.*eC1qg_1*gamma2gq_1[ii1]*LQ + 4.*eC1qq_1*gamma2qqp_1[ii1]*LQ + 4.*gamma3qqp_1[ii1]*LQ + 4.*B1q*eC1qg_1*gamma1gq_1[ii1]*LF*LQ - 4.*beta0*eC1qg_1*gamma1gq_1[ii1]*LF*LQ + 8.*eC1qg_1*gamma1qq_2[ii2]*gamma1gq_1[ii1]*LF*LQ + 4.*eC1qg_1*gamma1gg_1[ii1]*gamma1gq_1[ii1]*LF* LQ + 4.*eC1qq_1*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF*LQ + 4.*eC1qg_1*gamma1gq_1[ii1]* gamma1qq_1[ii1]*LF*LQ + 4.*gamma1qg_1[ii1]*gamma2gq_1[ii1]*LF*LQ + 4.*gamma1gq_1[ii1]*gamma2qg_1[ii1]*LF*LQ + 4.*B1q*gamma2qqp_1[ii1]*LF*LQ + 8.*gamma1qq_2[ii2]*gamma2qqp_1[ii1]*LF*LQ + 8.*gamma1qq_1[ii1]*gamma2qqp_1[ii1]*LF*LQ + 2.*B1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF2*LQ + 2.*beta0*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF2* LQ + 6.*gamma1qq_2[ii2]*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF2*LQ - 4.*B1q*eC1qg_1*gamma1gq_1[ii1]*LQ2 + 6.*beta0*eC1qg_1*gamma1gq_1[ii1]*LQ2 - 4.*eC1qg_1*gamma1qq_2[ii2]*gamma1gq_1[ii1]*LQ2 - 4.*gamma1qg_1[ii1]*gamma2gq_1[ii1]*LQ2 - 4.*B1q*gamma2qqp_1[ii1]*LQ2 + 4.*beta0*gamma2qqp_1[ii1]*LQ2 - 4.*gamma1qq_2[ii2]*gamma2qqp_1[ii1]*LQ2 - 4.*gamma1qq_1[ii1]*gamma2qqp_1[ii1]*LQ2 + 2.*A1q*eC1qg_1*gamma1gq_1[ii1]*LF*LQ2 - 4.*B1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF*LQ2 + 2.*beta0*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF*LQ2 - 6.*gamma1qq_2[ii2]*gamma1gq_1[ii1]* gamma1qg_1[ii1]*LF*LQ2 + 2.*A1q*gamma2qqp_1[ii1]*LF*LQ2 + A1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF2*LQ2 - 2.*A1q*eC1qg_1*gamma1gq_1[ii1]*LQ3 + 2.*B1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQ3 - 2.*beta0*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQ3 + 2.*gamma1qq_2[ii2]*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQ3 - 2.*A1q*gamma2qqp_1[ii1]*LQ3 - 2.*A1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF*LQ3 + A1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQ4))/12.;
              //complex <double> H3st_qqbp_2 = (2.*gamma1gq_2[ii2]*gamma1qg_2[ii2]*(gamma1gg_2[ii2] + 2.*gamma1qq_2[ii2])*LQF3 + 6.*gamma1gq_2[ii2]*(eC1qq_2*gamma1qg_2[ii2] + eC1qg_2*(gamma1gg_2[ii2] + gamma1qq_2[ii2]))* (LF2 + LQ2) - 3.*(-4.*eC2qg_2*gamma1gq_2[ii2]*LF - 4.*eC1qg_2*gamma2gq_2[ii2]*LF - 4.*eC1qq_2*gamma2qqp_2[ii2]*LF - 4.*gamma3qqp_2[ii2]*LF - 2.*beta0*eC1qg_2*gamma1gq_2[ii2]* LF2 - 4.*eC1qg_2*gamma1qq_1[ii1]*gamma1gq_2[ii2]*LF2 - 4.*gamma1qg_2[ii2]*gamma2gq_2[ii2]* LF2 - 4.*beta0*gamma2qqp_2[ii2]*LF2 - 4.*gamma1qq_1[ii1]*gamma2qqp_2[ii2]*LF2 - 4.*gamma1qq_2[ii2]*gamma2qqp_2[ii2]*LF2 - 2.*beta0*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF3 - 2.*gamma1qq_1[ii1]*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF3 - 4.*eC2qqp_2*(gamma1qq_1[ii1] + gamma1qq_2[ii2])*LQF - 2.*eC1qq_1*(2.*eC1qg_2*gamma1gq_2[ii2] + 2.*gamma2qqp_2[ii2] + gamma1gq_2[ii2]*gamma1qg_2[ii2]* LQF)*LQF + 4.*eC2qg_2*gamma1gq_2[ii2]*LQ + 4.*eC1qg_2*gamma2gq_2[ii2]*LQ + 4.*eC1qq_2*gamma2qqp_2[ii2]*LQ + 4.*gamma3qqp_2[ii2]*LQ + 4.*B1q*eC1qg_2*gamma1gq_2[ii2]*LF*LQ - 4.*beta0*eC1qg_2*gamma1gq_2[ii2]*LF*LQ + 8.*eC1qg_2*gamma1qq_1[ii1]*gamma1gq_2[ii2]*LF*LQ + 4.*eC1qg_2*gamma1gg_2[ii2]*gamma1gq_2[ii2]*LF* LQ + 4.*eC1qq_2*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF*LQ + 4.*eC1qg_2*gamma1gq_2[ii2]* gamma1qq_2[ii2]*LF*LQ + 4.*gamma1qg_2[ii2]*gamma2gq_2[ii2]*LF*LQ + 4.*gamma1gq_2[ii2]*gamma2qg_2[ii2]*LF*LQ + 4.*B1q*gamma2qqp_2[ii2]*LF*LQ + 8.*gamma1qq_1[ii1]*gamma2qqp_2[ii2]*LF*LQ + 8.*gamma1qq_2[ii2]*gamma2qqp_2[ii2]*LF*LQ + 2.*B1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF2*LQ + 2.*beta0*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF2* LQ + 6.*gamma1qq_1[ii1]*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF2*LQ - 4.*B1q*eC1qg_2*gamma1gq_2[ii2]*LQ2 + 6.*beta0*eC1qg_2*gamma1gq_2[ii2]*LQ2 - 4.*eC1qg_2*gamma1qq_1[ii1]*gamma1gq_2[ii2]*LQ2 - 4.*gamma1qg_2[ii2]*gamma2gq_2[ii2]*LQ2 - 4.*B1q*gamma2qqp_2[ii2]*LQ2 + 4.*beta0*gamma2qqp_2[ii2]*LQ2 - 4.*gamma1qq_1[ii1]*gamma2qqp_2[ii2]*LQ2 - 4.*gamma1qq_2[ii2]*gamma2qqp_2[ii2]*LQ2 + 2.*A1q*eC1qg_2*gamma1gq_2[ii2]*LF*LQ2 - 4.*B1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF*LQ2 + 2.*beta0*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF*LQ2 - 6.*gamma1qq_1[ii1]*gamma1gq_2[ii2]* gamma1qg_2[ii2]*LF*LQ2 + 2.*A1q*gamma2qqp_2[ii2]*LF*LQ2 + A1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF2*LQ2 - 2.*A1q*eC1qg_2*gamma1gq_2[ii2]*LQ3 + 2.*B1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQ3 - 2.*beta0*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQ3 + 2.*gamma1qq_1[ii1]*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQ3 - 2.*A1q*gamma2qqp_2[ii2]*LQ3 - 2.*A1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF*LQ3 + A1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQ4))/12.;
              //fix
              complex <double> H3st_qqb    = (LQF*(24.*eC1qg_2*eC1qq_1*gamma1gq_2[ii2] + 24.*eC2qg_2*gamma1gq_2[ii2] + 24.*eC1qq_2*eC1qq_1*gamma1qq_2[ii2] + 24.*eC2qq_2*gamma1qq_2[ii2] + 24.*eC2qq_1*gamma1qq_2[ii2] + 24.*eC1qq_2*eC1qg_1*gamma1gq_1[ii1] + 24.*eC2qg_1*gamma1gq_1[ii1] + 24.*eC1qq_2*eC1qq_1*gamma1qq_1[ii1] + 24.*eC2qq_2*gamma1qq_1[ii1] + 24.*eC2qq_1*gamma1qq_1[ii1] + 24.*eC1qg_2*gamma2gq_2[ii2] + 24.*eC1qq_2*gamma2qq_2[ii2] + 24.*eC1qq_1*gamma2qq_2[ii2] + 24.*eC1qg_1*gamma2gq_1[ii1] + 24.*eC1qq_2*gamma2qq_1[ii1] + 24.*eC1qq_1*gamma2qq_1[ii1] + 24.*gamma3qq_2[ii2] + 24.*gamma3qq_1[ii1] + 12.*beta0*eC1qg_2*gamma1gq_2[ii2]*LF + 12.*eC1qg_2*gamma1gg_2[ii2]*gamma1gq_2[ii2]*LF + 12.*eC1qq_2*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF + 12.*eC1qq_1*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF + 12.*beta1*gamma1qq_2[ii2]*LF + 12.*beta0*eC1qq_2*gamma1qq_2[ii2]*LF + 12.*beta0*eC1qq_1*gamma1qq_2[ii2]*LF + 12.*eC1qg_2*gamma1gq_2[ii2]*gamma1qq_2[ii2]*LF + 12.*eC1qq_2*pow(gamma1qq_2[ii2],2)*LF + 12.*eC1qq_1*pow(gamma1qq_2[ii2],2)*LF + 12.*beta0*eC1qg_1*gamma1gq_1[ii1]*LF + 24.*eC1qg_1*gamma1qq_2[ii2]*gamma1gq_1[ii1]*LF + 12.*eC1qg_1*gamma1gg_1[ii1]*gamma1gq_1[ii1]*LF + 12.*eC1qq_2*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF + 12.*eC1qq_1*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF + 12.*beta1*gamma1qq_1[ii1]*LF + 12.*beta0*eC1qq_2*gamma1qq_1[ii1]*LF + 12.*beta0*eC1qq_1*gamma1qq_1[ii1]*LF + 24.*eC1qg_2*gamma1gq_2[ii2]*gamma1qq_1[ii1]*LF + 24.*eC1qq_2*gamma1qq_2[ii2]*gamma1qq_1[ii1]*LF + 24.*eC1qq_1*gamma1qq_2[ii2]*gamma1qq_1[ii1]*LF + 12.*eC1qg_1*gamma1gq_1[ii1]*gamma1qq_1[ii1]*LF + 12.*eC1qq_2*pow(gamma1qq_1[ii1],2)*LF + 12.*eC1qq_1*pow(gamma1qq_1[ii1],2)*LF + 12.*gamma1qg_2[ii2]*gamma2gq_2[ii2]*LF + 12.*gamma1gq_2[ii2]*gamma2qg_2[ii2]*LF + 24.*beta0*gamma2qq_2[ii2]*LF + 24.*gamma1qq_2[ii2]*gamma2qq_2[ii2]*LF + 24.*gamma1qq_1[ii1]*gamma2qq_2[ii2]*LF + 12.*gamma1qg_1[ii1]*gamma2gq_1[ii1]*LF + 12.*gamma1gq_1[ii1]*gamma2qg_1[ii1]*LF + 24.*beta0*gamma2qq_1[ii1]*LF + 24.*gamma1qq_2[ii2]*gamma2qq_1[ii1]*LF + 24.*gamma1qq_1[ii1]*gamma2qq_1[ii1]*LF + 12.*beta0*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF2 + 8.*pow(beta0,2)*gamma1qq_2[ii2]*LF2 + 12.*beta0*pow(gamma1qq_2[ii2],2)*LF2 + 12.*beta0*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF2 + 12.*gamma1qq_2[ii2]*gamma1gq_1[ii1]*gamma1qg_1[ii1]* LF2 + 8.*pow(beta0,2)*gamma1qq_1[ii1]*LF2 + 12.*gamma1gq_2[ii2]*gamma1qg_2[ii2]*gamma1qq_1[ii1]* LF2 + 24.*beta0*gamma1qq_2[ii2]*gamma1qq_1[ii1]*LF2 + 12.*pow(gamma1qq_2[ii2],2)*gamma1qq_1[ii1]* LF2 + 12.*beta0*pow(gamma1qq_1[ii1],2)*LF2 + 12.*gamma1qq_2[ii2]*pow(gamma1qq_1[ii1],2)*LF2 + 4.*(gamma1gg_2[ii2]*gamma1gq_2[ii2]*gamma1qg_2[ii2] + 2.*gamma1gq_2[ii2]*gamma1qg_2[ii2]*gamma1qq_2[ii2] + pow(gamma1qq_2[ii2],3))*LQF2 + 4.*(gamma1gg_1[ii1]*gamma1gq_1[ii1]*gamma1qg_1[ii1] + 2.*gamma1gq_1[ii1]*gamma1qg_1[ii1]*gamma1qq_1[ii1] + pow(gamma1qq_1[ii1],3))*LQF2 - 24.*B1q*eC1qg_2*gamma1gq_2[ii2]*LQ + 36.*beta0*eC1qg_2*gamma1gq_2[ii2]*LQ - 12.*eC1qg_2*gamma1gg_2[ii2]*gamma1gq_2[ii2]*LQ - 12.*eC1qq_2*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQ - 12.*eC1qq_1*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQ - 24.*B2q*gamma1qq_2[ii2]*LQ + 12.*beta1*gamma1qq_2[ii2]*LQ - 24.*B1q*eC1qq_2*gamma1qq_2[ii2]*LQ + 36.*beta0*eC1qq_2*gamma1qq_2[ii2]*LQ - 24.*B1q*eC1qq_1*gamma1qq_2[ii2]*LQ + 36.*beta0*eC1qq_1*gamma1qq_2[ii2]*LQ - 12.*eC1qg_2*gamma1gq_2[ii2]*gamma1qq_2[ii2]*LQ - 12.*eC1qq_2*pow(gamma1qq_2[ii2],2)*LQ - 12.*eC1qq_1*pow(gamma1qq_2[ii2],2)*LQ - 24.*B1q*eC1qg_1*gamma1gq_1[ii1]*LQ + 36.*beta0*eC1qg_1*gamma1gq_1[ii1]*LQ - 24.*eC1qg_1*gamma1qq_2[ii2]*gamma1gq_1[ii1]*LQ - 12.*eC1qg_1*gamma1gg_1[ii1]*gamma1gq_1[ii1]*LQ - 12.*eC1qq_2*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQ - 12.*eC1qq_1*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQ - 24.*B2q*gamma1qq_1[ii1]*LQ + 12.*beta1*gamma1qq_1[ii1]*LQ - 24.*B1q*eC1qq_2*gamma1qq_1[ii1]*LQ + 36.*beta0*eC1qq_2*gamma1qq_1[ii1]*LQ - 24.*B1q*eC1qq_1*gamma1qq_1[ii1]*LQ + 36.*beta0*eC1qq_1*gamma1qq_1[ii1]*LQ - 24.*eC1qg_2*gamma1gq_2[ii2]*gamma1qq_1[ii1]*LQ - 24.*eC1qq_2*gamma1qq_2[ii2]*gamma1qq_1[ii1]*LQ - 24.*eC1qq_1*gamma1qq_2[ii2]*gamma1qq_1[ii1]*LQ - 12.*eC1qg_1*gamma1gq_1[ii1]*gamma1qq_1[ii1]*LQ - 12.*eC1qq_2*pow(gamma1qq_1[ii1],2)*LQ - 12.*eC1qq_1*pow(gamma1qq_1[ii1],2)*LQ - 12.*gamma1qg_2[ii2]*gamma2gq_2[ii2]*LQ - 12.*gamma1gq_2[ii2]*gamma2qg_2[ii2]*LQ - 24.*B1q*gamma2qq_2[ii2]*LQ + 24.*beta0*gamma2qq_2[ii2]*LQ - 24.*gamma1qq_2[ii2]*gamma2qq_2[ii2]*LQ - 24.*gamma1qq_1[ii1]*gamma2qq_2[ii2]*LQ - 12.*gamma1qg_1[ii1]*gamma2gq_1[ii1]*LQ - 12.*gamma1gq_1[ii1]*gamma2qg_1[ii1]*LQ - 24.*B1q*gamma2qq_1[ii1]*LQ + 24.*beta0*gamma2qq_1[ii1]*LQ - 24.*gamma1qq_2[ii2]*gamma2qq_1[ii1]*LQ - 24.*gamma1qq_1[ii1]*gamma2qq_1[ii1]*LQ - 12.*B1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF*LQ - 12.*B1q*beta0*gamma1qq_2[ii2]*LF*LQ + 8.*pow(beta0,2)*gamma1qq_2[ii2]*LF*LQ - 12.*B1q*pow(gamma1qq_2[ii2],2)*LF*LQ - 12.*B1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF*LQ - 24.*gamma1qq_2[ii2]*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF* LQ - 12.*B1q*beta0*gamma1qq_1[ii1]*LF*LQ + 8.*pow(beta0,2)*gamma1qq_1[ii1]*LF*LQ - 24.*gamma1gq_2[ii2]*gamma1qg_2[ii2]*gamma1qq_1[ii1]*LF*LQ - 24.*B1q*gamma1qq_2[ii2]*gamma1qq_1[ii1]*LF* LQ - 24.*pow(gamma1qq_2[ii2],2)*gamma1qq_1[ii1]*LF*LQ - 12.*B1q*pow(gamma1qq_1[ii1],2)*LF*LQ - 24.*gamma1qq_2[ii2]*pow(gamma1qq_1[ii1],2)*LF*LQ - 12.*A1q*eC1qg_2*gamma1gq_2[ii2]*LQ2 + 12.*B1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQ2 - 12.*beta0*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQ2 - 12.*A2q*gamma1qq_2[ii2]*LQ2 + 12.*pow(B1q,2)*gamma1qq_2[ii2]*LQ2 - 24.*B1q*beta0*gamma1qq_2[ii2]*LQ2 + 8.*pow(beta0,2)*gamma1qq_2[ii2]*LQ2 - 12.*A1q*eC1qq_2*gamma1qq_2[ii2]*LQ2 - 12.*A1q*eC1qq_1*gamma1qq_2[ii2]*LQ2 + 12.*B1q*pow(gamma1qq_2[ii2],2)*LQ2 - 12.*beta0*pow(gamma1qq_2[ii2],2)*LQ2 - 12.*A1q*eC1qg_1*gamma1gq_1[ii1]*LQ2 + 12.*B1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQ2 - 12.*beta0*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQ2 + 12.*gamma1qq_2[ii2]*gamma1gq_1[ii1]*gamma1qg_1[ii1]* LQ2 - 12.*A2q*gamma1qq_1[ii1]*LQ2 + 12.*pow(B1q,2)*gamma1qq_1[ii1]*LQ2 - 24.*B1q*beta0*gamma1qq_1[ii1]*LQ2 + 8.*pow(beta0,2)*gamma1qq_1[ii1]*LQ2 - 12.*A1q*eC1qq_2*gamma1qq_1[ii1]*LQ2 - 12.*A1q*eC1qq_1*gamma1qq_1[ii1]*LQ2 + 12.*gamma1gq_2[ii2]*gamma1qg_2[ii2]*gamma1qq_1[ii1]*LQ2 + 24.*B1q*gamma1qq_2[ii2]*gamma1qq_1[ii1]* LQ2 - 24.*beta0*gamma1qq_2[ii2]*gamma1qq_1[ii1]*LQ2 + 12.*pow(gamma1qq_2[ii2],2)*gamma1qq_1[ii1]* LQ2 + 12.*B1q*pow(gamma1qq_1[ii1],2)*LQ2 - 12.*beta0*pow(gamma1qq_1[ii1],2)*LQ2 + 12.*gamma1qq_2[ii2]*pow(gamma1qq_1[ii1],2)*LQ2 - 12.*A1q*gamma2qq_2[ii2]*LQ2 - 12.*A1q*gamma2qq_1[ii1]*LQ2 - 6.*A1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF*LQ2 - 6.*A1q*beta0*gamma1qq_2[ii2]*LF*LQ2 - 6.*A1q*pow(gamma1qq_2[ii2],2)*LF*LQ2 - 6.*A1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF*LQ2 - 6.*A1q*beta0*gamma1qq_1[ii1]*LF*LQ2 - 12.*A1q*gamma1qq_2[ii2]*gamma1qq_1[ii1]*LF*LQ2 - 6.*A1q*pow(gamma1qq_1[ii1],2)*LF*LQ2 + 6.*A1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQ3 + 12.*A1q*B1q*gamma1qq_2[ii2]*LQ3 - 14.*A1q*beta0*gamma1qq_2[ii2]*LQ3 + 6.*A1q*pow(gamma1qq_2[ii2],2)*LQ3 + 6.*A1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQ3 + 12.*A1q*B1q*gamma1qq_1[ii1]*LQ3 - 14.*A1q*beta0*gamma1qq_1[ii1]*LQ3 + 12.*A1q*gamma1qq_2[ii2]*gamma1qq_1[ii1]*LQ3 + 6.*A1q*pow(gamma1qq_1[ii1],2)*LQ3 + 3.*pow(A1q,2)*gamma1qq_2[ii2]*LQ4 + 3.*pow(A1q,2)*gamma1qq_1[ii1]*LQ4))/24.;
              complex <double> H3st_qq_1   = (LQF*(12.*eC2qg_1*gamma1gq_1[ii1] + 12.*eC2qqb_1*(gamma1qq_2[ii2] + gamma1qq_1[ii1]) + 12.*eC1qg_1*gamma2gq_1[ii1] + 12.*eC1qq_1*gamma2qqb_1[ii1] + 12.*gamma3qqb_1[ii1] + 6.*beta0*eC1qg_1*gamma1gq_1[ii1]*LF + 12.*eC1qg_1*gamma1qq_2[ii2]*gamma1gq_1[ii1]*LF + 6.*eC1qg_1*gamma1gg_1[ii1]*gamma1gq_1[ii1]*LF + 6.*eC1qq_1*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF + 6.*eC1qg_1*gamma1gq_1[ii1]*gamma1qq_1[ii1]*LF + 6.*gamma1qg_1[ii1]*gamma2gq_1[ii1]*LF + 6.*gamma1gq_1[ii1]*gamma2qg_1[ii1]*LF + 12.*beta0*gamma2qqb_1[ii1]*LF + 12.*gamma1qq_2[ii2]*gamma2qqb_1[ii1]*LF + 12.*gamma1qq_1[ii1]*gamma2qqb_1[ii1]*LF + 6.*beta0*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF2 + 6.*gamma1qq_2[ii2]*gamma1gq_1[ii1]*gamma1qg_1[ii1]* LF2 + 2.*gamma1gg_1[ii1]*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF2 + 4.*gamma1gq_1[ii1]*gamma1qg_1[ii1]*gamma1qq_1[ii1]*LF2 + 6.*eC1qq_2*(2.*eC1qg_1*gamma1gq_1[ii1] + 2.*gamma2qqb_1[ii1] + gamma1gq_1[ii1]*gamma1qg_1[ii1]* LQF) - 12.*B1q*eC1qg_1*gamma1gq_1[ii1]*LQ + 18.*beta0*eC1qg_1*gamma1gq_1[ii1]* LQ - 12.*eC1qg_1*gamma1qq_2[ii2]*gamma1gq_1[ii1]*LQ - 6.*eC1qg_1*gamma1gg_1[ii1]*gamma1gq_1[ii1]* LQ - 6.*eC1qq_1*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQ - 6.*eC1qg_1*gamma1gq_1[ii1]*gamma1qq_1[ii1]* LQ - 6.*gamma1qg_1[ii1]*gamma2gq_1[ii1]*LQ - 6.*gamma1gq_1[ii1]*gamma2qg_1[ii1]*LQ - 12.*B1q*gamma2qqb_1[ii1]*LQ + 12.*beta0*gamma2qqb_1[ii1]*LQ - 12.*gamma1qq_2[ii2]*gamma2qqb_1[ii1]*LQ - 12.*gamma1qq_1[ii1]*gamma2qqb_1[ii1]*LQ - 6.*B1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF*LQ - 12.*gamma1qq_2[ii2]*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF* LQ - 4.*gamma1gg_1[ii1]*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF*LQ - 8.*gamma1gq_1[ii1]*gamma1qg_1[ii1]*gamma1qq_1[ii1]*LF*LQ - 6.*A1q*eC1qg_1*gamma1gq_1[ii1]*LQ2 + 6.*B1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQ2 - 6.*beta0*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQ2 + 6.*gamma1qq_2[ii2]*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQ2 + 2.*gamma1gg_1[ii1]*gamma1gq_1[ii1]* gamma1qg_1[ii1]*LQ2 + 4.*gamma1gq_1[ii1]*gamma1qg_1[ii1]*gamma1qq_1[ii1]*LQ2 - 6.*A1q*gamma2qqb_1[ii1]*LQ2 - 3.*A1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF*LQ2 + 3.*A1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQ3))/12.;
              complex <double> H3st_qqp_1  = (LQF*(12.*eC2qg_1*gamma1gq_1[ii1] + 12.*eC2qqbp_1*(gamma1qq_2[ii2] + gamma1qq_1[ii1]) + 12.*eC1qg_1*gamma2gq_1[ii1] + 12.*eC1qq_1*gamma2qqbp_1[ii1] + 12.*gamma3qqbp_1[ii1] + 6.*beta0*eC1qg_1*gamma1gq_1[ii1]*LF + 12.*eC1qg_1*gamma1qq_2[ii2]*gamma1gq_1[ii1]*LF + 6.*eC1qg_1*gamma1gg_1[ii1]*gamma1gq_1[ii1]*LF + 6.*eC1qq_1*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF + 6.*eC1qg_1*gamma1gq_1[ii1]*gamma1qq_1[ii1]*LF + 6.*gamma1qg_1[ii1]*gamma2gq_1[ii1]*LF + 6.*gamma1gq_1[ii1]*gamma2qg_1[ii1]*LF + 12.*beta0*gamma2qqbp_1[ii1]*LF + 12.*gamma1qq_2[ii2]*gamma2qqbp_1[ii1]*LF + 12.*gamma1qq_1[ii1]*gamma2qqbp_1[ii1]*LF + 6.*beta0*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF2 + 6.*gamma1qq_2[ii2]*gamma1gq_1[ii1]*gamma1qg_1[ii1]* LF2 + 2.*gamma1gg_1[ii1]*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF2 + 4.*gamma1gq_1[ii1]*gamma1qg_1[ii1]*gamma1qq_1[ii1]*LF2 + 6.*eC1qq_2*(2.*eC1qg_1*gamma1gq_1[ii1] + 2.*gamma2qqbp_1[ii1] + gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQF) - 12.*B1q*eC1qg_1*gamma1gq_1[ii1]*LQ + 18.*beta0*eC1qg_1*gamma1gq_1[ii1]*LQ - 12.*eC1qg_1*gamma1qq_2[ii2]*gamma1gq_1[ii1]*LQ - 6.*eC1qg_1*gamma1gg_1[ii1]*gamma1gq_1[ii1]*LQ - 6.*eC1qq_1*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQ - 6.*eC1qg_1*gamma1gq_1[ii1]*gamma1qq_1[ii1]*LQ - 6.*gamma1qg_1[ii1]*gamma2gq_1[ii1]*LQ - 6.*gamma1gq_1[ii1]*gamma2qg_1[ii1]*LQ - 12.*B1q*gamma2qqbp_1[ii1]*LQ + 12.*beta0*gamma2qqbp_1[ii1]*LQ - 12.*gamma1qq_2[ii2]*gamma2qqbp_1[ii1]*LQ - 12.*gamma1qq_1[ii1]*gamma2qqbp_1[ii1]*LQ - 6.*B1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF*LQ - 12.*gamma1qq_2[ii2]*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF*LQ - 4.*gamma1gg_1[ii1]*gamma1gq_1[ii1]* gamma1qg_1[ii1]*LF*LQ - 8.*gamma1gq_1[ii1]*gamma1qg_1[ii1]*gamma1qq_1[ii1]*LF*LQ - 6.*A1q*eC1qg_1*gamma1gq_1[ii1]*LQ2 + 6.*B1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQ2 - 6.*beta0*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQ2 + 6.*gamma1qq_2[ii2]*gamma1gq_1[ii1]*gamma1qg_1[ii1]* LQ2 + 2.*gamma1gg_1[ii1]*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQ2 + 4.*gamma1gq_1[ii1]*gamma1qg_1[ii1]*gamma1qq_1[ii1]*LQ2 - 6.*A1q*gamma2qqbp_1[ii1]*LQ2 - 3.*A1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF*LQ2 + 3.*A1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQ3))/12.;
              complex <double> H3st_qqbp_1 = (LQF*(12.*eC2qg_1*gamma1gq_1[ii1] + 12.*eC2qqp_1*(gamma1qq_2[ii2] + gamma1qq_1[ii1]) + 12.*eC1qg_1*gamma2gq_1[ii1] + 12.*eC1qq_1*gamma2qqp_1[ii1] + 12.*gamma3qqp_1[ii1] + 6.*beta0*eC1qg_1*gamma1gq_1[ii1]*LF + 12.*eC1qg_1*gamma1qq_2[ii2]*gamma1gq_1[ii1]*LF + 6.*eC1qg_1*gamma1gg_1[ii1]*gamma1gq_1[ii1]*LF + 6.*eC1qq_1*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF + 6.*eC1qg_1*gamma1gq_1[ii1]*gamma1qq_1[ii1]*LF + 6.*gamma1qg_1[ii1]*gamma2gq_1[ii1]*LF + 6.*gamma1gq_1[ii1]*gamma2qg_1[ii1]*LF + 12.*beta0*gamma2qqp_1[ii1]*LF + 12.*gamma1qq_2[ii2]*gamma2qqp_1[ii1]*LF + 12.*gamma1qq_1[ii1]*gamma2qqp_1[ii1]*LF + 6.*beta0*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF2 + 6.*gamma1qq_2[ii2]*gamma1gq_1[ii1]*gamma1qg_1[ii1]* LF2 + 2.*gamma1gg_1[ii1]*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF2 + 4.*gamma1gq_1[ii1]*gamma1qg_1[ii1]*gamma1qq_1[ii1]*LF2 + 6.*eC1qq_2*(2.*eC1qg_1*gamma1gq_1[ii1] + 2.*gamma2qqp_1[ii1] + gamma1gq_1[ii1]*gamma1qg_1[ii1]* LQF) - 12.*B1q*eC1qg_1*gamma1gq_1[ii1]*LQ + 18.*beta0*eC1qg_1*gamma1gq_1[ii1]* LQ - 12.*eC1qg_1*gamma1qq_2[ii2]*gamma1gq_1[ii1]*LQ - 6.*eC1qg_1*gamma1gg_1[ii1]*gamma1gq_1[ii1]* LQ - 6.*eC1qq_1*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQ - 6.*eC1qg_1*gamma1gq_1[ii1]*gamma1qq_1[ii1]* LQ - 6.*gamma1qg_1[ii1]*gamma2gq_1[ii1]*LQ - 6.*gamma1gq_1[ii1]*gamma2qg_1[ii1]*LQ - 12.*B1q*gamma2qqp_1[ii1]*LQ + 12.*beta0*gamma2qqp_1[ii1]*LQ - 12.*gamma1qq_2[ii2]*gamma2qqp_1[ii1]*LQ - 12.*gamma1qq_1[ii1]*gamma2qqp_1[ii1]*LQ - 6.*B1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF*LQ - 12.*gamma1qq_2[ii2]*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF* LQ - 4.*gamma1gg_1[ii1]*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF*LQ - 8.*gamma1gq_1[ii1]*gamma1qg_1[ii1]*gamma1qq_1[ii1]*LF*LQ - 6.*A1q*eC1qg_1*gamma1gq_1[ii1]*LQ2 + 6.*B1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQ2 - 6.*beta0*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQ2 + 6.*gamma1qq_2[ii2]*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQ2 + 2.*gamma1gg_1[ii1]*gamma1gq_1[ii1]* gamma1qg_1[ii1]*LQ2 + 4.*gamma1gq_1[ii1]*gamma1qg_1[ii1]*gamma1qq_1[ii1]*LQ2 - 6.*A1q*gamma2qqp_1[ii1]*LQ2 - 3.*A1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF*LQ2 + 3.*A1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQ3))/12.;
              complex <double> H3st_qg_1   = (LQF*(24.*eC2qg_1*gamma1qq_2[ii2] + 24.*eC2qg_1*gamma1gg_1[ii1] + 24.*eC2qq_2*gamma1qg_1[ii1] + 24.*eC2qq_1*gamma1qg_1[ii1] + 24.*eC2qqb_1*gamma1qg_1[ii1] - 24.*eC2qqbp_1*gamma1qg_1[ii1] - 24.*eC2qqp_1*gamma1qg_1[ii1] + 24.*eC1qg_1*gamma2qq_2[ii2] + 24.*eC1qg_1*gamma2gg_1[ii1] + 24.*eC1qq_1*gamma2qg_1[ii1] + 24.*gamma3qg_1[ii1] + 12.*eC1qg_1*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF + 12.*beta0*eC1qg_1*gamma1qq_2[ii2]*LF + 12.*eC1qg_1*pow(gamma1qq_2[ii2],2)*LF + 12.*beta0*eC1qg_1*gamma1gg_1[ii1]*LF + 24.*eC1qg_1*gamma1qq_2[ii2]*gamma1gg_1[ii1]*LF + 12.*eC1qg_1*pow(gamma1gg_1[ii1],2)*LF + 12.*beta1*gamma1qg_1[ii1]*LF + 12.*beta0*eC1qq_1*gamma1qg_1[ii1]*LF + 24.*eC1qq_1*gamma1qq_2[ii2]*gamma1qg_1[ii1]*LF + 12.*eC1qq_1*gamma1gg_1[ii1]*gamma1qg_1[ii1]*LF + 12.*eC1qq_1*gamma1qg_1[ii1]*gamma1qq_1[ii1]*LF + 24.*gamma1qg_1[ii1]*gamma2qq_2[ii2]*LF + 12.*gamma1qg_1[ii1]*gamma2gg_1[ii1]*LF + 24.*beta0*gamma2qg_1[ii1]*LF + 24.*gamma1qq_2[ii2]*gamma2qg_1[ii1]*LF + 12.*gamma1gg_1[ii1]*gamma2qg_1[ii1]*LF + 12.*gamma1qq_1[ii1]*gamma2qg_1[ii1]*LF + 12.*gamma1qg_1[ii1]*gamma2qq_1[ii1]*LF + 12.*gamma1qg_1[ii1]*gamma2qqb_1[ii1]*LF - 12.*gamma1qg_1[ii1]*gamma2qqbp_1[ii1]*LF - 12.*gamma1qg_1[ii1]*gamma2qqp_1[ii1]*LF + 8.*pow(beta0,2)*gamma1qg_1[ii1]*LF2 + 12.*gamma1gq_2[ii2]*gamma1qg_2[ii2]*gamma1qg_1[ii1]*LF2 + 24.*beta0*gamma1qq_2[ii2]*gamma1qg_1[ii1]* LF2 + 12.*pow(gamma1qq_2[ii2],2)*gamma1qg_1[ii1]*LF2 + 12.*beta0*gamma1gg_1[ii1]*gamma1qg_1[ii1]* LF2 + 12.*gamma1qq_2[ii2]*gamma1gg_1[ii1]*gamma1qg_1[ii1]*LF2 + 4.*pow(gamma1gg_1[ii1],2)*gamma1qg_1[ii1]*LF2 + 12.*beta0*gamma1qg_1[ii1]*gamma1qq_1[ii1]*LF2 + 12.*gamma1qq_2[ii2]*gamma1qg_1[ii1]*gamma1qq_1[ii1]*LF2 + 4.*gamma1gg_1[ii1]*gamma1qg_1[ii1]* gamma1qq_1[ii1]*LF2 + 4.*gamma1qg_1[ii1]*pow(gamma1qq_1[ii1],2)*LF2 + 24.*eC1qg_2*gamma1gq_2[ii2]*(eC1qg_1 + gamma1qg_1[ii1]*LQF) - 12.*eC1qg_1*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQ - 24.*B1q*eC1qg_1*gamma1qq_2[ii2]*LQ + 36.*beta0*eC1qg_1*gamma1qq_2[ii2]*LQ - 12.*eC1qg_1*pow(gamma1qq_2[ii2],2)*LQ - 24.*B1q*eC1qg_1*gamma1gg_1[ii1]*LQ + 36.*beta0*eC1qg_1*gamma1gg_1[ii1]*LQ - 24.*eC1qg_1*gamma1qq_2[ii2]*gamma1gg_1[ii1]*LQ - 12.*eC1qg_1*pow(gamma1gg_1[ii1],2)*LQ - 24.*B2q*gamma1qg_1[ii1]*LQ + 12.*beta1*gamma1qg_1[ii1]*LQ - 24.*B1q*eC1qq_1*gamma1qg_1[ii1]*LQ + 36.*beta0*eC1qq_1*gamma1qg_1[ii1]*LQ - 24.*eC1qq_1*gamma1qq_2[ii2]*gamma1qg_1[ii1]*LQ - 12.*eC1qq_1*gamma1gg_1[ii1]*gamma1qg_1[ii1]*LQ - 12.*eC1qq_1*gamma1qg_1[ii1]*gamma1qq_1[ii1]*LQ - 24.*gamma1qg_1[ii1]*gamma2qq_2[ii2]*LQ - 12.*gamma1qg_1[ii1]*gamma2gg_1[ii1]*LQ - 24.*B1q*gamma2qg_1[ii1]*LQ + 24.*beta0*gamma2qg_1[ii1]*LQ - 24.*gamma1qq_2[ii2]*gamma2qg_1[ii1]*LQ - 12.*gamma1gg_1[ii1]*gamma2qg_1[ii1]*LQ - 12.*gamma1qq_1[ii1]*gamma2qg_1[ii1]*LQ - 12.*gamma1qg_1[ii1]*gamma2qq_1[ii1]*LQ - 12.*gamma1qg_1[ii1]*gamma2qqb_1[ii1]*LQ + 12.*gamma1qg_1[ii1]*gamma2qqbp_1[ii1]*LQ + 12.*gamma1qg_1[ii1]*gamma2qqp_1[ii1]*LQ - 12.*B1q*beta0*gamma1qg_1[ii1]*LF*LQ + 8.*pow(beta0,2)*gamma1qg_1[ii1]*LF*LQ - 24.*gamma1gq_2[ii2]*gamma1qg_2[ii2]*gamma1qg_1[ii1]*LF*LQ - 24.*B1q*gamma1qq_2[ii2]*gamma1qg_1[ii1]*LF* LQ - 24.*pow(gamma1qq_2[ii2],2)*gamma1qg_1[ii1]*LF*LQ - 12.*B1q*gamma1gg_1[ii1]*gamma1qg_1[ii1]*LF* LQ - 24.*gamma1qq_2[ii2]*gamma1gg_1[ii1]*gamma1qg_1[ii1]*LF*LQ - 8.*pow(gamma1gg_1[ii1],2)*gamma1qg_1[ii1]*LF*LQ - 12.*B1q*gamma1qg_1[ii1]*gamma1qq_1[ii1]*LF*LQ - 24.*gamma1qq_2[ii2]*gamma1qg_1[ii1]*gamma1qq_1[ii1]*LF*LQ - 8.*gamma1gg_1[ii1]*gamma1qg_1[ii1]* gamma1qq_1[ii1]*LF*LQ - 8.*gamma1qg_1[ii1]*pow(gamma1qq_1[ii1],2)*LF*LQ - 12.*A1q*eC1qg_1*gamma1qq_2[ii2]*LQ2 - 12.*A1q*eC1qg_1*gamma1gg_1[ii1]*LQ2 - 12.*A2q*gamma1qg_1[ii1]*LQ2 + 12.*pow(B1q,2)*gamma1qg_1[ii1]*LQ2 - 24.*B1q*beta0*gamma1qg_1[ii1]*LQ2 + 8.*pow(beta0,2)*gamma1qg_1[ii1]*LQ2 - 12.*A1q*eC1qq_1*gamma1qg_1[ii1]*LQ2 + 12.*gamma1gq_2[ii2]*gamma1qg_2[ii2]*gamma1qg_1[ii1]*LQ2 + 24.*B1q*gamma1qq_2[ii2]*gamma1qg_1[ii1]*LQ2 - 24.*beta0*gamma1qq_2[ii2]*gamma1qg_1[ii1]*LQ2 + 12.*pow(gamma1qq_2[ii2],2)*gamma1qg_1[ii1]*LQ2 + 12.*B1q*gamma1gg_1[ii1]*gamma1qg_1[ii1]*LQ2 - 12.*beta0*gamma1gg_1[ii1]*gamma1qg_1[ii1]*LQ2 + 12.*gamma1qq_2[ii2]*gamma1gg_1[ii1]*gamma1qg_1[ii1]* LQ2 + 4.*pow(gamma1gg_1[ii1],2)*gamma1qg_1[ii1]*LQ2 + 12.*B1q*gamma1qg_1[ii1]*gamma1qq_1[ii1]* LQ2 - 12.*beta0*gamma1qg_1[ii1]*gamma1qq_1[ii1]*LQ2 + 12.*gamma1qq_2[ii2]*gamma1qg_1[ii1]*gamma1qq_1[ii1]*LQ2 + 4.*gamma1gg_1[ii1]*gamma1qg_1[ii1]* gamma1qq_1[ii1]*LQ2 + 4.*gamma1qg_1[ii1]*pow(gamma1qq_1[ii1],2)*LQ2 - 12.*A1q*gamma2qg_1[ii1]*LQ2 - 6.*A1q*beta0*gamma1qg_1[ii1]*LF*LQ2 - 12.*A1q*gamma1qq_2[ii2]*gamma1qg_1[ii1]*LF*LQ2 - 6.*A1q*gamma1gg_1[ii1]*gamma1qg_1[ii1]*LF*LQ2 - 6.*A1q*gamma1qg_1[ii1]*gamma1qq_1[ii1]*LF*LQ2 + 12.*A1q*B1q*gamma1qg_1[ii1]*LQ3 - 14.*A1q*beta0*gamma1qg_1[ii1]*LQ3 + 12.*A1q*gamma1qq_2[ii2]*gamma1qg_1[ii1]*LQ3 + 6.*A1q*gamma1gg_1[ii1]*gamma1qg_1[ii1]*LQ3 + 6.*A1q*gamma1qg_1[ii1]*gamma1qq_1[ii1]*LQ3 + 3.*pow(A1q,2)*gamma1qg_1[ii1]*LQ4 + 12.*eC1qq_2*(2.*eC1qg_1*(gamma1qq_2[ii2] + gamma1gg_1[ii1]) + 2.*eC1qq_1*gamma1qg_1[ii1] + 2.*gamma2qg_1[ii1] + beta0*gamma1qg_1[ii1]*LF + 2.*gamma1qq_2[ii2]*gamma1qg_1[ii1]*LF + gamma1gg_1[ii1]*gamma1qg_1[ii1]*LF + gamma1qg_1[ii1]*gamma1qq_1[ii1]*LF - 2.*B1q*gamma1qg_1[ii1]*LQ + 3.*beta0*gamma1qg_1[ii1]*LQ - 2.*gamma1qq_2[ii2]*gamma1qg_1[ii1]*LQ - gamma1gg_1[ii1]*gamma1qg_1[ii1]*LQ - gamma1qg_1[ii1]*gamma1qq_1[ii1]*LQ - A1q*gamma1qg_1[ii1]*LQ2) + 24.*eC2qqbp_1*gamma1qg_1[ii1]*nf + 24.*eC2qqp_1*gamma1qg_1[ii1]*nf + 24.*eC1qg_1*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF*nf + 12.*gamma1qg_1[ii1]*gamma2qqbp_1[ii1]*LF*nf + 12.*gamma1qg_1[ii1]*gamma2qqp_1[ii1]*LF*nf + 8.*gamma1gq_1[ii1]*pow(gamma1qg_1[ii1],2)*LF2*nf - 24.*eC1qg_1*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQ*nf - 12.*gamma1qg_1[ii1]*gamma2qqbp_1[ii1]*LQ*nf - 12.*gamma1qg_1[ii1]*gamma2qqp_1[ii1]*LQ*nf - 16.*gamma1gq_1[ii1]*pow(gamma1qg_1[ii1],2)*LF*LQ*nf + 8.*gamma1gq_1[ii1]*pow(gamma1qg_1[ii1],2)*LQ2*nf))/24.;
              complex <double> H3st_qq_2   = (LQF*(12.*eC2qg_2*gamma1gq_2[ii2] + 12.*eC2qqb_2*(gamma1qq_1[ii1] + gamma1qq_2[ii2]) + 12.*eC1qg_2*gamma2gq_2[ii2] + 12.*eC1qq_2*gamma2qqb_2[ii2] + 12.*gamma3qqb_2[ii2] + 6.*beta0*eC1qg_2*gamma1gq_2[ii2]*LF + 12.*eC1qg_2*gamma1qq_1[ii1]*gamma1gq_2[ii2]*LF + 6.*eC1qg_2*gamma1gg_2[ii2]*gamma1gq_2[ii2]*LF + 6.*eC1qq_2*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF + 6.*eC1qg_2*gamma1gq_2[ii2]*gamma1qq_2[ii2]*LF + 6.*gamma1qg_2[ii2]*gamma2gq_2[ii2]*LF + 6.*gamma1gq_2[ii2]*gamma2qg_2[ii2]*LF + 12.*beta0*gamma2qqb_2[ii2]*LF + 12.*gamma1qq_1[ii1]*gamma2qqb_2[ii2]*LF + 12.*gamma1qq_2[ii2]*gamma2qqb_2[ii2]*LF + 6.*beta0*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF2 + 6.*gamma1qq_1[ii1]*gamma1gq_2[ii2]*gamma1qg_2[ii2]* LF2 + 2.*gamma1gg_2[ii2]*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF2 + 4.*gamma1gq_2[ii2]*gamma1qg_2[ii2]*gamma1qq_2[ii2]*LF2 + 6.*eC1qq_1*(2.*eC1qg_2*gamma1gq_2[ii2] + 2.*gamma2qqb_2[ii2] + gamma1gq_2[ii2]*gamma1qg_2[ii2]* LQF) - 12.*B1q*eC1qg_2*gamma1gq_2[ii2]*LQ + 18.*beta0*eC1qg_2*gamma1gq_2[ii2]* LQ - 12.*eC1qg_2*gamma1qq_1[ii1]*gamma1gq_2[ii2]*LQ - 6.*eC1qg_2*gamma1gg_2[ii2]*gamma1gq_2[ii2]* LQ - 6.*eC1qq_2*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQ - 6.*eC1qg_2*gamma1gq_2[ii2]*gamma1qq_2[ii2]* LQ - 6.*gamma1qg_2[ii2]*gamma2gq_2[ii2]*LQ - 6.*gamma1gq_2[ii2]*gamma2qg_2[ii2]*LQ - 12.*B1q*gamma2qqb_2[ii2]*LQ + 12.*beta0*gamma2qqb_2[ii2]*LQ - 12.*gamma1qq_1[ii1]*gamma2qqb_2[ii2]*LQ - 12.*gamma1qq_2[ii2]*gamma2qqb_2[ii2]*LQ - 6.*B1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF*LQ - 12.*gamma1qq_1[ii1]*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF* LQ - 4.*gamma1gg_2[ii2]*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF*LQ - 8.*gamma1gq_2[ii2]*gamma1qg_2[ii2]*gamma1qq_2[ii2]*LF*LQ - 6.*A1q*eC1qg_2*gamma1gq_2[ii2]*LQ2 + 6.*B1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQ2 - 6.*beta0*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQ2 + 6.*gamma1qq_1[ii1]*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQ2 + 2.*gamma1gg_2[ii2]*gamma1gq_2[ii2]* gamma1qg_2[ii2]*LQ2 + 4.*gamma1gq_2[ii2]*gamma1qg_2[ii2]*gamma1qq_2[ii2]*LQ2 - 6.*A1q*gamma2qqb_2[ii2]*LQ2 - 3.*A1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF*LQ2 + 3.*A1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQ3))/12.;
              complex <double> H3st_qqp_2  = (LQF*(12.*eC2qg_2*gamma1gq_2[ii2] + 12.*eC2qqbp_2*(gamma1qq_1[ii1] + gamma1qq_2[ii2]) + 12.*eC1qg_2*gamma2gq_2[ii2] + 12.*eC1qq_2*gamma2qqbp_2[ii2] + 12.*gamma3qqbp_2[ii2] + 6.*beta0*eC1qg_2*gamma1gq_2[ii2]*LF + 12.*eC1qg_2*gamma1qq_1[ii1]*gamma1gq_2[ii2]*LF + 6.*eC1qg_2*gamma1gg_2[ii2]*gamma1gq_2[ii2]*LF + 6.*eC1qq_2*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF + 6.*eC1qg_2*gamma1gq_2[ii2]*gamma1qq_2[ii2]*LF + 6.*gamma1qg_2[ii2]*gamma2gq_2[ii2]*LF + 6.*gamma1gq_2[ii2]*gamma2qg_2[ii2]*LF + 12.*beta0*gamma2qqbp_2[ii2]*LF + 12.*gamma1qq_1[ii1]*gamma2qqbp_2[ii2]*LF + 12.*gamma1qq_2[ii2]*gamma2qqbp_2[ii2]*LF + 6.*beta0*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF2 + 6.*gamma1qq_1[ii1]*gamma1gq_2[ii2]*gamma1qg_2[ii2]* LF2 + 2.*gamma1gg_2[ii2]*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF2 + 4.*gamma1gq_2[ii2]*gamma1qg_2[ii2]*gamma1qq_2[ii2]*LF2 + 6.*eC1qq_1*(2.*eC1qg_2*gamma1gq_2[ii2] + 2.*gamma2qqbp_2[ii2] + gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQF) - 12.*B1q*eC1qg_2*gamma1gq_2[ii2]*LQ + 18.*beta0*eC1qg_2*gamma1gq_2[ii2]*LQ - 12.*eC1qg_2*gamma1qq_1[ii1]*gamma1gq_2[ii2]*LQ - 6.*eC1qg_2*gamma1gg_2[ii2]*gamma1gq_2[ii2]*LQ - 6.*eC1qq_2*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQ - 6.*eC1qg_2*gamma1gq_2[ii2]*gamma1qq_2[ii2]*LQ - 6.*gamma1qg_2[ii2]*gamma2gq_2[ii2]*LQ - 6.*gamma1gq_2[ii2]*gamma2qg_2[ii2]*LQ - 12.*B1q*gamma2qqbp_2[ii2]*LQ + 12.*beta0*gamma2qqbp_2[ii2]*LQ - 12.*gamma1qq_1[ii1]*gamma2qqbp_2[ii2]*LQ - 12.*gamma1qq_2[ii2]*gamma2qqbp_2[ii2]*LQ - 6.*B1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF*LQ - 12.*gamma1qq_1[ii1]*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF*LQ - 4.*gamma1gg_2[ii2]*gamma1gq_2[ii2]* gamma1qg_2[ii2]*LF*LQ - 8.*gamma1gq_2[ii2]*gamma1qg_2[ii2]*gamma1qq_2[ii2]*LF*LQ - 6.*A1q*eC1qg_2*gamma1gq_2[ii2]*LQ2 + 6.*B1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQ2 - 6.*beta0*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQ2 + 6.*gamma1qq_1[ii1]*gamma1gq_2[ii2]*gamma1qg_2[ii2]* LQ2 + 2.*gamma1gg_2[ii2]*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQ2 + 4.*gamma1gq_2[ii2]*gamma1qg_2[ii2]*gamma1qq_2[ii2]*LQ2 - 6.*A1q*gamma2qqbp_2[ii2]*LQ2 - 3.*A1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF*LQ2 + 3.*A1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQ3))/12.;
              complex <double> H3st_qqbp_2 = (LQF*(12.*eC2qg_2*gamma1gq_2[ii2] + 12.*eC2qqp_2*(gamma1qq_1[ii1] + gamma1qq_2[ii2]) + 12.*eC1qg_2*gamma2gq_2[ii2] + 12.*eC1qq_2*gamma2qqp_2[ii2] + 12.*gamma3qqp_2[ii2] + 6.*beta0*eC1qg_2*gamma1gq_2[ii2]*LF + 12.*eC1qg_2*gamma1qq_1[ii1]*gamma1gq_2[ii2]*LF + 6.*eC1qg_2*gamma1gg_2[ii2]*gamma1gq_2[ii2]*LF + 6.*eC1qq_2*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF + 6.*eC1qg_2*gamma1gq_2[ii2]*gamma1qq_2[ii2]*LF + 6.*gamma1qg_2[ii2]*gamma2gq_2[ii2]*LF + 6.*gamma1gq_2[ii2]*gamma2qg_2[ii2]*LF + 12.*beta0*gamma2qqp_2[ii2]*LF + 12.*gamma1qq_1[ii1]*gamma2qqp_2[ii2]*LF + 12.*gamma1qq_2[ii2]*gamma2qqp_2[ii2]*LF + 6.*beta0*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF2 + 6.*gamma1qq_1[ii1]*gamma1gq_2[ii2]*gamma1qg_2[ii2]* LF2 + 2.*gamma1gg_2[ii2]*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF2 + 4.*gamma1gq_2[ii2]*gamma1qg_2[ii2]*gamma1qq_2[ii2]*LF2 + 6.*eC1qq_1*(2.*eC1qg_2*gamma1gq_2[ii2] + 2.*gamma2qqp_2[ii2] + gamma1gq_2[ii2]*gamma1qg_2[ii2]* LQF) - 12.*B1q*eC1qg_2*gamma1gq_2[ii2]*LQ + 18.*beta0*eC1qg_2*gamma1gq_2[ii2]* LQ - 12.*eC1qg_2*gamma1qq_1[ii1]*gamma1gq_2[ii2]*LQ - 6.*eC1qg_2*gamma1gg_2[ii2]*gamma1gq_2[ii2]* LQ - 6.*eC1qq_2*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQ - 6.*eC1qg_2*gamma1gq_2[ii2]*gamma1qq_2[ii2]* LQ - 6.*gamma1qg_2[ii2]*gamma2gq_2[ii2]*LQ - 6.*gamma1gq_2[ii2]*gamma2qg_2[ii2]*LQ - 12.*B1q*gamma2qqp_2[ii2]*LQ + 12.*beta0*gamma2qqp_2[ii2]*LQ - 12.*gamma1qq_1[ii1]*gamma2qqp_2[ii2]*LQ - 12.*gamma1qq_2[ii2]*gamma2qqp_2[ii2]*LQ - 6.*B1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF*LQ - 12.*gamma1qq_1[ii1]*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF* LQ - 4.*gamma1gg_2[ii2]*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF*LQ - 8.*gamma1gq_2[ii2]*gamma1qg_2[ii2]*gamma1qq_2[ii2]*LF*LQ - 6.*A1q*eC1qg_2*gamma1gq_2[ii2]*LQ2 + 6.*B1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQ2 - 6.*beta0*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQ2 + 6.*gamma1qq_1[ii1]*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQ2 + 2.*gamma1gg_2[ii2]*gamma1gq_2[ii2]* gamma1qg_2[ii2]*LQ2 + 4.*gamma1gq_2[ii2]*gamma1qg_2[ii2]*gamma1qq_2[ii2]*LQ2 - 6.*A1q*gamma2qqp_2[ii2]*LQ2 - 3.*A1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF*LQ2 + 3.*A1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQ3))/12.;
              complex <double> H3st_qg_2   = (LQF*(24.*eC2qg_2*gamma1qq_1[ii1] + 24.*eC2qg_2*gamma1gg_2[ii2] + 24.*eC2qq_1*gamma1qg_2[ii2] + 24.*eC2qq_2*gamma1qg_2[ii2] + 24.*eC2qqb_2*gamma1qg_2[ii2] - 24.*eC2qqbp_2*gamma1qg_2[ii2] - 24.*eC2qqp_2*gamma1qg_2[ii2] + 24.*eC1qg_2*gamma2qq_1[ii1] + 24.*eC1qg_2*gamma2gg_2[ii2] + 24.*eC1qq_2*gamma2qg_2[ii2] + 24.*gamma3qg_2[ii2] + 12.*eC1qg_2*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF + 12.*beta0*eC1qg_2*gamma1qq_1[ii1]*LF + 12.*eC1qg_2*pow(gamma1qq_1[ii1],2)*LF + 12.*beta0*eC1qg_2*gamma1gg_2[ii2]*LF + 24.*eC1qg_2*gamma1qq_1[ii1]*gamma1gg_2[ii2]*LF + 12.*eC1qg_2*pow(gamma1gg_2[ii2],2)*LF + 12.*beta1*gamma1qg_2[ii2]*LF + 12.*beta0*eC1qq_2*gamma1qg_2[ii2]*LF + 24.*eC1qq_2*gamma1qq_1[ii1]*gamma1qg_2[ii2]*LF + 12.*eC1qq_2*gamma1gg_2[ii2]*gamma1qg_2[ii2]*LF + 12.*eC1qq_2*gamma1qg_2[ii2]*gamma1qq_2[ii2]*LF + 24.*gamma1qg_2[ii2]*gamma2qq_1[ii1]*LF + 12.*gamma1qg_2[ii2]*gamma2gg_2[ii2]*LF + 24.*beta0*gamma2qg_2[ii2]*LF + 24.*gamma1qq_1[ii1]*gamma2qg_2[ii2]*LF + 12.*gamma1gg_2[ii2]*gamma2qg_2[ii2]*LF + 12.*gamma1qq_2[ii2]*gamma2qg_2[ii2]*LF + 12.*gamma1qg_2[ii2]*gamma2qq_2[ii2]*LF + 12.*gamma1qg_2[ii2]*gamma2qqb_2[ii2]*LF - 12.*gamma1qg_2[ii2]*gamma2qqbp_2[ii2]*LF - 12.*gamma1qg_2[ii2]*gamma2qqp_2[ii2]*LF + 8.*pow(beta0,2)*gamma1qg_2[ii2]*LF2 + 12.*gamma1gq_1[ii1]*gamma1qg_1[ii1]*gamma1qg_2[ii2]*LF2 + 24.*beta0*gamma1qq_1[ii1]*gamma1qg_2[ii2]* LF2 + 12.*pow(gamma1qq_1[ii1],2)*gamma1qg_2[ii2]*LF2 + 12.*beta0*gamma1gg_2[ii2]*gamma1qg_2[ii2]* LF2 + 12.*gamma1qq_1[ii1]*gamma1gg_2[ii2]*gamma1qg_2[ii2]*LF2 + 4.*pow(gamma1gg_2[ii2],2)*gamma1qg_2[ii2]*LF2 + 12.*beta0*gamma1qg_2[ii2]*gamma1qq_2[ii2]*LF2 + 12.*gamma1qq_1[ii1]*gamma1qg_2[ii2]*gamma1qq_2[ii2]*LF2 + 4.*gamma1gg_2[ii2]*gamma1qg_2[ii2]* gamma1qq_2[ii2]*LF2 + 4.*gamma1qg_2[ii2]*pow(gamma1qq_2[ii2],2)*LF2 + 24.*eC1qg_1*gamma1gq_1[ii1]*(eC1qg_2 + gamma1qg_2[ii2]*LQF) - 12.*eC1qg_2*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQ - 24.*B1q*eC1qg_2*gamma1qq_1[ii1]*LQ + 36.*beta0*eC1qg_2*gamma1qq_1[ii1]*LQ - 12.*eC1qg_2*pow(gamma1qq_1[ii1],2)*LQ - 24.*B1q*eC1qg_2*gamma1gg_2[ii2]*LQ + 36.*beta0*eC1qg_2*gamma1gg_2[ii2]*LQ - 24.*eC1qg_2*gamma1qq_1[ii1]*gamma1gg_2[ii2]*LQ - 12.*eC1qg_2*pow(gamma1gg_2[ii2],2)*LQ - 24.*B2q*gamma1qg_2[ii2]*LQ + 12.*beta1*gamma1qg_2[ii2]*LQ - 24.*B1q*eC1qq_2*gamma1qg_2[ii2]*LQ + 36.*beta0*eC1qq_2*gamma1qg_2[ii2]*LQ - 24.*eC1qq_2*gamma1qq_1[ii1]*gamma1qg_2[ii2]*LQ - 12.*eC1qq_2*gamma1gg_2[ii2]*gamma1qg_2[ii2]*LQ - 12.*eC1qq_2*gamma1qg_2[ii2]*gamma1qq_2[ii2]*LQ - 24.*gamma1qg_2[ii2]*gamma2qq_1[ii1]*LQ - 12.*gamma1qg_2[ii2]*gamma2gg_2[ii2]*LQ - 24.*B1q*gamma2qg_2[ii2]*LQ + 24.*beta0*gamma2qg_2[ii2]*LQ - 24.*gamma1qq_1[ii1]*gamma2qg_2[ii2]*LQ - 12.*gamma1gg_2[ii2]*gamma2qg_2[ii2]*LQ - 12.*gamma1qq_2[ii2]*gamma2qg_2[ii2]*LQ - 12.*gamma1qg_2[ii2]*gamma2qq_2[ii2]*LQ - 12.*gamma1qg_2[ii2]*gamma2qqb_2[ii2]*LQ + 12.*gamma1qg_2[ii2]*gamma2qqbp_2[ii2]*LQ + 12.*gamma1qg_2[ii2]*gamma2qqp_2[ii2]*LQ - 12.*B1q*beta0*gamma1qg_2[ii2]*LF*LQ + 8.*pow(beta0,2)*gamma1qg_2[ii2]*LF*LQ - 24.*gamma1gq_1[ii1]*gamma1qg_1[ii1]*gamma1qg_2[ii2]*LF*LQ - 24.*B1q*gamma1qq_1[ii1]*gamma1qg_2[ii2]*LF* LQ - 24.*pow(gamma1qq_1[ii1],2)*gamma1qg_2[ii2]*LF*LQ - 12.*B1q*gamma1gg_2[ii2]*gamma1qg_2[ii2]*LF* LQ - 24.*gamma1qq_1[ii1]*gamma1gg_2[ii2]*gamma1qg_2[ii2]*LF*LQ - 8.*pow(gamma1gg_2[ii2],2)*gamma1qg_2[ii2]*LF*LQ - 12.*B1q*gamma1qg_2[ii2]*gamma1qq_2[ii2]*LF*LQ - 24.*gamma1qq_1[ii1]*gamma1qg_2[ii2]*gamma1qq_2[ii2]*LF*LQ - 8.*gamma1gg_2[ii2]*gamma1qg_2[ii2]* gamma1qq_2[ii2]*LF*LQ - 8.*gamma1qg_2[ii2]*pow(gamma1qq_2[ii2],2)*LF*LQ - 12.*A1q*eC1qg_2*gamma1qq_1[ii1]*LQ2 - 12.*A1q*eC1qg_2*gamma1gg_2[ii2]*LQ2 - 12.*A2q*gamma1qg_2[ii2]*LQ2 + 12.*pow(B1q,2)*gamma1qg_2[ii2]*LQ2 - 24.*B1q*beta0*gamma1qg_2[ii2]*LQ2 + 8.*pow(beta0,2)*gamma1qg_2[ii2]*LQ2 - 12.*A1q*eC1qq_2*gamma1qg_2[ii2]*LQ2 + 12.*gamma1gq_1[ii1]*gamma1qg_1[ii1]*gamma1qg_2[ii2]*LQ2 + 24.*B1q*gamma1qq_1[ii1]*gamma1qg_2[ii2]*LQ2 - 24.*beta0*gamma1qq_1[ii1]*gamma1qg_2[ii2]*LQ2 + 12.*pow(gamma1qq_1[ii1],2)*gamma1qg_2[ii2]*LQ2 + 12.*B1q*gamma1gg_2[ii2]*gamma1qg_2[ii2]*LQ2 - 12.*beta0*gamma1gg_2[ii2]*gamma1qg_2[ii2]*LQ2 + 12.*gamma1qq_1[ii1]*gamma1gg_2[ii2]*gamma1qg_2[ii2]* LQ2 + 4.*pow(gamma1gg_2[ii2],2)*gamma1qg_2[ii2]*LQ2 + 12.*B1q*gamma1qg_2[ii2]*gamma1qq_2[ii2]* LQ2 - 12.*beta0*gamma1qg_2[ii2]*gamma1qq_2[ii2]*LQ2 + 12.*gamma1qq_1[ii1]*gamma1qg_2[ii2]*gamma1qq_2[ii2]*LQ2 + 4.*gamma1gg_2[ii2]*gamma1qg_2[ii2]* gamma1qq_2[ii2]*LQ2 + 4.*gamma1qg_2[ii2]*pow(gamma1qq_2[ii2],2)*LQ2 - 12.*A1q*gamma2qg_2[ii2]*LQ2 - 6.*A1q*beta0*gamma1qg_2[ii2]*LF*LQ2 - 12.*A1q*gamma1qq_1[ii1]*gamma1qg_2[ii2]*LF*LQ2 - 6.*A1q*gamma1gg_2[ii2]*gamma1qg_2[ii2]*LF*LQ2 - 6.*A1q*gamma1qg_2[ii2]*gamma1qq_2[ii2]*LF*LQ2 + 12.*A1q*B1q*gamma1qg_2[ii2]*LQ3 - 14.*A1q*beta0*gamma1qg_2[ii2]*LQ3 + 12.*A1q*gamma1qq_1[ii1]*gamma1qg_2[ii2]*LQ3 + 6.*A1q*gamma1gg_2[ii2]*gamma1qg_2[ii2]*LQ3 + 6.*A1q*gamma1qg_2[ii2]*gamma1qq_2[ii2]*LQ3 + 3.*pow(A1q,2)*gamma1qg_2[ii2]*LQ4 + 12.*eC1qq_1*(2.*eC1qg_2*(gamma1qq_1[ii1] + gamma1gg_2[ii2]) + 2.*eC1qq_2*gamma1qg_2[ii2] + 2.*gamma2qg_2[ii2] + beta0*gamma1qg_2[ii2]*LF + 2.*gamma1qq_1[ii1]*gamma1qg_2[ii2]*LF + gamma1gg_2[ii2]*gamma1qg_2[ii2]*LF + gamma1qg_2[ii2]*gamma1qq_2[ii2]*LF - 2.*B1q*gamma1qg_2[ii2]*LQ + 3.*beta0*gamma1qg_2[ii2]*LQ - 2.*gamma1qq_1[ii1]*gamma1qg_2[ii2]*LQ - gamma1gg_2[ii2]*gamma1qg_2[ii2]*LQ - gamma1qg_2[ii2]*gamma1qq_2[ii2]*LQ - A1q*gamma1qg_2[ii2]*LQ2) + 24.*eC2qqbp_2*gamma1qg_2[ii2]*nf + 24.*eC2qqp_2*gamma1qg_2[ii2]*nf + 24.*eC1qg_2*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF*nf + 12.*gamma1qg_2[ii2]*gamma2qqbp_2[ii2]*LF*nf + 12.*gamma1qg_2[ii2]*gamma2qqp_2[ii2]*LF*nf + 8.*gamma1gq_2[ii2]*pow(gamma1qg_2[ii2],2)*LF2*nf - 24.*eC1qg_2*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQ*nf - 12.*gamma1qg_2[ii2]*gamma2qqbp_2[ii2]*LQ*nf - 12.*gamma1qg_2[ii2]*gamma2qqp_2[ii2]*LQ*nf - 16.*gamma1gq_2[ii2]*pow(gamma1qg_2[ii2],2)*LF*LQ*nf + 8.*gamma1gq_2[ii2]*pow(gamma1qg_2[ii2],2)*LQ2*nf))/24.;

              complex <double> H3st_gg     = (LQF*(2.*eC2qg_2*gamma1qg_1[ii1] + 2.*eC2qg_1*gamma1qg_2[ii2] + 2.*eC1qg_2*gamma2qg_1[ii1] + beta0*eC1qg_2*gamma1qg_1[ii1]*LF + eC1qg_2*gamma1gg_1[ii1]*gamma1qg_1[ii1]*LF + eC1qg_2*gamma1qg_1[ii1]*gamma1qq_1[ii1]*LF + 2.*eC1qg_2*gamma1qg_1[ii1]*gamma1gg_2[ii2]*LF + 2.*eC1qq_2*gamma1qg_1[ii1]*gamma1qg_2[ii2]*LF + 2.*gamma1qg_2[ii2]*gamma2qg_1[ii1]*LF + 2.*gamma1qg_1[ii1]*gamma2qg_2[ii2]*LF + 2.*beta0*gamma1qg_1[ii1]*gamma1qg_2[ii2]*LF2 + gamma1gg_1[ii1]*gamma1qg_1[ii1]*gamma1qg_2[ii2]*LF2 + gamma1qg_1[ii1]*gamma1qq_1[ii1]*gamma1qg_2[ii2]* LF2 + gamma1qg_1[ii1]*gamma1gg_2[ii2]*gamma1qg_2[ii2]*LF2 + gamma1qg_1[ii1]*gamma1qg_2[ii2]*gamma1qq_2[ii2]*LF2 - 2.*B1q*eC1qg_2*gamma1qg_1[ii1]*LQ + 3.*beta0*eC1qg_2*gamma1qg_1[ii1]*LQ - eC1qg_2*gamma1gg_1[ii1]*gamma1qg_1[ii1]*LQ - eC1qg_2*gamma1qg_1[ii1]*gamma1qq_1[ii1]*LQ - 2.*eC1qg_2*gamma1qg_1[ii1]*gamma1gg_2[ii2]*LQ - 2.*eC1qq_2*gamma1qg_1[ii1]*gamma1qg_2[ii2]*LQ - 2.*gamma1qg_2[ii2]*gamma2qg_1[ii1]*LQ - 2.*gamma1qg_1[ii1]*gamma2qg_2[ii2]*LQ - 2.*B1q*gamma1qg_1[ii1]*gamma1qg_2[ii2]*LF*LQ - 2.*gamma1gg_1[ii1]*gamma1qg_1[ii1]*gamma1qg_2[ii2]*LF*LQ - 2.*gamma1qg_1[ii1]*gamma1qq_1[ii1]* gamma1qg_2[ii2]*LF*LQ - 2.*gamma1qg_1[ii1]*gamma1gg_2[ii2]*gamma1qg_2[ii2]*LF*LQ - 2.*gamma1qg_1[ii1]*gamma1qg_2[ii2]*gamma1qq_2[ii2]*LF*LQ - A1q*eC1qg_2*gamma1qg_1[ii1]*LQ2 + 2.*B1q*gamma1qg_1[ii1]*gamma1qg_2[ii2]*LQ2 - 2.*beta0*gamma1qg_1[ii1]*gamma1qg_2[ii2]*LQ2 + gamma1gg_1[ii1]*gamma1qg_1[ii1]*gamma1qg_2[ii2]*LQ2 + gamma1qg_1[ii1]*gamma1qq_1[ii1]*gamma1qg_2[ii2]* LQ2 + gamma1qg_1[ii1]*gamma1gg_2[ii2]*gamma1qg_2[ii2]*LQ2 + gamma1qg_1[ii1]*gamma1qg_2[ii2]*gamma1qq_2[ii2]*LQ2 - A1q*gamma1qg_1[ii1]*gamma1qg_2[ii2]*LF*LQ2 + A1q*gamma1qg_1[ii1]*gamma1qg_2[ii2]*LQ3 + 2.*eC1qq_1*gamma1qg_1[ii1]* (eC1qg_2 + gamma1qg_2[ii2]*LF - gamma1qg_2[ii2]*LQ) + eC1qg_1*(2.*eC1qg_2*(gamma1gg_1[ii1] + gamma1gg_2[ii2]) + 2.*eC1qq_2*gamma1qg_2[ii2] + 2.*gamma2qg_2[ii2] + beta0*gamma1qg_2[ii2]*LF + 2.*gamma1gg_1[ii1]*gamma1qg_2[ii2]*LF + gamma1gg_2[ii2]*gamma1qg_2[ii2]*LF + gamma1qg_2[ii2]*gamma1qq_2[ii2]*LF - 2.*B1q*gamma1qg_2[ii2]*LQ + 3.*beta0*gamma1qg_2[ii2]*LQ - 2.*gamma1gg_1[ii1]*gamma1qg_2[ii2]* LQ - gamma1gg_2[ii2]*gamma1qg_2[ii2]*LQ - gamma1qg_2[ii2]*gamma1qq_2[ii2]*LQ - A1q*gamma1qg_2[ii2]*LQ2)))/2.;
              complex <double> H3st_qbg_1  = ((2.*eC2qqb_2*gamma1qg_1[ii1] + (2.*eC1qg_2*gamma1gq_2[ii2] + 2.*gamma2qqb_2[ii2] + gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQF)*(eC1qg_1 + gamma1qg_1[ii1]*LQF))* LQF)/2.;
              complex <double> H3st_qbg_2  = ((2.*eC2qqb_1*gamma1qg_2[ii2] + (2.*eC1qg_1*gamma1gq_1[ii1] + 2.*gamma2qqb_1[ii1] + gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQF)*(eC1qg_2 + gamma1qg_2[ii2]*LQF))* LQF)/2.;
              complex <double> H3st_qpg_1  = ((2.*eC2qqp_2*gamma1qg_1[ii1] + (2.*eC1qg_2*gamma1gq_2[ii2] + 2.*gamma2qqp_2[ii2] + gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQF)*(eC1qg_1 + gamma1qg_1[ii1]*LQF))* LQF)/2.;
              complex <double> H3st_qpg_2  = ((2.*eC2qqp_1*gamma1qg_2[ii2] + (2.*eC1qg_1*gamma1gq_1[ii1] + 2.*gamma2qqp_1[ii1] + gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQF)*(eC1qg_2 + gamma1qg_2[ii2]*LQF))* LQF)/2.;
              complex <double> H3st_qbpg_1 = ((2.*eC2qqbp_2*gamma1qg_1[ii1] + (2.*eC1qg_2*gamma1gq_2[ii2] + 2.*gamma2qqbp_2[ii2] + gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQF)*(eC1qg_1 + gamma1qg_1[ii1]*LQF))* LQF)/2.;
              complex <double> H3st_qbpg_2 = ((2.*eC2qqbp_1*gamma1qg_2[ii2] + (2.*eC1qg_1*gamma1gq_1[ii1] + 2.*gamma2qqbp_1[ii1] + gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQF)*(eC1qg_2 + gamma1qg_2[ii2]*LQF))* LQF)/2.;

              //Mixed LR*LQF variations
              H3st_qqb    += -2*beta0*H2st_qqb*LR + (-beta1*LR - beta0*LR2)*H1st_qqb;
              H3st_qq_1   += -2*beta0*H2st_qq_1*LR;
              H3st_qq_2   += -2*beta0*H2st_qq_2*LR;
              H3st_qqp_1  += -2*beta0*H2st_qqp_1*LR;
              H3st_qqp_2  += -2*beta0*H2st_qqp_1*LR;
              H3st_qqbp_1 += -2*beta0*H2st_qqbp_1*LR;
              H3st_qqbp_2 += -2*beta0*H2st_qqbp_2*LR;
              H3st_qg_1   += -2*beta0*H2st_qg_1*LR + (-beta1*LR - beta0*LR2)*H1st_qg_1;
              H3st_qg_2   += -2*beta0*H2st_qg_2*LR + (-beta1*LR - beta0*LR2)*H1st_qg_2;
              H3st_gg     += -2*beta0*H2st_gg*LR;

              qqb[idx]    += as3*H3st_qqb;
              qg_1[idx]   += as3*H3st_qg_1;
              qg_2[idx]   += as3*H3st_qg_2;
              qq_1[idx]   += as3*H3st_qq_1;
              qq_2[idx]   += as3*H3st_qq_2;
              qqp_1[idx]  += as3*H3st_qqp_1;
              qqp_2[idx]  += as3*H3st_qqp_2;
              qqbp_1[idx] += as3*H3st_qqbp_1;
              qqbp_2[idx] += as3*H3st_qqbp_2;
              gg[idx]     += as3*H3st_gg;
              qbg_1[idx]  += as3*H3st_qbg_1;
              qbg_2[idx]  += as3*H3st_qbg_2;
              qpg_1[idx]  += as3*H3st_qpg_1;
              qpg_2[idx]  += as3*H3st_qpg_2;
              qbpg_1[idx] += as3*H3st_qbpg_1;
              qbpg_2[idx] += as3*H3st_qbpg_2;
            }
}

