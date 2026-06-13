#include "Integral_class.h"
#include "hep-rng.h"
#include "rng-emp.h"
#include "Evento.h"
#include "hhc2ph-event.h"
#include "hhc2ph-phasespace.h"
#include "hhc2ph-process.h"
#include "hhc2ph1jet.h"
#include "Inter_params.h"

// DYturbo includes
#include "settings.h"
#include "coupling.h"
#include "phasespace.h"

#include <cmath>

// --------------------------------------------------------------------------
//// REAL NLO qgINTg
// --------------------------------------------------------------------------
double Integral_class::Integrando_qgintg(const double x[],event_hhc2ph &pmed)
{  
  _Evento Evnt;
  // -----------------------------------------------------------------
  // -- Unused parameters  
  // -----------------------------------------------------------------
  unsigned int nu=2,nd=3;
  rng_emp rng( -238,8);
  process_hhc2ph *_M_proc = new hhc2ph1jet(rng, true, nu, nd, 0.1);
  // -----------------------------------------------------------------
  event_hhc2ph p2(1);
  Evnt.Set_ps(Evnt);
  double s=Evnt._M_s,q2min=Evnt.q2min,q2max=Evnt.q2max;
  basic_phasespace_hhc2ph * _M_ps = new basic_phasespace_hhc2ph(&rng, s);
  _M_ps -> important_sampling(2.0, 0.1);
  _M_ps -> phasespace_cuts(q2min, q2max);
  // -----------------------------------------------------------------
  double weightqgINTg=0.0,weight2=0.0;
  int flag,flagr;double smearing,Q2,fglR,FacNew;  
  smearing = opts.SigmaSmearing; 
  bool flagrmed = false;
  double weightmed = 0.0;
  try {
    // double operator()(_Event&,double&,const double x[]);
    // weight2 = _M_ps->operator()(p2,smearing,x); // -- Samplimg Ok for Mgamgam also // qT behaviour OK
    if (opts.withsmearing) 
      {
	// weight2   = _M_ps->operator()(p2,Q2,flagrmed,x); // -- Sampling Ok for Mgam gam -- PhaseSpece Daniel
	weight2   = _M_ps->operator()(p2,flagrmed,x); // -- Sampling Ok for Mgam gam -- MIO
        // -- flag med activated
        flagrmed  = true;
	// weightmed = _M_ps->operator()(pmed,Q2,flagrmed ,x); // -- Sampling Ok for Mgam gam -- PhaseSpece Daniel    
	weightmed = _M_ps->operator()(pmed,flagrmed ,x); // -- Sampling Ok for Mgam gam -- MIO
      }
    else
      {
	//	weight2 =_M_ps->operator()(p2,Q2,flagrmed ,x); // -- Sampling Ok for Mgam gam -- PhaseSpece Daniel 
	weight2 =_M_ps->operator()(p2,flagrmed ,x); // -- Sampling Ok for Mgam gam -- MIO 
        pmed = p2;
      }
    //    weight2 = _M_ps->operator()(p2,Q2,smearing,x); // -- Sampling Ok for Mgam gam
    // weight2 = _M_ps->operator()(p2,fglR,FacNew,x);
  } catch(numeric_error) {
    delete _M_ps;
    delete _M_proc;
    return 0.0;
  } catch(fp_exception) {
    delete _M_ps;
    delete _M_proc;
    return 0.0;
  } catch(...) {delete _M_ps;
    delete _M_proc;
    std::cerr<<"integral<..., 2, 0> : It was an unexpected error\n"
	     <<"while the phase space was being generated."<<std::endl;
    throw;
  }
  amplitude_hhc2ph amp(_M_proc,amplitude_hhc2ph::lo);  
  //----- real term NLO -----------------
  amp.qgINTg_contrib(weight2, p2);
  userfunc(p2,weightqgINTg, amp,4.0,4.0);
  //--------------------------------------
  delete _M_ps;
  delete _M_proc;
  return weightqgINTg;  
}

//--------------------------------------------------------------------------
//----    CT NLO ggINTg
//--------------------------------------------------------------------------

double Integral_class::Integrando_CTgg(const double x[],event_hhc2ph &pmed)
{
  _Evento Evnt;
  // -----------------------------------------------------------------
  // -- Unused parameters  
  // -----------------------------------------------------------------
  unsigned int nu=2,nd=3;
  rng_emp rng( -238,8);
  process_hhc2ph *_M_proc = new hhc2ph1jet(rng, true, nu, nd, 0.1);
  // -----------------------------------------------------------------
  event_hhc2ph p1(1);
  Evnt.Set_ps(Evnt);
  double s=Evnt._M_s,q2min=Evnt.q2min,q2max=Evnt.q2max;
  basic_phasespace_hhc2ph * _M_ps = new basic_phasespace_hhc2ph(&rng, s);
  _M_ps -> important_sampling(2.0, 0.1);
  _M_ps -> phasespace_cuts(q2min, q2max);
// -----------------------------------------------------------------
  double weight1=0.0,weightCTgg=0.0,weightH1=0.0,weightH1f=0.0;
  double smearing;
  bool flagrmed = false;
  double Q2,z1,z2;
  double weightmul = 0.0;
  try {  
    if(opts.withsmearing)
      { 
        flagrmed = true;
	weightmul = _M_ps->operator()(pmed,Q2,smearing,flagrmed,z1,z2,x);    
	flagrmed = false;
	weight1   = _M_ps->operator()(p1,Q2,smearing,flagrmed,z1,z2,x);    
      }
    else
      {
	weight1   = _M_ps->operator()(p1,Q2,smearing,flagrmed,z1,z2,x);
	pmed = p1;
      }

  } catch(numeric_error) {   
    delete _M_ps;
    delete _M_proc;
    return 0.0;
  } catch(fp_exception) {
    delete _M_ps;
    delete _M_proc;
    return 0.0;
  } catch(...) {delete _M_ps;
    delete _M_proc;
    std::cerr<<"integral<..., 2, 0> : It was an unexpected error\n"
	     <<"while the phase space was being generated."<<std::endl;
    throw;
  }  
  // ---- FLUJO y promediado !!
  weight1 *=1.0;
  amplitude_hhc2ph amp(_M_proc,amplitude_hhc2ph::lo);
  // Following commented lines witten in phasespace routine
  // // EXTRA TERM Espacio de fases
  // double tmp;
  // double pTmax = sqrt(s)/2.;
  // double pT2 = pTmax*pTmax*pow((tmp=x[6]),3.0);
  // // pT2 = 100.;
  // double pTJac = 3.0*pTmax*pTmax*pow(tmp,2.0);
  // // pTJac = 1.;
  // p1[1] = lorentzvector<double>(0.0,sqrt(pT2),0.0,sqrt(pT2));
  double pTJac = 1.;
  ///////////  CONTRATERMINO NLO ///////////////////////////
  int  order=1;
  amp.CTgg_contrib(weight1*pTJac,p1,z1,z2,order);
  userfunc(p1,weightCTgg,amp,2,20);
  //////////////////////////////////////////////////////////
  delete _M_ps;
  delete _M_proc;
  return weightCTgg;
}


// --------------------------------------------------------------------------
//// REAL NLO ggINTg
// --------------------------------------------------------------------------
double Integral_class::Integrando_ggintg(const double x[],event_hhc2ph &pmed)
{  
  _Evento Evnt;
  // -----------------------------------------------------------------
  // -- Unused parameters  
  // -----------------------------------------------------------------
  unsigned int nu=2,nd=3;
  rng_emp rng( -238,8);
  process_hhc2ph *_M_proc = new hhc2ph1jet(rng, true, nu, nd, 0.1);
  // -----------------------------------------------------------------
  Evnt.Set_ps(Evnt);
  event_hhc2ph p2(1);
  double s=Evnt._M_s,q2min=Evnt.q2min,q2max=Evnt.q2max;
  basic_phasespace_hhc2ph * _M_ps = new basic_phasespace_hhc2ph(&rng, s);
  _M_ps -> important_sampling(2.0, 0.1);
  _M_ps -> phasespace_cuts(q2min, q2max);
  // -----------------------------------------------------------------
  double weightggINTg=0.0,weight2=0.0,weightmul=0.0;
  bool flagrmed=false;double FacNew=11.;  
  try {
    if(opts.withsmearing)
      {
        flagrmed = true;
	weightmul = _M_ps->operator()(pmed,flagrmed,FacNew,x);
	flagrmed = false;
	weight2 = _M_ps->operator()(p2,flagrmed,FacNew,x);
      }
    else
      {
	weight2 = _M_ps->operator()(p2,flagrmed,FacNew,x);
        pmed = p2;
      }
  } catch(numeric_error) {
    delete _M_ps;
    delete _M_proc;
    return 0.0;
  } catch(fp_exception) {
    delete _M_ps;
    delete _M_proc;
    return 0.0;
  } catch(...) {delete _M_ps;
    delete _M_proc;
    std::cerr<<"integral<..., 2, 0> : It was an unexpected error\n"
	     <<"while the phase space was being generated."<<std::endl;
    throw;
  }
  amplitude_hhc2ph amp(_M_proc,amplitude_hhc2ph::lo);  
  //----- real term NLO -----------------
  amp.ggINTg_contrib(weight2, p2);
  userfunc(p2,weightggINTg, amp,4.0,4.0);
  //--------------------------------------
  delete _M_ps;
  delete _M_proc;
  return weightggINTg;  
}

// ----------------------------------------------------------------------


/////////////////////////////////////////
//Interference LO QG and QQBAR channels//
/////////////////////////////////////////

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
//// Signal and interference LO------------ LEANDRO---------------------
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------

double Integral_class::Integrando_qqbarhyy(const double x[],event_hhc2ph &pmed)
{
  _Evento Evnt;
  // -----------------------------------------------------------------
  // -- Unused parameters  
  // -----------------------------------------------------------------
  unsigned int nu=2,nd=3;
  rng_emp rng(-238,8);
  process_hhc2ph *_M_proc = new hhc2ph1jet(rng, true, nu, nd, 0.1);
  // -----------------------------------------------------------------
  event_hhc2ph p1(1);
  event_hhc2ph psolo(1);
  Evnt.Set_ps(Evnt);
  // Evnt.q2min = phasespace::mmin*phasespace::mmin;
  // Evnt.q2max = phasespace::mmax*phasespace::mmax;
  // Evnt._M_s  = opts.sroot*opts.sroot;
  double s=Evnt._M_s,q2min=Evnt.q2min,q2max=Evnt.q2max;
  //  cout << "s : " << s << "  ;  " << " q2min : " << q2min << "  ;  " << " q2max : " << q2max << endl; 
  basic_phasespace_hhc2ph * _M_ps = new basic_phasespace_hhc2ph(&rng, s);
  _M_ps -> important_sampling(2.0, 0.1);
  _M_ps -> phasespace_cuts(q2min, q2max);
  // ----------------------------------------------------------------- 
  double weight1=0.0,wghtfIntLOqq=0.0;
  double weightmed=0.0;
  double Q2,z1,z2,smearing;
  smearing = opts.SigmaSmearing;
  int flag;
  bool flagrmed = false;
  // ----------------------------------------------------------------- 
  try {

    if (opts.withsmearing) 
      {
	weight1 = _M_ps->operator()(p1,Q2,smearing,flag,flagrmed,x);    
        // -- flag med activated
        flagrmed = true;
	weightmed = _M_ps->operator()(pmed,Q2,smearing,flag,flagrmed,x);
      }
    else
      {
	weight1 = _M_ps->operator()(p1,Q2,smearing,flag,flagrmed,x);   
        pmed = p1;
      }  
  } catch(numeric_error) {
    
    delete _M_ps;
    delete _M_proc;
    return 0.0;
  } catch(fp_exception) {
    
    
    delete _M_ps;
    delete _M_proc;
    return 0.0;
  } catch(...) {delete _M_ps;
    delete _M_proc;
    std::cerr<<"integral<..., 2, 0> : It was an unexpected error\n"
	     <<"while the phase space was being generated."<<std::endl;
    throw;
  }
  // ----------------------------------------------------------------- 
  // ---- FLUJO y promediado incluido en hhc2phjet.cc!!
  //weight1 *=1.0;
  // ----------------------------------------------------------------- 
  amplitude_hhc2ph amp(_M_proc,amplitude_hhc2ph::lo);
  // ----------------------------------------------------------------- 
  //////// CONTRIBUCION Interferencia total  qqbar+qg+gq!!!!///////
  
  amp.MqqbarHyy_contrib(weight1,p1);
  userfunc(p1,wghtfIntLOqq,amp,19,19);
  
  // cout << "pT from med : " << (pmed[-2]+pmed[-3]).perp() << endl;
  // cout << "pT from p1  : " << (p1[-2]+p1[-3]).perp() << endl;

  //////////////////////////////////////////////////////
 delete _M_ps;
 delete _M_proc;
 return wghtfIntLOqq;
}










