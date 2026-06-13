#include "pdf.h"
#include "dyres_interface.h"
#include "interface.h"
#include "scales.h"
#include "settings.h"
#include "coupling.h"
#include "alphas.h"
#include "isnan.h"
#include "vjint.h"
#include "mcfm_interface.h"
#include "as3vjet_interface.h"
#include "resconst.h"
#include "constants.h"
#include "parton.h"
#include "config.h"

#ifdef USENEOPDF
#include <NeoPDF.hpp>
#endif
#include <LHAPDF/LHAPDF.h>
#include <LHAPDF/LHAGlue.h>
#include <math.h>

using namespace constants;

bool pdf::secondbeam = false;
LHAPDF::PDF* pdf::lhapdf = 0;

#ifdef USENEOPDF
neopdf::NeoPDF* pdf::neopdf = 0;
#endif

int pdf::order;
double pdf::xmin;
double pdf::qmin = 0.;
double pdf::mc;
double pdf::mb;
double pdf::mt;
double pdf::g;
double pdf::g1;
double pdf::g1c;
double pdf::g1m;
double pdf::g1f;
double pdf::g2;
double pdf::e;
double pdf::g0;
double pdf::q;
double pdf::blim;
double pdf::lambda;
bool pdf::photon = false;

std::vector<int32_t> pids = { -6, -5, -4, -3, -2, -1, 21, 1, 2, 3, 4, 5, 6 };

void (*pdf::xfxq)(const double &x, const double &Q, double *fPDF) = 0;
void (*pdf::xfxqphot)(const double &x, const double &Q, double &fotfx) = 0;
double (*pdf::alphas)(const double &Q) = 0;

void pdf::lhaxfxq(const double &x, const double &Q, double *r)
{
  vector <double> fPDF; fPDF.resize(13);
#ifdef USENEOPDF
  //fPDF[6-5] = neopdf_beam1->xfxQ2(-5,x,Q*Q);
  //fPDF[6-4] = neopdf_beam1->xfxQ2(-4,x,Q*Q);
  //fPDF[6-3] = neopdf_beam1->xfxQ2(-3,x,Q*Q);
  //fPDF[6-2] = neopdf_beam1->xfxQ2(-2,x,Q*Q);
  //fPDF[6-1] = neopdf_beam1->xfxQ2(-1,x,Q*Q);
  //fPDF[6+0] = neopdf_beam1->xfxQ2(21,x,Q*Q);
  //fPDF[6+1] = neopdf_beam1->xfxQ2(1 ,x,Q*Q);
  //fPDF[6+2] = neopdf_beam1->xfxQ2(2 ,x,Q*Q);
  //fPDF[6+3] = neopdf_beam1->xfxQ2(3 ,x,Q*Q);
  //fPDF[6+4] = neopdf_beam1->xfxQ2(4 ,x,Q*Q);
  //fPDF[6+5] = neopdf_beam1->xfxQ2(5 ,x,Q*Q);
  std::vector<double> kins = { x, Q*Q };
  //std::vector<double> neo_xfs;
  fPDF = neopdf->xfxQ2_pids(pids, kins);
#else
  lhapdf->xfxQ(x,Q,fPDF);
#endif
  copy(fPDF.begin(),fPDF.end(),r);
}

void pdf::lhaxfxqphot_beam1(const double &x, const double &Q, double &fotfx)
{
  //!!! Why here the PDF is divided by x???
#ifdef USENEOPDF
  fotfx = neopdf->xfxQ2(22,x,Q*Q)/x;
#else
  fotfx = pdf::lhapdf->xfxQ(22,x,Q)/x;
  //fotfx = lhapdf_beam1->xfxQ(22,x,Q)/x;
#endif
}

void pdf::lhaxfxqphot_beam2(const double &x, const double &Q, double &fotfx)
{
  //!!! Why here the PDF is divided by x???
#ifdef USENEOPDF
  fotfx = neopdf->xfxQ2(22,x,Q*Q)/x;
#else
  fotfx = lhapdf->xfxQ(22,x,Q)/x;
#endif
}

void pdf::lhaxfxqphot_empty(const double &x, const double &Q, double &fotfx)
{
  fotfx = 0.;
}

double pdf::lhaalphas(const double &Q)
{
  //Old LHAPDF interface
  //return LHAPDF::alphasPDF(Q);

  //New LHAPDF interface
  return lhapdf->alphasQ(Q);
}

double pdf::rgktalphas(double q)
{
  //Interface to the Runge-Kutta alphas solution, with order matched to opts.order
  //This is always the MSbar alphas

  //Take alphas(mZ) from the external (lhapdf) alphas
  double q0 = coupling::zmass;
  double as0 = lhapdf->alphasQ(q0);
  
  //Evolve from mZ to mu with Runge-Kutta (MSbar evolution)
  alphas::rgkt(q,q0,as0);
  switch (opts.order)
    {
    case 0: return real(alphas::asLO);     break;
    case 1: return real(alphas::asNLO);    break;
    case 2: return real(alphas::asNNLO);   break;
    case 3: return real(alphas::asNNNLO);  break;
    case 4: return real(alphas::asNNNNLO); break;
    default: return 0.;
    }
}

void pdf::init()
{
  if (!opts.externalpdf)
    {
      //printf(" ==Initialize PDF set from LHAPDF==\n\n");
      //printf("\n");
      LHAPDF::Info& cfg = LHAPDF::getConfig();

      cfg.set_entry("Verbosity"             , 0              );
      cfg.set_entry("Interpolator"          , "logcubic"     );
      cfg.set_entry("Extrapolator"          , "continuation" );
      cfg.set_entry("ForcePositive"         , 0              );
      cfg.set_entry("AlphaS_Type"           , "analytic"     );
      cfg.set_entry("MZ"                    , 91.1876        );
      cfg.set_entry("MUp"                   , 0.002          );
      cfg.set_entry("MDown"                 , 0.005          );
      cfg.set_entry("MStrange"              , 0.10           );
      cfg.set_entry("MCharm"                , 1.29           );
      cfg.set_entry("MBottom"               , 4.19           );
      cfg.set_entry("MTop"                  , 172.9          );
      cfg.set_entry("Pythia6LambdaV5Compat" , true           );

      //Old interface
      LHAPDF::initPDFSet(opts.LHAPDFset); //LHAPDF::initPDFSetByName(opts.LHAPDFset);
      LHAPDF::initPDF(opts.LHAPDFmember);

      //New interface
      lhapdf = LHAPDF::mkPDF(opts.LHAPDFset, opts.LHAPDFmember);
#ifdef USENEOPDF
      neopdf = new neopdf::NeoPDF(opts.LHAPDFset.c_str(), opts.LHAPDFmember);
#endif
      if (opts.LHAPDFset_beam2 != "none")
	{
#ifdef USENEOPDF
	  neopdf = new neopdf::NeoPDF(opts.LHAPDFset_beam2.c_str(), opts.LHAPDFmember_beam2);
#else
	  lhapdf = LHAPDF::mkPDF(opts.LHAPDFset_beam2, opts.LHAPDFmember_beam2);
#endif
	  secondbeam = true;
	}
      else
#ifdef USENEOPDF
	neopdf = new neopdf::NeoPDF(opts.LHAPDFset.c_str(), opts.LHAPDFmember);
#else
	lhapdf = LHAPDF::mkPDF(opts.LHAPDFset, opts.LHAPDFmember);
#endif
      
      xfxq = lhaxfxq;

      //if (LHAPDF::hasPhoton())
      if (lhapdf->hasFlavor(22))
      {
        photon = true;
        xfxqphot = lhaxfxqphot_beam1;
      }
      else
      {
        photon = false;
        xfxqphot = lhaxfxqphot_empty;
      }
      
      //if (opts.lhaalphas)
      alphas = lhaalphas;
      //else
      //alphas = rgktalphas

      if (opts.PDFerrors && LHAPDF::numberPDF() > 1)
      {
        opts.totpdf = LHAPDF::numberPDF()+1;
        pdferropts_.pdferr_ = true;
        pdferropts_.totpdf_ = LHAPDF::numberPDF()+1;
      }
      else
      {
        opts.totpdf = 1;
        pdferropts_.pdferr_ = false;
        pdferropts_.totpdf_ = 1;
      }

      order = LHAPDF::getOrderPDF(); //order of evolution

      LHAPDF::PDFInfo info(opts.LHAPDFset, opts.LHAPDFmember);
      qmin = info.get_entry_as<double>("QMin", -1);
      xmin = info.get_entry_as<double>("XMin", -1);

      mc = LHAPDF::getThreshold(4);
      mb = LHAPDF::getThreshold(5);
      mt = LHAPDF::getThreshold(6);
    }
    else
    {
      opts.totpdf = 1;
      pdferropts_.pdferr_ = false;
      pdferropts_.totpdf_ = 1;
    }

    // initialization of alphas
    setalphas();

    //read g from the PDF
    setg();

    //take the cmass and b mass from the PDF
    //      cmass=dsqrt(mcsq)
    //      bmass=dsqrt(mbsq)
}
/*
//run alphas
double pdf::alphas(double q, bool cmw)
{
  // **** obsolete running used in DYRes **** /
  //if (opts_.approxpdf_ == 1)
  //  {
  //    double asmz = extalphas(coupling::zmass);
  //    int nloop = 3;
  //    as = dyalphas_mcfm_(scales::ren,asmz,nloop);
  //  }
  // **************************************** /

  //MSbar evolution
  double as;
  if (opts.alphaslha)
    as = extalphas(q);
  else
    as = rgktalphas(q);

  //Conversion to CMW
  if (cmw)
    switch (opts.order)
      {
      case 0: as = as; break;
      case 1: as = as*(1 + as/pi*resconst::K1 ); break;
      case 2: as = as*(1 + as/pi*resconst::K1 + pow(as/pi,2)*resconst::K2 ); break;
      case 3: as = as*(1 + as/pi*resconst::K1 + pow(as/pi,2)*resconst::K2 + pow(as/pi,3)*resconst::K3 ); break;
      case 4: as = as*(1 + as/pi*resconst::K1 + pow(as/pi,2)*resconst::K2 + pow(as/pi,3)*resconst::K3 + pow(as/pi,4)*resconst::K4 ); break;
      }

  //cout << "alphas(q) at q = " << q << "  " << as << endl;
  //cout << resconst::K1 << endl;
  
  return as;
}
*/

//set value of alphas at the renormalization scale in the fortran common blocks
void pdf::setalphas()
{
  //couple_.amz_ = alphas(dymasses_.zmass_);
  double asmz = alphas(coupling::zmass);

  //run alphas
  double as;
  if (opts_.approxpdf_ == 1)
  {
    int nloop = 3;
    as = dyalphas_mcfm_(scales::ren,couple_.amz_,nloop);
  }
  else
    //as = dyalphas_lhapdf_(scales::ren);
    //as = LHAPDF::alphasPDF(scales::ren); //Old LHAPDF interface
    //as = lhapdf->alphasQ(scales::ren);   //New LHAPDF interface
    if (opts.alphaslha)
      as = lhapdf->alphasQ(scales::ren);		  //Allows external alphas
      else
        as = rgktalphas(scales::ren);

  //MCFM
  qcdcouple_.as_ = as;
  qcdcouple_.ason2pi_ = as/(2*M_PI);
  qcdcouple_.ason4pi_ = as/(4*M_PI);
  qcdcouple_.gsq_= 4*M_PI*as;
  //cout <<" as " << as << endl;

  //as3vjet
  as3vj_qcdcouple_.as_ = as;
  as3vj_qcdcouple_.ason2pi_ = as/(2*M_PI);
  as3vj_qcdcouple_.ason4pi_ = as/(4*M_PI);
  as3vj_qcdcouple_.gsq_= 4*M_PI*as;
  
  //vjet analytical calculation
  asnew_.as_ = as/M_PI;
  asp_.asp_ = as;
}

//set the value of the g-parameter of the non perturbative form factor
void pdf::setg()
{
  if (!opts.externalpdf)
    {
      LHAPDF::PDFInfo info(opts.LHAPDFset, opts.LHAPDFmember);
      g = info.get_entry_as<double>("g", -1);
      g1 = info.get_entry_as<double>("g1", -1);
      g1c = info.get_entry_as<double>("g1c", -1);
      g1m = info.get_entry_as<double>("g1m", -1);
      g1f = info.get_entry_as<double>("g1f", -1);
      g2 = info.get_entry_as<double>("g2", -1);
      e  = info.get_entry_as<double>("e", -1);
      g0  = info.get_entry_as<double>("g0", -1);
      q  = info.get_entry_as<double>("q", -1);
      blim  = info.get_entry_as<double>("blim", -1);
      lambda  = info.get_entry_as<double>("lambda", -1);
    }
  if (g >= 0)
    {
      cout << "g form factor: input from PDF member: " << g << endl;
      opts.g1 = g;
      g_param_.g_param_ = g;
      np_.g_ = g;
    }
  if (g1c >= 0)
    opts.g1c = g1c;
  if (g1m >= 0)
    opts.g1m = g1m;
  if (g1f >= 0)
    opts.g1f = g1f;
  if (g2 >= 0)
    opts.g2 = g2;
  if (e >= 0)
    opts.e = e;
  if (g0 >= 0)
    opts.g0 = g0;
  if (q > -1)
    opts.q = q;
  if (blim >= 0)
    {
      opts.blim = blim;
      opts.blim_pdf     = blim;
      opts.blim_sudakov = blim;
      opts.blim_expc    = blim;
    }
  if (lambda > -1)
    opts.lambda = lambda;
  
  if (g1 >= 0 || g2 >= 0)
    cout << "g form factor: input from PDF member: g1 = " << g1 << "; g2 = " << g2 << endl;
  if (g1c >= 0 || g1m >= 0 || g1f >= 0)
    cout << "g y-depedent form factor: input from PDF member: g1c = " << g1c << "; g1m = " << g1m << "; g1f = " << g1f << endl;
  if (e >= 0)
    cout << "e form factor: input from PDF member: e = " << e << endl;
  if (g0 >= 0)
    cout << "g0 form factor: input from PDF member: g0 = " << g0 << endl;
  if (q > -1)
    cout << "q form factor: input from PDF member: q = " << q << endl;
  if (blim >= 0)
    cout << "blim: input from PDF member: blim = " << blim << endl;
  if (lambda > -1)
    cout << "lambda form factor: input from PDF member: lambda = " << lambda << endl;
}


void dysetpdf_(int& member)
{
  if (opts.externalpdf) return;
  
  if (member == 0)
    {
      if (opts.PDFerrors && opts.totpdf > 1)
	LHAPDF::initPDF(0);
      else
	LHAPDF::initPDF(opts.LHAPDFmember);
    }
  else
    LHAPDF::initPDF(member);
  
  pdf::setalphas();
  //  setg(); //set g separately when setting a different PDF in the resummed part
}


void setmellinpdf_(int &member){
    // if member is still same than dont do anything
    //if (lastMember==member) return;
    //if (v_mellinpdf.size()<member) v_mellinpdf.
    // test current flag
    //if (v_mellinpdf[lastMember].isInitialized)
    // weights
    // moms
    // set init flag
    //lastMember=member;
}


void fdist_(int& ih, double& x, double& xmu, double fx[2*MAXNF+1])
{
  //set to zero if x out of range
  //check for nans
  if (isnan_ofast(x))
    {
      cout << "Bjorken-x is nan in fdist_" << endl;
      for (int i = -MAXNF; i <= MAXNF; i++)
	fx[MAXNF+i]=0.;
      return;
    }
  if (x > 1. || x <= 0.)
    {
      for (int i = -MAXNF; i <= MAXNF; i++)
	fx[MAXNF+i]=0.;
      return;
    }

  //Old LHAPDF interface
  //double fPDF[13];
  //LHAPDF::xfx(x,xmu,fPDF);

  //New LHAPDF interface
  //vector <double> fPDF; fPDF.resize(13);
  //pdf::lhapdf->xfxQ(x,xmu,fPDF);

  //Allows external PDF
  double fPDF[13];
  pdf::xfxq(x,xmu,fPDF);

  //cout << "fdist " << setprecision(16) << fPDF[6] << endl;
  //vector<int> pids = pdf::lhapdf->flavors();
  //for (int i = 0; i < pids.size(); i++)
  //cout << pids[i] << endl;

  //cout << " fdist x:" << x << " q: " << xmu << " - ";
  //for (int i = -MAXNF; i <= MAXNF; i++)
  //  cout << fPDF[6+i]/x << "  ";
  //cout << endl;
  
  //xmu *= xmu;
  //fPDF[6-5] = pdf::lhapdf->xfxQ(-5,x,xmu);
  //fPDF[6-4] = pdf::lhapdf->xfxQ(-4,x,xmu);
  //fPDF[6-3] = pdf::lhapdf->xfxQ(-3,x,xmu);
  //fPDF[6-2] = pdf::lhapdf->xfxQ(-2,x,xmu);
  //fPDF[6-1] = pdf::lhapdf->xfxQ(-1,x,xmu);
  //fPDF[6+0] = pdf::lhapdf->xfxQ(21,x,xmu);
  //fPDF[6+1] = pdf::lhapdf->xfxQ(1 ,x,xmu);
  //fPDF[6+2] = pdf::lhapdf->xfxQ(2 ,x,xmu);
  //fPDF[6+3] = pdf::lhapdf->xfxQ(3 ,x,xmu);
  //fPDF[6+4] = pdf::lhapdf->xfxQ(4 ,x,xmu);
  //fPDF[6+5] = pdf::lhapdf->xfxQ(5 ,x,xmu);

  //  for (int i = -MAXNF; i <= MAXNF; i++)
  //    cout << fPDF[6+i]/x << "  ";
  //  cout << endl;
  

  if (ih == 1) //proton
    for (int i = -MAXNF; i <= MAXNF; i++)
      fx[MAXNF+i]=fPDF[6+i]/x;
  else if (ih == -1) //antiproton
    for (int i = -MAXNF; i <= MAXNF; i++)
      fx[MAXNF+i]=fPDF[6-i]/x;

  //switch off flavours
  //fx[MAXNF-5]=0.; //bbar
  //fx[MAXNF-4]=0.; //cbar
  //fx[MAXNF-3]=0.; //sbar
  //fx[MAXNF-2]=0.;
  //fx[MAXNF-1]=0.;
  //fx[MAXNF+0]=0.; //gluon
  //fx[MAXNF+1]=0.;
  //fx[MAXNF+2]=0.;
  //fx[MAXNF+3]=0.; //s
  //fx[MAXNF+4]=0.; //c
  //fx[MAXNF+5]=0.; //b

  //impose positivity
  //fx[MAXNF-5]=max(0.,fx[MAXNF-5]);
  //fx[MAXNF-4]=max(0.,fx[MAXNF-4]);
  //fx[MAXNF-3]=max(0.,fx[MAXNF-3]);
  //fx[MAXNF-2]=max(0.,fx[MAXNF-2]);
  //fx[MAXNF-1]=max(0.,fx[MAXNF-1]);
  //fx[MAXNF+0]=max(0.,fx[MAXNF+0]);
  //fx[MAXNF+1]=max(0.,fx[MAXNF+1]);
  //fx[MAXNF+2]=max(0.,fx[MAXNF+2]);
  //fx[MAXNF+3]=max(0.,fx[MAXNF+3]);
  //fx[MAXNF+4]=max(0.,fx[MAXNF+4]);
  //fx[MAXNF+5]=max(0.,fx[MAXNF+5]);

  //make u = d and ubar = dbar
  //fx[MAXNF-1]=fx[MAXNF-2];
  //fx[MAXNF+1]=fx[MAXNF+2];

  /*
  if (x < 90./13000.*exp(-1.))
    {
      fx[MAXNF-5]=0.;
      fx[MAXNF-4]=0.;
      fx[MAXNF-3]=0.;
      fx[MAXNF-2]=0.;
      fx[MAXNF-1]=0.;
      fx[MAXNF+0]=0.;
      fx[MAXNF+1]=0.;
      fx[MAXNF+2]=0.;
      fx[MAXNF+3]=0.;
      fx[MAXNF+4]=0.;
      fx[MAXNF+5]=0.;
    }
  */

  //flat PDFs
  //fx[MAXNF-5]=1.; //bbar
  //fx[MAXNF-4]=1.; //cbar
  //fx[MAXNF-3]=1.; //sbar
  //fx[MAXNF-2]=1.;
  //fx[MAXNF-1]=1.;
  //fx[MAXNF+0]=1.; //gluon
  //fx[MAXNF+1]=1.;
  //fx[MAXNF+2]=1.;
  //fx[MAXNF+3]=1.; //s
  //fx[MAXNF+4]=1.; //c
  //fx[MAXNF+5]=1.; //b
}

// --> instead of overloading the fdist function, add a new function which retrieves only the photon PDF
void fdistphot_(int& beam, double& x, double& xmu, double &fotfx)
{
  pdf::fdistphot(beam, x, xmu, fotfx);
}

void pdf::fdistphot(int beam, double x, double mu, double& fotfx)
{
  pdf::xfxqphot(x,mu,fotfx);
  if (opts.nuclearpdf)
    {
      double a = 0.;
      if (beam == 1)
	a = opts.A1;
      if (beam == 2)
	a = opts.A2;
      fotfx *= a;
    }
}

void pdf_fdist_(int& beam, double& x, double& mu, double fx[2*MAXNF+1])
{
  pdf::fdist(beam, x, mu, fx);
}

//nuclear PDFs
void pdf::fdist(int beam, double x, double mu, double fx[2*MAXNF+1])
{
  //set to zero if x out of range
  if (x > 1. || x <= 0.)
    {
      for (int i = -MAXNF; i <= MAXNF; i++)
	fx[MAXNF+i]=0.;
      return;
    }

  //Allows external PDF
  double fPDF[14] = {0.};
  if (beam == 2 && secondbeam)
    lhaxfxq(x,mu,fPDF);
  else
    xfxq(x,mu,fPDF);

  if (!opts.nuclearpdf)
    {
      //proton, antiproton or pions
      int ih;
      if (beam == 1)
	ih = opts.ih1;
      if (beam == 2)
	ih = opts.ih2;

      if (ih == 1 || ih == -211 || fabs(ih) > 2212) //proton or negative pion (xFitter or Jam)
	for (int i = -MAXNF; i <= MAXNF; i++)
	  fx[MAXNF+i]=fPDF[6+i]/x;
      else if (ih == -1 || ih == 211) //antiproton
	for (int i = -MAXNF; i <= MAXNF; i++)
	  fx[MAXNF+i]=fPDF[6-i]/x;
      else if (ih == -221)// pi- a la GRV
	{
	  for (int i = -MAXNF; i <= MAXNF; i++)
	    fx[MAXNF+i]=fPDF[6+i]/x;
	  fx[MAXNF-2]=fx[MAXNF+2];//ubar is becoming valence
	  fx[MAXNF+2]=fx[MAXNF-1];//u is becoming sea contribution
	  //do I have to add sea compontents to valence?
	}
      else if (ih == 221)//pi+ a la GRV
	{
	  for (int i = -MAXNF; i <= MAXNF; i++)
	    fx[MAXNF+i]=fPDF[6+i]/x;
	  fx[MAXNF-1]=fx[MAXNF+1];//dbar is becoming valence
	  fx[MAXNF+1]=fx[MAXNF-2];//d is becoming sea contribution
	  //do I have to add sea components to valence?
	}
    }      
  else
    {
      //number of protons
      double p;
      if (beam == 1)
	p = opts.Z1;
      if (beam == 2)
	p = opts.Z2;

      //number of neutrons
      double n = 0;
      if (beam == 1)
	n = opts.A1-opts.Z1;
      if (beam == 2)
	n = opts.A2-opts.Z2;
  
      double fxp[2*MAXNF+1];
      double fxn[2*MAXNF+1];
      for (int i = -MAXNF; i <= MAXNF; i++)
	{
	  fxp[MAXNF+i]=fPDF[6+i]/x;
	  fxn[MAXNF+i]=fPDF[6+i]/x;
	}
      fxn[parton::U]  = fxp[parton::D];
      fxn[parton::D]  = fxp[parton::U];
      fxn[parton::Ub] = fxp[parton::Db];
      fxn[parton::Db] = fxp[parton::Ub];

      for (int i = -MAXNF; i <= MAXNF; i++)
	fx[MAXNF+i]=p*fxp[MAXNF+i]+n*fxn[MAXNF+i];
      //      for (int i = -MAXNF; i <= MAXNF; i++)
      //	cout << fx[MAXNF+i] << endl;
    }

  //switch off flavours
  //fx[parton::Bb]=0.; //bbar
  //fx[parton::Cb]=0.; //cbar
  //fx[parton::Sb]=0.; //sbar
  //fx[parton::Ub]=0.; //ubar
  //fx[parton::Db]=0.; //dbar
  //fx[parton::G ]=0.; //gluon
  //fx[parton::D ]=0.; //d
  //fx[parton::U ]=0.; //u
  //fx[parton::S ]=0.; //s
  //fx[parton::C ]=0.; //c
  //fx[parton::B ]=0.; //b

  //flat PDFs
  //fx[parton::Bb]=1.; //bbar
  //fx[parton::Cb]=1.; //cbar
  //fx[parton::Sb]=1.; //sbar
  //fx[parton::Ub]=1.;
  //fx[parton::Db]=1.;
  //fx[parton::G ]=1.; //gluon
  //fx[parton::D ]=1.;
  //fx[parton::U ]=1.;
  //fx[parton::S ]=1.; //s
  //fx[parton::C ]=1.; //c
  //fx[parton::B ]=1.; //b
}
