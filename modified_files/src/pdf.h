#ifndef pdf_h
#define pdf_h

#include "parton.h"
#include "config.h"

#ifdef USENEOPDF
#include <NeoPDF.hpp>
#endif
#include <LHAPDF/LHAPDF.h>

namespace pdf
{
  extern bool secondbeam;
  extern LHAPDF::PDF* lhapdf;
#ifdef USENEOPDF
  extern neopdf::NeoPDF* neopdf
#endif

  extern void init();
  extern void setalphas();
  extern void setg();

  extern void (*xfxq)(const double &x, const double &Q, double *fPDF);
  extern void (*xfxqphot)(const double &x, const double &Q, double &fotfx);
  extern void lhaxfxq(const double &x, const double &Q, double *fPDF);
  extern void lhaxfxqphot_beam1(const double &x, const double &Q, double &fotfx);
  extern void lhaxfxqphot_beam2(const double &x, const double &Q, double &fotfx);
  extern void lhaxfxqphot_empty(const double &x, const double &Q, double &fotfx);

  extern void fdist(int beam, double x, double mu, double fx[2*MAXNF+1]);
  extern void fdistphot(int beam, double x, double mu, double &fotfx);
  
  extern double (*alphas)(const double &Q); //external alphas routine (usually linked to lhapdf alphas)
  extern double lhaalphas(const double &Q);    //wrapper to lhapdf alphas
  extern double rgktalphas(double q);          //wrapper to internal Runge-Kutta alphas
  
  extern int order;
  extern double xmin;
  extern double qmin;
  extern double mc;
  extern double mb;
  extern double mt;
  extern double g;
  extern double g1;
  extern double g1c;
  extern double g1m;
  extern double g1f;
  extern double g2;
  extern double e;
  extern double g0;
  extern double q;
  extern double blim;
  extern double lambda;
  extern bool photon;
}

extern "C" {
  void fdist_(int& ih, double& x, double& xmu, double fx[2*MAXNF+1]);
  void pdf_fdist_(int& beam, double& x, double& mu, double fx[2*MAXNF+1]);
  void fdistphot_(int& beam, double& x, double& xmu, double &fotfx);
  void dysetpdf_(int& member);
  void setmellinpdf_(int& member);
  double dyalphas_mcfm_(double &q, double &amz, int &nloop);
  //  double dyalphas_lhapdf_(double &q); //no need for this, as there is a C++ LHAPDF function
}



#endif
