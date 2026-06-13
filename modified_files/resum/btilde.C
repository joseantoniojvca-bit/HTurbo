#include "btilde.h"
#include "settings.h"
#include "resconst.h"
#include "phasespace.h"
#include "gaussrules.h"
#include "besselint.h"
//#include "TMath.h"

#include <iostream>

using namespace std;

const int btilde::p = 4;

double x1 = 0.3;
double x2 = 0.6;
double x3 = 0.9;

int fact(int n)
{
  int p = 1;
  for (int i = 1; i <= n; i++)
    p *= i;
  return p;
}

double int1mj0ox(double bq)
{
  //Compute int_0^bq {dx/x * (1-J0(x))}
  int qrule = 100;
  double res = 0.;

  //quadrule loop
  double cc = 0.5*(bq);
  double mm = 0.5*(bq);
  for (int k = 0; k < qrule; k++)
    {
      double x = cc+mm*gr::xxx[qrule-1][k];
      double fac = gr::www[qrule-1][k]*mm;
      res += (1-fort_besj0_(x))/x * fac;
    }
  
  return res;
}

complex <double> btilde::mubtl(complex <double> b, double Q)
{
  complex <double> mub = resconst::b0/b;
  complex <double> mubtilde;
  if (!opts.modlog)
    mubtilde = mub;
  else if (opts.p == 1)
    mubtilde = mub * Q / sqrt((pow(mub,2) + pow(Q,2)));
  else
    mubtilde = mub * Q / pow((pow(mub,2*opts.p) + pow(Q,2*opts.p)),0.5/opts.p);

  //cout << b << "  " << mub << "  " << mubtilde << endl;
  return mubtilde;
}

complex <double> btilde::calc(complex <double> b, double Q)
{
  complex <double> mubtilde = mubtl(b,Q);
  return resconst::b0/mubtilde;
}

complex <double> btilde::blog(complex <double> b, double Q)
{
  if (opts.modlog == 0)
    return 2.*log(Q*b/resconst::b0);
  
  complex <double> blog = 2.*log(Q/mubtl(b,Q));
  //if (opts.modlog && real(blog) < 0.)
  //blog = max(0.,real(blog));

  //cout << " Q " << Q << "  mubtl " << mubtl(b,Q) << " blog " << blog << "  " << " qt " << phasespace::qt << endl;
  
  return blog;
}
