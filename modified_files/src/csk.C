#include "dyturbo.h"
#include "banner.h"
#include "settings.h"

#include "npff.h"
#include "blim.h"
#include "pdfevol.h"
#include "rad.h"
#include "cusp.h"
#include "pdf.h"
#include "alphas.h"
#include "scales.h"
#include "gaussrules.h"

#include "TError.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TAxis.h"
#include "TFile.h"
#include "TGraph.h"

#include <iostream>

using namespace std;

double csk(double *x, double *par)
{
  bool vfn = false;

  //b is the variable to plot along the x axis
  complex <double> b = x[0];
  
  //Set scales
  blim::set();
  pdfevol::scales(b);
  scales::set(opts.zmass);

  //Retrieve gk from the npff module
  double gk = -real(npff::gk(b));

  //Compute as(mubstar)
  double as;
  if (vfn)
    as = pdf::alphas(real(pdfevol::mubstar));
  else
    as = real(alphas::calc(pdfevol::mubstar, opts.order_sudak+1)) * M_PI;

  cout << "as analytic             " << gk << "    " << endl;
  //cout << "ren             " << scales::ren << "    " << endl;

  //The rapidity anomalous dimension is in as/4/pi normalisation
  double aso4pi = as / 4. / M_PI;

  //Compute the contribution from the rapidity anomalous dimension
  double KB = 0;
  if (opts.order_sudak >= 2)
    KB += pow(aso4pi,2) * rad::K1;
  if (opts.order_sudak >= 3)
    KB += pow(aso4pi,3) * rad::K2;
  if (opts.order_sudak >= 4)
    KB += pow(aso4pi,4) * rad::K3;

  //Not sure about additional factors to be applied to KB
  //KB = KB/2.;
  //KB = -KB;
  //KB = 2.*KB;

  //Compute the integrals  int_mubstar^mu=2 dmu'/mu' * as(mu')^k/pi
  int rule = 200;
  double q0 = real(pdfevol::mubstar);
  double q1 = 2;
  double cc = 0.5*(q1+q0);
  double mm = 0.5*(q1-q0);
  double as1 = 0;
  double as2 = 0;
  double as3 = 0;
  double as4 = 0;
  double as5 = 0;
  for (int i = 0; i < rule; i++)
    {
      double mu = cc+mm*gr::xxx[rule-1][i];
      //The Cusp anomalous dimension is in as/pi normalisation
      double asopi;
      if (vfn)
	asopi = pdf::alphas(real(pdfevol::mubstar))/M_PI;
      else
	asopi = real(alphas::calc(mu, opts.order_sudak+1));
      double fac = gr::www[rule-1][i]*mm;
      as1 += fac*asopi/mu;
      as2 += fac*pow(asopi,2)/mu;
      as3 += fac*pow(asopi,3)/mu;
      as4 += fac*pow(asopi,4)/mu;
      as5 += fac*pow(asopi,5)/mu;
      //cout << mu << "  " << fac << "  " << pdf::alphas(mu,opts.cmw) << endl;
    }

  //Compute the contribution from the Cusp anomalous dimension,
  double KA = 0;
  if (opts.order_sudak >= 0)
    KA += as1 * 2.*cusp::Gamcusp1;
  if (opts.order_sudak >= 1)
    KA += as2 * 2.*cusp::Gamcusp2;
  if (opts.order_sudak >= 2)
    KA += as3 * 2.*cusp::Gamcusp3;
  if (opts.order_sudak >= 3)
    KA += as4 * 2.*cusp::Gamcusp4;
  if (opts.order_sudak >= 4)
    KA += as5 * 2.*cusp::Gamcusp5;
  
  //Sum up all components
  double K = KB - KA;
  double f = K + gk;
  //cout << setw(15) << real(b) << setw(15) << real(pdfevol::mubstar) << setw(15) << as << setw(15) << -0.5*(-KA) << setw(15) << -0.5*KB << setw(15) << -0.5*K << setw(15) << -0.5*gk << setw(15) << -0.5*f << endl;
  return -f/2.; //Because we plots -1/2 K(b,mu=2)
}

int main(int argc, char * argv[])
{
  banner();  

  //Init constants
  DYTurbo::init_const();

  opts.parse_options(argc,argv);                      //parse options from command line (and parse input file)

  DYTurbo::init_params();                                      //Initialisation which depends on input file settings
  
  if (opts.verbose) opts.dumpAll();
  opts.silent = true;
  //  DYTurbo::Init(argc,argv);              //Init, read config file

  cout << "Compute Collins-Soper kernel" << endl;

  // Set the verbosity to show only errors and above
  gErrorIgnoreLevel = kError;

  string csk_filename="results";
  if (opts.output_filename!="") csk_filename = opts.output_filename;
  csk_filename  = csk_filename + "-csk";
  string pdf_filename  = csk_filename + ".pdf";
  
  TCanvas *c = new TCanvas();
  c->cd();

  int npts = 100;
  TF1 *fcsk = new TF1("csk", csk, 5./npts, 5);
  fcsk ->SetNpx(npts-1);
  fcsk->GetXaxis()->SetTitle("b");
  fcsk->GetYaxis()->SetTitle("-#frac{1}{2}K(b,#mu=2 GeV)");
  fcsk->GetYaxis()->SetRangeUser(-0.2,0.8);
  fcsk->SetTitle("");
  fcsk->Draw();
  c->Print(pdf_filename.c_str());
  cout << "CSK kernel plotted in " << pdf_filename << endl;

  double xmin = fcsk->GetXmin();
  double xmax = fcsk->GetXmax();

  double x[npts];
  double y[npts];

  for (int i = 0; i < npts; ++i)
    {
        x[i] = xmin + (xmax - xmin) * i / (npts - 1);
        y[i] = fcsk->Eval(x[i]);
    }

  TGraph *g = new TGraph(npts, x, y);
  g->SetName("csk");
  g->SetTitle("");
  
  string root_filename  = csk_filename + ".root";
  TFile *fout = new TFile(root_filename.c_str(),"RECREATE");

  fcsk->Write();
  //  g->Write();

  fout->Close();
  cout << "Results saved in " << root_filename << endl;

  //cout << setw(15) << "b" << setw(15) << "mubstar" << setw(15) << "as" << setw(15) << "-0.5*KA" << setw(15) << "-0.5*KB" << setw(15) << "-0.5*Kpert" << setw(15) << "-0.5*gk" << setw(15) << "-0.5*Ktot" << endl;
  //double par[1];
  //double xx[1];
  //for (double b = 0.1; b <= 5.0; b += 0.1)
  //  {
  //    xx[0] = b;
  //    csk(xx, par);
  //  }
    
}
