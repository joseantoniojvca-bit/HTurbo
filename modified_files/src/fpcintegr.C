#include "fpcintegr.h"

#include "dyres_interface.h"
#include "interface.h"
#include "omegaintegr.h"
#include "phasespace.h"
#include "settings.h"
#include "cubacall.h"
#include "isnan.h"
#include "qtint.h"
#include "ctint.h"
#include "ctintmc.h"

#include "Kinematics.h"

#include <iostream>
#include <iomanip>
//!!! also FPC should have the splitting of the as3 contribution, because the qtcut is different for the as3 term !!!
int fpcintegrand2d_cubature_v(unsigned ndim, long unsigned npts, const double x[], void *data, unsigned ncomp, double f[])
{
  tell_to_grid_we_are_alive(1000);
  //  cout << "parallel " << npts << endl;
#pragma omp parallel for num_threads(opts.cubacores) copyin(a_param_,scale_,facscale_,qcdcouple_)
  for (unsigned i = 0; i < npts; i++)
    {
      // evaluate the integrand for npts points
      double xi[ndim];
      double fi[ncomp];
      for (unsigned j = 0; j < ndim; j++)
	xi[j] = x[i*ndim + j];

      fpcintegrand2d(ndim, xi, ncomp, fi);
      
      for (unsigned k = 0; k < ncomp; ++k)
	f[i*ncomp + k] = fi[k];
    }
  return 0;
}

int fpcintegrand2d_cubature(unsigned ndim, const double x[], void *data, unsigned ncomp, double f[])
{
  tell_to_grid_we_are_alive(1000);
  fpcintegrand2d(ndim, x, ncomp, f);
  return 0;
}

integrand_t fpcintegrand2d(const int &ndim, const double x[], const int &ncomp, double f[])
//Calculates the fpc integrand as a function of m, y
//The integration in qt is factorised in LL1, LL2, LL3 and LL4 large logs
//dOmega integration is factorised in the costh moments
//The integration in alpha and beta is performed inside countdy
{
  tell_to_grid_we_are_alive(1000);
  clock_t begin_time, end_time;

  begin_time = clock();

  if (!opts.fixedorder || phasespace::qtmin > 0)
    {
      f[0]=0.;
      return 0;
    }
  
//  if (opts.PDFerrors)
//    for (int i = 1; i < opts.totpdf; i++)
//      f[i] = 0.;

  //Jacobian of the change of variables from the unitary square x[2] to the m, y boundaries
  double jac = 1.;
  bool status = true;

  double r2[2] = {x[0], x[1]};
  status = phasespace::gen_my(r2, jac, !opts.fixedorder, !opts.fixedorder);  //qtcut = !opts.fixedorder, qtswitching = !opts.fixedorder
  if (!status)
    {
      f[0] = 0.;
      return 0;
    }

  //integrate from 0 to qtcut
  double qtcut = max(opts.qtcut,opts.xqtcut*phasespace::m);
  double qtmn, qtmx;
  qtmn = 0.;
  qtmx = qtcut;
  phasespace::set_qt(0.); //Should fill bins with qt=0
  phasespace::set_phiV(0.);
  double m = phasespace::m;
  double y = phasespace::y;

  clock_t qtbt, qtet;
  qtbt = clock();
  qtint::calc(m,qtmn,qtmx,2);
  qtet = clock();
  
  //generate boson 4-momentum, with m, qt, y and phi=0 (not actually needed)
  omegaintegr::genV4p();

  //In this point of phase space (m, qt, y) the costh integration is performed by 
  //calculating the 0, 1 and 2 moments of costh
  //that are the integrals int(dcosth dphi1 dphi2), int(costh dcosth dphi1 dphi2), and int(costh^2 dcosth dphi1 dphi2) convoluted with cuts
  //Then the epxressions 1, costh and costh^2 in sigmaij are substituted by these costh moments
  double costh = 0;
  int mode = 2;
  dofill_.doFill_ = 1;

  //evaluate the fixed order expansion of the resummed cross section
  double qt = (qtmn+qtmx)/2.;
  clock_t cbt = clock();
  bool fpc = true;
  ctint::calc(costh,m,qt,y,mode,f,opts.order_ct,fpc);
  clock_t cet = clock();

  //avoid nans
  if (isnan_ofast(f[0]))
    {
      cout << "nan in fpcintegr 2d " << endl;
      f[0]=0.;
//      if (opts.PDFerrors)
//	for (int i = 1; i < opts.totpdf; i++)
//	  f[i] = 0.;
      return 0;
    }
	   
  f[0] = f[0]*jac;
//  if (opts.PDFerrors)
//    for (int i = 1; i < opts.totpdf; i++)
//      f[i] = f[i]*jac;

  end_time = clock();
  if (opts.timeprofile)
    cout << setw (3) << "m" << setw(10) << m << setw(4) << "y" << setw(10) <<  y
	 << setw(8) << "result" << setw(10) << f[0]
	 << setw(10) << "tot time" << setw(10) << float( end_time - begin_time ) /  CLOCKS_PER_SEC
      	 << setw(10) << "qtint"  << setw(10) << float( qtet - qtbt ) /  CLOCKS_PER_SEC
	 << setw(10) << "ctint"  << setw(10) << float( cet - cbt ) /  CLOCKS_PER_SEC
	 << endl;
  return 0;
}

integrand_t fpcintegrandMC(const int &ndim, const double x[], const int &ncomp, double f[],
			  void* userdata, const int &nvec, const int &core,
			  double &weight, const int &iter)
{
  tell_to_grid_we_are_alive(10000);
  clock_t begin_time, end_time;

  begin_time = clock();

  if (!opts.fixedorder || phasespace::qtmin > 0)
    {
      f[0]=0.;
      return 0;
    }
  
//  if (opts.PDFerrors)
//    for (int i = 1; i < opts.totpdf; i++)
//      f[i] = 0.;
    
  //Call the counterterm
  dofill_.doFill_ = (opts.vegascollect || int(iter==cubacall::last_iter));
  double r[8];
  for (int i = 0; i < ndim; i++)
    r[i]=x[i];
  bool fpc = true;
  ctintmc::calc(r, f, opts.order_ct, fpc);
  
  //avoid nans
  if (isnan_ofast(f[0]))
    {
      f[0]=0.;
//      if (opts.PDFerrors)
//	for (int i = 1; i < opts.totpdf; i++)
//	  f[i] = 0.;
      return 0;
    }
  
  //fill histograms
  if (opts.vegascollect || (iter==cubacall::last_iter))
    for (int i = 0; i < opts.totpdf; i++)
      {
	double wt = weight*f[i];
	//if (opts.vegascollect) wt *= 1./double(cubacall::last_iter);
	hists_setpdf_(&i);
	hists_fill_(phasespace::p3, phasespace::p4, &wt);
      }

  end_time = clock();
  if (opts.timeprofile)
    cout << setw (3) << "m" << setw(10) << phasespace::m << setw(4) << "qt" << setw(10) <<  phasespace::qt
	 << setw (3) << "y" << setw(10) << phasespace::y << setw(4) << "costh" << setw(10) <<  phasespace::costh
	 << setw(8) << "result" << setw(10) << f[0]
	 << setw(10) << "tot time" << setw(10) << float( end_time - begin_time ) /  CLOCKS_PER_SEC
	 << endl;
  return 0;
}
