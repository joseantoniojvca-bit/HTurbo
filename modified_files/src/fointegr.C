#include "fointegr.h"

#include "dyres_interface.h"
#include "interface.h"
#include "omegaintegr.h"
#include "phasespace.h"
#include "settings.h"
#include "cubacall.h"
#include "isnan.h"
#include "qtint.h"
#include "ctint.h"
#include "ctmellin.h"
#include "vjint.h"
#include "vjloint.h"
#include "switch.h"
#include "sparserules.h"
#include "flush.h"
#include "evolnative.h"
#include "pegasus.h"
#include "pmom.h"
#include "ccoeff.h"

#include <iostream>
#include <iomanip>

#include <unistd.h>
#include <signal.h>
#include <sys/wait.h>

double fointegrand0d()
{
  double res = 0.;
  if (opts.cubacores == 0)
    {
      for (unsigned i = 0; i < sr::n3; i++)
	{
	  double x[3];
	  x[0] = {0.5+0.5*sr::x3[0+i*3]};
	  x[1] = {0.5+0.5*sr::x3[1+i*3]};
	  x[2] = {0.5+0.5*sr::x3[2+i*3]};
	  double f[1];
	  fointegrand3d(1, x, 1, f);
	  res += 0.5*0.5*0.5 * sr::w3[i] * f[0];
	}
      return res;
    }

  //openmp multithreading
  if (opts.threading == 0)
    {
#pragma omp parallel for reduction(+:res) num_threads(opts.cubacores) \
  copyin(a_param_,scale_,facscale_,qcdcouple_,vjloint::mur,vjloint::muf, \
	 mellinint::wn,mellinint::Np,mellinint::Nm,mellinint::wn_1,mellinint::Np_1,mellinint::Nm_1,mellinint::wn_2,mellinint::Np_2,mellinint::Nm_2,\
	 evolnative::ans,evolnative::am,evolnative::ap,evolnative::al,evolnative::be,evolnative::ab,\
	 evolnative::rmin,evolnative::rplus,evolnative::rqq,evolnative::rqg,evolnative::rgq,evolnative::rgg,\
	 evolnative::rmmqq,evolnative::rmmqg,evolnative::rmmgq,evolnative::rmmgg,evolnative::rmpqq,evolnative::rmpqg,evolnative::rmpgq,evolnative::rmpgg,\
	 evolnative::rpmqq,evolnative::rpmqg,evolnative::rpmgq,evolnative::rpmgg,evolnative::rppqq,evolnative::rppqg,evolnative::rppgq,evolnative::rppgg,\
	 pegasus::gli,pegasus::vai,pegasus::m3i,pegasus::m8i,pegasus::m15i,pegasus::m24i,pegasus::sgi,pegasus::p3i,pegasus::p8i,pegasus::p15i,pegasus::p24i,\
	 order_,painp_,hfpainp_,pacthr_,pabthr_,patthr_,moms_,hsums_,pns0_,psg0_,pns1_,psg1_,pns2_,psg2_,asg2_,lsg_,spsums_,u1sg_,r1sg_,u1hsg_,u2sg_,r2sg_,u2hsg_,u2ns_, \
	 ccoeff::C1qg,ccoeff::C1qq,ccoeff::C1qqb,ccoeff::C1qqp,ccoeff::C1qqbp,ccoeff::C2qg,ccoeff::C2qq,ccoeff::C2qqb,ccoeff::C2qqp,ccoeff::C2qqbp,ccoeff::C3qg,ccoeff::C3qq,ccoeff::C3qqb,ccoeff::C3qqp,ccoeff::C3qqbp,ccoeff::C1qg_reg,\
	 ccoeff::C1qg_1,ccoeff::C1qq_1,ccoeff::C1qqb_1,ccoeff::C1qqp_1,ccoeff::C1qqbp_1,ccoeff::C2qg_1,ccoeff::C2qq_1,ccoeff::C2qqb_1,ccoeff::C2qqp_1,ccoeff::C2qqbp_1,ccoeff::C3qg_1,ccoeff::C3qq_1,ccoeff::C3qqb_1,ccoeff::C3qqp_1,ccoeff::C3qqbp_1,ccoeff::C1qg_reg_1,\
	 ccoeff::C1qg_2,ccoeff::C1qq_2,ccoeff::C1qqb_2,ccoeff::C1qqp_2,ccoeff::C1qqbp_2,ccoeff::C2qg_2,ccoeff::C2qq_2,ccoeff::C2qqb_2,ccoeff::C2qqp_2,ccoeff::C2qqbp_2,ccoeff::C3qg_2,ccoeff::C3qq_2,ccoeff::C3qqb_2,ccoeff::C3qqp_2,ccoeff::C3qqbp_2,ccoeff::C1qg_reg_2,\
	 pmom::gamma1qq,pmom::gamma1qqb,pmom::gamma1qqp,pmom::gamma1qqbp,pmom::gamma1qg,pmom::gamma1gq,pmom::gamma1gg,pmom::gamma2qq,pmom::gamma2qqb,pmom::gamma2qqp,pmom::gamma2qqbp,pmom::gamma2qg,pmom::gamma2gq,pmom::gamma2gg,pmom::gamma3qq,pmom::gamma3qqb,pmom::gamma3qqp,pmom::gamma3qqbp,pmom::gamma3qg,pmom::gamma3gq,pmom::gamma3gg,\
	 pmom::gamma1qq_1,pmom::gamma1qqb_1,pmom::gamma1qqp_1,pmom::gamma1qqbp_1,pmom::gamma1qg_1,pmom::gamma1gq_1,pmom::gamma1gg_1,pmom::gamma2qq_1,pmom::gamma2qqb_1,pmom::gamma2qqp_1,pmom::gamma2qqbp_1,pmom::gamma2qg_1,pmom::gamma2gq_1,pmom::gamma2gg_1,pmom::gamma3qq_1,pmom::gamma3qqb_1,pmom::gamma3qqp_1,pmom::gamma3qqbp_1,pmom::gamma3qg_1,pmom::gamma3gq_1,pmom::gamma3gg_1,\
	 pmom::gamma1qq_2,pmom::gamma1qqb_2,pmom::gamma1qqp_2,pmom::gamma1qqbp_2,pmom::gamma1qg_2,pmom::gamma1gq_2,pmom::gamma1gg_2,pmom::gamma2qq_2,pmom::gamma2qqb_2,pmom::gamma2qqp_2,pmom::gamma2qqbp_2,pmom::gamma2qg_2,pmom::gamma2gq_2,pmom::gamma2gg_2,pmom::gamma3qq_2,pmom::gamma3qqb_2,pmom::gamma3qqp_2,pmom::gamma3qqbp_2,pmom::gamma3qg_2,pmom::gamma3gq_2,pmom::gamma3gg_2)
      for (unsigned i = 0; i < sr::n3; i++)
	{
	  double x[3];
	  x[0] = {0.5+0.5*sr::x3[0+i*3]};
	  x[1] = {0.5+0.5*sr::x3[1+i*3]};
	  x[2] = {0.5+0.5*sr::x3[2+i*3]};
	  double f[1];
	  fointegrand3d(1, x, 1, f);
	  res += 0.5*0.5*0.5 * sr::w3[i] * f[0];
	}
      return res;
    }

  //fork wait multithreading
  if (opts.threading == 1)
    {
      flushall();
      int fd[2];
      if (pipe(fd) < 0)
	{
	  cout << "Error in fointegrand0d fork/wait: could not create pipe" << endl;
	  exit(-1);
	}
      
      for (unsigned i = 0; i < sr::n3; i++)
	{
	  if (i >= opts.cubacores)
	    wait(NULL);
      
	  pid_t id = fork();
	  if (id == 0)
	    {
	      close(fd[0]);
	      double x[3];
	      x[0] = {0.5+0.5*sr::x3[0+i*3]};
	      x[1] = {0.5+0.5*sr::x3[1+i*3]};
	      x[2] = {0.5+0.5*sr::x3[2+i*3]};
	      double f[1];
	      fointegrand3d(1, x, 1, f);
	      double partres = 0.5*0.5*0.5 * sr::w3[i] * f[0];
	      write(fd[1], &partres, sizeof partres);
	      exit(0);
	    }
	  else if (id < 0)
	    {
	      cout << "Error in fointegrand0d fork failed for process " << i << " pid " << id << endl;
	      exit (-1);
	    }
	}

      //wait for all children to finish
      int status;
      pid_t wpid;
      while ((wpid = wait(&status)) > 0)
        if (status < 0)
	  {
	    cout << "Process " << wpid << " terminated with status " << status << endl;
	    exit(-1);
	  }
  
      //close the write end of the pipe
      close(fd[1]);
      for(int i = 0; i < sr::n3; i++)
	{
	  double partres;
	  int nbytes = read(fd[0], &partres, sizeof partres);
	  if (nbytes > 0)
	      res += partres;
	  else
	    {
	      cout << "Error in fointegrand0d fork/wait: nothing on the pipe" << endl;
	      exit(-1);
	    }
	}
      //close the read end of the pipe
      close(fd[0]);
    }
  return res;
}

int fointegrand3d_cubature_v(unsigned ndim, long unsigned npts, const double x[], void *data, unsigned ncomp, double f[])
{
  //  cout << "parallel " << npts << endl;
#pragma omp parallel for num_threads(opts.cubacores) copyin(a_param_,scale_,facscale_,qcdcouple_)
  for (unsigned i = 0; i < npts; i++)
    {
      // evaluate the integrand for npts points
      double xi[ndim];
      double fi[ncomp];
      for (unsigned j = 0; j < ndim; j++)
	xi[j] = x[i*ndim + j];

      fointegrand3d(ndim, xi, ncomp, fi);
      
      for (unsigned k = 0; k < ncomp; ++k)
	f[i*ncomp + k] = fi[k];
    }
  tell_to_grid_we_are_alive(1000);
  return 0;
}

int fointegrand3d_cubature(unsigned ndim, const double x[], void *data, unsigned ncomp, double f[])
{
  fointegrand3d(ndim, x, ncomp, f);
  tell_to_grid_we_are_alive(1000);
  return 0;
}
integrand_t fointegrand3d(const int &ndim, const double x[], const int &ncomp, double f[])
//Generates the phase space 4 vectors
//Calculates the V+j LO/NLO and CT integrands as a function of m, qt, and y. Only for resummation, this routine is never called for fixedorder = true
{
  clock_t begin_time, end_time;

  begin_time = clock();

  if (opts.fixedorder)
    {
      f[0]=0.;
      return 0;
    }
  
  //Jacobian of the change of variables from the unitary square x[2] to the m, qt boundaries
  double jac = 1.;
  bool status = true;
  
  double r2[2] = {x[0], x[1]};
  status = phasespace::gen_my(r2, jac, true, false);   //qtcut = true, qtswitching = false
  if (!status)
    {
      f[0] = 0.;
      return 0;
    }
  double r = x[2];

  //Generate the boson transverse momentum between the integration boundaries
  double qtcut = max(opts.qtcut,opts.xqtcut*phasespace::m);
  double qtmn = max(qtcut, phasespace::qtmin);
  //phasespace::calcexpy();
  //double cosh2y=pow((phasespace::exppy+phasespace::expmy)*0.5,2);
  //double kinqtlim = sqrt(max(0.,pow(pow(opts.sroot,2)+phasespace::m2,2)/(4*pow(opts.sroot,2)*cosh2y)-phasespace::m2)); //introduced max to avoid neqative argument of sqrt when y=ymax
  //double kinqtlim = sqrt(max(0.,pow(pow(opts.sroot,2)+phasespace::m2,2)/(4*pow(opts.sroot,2))-phasespace::m2)); //introduced max to avoid negative argument of sqrt when y=ymax
  double kinqtlim = 1e10; //There should not be any kinematic limit on qt, since the counterterm is evaluated with born level kinematic
  double switchqtlim = switching::qtlimit(phasespace::m);
  //double qtlim = min(kinqtlim, switchqtlim);
  double qtlim = kinqtlim; //There should not be any switching limit on qt, since the v+jet is allowed to go to the kinematic limit
  double qtmx = min(qtlim, phasespace::qtmax);
  if (qtmn >= qtmx)
    {
      f[0] = 0.;
      return 0;
    }
  status = phasespace::gen_qt(x[2], jac, qtlim, true);
  if (!status)
    {
      f[0] = 0.;
      return 0;
    }
  //jac = jac *2.*phasespace::qt;

  //double qtmn = 1e-2;
  //double qtmx = max(opts.qtcut,opts.xqtcut*phasespace::m);
  //qtweight_(r,qtmn,qtmx,phasespace::qt,jac);
  //phasespace::set_qt(phasespace::qt);

  double m = phasespace::m;
  double qt = phasespace::qt;
  double y = phasespace::y;
  phasespace::calcexpy();
  //phasespace::calcmt();
  phasespace::set_phiV(0.);
  
  //evaluate the Vj (N)LO cross section and the CT in full-lepton phase space

  omegaintegr::genV4p();
  double vj = vjint::vint(m,qt,y);

  qtint::calc(m,qt,0,1);
  omegaintegr::genV4p();
  double costh = 0;
  int mode = 1;
  dofill_.doFill_ = 1;
  if (opts.ctmellin)
    ctmellin::calc(m,qt,y,mode,f);
  else
    ctint::calc(costh,m,qt,y,mode,f,min(2,opts.order_ct));
  double ctfull = f[0];

  //cout << " vj " << vj << " ct " << ctfull*2*qt << " diff " << vj+ctfull*2*qt << endl;
  f[0] =  vj+ctfull*2*qt;
  
  f[0] = f[0]*jac;

  if (isnan_ofast(f[0]))
    {
      cout << m << " " << qt << " " << y << " nan in fointegr:vint() " << endl;
      f[0] = 0.;  //avoid nans
    }
	   
  if (isnan_ofast(jac))
    {
      cout << m << " " << qt << " " << y << ": nan in gen_my() " << endl;
      f[0] = 0.;  //avoid nans
    }
  else

  end_time = clock();
  if (opts.timeprofile)
    cout << setw (3) << "m" << setw(10) << m << setw(4) << "qt" << setw(10) <<  qt << setw(4) << "y" << setw(10) <<  y
	 << setw(8) << "result" << setw(12) << f[0]
	 << setw(10) << "tot time" << setw(10) << float( end_time - begin_time ) /  CLOCKS_PER_SEC
	 << endl;
  
  tell_to_grid_we_are_alive(1000);
  return 0;
}

double fointegrand3d_smolyak(int ndim, double x[])
{
  int ncomp = 1;
  double f[ncomp];
  fointegrand3d(ndim, x, ncomp, f);
  return f[0];
}
