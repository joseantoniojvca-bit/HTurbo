#include "evoltrunc.h"

#include "hell.h"
#include "parton.h"
#include "alphas.h"
#include "alphaqed.h"
#include "mellinpdf.h"
#include "settings.h"
#include "scales.h"
#include "pdf.h"
#include "phasespace.h"
#include "mesq.h"

using namespace parton;

int evoltrunc::dim;
complex <double> *evoltrunc::gli;
complex <double> *evoltrunc::vai;
complex <double> *evoltrunc::m3i;
complex <double> *evoltrunc::m8i;
complex <double> *evoltrunc::m15i;
complex <double> *evoltrunc::m24i;
complex <double> *evoltrunc::sgi;
complex <double> *evoltrunc::p3i;
complex <double> *evoltrunc::p8i;
complex <double> *evoltrunc::p15i;
complex <double> *evoltrunc::p24i;
complex <double> *evoltrunc::fti;
complex <double> *evoltrunc::cxn;

complex <double> *evoltrunc::gli_muf;
complex <double> *evoltrunc::vai_muf;
complex <double> *evoltrunc::m3i_muf;
complex <double> *evoltrunc::m8i_muf;
complex <double> *evoltrunc::m15i_muf;
complex <double> *evoltrunc::m24i_muf;
complex <double> *evoltrunc::sgi_muf;
complex <double> *evoltrunc::p3i_muf;
complex <double> *evoltrunc::p8i_muf;
complex <double> *evoltrunc::p15i_muf;
complex <double> *evoltrunc::p24i_muf;
complex <double> *evoltrunc::fti_muf;

complex <double> *evoltrunc::rp;
complex <double> *evoltrunc::rm;
cmatrix *evoltrunc::Ep;
cmatrix *evoltrunc::Em;
cmatrix *evoltrunc::U1;
cmatrix *evoltrunc::U2;
cmatrix *evoltrunc::U3;
cmatrix *evoltrunc::U;

complex <double> *evoltrunc::R0NS;
complex <double> *evoltrunc::U1NSp;
complex <double> *evoltrunc::U1NSm;
complex <double> *evoltrunc::U1NSv;
complex <double> *evoltrunc::U2NSp;
complex <double> *evoltrunc::U2NSm;
complex <double> *evoltrunc::U2NSv;
complex <double> *evoltrunc::U3NSp;
complex <double> *evoltrunc::U3NSm;
complex <double> *evoltrunc::U3NSv;
complex <double> *evoltrunc::UNSp;
complex <double> *evoltrunc::UNSm;
complex <double> *evoltrunc::UNSv;

const int evoltrunc::maxord = 20;

const cmatrix evoltrunc::I = cmatrix(1.,0.,0.,1.);
const cmatrix evoltrunc::N = cmatrix(0.,0.,0.,0.);

double evoltrunc::beta0;
double evoltrunc::beta1;
double evoltrunc::beta2;
double evoltrunc::beta3;

double evoltrunc::b1;
double evoltrunc::b2;
double evoltrunc::b3;

double evoltrunc::beta(double as, int o)
{
  double a2  = pow(as,2);
  double a3  = pow(as,3);
  double a4  = pow(as,4);
  double a5  = pow(as,5);
  double a6  = pow(as,6);
  double a7  = pow(as,7);
  
  double b0 = beta0;
  double b1 = beta1;
  double b2 = beta2;
  double b3 = beta3;
  double b4 = beta4;
  //if (opts.rge)
  //  {
  //    double csi;
  //    if (opts.commute)
  //	csi = log(scales::ren/scales::fac);
  //    else
  //      csi = log(scales::ren/scales::res);
  //    //csi = log(scales::res/scales::fac);
  //
  //    //return real(alphas::fbetabar(as,csi,o+1));
  //
  //    double lk = csi;
  //    double lk2 = lk * lk;
  //    double lk3 = lk * lk2;
  //    double lk4 = lk * lk3;
  //
  //    int n = o+1;
  //    if (n == 1) return -a2*b0;
  //    if (n == 2) return -a2*b0 - a3*b1 + 5.*a4*b0*b1*lk + 10.*a4*pow(b0,3)*lk2;
  //    if (n == 3) return -a2*b0 - a3*b1 - a4*b2 + 3.*a5*pow(b1,2)*lk + 6.*a5*b0*b2*lk + 29.*a5*pow(b0,2)*b1*lk2 + 28.*a5*pow(b0,4)*lk3;
  //    if (n == 4) return -a2*b0 - a3*b1 - a4*b2 - a5*b3 + 7.*a6*b1*b2*lk + 7.*a6*b0*b3*lk + 77./2.*a6*b0*pow(b1,2)*lk2 + 38.*a6*pow(b0,2)*b2*lk2 + 563./3.*a6*pow(b0,3)*b1*lk3 + 532./3.*a6*pow(b0,5)*lk4;
  //    if (n == 5) return -a2*b0 - a3*b1 - a4*b2 - a5*b3 - a6*b4 + 4.*pow(b2,2)*a7*lk + 8.*b1*b3*a7*lk + 33./2.*pow(b1,3)*a7*lk2 + 8.*b0*b4*a7*lk + 98.*b0*b1*b2*a7*lk2 + 48.*pow(b0,2)*b3*a7*lk2 + 1193./3.*pow(b0,2)*pow(b1,2)*a7*lk3 + 260*pow(b0,3)*b2*a7*lk3 + 2633./2.*pow(b0,4)*b1*a7*lk4;
  //  }
  //else
    {
      int n = o+1;
      if (n == 1) return -a2*b0;
      if (n == 2) return -a2*b0 -a3*b1;
      if (n == 3) return -a2*b0 -a3*b1 -a4*b2;
      if (n == 4) return -a2*b0 -a3*b1 -a4*b2 -a5*b3;
      if (n == 5) return -a2*b0 -a3*b1 -a4*b2 -a5*b3 -a6*b4;
    }
  return 0;
}

void evoltrunc::init()
{
  //if (opts.mellin1d)
  //  dim = mellinint::mdim;
  //else
  dim = 2*mellinint::mdim;

  //Evolve the PDFs in msbar scheme, because PDFs in LHAPDF are in the msbar scheme
  beta0 = alphas::beta0_msbar/M_PI;
  beta1 = alphas::beta1_msbar/pi2;
  beta2 = alphas::beta2_msbar/pi3;
  beta3 = alphas::beta3_msbar/pi4;
  //beta0 = resconst::beta0/M_PI;
  //beta1 = resconst::beta1/pi2;
  //beta2 = resconst::beta2/pi3;
  //beta3 = resconst::beta3/pi4;
  
  b1 = beta1/beta0;
  b2 = beta2/beta0;
  b3 = beta3/beta0;

  //if (opts.melup <= 1)
  //  {
  //    allocate();
  //    //Cannot update on init, because it requires pmom::init, which requires pegasus psg common block
  //    update();
  //  }
}
void evoltrunc::release()
{
  //if (opts.melup <= 1)
  //free();
}

//Allocate PDF moments at the starting scale
void evoltrunc::allocate()
{
  gli  = new complex <double> [2*dim];
  vai  = new complex <double> [2*dim];
  m3i  = new complex <double> [2*dim];
  m8i  = new complex <double> [2*dim];
  m15i = new complex <double> [2*dim];
  m24i = new complex <double> [2*dim];
  sgi  = new complex <double> [2*dim];
  p3i  = new complex <double> [2*dim];
  p8i  = new complex <double> [2*dim];
  p15i = new complex <double> [2*dim];
  p24i = new complex <double> [2*dim];
  fti  = new complex <double> [2*dim];
  cxn = new complex <double> [2*dim];

  gli_muf  = new complex <double> [2*dim];
  vai_muf  = new complex <double> [2*dim];
  m3i_muf  = new complex <double> [2*dim];
  m8i_muf  = new complex <double> [2*dim];
  m15i_muf = new complex <double> [2*dim];
  m24i_muf = new complex <double> [2*dim];
  sgi_muf  = new complex <double> [2*dim];
  p3i_muf  = new complex <double> [2*dim];
  p8i_muf  = new complex <double> [2*dim];
  p15i_muf = new complex <double> [2*dim];
  p24i_muf = new complex <double> [2*dim];
  fti_muf  = new complex <double> [2*dim];
  
  rp = new complex <double> [2*dim];
  rm = new complex <double> [2*dim];
  Ep = new cmatrix [2*dim];
  Em = new cmatrix [2*dim];
  U1 = new cmatrix [2*dim];
  U2 = new cmatrix [2*dim];
  U3 = new cmatrix [2*dim];
  U = new cmatrix [2*dim*maxord];

  R0NS = new complex <double> [2*dim];
  U1NSp = new complex <double> [2*dim];
  U1NSm = new complex <double> [2*dim];
  U1NSv = new complex <double> [2*dim];
  U2NSp = new complex <double> [2*dim];
  U2NSm = new complex <double> [2*dim];
  U2NSv = new complex <double> [2*dim];
  U3NSp = new complex <double> [2*dim];
  U3NSm = new complex <double> [2*dim];
  U3NSv = new complex <double> [2*dim];
  UNSp = new complex <double> [2*dim*maxord];
  UNSm = new complex <double> [2*dim*maxord];
  UNSv = new complex <double> [2*dim*maxord];
}

void evoltrunc::free()
{
  delete[] gli;
  delete[] vai;
  delete[] m3i;
  delete[] m8i;
  delete[] m15i;
  delete[] m24i;
  delete[] sgi;
  delete[] p3i;
  delete[] p8i;
  delete[] p15i;
  delete[] p24i;
  delete[] fti;
  delete[] cxn;

  delete[] gli_muf;
  delete[] vai_muf;
  delete[] m3i_muf;
  delete[] m8i_muf;
  delete[] m15i_muf;
  delete[] m24i_muf;
  delete[] sgi_muf;
  delete[] p3i_muf;
  delete[] p8i_muf;
  delete[] p15i_muf;
  delete[] p24i_muf;
  delete[] fti_muf;
  
  delete[] rp;
  delete[] rm;
  delete[] Ep;
  delete[] Em;
  delete[] U1;
  delete[] U2;
  delete[] U3;
  delete[] U;

  delete[] R0NS;
  delete[] U1NSp;
  delete[] U1NSm;
  delete[] U1NSv;
  delete[] U2NSp;
  delete[] U2NSm;
  delete[] U2NSv;
  delete[] U3NSp;
  delete[] U3NSm;
  delete[] U3NSv;
  delete[] UNSp;
  delete[] UNSm;
  delete[] UNSv;
}

void evoltrunc::update()
{
  update_pdf();
  update_RU();
}

//update PDFs at the factorisation scale
void evoltrunc::update_pdf()
{
  double facscale = scales::fac;
  
  mellinpdf::allocate();
  mellinpdf::evalpdfs(facscale, phasespace::m, phasespace::y);
  mellinpdf::update_mellin();
  mellinpdf::transform();

  complex <double> qp[5];
  complex <double> qm[5];
  complex <double> gl;
  for (int sign = mesq::positive; sign <= mesq::negative; sign++)
    for (int n = 0; n < dim; n++)
      {
	int i = index(n,sign);

	int nn;
	//if (opts.mellin1d)
	//  nn = n;
	//else
	if (n < mellinint::mdim)
	  nn = n;
	else
	  nn = n-mellinint::mdim;
	
	//if (opts.mellin1d)
	//  if (sign == mesq::positive)
	//    {
	//      gl    = mellinpdf::GL[nn];
	//      qp[0] = mellinpdf::UP[nn] + mellinpdf::UB[nn];
	//      qp[1] = mellinpdf::DO[nn] + mellinpdf::DB[nn];
	//      qp[2] = mellinpdf::ST[nn] + mellinpdf::SB[nn];
	//      qp[3] = mellinpdf::CH[nn] + mellinpdf::CB[nn];
	//      qp[4] = mellinpdf::BO[nn] + mellinpdf::BB[nn];
	//      qm[0] = mellinpdf::UP[nn] - mellinpdf::UB[nn];
	//      qm[1] = mellinpdf::DO[nn] - mellinpdf::DB[nn];
	//      qm[2] = mellinpdf::ST[nn] - mellinpdf::SB[nn];
	//      qm[3] = mellinpdf::CH[nn] - mellinpdf::CB[nn];
	//      qm[4] = mellinpdf::BO[nn] - mellinpdf::BB[nn];
	//      fti_muf[i]  = mellinpdf::FT[nn];
	//    }
	//  else
	//    {
	//      gl    = conj(mellinpdf::GL[nn]);
	//      qp[0] = conj(mellinpdf::UP[nn] + mellinpdf::UB[nn]);
	//      qp[1] = conj(mellinpdf::DO[nn] + mellinpdf::DB[nn]);
	//      qp[2] = conj(mellinpdf::ST[nn] + mellinpdf::SB[nn]);
	//      qp[3] = conj(mellinpdf::CH[nn] + mellinpdf::CB[nn]);
	//      qp[4] = conj(mellinpdf::BO[nn] + mellinpdf::BB[nn]);
	//      qm[0] = conj(mellinpdf::UP[nn] - mellinpdf::UB[nn]);
	//      qm[1] = conj(mellinpdf::DO[nn] - mellinpdf::DB[nn]);
	//      qm[2] = conj(mellinpdf::ST[nn] - mellinpdf::SB[nn]);
	//      qm[3] = conj(mellinpdf::CH[nn] - mellinpdf::CB[nn]);
	//      qm[4] = conj(mellinpdf::BO[nn] - mellinpdf::BB[nn]);
	//      fti_muf[i]  = conj(mellinpdf::FT[nn]);
	//    }
	//else
	  if (n < mellinint::mdim)
	    if (sign == mesq::positive)
	      {
		gl    = mellinpdf::GL_1[nn];
		qp[0] = mellinpdf::UP_1[nn] + mellinpdf::UB_1[nn];
		qp[1] = mellinpdf::DO_1[nn] + mellinpdf::DB_1[nn];
		qp[2] = mellinpdf::ST_1[nn] + mellinpdf::SB_1[nn];
		qp[3] = mellinpdf::CH_1[nn] + mellinpdf::CB_1[nn];
		qp[4] = mellinpdf::BO_1[nn] + mellinpdf::BB_1[nn];
		qm[0] = mellinpdf::UP_1[nn] - mellinpdf::UB_1[nn];
		qm[1] = mellinpdf::DO_1[nn] - mellinpdf::DB_1[nn];
		qm[2] = mellinpdf::ST_1[nn] - mellinpdf::SB_1[nn];
		qm[3] = mellinpdf::CH_1[nn] - mellinpdf::CB_1[nn];
		qm[4] = mellinpdf::BO_1[nn] - mellinpdf::BB_1[nn];
		fti_muf[i]  = mellinpdf::FT_1[nn];
	      }
	    else
	      {
		gl    = conj(mellinpdf::GL_1[nn]);
		qp[0] = conj(mellinpdf::UP_1[nn] + mellinpdf::UB_1[nn]);
		qp[1] = conj(mellinpdf::DO_1[nn] + mellinpdf::DB_1[nn]);
		qp[2] = conj(mellinpdf::ST_1[nn] + mellinpdf::SB_1[nn]);
		qp[3] = conj(mellinpdf::CH_1[nn] + mellinpdf::CB_1[nn]);
		qp[4] = conj(mellinpdf::BO_1[nn] + mellinpdf::BB_1[nn]);
		qm[0] = conj(mellinpdf::UP_1[nn] - mellinpdf::UB_1[nn]);
		qm[1] = conj(mellinpdf::DO_1[nn] - mellinpdf::DB_1[nn]);
		qm[2] = conj(mellinpdf::ST_1[nn] - mellinpdf::SB_1[nn]);
		qm[3] = conj(mellinpdf::CH_1[nn] - mellinpdf::CB_1[nn]);
		qm[4] = conj(mellinpdf::BO_1[nn] - mellinpdf::BB_1[nn]);
		fti_muf[i]  = conj(mellinpdf::FT_1[nn]);
	      }
	  else
	    if (sign == mesq::positive)
	      {
		gl    = mellinpdf::GL_2[nn];
		qp[0] = mellinpdf::UP_2[nn] + mellinpdf::UB_2[nn];
		qp[1] = mellinpdf::DO_2[nn] + mellinpdf::DB_2[nn];
		qp[2] = mellinpdf::ST_2[nn] + mellinpdf::SB_2[nn];
		qp[3] = mellinpdf::CH_2[nn] + mellinpdf::CB_2[nn];
		qp[4] = mellinpdf::BO_2[nn] + mellinpdf::BB_2[nn];
		qm[0] = mellinpdf::UP_2[nn] - mellinpdf::UB_2[nn];
		qm[1] = mellinpdf::DO_2[nn] - mellinpdf::DB_2[nn];
		qm[2] = mellinpdf::ST_2[nn] - mellinpdf::SB_2[nn];
		qm[3] = mellinpdf::CH_2[nn] - mellinpdf::CB_2[nn];
		qm[4] = mellinpdf::BO_2[nn] - mellinpdf::BB_2[nn];
		fti_muf[i]  = mellinpdf::FT_2[nn];
	      }
	    else
	      {
		gl    = conj(mellinpdf::GL_2[nn]);
		qp[0] = conj(mellinpdf::UP_2[nn] + mellinpdf::UB_2[nn]);
		qp[1] = conj(mellinpdf::DO_2[nn] + mellinpdf::DB_2[nn]);
		qp[2] = conj(mellinpdf::ST_2[nn] + mellinpdf::SB_2[nn]);
		qp[3] = conj(mellinpdf::CH_2[nn] + mellinpdf::CB_2[nn]);
		qp[4] = conj(mellinpdf::BO_2[nn] + mellinpdf::BB_2[nn]);
		qm[0] = conj(mellinpdf::UP_2[nn] - mellinpdf::UB_2[nn]);
		qm[1] = conj(mellinpdf::DO_2[nn] - mellinpdf::DB_2[nn]);
		qm[2] = conj(mellinpdf::ST_2[nn] - mellinpdf::SB_2[nn]);
		qm[3] = conj(mellinpdf::CH_2[nn] - mellinpdf::CB_2[nn]);
		qm[4] = conj(mellinpdf::BO_2[nn] - mellinpdf::BB_2[nn]);
		fti_muf[i]  = conj(mellinpdf::FT_2[nn]);
	      }
	
	//arrays of non-singlet and singlet quark combinations for N_f = 3 (and 4 and 5)
	//defined as in Eq. (2.16) of hep-ph/0408244
	gli_muf[i]  = gl;	
	vai_muf[i]  = qm[0]+qm[1]+qm[2]+qm[3]+qm[4];
	m3i_muf[i]  = qm[0]-qm[1];
	m8i_muf[i]  = qm[0]+qm[1]-2.*qm[2];
	m15i_muf[i] = qm[0]+qm[1]+qm[2]-3.*qm[3];
	m24i_muf[i] = qm[0]+qm[1]+qm[2]+qm[3]-4.*qm[4];
	
	sgi_muf[i]  = qp[0]+qp[1]+qp[2]+qp[3]+qp[4];
	p3i_muf[i]  = qp[0]-qp[1];
	p8i_muf[i]  = qp[0]+qp[1]-2.*qp[2];
	p15i_muf[i] = qp[0]+qp[1]+qp[2]-3.*qp[3];
	p24i_muf[i] = qp[0]+qp[1]+qp[2]+qp[3]-4.*qp[4];
      }

  mellinpdf::free();

  //if (opts.mellin1d)
  //  for (int sign = mesq::positive; sign <= mesq::negative; sign++)
  //    for (int i = 0; i < mellinint::mdim; i++)
  //	{
  //	  int idx = index(i,sign);
  //	  if (sign == mesq::positive)
  //	    cxn[idx] = mellinint::Np[i];
  //	  else
  //	    cxn[idx] = mellinint::Nm[i];
  //	}
  //else
    for (int sign = mesq::positive; sign <= mesq::negative; sign++)
      for (int beam = 1; beam <= 2; beam++)
	for (int i = 0; i < mellinint::mdim; i++)
	  {
	    int idx = index(i,beam,sign);
	    if (beam == 1)
	      if (sign == mesq::positive)
		cxn[idx] = mellinint::Np_1[i];
	      else
		cxn[idx] = mellinint::Nm_1[i];
	    else
	      if (sign == mesq::positive)
		cxn[idx] = mellinint::Np_2[i];
	      else
		cxn[idx] = mellinint::Nm_2[i];
	  }
}

//retrieve evolved PDFs at the scale mub
void evoltrunc::retrieve_pdf_mub()
{
  //cout << "retrieve_pdf_mub" << endl;
  for (int sign = mesq::positive; sign <= mesq::negative; sign++)
    for (int n = 0; n < dim; n++)
      {

	int nn;
	//if (opts.mellin1d)
	//  nn = n;
	//else
	if (n < mellinint::mdim)
	  nn = n;
	else
	  nn = n-mellinint::mdim;
	
	//if (opts.mellin1d)
	//  if (sign == mesq::positive)
	//    pdfevol::retrieve1d_pos(nn);
	//  else
	//    pdfevol::retrieve1d_pos(nn+mellinint::mdim);
	//else
	if (n < mellinint::mdim)
	  if (sign == mesq::positive)
	    pdfevol::retrieve_beam1(nn);
	  else
	    pdfevol::retrieve_beam1(nn+mellinint::mdim);
	else
	  if (sign == mesq::positive)
	    pdfevol::retrieve_beam2_pos(nn);
	  else
	    pdfevol::retrieve_beam2_pos(nn+mellinint::mdim);

	complex<double>* fn,ft;
	if (n < mellinint::mdim)
	  {
	    fn = pdfevol::fn1;
	    ft = pdfevol::fotn1;
	  }
	else
	  {
	    fn = pdfevol::fn2;
	    ft = pdfevol::fotn2;
	  }
	
	complex <double> qp[5];
	complex <double> qm[5];

	//here flip back antiprotons and nuclear PDFs
	complex <double> fnn[2*MAXNF+1];
	copy(fn, fn+2*MAXNF+1, fnn);
	//if (n < mellinint::mdim) //beam 1
	//  {
	//    if (opts.ih1 == -1)
	//      {
	//	fn[bb] = fnn[b ];
	//	fn[cb] = fnn[c ];
	//	fn[sb] = fnn[s ];
	//	fn[db] = fnn[d ];
	//	fn[ub] = fnn[u ];
	//	fn[u ] = fnn[ub];
	//	fn[d ] = fnn[db];
	//	fn[s ] = fnn[sb];
	//	fn[c ] = fnn[cb];
	//	fn[b ] = fnn[bb];
	//      }
	//  }
	//else //beam 2
	//  {
	//    if (opts.ih2 == -1)
	//      {
	//	fn[bb] = fnn[b ];
	//	fn[cb] = fnn[c ];
	//	fn[sb] = fnn[s ];
	//	fn[db] = fnn[d ];
	//	fn[ub] = fnn[u ];
	//	fn[u ] = fnn[ub];
	//	fn[d ] = fnn[db];
	//	fn[s ] = fnn[sb];
	//	fn[c ] = fnn[cb];
	//	fn[b ] = fnn[bb];
	//      }
	//  }
	
//	if (n < mellinint::mdim || opts.ih2 == 1)
//	  {
	    qp[0] = fn[u] + fn[ub];
	    qp[1] = fn[d] + fn[db];
	    qp[2] = fn[s] + fn[sb];
	    qp[3] = fn[c] + fn[cb];
	    qp[4] = fn[b] + fn[bb];

	    qm[0] = fn[u] -  fn[ub];
	    qm[1] = fn[d] -  fn[db];
	    qm[2] = fn[s] -  fn[sb];
	    qm[3] = fn[c] -  fn[cb];
	    qm[4] = fn[b] -  fn[bb];
//	  }
//	else //flip back pbar -> p
//	  {
//	    qp[0] = fn[ub] + fn[u];
//	    qp[1] = fn[db] + fn[d];
//	    qp[2] = fn[sb] + fn[s];
//	    qp[3] = fn[cb] + fn[c];
//	    qp[4] = fn[bb] + fn[b];
//
//	    qm[0] = fn[ub] -  fn[u];
//	    qm[1] = fn[db] -  fn[d];
//	    qm[2] = fn[sb] -  fn[s];
//	    qm[3] = fn[cb] -  fn[c];
//	    qm[4] = fn[bb] -  fn[b];
//	  }

	//arrays of non-singlet and singlet quark combinations for N_f = 3 (and 4 and 5)
	//defined as in Eq. (2.16) of hep-ph/0408244
	int i = index(n,sign);
	//cout << "evoltrunc " << i << "  " << fn[g] << endl;
	//cout << i << "  " << gli[i] << "  " << fn[g] << endl;
	gli[i]  = fn[g];
	vai[i]  = qm[0]+qm[1]+qm[2]+qm[3]+qm[4];
	m3i[i]  = qm[0]-qm[1];
	m8i[i]  = qm[0]+qm[1]-2.*qm[2];
	m15i[i] = qm[0]+qm[1]+qm[2]-3.*qm[3];
	m24i[i] = qm[0]+qm[1]+qm[2]+qm[3]-4.*qm[4];
	
	sgi[i]  = qp[0]+qp[1]+qp[2]+qp[3]+qp[4];
	p3i[i]  = qp[0]-qp[1];
	p8i[i]  = qp[0]+qp[1]-2.*qp[2];
	p15i[i] = qp[0]+qp[1]+qp[2]-3.*qp[3];
	p24i[i] = qp[0]+qp[1]+qp[2]+qp[3]-4.*qp[4];
	fti[i]  = ft;
      }
}

//retrieve evolved PDFs at the scale muf
void evoltrunc::retrieve_pdf_muf()
{
  for (int i = 0; i < dim*2; i++)
    {
	gli[i]  = gli_muf[i];  
	vai[i]  = vai_muf[i];  
	m3i[i]  = m3i_muf[i];  
	m8i[i]  = m8i_muf[i];  
	m15i[i] = m15i_muf[i]; 
	m24i[i] = m24i_muf[i]; 
	sgi[i]  = sgi_muf[i];  
	p3i[i]  = p3i_muf[i];  
	p8i[i]  = p8i_muf[i];  
	p15i[i] = p15i_muf[i]; 
	p24i[i] = p24i_muf[i]; 
	fti[i]  = fti_muf[i]; 
    }
}

void evoltrunc::update_RU()
{
  int n,nn,sign;
  for (int i = 0; i < dim*2; i++)
    {
      n    = (i<dim)? i : i-dim;
      sign = (i<dim)? mesq::positive : mesq::negative;
      cmatrix g1sg,g2sg,g3sg,g4sg;
      complex <double> g1nsp,g1nsm,g1nsv;
      complex <double> g2nsp,g2nsm,g2nsv;
      complex <double> g3nsp,g3nsm,g3nsv;
      complex <double> g4nsp,g4nsm,g4nsv;
      if (opts.mellin1d)
	{
	  if (n < mellinint::mdim)
	    nn = n;
	  else
	    nn = n-mellinint::mdim;
	  
	  //Singlet
	  g1sg = gamma_LO(pmom::index(nn,sign))   ;
	  g2sg = gamma_NLO(pmom::index(nn,sign))  ;
	  g3sg = gamma_NNLO(pmom::index(nn,sign)) ;
	  g4sg = gamma_NNNLO(pmom::index(nn,sign));

	  //Non singlet
	  g1nsp = gammaNSp_LO(pmom::index(nn,sign));
	  g1nsm = gammaNSm_LO(pmom::index(nn,sign));
	  g1nsv = gammaNSv_LO(pmom::index(nn,sign));
	  g2nsp = gammaNSp_NLO(pmom::index(nn,sign));
	  g2nsm = gammaNSm_NLO(pmom::index(nn,sign));
	  g2nsv = gammaNSv_NLO(pmom::index(nn,sign));
	  g3nsp = gammaNSp_NNLO(pmom::index(nn,sign));
	  g3nsm = gammaNSm_NNLO(pmom::index(nn,sign));
	  g3nsv = gammaNSv_NNLO(pmom::index(nn,sign));
	  g4nsp = gammaNSp_NNNLO(pmom::index(nn,sign));
	  g4nsm = gammaNSm_NNNLO(pmom::index(nn,sign));
	  g4nsv = gammaNSv_NNNLO(pmom::index(nn,sign));
	}
      else
	if (n < mellinint::mdim)
	  {
	    nn = n;
	    g1sg = gamma_LO_1(pmom::index(nn,sign))   ;
	    g2sg = gamma_NLO_1(pmom::index(nn,sign))  ;
	    g3sg = gamma_NNLO_1(pmom::index(nn,sign)) ;
	    g4sg = gamma_NNNLO_1(pmom::index(nn,sign));
	    g1nsp = gammaNSp_LO_1(pmom::index(nn,sign));
	    g1nsm = gammaNSm_LO_1(pmom::index(nn,sign));
	    g1nsv = gammaNSv_LO_1(pmom::index(nn,sign));
	    g2nsp = gammaNSp_NLO_1(pmom::index(nn,sign));
	    g2nsm = gammaNSm_NLO_1(pmom::index(nn,sign));
	    g2nsv = gammaNSv_NLO_1(pmom::index(nn,sign));
	    g3nsp = gammaNSp_NNLO_1(pmom::index(nn,sign));
	    g3nsm = gammaNSm_NNLO_1(pmom::index(nn,sign));
	    g3nsv = gammaNSv_NNLO_1(pmom::index(nn,sign));
	    g4nsp = gammaNSp_NNNLO_1(pmom::index(nn,sign));
	    g4nsm = gammaNSm_NNNLO_1(pmom::index(nn,sign));
	    g4nsv = gammaNSv_NNNLO_1(pmom::index(nn,sign));
	  }
	else
	  {
	    nn = n-mellinint::mdim;
	    g1sg = gamma_LO_2(pmom::index(nn,sign))   ;
	    g2sg = gamma_NLO_2(pmom::index(nn,sign))  ;
	    g3sg = gamma_NNLO_2(pmom::index(nn,sign)) ;
	    g4sg = gamma_NNNLO_2(pmom::index(nn,sign));
	    g1nsp = gammaNSp_LO_2(pmom::index(nn,sign));
	    g1nsm = gammaNSm_LO_2(pmom::index(nn,sign));
	    g1nsv = gammaNSv_LO_2(pmom::index(nn,sign));
	    g2nsp = gammaNSp_NLO_2(pmom::index(nn,sign));
	    g2nsm = gammaNSm_NLO_2(pmom::index(nn,sign));
	    g2nsv = gammaNSv_NLO_2(pmom::index(nn,sign));
	    g3nsp = gammaNSp_NNLO_2(pmom::index(nn,sign));
	    g3nsm = gammaNSm_NNLO_2(pmom::index(nn,sign));
	    g3nsv = gammaNSv_NNLO_2(pmom::index(nn,sign));
	    g4nsp = gammaNSp_NNNLO_2(pmom::index(nn,sign));
	    g4nsm = gammaNSm_NNNLO_2(pmom::index(nn,sign));
	    g4nsv = gammaNSv_NNNLO_2(pmom::index(nn,sign));
	  }

      //Compute the Rk matrices with Eq. (2.21) of https://arxiv.org/pdf/hep-ph/0408244.pdf
      cmatrix R0 = g1sg/beta0;
      cmatrix R1 = g2sg/beta0 - b1*R0;
      cmatrix R2 = g3sg/beta0 - b1*R1 - b2*R0;
      cmatrix R3 = g4sg/beta0 - b1*R2 - b2*R1 - b3*R0;

      //Compute eigenvalues of the R0 matrix (See Eq. (2.27) of https://arxiv.org/pdf/hep-ph/0408244.pdf)
      complex<double> htr = R0.trace()/2.;
      complex<double> d = R0.det();
      rp[i] = htr + sqrt(htr*htr - d); 
      rm[i] = htr - sqrt(htr*htr - d);

      //Compute eigenvector projections (See Eq. (2.28) of https://arxiv.org/pdf/hep-ph/0408244.pdf)
      Ep[i] = 1./(rp[i]-rm[i])*(R0-rm[i]*I);
      Em[i] = 1./(rm[i]-rp[i])*(R0-rp[i]*I);

      //Compute Rktilde and Uk matrices (Eq. (2.25) and (2.31) of https://arxiv.org/pdf/hep-ph/0408244.pdf)
      cmatrix RTL1 = R1;
      U1[i] = -Em[i]*RTL1*Em[i]/1. - Ep[i]*RTL1*Ep[i]/1. + Ep[i]*RTL1*Em[i]/(rm[i]-rp[i]-1.) + Em[i]*RTL1*Ep[i]/(rp[i]-rm[i]-1.); //Eq. (2.31) of https://arxiv.org/pdf/hep-ph/0408244.pdf
      cmatrix RTL2 = R2+R1*U1[i];
      U2[i] = -Em[i]*RTL2*Em[i]/2. - Ep[i]*RTL2*Ep[i]/2. + Ep[i]*RTL2*Em[i]/(rm[i]-rp[i]-2.) + Em[i]*RTL2*Ep[i]/(rp[i]-rm[i]-2.); //Eq. (2.31) of https://arxiv.org/pdf/hep-ph/0408244.pdf
      cmatrix RTL3 = R3+R2*U1[i]+R1*U2[i];
      U3[i] = -Em[i]*RTL3*Em[i]/3. - Ep[i]*RTL3*Ep[i]/3. + Ep[i]*RTL3*Em[i]/(rm[i]-rp[i]-3.) + Em[i]*RTL3*Ep[i]/(rp[i]-rm[i]-3.); //Eq. (2.31) of https://arxiv.org/pdf/hep-ph/0408244.pdf

      //recursive relation
      if (opts.iterative)
	{
	  cmatrix R[maxord];
	  for (int k = 0; k < maxord; k++)
	    R[k] = N;
	      
	  if (opts.order_evol >= 1)
	    R[0] = R0;
	  if (opts.order_evol >= 2)
	    R[1] = R1;
	  if (opts.order_evol >= 3)
	    R[2] = R2;
	  if (opts.order_evol >= 4)
	    R[3] = R3;

	  if (opts.order_evol == 2)
	    for (int k = 2; k < maxord; k++)
	      R[k] += -b1*R[k-1];
	  else if (opts.order_evol == 3)
	    for (int k = 3; k < maxord; k++)
	      R[k] += -b1*R[k-1] -b2*R[k-2];
	  if (opts.order_evol >= 4)
	    for (int k = 4; k < maxord; k++)
	      R[k] += -b1*R[k-1] -b2*R[k-2] -b3*R[k-3];
	      
	  //Recursive relation to define RTL and U at all orders
	  cmatrix RTL[maxord];
	  for (int k = 1; k < maxord; k++)
	    {
	      RTL[k] = R[k];
	      for (int j = 1; j <= k-1; j++)
		RTL[k] += R[k-j]*U[idx(j,i)];
	      U[idx(k,i)] = -Em[i]*RTL[k]*Em[i]/double(k) - Ep[i]*RTL[k]*Ep[i]/double(k) + Ep[i]*RTL[k]*Em[i]/(rm[i]-rp[i]-double(k)) + Em[i]*RTL[k]*Ep[i]/(rp[i]-rm[i]-double(k)); //Eq. (2.31) of https://arxiv.org/pdf/hep-ph/0408244.pdf
	    }
	}

      //Non singlet (Eq. (2.21) of https://arxiv.org/pdf/hep-ph/0408244.pdf)
      complex <double> R0NSp = g1nsp/beta0;
      complex <double> R0NSm = g1nsm/beta0;
      complex <double> R0NSv = g1nsv/beta0;
      complex <double> R1NSp = g2nsp/beta0  - b1*R0NSp;
      complex <double> R1NSm = g2nsm/beta0  - b1*R0NSm;
      complex <double> R1NSv = g2nsv/beta0  - b1*R0NSv;
      complex <double> R2NSp = g3nsp/beta0  - b1*R1NSp - b2*R0NSp;
      complex <double> R2NSm = g3nsm/beta0  - b1*R1NSm - b2*R0NSm;
      complex <double> R2NSv = g3nsv/beta0  - b1*R1NSv - b2*R0NSv;
      complex <double> R3NSp = g4nsp/beta0  - b1*R2NSp - b2*R1NSp - b3*R0NSp;
      complex <double> R3NSm = g4nsm/beta0  - b1*R2NSm - b2*R1NSm - b3*R0NSm;
      complex <double> R3NSv = g4nsv/beta0  - b1*R2NSv - b2*R1NSv - b3*R0NSv;

      R0NS[i] = R0NSp;
	  
      U1NSp[i] = -R1NSp;
      U1NSm[i] = -R1NSm;
      U1NSv[i] = -R1NSv;

      complex <double> R2NSTLp = R2NSp+R1NSp*U1NSp[i];
      complex <double> R2NSTLm = R2NSm+R1NSm*U1NSm[i];
      complex <double> R2NSTLv = R2NSv+R1NSv*U1NSv[i];
      
      U2NSp[i] = -R2NSTLp/2.;
      U2NSm[i] = -R2NSTLm/2.;
      U2NSv[i] = -R2NSTLv/2.;

      complex <double> R3NSTLp = R3NSp+R2NSp*U1NSp[i]+R1NSp*U2NSp[i];
      complex <double> R3NSTLm = R3NSm+R2NSm*U1NSm[i]+R1NSm*U2NSm[i];
      complex <double> R3NSTLv = R3NSv+R2NSv*U1NSv[i]+R1NSv*U2NSv[i];

      U3NSp[i] = -R3NSTLp/3.;
      U3NSm[i] = -R3NSTLm/3.;
      U3NSv[i] = -R3NSTLv/3.;

      //recursive relation
      if (opts.iterative)
	{
	  complex <double> RNSp[maxord];
	  complex <double> RNSm[maxord];
	  complex <double> RNSv[maxord];
	  for (int k = 0; k < maxord; k++)
	    {
	      RNSp[k] = 0.;
	      RNSm[k] = 0.;
	      RNSv[k] = 0.;
	    }
	  if (opts.order_evol >= 1)
	    {
	      RNSp[0] = R0NSp;
	      RNSm[0] = R0NSm;
	      RNSv[0] = R0NSv;
	    }
	  if (opts.order_evol >= 2)
	    {
	      RNSp[1] = R1NSp;
	      RNSm[1] = R1NSm;
	      RNSv[1] = R1NSv;
	    }
	  if (opts.order_evol >= 3)
	    {
	      RNSp[2] = R2NSp;
	      RNSm[2] = R2NSm;
	      RNSv[2] = R2NSv;
	    }
	  if (opts.order_evol >= 4)
	    {
	      RNSp[3] = R3NSp;
	      RNSm[3] = R3NSm;
	      RNSv[3] = R3NSv;
	    }

	  if (opts.order_evol == 2)
	    for (int k = 2; k < maxord; k++)
	      {
		RNSp[k] += -b1*RNSp[k-1];
		RNSm[k] += -b1*RNSm[k-1];
		RNSv[k] += -b1*RNSv[k-1];
	      }
	  else if (opts.order_evol == 3)
	    for (int k = 3; k < maxord; k++)
	      {
		RNSp[k] += -b1*RNSp[k-1]-b2*RNSp[k-2];
		RNSm[k] += -b1*RNSm[k-1]-b2*RNSm[k-2];
		RNSv[k] += -b1*RNSv[k-1]-b2*RNSv[k-2];
	      }
	  else if (opts.order_evol >= 4)
	    for (int k = 4; k < maxord; k++)
	      {
		RNSp[k] += -b1*RNSp[k-1]-b2*RNSp[k-2]-b3*RNSp[k-3];
		RNSm[k] += -b1*RNSm[k-1]-b2*RNSm[k-2]-b3*RNSm[k-3];
		RNSv[k] += -b1*RNSv[k-1]-b2*RNSv[k-2]-b3*RNSv[k-3];
	      }
	      
	  //Recursive relation to define RTL and U at all orders
	  complex <double> RTLNSp[maxord];
	  complex <double> RTLNSm[maxord];
	  complex <double> RTLNSv[maxord];
	  for (int k = 1; k < maxord; k++)
	    {
	      RTLNSp[k] = RNSp[k];
	      RTLNSm[k] = RNSm[k];
	      RTLNSv[k] = RNSv[k];
	      for (int j = 1; j <= k-1; j++)
		{
		  RTLNSp[k] += RNSp[k-j]*UNSp[idx(j,i)];
		  RTLNSm[k] += RNSm[k-j]*UNSm[idx(j,i)];
		  RTLNSv[k] += RNSv[k-j]*UNSv[idx(j,i)];
		}
	      UNSp[idx(k,i)] = -RTLNSp[k]/double(k); //Eq. (2.31) of https://arxiv.org/pdf/hep-ph/0408244.pdf
	      UNSm[idx(k,i)] = -RTLNSm[k]/double(k); //Eq. (2.31) of https://arxiv.org/pdf/hep-ph/0408244.pdf
	      UNSv[idx(k,i)] = -RTLNSv[k]/double(k); //Eq. (2.31) of https://arxiv.org/pdf/hep-ph/0408244.pdf
	    }
	  //cout << U3NSp[i] << "  " << UNSp[idx(3,i)] << endl;
	  //cout << R3NSTLp << "  " << RTLNSp[3] << endl;
	}
    }
}

void evoltrunc::evolve()
{
  //Evolve from Q/muF to mub
  retrieve_pdf_muf();

  //Compute alphas for the evolution
  double asi;
  //if (opts.commute)
  //  asi = resint::alpqfac*4.*M_PI;
  //else
  //  asi = resint::alpqres*4.*M_PI;
  //Force msbar scheme for the PDF evolution (!!! but need the pmom::gamma in msbar)
  if (opts.commute)
    asi = pdf::alphas(scales::fac);
  else
    asi = pdf::alphas(scales::res);
  complex <double> asf = pdfevol::asl*asi;

  //cout << "asi " << asi << " asf " << asf << " asl " << pdfevol::asl << endl;
  
  bool hfswitch = false;
  if (opts.hfmode == 1)
    hfswitch = true;

  evolve(asi,asf, opts.order_evol, hfswitch, 0, opts.smallx);
  
  if (opts.hfmode == 2)
    {
      complex <double> asf_b = pdfevol::asl_b*asi;
      evolve(asi,asf_b, opts.order_evol, false, 1, opts.smallx);
      complex <double> asf_c = pdfevol::asl_c*asi;
      evolve(asi,asf_c, opts.order_evol, false, 2, opts.smallx);
    }
  else if (opts.hfmode == 3)
    {
      complex <double> asf_b = pdfevol::asl_b*asi;
      evolve(asi,asf_b, opts.order_evol, false, 1, false);
      complex <double> asf_c = pdfevol::asl_c*asi;
      evolve(asi,asf_c, opts.order_evol, false, 2, false);
    }
}

void evoltrunc::commute()
{
  if (!opts.commute)
    return;
  if (opts.kmufac == opts.kmures)
    return;

  //Evolve from Q to muF
  retrieve_pdf_mub();

  //should this piece be done in MSBAR or in CMW? --> the gamma anomalous dimension in the LQF*gamma terms in the H coefficients, the C coefficients and the beta function are in CMW
  //--> could do this piece in CMW to fully match the Hst in CMW
  double beta0 = alphas::beta0_msbar;
  double beta1 = alphas::beta1_msbar;
  double beta2 = alphas::beta2_msbar;
  double beta3 = alphas::beta3_msbar;
  //double beta0 = resconst::beta0;
  //double beta1 = resconst::beta1;
  //double beta2 = resconst::beta2;
  //double beta3 = resconst::beta3;
  
  //compute asi and asf
  double asi;
  complex <double> asf;
  double blog = log(pow(scales::res/scales::fac,2));
  double as = resint::aass;
  double LQR = resint::LR-resint::LQ;
  double as2 = pow(as,2);
  double as3 = pow(as,3);
  double xlambda = beta0*as*blog;
  double log1xlambda = log(1.-xlambda);
  double logas = 0.;

  //The order of this piece should match the order of the gamma terms in the Hst coefficients
  //int order = opts.order_evol;
  int order = opts.order_hcoef;
  if (order >= 1)
    logas += log1xlambda;
  if (order >= 2)
    logas += as* beta1/beta0*log1xlambda/(1.-xlambda);
  if (order >= 3)
    logas += as2* ((pow(beta1/beta0,2)-beta2/beta0) *xlambda/pow(1.-xlambda,2)
		   + pow(beta1/beta0,2)             *log1xlambda/pow(1.-xlambda,2)
		   - pow(beta1/beta0,2)             *pow(log1xlambda,2)/(2.*pow(1.-xlambda,2)));
  if (order >= 4)
    logas += as3*
      (3.*xlambda*(pow(beta0,2)*beta3*(-2. + xlambda)
		   - 2.*beta0*beta1*beta2*(-1. + xlambda)
		   + pow(beta1,3)*xlambda)
       + 6.*(-(pow(beta1,3)*xlambda)
	     + beta0*beta1*beta2*(1. + xlambda))*log1xlambda
       - 9.*pow(beta1,3)*pow(log1xlambda,2)
       + 2.*pow(beta1,3)*pow(log1xlambda,3)
       )/(6.*pow(beta0,3)*pow(1. - xlambda,3));
  if (order >= 2)
    logas += as*(LQR)
      *beta0*xlambda/(1.-xlambda);
  if (order >= 3)
    logas += as2*(+LQR*beta1                   *(xlambda-log1xlambda)/pow(1.-xlambda,2)
		  +LQR*beta1                   *xlambda/(1.-xlambda)                      //missing piece
		  +0.5*pow(LQR,2)*pow(beta0,2) *xlambda*(xlambda-2.)/pow(1.-xlambda,2));  //missing piece
  if (order >= 4)
    logas += as3*
      (LQR*(6.*pow(beta1,2)*pow(log1xlambda,2)
	    + 6.*beta1*log1xlambda*(-3.*beta1 + pow(beta0,2)*LQR)
	    + xlambda*(-6.*pow(beta1,2) + 6.*beta0*beta2*(4. - 3.*xlambda + pow(xlambda,2))
		       - 9.*pow(beta0,2)*beta1*(3. - 3.*xlambda + pow(xlambda,2))*LQR
		       + 2.*pow(beta0,4)*(3. - 3.*xlambda + pow(xlambda,2))*pow(LQR,2))))/(6.*beta0*pow(1. - xlambda,3));
  asi = resint::alpqres*4.*M_PI;
  asf = exp(-logas)*asi;

  evolve(asi,asf, order);
}

//Computed truncated solution of evolution equations
void evoltrunc::evolve(double asi, complex <double> asf, int order, bool hfswitch, int hfl, bool smallx)
{
  //Mellin arrays of evolution matrices
  cmatrix usg[dim*2];
  complex <double> unsp[dim*2];
  complex <double> unsm[dim*2];
  complex <double> unsv[dim*2];
  fill(usg,usg+dim*2, I);
  fill(unsp,unsp+dim*2, 1.);
  fill(unsm,unsm+dim*2, 1.);
  fill(unsv,unsv+dim*2, 1.);

  complex <double> asl = -log(asf/asi);

  //Singlet
  if (smallx)
    {
#ifdef USEHELLX
      hell::calcusg(asi, asf, order);
      copy(hell::usg, hell::usg+2*dim, usg);
#endif
    }
  else
    for (int i = 0; i < dim*2; i++)
      {
	//LO evolution matrix
	cmatrix ULO = Em[i]*exp(asl*rm[i]) + Ep[i]*exp(asl*rp[i]); //Eq. (2.29) of https://arxiv.org/pdf/hep-ph/0408244.pdf
      
	//The all-order iterative U and U^-1 matrices
	if (opts.iterative)
	  {
	    cmatrix UF = I;
	    cmatrix UI = I;
	    for (int k = 1; k < maxord; k++)
	      {
		UF += pow(asf,k)*U[idx(k,i)];
		UI += pow(asi,k)*U[idx(k,i)];
	      }
	  
	    //The full inverse UIm1 of UI
	    cmatrix UIm1 = UI.inverse();
	    if (opts.expc == 10)
	      usg[i] = ULO*UIm1;
	    else
	      usg[i] = UF*ULO*UIm1;
	  }
	else
	  {
	    //Eq. (2.24) of https://arxiv.org/pdf/hep-ph/0408244.pdf
	    //LO evolution
	    if (order >= 1)
	      usg[i] = ULO;
	      
	    //NLO evolution
	    if (order >= 2)
	      usg[i] += asf*U1[i]*ULO - asi*ULO*U1[i];
	      
	    //NNLO evolution
	    if (order >= 3)
	      usg[i] += asf*asf*U2[i]*ULO - asf*asi*U1[i]*ULO*U1[i] + asi*asi*ULO*(U1[i]*U1[i]-U2[i]);
	      
	    //NNNLO evolution
	    if (order >= 4)
	      usg[i] += asf*asf*asf*U3[i]*ULO - asf*asf*asi*U2[i]*ULO*U1[i] + asf*asi*asi*U1[i]*ULO*(U1[i]*U1[i]-U2[i])
		- asi*asi*asi*ULO*(U1[i]*U1[i]*U1[i]-U1[i]*U2[i]-U2[i]*U1[i]+U3[i]); //typo in Eq. (2.24) of https://arxiv.org/pdf/hep-ph/0408244.pdf, see Eq. (34) of https://doi.org/10.1016/0550-3213(94)90076-0
	  }
      }

  //Non singlet
  for (int i = 0; i < dim*2; i++)
    {
      complex<double> UNSLO = exp(R0NS[i]*asl);

      //The all-order iterative U and U^-1 matrices
      if (opts.iterative)
	{
	  complex <double> UNSpF = 1.;
	  complex <double> UNSmF = 1.;
	  complex <double> UNSvF = 1.;
	  complex <double> UNSpI = 1.;
	  complex <double> UNSmI = 1.;
	  complex <double> UNSvI = 1.;
	  for (int k = 1; k < maxord; k++)
	    {
	      UNSpF += pow(asf,k)*UNSp[idx(k,i)];
	      UNSmF += pow(asf,k)*UNSm[idx(k,i)];
	      UNSvF += pow(asf,k)*UNSv[idx(k,i)];
	      UNSpI += pow(asi,k)*UNSp[idx(k,i)];
	      UNSmI += pow(asi,k)*UNSm[idx(k,i)];
	      UNSvI += pow(asi,k)*UNSv[idx(k,i)];
	    }
	  if (opts.expc == 10)
	    {
	      unsp[i] = 1./UNSpI*UNSLO;
	      unsm[i] = 1./UNSmI*UNSLO;
	      unsv[i] = 1./UNSvI*UNSLO;
	    }
	  else
	    {
	      unsp[i] = UNSpF/UNSpI*UNSLO;
	      unsm[i] = UNSmF/UNSmI*UNSLO;
	      unsv[i] = UNSvF/UNSvI*UNSLO;
	    }
	}
      else
	{
	  if (order >= 1)
	    {
	      unsp[i] = UNSLO;
	      unsm[i] = UNSLO;
	      unsv[i] = UNSLO;
	    }
      
	  //NLO evolution
	  if (order >= 2)
	    {
	      unsp[i] += (asf-asi)*U1NSp[i]*UNSLO;
	      unsm[i] += (asf-asi)*U1NSm[i]*UNSLO;
	      unsv[i] += (asf-asi)*U1NSv[i]*UNSLO;
	    }
	      
	  //NNLO evolution
	  if (order >= 3)
	    {
	      unsp[i] += asf*asf*U2NSp[i]*UNSLO - asf*asi*U1NSp[i]*U1NSp[i]*UNSLO + asi*asi*(U1NSp[i]*U1NSp[i]-U2NSp[i])*UNSLO;
	      unsm[i] += asf*asf*U2NSm[i]*UNSLO - asf*asi*U1NSm[i]*U1NSm[i]*UNSLO + asi*asi*(U1NSm[i]*U1NSm[i]-U2NSm[i])*UNSLO;
	      unsv[i] += asf*asf*U2NSv[i]*UNSLO - asf*asi*U1NSv[i]*U1NSv[i]*UNSLO + asi*asi*(U1NSv[i]*U1NSv[i]-U2NSv[i])*UNSLO;
	    }

	  //NNNLO evolution
	  if (order >= 4)
	    {
	      unsp[i] += asf*asf*asf*U3NSp[i]*UNSLO - asf*asf*asi*U2NSp[i]*U1NSp[i]*UNSLO + asf*asi*asi*U1NSp[i]*UNSLO*(U1NSp[i]*U1NSp[i]-U2NSp[i]) - asi*asi*asi*UNSLO*(U1NSp[i]*U1NSp[i]*U1NSp[i]-2.*U1NSp[i]*U2NSp[i]+U3NSp[i]);
	      unsm[i] += asf*asf*asf*U3NSm[i]*UNSLO - asf*asf*asi*U2NSm[i]*U1NSm[i]*UNSLO + asf*asi*asi*U1NSm[i]*UNSLO*(U1NSm[i]*U1NSm[i]-U2NSm[i]) - asi*asi*asi*UNSLO*(U1NSm[i]*U1NSm[i]*U1NSm[i]-2.*U1NSm[i]*U2NSm[i]+U3NSm[i]);
	      unsv[i] += asf*asf*asf*U3NSv[i]*UNSLO - asf*asf*asi*U2NSv[i]*U1NSv[i]*UNSLO + asf*asi*asi*U1NSv[i]*UNSLO*(U1NSv[i]*U1NSv[i]-U2NSv[i]) - asi*asi*asi*UNSLO*(U1NSv[i]*U1NSv[i]*U1NSv[i]-2.*U1NSv[i]*U2NSv[i]+U3NSv[i]);
	    }
	}
      //cout << "analytic usg " << usg[0] << endl;
      //cout << "analytic unsp " << unsp[0] << endl;
      //cout << "analytic unsm " << unsm[0] << endl;
      //cout << "analytic unsv " << unsv[0] << endl;
    }

  //QED evolution before QCD evolution
  if (opts.qed && opts.qedorder_evol >= 1 && pdf::photon && hfl == 0)
    {
      int n,nn,sign;
      for (int i = 0; i < dim*2; i++)
	{
	  n    = (i<dim)? i : i-dim;
	  sign = (i<dim)? mesq::positive : mesq::negative;

	  //Flavour decomposition of the minus sector:  valence distributions
	  complex <double> bv = (vai[i] - m24i[i]) * 0.2;
	  complex <double> cv = bv + (m24i[i] - m15i[i]) * 0.25;
	  complex <double> sv = cv + (m15i[i] - m8i[i]) /3.;
	  complex <double> dv = sv + (m8i[i] - m3i[i]) * 0.5;
	  complex <double> uv = sv + (m8i[i] + m3i[i]) * 0.5;
	  //Flavour decomposition of the plus sector: antiquark distributions
	  complex <double> bp = (sgi[i] - p24i[i]) * 0.2;
	  complex <double> cp = bp + (p24i[i] - p15i[i]) * 0.25;
	  complex <double> sp = cp + (p15i[i] - p8i[i]) /3.;
	  complex <double> dp = sp + (p8i[i] - p3i[i]) * 0.5;
	  complex <double> up = sp + (p8i[i] + p3i[i]) * 0.5;
	  //sea distributions
	  complex <double> bs = 0.5*(bp - bv);
	  complex <double> cs = 0.5*(cp - cv);
	  complex <double> ss = 0.5*(sp - sv);
	  complex <double> ds = 0.5*(dp - dv);
	  complex <double> us = 0.5*(up - uv);
	  complex <double> ft = fti[i];
	  //QED non singlet
	  complex <double> NST1up = uv + 2.*us - cv - 2.*cs;
	  complex <double> NST1dn = dv + 2.*ds - sv - 2.*ss;
	  complex <double> NST3dn = dv + 2.*ds + sv + 2.*ss - 2.*bv - 4.*bs;
	  complex <double> SINup  = uv + 2.*us + cv + 2.*cs;
	  complex <double> SINdn  = dv + 2.*ds + sv + 2.*ss + bv + 2.*bs;

	  complex <double> Sqed   = log(1./pdfevol::aqedl);

	  double B0QED  =-4.*alphaqed::NF2/3.;
	  complex <double> QQqed;//=(3./4.)*(-2.*pmom::gamma1qq[i]); //Remove CF, add squared charge

	  if (n < mellinint::mdim)
	    nn = n;
	  else
	    nn = n-mellinint::mdim;
	  if (opts.mellin1d)
	    QQqed=(3./4.)*(-2.*pmom::gamma1qq[pmom::index(nn,sign)]); //Remove CF, add squared charge
	  else
	    if (n < mellinint::mdim)
	      QQqed=(3./4.)*(-2.*pmom::gamma1qq_1[pmom::index(nn,sign)]); //Remove CF, add squared charge
	    else
	      QQqed=(3./4.)*(-2.*pmom::gamma1qq_2[pmom::index(nn,sign)]); //Remove CF, add squared charge

	  //complex <double> QAqed= 2.*3.* QGF;                   //Es proporcional a la carga al cuadrado, incluye NC            
	  //complex <double> AQqed= (3./4.)*GQI;                  //Removemos el factor CF, es proporcional a la carga al cuadrado
	  complex <double> AAqed= 2. * (4./3.) * alphaqed::NF2;   //Only piece with summmed charges
	  complex <double> ansqedup = alphaqed::equarkup*4.*QQqed / (2.*B0QED);
	  complex <double> ansqeddn = alphaqed::equarkdn*4.*QQqed / (2.*B0QED);
	  complex <double> ansqedft = AAqed / (2.*B0QED);
	  
	  complex <double> ENSQEDup = exp(-ansqedup*Sqed);
	  complex <double> ENSQEDdn = exp(-ansqeddn*Sqed);
	  complex <double> ENSQEDft = exp(-ansqedft*Sqed);
      
	  //Evolve
	  uv     *= ENSQEDup;
	  dv     *= ENSQEDdn;
	  sv     *= ENSQEDdn;
	  cv     *= ENSQEDdn;
	  bv     *= ENSQEDdn;
	  NST1up *= ENSQEDup;
	  NST1dn *= ENSQEDdn;
	  NST3dn *= ENSQEDdn;
	  SINup  *= ENSQEDup;
	  SINdn  *= ENSQEDdn;
	  ft     *= ENSQEDft;
       
	  //Back to QCD evolution basis
	  complex <double> NS3N  = (-3.*NST1dn + 3.*NST1up - NST3dn - 2.*SINdn + 3.*SINup)/6.;
	  complex <double> NS8N  = ( 9.*NST1dn + 3.*NST1up - NST3dn - 2.*SINdn + 3.*SINup)/6. ;
	  complex <double> NS15N = ( 6.*NST1up + NST3dn + 2.*SINdn - 3.*SINup)/3.;
	  complex <double> NS24N = ( 5.*NST3dn)/3. - (2.*SINdn)/3. + SINup;
	  complex <double> SIN   = SINup+SINdn;
	  complex <double> NS35N = SIN;

	  ss = (10.* SIN + 2.* NS35N + 3.* NS24N + 5.* NS15N - 20.* NS8N - 60.* sv) / 120.;
	  ds = (10.* SIN + 2.* NS35N + 3.* NS24N + 5.* NS15N + 10.* NS8N - 30.* NS3N - 60.* dv) / 120.;
	  us = (10.* SIN + 2.* NS35N + 3.* NS24N + 5.* NS15N + 10.* NS8N + 30.* NS3N - 60.* uv) / 120.;
	  cs = (10.* SIN + 2. *NS35N + 3.* NS24N - 15.* NS15N - 60.* cv) / 120.;
	  bs = (10.* SIN + 2. *NS35N - 12.* NS24N - 60.* bv) / 120.;
	  
	  complex <double> fx[11];
	  fx[u ] = uv + us;
	  fx[d ] = dv + ds;
	  fx[s ] = sv + ss;
	  fx[c ] = cv + cs;
	  fx[b ] = bv + bs;
	  fx[g ] = gli[i];
	  fx[ub] = us;
	  fx[db] = ds;
	  fx[sb] = ss;    
	  fx[cb] = cs;	     
	  fx[bb] = bs;	     

	  if (opts.mellin1d)
	    pdfevol::storemoments(n, sign, fx);
	  else
	  if (n < mellinint::mdim)
	    pdfevol::storemoments_1(n, sign, fx);
	  else
	    pdfevol::storemoments_2(n-mellinint::mdim, sign, fx);
	  
	}
      //retrieve evolved PDFs into the input PDFs for the QCD evolution
      retrieve_pdf_mub();
    }
  
  //Loop on Mellin moments
  //cout << endl;
  //cout << "alphas " << asi << "  " << asf << endl;
  int n,nn,sign;
  for (int i = 0; i < dim*2; i++)
    {
      n    = (i<dim)? i : i-dim;
      sign = (i<dim)? mesq::positive : mesq::negative;

      //PDFs at the starting scale
      dytmatrix::vec2<complex<double>> f0(gli[i],sgi[i]); //(gluon,singlet)

      //Evolved PDFs

      //singlet
      dytmatrix::vec2<complex<double>> f1 = usg[i]*f0;
      complex <double> glf = f1.entry1();
      complex <double> sgf = f1.entry2();
      
      //non singlet
      complex <double> vaf  = unsv[i]*vai[i];
      complex <double> m3f  = unsm[i]*m3i[i];
      complex <double> m8f  = unsm[i]*m8i[i];
      complex <double> m15f = unsm[i]*m15i[i];
      complex <double> m24f = unsm[i]*m24i[i];
      complex <double> p3f  = unsp[i]*p3i[i];
      complex <double> p8f  = unsp[i]*p8i[i];
      complex <double> p15f = unsp[i]*p15i[i];
      complex <double> p24f = unsp[i]*p24i[i];
      
      //Flavour decomposition of the minus sector:  valence distributions
      complex <double> bv = (vaf - m24f) * 0.2;
      complex <double> cv = bv + (m24f - m15f) * 0.25;
      complex <double> sv = cv + (m15f - m8f) /3.;
      complex <double> dv = sv + (m8f - m3f) * 0.5;
      complex <double> uv = sv + (m8f + m3f) * 0.5;
      //Flavour decomposition of the plus sector: antiquark distributions
      complex <double> bp = (sgf - p24f) * 0.2;
      complex <double> cp = bp + (p24f - p15f) * 0.25;
      complex <double> sp = cp + (p15f - p8f) /3.;
      complex <double> dp = sp + (p8f - p3f) * 0.5;
      complex <double> up = sp + (p8f + p3f) * 0.5;

      complex <double> fx[11];
      fx[u ] = 0.5 * (uv + up);
      fx[d ] = 0.5 * (dv + dp);
      fx[s ] = 0.5 * (sv + sp);
      fx[c ] = 0.5 * (cv + cp);
      fx[b ] = 0.5 * (bv + bp);
      fx[g ] = glf;
      fx[ub] = 0.5 * (up - uv); 
      fx[db] = 0.5 * (dp - dv);
      fx[sb] = 0.5 * (sp - sv);
      fx[cb] = 0.5 * (cp - cv);
      fx[bb] = 0.5 * (bp - bv);

      //switching of HF pdfs
      if (hfswitch)
	{
	  double mub = real(pdfevol::mubstartilde);
	  //double mub = real(pdfevol::mub);
	  double delta = 1.;
	  double k = 2.;
	  if (mub < k*pdf::mb)
	    {
	      double swtch = exp(-pow((k*pdf::mb - mub),2)/pow((pdf::mb*delta),2));
	      fx[b ] *= swtch;
	      fx[bb] *= swtch;
	    }
	  if (mub < k*pdf::mc)
	    {
	      double swtch = exp(-pow((k*pdf::mc - mub),2)/pow((pdf::mc*delta),2));
	      fx[c ] *= swtch;
	      fx[cb] *= swtch;
	    }
	}
      
      //cout << "evoltrunc " << asf << "  "
      //	   << i << "  " << cxn[i] << "  "
      //	   << fx[g] << "  " << uv << "  " << dv << "  " << fx[ub] << "  " << fx[db] << "  " << fx[s] << "  " << fx[c] << "  " << fx[b]  << endl;

      if (hfl == 0)
	{
	  if (opts.mellin1d)
	    pdfevol::storemoments(n, sign, fx);
	  else
	  if (n < mellinint::mdim)
	    pdfevol::storemoments_1(n, sign, fx);
	  else
	    pdfevol::storemoments_2(n-mellinint::mdim, sign, fx);
	}
      else if (hfl == 1)
	{
	  if (opts.mellin1d)
	    pdfevol::storemoments_b(n, sign, fx);
	  else
	  if (n < mellinint::mdim)
	    pdfevol::storemoments_b(n, sign, fx);
	  else
	    pdfevol::storemoments_b(n-mellinint::mdim, sign, fx);
	}
      else if (hfl == 2)
	{
	  if (opts.mellin1d)
	    pdfevol::storemoments_c(n, sign, fx);
	  else
	  if (n < mellinint::mdim)
	    pdfevol::storemoments_c(n, sign, fx);
	  else
	    pdfevol::storemoments_c(n-mellinint::mdim, sign, fx);
	}

    }
    

  //evolnum::calculate();
  //for (int i = 0; i < mellinint::mdim; i++)
  //  {
  //    //int i = 0;
  //    //cout << "direct " << fabs(pdfevol::bstarscale) << "  "
  //    cout << "direct " << real(pdfevol::mubstartilde) << "  "
  //	   << i << "  " << mellinint::Np[i] << "  "
  //	   << pdfevol::fx1[i*11+(0+MAXNF )] << "  "
  //	   << pdfevol::fx1[i*11+(1+MAXNF )]-pdfevol::fx1[i*11+(-1+MAXNF)] << "  "
  //	   << pdfevol::fx1[i*11+(2+MAXNF )]-pdfevol::fx1[i*11+(-2+MAXNF)] << "  "
  //	   << pdfevol::fx1[i*11+(-1+MAXNF)] << "  "
  //	   << pdfevol::fx1[i*11+(-2+MAXNF)] << "  "
  //	   << pdfevol::fx1[i*11+(3+MAXNF )] << "  "
  //	   << pdfevol::fx1[i*11+(4+MAXNF )] << "  "
  //	   << pdfevol::fx1[i*11+(5+MAXNF )] << endl;
  //  }
  
  return;
}
