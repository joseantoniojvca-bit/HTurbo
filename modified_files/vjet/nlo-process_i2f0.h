//  Copyright (C) 2002 Zoltan Nagy
//                2012 Leandro Cieri
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#ifndef __NLO_PROCESS_I2F0_H__
#define __NLO_PROCESS_I2F0_H__ 1

//   nlojet++ includes
#include "nlo-dipole_mom.h"
#include "nlo-process.h"
//#include "H1.h"
#include "mstwpdf.h"
#include <math.h>
#include "Inter_params.h"

// DYTurbo includes
#include "coupling.h"
#include "parton.h"
#include "settings.h"
#include "../dyres/dyres_interface.h"
#include "pdf.h"
#include "interface.h"
#include "scales.h"

/* extern "C" { */
/*     double alphass_(double * mu2,int *Iset);    } */


/*  inline double ALPHASS(double mu2,int Iset)  { */
/*   return alphass_(  &mu2,&Iset); */
/*  } */

weight_hhc2ph GeT_pdf(double, double, double, unsigned int, unsigned int);
weight_hhc2ph GeT_pdfCT(double, double, double, unsigned int, unsigned int);
weight_hhc2ph GeT_pdfCT(double, double, double, unsigned int, unsigned int,double (&)[13] ,double (&)[13]);


static const double twopi = 6.28318530717958647692;
//double __1alfaQED2 =  1.0/alpha_qed/ alpha_qed ; 

  //
  //  Declaraton of the abstract class process in the
  //  case of hadron-hadron collision.
  template<typename _Weight, class _Event, class _EvenTraits>
  class process<_Weight, _Event, _EvenTraits, 2U, 0U>
  {
  public:
    //   types
    typedef _Weight weight_type;
    typedef _Event event_type;

    //   member access
    unsigned int npar() const { return _M_npar;}
    unsigned int npow() const { return _M_npow;}
    double alpha() const { return _M_alpha;}

    //  set the alpha parameter
    void alpha(double __new_alpha) { _M_alpha = __new_alpha;}

    //  number of quark flavours
    unsigned int nu() const { return _M_nu;}
    unsigned int nd() const { return _M_nd;}
    unsigned int nf() const { return _M_nu + _M_nd;}

    //   destructor
    virtual ~process() {}

        //   born contributions
    virtual void LO_term(const _Event&, _Weight&) = 0;

       // BOX term
    virtual void BOX_term(const _Event&, _Weight&) = 0;

//---------------------------------------------------------------------
//---------------------------------------------------------------------
//---------------------------------------------------------------------
//------- Interferencia Reales ---- Neri ------------------------------
//---------------------------------------------------------------------
//---------------------------------------------------------------------
//---------------------------------------------------------------------

    // ggffg term
    virtual void ggffg_term(const _Event&, _Weight&) = 0;
    // ggHg term
    virtual void ggHg_term(const _Event&, _Weight&) = 0;
    // ggINTg term
    virtual void ggINTg_term(const _Event&, _Weight&) = 0;
    // qgINTg term
    virtual void qgINTg_term(const _Event&, _Weight&) = 0;

//---------------------------------------------------------------------
//---------------------------------------------------------------------
//---------------------------------------------------------------------
// Interferencia LO ---------------- LEAN -----------------------------
//---------------------------------------------------------------------
//---------------------------------------------------------------------
//---------------------------------------------------------------------    

    // Señal
    virtual void MggHyy_term(const _Event&, _Weight&) = 0;
    virtual void MqqbarHyy_term(const _Event&, _Weight&) = 0;
    // CTgg Term
    virtual void CTgg_term(const _Event&, _Weight&,double,double,unsigned int,double,double,int) = 0;

    
//---------------------------------------------------------------------
//---------------------------------------------------------------------
//---------------------------------------------------------------------
//---------------------------------------------------------------------
//---------------------------------------------------------------------
//---------------------------------------------------------------------


         // H1f Term
    virtual void H1f_term(const _Event&, _Weight&) = 0;

    // H1 Term
    virtual void H1_term(const _Event&, _Weight&,double,double,unsigned int,double,double,int) = 0;

    // CTNLO Term
    virtual void CTNLO_term(const _Event&, _Weight&,double,double,unsigned int,double,double,int) = 0;

       // H2f Term
    virtual void H2f_term(const _Event&, _Weight&,double,unsigned int) = 0;

    //   born contributions
    virtual void born_term(const _Event&, _Weight&) = 0;

    //   real contributions
    virtual void real_term(const _Event&, _Weight&) = 0;

    //  finite contributions (1-loop,...)
    virtual void fini_term(double, double, double, double, const _Event&, _Weight *) = 0;

    //  finite contributions (1-loop,...)
    virtual void fini_term(double, double, double, double, const _Event&, _Weight *,double,double) = 0;

    //   dipole contributions
    virtual void dipole_term(const _Event&, const _Event&, int, int, int, _Weight&) = 0;

    //   select the dipoles (extra contraints for the dipole indices
    //     e.g: in photoproduction k > -1 no photon spectator)
    virtual bool dipole_index(int, int, int) {
      return true;
    }

    //   generate the dipole momenta
    bool dipole_mom(const _Event&, int, int, int, _Event&);

  protected:
    //  constructors
    explicit process(unsigned int np, unsigned int nw, unsigned int nu,
		     unsigned int nd, double al = 1.0)
      : _M_npar(np), _M_npow(nw), _M_nu(nu), _M_nd(nd), _M_alpha(al) {}

  private:
    //   data members
    unsigned int _M_npar, _M_npow, _M_nu, _M_nd;
    double _M_alpha;
  };

  template<typename _Weight, class _Event, class _EvenTraits>
  bool process<_Weight, _Event, _EvenTraits, 2U, 0U>::
  dipole_mom(const _Event& p, int i, int j, int k, _Event& q)
  {
    q[hadron(-1)] = p[hadron(-1)];
    q[hadron( 0)] = p[hadron( 0)];
    if(i > 0 && j > 0 && k > 0)
      return dipole_mom_fff<_Event>(_M_alpha, p, i, j, k, q);
    else if(i > 0 && j > 0 && k <= 0)
      return dipole_mom_ffi<_Event>(_M_alpha, p, i, j, k, q);
    else if(i <= 0 && j > 0 && k > 0)
      return dipole_mom_iff<_Event>(_M_alpha, p, i, j, k, q);
    else if(i <= 0 && j > 0 && k <= 0)
      return dipole_mom_ifi<_Event>(_M_alpha, p, i, j, k, q);
    else return false;
  }

//  template<class _Weight>
//  struct pdf_and_coupling<_Weight, 2U, 0U>
//  {
    //   Destructor
//    virtual ~pdf_and_coupling() {}

    //   The QCD coupling
 //   virtual double alpha_qcd(unsigned int, double) = 0;

    //   the parton distribution function
//    virtual _Weight pdf(double, double, double, unsigned int=2U, unsigned int=3U) = 0;
//  };

  template<typename _Weight, class _Event, class _EvenTraits>
  class amplitude<process<_Weight, _Event, _EvenTraits, 2U, 0U> >
  {
  public:
    //   public types
    typedef process<_Weight, _Event, _EvenTraits, 2U, 0U> process_type;
    typedef typename process_type::event_type event_type;
    typedef typename process_type::weight_type weight_type;
  //  typedef pdf_and_coupling<_Weight, 2U, 0U> pdf_type;

    //   type of the contributions
    enum contrib_type { notdef = -1, born = 0, real, sub, fini,LO,H1f,BOX,H2f,H1,CTNLO,ggffg,ggHg,ggINTg,qgINTg,LOIntHyy,CTgg,LOIntHyyqqbar};
    enum integral_type { lo = 0, nlo};

    //   constructor
    explicit amplitude(process_type *proc, integral_type itype)
      : _M_proc(proc), _M_itype(itype) {}

    //   leading  order contribution LO POSTA
    void LO_contrib(double w, const _Event& p) {
      _M_p = &p; _M_weight = w; _M_calc = false; _M_contr = LO;
    }

    //   leading  H1f order contribution NLO POSTA
    void H1f_contrib(double w, const _Event& p) {
      _M_p = &p; _M_weight = w; _M_calc = false; _M_contr = H1f;
    }

    //   leading  H1 order contribution NLO POSTA
    void H1_contrib(double w, const _Event& p,double z1, double z2,int orderH) {
      _M_p = &p; _M_weight = w; _M_calc = false; _M_contr = H1;
      _M_H1.order = orderH;
      _M_H1.z1 = z1; _M_H1.z2=z2;
    }

    //   leading  CTNLO order contribution NLO POSTA
    void CTNLO_contrib(double w, const _Event& p,double z1, double z2,int orderct) {
      _M_p = &p; _M_weight = w; _M_calc = false; _M_contr = CTNLO;
       _M_CTNLO.order = orderct;
      _M_CTNLO.z1 = z1; _M_CTNLO.z2=z2;
    }

    //   leading  order contribution LO POSTA
    void H2f_contrib(double w, const _Event& p) {
      _M_p = &p; _M_weight = w; _M_calc = false; _M_contr = H2f;
    }

    //   leading  order contribution BOX POSTA
    void BOX_contrib(double w, const _Event& p) {
      _M_p = &p; _M_weight = w; _M_calc = false; _M_contr = BOX;
    }
//---------------------------------------------------------------------
//---------------------------------------------------------------------
//---------------------------------------------------------------------
//---------------------------------------------------------------------
//---------------------------------------------------------------------
//---------------------------------------------------------------------
//---------------------------------------------------------------------
//---------------------------------------------------------------------
//---------------------------------------------------------------------


    //   leading  order contribution ggffg POSTA
    void ggffg_contrib(double w, const _Event& p) {
      _M_p = &p; _M_weight = w; _M_calc = false; _M_contr = ggffg;
    }
    //   leading  order contribution ggHg POSTA
    void ggHg_contrib(double w, const _Event& p) {
      _M_p = &p; _M_weight = w; _M_calc = false; _M_contr = ggHg;
    }
    //   leading  order contribution ggINTg POSTA
    void ggINTg_contrib(double w, const _Event& p) {
      _M_p = &p; _M_weight = w; _M_calc = false; _M_contr = ggINTg;
    }
    //   leading  order contribution qgINTg POSTA
    void qgINTg_contrib(double w, const _Event& p) {
      _M_p = &p; _M_weight = w; _M_calc = false; _M_contr = qgINTg;
    }

//---------------------------------------------------------------------
//---------------------------------------------------------------------
//---------------------------------------------------------------------
//---------------------------------------------------------------------
//---------------------------------------------------------------------
//---------------------------------------------------------------------
//---------------------------------------------------------------------
//---------------------------------------------------------------------
//---------------------------------------------------------------------
   
    //   leading  order contribution señal LO
    void MggHyy_contrib(double w, const _Event& p) {
      _M_p = &p; _M_weight = w; _M_calc = false; _M_contr = LOIntHyy;
    }
    //   leading  order contribution señal LO QG AND QQBAR
    void MqqbarHyy_contrib(double w, const _Event& p) {
      _M_p = &p; _M_weight = w; _M_calc = false; _M_contr = LOIntHyyqqbar;
    }



    //   leading  CTNLO order contribution NLO POSTA
    void CTgg_contrib(double w, const _Event& p,double z1, double z2,int orderct) {
      _M_p = &p; _M_weight = w; _M_calc = false; _M_contr = CTgg;
       _M_CTgg.order = orderct;
      _M_CTgg.z1 = z1; _M_CTgg.z2=z2;
    }
   
//---------------------------------------------------------------------
//---------------------------------------------------------------------
//---------------------------------------------------------------------
//---------------------------------------------------------------------
//---------------------------------------------------------------------
//---------------------------------------------------------------------
//---------------------------------------------------------------------
//---------------------------------------------------------------------
//---------------------------------------------------------------------
   
   
   

    //   leading  order contribution
    void leading_order(double w, const _Event& p) {
      _M_p = &p; _M_weight = w; _M_calc = false; _M_contr = born;
    }

    //    next-to-leading order contributions
    void next_to_leading_order_real(double w, const _Event& p) {
      _M_p = &p; _M_weight = w; _M_calc = false; _M_contr = real;
    }

    void next_to_leading_order_fini(double x1, double xjac1, double x2, double xjac2, double w, const _Event& p)
    {
      _M_p = &p; _M_weight = w; _M_calc = false; _M_contr = fini;

      //   fini specific options
      _M_fini.x1 = x1; _M_fini.xjac1 = xjac1;
      _M_fini.x2 = x2; _M_fini.xjac2 = xjac2;
      _M_fini.mode = 0;
    }

    void next_to_leading_order_finix1() { _M_fini.mode = 0;}
    void next_to_leading_order_finix2() { _M_fini.mode = 1;}
    void next_to_leading_order_fini1()  { _M_fini.mode = 2;}

    void next_to_leading_order_sub(double w, const _Event& p, const _Event& dp, int i, int j, int k)
    {
      _M_p = &p; _M_weight = w; _M_calc = false; _M_contr = sub;

      //   dipole specific options
      _M_dipole.dp = &dp; _M_dipole.i = i;
      _M_dipole.j = j; _M_dipole.k = k;
    }

    //    get the value of the current amplitude
    _Weight operator()( double, double, double=1.0) const;
 //   _Weight operator()(c_mstwpdf& , double, double, double=1.0) const;

     _Weight operator()( double, double, double=1.0,double=1.0,double=1.0) const;
//    _Weight operator()(c_mstwpdf& , double, double, double=1.0,double=1.0,double=1.0) const;

    //    get number of the jets
    unsigned int npar() const { return _M_proc -> npar();}
    unsigned int npow() const { return _M_proc -> npow();}

    //  number of quark flavours
    unsigned int nu() const { return _M_proc -> nu();}
    unsigned int nd() const { return _M_proc -> nd();}
    unsigned int nf() const { return _M_proc -> nf();}

    //    get the contributionb type of the current amplitude
    contrib_type contrib() const { return _M_contr;}
    integral_type integral() const { return _M_itype;}

  private:
    //   pointer to the process
    process_type *_M_proc;

    //   store the current process type and function
    contrib_type _M_contr;
    integral_type _M_itype;

    //   store the event
    const _Event *_M_p;
    double _M_weight;

    //   value of the current amp
    mutable bool _M_calc;
    mutable _Weight _M_amp;

    //   contribution specific variables
    //   dipole contribution
    struct {
      //   dipole event
      const _Event *dp;

      //   dipole indices
      int i, j, k;

      void operator()(process_type *proc, const _Event *p, _Weight& res) const {
	proc -> dipole_term(*p, *dp, i, j, k, res);
      }
    } _M_dipole;

    //   finite contributions
    struct {
      //  finite type
      unsigned int mode;

      //  the x1, x2 integral
      double x1, xjac1, x2, xjac2;

      //   log of the scales
      mutable double lxr, lxf;

      //  scale independent decomposation
      mutable _Weight amp[7];

      void operator()(process_type *proc, const _Event *p) const {
	proc -> fini_term(x1, xjac1, x2, xjac2, *p, amp);
      }

      void operator()(process_type *proc, const _Event *p,double rng1,double rng2) const {
	proc -> fini_term(x1, xjac1, x2, xjac2, *p, amp, rng1, rng2);
      }

      void finix1(_Weight& res) const {
	res = amp[0] + amp[3]*lxf;
      }

      void finix2(_Weight& res) const {
	res = amp[1] + amp[4]*lxf;
      }

      void fini1(_Weight& res) const {
	res = amp[2] + amp[5]*lxf + amp[6]*lxr;
      }
    } _M_fini;





    //   finite contributions
    struct {

      //   log of the scales
      mutable double mur2;
      mutable unsigned int Nf;


      void operator()(process_type *proc, const _Event *p, _Weight& res) const {
	proc -> H2f_term(*p, res, mur2,Nf);
      }


    } _M_H2f;


    //   finite contributions
    struct {

      //   log of the scales
      mutable double muf2;
      mutable double mur2;
      mutable unsigned int Nf;
      int order;
      mutable double z1,z2;

      void operator()(process_type *proc, const _Event *p, _Weight& res) const {
	proc ->  H1_term(*p, res, muf2,mur2,Nf,z1,z2,order);
      }


    } _M_H1;


     //   CTNLO contributions
    struct {

      //   log of the scales
      mutable double muf2;
       mutable double mur2;
      mutable unsigned int Nf;
      int order;
      mutable double z1,z2;

      void operator()(process_type *proc, const _Event *p, _Weight& res) const {
	proc ->  CTNLO_term(*p, res, muf2,mur2,Nf,z1,z2,order);
      }


    } _M_CTNLO;



    //   private member functions
    _Weight update_scales(double, double, double) const;

    //   contributions
    void amp_LO() const;
    void amp_BOX() const;
//-----------------------------------
//-- Interferencia ---- Reals -- Neri    
//-----------------------------------
    void amp_ggffg() const;
    void amp_ggHg() const;
    void amp_ggINTg() const;
    void amp_qgINTg() const;
//-----------------------------------
// Interferencia LO ---- LEAN
//-----------------------------------    
    void amp_MggHyy() const;
    void amp_MqqbarHyy() const;
    

     //   CTgg contributions
    struct {

      //   log of the scales
      mutable double muf2;
       mutable double mur2;
      mutable unsigned int Nf;
      int order;
      mutable double z1,z2;

      void operator()(process_type *proc, const _Event *p, _Weight& res) const {
	proc ->  CTgg_term(*p, res, muf2,mur2,Nf,z1,z2,order);
      }


    } _M_CTgg;
   
//-----------------------------------
//-----------------------------------
    void amp_H1f() const;
    void amp_H1() const;
    void amp_CTNLO() const;
    void amp_CTgg() const;
    void amp_H2f() const;
    void amp_born() const;
    void amp_real() const;
    void amp_dipole() const;
    void amp_fini() const;
    void amp_fini(double,double) const;
  };






 template<typename _Weight, class _Event, class _EvenTraits>
  _Weight amplitude<process<_Weight, _Event, _EvenTraits, 2U, 0U> >::
  update_scales( double mr2, double mf2, double coef) const
  {

	 int iset;
	 iset=parpdfs2gamma_.iset2g;

	 // -- Setting scales for PDF in DYTurbo
     scales::set(sqrt(mr2));
     scales::mcfm();
     pdf::setalphas();



    if(_M_contr == fini) {
      _M_fini.lxr = std::log(mr2);
      _M_fini.lxf = std::log(mf2);
    }

    if(_M_contr == H2f) {
      _M_H2f.mur2 = mr2;
      _M_H2f.Nf = this->nf() ;
          }

        if(_M_contr == H1) {
			 _Weight muleto;
			 for(int imul=0;imul<17;imul++) muleto[imul]=1.0;
      _M_H1.muf2 = mf2;
      _M_H1.mur2 = mr2;
      _M_H1.Nf = this->nf() ;
        //      double xmur = sqrt(mr2);
        // double xalpha_s =AlphaS(xmur)/twopi;
      double xalpha_s = (qcdcouple_.as_)/2./__DEB_PI;

      if(_M_H1.order==2)
	{
	  //  return coef*muleto*__1alfaQED2*xalpha_s*xalpha_s;
	  return coef*muleto*xalpha_s*xalpha_s;
	}


      return coef*muleto*xalpha_s;
          }

    if(_M_contr == CTNLO) {

		_Weight cmuleto;
			 for(int cimul=0;cimul<17;cimul++) cmuleto[cimul]=1.0;


	  _M_CTNLO.muf2 = mf2;
      _M_CTNLO.mur2 = mr2;
      _M_CTNLO.Nf = this->nf() ;


			 if(_M_CTNLO.order==2){
				double cTxalpha_s =  (qcdcouple_.as_)/2./__DEB_PI;
				return coef*cmuleto*cTxalpha_s*cTxalpha_s;

			 }

	  //	 double cxmur = sqrt(mr2);
      //     double cxalpha_s =AlphaS(cxmur)/twopi;
         double cxalpha_s =  (qcdcouple_.as_)/2./__DEB_PI;
         return coef*cmuleto*cxalpha_s;
          }

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
// CTgg
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
    if(_M_contr == CTgg) 
      {
	
	_Weight cmuleto;
	for(int cimul=0;cimul<17;cimul++) cmuleto[cimul]=1.0;
	
	
	_M_CTgg.muf2 = mf2;
	_M_CTgg.mur2 = mr2;
	_M_CTgg.Nf = this->nf() ;
	
	if(_M_CTgg.order==2)
	  {
	    double cTxalpha_s = (qcdcouple_.as_)/2./__DEB_PI;
	    return coef*cmuleto*cTxalpha_s*cTxalpha_s;   
	  }	
	double cxalpha_s =  (qcdcouple_.as_)/__DEB_PI;
	double alphas_porder =  (qcdcouple_.as_)/2./__DEB_PI;
	return cmuleto*coef/( opts.aemmz )/( opts.aemmz)*cxalpha_s*std::pow(alphas_porder *twopi, (int) 2);
      }
//////////////////////////////////////////////////////////////////////////

    //----- calculating the pdf -----
    const _Event& p = *_M_p;
    double s = p[hadron(-1)]*p[hadron(0)];
    double eta1 = p[-1]*p[hadron(0)]/s;
    double eta2 = p[0]*p[hadron(-1)]/s;

    if(_M_contr == fini && _M_fini.mode == 0) eta1 /= _M_fini.x1;
    if(_M_contr == fini && _M_fini.mode == 1) eta2 /= _M_fini.x2;

    // Reals interference ggffg,ggHg,ggINTg
    if(_M_contr == ggINTg ||_M_contr == ggHg || _M_contr ==  ggffg   || _M_contr == qgINTg ){

    eta1 = p[-1]*p[hadron(0)]/s;
    eta2 = p[0]*p[hadron(-1)]/s;

  };


  _Weight retval = GeT_pdf(eta1, eta2, mf2,this->nu(), this->nd());


    //----- calculating alpha_s ------
    unsigned int aspow = _M_proc->npow();
   // double as = pdf->alpha_qcd(this->nf(), mr2);

     if(_M_contr != born) ++aspow;
     if(_M_contr == (LO) ) aspow=0;
     if(_M_contr == (H1f) ) aspow=1;
   //  if(_M_contr == (ggffg) ) aspow=3;
     

    // double mur = sqrt(mr2);
    // CALCULA EL ALPHA_S HAY QUE PONERLA AL CUADRADO SI ES A NNLO
    // double alpha_s =AlphaS(mur)/twopi;

     //Set scales                                                                                                                                                                                    
   
     double alpha_s =  (qcdcouple_.as_)/2./__DEB_PI;

     if(_M_contr != (LOIntHyy) && _M_contr != (LOIntHyyqqbar) && _M_contr != (ggffg) && _M_contr != (ggHg) && _M_contr != (ggINTg) && _M_contr != (qgINTg) )  retval *= coef*std::pow(alpha_s, (int) aspow);
     
     if(_M_contr == (ggffg) || _M_contr == (ggHg) || _M_contr == (ggINTg) || _M_contr == (qgINTg) ){
		  aspow=3.0;
		 // std::cout << "alpha_s   : "  << alpha_s << std::endl;
		  retval *= coef*std::pow(alpha_s, (int) aspow);
	//	 _Weight cmuleto;
	//		 for(int cimul=0;cimul<17;cimul++) cmuleto[cimul]=0.0;
//cmuleto[0]=1.0;
	//	  retval = cmuleto;
		  }
     
     //     cout << " _M_contr : " << _M_contr << endl;

     if(_M_contr == (LOIntHyy) ){ aspow=2; retval *= coef/( opts.aemmz )/( opts.aemmz )*std::pow(alpha_s*twopi, (int) aspow); 
       //retval *= coef/( alpha_qed )/( alpha_qed )*std::pow(0.112, (int) aspow); 
      
       // CANCELO LA INCLUSION DE ALPHA_QED^2 QUE PUSO DESDE USER_2FOT
       //         cout << " alphass nlo-process : " << alpha_s << endl;
       //  cout << " coef nlo-process : " << coef << endl;
       //  cout << " alpha_qed nlo-process : " << alpha_qed << endl;
};
	 if(_M_contr == (LOIntHyyqqbar) ){ aspow=2; retval *= coef/(opts.aemmz )/( opts.aemmz )*std::pow((alpha_s*twopi)*4.0*__DEB_PI, (int) aspow);
           /* cout << " pepa " << coef/(opts.aemmz )/( opts.aemmz )*std::pow((alpha_s*twopi)*4.0*__DEB_PI, (int) aspow) << endl; */
	   /*   cout<< " coef : " << coef << endl; */
           /* cout << "alpha  : " << opts.aemmz << endl; */
	   /* cout << " std::pow((alpha_s*twopi)*4.0*__DEB_PI, (int) aspow) : " << std::pow((alpha_s*twopi)*4.0*__DEB_PI, (int) aspow) << endl; */
	   //cout << " alphass : " << alpha_s << endl;
  };
	// if(_M_contr == (LOIntHyyqqbar) ){ aspow=2; retval *= coef/( alpha_qed )/( alpha_qed );  };
		
		 
		
		
   

    return retval;
  }



  template<typename _Weight, class _Event, class _EvenTraits>
  void amplitude<process<_Weight, _Event, _EvenTraits, 2U, 0U> >::amp_LO() const
  {
    if(!_M_calc) {
      //----- calculate the matrix element squares -----
      _M_proc -> LO_term(*_M_p, _M_amp);

       //----- phase space weight -----
      _M_amp *= _M_weight;

      _M_calc = true;
    }
  }

  template<typename _Weight, class _Event, class _EvenTraits>
  void amplitude<process<_Weight, _Event, _EvenTraits, 2U, 0U> >::amp_H2f() const
  {
    if(!_M_calc) {
      //----- calculate the matrix element squares -----
      _M_H2f(_M_proc,_M_p, _M_amp);

      //----- phase space weight -----
      _M_amp *= _M_weight;

      _M_calc = true;
    }
  }

  template<typename _Weight, class _Event, class _EvenTraits>
  void amplitude<process<_Weight, _Event, _EvenTraits, 2U, 0U> >::amp_H1() const
  {
    if(!_M_calc) {
      //----- calculate the matrix element squares -----
      _M_H1(_M_proc,_M_p, _M_amp);

      //----- phase space weight -----
      _M_amp *= _M_weight;

      _M_calc = true;
    }
  }

  template<typename _Weight, class _Event, class _EvenTraits>
  void amplitude<process<_Weight, _Event, _EvenTraits, 2U, 0U> >::amp_CTNLO() const
  {
    if(!_M_calc) {
      //----- calculate the matrix element squares -----
      _M_CTNLO(_M_proc,_M_p, _M_amp);

      //----- phase space weight -----
      _M_amp *= _M_weight;

      _M_calc = true;
    }
  }

 template<typename _Weight, class _Event, class _EvenTraits>
  void amplitude<process<_Weight, _Event, _EvenTraits, 2U, 0U> >::amp_BOX() const
  {
    if(!_M_calc) {
      //----- calculate the matrix element squares -----
      _M_proc -> BOX_term(*_M_p, _M_amp);

      //----- phase space weight -----
      _M_amp *= _M_weight;

      _M_calc = true;
    }
  }

//-----------------------------------------------------------------------------  
//-----------------------------------------------------------------------------  
//-----------------------------------------------------------------------------  
//------- Intefrencia Reales Neri --- 
//-----------------------------------------------------------------------------  
//-----------------------------------------------------------------------------  
//-----------------------------------------------------------------------------  
  

template<typename _Weight, class _Event, class _EvenTraits>
   void amplitude<process<_Weight, _Event, _EvenTraits, 2U, 0U> >::amp_ggffg() const
   {
    if(!_M_calc) {
      //----- calculate the matrix element squares -----
      _M_proc -> ggffg_term(*_M_p, _M_amp);

   //   std::cout << " from amplitude<process<_Weight, _Event, _EvenTraits, 2U, 0U> >::amp_ggffg() const " << std::endl;

      //----- phase space weight -----
      _M_amp *= _M_weight;

      _M_calc = true;
    }
   }
template<typename _Weight, class _Event, class _EvenTraits>
   void amplitude<process<_Weight, _Event, _EvenTraits, 2U, 0U> >::amp_ggHg() const
   {
    if(!_M_calc) {
      //----- calculate the matrix element squares -----
      _M_proc -> ggHg_term(*_M_p, _M_amp);
      
  //     std::cout << " from amplitude<process<_Weight, _Event, _EvenTraits, 2U, 0U> >::amp_ggHg() const" << std::endl;

      //----- phase space weight -----
      _M_amp *= _M_weight;

      _M_calc = true;
    }
   }
// -- Interference Lean 24/0//2024
template<typename _Weight, class _Event, class _EvenTraits>
   void amplitude<process<_Weight, _Event, _EvenTraits, 2U, 0U> >::amp_qgINTg() const
   {
    if(!_M_calc) {
      //----- calculate the matrix element squares -----
      _M_proc -> qgINTg_term(*_M_p, _M_amp);
      
   //    std::cout << " from amplitude<process<_Weight, _Event, _EvenTraits, 2U, 0U> >::amp_ggINTg() const" << std::endl;

      //----- phase space weight -----
      _M_amp *= _M_weight;

      _M_calc = true;
    }
   }

template<typename _Weight, class _Event, class _EvenTraits>
   void amplitude<process<_Weight, _Event, _EvenTraits, 2U, 0U> >::amp_ggINTg() const
   {
    if(!_M_calc) {
      //----- calculate the matrix element squares -----
      _M_proc -> ggINTg_term(*_M_p, _M_amp);
      
   //    std::cout << " from amplitude<process<_Weight, _Event, _EvenTraits, 2U, 0U> >::amp_ggINTg() const" << std::endl;

      //----- phase space weight -----
      _M_amp *= _M_weight;

      _M_calc = true;
    }
   }

//-----------------------------------------------------------------------------  
//-----------------------------------------------------------------------------  
//-----------------------------------------------------------------------------  
//------- Intefrencia LO --- LEAN 
//-----------------------------------------------------------------------------  
//-----------------------------------------------------------------------------  
//-----------------------------------------------------------------------------   

template<typename _Weight, class _Event, class _EvenTraits>
   void amplitude<process<_Weight, _Event, _EvenTraits, 2U, 0U> >::amp_MggHyy() const
   {
    if(!_M_calc) {
      //----- calculate the matrix element squares -----
      _M_proc -> MggHyy_term(*_M_p, _M_amp);

      //----- phase space weight -----
      _M_amp *= _M_weight;

      _M_calc = true;
    }
   }
   
   template<typename _Weight, class _Event, class _EvenTraits>
   void amplitude<process<_Weight, _Event, _EvenTraits, 2U, 0U> >::amp_MqqbarHyy() const
   {
    if(!_M_calc) {
      //----- calculate the matrix element squares -----
      _M_proc -> MqqbarHyy_term(*_M_p, _M_amp);

      //----- phase space weight -----
      _M_amp *= _M_weight;

      _M_calc = true;
    }
   }
   
     template<typename _Weight, class _Event, class _EvenTraits>
  void amplitude<process<_Weight, _Event, _EvenTraits, 2U, 0U> >::amp_CTgg() const
  {
    if(!_M_calc) {
      //----- calculate the matrix element squares -----
      _M_CTgg(_M_proc,_M_p, _M_amp);

      //----- phase space weight -----
      _M_amp *= _M_weight;

      _M_calc = true;
    }
  }
  
  
   
//-----------------------------------------------------------------------------  
//-----------------------------------------------------------------------------  
//-----------------------------------------------------------------------------  
//-----------------------------------------------------------------------------  
//-----------------------------------------------------------------------------  
//-----------------------------------------------------------------------------     
   

 template<typename _Weight, class _Event, class _EvenTraits>
  void amplitude<process<_Weight, _Event, _EvenTraits, 2U, 0U> >::amp_H1f() const
  {
    if(!_M_calc) {
      //----- calculate the matrix element squares -----

      _M_proc -> H1f_term(*_M_p, _M_amp);

      //----- phase space weight -----
      _M_amp *= _M_weight;

      _M_calc = true;
    }
  }





  template<typename _Weight, class _Event, class _EvenTraits>
  void amplitude<process<_Weight, _Event, _EvenTraits, 2U, 0U> >::amp_born() const
  {
    if(!_M_calc) {
      //----- calculate the matrix element squares -----
      _M_proc -> born_term(*_M_p, _M_amp);


      //----- phase space weight -----
      _M_amp *= _M_weight;

      _M_calc = true;
    }
  }

  template<typename _Weight, class _Event, class _EvenTraits>
  void amplitude<process<_Weight, _Event, _EvenTraits, 2U, 0U> >::amp_real() const
  {
    if(!_M_calc) {
      //----- calculate the matrix element squares -----
      _M_proc -> real_term(*_M_p, _M_amp);


      //----- phase space weight -----
      _M_amp *= _M_weight;


      _M_calc = true;
    }
  }

  template<typename _Weight, class _Event, class _EvenTraits>
  void amplitude<process<_Weight, _Event, _EvenTraits, 2U, 0U> >::amp_dipole() const
  {
    if(!_M_calc) {
      //----- dipole contributions -----
      _M_dipole(_M_proc, _M_p, _M_amp);

      //----- phase space weight -----


      _M_amp *= _M_weight;
      _M_calc = true;
    }
  }

  template<typename _Weight, class _Event, class _EvenTraits>
  void amplitude<process<_Weight, _Event, _EvenTraits, 2U, 0U> >::amp_fini() const
  {
    if(!_M_calc) {
      //----- calculate the matrix element squares -----
      _M_fini(_M_proc, _M_p);

      //----- phase space weight -----
      for(unsigned int i = 0; i < 7; i++)
	_M_fini.amp[i] *= _M_weight;
      _M_calc = true;
    }

    switch(_M_fini.mode) {
    case 0: _M_fini.finix1(_M_amp); break;
    case 1: _M_fini.finix2(_M_amp); break;
    case 2: _M_fini.fini1(_M_amp);  break;
    }
  }

template<typename _Weight, class _Event, class _EvenTraits>
  void amplitude<process<_Weight, _Event, _EvenTraits, 2U, 0U> >::amp_fini(double rng1,double rng2) const
  {
    if(!_M_calc) {
      //----- calculate the matrix element squares -----
      _M_fini(_M_proc, _M_p,rng1,rng2);

      //----- phase space weight -----
      for(unsigned int i = 0; i < 7; i++)
	_M_fini.amp[i] *= _M_weight;
      _M_calc = true;
    }

    switch(_M_fini.mode) {
    case 0: _M_fini.finix1(_M_amp); break;
    case 1: _M_fini.finix2(_M_amp); break;
    case 2: _M_fini.fini1(_M_amp);  break;
    }
  }

  template<typename _Weight, class _Event, class _EvenTraits>
  _Weight amplitude<process<_Weight, _Event, _EvenTraits, 2U, 0U> >::
  operator()(double mr2, double mf2, double coef) const
  {
    _Weight fac = update_scales(mr2, mf2, coef);


    switch(_M_contr) {
    case LO: this -> amp_LO();   break;
    case H1: this -> amp_H1();   break;
    case H2f: this -> amp_H2f();   break;
    case BOX: this -> amp_BOX();   break;
    case ggffg: this -> amp_ggffg();   break;
    case ggHg: this -> amp_ggHg();   break;
    case ggINTg: this -> amp_ggINTg();   break;
    case qgINTg: this -> amp_qgINTg();   break;
    case H1f: this -> amp_H1f();   break;
    case CTNLO: this -> amp_CTNLO();   break;
    case CTgg: this -> amp_CTgg();   break;
    case born: this -> amp_born();   break;
    case real: this -> amp_real();   break;
    case sub:  this -> amp_dipole(); break;
    case fini: this -> amp_fini();   break;
    case LOIntHyy: this -> amp_MggHyy(); break;
    case LOIntHyyqqbar: this -> amp_MqqbarHyy(); break;
    
    default: throw "amplitude<..., 2, 0> : no contribution defined"; break;
    }


    return fac*_M_amp;
  }




template<typename _Weight, class _Event, class _EvenTraits>
  _Weight amplitude<process<_Weight, _Event, _EvenTraits, 2U, 0U> >::
  operator()(double mr2, double mf2, double coef,double rng1,double rng2) const
  {

    _Weight fac = update_scales(mr2, mf2, coef);

   // std::cout << " from operator()(double mr2, double mf2, double coef,double rng1,double rng2) const " << std::endl;
   
  // std::cout << _M_contr << std::endl;

    switch(_M_contr) {
    case LO:     this -> amp_LO();   break;
    case H1:     this -> amp_H1();   break;
    case H2f:    this -> amp_H2f();   break;
    case BOX:    this -> amp_BOX();   break;
    case ggffg:  this -> amp_ggffg();   break;
    case ggHg:   this -> amp_ggHg();   break;
    case ggINTg: this -> amp_ggINTg();   break;
    case qgINTg: this -> amp_qgINTg();   break;
    case H1f:    this -> amp_H1f();   break;
    case CTNLO:  this -> amp_CTNLO();   break;
    case CTgg:  this -> amp_CTgg();   break;
    case born:   this -> amp_born();   break;
    case real:   this -> amp_real();   break;
    case sub:    this -> amp_dipole(); break;
    case fini:   this -> amp_fini(rng1,rng2);   break;
    case LOIntHyy: this -> amp_MggHyy(); break;
    case LOIntHyyqqbar: this -> amp_MqqbarHyy(); break;
      //
    default: throw "amplitude<..., 2, 0> : no contribution defined"; break;
    }

    return fac*_M_amp;
  }



#endif
