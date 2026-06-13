//  Copyright (C) 2002 Zoltan Nagy
//
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
#ifndef __AMPQ2G1P2_H__
#define __AMPQ2G1P2_H__ 1

#include "ampq2g3.h"

  class ampq2g1p2 : private ampq2g3
  {
    //   private types
    typedef std::complex<double> _ComplexD;
    typedef std::pair<double, _ComplexD> _Pair;

  public:
    //   constructor
  ampq2g1p2(const innerprod_type& __x, const random_generator& rng)
    : ampq2g3(__x, rng) {}
    
    
    //   tree level matrix element squared
    double LO_tree(int, int, int, int, int) const;
    
    // H2f

    double H2f(int , int , int ,int ,int ,double ,unsigned int,double) const;
    
    //   tree level matrix element squared H1f
    double LO_H1f(int, int, int, int, int) const;
    
    //   BOX
    double BOXMTR(int, int, int, int, int) const;
    
    //   tree level helicity amplitudes
    
    //   tree level matrix element squared
    double su3_tree(int, int, int, int, int) const;
    
    //  color corraleted amplitudes
    _Pair su3_cc(int, int, int, int, int, int, int) const;
    
    //   amplitudes for the finite part
    void su3_kp(unsigned int, int, int, int, int, int, int, int, su3_kp_i2 *, double=1.0) const;
    
    // ---------------------------------------------------------------------
    // SIGNAL - BACKGROUND INTERFERENCE
    // ---------------------------------------------------------------------
    void ggg_background(int,int,int,int,int, std::vector<std::complex<double> >&) const;
    void ggg_signal(int,int,int,int,int, std::vector<std::complex<double> >&) const;
    std::complex<double> Higgs_prop(double ss) const;
    double wgtR(int, int, int, int, int ) const;
    double wgtHsig(int, int, int, int, int ) const;
    double wgtHint(int, int, int, int, int ) const;
    double AmpIntqgNLO(int, int, int, int, int ) const; // New -- qg channel
    
    std::complex<double> A5yy(std::complex<double>(ampq2g1p2::*)(int,int,int,int,int)const,
                              int,int,int,int,int) const;
    
    /*
      std::complex<double> Mppppp(int , int , int , int , int ) const;
      std::complex<double> Mmmmmm(int , int , int , int , int ) const;
      std::complex<double> Mmpppp(int , int , int , int , int ) const;
      std::complex<double> Mpmmmm(int , int , int , int , int ) const;
      std::complex<double> Mmmppp(int , int , int , int , int ) const;
      std::complex<double> Mppmmm(int , int , int , int , int ) const;
      std::complex<double> Mmpmpp(int , int , int , int , int ) const;
      std::complex<double> Mpmpmm(int , int , int , int , int ) const;
      * */
    // ---------------------------------------------------------------------
    // SIGNAL - BACKGROUND INTERFERENCE LO - LEANDRO
    // ---------------------------------------------------------------------
    std::complex<double> LO_MggH(int, int, int, int, int) const;
    std::complex<double> LO_MHyy(int, int, int, int, int) const;
    double LO_Mggyy(int, int, int, int, int,double &) const;
    double LO_MqqbarH(int, int, int, int, int) const;
    double LO_MqgH(int, int, int, int, int) const;
    
    // ---------------------------------------------------------------------
    // ---------------------------------------------------------------------
    // ---------------------------------------------------------------------

    struct amp_1loop {
      _ComplexD U0, D0, U1, D1;
    };
    
    //   1-loop matrix elements (quark electric charges are included!)
    void matrix_1loop_mpmpp(unsigned int, unsigned int, int, int, int, int, int, amp_1loop *) const;
    void matrix_1loop_mppmp(unsigned int, unsigned int, int, int, int, int, int, amp_1loop *) const;
    void matrix_1loop_mpppm(unsigned int, unsigned int, int, int, int, int, int, amp_1loop *) const;
    void matrix_1loop_mppmm(unsigned int, unsigned int, int, int, int, int, int, amp_1loop *) const;
    void matrix_1loop_mpmpm(unsigned int, unsigned int, int, int, int, int, int, amp_1loop *) const;
    void matrix_1loop_mpmmp(unsigned int, unsigned int, int, int, int, int, int, amp_1loop *) const;

    //   1-loop level matrix element squared (quark electric charges are included!)
    void su3_1loop(unsigned int, unsigned int, int, int, int, int, int, double *) const;
    void su3_1loop_mch(unsigned int, unsigned int, int, int, int, int, int, double *) const;
    void su3_1loop_mch(unsigned int, unsigned int, int, int, int, int, int, double *,double,double) const;

  private:
    //   private members
    _ComplexD amphtree(int, int, int, int, int) const;
    static void matrix_1loop(unsigned int, unsigned int, const _AmpPrim *, const _AmpPrim *, amp_1loop *);
  };


// Funciones para H -> gam gam

std::complex<double> ff( double );
std::complex<double> F1( double );
std::complex<double> F1on2( double );
std::complex<double> CorcheteMartin(double );
std::complex<double> SumaMartin(double );  

double LLi(const double &x);

// Li is the dilogarithm.  It assumes 0 < x < 1.

double LLi2(double r, double theta);

// Li2 is the Real part of the dilog function of
// the complex variable Z = r * exp(I*theta), 
// with r>0 and -PI < theta < PI.  
// Li2 = - Real( int(0..Z) log(1.0-z)/z ).
//
// (If theta=0 and r < 1, it reduces to Li(r).)

double LLi2(double x);

// This is the same as Li2(r, theta) with (r=x,theta=0) for x>0 an
// (r=-x,theta=PI) for x<0.  (Note that this only gives the Real part.)

double LLi3(double x);

double LLi4(double x);

complex<double> LL0(double r1, double r2);

complex<double> LL1(double r1, double r2);

complex<double> LL2(double r1, double r2);

complex<double> LLs1(double r1, double r2, double r3, double r4);

//complex<double> Mpppp,Mmppp,Mmmpp,Mmpmp,Mmppm;

     complex<double> Mppppp(int i1, int i2, int i3, int i4, int i5);
     complex<double> Mmpppp(int i1, int i2, int i3, int i4, int i5);
     complex<double> Mmmppp(int i1, int i2, int i3, int i4, int i5);
     complex<double> Mmpmpp(int i1, int i2, int i3, int i4, int i5);
     complex<double> Mmmmmm(int i1, int i2, int i3, int i4, int i5);
     complex<double> Mpmmmm(int i1, int i2, int i3, int i4, int i5);
     complex<double> Mppmmm(int i1, int i2, int i3, int i4, int i5);
     complex<double> Mpmpmm(int i1, int i2, int i3, int i4, int i5);

 void Lggg_background(std::vector<std::complex<double> >& helamps) ;

#endif
