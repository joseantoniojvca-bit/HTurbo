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
#ifndef __NLO_PSG_PHASESPACE_N2I2F0_H__
#define __NLO_PSG_PHASESPACE_N2I2F0_H__ 1

#include "nlo-event.h"
#include "nlo-phasespace.h"
#include "psg-phasespace.h"
#include "psg-psgen.h"
//#include "psg-rambo.h"
#include "psg-dipole.h"

template<>
class basic_phasespace<hadronic_event<lorentzvector<double>, hadronic_event_traits<2U,2U,0U> > >
  :  public phasespacel<hadronic_event<lorentzvector<double>, hadronic_event_traits<2U,2U,0U> > >
{
  //  private types
  //   abstract base type
  typedef phasespacel<hadronic_event<lorentzvector<double>, hadronic_event_traits<2U,2U,0U> > > _Base;
  
  //   event_type related type aliases 
  typedef _Base::event_type _Event;
  typedef event_type::iterator _Iterator;
  
  //   phase space generator related aliases
  typedef phasespace_generator<lorentzvector<double>, _Iterator> _PhaseSpaceGen;
  //   typedef rambo<_Iterator> _Rambo;
  
 public:
  typedef _Base::event_type event_type;
  
  //  constructor
 basic_phasespace(const random_generator *rng, double s)
   : _Base(rng), _M_s(s), _M_beta(2.0), _M_eps(1.0), 
    _M_intpsgen(true) {}
  
 basic_phasespace(const random_generator *rng, _PhaseSpaceGen *psg, double s)
   : _Base(rng), _M_s(s), _M_beta(2.0), _M_eps(1.0), _M_q2min(0.0), _M_q2max(s),
    _M_intpsgen(false) {}
  
  ~basic_phasespace() {
    
  }
  
  //   set the important sampling parameters
  void important_sampling(double beta, double eps) {
    _M_beta = beta; _M_eps = eps;
  }
  
  //   set the phase space cuts
  void phasespace_cuts(double, double);
  
  
  
  //   generate the events
  double operator()(_Event&,double&,double&,double&);
  double operator()(_Event&,double,double&,double&,double,double);
  
  double operator()(_Event&,double&,double&);
  double operator()(_Event&);
  double operator()(const _Event&, _Event&);
  
  //-- Phase space -- Interference LO
  double operator()(_Event&,double&,double&,int&);
  double operator()(_Event&,double&,double&,int&,int&,const double x[]);
  double operator()(_Event&,double&,double&,int&,bool,const double x[]);
  double operator()(_Event&,double&,const double x[]);
  double operator()(_Event&,double&,bool,const double x[]); // New qg NLO
  double operator()(_Event& p,bool flagrmed,const double x[]); // Another qg NLO 
  //-- Phase space -- Interference LO - Counter term gg
  double operator()(event_type& p,double &Q2,double &smearing, bool flagrmed,double &z1, double &z2,const double x[]);
  //         double operator()(_Event&,double&,double&,int&,double&,double&,const double x[]);
  // -- Phase space 2->3 all outgoing 
  double operator()(_Event&, bool, double&, const double x[]);
  double operator()(_Event&,bool&); 
 private:
  //    data members
  //   energy of the incomings
  double _M_s, _M_beta, _M_eps;
  double _M_q2min, _M_q2max;
  
  //   base phase space generator (e.g. rambo)
  bool _M_intpsgen;
  // _PhaseSpaceGen *_M_psgen;
  
  //   private static members
    static void _S_safety_cut(const event_type&);
  };



#endif
