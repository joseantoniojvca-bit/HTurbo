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
#ifndef __H1_NLO_H__
#define __H1_NLO_H__ 1



#include "qt_splitt.h"



#include "nlo-process.h"

#include "mstwpdf.h"

#include "nlo-color.h"
#include "hhc2ph-event.h"

#include "hhc2ph-process.h"

double  ampBORN( event_hhc2ph );
  
  double  H1( event_hhc2ph );
  
  weight_hhc2ph  procH1(c_mstwpdf *, event_hhc2ph ,double ,double , double , double ,int ,int ,int ,double );

//CTNLO gg Leandro/////////////////////////////////////////////////////
 
  std::complex<double> LO_MggHB( event_hhc2ph );
 
  std::complex<double> LO_MHyyB( event_hhc2ph );
  
 double LO_MggyyB( event_hhc2ph , double &);
  
 double _amp_MggHyyB( event_hhc2ph );

// Funciones para H -> gam gam


std::complex<double> CorcheteMartin(double );
std::complex<double> SumaMartin(double );  
 

#endif





  
  
