//      weight.cpp
//      
//      Copyright 2010 Leandro <leandro@ubuntu>
//      
//      This program is free software; you can redistribute it and/or modify
//      it under the terms of the GNU General Public License as published by
//      the Free Software Foundation; either version 2 of the License, or
//      (at your option) any later version.
//      
//      This program is distributed in the hope that it will be useful,
//      but WITHOUT ANY WARRANTY; without even the implied warranty of
//      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//      GNU General Public License for more details.
//      
//      You should have received a copy of the GNU General Public License
//      along with this program; if not, write to the Free Software
//      Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
//      MA 02110-1301, USA.


#include <iostream>
#include <iomanip>
#include <algorithm>

#include <fstream>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>



#include "wadp.h"









void wadp::asignawadp(double v1,double v2,double v3)
{
	
	
	RESULT=v1;
	ERR=v2;
	FWW=v3;
		
	
}

wadp wadp::operator+(wadp wt)
{ 

	return wadp(RESULT + wt.RESULT, ERR + wt.ERR , FWW + wt.FWW);
	
	
 }

wadp wadp::operator-(wadp wt)
{ 

	return wadp(RESULT - wt.RESULT,ERR-wt.ERR,FWW - wt.FWW);
	
	
 }

wadp wadp::operator*(wadp wt)
{ 

	return wadp(RESULT * wt.RESULT,ERR*wt.ERR,FWW * wt.FWW);
	
 }
 
 wadp realwadpprod(double b, wadp v)
{return 	wadp(b*(v.WRE()),b*(v.WER()),b*(v.WFW()));
}

double sumawadp(wadp wtt)
{
	
	return wtt.WRE() + wtt.WER() + wtt.WFW();
	
}

void coutwadp(const wadp& wt)	
{

   
   std::cout <<"[ " ;
   std::cout << std::setw(7) <<  wt.WRE() << " , " ;
   std::cout << std::setw(7) <<  wt.WER() << " , " ;
   std::cout << std::setw(7) <<  wt.WFW() << "  ";
   std::cout << " ]";
   std::cout << "\n\n";
	
	
   
	
	
} 



