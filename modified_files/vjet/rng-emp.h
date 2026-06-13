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
#ifndef __RNG_EMP_H__
#define __RNG_EMP_H__ 1

#include "hep-rng.h"
//#include "Integral_class.h"






extern "C" 
{
extern struct {
// La dimension de estas variables tienen que ser las definidas en
// el archivo fortran
          
         double xxextern[30];
   
} xxrands_;
}



  class rng_emp :  public random_generator
  {
  public:
    //   store the state of the generator
    struct state_type {
           int mti;
     unsigned long flg;
      unsigned long flg2;
    };
    
    
   void integrandob( double &,double &,double &,
  double &,double &,double &,double &,double &,double &);
  
  //  int const Intdim;
   
    double xx[50]; // Tiene la maxima dimension que puede tener la Integral - Si se reajusta desde el mail cambiar
    

    //   constructors
    explicit rng_emp (int s = -238, int Intdim=0) 
      : random_generator(0.0, 1.0) { this -> set(s);
      
     
           
      for(int i=0;i<Intdim;i++) xx[i]=1.;//xrands_.xxextern[i];              
                    
      
      } // - Esto en realidad no setea lo limites es empty
    
    //   seeding the generator
    void set(int); // Este deberia ser virtual - Si lo es
    
    //   generate unsigned long and double random number
   // unsigned long int get() const;
    double get_double() const;
    
   
    
  //  void setdim(int numdim){ Intdim = numdim;};
    
  private:
    mutable state_type _M_state;
  };
  




#endif
