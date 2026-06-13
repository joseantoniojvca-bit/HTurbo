// Solo para llamar a mis Randoms en la Forma de Nagy

#include "rng-emp.h"


  
#define N 40 // Numero máximo de dimensiones

  
  
  
  
  
  double rng_emp::get_double() const{
    
    if(_M_state.mti >= N) {
     
    
      std::cout << " Estás requiriendo más randoms que la dimensión de tu integral - Ajustar dimensión " << endl;
      
     _M_state.mti = 0;
    }
     
 //     cout <<  "ANtes  " << _M_state.mti << endl;
     
   
   
 //   if( _M_state.mti > 16 && (_M_state.flg != 31415985742320053281282046584520) ) {
   
//   _M_state.mti=17; _M_state.flg= 31415985742320053281282046584520;
   
//    }
    
    
    int rng_iden = _M_state.mti;
    
    
    //cout << "Despues   :"<< _M_state.mti << endl;
    
     _M_state.mti++;
    
   //  _M_state_ext = _M_state.mti;
        
    return xx[rng_iden] ;
    
    
  }
  

  
  void rng_emp::set( int s)
  {
    if(s == 0) s = -238;
    
   
    
    _M_state.mti = 0; // Setea el uso del primer random
  }

