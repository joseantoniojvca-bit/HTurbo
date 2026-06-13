// wadp.h 
#if !defined( _WADP_H_ )
#define _WADP_H_



class wadp
{
       private:
           
            double RESULT,ERR,FWW;
      
                    
           
       public:
           
           wadp(double res=0.0, double error=0.0, double fww=0.0):
             RESULT(res), ERR(error),FWW(fww){

                                    }

          
           
           void asignawadp(double a, double b,double c);
          
           double WRE() const {return RESULT;}
           double WER() const {return ERR;}
           double WFW() const {return FWW;}
           
           
           wadp operator+( wadp  v);
           
           wadp operator-( wadp  v);
      
           wadp operator*( wadp  v);
           

 
           
        
        



};

#endif // _WADP_H_
             
                 
