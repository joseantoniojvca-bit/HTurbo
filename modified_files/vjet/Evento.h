#ifndef EVENTOS_H
#define EVENTOS_H 1

double U_radius(const lorentzvector<double>& , const lorentzvector<double>& );
bool U_photon_izolation(const event_hhc2ph& , double , double , double );

class _Evento { 
 
public:
    double _M_s;
    double q2max;
    double q2min;
    void Set_ps(_Evento & );
    void Set_VEGAS(_Evento & );
// mutable event_hhc2ph p1;
// mutable event_hhc2ph p2;
    
    int i0,inn,ifnn;
    int n0,nn0,fnn0;
    int ncl0,nncl,fnncl;
    int iifborn,iifnnlo,iinnlo,iihyylo,iirhyy;
    char *cppname[25];
  
} ;


#endif


