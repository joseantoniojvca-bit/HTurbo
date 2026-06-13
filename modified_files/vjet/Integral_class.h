#ifndef INTEGRAL_CLASS
#define INTEGRAL_CLASS


#include <iostream>
#include "hhc2ph-event.h"
#include "hhc2ph-phasespace.h"
#include "hhc2ph-process.h"
#include "mstwpdf.h"


class Integral_class  {

  public:
    double vweght;
    //  phasespace_hhc2ph *_M_ps;
    process_hhc2ph *_M_proc;
    //     c_mstwpdf *_M_pdf;
    //     ~Integral_class() {
    //     delete _M_pdf;
    //  }

    double calcula_nlo(); // Calcula Real NLO
    double calcula_LO(); // Calcula born
    double calcula_nnlo(); // Calcula nnlo
    void Cuts(event_hhc2ph , bool& ,double,double);
    double Integrando_Fborn(double );
    double Integrando_boxh(double );
    double Integrando_h2(double );
    double Integrando_nnlo(double );
    double Integrando_fini(double );
    double Integrando_Full_nlo(double );
    double Integrando_Full_nnlo(double );
    double Integrando_Full_nlo_CT(double );
    double Integrando_Full_nnlo_CT(double );
// Interferencia NLO -- LEAN    
    double Integrando_ggffg(double);
    double Integrando_gghg(double);
    double Integrando_ggintg(const double x[],event_hhc2ph & p);
    double Integrando_qgintg(const double x[],event_hhc2ph & p);
// Interferencia LO  -- LEAN
    double Integrando_gghyy(double);
    double Integrando_qqbarhyy(const double x[],event_hhc2ph & p);
    double Integrando_CTgg(const double x[],event_hhc2ph & p);
//  weight_hhc2ph  Integral_class::GeT_pdf(c_mstwpdf *,double , double, double , unsigned int , unsigned int )
    void userfunc(event_hhc2ph,double & ,amplitude_hhc2ph& ,double,double);

};







#endif


