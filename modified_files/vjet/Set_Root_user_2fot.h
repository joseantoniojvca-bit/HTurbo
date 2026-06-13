#ifndef __SET_ROOT_USER_2FOT__
#define __SET_ROOT_USER_2FOT__ 1


#include "Integral_class.h"
#include "hhc2ph-event.h"
#include "hhc2ph-phasespace.h"
#include "hhc2ph-process.h"
#include "hep-rng.h"
#include "rng-emp.h"
#include "Evento.h"

#include <algorithm>
// Root
/* #include <Riostream.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1D.h>
#include <TNtuple.h>
*/
extern double contrib;
extern unsigned int contpt;

/*
extern "C" {
    void topfill_(int* hident,double* var,double* weight);
}

inline void fillhist(int hident, double var, double weight) {
    topfill_(&hident, &var, &weight);
}

extern bool fillntuple;
extern int nhists;

// declare here The Root Histograms, NTuples and TTres
extern TFile *F;
extern TNtuple *T;
extern TH1D* hists[50];
extern string histvar[50];
*/

#endif
