#ifndef __SET_ROOT_WORKSPACE__
#define __SET_ROOT_WORKSPACE__ 1

#include <iostream>
#include "hep-rng.h"
#include "rng-emp.h"
#include "Integral_class.h"
#include "Evento.h"
#include "hhc2ph-process.h"
#include "hhc2ph1jet.h"
#include "mstwpdf.h"
#include "SetVEGAS.h"


#include <Riostream.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1D.h>
#include <TNtuple.h>

#include <string>
using std::string;

double contrib;
unsigned int contpt;

bool fillntuple;
int nhists;

// declare here The Root Histograms, NTuples and TTres
TFile *F;
TNtuple *T;
TH1D* hists[50];
string histvar[50];

#endif
