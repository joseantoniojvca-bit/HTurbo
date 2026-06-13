#ifndef _params_
#define _params_
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>
//
// cm energy
//
extern double ecm;

//
// acceptance cuts
//
// photon kinematic cuts
extern double ph1_ptmin;
extern double ph2_ptmin;
extern double ph_eta_max;
extern double crack_eta_min;
extern double crack_eta_max;
extern double ph_deltaR_min;
// photon isolation cuts
extern double ph_iso_R0;
extern double ph_iso_E0;
// diphoton invariant mass cuts
extern double mgg_min;
extern double mgg_max;
extern double muR;
extern double muF;
extern std::string Scale;


//
// output options
//
extern bool SaveTopdHistos;
extern bool SaveRootHistos;
extern bool SaveRootNtuple;

//
// 2gNNLO options
//
// contributions to switch on/off (1/0)
extern int _iifborn; // BORN
extern int _iifnnlo; // NLO 
extern int _iinnlo; // NNLO
//extern int _iilean; // NNLO
extern int _iihyylo; // Interference LO
extern int _iirhyy; // Interference Reals NLO
// run option
//     1: RESTART
//     2: CONTINUE
//     3: RESTART WITH PREVIOUS GRIDS
extern int _i0  ; // BORN
extern int _ifnn; // NLO 
extern int _inn ; // NNLO
// Number of iterations 
// With the negative number, 
// you are forcing to the integrator to reach the desired number of Iteration
// In the case of Root histograms use a positive number of Iterations ONLY!!!!!!!!!!!!!!!
// Because in this case, the code do not accumulate weights in Root.
extern int _n0  ; // BORN
extern int _fnn0; // NLO
extern int _nn0 ; // NNLO
// Number of calls per iteration
extern int _ncl0 ;  // BORN
extern int _fnncl; // NLO
extern int _nncl ;  // NNLO



void ReadParams(const char* inputfile);
void DumpParams();
#endif
