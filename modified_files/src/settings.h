#ifndef settings_h
#define settings_h

#include <vector>
using namespace std;

#include <algorithm>
#include <map>
#include <string>
#include <fstream>

inline bool file_exists (string name)
{
  ifstream f(name.c_str());
  return f.good();
}

class InputParser {
    public:
        // constructor
        InputParser( string _filename = "", string _charset="#=[ ]", string _white=" \t");
        ~InputParser();
        // getters
        double GetNumber(string name);
        string GetString(string name);
        bool GetBool(string name);
        void GetVectorDouble(string name, vector<double> &vec);
        void parse_file(const string fname);
    private :
        // functions
        void trim(string & str);
        void has_key(const string key);
        // data members
        map<string,string> data; ///< Data table.
        string filename; ///< Input file.
        // parser characters
        char Ccommnt; 
        char Cassign;
        char CopenAr;
        char CclosAr;
        char CdeliAr;
        string Swhite;
};

// Simple program quiting exception
struct QuitProgram : public std::runtime_error { QuitProgram(string msg) : std::runtime_error (msg) {}; };

// Forward declaration
namespace cxxopts{
    class Options;
}
namespace po=cxxopts; // inspired by po = boost::program_options

class settings
{
public:
  settings() {};
  void parse_options(int argc, char * argv[]);
  void readfromfile(const string fname);
  void check_consistency();
  void parse_binning(string name, vector<double> &vec, po::Options &args);

  // private:
  void dumpAll();
  void dumpI(string var, int    val );
  void dumpD(string var, double val );
  void dumpS(string var, string val );
  void dumpB(string var, bool   val );

  // string helpers
  void ToLower(string &val){std::transform(val.begin(), val.end(), val.begin(), ::tolower);}
  void ToUpper(string &val){std::transform(val.begin(), val.end(), val.begin(), ::toupper);}
  vector<string> Tokenize(string val, char Delim=',');
  // is number: http://stackoverflow.com/a/16575025
  bool IsNumber(const string &s);

  //process settings
  double sroot;
  int    ih1;
  int    ih2;
  int    nproc;

  //nuclear beams
  bool nuclearpdf;
  double    Z1;
  double    Z2;
  double    A1;
  double    A2;
  
  //resummation or fixed order switch
  bool fixedorder;
  bool primed;
  
  //order
  int    order;       //Main order
  int    order_sudak; //Order of the Sudakov
  int    order_hcoef; //Order of the H coefficients
  int    order_evol;  //Order of the PDF evolution
  int    order_expc;  //Order of the C exponentiation
  int    order_ct;    //Order of the counter term
  int    order_vjet;  //Order of the V+jet
  //  int    order_alphas; //Order of the alphas running

  int    order_ct_temp = 0;    //Internal variable used to split CT into as1+as2 and as3 pieces
  
  // QED and Mixed QED+QCD corrections
  bool qed;             // ON/OFF
  int qedorder;         //Main order of the QED corrections 0,1,2
  int qedorder_sudak; //QED order of the Sudakov		 
  int qedorder_hcoef;	//QED order of the H coefficients	 
  int qedorder_evol;	//QED order of the PDF evolution	 
  int qedorder_expc;	//QED order of the C exponentiation	 
  int qedorder_ct;	//QED order of the counter term	 
  int qedorder_vjet;	//QED order of the V+jet
  bool alpha0;          //use alphaQED(0) instead of alphaQED(mur)
  
  bool qbox;
  bool qtriang;

  //N4LL uncertainties
  double A5;
  bool a5lc;
  double H4;
  double gam4NSp,gam4NSm,gam4NSv;
  double gam4SGqq,gam4SGqg,gam4SGgq,gam4SGgg;
  double C4qq,C4qg,C4qqb,C4qqp,C4qqbp;
  int sgapprox;
  bool nsapprox;
  
  //Non-perturbative form factor
  int npff;
  double g1,g2,g3;     //Gaussian
  double g1x,g1a,g1b;  //x-dependent g1
  double g1c,g1m,g1f;  //y-dependent g1
  bool gjxnspace;      //x-dependent g1 in Mellin space
  double e;            //Exponential
  double q;            //Quartic
  double sx;           //Sextic
  double g0;           //Collins-Rogers
  double bref;         //Collins-Rogers reference blim
  double qCS;          //CS kernel quartic term
  double Q0;           //reference mass
  double a2,a2p;       //Dokshitzer, Marchesini, Webber
  double lambda,g2A,g2B,NA,sigmaA,alphaA,NB,sigmaB,alphaB; //Pavia19
  double Linf,L2,L4,linf,l2,l4; //SCETlib
  
  //Flavour dependent g1
  bool flavour_kt;
  double g1_uv = 0.5;
  double g1_us = 0.5;
  double g1_dv = 0.5;
  double g1_ds = 0.5;
  double g1_ss = 0.5;
  double g1_ch = 0.5;
  double g1_bo = 0.5;
  double g1_gl = 0.5;

  //PDF settings
  string LHAPDFset    ;
  int    LHAPDFmember ;
  string LHAPDFset_beam2;
  int    LHAPDFmember_beam2;
  bool   externalpdf;

  //alphas running
  bool alphaslha;

  //set the order of the evolution not smaller than the LHAPDF order
  bool orderlha;

  //include mur scale variations in the PDF evolution
  bool murpdf;
  
  //set blim_pdf not larger than b0/Qmin
  bool blimlha;
  
  //functional forms of the QCD scales
  int fmures;
  int fmuren;
  int fmufac;
  
  //scale factors for the QCD scales
  double kmures;
  double kmuren;
  double kmufac;

  //scale factors for the matching scales
  double kmuc;
  double kmub;
  double kmut;
  
  //  double a_param;
  //series acceleration
  int accel;
  
  //EW parameters
  int ewscheme;
  double Gf, hmass, zmass, wmass;
  double xw, aemmz;
  double hwidth, zwidth, wwidth;

  //Running width propagator
  bool runningwidth;

  //Running width to fixed width translation
  bool conv2fixw;
  
  //CKM matrix
  double Vud, Vus, Vub;
  double Vcd, Vcs, Vcb;

  //Z/gamma* coupling
  double Zuu, Zdd, Zss, Zcc, Zbb;

  //resonance mass and width (used for breit wigner unweighting)
  double rmass, rwidth;

  // photon,Z switches
  bool useGamma, useZ;

  //integration boundaries
  //double ylow;
  //double yhigh;
  //double mlow;
  //double mhigh;

  //Resummation damping
  bool damp;
  double dampk, dampdelta;
  int dampmode;

  //Kinematic limit in qt
  bool qtlim;
  
  //Resummation scheme
  int resscheme;
  
  //Resummation cutoff
  double qtcutoff;
  double xqtcutoff;

  //cutoff on invariant mass between emitted and radiator fot V+J
  double mcutoff;
  double mcutoffas3;

  //Nagy parameters
  double aii,aif,afi,aff;
  
  //Modified logarithms
  bool modlog;
  int p;
  
  //qtcut
  double xqtcut, qtcut, qtfpc;
  double xqtcutas3, qtcutas3;

  //taucut
  double taucut;
  int ftau;
  double ptjetmin;
  bool tauboost;
  
  //integration settings
  int rseed;

  //dimension of integration for the resummed part
  int intDimRes;
  bool resint0d, resint1d, resint2d, resint3d, resintvegas;

  //dimension of integration for the born configuration
  int intDimBorn;
  bool bornint1d, bornint2d, bornintvegas4d, bornintvegas6d;
  
  //type of integration for the counterterm
  int intDimCT;
  bool ctint0d, ctint1d, ctint2d, ctint3d, ctintvegas6d, ctintvegas8d;

  //type of integration for the finite order
  int intDimFO;
  bool foint0d, foint3d;
  
  //type of integration for the V+j at LO
  int intDimVJ;
  bool vjint0d,vjint3d, vjint5d, vjintvegas7d;

  //type of integration for the FPC term
  int intDimFPC;
  bool fpcint2d, fpcintvegas8d;
  
  //type of integration, automatic selector
  bool BORNquad, CTquad, VJquad, FPCquad;
  
  //term switches
  bool doBORN;
  bool doCT;
  bool doVJ, doVJREAL, doVJVIRT;
  bool doCTAS3;
  bool doVJAS3, doVJVV, doVJVR, doVJRR;
  bool doFPC;

  //compute CT+VJ
  bool finord;
  
  //multithreading model
  int threading;
  
  //Cuba settings
  int cubaverbosity;
  int cubacores;
  int cubanbatch;
  int niterBORN;
  int niterCT;
  int niterVJ;
  int vegasncallsBORN  ;
  int vegasncallsCT   ;
  int vegasncallsVJLO   ;
  long long int vegasncallsVJREAL ;
  int vegasncallsVJVIRT ;
  int vegasncallsVJVV ;
  int vegasncallsVJVR ;
  long long int vegasncallsVJRR ;
  bool vegascollect;
  bool vegascorr;
  
  //cubature settings
  bool pcubature;
  double relaccuracy;
  double absaccuracy;

  //smolyak
  int level;
  
  //costh boundaries
  double costhmin, costhmax;

  //xF cuts
  double xfmin, xfmax;
  
  //lepton fiducial cuts
  bool makecuts;
  double lptcut, lycut; //charged leptons
  double mtcut, etmisscut;
  double lepptcut, lepycut, alpptcut, alpycut; //lepton and antilepton
  double lcptcut, lcymin, lcymax, lfptcut, lfymin, lfymax; //lc and lf are absolute-rapidity-ordered leptons
  double cthCSmin, cthCSmax;

  //integration types and settings for costh phi_lep phase space
  bool cubaint;
  int suavepoints;

  bool trapezint;
  int nphitrape;
  int ncstart;

  //quadrature rule in phi_lep
  bool quadint;
  int phiintervals;
  int phirule;

  //settings for Bessel integration
  double bintaccuracy;

  //numerical integration of the Sudakov
  bool numsud;
  //use a numerical solution for the running of alphas
  bool asrgkt;
  //running coupling in the sudakov
  int asrun;
  //use the CMW scheme
  bool cmw, cmwfo;

  //HF 
  bool deltaismb, deltaismc, deltafsmb, deltafsmc;

  //ZMVFNs
  bool vfncqg;
  
  //numerical integration of the C exponentiation
  bool numexpc;
  
  //b-space prescription
  int bprescription;
  int fbstar;
  
  //blim parameter of the bstar prescription (acts as an IR cut-off)
  double blim;
  double blim_pdf, blim_sudakov, blim_expc;

  //force bstar prescription
  bool bstar_pdf, bstar_sudakov, bstar_expc;
  
  //arg(z) in the complex plane for the minimal prescription
  double phibr;
  //select the point bc, where the integration contour is bended in the complex plane, as a fraction of b_L = ... (Landau singularity)
  double bcf;
  
  //settings for the Mellin transform
  int mellintr;

  //settings for the Mellin inversion
  int mellininv;
  int mellinintervals;
  int mellinrule;
  double zmax;
  int ncycle;
  double cpoint;
  double cshift;
  double phi;
  int mellincores;
  bool mellin1d;
  int melup;
  bool xspace;
  bool ctmellin;
  
  //settings for x-to-N Mellin transform
  int pdfrule;

  //Invariant mass rule in 0d integration
  int mrule;
  
  //settings for rapidity integration in 2D resummed piece
  int yintervals;
  int yrule;

  //settings for qt integration in 2D counter term
  int qtintervals;
  int qtrule;

  //settings for alfa beta scaled-PDF integration in counter term and born fixed order
  int abintervals;
  int abrule;

  //settings for the phi integration in the V+J 5d LO term when makecuts is false
  int vjphirule;

  //settings for the z1, z2 integration in the V+J 3d NLO singular term
  int zrule;

  //settings for the x integration in the V+J 3d delta term
  int xrule;
  
  //qt-recoil prescriptions
  bool qtrec_naive, qtrec_cs, qtrec_kt0;

  // PDF errors
  bool PDFerrors;
  int totpdf;
  
  //debug settings
  bool timeprofile;
  bool verbose;
  bool gridverbose;

  //output settings
  bool texttable;
  bool redirect;
  bool unicode;
  bool silent;
  bool makehistos;
  string output_filename; // Output Filenames
  
  //resummed code in C++
  bool resumcpp = true; //use C++ code for resummation     

  //counter term code in C++
  bool ctcpp = true;   //use C++ code for the counter term
  
  //dyres, pegasus,lhapdf PDF evolution
  int evolmode;

  //Do not factorize muf variations with expc
  bool mufvar;

  //iterative or truncated solution of the PDF evolution
  bool iterative;

  //in backward evolution mode switch off heavy flavour PDFs below threshold
  int hfmode;
  
  //smallx resummation in the PDF evolution
  bool smallx;
  //number of steps in the path integral
  int nsub;
  
  //Evolve PDFs from mufac instead of mures
  //bool mufevol;

  //Comute PDF evolution U(mub,Q) = U(muF,Q)*U(mub,muF)
  /*
  pdf(b) = U(b,F)*pdf(F)
    = U(b,Q)*U(Q,F)*pdf(F)
    = (U(Q,F)*U(b,Q)-[U(b,Q),U(Q,F)])*pdf(F)
    = (U(Q,F)*U(b,F)*U(F,Q)-[U(b,Q),U(Q,F)])*pdf(F)
    = (U(Q,F)*(U(F,Q)*U(b,F) - [U(b,F),U(F,Q)])-[U(b,Q),U(Q,F)])*pdf(F)
    = (U(Q,F)*U(F,Q)*U(b,F) - U(Q,F)*[U(b,F),U(F,Q)] - [U(b,Q),U(Q,F)])*pdf(F)
    = (U(Q,F)*U(F,Q)*U(b,F) - U(Q,F)*U(b,F)*U(F,Q) + U(Q,F)*U(F,Q)*U(b,F) - U(b,Q)*U(Q,F) + U(Q,F)*U(b,Q))*pdf(F)
    = (U(Q,F)*U(F,Q)*U(b,F) - U(Q,F)*U(b,Q) + U(b,F) - U(b,F) + U(Q,F)*U(b,Q))*pdf(F)
    = U(Q,F)*U(F,Q)*U(b,F)*pdf(F)
  */
  bool commute;

  //switches for C exponentiation
  int expc, ntaylor;
  double expcreg;

  //sum all logs before exponentiation
  bool sumlogs;
  
  //bin width normalisation
  bool ptbinwidth, ybinwidth, mbinwidth;

  //ds/dqt^2
  bool dsdqt2;
  bool dsdxf;
  bool edsdp3;
  bool m3dsdm;

  bool exactkin;
  
  // Force to loop over all bins even you have all Vegas integrands
  bool force_binsampling = false;

  bool lepbounds;
  
  // Calculate helicity cross sections
  int helicity;

  // New implementation LEAN
  // -- Photons
  double mtop, mbottom, mcharm, mtau;
  bool removebr = false;
  string calc_type = "sig";
  string decay_type = "full";
  bool apply_crack = false;
  double phpt_cut_h,phpt_cut_s,phycut,crack_eta_max,crack_eta_min;
  double ph_iso_R0, ph_iso_E0 , n_isol;
  double SigmaSmearing;
  int vegasncallsqqbarLO;
  int vegasncallsggNLO;
  bool LO_qg_ON, LO_qqbar_ON,gg_only,qg_only;
  bool withsmearing,prodmtinf;
//  double Rgamgam;
  double eta_jet;
  bool Dixon_eta_jet_cut;

  //optional default histograms
  bool h2d, h3d;
  bool hetal, hptl, hptnu, hmt;
  bool haiqt, haiy, haim;
};

class binning
{
 public:
  binning() {};
  void readfromfile(const string fname);
  void GetBins(string name,vector<double> &bins);
  // private:
  vector <double> qtbins;
  vector <double> ybins;
  vector <double> mbins;
  vector <double> ptlbins;
  vector <double> etalbins;
  vector <double> hist_qt_bins;
  vector <double> hist_y_bins;
  vector <double> hist_m_bins;
 private:
  InputParser in;
};


extern settings opts;
extern binning bins;

#endif
