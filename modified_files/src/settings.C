#include <cmath>
#include <cstring>
#include <iostream>
#include <stdexcept>

#include "settings.h"
#include "dyres_interface.h"
#include "interface.h"
#include "phasespace.h"
#include "gaussrules.h"
#include "pdf.h"
#include "resconst.h"
#include "HistoHandler.h"

// CXX option parser: https://raw.githubusercontent.com/jarro2783/cxxopts/master/src/cxxopts.hpp
#include "cxxopts.hpp"

#include "handy_typdefs.h"


settings opts;
binning bins;

void settings::parse_options(int argc, char* argv[]){
    // Declare the supported options.
    po::Options args(argv[0], " [config.in] \n\n"
            " Fast Drell-Yan Monte Carlo and quadrature integrator. \n\n"
            " NOTE: Command line options are overiding the default and config file settings."
            );
    // Hidden arguments
    args.add_options("Hidden")
        ("conf_file"     , "Name of output file. ", po::value<string>()->default_value("")  )
    ;
    // Program options
    args.add_options("")
        ("h,help"            , "Print this help and die."       )
        ("v,verbose"         , "Be verbose"                     )
      //        ("q,small-stat"      , "Set quick run with small stat." )
        ("p,proc"            , "Set process [z0/wp/wm]"                              , po::value<string>() )
        ("c,collider"        , "Set beam conditions [tev2/lhc5/lhc7/lhc8/lhc13/lhc14]"    , po::value<string>() )
        ("o,order"           , "Set order [0:LO, 1:NLL+NLO, 2:NNLL+NNLO, 3:N3LL+N3LO]"            , po::value<int>() )
        ("f,fixedorder"      , "Set fixed order only"              )
        ("e,resummation"     , "Set resummation"                )
        ("P,primed"          , "primed orders NLL', NNLL', N3LL', N4LL'")
        ("U,unprimed"        , "unprimed orders NLL, NNLL, N3LL, N4LL")
        ("t,term"            , "Set term [BORN,CT,VJ,VJVIRT,VJREAL,CTAS3,VJAS3,VJVV,VJVR,VJRR,FPC,ALL]"   , po::value<string>() )
        ("r,seed"            , "Set random seed [integer]"                           , po::value<int>()    )
        ("s,pdfset"          , "Set PDF set [LHAPDF name]"                           , po::value<string>() )
        ("m,pdfvar"          , "Set PDF member [integer/all]"                        , po::value<string>() )
        ("u,evmode"          , "Set PDF evolution mode [integer]"                    , po::value<int>() )
        ("g,gpar"            , "Set non-perturbative Sudakov parameter g1"           , po::value<double>() )
        ("b,blim"            , "Set the blim parameter of the bstar prescription"    , po::value<double>() )
        ("l,bpre"            , "Set the prescription for dealing with the Landau pole in b-space"                               , po::value<int>()    )
        ("q,qtcut"           , "Set fixed IR cut-off for the qt subtraction: qt > qtcut"   , po::value<double>() )
        ("x,xqtcut"          , "Set dynamic IR cut-off for the qt subtraction: qt/m > xqtcut"              , po::value<double>() )
        ("j,taucut"          , "Set taucut for the jettiness subtraction"           , po::value<double>() )
        ("R,kmuren"          , "Set relative renormalization scale"                 , po::value<double>() )
        ("F,kmufac"          , "Set relative factorization scale"                   , po::value<double>() )
        ("Q,kmures"          , "Set relative resummation scale"                     , po::value<double>() )
        ("npff"              , "Set the non perturbative form factor functional form"    , po::value<int>() )
        ("g0"                , "Set the non perturbative CS kernel parameter g0"    , po::value<double>() )
        ("qpar"              , "Set non-perturbative Sudakov parameter q"           , po::value<double>() )
        ("mkcuts"            , "Apply lepton fiducial cuts"                )
        ("nocuts"            , "Full-lepton phase space"                )
        ("cores"             , "Set number of multi-threading cores [integer]"      , po::value<int>()    )
        ("qtbins"            , "Set equidistant binning for mass [N,lo,hi]"            , po::value<string>() )
        ("ybins"             , "Set equidistant binning for qt [N,lo,hi]"              , po::value<string>() )
        ("mbins"             , "Set equidistant binning for y [N,lo,hi]"               , po::value<string>() )
        ("grid"              , "Be verbose to avoid jobs being killed")
        ("n,ofname"          , "Set output file name [filename]"                      , po::value<string>() )
        ("qedON"             , "Turn ON QED corrections"               )          
	("qedOFF"            , "Turn OFF QED corrections"          )                 
        ("qedorder"          , "Set the order of the QED corrections [integer]"      , po::value<int>()    )

    ;
    // Parse
    try {
        args.parse_positional( std::vector<string>({"conf_file"}) );
        args.parse(argc,argv);
    }
    catch (po::OptionException &e){
        printf("%s\n", args.help().c_str());
        printf("Bad arguments: %s \n",e.what());
        throw e;
    }
    // Print help and die
    if (args.count("help")) {
      //throw QuitProgram(args.help().c_str());
      //Exit gently
      cout << args.help() << endl;
      exit (-1);
    }

    // load config file (or default settings)
    readfromfile      ( args["conf_file"].as<string>() );
    bins.readfromfile ( args["conf_file"].as<string>() );

    externalpdf = false;
    
    // NOTE: Command line options are overiding the default and config file settings.
    // verbose
    if (args.count("ofname")) output_filename=args["ofname"].as<string>();
    if (args.count("verbose"))
      {
        verbose=true;
        cubaverbosity=3;
      }
    if (args.count("grid"))
      {
        gridverbose=true;
	redirect=true;
        //cubaverbosity=3;
      }
    if (opts.redirect)
      {
	string logfile = output_filename + ".log";
	freopen(logfile.c_str(),"w",stdout);
      }
    if (opts.output_filename!="") HistoHandler::result_filename = opts.output_filename;
    
    // rseed
    if (args.count("seed")) rseed=args["seed"].as<int>();
    // order
    if (args.count("order")) order=args["order"].as<int>();
    // fixed order
    if (args.count("fixedorder"))
      fixedorder          = true;
    if (args.count("resummation"))
      fixedorder          = false;
    if (args.count("primed"))
      primed          = true;
    if (args.count("unprimed"))
      primed          = false;
    
    // small stat
    if (args.count("small-stat")){
        niterBORN          = 1;
        niterCT            = 1;
        niterVJ            = 1;
        vegasncallsBORN    = (fixedorder ? 1e5 : 1e3);
        vegasncallsCT      = 1e4;
        vegasncallsVJLO    = 1e5;
        vegasncallsVJREAL  = 1e5;
        vegasncallsVJVIRT  = 1e5;
        vegasncallsVJVV    = 1e5;
        vegasncallsVJVR    = 1e5;
        vegasncallsVJRR    = 1e5;
        relaccuracy       = 0.1;
    }

    // proc
    if (args.count("proc")) {
        string val=args["proc"].as<string>(); 
        ToLower(val);
        if        (val == "z0"){ nproc=3;
        } else if (val == "wp"){ nproc=1;
        } else if (val == "wm"){ nproc=2;
        } else if (val == "h0"){ nproc=4; //CJU
        } else {
            throw QuitProgram("Unsupported value of proc: "+val);
        }
    }

    // PDF
    if (args.count("pdfset")) LHAPDFset=args["pdfset"].as<string>();
    if (args.count("pdfvar")) {
        // set default values
        //PDFerrors=false;
        LHAPDFmember=0;
        string val=args["pdfvar"].as<string>();
        ToLower(val);
        //if   (val=="all") { PDFerrors=true; }
	if (IsNumber(val))
	  LHAPDFmember=stod(val);
	else
	  throw QuitProgram("Unsupported value of pdfvar: "+val);
    }

    if (args.count("evmode")) {
      opts.evolmode=args["evmode"].as<int>();
    }
    //Gpar
    if (args.count("gpar")) {
        opts.g1 = args["gpar"].as<double>();
    }
    if (args.count("qpar")) {
        opts.q = args["qpar"].as<double>();
    }
    if (args.count("blim")) {
        opts.blim = args["blim"].as<double>();
    }
    if (args.count("bpre")) {
      opts.bprescription=args["bpre"].as<int>();
    }
    if (args.count("npff")) {
      opts.npff=args["npff"].as<int>();
    }
    if (args.count("g0")) {
        opts.g0 = args["g0"].as<double>();
    }
    if (args.count("qtcut")) {
        opts.qtcut = args["qtcut"].as<double>();
	opts.xqtcut = 0.;
    }
    if (args.count("xqtcut")) {
        opts.xqtcut = args["xqtcut"].as<double>();
	opts.qtcut = 0.;
    }
    if (args.count("taucut")) {
        opts.taucut = args["taucut"].as<double>();
    }
    if (args.count("kmuren")) {
        opts.kmuren = args["kmuren"].as<double>();
    }
    if (args.count("kmufac")) {
        opts.kmufac = args["kmufac"].as<double>();
    }
    if (args.count("kmures")) {
        opts.kmures = args["kmures"].as<double>();
    }

    // Collider
    if (args.count("collider")) {
        string val=args["collider"].as<string>();
        ToLower(val);
        if        (val == "tev1"  ){ sroot=1.80e3; ih1=1; ih2=-1;
        } else if (val == "tev2"  ){ sroot=1.96e3; ih1=1; ih2=-1;
        } else if (val == "lhc5"  ){ sroot=5.00e3; ih1=1; ih2=1;
        } else if (val == "lhc7"  ){ sroot=7.00e3; ih1=1; ih2=1;
        } else if (val == "lhc8"  ){ sroot=8.00e3; ih1=1; ih2=1;
        } else if (val == "lhc13" ){ sroot=13.0e3; ih1=1; ih2=1;
        } else if (val == "lhc14" ){ sroot=14.0e3; ih1=1; ih2=1;
        } else {
            throw QuitProgram("Unsupported value of collider: "+val);
        }
    }

    // term
    if (args.count("term")) {
        // first turn off all terms
      doBORN = doCT = doVJ = doVJREAL = doVJVIRT = doCTAS3 = doVJAS3 = doVJVV = doVJVR = doVJRR = doFPC = false; //doQPC = false ;
        string val=args["term"].as<string>();
        ToUpper(val);
        for (auto piece : Tokenize(val)) {
            if        ( piece == "BORN"    ) { doBORN  = true;
            } else if ( piece == "CT"      ) { doCT    = true;
            } else if ( piece == "VJ"      ) { doVJ    = true; doVJVIRT = true; doVJREAL = true;
            } else if ( piece == "VJVIRT"  ) { doVJ    = true; doVJVIRT = true;
	    } else if ( piece == "VJREAL"  ) { doVJ    = true; doVJREAL = true;
            } else if ( piece == "CTAS3"   ) { doCTAS3 = true;
            } else if ( piece == "VJAS3"   ) { doVJAS3 = true; doVJVV = true; doVJVR = true; doVJRR = true;
	    } else if ( piece == "VJVV"    ) { doVJAS3 = true; doVJVV = true;
	    } else if ( piece == "VJVR"    ) { doVJAS3 = true; doVJVR = true;
	    } else if ( piece == "VJRR"    ) { doVJAS3 = true; doVJRR = true;
	    } else if ( piece == "FPC"     ) { doFPC   = true;
	      //	    } else if ( piece == "QPC"     ) { doQPC   = true;
	    } else if ( piece == "ALL"     ) { doBORN = doCT = doVJ = doVJREAL = doVJVIRT = doCTAS3 = doVJAS3 = doVJVV = doVJVR = doVJRR = doFPC = true; //doQPC = true;
            } else {
                throw QuitProgram("Unsupported value of term : "+piece);
            }
        }
    }

    // cores
    if (args.count("cores")) cubacores=args["cores"].as<int>();

    if (args.count("mkcuts"))
      makecuts          = true;
    if (args.count("nocuts"))
      makecuts          = false;
    // QED and Mixed QCD-QED corrections   
    if (args.count("qedON"))
      qed          = true;
    if (args.count("qedOFF"))
      qed          = false;

    if(args.count("qedorder")) qedorder=args["qedorder"].as<int>();
    
    //binning
    parse_binning("qtbins" , bins.qtbins ,args);
    parse_binning("ybins"  , bins.ybins  ,args);
    parse_binning("mbins"  , bins.mbins  ,args);
    parse_binning("ptlbins"   , bins.ptlbins  ,args);
    parse_binning("etalbins"  , bins.etalbins  ,args);


    // check consistency of settings --> move to after PDF init
    check_consistency();
}

// PROTOTYPE OF VARIABLE REGISTRATION
struct  var_def{
    string name = "";
    void * ptr = NULL;
    char type = '0'; // 0-group name, b-bool, s-string, i-int, d-doubl, v-VecDbl

    string dump(){
        SStream tmp;
        if (type == '0') {
            tmp << "# " << name;
        } else {
            tmp << name << " = ";
            if (type != 'b'){
                bool *val = (bool*) ptr;
                tmp << ((*val) ? "true" : "false");
            }
            if (type != 'd'){
                double *val = (double*) ptr;
                tmp << *val;
            }
            if (type != 'i'){
                int *val = (int*) ptr;
                tmp << *val;
            }
            if (type != 's'){
                string *val = (string*) ptr;
                tmp << *val;
            }
            if (type != 'v'){
                VecDbl *val = (VecDbl*) ptr;
                tmp << "[ ";
                for (double v: *val) tmp << v << " ";
                tmp <<  " ]";
            }
        }
        tmp << "\n";
        return tmp.str();
    }
};
Vec<var_def> var_reg;

void GroupName(string name){ var_def add; add.name=name; var_reg.push_back(add); };
void RegVar(string name, bool   *var){ var_def add; add.name=name; add.ptr=var; add.type='b'; var_reg.push_back(add); };
void RegVar(string name, double *var){ var_def add; add.name=name; add.ptr=var; add.type='d'; var_reg.push_back(add); };
void RegVar(string name, int    *var){ var_def add; add.name=name; add.ptr=var; add.type='i'; var_reg.push_back(add); };
void RegVar(string name, string *var){ var_def add; add.name=name; add.ptr=var; add.type='s'; var_reg.push_back(add); };
void RegVar(string name, VecDbl *var){ var_def add; add.name=name; add.ptr=var; add.type='v'; var_reg.push_back(add); };

void GetVar(string name, bool   &var, InputParser &in){ var = in.GetBool(name)  ; };
void GetVar(string name, double &var, InputParser &in){ var = in.GetNumber(name); };
void GetVar(string name, int    &var, InputParser &in){ var = in.GetNumber(name); };
void GetVar(string name, string &var, InputParser &in){ var = in.GetString(name); };
void GetVar(string name, VecDbl &var, InputParser &in){ var.clear(); in.GetVectorDouble(name,var); };

template <class Tvar>
void addVariable(string name, Tvar& var, InputParser &in){
    GetVar(name,var,in);
    AddToDump(name,&var);
};

// #define GETVARIABLE(var) addVariable( ##var, var, in);

// END OF VARIABLE REGISTRATION



void settings::readfromfile(const string fname){
    //read input settings from file
    InputParser in(fname);
    sroot          = in.GetNumber ( "sroot"      );
    ih1            = in.GetNumber ( "ih1"        );
    ih2            = in.GetNumber ( "ih2"        );
    nuclearpdf     = in.GetBool ( "nuclearpdf"        );
    Z1            = in.GetNumber ( "Z1"        );
    Z2            = in.GetNumber ( "Z2"        );
    A1            = in.GetNumber ( "A1"        );
    A2            = in.GetNumber ( "A2"        );
    nproc          = in.GetNumber ( "nproc"      );
    alphaslha      = in.GetBool   ( "alphaslha"  );
    orderlha       = in.GetBool   ( "orderlha"  );
    murpdf         = in.GetBool   ( "murpdf"  );
    blimlha        = in.GetBool   ( "blimlha"  );
    kmuren         = in.GetNumber ( "kmuren"     );
    kmufac         = in.GetNumber ( "kmufac"     );
    kmures         = in.GetNumber ( "kmures"     );
    fmuren         = in.GetNumber ( "fmuren"     );
    fmufac         = in.GetNumber ( "fmufac"     );
    fmures         = in.GetNumber ( "fmures"     );
    kmuc           = in.GetNumber ( "kmuc"       );
    kmub           = in.GetNumber ( "kmub"       );
    kmut           = in.GetNumber ( "kmut"       );
    accel          = in.GetNumber ( "accel"       );
    npff           = in.GetNumber ( "npff"        );
    g1             = in.GetNumber ( "g1"          );
    g2             = in.GetNumber ( "g2"          );
    g3             = in.GetNumber ( "g3"          );
    g1x            = in.GetNumber ( "g1x"          );
    g1a            = in.GetNumber ( "g1a"          );
    g1b            = in.GetNumber ( "g1b"          );
    gjxnspace      = in.GetBool   ( "gjxnspace"          );
    g1c            = in.GetNumber ( "g1c"          );
    g1m            = in.GetNumber ( "g1m"          );
    g1f            = in.GetNumber ( "g1f"          );
    e              = in.GetNumber ( "e"          );
    q              = in.GetNumber ( "q"          );
    sx             = in.GetNumber ( "sx"          );
    g0             = in.GetNumber ( "g0"          );
    bref           = in.GetNumber ( "bref"        );
    qCS            = in.GetNumber ( "qCS"          );
    Q0             = in.GetNumber ( "Q0"          );
    a2             = in.GetNumber ( "a2"          );
    a2p            = in.GetNumber ( "a2p"         );
    lambda         = in.GetNumber ( "lambda"      );
    g2A		   = in.GetNumber ( "g2A"         );
    g2B		   = in.GetNumber ( "g2B"         );
    NA		   = in.GetNumber ( "NA"          );
    sigmaA	   = in.GetNumber ( "sigmaA"       );
    alphaA	   = in.GetNumber ( "alphaA"       );
    NB		   = in.GetNumber ( "NB"         );
    sigmaB	   = in.GetNumber ( "sigmaB"      );
    alphaB         = in.GetNumber ( "alphaB"      );
    Linf           = in.GetNumber ( "Linf"      );
    L2             = in.GetNumber ( "L2"      );
    L4             = in.GetNumber ( "L4"      );
    linf           = in.GetNumber ( "linf"      );
    l2             = in.GetNumber ( "l2"      );
    l4             = in.GetNumber ( "l4"      );
    flavour_kt     = in.GetBool   ( "flavour_kt"     );
    g1_uv          = in.GetNumber ( "g1_uv"          );
    g1_us          = in.GetNumber ( "g1_us"          );
    g1_dv          = in.GetNumber ( "g1_dv"          );
    g1_ds          = in.GetNumber ( "g1_ds"          );
    g1_ss          = in.GetNumber ( "g1_ss"          );
    g1_ch          = in.GetNumber ( "g1_ch"          );
    g1_bo          = in.GetNumber ( "g1_bo"          );
    g1_gl          = in.GetNumber ( "g1_gl"          );
    order          = in.GetNumber ( "order"          );
    order_sudak    = in.GetNumber ( "order_sudak"    );
    order_hcoef    = in.GetNumber ( "order_hcoef"    );
    order_evol     = in.GetNumber ( "order_evol"    );
    order_expc     = in.GetNumber ( "order_expc"    );
    order_ct       = in.GetNumber ( "order_ct"    );
    order_vjet     = in.GetNumber ( "order_vjet"    );
    qbox           = in.GetBool   ( "qbox"    );
    qtriang        = in.GetBool   ( "qtriang"    );
    A5             = in.GetNumber ( "A5"    );
    a5lc           = in.GetBool ( "a5lc"    );
    H4             = in.GetNumber ( "H4"    );
    gam4NSp         = in.GetNumber ( "gam4NSp"    );
    gam4NSm         = in.GetNumber ( "gam4NSm"    );
    gam4NSv         = in.GetNumber ( "gam4NSv"    );
    gam4SGqq         = in.GetNumber ( "gam4SGqq"    );
    gam4SGqg         = in.GetNumber ( "gam4SGqg"    );
    gam4SGgq         = in.GetNumber ( "gam4SGgq"    );
    gam4SGgg         = in.GetNumber ( "gam4SGgg"    );
    C4qq             = in.GetNumber ( "C4qq"    );
    C4qg             = in.GetNumber ( "C4qg"    );
    C4qqb             = in.GetNumber ( "C4qqb"    );
    C4qqp             = in.GetNumber ( "C4qqp"    );
    C4qqbp             = in.GetNumber ( "C4qqbp"    );
    sgapprox       = in.GetNumber ( "sgapprox"    );
    nsapprox       = in.GetBool ( "nsapprox"    );
    runningwidth   = in.GetBool   ( "runningwidth"   );
    conv2fixw      = in.GetBool   ( "conv2fixw"   );
    rseed          = in.GetNumber ( "rseed"          );
    blim           = in.GetNumber ( "blim"           );
    blim_pdf       = in.GetNumber ( "blim_pdf"       );
    blim_sudakov   = in.GetNumber ( "blim_sudakov"   );
    blim_expc      = in.GetNumber ( "blim_expc"      );
    bstar_pdf      = in.GetBool   ( "bstar_pdf"       );
    bstar_sudakov  = in.GetBool   ( "bstar_sudakov"   );
    bstar_expc     = in.GetBool   ( "bstar_expc"      );
    LHAPDFset      = in.GetString ( "LHAPDFset"      );
    LHAPDFmember   = in.GetNumber ( "LHAPDFmember"   );
    LHAPDFset_beam2      = in.GetString ( "LHAPDFset_beam2"      );
    LHAPDFmember_beam2   = in.GetNumber ( "LHAPDFmember_beam2"   );
    ewscheme       = in.GetNumber ( "ewscheme"       );
    Gf             = in.GetNumber ( "Gf"           );
    zmass          = in.GetNumber ( "zmass"        );
    wmass          = in.GetNumber ( "wmass"        );
    xw             = in.GetNumber ( "xw"           );
    aemmz          = in.GetNumber ( "aemmz"        );
    zwidth         = in.GetNumber ( "zwidth"       );
    wwidth         = in.GetNumber ( "wwidth"       );
    Vud            = in.GetNumber ( "Vud"        );
    Vus            = in.GetNumber ( "Vus"        );
    Vub            = in.GetNumber ( "Vub"        );
    Vcd            = in.GetNumber ( "Vcd"        );
    Vcs            = in.GetNumber ( "Vcs"        );
    Vcb            = in.GetNumber ( "Vcb"        );
    Zuu            = in.GetNumber ( "Zuu"        );
    Zdd            = in.GetNumber ( "Zdd"        );
    Zcc            = in.GetNumber ( "Zcc"        );
    Zss            = in.GetNumber ( "Zss"        );
    Zbb            = in.GetNumber ( "Zbb"        );
    damp           = in.GetBool   ( "damp"          );
    dampk          = in.GetNumber ( "dampk"          );
    dampdelta      = in.GetNumber ( "dampdelta"      );
    dampmode       = in.GetNumber ( "dampmode"       );
    qtlim          = in.GetBool   ( "qtlim"          );
    resscheme      = in.GetNumber ( "resscheme"       );
    qtcutoff       = in.GetNumber ( "qtcutoff"       );
    xqtcutoff      = in.GetNumber ( "xqtcutoff"       );
    mcutoff        = in.GetNumber ( "mcutoff"       );
    mcutoffas3     = in.GetNumber ( "mcutoffas3"    );
    aii            = in.GetNumber ( "aii"       );
    aif            = in.GetNumber ( "aif"       );
    afi            = in.GetNumber ( "afi"       );
    aff            = in.GetNumber ( "aff"       );
    modlog         = in.GetBool ( "modlog"       );
    p              = in.GetNumber ( "p"       );
    xqtcut         = in.GetNumber ( "xqtcut"         );
    qtcut          = in.GetNumber ( "qtcut"          );
    qtfpc          = in.GetNumber ( "qtfpc"          );
    xqtcutas3      = in.GetNumber ( "xqtcutas3"         );
    qtcutas3       = in.GetNumber ( "qtcutas3"          );
    taucut         = in.GetNumber ( "taucut"          );
    ftau           = in.GetNumber ( "ftau"          );
    ptjetmin       = in.GetNumber ( "ptjetmin"        );
    tauboost       = in.GetBool   ( "tauboost"          );
    intDimRes      = in.GetNumber ( "intDimRes"      );
    intDimBorn     = in.GetNumber ( "intDimBorn"     );
    intDimCT       = in.GetNumber ( "intDimCT"       );
    intDimVJ       = in.GetNumber ( "intDimVJ"        );
    intDimFO       = in.GetNumber ( "intDimFO"       );
    intDimFPC      = in.GetNumber ( "intDimFPC"        );
    BORNquad       = in.GetBool   ( "BORNquad"          );
    CTquad         = in.GetBool   ( "CTquad"            );
    VJquad         = in.GetBool   ( "VJquad"            );
    FPCquad        = in.GetBool   ( "FPCquad"            );
    fixedorder         = in.GetBool   ( "fixedorder_only" );
    primed             = in.GetBool   ( "primed" );
    doBORN             = in.GetBool   ( "doBORN"          );
    doCT               = in.GetBool   ( "doCT"            );
    doVJ               = in.GetBool   ( "doVJ"            );
    doVJREAL           = in.GetBool   ( "doVJREAL"        );
    doVJVIRT           = in.GetBool   ( "doVJVIRT"        );
    doCTAS3            = in.GetBool   ( "doCTAS3"        );
    doVJAS3            = in.GetBool   ( "doVJAS3"        );
    doVJVV             = in.GetBool   ( "doVJVV"        );
    doVJVR             = in.GetBool   ( "doVJVR"        );
    doVJRR             = in.GetBool   ( "doVJRR"        );
    doFPC              = in.GetBool   ( "doFPC"           );
    finord             = in.GetBool   ( "finord"            );
    threading          = in.GetNumber ( "threading"       );
    cubaverbosity      = in.GetNumber ( "cubaverbosity"   );
    cubacores          = in.GetNumber ( "cores"       );
    cubanbatch         = in.GetNumber ( "cubanbatch"      );
    niterBORN          = in.GetNumber ( "niterBORN"       );
    niterCT            = in.GetNumber ( "niterCT"         );
    niterVJ            = in.GetNumber ( "niterVJ"         );
    vegasncallsBORN    = in.GetNumber ( "vegasncallsBORN" );
    vegasncallsCT      = in.GetNumber ( "vegasncallsCT"   );
    vegasncallsVJLO    = in.GetNumber ( "vegasncallsVJLO"   );
    vegasncallsVJREAL  = in.GetNumber ( "vegasncallsVJREAL" );
    vegasncallsVJVIRT  = in.GetNumber ( "vegasncallsVJVIRT" );
    vegasncallsVJVV    = in.GetNumber ( "vegasncallsVJVV" );
    vegasncallsVJVR    = in.GetNumber ( "vegasncallsVJVR" );
    vegasncallsVJRR    = in.GetNumber ( "vegasncallsVJRR" );
    vegascollect       = in.GetBool ( "vegascollect" );
    vegascorr          = in.GetBool ( "vegascorr" );
    pcubature          = in.GetBool   ( "pcubature" );
    relaccuracy        = in.GetNumber ( "relaccuracy" );
    absaccuracy        = in.GetNumber ( "absaccuracy" );
    level              = in.GetNumber ( "level" );
    costhmin           = in.GetNumber ( "costhmin"     );
    costhmax           = in.GetNumber ( "costhmax"     );
    xfmin           = in.GetNumber ( "xfmin"     );
    xfmax           = in.GetNumber ( "xfmax"     );
    makecuts           = in.GetBool   ( "makecuts"     );
    qed                = in.GetBool   ( "qed"          );
    qedorder           = in.GetNumber ( "qedorder"          );
    qedorder_sudak    = in.GetNumber ( "qedorder_sudak"    );
    qedorder_hcoef    = in.GetNumber ( "qedorder_hcoef"    );
    qedorder_evol     = in.GetNumber ( "qedorder_evol"    );
    qedorder_expc     = in.GetNumber ( "qedorder_expc"    );
    qedorder_ct       = in.GetNumber ( "qedorder_ct"    );
    qedorder_vjet     = in.GetNumber ( "qedorder_vjet"    );
    alpha0            = in.GetBool ( "alpha0"    );
    lptcut             = in.GetNumber ( "lptcut"         );
    lycut              = in.GetNumber ( "lycut"          );
    mtcut              = in.GetNumber ( "mtcut"          );
    etmisscut          = in.GetNumber ( "etmisscut"      );
    lepptcut            = in.GetNumber ( "lepptcut"        );
    lepycut             = in.GetNumber ( "lepycut"         );
    alpptcut            = in.GetNumber ( "alpptcut"        );
    alpycut             = in.GetNumber ( "alpycut"          );
    lcptcut            = in.GetNumber ( "lcptcut"        );
    lcymin             = in.GetNumber ( "lcymin"         );
    lcymax             = in.GetNumber ( "lcymax"         );
    lfptcut            = in.GetNumber ( "lfptcut"        );
    lfymin             = in.GetNumber ( "lfymin"         );
    lfymax             = in.GetNumber ( "lfymax"         );
    cthCSmin           = in.GetNumber ( "cthCSmin"         );
    cthCSmax           = in.GetNumber ( "cthCSmax"         );
    cubaint            = in.GetBool   ( "cubaint"         ); //true    # integration with     Cuba       Suave
    trapezint          = in.GetBool   ( "trapezint"       ); //false   # trapezoidal rule     for        the     phi_lep     integration     and         semi-analytical for         costh
    quadint            = in.GetBool   ( "quadint"         ); //false   # quadrature  rule     for        the     phi_lep     integration     and         semi-analytical for         costh
    suavepoints        = in.GetNumber ( "suavepoints"     ); //1000000 # number      of       points     for     suave       integration,    newpoints   is              set         to    suavepoints/10;
    nphitrape          = in.GetNumber ( "nphitrape"       ); //1000    # number      of       steps      for     trapezoidal rule            of          phi_lep         integration
    phiintervals       = in.GetNumber ( "phiintervals"    );
    phirule            = in.GetNumber ( "phirule"         );
    ncstart            = in.GetNumber ( "ncstart"         ); //1000    # starting    sampling for        the     costh       semi-analytical integration (common         settings    for   the             trapezoidal and quadrature rules)
    qtrec_naive        = in.GetBool   ( "qtrec_naive"     ); //false
    qtrec_cs           = in.GetBool   ( "qtrec_cs"        ); //false
    qtrec_kt0          = in.GetBool   ( "qtrec_kt0"       ); //true
    timeprofile        = in.GetBool   ( "timeprofile"     ); //false   # debug       and      time       profile resummation integration
    verbose            = in.GetBool   ( "verbose"         ); //false   # debug       and      time       profile costh       phi_lep         integration
    gridverbose        = in.GetBool   ( "gridverbose"     ); //false   # debug       and      time       profile costh       phi_lep         integration
    texttable          = in.GetBool   ( "texttable"       ); //
    redirect           = in.GetBool   ( "redirect"       ); //
    unicode            = in.GetBool   ( "unicode"         ); //
    silent             = in.GetBool   ( "silent"         ); //
    makehistos         = in.GetBool   ( "makehistos"      ); //
    useGamma           = in.GetBool ( "useGamma" );//
    useZ               = in.GetBool ( "useZ" );//
    //    PDFerrors           = in.GetBool ( "PDFerrors" );//
    opts_.approxpdf_    = in.GetNumber ( "opts_approxpdf" ); //0
    opts_.pdfintervals_ = in.GetNumber ( "opts_pdfintervals" ); //100
    pdfrule            = in.GetNumber ( "pdfrule" );
    evolmode           = in.GetNumber  ("evolmode");
    iterative          = in.GetBool  ("iterative");
    hfmode             = in.GetNumber  ("hfmode");
    smallx             = in.GetBool  ("smallx");
    nsub               = in.GetNumber ("nsub");
    //mufevol            = in.GetBool  ("mufevol");
    mufvar             = in.GetBool  ("mufvar");
    commute            = in.GetBool  ("commute");
    expc               = in.GetNumber  ("expc");
    ntaylor            = in.GetNumber  ("ntaylor");
    expcreg            = in.GetNumber  ("expcreg");
    sumlogs            = in.GetBool  ("sumlogs");
    numsud             = in.GetBool   ("numsud");
    asrgkt             = in.GetBool   ("asrgkt");
    deltaismb          = in.GetBool   ("deltaismb");
    deltafsmb          = in.GetBool   ("deltafsmb");
    deltaismc          = in.GetBool   ("deltaismc");
    deltafsmc          = in.GetBool   ("deltafsmc");
    vfncqg             = in.GetBool   ("vfncqg");
    asrun              = in.GetNumber ("asrun");
    cmw                = in.GetBool   ("cmw");
    cmwfo              = in.GetBool   ("cmwfo");
    numexpc            = in.GetBool   ("numexpc");
    bprescription      = in.GetNumber   ("bprescription");
    fbstar             = in.GetNumber   ("fbstar");
    bintaccuracy       = in.GetNumber ( "bintaccuracy" );
    phibr              = in.GetNumber ( "phibr" );
    bcf                = in.GetNumber ( "bcf" );
    mellintr           = in.GetNumber ( "mellintr" );
    mellininv          = in.GetNumber ( "mellininv" );
    mellinintervals    = in.GetNumber ( "mellinintervals" );
    mellinrule         = in.GetNumber ( "mellinrule" );
    zmax               = in.GetNumber ( "zmax" );
    ncycle             = in.GetNumber ( "ncycle" );
    cpoint             = in.GetNumber ( "cpoint" );
    cshift             = in.GetNumber ( "cshift" );
    phi                = in.GetNumber ( "phi" );
    mellincores        = in.GetNumber ( "mellincores" );
    mellin1d           = in.GetBool   ( "mellin1d" );
    melup              = in.GetNumber ( "melup" );
    xspace             = in.GetBool   ( "xspace" );
    ctmellin           = in.GetBool   ( "ctmellin" );
    mrule              = in.GetNumber ( "mrule" );
    yintervals         = in.GetNumber ( "yintervals" );
    yrule              = in.GetNumber ( "yrule" );
    qtintervals        = in.GetNumber ( "qtintervals" );
    qtrule             = in.GetNumber ( "qtrule" );
    abintervals        = in.GetNumber ( "abintervals" );
    abrule             = in.GetNumber ( "abrule" );
    vjphirule          = in.GetNumber ( "vjphirule" );
    zrule              = in.GetNumber ( "zrule" );
    xrule              = in.GetNumber ( "xrule" );
    ptbinwidth         = in.GetBool ( "ptbinwidth" );
    ybinwidth          = in.GetBool ( "ybinwidth" );
    mbinwidth          = in.GetBool ( "mbinwidth" );
    dsdqt2             = in.GetBool ( "dsdqt2" );
    dsdxf              = in.GetBool ( "dsdxf" );
    edsdp3             = in.GetBool ( "edsdp3" );
    m3dsdm             = in.GetBool ( "m3dsdm" );
    exactkin           = in.GetBool ( "exactkin" );
    force_binsampling  = in.GetBool ( "force_binsampling" );
    lepbounds          = in.GetBool ( "lepbounds" );
    helicity           = in.GetNumber ( "helicity" );
    h2d                = in.GetBool ( "h2d" );
    h3d                = in.GetBool ( "h2d" );
    hetal              = in.GetBool ( "hetal" );
    hptl               = in.GetBool ( "hptl" );
    hptnu              = in.GetBool ( "hptnu" );
    hmt                = in.GetBool ( "hmt" );
    haiqt              = in.GetBool ( "haiqt" );
    haiy               = in.GetBool ( "haiy" );
    haim               = in.GetBool ( "haim" );
    output_filename    = in.GetString ( "output_filename" );
    calc_type          = in.GetString ( "calc_type" );
    decay_type         = in.GetString ( "decay_type" );
    removebr           = in.GetBool   ( "removebr" );
    hmass          = in.GetNumber ( "hmass"        );
    zmass          = in.GetNumber ( "zmass"        );
    hwidth         = in.GetNumber ( "hwidth"       );
    mtop               = in.GetNumber ( "mtop" );
    mbottom            = in.GetNumber ( "mbottom" );
    mcharm             = in.GetNumber ( "mcharm" );
    mtau               = in.GetNumber ( "mtau" );
    phpt_cut_h         = in.GetNumber ( "phpt_cut_h" );
    phpt_cut_s         = in.GetNumber ( "phpt_cut_s" );
    phycut             = in.GetNumber ( "phycut" );
    apply_crack           = in.GetBool   ( "apply_crack" );
    crack_eta_max      = in.GetNumber ( "crack_eta_max" );
    crack_eta_min      = in.GetNumber ( "crack_eta_min" );
    ph_iso_R0          = in.GetNumber ( "ph_iso_R0"    );
    ph_iso_E0          = in.GetNumber ( "ph_iso_E0"   );
    n_isol             = in.GetNumber ( "n_isol"      );
//    Rgamgam            = in.GetNumber ( "Rgamgam" );
    vegasncallsqqbarLO = in.GetNumber ( "vegasncallsqqbarLO" );
    vegasncallsggNLO   = in.GetNumber ( "vegasncallsggNLO" );
    LO_qg_ON           = in.GetBool   ( "LO_qg_ON" );
    LO_qqbar_ON        = in.GetBool   ( "LO_qqbar_ON" );
    gg_only            = in.GetBool   ( "gg_only" );
    qg_only            = in.GetBool   ( "qg_only" );
    withsmearing       = in.GetBool   ( "withsmearing" );
    SigmaSmearing      = in.GetNumber ( "SigmaSmearing" );
    prodmtinf          = in.GetBool   ( "prodmtinf" );
    eta_jet            = in.GetNumber ( "eta_jet" );
    Dixon_eta_jet_cut  = in.GetBool   ( "Dixon_eta_jet_cut" );
    return ;
}

void settings::check_consistency(){

    // additional conditions
    //if (order != 0 && order != 1 && order != 2 && order != 3)
    //throw invalid_argument("Invalid order, please select 0 (LL) 1 (NLL) 2 (NNLL) or 3 (NNNLL)");

    if (nproc != 1 && nproc != 2 && nproc != 3 && nproc != 4)
      throw invalid_argument("Wrong process, please select nproc = 1 (W+), 2 (W-), 3(Z), or 4(H)");

    if (order_sudak < 0)
      order_sudak = order;

    if (order_hcoef < 0)
      if (primed)
	order_hcoef = order;
      else
	order_hcoef = max(0,order-1);

    if (order_evol < 0)
      order_evol = order;

    if (order_expc < 0)
      order_expc = order;

    if (order_ct < 0)
      {
	if (primed)
	  order_ct = order;
	else
	  order_ct = max(0,order-1);
	//order_ct = min(2,order_ct);
      }

    if (order_vjet < 0)
      {
	if (primed)
	  order_vjet = order;
	else
	  order_vjet = max(0,order-1);
	order_vjet = min(3,order_vjet);
      }

    if (opts.externalpdf)
      if (opts.orderlha)
	//order_evol = max(order_evol,pdf::order+1);
	order_evol = pdf::order+1;
    
    //order = max(order,order_sudak);
    //order = max(order,order_hcoef);
    //order = max(order,order_evol);
    //order = max(order,order_expc);
    
    if (fixedorder == false || !opts.makecuts || order == 0)
      doFPC = false;

    //    if (fixedorder == false || order == 0)
    //      doQPC = false;
    
    if (order_ct == 0)
      doCT = false;

    if (order_ct < 3)
      doCTAS3 = false;
    
    if (order_vjet == 0)
      doVJ = false;

    if (order_vjet < 2)
      {
	doVJREAL = false;
	doVJVIRT = false;
      }

    if (order_vjet < 3)
      {
	doVJAS3 = false;
	doVJVV  = false;
	doVJVR  = false;
	doVJRR  = false;
      }
    
    if (asrun != 0 && asrun != 1 && asrun != 2 && asrun != 3 && asrun != 4 && asrun != 5)
      {
	cout << "Invalid option for asrun, should be between 0 and 5" << endl;
	exit (-1);
      }	

    //numerical or iterative alphas or non truncated Sudakov require numerical integration of the sudakov
    if (asrgkt || order >= 4)
      {
	numsud = true;
	numexpc = true;
      }
    
    //In fixed order mode, a_param must be one
    if (fixedorder == true)
      {
	//cout << "Asked for fixed order only predictions, enforce kmures = 1.0" << endl;
	kmures = 1.0;
	fmures = 1;
      }

    if (npff != 0 && npff != 1 && npff != 2 && npff != 3 && npff != 4 && npff != 5 && npff != 6 && npff != 7 && npff != 8  && npff != 9  && npff != 10)
      {
	cout << "Invalid option for npff, should be between 0 and 10" << endl;
	exit (-1);
      }	

    //convert g at the Z mass into a2 and a2p, assuming g1 = 0
    //if (npff == 3)
    //  {
    //	a2 = g1/log(zmass)/2.; //(g-g1)/log(91.)/2; a2 = g2/2
    //	a2p = a2;//a2 - g1; a2p = g2/2 - g1
    //  }

    ////Force Gaussian NP parameters to be positive
    //g1 = max(0.,g1);
    //g2 = max(0.,g2);
    //g3 = max(0.,g3);
    //g1a = max(0.,g1a);
    //g1b = max(0.,g1b);
    
    //If default values of blim for pdfs, sudakov and expc are not set, take blim as default
    if (blim_pdf == 0)
      blim_pdf = blim;

    if (blim_sudakov == 0)
      blim_sudakov = blim;

    if (blim_expc == 0)
      blim_expc = blim;

    //local bstar requires numerical integration of the sudakov
    if (bprescription == 4)
      {
	numsud = true;
	numexpc = true;
      }

    //QED corrections not implemented in hcoeff::calcb() required by some expc modes
    if (opts.expc >= 8 && opts.expc <= 10)
      qed = false;
    
    //if (bcf < 0 || bcf > 1)
    if (bcf > 1)
      {
	cout << "bcf should be between 0 and 1 (or negative)" << endl;
	exit (-1);
      }
      
    if (evolmode > 5 || evolmode < 0)
      {
	cout << "wrong value for evolmode: available evolmodes: 0,1,2,3,4,5" << endl;
	exit (-1);
      }

    //reset evolmode to zero at leading log
    if (evolmode == 1 && order_evol == 0)
      evolmode = 0;

    //use evolmode 1 at N3LL (evolmode 0 not implemented)
    if (evolmode == 0 && order_evol >= 3)
      evolmode = 1;

    //use evolmode 5 at N4LL (evolmode 0 and 1 not implemented)
    if (evolmode <= 1 && order_evol >= 4)
      evolmode = 5;

#ifndef USEHELLX
  if (opts.evolmode == 4 || opts.smallx)
    {
      cout << "For evolmode = 4 and smallx HELLx must be linked" << endl;
      exit (-1);
    }
#endif
	
    if (qedorder_sudak < 0)
      qedorder_sudak = qedorder;

    if (qedorder_hcoef < 0)
      qedorder_hcoef = qedorder;

    if (qedorder_evol < 0)
      qedorder_evol = qedorder;

    if (qedorder_expc < 0)
      qedorder_expc = qedorder;

    if (qedorder_ct < 0)
      qedorder_ct = qedorder;

    if (qedorder_vjet < 0)
      qedorder_vjet = qedorder;
    
    if (qedorder_evol >= 1 && pdf::photon && opts.doBORN && fixedorder == false && (evolmode == 1 || evolmode == 3))
      {
	cout << "Cannot perform QED evolution with evolmode 1 or 3" << endl;
	exit (-1);
      }
    
    //Force bstar prescription for PDFs for evolmode 2
    if (evolmode == 2)
	bstar_pdf = true;
      
    if (blimlha && evolmode == 2)
      if (pdf::qmin > 0)
	blim_pdf = min(resconst::b0/pdf::qmin, blim_pdf);

    
    //if (fmufac > 0 && evolmode < 3 && order > 0 && fixedorder == false)
    //  {
    //	//cannot use a dynamic muren, mufac, when the PDFs are converted from x-space to N-space at the factorisation scale
    //	cout << "At NLL and NNLL mufac = mll is possible only with evolmode = 3 or 4" << endl;
    //	exit (-1);
    //  }

    //minimal b prescription in the complex plain available only for evolmode 0
    //real axix minimal presxription works only at LL

    //with modlog = false can run only the resummation term

    //cannot take PDFs in xspace when Mellin convolution is needed
    if (xspace && (order_hcoef > 0 || order_expc > 0))
      xspace = false;
    
//    if (PDFerrors == true && LHAPDFmember != 0)
//      {
//	cout << "Asked for PDFerrors, enforce LHAPDFmember  = 0" << endl;
//	LHAPDFmember = 0;
//      }
    PDFerrors = false;
    
    if (qtcut <= 0 && xqtcut <= 0)
      {
	cout << "At least one between qtcut and xqtcut must be > 0" << endl;
	exit (-1);
      }

    if (qtcutas3 <= 0 && xqtcutas3 <= 0)
      {
	cout << "At least one between qtcutas3 and xqtcutas3 must be > 0" << endl;
	exit (-1);
      }
    
    //Automatic selector of integration type
    if (intDimBorn < 0)
      if (BORNquad)
	{
	  intDimBorn = 2;
	  if (order >= 3)
	    intDimBorn = 1;
	}
      else
	intDimBorn = 4; //6;

    if (intDimRes < 0)
      if (BORNquad)
	if (opts.makecuts)
	  intDimRes = 2;
	else
	  intDimRes = 1; // -> Check this works also when not fully integrated in rapidity!!!
      else
	intDimRes = 4;

    //Add some safety cuts on mll and yll, when makecuts is true
    //if (opts.makecuts)
    //  {
    //	if (opts.proc == 1 || opts.proc == 2)
    //	  mmin = max(mmin,opts.mtcut);
    // also, mmin should not be zero!
    //	ymax = max(ymax,eta1,eta2);
    //	ymin = min(ymin,-eta1,-eta2);
    //	    
    //  }
    
    //Determine wether the integration range in rapidity is larger than [-ymax,ymax]
    bool fullrap = true;
    double ylim = log(opts.sroot/bins.mbins.front());
    if (bins.ybins.front() >= -ylim || bins.ybins.back() <= ylim)
      fullrap = false;
    if (bins.ybins.size() > 2)
      fullrap = false;
    if (dsdxf)
      fullrap = false;
    
    if (bins.ybins.size() > 2 && dsdxf)
      {
	cout << "dsdxf option should be used only with one large rapidity bin" << endl;
	exit (-1);
      }
    
    //Cannot use mellin1d option if yrange is below ymax or makecuts is true
    if (opts.makecuts || !fullrap)
      mellin1d = false;

    //Need to compute both beams for fixed target experiments with nuclear PDFs
    if (nuclearpdf) 
      mellin1d = false;
    
    if (intDimCT < 0)
      if (CTquad)
	if (opts.makecuts || !fullrap || !mellin1d || order_ct >= 3)
	  intDimCT = 2;
	else
	  intDimCT = 1;
      else
	intDimCT = 6;

    if (intDimVJ < 0)
      if (VJquad)
	if (opts.makecuts)
	  intDimVJ = 5;
	else
	  intDimVJ = 3;
      else
	intDimVJ = 7;

    if (intDimFPC < 0)
      if (FPCquad)
	intDimFPC = 2;
      else
	intDimFPC = 8;
    
    //cannot update mellin support points at every phase-space point when rapidity integrals are calculated numerically
    if (makecuts && (intDimRes <= 2 || intDimBorn == 1))
      if (melup == 2)
	melup = 1;
          
    if (makecuts && doVJ && intDimVJ < 7 && order_vjet >= 2)
      {

	cout << "cannot perform quadrature integration for V+jet at NNLO with cuts on the leptons" << endl;
	exit (-1);
      }

    // resummation term integration dimension
    if (intDimRes<4 && intDimRes>=0){
        resint0d = (intDimRes == 0);
        resint1d = (intDimRes == 1);
        resint2d = (intDimRes == 2);
        resint3d = (intDimRes == 3);
        resintvegas = false;
    } else {
        resint2d = false;
        resint3d = false;
        resintvegas = true;
    }

    //Cannot use mellin1d option for vegas integration
    if (doBORN && resintvegas)
      mellin1d = false;

    // born term integration dimension
    if (intDimBorn < 4 && intDimBorn>0){
      bornint1d      = (intDimBorn == 1);
      bornint2d      = (intDimBorn == 2);
      bornintvegas4d = false;
      bornintvegas6d = false;
    } else {
      bornint1d      = false;
      bornint2d      = false;
      bornintvegas4d = (intDimBorn == 4);
      bornintvegas6d = (intDimBorn >  5);
    }

    // counter term integration dimension
    if (intDimCT<4 && intDimCT>=0){
        ctint0d = (intDimCT == 0);
        ctint1d = (intDimCT == 1);
        ctint2d = (intDimCT == 2);
        ctint3d = (intDimCT == 3);
        ctintvegas6d = false;
	ctintvegas8d = false;
    } else {
        ctint2d = false;
        ctint3d = false;
        ctintvegas6d = (intDimCT <= 6);
	ctintvegas8d = (intDimCT > 6);
    }

    // counter term integration dimension
    if (intDimFO < 4 && intDimFO >= 0)
      {
        foint0d = (intDimFO == 0);
        foint3d = (intDimFO >= 1);
      }
    else
      {
        foint0d = false;
        foint3d = true;
      }
    
    //Cannot use mellin1d option for vegas integration
    if (doCT && (ctintvegas8d || ctintvegas6d))
      mellin1d = false;
    
    // V+J integration dimension
    if (intDimVJ < 7 && intDimVJ >= 0)
      {
        vjint0d = (intDimVJ == 0);
        vjint3d = (intDimVJ == 3);
        vjint5d = (intDimVJ == 5);
        vjintvegas7d = false;
      }
    else
      {
        vjint0d = false;
        vjint3d = false;
	vjint5d = false;
        vjintvegas7d = true;
      }

    // FPC term integration dimension
    if (intDimFPC<4 && intDimFPC>0)
      {
	fpcint2d = (intDimFPC == 2);
	fpcintvegas8d = false;
      }
    else
      {
        fpcint2d = false;
	fpcintvegas8d = (intDimFPC >= 4);
    }

    // Finite order
    finord = (finord && doCT && doVJ && !makecuts && !fixedorder && (ctint0d||ctint1d||ctint2d||ctint3d));

    if (finord && ctmellin)
      mellin1d = false;
    
    if (doBORN && mellin1d && (!(resint2d || resint1d || resint0d) || makecuts))
      {
	cout << "mellin1d option is possible only for 0d, 1d or 2d integration of resummed piece, no cuts on leptons, and integration between [-ymin,ymax]" << endl;
	exit (-1);
      }

    if (mellinintervals*mellinrule >= 512) // >= ndim //(from pegasus.h)
      {
	cout << "mellinintervals*mellinrule  should be less than 512 " << endl;
	exit(-1);
      }
    
    //adjust zmax and cpoint between mellin1d and mellin2d !!! --> This is very empirical, why should zmax and cpoint be different between 1d and 2d???
    //The scaling from 1d to 2d should account for the change in the inverse transform of [F(N)]^2 -> F(N)
    if (!mellin1d)
      {
	//zmax   = (zmax+1)*2-1;
	//cpoint = (cpoint+1)*2-1;
      }
    
    //Check that no rule is bigger than the maximum value computed
    if (pdfrule       > GRNMAX
	|| phirule    > GRNMAX
	|| mellinrule > GRNMAX
	|| mrule      > GRNMAX
	|| yrule      > GRNMAX
	|| qtrule     > GRNMAX
	|| abrule     > GRNMAX
	|| vjphirule  > GRNMAX
	|| zrule      > GRNMAX
	|| xrule      > GRNMAX)
      {
	cout << "Check that no quadrature rule is larger than GRNMAX = " << GRNMAX << endl;
	exit (-1);
      }

    if (mellininv > 1 || mellininv < 0)
      {
	cout << "Not valid option for mellininv (should be 0 or 1) " << endl;
	exit (-1);
      }
    
    if (makecuts && vjint3d)
      {
	cout << "Required cuts on the final state leptons, enforce 5D integration for V+J fixed order cross section" << endl;
	vjint3d = false;
        vjint5d = true;
      }

    if (makecuts && helicity >= 0)
      {
	cout << "Required cuts on the final state leptons, cannot calculate helicity cross sections, enforce helicity = -1" << endl;
	helicity = -1;
      }
    
    if (opts_.approxpdf_ == 1)
      {
	cout << "DYRES-style approximate PDF requested, enforce vegas integration for the resummed cross section" << endl;
	resint2d = false;
	resint3d = false;
	resintvegas = true;
      }

    //disable ptbinwidth, mbinwidth, ybinwidth options if there non-integrable observables
    
    // -- binning
    // check bins size
    if ( bins.qtbins .size() < 2) throw QuitProgram("Option `qtbins` needs at least 2 items ");
    if ( bins.ybins  .size() < 2) throw QuitProgram("Option `ybins`  needs at least 2 items ");
    if ( bins.mbins  .size() < 2) throw QuitProgram("Option `mbins`  needs at least 2 items ");
    // check sorting
    sort( bins.qtbins .begin () , bins.qtbins .end () ) ;
    sort( bins.ybins  .begin () , bins.ybins  .end () ) ;
    sort( bins.mbins  .begin () , bins.mbins  .end () ) ;
    // set histogram bins
    bins.hist_qt_bins = bins.qtbins ;
    bins.hist_y_bins  = bins.ybins ;
    bins.hist_m_bins  = bins.mbins ;
    // integration boundaries
    phasespace::setbounds(
            bins. mbins  .front() ,
            bins. mbins  .back()  ,
            bins. qtbins .front() ,
            bins. qtbins .back()  ,
            bins. ybins  .front() ,
            bins. ybins  .back()
            );

    return ;
}


void settings::dumpAll(){
    printf("==Listing settings==\n");
    bool print_inputs = true;
    bool print_process_inputs = true;
    bool print_masses = true;

    if (print_process_inputs) {
        printf("Input process settings:\n");
        dumpD ( "sroot       ",  energy_      . sroot_      ) ;
        dumpI ( "ih1         ",  density_     . ih1_        ) ;
        dumpI ( "ih2         ",  density_     . ih2_        ) ;
        dumpI ( "Z1         ",  Z1        ) ;
        dumpI ( "Z2         ",  Z2        ) ;
        dumpI ( "A1         ",  A1        ) ;
        dumpI ( "A2         ",  A2        ) ;
        dumpI ( "nproc       ",  nproc_       . nproc_      ) ;
    }

    if (print_inputs) {
        printf("Input settings:\n");
        dumpD ("g              ",  g_param_     . g_param_    ) ;
	dumpI ("npff           ",npff               );
	dumpD ("g1             ",g1               );
	dumpD ("g2             ",g2               );
	dumpD ("g3             ",g3               );
	dumpD ("g1x            ",g1x              );
	dumpD ("g1a            ",g1a              );
	dumpD ("g1b            ",g1b              );
	dumpB ("gjxnspace      ",gjxnspace        );
	dumpD ("g1c            ",g1c              );
	dumpD ("g1m            ",g1m              );
	dumpD ("g1f            ",g1f              );
	dumpD ("e              ",e                );
	dumpD ("q              ",q               );
	dumpD ("sx             ",sx               );
	dumpD ("g0             ",g0               );
	dumpD ("bref           ",bref              );
	dumpD ("qCS            ",qCS               );
	dumpD ("Q0             ",Q0               );
	dumpD ("a2             ",a2               );
	dumpD ("a2p            ",a2p              );
        dumpD ("lambda         ", lambda         );
        dumpD ("g2A	       ", g2A	       );
        dumpD ("g2B	       ", g2B	       );
        dumpD ("NA	       ", NA	       );
        dumpD ("sigmaA	       ", sigmaA	       );
        dumpD ("alphaA	       ", alphaA	       );
        dumpD ("NB	       ", NB	       );
        dumpD ("sigmaB	       ", sigmaB	       );
        dumpD ("alphaB         ", alphaB         );
        dumpD ("Linf           ", Linf         );
        dumpD ("L2             ", L2         );
        dumpD ("L4             ", L4         );
        dumpD ("linf           ", linf         );
        dumpD ("l2             ", l2         );
        dumpD ("l4             ", l4         );
	dumpB ("flavour_kt     ",flavour_kt       );
	dumpD ("g1_uv          ",g1_uv            );
	dumpD ("g1_us          ",g1_us            );
	dumpD ("g1_dv          ",g1_dv            );
	dumpD ("g1_ds          ",g1_ds            );
	dumpD ("g1_ss          ",g1_ss            );
	dumpD ("g1_ch          ",g1_ch            );
	dumpD ("g1_bo          ",g1_bo            );
	dumpD ("g1_gl          ",g1_gl            );
	dumpI ( "order       ",  nnlo_        . order_      ) ;
	dumpI ( "order_sudak ",  order_sudak      ) ;
	dumpI ( "order_hcoef ",  order_hcoef      ) ;
	dumpI ( "order_evol ",   order_evol      ) ;
	dumpI ( "order_expc ",   order_expc      ) ;
	dumpI ( "order_ct ",     order_ct      ) ;
	dumpI ( "order_vjet ",   order_vjet      ) ;
	dumpB ( "qbox ",   qbox      ) ;
	dumpB ( "qtriang ",   qtriang      ) ;
	dumpD ( "A5 ",   A5      ) ;
	dumpB ( "a5lc ",   a5lc      ) ;
	dumpD ( "H4 ",   H4      ) ;
	dumpD ( "gam4NSp ",   gam4NSp      ) ;
	dumpD ( "gam4NSm ",   gam4NSm      ) ;
	dumpD ( "gam4NSv ",   gam4NSv      ) ;
	dumpD ( "gam4SGqq ",   gam4SGqq      ) ;
	dumpD ( "gam4SGqg ",   gam4SGqg      ) ;
	dumpD ( "gam4SGgq ",   gam4SGgq      ) ;
	dumpD ( "gam4SGgg ",   gam4SGgg      ) ;
	dumpD ( "C4qq ",   C4qq      ) ;
	dumpD ( "C4qg ",   C4qg      ) ;
	dumpD ( "C4qqb ",   C4qqb      ) ;
	dumpD ( "C4qqp ",   C4qqp      ) ;
	dumpD ( "C4qqbp ",   C4qqbp      ) ;
	dumpI ( "sgapprox ", sgapprox      ) ;
	dumpB ( "nsapprox ", nsapprox      ) ;
        dumpB ( "alphaslha   ",  alphaslha                  ) ;
        dumpB ( "orderlha   ",  orderlha                  ) ;
        dumpB ( "murpdf   ",  murpdf                  ) ;
        dumpB ( "blimlha   ",  blimlha                  ) ;
        dumpD ( "kmuren      ",  kmuren                     ) ;
        dumpD ( "kmufac      ",  kmufac                     ) ;
	dumpD ( "kmures       ",  kmures   ) ;
        dumpI ( "fmuren      ",  fmuren                     ) ;
        dumpI ( "fmufac      ",  fmufac                     ) ;
	dumpI ( "fmures       ",  fmures   ) ;
        dumpD ( "kmuc      ",  kmuc                     ) ;
        dumpD ( "kmub      ",  kmub                     ) ;
	dumpD ( "kmut      ",  kmut   ) ;
	dumpI ( "accel      ",  accel   ) ;
        dumpD( "blim              ",  blim    ) ;
        dumpD( "blim_pdf               ",  blim_pdf    ) ;
        dumpD( "blim_sudakov           ",  blim_sudakov    ) ;
        dumpD( "blim_expc              ",  blim_expc    ) ;
        dumpB( "bstar_pdf              ",  bstar_pdf    ) ;
        dumpB( "bstar_sudakov          ",  bstar_sudakov    ) ;
        dumpB( "bstar_expc             ",  bstar_expc    ) ;
        dumpB("nuclearpdf          ", nuclearpdf           );
        dumpS("LHAPDFset          ", LHAPDFset           );
        dumpI("LHAPDFmember       ", LHAPDFmember        );
        dumpI("LHAPDFmember_beam2       ", LHAPDFmember_beam2        );
        dumpS("LHAPDFset_beam2          ", LHAPDFset_beam2           );
        dumpI("rseed              ", rseed               );
	dumpI("ewscheme           ", ewscheme );
	dumpD("Gf"                 , Gf);
	dumpD("zmass"              , zmass);
	dumpD("wmass"              , wmass   );
	dumpD("xw"                 , xw);
	dumpD("aemmz"              , aemmz);
	dumpD("zwidth"             , zwidth);
	dumpD("wwidth"             , wwidth);
	dumpB("runningwidth"       , runningwidth);
	dumpB("conv2fixw"       , conv2fixw);
	dumpD( "Vud",        Vud);
	dumpD( "Vus",        Vus);
	dumpD( "Vub",        Vub);
	dumpD( "Vcd",        Vcd);
	dumpD( "Vcs",        Vcs);
	dumpD( "Vcb",        Vcb);
	dumpD( "Zuu",        Zuu);
	dumpD( "Zdd",        Zdd);
	dumpD( "Zss",        Zss);
	dumpD( "Zcc",        Zcc);
	dumpD( "Zbb",        Zbb);
        //dumpD("ylow               ", ylow                );
        //dumpD("yhigh              ", yhigh               );
        //dumpD("mlow               ", mlow                );
        //dumpD("mhigh              ", mhigh               );
        dumpB("damp",             damp       );
	dumpD("dampk",             dampk       );
        dumpD("dampdelta",     dampdelta      );
        dumpI("dampmode",       dampmode     );
        dumpB("qtlim",          qtlim       );
        dumpI("resscheme",      resscheme     );
	dumpD("qtcutoff",       qtcutoff     );
	dumpD("xqtcutoff",      xqtcutoff     );
	dumpD("mcutoff",        mcutoff     );
	dumpD("mcutoffas3",     mcutoffas3     );
	dumpD("aii",        aii     );
	dumpD("aif",        aif     );
	dumpD("afi",        afi     );
	dumpD("aff",        aff     );
	dumpB("modlog",         modlog     );
    dumpI("p",         p     );
	dumpD("xqtcut",       xqtcut     );
	dumpD("qtcut",       qtcut     );
	dumpD("qtfpc",       qtfpc     );
	dumpD("xqtcutas3",       xqtcutas3     );
	dumpD("qtcutas3",       qtcutas3     );
	dumpD("taucut",      taucut     );
	dumpI("ftau",        ftau     );
	dumpD("ptjetmin",      ptjetmin     );
        dumpB("tauboost          ", tauboost            );
        dumpB("useGamma          ", useGamma            );
        dumpB("useZ              ", useZ            );
        dumpI("intDimRes         ", intDimRes           );
        dumpB("resint0d          ", resint0d            );
        dumpB("resint1d          ", resint1d            );
        dumpB("resint2d          ", resint2d            );
        dumpB("resint3d          ", resint3d            );
        dumpB("resintvegas       ", resintvegas         );
        dumpI("intDimBorn        ", intDimBorn          );
        dumpB("bornint1d         ", bornint1d           );
        dumpB("bornint2d         ", bornint2d           );
        dumpB("bornintvegas4d    ", bornintvegas4d      );
        dumpB("bornintvegas6d    ", bornintvegas6d      );
        dumpI("intDimCT          ", intDimCT            );
        dumpB("ctint2d           ", ctint2d             );
        dumpB("ctint3d           ", ctint3d             );
        dumpB("ctintvegas6d      ", ctintvegas6d        );
        dumpB("ctintvegas8d      ", ctintvegas8d        );
        dumpI("intDimVJ          ", intDimVJ            );
        dumpB("vjint3d           ", vjint3d             );
        dumpB("vjint5d           ", vjint5d             );
        dumpB("vjintvegas7d      ", vjintvegas7d        );
        dumpI("intDimFO          ", intDimFO            );
        dumpB("foint0d           ", foint0d             );
        dumpB("foint3d           ", foint3d             );
        dumpI("intDimFPC          ", intDimFPC            );
        dumpB("fpcint2d           ", fpcint2d             );
        dumpB("fpcintvegas8d      ", fpcintvegas8d        );
        dumpB("fixedorder_only   ", fixedorder          );
        dumpB("primed   ", primed          );
        dumpB("BORNquad          ", BORNquad              );
        dumpB("CTquad            ", CTquad                );
        dumpB("VJquad            ", VJquad                );
        dumpB("doBORN            ", doBORN              );
        dumpB("doCT              ", doCT                );
        dumpB("doVJ              ", doVJ                );
        dumpB("doVJREAL          ", doVJREAL            );
        dumpB("doVJVIRT          ", doVJVIRT            );
        dumpB("doCTAS3           ", doCTAS3             );
        dumpB("doVJAS3           ", doVJAS3             );
        dumpB("doVJVV           ", doVJVV            );
        dumpB("doVJVR           ", doVJVR            );
        dumpB("doVJRR           ", doVJRR            );
        dumpB("doFPC             ", doFPC               );
        dumpB("finord            ", finord           );
        dumpI("threading         ", threading           );
        dumpI("cubaverbosity     ", cubaverbosity       );
        dumpI("cores             ", cubacores           );
        dumpI("cubanbatch        ", cubanbatch          );
        dumpI("niterBORN         ", niterBORN           );
        dumpI("niterCT           ", niterCT             );
        dumpI("niterVJ           ", niterVJ             );
        dumpD("vegasncallsBORN   ", vegasncallsBORN     );
        dumpD("vegasncallsCT     ", vegasncallsCT       );
        dumpD("vegasncallsVJLO   ", vegasncallsVJLO     );
        dumpD("vegasncallsVJREAL ", vegasncallsVJREAL   );
        dumpD("vegasncallsVJVIRT ", vegasncallsVJVIRT   );
        dumpD("vegasncallsVJVV   ", vegasncallsVJVV    );
        dumpD("vegasncallsVJVR   ", vegasncallsVJVR    );
        dumpD("vegasncallsVJRR   ", vegasncallsVJRR    );
	dumpB("vegascollect      ", vegascollect        );
	dumpB("vegascorr         ", vegascorr           );
        dumpB("pcubature         ", pcubature           );
        dumpD("relaccuracy       ", relaccuracy         );
        dumpD("absaccuracy       ", absaccuracy         );	
        dumpI("level             ", level               );
        dumpD("costhmin          ", costhmin            );
        dumpD("costhmax          ", costhmax            );
        dumpD("xfmin          ", xfmin            );
        dumpD("xfmax          ", xfmax            );
        dumpB("makecuts          ", makecuts            );
	dumpB("qed               ", qed                 );
	dumpB("qedorder          ", qedorder            );
        dumpD("lptcut            ", lptcut              );
        dumpD("lycut             ", lycut               );
        dumpD("mtcut             ", mtcut               );
	dumpD("cthCSmin          ", cthCSmin            );
	dumpD("cthCSmax          ", cthCSmax            );
        dumpD("etmisscut         ", etmisscut           );
        dumpB("cubaint           ", cubaint             );
        dumpB("trapezint         ", trapezint           );
        dumpB("quadint           ", quadint             );
        dumpI("suavepoints       ", suavepoints         );
        dumpI("nphitrape         ", nphitrape           );
        dumpI("phiintervals      ", phiintervals        );
        dumpI("phirule           ", phirule             );
        dumpI("ncstart           ", ncstart             );
        dumpB("qtrec_naive       ", qtrec_naive         );
        dumpB("qtrec_cs          ", qtrec_cs            );
        dumpB("qtrec_kt0         ", qtrec_kt0           );
        dumpB("timeprofile       ", timeprofile         );
        dumpB("verbose           ", verbose             );
        dumpB("gridverbose       ", gridverbose         );	
        dumpB("unicode           ", unicode         );	
        dumpB("silent            ", silent          );	
        dumpB("makehistos        ", makehistos         );	
        dumpB("texttable         ", texttable           );
        dumpB("resumcpp          ", resumcpp            );
	dumpB("ctcpp             ", ctcpp               );
        dumpI("approxpdf         ", opts_.approxpdf_    );
        dumpI("pdfintervals      ", opts_.pdfintervals_ );
        dumpI("evolmode          ", evolmode            );
        dumpB("iterative         ", iterative            );
        dumpB("hfmode            ", hfmode            );
        dumpB("smallx            ", smallx            );
        dumpI("nsub              ", nsub            );
        //dumpB("mufevol           ", mufevol            );
        dumpB("mufvar            ", mufvar            );
        dumpB("commute           ", commute            );
        dumpI("expc              ", expc             );
        dumpI("ntaylor           ", ntaylor             );
        dumpD("expcreg           ", expcreg             );
	dumpB("sumlogs           ", sumlogs       );
	dumpB("numsud            ", numsud       );
    dumpB("asrgkt            ", asrgkt       );
	dumpB("deltaismb         ", deltaismb       );
	dumpB("deltafsmb         ", deltafsmb       );
	dumpB("deltaismc         ", deltaismc       );
	dumpB("deltafsmc         ", deltafsmc       );
	dumpB("vfncqg            ", vfncqg       );
	dumpI("asrun             ", asrun       );
	dumpB("cmw               ", cmw       );
	dumpB("cmwfo             ", cmwfo     );
	dumpB("numexpc           ", numexpc      );
	dumpI("bprescription     ", bprescription       );
	dumpI("fbstar     ", fbstar       );
	dumpB("phibr             ", phibr       );
	dumpB("bcf               ", bcf       );
	//        dumpB("PDFerrors         ", PDFerrors           );
        dumpI("mellintr          ", mellintr     );
        dumpI("mellininv         ", mellininv     );
	dumpI("mellinintervals   ", mellinintervals     );
        dumpI("mellinrule        ", mellinrule          );
        dumpD("zmax              ", zmax                );
        dumpI("ncycle            ", ncycle                );
	dumpD("cpoint            ", cpoint              );
	dumpD("cshift            ", cshift              );
	dumpD("phi               ", phi              );
        dumpD("mellincores       ", mellincores         );
	dumpB("mellin1d          ", mellin1d            );
	dumpI("melup             ", melup            );
	dumpB("xspace            ", xspace              );
	dumpB("ctmellin          ", ctmellin              );
        dumpI("mrule             ", mrule               );
        dumpI("yintervals        ", yintervals          );
        dumpI("yrule             ", yrule               );
        dumpI("qtintervals       ", qtintervals          );
        dumpI("qtrule            ", qtrule               );
        dumpI("abintervals       ", abintervals          );
        dumpI("abrule            ", abrule               );
	dumpI("vjphirule         ", vjphirule            );
	dumpI("zrule             ", zrule               );
	dumpI("xrule             ", xrule               );
        dumpB("ptbinwidth        ", ptbinwidth          );
        dumpB("ybinwidth         ", ybinwidth           );
        dumpB("mbinwidth         ", mbinwidth           );
        dumpB("dsdqt2            ", dsdqt2           );
        dumpB("dsdxf             ", dsdxf           );
        dumpB("edsdp3            ", edsdp3           );
        dumpB("m3dsdm            ", m3dsdm           );
        dumpB("exactkin          ", exactkin           );
        dumpB("force_binsampling ", force_binsampling   );
        dumpB("lepbounds "        , lepbounds   );
        dumpI("helicity ",          helicity   );
        dumpB("h2d"               , h2d   );
        dumpB("h3d"               , h3d   );
        dumpB("hetal"             , hetal   );
        dumpB("hptl"              , hptl   );
        dumpB("hptnu"             , hptnu   );
        dumpB("hmt"               , hmt   );
        dumpB("haiqt"             , haiqt   );
        dumpB("haiy"              , haiy   );
        dumpB("haim"              , haim   );
        dumpS("output_filename ",   output_filename   );
    }

    if (print_masses){
        printf("Masses and EW constants:\n");
        dumpD ( "md        " , dymasses_   . md_        );
        dumpD ( "mu        " , dymasses_   . mu_        );
        dumpD ( "ms        " , dymasses_   . ms_        );
        dumpD ( "mc        " , dymasses_   . mc_        );
        dumpD ( "mb        " , dymasses_   . mb_        );
        dumpD ( "mt        " , dymasses_   . mt_        );
        dumpD ( "mel       " , dymasses_   . mel_       );
        dumpD ( "mmu       " , dymasses_   . mmu_       );
        dumpD ( "mtau      " , dymasses_   . mtau_      );
        dumpD ( "hmass     " , dymasses_   . hmass_     );
        dumpD ( "hwidth    " , dymasses_   . hwidth_    );
        dumpD ( "wmass     " , dymasses_   . wmass_     );
        dumpD ( "wwidth    " , dymasses_   . wwidth_    );
        dumpD ( "zmass     " , dymasses_   . zmass_     );
        dumpD ( "zwidth    " , dymasses_   . zwidth_    );
        dumpD ( "twidth    " , dymasses_   . twidth_    );
        dumpD ( "mtausq    " , dymasses_   . mtausq_    );
        dumpD ( "mcsq      " , dymasses_   . mcsq_      );
        dumpD ( "mbsq      " , dymasses_   . mbsq_      );
	/*
        dumpD ( "Gf_inp    " , ewinput_  . Gf_inp_    );
        dumpD ( "aemmz_inp " , ewinput_  . aemmz_inp_ );
        dumpD ( "xw_inp    " , ewinput_  . xw_inp_    );
        dumpD ( "wmass_inp " , ewinput_  . wmass_inp_ );
        dumpD ( "zmass_inp " , ewinput_  . zmass_inp_ );
	*/
        dumpD ( "Vud       " , cabib_    . Vud_       );
        dumpD ( "Vus       " , cabib_    . Vus_       );
        dumpD ( "Vub       " , cabib_    . Vub_       );
        dumpD ( "Vcd       " , cabib_    . Vcd_       );
        dumpD ( "Vcs       " , cabib_    . Vcs_       );
        dumpD ( "Vcb       " , cabib_    . Vcb_       );
        dumpD ( "epinv     " , mcfm_epinv_    . epinv_     );
        dumpD ( "epinv2    " , mcfm_epinv2_   . epinv2_    );
        dumpD ( "hmass     " , dymasses_   . hmass_     );
        dumpD ( "hwidth    " , dymasses_   . hwidth_    );
        dumpD ( "zmass     " , dymasses_   . zmass_     );
        dumpD("hwidth"             , hwidth);
    }

    printf("end of setting dump.\n\n");
}

void settings::dumpI( string var,int    val){
    printf( " %s = %d\n", var.c_str(), val);
}
void settings::dumpD( string var,double val){
    printf( " %s = %f\n", var.c_str(), val);
}
void settings::dumpS( string var,string val){
    printf( " %s = %s\n", var.c_str(), val.c_str());
}
void settings::dumpB( string var,bool val){
    printf( " %s = %s\n", var.c_str(), val ? "true" : "false" );
}

vector<string> settings::Tokenize(string val,char Delim){
    vector<string> vec;
    size_t pos = 0;
    string tmp;
    while (!val.empty()){
        // find delim
        pos = val.find_first_of(Delim);
        if (pos==string::npos){
            // last item
            tmp = val;
            val.clear();
        } else {
            // get substring
            tmp = val.substr(0,pos);
            val = val.substr(pos+1);
        }
        // add to vector
        if (!tmp.empty()) vec.push_back(tmp);
    }
    return vec;
}

bool settings::IsNumber(const string &s) {
    double dummy;
    try {
        dummy = stod(s);
        return true;
    } catch (const std::exception &e){
        return false;
    }
}

void settings::parse_binning(string name, vector<double> &bins, po::Options &args){
    if (args.count(name)) {
        string e("Unsupported value of "+name+" : need 3 numbers seperated by comma: 'N,lo,hi' ");
        string val=args[name.c_str()].as<string>();
        ToLower(val);
        vector<string> vec = Tokenize(val);
        // check value
        if (vec.size()!=3) throw QuitProgram(e+" not 3 numbers.");
        for (auto s : vec) if (!IsNumber(s))  throw QuitProgram(e+" not numbers.");
        // retrieve N,lo,hi
        int N = stod(vec[0]);
        double lo = stod(vec[1]);
        double hi = stod(vec[2]);
        if (lo>hi)  throw QuitProgram(e+" lo is more then hi.");
        if (N<1)  throw QuitProgram(e+" N is not at least 1.");
        // make binning
        bins.clear();
	//loop with double has problems for equality test (try --ybins 25,0,5)
	//for (double loedge=lo; loedge<=hi; loedge+=(hi-lo)/double(N)) bins.push_back(loedge);
	for (int i=0; i <= N; i++) bins.push_back(lo+i*(hi-lo)/double(N));
    }
}

void binning::GetBins(string name, vector<double> &vec){
    vec.clear();
    // CLI defined
    if (name == "qt" && qtbins .size()!=0) { vec = qtbins; return; }
    if (name == "y"  && ybins  .size()!=0) { vec = ybins;  return; }
    if (name == "m"  && mbins  .size()!=0) { vec = mbins;  return; }
    if (name == "ptl"  && ptlbins  .size()!=0) { vec = ptlbins;  return; }
    if (name == "etal" && etalbins .size()!=0) { vec = etalbins;  return; }
    // input file defined
    name+="_bins";
    in.GetVectorDouble(name.c_str(), vec);
    if ( vec.size() < 2) throw QuitProgram( "Option `" +name+ "` needs at least 2 items ");
}

void binning::readfromfile(const string fname){
    in.parse_file(fname);
    in.GetVectorDouble("qt_bins" , qtbins );
    in.GetVectorDouble("y_bins"  , ybins  );
    in.GetVectorDouble("m_bins"  , mbins  );
    //lepton bins are possible only with cuts on leptons, and when the charged lepton is uniquely defined, i.e. in W processes
    //if (opts.makecuts && opts.nproc != 3)
    //{
	in.GetVectorDouble("ptl_bins"  , ptlbins  );
	in.GetVectorDouble("etal_bins"  , etalbins  );
	//}
    hist_qt_bins .clear();
    hist_y_bins  .clear();
    hist_m_bins  .clear();
    return;
}

// InputParser definitions
//
InputParser::InputParser( string _filename, string _charset, string _white):
    filename ( _filename   ),
    Ccommnt  ( _charset[0] ),
    Cassign  ( _charset[1] ),
    CopenAr  ( _charset[2] ),
    CdeliAr  ( _charset[3] ),
    CclosAr  ( _charset[4] ),
    Swhite   ( _white      )
{
  // parse default
  string defname;
  char* dytpath = getenv ("DYTPATH");
  if (dytpath != NULL && file_exists(string(dytpath)+"/default.in"))
    defname = string(dytpath)+"/default.in";
  else if (file_exists((string)SHAREDIR+"/default.in"))
    defname = (string)SHAREDIR+"/default.in";
  else if (file_exists("share/dyturbo/default.in"))
    defname = "share/dyturbo/default.in";
  
  try {parse_file(defname);}
  catch (invalid_argument &e1)
    {
      throw invalid_argument( string( "Can not load default settings:\n")
			      + string(e1.what()) + string("\n")
			      );
    }

  // parse custom settings
  if (filename.compare("")!=0) parse_file(filename);
}


InputParser::~InputParser(){
}

double InputParser::GetNumber(string name){
    has_key(name);
    string val = data[name];
    try { return stod(val);
    } catch (const std::exception &e){
        printf("Cannot read option '%s' with value '%s' as number.\n",name.c_str(), val.c_str() );
        throw e;
    }
}

string InputParser::GetString(string name){
    has_key(name);
    string val = data[name];
    trim(val);
    return val;
}


bool InputParser::GetBool(string name){
    has_key(name);
    string val = data[name];
    // lower case
    std::transform(val.begin(), val.end(), val.begin(), ::tolower);
    bool result = (val.compare(0,4,"true") == 0);
    //printf(" InputParser::GetBool processing key `%s` with value `%s` as result `%d`\n", name.c_str(), val.c_str(), result);
    return result;
}


void InputParser::GetVectorDouble(string name, vector<double> &vec){
    has_key(name);
    vec.clear();
    string val = data[name];
    // has open/close array
    size_t strBegin = val.find(CopenAr,0);
    size_t strEnd   = val.find(CclosAr,0);
    // we need both open and close otherwise is something wrong
    if (strBegin==string::npos || strEnd==string::npos) throw invalid_argument("Missing open/close character.");
    // parse what is between them
    val = val.substr(strBegin+1,strEnd-strBegin-1);
    size_t pos = 0;
    string tmp;
    while (!val.empty()){
        // find delim
        pos = val.find(CdeliAr);
        if (pos==string::npos){
            // last item
            tmp = val;
            val.clear();
        } else {
            // get substring
            tmp = val.substr(0,pos);
            val = val.substr(pos+1);
        }
        // get substring
        trim(tmp);
        // expecting number otherwise exception
        if (!tmp.empty()) vec.push_back(stod(tmp));
    }
    return;
}


void InputParser::parse_file(const string fname){
    ifstream fstrm(fname.c_str());
    if ( ! fstrm.good() ) throw invalid_argument( string("Uknown file name ") + fname );
    string line;
    // line by line
    while(getline(fstrm,line)) {
        size_t pos = 0;
        // check for comments
        pos = line.find(Ccommnt,pos);
        if (pos != string::npos) line = line.substr(0,pos);
        // split string by delim
        pos = line.find(Cassign,0);
        if (pos == string::npos || line.empty()) continue;
        string key  = line.substr(0,pos);
        string val = line.substr(pos+1);
        // remove leading / trailing whitespace
        trim(key);
        trim(val);
        /// @todo: Multiline setting. If has CopenAr load until CclosAr
        // save to map
        data[key]=val;
    }
    fstrm.close();
    return;
}


void InputParser::trim(string & str){
    size_t strBegin = str.find_first_not_of(Swhite);
    if (strBegin == std::string::npos){
        str="";
        return ; // no content
    }
    size_t strEnd = str.find_last_not_of(Swhite);
    size_t strRange = strEnd - strBegin + 1;
    str = str.substr(strBegin, strRange);
    return;
}


void InputParser::has_key(const string key){
    if ( data.count(key) == 0 ){
        string msg = "Missing setting with name '";
        msg += key;
        msg += "'";
        throw invalid_argument(msg.c_str());
    }
    return;
}
