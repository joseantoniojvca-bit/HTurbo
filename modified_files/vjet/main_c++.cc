//       main_c++.cc
//
//      Copyright 2011 Leandro Cieri
//
//      This program is free software; you can redistribute it and/or modify
//      it under the terms of the GNU General Public License as published by
//      the Free Software Foundation; either version 2 of the License, or
//      (at your option) any later version.
//
//      This program is distributed in the hope that it will be useful,
//      but WITHOUT ANY WARRANTY; without even the implied warranty of
//      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//      GNU General Public License for more details.
//
//      You should have received a copy of the GNU General Public License
//      along with this program; if not, write to the Free Software
//      Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
//      MA 02110-1301, USA.
using namespace std;

#include <fstream>
#include <iostream>
#include <sstream>
#include <cstring>
#include <stdio.h>
#include "Set_Root_workspace.h"
//c_mstwpdf *mstpdf;

void BookRootObjects();
void SaveRootObjects();
#include "params.h"

extern "C" {

	void bookmbobjects_();

}

extern "C" {
    void bookplot_(int* hident,const char* nameHist, int* nbins , double* xmin,
        double* xmax,const char* var,const char* unit,const char* yscale,
        int* nunmber,int len_a,int len_b,int len_c,int len_d);
}

inline void bookplotc(int hident,const char *nameHist, int nbins , double xmin,
    double xmax, char* var,char* unit,char* yscale,int nunmber) {
    bookplot_(&hident,nameHist, &nbins , &xmin,
    &xmax, var, unit, yscale,&nunmber,strlen(nameHist),strlen(var),strlen(unit),strlen(yscale));
}



/* Wrapper a Fortran */
extern "C" int vegasmain_();
int main()
{

  
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;
 std::cout <<"CCCCCCCCCCCCCCC   2gNNLO, Version 2.1   CCCCCCCCCCCCC"<<endl;
 std::cout <<"C                                                   C"<<endl;
 std::cout <<"C  Written by S. Catani, D. de Florian, L. Cieri    C"<<endl;
 std::cout <<"C             G. Ferrera and M. Grazzini            C"<<endl;
 std::cout <<"C                                                   C"<<endl;
 std::cout <<"C  Please refer to:                                 C"<<endl;
 std::cout <<"C       S. Catani, D. de Florian, L. Cieri,         C"<<endl;
 std::cout <<"C           G. Ferrera and M. Grazzini,             C"<<endl;
 std::cout <<"C        Phys. Rev.Lett. 108 (2012) 072001          C"<<endl;
 std::cout <<"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"<<endl;






    // read and dump to screen user configuration options
    ReadParams("input.txt");
    DumpParams();

  // Puntero a PDF's
  // char filename[100];
  // char prefix[] = "Grids/mstw2008nnlo"; // prefix for the grid files

  // Consider only the central PDF set to start with.
  //  sprintf(filename,"%s.%2.2d.dat",prefix,0);
  //  mstpdf = new c_mstwpdf(filename); // default: warn=false, fatal=true




    Integral_class Int_evnt;

    _Evento evnt;

    evnt.Set_VEGAS(evnt);

    charge_VEGAS(evnt);

    BookRootObjects();

    vegasmain_();

    SaveRootObjects();

    return 0;
}

string Unit(string var) {
    if (var=="Mgg" || var=="pTgg" || var=="pT") return "GeV";
    else if (var=="dphi") return "rad";
    else return "";
}

string Label(string var) {
    if (var=="Mgg") return "m_{#gamma#gamma}";
    else if (var=="eta") return "#eta_{#gamma#gamma}";
    else if (var=="pT") return "p_{T}"; // what's this??
    else if (var=="pTgg") return "p_{T,#gamma#gamma}";
    else if (var=="dphi") return "#Delta#phi";
    else if (var=="CosP") return "cos#theta*";
    else if (var=="CosA") return "cos#theta*";
    else return "";
}

void BookRootObjects() {
    if (!SaveRootNtuple && !SaveRootHistos) return;

    // open file for output
    F = new TFile("ReRoot.root","RECREATE");

    nhists=0;
    T=0;
    // process configuration file
    ifstream infile("input.txt");
    if (!infile.good()) {
	exit(-1);
    }
    while (!infile.eof()) {
	char buf[500];
	infile.getline(buf,500,'\n');
	if (buf[0]=='#') continue;
	stringstream str(buf);
	string token;
	str >> token;
	if( !infile.good()  ) break;
	if (token=="hist") {
	    // add new histogram
	    int ihist, nbins;
	    string var;
	    double xmin, xmax;
	    string yscale;
	    str >> ihist >> var >> nbins >> xmin >> xmax >> yscale;

	    char name[50];
	    sprintf(name,"hist%d",ihist);
	    //int idhist = ihist-1;
	    int idhist = nhists;
	    hists[idhist] = new TH1D(name,name,nbins,xmin,xmax);
	    string unit = Unit(var);
	    string label = Label(var);
	    if (unit!="") {
		sprintf(name,"%s [%s]",label.c_str(),unit.c_str());
		hists[idhist]->GetXaxis()->SetTitle(name);
		sprintf(name,"d#sigma/d%s [fb/%s]",label.c_str(),unit.c_str());
		hists[idhist]->GetYaxis()->SetTitle(name);
	    }
	    else {
		sprintf(name,"%s",label.c_str());
		hists[idhist]->GetXaxis()->SetTitle(name);
		sprintf(name,"d#sigma/d%s [fb]",label.c_str());
		hists[idhist]->GetYaxis()->SetTitle(name);
	    }
	    histvar[idhist] = var;
	    nhists++;
	}
    }

    // create output ntuple, if desired
    if (SaveRootNtuple)
	T = new TNtuple("leantuple","datacoment","Mff:weight");
}

void SaveRootObjects() {
    if (!SaveRootNtuple && !SaveRootHistos) return;

    // save ntuple
    if (SaveRootNtuple)
	T-> Write();

    // save histograms
    if (SaveRootHistos) {
	for (int i=0; i<nhists; i++)
	    hists[i]->Write();
    }

    // close output file
    F->Close();
}

extern "C"
{  // PROCESO TODO PROCESO CINEMATICA BORN
  void integrandob_(double &winteb, double weght) {

  contrib = 1.0;
  contpt = 1U;

  winteb = 0.0;

  Integral_class Inte;

  Inte.vweght=weght;

  winteb=Inte.Integrando_Fborn(weght);


 return;
   }
  }



  extern "C"
{  // PROCESO TODO PROCESO CINEMATICA BORN
  void integrandoboxh_( double &winteb,double weght) {

  contrib = 1.0;
  contpt = 1U;

  winteb = 0.0;

  Integral_class Inte;

  Inte.vweght=weght;

  winteb=Inte.Integrando_boxh(weght);


 return;
   }
  }

extern "C"
{  // PROCESO TODO PROCESO CINEMATICA BORN
  void integrandoctnlo_( double &wwct,double weght) {

  contrib = 2.0;
  contpt = 2U;


  wwct = 0.0;

  Integral_class Inte;

  Inte.vweght=weght;

  wwct=Inte.Integrando_Full_nlo_CT(weght);
  //wwct=Inte.Integrando_CTgg(weght);


 return;
   }
  }


  extern "C"
{  // PROCESO TODO PROCESO CINEMATICA BORN
  void integrandoctnnlo_( double &wwctnn,double weght) {

  contrib = 2.0;
  contpt = 2U;

  wwctnn = 0.0;

  Integral_class Inte;

  Inte.vweght=weght;

  wwctnn=Inte.Integrando_Full_nnlo_CT(weght);


 return;
   }
  }


   extern "C"
{  // PROCESO TODO PROCESO CINEMATICA BORN
  void integrandoh2_( double &wwh2,double weght) {

  contrib = 2.0;
  contpt = 1U;

  wwh2 = 0.0;

  Integral_class Inte;

  Inte.vweght=weght;

  wwh2=Inte.Integrando_h2(weght);


 return;
   }
  }


extern "C"
{  // PROCESO  CINEMATICA Full NNLO: 2Real + Dip + fini y Real NLO
  void integrandofullnlo_( double &wfinten,double weght) {

  contrib = 2.0;
  contpt = 1U;

  wfinten = 0.0;



  Integral_class Inte;

  Inte.vweght=weght;

  wfinten=Inte.Integrando_Full_nlo(weght);




 return;
   }
  }

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------


extern "C"
{  // PROCESO  CINEMATICA Full NLO: N^3LO ggffg
  void integgffg_( double &wfinten,double weght) {

  contrib = 2.0;
  contpt = 1U;

  wfinten = 0.0;



  Integral_class Inte;

  Inte.vweght=weght;

  wfinten=Inte.Integrando_ggffg(weght);




 return;
   }
  }
  
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
 
  
extern "C"
{  // PROCESO  CINEMATICA Full NLO: N^3LO ggHg
  void integghg_( double &wfinten,double weght) {

  contrib = 2.0;
  contpt = 1U;

  wfinten = 0.0;

 // std::cout << " desde integghg_ " << std::endl;

  Integral_class Inte;

  Inte.vweght=weght;

  wfinten=Inte.Integrando_gghg(weght);




 return;
   }
  }
  
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
  
extern "C"
{  // PROCESO  CINEMATICA Full NLO: N^3LO ggINTg
  void integgintg_( double &wfinten,double weght) {

  contrib = 2.0;
  contpt = 1U;

  wfinten = 0.0;



  Integral_class Inte;

  Inte.vweght=weght;

  wfinten=Inte.Integrando_ggintg(weght);




 return;
   }
  }

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// Interferencia LO --- LEAN
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
  
extern "C"
{  // PROCESO  CINEMATICA Full NLO: N^3LO ggINTg
  void integghyy_( double &wgghyy,double weght) {

  contrib = 1.0;
  contpt = 1U;

  wgghyy = 0.0;

  Integral_class Inte;

  Inte.vweght=weght;

  wgghyy=Inte.Integrando_gghyy(weght);




 return;
   }
  }
  
extern "C"
{  // PROCESO  CINEMATICA Full NLO: N^3LO ggINTg QG AND QQBAR
  void inteqqbarhyy_( double &wqqbarhyy,double weght) {

  contrib = 2.0;
  contpt = 1U;

  wqqbarhyy = 0.0;

  Integral_class Inte;

  Inte.vweght=weght;

  wqqbarhyy=Inte.Integrando_qqbarhyy(weght);


 return;
   }
  }
    
  
  // CT gg NLO
  
extern "C"
{  // PROCESO TODO PROCESO CINEMATICA BORN
  void integrandoctggnlo_( double &wwct,double weght) {

  contrib = 2.0;
  contpt = 2U;


  wwct = 0.0;

  Integral_class Inte;

  Inte.vweght=weght;

  //wwct=Inte.Integrando_Full_nlo_CT(weght);
  wwct=Inte.Integrando_CTgg(weght);


 return;
   }
  }

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------

extern "C"
{  // PROCESO TODO PROCESE CINEMATICA NNLO 2Real + Dip + Fini
  void integrandon_( double &winten,double weght) {

  contrib = 3.0;
  contpt = 1U;

  winten = 0.0;

  Integral_class Inte;

  Inte.vweght=weght;

  winten=Inte.Integrando_Full_nnlo(weght);


 return;
   }
  }




 extern "C"
{  // PROCESO TODO PROCESE CINEMATICA NNLO 2Real + Dip + Fini
  void integrandofini_( double &winten,double weght) {

  contrib = 2.0;
  contpt = 1U;

  winten = 0.0;

  Integral_class Inte;

  Inte.vweght=weght;

  winten=Inte.Integrando_fini(weght);


 return;
   }
  }


void  bookmbobjects_() {
    if (!SaveTopdHistos) return;



    nhists=0;
    // process configuration file
    ifstream infile("input.txt");
    if (!infile.good()) {
	exit(-1);
    }
    while (!infile.eof()) {
	char buf[500];
	infile.getline(buf,500,'\n');
	if (buf[0]=='#') continue;
	stringstream str(buf);
	string token;
	str >> token;
	if( !infile.good()  ) break;
	if (token=="hist") {
	    // add new histogram
	    int ihist, nbins;
	    string var;
	    double xmin, xmax;
	    string yscale;
	    str >> ihist >> var >> nbins >> xmin >> xmax >> yscale;

        char *escala=new char[yscale.size()+1];
        escala[yscale.size()]=0;
        memcpy(escala,yscale.c_str(),yscale.size());

	    char *vari=new char[var.size()+1];
        vari[var.size()]=0;
        memcpy(vari,var.c_str(),var.size());



	    char name[50];
	    sprintf(name,"hist%d",ihist);
	    //int idhist = ihist-1;
	    int idhist = nhists;
//	    hists[idhist] = new TH1D(name,name,nbins,xmin,xmax);
	    bookplotc(idhist+1,name, nbins ,  xmin,
        xmax, vari,"[fb]", escala,1);
	    histvar[idhist] = var;
	    nhists++;
	}
    }


}
