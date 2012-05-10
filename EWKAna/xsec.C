#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <sstream>                  // class for parsing strings
#endif

void xsec(const TString infilename="input.txt")
{
  Double_t nWmp, nWmm, nWm, nZmm, nWep, nWem, nWe, nZee;
  Double_t accWmp, accWmm, accWm, accZmm, accWep, accWem, accWe, accZee;
  Double_t lumi, lumiErr;
  
  Double_t errWmp=0, errWmm=0, errWm=0, errZmm=0, errWmpm=0, errWZm=0;
  Double_t errWep=0, errWem=0, errWe=0, errZee=0, errWepm=0, errWZe=0; 
  
  ifstream ifs;
  ifs.open(infilename.Data());
  assert(ifs.is_open());
  string line;
  Int_t state=0;
  while(getline(ifs,line)) {
    if(line[0]=='#') continue;
    
    if(state==0) {
      stringstream ss(line);
      ss >> lumi >> lumiErr;
      state++;
      
    } else if(state==1) {
      stringstream ss(line);
      string label;
      ss >> label >> nWmp >> nWmm >> nWm >> nZmm >> nWep >> nWem >> nWe >> nZee;
      state++;
      
    } else if(state==2) {
      stringstream ss(line);
      string label;
      ss >> label >> accWmp >> accWmm >> accWm >> accZmm >> accWep >> accWem >> accWe >> accZee;
      state++;
    
    } else if(state==3) {  
      stringstream ss(line);
      string label;
      state++;
    
    } else {
      stringstream ss(line);
      string label;
      Double_t err[12];
      ss >> label >> err[0] >> err[1] >> err[2] >> err[3] >> err[4] >> err[5] >> err[6] >> err[7] >> err[8] >> err[9] >> err[10] >> err[11];
      errWmp  += err[0] *err[0];
      errWmm  += err[1] *err[1];
      errWm   += err[2] *err[2];
      errZmm  += err[3] *err[3];
      errWmpm += err[4] *err[4];
      errWZm  += err[5] *err[5];
      errWep  += err[6] *err[6];
      errWem  += err[7] *err[7];
      errWe   += err[8] *err[8];
      errZee  += err[9] *err[9];
      errWepm += err[10]*err[10];
      errWZe  += err[11]*err[11];
    }
  }
  ifs.close();
  
  cout << endl;   
  cout << setprecision(2) << fixed; 
  
  Double_t xsWmp = 1.005*nWmp/accWmp/lumi/1000.;
  cout << "    W+(mu): " << setw(8) << xsWmp  << " +/- " << setw(7) << xsWmp/sqrt(nWmp)             << " (stat) +/- " << setw(7) << sqrt(errWmp)*xsWmp   << " (syst) +/- " << setw(7) << xsWmp*lumiErr << " (lumi) pb" << endl;
  
  Double_t xsWmm = 1.005*nWmm/accWmm/lumi/1000.;
  cout << "    W-(mu): " << setw(8) << xsWmm  << " +/- " << setw(7) << xsWmm/sqrt(nWmm)             << " (stat) +/- " << setw(7) << sqrt(errWmm)*xsWmm   << " (syst) +/- " << setw(7) << xsWmm*lumiErr << " (lumi) pb" << endl;
  
  Double_t xsWm = 1.005*nWm/accWm/lumi/1000.;
  cout << "     W(mu): " << setw(8) << xsWm   << " +/- " << setw(7) << xsWm/sqrt(nWm)               << " (stat) +/- " << setw(7) << sqrt(errWm)*xsWm     << " (syst) +/- " << setw(7) <<  xsWm*lumiErr  << " (lumi) pb" <<  endl;
  
  Double_t xsZmm  = 1.01*nZmm/accZmm/lumi/1000.;
  cout << "     Z(mu): " << setw(8) << xsZmm  << " +/- " << setw(7) << xsZmm/sqrt(nZmm)             << " (stat) +/- " << setw(7) << sqrt(errZmm)*xsZmm   << " (syst) +/- " << setw(7) << xsZmm*lumiErr << " (lumi) pb" << endl;
  
  cout << setprecision(4) << fixed;
  
  Double_t xsWmpm = nWmp/nWmm/accWmp*accWmm;
  cout << " W+/W-(mu): " << setw(8) << xsWmpm << " +/- " << setw(7) << sqrt(1./nWmp+1./nWmm)*xsWmpm << " (stat) +/- " << setw(7) << sqrt(errWmpm)*xsWmpm << " (syst)" << endl;
  
  Double_t xsWZm = nWm/nZmm/accWm*accZmm;
  cout << "   W/Z(mu): " << setw(8) << xsWZm  << " +/- " << setw(7) << sqrt(1./nWm+1./nZmm)*xsWZm   << " (stat) +/- " << setw(7) << sqrt(errWZm)*xsWZm   << " (syst)" << endl;

  cout << endl;
  cout << setprecision(2) << fixed; 
  
  Double_t xsWep = nWep/accWep/lumi/1000.;
  cout << "     W+(e): " << setw(8) << xsWep  << " +/- " << setw(7) << xsWep/sqrt(nWep)             << " (stat) +/- " << setw(7) << sqrt(errWep)*xsWep   << " (syst) +/- " << setw(7) << xsWep*lumiErr << " (lumi) pb" << endl;
  
  Double_t xsWem = nWem/accWem/lumi/1000.;
  cout << "     W-(e): " << setw(8) << xsWem  << " +/- " << setw(7) << xsWem/sqrt(nWem)             << " (stat) +/- " << setw(7) << sqrt(errWem)*xsWem   << " (syst) +/- " << setw(7) << xsWem*lumiErr << " (lumi) pb" << endl;
  
  Double_t xsWe = nWe/accWe/lumi/1000.;
  cout << "      W(e): " << setw(8) << xsWe   << " +/- " << setw(7) << xsWe/sqrt(nWe)               << " (stat) +/- " << setw(7) << sqrt(errWe)*xsWe     << " (syst) +/- " << setw(7) << xsWe*lumiErr  << " (lumi) pb" <<  endl;
  
  Double_t xsZee = nZee/accZee/lumi/1000.;
  cout << "      Z(e): " << setw(8) << xsZee  << " +/- " << setw(7) << xsZee/sqrt(nZee)             << " (stat) +/- " << setw(7) << sqrt(errZee)*xsZee   << " (syst) +/- " << setw(7) << xsZee*lumiErr << " (lumi) pb" << endl;
  
  cout << setprecision(4) << fixed;
  
  Double_t xsWepm = nWep/nWem/accWep*accWem;
  cout << "  W+/W-(e): " << setw(8) << xsWepm << " +/- " << setw(7) << sqrt(1./nWep+1./nWem)*xsWepm << " (stat) +/- " << setw(7) << sqrt(errWepm)*xsWepm << " (syst)" << endl;
  
  Double_t xsWZe = nWe/nZee/accWe*accZee;
  cout << "    W/Z(e): " << setw(8) << xsWZe  << " +/- " << setw(7) << sqrt(1./nWe+1./nZee)*xsWZe   << " (stat) +/- " << setw(7) << sqrt(errWZe)*xsWZe   << " (syst)" << endl;

  cout << endl;
}
