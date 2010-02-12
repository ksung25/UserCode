#ifndef MYTOOLS_HH
#define MYTOOLS_HH

#include <TGraphAsymmErrors.h>

#ifdef USE_ROOSTATSCMS
#include "PhysicsTools/RooStatsCms/interface/ClopperPearsonBinomialInterval.h"
#include "PhysicsTools/RooStatsCms/interface/FeldmanCousinsBinomialInterval.h"
#endif

namespace toolbox 
{
Double_t calcEff(Int_t pass, Int_t total, Double_t *errl=0, Double_t *errh=0, Int_t method=0);
Double_t calcEff(Double_t pass, Double_t total, Double_t *errl=0, Double_t *errh=0, Int_t method=0);

TGraphAsymmErrors* makeEffGraph(UInt_t n, Double_t *binEdges, Int_t *numer, Int_t *denom, Int_t method=0);
TGraphAsymmErrors* makeEffGraph(UInt_t n, Double_t *binEdges, Double_t *numer, Double_t *denom, Int_t method=0);

Double_t deltaR(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2);

Double_t deltaPhi(Double_t phi1, Double_t phi2);

Int_t roundToInt(Double_t x);
}

//------------------------------------------------------------------------------------------------------------------------
Double_t toolbox::calcEff(Int_t pass, Int_t total, Double_t *errl, Double_t *errh, Int_t method)
{
  // method: 0 -> Bayes Divide
  //         1 -> Feldman-Cousins 
  //         2 -> Clopper-Pearson
  
  Double_t r = (total>0) ? (Double_t)pass/(Double_t)total : 0;
  if(errl) *errl = 0;
  if(errh) *errh = 0;
    
  if(method==0) {
    const Double_t conf = 0.68269;
    TGraphAsymmErrors g;
    Double_t mode, low, high;
    g.Efficiency(pass,total,conf,mode,low,high);
    if(errl) *errl = mode - low;
    if(errh) *errh = high - mode;
  }
#ifdef USE_ROOSTATSCMS
  const Double_t alpha = (1-0.68269);
  if(method==1) {
    FeldmanCousinsBinomialInterval fc;    
    fc.init(alpha);
    fc.calculate(pass,total);
    if(errl) *errl = r - fc.lower();
    if(errh) *errh = fc.upper() - r;
  }

  if(method==2) {
    ClopperPearsonBinomialInterval cp;    
    cp.init(alpha);
    cp.calculate(pass,total);
    if(errl) *errl = r - cp.lower();
    if(errh) *errh = cp.upper() - r;
  }
#endif  
  return r;
}

Double_t toolbox::calcEff(Double_t pass, Double_t total, Double_t *errl, Double_t *errh, Int_t method) 
{
  // Round values to whole numbers first
  return calcEff(roundToInt(pass),roundToInt(total),errl,errh,method);
}

//------------------------------------------------------------------------------------------------------------------------
TGraphAsymmErrors* toolbox::makeEffGraph(UInt_t n, Double_t *binEdges, Int_t *numer, Int_t *denom, Int_t method)
{
  // method: 0 -> Bayes Divide
  //         1 -> Feldman-Cousins 
  //         2 -> Clopper-Pearson
  
  if(method==0) {
    TGraphAsymmErrors *g = new TGraphAsymmErrors();
    TH1F hNumer("hNumer","",n,binEdges); 
    for(UInt_t ibin=0; ibin<n; ibin++) { hNumer.SetBinContent(ibin+1,numer[ibin]); }
    TH1F hDenom("hDenom","",n,binEdges); 
    for(UInt_t ibin=0; ibin<n; ibin++) { hDenom.SetBinContent(ibin+1,denom[ibin]); }
    g->BayesDivide(&hNumer,&hDenom);
    return g;
  }
#ifdef USE_ROOSTATSCMS
cout << "Using RooStatsCMS" << endl;
  Double_t xval[n], yval[n], xerr[n], yerrl[n], yerrh[n];
  const Double_t alpha = (1-0.68269);
  
  if(method==1) {
    FeldmanCousinsBinomialInterval fc;    
    fc.init(alpha);
    for(UInt_t ibin=0; ibin<n; ibin++) {
      fc.calculate(numer[ibin], denom[ibin]);
      xval[ibin]  = 0.5*(binEdges[ibin+1]+binEdges[ibin]);
      xerr[ibin]  = 0.5*(binEdges[ibin+1]-binEdges[ibin]);
      yval[ibin]  = (denom[ibin]>0) ? (Double_t)numer[ibin]/(Double_t)denom[ibin] : 0;
      yerrl[ibin] = yval[ibin] - fc.lower();
      yerrh[ibin] = fc.upper() - yval[ibin];
    }
    return new TGraphAsymmErrors(n,xval,yval,xerr,xerr,yerrl,yerrh); 
  }
  
  if(method==2) {
    ClopperPearsonBinomialInterval cp;    
    cp.init(alpha);
    for(UInt_t ibin=0; ibin<n; ibin++) {
      cp.calculate(numer[ibin], denom[ibin]);
      xval[ibin]  = 0.5*(binEdges[ibin+1]+binEdges[ibin]);
      xerr[ibin]  = 0.5*(binEdges[ibin+1]-binEdges[ibin]);
      yval[ibin]  = (denom[ibin]>0) ? (Double_t)numer[ibin]/(Double_t)denom[ibin] : 0;
      yerrl[ibin] = yval[ibin] - cp.lower();
      yerrh[ibin] = cp.upper() - yval[ibin];
    }
    return new TGraphAsymmErrors(n,xval,yval,xerr,xerr,yerrl,yerrh); 
  }    
#endif
  return 0;
}

TGraphAsymmErrors* toolbox::makeEffGraph(UInt_t n, Double_t *binEdges, Double_t *numer, Double_t *denom, Int_t method)
{
  Int_t v1[n], v2[n];
  for(UInt_t i=0; i<n; i++) {
    v1[i] = roundToInt(numer[i]);
    v2[i] = roundToInt(denom[i]);
  }
  
  return makeEffGraph(n,binEdges,v1,v2,method);
}

//------------------------------------------------------------------------------------------------------------------------
Double_t toolbox::deltaR(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2)
{
  const Double_t pi = 3.14159265358979323846;
  Double_t dphi = fabs(phi1-phi2);
  while (dphi>pi)
    dphi = fabs(dphi - 2.0*pi);
    
  Double_t deta = eta1-eta2;
  
  return sqrt(dphi*dphi + deta*deta);
}

//------------------------------------------------------------------------------------------------------------------------
Double_t toolbox::deltaPhi(Double_t phi1, Double_t phi2) 
{
  // Compute dPhi between two given angles. Results is in [-pi/2,pi/2].
  const Double_t pi = 3.14159265358979323846;
  Double_t dphi = fabs(phi1-phi2);
  while (dphi>pi)
    dphi = fabs(dphi - 2.0*pi);
  return(dphi);
}

//------------------------------------------------------------------------------------------------------------------------
Int_t toolbox::roundToInt(Double_t x)
{
  if(x>0)
    return ((x-floor(x)) < (ceil(x)-x)) ? (Int_t)floor(x) : (Int_t)ceil(x);
  else
    return ((x-floor(x)) < (ceil(x)-x)) ? (Int_t)ceil(x) : (Int_t)floor(x);
}

#endif
