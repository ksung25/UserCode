//================================================================================================
//
// Plot
//
//  * 
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                        // access to gROOT, entry point to ROOT system
#include <TSystem.h>                      // interface to OS
#include <TFile.h>                        // file handle class
#include <TTree.h>                        // class to access ntuples
#include <TBenchmark.h>                   // class to track macro running statistics
#include <TH1D.h>                         // histogram class
#include <vector>                         // STL vector class
#include <iostream>                       // standard I/O
#include <iomanip>                        // functions to format standard I/O
#include <fstream>                        // functions for file I/O
#include <string>                         // C++ string class
#include <sstream>                        // class for parsing strings
#include "Math/LorentzVector.h"           // 4-vector class

#include "../Utils/MyTools.hh"	          // various helper functions
#include "../Utils/CPlot.hh"	          // helper class for plots
#include "../Utils/MitStyleRemix.hh"      // style settings for drawing

#include "../Utils/ZSignals.hh"           // define models for Z signal PDFs
#include "../Utils/ZBackgrounds.hh"       // define models for background PDFs

// RooFit headers
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooCategory.h"
#include "RooArgList.h"
#include "RooDataHist.h"
#include "RooFormulaVar.h"
#include "RooHistPdf.h"
#include "RooGenericPdf.h"
#include "RooAddPdf.h"
#include "RooSimultaneous.h"
#include "RooGaussian.h"
#include "RooNLLVar.h"
#include "RooConstVar.h"
#include "RooMinuit.h"
#include "RooFitResult.h"
#endif

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVector;

enum { eCount, eBWxCB, eMCxGaus };
enum { eNone, eExp, eErfcExp, eDblExp, eLinExp, eQuadExp };


//=== FUNCTION DECLARATIONS ======================================================================================

// make webpage
void makeHTML(const TString outDir);

// obtain a signal+background PDF
RooAbsPdf* getModel(Int_t sigType, Int_t bkgType, const char *label, TH1D *histTemplate,
                    RooRealVar &m, RooFormulaVar &NfitSig, RooRealVar &NfitBkg);


//=== MAIN MACRO ================================================================================================= 

void fitZmm(const TString infilename="/data/blue/ksung/EWKAna/test/Selection/Zmumu/ntuples/data_select.root",      // input ntuple
            const TString outputDir="test"        // output directory
) {
  gBenchmark->Start("fitZmm");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================   
  
  const Bool_t doBinned = kTRUE;
  
  const Double_t MASS_LOW  = 60;
  const Double_t MASS_HIGH = 120;
  const Int_t    NBINS     = 30;
  
  const Int_t typeMuMuNoSelSig = eBWxCB;
  const Int_t typeMuMuNoSelBkg = eExp;
  const Int_t typeMuStaSig     = eBWxCB;
  const Int_t typeMuStaBkg     = eExp;
  const Int_t typeMuTrkSig     = eBWxCB;
  const Int_t typeMuTrkBkg     = eExp;
  
  const TString format("png");

  
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  

  enum { eMuMu2HLT=1, eMuMu1HLT, eMuMuNoSel, eMuSta, eMuTrk };
  
  // Create output directory
  gSystem->mkdir(outputDir,kTRUE);
  CPlot::sOutDir = outputDir;  
  
  Float_t mass;
  
  TTree *treeMuMu2HLT = new TTree("treeMuMu2HLT","treeMuMu2HLT");
  treeMuMu2HLT->Branch("m",&mass,"m/F");
  treeMuMu2HLT->SetDirectory(0);
  TH1D *hMuMu2HLT = new TH1D("hMuMu2HLT","",NBINS,MASS_LOW,MASS_HIGH);
  hMuMu2HLT->Sumw2();
  
  TTree *treeMuMu1HLT = new TTree("treeMuMu1HLT","treeMuMu1HLT");
  treeMuMu1HLT->Branch("m",&mass,"m/F");
  treeMuMu1HLT->SetDirectory(0);
  TH1D *hMuMu1HLT = new TH1D("hMuMu1HLT","",NBINS,MASS_LOW,MASS_HIGH);
  hMuMu1HLT->Sumw2();
  
  TTree *treeMuMuNoSel = new TTree("treeMuMuNoSel","treeMuMuNoSel");
  treeMuMuNoSel->Branch("m",&mass,"m/F");
  treeMuMuNoSel->SetDirectory(0);
  TH1D *hMuMuNoSel = new TH1D("hMuMuNoSel","",NBINS,MASS_LOW,MASS_HIGH);
  hMuMuNoSel->Sumw2();
  
  TTree *treeMuSta = new TTree("treeMuSta","treeMuSta");
  treeMuSta->Branch("m",&mass,"m/F");
  treeMuSta->SetDirectory(0);
  TH1D *hMuSta = new TH1D("hMuSta","",NBINS,MASS_LOW,MASS_HIGH);
  hMuSta->Sumw2();
  
  TTree *treeMuTrk = new TTree("treeMuTrk","treeMuTrk");
  treeMuTrk->Branch("m",&mass,"m/F");
  treeMuTrk->SetDirectory(0);
  TH1D *hMuTrk = new TH1D("hMuTrk","",NBINS,MASS_LOW,MASS_HIGH); 
  hMuTrk->Sumw2();
  
  //
  // Declare variables to read in ntuple
  //
  UInt_t  runNum, lumiSec, evtNum;
  UInt_t  category;
  UInt_t  npv, npu;
  Float_t scale1fb;
  Float_t met, metPhi, sumEt, u1, u2;
  Int_t   q1, q2;
  LorentzVector *dilep=0, *lep1=0, *lep2=0;

  TFile *infile=0;
  TTree *intree=0;

  // Read input file and get the TTrees
  cout << "Processing " << infilename << "..." << endl;
  infile = new TFile(infilename);	  assert(infile);
  intree = (TTree*)infile->Get("Events"); assert(intree);

  intree->SetBranchAddress("runNum",   &runNum);     // event run number
  intree->SetBranchAddress("lumiSec",  &lumiSec);    // event lumi section
  intree->SetBranchAddress("evtNum",   &evtNum);     // event number
  intree->SetBranchAddress("category", &category);   // dilepton category
  intree->SetBranchAddress("npv",      &npv);	     // number of primary vertices
  intree->SetBranchAddress("npu",      &npu);	     // number of in-time PU events (MC)
  intree->SetBranchAddress("scale1fb", &scale1fb);   // event weight per 1/fb (MC)
  intree->SetBranchAddress("met",      &met);	     // MET
  intree->SetBranchAddress("metPhi",   &metPhi);     // phi(MET)
  intree->SetBranchAddress("sumEt",    &sumEt);      // Sum ET
  intree->SetBranchAddress("u1",       &u1);	     // parallel component of recoil
  intree->SetBranchAddress("u2",       &u2);	     // perpendicular component of recoil
  intree->SetBranchAddress("q1",       &q1);	     // charge of tag lepton
  intree->SetBranchAddress("q2",       &q2);	     // charge of probe lepton
  intree->SetBranchAddress("dilep",    &dilep);      // dilepton 4-vector
  intree->SetBranchAddress("lep1",     &lep1);       // tag lepton 4-vector
  intree->SetBranchAddress("lep2",     &lep2);       // probe lepton 4-vector
  
  //
  // loop over events
  //
  for(UInt_t ientry=0; ientry<intree->GetEntries(); ientry++) {
    intree->GetEntry(ientry);
   
    mass = dilep->M();
    
    if     (category == eMuMu2HLT)  { treeMuMu2HLT->Fill();  hMuMu2HLT->Fill(mass); }
    else if(category == eMuMu1HLT)  { treeMuMu1HLT->Fill();  hMuMu1HLT->Fill(mass); }
    else if(category == eMuMuNoSel) { treeMuMuNoSel->Fill(); hMuMuNoSel->Fill(mass); }
    else if(category == eMuSta)     { treeMuSta->Fill();     hMuSta->Fill(mass); }
    else if(category == eMuTrk)     { treeMuTrk->Fill();     hMuTrk->Fill(mass); }
  }
  
  UInt_t nMuMu2HLT  = treeMuMu2HLT->GetEntries();
  UInt_t nMuMu1HLT  = treeMuMu1HLT->GetEntries();
  UInt_t nMuMuNoSel = treeMuMuNoSel->GetEntries();
  UInt_t nMuSta     = treeMuSta->GetEntries();
  UInt_t nMuTrk     = treeMuTrk->GetEntries();
      
  RooRealVar m("m","m",MASS_LOW,MASS_HIGH);
  m.setBins(NBINS);

  RooAbsData *dataMuMuNoSel=0;
  RooAbsData *dataMuSta=0;
  RooAbsData *dataMuTrk=0;
  RooAbsData *dataNonGolden=0;

  
  // Combine low purity categories into one dataset for simultaneous fitting
  // The golden categories 2HLT and 1HLT are not fit: they enter minimization
  // via the constraint terms.
  RooCategory rooCat("rooCat","rooCat");
  rooCat.defineType("MuMuNoSel");
  rooCat.defineType("MuSta");
  rooCat.defineType("MuTrk");
    
  if(doBinned) {
    dataMuMuNoSel = new RooDataHist("dataMuMuNoSel","dataMuMuNoSel",RooArgSet(m),hMuMuNoSel);
    dataMuSta	  = new RooDataHist("dataMuSta",    "dataMuSta",    RooArgSet(m),hMuSta);
    dataMuTrk	  = new RooDataHist("dataMuTrk",    "dataMuTrk",    RooArgSet(m),hMuTrk);
    dataNonGolden = new RooDataHist("dataNonGolden","dataNonGolden", RooArgList(m), Index(rooCat),
                                    Import("MuMuNoSel", *((RooDataHist*)dataMuMuNoSel)),
				    Import("MuSta",     *((RooDataHist*)dataMuSta)),
                                    Import("MuTrk",     *((RooDataHist*)dataMuTrk)));
  } else {  
    dataMuMuNoSel = new RooDataSet("dataMuMuNoSel","dataMuMuNoSel",treeMuMuNoSel,RooArgSet(m));
    dataMuSta	  = new RooDataSet("dataMuSta",	   "dataMuSta",	   treeMuSta,    RooArgSet(m));
    dataMuTrk	  = new RooDataSet("dataMuTrk",	   "dataMuTrk",	   treeMuTrk,    RooArgSet(m));
    dataNonGolden = new RooDataSet("dataNonGolden","dataNonGolden", RooArgList(m), Index(rooCat),
                                   Import("MuMuNoSel", *((RooDataSet*)dataMuMuNoSel)),
				   Import("MuSta",     *((RooDataSet*)dataMuSta)),
                                   Import("MuTrk",     *((RooDataSet*)dataMuTrk)));
  }
  
  //
  // Create PDFs
  //
  
  // Primary parameters
  UInt_t NzMax = 2*(nMuMu2HLT + nMuMu1HLT + nMuMuNoSel + nMuSta + nMuTrk);
  RooRealVar Nz("Nz","Nz",nMuMu2HLT+nMuMu1HLT,0,NzMax);
  RooRealVar effHLT("effHLT","effHLT",0.90,0.80,1.0);
  RooRealVar effSel("effSel","effSel",0.95,0.80,1.0);
  RooRealVar effTrk("effTrk","effTrk",0.98,0.90,1.0);
  RooRealVar effSta("effSta","effSta",0.97,0.90,1.0);    
  
  // The expected background yields in the non-golden samples for extended likelihood
  RooRealVar NfitBkgMuMuNoSel("NfitBkgMuMuNoSel","NfitBkgMuMuNoSel",0.2*nMuMuNoSel,0,nMuMuNoSel);
  RooRealVar NfitBkgMuSta    ("NfitBkgMuSta",    "NfitBkgMuSta",    0.2*nMuSta,    0,nMuSta);
  RooRealVar NfitBkgMuTrk    ("NfitBkgMuTrk",    "NfitBkgMuTrk",    0.2*nMuTrk,    0,nMuTrk);  
  
  // The expected numbers of signal events in each sample
  RooFormulaVar NfitMuMu2HLT("NfitMuMu2HLT","NfitMuMu2HLT",
                             "Nz*effHLT*effHLT*effTrk*effTrk*effSta*effSta*effSel*effSel",
                             RooArgList(Nz,effHLT,effTrk,effSta,effSel));
  RooFormulaVar NfitMuMu1HLT("NfitMuMu1HLT","NfitMuMu1HLT",
                             "2*Nz*effHLT*(1-effHLT)*effTrk*effTrk*effSta*effSta*effSel*effSel",
                             RooArgList(Nz,effHLT,effTrk,effSta,effSel));
  RooFormulaVar NfitMuMuNoSel("NfitMuMuNoSel","NfitMuMuNoSel",
                              "2*Nz*effHLT*effTrk*effTrk*effSta*effSta*effSel*(1-effSel)",
                              RooArgList(Nz,effHLT,effTrk,effSta,effSel)); 
  RooFormulaVar NfitMuSta("NfitMuSta","NfitMuSta",
                          "2*Nz*effHLT*effTrk*(1-effTrk)*effSta*effSta*effSel",
                          RooArgList(Nz,effHLT,effTrk,effSta,effSel));
  RooFormulaVar NfitMuTrk("NfitMuTrk","NfitMuTrk",
                          "2*Nz*effHLT*effTrk*effTrk*effSta*(1-effSta)*effSel",
                          RooArgList(Nz,effHLT,effTrk,effSta,effSel));
  
  TH1D *tmpltMuMuNoSel=0;
  TH1D *tmpltMuSta=0;
  TH1D *tmpltMuTrk=0;
  
  //
  // Put together total PDFs
  //
  RooAbsPdf *pdfMuMuNoSel = getModel(typeMuMuNoSelSig, typeMuMuNoSelBkg, "MuMuNoSel", tmpltMuMuNoSel, m, NfitMuMuNoSel, NfitBkgMuMuNoSel);
  RooAbsPdf *pdfMuSta     = getModel(typeMuStaSig,     typeMuStaBkg,     "MuSta",     tmpltMuSta,     m, NfitMuSta,     NfitBkgMuSta);
  RooAbsPdf *pdfMuTrk     = getModel(typeMuTrkSig,     typeMuTrkBkg,     "MuTrk",     tmpltMuTrk,     m, NfitMuTrk,     NfitBkgMuTrk); 

  // PDF for simultaneous fit
  RooSimultaneous pdfTotal("pdfTotal","pdfTotal",rooCat);
  pdfTotal.addPdf(*pdfMuMuNoSel,"MuMuNoSel");
  pdfTotal.addPdf(*pdfMuSta,    "MuSta");  
  pdfTotal.addPdf(*pdfMuTrk,    "MuTrk");			     
                             
  //
  // Define likelihood, add constraints, and run the fit
  //
  
  // Extra terms to likelihood
  RooGaussian constraintMuMu2HLT("constraintMuMu2HLT","constraintMuMu2HLT", NfitMuMu2HLT, RooConst(nMuMu2HLT), RooConst(sqrt(nMuMu2HLT)));
  RooGaussian constraintMuMu1HLT("constraintMuMu1HLT","constraintMuMu1HLT", NfitMuMu1HLT, RooConst(nMuMu1HLT), RooConst(sqrt(nMuMu1HLT)));
  RooGaussian constraintMuMuNoSel("constraintMuMuNoSel","constraintMuMuNoSel", NfitMuMuNoSel, RooConst(nMuMuNoSel), RooConst(sqrt(nMuMuNoSel)));
  RooGaussian constraintMuSta("constraintMuSta","constraintMuSta", NfitMuSta, RooConst(nMuSta), RooConst(sqrt(nMuSta)));
  RooGaussian constraintMuTrk("constraintMuTrk","constraintMuTrk", NfitMuTrk, RooConst(nMuTrk), RooConst(sqrt(nMuTrk)));
  
//  RooPoisson constraintMuMu2HLT("constraintMuMu2HLT","constraintMuMu2HLT", NfitMuMu2HLT, RooConst(nMuMu2HLT));
//  RooPoisson constraintMuMu1HLT("constraintMuMu1HLT","constraintMuMu1HLT", NfitMuMu1HLT, RooConst(nMuMu1HLT));
//  RooPoisson constraintMuMuNoSel("constraintMuMuNoSel","constraintMuMuNoSel", NfitMuMuNoSel, RooConst(nMuMuNoSel));
//  RooPoisson constraintMuSta("constraintMuSta","constraintMuSta", NfitMuSta, RooConst(nMuSta));
//  RooPoisson constraintMuTrk("constraintMuTrk","constraintMuTrk", NfitMuTrk, RooConst(nMuTrk));

  // Define goodness of fit including the constraints
  RooArgList fitConstraints;
  fitConstraints.add(constraintMuMu2HLT);
  fitConstraints.add(constraintMuMu1HLT);
  if(typeMuMuNoSelSig==eCount) fitConstraints.add(constraintMuMuNoSel);
  if(typeMuStaSig==eCount)     fitConstraints.add(constraintMuSta);
  if(typeMuTrkSig==eCount)     fitConstraints.add(constraintMuTrk);
/*  
  RooNLLVar *nll = (RooNLLVar*)pdfTotal.createNLL(*dataNonGolden, Extended(kTRUE), ExternalConstraints(fitConstraints));

  RooMinuit rminuit(*nll);
  cout << "Explicitly telling Minuit to use error level 0.5 for likelihood" << endl;
  rminuit.setErrorLevel(0.5);
  RooFitResult *result = rminuit.fit("smhr");
*/
  RooFitResult *result = pdfTotal.fitTo(*dataNonGolden, Extended(kTRUE), ExternalConstraints(fitConstraints), Save(kTRUE));


  //--------------------------------------------------------------------------------------------------------------
  // Make plots 
  //==============================================================================================================  
  TCanvas *c = MakeCanvas("c","c",800,600);
  
  char ylabel[100];
  char cattext[100];
  
  sprintf(ylabel,"Events / %.1f GeV/c^{2}",hMuMu2HLT->GetBinWidth(1));
  sprintf(cattext,"#mu+#mu (2HLT)");
  CPlot plotMuMu2HLT("mumu2hlt","","mass [GeV/c^{2}]",ylabel);
  plotMuMu2HLT.AddHist1D(hMuMu2HLT,"E");
  plotMuMu2HLT.AddTextBox(cattext,0.70,0.85,0.90,0.79,0);
  plotMuMu2HLT.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events / %.1f GeV/c^{2}",hMuMu1HLT->GetBinWidth(1));
  sprintf(cattext,"#mu+#mu (1HLT)");
  CPlot plotMuMu1HLT("mumu1hlt","","mass [GeV/c^{2}]",ylabel);
  plotMuMu1HLT.AddHist1D(hMuMu1HLT,"E");
  plotMuMu1HLT.AddTextBox(cattext,0.70,0.85,0.90,0.79,0);
  plotMuMu1HLT.Draw(c,kTRUE,format);
  
  sprintf(ylabel,"Events / %.1f GeV/c^{2}",hMuMuNoSel->GetBinWidth(1));
  sprintf(cattext,"#mu+#mu (No Sel)");
  CPlot plotMuMuNoSel("mumunosel","","mass [GeV/c^{2}]",ylabel);
  plotMuMuNoSel.AddHist1D(hMuMuNoSel,"E");
  plotMuMuNoSel.AddTextBox(cattext,0.70,0.85,0.90,0.79,0);
  plotMuMuNoSel.Draw(c,kTRUE,format);
  
  sprintf(ylabel,"Events / %.1f GeV/c^{2}",hMuSta->GetBinWidth(1));
  sprintf(cattext,"#mu+#mu_{STA}");
  CPlot plotMuSta("musta","","mass [GeV/c^{2}]",ylabel);
  plotMuSta.AddHist1D(hMuSta,"E");
  plotMuSta.AddTextBox(cattext,0.70,0.85,0.90,0.79,0);
  plotMuSta.Draw(c,kTRUE,format);
  
  sprintf(ylabel,"Events / %.1f GeV/c^{2}",hMuTrk->GetBinWidth(1));
  sprintf(cattext,"#mu+track");
  CPlot plotMuTrk("mutrk","","mass [GeV/c^{2}]",ylabel);
  plotMuTrk.AddHist1D(hMuTrk,"E");
  plotMuTrk.AddTextBox(cattext,0.70,0.85,0.90,0.79,0);
  plotMuTrk.Draw(c,kTRUE,format);  
  
    
  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
   
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;  
  cout << endl;
  cout << "   MuMu2HLT: " << setw(8) << nMuMu2HLT << " events" << endl; 
  cout << "   MuMu1HLT: " << setw(8) << nMuMu1HLT << " events" << endl; 
  cout << "  MuMuNoSel: " << setw(8) << nMuMuNoSel << " events" << endl;
  cout << "      MuSta: " << setw(8) << nMuSta << " events" << endl;    
  cout << "      MuTrk: " << setw(8) << nMuTrk << " events" << endl;      
  cout << endl;
  cout << "       N_Z: " << Nz.getVal()     << " +/- " << Nz.getPropagatedError(*result) << endl;
  cout << "  eff(HLT): " << effHLT.getVal() << " +/- " << effHLT.getPropagatedError(*result) << endl;
  cout << "  eff(Sel): " << effSel.getVal() << " +/- " << effSel.getPropagatedError(*result) << endl;
  cout << "  eff(Trk): " << effTrk.getVal() << " +/- " << effTrk.getPropagatedError(*result) << endl;
  cout << "  eff(Sta): " << effSta.getVal() << " +/- " << effSta.getPropagatedError(*result) << endl;
  cout << endl;
  
  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;     
  
  gBenchmark->Show("fitZmm");
}


void makeHTML(const TString outDir)
{
}


RooAbsPdf* getModel(Int_t sigType, Int_t bkgType, const char *label, TH1D *histTemplate,
                    RooRealVar &m, RooFormulaVar &NfitSig, RooRealVar &NfitBkg)
{
  char name[100];
  
  sprintf(name,"sig%s",label);
  CSignalModel *sigModel=0;  
  if     (sigType==eCount)   { sigModel = new CBreitWignerConvCrystalBall(name, m); } 
  else if(sigType==eBWxCB)   { sigModel = new CBreitWignerConvCrystalBall(name, m); } 
  else if(sigType==eMCxGaus) { sigModel = new CMCTemplateConvGaussian(name, m, histTemplate); } 
  else { 
    cout << "Not valid background model choice for " << label << endl;
    assert(0);
  }
  
  sprintf(name,"bkg%s",label);  
  CBackgroundModel *bkgModel=0;
  if     (bkgType==eNone)    { NfitBkg.setVal(0); } 
  else if(bkgType==eExp)     { bkgModel = new CExponential(name, m); } 
  else if(bkgType==eErfcExp) { bkgModel = new CErfExpo(name, m); } 
  else if(bkgType==eDblExp)  { bkgModel = new CDoubleExp(name, m); } 
  else if(bkgType==eLinExp)  { bkgModel = new CLinearExp(name, m); } 
  else if(bkgType==eQuadExp) { bkgModel = new CQuadraticExp(name, m); } 
  else {    
    cout << "Not valid signal model choice for MuSta!" << endl;
    assert(0);
  }
  
  sprintf(name,"pdf%s",label);
  return new RooAddPdf(name,name,
                       (bkgType==eNone) ? RooArgList(*(sigModel->model)) : RooArgList(*(sigModel->model),*(bkgModel->model)),
		       (bkgType==eNone) ? RooArgList(NfitSig) : RooArgList(NfitSig,NfitBkg));  
}
