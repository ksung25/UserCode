//================================================================================================
//
// Study isolation efficiency on Z sample and QCD sample
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <TFile.h>
#include <TRFIOFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TVirtualFitter.h>
#include <TBenchmark.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <vector>
#include <iostream>
#include <iomanip>
#endif

#include "CPlot.hh"
#include "MitStyleRemix.hh"

#include "ZAnaStructDefs.hh"     // define structures to read in ntuple

//#define USE_ROOSTATSCMS
#include "MyTools.hh"

void plotIsoEffMC() {
  
  gBenchmark->Start("plotIsoEffMC");
  gSystem->Load("libRFIO");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================
  
  Bool_t doSave  = false;   // save plots?
  CPlot::sOutDir = ".";    // output directory
  TString format = "png";  // output file format
  Int_t method   = 0;      // efficiency errors calculation method
  
  // signal sample
  //---------------
  TString fnameSig("/castor/cern.ch/user/k/ksung/ntuples/00/s09-zmm-7-mc3_ntuple.root");
//  TString fnameSig("/castor/cern.ch/user/k/ksung/ntuples/00/s09-zmm-mc3_ntuple.root");
  TString labelSig("Z #rightarrow #mu#mu");  
  
  // QCD background sample
  //-----------------------  
  TString fnameBkg("/castor/cern.ch/user/k/ksung/ntuples/00/s09-incmu15-7-mc3_ntuple.root");
//  TString fnameBkg("/castor/cern.ch/user/k/ksung/ntuples/00/s09-incmu15-mc3_ntuple.root");
  TString labelBkg("QCD"); 

  //
  // Scan settings
  //
  Double_t iso03Low = 2.0, iso03High = 12.0, iso03Step = 1.0;
  Double_t iso05Low = 3.0, iso05High = 21.0, iso05Step = 2.0;
  Double_t relIso03Low = 0.03, relIso03High = 0.31, relIso03Step = 0.02;
  Double_t relIso05Low = 0.10, relIso05High = 0.40, relIso05Step = 0.02;
  
  //
  // Cuts and requirements
  //
  const Double_t muPtCut        = 20;
  const Double_t muEtaCut       = 2;
  const Double_t massMin        = 40;
  const Double_t trkIso03Cut    = 3;
  const Double_t trkIso05Cut    = 6;
  const Double_t trkRelIso03Cut = 0.09;
  const Double_t trkRelIso05Cut = 0.18;
  
  //
  // Canvas dimensions
  //
  Int_t canw=800, canh=600;

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================
  
  //
  // Variables for efficiency calculation in isolation scan
  //  
  Int_t denomSig = 0, denomBkg = 0;
  
  const Int_t nIso03Steps = Int_t((iso03High-iso03Low)/iso03Step)+1;
  const Int_t nIso05Steps = Int_t((iso05High-iso05Low)/iso05Step)+1;
  Int_t numerIso03Sig[nIso03Steps];
  Int_t numerIso03Bkg[nIso03Steps];
  for(Int_t i=0; i<nIso03Steps; i++) {
    numerIso03Sig[i] = 0;
    numerIso03Bkg[i] = 0;
  }
  Int_t numerIso05Sig[nIso05Steps];
  Int_t numerIso05Bkg[nIso05Steps];
  for(Int_t i=0; i<nIso05Steps; i++) {
    numerIso05Sig[i] = 0;
    numerIso05Bkg[i] = 0;
  }
  Double_t iso03EffSig[nIso03Steps]; Double_t iso03ErrlSig[nIso03Steps]; Double_t iso03ErrhSig[nIso03Steps];
  Double_t iso03EffBkg[nIso03Steps]; Double_t iso03ErrlBkg[nIso03Steps]; Double_t iso03ErrhBkg[nIso03Steps];
  Double_t iso05EffSig[nIso05Steps]; Double_t iso05ErrlSig[nIso05Steps]; Double_t iso05ErrhSig[nIso05Steps];
  Double_t iso05EffBkg[nIso05Steps]; Double_t iso05ErrlBkg[nIso05Steps]; Double_t iso05ErrhBkg[nIso05Steps];
    
  const Int_t nRelIso03Steps = Int_t(relIso03High/relIso03Step-relIso03Low/relIso03Step)+1;
  const Int_t nRelIso05Steps = Int_t(relIso05High/relIso05Step-relIso05Low/relIso05Step)+1;
  Int_t numerRelIso03Sig[nRelIso03Steps];
  Int_t numerRelIso03Bkg[nRelIso03Steps];
  for(Int_t i=0; i<nRelIso03Steps; i++) {
    numerRelIso03Sig[i] = 0;
    numerRelIso03Bkg[i] = 0;
  }
  Int_t *numerRelIso05Sig = new Int_t[nRelIso05Steps];
  Int_t *numerRelIso05Bkg = new Int_t[nRelIso05Steps];
  for(Int_t i=0; i<nRelIso05Steps; i++) {
    numerRelIso05Sig[i] = 0;
    numerRelIso05Bkg[i] = 0;
  }
  Double_t relIso03EffSig[nRelIso03Steps]; Double_t relIso03ErrlSig[nRelIso03Steps]; Double_t relIso03ErrhSig[nRelIso03Steps];
  Double_t relIso03EffBkg[nRelIso03Steps]; Double_t relIso03ErrlBkg[nRelIso03Steps]; Double_t relIso03ErrhBkg[nRelIso03Steps];
  Double_t relIso05EffSig[nRelIso05Steps]; Double_t relIso05ErrlSig[nRelIso05Steps]; Double_t relIso05ErrhSig[nRelIso05Steps];
  Double_t relIso05EffBkg[nRelIso05Steps]; Double_t relIso05ErrlBkg[nRelIso05Steps]; Double_t relIso05ErrhBkg[nRelIso05Steps];
    
  Int_t numerIso03CutSig=0, numerIso03CutBkg=0;
  Int_t numerIso05CutSig=0, numerIso05CutBkg=0;
  Int_t numerRelIso03CutSig=0, numerRelIso03CutBkg=0;
  Int_t numerRelIso05CutSig=0, numerRelIso05CutBkg=0;  
  Double_t iso03CutEffSig, iso05CutEffSig, relIso03CutEffSig, relIso05CutEffSig; 
  Double_t iso03CutEffBkg, iso05CutEffBkg, relIso03CutEffBkg, relIso05CutEffBkg;
  
  //
  // Variables for efficiency calculation vs. pT, eta, phi
  //  
  Double_t ptBinEdges[] = { 20, 24, 26, 28, 30, 32, 34, 36, 38, 40, 
                            44, 48, 52, 60, 70, 90 };
  const UInt_t nPtBins = sizeof(ptBinEdges)/sizeof(Double_t) - 1;			    		    			    
  Int_t denomPtSig[nPtBins];
  Int_t denomPtBkg[nPtBins];
  Int_t numerIso03vPtSig[nPtBins];
  Int_t numerIso03vPtBkg[nPtBins]; 
  Int_t numerIso05vPtSig[nPtBins];
  Int_t numerIso05vPtBkg[nPtBins];
  Int_t numerRelIso03vPtSig[nPtBins];
  Int_t numerRelIso03vPtBkg[nPtBins]; 
  Int_t numerRelIso05vPtSig[nPtBins];
  Int_t numerRelIso05vPtBkg[nPtBins];
  for(UInt_t i=0; i<nPtBins; i++) {
    denomPtSig[i] = 0;
    denomPtBkg[i] = 0;
    numerIso03vPtSig[i] = 0;
    numerIso03vPtBkg[i] = 0;
    numerIso05vPtSig[i] = 0;
    numerIso05vPtBkg[i] = 0;
    numerRelIso03vPtSig[i] = 0;
    numerRelIso03vPtBkg[i] = 0;
    numerRelIso05vPtSig[i] = 0;
    numerRelIso05vPtBkg[i] = 0;
  }
  
  Double_t etaBinEdges[21];
  const UInt_t nEtaBins = sizeof(etaBinEdges)/sizeof(Double_t) - 1; 
  for(UInt_t i=0; i<nEtaBins+1; i++) { etaBinEdges[i] = -2.0+i*0.2; }
  Int_t denomEtaSig[nEtaBins];
  Int_t denomEtaBkg[nEtaBins];
  Int_t numerIso03vEtaSig[nEtaBins];
  Int_t numerIso03vEtaBkg[nEtaBins]; 
  Int_t numerIso05vEtaSig[nEtaBins];
  Int_t numerIso05vEtaBkg[nEtaBins];
  Int_t numerRelIso03vEtaSig[nEtaBins];
  Int_t numerRelIso03vEtaBkg[nEtaBins]; 
  Int_t numerRelIso05vEtaSig[nEtaBins];
  Int_t numerRelIso05vEtaBkg[nEtaBins];
  for(UInt_t i=0; i<nEtaBins; i++) {
    denomEtaSig[i] = 0;
    denomEtaBkg[i] = 0;
    numerIso03vEtaSig[i] = 0;
    numerIso03vEtaBkg[i] = 0;
    numerIso05vEtaSig[i] = 0;
    numerIso05vEtaBkg[i] = 0;
    numerRelIso03vEtaSig[i] = 0;
    numerRelIso03vEtaBkg[i] = 0;
    numerRelIso05vEtaSig[i] = 0;
    numerRelIso05vEtaBkg[i] = 0;
  }

  Double_t phiBinEdges[33];
  const UInt_t nPhiBins = sizeof(phiBinEdges)/sizeof(Double_t) - 1; 
  for(UInt_t i=0; i<nPhiBins+1; i++) { phiBinEdges[i] = -3.2+i*0.2; }
  Int_t denomPhiSig[nPhiBins];
  Int_t denomPhiBkg[nPhiBins];
  Int_t numerIso03vPhiSig[nPhiBins];
  Int_t numerIso03vPhiBkg[nPhiBins]; 
  Int_t numerIso05vPhiSig[nPhiBins];
  Int_t numerIso05vPhiBkg[nPhiBins];
  Int_t numerRelIso03vPhiSig[nPhiBins];
  Int_t numerRelIso03vPhiBkg[nPhiBins]; 
  Int_t numerRelIso05vPhiSig[nPhiBins];
  Int_t numerRelIso05vPhiBkg[nPhiBins];
  for(UInt_t i=0; i<nPhiBins; i++) {
    denomPhiSig[i] = 0;
    denomPhiBkg[i] = 0;
    numerIso03vPhiSig[i] = 0;
    numerIso03vPhiBkg[i] = 0;
    numerIso05vPhiSig[i] = 0;
    numerIso05vPhiBkg[i] = 0;
    numerRelIso03vPhiSig[i] = 0;
    numerRelIso03vPhiBkg[i] = 0;
    numerRelIso05vPhiSig[i] = 0;
    numerRelIso05vPhiBkg[i] = 0;
  }
     
  //
  // Perform scan and calculation
  //
  TFile *infile=0;
  TTree *infoTree=0, *genTree=0, *dimuonTree=0, *muonTree=0;
    
  // Data structures to store info from TTrees
  TEventInfo eventInfo;
  TGenInfo   genInfo;
  TDimuon    dimuon;
  TMuon      muon;
  
  // signal sample
  if(fnameSig.Contains("/castor/cern.ch",TString::kIgnoreCase))
    infile = new TRFIOFile(fnameSig);
  else
    infile = new TFile(fnameSig); 
  assert(infile);
    
  infoTree   = (TTree*)infile->FindObjectAny("EventInfo"); assert(infoTree);
  genTree    = (TTree*)infile->FindObjectAny("GenInfo");   assert(genTree);
  dimuonTree = (TTree*)infile->FindObjectAny("Dimuon");    assert(dimuonTree);
  muonTree   = (TTree*)infile->FindObjectAny("Muon");      assert(muonTree);

  // Set branch address to structures that will store the info  
  infoTree  ->SetBranchAddress("EventInfo",&eventInfo);
  genTree   ->SetBranchAddress("GenInfo",  &genInfo);
  dimuonTree->SetBranchAddress("Dimuon",   &dimuon);
  muonTree  ->SetBranchAddress("Muon",     &muon);
 
  // loop over events
  for(UInt_t ientry=0; ientry<dimuonTree->GetEntries(); ientry++) {
    dimuonTree->GetEntry(ientry);
      
    assert(dimuon.eventId <= infoTree->GetEntries());
    infoTree->GetEntry(dimuon.eventId);
            
    //if(eventInfo.HLT_Mu15==0)                                          continue;  // trigger accept                                      
    if(dimuon.mass < massMin)                                          continue;  // mass region
    if(dimuon.pt_1 < muPtCut || dimuon.pt_2 < muPtCut)                 continue;  // muon pT cut
    if(fabs(dimuon.eta_1) > muEtaCut || fabs(dimuon.eta_2) > muEtaCut) continue;  // muon eta cut 
    if(dimuon.q_1 == dimuon.q_2)                                       continue;  // require opposite charge 
    
    // average efficiency for specific cut values
    if(dimuon.trkIso03_1 < trkIso03Cut) numerIso03CutSig++;
    if(dimuon.trkIso03_2 < trkIso03Cut) numerIso03CutSig++;
    if(dimuon.trkIso03_1/dimuon.pt_1 < trkRelIso03Cut) numerRelIso03CutSig++;
    if(dimuon.trkIso03_2/dimuon.pt_2 < trkRelIso03Cut) numerRelIso03CutSig++;
    
    if(dimuon.trkIso05_1 < trkIso05Cut) numerIso05CutSig++;
    if(dimuon.trkIso05_2 < trkIso05Cut) numerIso05CutSig++;
    if(dimuon.trkIso05_1/dimuon.pt_1 < trkRelIso05Cut) numerRelIso05CutSig++;
    if(dimuon.trkIso05_2/dimuon.pt_2 < trkRelIso05Cut) numerRelIso05CutSig++;
    
    // scan
    denomSig+=2;            
    for(Int_t istep=0; istep<nIso03Steps; istep++) {
      if(dimuon.trkIso03_1 < iso03Low+istep*iso03Step) numerIso03Sig[istep]++;
      if(dimuon.trkIso03_2 < iso03Low+istep*iso03Step) numerIso03Sig[istep]++;
    } 
    
    for(Int_t istep=0; istep<nIso05Steps; istep++) {
      if(dimuon.trkIso05_1 < iso05Low+istep*iso05Step) numerIso05Sig[istep]++;
      if(dimuon.trkIso05_2 < iso05Low+istep*iso05Step) numerIso05Sig[istep]++;
    }
    for(Int_t istep=0; istep<nRelIso03Steps; istep++) {
      if(dimuon.trkIso03_1/dimuon.pt_1 < relIso03Low+istep*relIso03Step) numerRelIso03Sig[istep]++;
      if(dimuon.trkIso03_2/dimuon.pt_2 < relIso03Low+istep*relIso03Step) numerRelIso03Sig[istep]++;
    }
    for(Int_t istep=0; istep<nRelIso05Steps; istep++) {
      if(dimuon.trkIso05_1/dimuon.pt_1 < relIso05Low+istep*relIso05Step) numerRelIso05Sig[istep]++;
      if(dimuon.trkIso05_2/dimuon.pt_2 < relIso05Low+istep*relIso05Step) numerRelIso05Sig[istep]++;
    }
    
        
    // pT dependence
    for(UInt_t ibin=0; ibin<nPtBins; ibin++) {
      if(dimuon.pt_1 > ptBinEdges[ibin] && dimuon.pt_1 < ptBinEdges[ibin+1]) {
        denomPtSig[ibin]++;
	if(dimuon.trkIso03_1 < trkIso03Cut) numerIso03vPtSig[ibin]++;
	if(dimuon.trkIso05_1 < trkIso05Cut) numerIso05vPtSig[ibin]++;
	if(dimuon.trkIso03_1/dimuon.pt_1 < trkRelIso03Cut) numerRelIso03vPtSig[ibin]++;
	if(dimuon.trkIso05_1/dimuon.pt_1 < trkRelIso05Cut) numerRelIso05vPtSig[ibin]++;
      }
      if(dimuon.pt_2 > ptBinEdges[ibin] && dimuon.pt_2 < ptBinEdges[ibin+1]) {
        denomPtSig[ibin]++;
	if(dimuon.trkIso03_2 < trkIso03Cut) numerIso03vPtSig[ibin]++;
	if(dimuon.trkIso05_2 < trkIso05Cut) numerIso05vPtSig[ibin]++;
	if(dimuon.trkIso03_2/dimuon.pt_2 < trkRelIso03Cut) numerRelIso03vPtSig[ibin]++;
	if(dimuon.trkIso05_2/dimuon.pt_2 < trkRelIso05Cut) numerRelIso05vPtSig[ibin]++;
      }
    }
    
    // eta dependence
    for(UInt_t ibin=0; ibin<nEtaBins; ibin++) {
      if(dimuon.eta_1 > etaBinEdges[ibin] && dimuon.eta_1 < etaBinEdges[ibin+1]) {
        denomEtaSig[ibin]++;
	if(dimuon.trkIso03_1 < trkIso03Cut) numerIso03vEtaSig[ibin]++;
	if(dimuon.trkIso05_1 < trkIso05Cut) numerIso05vEtaSig[ibin]++;
	if(dimuon.trkIso03_1/dimuon.pt_1 < trkRelIso03Cut) numerRelIso03vEtaSig[ibin]++;
	if(dimuon.trkIso05_1/dimuon.pt_1 < trkRelIso05Cut) numerRelIso05vEtaSig[ibin]++;
      }
      if(dimuon.eta_2 > etaBinEdges[ibin] && dimuon.eta_2 < etaBinEdges[ibin+1]) {
        denomEtaSig[ibin]++;
	if(dimuon.trkIso03_2 < trkIso03Cut) numerIso03vEtaSig[ibin]++;
	if(dimuon.trkIso05_2 < trkIso05Cut) numerIso05vEtaSig[ibin]++;
	if(dimuon.trkIso03_2/dimuon.pt_2 < trkRelIso03Cut) numerRelIso03vEtaSig[ibin]++;
	if(dimuon.trkIso05_2/dimuon.pt_2 < trkRelIso05Cut) numerRelIso05vEtaSig[ibin]++;
      }
    }
    
    // phi dependence
    for(UInt_t ibin=0; ibin<nPhiBins; ibin++) {
      if(dimuon.phi_1 > phiBinEdges[ibin] && dimuon.phi_1 < phiBinEdges[ibin+1]) {
        denomPhiSig[ibin]++;
	if(dimuon.trkIso03_1 < trkIso03Cut) numerIso03vPhiSig[ibin]++;
	if(dimuon.trkIso05_1 < trkIso05Cut) numerIso05vPhiSig[ibin]++;
	if(dimuon.trkIso03_1/dimuon.pt_1 < trkRelIso03Cut) numerRelIso03vPhiSig[ibin]++;
	if(dimuon.trkIso05_1/dimuon.pt_1 < trkRelIso05Cut) numerRelIso05vPhiSig[ibin]++;
      }
      if(dimuon.phi_2 > phiBinEdges[ibin] && dimuon.phi_2 < phiBinEdges[ibin+1]) {
        denomPhiSig[ibin]++;
	if(dimuon.trkIso03_2 < trkIso03Cut) numerIso03vPhiSig[ibin]++;
	if(dimuon.trkIso05_2 < trkIso05Cut) numerIso05vPhiSig[ibin]++;
	if(dimuon.trkIso03_2/dimuon.pt_2 < trkRelIso03Cut) numerRelIso03vPhiSig[ibin]++;
	if(dimuon.trkIso05_2/dimuon.pt_2 < trkRelIso05Cut) numerRelIso05vPhiSig[ibin]++;
      }
    }
  }
  delete infile;
  infile=0, infoTree=0, genTree=0, dimuonTree=0, muonTree=0; 

  // background sample
  if(fnameBkg.Contains("/castor/cern.ch",TString::kIgnoreCase))
    infile = new TRFIOFile(fnameBkg);
  else
    infile = new TFile(fnameBkg); 
  assert(infile);
    
  infoTree   = (TTree*)infile->FindObjectAny("EventInfo"); assert(infoTree);
  genTree    = (TTree*)infile->FindObjectAny("GenInfo");   assert(genTree);
  dimuonTree = (TTree*)infile->FindObjectAny("Dimuon");    assert(dimuonTree);
  muonTree   = (TTree*)infile->FindObjectAny("Muon");      assert(muonTree);

  // Set branch address to structures that will store the info  
  infoTree  ->SetBranchAddress("EventInfo",&eventInfo);
  genTree   ->SetBranchAddress("GenInfo",  &genInfo);
  dimuonTree->SetBranchAddress("Dimuon",   &dimuon);
  muonTree  ->SetBranchAddress("Muon",     &muon);
 
  // loop over events
  for(UInt_t ientry=0; ientry<dimuonTree->GetEntries(); ientry++) {
    dimuonTree->GetEntry(ientry);
      
    assert(dimuon.eventId <= infoTree->GetEntries());
    infoTree->GetEntry(dimuon.eventId);
            
    //if(eventInfo.HLT_Mu15==0)                                          continue;  // trigger accept                                      
    if(dimuon.mass < massMin)                                          continue;  // mass region
    if(dimuon.pt_1 < muPtCut || dimuon.pt_2 < muPtCut)                 continue;  // muon pT cut
    if(fabs(dimuon.eta_1) > muEtaCut || fabs(dimuon.eta_2) > muEtaCut) continue;  // muon eta cut
    if(dimuon.q_1 == dimuon.q_2)                                       continue;  // require opposite charge
    
    // average efficiency for specific cut values
    if(dimuon.trkIso03_1 < trkIso03Cut) numerIso03CutBkg++;
    if(dimuon.trkIso03_2 < trkIso03Cut) numerIso03CutBkg++;
    if(dimuon.trkIso03_1/dimuon.pt_1 < trkRelIso03Cut) numerRelIso03CutBkg++;
    if(dimuon.trkIso03_2/dimuon.pt_2 < trkRelIso03Cut) numerRelIso03CutBkg++;
    
    if(dimuon.trkIso05_1 < trkIso05Cut) numerIso05CutBkg++;
    if(dimuon.trkIso05_2 < trkIso05Cut) numerIso05CutBkg++;
    if(dimuon.trkIso05_1/dimuon.pt_1 < trkRelIso05Cut) numerRelIso05CutBkg++;
    if(dimuon.trkIso05_2/dimuon.pt_2 < trkRelIso05Cut) numerRelIso05CutBkg++;
    
    // scan isolation cut values
    denomBkg+=2;        
    for(Int_t istep=0; istep<nIso03Steps; istep++) {
      if(dimuon.trkIso03_1 < iso03Low+istep*iso03Step) numerIso03Bkg[istep]++;
      if(dimuon.trkIso03_2 < iso03Low+istep*iso03Step) numerIso03Bkg[istep]++;
    }
    for(Int_t istep=0; istep<nIso05Steps; istep++) {
      if(dimuon.trkIso05_1 < iso05Low+istep*iso05Step) numerIso05Bkg[istep]++;
      if(dimuon.trkIso05_2 < iso05Low+istep*iso05Step) numerIso05Bkg[istep]++;
    }
    for(Int_t istep=0; istep<nRelIso03Steps; istep++) {
      if(dimuon.trkIso03_1/dimuon.pt_1 < relIso03Low+istep*relIso03Step) numerRelIso03Bkg[istep]++;
      if(dimuon.trkIso03_2/dimuon.pt_2 < relIso03Low+istep*relIso03Step) numerRelIso03Bkg[istep]++;
    }
    for(Int_t istep=0; istep<nRelIso05Steps; istep++) {
      if(dimuon.trkIso05_1/dimuon.pt_1 < relIso05Low+istep*relIso05Step) numerRelIso05Bkg[istep]++;
      if(dimuon.trkIso05_2/dimuon.pt_2 < relIso05Low+istep*relIso05Step) numerRelIso05Bkg[istep]++;
    }
    
    // pT dependence
    for(UInt_t ibin=0; ibin<nPtBins; ibin++) {
      if(dimuon.pt_1 > ptBinEdges[ibin] && dimuon.pt_1 < ptBinEdges[ibin+1]) {
        denomPtBkg[ibin]++;
	if(dimuon.trkIso03_1 < trkIso03Cut) numerIso03vPtBkg[ibin]++;
	if(dimuon.trkIso05_1 < trkIso05Cut) numerIso05vPtBkg[ibin]++;
	if(dimuon.trkIso03_1/dimuon.pt_1 < trkRelIso03Cut) numerRelIso03vPtBkg[ibin]++;
	if(dimuon.trkIso05_1/dimuon.pt_1 < trkRelIso05Cut) numerRelIso05vPtBkg[ibin]++;
      }
      if(dimuon.pt_2 > ptBinEdges[ibin] && dimuon.pt_2 < ptBinEdges[ibin+1]) {
        denomPtBkg[ibin]++;
	if(dimuon.trkIso03_2 < trkIso03Cut) numerIso03vPtBkg[ibin]++;
	if(dimuon.trkIso05_2 < trkIso05Cut) numerIso05vPtBkg[ibin]++;
	if(dimuon.trkIso03_2/dimuon.pt_2 < trkRelIso03Cut) numerRelIso03vPtBkg[ibin]++;
	if(dimuon.trkIso05_2/dimuon.pt_2 < trkRelIso05Cut) numerRelIso05vPtBkg[ibin]++;
      }
    }
    
    // eta dependence
    for(UInt_t ibin=0; ibin<nEtaBins; ibin++) {
      if(dimuon.eta_1 > etaBinEdges[ibin] && dimuon.eta_1 < etaBinEdges[ibin+1]) {
        denomEtaBkg[ibin]++;
	if(dimuon.trkIso03_1 < trkIso03Cut) numerIso03vEtaBkg[ibin]++;
	if(dimuon.trkIso05_1 < trkIso05Cut) numerIso05vEtaBkg[ibin]++;
	if(dimuon.trkIso03_1/dimuon.pt_1 < trkRelIso03Cut) numerRelIso03vEtaBkg[ibin]++;
	if(dimuon.trkIso05_1/dimuon.pt_1 < trkRelIso05Cut) numerRelIso05vEtaBkg[ibin]++;
      }
      if(dimuon.eta_2 > etaBinEdges[ibin] && dimuon.eta_2 < etaBinEdges[ibin+1]) {
        denomEtaBkg[ibin]++;
	if(dimuon.trkIso03_2 < trkIso03Cut) numerIso03vEtaBkg[ibin]++;
	if(dimuon.trkIso05_2 < trkIso05Cut) numerIso05vEtaBkg[ibin]++;
	if(dimuon.trkIso03_2/dimuon.pt_2 < trkRelIso03Cut) numerRelIso03vEtaBkg[ibin]++;
	if(dimuon.trkIso05_2/dimuon.pt_2 < trkRelIso05Cut) numerRelIso05vEtaBkg[ibin]++;
      }
    }
    
    // phi dependence
    for(UInt_t ibin=0; ibin<nPhiBins; ibin++) {
      if(dimuon.phi_1 > phiBinEdges[ibin] && dimuon.phi_1 < phiBinEdges[ibin+1]) {
        denomPhiBkg[ibin]++;
	if(dimuon.trkIso03_1 < trkIso03Cut) numerIso03vPhiBkg[ibin]++;
	if(dimuon.trkIso05_1 < trkIso05Cut) numerIso05vPhiBkg[ibin]++;
	if(dimuon.trkIso03_1/dimuon.pt_1 < trkRelIso03Cut) numerRelIso03vPhiBkg[ibin]++;
	if(dimuon.trkIso05_1/dimuon.pt_1 < trkRelIso05Cut) numerRelIso05vPhiBkg[ibin]++;
      }
      if(dimuon.phi_2 > phiBinEdges[ibin] && dimuon.phi_2 < phiBinEdges[ibin+1]) {
        denomPhiBkg[ibin]++;
	if(dimuon.trkIso03_2 < trkIso03Cut) numerIso03vPhiBkg[ibin]++;
	if(dimuon.trkIso05_2 < trkIso05Cut) numerIso05vPhiBkg[ibin]++;
	if(dimuon.trkIso03_2/dimuon.pt_2 < trkRelIso03Cut) numerRelIso03vPhiBkg[ibin]++;
	if(dimuon.trkIso05_2/dimuon.pt_2 < trkRelIso05Cut) numerRelIso05vPhiBkg[ibin]++;
      }
    }
  }
  delete infile;
  infile=0, infoTree=0, genTree=0, dimuonTree=0, muonTree=0; 
  
  // generate efficiency graphs scan of isolation cut values
  for(Int_t istep=0; istep<nIso03Steps; istep++) {
    iso03EffSig[istep] = toolbox::calcEff(numerIso03Sig[istep],denomSig,&(iso03ErrlSig[istep]),&(iso03ErrhSig[istep]),method);
    iso03EffBkg[istep] = toolbox::calcEff(numerIso03Bkg[istep],denomBkg,&(iso03ErrlBkg[istep]),&(iso03ErrhBkg[istep]),method);
  }
  for(Int_t istep=0; istep<nIso05Steps; istep++) {
    iso05EffSig[istep] = toolbox::calcEff(numerIso05Sig[istep],denomSig,&(iso05ErrlSig[istep]),&(iso05ErrhSig[istep]),method);
    iso05EffBkg[istep] = toolbox::calcEff(numerIso05Bkg[istep],denomBkg,&(iso05ErrlBkg[istep]),&(iso05ErrhBkg[istep]),method);
  }
  TGraph *grIso03Scan = new TGraphErrors(nIso03Steps,iso03EffBkg,iso03EffSig);
  TGraph *grIso05Scan = new TGraphErrors(nIso05Steps,iso05EffBkg,iso05EffSig);
  
  for(Int_t istep=0; istep<nRelIso03Steps; istep++) {
    relIso03EffSig[istep] = toolbox::calcEff(numerRelIso03Sig[istep],denomSig,&(relIso03ErrlSig[istep]),&(relIso03ErrhSig[istep]),method);
    relIso03EffBkg[istep] = toolbox::calcEff(numerRelIso03Bkg[istep],denomBkg,&(relIso03ErrlBkg[istep]),&(relIso03ErrhBkg[istep]),method);
  }
  for(Int_t istep=0; istep<nRelIso05Steps; istep++) {
    relIso05EffSig[istep] = toolbox::calcEff(numerRelIso05Sig[istep],denomSig,&(relIso05ErrlSig[istep]),&(relIso05ErrhSig[istep]),method);
    relIso05EffBkg[istep] = toolbox::calcEff(numerRelIso05Bkg[istep],denomBkg,&(relIso05ErrlBkg[istep]),&(relIso05ErrhBkg[istep]),method);
  }
  TGraphErrors *grRelIso03Scan = new TGraphErrors(nRelIso03Steps,relIso03EffBkg,relIso03EffSig);
  TGraphErrors *grRelIso05Scan = new TGraphErrors(nRelIso05Steps,relIso05EffBkg,relIso05EffSig);
  
  iso03CutEffSig = toolbox::calcEff(numerIso03CutSig,denomSig);
  iso03CutEffBkg = toolbox::calcEff(numerIso03CutBkg,denomBkg);
  TGraph *grIso03Point = new TGraphErrors(1,&iso03CutEffBkg,&iso03CutEffSig);
  iso05CutEffSig = toolbox::calcEff(numerIso05CutSig,denomSig);
  iso05CutEffBkg = toolbox::calcEff(numerIso05CutBkg,denomBkg);
  TGraph *grIso05Point = new TGraphErrors(1,&iso05CutEffBkg,&iso05CutEffSig);
  relIso03CutEffSig = toolbox::calcEff(numerRelIso03CutSig,denomSig);
  relIso03CutEffBkg = toolbox::calcEff(numerRelIso03CutBkg,denomBkg);
  TGraph *grRelIso03Point = new TGraphErrors(1,&relIso03CutEffBkg,&relIso03CutEffSig);
  relIso05CutEffSig = toolbox::calcEff(numerRelIso05CutSig,denomSig);
  relIso05CutEffBkg = toolbox::calcEff(numerRelIso05CutBkg,denomBkg);
  TGraph *grRelIso05Point = new TGraphErrors(1,&relIso05CutEffBkg,&relIso05CutEffSig);
 
  //--------------------------------------------------------------------------------------------------------------
  // Make plots 
  //==============================================================================================================  
  TCanvas *c = MakeCanvas("c","c",canw,canh);
  char label[50];
    
  //
  // Plots for isolation cut scan 
  //    
  CPlot plotScanIso("scanIso","","QCD muon isolation efficiency","Z#rightarrow#mu#mu muon isolation efficiency");
  plotScanIso.AddGraph(grIso03Scan,"#DeltaR < 0.3","C");
  plotScanIso.AddGraph(grIso05Scan,"#DeltaR < 0.5","C",kGreen+2,kFullTriangleDown);
  sprintf(label,"#DeltaR < 0.3, #Sigma^{}p_{T} < %i GeV/c",(Int_t)trkIso03Cut);
  plotScanIso.AddGraph(grIso03Point,label,"",kRed,kFullSquare);
  sprintf(label,"#DeltaR < 0.5, #Sigma^{}p_{T} < %i GeV/c",(Int_t)trkIso05Cut);
  plotScanIso.AddGraph(grIso05Point,label,"",kBlue,kFullSquare);
  plotScanIso.AddTextBox("#Sigma^{}p_{T}",0.21,0.87,0.35,0.8,0);
  plotScanIso.SetYRange(0.9,1.0);
  plotScanIso.SetXRange(0.05,0.5);
  plotScanIso.TransLegend(0,-0.4);
  plotScanIso.Draw(c,doSave,format);
  
  CPlot plotScanRelIso("scanRelIso","","QCD muon isolation efficiency","Z#rightarrow#mu#mu muon isolation efficiency");
  plotScanRelIso.AddGraph(grRelIso03Scan,"#DeltaR < 0.3","C",kCyan+2,kFullTriangleUp);
  plotScanRelIso.AddGraph(grRelIso05Scan,"#DeltaR < 0.5","C",kMagenta+2,kOpenDiamond);
  sprintf(label,"#DeltaR < 0.3, #Sigma^{}p_{T}/^{}p^{#mu}_{T} < %.2f",trkRelIso03Cut);
  plotScanRelIso.AddGraph(grRelIso03Point,label,"",kOrange+7,kFullSquare);
  sprintf(label,"#DeltaR < 0.5, #Sigma^{}p_{T}/^{}p^{#mu}_{T} < %.2f",trkRelIso05Cut);
  plotScanRelIso.AddGraph(grRelIso05Point,label,"",kGray+3,kFullSquare);
  plotScanRelIso.AddTextBox("#Sigma^{}p_{T}/^{}p^{#mu}_{T}",0.21,0.87,0.35,0.8,0);
  plotScanRelIso.SetYRange(0.9,1.0);
  plotScanRelIso.SetXRange(0.05,0.5);
  plotScanRelIso.TransLegend(0,-0.4);
  plotScanRelIso.Draw(c,doSave,format);
  
  CPlot plotScanIso03("scanIso03","","QCD muon isolation efficiency","Z#rightarrow#mu#mu muon isolation efficiency");
  plotScanIso03.AddGraph(grIso03Scan,"#Sigma^{}p_{T}","C");
  plotScanIso03.AddGraph(grRelIso03Scan,"#Sigma^{}p_{T}/^{}p^{#mu}_{T}","C",kCyan+2,kFullTriangleUp);
  plotScanIso03.AddTextBox("#DeltaR < 0.3",0.21,0.87,0.35,0.8,0);
  sprintf(label,"#DeltaR < 0.3, #Sigma^{}p_{T} < %i GeV/c",(Int_t)trkIso03Cut);
  plotScanIso03.AddGraph(grIso03Point,label,"",kRed,kFullSquare);
  sprintf(label,"#DeltaR < 0.3, #Sigma^{}p_{T}/^{}p^{#mu}_{T} < %.2f",trkRelIso03Cut);
  plotScanIso03.AddGraph(grRelIso03Point,label,"",kOrange+7,kFullSquare);
  plotScanIso03.SetYRange(0.9,1.0);
  plotScanIso03.SetXRange(0.05,0.5);
  plotScanIso03.TransLegend(0,-0.4);
  plotScanIso03.Draw(c,doSave,format);
  
  CPlot plotScanIso05("scanIso05","","QCD muon isolation efficiency","Z#rightarrow#mu#mu muon isolation efficiency");
  plotScanIso05.AddGraph(grIso05Scan,"#Sigma^{}p_{T}","C",kGreen+2,kFullTriangleDown);
  plotScanIso05.AddGraph(grRelIso05Scan,"#Sigma^{}p_{T}/^{}p^{#mu}_{T}","C",kMagenta+2,kOpenDiamond);
  plotScanIso05.AddTextBox("#DeltaR < 0.5",0.21,0.87,0.35,0.8,0);
  sprintf(label,"#DeltaR < 0.5, #Sigma^{}p_{T} < %i GeV/c",(Int_t)trkIso05Cut);
  plotScanIso05.AddGraph(grIso05Point,label,"",kBlue,kFullSquare);
  sprintf(label,"#DeltaR < 0.5, #Sigma^{}p_{T}/^{}p^{#mu}_{T} < %.2f",trkRelIso05Cut);
  plotScanIso05.AddGraph(grRelIso05Point,label,"",kGray+3,kFullSquare);
  plotScanIso05.SetYRange(0.9,1.0);
  plotScanIso05.SetXRange(0.05,0.5);
  plotScanIso05.TransLegend(0,-0.4);
  plotScanIso05.Draw(c,doSave,format);
  
  CPlot plotScanIsoAll("scanIsoAll","","QCD muon isolation efficiency","Z#rightarrow#mu#mu muon isolation efficiency");
  plotScanIsoAll.AddGraph(grIso03Scan,"#DeltaR < 0.3, #Sigma^{}p_{T}","C");
  plotScanIsoAll.AddGraph(grIso05Scan,"#DeltaR < 0.5, #Sigma^{}p_{T}","C",kGreen+2,kFullTriangleDown);
  plotScanIsoAll.AddGraph(grRelIso03Scan,"#DeltaR < 0.3, #Sigma^{}p_{T}/^{}p^{#mu}_{T}","C",kCyan+2,kFullTriangleUp);
  plotScanIsoAll.AddGraph(grRelIso05Scan,"#DeltaR < 0.5, #Sigma^{}p_{T}/^{}p^{#mu}_{T}","C",kMagenta+2,kOpenDiamond);
  sprintf(label,"#DeltaR < 0.3, #Sigma^{}p_{T} < %i GeV/c",(Int_t)trkIso03Cut);
  plotScanIsoAll.AddGraph(grIso03Point,label,"",kRed,kFullSquare);
  plotScanIsoAll.SetYRange(0.9,1.0);
  plotScanIsoAll.SetXRange(0.05,0.5);
  plotScanIsoAll.TransLegend(0,-0.4);
  plotScanIsoAll.Draw(c,doSave,format);  
  
  //
  // Generate graphs of efficiency vs. muon pT
  //
  TGraphAsymmErrors *grIso03vPtSig = toolbox::makeEffGraph(nPtBins,ptBinEdges,numerIso03vPtSig,denomPtSig,method);
  TGraphAsymmErrors *grIso03vPtBkg = toolbox::makeEffGraph(nPtBins,ptBinEdges,numerIso03vPtBkg,denomPtBkg,method);
  TGraphAsymmErrors *grIso05vPtSig = toolbox::makeEffGraph(nPtBins,ptBinEdges,numerIso05vPtSig,denomPtSig,method);
  TGraphAsymmErrors *grIso05vPtBkg = toolbox::makeEffGraph(nPtBins,ptBinEdges,numerIso05vPtBkg,denomPtBkg,method);

  TGraphAsymmErrors *grRelIso03vPtSig = toolbox::makeEffGraph(nPtBins,ptBinEdges,numerRelIso03vPtSig,denomPtSig,method);
  TGraphAsymmErrors *grRelIso03vPtBkg = toolbox::makeEffGraph(nPtBins,ptBinEdges,numerRelIso03vPtBkg,denomPtBkg,method);
  TGraphAsymmErrors *grRelIso05vPtSig = toolbox::makeEffGraph(nPtBins,ptBinEdges,numerRelIso05vPtSig,denomPtSig,method);
  TGraphAsymmErrors *grRelIso05vPtBkg = toolbox::makeEffGraph(nPtBins,ptBinEdges,numerRelIso05vPtBkg,denomPtBkg,method);
  
  CPlot plotIsoEffvPt("iso03EffvPt","#DeltaR < 0.3","muon p_{T} [GeV/c]","isolation efficiency");
  plotIsoEffvPt.AddGraph(grIso03vPtSig,"#Sigma^{}p_{T} < 3 GeV/c, Z#rightarrow#mu#mu","",kBlack);
  plotIsoEffvPt.AddGraph(grRelIso03vPtSig,"#Sigma^{}p_{T}/^{}p^{#mu}_{T} < 0.09, Z#rightarrow#mu#mu","",kRed,kFullTriangleUp);
  plotIsoEffvPt.AddGraph(grIso03vPtBkg,"#Sigma^{}p_{T} < 3 GeV/c, QCD","",kBlue,kFullSquare);
  plotIsoEffvPt.AddGraph(grRelIso03vPtBkg,"#Sigma^{}p_{T}/^{}p^{#mu}_{T} < 0.09, QCD","",kGreen+2,kFullTriangleDown);
  plotIsoEffvPt.SetYRange(0,1.01);
  plotIsoEffvPt.TransLegend(0,-0.1);
  plotIsoEffvPt.Draw(c,doSave,format);
  
  CPlot plotIso03EffvPtSig("iso03EffvPtSig","#DeltaR < 0.3","muon p_{T} [GeV/c]","isolation efficiency");
  plotIso03EffvPtSig.AddGraph(grIso03vPtSig,"#Sigma^{}p_{T} < 3 GeV/c, Z#rightarrow#mu#mu","",kBlack);
  plotIso03EffvPtSig.AddGraph(grRelIso03vPtSig,"#Sigma^{}p_{T}/^{}p^{#mu}_{T} < 0.09, Z#rightarrow#mu#mu","",kRed,kFullTriangleUp);
  plotIso03EffvPtSig.SetYRange(0.93,1.0);
  plotIso03EffvPtSig.TransLegend(0,-0.6);
  plotIso03EffvPtSig.Draw(c,doSave,format);  
plotIso03EffvPtSig.Draw(c,true,"C");

  CPlot plotIso03EffvPtBkg("iso03EffvPtBkg","#DeltaR < 0.3","muon p_{T} [GeV/c]","isolation efficiency");
  plotIso03EffvPtBkg.AddGraph(grIso03vPtBkg,"#Sigma^{}p_{T} < 3 GeV/c, QCD","",kBlue,kFullSquare);
  plotIso03EffvPtBkg.AddGraph(grRelIso03vPtBkg,"#Sigma^{}p_{T}/^{}p^{#mu}_{T} < 0.09, QCD","",kGreen+2,kFullTriangleDown);
  plotIso03EffvPtBkg.SetYRange(0,0.5);
  plotIso03EffvPtBkg.TransLegend(0,-0.6);
  plotIso03EffvPtBkg.Draw(c,doSave,format);

  CPlot plotIso05EffvPt("iso05EffvPt","#DeltaR < 0.5","muon p_{T} [GeV/c]","isolation efficiency");
  plotIso05EffvPt.AddGraph(grIso05vPtSig,"#Sigma^{}p_{T} < 3 GeV/c, Z#rightarrow#mu#mu","",kBlack);
  plotIso05EffvPt.AddGraph(grRelIso05vPtSig,"#Sigma^{}p_{T}/^{}p^{#mu}_{T} < 0.09, Z#rightarrow#mu#mu","",kRed,kFullTriangleUp);
  plotIso05EffvPt.AddGraph(grIso05vPtBkg,"#Sigma^{}p_{T} < 3 GeV/c, QCD","",kBlue,kFullSquare);
  plotIso05EffvPt.AddGraph(grRelIso05vPtBkg,"#Sigma^{}p_{T}/^{}p^{#mu}_{T} < 0.09, QCD","",kGreen+2,kFullTriangleDown);
  plotIso05EffvPt.SetYRange(0,1.01);
  plotIso05EffvPt.TransLegend(0,-0.1);
  plotIso05EffvPt.Draw(c,doSave,format);
  
  CPlot plotIso05EffvPtSig("iso05EffvPtSig","#DeltaR < 0.5","muon p_{T} [GeV/c]","isolation efficiency");
  plotIso05EffvPtSig.AddGraph(grIso05vPtSig,"#Sigma^{}p_{T} < 3 GeV/c, Z#rightarrow#mu#mu","",kBlack);
  plotIso05EffvPtSig.AddGraph(grRelIso05vPtSig,"#Sigma^{}p_{T}/^{}p^{#mu}_{T} < 0.09, Z#rightarrow#mu#mu","",kRed,kFullTriangleUp);
  plotIso05EffvPtSig.SetYRange(0.93,1.0);
  plotIso05EffvPtSig.TransLegend(0,-0.6);
  plotIso05EffvPtSig.Draw(c,doSave,format);  

  CPlot plotIso05EffvPtBkg("iso05EffvPtBkg","#DeltaR < 0.5","muon p_{T} [GeV/c]","isolation efficiency");
  plotIso05EffvPtBkg.AddGraph(grIso05vPtBkg,"#Sigma^{}p_{T} < 3 GeV/c, QCD","",kBlue,kFullSquare);
  plotIso05EffvPtBkg.AddGraph(grRelIso05vPtBkg,"#Sigma^{}p_{T}/^{}p^{#mu}_{T} < 0.09, QCD","",kGreen+2,kFullTriangleDown);
  plotIso05EffvPtBkg.SetYRange(0,0.5);
  plotIso05EffvPtBkg.TransLegend(0,-0.6);
  plotIso05EffvPtBkg.Draw(c,doSave,format);
    
  //
  // Generate graphs of efficiency vs. muon eta
  //
  TGraphAsymmErrors *grIso03vEtaSig = toolbox::makeEffGraph(nEtaBins,etaBinEdges,numerIso03vEtaSig,denomEtaSig,method);
  TGraphAsymmErrors *grIso03vEtaBkg = toolbox::makeEffGraph(nEtaBins,etaBinEdges,numerIso03vEtaBkg,denomEtaBkg,method);
  TGraphAsymmErrors *grIso05vEtaSig = toolbox::makeEffGraph(nEtaBins,etaBinEdges,numerIso05vEtaSig,denomEtaSig,method);
  TGraphAsymmErrors *grIso05vEtaBkg = toolbox::makeEffGraph(nEtaBins,etaBinEdges,numerIso05vEtaBkg,denomEtaBkg,method);

  TGraphAsymmErrors *grRelIso03vEtaSig = toolbox::makeEffGraph(nEtaBins,etaBinEdges,numerRelIso03vEtaSig,denomEtaSig,method);
  TGraphAsymmErrors *grRelIso03vEtaBkg = toolbox::makeEffGraph(nEtaBins,etaBinEdges,numerRelIso03vEtaBkg,denomEtaBkg,method);
  TGraphAsymmErrors *grRelIso05vEtaSig = toolbox::makeEffGraph(nEtaBins,etaBinEdges,numerRelIso05vEtaSig,denomEtaSig,method);
  TGraphAsymmErrors *grRelIso05vEtaBkg = toolbox::makeEffGraph(nEtaBins,etaBinEdges,numerRelIso05vEtaBkg,denomEtaBkg,method);
  
  CPlot plotIso03EffvEta("iso03EffvEta","#DeltaR < 0.3","muon #eta","isolation efficiency");
  plotIso03EffvEta.AddGraph(grIso03vEtaSig,"#Sigma^{}p_{T} < 3 GeV/c, Z#rightarrow#mu#mu","",kBlack);
  plotIso03EffvEta.AddGraph(grRelIso03vEtaSig,"#Sigma^{}p_{T}/^{}p^{#mu}_{T} < 0.09, Z#rightarrow#mu#mu","",kRed,kFullTriangleUp);
  plotIso03EffvEta.AddGraph(grIso03vEtaBkg,"#Sigma^{}p_{T} < 3 GeV/c, QCD","",kBlue,kFullSquare);
  plotIso03EffvEta.AddGraph(grRelIso03vEtaBkg,"#Sigma^{}p_{T}/^{}p^{#mu}_{T} < 0.09, QCD","",kGreen+2,kFullTriangleDown);
  plotIso03EffvEta.SetYRange(0,1.01);
  plotIso03EffvEta.TransLegend(0,-0.1);
  plotIso03EffvEta.Draw(c,doSave,format);
  
  CPlot plotIso03EffvEtaSig("iso03EffvEtaSig","#DeltaR < 0.3","muon #eta","isolation efficiency");
  plotIso03EffvEtaSig.AddGraph(grIso03vEtaSig,"#Sigma^{}p_{T} < 3 GeV/c, Z#rightarrow#mu#mu","",kBlack);
  plotIso03EffvEtaSig.AddGraph(grRelIso03vEtaSig,"#Sigma^{}p_{T}/^{}p^{#mu}_{T} < 0.09, Z#rightarrow#mu#mu","",kRed,kFullTriangleUp);
  plotIso03EffvEtaSig.SetYRange(0.93,1.0);
  plotIso03EffvEtaSig.Draw(c,doSave,format);  
plotIso03EffvEtaSig.Draw(c,true,"C");

  CPlot plotIso03EffvEtaBkg("iso03EffvEtaBkg","#DeltaR < 0.3","muon #eta","isolation efficiency");
  plotIso03EffvEtaBkg.AddGraph(grIso03vEtaBkg,"#Sigma^{}p_{T} < 3 GeV/c, QCD","",kBlue,kFullSquare);
  plotIso03EffvEtaBkg.AddGraph(grRelIso03vEtaBkg,"#Sigma^{}p_{T}/^{}p^{#mu}_{T} < 0.09, QCD","",kGreen+2,kFullTriangleDown);
  plotIso03EffvEtaBkg.SetYRange(0,0.5);
  plotIso03EffvEtaBkg.Draw(c,doSave,format);

  CPlot plotIso05EffvEta("iso05EffvEta","#DeltaR < 0.5","muon #eta","isolation efficiency");
  plotIso05EffvEta.AddGraph(grIso05vEtaSig,"#Sigma^{}p_{T} < 3 GeV/c, Z#rightarrow#mu#mu","",kBlack);
  plotIso05EffvEta.AddGraph(grRelIso05vEtaSig,"#Sigma^{}p_{T}/^{}p^{#mu}_{T} < 0.09, Z#rightarrow#mu#mu","",kRed,kFullTriangleUp);
  plotIso05EffvEta.AddGraph(grIso05vEtaBkg,"#Sigma^{}p_{T} < 3 GeV/c, QCD","",kBlue,kFullSquare);
  plotIso05EffvEta.AddGraph(grRelIso05vEtaBkg,"#Sigma^{}p_{T}/^{}p^{#mu}_{T} < 0.09, QCD","",kGreen+2,kFullTriangleDown);
  plotIso05EffvEta.SetYRange(0,1.01);
  plotIso05EffvEta.TransLegend(0,-0.1);
  plotIso05EffvEta.Draw(c,doSave,format);
  
  CPlot plotIso05EffvEtaSig("iso05EffvEtaSig","#DeltaR < 0.5","muon #eta","isolation efficiency");
  plotIso05EffvEtaSig.AddGraph(grIso05vEtaSig,"#Sigma^{}p_{T} < 3 GeV/c, Z#rightarrow#mu#mu","",kBlack);
  plotIso05EffvEtaSig.AddGraph(grRelIso05vEtaSig,"#Sigma^{}p_{T}/^{}p^{#mu}_{T} < 0.09, Z#rightarrow#mu#mu","",kRed,kFullTriangleUp);
  plotIso05EffvEtaSig.SetYRange(0.93,1.0);
  plotIso05EffvEtaSig.Draw(c,doSave,format);  

  CPlot plotIso05EffvEtaBkg("iso05EffvEtaBkg","#DeltaR < 0.5","muon #eta","isolation efficiency");
  plotIso05EffvEtaBkg.AddGraph(grIso05vEtaBkg,"#Sigma^{}p_{T} < 3 GeV/c, QCD","",kBlue,kFullSquare);
  plotIso05EffvEtaBkg.AddGraph(grRelIso05vEtaBkg,"#Sigma^{}p_{T}/^{}p^{#mu}_{T} < 0.09, QCD","",kGreen+2,kFullTriangleDown);
  plotIso05EffvEtaBkg.SetYRange(0,0.5);
  plotIso05EffvEtaBkg.Draw(c,doSave,format);
  
  //
  // Generate graphs of efficiency vs. muon eta
  //
  TGraphAsymmErrors *grIso03vPhiSig = toolbox::makeEffGraph(nPhiBins,phiBinEdges,numerIso03vPhiSig,denomPhiSig,method);
  TGraphAsymmErrors *grIso03vPhiBkg = toolbox::makeEffGraph(nPhiBins,phiBinEdges,numerIso03vPhiBkg,denomPhiBkg,method);
  TGraphAsymmErrors *grIso05vPhiSig = toolbox::makeEffGraph(nPhiBins,phiBinEdges,numerIso05vPhiSig,denomPhiSig,method);
  TGraphAsymmErrors *grIso05vPhiBkg = toolbox::makeEffGraph(nPhiBins,phiBinEdges,numerIso05vPhiBkg,denomPhiBkg,method);

  TGraphAsymmErrors *grRelIso03vPhiSig = toolbox::makeEffGraph(nPhiBins,phiBinEdges,numerRelIso03vPhiSig,denomPhiSig,method);
  TGraphAsymmErrors *grRelIso03vPhiBkg = toolbox::makeEffGraph(nPhiBins,phiBinEdges,numerRelIso03vPhiBkg,denomPhiBkg,method);
  TGraphAsymmErrors *grRelIso05vPhiSig = toolbox::makeEffGraph(nPhiBins,phiBinEdges,numerRelIso05vPhiSig,denomPhiSig,method);
  TGraphAsymmErrors *grRelIso05vPhiBkg = toolbox::makeEffGraph(nPhiBins,phiBinEdges,numerRelIso05vPhiBkg,denomPhiBkg,method);
      
  CPlot plotIso03EffvPhi("iso03EffvPhi","#DeltaR < 0.3","muon #phi","isolation efficiency");
  plotIso03EffvPhi.AddGraph(grIso03vPhiSig,"#Sigma^{}p_{T} < 3 GeV/c, Z#rightarrow#mu#mu","",kBlack);
  plotIso03EffvPhi.AddGraph(grRelIso03vPhiSig,"#Sigma^{}p_{T}/^{}p^{#mu}_{T} < 0.09, Z#rightarrow#mu#mu","",kRed,kFullTriangleUp);
  plotIso03EffvPhi.AddGraph(grIso03vPhiBkg,"#Sigma^{}p_{T} < 3 GeV/c, QCD","",kBlue,kFullSquare);
  plotIso03EffvPhi.AddGraph(grRelIso03vPhiBkg,"#Sigma^{}p_{T}/^{}p^{#mu}_{T} < 0.09, QCD","",kGreen+2,kFullTriangleDown);
  plotIso03EffvPhi.SetYRange(0,1.01);
  plotIso03EffvPhi.TransLegend(0,-0.1);
  plotIso03EffvPhi.Draw(c,doSave,format);  

  CPlot plotIso03EffvPhiSig("iso03EffvPhiSig","#DeltaR < 0.3","muon #phi","isolation efficiency");
  plotIso03EffvPhiSig.AddGraph(grIso03vPhiSig,"#Sigma^{}p_{T} < 3 GeV/c, Z#rightarrow#mu#mu","",kBlack);
  plotIso03EffvPhiSig.AddGraph(grRelIso03vPhiSig,"#Sigma^{}p_{T}/^{}p^{#mu}_{T} < 0.09, Z#rightarrow#mu#mu","",kRed,kFullTriangleUp);
  plotIso03EffvPhiSig.SetYRange(0.93,1.0);
  plotIso03EffvPhiSig.Draw(c,doSave,format);  
plotIso03EffvPhiSig.Draw(c,true,"C");

  CPlot plotIso03EffvPhiBkg("iso03EffvPhiBkg","#DeltaR < 0.3","muon #phi","isolation efficiency");
  plotIso03EffvPhiBkg.AddGraph(grIso03vPhiBkg,"#Sigma^{}p_{T} < 3 GeV/c, QCD","",kBlue,kFullSquare);
  plotIso03EffvPhiBkg.AddGraph(grRelIso03vPhiBkg,"#Sigma^{}p_{T}/^{}p^{#mu}_{T} < 0.09, QCD","",kGreen+2,kFullTriangleDown);
  plotIso03EffvPhiBkg.SetYRange(0,0.5);
  plotIso03EffvPhiBkg.Draw(c,doSave,format);

  CPlot plotIso05EffvPhi("iso05EffvPhi","#DeltaR < 0.5","muon #phi","isolation efficiency");
  plotIso05EffvPhi.AddGraph(grIso05vPhiSig,"#Sigma^{}p_{T} < 3 GeV/c, Z#rightarrow#mu#mu","",kBlack);
  plotIso05EffvPhi.AddGraph(grRelIso05vPhiSig,"#Sigma^{}p_{T}/^{}p^{#mu}_{T} < 0.09, Z#rightarrow#mu#mu","",kRed,kFullTriangleUp);
  plotIso05EffvPhi.AddGraph(grIso05vPhiBkg,"#Sigma^{}p_{T} < 3 GeV/c, QCD","",kBlue,kFullSquare);
  plotIso05EffvPhi.AddGraph(grRelIso05vPhiBkg,"#Sigma^{}p_{T}/^{}p^{#mu}_{T} < 0.09, QCD","",kGreen+2,kFullTriangleDown);
  plotIso05EffvPhi.SetYRange(0,1.01);
  plotIso05EffvPhi.TransLegend(0,-0.1);
  plotIso05EffvPhi.Draw(c,doSave,format);  

  CPlot plotIso05EffvPhiSig("iso05EffvPhiSig","#DeltaR < 0.5","muon #phi","isolation efficiency");
  plotIso05EffvPhiSig.AddGraph(grIso05vPhiSig,"#Sigma^{}p_{T} < 3 GeV/c, Z#rightarrow#mu#mu","",kBlack);
  plotIso05EffvPhiSig.AddGraph(grRelIso05vPhiSig,"#Sigma^{}p_{T}/^{}p^{#mu}_{T} < 0.09, Z#rightarrow#mu#mu","",kRed,kFullTriangleUp);
  plotIso05EffvPhiSig.SetYRange(0.93,1.0);
  plotIso05EffvPhiSig.Draw(c,doSave,format);  

  CPlot plotIso05EffvPhiBkg("iso05EffvPhiBkg","#DeltaR < 0.5","muon #phi","isolation efficiency");
  plotIso05EffvPhiBkg.AddGraph(grIso05vPhiBkg,"#Sigma^{}p_{T} < 3 GeV/c, QCD","",kBlue,kFullSquare);
  plotIso05EffvPhiBkg.AddGraph(grRelIso05vPhiBkg,"#Sigma^{}p_{T}/^{}p^{#mu}_{T} < 0.09, QCD","",kGreen+2,kFullTriangleDown);
  plotIso05EffvPhiBkg.SetYRange(0,0.5);
  plotIso05EffvPhiBkg.Draw(c,doSave,format);
    
  //--------------------------------------------------------------------------------------------------------------
  // Summary print out
  //==============================================================================================================
  cout << endl;
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  cout << endl;
  
  cout << "  Signal sample: " << denomSig << " events considered" << endl;
  cout << "  Background sample: " << denomBkg << " events considered" << endl;
  cout << endl;
  
  cout << " dR = 0.3 " << endl;
  cout << "----------" << endl;
  for(Int_t istep=0; istep<nIso03Steps; istep++) {
    cout << "dR = 0.3, iso cut = " << setprecision(1) << fixed << iso03Low+istep*iso03Step << ": ";
    cout << "  Sig Eff = " << setprecision(4) << fixed << iso03EffSig[istep]; 
    cout << " [-" << setprecision(4) << fixed << iso03ErrlSig[istep];
    cout << ", +" << setprecision(4) << fixed << iso03ErrhSig[istep] << "]; ";
    cout << "  Bkg Eff = " << setprecision(4) << fixed << iso03EffBkg[istep];
    cout << " [-" << setprecision(4) << fixed << iso03ErrlBkg[istep];
    cout << ", +" << setprecision(4) << fixed << iso03ErrhBkg[istep] << "]" << endl;
  }
  cout << endl;
  
  cout << " dR = 0.5 " << endl;
  cout << "----------" << endl;
  for(Int_t istep=0; istep<nIso05Steps; istep++) {
    cout << "iso cut = " << setprecision(1) << fixed << iso05Low+istep*iso05Step << ": ";
    cout << "  Sig Eff = " << setprecision(4) << fixed << iso05EffSig[istep]; 
    cout << " [-" << setprecision(4) << fixed << iso05ErrlSig[istep];
    cout << ", +" << setprecision(4) << fixed << iso05ErrhSig[istep] << "]; ";
    cout << "  Bkg Eff = " << setprecision(4) << fixed << iso05EffBkg[istep];
    cout << " [-" << setprecision(4) << fixed << iso05ErrlBkg[istep];
    cout << ", +" << setprecision(4) << fixed << iso05ErrhBkg[istep] << "]" << endl;
  }
  cout << endl;
  
  cout << " dR = 0.3 " << endl;
  cout << "----------" << endl;    
  for(Int_t istep=0; istep<nRelIso03Steps; istep++) {
    cout << "rel iso cut = " << setprecision(4) << relIso03Low+istep*relIso03Step << ": ";
    cout << "  Sig Eff = " << setprecision(4) << fixed << relIso03EffSig[istep]; 
    cout << " [-" << setprecision(4) << fixed << relIso03ErrlSig[istep];
    cout << ", +" << setprecision(4) << fixed << relIso03ErrhSig[istep] << "]; ";
    cout << "  Bkg Eff = " << setprecision(4) << fixed << relIso03EffBkg[istep];
    cout << " [-" << setprecision(4) << fixed << relIso03ErrlBkg[istep];
    cout << ", +" << setprecision(4) << fixed << relIso03ErrhBkg[istep] << "]" << endl;
  }
  cout << endl;

  cout << " dR = 0.5 " << endl;
  cout << "----------" << endl;  
  for(Int_t istep=0; istep<nRelIso05Steps; istep++) {
    cout << "rel iso cut = " << setprecision(4) << relIso05Low+istep*relIso05Step << ": ";
    cout << "  Sig Eff = " << setprecision(4) << fixed << relIso05EffSig[istep]; 
    cout << " [-" << setprecision(4) << fixed << relIso05ErrlSig[istep];
    cout << ", +" << setprecision(4) << fixed << relIso05ErrhSig[istep] << "]; ";
    cout << "  Bkg Eff = " << setprecision(4) << fixed << relIso05EffBkg[istep];
    cout << " [-" << setprecision(4) << fixed << relIso05ErrlBkg[istep];
    cout << ", +" << setprecision(4) << fixed << relIso05ErrhBkg[istep] << "]" << endl;
  }
  cout << endl;
  
  gBenchmark->Show("plotIsoEffMC");
}

