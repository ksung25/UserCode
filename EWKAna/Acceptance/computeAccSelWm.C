//================================================================================================
//
// Compute W->munu acceptance at full selection level
//
//  * outputs results summary text file
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TClonesArray.h>           // ROOT array class
#include <TBenchmark.h>             // class to track macro running statistics
#include <TH1D.h>                   // histogram class
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings
#include "Math/LorentzVector.h"     // 4-vector class

// define structures to read in ntuple
#include "EWKAna/Ntupler/interface/EWKAnaDefs.hh"
#include "EWKAna/Ntupler/interface/TEventInfo.hh"
#include "EWKAna/Ntupler/interface/TGenInfo.hh"
#include "EWKAna/Ntupler/interface/TMuon.hh"

// helper functions for lepton ID selection
#include "EWKAna/Utils/LeptonIDCuts.hh"

// helper class to handle efficiency tables
#include "CEffUser2D.hh"
#endif

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVector;


//=== MAIN MACRO ================================================================================================= 

void computeAccSelWm(const TString conf,       // input file
                     const TString outputDir,  // output directory
		     const Int_t   charge,     // 0 = inclusive, +1 = W+, -1 = W-
		     const Int_t   doPU=0      // PU-reweight scheme (0 = no reweight, 1 = 2011A, 2 = 2011B, 3 = All 2011)
) {
  gBenchmark->Start("computeAccSelWm");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 

  const Double_t PT_CUT     = 25;
  const Double_t ETA_CUT    = 2.1;
  const Double_t ETA_BARREL = 1.2;
  const Double_t ETA_ENDCAP = 1.2;
  
  // pile-up weight files
  TString pufname("");
  if(doPU==1) pufname = "/data/blue/ksung/EWKAna/test/Utils/PileupReweighting.Summer11DYmm_To_Run2011A.root";
  if(doPU==2) pufname = "/data/blue/ksung/EWKAna/test/Utils/PileupReweighting.Summer11DYmm_To_Run2011B.root";
  if(doPU==3) pufname = "/data/blue/ksung/EWKAna/test/Utils/PileupReweighting.Summer11DYmm_To_Full2011.root";
  
  // efficiency files
  TString dataHLTEffName("../Efficiency/Data_MuHLTEff/analysis/eff.root");
  TString zmmHLTEffName("../Efficiency/Zmm_MuHLTEff/analysis/eff.root");
  TString dataSelEffName("../Efficiency/Data_MuSelEff/analysis/eff.root");
  TString zmmSelEffName("../Efficiency/Zmm_MuSelEff/analysis/eff.root");
  TString dataTrkEffName("../Efficiency/Data_MuTrkEff/analysis/eff.root");
  TString zmmTrkEffName("../Efficiency/Zmm_MuTrkEff/analysis/eff.root");
  TString dataStaEffName("../Efficiency/Data_MuStaEff/analysis/eff.root");
  TString zmmStaEffName("../Efficiency/Zmm_MuStaEff/analysis/eff.root");
  

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  

  vector<TString> fnamev;  // file name per input file
  vector<TString> labelv;  // TLegend label per input file
  vector<Int_t>   colorv;  // plot color per input file
  vector<Int_t>   linev;   // plot line style per input file

  //
  // parse .conf file
  //
  ifstream ifs;
  ifs.open(conf.Data());
  assert(ifs.is_open());
  string line;
  while(getline(ifs,line)) {
    if(line[0]=='#') continue;
    
    string fname;
    Int_t color, linesty;
    stringstream ss(line);
    ss >> fname >> color >> linesty;
    string label = line.substr(line.find('@')+1);
    fnamev.push_back(fname);
    labelv.push_back(label);
    colorv.push_back(color);
    linev.push_back(linesty);
  }
  ifs.close();

  // Create output directory
  gSystem->mkdir(outputDir,kTRUE);
  
  
  // Get pile-up weights
  TFile *pufile=0;
  TH1D  *puWeights=0;
  if(doPU>0) {
    pufile    = new TFile(pufname);	         assert(pufile);
    puWeights = (TH1D*)pufile->Get("puWeights"); assert(puWeights);
  }
  
  //
  // Get efficiency
  //
  TFile *dataHLTEffFile = new TFile(dataHLTEffName);
  CEffUser2D dataHLTEff;
  if(dataHLTEffName) {    
    dataHLTEff.loadEff((TH2D*)dataHLTEffFile->Get("hEffEtaPt"), 
                       (TH2D*)dataHLTEffFile->Get("hErrlEtaPt"),
		       (TH2D*)dataHLTEffFile->Get("hErrhEtaPt"));
  }
  
  TFile *zmmHLTEffFile = new TFile(zmmHLTEffName);
  CEffUser2D zmmHLTEff;
  if(zmmHLTEffName) {
    zmmHLTEff.loadEff((TH2D*)zmmHLTEffFile->Get("hEffEtaPt"),
                      (TH2D*)zmmHLTEffFile->Get("hErrlEtaPt"),
                      (TH2D*)zmmHLTEffFile->Get("hErrhEtaPt"));
  }
  
  TFile *dataSelEffFile = new TFile(dataSelEffName);
  CEffUser2D dataSelEff;
  if(dataSelEffName) {
    dataSelEff.loadEff((TH2D*)dataSelEffFile->Get("hEffEtaPt"),
                       (TH2D*)dataSelEffFile->Get("hErrlEtaPt"),
                       (TH2D*)dataSelEffFile->Get("hErrhEtaPt"));
  }
  
  TFile *zmmSelEffFile = new TFile(zmmSelEffName);
  CEffUser2D zmmSelEff;
  if(zmmSelEffName) {
    zmmSelEff.loadEff((TH2D*)zmmSelEffFile->Get("hEffEtaPt"),
                      (TH2D*)zmmSelEffFile->Get("hErrlEtaPt"),
                      (TH2D*)zmmSelEffFile->Get("hErrhEtaPt"));
  }
  
  TFile *dataTrkEffFile = new TFile(dataTrkEffName);
  CEffUser2D dataTrkEff;
  if(dataTrkEffName) {
    dataTrkEff.loadEff((TH2D*)dataTrkEffFile->Get("hEffEtaPt"),
                       (TH2D*)dataTrkEffFile->Get("hErrlEtaPt"),
                       (TH2D*)dataTrkEffFile->Get("hErrhEtaPt"));
  }
  
  TFile *zmmTrkEffFile = new TFile(zmmTrkEffName);
  CEffUser2D zmmTrkEff;
  if(zmmTrkEffName) {
    zmmTrkEff.loadEff((TH2D*)zmmTrkEffFile->Get("hEffEtaPt"),
                      (TH2D*)zmmTrkEffFile->Get("hErrlEtaPt"),
                      (TH2D*)zmmTrkEffFile->Get("hErrhEtaPt"));
  }
  
  TFile *dataStaEffFile = new TFile(dataStaEffName);
  CEffUser2D dataStaEff;
  if(dataStaEffName) {
    dataStaEff.loadEff((TH2D*)dataStaEffFile->Get("hEffEtaPt"),
                       (TH2D*)dataStaEffFile->Get("hErrlEtaPt"),
                       (TH2D*)dataStaEffFile->Get("hErrhEtaPt"));
  }
  
  TFile *zmmStaEffFile = new TFile(zmmStaEffName);
  CEffUser2D zmmStaEff;
  if(zmmStaEffName) {
    zmmStaEff.loadEff((TH2D*)zmmStaEffFile->Get("hEffEtaPt"),
                      (TH2D*)zmmStaEffFile->Get("hErrlEtaPt"),
                      (TH2D*)zmmStaEffFile->Get("hErrhEtaPt"));
  }
  
  // Data structures to store info from TTrees
  mithep::TEventInfo *info = new mithep::TEventInfo();
  mithep::TGenInfo   *gen  = new mithep::TGenInfo();
  TClonesArray *muonArr    = new TClonesArray("mithep::TMuon");
  
  TFile *infile=0;
  TTree *eventTree=0;

  // Variables to store acceptances and uncertainties (per input file)
  vector<Double_t> nEvtsv, nSelv, nSelBv, nSelEv;
  vector<Double_t> accv, accBv, accEv;
  vector<Double_t> accErrv, accErrBv, accErrEv;
  vector<Double_t> nSelCorrv, nSelBCorrv, nSelECorrv;
  vector<Double_t> accCorrv, accBCorrv, accECorrv;
  vector<Double_t> accErrCorrv, accErrBCorrv, accErrECorrv;
  
  //
  // loop through files
  //
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {  

    // Read input file and get the TTrees
    cout << "Processing " << fnamev[ifile] << " ..." << endl;
    infile = new TFile(fnamev[ifile]); 
    assert(infile);
  
    eventTree = (TTree*)infile->Get("Events"); assert(eventTree);  
    eventTree->SetBranchAddress("Info", &info);    TBranch *infoBr = eventTree->GetBranch("Info");
    eventTree->SetBranchAddress("Gen",  &gen);     TBranch *genBr  = eventTree->GetBranch("Gen");
    eventTree->SetBranchAddress("Muon", &muonArr); TBranch *muonBr = eventTree->GetBranch("Muon"); 

    nEvtsv.push_back(0);
    nSelv.push_back(0);
    nSelBv.push_back(0);
    nSelEv.push_back(0);
    nSelCorrv.push_back(0);
    nSelBCorrv.push_back(0);
    nSelECorrv.push_back(0);
    
    //
    // loop over events
    //    
    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
      genBr->GetEntry(ientry);
      if(charge==-1 && gen->id_1!= EGenType::kMuon) continue;  // check for W-
      if(charge== 1 && gen->id_2!=-EGenType::kMuon) continue;  // check for W+
      infoBr->GetEntry(ientry);     
    
      Double_t weight=1;
      if(doPU>0) weight *= puWeights->GetBinContent(info->nPU+1);
      nEvtsv[ifile]+=weight;
      
      // trigger requirement               
      ULong_t trigger = // MC trigger?  kHLT_Mu15;
      ULong_t trigObj = // MC trigger?  kHLT_Mu15_MuObj;   
      if(!(info->triggerBits & trigger)) continue;  
      
      // good vertex requirement
      if(!(info->hasGoodPV)) continue;
    
      muonArr->Clear();
      muonBr->GetEntry(ientry);
      Int_t nLooseLep=0;
      const mithep::TMuon *goodMuon=0;
      Bool_t passSel=kFALSE;
      for(Int_t i=0; i<muonArr->GetEntriesFast(); i++) {
  	const mithep::TMuon *mu = (mithep::TMuon*)((*muonArr)[i]);

        if(fabs(mu->eta) > ETA_CUT) continue;  // lepton |eta| cut
        if(mu->pt	 < 10)      continue;  // loose lepton pT cut
        if(passMuonLooseID(mu)) nLooseLep++;   // loose lepton selection
        if(nLooseLep>1) {  // extra lepton veto
          passSel=kFALSE;
          break;
        }
        
        if(mu->pt < PT_CUT)		  continue;  // lepton pT cut	
        if(!passMuonID(mu))		  continue;  // lepton selection
        if(!(mu->hltMatchBits & trigObj)) continue;  // check trigger matching
	
	passSel=kTRUE;
	goodMuon=mu;
      }
      
      if(passSel) {
        
	/******** We have a W candidate! HURRAY! ********/
        
	Bool_t isBarrel = (fabs(goodMuon->eta)<ETA_BARREL) ? kTRUE : kFALSE;
        
	Double_t corr=1;
	if(dataHLTEffFile && zmmHLTEffFile) {
	  Double_t effdata = dataHLTEff.getEff(fabs(goodMuon->eta), goodMuon->pt);
	  Double_t effmc   = zmmHLTEff.getEff(fabs(goodMuon->eta), goodMuon->pt);
	  corr *= (1.-(1.-effdata)*(1.-effdata))/(1.-(1.-effmc)*(1.-effmc));
	}
	if(dataSelEffFile && zmmSelEffFile) {
	  Double_t effdata = dataSelEff.getEff(fabs(goodMuon->eta), goodMuon->pt);
	  Double_t effmc   = zmmSelEff.getEff(fabs(goodMuon->eta), goodMuon->pt);
	  corr *= effdata/effmc;
	}
	if(dataTrkEffFile && zmmTrkEffFile) {
	  Double_t effdata = dataTrkEff.getEff(fabs(goodMuon->eta), goodMuon->pt);
	  Double_t effmc   = zmmTrkEff.getEff(fabs(goodMuon->eta), goodMuon->pt);
	  corr *= effdata/effmc;
	}
	if(dataStaEffFile && zmmStaEffFile) {
	  Double_t effdata = dataStaEff.getEff(fabs(goodMuon->eta), goodMuon->pt);
	  Double_t effmc   = zmmStaEff.getEff(fabs(goodMuon->eta), goodMuon->pt);
	  corr *= effdata/effmc;
	}
	
	nSelv[ifile]+=weight;
	nSelCorrv[ifile]+=weight*corr;
  	if(isBarrel) { nSelBv[ifile]+=weight; nSelBCorrv[ifile]+=weight*corr; }
	else	     { nSelEv[ifile]+=weight; nSelECorrv[ifile]+=weight*corr; }
      }
    }
    
    // compute acceptances
    accv.push_back(nSelv[ifile]/nEvtsv[ifile]);   accErrv.push_back(accv[ifile]*sqrt((1.-accv[ifile])/nEvtsv[ifile]));
    accBv.push_back(nSelBv[ifile]/nEvtsv[ifile]); accErrBv.push_back(accBv[ifile]*sqrt((1.-accBv[ifile])/nEvtsv[ifile]));
    accEv.push_back(nSelEv[ifile]/nEvtsv[ifile]); accErrEv.push_back(accEv[ifile]*sqrt((1.-accEv[ifile])/nEvtsv[ifile]));
    
    accCorrv.push_back(nSelCorrv[ifile]/nEvtsv[ifile]);   accErrCorrv.push_back(accCorrv[ifile]*sqrt((1.-accCorrv[ifile])/nEvtsv[ifile]));
    accBCorrv.push_back(nSelBCorrv[ifile]/nEvtsv[ifile]); accErrBCorrv.push_back(accBCorrv[ifile]*sqrt((1.-accBCorrv[ifile])/nEvtsv[ifile]));
    accECorrv.push_back(nSelECorrv[ifile]/nEvtsv[ifile]); accErrECorrv.push_back(accECorrv[ifile]*sqrt((1.-accECorrv[ifile])/nEvtsv[ifile]));
    
    delete infile;
    infile=0, eventTree=0;  
  }
  delete info;
  delete gen;
  delete muonArr;
  
    
  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
   
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  if(charge== 0) cout << " W -> mu nu"  << endl;
  if(charge==-1) cout << " W- -> mu nu" << endl;
  if(charge== 1) cout << " W+ -> mu nu" << endl;
  cout << "  pT > " << PT_CUT << endl;
  cout << "  |eta| < " << ETA_CUT << endl;
  cout << "  Barrel definition: |eta| < " << ETA_BARREL << endl;
  cout << "  Endcap definition: |eta| > " << ETA_ENDCAP << endl;
  cout << endl;
  
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
    cout << "   ================================================" << endl;
    cout << "    Label: " << labelv[ifile] << endl;
    cout << "     File: " << fnamev[ifile] << endl;
    cout << endl;
    cout << "    *** Acceptance ***" << endl;
    cout << "     barrel: " << setw(12) << nSelBv[ifile] << " / " << nEvtsv[ifile] << " = " << accBv[ifile] << " +/- " << accErrBv[ifile];
    cout << "  ==eff corr==> " << accBCorrv[ifile] << " +/- " << accErrBCorrv[ifile] << endl;
    cout << "     endcap: " << setw(12) << nSelEv[ifile] << " / " << nEvtsv[ifile] << " = " << accEv[ifile] << " +/- " << accErrEv[ifile];
    cout << "  ==eff corr==> " << accECorrv[ifile] << " +/- " << accErrECorrv[ifile] << endl;
    cout << "      total: " << setw(12) << nSelv[ifile]  << " / " << nEvtsv[ifile] << " = " << accv[ifile]  << " +/- " << accErrv[ifile];
    cout << "  ==eff corr==> " << accCorrv[ifile]  << " +/- " << accErrCorrv[ifile] << endl;
    cout << endl;
  }
  
  char txtfname[100];
  sprintf(txtfname,"%s/sel.txt",outputDir.Data());
  ofstream txtfile;
  txtfile.open(txtfname);
  txtfile << "*" << endl;
  txtfile << "* SUMMARY" << endl;
  txtfile << "*--------------------------------------------------" << endl;
  if(charge== 0) txtfile << " W -> mu nu"  << endl;
  if(charge==-1) txtfile << " W- -> mu nu" << endl;
  if(charge== 1) txtfile << " W+ -> mu nu" << endl;
  txtfile << "  pT > " << PT_CUT << endl;
  txtfile << "  |eta| < " << ETA_CUT << endl;
  txtfile << "  Barrel definition: |eta| < " << ETA_BARREL << endl;
  txtfile << "  Endcap definition: |eta| > " << ETA_ENDCAP << endl;
  txtfile << endl;
  
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
    txtfile << "   ================================================" << endl;
    txtfile << "    Label: " << labelv[ifile] << endl;
    txtfile << "     File: " << fnamev[ifile] << endl;
    txtfile << endl;
    txtfile << "    *** Acceptance ***" << endl;
    txtfile << "     barrel: " << setw(12) << nSelBv[ifile] << " / " << nEvtsv[ifile] << " = " << accBv[ifile] << " +/- " << accErrBv[ifile];
    txtfile << "  ==eff corr==> " << accBCorrv[ifile] << " +/- " << accErrBCorrv[ifile] << endl;
    txtfile << "     endcap: " << setw(12) << nSelEv[ifile] << " / " << nEvtsv[ifile] << " = " << accEv[ifile] << " +/- " << accErrBv[ifile];
    txtfile << "  ==eff corr==> " << accECorrv[ifile] << " +/- " << accErrECorrv[ifile] << endl;
    txtfile << "      total: " << setw(12) << nSelv[ifile]  << " / " << nEvtsv[ifile] << " = " << accv[ifile]  << " +/- " << accErrv[ifile];
    txtfile << "  ==eff corr==> " << accCorrv[ifile]  << " +/- " << accErrCorrv[ifile] << endl;
    txtfile << endl;
  }
  txtfile.close();  
  
  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;  
      
  gBenchmark->Show("computeAccSelWm"); 
}
