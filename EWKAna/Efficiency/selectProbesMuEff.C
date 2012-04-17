//================================================================================================
//
// Select probes for muon efficiencies with Tag&Probe method
//
//  * outputs ROOT file with a TTree of probes
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                        // access to gROOT, entry point to ROOT system
#include <TSystem.h>                      // interface to OS
#include <TFile.h>                        // file handle class
#include <TTree.h>                        // class to access ntuples
#include <TBenchmark.h>                   // class to track macro running statistics
#include <vector>                         // STL vector class
#include <iostream>                       // standard I/O
#include <iomanip>                        // functions to format standard I/O
#include "Math/LorentzVector.h"           // 4-vector class

// structure for output ntuple
#include "EffData.hh" 
#endif

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVector;


//=== MAIN MACRO ================================================================================================= 

void selectProbesMuEff(const TString infilename,          // input ntuple
                       const TString outputDir, 	  // output directory
		       const Int_t   effType, 	          // type of efficiency to compute
		       const Bool_t  doGenMatch = kFALSE  // match to generator leptons
) {
  gBenchmark->Start("selectProbesMuEff");

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
  
  enum { eHLTEff, eSelEff, eTrkEff, eStaEff };  // event category enum
  if(effType > eStaEff) {
    cout << "Invalid effType option! Exiting..." << endl;
    return;
  }
  
  enum { eMuMu2HLT=1, eMuMu1HLT, eMuMuNoSel, eMuSta, eMuTrk };
  
  Double_t nProbes = 0;
  
  //
  // Set up output ntuple
  //
  gSystem->mkdir(outputDir,kTRUE);
  TFile *outFile = new TFile(outputDir+TString("/probes.root"),"RECREATE"); 
  TTree *outTree = new TTree("Events","Events");
  EffData data;
  outTree->Branch("Events",&data.mass,"mass/F:pt:eta:phi:weight:q/I:npv/i:npu:pass:runNum:lumiSec:evtNum");

  //
  // Declare output ntuple variables
  //
  UInt_t  runNum, lumiSec, evtNum;
  UInt_t  matchGen;
  UInt_t  category;
  UInt_t  npv, npu;
  Float_t scale1fb;
  Float_t met, metPhi, sumEt, u1, u2;
  Int_t   q1, q2;
  LorentzVector *dilep=0, *lep1=0, *lep2=0;
  LorentzVector *sta1=0, *sta2=0;
  
  // Read input file and get the TTrees
  cout << "Processing " << infilename << "..." << endl;
  TFile *infile = new TFile(infilename);	 assert(infile);
  TTree *intree = (TTree*)infile->Get("Events"); assert(intree);

  intree->SetBranchAddress("runNum",   &runNum);     // event run number
  intree->SetBranchAddress("lumiSec",  &lumiSec);    // event lumi section
  intree->SetBranchAddress("evtNum",   &evtNum);     // event number
  intree->SetBranchAddress("matchGen", &matchGen);   // event has both leptons matched to MC Z->ll
  intree->SetBranchAddress("category", &category);   // dilepton category
  intree->SetBranchAddress("npv",      &npv);	     // number of primary vertices
  intree->SetBranchAddress("npu",      &npu);	     // number of in-time PU events (MC)
  intree->SetBranchAddress("scale1fb", &scale1fb);   // event weight per 1/fb (MC)
  intree->SetBranchAddress("met",      &met);	     // MET
  intree->SetBranchAddress("metPhi",   &metPhi);     // phi(MET)
  intree->SetBranchAddress("sumEt",    &sumEt);	     // Sum ET
  intree->SetBranchAddress("u1",       &u1);	     // parallel component of recoil
  intree->SetBranchAddress("u2",       &u2);	     // perpendicular component of recoil
  intree->SetBranchAddress("q1",       &q1);	     // charge of tag lepton
  intree->SetBranchAddress("q2",       &q2);	     // charge of probe lepton
  intree->SetBranchAddress("dilep",    &dilep);      // dilepton 4-vector
  intree->SetBranchAddress("lep1",     &lep1);       // tag lepton 4-vector
  intree->SetBranchAddress("lep2",     &lep2);       // probe lepton 4-vector
  intree->SetBranchAddress("sta1",     &sta1);       // tag STA muon 4-vector
  intree->SetBranchAddress("sta2",     &sta2);       // probe STA muon 4-vector 
  
  //
  // loop over events
  //
  for(UInt_t ientry=0; ientry<intree->GetEntries(); ientry++) {
    intree->GetEntry(ientry);
    
    // check GEN match if necessary
    if(doGenMatch && !matchGen) continue;
    
    Bool_t  pass = kFALSE;
    Float_t mass = 0;
    
    if(effType==eHLTEff) {
      //
      // probe = muon passing selection
      // pass  = matched to HLT
      // * MuMu2HLT event means a passing probe, MuMu1HLT event means a failing probe
      //   all other categories do not satisfy probe requirements
      //    
      if     (category==eMuMu2HLT)  { pass=kTRUE; }
      else if(category==eMuMu1HLT)  { pass=kFALSE; }
      else if(category==eMuMuNoSel) { continue; }
      else if(category==eMuSta)     { continue; }
      else                          { continue; }
      
      mass = dilep->M();
    
    } else if(effType==eSelEff) {
      //
      // probe = GLB muon
      // pass  = passing selection
      // * MuMu2HLT, MuMu1HLT event means a passing probe, MuMuNoSel event means a failing probe,
      //   all other categories do not satisfy probe requirements
      //    
      if     (category==eMuMu2HLT)  { pass=kTRUE; }
      else if(category==eMuMu1HLT)  { pass=kTRUE; }
      else if(category==eMuMuNoSel) { pass=kFALSE; }
      else if(category==eMuSta)     { continue; }
      else                          { continue; }
      
      mass = dilep->M();
    
    } else if(effType==eTrkEff) {
      //
      // probe = STA muon
      // pass  = is also a GLB muon
      // * MuMu2HLT, MuMu1HLT, MuMuNoSel event means a passing probe, MuSta event means a failing probe, 
      //   MuTrk event does not satisfy probe requirements
      //    
      if     (category==eMuMu2HLT)  { pass=kTRUE; }
      else if(category==eMuMu1HLT)  { pass=kTRUE; }
      else if(category==eMuMuNoSel) { pass=kTRUE; }
      else if(category==eMuSta)     { pass=kFALSE; }
      else                          { continue; }
      
      // compute mass using probe STA muon pT
      LorentzVector tp = *lep1 + *sta2;
      mass = tp.M();    
    
    } else if(effType==eStaEff) {
      //
      // probe = tracker track
      // pass  = is also a GLB muon
      // * MuMu2HLT, MuMu1HLT, MuMuNoSel event means a passing probe, MuTrk event means a failing probe, 
      //   MuSta event does not satisfy probe requirements
      //    
      if     (category==eMuMu2HLT)  { pass=kTRUE; }
      else if(category==eMuMu1HLT)  { pass=kTRUE; }
      else if(category==eMuMuNoSel) { pass=kTRUE; }
      else if(category==eMuSta)     { continue; }
      else                          { pass=kFALSE; }
      
      mass = dilep->M();
    }
    
    nProbes += scale1fb;

    // Fill tree
    data.mass	 = mass;
    data.pt	 = (effType==eTrkEff) ? sta2->Pt() : lep2->Pt();
    data.eta	 = (effType==eTrkEff) ? sta2->Eta() : lep2->Eta();
    data.phi	 = (effType==eTrkEff) ? sta2->Phi() : lep2->Phi();
    data.weight  = scale1fb;
    data.q	 = q2;
    data.npv	 = npv;
    data.npu	 = npu;
    data.pass	 = (pass) ? 1 : 0;
    data.runNum  = runNum;
    data.lumiSec = lumiSec;
    data.evtNum  = evtNum;
    outTree->Fill(); 
  }  
  delete infile;
  infile=0, intree=0;	   


  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
   
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  cout << endl;

  cout << " Number of probes selected: " << nProbes << endl;
  
  outFile->Write();
  outFile->Close();
  delete outFile;
  
  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;  
      
  gBenchmark->Show("selectProbesMuEff"); 
}
