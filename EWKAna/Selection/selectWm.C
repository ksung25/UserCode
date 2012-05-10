//================================================================================================
//
// Select W->munu candidates
//
//  * outputs ROOT files of events passing selection
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TClonesArray.h>           // ROOT array class
#include <TBenchmark.h>             // class to track macro running statistics
#include <TVector2.h>               // 2D vector class
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include "Math/LorentzVector.h"     // 4-vector class

#include "ConfParse.hh"             // input conf file parser
#include "EWKAna/Utils/CSample.hh"  // helper class to handle samples
#include "EWKAna/Utils/MyTools.hh"  // various helper functions

// define structures to read in ntuple
#include "EWKAna/Ntupler/interface/EWKAnaDefs.hh"
#include "EWKAna/Ntupler/interface/TEventInfo.hh"
#include "EWKAna/Ntupler/interface/TGenInfo.hh"
#include "EWKAna/Ntupler/interface/TMuon.hh"
#include "EWKAna/Ntupler/interface/TVertex.hh"

// lumi section selection with JSON files
#include "MitAna/DataCont/interface/RunLumiRangeMap.h"

// helper functions for lepton ID selection
#include "EWKAna/Utils/LeptonIDCuts.hh"
#endif

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVector;


//=== MAIN MACRO ================================================================================================= 

void selectWm(const TString conf,      // input file
              const TString outputDir  // output directory
) {
  gBenchmark->Start("selectWm");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 

  const Double_t PT_CUT    = 20;
  const Double_t ETA_CUT   = 2.4;
  const Double_t MUON_MASS = 0.105658369;


  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  

  vector<TString>  snamev;      // sample name (for output files)  
  vector<CSample*> samplev;     // data/MC samples

  //
  // parse .conf file
  //
  confParse(conf, snamev, samplev);
  const Bool_t hasData = (samplev[0]->fnamev.size()>0);

  // Create output directory
  gSystem->mkdir(outputDir,kTRUE);
  const TString ntupDir = outputDir + TString("/ntuples");
  gSystem->mkdir(ntupDir,kTRUE);
  
  //
  // Declare output ntuple variables
  //
  UInt_t  runNum, lumiSec, evtNum;
  UInt_t  npv, npu;
  Float_t genWPt, genWPhi;
  Float_t scale1fb;
  Float_t met, metPhi, sumEt, mt, u1, u2;
  Int_t   q;
  LorentzVector *lep=0;
  ///// muon specific /////
  Float_t trkIso, emIso, hadIso;
  UInt_t nPixHits, nTkHits, nValidHits, nMatch;  
  
  // Data structures to store info from TTrees
  mithep::TEventInfo *info = new mithep::TEventInfo();
  mithep::TGenInfo   *gen  = new mithep::TGenInfo();
  TClonesArray *muonArr    = new TClonesArray("mithep::TMuon");
  TClonesArray *pvArr      = new TClonesArray("mithep::TVertex");
  
  TFile *infile=0;
  TTree *eventTree=0;
  
  //
  // loop over samples
  //  
  for(UInt_t isam=0; isam<samplev.size(); isam++) {
    
    // Assume data sample is first sample in .conf file
    // If sample is empty (i.e. contains no ntuple files), skip to next sample
    if(isam==0 && !hasData) continue;
    
    // Assume signal sample is given name "wm"
    // If it's the signal sample, toggle flag to do GEN matching
    Bool_t isSignal = (snamev[isam].CompareTo("wm",TString::kIgnoreCase)==0);
  
    CSample* samp = samplev[isam];
  
    //
    // Set up output ntuple
    //
    TString outfilename = ntupDir + TString("/") + snamev[isam] + TString("_select.root");
    TFile *outFile = new TFile(outfilename,"RECREATE"); 
    TTree *outTree = new TTree("Events","Events");

    outTree->Branch("runNum",   &runNum,   "runNum/i");     // event run number
    outTree->Branch("lumiSec",  &lumiSec,  "lumiSec/i");    // event lumi section
    outTree->Branch("evtNum",   &evtNum,   "evtNum/i");     // event number
    outTree->Branch("npv",      &npv,      "npv/i");        // number of primary vertices
    outTree->Branch("npu",      &npu,      "npu/i");        // number of in-time PU events (MC)
    outTree->Branch("genWPt",   &genWPt,   "genWPt/F");     // GEN W boson pT (signal MC)
    outTree->Branch("genWPhi",  &genWPhi,  "genWPhi/F");    // GEN W boson phi (signal MC)
    outTree->Branch("scale1fb", &scale1fb, "scale1fb/F");   // event weight per 1/fb (MC)
    outTree->Branch("met",      &met,      "met/F");        // MET
    outTree->Branch("metPhi",   &metPhi,   "metPhi/F");     // phi(MET)
    outTree->Branch("sumEt",    &sumEt,    "sumEt/F");      // Sum ET
    outTree->Branch("mt",       &mt,       "mt/F");         // transverse mass
    outTree->Branch("u1",       &u1,       "u1/F");         // parallel component of recoil
    outTree->Branch("u2",       &u2,       "u2/F");         // perpendicular component of recoil
    outTree->Branch("q",        &q,        "q/I");          // lepton charge
    outTree->Branch("lep", "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &lep);   // lepton 4-vector
    ///// muon specific /////
    outTree->Branch("trkIso",     &trkIso,     "trkIso/F");       // track isolation of lepton
    outTree->Branch("emIso",      &emIso,      "emIso/F");        // ECAL isolation of lepton
    outTree->Branch("hadIso",     &hadIso,     "hadIso/F");       // HCAL isolation of lepton
    outTree->Branch("nPixHits",   &nPixHits,   "nPixHits/i");	  // number of pixel hits of muon
    outTree->Branch("nTkHits",    &nTkHits,    "nTkHits/i");	  // number of tracker hits of muon
    outTree->Branch("nMatch",     &nMatch,     "nMatch/i");	  // number of matched segments of muon	 
    outTree->Branch("nValidHits", &nValidHits, "nValidHits/i");   // number of valid muon hits of muon 
    
    //
    // loop through files
    //
    const UInt_t nfiles = samp->fnamev.size();
    for(UInt_t ifile=0; ifile<nfiles; ifile++) {  

      // Read input file and get the TTrees
      cout << "Processing " << samp->fnamev[ifile] << " [xsec = " << samp->xsecv[ifile] << " pb] ... "; cout.flush();
      infile = new TFile(samp->fnamev[ifile]); 
      assert(infile);
      
      Bool_t hasJSON = kFALSE;
      mithep::RunLumiRangeMap rlrm;
      if(samp->jsonv[ifile].CompareTo("NONE")!=0) { 
        hasJSON = kTRUE;
        rlrm.AddJSONFile(samp->jsonv[ifile].Data()); 
      }
  
      eventTree = (TTree*)infile->Get("Events");
      assert(eventTree);  
      eventTree->SetBranchAddress("Info", &info);    TBranch *infoBr = eventTree->GetBranch("Info");
      eventTree->SetBranchAddress("Muon", &muonArr); TBranch *muonBr = eventTree->GetBranch("Muon");
      eventTree->SetBranchAddress("PV",   &pvArr);   TBranch *pvBr   = eventTree->GetBranch("PV");
      TBranch *genBr=0;
      if(isSignal) {
        eventTree->SetBranchAddress("Gen", &gen);
	genBr = eventTree->GetBranch("Gen");
      }
    
      // Compute MC event weight per 1/fb
      Double_t weight = 1;
      const Double_t xsec = samp->xsecv[ifile];
      if(xsec>0) weight = 1000.*xsec/(Double_t)eventTree->GetEntries();     

      //
      // loop over events
      //
      Double_t nsel=0, nselvar=0;
      for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
        infoBr->GetEntry(ientry);
	
	if(genBr) genBr->GetEntry(ientry);
        
	// check for certified lumi (if applicable)
        mithep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);      
        if(hasJSON && !rlrm.HasRunLumi(rl)) continue;  

        // trigger requirement               
        ULong64_t trigger = kHLT_Mu15_eta2p1;
	ULong64_t trigObj = kHLT_Mu15_eta2p1_MuObj;   
        if(!(info->triggerBits & trigger)) continue;      
      
        // good vertex requirement
        if(!(info->hasGoodPV)) continue;
        pvArr->Clear();
        pvBr->GetEntry(ientry);
      
        //
	// SELECTION PROCEDURE:
	//  (1) Look for 1 good muon matched to trigger
	//  (2) Reject event if another muon is present passing looser cuts
	//
	muonArr->Clear();
        muonBr->GetEntry(ientry);
	Int_t nLooseLep=0;
	const mithep::TMuon *goodMuon=0;
	Bool_t passSel=kFALSE;
        for(Int_t i=0; i<muonArr->GetEntriesFast(); i++) {
          const mithep::TMuon *mu = (mithep::TMuon*)((*muonArr)[i]);
	
	  if(passMuonLooseID(mu)) nLooseLep++;
	  if(nLooseLep>1) {  // extra lepton veto
	    passSel=kFALSE;
	    break;
	  }
	  
	  if(mu->pt        < PT_CUT)        continue;  // lepton pT cut
	  if(fabs(mu->eta) > ETA_CUT)       continue;  // lepton |eta| cut
	  if(!passMuonID(mu))               continue;  // lepton selection
	  if(!(mu->hltMatchBits & trigObj)) continue;  // check trigger matching
	  
	  passSel=kTRUE;
	  goodMuon = mu;
	}
	
	if(passSel) {
	  
	  /******** We have a W candidate! HURRAY! ********/
	    
	  nsel+=weight;
          nselvar+=weight*weight;
	  
	  LorentzVector vLep(goodMuon->pt, goodMuon->eta, goodMuon->phi, MUON_MASS);  
	  	  
	  //
	  // Fill tree
	  //
	  runNum   = info->runNum;
	  lumiSec  = info->lumiSec;
	  evtNum   = info->evtNum;
	  npv	   = pvArr->GetEntriesFast();
	  npu	   = info->nPU;
	  genWPt   = 0;
	  genWPhi  = 0;
	  u1       = 0;
	  u2       = 0;
	  if(isSignal) {
	    genWPt   = gen->vpt;
            genWPhi  = gen->vphi;
	    TVector2 vWPt((gen->vpt)*cos(gen->vphi),(gen->vpt)*sin(gen->vphi));
	    TVector2 vLepPt(vLep.Px(),vLep.Py());      
            TVector2 vMet((info->pfMET)*cos(info->pfMETphi), (info->pfMET)*sin(info->pfMETphi));        
            TVector2 vU = vMet+vLepPt;
            u1 = -((vWPt.Px())*(vU.Px()) + (vWPt.Py())*(vU.Py()))/(gen->vpt);  // u1 = -(pT . u)/|pT|
            u2 =  ((vWPt.Px())*(vU.Py()) - (vWPt.Py())*(vU.Px()))/(gen->vpt);  // u2 =  (pT x u)/|pT|
	  }
	  scale1fb = weight;
	  met	   = info->pfMET;
	  metPhi   = info->pfMETphi;
	  sumEt    = info->pfSumET;
	  mt       = sqrt( 2.0 * (vLep.Pt()) * (info->pfMET) * (1.0-cos(toolbox::deltaPhi(vLep.Phi(),info->pfMETphi))) );
	  q        = goodMuon->q;	  
	  lep      = &vLep;	  

	  ///// muon specific /////
	  trkIso     = goodMuon->trkIso03;
	  emIso      = goodMuon->emIso03;
	  hadIso     = goodMuon->hadIso03;
	  nPixHits   = goodMuon->nPixHits;
	  nTkHits    = goodMuon->nTkHits;
	  nMatch     = goodMuon->nMatch;
	  nValidHits = goodMuon->nValidHits;
	  
	  outTree->Fill();
        }
      }
      delete infile;
      infile=0, eventTree=0;    

      cout << nsel  << " +/- " << sqrt(nselvar) << endl;
    }
    outFile->Write();
    outFile->Close();
  }
  delete info;
  delete gen;
  delete muonArr;
  delete pvArr;
  
    
  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
   
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  cout << " W -> mu nu" << endl;
  cout << "  pT > " << PT_CUT << endl;
  cout << "  |eta| < " << ETA_CUT << endl;
  cout << endl;
  
  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;  
      
  gBenchmark->Show("selectWm"); 
}
