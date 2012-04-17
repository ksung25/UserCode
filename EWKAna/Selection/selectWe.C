//================================================================================================
//
// Select W->enu candidates
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
#include "EWKAna/Ntupler/interface/TElectron.hh"
#include "EWKAna/Ntupler/interface/TVertex.hh"

// lumi section selection with JSON files
#include "MitAna/DataCont/interface/RunLumiRangeMap.h"

// helper functions for lepton ID selection
#include "EWKAna/Utils/LeptonIDCuts.hh"
#endif

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVector;


//=== MAIN MACRO ================================================================================================= 

void selectWe(const TString conf,       // input file
              const TString outputDir   // output directory
) {
  gBenchmark->Start("selectWe");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 

  const Double_t PT_CUT   = 25;
  const Double_t ETA_CUT  = 2.5;
  const Double_t ELE_MASS = 0.000511;
  
  const Double_t ECAL_GAP_LOW  = 1.4442;
  const Double_t ECAL_GAP_HIGH = 1.566;
  
  
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
  Float_t met, metPhi, sumEt, mt;
  Int_t   q;
  LorentzVector *lep=0;
  ///// electron specific /////
  Float_t trkIso, emIso, hadIso;
  Float_t sigieie, hovere, eoverp, fbrem;
  Float_t dphi, deta;
  LorentzVector *sc=0;
  
  // Data structures to store info from TTrees
  mithep::TEventInfo *info  = new mithep::TEventInfo();
  mithep::TGenInfo   *gen   = new mithep::TGenInfo();
  TClonesArray *electronArr = new TClonesArray("mithep::TElectron");
  TClonesArray *pvArr       = new TClonesArray("mithep::TVertex");
  
  TFile *infile=0;
  TTree *eventTree=0;
  
  //
  // loop over samples
  //  
  for(UInt_t isam=0; isam<samplev.size(); isam++) {
    
    // Assume data sample is first sample in .conf file
    // If sample is empty (i.e. contains no ntuple files), skip to next sample
    if(isam==0 && !hasData) continue;
    
    // Assume signal sample is given name "we"
    // If it's the signal sample, toggle flag to store GEN W kinematics
    Bool_t isSignal = (snamev[isam].CompareTo("we",TString::kIgnoreCase)==0);
  
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
    outTree->Branch("scale1fb", &scale1fb, "scale1fb/F");   // event weight per 1/fb (MC)
    outTree->Branch("met",      &met,      "met/F");        // MET
    outTree->Branch("metPhi",   &metPhi,   "metPhi/F");     // phi(MET)
    outTree->Branch("sumEt",    &sumEt,    "sumEt/F");      // Sum ET
    outTree->Branch("mt",       &mt,       "mt/F");         // transverse mass
    outTree->Branch("q",        &q,        "q/I");          // lepton charge
    outTree->Branch("lep", "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &lep);  // lepton 4-vector
    ///// electron specific /////
    outTree->Branch("trkIso",  &trkIso,  "trkIso/F");   // track isolation of tag lepton
    outTree->Branch("emIso",   &emIso,   "emIso/F");    // ECAL isolation of tag lepton
    outTree->Branch("hadIso",  &hadIso,  "hadIso/F");   // HCAL isolation of tag lepton
    outTree->Branch("sigieie", &sigieie, "sigieie/F");  // sigma-ieta-ieta of electron
    outTree->Branch("hovere",  &hovere,  "hovere/F");   // H/E of electron
    outTree->Branch("eoverp",  &eoverp,	 "eoverp/F");   // E/p of electron
    outTree->Branch("fbrem",   &fbrem,   "fbrem/F");    // brem fraction of electron
    outTree->Branch("dphi",    &dphi,	 "dphi/F");     // GSF track - ECAL dphi of electron
    outTree->Branch("deta",    &deta,    "deta/F");     // GSF track - ECAL deta of electron
    outTree->Branch("sc",  "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sc);   // electron Supercluster 4-vector
    
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
      eventTree->SetBranchAddress("Info",     &info);        TBranch *infoBr     = eventTree->GetBranch("Info");
      eventTree->SetBranchAddress("Electron", &electronArr); TBranch *electronBr = eventTree->GetBranch("Electron");
      eventTree->SetBranchAddress("PV",       &pvArr);       TBranch *pvBr       = eventTree->GetBranch("PV");
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
        ULong64_t trigger = (isam==0) ? kHLT_Ele22_CaloIdL_CaloIsoVL        : kHLT_Ele17_CaloIdL_CaloIsoVL;
	ULong64_t trigObj = (isam==0) ? kHLT_Ele22_CaloIdL_CaloIsoVL_EleObj : kHLT_Ele17_CaloIdL_CaloIsoVL_EleObj;  
        if(!(info->triggerBits & trigger)) continue;      
      
        // good vertex requirement
        if(!(info->hasGoodPV)) continue;
        pvArr->Clear();
        pvBr->GetEntry(ientry);
      
        //
	// SELECTION PROCEDURE:
	//  (1) Look for 1 good electron matched to trigger
	//  (2) Reject event if another electron is present passing looser cuts
	//
	electronArr->Clear();
        electronBr->GetEntry(ientry);
	Int_t nLooseLep=0;
	const mithep::TElectron *goodEle=0;
	Bool_t passSel=kFALSE;
        for(Int_t i=0; i<electronArr->GetEntriesFast(); i++) {
          const mithep::TElectron *ele = (mithep::TElectron*)((*electronArr)[i]);
	  
	  // check ECAL gap
	  if(fabs(ele->scEta)>=ECAL_GAP_LOW && fabs(ele->scEta)<=ECAL_GAP_HIGH) continue;
	  
	  if(passEleLooseID(ele)) nLooseLep++;
	  if(nLooseLep>1) {  // extra lepton veto
	    passSel=kFALSE;
	    break;
	  }
	  
	  if(ele->scEt        < PT_CUT)      continue;  // lepton pT cut
	  if(fabs(ele->scEta) > ETA_CUT)     continue;  // lepton |eta| cut
	  if(!passEleID(ele))                continue;  // lepton selection
	  if(!(ele->hltMatchBits & trigObj)) continue;  // check trigger matching
	  
	  passSel=kTRUE;
	  goodEle = ele;  
	}
	
	if(passSel) {
	  
	  /******** We have a W candidate! HURRAY! ********/
	    
	  nsel+=weight;
          nselvar+=weight*weight;
	  
	  LorentzVector vLep(goodEle->pt, goodEle->eta, goodEle->phi, ELE_MASS);  
	  LorentzVector vSC(goodEle->scEt, goodEle->scEta, goodEle->scPhi, ELE_MASS); 	  
	  
	  //
	  // Fill tree
	  //
	  runNum   = info->runNum;
	  lumiSec  = info->lumiSec;
	  evtNum   = info->evtNum;
	  npv	   = pvArr->GetEntriesFast();
	  npu	   = info->nPU;
	  genWPt   = (isSignal) ? gen->vpt : 0;
	  genWPhi  = (isSignal) ? gen->vphi : 0;
	  scale1fb = weight;
	  met	   = info->pfMET;
	  metPhi   = info->pfMETphi;
	  sumEt    = info->pfSumET;
	  mt       = sqrt( 2.0 * (vLep.Pt()) * (info->pfMET) * (1.0-cos(toolbox::deltaPhi(vLep.Phi(),info->pfMETphi))) );
	  q        = goodEle->q;	  
	  lep      = &vLep;	  
	  
	  ///// electron specific /////
	  sc	  = &vSC;
	  trkIso  = goodEle->trkIso03;
	  emIso   = goodEle->emIso03;
	  hadIso  = goodEle->hadIso03;
	  sigieie = goodEle->sigiEtaiEta;
	  hovere  = goodEle->HoverE;
	  eoverp  = goodEle->EoverP;
	  fbrem   = goodEle->fBrem;
	  dphi    = goodEle->deltaPhiIn;
	  deta    = goodEle->deltaEtaIn;
	  
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
  delete electronArr;
  delete pvArr;
  
    
  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
   
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  cout << " W -> e nu" << endl;
  cout << "  pT > " << PT_CUT << endl;
  cout << "  |eta| < " << ETA_CUT << endl;
  cout << endl;

  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;  
      
  gBenchmark->Show("selectWe"); 
}
