//================================================================================================
//
// Select Z->mumu candidates
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

void selectZmm(const TString conf,      // input file
               const TString outputDir  // output directory
) {
  gBenchmark->Start("selectZmm");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 

  const Double_t MASS_LOW  = 40;
  const Double_t MASS_HIGH = 200;
  const Double_t PT_CUT    = 20;
  const Double_t ETA_CUT   = 2.4;
  const Double_t MUON_MASS = 0.105658369;


  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  

  enum { eMuMu2HLT=1, eMuMu1HLT, eMuMuNoSel, eMuSta, eMuTrk };  // event category enum
  
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
  UInt_t  matchGen;
  UInt_t  category;
  UInt_t  npv, npu;
  Float_t genZPt, genZPhi;
  Float_t scale1fb;
  Float_t met, metPhi, sumEt, u1, u2;
  Int_t   q1, q2;
  LorentzVector *dilep=0, *lep1=0, *lep2=0;
  ///// muon specific /////
  Float_t trkIso1, emIso1, hadIso1, trkIso2, emIso2, hadIso2;
  UInt_t nPixHits1, nTkHits1, nPixHits2, nTkHits2;
  UInt_t nValidHits1, nMatch1, nValidHits2, nMatch2;
  LorentzVector *sta1=0, *sta2=0;
  
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
    
    // Assume signal sample is given name "zmm"
    // If it's the signal sample, toggle flag to do GEN matching
    Bool_t isSignal = (snamev[isam].CompareTo("zmm",TString::kIgnoreCase)==0);
    
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
    outTree->Branch("matchGen", &matchGen, "matchGen/i");   // event has both leptons matched to MC Z->ll
    outTree->Branch("category", &category, "category/i");   // dilepton category
    outTree->Branch("npv",      &npv,      "npv/i");        // number of primary vertices
    outTree->Branch("npu",      &npu,      "npu/i");        // number of in-time PU events (MC)
    outTree->Branch("genZPt",   &genZPt,   "genZPt/F");     // GEN Z boson pT (signal MC)
    outTree->Branch("genZPhi",  &genZPhi,  "genZPhi/F");    // GEN Z boson phi (signal MC)
    outTree->Branch("scale1fb", &scale1fb, "scale1fb/F");   // event weight per 1/fb (MC)
    outTree->Branch("met",      &met,      "met/F");        // MET
    outTree->Branch("metPhi",   &metPhi,   "metPhi/F");     // phi(MET)
    outTree->Branch("sumEt",    &sumEt,    "sumEt/F");      // Sum ET
    outTree->Branch("u1",       &u1,       "u1/F");         // parallel component of recoil
    outTree->Branch("u2",       &u2,       "u2/F");         // perpendicular component of recoil
    outTree->Branch("q1",       &q1,       "q1/I");         // charge of tag lepton
    outTree->Branch("q2",       &q2,       "q2/I");         // charge of probe lepton
    outTree->Branch("dilep", "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &dilep);  // dilepton 4-vector
    outTree->Branch("lep1",  "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &lep1);   // tag lepton 4-vector
    outTree->Branch("lep2",  "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &lep2);   // probe lepton 4-vector    
    ///// muon specific /////
    outTree->Branch("trkIso1",     &trkIso1,     "trkIso1/F");       // track isolation of tag lepton
    outTree->Branch("trkIso2",     &trkIso2,     "trkIso2/F");       // track isolation of probe lepton
    outTree->Branch("emIso1",      &emIso1,      "emIso1/F");        // ECAL isolation of tag lepton
    outTree->Branch("emIso2",      &emIso2,      "emIso2/F");        // ECAL isolation of probe lepton
    outTree->Branch("hadIso1",     &hadIso1,     "hadIso1/F");       // HCAL isolation of tag lepton
    outTree->Branch("hadIso2",     &hadIso2,     "hadIso2/F");       // HCAL isolation of probe lepton
    outTree->Branch("nPixHits1",   &nPixHits1,	 "nPixHits1/i");     // number of pixel hits of tag muon
    outTree->Branch("nPixHits2",   &nPixHits2,	 "nPixHits2/i");     // number of pixel hits of probe muon
    outTree->Branch("nTkHits1",    &nTkHits1,	 "nTkHits1/i");      // number of tracker hits of tag muon
    outTree->Branch("nTkHits2",    &nTkHits2,	 "nTkHits2/i");      // number of tracker hits of probe muon
    outTree->Branch("nMatch1",     &nMatch1,	 "nMatch1/i");       // number of matched segments of tag muon
    outTree->Branch("nMatch2",     &nMatch2,	 "nMatch2/i");       // number of matched segments of probe muon    
    outTree->Branch("nValidHits1", &nValidHits1, "nValidHits1/i");   // number of valid muon hits of tag muon
    outTree->Branch("nValidHits2", &nValidHits2, "nValidHits2/i");   // number of valid muon hits of probe muon
    outTree->Branch("sta1",  "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sta1);   // tag STA muon 4-vector
    outTree->Branch("sta2",  "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sta2);   // probe STA muon 4-vector 
    
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
  
      eventTree = (TTree*)infile->Get("Events"); assert(eventTree);  
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
	//  (1) Find a good muon matched to trigger -> this will be the "tag"
	//  (2) Pair the tag with various probe types which form a tag+probe mass inside 
	//      the Z window and divide candidates into exclusive categories as follows:
	//      (a) if probe is a good muon matched to trigger                  -> MuMu2HLT category
	//      (b) if probe is a good muon not matched to trigger              -> MuMu1HLT category
	//      (c) if probe is a muon failing selection cuts                   -> MuMuNoSel category
	//      (d) if probe is a standalone muon but not global                -> MuSta category
	//      (e) if probe is a tracker muon or non-muon track but not global -> MuTrk category
	//
	muonArr->Clear();
        muonBr->GetEntry(ientry);
        for(Int_t i1=0; i1<muonArr->GetEntriesFast(); i1++) {
          const mithep::TMuon *tag = (mithep::TMuon*)((*muonArr)[i1]);
	
	  if(tag->pt        < PT_CUT)        continue;  // lepton pT cut
	  if(fabs(tag->eta) > ETA_CUT)       continue;  // lepton |eta| cut
	  if(!passMuonID(tag))               continue;  // lepton selection
	  if(!(tag->hltMatchBits & trigObj)) continue;  // check trigger matching
	
	  LorentzVector vTag(tag->pt, tag->eta, tag->phi, MUON_MASS);
	  LorentzVector vTagSta(tag->staPt, tag->staEta, tag->staPhi, MUON_MASS);
	
	  for(Int_t i2=0; i2<muonArr->GetEntriesFast(); i2++) {
	    if(i1==i2) continue;
	    const mithep::TMuon *probe = (mithep::TMuon*)((*muonArr)[i2]);
	  
	    if(tag->q == probe->q)         continue;  // opposite charge requirement
	    if(probe->pt        < PT_CUT)  continue;  // lepton pT cut
	    if(fabs(probe->eta) > ETA_CUT) continue;  // lepton |eta| cut
	    
	    LorentzVector vProbe(probe->pt, probe->eta, probe->phi, MUON_MASS);
	    LorentzVector vProbeSta(0,0,0,0);
	    if(probe->typeBits & kStandalone)
	      vProbeSta.SetCoordinates(probe->staPt, probe->staEta, probe->staPhi, MUON_MASS);

	    // mass window
	    LorentzVector vDilep = vTag + vProbe;
	    if((vDilep.M()<MASS_LOW) || (vDilep.M()>MASS_HIGH)) continue;
	    
	    // determine event category
	    UInt_t icat=0;
	    if(passMuonID(probe)) {
	      if(probe->hltMatchBits & trigObj) {
	        if(i1>i2) continue;  // make sure we don't double count MuMu2HLT category
		icat=eMuMu2HLT;
		
	      } else {
	        icat=eMuMu1HLT;
	      }
	    } 
	    else if(probe->typeBits & kGlobal)    { icat=eMuMuNoSel; } 
	    else if(probe->typeBits==kStandalone) { icat=eMuSta; } 
	    else                                  { icat=eMuTrk; }
	    if(icat==0) continue;
	    
	    
	    /******** We have a Z candidate! HURRAY! ********/
	    
	    nsel+=weight;
            nselvar+=weight*weight;
	    
	    // Perform matching of dileptons to GEN leptons from Z decay
	    Bool_t hasGenMatch = kFALSE;
	    if(isSignal) {
	      Bool_t match1 = ( (abs(gen->id_1)==EGenType::kMuon) && ((toolbox::deltaR(tag->eta, tag->phi, gen->eta_1, gen->phi_1) < 0.5)) )
	                      || ( (abs(gen->id_2)==EGenType::kMuon) && ((toolbox::deltaR(tag->eta, tag->phi, gen->eta_2, gen->phi_2) < 0.5)) );
	      Bool_t match2 = ( (abs(gen->id_1)==EGenType::kMuon) && ((toolbox::deltaR(probe->eta, probe->phi, gen->eta_1, gen->phi_1) < 0.5)) )
	                      || ( (abs(gen->id_2)==EGenType::kMuon) && ((toolbox::deltaR(probe->eta, probe->phi, gen->eta_2, gen->phi_2) < 0.5)) );
	      if(match1 && match2) hasGenMatch = kTRUE;
	    };
	    
	    //
	    // Fill tree
	    //
	    runNum   = info->runNum;
	    lumiSec  = info->lumiSec;
	    evtNum   = info->evtNum;
	    matchGen = hasGenMatch ? 1 : 0;
	    category = icat;
	    npv      = pvArr->GetEntriesFast();
	    npu      = info->nPU;
	    genZPt   = (isSignal) ? gen->vpt : 0;
	    genZPhi  = (isSignal) ? gen->vphi : 0;
	    scale1fb = weight;
	    met      = info->pfMET;
	    metPhi   = info->pfMETphi;
	    sumEt    = info->pfSumET;
	    
	    lep1 = &vTag;
	    q1   = tag->q;
	    
	    lep2 = &vProbe;
	    q2   = probe->q;
	    
	    dilep = &vDilep;
	    
	    TVector2 vZPt((vDilep.Pt())*cos(vDilep.Phi()),(vDilep.Pt())*sin(vDilep.Phi()));        
            TVector2 vMet((info->pfMET)*cos(info->pfMETphi), (info->pfMET)*sin(info->pfMETphi));        
            TVector2 vU = vMet+vZPt;
            u1 = -((vDilep.Px())*(vU.Px()) + (vDilep.Py())*(vU.Py()))/(vDilep.Pt());  // u1 = -(pT . u)/|pT|
            u2 =  ((vDilep.Px())*(vU.Py()) - (vDilep.Py())*(vU.Px()))/(vDilep.Pt());  // u2 =  (pT x u)/|pT|
	  
	    ///// muon specific /////
	    sta1        = &vTagSta;
	    trkIso1     = tag->trkIso03;
	    emIso1      = tag->emIso03;
	    hadIso1     = tag->hadIso03;
	    nPixHits1   = tag->nPixHits;
	    nTkHits1    = tag->nTkHits;
	    nMatch1     = tag->nMatch;
	    nValidHits1 = tag->nValidHits;
	    
	    sta2        = &vProbeSta;
	    trkIso2     = probe->trkIso03;
	    emIso2      = probe->emIso03;
	    hadIso2     = probe->hadIso03;
	    nPixHits2   = probe->nPixHits;
	    nTkHits2    = probe->nTkHits;
	    nMatch2     = probe->nMatch;
	    nValidHits2 = probe->nValidHits;
	    
	    outTree->Fill();
	  }
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
  cout << " Z -> mu mu" << endl;
  cout << "  Mass window: [" << MASS_LOW << ", " << MASS_HIGH << "]" << endl;
  cout << "  pT > " << PT_CUT << endl;
  cout << "  |eta| < " << ETA_CUT << endl;
  cout << endl;
  
  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;  
      
  gBenchmark->Show("selectZmm"); 
}
