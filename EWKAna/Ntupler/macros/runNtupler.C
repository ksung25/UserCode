#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/Catalog/interface/Catalog.h"
#include "MitAna/TreeMod/interface/Analysis.h"
#include "EWKAna/Ntupler/interface/NtuplerMod.hh"
#endif

using namespace mithep;

/*
 *  useGen options:
 *  0: ignore GEN info
 *  1: H->WW/ZZ->2l2nu
 *  2: Z->ll
 *  3: W->lnu
 *  4: WW->2l2nu
 *  5: VZ->X+ll
 *  6: W+jets (flag for Wgamma)
 *  7: Z+jets (flag for Wgamma)
 */
   
//==================================================================================================
/*
 * Run on a BAMBU fileset
 *
 * Example usage:
 *   root -l -q -b runNtupler.C+\(\"0000\",\"f11-h250zz2l-gf-v14b-pu\",\"cern/filefi/025\",\"/home/mitprod/catalog\",0,1,0,-1,0\)
 *
 * Output file name has standard format: <dataset>_<fileset>_ntuple.root
 *
 */
void runNtupler(
    const char   *fileset,      // "4-digit" string that labels a group of files
    const char   *dataset,      // BAMBU dataset name
    const char   *book,         // BAMBU book containing the dataset
    const char   *catalogDir,   // BAMBU catalog directory
    const Bool_t  isData,       // flag to indicate processing of collision data
    const Int_t   useGen,       // which MC process?
    const Int_t   fsrmode,      // indicates how generator level FSR effects are stored
    const Int_t   nevents,      // number of events to process
    const Bool_t  skipHLTFail   // skip events if no HLT accept
  )
{
  gDebugMask  = Debug::kAnalysis;  // debug message category
  gDebugLevel = 1;                 // higher level allows more messages to print
 
  char output[100];
  sprintf(output,"%s_%s_ntuple.root",dataset,fileset); 
  
  // muon kinematics
  const Double_t muPtMin  = 0;
  const Double_t muPtMax  = 7000;
  const Double_t muEtaMin = -3;
  const Double_t muEtaMax =  3;

  // electron kinematics
  const Double_t eleEtMin  = 9;
  const Double_t eleEtMax  = 7000;
  const Double_t eleEtaMin = -3;
  const Double_t eleEtaMax =  3;
  
  // jet requirements
  const Double_t jetPtMin = 10;

  // photon requirements
  const Double_t photonEtMin = 9;
      
  // good PV requirements
  const UInt_t   minNTracksFit = 0;
  const Double_t minNdof       = 4;
  const Double_t maxAbsZ       = 24;
  const Double_t maxRho        = 2;
  
  //
  // setup analysis object
  //
  Analysis *ana = new Analysis;
  ana->SetUseHLT(kTRUE);
  if(nevents>0) 
    ana->SetProcessNEvents(nevents);
  printf("\nRely on Catalog: %s\n",catalogDir);
  printf("  -> Book: %s  Dataset: %s  Fileset: %s <-\n\n",book,dataset,fileset);
  Catalog *c = new Catalog(catalogDir);
  Dataset *d = NULL;
  d = c->FindDataset(book,dataset,fileset);
  ana->AddDataset(d);
    
  //
  // setup ntupler module
  //
  NtuplerMod *mymod = new NtuplerMod;
  mymod->SetOutputName(output);          // output ntuple file name
  mymod->SetIsData(isData);              // toggle data specific or MC specific procedures
  mymod->SetUseGen(useGen);              // use generator info
  mymod->SetSkipIfHLTFail(skipHLTFail);  // skip to next event if no HLT accept
  mymod->SetFSRMode(fsrmode);
  mymod->SetMuPtMin(muPtMin);
  mymod->SetMuPtMax(muPtMax);
  mymod->SetMuEtaMin(muEtaMin);
  mymod->SetMuEtaMax(muEtaMax);
  mymod->SetEleEtMin(eleEtMin);
  mymod->SetEleEtMax(eleEtMax);
  mymod->SetEleEtaMin(eleEtaMin);
  mymod->SetEleEtaMax(eleEtaMax);
  mymod->SetJetPtMin(jetPtMin);
  mymod->SetPhotonEtMin(photonEtMin);
  mymod->SetMinNTracksFit(minNTracksFit);
  mymod->SetMinNdof(minNdof);
  mymod->SetMaxAbsZ(maxAbsZ);
  mymod->SetMaxRho(maxRho);

  // Jet corrections
  mymod->AddJetCorr("START50_V15_L1FastJet_AK5PF.txt");
  mymod->AddJetCorr("START50_V15_L2Relative_AK5PF.txt");
  mymod->AddJetCorr("START50_V15_L3Absolute_AK5PF.txt");
  if(isData)
    mymod->AddJetCorr("START50_V15_L2L3Residual_AK5PF.txt");
  
  //
  // SingleMu
  //
  mymod->AddTrigger("HLT_Mu15_v14",kHLT_Mu15,"hltSingleMu15L3Filtered15",kHLT_Mu15_MuObj);
  mymod->AddTrigger("HLT_Mu15_eta2p1_v3",kHLT_Mu15_eta2p1,"hltL3fL1sMu7L1fEta2p1L2fEta2p1f7L3FilteredEta2p1Filtered15",kHLT_Mu15_eta2p1_MuObj);
  
  //
  // SingleElectron
  //
  mymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_v11",kHLT_Ele17_CaloIdL_CaloIsoVL,"hltEle17CaloIdLCaloIsoVLPixelMatchFilter",kHLT_Ele17_CaloIdL_CaloIsoVL_EleObj,22);
  mymod->AddTrigger("HLT_Ele22_CaloIdL_CaloIsoVL_v3",kHLT_Ele22_CaloIdL_CaloIsoVL,"hltEle22CaloIdLCaloIsoVLPixelMatchFilter",kHLT_Ele22_CaloIdL_CaloIsoVL_EleObj);
  mymod->AddTrigger("HLT_Ele22_CaloIdL_CaloIsoVL_v4",kHLT_Ele22_CaloIdL_CaloIsoVL,"hltEle22CaloIdLCaloIsoVLPixelMatchFilter",kHLT_Ele22_CaloIdL_CaloIsoVL_EleObj);
  mymod->AddTrigger("HLT_Ele22_CaloIdL_CaloIsoVL_v5",kHLT_Ele22_CaloIdL_CaloIsoVL,"hltEle22CaloIdLCaloIsoVLPixelMatchFilter",kHLT_Ele22_CaloIdL_CaloIsoVL_EleObj);
  
  
  mymod->SetPrintHLT(kFALSE); // print HLT table at start of analysis?
  
  ana->AddSuperModule(mymod); 
    
  //
  // run analysis after successful initialisation
  //
  ana->Run(!gROOT->IsBatch());
}

//==================================================================================================
/*
 * Run on a single BAMBU file (mainly for testing purposes)
 *
 */
void runNtupler(
    const char *file   = "/castor/cern.ch/user/p/paus/filefi/026/s12-zmmm20-v15/009E6E12-F36E-E111-9AED-003048678FA0.root",
    const char *output = "ntuple.root",
    Bool_t isData      = kFALSE,
    Int_t  useGen      = 2,
    Int_t  fsrmode     = 0,
    Int_t  nevents     = -1,
    Bool_t skipHLTFail = kFALSE
  )
{
  gDebugMask  = Debug::kAnalysis;  // debug message category
  gDebugLevel = 1;                 // higher level allows more messages to print
  
  // muon kinematics
  const Double_t muPtMin  = 0;
  const Double_t muPtMax  = 7000;
  const Double_t muEtaMin = -3;
  const Double_t muEtaMax =  3;

  // electron kinematics
  const Double_t eleEtMin  = 9;
  const Double_t eleEtMax  = 7000;
  const Double_t eleEtaMin = -3;
  const Double_t eleEtaMax =  3;
  
  // jet requirements
  const Double_t jetPtMin = 10;

  // photon requirements
  const Double_t photonEtMin = 9;
      
  // good PV requirements
  const UInt_t   minNTracksFit = 0;
  const Double_t minNdof       = 4;
  const Double_t maxAbsZ       = 24;
  const Double_t maxRho        = 2;
  
  //
  // setup analysis object
  //
  Analysis *ana = new Analysis;
  ana->SetUseHLT(kTRUE);
  if(nevents>0) 
    ana->SetProcessNEvents(nevents);
  ana->AddFile(file);
    
  //
  // setup ntupler module
  //
  NtuplerMod *mymod = new NtuplerMod;
  mymod->SetOutputName(output);          // output ntuple file name
  mymod->SetIsData(isData);              // toggle data specific or MC specific procedures
  mymod->SetUseGen(useGen);              // use generator info
  mymod->SetSkipIfHLTFail(skipHLTFail);  // skip to next event if no HLT accept
  mymod->SetFSRMode(fsrmode);
  mymod->SetMuPtMin(muPtMin);
  mymod->SetMuPtMax(muPtMax);
  mymod->SetMuEtaMin(muEtaMin);
  mymod->SetMuEtaMax(muEtaMax);
  mymod->SetEleEtMin(eleEtMin);
  mymod->SetEleEtMax(eleEtMax);
  mymod->SetEleEtaMin(eleEtaMin);
  mymod->SetEleEtaMax(eleEtaMax);
  mymod->SetJetPtMin(jetPtMin);
  mymod->SetPhotonEtMin(photonEtMin);
  mymod->SetMinNTracksFit(minNTracksFit);
  mymod->SetMinNdof(minNdof);
  mymod->SetMaxAbsZ(maxAbsZ);
  mymod->SetMaxRho(maxRho);

  // Jet corrections
  mymod->AddJetCorr("START50_V15_L1FastJet_AK5PF.txt");
  mymod->AddJetCorr("START50_V15_L2Relative_AK5PF.txt");
  mymod->AddJetCorr("START50_V15_L3Absolute_AK5PF.txt");
  if(isData)
    mymod->AddJetCorr("START50_V15_L2L3Residual_AK5PF.txt");

  //
  // SingleMu
  //
  mymod->AddTrigger("HLT_Mu15_v14",kHLT_Mu15,"hltSingleMu15L3Filtered15",kHLT_Mu15_MuObj);
  mymod->AddTrigger("HLT_Mu15_eta2p1_v3",kHLT_Mu15_eta2p1,"hltL3fL1sMu7L1fEta2p1L2fEta2p1f7L3FilteredEta2p1Filtered15",kHLT_Mu15_eta2p1_MuObj);
  
  //
  // SingleElectron
  //
  mymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_v11",kHLT_Ele17_CaloIdL_CaloIsoVL,"hltEle17CaloIdLCaloIsoVLPixelMatchFilter",kHLT_Ele17_CaloIdL_CaloIsoVL_EleObj,22);
  mymod->AddTrigger("HLT_Ele22_CaloIdL_CaloIsoVL_v3",kHLT_Ele22_CaloIdL_CaloIsoVL,"hltEle22CaloIdLCaloIsoVLPixelMatchFilter",kHLT_Ele22_CaloIdL_CaloIsoVL_EleObj);
  mymod->AddTrigger("HLT_Ele22_CaloIdL_CaloIsoVL_v4",kHLT_Ele22_CaloIdL_CaloIsoVL,"hltEle22CaloIdLCaloIsoVLPixelMatchFilter",kHLT_Ele22_CaloIdL_CaloIsoVL_EleObj);
  mymod->AddTrigger("HLT_Ele22_CaloIdL_CaloIsoVL_v5",kHLT_Ele22_CaloIdL_CaloIsoVL,"hltEle22CaloIdLCaloIsoVLPixelMatchFilter",kHLT_Ele22_CaloIdL_CaloIsoVL_EleObj);

  
  mymod->SetPrintHLT(kFALSE); // print HLT table at start of analysis?
  
  ana->AddSuperModule(mymod); 
  
  //
  // run analysis after successful initialisation
  //
  ana->Run(!gROOT->IsBatch());
}
