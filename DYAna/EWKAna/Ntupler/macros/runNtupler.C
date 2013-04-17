#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/Catalog/interface/Catalog.h"
#include "MitAna/TreeMod/interface/Analysis.h"
#include "EWKAna/Ntupler/interface/NtuplerMod.hh"
#endif

using namespace mithep;
   
//==================================================================================================
/*
 * Run on a BAMBU fileset
 *
 * Example usage:
 *   root -l -q -b runNtupler.C+\(\"0000\",\"s11-zeem20-v11-pu\",\"cern/filefi/021\",\"/home/mitprod/catalog\",0,1,0,-1,0\)
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
    const Bool_t  useGen,       // flag to look at generator info
    const Int_t   fsrmode,      // indicates how generator level FSR effects are stored
    const Int_t   nevents,      // number of events to process
    const Bool_t  skipHLTFail,  // skip events if no HLT accept
    const Bool_t  is2011        // is this 2011 data or MC?
  )
{
  gDebugMask  = Debug::kAnalysis;  // debug message category
  gDebugLevel = 1;                 // higher level allows more messages to print
 
  char output[100];
  sprintf(output,"%s_%s_ntuple.root",dataset,fileset); 
  
  // muon kinematics
  const Double_t muPtMin  = 10;
  const Double_t muPtMax  = 7000;
  const Double_t muEtaMin = -3;
  const Double_t muEtaMax =  3;

  // electron kinematics
  const Double_t eleEtMin  = 9;
  const Double_t eleEtMax  = 7000;
  const Double_t eleEtaMin = -3;
  const Double_t eleEtaMax =  3;

  // mass region
  const Double_t massMin = 10;
  const Double_t massMax = 7000;
  
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
  mymod->SetMassMin(massMin);
  mymod->SetMassMax(massMax);
  mymod->SetJetPtMin(jetPtMin);
  mymod->SetPhotonEtMin(photonEtMin);
  mymod->SetMinNTracksFit(minNTracksFit);
  mymod->SetMinNdof(minNdof);
  mymod->SetMaxAbsZ(maxAbsZ);
  mymod->SetMaxRho(maxRho);

  // Jet corrections
  if(is2011) {
    mymod->AddJetCorr("START42_V12_AK5PF_L1FastJet.txt");
    mymod->AddJetCorr("START42_V12_AK5PF_L2Relative.txt");
    mymod->AddJetCorr("START42_V12_AK5PF_L3Absolute.txt");
    if(isData)
      mymod->AddJetCorr("START42_V12_AK5PF_L2L3Residual.txt");
  } else {
    mymod->AddJetCorr("START53_V7F_L1FastJet_AK5PF.txt");
    mymod->AddJetCorr("START53_V7F_L2Relative_AK5PF.txt");
    mymod->AddJetCorr("START53_V7F_L3Absolute_AK5PF.txt");
    if(isData)
      mymod->AddJetCorr("START53_V7F_L2L3Residual_AK5PF.txt");  
  }
  
  if(is2011) {
#include "HLT2011.hh"
  } else {
#include "HLT2012.hh"  
  }
  
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
    const char *file   = "/castor/cern.ch/user/p/pharris/test/XX-MITDATASET-XX_000_10_1_5jt.root",
    const char *output = "ntuple.root",
    Bool_t isData      = kFALSE,
    Bool_t useGen      = kTRUE,
    Int_t  fsrmode     = 0,
    Int_t  nevents     = -1,
    Bool_t skipHLTFail = kFALSE,
    Bool_t is2011      = kTRUE
  )
{
  gDebugMask  = Debug::kAnalysis;  // debug message category
  gDebugLevel = 1;                 // higher level allows more messages to print
  
  // muon kinematics
  const Double_t muPtMin  = 10;
  const Double_t muPtMax  = 7000;
  const Double_t muEtaMin = -3;
  const Double_t muEtaMax =  3;

  // electron kinematics
  const Double_t eleEtMin  = 9;
  const Double_t eleEtMax  = 7000;
  const Double_t eleEtaMin = -3;
  const Double_t eleEtaMax =  3;

  // mass region
  const Double_t massMin = 10;
  const Double_t massMax = 7000;
    
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
  mymod->SetMassMin(massMin);
  mymod->SetMassMax(massMax);
  mymod->SetJetPtMin(jetPtMin);
  mymod->SetPhotonEtMin(photonEtMin);
  mymod->SetMinNTracksFit(minNTracksFit);
  mymod->SetMinNdof(minNdof);
  mymod->SetMaxAbsZ(maxAbsZ);
  mymod->SetMaxRho(maxRho);

  // Jet corrections
  if(is2011) {
    mymod->AddJetCorr("START42_V12_AK5PF_L1FastJet.txt");
    mymod->AddJetCorr("START42_V12_AK5PF_L2Relative.txt");
    mymod->AddJetCorr("START42_V12_AK5PF_L3Absolute.txt");
    if(isData)
      mymod->AddJetCorr("START42_V12_AK5PF_L2L3Residual.txt");
  } else {
    mymod->AddJetCorr("START53_V7F_L1FastJet_AK5PF.txt");
    mymod->AddJetCorr("START53_V7F_L2Relative_AK5PF.txt");
    mymod->AddJetCorr("START53_V7F_L3Absolute_AK5PF.txt");
    if(isData)
      mymod->AddJetCorr("START53_V7F_L2L3Residual_AK5PF.txt");  
  }

  if(is2011) {
#include "HLT2011.hh"
  } else {
#include "HLT2012.hh"  
  }  
  
  mymod->SetPrintHLT(kFALSE); // print HLT table at start of analysis?
  
  ana->AddSuperModule(mymod); 
  
  //
  // run analysis after successful initialisation
  //
  ana->Run(!gROOT->IsBatch());
}
