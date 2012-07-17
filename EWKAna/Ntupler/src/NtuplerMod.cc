#include "EWKAna/Ntupler/interface/NtuplerMod.hh"
#include "MitAna/DataTree/interface/Names.h"
#include "MitAna/DataTree/interface/MCParticleCol.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/Track.h"
#include "MitAna/DataTree/interface/Vertex.h"
#include "MitAna/DataTree/interface/BeamSpotCol.h"
#include "MitAna/DataTree/interface/PileupInfoCol.h"
#include "MitAna/DataTree/interface/PileupEnergyDensityCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/PFJet.h"
#include "MitAna/DataTree/interface/Photon.h"
#include "MitAna/DataTree/interface/TriggerTable.h"
#include "MitAna/DataTree/interface/TriggerObjectsTable.h"
#include "MitAna/DataTree/interface/TriggerName.h"
#include "MitAna/DataTree/interface/Met.h"
#include "MitAna/DataTree/interface/PFMet.h"
#include "MitAna/DataTree/interface/DecayParticle.h"
#include "MitAna/DataTree/interface/StableData.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitPhysics/Utils/interface/ElectronTools.h"
#include "MitPhysics/Utils/interface/IsolationTools.h"
#include "MitPhysics/Utils/interface/PhotonTools.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLorentzVector.h>
#include <vector>

using namespace mithep;

ClassImp(mithep::NtuplerMod)

NtuplerMod::NtuplerMod(const char *name, const char *title):
  BaseMod        (name,title),
  fOutputFile    (0),
  fOutputName    ("ntuple.root"),
  fPartName      (Names::gkMCPartBrn),
  fMCEvtInfoName (Names::gkMCEvtInfoBrn),
  fMuonName      (Names::gkMuonBrn),
  fElectronName  (Names::gkElectronBrn),
  fPrimVtxName   (Names::gkPVBrn),
  fBeamSpotName  (Names::gkBeamSpotBrn),
  fPFJetName     (Names::gkPFJetBrn),
  fPhotonName    (Names::gkPhotonBrn),
  fTrigMaskName  (Names::gkHltBitBrn),
  fPFMetName     ("PFMet"),
  fConversionName(Names::gkMvfConversionBrn),
  fPileupName    (Names::gkPileupInfoBrn),
  fPUEnergyDensityName(Names::gkPileupEnergyDensityBrn),
  fPFCandidateName(Names::gkPFCandidatesBrn),
  fTracksName    (Names::gkTrackBrn),
  fParticles     (0),
  fMCEvtInfo     (0),
  fMuons         (0),
  fElectrons     (0),
  fPrimVerts     (0),
  fBeamSpot      (0),
  fPFJets        (0),
  fPhotons       (0),  
  fTrigMask      (0),
  fPFMet         (0),
  fConversions   (0),
  fPileup        (0),
  fPUEnergyDensity(0),
  fPFCandidates  (0),
  fPFPileUp      (0),
  fPFNoPileUp    (0),
  fTracks        (0),
  fIsData        (kFALSE),
  fUseGen        (0),
  fPrintTable    (kFALSE),
  fSkipIfHLTFail (kFALSE),
  fMuPtMin       (15),
  fMuPtMax       (1000),
  fMuEtaMin      (-3),
  fMuEtaMax      (3),
  fEleEtMin      (15),
  fEleEtMax      (1000),
  fEleEtaMin     (-3),
  fEleEtaMax     (3),
  fJetPtMin      (15),
  fPhotonEtMin   (10),
  fMinNTracksFit (0),
  fMinNdof       (4),
  fMaxAbsZ       (24),
  fMaxRho        (2),
  fEventTree     (0),
  fFSRMode       (0),
  fJetCorrector  (0),
  fEleMVA        (0)
{
  // Constructor
  
  // Don't write TObject part of the objects
  TEventInfo::Class()->IgnoreTObjectStreamer();
  TGenInfo::Class()->IgnoreTObjectStreamer();
  TMuon::Class()->IgnoreTObjectStreamer();
  TElectron::Class()->IgnoreTObjectStreamer();
  TJet::Class()->IgnoreTObjectStreamer();
  TPhoton::Class()->IgnoreTObjectStreamer();
  TVertex::Class()->IgnoreTObjectStreamer();
}

//--------------------------------------------------------------------------------------------------
NtuplerMod::~NtuplerMod()
{
  // Destructor
}	

//--------------------------------------------------------------------------------------------------      
void NtuplerMod::Begin()
{
}

//--------------------------------------------------------------------------------------------------
void NtuplerMod::SlaveBegin()
{
  //
  // Request BAMBU branches
  //
  ReqBranch(fPartName,            fParticles); 
  ReqBranch(fMCEvtInfoName,       fMCEvtInfo);
  ReqBranch(fMuonName,            fMuons);
  ReqBranch(fElectronName,        fElectrons);
  ReqBranch(fPrimVtxName,         fPrimVerts);
  ReqBranch(fBeamSpotName,        fBeamSpot);
  ReqBranch(fPFJetName,           fPFJets);
  ReqBranch(fTrigMaskName,        fTrigMask);
  ReqBranch(fPFMetName,           fPFMet);
  ReqBranch(fPhotonName,          fPhotons);
  ReqBranch(fConversionName,      fConversions);
  ReqBranch(fPileupName,          fPileup);
  ReqBranch(fPUEnergyDensityName, fPUEnergyDensity);
  ReqBranch(fPFCandidateName,     fPFCandidates); 
  ReqBranch(fTracksName,          fTracks);

  // Pileup and NoPileup collections of PFCandidates
  fPFPileUp   = new PFCandidateOArr;
  fPFNoPileUp = new PFCandidateOArr;
   
  //
  // Set up arrays
  //
  fMuonArr     = new TClonesArray("mithep::TMuon");	assert(fMuonArr);  
  fElectronArr = new TClonesArray("mithep::TElectron"); assert(fElectronArr);
  fPFJetArr    = new TClonesArray("mithep::TJet");	assert(fPFJetArr);
  fPhotonArr   = new TClonesArray("mithep::TPhoton");	assert(fPhotonArr);
  fPVArr       = new TClonesArray("mithep::TVertex");	assert(fPVArr);
  
  //
  // Create output file
  //
  fOutputFile = new TFile(fOutputName, "RECREATE");
  
  //
  // Initialize data trees and structs
  // 
  fEventTree = new TTree("Events","Events");

  fEventTree->Branch("Info",&fEventInfo);
  if(!fIsData && fUseGen)
    fEventTree->Branch("Gen",&fGenInfo);
  
  fEventTree->Branch("Muon",    &fMuonArr);
  fEventTree->Branch("Electron",&fElectronArr);
  fEventTree->Branch("PFJet",   &fPFJetArr);
  fEventTree->Branch("Photon",  &fPhotonArr);
  fEventTree->Branch("PV",      &fPVArr);
  
  //
  // Set up jet corrections for PF jets
  //
  std::vector<JetCorrectorParameters> correctionParameters;
  for(UInt_t icorr=0; icorr<fJetCorrParsv.size(); icorr++)
    correctionParameters.push_back(JetCorrectorParameters(fJetCorrParsv[icorr].Data()));
    
  // initialize jet corrector class
  fJetCorrector = new FactorizedJetCorrector(correctionParameters);
  
  // initialize electron MVA
  fEleMVA = new ElectronIDMVA();

  fEleMVA->Initialize("BDTG method",
                      "Subdet0LowPt_WithIPInfo_BDTG.weights.xml",   // 0     < |eta| < 1,     10 < pT < 20
		      "Subdet1LowPt_WithIPInfo_BDTG.weights.xml",   // 1     < |eta| < 1.479, 10 < pT < 20
		      "Subdet2LowPt_WithIPInfo_BDTG.weights.xml",   // 1.479 < |eta| < 2.5,   10 < pT < 20
		      "Subdet0HighPt_WithIPInfo_BDTG.weights.xml",  // 0     < |eta| < 1,     pT > 20
		      "Subdet1HighPt_WithIPInfo_BDTG.weights.xml",  // 1     < |eta| < 1.479, pT > 20
		      "Subdet2HighPt_WithIPInfo_BDTG.weights.xml",  // 1.479 < |eta| < 2.5,   pT > 20
		      ElectronIDMVA::kWithIPInfo);
 
  fLastRunLumi = RunLumiRangeMap::RunLumiPairType(0,0);
}

//--------------------------------------------------------------------------------------------------
void NtuplerMod::SlaveTerminate()
{
  //
  // Save to ROOT file
  //
  fEventTree->Print();

  fOutputFile->Write();
  fOutputFile->Close();
  
  delete fPFPileUp;
  delete fPFNoPileUp;
  
  delete fMuonArr;
  delete fElectronArr;
  delete fPFJetArr;
  delete fPhotonArr;
  delete fPVArr;
 
  delete fJetCorrector;
  delete fEleMVA;
   
  //
  // Dump JSON file
  //
  TString jsonfname = fOutputName+TString(".json");
  if(fIsData)
    fRunLumiSet.DumpJSONFile(jsonfname.Data());
}

//--------------------------------------------------------------------------------------------------
void NtuplerMod::Terminate()
{
}

//--------------------------------------------------------------------------------------------------
void NtuplerMod::BeginRun()
{
  if(HasHLTInfo() && fPrintTable) { GetHLTTable()->Print(); }
}

//--------------------------------------------------------------------------------------------------
void NtuplerMod::EndRun()
{
}

//--------------------------------------------------------------------------------------------------
void NtuplerMod::Process()
{
  RunLumiRangeMap::RunLumiPairType rl(GetEventHeader()->RunNum(), GetEventHeader()->LumiSec());
  if(rl!=fLastRunLumi) {
    fLastRunLumi = rl;
    fRunLumiSet.Add(rl);
  }

  //
  // Load branches
  //
  if(!fIsData && fUseGen) LoadBranch(fPartName);
  if(!fIsData && fUseGen) LoadBranch(fMCEvtInfoName);
  LoadBranch(fMuonName);
  LoadBranch(fElectronName);
  LoadBranch(fPrimVtxName);
  LoadBranch(fBeamSpotName);
  LoadBranch(fPFJetName);
  LoadBranch(fTrigMaskName);
  LoadBranch(fPFMetName); 
  LoadBranch(fPhotonName);
  LoadBranch(fConversionName);
  LoadBranch(fPUEnergyDensityName);
  LoadBranch(fPFCandidateName);
  LoadBranch(fTracksName);
  if(!fIsData)
    LoadBranch(fPileupName);

  //
  // Scan generator info
  //
  if(fUseGen) {
    if(fUseGen==1) FillGenH();
    if(fUseGen==2) FillGenZ();
    if(fUseGen==3) FillGenW();
    if(fUseGen==4) FillGenWW();
    if(fUseGen==5) FillGenVZ();
    if(fUseGen==6) FillGenWjets();
    if(fUseGen==7) FillGenZjets();
  }
  
  //
  // Get HLT info. Trigger objects can be matched by name to the corresponding trigger that passed.
  //
  ULong64_t trigbits=0;
  if(HasHLTInfo()) {
    const TriggerTable *hltTable = GetHLTTable();
    assert(hltTable);
    for(UInt_t itrig=0; itrig<fTriggerNamesv.size(); itrig++) {
      const TriggerName *trigname = hltTable->Get(fTriggerNamesv[itrig].Data());
      if(!trigname) continue;
      if(fTrigMask->At(trigname->Id())) { trigbits |= fTriggerIdsv[itrig]; }
    }  
  }
  if(fSkipIfHLTFail && (trigbits==0))
    return;
  
  IncNEventsProcessed();
  
  fMuonArr->Clear();
  fElectronArr->Clear();
  fPFJetArr->Clear();
  fPhotonArr->Clear();
  fPVArr->Clear();


  //
  // Get beam spot. If no beam spot information is available, default the coordinates to 99999
  //
  Double_t bsx=99999, bsy=99999, bsz=99999;
  if(fBeamSpot) {
    if(fBeamSpot->GetEntries() > 1) 
      std::cout << "********** More than 1 beam spot! **********" << std::endl;
    const BeamSpot *bs = fBeamSpot->At(0);
    bsx = bs->X();
    bsy = bs->Y();
    bsz = bs->Z();
  }

  //
  // Get primary vertices
  // Assumes primary vertices are ordered by sum-pT^2 (as should be in CMSSW)
  // NOTE: if no PV is found from fitting tracks, the beamspot is used
  //
  fVertex = 0;
  Bool_t hasGoodPV = kFALSE;  
  for(UInt_t i=0; i<fPrimVerts->GetEntries(); ++i) {
    const Vertex *pv = fPrimVerts->At(i);
    
    // Select best PV for corrected d0; if no PV passing cuts, the first PV in the collection will be used
    //if(!pv->IsValid()) continue;
    if(pv->NTracksFit()     < fMinNTracksFit) continue;
    if(pv->Ndof()	    < fMinNdof)	      continue;
    if(fabs(pv->Z())	    > fMaxAbsZ)	      continue;
    if(pv->Position().Rho() > fMaxRho)	      continue;    
    hasGoodPV = kTRUE;
    
    FillPV(pv);
    if(!fVertex) fVertex = pv;
  }
  if(!fVertex) fVertex = fPrimVerts->At(0);
   
  //
  // Separate PF candidates into those associated
  // with PV and those associated with PU
  //
  separatePileUp(kTRUE);
  
  //
  // Loop through muons (and general tracks if desired).
  //
  vector<const Muon*> muonv;    // array of pointers to preselected muons ... 
  
  assert(fMuons);
  for(UInt_t i=0; i<fMuons->GetEntries(); ++i) {
    const Muon *mu = fMuons->At(i); 
    if(!mu->HasTrk()) continue; 
    
    // Use tracker tracks for kinematics when available
    const Track *muTrk=0;
    if(mu->HasTrackerTrk())     { muTrk = mu->TrackerTrk(); }
    else if(mu->HasGlobalTrk()) { muTrk = mu->GlobalTrk(); }
    else                        { muTrk = mu->StandaloneTrk(); } 
              
    if((muTrk->Eta() < fMuEtaMin) || (muTrk->Eta() > fMuEtaMax)) continue;   
    if((muTrk->Pt()  > fMuPtMax)  || (muTrk->Pt()  < fMuPtMin))  continue;
    
    FillMuon(mu);  
  }
       
  assert(fTracks);
  for(UInt_t i=0; i<fTracks->GetEntries(); ++i) {
    const Track *track = fTracks->At(i);
  
    // Check that the track is not associated with a muon.
    // If it is, skip to next track...
    Bool_t isMuon = kFALSE;
    for(UInt_t j=0; j<fMuons->GetEntries(); ++j) {
      if(track == (fMuons->At(j)->TrackerTrk())) isMuon = kTRUE;
    }
    if(isMuon) continue;
  
    if(track->Pt() < 20) continue;  // pT cut
    if((track->Eta() < fMuEtaMin) || (track->Eta() > fMuEtaMax)) continue;  // eta cut
  
    FillMuon(track);
  }


  //
  // Loop through electrons.
  //
  vector<const Electron*> elev;  // array of pointers to preselected electrons ... 
  ElectronTools eleTools;        // helper class for electron ID decisions
  
  assert(fElectrons);                
  for(UInt_t i=0; i<fElectrons->GetEntries(); ++i) {
    const Electron *ele = fElectrons->At(i);  
   
    if((ele->Pt()  < fEleEtMin)  || (ele->Pt()  > fEleEtMax))  continue;  // electron pT cut
    if((ele->Eta() < fEleEtaMin) || (ele->Eta() > fEleEtaMax)) continue;  // electron eta cut
    if(!eleTools.PassSpikeRemovalFilter(ele))                  continue;  // spike cleaning
        
    FillElectron(ele);  // fill electron data object
  }
  
  //
  // Loop through jets
  //
  assert(fPFJets);
  for(UInt_t i=0; i<fPFJets->GetEntries(); ++i) {
    const PFJet *jet = fPFJets->At(i);
    
    const FourVectorM rawMom = jet->RawMom();
    fJetCorrector->setJetEta(rawMom.Eta());
    fJetCorrector->setJetPt(rawMom.Pt());
    fJetCorrector->setJetPhi(rawMom.Phi());
    fJetCorrector->setJetE(rawMom.E());
    fJetCorrector->setRho(fPUEnergyDensity->At(0)->RhoHighEta());
    fJetCorrector->setJetA(jet->JetArea());
    fJetCorrector->setJetEMF(-99.0);     
    
    if( ((jet->RawMom().Pt())*(fJetCorrector->getCorrection()) > fJetPtMin)
         || (jet->TrackCountingHighEffBJetTagsDisc()!=-100) ) {
       FillJet(jet);
    }
  }  

    
  //
  // Loop through photons
  //
  assert(fPhotons);
  for(UInt_t i=0; i<fPhotons->GetEntries(); ++i) {
    const Photon *pho= fPhotons->At(i);
    if(pho->SCluster()->Et() > fPhotonEtMin) { FillPhoton(pho); }
  } 
  
  //
  // Compute MET
  //
  TLorentzVector pfmet; pfmet.SetPxPyPzE(fPFMet->At(0)->Mex(),fPFMet->At(0)->Mey(),0,0);
  Double_t trkMetx=0, trkMety=0, trkSumET=0;
  
  assert(fPFCandidates);
  for(UInt_t i=0; i<fPFCandidates->GetEntries(); ++i) {
    const Double_t trkDzCut  = 0.1;
    
    const PFCandidate *pfcand = fPFCandidates->At(i);
    if( (pfcand->HasTrackerTrk() && (fabs(pfcand->TrackerTrk()->DzCorrected(*fVertex))<trkDzCut)) ||
        (pfcand->HasGsfTrk()     && (fabs(pfcand->GsfTrk()->DzCorrected(*fVertex))<trkDzCut)) ) {
      
      trkMetx  -= pfcand->Px();
      trkMety  -= pfcand->Py();
      trkSumET += pfcand->Pt();
    }
  }
  TLorentzVector trkmet; trkmet.SetPxPyPzE(trkMetx,trkMety,0,0);
        
  //
  // Identify bunch crossing 0 in PU info
  //
  Int_t ibx0=-1, ibxm=-1, ibxp=-1;
  if(!fIsData) {
    assert(fPileup);    
    for(UInt_t i=0; i<fPileup->GetEntries(); ++i) {
      if(fPileup->At(i)->GetBunchCrossing() == 0) ibx0=i;
      if(fPileup->At(i)->GetBunchCrossing() ==-1) ibxm=i;
      if(fPileup->At(i)->GetBunchCrossing() == 1) ibxp=i;
    }
  }
 
  //
  // Fill event info tree
  //    
  fEventInfo.runNum       = GetEventHeader()->RunNum();
  fEventInfo.evtNum       = GetEventHeader()->EvtNum();
  fEventInfo.lumiSec      = GetEventHeader()->LumiSec();
  fEventInfo.nPU          = (ibx0>-1) ? fPileup->At(ibx0)->GetPU_NumInteractions() : 0;
  fEventInfo.nPUminus     = (ibxm>-1) ? fPileup->At(ibxm)->GetPU_NumInteractions() : 0;
  fEventInfo.nPUplus      = (ibxp>-1) ? fPileup->At(ibxp)->GetPU_NumInteractions() : 0;
  fEventInfo.triggerBits  = trigbits;
  fEventInfo.pvx          = fVertex->X();
  fEventInfo.pvy          = fVertex->Y();
  fEventInfo.pvz          = fVertex->Z();
  fEventInfo.bsx          = bsx;
  fEventInfo.bsy          = bsy;
  fEventInfo.bsz          = bsz;
  fEventInfo.pfMET        = pfmet.Pt();
  fEventInfo.pfMETphi     = pfmet.Phi();
  fEventInfo.pfSumET      = fPFMet->At(0)->SumEt();
  fEventInfo.trkMET       = trkmet.Pt();
  fEventInfo.trkMETphi    = trkmet.Phi();
  fEventInfo.trkSumET     = trkSumET;
  fEventInfo.rhoLowEta    = fPUEnergyDensity->At(0)->RhoLowEta();
  fEventInfo.rhoHighEta   = fPUEnergyDensity->At(0)->RhoHighEta();
  fEventInfo.hasGoodPV    = hasGoodPV;
  
  // Fill the tree
  fEventTree->Fill();
}

//--------------------------------------------------------------------------------------------------
void NtuplerMod::FillMuon(const Muon *mu)
{
  assert(mu);

  TClonesArray &rMuonArr = *fMuonArr;
  assert(rMuonArr.GetEntries() < rMuonArr.GetSize());
  const Int_t index = rMuonArr.GetEntries();  
  new(rMuonArr[index]) TMuon();
  TMuon *pMuon = (TMuon*)rMuonArr[index];
  
  // Use tracker track when available
  const Track *muTrk=0;
  if(mu->HasTrackerTrk())         { muTrk = mu->TrackerTrk(); }
  else if(mu->HasGlobalTrk())     { muTrk = mu->GlobalTrk(); }
  else if(mu->HasStandaloneTrk()) { muTrk = mu->StandaloneTrk(); }
  assert(muTrk);                  
  
  pMuon->pt         = muTrk->Pt();
  pMuon->ptErr      = muTrk->PtErr();
  pMuon->eta        = muTrk->Eta();
  pMuon->phi        = muTrk->Phi();
  pMuon->trkIso03   = mu->IsoR03SumPt();
  pMuon->emIso03    = mu->IsoR03EmEt();
  pMuon->hadIso03   = mu->IsoR03HadEt();  
  pMuon->pfChIso03  = computePFIso(muTrk, 0.0, 0.3, 0.0, 4);
  pMuon->pfNeuIso03 = computePFIso(muTrk, 0.5, 0.3, 0.0, 2);
  pMuon->pfGamIso03 = computePFIso(muTrk, 0.5, 0.3, 0.0, 3);
  pMuon->puIso03    = computePUIso(muTrk, 0.0, 0.3, 0);
  pMuon->pfChIso04  = computePFIso(muTrk, 0.0, 0.4, 0.0, 4);
  pMuon->pfNeuIso04 = computePFIso(muTrk, 0.5, 0.4, 0.0, 2);
  pMuon->pfGamIso04 = computePFIso(muTrk, 0.5, 0.4, 0.0, 3);
  pMuon->puIso04    = computePUIso(muTrk, 0.0, 0.4, 0); 
  pMuon->d0         = muTrk->D0Corrected(*fVertex);
  pMuon->dz         = muTrk->DzCorrected(*fVertex);
  pMuon->tkNchi2    = (mu->HasTrackerTrk()) ? mu->TrackerTrk()->RChi2() : 0;
  
  if(mu->HasGlobalTrk())          { pMuon->muNchi2 = mu->GlobalTrk()->RChi2();     }
  else if(mu->HasStandaloneTrk()) { pMuon->muNchi2 = mu->StandaloneTrk()->RChi2(); }
  else if(mu->HasTrackerTrk())    { pMuon->muNchi2 = mu->TrackerTrk()->RChi2();    }
  
  pMuon->trkKink    = mu->TrkKink();
  pMuon->gblKink    = mu->GlbKink();
  pMuon->mva        = 0;
  pMuon->q          = muTrk->Charge();
  pMuon->nValidHits = mu->NValidHits();
  
  pMuon->qualityBits = mu->Quality().QualityMask().Mask();

  //
  // NOTE:
  // It is possible for a muon to be TK+SA. The muon reco associates a TK with a SA if
  // chamber matches for the TK and hits for the SA share DetIDs
  //   (see hypernews thread: https://hypernews.cern.ch/HyperNews/CMS/get/csa08-muons/57/2/1/1/1.html)
  //	      
  pMuon->typeBits = 0;
  if(mu->IsGlobalMuon())     { pMuon->typeBits |= kGlobal; }
  if(mu->IsTrackerMuon())    { pMuon->typeBits |= kTracker; }
  if(mu->IsStandaloneMuon()) { pMuon->typeBits |= kStandalone; }
  if(mu->IsPFMuon())         { pMuon->typeBits |= kPFMuon; }
  
  pMuon->nTkHits      = (mu->HasTrackerTrk()) ? mu->TrackerTrk()->NHits() : 0;
  pMuon->nPixHits     = muTrk->NPixelHits();
  pMuon->nTkLayers    = mu->NTrkLayersHit();
  pMuon->nPixLayers   = mu->NPxlLayersHit();
  pMuon->nSeg         = mu->NSegments();
  pMuon->nMatch       = mu->NMatches();
  pMuon->hltMatchBits = MatchHLT(muTrk->Eta(),muTrk->Phi());
  
  pMuon->staPt  = (mu->HasStandaloneTrk()) ? mu->StandaloneTrk()->Pt()  : -1;  
  pMuon->staEta = (mu->HasStandaloneTrk()) ? mu->StandaloneTrk()->Eta() : -999;
  pMuon->staPhi = (mu->HasStandaloneTrk()) ? mu->StandaloneTrk()->Phi() : -999;
  pMuon->trkID  = (mu->HasTrackerTrk()) ? mu->TrackerTrk()->GetUniqueID() : 0;
  
  pMuon->pfPx=0;
  pMuon->pfPy=0;
  assert(fPFCandidates);
  for(UInt_t i=0; i<fPFCandidates->GetEntries(); ++i) {    
    const PFCandidate *pfcand = fPFCandidates->At(i);
    if(mu->HasTrackerTrk() && mu->TrackerTrk() == pfcand->TrackerTrk()) {
      pMuon->pfPx = pfcand->Px();
      pMuon->pfPy = pfcand->Py();
      break;
    }
  }
}

void NtuplerMod::FillMuon(const Track *mu)
{
  assert(mu);

  TClonesArray &rMuonArr = *fMuonArr;
  assert(rMuonArr.GetEntries() < rMuonArr.GetSize());
  const Int_t index = rMuonArr.GetEntries();  
  new(rMuonArr[index]) TMuon();
  TMuon *pMuon = (TMuon*)rMuonArr[index];
  
  // Use tracker track when available
  const Track *muTrk=mu;
  assert(muTrk);                  
  
  pMuon->pt           = muTrk->Pt();
  pMuon->ptErr        = muTrk->PtErr();
  pMuon->eta          = muTrk->Eta();
  pMuon->phi          = muTrk->Phi();
  pMuon->trkIso03     = 0;
  pMuon->emIso03      = 0;
  pMuon->hadIso03     = 0;
  pMuon->pfChIso03    = computePFIso(mu, 0.0, 0.3, 0.0, 4);
  pMuon->pfNeuIso03   = computePFIso(mu, 0.5, 0.3, 0.0, 2);
  pMuon->pfGamIso03   = computePFIso(mu, 0.5, 0.3, 0.0, 3);
  pMuon->puIso03      = computePUIso(mu, 0.0, 0.3, 0);
  pMuon->pfChIso04    = computePFIso(mu, 0.0, 0.4, 0.0, 4);
  pMuon->pfNeuIso04   = computePFIso(mu, 0.5, 0.4, 0.0, 2);
  pMuon->pfGamIso04   = computePFIso(mu, 0.5, 0.4, 0.0, 3);
  pMuon->puIso04      = computePUIso(mu, 0.0, 0.4, 0); 
  pMuon->d0           = muTrk->D0Corrected(*fVertex);
  pMuon->dz           = muTrk->DzCorrected(*fVertex);
  pMuon->tkNchi2      = muTrk->RChi2();
  pMuon->muNchi2      = muTrk->RChi2();  
  pMuon->trkKink      = 0;
  pMuon->gblKink      = 0;
  pMuon->mva          = 0;
  pMuon->q            = muTrk->Charge();
  pMuon->nValidHits   = 0;  
  pMuon->qualityBits  = 0;            
  pMuon->typeBits     = 0;  
  pMuon->nTkHits      = muTrk->NHits();
  pMuon->nPixHits     = muTrk->NPixelHits();
  pMuon->nSeg         = 0;
  pMuon->nMatch       = 0;
  pMuon->hltMatchBits = MatchHLT(muTrk->Eta(),muTrk->Phi());  
  pMuon->staPt        = -1;  
  pMuon->staEta       = -999;
  pMuon->staPhi       = -999;
  pMuon->trkID        = muTrk->GetUniqueID();
  
  pMuon->nTkLayers=0;
  pMuon->nPixLayers=0;
  if(mu->Hit(Track::PXB1)) { pMuon->nTkLayers++; pMuon->nPixLayers++; } 
  if(mu->Hit(Track::PXB2)) { pMuon->nTkLayers++; pMuon->nPixLayers++; }
  if(mu->Hit(Track::PXB3)) { pMuon->nTkLayers++; pMuon->nPixLayers++; }
  if(mu->Hit(Track::PXF1)) { pMuon->nTkLayers++; pMuon->nPixLayers++; }
  if(mu->Hit(Track::PXF2)) { pMuon->nTkLayers++; pMuon->nPixLayers++; }
  if(mu->Hit(Track::TIB1) || mu->Hit(Track::TIB1S)) pMuon->nTkLayers++;
  if(mu->Hit(Track::TIB2) || mu->Hit(Track::TIB2S)) pMuon->nTkLayers++;
  if(mu->Hit(Track::TIB3)) pMuon->nTkLayers++;
  if(mu->Hit(Track::TIB4)) pMuon->nTkLayers++;
  if(mu->Hit(Track::TID1) || mu->Hit(Track::TID1S)) pMuon->nTkLayers++;
  if(mu->Hit(Track::TID2) || mu->Hit(Track::TID2S)) pMuon->nTkLayers++;
  if(mu->Hit(Track::TID3) || mu->Hit(Track::TID3S)) pMuon->nTkLayers++;
  if(mu->Hit(Track::TOB1) || mu->Hit(Track::TOB1S)) pMuon->nTkLayers++;
  if(mu->Hit(Track::TOB2) || mu->Hit(Track::TOB2S)) pMuon->nTkLayers++;
  if(mu->Hit(Track::TOB3)) pMuon->nTkLayers++;
  if(mu->Hit(Track::TOB4)) pMuon->nTkLayers++;
  if(mu->Hit(Track::TOB5)) pMuon->nTkLayers++;
  if(mu->Hit(Track::TOB6)) pMuon->nTkLayers++;
  if(mu->Hit(Track::TEC1) || mu->Hit(Track::TEC1S)) pMuon->nTkLayers++;
  if(mu->Hit(Track::TEC2) || mu->Hit(Track::TEC2S)) pMuon->nTkLayers++;
  if(mu->Hit(Track::TEC3) || mu->Hit(Track::TEC3S)) pMuon->nTkLayers++;
  if(mu->Hit(Track::TEC4) || mu->Hit(Track::TEC4S)) pMuon->nTkLayers++;
  if(mu->Hit(Track::TEC5) || mu->Hit(Track::TEC5S)) pMuon->nTkLayers++;
  if(mu->Hit(Track::TEC6) || mu->Hit(Track::TEC6S)) pMuon->nTkLayers++;
  if(mu->Hit(Track::TEC7) || mu->Hit(Track::TEC7S)) pMuon->nTkLayers++;
  if(mu->Hit(Track::TEC8) || mu->Hit(Track::TEC8S)) pMuon->nTkLayers++;
  if(mu->Hit(Track::TEC9) || mu->Hit(Track::TEC9S)) pMuon->nTkLayers++;
  
  pMuon->pfPx=0;
  pMuon->pfPy=0;
  assert(fPFCandidates);
  for(UInt_t i=0; i<fPFCandidates->GetEntries(); ++i) {    
    const PFCandidate *pfcand = fPFCandidates->At(i);
    if(muTrk == pfcand->TrackerTrk()) {
      pMuon->pfPx = pfcand->Px();
      pMuon->pfPy = pfcand->Py();
      break;
    }
  }
}

//--------------------------------------------------------------------------------------------------
void NtuplerMod::FillElectron(const Electron *ele)
{
  assert(ele);
 
  ElectronTools eleTools;
  
  TClonesArray &rElectronArr = *fElectronArr;
  assert(rElectronArr.GetEntries() < rElectronArr.GetSize());
  const Int_t index = rElectronArr.GetEntries();  
  new(rElectronArr[index]) TElectron();
  TElectron *pElectron = (TElectron*)rElectronArr[index];
                    
  pElectron->pt              = ele->Pt();
  pElectron->eta             = ele->Eta();
  pElectron->phi             = ele->Phi();
  pElectron->trkIso03        = ele->TrackIsolationDr03();
  pElectron->emIso03         = ele->EcalRecHitIsoDr03();
  pElectron->hadIso03        = ele->HcalTowerSumEtDr03();
  pElectron->pfChIso03       = computePFIso(ele, 0, 0.3, 0.015, 4);
  pElectron->pfNeuIso03      = computePFIso(ele, 0, 0.3, 0,     2);
  pElectron->pfGamIso03      = computePFIso(ele, 0, 0.3, 0.08,  3);
  pElectron->puIso03         = computePUIso(ele, 0, 0.3, 0);
  pElectron->pfChIso04       = computePFIso(ele, 0, 0.4, 0.015, 4);
  pElectron->pfNeuIso04      = computePFIso(ele, 0, 0.4, 0,     2);
  pElectron->pfGamIso04      = computePFIso(ele, 0, 0.4, 0.08,  3);
  pElectron->puIso04         = computePUIso(ele, 0, 0.4, 0);
  pElectron->d0              = ele->BestTrk()->D0Corrected(*fVertex);
  pElectron->dz              = ele->BestTrk()->DzCorrected(*fVertex);  
  pElectron->scEt            = (ele->SCluster()->Energy())*(ele->Pt())/(ele->P());
  pElectron->scEta           = ele->SCluster()->Eta();
  pElectron->scPhi           = ele->SCluster()->Phi();
  pElectron->ecalE           = ele->EcalEnergy();
  pElectron->HoverE          = ele->HadronicOverEm();
  pElectron->EoverP          = ele->ESuperClusterOverP();
  pElectron->fBrem           = ele->FBrem();
  pElectron->deltaEtaIn      = ele->DeltaEtaSuperClusterTrackAtVtx();
  pElectron->deltaPhiIn      = ele->DeltaPhiSuperClusterTrackAtVtx();
  pElectron->sigiEtaiEta     = ele->CoviEtaiEta();
  pElectron->nExpHitsInner   = ele->BestTrk()->NExpectedHitsInner();
  pElectron->partnerDeltaCot = ele->ConvPartnerDCotTheta();
  pElectron->partnerDist     = ele->ConvPartnerDist();
  pElectron->q               = ele->Charge(); 
  
  pElectron->hltMatchBits    = MatchHLT(ele->SCluster()->Eta(),ele->SCluster()->Phi());
  pElectron->scID            = ele->SCluster()->GetUniqueID();
  pElectron->trkID           = (ele->HasTrackerTrk()) ? ele->TrackerTrk()->GetUniqueID() : 0;
  pElectron->isConv          = IsConversion(ele);
  
  pElectron->typeBits = 0;
  if(ele->IsEcalDriven())    { pElectron->typeBits |= kEcalDriven; }
  if(ele->IsTrackerDriven()) { pElectron->typeBits |= kTrackerDriven; }
  
  pElectron->pfPx=0;
  pElectron->pfPy=0;
  assert(fPFCandidates);
  for(UInt_t i=0; i<fPFCandidates->GetEntries(); ++i) {
    const PFCandidate *pfcand = fPFCandidates->At(i);
 
    if( (pfcand->HasTrackerTrk() && ele->TrackerTrk() == pfcand->TrackerTrk()) ||
        (pfcand->HasGsfTrk() && ele->GsfTrk() == pfcand->GsfTrk()) ) {
          pElectron->pfPx = pfcand->Px();
          pElectron->pfPy = pfcand->Py();
          break;
    }	 
  }  
  
  pElectron->mva = fEleMVA->MVAValue(ele, fVertex);
}

//--------------------------------------------------------------------------------------------------
void NtuplerMod::FillJet(const PFJet *jet)
{
  TClonesArray &rPFJetArr = *fPFJetArr;
  assert(rPFJetArr.GetEntries() < rPFJetArr.GetSize());
  const Int_t index = rPFJetArr.GetEntries();  
  new(rPFJetArr[index]) TJet();
  TJet *pPFJet = (TJet*)rPFJetArr[index]; 
  
  const FourVectorM rawMom = jet->RawMom();
  fJetCorrector->setJetEta(rawMom.Eta());
  fJetCorrector->setJetPt(rawMom.Pt());
  fJetCorrector->setJetPhi(rawMom.Phi());
  fJetCorrector->setJetE(rawMom.E());
  fJetCorrector->setRho(fPUEnergyDensity->At(0)->RhoHighEta());
  fJetCorrector->setJetA(jet->JetArea());
  fJetCorrector->setJetEMF(-99.0);

  pPFJet->pt    = (rawMom.Pt())*(fJetCorrector->getCorrection());
  pPFJet->eta   = rawMom.Eta();
  pPFJet->phi   = rawMom.Phi();
  pPFJet->mass  = jet->Mass();
  pPFJet->rawPt = rawMom.Pt();
  pPFJet->area  = jet->JetArea();
  
  pPFJet->tche = jet->TrackCountingHighEffBJetTagsDisc();
  pPFJet->tchp = jet->TrackCountingHighPurBJetTagsDisc();
  
  pPFJet->mcFlavor = jet->MatchedMCFlavor();
  
  pPFJet->hltMatchBits = MatchHLT(jet->Eta(),jet->Phi());
  
  Double_t numer=0, denom=0;
  for(UInt_t ipf=0; ipf<jet->NPFCands(); ipf++) {
    const PFCandidate *pfcand = jet->PFCand(ipf);
    if(!pfcand->HasTrk()) continue;
    numer += (pfcand->Pt())*(pfcand->Pt())*(pfcand->BestTrk()->DzCorrected(*fVertex));
    denom += (pfcand->Pt())*(pfcand->Pt());
  }
  pPFJet->dz = (denom>0) ? numer/denom : 999;
}

      
//--------------------------------------------------------------------------------------------------
void NtuplerMod::FillPhoton(const Photon *pho)
{
  TClonesArray &rPhotonArr = *fPhotonArr;
  assert(rPhotonArr.GetEntries() < rPhotonArr.GetSize());
  const Int_t index = rPhotonArr.GetEntries();  
  new(rPhotonArr[index]) TPhoton();
  TPhoton *pPhoton = (TPhoton*)rPhotonArr[index];
  
  pPhoton->pt           = pho->Pt(); 
  pPhoton->eta          = pho->Eta();
  pPhoton->phi          = pho->Phi();
  pPhoton->scEt         = pho->SCluster()->Et(); 
  pPhoton->scEta        = pho->SCluster()->Eta();
  pPhoton->scPhi        = pho->SCluster()->Phi();
  pPhoton->trkIso04     = pho->HollowConeTrkIsoDr04();
  pPhoton->trkIso04NoPV = IsolationTools::TrackIsolationNoPV(pho, fBeamSpot->At(0), 0.4, 0.04, 0.0, 0.015, 0.1, TrackQuality::highPurity, fTracks);
  pPhoton->emIso04      = pho->EcalRecHitIsoDr04();
  pPhoton->hadIso04     = pho->HcalTowerSumEtDr04();
  pPhoton->HoverE       = pho->HadOverEm();
  pPhoton->R9           = pho->R9();
  pPhoton->sigiEtaiEta  = pho->CoviEtaiEta();
  pPhoton->sigiPhiiPhi  = sqrt(pho->SCluster()->Seed()->CoviPhiiPhi());
  pPhoton->hltMatchBits = MatchHLT(pho->SCluster()->Eta(),pho->SCluster()->Phi());
  pPhoton->scID         = pho->SCluster()->GetUniqueID();
  pPhoton->hasPixelSeed = pho->HasPixelSeed();
  
  pPhoton->passEleVetoConvRec = PhotonTools::PassElectronVetoConvRecovery(pho,fElectrons,fConversions,fBeamSpot->At(0));
  pPhoton->passConvId         = PhotonTools::PassConversionId(pho,PhotonTools::MatchedConversion(pho,fConversions,fBeamSpot->At(0)));
  
  Bool_t passFilter=kTRUE;
  if((pho->SCluster()->Seed()->Energy() > 5.0) && 
     (pho->SCluster()->Seed()->EMax() / pho->SCluster()->Seed()->E3x3() > 0.95))
     passFilter=kFALSE;
  
  pPhoton->passSpikeFilter = passFilter;
}

//--------------------------------------------------------------------------------------------------
void NtuplerMod::FillPV(const Vertex *pv) 
{
  TClonesArray &rPVArr = *fPVArr;
  assert(rPVArr.GetEntries() < rPVArr.GetSize());
  const Int_t index = rPVArr.GetEntries();  
  new(rPVArr[index]) TVertex();
  TVertex *pVertex = (TVertex*)rPVArr[index];
  
  pVertex->nTracksFit = pv->NTracksFit();
  pVertex->ndof       = pv->Ndof();      
  pVertex->chi2       = pv->Chi2();  
  pVertex->x          = pv->X();
  pVertex->y          = pv->Y();
  pVertex->z          = pv->Z();
  
  pVertex->sumPt=0;
  for(UInt_t itrk=0; itrk<pv->NTracks(); itrk++)
    pVertex->sumPt += pv->Trk(itrk)->Pt();				   
}
      
//--------------------------------------------------------------------------------------------------
ULong64_t NtuplerMod::MatchHLT(const Double_t eta, const Double_t phi)
{
  ULong64_t bits = 0;
  
  const Double_t hltMatchR = 0.2;
  
  if(HasHLTInfo()) {
    const TriggerTable *hltTable = GetHLTTable();
    assert(hltTable);
    for(UInt_t itrig=0; itrig<fTriggerNamesv.size(); itrig++) {
      const TriggerName *trigname = hltTable->Get(fTriggerNamesv[itrig].Data());
      if(!trigname) continue;
      
      const TList *list = GetHLTObjectsTable()->GetList(fTriggerNamesv[itrig].Data());
      if(!list) continue;
      TIter iter(list->MakeIterator());
      const TriggerObject *to = dynamic_cast<const TriggerObject*>(iter.Next()); 
   
      while(to) {             
        if(to->IsHLT()) {          
	  
	  if(fTriggerObjNames1v[itrig].Length()>0 && fTriggerObjNames1v[itrig].CompareTo(to->ModuleName())==0) {
	    Bool_t match = kTRUE;
	    if(to->Pt() < fTriggerObjMinPt1v[itrig])                       match=kFALSE;  // minimum pT threshold on trigger object
	    if(MathUtils::DeltaR(phi,eta,to->Phi(),to->Eta()) > hltMatchR) match=kFALSE;  // eta-phi matching
	    if(match) bits |= fTriggerObjIds1v[itrig];
	  }
	  
	  if(fTriggerObjNames2v[itrig].Length()>0 && fTriggerObjNames2v[itrig].CompareTo(to->ModuleName())==0) {
	    Bool_t match = kTRUE;
	    if(to->Pt() < fTriggerObjMinPt2v[itrig])                       match=kFALSE;  // minimum pT threshold on trigger object
	    if(MathUtils::DeltaR(phi,eta,to->Phi(),to->Eta()) > hltMatchR) match=kFALSE;  // eta-phi matching
	    if(match) bits |= fTriggerObjIds2v[itrig];
	  }
	  
	  if(fTriggerObjNames1v[itrig].Length()==0 && fTriggerObjNames2v[itrig].Length()==0) {
	    Bool_t match = kTRUE;
	    if(to->Pt() < fTriggerObjMinPt1v[itrig])                       match=kFALSE;  // minimum pT threshold on trigger object
	    if(MathUtils::DeltaR(phi,eta,to->Phi(),to->Eta()) > hltMatchR) match=kFALSE;  // eta-phi matching
	    if(match) bits |= fTriggerObjIds1v[itrig];
          }
        } 
        to = dynamic_cast<const TriggerObject*>(iter.Next());
      }    
    }
  }
  
  return bits;
}

ULong64_t NtuplerMod::MatchHLT(const Double_t pt, const Double_t eta, const Double_t phi)
{
  ULong64_t bits = 0;
  
  const Double_t hltMatchR = 0.2;
  const Double_t hltMatchPtFrac = 1;
  
  if(HasHLTInfo()) {
    const TriggerTable *hltTable = GetHLTTable();
    assert(hltTable);
    for(UInt_t itrig=0; itrig<fTriggerNamesv.size(); itrig++) {
      const TriggerName *trigname = hltTable->Get(fTriggerNamesv[itrig].Data());
      if(!trigname) continue;
      
      const TList *list = GetHLTObjectsTable()->GetList(fTriggerNamesv[itrig].Data());
      if(!list) continue;
      TIter iter(list->MakeIterator());
      const TriggerObject *to = dynamic_cast<const TriggerObject*>(iter.Next()); 
    
      while(to) {         
        if(to->IsHLT()) {
	  
	  if(fTriggerObjNames1v[itrig].Length()>0 && fTriggerObjNames1v[itrig].CompareTo(to->ModuleName())==0) {
	    Bool_t match = kTRUE;
	    if(to->Pt() < fTriggerObjMinPt1v[itrig])                       match=kFALSE;  // minimum pT threshold on trigger object
	    if(MathUtils::DeltaR(phi,eta,to->Phi(),to->Eta()) > hltMatchR) match=kFALSE;  // eta-phi matching
	    if(fabs(pt - to->Pt())>hltMatchPtFrac*(to->Pt()))              match=kFALSE;  // pT matching
	    if(match) bits |= fTriggerObjIds1v[itrig];
	  }
	  
	  if(fTriggerObjNames2v[itrig].Length()>0 && fTriggerObjNames2v[itrig].CompareTo(to->ModuleName())==0) {
	    Bool_t match = kTRUE;
	    if(to->Pt() < fTriggerObjMinPt2v[itrig])                       match=kFALSE;  // minimum pT threshold on trigger object
	    if(MathUtils::DeltaR(phi,eta,to->Phi(),to->Eta()) > hltMatchR) match=kFALSE;  // eta-phi matching
	    if(fabs(pt - to->Pt())>hltMatchPtFrac*(to->Pt()))              match=kFALSE;  // pT matching
	    if(match) bits |= fTriggerObjIds2v[itrig];
	  }
	  
	  if(fTriggerObjNames1v[itrig].Length()==0 && fTriggerObjNames2v[itrig].Length()==0) {
	    Bool_t match = kTRUE;
	    if(to->Pt() < fTriggerObjMinPt1v[itrig])                       match=kFALSE;  // minimum pT threshold on trigger object
	    if(MathUtils::DeltaR(phi,eta,to->Phi(),to->Eta()) > hltMatchR) match=kFALSE;  // eta-phi matching
	    if(fabs(pt - to->Pt())>hltMatchPtFrac*(to->Pt()))              match=kFALSE;  // pT matching
	    if(match) bits |= fTriggerObjIds1v[itrig];
          }
        }
        to = dynamic_cast<const TriggerObject*>(iter.Next());
      }    
    }
  }
  
  return bits;
}

//--------------------------------------------------------------------------------------------------
Bool_t NtuplerMod::IsConversion(const Electron *ele) 
{
  Bool_t isGoodConversion = kFALSE;
  
  const UInt_t   nWrongHitsMax = 0;
  const Double_t probMin       = 1e-6;
  const Double_t lxyMin        = 2.0;
  const Bool_t   matchCkf      = kTRUE;
  const Bool_t   requireArbitratedMerged = kFALSE;
  
  for (UInt_t ifc=0; ifc<fConversions->GetEntries(); ifc++) {
    Bool_t ConversionMatchFound = kFALSE;
    for (UInt_t d=0; d<fConversions->At(ifc)->NDaughters(); d++) {
      const Track *trk = dynamic_cast<const ChargedParticle*>
        (fConversions->At(ifc)->Daughter(d))->Trk();
      if(ele->GsfTrk() == trk || (matchCkf && ele->TrackerTrk()==trk)) {
        ConversionMatchFound = kTRUE;
        break;
      }
    }

    // if match between the e-track and one of the conversion legs
    if (ConversionMatchFound == kTRUE){
      isGoodConversion =  (fConversions->At(ifc)->Prob() > probMin) &&
        (!requireArbitratedMerged || fConversions->At(ifc)->Quality().Quality(ConversionQuality::arbitratedMerged)) &&
        (fConversions->At(ifc)->LxyCorrected(fVertex) > lxyMin);

      if (isGoodConversion == kTRUE) {
        for (UInt_t d=0; d<fConversions->At(ifc)->NDaughters(); d++) {
          const Track *trk = dynamic_cast<const ChargedParticle*>
            (fConversions->At(ifc)->Daughter(d))->Trk();
          if (trk) {
            const StableData *sd = dynamic_cast<const StableData*>
              (fConversions->At(ifc)->DaughterDat(d));
            if (sd->NWrongHits() > nWrongHitsMax)
              isGoodConversion = kFALSE;
          } else {
            isGoodConversion = kFALSE;
          }
        }
      }
    }

    if(isGoodConversion == kTRUE) break;
    
  } // loop over all fConversions 

  return isGoodConversion;
}

//--------------------------------------------------------------------------------------------------
Float_t NtuplerMod::computePFIso(const Electron *ele, const Double_t ptMin,
                                 const Double_t extRadius, const Double_t intRadius,
                                 const Int_t isoType)
{
  assert(ele);

  Float_t iso = 0.0;
  for(UInt_t i=0; i<fPFNoPileUp->GetEntries(); i++)
  {
    const PFCandidate *pf = fPFNoPileUp->At(i);
    assert(pf);

    PFCandidate::EPFType pfType = pf->PFType();

    if((isoType == 1 && (pfType == PFCandidate::eHadron   ||
                         pfType == PFCandidate::eElectron ||
                         pfType == PFCandidate::eMuon))        ||
       (isoType == 2 && pfType == PFCandidate::eNeutralHadron) ||
       (isoType == 3 && pfType == PFCandidate::eGamma)         ||
       (isoType == 4 && pfType == PFCandidate::eHadron))
    {
      if(pf->Pt() >= ptMin &&
         !(pf->HasTrackerTrk() && ele->HasTrackerTrk() && pf->TrackerTrk() == ele->TrackerTrk()) &&
	 !(pf->HasGsfTrk()     && ele->HasGsfTrk()     && pf->GsfTrk()     == ele->GsfTrk()))
      {
        // Add p_T to running sum if PFCandidate is close enough
        Double_t dr = MathUtils::DeltaR(ele->Mom(), pf->Mom());
        if(ele->IsEB()) {// && ele->Mva()>=-0.1) {
	  if(dr < extRadius) iso += pf->Pt();
	} else {
	  if(dr < extRadius && dr >= intRadius) iso += pf->Pt();
	}
      }
    }
  }

  return iso;
}

Float_t NtuplerMod::computePFIso(const Track *track, const Double_t ptMin, 
				 const Double_t extRadius, const Double_t intRadius,
				 const Int_t isoType)
{
  assert(track);

  Float_t iso = 0.0;
  for(UInt_t i=0; i<fPFNoPileUp->GetEntries(); i++)
  {
    const PFCandidate *pf = fPFNoPileUp->At(i);
    assert(pf);

    PFCandidate::EPFType pfType = pf->PFType();

    if((isoType == 1 && (pfType == PFCandidate::eHadron   ||
                         pfType == PFCandidate::eElectron ||
                         pfType == PFCandidate::eMuon))        ||
       (isoType == 2 && pfType == PFCandidate::eNeutralHadron) ||
       (isoType == 3 && pfType == PFCandidate::eGamma)         ||
       (isoType == 4 && pfType == PFCandidate::eHadron))
    {
      if(pf->Pt() >= ptMin &&
         !(pf->TrackerTrk() && pf->TrackerTrk() == track))
      {
        // Add p_T to running sum if PFCandidate is close enough
        Double_t dr = MathUtils::DeltaR(track->Mom(), pf->Mom());
        if(dr < extRadius && dr >= intRadius) iso += pf->Pt();
      }
    }
  }

  return iso;
}

//--------------------------------------------------------------------------------------------------
Float_t NtuplerMod::computePUIso(const Electron *ele, const Double_t ptMin,
                                 const Double_t extRadius, const Double_t intRadius)
{
  assert(fPFPileUp);

  Float_t iso = 0.0;

  for(UInt_t i=0; i<fPFPileUp->GetEntries(); i++) {
    const PFCandidate *pf = fPFPileUp->At(i);
    assert(pf);

    if(pf->Pt() >= ptMin) {
      // Add p_T to running sum if PFCandidate is within isolation cone
      Double_t dr = MathUtils::DeltaR(ele->Mom(), pf->Mom());
      if(dr < extRadius && dr >= intRadius) iso += pf->Pt();
    }
  }

  return iso;
}

Float_t NtuplerMod::computePUIso(const Track *track, const Double_t ptMin,
                                 const Double_t extRadius, const Double_t intRadius)
{
  assert(fPFPileUp);

  Float_t iso = 0.0;

  for(UInt_t i=0; i<fPFPileUp->GetEntries(); i++) {
    const PFCandidate *pf = fPFPileUp->At(i);
    assert(pf);

    if(pf->Pt() >= ptMin) {
      // Add p_T to running sum if PFCandidate is within isolation cone
      Double_t dr = MathUtils::DeltaR(track->Mom(), pf->Mom());
      if(dr < extRadius && dr >= intRadius) iso += pf->Pt();
    }
  }

  return iso;
}

//--------------------------------------------------------------------------------------------------
void NtuplerMod::separatePileUp(const Bool_t checkClosestZVertex)
{
  assert(fPFPileUp);
  assert(fPFNoPileUp);

  fPFPileUp->Reset();
  fPFNoPileUp->Reset();

  for(UInt_t i=0; i<fPFCandidates->GetEntries(); i++) {
    const PFCandidate *pf = fPFCandidates->At(i);
    assert(pf);

    if(pf->PFType() == PFCandidate::eHadron) {
      if(pf->HasTrackerTrk() &&
         fVertex->HasTrack(pf->TrackerTrk()) &&
         fVertex->TrackWeight(pf->TrackerTrk()) > 0)
      {
        fPFNoPileUp->Add(pf);
      
      } else {
        Bool_t vertexFound = kFALSE;
        const Vertex *closestVtx = 0;
        Double_t dzmin = 10000;

        for(UInt_t j=0; j<fPrimVerts->GetEntries(); j++) {
          const Vertex *vtx = fPrimVerts->At(j);
          assert(vtx);

          if(pf->HasTrackerTrk() &&
             vtx->HasTrack(pf->TrackerTrk()) &&
             vtx->TrackWeight(pf->TrackerTrk()) > 0)
          {
            vertexFound = kTRUE;
            closestVtx = vtx;
            break;
          }

          Double_t dz = fabs(pf->SourceVertex().Z() - vtx->Z());
          if(dz < dzmin) {
            closestVtx = vtx;
            dzmin = dz;
          }
        }

        if(checkClosestZVertex) {
          // Fallback: if track is not associated with any vertex,
          // associate it with the vertex closest in z
          if(vertexFound || closestVtx != fVertex)
            fPFPileUp->Add(pf);
          else
            fPFNoPileUp->Add(pf);
        
	} else {
          if(vertexFound && closestVtx != fVertex)
            fPFPileUp->Add(pf);
          else
            fPFNoPileUp->Add(pf); // Ridiculous but that's how it is
        }
      }
    } else {
      fPFNoPileUp->Add(pf);
    }
  }
}

//--------------------------------------------------------------------------------------------------
void NtuplerMod::FillGenH() 
{
  assert(fParticles);
  assert(fMCEvtInfo);
 
  const MCParticle *boson=0, *dau1=0, *dau2=0;
  
  Int_t id1=0, id2=0;
  
  for(UInt_t i=0; i<fParticles->GetEntries(); ++i) {
    const MCParticle *p = fParticles->At(i);
      
    if( (p->PdgId() == 25) && (p->Status() == 3) ) {
      boson = p;
      
      // loop through daughters and look for Ws
      for(UInt_t ii=0; ii<boson->NDaughters(); ii++) {
        const MCParticle *bosonv = boson->Daughter(ii); 
	if(abs(bosonv->PdgId())==24 || abs(bosonv->PdgId())==23) {  
	
          // Loop through daughters and look for leptons  
          for(UInt_t j=0; j<bosonv->NDaughters(); j++) {
            const MCParticle *d = bosonv->Daughter(j);
            if((abs(d->PdgId())==11) || (abs(d->PdgId())==13) || (abs(d->PdgId())==15)) {
          
	      // traverse down daughter muon tree
              while(d->HasDaughter(d->PdgId(),kTRUE) && (d->Status()!=1))
                d = d->FindDaughter(d->PdgId(),kTRUE);	  

              if(d->PdgId()== 11) id1 =  EGenType::kElectron;
	      if(d->PdgId()==-11) id2 = -EGenType::kElectron;
	      if(d->PdgId()== 13) id1 =  EGenType::kMuon;
	      if(d->PdgId()==-13) id2 = -EGenType::kMuon;
	      if(d->PdgId()== 15) id1 =  EGenType::kTau;
	      if(d->PdgId()==-15) id2 = -EGenType::kTau;
	      
	      if(abs(d->PdgId())==15) {
	        if(d->HasDaughter(11)) { 
	          d = d->FindDaughter(11); 
	          if(d->PdgId()== 11) id1 =  EGenType::kTauElectron;
	          if(d->PdgId()==-11) id2 = -EGenType::kTauElectron;
	        }
	        if(d->HasDaughter(13)) { 
	          d = d->FindDaughter(13); 
	          if(d->PdgId()== 13) id1 =  EGenType::kTauMuon;
	          if(d->PdgId()==-13) id2 = -EGenType::kTauMuon;
	        }
	      }
	      
	      if(d->PdgId()>0) dau1 = d;
              if(d->PdgId()<0) dau2 = d;
            }
          }	  
	}	
      }      		      
    }                            
  }
  
  assert(boson);
    
  if(!dau1 || !dau2) {
    for(UInt_t j=0; j<boson->NDaughters(); j++) {
      const MCParticle *d = boson->Daughter(j);
      if(d->PdgId()==boson->PdgId()) continue;
      if(d->PdgId()>0) dau1 = d;
      if(d->PdgId()<0) dau2 = d;
    }
    //cout << "Weird: " << boson->PdgId() << " --> " << dau1->PdgId() << " + " << dau2->PdgId() << endl;
  }
    
  assert(dau1);
  assert(dau2);  
    
  // Get the status=3 daughters. Not applicable for HORACE or PHOTOS modes.
  const MCParticle *vdau1=dau1;
  const MCParticle *vdau2=dau2;
  if(fFSRMode==0) {
    while(vdau1->Status()!=3)
      vdau1 = vdau1->FindMother(vdau1->PdgId(),kTRUE);
    
    while(vdau2->Status()!=3)
      vdau2 = vdau2->FindMother(vdau2->PdgId(),kTRUE);
  }
  
  FourVectorM vDilep = dau1->Mom() + dau2->Mom();
  
  fGenInfo.pid_1    = fMCEvtInfo->Id1();
  fGenInfo.pid_2    = fMCEvtInfo->Id2();
  fGenInfo.procId   = fMCEvtInfo->ProcessId();
  fGenInfo.x_1      = fMCEvtInfo->X1();
  fGenInfo.x_2      = fMCEvtInfo->X2();
  fGenInfo.scalePDF = fMCEvtInfo->ScalePdf();
  fGenInfo.weight   = fMCEvtInfo->Weight();
  fGenInfo.vmass    = boson->Mass();
  fGenInfo.vpt      = boson->Pt();
  fGenInfo.vy       = boson->Rapidity();
  fGenInfo.vphi     = boson->Phi();
  fGenInfo.mass     = vDilep.M();
  fGenInfo.pt       = vDilep.Pt(); 
  fGenInfo.y        = vDilep.Rapidity(); 
  fGenInfo.phi      = vDilep.Phi(); 
  fGenInfo.vpt_1    = vdau1->Pt(); 
  fGenInfo.veta_1   = vdau1->Eta(); 
  fGenInfo.vphi_1   = vdau1->Phi();
  fGenInfo.vpt_2    = vdau2->Pt();
  fGenInfo.veta_2   = vdau2->Eta(); 
  fGenInfo.vphi_2   = vdau2->Phi(); 
  fGenInfo.id       = EGenType::kHiggs;
  fGenInfo.pt_1     = dau1->Pt(); 
  fGenInfo.eta_1    = dau1->Eta(); 
  fGenInfo.phi_1    = dau1->Phi();
  fGenInfo.id_1     = id1;
  fGenInfo.pt_2     = dau2->Pt();
  fGenInfo.eta_2    = dau2->Eta(); 
  fGenInfo.phi_2    = dau2->Phi(); 
  fGenInfo.id_2     = id2;
  fGenInfo.npho     = 0;
  fGenInfo.phopt    = 0;
  fGenInfo.phoeta   = 0;
  fGenInfo.phophi   = 0;
  fGenInfo.decx     = boson->DecayVertex().X();
  fGenInfo.decy     = boson->DecayVertex().Y(); 
  fGenInfo.decz     = boson->DecayVertex().Z();
}

//--------------------------------------------------------------------------------------------------
void NtuplerMod::FillGenZ() 
{
  assert(fParticles);
  assert(fMCEvtInfo);
  
  const MCParticle *boson=0, *dau1=0, *dau2=0;
  
  Int_t id1=0, id2=0;
  
  for(UInt_t i=0; i<fParticles->GetEntries(); ++i) {
    const MCParticle *p = fParticles->At(i);
      
    if( (p->PdgId() == 23) && (p->Status() == 3) ) {
      boson = p;

      // Loop through daughters and look for leptons  
      for(UInt_t j=0; j<boson->NDaughters(); j++) {
        const MCParticle *d = boson->Daughter(j);
        if((abs(d->PdgId())==11) || (abs(d->PdgId())==13) || (abs(d->PdgId())==15) ||
	   (abs(d->PdgId())==12) || (abs(d->PdgId())==14) || (abs(d->PdgId())==16)) {
          
          // traverse down daughter muon tree
          while(d->HasDaughter(d->PdgId(),kTRUE) && (d->Status()!=1))
            d = d->FindDaughter(d->PdgId(),kTRUE);	  

          if(d->PdgId()== 11) id1 =  EGenType::kElectron;
	  if(d->PdgId()==-11) id2 = -EGenType::kElectron;
	  if(d->PdgId()== 13) id1 =  EGenType::kMuon;
	  if(d->PdgId()==-13) id2 = -EGenType::kMuon;
	  if(d->PdgId()== 15) id1 =  EGenType::kTau;
	  if(d->PdgId()==-15) id2 = -EGenType::kTau;
	  
	  if(abs(d->PdgId())==15) {
	    if(d->HasDaughter(11)) { 
	      d = d->FindDaughter(11); 
	      if(d->PdgId()== 11) id1 =  EGenType::kTauElectron;
	      if(d->PdgId()==-11) id2 = -EGenType::kTauElectron;
	    }
	    if(d->HasDaughter(13)) { 
	      d = d->FindDaughter(13); 
	      if(d->PdgId()== 13) id1 =  EGenType::kTauMuon;
	      if(d->PdgId()==-13) id2 = -EGenType::kTauMuon;
	    }
	  }
	  
	  if(d->PdgId()>0) dau1 = d;
          if(d->PdgId()<0) dau2 = d;
        }  
      }		      
    }                            
  }
  
  assert(boson);
  assert(dau1);
  assert(dau2);
    
  // Get the status=3 daughters. Not applicable for HORACE or PHOTOS modes.
  const MCParticle *vdau1=dau1;
  const MCParticle *vdau2=dau2;
  if(fFSRMode==0) {
    const MCParticle *tmp=0;
    while(vdau1->Status()!=3) {
      tmp = vdau1->FindMother(vdau1->PdgId(),kTRUE);
      if(!tmp) tmp = vdau1->FindMother(15,kFALSE);
      assert(tmp);
      vdau1=tmp;
    }
    
    tmp=0;
    while(vdau2->Status()!=3) {
      tmp = vdau2->FindMother(vdau2->PdgId(),kTRUE);
      if(!tmp) tmp = vdau2->FindMother(15,kFALSE);
      assert(tmp);
      vdau2=tmp;
    }
    assert(vdau2);
  }
  
  FourVectorM vDilep = dau1->Mom() + dau2->Mom();
  
  fGenInfo.pid_1    = fMCEvtInfo->Id1();
  fGenInfo.pid_2    = fMCEvtInfo->Id2();
  fGenInfo.procId   = fMCEvtInfo->ProcessId();
  fGenInfo.x_1      = fMCEvtInfo->X1();
  fGenInfo.x_2      = fMCEvtInfo->X2();
  fGenInfo.scalePDF = fMCEvtInfo->ScalePdf();
  fGenInfo.weight   = fMCEvtInfo->Weight();
  fGenInfo.vmass    = boson->Mass();
  fGenInfo.vpt      = boson->Pt();
  fGenInfo.vy       = boson->Rapidity();
  fGenInfo.vphi     = boson->Phi();
  fGenInfo.vpt_1    = vdau1->Pt(); 
  fGenInfo.veta_1   = vdau1->Eta(); 
  fGenInfo.vphi_1   = vdau1->Phi();
  fGenInfo.vpt_2    = vdau2->Pt();
  fGenInfo.veta_2   = vdau2->Eta(); 
  fGenInfo.vphi_2   = vdau2->Phi(); 
  fGenInfo.mass     = vDilep.M();
  fGenInfo.pt       = vDilep.Pt(); 
  fGenInfo.y        = vDilep.Rapidity(); 
  fGenInfo.phi      = vDilep.Phi(); 
  fGenInfo.id       = EGenType::kZ;
  fGenInfo.pt_1     = dau1->Pt(); 
  fGenInfo.eta_1    = dau1->Eta(); 
  fGenInfo.phi_1    = dau1->Phi();
  fGenInfo.id_1     = id1;
  fGenInfo.pt_2     = dau2->Pt();
  fGenInfo.eta_2    = dau2->Eta(); 
  fGenInfo.phi_2    = dau2->Phi(); 
  fGenInfo.id_2     = id2;
  fGenInfo.npho     = 0;
  fGenInfo.phopt    = 0;
  fGenInfo.phoeta   = 0;
  fGenInfo.phophi   = 0;
  fGenInfo.decx     = boson->DecayVertex().X();
  fGenInfo.decy     = boson->DecayVertex().Y(); 
  fGenInfo.decz     = boson->DecayVertex().Z();
}

//--------------------------------------------------------------------------------------------------
void NtuplerMod::FillGenW() 
{
  assert(fParticles);
  assert(fMCEvtInfo);
  
  const MCParticle *boson=0, *dau1=0, *dau2=0;
  
  Int_t id1=0, id2=0;
  
  for(UInt_t i=0; i<fParticles->GetEntries(); ++i) {
    const MCParticle *p = fParticles->At(i);
      
    if( (abs(p->PdgId()) == 24) && (p->Status() == 3) ) {
      boson = p;

      // Loop through daughters and look for leptons  
      for(UInt_t j=0; j<boson->NDaughters(); j++) {
        const MCParticle *d = boson->Daughter(j);
        if(d->PdgId() == boson->PdgId()) continue;  
        
	// traverse down daughter lepton tree
        while(d->HasDaughter(d->PdgId(),kTRUE) && (d->Status()!=1))
          d = d->FindDaughter(d->PdgId(),kTRUE);	  

        if(d->PdgId()== 11) id1 =  EGenType::kElectron;
	if(d->PdgId()==-11) id2 = -EGenType::kElectron;
	if(d->PdgId()== 13) id1 =  EGenType::kMuon;
	if(d->PdgId()==-13) id2 = -EGenType::kMuon;
	if(d->PdgId()== 15) id1 =  EGenType::kTau;
	if(d->PdgId()==-15) id2 = -EGenType::kTau;
	
	if(abs(d->PdgId())==15) {
	  if(d->HasDaughter(11)) { 
	    d = d->FindDaughter(11);
	    if(d->PdgId()== 11) id1 =  EGenType::kTauElectron;
	    if(d->PdgId()==-11) id2 = -EGenType::kTauElectron; 
	  }
	  if(d->HasDaughter(13)) { 
	    d = d->FindDaughter(13);
	    if(d->PdgId()== 13) id1 =  EGenType::kTauMuon;
	    if(d->PdgId()==-13) id2 = -EGenType::kTauMuon; 
	  }
	}
	
	if(d->PdgId()>0) dau1 = d;
        if(d->PdgId()<0) dau2 = d; 
      }		      
    }                            
  }
  
  assert(boson);
  if(!dau1 || !dau2) {
    while(boson->HasDaughter(boson->PdgId(),kTRUE))
      boson = boson->FindDaughter(boson->PdgId(),kTRUE);
    
    // Loop through daughters and look for leptons  
    for(UInt_t j=0; j<boson->NDaughters(); j++) {
      const MCParticle *d = boson->Daughter(j);
      if(d->PdgId() == boson->PdgId()) continue;  
      
      // traverse down daughter lepton tree
      while(d->HasDaughter(d->PdgId(),kTRUE) && (d->Status()!=1))
    	d = d->FindDaughter(d->PdgId(),kTRUE);  	

      if(d->PdgId()== 11) id1 =  EGenType::kElectron;
      if(d->PdgId()==-11) id2 = -EGenType::kElectron;
      if(d->PdgId()== 13) id1 =  EGenType::kMuon;
      if(d->PdgId()==-13) id2 = -EGenType::kMuon;
      if(d->PdgId()== 15) id1 =  EGenType::kTau;
      if(d->PdgId()==-15) id2 = -EGenType::kTau;

      if(abs(d->PdgId())==15) {
        if(d->HasDaughter(11)) { 
          d = d->FindDaughter(11);
          if(d->PdgId()== 11) id1 =  EGenType::kTauElectron;
          if(d->PdgId()==-11) id2 = -EGenType::kTauElectron; 
        }
        if(d->HasDaughter(13)) { 
          d = d->FindDaughter(13);
          if(d->PdgId()== 13) id1 =  EGenType::kTauMuon;
          if(d->PdgId()==-13) id2 = -EGenType::kTauMuon; 
        }
      }

      if(d->PdgId()>0) dau1 = d;
      if(d->PdgId()<0) dau2 = d; 
    }
  }
  assert(dau1);
  assert(dau2);
    
  // Get the status=3 daughters. Not applicable for HORACE or PHOTOS modes.
  const MCParticle *vdau1=dau1;
  const MCParticle *vdau2=dau2;
  if(fFSRMode==0) {
    const MCParticle *tmp=0;
    while(vdau1->Status()!=3) {
      tmp = vdau1->FindMother(vdau1->PdgId(),kTRUE);
      if(!tmp) tmp = vdau1->FindMother(15,kFALSE);
      if(!tmp) break;
      vdau1=tmp;
    }
    
    tmp=0;
    while(vdau2->Status()!=3) {
      tmp = vdau2->FindMother(vdau2->PdgId(),kTRUE);
      if(!tmp) tmp = vdau2->FindMother(15,kFALSE);
      if(!tmp) break;
      vdau2=tmp;
    }
  }
  
  FourVectorM vDilep = dau1->Mom() + dau2->Mom();
  
  fGenInfo.pid_1    = fMCEvtInfo->Id1();
  fGenInfo.pid_2    = fMCEvtInfo->Id2();
  fGenInfo.procId   = fMCEvtInfo->ProcessId();
  fGenInfo.x_1      = fMCEvtInfo->X1();
  fGenInfo.x_2      = fMCEvtInfo->X2();
  fGenInfo.scalePDF = fMCEvtInfo->ScalePdf();
  fGenInfo.weight   = fMCEvtInfo->Weight();
  fGenInfo.vmass    = boson->Mass();
  fGenInfo.vpt      = boson->Pt();
  fGenInfo.vy       = boson->Rapidity();
  fGenInfo.vphi     = boson->Phi();
  fGenInfo.vpt_1    = vdau1->Pt(); 
  fGenInfo.veta_1   = vdau1->Eta(); 
  fGenInfo.vphi_1   = vdau1->Phi();
  fGenInfo.vpt_2    = vdau2->Pt();
  fGenInfo.veta_2   = vdau2->Eta(); 
  fGenInfo.vphi_2   = vdau2->Phi(); 
  fGenInfo.mass     = vDilep.M();
  fGenInfo.pt       = vDilep.Pt(); 
  fGenInfo.y        = vDilep.Rapidity(); 
  fGenInfo.phi      = vDilep.Phi(); 
  fGenInfo.id       = boson->PdgId();
  fGenInfo.pt_1     = dau1->Pt(); 
  fGenInfo.eta_1    = dau1->Eta(); 
  fGenInfo.phi_1    = dau1->Phi();
  fGenInfo.id_1     = id1;
  fGenInfo.pt_2     = dau2->Pt();
  fGenInfo.eta_2    = dau2->Eta(); 
  fGenInfo.phi_2    = dau2->Phi(); 
  fGenInfo.id_2     = id2;
  fGenInfo.npho     = 0;
  fGenInfo.phopt    = 0;
  fGenInfo.phoeta   = 0;
  fGenInfo.phophi   = 0;
  fGenInfo.decx     = boson->DecayVertex().X();
  fGenInfo.decy     = boson->DecayVertex().Y(); 
  fGenInfo.decz     = boson->DecayVertex().Z();
}

//--------------------------------------------------------------------------------------------------
void NtuplerMod::FillGenWW() 
{
  assert(fParticles);
  assert(fMCEvtInfo);
  
  const MCParticle *boson=0, *dau1=0, *dau2=0;
  
  Int_t id=0, id1=0, id2=0;
  Int_t nW=0;
  
  for(UInt_t i=0; i<fParticles->GetEntries(); ++i) {
    const MCParticle *p = fParticles->At(i);
    
    if( (abs(p->PdgId()) == 24) && (p->Status() == 3) ) {
      boson = p;
      nW++;
      if(nW==2) id = EGenType::kWW;

      // Loop through daughters and look for leptons  
      for(UInt_t j=0; j<boson->NDaughters(); j++) {
        const MCParticle *d = boson->Daughter(j);
        if(d->PdgId() == boson->PdgId()) continue;  
        
	// traverse down daughter lepton tree
        while(d->HasDaughter(d->PdgId(),kTRUE) && (d->Status()!=1))
          d = d->FindDaughter(d->PdgId(),kTRUE);	  

        // ignore neutrinos
	if(abs(d->PdgId())==12) continue;
	if(abs(d->PdgId())==14) continue;
	if(abs(d->PdgId())==16) continue;
	
	if(d->PdgId()== 11) id1 =  EGenType::kElectron;
	if(d->PdgId()==-11) id2 = -EGenType::kElectron;
	if(d->PdgId()== 13) id1 =  EGenType::kMuon;
	if(d->PdgId()==-13) id2 = -EGenType::kMuon;
	if(d->PdgId()== 15) id1 =  EGenType::kTau;
	if(d->PdgId()==-15) id2 = -EGenType::kTau;
	
	if(abs(d->PdgId())==15) {
	  if(d->HasDaughter(11)) { 
	    d = d->FindDaughter(11); 
	    if(d->PdgId()== 11) id1 =  EGenType::kTauElectron;
	    if(d->PdgId()==-11) id2 = -EGenType::kTauElectron; 
	  }
	  if(d->HasDaughter(13)) { 
	    d = d->FindDaughter(13); 
	    if(d->PdgId()== 13) id1 =  EGenType::kTauMuon;
	    if(d->PdgId()==-13) id2 = -EGenType::kTauMuon;
	  }
	}
	
	if(d->PdgId()>0) dau1 = d;
        if(d->PdgId()<0) dau2 = d; 
      }		      
    }                            
  }
  
  FourVectorM vDilep;
  const MCParticle *vdau1=dau1;
  const MCParticle *vdau2=dau2;
  if(dau1 && dau2) {
    vDilep = dau1->Mom() + dau2->Mom();
    
    // Get the status=3 daughters. Not applicable for HORACE or PHOTOS modes.    
    if(fFSRMode==0) {
      while(vdau1->Status()!=3)
        vdau1 = vdau1->FindMother(vdau1->PdgId(),kTRUE);
    
      while(vdau2->Status()!=3)
        vdau2 = vdau2->FindMother(vdau2->PdgId(),kTRUE);
    }    
  }
  
  fGenInfo.pid_1    = fMCEvtInfo->Id1();
  fGenInfo.pid_2    = fMCEvtInfo->Id2();
  fGenInfo.procId   = fMCEvtInfo->ProcessId();
  fGenInfo.x_1      = fMCEvtInfo->X1();
  fGenInfo.x_2      = fMCEvtInfo->X2();
  fGenInfo.scalePDF = fMCEvtInfo->ScalePdf();
  fGenInfo.weight   = fMCEvtInfo->Weight();
  fGenInfo.vmass    = 0;
  fGenInfo.vpt      = 0;
  fGenInfo.vy       = 0;
  fGenInfo.vphi     = 0;
  
  fGenInfo.vpt_1   = vdau1 ? vdau1->Pt()  : 0; 
  fGenInfo.veta_1  = vdau1 ? vdau1->Eta() : 0; 
  fGenInfo.vphi_1  = vdau1 ? vdau1->Phi() : 0;
  
  fGenInfo.vpt_2   = vdau2 ? vdau2->Pt()  : 0;
  fGenInfo.veta_2  = vdau2 ? vdau2->Eta() : 0; 
  fGenInfo.vphi_2  = vdau2 ? vdau2->Phi() : 0; 
      
  fGenInfo.mass   = (dau1 && dau2) ? vDilep.M()        : 0;
  fGenInfo.pt     = (dau1 && dau2) ? vDilep.Pt()       : 0; 
  fGenInfo.y      = (dau1 && dau2) ? vDilep.Rapidity() : 0; 
  fGenInfo.phi    = (dau1 && dau2) ? vDilep.Phi()      : 0; 
  fGenInfo.id     = id;
  
  fGenInfo.pt_1   = dau1 ? dau1->Pt()  : 0; 
  fGenInfo.eta_1  = dau1 ? dau1->Eta() : 0; 
  fGenInfo.phi_1  = dau1 ? dau1->Phi() : 0;
  fGenInfo.id_1   = id1;
  
  fGenInfo.pt_2   = dau2 ? dau2->Pt()  : 0;
  fGenInfo.eta_2  = dau2 ? dau2->Eta() : 0; 
  fGenInfo.phi_2  = dau2 ? dau2->Phi() : 0; 
  fGenInfo.id_2   = id2;
  
  fGenInfo.npho     = 0;
  fGenInfo.phopt    = 0;
  fGenInfo.phoeta   = 0;
  fGenInfo.phophi   = 0;
  
  fGenInfo.decx   = boson ? boson->DecayVertex().X() : -999;
  fGenInfo.decy   = boson ? boson->DecayVertex().Y() : -999; 
  fGenInfo.decz   = boson ? boson->DecayVertex().Z() : -999;
}

//--------------------------------------------------------------------------------------------------
void NtuplerMod::FillGenVZ() 
{
  assert(fParticles);
  assert(fMCEvtInfo);
  
  const MCParticle *boson=0, *dau1=0, *dau2=0;
  
  Int_t id1=0, id2=0;
  
  for(UInt_t i=0; i<fParticles->GetEntries(); ++i) {
    const MCParticle *p = fParticles->At(i);
      
    if( (p->PdgId() == 23) && (p->Status() == 3) ) {

      // Loop through daughters and look for leptons  
      for(UInt_t j=0; j<p->NDaughters(); j++) {
        const MCParticle *d = p->Daughter(j);
        if((abs(d->PdgId())==11) || (abs(d->PdgId())==13) || (abs(d->PdgId())==15)) {
          
          // traverse down daughter muon tree
          while(d->HasDaughter(d->PdgId(),kTRUE) && (d->Status()!=1))
            d = d->FindDaughter(d->PdgId(),kTRUE);	  

          if(d->PdgId()== 11) id1 =  EGenType::kElectron;
	  if(d->PdgId()==-11) id2 = -EGenType::kElectron;
	  if(d->PdgId()== 13) id1 =  EGenType::kMuon;
	  if(d->PdgId()==-13) id2 = -EGenType::kMuon;
	  if(d->PdgId()== 15) id1 =  EGenType::kTau;
	  if(d->PdgId()==-15) id2 = -EGenType::kTau;
	  
	  if(abs(d->PdgId())==15) {
	    if(d->HasDaughter(11)) { 
	      d = d->FindDaughter(11); 
	      if(d->PdgId()== 11) id1 =  EGenType::kTauElectron;
	      if(d->PdgId()==-11) id2 = -EGenType::kTauElectron;
	    }
	    if(d->HasDaughter(13)) { 
	      d = d->FindDaughter(13); 
	      if(d->PdgId()== 13) id1 =  EGenType::kTauMuon;
	      if(d->PdgId()==-13) id2 = -EGenType::kTauMuon;
	    }
	  }
	  
	  if(d->PdgId()>0) dau1 = d;
          if(d->PdgId()<0) dau2 = d;
	  if(dau1 && dau2) boson = p;
        }  
      }		      
    }                            
  }
  
  Bool_t hasZll = (boson && dau1 && dau2);
  
  FourVectorM vDilep;
  const MCParticle *vdau1=dau1;
  const MCParticle *vdau2=dau2;
  if(hasZll) {
    vDilep = dau1->Mom() + dau2->Mom();
    
    // Get the status=3 daughters. Not applicable for HORACE or PHOTOS modes.    
    if(fFSRMode==0) {
      while(vdau1->Status()!=3)
        vdau1 = vdau1->FindMother(vdau1->PdgId(),kTRUE);
    
      while(vdau2->Status()!=3)
        vdau2 = vdau2->FindMother(vdau2->PdgId(),kTRUE);
    } 
  }
  
  fGenInfo.pid_1    = fMCEvtInfo->Id1();
  fGenInfo.pid_2    = fMCEvtInfo->Id2();
  fGenInfo.procId   = fMCEvtInfo->ProcessId();
  fGenInfo.x_1      = fMCEvtInfo->X1();
  fGenInfo.x_2      = fMCEvtInfo->X2();
  fGenInfo.scalePDF = fMCEvtInfo->ScalePdf();
  fGenInfo.weight   = fMCEvtInfo->Weight();
  fGenInfo.vmass    = (hasZll) ? boson->Mass() : 0;
  fGenInfo.vpt      = (hasZll) ? boson->Pt() : 0;
  fGenInfo.vy       = (hasZll) ? boson->Rapidity() : 0;
  fGenInfo.vphi     = (hasZll) ? boson->Phi() : 0;
  fGenInfo.vpt_1    = (hasZll) ? vdau1->Pt() : 0; 
  fGenInfo.veta_1   = (hasZll) ? vdau1->Eta() : 0; 
  fGenInfo.vphi_1   = (hasZll) ? vdau1->Phi() : 0;
  fGenInfo.vpt_2    = (hasZll) ? vdau2->Pt() : 0;
  fGenInfo.veta_2   = (hasZll) ? vdau2->Eta() : 0; 
  fGenInfo.vphi_2   = (hasZll) ? vdau2->Phi() : 0;
  fGenInfo.mass     = (hasZll) ? vDilep.M() : 0;
  fGenInfo.pt       = (hasZll) ? vDilep.Pt() : 0; 
  fGenInfo.y        = (hasZll) ? vDilep.Rapidity() : 0; 
  fGenInfo.phi      = (hasZll) ? vDilep.Phi() : 0; 
  fGenInfo.id       = (hasZll) ? EGenType::kVZ : 0;
  fGenInfo.pt_1     = (hasZll) ? dau1->Pt() : 0; 
  fGenInfo.eta_1    = (hasZll) ? dau1->Eta() : 0; 
  fGenInfo.phi_1    = (hasZll) ? dau1->Phi() : 0;
  fGenInfo.id_1     = (hasZll) ? id1 : 0;
  fGenInfo.pt_2     = (hasZll) ? dau2->Pt() : 0;
  fGenInfo.eta_2    = (hasZll) ? dau2->Eta() : 0; 
  fGenInfo.phi_2    = (hasZll) ? dau2->Phi() : 0; 
  fGenInfo.id_2     = (hasZll) ? id2 : 0;
  fGenInfo.npho     = 0;
  fGenInfo.phopt    = 0;
  fGenInfo.phoeta   = 0;
  fGenInfo.phophi   = 0;
  fGenInfo.decx     = (hasZll) ? boson->DecayVertex().X() : 0;
  fGenInfo.decy     = (hasZll) ? boson->DecayVertex().Y() : 0; 
  fGenInfo.decz     = (hasZll) ? boson->DecayVertex().Z() : 0;
}

//--------------------------------------------------------------------------------------------------
void NtuplerMod::FillGenWjets() 
{
  assert(fParticles);
  assert(fMCEvtInfo);
  
  const MCParticle *boson=0, *dau1=0, *dau2=0;
  
  Int_t id1=0, id2=0;
  
  Bool_t isWgamma=kFALSE;
  
  for(UInt_t i=0; i<fParticles->GetEntries(); ++i) {
    const MCParticle *p = fParticles->At(i);
    
    if( (p->PdgId() == 22) && (p->Status() == 1) ) {
       
       // Radiated photon
       if(p->DistinctMother() && (p->DistinctMother()->Status() == 3)
          && (p->DistinctMother()->Is(MCParticle::kEl)  || p->DistinctMother()->Is(MCParticle::kMu) ||
              p->DistinctMother()->Is(MCParticle::kTau) || p->DistinctMother()->Is(MCParticle::kW)) ) {
         isWgamma=kTRUE;
       }
       
       // ISR photon
       if(p->DistinctMother() && p->DistinctMother()->IsParton()) isWgamma=kTRUE;
    }  
    
    if( (abs(p->PdgId()) == 24) && (p->Status() == 3) ) {
      boson = p;

      // Loop through daughters and look for leptons  
      for(UInt_t j=0; j<boson->NDaughters(); j++) {
        const MCParticle *d = boson->Daughter(j);
        if(d->PdgId() == boson->PdgId()) continue;  
        
        // traverse down daughter lepton tree
        while(d->HasDaughter(d->PdgId(),kTRUE) && (d->Status()!=1))
          d = d->FindDaughter(d->PdgId(),kTRUE);          

        if(d->PdgId()== 11) id1 =  EGenType::kElectron;
        if(d->PdgId()==-11) id2 = -EGenType::kElectron;
        if(d->PdgId()== 13) id1 =  EGenType::kMuon;
        if(d->PdgId()==-13) id2 = -EGenType::kMuon;
        if(d->PdgId()== 15) id1 =  EGenType::kTau;
        if(d->PdgId()==-15) id2 = -EGenType::kTau;
        
        if(abs(d->PdgId())==15) {
          if(d->HasDaughter(11)) { 
            d = d->FindDaughter(11);
            if(d->PdgId()== 11) id1 =  EGenType::kTauElectron;
            if(d->PdgId()==-11) id2 = -EGenType::kTauElectron; 
          }
          if(d->HasDaughter(13)) { 
            d = d->FindDaughter(13);
            if(d->PdgId()== 13) id1 =  EGenType::kTauMuon;
            if(d->PdgId()==-13) id2 = -EGenType::kTauMuon; 
          }
        }
        
        if(d->PdgId()>0) dau1 = d;
        if(d->PdgId()<0) dau2 = d; 
      }               
    }                            
  }
  
  assert(boson);
  if(!dau1 || !dau2) {
    while(boson->HasDaughter(boson->PdgId(),kTRUE))
      boson = boson->FindDaughter(boson->PdgId(),kTRUE);
    
    // Loop through daughters and look for leptons  
    for(UInt_t j=0; j<boson->NDaughters(); j++) {
      const MCParticle *d = boson->Daughter(j);
      if(d->PdgId() == boson->PdgId()) continue;  
      
      // traverse down daughter lepton tree
      while(d->HasDaughter(d->PdgId(),kTRUE) && (d->Status()!=1))
        d = d->FindDaughter(d->PdgId(),kTRUE);          

      if(d->PdgId()== 11) id1 =  EGenType::kElectron;
      if(d->PdgId()==-11) id2 = -EGenType::kElectron;
      if(d->PdgId()== 13) id1 =  EGenType::kMuon;
      if(d->PdgId()==-13) id2 = -EGenType::kMuon;
      if(d->PdgId()== 15) id1 =  EGenType::kTau;
      if(d->PdgId()==-15) id2 = -EGenType::kTau;

      if(abs(d->PdgId())==15) {
        if(d->HasDaughter(11)) { 
          d = d->FindDaughter(11);
          if(d->PdgId()== 11) id1 =  EGenType::kTauElectron;
          if(d->PdgId()==-11) id2 = -EGenType::kTauElectron; 
        }
        if(d->HasDaughter(13)) { 
          d = d->FindDaughter(13);
          if(d->PdgId()== 13) id1 =  EGenType::kTauMuon;
          if(d->PdgId()==-13) id2 = -EGenType::kTauMuon; 
        }
      }

      if(d->PdgId()>0) dau1 = d;
      if(d->PdgId()<0) dau2 = d; 
    }
  }
  assert(dau1);
  assert(dau2);
  
  // Get the status=3 daughters. Not applicable for HORACE or PHOTOS modes.
  const MCParticle *vdau1=dau1;
  const MCParticle *vdau2=dau2;
  if(fFSRMode==0) {
    while(vdau1->Status()!=3)
      vdau1 = vdau1->FindMother(vdau1->PdgId(),kTRUE);
  
    while(vdau2->Status()!=3)
      vdau2 = vdau2->FindMother(vdau2->PdgId(),kTRUE);
  }  
  
  FourVectorM vDilep = dau1->Mom() + dau2->Mom();
  
  fGenInfo.pid_1    = fMCEvtInfo->Id1();
  fGenInfo.pid_2    = fMCEvtInfo->Id2();
  fGenInfo.x_1      = fMCEvtInfo->X1();
  fGenInfo.x_2      = fMCEvtInfo->X2();
  fGenInfo.scalePDF = fMCEvtInfo->ScalePdf();
  fGenInfo.weight   = fMCEvtInfo->Weight();
  fGenInfo.vmass    = boson->Mass();
  fGenInfo.vpt      = boson->Pt();
  fGenInfo.vy       = boson->Rapidity();
  fGenInfo.vphi     = boson->Phi();
  fGenInfo.vpt_1    = vdau1->Pt(); 
  fGenInfo.veta_1   = vdau1->Eta(); 
  fGenInfo.vphi_1   = vdau1->Phi();
  fGenInfo.vpt_2    = vdau2->Pt();
  fGenInfo.veta_2   = vdau2->Eta(); 
  fGenInfo.vphi_2   = vdau2->Phi();
  fGenInfo.mass     = vDilep.M();
  fGenInfo.pt       = vDilep.Pt(); 
  fGenInfo.y        = vDilep.Rapidity(); 
  fGenInfo.phi      = vDilep.Phi();
  if(isWgamma) fGenInfo.id = (boson->PdgId()>0) ? EGenType::kWgamma : -EGenType::kWgamma;
  else         fGenInfo.id = boson->PdgId();
  fGenInfo.pt_1     = dau1->Pt(); 
  fGenInfo.eta_1    = dau1->Eta(); 
  fGenInfo.phi_1    = dau1->Phi();
  fGenInfo.id_1     = id1;
  fGenInfo.pt_2     = dau2->Pt();
  fGenInfo.eta_2    = dau2->Eta(); 
  fGenInfo.phi_2    = dau2->Phi(); 
  fGenInfo.id_2     = id2;
  fGenInfo.npho     = 0;
  fGenInfo.phopt    = 0;
  fGenInfo.phoeta   = 0;
  fGenInfo.phophi   = 0;
  fGenInfo.decx     = boson->DecayVertex().X();
  fGenInfo.decy     = boson->DecayVertex().Y(); 
  fGenInfo.decz     = boson->DecayVertex().Z();
}

//--------------------------------------------------------------------------------------------------
void NtuplerMod::FillGenZjets() 
{
  assert(fParticles);
  assert(fMCEvtInfo);
  
  const MCParticle *boson=0, *dau1=0, *dau2=0;
  
  Int_t id1=0, id2=0;
  
  Bool_t isZgamma=kFALSE;
  
  for(UInt_t i=0; i<fParticles->GetEntries(); ++i) {
    const MCParticle *p = fParticles->At(i);

    if( (p->PdgId() == 22) && (p->Status() == 1) ) {
       
       // Radiated photon
       if(p->DistinctMother() && (p->DistinctMother()->Status() == 3)
          && (p->DistinctMother()->Is(MCParticle::kEl)  || p->DistinctMother()->Is(MCParticle::kMu) ||
              p->DistinctMother()->Is(MCParticle::kTau)) ) {
         isZgamma=kTRUE;
       }
       
       // ISR photon
       if(p->DistinctMother() && p->DistinctMother()->IsParton()) isZgamma=kTRUE;
    }
          
    if( (p->PdgId() == 23) && (p->Status() == 3) ) {
      boson = p;

      // Loop through daughters and look for leptons  
      for(UInt_t j=0; j<boson->NDaughters(); j++) {
        const MCParticle *d = boson->Daughter(j);
        if((abs(d->PdgId())==11) || (abs(d->PdgId())==13) || (abs(d->PdgId())==15)) {
          
          // traverse down daughter muon tree
          while(d->HasDaughter(d->PdgId(),kTRUE) && (d->Status()!=1))
            d = d->FindDaughter(d->PdgId(),kTRUE);	  

          if(d->PdgId()== 11) id1 =  EGenType::kElectron;
	  if(d->PdgId()==-11) id2 = -EGenType::kElectron;
	  if(d->PdgId()== 13) id1 =  EGenType::kMuon;
	  if(d->PdgId()==-13) id2 = -EGenType::kMuon;
	  if(d->PdgId()== 15) id1 =  EGenType::kTau;
	  if(d->PdgId()==-15) id2 = -EGenType::kTau;
	  
	  if(abs(d->PdgId())==15) {
	    if(d->HasDaughter(11)) { 
	      d = d->FindDaughter(11); 
	      if(d->PdgId()== 11) id1 =  EGenType::kTauElectron;
	      if(d->PdgId()==-11) id2 = -EGenType::kTauElectron;
	    }
	    if(d->HasDaughter(13)) { 
	      d = d->FindDaughter(13); 
	      if(d->PdgId()== 13) id1 =  EGenType::kTauMuon;
	      if(d->PdgId()==-13) id2 = -EGenType::kTauMuon;
	    }
	  }
	  
	  if(d->PdgId()>0) dau1 = d;
          if(d->PdgId()<0) dau2 = d;
        }  
      }		      
    }                            
  }
  
  assert(boson);
  assert(dau1);
  assert(dau2);
  
  // Get the status=3 daughters. Not applicable for HORACE or PHOTOS modes.
  const MCParticle *vdau1=dau1;
  const MCParticle *vdau2=dau2;
  if(fFSRMode==0) {
    while(vdau1->Status()!=3)
      vdau1 = vdau1->FindMother(vdau1->PdgId(),kTRUE);
  
    while(vdau2->Status()!=3)
      vdau2 = vdau2->FindMother(vdau2->PdgId(),kTRUE);
  }
    
  FourVectorM vDilep = dau1->Mom() + dau2->Mom();
  
  fGenInfo.pid_1    = fMCEvtInfo->Id1();
  fGenInfo.pid_2    = fMCEvtInfo->Id2();
  fGenInfo.procId   = fMCEvtInfo->ProcessId();
  fGenInfo.x_1      = fMCEvtInfo->X1();
  fGenInfo.x_2      = fMCEvtInfo->X2();
  fGenInfo.scalePDF = fMCEvtInfo->ScalePdf();
  fGenInfo.weight   = fMCEvtInfo->Weight();
  fGenInfo.vmass    = boson->Mass();
  fGenInfo.vpt      = boson->Pt();
  fGenInfo.vy       = boson->Rapidity();
  fGenInfo.vphi     = boson->Phi();
  fGenInfo.vpt_1    = vdau1->Pt(); 
  fGenInfo.veta_1   = vdau1->Eta(); 
  fGenInfo.vphi_1   = vdau1->Phi();
  fGenInfo.vpt_2    = vdau2->Pt();
  fGenInfo.veta_2   = vdau2->Eta(); 
  fGenInfo.vphi_2   = vdau2->Phi();
  fGenInfo.mass     = vDilep.M();
  fGenInfo.pt       = vDilep.Pt(); 
  fGenInfo.y        = vDilep.Rapidity(); 
  fGenInfo.phi      = vDilep.Phi(); 
  fGenInfo.id       = isZgamma ? EGenType::kZgamma : EGenType::kZ;
  fGenInfo.pt_1     = dau1->Pt(); 
  fGenInfo.eta_1    = dau1->Eta(); 
  fGenInfo.phi_1    = dau1->Phi();
  fGenInfo.id_1     = id1;
  fGenInfo.pt_2     = dau2->Pt();
  fGenInfo.eta_2    = dau2->Eta(); 
  fGenInfo.phi_2    = dau2->Phi(); 
  fGenInfo.id_2     = id2;
  fGenInfo.npho     = 0;
  fGenInfo.phopt    = 0;
  fGenInfo.phoeta   = 0;
  fGenInfo.phophi   = 0;
  fGenInfo.decx     = boson->DecayVertex().X();
  fGenInfo.decy     = boson->DecayVertex().Y(); 
  fGenInfo.decz     = boson->DecayVertex().Z();
}
