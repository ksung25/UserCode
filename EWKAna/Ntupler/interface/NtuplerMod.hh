#ifndef EWKANA_NTUPLER_NTUPLERMOD_H
#define EWKANA_NTUPLER_NTUPLERMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h"
#include "MitAna/DataTree/interface/MCEventInfo.h"
#include "MitAna/DataTree/interface/CollectionsFwd.h"
#include "MitAna/DataTree/interface/PileupInfoFwd.h"
#include "MitAna/DataTree/interface/PileupEnergyDensityFwd.h"
#include "MitAna/DataTree/interface/BaseVertex.h"
#include "MitAna/DataTree/interface/TriggerMask.h"
#include "MitAna/DataCont/interface/RunLumiRangeMap.h"
#include "MitAna/DataCont/interface/RunLumiSet.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "MitPhysics/Utils/interface/ElectronIDMVA.h"

#include "EWKAnaDefs.hh"
#include "TEventInfo.hh"
#include "TGenInfo.hh"
#include <TClonesArray.h>
#include "TMuon.hh"
#include "TElectron.hh"
#include "TJet.hh"
#include "TPhoton.hh"
#include "TVertex.hh"

#include <vector>

class TTree;
class TFile;
class TString;

namespace mithep
{  
  class NtuplerMod : public BaseMod
  {    
    public:
      NtuplerMod(const char *name="NtuplerMod", const char *title="BAMBU to ntuple");
      ~NtuplerMod();	    
      
      void SetOutputName(const char *f)         { fOutputName = f; }
       
      void SetIsData(const Bool_t flag)         { fIsData = flag; }
      void SetUseGen(const Int_t useGen)        { fUseGen = useGen; }
      void SetSkipIfHLTFail(const Bool_t flag)  { fSkipIfHLTFail = flag; } 
          
      void SetMuPtMin(const Double_t pt)        { fMuPtMin = pt; }
      void SetMuPtMax(const Double_t pt)        { fMuPtMax = pt; }
      void SetMuEtaMin(const Double_t eta)      { fMuEtaMin = eta; }
      void SetMuEtaMax(const Double_t eta)      { fMuEtaMax = eta; }
      void SetEleEtMin(const Double_t pt)       { fEleEtMin = pt; }
      void SetEleEtMax(const Double_t pt)       { fEleEtMax = pt; }
      void SetEleEtaMin(const Double_t eta)     { fEleEtaMin = eta; }
      void SetEleEtaMax(const Double_t eta)     { fEleEtaMax = eta; }
      void SetJetPtMin(const Double_t et)       { fJetPtMin = et; }
      void SetPhotonEtMin(const Double_t et)    { fPhotonEtMin = et; }
      void SetMinNTracksFit(const UInt_t ntrks) { fMinNTracksFit = ntrks; }
      void SetMinNdof(const Double_t ndof)      { fMinNdof = ndof; }
      void SetMaxAbsZ(const Double_t z)         { fMaxAbsZ = z; }
      void SetMaxRho(const Double_t rho)        { fMaxRho = rho; }
      void SetPrintHLT(const Bool_t flag)       { fPrintTable = flag; }
      void SetFSRMode(const Int_t mode)         { fFSRMode = mode; }
      
      void AddTrigger(const char* name, const ULong64_t id,
                      const char* objName1="", const ULong64_t objId1=0, const Double_t minPt1=0,
		      const char* objName2="", const ULong64_t objId2=0, const Double_t minPt2=0) {
        fTriggerNamesv.push_back(name);
	fTriggerIdsv.push_back(id);
	
	fTriggerObjNames1v.push_back(objName1);
        fTriggerObjIds1v.push_back(objId1);
	fTriggerObjMinPt1v.push_back(minPt1);
	
	fTriggerObjNames2v.push_back(objName2);
        fTriggerObjIds2v.push_back(objId2);
	fTriggerObjMinPt2v.push_back(minPt2);
      }
      
      void AddJetCorr(const char* name) {
        fJetCorrParsv.push_back(name);
      }
      
    protected:
      void Begin();
      void BeginRun();
      void EndRun();
      void SlaveBegin();
      void SlaveTerminate();
      void Terminate();
      void Process();
      
      // Fill MC info
      void FillGenZ();
      void FillGenH();
      void FillGenW();
      void FillGenWW();
      void FillGenVZ();
      void FillGenWjets();
      void FillGenZjets();
      
      // Fill muon data object
      void FillMuon(const Muon *mu);
      void FillMuon(const Track *mu);
      
      // Fill electron data object
      void FillElectron(const Electron *ele);
      
      // Fill jet data object
      void FillJet(const PFJet *jet);
      
      // Fill photon data object
      void FillPhoton(const Photon *pho);
      
      // Fill vertex data object
      void FillPV(const Vertex *pv);
      
      // Match muon to HLT primitive
      ULong64_t MatchHLT(const Double_t eta, const Double_t phi);
      ULong64_t MatchHLT(const Double_t pt, const Double_t eta, const Double_t phi);
      
      // Check for conversion with MVF
      Bool_t IsConversion(const Electron *ele);
      
      // Compute PF isolation
      Float_t computePFIso(const Electron *ele, const Double_t ptMin,
                           const Double_t extRadius, const Double_t intRadius,
                           const Int_t isoType);
      Float_t computePFIso(const Track *track, const Double_t ptMin,
                           const Double_t extRadius, const Double_t intRadius,
                           const Int_t isoType);		   
      Float_t computePUIso(const Electron *ele, const Double_t ptMin, 
                           const Double_t extRadius, const Double_t intRadius);
      Float_t computePUIso(const Track *track, const Double_t ptMin, 
                           const Double_t extRadius, const Double_t intRadius);
      void separatePileUp(Bool_t checkClosestZVertex);
      
      TFile                        *fOutputFile;           // output file handle
      TString                       fOutputName;           // output file name
      
      TString                       fPartName;             // MC particle collection name
      TString                       fMCEvtInfoName;        // MC event info name
      TString                       fMuonName;             // muon collection name
      TString                       fElectronName;         // electron collection name
      TString                       fPrimVtxName;          // primary vertex collection name
      TString                       fBeamSpotName;         // pointer to beam spot branch
      TString                       fPFJetName;            // particle flow jet collection name
      TString                       fPhotonName;           // photon collection name
      TString                       fTrigMaskName;         // trigger mask name
      TString                       fPFMetName;            // particle flow MET collection name
      TString                       fConversionName;       // conversion collection name
      TString                       fPileupName;           // pile-up info name
      TString                       fPUEnergyDensityName;  // Fastjet correction info name
      TString                       fPFCandidateName;      // particle flow candidates collection name
      TString                       fTracksName;           // track collection name           
      
      const MCParticleCol          *fParticles;       // MC particle collection handle
      const MCEventInfo            *fMCEvtInfo;       // MC event info handle
      const MuonCol                *fMuons;           // muon collection handle
      const ElectronCol            *fElectrons;       // electron collection handle
      const VertexCol              *fPrimVerts;       // primary vertex collection handle
      const BeamSpotCol            *fBeamSpot;        // pointer to beam spot branch
      const PFJetCol               *fPFJets;          // particle flow jet collection handle
      const PhotonCol              *fPhotons;         // photon collection handle
      const TriggerMask            *fTrigMask;        // trigger mask handle
      const MetCol                 *fTCMet;           // track-corrected MET handle
      const PFMetCol               *fPFMet;           // particle flow MET handle
      const DecayParticleCol       *fConversions;     // conversion collection handle
      const PileupInfoCol          *fPileup;          // pile-up info handle
      const PileupEnergyDensityCol *fPUEnergyDensity; // Fastjet correction info handle
      const PFCandidateCol         *fPFCandidates;    // particle flow candidates collection handle
      PFCandidateOArr              *fPFPileUp;
      PFCandidateOArr              *fPFNoPileUp;
      const TrackCol               *fTracks;          // tracks collection handle
      
      Bool_t                  fIsData;          // flag to indicate if processing collision data
      Int_t                   fUseGen;          // flag whether to look at generator info
      Bool_t                  fPrintTable;      // flag whether to print out HLT table
      Bool_t                  fSkipIfHLTFail;   // flag whether to skip event processing if HLT does not accept
       
      Double_t                fMuPtMin;         // minimum reco muon pT
      Double_t                fMuPtMax;         // maximum reco muon pT
      Double_t                fMuEtaMin;        // minimum reco muon eta
      Double_t                fMuEtaMax;        // maximum reco muon eta
      Double_t                fEleEtMin;        // minimum electron supercluster ET
      Double_t                fEleEtMax;        // maximum electron supercluster ET
      Double_t                fEleEtaMin;       // minimum electron supercluster eta
      Double_t                fEleEtaMax;       // maximum electron supercluster eta
      Double_t                fJetPtMin;        // minimum jet ET
      Double_t                fPhotonEtMin;     // minimum photon supercluster ET
      UInt_t                  fMinNTracksFit;   // minimum number of tracks used for a good primary vertex
      Double_t                fMinNdof;         // minimum degrees of freedom for a good primary vertex
      Double_t                fMaxAbsZ;         // maximum z displacement for a good primary vertex
      Double_t                fMaxRho;          // maximum transverse displacement for a good primary vertex 
      
      TTree*                  fEventTree;       // event tree

      const Vertex*           fVertex;          // best primary vertex in the event
            
      TEventInfo              fEventInfo;       // general event information
      TGenInfo                fGenInfo;         // generator information
      TClonesArray           *fMuonArr;         // muon array
      TClonesArray           *fElectronArr;     // electron array
      TClonesArray           *fPFJetArr;        // particle flow jet array
      TClonesArray           *fPhotonArr;       // photon array
      TClonesArray           *fPVArr;           // valid primary vertex array
      
      vector<TString>         fTriggerNamesv;   // names of triggers we're interested in 
      vector<ULong64_t>       fTriggerIdsv;     // corresponding ETriggerBit value
      
      vector<TString>         fTriggerObjNames1v;
      vector<ULong64_t>       fTriggerObjIds1v;
      vector<Double_t>        fTriggerObjMinPt1v;

      vector<TString>         fTriggerObjNames2v;
      vector<ULong64_t>       fTriggerObjIds2v;
      vector<Double_t>        fTriggerObjMinPt2v;
      
      Int_t                   fFSRMode;         // flag to indicate how to parse GEN-level tree (depends on MC tool)
                                                // 0 -> PYTHIA: traverse down tree from status=3 boson
                                                // 1 -> HORACE: look at immediate daughters of status=3 boson
                                                // 2 -> PHOTOS: look at immediate daughters of status=2 boson 

      vector<TString>         fJetCorrParsv;
      FactorizedJetCorrector *fJetCorrector;    // CMSSW class to handle jet corrections
      
      ElectronIDMVA          *fEleMVA;      
      
      RunLumiSet                       fRunLumiSet;
      RunLumiRangeMap::RunLumiPairType fLastRunLumi;
      
    ClassDef(NtuplerMod,1)
  };
}

#endif
