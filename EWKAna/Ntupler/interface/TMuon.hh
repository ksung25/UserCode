#ifndef EWKANA_NTUPLER_TMUON_HH
#define EWKANA_NTUPLER_TMUON_HH

#include <TObject.h>

namespace mithep 
{
  class TMuon : public TObject
  {
    public:
      TMuon(){}
      ~TMuon(){} 
  
      Float_t   pt, ptErr, eta, phi;    // kinematics
      Float_t   staPt, staEta, staPhi;  // standalone muon measurements
      Float_t   trkIso03;	        // track isolation
      Float_t   emIso03;	        // ECAL-based isolation
      Float_t   hadIso03;	        // HCAL-based isolation
      Float_t   pfChIso03, pfChIso04;   // Particle Flow charged isolation
      Float_t   pfNeuIso03, pfNeuIso04; // Particle Flow neutral hadron isolation
      Float_t   pfGamIso03, pfGamIso04; // Particle Flow gamma isolation
      Float_t   puIso03, puIso04;       // Particle Flow isolation from PU
      Float_t   pfPx, pfPy;             // Matching Particle Flow candidate (px,py)
      Float_t   d0, dz;                 // impact parameter
      Float_t   tkNchi2;	        // track chi^2/ndf 
      Float_t   muNchi2;	        // global muon chi^2/ndf
      Float_t   trkKink;                // kink of tracker track
      Float_t   gblKink;                // kink of global track
      Float_t   mva;                    // MVA muon ID (not filled yet)
      Int_t     q;		        // charge
      Int_t     nValidHits;	        // number of valid hits in muon system
      UInt_t    qualityBits;            // bits for various muon quality criteria
      UInt_t    typeBits;	        // global muon, tracker muon, or standalone muon
      UInt_t    nTkHits;	        // number of inner tracker hits
      UInt_t    nPixHits;	        // number of pixel hits
      UInt_t    nTkLayers;	        // number of inner tracker layers
      UInt_t    nPixLayers;	        // number of pixel layers
      UInt_t    nSeg;  	                // number of muon segments
      UInt_t    nMatch;                 // number of muon chambers matched to segments      
      UInt_t    trkID;                  // tracker track ID
      ULong64_t hltMatchBits;           // bits for matching with HLT primitives 

    ClassDef(TMuon,1)
  };  
}
#endif
