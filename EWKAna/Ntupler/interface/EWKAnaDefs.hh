#ifndef EWKANA_NTUPLER_EWKANADEFS_HH 
#define EWKANA_NTUPLER_EWKANADEFS_HH

namespace EGenType {
enum {
  kMuon        = 1,
  kElectron    = 2,
  kTau         = 3,
  kTauMuon     = 4,
  kTauElectron = 5,
  kW           = 6,
  kZ           = 7,
  kWW          = 8,
  kHiggs       = 9,
  kVZ          = 10,
  kWgamma      = 11,
  kZgamma      = 12
};
}

enum EMuType 
{ 
  kGlobal     = 1, 
  kTracker    = 2, 
  kStandalone = 4
};

enum EEleType
{
  kEcalDriven    = 1,
  kTrackerDriven = 2
};

enum EQualityBit
{ 
  // descriptions from DataFormats/MuonReco/interface/MuonSelectors.h
  kAll  			    = 0x000001,  // dummy options - always true
  kAllGlobalMuons		    = 0x000002,  // checks isGlobalMuon flag
  kAllStandAloneMuons		    = 0x000004,  // checks isStandAloneMuon flag
  kAllTrackerMuons		    = 0x000008,  // checks isTrackerMuon flag
  kTrackerMuonArbitrated	    = 0x000010,  // resolve ambiguity of sharing segments
  kAllArbitrated		    = 0x000020,  // all muons with the tracker muon arbitrated
  kGlobalMuonPromptTight	    = 0x000040,  // global muons with tighter fit requirements
  kTMLastStationLoose		    = 0x000080,  // penetration depth loose selector
  kTMLastStationTight		    = 0x000100,  // penetration depth tight selector
  kTM2DCompatibilityLoose	    = 0x000200,  // likelihood based loose selector
  kTM2DCompatibilityTight	    = 0x000400,  // likelihood based tight selector
  kTMOneStationLoose		    = 0x000800,  // require one well matched segment
  kTMOneStationTight		    = 0x001000,  // require one well matched segment
  kTMLastStationOptimizedLowPtLoose = 0x002000,  // combination of TMLastStation and TMOneStation
  kTMLastStationOptimizedLowPtTight = 0x004000,  // combination of TMLastStation and TMOneStation
  kGMTkChiCompatibility 	    = 0x008000,  // require tk stub have good chi2 relative to glb track
  kGMStaChiCompatibility	    = 0x010000,  // require sta stub have good chi2 compatibility relative to glb track
  kGMTkKinkTight		    = 0x020000,  // require a small kink value in the tracker stub
  kTMLastStationAngLoose	    = 0x040000,  // TMLastStationLoose with additional angular cuts
  kTMLastStationAngTight	    = 0x080000,  // TMLastStationTight with additional angular cuts
  kTMOneStationAngLoose 	    = 0x100000,  // TMOneStationLoose with additional angular cuts
  kTMOneStationAngTight 	    = 0x200000,  // TMOneStationTight with additional angular cuts
  //The two algorithms that follow are identical to what were known as
  //TMLastStationOptimizedLowPt* (sans the Barrel) as late as revision
  //1.7 of this file. The names were changed because indeed the low pt
  //optimization applies only to the barrel region, whereas the sel-
  //ectors above are more efficient at low pt in the endcaps, which is
  //what we feel is more suggestive of the algorithm name. This will be
  //less confusing for future generations of CMS members, I hope...
  kTMLastStationOptimizedBarrelLowPtLoose = 0x400000,  // combination of TMLastStation and TMOneStation but with low pT optimization in barrel only
  kTMLastStationOptimizedBarrelLowPtTight = 0x800000   // combination of TMLastStation and TMOneStation but with low pT optimization in barrel only
}; 

enum ETriggerBit
{  
  // SingleMu
  kHLT_Mu15        = 1UL<<0,
  kHLT_Mu15_eta2p1 = 1UL<<1,
  
  // SingleElectron
  kHLT_Ele17_CaloIdL_CaloIsoVL = 1UL<<2,
  kHLT_Ele22_CaloIdL_CaloIsoVL = 1UL<<3
};

enum ETriggerObjBit
{ 
  // SingleMu
  kHLT_Mu15_MuObj        = 1UL<<0,
  kHLT_Mu15_eta2p1_MuObj = 1UL<<1,
  
  // SingleElectron
  kHLT_Ele17_CaloIdL_CaloIsoVL_EleObj = 1UL<<2,
  kHLT_Ele22_CaloIdL_CaloIsoVL_EleObj = 1UL<<3
};
#endif
