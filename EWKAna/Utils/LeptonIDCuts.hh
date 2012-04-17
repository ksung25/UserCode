#ifndef EWKANA_UTILS_LEPTONIDCUTS_HH
#define EWKANA_UTILS_LEPTONIDCUTS_HH

#include "EWKAna/Ntupler/interface/TElectron.hh"
#include "EWKAna/Ntupler/interface/TMuon.hh"
#include <TMath.h>
#include <cassert>
		   
Bool_t passMuonID(const mithep::TMuon *muon, const Double_t rho=0);
Bool_t passMuonLooseID(const mithep::TMuon *muon, const Double_t rho=0);
Bool_t passEleID(const mithep::TElectron *electron, const Double_t rho=0);
Bool_t passEleLooseID(const mithep::TElectron *electron, const Double_t rho=0);


//=== FUNCTION DEFINITIONS ======================================================================================

//--------------------------------------------------------------------------------------------------
Bool_t passMuonID(const mithep::TMuon *muon, const Double_t rho)
{
  if(muon->nTkHits    < 11)        return kFALSE;
  if(muon->nPixHits   < 1)         return kFALSE;
  if(fabs(muon->d0)   > 0.1)       return kFALSE;
  if(fabs(muon->dz)   > 0.1)       return kFALSE;
  if(muon->muNchi2    > 10)        return kFALSE;
  if(muon->nMatch     < 2)         return kFALSE;
  if(muon->nValidHits < 1)         return kFALSE;
  if(!(muon->typeBits & kGlobal))  return kFALSE;
  if(!(muon->typeBits & kTracker)) return kFALSE;
  
  Double_t iso = muon->trkIso03 + muon->emIso03 + muon->hadIso03 - rho*0.3*0.3*TMath::Pi();
  if(iso > 0.15*(muon->pt)) return kFALSE;

  return kTRUE;
}

//--------------------------------------------------------------------------------------------------
Bool_t passMuonLooseID(const mithep::TMuon *muon, const Double_t rho)
{
  if(muon->nPixHits   < 1)         return kFALSE;
  if(muon->muNchi2    > 10)        return kFALSE;
  if(muon->nValidHits < 1)         return kFALSE;
  if(!(muon->typeBits & kGlobal))  return kFALSE;
  if(!(muon->typeBits & kTracker)) return kFALSE;

  return kTRUE;
}

//--------------------------------------------------------------------------------------------------
Bool_t passEleID(const mithep::TElectron *electron, const Double_t rho)
{
  const Double_t ECAL_GAP_LOW  = 1.4442;
  const Double_t ECAL_GAP_HIGH = 1.566;
  
  if((fabs(electron->scEta)>ECAL_GAP_LOW) && (fabs(electron->scEta)<ECAL_GAP_HIGH)) return kFALSE;
  
  if(!(electron->typeBits & kEcalDriven)) return kFALSE;
  
  if(fabs(electron->d0) > 0.1) return kFALSE;
  if(fabs(electron->dz) > 0.1) return kFALSE;
  
  // conversion rejection
  if(electron->nExpHitsInner > 0) return kFALSE;
  if(electron->isConv)            return kFALSE;
     
  // barrel/endcap dependent requirments      
  if(fabs(electron->scEta)<=ECAL_GAP_LOW) {
    // barrel
    Double_t iso = electron->trkIso03 + TMath::Max(electron->emIso03-1,Float_t(0)) + electron->hadIso03 - rho*0.3*0.3*TMath::Pi();
    if(iso > 0.07*(electron->pt)) return kFALSE;
     
    if(electron->sigiEtaiEta	  > 0.01)  return kFALSE;
    if(fabs(electron->deltaPhiIn) > 0.06)  return kFALSE;
    if(fabs(electron->deltaEtaIn) > 0.004) return kFALSE;
    if(electron->HoverE 	  > 0.04)  return kFALSE;
  
  } else {
    // endcap
    Double_t iso = electron->trkIso03 + electron->emIso03 + electron->hadIso03 - rho*0.3*0.3*TMath::Pi();
    if(iso > 0.06*(electron->pt)) return kFALSE;
     
    if(electron->sigiEtaiEta	  > 0.03)  return kFALSE;
    if(fabs(electron->deltaPhiIn) > 0.03)  return kFALSE;
    if(fabs(electron->deltaEtaIn) > 0.007) return kFALSE;
    if(electron->HoverE 	  > 0.10)  return kFALSE;
  }

  return kTRUE;
}

//--------------------------------------------------------------------------------------------------
Bool_t passEleLooseID(const mithep::TElectron *electron, const Double_t rho)
{
  const Double_t ECAL_GAP_LOW  = 1.4442;
  const Double_t ECAL_GAP_HIGH = 1.566;
  
  if((fabs(electron->scEta)>ECAL_GAP_LOW) && (fabs(electron->scEta)<ECAL_GAP_HIGH)) return kFALSE;
  
  if(!(electron->typeBits & kEcalDriven)) return kFALSE;
  
  // conversion rejection
  if(electron->nExpHitsInner > 1) return kFALSE;
  if(electron->isConv)            return kFALSE;
     
  // barrel/endcap dependent requirments      
  if(fabs(electron->scEta)<=ECAL_GAP_LOW) {
    // barrel
    Double_t iso = electron->trkIso03 + TMath::Max(electron->emIso03-1,Float_t(0)) + electron->hadIso03 - rho*0.3*0.3*TMath::Pi();
    if(iso > 0.15*(electron->pt)) return kFALSE;
     
    if(electron->sigiEtaiEta	  > 0.01)  return kFALSE;
    if(fabs(electron->deltaPhiIn) > 0.8)   return kFALSE;
    if(fabs(electron->deltaEtaIn) > 0.007) return kFALSE;
    if(electron->HoverE 	  > 0.15)  return kFALSE;
  
  } else {
    // endcap
    Double_t iso = electron->trkIso03 + electron->emIso03 + electron->hadIso03 - rho*0.3*0.3*TMath::Pi();
    if(iso > 0.10*(electron->pt)) return kFALSE;
     
    if(electron->sigiEtaiEta	  > 0.03) return kFALSE;
    if(fabs(electron->deltaPhiIn) > 0.7)  return kFALSE;
    if(fabs(electron->deltaEtaIn) > 0.01) return kFALSE;
    if(electron->HoverE 	  > 0.10) return kFALSE;
  }

  return kTRUE;
}
#endif

