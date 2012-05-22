#ifndef EWKANA_NTUPLER_TGENINFO_HH
#define EWKANA_NTUPLER_TGENINFO_HH

#include <TObject.h>

namespace mithep 
{
  // Generator level info data object
  class TGenInfo : public TObject
  {
    public:
      TGenInfo(){}
      ~TGenInfo(){}
      
      UInt_t  npho;                      // number of FSR photons (not filled yet)
      Int_t   pid_1, pid_2;	         // parton ID
      Int_t   id;                        // parent boson ID
      Int_t   id_1, id_2;                // lepton ID
      Int_t   procId;                    // MC process ID
      Float_t x_1, x_2;		         // parton momentum fraction
      Float_t scalePDF;                  // Q-value used for PDF evaluation
      Float_t weight;		         // event weight
      Float_t vmass, vpt, vy, vphi;      // boson info
      Float_t vpt_1, veta_1, vphi_1;     // lepton info (pre FSR)
      Float_t vpt_2, veta_2, vphi_2;
      Float_t mass, pt, y, phi;          // dilepton info
      Float_t pt_1, eta_1, phi_1;        // lepton info
      Float_t pt_2, eta_2, phi_2;  
      Float_t phopt, phoeta, phophi;     // leading photon kinematics (not filled yet)
      Float_t decx, decy, decz;	         // boson decay vertex
      	  
    ClassDef(TGenInfo,1)
  };
}
#endif
