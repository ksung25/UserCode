
//
// Mode
//

enum Mode {
    MuonTagAndProbe         =1,
    MuonFakeRate            =2,
    ElectronTagAndProbe     =3,
    ElectronFakeRate        =4,
    MuonTagAndProbeMC       =5,
    ElectronTagAndProbeMC   =6
};

//
// tag and probe / FR options
//

enum Option {

    //
    // Muon
    //

    MuonIsoDenominator2012=10,
    MuonIDDenominator2012=11,
    MuonFO=12,
    MuonIDIsoDenominator=13,
    MuonIDIsoRndTagDenominator=14,

    MuonIsoNumerator2011=200,
    MuonIDNumerator2011=201,
    MuonIDIsoNumerator2011=202,
    MuonIsoNumerator2012=20,
    MuonIDNumerator2012=21,
    MuonIDIsoNumerator2012=22,
    MuonTrigSglNumerator2012=23,
    MuonTrigDblLeadNumerator2012=24,
    MuonTrigDblTrailNumerator2012=25,

    //
    // Electron
    //

    ElectronIsoDenominator2012=30,
    ElectronIDDenominator2012=31,
    ElectronFO=32,
    ElectronIDIsoDenominator=33,
    ElectronIDIsoRndTagDenominator=34,

    ElectronIsoNumerator2011=400,
    ElectronIDNumerator2011=401,
    ElectronIDIsoNumerator2011=402,
    ElectronIsoNumerator2012=40,
    ElectronIDNumerator2012=41,
    ElectronIDIsoNumerator2012=42,
    ElectronTrigSglNumerator2012=43,
    ElectronTrigDblLeadNumerator2012=44,
    ElectronTrigDblTrailNumerator2012=45,

};

//
// triggers
//

struct TriggerResults {
    UInt_t HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_tag_;
    UInt_t HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_probe_;
    UInt_t HLT_Mu17_Mu8_TrailingLeg_tag_;
    UInt_t HLT_Mu17_Mu8_TrailingLeg_probe_;
    UInt_t HLT_Mu17_Mu8_LeadingLeg_tag_;
    UInt_t HLT_Mu17_Mu8_LeadingLeg_probe_;
    UInt_t HLT_Mu17_Mu8_tag_;
    UInt_t HLT_Mu17_Mu8_probe_;
    UInt_t HLT_Mu17_TkMu8_TrailingLeg_tag_;
    UInt_t HLT_Mu17_TkMu8_TrailingLeg_probe_;
    UInt_t HLT_Mu17_TkMu8_LeadingLeg_tag_;
    UInt_t HLT_Mu17_TkMu8_LeadingLeg_probe_;
    UInt_t HLT_Mu17_TkMu8_TrailingLegTrkFiltered_tag_;
    UInt_t HLT_Mu17_TkMu8_TrailingLegTrkFiltered_probe_;
    UInt_t HLT_Mu17_TkMu8_tag_;
    UInt_t HLT_Mu17_TkMu8_probe_;
    UInt_t HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_tag_;
    UInt_t HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_probe_;
    UInt_t HLT_IsoMu24_eta2p1_tag_;
    UInt_t HLT_IsoMu24_eta2p1_probe_;
    UInt_t HLT_Ele17_Ele8_L1sL1DoubleEG137_tag_;
    UInt_t HLT_Ele17_Ele8_L1sL1DoubleEG137_probe_;
    UInt_t HLT_Ele17_Ele8_LeadingLeg_tag_;
    UInt_t HLT_Ele17_Ele8_LeadingLeg_probe_;
    UInt_t HLT_Ele17_Ele8_TrailingLeg_tag_;
    UInt_t HLT_Ele17_Ele8_TrailingLeg_probe_;
    UInt_t HLT_Ele17_Ele8_tag_;
    UInt_t HLT_Ele17_Ele8_probe_;
    UInt_t HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_tag_;
    UInt_t HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_probe_;
    UInt_t HLT_Ele27_WP80_tag_;
    UInt_t HLT_Ele27_WP80_probe_;
    UInt_t HLT_Ele17_Ele8_Mass50_LeadingLeg_tag_;
    UInt_t HLT_Ele17_Ele8_Mass50_LeadingLeg_probe_;
    UInt_t HLT_Ele17_Ele8_Mass50_TrailingLeg_tag_;
    UInt_t HLT_Ele17_Ele8_Mass50_TrailingLeg_probe_;
    UInt_t HLT_Ele20_SC4_Mass50_LeadingLeg_tag_;
    UInt_t HLT_Ele20_SC4_Mass50_LeadingLeg_probe_;
    UInt_t HLT_Ele20_SC4_Mass50_TrailingLeg_tag_;
    UInt_t HLT_Ele20_SC4_Mass50_TrailingLeg_probe_;
    UInt_t HLT_Ele32_SC17_Mass50_LeadingLeg_tag_;
    UInt_t HLT_Ele32_SC17_Mass50_LeadingLeg_probe_;
    UInt_t HLT_Ele32_SC17_Mass50_TrailingLeg_tag_;
    UInt_t HLT_Ele32_SC17_Mass50_TrailingLeg_probe_;
    UInt_t HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_probe_;
    UInt_t HLT_Mu8_probe_;
};

//
// map branches
//

void mapExtraBranches(LeptonTree *leptonTree, TriggerResults &triggerResults) 
{
    leptonTree->tree_->SetBranchAddress("HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_tag", &triggerResults.HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_tag_);
    leptonTree->tree_->SetBranchAddress("HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_probe", &triggerResults.HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_probe_);
    leptonTree->tree_->SetBranchAddress("HLT_Mu17_Mu8_TrailingLeg_tag", &triggerResults.HLT_Mu17_Mu8_TrailingLeg_tag_);
    leptonTree->tree_->SetBranchAddress("HLT_Mu17_Mu8_TrailingLeg_probe", &triggerResults.HLT_Mu17_Mu8_TrailingLeg_probe_);
    leptonTree->tree_->SetBranchAddress("HLT_Mu17_Mu8_LeadingLeg_tag", &triggerResults.HLT_Mu17_Mu8_LeadingLeg_tag_);
    leptonTree->tree_->SetBranchAddress("HLT_Mu17_Mu8_LeadingLeg_probe", &triggerResults.HLT_Mu17_Mu8_LeadingLeg_probe_);
    leptonTree->tree_->SetBranchAddress("HLT_Mu17_Mu8_tag", &triggerResults.HLT_Mu17_Mu8_tag_);
    leptonTree->tree_->SetBranchAddress("HLT_Mu17_Mu8_probe", &triggerResults.HLT_Mu17_Mu8_probe_);
    leptonTree->tree_->SetBranchAddress("HLT_Mu17_TkMu8_TrailingLeg_tag", &triggerResults.HLT_Mu17_TkMu8_TrailingLeg_tag_);
    leptonTree->tree_->SetBranchAddress("HLT_Mu17_TkMu8_TrailingLeg_probe", &triggerResults.HLT_Mu17_TkMu8_TrailingLeg_probe_);
    leptonTree->tree_->SetBranchAddress("HLT_Mu17_TkMu8_LeadingLeg_tag", &triggerResults.HLT_Mu17_TkMu8_LeadingLeg_tag_);
    leptonTree->tree_->SetBranchAddress("HLT_Mu17_TkMu8_LeadingLeg_probe", &triggerResults.HLT_Mu17_TkMu8_LeadingLeg_probe_);
    leptonTree->tree_->SetBranchAddress("HLT_Mu17_TkMu8_TrailingLegTrkFiltered_tag", &triggerResults.HLT_Mu17_TkMu8_TrailingLegTrkFiltered_tag_);
    leptonTree->tree_->SetBranchAddress("HLT_Mu17_TkMu8_TrailingLegTrkFiltered_probe", &triggerResults.HLT_Mu17_TkMu8_TrailingLegTrkFiltered_probe_);
    leptonTree->tree_->SetBranchAddress("HLT_Mu17_TkMu8_tag", &triggerResults.HLT_Mu17_TkMu8_tag_);
    leptonTree->tree_->SetBranchAddress("HLT_Mu17_TkMu8_probe", &triggerResults.HLT_Mu17_TkMu8_probe_);
    leptonTree->tree_->SetBranchAddress("HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_tag", &triggerResults.HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_tag_);
    leptonTree->tree_->SetBranchAddress("HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_probe", &triggerResults.HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_probe_);
    leptonTree->tree_->SetBranchAddress("HLT_IsoMu24_eta2p1_tag", &triggerResults.HLT_IsoMu24_eta2p1_tag_);
    leptonTree->tree_->SetBranchAddress("HLT_IsoMu24_eta2p1_probe", &triggerResults.HLT_IsoMu24_eta2p1_probe_);
    leptonTree->tree_->SetBranchAddress("HLT_Ele17_Ele8_L1sL1DoubleEG137_tag", &triggerResults.HLT_Ele17_Ele8_L1sL1DoubleEG137_tag_);
    leptonTree->tree_->SetBranchAddress("HLT_Ele17_Ele8_L1sL1DoubleEG137_probe", &triggerResults.HLT_Ele17_Ele8_L1sL1DoubleEG137_probe_);
    leptonTree->tree_->SetBranchAddress("HLT_Ele17_Ele8_LeadingLeg_tag", &triggerResults.HLT_Ele17_Ele8_LeadingLeg_tag_);
    leptonTree->tree_->SetBranchAddress("HLT_Ele17_Ele8_LeadingLeg_probe", &triggerResults.HLT_Ele17_Ele8_LeadingLeg_probe_);
    leptonTree->tree_->SetBranchAddress("HLT_Ele17_Ele8_TrailingLeg_tag", &triggerResults.HLT_Ele17_Ele8_TrailingLeg_tag_);
    leptonTree->tree_->SetBranchAddress("HLT_Ele17_Ele8_TrailingLeg_probe", &triggerResults.HLT_Ele17_Ele8_TrailingLeg_probe_);
    leptonTree->tree_->SetBranchAddress("HLT_Ele17_Ele8_tag", &triggerResults.HLT_Ele17_Ele8_tag_);
    leptonTree->tree_->SetBranchAddress("HLT_Ele17_Ele8_probe", &triggerResults.HLT_Ele17_Ele8_probe_);
    leptonTree->tree_->SetBranchAddress("HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_tag", &triggerResults.HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_tag_);
    leptonTree->tree_->SetBranchAddress("HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_probe", &triggerResults.HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_probe_);
    leptonTree->tree_->SetBranchAddress("HLT_Ele27_WP80_tag", &triggerResults.HLT_Ele27_WP80_tag_);
    leptonTree->tree_->SetBranchAddress("HLT_Ele27_WP80_probe", &triggerResults.HLT_Ele27_WP80_probe_);
    leptonTree->tree_->SetBranchAddress("HLT_Ele17_Ele8_Mass50_LeadingLeg_tag", &triggerResults.HLT_Ele17_Ele8_Mass50_LeadingLeg_tag_);
    leptonTree->tree_->SetBranchAddress("HLT_Ele17_Ele8_Mass50_LeadingLeg_probe", &triggerResults.HLT_Ele17_Ele8_Mass50_LeadingLeg_probe_);
    leptonTree->tree_->SetBranchAddress("HLT_Ele17_Ele8_Mass50_TrailingLeg_tag", &triggerResults.HLT_Ele17_Ele8_Mass50_TrailingLeg_tag_);
    leptonTree->tree_->SetBranchAddress("HLT_Ele17_Ele8_Mass50_TrailingLeg_probe", &triggerResults.HLT_Ele17_Ele8_Mass50_TrailingLeg_probe_);
    leptonTree->tree_->SetBranchAddress("HLT_Ele20_SC4_Mass50_LeadingLeg_tag", &triggerResults.HLT_Ele20_SC4_Mass50_LeadingLeg_tag_);
    leptonTree->tree_->SetBranchAddress("HLT_Ele20_SC4_Mass50_LeadingLeg_probe", &triggerResults.HLT_Ele20_SC4_Mass50_LeadingLeg_probe_);
    leptonTree->tree_->SetBranchAddress("HLT_Ele20_SC4_Mass50_TrailingLeg_tag", &triggerResults.HLT_Ele20_SC4_Mass50_TrailingLeg_tag_);
    leptonTree->tree_->SetBranchAddress("HLT_Ele20_SC4_Mass50_TrailingLeg_probe", &triggerResults.HLT_Ele20_SC4_Mass50_TrailingLeg_probe_);
    leptonTree->tree_->SetBranchAddress("HLT_Ele32_SC17_Mass50_LeadingLeg_tag", &triggerResults.HLT_Ele32_SC17_Mass50_LeadingLeg_tag_);
    leptonTree->tree_->SetBranchAddress("HLT_Ele32_SC17_Mass50_LeadingLeg_probe", &triggerResults.HLT_Ele32_SC17_Mass50_LeadingLeg_probe_);
    leptonTree->tree_->SetBranchAddress("HLT_Ele32_SC17_Mass50_TrailingLeg_tag", &triggerResults.HLT_Ele32_SC17_Mass50_TrailingLeg_tag_);
    leptonTree->tree_->SetBranchAddress("HLT_Ele32_SC17_Mass50_TrailingLeg_probe", &triggerResults.HLT_Ele32_SC17_Mass50_TrailingLeg_probe_);
    leptonTree->tree_->SetBranchAddress("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_probe", &triggerResults.HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_probe_);
    leptonTree->tree_->SetBranchAddress("HLT_Mu8_probe", &triggerResults.HLT_Mu8_probe_);

}

//
// probes
//

bool isProbe(unsigned int option, const LeptonTree &leptonTree, const TriggerResults &triggerResults, const Mode mode)
{

    bool isMC = false;
    if (mode == ElectronTagAndProbeMC || mode == MuonTagAndProbeMC) isMC = true;

    //
    // electrons
    //

    if (option == ElectronFO) {
        if ((leptonTree.leptonSelection_ & LeptonTree::PassEleFOICHEP2012) != LeptonTree::PassEleFOICHEP2012)       return false;
        if (triggerResults.HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_probe_ == 0)                                return false;
    }

    if (option == ElectronIsoDenominator2012) {
        if (leptonTree.tag_.Pt() < 32.0)                                            return false;
        if (triggerResults.HLT_Ele27_WP80_tag_ == 0 && !isMC)                                    return false;
        if ((leptonTree.leptonSelection_ & LeptonTree::PassEleIsoICHEP2012) 
                != LeptonTree::PassEleIsoICHEP2012)                                 return false;
    }

    if (option == ElectronIDDenominator2012) {
        if (leptonTree.tag_.Pt() < 32.0)                                            return false;
        if (triggerResults.HLT_Ele27_WP80_tag_ == 0 && !isMC)                                    return false;
        if ((leptonTree.leptonSelection_ & LeptonTree::PassEleIDICHEP2012) 
                != LeptonTree::PassEleIDICHEP2012)                                  return false;
    }

    if (option == ElectronIDIsoDenominator) {
        if (leptonTree.tag_.Pt() < 32.0)                                           return false;
        if (triggerResults.HLT_Ele27_WP80_tag_ == 0 && !isMC)                      return false;
        if ((leptonTree.leptonSelection_ & LeptonTree::PassEleIsoICHEP2012)
                != LeptonTree::PassEleIsoICHEP2012)                                return false;
        if ((leptonTree.leptonSelection_ & LeptonTree::PassEleIDICHEP2012)
                != LeptonTree::PassEleIDICHEP2012)                                 return false;
    }

    if (option == ElectronIDIsoRndTagDenominator) {
        if (leptonTree.tag_.Pt() < 32.0)                                           return false;
        if (triggerResults.HLT_Ele27_WP80_tag_ == 0 && !isMC)                      return false;
        if ((leptonTree.leptonSelection_ & LeptonTree::PassEleIsoICHEP2012)
                != LeptonTree::PassEleIsoICHEP2012)                                return false;
        if ((leptonTree.leptonSelection_ & LeptonTree::PassEleIDICHEP2012)
                != LeptonTree::PassEleIDICHEP2012)                                 return false;
        if ((leptonTree.qTag_ < 0 && leptonTree.rnd_ < 0.5) ||
            (leptonTree.qTag_ > 0 && leptonTree.rnd_ >= 0.5))                      return false;
    }

    //
    // muons
    //

    if (option == MuonFO) {
        if ((leptonTree.leptonSelection_ & LeptonTree::PassMuFOICHEP2012) != LeptonTree::PassMuFOICHEP2012)       return false;
        if (triggerResults.HLT_Mu8_probe_ == 0) return false;
    }

    if (option == MuonIsoDenominator2012) {
        if (leptonTree.tag_.Pt() < 30.0)                                            return false;
        if (triggerResults.HLT_IsoMu24_eta2p1_tag_ == 0 && !isMC)                                    return false;
        if ((leptonTree.leptonSelection_ & LeptonTree::PassMuIsoICHEP2012)
                != LeptonTree::PassMuIsoICHEP2012)                                 return false;
    }

    if (option == MuonIDDenominator2012) {
        if (leptonTree.tag_.Pt() < 30.0)                                            return false;
        if (triggerResults.HLT_IsoMu24_eta2p1_tag_ == 0 && !isMC)                                    return false;
        if ((leptonTree.leptonSelection_ & LeptonTree::PassMuIDICHEP2012)
                != LeptonTree::PassMuIDICHEP2012)                                  return false;
    }

    if (option == MuonIDIsoDenominator) {
        if (leptonTree.tag_.Pt() < 30.0)                                           return false;
        if (triggerResults.HLT_IsoMu24_eta2p1_tag_ == 0 && !isMC)                  return false;
        if ((leptonTree.leptonSelection_ & LeptonTree::PassMuIsoICHEP2012)
                != LeptonTree::PassMuIsoICHEP2012)                                 return false;
        if ((leptonTree.leptonSelection_ & LeptonTree::PassMuIDICHEP2012)
                != LeptonTree::PassMuIDICHEP2012)                                  return false;
    }

    if (option == MuonIDIsoRndTagDenominator) {
        if (leptonTree.tag_.Pt() < 30.0)                                           return false;
        if (triggerResults.HLT_IsoMu24_eta2p1_tag_ == 0 && !isMC)                  return false;
        if ((leptonTree.leptonSelection_ & LeptonTree::PassMuIsoICHEP2012)
                != LeptonTree::PassMuIsoICHEP2012)                                 return false;
        if ((leptonTree.leptonSelection_ & LeptonTree::PassMuIDICHEP2012)
                != LeptonTree::PassMuIDICHEP2012)                                  return false;
        if ((leptonTree.qTag_ < 0 && leptonTree.rnd_ < 0.5) ||
            (leptonTree.qTag_ > 0 && leptonTree.rnd_ >= 0.5))                      return false;
    }

    return true;
}

bool isPassProbe(unsigned int option, const LeptonTree &leptonTree, const TriggerResults &triggerResults, const Mode mode)
{

    bool isMC = false;
    if (mode == ElectronTagAndProbeMC || mode == MuonTagAndProbeMC) isMC = true;

    //
    // Electrons
    //

    if (option == ElectronIsoNumerator2012) {
        if ((leptonTree.leptonSelection_ & LeptonTree::PassEleIsoICHEP2012) != LeptonTree::PassEleIsoICHEP2012)       return false;
    }

    if (option == ElectronIDNumerator2012) {
        if ((leptonTree.leptonSelection_ & LeptonTree::PassEleIDICHEP2012) != LeptonTree::PassEleIDICHEP2012)       return false;
    }

    if (option == ElectronIDIsoNumerator2012) {
        unsigned int IDIso = LeptonTree::PassEleIsoICHEP2012 | LeptonTree::PassEleIDICHEP2012;
        if ((leptonTree.leptonSelection_ & IDIso) != IDIso)       return false;
    }

    if (option == ElectronIsoNumerator2011) {
        if ((leptonTree.leptonSelection_ & LeptonTree::PassEleIso) != LeptonTree::PassEleIso)       return false;
    }

    if (option == ElectronIDNumerator2011) {
        if ((leptonTree.leptonSelection_ & LeptonTree::PassEleID) != LeptonTree::PassEleID)       return false;
    }

    if (option == ElectronIDIsoNumerator2011) {
        unsigned int IDIso = LeptonTree::PassEleIso | LeptonTree::PassEleID;
        if ((leptonTree.leptonSelection_ & IDIso) != IDIso)       return false;
    }

    if (option == ElectronTrigSglNumerator2012) {
        if (triggerResults.HLT_Ele27_WP80_probe_ == 0 && !isMC)   return false;
    }   
        
    if (option == ElectronTrigDblLeadNumerator2012) {
        if (triggerResults.HLT_Ele17_Ele8_LeadingLeg_probe_ == 0 && !isMC) return false;
    }   
        
    if (option == ElectronTrigDblTrailNumerator2012) {
        if (triggerResults.HLT_Ele17_Ele8_TrailingLeg_probe_ == 0 && !isMC) return false;
    }

    //
    // Muons
    //

    if (option == MuonIsoNumerator2012) {
        if ((leptonTree.leptonSelection_ & LeptonTree::PassMuIsoICHEP2012) != LeptonTree::PassMuIsoICHEP2012)       return false;
    }

    if (option == MuonIDNumerator2012) {
        if ((leptonTree.leptonSelection_ & LeptonTree::PassMuIDICHEP2012) != LeptonTree::PassMuIDICHEP2012)       return false;
    }

    if (option == MuonIDIsoNumerator2012) {
        unsigned int IDIso = LeptonTree::PassMuIsoICHEP2012 | LeptonTree::PassMuIDICHEP2012;
        if ((leptonTree.leptonSelection_ & IDIso) != IDIso)       return false;
    }

    if (option == MuonIsoNumerator2011) {
        if ((leptonTree.leptonSelection_ & LeptonTree::PassMuIso) != LeptonTree::PassMuIso)       return false;
    }

    if (option == MuonIDNumerator2011) {
        if ((leptonTree.leptonSelection_ & LeptonTree::PassMuID) != LeptonTree::PassMuID)       return false;
    }

    if (option == MuonIDIsoNumerator2011) {
        unsigned int IDIso = LeptonTree::PassMuIso | LeptonTree::PassMuID;
        if ((leptonTree.leptonSelection_ & IDIso) != IDIso)       return false;
    }

    if (option == MuonTrigSglNumerator2012) {
        if (triggerResults.HLT_IsoMu24_eta2p1_probe_ == 0 && !isMC)   return false;
    }

    if (option == MuonTrigDblLeadNumerator2012) {
        if ((triggerResults.HLT_Mu17_TkMu8_LeadingLeg_probe_ == 0
                && triggerResults.HLT_Mu17_Mu8_LeadingLeg_probe_ == 0) && !isMC) return false;
    }

    if (option == MuonTrigDblTrailNumerator2012) {
        if ((triggerResults.HLT_Mu17_TkMu8_TrailingLeg_probe_ == 0
                && triggerResults.HLT_Mu17_Mu8_TrailingLeg_probe_ == 0) && !isMC) return false;
    }

    return true;
}

