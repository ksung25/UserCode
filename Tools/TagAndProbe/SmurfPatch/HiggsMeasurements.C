
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
// Egamma cut based
//

    enum CutType {
        DETAIN          = (1<<0),
        DPHIIN          = (1<<1),
        SIGMAIETAIETA   = (1<<2),
        HOE             = (1<<3),
        OOEMOOP         = (1<<4),
        D0VTX           = (1<<5),
        DZVTX           = (1<<6),
        ISO             = (1<<7),
        VTXFIT          = (1<<8),
        MHITS           = (1<<9)
    };

// all possible cuts pass
static const unsigned int PassAll         = DETAIN | DPHIIN | SIGMAIETAIETA | HOE | OOEMOOP | D0VTX | DZVTX | ISO | VTXFIT | MHITS;

//
// tag and probe / FR options
//

enum Option {

    // single electron triggers
    DenominatorEleSingle27WP80  = 0,
    NumeratorEleSingle27WP80    = 1,

    // single muon triggers
    DenominatorMuSingleIso24    = 100,
    NumeratorMuSingleIso24      = 101,

    //
    // double electron triggers
    //

    DenominatorEle17d8Trail   = 200,
    NumeratorEle17d8Trail     = 201,
    DenominatorEle17d8Lead   = 202,
    NumeratorEle17d8Lead     = 203,
    DenominatorEle17d8Dz      = 204,
    NumeratorEle17d8Dz        = 205,

    //
    // double muon triggers
    //

    DenominatorMu17dTk8Trail  =300,
    NumeratorMu17dTk8Trail    =301,
    DenominatorMu17dTk8Lead   =302,
    NumeratorMu17dTk8Lead     =303,
    DenominatorMu17dTk8Dz     =304,
    NumeratorMu17dTk8Dz       =305,

    DenominatorMu17d8Trail  =310,
    NumeratorMu17d8Trail    =311,
    DenominatorMu17d8Lead   =312,
    NumeratorMu17d8Lead     =313,
    DenominatorMu17d8Dz     =314,
    NumeratorMu17d8Dz       =315,

    // offline denominators
    // and numerators
    //

    DenominatorMu                   = 5000,
    DenominatorMuHWW2011ID          = 5001,
    DenominatorMuHWW2011Iso         = 5002,

    NumeratorMuHWW2011ID            = 5003,
    NumeratorMuHWW2011Iso           = 5004,
    NumeratorMuHWW2011IDIso         = 5005,
    NumeratorMuHWW2011IDPlusPFMu    = 5006,

    //
    // offline numerators
    //

    DenominatorEle                   = 6000,
    DenominatorEleHWW2011ID          = 6001,
    DenominatorEleHWW2011Iso         = 6002,

    NumeratorEleHWW2011ID                = 6003,
    NumeratorEleHWW2011Iso           = 6004,
    NumeratorEleHWW2011IDIso         = 6005,
    NumeratorElePOG2012Iso          = 6006,
    NumeratorEleHWW2011IDIsoMVA     = 6007,
    NumeratorElePOG2012IDIso        = 6008,

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

// Define probe and passing probe condition
bool isProbe(unsigned int option, const LeptonTree &leptonTree, const TriggerResults &triggerResults) {

    //
    // electron triggers
    //

    if (option == DenominatorEleSingle27WP80) {
        if (leptonTree.tag_.Pt() < 30.0)                                            return false;
        if (triggerResults.HLT_Ele27_WP80_tag_ == 0)                        return false;
        if ((leptonTree.leptonSelection_ & LeptonTree::PassEleID) != LeptonTree::PassEleID)     return false;
        if ((leptonTree.leptonSelection_ & LeptonTree::PassEleIso) != LeptonTree::PassEleIso)   return false;
    }

    if (option == DenominatorEle17d8Trail) {
        if ((leptonTree.qTag_ < 0 && leptonTree.rnd_ < 0.5) ||
                (leptonTree.qTag_ > 0 && leptonTree.rnd_ >= 0.5))                        return false;
        if (leptonTree.tag_.Pt() < 30.0)                                            return false;
        if (triggerResults.HLT_Ele27_WP80_tag_ == 0)                        return false;
        if (triggerResults.HLT_Ele17_Ele8_LeadingLeg_tag_ == 0)             return false;
        if ((leptonTree.leptonSelection_ & LeptonTree::PassEleID) != LeptonTree::PassEleID)     return false;
        if ((leptonTree.leptonSelection_ & LeptonTree::PassEleIso) != LeptonTree::PassEleIso)   return false;
    }

    if (option == DenominatorEle17d8Lead) {
        if ((leptonTree.qTag_ < 0 && leptonTree.rnd_ < 0.5) ||
                (leptonTree.qTag_ > 0 && leptonTree.rnd_ >= 0.5))                        return false;
        if (leptonTree.tag_.Pt() < 30.0)                                            return false;
        if (triggerResults.HLT_Ele27_WP80_tag_ == 0)                        return false;
        if (triggerResults.HLT_Ele17_Ele8_TrailingLeg_tag_ == 0)             return false;
        if ((leptonTree.leptonSelection_ & LeptonTree::PassEleID) != LeptonTree::PassEleID)     return false;
        if ((leptonTree.leptonSelection_ & LeptonTree::PassEleIso) != LeptonTree::PassEleIso)   return false;
    }

    if (option == DenominatorEle17d8Dz) {
        if ((leptonTree.qTag_ < 0 && leptonTree.rnd_ < 0.5) ||
                (leptonTree.qTag_ > 0 && leptonTree.rnd_ >= 0.5))                        return false;
        if (leptonTree.tag_.Pt() < 30.0)                                            return false;
        if (triggerResults.HLT_Ele27_WP80_tag_ == 0)                        return false;
        if (triggerResults.HLT_Ele17_Ele8_LeadingLeg_tag_ == 0)             return false;
        if ((leptonTree.leptonSelection_ & LeptonTree::PassEleID) != LeptonTree::PassEleID)     return false;
        if ((leptonTree.leptonSelection_ & LeptonTree::PassEleIso) != LeptonTree::PassEleIso)   return false;
        if (triggerResults.HLT_Ele17_Ele8_TrailingLeg_probe_ == 0)          return false;
    }

    //
    // muon triggers
    //

    if (option == DenominatorMuSingleIso24) {
        if (leptonTree.tag_.Pt() < 30.0)                                            return false;
        if (triggerResults.HLT_IsoMu24_eta2p1_tag_ == 0)                  return false;
        if ((leptonTree.leptonSelection_ & LeptonTree::PassMuID) != LeptonTree::PassMuID)     return false;
        if ((leptonTree.leptonSelection_ & LeptonTree::PassMuIso) != LeptonTree::PassMuIso)   return false;
    }

    if (option == DenominatorMu17dTk8Trail) {
        if ((leptonTree.qTag_ < 0 && leptonTree.rnd_ < 0.5) ||
                (leptonTree.qTag_ > 0 && leptonTree.rnd_ >= 0.5))                        return false;
        if (leptonTree.tag_.Pt() < 30.0)                                            return false;
        if (triggerResults.HLT_IsoMu24_eta2p1_tag_ == 0)                  return false;
        if (triggerResults.HLT_Mu17_TkMu8_LeadingLeg_tag_ == 0)           return false;
        if ((leptonTree.leptonSelection_ & LeptonTree::PassMuID) != LeptonTree::PassMuID)     return false;
        if ((leptonTree.leptonSelection_ & LeptonTree::PassMuIso) != LeptonTree::PassMuIso)   return false;
    }

    if (option == DenominatorMu17dTk8Lead) {
        if ((leptonTree.qTag_ < 0 && leptonTree.rnd_ < 0.5) ||
                (leptonTree.qTag_ > 0 && leptonTree.rnd_ >= 0.5))                        return false;
        if (leptonTree.tag_.Pt() < 30.0)                                            return false;
        if (triggerResults.HLT_IsoMu24_eta2p1_tag_ == 0)                    return false;
        if (triggerResults.HLT_Mu17_TkMu8_TrailingLegTrkFiltered_tag_ == 0) return false;
        if ((leptonTree.leptonSelection_ & LeptonTree::PassMuID) != LeptonTree::PassMuID)       return false;
        if ((leptonTree.leptonSelection_ & LeptonTree::PassMuIso) != LeptonTree::PassMuIso)     return false;
    }

    if (option == DenominatorMu17dTk8Dz) {
        if ((leptonTree.qTag_ < 0 && leptonTree.rnd_ < 0.5) ||
                (leptonTree.qTag_ > 0 && leptonTree.rnd_ >= 0.5))                        return false;
        if (leptonTree.tag_.Pt() < 30.0)                                            return false;
        if (triggerResults.HLT_IsoMu24_eta2p1_tag_ == 0)                  return false;
        if (triggerResults.HLT_Mu17_TkMu8_LeadingLeg_tag_ == 0)           return false;
        if ((leptonTree.leptonSelection_ & LeptonTree::PassMuID) != LeptonTree::PassMuID)     return false;
        if ((leptonTree.leptonSelection_ & LeptonTree::PassMuIso) != LeptonTree::PassMuIso)   return false;
        if (triggerResults.HLT_Mu17_TkMu8_TrailingLegTrkFiltered_probe_ == 0) return false;
    }

    //

    if (option == DenominatorMu17d8Trail) {
        if ((leptonTree.qTag_ < 0 && leptonTree.rnd_ < 0.5) ||
                (leptonTree.qTag_ > 0 && leptonTree.rnd_ >= 0.5))                        return false;
        if (leptonTree.tag_.Pt() < 30.0)                                            return false;
        if (triggerResults.HLT_IsoMu24_eta2p1_tag_ == 0)                  return false;
        if (triggerResults.HLT_Mu17_Mu8_LeadingLeg_tag_ == 0)           return false;
        if ((leptonTree.leptonSelection_ & LeptonTree::PassMuID) != LeptonTree::PassMuID)     return false;
        if ((leptonTree.leptonSelection_ & LeptonTree::PassMuIso) != LeptonTree::PassMuIso)   return false;
    }

    if (option == DenominatorMu17d8Lead) {
        if ((leptonTree.qTag_ < 0 && leptonTree.rnd_ < 0.5) ||
                (leptonTree.qTag_ > 0 && leptonTree.rnd_ >= 0.5))                        return false;
        if (leptonTree.tag_.Pt() < 30.0)                                            return false;
        if (triggerResults.HLT_IsoMu24_eta2p1_tag_ == 0)                    return false;
        if (triggerResults.HLT_Mu17_Mu8_TrailingLeg_tag_ == 0) return false;
        if ((leptonTree.leptonSelection_ & LeptonTree::PassMuID) != LeptonTree::PassMuID)       return false;
        if ((leptonTree.leptonSelection_ & LeptonTree::PassMuIso) != LeptonTree::PassMuIso)     return false;
    }

    if (option == DenominatorMu17d8Dz) {
        if ((leptonTree.qTag_ < 0 && leptonTree.rnd_ < 0.5) ||
                (leptonTree.qTag_ > 0 && leptonTree.rnd_ >= 0.5))                        return false;
        if (leptonTree.tag_.Pt() < 30.0)                                            return false;
        if (triggerResults.HLT_IsoMu24_eta2p1_tag_ == 0)                  return false;
        if (triggerResults.HLT_Mu17_Mu8_LeadingLeg_tag_ == 0)           return false;
        if ((leptonTree.leptonSelection_ & LeptonTree::PassMuID) != LeptonTree::PassMuID)     return false;
        if ((leptonTree.leptonSelection_ & LeptonTree::PassMuIso) != LeptonTree::PassMuIso)   return false;
        if (triggerResults.HLT_Mu17_Mu8_TrailingLeg_probe_ == 0) return false;
    }

    //
    // offline
    //

    if (option == DenominatorEle) {
        if ((leptonTree.leptonSelection_ & LeptonTree::PassEleFO) != LeptonTree::PassEleFO)     return false;
        if (abs(leptonTree.sceta_) < 1.5) {
            if (leptonTree.sieie_ > 0.011) return false;
            if (leptonTree.hoe_ > 0.10) return false;
        } else {
            if (leptonTree.sieie_ > 0.031) return false;
            if (leptonTree.hoe_ > 0.075) return false;
        }
    }


    if (option == DenominatorEleHWW2011ID) {
        if ((leptonTree.leptonSelection_ & LeptonTree::PassEleID) != LeptonTree::PassEleID)     return false;
    }

    if (option == DenominatorEleHWW2011Iso) {
        if ((leptonTree.leptonSelection_ & LeptonTree::PassEleIso) != LeptonTree::PassEleIso)   return false;
    }


    if (option == DenominatorMu) {
        // nothing
    }

    if (option == DenominatorMuHWW2011ID) {
        if ((leptonTree.leptonSelection_ & LeptonTree::PassMuID) != LeptonTree::PassMuID)     return false;
    }

    if (option == DenominatorMuHWW2011Iso) {
        if ((leptonTree.leptonSelection_ & LeptonTree::PassMuIso) != LeptonTree::PassMuIso)   return false;
    }


    return true;

}

bool isPassProbe(unsigned int option, const LeptonTree &leptonTree, const TriggerResults &triggerResults) {

    //
    // electron triggers
    //

    if (option == NumeratorEleSingle27WP80) {
        if (triggerResults.HLT_Ele27_WP80_probe_ == 0)                  return false;
    }

    if (option == NumeratorEle17d8Trail) {
        if (triggerResults.HLT_Ele17_Ele8_TrailingLeg_probe_ == 0)          return false;
    }

    if (option == NumeratorEle17d8Lead) {
        if (triggerResults.HLT_Ele17_Ele8_LeadingLeg_probe_ == 0)          return false;
    }

    if (option == NumeratorEle17d8Dz) {
        if (triggerResults.HLT_Ele17_Ele8_probe_ == 0)          return false;
    }

    //
    // muon triggers
    //

    if (option == NumeratorMuSingleIso24) {
        if (triggerResults.HLT_IsoMu24_eta2p1_probe_ == 0)              return false;
    }
    if (option == NumeratorMu17dTk8Trail) {
        if (triggerResults.HLT_Mu17_TkMu8_TrailingLegTrkFiltered_probe_ == 0) return false;
    }
    if (option == NumeratorMu17dTk8Lead) {
        if (triggerResults.HLT_Mu17_TkMu8_LeadingLeg_probe_ == 0) return false;
    }
    if (option == NumeratorMu17dTk8Dz) {
        if (triggerResults.HLT_Mu17_TkMu8_probe_ == 0) return false;
    }

    //

    if (option == NumeratorMu17d8Trail) {
        if (triggerResults.HLT_Mu17_Mu8_TrailingLeg_probe_ == 0) return false;
    }
    if (option == NumeratorMu17d8Lead) {
        if (triggerResults.HLT_Mu17_Mu8_LeadingLeg_probe_ == 0) return false;
    }
    if (option == NumeratorMu17d8Dz) {
        if (triggerResults.HLT_Mu17_Mu8_probe_ == 0) return false;
    }

    //
    // offline
    //

    if (option == NumeratorEleHWW2011ID) {
        if ((leptonTree.leptonSelection_ & LeptonTree::PassEleID) != LeptonTree::PassEleID)     return false;
    }
    if (option == NumeratorEleHWW2011Iso) {
        if ((leptonTree.leptonSelection_ & LeptonTree::PassEleIso) != LeptonTree::PassEleIso)   return false;
    }
    if (option == NumeratorEleHWW2011IDIso) {
        if ((leptonTree.leptonSelection_ & LeptonTree::PassEleID) != LeptonTree::PassEleID)   return false;
        if ((leptonTree.leptonSelection_ & LeptonTree::PassEleIso) != LeptonTree::PassEleIso)   return false;
    }

    if (option == NumeratorElePOG2012Iso) {
        float egammaPOGNeutral = TMath::Max(float(0.0), leptonTree.el_pfnhiso04_ + leptonTree.el_pfemiso04_ - leptonTree.el_ea04_ * TMath::Max(float(0.0), leptonTree.rhoIsoAll_));
        float egammaPOGIso = (leptonTree.el_pfchiso04_ + egammaPOGNeutral) / leptonTree.probe_.Pt();
        if (egammaPOGIso > 0.15) return false;
    }


    if (option == NumeratorEleHWW2011IDIsoMVA) {

        // define kinematic bin
        Int_t subdet = 0;
        if (fabs(leptonTree.sceta_) < 1.0) subdet = 0;
        else if (fabs(leptonTree.sceta_) < 1.479) subdet = 1;
        else subdet = 2;
        Int_t ptBin = 0;
        if (leptonTree.probe_.Pt() > 20.0) ptBin = 1;
        Int_t MVABin = -1;
        if      (subdet == 0 && ptBin == 0) MVABin = 0;
        else if (subdet == 1 && ptBin == 0) MVABin = 1;
        else if (subdet == 2 && ptBin == 0) MVABin = 2;
        else if (subdet == 0 && ptBin == 1) MVABin = 3;
        else if (subdet == 1 && ptBin == 1) MVABin = 4;
        else if (subdet == 2 && ptBin == 1) MVABin = 5;

        // define cut on MVA value
        Double_t MVACut = -999.;
        if      (MVABin == 0) MVACut = 0.4202;
        else if (MVABin == 1) MVACut = 0.6206;
        else if (MVABin == 2) MVACut = 0.6190;
        else if (MVABin == 3) MVACut = 0.9590;
        else if (MVABin == 4) MVACut = 0.9586;
        else if (MVABin == 5) MVACut = 0.9278;

        if (leptonTree.electronHWW2011IDIsoMVA_ < MVACut) return false;
    }

    if (option == NumeratorElePOG2012IDIso) {

        float egammaPOGNeutral = TMath::Max(float(0.0), leptonTree.el_pfnhiso04_ + leptonTree.el_pfemiso04_ - leptonTree.el_ea04_ * TMath::Max(float(0.0), leptonTree.rhoIsoAll_));
        float egammaPOGIso = (leptonTree.el_pfchiso04_ + egammaPOGNeutral) / leptonTree.probe_.Pt();

        float eta = fabs(leptonTree.probe_.Eta());
        unsigned int bin = 0;
        if (leptonTree.probe_.Pt() > 20.0) {
            if      (eta < 0.8)                 bin = 0;
            else if      (eta >= 0.8 && eta < 1.479) bin = 1;
            else                                bin = 2;
        } else {
            if      (eta < 0.8)                 bin = 3;
            else if      (eta >= 0.8 && eta < 1.479) bin = 4;
            else                                bin = 5;
        }

        if (bin == 2 || bin == 5) {
            if (egammaPOGIso > 0.10)                                return false;
        } else {
            if (egammaPOGIso > 0.15)                                return false;
        }

        if      (bin == 0 && leptonTree.egammaPOG2012MVA_ < 0.94)    return false;
        if      (bin == 1 && leptonTree.egammaPOG2012MVA_ < 0.85)    return false;
        if      (bin == 2 && leptonTree.egammaPOG2012MVA_ < 0.92)    return false;
        if      (bin == 3 && leptonTree.egammaPOG2012MVA_ < 0.00)    return false;
        if      (bin == 4 && leptonTree.egammaPOG2012MVA_ < 0.10)    return false;
        if      (bin == 5 && leptonTree.egammaPOG2012MVA_ < 0.62)    return false;

        return true;


    }


    if (option == NumeratorMuHWW2011ID) {
        if ((leptonTree.leptonSelection_ & LeptonTree::PassMuID) != LeptonTree::PassMuID)     return false;
    }
    if (option == NumeratorMuHWW2011IDPlusPFMu) {
        if ((leptonTree.leptonSelection_ & LeptonTree::PassMuID) != LeptonTree::PassMuID)       return false;
        if ((leptonTree.leptonSelection_ & LeptonTree::PassMuIsPF) != LeptonTree::PassMuIsPF)   return false;
    }

    if (option == NumeratorMuHWW2011Iso) {
        if ((leptonTree.leptonSelection_ & LeptonTree::PassMuIso) != LeptonTree::PassMuIso)   return false;
    }
    if (option == NumeratorMuHWW2011IDIso) {
        if ((leptonTree.leptonSelection_ & LeptonTree::PassMuID) != LeptonTree::PassMuID)   return false;
        if ((leptonTree.leptonSelection_ & LeptonTree::PassMuIso) != LeptonTree::PassMuIso)   return false;
    }

    return true;
}







