//================================================================================================
//
// Tag & Probe analysis of isolation efficiency
//
// Use the macro plotIsoEffMC.C to get isolation efficiency from Monte Carlo...
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <TFile.h>
#include <TRFIOFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TVirtualFitter.h>
#include <TBenchmark.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <vector>
#include <iostream>
#include <iomanip>
#endif

#include "CPlot.hh"
#include "MitStyleRemix.hh"

#include "ZAnaStructDefs.hh"     // define structures to read in ntuple

//#define USE_ROOSTATSCMS
#include "MyTools.hh"

void plotIsoEffTP() {
  
  gBenchmark->Start("plotIsoEffTP");
  gSystem->Load("libRFIO");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================
  
  Bool_t doSave  = true;  // save plots?
  CPlot::sOutDir = ".";    // output directory
  TString format = "png";  // output file format
  Double_t lumi  = 10;     // luminosity normalization [pb^-1]
  Int_t method   = 0;      // efficiency errors calculation method
  
  // 
  // Values for data samples (signal sample should be first!)
  //   * For cross section values, look at /UserCode/MitPhysics/data/xs.dat in CVS
  //
  vector<TString>  fnamev;   // file names
  vector<Double_t> xsecv;    // cross section 
  vector<Double_t> weightv;  // sample weight (to be computed)
  
  // signal sample
  //---------------
  fnamev.push_back("/castor/cern.ch/user/k/ksung/ntuples/00/s09-zmm-7-mc3_ntuple.root"); xsecv.push_back(1606.6);
//  fnamev.push_back("/castor/cern.ch/user/k/ksung/ntuples/00/s09-zmm-mc3_ntuple.root"); xsecv.push_back(2323.6);
  
  // background samples
  //--------------------
  fnamev.push_back("/castor/cern.ch/user/k/ksung/ntuples/00/s09-wm-7-mc3_ntuple.root"); xsecv.push_back(9679.9*0.742);
//  fnamev.push_back("/castor/cern.ch/user/k/ksung/ntuples/00/s09-wm-mc3_ntuple.root"); xsecv.push_back(14253.7*0.691);
  
  fnamev.push_back("/castor/cern.ch/user/k/ksung/ntuples/00/s09-ttbar-7-mc3_ntuple.root"); xsecv.push_back(165.0);
//  fnamev.push_back("/castor/cern.ch/user/k/ksung/ntuples/00/s09-ttbar-mc3_ntuple.root"); xsecv.push_back(415.0);
  
  fnamev.push_back("/castor/cern.ch/user/k/ksung/ntuples/00/s09-ztt-7-mc3_ntuple.root"); xsecv.push_back(1606.6);
//  fnamev.push_back("/castor/cern.ch/user/k/ksung/ntuples/00/s09-ztt-mc3_ntuple.root"); xsecv.push_back(2323.6);
  
  fnamev.push_back("/castor/cern.ch/user/k/ksung/ntuples/00/s09-incmu15-7-mc3_ntuple.root"); xsecv.push_back(0.2969*0.00037*1e+09);
//  fnamev.push_back("/castor/cern.ch/user/k/ksung/ntuples/00/s09-incmu15-mc3_ntuple.root"); xsecv.push_back(0.5091*0.0002881*1000000000.);
  
  //
  // Cuts and requirements
  //
  //  Preselection: 70 < mass < 110
  //      TAG => isolated global muon, pT > 20, |eta| < 2.0 
  //    PROBE => global muon, pT > 20, |eta| < 2.0
  //
  const Double_t tagPtCut    = 20;
  const Double_t tagEtaCut   = 2;
  const UInt_t   tagType     = kGlobal;
  const Double_t probePtCut  = 20;
  const Double_t probeEtaCut = 2;
  const UInt_t   probeType   = kGlobal;
  const Double_t massMin     = 70;
  const Double_t massMax     = 110;
  const Double_t trkIsoCut   = 3;
  
  //
  // Canvas dimensions
  //
  Int_t canw=800, canh=600;

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================
  
  //
  // Variables for efficiency calculations
  //
  Double_t numer=0, denom=0;
  
  Double_t ptBinEdges[] = { 20, 24, 26, 28, 30, 32, 34, 36, 38, 40, 
                            44, 48, 52, 60, 70, 90 };
  const UInt_t nPtBins = sizeof(ptBinEdges)/sizeof(Double_t) - 1;  
  Double_t numerPt[nPtBins];
  Double_t denomPt[nPtBins];
  for(UInt_t ibin=0; ibin<nPtBins; ibin++) { numerPt[ibin] = denomPt[ibin] = 0; }
  
  Double_t etaBinEdges[21];
  const UInt_t nEtaBins = sizeof(etaBinEdges)/sizeof(Double_t) - 1; 
  for(UInt_t i=0; i<nEtaBins+1; i++) { etaBinEdges[i] = -2.0+i*0.2; }
  Double_t numerEta[nEtaBins];
  Double_t denomEta[nEtaBins];
  for(UInt_t ibin=0; ibin<nEtaBins; ibin++) { numerEta[ibin] = denomEta[ibin] = 0; }
  
  Double_t phiBinEdges[17];
  const UInt_t nPhiBins = sizeof(phiBinEdges)/sizeof(Double_t) - 1;
  for(UInt_t i=0; i<nPhiBins+1; i++) { phiBinEdges[i] = -3.2+i*0.4; }
  Double_t numerPhi[nPhiBins];
  Double_t denomPhi[nPhiBins];
  for(UInt_t ibin=0; ibin<nPhiBins; ibin++) { numerPhi[ibin] = denomPhi[ibin] = 0; }
  
  TFile *infile=0;
  TTree *infoTree=0, *genTree=0, *dimuonTree=0, *muonTree=0;
    
  // Data structures to store info from TTrees
  TEventInfo eventInfo;
  TGenInfo   genInfo;
  TDimuon    dimuon;
  TMuon      muon;

  //
  // Loop through samples
  //
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {  
  
    // Read input file and get the TTrees
    if(fnamev[ifile].Contains("/castor/cern.ch",TString::kIgnoreCase))
      infile = new TRFIOFile(fnamev[ifile]);
    else
      infile = new TFile(fnamev[ifile]); 
    assert(infile);
    
    infoTree   = (TTree*)infile->FindObjectAny("EventInfo"); assert(infoTree);
    genTree    = (TTree*)infile->FindObjectAny("GenInfo");   assert(genTree);
    dimuonTree = (TTree*)infile->FindObjectAny("Dimuon");    assert(dimuonTree);
    muonTree   = (TTree*)infile->FindObjectAny("Muon");      assert(muonTree);

    // Set branch address to structures that will store the info  
    infoTree  ->SetBranchAddress("EventInfo",&eventInfo);
    genTree   ->SetBranchAddress("GenInfo",  &genInfo);
    dimuonTree->SetBranchAddress("Dimuon",   &dimuon);
    muonTree  ->SetBranchAddress("Muon",     &muon);
    
    weightv.push_back(lumi*xsecv[ifile]/(Double_t)infoTree->GetEntries());  
  
    // loop over events
    for(UInt_t ientry=0; ientry<dimuonTree->GetEntries(); ientry++) {
      dimuonTree->GetEntry(ientry);
      
      assert(dimuon.eventId <= infoTree->GetEntries());
      infoTree->GetEntry(dimuon.eventId);
            
      if(eventInfo.HLT_Mu15==0)                              continue;  // trigger accept                                      
      if((dimuon.mass < massMin) || (dimuon.mass > massMax)) continue;  // mass region
      if(dimuon.q_1 == dimuon.q_2)                           continue;  // require opposite charge 
      
      // muon 1 is tag
      if((dimuon.pt_1 > tagPtCut) && (fabs(dimuon.eta_1) < tagEtaCut) && (dimuon.trkIso03_1 < trkIsoCut) && (dimuon.type_1 == tagType)) {
	// muon 2 is probe
	if((dimuon.pt_2 > probePtCut) && (fabs(dimuon.eta_2) < probeEtaCut) && (dimuon.type_2 == probeType)) {	  
	  Bool_t pass = dimuon.trkIso03_2 < trkIsoCut;
	  denom += weightv[ifile];
	  if(pass) { numer += weightv[ifile]; }
	  for(UInt_t ibin=0; ibin<nPtBins; ibin++) { 
	    if((dimuon.pt_2 > ptBinEdges[ibin]) && (dimuon.pt_2 < ptBinEdges[ibin+1])) {
	      denomPt[ibin] += weightv[ifile];  
	      if(pass) { numerPt[ibin] += weightv[ifile]; }
	      break; 
	    }
	  }
	  for(UInt_t ibin=0; ibin<nEtaBins; ibin++) { 
	    if((dimuon.eta_2 > etaBinEdges[ibin]) && (dimuon.eta_2 < etaBinEdges[ibin+1])) {
	      denomEta[ibin] += weightv[ifile]; 
	      if(pass) { numerEta[ibin] += weightv[ifile]; } 
	      break;
	    }
	  }
	  for(UInt_t ibin=0; ibin<nPhiBins; ibin++) { 
	    if((dimuon.phi_2 > phiBinEdges[ibin]) && (dimuon.phi_2 < phiBinEdges[ibin+1])) {
	      denomPhi[ibin] += weightv[ifile]; 
	      if(pass) { numerPhi[ibin] += weightv[ifile]; } 
	      break;
	    }
	  }
	}
      }
      
      // muon 2 is tag
      if((dimuon.pt_2 > tagPtCut) && (fabs(dimuon.eta_2) < tagEtaCut) && (dimuon.trkIso03_2 < trkIsoCut) && (dimuon.type_2 == tagType)) {
	// muon 1 is probe
	if((dimuon.pt_1 > probePtCut) && (fabs(dimuon.eta_1) < probeEtaCut) && (dimuon.type_1 == probeType)) {
	  Bool_t pass = dimuon.trkIso03_1 < trkIsoCut;
	  denom += weightv[ifile];
	  if(pass) { numer += weightv[ifile]; }
	  for(UInt_t ibin=0; ibin<nPtBins; ibin++) { 
	    if((dimuon.pt_1 > ptBinEdges[ibin]) && (dimuon.pt_1 < ptBinEdges[ibin+1])) {
	      denomPt[ibin] += weightv[ifile];  
	      if(pass) { numerPt[ibin] += weightv[ifile]; }
	      break; 
	    }
	  }
	  for(UInt_t ibin=0; ibin<nEtaBins; ibin++) { 
	    if((dimuon.eta_1 > etaBinEdges[ibin]) && (dimuon.eta_1 < etaBinEdges[ibin+1])) {
	      denomEta[ibin] += weightv[ifile]; 
	      if(pass) { numerEta[ibin] += weightv[ifile]; } 
	      break;
	    }
	  }
	  for(UInt_t ibin=0; ibin<nPhiBins; ibin++) { 
	    if((dimuon.phi_1 > phiBinEdges[ibin]) && (dimuon.phi_1 < phiBinEdges[ibin+1])) {
	      denomPhi[ibin] += weightv[ifile]; 
	      if(pass) { numerPhi[ibin] += weightv[ifile]; } 
	      break;
	    }
	  }
	}
      }
    }

    delete infile;
    infile=0, infoTree=0, genTree=0, dimuonTree=0, muonTree=0;    
  }
  
  //--------------------------------------------------------------------------------------------------------------
  // Make plots 
  //==============================================================================================================  
  TCanvas *c = MakeCanvas("c","c",canw,canh);
  char label[50];
  sprintf(label,"Tag & Probe (%i pb^{-1})",(Int_t)lumi);

  //
  // pT dependence
  //
   TGraphAsymmErrors *grEffPtMC = new TGraphAsymmErrors(15);  // copy and paste from iso03EffvPtSig.C
   grEffPtMC->SetName("iso03EffvPtSig_gr_0");
   grEffPtMC->SetTitle("#DeltaR < 0.3");
   grEffPtMC->SetFillColor(1);
   grEffPtMC->SetLineWidth(2);
   grEffPtMC->SetMarkerStyle(8);
   grEffPtMC->SetMarkerSize(0.9);
   grEffPtMC->SetPoint(0,22,0.9627918);
   grEffPtMC->SetPointError(0,2,2,0.00133598,0.001305926);
   grEffPtMC->SetPoint(1,25,0.9631563);
   grEffPtMC->SetPointError(1,1,1,0.001663158,0.001616369);
   grEffPtMC->SetPoint(2,27,0.9662244);
   grEffPtMC->SetPointError(2,1,1,0.001500864,0.001459145);
   grEffPtMC->SetPoint(3,29,0.9672292);
   grEffPtMC->SetPointError(3,1,1,0.001372564,0.001336507);
   grEffPtMC->SetPoint(4,31,0.9673415);
   grEffPtMC->SetPointError(4,1,1,0.001292708,0.001260563);
   grEffPtMC->SetPoint(5,33,0.9685371);
   grEffPtMC->SetPointError(5,1,1,0.001195673,0.001167064);
   grEffPtMC->SetPoint(6,35,0.9703608);
   grEffPtMC->SetPointError(6,1,1,0.001092654,0.001067224);
   grEffPtMC->SetPoint(7,37,0.9711445);
   grEffPtMC->SetPointError(7,1,1,0.001016808,0.0009941471);
   grEffPtMC->SetPoint(8,39,0.9763942);
   grEffPtMC->SetPointError(8,1,1,0.0008652992,0.0008451448);
   grEffPtMC->SetPoint(9,42,0.9771184);
   grEffPtMC->SetPointError(9,2,2,0.0005590694,0.0005503147);
   grEffPtMC->SetPoint(10,46,0.9795804);
   grEffPtMC->SetPointError(10,2,2,0.0005786637,0.0005681531);
   grEffPtMC->SetPoint(11,50,0.9805051);
   grEffPtMC->SetPointError(11,2,2,0.0008445414,0.0008212922);
   grEffPtMC->SetPoint(12,56,0.9791224);
   grEffPtMC->SetPointError(12,4,4,0.0009763667,0.0009474545);
   grEffPtMC->SetPoint(13,65,0.9779054);
   grEffPtMC->SetPointError(13,5,5,0.001482698,0.001420568);
   grEffPtMC->SetPoint(14,80,0.9770187);
   grEffPtMC->SetPointError(14,10,10,0.001932364,0.00183205);

  TGraphAsymmErrors *grEffPt = toolbox::makeEffGraph(nPtBins, ptBinEdges, numerPt, denomPt, method);
  CPlot plotIsoEffvPt("isoEffvPt","","p_{T} [GeV/c]","isolation efficiency");
  plotIsoEffvPt.AddGraph(grEffPt,label,"");
  plotIsoEffvPt.AddGraph(grEffPtMC,"MC","",kRed,kFullTriangleDown);
  plotIsoEffvPt.SetYRange(0.92,1.01);
  plotIsoEffvPt.Draw(c,doSave,format);

  //
  // eta dependence
  //
   TGraphAsymmErrors *grEffEtaMC = new TGraphAsymmErrors(20);  // copy and paste from iso03EffvEtaSig.C
   grEffEtaMC->SetPoint(0,-1.9,0.9770189);
   grEffEtaMC->SetPointError(0,0.1,0.1,0.001287716,0.001242391);
   grEffEtaMC->SetPoint(1,-1.7,0.9766598);
   grEffEtaMC->SetPointError(1,0.1,0.1,0.001237384,0.001196115);
   grEffEtaMC->SetPoint(2,-1.5,0.9744966);
   grEffEtaMC->SetPointError(2,0.1,0.1,0.001244215,0.001206013);
   grEffEtaMC->SetPoint(3,-1.3,0.9758192);
   grEffEtaMC->SetPointError(3,0.1,0.1,0.001150529,0.001116003);
   grEffEtaMC->SetPoint(4,-1.1,0.9749559);
   grEffEtaMC->SetPointError(4,0.1,0.1,0.00112609,0.001094131);
   grEffEtaMC->SetPoint(5,-0.9,0.972916);
   grEffEtaMC->SetPointError(5,0.1,0.1,0.001139904,0.001109633);
   grEffEtaMC->SetPoint(6,-0.7,0.9732467);
   grEffEtaMC->SetPointError(6,0.1,0.1,0.001102692,0.00107399);
   grEffEtaMC->SetPoint(7,-0.5,0.9719888);
   grEffEtaMC->SetPointError(7,0.1,0.1,0.001106241,0.001078657);
   grEffEtaMC->SetPoint(8,-0.3,0.9716777);
   grEffEtaMC->SetPointError(8,0.1,0.1,0.001104444,0.001077253);
   grEffEtaMC->SetPoint(9,-0.1,0.9722481);
   grEffEtaMC->SetPointError(9,0.1,0.1,0.001083369,0.001056653);
   grEffEtaMC->SetPoint(10,0.1,0.9724435);
   grEffEtaMC->SetPointError(10,0.1,0.1,0.001082571,0.001055703);
   grEffEtaMC->SetPoint(11,0.3,0.9728172);
   grEffEtaMC->SetPointError(11,0.1,0.1,0.001080864,0.001053711);
   grEffEtaMC->SetPoint(12,0.5,0.9722489);
   grEffEtaMC->SetPointError(12,0.1,0.1,0.001099607,0.001072093);
   grEffEtaMC->SetPoint(13,0.7,0.9719029);
   grEffEtaMC->SetPointError(13,0.1,0.1,0.001131266,0.001102526);
   grEffEtaMC->SetPoint(14,0.9,0.9746468);
   grEffEtaMC->SetPointError(14,0.1,0.1,0.001100085,0.001069938);
   grEffEtaMC->SetPoint(15,1.1,0.9758824);
   grEffEtaMC->SetPointError(15,0.1,0.1,0.001105469,0.001073473);
   grEffEtaMC->SetPoint(16,1.3,0.9734979);
   grEffEtaMC->SetPointError(16,0.1,0.1,0.001208095,0.001173404);
   grEffEtaMC->SetPoint(17,1.5,0.9754444);
   grEffEtaMC->SetPointError(17,0.1,0.1,0.00121324,0.001175492);
   grEffEtaMC->SetPoint(18,1.7,0.976848);
   grEffEtaMC->SetPointError(18,0.1,0.1,0.001222359,0.001181744);
   grEffEtaMC->SetPoint(19,1.9,0.9760142);
   grEffEtaMC->SetPointError(19,0.1,0.1,0.001314578,0.001269332);
  
  TGraphAsymmErrors *grEffEta = toolbox::makeEffGraph(nEtaBins, etaBinEdges, numerEta, denomEta, method);
  CPlot plotIsoEffvEta("isoEffvEta","","#eta","isolation efficiency");
  plotIsoEffvEta.AddGraph(grEffEta,label,"");
  plotIsoEffvEta.AddGraph(grEffEtaMC,"MC Z#rightarrow#mu#mu","",kRed,kFullTriangleDown);
  plotIsoEffvEta.SetYRange(0.92,1.01);
  plotIsoEffvEta.Draw(c,doSave,format);
  
  //
  // phi dependence
  //  
   TGraphAsymmErrors *grEffPhiMC = new TGraphAsymmErrors(32);  // copy and paste from iso03EffvPhiSig.C
   grEffPhiMC->SetPoint(0,-3.1,0.9766281);
   grEffPhiMC->SetPointError(0,0.1,0.1,0.00164703,0.001574799);
   grEffPhiMC->SetPoint(1,-2.9,0.9716119);
   grEffPhiMC->SetPointError(1,0.1,0.1,0.001513758,0.001463252);
   grEffPhiMC->SetPoint(2,-2.7,0.9746968);
   grEffPhiMC->SetPointError(2,0.1,0.1,0.001434234,0.001383312);
   grEffPhiMC->SetPoint(3,-2.5,0.9743509);
   grEffPhiMC->SetPointError(3,0.1,0.1,0.001430623,0.001380631);
   grEffPhiMC->SetPoint(4,-2.3,0.9732208);
   grEffPhiMC->SetPointError(4,0.1,0.1,0.001481182,0.001429899);
   grEffPhiMC->SetPoint(5,-2.1,0.9745715);
   grEffPhiMC->SetPointError(5,0.1,0.1,0.001438912,0.001387913);
   grEffPhiMC->SetPoint(6,-1.9,0.9737388);
   grEffPhiMC->SetPointError(6,0.1,0.1,0.001455186,0.001404687);
   grEffPhiMC->SetPoint(6,-1.9,0.9737388);
   grEffPhiMC->SetPointError(6,0.1,0.1,0.001455186,0.001404687);
   grEffPhiMC->SetPoint(7,-1.7,0.9733765);
   grEffPhiMC->SetPointError(7,0.1,0.1,0.001457048,0.001407107);
   grEffPhiMC->SetPoint(8,-1.5,0.9762056);
   grEffPhiMC->SetPointError(8,0.1,0.1,0.001383781,0.001333339);
   grEffPhiMC->SetPoint(9,-1.3,0.9744484);
   grEffPhiMC->SetPointError(9,0.1,0.1,0.001448116,0.00139672);
   grEffPhiMC->SetPoint(10,-1.1,0.9736108);
   grEffPhiMC->SetPointError(10,0.1,0.1,0.001457676,0.001407252);
   grEffPhiMC->SetPoint(11,-0.9,0.9700821);
   grEffPhiMC->SetPointError(11,0.1,0.1,0.001539649,0.001490097);
   grEffPhiMC->SetPoint(12,-0.7,0.9749822);
   grEffPhiMC->SetPointError(12,0.1,0.1,0.001413711,0.00136365);
   grEffPhiMC->SetPoint(13,-0.5,0.9741379);
   grEffPhiMC->SetPointError(13,0.1,0.1,0.001431171,0.001381547);
   grEffPhiMC->SetPoint(14,-0.3,0.974962);
   grEffPhiMC->SetPointError(14,0.1,0.1,0.001424016,0.001373277);
   grEffPhiMC->SetPoint(15,-0.1,0.977269);
   grEffPhiMC->SetPointError(15,0.1,0.1,0.001372532,0.0013206);
   grEffPhiMC->SetPoint(16,0.1,0.9751513);
   grEffPhiMC->SetPointError(16,0.1,0.1,0.001425035,0.001373841);
   grEffPhiMC->SetPoint(17,0.3,0.9729943);
   grEffPhiMC->SetPointError(17,0.1,0.1,0.001466627,0.00141675);
   grEffPhiMC->SetPoint(18,0.5,0.9750139);
   grEffPhiMC->SetPointError(18,0.1,0.1,0.001418797,0.001368319);
   grEffPhiMC->SetPoint(19,0.7,0.9707513);
   grEffPhiMC->SetPointError(19,0.1,0.1,0.001536987,0.001486475);
   grEffPhiMC->SetPoint(20,0.9,0.974361);
   grEffPhiMC->SetPointError(20,0.1,0.1,0.001432312,0.001382182);
   grEffPhiMC->SetPoint(21,1.1,0.9742814);
   grEffPhiMC->SetPointError(21,0.1,0.1,0.001438959,0.001388527);
   grEffPhiMC->SetPoint(22,1.3,0.9728129);
   grEffPhiMC->SetPointError(22,0.1,0.1,0.001491951,0.001440709);
   grEffPhiMC->SetPoint(23,1.5,0.9770097);
   grEffPhiMC->SetPointError(23,0.1,0.1,0.001365718,0.001314862);
   grEffPhiMC->SetPoint(24,1.7,0.9733514);
   grEffPhiMC->SetPointError(24,0.1,0.1,0.001462829,0.001412544);
   grEffPhiMC->SetPoint(25,1.9,0.9729664);
   grEffPhiMC->SetPointError(25,0.1,0.1,0.001483648,0.001432678);
   grEffPhiMC->SetPoint(26,2.1,0.9717292);
   grEffPhiMC->SetPointError(26,0.1,0.1,0.001494765,0.001445296);
   grEffPhiMC->SetPoint(27,2.3,0.9735963);
   grEffPhiMC->SetPointError(27,0.1,0.1,0.001456227,0.001405928);
   grEffPhiMC->SetPoint(28,2.5,0.9746269);
   grEffPhiMC->SetPointError(28,0.1,0.1,0.001419897,0.001370106);
   grEffPhiMC->SetPoint(29,2.7,0.971483);
   grEffPhiMC->SetPointError(29,0.1,0.1,0.001509714,0.001459701);
   grEffPhiMC->SetPoint(30,2.9,0.9736262);
   grEffPhiMC->SetPointError(30,0.1,0.1,0.001447959,0.001398164);
   grEffPhiMC->SetPoint(31,3.1,0.973143);
   grEffPhiMC->SetPointError(31,0.1,0.1,0.001750891,0.001679881);

  TGraphAsymmErrors *grEffPhi = toolbox::makeEffGraph(nPhiBins, phiBinEdges, numerPhi, denomPhi, method);
  CPlot plotIsoEffvPhi("isoEffvPhi","","#phi","isolation efficiency");
  plotIsoEffvPhi.AddGraph(grEffPhi,label,"");
  plotIsoEffvPhi.AddGraph(grEffPhiMC,"MC","",kRed,kFullTriangleDown);
  plotIsoEffvPhi.SetYRange(0.92,1.01);
  plotIsoEffvPhi.Draw(c,doSave,format);  
    
  //--------------------------------------------------------------------------------------------------------------
  // Summary print out
  //==============================================================================================================
  cout << endl;
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  cout << "  Luminosity = " << lumi << " / pb" << endl;
  cout << endl;

  Double_t errl, errh;
  Double_t eff = toolbox::calcEff(numer,denom,&errl,&errh,method);  
  cout << "  Number of probes: " << setprecision(1) << fixed << denom << endl;
  cout << "  Average Efficiency = " << setprecision(4) << fixed << eff;
  cout << " [-" << setprecision(4) << fixed << errl;
  cout << ", +" << setprecision(4) << fixed << errh << "]" << endl;
  cout << endl;
  
  gBenchmark->Show("plotIsoEffTP");
}

