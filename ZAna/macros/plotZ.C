//================================================================================================
//
// Z candidate categories:
//   1. mu-mu-2HLT  : 2 isolated global muons both matched to HLT triggering muon
//   2. mu-mu-1HLT  : 2 isolated global muons only one compatible with an HLT muon
//   3. mu-sa       : 1 isolated global muon + 1 standalone muon
//   4. mu-tk       : 1 isolated global muon + 1 tracker muon
//   5. mu-mu-noIso : 2 global muons, at least one compatible with HLT, at least one FAILS isolation
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
#include <TCut.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include "RooGlobalFunc.h"
#endif

#include "CPlot.hh"
#include "MitStyleRemix.hh"

#include "ZAnaStructDefs.hh"     // define structures to read in ntuple

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooPolynomial.h"
#include "RooPlot.h"

void plotZ() {
  
  gBenchmark->Start("plotZ");
  gSystem->Load("libRFIO");
  using namespace RooFit ;
  
  Int_t category = 0;
  TString title = "";
  if(category==1)      title = "mu-mu-2HLT";
  else if(category==2) title = "mu-mu-1HLT";
  else if(category==3) title = "mu-sa";
  else if(category==4) title = "mu-tk";
  else if(category==5) title = "mu-mu-noIso";
  
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================
  
  Bool_t doSave  = true;   // save plots?
  CPlot::sOutDir = ".";    // output directory
  TString format = "png";  // output file format
  Double_t lumi  = 10;     // luminosity normalization [pb^-1]
  
  // 
  // Values for data samples (signal sample should be first!)
  //   * For cross section values, look at /UserCode/MitPhysics/data/xs.dat in CVS
  //
  vector<TString>  fnamev;   // file names
  vector<Double_t> xsecv;    // cross section 
  vector<TString>  labelv;   // legend label
  vector<Int_t>    colorv;   // color in plots
  vector<Double_t> weightv;  // sample weight (to be computed)
  
  // signal sample
  //---------------
  fnamev.push_back("/castor/cern.ch/user/k/ksung/ntuples/00/s09-zmm-7-mc3_ntuple.root"); xsecv.push_back(1606.6);
//  fnamev.push_back("/castor/cern.ch/user/k/ksung/ntuples/00/s09-zmm-mc3_ntuple.root"); xsecv.push_back(2323.6);
  labelv.push_back("Z #rightarrow #mu#mu");
  colorv.push_back(kBlue);  
  
  // background samples
  //--------------------
  fnamev.push_back("/castor/cern.ch/user/k/ksung/ntuples/00/s09-wm-7-mc3_ntuple.root"); xsecv.push_back(9679.9*0.742);
//  fnamev.push_back("/castor/cern.ch/user/k/ksung/ntuples/00/s09-wm-mc3_ntuple.root"); xsecv.push_back(14253.7*0.691);
  labelv.push_back("W #rightarrow #mu#nu");
  colorv.push_back(kOrange+7);  
  
  fnamev.push_back("/castor/cern.ch/user/k/ksung/ntuples/00/s09-ttbar-7-mc3_ntuple.root"); xsecv.push_back(165.0);
//  fnamev.push_back("/castor/cern.ch/user/k/ksung/ntuples/00/s09-ttbar-mc3_ntuple.root"); xsecv.push_back(415.0);
  labelv.push_back("t#bar{t}");
  colorv.push_back(kGreen+2);  
  
  fnamev.push_back("/castor/cern.ch/user/k/ksung/ntuples/00/s09-ztt-7-mc3_ntuple.root"); xsecv.push_back(1606.6);
//  fnamev.push_back("/castor/cern.ch/user/k/ksung/ntuples/00/s09-ztt-mc3_ntuple.root"); xsecv.push_back(2323.6);
  labelv.push_back("Z #rightarrow #tau#tau");     
  colorv.push_back(kMagenta+2);
  
  fnamev.push_back("/castor/cern.ch/user/k/ksung/ntuples/00/s09-incmu15-7-mc3_ntuple.root"); xsecv.push_back(0.2969*0.00037*1e+09);
//  fnamev.push_back("/castor/cern.ch/user/k/ksung/ntuples/00/s09-incmu15-mc3_ntuple.root"); xsecv.push_back(0.5091*0.0002881*1000000000.);
  labelv.push_back("QCD");    
  colorv.push_back(kRed);  
  
  //
  // Cuts and requirements
  //
  const Double_t muPtCut   = 20;
  const Double_t muEtaCut  = 2.5;
  const Double_t massMin   = 20;
  const Double_t trkIsoCut = 3;
  //const Double_t d0Cut     = 0.2;
  
  //
  // Canvas dimensions
  //
  Int_t canw=800, canh=600;


  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================
    
  //
  // Set up histograms
  //  
  TH1F *hZMass = new TH1F("hZMass","",180,20,200);  // mass (sig + bkg)
  hZMass->Sumw2();
  
  // mass histogram for each sample
  vector <TH1F*> hZMassv;
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
    Int_t nbins;
    Double_t xlow, xhigh;
    char hname[100];
    
    nbins = hZMass->GetNbinsX();
    xlow  = hZMass->GetXaxis()->GetBinLowEdge(1);
    xhigh = hZMass->GetXaxis()->GetBinUpEdge(nbins);
    sprintf(hname,"%s_%i",hZMass->GetName(),ifile);
    hZMassv.push_back(new TH1F(hname,"",nbins,xlow,xhigh)); 
    hZMassv[ifile]->Sumw2();
  }
 
  Double_t brMass, brWeight;
  TTree massTree("massTree","Mass Tree");
  massTree.SetDirectory(0);  // force tree to be memory-resident
  massTree.Branch("mass",&brMass,"mass/D");
  massTree.Branch("weight",&brWeight,"weight/D");
      
  //
  // Fill histograms
  //
  TFile *infile=0;
  TTree *infoTree=0, *genTree=0, *dimuonTree=0, *muonTree=0;
    
  // Data structures to store info from TTrees
  TEventInfo eventInfo;
  TGenInfo   genInfo;
  TDimuon    dimuon;
  TMuon      muon;
  
  // loop through samples
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
            
      if(eventInfo.HLT_Mu15==0) continue;                                           // trigger accept                                      
      if(dimuon.mass < massMin) continue;                                           // mass region
      if(dimuon.pt_1 < muPtCut || dimuon.pt_2 < muPtCut) continue;                  // muon pT cut
      if(fabs(dimuon.eta_1) > muEtaCut || fabs(dimuon.eta_2) > muEtaCut) continue;  // muon eta cut
      
      if(category==1) {
        if((dimuon.type_1 != kGlobal) || (dimuon.type_2 != kGlobal))           continue;  // 2 global muons
	if((dimuon.trkIso03_1 > trkIsoCut) || (dimuon.trkIso03_2 > trkIsoCut)) continue;  // track isolation
	if(dimuon.pt_1 < 15)                                                   continue;  // satisfy HLT requirements 
      } else if(category==2) {
        if((dimuon.type_1 != kGlobal) || (dimuon.type_2 != kGlobal))           continue;  // 2 global muons
	if((dimuon.trkIso03_1 > trkIsoCut) || (dimuon.trkIso03_2 > trkIsoCut)) continue;  // track isolation
	if((dimuon.pt_1 < 15) || (dimuon.pt_2 > 15))                           continue;  // satisfy HLT requirements
      } else if(category==3) {
        if((dimuon.type_1 == kGlobal)                  // one global muon...
	   && (dimuon.trkIso03_1 < trkIsoCut)          // and isolated...
	   && (dimuon.pt_1 > 15)) {                    // and satisfy HLT requirements
	  if(dimuon.type_2 != kStandalone)             // one standalone muon
	    continue;
	} else if((dimuon.type_2 == kGlobal)           // one global muon...
	           && (dimuon.trkIso03_2 < trkIsoCut)  // and isolated...
		   && (dimuon.pt_2 > 15)) {            // and satisfy HLT requirements
	  if(dimuon.type_1 != kStandalone)             // one standalone muon
	    continue;
	} else {
	  continue;
	}
      } else if(category==4) {     
        if((dimuon.type_1 == kGlobal)                  // one global muon
	   && (dimuon.trkIso03_1 < trkIsoCut)          // and isolated...
	   && (dimuon.pt_1 > 15)) {                    // and satisfy HLT requirements
	  if(dimuon.type_2 != kTracker)                // one tracker muon
	    continue;
	} else if((dimuon.type_2 == kGlobal)           // one global muon
	           && (dimuon.trkIso03_2 < trkIsoCut)  // and isolated...
		   && (dimuon.pt_2 > 15)) {            // and satisfy HLT requirements
	  if(dimuon.type_1 != kTracker)                // one tracker muon
	    continue;
	} else {
	  continue;
	}
      } else if(category==5) { 
        if((dimuon.type_1 != kGlobal) || (dimuon.type_2 != kGlobal))           continue;  // 2 global muons
	if((dimuon.trkIso03_1 < trkIsoCut) && (dimuon.trkIso03_2 < trkIsoCut)) continue;  // at least one fails isolation
	if(dimuon.pt_1 < 15)                                                   continue;  // at least one satisfy HLT requirements
      } else {
        if((dimuon.trkIso03_1 > trkIsoCut) || (dimuon.trkIso03_2 > trkIsoCut)) continue;  // track isolation
        //if(fabs(dimuon.d0_1) > d0Cut || fabs(dimuon.d0_2) > d0Cut)             continue;  // d0 cut
      }
      
      hZMass->Fill(dimuon.mass,weightv[ifile]);
      hZMassv[ifile]->Fill(dimuon.mass,weightv[ifile]);
      
      brMass = dimuon.mass;
      brWeight = weightv[ifile];
      massTree.Fill();
    }

    delete infile;
    infile=0, infoTree=0, genTree=0, dimuonTree=0, muonTree=0; 
  }
  
  //
  // Make plots
  //
  TCanvas *c = MakeCanvas("c","c",canw,canh);
  char name[50];
  
  sprintf(name,"zmass%i",category);
  CPlot plotZMass(name,title,"m(#mu#mu) [GeV/c^{2}]","events / 1.0 GeV/c^{2}");
  plotZMass.AddHist1D(hZMass,"Sig + Bkg","E");
  for(UInt_t ibg=1; ibg<hZMassv.size(); ibg++)
    plotZMass.AddToStack(hZMassv[ibg],labelv[ibg],colorv[ibg]);
  plotZMass.AddToStack(hZMassv[0],labelv[0],colorv[0]);
  plotZMass.SetLogy();
  plotZMass.AddTextBox("#int#font[12]{L}dt = 10 pb^{-1}",0.21,0.85,0.41,0.8,0);
  plotZMass.Draw(c,doSave,format);  

/* WORK IN PROGRESS :P
  //
  // RooFit
  //  
  RooRealVar mass("mass","dimuon mass",40.,200.);
  RooRealVar weight("weight","event weight",0,1000);
  RooRealVar sigMean("sigMean","signal mean",91.,70.,110.);
  RooRealVar sigWidth("sigWidth","signal width",20.,0.,50.);
  RooGaussian signal("signal","signal",mass,sigMean,sigWidth);

//  RooRealVar nsig("nsig","signal events",2000,0,3000);
//  RooRealVar nbkg("nbkg","background events",40,0,3000);
//  RooAddPdf model("model","signal + background",RooArgList(signal,background),RooArgList(nsig,nbkg));
  
//  RooDataHist data("data","MC data",mass,hZMass);
RooDataSet data("data","MC data",&massTree,mass);
//  RooDataSet data("data","MC data",&massTree,RooArgSet(mass,weight),0,"weight");
  signal.fitTo(data);
  RooPlot *frame = mass.frame();
  data.plotOn(frame);
  signal.plotOn(frame,LineColor(kRed));
  sigMean.Print();
  sigWidth.Print();
  frame->Draw();
*/  

  //--------------------------------------------------------------------------------------------------------------
  // Summary print out
  //==============================================================================================================
  cout << endl;
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  cout << "  Luminosity = " << lumi << " / pb" << endl;
  cout << endl;
  
  cout << "  Signal file: " << fnamev[0] << endl;
  for(UInt_t ibg=1; ibg<hZMassv.size(); ibg++) {
    cout << "  Bkgd " << ibg << " file: " << fnamev[ibg] << endl;
  }
  cout << endl;
  
  cout << setw(20) << "Total events =";
  cout << setw(8) << setprecision(1) << fixed << hZMass->Integral() << endl;
  
  cout << setw(20) << "Signal =";
  cout << setw(8) << setprecision(1) << fixed << hZMassv[0]->Integral();
  cout << " \u00B1 "; 
  cout << setw(5) << setprecision(2) << sqrt(hZMassv[0]->GetEffectiveEntries())*weightv[0] << endl;; 
  
  for(UInt_t ibg=1; ibg<hZMassv.size(); ibg++) {    
    cout << setw(17) << "Bkgd " << ibg << " =";
    cout << setw(8) << setprecision(1) << fixed << hZMassv[ibg]->Integral();
    cout << " \u00B1 ";
    cout << setw(5) << setprecision(2) << sqrt(hZMassv[ibg]->GetEffectiveEntries())*weightv[ibg] << endl;     
  }
  cout << endl;
  
  gBenchmark->Show("plotZ");
}
