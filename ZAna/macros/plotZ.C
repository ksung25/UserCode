//================================================================================================
//
// Plots for reconstructed Z candidates
//
// Z candidate categories:
//   0. standard selection: 2 isolated global muons at least one matched to an HLT primitive
//   1. mu-mu-2HLT  : 2 isolated global muons both matched to HLT primitives
//   2. mu-mu-1HLT  : 2 isolated global muons only one compatible with an HLT primitive
//   3. mu-tk       : 1 isolated global muon + 1 isolated tracker muon
//   4. mu-sa       : 1 isolated global muon + 1 isolated standalone muon
//   5. mu-mu-noIso : 2 global muons, at least one compatible with HLT, at least one FAILS isolation
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>               // access to gROOT, entry point to ROOT system
#include <TFile.h>               // file handle class
#include <TRFIOFile.h>           // handles RFIO (read from CASTOR)
#include <TTree.h>               // class to access ntuples
#include <TCanvas.h>             // class for drawing
#include <TH1F.h>                // 1D histograms
#include <TF1.h>                 // 1D functions
#include <TVirtualFitter.h>      // abstract base class for fitting in ROOT
#include <TBenchmark.h>          // class to track macro running statistics
#include <vector>                // STL vector class
#include <iostream>              // standard I/O
#include <iomanip>               // functions to format standard I/O
#include "RooGlobalFunc.h"
#endif

#include "CPlot.hh"              // helper class for plots
#include "MitStyleRemix.hh"      // style settings for drawing

#include "ZAnaStructDefs.hh"     // define structures to read in ntuple

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooPolynomial.h"
#include "RooPlot.h"

void plotZ(Int_t category=0) 
{  
  gBenchmark->Start("plotZ");
  gSystem->Load("libRFIO");
  using namespace RooFit ;

    
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================
  
  Bool_t doSave  = true;      // save plots?
  CPlot::sOutDir = "ZPlots";  // output directory
  TString format = "png";     // output file format
  Double_t lumi  = 10;        // luminosity normalization [pb^-1]
  
  // 
  // Values for samples (data sample first, then MC signal sample, then MC background sample)
  //   * For cross section values, look at /UserCode/MitPhysics/data/xs.dat in CVS
  //
  vector<TString>  fnamev;   // file names
  vector<Double_t> xsecv;    // cross section 
  vector<TString>  labelv;   // legend label
  vector<Int_t>    colorv;   // color in plots
  vector<Double_t> weightv;  // sample weight (to be computed)
  
  // data (Use blank ("") name if no data sample to process)
  //---------------------------------------------------------
  fnamev.push_back("");  
  xsecv.push_back(0);
  labelv.push_back("data");
  colorv.push_back(kBlack);
  weightv.push_back(1);
  
  // signal MC sample
  //------------------
  fnamev.push_back("/desktop/08a/ksung/ZAna/s09-zmm-7-mc3_ntuple.root"); xsecv.push_back(1606.6);
//  fnamev.push_back("/desktop/08a/ksung/ZAna/s09-zmm-mc3_ntuple.root"); xsecv.push_back(2323.6);
  labelv.push_back("Z #rightarrow #mu#mu");
  colorv.push_back(kBlue);  
  
  // background MC samples
  //-----------------------
  fnamev.push_back("/desktop/08a/ksung/ZAna/s09-wm-7-mc3_ntuple.root"); xsecv.push_back(9679.9*0.742);
//  fnamev.push_back("/desktop/08a/ksung/ZAna/s09-wm-mc3_ntuple.root"); xsecv.push_back(14253.7*0.691);
  labelv.push_back("W #rightarrow #mu#nu");
  colorv.push_back(kOrange+7);  
  
  fnamev.push_back("/desktop/08a/ksung/ZAna/s09-ttbar-7-mc3_ntuple.root"); xsecv.push_back(165.0);
//  fnamev.push_back("/desktop/08a/ksung/ZAna/s09-ttbar-mc3_ntuple.root"); xsecv.push_back(415.0);
  labelv.push_back("t#bar{t}");
  colorv.push_back(kGreen+2);  
  
  fnamev.push_back("/desktop/08a/ksung/ZAna/s09-ztt-7-mc3_ntuple.root"); xsecv.push_back(1606.6);
//  fnamev.push_back("/desktop/08a/ksung/ZAna/s09-ztt-mc3_ntuple.root"); xsecv.push_back(2323.6);
  labelv.push_back("Z #rightarrow #tau#tau");     
  colorv.push_back(kMagenta+2);
  
  fnamev.push_back("/desktop/08a/ksung/ZAna/s09-incmu15-7-mc3_ntuple.root"); xsecv.push_back(0.2969*0.00037*1e+09);
//  fnamev.push_back("/desktop/08a/ksung/ZAna/s09-incmu15-mc3_ntuple.root"); xsecv.push_back(0.5091*0.0002881*1000000000.);
  labelv.push_back("QCD");    
  colorv.push_back(kRed);  
  
  fnamev.push_back("/desktop/08a/ksung/ZAna/s09-ww-7-mc3_ntuple.root"); xsecv.push_back(42.9);
//  fnamev.push_back("/desktop/08a/ksung/ZAna/s09-ww-mc3_ntuple.root"); xsecv.push_back(71.4);
  labelv.push_back("WW");     
  colorv.push_back(kYellow+2);
    
  //
  // Cuts and requirements
  //
  const Double_t muPtCut   = 20;
  const Double_t muEtaCut  = 3;
  const Double_t massMin   = 40;
  const Double_t trkIsoCut = 3;
  const UInt_t trigger     = kHLT_Mu15;
  
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
  vector <TH1F*> hZMassv;
  vector <TH1F*> hZPtv;
  vector <TH1F*> hZyv;
  vector <TH1F*> hZPhiv;
  
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
    char hname[100];    
    sprintf(hname,"hZMass_%i",ifile); hZMassv.push_back(new TH1F(hname,"",180,20,200)); hZMassv[ifile]->Sumw2();
    sprintf(hname,"hZPt_%i",ifile);   hZPtv.push_back(new TH1F(hname,"",200,0,200));    hZPtv[ifile]->Sumw2();
    sprintf(hname,"hZy_%i",ifile);    hZyv.push_back(new TH1F(hname,"",50,-5,5));       hZyv[ifile]->Sumw2();
    sprintf(hname,"hZPhi_%i",ifile);  hZPhiv.push_back(new TH1F(hname,"",50,-5,5));     hZPhiv[ifile]->Sumw2(); 
  }
 
  //
  // Set up TTree for unbinned likelihood fit
  //
  Double_t brMass, brWeight;
  TTree massTree("massTree","Mass Tree");
  massTree.SetDirectory(0);  // force tree to be memory-resident
  massTree.Branch("mass",&brMass,"mass/D");
  massTree.Branch("weight",&brWeight,"weight/D");
      
  //
  // Access samples and fill histograms
  //
  TFile *infile=0;
  TTree *infoTree=0, *genTree=0, *dimuonTree=0, *muonTree=0;
    
  // Data structures to store info from TTrees
  TEventInfo eventInfo;
  TGenInfo   genInfo;
  TDimuon    dimuon;
  TMuon      muon;
  
  // flag for if there is a data sample to process
  Bool_t hasData = (fnamev[0].Length() > 0);
  
  // loop over samples
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {  
    
    // if there is no data sample, skip to MC samples
    if((ifile==0) && !hasData) {
      labelv[0] = "MC data";
      continue;
    }
    
    // Read input file
    if(fnamev[ifile].Contains("/castor/cern.ch",TString::kIgnoreCase))
      infile = new TRFIOFile(fnamev[ifile]);
    else
      infile = new TFile(fnamev[ifile]); 
    assert(infile);
    
    // Get the TTrees
    infoTree   = (TTree*)infile->FindObjectAny("EventInfo"); assert(infoTree);
    genTree    = (TTree*)infile->FindObjectAny("GenInfo");   assert(genTree);
    dimuonTree = (TTree*)infile->FindObjectAny("Dimuon");    assert(dimuonTree);
    muonTree   = (TTree*)infile->FindObjectAny("Muon");      assert(muonTree);

    // Set branch address to structures that will store the info  
    infoTree  ->SetBranchAddress("EventInfo",&eventInfo);
    genTree   ->SetBranchAddress("GenInfo",  &genInfo);
    dimuonTree->SetBranchAddress("Dimuon",   &dimuon);
    muonTree  ->SetBranchAddress("Muon",     &muon);
    
    // Weight MC samples to the requested luminosity
    if(ifile>0)
      weightv.push_back(lumi*xsecv[ifile]/(Double_t)infoTree->GetEntries());      
   
    // loop over events
    for(UInt_t ientry=0; ientry<dimuonTree->GetEntries(); ientry++) {
      dimuonTree->GetEntry(ientry);
      
      // get corresponding event info
      assert(dimuon.eventId <= infoTree->GetEntries());
      infoTree->GetEntry(dimuon.eventId-1);
      assert(eventInfo.eventId == dimuon.eventId);
            
      if(!(eventInfo.triggerBits & trigger))                                 continue;  // no trigger accept? Skip to next event...                                     
      if(dimuon.mass < massMin)                                              continue;  // outside mass region? Skip to next event...
      if((dimuon.pt_1 < muPtCut) || (dimuon.pt_2 < muPtCut))                 continue;  // below muon pT cut? Skip to next event...
      if((fabs(dimuon.eta_1) > muEtaCut) && (fabs(dimuon.eta_2) > muEtaCut)) continue;  // outside eta range? Skip to next event...
      if(dimuon.q_1 == dimuon.q_2)                                           continue;  // same charge? Skip to next event...    
      
      // Check if the Z candidate satisfies a category and if that category is 
      // the one we're interested in. Check categories 1 ~ 5 in sequence:
      // (check category 1, if fail check category 2, if fail...etc.)
      if((category==0) &&
         ((dimuon.typeBits_1 & kGlobal)     && (dimuon.typeBits_2 & kGlobal))   &&           // 2 global muons
         ((dimuon.trkIso03_1 < trkIsoCut)   && (dimuon.trkIso03_2 < trkIsoCut)) &&           // isolation
         ((dimuon.hltMatchBits_1 & trigger) || (dimuon.hltMatchBits_2 & trigger))) {         // at least one match to HLT
	 // Event passes standard selection
	  	 	 
      } else if(((dimuon.typeBits_1 & kGlobal)     && (dimuon.typeBits_2 & kGlobal))   &&    // 2 global muons
                ((dimuon.trkIso03_1 < trkIsoCut)   && (dimuon.trkIso03_2 < trkIsoCut)) &&    // isolation
	        ((dimuon.hltMatchBits_1 & trigger) && (dimuon.hltMatchBits_2 & trigger))) {  // both match to HLT
	// mu-mu-2HLT event	
	if(category!=1) continue;  // Not the category we want. Skip to next event... 
      
      } else if(((dimuon.typeBits_1 & kGlobal)     && (dimuon.typeBits_2 & kGlobal))   &&    // 2 global muons                
		((dimuon.trkIso03_1 < trkIsoCut)   && (dimuon.trkIso03_2 < trkIsoCut)) &&    // isolation
		((dimuon.hltMatchBits_1 & trigger) != (dimuon.hltMatchBits_2 & trigger))) {  // only one match to HLT
	// mu-mu-1HLT event
	if(category!=2) continue;  // Not the category we want. Skip to next event... 
      
      } else if(( ((dimuon.typeBits_1 & kGlobal) && ((dimuon.typeBits_2 & (kGlobal+kTracker)) == kTracker)) ||    // 1 global + 1 tracker OR...
                  ((dimuon.typeBits_2 & kGlobal) && ((dimuon.typeBits_1 & (kGlobal+kTracker)) == kTracker)) ) &&  // 1 tracker + 1 global
		  ((dimuon.trkIso03_1 < trkIsoCut) && (dimuon.trkIso03_2 < trkIsoCut))) {                         // isolation
	// mu-tk event
	if(category!=3) continue;  // Not the category we want. Skip to next event... 
	 	
      } else if(( ((dimuon.typeBits_1 & kGlobal) && (dimuon.typeBits_2 == kStandalone)) ||    // 1 global + 1 standalone OR...
                  ((dimuon.typeBits_2 & kGlobal) && (dimuon.typeBits_1 == kStandalone)) ) &&  // 1 standalone + 1 global
                ((dimuon.trkIso03_1 < trkIsoCut) && (dimuon.trkIso03_2 < trkIsoCut))) {       // isolation
	// mu-sa event
	if(category!=4) continue;  // Not the category we want. Skip to next event...
	
      } else if(((dimuon.typeBits_1 & kGlobal)     && (dimuon.typeBits_2 & kGlobal))   &&    // 2 global muons
                ((dimuon.trkIso03_1 < trkIsoCut)   != (dimuon.trkIso03_2 < trkIsoCut)) &&    // at least one fails isolation
                ((dimuon.hltMatchBits_1 & trigger) || (dimuon.hltMatchBits_2 & trigger))) {  // at least one match to HLT
        // mu-mu-noIso event
	if(category!=5) continue;  // Not the category we want. Skip to next event...
	
      } else {
	// Event does not pass selection cut of any category. Skip to next event...
	continue;
      }
      
      // Selection cuts passed. Fill histograms.
      hZMassv[ifile]->Fill(dimuon.mass,weightv[ifile]);
      hZPtv[ifile]->Fill(dimuon.pt,weightv[ifile]);
      hZyv[ifile]->Fill(dimuon.y,weightv[ifile]);
      hZPhiv[ifile]->Fill(dimuon.phi,weightv[ifile]);
      
      // If MC only, then the 'data' histograms will accumulate contributions from each sample
      if(!hasData) {
        hZMassv[0]->Fill(dimuon.mass,weightv[ifile]);
        hZPtv[0]->Fill(dimuon.pt,weightv[ifile]);	 
        hZyv[0]->Fill(dimuon.y,weightv[ifile]);	 
        hZPhiv[0]->Fill(dimuon.phi,weightv[ifile]);  
      }
      
      // Fill TTree for RooFit
      if(!hasData || (ifile==0)) {
        brMass = dimuon.mass;
        brWeight = weightv[ifile];
        massTree.Fill();    
      }   
    }

    delete infile;
    infile=0, infoTree=0, genTree=0, dimuonTree=0, muonTree=0; 
  }
  
  //
  // Make plots
  //
  TCanvas *c = MakeCanvas("c","c",canw,canh);
  
  // string buffers
  char pname[50];    // plot name
  char ylabel[50];   // y-axis label
  char text[100];    // text box string
  
  // plot title
  TString title = "";
  if(category==1)      title = "#mu#mu-2HLT";
  else if(category==2) title = "#mu#mu-1HLT";
  else if(category==3) title = "#mut";
  else if(category==4) title = "#mus";
  else if(category==5) title = "#mu#mu-noIso";  

  // set how many decimal places to print depending on "lumi" variable
  Int_t precision = 2;
  if(lumi >= 1) precision = 1;
  if(lumi >= 10) precision = 0;
  
  // dimuon mass
  sprintf(pname,"zmass%i",category);
  sprintf(ylabel,"events / %.1f GeV/c^{2}",hZMassv[0]->GetBinWidth(1));
  CPlot plotZMass(pname,title,"m(#mu#mu) [GeV/c^{2}]",ylabel);
  plotZMass.AddHist1D(hZMassv[0],labelv[0],"E",colorv[0]);
  for(UInt_t ibg=2; ibg<hZMassv.size(); ibg++)
    plotZMass.AddToStack(hZMassv[ibg],labelv[ibg],colorv[ibg]);
  plotZMass.AddToStack(hZMassv[1],labelv[1],colorv[1]);  
  sprintf(text,"#int#font[12]{L}dt = %.*f pb^{-1}",precision,lumi);
  plotZMass.AddTextBox(text,0.21,0.85,0.41,0.8,0);
  plotZMass.SetLogy();
  plotZMass.Draw(c,doSave,format);    

  // dimuon pT
  sprintf(pname,"zpt%i",category);
  sprintf(ylabel,"events / %.1f GeV/c",hZPtv[0]->GetBinWidth(1));
  CPlot plotZPt(pname,title,"p_{T}(#mu#mu) [GeV/c]",ylabel);
  plotZPt.AddHist1D(hZPtv[0],labelv[0],"E",colorv[0]);
  for(UInt_t ibg=2; ibg<hZPtv.size(); ibg++)
    plotZPt.AddToStack(hZPtv[ibg],labelv[ibg],colorv[ibg]);
  plotZPt.AddToStack(hZPtv[1],labelv[1],colorv[1]);  
  sprintf(text,"#int#font[12]{L}dt = %.*f pb^{-1}",precision,lumi);
  plotZPt.AddTextBox(text,0.28,0.85,0.48,0.8,0);
  plotZPt.SetLogy();
  plotZPt.Draw(c,doSave,format); 
  
  // dimuon rapidity
  sprintf(pname,"zy%i",category);
  sprintf(ylabel,"events / %.1f",hZyv[0]->GetBinWidth(1));
  CPlot plotZy(pname,title,"y(#mu#mu)",ylabel);
  plotZy.AddHist1D(hZyv[0],labelv[0],"E",colorv[0]);
  for(UInt_t ibg=2; ibg<hZyv.size(); ibg++)
    plotZy.AddToStack(hZyv[ibg],labelv[ibg],colorv[ibg]);
  plotZy.AddToStack(hZyv[1],labelv[1],colorv[1]);  
  sprintf(text,"#int#font[12]{L}dt = %.*f pb^{-1}",precision,lumi);
  plotZy.AddTextBox(text,0.21,0.85,0.41,0.8,0);
  plotZy.TransLegend(0.1,0);
  //plotZy.SetLogy();
  plotZy.Draw(c,doSave,format); 
  
  // dimuon phi
  sprintf(pname,"zphi%i",category);
  sprintf(ylabel,"events / %.1f",hZPhiv[0]->GetBinWidth(1));
  CPlot plotZPhi(pname,title,"#phi(#mu#mu)",ylabel);
  plotZPhi.AddHist1D(hZPhiv[0],labelv[0],"E",colorv[0]);
  for(UInt_t ibg=2; ibg<hZPhiv.size(); ibg++)
    plotZPhi.AddToStack(hZPhiv[ibg],labelv[ibg],colorv[ibg]);
  plotZPhi.AddToStack(hZPhiv[1],labelv[1],colorv[1]);  
  sprintf(text,"#int#font[12]{L}dt = %.*f pb^{-1}",precision,lumi);
  plotZPhi.AddTextBox(text,0.21,0.85,0.41,0.8,0);
  plotZPhi.TransLegend(0.1,0);
  //plotZPhi.SetLogy();
  plotZPhi.Draw(c,doSave,format); 
    
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
  cout << "  Luminosity = " << setprecision(2) << lumi << " / pb" << endl;
  cout << endl;
  
  cout << "       Data file: " << fnamev[0] << endl;
  
  if(category==1)      { cout << " *** Category: mu-mu-2HLT" << endl; }
  else if(category==2) { cout << " *** Category: mu-mu-1HLT" << endl; }
  else if(category==3) { cout << " *** Category: mu-tk" << endl; }
  else if(category==4) { cout << " *** Category: mu-sa" << endl; }
  else if(category==5) { cout << " *** Category: mu-mu-noIso" << endl; }
  
  cout << setw(20) << "Total events =";
  cout << setw(8) << setprecision(1) << fixed << hZMassv[0]->Integral() << endl;
  
  cout << setw(20) << "MC Signal =";
  cout << setw(8) << setprecision(1) << fixed << hZMassv[1]->Integral();
  cout << " \u00B1 "; 
  cout << setw(4) << setprecision(1) << sqrt(hZMassv[1]->GetEffectiveEntries())*weightv[1];
  cout << ",  xsec = " << setw(9) << setprecision(1) << fixed << xsecv[1] << " pb";
  cout << ",  file: " << fnamev[1] << endl; 
  
  for(UInt_t ibg=2; ibg<hZMassv.size(); ibg++) {    
    cout << setw(17) << "MC Bkgd " << ibg-1 << " =";
    cout << setw(8) << setprecision(1) << fixed << hZMassv[ibg]->Integral();
    cout << " \u00B1 ";
    cout << setw(4) << setprecision(1) << sqrt(hZMassv[ibg]->GetEffectiveEntries())*weightv[ibg];
    cout << ",  xsec = " << setw(9) << setprecision(1) << fixed << xsecv[ibg] << " pb"; 
    cout << ",  file: " << fnamev[ibg] << endl;     
  }
  cout << endl;  
  
  if(doSave) {
    cout << " <> Output save in " << CPlot::sOutDir << "/" << endl;
    cout << endl;
  }
  
  gBenchmark->Show("plotZ");
}
