//================================================================================================
//
// Plot muon related distributions
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
#endif

#include "CPlot.hh"
#include "MitStyleRemix.hh"

#include "ZAnaStructDefs.hh"     // define structures to read in ntuple

void plotMuons() {
  
  gBenchmark->Start("plotZ");
  gSystem->Load("libRFIO");
  
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
  const Double_t muPtCut  = 20;
  const Double_t muEtaCut = 2;
  
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
  vector<TH1F*> hPtv;
  vector<TH1F*> hEtav;
  vector<TH1F*> hPhiv;
  vector<TH1F*> hTrkIso03v;
  vector<TH1F*> hTrkIso05v;
  vector<TH1F*> hRelTrkIso03v;
  vector<TH1F*> hRelTrkIso05v;
  vector<TH1F*> hD0v;
  vector<TH1F*> hD0Sigv;
  vector<TH1F*> hCaloCompv;
  vector<TH1F*> hSegCompv;
  vector<TH1F*> hNTkHitsv;
  vector<TH1F*> hNChi2v;

  char hname[100];
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
    sprintf(hname,"hPt_%i",ifile);          hPtv.push_back(new TH1F(hname,"",100,0,100));
    sprintf(hname,"hEta_%i",ifile);         hEtav.push_back(new TH1F(hname,"",150,-3,3));
    sprintf(hname,"hPhi_%i",ifile);         hPhiv.push_back(new TH1F(hname,"",40,-3.2,3.2));
    sprintf(hname,"hTrkIso03_%i",ifile);    hTrkIso03v.push_back(new TH1F(hname,"",100,0,10));
    sprintf(hname,"hTrkIso05_%i",ifile);    hTrkIso05v.push_back(new TH1F(hname,"",100,0,10));
    sprintf(hname,"hRelTrkIso03_%i",ifile); hRelTrkIso03v.push_back(new TH1F(hname,"",100,0,0.5));
    sprintf(hname,"hRelTrkIso05_%i",ifile); hRelTrkIso05v.push_back(new TH1F(hname,"",100,0,0.5));
    sprintf(hname,"hD0_%i",ifile);          hD0v.push_back(new TH1F(hname,"",100,-2,2));
    sprintf(hname,"hD0Sig_%i",ifile);       hD0Sigv.push_back(new TH1F(hname,"",100,-5,5));
    sprintf(hname,"hCaloComp_%i",ifile);    hCaloCompv.push_back(new TH1F(hname,"",100,0,1));    
    sprintf(hname,"hSegComp_%i",ifile);     hSegCompv.push_back(new TH1F(hname,"",100,0,1));
    sprintf(hname,"hNTkHits_%i",ifile);     hNTkHitsv.push_back(new TH1F(hname,"",30,-0.5,29.5));
    sprintf(hname,"hNChi2_%i",ifile);       hNChi2v.push_back(new TH1F(hname,"",100,0,5));        
  }
  
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
    
    for(UInt_t ientry=0; ientry<muonTree->GetEntries(); ientry++) {
      muonTree->GetEntry(ientry);
      
      if(fabs(muon.eta) < muEtaCut) hPtv[ifile]->Fill(muon.pt);
      
      if(muon.pt > muPtCut) hEtav[ifile]->Fill(muon.eta);
      
      if((fabs(muon.eta) > muEtaCut) || (muon.pt < muPtCut)) continue;
  
      hPhiv[ifile]        ->Fill(muon.phi);
      hTrkIso03v[ifile]   ->Fill(muon.trkIso03);
      hTrkIso05v[ifile]   ->Fill(muon.trkIso05);
      hRelTrkIso03v[ifile]->Fill(muon.trkIso03/muon.pt);
      hRelTrkIso05v[ifile]->Fill(muon.trkIso05/muon.pt);
      hD0v[ifile]         ->Fill(muon.d0);
      hD0Sigv[ifile]      ->Fill(muon.d0/muon.d0Err);
      hCaloCompv[ifile]   ->Fill(muon.caloComp);
      hSegCompv[ifile]    ->Fill(muon.segComp);
      hNTkHitsv[ifile]    ->Fill(muon.nTkHits);
      hNChi2v[ifile]      ->Fill(muon.nchi2);           
    }

    delete infile;
    infile=0, infoTree=0, genTree=0, dimuonTree=0, muonTree=0; 
  }
  
  //
  // Make plots
  //
  TCanvas *c = MakeCanvas("c","c",canw,canh);
  char binstr[20];  
  char text[100];
  
  // muon pT
  sprintf(binstr,"fraction / %.1f GeV/c",hPtv[0]->GetBinWidth(1));
  CPlot plotPt("pt","","muon p_{T} [GeV/c]","fraction / 1.0 GeV/c");
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hPtv[i]->Scale(1.0/hPtv[i]->Integral());
    plotPt.AddHist1D(hPtv[i],labelv[i],"hist",colorv[i]); 
  }
  sprintf(text,"|#eta| < %.1f",muEtaCut);
  plotPt.AddTextBox(text,0.63,0.5,0.76,0.55,0);
  plotPt.Draw(c,doSave,format);
  
  // muon eta
  sprintf(binstr,"fraction / %.2f",hEtav[0]->GetBinWidth(1)); 
  CPlot plotEta("eta","","muon #eta",binstr);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hEtav[i]->Scale(1.0/hEtav[i]->Integral());
    plotEta.AddHist1D(hEtav[i],labelv[i],"hist",colorv[i]); 
  }
  sprintf(text,"p_{T} > %i GeV/c",(Int_t)muPtCut);
  plotEta.AddTextBox(text,0.23,0.84,0.45,0.9,0);
  plotEta.TransLegend(0.12,0);
  plotEta.Draw(c,doSave,format);
  
  // muon phi
  sprintf(binstr,"fraction / %.1f radians",hPhiv[0]->GetBinWidth(1));
  CPlot plotPhi("phi","","muon #phi",binstr);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hPhiv[i]->Scale(1.0/hPhiv[i]->Integral());
    plotPhi.AddHist1D(hPhiv[i],labelv[i],"hist",colorv[i]); 
  }
  sprintf(text,"p_{T} > %i GeV/c, |#eta| < %.1f",(Int_t)muPtCut,muEtaCut);  
  plotPhi.AddTextBox(text,0.23,0.84,0.55,0.9,0);
  plotPhi.TransLegend(-0.2,-0.4);  
  plotPhi.Draw(c,doSave,format);
  
  // track isolation 0.3
  sprintf(binstr,"fraction / %.1f GeV/c",hTrkIso03v[0]->GetBinWidth(1));
  CPlot plotTrkIso03("trkIso03","","track isolation (#DeltaR<0.3) [GeV/c]",binstr);
  for(UInt_t i=0; i<fnamev.size(); i++) {
    hTrkIso03v[i]->Scale(1.0/hTrkIso03v[i]->Integral()); 
    plotTrkIso03.AddHist1D(hTrkIso03v[i],labelv[i],"hist",colorv[i]); 
  }
  sprintf(text,"p_{T} > %i GeV/c, |#eta| < %.1f",(Int_t)muPtCut,muEtaCut);  
  plotTrkIso03.AddTextBox(text,0.23,0.84,0.55,0.9,0);
  plotTrkIso03.SetLogy();
  plotTrkIso03.Draw(c,doSave,format);
  
  // track isolation 0.5
  sprintf(binstr,"fraction / %.2f GeV/c",hTrkIso05v[0]->GetBinWidth(1));
  CPlot plotTrkIso05("trkIso05","","track isolation (#DeltaR<0.5) [GeV/c]",binstr);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hTrkIso05v[i]->Scale(1.0/hTrkIso05v[i]->Integral());
    plotTrkIso05.AddHist1D(hTrkIso05v[i],labelv[i],"hist",colorv[i]); 
  }
  sprintf(text,"p_{T} > %i GeV/c, |#eta| < %.1f",(Int_t)muPtCut,muEtaCut);  
  plotTrkIso05.AddTextBox(text,0.23,0.84,0.55,0.9,0);
  plotTrkIso05.SetLogy();
  plotTrkIso05.Draw(c,doSave,format);
  
  // relative track isolation 0.3
  sprintf(binstr,"fraction / %.3f",hRelTrkIso03v[0]->GetBinWidth(1));
  CPlot plotRelTrkIso03("relTrkIso03","","relative track isolation (#DeltaR<0.3)",binstr);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hRelTrkIso03v[i]->Scale(1.0/hRelTrkIso03v[i]->Integral());
    plotRelTrkIso03.AddHist1D(hRelTrkIso03v[i],labelv[i],"hist",colorv[i]); 
  }
  sprintf(text,"p_{T} > %i GeV/c, |#eta| < %.1f",(Int_t)muPtCut,muEtaCut);
  plotRelTrkIso03.AddTextBox(text,0.23,0.84,0.55,0.9,0);
  plotRelTrkIso03.SetLogy();
  plotRelTrkIso03.Draw(c,doSave,format);
  
  // relative track isolation 0.5
  sprintf(binstr,"fraction / %.3f",hRelTrkIso05v[0]->GetBinWidth(1));
  CPlot plotRelTrkIso05("relTrkIso05","","relative track isolation (#DeltaR<0.5)",binstr);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hRelTrkIso05v[i]->Scale(1.0/hRelTrkIso05v[i]->Integral());
    plotRelTrkIso05.AddHist1D(hRelTrkIso05v[i],labelv[i],"hist",colorv[i]); 
  }
  sprintf(text,"p_{T} > %i GeV/c, |#eta| < %.1f",(Int_t)muPtCut,muEtaCut);
  plotRelTrkIso05.AddTextBox(text,0.23,0.84,0.55,0.9,0);
  plotRelTrkIso05.SetLogy();
  plotRelTrkIso05.Draw(c,doSave,format);
  
  // impact parameter
  sprintf(binstr,"fraction / %.2f cm",hD0v[0]->GetBinWidth(1));
  CPlot plotD0("d0","","d_{0} [cm]",binstr);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hD0v[i]->Scale(1.0/hD0v[i]->Integral());
    plotD0.AddHist1D(hD0v[i],labelv[i],"hist",colorv[i]); 
  }
  sprintf(text,"p_{T} > %i GeV/c, |#eta| < %.1f",(Int_t)muPtCut,muEtaCut);
  plotD0.AddTextBox(text,0.23,0.84,0.55,0.9,0);
  plotD0.SetLogy();
  plotD0.Draw(c,doSave,format);
  
  // impact parameter significance
  sprintf(binstr,"fraction / %.1f",hD0Sigv[0]->GetBinWidth(1));
  CPlot plotD0Sig("d0Sig","","d_{0}/#sigma",binstr);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hD0Sigv[i]->Scale(1.0/hD0Sigv[i]->Integral());
    plotD0Sig.AddHist1D(hD0Sigv[i],labelv[i],"hist",colorv[i]); 
  }
  sprintf(text,"p_{T} > %i GeV/c, |#eta| < %.1f",(Int_t)muPtCut,muEtaCut);
  plotD0Sig.AddTextBox(text,0.23,0.84,0.55,0.9,0);
  plotD0Sig.Draw(c,doSave,format);
  
  // calorimeter compatability
  sprintf(binstr,"fraction / %.2f",hCaloCompv[0]->GetBinWidth(1));
  CPlot plotCaloComp("caloComp","","calorimeter compatibility",binstr);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hCaloCompv[i]->Scale(1.0/hCaloCompv[i]->Integral());
    plotCaloComp.AddHist1D(hCaloCompv[i],labelv[i],"hist",colorv[i]); 
  }
  sprintf(text,"p_{T} > %i GeV/c, |#eta| < %.1f",(Int_t)muPtCut,muEtaCut);
  plotCaloComp.AddTextBox(text,0.23,0.84,0.55,0.9,0);
  plotCaloComp.Draw(c,doSave,format);
  
  // segment compatibility
  sprintf(binstr,"fraction / %.2f",hSegCompv[0]->GetBinWidth(1));
  CPlot plotSegComp("segComp","","segment compatibility",binstr);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hSegCompv[i]->Scale(1.0/hSegCompv[i]->Integral());
    plotSegComp.AddHist1D(hSegCompv[i],labelv[i],"hist",colorv[i]); 
  }
  sprintf(text,"p_{T} > %i GeV/c, |#eta| < %.1f",(Int_t)muPtCut,muEtaCut);
  plotSegComp.AddTextBox(text,0.23,0.84,0.55,0.9,0);
  plotSegComp.Draw(c,doSave,format);
  
  // number of tracker hits
  CPlot plotNTkHits("nTkHits","","number of tracker hits","fraction");
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hNTkHitsv[i]->Scale(1.0/hNTkHitsv[i]->Integral());
    plotNTkHits.AddHist1D(hNTkHitsv[i],labelv[i],"hist",colorv[i]); 
  }
  sprintf(text,"p_{T} > %i GeV/c, |#eta| < %.1f",(Int_t)muPtCut,muEtaCut);
  plotNTkHits.AddTextBox(text,0.23,0.84,0.55,0.9,0);
  plotNTkHits.TransLegend(0.1,0);
  plotNTkHits.Draw(c,doSave,format);
  
  // normalized chi^2 for muon track
    sprintf(binstr,"fraction / %.2f",hNChi2v[0]->GetBinWidth(1));
  CPlot plotNChi2("nchi2","","#chi^{2}/NDF",binstr);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hNChi2v[i]->Scale(1.0/hNChi2v[i]->Integral());
    plotNChi2.AddHist1D(hNChi2v[i],labelv[i],"hist",colorv[i]); 
  }
  sprintf(text,"p_{T} > %i GeV/c, |#eta| < %.1f",(Int_t)muPtCut,muEtaCut);
  plotNChi2.AddTextBox(text,0.55,0.44,0.87,0.5,0);
  plotNChi2.Draw(c,doSave,format);
  
  gBenchmark->Show("plotZ");
}
