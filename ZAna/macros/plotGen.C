//================================================================================================
//
//  Plot generator level distributions for signal sample
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <TFile.h>
#include <TRFIOFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TBenchmark.h>
#include <iostream>
#include <iomanip>
#endif

#include "CPlot.hh"
#include "MitStyleRemix.hh"

#include "ZAnaStructDefs.hh"     // define structures to read in ntuple

void plotGen() 
{
  gBenchmark->Start("plotGen");
  gSystem->Load("libRFIO");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================
  
  Bool_t doSave  = true;   // save plots?
  CPlot::sOutDir = ".";    // output directory
  TString format = "png";  // output file format
  
  // signal sample
  //---------------
  TString fnameSig("/castor/cern.ch/user/k/ksung/ntuples/00/s09-zmm-7-mc3_ntuple.root");
//  TString fnameSig("/castor/cern.ch/user/k/ksung/ntuples/00/s09-zmm-mc3_ntuple.root");
  TString labelSig("Z #rightarrow #mu#mu"); 

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
  TH1F *hZMass1 = new TH1F("ZMass1","",200,0,1000);
  TH1F *hZMass2 = new TH1F("ZMass2","",200,0,200);
  TH1F *hZPt    = new TH1F("ZPt","",200,0,500);
  TH1F *hZy     = new TH1F("Zy","",200,-8,8);
  TH1F *hZPhi   = new TH1F("ZPhi","",160,-3.2,3.2);
  
  TH1F *hDimuonMass1 = new TH1F("DimuonMass1","",200,0,1000);
  TH1F *hDimuonMass2 = new TH1F("DimuonMass2","",200,0,200);
  TH1F *hDimuonPt    = new TH1F("DimuonPt","",200,0,500);
  TH1F *hDimuony     = new TH1F("Dimuony","",200,-8,8);
  TH1F *hDimuonPhi   = new TH1F("DimuonPhi","",160,-3.2,3.2);
  
  TH1F *hMuonPt  = new TH1F("MuonPt","",250,0,500);
  TH1F *hMuonEta = new TH1F("MuonEta","",200,-8,8);
  TH1F *hMuonPhi = new TH1F("MuonPhi","",160,-3.2,3.2);
  
  // Data structures to store info from TTrees
  TEventInfo eventInfo;
  TGenInfo   genInfo;
  
  // signal sample
  TFile *infile;
  if(fnameSig.Contains("/castor/cern.ch",TString::kIgnoreCase))
    infile = new TRFIOFile(fnameSig);
  else
    infile = new TFile(fnameSig); 
  assert(infile);
    
  TTree *infoTree = (TTree*)infile->FindObjectAny("EventInfo"); assert(infoTree);
  TTree *genTree  = (TTree*)infile->FindObjectAny("GenInfo");   assert(genTree);

  // Set branch address to structures that will store the info  
  infoTree  ->SetBranchAddress("EventInfo",&eventInfo);
  genTree   ->SetBranchAddress("GenInfo",  &genInfo);  
  
  // loop over events
  for(UInt_t ientry=0; ientry<genTree->GetEntries(); ientry++) {
    genTree->GetEntry(ientry);
    
    hZMass1->Fill(genInfo.zmass);
    hZMass2->Fill(genInfo.zmass);
    hZPt   ->Fill(genInfo.zpt);
    hZy    ->Fill(genInfo.zy);
    hZPhi  ->Fill(genInfo.zphi);
    
    hDimuonMass1->Fill(genInfo.mass);
    hDimuonMass2->Fill(genInfo.mass);
    hDimuonPt   ->Fill(genInfo.pt);
    hDimuony    ->Fill(genInfo.y);
    hDimuonPhi  ->Fill(genInfo.phi);
    
    hMuonPt->Fill(genInfo.pt_1);   hMuonPt->Fill(genInfo.pt_2);
    hMuonEta->Fill(genInfo.eta_1); hMuonEta->Fill(genInfo.eta_2);
    hMuonPhi->Fill(genInfo.phi_1); hMuonPhi->Fill(genInfo.phi_2);
  }
  
  //
  // Make plots
  //
  TCanvas *c = MakeCanvas("c","c",canw,canh);
  char binstr[20];
  char text[100];  
  
  sprintf(binstr,"events / %.1f GeV/c^{2}",hZMass1->GetBinWidth(1));
  CPlot plotGenMass1("GenMass1","","mass [GeV/c^{2}]",binstr);
  plotGenMass1.AddHist1D(hZMass1,"");
  plotGenMass1.AddHist1D(hDimuonMass1,"",kRed,5);
  sprintf(text,"%i Z events",(Int_t)hZMass1->Integral()); 
  plotGenMass1.AddTextBox(text,0.65,0.75,0.9,0.85,0);
  sprintf(text,"%i #mu^{+}#mu^{-} events",(Int_t)hDimuonMass1->Integral());
  plotGenMass1.AddTextBox(text,0.65,0.65,0.92,0.75,0,kRed);
  plotGenMass1.SetLogy();
  plotGenMass1.Draw(c,doSave,format);
  
  sprintf(binstr,"events / %.1f GeV/c^{2}",hZMass2->GetBinWidth(1));
  CPlot plotGenMass2("GenMass2","","mass [GeV/c^{2}]",binstr);
  plotGenMass2.AddHist1D(hZMass2,"");
  plotGenMass2.AddHist1D(hDimuonMass2,"",kRed,5);
  sprintf(text,"%i Z events",(Int_t)hZMass2->Integral());
  plotGenMass2.AddTextBox(text,0.65,0.75,0.9,0.85,0);
  sprintf(text,"%i #mu^{+}#mu^{-} events",(Int_t)hDimuonMass2->Integral());
  plotGenMass2.AddTextBox(text,0.65,0.65,0.92,0.75,0,kRed);
  plotGenMass2.SetLogy();
  plotGenMass2.Draw(c,doSave,format);
  
  sprintf(binstr,"events / %.1f GeV/c",hZPt->GetBinWidth(1));
  CPlot plotGenPt("GenPt","","p_{T} [GeV/c]",binstr);
  plotGenPt.AddHist1D(hZPt,"");
  plotGenPt.AddHist1D(hDimuonPt,"",kRed,5);
  sprintf(text,"%i Z events",(Int_t)hZPt->Integral());
  plotGenPt.AddTextBox(text,0.65,0.75,0.9,0.85,0);
  sprintf(text,"%i #mu^{+}#mu^{-} events",(Int_t)hDimuonPt->Integral());
  plotGenPt.AddTextBox(text,0.65,0.65,0.92,0.75,0,kRed);
  plotGenPt.SetLogy();
  plotGenPt.Draw(c,doSave,format);
  
  sprintf(binstr,"events / %.2f",hZy->GetBinWidth(1));
  CPlot plotGeny("Geny","","y",binstr);
  plotGeny.AddHist1D(hZy,"");
  plotGeny.AddHist1D(hDimuony,"",kRed,5);
  sprintf(text,"%i Z events",(Int_t)hZy->Integral());
  plotGeny.AddTextBox(text,0.45,0.25,0.7,0.35,0);
  sprintf(text,"%i #mu^{+}#mu^{-} events",(Int_t)hDimuony->Integral());
  plotGeny.AddTextBox(text,0.45,0.15,0.72,0.25,0,kRed);
  plotGeny.Draw(c,doSave,format);
  
  sprintf(binstr,"events / %.2f radians",hZPhi->GetBinWidth(1));
  CPlot plotGenPhi("GenPhi","","#phi",binstr);
  plotGenPhi.AddHist1D(hZPhi,"");
  plotGenPhi.AddHist1D(hDimuonPhi,"",kRed,5);
  sprintf(text,"%i Z events",(Int_t)hZPhi->Integral());
  plotGenPhi.AddTextBox(text,0.65,0.25,0.9,0.35,0);
  sprintf(text,"%i #mu^{+}#mu^{-} events",(Int_t)hDimuonPhi->Integral());
  plotGenPhi.AddTextBox(text,0.65,0.15,0.93,0.25,0,kRed);
  plotGenPhi.Draw(c,doSave,format);
  
  sprintf(binstr,"events / %.1f GeV/c",hMuonPt->GetBinWidth(1));
  CPlot plotGenMuonPt("GenMuonPt","","muon p_{T} [GeV/c]",binstr);
  plotGenMuonPt.AddHist1D(hMuonPt);
  sprintf(text,"%i events",(Int_t)hMuonPt->Integral());
  plotGenMuonPt.AddTextBox(text,0.65,0.8,0.9,0.9,0);
  plotGenMuonPt.SetLogy();
  plotGenMuonPt.Draw(c,doSave,format);
  
  sprintf(binstr,"events / %.2f",hMuonEta->GetBinWidth(1));  
  CPlot plotGenMuonEta("GenMuonEta","","muon #eta",binstr);
  plotGenMuonEta.AddHist1D(hMuonEta);
  sprintf(text,"%i events",(Int_t)hMuonPhi->Integral());
  plotGenMuonEta.AddTextBox(text,0.22,0.8,0.47,0.9,0);
  plotGenMuonEta.Draw(c,doSave,format);
  
  sprintf(binstr,"events / %.2f radians",hMuonPhi->GetBinWidth(1));  
  CPlot plotGenMuonPhi("GenMuonPhi","","muon #phi",binstr); 
  plotGenMuonPhi.AddHist1D(hMuonPhi);
  sprintf(text,"%i events",(Int_t)hMuonPhi->Integral());
  plotGenMuonPhi.AddTextBox(text,0.65,0.15,0.9,0.25,0);
  plotGenMuonPhi.Draw(c,doSave,format); 
    
  //
  // Summary print out
  //
  cout << endl;
  cout << "*" << endl; 
  cout << "* Number of generated Z: " << genTree->GetEntries() << endl; 
  cout << "*" << endl; 
  cout << endl;
  
  gBenchmark->Show("plotGen");
}
