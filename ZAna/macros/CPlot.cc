#include "CPlot.hh"
#include <TLatex.h>

CPlot::CPlot()
{
  TString name = "plot"; 
  name += sCount;
  CPlot(name,"","",""); 
}

CPlot::CPlot(TString name, TString title, TString xtitle, TString ytitle):
fStack(0),
fName(name),
fTitle(title),
fXTitle(xtitle),
fYTitle(ytitle),
fXmin(0),
fXmax(0),
fYmin(0),
fYmax(0),
fLogx(0),
fLogy(0),
fGridx(0),
fGridy(0),
fRebin(1),
fLeg(0),
fShowStats(0),
fStatsX(0.63),
fStatsY(0.90)
{
  sCount++;
}

//--------------------------------------------------------------------------------------------------
void CPlot::AddHist1D(TH1F *h, TString drawopt, int color, int linesty, int fillsty)
{
  if(!h)
    return;
   
  h->SetLineColor(color);
  h->SetLineStyle(linesty);
  h->SetFillColor(color);
  h->SetFillStyle(fillsty);  
  
  CPlotItem item;
  item.hist1D = h;
  item.drawopt = drawopt;
  fItems.push_back(item);
}

void CPlot::AddHist1D(TH1F *h, TString label, TString drawopt, int color, int linesty, int fillsty)
{
  if(!h)
    return;

  if(!fLeg)
    fLeg = new TLegend(0.6,0.84,0.93,0.9);
  else
    fLeg->SetY1(fLeg->GetY1()-0.06);
  fLeg->AddEntry(h,label,"L"); 
  
  AddHist1D(h,drawopt,color,linesty,fillsty);
}

void CPlot::AddHist1D(TFile *f, TString histName, TString drawopt, int color, int linesty, int fillsty)
{
  if(!f)
    return;
  
  TH1F *h = (TH1F*)f->FindObjectAny(histName);
  AddHist1D(h,drawopt,color,linesty,fillsty);
}

void CPlot::AddHist1D(TFile *f, TString histName, TString label, TString drawopt, int color, int linesty, int fillsty)
{
  if(!f)
    return;
  
  TH1F *h = (TH1F*)f->FindObjectAny(histName);
  AddHist1D(h,label,drawopt,color,linesty,fillsty);
}

//--------------------------------------------------------------------------------------------------
void CPlot::AddToStack(TH1F *h, int color)
{
  if(!h)
    return;
    
  if(!fStack)
    fStack = new THStack();
  
  fStack->Add(h);
  AddHist1D(h,"",color,1,1001);
}

void CPlot::AddToStack(TH1F *h, TString label, int color)
{
  if(!h)
    return;
  
  if(!fStack)
    fStack = new THStack();
  
  fStack->Add(h);
  AddHist1D(h,label,"",color,1,1001);
}

void CPlot::AddToStack(TFile *f, TString histName, int color)
{
  if(!f)
    return;
  
  TH1F *h = (TH1F*)f->FindObjectAny(histName);
  AddToStack(h,color);
}

void CPlot::AddToStack(TFile *f, TString histName, TString label, int color)
{
  if(!f)
    return;
  
  TH1F *h = (TH1F*)f->FindObjectAny(histName);
  AddToStack(h,label,color);
}  
  
//--------------------------------------------------------------------------------------------------
void CPlot::AddHist2D(TH2F *h, TString drawopt, int fillcolor, int linecolor)
{
  if(!h)
    return;
  
  h->SetLineColor(linecolor);
  h->SetFillColor(fillcolor);
  h->SetMarkerStyle(kFullDotMedium); 
  
  CPlotItem item;
  item.hist2D = h;
  item.drawopt = drawopt;
  fItems.push_back(item);
}

void CPlot::AddHist2D(TFile *f, TString histName, TString drawopt, int fillcolor, int linecolor)
{
  if(!f)
    return;
  
  TH2F *h = (TH2F*)f->FindObjectAny(histName);
  AddHist2D(h,drawopt,linecolor,fillcolor);
}

//--------------------------------------------------------------------------------------------------
void CPlot::AddGraph(TGraph *gr, TString drawopt, int color, int marksty)
{
  if(!gr)
    return;
  
  gr->SetMarkerColor(color);
  gr->SetLineColor(color);
  gr->SetLineWidth(2);
  gr->SetMarkerStyle(marksty);
  gr->SetMarkerSize(0.9);
  
  CPlotItem item;
  item.graph = gr;
  item.drawopt = drawopt;  
  fItems.push_back(item);
}

void CPlot::AddGraph(TGraph *gr, TString label, TString drawopt, int color, int marksty)
{
  if(!gr)
    return;

  if(!fLeg)
    fLeg = new TLegend(0.6,0.84,0.93,0.9);
  else
    fLeg->SetY1(fLeg->GetY1()-0.06);
  fLeg->AddEntry(gr,label,"P"); 
  
  AddGraph(gr,drawopt,color,marksty);
}

void CPlot::AddGraph(TFile *f, TString grName, TString drawopt, int color, int marksty)
{
  if(!f)
    return;
  
  TGraph *gr = (TGraph*)f->FindObjectAny(grName);
  AddGraph(gr,drawopt,color,marksty);
}

void CPlot::AddGraph(TFile *f, TString grName, TString label, TString drawopt, int color, int marksty)
{
  if(!f)
    return;
  
  TGraph *gr = (TGraph*)f->FindObjectAny(grName);
  AddGraph(gr,label,drawopt,color,marksty);
}

//--------------------------------------------------------------------------------------------------
void CPlot::AddTextBox(TString text, double x1, double y1, double x2, double y2,
                       int bordersize, int textcolor, int fillcolor)
{
  TPaveText *tb = new TPaveText(x1,y1,x2,y2,"NDC");
  tb->SetTextColor(textcolor);
  if(fillcolor==-1)
    tb->SetFillStyle(0);
  else
    tb->SetFillColor(fillcolor);
  tb->SetBorderSize(bordersize);
  tb->AddText(text);
  fTextBoxes.push_back(tb);
}

void CPlot::AddTextBox(double x1, double y1, double x2, double y2, 
                       int bordersize, int textcolor, int fillcolor, int nlines,...)
{
  TPaveText *tb = new TPaveText(x1,y1,x2,y2,"NDC");
  tb->SetTextColor(textcolor);
  if(fillcolor==-1)
    tb->SetFillStyle(0);
  else
    tb->SetFillColor(fillcolor);
  tb->SetBorderSize(bordersize);
  tb->SetTextAlign(12);
  
  va_list ap;
  va_start(ap,nlines);
  for(int i=0; i<nlines; i++) {
    TString textline(va_arg(ap,char*));
    tb->AddText(textline);
  }
  va_end(ap);
  
  fTextBoxes.push_back(tb);  
}

//--------------------------------------------------------------------------------------------------
void CPlot::AddLine(double x1, double y1, double x2, double y2, int color, int style)
{
  TLine *line = new TLine(x1,y1,x2,y2);
  line->SetLineColor(color);
  line->SetLineStyle(style);
  line->SetLineWidth(2);
  fLines.push_back(line);
}

void CPlot::AddLine(double x1, double y1, double x2, double y2, 
                    int color, int style, TString label)
{
  TLine *line = new TLine(x1,y1,x2,y2);
  line->SetLineColor(color);
  line->SetLineStyle(style);
  line->SetLineWidth(2);
  fLines.push_back(line);

  if(!fLeg)
    fLeg = new TLegend(0.6,0.84,0.93,0.9);
  else
    fLeg->SetY1(fLeg->GetY1()-0.06);
  fLeg->AddEntry(line,label,"L");
}

//--------------------------------------------------------------------------------------------------
void CPlot::AddFcn(TF1* fcn, int color, int linesty)
{
  if(!fcn)
    return;
  
  fcn->SetLineColor(color);
  fcn->SetLineStyle(linesty);
  fFcns.push_back(fcn);
}

void CPlot::AddFcn(TF1* fcn, TString label, int color, int linesty)
{
  if(!fcn)
    return;
    
  if(!fLeg)
    fLeg = new TLegend(0.6,0.84,0.93,0.9);
  else
    fLeg->SetY1(fLeg->GetY1()-0.06);
  fLeg->AddEntry(fcn,label,"L");
  
  AddFcn(fcn,color,linesty);
} 

//--------------------------------------------------------------------------------------------------
void CPlot::Draw(TCanvas *c, bool doSave, TString format)
{
  if(!fItems.size())
    return;    
  
  int nHist1D=0, nHist2D=0, nGraph=0;
  for(uint i=0; i<fItems.size(); i++) {
    if(fItems[i].hist1D != 0) nHist1D++;
    if(fItems[i].hist2D != 0) nHist2D++;
    if(fItems[i].graph  != 0) nGraph++;
  }
  
  //
  // Draw 2D histogram, save if necessary, then exit
  //   Suggested options for:
  //     contour plot -> "CONT4Z"
  //     lego plot    -> "LEGO1 0"
  //     color plot   -> "COLZ"
  //   Default is scatter plot
  //  
  if(nHist2D>0) {
    for(uint i=0; i<fItems.size(); i++) {
      if(fItems[i].hist2D==0) continue;
      
      fItems[i].hist2D->Draw(fItems[i].drawopt);
      fItems[i].hist2D->SetTitle(fTitle);
      fItems[i].hist2D->GetXaxis()->SetTitle(fXTitle);
      fItems[i].hist2D->GetYaxis()->SetTitle(fYTitle);
    
      for(uint j=0; j<fTextBoxes.size(); j++)
        fTextBoxes[j]->Draw();
    
      gSystem->mkdir(sOutDir,true);
      TString outname = sOutDir+TString("/")+fName+TString(".");        
      if(doSave) {
        if(format.CompareTo("all",TString::kIgnoreCase)==0) {
	  c->SaveAs(outname+TString("png"));
	  c->SaveAs(outname+TString("eps"));
	  c->SaveAs(outname+TString("C"));
	} else {
	  c->SaveAs(outname+format);
	}
      }
      
      return;
    }
  }    
  
  // 
  // Draw 1D histograms
  //   Histograms are cloned so that content and properties 
  //   of the original histograms are not changed
  //
  std::vector<TH1F*> vHists;
  std::vector<TString> vHistOpts;
  if(nHist1D>0) {   
    
    double ymax=0;
    uint ifirst=0;
    
    for(uint i=0; i<fItems.size(); i++) {
      if(fItems[i].hist1D==0) continue;
      if(fStack && fStack->GetHists()->Contains(fItems[i].hist1D)) continue;
      
      TString hname = fName;
      hname += "_h_";
      hname += i;
    
      TH1F *h;     
      if(fRebin>1)      
        h = (TH1F*)fItems[i].hist1D->Rebin(fRebin,hname);
      else
        h = (TH1F*)fItems[i].hist1D->Clone(hname);
  
      if(fXmin < fXmax)
        h->GetXaxis()->SetRangeUser(fXmin,fXmax);
      
      if(fYmin < fYmax) {
        h->GetYaxis()->SetRangeUser(fYmin,fYmax);
      } else {
        if(ymax < h->GetMaximum()) {
	  ymax = h->GetMaximum();
	  ifirst = vHists.size();
	}
      }
      
      vHists.push_back(h);
      vHistOpts.push_back(fItems[i].drawopt);
    }
    
    vHists[ifirst]->SetTitle(fTitle);
    vHists[ifirst]->SetXTitle(fXTitle);
    vHists[ifirst]->SetYTitle(fYTitle);
    vHists[ifirst]->SetLineWidth(2);
    vHists[ifirst]->Draw(vHistOpts[ifirst].Data());
    
    //
    // Draw histogram stack
    //
    if(fStack)
      fStack->Draw("hist same");
        
    for(uint i=0; i<vHists.size(); i++) {
      TH1F *h = vHists[i];              
      h->SetLineWidth(2);
      char opt[100];
      sprintf(opt,"same%s",vHistOpts[i].Data());    
      h->Draw(opt);
    }
  }  
  
  //
  // Draw graphs
  //
  std::vector<TGraph*> vGraphs;
  std::vector<TString> vGraphOpts;
  if(nGraph>0) {    
    for(uint i=0; i<fItems.size(); i++) {
      if(fItems[i].graph==0) continue;
    
      TString grName = fName;
      grName += "_gr_";
      grName += i;
      
      TGraph *gr = (TGraph*)fItems[i].graph->Clone(grName);
      
      if(fXmin < fXmax)
        gr->GetXaxis()->SetLimits(fXmin,fXmax);
//        gr->GetXaxis()->SetRangeUser(fXmin,fXmax);
    
      if(fYmin < fYmax)
        gr->GetYaxis()->SetRangeUser(fYmin,fYmax);
	
      vGraphs.push_back(gr);
      vGraphOpts.push_back(fItems[i].drawopt);
    }
    
    if(vHists.size()==0) {
      vGraphs[0]->SetTitle(fTitle);
      vGraphs[0]->GetXaxis()->SetTitle(fXTitle);
      vGraphs[0]->GetYaxis()->SetTitle(fYTitle);
    }
    
    for(uint i=0; i<vGraphs.size(); i++) {
      TGraph *gr = vGraphs[i];
      char opt[100];
      (i==0 && nHist1D==0) ? sprintf(opt,"AP%s",vGraphOpts[i].Data()) : sprintf(opt,"P%s",vGraphOpts[i].Data());
      gr->Draw(opt);
    }
  }
    
  //
  // Draw legend
  //
  if(fLeg) {
    fLeg->SetFillStyle(0);
    fLeg->SetBorderSize(0);
    fLeg->Draw();
  }
  
  //
  // Draw statistics box
  //
  TLatex *stat=0, *sval=0;
  if(fShowStats) {
    char buffer[20];
    stat = new TLatex[3*vHists.size()];
    sval = new TLatex[3*vHists.size()];
    for(uint i=0; i<vHists.size(); i++) {
      int x = fShowStats;
      
      // number of entries
      if(x / 100) {
        stat[3*i].SetNDC(); stat[3*i].SetTextAlign(13); stat[3*i].SetTextSize(0.03);
        stat[3*i].SetText(fStatsX,fStatsY-0.04*(3*i)-0.005*i,"Entries");
	stat[3*i].SetTextColor(vHists[i]->GetLineColor());
        stat[3*i].Draw();   
        sprintf(buffer,"%i",int(vHists[i]->Integral()));
        sval[3*i].SetNDC(); sval[3*i].SetTextAlign(33); sval[3*i].SetTextSize(0.03);
        sval[3*i].SetText(fStatsX+0.3,fStatsY-0.04*(3*i)-0.005*i,buffer);
	sval[3*i].SetTextColor(vHists[i]->GetLineColor());
        sval[3*i].Draw();
      }
      
      // mean
      x = x % 100;
      if(x / 10) {
        stat[3*i+1].SetNDC(); stat[3*i+1].SetTextAlign(13); stat[3*i+1].SetTextSize(0.03);
        stat[3*i+1].SetText(fStatsX,fStatsY-0.04*(3*i+1)-0.005*i,"Mean");
	stat[3*i+1].SetTextColor(vHists[i]->GetLineColor());
        stat[3*i+1].Draw();   
        sprintf(buffer,"%g",vHists[i]->GetMean());
        sval[3*i+1].SetNDC(); sval[3*i+1].SetTextAlign(33); sval[3*i+1].SetTextSize(0.03);
        sval[3*i+1].SetText(fStatsX+0.3,fStatsY-0.04*(3*i+1)-0.005*i,buffer);
	sval[3*i+1].SetTextColor(vHists[i]->GetLineColor());
        sval[3*i+1].Draw();
      }
      
      // RMS
      x = x % 10;
      if(x) {
        stat[3*i+2].SetNDC(); stat[3*i+2].SetTextAlign(13); stat[3*i+2].SetTextSize(0.03);
        stat[3*i+2].SetText(fStatsX,fStatsY-0.04*(3*i+2)-0.005*i,"RMS");
	stat[3*i+2].SetTextColor(vHists[i]->GetLineColor());
        stat[3*i+2].Draw();   
        sprintf(buffer,"%g",vHists[i]->GetRMS());
        sval[3*i+2].SetNDC(); sval[3*i+2].SetTextAlign(33); sval[3*i+2].SetTextSize(0.03);
        sval[3*i+2].SetText(fStatsX+0.3,fStatsY-0.04*(3*i+2)-0.005*i,buffer);
	sval[3*i+2].SetTextColor(vHists[i]->GetLineColor());
        sval[3*i+2].Draw();
      }
    }
  }
  
  //
  // Draw functions
  //
  for(uint i=0; i<fFcns.size(); i++)
    (i==0 && vHists.size()==0 && vGraphs.size()==0) ? fFcns[i]->Draw() : fFcns[i]->Draw("sameC");
  
  //
  // Draw lines
  //
  for(uint i=0; i<fLines.size(); i++)
    fLines[i]->Draw();
  
  //
  // Draw textboxes
  //
  for(uint i=0; i<fTextBoxes.size(); i++)
    fTextBoxes[i]->Draw();    
  
  //
  // Set log scale if necessary
  // 
  c->SetLogx(fLogx);
  c->SetLogy(fLogy);
  
  //
  // Set grid lines if necessary
  //
  c->SetGridx(fGridx);
  c->SetGridy(fGridy);
  
  //
  // Save plot if necessary
  //
  gSystem->mkdir(sOutDir,true);
  TString outname = sOutDir+TString("/")+fName+TString(".");
  if(doSave) {
    if(format.CompareTo("all",TString::kIgnoreCase)==0) {
      c->SaveAs(outname+TString("png"));
      c->SaveAs(outname+TString("eps"));
      c->SaveAs(outname+TString("C"));
    } else {
      c->SaveAs(outname+format);
    }
    
    delete [] stat;
    delete [] sval;
//    for(uint i=0; i<vHists.size(); i++)
//      delete vHists[i];
  }
}
