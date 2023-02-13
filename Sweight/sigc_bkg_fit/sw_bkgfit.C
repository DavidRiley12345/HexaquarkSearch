#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDstD0BG.h"
#include "RooAddPdf.h"
#include "RooWorkspace.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLegend.h"
#include "RooFitResult.h"

void AddModel(RooWorkspace* ws);
void AddData(RooWorkspace* ws);
void plot_weighted();

//main function
void sigc_bkg_fit(){
   RooWorkspace ws("ws");
   AddModel(&ws);
   AddData(&ws);
   
   RooFitResult* fitResult = ws.pdf("model")->fitTo(ws.data("data"), RooFit::Save(true));
   ws.import(*fitResult);

   plot_weighted();
} 

//create a background model for the bkg
void AddModel(RooWorkspace* ws){
   RooRealVar Sigmac_M("Sigmac_M","",2390.0,2700.0);
   RooRealVar dm0("dm0","",2420);
   RooRealVar a("a","",100000,0,10000000);
   RooRealVar b("b","",-0.005,-4,0.);
   RooRealVar c("c","",-0.005,-4,0.);
   RooDstD0BG bkgPDF("bkgPDF","",Sigmac_M,dm0,a,b,c);
   
   RooRealVar nbkg("nbkg","",1,100000);

   RooAddPdf model("model","",RooArgList(bkgPDF),RooArgList(nbkg));

   ws->import(model);
}

void AddData(RooWorkspace* ws){
   TFile* fin = TFile::Open("/home/ppe/d/driley/git_HexaquarkSummer/Sweight/allvars/weighted.root"); 
   TTree* tree = (TTree*)fin->Get("weightedtree");  
   RooRealVar* Sigmac_M = ws->var("Sigmac_M");
   RooDataSet *datain = new RooDataSet("data","",tree,RooArgSet(*Sigmac_M));
   ws->import(*datain);
};

void plot_weighted()
{ 
   TFile* fin = TFile::Open("/home/ppe/d/driley/git_HexaquarkSummer/Sweight/allvars/weighted.root"); 
   TTree* tree = (TTree*)fin->Get("weightedtree");

   // Variables and their ranges
   const int nVars = 1;
   string vars[nVars] = {"Sigmac_M"};
   double minVals[nVars] = {2350};
   double maxVals[nVars] = {3000};

   // Loop over variables
   for (int i = 0; i < nVars; i++) {
      TH1D *h_sig = new TH1D("h_sig", ("Sweighted " + vars[i]).c_str(), 80, minVals[i], maxVals[i]);
      TH1D *h_bkg = new TH1D("h_bkg", ("Bweighted " + vars[i]).c_str(), 80, minVals[i], maxVals[i]);

      TCanvas* cdata = new TCanvas("","",800,800);
      cdata->cd(1);
      tree->Draw((vars[i] + ">>h_sig").c_str(), "nsig_sw","norm");
      tree->Draw((vars[i] + ">>h_bkg").c_str(), "nbkg_sw","norm");

      h_sig->SetLineColor(kRed);
      h_bkg->SetLineColor(kBlue);

      h_sig->Draw();
      h_bkg->Draw("same");
    
      TLegend *legend = new TLegend(0.7,0.7,0.9,0.9);
      legend->AddEntry(h_sig,"Signal");
      legend->AddEntry(h_bkg,"Background");
      legend->Draw();

      cdata->SaveAs((vars[i] + ".gif").c_str());
  }
}