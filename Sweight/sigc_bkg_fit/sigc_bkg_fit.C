#include "RooRealVar.h"
#include "RooStats/SPlot.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooAddition.h"
#include "RooProduct.h"
#include "TCanvas.h"
#include "RooAbsPdf.h"
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooConstVar.h"
#include "ROOT/RDataFrame.hxx"
#include "TTree.h"
// use this order for safety on library loading
using namespace RooFit;
using namespace RooStats;

void AddModel(RooWorkspace*);
void AddData(RooWorkspace*);
void plot_weighted(RooWorkspace*);
//main function
void sigc_bkg_fit(){
   RooWorkspace* wspace = new RooWorkspace("myWS");
   AddModel(wspace);
   AddData(wspace);
   plot_weighted(wspace);
   delete wspace;
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
   RooRealVar nbkg_sw("nbkg_sw","",-10,10);
   RooAddPdf model("model","",RooArgList(bkgPDF),RooArgList(nbkg));
   ws->import(nbkg_sw);
   ws->import(model);
}

void AddData(RooWorkspace* ws){
   TFile* fin = TFile::Open("/home/ppe/d/driley/git_HexaquarkSummer/Sweight/allvars/weighted.root"); 
   TTree* tree = (TTree*)fin->Get("weightedtree");  
   RooRealVar* Sigmac_M = ws->var("Sigmac_M");
   RooRealVar* nbkg_sw = ws->var("nbkg_sw");
   RooArgSet argset(*Sigmac_M, *nbkg_sw);

   RooDataSet *datain = new RooDataSet("","",tree,argset,wgtVarName("nbkg_sw");
   ws->import(*datain, Rename("data"));
};


void plot_weighted(RooWorkspace* ws){ 
   TFile* fin = TFile::Open("/home/ppe/d/driley/git_HexaquarkSummer/Sweight/allvars/weighted.root"); 
   TTree* tree = (TTree*)fin->Get("weightedtree");

   RooAbsPdf* model = ws->("model");
   RooAbsPdf* bkgPDF = ws->("bkgPDF");
   RooRealVar* Sigmac_M = ws->var("Sigmac_M");   


   RooDataSet* data = (RooDataSet*) ws->data("data");
   
   model->fitTo(*data, Extended());
 
   TCanvas* cdata = new TCanvas("","",800,800);
   cdata->cd(1);
   
   RooPlot* frame = Sigmac_M->frame();

   data->plotOn(frame);
   model->plotOn(frame);   
   TLegend *legend = new TLegend(0.7,0.7,0.9,0.9);
   legend->Draw();
   cdata->SaveAs("sigma_bkg_fit.gif");
}

