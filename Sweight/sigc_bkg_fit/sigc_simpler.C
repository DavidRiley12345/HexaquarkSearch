#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDstD0BG.h"
#include "RooAddPdf.h"
#include "TCanvas.h"
#include "RooFit.h"
#include "TTree.h"

using namespace RooFit;

int sigc_simpler() {
   // Create a RooRealVar for the mass variable
   RooRealVar sig_M("sig_M","",2350.0,3000.0);

   // Create background PDF variables
   RooRealVar dm0("dm0","",2350,2350,2500);
   RooRealVar a("a","",100000,0,10000000);
   RooRealVar b("b","",-0.005,-4,4);
   RooRealVar c("c","",-0.005,-4,4);

   // Create the background PDF
   RooDstD0BG bkgPDF("bkgPDF","",sig_M,dm0,a,b,c);

   // Create the normalization variable for the background PDF
   RooRealVar nbkg("nbkg","",1,100000);

   // Create the total PDF
   RooAddPdf model("model","",RooArgList(bkgPDF),RooArgList(nbkg));

   // Open the data file and get the tree
   TFile* fin = TFile::Open("/home/ppe/d/driley/git_HexaquarkSummer/Sweight/allvars/weighted.root"); 
   TTree* tree = (TTree*)fin->Get("weightedtree");

   // Create the RooDataSet from the tree
   RooRealVar nbkg_sw("nbkg_sw","",-10,10);
   RooArgSet argset(sig_M,nbkg_sw);
   RooDataSet data("data", "",argset,Import(*tree),WeightVar("nbkg_sw"));

   // Fit the model to the data
   model.fitTo(data, Extended());

   // Plot the data and the fit
   TCanvas *c1 = new TCanvas("c1","",800,800);
   RooPlot* frame = sig_M.frame();
   data.plotOn(frame);
   model.plotOn(frame);
   frame->Draw();
   c1->SaveAs("sigma_bkg_fit.gif");
}
