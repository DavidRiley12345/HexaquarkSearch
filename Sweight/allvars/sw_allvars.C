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

void sw_allvar_makesig_M(){
   ROOT::RDataFrame d("weightedtree","/home/ppe/d/driley/git_HexaquarkSummer/Sweight/weighted.root");
   auto d_new = d.Define("sig_M","sqrt(pow(pi_PE + Lambdacp_PE,2) - pow(pi_PX + Lambdacp_PX,2) - pow(pi_PY + Lambdacp_PY,2) - pow(pi_PZ + Lambdacp_PZ,2))");
   d_new.Snapshot("weightedtree","/home/ppe/d/driley/git_HexaquarkSummer/Sweight/allvars/weighted.root");
} 

void sw_allvars()
{ 
  sw_allvar_makesig_M(); 

  TFile* fin = TFile::Open("/home/ppe/d/driley/git_HexaquarkSummer/Sweight/allvars/weighted.root"); 
  TTree* tree = (TTree*)fin->Get("weightedtree");

  // Variables and their ranges
  const int nVars = 19;
  string vars[nVars] = {"sig_M", "Lambdacp_M", "Lambdacp_P", "Lambdacp_PT", "Lambdacp_p_P", "Lambdacp_p_PT","Lambdacp_K_P", "Lambdacp_K_PT", "Lambdacp_pi_P", "Lambdacp_pi_PT", "pi_P","pi_PT","pi_PX", "pi_PY", "pi_PZ", "Lambdacp_PX", "Lambdacp_PY", "Lambdacp_PZ","totCandidates"};
  double minVals[nVars] = {2350,         2220,            0,             0,              0,               0,             0,               0,               0,                0,      0,      0, -10000,  -10000,       0,        -30000,        -30000,     0,0 };
  double maxVals[nVars] = {3000,         2360,       250000,         20000,         120000,           10000,        100000,           6000,          60000,           6000, 30000,   3000,  10000,   10000,  220000,         30000,         30000,    450000,1000 };

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

