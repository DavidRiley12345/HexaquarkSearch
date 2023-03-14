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

// use this order for safety on library loading
using namespace RooFit;
using namespace RooStats;

// see below for implementation

void gensweight_simpler(){
    ROOT::EnableImplicitMT();

    // Create the model of the signal and background distributions of the sigmac mass
    // using the DstDOBG and breitwigner functions

    // Create a RooRealVar for the mass variable

    RooRealVar sig_M("sig_M","",2350.0,3000.0);

    // create signal pdfs

    RooRealVar sig1_mean("sig1_mean","",2455,2450,2460);
    RooRealVar sig1_width("sig1_width","",5,0.1,15);
    RooBreitWigner sig1("sig1","",sig_M,sig1_mean,sig1_width);

    RooRealVar sig2_mean("sig2_mean","",2520,2500,2540);
    RooRealVar sig2_width("sig2_width","",5,0.1,25);
    RooBreitWigner sig2("sig2","",sig_M,sig2_mean,sig2_width);

    // combine the signal pdf's

    RooRealVar sig1_frac("sig1_frac","",0.5,0,1);
    RooRealVar sig2_frac("sig2_frac","",0.1,0,1);
    RooAddPdf sigPDF("sigPDF","",RooArgList(sig1,sig2),RooArgList(sig1_frac,sig2_frac));

    // Now we want to create the background PDF but we will do this by using SPlot weighting on the Lambda_C variable
    // and using the background weighted Lambda_C variable to create a PDF

    //first we import the sweighted data

    TFile* fin = TFile::Open("/home/ppe/d/driley/git_HexaquarkSummer/Sweight/allvars/weighted.root");
    TTree* tree = (TTree*)fin->Get("weightedtree");

    // Create the RooDataSet from the tree
    RooRealVar nbkg_sw("nbkg_sw","",-10,10);
    RooRealVar nsig_sw("nsig_sw","",-10,10);
    RooArgSet argset1(sig_M,nbkg_sw);
    RooArgSet argset2(sig_M,nsig_sw);
    RooDataSet bkgdata("bkgdata","",argset1,Import(*tree),WeightVar("nbkg_sw"));
    RooDataSet sigdata("sigdata","",argset2,Import(*tree),WeightVar("nsig_sw"));

    // We now have signal and background data weighted in these RooDataSets
    // We want to use RooHistPdf to create a PDF from the weighted data

    // TH1* hist = bkgdata.createHistogram("hist",sig_M,Binning(30));

    // RooDataHist *bkgHist = new RooDataHist("bkgHist","",RooArgList(sig_M),hist);

    // RooHistPdf bkgPDF("bkgPDF","",RooArgSet(sig_M),*bkgHist,2.0);

    // We might prefer to use RooKeysPdf to create a PDF from the weighted data

    RooKeysPdf bkgPDF("bkgPDF2","",sig_M,bkgdata,RooKeysPdf::MirrorBoth,2.0);   

    // lets output them to make sure they look ok

    TCanvas* c1 = new TCanvas("c1","c1",800,600);
    RooPlot* frame = sig_M.frame();
    bkgdata.plotOn(frame,MarkerColor(kRed),Binning(40));
    sigdata.plotOn(frame,MarkerColor(kGreen),Binning(40));
    bkgPDF.plotOn(frame);
    frame->Draw();
    c1->SaveAs("bkgPDFfitting.png");
  
    // Now we have the signal and background PDFs we can create the model for our distribution
    // first lets create our signal and bkg yeilds

    RooRealVar sig_yield("sig_yield","",1000,0,100000);
    RooRealVar bkg_yield("bkg_yield","",100000,0,100000);

    RooAddPdf model("model","",RooArgList(sigPDF,bkgPDF),RooArgList(sig_yield,bkg_yield));

    // Lets import our data

    TFile* fin2 = TFile::Open("/home/ppe/d/driley/git_HexaquarkSummer/Sweight/sigmac_sw/sigmac.root");
    TTree* tree2 = (TTree*)fin2->Get("sigmac");

    RooDataSet datain("datain","",tree2,RooArgSet(sig_M));

    // Now we can fit the model to the data

    model.fitTo(datain,Extended());

    // check the output

    TCanvas* c2 = new TCanvas("c2","c2",800,600);
    RooPlot* frame2 = sig_M.frame();
    datain.plotOn(frame2);
    model.plotOn(frame2);
    model.plotOn(frame2,Components(bkgPDF),LineStyle(kDashed),LineColor(kRed));
    model.plotOn(frame2,Components(sigPDF),LineStyle(kDashed),LineColor(kGreen));
    frame2->Draw();
    c2->SaveAs("sigc_M_fit.png");
}