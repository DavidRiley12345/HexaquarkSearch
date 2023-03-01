#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDstD0BG.h"
#include "RooAddPdf.h"
#include "TCanvas.h"
#include "RooFit.h"
#include "TTree.h"

using namespace RooFit;

void sigc_keyspdf() {
    ROOT::EnableImplicitMT();
    // Create a RooRealVar for the mass variable
    RooRealVar sig_M("sig_M","",2350.0,3000.0);

    // Create background PDF variables
    RooRealVar dm0("dm0","",2360,2350,2500);
    RooRealVar a("a","",800,0,10000);
    RooRealVar b("b","",-8,-20,4);
    RooRealVar c("c","",-0.005,-6,4);

    // Create the background PDF
    RooDstD0BG bkgPDF("bkgPDF","",sig_M,dm0,a,b,c);
   
    // Create the normalization variable for the background PDF
    RooRealVar nbkg("nbkg","",1,1000000);

    // Create the total PDF
    RooAddPdf model("model","",RooArgList(bkgPDF),RooArgList(nbkg));

    // Open the data file and get the tree
    TFile* fin = TFile::Open("/home/ppe/d/driley/git_HexaquarkSummer/Sweight/allvars/weighted.root"); 
    TTree* tree = (TTree*)fin->Get("weightedtree");

    // Create the RooDataSet from the tree
    RooRealVar nbkg_sw("nbkg_sw","",-10,10);
    RooRealVar nsig_sw("nsig_sw","",-10,10);
    RooArgSet argset(sig_M,nbkg_sw);
    RooArgSet argset2(sig_M,nsig_sw);
    RooDataSet data("data", "Sigc bkg weighted",argset,Import(*tree),WeightVar("nbkg_sw"));
    RooDataSet sigdata("sigdata","Sigc sig weighted",argset2,Import(*tree),WeightVar("nsig_sw"));
    // Fit the model to the data
    model.fitTo(data, Extended());
   
    RooKeysPdf keysPDF("keysPDF","",sig_M,data);

    RooWorkspace workspace("workspace");
    workspace.import(keysPDF);
   
    TFile* outfile = new TFile("/home/ppe/d/driley/git_HexaquarkSummer/Sweight/sigc_bkg_fit/pdf.root","RECREATE");
    workspace.Write();
    outfile->Close();

    TH1* hist = data.createHistogram("hist",sig_M,Binning(50));
    hist->Scale(0.5/ hist->Integral());

    //TFile* outhist = new TFile("/home/ppe/d/driley/git_HexaquarkSummer/Sweight/sigc_bkg_fit/hist.root","RECREATE"); 
    data.write("hist.root");
    //outhist->Close(); 

    // Plot the data and the fit
    TCanvas *c1 = new TCanvas("c1","",800,600);
    RooPlot* frame = sig_M.frame();
    //data.plotOn(frame,MarkerColor(kRed),Normalization(1));
    //sigdata.plotOn(frame,MarkerColor(kGreen),Normalization(1));
    model.plotOn(frame);
        
    keysPDF.plotOn(frame);
    frame->addTH1(hist);
    frame->Draw();
    frame->SetTitle("sigmac bkg fit");
    frame->SetXTitle("sigmac invariant mass");
    c1->SaveAs("sigma_bkg_fit.gif");
}
