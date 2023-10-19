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
#include "TLorentzVector.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooConstVar.h"
#include "ROOT/RDataFrame.hxx"
#include "TTree.h"

// this file will output the two and three particle backgrounds of the decay products of the Lc's

void particleBGs(){

    ROOT::EnableImplicitMT();

    // define names for colours
    Int_t black  = 1;
    Int_t red    = 2;
    Int_t green  = 3;
    Int_t blue   = 4;
    Int_t yellow = 5; 
    Int_t magenta= 6;
    Int_t cyan   = 7;
    Int_t purple = 9;
    
    // Use times new roman, precision 2 
    Int_t lhcbFont        = 132;  // Old LHCb style: 62;
    // Line thickness
    Double_t lhcbWidth    = 2.00; // Old LHCb style: 3.00;
    // Text size
    Double_t lhcbTSize    = 0.06; 
    
    // use plain black on white colors
    gROOT->SetStyle("Plain"); 
    TStyle *lhcbStyle= new TStyle("lhcbStyle","LHCb plots style");
    
    //lhcbStyle->SetErrorX(0); //  don't suppress the error bar along X

    lhcbStyle->SetFillColor(1);
    lhcbStyle->SetFillStyle(1001);   // solid
    lhcbStyle->SetFrameFillColor(0);
    lhcbStyle->SetFrameBorderMode(0);
    lhcbStyle->SetPadBorderMode(0);
    lhcbStyle->SetPadColor(0);
    lhcbStyle->SetCanvasBorderMode(0);
    lhcbStyle->SetCanvasColor(0);
    lhcbStyle->SetStatColor(0);
    lhcbStyle->SetLegendBorderSize(0);

    // If you want the usual gradient palette (blue -> red)
    lhcbStyle->SetPalette(1);
    // If you want colors that correspond to gray scale in black and white:
    //int colors[8] = {0,5,7,3,6,2,4,1};
    //lhcbStyle->SetPalette(8,colors);
    //lhcbStyle->SetPalette(1);

    // set the paper & margin sizes
    lhcbStyle->SetPaperSize(20,26);
    lhcbStyle->SetPadTopMargin(0.07);
    lhcbStyle->SetPadRightMargin(0.05); // increase for colz plots
    lhcbStyle->SetPadBottomMargin(0.16);
    lhcbStyle->SetPadLeftMargin(0.18);
    
    // use large fonts
    lhcbStyle->SetLegendFont(lhcbFont);
    lhcbStyle->SetTextFont(lhcbFont);
    lhcbStyle->SetTextSize(lhcbTSize);
    lhcbStyle->SetLabelFont(lhcbFont,"x");
    lhcbStyle->SetLabelFont(lhcbFont,"y");
    lhcbStyle->SetLabelFont(lhcbFont,"z");
    lhcbStyle->SetLabelSize(lhcbTSize,"x");
    lhcbStyle->SetLabelSize(lhcbTSize,"y");
    lhcbStyle->SetLabelSize(lhcbTSize,"z");
    lhcbStyle->SetTitleFont(lhcbFont);
    lhcbStyle->SetTitleFont(lhcbFont,"x");
    lhcbStyle->SetTitleFont(lhcbFont,"y");
    lhcbStyle->SetTitleFont(lhcbFont,"z");
    lhcbStyle->SetTitleSize(1.2*lhcbTSize,"x");
    lhcbStyle->SetTitleSize(1.2*lhcbTSize,"y");
    lhcbStyle->SetTitleSize(1.2*lhcbTSize,"z");

    // use medium bold lines and thick markers
    lhcbStyle->SetLineWidth(lhcbWidth);
    lhcbStyle->SetFrameLineWidth(lhcbWidth);
    lhcbStyle->SetHistLineWidth(lhcbWidth);
    lhcbStyle->SetFuncWidth(lhcbWidth);
    lhcbStyle->SetGridWidth(lhcbWidth);
    lhcbStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
    lhcbStyle->SetMarkerStyle(20);
    lhcbStyle->SetMarkerSize(1.0);

    // label offsets
    lhcbStyle->SetLabelOffset(0.010,"X");
    lhcbStyle->SetLabelOffset(0.010,"Y");

    // by default, do not display histogram decorations:
    lhcbStyle->SetOptStat(0);  
    //lhcbStyle->SetOptStat("emr");  // show only nent -e , mean - m , rms -r
    // full opts at http://root.cern.ch/root/html/TStyle.html#TStyle:SetOptStat
    lhcbStyle->SetStatFormat("6.3g"); // specified as c printf options
    lhcbStyle->SetOptTitle(1);
    lhcbStyle->SetOptFit(0);
    //lhcbStyle->SetOptFit(1011); // order is probability, Chi2, errors, parameters
    //titles
    lhcbStyle->SetTitleOffset(0.95,"X");
    lhcbStyle->SetTitleOffset(1.3,"Y");
    lhcbStyle->SetTitleOffset(1.2,"Z");
    lhcbStyle->SetTitleFillColor(0);
    lhcbStyle->SetTitleStyle(0);
    lhcbStyle->SetTitleBorderSize(0);
    lhcbStyle->SetTitleFont(lhcbFont,"title");
    lhcbStyle->SetTitleX(0.0);
    lhcbStyle->SetTitleY(1.0); 
    lhcbStyle->SetTitleW(1.0);
    lhcbStyle->SetTitleH(0.05);
    
    // look of the statistics box:
    lhcbStyle->SetStatBorderSize(0);
    lhcbStyle->SetStatFont(lhcbFont);
    lhcbStyle->SetStatFontSize(0.05);
    lhcbStyle->SetStatX(0.9);
    lhcbStyle->SetStatY(0.9);
    lhcbStyle->SetStatW(0.25);
    lhcbStyle->SetStatH(0.15);

    // put tick marks on top and RHS of plots
    lhcbStyle->SetPadTickX(1);
    lhcbStyle->SetPadTickY(1);

    // histogram divisions: only 5 in x to avoid label overlaps
    lhcbStyle->SetNdivisions(505,"x");
    lhcbStyle->SetNdivisions(510,"y");
    
    gROOT->SetStyle("lhcbStyle");
    gROOT->ForceStyle();

    auto d = ROOT::RDataFrame("B2LcLcpiSS/DecayTree", {"/data/lhcb01/mwhitehead/LcLcpi_2018_MU.root","/data/lhcb01/mwhitehead/LcLcpi_2018_MD.root","/data/lhcb01/mwhitehead/LcLcpi_2017_MU.root","/data/lhcb01/mwhitehead/LcLcpi_2017_MD.root"});
    auto df = d.Filter("Lambdacp_p_ProbNNp > 0.7 && Lambdacp_K_ProbNNk > 0.6 && Lambdacp_pi_ProbNNpi > 0.4 && Lambdacm_p_ProbNNp > 0.7 && Lambdacm_K_ProbNNk > 0.6 && Lambdacm_pi_ProbNNpi > 0.4");
    auto df1 = df.Define("protonp",    "TLorentzVector(Lambdacp_p_PX, Lambdacp_p_PY, Lambdacp_p_PZ, Lambdacp_p_PE)")
                 .Define("kaonp",    "TLorentzVector(Lambdacp_K_PX, Lambdacp_K_PY, Lambdacp_K_PZ, Lambdacp_K_PE)")
                 .Define("pionp",  "TLorentzVector(Lambdacp_pi_PX, Lambdacp_pi_PY, Lambdacp_pi_PZ, Lambdacp_pi_PE)")
                 .Define("protonm",    "TLorentzVector(Lambdacm_p_PX, Lambdacm_p_PY, Lambdacm_p_PZ, Lambdacm_p_PE)")
                 .Define("kaonm",    "TLorentzVector(Lambdacm_K_PX, Lambdacm_K_PY, Lambdacm_K_PZ, Lambdacm_K_PE)")
                 .Define("pionm",  "TLorentzVector(Lambdacm_pi_PX, Lambdacm_pi_PY, Lambdacm_pi_PZ, Lambdacm_pi_PE)");

    // TH1F* h_pion_pion = new TH1F("h_pion_pion", "Invariant mass of pion-pion pairs", 100, 0, 5);
    // TH1F* h_proton_kaon = new TH1F("h_proton_kaon", "Invariant mass of proton-kaon pairs", 100, 0, 5);
    // TH1F* h_proton_pion = new TH1F("h_proton_pion", "Invariant mass of proton-pion pairs", 100, 0, 5);
    // TH1F* h_kaon_pion = new TH1F("h_kaon_pion", "Invariant mass of kaon-pion pairs", 100, 0, 5);
    // TH1F* h_kaonm_kaonp = new TH1F("h_kaonm_kaonp", "Invariant mass of kaonm-kaonp pairs", 100, 0, 5);

    // Calculate invariant masses for all sub-particle pairs
    auto df2 = df1.Define("inv_mass_pion_pion","(pionp+pionm).M()")
                  .Define("inv_mass_proton_kaon","(protonp+kaonm).M()")
                  .Define("inv_mass_proton_pion","(protonp+pionm).M()")
                  .Define("inv_mass_kaon_pion","(kaonp+pionm).M()")
                  .Define("inv_mass_kaonm_kaonp","(kaonm+kaonp).M()")
                  .Define("inv_mass_pKmpim","(protonp+kaonm+pionm).M()");

    // Fill histograms with invariant mass values from dataframe
    auto h1 = df2.Histo1D({"pion-pion background", "pion-pion background", 25, 200, 3000}, "inv_mass_pion_pion");
    auto h2 = df2.Histo1D({"proton-kaon background", "proton-kaon background", 25, 1400, 3000}, "inv_mass_proton_kaon");
    auto h3 = df2.Histo1D({"proton-pion background", "proton-pion background", 25, 1000, 3000}, "inv_mass_proton_pion");
    auto h4 = df2.Histo1D({"kaon-pion background", "kaon-pion background", 25, 600, 3000}, "inv_mass_kaon_pion");
    auto h5 = df2.Histo1D({"kaonm-kaonp background", "kaonm-kaonp background", 25, 900, 2500}, "inv_mass_kaonm_kaonp");
    auto h6 = df2.Histo1D({"pKmpim background", "pKmpim background", 25, 1800, 4000}, "inv_mass_pKmpim");

    TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
    h1->SetTitle("Invariant mass of pion-pion pairs (SS)");
    h1->SetXTitle("m(#pi^{+} #pi^{-}) [MeV]");
    h1->SetYTitle(("Candidates / " + std::to_string(static_cast<int>(h1->GetBinWidth(1)*1000)/1000) + " MeV").c_str());
    h1->Draw("pe");
    c1->SaveAs("pion_pion.png");

    TCanvas* c2 = new TCanvas("c2", "c2", 800, 600);
    h2->SetTitle("Invariant mass of proton-kaon pairs (SS)");
    h2->SetXTitle("m(pK) [MeV]");
    h2->SetYTitle(("Candidates / " + std::to_string(static_cast<int>(h2->GetBinWidth(1)*1000)/1000) + " MeV").c_str());
    h2->Draw("pe");
    c2->SaveAs("proton_kaon.png");

    TCanvas* c3 = new TCanvas("c3", "c3", 800, 600);
    h3->SetTitle("Invariant mass of proton-pion pairs (SS)");
    h3->SetXTitle("m(p #pi) [MeV]");
    h3->SetYTitle(("Candidates / " + std::to_string(static_cast<int>(h3->GetBinWidth(1)*1000)/1000) + " MeV").c_str());
    h3->Draw("pe");
    c3->SaveAs("proton_pion.png");

    TCanvas* c4 = new TCanvas("c4", "c4", 800, 600);
    h4->SetTitle("Invariant mass of kaon-pion pairs (SS)");
    h4->SetYTitle(("Candidates / " + std::to_string(static_cast<int>(h4->GetBinWidth(1)*1000)/1000) + " MeV").c_str());
    h4->SetXTitle("m(K #pi) [MeV]");
    h4->Draw("pe");
    c4->SaveAs("kaon_pion.png");

    TCanvas* c5 = new TCanvas("c5", "c5", 800, 600);
    h5->SetTitle("Invariant mass of kaonm-kaonp pairs (SS)");
    h5->SetYTitle(("Candidates / " + std::to_string(static_cast<int>(h5->GetBinWidth(1)*1000)/1000) + " MeV").c_str());
    h5->SetXTitle("m(K^{+}K^{-}) [MeV]");
    h5->Draw("pe");
    c5->SaveAs("kaonm_kaonp.png");

    TCanvas* c6 = new TCanvas("c6", "c6", 800, 600);
    h6->SetTitle("Invariant mass of pKmpim pairs (SS)");
    h6->SetYTitle(("Candidates / " + std::to_string(static_cast<int>(h6->GetBinWidth(1)*1000)/1000) + " MeV").c_str());
    h6->SetXTitle("m(pK^{-}#pi^{-}) [MeV]");
    h6->Draw("pe");
    c6->SaveAs("pKmpim.png");
};

