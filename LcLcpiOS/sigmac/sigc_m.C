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

void sigc_m(){

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
    
    // plotting the invariant mass of the sigmac with the full set of cuts applied

    ROOT::EnableImplicitMT();

    ROOT::RDataFrame df("B2LcLcpiOS/DecayTree", {"/data/lhcb01/mwhitehead/LcLcpi_2018_MU.root","/data/lhcb01/mwhitehead/LcLcpi_2018_MD.root","/data/lhcb01/mwhitehead/LcLcpi_2017_MU.root","/data/lhcb01/mwhitehead/LcLcpi_2017_MD.root"});
    auto df1 = df.Filter("Lambdacp_M < 2305 && Lambdacp_M > 2270")
    .Filter("Lambdacp_p_ProbNNp > 0.7 && Lambdacp_K_ProbNNk > 0.6 && Lambdacp_pi_ProbNNpi > 0.4")
    .Define("Angle_pp_pip", "(180/3.14)*(acos((Lambdacp_pi_PX*Lambdacp_p_PX + Lambdacp_pi_PY*Lambdacp_p_PY + Lambdacp_pi_PZ*Lambdacp_p_PZ) / (sqrt(pow(Lambdacp_pi_PX,2)+pow(Lambdacp_pi_PY,2)+pow(Lambdacp_pi_PZ,2)) * sqrt(pow(Lambdacp_p_PX,2)+pow(Lambdacp_p_PY,2)+pow(Lambdacp_p_PZ,2)))))")
    .Define("Angle_pm_pim", "(180/3.14)*(acos((Lambdacm_pi_PX*Lambdacm_p_PX + Lambdacm_pi_PY*Lambdacm_p_PY + Lambdacm_pi_PZ*Lambdacm_p_PZ) / (sqrt(pow(Lambdacm_pi_PX,2)+pow(Lambdacm_pi_PY,2)+pow(Lambdacm_pi_PZ,2)) * sqrt(pow(Lambdacm_p_PX,2)+pow(Lambdacm_p_PY,2)+pow(Lambdacm_p_PZ,2)))))")
    .Define("Angle_kp_pp", "(180/3.14)*(acos((Lambdacm_K_PX*Lambdacp_p_PX + Lambdacm_K_PY*Lambdacp_p_PY + Lambdacm_K_PZ*Lambdacp_p_PZ) / (sqrt(pow(Lambdacm_K_PX,2)+pow(Lambdacm_K_PY,2)+pow(Lambdacm_K_PZ,2)) * sqrt(pow(Lambdacp_p_PX,2)+pow(Lambdacp_p_PY,2)+pow(Lambdacp_p_PZ,2)))))")
    .Define("Angle_km_pm", "(180/3.14)*(acos((Lambdacp_K_PX*Lambdacm_p_PX + Lambdacp_K_PY*Lambdacm_p_PY + Lambdacp_K_PZ*Lambdacm_p_PZ) / (sqrt(pow(Lambdacp_K_PX,2)+pow(Lambdacp_K_PY,2)+pow(Lambdacp_K_PZ,2)) * sqrt(pow(Lambdacm_p_PX,2)+pow(Lambdacm_p_PY,2)+pow(Lambdacm_p_PZ,2)))))")
    .Define("Angle_kp_pip", "(180/3.14)*(acos((Lambdacm_K_PX*Lambdacp_pi_PX + Lambdacm_K_PY*Lambdacp_pi_PY + Lambdacm_K_PZ*Lambdacp_pi_PZ) / (sqrt(pow(Lambdacm_K_PX,2)+pow(Lambdacm_K_PY,2)+pow(Lambdacm_K_PZ,2)) * sqrt(pow(Lambdacp_pi_PX,2)+pow(Lambdacp_pi_PY,2)+pow(Lambdacp_pi_PZ,2)))))")
    .Filter("Angle_pp_pip >0.05 && Angle_pm_pim >0.05 && Angle_kp_pp >0.05 && Angle_km_pm >0.05 && Angle_kp_pip >0.05")
    .Filter("Lambdacp_K_ORIVX_CHI2 < 10 && Lambdacp_p_ORIVX_CHI2 < 10 && Lambdacp_pi_ORIVX_CHI2 < 10")
    .Filter("Lambdacm_K_ORIVX_CHI2 < 10 && Lambdacm_p_ORIVX_CHI2 < 10 && Lambdacm_pi_ORIVX_CHI2 < 10");
    // .Filter("Lambdacm_M < 2305 && Lambdacm_M > 2270")
    // .Filter("Lambdacm_p_ProbNNp > 0.7 && Lambdacm_K_ProbNNk > 0.6 && Lambdacm_pi_ProbNNpi > 0.4")
    // .Filter("LcLcpiOS_ENDVERTEX_CHI2 < 9")

    auto df2 = df1.Define("sigc_M","sqrt(pow(Lambdacp_PE + pi_PE,2) - pow(Lambdacp_PX + pi_PX,2) - pow(Lambdacp_PY + pi_PY,2) - pow(Lambdacp_PZ + pi_PZ,2))");

    auto h1 = df1.Histo1D({"h1", "h1", 100, 2200, 2320}, "Lambdacp_M");
    auto h2 = df2.Histo1D({"#Sigma_{c} Invariant Mass Plot", "#Sigma_{c} Invariant Mass Plot", 50, 2400, 3000}, "sigc_M");

    // plot the histogram

    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    h1->Draw();
    c1->SaveAs("lc_m.png");

    // plot the second histogram

    TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
    
    // set the yaxis title to show events per bin width calculating bin width

    h2->SetYTitle(("Candidates / " + std::to_string(static_cast<int>(h2->GetBinWidth(1)*1000)/1000) + " MeV").c_str());
    h2->SetXTitle("m(#Lambda_{c}#pi) [MeV/c^{2}]");
    h2->Draw();
    c2->SaveAs("sigc_m.png");

};