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

// this file will output the PID of the proton, Kaon, and pion

void mom_angle() 
{

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
    lhcbStyle->SetPadTopMargin(0.1);
    lhcbStyle->SetPadRightMargin(0.05); // increase for colz plots
    lhcbStyle->SetPadBottomMargin(0.18);
    lhcbStyle->SetPadLeftMargin(0.19);
    
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
    lhcbStyle->SetLabelOffset(0.01,"Y");

    // by default, do not display histogram decorations:
    lhcbStyle->SetOptStat(0);  
    //lhcbStyle->SetOptStat("emr");  // show only nent -e , mean - m , rms -r
    // full opts at http://root.cern.ch/root/html/TStyle.html#TStyle:SetOptStat
    lhcbStyle->SetStatFormat("6.3g"); // specified as c printf options
    lhcbStyle->SetOptTitle(1);
    lhcbStyle->SetOptFit(0);
    //lhcbStyle->SetOptFit(1011); // order is probability, Chi2, errors, parameters
    //titles
    lhcbStyle->SetTitleOffset(0.5,"T");
    lhcbStyle->SetTitleOffset(0.95,"X");
    lhcbStyle->SetTitleOffset(1.2,"Y");
    lhcbStyle->SetTitleOffset(1.2,"Z");
    lhcbStyle->SetTitleFillColor(0);
    lhcbStyle->SetTitleStyle(0);
    lhcbStyle->SetTitleBorderSize(0);
    lhcbStyle->SetTitleFont(lhcbFont,"title");
    lhcbStyle->SetTitleFontSize(10);
    lhcbStyle->SetTitleX(0);
    lhcbStyle->SetTitleY(0.9); 
    lhcbStyle->SetTitleW(1.0);
    lhcbStyle->SetTitleH(0.1);
    
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

    lhcbStyle->SetTitleFontSize(0.1);
    
    gROOT->SetStyle("lhcbStyle");
    gROOT->ForceStyle();
      

    // create a dataframe

    ROOT::RDataFrame df("B2LcLcpiOS/DecayTree", {"/data/lhcb01/mwhitehead/LcLcpi_2018_MU.root","/data/lhcb01/mwhitehead/LcLcpi_2018_MD.root","/data/lhcb01/mwhitehead/LcLcpi_2017_MU.root","/data/lhcb01/mwhitehead/LcLcpi_2017_MD.root"});

    auto dffiltered = df.Filter("Lambdacp_p_ProbNNp > 0.7 && Lambdacp_K_ProbNNk > 0.6 && Lambdacp_pi_ProbNNpi > 0.4 && Lambdacm_p_ProbNNp > 0.7 && Lambdacm_K_ProbNNk > 0.6 && Lambdacm_pi_ProbNNpi > 0.4");

    //this script is going to calculate the angle between the momentum vector of the decay products of the Lambdac's
    //first we will define columns in the array which calculate the total angle between the 

    auto df1 = dffiltered.Define("angle_pp_pm",  "(180/3.14)*(acos((Lambdacp_p_PX*Lambdacm_p_PX   + Lambdacp_p_PY*Lambdacm_p_PY   + Lambdacp_p_PZ*Lambdacm_p_PZ)   / (sqrt(pow(Lambdacp_p_PX,2)+pow(Lambdacp_p_PY,2)+pow(Lambdacp_p_PZ,2))    * sqrt(pow(Lambdacm_p_PX,2)+pow(Lambdacm_p_PY,2)+pow(Lambdacm_p_PZ,2   )))))")
                 .Define("angle_pp_Km",  "(180/3.14)*(acos((Lambdacp_p_PX*Lambdacm_K_PX   + Lambdacp_p_PY*Lambdacm_K_PY   + Lambdacp_p_PZ*Lambdacm_K_PZ)   / (sqrt(pow(Lambdacp_p_PX,2)+pow(Lambdacp_p_PY,2)+pow(Lambdacp_p_PZ,2))    * sqrt(pow(Lambdacm_K_PX,2)+pow(Lambdacm_K_PY,2)+pow(Lambdacm_K_PZ,2   )))))")
                 .Define("angle_pp_pim", "(180/3.14)*(acos((Lambdacp_p_PX*Lambdacm_pi_PX  + Lambdacp_p_PY*Lambdacm_pi_PY  + Lambdacp_p_PZ*Lambdacm_pi_PZ)  / (sqrt(pow(Lambdacp_p_PX,2)+pow(Lambdacp_p_PY,2)+pow(Lambdacp_p_PZ,2))    * sqrt(pow(Lambdacm_pi_PX,2)+pow(Lambdacm_pi_PY,2)+pow(Lambdacm_pi_PZ,2)))))")
                 .Define("angle_pp_pip", "(180/3.14)*(acos((Lambdacp_p_PX*Lambdacp_pi_PX  + Lambdacp_p_PY*Lambdacp_pi_PY  + Lambdacp_p_PZ*Lambdacp_pi_PZ)  / (sqrt(pow(Lambdacp_p_PX,2)+pow(Lambdacp_p_PY,2)+pow(Lambdacp_p_PZ,2))    * sqrt(pow(Lambdacp_pi_PX,2)+pow(Lambdacp_pi_PY,2)+pow(Lambdacp_pi_PZ,2)))))")
                 .Define("angle_pp_Kp",  "(180/3.14)*(acos((Lambdacp_p_PX*Lambdacp_K_PX   + Lambdacp_p_PY*Lambdacp_K_PY   + Lambdacp_p_PZ*Lambdacp_K_PZ)   / (sqrt(pow(Lambdacp_p_PX,2)+pow(Lambdacp_p_PY,2)+pow(Lambdacp_p_PZ,2))    * sqrt(pow(Lambdacp_K_PX,2)+pow(Lambdacp_K_PY,2)+pow(Lambdacp_K_PZ,2   )))))")
                 .Define("angle_Kp_pm",  "(180/3.14)*(acos((Lambdacp_K_PX*Lambdacm_p_PX   + Lambdacp_K_PY*Lambdacm_p_PY   + Lambdacp_K_PZ*Lambdacm_p_PZ)   / (sqrt(pow(Lambdacp_K_PX,2)+pow(Lambdacp_K_PY,2)+pow(Lambdacp_K_PZ,2))    * sqrt(pow(Lambdacm_p_PX,2)+pow(Lambdacm_p_PY,2)+pow(Lambdacm_p_PZ,2   )))))")
                 .Define("angle_Kp_Km",  "(180/3.14)*(acos((Lambdacp_K_PX*Lambdacm_K_PX   + Lambdacp_K_PY*Lambdacm_K_PY   + Lambdacp_K_PZ*Lambdacm_K_PZ)   / (sqrt(pow(Lambdacp_K_PX,2)+pow(Lambdacp_K_PY,2)+pow(Lambdacp_K_PZ,2))    * sqrt(pow(Lambdacm_K_PX,2)+pow(Lambdacm_K_PY,2)+pow(Lambdacm_K_PZ,2   )))))")
                 .Define("angle_Kp_pim", "(180/3.14)*(acos((Lambdacp_K_PX*Lambdacm_pi_PX  + Lambdacp_K_PY*Lambdacm_pi_PY  + Lambdacp_K_PZ*Lambdacm_pi_PZ)  / (sqrt(pow(Lambdacp_K_PX,2)+pow(Lambdacp_K_PY,2)+pow(Lambdacp_K_PZ,2))    * sqrt(pow(Lambdacm_pi_PX,2)+pow(Lambdacm_pi_PY,2)+pow(Lambdacm_pi_PZ,2)))))")
                 .Define("angle_Kp_pip", "(180/3.14)*(acos((Lambdacp_K_PX*Lambdacp_pi_PX  + Lambdacp_K_PY*Lambdacp_pi_PY  + Lambdacp_K_PZ*Lambdacp_pi_PZ)  / (sqrt(pow(Lambdacp_K_PX,2)+pow(Lambdacp_K_PY,2)+pow(Lambdacp_K_PZ,2))    * sqrt(pow(Lambdacp_pi_PX,2)+pow(Lambdacp_pi_PY,2)+pow(Lambdacp_pi_PZ,2)))))")
                 .Define("angle_pip_pm", "(180/3.14)*(acos((Lambdacp_pi_PX*Lambdacm_p_PX  + Lambdacp_pi_PY*Lambdacm_p_PY  + Lambdacp_pi_PZ*Lambdacm_p_PZ)  / (sqrt(pow(Lambdacp_pi_PX,2)+pow(Lambdacp_pi_PY,2)+pow(Lambdacp_pi_PZ,2)) * sqrt(pow(Lambdacm_p_PX,2)+pow(Lambdacm_p_PY,2)+pow(Lambdacm_p_PZ,2   )))))")
                 .Define("angle_pip_Km", "(180/3.14)*(acos((Lambdacp_pi_PX*Lambdacm_K_PX  + Lambdacp_pi_PY*Lambdacm_K_PY  + Lambdacp_pi_PZ*Lambdacm_K_PZ)  / (sqrt(pow(Lambdacp_pi_PX,2)+pow(Lambdacp_pi_PY,2)+pow(Lambdacp_pi_PZ,2)) * sqrt(pow(Lambdacm_K_PX,2)+pow(Lambdacm_K_PY,2)+pow(Lambdacm_K_PZ,2   )))))")
                 .Define("angle_pip_pim","(180/3.14)*(acos((Lambdacp_pi_PX*Lambdacm_pi_PX + Lambdacp_pi_PY*Lambdacm_pi_PY + Lambdacp_pi_PZ*Lambdacm_pi_PZ) / (sqrt(pow(Lambdacp_pi_PX,2)+pow(Lambdacp_pi_PY,2)+pow(Lambdacp_pi_PZ,2)) * sqrt(pow(Lambdacm_pi_PX,2)+pow(Lambdacm_pi_PY,2)+pow(Lambdacm_pi_PZ,2)))))")
                 .Define("angle_pim_pm", "(180/3.14)*(acos((Lambdacm_pi_PX*Lambdacm_p_PX  + Lambdacm_pi_PY*Lambdacm_p_PY  + Lambdacm_pi_PZ*Lambdacm_p_PZ)  / (sqrt(pow(Lambdacm_pi_PX,2)+pow(Lambdacm_pi_PY,2)+pow(Lambdacm_pi_PZ,2)) * sqrt(pow(Lambdacm_p_PX,2)+pow(Lambdacm_p_PY,2)+pow(Lambdacm_p_PZ,2   )))))")
                 .Define("angle_pim_Km", "(180/3.14)*(acos((Lambdacm_pi_PX*Lambdacm_K_PX  + Lambdacm_pi_PY*Lambdacm_K_PY  + Lambdacm_pi_PZ*Lambdacm_K_PZ)  / (sqrt(pow(Lambdacm_pi_PX,2)+pow(Lambdacm_pi_PY,2)+pow(Lambdacm_pi_PZ,2)) * sqrt(pow(Lambdacm_K_PX,2)+pow(Lambdacm_K_PY,2)+pow(Lambdacm_K_PZ,2   )))))")
                 .Define("angle_Km_pm",  "(180/3.14)*(acos((Lambdacm_K_PX*Lambdacm_p_PX   + Lambdacm_K_PY*Lambdacm_p_PY   + Lambdacm_K_PZ*Lambdacm_p_PZ)   / (sqrt(pow(Lambdacm_K_PX,2)+pow(Lambdacm_K_PY,2)+pow(Lambdacm_K_PZ,2))    * sqrt(pow(Lambdacm_p_PX,2)+pow(Lambdacm_p_PY,2)+pow(Lambdacm_p_PZ,2   )))))");

    // now we want to plot each of these angles using the Histo1D function

    auto h_angle_pp_pm = df1.Histo1D({"h_angle_pp_pm","h_angle_pp_pm",20,0,0.3},"angle_pp_pm");
    auto h_angle_pp_Km = df1.Histo1D({"h_angle_pp_Km","h_angle_pp_Km",20,0,0.3},"angle_pp_Km");
    auto h_angle_pp_pim = df1.Histo1D({"h_angle_pp_pim","h_angle_pp_pim",20,0,0.3},"angle_pp_pim");
    auto h_angle_pp_pip = df1.Histo1D({"h_angle_pp_pip","h_angle_pp_pip",20,0,0.3},"angle_pp_pip");
    auto h_angle_pp_Kp = df1.Histo1D({"h_angle_pp_Kp","h_angle_pp_Kp",20,0,0.3},"angle_pp_Kp");
    auto h_angle_Kp_pm = df1.Histo1D({"h_angle_Kp_pm","h_angle_Kp_pm",20,0,0.3},"angle_Kp_pm");
    auto h_angle_Kp_Km = df1.Histo1D({"h_angle_Kp_Km","h_angle_Kp_Km",20,0,0.3},"angle_Kp_Km");
    auto h_angle_Kp_pim = df1.Histo1D({"h_angle_Kp_pim","h_angle_Kp_pim",20,0,0.3},"angle_Kp_pim");
    auto h_angle_Kp_pip = df1.Histo1D({"h_angle_Kp_pip","h_angle_Kp_pip",20,0,0.3},"angle_Kp_pip");
    auto h_angle_pip_pm = df1.Histo1D({"h_angle_pip_pm","h_angle_pip_pm",20,0,0.3},"angle_pip_pm");
    auto h_angle_pip_Km = df1.Histo1D({"h_angle_pip_Km","h_angle_pip_Km",20,0,0.3},"angle_pip_Km");
    auto h_angle_pip_pim = df1.Histo1D({"h_angle_pip_pim","h_angle_pip_pim",20,0,0.3},"angle_pip_pim");
    auto h_angle_pim_pm = df1.Histo1D({"h_angle_pim_pm","h_angle_pim_pm",20,0,0.3},"angle_pim_pm");
    auto h_angle_pim_Km = df1.Histo1D({"h_angle_pim_Km","h_angle_pim_Km",20,0,0.3},"angle_pim_Km");
    auto h_angle_Km_pm = df1.Histo1D({"h_angle_Km_pm","h_angle_Km_pm",20,0,0.3},"angle_Km_pm");

    // now we create 3*15 canvas and draw each one

    TCanvas *c1 = new TCanvas("c1","c1",3000,4000);
    
    
    // Create a TPad to hold the title
    TPad *titlePad = new TPad("titlePad","titlePad",0,0.93,1,1);
    titlePad->SetFillColor(0);
    titlePad->SetBorderMode(0);
    titlePad->SetBorderSize(0);
    titlePad->SetGridx();
    titlePad->SetGridy();
    titlePad->Draw();

    titlePad->cd();
    // Create the TPaveText object and add it to the TPad
    TPaveText *pt = new TPaveText(0,0,1,1,"brndc");
    pt->AddText("Momentum angles between particles");
    pt->SetTextAlign(22);
    pt->SetFillStyle(0);
    pt->SetBorderSize(0);
    pt->SetTextSize(0.5);
    pt->Draw();

    c1->cd();
    // Create the TPad for the histograms and divide it
    TPad *histPad = new TPad("histPad","histPad",0,0,1,0.93);
    histPad->SetFillColor(0);
    histPad->SetBorderMode(0);
    histPad->SetBorderSize(0);
    histPad->SetGridx();
    histPad->SetGridy();
    histPad->Divide(3,5);
    histPad->Draw();

    lhcbStyle->SetTitleY(1.03); 
    
    // Draw the histograms
    histPad->cd(1);
    h_angle_pp_pm->SetTitle("#font[40]{p+ and p-}; Angle (radians); Entries");
    h_angle_pp_pm->Draw();
    histPad->cd(2);
    h_angle_pp_Km->SetTitle("#font[40]{p+ and K-}; Angle (radians); Entries");
    h_angle_pp_Km->Draw();
    histPad->cd(3);
    h_angle_pp_pim->SetTitle("#font[40]{p+ and pi-}; Angle (radians); Entries");
    h_angle_pp_pim->Draw();
    histPad->cd(4);
    h_angle_pp_pip->SetTitle("#font[40]{p+ and pi+}; Angle (radians); Entries");
    h_angle_pp_pip->Draw();
    histPad->cd(5);
    h_angle_pp_Kp->SetTitle("#font[40]{p+ and K+}; Angle (radians); Entries");
    h_angle_pp_Kp->Draw();
    histPad->cd(6);
    h_angle_Kp_pm->SetTitle("#font[40]{K+ and p-}; Angle (radians); Entries");
    h_angle_Kp_pm->Draw();
    histPad->cd(7);
    h_angle_Kp_Km->SetTitle("#font[40]{K+ and K-}; Angle (radians); Entries");
    h_angle_Kp_Km->Draw();
    histPad->cd(8);
    h_angle_Kp_pim->SetTitle("#font[40]{K+ and pi-}; Angle (radians); Entries");
    h_angle_Kp_pim->Draw();
    histPad->cd(9);
    h_angle_Kp_pip->SetTitle("#font[40]{K+ and pi+}; Angle (radians); Entries");
    h_angle_Kp_pip->Draw();
    histPad->cd(10);
    h_angle_pip_pm->SetTitle("#font[40]{pi+ and p-}; Angle (radians); Entries");
    h_angle_pip_pm->Draw();
    histPad->cd(11);
    h_angle_pip_Km->SetTitle("#font[40]{pi+ and K-}; Angle (radians); Entries");
    h_angle_pip_Km->Draw();
    histPad->cd(12);
    h_angle_pip_pim->SetTitle("#font[40]{pi+ and pi-}; Angle (radians); Entries");
    h_angle_pip_pim->Draw();
    histPad->cd(13);
    h_angle_pim_pm->SetTitle("#font[40]{pi- and p-}; Angle (radians); Entries");
    h_angle_pim_pm->Draw();
    histPad->cd(14);
    h_angle_pim_Km->SetTitle("#font[40]{pi- and K-}; Angle (radians); Entries");
    h_angle_pim_Km->Draw();
    histPad->cd(15);
    h_angle_Km_pm->SetTitle("#font[40]{K- and p-}; Angle (radians); Entries");
    h_angle_Km_pm->Draw();

    c1->SaveAs("angle_plots.png");

    lhcbStyle->SetTitleY(1.01); 

    // create a new TCanvas and draw histograms for the ones we cut on
    TCanvas *c2 = new TCanvas("c2","c2",1800,500);
    c2->Divide(5,1);
    c2->cd(1);
    h_angle_pp_pip->Draw();
    c2->cd(2);
    h_angle_pim_pm->Draw();
    c2->cd(3);
    h_angle_pp_Kp->Draw();
    c2->cd(4);
    h_angle_Km_pm->Draw();
    c2->cd(5);
    h_angle_Kp_pip->Draw();
    c2->SaveAs("angle_plots_cut.png");


};

