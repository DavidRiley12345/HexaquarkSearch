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
#include "TKDE.h"
#include "RooAbsPdf.h"
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooConstVar.h"

// use this order for safety on library loading
using namespace RooFit;
using namespace RooStats;

// see below for implementation

void bg_fit(){

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

    ROOT::EnableImplicitMT();

    auto df = ROOT::RDataFrame("B2LcLcpiOS/DecayTree",{"/data/lhcb01/mwhitehead/LcLcpi_2018_MU.root","/data/lhcb01/mwhitehead/LcLcpi_2018_MD.root","/data/lhcb01/mwhitehead/LcLcpi_2017_MU.root","/data/lhcb01/mwhitehead/LcLcpi_2017_MD.root"});

    auto df1 = df.Filter("Lambdacp_p_ProbNNp > 0.7 && Lambdacp_K_ProbNNk > 0.6 && Lambdacp_pi_ProbNNpi > 0.4 && Lambdacm_p_ProbNNp > 0.7 && Lambdacm_K_ProbNNk > 0.6 && Lambdacm_pi_ProbNNpi > 0.4 && pi_ProbNNpi > 0.4");
    auto df2 = df1.Define("Sigmac_M","sqrt(pow(Lambdacp_PE + pi_PE,2) - pow(Lambdacp_PX + pi_PX,2) - pow(Lambdacp_PY + pi_PY,2) - pow(Lambdacp_PZ + pi_PZ,2))");
        
    df2.Snapshot("DecayTree","/home/ppe/d/driley/git_HexaquarkSummer/LcLcpiOS/sweight_bg_fit/cut.root",{"Lambdacp_M","Lambdacm_M","Sigmac_M","LcLcpiOS_M"});

    // we want to load that data back into a roodataset

    TFile* fin = TFile::Open("/home/ppe/d/driley/git_HexaquarkSummer/LcLcpiOS/sweight_bg_fit/cut.root");
    TTree* tree = (TTree*)fin->Get("DecayTree");

    std::cout << "\033[1;31m"; 
    std::cout << "Tree Imported" << std::endl;
    std::cout << "\033[0m";

    // Create the RooDataSet from the tree

    RooRealVar Lambdacp_M("Lambdacp_M","",2220.0,2350.0);
    RooRealVar Lambdacm_M("Lambdacm_M","",2220.0,2350.0);
    RooRealVar Sigmac_M("Sigmac_M","",2350.0,3000.0);
    RooRealVar LcLcpiOS_M("LcLcpiOS_M","",2400.0,7000.0);

    RooDataSet *data = new RooDataSet("datain","",tree,RooArgSet(Lambdacp_M,Sigmac_M,Lambdacm_M,LcLcpiOS_M));

    // we now want to create the sweights on the lambdac mass

    // output message confirming data loaded to screen

    std::cout << "Data loaded" << std::endl;

    // Create the RooRealVar for the mass variable
    
    RooRealVar nsig_lc("nsig_lc","nsig_lc",10000,0,10000000);
    RooRealVar nbkg_lc("nbkg_lc","nbkg_lc",10000,0,10000000);

    // create signal pdfs as a double gaussian
    RooRealVar sig1_lc_mean("sig1_lc_mean","sig1_lc_mean",2287,2220,2350);
    RooRealVar sig1_lc_width("sig1_lc_width","sig1_lc_width",5,0.1,15);
    RooRealVar sig2_lc_width("sig2_lc_width","sig2_lc_width",5,0.1,15);
    RooGaussian sig1_lc("sig1_lc","sig1_lc",Lambdacp_M,sig1_lc_mean,sig1_lc_width);
    RooGaussian sig2_lc("sig2_lc","sig2_lc",Lambdacp_M,sig1_lc_mean,sig2_lc_width);
    RooRealVar sig1_lc_frac("sig1_lc_frac","sig1_lc_frac",0.5,0,1);

    // combine the signal pdf's

    RooAddPdf sig_lc("sig_lc","sig_lc",RooArgList(sig1_lc,sig2_lc),RooArgList(sig1_lc_frac));

    //create exponential background pdf
    RooRealVar bkg1_slope("bkg1_slope","bkg1_slope",-0.000238,-1,1);
    RooExponential bkg_lc("bkg_lc","bkg_lc",Lambdacp_M,bkg1_slope);

    // combine the signal and background pdf's
    RooAddPdf model_lc("model_lc","model_lc",RooArgList(sig_lc,bkg_lc),RooArgList(nsig_lc,nbkg_lc));

    // fit the model to the data

    model_lc.fitTo(*data,Extended());

    std::cout << "\033[1;31m"; 
    std::cout << "Lc model fitted" << std::endl;
    std::cout << "\033[0m";


    // lets plot the fit to see if it looks good

    TCanvas* c = new TCanvas("c","c",800,600);

    TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
    pad1->Draw();

    pad1->cd();
    RooPlot* frame = Lambdacp_M.frame();
    data->plotOn(frame);
    model_lc.plotOn(frame,LineColor(kBlack));
    model_lc.plotOn(frame,Components(bkg_lc),LineStyle(kDashed),LineColor(kRed));
    model_lc.plotOn(frame,Components(sig1_lc),LineStyle(kDashed),LineColor(kGreen));
    model_lc.plotOn(frame,Components(sig2_lc),LineStyle(kDashed),LineColor(kBlue));
    frame->SetTitle("#Lambda_{c} Mass Fit");
    frame->SetYTitle("Events");
    frame->SetXTitle("Mass (MeV)");
    frame->Draw();
    


    // lets sweight the data

    RooStats::SPlot* sData1 = new RooStats::SPlot("sData1","An SPlot",*data,&model_lc,RooArgList(nsig_lc,nbkg_lc));

    std::cout << "\033[1;31m"; 
    std::cout << "SPlot on Lc done" << std::endl;
    std::cout << "\033[0m";

    sData1->Print("v");

    // now lets plot the sweights themselves from the columns names sw_sig_lc and sw_bkg_lc

    TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
    c->cd();
    pad2->Draw();

    pad2->cd();

    c->SaveAs("lc_fit_sb.png");

    std::cout << "\033[1;31m"; 
    std::cout << "lc_fit_sb created" << std::endl;
    std::cout << "\033[0m";


    // now that we have the parameters we need for the model, we will fix these paramters allowing only the signal and background values to vary

    sig1_lc_mean.setConstant(); 
    sig1_lc_width.setConstant();
    sig2_lc_width.setConstant();
    sig1_lc_frac.setConstant();
    bkg1_slope.setConstant();
    // nsig_lc.setConstant(0);
    // nbkg_lc.setConstant(0);

    // now we want to apply the full mass window cut on the data

    // RooDataSet* data = (RooDataSet*)datain->reduce("Lambdacp_M > 2220 && Lambdacp_M < 2350");

    // // now lets fit the model to the data again with the mass window applied

    // we need to create a whole new model with the new parameters

    

    Lambdacp_M.setRange("signal",2220,2350);

    RooGaussian sig1_lc2("sig1_lc2","sig1_lc2",Lambdacp_M,sig1_lc_mean,sig1_lc_width);
    RooGaussian sig2_lc2("sig2_lc2","sig2_lc2",Lambdacp_M,sig1_lc_mean,sig2_lc_width);
    RooAddPdf sig_lc2("sig_lc2","sig_lc2",RooArgList(sig1_lc2,sig2_lc2),RooArgList(sig1_lc_frac));
    RooExponential bkg_lc2("bkg_lc2","bkg_lc2",Lambdacp_M,bkg1_slope);

    RooRealVar nsig_lc2("nsig_lc2","nsig_lc2",10000,0,10000000);
    RooRealVar nbkg_lc2("nbkg_lc2","nbkg_lc2",10000,0,10000000);

    RooAddPdf model_lc2("model_lc2","model_lc2",RooArgList(sig_lc2,bkg_lc2),RooArgList(nsig_lc2,nbkg_lc2));

    model_lc2.fitTo(*data,Extended(),Range("signal"));

    std::cout << "\033[1;31m"; 
    std::cout << "no sb model fit" << std::endl;
    std::cout << "\033[0m";


    // lets plot this new fit to see if it looks resonable

    TCanvas* c2 = new TCanvas("c2","c2",800,600);
    c2->cd();
    RooPlot* frame2 = Lambdacp_M.frame();
    data->plotOn(frame2);
    model_lc2.plotOn(frame2,LineColor(kBlack));
    model_lc2.plotOn(frame2,Components(bkg_lc2),LineStyle(kDashed),LineColor(kRed));
    model_lc2.plotOn(frame2,Components(sig_lc2),LineStyle(kDashed),LineColor(kGreen));
    frame2->Draw();
    c2->SaveAs("lc_fit_nosb.png");

    std::cout << "\033[1;31m"; 
    std::cout << "no sb figure created" << std::endl;
    std::cout << "\033[0m";


    // Lambdacp_M.setMin(2270);
    // Lambdacp_M.setMax(2305);

    // now we want to create the splot

    RooMsgService::instance().setSilentMode(true);

    SPlot* sData = new SPlot("sData","An SPlot",*data,&model_lc2,RooArgList(nsig_lc2,nbkg_lc2));

    std::cout << "\033[1;31m"; 
    std::cout << "no sb sweights created" << std::endl;
    std::cout << "\033[0m";

    
    // data should now contain the sweights lets output the variables in data to check

    // Check that our weights have the desired properties
   std::cout << "Check SWeights:" << std::endl;
   std::cout << std::endl <<  "Yield of signal is "
               << nsig_lc.getVal() << ".  From sWeights it is "
               << sData->GetYieldFromSWeight("nsig_lc2_sw") << std::endl;
   std::cout << "Yield of bkg is "
               << nbkg_lc.getVal() << ".  From sWeights it is "
               << sData->GetYieldFromSWeight("nbkg_lc2_sw") << std::endl
               << std::endl;
   for(Int_t i=0; i < 10; i++)
      {
      std::cout << "sig   " << sData->GetSWeight(i,"nsig_lc2")
                  << "   bkg Weight   " << sData->GetSWeight(i,"nbkg_lc2")
                  << "  Total Weight   " << sData->GetSumOfEventSWeight(i)
                  << std::endl;
      }
   std::cout << std::endl;

    sData->Print();
    data->Print();
    
    // lets get some information 

    std::cout << "sum of sweights" << sData->GetSumOfEventSWeight(1) << std::endl;
    std::cout << "sum of sweights" << sData->GetSumOfEventSWeight(1) << std::endl;
    
    // create RooRealVar for sweights

    RooRealVar nsig_lc2_sw("nsig_lc2_sw","nsig_lc2_sw",-10,10);
    RooRealVar nbkg_lc2_sw("nbkg_lc2_sw","nbkg_lc2_sw",-10,10);
    
    // now lets plot the sweights

    TCanvas* c11 = new TCanvas("c11","c11",800,800);
    c11->cd();
    // add tpad

    lhcbStyle->SetTitleH(0.1);

    TPad* lc_fit_pad= new TPad("lc_fit_pad","lc_fit_pad",0,0.3,1,1);
    lc_fit_pad->SetTopMargin(0.1);
    lc_fit_pad->SetBottomMargin(0.18);
    lc_fit_pad->Draw();
    lc_fit_pad->cd();
    frame2->SetTitle("#Lambda_{c} mass with sWeights");
    frame2->GetXaxis()->SetTitle("m(#Lambda_{c}) [MeV]");
    frame2->Draw();

    // TLegend 

    TLegend* leg = new TLegend(0.7,0.65,0.9,0.85);
    leg->AddEntry(frame2->getObject(0),"Data","lep");
    leg->AddEntry(frame2->getObject(1),"Total Fit","l");
    leg->AddEntry(frame2->getObject(2),"Bkg Fit","l");
    leg->AddEntry(frame2->getObject(3),"Sig Fit","l");
    leg->Draw();
    
    c11->cd();

    TPad* sWeight_pad = new TPad("sWeight_pad","sWeight_pad",0,0,1,0.3);
    sWeight_pad->SetBottomMargin(0.2);
    sWeight_pad->SetTopMargin(0.01);
    sWeight_pad->Draw();
    sWeight_pad->cd();


    RooPlot* frame11 = Lambdacp_M.frame();
    data->plotOnXY(frame11,YVar(nsig_lc2_sw),MarkerColor(kGreen));
    data->plotOnXY(frame11,YVar(nbkg_lc2_sw),MarkerColor(kRed));
    frame11->SetTitle(0);
    frame11->GetYaxis()->SetTitleSize(0.14);
    frame11->GetYaxis()->SetTitleOffset(0.6);
    frame11->SetYTitle("SWeights");
    frame11->SetXTitle("");
    frame11->SetLabelSize(0,"X");
    frame11->SetLabelSize(0.1,"Y");
    frame11->Draw();

    TLegend* leg2 = new TLegend(0.7,0.8,0.95,0.95);
    leg2->AddEntry(frame11->getObject(0),"Sig","lep");
    leg2->AddEntry(frame11->getObject(1),"Bkg","lep");
    leg2->Draw();

    c11->SaveAs("lc_sweights.png");


    // c11->SaveAs("lc_sweights.png");

    // lets fill histograms with Sigmac_M wieghted by the sweights

    TH1F* h_sig = new TH1F("h_sig","h_sig",25,2350,3000);
    TH1F* h_bkg = new TH1F("h_bkg","h_bkg",25,2350,3000);

    // fill histograms using fillHistogram and weighting the Sigmac_M mass with each sweight

    // data->fillHistogram(h_sig,RooArgList(Sigmac_M),nullptr,"nsig_lc_sw");
    // data->fillHistogram(h_bkg,RooArgList(Sigmac_M),nullptr,"nbkg_lc_sw");
    
    for (int i=0; i<data->numEntries(); i++) {
        // Get values and weights for this entry
        double sigmac_m = data->get(i)->getRealValue("Sigmac_M");
        double weight1 = data->get(i)->getRealValue("nsig_lc2_sw");
        double weight2 = data->get(i)->getRealValue("nbkg_lc2_sw");

        // Fill histograms with weighted entries
        h_sig->Fill(sigmac_m, weight1);
        h_bkg->Fill(sigmac_m, weight2);
    }



    // plot those histograms on the same canvas, with the different hist in a separate canvas

    TCanvas* c10 = new TCanvas("c10","c10",800,600);
    c10->cd();
    h_sig->SetMarkerColor(kGreen);
    h_sig->SetMarkerStyle(20);
    h_sig->SetTitle("#Sigma_{c} invariant mass");
    h_sig->SetYTitle("Events");
    h_sig->SetXTitle("m(pK#pi#pi) [MeV/c^{2}]");
    h_bkg->SetMarkerColor(kRed);
    h_bkg->SetMarkerStyle(20); 
    h_sig->Draw("PE");
    h_bkg->Draw("PE SAME");

    TLegend *legend = new TLegend(0.7,0.7,0.9,0.9);
    legend->AddEntry(h_sig,"Signal");
    legend->AddEntry(h_bkg,"Background");
    legend->Draw();

    c10->SaveAs("sweighted_sigmac.png");
       
    std::cout << "\033[1;31m"; 
    std::cout << "sweighted sig_c figure created" << std::endl;
    std::cout << "\033[0m";


    // we have bkg pdf stored in h_bkg lets fit this with the RooDstD0BG

    // create a datahist from the histogram 
    RooDataSet *dataSet = new RooDataSet("dataSet", "Data Set", RooArgSet(Sigmac_M));
    for (int i = 0; i < h_bkg->GetNbinsX(); ++i) {
        Sigmac_M.setVal(h_bkg->GetBinCenter(i + 1));
        dataSet->add(RooArgSet(Sigmac_M), h_bkg->GetBinContent(i + 1));
    }
    // create the RooKeysPdf

    RooKeysPdf bkgPDF("bkgPDF","bkgPDF",Sigmac_M,*dataSet,RooKeysPdf::MirrorBoth,1.5);

    // // create the model using DstD0BG

    // RooRealVar dm0("dm0","dm0", 2350,2350,2500);
    // RooRealVar b1("b1","b1",600, 200, 1000);
    // RooRealVar c1("c1","c1",-0.001,-0.5,0.5);
    // RooRealVar d1("d1","d1",-0.001,-0.5,0.5);
    // RooDstD0BG bkgPDF("bkgPDF","bkgPDF",Sigmac_M,dm0,b1,c1,d1);

    // bkgPDF.fitTo(h_bkg_datahist);

    // // lets set those as constant

    // dm0.setConstant();
    // b1.setConstant();
    // c1.setConstant();
    // d1.setConstant();

    // lets plot the bkgPDF and the datahist

    // TCanvas* c12 = new TCanvas("c12","c12",800,600);
    // RooPlot* frame123 = Sigmac_M.frame();
    // bkgPDF.plotOn(frame123);
    // frame123->Draw();
    // c12->SaveAs("bkgpdf.png");

    // RooKeysPdf bkgPDF("bkgPDF","bkgPDF",Sigmac_M,h_bkg_datahist,RooKeysPdf::MirrorBoth,1.5);

    // now we can now try and fit the Sigmac_M with this new bkgPDF

    // create a new model

    RooRealVar nsig("nsig","nsig",0,100000);
    RooRealVar nbkg("nbkg","nbkg",0,10000000);

    //create two breitwigner functions

    RooRealVar mean1("mean1","mean1",2455,2400,2500);
    RooRealVar width1("width1","width1",10,0,30);
    RooBreitWigner bw1("bw1","bw1",Sigmac_M,mean1,width1);

    RooRealVar mean2("mean2","mean2",2520,2500,2550);
    RooRealVar width2("width2","width2",10,0,30);
    RooGaussian bw2("bw2","bw2",Sigmac_M,mean2,width2);

    // add a ratio for these two breitwigners

    RooRealVar ratio("ratio","ratio",0.5,-1,1);

    // sum these two pdfs

    RooAddPdf sigc_total("sigc_total","sigc_total",RooArgList(bw1,bw2),RooArgList(ratio));

    // create a total pdf

    RooAddPdf model("model","model",RooArgList(sigc_total,bkgPDF),RooArgList(nsig,nbkg));

    // fit the model to the data

    model.fitTo(*data);

    data->Print();

    // plot the model

    TCanvas* c12 = new TCanvas("c12","c12",800,600);
    c12->cd();
    RooPlot* frame12 = Sigmac_M.frame();
    data->plotOn(frame12);
    model.plotOn(frame12);
    model.plotOn(frame12,Components(bkgPDF),LineStyle(kDashed));
    model.plotOn(frame12,Components(sigc_total),LineStyle(kDashed));
    frame12->Draw();
    c12->SaveAs("sigmac_fit.png");

    // now do sweights on the Sigmac_M fit

    // RooStats::SPlot* sData2 = new RooStats::SPlot("sData2","An SPlot",*data,&model,RooArgList(nsig,nbkg));

    // // create new histograms for the sweighted LcLcpiOS_M

    // TH1F* h_sig2 = new TH1F("h_sig2","h_sig2",30,4500,8000);
    // TH1F* h_bkg2 = new TH1F("h_bkg2","h_bkg2",30,4500,8000);

    // // fill histograms using fillHistogram and weighting the Sigmac_M mass with each sweight

    // for (int i=0; i<data->numEntries(); i++) {
    //     // Get values and weights for this entry
    //     double lclcpi_m = data->get(i)->getRealValue("LcLcpiOS_M");
    //     double weight1 = data->get(i)->getRealValue("nsig_sw");
    //     double weight2 = data->get(i)->getRealValue("nbkg_sw");

    //     // Fill histograms with weighted entries

    //     h_sig2->Fill(lclcpi_m, weight1);
    //     h_bkg2->Fill(lclcpi_m, weight2);
    // }

    // // we have to add a cut on Lambdacm_M to plot only the signal region

    // data->reduce("Lambdacm_M > 2270 && Lambdacm_M < 2305");

    // // plot those histograms

    // TCanvas* c13 = new TCanvas("c13","c13",800,1200);
    // c13->cd(1);
    // h_sig2->SetFillColor(kGreen);
    // h_bkg2->SetFillColor(kRed);
    // h_sig2->Draw();
    // h_bkg2->Draw("SAME");
    // c13->SaveAs("sweighted_lclcpi.png");


}
    