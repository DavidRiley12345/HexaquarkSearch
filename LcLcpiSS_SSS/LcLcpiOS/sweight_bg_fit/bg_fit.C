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

void bg_fit(){
    ROOT::EnableImplicitMT();

    auto df = ROOT::RDataFrame("B2LcLcpiOS/DecayTree",{"/data/lhcb01/mwhitehead/LcLcpi_2018_MU.root","/data/lhcb01/mwhitehead/LcLcpi_2018_MD.root","/data/lhcb01/mwhitehead/LcLcpi_2017_MU.root","/data/lhcb01/mwhitehead/LcLcpi_2017_MD.root"});

    auto df1 = df.Filter("Lambdacp_p_ProbNNp > 0.7 && Lambdacp_K_ProbNNk > 0.6 && Lambdacp_pi_ProbNNpi > 0.4 && Lambdacm_p_ProbNNp > 0.7 && Lambdacm_K_ProbNNk > 0.6 && Lambdacm_pi_ProbNNpi > 0.4");
    auto df2 = df1.Define("Sigmac_M","sqrt(pow(Lambdacp_PE + pi_PE,2) - pow(Lambdacp_PX + pi_PX,2) - pow(Lambdacp_PY + pi_PY,2) - pow(Lambdacp_PZ + pi_PZ,2))");
        
    df2.Snapshot("DecayTree","/home/ppe/d/driley/git_HexaquarkSummer/LcLcpiOS/sweight_bg_fit/cut.root",{"Lambdacp_M","Lambdacm_M","Sigmac_M","LcLcpiOS_M"});

    // we want to load that data back into a roodataset

    TFile* fin = TFile::Open("/home/ppe/d/driley/git_HexaquarkSummer/LcLcpiOS/sweight_bg_fit/cut.root");
    TTree* tree = (TTree*)fin->Get("DecayTree");

    // Create the RooDataSet from the tree

    RooRealVar Lambdacp_M("Lambdacp_M","",2220.0,2350.0);
    RooRealVar Lambdacm_M("Lambdacm_M","",2220.0,2350.0);
    RooRealVar Sigmac_M("Sigmac_M","",2400.0,3000.0);
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
    RooGaussian sig_lc("sig_lc","sig_lc",Lambdacp_M,sig1_lc_mean,sig1_lc_width);

    //create exponential background pdf
    RooRealVar bkg1_slope("bkg1_slope","bkg1_slope",-0.000238);
    RooExponential bkg_lc("bkg_lc","bkg_lc",Lambdacp_M,bkg1_slope);

    // combine the signal and background pdf's
    RooAddPdf model_lc("model_lc","model_lc",RooArgList(sig_lc,bkg_lc),RooArgList(nsig_lc,nbkg_lc));

    // fit the model to the data

    model_lc.fitTo(*data,Extended());

    // lets plot the fit to see if it looks good

    TCanvas* c = new TCanvas("c","c",800,600);
    c->cd();
    RooPlot* frame = Lambdacp_M.frame();
    data->plotOn(frame);
    model_lc.plotOn(frame);
    model_lc.plotOn(frame,Components(bkg_lc),LineStyle(kDashed));
    model_lc.plotOn(frame,Components(sig_lc),LineStyle(kDashed),LineColor(kRed));
    frame->Draw();
    c->SaveAs("lc_fit_sb.png");
    
    // now that we have the parameters we need for the model, we will fix these paramters allowing only the signal and background values to vary

    // sig1_lc_mean.setConstant(); 
    // sig1_lc_width.setConstant();
    // bkg1_slope.setConstant();
    // nsig_lc.setConstant(0);
    // nbkg_lc.setConstant(0);

    // now we want to apply the full mass window cut on the data

    // RooDataSet* data = (RooDataSet*)datain->reduce("Lambdacp_M > 2220 && Lambdacp_M < 2350");

    // // now lets fit the model to the data again with the mass window applied

    // we need to create a whole new model with the new parameters

    

    Lambdacp_M.setRange("signal",2270,2305);

    RooGaussian sig_lc2("sig_lc2","sig_lc2",Lambdacp_M,sig1_lc_mean,sig1_lc_width);
    RooExponential bkg_lc2("bkg_lc2","bkg_lc2",Lambdacp_M,bkg1_slope);

    RooRealVar nsig_lc2("nsig_lc2","nsig_lc2",10000,0,10000000);
    RooRealVar nbkg_lc2("nbkg_lc2","nbkg_lc2",10000,0,10000000);

    RooAddPdf model_lc2("model_lc2","model_lc2",RooArgList(sig_lc2,bkg_lc2),RooArgList(nsig_lc2,nbkg_lc2));

    model_lc2.fitTo(*data,Extended(),Range("signal"));

    // lets plot this new fit to see if it looks resonable

    TCanvas* c2 = new TCanvas("c2","c2",800,600);
    c2->cd();
    RooPlot* frame2 = Lambdacp_M.frame();
    data->plotOn(frame2);
    model_lc2.plotOn(frame2);
    model_lc2.plotOn(frame2,Components(bkg_lc2),LineStyle(kDashed));
    model_lc2.plotOn(frame2,Components(sig_lc2),LineStyle(kDashed),LineColor(kRed));
    frame2->Draw();
    c2->SaveAs("lc_fit_nosb.png");

    Lambdacp_M.setMin(2270);
    Lambdacp_M.setMax(2305);

    // now we want to create the splot

    SPlot* sData = new SPlot("sData","An SPlot",*data,&model_lc2,RooArgList(nsig_lc2,nbkg_lc2));
    
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

    // TCanvas* c11 = new TCanvas("c11","c11",800,600);
    // c11->cd();
    // RooPlot* frame11 = Lambdacp_M.frame();
    // data->plotOnXY(frame11,YVar(nsig_lc_sw),MarkerColor(kRed));
    // data->plotOnXY(frame11,YVar(nbkg_lc_sw),MarkerColor(kBlue));
    // frame11->Draw();
    // c11->SaveAs("lc_sweights.png");

    // lets fill histograms with Sigmac_M wieghted by the sweights

    TH1F* h_sig = new TH1F("h_sig","h_sig",30,2400,3000);
    TH1F* h_bkg = new TH1F("h_bkg","h_bkg",30,2400,3000);

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

    // create a histogram which shows the difference between these two

    TH1F* h_diff = new TH1F("h_diff","h_diff",40,2400,3000);

    for (int i = 1; i < data->numEntries(); i++) {
        h_diff->SetBinContent(i,h_sig->GetBinContent(i)-h_bkg->GetBinContent(i));
    }

    // plot those histograms on the same canvas, with the different hist in a separate canvas

    TCanvas* c10 = new TCanvas("c10","c10",800,1200);
    c10->Divide(1,2);
    c10->cd(1);
    h_sig->SetFillColor(kGreen);
    h_bkg->SetFillColor(kRed);
    h_sig->Draw();
    h_bkg->Draw("SAME");

    c10->cd(2);
    h_diff->SetMarkerColor(kBlue);
    h_diff->Draw("pe");
    c10->SaveAs("sweighted_sigmac.png");
       
    // we have the bkgPDF stored in h_bkg, lets use this to create a roodatahist

    RooDataHist h_bkg_datahist("bkgPDF","bkgPDF",RooArgList(Sigmac_M),h_bkg);

    RooHistPdf bkgPDF("bkgPDF","bkgPDF",RooArgSet(Sigmac_M),h_bkg_datahist,2);

    // lets use the data hist to create a rookeyspdf

    // RooKeysPdf bkgPDF("bkgPDF","bkgPDF",Sigmac_M,h_bkg_datahist,RooKeysPdf::MirrorBoth,1.5);

    // now we can now try and fit the Sigmac_M with this new bkgPDF

    // create a new model

    RooRealVar nsig("nsig","nsig",0,100000);
    RooRealVar nbkg("nbkg","nbkg",0,10000000);

    //create two breitwigner functions

    RooRealVar mean1("mean1","mean1",2455,2400,2500);
    RooRealVar width1("width1","width1",10,0,100);
    RooBreitWigner bw1("bw1","bw1",Sigmac_M,mean1,width1);

    RooRealVar mean2("mean2","mean2",2520,2500,2550);
    RooRealVar width2("width2","width2",10,0,100);
    RooBreitWigner bw2("bw2","bw2",Sigmac_M,mean2,width2);

    // add a ratio for these two breitwigners

    RooRealVar ratio("ratio","ratio",0.5,-10,10);

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

    RooStats::SPlot* sData2 = new RooStats::SPlot("sData2","An SPlot",*data,&model,RooArgList(nsig,nbkg));

    // create new histograms for the sweighted LcLcpiOS_M

    TH1F* h_sig2 = new TH1F("h_sig2","h_sig2",30,4500,8000);
    TH1F* h_bkg2 = new TH1F("h_bkg2","h_bkg2",30,4500,8000);

    // fill histograms using fillHistogram and weighting the Sigmac_M mass with each sweight

    for (int i=0; i<data->numEntries(); i++) {
        // Get values and weights for this entry
        double lclcpi_m = data->get(i)->getRealValue("LcLcpiOS_M");
        double weight1 = data->get(i)->getRealValue("nsig_sw");
        double weight2 = data->get(i)->getRealValue("nbkg_sw");

        // Fill histograms with weighted entries

        h_sig2->Fill(lclcpi_m, weight1);
        h_bkg2->Fill(lclcpi_m, weight2);
    }

    // we have to add a cut on Lambdacm_M to plot only the signal region

    data->reduce("Lambdacm_M > 2270 && Lambdacm_M < 2305");

    // plot those histograms

    TCanvas* c13 = new TCanvas("c13","c13",800,1200);
    c13->cd(1);
    h_sig2->SetFillColor(kGreen);
    h_bkg2->SetFillColor(kRed);
    h_sig2->Draw();
    h_bkg2->Draw("SAME");
    c13->SaveAs("sweighted_lclcpi.png");


}
    