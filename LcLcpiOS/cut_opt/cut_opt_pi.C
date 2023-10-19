#include <TFile.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TF1.h>
#include <TGraph.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TMath.h>
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooGaussian.h>
#include <RooExponential.h>
#include <RooAddPdf.h>
#include <RooArgList.h>
#include <RooFitResult.h>
#include <iostream>
#include <array>

void fit(const double cutval, std::array<double, 19>& purities, std::array<double, 19>& sigs, std::array<double, 19>& total)
{
    TFile* fin = TFile::Open("/data/lhcb01/mwhitehead/LcLcpi_2018_MD.root");
    TTree* tin = dynamic_cast<TTree*>(fin->Get("B2LcLcpiOS/DecayTree"));

    const double lowerLimit = 2220;
    const double upperLimit = 2350;
    const int nbins = 50;
    const double binwidth = (upperLimit - lowerLimit) / nbins;
    TH1D* h = new TH1D("h", "", nbins, lowerLimit, upperLimit);
    tin->Draw("Lambdacp_M>>h",  TString::Format("(Lambdacp_ProbNNp > %f) && (Lambdacp_ProbNNk > %f)", p_cut, k_cut));

    RooRealVar mass("mass", "", lowerLimit, upperLimit);
    RooDataHist data("data", "data", mass, h);
    RooRealVar mean("mean", "", 2287., 2270.0, 2290.0);
    RooRealVar sigma("sigma", "", 4., 0.0, 20.0);
    RooGaussian signal("signal", "", mass, mean, sigma);
    RooRealVar slope("slope", "", -1000, 1000);
    RooExponential bkgPDF("comboPDF", "", mass, slope);
    RooRealVar nsig("nsig", "", -1000, 1000000);
    RooRealVar nbkg("nbkg", "", -1000, 10000000);
    RooAddPdf totalPDF("totalPDF", "", RooArgList(signal, bkgPDF), RooArgList(nsig, nbkg));
    RooFitResult* fitresult = totalPDF.fitTo(data, RooFit::Extended());

    double sig = nsig.getValV();
    double bkg = nbkg.getValV();
    double purity = sig / (sig+bkg);
    double significance = sig / TMath::Sqrt(sig + bkg);

    purities.push_back(purity);
    sigs.push_back(significance);
    total.push_back(sig+bkg);

    std::cout << "cut val: " << cutval << std::endl;
    std::cout << "Signal : " << sig << std::endl;
    std::cout << "Backgr : " << bkg << std::endl;
    std::cout << "Purity : " << purity << std::endl;
    std::cout << "Signif : " << significance << std::endl;

    TCanvas* c3 = new TCanvas("c3", TString::Format("Fit for p cut at %f", cutval), 800, 600);
    c3->SetLeftMargin(0.12);
    c3->SetRightMargin(0.12);
}