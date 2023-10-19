#include <ROOT/RDataFrame.hxx>
#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooGaussian.h>
#include <RooExponential.h>
#include <RooAddPdf.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <TMath.h>

void pKpi_SG() {

    TFile *fin = TFile::Open("/data/lhcb01/mwhitehead/LcLcpi_2018_MD.root");
    ROOT::RDataFrame din("B2LcLcpiOS/DecayTree", fin);

    auto d = din.Filter("Lambdacp_p_ProbNNp > 0.7 && Lambdacp_K_ProbNNk > 0.6 && Lambdacp_pi_ProbNNpi > 0.4");

    // we want to fit the invariant mass of the pKpi system using a double gaussian and exponential background

    // define the signal variables

    RooRealVar m("m", "m", 2220, 2350);
    RooRealVar mean1("mean1", "mean1", 2287, 2270, 2290);
    RooRealVar sigma1("sigma1", "sigma1", 0, 20);
    RooRealVar frac("frac", "frac", -10, 10);
    RooRealVar nsig("nsig", "nsig", 0, 1000000);

    // define the background variables

    RooRealVar a("a", "a", -0.1, 0.1);
    RooRealVar nbkg("nbkg", "nbkg", 1, 10000000);

    // define the pdfs

    RooGaussian sg("sg", "sg", m, mean1, sigma1);
    RooExponential bkg("bkg", "bkg", m, a);

    // define the total pdf

    RooAddPdf model("model", "model", RooArgList(sg, bkg), RooArgList(nsig, nbkg));

    // define the data

    auto h = d.Histo1D({"h", "h", 40, 2220, 2350}, "Lambdacp_M");

    TH1D* th1d = static_cast<TH1D*>(h.GetValue().Clone());

    RooDataHist data("data", "data", RooArgList(m), th1d);

    // fit the data

    RooFitResult *r = model.fitTo(data, RooFit::Save());

    // create a pullplot of the fit calculating the difference between the data and the fit divided by the error

    // plot the data and the fit

    RooPlot *frame = m.frame();
    frame->SetTitle("pKpi Gaussian Fit; m{pK#pi} (MeV/c^{2}); Entries");
    data.plotOn(frame); 
    model.plotOn(frame, RooFit::Components("bkg"), RooFit::LineStyle(kDashed), RooFit::LineColor(kRed));
    model.plotOn(frame, RooFit::Components("sg"), RooFit::LineStyle(kDashed), RooFit::LineColor(kGreen));
    model.plotOn(frame, RooFit::LineColor(kBlack));


    TCanvas *c = new TCanvas("c", "c", 800, 800);
    c->cd();
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1);
    pad1->SetBottomMargin(0.1);
    pad1->SetLeftMargin(0.15);
    pad1->Draw();
    pad1->cd();
    frame->Draw();

    // add a legend with the chi2 and the number of entries

    TLegend *leg = new TLegend(0.6, 0.8, 0.9, 0.9);
    leg->AddEntry(frame->getObject(3), Form("#chi^{2}: %.2f", frame->chiSquare()), "l");
    leg->Draw();
    
    // plot the pull plot below the main plot
    c->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.3);
    pad2->SetTopMargin(0.1);
    pad2->SetBottomMargin(0.4);
    pad2->SetLeftMargin(0.15);
    pad2->Draw();
    pad2->cd();

    RooHist *pullHist = frame->pullHist();
    RooPlot *pullPlot = m.frame();
    pullPlot->SetTitle("");
    pullPlot->addPlotable(pullHist, "P");
    pullPlot->SetYTitle("Pull");
    pullPlot->SetTitleSize(0.1, "XYZ");
    pullPlot->SetLabelSize(0.1, "XYZ");
    pullPlot->SetTitleOffset(0.5, "Y");
    pullPlot->SetNdivisions(505, "Y");
    pullPlot->Draw();

    c->SaveAs("pKpifit_sg.png");

    // print fit results

    r->Print();

    // calc the chi2

    double chi2 = frame->chiSquare();
    std::cout << "Chi2: " << chi2 << std::endl;
    
}