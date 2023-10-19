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

void Optimize2DCut() {

    // get the current TStyle

    TStyle *currentStyle = gStyle;

    // set the TStyle statbox to zero

    currentStyle->SetOptStat(0);
    
    TFile* fin = TFile::Open("/data/lhcb01/mwhitehead/LcLcpi_2018_MD.root");
    ROOT::RDataFrame d("B2LcLcpiOS/DecayTree", fin);

    const double lowerLimit = 2220;
    const double upperLimit = 2350;
    const int nbins = 50;

    TH2D* hSignificance = new TH2D("hSignificance", "Significance;Proton cut value;Kaon cut value", 10, 0, 1, 10, 0, 1);
    TH2D* hPurity = new TH2D("hPurity", "Purity;Proton cut value ;Kaon cut value", 10, 0, 1, 10, 0, 1);

    TCanvas* cFits = new TCanvas("cFits", "Fit plots", 2400, 2400);
    cFits->Divide(10, 10);

    // Loop over proton and kaon cut values
    int padIndex = 1;
    for (double p_cut = 0.05; p_cut <= 0.95; p_cut += 0.1) {
        for (double k_cut = 0.05; k_cut <= 0.95; k_cut += 0.1) {
            auto d_cut = d.Filter(TString::Format("(Lambdacp_p_ProbNNp > %f) && (Lambdacp_K_ProbNNk > %f)", p_cut, k_cut).Data());
            auto h = d_cut.Histo1D({"h", "", nbins, lowerLimit, upperLimit}, "Lambdacp_M");

            RooRealVar mass("mass", "", lowerLimit, upperLimit);

            TH1D* th1d = static_cast<TH1D*>(h.GetValue().Clone());

            RooDataHist *data = new RooDataHist("h2", "h2", RooArgList(mass), th1d);

            // Signal and background models
            RooRealVar mean("mean", "", 2287., 2282.0, 2290.0);
            RooRealVar sigma("sigma", "", 6, 5, 7);
            RooGaussian signal("signal", "", mass, mean, sigma);
            RooRealVar slope("slope", "", -0.01, 0.01);
            RooExponential bkgPDF("comboPDF", "", mass, slope);
            RooRealVar nsig("nsig", "", 40000, 100, 10000000);
            RooRealVar nbkg("nbkg", "", 20000, 100, 10000000);
            RooAddPdf totalPDF("totalPDF", "", RooArgList(signal, bkgPDF), RooArgList(nsig, nbkg));

            // Fit and results
            RooFitResult* fitresult = totalPDF.fitTo(*data, RooFit::Extended());
            double sig = nsig.getValV();
            double bkg = nbkg.getValV();
            double purity = sig / (sig + bkg);
            double significance = sig / TMath::Sqrt(sig + bkg);

            // Output
            std::cout << "Proton cut: " << p_cut << ", Kaon cut: " << k_cut << std::endl;
            std::cout << "Signal : " << sig << std::endl;
            std::cout << "Backgr : " << bkg << std::endl;
            std::cout << "Purity : " << purity << std::endl;
            std::cout << "Signif : " << significance << std::endl;

            // Plotting
            // TCanvas* c3 = new TCanvas("c3", TString::Format("Fit for p cut at %f, k cut at %f", p_cut, k_cut), 800, 600);
            // c3->SetLeftMargin(0.12);
            // c3->SetRightMargin(0.12);

            // RooPlot* frame = mass.frame(RooFit::Title(""));
            // data->plotOn(frame);
            // totalPDF.plotOn(frame);
            // totalPDF.plotOn(frame, RooFit::Components("comboPDF"), RooFit::LineStyle(kDashed));
            // frame->Draw();

            // c3->SaveAs(TString::Format("Fit_p_cut_%f_k_cut_%f.png", p_cut, k_cut));


            cFits->cd(padIndex);
            RooPlot* frame1 = mass.frame(RooFit::Title(""));
            // data->plotOn(frame1);
            // totalPDF.plotOn(frame1);
            // totalPDF.plotOn(frame1, RooFit::Components("comboPDF"), RooFit::LineStyle(kDashed));
            // frame1->Draw();

            padIndex++;

            hSignificance->Fill(p_cut, k_cut, significance);
            hPurity->Fill(p_cut, k_cut, purity);

            // Clear memory for the next iteration
            // delete c3;
            // delete frame;
            data->reset();
            totalPDF.getComponents()->removeAll();
        }
    }

// Save histograms of significance and purity
TCanvas* cSignificance = new TCanvas("cSignificance", "Significance", 800, 800);
cSignificance->SetLeftMargin(0.12);
cSignificance->SetRightMargin(0.12);
hSignificance->Draw("COLZ");
hSignificance->SetTitleSize(0.05, "XYZ");
cSignificance->SaveAs("Significance_2D_plot.png");

TCanvas* cPurity = new TCanvas("cPurity", "Purity", 800, 800);
cPurity->SetLeftMargin(0.12);
cPurity->SetRightMargin(0.12);
hPurity->Draw("COLZ");
hPurity->SetTitleSize(0.05, "XYZ");
cPurity->SaveAs("Purity_2D_plot.png");

cFits->SaveAs("All_Fit_Plots.png");

// Clean up
delete cSignificance;
delete cPurity;
delete hSignificance;
delete hPurity;
delete fin;
}

int main() {
    Optimize2DCut();
    return 0;
}
