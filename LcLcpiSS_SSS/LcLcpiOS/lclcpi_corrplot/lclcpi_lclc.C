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



void lclcpi_lclc(){

    ROOT::EnableImplicitMT();

    std::vector<std::string> names = {"/data/lhcb01/mwhitehead/LcLcpi_2018_MU.root","/data/lhcb01/mwhitehead/LcLcpi_2018_MD.root","/data/lhcb01/mwhitehead/LcLcpi_2017_MU.root","/data/lhcb01/mwhitehead/LcLcpi_2017_MD.root"};

    ROOT::RDataFrame df("B2LcLcpiSS/DecayTree", names);
    
    auto df1 = df.Define("lclc","sqrt(pow(Lambdacp_PE + Lambdacm_PE,2) - pow(Lambdacp_PX + Lambdacm_PX,2) - pow(Lambdacp_PY + Lambdacm_PY,2) - pow(Lambdacp_PZ + Lambdacm_PZ,2))");

    auto df2 = df1.Filter("Lambdacp_p_ProbNNp > 0.7 && Lambdacm_p_ProbNNp > 0.7 && Lambdacp_K_ProbNNk > 0.6 && Lambdacm_K_ProbNNk > 0.6 && Lambdacp_pi_ProbNNpi > 0.4 && Lambdacm_pi_ProbNNpi > 0.4 && (Lambdacp_M > 2270 && Lambdacp_M < 2305) && (Lambdacm_M > 2270 && Lambdacm_M < 2305)");

    //create a model for our 2d histogram

    auto model =  ROOT::RDF::TH2DModel("model", "model", 20, 4600, 8000, 20, 4450, 6000);
    
    auto hist = df2.Histo2D(model, "LcLcpiSS_M", "lclc");
    double corr = hist->GetCorrelationFactor();

    TCanvas c1;
    hist->Draw("colz");
    hist->SetTitle("LcLcpi LcLc Correlation Plot");
    hist->GetXaxis()->SetTitle("LcLcpi Mass (MeV)");
    hist->GetYaxis()->SetTitle("LcLc Mass (MeV)");
    c1.SaveAs("lclcpi_lclc.png");

    std::cout << "Correlation Factor: " << corr << std::endl;

    //lets also output the lclc invariant mass distribution

    auto hist2 = df2.Histo1D({"LcLc invariant mass", "LcLc invariant mass", 50, 4450, 6000}, "lclc");
    TCanvas c2;
    hist2->Draw();
    hist2->SetTitle("LcLc Invariant Mass Distribution");
    hist2->GetXaxis()->SetTitle("LcLc Mass (MeV)");
    hist2->GetYaxis()->SetTitle(("Candidates / " + std::to_string(hist->GetBinWidth(1)) + " MeV").c_str());
    c2.SaveAs("lclc.png");

};