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


void lcpi_lc(){

    ROOT::EnableImplicitMT();

    std::vector<std::string> names = {"/data/lhcb01/mwhitehead/LcLcpi_2018_MU.root","/data/lhcb01/mwhitehead/LcLcpi_2018_MD.root","/data/lhcb01/mwhitehead/LcLcpi_2017_MU.root","/data/lhcb01/mwhitehead/LcLcpi_2017_MD.root"};

    ROOT::RDataFrame df("B2LcLcpiOS/DecayTree", names);

    //define the lambdac plus pion invariant mass

    auto df1 = df.Define("lcpi_M","sqrt(pow(Lambdacp_PE + pi_PE,2) - pow(Lambdacp_PX + pi_PX,2) - pow(Lambdacp_PY + pi_PY,2) - pow(Lambdacp_PZ + pi_PZ,2))");

    auto df2 = df1.Filter("Lambdacp_p_ProbNNp > 0.7 && Lambdacm_p_ProbNNp > 0.7 && Lambdacp_K_ProbNNk > 0.6 && Lambdacm_K_ProbNNk > 0.6 && Lambdacp_pi_ProbNNpi > 0.4 && Lambdacm_pi_ProbNNpi > 0.4 && (Lambdacp_M > 2220 && Lambdacp_M < 2350)");

    //create a model for our 2d histogram

    auto model =  ROOT::RDF::TH2DModel("model", "model", 30, 2350, 5000, 30, 2220, 2350);
    
    auto hist = df2.Histo2D(model, "lcpi_M", "Lambdacp_M");
    double corr = hist->GetCorrelationFactor();

    TCanvas c1;
    hist->Draw("colz");
    hist->SetStats(0);
    hist->SetTitle("LambdacM and Lcpi Correlation Plot");
    hist->GetXaxis()->SetTitle("Lc Mass (MeV)");
    hist->GetYaxis()->SetTitle("Lcpi Mass (MeV)");

    // Add text box with correlation factor

    TPaveText *textbox = new TPaveText(0.57, 0.78, 0.90, 0.90, "NDC");
    textbox->SetFillColor(0);
    textbox->SetTextFont(43);
    textbox->SetTextSize(20);
    textbox->SetBorderSize(1);
    textbox->SetLineColor(1);
    textbox->SetLineWidth(2);
    textbox->AddText(Form("Correlation Factor: %.3f", corr));
    textbox->Draw();
            
    
    c1.SaveAs("lcpi_lc.png");

    std::cout << "Correlation Factor: " << corr << std::endl;

};