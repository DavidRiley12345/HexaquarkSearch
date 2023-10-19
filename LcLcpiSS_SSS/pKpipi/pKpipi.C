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

void pKpipi() {

    ROOT::EnableImplicitMT();

    ROOT::RDataFrame din("B2LcLcpiSS/DecayTree", {"/data/lhcb01/mwhitehead/LcLcpi_2018_MD.root", "/data/lhcb01/mwhitehead/LcLcpi_2018_MU.root","/data/lhcb01/mwhitehead/LcLcpi_2017_MD.root","/data/lhcb01/mwhitehead/LcLcpi_2017_MU.root"});

    auto d = din.Filter("Lambdacp_p_ProbNNp > 0.7 && Lambdacp_K_ProbNNk > 0.6 && Lambdacp_pi_ProbNNpi > 0.4 && pi_ProbNNpi > 0.5")
                .Filter("Lambdacp_M > 2270 && Lambdacp_M < 2305")
                .Define("Angle_pp_pip", "(180/3.14)*(acos((Lambdacp_pi_PX*Lambdacp_p_PX + Lambdacp_pi_PY*Lambdacp_p_PY + Lambdacp_pi_PZ*Lambdacp_p_PZ) / (sqrt(pow(Lambdacp_pi_PX,2)+pow(Lambdacp_pi_PY,2)+pow(Lambdacp_pi_PZ,2)) * sqrt(pow(Lambdacp_p_PX,2)+pow(Lambdacp_p_PY,2)+pow(Lambdacp_p_PZ,2)))))")
                .Define("Angle_pm_pim", "(180/3.14)*(acos((Lambdacm_pi_PX*Lambdacm_p_PX + Lambdacm_pi_PY*Lambdacm_p_PY + Lambdacm_pi_PZ*Lambdacm_p_PZ) / (sqrt(pow(Lambdacm_pi_PX,2)+pow(Lambdacm_pi_PY,2)+pow(Lambdacm_pi_PZ,2)) * sqrt(pow(Lambdacm_p_PX,2)+pow(Lambdacm_p_PY,2)+pow(Lambdacm_p_PZ,2)))))")
                .Define("Angle_kp_pp", "(180/3.14)*(acos((Lambdacm_K_PX*Lambdacp_p_PX + Lambdacm_K_PY*Lambdacp_p_PY + Lambdacm_K_PZ*Lambdacp_p_PZ) / (sqrt(pow(Lambdacm_K_PX,2)+pow(Lambdacm_K_PY,2)+pow(Lambdacm_K_PZ,2)) * sqrt(pow(Lambdacp_p_PX,2)+pow(Lambdacp_p_PY,2)+pow(Lambdacp_p_PZ,2)))))")
                .Define("Angle_km_pm", "(180/3.14)*(acos((Lambdacp_K_PX*Lambdacm_p_PX + Lambdacp_K_PY*Lambdacm_p_PY + Lambdacp_K_PZ*Lambdacm_p_PZ) / (sqrt(pow(Lambdacp_K_PX,2)+pow(Lambdacp_K_PY,2)+pow(Lambdacp_K_PZ,2)) * sqrt(pow(Lambdacm_p_PX,2)+pow(Lambdacm_p_PY,2)+pow(Lambdacm_p_PZ,2)))))")
                .Define("Angle_kp_pip", "(180/3.14)*(acos((Lambdacm_K_PX*Lambdacp_pi_PX + Lambdacm_K_PY*Lambdacp_pi_PY + Lambdacm_K_PZ*Lambdacp_pi_PZ) / (sqrt(pow(Lambdacm_K_PX,2)+pow(Lambdacm_K_PY,2)+pow(Lambdacm_K_PZ,2)) * sqrt(pow(Lambdacp_pi_PX,2)+pow(Lambdacp_pi_PY,2)+pow(Lambdacp_pi_PZ,2)))))")
                .Filter("Angle_pp_pip >0.05 && Angle_pm_pim >0.05 && Angle_kp_pp >0.05 && Angle_km_pm >0.05 && Angle_kp_pip >0.05")
                .Filter("Lambdacp_K_ORIVX_CHI2 < 10 && Lambdacp_p_ORIVX_CHI2 < 10 && Lambdacp_pi_ORIVX_CHI2 < 10")
                .Define("Sigmac_M", "sqrt(pow(Lambdacp_p_PE + pi_PE,2) - pow(Lambdacp_p_PX + pi_PX,2) - pow(Lambdacp_p_PY + pi_PY,2) - pow(Lambdacp_p_PZ + pi_PZ,2))");
    // we want to fit the invariant mass of the pKpi system using a double gaussian and exponential background

    // define the data

    auto h = d.Histo1D({"h", "h",40, 1000, 6000}, "Sigmac_M");

    TH1D* th1d = static_cast<TH1D*>(h.GetValue().Clone());

    RooRealVar m("m", "m", 1000, 6000);

    RooDataHist data("data", "data", RooArgList(m), th1d);

    // plot this data on a canvas

    RooPlot* frame = m.frame();
    
    data.plotOn(frame);

    TCanvas* c = new TCanvas("c", "c", 800, 600);

    frame->Draw();

    c->SaveAs("pKpipi.png");

}