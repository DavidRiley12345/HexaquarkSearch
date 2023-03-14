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



void allplots(){

    ROOT::EnableImplicitMT(); 

    // read in the LcLcpiSS data

    ROOT::RDataFrame df("B2LcLcpiSS/DecayTree", {"/data/lhcb01/mwhitehead/LcLcpi_2018_MU.root","/data/lhcb01/mwhitehead/LcLcpi_2018_MD.root","/data/lhcb01/mwhitehead/LcLcpi_2017_MU.root","/data/lhcb01/mwhitehead/LcLcpi_2017_MD.root"});

    // do the usual cuts on the LcLcpiSS data

    auto SS_df = df.Filter("Lambdacp_K_ProbNNk > 0.6 && Lambdacm_K_ProbNNk > 0.6 && Lambdacp_p_ProbNNp > 0.7 && Lambdacm_p_ProbNNp > 0.7 && Lambdacp_pi_ProbNNpi > 0.4 && Lambdacm_pi_ProbNNpi > 0.4")
    .Filter("LcLcpiSS_ENDVERTEX_CHI2 < 9")
    .Define("Angle_pp_pip", "(180/3.14)*(acos((Lambdacp_pi_PX*Lambdacp_p_PX + Lambdacp_pi_PY*Lambdacp_p_PY + Lambdacp_pi_PZ*Lambdacp_p_PZ) / (sqrt(pow(Lambdacp_pi_PX,2)+pow(Lambdacp_pi_PY,2)+pow(Lambdacp_pi_PZ,2)) * sqrt(pow(Lambdacp_p_PX,2)+pow(Lambdacp_p_PY,2)+pow(Lambdacp_p_PZ,2)))))")
    .Define("Angle_pm_pim", "(180/3.14)*(acos((Lambdacm_pi_PX*Lambdacm_p_PX + Lambdacm_pi_PY*Lambdacm_p_PY + Lambdacm_pi_PZ*Lambdacm_p_PZ) / (sqrt(pow(Lambdacm_pi_PX,2)+pow(Lambdacm_pi_PY,2)+pow(Lambdacm_pi_PZ,2)) * sqrt(pow(Lambdacm_p_PX,2)+pow(Lambdacm_p_PY,2)+pow(Lambdacm_p_PZ,2)))))")
    .Define("Angle_kp_pp", "(180/3.14)*(acos((Lambdacm_K_PX*Lambdacp_p_PX + Lambdacm_K_PY*Lambdacp_p_PY + Lambdacm_K_PZ*Lambdacp_p_PZ) / (sqrt(pow(Lambdacm_K_PX,2)+pow(Lambdacm_K_PY,2)+pow(Lambdacm_K_PZ,2)) * sqrt(pow(Lambdacp_p_PX,2)+pow(Lambdacp_p_PY,2)+pow(Lambdacp_p_PZ,2)))))")
    .Define("Angle_km_pm", "(180/3.14)*(acos((Lambdacp_K_PX*Lambdacm_p_PX + Lambdacp_K_PY*Lambdacm_p_PY + Lambdacp_K_PZ*Lambdacm_p_PZ) / (sqrt(pow(Lambdacp_K_PX,2)+pow(Lambdacp_K_PY,2)+pow(Lambdacp_K_PZ,2)) * sqrt(pow(Lambdacm_p_PX,2)+pow(Lambdacm_p_PY,2)+pow(Lambdacm_p_PZ,2)))))")
    .Define("Angle_kp_pip", "(180/3.14)*(acos((Lambdacm_K_PX*Lambdacp_pi_PX + Lambdacm_K_PY*Lambdacp_pi_PY + Lambdacm_K_PZ*Lambdacp_pi_PZ) / (sqrt(pow(Lambdacm_K_PX,2)+pow(Lambdacm_K_PY,2)+pow(Lambdacm_K_PZ,2)) * sqrt(pow(Lambdacp_pi_PX,2)+pow(Lambdacp_pi_PY,2)+pow(Lambdacp_pi_PZ,2)))))")
    .Filter("Angle_pp_pip >0.05 && Angle_pm_pim >0.05 && Angle_kp_pp >0.05 && Angle_km_pm >0.05 && Angle_kp_pip >0.05")
    .Filter("pi_ProbNNpi > 0.5")
    .Filter("pi_ProbNNghost < 0.1 && Lambdacp_K_ProbNNghost < 0.1 && Lambdacp_p_ProbNNghost < 0.1 && Lambdacp_pi_ProbNNghost < 0.1 && Lambdacm_K_ProbNNghost < 0.1 && Lambdacm_p_ProbNNghost < 0.1 && Lambdacm_pi_ProbNNghost < 0.1 ")
    .Filter("Lambdacp_K_ORIVX_CHI2 < 10 && Lambdacp_p_ORIVX_CHI2 < 10 && Lambdacp_pi_ORIVX_CHI2 < 10 && Lambdacm_K_ORIVX_CHI2 < 10 && Lambdacm_p_ORIVX_CHI2 < 10 && Lambdacm_pi_ORIVX_CHI2 < 10 ");
    
    ROOT::RDataFrame sdf("B2LcLcpiSSS/DecayTree", {"/data/lhcb01/mwhitehead/LcLcpi_2018_MU.root","/data/lhcb01/mwhitehead/LcLcpi_2018_MD.root","/data/lhcb01/mwhitehead/LcLcpi_2017_MU.root","/data/lhcb01/mwhitehead/LcLcpi_2017_MD.root"});

    // do the usual cuts on the LcLcpiSS data

    auto SSS_df = sdf.Filter("Lambdacp_K_ProbNNk > 0.6 && Lambdacm_K_ProbNNk > 0.6 && Lambdacp_p_ProbNNp > 0.7 && Lambdacm_p_ProbNNp > 0.7 && Lambdacp_pi_ProbNNpi > 0.4 && Lambdacm_pi_ProbNNpi > 0.4")
    .Filter("LcLcpiSSS_ENDVERTEX_CHI2 < 9")
    .Define("Angle_pp_pip", "(180/3.14)*(acos((Lambdacp_pi_PX*Lambdacp_p_PX + Lambdacp_pi_PY*Lambdacp_p_PY + Lambdacp_pi_PZ*Lambdacp_p_PZ) / (sqrt(pow(Lambdacp_pi_PX,2)+pow(Lambdacp_pi_PY,2)+pow(Lambdacp_pi_PZ,2)) * sqrt(pow(Lambdacp_p_PX,2)+pow(Lambdacp_p_PY,2)+pow(Lambdacp_p_PZ,2)))))")
    .Define("Angle_pm_pim", "(180/3.14)*(acos((Lambdacm_pi_PX*Lambdacm_p_PX + Lambdacm_pi_PY*Lambdacm_p_PY + Lambdacm_pi_PZ*Lambdacm_p_PZ) / (sqrt(pow(Lambdacm_pi_PX,2)+pow(Lambdacm_pi_PY,2)+pow(Lambdacm_pi_PZ,2)) * sqrt(pow(Lambdacm_p_PX,2)+pow(Lambdacm_p_PY,2)+pow(Lambdacm_p_PZ,2)))))")
    .Define("Angle_kp_pp", "(180/3.14)*(acos((Lambdacm_K_PX*Lambdacp_p_PX + Lambdacm_K_PY*Lambdacp_p_PY + Lambdacm_K_PZ*Lambdacp_p_PZ) / (sqrt(pow(Lambdacm_K_PX,2)+pow(Lambdacm_K_PY,2)+pow(Lambdacm_K_PZ,2)) * sqrt(pow(Lambdacp_p_PX,2)+pow(Lambdacp_p_PY,2)+pow(Lambdacp_p_PZ,2)))))")
    .Define("Angle_km_pm", "(180/3.14)*(acos((Lambdacp_K_PX*Lambdacm_p_PX + Lambdacp_K_PY*Lambdacm_p_PY + Lambdacp_K_PZ*Lambdacm_p_PZ) / (sqrt(pow(Lambdacp_K_PX,2)+pow(Lambdacp_K_PY,2)+pow(Lambdacp_K_PZ,2)) * sqrt(pow(Lambdacm_p_PX,2)+pow(Lambdacm_p_PY,2)+pow(Lambdacm_p_PZ,2)))))")
    .Define("Angle_kp_pip", "(180/3.14)*(acos((Lambdacm_K_PX*Lambdacp_pi_PX + Lambdacm_K_PY*Lambdacp_pi_PY + Lambdacm_K_PZ*Lambdacp_pi_PZ) / (sqrt(pow(Lambdacm_K_PX,2)+pow(Lambdacm_K_PY,2)+pow(Lambdacm_K_PZ,2)) * sqrt(pow(Lambdacp_pi_PX,2)+pow(Lambdacp_pi_PY,2)+pow(Lambdacp_pi_PZ,2)))))")
    .Filter("Angle_pp_pip >0.05 && Angle_pm_pim >0.05 && Angle_kp_pp >0.05 && Angle_km_pm >0.05 && Angle_kp_pip >0.05")
    .Filter("pi_ProbNNpi > 0.5")
    .Filter("pi_ProbNNghost < 0.1 && Lambdacp_K_ProbNNghost < 0.1 && Lambdacp_p_ProbNNghost < 0.1 && Lambdacp_pi_ProbNNghost < 0.1 && Lambdacm_K_ProbNNghost < 0.1 && Lambdacm_p_ProbNNghost < 0.1 && Lambdacm_pi_ProbNNghost < 0.1 ")
    .Filter("Lambdacp_K_ORIVX_CHI2 < 10 && Lambdacp_p_ORIVX_CHI2 < 10 && Lambdacp_pi_ORIVX_CHI2 < 10 && Lambdacm_K_ORIVX_CHI2 < 10 && Lambdacm_p_ORIVX_CHI2 < 10 && Lambdacm_pi_ORIVX_CHI2 < 10 ");

    // make a histogram of of the lcp invariant mass from each of the two decay modes

    TH1D h11 = *SSS_df.Histo1D({"Lambdacp_M_SS","Lambdacp_M_SS",40,2220,2350},"Lambdacp_M");
    TH1D h12 = *SS_df.Histo1D({"Lambdacp_M_SS","Lambdacp_M_SS",40,2220,2350},"Lambdacp_M");

    TCanvas *c1 = new TCanvas("c1","c1",800,600);
    h11.Add(&h11,&h12);
    h11.Draw("pe");
    c1->SaveAs("Lambdacp_M.png");

    // now create the lcpi invariant mass

    auto df2 = SS_df.Define("lcpip_M","sqrt(pow(Lambdacp_PE+pi_PE,2)-pow(Lambdacp_PX+pi_PX,2)-pow(Lambdacp_PY+pi_PY,2)-pow(Lambdacp_PZ+pi_PZ,2))")
                   .Define("lcpim_M","sqrt(pow(Lambdacm_PE+pi_PE,2)-pow(Lambdacm_PX+pi_PX,2)-pow(Lambdacm_PY+pi_PY,2)-pow(Lambdacm_PZ+pi_PZ,2))")
                   .Filter("Lambdacp_M > 2270 && Lambdacp_M < 2305")
                   .Filter("Lambdacm_M > 2270 && Lambdacm_M < 2305");
    TH1D h2 = *df2.Histo1D({"lcpi_M_SS","lcpi_M_SS",50,2400,3000},"lcpip_M");
    TH1D h2m= *df2.Histo1D({"lcpi_M_SS","lcpi_M_SS",50,2400,3000},"lcpim_M");

    auto sdf2 = SSS_df.Define("lcpip_M","sqrt(pow(Lambdacp_PE+pi_PE,2)-pow(Lambdacp_PX+pi_PX,2)-pow(Lambdacp_PY+pi_PY,2)-pow(Lambdacp_PZ+pi_PZ,2))")
                   .Define("lcpim_M","sqrt(pow(Lambdacm_PE+pi_PE,2)-pow(Lambdacm_PX+pi_PX,2)-pow(Lambdacm_PY+pi_PY,2)-pow(Lambdacm_PZ+pi_PZ,2))")
                   .Filter("Lambdacp_M > 2270 && Lambdacp_M < 2305")
                   .Filter("Lambdacm_M > 2270 && Lambdacm_M < 2305");
    TH1D sh2 = *sdf2.Histo1D({"lcpi_M_SS","lcpi_M_SS",50,2400,3000},"lcpip_M");
    TH1D sh2m= *sdf2.Histo1D({"lcpi_M_SS","lcpi_M_SS",50,2400,3000},"lcpim_M");

    // plot that histogram

    TCanvas *c2 = new TCanvas("c2","c2",800,600);
    h2.Add(&h2,&h2m);
    h2.Add(&h2,&sh2);
    h2.Add(&h2,&sh2m);
    h2.Draw("pe");
    c2->SaveAs("lcpi_M.png");

    // now finall plot the LcLcpiSS invariant mass after adding the mass window to the other Lcm

    auto df4 = df2.Filter("Lambdacm_M > 2270 && Lambdacm_M < 2305");
    auto h3 = *df4.Histo1D({"LcLcpiSS_M","LcLcpiSS_M",35,4500,20000},"LcLcpiSS_LcLcpiSS_DTF_M_Lambdacp_Lambdam_PV");

    auto sdf4 = sdf2.Filter("Lambdacm_M > 2270 && Lambdacm_M < 2305");
    auto sh3 = *sdf4.Histo1D({"LcLcpiSS_M","LcLcpiSS_M",35,4500,20000},"LcLcpiSSS_LcLcpiSSS_DTF_M_Lambdacp_Lambdam_PV");
    
    // plot that histogram

    TCanvas *c3 = new TCanvas("c3","c3",800,600);
    h3.Add(&h3,&sh3);
    h3.Draw("pe");
    c3->SaveAs("LcLcpiSS_M.png");

};