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

void AddModel(RooWorkspace*);
void AddData(RooWorkspace*);
void DoSPlot(RooWorkspace*);
void MakePlots(RooWorkspace*);
void sigmac_gensweight()
{
   // Create a new workspace to manage the project.
   RooWorkspace* wspace = new RooWorkspace("myWS");
   std::cout << "Wkspace made" << std::endl; 
   // add the signal and background models to the workspace.
   // Inside this function you will find a description our model.
   AddModel(wspace);
   std::cout << "model added" << std::endl;
   // add some toy data to the workspace
   AddData(wspace);
   std::cout << "data added" << std::endl;
   // inspect the workspace if you wish
   //  wspace->Print();
   // do sPlot.
   //This wil make a new dataset with sWeights added for every event.
   DoSPlot(wspace);
   // Make some plots showing the discriminating variable and
   // the control variable after unfolding.
   MakePlots(wspace);
   // cleanup
   delete wspace;
}
//____________________________________
void AddModel(RooWorkspace* ws){

   RooRealVar Sigmac_M("Sigmac_M","",2390.0,2700.0);
   RooRealVar nsig("nsig","",10000,1,1000000);
   RooRealVar nbkg("nbkg","",100000,1,10000000);

   //create a two gaussian PDF, one for the each of the two sigmac masses we observed

   RooRealVar  gauss1Mean("gauss1Mean","",2453.75,2453.,2455.);
   RooRealVar  gauss1Std("gauss1Std","",0.1,15.0);
   RooGaussian gaussian1Sig("gaussian1Sig","",Sigmac_M,gauss1Mean,gauss1Std);

   RooRealVar gauss2Mean("gauss2Mean","",2518.4,2500.,2525);
   RooRealVar gauss2Std("gauss2Std","",0.1,30);
   RooGaussian gaussian2Sig("gaussian2Sig","",Sigmac_M,gauss2Mean,gauss2Std);

   //mixing parameters for the two gaussians
   
   RooRealVar mix1("mix1","",1.,10000000);
   RooRealVar mix2("mix2","",1.,10000000);

   //combine the two gaussians
   
   RooAddPdf signalPDF("signalPDF","",RooArgList(gaussian1Sig,gaussian2Sig),RooArgList(mix1,mix2));

   //add a background model, we found this function seemed to fit well, in future could be a good idea to use the 
   //sweighted background in the sigma plots to model background, but this shape will do here
   RooRealVar dm0("dm0","",2420); 
   RooRealVar a("a","",100000,0,10000000);
   RooRealVar b("b","",-0.005,-4,0.);
   RooRealVar c("c","",-0.005,-4,0.);
   RooDstD0BG bkgPDF("bkgPDF","",Sigmac_M,dm0,a,b,c);

   RooAddPdf model("model","",RooArgList(signalPDF,bkgPDF),RooArgList(nsig,nbkg));	 
   RooRealVar pi_PT("pi_PT","",0,25000);
   ws->import(pi_PT);
   ws->import(model);
}
//____________________________________
void AddData(RooWorkspace* ws){
   // Add a toy dataset
   // how many events do we want?
   //TFile* fin = TFile :: Open("/home/ppe/d/driley/LcLcpiOS/Sweight/cut.root");
   //TTree* tin = (TTree*)fin->Get("Lcm_cut");
   //TFile* fin = TFile :: Open("/data/lhcb01/mwhitehead/LcLcpi_2018_MD.root");
   //TTree* tin = (TTree*)fin->Get("B2LcLcpiOS/DecayTree");
   std::cout << "adding sigmac_M to tree:" << std::endl;
   ROOT::RDataFrame df1("B2LcLcpiOS/DecayTree","/data/lhcb01/mwhitehead/LcLcpi_2018_MD.root");
   std::cout << "opening tree as RDataFrame to define and filter" << std::endl;
   auto df2 = df1.Define("Sigmac_M","sqrt(pow(pi_PE + Lambdacp_PE,2) - pow(pi_PX + Lambdacp_PX,2) - pow(pi_PY + Lambdacp_PY,2) - pow(pi_PZ + Lambdacp_PZ,2))");
   std::cout << "defined: filtering ..." << std::endl;
   auto df3 = df2.Filter("Lambdacp_K_ProbNNk > 0.4 && Lambdacm_K_ProbNNk > 0.4 && Lambdacp_p_ProbNNp > 0.4 && Lambdacm_p_ProbNNp > 0.4 && Lambdacp_pi_ProbNNpi > 0.4 && Lambdacm_pi_ProbNNpi > 0.4 && (Lambdacp_M < 2350 && Lambdacp_M > 2220)");
   df3.Snapshot("sigmac","/home/ppe/d/driley/git_HexaquarkSummer/Sweight/sigmac_sw/sigmac.root",{"Sigmac_M","pi_PT"});
   std::cout << "snapshot saved" << std::endl;
   TFile* fin = TFile :: Open("/home/ppe/d/driley/git_HexaquarkSummer/Sweight/sigmac_sw/sigmac.root");
   TTree* tin = (TTree*)fin->Get("sigmac");
   RooRealVar* Sigmac_M = ws->var("Sigmac_M");   
   RooRealVar* pi_PT = ws->var("pi_PT"); 
   std::cout << ws->allVars() << std::endl;
   RooDataSet *datain = new RooDataSet("","",tin,RooArgSet(*Sigmac_M,*pi_PT));
   datain->Print("v");
   std::cout <<"---------------------------------" << std::endl;
 
   // import data into workspace
   ws->import(*datain, Rename("data"));
}
//____________________________________
void DoSPlot(RooWorkspace* ws){
   std::cout << "Calculate sWeights" << std::endl;
   // get what we need out of the workspace to do the fit
   RooAbsPdf* model = ws->pdf("model");
   RooRealVar* nsig = ws->var("nsig");
   RooRealVar* nbkg = ws->var("nbkg");
   RooDataSet* data = (RooDataSet*) ws->data("data");
   // fit the model to the data.
   model->fitTo(*data, Extended() );
   // The sPlot technique requires that we fix the parameters
   // of the model that are not yields after doing the fit.
   //
   // This *could* be done with the lines below, however this is taken care of
   // by the RooStats::SPlot constructor (or more precisely the AddSWeight
   // method).
   //
   //RooRealVar* sigmaZ = ws->var("sigmaZ");
   //RooRealVar* qcdMassDecayConst = ws->var("qcdMassDecayConst");
   //sigmaZ->setConstant();
   //qcdMassDecayConst->setConstant();
   RooMsgService::instance().setSilentMode(true);
   // Now we use the SPlot class to add SWeights to our data set
   // based on our model and our yield variables
   RooStats::SPlot* sData = new RooStats::SPlot("sData","An SPlot",
                                             *data, model, RooArgList(*nsig,*nbkg) );
   // Check that our weights have the desired properties
   std::cout << "Check SWeights:" << std::endl;
   std::cout << std::endl <<  "Yield of signal is "
               << nsig->getVal() << ".  From sWeights it is "
               << sData->GetYieldFromSWeight("nsig") << std::endl;
   std::cout << "Yield of bkg is "
               << nbkg->getVal() << ".  From sWeights it is "
               << sData->GetYieldFromSWeight("nbkg") << std::endl
               << std::endl;
   for(Int_t i=0; i < 10; i++)
      {
      std::cout << "sig   " << sData->GetSWeight(i,"nsig")
                  << "   bkg Weight   " << sData->GetSWeight(i,"nbkg")
                  << "  Total Weight   " << sData->GetSumOfEventSWeight(i)
                  << std::endl;
      }
   std::cout << std::endl;
   // import this new dataset with sWeights
   
   std::cout << "------------------------------------import new dataset with sWeights" << std::endl;
   data->Print("v");
   ws->import(*data, Rename("dataWithSWeights"));
}
void MakePlots(RooWorkspace* ws){
   // Here we make plots of the discriminating variable (invMass) after the fit
   // and of the control variable (isolation) after unfolding with sPlot.
   std::cout << "make plots" << std::endl;
   // make our canvas
   TCanvas* cdata = new TCanvas("sPlot","sPlot demo", 800, 800);
   cdata->Divide(2,2);
   // get what we need out of the workspace
   RooAbsPdf* model = ws->pdf("model");
   RooAbsPdf* signalPDF = ws->pdf("signalPDF");
   RooAbsPdf* bkgPDF = ws->pdf("bkgPDF");
   RooRealVar* nsig_sw = ws->var("nsig_sw");
   RooRealVar* nbkg_sw = ws->var("nbkg_sw");
   RooRealVar* pi_PT=ws->var("pi_PT");
   RooRealVar* Sigmac_M = ws->var("Sigmac_M");
   // note, we get the dataset with sWeights
   RooDataSet* data = (RooDataSet*) ws->data("dataWithSWeights");
   // this shouldn't be necessary, need to fix something with workspace
   // do this to set parameters back to their fitted values.
   model->fitTo(*data, Extended() );
   //plot invMass for data with full model and individual components overlaid
   //  TCanvas* cdata = new TCanvas();
   cdata->cd(1);
   RooPlot* frame = Sigmac_M->frame() ;
   data->plotOn(frame) ;
   model->plotOn(frame) ;
   model->plotOn(frame,Components(*signalPDF),LineStyle(kDashed), LineColor(kGreen)) ;
   model->plotOn(frame,Components(*bkgPDF),LineStyle(kDashed),LineColor(kRed)) ;
   frame->SetTitle("Sigmac_M (Gaussian signal, exp bkg");
   frame->Draw();
   
   cdata->cd(3);
   data->get()->Print("v");
   RooPlot* frame1 = Sigmac_M->frame();
   data->plotOnXY(frame1, YVar(*nsig_sw),MarkerColor(kGreen));
   data->plotOnXY(frame1, YVar(*nbkg_sw),MarkerColor(kRed));
   frame1->SetTitle("Plot of signal and bkg sweights");
   frame1->Draw();
   
   cdata->cd(2);
   RooPlot* frame2 = pi_PT->frame();
   RooDataSet wdata_sigsw("","",data,*data->get(),nullptr,"nsig_sw");
   wdata_sigsw.plotOn(frame2);
   frame2->SetTitle("Transverse mom: signal weight");
   frame2->Draw();
 
   cdata->cd(4);
   RooPlot* frame3 = pi_PT->frame();
   RooDataSet wdata_bkgsw("","",data,*data->get(),nullptr,"nbkg_sw");
   wdata_bkgsw.plotOn(frame3);
   frame3->SetTitle("Transverse mom: bkg weight");
   frame3->Draw();   
   cdata->SaveAs("SPlot.gif");
}

