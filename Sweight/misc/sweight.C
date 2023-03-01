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
void sweight()
{
   ROOT::EnableImplicitMT(); 
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

   RooRealVar Lambdacp_M("Lambdacp_M","",2220.0,2350.0);
   RooRealVar nsig("nsig","",10000,-1000,1000000);
   RooRealVar nbkg("nbkg","",100000,-1000,10000000);

   RooRealVar gaussMean("gaussMean","",2250.0,2300.0);
   RooRealVar gaussSig("gaussSig","",0.1,15.0);
   RooGaussian gaussianSig("gaussianSig","",Lambdacp_M,gaussMean,gaussSig);

   RooRealVar expConst("expConst","",-0.002,-10,10);
   RooExponential bkgPDF("bkgPDF","",Lambdacp_M,expConst);

   RooAddPdf model("model","",RooArgList(gaussianSig,bkgPDF),RooArgList(nsig,nbkg));	 
   RooRealVar Lambdacp_PT("Lambdacp_PT","",0,25000);
   ws->import(Lambdacp_PT);
   ws->import(model);
}
//____________________________________
void AddData(RooWorkspace* ws){
   // Add a toy dataset
   // how many events do we want?
   TFile* fin = TFile :: Open("/home/ppe/d/driley/git_HexaquarkSummer/Sweight/cut.root");
   TTree* tin = (TTree*)fin->Get("Lcm");
   //TFile* fin = TFile :: Open("/data/lhcb01/mwhitehead/LcLcpi_2018_MD.root");
   //TTree* tin = (TTree*)fin->Get("B2LcLcpiOS/DecayTree");
   RooRealVar* Lambdacp_M = ws->var("Lambdacp_M");   
   RooRealVar* Lambdacp_PT = ws->var("Lambdacp_PT"); 
   std::cout << ws->allVars() << std::endl;
   RooDataSet *datain = new RooDataSet("","",tin,RooArgSet(*Lambdacp_M,*Lambdacp_PT));
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
   cdata->Divide(1,3);
   // get what we need out of the workspace
   RooAbsPdf* model = ws->pdf("model");
   RooAbsPdf* gaussianSig = ws->pdf("gaussianSig");
   RooAbsPdf* bkgPDF = ws->pdf("bkgPDF");
   RooRealVar* nsig_sw = ws->var("nsig_sw");
   RooRealVar* nbkg_sw = ws->var("nbkg_sw");
   RooRealVar* Lambdacp_PT=ws->var("Lambdacp_PT");
   RooRealVar* Lambdacp_M = ws->var("Lambdacp_M");
   // note, we get the dataset with sWeights
   RooDataSet* data = (RooDataSet*) ws->data("dataWithSWeights");
   // this shouldn't be necessary, need to fix something with workspace
   // do this to set parameters back to their fitted values.
   model->fitTo(*data, Extended() );
   //plot invMass for data with full model and individual components overlaid
   //  TCanvas* cdata = new TCanvas();
   cdata->cd(1);
   RooPlot* frame = Lambdacp_M->frame() ;
   data->plotOn(frame) ;
   model->plotOn(frame) ;
   model->plotOn(frame,Components(*gaussianSig),LineStyle(kDashed), LineColor(kGreen)) ;
   model->plotOn(frame,Components(*bkgPDF),LineStyle(kDashed),LineColor(kRed)) ;
   frame->SetTitle("Lambdacp_M (Gaussian signal, exp bkg");
   frame->Draw();
   
   cdata->cd(2);
   data->get()->Print("v");
   RooPlot* frame1 = Lambdacp_M->frame();
   data->plotOnXY(frame1, YVar(*nsig_sw),MarkerColor(kGreen));
   data->plotOnXY(frame1, YVar(*nbkg_sw),MarkerColor(kRed));
   frame1->SetTitle("Plot of signal and bkg sweights");
   frame1->Draw();
   
   cdata->cd(3);
   RooPlot* frame2 = Lambdacp_PT->frame();
   RooDataSet wdata_bkgsw("","",data,*data->get(),nullptr,"nbkg_sw");
   RooDataSet wdata_sigsw("","",data,*data->get(),nullptr,"nsig_sw");
   wdata_bkgsw.plotOn(frame2,MarkerColor(kRed));
   wdata_sigsw.plotOn(frame2,MarkerColor(kGreen));
   frame2->SetTitle("Lcm_PT sig vs bkg");
   frame2->Draw();
   cdata->SaveAs("SPlot.gif");
}

