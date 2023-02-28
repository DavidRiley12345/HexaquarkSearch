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
void gensweight()
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
   
   RooRealVar Lambdacm_M("Lambdacm_M","",2220.0,2350.0);
   RooRealVar Lamndacm_P("Lambdacm_P","Lambdacm_P",0,450000);
   RooRealVar Lambdacm_PT("Lambdacm_PT","Lambdacm_PT",0,30000); 
   RooRealVar Lambdacp_P("Lambdacp_P","Lambdacp_P",0,450000);
   RooRealVar Lambdacp_PT("Lambdacp_PT","Lambdacp_PT",0,30000);
   RooRealVar Lambdacp_p_P("Lambdacp_p_P","Lambdacp_p_P",0,220000);
   RooRealVar Lambdacp_p_PT("Lambdacp_p_PT","Lambdacp_p_PT",0,16000);
   RooRealVar Lambdacp_K_P("Lambdacp_K_P","Lambdacp_K_P",0,125000);
   RooRealVar Lambdacp_K_PT("Lambdacp_K_PT","Lambdacp_K_PT",0,13000);
   RooRealVar Lambdacp_pi_P("Lambdacp_pi_P","Lambdacp_pi_P",0,140000);
   RooRealVar Lambdacp_pi_PT("Lambdacp_pi_PT","Lambdacp_pi_PT",0,140000);
   RooRealVar pi_P("pi_P","pi_P",0,900000);
   RooRealVar pi_PT("pi_PT","pi_PT",0,16000);
   RooRealVar pi_PE("pi_PE","pi_PE",0,220000);
   RooRealVar pi_PX("pi_PX","pi_PX",-10000,10000);
   RooRealVar pi_PY("pi_PY","pi_PY",-10000,10000);
   RooRealVar pi_PZ("pi_PZ","pi_PZ",0,220000);
   RooRealVar Lambdacm_PE("Lambdacm_PE","Lambdacm_PE",0,450000);
   RooRealVar Lambdacm_PX("Lambdacm_PX","Lambdacm_PX",-30000,30000);
   RooRealVar Lambdacm_PY("Lambdacm_PY","Lambdacm_PY",-30000,30000);
   RooRealVar Lambdacm_PZ("Lambdacm_PZ","Lambdacm_PZ",0,450000); 
   RooRealVar Lambdacp_PE("Lambdacp_PE","Lambdacp_PE",0,450000);
   RooRealVar Lambdacp_PX("Lambdacp_PX","Lambdacp_PX",-30000,30000);
   RooRealVar Lambdacp_PY("Lambdacp_PY","Lambdacp_PY",-30000,30000);
   RooRealVar Lambdacp_PZ("Lambdacp_PZ","Lambdacp_PZ",0,450000);
   RooRealVar totCandidates("totCandidates","totCandidates",0,7000);  
   ws->import(Lambdacp_P);
   ws->import(Lambdacp_PT);
   ws->import(Lambdacp_p_P);
   ws->import(Lambdacp_p_PT);
   ws->import(Lambdacp_K_P);
   ws->import(Lambdacp_K_PT);
   ws->import(Lambdacp_pi_P);
   ws->import(Lambdacp_pi_PT);
   ws->import(pi_P);
   ws->import(pi_PT);
   ws->import(pi_PE);
   ws->import(pi_PX);
   ws->import(pi_PY);
   ws->import(pi_PZ);
   ws->import(Lambdacp_PE);
   ws->import(Lambdacp_PX);
   ws->import(Lambdacp_PY);
   ws->import(Lambdacp_PZ);
   ws->import(Lambdacm_PE);
   ws->import(Lambdacm_PX);
   ws->import(Lambdacm_PY);
   ws->import(Lambdacm_PZ);
   ws->import(totCandidates);
   ws->import(model);
}
//____________________________________
void AddData(RooWorkspace* ws){
   // Add a toy dataset
   // how many events do we want?
   TFile* fin = TFile :: Open("/home/ppe/d/driley/git_HexaquarkSummer/Sweight/cut.root");
   TTree* tin = (TTree*)fin->Get("Lcm");
   std::cout << "tfile opened" << std::endl;
   //TFile* fin = TFile :: Open("/data/lhcb01/mwhitehead/LcLcpi_2018_MD.root");
   //TTree* tin = (TTree*)fin->Get("B2LcLcpiOS/DecayTree");
   RooRealVar* Lambdacp_M = ws->var("Lambdacp_M");   
   RooRealVar* Lambdacp_P = ws->var("Lambdacp_P");
   RooRealVar* Lambdacp_PT = ws->var("Lambdacp_PT");
   RooRealVar* Lambdacp_p_P = ws->var("Lambdacp_p_P");
   RooRealVar* Lambdacp_p_PT = ws->var("Lambdacp_p_PT");
   RooRealVar* Lambdacp_K_P = ws->var("Lambdacp_K_P");
   RooRealVar* Lambdacp_K_PT = ws->var("Lambdacp_K_PT");
   RooRealVar* Lambdacp_pi_P = ws->var("Lambdacp_pi_P");
   RooRealVar* Lambdacp_pi_PT = ws->var("Lambdacp_pi_PT");
   RooRealVar* pi_P = ws->var("pi_P");
   RooRealVar* pi_PT = ws->var("pi_PT");
   RooRealVar* pi_PE = ws->var("pi_PE");
   RooRealVar* pi_PX = ws->var("pi_PX");
   RooRealVar* pi_PY = ws->var("pi_PY");
   RooRealVar* pi_PZ = ws->var("pi_PZ");
   RooRealVar* Lambdacp_PE = ws->var("Lambdacp_PE");
   RooRealVar* Lambdacp_PX = ws->var("Lambdacp_PX");
   RooRealVar* Lambdacp_PY = ws->var("Lambdacp_PY");
   RooRealVar* Lambdacp_PZ = ws->var("Lambdacp_PZ"); 
   RooRealVar* Lambdacm_PE = ws->var("Lambdacm_PE");
   RooRealVar* Lambdacm_PX = ws->var("Lambdacm_PX");
   RooRealVar* Lambdacm_PY = ws->var("Lambdacm_PY");
   RooRealVar* Lambdacm_PZ = ws->var("Lambdacm_PZ");
   RooRealVar* totCandidates = ws->var("totCandidates");
   std::cout << ws->allVars() << std::endl;
   auto variables = new RooArgSet();
   variables->add(*Lambdacp_M);
   variables->add(*Lambdacp_P);
   variables->add(*Lambdacp_PT);
   variables->add(*Lambdacp_p_P);
   variables->add(*Lambdacp_p_PT);
   variables->add(*Lambdacp_K_P);
   variables->add(*Lambdacp_K_PT);
   variables->add(*Lambdacp_pi_P);
   variables->add(*Lambdacp_pi_PT);
   variables->add(*pi_P);
   variables->add(*pi_PT);
   variables->add(*pi_PE);
   variables->add(*pi_PX);
   variables->add(*pi_PY);
   variables->add(*pi_PZ);
   variables->add(*Lambdacp_PE);
   variables->add(*Lambdacp_PX);
   variables->add(*Lambdacp_PY);
   variables->add(*Lambdacp_PZ);
   variables->add(*Lambdacm_PE);
   variables->add(*Lambdacm_PX);
   variables->add(*Lambdacm_PY);
   variables->add(*Lambdacm_PZ);
   variables->add(*totCandidates);
   std::cout << "adding data to RooDataSet" << std::endl;
   RooDataSet *datain = new RooDataSet("","",tin,*variables);
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
   RooAbsData::setDefaultStorageType(RooAbsData::Tree);
   RooDataSet* dataNew = new RooDataSet("dataNew","dataNew",data,*data->get());
   dataNew->convertToTreeStore();
   const TTree* tree = dataNew->tree();

   TFile* outputFile = TFile::Open("weighted.root","RECREATE");
   outputFile->WriteObject(tree,"weightedtree"); 

}


