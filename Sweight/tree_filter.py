import ROOT as r
import math
import sys
from array import array

filename = __file__.replace('.py','')

#fin = r.TFile.Open('/data/lhcb01/mwhitehead/LcLcpi_2018_MD.root')
#tin = fin.Get('B2LcLcpiOS/DecayTree')

names = r.std.vector("string")()

for n in ['/data/lhcb01/mwhitehead/LcLcpi_2018_MU.root','/data/lhcb01/mwhitehead/LcLcpi_2018_MD.root','/data/lhcb01/mwhitehead/LcLcpi_2017_MU.root','/data/lhcb01/mwhitehead/LcLcpi_2017_MD.root']: names.push_back(n)

df0 = r.RDataFrame("B2LcLcpiOS/DecayTree",names)
df1 = df0.Define("Lambdacp_invMass","Lambdacp_M")
df11= df1.Filter("Lambdacp_M < 2350 && Lambdacp_M > 2220")
df2 = df11.Filter("Lambdacp_K_ProbNNk > 0.4 && Lambdacm_K_ProbNNk > 0.4 && Lambdacp_p_ProbNNp > 0.4 && Lambdacm_p_ProbNNp > 0.4 && Lambdacp_pi_ProbNNpi > 0.4 && Lambdacm_pi_ProbNNpi > 0.4")

model_LcM = r.RDF.TH1DModel("model_LcM","Lambdacp_M", 50, 2220,2350)

# create snapshot of dataset once its been cut

col_list = r.std.vector("std::string")()

col_list.push_back('Lambdacp_M')
col_list.push_back('Lambdacp_P')
col_list.push_back('Lambdacp_PT')

col_list.push_back('Lambdacm_M')
col_list.push_back('Lambdacm_P')
col_list.push_back('Lambdacm_PT')

col_list.push_back('Lambdacp_p_P')
col_list.push_back('Lambdacp_p_PT')

col_list.push_back('Lambdacp_K_P')
col_list.push_back('Lambdacp_K_PT')

col_list.push_back('Lambdacp_pi_P')
col_list.push_back('Lambdacp_pi_PT')

col_list.push_back('pi_P')
col_list.push_back('pi_PT')
col_list.push_back('pi_PE')
col_list.push_back('pi_PX')
col_list.push_back('pi_PY')
col_list.push_back('pi_PZ')

col_list.push_back('Lambdacp_PE')
col_list.push_back('Lambdacp_PX')
col_list.push_back('Lambdacp_PY')
col_list.push_back('Lambdacp_PZ')

col_list.push_back('Lambdacm_PE')
col_list.push_back('Lambdacm_PX')
col_list.push_back('Lambdacm_PY')
col_list.push_back('Lambdacm_PZ')

col_list.push_back('totCandidates')

df2.Snapshot("Lcm","cut.root",col_list)
'''
mass = r.RooRealVar("mass","",2220,2350)
arglist = r.RooArgList(mass)

# read datafrom snapshot into a histogram for fitting and a RooDataSet for Sweighting

filein = r.TFile.Open('cut.root')
treein = filein.Get("Lcm_cut")
rdfin =  r.RDataFrame(treein)
hin = rdfin.Histo1D(model_LcM,"Lambdacp_M")
hinclone = hin.Clone()
cuthist = r.RooDataHist("hist","hist",arglist,hinclone)

x = r.RooRealVar("Lambdacp_M","Lambdacp_M",2220,2350)
massSet = r.RooArgSet(x)

dataset = r.RooDataSet("","",treein,massSet)

# set up variables for RooFit

mean = r.RooRealVar("mean","",2286,2250.0,2300.0)
sigma = r.RooRealVar("sigma","",5,0.1,15.0)
Gauss = r.RooGaussian("signal","",mass,mean,sigma)

slope =r.RooRealVar("slope","",-0.002,-10,10)
bkgPDF = r.RooExponential("comboPDF","",mass,slope)

nsig = r.RooRealVar("nsig","",10000 ,-1000,1000000)
nbkg = r.RooRealVar("nbkg","",100000,-1000,1000000)

totalPDF = r.RooAddPdf("totalPDF","",r.RooArgList(Gauss,bkgPDF),r.RooArgList(nsig,nbkg))
totalPDF.fitTo(dataset,r.RooFit.Extended());

c3 = r.TCanvas("c3","Fit",800,600);
c3.Divide(2,1)
c3.cd(1)
plot1 = mass.frame(50);

#dataset.plotOn(plot1);
cuthist.plotOn(plot1);
totalPDF.plotOn(plot1);

plot1.GetYaxis().SetTitle("Candidates / X MeV/c^{2}");
plot1.GetXaxis().SetTitle("m(#Lambda_{c}) (MeV/c^{2})");
plot1.Draw();
c3.SaveAs("fitresult.C");
c3.SaveAs("fitresult.png");

print(nsig.getValV(),nbkg.getValV())

mean.setConstant()
sigma.setConstant()
slope.setConstant()

r.RooMsgService.instance().setSilentMode(True)

sData = r.RooStats.SPlot("splot","splot",dataset,totalPDF,r.RooArgList(nsig,nbkg))
sData.GetSWeightVars().Print()
print("yield from sweight", sData.GetYieldFromSWeight("Lambdacp_M"))

print(nsig.getValV(),nbkg.getValV())

c3.cd(2)
plot2 = mass.frame(50);
#sig_sw = r.RooDataSet("sw_name","sw_title", dataset, dataset.get(), 0, "nsig_sw")

'''
