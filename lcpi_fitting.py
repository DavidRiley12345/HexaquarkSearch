import ROOT as r
import math
import sys
from array import array
import numpy as np

filename = __file__.replace('.py','')

r.EnableImplicitMT()

canvas = 1
savedraw = 1

names = r.std.vector('string')()
for n in ['/data/lhcb01/mwhitehead/LcLcpi_2018_MU.root','/data/lhcb01/mwhitehead/LcLcpi_2018_MD.root','/data/lhcb01/mwhitehead/LcLcpi_2017_MU.root','/data/lhcb01/mwhitehead/LcLcpi_2017_MD.root']: names.push_back(n)

df0 = r.RDataFrame('B2LcLcpiOS/DecayTree',names)
'''
fin = r.TFile.Open('/data/lhcb01/mwhitehead/LcLcpi_2018_MD.root')
tin = fin.Get('B2LcLcpiOS/DecayTree')
df0 = r.RDataFrame(tin)
'''
df1 = df0.Define("Lambdacp_invMass", "sqrt(pow(Lambdacp_PE+pi_PE,2)-pow(Lambdacp_PX+pi_PX,2)-pow(Lambdacp_PY+pi_PY,2)-pow(Lambdacp_PZ+pi_PZ,2))")
df2 = df1.Filter("Lambdacp_p_ProbNNp > 0.7 && Lambdacp_K_ProbNNk > 0.6 && Lambdacp_pi_ProbNNpi > 0.4 && Lambdacm_p_ProbNNp > 0.7 && Lambdacm_K_ProbNNk > 0.6 && Lambdacm_pi_ProbNNpi > 0.4")

mass = r.RooRealVar("mass","",2390.,3000.)
arglist = r.RooArgList(mass)

mean = r.RooRealVar("mean","",2453.75,2453.,2455.)
sigma = r.RooRealVar("sigma","",1.,15)
Gauss = r.RooGaussian("signal","",mass,mean,sigma)

mean2 = r.RooRealVar("mean2","",2518.4,2500,2525)
sigma2 = r.RooRealVar("sigma2","",0.1,30.0)
Gauss2 = r.RooGaussian("signal2","",mass,mean2,sigma2)

mix = r.RooRealVar('mix','',1.,10000000.)
mix2 = r.RooRealVar('mix2','',1.,10000000.)

totalGauss = r.RooAddPdf("totalGauss","",r.RooArgList(Gauss,Gauss2),r.RooArgList(mix,mix2))

#slope =r.RooRealVar("slope","",-1000.,1000.)
#bkgPDF = r.RooExponential("comboPDF","",mass,slope)
#slopearglist = r.RooArgList(slope)
#bkgPDF = r.RooPolynomial("comboPDF","",mass,slopearglist,1)

dm0= r.RooRealVar('dm0','',2390)

a = r.RooRealVar('a','a',100000.,0.,10000000.)
b = r.RooRealVar('b','b',-.005,-0.9,.0)
c = r.RooRealVar('c','c',-.005,-0.9,.0)
bg = r.RooDstD0BG('bg','',mass,dm0,a,b,c)

m0 = r.RooRealVar('m0','',8000)
argparam = r.RooRealVar('argparam','',0.,10.)
argus = r.RooArgusBG('bg2','',mass,m0,argparam)

atest = r.RooRealVar('atest','atest',10.)
btest = r.RooRealVar('btest','btest',1.)
ctest = r.RooRealVar('ctest','ctest',1.)
test = r.RooDstD0BG('test','',mass,dm0,atest,btest,ctest)

nsig = r.RooRealVar("nsig","",-1000,10000000)
nbkg = r.RooRealVar("nbkg","",-1000,100000000)

totalPDF = r.RooAddPdf("totalPDF","",r.RooArgList(totalGauss,bg),r.RooArgList(nsig,nbkg))

model_LcM = r.RDF.TH1DModel("model_LcM","Lambdacp_M", 100, 2390.,3000.)
h = df2.Histo1D(model_LcM, "Lambdacp_invMass")
h3 = h.Clone()
	        
data = r.RooDataHist("data","data",arglist,h3)	
totalPDF.fitTo(data,r.RooFit.Extended());

if canvas == True:
	c3 = r.TCanvas("c3","Fit",800,600);
	c3.SetLeftMargin(0.12)
	c3.SetRightMargin(0.12)
	plot2 = mass.frame(50);
	data.plotOn(plot2);
	totalPDF.plotOn(plot2,r.RooFit.LineColor(1));
	totalPDF.plotOn(plot2,r.RooFit.Components('signal'),r.RooFit.LineColor(2));
	totalPDF.plotOn(plot2,r.RooFit.Components('signal2'),r.RooFit.LineColor(3));
	totalPDF.plotOn(plot2,r.RooFit.Components('bg'),r.RooFit.LineColor(4));
	plot2.SetTitle('LcPi Invariant Mass Fit')
	plot2.GetYaxis().SetTitle("Candidates / X MeV/c^{2}");
	plot2.GetXaxis().SetTitle("m(#Lambda_{c}#pi) (MeV/c^{2})");
	plot2.Draw();
	if savedraw == True:
		#c3.SaveAs("{}_fitresult.C".format(cutval));
		c3.SaveAs("Fit.png");


'''
rangetest = r.RooRealVar('testmass','',0,100000)
atest = r.RooRealVar('atest','atest',100)
btest = r.RooRealVar('btest','btest',-0.1)
ctest = r.RooRealVar('ctest','ctest',-0.05)
test = r.RooDstD0BG('test','',rangetest,dm0,atest,btest,ctest)
c1 = r.TCanvas('c1','',800,600)
plot1 = rangetest.frame(50)
test.plotOn(plot1)
plot1.Draw()
c1.SaveAs("test.png")
'''
	
