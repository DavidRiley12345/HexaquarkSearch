import ROOT as r
import math
import sys
from array import array
import numpy as np
filename = __file__.replace('.py','')



start = 0.0
end = 1.
numb = 9

#cuts = np.linspace(start,end,num=numb)
#step = np.linspace(start,end,num=numb,retstep = True)[1]
#cuts = array('d',[0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95])
cuts = array('d',[0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95])
#cuts= array('d',[0.05,0.1])
nbins = 10
cutmax = max(cuts)
cutmin = min(cuts)

canvas = 1
savedraw = 1

fin = r.TFile.Open('/data/lhcb01/mwhitehead/LcLcpi_2018_MD.root')
tin = fin.Get('B2LcLcpiOS/DecayTree')
df0 = r.RDataFrame(tin)
df1 = df0.Define("Lambdacp_invMass", "Lambdacp_M")

Hist2DSig = r.TH2D('','',nbins,0,1,nbins,0,1);
Hist2DPurity = r.TH2D('','',nbins,0,1,nbins,0,1);

mass = r.RooRealVar("mass","",2220,2350)
arglist = r.RooArgList(mass)

mean = r.RooRealVar("mean","",2287.,2270.0,2300.0)
sigma = r.RooRealVar("sigma","",4.,0.5,20.0)
Gauss = r.RooGaussian("signal","",mass,mean,sigma)

mean2 = r.RooRealVar("mean2","",2287.,2270.0,2300.0)
sigma2 = r.RooRealVar("sigma2","",4.,0.5,20.0)
Gauss2 = r.RooGaussian("signal2","",mass,mean,sigma2)

mix = r.RooRealVar('mix','',-1,1)
mix2 = r.RooRealVar('mix2','',10000,-1000.,10000000.)
totalGauss = r.RooAddPdf("totalGauss","",r.RooArgList(Gauss,Gauss2),r.RooArgList(mix))

slope =r.RooRealVar("slope","",-1000.,1000.)
bkgPDF = r.RooExponential("comboPDF","",mass,slope)
#slopearglist = r.RooArgList(slope)
#bkgPDF = r.RooPolynomial("comboPDF","",mass,slopearglist,1)

nsig = r.RooRealVar("nsig","",-1000,1000000)
nbkg = r.RooRealVar("nbkg","",-1000,1000000)

totalPDF = r.RooAddPdf("totalPDF","",r.RooArgList(totalGauss,bkgPDF),r.RooArgList(nsig,nbkg))

for pcut in cuts:
	for kcut in cuts:
		
		pcutval = pcut
		kcutval = kcut
		
		currBinSig = Hist2DSig.FindBin(pcutval,kcutval)
		currBinPurity = Hist2DPurity.FindBin(pcutval,kcutval)		

		df2 = df1.Filter("Lambdacp_p_ProbNNp > {}".format(pcutval)) 
		df3 = df2.Filter("Lambdacp_K_ProbNNk > {}".format(kcutval))

		model_LcM = r.RDF.TH1DModel("model_LcM","Lambdacp_M", 100, 2220,2350)

		h = df3.Histo1D(model_LcM, "Lambdacp_invMass")
		h3 = h.Clone()
	        
		data = r.RooDataHist("data","data",arglist,h3)	
		totalPDF.fitTo(data,r.RooFit.Extended());
		
		sig = nsig.getValV()
		bkg = nbkg.getValV()
		purity = sig / (sig+bkg)
		significance = sig / math.sqrt(sig + bkg)
			
		Hist2DSig.SetBinContent(currBinSig,significance)
		Hist2DPurity.SetBinContent(currBinPurity,purity)
		print("pcut val: {0} ,kcut val: {1}".format(pcutval,kcutval))
		print("")
		
		print("Signal : {}".format(sig))
		print("Backgr : {}".format(bkg))
		print("")
		print("Purity : {}".format(purity))
		print("Signif : {}".format(significance))
		
		if canvas == True:
			c3 = r.TCanvas("c3","Fit for p cut at {0}{1}".format(pcutval,kcutval),800,600);
			c3.SetLeftMargin(0.12)
			c3.SetRightMargin(0.12)
			plot2 = mass.frame(50);
			data.plotOn(plot2);
			totalPDF.plotOn(plot2,r.RooFit.Components("signal"),r.RooFit.LineColor(2))
			totalPDF.plotOn(plot2,r.RooFit.Components("signal2"),r.RooFit.LineColor(3))
			totalPDF.plotOn(plot2,r.RooFit.Components("comboPDF"),r.RooFit.LineColor(4))
			totalPDF.plotOn(plot2,r.RooFit.LineColor(1));
			plot2.SetTitle('Plot of PID cuts: {0} p {1} K'.format(pcutval,kcutval))
			plot2.GetYaxis().SetTitle("Candidates / X MeV/c^{2}");
			plot2.GetXaxis().SetTitle("m(#Lambda_{c}) (MeV/c^{2})");
			plot2.Draw();
			if savedraw == True:
				#c3.SaveAs("{}_fitresult.C".format(cutval));
				c3.SaveAs("/home/ppe/d/driley/git_HexaquarkSummer/LcLcpiOS/pKpi/cut_opt/pk_plots/p{0}k{1}_fitresult.png".format(pcutval,kcutval));

c = r.TCanvas('c','c', 650,450)
c.Divide(2,1)
c.cd(1)
Hist2DSig.SetStats(False)
Hist2DSig.SetTitle('pK significance plot;pcut;kcut')
Hist2DSig.Draw('COLZ')

c.cd(2)
Hist2DPurity.SetStats(False)
Hist2DPurity.SetTitle('pK purity plot;pcut;kcut')
Hist2DPurity.Draw('COLZ')


c.SaveAs('2D.png')
c.SaveAs('2D.C') 


 	
'''
c1 = r.TCanvas('c1','Purity with cut on Prob pk',200,10,700,500)
graph1 = r.TGraph2D(n,x,y,sigs)
graph1.SetTitle("Purity with Prob pK cut;Cut Value;Purity")
graph1.Draw('surf1')
c1.SaveAs("pK_fitplotpurity_linear.png")


c2 = r.TCanvas('c2','Significance With Cut on Prob pK',200,10,700,500)
graph2 = r.TGraph(n,cuts,sigs)
graph2.SetTitle("Significance with Prob pK cut;Cut Value;Significance")
graph2.Draw()
c2.SaveAs("pK_fitplotsignificance_linear.png")
'''
