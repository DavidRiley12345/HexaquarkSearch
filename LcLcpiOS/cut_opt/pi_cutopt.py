import ROOT as r
import math
from array import array

#use the lhcbstyle guide - comment out if not needed
r.gROOT.ProcessLine(".x /home/ppe/d/driley/lhcbStyle.C")
r.gROOT.ForceStyle()

currstyle = r.gROOT.GetStyle("lhcbStyle")
currstyle.SetOptTitle(1)
currstyle.SetTitleH(0.08)

cuts = [0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,0.99]
pid = 'pi'
purities = []
sigs = []

fin = r.TFile.Open('/data/lhcb01/mwhitehead/LcLcpi_2018_MD.root')
tin = fin.Get('B2LcLcpiOS/DecayTree')
df0 = r.RDataFrame(tin)
df1 = df0.Define("Lambdacp_invMass", "Lambdacp_M")

model_LcM = r.RDF.TH1DModel("model_LcM","Lambdacp_M", 50, 2220,2350)
mass = r.RooRealVar("mass","",2220,2350)
mean = r.RooRealVar("mean","",2287.,2270.0,2290.0)
sigma = r.RooRealVar("sigma","",4.,0.0,20.0)
Gauss = r.RooGaussian("signal","",mass,mean,sigma)
slope =r.RooRealVar("slope","",-1000,1000)
bkgPDF = r.RooExponential("comboPDF","",mass,slope)
nsig = r.RooRealVar("nsig","",-1000,1000000)
nbkg = r.RooRealVar("nbkg","",-1000,10000000)
totalPDF = r.RooAddPdf("totalPDF","",r.RooArgList(Gauss,bkgPDF),r.RooArgList(nsig,nbkg))

for cut in cuts:
    df2 = df1.Filter("Lambdacp_{0}_ProbNN{0} > {1}".format(pid,cut))
    h = df2.Histo1D(model_LcM, "Lambdacp_invMass")
    h3 = h.Clone()
    data = r.RooDataHist("data","data",r.RooArgList(mass), h3)
    totalPDF.fitTo(data,r.RooFit.Extended())
    sig = nsig.getValV()
    bkg = nbkg.getValV()
    purity = sig / (sig+bkg)
    significance = sig / math.sqrt(sig + bkg)
    purities.append(purity)
    sigs.append(significance)

c1 = r.TCanvas('c1','Purity with cut on Prob {}'.format(pid),200,10,700,500)
graph = r.TGraph(len(cuts), array('d',cuts), array('d',purities))
graph.SetTitle("Purity with Prob {} cut;Cut Value;Purity".format(pid))
graph.SetMarkerStyle(0)
graph.Draw()
c1.SaveAs("{}_fitplotpurity_exp.png".format(pid))

c2 = r.TCanvas('c2','Significance With Cut on Prob {}'.format(pid),200,10,700,500)
graph = r.TGraph(len(cuts), array('d',cuts), array('d',sigs))
graph.SetTitle("Significance with Prob {} cut;Cut Value;Significance".format(pid))
graph.SetMarkerStyle(0)
graph.Draw()
c2.SaveAs("{}_fitplotsignificance_exp.png".format(pid))