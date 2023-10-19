import ROOT as r
import math
from array import array

cuts = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
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
    df2 = df1.Filter("Lambdacp_IPCHI2_OWNPV < {}".format(cut))
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

c1 = r.TCanvas('c1','Purity with cut on Lambdacp_IPCHI2',200,10,700,500)
graph = r.TGraph(len(cuts), array('d',cuts), array('d',purities))
graph.SetTitle("Purity with cut on Lambdacp_IPCHI2;Cut Value;Purity")
graph.Draw()
c1.SaveAs("IPCHI_fitplotpurity.png")

c2 = r.TCanvas('c2','Significance with cut on Lambdacp_IPCHI2',200,10,700,500)
graph = r.TGraph(len(cuts), array('d',cuts), array('d',sigs))
graph.SetTitle("Significance cut on Lambdacp_IPCHI2;Cut Value;Significance")
graph.Draw()
c2.SaveAs("FDCHI_fitplotsignificance.png")