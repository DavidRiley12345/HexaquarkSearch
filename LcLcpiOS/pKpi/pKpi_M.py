import ROOT as r
import math
import sys

filename = __file__.strip('.py')
'''
r.gROOT.ProcessLine(".x lhcbStyle.C")
r.gROOT.ForceStyle()
'''
#fin = r.TFile.Open('/data/lhcb01/mwhitehead/LcLcpi_2018_MD.root' , '/data/lhcb01/mwhitehead/LcLcpi_2018_MU.root')
#tin = fin.Get('B2LcLcpiOS/DecayTree')
#df0 = r.RDataFrame(tin)

names = r.std.vector('string')()
for n in ['/data/lhcb01/mwhitehead/LcLcpi_2018_MU.root','/data/lhcb01/mwhitehead/LcLcpi_2018_MD.root']: names.push_back(n)

df0 = r.RDataFrame('B2LcLcpiOS/DecayTree',names)

df1 = df0.Define("Lambdacp_invMass", "sqrt(pow(Lambdacp_p_PE+Lambdacp_K_PE+Lambdacp_pi_PE,2)-pow(Lambdacp_p_PX+Lambdacp_K_PX+Lambdacp_pi_PX,2)-pow(Lambdacp_p_PY+Lambdacp_K_PY+Lambdacp_pi_PY,2)-pow(Lambdacp_p_PZ+Lambdacp_K_PZ+Lambdacp_pi_PZ,2))")

df2 = df1.Filter("Lambdacp_p_ProbNNp > 0.6")
df3 = df2.Filter("Lambdacp_K_ProbNNk > 0.6")
df4 = df3.Filter("Lambdacp_pi_ProbNNpi > 0.6")



model_M_K_P = r.RDF.TH1DModel("{}".format(filename), "{}".format(filename), 150,2150.,2430.);

h = df1.Histo1D(model_M_K_P, "Lambdacp_invMass");
ProbNNcuts = df4.Histo1D(model_M_K_P, "Lambdacp_invMass");


hist=h.Clone();
hist.SetName("hist");

hist.GetXaxis().SetTitle("#font[12]{P(K)}[MeV]")

filename = __file__.strip('.py')

c = r.TCanvas("c", "c", 650, 450)
c.SetLeftMargin(0.12)
c.SetRightMargin(0.12)
c.cd()
hist.Draw("h")
ProbNNcuts.Draw("SAME")

c.SaveAs("{}.C".format(filename))
c.SaveAs("{}.png".format(filename))
