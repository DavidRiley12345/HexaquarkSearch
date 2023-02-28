import ROOT as r
import math
import sys



filename = __file__.replace('.py','')

#use the lhcbstyle guide - comment out if not needed
r.gROOT.ProcessLine(".x /home/ppe/d/driley/lhcbStyle.C")
r.gROOT.ForceStyle()

currstyle = r.gROOT.GetStyle("lhcbStyle")
currstyle.SetOptTitle(1)

#load file names into a string to be read into dataframe
names = r.std.vector('string')()
for n in ['/data/lhcb01/mwhitehead/LcLcpi_2018_MU.root','/data/lhcb01/mwhitehead/LcLcpi_2018_MD.root','/data/lhcb01/mwhitehead/LcLcpi_2017_MU.root','/data/lhcb01/mwhitehead/LcLcpi_2017_MD.root']: names.push_back(n)

df0 = r.RDataFrame('B2LcLcpiOS/DecayTree',names)

#define the invariant mass as a new column in the tree
df1 = df0.Define("Lambdacp_invMass", "sqrt(pow(Lambdacp_p_PE+Lambdacp_K_PE+Lambdacp_pi_PE,2)-pow(Lambdacp_p_PX+Lambdacp_K_PX+Lambdacp_pi_PX,2)-pow(Lambdacp_p_PY+Lambdacp_K_PY+Lambdacp_pi_PY,2)-pow(Lambdacp_p_PZ+Lambdacp_K_PZ+Lambdacp_pi_PZ,2))")

#filter the daughter nuclei based on the Prob variables
df2 = df1.Filter("Lambdacp_p_ProbNNp > 0.7")
df3 = df2.Filter("Lambdacp_K_ProbNNk > 0.6")
df4 = df3.Filter("Lambdacp_pi_ProbNNpi > 0.4")


#create histogram model and fill using the data from the desired column
model_M_K_P = r.RDF.TH1DModel("{}".format(filename), "{}".format(filename), 100,2170.,2400.);

h = df4.Histo1D(model_M_K_P, "Lambdacp_invMass");

hist=h.Clone();
hist.SetTitle("#font[12]{pK#pi Invariant Mass plot (OS_Lambdacp)}")

upperline = r.TLine(2305,0,2305,225000)
upperline.SetLineStyle(2)

lowerline = r.TLine(2270,0,2270,225000)
lowerline.SetLineStyle(2)



hist.GetYaxis().SetTitle("#font[12]{candidates}")
hist.GetXaxis().SetTitle("#font[12]{m(pK#pi)}[MeV]") 
c = r.TCanvas("c", "c", 650, 450)
c.SetLeftMargin(0.13)
c.SetRightMargin(0.12)
c.SetTopMargin(0.12)
c.cd()

hist.Draw("h")
upperline.Draw()
lowerline.Draw()

#save as both png and ROOT macro
c.SaveAs("{}.C".format(filename))
c.SaveAs("{}.png".format(filename))
