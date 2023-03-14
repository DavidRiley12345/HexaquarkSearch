import ROOT as r
import math
import sys



filename = __file__.replace('.py','')

#use the lhcbstyle guide - comment out if not needed
r.gROOT.ProcessLine(".x /home/ppe/d/driley/lhcbStyle.C")
r.gROOT.ForceStyle()

currstyle = r.gROOT.GetStyle("lhcbStyle")
currstyle.SetOptTitle(1)
#currstyle.SetTitleW(0.)
#currstyle.SetTitleH(0.1)
 

#load file names into a string to be read into dataframe
names = r.std.vector('string')()
for n in ['/data/lhcb01/mwhitehead/LcLcpi_2018_MU.root','/data/lhcb01/mwhitehead/LcLcpi_2018_MD.root','/data/lhcb01/mwhitehead/LcLcpi_2017_MU.root','/data/lhcb01/mwhitehead/LcLcpi_2017_MD.root']: names.push_back(n)

df0 = r.RDataFrame('B2LcLcpiOS/DecayTree',names)

#define the invariant mass as a new column in the tree
df1 = df0.Define("Lambdacp_invMass", "sqrt(pow(Lambdacp_p_PE+Lambdacp_K_PE+Lambdacp_pi_PE,2)-pow(Lambdacp_p_PX+Lambdacp_K_PX+Lambdacp_pi_PX,2)-pow(Lambdacp_p_PY+Lambdacp_K_PY+Lambdacp_pi_PY,2)-pow(Lambdacp_p_PZ+Lambdacp_K_PZ+Lambdacp_pi_PZ,2))")

#create histogram model and fill using the data from the desired column
model_M_K_P = r.RDF.TH1DModel("{}".format(filename), "{}".format(filename), 100,2170.,2400.);

h = df1.Histo1D(model_M_K_P, "Lambdacp_invMass");

hist=h.Clone();
hist.SetTitle("#font[12]{#scale[1.5]{pK#pi Invariant Mass plot (OS_Lambdacp)}}")

upperline = r.TLine(2300,0,2300,225000)
upperline.SetLineStyle(2)

lowerline = r.TLine(2270,0,2270,225000)
lowerline.SetLineStyle(2)

hist.GetYaxis().SetTitle("#font[12]{candidates}")
hist.GetXaxis().SetTitle("#font[12]{m(pK#pi)}[MeV]")
hist.GetYaxis().SetTitleOffset(1.2) 
c = r.TCanvas("c", "c", 1300, 900)
c.SetLeftMargin(0.16)
c.SetRightMargin(0.12)
c.SetTopMargin(0.12)
c.cd()

hist.Draw("h")
#save as both png and ROOT macro
c.SaveAs("{}.C".format(filename))
c.SaveAs("{}.png".format(filename))
