import ROOT as r
import math
import sys



filename = __file__.replace('.py','')

#use the lhcbstyle guide - comment out if not needed
r.gROOT.ProcessLine(".x /home/ppe/d/driley/lhcbStyle.C")
r.gROOT.ForceStyle()

currstyle = r.gROOT.GetStyle("lhcbStyle")
currstyle.SetOptTitle(1)
currstyle.SetTitleH(0.1)

# set the y title margin to be wider

currstyle.SetTitleYOffset(1.3)

#load file names into a string to be read into dataframe
names = r.std.vector('string')()
for n in ['/data/lhcb01/mwhitehead/LcLcpi_2018_MU.root','/data/lhcb01/mwhitehead/LcLcpi_2018_MD.root','/data/lhcb01/mwhitehead/LcLcpi_2017_MU.root','/data/lhcb01/mwhitehead/LcLcpi_2017_MD.root']: names.push_back(n)

df0 = r.RDataFrame('B2LcLcpiOS/DecayTree',names)

#define the invariant mass as a new column in the tree
df1 = df0.Define("Lambdacp_invMass", "sqrt(pow(Lambdacp_p_PE+Lambdacp_K_PE+Lambdacp_pi_PE,2)-pow(Lambdacp_p_PX+Lambdacp_K_PX+Lambdacp_pi_PX,2)-pow(Lambdacp_p_PY+Lambdacp_K_PY+Lambdacp_pi_PY,2)-pow(Lambdacp_p_PZ+Lambdacp_K_PZ+Lambdacp_pi_PZ,2))")

#filter the daughter nuclei based on the Prob variables
df2 = df1.Filter("Lambdacp_p_ProbNNp > 0.7 && Lambdacp_K_ProbNNk > 0.6 && Lambdacp_pi_ProbNNpi > 0.4").Define("Angle_pp_pip", "(180/3.14)*(acos((Lambdacp_pi_PX*Lambdacp_p_PX + Lambdacp_pi_PY*Lambdacp_p_PY + Lambdacp_pi_PZ*Lambdacp_p_PZ) / (sqrt(pow(Lambdacp_pi_PX,2)+pow(Lambdacp_pi_PY,2)+pow(Lambdacp_pi_PZ,2)) * sqrt(pow(Lambdacp_p_PX,2)+pow(Lambdacp_p_PY,2)+pow(Lambdacp_p_PZ,2)))))").Define("Angle_pm_pim", "(180/3.14)*(acos((Lambdacm_pi_PX*Lambdacm_p_PX + Lambdacm_pi_PY*Lambdacm_p_PY + Lambdacm_pi_PZ*Lambdacm_p_PZ) / (sqrt(pow(Lambdacm_pi_PX,2)+pow(Lambdacm_pi_PY,2)+pow(Lambdacm_pi_PZ,2)) * sqrt(pow(Lambdacm_p_PX,2)+pow(Lambdacm_p_PY,2)+pow(Lambdacm_p_PZ,2)))))").Define("Angle_kp_pp", "(180/3.14)*(acos((Lambdacm_K_PX*Lambdacp_p_PX + Lambdacm_K_PY*Lambdacp_p_PY + Lambdacm_K_PZ*Lambdacp_p_PZ) / (sqrt(pow(Lambdacm_K_PX,2)+pow(Lambdacm_K_PY,2)+pow(Lambdacm_K_PZ,2)) * sqrt(pow(Lambdacp_p_PX,2)+pow(Lambdacp_p_PY,2)+pow(Lambdacp_p_PZ,2)))))").Define("Angle_km_pm", "(180/3.14)*(acos((Lambdacp_K_PX*Lambdacm_p_PX + Lambdacp_K_PY*Lambdacm_p_PY + Lambdacp_K_PZ*Lambdacm_p_PZ) / (sqrt(pow(Lambdacp_K_PX,2)+pow(Lambdacp_K_PY,2)+pow(Lambdacp_K_PZ,2)) * sqrt(pow(Lambdacm_p_PX,2)+pow(Lambdacm_p_PY,2)+pow(Lambdacm_p_PZ,2)))))").Define("Angle_kp_pip", "(180/3.14)*(acos((Lambdacm_K_PX*Lambdacp_pi_PX + Lambdacm_K_PY*Lambdacp_pi_PY + Lambdacm_K_PZ*Lambdacp_pi_PZ) / (sqrt(pow(Lambdacm_K_PX,2)+pow(Lambdacm_K_PY,2)+pow(Lambdacm_K_PZ,2)) * sqrt(pow(Lambdacp_pi_PX,2)+pow(Lambdacp_pi_PY,2)+pow(Lambdacp_pi_PZ,2)))))").Filter("Angle_pp_pip >0.05 && Angle_pm_pim >0.05 && Angle_kp_pp >0.05 && Angle_km_pm >0.05 && Angle_kp_pip >0.05").Filter("Lambdacp_K_ORIVX_CHI2 < 10 && Lambdacp_p_ORIVX_CHI2 < 10 && Lambdacp_pi_ORIVX_CHI2 < 10").Filter("LcLcpiOS_ENDVERTEX_CHI2 < 9");


#create histogram model and fill using the data from the desired column
model_M_K_P = r.RDF.TH1DModel("{}".format(filename), "{}".format(filename), 100,2220.,2350.);

h = df2.Histo1D(model_M_K_P, "Lambdacp_invMass");

hist=h.Clone();
hist.SetTitle("#font[12]{pK#pi Invariant Mass plot #Lambda_{c}^{+}}")

upperline = r.TLine(2305,0,2305,60000)
upperline.SetLineStyle(2)

lowerline = r.TLine(2270,0,2270,60000)
lowerline.SetLineStyle(2)



hist.GetYaxis().SetTitle("candidates / {0} MeV".format(int(hist.GetBinWidth(1))))
hist.GetXaxis().SetTitle("#font[12]{m(p^{+}K^{-}#pi^{+})}[MeV]") 
c = r.TCanvas("c", "c", 650, 450)
c.SetLeftMargin(0.17)
c.SetRightMargin(0.12)
c.SetTopMargin(0.12)
c.cd()

hist.Draw("h")
upperline.Draw()
lowerline.Draw()

#save as both png and ROOT macro
c.SaveAs("{}.C".format(filename))
c.SaveAs("{}.png".format(filename))
