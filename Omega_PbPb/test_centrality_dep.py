import ROOT
import numpy as np

cent = [[0,5],[5,10],[30,50]]

f_ratios = ROOT.TFile("ratio_sys.root")
g = ROOT.TGraphErrors()
for i_c, c in enumerate(cent):
    ratio = f_ratios.Get(f"h_ratio_{c[0]}_{c[1]}")
    sys = f_ratios.Get(f"h_sys;{i_c+1}")
    ratio_val = ratio.GetFunction("pol0").GetParameter(0)
    ratio_err = ratio.GetFunction("pol0").GetParError(0)
    ratio_sys = sys.GetStdDev()
    g.AddPoint(0.5*(c[0]+c[1]),ratio_val)
    g.SetPointError(i_c,0,np.sqrt(ratio_err*ratio_err+ratio_sys*ratio_sys))
g.Fit("pol0")
g.GetXaxis().SetTitle("Centrality (%)")
g.GetYaxis().SetTitle("#bar{#Omega}^{+} / #Omega^{-}")
c = ROOT.TCanvas("c","c")
g.Draw()
c.Print("c.png")