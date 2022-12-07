import ROOT

staggering = False

ROOT.gStyle.SetOptStat(0)
cent = [[0,5],[5,10],[10,30],[30,50],[50,90]]

f_names = ["FinalPlot3D_new.root","FinalPlot3D_new_fixmuQ_nopions.root"]
file = [ROOT.TFile(f) for f in f_names]
cvs = [f.Get("cMuBCent") for f in file]
g = [c.FindObject("MuBCent") for c in cvs]
gCorr = [c.FindObject("gMuBCentCorrUnc") for c in cvs]
n = cvs[0].FindObject("NaturePoint")

c = ROOT.TCanvas()
c.cd()
frame = ROOT.TH2D("frame",";Centrality (%);#mu_{#it{B}} (MeV)",1,0,90,1,-3.6,5.)
frame.Draw()
n.Draw("samepe5")
#gCorr[0].Draw("samepe5")
g[0].Draw("pesame")
for p in range(gCorr[1].GetN()):
    if staggering:
        gCorr[1].SetPointX(p,gCorr[1].GetPointX(p)+2)
        g[1].SetPointX(p,g[1].GetPointX(p)+2)
    for gg in g:
        gg.SetPointError(p,0.,gg.GetErrorY(p))
gCorr[1].SetFillColor(ROOT.kBlue-9)
#gCorr[1].Draw("samepe5")
g[1].SetLineColor(ROOT.kBlue)
g[1].SetMarkerColor(ROOT.kBlue)
g[1].Draw("samepe")

l = ROOT.TLegend(0.26817,0.175652,0.452381,0.325217)
l.SetTextFont(44)
l.SetTextSize(22)
l.AddEntry(g[0],"Thermal-FIST, #it{T}_{ch}=155 MeV, #mu_{#it{S}} constrained","pe")
l.AddEntry(g[1],"Thermal-FIST, #it{T}_{ch}=155 MeV, #mu_{#it{S}} and #mu_{#it{Q}} constrained","pe")
l.AddEntry(n,"SHM fit, #it{Nature} #bf{561}, 321-330 (2018)")
l.Draw("same")

t = ROOT.TLatex()
t.SetTextFont(44)
t.SetTextSize(40)
t.DrawLatex(30.8736,3.91661,"ALICE")
t.SetTextSize(35)
t.DrawLatex(30.8736,2.91661,"Pb-Pb #sqrt{#it{s}_{NN}}=5.02 TeV")
t.SetTextSize(25)
t.DrawLatex(17.,-1.,"0.5 MeV correlated uncertainty not shown")

o = ROOT.TFile("muBvsCent.root","recreate")
c.Write()
c.Print("MuBvsCent.pdf")