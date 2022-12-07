import ROOT

def InvertPoints(g,a,b):
    buffer = []
    buffer.append(g.GetPointX(a))
    buffer.append(g.GetErrorX(a))
    buffer.append(g.GetPointY(a))
    buffer.append(g.GetErrorY(a))
    g.SetPoint(a,g.GetPointX(b),g.GetPointY(b))
    g.SetPointError(a,g.GetErrorX(b),g.GetErrorY(b))
    g.SetPoint(b,buffer[0],buffer[2])
    g.SetPointError(b,buffer[1],buffer[3])
    return

staggering = False

ROOT.gStyle.SetOptStat(0)
cent = [[0,5],[5,10],[10,30],[30,50],[50,90]]

f_names = ["FinalPlot3D_new.root","FinalPlot3D_new_fixmuQ_nopions.root"]
file = [ROOT.TFile(f) for f in f_names]
cvs = [f.Get("cMuQuncorr") for f in file]
g = [c.FindObject("MuQCent") for c in cvs]

c = ROOT.TCanvas()
c.cd()
frame = ROOT.TH2D("frame",";Centrality (%);#mu_{#it{Q}} (MeV)",1,0,90,1,-1,.6)
frame.Draw()
#gCorr[0].Draw("samepe5")
for p in range(g[1].GetN()):
    if staggering:
        g[1].SetPointX(p,g[1].GetPointX(p)+2)
    for gg in g:
        gg.SetPointError(p,0.,gg.GetErrorY(p))
InvertPoints(g[1],2,3)
#gCorr[1].Draw("samepe5")
g[0].Draw("pesame")
g[0].GetFunction("pol0").Delete()
g[1].SetLineColor(ROOT.kBlue)
g[1].SetFillColor(ROOT.kBlue)
g[1].SetMarkerColor(ROOT.kBlue)
g[1].Draw("e3same")
g[1].GetFunction("pol0").Delete()

l = ROOT.TLegend(0.256892,0.151304,0.442356,0.257391)
l.SetTextFont(44)
l.SetTextSize(22)
l.AddEntry(g[0],"Thermal-FIST, #it{T}_{ch}=155 MeV, #mu_{#it{S}} constrained","pe")
l.AddEntry(g[1],"Thermal-FIST, #it{T}_{ch}=155 MeV, #mu_{#it{S}} and #mu_{#it{Q}} constrained","l")
l.Draw("same")

t = ROOT.TLatex()
t.SetTextFont(44)
t.SetTextSize(40)
t.DrawLatex(30.8736,0.4,"ALICE")
t.SetTextSize(35)
t.DrawLatex(30.8736,0.2,"Pb-Pb #sqrt{#it{s}_{NN}}=5.02 TeV")
t.SetTextSize(25)
t.DrawLatex(17.,-.671795,"0.9 MeV correlated uncertainty not shown")

o = ROOT.TFile("muQvsCent.root","recreate")
c.Write()
c.Print("MuQvsCent.pdf")