import ROOT
import numpy as np

n_part = [383.4,331.2,109]
n_part_err = [17.8,19.6,26.6]

CENTRALITY_CLASSES = [[0,5], [5,10], [30,50]]
centrality_colors = [ROOT.kOrange+7, ROOT.kAzure+4, ROOT.kTeal+4]
centrality_colors_area = [ROOT.kOrange-9, ROOT.kAzure-9, ROOT.kGreen-8]

Z_PB208 = 82
A_PB208 = 208

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)


in_file = ROOT.TFile("FinalPlot3D.root")
in_file_glauber_std = ROOT.TFile("GlauberMcCurves/standard.root")
in_file_glauber_mod = ROOT.TFile("GlauberMcCurves/modified.root")
cMuB = in_file.Get("cMuBCent")
cMuI = in_file.Get("cMuICent")

out_file = ROOT.TFile("MuQtoMuBplot.root","recreate")
canv = ROOT.TCanvas("canv","canv",600,500)
h_frame = ROOT.TH2F("h_frame"," ",1,0.,450,1,-0.3,.7)
h_frame.GetXaxis().SetNdivisions(10)
h_frame.GetYaxis().SetNdivisions(10)

graph = ROOT.TGraphErrors()
h_frame.GetYaxis().SetTitle("#mu_{Q}/#mu_{B}")
h_frame.GetXaxis().SetTitle("#LT #it{N}_{part} #GT")
h_frame.Draw()
graph_corr = ROOT.TGraphErrors()
points = []

gMuB = cMuB.GetPrimitive("MuBCent")
gMuI = cMuI.GetPrimitive("MuICent")
gMuBCorr = cMuB.GetPrimitive("gMuBCentCorrUnc")
gMuICorr = cMuI.GetPrimitive("gMuICentCorrUnc")

for i_cent, cent in enumerate(CENTRALITY_CLASSES):
    muB = gMuB.GetPointY(i_cent)
    muB_err = gMuB.GetErrorY(i_cent)
    muB_err_corr = gMuBCorr.GetErrorY(i_cent)
    muI = gMuI.GetPointY(i_cent)
    muI_err = gMuI.GetErrorY(i_cent)
    muI_err_corr = gMuICorr.GetErrorY(i_cent)
    corr_mat = in_file.Get(f"CovMat_{cent[0]}_{cent[1]}")
    corr = corr_mat.GetBinContent(2,1)
    cov = corr*muB_err*muI_err
    graph.AddPoint(n_part[i_cent],muI/muB)
    graph.SetPointError(i_cent,0.,np.sqrt(((muI/(muB**2))**2)*(muB_err**2)+((1/muB)**2)*(muI_err**2)-(2.*muI/(muB**3))*cov)) # n_part_err[i_cent]
    print(f"E = {np.sqrt(((muI/(muB**2))**2)*(muB_err**2)+((1/muB)**2)*(muI_err**2)-(2.*muI/(muB**3))*cov)}, muB = {muB}, errmuB = {muB_err}, muQ = {muI}, errmuQ = {muI_err}, cov = {cov}")

for i_cent, cent in enumerate(CENTRALITY_CLASSES):
    muB = gMuB.GetPointY(i_cent)
    muB_err = gMuB.GetErrorY(i_cent)
    muB_err_corr = gMuBCorr.GetErrorY(i_cent)
    muI = gMuI.GetPointY(i_cent)
    muI_err = gMuI.GetErrorY(i_cent)
    muI_err_corr = gMuICorr.GetErrorY(i_cent)
    corr_mat = in_file.Get(f"CovMat_{cent[0]}_{cent[1]}")
    corr = corr_mat.GetBinContent(2,1)
    cov = corr*muB_err_corr*muI_err_corr
    graph_corr.AddPoint(n_part[i_cent],muI/muB)
    graph_corr.SetPointError(i_cent,n_part_err[i_cent],np.sqrt(((muI/(muB**2))**2)*(muB_err_corr**2)+((1/muB)**2)*(muI_err_corr**2)-(2.*muI/(muB**3))*cov))
    print(f"E = {np.sqrt(((muI/(muB**2))**2)*(muB_err_corr**2)+((1/muB)**2)*(muI_err_corr**2)-(2.*muI/(muB**3))*cov)}, muB = {muB}, errmuB = {muB_err}, muQ = {muI}, errmuQ = {muI_err}, cov = {cov}")

c_std = in_file_glauber_std.Get("Canvas_1")
c_mod = in_file_glauber_mod.Get("Canvas_1")
h_std = c_std.GetPrimitive("htemp_1")
h_mod = c_mod.GetPrimitive("htemp_1")
g_std = ROOT.TGraphErrors(h_std)
g_mod = ROOT.TGraphErrors(h_mod)
g_std.SetFillColor(ROOT.kBlue)
g_std.SetMarkerSize(0.1)
g_std.SetMarkerColor(ROOT.kBlue)
g_mod.SetFillColor(ROOT.kRed)
g_mod.SetMarkerSize(0.1)
g_mod.SetMarkerColor(ROOT.kRed)

graph.SetMarkerStyle(20)
graph.SetMarkerSize(1.2)
graph.SetMarkerColor(ROOT.kRed)
graph.SetLineColor(ROOT.kRed)
graph.SetFillColor(0)
graph.SetFillStyle(0)
graph_corr.SetMarkerColor(ROOT.kRed-9)
graph_corr.SetLineWidth(0)
graph_corr.SetMarkerSize(0)
graph_corr.SetLineColor(0)
graph_corr.SetFillColor(ROOT.kRed-9)

# f = ROOT.TF1("f","pol0")
# f.FixParameter(0,0.4)
# graph.Fit(f)
# chi2 = f.GetChisquare()


#print(f"chi2 = {chi2}")
canv.cd()

# add labels
text = ROOT.TLatex()
text.SetTextFont(43)
text.SetTextSize(28)
text.DrawLatex(36, 0.58 ,"ALICE")
text.SetTextSize(24)
text.DrawLatex(36, 0.48,"Pb-Pb #sqrt{s_{NN}}=5.02 TeV")

# legend
leg = ROOT.TLegend(0.399666,0.191579,0.623746,0.452632)
leg.SetTextFont(43)
leg.SetTextSize(22)
leg.AddEntry(graph,"Uncorr. uncert.","pe")
leg.AddEntry(graph_corr,"Corr. uncert.","f")
leg.AddEntry(g_std,"TGlauberMC, std. neutron skin","l")
leg.AddEntry(g_mod,"TGlauberMC, inc. neutron skin","l")
leg.Draw("same")

graph_corr.Draw("pe5same")
graph.Draw("pesame")
g_std.Draw("e3same")
g_mod.Draw("e3same")
out_file.cd()
canv.Write()
canv.Print("MuQtoMuB.pdf")
out_file.Close()