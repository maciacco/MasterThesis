import ROOT
import numpy as np

staggering = False

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

def chi2(par,g,gCorr,gB,print_flag):
    val = []
    uncorr = []
    corr = []
    polarity = []
    cov_m = []
    par_vec = []
    for i_c in range(5):
        val.append(g.GetPointY(i_c))
        uncorr.append(np.abs(g.GetErrorY(i_c)))
        corr.append(np.abs(gCorr.GetErrorY(i_c)))
        polarity.append(np.abs(gB.GetErrorY(i_c)))
        par_vec.append(par)
    polarity_min = min(polarity)
    corr_min = min(corr)
    if print_flag:
        #print(corr)
        pass
    for i in range(5):
        cov_m.append([(polarity_min**2 + corr[i]**2) for _ in range(5)])
        for j in range(5):
            if j == i:
                continue
            cov_m[i][j] = corr[i]*corr[j] + polarity[i]*polarity[j]
        cov_m[i][i] = uncorr[i]**2 + corr[i]**2 + polarity[i]**2
    cov_m = np.array(cov_m)
    val = np.array(val)
    par_vec = np.array(par_vec)
    cov_m_inv = np.linalg.inv(cov_m)
    diff_val_par = val-par_vec
    diff_val_par_T = np.matrix.transpose(diff_val_par)
    tmp = np.matmul(cov_m_inv,diff_val_par)
    chisq = np.matmul(diff_val_par_T,tmp)
    if print_flag:
        print(cov_m)
        print(diff_val_par)
    return chisq

ROOT.gStyle.SetOptStat(0)
cent = [[0,5],[5,10],[10,30],[30,50],[50,90]]

f_names = ["FinalPlot3D_new_2.root","FinalPlot3D_new_2fixmuQ_nopions.root"]
file = [ROOT.TFile(f) for f in f_names]
cvs = [f.Get("cMuBCent") for f in file]
g = [c.FindObject("MuBCent") for c in cvs]
gB = [c.FindObject("gMuBCentPolarityUnc") for c in cvs]
gCorr = [c.FindObject("gMuBCentCorrUnc") for c in cvs]
gCorrPlot = [ROOT.TGraphErrors(g) for g in gCorr]
n = cvs[0].FindObject("NaturePoint")

c = ROOT.TCanvas()
c.cd()
c.SetRightMargin(0.02)
c.SetTopMargin(0.03)
frame = ROOT.TH2D("frame",";Centrality (%);#mu_{#it{B}} (MeV)",1,0,90,1,-0.1,1.3)
frame.Draw()
n.SetFillColor(ROOT.kGray+1)
#n.Draw("samee5")
#gCorr[0].Draw("samepe5")
for p in range(gCorr[1].GetN()):
    if staggering:
        gCorr[1].SetPointX(p,gCorr[1].GetPointX(p)+2)
        g[1].SetPointX(p,g[1].GetPointX(p)+2)
    for gg in g:
        gg.SetPointError(p,0.,gg.GetErrorY(p))
        gg.SetLineWidth(2)
        gg.SetMarkerSize(1.5)
    for polPlot, corrPlot in zip(gB, gCorrPlot):
        corrPlot.SetPointError(p,0.,np.sqrt(corrPlot.GetErrorY(p)**2+polPlot.GetErrorY(p)**2))
InvertPoints(g[0],2,3)
InvertPoints(g[1],2,3)
InvertPoints(gB[0],2,3)
InvertPoints(gB[1],2,3)
InvertPoints(gCorr[0],2,3)
InvertPoints(gCorr[1],2,3)
InvertPoints(gCorrPlot[0],2,3)
InvertPoints(gCorrPlot[1],2,3)
gCorrPlot[1].SetFillColor(ROOT.kBlue)
gCorrPlot[1].SetFillStyle(3345)
#gCorrPlot[1].Draw("samepe3")
gCorrPlot[0].SetFillColor(ROOT.kRed)
gCorrPlot[0].SetFillStyle(3354)
#gCorrPlot[0].Draw("samepe3")
# gB[0].SetFillColor(ROOT.kRed-10)
# gB[1].SetFillColor(ROOT.kBlue-10)
# gB[1].Draw("e3same")
# gB[0].Draw("e3same")
g[0].Draw("pesame")
#InvertPoints(gB[1],2,3)
g[1].SetLineColor(ROOT.kBlue)
g[1].SetMarkerColor(ROOT.kBlue)
g[1].Draw("samepe")

print(f"chisq = {chi2(0.8,g[1],gCorr[1],gB[1],True)}")

cc = ROOT.TCanvas()
cc.cd()
h = ROOT.TH1D("h",";#mu_{B} (MeV);#chi^{2}",150,0,1.5)
for i_p in range(150):
    h.SetBinContent(i_p+1,chi2(i_p*0.01,g[1],gCorr[1],gB[1],False))

h2 = ROOT.TH1D("h2",";#mu_{B} (MeV);#chi^{2}",150,0,1.5)
for i_p in range(150):
    h2.SetBinContent(i_p+1,chi2(i_p*0.01,g[0],gCorr[0],gB[0],False))


c.cd()
l = ROOT.TLegend(0.26817,0.175652,0.452381,0.295217)
l.SetTextFont(44)
l.SetTextSize(22)
l.AddEntry(g[0],"Thermal-FIST, #it{T}_{ch}=155 MeV, #mu_{#it{S}} constrained","pe")
l.AddEntry(g[1],"Thermal-FIST, #it{T}_{ch}=155 MeV, #mu_{#it{S}} and #mu_{#it{Q}} constrained","pe")
#l.AddEntry(n,"SHM fit, #it{Nature} #bf{561}, 321-330 (2018)")
l.Draw("same")

t = ROOT.TLatex()
t.SetTextFont(44)
t.SetTextSize(40)
t.DrawLatex(40.,1.12,"ALICE")
t.SetTextSize(35)
t.DrawLatex(40.,0.95,"Pb-Pb #sqrt{#it{s}_{NN}}=5.02 TeV")
t.SetTextSize(25)
t.DrawLatex(20.,0.25,"#pm0.15 MeV corr. sys. not shown")
#t.DrawLatex(17.,-1.,"0.5 MeV correlated uncertainty not shown")

o = ROOT.TFile("muBvsCent_1.root","recreate")
cc.cd()
h.Fit("pol2")
h2.Fit("pol2")
c.Print("MuBvsCent.pdf")

p0 = h.GetFunction("pol2").GetParameter(0)
p1 = h.GetFunction("pol2").GetParameter(1)
p2 = h.GetFunction("pol2").GetParameter(2)
val = -p1/2./p2
val_err = np.sqrt(p1**2-4*p2*(p1**2/4./p2-1))/2./p2
format_val = "{:.4f}".format(val)
format_val_err = "{:.4f}".format(val_err)
format_chi2 = "{:.2f}".format(h.GetFunction("pol2").Eval(val))
print(f"muB = {format_val} +/- {format_val_err} MeV, chi2 = {format_chi2}/4")

p0 = h2.GetFunction("pol2").GetParameter(0)
p1 = h2.GetFunction("pol2").GetParameter(1)
p2 = h2.GetFunction("pol2").GetParameter(2)
val = -p1/2./p2
val_err = np.sqrt(p1**2-4*p2*(p1**2/4./p2-1))/2./p2
format_val = "{:.4f}".format(val)
format_val_err = "{:.4f}".format(val_err)
format_chi2 = "{:.2f}".format(h2.GetFunction("pol2").Eval(val))
print(f"muB = {format_val} +/- {format_val_err} MeV, chi2 = {format_chi2}/4")

h.Write()
h2.Write()
c.Write()