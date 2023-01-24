import ROOT
import numpy as np

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

def chi2(par,g,gCorr,gB,print_cov=False):
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
    if print_cov:
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
    if print_cov:
        print(cov_m)
    return chisq

staggering = False

ROOT.gStyle.SetOptStat(0)
cent = [[0,5],[5,10],[10,30],[30,50],[50,90]]

f_names = ["FinalPlot3D_new.root","FinalPlot3D_new_fixmuQ_nopions.root"]
file = [ROOT.TFile(f) for f in f_names]
cvs = [f.Get("cMuQuncorr") for f in file]
g = [c.FindObject("MuQCent") for c in cvs]
gB = cvs[0].FindObject("MuQCentPolarity")
gCorr = cvs[0].FindObject("MuQCentCorrUnc")
gCorrPlot = ROOT.TGraphErrors(gCorr)

c = ROOT.TCanvas()
c.cd()
c.SetRightMargin(0.02)
c.SetTopMargin(0.03)
frame = ROOT.TH2D("frame",";Centrality (%);#mu_{#it{Q}} (MeV)",1,0,90,1,-2.,1.8)
frame.Draw()
#gCorr[0].Draw("samepe5")
for p in range(g[1].GetN()):
    if staggering:
        g[1].SetPointX(p,g[1].GetPointX(p)+2)
    for gg in g:
        gg.SetPointError(p,0.,gg.GetErrorY(p))
        gg.SetLineWidth(2)
        gg.SetMarkerSize(1.5)
        gCorr.SetPointError(p,0.,gCorr.GetErrorY(p))
    gCorrPlot.SetPointError(p,0.,np.sqrt(gCorrPlot.GetErrorY(p)**2+gB.GetErrorY(p)**2))
InvertPoints(g[0],2,3)
InvertPoints(g[1],2,3)
InvertPoints(gCorr,2,3)
InvertPoints(gB,2,3)
InvertPoints(gCorrPlot,2,3)
gCorrPlot.SetFillColor(ROOT.kRed)
gCorrPlot.SetFillStyle(3354)
gCorrPlot.Draw("e3")
# gB.SetLineColor(ROOT.kRed)
# gB.SetFillColor(ROOT.kRed-10)
# gB.Draw("e3same")
g[0].Draw("pesame")
g[0].GetFunction("pol0").Delete()
g[1].SetLineColor(ROOT.kBlue)
g[1].SetFillColor(ROOT.kBlue)
g[1].SetMarkerColor(ROOT.kBlue)
g[1].Draw("e3same")
g[1].GetFunction("pol0").Delete()

print(f"chisq = {chi2(0.1,g[0],gCorr,gB,True)}")

cc = ROOT.TCanvas()
cc.cd()
h = ROOT.TH1D("h",";#mu_{Q} (MeV);#chi^{2}",300,-1.5,1.5)
for i_p in range(300):
    h.SetBinContent(i_p+1,chi2(-1.5+i_p*0.01,g[0],gCorr,gB))

c.cd()
l = ROOT.TLegend(0.26817,0.175652,0.452381,0.295217)
l.SetTextFont(44)
l.SetTextSize(22)
l.AddEntry(g[0],"Thermal-FIST, #it{T}_{ch}=155 MeV, #mu_{#it{S}} constrained","pe")
l.AddEntry(g[1],"Thermal-FIST, #it{T}_{ch}=155 MeV, #mu_{#it{S}} and #mu_{#it{Q}} constrained","l")
l.Draw("same")

t = ROOT.TLatex()
t.SetTextFont(44)
t.SetTextSize(40)
t.DrawLatex(40.,1.3,"ALICE")
t.SetTextSize(35)
t.DrawLatex(40.,.9,"Pb-Pb #sqrt{#it{s}_{NN}}=5.02 TeV")
t.SetTextSize(25)
#t.DrawLatex(15.7949,-0.646823,"0.9 MeV correlated uncertainty not shown")

o = ROOT.TFile("muQvsCent.root","recreate")
c.Write()
cc.cd()
h.Fit("pol2")
c.Print("MuQvsCent.pdf")
p0 = h.GetFunction("pol2").GetParameter(0)
p1 = h.GetFunction("pol2").GetParameter(1)
p2 = h.GetFunction("pol2").GetParameter(2)
val = -p1/2./p2
val_err = np.sqrt(p1**2-4*p2*(p1**2/4./p2-1))/2./p2
format_val = "{:.2f}".format(val)
format_val_err = "{:.2f}".format(val_err)
format_chi2 = "{:.2f}".format(h.GetFunction("pol2").Eval(val))
print(f"muQ = {format_val} +/- {format_val_err} MeV, chi2 = {format_chi2}/4")
h.Write()