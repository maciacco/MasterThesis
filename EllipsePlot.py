import ROOT
import numpy as np


CENTRALITY_CLASSES = [[0,5], [5,10], [30,50], [10, 30], [50,90]]
centrality_colors = [ROOT.kRed, ROOT.kOrange-3, ROOT.kAzure+4,ROOT.kGreen+2,ROOT.kMagenta+2]
centrality_colors_area = [ROOT.kRed, ROOT.kOrange-3, ROOT.kAzure+4,ROOT.kGreen+2,ROOT.kMagenta+2]

Z_PB208 = 82
A_PB208 = 208

# compute 2-dimensional chi2 with correlations
def chi2_2d(in_file):
    in_file = ROOT.TFile("FinalPlot3D.root")
    cMuB = in_file.Get("cMuBCent")
    cMuI = in_file.Get("cMuICent")
    gMuB = cMuB.GetPrimitive("MuBCent")
    gMuI = cMuI.GetPrimitive("MuICent")
    gMuBCorr = cMuB.GetPrimitive("gMuBCentCorrUnc")
    gMuICorr = cMuI.GetPrimitive("gMuICentCorrUnc")

    chi2 = 0
    for i_cent, cent in enumerate(CENTRALITY_CLASSES):
        muB = gMuB.GetPointY(i_cent)
        muB_err = gMuB.GetErrorY(i_cent)
        muI = gMuI.GetPointY(i_cent)
        muI_err = gMuI.GetErrorY(i_cent)
        corr_mat = in_file.Get(f"CovMat_{cent[0]}_{cent[1]}")
        corr = corr_mat.GetBinContent(2,1)
        cov = corr*muB_err*muI_err
        A = Z_PB208/A_PB208
        distance = np.abs(A*muB-muI)/np.sqrt(1+A**2)
        #var = (((muB-(A*muI+muB)/(A**2+1))/distance)**2)*(muB_err**2)+(((muI-A*(A*muI+muB)/(A**2+1))/distance)**2)*(muI_err**2)+2*((muI-A*(A*muI+muB)/(A**2+1))*(muB-(A*muI+muB)/(A**2+1))/distance/distance)*cov
        pre_factor = 1/(A**2+1)/2/np.sqrt((A*muI-A*A*muB)**2+(A*muB-muI)**2)
        dLdx = pre_factor*(-2*(A*muI-A*A*muB)*A*A+2*A*(A*muB-muI))
        dLdy = pre_factor*(2*A*(A*muI-A*A*muB)-2*(A*muB-muI))
        var = dLdx*dLdx*muB_err*muB_err+dLdy*dLdy*muI_err*muI_err+2*dLdx*dLdy*cov
        print(f"cent = {cent[0]}-{cent[1]}%, delta = {distance}, error = {np.sqrt(var)}")
        print(f"muB = {muB}, errmuB = {muB_err}, muQ = {muI}, errmuQ = {muI_err}, cov = {cov}")
        chi2 = chi2 + distance*distance/var
    return chi2


ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)


in_file = ROOT.TFile("FinalPlot3D.root")
cMuB = in_file.Get("cMuBCent")
cMuI = in_file.Get("cMuICent")

out_file = ROOT.TFile("EllipsePlot.root","recreate")
canv = ROOT.TCanvas("canv","canv",600,500)
h_frame = ROOT.TH2F("h_frame"," ",1,0.,1.3,1,-0.45,.55)
h_frame.GetXaxis().SetTitle("#mu_{#it{B}} (MeV)")
h_frame.GetYaxis().SetTitle("#mu_{#it{Q}} (MeV)")
h_frame.GetXaxis().SetNdivisions(10)
h_frame.GetYaxis().SetNdivisions(5)

h_frame.Draw()
ellipses = []
ellipses_corr = []
points = []

for i_cent, cent in enumerate(CENTRALITY_CLASSES):
    gMuB = cMuB.GetPrimitive("MuBCent")
    gMuI = cMuI.GetPrimitive("MuICent")
    gMuBCorr = cMuB.GetPrimitive("gMuBCentCorrUnc")
    gMuICorr = cMuI.GetPrimitive("gMuICentCorrUnc")
    muB = gMuB.GetPointY(i_cent)
    muB_err_corr = gMuBCorr.GetErrorY(i_cent)
    muI = gMuI.GetPointY(i_cent)
    muI_err_corr = gMuICorr.GetErrorY(i_cent)
    corr_mat = in_file.Get(f"CovMat_{cent[0]}_{cent[1]}")
    corr = corr_mat.GetBinContent(2,1)
    ellipses_corr.append(ROOT.RooEllipse(f"corrEllipse_{cent[0]}_{cent[1]}",muB,muI,muB_err_corr,muI_err_corr,corr))
    ellipses_corr[i_cent].SetLineColor(centrality_colors[i_cent])
    ellipses_corr[i_cent].SetFillStyle(3345)
    ellipses_corr[i_cent].SetFillColor(centrality_colors_area[i_cent])
    ellipses_corr[i_cent].Draw("lfsame")

ellipses_corr.append(ROOT.RooEllipse("dummy_corr",-100,-100,1,1,0))
ellipses_corr[5].SetFillStyle(3345)
ellipses_corr[5].SetFillColor(ROOT.kBlack)
ellipses_corr[5].SetLineColor(ROOT.kBlack)

for i_cent, cent in enumerate(CENTRALITY_CLASSES):
    gMuB = cMuB.GetPrimitive("MuBCent")
    gMuI = cMuI.GetPrimitive("MuICent")
    muB = gMuB.GetPointY(i_cent)
    muB_err = gMuB.GetErrorY(i_cent)
    muI = gMuI.GetPointY(i_cent)
    muI_err = gMuI.GetErrorY(i_cent)
    corr_mat = in_file.Get(f"CovMat_{cent[0]}_{cent[1]}")
    corr = corr_mat.GetBinContent(2,1)
    ellipses.append(ROOT.RooEllipse(f"Ellipse_{cent[0]}_{cent[1]}",muB,muI,muB_err,muI_err,corr))
    ellipses[i_cent].SetLineColor(centrality_colors[i_cent])

    ellipses[i_cent].Draw("lsame")
    points.append(ROOT.TGraph())
    points[i_cent].AddPoint(muB,muI)
    points[i_cent].SetMarkerColor(centrality_colors[i_cent])
    points[i_cent].SetMarkerStyle(20)
    points[i_cent].SetMarkerSize(1)
    points[i_cent].Draw("psame")

ellipses.append(ROOT.RooEllipse("dummy",-100,-100,1,1,0))
ellipses[5].SetLineColor(ROOT.kBlack)

# plot Pb208 line
# f_Pb208 = ROOT.TF1("f_Pb208","[0]*x",0,100)
# f_Pb208.SetParameter(0,Z_PB208/A_PB208)
# f_Pb208.SetLineColor(ROOT.kRed-3)
# f_Pb208.Draw("same")

# add labels
text = ROOT.TLatex()
text.SetTextFont(43)
text.SetTextSize(28)
text.DrawLatex(0.0634404,0.454049 ,"ALICE")
text.SetTextSize(24)
text.DrawLatex(0.0634404,0.378745 ,"Pb-Pb #sqrt{s_{NN}}=5.02 TeV")

# legend
leg = ROOT.TLegend(0.168896,0.423158,0.240803,0.715789)
leg.SetTextFont(43)
leg.SetTextSize(22)
leg.AddEntry(points[0],"0-5%","p")
leg.AddEntry(points[1],"5-10%","p")
leg.AddEntry(points[3],"10-30%","p")
leg.AddEntry(points[2],"30-50%","p")
leg.AddEntry(points[4],"50-90%","p")
leg.Draw("same")

legErr = ROOT.TLegend(0.153846,0.151579,0.342809,0.265263)
legErr.SetNColumns(1)
legErr.SetTextFont(43)
legErr.SetTextSize(22)
legErr.AddEntry(ellipses[5],"Uncorr. uncert.", "f")
legErr.AddEntry(ellipses_corr[5],"Corr. uncert.", "f")
legErr.Draw("same")

# chi2 = chi2_2d(in_file)
# print(f"Chi2 = {chi2}")

# canv.cd()
# legLine = ROOT.TLegend(0.637124,0.176842,0.924749,0.231579)
# legLine.SetNColumns(1)
# legLine.SetTextFont(43)
# legLine.SetTextSize(22)
# legLine.AddEntry(f_Pb208,"#frac{#mu_{#it{Q}}}{#mu_{#it{B}}} = #frac{Z}{A}({}^{208}Pb)","l")
# legLine.Draw("same")

#format_chi2 = "{:.1f}".format(chi2)
#text.DrawLatex(0.737978,0.501619,"#chi^{2}_{TLS}/ndf = "+format_chi2+"/3")

out_file.cd()
canv.Write()
canv.Print("ellipsePlot.pdf")
out_file.Close()