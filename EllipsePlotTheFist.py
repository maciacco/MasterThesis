import ROOT
import numpy as np

fit_res = []
fit_res_fixMuQ = []
fit_res_fixMuQ_nopions = []
fit_res.append([0.712761, 0.363998, 0.159984, 0.123506, -0.791179]) # chi2 = 1.0596
fit_res_fixMuQ.append([1.09433, 0.100112]) # chi2 = 10.7356
fit_res_fixMuQ_nopions.append([1.09939, 0.100178]) # chi2 = 0.695672
fit_res.append([0.251577, 1.2497, 0.161502, 0.123836, -0.793272]) # chi2 = 5.16327
fit_res_fixMuQ.append([1.5397, 0.099934]) # chi2 = 112.439
fit_res_fixMuQ_nopions.append([1.55914, 0.100142]) # chi2 = 2.02058
fit_res.append([1.17425, -0.527279, 0.158317, 0.123128, -0.788817]) # chi2 = 0.993692
fit_res_fixMuQ.append([0.649667, 0.100177]) # chi2 = 18.9234
fit_res_fixMuQ_nopions.append([0.63977, 0.100189]) # chi2 = 1.45573
fit_res.append([0.724201, 0.299973, 0.150082, 0.111395, -0.761083]) # chi2 = 4.77823
fit_res_fixMuQ.append([1.03859, 0.0995633]) # chi2 = 12.7981
fit_res_fixMuQ_nopions.append([1.04418, 0.0996003]) # chi2 = 4.0161
fit_res.append([0.272802, 1.18598, 0.151581, 0.111655, -0.763454]) # chi2 = 6.79138
fit_res_fixMuQ.append([1.49239, 0.0993922]) # chi2 = 125.647
fit_res_fixMuQ_nopions.append([1.51469, 0.0995805]) # chi2 = 2.74022
fit_res.append([1.17451, -0.591067, 0.1488, 0.111232, -0.759169]) # chi2 = 7.52173
fit_res_fixMuQ.append([0.585884, 0.0996545]) # chi2 = 34.8721
fit_res_fixMuQ_nopions.append([0.576677, 0.0987702]) # chi2 = 5
fit_res.append([0.856205, -0.32537, 0.138354, 0.120137, -0.893478]) # chi2 = 3.71918
fit_res_fixMuQ.append([0.524201, 0.0638392]) # chi2 = 10.4624
fit_res_fixMuQ_nopions.append([0.521517, 0.0638321]) # chi2 = 3.03759
fit_res.append([0.413393, 0.57488, 0.139981, 0.120985, -0.89536]) # chi2 = 2.23355
fit_res_fixMuQ.append([1.01924, 0.0638172]) # chi2 = 25.7678
fit_res_fixMuQ_nopions.append([1.02183, 0.0638299]) # chi2 = 1.56665
fit_res.append([1.29588, -1.22711, 0.136053, 0.118802, -0.89064]) # chi2 = 9
fit_res_fixMuQ.append([0.0283267, 0.0639664]) # chi2 = 115.876
fit_res_fixMuQ_nopions.append([0.0233433, 0.0636091]) # chi2 = 5
fit_res.append([0.828826, 0.140898, 0.09488, 0.0734115, -0.800945]) # chi2 = 3.40207
fit_res_fixMuQ.append([0.984819, 0.0582182]) # chi2 = 7.41971
fit_res_fixMuQ_nopions.append([0.986853, 0.0582288]) # chi2 = 2.90952
fit_res.append([0.333476, 1.04138, 0.0957711, 0.0736439, -0.803095]) # chi2 = 7.02403
fit_res_fixMuQ.append([1.4233, 0.0581529]) # chi2 = 217.522
fit_res_fixMuQ_nopions.append([1.43719, 0.0581324]) # chi2 = 5
fit_res.append([1.3209, -0.762161, 0.0938943, 0.0731459, -0.798547]) # chi2 = 2.33157
fit_res_fixMuQ.append([0.549657, 0.0582763]) # chi2 = 109.282
fit_res_fixMuQ_nopions.append([0.535755, 0.0582398]) # chi2 = 1.03048
fit_res.append([0.79105, -0.343202, 0.205461, 0.179966, -0.903472]) # chi2 = 1.77357
fit_res_fixMuQ.append([0.437958, 0.0904812]) # chi2 = 5.14968
fit_res_fixMuQ_nopions.append([0.435647, 0.0904559]) # chi2 = 1.26367
fit_res.append([0.302007, 0.565742, 0.208151, 0.181488, -0.905428]) # chi2 = 0.440131
fit_res_fixMuQ.append([0.898153, 0.090461]) # chi2 = 10.3317
fit_res_fixMuQ_nopions.append([0.900612, 0.0904335]) # chi2 = 0.292144
fit_res.append([1.27559, -1.25227, 0.203589, 0.179062, -0.902132]) # chi2 = 4.05683
fit_res_fixMuQ.append([-0.0238193, 0.0906695]) # chi2 = 52.7022
fit_res_fixMuQ_nopions.append([-0.0300183, 0.0904955]) # chi2 = 2.92666

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

out_file = ROOT.TFile("EllipsePlot.root","recreate")
canv = ROOT.TCanvas("canv","canv",600,600)
h_frame = ROOT.TH2F("h_frame"," ",1,.1,1.5,1,-1.6,1.6)
h_frame.GetXaxis().SetTitle("#mu_{#it{B}} (MeV)")
h_frame.GetYaxis().SetTitle("#mu_{#it{Q}} (MeV)")
h_frame.GetXaxis().SetNdivisions(10)
h_frame.GetYaxis().SetNdivisions(10)

h_frame.Draw()
ellipses = []
ellipses_corr = []
points = []

for i_cent, cent in enumerate(CENTRALITY_CLASSES):
    muB = fit_res[int(3*i_cent)][0]
    muB_err_corr = np.abs(fit_res[int(3*i_cent+2)][0]-fit_res[int(3*i_cent+1)][0])*0.5
    muI = fit_res[int(3*i_cent)][1]
    muI_err_corr = np.abs(fit_res[int(3*i_cent+2)][1]-fit_res[int(3*i_cent+1)][1])*0.5
    corr = fit_res[int(3*i_cent)][4]
    ellipses_corr.append(ROOT.RooEllipse(f"corrEllipse_{cent[0]}_{cent[1]}",muB,muI,muB_err_corr,muI_err_corr,corr))
    ellipses_corr[i_cent].SetLineColor(centrality_colors[i_cent])
    ellipses_corr[i_cent].SetLineStyle(ROOT.kDashed)
    ellipses_corr[i_cent].SetLineWidth(2)
    ellipses_corr[i_cent].Draw("lsame")

ellipses_corr.append(ROOT.RooEllipse("dummy_corr",-100,-100,1,1,0))
ellipses_corr[5].SetLineStyle(ROOT.kDashed)
ellipses_corr[5].SetLineWidth(2)
ellipses_corr[5].SetLineColor(ROOT.kBlack)

for i_cent, cent in enumerate(CENTRALITY_CLASSES):
    muB = fit_res[int(3*i_cent)][0]
    muB_err = fit_res[int(3*i_cent)][2]
    muI = fit_res[int(3*i_cent)][1]
    muI_err = fit_res[int(3*i_cent)][3]
    corr = fit_res[int(3*i_cent)][4]
    ellipses.append(ROOT.RooEllipse(f"Ellipse_{cent[0]}_{cent[1]}",muB,muI,muB_err,muI_err,corr))
    ellipses[i_cent].SetLineColor(centrality_colors[i_cent])
    ellipses[i_cent].SetLineWidth(2)
    print(f"muB = {muB}, errmuB = {muB_err}, muQ = {muI}, errmuQ = {muI_err}, corr = {corr}")
    # ellipses[i_cent].SetFillStyle(3345)
    # ellipses[i_cent].SetFillColor(centrality_colors_area[i_cent])

    ellipses[i_cent].Draw("lsame")
    points.append(ROOT.TGraph())
    points[i_cent].AddPoint(muB,muI)
    points[i_cent].SetMarkerColor(centrality_colors[i_cent])
    points[i_cent].SetMarkerStyle(20)
    points[i_cent].SetMarkerSize(1.2)
    points[i_cent].Draw("psame")

ellipses.append(ROOT.RooEllipse("dummy",-100,-100,1,1,0))
ellipses[5].SetLineWidth(2)
#ellipses[5].SetFillColor(ROOT.kBlack)
ellipses[5].SetLineColor(ROOT.kBlack)

# plot Pb208 line
# f_Pb208 = ROOT.TF1("f_Pb208","[0]*x",0,100)
# f_Pb208.SetParameter(0,Z_PB208/A_PB208)
# f_Pb208.SetLineColor(ROOT.kRed-3)
# f_Pb208.Draw("same")

# add labels
text = ROOT.TLatex()
text.SetTextFont(43)
text.SetTextSize(32)
text.DrawLatex(1.20355,1.19636,"ALICE")
text.SetTextSize(28)
text.DrawLatex(0.740964,0.947157,"Pb-Pb #sqrt{s_{NN}}=5.02 TeV")

# legend
leg = ROOT.TLegend(0.152174,0.292969,0.307692,0.515625)
leg.SetTextFont(43)
leg.SetTextSize(22)
leg.AddEntry(points[0],"0-5%","p")
leg.AddEntry(points[1],"5-10%","p")
leg.AddEntry(points[3],"10-30%","p")
leg.AddEntry(points[2],"30-50%","p")
leg.AddEntry(points[4],"50-90%","p")
leg.Draw("same")

legErr = ROOT.TLegend(0.152174,0.164211,0.342809,0.277895)
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
canv.Print("ellipsePlot_fist.pdf")
out_file.Close()