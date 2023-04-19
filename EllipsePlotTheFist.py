import ROOT
import numpy as np

chi2 = []
chi2_fixmuQ = []
chi2_fixmuQ_nopions = []
fit_res = []
fit_res_fixMuQ = []
fit_res_fixMuQ_nopions = []

fit_res.append([0.646566, 0.459121, 0.160271, 0.123529, -0.792434])
chi2.append(0.970221)
fit_res_fixMuQ.append([1.12586, 0.0999647, -0.0283059])
chi2_fixmuQ.append(15.9348)
fit_res_fixMuQ_nopions.append([1.13196, 0.100004, -0.0284593])
chi2_fixmuQ_nopions.append(0.515368)
fit_res.append([0.179215, 1.3487, 0.161486, 0.123803, -0.794266])
chi2.append(3.95279)
fit_res_fixMuQ.append([1.57088, 0.099782, -0.0394945])
chi2_fixmuQ.append(128.816)
fit_res_fixMuQ_nopions.append([1.59169, 0.0999803, -0.0400177])
chi2_fixmuQ_nopions.append(1.89386)
fit_res.append([1.11276, -0.434212, 0.158374, 0.122998, -0.789521])
chi2.append(1.10316)
fit_res_fixMuQ.append([0.680615, 0.100026, -0.0171117])
chi2_fixmuQ.append(13.1512)
fit_res_fixMuQ_nopions.append([0.672583, 0.100017, -0.0169098])
chi2_fixmuQ_nopions.append(1.26066)
fit_res.append([0.50497, 0.824122, 0.161033, 0.123783, -0.793697])
chi2.append(0.846032)
fit_res_fixMuQ.append([1.36011, 0.099905, -0.0341952])
chi2_fixmuQ.append(47.9625)
fit_res_fixMuQ_nopions.append([1.37137, 0.100009, -0.0344785])
chi2_fixmuQ_nopions.append(0.165534)
fit_res.append([0.787958, 0.0935532, 0.159421, 0.123239, -0.791019])
chi2.append(1.28637)
fit_res_fixMuQ.append([0.890967, 0.0999835, -0.0224003])
chi2_fixmuQ.append(2.00803)
fit_res_fixMuQ_nopions.append([0.892611, 0.100003, -0.0224416])
chi2_fixmuQ_nopions.append(1.13344)

fit_res.append([0.781426, 0.172088, 0.129462, 0.120502, -0.575297])
chi2.append(0.132713)
fit_res_fixMuQ.append([1.04704, 0.0997645, -0.0263242])
chi2_fixmuQ.append(15.8296)
fit_res_fixMuQ_nopions.append([1.05344, 0.0997956, -0.026485])
chi2_fixmuQ_nopions.append(3.96268)
fit_res.append([0.20045, 1.2615, 0.154049, 0.114629, -0.771964])
chi2.append(4.91523)
fit_res_fixMuQ.append([1.50106, 0.0995939, -0.0377391])
chi2_fixmuQ.append(132.243)
fit_res_fixMuQ_nopions.append([1.5235, 0.0997715, -0.0383032])
chi2_fixmuQ_nopions.append(2.67027)
fit_res.append([1.12138, -0.527379, 0.15115, 0.113992, -0.767127])
chi2.append(7.56368)
fit_res_fixMuQ.append([0.593652, 0.0998487, -0.0149253])
chi2_fixmuQ.append(28.3695)
fit_res_fixMuQ_nopions.append([0.585925, 0.0989956, -0.014731])
chi2_fixmuQ_nopions.append(5)
fit_res.append([0.520498, 0.733978, 0.153444, 0.114551, -0.771052])
chi2.append(4.90904)
fit_res_fixMuQ.append([1.28119, 0.0997111, -0.0322112])
chi2_fixmuQ.append(48.6687)
fit_res_fixMuQ_nopions.append([1.29367, 0.099786, -0.032525])
chi2_fixmuQ_nopions.append(4.07446)
fit_res.append([0.80296, 0.00311398, 0.152021, 0.114138, -0.768449])
chi2.append(4.1706)
fit_res_fixMuQ.append([0.812729, 0.0997815, -0.0204332])
chi2_fixmuQ.append(4.11037)
fit_res_fixMuQ_nopions.append([0.813155, 0.0998048, -0.0204439])
chi2_fixmuQ_nopions.append(4.07014)

fit_res.append([0.697399, -0.178802, 0.141653, 0.118879, -0.864176])
chi2.append(3.39079)
fit_res_fixMuQ.append([0.51438, 0.0732017, -0.0129323])
chi2_fixmuQ.append(5.36281)
fit_res_fixMuQ_nopions.append([0.512595, 0.0732049, -0.0128874])
chi2_fixmuQ_nopions.append(3.13772)
fit_res.append([0.256138, 0.717743, 0.143405, 0.119689, -0.866585])
chi2.append(2.2324)
fit_res_fixMuQ.append([1.00781, 0.0731757, -0.0253379])
chi2_fixmuQ.append(39.7421)
fit_res_fixMuQ_nopions.append([1.01276, 0.0731883, -0.0254624])
chi2_fixmuQ_nopions.append(1.54537)
fit_res.append([1.13366, -1.0757, 0.139298, 0.117667, -0.860728])
chi2.append(9)
fit_res_fixMuQ.append([0.0208937, 0.0733509, -0.000525299])
chi2_fixmuQ.append(92.8789)
fit_res_fixMuQ_nopions.append([0.0142132, 0.0728957, -0.00035734])
chi2_fixmuQ_nopions.append(5)
fit_res.append([0.697399, -0.178802, 0.141653, 0.118879, -0.864176])
chi2.append(3.39079)
fit_res_fixMuQ.append([0.51438, 0.0732017, -0.0129323])
chi2_fixmuQ.append(5.36281)
fit_res_fixMuQ_nopions.append([0.512595, 0.0732049, -0.0128874])
chi2_fixmuQ_nopions.append(3.13772)
fit_res.append([0.697399, -0.178802, 0.141653, 0.118879, -0.864176])
chi2.append(3.39079)
fit_res_fixMuQ.append([0.51438, 0.0732017, -0.0129323])
chi2_fixmuQ.append(5.36281)
fit_res_fixMuQ_nopions.append([0.512595, 0.0732049, -0.0128874])
chi2_fixmuQ_nopions.append(3.13772)

fit_res.append([0.800627, 0.160439, 0.0958831, 0.0743065, -0.80318])
chi2.append(2.93535)
fit_res_fixMuQ.append([0.975953, 0.0585538, -0.024537])
chi2_fixmuQ.append(8.5143)
fit_res_fixMuQ_nopions.append([0.978234, 0.0585657, -0.0245943])
chi2_fixmuQ_nopions.append(2.98952)
fit_res.append([0.304551, 1.05982, 0.0968229, 0.074562, -0.805448])
chi2.append(6.38844)
fit_res_fixMuQ.append([1.41439, 0.0584859, -0.03556])
chi2_fixmuQ.append(219.35)
fit_res_fixMuQ_nopions.append([1.42864, 0.0584718, -0.0359184])
chi2_fixmuQ_nopions.append(5)
fit_res.append([1.29288, -0.741484, 0.0947708, 0.0739918, -0.800457])
chi2.append(1.91356)
fit_res_fixMuQ.append([0.540445, 0.058606, -0.0135876])
chi2_fixmuQ.append(101.153)
fit_res_fixMuQ_nopions.append([0.527268, 0.0585711, -0.0132563])
chi2_fixmuQ_nopions.append(1.09386)
fit_res.append([0.705583, 0.404583, 0.0962231, 0.0744291, -0.804094])
chi2.append(2.57278)
fit_res_fixMuQ.append([1.13538, 0.0585435, -0.0285453])
chi2_fixmuQ.append(34.9773)
fit_res_fixMuQ_nopions.append([1.1406, 0.0585632, -0.0286765])
chi2_fixmuQ_nopions.append(2.79656)
fit_res.append([0.895439, -0.0838715, 0.0955375, 0.0741821, -0.802246])
chi2.append(3.36936)
fit_res_fixMuQ.append([0.817058, 0.0585581, -0.0205421])
chi2_fixmuQ.append(4.23893)
fit_res_fixMuQ_nopions.append([0.815846, 0.0585684, -0.0205116])
chi2_fixmuQ_nopions.append(3.23815)

fit_res.append([0.649342, -0.208717, 0.19268, 0.164762, -0.883546])
chi2.append(1.59046)
fit_res_fixMuQ.append([0.432747, 0.0926938, -0.0108799])
chi2_fixmuQ.append(3.04905)
fit_res_fixMuQ_nopions.append([0.430914, 0.0926928, -0.0108338])
chi2_fixmuQ_nopions.append(1.34692)
fit_res.append([0.163009, 0.694674, 0.195181, 0.166011, -0.885776])
chi2.append(0.272443)
fit_res_fixMuQ.append([0.892224, 0.0926737, -0.0224319])
chi2_fixmuQ.append(18.4048)
fit_res_fixMuQ_nopions.append([0.895879, 0.0926761, -0.0225238])
chi2_fixmuQ_nopions.append(0.416871)
fit_res.append([1.13111, -1.11327, 0.190841, 0.16404, -0.881969])
chi2.append(3.97182)
fit_res_fixMuQ.append([-0.0271187, 0.0928731, 0.000681804])
chi2_fixmuQ.append(49.7348)
fit_res_fixMuQ_nopions.append([-0.034709, 0.0927288, 0.000872637])
chi2_fixmuQ_nopions.append(2.9769)
fit_res.append([0.649342, -0.208717, 0.19268, 0.164762, -0.883546])
chi2.append(1.59046)
fit_res_fixMuQ.append([0.432747, 0.0926938, -0.0108799])
chi2_fixmuQ.append(3.04905)
fit_res_fixMuQ_nopions.append([0.430914, 0.0926928, -0.0108338])
chi2_fixmuQ_nopions.append(1.34692)
fit_res.append([0.649342, -0.208717, 0.19268, 0.164762, -0.883546])
chi2.append(1.59046)
fit_res_fixMuQ.append([0.432747, 0.0926938, -0.0108799])
chi2_fixmuQ.append(3.04905)
fit_res_fixMuQ_nopions.append([0.430914, 0.0926928, -0.0108338])
chi2_fixmuQ_nopions.append(1.34692)


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
h_frame = ROOT.TH2F("h_frame"," ",1,.1,1.35,1,-1.35,1.6)
h_frame.GetXaxis().SetTitle("#mu_{#it{B}} (MeV)")
h_frame.GetYaxis().SetTitle("#mu_{#it{Q}} (MeV)")
h_frame.GetXaxis().SetNdivisions(10)
h_frame.GetYaxis().SetNdivisions(10)

h_frame.Draw()
ellipses = []
ellipses_corr = []
ellipses_corr_ = []
ellipses_polarity = []
points = []
ROOT.gStyle.SetEndErrorSize(5)

for i_cent, cent in enumerate(CENTRALITY_CLASSES):
    muB = fit_res[int(5*i_cent)][0]
    muB_err_corr = np.abs(fit_res[int(5*i_cent+2)][0]-fit_res[int(5*i_cent+1)][0])*0.5
    muI = fit_res[int(5*i_cent)][1]
    muI_err_corr = np.abs(fit_res[int(5*i_cent+2)][1]-fit_res[int(5*i_cent+1)][1])*0.5
    corr_ = fit_res[int(5*i_cent)][4]
    muB_err_polarity = np.abs(fit_res[int(5*i_cent+4)][0]-fit_res[int(5*i_cent+3)][0])*0.5
    muI_err_polarity = np.abs(fit_res[int(5*i_cent+4)][1]-fit_res[int(5*i_cent+3)][1])*0.5
    #corr = -(muB_err_corr*muI_err_corr+muB_err_polarity*muI_err_polarity)/np.sqrt((muB_err_corr**2+muB_err_polarity**2)*(muI_err_corr**2+muI_err_polarity**2))
    #ellipses_corr.append(ROOT.RooEllipse(f"corrEllipse_{cent[0]}_{cent[1]}",muB,muI,np.sqrt(muB_err_corr**2+muB_err_polarity**2),np.sqrt(muI_err_corr**2+muI_err_polarity**2),corr,100000))
    ellipses_corr_.append(ROOT.TGraphErrors(1,np.array([muB]),np.array([muI]),np.array([np.sqrt(muB_err_corr**2+muB_err_polarity**2)]),np.array([np.sqrt(muI_err_corr**2+muI_err_polarity**2)])))
    ellipses_corr_[i_cent].SetLineColor(centrality_colors[i_cent])
    #ellipses_corr_[i_cent].SetLineStyle(ROOT.kDashed)
    ellipses_corr_[i_cent].SetFillColor(centrality_colors[i_cent])
    ellipses_corr_[i_cent].SetFillStyle(3345)
    ellipses_corr_[i_cent].SetLineWidth(2)
    ellipses_corr_[i_cent].Draw("e[]same")
    #print(f"muB = {muB}, corr = 1")

for i_cent, cent in enumerate(CENTRALITY_CLASSES):
    muB = fit_res[int(5*i_cent)][0]
    muB_err_corr = np.abs(fit_res[int(5*i_cent+2)][0]-fit_res[int(5*i_cent+1)][0])*0.5
    muI = fit_res[int(5*i_cent)][1]
    muI_err_corr = np.abs(fit_res[int(5*i_cent+2)][1]-fit_res[int(5*i_cent+1)][1])*0.5
    corr_ = fit_res[int(5*i_cent)][4]
    muB_err_polarity = np.abs(fit_res[int(5*i_cent+4)][0]-fit_res[int(5*i_cent+3)][0])*0.5
    muI_err_polarity = np.abs(fit_res[int(5*i_cent+4)][1]-fit_res[int(5*i_cent+3)][1])*0.5
    #corr = -(muB_err_corr*muI_err_corr+muB_err_polarity*muI_err_polarity)/np.sqrt((muB_err_corr**2+muB_err_polarity**2)*(muI_err_corr**2+muI_err_polarity**2))
    #ellipses_corr.append(ROOT.RooEllipse(f"corrEllipse_{cent[0]}_{cent[1]}",muB,muI,np.sqrt(muB_err_corr**2+muB_err_polarity**2),np.sqrt(muI_err_corr**2+muI_err_polarity**2),corr,100000))
    ellipses_corr.append(ROOT.TGraphErrors(1,np.array([muB]),np.array([muI]),np.array([np.sqrt(muB_err_corr**2+muB_err_polarity**2)]),np.array([np.sqrt(muI_err_corr**2+muI_err_polarity**2)])))
    ellipses_corr[i_cent].SetLineColor(centrality_colors[i_cent])
    ellipses_corr[i_cent].SetLineStyle(4)
    ellipses_corr[i_cent].SetFillColor(centrality_colors[i_cent])
    ellipses_corr[i_cent].SetFillStyle(3345)
    ellipses_corr[i_cent].SetLineWidth(2)
    ellipses_corr[i_cent].Draw("esame")
    print(f"muB = {muB}, corr = 1")

ellipses_corr.append(ROOT.RooEllipse("dummy_corr",-100,-100,1,1,0))
ellipses_corr[5].SetLineStyle(4)
ellipses_corr[5].SetLineWidth(2)
ellipses_corr[5].SetLineColor(ROOT.kBlack)

# for i_cent, cent in enumerate(CENTRALITY_CLASSES):
#     muB = fit_res[int(5*i_cent)][0]
#     muB_err_polarity = np.abs(fit_res[int(5*i_cent+4)][0]-fit_res[int(5*i_cent+3)][0])*0.5
#     muI = fit_res[int(5*i_cent)][1]
#     muI_err_polarity = np.abs(fit_res[int(5*i_cent+4)][1]-fit_res[int(5*i_cent+3)][1])*0.5
#     polarity = fit_res[int(5*i_cent)][4]
#     ellipses_polarity.append(ROOT.RooEllipse(f"polarityEllipse_{cent[0]}_{cent[1]}",muB,muI,muB_err_polarity,muI_err_polarity,polarity))
#     ellipses_polarity[i_cent].SetLineColor(centrality_colors[i_cent])
#     ellipses_polarity[i_cent].SetLineStyle(ROOT.kDotted)
#     ellipses_polarity[i_cent].SetLineWidth(2)
#     ellipses_polarity[i_cent].Draw("lsame")

# ellipses_polarity.append(ROOT.RooEllipse("dummy_polarity",-100,-100,1,1,0))
# ellipses_polarity[5].SetLineStyle(ROOT.kDotted)
# ellipses_polarity[5].SetLineWidth(2)
# ellipses_polarity[5].SetLineColor(ROOT.kBlack)

for i_cent, cent in enumerate(CENTRALITY_CLASSES):
    muB = fit_res[int(5*i_cent)][0]
    muB_err = fit_res[int(5*i_cent)][2]
    muI = fit_res[int(5*i_cent)][1]
    muI_err = fit_res[int(5*i_cent)][3]
    corr = fit_res[int(5*i_cent)][4]
    ellipses.append(ROOT.RooEllipse(f"Ellipse_{cent[0]}_{cent[1]}",muB,muI,muB_err,muI_err,corr,100000))
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
text.SetTextSize(29)
text.DrawLatex(0.18,1.3,"ALICE")
text.SetTextSize(24)
text.DrawLatex(0.18,1.1,"Pb-Pb")
text.DrawLatex(0.18,0.9,"#sqrt{s_{NN}}=5.02 TeV")
text.DrawLatex(0.9,1.3,"Thermal-FIST")
text.DrawLatex(0.9,1.1,"#it{T}_{ch}=155#pm2 MeV")
text.DrawLatex(0.9,0.9,"#mu_{#it{S}} constrained")
#text.SetTextSize(28)
#text.DrawLatex(1.19,0.69795400,"|y| < 0.5")

# legend
#leg = ROOT.TLegend(0.152174,0.292969,0.307692,0.515625)
leg = ROOT.TLegend(0.697324,0.158261,0.851171,0.38087)
leg.SetTextFont(43)
leg.SetTextSize(22)
leg.AddEntry(points[0],"0-5%","p")
leg.AddEntry(points[1],"5-10%","p")
leg.AddEntry(points[3],"10-30%","p")
leg.AddEntry(points[2],"30-50%","p")
leg.AddEntry(points[4],"50-90%","p")
leg.Draw("same")

legErr = ROOT.TLegend(0.152174+0.03,0.164211,0.342809+0.05,0.277895)
legErr.SetNColumns(1)
legErr.SetTextFont(43)
legErr.SetTextSize(22)
legErr.AddEntry(ellipses[5],"Uncorr. uncert.", "f")
legErr.AddEntry(ellipses_corr[5],"Corr. uncert.", "l")
#legErr.AddEntry(ellipses_polarity[5],"B field uncert.", "f")
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
canv.SetRightMargin(0.02)
canv.SetLeftMargin(0.13)
canv.SetTopMargin(0.02)
canv.Write()
canv.Print("ellipsePlot_fist.pdf")
out_file.Close()