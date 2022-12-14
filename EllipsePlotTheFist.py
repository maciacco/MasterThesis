import ROOT
import numpy as np

chi2 = []
chi2_fixmuQ = []
chi2_fixmuQ_nopions = []
fit_res = []
fit_res_fixMuQ = []
fit_res_fixMuQ_nopions = []

fit_res.append([0.708006, 0.368355, 0.160412, 0.123792, -0.792493])
chi2.append(0.962495)
fit_res_fixMuQ.append([1.09442, 0.100108])
chi2_fixmuQ.append(10.7233)
fit_res_fixMuQ_nopions.append([1.0995, 0.100145])
chi2_fixmuQ_nopions.append(0.60714)
fit_res.append([0.239212, 1.25969, 0.161754, 0.124116, -0.794516])
chi2.append(3.90689)
fit_res_fixMuQ.append([1.53987, 0.0999371])
chi2_fixmuQ.append(112.435)
fit_res_fixMuQ_nopions.append([1.55907, 0.100116])
chi2_fixmuQ_nopions.append(1.94999)
fit_res.append([1.17473, -0.526276, 0.158512, 0.123263, -0.789586])
chi2.append(1.1251)
fit_res_fixMuQ.append([0.649775, 0.100175])
chi2_fixmuQ.append(18.9201)
fit_res_fixMuQ_nopions.append([0.640038, 0.10016])
chi2_fixmuQ_nopions.append(1.38545)
fit_res.append([0.566466, 0.733513, 0.161216, 0.124064, -0.793823])
chi2.append(0.833115)
fit_res_fixMuQ.append([1.32886, 0.100055])
chi2_fixmuQ.append(38.1696)
fit_res_fixMuQ_nopions.append([1.33895, 0.100139])
chi2_fixmuQ_nopions.append(0.187315)
fit_res.append([0.849185, 0.00272185, 0.159554, 0.123501, -0.79107])
chi2.append(1.29216)
fit_res_fixMuQ.append([0.859691, 0.100126])
chi2_fixmuQ.append(1.29994)
fit_res_fixMuQ_nopions.append([0.860066, 0.100146])
chi2_fixmuQ_nopions.append(1.26343)

fit_res.append([0.716453, 0.307134, 0.150533, 0.111776, -0.763182])
chi2.append(4.39531)
fit_res_fixMuQ.append([1.0387, 0.099562])
chi2_fixmuQ.append(12.7941)
fit_res_fixMuQ_nopions.append([1.04439, 0.0995929])
chi2_fixmuQ_nopions.append(3.99529)
fit_res.append([0.255142, 1.20029, 0.151855, 0.112048, -0.76545])
chi2.append(4.78801)
fit_res_fixMuQ.append([1.49255, 0.099395])
chi2_fixmuQ.append(125.642)
fit_res_fixMuQ_nopions.append([1.51484, 0.0995692])
chi2_fixmuQ_nopions.append(2.69494)
fit_res.append([1.17545, -0.589442, 0.14894, 0.111434, -0.76047])
chi2.append(7.56329)
fit_res_fixMuQ.append([0.585999, 0.099653])
chi2_fixmuQ.append(34.873)
fit_res_fixMuQ_nopions.append([0.576592, 0.0987852])
chi2_fixmuQ_nopions.append(5)
fit_res.append([0.575151, 0.672392, 0.151215, 0.111964, -0.764459])
chi2.append(4.84168)
fit_res_fixMuQ.append([1.27266, 0.0995114])
chi2_fixmuQ.append(43.4741)
fit_res_fixMuQ_nopions.append([1.28465, 0.0995833])
chi2_fixmuQ_nopions.append(4.09971)
fit_res.append([0.857477, -0.0586904, 0.149788, 0.111565, -0.761784])
chi2.append(4.12635)
fit_res_fixMuQ.append([0.80487, 0.0995784])
chi2_fixmuQ.append(4.23564)
fit_res_fixMuQ_nopions.append([0.804083, 0.0996017])
chi2_fixmuQ_nopions.append(4.10866)

fit_res.append([0.857839, -0.325634, 0.138105, 0.119763, -0.893128])
chi2.append(3.47322)
fit_res_fixMuQ.append([0.524393, 0.0638413])
chi2_fixmuQ.append(10.5192)
fit_res_fixMuQ_nopions.append([0.521816, 0.0638337])
chi2_fixmuQ_nopions.append(3.1146)
fit_res.append([0.414803, 0.573214, 0.139981, 0.120769, -0.895342])
chi2.append(2.02227)
fit_res_fixMuQ.append([1.01937, 0.0638168])
chi2_fixmuQ.append(25.7632)
fit_res_fixMuQ_nopions.append([1.02206, 0.0638276])
chi2_fixmuQ_nopions.append(1.5424)
fit_res.append([1.29707, -1.22584, 0.135825, 0.118479, -0.89032])
chi2.append(9)
fit_res_fixMuQ.append([0.0285132, 0.0639634])
chi2_fixmuQ.append(115.891)
fit_res_fixMuQ_nopions.append([0.022766, 0.0636249])
chi2_fixmuQ_nopions.append(5)
fit_res.append([0.810884, -0.203695, 0.138361, 0.119904, -0.893444])
chi2.append(3.18565)
fit_res_fixMuQ.append([0.604379, 0.0638352])
chi2_fixmuQ.append(5.77991)
fit_res_fixMuQ_nopions.append([0.602757, 0.0638383])
chi2_fixmuQ_nopions.append(2.97751)
fit_res.append([0.90485, -0.44765, 0.137778, 0.119566, -0.892716])
chi2.append(3.66621)
fit_res_fixMuQ.append([0.444308, 0.0638519])
chi2_fixmuQ.append(17.433)
fit_res_fixMuQ_nopions.append([0.440695, 0.0638423])
chi2_fixmuQ_nopions.append(3.38683)

fit_res.append([0.829428, 0.141044, 0.0948749, 0.0733367, -0.801046])
chi2.append(2.96954)
fit_res_fixMuQ.append([0.985061, 0.0582218])
chi2_fixmuQ.append(7.46993)
fit_res_fixMuQ_nopions.append([0.987111, 0.0582322])
chi2_fixmuQ_nopions.append(2.95776)
fit_res.append([0.333241, 1.04078, 0.0958088, 0.0735831, -0.803324])
chi2.append(6.33835)
fit_res_fixMuQ.append([1.42343, 0.058154])
chi2_fixmuQ.append(217.516)
fit_res_fixMuQ_nopions.append([1.43765, 0.0581413])
chi2_fixmuQ_nopions.append(5)
fit_res.append([1.32183, -0.761355, 0.0937961, 0.0730418, -0.798378])
chi2.append(2.06972)
fit_res_fixMuQ.append([0.549768, 0.0582754])
chi2_fixmuQ.append(109.286)
fit_res_fixMuQ_nopions.append([0.536124, 0.0582323])
chi2_fixmuQ_nopions.append(0.99841)
fit_res.append([0.734428, 0.385186, 0.0952013, 0.0734494, -0.801931])
chi2.append(2.53909)
fit_res_fixMuQ.append([1.14445, 0.0582122])
chi2_fixmuQ.append(32.8953)
fit_res_fixMuQ_nopions.append([1.14956, 0.0582286])
chi2_fixmuQ_nopions.append(2.74959)
fit_res.append([0.924239, -0.103301, 0.0945198, 0.073207, -0.800065])
chi2.append(3.34867)
fit_res_fixMuQ.append([0.826315, 0.058229])
chi2_fixmuQ.append(4.88989)
fit_res_fixMuQ_nopions.append([0.82472, 0.0582401])
chi2_fixmuQ_nopions.append(3.27125)

fit_res.append([0.793914, -0.34432, 0.204878, 0.179192, -0.90296])
chi2.append(1.61171)
fit_res_fixMuQ.append([0.437895, 0.0904968])
chi2_fixmuQ.append(5.20226)
fit_res_fixMuQ_nopions.append([0.435359, 0.0904916])
chi2_fixmuQ_nopions.append(1.3741)
fit_res.append([0.307313, 0.560639, 0.207649, 0.180695, -0.90497])
chi2.append(0.239696)
fit_res_fixMuQ.append([0.898222, 0.0904673])
chi2_fixmuQ.append(10.3544)
fit_res_fixMuQ_nopions.append([0.900629, 0.0904704])
chi2_fixmuQ_nopions.append(0.423425)
fit_res.append([1.27693, -1.25102, 0.203126, 0.178464, -0.90174])
chi2.append(4.03199)
fit_res_fixMuQ.append([-0.0236317, 0.0906651])
chi2_fixmuQ.append(52.7105)
fit_res_fixMuQ_nopions.append([-0.030415, 0.090523])
chi2_fixmuQ_nopions.append(2.98998)
fit_res.append([0.746567, -0.222153, 0.20541, 0.179533, -0.903379])
chi2.append(1.55846)
fit_res_fixMuQ.append([0.518375, 0.090482])
chi2_fixmuQ.append(2.93424)
fit_res_fixMuQ_nopions.append([0.516817, 0.0904774])
chi2_fixmuQ_nopions.append(1.31862)
fit_res.append([0.841117, -0.466472, 0.204543, 0.179011, -0.902705])
chi2.append(1.74637)
fit_res_fixMuQ.append([0.357365, 0.0905119])
chi2_fixmuQ.append(8.41588)
fit_res_fixMuQ_nopions.append([0.353942, 0.0904871])
chi2_fixmuQ_nopions.append(1.35242)

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
ellipses_polarity = []
points = []

for i_cent, cent in enumerate(CENTRALITY_CLASSES):
    muB = fit_res[int(5*i_cent)][0]
    muB_err_corr = np.abs(fit_res[int(5*i_cent+2)][0]-fit_res[int(5*i_cent+1)][0])*0.5
    muI = fit_res[int(5*i_cent)][1]
    muI_err_corr = np.abs(fit_res[int(5*i_cent+2)][1]-fit_res[int(5*i_cent+1)][1])*0.5
    corr = fit_res[int(5*i_cent)][4]
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
    muB = fit_res[int(5*i_cent)][0]
    muB_err_polarity = np.abs(fit_res[int(5*i_cent+4)][0]-fit_res[int(5*i_cent+3)][0])*0.5
    muI = fit_res[int(5*i_cent)][1]
    muI_err_polarity = np.abs(fit_res[int(5*i_cent+4)][1]-fit_res[int(5*i_cent+3)][1])*0.5
    polarity = fit_res[int(5*i_cent)][4]
    ellipses_polarity.append(ROOT.RooEllipse(f"polarityEllipse_{cent[0]}_{cent[1]}",muB,muI,muB_err_polarity,muI_err_polarity,polarity))
    ellipses_polarity[i_cent].SetLineColor(centrality_colors[i_cent])
    ellipses_polarity[i_cent].SetLineStyle(ROOT.kDotted)
    ellipses_polarity[i_cent].SetLineWidth(2)
    ellipses_polarity[i_cent].Draw("lsame")

ellipses_polarity.append(ROOT.RooEllipse("dummy_polarity",-100,-100,1,1,0))
ellipses_polarity[5].SetLineStyle(ROOT.kDotted)
ellipses_polarity[5].SetLineWidth(2)
ellipses_polarity[5].SetLineColor(ROOT.kBlack)

for i_cent, cent in enumerate(CENTRALITY_CLASSES):
    muB = fit_res[int(5*i_cent)][0]
    muB_err = fit_res[int(5*i_cent)][2]
    muI = fit_res[int(5*i_cent)][1]
    muI_err = fit_res[int(5*i_cent)][3]
    corr = fit_res[int(5*i_cent)][4]
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
text.DrawLatex(1.20355,1.25,"ALICE")
text.SetTextSize(28)
text.DrawLatex(0.740964,0.947157,"Pb-Pb #sqrt{s_{NN}}=5.02 TeV")
#text.SetTextSize(28)
#text.DrawLatex(1.19,0.69795400,"|y| < 0.5")

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
legErr.AddEntry(ellipses_polarity[5],"B field uncert.", "f")
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