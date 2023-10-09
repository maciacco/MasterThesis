import ROOT
import numpy as np

marker = [34, 45, 47, 33, 29]

chi2 = []
chi2_fixmuQ = []
chi2_fixmuQ_nopions = []
fit_res = []
fit_res_fixMuQ = []
fit_res_fixMuQ_nopions = []

fit_res.append([0.732082, 0.375965, 0.157341, 0.12004, -0.783795])
chi2.append(0.853811)
fit_res_fixMuQ.append([1.12657, 0.0999665, -0.0283236])
chi2_fixmuQ.append(11.6751)
fit_res_fixMuQ_nopions.append([1.13196, 0.100004, -0.0284593])
chi2_fixmuQ_nopions.append(0.515368)
fit_res.append([0.263754, 1.26683, 0.158628, 0.12032, -0.785797])
chi2.append(3.71891)
fit_res_fixMuQ.append([1.57095, 0.0997854, -0.0394962])
chi2_fixmuQ.append(120.671)
fit_res_fixMuQ_nopions.append([1.59169, 0.0999803, -0.0400177])
chi2_fixmuQ_nopions.append(1.89386)
fit_res.append([1.19837, -0.518325, 0.155481, 0.119549, -0.780839])
chi2.append(1.08175)
fit_res_fixMuQ.append([0.682777, 0.100033, -0.017166])
chi2_fixmuQ.append(19.358)
fit_res_fixMuQ_nopions.append([0.672583, 0.100017, -0.0169098])
chi2_fixmuQ_nopions.append(1.26066)
fit_res.append([0.590557, 0.741128, 0.158119, 0.120289, -0.785134])
chi2.append(0.716383)
fit_res_fixMuQ.append([1.36059, 0.0999093, -0.0342073])
chi2_fixmuQ.append(41.3203)
fit_res_fixMuQ_nopions.append([1.37137, 0.100009, -0.0344785])
chi2_fixmuQ_nopions.append(0.165534)
fit_res.append([0.873264, 0.0102989, 0.156495, 0.119764, -0.782333])
chi2.append(1.18213)
fit_res_fixMuQ.append([0.892101, 0.0999822, -0.0224288])
chi2_fixmuQ.append(1.19981)
fit_res_fixMuQ_nopions.append([0.892611, 0.100003, -0.0224416])
chi2_fixmuQ_nopions.append(1.13344)

fit_res.append([0.730364, 0.302259, 0.150009, 0.11095, -0.760148])
chi2.append(4.34376)
fit_res_fixMuQ.append([1.04772, 0.0997644, -0.0263413])
chi2_fixmuQ.append(12.6284)
fit_res_fixMuQ_nopions.append([1.05344, 0.0997956, -0.026485])
chi2_fixmuQ_nopions.append(3.96268)
fit_res.append([0.268336, 1.19585, 0.15133, 0.111216, -0.762433])
chi2.append(4.7174)
fit_res_fixMuQ.append([1.50084, 0.0995951, -0.0377335])
chi2_fixmuQ.append(126.55)
fit_res_fixMuQ_nopions.append([1.5235, 0.0997715, -0.0383032])
chi2_fixmuQ_nopions.append(2.67027)
fit_res.append([1.19, -0.594735, 0.148423, 0.110616, -0.757425])
chi2.append(7.53303)
fit_res_fixMuQ.append([0.595799, 0.0998555, -0.0149793])
chi2_fixmuQ.append(35.7424)
fit_res_fixMuQ_nopions.append([0.585925, 0.0989956, -0.014731])
chi2_fixmuQ_nopions.append(5)
fit_res.append([0.589086, 0.66753, 0.150689, 0.111135, -0.761433])
chi2.append(4.78684)
fit_res_fixMuQ.append([1.28152, 0.0997129, -0.0322195])
chi2_fixmuQ.append(43.4794)
fit_res_fixMuQ_nopions.append([1.29367, 0.099786, -0.032525])
chi2_fixmuQ_nopions.append(4.07446)
fit_res.append([0.871353, -0.0635727, 0.149266, 0.110743, -0.758742])
chi2.append(4.07776)
fit_res_fixMuQ.append([0.81406, 0.099781, -0.0204667])
chi2_fixmuQ.append(4.23337)
fit_res_fixMuQ_nopions.append([0.813155, 0.0998048, -0.0204439])
chi2_fixmuQ_nopions.append(4.07014)

fit_res.append([0.759228, -0.238615, 0.138515, 0.115424, -0.857544])
chi2.append(3.40606)
fit_res_fixMuQ.append([0.515186, 0.0732047, -0.0129526])
chi2_fixmuQ.append(7.37409)
fit_res_fixMuQ_nopions.append([0.512595, 0.0732049, -0.0128874])
chi2_fixmuQ_nopions.append(3.13772)
fit_res.append([0.317167, 0.659052, 0.140228, 0.116181, -0.860022])
chi2.append(2.15119)
fit_res_fixMuQ.append([1.00787, 0.0731753, -0.0253395])
chi2_fixmuQ.append(35.8421)
fit_res_fixMuQ_nopions.append([1.01276, 0.0731883, -0.0254624])
chi2_fixmuQ_nopions.append(1.54537)
fit_res.append([1.19678, -1.13738, 0.136317, 0.114356, -0.854184])
chi2.append(9)
fit_res_fixMuQ.append([0.0219876, 0.0733615, -0.000552801])
chi2_fixmuQ.append(108.149)
fit_res_fixMuQ_nopions.append([0.0142132, 0.0728957, -0.00035734])
chi2_fixmuQ_nopions.append(5)
fit_res.append([0.759228, -0.238615, 0.138515, 0.115424, -0.857544])
chi2.append(3.40606)
fit_res_fixMuQ.append([0.515186, 0.0732047, -0.0129526])
chi2_fixmuQ.append(7.37409)
fit_res_fixMuQ_nopions.append([0.512595, 0.0732049, -0.0128874])
chi2_fixmuQ_nopions.append(3.13772)
fit_res.append([0.759228, -0.238615, 0.138515, 0.115424, -0.857544])
chi2.append(3.40606)
fit_res_fixMuQ.append([0.515186, 0.0732047, -0.0129526])
chi2_fixmuQ.append(7.37409)
fit_res_fixMuQ_nopions.append([0.512595, 0.0732049, -0.0128874])
chi2_fixmuQ_nopions.append(3.13772)

fit_res.append([0.825807, 0.136206, 0.0942823, 0.0723968, -0.795644])
chi2.append(2.93335)
fit_res_fixMuQ.append([0.976131, 0.0585533, -0.0245414])
chi2_fixmuQ.append(7.30645)
fit_res_fixMuQ_nopions.append([0.978234, 0.0585657, -0.0245943])
chi2_fixmuQ_nopions.append(2.98952)
fit_res.append([0.329507, 1.0359, 0.0952184, 0.0726401, -0.797982])
chi2.append(6.34478)
fit_res_fixMuQ.append([1.41389, 0.058485, -0.0355474])
chi2_fixmuQ.append(220.935)
fit_res_fixMuQ_nopions.append([1.42864, 0.0584718, -0.0359184])
chi2_fixmuQ_nopions.append(5)
fit_res.append([1.31811, -0.766036, 0.0931977, 0.0721053, -0.792898])
chi2.append(1.946)
fit_res_fixMuQ.append([0.541527, 0.0586082, -0.0136148])
chi2_fixmuQ.append(113.489)
fit_res_fixMuQ_nopions.append([0.527268, 0.0585711, -0.0132563])
chi2_fixmuQ_nopions.append(1.09386)
fit_res.append([0.730795, 0.380368, 0.0946166, 0.0725128, -0.796577])
chi2.append(2.56288)
fit_res_fixMuQ.append([1.13538, 0.0585432, -0.0285452])
chi2_fixmuQ.append(32.922)
fit_res_fixMuQ_nopions.append([1.1406, 0.0585632, -0.0286765])
chi2_fixmuQ_nopions.append(2.79656)
fit_res.append([0.920577, -0.108123, 0.0939439, 0.07228, -0.794696])
chi2.append(3.37534)
fit_res_fixMuQ.append([0.817574, 0.0585581, -0.0205551])
chi2_fixmuQ.append(5.09777)
fit_res_fixMuQ_nopions.append([0.815846, 0.0585684, -0.0205116])
chi2_fixmuQ_nopions.append(3.23815)

fit_res.append([0.714549, -0.271576, 0.188885, 0.160705, -0.878581])
chi2.append(1.62153)
fit_res_fixMuQ.append([0.433449, 0.0926987, -0.0108976])
chi2_fixmuQ.append(4.33548)
fit_res_fixMuQ_nopions.append([0.430914, 0.0926928, -0.0108338])
chi2_fixmuQ_nopions.append(1.34692)
fit_res.append([0.228141, 0.632279, 0.191322, 0.161876, -0.880854])
chi2.append(0.281467)
fit_res_fixMuQ.append([0.892337, 0.0926727, -0.0224347])
chi2_fixmuQ.append(16.1368)
fit_res_fixMuQ_nopions.append([0.895879, 0.0926761, -0.0225238])
chi2_fixmuQ_nopions.append(0.416871)
fit_res.append([1.19689, -1.17718, 0.187196, 0.160096, -0.87708])
chi2.append(4.01259)
fit_res_fixMuQ.append([-0.0262955, 0.0928863, 0.000661107])
chi2_fixmuQ.append(57.7456)
fit_res_fixMuQ_nopions.append([-0.034709, 0.0927288, 0.000872637])
chi2_fixmuQ_nopions.append(2.9769)
fit_res.append([0.714549, -0.271576, 0.188885, 0.160705, -0.878581])
chi2.append(1.62153)
fit_res_fixMuQ.append([0.433449, 0.0926987, -0.0108976])
chi2_fixmuQ.append(4.33548)
fit_res_fixMuQ_nopions.append([0.430914, 0.0926928, -0.0108338])
chi2_fixmuQ_nopions.append(1.34692)
fit_res.append([0.714549, -0.271576, 0.188885, 0.160705, -0.878581])
chi2.append(1.62153)
fit_res_fixMuQ.append([0.433449, 0.0926987, -0.0108976])
chi2_fixmuQ.append(4.33548)
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
h_frame = ROOT.TH2F("h_frame"," ",1,-1.65,1.65,1,-1.65,1.65)
h_frame.GetXaxis().SetTitle("#it{#mu}_{B} (MeV)")
h_frame.GetYaxis().SetTitle("#it{#mu}_{Q} (MeV)")
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
    points[i_cent].SetMarkerStyle(marker[i_cent])
    points[i_cent].SetMarkerSize(2.0)
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
text.SetNDC()
text.DrawLatex(0.2,0.9,"ALICE")
text.SetTextSize(24)
text.DrawLatex(0.2,0.85,"Pb#minusPb, #sqrt{#it{s}_{NN}} = 5.02 TeV")
text.DrawLatex(0.2,0.45,"Thermal-FIST")
text.DrawLatex(0.2,0.4,"#it{T}_{ch} = 155 #pm 2 MeV")
text.DrawLatex(0.2,0.35,"#it{#mu}_{S} constrained")
#text.SetTextSize(28)
#text.DrawLatex(1.19,0.69795400,"|y| < 0.5")

# legend
#leg = ROOT.TLegend(0.152174,0.292969,0.307692,0.515625)
leg = ROOT.TLegend(0.2,0.55,0.4,0.8)
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

# axes
l = ROOT.TLine()

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