#!/usr/bin/env python3

import ROOT
import numpy as np
ROOT.gStyle.SetHistTopMargin(0)
path_he3 = './He3_PbPb/out'
path_triton = './Triton_PbPb/out'
path_hyp = './Hypertriton_PbPb'
path_proton = './Proton_PbPb/out'
path_pion = './Pion_PbPb/out'
path_omega = './Omega_PbPb'
centrality_classes = [[0, 5], [5, 10], [30, 50], [10, 30], [50,90]]
centrality_colors = [ROOT.kRed, ROOT.kOrange-3,ROOT.kAzure+4,ROOT.kGreen+2, ROOT.kMagenta+2]
particle_ratios = ["#bar{#Omega}^{+} / #Omega^{-}","#pi^{-} / #pi^{+}","#bar{p} / p","{}_{#bar{#Lambda}}^{3}#bar{H} / ^{3}_{#Lambda}H","^{3}#bar{H} / ^{3}H","^{3}#bar{He} / ^{3}He"]
n_part = [2.5,7.5,40,20,70]
n_part_err = [1.,1.,1.,1.,1.]

fixmuQ = False
nopions = False

chi2 = [1.0596,4.77823,3.71918,3.40207,1.77357]
chi2_fixmuQ = [10.7356,12.7981,10.4624,7.41971,5.14968]
chi2_fixmuQ_nopions = [0.695672,4.0161,3.03759,2.90952,1.26367]

ndf = [4,4,4,4,3]

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

fit_val=[[0.997025,0.986166,0.963691,0.968228,0.968138,0.995232],
[0.997547,0.986858,0.964858,0.968599,0.968861,0.995346],
[1.00266,0.993391,0.975557,0.97147,0.975231,0.996352],
[0.998845,0.987619,0.964911,0.966667,0.967986,0.995227],
[1.00281,0.994454,0.978471,0.974147,-1,0.996745]]

fit_val_fixmuQ=[[1.00022,0.986473,0.959204,0.958864,0.961662,0.994283],
[1.00021,0.987157,0.961242,0.960918,0.963579,0.994573],
[1.00011,0.993497,0.980246,0.980079,0.981448,0.997257],
[1.0002,0.987818,0.963211,0.962903,0.965431,0.994854],
[1.00009,0.994564,0.983469,0.983329,-1,0.997708]]

fit_val_fixmuQ_nopions=[[1.00022,0.98641,0.95902,0.958678,0.961488,0.994257],
[1.00021,0.987088,0.961037,0.960711,0.963386,0.994544],
[1.00011,0.99353,0.980346,0.98018,0.981543,0.997271],
[1.0002,0.987793,0.963136,0.962828,0.965361,0.994843],
[1.00009,0.994593,0.983556,0.983417,-1,0.99772]]

def muBmuQ(fixmuQ,nopions,i_cent):
    fit_r = fit_res
    if fixmuQ:
        fit_r = fit_res_fixMuQ
        if nopions:
            fit_r = fit_res_fixMuQ_nopions
    if not fixmuQ and not nopions:
        muB = fit_r[i_cent*3][0]
        muB_err = fit_r[i_cent*3][2]
        muB_corr = np.abs(fit_r[(i_cent)*3+1][0]-fit_r[(i_cent)*3+2][0])*0.5
        muQ = fit_r[i_cent*3][1]
        muQ_err = fit_r[i_cent*3][3]
        muQ_corr = np.abs(fit_r[(i_cent)*3+1][1]-fit_r[(i_cent)*3+2][1])*0.5
    else:
        muB = fit_r[i_cent*3][0]
        muB_err = fit_r[i_cent*3][1]
        muB_corr = np.abs(fit_r[(i_cent)*3+1][0]-fit_r[(i_cent)*3+2][0])*0.5
        muQ = -muB/50.
        muQ_err = muB_err/50.
        muQ_corr = 0
    return [muB,muB_err,muB_corr,muQ,muQ_err,muQ_corr]

apply_scaling_radius = False

TLATEX_TEXT_SIZE = 28
CENT_LIMIT_NUCLEI = 4

ROOT.gStyle.SetPadTickY(1)
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)
ROOT.gStyle.SetTextFont(44)

file_out = ROOT.TFile.Open('FinalPlot3D_new.root', 'recreate')
if apply_scaling_radius:
    file_out = ROOT.TFile.Open('FinalPlot3D_LHC22b9.root', 'recreate')
file_hijing = ROOT.TFile.Open('HIJINGRatios.root')

hMuBCentFrame = ROOT.TH1D("frame",";Centrality (%);#mu_{#it{B}} (MeV)",1,0,90)
hMuQCentFrame = ROOT.TH1D("frame",";#LT#it{N}_{part}#GT;#mu_{#it{Q}} (MeV)",1,0,450)
gMuBCent = ROOT.TGraphErrors()
gMuBCentCorrUnc = ROOT.TGraphErrors()
gMuQCent = ROOT.TGraphErrors()
gMuQCentCorrUnc = ROOT.TGraphErrors()

for i_cent, cent in enumerate(centrality_classes):
    fixmuQ_nopions = ''
    if fixmuQ:
        fixmuQ_nopions = 'fixmuQ'
    if fixmuQ and nopions:
        fixmuQ_nopions = 'fixmuQ_nopions'

    # get histograms
    file_material_budget = ROOT.TFile.Open("MaterialBudgetUncertainty.root")
    file_material_budget_hyp = ROOT.TFile.Open("MaterialBudgetUncertaintyHyp.root")
    file_he3 = ROOT.TFile.Open(path_he3 + '/SpectraHe3_kINT7.root')
    file_triton = ROOT.TFile.Open(path_triton + '/SpectraHe3.root')
    file_hyp = ROOT.TFile.Open(path_hyp + '/Ratio.root')
    file_proton = ROOT.TFile.Open(path_proton + '/out_tof_extend_4_old/SystematicsAllEPtNotCombinedTOF_extend_4.root')
    file_pion = ROOT.TFile.Open(path_pion + '/SystematicsAllEPtNotCombined_extend2.root')
    file_omega = ROOT.TFile.Open(path_omega + '/ratio_cutCompetingMass-3.root')
    file_he3_syst = ROOT.TFile.Open(path_he3 + '/SystematicsAll_extend2.root')
    file_triton_syst = ROOT.TFile.Open(path_triton + '/SystematicsAll_extend.root')
    file_he3_syst_eff_prim = ROOT.TFile.Open(path_he3 + '/SystematicsEfficiencyPrimary.root')
    file_hyp_syst = ROOT.TFile.Open(path_hyp + '/SystematicsRatio.root')
    # file_he3_scaling = ROOT.TFile.Open('He3_PbPb/MaterialRadiusCheck.root')
    # file_proton_scaling = ROOT.TFile.Open('Proton_PbPb/MaterialRadiusCheck.root')
    # file_pion_scaling = ROOT.TFile.Open('Pion_PbPb/MaterialRadiusCheck.root')
    file_out.cd()
    ratio_he3 = file_he3.Get(f'1.0_89_0.1_2.5_1_1_1/fRatio_{cent[0]}_{cent[1]}')
    ratio_proton = file_proton.Get(f'fRatio_{cent[0]}_{cent[1]}')
    ratio_pion = file_pion.Get(f'fRatio_{cent[0]}_{cent[1]}')
    ratio_omega = file_omega.Get(f'h_ratio_{cent[0]}_{cent[1]}')
    ratio_triton = ROOT.TH1D()
    ratio_hyp = ROOT.TH1D()
    ratio_triton_distribution = ROOT.TH1D()
    ratio_hyp_distribution = ROOT.TH1D()
    if i_cent < CENT_LIMIT_NUCLEI:
        ratio_hyp = file_hyp.Get(f'fRatio_{cent[0]}_{cent[1]}')
    ratio_triton = file_triton.Get(f'1.0_69_0.1_2.5_1_1_1/fRatio_{cent[0]}_{cent[1]}')
    ratio_triton_distribution = file_triton_syst.Get(f'hist/fFitPar_{cent[0]}_{cent[1]}')
    ratio_hyp_distribution = file_hyp_syst.Get(f'fParameterDistribution_{cent[0]}_{cent[1]}')
    ratio_he3_distribution = file_he3_syst.Get(f'hist/fFitPar_{cent[0]}_{cent[1]}')
    ratio_he3_distribution_eff_prim = file_he3_syst_eff_prim.Get(f'fRatioDistribution_{cent[0]}_{cent[1]}')
    ratio_proton_pt_correlated = file_proton.Get(f'fRatioDistributionTrials_{cent[0]}_{cent[1]}')
    ratio_pion_pt_correlated = file_pion.Get(f'fRatioDistributionTrials_{cent[0]}_{cent[1]}')
    ratio_omega_distribution = file_omega.Get(f'h_sys;{1+2*i_cent}')

    ratio_he3.Fit("pol0","QR","",2.,10.)
    ratio_he3.Write()

    # triton
    ratio_triton.Fit("pol0","QR","",1.6,3.)
    ratio_triton.Write()

    if i_cent < CENT_LIMIT_NUCLEI:
        ratio_hyp.Fit("pol0","QR","",0.,35.)
        ratio_hyp.Write()

    # proton
    ratio_proton.Fit("pol0","QR","",0.5,3.)
    ratio_proton.Write()

    # pion
    ratio_pion.Fit("pol0","QR","",0.7,1.6)
    ratio_pion.Write()

    # omega
    ratio_omega.Fit("pol0","QR","",1,10)
    ratio_omega.Write()

    # get fit functions
    fit_he3 = ratio_he3.GetFunction("pol0")
    fit_proton = ratio_proton.GetFunction("pol0")
    fit_pion = ratio_pion.GetFunction("pol0")
    fit_omega = ratio_omega.GetFunction("pol0")
    fit_triton = ratio_triton.GetFunction("pol0")
    fit_hyp = ROOT.TF1()
    if i_cent < CENT_LIMIT_NUCLEI:
        fit_hyp = ratio_hyp.GetFunction("pol0")

    # get ratios and errors
    ratio_he3 = fit_he3.GetParameter(0)
    ratio_he3_err = fit_he3.GetParError(0)/ratio_he3

    ratio_triton_err = 0.
    ratio_hyp_err = 0.
    ratio_triton = fit_triton.GetParameter(0)
    ratio_triton_err = fit_triton.GetParError(0)/ratio_triton
    if i_cent < CENT_LIMIT_NUCLEI:
        ratio_hyp = fit_hyp.GetParameter(0)
        ratio_hyp_err = fit_hyp.GetParError(0)/ratio_hyp
    
    ratio_proton = fit_proton.GetParameter(0)
    ratio_proton_err = 0#fit_proton.GetParError(0)

    ratio_pion = fit_pion.GetParameter(0)
    ratio_pion_err = 0#fit_pion.GetParError(0)

    ratio_omega = fit_omega.GetParameter(0)
    ratio_omega_err = fit_omega.GetParError(0)

    # systematic error
    syst_he3 = ratio_he3_distribution.GetRMS()/ratio_he3
    syst_he3_eff_prim = ratio_he3_distribution_eff_prim.GetRMS()/ratio_he3
    syst_omega = ratio_omega_distribution.GetRMS()
    syst_triton = 0.
    syst_hyp = 0.
    syst_triton = ratio_triton_distribution.GetRMS()/ratio_triton
    if i_cent < CENT_LIMIT_NUCLEI:
        syst_hyp = ratio_hyp_distribution.GetRMS()/ratio_hyp

    # uncorrelated sys
    syst_he3 = syst_he3*ratio_he3
    syst_he3_eff_prim = syst_he3_eff_prim*ratio_he3
    syst_he3 = np.sqrt(syst_he3*syst_he3+syst_he3_eff_prim*syst_he3_eff_prim)
    syst_triton = syst_triton*ratio_triton
    syst_triton = np.sqrt(syst_triton*syst_triton)
    if i_cent < CENT_LIMIT_NUCLEI:
        syst_hyp = syst_hyp*ratio_hyp
        syst_hyp = np.sqrt(syst_hyp*syst_hyp)
    syst_proton = fit_proton.GetParError(0)
    syst_pion = fit_pion.GetParError(0)
    
    # final plot
    stat_proton = ratio_proton_err
    stat_pion = ratio_pion_err
    stat_he3 = ratio_he3_err
    stat_hyp = 0
    stat_triton = ratio_triton_err
    if i_cent < CENT_LIMIT_NUCLEI:
        stat_hyp = ratio_hyp_err
    ratios_vs_b = ROOT.TH2D(f'fRatio_vs_b_{cent[0]}_{cent[1]}', ';#it{B}+#it{S}/3;#it{Q}-#it{S}/3;Antimatter / Matter', 10, -0.5, 9.5, 7, -0.5, 6.5)
    ratios_vs_b.SetBinContent(10, 7, ratio_he3)
    ratios_vs_b.SetBinError(10, 7, np.sqrt(syst_he3*syst_he3+stat_he3*stat_he3))
    ratios_vs_b.SetBinContent(4, 4, ratio_proton)
    ratios_vs_b.SetBinError(4, 4, np.sqrt(syst_proton*syst_proton+stat_proton*stat_proton))
    ratios_vs_b.SetBinContent(1, 4, ratio_pion)
    ratios_vs_b.SetBinError(1, 4, np.sqrt(syst_pion*syst_pion+stat_pion*stat_pion))
    ratios_vs_b.SetBinContent(1, 1, ratio_omega)
    ratios_vs_b.SetBinError(1, 1, np.sqrt(syst_omega*syst_omega+ratio_omega_err*ratio_omega_err))
    ratios_vs_b.SetBinContent(10, 4, ratio_triton)
    ratios_vs_b.SetBinError(10, 4, np.sqrt(syst_triton*syst_triton+stat_triton*stat_triton))
    if i_cent < CENT_LIMIT_NUCLEI:
        ratios_vs_b.SetBinContent(9, 5, ratio_hyp)
        ratios_vs_b.SetBinError(9, 5, np.sqrt(syst_hyp*syst_hyp+stat_hyp*stat_hyp))

    # set labels
    a = ratios_vs_b.GetXaxis()
    a.SetBit(ROOT.TAxis.kLabelsHori)
    a.SetBinLabel(1, '0')
    a.SetBinLabel(4, '1')
    a.SetBinLabel(7, '2')
    a.SetBinLabel(10, '3')
    a.SetLabelSize(0.062)

    aa = ratios_vs_b.GetYaxis()
    aa.SetBit(ROOT.TAxis.kLabelsHori)
    aa.SetBinLabel(1, '0')
    aa.SetBinLabel(4, '1')
    aa.SetBinLabel(7, '2')
    aa.SetLabelSize(0.062)

    # chi2 text
    formatted_chi2 = "{:.2f}".format(chi2[i_cent])
    if fixmuQ:
        formatted_chi2 = "{:.2f}".format(chi2_fixmuQ[i_cent])
        if nopions:
            formatted_chi2 = "{:.2f}".format(chi2_fixmuQ_nopions[i_cent])
    dof = ndf[i_cent]
    if fixmuQ and not nopions:
        dof = ndf[i_cent] + 1
    text_chi2 = ROOT.TLatex(-0.27, -0.44, "#chi^{2}/NDF = "+formatted_chi2+"/"+str(dof))
    text_chi2.SetTextSize(TLATEX_TEXT_SIZE)
    text_chi2.SetTextColor(ROOT.kBlack)

    # write to file
    ratios_vs_b.Write()

    # make final plot for approval (preliminary)
    leg_ratios_particle = ROOT.TLegend(0.164659,0.291729,0.359438,0.44812)
    leg_ratios_particle.SetTextFont(43)
    leg_ratios_particle.SetTextSize(22)
    cRatiosParticle = ROOT.TCanvas(f"cRatiosParticle_{cent[0]}_{cent[1]}",f"cRatiosParticle_{cent[0]}_{cent[1]}",500,500)
    pad1 = ROOT.TPad("pad1","pad1",0.0,0.3,1.0,1.0,0)
    pad1.SetFillColor(0)
    pad1.SetFrameBorderMode(0)
    pad1.SetBottomMargin(0.0)
    pad1.Draw()
    pad2 = ROOT.TPad("pad2", "pad2", 0.0, 0.0, 1.0, 0.3, 0)
    pad2.SetFillColor(0)
    pad2.SetFrameBorderMode(0)
    pad2.SetFillColor(0)
    pad2.SetFrameBorderMode(0)
    pad2.SetTopMargin(0.)
    pad2.SetBottomMargin(0.27)
    pad2.SetGridy()
    pad2.Draw()
    pad1.cd()
    hRatiosParticle = ROOT.TH1D(f"hRatiosParticle_{cent[0]}_{cent[1]}",f" ",6,0,6)
    hNSigmaRatioFitParticle = ROOT.TH1D(f"hNSigmaRatioFitParticle_{cent[0]}_{cent[1]}",f"{cent[0]}-{cent[1]}%",6,0,6)
    hRatiosParticleFit = ROOT.TH1D(f"hRatiosParticleFit_{cent[0]}_{cent[1]}",f"{cent[0]}-{cent[1]}%",6,0,6)
    hRatiosParticle.GetYaxis().SetTitle("Ratio")
    hRatiosParticle.GetYaxis().CenterTitle()
    hRatiosParticle.GetXaxis().SetTickSize(0.07*0.4)
    hRatiosParticle.GetYaxis().SetTitleSize(0.13*0.45)
    hRatiosParticle.GetYaxis().SetLabelSize(0.12*0.42)
    hRatiosParticle.GetYaxis().SetTitleOffset(0.9)
    for i_part in range(0,6):
        hRatiosParticle.GetXaxis().SetBinLabel(i_part+1,particle_ratios[i_part])
        hRatiosParticle.GetXaxis().SetLabelSize(0.06)
        ratio = -9999
        ratio_err = 0
        fit = -999
        if i_part == 0:
            ratio = ratio_omega
            ratio_err = np.sqrt(syst_omega*syst_omega+ratio_omega_err*ratio_omega_err)
            fit = fit_val[i_cent][5]
            if fixmuQ:
                fit = fit_val_fixmuQ[i_cent][5]
            if fixmuQ and nopions:
                fit = fit_val_fixmuQ_nopions[i_cent][5]
        if i_part == 1:
            ratio = ratio_pion
            ratio_err = np.sqrt(syst_pion*syst_pion+stat_pion*stat_pion)
            fit = fit_val[i_cent][0]
            if fixmuQ:
                fit = fit_val_fixmuQ[i_cent][0]
            if fixmuQ and nopions:
                fit = fit_val_fixmuQ_nopions[i_cent][0]
        if i_part == 2:
            ratio = ratio_proton
            ratio_err = np.sqrt(syst_proton*syst_proton+stat_proton*stat_proton)
            fit = fit_val[i_cent][1]
            if fixmuQ:
                fit = fit_val_fixmuQ[i_cent][1]
            if fixmuQ and nopions:
                fit = fit_val_fixmuQ_nopions[i_cent][1]
        if i_part == 3 and i_cent < CENT_LIMIT_NUCLEI:
            ratio = ratio_hyp
            ratio_err = np.sqrt(syst_hyp*syst_hyp+stat_hyp*stat_hyp)
            fit = fit_val[i_cent][4]
            if fixmuQ:
                fit = fit_val_fixmuQ[i_cent][4]
            if fixmuQ and nopions:
                fit = fit_val_fixmuQ_nopions[i_cent][4]
        if i_part == 4:
            ratio = ratio_triton
            ratio_err = np.sqrt(syst_triton*syst_triton+stat_triton*stat_triton)
            fit = fit_val[i_cent][3]
            if fixmuQ:
                fit = fit_val_fixmuQ[i_cent][3]
            if fixmuQ and nopions:
                fit = fit_val_fixmuQ_nopions[i_cent][3]
        if i_part == 5:
            ratio = ratio_he3
            ratio_err = np.sqrt(syst_he3*syst_he3+stat_he3*stat_he3)
            fit = fit_val[i_cent][2]
            if fixmuQ:
                fit = fit_val_fixmuQ[i_cent][2]
            if fixmuQ and nopions:
                fit = fit_val_fixmuQ_nopions[i_cent][2]
        hRatiosParticle.SetBinContent(i_part+1,ratio)
        hRatiosParticle.SetBinError(i_part+1,ratio_err)
        hRatiosParticleFit.SetBinContent(i_part+1,fit)
        hRatiosParticleFit.SetBinError(i_part+1,0)
        print(f"fit = {fit}")
    hRatiosParticle.GetYaxis().SetNdivisions(10)
    hRatiosParticle.GetYaxis().SetRangeUser(0.67,1.18)
    gRatiosParticle = ROOT.TGraphErrors(hRatiosParticle)
    gRatiosParticle.SetName(f"gRatiosParticle_{cent[0]}_{cent[1]}")
    gRatiosParticleFit = ROOT.TGraphErrors(hRatiosParticleFit)
    gRatiosParticleFit.SetName(f"gRatiosParticleFit_{cent[0]}_{cent[1]}")
    for i_part in range(0,6):
        gRatiosParticle.SetPointError(i_part,0,hRatiosParticle.GetBinError(i_part+1))
        gRatiosParticleFit.SetPointError(i_part,0.3,0)
    
    # data to fit ratio
    for i_ticks in range(6):
        hNSigmaRatioFitParticle.GetYaxis().SetNdivisions(5)
    for i_part in range(0,6):
        hNSigmaRatioFitParticle.GetXaxis().SetBinLabel(i_part+1,particle_ratios[i_part])
        hNSigmaRatioFitParticle.GetXaxis().SetLabelSize(0.2)
        ratio = -999.
        ratio_err = -1.
        fit = -9999.
        if i_part == 0:
            ratio = ratio_omega
            ratio_err = np.sqrt(syst_omega*syst_omega+ratio_omega_err*ratio_omega_err)
            fit = fit_val[i_cent][5]
            if fixmuQ:
                fit = fit_val_fixmuQ[i_cent][5]
            if fixmuQ and nopions:
                fit = fit_val_fixmuQ_nopions[i_cent][5]
        if i_part == 1:
            ratio = ratio_pion
            ratio_err = np.sqrt(syst_pion*syst_pion+stat_pion*stat_pion)
            fit = fit_val[i_cent][0]
            if fixmuQ:
                fit = fit_val_fixmuQ[i_cent][0]
            if fixmuQ and nopions:
                fit = fit_val_fixmuQ_nopions[i_cent][0]
        if i_part == 2:
            ratio = ratio_proton
            ratio_err = np.sqrt(syst_proton*syst_proton+stat_proton*stat_proton)
            fit = fit_val[i_cent][1]
            if fixmuQ:
                fit = fit_val_fixmuQ[i_cent][1]
            if fixmuQ and nopions:
                fit = fit_val_fixmuQ_nopions[i_cent][1]
        if i_part == 3 and i_cent < CENT_LIMIT_NUCLEI:
            ratio = ratio_hyp
            ratio_err = np.sqrt(syst_hyp*syst_hyp+stat_hyp*stat_hyp)
            fit = fit_val[i_cent][4]
            if fixmuQ:
                fit = fit_val_fixmuQ[i_cent][4]
            if fixmuQ and nopions:
                fit = fit_val_fixmuQ_nopions[i_cent][4]
        if i_part == 4:
            ratio = ratio_triton
            ratio_err = np.sqrt(syst_triton*syst_triton+stat_triton*stat_triton)
            fit = fit_val[i_cent][3]
            if fixmuQ:
                fit = fit_val_fixmuQ[i_cent][3]
            if fixmuQ and nopions:
                fit = fit_val_fixmuQ_nopions[i_cent][3]
        if i_part == 5:
            ratio = ratio_he3
            ratio_err = np.sqrt(syst_he3*syst_he3+stat_he3*stat_he3)
            fit = fit_val[i_cent][2]
            if fixmuQ:
                fit = fit_val_fixmuQ[i_cent][2]
            if fixmuQ and nopions:
                fit = fit_val_fixmuQ_nopions[i_cent][2]
        if ratio > 0:
            hNSigmaRatioFitParticle.SetBinContent(i_part+1,(ratio-fit)/ratio_err)
            hNSigmaRatioFitParticle.SetBinError(i_part+1,0)
        else:
            hNSigmaRatioFitParticle.SetBinContent(i_part+1,-99999.)
            hNSigmaRatioFitParticle.SetBinError(i_part+1,0)
    hNSigmaRatioFitParticle.SetLineColor(centrality_colors[i_cent])
    hNSigmaRatioFitParticle.SetMarkerColor(centrality_colors[i_cent])
    hNSigmaRatioFitParticle.SetTitle("")
    hNSigmaRatioFitParticle.GetXaxis().SetTickSize(0.07)
    hNSigmaRatioFitParticle.GetYaxis().SetTickSize(0.035)
    hNSigmaRatioFitParticle.GetYaxis().SetTitleSize(0.13)
    hNSigmaRatioFitParticle.GetYaxis().SetLabelSize(0.12)
    hNSigmaRatioFitParticle.GetYaxis().SetTitleOffset(0.4)
    hNSigmaRatioFitParticle.GetXaxis().SetLabelOffset(0.01)
    # hNSigmaRatioFitParticle.SetLineWidth(1)
    hNSigmaRatioFitParticle.GetYaxis().SetTitle("pull")#"#frac{data - fit}{#sigma_{data}}")
    hNSigmaRatioFitParticle.GetYaxis().CenterTitle()
    hNSigmaRatioFitParticle.GetYaxis().SetRangeUser(-4.7,4.7)
    # gHijingNsigma = ROOT.TGraphErrors(gHijingPredictions)
    # gHijingNsigma.SetPointY(0,(gHijingNsigma.GetPointY(0)-fit_expo.Eval(0,1))/np.sqrt(syst_pion*syst_pion+stat_pion*stat_pion))
    # gHijingNsigma.SetPointY(1,(gHijingNsigma.GetPointY(1)-fit_expo.Eval(3,0.5))/np.sqrt(syst_proton*syst_proton+stat_proton*stat_proton))
    # hLineZero = ROOT.TLine(0.,0.,6.,0.)
    # hLineZero.SetLineColor(ROOT.kBlack)
    # hLineZero.SetLineStyle(ROOT.kDashed)
    # hLineOnePos = ROOT.TLine(0.,1.,6.,1.)
    # hLineOnePos.SetLineColor(ROOT.kBlack)
    # hLineOnePos.SetLineStyle(ROOT.kDashed)
    # hLineOneNeg = ROOT.TLine(0.,-1.,6.,-1.)
    # hLineOneNeg.SetLineColor(ROOT.kBlack)
    # hLineOneNeg.SetLineStyle(ROOT.kDashed)
    # hLineTwoPos = ROOT.TLine(0.,2.,6.,2.)
    # hLineTwoPos.SetLineColor(ROOT.kBlack)
    # hLineTwoPos.SetLineStyle(ROOT.kDashed)
    # hLineTwoNeg = ROOT.TLine(0.,-2.,6.,-2.)
    # hLineTwoNeg.SetLineColor(ROOT.kBlack)
    # hLineTwoNeg.SetLineStyle(ROOT.kDashed)
    hLineZero = ROOT.TLine(0.,0.,6.,0.)
    hLineZero.SetLineColor(ROOT.kBlack)
    hLineZero.SetLineStyle(ROOT.kDashed)
    hLineOnePos = ROOT.TLine(0.,2.,6.,2.)
    hLineOnePos.SetLineColor(ROOT.kBlack)
    hLineOnePos.SetLineStyle(ROOT.kDashed)
    hLineOneNeg = ROOT.TLine(0.,-2.,6.,-2.)
    hLineOneNeg.SetLineColor(ROOT.kBlack)
    hLineOneNeg.SetLineStyle(ROOT.kDashed)
    hLineTwoPos = ROOT.TLine(0.,4.,6.,4.)
    hLineTwoPos.SetLineColor(ROOT.kBlack)
    hLineTwoPos.SetLineStyle(ROOT.kDashed)
    hLineTwoNeg = ROOT.TLine(0.,-4.,6.,-4.)
    hLineTwoNeg.SetLineColor(ROOT.kBlack)
    hLineTwoNeg.SetLineStyle(ROOT.kDashed)

    hRatiosParticle.GetYaxis().SetRangeUser(0.45,1.3)
    text_alice = ROOT.TLatex(0.174185,0.17764,"ALICE")
    text_alice.SetNDC(True)
    text_alice.SetTextFont(43)
    text_alice.SetTextSize(28)
    text_energy = ROOT.TLatex(0.175439,0.0732919,"Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV")
    text_energy.SetNDC(True)
    text_energy.SetTextFont(43)
    text_energy.SetTextSize(25)
    text_centrality = ROOT.TLatex(0.18,0.8,f"{cent[0]}-{cent[1]}%")
    text_centrality.SetNDC(True)
    text_centrality.SetTextFont(43)
    text_centrality.SetTextSize(25)
    gRatiosParticle = ROOT.TGraphErrors(hRatiosParticle)
    gRatiosParticleFit = ROOT.TGraphErrors(hRatiosParticleFit)
    for i_part in range(0,6):
        gRatiosParticle.SetPointError(i_part,0,hRatiosParticle.GetBinError(i_part+1))
        gRatiosParticleFit.SetPointError(i_part,0.3,0)
    hRatiosParticle.SetLineColor(ROOT.kWhite)
    hRatiosParticle.SetMarkerColor(ROOT.kWhite)
    hRatiosParticle.Draw("")
    hRatiosParticle.GetYaxis().SetDecimals()
    hRatiosParticle.GetXaxis().SetDecimals()
    gRatiosParticle.SetMarkerStyle(20)
    gRatiosParticle.SetMarkerSize(1.4)
    gRatiosParticle.SetLineWidth(1)
    gRatiosParticle.SetLineColor(centrality_colors[i_cent])
    gRatiosParticle.SetMarkerColor(centrality_colors[i_cent])
    gRatiosParticleFit.SetMarkerSize(0)
    gRatiosParticleFit.SetLineWidth(2)
    gRatiosParticleFit.SetLineColor(ROOT.kBlack)
    # gHijingPredictions.SetLineColor(ROOT.kGreen+1)
    # gHijingPredictions.SetLineWidth(2)
    gRatiosParticleFit.Draw("e same")
    # gHijingPredictions.Draw("pe same")
    gRatiosParticle.Draw("pe same")
    leg_ratios_particle.AddEntry(gRatiosParticle,"data","pe")
    leg_ratios_particle.AddEntry(gRatiosParticleFit,"fit")

    muBmuQ_list = muBmuQ(fixmuQ,nopions,i_cent)
    text_mu_b = ROOT.TLatex(1,0,"#mu_{B}="+"{:.2f}".format(muBmuQ_list[0])+"#pm"+"{:.2f}".format(muBmuQ_list[1])+"(uncorr.)#pm"+"{:.2f}".format(muBmuQ_list[2])+"(corr.) MeV")
    text_mu_b.SetName("mu_b")
    text_mu_q = ROOT.TLatex(1.2,0,"#mu_{Q}="+"{:.2f}".format(muBmuQ_list[3])+"#pm"+"{:.2f}".format(muBmuQ_list[4])+"(uncorr.)#pm"+"{:.2f}".format(muBmuQ_list[5])+"(corr.) MeV")
    if fixmuQ:
        text_mu_q = ROOT.TLatex(1.2,0,"#mu_{Q}="+"{:.2f}".format(muBmuQ_list[3])+" MeV")
    text_mu_q.SetName("mu_q")
    text_chi2.SetName("chi2")
    text_centrality.SetName("centrality")
    # leg_ratios_particle.AddEntry(gHijingPredictions,"HIJING")
    text_mu_b.SetTextFont(43)
    text_mu_q.SetTextFont(43)
    text_chi2.SetTextFont(43)
    text_mu_b.SetTextSize(20)
    text_mu_q.SetTextSize(20)
    text_chi2.SetTextSize(20)
    text_mu_b.SetNDC(True)
    text_mu_q.SetNDC(True)
    text_chi2.SetNDC(True)
    text_mu_b.SetX(0.18)
    text_mu_b.SetY(0.16)
    text_mu_q.SetX(0.18)
    text_mu_q.SetY(0.07)
    text_chi2.SetX(0.18)
    text_chi2.SetY(0.24)
    text_mu_b.Draw("same")
    text_mu_q.Draw("same")
    text_chi2.Draw("same")
    #text_alice.Draw("same")
    #text_energy.Draw("same")
    text_centrality.Draw("same")
    leg_ratios_particle.Draw("same")
    pad2.cd()
    hNSigmaRatioFitParticle.SetMarkerStyle(20)
    hNSigmaRatioFitParticle.SetMarkerSize(1.4)
    hNSigmaRatioFitParticle.Draw("axis")
    hLineZero.Draw("same")
    hLineOnePos.Draw("same")
    hLineOneNeg.Draw("same")
    hLineTwoPos.Draw("same")
    hLineTwoNeg.Draw("same")
    # gHijingNsigma.SetLineColor(ROOT.kGreen+1)
    # gHijingNsigma.SetLineWidth(2)
    hNSigmaRatioFitParticle.Draw("psame")
    # gHijingNsigma.Draw("pe same")
    cRatiosParticle.Write()
    cRatiosParticle.Print(f"SHMfitToRatios_{cent[0]}_{cent[1]}{fixmuQ_nopions}.pdf")
    gRatiosParticle.Write()
    gRatiosParticleFit.Write()
    hNSigmaRatioFitParticle.Write()

    gMuBCent.AddPoint(n_part[i_cent],muBmuQ_list[0])
    gMuBCent.SetPointError(gMuBCent.GetN()-1,n_part_err[i_cent],muBmuQ_list[1])
    gMuBCentCorrUnc.AddPoint(n_part[i_cent],muBmuQ_list[0])
    gMuBCentCorrUnc.SetPointError(gMuBCent.GetN()-1,.5,muBmuQ_list[2])
    gMuQCent.AddPoint(n_part[i_cent],muBmuQ_list[3])
    gMuQCent.SetPointError(gMuBCent.GetN()-1,0.,muBmuQ_list[4])
    gMuQCentCorrUnc.AddPoint(n_part[i_cent],muBmuQ_list[3])
    gMuQCentCorrUnc.SetPointError(gMuQCent.GetN()-1,10,muBmuQ_list[5])
    # print(f'mu_I_error = {mu_I_error}')

gMuBCent.SetName("MuBCent")
gMuBCentCorrUnc.SetName("gMuBCentCorrUnc")

legMuBCent = ROOT.TLegend(0.18,0.15,0.60,0.28)
gPublishedResult = ROOT.TGraphErrors()
gPublishedResult.SetName("NaturePoint")
gPublishedResult.AddPoint(5,0.7)
gPublishedResult.SetPointError(0,5,3.8)
gPublishedResult.SetMarkerStyle(0)
gPublishedResult.SetMarkerSize(0)
gPublishedResult.SetLineWidth(0)
gPublishedResult.SetFillStyle(3345)
gPublishedResult.SetFillColor(ROOT.kBlack)
text_alice = ROOT.TLatex(0.18,0.8,"ALICE")
text_alice.SetTextFont(43)
text_alice.SetTextSize(28)
text_alice.SetNDC()
text_energy = ROOT.TLatex(0.18,0.73,"Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV")
text_energy.SetTextFont(43)
text_energy.SetTextSize(25)
text_energy.SetNDC()
text_uncorr_uncertainty = ROOT.TLatex(0.18,0.66,"Uncorrelated uncertainties only")
text_uncorr_uncertainty.SetTextFont(43)
text_uncorr_uncertainty.SetNDC()
cMuBCent = ROOT.TCanvas("cMuBCent","cMuBCent")
cMuBCent.cd()
gMuBCent.SetMarkerStyle(20)
gMuBCent.SetMarkerSize(1.2)
gMuBCent.SetMarkerColor(ROOT.kRed)
gMuBCent.SetLineColor(ROOT.kRed)
gMuBCent.SetFillColor(0)
gMuBCent.SetFillStyle(0)
gMuBCentCorrUnc.SetMarkerColor(ROOT.kRed-9)
gMuBCentCorrUnc.SetLineWidth(0)
gMuBCentCorrUnc.SetMarkerSize(0)
gMuBCentCorrUnc.SetLineColor(0)
gMuBCentCorrUnc.SetFillColor(ROOT.kRed-9)
hMuBCentFrame.GetYaxis().SetRangeUser(-3.6,5.0)
hMuBCentFrame.GetYaxis().SetNdivisions(510)
hMuBCentFrame.GetYaxis().SetDecimals()
hMuBCentFrame.GetXaxis().SetDecimals()
hMuBCentFrame.Draw("axis")
gPublishedResult.Draw("pe5same")
gMuBCentCorrUnc.Draw("pe5same")
gMuBCent.Draw("pe5same")
legMuBCent.SetNColumns(2)
legMuBCent.SetTextFont(43)
legMuBCent.SetTextSize(22)
legMuBCent.AddEntry(gMuBCent,"Uncorr. uncert.", "pf")
legMuBCent.AddEntry(gMuBCentCorrUnc,"Corr. uncert.")
legMuBCent.AddEntry(gPublishedResult,"SHM fit, #it{Nature} #bf{561}, 321-330 (2018)")
legMuBCent.Draw("same")
text_alice.Draw("same")
text_energy.Draw("same")
cMuBCent.Write()
cMuBCent.Print("MuBvsNpart.eps")

cB = ROOT.TCanvas("cMuBuncorr","cMuBuncorr")
cB.cd()
gMuBCent.Fit("pol0")
gMuBCent.GetYaxis().SetTitle("#mu_{B} (MeV)")
gMuBCent.GetXaxis().SetTitle("Centrality (%)")
gMuBCent.Draw("ape")


cQ = ROOT.TCanvas("cMuQuncorr","cMuQuncorr")
cQ.cd()
gMuQCent.SetName("MuQCent")
gMuQCent.SetMarkerStyle(20)
gMuQCent.SetMarkerSize(1.2)
gMuQCent.SetMarkerColor(ROOT.kRed)
gMuQCent.SetLineColor(ROOT.kRed)
gMuQCent.SetFillColor(0)
gMuQCent.SetFillStyle(0)
gMuQCent.Fit("pol0")
gMuQCent.GetYaxis().SetTitle("#mu_{Q} (MeV)")
gMuQCent.GetXaxis().SetTitle("Centrality (%)")
gMuQCent.Draw("ape")

cB.Write()
cQ.Write()
gMuBCent.Write()


file_out.Close()