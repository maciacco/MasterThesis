#!/usr/bin/env python3

import ROOT
import numpy as np
ROOT.gStyle.SetHistTopMargin(0)
path_he3 = './He3_PbPb/out'
path_hyp = './Hypertriton_PbPb'
path_proton = './Proton_PbPb/out'
path_pion = './Pion_PbPb/out'
centrality_classes = [[0, 5], [5, 10], [30, 50]]
centrality_colors = [ROOT.kOrange+7, ROOT.kAzure+4, ROOT.kTeal+4]

TLATEX_TEXT_SIZE = 28

ROOT.gStyle.SetPadTickY(1)
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)
ROOT.gStyle.SetTextFont(44)

file_he3 = ROOT.TFile.Open(path_he3 + '/SpectraHe3.root')
file_hyp = ROOT.TFile.Open(path_hyp + '/Ratio.root')
file_proton = ROOT.TFile.Open(path_proton + '/SystematicsAllEPtNotCombined.root')
file_pion = ROOT.TFile.Open(path_pion + '/SystematicsAllEPtNotCombined.root')
file_he3_syst = ROOT.TFile.Open(path_he3 + '/SystematicsAll.root')
file_hyp_syst = ROOT.TFile.Open(path_hyp + '/Systematics.root')

file_out = ROOT.TFile.Open('FinalPlot3D.root', 'recreate')

for i_cent, cent in enumerate(centrality_classes):

    # get histograms
    ratio_he3 = file_he3.Get(f'1.0_89_1_1_1/fRatio_{cent[0]}_{cent[1]}')
    ratio_hyp = file_hyp.Get(f'fRatio_{cent[0]}_{cent[1]}')
    ratio_proton = file_proton.Get(f'fRatio_{cent[0]}_{cent[1]}')
    ratio_pion = file_pion.Get(f'fRatio_{cent[0]}_{cent[1]}')
    ratio_he3_distribution = file_he3_syst.Get(f'hist/fFitPar_{cent[0]}_{cent[1]}')
    ratio_hyp_distribution = file_hyp_syst.Get(f'fParameterDistribution_{cent[0]}_{cent[1]}')
    ratio_proton_pt_correlated = file_proton.Get(f'fRatioDistributionTrials_{cent[0]}_{cent[1]}')
    ratio_pion_pt_correlated = file_pion.Get(f'fRatioDistributionTrials_{cent[0]}_{cent[1]}')

    # get fit functions
    fit_he3 = ratio_he3.GetFunction("pol0")
    fit_hyp = ratio_hyp.GetFunction("pol0")
    fit_proton = ratio_proton.GetFunction("pol0")
    fit_pion = ratio_pion.GetFunction("pol0")

    # get ratios and errors
    ratio_he3 = fit_he3.GetParameter(0)
    ratio_he3_err = fit_he3.GetParError(0)

    ratio_hyp = fit_hyp.GetParameter(0)
    ratio_hyp_err = fit_hyp.GetParError(0)
    
    ratio_proton = fit_proton.GetParameter(0)
    ratio_proton_err = 0#fit_proton.GetParError(0)

    ratio_pion = fit_pion.GetParameter(0)
    ratio_pion_err = 0#fit_proton.GetParError(0)

    # systematic error
    syst_he3 = ratio_he3_distribution.GetRMS()
    syst_he3_abs = np.sqrt(0.00861352*0.00861352+0.00473124*0.00473124)*ratio_he3 # from absorption cross section variation
    syst_he3 = np.sqrt(syst_he3*syst_he3+syst_he3_abs*syst_he3_abs)
    syst_hyp = ratio_hyp_distribution.GetRMS()
    syst_hyp_abs = np.sqrt(0.00861352*0.00861352+0.00473124*0.00473124)*ratio_hyp # from absorption cross section variation
    syst_hyp = np.sqrt(syst_hyp*syst_hyp+syst_hyp_abs*syst_hyp_abs)
    syst_proton = fit_proton.GetParError(0)
    syst_pion = fit_pion.GetParError(0)
    syst_proton_pt_correlated = ratio_proton_pt_correlated.GetRMS()
    syst_proton_abs = np.sqrt(0.00227815*0.00227815+0.000943191*0.000943191)*ratio_proton # from absorption cross section variation
    syst_proton = np.sqrt(syst_proton*syst_proton+syst_proton_pt_correlated*syst_proton_pt_correlated+syst_proton_abs*syst_proton_abs)
    syst_pion_pt_correlated = ratio_pion_pt_correlated.GetRMS()
    syst_pion = np.sqrt(syst_pion*syst_pion+syst_pion_pt_correlated*syst_pion_pt_correlated)

    # final plot
    stat_proton = ratio_proton_err
    stat_pion = ratio_pion_err
    stat_hyp = ratio_hyp_err
    stat_he3 = ratio_he3_err
    ratios_vs_b = ROOT.TH2D(f'fRatio_vs_b_{cent[0]}_{cent[1]}', ';#it{B}+#it{S}/3;#it{I}_{3};Antimatter / Matter', 10, -0.5, 9.5, 11, -0.05, 1.05)
    ratios_vs_b.SetBinContent(10, 6, ratio_he3)
    ratios_vs_b.SetBinError(10, 6, np.sqrt(syst_he3*syst_he3+stat_he3*stat_he3))
    ratios_vs_b.SetBinContent(9, 1, ratio_hyp)
    ratios_vs_b.SetBinError(9, 1, np.sqrt(syst_hyp*syst_hyp+stat_hyp*stat_hyp))
    ratios_vs_b.SetBinContent(4, 6, ratio_proton)
    ratios_vs_b.SetBinError(4, 6, np.sqrt(syst_proton*syst_proton+stat_proton*stat_proton))
    ratios_vs_b.SetBinContent(1, 11, ratio_pion)
    ratios_vs_b.SetBinError(1, 11, np.sqrt(syst_pion*syst_pion+stat_pion*stat_pion))

    # set labels
    a = ratios_vs_b.GetXaxis()
    a.SetBit(ROOT.TAxis.kLabelsHori)
    a.SetBinLabel(1, '0')
    a.SetBinLabel(4, '1')
    a.SetBinLabel(7, '2')
    a.SetBinLabel(10, '3')
    a.SetLabelSize(0.062)

    # fit to data (exponential)
    fit_expo = ROOT.TF2(f"fit_expo_{cent[0]}_{cent[1]}", "TMath::Exp(-2./3.*[0]*x-2.*[1]*y)", -0.5, 9.5,-0.05,1.05,"")
    ratios_vs_b.Fit(f"fit_expo_{cent[0]}_{cent[1]}","SN")
    
    # chi2 and fit parameter
    formatted_chi2 = "{:.2f}".format(fit_expo.GetChisquare())
    print(f"chi2 = {formatted_chi2}")
    print(f"prob = {fit_expo.GetProb()}")
    fit_parameter_0 = fit_expo.GetParameter(0)
    fit_parameter_error_0 = fit_expo.GetParError(0)
    fit_parameter_1 = fit_expo.GetParameter(1)
    fit_parameter_error_1 = fit_expo.GetParError(1)
    print(f"mu_B / T = {fit_parameter_0} +/- {fit_parameter_error_0}; mu_I3 / T =  {fit_parameter_1} +/- {fit_parameter_error_1}")
    temperature = 155. # MeV
    print(f"mu_B (T = 155 MeV) = {fit_parameter_0*155} +/- {fit_parameter_error_0*155} +/- {fit_parameter_0*2} MeV; mu_I3 (T = 155 MeV) = {fit_parameter_1*155} +/- {fit_parameter_error_1*155} +/- {fit_parameter_1*2} MeV")
    
    # # format fit parameter and baryon chemical potential values and errors
    # formatted_fit_parameter = "{:.4f}".format(fit_parameter)
    # formatted_fit_parameter_error = "{:.4f}".format(fit_parameter_error)
    # formatted_mu_b = "{:.2f}".format(fit_parameter*155)
    # mu_b_error = np.sqrt(fit_parameter_error*fit_parameter_error/fit_parameter/fit_parameter)*fit_parameter*155
    # formatted_mu_b_error = "{:.2f}".format(fit_parameter_error*155)
    
    # # chi2 text
    # text_chi2 = ROOT.TLatex(0.5, 0.85, "#chi^{2}/NDF = "+formatted_chi2+"/"+str(fit_expo.GetNDF()))
    # text_chi2.SetTextSize(TLATEX_TEXT_SIZE)
    # text_chi2.SetTextColor(ROOT.kBlack)
    
    # # mu_b / T ratio
    # text_mu_b_over_T = ROOT.TLatex(0.5, 0.8, "#mu_{#it{B}}/#it{T} = "+formatted_fit_parameter+" #pm "+formatted_fit_parameter_error)
    # text_mu_b_over_T.SetTextSize(TLATEX_TEXT_SIZE)
    # text_mu_b_over_T.SetTextColor(ROOT.kBlack)
    
    # # mu_b at T = 155 MeV
    # formatted_temperature_error = "{:.2f}".format(fit_parameter*2)
    # text_mu_b = ROOT.TLatex(0.5, 0.75, "#mu_{#it{B}} = "+formatted_mu_b+" #pm "+formatted_mu_b_error+" #pm "+formatted_temperature_error+" MeV")
    # text_mu_b.SetTextSize(TLATEX_TEXT_SIZE)
    # text_mu_b.SetTextColor(ROOT.kBlack)

    # # T = 155 +/- 2 MeV
    # text_T = ROOT.TLatex(0.5, 0.7, "#it{T} = 155 #pm 2 MeV")
    # text_T.SetTextSize(TLATEX_TEXT_SIZE)
    # text_T.SetTextColor(ROOT.kBlack)

    # write to file
    ratios_vs_b.Write()

    # save plot
    c = ROOT.TCanvas(f"c_{cent[0]}_{cent[1]}", "c")
    c.cd()
    ratios_vs_b.SetMarkerStyle(20)
    ratios_vs_b.SetMarkerSize(0.8)
    ratios_vs_b.SetTitle(f"{cent[0]}-{cent[1]}%")
    fit_expo.SetMinimum(0.6)
    fit_expo.SetMaximum(1.23)
    fit_expo.SetNpx(10)
    fit_expo.SetNpy(10)
    #ratios_vs_b.SetMinimum(0.6)
    ratios_vs_b.SetMaximum(1.2)
    ratios_vs_b.Draw("pe")
    ratios_vs_b.GetZaxis().SetRangeUser(0.6,1.2)
    fit_expo.Draw("surf0 same")
    # ratios_vs_b.Draw("same")
    # ratios_vs_b_graph.Draw("P5 same")
    # text_chi2.Draw("same")
    # text_mu_b_over_T.Draw("same")
    # text_mu_b.Draw("same")
    # text_T.Draw("same")
    c.Write()
    fit_expo.Write()
    # c.Print(f"Ratios_{cent[0]}_{cent[1]}.pdf")

file_out.Close()
