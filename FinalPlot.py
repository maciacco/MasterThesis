#!/usr/bin/env python3
import ROOT
import numpy as np

path_he3 = './He3_PbPb/out'
path_hyp = './Hypertriton_PbPb'
path_proton = './Proton_PbPb/out'
centrality_classes = [[0, 5], [5, 10], [30, 50]]
centrality_colors = [ROOT.kOrange+7, ROOT.kAzure+4, ROOT.kTeal+4]

TLATEX_TEXT_SIZE = 28

ROOT.gStyle.SetPadTickY(1)
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetTextFont(44)

file_he3 = ROOT.TFile.Open(path_he3 + '/SpectraHe3.root')
file_hyp = ROOT.TFile.Open(path_hyp + '/Ratio.root')
file_proton = ROOT.TFile.Open(path_proton + '/SpectraProton_MC21l5_raw.root')
file_he3_syst = ROOT.TFile.Open(path_he3 + '/SystematicsAll.root')
file_he3_syst_abs = ROOT.TFile.Open(path_he3 + '/AbsError.root')
file_hyp_syst = ROOT.TFile.Open(path_hyp + '/Systematics.root')
file_hyp_syst_abs = ROOT.TFile.Open(path_hyp + '/AbsError.root')
file_proton_syst = ROOT.TFile.Open(path_proton + '/SystematicsAll.root')
file_proton_syst_abs = ROOT.TFile.Open(path_proton + '/AbsErrorMCorrection.root')

file_out = ROOT.TFile.Open('FinalPlot.root', 'recreate')

for i_cent, cent in enumerate(centrality_classes):

    # get histograms
    ratio_he3 = file_he3.Get(f'1.0_89_1_1_1/fRatio_{cent[0]}_{cent[1]}')
    #ratio_he3_abs = file_he3.Get(f'1.0_89_1_1_1/fRatio_{cent[0]}_{cent[1]}')
    ratio_hyp = file_hyp.Get(f'fRatio_{cent[0]}_{cent[1]}')
    ratio_proton = file_proton.Get(f'fRatio_{cent[0]}_{cent[1]}')
    ratio_he3_distribution = file_he3_syst.Get(f'hist/fFitPar_{cent[0]}_{cent[1]}')
    ratio_he3_distribution_abs = file_he3_syst_abs.Get(f'fFitPar_{cent[0]}_{cent[1]}')
    ratio_hyp_distribution = file_hyp_syst.Get(f'fParameterDistribution_{cent[0]}_{cent[1]}')
    ratio_hyp_distribution_abs = file_hyp_syst_abs.Get(f'fParameterDistribution_{cent[0]}_{cent[1]}')
    ratio_proton_distribution = file_proton_syst.Get(f'hist/fFitPar_{cent[0]}_{cent[1]}')
    ratio_proton_distribution_abs = file_proton_syst_abs.Get(f'fFitPar_{cent[0]}_{cent[1]}')

    # get fit functions
    fit_he3 = ratio_he3.GetFunction("pol0")
    fit_hyp = ratio_hyp.GetFunction("pol0")
    fit_proton = ratio_proton.GetFunction("pol0")

    # get ratios and errors
    ratio_he3 = fit_he3.GetParameter(0)
    ratio_he3_err = fit_he3.GetParError(0)

    ratio_hyp = fit_hyp.GetParameter(0)
    ratio_hyp_err = fit_hyp.GetParError(0)
    
    ratio_proton = fit_proton.GetParameter(0)
    ratio_proton_err = 0#fit_proton.GetParError(0)

    # systematic error
    syst_he3 = ratio_he3_distribution.GetRMS()
    syst_he3_abs = ratio_he3_distribution_abs.GetRMS()
    #syst_he3 = np.sqrt(syst_he3*syst_he3+syst_he3_abs*syst_he3_abs)
    syst_hyp = ratio_hyp_distribution.GetRMS()
    syst_hyp_abs = ratio_hyp_distribution_abs.GetRMS()
    #syst_hyp = np.sqrt(syst_hyp*syst_hyp+syst_hyp_abs*syst_hyp_abs)
    syst_proton = fit_proton.GetParError(0)#ratio_proton_distribution.GetRMS()
    syst_proton_abs = ratio_proton_distribution_abs.GetRMS()
    #syst_proton = np.sqrt(syst_proton*syst_proton+syst_proton_abs*syst_proton_abs)

    # final plot
    ratios_vs_b = ROOT.TH1D(f'fRatio_vs_b_{cent[0]}_{cent[1]}', ';B+S/3; Antimatter / Matter', 10, -0.5, 9.5)
    ratios_vs_b.SetBinContent(10, ratio_he3)
    ratios_vs_b.SetBinError(10, ratio_he3_err)
    ratios_vs_b.SetBinContent(9, ratio_hyp)
    ratios_vs_b.SetBinError(9, ratio_hyp_err)
    ratios_vs_b.SetBinContent(4, ratio_proton)
    ratios_vs_b.SetBinError(4, ratio_proton_err)

    # set labels
    a = ratios_vs_b.GetXaxis()
    a.SetBit(ROOT.TAxis.kLabelsHori)
    a.SetBinLabel(1, '0')
    a.SetBinLabel(4, '1')
    a.SetBinLabel(7, '2')
    a.SetBinLabel(10, '3')
    a.SetLabelSize(0.062)

    # set histogram style
    ratios_vs_b.SetTitle(f'{cent[0]}-{cent[1]}%')
    ratios_vs_b.SetMarkerColor(centrality_colors[i_cent])
    ratios_vs_b.SetLineColor(centrality_colors[i_cent])
    ratios_vs_b.SetMarkerStyle(20)
    ratios_vs_b.SetMarkerSize(0.8)
    ratios_vs_b.GetYaxis().SetRangeUser(0.6, 1.2)

    ratios_vs_b_graph = ROOT.TGraphErrors(ratios_vs_b)
    for _ in range(3):
        ratios_vs_b_graph.RemovePoint(0)
    for _ in range(4):
        ratios_vs_b_graph.RemovePoint(1)

    # set systematic errors
    ratios_vs_b_graph.SetPointError(1, 0.5, syst_hyp)
    ratios_vs_b_graph.SetPointError(2, 0.5, syst_he3)
    ratios_vs_b_graph.SetPointError(0, 0.5, syst_proton)

    ratios_vs_b_graph.SetFillColor(ROOT.kWhite)
    ratios_vs_b_graph.SetFillStyle(3000)
    ratios_vs_b_graph.SetLineColor(centrality_colors[i_cent])
    ratios_vs_b_graph.SetDrawOption("AP5")

    # define histogram for fit
    ratios_vs_b_fit = ROOT.TH1D(ratios_vs_b)
    stat_proton = ratios_vs_b.GetBinError(4)
    stat_hyp = ratios_vs_b.GetBinError(9)
    stat_he3 = ratios_vs_b.GetBinError(10)
    ratios_vs_b_fit.SetBinError(4, np.sqrt(syst_proton*syst_proton+stat_proton*stat_proton))
    ratios_vs_b_fit.SetBinError(9, np.sqrt(syst_hyp*syst_hyp+stat_hyp*stat_hyp))
    ratios_vs_b_fit.SetBinError(10, np.sqrt(syst_he3*syst_he3+stat_he3*stat_he3))

    # fit to data (exponential)
    fit_expo = ROOT.TF1(f"fit_expo_{cent[0]}_{cent[1]}", "TMath::Exp(-2./3.*[0]*x)", -0.5, 9.5)
    fit_expo.SetNpx(10000)
    ratios_vs_b_fit.Fit(f"fit_expo_{cent[0]}_{cent[1]}")#,"R","",6.5,9.5)
    
    # chi2 and fit parameter
    formatted_chi2 = "{:.2f}".format(fit_expo.GetChisquare())
    print(f"chi2 = {formatted_chi2}")
    print(f"prob = {fit_expo.GetProb()}")
    fit_parameter = fit_expo.GetParameter(0)
    fit_parameter_error = fit_expo.GetParError(0)
    print(f"mu_B / T = {fit_parameter} +/- {fit_parameter_error}")
    temperature = 155. # MeV
    print(f"mu_B (T = 155 MeV) = {fit_parameter*155} +/- {fit_parameter_error*155}")
    
    # format fit parameter and baryon chemical potential values and errors
    formatted_fit_parameter = "{:.4f}".format(fit_parameter)
    formatted_fit_parameter_error = "{:.4f}".format(fit_parameter_error)
    formatted_mu_b = "{:.2f}".format(fit_parameter*155)
    mu_b_error = np.sqrt(fit_parameter_error*fit_parameter_error/fit_parameter/fit_parameter+2*2/155/155)*fit_parameter*155
    formatted_mu_b_error = "{:.2f}".format(fit_parameter_error*155)
    
    # chi2 text
    text_chi2 = ROOT.TLatex(0.5, 0.85, "#chi^{2}/NDF = "+formatted_chi2+"/"+str(fit_expo.GetNDF()))
    text_chi2.SetTextSize(TLATEX_TEXT_SIZE)
    text_chi2.SetTextColor(ROOT.kBlack)
    
    # mu_b / T ratio
    text_mu_b_over_T = ROOT.TLatex(0.5, 0.8, "#mu_{#it{B}}/#it{T} = "+formatted_fit_parameter+" #pm "+formatted_fit_parameter_error)
    text_mu_b_over_T.SetTextSize(TLATEX_TEXT_SIZE)
    text_mu_b_over_T.SetTextColor(ROOT.kBlack)
    
    # mu_b at T = 155 MeV
    text_mu_b = ROOT.TLatex(0.5, 0.75, "#mu_{#it{B}} = "+formatted_mu_b+" #pm "+formatted_mu_b_error+" MeV, #it{T} = 155 MeV")
    text_mu_b.SetTextSize(TLATEX_TEXT_SIZE)
    text_mu_b.SetTextColor(ROOT.kBlack)

    # write to file
    ratios_vs_b.Write()

    # save plot
    c = ROOT.TCanvas(f"c_{cent[0]}_{cent[1]}", "c")
    c.cd()
    ratios_vs_b.Draw("pe")
    fit_expo.Draw("same")
    ratios_vs_b.Draw("same")
    ratios_vs_b_graph.Draw("P5 same")
    text_chi2.Draw("same")
    text_mu_b_over_T.Draw("same")
    text_mu_b.Draw("same")
    c.Write()
    fit_expo.Write()
    c.Print(f"Ratios_{cent[0]}_{cent[1]}.pdf")

file_out.Close()
