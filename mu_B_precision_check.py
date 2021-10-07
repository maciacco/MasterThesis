#!/usr/bin/env python3
import ROOT
import numpy as np

N_TRIALS = 1000
N_UNCERTAINTIES = 150

path_he3 = './He3_PbPb/out'
path_hyp = './Hypertriton_PbPb'
centrality_classes = [[0, 5], [5, 10], [30, 50]]
centrality_colors = [ROOT.kOrange+7, ROOT.kAzure+7, ROOT.kTeal+4]

ROOT.gStyle.SetPadTickY(1)
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(1111)
ROOT.gStyle.SetTextFont(4)

file_he3 = ROOT.TFile.Open(path_he3 + '/SpectraHe3.root')
file_hyp = ROOT.TFile.Open(path_hyp + '/Ratio.root')
file_he3_syst = ROOT.TFile.Open(path_he3 + '/SystematicsAll.root')
file_hyp_syst = ROOT.TFile.Open(path_hyp + '/Systematics.root')

file_out = ROOT.TFile.Open('MuBPrecisionCheck.root', 'recreate')

for i_cent, cent in enumerate(centrality_classes):

    # get histograms
    ratio_he3 = file_he3.Get(f'1.0_89_1_1_1/fRatio_{cent[0]}_{cent[1]}')
    ratio_hyp = file_hyp.Get(f'fRatio_{cent[0]}_{cent[1]}')
    ratio_he3_distribution = file_he3_syst.Get(f'hist/fFitPar_{cent[0]}_{cent[1]}')
    ratio_hyp_distribution = file_hyp_syst.Get(f'fParameterDistribution_{cent[0]}_{cent[1]}')

    # get fit functions
    fit_he3 = ratio_he3.GetFunction("pol0")
    fit_hyp = ratio_hyp.GetFunction("pol0")

    # get ratios and errors
    ratio_he3 = fit_he3.GetParameter(0)
    ratio_he3_err = fit_he3.GetParError(0)

    ratio_hyp = fit_hyp.GetParameter(0)
    ratio_hyp_err = fit_hyp.GetParError(0)

    # systematic error
    syst_he3 = ratio_he3_distribution.GetRMS()
    syst_hyp = ratio_hyp_distribution.GetRMS()

    # final plot
    ratios_vs_b = ROOT.TH1D(f'fRatio_vs_b_{cent[0]}_{cent[1]}', ';B+S/3; Antimatter / Matter', 10, -0.5, 9.5)
    ratios_vs_b.SetBinContent(10, ratio_he3)
    ratios_vs_b.SetBinError(10, ratio_he3_err)
    ratios_vs_b.SetBinContent(9, ratio_hyp)
    ratios_vs_b.SetBinError(9, ratio_hyp_err)

    # set labels
    a = ratios_vs_b.GetXaxis()
    a.SetBit(ROOT.TAxis.kLabelsHori)
    a.SetBinLabel(1, '0')
    a.SetBinLabel(4, '1')
    a.SetBinLabel(7, '2')
    a.SetBinLabel(10, '3')
    a.SetLabelSize(0.062)

    # define histogram for fit
    ratios_vs_b_fit = ROOT.TH1D(ratios_vs_b)
    stat_hyp = ratios_vs_b.GetBinError(9)
    stat_he3 = ratios_vs_b.GetBinError(10)
    ratios_vs_b_fit.SetBinError(9, np.sqrt(syst_hyp*syst_hyp+stat_hyp*stat_hyp))
    ratios_vs_b_fit.SetBinError(10, np.sqrt(syst_he3*syst_he3+stat_he3*stat_he3))

    # fit to data (exponential)
    fit_expo = ROOT.TF1(f"fit_expo_{cent[0]}_{cent[1]}", "TMath::Exp(-2./3.*[0]*x)", -0.5, 9.5)
    fit_expo.SetNpx(10000)
    ratios_vs_b_fit.Fit(f"fit_expo_{cent[0]}_{cent[1]}", "Q")
    
    # chi2 and fit parameter
    formatted_chi2 = "{:.2f}".format(fit_expo.GetChisquare())
    print(f"chi2 = {formatted_chi2}")
    fit_parameter = fit_expo.GetParameter(0)
    fit_parameter_error = fit_expo.GetParError(0)
    print(f"mu_B / T = {fit_parameter} +/- {fit_parameter_error}")
    temperature = 155. # MeV
    print(f"mu_B (T = 155 MeV) = {fit_parameter*155} +/- {fit_parameter_error*155}")

    # get proton ratio expected from exponential fit
    proton_ratio_fit = fit_expo.Eval(ratios_vs_b_fit.GetBinCenter(4))
    print(f"Proton ratio from fit = {proton_ratio_fit}")

    # mu_b histogram
    h_mu_b = ROOT.TH1D(f"mu_b_{cent[0]}_{cent[1]}", ";#sigma(^{3} #bar{He}/^{3} He) / #sigma(#bar{p}/p);#mu_{B} (MeV)",150, 0.5, 75.5)
    h_mu_b.SetTitle(f"{cent[0]}-{cent[1]}%")
    
    # mu_b uncertainty histogram
    h_mu_b_sigma = ROOT.TH1D(f"mu_b_sigma_{cent[0]}_{cent[1]}", ";#sigma(^{3} #bar{He}/^{3} He) / #sigma(#bar{p}/p);#sigma(#mu_{B}) (MeV)",150, 0.5, 75.5)
    h_mu_b_sigma.SetTitle(f"{cent[0]}-{cent[1]}%")

    # mu_b uncertainty histogram
    h_mu_b_uncertainty = ROOT.TH1D(f"mu_b_uncertainty_{cent[0]}_{cent[1]}", ";#sigma(^{3}#bar{He}/^{3}He) / #sigma(#bar{p}/p);#sigma(#mu_{B}) (MeV)",150, 0.5, 75.5)
    h_mu_b_uncertainty.SetTitle(f"{cent[0]}-{cent[1]}%")

    for i_uncertainties in range(N_UNCERTAINTIES):

        # define mu_b histogram
        h_mu_b_tmp = ROOT.TH1D(f"mu_b_tmp_n_{i_uncertainties+1}_{cent[0]}_{cent[1]}", ";#mu_{B} (MeV);Entries",400,-4.,4.)
        h_mu_b_tmp.SetDrawOption("pe")
        h_mu_b_tmp_error = ROOT.TH1D(f"mu_b_tmp_error_n_{i_uncertainties+1}_{cent[0]}_{cent[1]}", ";#sigma(#mu_{B}) (MeV);Entries",24000,0.,0.8)
        h_mu_b_tmp_error.SetDrawOption("pe")
        
        # deuteron ratio uncertainty (starting from that of he3)
        proton_ratio_sigma = ratios_vs_b_fit.GetBinError(10)*2./(i_uncertainties+1)
        print(f"Proton ratio uncertainty = {proton_ratio_sigma}; scale factor = {i_uncertainties}")
        
        # make direcotry in output file
        file_out.mkdir(f"{i_uncertainties+1}_scale_factor")
        file_out.cd(f"{i_uncertainties+1}_scale_factor")
        
        for i_trial in range(N_TRIALS):
            ratios_vs_b_fit.SetName(f"ratio_vs_b_fit_trial_{i_trial}_{cent[0]}_{cent[1]}")
            
            # generate deuteron ratio measurement
            proton_ratio_gen = ROOT.gRandom.Gaus(proton_ratio_fit, proton_ratio_sigma)
            # print(f"generated deuteron ratio = {proton_ratio_gen}")

            # set generated points in histogram
            ratios_vs_b_fit.SetBinContent(4, proton_ratio_gen)
            ratios_vs_b_fit.SetBinError(4, proton_ratio_sigma)
            
            # fit to data
            fit_expo_tmp = ROOT.TF1(f"fit_expo_tmp_{cent[0]}_{cent[1]}", "TMath::Exp(-2./3.*[0]*x)", -0.5, 9.5)
            ratios_vs_b_fit.Fit(f"fit_expo_tmp_{cent[0]}_{cent[1]}","Q")
            
            # measure mu_b
            fit_parameter_tmp = fit_expo_tmp.GetParameter(0)
            fit_parameter_tmp_error = fit_expo_tmp.GetParError(0)
            mu_b_tmp = fit_parameter_tmp*155
            mu_b_tmp_error = fit_parameter_tmp_error*155
            # print(f"mu_B (T = 155 MeV) = {mu_b_tmp} +/- {fit_parameter_tmp_error*155}")
        
            # fill mu_b histogram
            if fit_expo_tmp.GetProb() < 0.975 and fit_expo_tmp.GetProb() > 0.025:
                h_mu_b_tmp.Fill(mu_b_tmp)
                h_mu_b_tmp_error.Fill(mu_b_tmp_error)
                if (i_trial % 1000 == 0) :
                    ratios_vs_b_fit.Write()
        
        # get mean and RMS from gaussian fit
        h_mu_b_tmp.Fit("gaus","Q")
        mu_b_mean = h_mu_b_tmp.GetFunction("gaus").GetParameter(1)
        mu_b_sigma = h_mu_b_tmp.GetFunction("gaus").GetParameter(2)
        mu_b_sigma_error = h_mu_b_tmp.GetFunction("gaus").GetParError(2)
        h_mu_b_tmp_error.Fit("gaus","Q")
        mu_b_error = h_mu_b_tmp_error.GetFunction("gaus").GetParameter(1)
        mu_b_error_error = h_mu_b_tmp_error.GetFunction("gaus").GetParError(1)
        
        h_mu_b.SetBinContent(i_uncertainties+1, mu_b_mean)
        h_mu_b.SetBinError(i_uncertainties+1, mu_b_error)
        h_mu_b_sigma.SetBinContent(i_uncertainties+1, mu_b_sigma)
        h_mu_b_sigma.SetBinError(i_uncertainties+1, mu_b_sigma_error)
        h_mu_b_uncertainty.SetBinContent(i_uncertainties+1, mu_b_error)
        h_mu_b_uncertainty.SetBinError(i_uncertainties+1, mu_b_error_error)
        
        h_mu_b_tmp.Write()
        h_mu_b_tmp_error.Write()
        
    file_out.cd()    
    h_mu_b.Write()
    h_mu_b_sigma.Write()
    h_mu_b_uncertainty.Write()
    
    # write to PNG
    
    # write png (uncertainty)
    c = ROOT.TCanvas("c", "c")
    h_mu_b_uncertainty.SetMarkerStyle(20)
    h_mu_b_uncertainty.SetMarkerSize(0.8)
    h_mu_b_uncertainty.Draw()
    print_name = "./"+h_mu_b_uncertainty.GetName()+".pdf"
    c.Print(print_name)   

file_out.Close()