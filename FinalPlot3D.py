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
centrality_colors = [ROOT.kOrange+7, ROOT.kAzure+4, ROOT.kTeal+4, ROOT.kCyan+1, ROOT.kMagenta-3]
particle_ratios = ["#bar{#Omega}^{+} / #Omega^{-}","#pi^{-} / #pi^{+}","#bar{p} / p","{}_{#bar{#Lambda}}^{3}#bar{H} / ^{3}_{#Lambda}H","^{3}#bar{H} / ^{3}H","^{3}#bar{He} / ^{3}He"]
n_part = [383.4,331.2,109,180,50]
n_part_err = [17.8,19.6,26.6,20,40]

apply_scaling_radius = False

TLATEX_TEXT_SIZE = 28
CENT_LIMIT_NUCLEI = 4

ROOT.gStyle.SetPadTickY(1)
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)
ROOT.gStyle.SetTextFont(44)

file_out = ROOT.TFile.Open('FinalPlot3D.root', 'recreate')
if apply_scaling_radius:
    file_out = ROOT.TFile.Open('FinalPlot3D_LHC22b9.root', 'recreate')
file_hijing = ROOT.TFile.Open('HIJINGRatios.root')

hMuBCentFrame = ROOT.TH1D("frame",";#LT#it{N}_{part}#GT;#mu_{#it{B}} (MeV)",1,0,450)
hMuICentFrame = ROOT.TH1D("frame",";#LT#it{N}_{part}#GT;#mu_{#it{Q}} (MeV)",1,0,450)
gMuBCent = ROOT.TGraphErrors()
gMuBCentCorrUnc = ROOT.TGraphErrors()
gMuICent = ROOT.TGraphErrors()
gMuICentCorrUnc = ROOT.TGraphErrors()

for i_cent, cent in enumerate(centrality_classes):
    results_muB = []
    results_muI = []
    for sign in ["Minus","Plus","STD"]:
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

        # apply material budget variation
        f_material_budget_he3 = file_material_budget.Get(f"fMHe3EffRatio{sign}_0_5")
        f_material_budget_antihe3 = file_material_budget.Get(f"fAHe3EffRatio{sign}_0_5")
        f_material_budget_proton = file_material_budget.Get(f"fMProtonEffRatio{sign}_0_5")
        f_material_budget_antiproton = file_material_budget.Get(f"fAProtonEffRatio{sign}_0_5")
        f_material_budget_proton_lowPT = file_material_budget.Get(f"fMProtonLowPTEffRatio{sign}_0_5")
        f_material_budget_antiproton_lowPT = file_material_budget.Get(f"fAProtonLowPTEffRatio{sign}_0_5")
        f_material_budget_pion = file_material_budget.Get(f"fMPionEffRatio{sign}_0_5")
        f_material_budget_antipion = file_material_budget.Get(f"fAPionEffRatio{sign}_0_5")
        f_material_budget_hyp = ROOT.TH1D()
        f_material_budget_antihyp = ROOT.TH1D()
        f_material_budget_triton = file_material_budget.Get(f"fMTritonEffRatio{sign}_0_5")
        f_material_budget_antitriton = file_material_budget.Get(f"fATritonEffRatio{sign}_0_5")
        if i_cent < CENT_LIMIT_NUCLEI:
            f_material_budget_hyp = file_material_budget_hyp.Get(f"{cent[0]}_{cent[1]}/BGBW/fCorrection_matter_{cent[0]}_{cent[1]}_BGBW")
            f_material_budget_antihyp = file_material_budget_hyp.Get(f"{cent[0]}_{cent[1]}/BGBW/fCorrection_antimatter_{cent[0]}_{cent[1]}_BGBW")
        
        # he3
        if sign!= "STD":
            for i_bins in range(1,ratio_he3.GetNbinsX()+1):
                tmp_bin_pt = ratio_he3.GetBinCenter(i_bins)
                eff_ratio_mat_he3 = f_material_budget_he3.GetBinContent(i_bins)
                eff_ratio_mat_antihe3 = f_material_budget_antihe3.GetBinContent(i_bins)
                if eff_ratio_mat_antihe3 < 1.e-6 or eff_ratio_mat_he3 < 1.e-6:
                    continue
                print(f"eff_ratio_mat_he3 = {eff_ratio_mat_he3}, eff_ratio_mat_antihe3 = {eff_ratio_mat_antihe3}")
                tmp_bin_content = ratio_he3.GetBinContent(i_bins)
                ratio_he3.SetBinContent(i_bins,tmp_bin_content*(eff_ratio_mat_he3/eff_ratio_mat_antihe3))
        ratio_he3.Fit("pol0","QR","",2.,10.)
        ratio_he3.Write()

        # triton
        if sign!= "STD":
            for i_bins in range(1,ratio_triton.GetNbinsX()+1):
                tmp_bin_pt = ratio_triton.GetBinCenter(i_bins)
                eff_ratio_mat_triton = f_material_budget_triton.GetBinContent(i_bins)
                eff_ratio_mat_antitriton = f_material_budget_antitriton.GetBinContent(i_bins)
                if eff_ratio_mat_antitriton < 1.e-6 or eff_ratio_mat_triton < 1.e-6:
                    continue
                print(f"eff_ratio_mat_triton = {eff_ratio_mat_triton}, eff_ratio_mat_antitriton = {eff_ratio_mat_antitriton}")
                tmp_bin_content = ratio_triton.GetBinContent(i_bins)
                ratio_triton.SetBinContent(i_bins,tmp_bin_content*(eff_ratio_mat_triton/eff_ratio_mat_antitriton))
        ratio_triton.Fit("pol0","QR","",1.6,3.)
        ratio_triton.Write()

        if i_cent < CENT_LIMIT_NUCLEI:
            # hyp
            # if sign!= "STD":
            #     for i_bins in range(1,ratio_hyp.GetNbinsX()+1):
            #         tmp_bin_ct = ratio_hyp.GetBinCenter(i_bins)
            #         eff_ratio_mat_hyp = f_material_budget_hyp.GetBinContent(f_material_budget_hyp.FindBin(tmp_bin_ct))
            #         eff_ratio_mat_antihyp = f_material_budget_antihyp.GetBinContent(f_material_budget_hyp.FindBin(tmp_bin_ct))
            #         tmp_bin_content = ratio_hyp.GetBinContent(i_bins)
            #         ratio_hyp.SetBinContent(i_bins,tmp_bin_content*(eff_ratio_mat_hyp/eff_ratio_mat_antihyp))
            ratio_hyp.Fit("pol0","QR","",0.,35.)
            ratio_hyp.Write()

        # proton
        if sign!= "STD":
            for i_bins in range(5,ratio_proton.GetNbinsX()+1):
                tmp_bin_pt = ratio_proton.GetBinCenter(i_bins)
                eff_ratio_mat_proton = f_material_budget_proton.GetBinContent(i_bins)
                eff_ratio_mat_antiproton = f_material_budget_antiproton.GetBinContent(i_bins)
                if eff_ratio_mat_antiproton < 1.e-6 or eff_ratio_mat_proton < 1.e-6:
                    continue
                tmp_bin_content = ratio_proton.GetBinContent(i_bins)
                print(f"eff_ratio_mat_proton = {eff_ratio_mat_proton}, eff_ratio_mat_antiproton = {eff_ratio_mat_antitriton}")
                if tmp_bin_pt<0.81:
                    eff_ratio_mat_proton = f_material_budget_proton_lowPT.GetBinContent(i_bins)
                    eff_ratio_mat_antiproton = f_material_budget_antiproton_lowPT.GetBinContent(i_bins)
                ratio_proton.SetBinContent(i_bins,tmp_bin_content*(eff_ratio_mat_proton/eff_ratio_mat_antiproton))
        ratio_proton.Fit("pol0","QR","",0.5,3.)
        ratio_proton.Write()

        # pion
        if sign!= "STD":
            for i_bins in range(1,ratio_pion.GetNbinsX()+1):
                tmp_bin_pt = ratio_pion.GetBinCenter(i_bins)
                eff_ratio_mat_pion = f_material_budget_pion.GetBinContent(i_bins)
                eff_ratio_mat_antipion = f_material_budget_antipion.GetBinContent(i_bins)
                if eff_ratio_mat_antipion < 1.e-6 or eff_ratio_mat_pion < 1.e-6:
                    continue
                tmp_bin_content = ratio_pion.GetBinContent(i_bins)
                ratio_pion.SetBinContent(i_bins,tmp_bin_content*(eff_ratio_mat_pion/eff_ratio_mat_antipion))
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

        # if apply_scaling_radius:
            # a_scaling = file_he3_scaling.Get(f"GraphA")
            # m_scaling = file_he3_scaling.Get(f"GraphM")
            # a_scaling_val = a_scaling.GetFunction("pol0").GetParameter(0)
            # m_scaling_val = m_scaling.GetFunction("pol0").GetParameter(0)
            # ratio_he3 = ratio_he3*(1+m_scaling_val-a_scaling_val)
            # ratio_hyp = 0.
            # if i_cent < CENT_LIMIT_NUCLEI:
            #     ratio_hyp = ratio_hyp*(1+m_scaling_val-a_scaling_val)
            # a_scaling = file_proton_scaling.Get(f"GraphA")
            # m_scaling = file_proton_scaling.Get(f"GraphM")
            # a_scaling_val = a_scaling.GetFunction("pol0").GetParameter(0)
            # m_scaling_val = m_scaling.GetFunction("pol0").GetParameter(0)
            # ratio_proton = ratio_proton*(1+m_scaling_val-a_scaling_val)
            # a_scaling = file_pion_scaling.Get(f"GraphA")
            # m_scaling = file_pion_scaling.Get(f"GraphM")
            # a_scaling_val = a_scaling.GetFunction("pol0").GetParameter(0)
            # m_scaling_val = m_scaling.GetFunction("pol0").GetParameter(0)
            # ratio_pion = ratio_pion*(1+m_scaling_val-a_scaling_val)

        ratio_he3_err = ratio_he3_err*ratio_he3
        ratio_triton_err = ratio_triton_err*ratio_triton
        if i_cent < CENT_LIMIT_NUCLEI:
            ratio_hyp_err = ratio_hyp_err*ratio_hyp

        correlated_err_muB = 0.
        correlated_err_muI = 0.
        if sign=="STD":
            # correlated syst
            syst_he3 = syst_he3*ratio_he3
            syst_he3_eff_prim = syst_he3_eff_prim*ratio_he3
            syst_he3_abs = np.sqrt(0.00905074*0.00905074+0.00487767*0.00487767)*ratio_he3 # from absorption cross section variation
            syst_he3 = np.sqrt(syst_he3_abs*syst_he3_abs)
            syst_triton = syst_triton*ratio_triton
            syst_triton_abs = np.sqrt(0.107459*0.107459+0.0177524*0.0177524)*ratio_triton # from absorption cross section variation
            syst_triton = np.sqrt(syst_triton_abs*syst_triton_abs)
            if i_cent < CENT_LIMIT_NUCLEI:
                syst_hyp = syst_hyp*ratio_hyp
                syst_hyp_abs = np.sqrt(0.00905074*0.00905074+0.00487767*0.00487767)*ratio_hyp # from absorption cross section variation
                syst_hyp = np.sqrt(syst_hyp_abs*syst_hyp_abs)
            syst_proton = fit_proton.GetParError(0)
            syst_pion = fit_pion.GetParError(0)
            syst_proton_pt_correlated = ratio_proton_pt_correlated.GetRMS()
            syst_proton_abs = np.sqrt(0.00129523*0.00129523+0.000173974*0.000173974)*ratio_proton # from absorption cross section variation
            syst_proton = np.sqrt(syst_proton_abs*syst_proton_abs+syst_proton_pt_correlated*syst_proton_pt_correlated)
            syst_pion_pt_correlated = ratio_pion_pt_correlated.GetRMS()
            syst_pion_abs = 0.0073795806*ratio_pion
            syst_pion = np.sqrt(syst_pion_pt_correlated*syst_pion_pt_correlated+syst_pion_abs*syst_pion_abs)

            # final plot
            stat_proton = ratio_proton_err
            stat_pion = ratio_pion_err
            stat_he3 = ratio_he3_err
            stat_hyp = 0.
            stat_triton = 0.
            stat_triton = ratio_triton_err
            if i_cent < CENT_LIMIT_NUCLEI:
                stat_hyp = ratio_hyp_err
            ratios_vs_b = ROOT.TH2D(f'fRatio_vs_b_{cent[0]}_{cent[1]}', ';#it{B}+#it{S}/3;#it{Q}-#it{S}/3;Antimatter / Matter', 10, -0.5, 9.5, 7, -0.5, 6.5)
            ratios_vs_b.SetBinContent(10, 7, ratio_he3)
            ratios_vs_b.SetBinError(10, 7, np.sqrt(syst_he3*syst_he3))
            ratios_vs_b.SetBinContent(4, 4, ratio_proton)
            ratios_vs_b.SetBinError(4, 4, np.sqrt(syst_proton*syst_proton))
            ratios_vs_b.SetBinContent(1, 4, ratio_pion)
            ratios_vs_b.SetBinError(1, 4, np.sqrt(syst_pion*syst_pion))
            ratios_vs_b.SetBinContent(10, 4, ratio_triton)
            ratios_vs_b.SetBinError(10, 4, np.sqrt(syst_triton*syst_triton))
            if i_cent < CENT_LIMIT_NUCLEI:
                ratios_vs_b.SetBinContent(9, 5, ratio_hyp)
                ratios_vs_b.SetBinError(9, 5, np.sqrt(syst_hyp*syst_hyp))

            # fit to data (exponential)
            fit_expo = ROOT.TF2(f"fit_expo_{cent[0]}_{cent[1]}", "TMath::Exp(-2./3.*[0]*x-2./3.*[1]*y)", -0.5, 9.5,0.5,6.5,"")
            r = ratios_vs_b.Fit(f"fit_expo_{cent[0]}_{cent[1]}","SN")
            fit_parameter_error_0 = fit_expo.GetParError(0)
            fit_parameter_error_1 = fit_expo.GetParError(1)
            correlated_err_muB = fit_parameter_error_0*156.2
            correlated_err_muI = fit_parameter_error_1*156.2
        
        # uncorrelated sys
        syst_he3 = syst_he3*ratio_he3
        syst_he3_eff_prim = syst_he3_eff_prim*ratio_he3
        syst_he3_abs = np.sqrt(0.00905074*0.00905074+0.00487767*0.00487767)*ratio_he3 # from absorption cross section variation
        syst_he3 = np.sqrt(syst_he3*syst_he3+syst_he3_eff_prim*syst_he3_eff_prim)
        syst_triton = syst_triton*ratio_triton
        syst_triton_abs = np.sqrt(0.107459*0.107459+0.0177524*0.0177524)*ratio_triton # from absorption cross section variation
        syst_triton = np.sqrt(syst_triton*syst_triton)
        if i_cent < CENT_LIMIT_NUCLEI:
            syst_hyp = syst_hyp*ratio_hyp
            syst_hyp_abs = np.sqrt(0.00905074*0.00905074+0.00487767*0.00487767)*ratio_hyp # from absorption cross section variation
            syst_hyp = np.sqrt(syst_hyp*syst_hyp)
        syst_proton = fit_proton.GetParError(0)
        syst_pion = fit_pion.GetParError(0)
        syst_proton_pt_correlated = ratio_proton_pt_correlated.GetRMS()
        syst_proton_abs = np.sqrt(0.00129523*0.00129523+0.000173974*0.000173974)*ratio_proton # from absorption cross section variation
        syst_proton = np.sqrt(syst_proton*syst_proton)
        syst_pion_pt_correlated = ratio_pion_pt_correlated.GetRMS()
        syst_pion = np.sqrt(syst_pion*syst_pion)

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
        aa.SetBinLabel(10, '3')
        aa.SetLabelSize(0.062)

        # fit to data (exponential)
        fit_expo = ROOT.TF2(f"fit_expo_{cent[0]}_{cent[1]}", "TMath::Exp(-2./3.*[0]*x-2./3.*[1]*y)", -0.5, 9.5,-0.5,6.5,"")
        r = ratios_vs_b.Fit(f"fit_expo_{cent[0]}_{cent[1]}","SN")
        
        # chi2 and fit parameter
        formatted_chi2 = "{:.2f}".format(fit_expo.GetChisquare())
        print(f"chi2 = {formatted_chi2}")
        print(f"prob = {fit_expo.GetProb()}")
        fit_parameter_0 = fit_expo.GetParameter(0)
        fit_parameter_error_0 = fit_expo.GetParError(0)
        fit_parameter_1 = fit_expo.GetParameter(1)
        fit_parameter_error_1 = fit_expo.GetParError(1)
        print(f"mu_B / T = {fit_parameter_0} +/- {fit_parameter_error_0}; mu_Q / T =  {fit_parameter_1} +/- {fit_parameter_error_1}")
        temperature = 156.2 # MeV
        print(f"mu_B (T = 156.2 MeV) = {fit_parameter_0*156.2} +/- {fit_parameter_error_0*156.2} +/- {fit_parameter_0*1.5} MeV; mu_Q (T = 156.2 MeV) = {fit_parameter_1*156.2} +/- {fit_parameter_error_1*156.2} +/- {fit_parameter_1*1.5} MeV")
        
        # format fit parameter and baryon chemical potential values and errors
        formatted_fit_parameter_0 = "{:.4f}".format(fit_parameter_0)
        formatted_fit_parameter_error_0 = "{:.4f}".format(fit_parameter_error_0)
        formatted_mu_b = "{:.2f}".format(fit_parameter_0*156.2)
        mu_b_error = np.sqrt(fit_parameter_error_0*fit_parameter_error_0/fit_parameter_0/fit_parameter_0)*fit_parameter_0*156.2
        formatted_mu_b_error = "{:.2f}".format(fit_parameter_error_0*156.2)

        formatted_mu_I = "{:.2f}".format(fit_parameter_1*156.2)
        mu_I_error = np.sqrt(fit_parameter_error_1*fit_parameter_error_1/fit_parameter_1/fit_parameter_1)*np.abs(fit_parameter_1)*156.2
        formatted_mu_I_error = "{:.2f}".format(np.abs(fit_parameter_error_1*156.2))
        
        # chi2 text
        text_chi2 = ROOT.TLatex(-0.27, -0.44, "#chi^{2}/NDF = "+formatted_chi2+"/"+str(fit_expo.GetNDF()))
        text_chi2.SetTextSize(TLATEX_TEXT_SIZE)
        text_chi2.SetTextColor(ROOT.kBlack)
        
        # # mu_b / T ratio
        # text_mu_b_over_T = ROOT.TLatex(0.5, 0.8, "#mu_{#it{B}}/#it{T} = "+formatted_fit_parameter+" #pm "+formatted_fit_parameter_error)
        # text_mu_b_over_T.SetTextSize(TLATEX_TEXT_SIZE)
        # text_mu_b_over_T.SetTextColor(ROOT.kBlack)

        # # T = 156.2 +/- 2 MeV
        # text_T = ROOT.TLatex(0.5, 0.7, "#it{T} = 156.2 #pm 2 MeV")
        # text_T.SetTextSize(TLATEX_TEXT_SIZE)
        # text_T.SetTextColor(ROOT.kBlack)

        # write to file
        ratios_vs_b.Write()

        if sign!="STD":
            results_muB.append(fit_parameter_0*156.2)
            results_muI.append(fit_parameter_1*156.2)
            continue

        # save plots
        error_mb = 0.5*(np.abs(results_muB[0]-results_muB[1]))
        error_mb = np.sqrt(error_mb*error_mb+correlated_err_muB*correlated_err_muB)
        error_mi = 0.5*(np.abs(results_muI[0]-results_muI[1]))
        error_mi = np.sqrt(error_mi*error_mi+correlated_err_muI*correlated_err_muI)
        print(f"* * * * * * E R R O R * * * * * *")
        print(f"MB Error = {error_mb}")
        print(f"* * * * * * * * * * * * * * * * *")
        formatted_mb_error = "{:.2f}".format(error_mb)
        formatted_mi_error = "{:.2f}".format(error_mi)
        c = ROOT.TCanvas(f"c_{cent[0]}_{cent[1]}", "c",500,500)
        m = r.GetCorrelationMatrix()
        print(f'* * * * * * * * * * * * * * * * *')
        print(f'muI - muB correlation = {m[0][1]}')
        print(f'* * * * * * * * * * * * * * * * *')
        hMat = ROOT.TH2D(f"CovMat_{cent[0]}_{cent[1]}",f"CovMat_{cent[0]}_{cent[1]}",2,0,2,2,0,2)
        for i in range(2):
            hMat.SetBinContent(i+1,i+1,m[i][i])
            print(f"i = {i},j = {i}; m = {m[i][i]}")
            hMat.SetBinContent(2-i,i+1,m[1-i][i])
            print(f"i = {1-i}, j = {i}; m = {m[1-i][i]}")
        file_out.cd()
        hMat.Write()
        
        # mu_b at T = 156.2 MeV
        formatted_temperature_error = "{:.2f}".format(fit_parameter_0*1.5)
        text_mu_b = ROOT.TLatex(-0.816381, -1.18392, "#mu_{#it{B}} = "+formatted_mu_b+" #pm "+formatted_mu_b_error+" #pm "+formatted_mb_error+" #pm "+formatted_temperature_error+" MeV")
        text_mu_b.SetTextSize(TLATEX_TEXT_SIZE)
        text_mu_b.SetTextColor(ROOT.kBlack)

        # mu_I3 at T = 156.2 MeV
        formatted_temperature_error = "{:.2f}".format(np.abs(fit_parameter_1)*1.5)
        text_mu_I = ROOT.TLatex(-0.27, -0.59, "#mu_{#it{Q}} = "+formatted_mu_I+" #pm "+formatted_mu_I_error+" #pm "+formatted_mi_error+" #pm "+formatted_temperature_error+" MeV")
        text_mu_I.SetTextSize(TLATEX_TEXT_SIZE)
        text_mu_I.SetTextColor(ROOT.kBlack)

        c.cd()
        ROOT.gPad.SetTheta(40)
        ROOT.gPad.SetPhi(40)
        ratios_vs_b.SetMarkerStyle(20)
        ratios_vs_b.SetMarkerSize(1.1)
        ratios_vs_b.GetZaxis().SetNdivisions(3)
        ratios_vs_b.GetYaxis().SetNdivisions(5)
        ratios_vs_b.SetLineColor(centrality_colors[i_cent])
        ratios_vs_b.SetLineWidth(2)
        ratios_vs_b.SetMarkerColor(centrality_colors[i_cent])
        ratios_vs_b.SetTitle(f"{cent[0]}-{cent[1]}%")
        fit_expo.SetMinimum(0.67)
        fit_expo.SetMaximum(1.13) # 1.23 to correctly visualise in file
        fit_expo.SetNpx(10)
        fit_expo.SetNpy(10)
        fit_expo.SetLineColor(ROOT.kGray+2)
        #ratios_vs_b.SetMinimum(0.6)
        ratios_vs_b.SetMaximum(1.2)
        ratios_vs_b.Draw("pe")
        ratios_vs_b.GetZaxis().SetRangeUser(0.67,1.13)
        fit_expo.Draw("surf0 same")
        ratios_vs_b.Draw("pe same")
        # ratios_vs_b_graph.Draw("P5 same")
        text_chi2.Draw("same")
        # text_mu_b_over_T.Draw("same")
        text_mu_b.Draw("same")
        # text_T.Draw("same")
        c.Write()
        fit_expo.Write()
        c.Print(f"Ratios_{cent[0]}_{cent[1]}_3D.pdf")

        # make final plot for approval (preliminary)
        leg_ratios_particle = ROOT.TLegend(0.72,0.1,0.92,0.3)
        leg_ratios_particle.SetTextFont(43)
        leg_ratios_particle.SetTextSize(22)
        cRatiosParticle = ROOT.TCanvas(f"cRatiosParticle_{cent[0]}_{cent[1]}",f"cRatiosParticle_{cent[0]}_{cent[1]}")
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
            ratio = 1
            ratio_err = 0
            fit = 1
            if i_part == 0:
                ratio = ratio_omega
                ratio_err = np.sqrt(syst_omega*syst_omega+ratio_omega_err*ratio_omega_err)
                fit = fit_expo.Eval(0,0)
            if i_part == 1:
                ratio = ratio_pion
                ratio_err = np.sqrt(syst_pion*syst_pion+stat_pion*stat_pion)
                fit = fit_expo.Eval(0,3)
            if i_part == 2:
                ratio = ratio_proton
                ratio_err = np.sqrt(syst_proton*syst_proton+stat_proton*stat_proton)
                fit = fit_expo.Eval(3,3)
            if i_part == 3 and i_cent < CENT_LIMIT_NUCLEI:
                ratio = ratio_hyp
                ratio_err = np.sqrt(syst_hyp*syst_hyp+stat_hyp*stat_hyp)
                fit = fit_expo.Eval(8,4)
            if i_part == 4:
                ratio = ratio_triton
                ratio_err = np.sqrt(syst_triton*syst_triton+stat_triton*stat_triton)
                fit = fit_expo.Eval(9,3)
            if i_part == 5:
                ratio = ratio_he3
                ratio_err = np.sqrt(syst_he3*syst_he3+stat_he3*stat_he3)
                fit = fit_expo.Eval(9,6)
            hRatiosParticle.SetBinContent(i_part+1,ratio)
            hRatiosParticle.SetBinError(i_part+1,ratio_err)
            hRatiosParticleFit.SetBinContent(i_part+1,fit)
            hRatiosParticleFit.SetBinError(i_part+1,0)
            print(f"fit = {fit}")
        hRatiosParticle.GetYaxis().SetRangeUser(0.67,1.18)
        gRatiosParticle = ROOT.TGraphErrors(hRatiosParticle)
        gRatiosParticle.SetName(f"gRatiosParticle_{cent[0]}_{cent[1]}")
        gRatiosParticleFit = ROOT.TGraphErrors(hRatiosParticleFit)
        gRatiosParticleFit.SetName(f"gRatiosParticleFit_{cent[0]}_{cent[1]}")
        for i_part in range(0,6):
            gRatiosParticle.SetPointError(i_part,0,hRatiosParticle.GetBinError(i_part+1))
            gRatiosParticleFit.SetPointError(i_part,0.3,0)
        # hijing_predictions = file_hijing.Get(f"fRatios_{cent[0]}_{cent[1]}")
        # gHijingPredictions = ROOT.TGraphErrors()
        # gHijingPredictions.AddPoint(gRatiosParticle.GetPointX(0),hijing_predictions.GetBinContent(1,21)) # pion
        # gHijingPredictions.SetPointError(0,0.25,hijing_predictions.GetBinError(1,21)) # pion
        # gHijingPredictions.AddPoint(gRatiosParticle.GetPointX(1),hijing_predictions.GetBinContent(4,16)) # proton
        # gHijingPredictions.SetPointError(1,0.25,hijing_predictions.GetBinError(4,16)) # pion

        # data to fit ratio
        for i_ticks in range(6):
            hNSigmaRatioFitParticle.GetYaxis().SetNdivisions(5)
        for i_part in range(0,6):
            hNSigmaRatioFitParticle.GetXaxis().SetBinLabel(i_part+1,particle_ratios[i_part])
            hNSigmaRatioFitParticle.GetXaxis().SetLabelSize(0.2)
            ratio = -999.
            ratio_err = -999.
            fit = -999.
            if i_part == 0:
                ratio = ratio_omega
                ratio_err = np.sqrt(syst_omega*syst_omega+ratio_omega_err*ratio_omega_err)
                fit = fit_expo.Eval(0,0)
            if i_part == 1:
                ratio = ratio_pion
                ratio_err = np.sqrt(syst_pion*syst_pion+stat_pion*stat_pion)
                fit = fit_expo.Eval(0,3)
            if i_part == 2:
                ratio = ratio_proton
                ratio_err = np.sqrt(syst_proton*syst_proton+stat_proton*stat_proton)
                fit = fit_expo.Eval(3,3)
            if i_part == 3 and i_cent < CENT_LIMIT_NUCLEI:
                ratio = ratio_hyp
                ratio_err = np.sqrt(syst_hyp*syst_hyp+stat_hyp*stat_hyp)
                fit = fit_expo.Eval(8,4)
            if i_part == 4:
                ratio = ratio_triton
                ratio_err = np.sqrt(syst_triton*syst_triton+stat_triton*stat_triton)
                fit = fit_expo.Eval(9,3)
            if i_part == 5:
                ratio = ratio_he3
                ratio_err = np.sqrt(syst_he3*syst_he3+stat_he3*stat_he3)
                fit = fit_expo.Eval(9,6)
            if ratio > 0:
                hNSigmaRatioFitParticle.SetBinContent(i_part+1,(ratio-fit)/ratio_err)
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
        hNSigmaRatioFitParticle.GetYaxis().SetRangeUser(-2.3,2.3)
        # gHijingNsigma = ROOT.TGraphErrors(gHijingPredictions)
        # gHijingNsigma.SetPointY(0,(gHijingNsigma.GetPointY(0)-fit_expo.Eval(0,1))/np.sqrt(syst_pion*syst_pion+stat_pion*stat_pion))
        # gHijingNsigma.SetPointY(1,(gHijingNsigma.GetPointY(1)-fit_expo.Eval(3,0.5))/np.sqrt(syst_proton*syst_proton+stat_proton*stat_proton))
        hLineZero = ROOT.TLine(0.,0.,5.,0.)
        hLineZero.SetLineColor(ROOT.kBlack)
        hLineZero.SetLineStyle(ROOT.kDashed)
        hLineOnePos = ROOT.TLine(0.,1.,5.,1.)
        hLineOnePos.SetLineColor(ROOT.kBlack)
        hLineOnePos.SetLineStyle(ROOT.kDashed)
        hLineOneNeg = ROOT.TLine(0.,-1.,5.,-1.)
        hLineOneNeg.SetLineColor(ROOT.kBlack)
        hLineOneNeg.SetLineStyle(ROOT.kDashed)
        hLineTwoPos = ROOT.TLine(0.,2.,5.,2.)
        hLineTwoPos.SetLineColor(ROOT.kBlack)
        hLineTwoPos.SetLineStyle(ROOT.kDashed)
        hLineTwoNeg = ROOT.TLine(0.,-2.,5.,-2.)
        hLineTwoNeg.SetLineColor(ROOT.kBlack)
        hLineTwoNeg.SetLineStyle(ROOT.kDashed)

        hRatiosParticle.GetYaxis().SetRangeUser(0.67,1.18)
        text_alice = ROOT.TLatex(0.18,0.3,"ALICE")
        text_alice.SetNDC(True)
        text_alice.SetTextFont(43)
        text_alice.SetTextSize(28)
        text_energy = ROOT.TLatex(0.18,0.2,"Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV")
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
        # leg_ratios_particle.AddEntry(gHijingPredictions,"HIJING")
        text_mu_b.SetTextFont(43)
        text_mu_I.SetTextFont(43)
        text_chi2.SetTextFont(43)
        text_mu_b.SetTextSize(20)
        text_mu_I.SetTextSize(20)
        text_chi2.SetTextSize(22)
        text_mu_b.SetNDC(True)
        text_mu_I.SetNDC(True)
        text_chi2.SetNDC(True)
        text_mu_b.SetX(0.52)
        text_mu_b.SetY(0.68)
        text_mu_I.SetX(0.52)
        text_mu_I.SetY(0.75)
        text_chi2.SetX(0.18)
        text_chi2.SetY(0.1)
        # text_mu_b.Draw("same")
        # text_mu_I.Draw("same")
        text_chi2.Draw("same")
        text_alice.Draw("same")
        text_energy.Draw("same")
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
        cRatiosParticle.Print(f"SHMfitToRatios_{cent[0]}_{cent[1]}.eps")
        gRatiosParticle.Write()
        gRatiosParticleFit.Write()
        hNSigmaRatioFitParticle.Write()

        gMuBCent.AddPoint(n_part[i_cent],fit_parameter_0*156.2)
        gMuBCent.SetPointError(gMuBCent.GetN()-1,n_part_err[i_cent],mu_b_error)
        gMuBCentCorrUnc.AddPoint(n_part[i_cent],fit_parameter_0*156.2)
        gMuBCentCorrUnc.SetPointError(gMuBCent.GetN()-1,10,np.abs(error_mb))
        gMuICent.AddPoint(n_part[i_cent],fit_parameter_1*156.2)
        gMuICent.SetPointError(gMuBCent.GetN()-1,0.,mu_I_error)
        gMuICentCorrUnc.AddPoint(n_part[i_cent],fit_parameter_1*156.2)
        gMuICentCorrUnc.SetPointError(gMuICent.GetN()-1,10,np.abs(error_mi))
        print(f'mu_I_error = {mu_I_error}')

gMuBCent.SetName("MuBCent")
gMuICent.SetName("MuICent")
gMuBCentCorrUnc.SetName("gMuBCentCorrUnc")
gMuICentCorrUnc.SetName("gMuICentCorrUnc")

legMuBCent = ROOT.TLegend(0.18,0.15,0.60,0.28)
gPublishedResult = ROOT.TGraphErrors()
gPublishedResult.AddPoint(354.7,0.7)
gPublishedResult.SetPointError(0,33,3.8)
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
gHijingMuB = file_hijing.Get("gMuBCent")
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
# gMuBCent.Fit("pol0")
# gMuBCent.GetFunction("pol0").SetLineColor(ROOT.kGray+2)
# gMuBCent.GetFunction("pol0").SetLineStyle(ROOT.kDashed)
gPublishedResult.Draw("pe5same")
gMuBCentCorrUnc.Draw("pe5same")
gMuBCent.Draw("pe5same")
#gHijingMuB.Draw("e3same")
legMuBCent.SetNColumns(2)
legMuBCent.SetTextFont(43)
legMuBCent.SetTextSize(22)
legMuBCent.AddEntry(gMuBCent,"Uncorr. uncert.", "pf")
legMuBCent.AddEntry(gMuBCentCorrUnc,"Corr. uncert.")
legMuBCent.AddEntry(gPublishedResult,"SHM fit, #it{Nature} #bf{561}, 321-330 (2018)")
legMuBCent.Draw("same")
text_alice.Draw("same")
text_energy.Draw("same")
#text_uncorr_uncertainty.Draw("same")
cMuBCent.Write()
cMuBCent.Print("MuBvsNpart.eps")

legMuICent = ROOT.TLegend(0.18, 0.165, 0.6, 0.228)
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
cMuICent = ROOT.TCanvas("cMuICent","cMuICent")
gHijingMuI = file_hijing.Get("gMuICent")
cMuICent.cd()
gMuICent.SetMarkerStyle(20)
gMuICent.SetMarkerSize(1.2)
gMuICent.SetMarkerColor(ROOT.kRed)
gMuICent.SetLineColor(ROOT.kRed)
gMuICent.SetFillColor(0)
gMuICent.SetFillStyle(0)
gMuICentCorrUnc.SetMarkerColor(ROOT.kRed-9)
gMuICentCorrUnc.SetLineWidth(0)
gMuICentCorrUnc.SetMarkerSize(0)
gMuICentCorrUnc.SetLineColor(0)
gMuICentCorrUnc.SetFillColor(ROOT.kRed-9)
hMuICentFrame.GetYaxis().SetRangeUser(-0.3,0.5)
hMuICentFrame.GetYaxis().SetNdivisions(510)
hMuICentFrame.GetYaxis().SetDecimals()
hMuICentFrame.GetXaxis().SetDecimals()
hMuICentFrame.Draw("axis")
# gMuICent.Fit("pol0")
# gMuICent.GetFunction("pol0").SetLineColor(ROOT.kGray+2)
# gMuICent.GetFunction("pol0").SetLineStyle(ROOT.kDashed)
gMuICentCorrUnc.Draw("pe5same")
gMuICent.Draw("pe5same")
#gHijingMuI.Draw("e3same")
legMuICent.SetNColumns(2)
legMuICent.SetTextFont(43)
legMuICent.SetTextSize(22)
legMuICent.AddEntry(gMuICent,"Uncorr. uncert.", "pf")
legMuICent.AddEntry(gMuICentCorrUnc,"Corr. uncert.")
legMuICent.Draw("same")
text_alice.Draw("same")
text_energy.Draw("same")
#legMuICent.Draw("same")
cMuICent.Write()
cMuICent.Print("MuIvsCent.pdf")

gMuBCent.Write()
gMuICent.Write()

file_out.Close()
