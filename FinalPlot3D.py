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
particle_ratios = ["#pi^{-} / #pi^{+}","#bar{p} / p","{}_{#bar{#Lambda}}^{3}#bar{H} / ^{3}_{#Lambda}H","^{3}#bar{He} / ^{3}He"]
n_part = [383.4,331.2,109]
n_part_err = [17.8,19.6,26.6]

apply_scaling_radius = False

TLATEX_TEXT_SIZE = 28

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
hMuICentFrame = ROOT.TH1D("frame",";Centrality (%);#mu_{#it{I}_{3}} (MeV)",1,0,50)
gMuBCent = ROOT.TGraphErrors()
gMuBCentCorrUnc = ROOT.TGraphErrors()
gMuICent = ROOT.TGraphErrors()

for i_cent, cent in enumerate(centrality_classes):
    results_muB = []
    results_muI = []
    for sign in [-1.,1.,0.]:
        # get histograms
        file_material_budget = ROOT.TFile.Open("MaterialBudgetUncertainty.root")
        file_material_budget_hyp = ROOT.TFile.Open("MaterialBudgetUncertaintyHyp.root")
        file_he3 = ROOT.TFile.Open(path_he3 + '/SpectraHe3.root')
        file_hyp = ROOT.TFile.Open(path_hyp + '/Ratio.root')
        file_proton = ROOT.TFile.Open(path_proton + '/SystematicsAllEPtNotCombinedTOFTEST.root')
        file_pion = ROOT.TFile.Open(path_pion + '/SystematicsAllEPtNotCombined.root')
        file_he3_syst = ROOT.TFile.Open(path_he3 + '/SystematicsAll.root')
        file_he3_syst_eff_prim = ROOT.TFile.Open(path_he3 + '/SystematicsEfficiencyPrimary.root')
        file_hyp_syst = ROOT.TFile.Open(path_hyp + '/Systematics.root')
        file_he3_scaling = ROOT.TFile.Open('He3_PbPb/MaterialRadiusCheck.root')
        file_proton_scaling = ROOT.TFile.Open('Proton_PbPb/MaterialRadiusCheck.root')
        file_pion_scaling = ROOT.TFile.Open('Pion_PbPb/MaterialRadiusCheck.root')

        file_out.cd()
        ratio_he3 = file_he3.Get(f'1.0_89_0.1_2.5_1_1_1/fRatio_{cent[0]}_{cent[1]}')
        ratio_hyp = file_hyp.Get(f'fRatio_{cent[0]}_{cent[1]}')
        ratio_proton = file_proton.Get(f'fRatio_{cent[0]}_{cent[1]}')
        ratio_pion = file_pion.Get(f'fRatio_{cent[0]}_{cent[1]}')
        ratio_he3_distribution = file_he3_syst.Get(f'hist/fFitPar_{cent[0]}_{cent[1]}')
        ratio_he3_distribution_eff_prim = file_he3_syst_eff_prim.Get(f'fRatioDistribution_{cent[0]}_{cent[1]}')
        ratio_hyp_distribution = file_hyp_syst.Get(f'fParameterDistribution_{cent[0]}_{cent[1]}')
        ratio_proton_pt_correlated = file_proton.Get(f'fRatioDistributionTrials_{cent[0]}_{cent[1]}')
        ratio_pion_pt_correlated = file_pion.Get(f'fRatioDistributionTrials_{cent[0]}_{cent[1]}')

        # apply material budget variation
        f_material_budget_he3 = file_material_budget.Get(f"fMHe3EffError_{cent[0]}_{cent[1]}")
        f_material_budget_antihe3 = file_material_budget.Get(f"fAHe3EffError_{cent[0]}_{cent[1]}")
        f_material_budget_hyp = file_material_budget_hyp.Get(f"{cent[0]}_{cent[1]}/BGBW/fCorrection_matter_{cent[0]}_{cent[1]}_BGBW")
        f_material_budget_antihyp = file_material_budget_hyp.Get(f"{cent[0]}_{cent[1]}/BGBW/fCorrection_antimatter_{cent[0]}_{cent[1]}_BGBW")
        f_material_budget_proton = file_material_budget.Get(f"fMProtonEffError_{cent[0]}_{cent[1]}")
        f_material_budget_antiproton = file_material_budget.Get(f"fAProtonEffError_{cent[0]}_{cent[1]}")
        f_material_budget_proton_lowPT = file_material_budget.Get(f"fMProtonLowPTEffError_{cent[0]}_{cent[1]}")
        f_material_budget_antiproton_lowPT = file_material_budget.Get(f"fAProtonLowPTEffError_{cent[0]}_{cent[1]}")
        f_material_budget_pion = file_material_budget.Get(f"fMPionEffError_{cent[0]}_{cent[1]}")
        f_material_budget_antipion = file_material_budget.Get(f"fAPionEffError_{cent[0]}_{cent[1]}")
        for i_bins in range(ratio_he3.GetNbinsX()):
            tmp_bin_pt = ratio_he3.GetBinCenter(i_bins)
            tmp_error_he3 = f_material_budget_he3.GetFunction("fitFunction").Eval(tmp_bin_pt)
            tmp_error_antihe3 = f_material_budget_antihe3.GetFunction("fitFunction").Eval(tmp_bin_pt)
            tmp_bin_content = ratio_he3.GetBinContent(i_bins)
            ratio_he3.SetBinContent(i_bins,tmp_bin_content*(1.+sign*np.sqrt(tmp_error_he3*tmp_error_he3+tmp_error_antihe3*tmp_error_antihe3)))
        ratio_he3.Fit("pol0","QR","",2.,10.)
        ratio_he3.Write()
        for i_bins in range(ratio_hyp.GetNbinsX()):
            tmp_bin_ct = ratio_hyp.GetBinCenter(i_bins)
            tmp_error_hyp = f_material_budget_hyp.GetBinContent(f_material_budget_hyp.FindBin(tmp_bin_ct))
            tmp_error_antihyp = f_material_budget_antihyp.GetBinContent(f_material_budget_hyp.FindBin(tmp_bin_ct))
            tmp_bin_content = ratio_hyp.GetBinContent(i_bins)
            ratio_hyp.SetBinContent(i_bins,tmp_bin_content*(1.+sign*np.sqrt(tmp_error_hyp*tmp_error_hyp+tmp_error_antihyp*tmp_error_antihyp)))
        ratio_hyp.Fit("pol0","QR","",0.,35.)
        ratio_hyp.Write()
        for i_bins in range(ratio_proton.GetNbinsX()):
            tmp_bin_pt = ratio_proton.GetBinCenter(i_bins)
            tmp_error_proton = f_material_budget_proton.GetFunction("fitFunction").Eval(tmp_bin_pt)
            tmp_error_antiproton = f_material_budget_antiproton.GetFunction("fitFunction").Eval(tmp_bin_pt)
            tmp_bin_content = ratio_proton.GetBinContent(i_bins)
            if tmp_bin_pt<0.81:
                tmp_error_proton = f_material_budget_proton_lowPT.GetFunction("fitFunction").Eval(tmp_bin_pt)
                tmp_error_antiproton = f_material_budget_antiproton_lowPT.GetFunction("fitFunction").Eval(tmp_bin_pt)
            ratio_proton.SetBinContent(i_bins,tmp_bin_content*(1.+sign*np.sqrt(tmp_error_proton*tmp_error_proton+tmp_error_antiproton*tmp_error_antiproton)))
        ratio_proton.Fit("pol0","QR","",0.5,3.)
        ratio_proton.Write()
        for i_bins in range(ratio_pion.GetNbinsX()):
            tmp_bin_pt = ratio_pion.GetBinCenter(i_bins)
            tmp_error_pion = f_material_budget_pion.GetFunction("fitFunction").Eval(tmp_bin_pt)
            tmp_error_antipion = f_material_budget_antipion.GetFunction("fitFunction").Eval(tmp_bin_pt)
            tmp_bin_content = ratio_pion.GetBinContent(i_bins)
            ratio_pion.SetBinContent(i_bins,tmp_bin_content*(1.+sign*np.sqrt(tmp_error_pion*tmp_error_pion+tmp_error_antipion*tmp_error_antipion)))
        ratio_pion.Fit("pol0","QR","",0.7,1.6)
        ratio_pion.Write()

        # get fit functions
        fit_he3 = ratio_he3.GetFunction("pol0")
        fit_hyp = ratio_hyp.GetFunction("pol0")
        fit_proton = ratio_proton.GetFunction("pol0")
        fit_pion = ratio_pion.GetFunction("pol0")

        # get ratios and errors
        ratio_he3 = fit_he3.GetParameter(0)
        ratio_he3_err = fit_he3.GetParError(0)/ratio_he3

        ratio_hyp = fit_hyp.GetParameter(0)
        ratio_hyp_err = fit_hyp.GetParError(0)/ratio_hyp
        
        ratio_proton = fit_proton.GetParameter(0)
        ratio_proton_err = 0#fit_proton.GetParError(0)

        ratio_pion = fit_pion.GetParameter(0)
        ratio_pion_err = 0#fit_proton.GetParError(0)

        # systematic error
        syst_he3 = ratio_he3_distribution.GetRMS()/ratio_he3
        syst_he3_eff_prim = ratio_he3_distribution_eff_prim.GetRMS()/ratio_he3
        syst_hyp = ratio_hyp_distribution.GetRMS()/ratio_hyp

        if apply_scaling_radius:
            a_scaling = file_he3_scaling.Get(f"GraphA")
            m_scaling = file_he3_scaling.Get(f"GraphM")
            a_scaling_val = a_scaling.GetFunction("pol0").GetParameter(0)
            m_scaling_val = m_scaling.GetFunction("pol0").GetParameter(0)
            ratio_he3 = ratio_he3*(1+m_scaling_val-a_scaling_val)
            ratio_hyp = ratio_hyp*(1+m_scaling_val-a_scaling_val)
            a_scaling = file_proton_scaling.Get(f"GraphA")
            m_scaling = file_proton_scaling.Get(f"GraphM")
            a_scaling_val = a_scaling.GetFunction("pol0").GetParameter(0)
            m_scaling_val = m_scaling.GetFunction("pol0").GetParameter(0)
            ratio_proton = ratio_proton*(1+m_scaling_val-a_scaling_val)
            a_scaling = file_pion_scaling.Get(f"GraphA")
            m_scaling = file_pion_scaling.Get(f"GraphM")
            a_scaling_val = a_scaling.GetFunction("pol0").GetParameter(0)
            m_scaling_val = m_scaling.GetFunction("pol0").GetParameter(0)
            ratio_pion = ratio_pion*(1+m_scaling_val-a_scaling_val)

        ratio_he3_err = ratio_he3_err*ratio_he3
        ratio_hyp_err = ratio_hyp_err*ratio_hyp

        syst_he3 = syst_he3*ratio_he3
        syst_he3_eff_prim = syst_he3_eff_prim*ratio_he3
        syst_he3_abs = np.sqrt(0.00861352*0.00861352+0.00473124*0.00473124)*ratio_he3 # from absorption cross section variation
        syst_he3 = np.sqrt(syst_he3*syst_he3+syst_he3_abs*syst_he3_abs+syst_he3_eff_prim*syst_he3_eff_prim)
        syst_hyp = syst_hyp*ratio_hyp
        syst_hyp_abs = np.sqrt(0.00861352*0.00861352+0.00473124*0.00473124)*ratio_hyp # from absorption cross section variation
        syst_hyp = np.sqrt(syst_hyp*syst_hyp+syst_hyp_abs*syst_hyp_abs)
        syst_proton = fit_proton.GetParError(0)
        syst_pion = fit_pion.GetParError(0)
        syst_proton_pt_correlated = ratio_proton_pt_correlated.GetRMS()
        syst_proton_abs = np.sqrt(0.00115085*0.00115085+0.000172469*0.000172469)*ratio_proton # from absorption cross section variation
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
        
        # format fit parameter and baryon chemical potential values and errors
        formatted_fit_parameter_0 = "{:.4f}".format(fit_parameter_0)
        formatted_fit_parameter_error_0 = "{:.4f}".format(fit_parameter_error_0)
        formatted_mu_b = "{:.2f}".format(fit_parameter_0*155)
        mu_b_error = np.sqrt(fit_parameter_error_0*fit_parameter_error_0/fit_parameter_0/fit_parameter_0)*fit_parameter_0*155
        formatted_mu_b_error = "{:.2f}".format(fit_parameter_error_0*155)

        formatted_mu_I = "{:.2f}".format(fit_parameter_1*155)
        mu_I_error = np.sqrt(fit_parameter_error_1*fit_parameter_error_1/fit_parameter_1/fit_parameter_1)*np.abs(fit_parameter_1)*155
        formatted_mu_I_error = "{:.2f}".format(np.abs(fit_parameter_error_1*155))
        
        # chi2 text
        text_chi2 = ROOT.TLatex(-0.27, -0.44, "#chi^{2}/NDF = "+formatted_chi2+"/"+str(fit_expo.GetNDF()))
        text_chi2.SetTextSize(TLATEX_TEXT_SIZE)
        text_chi2.SetTextColor(ROOT.kBlack)
        
        # # mu_b / T ratio
        # text_mu_b_over_T = ROOT.TLatex(0.5, 0.8, "#mu_{#it{B}}/#it{T} = "+formatted_fit_parameter+" #pm "+formatted_fit_parameter_error)
        # text_mu_b_over_T.SetTextSize(TLATEX_TEXT_SIZE)
        # text_mu_b_over_T.SetTextColor(ROOT.kBlack)

        # # T = 155 +/- 2 MeV
        # text_T = ROOT.TLatex(0.5, 0.7, "#it{T} = 155 #pm 2 MeV")
        # text_T.SetTextSize(TLATEX_TEXT_SIZE)
        # text_T.SetTextColor(ROOT.kBlack)

        # write to file
        ratios_vs_b.Write()

        if sign!=0:
            results_muB.append(fit_parameter_0*155)
            results_muI.append(fit_parameter_1*155)
            continue

        # save plots
        error_mb = 0.5*(results_muB[0]-results_muB[1])
        error_mi = 0.5*(results_muI[0]-results_muI[1])
        print(f"* * * * * * E R R O R * * * * * *")
        print(f"MB Error = {error_mb}")
        print(f"* * * * * * * * * * * * * * * * *")
        formatted_mb_error = "{:.2f}".format(error_mb)
        formatted_mi_error = "{:.2f}".format(error_mi)
        c = ROOT.TCanvas(f"c_{cent[0]}_{cent[1]}", "c")
        
        # mu_b at T = 155 MeV
        formatted_temperature_error = "{:.2f}".format(fit_parameter_0*2)
        text_mu_b = ROOT.TLatex(-0.27, -0.59, "#mu_{#it{B}} = "+formatted_mu_b+" #pm "+formatted_mu_b_error+" #pm "+formatted_mb_error+" #pm "+formatted_temperature_error+" MeV")
        text_mu_b.SetTextSize(TLATEX_TEXT_SIZE)
        text_mu_b.SetTextColor(ROOT.kBlack)

        # mu_I3 at T = 155 MeV
        formatted_temperature_error = "{:.2f}".format(np.abs(fit_parameter_1)*2)
        text_mu_I = ROOT.TLatex(-0.27, -0.59, "#mu_{#it{I}_{3}} = "+formatted_mu_I+" #pm "+formatted_mu_I_error+" #pm "+formatted_mi_error+" #pm "+formatted_temperature_error+" MeV")
        text_mu_I.SetTextSize(TLATEX_TEXT_SIZE)
        text_mu_I.SetTextColor(ROOT.kBlack)

        c.cd()
        ratios_vs_b.SetMarkerStyle(20)
        ratios_vs_b.SetMarkerSize(1.1)
        ratios_vs_b.SetLineColor(centrality_colors[i_cent])
        ratios_vs_b.SetLineWidth(2)
        ratios_vs_b.SetMarkerColor(centrality_colors[i_cent])
        ratios_vs_b.SetTitle(f"{cent[0]}-{cent[1]}%")
        fit_expo.SetMinimum(0.6)
        fit_expo.SetMaximum(1.2) # 1.23 to correctly visualise in file
        fit_expo.SetNpx(10)
        fit_expo.SetNpy(10)
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

        # make final plot for approval
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
        hRatiosParticle = ROOT.TH1D(f"hRatiosParticle_{cent[0]}_{cent[1]}",f" ",4,0,4)
        hNSigmaRatioFitParticle = ROOT.TH1D(f"hNSigmaRatioFitParticle_{cent[0]}_{cent[1]}",f"{cent[0]}-{cent[1]}%",4,0,4)
        hRatiosParticleFit = ROOT.TH1D(f"hRatiosParticleFit_{cent[0]}_{cent[1]}",f"{cent[0]}-{cent[1]}%",4,0,4)
        hRatiosParticle.GetYaxis().SetTitle("Ratio")
        hRatiosParticle.GetYaxis().CenterTitle()
        hRatiosParticle.GetXaxis().SetTickSize(0.07*0.4)
        hRatiosParticle.GetYaxis().SetTitleSize(0.13*0.45)
        hRatiosParticle.GetYaxis().SetLabelSize(0.12*0.42)
        hRatiosParticle.GetYaxis().SetTitleOffset(0.9)
        for i_part in range(0,4):
            hRatiosParticle.GetXaxis().SetBinLabel(i_part+1,particle_ratios[i_part])
            hRatiosParticle.GetXaxis().SetLabelSize(0.06)
            ratio = 1
            ratio_err = 0
            fit = 1
            if i_part == 0:
                ratio = ratio_pion
                ratio_err = np.sqrt(syst_pion*syst_pion+stat_pion*stat_pion)
                fit = fit_expo.Eval(0,1)
            if i_part == 1:
                ratio = ratio_proton
                ratio_err = np.sqrt(syst_proton*syst_proton+stat_proton*stat_proton)
                fit = fit_expo.Eval(3,0.5)
            if i_part == 2:
                ratio = ratio_hyp
                ratio_err = np.sqrt(syst_hyp*syst_hyp+stat_hyp*stat_hyp)
                fit = fit_expo.Eval(8,0)
            if i_part == 3:
                ratio = ratio_he3
                ratio_err = np.sqrt(syst_he3*syst_he3+stat_he3*stat_he3)
                fit = fit_expo.Eval(9,0.5)
            hRatiosParticle.SetBinContent(i_part+1,ratio)
            hRatiosParticle.SetBinError(i_part+1,ratio_err)
            hRatiosParticleFit.SetBinContent(i_part+1,fit)
            hRatiosParticleFit.SetBinError(i_part+1,0)
            print(f"fit = {fit}")
        hRatiosParticle.GetYaxis().SetRangeUser(0.67,1.13)
        gRatiosParticle = ROOT.TGraphErrors(hRatiosParticle)
        gRatiosParticle.SetName(f"gRatiosParticle_{cent[0]}_{cent[1]}")
        gRatiosParticleFit = ROOT.TGraphErrors(hRatiosParticleFit)
        gRatiosParticleFit.SetName(f"gRatiosParticleFit_{cent[0]}_{cent[1]}")
        for i_part in range(0,4):
            gRatiosParticle.SetPointError(i_part,0,hRatiosParticle.GetBinError(i_part+1))
            gRatiosParticleFit.SetPointError(i_part,0.3,0)
        hijing_predictions = file_hijing.Get(f"fRatios_{cent[0]}_{cent[1]}")
        gHijingPredictions = ROOT.TGraphErrors()
        gHijingPredictions.AddPoint(gRatiosParticle.GetPointX(0),hijing_predictions.GetBinContent(1,21)) # pion
        gHijingPredictions.SetPointError(0,0.25,hijing_predictions.GetBinError(1,21)) # pion
        gHijingPredictions.AddPoint(gRatiosParticle.GetPointX(1),hijing_predictions.GetBinContent(4,16)) # proton
        gHijingPredictions.SetPointError(1,0.25,hijing_predictions.GetBinError(4,16)) # pion

        # data to fit ratio
        for i_ticks in range(5):
            hNSigmaRatioFitParticle.GetYaxis().SetNdivisions(5)
        for i_part in range(0,4):
            hNSigmaRatioFitParticle.GetXaxis().SetBinLabel(i_part+1,particle_ratios[i_part])
            hNSigmaRatioFitParticle.GetXaxis().SetLabelSize(0.2)
            ratio = 1
            ratio_err = 0
            fit = 1
            if i_part == 0:
                ratio = ratio_pion
                ratio_err = np.sqrt(syst_pion*syst_pion+stat_pion*stat_pion)
                fit = fit_expo.Eval(0,1)
            if i_part == 1:
                ratio = ratio_proton
                ratio_err = np.sqrt(syst_proton*syst_proton+stat_proton*stat_proton)
                fit = fit_expo.Eval(3,0.5)
            if i_part == 2:
                ratio = ratio_hyp
                ratio_err = np.sqrt(syst_hyp*syst_hyp+stat_hyp*stat_hyp)
                fit = fit_expo.Eval(8,0)
            if i_part == 3:
                ratio = ratio_he3
                ratio_err = np.sqrt(syst_he3*syst_he3+stat_he3*stat_he3)
                fit = fit_expo.Eval(9,0.5)
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
        hNSigmaRatioFitParticle.GetYaxis().SetTitle("#frac{data - fit}{#sigma_{data}}")
        hNSigmaRatioFitParticle.GetYaxis().CenterTitle()
        hNSigmaRatioFitParticle.GetYaxis().SetRangeUser(-2.3,2.3)
        gHijingNsigma = ROOT.TGraphErrors(gHijingPredictions)
        gHijingNsigma.SetPointY(0,(gHijingNsigma.GetPointY(0)-fit_expo.Eval(0,1))/np.sqrt(syst_pion*syst_pion+stat_pion*stat_pion))
        gHijingNsigma.SetPointY(1,(gHijingNsigma.GetPointY(1)-fit_expo.Eval(3,0.5))/np.sqrt(syst_proton*syst_proton+stat_proton*stat_proton))
        hLineZero = ROOT.TLine(0.,0.,4.,0.)
        hLineZero.SetLineColor(ROOT.kBlack)
        hLineZero.SetLineStyle(ROOT.kDashed)
        hLineOnePos = ROOT.TLine(0.,1.,4.,1.)
        hLineOnePos.SetLineColor(ROOT.kBlack)
        hLineOnePos.SetLineStyle(ROOT.kDashed)
        hLineOneNeg = ROOT.TLine(0.,-1.,4.,-1.)
        hLineOneNeg.SetLineColor(ROOT.kBlack)
        hLineOneNeg.SetLineStyle(ROOT.kDashed)
        hLineTwoPos = ROOT.TLine(0.,2.,4.,2.)
        hLineTwoPos.SetLineColor(ROOT.kBlack)
        hLineTwoPos.SetLineStyle(ROOT.kDashed)
        hLineTwoNeg = ROOT.TLine(0.,-2.,4.,-2.)
        hLineTwoNeg.SetLineColor(ROOT.kBlack)
        hLineTwoNeg.SetLineStyle(ROOT.kDashed)

        hRatiosParticle.GetYaxis().SetRangeUser(0.72,1.08)
        text_alice = ROOT.TLatex(0.18,0.3,"ALICE Preliminary")
        text_alice.SetNDC(True)
        text_alice.SetTextFont(63)
        text_alice.SetTextSize(25)
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
        for i_part in range(0,4):
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
        gHijingPredictions.SetLineColor(ROOT.kGreen+1)
        gHijingPredictions.SetLineWidth(2)
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
        gHijingNsigma.SetLineColor(ROOT.kGreen+1)
        gHijingNsigma.SetLineWidth(2)
        hNSigmaRatioFitParticle.Draw("psame")
        # gHijingNsigma.Draw("pe same")
        cRatiosParticle.Write()
        cRatiosParticle.Print(f"SHMfitToRatios_{cent[0]}_{cent[1]}.eps")
        gRatiosParticle.Write()
        gRatiosParticleFit.Write()
        hNSigmaRatioFitParticle.Write()

        gMuBCent.AddPoint(n_part[i_cent],fit_parameter_0*155)
        gMuBCent.SetPointError(gMuBCent.GetN()-1,n_part_err[i_cent],mu_b_error)
        gMuBCentCorrUnc.AddPoint(n_part[i_cent],fit_parameter_0*155)
        gMuBCentCorrUnc.SetPointError(gMuBCent.GetN()-1,10,error_mb)
        gMuICent.AddPoint(0.5*(cent[1]+cent[0]),fit_parameter_1*155)
        gMuICent.SetPointError(gMuBCent.GetN()-1,0.,mu_I_error)
        print(f'mu_I_error = {mu_I_error}')

gMuBCent.SetName("MuBCent")
gMuICent.SetName("MuICent")

legMuBCent = ROOT.TLegend(0.18,0.15,0.60,0.28)
gPublishedResult = ROOT.TGraphErrors()
gPublishedResult.AddPoint(354.7,0.7)
gPublishedResult.SetPointError(0,33,3.8)
gPublishedResult.SetMarkerStyle(0)
gPublishedResult.SetMarkerSize(0)
gPublishedResult.SetLineWidth(0)
gPublishedResult.SetFillStyle(3345)
gPublishedResult.SetFillColor(ROOT.kBlack)
text_alice = ROOT.TLatex(0.18,0.8,"ALICE Preliminary")
text_alice.SetTextFont(63)
text_alice.SetTextSize(25)
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

legMuICent = ROOT.TLegend(0.6,0.6,0.8,0.8)
cMuICent = ROOT.TCanvas("cMuICent","cMuICent")
gHijingMuI = file_hijing.Get("gMuICent")
cMuICent.cd()
gMuICent.SetMarkerStyle(20)
gMuICent.SetMarkerSize(1.0)
gMuICent.SetMarkerColor(ROOT.kRed)
gMuICent.SetLineColor(ROOT.kRed)
hMuICentFrame.GetYaxis().SetRangeUser(-0.3,0.5)
hMuICentFrame.Draw("pe")
gMuICent.Draw("pesame")
#gHijingMuI.Draw("e3same")
legMuICent.AddEntry(gMuICent,"Data")
legMuICent.AddEntry(gHijingMuI,"HIJING")
#legMuICent.Draw("same")
cMuICent.Print("MuIvsCent.pdf")
gMuBCent.Write()
gMuICent.Write()

file_out.Close()
