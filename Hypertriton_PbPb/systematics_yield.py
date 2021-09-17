#!/usr/bin/env python3
import os
import pickle
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import ROOT
import uproot
import yaml

SPEED_OF_LIGHT = 2.99792458
N_TRIALS = 10000
SPLIT = True
MAX_EFF = 0.99
THRESH_EFF = 0.90

# avoid pandas warning
warnings.simplefilter(action='ignore', category=FutureWarning)
ROOT.gROOT.SetBatch()

##################################################################
# read configuration file
##################################################################
config = 'config.yaml'
with open(os.path.expandvars(config), 'r') as stream:
    try:
        params = yaml.full_load(stream)
    except yaml.YAMLError as exc:
        print(exc)

ANALYSIS_RESULTS_PATH = params['ANALYSIS_RESULTS_PATH']
CT_BINS_CENT = params['CT_BINS_CENT']
CENTRALITY_LIST = params['CENTRALITY_LIST']
RANDOM_STATE = params['RANDOM_STATE']
##################################################################

# split matter/antimatter
SPLIT_LIST = ['']
if SPLIT:
    SPLIT_LIST = ['antimatter', 'matter']

eff_cut_dict = pickle.load(open("file_eff_cut_dict", "rb"))
presel_eff_file = uproot.open('PreselEff.root')
signal_extraction_file = ROOT.TFile.Open('SignalExtraction.root')
signal_extraction_keys = uproot.open('SignalExtraction.root').keys()

abs_correction_file = ROOT.TFile.Open('He3_abs.root')
systematics_file = ROOT.TFile.Open('SystematicsYield.root', 'recreate')

for i_cent_bins in range(len(CENTRALITY_LIST)):
    cent_bins = CENTRALITY_LIST[i_cent_bins]

    h_parameter_distribution = [ROOT.TH1D(f'fParameterDistribution_A_{cent_bins[0]}_{cent_bins[1]}', f'{cent_bins[0]}-{cent_bins[1]}%', 200, 0, 2), ROOT.TH1D(f'fParameterDistribution_M_{cent_bins[0]}_{cent_bins[1]}', f'{cent_bins[0]}-{cent_bins[1]}%', 200, 0, 2)]
    h_prob_distribution = [ROOT.TH1D(f'fProbDistribution_A_{cent_bins[0]}_{cent_bins[1]}', f'{cent_bins[0]}-{cent_bins[1]}%', 100, 0, 1),ROOT.TH1D(f'fProbDistribution_M_{cent_bins[0]}_{cent_bins[1]}', f'{cent_bins[0]}-{cent_bins[1]}%', 100, 0, 1)]

    systematics_file.mkdir(f'{cent_bins[0]}_{cent_bins[1]}')
    ####################################################
    # MULTITRIAL ESTIMATION OF SYSTEMATIC UNCERTAINTY
    ####################################################
    # sample 10000 combinations
    i_trial=0
    while i_trial < N_TRIALS:
        h_corrected_yields = [ROOT.TH1D(), ROOT.TH1D()]
        if (i_trial % 100) == 0:
            print(f'{cent_bins[0]}-{cent_bins[1]}%: i_trial = {i_trial}')
        for i_split, split in enumerate(SPLIT_LIST):

            # get preselection efficiency histogram
            presel_eff_counts, presel_eff_edges = presel_eff_file[
                f'fPreselEff_vs_ct_{split}_{cent_bins[0]}_{cent_bins[1]};1'].to_numpy()
            presel_eff_bin_centers = (presel_eff_edges[1:]+presel_eff_edges[:-1])/2

            func_name = 'BGBW'
            if cent_bins[1] == 90:
                func_name = 'BlastWave'
            g_abs_correction = ROOT.TGraphAsymmErrors()
            g_abs_correction = abs_correction_file.Get(f"{cent_bins[0]}_{cent_bins[1]}/{func_name}/fEffCt_{split}_{cent_bins[0]}_{cent_bins[1]}_{func_name}")

            # list of corrected yields
            ct_bins_tmp = [0]
            ct_bins_tmp += CT_BINS_CENT[i_cent_bins]
            bins = np.array(ct_bins_tmp, dtype=float)

            h_corrected_yields[i_split] = ROOT.TH1D(
                f'fYields_{split}_{cent_bins[0]}_{cent_bins[1]}', f'{split}, {cent_bins[0]}-{cent_bins[1]}%', len(bins)-1, bins)

            for ct_bins in zip(CT_BINS_CENT[i_cent_bins][:-1], CT_BINS_CENT[i_cent_bins][1:]):

                bin = f'{split}_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}'
                formatted_eff_cut = "{:.2f}".format(eff_cut_dict[bin])

                # look for plot with eff = eff_cut (or the nearest one)
                bkg_shape = 'pol1'
                if ROOT.gRandom.Rndm() > 0.5:
                    bkg_shape = 'expo'

                eff_cut_increment = 0
                eff_cut_sign = -1
                while signal_extraction_keys.count(f'{bin}_{bkg_shape}/fInvMass_{formatted_eff_cut};1') == 0:
                    if eff_cut_sign == -1:
                        eff_cut_increment += 0.01
                    eff_cut_sign *= -1
                    formatted_eff_cut = "{:.2f}".format(eff_cut_dict[bin]+eff_cut_increment*eff_cut_sign)

                # random sample of cut (upper and lower limits from significance scan)
                bin_range = f'{split}_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}_range'
                eff_cut_range = eff_cut_dict[bin_range]/100 - 0.01
                #print(f"BDT efficiency cut variation range: {eff_cut_range}")
                lower_limit = float(formatted_eff_cut) - eff_cut_range
                upper_limit = float(formatted_eff_cut) + eff_cut_range
                if float(formatted_eff_cut) > THRESH_EFF:
                    upper_limit = MAX_EFF
                random_cut = lower_limit + ROOT.gRandom.Rndm()*(upper_limit-lower_limit)
                formatted_eff_cut = "{:.2f}".format(random_cut)

                # get signal
                h_raw_yield = signal_extraction_file.Get(f'{bin}_{bkg_shape}/fRawYields;1')
                eff_index = h_raw_yield.FindBin(float(formatted_eff_cut))
                raw_yield = h_raw_yield.GetBinContent(eff_index)
                raw_yield_error = h_raw_yield.GetBinError(eff_index)
                if raw_yield < 0.1:
                    continue

                # apply corrections
                presel_eff_map = np.logical_and(
                    presel_eff_bin_centers > ct_bins[0],
                    presel_eff_bin_centers < ct_bins[1])
                presel_eff = presel_eff_counts[presel_eff_map]
                bdt_eff = float(formatted_eff_cut)
                eff = presel_eff * eff_cut_dict[bin]
                abs = g_abs_correction.GetPointY(CT_BINS_CENT[i_cent_bins].index(ct_bins[0]))
                ct_bin_index = h_corrected_yields[i_split].FindBin(ct_bins[0]+0.5)

                h_corrected_yields[i_split].SetBinContent(ct_bin_index, raw_yield/eff[0]/abs)
                h_corrected_yields[i_split].SetBinError(ct_bin_index, raw_yield_error/eff[0]/abs)

            # set labels
            h_corrected_yields[i_split].GetXaxis().SetTitle("#it{c}t (cm)")

            #h_corrected_yields[i_split].Scale(1., "width")
            for i_bin in range(len(bins))[2:]:
                bin_width = h_corrected_yields[i_split].GetBinWidth(i_bin)
                print(f"bin: {h_corrected_yields[i_split].GetBinLowEdge(i_bin)}; bin width: {bin_width}")
                bin_content = h_corrected_yields[i_split].GetBinContent(i_bin)
                bin_error = h_corrected_yields[i_split].GetBinError(i_bin)
                h_corrected_yields[i_split].SetBinContent(i_bin, bin_content/bin_width)
                h_corrected_yields[i_split].SetBinError(i_bin, bin_error/bin_width)
            h_corrected_yields[i_split].GetYaxis().SetRangeUser(1., 450.)
            h_corrected_yields[i_split].SetMarkerStyle(20)
            h_corrected_yields[i_split].SetMarkerSize(0.8)
            h_corrected_yields[i_split].Write()

            fit_function_expo = ROOT.TF1("expo", "expo", 2, 35)
            if cent_bins[0] == 30:
                fit_function_expo = ROOT.TF1("expo", "expo", 2, 14)
            elif cent_bins[1] == 90:
                fit_function_expo = ROOT.TF1("expo", "expo", 0, 35)
            res = h_corrected_yields[i_split].Fit(fit_function_expo, "RMIS+")

            # compute and display lifetime
            tau = -1/fit_function_expo.GetParameter(1)*100/SPEED_OF_LIGHT # ps
            tau_error = fit_function_expo.GetParError(1)*100/SPEED_OF_LIGHT/fit_function_expo.GetParameter(1)/fit_function_expo.GetParameter(1) # ps
            tau_text = ROOT.TLatex(20, 100, '#tau = ' + "{:.0f}".format(tau) + '#pm' + "{:.0f}".format(tau_error) + ' ps')
            tau_text.SetTextSize(0.05)
            
            # display chi2
            formatted_chi2_lifetime = "{:.2f}".format(fit_function_expo.GetChisquare())
            chi2_lifetime_text = ROOT.TLatex(20, 70, "#chi^{2}/NDF = "+formatted_chi2_lifetime+"/"+str(fit_function_expo.GetNDF()))
            chi2_lifetime_text.SetTextSize(0.05)

            # compute yield
            integral = float()
            integral_error = float()
            integral = fit_function_expo.IntegralError(0, 1.e9)
            integral_error = fit_function_expo.IntegralError(0, 1.e9)
            formatted_integral = "{:.7f}".format(integral/2./evts*1.e9)
            formatted_integral_error = "{:.7f}".format(integral_error/2./evts*1.e9)
            integral_yield = ROOT.TLatex(20, 40, "1/#it{N_{ev}} d#it{N}/d#it{y} = "+formatted_integral+" #pm "+formatted_integral_error)
            integral_yield.SetTextSize(0.05)

            systematics_file.cd(f'{cent_bins[0]}_{cent_bins[1]}')

            if fit_function_expo.GetProb() > 0.025 and fit_function_expo.GetProb() < 0.975 and res.Status() == 0 and res.Ndf() > 1:
                # h_ratio.Write()
                h_parameter_distribution[i_split].Fill(integral)
                h_prob_distribution[i_split].Fill(res.Prob())
                i_trial += 0.5

        del h_corrected_yields

    systematics_file.cd()
    for i_split in range(2):
        h_parameter_distribution[i_split].GetXaxis().SetRangeUser(0., 2.)
        h_parameter_distribution[i_split].GetYaxis().SetTitle("Entries")
        h_parameter_distribution[i_split].GetXaxis().SetTitle("yield")
        h_parameter_distribution[i_split].SetDrawOption("histo")
        h_parameter_distribution[i_split].SetLineWidth(2)
        h_parameter_distribution[i_split].SetFillStyle(3345)
        h_parameter_distribution[i_split].SetFillColor(ROOT.kBlue)
        h_parameter_distribution[i_split].Write()
        h_prob_distribution[i_split].GetXaxis().SetTitle("Prob")
        h_prob_distribution[i_split].GetYaxis().SetTitle("Entries")
        h_prob_distribution[i_split].Write()

        # plot systematics distribution
        c = ROOT.TCanvas("c", "c")
        ROOT.gStyle.SetOptStat(110001110)
        ROOT.gStyle.SetStatX(0.85);
        ROOT.gStyle.SetStatY(0.85);
        c.cd()
        c.SetTicks(1, 1)
        mean = h_parameter_distribution[i_split].GetMean()
        std_dev = h_parameter_distribution[i_split].GetRMS()
        h_parameter_distribution[i_split].Rebin(2)
        h_parameter_distribution[i_split].GetXaxis().SetRangeUser(mean-5*std_dev, mean+5*std_dev)
        h_parameter_distribution[i_split].Draw("histo")
        c.Print(f"plots/{h_parameter_distribution[i_split].GetName()}.png")

        del h_parameter_distribution[i_split]
        del h_prob_distribution[i_split]

systematics_file.Close()
