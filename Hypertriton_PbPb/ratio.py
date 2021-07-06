#!/usr/bin/env python3
import os
import pickle
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import norm
import ROOT
import uproot
import yaml

SPEED_OF_LIGHT = 2.99792458
SPLIT = True

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

if not os.path.isdir("plots/cpt_and_lifetime"):
    os.mkdir("plots/cpt_and_lifetime")

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)

eff_cut_dict = pickle.load(open("file_eff_cut_dict", "rb"))
presel_eff_file = uproot.open('PreselEff.root')
signal_extraction_file = ROOT.TFile.Open('SignalExtraction.root')
signal_extraction_keys = uproot.open('SignalExtraction.root').keys()

abs_correction_file = ROOT.TFile.Open('He3_abs.root')
ratio_file = ROOT.TFile.Open('Ratio.root', 'recreate')

for i_cent_bins in range(len(CENTRALITY_LIST)):
    cent_bins = CENTRALITY_LIST[i_cent_bins]

    h_corrected_yields = [ROOT.TH1D(), ROOT.TH1D()]
    for i_split, split in enumerate(SPLIT_LIST):
        print(f'{i_split} -> {split}')
        # get preselection efficiency and abs correction histograms
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
        if cent_bins[0]==30:
            ct_bins_tmp = [0, 2, 4, 7, 14, 35]
        #if cent_bins[1] == 90:
        #    ct_bin_tmp = CT_BINS_CENT[i_cent_bins]
        bins = np.array(ct_bins_tmp, dtype=float)
        # print(bins)
        h_corrected_yields[i_split] = ROOT.TH1D(
            f'fYields_{split}_{cent_bins[0]}_{cent_bins[1]}', f'{split}, {cent_bins[0]}-{cent_bins[1]}%', len(bins)-1, bins)

        for ct_bins in zip(CT_BINS_CENT[i_cent_bins][:-1], CT_BINS_CENT[i_cent_bins][1:]):

            bin = f'{split}_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}'
            formatted_eff_cut = "{:.2f}".format(eff_cut_dict[bin])

            # look for plot with eff = eff_cut (or the nearest one)
            bkg_shape = 'pol1'
            eff_cut_increment = 0
            eff_cut_sign = -1
            while signal_extraction_keys.count(f'{bin}_{bkg_shape}/fInvMass_{formatted_eff_cut};1') == 0:
                if eff_cut_sign == -1:
                    eff_cut_increment += 0.01
                eff_cut_sign *= -1
                formatted_eff_cut = "{:.2f}".format(eff_cut_dict[bin]+eff_cut_increment*eff_cut_sign)

            # get signal
            h_raw_yield = signal_extraction_file.Get(f'{bin}_{bkg_shape}/fRawYields;1')
            eff_index = h_raw_yield.FindBin(float(formatted_eff_cut))
            raw_yield = h_raw_yield.GetBinContent(eff_index)
            raw_yield_error = h_raw_yield.GetBinError(eff_index)
            print(f'eff_cut = {formatted_eff_cut}, raw_yield = {raw_yield}+{raw_yield_error}')

            # apply corrections
            # 1. efficiency (presel x BDT)
            presel_eff_map = np.logical_and(
                presel_eff_bin_centers > ct_bins[0],
                presel_eff_bin_centers < ct_bins[1])
            presel_eff = presel_eff_counts[presel_eff_map]
            bdt_eff = float(formatted_eff_cut)
            print(f'bin: {presel_eff_map}, presel_eff: {presel_eff}')
            eff = presel_eff * eff_cut_dict[bin]
            # 2. absorption correction
            abs = g_abs_correction.GetPointY(CT_BINS_CENT[i_cent_bins].index(ct_bins[0]))
            print(f"absorption correction for point {CT_BINS_CENT[i_cent_bins].index(ct_bins[0])}: {abs}")

            ct_bin_index = h_corrected_yields[i_split].FindBin(ct_bins[0]+0.5)
            h_corrected_yields[i_split].SetBinContent(ct_bin_index, raw_yield/eff[0]/abs)
            h_corrected_yields[i_split].SetBinError(ct_bin_index, raw_yield_error/eff[0]/abs)

        # set labels
        h_corrected_yields[i_split].GetXaxis().SetTitle("#it{c}t (cm)")
        h_corrected_yields[i_split].GetYaxis().SetTitle("d#it{N}/d(#it{c}t) (cm^{-1})")
        #h_corrected_yields[i_split].Scale(1, "width")
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

        # fit with exponential pdf
        fit_function_expo = ROOT.TF1("expo", "expo", 2, 35)
        if cent_bins[0] == 30:
            fit_function_expo = ROOT.TF1("expo", "expo", 2, 14)
        elif cent_bins[1] == 90:
            fit_function_expo = ROOT.TF1("expo", "expo", 0, 35)
        h_corrected_yields[i_split].Fit(fit_function_expo, "RMIS+")

        # compute and display lifetime
        tau = -1/fit_function_expo.GetParameter(1)*100/SPEED_OF_LIGHT # ps
        tau_error = fit_function_expo.GetParError(1)*100/SPEED_OF_LIGHT/fit_function_expo.GetParameter(1)/fit_function_expo.GetParameter(1) # ps
        tau_text = ROOT.TLatex(20, 100, '#tau = ' + "{:.0f}".format(tau) + '#pm' + "{:.0f}".format(tau_error) + ' ps')
        tau_text.SetTextSize(0.035)
        
        # display chi2
        formatted_chi2_lifetime = "{:.2f}".format(fit_function_expo.GetChisquare())
        chi2_lifetime_text = ROOT.TLatex(20, 70, "#chi^{2}/NDF = "+formatted_chi2_lifetime+"/"+str(fit_function_expo.GetNDF()))
        chi2_lifetime_text.SetTextSize(0.035)

        # draw on canvas
        canv = ROOT.TCanvas()
        ROOT.gStyle.SetOptStat(0)
        canv.SetTicks(1, 1)
        canv.SetName(h_corrected_yields[i_split].GetName())
        canv.SetTitle(h_corrected_yields[i_split].GetTitle())
        canv.cd()
        h_corrected_yields[i_split].Draw("")
        tau_text.Draw("same")
        chi2_lifetime_text.Draw("same")
        canv.SetLogy()
        canv.Write() # write to file
        canv.Print(f'plots/cpt_and_lifetime/{h_corrected_yields[i_split].GetName()}.png')

    # ratios
    ROOT.gStyle.SetOptStat(0)
    h_ratio = ROOT.TH1D(h_corrected_yields[0])
    h_ratio.SetName(f'fRatio_{cent_bins[0]}_{cent_bins[1]}')
    h_ratio.SetTitle(f'{cent_bins[0]}-{cent_bins[1]}%')
    h_ratio.Divide(h_corrected_yields[0], h_corrected_yields[1], 1, 1)
    h_ratio.GetYaxis().SetTitle("ratio ^{3}_{#bar{#Lambda}}#bar{H} / ^{3}_{#Lambda}H")
    h_ratio.GetYaxis().SetRangeUser(0., 1.8)
    h_ratio.GetXaxis().SetRangeUser(0., 35.)
    h_ratio.SetMarkerStyle(20)
    h_ratio.SetMarkerSize(0.8)
    h_ratio.Fit("pol0")
    h_ratio.Write()

    # plot ratios
    c = ROOT.TCanvas("c", "c")
    c.SetTicks(1, 1)
    c.cd()
    h_ratio.Draw()
    formatted_ratio = "{:.2f}".format(h_ratio.GetFunction("pol0").GetParameter(0))
    formatted_ratio_error = "{:.2f}".format(h_ratio.GetFunction("pol0").GetParError(0))
    text_x_position = 20
    ratio_text = ROOT.TLatex(text_x_position, 1.6, f"R = {formatted_ratio} #pm {formatted_ratio_error}")
    ratio_text.SetTextFont(44)
    ratio_text.SetTextSize(28)
    ratio_text.Draw("same")
    formatted_chi2 = "{:.2f}".format(h_ratio.GetFunction("pol0").GetChisquare())
    chi2_text = ROOT.TLatex(text_x_position, 1.45, "#chi^{2}/NDF = "+formatted_chi2+"/"+str(h_ratio.GetFunction("pol0").GetNDF()))
    chi2_text.SetTextFont(44)
    chi2_text.SetTextSize(28)
    chi2_text.Draw("same")
    c.Print(f"plots/{h_ratio.GetName()}.png")

    del h_corrected_yields
    del h_ratio

ratio_file.Close()
