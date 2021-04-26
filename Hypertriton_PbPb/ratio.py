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

eff_cut_dict = pickle.load(open("file_eff_cut_dict", "rb"))
presel_eff_file = uproot.open('PreselEff.root')
signal_extraction_file = ROOT.TFile.Open('SignalExtraction.root')
signal_extraction_keys = uproot.open('SignalExtraction.root').keys()

ratio_file = ROOT.TFile.Open('Ratio.root', 'recreate')

for i_cent_bins in range(len(CENTRALITY_LIST)):
    cent_bins = CENTRALITY_LIST[i_cent_bins]

    h_corrected_yields = [ROOT.TH1D(), ROOT.TH1D()]
    for i_split, split in enumerate(SPLIT_LIST):
        print(f'{i_split} -> {split}')
        # get preselection efficiency histogram
        presel_eff_counts, presel_eff_edges = presel_eff_file[
            f'fPreselEff_vs_ct_{split}_{cent_bins[0]}_{cent_bins[1]};1'].to_numpy()
        presel_eff_bin_centers = (presel_eff_edges[1:]+presel_eff_edges[:-1])/2

        # list of corrected yields
        ct_bins_tmp = [0]
        ct_bins_tmp += CT_BINS_CENT[i_cent_bins]
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
            presel_eff_map = np.logical_and(
                presel_eff_bin_centers > ct_bins[0],
                presel_eff_bin_centers < ct_bins[1])
            presel_eff = presel_eff_counts[presel_eff_map]
            bdt_eff = float(formatted_eff_cut)
            print(f'bin: {presel_eff_map}, presel_eff: {presel_eff}')
            eff = presel_eff * eff_cut_dict[bin]
            # print(f'corrected yield = {raw_yield/eff[0]}')
            ct_bin_index = h_corrected_yields[i_split].FindBin(ct_bins[0]+0.5)
            # print(f'index = {ct_bin_index}')
            h_corrected_yields[i_split].SetBinContent(ct_bin_index, raw_yield/eff[0])
            h_corrected_yields[i_split].SetBinError(ct_bin_index, raw_yield_error/eff[0])

        # set labels
        h_corrected_yields[i_split].GetXaxis().SetTitle("#it{c}t (cm)")
        h_corrected_yields[i_split].GetYaxis().SetTitle("d#it{N}/d(#it{c}t) (cm^{-1})")
        h_corrected_yields[i_split].Scale(1., "width")
        h_corrected_yields[i_split].Write()

        # fit with exponential pdf
        fit_function_expo = ROOT.TF1("expo", "expo", 2, 35)
        if cent_bins[0] == 30:
            fit_function_expo = ROOT.TF1("expo", "expo", 2, 14)
        elif cent_bins[1] == 90:
            fit_function_expo = ROOT.TF1("expo", "expo", 2, 35)
        h_corrected_yields[i_split].Fit(fit_function_expo, "RMLS+")

        # compute lifetime
        tau = -1/fit_function_expo.GetParameter(1)*100/SPEED_OF_LIGHT # ps
        tau_error = fit_function_expo.GetParError(1)*100/SPEED_OF_LIGHT/fit_function_expo.GetParameter(1)/fit_function_expo.GetParameter(1) # ps
        tau_text = ROOT.TLatex(4, 0.9*h_corrected_yields[i_split].GetMaximum(), '#tau = ' + "{:.2f}".format(tau) + '#pm' + "{:.2f}".format(tau_error) + ' ps')
        tau_text.SetTextSize(0.035)

        # draw on canvas
        canv = ROOT.TCanvas()
        canv.SetName(h_corrected_yields[i_split].GetName())
        canv.SetTitle(h_corrected_yields[i_split].GetTitle())
        canv.cd()
        h_corrected_yields[i_split].Draw("")
        tau_text.Draw("same")
        canv.SetLogy()
        canv.Write() # write to file

    # ratios
    h_ratio = ROOT.TH1D(h_corrected_yields[0])
    h_ratio.SetName(f'fRatio_{cent_bins[0]}_{cent_bins[1]}')
    h_ratio.SetTitle(f'{cent_bins[0]}-{cent_bins[1]}%')
    h_ratio.Divide(h_corrected_yields[0], h_corrected_yields[1], 1, 1)
    h_ratio.GetYaxis().SetTitle("ratio ^{3}_{#bar{#Lambda}}#bar{H} / ^{3}_{#Lambda}H")
    h_ratio.SetMarkerStyle(20)
    h_ratio.SetMarkerSize(0.8)
    h_ratio.Fit("pol0")
    h_ratio.Write()

    del h_corrected_yields
    del h_ratio

ratio_file.Close()
