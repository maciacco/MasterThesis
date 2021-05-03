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

systematics_file = ROOT.TFile.Open('SystematicsLifetime.root', 'recreate')

for i_cent_bins in range(len(CENTRALITY_LIST)):
    cent_bins = CENTRALITY_LIST[i_cent_bins]
    if cent_bins[1] < 80:
        continue

    h_asymmetry_distribution = ROOT.TH1D(f'fAsymmetryDistribution_{cent_bins[0]}_{cent_bins[1]}', f'{cent_bins[0]}-{cent_bins[1]}%', 400, -200, 200)
    h_prob_distribution = ROOT.TH1D(f'fProbDistribution_{cent_bins[0]}_{cent_bins[1]}', f'{cent_bins[0]}-{cent_bins[1]}%', 100, 0, 1)
    h_lifetime = [ROOT.TH1D(), ROOT.TH1D()]
    h_lifetime[0] = ROOT.TH1D(f'fLifetimeAntimatterDistribution_{cent_bins[0]}_{cent_bins[1]}', f'{cent_bins[0]}-{cent_bins[1]}%', 1000, 0, 1000)
    h_lifetime[1] = ROOT.TH1D(f'fLifetimeMatterDistribution_{cent_bins[0]}_{cent_bins[1]}', f'{cent_bins[0]}-{cent_bins[1]}%', 1000, 0, 1000)
    h_fit_status = ROOT.TH1D(f'fFitStatus_{cent_bins[0]}_{cent_bins[1]}', f'{cent_bins[0]}-{cent_bins[1]}%', 100, 0, 99)

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

        lifetime_tmp = [-9999., -9999.]
        for i_split, split in enumerate(SPLIT_LIST):

            # get preselection efficiency histogram
            presel_eff_counts, presel_eff_edges = presel_eff_file[
                f'fPreselEff_vs_ct_{split}_{cent_bins[0]}_{cent_bins[1]};1'].to_numpy()
            presel_eff_bin_centers = (presel_eff_edges[1:]+presel_eff_edges[:-1])/2

            # list of corrected yields
            ct_bins_tmp = [0]
            ct_bins_tmp += CT_BINS_CENT[i_cent_bins]
            #if cent_bins[1] == 90:
            #    ct_bin_tmp = CT_BINS_CENT[i_cent_bins]
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

                # random sample of cut
                lower_limit = float(formatted_eff_cut) - 0.10
                upper_limit = float(formatted_eff_cut) + 0.10
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
                ct_bin_index = h_corrected_yields[i_split].FindBin(ct_bins[0]+0.5)

                h_corrected_yields[i_split].SetBinContent(ct_bin_index, raw_yield/eff[0])
                h_corrected_yields[i_split].SetBinError(ct_bin_index, raw_yield_error/eff[0])

            # set labels
            h_corrected_yields[i_split].GetXaxis().SetTitle("#it{c}t (cm)")
            h_corrected_yields[i_split].GetYaxis().SetTitle("d#it{N}/d(#it{c}t) (cm^{-1})")
            for i_bin in range(len(bins))[2:]:
                bin_width = h_corrected_yields[i_split].GetBinWidth(i_bin)
                #print(f"bin: {h_corrected_yields[i_split].GetBinLowEdge(i_bin)}; bin width: {bin_width}")
                bin_content = h_corrected_yields[i_split].GetBinContent(i_bin)
                bin_error = h_corrected_yields[i_split].GetBinError(i_bin)
                h_corrected_yields[i_split].SetBinContent(i_bin, bin_content/bin_width)
                h_corrected_yields[i_split].SetBinError(i_bin, bin_error/bin_width)

            # fit with exponential pdf
            fit_function_expo = ROOT.TF1("expo", "expo", 2, 35)
            if cent_bins[0] == 30:
                fit_function_expo = ROOT.TF1("expo", "expo", 2, 14)
            elif cent_bins[1] == 90:
                fit_function_expo = ROOT.TF1("expo", "expo", 0, 35)
            res = h_corrected_yields[i_split].Fit(fit_function_expo, "QRMIS+")

            # compute lifetime
            tau = -1/fit_function_expo.GetParameter(1)*100/SPEED_OF_LIGHT # ps

            if (res.Prob() > 0.025) and (res.Prob() < 0.975) and (fit_function_expo.GetNDF() == 5):
                systematics_file.cd(f'{cent_bins[0]}_{cent_bins[1]}')
                # h_corrected_yields[i_split].Write()
                lifetime_tmp[i_split] = tau
                h_prob_distribution.Fill(res.Prob())
                h_fit_status.Fill(fit_function_expo.IsValid())

        if (lifetime_tmp[0] > 0) and (lifetime_tmp[1] > 0):
            h_lifetime[0].Fill(lifetime_tmp[0])
            h_lifetime[1].Fill(lifetime_tmp[1])
            h_asymmetry_distribution.Fill(lifetime_tmp[1]-lifetime_tmp[0])
            i_trial+=1
    systematics_file.cd()
    h_asymmetry_distribution.GetXaxis().SetTitle("Asymmetry (ps)")
    h_asymmetry_distribution.GetYaxis().SetTitle("Entries")
    h_asymmetry_distribution.Write()
    h_prob_distribution.GetXaxis().SetTitle("Prob")
    h_prob_distribution.GetYaxis().SetTitle("Entries")
    h_prob_distribution.Write()
    h_lifetime[0].GetXaxis().SetTitle("#tau (ps)")
    h_lifetime[1].GetXaxis().SetTitle("#tau (ps)")
    h_lifetime[0].Write()
    h_lifetime[1].Write()
    h_fit_status.Write()

    del h_asymmetry_distribution
    del h_prob_distribution
    del h_fit_status
    del h_lifetime

systematics_file.Close()
