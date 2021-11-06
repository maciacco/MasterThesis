#!/usr/bin/env python3
import os
import pickle
import warnings
import numpy as np

import ROOT
import uproot
import yaml

ROOT.gROOT.SetBatch()

SPEED_OF_LIGHT = 2.99792458
SPLIT = True
N_TRIALS = 10000

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
SPLIT_LIST = ['all']
if SPLIT:
    SPLIT_LIST = ['antimatter', 'matter']

analysis_results_file = uproot.open(os.path.expandvars(ANALYSIS_RESULTS_PATH))
ratio_file = ROOT.TFile.Open('Ratio.root')

# get centrality selected histogram
cent_counts, cent_edges = analysis_results_file['Centrality_selected;1'].to_numpy()
cent_bin_centers = (cent_edges[:-1]+cent_edges[1:])/2

eff_correction_file = ROOT.TFile.Open('EffAbsCorrection.root')
uncertainty_file = ROOT.TFile.Open('AbsError.root', 'recreate')

h_uncertainty_abs = [ROOT.TH1D(), ROOT.TH1D(), ROOT.TH1D()]

for i_cent_bins in range(len(CENTRALITY_LIST)):
    cent_bins = CENTRALITY_LIST[i_cent_bins]

    # declare uncertainty histogram
    h_uncertainty_abs = ROOT.TH1D(f"fParameterDistribution_{cent_bins[0]}_{cent_bins[1]}",";Ratio;Entries",2500,0,5)

    # get relative systematic error histogram
    func_name = 'BGBW'
    if (cent_bins[1] == 90) or (cent_bins[0] == 0 and cent_bins[1] == 10):
        func_name = 'BlastWave'
    h_relative_sys_antimatter = ROOT.TH1D()
    h_relative_sys_antimatter = eff_correction_file.Get(f"{cent_bins[0]}_{cent_bins[1]}/{func_name}/fRelativeSystematicUncertainty_antimatter_{cent_bins[0]}_{cent_bins[1]}_{func_name}")

    h_relative_sys_matter = ROOT.TH1D()
    h_relative_sys_matter = eff_correction_file.Get(f"{cent_bins[0]}_{cent_bins[1]}/{func_name}/fRelativeSystematicUncertainty_matter_{cent_bins[0]}_{cent_bins[1]}_{func_name}")

    # get number of events
    cent_range_map = np.logical_and(cent_bin_centers > cent_bins[0], cent_bin_centers < cent_bins[1])
    counts_cent_range = cent_counts[cent_range_map]
    evts = np.sum(counts_cent_range)
    print(f'Number of events: {evts}')

    for _ in range(N_TRIALS):
        h_corrected_yields = [ROOT.TH1D(), ROOT.TH1D()]
        for i_split, split in enumerate(SPLIT_LIST):
            #print(f'{i_split} -> {split}')

            nsigma = ROOT.gRandom.Gaus()

            func_name = 'BGBW'
            if (cent_bins[1] == 90) or (cent_bins[0] == 0 and cent_bins[1] == 10):
                func_name = 'BlastWave'

            # list of corrected yields
            ct_bins_tmp = [0]
            ct_bins_tmp += CT_BINS_CENT[i_cent_bins]
            if cent_bins[0]==30:
                ct_bins_tmp = [0, 2, 4, 7, 14, 35]
            #if cent_bins[1] == 90:
            #    ct_bin_tmp = CT_BINS_CENT[i_cent_bins]
            bins = np.array(ct_bins_tmp, dtype=float)
            # print(bins)
            h_corrected_yields[i_split] = ROOT.TH1D()
            h_corrected_yields[i_split] = ratio_file.Get(f'fYields_{split}_{cent_bins[0]}_{cent_bins[1]};1')

            for ct_bins in zip(CT_BINS_CENT[i_cent_bins][:-1], CT_BINS_CENT[i_cent_bins][1:]):
                bin = f'{split}_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}'

                ct_bin_index = h_corrected_yields[i_split].FindBin(ct_bins[0]+0.5)

                syst_abs = 0
                if i_split == 0:
                    syst_abs = h_relative_sys_antimatter.GetBinContent(ct_bin_index)
                else:
                    syst_abs = h_relative_sys_matter.GetBinContent(ct_bin_index)
                corrected_yield = h_corrected_yields[i_split].GetBinContent(ct_bin_index)
                corrected_yield_error = h_corrected_yields[i_split].GetBinError(ct_bin_index)
                correction = 1+nsigma*syst_abs
                #print(f"corrected_yields = {corrected_yield}")
                h_corrected_yields[i_split].SetBinContent(ct_bin_index, corrected_yield*correction)
                h_corrected_yields[i_split].SetBinError(ct_bin_index, corrected_yield_error*correction)

            # # fit with exponential pdf
            # fit_function_expo = ROOT.TF1("expo", "expo", 2, 35)
            # if cent_bins[0] == 30:
            #     fit_function_expo = ROOT.TF1("expo", "expo", 2, 14)
            # elif cent_bins[1] == 90:
            #     fit_function_expo = ROOT.TF1("expo", "expo", 0, 35)
            # res_expo = h_corrected_yields[i_split].Fit(fit_function_expo, "RMIS+")

            # # compute lifetime
            # tau = -1/fit_function_expo.GetParameter(1)*100/SPEED_OF_LIGHT # ps
            # tau_error = fit_function_expo.GetParError(1)*100/SPEED_OF_LIGHT/fit_function_expo.GetParameter(1)/fit_function_expo.GetParameter(1) # ps
            
            # # compute chi2
            # formatted_chi2_lifetime = "{:.2f}".format(fit_function_expo.GetChisquare())

            # # compute yield
            # integral = float()
            # integral_error = float()
            # integral = fit_function_expo.Integral(0, 1.e9)
            # integral_error = fit_function_expo.IntegralError(0, 1.e9, res_expo.GetParams(), res_expo.GetCovarianceMatrix().GetMatrixArray())
            # formatted_integral = "{:.10f}".format(integral/2./evts)
            # formatted_integral_error = "{:.10f}".format(integral_error/evts/2.)

        # ratios
        ROOT.gStyle.SetOptStat(0)
        h_ratio = ROOT.TH1D(h_corrected_yields[0])
        h_ratio.SetName(f'fRatio_{cent_bins[0]}_{cent_bins[1]}')
        h_ratio.Divide(h_corrected_yields[0], h_corrected_yields[1], 1, 1)
        h_ratio.Fit("pol0","Q")
        # h_ratio.Write()

        ratio_val = h_ratio.GetFunction("pol0").GetParameter(0)
        h_uncertainty_abs.Fill(ratio_val)

        del h_corrected_yields
        del h_ratio
    
    h_uncertainty_abs.Write()

ratio_file.Close()
