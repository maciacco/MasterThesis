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
from helpers import significance_error, expected_signal
from hipe4ml import analysis_utils

plt.rcParams.update({'font.size': 13})

SPLIT = True
MAX_EFF = 1.00

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
CT_BINS = params['CT_BINS']
CT_BINS_CENT = params['CT_BINS_CENT']
CENTRALITY_LIST = params['CENTRALITY_LIST']
RANDOM_STATE = params['RANDOM_STATE']
##################################################################

# split matter/antimatter
SPLIT_LIST = ['all']
if SPLIT:
    SPLIT_LIST = ['antimatter', 'matter']

score_eff_arrays_dict = pickle.load(open("file_score_eff_dict", "rb"))
eff_array = np.arange(0.40, MAX_EFF, 0.01)
analysis_results_file = uproot.open(os.path.expandvars(ANALYSIS_RESULTS_PATH))

# get centrality selected histogram
cent_counts, cent_edges = analysis_results_file['Centrality_selected;1'].to_numpy()
cent_bin_centers = (cent_edges[:-1]+cent_edges[1:])/2

# cut dictionary
eff_cut_dict = dict()

for split in SPLIT_LIST:
    split_ineq_sign = '> -999999999'
    if SPLIT:
        split_ineq_sign = '> 0.5'
        if split == 'antimatter':
            split_ineq_sign = '< 0.5'
    for i_cent_bins in range(len(CENTRALITY_LIST)):
        cent_bins = CENTRALITY_LIST[i_cent_bins]
        presel_eff_file = uproot.open(f'PreselEff_{cent_bins[0]}_{cent_bins[1]}_fineCt.root')

        # get number of events
        cent_range_map = np.logical_and(cent_bin_centers > cent_bins[0], cent_bin_centers < cent_bins[1])
        counts_cent_range = cent_counts[cent_range_map]
        evts = np.sum(counts_cent_range)
        print(f'Number of events: {evts}')

        # get preselection efficiency histogram
        presel_eff_counts, presel_eff_edges = presel_eff_file[
            f'fPreselEff_vs_ct_{split}_{cent_bins[0]}_{cent_bins[1]};2'].to_numpy()
        presel_eff_bin_centers = (presel_eff_edges[1:]+presel_eff_edges[:-1])/2
        
        for ct_ in CT_BINS:

            bin_df = f'{cent_bins[0]}_{cent_bins[1]}_{ct_[0]}_{ct_[1]}'
            df_data = pd.read_parquet(f'df_test/all_{bin_df}_data_apply.parquet.gzip')
            df_signal = pd.read_parquet(f'df_test/test_data_{bin_df}_predict.parquet.gzip')
        
            for ct_bins in zip(CT_BINS_CENT[i_cent_bins][:-1], CT_BINS_CENT[i_cent_bins][1:]):
                if (ct_bins[0] < (ct_[0]-0.01)) or (ct_bins[1] > (ct_[1]+0.01)) or ct_bins[0]<0.5:
                    continue
                bin = f'{split}_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}'
                # plot directory
                if not os.path.isdir('plots/significance_scan'):
                    os.mkdir('plots/significance_scan')

                # lists filled inside loop
                significance_list = []
                significance_err_list = []

                df_signal_ct = df_signal.query(f"pdg {split_ineq_sign} and ct >= {ct_bins[0]} and ct < {ct_bins[1]} and mass > 1.6424500 and mass < 1.7024500 and pt > 0.5 and pt < 4.5 and (pdg==3334 or pdg==-3334) and isReconstructed and tpcClV0Pi > 69 and bachBarCosPA < 0.99995 and tpcClV0Pr > 69 and tpcClBach > 69 and radius < 25 and radiusV0 < 25 and dcaV0prPV < 2.5 and dcaV0piPV < 2.5 and dcaV0PV < 2.5 and dcaBachPV < 2.5 and eta < 0.8 and eta > -0.8 and isOmega")
                
                # score-eff relation
                score_array = analysis_utils.score_from_efficiency_array(df_signal_ct["y_true"], df_signal_ct["model_output"], efficiency_selected=eff_array, keep_lower=False)
                delta_eff = 0
                if ct_bins[0] > 6.9:
                    delta_eff = 0.1
                for eff_score in zip(eff_array, score_array):
                    if (ct_bins[0] >0) and (eff_score[0] < (0.49+delta_eff)):
                        continue
                    formatted_eff = "{:.2f}".format(eff_score[0])
                    print(f'processing {bin}: eff = {eff_score[0]:.2f}, score = {eff_score[1]:.2f}...')

                    # select data in the sidebands
                    df_data_sel = df_data.query(f'matter {split_ineq_sign} and model_output > {eff_score[1]}')

                    # make histogram
                    hyp_mass = 1.6725
                    sigma = 0.0015
                    m_min = 1.66 #hyp_mass - 3*sigma
                    m_max = 1.685 #hyp_mass + 3*sigma
                    nBins = 55

                    counts, bin_edges = np.histogram(np.array(df_data_sel['mass']), nBins, range=(1.65,1.7))
                    bin_centers = (bin_edges[1:]+bin_edges[:-1])/2.
                    side_map = np.logical_or(bin_centers < m_min, bin_centers > m_max)
                    mass_map = np.logical_not(side_map)

                    side_bins = bin_centers[side_map]
                    side_counts = counts[side_map]
                    side_errors = np.sqrt(side_counts)

                    # polynomial fit to background
                    pol = np.polynomial.Polynomial.fit(side_bins, side_counts, deg=2)
                    xx_side, yy_side = pol.linspace()  # plot polynomial

                    # compute background
                    pol_integral = pol.integ()
                    bin_size = (1.7-1.65)/nBins
                    bkg = (pol_integral(m_max) - pol_integral(m_min))/bin_size

                    # compute eff = presel_eff * BDT_eff
                    presel_eff_map = np.logical_and(
                        presel_eff_bin_centers > ct_bins[0],
                        presel_eff_bin_centers < ct_bins[1])
                    presel_eff = presel_eff_counts[presel_eff_map]
                    eff = presel_eff * eff_score[0]

                    # compute expected signal
                    sig = expected_signal(cent_bins, ct_bins, eff, evts)[0]
                    if not SPLIT:
                        sig *= 2
                    mass_bins = bin_centers[mass_map]
                    mass_counts = norm.pdf(mass_bins, hyp_mass, sigma)
                    pol_offset = pol(mass_bins)
                    for i, counts in enumerate(mass_counts):
                        mass_counts[i] = float(mass_counts[i]*sig*bin_size+pol_offset[i])

                    xx_mass = np.linspace(m_min, m_max, 100)
                    yy_mass = norm.pdf(xx_mass, hyp_mass, sigma)
                    yy_offset = pol(xx_mass)  # plot polynomial
                    for i, counts in enumerate(yy_mass):
                        yy_mass[i] = float(yy_mass[i]*sig*bin_size+yy_offset[i])

                    # compute significance
                    significance_list.append(sig/np.sqrt(sig+bkg))
                    significance_err_list.append(significance_error(sig, bkg))

                    # plot histograms
                    if not os.path.isdir(f'plots/significance_scan/{bin}'):
                        os.mkdir(f'plots/significance_scan/{bin}')
                    fig = plt.figure()
                    mass_errors = []
                    for counts in mass_counts:
                        mass_errors.append(np.sqrt(counts))
                    plt.errorbar(side_bins, side_counts, side_errors, fmt='o', label='Data', color='blue', ecolor='black')
                    plt.errorbar(mass_bins, mass_counts, mass_errors, fmt='o', label='Pseudodata', color='red', ecolor='black')
                    plt.plot(xx_side, yy_side, label='Background fit', color='green')
                    plt.plot(xx_mass, yy_mass, label='Gaussian model', color='orange')
                    plt.xlabel(r'$M \ (^{3}He + \pi^{-}) \ (\mathrm{GeV}/\it{c}^{2})$')
                    plt.ylabel(r'$\mathrm{Entries}/(2.5 \mathrm{MeV}/\it{c}^{2})$')
                    handles, labels = fig.gca().get_legend_handles_labels()
                    order = [2, 3, 0, 1]
                    plt.legend([handles[idx] for idx in order], [labels[idx] for idx in order], loc='upper right')
                    plt.savefig(f'plots/significance_scan/{bin}/{formatted_eff}_{bin}.png')
                    plt.close('all')

                eff_array_reduced = eff_array[int(9+delta_eff*100):]
                if (ct_bins[0] == 0):
                    eff_array_reduced = eff_array
                significance_array = np.asarray(significance_list)
                significance_err_array = np.asarray(significance_err_list)

                significance_array = significance_array*(eff_array_reduced)
                significance_err_array = significance_err_array*(eff_array_reduced)

                low_limit = significance_array - significance_err_array
                up_limit = significance_array + significance_err_array

                # cut
                cut_eff_index = significance_array == significance_array.max()
                cut_eff = eff_array_reduced[cut_eff_index]

                fig = plt.figure()
                plt.plot(eff_array_reduced, significance_array, 'b', label='Expected Significance')
                plt.fill_between(eff_array_reduced, low_limit, up_limit,
                                facecolor='deepskyblue', label=r'$ \pm 1\sigma$', alpha=0.3)
                plt.vlines(cut_eff[0], 0, significance_array.max(), colors='red', linestyles='dashed',
                        label=f'cut = {"{:.2f}".format(cut_eff[0])}')

                handles, labels = fig.gca().get_legend_handles_labels()
                order = [0, 2, 1]
                plt.legend([handles[idx] for idx in order], [labels[idx] for idx in order], loc='lower left')

                plt.xlabel("BDT Efficiency")
                plt.ylabel("Significance x BDT Efficiency")
                plt.xlim(0.5, MAX_EFF-0.01)
                if ct_bins[0] == 0:
                    plt.xlim(0.1, MAX_EFF-0.01)
                plt.ylim(0.3, up_limit.max()+0.3)

                plt.savefig(f'plots/significance_scan/{bin}.png')
                plt.close('all')

                eff_cut_dict[bin] = cut_eff[0]
                
                # get the systematic variation range
                i_eff = 1
                while i_eff < 11:
                    tmp_eff_var_right = cut_eff[0] + i_eff/100
                    tmp_eff_var_left = cut_eff[0] - i_eff/100
                    
                    # get indices in efficiency array
                    i_eff_var_right = np.where(eff_array_reduced == tmp_eff_var_right)[0]
                    i_eff_var_left = np.where(eff_array_reduced == tmp_eff_var_left)[0]
                    
                    # find corresponding significance
                    tmp_signif_left = significance_array[i_eff_var_right] / tmp_eff_var_right
                    tmp_signif_right = significance_array[i_eff_var_left] / tmp_eff_var_left
                    
                    # if one of the two is smaller than 3, exit
                    if tmp_signif_left<3. or tmp_signif_right<3. :
                        break
                    
                    else:
                        i_eff += 1

                bin_range = f'{split}_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}_range'
                eff_cut_dict[bin_range] = i_eff
                print(f'BDT efficiency cut variation range: +/-{i_eff-1}%')

pickle.dump(eff_cut_dict, open("file_eff_cut_dict", "wb"))
