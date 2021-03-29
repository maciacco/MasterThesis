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

DATA_PATH = params['DATA_PATH']
CT_BINS = params['CT_BINS']
CENTRALITY_LIST = params['CENTRALITY_LIST']
RANDOM_STATE = params['RANDOM_STATE']
##################################################################

# split matter/antimatter
SPLIT_LIST = ['']
if SPLIT:
    SPLIT_LIST = ['antimatter', 'matter']

score_eff_arrays_dict = pickle.load(open("file_score_eff_dict", "rb"))
eff_array = np.arange(0.10, 0.91, 0.01)

for split in SPLIT_LIST:
    for i_cent_bins in range(len(CENTRALITY_LIST)):
        cent_bins = CENTRALITY_LIST[i_cent_bins]
        for ct_bins in zip(CT_BINS[i_cent_bins][:-1], CT_BINS[i_cent_bins][1:]):

            bin = f'{split}_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}'
            df_data = pd.read_parquet(f'df/{bin}')

            # plot directory
            if not os.path.isdir('plots/significance_scan'):
                os.mkdir('plots/significance_scan')
            for eff_score in zip(eff_array, score_eff_arrays_dict[bin]):
                if (ct_bins[0] > 0) and (eff_score[0] < 0.50):
                    continue
                formatted_eff = "{:.2f}".format(eff_score[0])
                print(f'processing {bin}: eff = {eff_score[0]:.2f}, score = {eff_score[1]:.2f}...')

                # select data in the sidebands
                df_data_sel = df_data.query(f'model_output > {eff_score[1]}')

                # make histogram
                hyp_mass = 2.991
                sigma = 0.0015
                m_min = hyp_mass - 3*sigma
                m_max = hyp_mass + 3*sigma
                nBins = 32

                counts, bin_edges = np.histogram(np.array(df_data_sel['m']), nBins, range=(2.96, 3.04))
                bin_centers = (bin_edges[1:]+bin_edges[:-1])/2.
                side_map = np.logical_or(bin_centers < m_min, bin_centers > m_max)
                mass_map = np.logical_not(side_map)

                side_bins = bin_centers[side_map]
                side_counts = counts[side_map]
                side_errors = np.sqrt(side_counts)

                # polynomial fit to background
                pol = np.polynomial.Polynomial.fit(side_bins, side_counts, deg=2)
                xx, yy = pol.linspace()

                # plot histograms
                if not os.path.isdir(f'plots/significance_scan/{bin}'):
                    os.mkdir(f'plots/significance_scan/{bin}')
                plt.errorbar(side_bins, side_counts, side_errors, fmt='o')
                plt.plot(xx, yy)
                plt.xlabel('Invariant mass (GeV/c^2)')
                plt.ylabel('Entries')
                plt.savefig(f'plots/significance_scan/{bin}/{formatted_eff}_{bin}.png')
                plt.close('all')
