#!/usr/bin/env python3
import os
import pickle
import warnings
import argparse

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import ROOT
import uproot
import xgboost as xgb
import yaml
from hipe4ml import analysis_utils, plot_utils
from hipe4ml.analysis_utils import train_test_generator
from hipe4ml.model_handler import ModelHandler
from hipe4ml.tree_handler import TreeHandler
from sklearn.model_selection import train_test_split


def presel_eff_hist(df_list, col_name, split, cent_bins, bins):
    # fill histograms (vs. ct and vs. pt)
    bins_array = np.asarray(bins, dtype=float)
    hist_eff = ROOT.TH1F(
        f'fPreselEff_vs_{col_name}_{split}_{cent_bins[0]}_{cent_bins[1]}',
        f'Preselection Efficiency, {split}, {cent_bins[0]}-{cent_bins[1]}%', len(bins) - 1, bins_array)
    hist_gen = ROOT.TH1F('fPreselGen_vs_{col_name}', 'Gen', len(bins)-1, bins_array)

    for val in df_list[0][col_name]:
        hist_eff.Fill(val)
    for val in df_list[1][col_name]:
        hist_gen.Fill(val)

    # compute efficiency and set properties
    hist_eff.Divide(hist_eff, hist_gen, 1, 1, "B")
    if col_name == 'ct':
        hist_eff.GetXaxis().SetTitle('#it{c}t (cm)')
    elif col_name == 'pt':
        hist_eff.GetXaxis().SetTitle('#it{p}_{T} (GeV/#it{c})')
    hist_eff.GetYaxis().SetTitle('Efficiency')
    hist_eff.SetMinimum(0)
    hist_eff.SetDrawOption("histo")
    hist_eff.SetLineWidth(2)

    # return histogram
    return hist_eff


parser = argparse.ArgumentParser(prog='ml_analysis', allow_abbrev=True)
parser.add_argument('-split', action='store_true')
parser.add_argument('-mergecentrality', action='store_true')
parser.add_argument('-eff', action='store_true')
parser.add_argument('-train', action='store_true')
parser.add_argument('-computescoreff', action='store_true')
parser.add_argument('-application', action='store_true')
args = parser.parse_args()

SPLIT = args.split
MAX_EFF = 1.00
DUMP_HYPERPARAMS = False

# training
TRAINING = not args.application
PLOT_DIR = 'plots'
MAKE_PRESELECTION_EFFICIENCY = args.eff
MAKE_FEATURES_PLOTS = False
MAKE_TRAIN_TEST_PLOT = args.train
OPTIMIZE = False
OPTIMIZED = False
TRAIN = args.train
COMPUTE_SCORES_FROM_EFF = args.computescoreff
MERGE_CENTRALITY = args.mergecentrality

# application
APPLICATION = args.application

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
MC_PATH = params['MC_SIGNAL_PATH']
BKG_PATH = params['LS_BACKGROUND_PATH']
CT_BINS = params['CT_BINS']
CT_BINS_CENT = params['CT_BINS_CENT']
PT_BINS = params['PT_BINS']
CENTRALITY_LIST = params['CENTRALITY_LIST']
TRAINING_COLUMNS_LIST = params['TRAINING_COLUMNS']
RANDOM_STATE = params['RANDOM_STATE']
HYPERPARAMS = params['HYPERPARAMS']
HYPERPARAMS_RANGES = params['HYPERPARAMS_RANGES']
##################################################################

# split matter/antimatter
SPLIT_LIST = ['all']
if SPLIT:
    SPLIT_LIST = ['antimatter', 'matter']

if TRAINING:

    df_signal = uproot.open(os.path.expandvars(MC_PATH))['SignalTable'].arrays(library="pd")
    df_background = uproot.open(os.path.expandvars(BKG_PATH))['DataTable'].arrays(library="pd")

    # make plot directory
    if not os.path.isdir(PLOT_DIR):
        os.mkdir(PLOT_DIR)

    # make dataframe directory
    if not os.path.isdir('df'):
        os.mkdir('df')

    for split in SPLIT_LIST:

        split_ineq_sign = '> -0.1'
        if SPLIT:
            split_ineq_sign = '> 0.5'
            if split == 'antimatter':
                split_ineq_sign = '< 0.5'

        for i_cent_bins in range(len(CENTRALITY_LIST)):
            cent_bins = CENTRALITY_LIST[i_cent_bins]

            if MAKE_PRESELECTION_EFFICIENCY and not MAKE_FEATURES_PLOTS and not MAKE_TRAIN_TEST_PLOT:
                ##############################################################
                # PRESELECTION EFFICIENCY
                ##############################################################
                df_generated = uproot.open(os.path.expandvars(MC_PATH))['GenTable'].arrays(library="pd")
                df_signal_cent = df_signal.query(
                    f'Matter {split_ineq_sign} and centrality > {cent_bins[0]} and centrality < {cent_bins[1]} and pt > 2 and pt < 10')
                df_generated_cent = df_generated.query(
                    f'matter {split_ineq_sign} and centrality > {cent_bins[0]} and centrality < {cent_bins[1]} and pt > 2 and pt < 10')
                del df_generated

                # fill histograms (vs. ct and vs. pt)
                hist_eff_ct = presel_eff_hist([df_signal_cent, df_generated_cent], 'ct', split, cent_bins, CT_BINS_CENT[i_cent_bins])
                hist_eff_pt = presel_eff_hist([df_signal_cent, df_generated_cent], 'pt', split, cent_bins, PT_BINS)

                # plot histograms
                if not os.path.isdir(f'{PLOT_DIR}/presel_eff'):
                    os.mkdir(f'{PLOT_DIR}/presel_eff')
                c1 = ROOT.TCanvas()
                ROOT.gStyle.SetOptStat(0)
                c1.cd()
                c1.SetTicks(1,1)
                hist_eff_ct.Draw("histo")
                c1.Print(f'{PLOT_DIR}/presel_eff/hPreselEffVsCt_{split}_{cent_bins[0]}_{cent_bins[1]}.pdf')
                hist_eff_pt.Draw("histo")
                c1.Print(f'{PLOT_DIR}/presel_eff/hPreselEffVsPt_{split}_{cent_bins[0]}_{cent_bins[1]}.pdf')

                root_file_presel_eff = ROOT.TFile("PreselEff.root", "update")
                hist_eff_ct.Write()
                hist_eff_pt.Write()
                root_file_presel_eff.Close()

                del df_signal_cent
                del df_generated_cent

        ######################################################################

    # second condition needed because of issue with Qt libraries
    if MAKE_FEATURES_PLOTS and not MAKE_PRESELECTION_EFFICIENCY and not TRAIN:
        ######################################################
        # PLOT FEATURES DISTRIBUTIONS AND CORRELATIONS
        ######################################################

        df_signal_ct = df_signal.query(f'pt > 2 and pt < 10')
        df_background_ct = df_background.query(f'pt > 2 and pt < 10')

        # define tree handlers
        signal_tree_handler = TreeHandler()
        background_tree_handler = TreeHandler()
        signal_tree_handler.set_data_frame(df_signal_ct)
        background_tree_handler.set_data_frame(df_background_ct)
        del df_signal_ct, df_background_ct

        if not os.path.isdir(f'{PLOT_DIR}/features'):
            os.mkdir(f'{PLOT_DIR}/features')

        leg_labels = ['background', 'signal']
        plot_utils.plot_distr(
            [background_tree_handler, signal_tree_handler],
            TRAINING_COLUMNS_LIST, bins=40, labels=leg_labels, log=True, density=True, figsize=(12, 7),
            alpha=0.3, grid=False)
        plt.subplots_adjust(left=0.06, bottom=0.06, right=0.99, top=0.96, hspace=0.55, wspace=0.55)
        plt.savefig(f'{PLOT_DIR}/features/FeaturePlots')
        plot_utils.plot_corr([background_tree_handler], TRAINING_COLUMNS_LIST, ['background'])
        plt.savefig(f'{PLOT_DIR}/features/BackgroundCorrelationMatrix')
        plot_utils.plot_corr([signal_tree_handler], TRAINING_COLUMNS_LIST, ['signal'])
        plt.savefig(f'{PLOT_DIR}/features/SignalCorrelationMatrix')
        plt.close('all')

        ###########################################################