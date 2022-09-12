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
    hist_eff.Write()
    for val in df_list[1][f"{col_name}MC"]:
        hist_gen.Fill(val)
    hist_gen.Write()

    # compute efficiency and set properties
    hist_eff.Divide(hist_eff, hist_gen, 1, 1, "B")
    if col_name == 'ct':
        hist_eff.GetXaxis().SetTitle('#it{c}t (cm)')
    elif col_name == 'pt':
        hist_eff.GetXaxis().SetTitle('#it{p}_{T} (GeV/#it{c})')
    hist_eff.GetYaxis().SetTitle('Efficiency')
    hist_eff.SetMinimum(0)
    hist_eff.SetDrawOption("pe")
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
MAKE_FEATURES_PLOTS = True
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

PSEUDODATA_PATH = params['PSEUDODATA_PATH']
BKG_PATH = params['BKG_PATH']
MC_SIGNAL_PATH = params['MC_SIGNAL_PATH']
MC_SIGNAL_PATH_GEN = params['MC_SIGNAL_PATH_GEN']
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

    # make plot directory
    if not os.path.isdir(PLOT_DIR):
        os.mkdir(PLOT_DIR)

    # make dataframe directory
    if not os.path.isdir('df'):
        os.mkdir('df')

    for split in SPLIT_LIST:

        split_ineq_sign = '> -999999999.'
        if SPLIT:
            split_ineq_sign = '> 0.5'
            if split == 'antimatter':
                split_ineq_sign = '< 0.5'

        for i_cent_bins in range(len(CENTRALITY_LIST)):
            cent_bins = CENTRALITY_LIST[i_cent_bins]

            df_signal = uproot.open(os.path.expandvars(f"AnalysisResults_{cent_bins[0]}_{cent_bins[1]}.root"))['XiOmegaTree'].arrays(library="pd")

            if MAKE_PRESELECTION_EFFICIENCY and not MAKE_FEATURES_PLOTS and not MAKE_TRAIN_TEST_PLOT:
                ##############################################################
                # PRESELECTION EFFICIENCY
                ##############################################################
                # df_generated = uproot.open(os.path.expandvars(MC_SIGNAL_PATH_GEN))['LambdaTree'].arrays(library="pd")
                
                root_file_presel_eff = ROOT.TFile(f"PreselEff_{cent_bins[0]}_{cent_bins[1]}_fineCt_plots.root", "update")
                df_signal_cent = df_signal.query(f'matter {split_ineq_sign} and (pdg==3334 or pdg==-3334) and mass > 1.645 and mass < 1.70 and pt > 0.5 and pt < 4.5 and isReconstructed and bachBarCosPA < 0.99995 and tpcClV0Pi > 69 and tpcClV0Pr > 69 and tpcClBach > 69 and radius < 25 and radiusV0 < 25 and dcaV0prPV < 2.5 and dcaV0piPV < 2.5 and dcaV0PV < 2.5 and dcaBachPV < 2.5 and eta < 0.8 and eta > -0.8 and isOmega and flag==1') # and (hasTOFhit or hasITSrefit)') # pt cut? (hasTOFhit or hasITSrefit) and 
                df_generated_cent = df_signal.query(
                    f'pdg {split_ineq_sign} and ptMC > 0.5 and ptMC < 4.5 and (pdg==3334 or pdg==-3334) and flag==1') # pt cut?
                #del df_generated

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
                hist_eff_ct.Draw("pe")
                c1.Print(f'{PLOT_DIR}/presel_eff/hPreselEffVsCt_{split}_{cent_bins[0]}_{cent_bins[1]}.pdf')
                hist_eff_pt.Draw("pe")
                c1.Print(f'{PLOT_DIR}/presel_eff/hPreselEffVsPt_{split}_{cent_bins[0]}_{cent_bins[1]}.pdf')

                hist_eff_ct.Write()
                hist_eff_pt.Write()
                root_file_presel_eff.Close()

                del df_signal_cent
                del df_generated_cent

        ######################################################################

    del df_signal

    df_signal = uproot.open(os.path.expandvars(f"AnalysisResults_0_5.root"))['XiOmegaTree'].arrays(library="pd")

    if MAKE_FEATURES_PLOTS and not MAKE_PRESELECTION_EFFICIENCY and not TRAIN:
        ######################################################
        # PLOT FEATURES DISTRIBUTIONS AND CORRELATIONS
        ######################################################

        df_background = uproot.open(os.path.expandvars("/data/mciacco/Omega_PbPb/AnalysisResults-data.root"))['XiOmegaTree'].arrays(library="pd")
        df_prompt_ct = df_signal.query(f'mass > 1.6424500 and mass < 1.7024500 and pt > 0.5 and pt < 4.5 and (pdg==3334 or pdg==-3334) and isReconstructed and tpcClV0Pi > 69 and bachBarCosPA < 0.99995 and tpcClV0Pr > 69 and tpcClBach > 69 and radius < 25 and radiusV0 < 25 and dcaV0prPV < 2.5 and dcaV0piPV < 2.5 and dcaV0PV < 2.5 and dcaBachPV < 2.5 and eta < 0.8 and eta > -0.8 and isOmega') # pt cut?
        df_background_ct = df_background.query(f'mass > 1.6424500 and mass < 1.7024500 and (mass < 1.66 or mass > 1.685) and pt > 0.5 and pt < 4.5 and bachBarCosPA < 0.99995 and tpcClV0Pi > 69 and tpcClV0Pr > 69 and tpcClBach > 69 and radius < 25 and radiusV0 < 25 and dcaV0prPV < 2.5 and dcaV0piPV < 2.5 and dcaV0PV < 2.5 and dcaBachPV < 2.5 and eta < 0.8 and eta > -0.8 and isOmega') # pt cut?
        #print(df_prompt_ct.keys())
        #print(df_background_ct.keys())

        # define tree handlers
        prompt_tree_handler = TreeHandler()
        background_tree_handler = TreeHandler()
        prompt_tree_handler.set_data_frame(df_prompt_ct)
        background_tree_handler.set_data_frame(df_background_ct)
        del df_prompt_ct, df_background_ct

        if not os.path.isdir(f'{PLOT_DIR}/features'):
            os.mkdir(f'{PLOT_DIR}/features')

        leg_labels = ['background', 'signal']
        bins = [40,40,29,40,40,29,40,40,40]
        for i_col, col in enumerate(TRAINING_COLUMNS_LIST):   
            plot_distr = plot_utils.plot_distr(
                [background_tree_handler, prompt_tree_handler],
            col, bins=bins[i_col], labels=leg_labels, log=True, density=True, figsize=(12, 12),
            alpha=0.5, grid=False)
            plt.savefig(f'{PLOT_DIR}/features/FeaturePlots_{col}_{bins[i_col]}.png')
        
        bkg_corr = plot_utils.plot_corr([background_tree_handler], TRAINING_COLUMNS_LIST, ['Background'])
        bkg_corr.set_size_inches(6,6)
        plt.subplots_adjust(left=0.1, bottom=0.06, right=0.99, top=0.96, hspace=0.55, wspace=0.55)
        plt.tight_layout()
        plt.savefig(f'{PLOT_DIR}/features/BackgroundCorrelationMatrix_fullMC.pdf')
        p_corr = plot_utils.plot_corr([prompt_tree_handler], TRAINING_COLUMNS_LIST, ['Signal'])
        p_corr.set_size_inches(6,6)
        plt.tight_layout()
        plt.savefig(f'{PLOT_DIR}/features/PromptCorrelationMatrix_fullMC.pdf')
        plt.close('all')

        ###########################################################
