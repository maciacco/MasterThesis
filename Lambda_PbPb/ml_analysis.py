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

parser = argparse.ArgumentParser(prog='ml_analysis', allow_abbrev=True)
parser.add_argument('-split', action='store_true')
parser.add_argument('-dotraining', action='store_true')
parser.add_argument('-generate', action='store_true')
parser.add_argument('-mergecentrality', action='store_true')
parser.add_argument('-eff', action='store_true')
parser.add_argument('-train', action='store_true')
parser.add_argument('-computescoreff', action='store_true')
parser.add_argument('-application', action='store_true')
args = parser.parse_args()

SPLIT = args.split
SAMPLE_SIDEBANDS = args.generate
MAX_EFF = 1.00
DUMP_HYPERPARAMS = True

# training
PLOT_DIR = 'plots'
MAKE_PRESELECTION_EFFICIENCY = args.eff
MAKE_TRAIN_TEST_PLOT = True
OPTIMIZE = False
OPTIMIZED = True
TRAIN = False #args.train
COMPUTE_SCORES_FROM_EFF = args.computescoreff
TRAINING = not args.application and (COMPUTE_SCORES_FROM_EFF or TRAIN)
MERGE_CENTRALITY = args.mergecentrality
CREATE_TRAIN_TEST = False

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
BKG_PATH = params['BKG_PATH']
MC_SIGNAL_PATH = params['MC_SIGNAL_PATH']
MC_SIGNAL_PATH_GEN = params['MC_SIGNAL_PATH_GEN']
CT_BINS = params['CT_BINS']
CT_BINS_APPLY = params['CT_BINS_APPLY']
CT_BINS_CENT = params['CT_BINS_CENT']
PT_BINS = params['PT_BINS']
CENTRALITY_LIST = params['CENTRALITY_LIST']
TRAINING_COLUMNS_LIST = params['TRAINING_COLUMNS']
RANDOM_STATE = params['RANDOM_STATE']
HYPERPARAMS = params['HYPERPARAMS']
HYPERPARAMS_RANGES = params['HYPERPARAMS_RANGES']
##################################################################

USE_PD = True

# split matter/antimatter
SPLIT_LIST = ['all']
if SPLIT:
    SPLIT_LIST = ['antimatter', 'matter']

if SAMPLE_SIDEBANDS and not TRAINING and not APPLICATION:

    if USE_PD:
        df_data = uproot.open(os.path.expandvars("../data/Lambda_PbPb/data.root"))['LambdaTree'].arrays(library="pd")
        df_background = df_data.sample(frac=0.05)
        df_data = df_data.drop(df_background.index)
        df_background.to_parquet('df/background_dataset', compression='gzip')
        df_data.to_parquet('df/data_dataset', compression='gzip')

    else:
        df_data = ROOT.RDataFrame("LambdaTree","../data/Lambda_PbPb/data.root")
        df_index = df_data.Define("index","gRandom->Rndm()+mass-mass")
        df_index.Filter("index < 0.05").Snapshot("LambdaTree","../data/Lambda_PbPb/dataBackground.root")
        df_index.Filter("index > 0.05").Snapshot("LambdaTree","../data/Lambda_PbPb/dataSample.root")


if TRAINING:

    df_signal = uproot.open(os.path.expandvars(MC_SIGNAL_PATH))['LambdaTree'].arrays(library="pd")
    df_background = pd.DataFrame()
    if USE_PD:
        df_background = pd.read_parquet('df/background_dataset')
    else:
        df_background = uproot.open(os.path.expandvars(f'{BKG_PATH}'))['LambdaTree'].arrays(library="pd")

    # make plot directory
    if not os.path.isdir(PLOT_DIR):
        os.mkdir(PLOT_DIR)

    # make dataframe directory
    if not os.path.isdir('df'):
        os.mkdir('df')

    score_eff_arrays_dict = dict()

    for ct_bins in CT_BINS:

        if CREATE_TRAIN_TEST and TRAIN:
            df_signal_ct = df_signal.query(f'ct > {ct_bins[0]} and ct < {ct_bins[1]} and pt > 0.5 and pt < 3 and isReconstructed and tpcClV0Pi > 69 and tpcClV0Pr > 69 and radius > 3')
            df_background_ct = df_background.query(f'ct > {ct_bins[0]} and ct < {ct_bins[1]} and pt > 0.5 and pt < 3 and ( mass < 1.105 or mass > 1.13 ) and tpcClV0Pi > 69 and tpcClV0Pr > 69 and radius > 3')

            # define tree handlers
            signal_tree_handler = TreeHandler()
            background_tree_handler = TreeHandler()
            signal_tree_handler.set_data_frame(df_signal_ct)
            background_tree_handler.set_data_frame(df_background_ct)
            del df_signal_ct, df_background_ct

            # split data into training and test set
            train_test_data = train_test_generator([signal_tree_handler, background_tree_handler], [
                1, 0], test_size=0.5, random_state=RANDOM_STATE)
            train_test_data[0]['y_true'] = train_test_data[1]
            train_test_data[2]['y_true'] = train_test_data[3]
            train_test_data[0].to_parquet(f'df/train_data_{ct_bins[0]}_{ct_bins[1]}.parquet.gzip',compression='gzip')
            train_test_data[2].to_parquet(f'df/test_data_{ct_bins[0]}_{ct_bins[1]}.parquet.gzip',compression='gzip')
            continue
        else:
            train_test_data = [pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame()]
            train_test_data[0] = pd.read_parquet(f'df/train_test/train_data_{ct_bins[0]}_{ct_bins[1]}.parquet.gzip')
            train_test_data[2] = pd.read_parquet(f'df/train_test/test_data_{ct_bins[0]}_{ct_bins[1]}.parquet.gzip')


        for split in SPLIT_LIST:
            split_ineq_sign = '> -0.1'
            if SPLIT:
                split_ineq_sign = '> 0.5'
                if split == 'antimatter':
                    split_ineq_sign = '< 0.5'

            for i_cent_bins in range(len(CENTRALITY_LIST)):
                cent_bins = CENTRALITY_LIST[i_cent_bins]
                
                bin = f'{split}_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}'
                ##############################################################
                # TRAINING AND TEST SET PREPARATION
                ##############################################################

                # features plot
                leg_labels = ['background', 'signal']

                model_clf = xgb.XGBClassifier(use_label_encoder=False)
                model_hdl = ModelHandler(model_clf, TRAINING_COLUMNS_LIST)
                model_hdl.set_model_params(HYPERPARAMS)

                # hyperparameters optimization and model training
                if not os.path.isdir('models'):
                    os.mkdir('models')
                bin_model = bin
                if MERGE_CENTRALITY:
                    bin_model = f'all_0_90_{ct_bins[0]}_{ct_bins[1]}'

                if OPTIMIZE and TRAIN:
                    model_hdl.optimize_params_bayes(train_test_data, HYPERPARAMS_RANGES,
                                                    'roc_auc', nfold=5, init_points=10, n_iter=10, njobs=-1)

                isModelTrained = os.path.isfile(f'models/{bin_model}_trained')
                print(f'isModelTrained {bin_model}: {isModelTrained}')
                if TRAIN and not isModelTrained:
                    print(
                    f'Number of candidates ({split}) for training in {ct_bins[0]} <= ct < {ct_bins[1]} cm: {len(train_test_data[0])}')
                    print(
                    f'signal candidates: {np.count_nonzero(train_test_data[1] == 1)}; background candidates: {np.count_nonzero(train_test_data[1] == 0)}; n_cand_bkg / n_cand_signal = {np.count_nonzero(train_test_data[1] == 0) / np.count_nonzero(train_test_data[1] == 1)}')
                    model_hdl.train_test_model(train_test_data, return_prediction=True)
                    model_file_name = str(f'models/{bin_model}_trained')
                    if OPTIMIZE:
                        model_file_name = str(f'models/{bin_model}_optimized_trained')
                    model_hdl.dump_model_handler(model_file_name)
                elif COMPUTE_SCORES_FROM_EFF and isModelTrained:
                    if OPTIMIZED:
                        model_hdl.load_model_handler(f'models/{bin_model}_trained')
                    else:
                        model_hdl.load_model_handler(f'models/{bin_model}_trained')
                else:
                    continue

                ct_bins_df_index = int(ct_bins[0]/5 -1)
                for ct_bins_df in zip(CT_BINS_APPLY[i_cent_bins][ct_bins_df_index][:-1], CT_BINS_APPLY[i_cent_bins][ct_bins_df_index][1:]):
                    bin_df = f'{split}_{cent_bins[0]}_{cent_bins[1]}_{ct_bins_df[0]}_{ct_bins_df[1]}'
                    
                    # get only centrality selected
                    train_test_data_cent = [pd.DataFrame(), [], pd.DataFrame(), []]
                    train_test_data_cent[0] = train_test_data[0].query(f'matter {split_ineq_sign} and centrality > {cent_bins[0]} and centrality < {cent_bins[1]} and ct >= {ct_bins_df[0]} and ct < {ct_bins_df[1]}')
                    train_test_data_cent[2] = train_test_data[2].query(f'matter {split_ineq_sign} and centrality > {cent_bins[0]} and centrality < {cent_bins[1]} and ct >= {ct_bins_df[0]} and ct < {ct_bins_df[1]}')
                    train_test_data_cent[1] = train_test_data_cent[0]['y_true']
                    train_test_data_cent[3] = train_test_data_cent[2]['y_true']

                    # get predictions for training and test sets
                    test_y_score = model_hdl.predict(train_test_data_cent[2])
                    train_y_score = model_hdl.predict(train_test_data_cent[0])

                    # second condition needed because of issue with Qt libraries
                    if MAKE_TRAIN_TEST_PLOT and not MAKE_PRESELECTION_EFFICIENCY:
                        if not os.path.isdir(f'{PLOT_DIR}/train_test_out'):
                            os.mkdir(f'{PLOT_DIR}/train_test_out')
                        plot_utils.plot_output_train_test(model_hdl, train_test_data_cent,
                                                        logscale=True, density=True, labels=leg_labels)
                        plt.savefig(f'{PLOT_DIR}/train_test_out/{bin_df}_out.pdf')

                        plot_utils.plot_feature_imp(train_test_data_cent[0], train_test_data_cent[1], model_hdl)
                        plt.savefig(f'{PLOT_DIR}/train_test_out/feature_imp_training_{bin_df}.pdf')
                        plot_utils.plot_roc_train_test(
                            train_test_data_cent[3],
                            test_y_score, train_test_data_cent[1],
                            train_y_score, labels=leg_labels)
                        plt.savefig(f'{PLOT_DIR}/train_test_out/roc_train_test_{bin_df}.pdf')
                        plt.close('all')

                    if COMPUTE_SCORES_FROM_EFF:
                        # get scores corresponding to BDT efficiencies using test set
                        eff_array = np.arange(0.10, MAX_EFF, 0.01)
                        score_array = analysis_utils.score_from_efficiency_array(
                            train_test_data_cent[3], test_y_score, efficiency_selected=eff_array, keep_lower=False)
                        score_eff_arrays_dict[bin_df] = score_array

                        # write test set data frame
                        train_test_data_cent[2]['model_output'] = test_y_score
                        train_test_data_cent[2]['y_true'] = train_test_data_cent[3]
                        train_test_data_cent_tmp = train_test_data_cent[2].query(f'y_true > 0.5 and ct >= {ct_bins_df[0]} and ct < {ct_bins_df[1]}')
                        train_test_data_cent_tmp.to_parquet(f'df/mc_{bin_df}', compression='gzip')

                # get the model hyperparameters
                if DUMP_HYPERPARAMS and TRAIN:
                    if not os.path.isdir('hyperparams'):
                        os.mkdir('hyperparams')
                    model_params_dict = model_hdl.get_model_params()
                    with open(f'hyperparams/model_params_{bin}.yml', 'w') as outfile:
                        yaml.dump(model_params_dict, outfile, default_flow_style=False)

                    # save roc-auc
                del train_test_data_cent
                ##############################################################

    if COMPUTE_SCORES_FROM_EFF:
        pickle.dump(score_eff_arrays_dict, open("file_score_eff_dict", "wb"))

# apply model to data
if APPLICATION:
    if not os.path.isdir('df'):
        os.mkdir('df')
    score_eff_arrays_dict = pickle.load(open("file_score_eff_dict", "rb"))

    for split in SPLIT_LIST:

        split_ineq_sign = '> -0.1'
        if SPLIT:
            split_ineq_sign = '> 0.5'
            if split == 'antimatter':
                split_ineq_sign = '< 0.5'

        for i_cent_bins in range(len(CENTRALITY_LIST)):
            cent_bins = CENTRALITY_LIST[i_cent_bins]

            for ct_bins in zip(CT_BINS_CENT[i_cent_bins][:-1], CT_BINS_CENT[i_cent_bins][1:]):
                # if (ct_bins[0] < 19.4):
                #     continue
                bin = f'{split}_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}'

                ct_bins_df_index = int(ct_bins[0]/5 -1)
                model_hdl = ModelHandler()
                bin_model = f'{split}_{cent_bins[0]}_{cent_bins[1]}_{CT_BINS_APPLY[i_cent_bins][ct_bins_df_index][0]}_{CT_BINS_APPLY[i_cent_bins][ct_bins_df_index][-1]}'
                if MERGE_CENTRALITY:
                    bin_model = f'all_0_90_{CT_BINS_APPLY[i_cent_bins][ct_bins_df_index][0]}_{CT_BINS_APPLY[i_cent_bins][ct_bins_df_index][-1]}'
                if OPTIMIZED:
                    model_hdl.load_model_handler(f'models/{bin_model}_trained') # *_optimized_trained
                else:
                    model_hdl.load_model_handler(f'models/{bin_model}_trained')

                eff_array = np.arange(0.10, MAX_EFF, 0.01)
                if USE_PD:
                    df_data = pd.read_parquet('df/data_dataset')
                    df_data_cent = df_data.query(
                    f'matter {split_ineq_sign} and centrality > {cent_bins[0]} and centrality < {cent_bins[1]} and pt > 0.5 and pt < 3 and ct > {ct_bins[0]} and ct < {ct_bins[1]} and tpcClV0Pi > 69 and tpcClV0Pr > 69 and radius > 3')
                    del df_data

                    data_y_score = model_hdl.predict(df_data_cent)
                    df_data_cent['model_output'] = data_y_score

                    df_data_cent = df_data_cent.query(f'model_output > {score_eff_arrays_dict[bin][len(eff_array)-1]}')
                    df_data_cent.to_parquet(f'df/{bin}.parquet.gzip', compression='gzip')
                else:
                    df_data = TreeHandler()
                    df_data.get_handler_from_large_file(DATA_PATH, "LambdaTree",
                        preselection=f'matter {split_ineq_sign} and centrality > {cent_bins[0]} and centrality < {cent_bins[1]} and pt > 0.5 and pt < 3 and ct > {ct_bins[0]} and ct < {ct_bins[1]}',
                        max_workers=8)

                    df_data.apply_model_handler(model_hdl)
                    df_data.apply_preselections(f'model_output > {score_eff_arrays_dict[bin][len(eff_array)-1]}')
                    df_data.write_df_to_parquet_files(bin,"df/")
