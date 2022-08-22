#!/usr/bin/env python3
import os
import pickle
import warnings
import argparse

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import ROOT
import uproot
import xgboost as xgb
import yaml
from hipe4ml import analysis_utils, plot_utils
from hipe4ml.analysis_utils import train_test_generator
from hipe4ml.model_handler import ModelHandler
from hipe4ml.tree_handler import TreeHandler
from sklearn.model_selection import train_test_split
from sklearn.utils.class_weight import compute_sample_weight

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
MAX_SCORE = 1.00
MAX_EFF = 1.00
USE_PD = True
DUMP_HYPERPARAMS = True
USE_REAL_DATA = False

PRODUCE_DATASETS = args.generate
TRAINING = False
use_gpu = False

# avoid pandas warning
warnings.simplefilter(action='ignore', category=FutureWarning)
ROOT.gROOT.SetBatch()
ROOT.EnableImplicitMT(10)

# training
PLOT_DIR = 'plots_test_exact'
MAKE_PRESELECTION_EFFICIENCY = args.eff
MAKE_TRAIN_TEST_PLOT = True
OPTIMIZE = False
OPTIMIZED = False
TRAIN = args.dotraining
COMPUTE_SCORES_FROM_EFF = args.computescoreff
TRAINING = args.train and (COMPUTE_SCORES_FROM_EFF or TRAIN)
MERGE_CENTRALITY = args.mergecentrality
CREATE_TRAIN_TEST = True

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
CT_BINS_APPLY = params['CT_BINS_APPLY']
CT_BINS_CENT = params['CT_BINS_CENT']
PT_BINS = params['PT_BINS']
CENTRALITY_LIST = params['CENTRALITY_LIST']
CENTRALITY_LIST_TRAINING = params['CENTRALITY_LIST_TRAINING']
TRAINING_COLUMNS_LIST = params['TRAINING_COLUMNS']
RANDOM_STATE = params['RANDOM_STATE']
HYPERPARAMS = params['HYPERPARAMS']
HYPERPARAMS_RANGES = params['HYPERPARAMS_RANGES']
##################################################################

ROOT.gInterpreter.ProcessLine(".L help.h+")
# split matter/antimatter
SPLIT_LIST = ['all']
if SPLIT:
    SPLIT_LIST = ['antimatter', 'matter']

if PRODUCE_DATASETS and not TRAIN:

    for i_cent_bins in range(len(CENTRALITY_LIST_TRAINING)):
        cent_bins = CENTRALITY_LIST_TRAINING[i_cent_bins]

        # get sidebands (both for training and as pseudo-data)
        df_data = ROOT.RDataFrame("LambdaTree","/data/mciacco/LambdaPrompt_PbPb/AnalysisResults_data.root")
        df_index = df_data.Define("index_1","gRandom->Rndm()+mass-mass")
        df_index.Filter(f"(mass < 1.105 || mass > 1.13) && centrality > {cent_bins[0]} && centrality < {cent_bins[1]}").Snapshot("LambdaTree",f"/data/mciacco/LambdaPrompt_PbPb/trainingBackground_{cent_bins[0]}_{cent_bins[1]}.root")

if TRAINING:

    # make plot directory
    if not os.path.isdir(PLOT_DIR):
        os.mkdir(PLOT_DIR)

    # make dataframe directory
    if not os.path.isdir('df_test'):
        os.mkdir('df_test')

    score_eff_arrays_dict = dict()

    for i_cent_bins in range(len(CENTRALITY_LIST_TRAINING)):
        cent_bins = CENTRALITY_LIST_TRAINING[i_cent_bins]

        df_signal = uproot.open(os.path.expandvars(f"/data/mciacco/LambdaPrompt_PbPb/AnalysisResults_reweight_BW_{cent_bins[0]}_{cent_bins[1]}.root"))['LambdaTree'].arrays(library="pd")
        df_signal_prompt = df_signal.query("flag==1")
        df_signal_non_prompt_omega = df_signal.query("flag==4")
        df_signal_non_prompt = uproot.open(os.path.expandvars(f"/data/mciacco/LambdaPrompt_PbPb/AnalysisResults_Xi0_reweight_BW_{cent_bins[0]}_{cent_bins[1]}.root"))['LambdaTree'].arrays(library="pd")
        df_background = uproot.open(os.path.expandvars(f"/data/mciacco/LambdaPrompt_PbPb/trainingBackground_{cent_bins[0]}_{cent_bins[1]}.root"))['LambdaTree'].arrays(library="pd")

        for ct_bins in CT_BINS:

            for split in SPLIT_LIST:
                split_ineq_sign = '> -0.1'
                if SPLIT:
                    split_ineq_sign = '> 0.5'
                    if split == 'antimatter':
                        split_ineq_sign = '< 0.5'

                bin = f'{split}_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}'

                train_test_data = [pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame()]
                if CREATE_TRAIN_TEST and (COMPUTE_SCORES_FROM_EFF or TRAIN):
                    df_prompt_ct = df_signal_prompt.query(f'ct > {ct_bins[0]} and ct < {ct_bins[1]} and pt > 0.5 and pt < 3.5 and isReconstructed and (flag==1) and tpcClV0Pi > 69 and tpcClV0Pr > 69 and radius > 3 and radius < 100 and dcaPrPV < 20 and dcaPiPV < 20 and dcaV0PV < 10 and eta < 0.8 and eta > -0.8')
                    df_nonprompt_ct_Xi = df_signal_non_prompt.query(f'ct > {ct_bins[0]} and ct < {ct_bins[1]} and pt > 0.5 and pt < 3.5 and isReconstructed and (flag==2) and tpcClV0Pi > 69 and tpcClV0Pr > 69 and radius > 3 and radius < 100 and dcaPrPV < 20 and dcaPiPV < 20 and dcaV0PV < 10 and eta < 0.8 and eta > -0.8')
                    df_nonprompt_ct_Omega = df_signal_non_prompt_omega.query(f'ct > {ct_bins[0]} and ct < {ct_bins[1]} and pt > 0.5 and pt < 3.5 and isReconstructed and (flag==4) and tpcClV0Pi > 69 and tpcClV0Pr > 69 and radius > 3 and radius < 100 and dcaPrPV < 20 and dcaPiPV < 20 and dcaV0PV < 10 and eta < 0.8 and eta > -0.8')
                    df_background_ct = df_background.query(f'ct > {ct_bins[0]} and ct < {ct_bins[1]} and (mass < 1.102 or mass > 1.13) and pt > 0.5 and pt < 3.5 and tpcClV0Pi > 69 and tpcClV0Pr > 69 and radius > 3 and radius < 100 and dcaPrPV < 20 and dcaPiPV < 20 and dcaV0PV < 10 and eta < 0.8 and eta > -0.8')

                    n_background = df_background_ct.shape[0]
                    dfs_nonprompt_ct = [df_nonprompt_ct_Xi,df_nonprompt_ct_Omega]
                    df_nonprompt_ct = pd.concat(dfs_nonprompt_ct)
                    n_non_prompt = df_nonprompt_ct.shape[0]
                    if 5*n_non_prompt/n_background < 1:
                        df_background_ct = df_background_ct.sample(frac=5*n_non_prompt/n_background)

                    # define tree handlers
                    prompt_tree_handler = TreeHandler()
                    nonprompt_tree_handler = TreeHandler()
                    background_tree_handler = TreeHandler()
                    prompt_tree_handler.set_data_frame(df_prompt_ct)
                    nonprompt_tree_handler.set_data_frame(df_nonprompt_ct)
                    background_tree_handler.set_data_frame(df_background_ct)
                    del df_prompt_ct, df_nonprompt_ct, df_background_ct, df_nonprompt_ct_Xi, df_nonprompt_ct_Omega

                    # split data into training and test set
                    train_test_data = train_test_generator([background_tree_handler, nonprompt_tree_handler, prompt_tree_handler], [
                        0, 1, 2], test_size=0.5, random_state=RANDOM_STATE)
                    train_test_data[0]['y_true'] = train_test_data[1]
                    train_test_data[2]['y_true'] = train_test_data[3]
                    train_test_data[0].to_parquet(f'df_test/train_data_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}.parquet.gzip',compression='gzip')
                    train_test_data[2].to_parquet(f'df_test/test_data_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}.parquet.gzip',compression='gzip')
                    # continue
                else:
                    train_test_data[0] = pd.read_parquet(f'df_test/train_data_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}.parquet.gzip')
                    train_test_data[2] = pd.read_parquet(f'df_test/test_data_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}.parquet.gzip')
                    train_test_data[1] = train_test_data[0]['y_true']
                    train_test_data[3] = train_test_data[2]['y_true']
                ##############################################################
                # TRAINING AND TEST SET PREPARATION
                ##############################################################

                # features plot
                leg_labels = ['background', 'non_prompt', 'prompt']
                
                model_clf = xgb.XGBClassifier()
                if use_gpu:
                    model_clf = xgb.XGBClassifier(use_label_encoder=False, tree_method="gpu_hist", gpu_id=0)
                else:
                    model_clf = xgb.XGBClassifier(use_label_encoder=False, n_jobs=2, tree_method="hist")
                model_hdl = ModelHandler(model_clf, TRAINING_COLUMNS_LIST)
                model_hdl.set_model_params(HYPERPARAMS)

                # hyperparameters optimization and model training
                print(
                f'Number of candidates ({split}) for training in {ct_bins[0]} <= ct < {ct_bins[1]} cm: {len(train_test_data[0])}')
                print(
                f'prompt candidates: {np.count_nonzero(train_test_data[1] == 2)}; non-prompt candidates: {np.count_nonzero(train_test_data[1] == 1)}; background candidates: {np.count_nonzero(train_test_data[1] == 0)}; n_cand_bkg / n_cand_signal = {np.count_nonzero(train_test_data[1] == 0) / np.count_nonzero(train_test_data[1] == 2)}')
                
                if not os.path.isdir('models_test'):
                    os.mkdir('models_test')
                bin_model = bin
                if MERGE_CENTRALITY:
                    bin_model = f'all__{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}'

                if OPTIMIZE and TRAIN:
                    model_hdl.optimize_params_optuna(train_test_data, HYPERPARAMS_RANGES,
                                                    'roc_auc_ovr', nfold=5, timeout=300)

                isModelTrained = os.path.isfile(f'models_test/{bin_model}_trained')
                print(f'isModelTrained {bin_model}: {isModelTrained}')
                if TRAIN and not isModelTrained:
                    # weights={0:2,1:6,2:3}
                    # sample_weights = compute_sample_weight(class_weight=weights,y=train_test_data[0]['y_true'])
                    model_hdl.train_test_model(train_test_data, multi_class_opt="ovr", return_prediction=True, output_margin=False) #, sample_weight=sample_weights)
                    model_file_name = str(f'models_test/{bin_model}_trained')
                    if OPTIMIZE:
                        model_file_name = str(f'models_test/{bin_model}_optimized_trained')
                    model_hdl.dump_model_handler(model_file_name)
                elif COMPUTE_SCORES_FROM_EFF and isModelTrained:
                    print('Model trained...')
                    if OPTIMIZED:
                        model_hdl.load_model_handler(f'models_test/{bin_model}_trained')
                    else:
                        model_hdl.load_model_handler(f'models_test/{bin_model}_trained')
                else:
                    continue
                model_file_name = str(f'models_test/{bin_model}.model')
                model_hdl.dump_original_model(model_file_name,True)

                ct_bins_df_index = int(ct_bins[0]/5)
                
                bin_df = f'{split}_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}'
                print(bin_df)
                # get only centrality selected
                train_test_data_cent = [pd.DataFrame(), [], pd.DataFrame(), []]
                train_test_data_cent[0] = train_test_data[0].query(f'matter {split_ineq_sign} and ct >= {ct_bins[0]} and ct < {ct_bins[1]}')
                train_test_data_cent[2] = train_test_data[2].query(f'matter {split_ineq_sign} and ct >= {ct_bins[0]} and ct < {ct_bins[1]}')
                train_test_data_cent[1] = train_test_data_cent[0]['y_true']
                train_test_data_cent[3] = train_test_data_cent[2]['y_true']

                # get predictions for training and test sets
                print(train_test_data_cent[2])
                if (ct_bins[0]==0):
                    continue
                test_y_score = model_hdl.predict(train_test_data_cent[2], output_margin=False)
                print(test_y_score)
                train_y_score = model_hdl.predict(train_test_data_cent[0], output_margin=False)
                train_test_data_cent[0].loc[:,'model_output_prompt'] = train_y_score[:,2]
                train_test_data_cent[2].loc[:,'model_output_prompt'] = test_y_score[:,2]
                train_test_data_cent[0].loc[:,'model_output_non_prompt'] = train_y_score[:,1]
                train_test_data_cent[2].loc[:,'model_output_non_prompt'] = test_y_score[:,1]
                train_test_data_cent[0].loc[:,'model_output_background'] = train_y_score[:,0]
                train_test_data_cent[2].loc[:,'model_output_background'] = test_y_score[:,0]

                # write
                train_test_data_cent[0].to_parquet(f'df_test/train_data_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}_predict.parquet.gzip',compression='gzip')
                train_test_data_cent[2].to_parquet(f'df_test/test_data_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}_predict.parquet.gzip',compression='gzip')
                
                # second condition needed because of issue with Qt libraries
                if MAKE_TRAIN_TEST_PLOT and not MAKE_PRESELECTION_EFFICIENCY:
                    if not os.path.isdir(f'{PLOT_DIR}/train_test_out'):
                        os.mkdir(f'{PLOT_DIR}/train_test_out')
                    out_figs = plot_utils.plot_output_train_test(model_hdl, train_test_data_cent, bins=50,
                                                    logscale=True, density=True, labels=leg_labels, output_margin=False)
                    for i_label, label in enumerate(leg_labels):
                        out_figs[i_label].savefig(f'{PLOT_DIR}/train_test_out/{bin_df}_out_{label}.pdf')

                    feat_imp = plot_utils.plot_feature_imp(train_test_data_cent[0], train_test_data_cent[1], model_hdl)
                    for i_label, label in enumerate(leg_labels):
                       feat_imp[i_label].savefig(f'{PLOT_DIR}/train_test_out/feature_imp_training_{bin_df}_{label}.pdf')
                    feat_imp[3].savefig(f'{PLOT_DIR}/train_test_out/feature_imp_training_{bin_df}_all.pdf')
                    plot_utils.plot_roc_train_test(
                        train_test_data_cent[3],
                        test_y_score, train_test_data_cent[1],
                        train_y_score, labels=leg_labels, multi_class_opt="ovr")
                    plt.savefig(f'{PLOT_DIR}/train_test_out/roc_train_test_{bin_df}.pdf')
                    plt.close('all')

                if COMPUTE_SCORES_FROM_EFF:
                    pass
                    # get scores corresponding to BDT prompt efficiencies using test set
                    eff_selected = np.arange(0.1, MAX_EFF, 0.01)
                    y_truth_tmp = train_test_data_cent[2]['y_true'].apply(lambda x : int(x/2))
                    score = analysis_utils.score_from_efficiency_array(
                        y_truth_tmp, test_y_score[:,0], eff_selected, keep_lower=False)
                    print(y_truth_tmp)
                    print(score)
                    score_eff_arrays_dict[bin_df] = score

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

    if COMPUTE_SCORES_FROM_EFF and TRAIN:
        pickle.dump(score_eff_arrays_dict, open("file_score_eff_dict", "wb"))


# apply model to data
if APPLICATION:
    if not os.path.isdir('df_test'):
        os.mkdir('df_test')
    score_eff_arrays_dict = pickle.load(open("file_score_eff_dict", "rb"))

    for split in SPLIT_LIST:

        split_ineq_sign = '> -0.1'
        if SPLIT:
            split_ineq_sign = '> 0.5'
            if split == 'antimatter':
                split_ineq_sign = '< 0.5'

        for i_cent_bins in range(len(CENTRALITY_LIST)):
            cent_bins = CENTRALITY_LIST[i_cent_bins]
            
            for ct_bins in CT_BINS:
                if (ct_bins[0] < 15 or ct_bins[1] > 40):
                    continue
                bin = f'{split}_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}_data_apply'

                ct_bins_index = int(ct_bins[0]/5 -1)
                model_hdl = ModelHandler()
                bin_model = f'{split}_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[-1]}'
                if MERGE_CENTRALITY:
                    bin_model = f'all_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[-1]}'
                if OPTIMIZED:
                    model_hdl.load_model_handler(f'models_test/{bin_model}_optimized_trained') # *_optimized_trained
                else:
                    print(bin_model)
                    model_hdl.load_model_handler(f'models_test/{bin_model}_trained')

                eff_array = np.arange(0.10, MAX_EFF, 0.01)
                if USE_REAL_DATA:
                    if USE_PD:
                        #df_data = pd.read_parquet('df_test/data_dataset')
                        df_data = uproot.open(f'/data/mciacco/LambdaPrompt_PbPb/test_real_data_0_5.root')['LambdaTree'].arrays(library="pd")
                        #df_data = df_data.append(df_data_r, ignore_index=True)
                        df_data_cent = df_data.query(
                        f'matter {split_ineq_sign} and centrality > {cent_bins[0]} and centrality < {cent_bins[1]} and pt > 0.5 and pt < 3.5 and ct > {ct_bins[0]} and ct < {ct_bins[1]} and tpcClV0Pi > 69 and tpcClV0Pr > 69 and radius > 3 and radius < 100 and dcaPrPV < 20 and dcaPiPV < 20 and dcaV0PV < 10 and eta < 0.8 and eta > -0.8')
                        del df_data

                        data_y_score = model_hdl.predict(df_data_cent, output_margin=False)
                        df_data_cent['model_output_0'] = data_y_score[:,0]
                        df_data_cent['model_output_1'] = data_y_score[:,1]
                        df_data_cent['model_output_2'] = data_y_score[:,2]

                        # df_data_cent = df_data_cent.query(f'model_output > {score_eff_arrays_dict[bin][len(eff_array)-1]}')
                        df_data_cent.to_parquet(f'df_test/{bin}.parquet.gzip', compression='gzip')
                    else:
                        df_data = TreeHandler()
                        df_data.get_handler_from_large_file(f"/data/mciacco/LambdaPrompt_PbPb/test_real_data_0_5.root", "LambdaTree",
                            preselection=f'matter {split_ineq_sign} and centrality > {cent_bins[0]} and centrality < {cent_bins[1]} and pt > 0.5 and pt < 3.5 and ct > {ct_bins[0]} and ct < {ct_bins[1]} and tpcClV0Pi > 69 and tpcClV0Pr > 69 and radius > 3 and radius < 100 and dcaPrPV < 20 and dcaPiPV < 20 and dcaV0PV < 10 and eta < 0.8 and eta > -0.8',
                            max_workers=4, model_handler=model_hdl, output_margin=False)


                        #df_data.apply_model_handler(model_hdl, output_margin=False)
                        #df_data.apply_preselections(f'model_output > {score_eff_arrays_dict[bin][len(eff_array)-1]}')
                        df_data.write_df_to_parquet_files(bin,"df_test/")
                else:
                    bin = f'{split}_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}_mc_apply'
                    df_data = TreeHandler()
                    df_data.get_handler_from_large_file(f"/data/mciacco/LambdaPrompt_PbPb/AnalysisResults_reweight_BW_0_5.root", "LambdaTree",
                        preselection=f'isReconstructed and flag<4 and pdg {split_ineq_sign} and pt > 0.5 and pt < 3.5 and ct > {ct_bins[0]} and ct < {ct_bins[1]} and tpcClV0Pi > 69 and tpcClV0Pr > 69 and radius > 3 and radius < 100 and dcaPrPV < 20 and dcaPiPV < 20 and dcaV0PV < 10 and eta < 0.8 and eta > -0.8',
                        max_workers=4, model_handler=model_hdl, output_margin=False)

                    df_data.write_df_to_parquet_files(bin,"df_test/")
