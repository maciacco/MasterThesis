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
USE_PD = False
DUMP_HYPERPARAMS = True

PRODUCE_DATASETS = args.generate
TRAINING = False

# avoid pandas warning
warnings.simplefilter(action='ignore', category=FutureWarning)
ROOT.gROOT.SetBatch()

# training
PLOT_DIR = 'plots'
MAKE_PRESELECTION_EFFICIENCY = args.eff
MAKE_TRAIN_TEST_PLOT = True
OPTIMIZE = False
OPTIMIZED = True
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
TRAINING_COLUMNS_LIST = params['TRAINING_COLUMNS']
RANDOM_STATE = params['RANDOM_STATE']
HYPERPARAMS = params['HYPERPARAMS']
HYPERPARAMS_RANGES = params['HYPERPARAMS_RANGES']
##################################################################

# split matter/antimatter
SPLIT_LIST = ['all']
if SPLIT:
    SPLIT_LIST = ['antimatter', 'matter']

if PRODUCE_DATASETS and not TRAIN:
    b_list = ["pt","ct","mass","matter","eta","radius","centrality","hasITSrefit","cosPA", "dcaV0tracks", "dcaPiPV", "dcaPrPV", "dcaV0PV", "tpcNsigmaPr", "tpcNsigmaPi", "tpcClV0Pr", "tpcClV0Pi"]
    branch_list = ROOT.std.vector('std::string')()
    for b in b_list:
        branch_list.push_back(b)
    b_list_mc = ["flag","pt","ct","ptMC","ctMC","matter","radius","centrality", "tpcClV0Pr", "tpcClV0Pi","isReconstructed","isPrimary"]
    branch_list_mc = ROOT.std.vector('std::string')()
    for b in b_list_mc:
        branch_list_mc.push_back(b)
    # split MC sample into pseudo-data + (training + testing)
    df_MC = ROOT.RDataFrame("LambdaTree","../data/LambdaPrompt_PbPb/AnalysisResults_LambdaMC_.root")
    df_index = df_MC.Define("index","gRandom->Rndm()+mass-mass")
    df_index.Filter("index < 0.5 && isReconstructed").Snapshot("LambdaTree","../data/LambdaPrompt_PbPb/pseudodataSample_small.root",branch_list)

    # get sidebands (both for training and as pseudo-data)
    df_data = ROOT.RDataFrame("LambdaTree","../data/old_Lambda_PbPb/data.root")
    df_index = df_data.Define("index","gRandom->Rndm()+mass-mass")
    df_index.Filter("index < 0.05 && (mass < 1.105 || mass > 1.13)").Snapshot("LambdaTree","../data/LambdaPrompt_PbPb/trainingBackground.root")
    df_index.Filter("index > 0.05 && index < 0.30 && (mass < 1.105 || mass > 1.13)").Snapshot("LambdaTree","../data/LambdaPrompt_PbPb/pseudodataBackground.root",branch_list)#,{"pt","ct","mass","matter","eta","radius","centrality","hasITSrefit","cosPA", "dcaV0tracks", "dcaPiPV", "dcaPrPV", "dcaV0PV", "tpcNsigmaPr", "tpcNsigmaPi", "tpcClV0Pr", "tpcClV0Pi", "radius"})

