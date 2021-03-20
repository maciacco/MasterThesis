#!/usr/bin/env python3
import os
import pickle
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import ROOT
import uproot
import xgboost as xgb
import yaml
from hipe4ml.model_handler import ModelHandler
from hipe4ml.tree_handler import TreeHandler

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
eff_array = np.arange(0.10, 0.51, 0.01)

for split in SPLIT_LIST:
    for cent_bins in CENTRALITY_LIST:
        for ct_bins in zip(CT_BINS[:-1], CT_BINS[1:]):

            bin = f'{split}_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}'
            df_data = pd.read_parquet(f'df/{bin}')

            for eff_score in zip(eff_array, score_eff_arrays_dict[bin]):
                # print(f'eff = {1-eff_score[0]}, score = {eff_score[1]}')
                df_data.query(f'model_output > {eff_score[1]}')
