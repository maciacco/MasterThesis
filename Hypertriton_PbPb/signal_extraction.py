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


def ndarray2roo(ndarray, var):
    if isinstance(ndarray, ROOT.RooDataSet):
        print('Already a RooDataSet')
        return ndarray

    assert isinstance(ndarray, np.ndarray), 'Did not receive NumPy array'
    assert len(ndarray.shape) == 1, 'Can only handle 1d array'

    name = var.GetName()
    x = np.zeros(1, dtype=np.float64)

    tree = ROOT.TTree('tree', 'tree')
    tree.Branch(f'{name}', x, f'{name}/D')

    for i in ndarray:
        x[0] = i
        tree.Fill()

    array_roo = ROOT.RooDataSet(
        'data', 'dataset from tree', tree, ROOT.RooArgSet(var))
    return array_roo


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
                formatted_eff="{:.2f}".format(1-eff_score[0])
                print(f'eff = {1-eff_score[0]:.2f}, score = {eff_score[1]:.2f}')
                df_data_sel = df_data.query(f'model_output > {eff_score[1]}')

                roo_m = ROOT.RooRealVar("m", "#it{M} (^{3}He + #pi^{-})", 2.96, 3.04, "GeV/#it{c}^{2}")
                roo_data = ndarray2roo(np.array(df_data_sel['m']), roo_m)
                frame = roo_m.frame(2.96, 3.04)
                frame.SetTitle(str(ct_bins[0])+'#leq #it{c}t<'+str(ct_bins[1])+' cm, '+str(cent_bins[0])+'-'+str(cent_bins[1])+'%, '+str(formatted_eff))
                roo_data.plotOn(frame)
                canv = ROOT.TCanvas()
                frame.Draw("")

                if not os.path.isdir('plots/signal_extraction'):
                    os.mkdir('plots/signal_extraction')
                canv.Print(f'plots/signal_extraction/{1-eff_score[0]:.2f}_{bin}.png')
