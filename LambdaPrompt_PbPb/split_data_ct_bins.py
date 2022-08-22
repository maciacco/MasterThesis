import ROOT
import os
import pickle
import warnings
import argparse
import yaml
import pandas as pd
import matplotlib.pyplot as plt
import helpers
import numpy as np
import uproot
from hipe4ml import analysis_utils, plot_utils

SPLIT = True
MAX_EFF = 1

ROOT.gInterpreter.ProcessLine(".L help.h+")
warnings.simplefilter(action='ignore', category=FutureWarning)
ROOT.EnableImplicitMT(40)
ROOT.gROOT.SetBatch()

ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)
ROOT.RooMsgService.instance().setSilentMode(ROOT.kTRUE)
ROOT.gErrorIgnoreLevel = ROOT.kError

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

score_eff_dict = pickle.load(open('file_score_eff_dict','rb'))
eff_array = np.arange(0.10, MAX_EFF, 0.01)

for i_cent_bins in range(len(CENTRALITY_LIST)):

    cent_bins = CENTRALITY_LIST[i_cent_bins]
    if not os.path.isdir(f"/data/mciacco/LambdaPrompt_PbPb/df_data_{cent_bins[0]}_{cent_bins[1]}"):
        os.mkdir(f"/data/mciacco/LambdaPrompt_PbPb/df_data_{cent_bins[0]}_{cent_bins[1]}")

    df_data_ = ROOT.RDataFrame("LambdaTreeBDTOut",f"/data/mciacco/LambdaPrompt_PbPb/AnalysisResults_lambda_{cent_bins[0]}_{cent_bins[1]}.root")
    
    for ct_ in CT_BINS:
        if ct_[0] < 10 or ct_[1] > 40:
            continue
        for ct_bins in zip(CT_BINS_CENT[i_cent_bins][:-1], CT_BINS_CENT[i_cent_bins][1:]):
            if ct_bins[0] < ct_[0] or ct_bins[1] > ct_[1]:
                continue
            # if ct_bins[0] < 32 or ct_bins[1] > ct_[1]:
            #     continue
        
            df_data_.Filter(f"pt > 0.5 && pt < 3.5 && ct > {ct_bins[0]} && ct < {ct_bins[1]}").Snapshot("LambdaTreeBDTOut",f"/data/mciacco/LambdaPrompt_PbPb/df_data_{cent_bins[0]}_{cent_bins[1]}/AnalysisResults_lambda_{cent_bins[0]}_{cent_bins[1]}_ct_{ct_bins[0]}_{ct_bins[1]}.root")
