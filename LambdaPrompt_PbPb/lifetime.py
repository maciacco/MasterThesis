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

SPLIT = False

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

raw_yields_file = ROOT.TFile('out.root')
score_eff_dict = pickle.load(open('file_score_eff_dict','rb'))
presel_eff_file = ROOT.TFile('PreselEff.root')
f = ROOT.TFile("f.root","recreate")

for split in SPLIT_LIST:
    split_ineq_sign = '> -0.1'
    if SPLIT:
        split_ineq_sign = '> 0.5'
        if split == 'antimatter':
            split_ineq_sign = '< 0.5'

    for i_cent_bins in range(len(CENTRALITY_LIST)):
        h = ROOT.TH1D("h","h",80,0,40)
        h_raw = raw_yields_file.Get("h")
        h_pres_eff = presel_eff_file.Get("fPreselEff_vs_ct_all_0_90;2")
        cent_bins = CENTRALITY_LIST[i_cent_bins]
        for ct_bins in zip(CT_BINS_CENT[i_cent_bins][:-1], CT_BINS_CENT[i_cent_bins][1:]):
            if ct_bins[0] < 5 or ct_bins[1] > 10:
                continue
            bdt_cut = 0.2
            bin = f'{split}_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}'
            h.SetBinContent(h.FindBin(0.5*(ct_bins[0]+ct_bins[1])),h_raw.GetBinContent(h_raw.FindBin(0.5*(ct_bins[0]+ct_bins[1])))/h_pres_eff.GetBinContent(h.FindBin(0.5*(ct_bins[0]+ct_bins[1])))/0.9)
        f.cd()
        h.Write()