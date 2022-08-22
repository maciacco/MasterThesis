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

ROOT.gInterpreter.ProcessLine("#include \"../utils/RooDSCBShape.h\"")
ROOT.gInterpreter.ProcessLine(".L ../utils/RooDSCBShape.cxx+")
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

if not os.path.isdir("plots/signal_extraction_"):
    os.mkdir("plots/signal_extraction_")

#df_data = uproot.open(os.path.expandvars(f"/data/mciacco/LambdaPrompt_PbPb/AnalysisResults_data_filtered_dcaV0PV.root"))['LambdaTree'].arrays(library="pd")

for i_cent_bins in range(len(CENTRALITY_LIST)):

    cent_bins = CENTRALITY_LIST[i_cent_bins]
    if not os.path.isdir(f"plots/signal_extraction_/{cent_bins[0]}_{cent_bins[1]}"):
        os.mkdir(f"plots/signal_extraction_/{cent_bins[0]}_{cent_bins[1]}")
    df_data_gen = uproot.open(os.path.expandvars(f"/data/mciacco/LambdaPrompt_PbPb/AnalysisResults_reweight_BW_{cent_bins[0]}_{cent_bins[1]}.root"))['LambdaTree'].arrays(library="pd")
    for ct_ in CT_BINS:
        if ct_[0] < 5 or ct_[1] > 40:
            continue
        for ct_bins in zip(CT_BINS_CENT[i_cent_bins][:-1], CT_BINS_CENT[i_cent_bins][1:]):
            if ct_bins[0] < ct_[0] or ct_bins[1] > ct_[1]:
                continue
        
            h_eff = [[],[]] # 0 -> antim; 1 -> m
            h_eff[0] = [ROOT.TH1D("fAccEffP","fAccEffP",100,0,1),ROOT.TH1D("fAccEffNP","fAccEffNP",100,0,1)]

            h_eff[1] = [ROOT.TH1D("fAccEffP","fAccEffP",100,0,1),ROOT.TH1D("fAccEffNP","fAccEffNP",100,0,1)]

            bin = f'all_{cent_bins[0]}_{cent_bins[1]}_{ct_[0]}_{ct_[1]}_mc_apply'
            df_data = pd.read_parquet(f'df_test/{bin}.parquet.gzip')
            df_data_ = df_data.query(f"mass > 1.09 and mass < 1.15 and ct >= {ct_bins[0]} and ct < {ct_bins[1]}")
            df_data_gen_ = df_data_gen.query(f"ctMC >= {ct_bins[0]} and ctMC < {ct_bins[1]} and ptMC >= 0.5 and ptMC < 3.5")

            for bdt_score in np.arange(0,1,0.05):
                print(f'processing {bin}: bkg cut = {bdt_score:.4f}')

                for i_split, split in enumerate(SPLIT_LIST):
                    split_ineq_sign = '> -0.1'
                    if SPLIT:
                        split_ineq_sign = '> 0.5'
                        if split == 'antimatter':
                            split_ineq_sign = '< 0.5'

                    for i_p_np, p_np in enumerate([1,2]):
                        # get invariant mass
                        df_data_cut = df_data_.query(f"flag == {p_np} and matter {split_ineq_sign} and model_output_0 < .2 and model_output_1 > {bdt_score} and pt < 3.5")
                        df_data_gen_cut = df_data_gen_.query(f"pdg {split_ineq_sign} and flag == {p_np}")
                        rec = df_data_cut.shape[0]
                        gen = df_data_gen_cut.shape[0]
                        h_eff[i_split][i_p_np].SetBinContent(h_eff[i_split][i_p_np].FindBin(bdt_score+0.0001),rec/gen)
                        h_eff[i_split][i_p_np].SetBinError(h_eff[i_split][i_p_np].FindBin(bdt_score+0.0001),np.sqrt((rec/gen)*(1-rec/gen)/gen))

            f = ROOT.TFile(f"efficiency_noRadius_nodcaV0tracksCut_lambda_splitCentInMC_saveHistos_{cent_bins[0]}_{cent_bins[1]}.root","update")
            for i_split, split in enumerate(SPLIT_LIST):
                split_ineq_sign = '> -0.1'
                if SPLIT:
                    split_ineq_sign = '> 0.5'
                    if split == 'antimatter':
                        split_ineq_sign = '< 0.5'
                for i_p_np, p_np in enumerate([1,2]):
                    if not f.cd(f'{split}_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}_{i_p_np}'):
                        f.mkdir(f'{split}_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}_{i_p_np}')
                    f.cd(f'{split}_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}_{i_p_np}')
                    h_eff[i_split][i_p_np].Write()
            f.Close()
