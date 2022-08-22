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
           
for i_cent_bins in range(len(CENTRALITY_LIST)):

    cent_bins = CENTRALITY_LIST[i_cent_bins]

    file_eff = ROOT.TFile.Open(f"efficiency_noRadius_nodcaV0tracksCut_lambda_splitCentInMC_saveHistos_{cent_bins[0]}_{cent_bins[1]}.root")
    file_raw = ROOT.TFile.Open(f"fit_yield_noRadius_nodcaV0tracksCut_lambda_splitCentInMC_saveHistos_{cent_bins[0]}_{cent_bins[1]}.root")

    if not os.path.isdir(f"plots/signal_extraction_/{cent_bins[0]}_{cent_bins[1]}"):
        os.mkdir(f"plots/signal_extraction_/{cent_bins[0]}_{cent_bins[1]}")
    for ct_ in CT_BINS:
        if ct_[0] < 5 or ct_[1] > 15:
            continue
        for ct_bins in zip(CT_BINS_CENT[i_cent_bins][:-1], CT_BINS_CENT[i_cent_bins][1:]):
            if ct_bins[0] < ct_[0] or ct_bins[1] > ct_[1]:
                continue
            
            h_corr_yields = []
            h_yields_p = []
            h_yields_np = []
            h_yields = []
            h_corr_yields.append(ROOT.TH1D(f"fACorrYields_{cent_bins[0]}_{cent_bins[1]}",f"fACorrYields_{cent_bins[0]}_{cent_bins[1]}",100,0,1))
            h_corr_yields.append(ROOT.TH1D(f"fMCorrYields_{cent_bins[0]}_{cent_bins[1]}",f"fMCorrYields_{cent_bins[0]}_{cent_bins[1]}",100,0,1))
            h_yields_p.append(ROOT.TH1D(f"fAYieldsP_{cent_bins[0]}_{cent_bins[1]}",f"fAYieldsP_{cent_bins[0]}_{cent_bins[1]}",100,0,1))
            h_yields_p.append(ROOT.TH1D(f"fMYieldsP_{cent_bins[0]}_{cent_bins[1]}",f"fMYieldsP_{cent_bins[0]}_{cent_bins[1]}",100,0,1))
            h_yields_np.append(ROOT.TH1D(f"fAYieldsNP_{cent_bins[0]}_{cent_bins[1]}",f"fAYieldsNP_{cent_bins[0]}_{cent_bins[1]}",100,0,1))
            h_yields_np.append(ROOT.TH1D(f"fMYieldsNP_{cent_bins[0]}_{cent_bins[1]}",f"fMYieldsNP_{cent_bins[0]}_{cent_bins[1]}",100,0,1))
            h_yields.append(ROOT.TH1D(f"fAYields_{cent_bins[0]}_{cent_bins[1]}",f"fAYields_{cent_bins[0]}_{cent_bins[1]}",100,0,1))
            h_yields.append(ROOT.TH1D(f"fMYields_{cent_bins[0]}_{cent_bins[1]}",f"fMYields_{cent_bins[0]}_{cent_bins[1]}",100,0,1))

            for i_split, split in enumerate(SPLIT_LIST):
                split_ineq_sign = '> -0.1'
                if SPLIT:
                    split_ineq_sign = '> 0.5'
                    if split == 'antimatter':
                        split_ineq_sign = '< 0.5'

                h_raw = file_raw.Get(f"{split}_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}_pol/fRawYields_pol")
                h_eff_p = file_eff.Get(f"{split}_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}_0/fAccEffP")
                h_eff_np = file_eff.Get(f"{split}_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}_1/fAccEffNP")
                l_raw = []
                l_raw_err = []
                l_eff_p = []
                l_eff_np = []
                l_eff_p_err = []
                l_eff_np_err = []
                l_score = []

                bin = f'all_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}'
                for bdt_score in np.arange(0,1,0.05):
                    print(f'processing {bin}: bkg cut = {bdt_score:.4f}')
                    tmp_raw = h_raw.GetBinContent(h_raw.FindBin(bdt_score+0.0001))
                    tmp_raw_err = h_raw.GetBinError(h_raw.FindBin(bdt_score+0.0001))
                    tmp_eff_p = h_eff_p.GetBinContent(h_eff_p.FindBin(bdt_score+0.0001))
                    tmp_eff_p_err = h_eff_p.GetBinError(h_eff_p.FindBin(bdt_score+0.0001))
                    tmp_eff_np = h_eff_np.GetBinContent(h_eff_np.FindBin(bdt_score+0.0001))
                    tmp_eff_np_err = h_eff_np.GetBinError(h_eff_np.FindBin(bdt_score+0.0001))
                    if tmp_raw < 1.e-6 or tmp_eff_p < 1.e-6 or tmp_eff_np < 1.e-6:
                        continue
                    l_raw.append(tmp_raw)
                    l_raw_err.append(tmp_raw_err)
                    l_eff_p.append(tmp_eff_p)
                    l_eff_p_err.append(tmp_eff_p_err)
                    l_eff_np.append(tmp_eff_np)
                    l_eff_np_err.append(tmp_eff_np_err)
                    l_score.append(bdt_score)
                corr_yield, cov, red_chi2, dict_of_mat = helpers.GetPromptFDYieldsAnalyticMinimisation(l_eff_p,l_eff_np,l_raw,l_eff_p_err,l_eff_np_err,l_raw_err,nMaxIter=10000)
                print(f"chi2 = {red_chi2}")
                for (raw, eff_p, eff_np, raw_err, eff_p_err, eff_np_err, bdt_score) in zip(l_raw, l_eff_p, l_eff_np, l_raw_err, l_eff_p_err, l_eff_np_err, l_score):
                    f_p = eff_p*corr_yield[0]/(eff_p*corr_yield[0]+eff_np*corr_yield[1])
                    corr_yield_p = f_p*raw/eff_p
                    corr_yield_p_err = corr_yield_p*np.sqrt(raw_err*raw_err/raw/raw+eff_p_err*eff_p_err/eff_p/eff_p)
                    h_corr_yields[i_split].SetBinContent(h_corr_yields[i_split].FindBin(bdt_score+0.0001),corr_yield_p)
                    h_corr_yields[i_split].SetBinError(h_corr_yields[i_split].FindBin(bdt_score+0.0001),corr_yield_p_err)
                    h_yields_p[i_split].SetBinContent(h_corr_yields[i_split].FindBin(bdt_score+0.0001),corr_yield[0]*eff_p)
                    h_yields_np[i_split].SetBinContent(h_corr_yields[i_split].FindBin(bdt_score+0.0001),corr_yield[1]*eff_np)
                    h_yields[i_split].SetBinContent(h_corr_yields[i_split].FindBin(bdt_score+0.0001),corr_yield[1]*eff_np+corr_yield[0]*eff_p)

            f = ROOT.TFile(f"corrYields_noRadius_nodcaV0tracksCut_lambda_splitCentInMC_saveHistos_{cent_bins[0]}_{cent_bins[1]}.root","update")
            for i_split, split in enumerate(SPLIT_LIST):
                split_ineq_sign = '> -0.1'
                if SPLIT:
                    split_ineq_sign = '> 0.5'
                    if split == 'antimatter':
                        split_ineq_sign = '< 0.5'
                if not f.cd(f'{split}_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}'):
                    f.mkdir(f'{split}_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}')
                f.cd(f'{split}_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}')
                h_corr_yields[i_split].Write()
                h_yields_p[i_split].Write()
                h_yields_np[i_split].Write()
                h_yields[i_split].Write()
            f.Close()
