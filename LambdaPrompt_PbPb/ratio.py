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

SPLIT = True
N_TRIALS = 1000
MAX_EFF = 1
speed_of_light = 0.0299792458

warnings.simplefilter(action='ignore', category=FutureWarning)
ROOT.gROOT.SetBatch()

bkg_function = ['pol','expo']

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

raw_yields_file = ROOT.TFile('fit_13_noRadius_nodcaV0tracksCut_lambda_splitCentInMC_saveHistos_0_5.root')
# score_eff_dict = pickle.load(open('second_round/file_score_eff_dict','rb'))
eff_array = np.arange(0.10, MAX_EFF, 0.01)
presel_eff_file = ROOT.TFile('PreselEff.root')
f = ROOT.TFile("ratio_lambda_splitCentInMC__0_5.root","recreate")

# for i_cent_bins in range(len(CENTRALITY_LIST)):
#     h = []
#     h.append(ROOT.TH1D(f"h_antimatter",f"h_antimatter",len(CT_BINS_CENT[i_cent_bins])-1,np.asarray(CT_BINS_CENT[i_cent_bins],dtype="float")))
#     h.append(ROOT.TH1D(f"h_matter",f"h_matter",len(CT_BINS_CENT[i_cent_bins])-1,np.asarray(CT_BINS_CENT[i_cent_bins],dtype="float")))

#     for i_split, split in enumerate(SPLIT_LIST):
#         split_ineq_sign = '> -0.1'
#         if SPLIT:
#             split_ineq_sign = '> 0.5'
#             if split == 'antimatter':
#                 split_ineq_sign = '< 0.5'

#         cent_bins = CENTRALITY_LIST[i_cent_bins]
#         h_pres_eff = presel_eff_file.Get(f"fPreselEff_vs_ct_{split}_{cent_bins[0]}_{cent_bins[1]};2")
#         for ct_bins in zip(CT_BINS_CENT[i_cent_bins][:-1], CT_BINS_CENT[i_cent_bins][1:]):
#             if ct_bins[0] < 5 or ct_bins[1] > 40:
#                 continue
#             h_raw = raw_yields_file.Get(f"{split}_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}_pol/fRawYields_pol")
#             if not h_raw:
#                 continue
#             bin = f'{split}_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}'
#             print(bin)

#             delta_t = ct_bins[1]-ct_bins[0]
#             raw_yield = h_raw.GetBinContent(h_raw.FindBin(0.9))
#             raw_yield_error = h_raw.GetBinError(h_raw.FindBin(0.9))
#             j = 0
#             i = np.power(-1,j)*0.01*int(j/2)
#             while raw_yield < 1.e-6 and (0.9+i < 0.9) and (0.9+i > 0.7):
#                 i = np.power(-1,j)*0.01*int(j/2)
#                 j = j + 1
#                 raw_yield = h_raw.GetBinContent(h_raw.FindBin(0.9+i))
#                 raw_yield_error = h_raw.GetBinError(h_raw.FindBin(0.9+i))
#             h[i_split].SetBinContent(h[i_split].FindBin(0.5*(ct_bins[0]+ct_bins[1])),raw_yield/delta_t/h_pres_eff.GetBinContent(h[i_split].FindBin(0.5*(ct_bins[0]+ct_bins[1]))))
#             h[i_split].SetBinError(h[i_split].FindBin(0.5*(ct_bins[0]+ct_bins[1])),raw_yield_error/delta_t/h_pres_eff.GetBinContent(h[i_split].FindBin(0.5*(ct_bins[0]+ct_bins[1]))))

#     h_ratio = ROOT.TH1D(f"h_ratio_{cent_bins[0]}_{cent_bins[1]}",f"h_ratio_{cent_bins[0]}_{cent_bins[1]}",len(CT_BINS_CENT[i_cent_bins])-1,np.asarray(CT_BINS_CENT[i_cent_bins],dtype="float"))
#     h_ratio.Divide(h[0],h[1],1,1)
#     h_ratio.GetXaxis().SetTitle("#it{c}t (cm)")
#     h_ratio.GetYaxis().SetTitle("Ratio #bar{#Lambda}/#Lambda")
#     h_ratio.Fit("pol0")
#     h_ratio.Write()

for i_cent_bins in range(len(CENTRALITY_LIST)):
    h = []
    h.append(ROOT.TH1D(f"h_antimatter",f"h_antimatter",len(CT_BINS_CENT[i_cent_bins])-1,np.asarray(CT_BINS_CENT[i_cent_bins],dtype="float")))
    h.append(ROOT.TH1D(f"h_matter",f"h_matter",len(CT_BINS_CENT[i_cent_bins])-1,np.asarray(CT_BINS_CENT[i_cent_bins],dtype="float")))

    cent_bins = CENTRALITY_LIST[i_cent_bins]
    h_pres_eff_a = presel_eff_file.Get(f"fPreselEff_vs_ct_antimatter_{cent_bins[0]}_{cent_bins[1]};2")
    h_pres_eff_m = presel_eff_file.Get(f"fPreselEff_vs_ct_matter_{cent_bins[0]}_{cent_bins[1]};2")
    for ct_bins in zip(CT_BINS_CENT[i_cent_bins][:-1], CT_BINS_CENT[i_cent_bins][1:]):
        if ct_bins[0] < 5 or ct_bins[1] > 40:
            continue
        h_raw_a = raw_yields_file.Get(f"antimatter_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}_expo/fRawYields_expo")
        h_raw_m = raw_yields_file.Get(f"matter_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}_expo/fRawYields_expo")
        if not h_raw_a or not h_raw_m:
            continue
        bin = f'all_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}'
        print(bin)

        delta_t = ct_bins[1]-ct_bins[0]
        raw_yield_a = h_raw_a.GetBinContent(h_raw_a.FindBin(0.8))
        raw_yield_m = h_raw_m.GetBinContent(h_raw_m.FindBin(0.8))
        raw_yield_error_a = h_raw_a.GetBinError(h_raw_a.FindBin(0.8))
        raw_yield_error_m = h_raw_m.GetBinError(h_raw_m.FindBin(0.8))
        j = 1
        i = np.power(-1,j)*0.01*int(j/2)
        while (raw_yield_a < 1.e-6 or raw_yield_m < 1.e-6) and (0.8+i < 0.8) and (0.8+i > 0.7):
            i = np.power(-1,j)*0.01*int(j/2)
            j = j + 1
            raw_yield_a = h_raw_a.GetBinContent(h_raw_a.FindBin(0.8+i))
            raw_yield_m = h_raw_m.GetBinContent(h_raw_m.FindBin(0.8+i))
            raw_yield_error_a = h_raw_a.GetBinError(h_raw_a.FindBin(0.8+i))
            raw_yield_error_m = h_raw_m.GetBinError(h_raw_m.FindBin(0.8+i))
        
        eff_a = h_pres_eff_a.GetBinContent(h[0].FindBin(0.5*(ct_bins[0]+ct_bins[1])))
        eff_error_a = h_pres_eff_a.GetBinError(h[0].FindBin(0.5*(ct_bins[0]+ct_bins[1])))
        eff_m = h_pres_eff_m.GetBinContent(h[0].FindBin(0.5*(ct_bins[0]+ct_bins[1])))
        eff_error_m = h_pres_eff_m.GetBinError(h[0].FindBin(0.5*(ct_bins[0]+ct_bins[1])))
        h[0].SetBinContent(h[0].FindBin(0.5*(ct_bins[0]+ct_bins[1])),raw_yield_a/delta_t/h_pres_eff_a.GetBinContent(h[0].FindBin(0.5*(ct_bins[0]+ct_bins[1]))))
        h[0].SetBinError(h[0].FindBin(0.5*(ct_bins[0]+ct_bins[1])),np.sqrt((raw_yield_error_a**2)/(eff_a**2)+(eff_error_a**2)/(eff_a**4)*(raw_yield_a**2)))
        h[1].SetBinContent(h[1].FindBin(0.5*(ct_bins[0]+ct_bins[1])),raw_yield_m/delta_t/h_pres_eff_m.GetBinContent(h[1].FindBin(0.5*(ct_bins[0]+ct_bins[1]))))
        h[1].SetBinError(h[1].FindBin(0.5*(ct_bins[0]+ct_bins[1])),np.sqrt((raw_yield_error_m**2)/(eff_m**2)+(eff_error_m**2)/(eff_m**4)*(raw_yield_m**2)))

    h_ratio = ROOT.TH1D(f"h_ratio",f"h_ratio",len(CT_BINS_CENT[i_cent_bins])-1,np.asarray(CT_BINS_CENT[i_cent_bins],dtype="float"))
    h_ratio.Divide(h[0],h[1],1,1)
    h_ratio.Fit("pol0")
    h_ratio.Write()


    # systematics
    # h_sys = ROOT.TH1D(f"h_sys",f"h_sys",100,0.85,1)
    # i=0
    # while i < N_TRIALS:
    #     h_tmp = []
    #     h_tmp.append(ROOT.TH1D("h_tmp_a","h_tmp_a",len(CT_BINS_CENT[i_cent_bins])-1,np.asarray(CT_BINS_CENT[i_cent_bins],dtype="float")))
    #     h_tmp.append(ROOT.TH1D("h_tmp_m","h_tmp_m",len(CT_BINS_CENT[i_cent_bins])-1,np.asarray(CT_BINS_CENT[i_cent_bins],dtype="float")))
    #     for ct_bins in zip(CT_BINS_CENT[i_cent_bins][:-1], CT_BINS_CENT[i_cent_bins][1:]):
    #         if ct_bins[0] < 5 or ct_bins[1] > 40:
    #             continue
    #         cut_rndm = int(ROOT.gRandom.Rndm()*20)*0.01+0.7
    #         bkg_index = int(ROOT.gRandom.Rndm())
    #         split_ineq_sign = '> -0.1'
    #         for i_split, split in enumerate(SPLIT_LIST):
    #             if SPLIT:
    #                 split_ineq_sign = '> 0.5'
    #                 if split == 'antimatter':
    #                     split_ineq_sign = '< 0.5'
    #             cent_bins = CENTRALITY_LIST[i_cent_bins]
    #             h_pres_eff = presel_eff_file.Get(f"fPreselEff_vs_ct_{split}_{cent_bins[0]}_{cent_bins[1]};2")
    #             h_raw = raw_yields_file.Get(f"{split}_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}_{bkg_function[bkg_index]}/fRawYields_{bkg_function[bkg_index]}")
    #             if not h_raw:
    #                 continue
    #             bin = f'{split}_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}_{cut_rndm}'
    #             print(bin)

    #             h_tmp[i_split].SetBinContent(h_tmp[i_split].FindBin(0.5*(ct_bins[0]+ct_bins[1])),h_raw.GetBinContent(h_raw.FindBin(cut_rndm+0.0001))/h_pres_eff.GetBinContent(h_tmp[i_split].FindBin(0.5*(ct_bins[0]+ct_bins[1]))))
    #             h_tmp[i_split].SetBinError(h_tmp[i_split].FindBin(0.5*(ct_bins[0]+ct_bins[1])),h_raw.GetBinError(h_raw.FindBin(cut_rndm+0.0001))/h_pres_eff.GetBinContent(h_tmp[i_split].FindBin(0.5*(ct_bins[0]+ct_bins[1]))))

    #     h_tmp_ratio = ROOT.TH1D("h_tmp_ratio","h_tmp_ratio",len(CT_BINS_CENT[i_cent_bins])-1,np.asarray(CT_BINS_CENT[i_cent_bins],dtype="float"))
    #     h_tmp_ratio.Divide(h_tmp[0],h_tmp[1])
    #     h_tmp_ratio.Fit("pol0")
    #     if (h_tmp_ratio.GetFunction("pol0").GetChisquare()/h_tmp_ratio.GetFunction("pol0").GetNDF()>2.):
    #         continue
    #     i = i+1
    #     h_sys.Fill(h_tmp_ratio.GetFunction("pol0").GetParameter(0))
    # h_sys.Write()