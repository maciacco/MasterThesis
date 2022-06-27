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
N_TRIALS = 10000
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

raw_yields_file = ROOT.TFile('out_1.root')
score_eff_dict = pickle.load(open('second_round/file_score_eff_dict','rb'))
eff_array = np.arange(0.10, MAX_EFF, 0.01)
presel_eff_file = ROOT.TFile('PreselEff.root')
f = ROOT.TFile("f.root","recreate")

for split in SPLIT_LIST:
    split_ineq_sign = '> -0.1'
    if SPLIT:
        split_ineq_sign = '> 0.5'
        if split == 'antimatter':
            split_ineq_sign = '< 0.5'

    for i_cent_bins in range(len(CENTRALITY_LIST)):
        h = ROOT.TH1D("h","h",len(CT_BINS_CENT[i_cent_bins])-1,np.asarray(CT_BINS_CENT[i_cent_bins],dtype="float"))
        h_pres_eff = presel_eff_file.Get("fPreselEff_vs_ct_all_0_90")
        cent_bins = CENTRALITY_LIST[i_cent_bins]
        for ct_bins in zip(CT_BINS_CENT[i_cent_bins][:-1], CT_BINS_CENT[i_cent_bins][1:]):
            if ct_bins[0] < 10 or ct_bins[1] > 40:
                continue
            h_raw = raw_yields_file.Get(f"all_0_0_{ct_bins[0]}_{ct_bins[1]}_pol/fRawYields")
            if not h_raw:
                continue
            bin = f'{split}_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}'
            print(bin)

            delta_t = ct_bins[1]-ct_bins[0]
            raw_yield = h_raw.GetBinContent(h_raw.FindBin(0.861))
            raw_yield_error = h_raw.GetBinError(h_raw.FindBin(0.861))
            i = 0.
            while raw_yield < 1.e-6 and (0.861+i < 0.95):
                i = i + 0.01
                raw_yield = h_raw.GetBinContent(h_raw.FindBin(0.861+i))
                raw_yield_error = h_raw.GetBinError(h_raw.FindBin(0.861+i))
            while raw_yield < 1.e-6 and (0.861+i > 0.76):
                i = i - 0.01
                raw_yield = h_raw.GetBinContent(h_raw.FindBin(0.861+i))
                raw_yield_error = h_raw.GetBinError(h_raw.FindBin(0.861+i))
            h.SetBinContent(h.FindBin(0.5*(ct_bins[0]+ct_bins[1])),raw_yield/delta_t/h_pres_eff.GetBinContent(h.FindBin(0.5*(ct_bins[0]+ct_bins[1]))))
            h.SetBinError(h.FindBin(0.5*(ct_bins[0]+ct_bins[1])),raw_yield_error/delta_t/h_pres_eff.GetBinContent(h.FindBin(0.5*(ct_bins[0]+ct_bins[1]))))
        h.Fit("expo","","I")
        lifetime = ROOT.TLatex(10,1000,"#it{c}#tau = " + str(-1./h.GetFunction("expo").GetParameter(1)/speed_of_light) + " +/- "  + str(h.GetFunction("expo").GetParError(1)/h.GetFunction("expo").GetParameter(1)/h.GetFunction("expo").GetParameter(1)/speed_of_light) + " ps")
        f.cd()
        c = ROOT.TCanvas("lifetime","lifetime")
        h.Draw()
        lifetime.Draw("same")
        h.Write()
        c.Write()

        # systematics
        h_sys = ROOT.TH1D("h_sys","h_sys",10000,0,1000)
        i=0
        while i < N_TRIALS:
            h_tmp = ROOT.TH1D("h_tmp","h_tmp",80,0,40)
            h_pres_eff = presel_eff_file.Get("fPreselEff_vs_ct_all_0_90")
            cent_bins = CENTRALITY_LIST[i_cent_bins]
            for ct_bins in zip(CT_BINS_CENT[i_cent_bins][:-1], CT_BINS_CENT[i_cent_bins][1:]):
                if ct_bins[0] < 10 or ct_bins[1] > 40:
                    continue
                cut_rndm = int(ROOT.gRandom.Rndm()*5)*0.04+0.76
                bkg_index = int(ROOT.gRandom.Rndm())
                h_raw = raw_yields_file.Get(f"all_0_0_{ct_bins[0]}_{ct_bins[1]}_{bkg_function[bkg_index]}/fRawYields")
                if not h_raw:
                    continue
                bin = f'{split}_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}'
                print(bin)

                h_tmp.SetBinContent(h_tmp.FindBin(0.5*(ct_bins[0]+ct_bins[1])),h_raw.GetBinContent(h_raw.FindBin(cut_rndm+0.001))/h_pres_eff.GetBinContent(h_tmp.FindBin(0.5*(ct_bins[0]+ct_bins[1]))))
                h_tmp.SetBinError(h_tmp.FindBin(0.5*(ct_bins[0]+ct_bins[1])),h_raw.GetBinError(h_raw.FindBin(cut_rndm+0.001))/h_pres_eff.GetBinContent(h_tmp.FindBin(0.5*(ct_bins[0]+ct_bins[1]))))
            
            h_tmp.Fit("expo")
            i = i+1
            if (h_tmp.GetFunction("expo").GetChisquare()/h_tmp.GetFunction("expo").GetNDF()>2.) or h_tmp.GetFunction("expo").GetNDF()<10:
                continue
            slope = h_tmp.GetFunction("expo").GetParameter(1)
            h_sys.Fill(-1./slope/speed_of_light)
        h_sys.Write()