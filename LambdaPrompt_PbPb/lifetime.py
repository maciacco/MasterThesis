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
N_TRIALS = 100
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

for split in SPLIT_LIST:
    split_ineq_sign = '> -0.1'
    if SPLIT:
        split_ineq_sign = '> 0.5'
        if split == 'antimatter':
            split_ineq_sign = '< 0.5'

    for i_cent_bins in range(len(CENTRALITY_LIST)):
        if i_cent_bins > 0:
            continue
        cent_bins = CENTRALITY_LIST[i_cent_bins]

        raw_yields_file = ROOT.TFile(f'fit_{cent_bins[0]}_{cent_bins[1]}.root')
        score_eff_dict = pickle.load(open('second_round/file_score_eff_dict','rb'))
        eff_array = np.arange(0.10, MAX_EFF, 0.01)
        presel_eff_file = ROOT.TFile(f'PreselEff_{cent_bins[0]}_{cent_bins[1]}.root')
        f = ROOT.TFile(f"lifetime_{cent_bins[0]}_{cent_bins[1]}.root","recreate")
        h = ROOT.TH1D(f"h_{split}",f"h_{split}",len(CT_BINS_CENT[i_cent_bins])-1,np.asarray(CT_BINS_CENT[i_cent_bins],dtype="float"))

        h_pres_eff = presel_eff_file.Get(f"fPreselEff_vs_ct_{split}_{cent_bins[0]}_{cent_bins[1]};1")
        for ct_bins in zip(CT_BINS_CENT[i_cent_bins][:-1], CT_BINS_CENT[i_cent_bins][1:]):
            if (ct_bins[0] < 10 or ct_bins[1] > 40): # or (ct_bins[0] > 32 and ct_bins[1] < 35):
                continue
            h_raw = raw_yields_file.Get(f"{split}_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}_pol/fRawYields_pol")
            if not h_raw:
                continue
            bin = f'{split}_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}'
            print(bin)

            delta_t = ct_bins[1]-ct_bins[0]
            raw_yield = h_raw.GetBinContent(h_raw.FindBin(0.8))
            raw_yield_error = h_raw.GetBinError(h_raw.FindBin(0.8))
            j = 0
            i = np.power(-1,j)*0.01*int(j/2)
            while raw_yield < 1.e-6 and (0.8+i < 0.9) and (0.8+i > 0.7):
                i = np.power(-1,j)*0.01*int(j/2)
                j = j + 1
                raw_yield = h_raw.GetBinContent(h_raw.FindBin(0.8+i))
                raw_yield_error = h_raw.GetBinError(h_raw.FindBin(0.8+i))
            h.SetBinContent(h.FindBin(0.5*(ct_bins[0]+ct_bins[1])),raw_yield/delta_t/h_pres_eff.GetBinContent(h.FindBin(0.5*(ct_bins[0]+ct_bins[1]))))
            h.SetBinError(h.FindBin(0.5*(ct_bins[0]+ct_bins[1])),raw_yield*h_pres_eff.GetBinError(h.FindBin(0.5*(ct_bins[0]+ct_bins[1])))/delta_t/h_pres_eff.GetBinContent(h.FindBin(0.5*(ct_bins[0]+ct_bins[1])))/h_pres_eff.GetBinContent(h.FindBin(0.5*(ct_bins[0]+ct_bins[1]))))
        h.Fit("expo","IMR+","",0,40)
        h.GetXaxis().SetTitle("#it{c}t (cm)")
        h.GetYaxis().SetTitle("d#it{N}/d(#it{c}t) (cm^{-1})")
        lifetime = ROOT.TLatex(5,1.e6,"#it{c}#tau = " + "{:.2f}".format(-1./h.GetFunction("expo").GetParameter(1)/speed_of_light) + " +/- "  + "{:.2f}".format(h.GetFunction("expo").GetParError(1)/h.GetFunction("expo").GetParameter(1)/h.GetFunction("expo").GetParameter(1)/speed_of_light) + " ps")
        f.cd()
        c = ROOT.TCanvas(f"lifetime_{split}",f"lifetime_{split}")
        h.Draw()
        lifetime.Draw("same")
        h.Write()
        c.SetLogy()
        c.Write()
        c.Print(f"plots/lifetime_{split}_{cent_bins[0]}_{cent_bins[1]}.png")

        # systematics
        # h_sys = ROOT.TH1D(f"h_sys_{split}",f"h_sys_{split}",10000,0,1000)
        # i=0
        # while i < N_TRIALS:
        #     h_tmp = ROOT.TH1D("h_tmp","h_tmp",len(CT_BINS_CENT[i_cent_bins])-1,np.asarray(CT_BINS_CENT[i_cent_bins],dtype="float"))
        #     cent_bins = CENTRALITY_LIST[i_cent_bins]
        #     h_pres_eff = presel_eff_file.Get(f"fPreselEff_vs_ct_{split}_{cent_bins[0]}_{cent_bins[1]}")
        #     for ct_bins in zip(CT_BINS_CENT[i_cent_bins][:-1], CT_BINS_CENT[i_cent_bins][1:]):
        #         if ct_bins[0] < 7 or ct_bins[1] > 40:
        #             continue
        #         cut_rndm = ROOT.gRandom.Rndm()*0.1+0.85
        #         bkg_index = int(ROOT.gRandom.Rndm())
        #         h_raw = raw_yields_file.Get(f"{split}_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}_{bkg_function[bkg_index]}/fRawYields_{bkg_function[bkg_index]}")
        #         if not h_raw:
        #             continue
        #         bin = f'{split}_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}'
        #         print(bin)

        #         h_tmp.SetBinContent(h_tmp.FindBin(0.5*(ct_bins[0]+ct_bins[1])),h_raw.GetBinContent(h_raw.FindBin(cut_rndm+0.001))/h_pres_eff.GetBinContent(h_tmp.FindBin(0.5*(ct_bins[0]+ct_bins[1]))))
        #         h_tmp.SetBinError(h_tmp.FindBin(0.5*(ct_bins[0]+ct_bins[1])),h_raw.GetBinError(h_raw.FindBin(cut_rndm+0.001))/h_pres_eff.GetBinContent(h_tmp.FindBin(0.5*(ct_bins[0]+ct_bins[1]))))
            
        #     h_tmp.Fit("expo")
        #     i = i+1
        #     if (h_tmp.GetFunction("expo").GetChisquare()/h_tmp.GetFunction("expo").GetNDF()>300.):
        #         continue
        #     slope = h_tmp.GetFunction("expo").GetParameter(1)
        #     h_sys.Fill(-1./slope/speed_of_light)
        # h_sys.Write()
