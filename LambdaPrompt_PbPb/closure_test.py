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
N_TRIALS = 10000
MAX_EFF = 1
speed_of_light = 0.0299792458

warnings.simplefilter(action='ignore', category=FutureWarning)
ROOT.EnableImplicitMT(40)
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

raw_yields_file = ROOT.TFile('out_1_mc_split.root')
score_eff_dict = pickle.load(open('second_round/file_score_eff_dict','rb'))
eff_array = np.arange(0.10, MAX_EFF, 0.01)
presel_eff_file = ROOT.TFile('PreselEff.root')
f = ROOT.TFile("f_split.root","recreate")

df_MC = ROOT.RDataFrame("LambdaTree","/data/mciacco/LambdaPrompt_PbPb/mc.root")

for split in SPLIT_LIST:
    split_ineq_sign = '> -999999999.'
    if SPLIT:
        split_ineq_sign = '> 0.5'
        if split == 'antimatter':
            split_ineq_sign = '< 0.5'

    for i_cent_bins in range(len(CENTRALITY_LIST)):
        h = ROOT.TH1D(f"hYield_{split}",f"hYield_{split}",len(CT_BINS_CENT[i_cent_bins])-1,np.asarray(CT_BINS_CENT[i_cent_bins],dtype="float"))
        h_pres_eff = presel_eff_file.Get(f"fPreselEff_vs_ct_{split}_0_90")
        cent_bins = CENTRALITY_LIST[i_cent_bins]
        for ct_bins in zip(CT_BINS_CENT[i_cent_bins][:-1], CT_BINS_CENT[i_cent_bins][1:]):
            if ct_bins[0] < 1 or ct_bins[1] > 40:
                continue
            h_raw = raw_yields_file.Get(f"{split}_0_0_{ct_bins[0]}_{ct_bins[1]}_bkg/fRawYields")
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
        # h.Fit("expo","","I")
        # lifetime = ROOT.TLatex(10,1000,"#it{c}#tau = " + str(-1./h.GetFunction("expo").GetParameter(1)/speed_of_light) + " +/- "  + str(h.GetFunction("expo").GetParError(1)/h.GetFunction("expo").GetParameter(1)/h.GetFunction("expo").GetParameter(1)/speed_of_light) + " ps")
        bins = np.asarray(CT_BINS_CENT[i_cent_bins],dtype="float")
        df_MC_gen = df_MC.Filter(f"pdg {split_ineq_sign} && centrality > {cent_bins[0]} && centrality < {cent_bins[1]} && ptMC > 0.5 && ptMC < 3.5 && index_1 < 0.3 && flag == 1")
        h_MC_gen_tmp = df_MC_gen.Histo1D((f"h_{split}",f"h_{split}",len(CT_BINS_CENT[i_cent_bins])-1,1,40),"ctMC")
        h_MC_gen_tmp_pt = df_MC_gen.Histo1D((f"h_{split}_pt",f"h_{split}_pt",1000,0,10),"ptMC")
        delta_t = CT_BINS_CENT[i_cent_bins][1]-CT_BINS_CENT[i_cent_bins][0]
        h_MC_gen = h_MC_gen_tmp.GetPtr()
        h_MC_gen.Scale(1./delta_t)
        f.cd()
        h_MC_gen_tmp_pt.Write()
        c = ROOT.TCanvas(f"fYield_{split}",f"yield {split}")
        c.SetLogy()
        h_MC_gen.Draw()
        h_MC_gen.GetXaxis().SetTitle("#it{c}t (cm)")
        h_MC_gen.GetYaxis().SetTitle("d#it{N}/d(#it{c}t) (cm^{-1})")
        leg = ROOT.TLegend()
        leg.SetX1NDC(0.3)
        leg.SetY1NDC(0.4)
        leg.AddEntry(h_MC_gen,"generated","l")
        leg.AddEntry(h,"reconstructed","l")
        h.SetLineColor(ROOT.kRed)
        h.Draw("same")
        leg.Draw("same")
        # lifetime.Draw("same")
        h.Write()
        h_MC_gen.Write()
        c.Write()
        del h, h_MC_gen
