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
ROOT.gROOT.SetBatch()

bkg_function = ['pol','expo']
cut_var = [0.1,0.1,0.15]
dof = [13,13,20]
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

raw_yields_file = ROOT.TFile('SignalExtraction-data_fineCt_lowEff.root')
eff_array = np.arange(0.10, MAX_EFF, 0.01)
f = ROOT.TFile("cut_selection.root","recreate")
for i_cent_bins in range(len(CENTRALITY_LIST)):

    cent_bins = CENTRALITY_LIST[i_cent_bins]
    h_cut_distr = ROOT.TH1D(f"h_distr_{cent_bins[0]}_{cent_bins[1]}",";cut;Entries",100,0,1)
    for split in SPLIT_LIST:
        for ct_bins in zip(CT_BINS_CENT[i_cent_bins][:-1], CT_BINS_CENT[i_cent_bins][1:]):
            if ct_bins[0] < 1 or ct_bins[1] > 40:
                continue
            h_compare_slope = ROOT.TH1D(f"compare_slope_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}",";cut BDT efficiency;slope",100,0,1)
            h_chi2 = ROOT.TH1D(f"chi2_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}",";cut BDT efficiency;chi2",100,0,1)
            h_raw = raw_yields_file.Get(f"{split}_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}_pol1/fRawYields")
            ref_pol1 = ROOT.TF1("ref_pol1","pol1")
            ref_pol1.FixParameter(0,0.)
            ref_pol1_1 = ROOT.TF1("ref_pol1_1","pol1")
            ref_pol1_1.FixParameter(0,0.)
            h_raw.Fit("ref_pol1","R","",0.9,1.00)
            h_raw.Fit("ref_pol1_1","R","",0.7,1)
            for cut_eff in eff_array:
                if cut_eff<0.3:
                    continue
                h_raw.Fit("pol1","R","",cut_eff-cut_var[i_cent_bins],cut_eff+cut_var[i_cent_bins])
                if h_raw.GetFunction("pol1").GetNDF()<dof[i_cent_bins]:
                    continue
                # if h_raw.GetFunction("pol1").GetProb()>0.975 or h_raw.GetFunction("pol1").GetProb()<0.025:
                #     continue
                if (h_raw.GetFunction("pol1").GetChisquare()/h_raw.GetFunction("pol1").GetNDF())>2:
                    continue
                h_compare_slope.SetBinContent(int(cut_eff*100)+1,h_raw.GetFunction("pol1").GetParameter(1))
                h_compare_slope.SetBinError(int(cut_eff*100)+1,h_raw.GetFunction("pol1").GetParError(1))
                h_chi2.SetBinContent(int(cut_eff*100)+1,h_raw.GetFunction("pol1").GetChisquare()/h_raw.GetFunction("pol1").GetNDF())
                f.mkdir(f"fits_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}")
                f.cd(f"fits_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}")
                c = ROOT.TCanvas(f"fit_{'{:.2f}'.format(cut_eff)}_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}")
                l_left = ROOT.TLine(cut_eff-cut_var[i_cent_bins],0,cut_eff-cut_var[i_cent_bins],h_raw.GetMaximum())
                l_right = ROOT.TLine(cut_eff+cut_var[i_cent_bins],0,cut_eff+cut_var[i_cent_bins],h_raw.GetMaximum())
                l_cent = ROOT.TLine(cut_eff,0,cut_eff,h_raw.GetMaximum())
                h_raw.GetFunction("pol1").SetRange(0,1)
                c.cd()
                h_raw.Draw()
                l_left.SetLineWidth(2)
                l_left.SetLineStyle(ROOT.kDashed)
                l_left.SetLineColor(ROOT.kBlue)
                l_right.SetLineWidth(2)
                l_right.SetLineStyle(ROOT.kDashed)
                l_right.SetLineColor(ROOT.kBlue)
                l_cent.SetLineWidth(2)
                l_cent.SetLineStyle(ROOT.kDashed)
                l_cent.SetLineColor(ROOT.kGreen)
                l_left.Draw("same")
                l_right.Draw("same")
                l_cent.Draw("same")
                c.Write()
                if (h_raw.GetFunction("pol1").GetParameter(1)-ref_pol1_1.GetParameter(1))>0 and (h_raw.GetFunction("pol1").GetParameter(1)-ref_pol1.GetParameter(1))<0:
                    h_cut_distr.Fill(cut_eff,cut_eff)
            c = ROOT.TCanvas(f"c_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}",f"c_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}")
            c.cd()
            h_compare_slope.Draw("pe")
            ref = ROOT.TF1("ref","pol0",0,1)
            ref_2 = ROOT.TF1("ref_2","pol0",0,1)
            ref.SetParameter(0,ref_pol1.GetParameter(1))
            ref_2.SetParameter(0,ref_pol1_1.GetParameter(1))
            h_chi2.Write()
            ref.Draw("same")
            ref_2.Draw("same")
            f.cd()
            c.Write()
    
    f.cd()
    h_cut_distr.Write()
f.Close()