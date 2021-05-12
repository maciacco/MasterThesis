#!/usr/bin/env python3
import os
import pickle
import warnings

import numpy as np
import ROOT
import yaml

HE_3_MASS = 2.809230089

##################################################################
# read configuration file
##################################################################
config = 'config.yaml'
with open(os.path.expandvars(config), 'r') as stream:
    try:
        params = yaml.full_load(stream)
    except yaml.YAMLError as exc:
        print(exc)

CT_BINS_CENT = params['CT_BINS_CENT']
CENTRALITY_LIST = params['CENTRALITY_LIST']
##################################################################

split_list = ['antimatter', 'matter']
cent_bins_MB = [[0, 10], [10, 40], [40, 90]]

# mc input file
mc_file = './AnalysisResults.root'
outfile = ROOT.TFile("He3_abs.root", "recreate")
centfile = ROOT.TFile("AnalysisResults_18.root")

# get event centrality distribution
cent_dist = centfile.Get("Centrality_selected")
cent_dist_max = cent_dist.GetMaximum()

# functions
n_fun = [1, 1, 1, 1]  # [5, 5, 5, 4]
func = [[] for _ in range(len(CENTRALITY_LIST)-1)]
funcMB = [[] for _ in range(3)]
func_max = [[] for _ in range(len(CENTRALITY_LIST)-1)]
func_max_MB = [[] for _ in range(3)]

# functions names
func_names = ["BGBW", "Boltzmann", "Mt-exp", "Pt-exp", "LevyTsallis"]
func_names_MB = ["BlastWave", "Boltzmann", "LevyTsallis", "Mt-exp"]

# functions input files
input_func_file = ROOT.TFile("Anti_fits.root")
input_func_file_MB = ROOT.TFile("BlastWaveFits.root")

# get functions and maxima from file
cent_index = [1, 2, 4]
for i_cent in range(len(CENTRALITY_LIST)-1):
    for i_fun in range(n_fun[i_cent]):
        func[i_cent].append(input_func_file.Get(
            f"{func_names[i_fun]}/{cent_index[i_cent]}/{func_names[i_fun]}{cent_index[i_cent]}"))
        func_max[i_cent].append(func[i_cent][i_fun].GetMaximum())
for i_cent in range(3):
    for i_fun in range(n_fun[-1]):
        funcMB[i_cent].append(input_func_file_MB.Get(f"{func_names_MB[i_fun]}/{func_names_MB[i_fun]}{i_cent}"))
        func_max_MB[i_cent].append(funcMB[i_cent][i_fun].GetMaximum())

# book histograms
cent_len = len(CENTRALITY_LIST)
h_abs_radius = [[] for _ in range(cent_len)]
h_abs_ct = [[] for _ in range(cent_len)]
h_gen_radius = [[] for _ in range(cent_len)]
h_gen_ct = [[] for _ in range(cent_len)]
h_gen_pt = [[] for _ in range(cent_len)]
h_rec_radius = [[] for _ in range(cent_len)]
h_rec_ct = [[] for _ in range(cent_len)]
h_rec_pt = [[] for _ in range(cent_len)]
for i_cent in range(cent_len):
    for i_fun in range(n_fun[i_cent]):
        ct_bins = np.asarray(CT_BINS_CENT[i_cent], dtype="float")
        cent_bins = CENTRALITY_LIST[i_cent]
        f_name = func_names[i_fun]
        if cent_bins[1] == 90:
            f_name = func_names_MB[i_fun]
        h_abs_radius[i_cent].append([])
        h_abs_ct[i_cent].append([])
        h_gen_radius[i_cent].append([])
        h_gen_ct[i_cent].append([])
        h_gen_pt[i_cent].append([])
        h_rec_radius[i_cent].append([])
        h_rec_ct[i_cent].append([])
        h_rec_pt[i_cent].append([])
        for split in split_list:
            h_abs_radius[i_cent][i_fun].append(ROOT.TH1D(
                f"fAbsRadius_{split}_{cent_bins[0]}_{cent_bins[1]}_{f_name}", ";#it{R}_{#it{abs}} (cm);Entries",
                1000, 0, 1000))
            h_abs_ct[i_cent][i_fun].append(ROOT.TH1D(
                f"fAbsCt_{split}_{cent_bins[0]}_{cent_bins[1]}_{f_name}", ";#it{c}t (cm);Entries",
                len(ct_bins)-1,ct_bins))
            h_gen_radius[i_cent][i_fun].append(ROOT.TH1D(
                f"fGenRadius_{split}_{cent_bins[0]}_{cent_bins[1]}_{f_name}", ";#it{R}_{#it{abs}} (cm);Entries",
                1000, 0, 1000))
            h_gen_ct[i_cent][i_fun].append(ROOT.TH1D(
                f"fGenCt_{split}_{cent_bins[0]}_{cent_bins[1]}_{f_name}", ";#it{c}t (cm);Entries",
                len(ct_bins)-1,ct_bins))
            h_gen_pt[i_cent][i_fun].append(ROOT.TH1D(
                f"fGenPt_{split}_{cent_bins[0]}_{cent_bins[1]}_{f_name}", ";#it{p}_{T} (GeV/#it{c});Entries",
                50, 0, 10))
            h_rec_radius[i_cent][i_fun].append(ROOT.TH1D(
                f"fGenRadius_{split}_{cent_bins[0]}_{cent_bins[1]}_{f_name}", ";#it{R}_{#it{abs}} (cm);Entries",
                1000, 0, 1000))
            h_rec_ct[i_cent][i_fun].append(ROOT.TH1D(
                f"fRecCt_{split}_{cent_bins[0]}_{cent_bins[1]}_{f_name}", ";#it{c}t (cm);Entries",
                len(ct_bins)-1,ct_bins))
            h_rec_pt[i_cent][i_fun].append(ROOT.TH1D(
                f"fRecPt_{split}_{cent_bins[0]}_{cent_bins[1]}_{f_name}", ";#it{p}_{T} (GeV/#it{c});Entries",
                50, 0, 10))

# read tree
data_frame_he3 = ROOT.RDataFrame('STree', mc_file)
data_frame_he3 = data_frame_he3.Filter('pt > 2. and pt < 10. and (flag & 1)==1')
np_he3 = data_frame_he3.AsNumpy(["pt", "pdg", "absCt", "eta"])

# analysis in centrality classes
for he3 in zip(np_he3['pt'], np_he3['pdg'], np_he3['absCt'], np_he3['eta']):
    i_matt = 0
    if he3[1] == 1000020030:
        i_matt = 1
    absCt = he3[2]

    # analysis in centrality classes
    for i_cent in range(len(CENTRALITY_LIST)-1):
        for i_fun in range(n_fun[i_cent]):
            # rejection sampling to reweight pt
            if ROOT.gRandom.Rndm()*func_max[i_cent][i_fun] > func[i_cent][i_fun].Eval(he3[0]):
                continue

            # sample decay ct and ckeck for absorption
            decCt = ROOT.gRandom.Exp(7.6)
            # polar angle from eta
            tmp = abs(he3[3])
            tmp = ROOT.TMath.Exp(tmp)
            theta = 2*ROOT.TMath.ATan(tmp) # eta = -log[tan(theta/2)]
            # momentum from transverse momentum and angle
            mom = he3[0]/ROOT.TMath.Sin(theta)
            # absorption radius
            abs_radius = absCt*mom/HE_3_MASS
            # decay radius
            dec_radius = decCt*mom/HE_3_MASS

            h_abs_ct[i_cent][i_fun][i_matt].Fill(absCt)
            h_abs_radius[i_cent][i_fun][i_matt].Fill(abs_radius)
            h_gen_radius[i_cent][i_fun][i_matt].Fill(dec_radius)
            h_gen_ct[i_cent][i_fun][i_matt].Fill(decCt)
            h_gen_pt[i_cent][i_fun][i_matt].Fill(he3[0])
            if not (decCt > absCt):  # decCt < absCt
                h_rec_radius[i_cent][i_fun][i_matt].Fill(dec_radius)
                h_rec_ct[i_cent][i_fun][i_matt].Fill(decCt)
                h_rec_pt[i_cent][i_fun][i_matt].Fill(he3[0])
            if (absCt < -0.5):  # decCt < absCt
                h_rec_radius[i_cent][i_fun][i_matt].Fill(dec_radius)
                h_rec_ct[i_cent][i_fun][i_matt].Fill(decCt)
                h_rec_pt[i_cent][i_fun][i_matt].Fill(he3[0])

    # MB analysis
    # sample centrality
    sample_centrality = ROOT.gRandom.Rndm()*90.
    uniform_event = ROOT.gRandom.Rndm()*cent_dist_max
    if uniform_event > cent_dist.GetBinContent(cent_dist.FindBin(sample_centrality)):
        continue
    i_cent_MB = 0
    if sample_centrality > 10. and sample_centrality < 40.:
        i_cent_MB = 1
    elif sample_centrality > 40.:
        i_cent_MB = 2
    for i_fun in range(n_fun[i_cent_MB]):
        # rejection sampling to reweight pt
        if ROOT.gRandom.Rndm()*func_max_MB[i_cent_MB][i_fun] > funcMB[i_cent_MB][i_fun].Eval(he3[0]):
            continue

        # sample decay ct and ckeck for absorption
        decCt = ROOT.gRandom.Exp(7.6)
        h_gen_ct[3][i_fun][i_matt].Fill(decCt)
        h_gen_pt[3][i_fun][i_matt].Fill(he3[0])
        if (decCt < absCt) or (absCt < -0.5):  # decCt < absCt
            h_rec_ct[3][i_fun][i_matt].Fill(decCt)
            h_rec_pt[3][i_fun][i_matt].Fill(he3[0])

# write histograms and compute efficiency
for i_cent in range(cent_len):
    cent_bins = CENTRALITY_LIST[i_cent]
    outfile.mkdir(f"{cent_bins[0]}_{cent_bins[1]}")
    outfile.cd(f"{cent_bins[0]}_{cent_bins[1]}")
    for i_fun in range(n_fun[i_cent]):
        func_name = func_names[i_fun]
        if cent_bins[1] == 90:
            func_name = func_names_MB[i_fun]
        outfile.mkdir(f"{cent_bins[0]}_{cent_bins[1]}/{func_name}")
        outfile.cd(f"{cent_bins[0]}_{cent_bins[1]}/{func_name}")
        for i_matt in [0, 1]:
            h_abs_radius[i_cent][i_fun][i_matt].Write()
            h_abs_ct[i_cent][i_fun][i_matt].Write()

            h_rec_radius[i_cent][i_fun][i_matt].Write()
            h_gen_radius[i_cent][i_fun][i_matt].Write()
            eff_radius = ROOT.TGraphAsymmErrors(h_rec_radius[i_cent][i_fun][i_matt], h_gen_radius[i_cent][i_fun][i_matt])
            eff_radius.GetXaxis().SetTitle("#it{R}_{#it{abs}} (cm)")
            eff_radius.GetYaxis().SetTitle("1 - #it{f}_{abs}")
            eff_radius.Write(f"fEffRadius_{split_list[i_matt]}_{cent_bins[0]}_{cent_bins[1]}_{func_name}")

            h_rec_ct[i_cent][i_fun][i_matt].Write()
            h_gen_ct[i_cent][i_fun][i_matt].Write()
            eff_ct = ROOT.TGraphAsymmErrors(h_rec_ct[i_cent][i_fun][i_matt], h_gen_ct[i_cent][i_fun][i_matt])
            eff_ct.GetXaxis().SetTitle("#it{c}t (cm)")
            eff_ct.GetYaxis().SetTitle("1 - #it{f}_{abs}")
            eff_ct.Write(f"fEffCt_{split_list[i_matt]}_{cent_bins[0]}_{cent_bins[1]}_{func_name}")

            h_rec_pt[i_cent][i_fun][i_matt].Write()
            h_gen_pt[i_cent][i_fun][i_matt].Write()
            eff_pt = ROOT.TGraphAsymmErrors(h_rec_pt[i_cent][i_fun][i_matt], h_gen_pt[i_cent][i_fun][i_matt])
            eff_pt.GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c})")
            eff_pt.GetYaxis().SetTitle("1 - #it{f}_{abs}")
            eff_pt.Write(f"fEffPt_{split_list[i_matt]}_{cent_bins[0]}_{cent_bins[1]}_{func_name}")

outfile.Close()
