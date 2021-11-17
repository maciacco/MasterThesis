#!/usr/bin/env python3
import os
import numpy as np
import ROOT
import yaml

def he3_correction_pt(i_matt, pt):
    if i_matt == 0:
        return 1.04948*ROOT.TMath.Power(pt,-0.01525)
    return 0.99274*ROOT.TMath.Power(pt,0.00143)
    
def he3_uncertainty_pt(i_matt, pt):
    if i_matt == 0:
        return 0.02088*ROOT.TMath.Power(pt,-0.48766)
    return 0.00294*ROOT.TMath.Power(pt,-0.19483)

TOY_MC_EFF = 0.8 # value used in the toy MC (no physical meaning, just to keep eff < 1)
TOY_MC_CT = 7.6  # value of the proper decay length used in the toy MC
N_TRIALS = 1e7   # number of iterations for the toy MC 

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

if not os.path.isdir(f'plots/efficiency_correction'):
    os.mkdir(f'plots/efficiency_correction')

split_list = ['antimatter', 'matter']
cent_bins_MB = [[0, 10], [10, 40], [40, 90]]

# mc input file
mc_file = './AnalysisResults.root'
outfile = ROOT.TFile("EffAbsCorrection.root", "recreate")
centfile = ROOT.TFile("../data/Hypertriton_PbPb/AnalysisResults_18.root")

# get event centrality distribution
cent_dist = centfile.Get("Centrality_selected")
cent_dist_max = cent_dist.GetMaximum()

# functions
n_fun = [1, 1, 1, 1, 1]  # [5, 5, 5, 4, 4]
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
for i_cent in range(len(CENTRALITY_LIST)-2):
    for i_fun in range(n_fun[i_cent]):
        func[i_cent].append(input_func_file.Get(
            f"{func_names[i_fun]}/{cent_index[i_cent]}/{func_names[i_fun]}{cent_index[i_cent]}"))
        func_max[i_cent].append(func[i_cent][i_fun].GetMaximum())

for i_cent in range(3):
    for i_fun in range(n_fun[-1]):
        funcMB[i_cent].append(input_func_file_MB.Get(f"{func_names_MB[i_fun]}/{func_names_MB[i_fun]}{i_cent}"))
        func_max_MB[i_cent].append(funcMB[i_cent][i_fun].GetMaximum())

# 0-10%
i_cent = 0
for i_fun in range(n_fun[-2]):
    func[3].append(input_func_file_MB.Get(f"{func_names_MB[i_fun]}/{func_names_MB[i_fun]}{i_cent}"))
    func_max[3].append(funcMB[i_cent][i_fun].GetMaximum())

# book histograms
cent_len = len(CENTRALITY_LIST)
h_eff_correction_ct = [[] for _ in range(cent_len)]
h_eff_correction_pt = [[] for _ in range(cent_len)] # check consistency
h_eff_correction_ct_syst = [[] for _ in range(cent_len)]
h_eff_correction_pt_syst = [[] for _ in range(cent_len)] # check consistency
h_gen_ct = [[] for _ in range(cent_len)]
h_gen_pt = [[] for _ in range(cent_len)]
h_rec_ct = [[] for _ in range(cent_len)]
h_rec_pt = [[] for _ in range(cent_len)]
h_rec_ct_corrected = [[] for _ in range(cent_len)]
h_rec_pt_corrected = [[] for _ in range(cent_len)]
h_rec_ct_syst = [[] for _ in range(cent_len)]
h_rec_pt_syst = [[] for _ in range(cent_len)]
for i_cent in range(cent_len):
    for i_fun in range(n_fun[i_cent]):
        ct_bins = np.asarray(CT_BINS_CENT[i_cent], dtype="float")
        cent_bins = CENTRALITY_LIST[i_cent]
        f_name = func_names[i_fun]
        if (cent_bins[1] == 90):
            f_name = func_names_MB[i_fun]
        if (cent_bins[0] == 0) and (cent_bins[1] == 10):
            f_name = func_names_MB[i_fun]
        h_eff_correction_ct[i_cent].append([])
        h_eff_correction_pt[i_cent].append([])
        h_eff_correction_ct_syst[i_cent].append([])
        h_eff_correction_pt_syst[i_cent].append([])
        h_gen_ct[i_cent].append([])
        h_gen_pt[i_cent].append([])
        h_rec_ct[i_cent].append([])
        h_rec_pt[i_cent].append([])
        h_rec_ct_corrected[i_cent].append([])
        h_rec_pt_corrected[i_cent].append([])
        h_rec_ct_syst[i_cent].append([])
        h_rec_pt_syst[i_cent].append([])
        for split in split_list:
            h_eff_correction_ct[i_cent][i_fun].append(ROOT.TH1D(
                f"fEffCorrectionCt_{split}_{cent_bins[0]}_{cent_bins[1]}_{f_name}", ";#it{c}t (cm);Entries",
                len(ct_bins)-1,ct_bins))
            h_eff_correction_pt[i_cent][i_fun].append(ROOT.TH1D(
                f"fEffCorrectionPt_{split}_{cent_bins[0]}_{cent_bins[1]}_{f_name}", ";#it{p}_{T} (GeV/#it{c});Entries",
                50, 0, 10))
            h_eff_correction_ct_syst[i_cent][i_fun].append(ROOT.TH1D(
                f"fEffCorrectionCtSyst_{split}_{cent_bins[0]}_{cent_bins[1]}_{f_name}", ";#it{c}t (cm);Entries",
                len(ct_bins)-1,ct_bins))
            h_eff_correction_pt_syst[i_cent][i_fun].append(ROOT.TH1D(
                f"fEffCorrectionPtSyst_{split}_{cent_bins[0]}_{cent_bins[1]}_{f_name}", ";#it{p}_{T} (GeV/#it{c});Entries",
                50, 0, 10))
            h_gen_ct[i_cent][i_fun].append(ROOT.TH1D(
                f"fGenCt_{split}_{cent_bins[0]}_{cent_bins[1]}_{f_name}", ";#it{c}t (cm);Entries",
                len(ct_bins)-1,ct_bins))
            h_gen_pt[i_cent][i_fun].append(ROOT.TH1D(
                f"fGenPt_{split}_{cent_bins[0]}_{cent_bins[1]}_{f_name}", ";#it{p}_{T} (GeV/#it{c});Entries",
                50, 0, 10))
            h_rec_ct[i_cent][i_fun].append(ROOT.TH1D(
                f"fRecCt_{split}_{cent_bins[0]}_{cent_bins[1]}_{f_name}", ";#it{c}t (cm);Entries",
                len(ct_bins)-1,ct_bins))
            h_rec_pt[i_cent][i_fun].append(ROOT.TH1D(
                f"fRecPt_{split}_{cent_bins[0]}_{cent_bins[1]}_{f_name}", ";#it{p}_{T} (GeV/#it{c});Entries",
                50, 0, 10))
            h_rec_ct_corrected[i_cent][i_fun].append(ROOT.TH1D(
                f"fRecCtCorrected_{split}_{cent_bins[0]}_{cent_bins[1]}_{f_name}", ";#it{c}t (cm);Entries",
                len(ct_bins)-1,ct_bins))
            h_rec_pt_corrected[i_cent][i_fun].append(ROOT.TH1D(
                f"fRecPtCorrected_{split}_{cent_bins[0]}_{cent_bins[1]}_{f_name}", ";#it{p}_{T} (GeV/#it{c});Entries",
                50, 0, 10))
            h_rec_ct_syst[i_cent][i_fun].append(ROOT.TH1D(
                f"fRecCtSyst_{split}_{cent_bins[0]}_{cent_bins[1]}_{f_name}", ";#it{c}t (cm);Entries",
                len(ct_bins)-1,ct_bins))
            h_rec_pt_syst[i_cent][i_fun].append(ROOT.TH1D(
                f"fRecPtSyst_{split}_{cent_bins[0]}_{cent_bins[1]}_{f_name}", ";#it{p}_{T} (GeV/#it{c});Entries",
                50, 0, 10))

# begin toy MC
trial = 0
while trial < N_TRIALS:

    for i_matt in [0, 1]:

        # analysis in centrality classes
        for i_cent in range(len(CENTRALITY_LIST)-1):
            for i_fun in range(n_fun[i_cent]):
                
                # sample pt in [2, 10] GeV/c
                pt = func[i_cent][i_fun].GetRandom(2, 10)

                # sample decay ct 
                dec_ct = ROOT.gRandom.Exp(TOY_MC_CT)

                # efficiency correction
                eff_corr = he3_correction_pt(i_matt,pt)
                
                # absorption systematic error
                abs_syst = he3_uncertainty_pt(i_matt,pt)

                h_gen_ct[i_cent][i_fun][i_matt].Fill(dec_ct)
                h_gen_pt[i_cent][i_fun][i_matt].Fill(pt)
                if (ROOT.gRandom.Rndm() < TOY_MC_EFF):
                    h_rec_ct[i_cent][i_fun][i_matt].Fill(dec_ct)
                    h_rec_pt[i_cent][i_fun][i_matt].Fill(pt)
                if (ROOT.gRandom.Rndm() < TOY_MC_EFF*eff_corr):
                    h_rec_ct_corrected[i_cent][i_fun][i_matt].Fill(dec_ct)
                    h_rec_pt_corrected[i_cent][i_fun][i_matt].Fill(pt)
                if (ROOT.gRandom.Rndm() < TOY_MC_EFF*(1+abs_syst)):
                    h_rec_ct_syst[i_cent][i_fun][i_matt].Fill(dec_ct)
                    h_rec_pt_syst[i_cent][i_fun][i_matt].Fill(pt)

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

            # sample pt in [2, 10] GeV/c
            pt = funcMB[i_cent_MB][i_fun].GetRandom(2, 10)

            # sample decay ct 
            dec_ct = ROOT.gRandom.Exp(TOY_MC_CT)

            # efficiency correction
            eff_corr = he3_correction_pt(i_matt,pt)

            h_gen_ct[4][i_fun][i_matt].Fill(dec_ct)
            h_gen_pt[4][i_fun][i_matt].Fill(pt)
            if (ROOT.gRandom.Rndm()*eff_corr < TOY_MC_EFF):
                h_rec_ct[4][i_fun][i_matt].Fill(dec_ct)
                h_rec_pt[4][i_fun][i_matt].Fill(pt)
            if (ROOT.gRandom.Rndm()*eff_corr < TOY_MC_EFF*eff_corr):
                h_rec_ct_corrected[4][i_fun][i_matt].Fill(dec_ct)
                h_rec_pt_corrected[4][i_fun][i_matt].Fill(pt)
                
    trial += 1
    
# write histograms and compute efficiency
for i_cent in range(cent_len):
    cent_bins = CENTRALITY_LIST[i_cent]
    outfile.mkdir(f"{cent_bins[0]}_{cent_bins[1]}")
    outfile.cd(f"{cent_bins[0]}_{cent_bins[1]}")
    for i_fun in range(n_fun[i_cent]):
        func_name = func_names[i_fun]
        if cent_bins[1] == 90:
            func_name = func_names_MB[i_fun]
        if (cent_bins[0] == 0) and  (cent_bins[1] == 10):
            func_name = func_names_MB[i_fun]
        outfile.mkdir(f"{cent_bins[0]}_{cent_bins[1]}/{func_name}")
        outfile.cd(f"{cent_bins[0]}_{cent_bins[1]}/{func_name}")
        for i_matt in [0, 1]:
            # define consant function
            f_const = ROOT.TF1("f","1",0,50)

            # ct
            h_rec_ct[i_cent][i_fun][i_matt].Write()
            h_gen_ct[i_cent][i_fun][i_matt].Write()
            eff_ct = ROOT.TH1D(h_rec_ct[i_cent][i_fun][i_matt])
            eff_ct.Divide(h_rec_ct[i_cent][i_fun][i_matt],h_gen_ct[i_cent][i_fun][i_matt],1,1,"B")
            eff_ct.GetXaxis().SetTitle("#it{c}t (cm)")
            eff_ct.GetYaxis().SetTitle("efficiency")
            eff_ct.SetTitle(f"{cent_bins[0]}-{cent_bins[1]}%")
            eff_ct.Write(f"fEffCt_{split_list[i_matt]}_{cent_bins[0]}_{cent_bins[1]}_{func_name}")
            
            h_rec_ct_corrected[i_cent][i_fun][i_matt].Write()
            eff_ct_corrected = ROOT.TH1D(h_rec_ct[i_cent][i_fun][i_matt])
            eff_ct_corrected.Divide(h_rec_ct_corrected[i_cent][i_fun][i_matt],h_gen_ct[i_cent][i_fun][i_matt],1,1,"B")
            eff_ct_corrected.GetXaxis().SetTitle("#it{c}t (cm)")
            eff_ct_corrected.GetYaxis().SetTitle("efficiency corrected")
            eff_ct_corrected.SetTitle(f"{cent_bins[0]}-{cent_bins[1]}%")
            eff_ct_corrected.Write(f"fEffCtCorrected_{split_list[i_matt]}_{cent_bins[0]}_{cent_bins[1]}_{func_name}")

            h_correction = ROOT.TH1D(h_rec_ct[i_cent][i_fun][i_matt])
            h_correction.Divide(h_rec_ct_corrected[i_cent][i_fun][i_matt],h_rec_ct[i_cent][i_fun][i_matt])
            h_correction.GetXaxis().SetTitle("#it{c}t (cm)")
            h_correction.GetYaxis().SetTitle("correction")
            h_correction.SetTitle(f"{cent_bins[0]}-{cent_bins[1]}%")
            h_correction.Write(f"fCorrection_{split_list[i_matt]}_{cent_bins[0]}_{cent_bins[1]}_{func_name}")
            
            h_correction_syst = ROOT.TH1D(h_rec_ct[i_cent][i_fun][i_matt])
            h_correction_syst.Divide(h_rec_ct_syst[i_cent][i_fun][i_matt],h_rec_ct[i_cent][i_fun][i_matt])
            h_correction_syst.Add(f_const,-1)
            h_correction_syst.GetXaxis().SetTitle("#it{c}t (cm)")
            h_correction_syst.GetYaxis().SetTitle("relative systematic uncertainty")
            h_correction_syst.SetTitle(f"{cent_bins[0]}-{cent_bins[1]}%")
            h_correction_syst.Write(f"fRelativeSystematicUncertainty_{split_list[i_matt]}_{cent_bins[0]}_{cent_bins[1]}_{func_name}")

            # plot the correction
            canv = ROOT.TCanvas("canv", "canv")
            canv.cd()
            canv.SetTicks(1, 1)
            h_correction.Draw()
            canv.Print(f"plots/efficiency_correction/fCorrection_{split_list[i_matt]}_{cent_bins[0]}_{cent_bins[1]}_{func_name}.pdf")

            # pt
            h_rec_pt[i_cent][i_fun][i_matt].Write()
            h_gen_pt[i_cent][i_fun][i_matt].Write()
            eff_pt = ROOT.TH1D(h_rec_pt[i_cent][i_fun][i_matt])
            eff_pt.Divide(h_rec_pt[i_cent][i_fun][i_matt],h_gen_pt[i_cent][i_fun][i_matt],1,1,"B")
            eff_pt.GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c})")
            eff_pt.GetYaxis().SetTitle("efficiency")
            eff_ct.SetTitle(f"{cent_bins[0]}-{cent_bins[1]}%")
            eff_pt.Write(f"fEffPt_{split_list[i_matt]}_{cent_bins[0]}_{cent_bins[1]}_{func_name}")

outfile.Close()