#!/usr/bin/env python3
import os
import pickle
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import ROOT
import uproot
import yaml

SPEED_OF_LIGHT = 2.99792458
N_TRIALS = 10000
SPLIT = False
MAX_EFF = 0.99
THRESH_EFF = 0.90

# avoid pandas warning
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

ANALYSIS_RESULTS_PATH = params['ANALYSIS_RESULTS_PATH']
CT_BINS_CENT = params['CT_BINS_CENT']
CENTRALITY_LIST = params['CENTRALITY_LIST']
RANDOM_STATE = params['RANDOM_STATE']
##################################################################

# split matter/antimatter
SPLIT_LIST = ['all']
if SPLIT:
    SPLIT_LIST = ['antimatter', 'matter']

eff_cut_dict = pickle.load(open("file_eff_cut_dict", "rb"))
presel_eff_file = uproot.open('PreselEff.root')
signal_extraction_file = ROOT.TFile.Open('SignalExtraction.root')
signal_extraction_keys = uproot.open('SignalExtraction.root').keys()

# abs_correction_file = ROOT.TFile.Open('He3_abs.root')
systematics_file = ROOT.TFile.Open('SystematicsDeltaMass.root', 'recreate')

for i_cent_bins in range(len(CENTRALITY_LIST)):
    cent_bins = CENTRALITY_LIST[i_cent_bins]
    if cent_bins[1] < 80:
        continue

    h_m_minus_mc_distribution = ROOT.TH1D(f'fDeltaMassDistribution_{cent_bins[0]}_{cent_bins[1]}', f'{cent_bins[0]}-{cent_bins[1]}%', 120, 0, 60)
    h_asymmetry_distribution = ROOT.TH1D(f'fAsymmetryDistribution_{cent_bins[0]}_{cent_bins[1]}', f'{cent_bins[0]}-{cent_bins[1]}%', 400, -200, 200)
    h_prob_distribution = ROOT.TH1D(f'fProbDistribution_{cent_bins[0]}_{cent_bins[1]}', f'{cent_bins[0]}-{cent_bins[1]}%', 100, 0, 1)
    h_lifetime = [ROOT.TH1D(), ROOT.TH1D()]
    h_lifetime[0] = ROOT.TH1D(f'fLifetimeAntimatterDistribution_{cent_bins[0]}_{cent_bins[1]}', f'{cent_bins[0]}-{cent_bins[1]}%', 1000, 0, 1000)
    h_lifetime[1] = ROOT.TH1D(f'fLifetimeMatterDistribution_{cent_bins[0]}_{cent_bins[1]}', f'{cent_bins[0]}-{cent_bins[1]}%', 1000, 0, 1000)
    h_fit_status = ROOT.TH1D(f'fFitStatus_{cent_bins[0]}_{cent_bins[1]}', f'{cent_bins[0]}-{cent_bins[1]}%', 100, 0, 99)

    systematics_file.mkdir(f'{cent_bins[0]}_{cent_bins[1]}')
    ####################################################
    # MULTITRIAL ESTIMATION OF SYSTEMATIC UNCERTAINTY
    ####################################################
    # sample 10000 combinations
    i_trial=0
    while i_trial < N_TRIALS:
        h_corrected_yields = [ROOT.TH1D(), ROOT.TH1D()]
        h_mass = [ROOT.TH1D(), ROOT.TH1D()]
        h_mass_mc = [ROOT.TH1D(), ROOT.TH1D()]
        h_mass_dat = [ROOT.TH1D(), ROOT.TH1D()]
        h_m_minus_mc = [ROOT.TH1D(), ROOT.TH1D()]
        h_m_minus_mc_distribution_reject = [ROOT.TH1D(), ROOT.TH1D()] # distribution to reject outliers
        if (i_trial % 100) == 0:
            print(f'{cent_bins[0]}-{cent_bins[1]}%: i_trial = {i_trial}')

        lifetime_tmp = [-9999., -9999.]
        for i_split, split in enumerate(SPLIT_LIST):

            # get preselection efficiency histogram
            presel_eff_counts, presel_eff_edges = presel_eff_file[
                f'fPreselEff_vs_ct_{split}_{cent_bins[0]}_{cent_bins[1]};1'].to_numpy()
            presel_eff_bin_centers = (presel_eff_edges[1:]+presel_eff_edges[:-1])/2

            # list of corrected yields
            ct_bins_tmp = [0]
            ct_bins_tmp += CT_BINS_CENT[i_cent_bins]
            #if cent_bins[1] == 90:
            #    ct_bin_tmp = CT_BINS_CENT[i_cent_bins]
            bins = np.array(ct_bins_tmp, dtype=float)

            h_corrected_yields[i_split] = ROOT.TH1D(
                f'fYields_{split}_{cent_bins[0]}_{cent_bins[1]}', f'{split}, {cent_bins[0]}-{cent_bins[1]}%', len(bins)-1, bins)
            h_mass[i_split] = ROOT.TH1D(
                f'fMass_{split}_{cent_bins[0]}_{cent_bins[1]}', f'{split}, {cent_bins[0]}-{cent_bins[1]}%', len(bins)-1, bins)
            h_mass_mc[i_split] = ROOT.TH1D(
                f'fMassMC_{split}_{cent_bins[0]}_{cent_bins[1]}', f'{split}, {cent_bins[0]}-{cent_bins[1]}%', len(bins)-1, bins)
            h_mass_dat[i_split] = ROOT.TH1D(
                f'fMassData_{split}_{cent_bins[0]}_{cent_bins[1]}', f'{split}, {cent_bins[0]}-{cent_bins[1]}%', len(bins)-1, bins)
            h_m_minus_mc[i_split] = ROOT.TH1D(
                f'fDeltaMassLambda_{split}_{cent_bins[0]}_{cent_bins[1]}', f'{split}, {cent_bins[0]}-{cent_bins[1]}%', len(bins)-1, bins)
            h_m_minus_mc_distribution_reject[i_split] = ROOT.TH1D(
                f'fDeltaMassLambdaDistribution_{split}_{cent_bins[0]}_{cent_bins[1]}', f'{split}, {cent_bins[0]}-{cent_bins[1]}%', 800, -400, 400)

            for ct_bins in zip(CT_BINS_CENT[i_cent_bins][:-1], CT_BINS_CENT[i_cent_bins][1:]):

                bin = f'{split}_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}'
                formatted_eff_cut = "0.80"#"{:.2f}".format(eff_cut_dict[bin])

                # look for plot with eff = eff_cut (or the nearest one)
                bkg_shape = 'pol1'
                # if ROOT.gRandom.Rndm() > 0.5:
                #     bkg_shape = 'expo'

                eff_cut_increment = 0
                eff_cut_sign = -1
                while signal_extraction_keys.count(f'{bin}_{bkg_shape}/fInvMass_{formatted_eff_cut};1') == 0:
                    if eff_cut_sign == -1:
                        eff_cut_increment += 0.01
                    eff_cut_sign *= -1
                    formatted_eff_cut = "{:.2f}".format(eff_cut_dict[bin]+eff_cut_increment*eff_cut_sign)

                # random sample of cut (upper and lower limits from significance scan)
                bin_range = f'{split}_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}_range'
                eff_cut_range = eff_cut_dict[bin_range]/100 - 0.01
                #print(f"BDT efficiency cut variation range: {eff_cut_range}")
                lower_limit = 0.70 #float(formatted_eff_cut) - eff_cut_range
                upper_limit = 0.90 #float(formatted_eff_cut) + eff_cut_range
                # if float(formatted_eff_cut) > THRESH_EFF:
                #     upper_limit = MAX_EFF
                random_cut = lower_limit + ROOT.gRandom.Rndm()*(upper_limit-lower_limit)
                formatted_eff_cut = "{:.2f}".format(random_cut)
                # print("Formatted eff cut = "+formatted_eff_cut)

                # get signal
                h_raw_yield = signal_extraction_file.Get(f'{bin}_{bkg_shape}/fRawYields;1')
                eff_index = h_raw_yield.FindBin(float(formatted_eff_cut))
                raw_yield = h_raw_yield.GetBinContent(eff_index)
                raw_yield_error = h_raw_yield.GetBinError(eff_index)
                if raw_yield < 0.1:
                    continue

                # get mass
                h_mass_eff = signal_extraction_file.Get(f'{bin}_{bkg_shape}/fMass;1')
                h_mass_mc_eff = signal_extraction_file.Get(f'{bin}_{bkg_shape}/fMassMC;1')
                h_mass_dat_eff = signal_extraction_file.Get(f'{bin}_{bkg_shape}/fMassData;1')
                eff_index = h_mass_eff.FindBin(float(formatted_eff_cut))
                mass = h_mass_eff.GetBinContent(eff_index)
                mass_error = h_mass_eff.GetBinError(eff_index)
                mass_mc = h_mass_mc_eff.GetBinContent(eff_index)
                mass_mc_error = h_mass_mc_eff.GetBinError(eff_index)
                mass_dat = h_mass_dat_eff.GetBinContent(eff_index)
                mass_dat_error = h_mass_dat_eff.GetBinError(eff_index)
                
                # apply corrections
                presel_eff_map = np.logical_and(
                    presel_eff_bin_centers > ct_bins[0],
                    presel_eff_bin_centers < ct_bins[1])
                presel_eff = presel_eff_counts[presel_eff_map]
                bdt_eff = random_cut
                eff = presel_eff * eff_cut_dict[bin]
                #abs = g_abs_correction.GetPointY(CT_BINS_CENT[i_cent_bins].index(ct_bins[0]))

                ct_bin_index = h_corrected_yields[i_split].FindBin(ct_bins[0]+0.05)

                h_corrected_yields[i_split].SetBinContent(ct_bin_index, raw_yield/eff[0])
                h_corrected_yields[i_split].SetBinError(ct_bin_index, raw_yield_error/eff[0])
                h_mass[i_split].SetBinContent(ct_bin_index, mass)
                h_mass[i_split].SetBinError(ct_bin_index, mass_error)
                h_mass_mc[i_split].SetBinContent(ct_bin_index, mass_mc)
                h_mass_mc[i_split].SetBinError(ct_bin_index, mass_mc_error)
                h_mass_dat[i_split].SetBinContent(ct_bin_index, mass_dat)
                h_mass_dat[i_split].SetBinError(ct_bin_index, mass_dat_error)
                # h_m_minus_mc_distribution_reject[i_split].Fill((mass_dat-mass_mc)*1e6)
                h_m_minus_mc[i_split].SetBinContent(ct_bin_index, (mass_dat - mass_mc)*1e6)
                h_m_minus_mc[i_split].SetBinError(ct_bin_index, 1e6*np.sqrt(mass_dat_error*mass_dat_error+mass_mc_error*mass_mc_error))

            mean_delta_mass = h_m_minus_mc_distribution_reject[i_split].GetMean()
            sigma_delta_mass = h_m_minus_mc_distribution_reject[i_split].GetRMS()
            n_entries = h_m_minus_mc_distribution_reject[i_split].GetEntries()
            prob_reject = 1/4/n_entries
            gaus_cdf = ROOT.TF1("gaus_cdf","ROOT::Math::normal_cdf(x)")
            limit = -gaus_cdf.GetX(prob_reject,-5,0)
            # print(F"LIMIT = {limit}")

            # for ct_bins in zip(CT_BINS_CENT[i_cent_bins][:-1], CT_BINS_CENT[i_cent_bins][1:]):
            #     ct_bin_index = h_corrected_yields[i_split].FindBin(ct_bins[0]+0.5)
            #     bin_content = h_m_minus_mc[i_split].GetBinContent(ct_bin_index)
            #     if np.abs(bin_content-mean_delta_mass) > (limit*sigma_delta_mass):
            #         # print("remove point")
            #         h_m_minus_mc[i_split].SetBinContent(ct_bin_index,0)
            #         h_m_minus_mc[i_split].SetBinError(ct_bin_index,0)

            # set labels
            h_corrected_yields[i_split].GetXaxis().SetTitle("#it{c}t (cm)")
            h_corrected_yields[i_split].GetYaxis().SetTitle("d#it{N}/d(#it{c}t) (cm^{-1})")
            h_mass[i_split].GetXaxis().SetTitle("#it{c}t (cm)")
            h_mass[i_split].GetYaxis().SetTitle("#it{m}_{#Lambda} (GeV/#it{c}^{2})")
            h_mass[i_split].Fit("pol0","Q")
            h_mass_mc[i_split].GetXaxis().SetTitle("#it{c}t (cm)")
            h_mass_mc[i_split].GetYaxis().SetTitle("#it{m}_{#Lambda} (GeV/#it{c}^{2})")
            h_mass_dat[i_split].GetXaxis().SetTitle("#it{c}t (cm)")
            h_mass_dat[i_split].GetYaxis().SetTitle("#it{m}_{#Lambda} (GeV/#it{c}^{2})")
            h_m_minus_mc[i_split].GetXaxis().SetTitle("#it{c}t (cm)")
            h_m_minus_mc[i_split].GetYaxis().SetTitle("#it{m}_{#Lambda}^{data}-#it{m}_{#Lambda}^{MC} (keV/#it{c}^{2})")
            res_deltaM = h_m_minus_mc[i_split].Fit("pol0","QS")

            for i_bin in range(len(bins))[2:]:
                bin_width = h_corrected_yields[i_split].GetBinWidth(i_bin)
                #print(f"bin: {h_corrected_yields[i_split].GetBinLowEdge(i_bin)}; bin width: {bin_width}")
                bin_content = h_corrected_yields[i_split].GetBinContent(i_bin)
                bin_error = h_corrected_yields[i_split].GetBinError(i_bin)
                h_corrected_yields[i_split].SetBinContent(i_bin, bin_content/bin_width)
                h_corrected_yields[i_split].SetBinError(i_bin, bin_error/bin_width)

            # fit with exponential pdf
            fit_function_expo = ROOT.TF1("expo", "expo", 2, 35)
            if cent_bins[0] == 30:
                fit_function_expo = ROOT.TF1("expo", "expo", 2, 14)
            elif cent_bins[1] == 90:
                fit_function_expo = ROOT.TF1("expo", "expo", 0, 35)
            res = h_corrected_yields[i_split].Fit(fit_function_expo, "QRMIS+")

            # compute lifetime
            tau = -1/fit_function_expo.GetParameter(1)*100/SPEED_OF_LIGHT # ps

            if h_m_minus_mc[i_split].GetEntries() > 0:
                if (h_m_minus_mc[i_split].GetFunction("pol0").GetChisquare()/h_m_minus_mc[i_split].GetFunction("pol0").GetNDF()) < 3:
                    systematics_file.cd(f'{cent_bins[0]}_{cent_bins[1]}')
                    # h_corrected_yields[i_split].Write()
                    #lifetime_tmp[i_split] = tau
                    h_m_minus_mc_distribution.Fill(h_m_minus_mc[i_split].GetFunction("pol0").GetParameter(0))
                    h_prob_distribution.Fill(res_deltaM.Prob())
                    #h_fit_status.Fill(fit_function_expo.IsValid())
                    i_trial+=1
        
    systematics_file.cd()

    #     if (lifetime_tmp[0] > 0) and (lifetime_tmp[1] > 0):
    #         h_lifetime[0].Fill(lifetime_tmp[0])
    #         h_lifetime[1].Fill(lifetime_tmp[1])
    #         h_asymmetry_distribution.Fill(lifetime_tmp[1]-lifetime_tmp[0])
    #         i_trial+=1
    # systematics_file.cd()
    # h_asymmetry_distribution.GetXaxis().SetTitle("Asymmetry (ps)")
    # h_asymmetry_distribution.GetYaxis().SetTitle("Entries")
    # h_asymmetry_distribution.Write()

    # save plots
    # h_asymmetry_distribution.Rebin(8)
    # c = ROOT.TCanvas("c", "c")
    # ROOT.gStyle.SetOptStat(110001110)
    # c.SetTicks(1, 1)
    # h_asymmetry_distribution.SetDrawOption("histo")
    # h_asymmetry_distribution.SetLineWidth(2)
    # h_asymmetry_distribution.SetFillStyle(3345)
    # h_asymmetry_distribution.SetFillColor(ROOT.kBlue)
    # h_asymmetry_distribution.Draw("histo")
    # c.Print(f"plots/{h_asymmetry_distribution.GetName()}.png")

    h_m_minus_mc_distribution.Write()
    h_m_minus_mc_distribution.GetXaxis().SetTitle("#it{m}_{#Lambda}^{data}-#it{m}_{#Lambda}^{MC} (keV/#it{c}^{2})")
    h_m_minus_mc_distribution.GetYaxis().SetTitle("Entries")
    h_prob_distribution.GetXaxis().SetTitle("Prob")
    h_prob_distribution.GetYaxis().SetTitle("Entries")
    h_prob_distribution.Write()
    # h_lifetime[0].GetXaxis().SetTitle("#tau (ps)")
    # h_lifetime[1].GetXaxis().SetTitle("#tau (ps)")
    # h_lifetime[0].Write()
    # h_lifetime[1].Write()
    h_fit_status.Write()

    del h_asymmetry_distribution
    del h_prob_distribution
    del h_fit_status
    del h_lifetime

systematics_file.Close()
