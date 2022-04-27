#!/usr/bin/env python3
import os
import pickle
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import norm
import ROOT
import uproot
import yaml

ROOT.gInterpreter.ProcessLine('#include "../utils/Utils.h"')
ROOT.gInterpreter.ProcessLine('using namespace utils;')

m_PDG = 1.115683

SPEED_OF_LIGHT = 2.99792458
SPLIT = False

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

if not os.path.isdir("plots/cpt_and_lifetime"):
    os.mkdir("plots/cpt_and_lifetime")

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)

eff_cut_dict = pickle.load(open("file_eff_cut_dict", "rb"))
presel_eff_file = uproot.open('PreselEff.root')
analysis_results_file = uproot.open(ANALYSIS_RESULTS_PATH)
signal_extraction_file = ROOT.TFile.Open('SignalExtraction.root')
signal_extraction_keys = uproot.open('SignalExtraction.root').keys()

# get centrality selected histogram
cent_counts, cent_edges = analysis_results_file['Centrality_selected;1'].to_numpy()
cent_bin_centers = (cent_edges[:-1]+cent_edges[1:])/2

ratio_file = ROOT.TFile.Open('Mass.root', 'recreate')

for i_cent_bins in range(len(CENTRALITY_LIST)):
    cent_bins = CENTRALITY_LIST[i_cent_bins]

    # get number of events
    cent_range_map = np.logical_and(cent_bin_centers > cent_bins[0], cent_bin_centers < cent_bins[1])
    counts_cent_range = cent_counts[cent_range_map]
    evts = np.sum(counts_cent_range)
    print(f'Number of events: {evts}')

    h_corrected_yields = [ROOT.TH1D(), ROOT.TH1D()]
    h_mass = [ROOT.TH1D(), ROOT.TH1D()]
    h_mass_mc = [ROOT.TH1D(), ROOT.TH1D()]
    h_mass_dat = [ROOT.TH1D(), ROOT.TH1D()]
    h_m_minus_mc = [ROOT.TH1D(), ROOT.TH1D()]
    h_m_minus_mc_distribution = [ROOT.TH1D(), ROOT.TH1D()] # distribution to reject outliers
    for i_split, split in enumerate(SPLIT_LIST):
        print(f'{i_split} -> {split}')
        # get preselection efficiency and abs correction histograms
        presel_eff_counts, presel_eff_edges = presel_eff_file[
            f'fPreselEff_vs_ct_{split}_{cent_bins[0]}_{cent_bins[1]};1'].to_numpy()
        presel_eff_bin_centers = (presel_eff_edges[1:]+presel_eff_edges[:-1])/2

        # list of corrected yields
        ct_bins_tmp = [0]
        ct_bins_tmp += CT_BINS_CENT[i_cent_bins]
        if cent_bins[0]==30:
            ct_bins_tmp = [0, 2, 4, 7, 14, 35]
        #if cent_bins[1] == 90:
        #    ct_bin_tmp = CT_BINS_CENT[i_cent_bins]
        bins = np.array(ct_bins_tmp, dtype=float)
        # print(bins)
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
        h_m_minus_mc_distribution[i_split] = ROOT.TH1D(
            f'fDeltaMassLambdaDistribution_{split}_{cent_bins[0]}_{cent_bins[1]}', f'{split}, {cent_bins[0]}-{cent_bins[1]}%', 800, -400, 400)

        for ct_bins in zip(CT_BINS_CENT[i_cent_bins][:-1], CT_BINS_CENT[i_cent_bins][1:]):

            bin = f'{split}_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}'
            formatted_eff_cut = "0.80"# "{:.2f}".format(eff_cut_dict[bin])

            # look for plot with eff = eff_cut (or the nearest one)
            bkg_shape = 'pol1'
            eff_cut_increment = 0
            eff_cut_sign = -1
            while signal_extraction_keys.count(f'{bin}_{bkg_shape}/fInvMass_{formatted_eff_cut};1') == 0:
                if eff_cut_sign == -1:
                    eff_cut_increment += 0.01
                eff_cut_sign *= -1
                formatted_eff_cut = "{:.2f}".format(eff_cut_dict[bin]+eff_cut_increment*eff_cut_sign)

            # get signal
            h_raw_yield = signal_extraction_file.Get(f'{bin}_{bkg_shape}/fRawYields;1')
            eff_index = h_raw_yield.FindBin(float(formatted_eff_cut))
            raw_yield = h_raw_yield.GetBinContent(eff_index)
            raw_yield_error = h_raw_yield.GetBinError(eff_index)
            print(f'eff_cut = {formatted_eff_cut}, raw_yield = {raw_yield}+{raw_yield_error}')

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
            print(f'eff_cut = {formatted_eff_cut}, mass = {mass}+{mass_error}')

            # apply corrections
            # 1. efficiency (presel x BDT)
            presel_eff_map = np.logical_and(
                presel_eff_bin_centers > ct_bins[0],
                presel_eff_bin_centers < ct_bins[1])
            presel_eff = presel_eff_counts[presel_eff_map]
            bdt_eff = float(formatted_eff_cut)
            print(f'bin: {presel_eff_map}, presel_eff: {presel_eff}')
            eff = presel_eff * eff_cut_dict[bin]
            # 2. absorption correction
            # abs = g_abs_correction.GetPointY(CT_BINS_CENT[i_cent_bins].index(ct_bins[0]))
            # print(f"absorption correction for point {CT_BINS_CENT[i_cent_bins].index(ct_bins[0])}: {abs}")
            # 3. efficiency correction
            # eff_correct = eff_abs_correction.GetBinContent(eff_abs_correction.FindBin(ct_bins[0]))

            ct_bin_index = h_corrected_yields[i_split].FindBin(ct_bins[0]+0.05)
            print(f"ct_index = {ct_bin_index}")
            h_corrected_yields[i_split].SetBinContent(ct_bin_index, raw_yield/eff[0])
            h_corrected_yields[i_split].SetBinError(ct_bin_index, raw_yield_error/eff[0])
            h_mass[i_split].SetBinContent(ct_bin_index, mass)
            h_mass[i_split].SetBinError(ct_bin_index, mass_error)
            h_mass_mc[i_split].SetBinContent(ct_bin_index, (mass_mc-m_PDG)*1e6)
            h_mass_mc[i_split].SetBinError(ct_bin_index, mass_mc_error*1e6)
            h_mass_dat[i_split].SetBinContent(ct_bin_index, (mass_dat-m_PDG)*1e6)
            h_mass_dat[i_split].SetBinError(ct_bin_index, mass_dat_error*1e6)
            h_m_minus_mc[i_split].SetBinContent(ct_bin_index, (mass_dat - mass_mc)*1e6)
            h_m_minus_mc[i_split].SetBinError(ct_bin_index, 1e6*np.sqrt(mass_dat_error*mass_dat_error+mass_mc_error*mass_mc_error))
            # h_m_minus_mc_distribution[i_split].Fill((mass_dat-mass_mc)*1e6)

        # mean_delta_mass = h_m_minus_mc_distribution[i_split].GetMean()
        # sigma_delta_mass = h_m_minus_mc_distribution[i_split].GetRMS()
        # n_entries = h_m_minus_mc_distribution[i_split].GetEntries()
        # prob_reject = 1/4/n_entries
        # gaus_cdf = ROOT.TF1("gaus_cdf","ROOT::Math::normal_cdf(x)")
        # limit = -gaus_cdf.GetX(prob_reject,-5,0)
        # print(F"LIMIT = {limit}")

        # for ct_bins in zip(CT_BINS_CENT[i_cent_bins][:-1], CT_BINS_CENT[i_cent_bins][1:]):
        #     ct_bin_index = h_corrected_yields[i_split].FindBin(ct_bins[0]+0.5)
        #     bin_content = h_m_minus_mc[i_split].GetBinContent(ct_bin_index)
        #     if np.abs(bin_content-mean_delta_mass) > (limit*sigma_delta_mass):
        #         print("remove point")
        #         h_m_minus_mc[i_split].SetBinContent(ct_bin_index,0)
        #         h_m_minus_mc[i_split].SetBinError(ct_bin_index,0)

        # set labels
        h_corrected_yields[i_split].GetXaxis().SetTitle("#it{c}t (cm)")
        h_corrected_yields[i_split].GetYaxis().SetTitle("d#it{N}/d(#it{c}t) (cm^{-1})")
        h_mass[i_split].GetXaxis().SetTitle("#it{c}t (cm)")
        h_mass[i_split].GetYaxis().SetTitle("#it{m}_{#Lambda} (GeV/#it{c}^{2})")
        h_mass[i_split].Fit("pol0")
        h_mass_mc[i_split].GetXaxis().SetTitle("#it{c}t (cm)")
        h_mass_mc[i_split].GetYaxis().SetTitle("#it{m}_{#Lambda}^{MC} - #it{m}_{PDG} (keV/#it{c}^{2})")
        h_mass_dat[i_split].GetXaxis().SetTitle("#it{c}t (cm)")
        h_mass_dat[i_split].GetYaxis().SetTitle("#it{m}_{#Lambda}^{data} - #it{m}_{PDG} (keV/#it{c}^{2})")
        h_m_minus_mc[i_split].GetXaxis().SetTitle("#it{c}t (cm)")
        h_m_minus_mc[i_split].GetYaxis().SetTitle("#it{m}_{#Lambda}^{data} - #it{m}_{#Lambda}^{MC} (keV/#it{c}^{2})")
        h_m_minus_mc[i_split].Fit("pol0")
        prob = h_m_minus_mc[i_split].GetFunction("pol0").GetProb()
        gaus = ROOT.TF1("gausn","gausn")
        gaus.SetParameter(0,1)
        gaus.SetParameter(1,0)
        gaus.SetParameter(2,1)
        integral = gaus.Integral(-2.64,2.64)
        #print(f"Probability = {prob}, integral = {(integral)}, n_entries_hist = {h_m_minus_mc_distribution[i_split].GetEntries()}")
        
        mass_fit = h_m_minus_mc[i_split].GetFunction("pol0").GetParameter(0)
        for ct_bins in zip(CT_BINS_CENT[i_cent_bins][:-1], CT_BINS_CENT[i_cent_bins][1:]):
            ct_bin_index = h_corrected_yields[i_split].FindBin(ct_bins[0]+0.05)
            bin_content = h_m_minus_mc[i_split].GetBinContent(ct_bin_index)
            bin_error = h_m_minus_mc[i_split].GetBinError(ct_bin_index)
            h_m_minus_mc_distribution[i_split].Fill((bin_content-mass_fit)/bin_error)
        # h_m_minus_mc_distribution[i_split].Fit("gaus","L")

        #h_corrected_yields[i_split].Scale(1, "width")
        for i_bin in range(len(bins))[2:]:
            bin_width = h_corrected_yields[i_split].GetBinWidth(i_bin)
            # print(f"bin: {h_corrected_yields[i_split].GetBinLowEdge(i_bin)}; bin width: {bin_width}")
            bin_content = h_corrected_yields[i_split].GetBinContent(i_bin)
            bin_error = h_corrected_yields[i_split].GetBinError(i_bin)
            h_corrected_yields[i_split].SetBinContent(i_bin, bin_content/bin_width)
            h_corrected_yields[i_split].SetBinError(i_bin, bin_error/bin_width)
        #h_corrected_yields[i_split].GetYaxis().SetRangeUser(1., 450.)
        h_corrected_yields[i_split].SetMarkerStyle(20)
        h_corrected_yields[i_split].SetMarkerSize(0.8)
        h_corrected_yields[i_split].Write()
        h_mass[i_split].SetMarkerStyle(20)
        h_mass[i_split].GetXaxis().SetRangeUser(5,35)
        h_mass[i_split].GetYaxis().SetRangeUser(1.115,1.117)
        h_mass[i_split].SetMarkerSize(0.8)
        h_mass[i_split].Write()
        h_mass_mc[i_split].SetMarkerStyle(20)
        h_mass_mc[i_split].GetXaxis().SetRangeUser(5,35)
        h_mass_mc[i_split].GetYaxis().SetRangeUser(0,1000)
        h_mass_mc[i_split].SetMarkerSize(0.8)
        h_mass_mc[i_split].Write()
        h_mass_dat[i_split].SetMarkerStyle(20)
        h_mass_dat[i_split].GetXaxis().SetRangeUser(5,35)
        h_mass_dat[i_split].GetYaxis().SetRangeUser(0,1000)
        h_mass_dat[i_split].SetMarkerSize(0.8)
        h_mass_dat[i_split].Write()
        h_m_minus_mc[i_split].SetMarkerStyle(20)
        h_m_minus_mc[i_split].GetXaxis().SetRangeUser(5,35)
        h_m_minus_mc[i_split].GetYaxis().SetRangeUser(-5e2,5e2)
        h_m_minus_mc[i_split].SetMarkerSize(0.8)
        h_m_minus_mc[i_split].Write()
        h_m_minus_mc_distribution[i_split].Write()

        # mass shift chi2 and value
        canv_delta_mass = ROOT.TCanvas("fDeltaMassCanvas","canv_delta_mass")
        h_m_minus_mc[i_split].Draw()
        delta_mass_kev = h_m_minus_mc[i_split].GetFunction("pol0").GetParameter(0)
        delta_mass_error_kev = h_m_minus_mc[i_split].GetFunction("pol0").GetParError(0)
        formatted_delta_mass = "{:.1f}".format(delta_mass_kev)
        formatted_delta_mass_error = "{:.1f}".format(delta_mass_error_kev)
        delta_mass_text = ROOT.TLatex(20, 300, "#delta#it{m} = "+formatted_delta_mass+" #pm "+formatted_delta_mass_error+" keV/#it{c}^{2}")
        delta_mass_text.SetTextFont(44)
        delta_mass_text.SetTextSize(28)
        delta_mass_text.Draw("same")
        formatted_chi2_delta_mass = "{:.2f}".format(h_m_minus_mc[i_split].GetFunction("pol0").GetChisquare())
        chi2_text_delta_mass = ROOT.TLatex(20, 400, "#chi^{2}/NDF = "+formatted_chi2_delta_mass+"/"+str(h_m_minus_mc[i_split].GetFunction("pol0").GetNDF()))
        chi2_text_delta_mass.SetTextFont(44)
        chi2_text_delta_mass.SetTextSize(28)
        chi2_text_delta_mass.Draw("same")
        canv_delta_mass.Write()

        # fit with exponential pdf
        fit_function_expo = ROOT.TF1("expo", "expo", 2, 35)
        if cent_bins[0] == 30:
            fit_function_expo = ROOT.TF1("expo", "expo", 2, 14)
        elif cent_bins[1] == 90:
            fit_function_expo = ROOT.TF1("expo", "expo", 0, 35)
        res_expo = h_corrected_yields[i_split].Fit(fit_function_expo, "RMIS+")

        # compute and display lifetime
        tau = -1/fit_function_expo.GetParameter(1)*100/SPEED_OF_LIGHT # ps
        tau_error = fit_function_expo.GetParError(1)*100/SPEED_OF_LIGHT/fit_function_expo.GetParameter(1)/fit_function_expo.GetParameter(1) # ps
        tau_text = ROOT.TLatex(20, 100, '#tau = ' + "{:.0f}".format(tau) + '#pm' + "{:.0f}".format(tau_error) + ' ps')
        tau_text.SetTextSize(0.05)
        
        # display chi2
        formatted_chi2_lifetime = "{:.2f}".format(fit_function_expo.GetChisquare())
        chi2_lifetime_text = ROOT.TLatex(20, 70, "#chi^{2}/NDF = "+formatted_chi2_lifetime+"/"+str(fit_function_expo.GetNDF()))
        chi2_lifetime_text.SetTextSize(0.05)

        # compute yield
        integral = float()
        integral_error = float()
        integral = fit_function_expo.Integral(0, 1.e9)
        integral_error = fit_function_expo.IntegralError(0, 1.e9, res_expo.GetParams(), res_expo.GetCovarianceMatrix().GetMatrixArray())
        formatted_integral = "{:.10f}".format(integral/2./evts)
        formatted_integral_error = "{:.10f}".format(integral_error/evts/2.)
        # print(f"N_ev = {evts}")
        integral_yield = ROOT.TLatex(20, 40, "1/#it{N_{ev}} #times BR #times d#it{N}/d#it{y} = "+formatted_integral+" #pm "+formatted_integral_error)
        integral_yield.SetTextSize(0.05)

        # draw on canvas
        canv = ROOT.TCanvas()
        ROOT.gStyle.SetOptStat(0)
        canv.SetTicks(1, 1)
        canv.SetName(h_corrected_yields[i_split].GetName())
        canv.SetTitle(h_corrected_yields[i_split].GetTitle())
        canv.cd()
        h_corrected_yields[i_split].Draw("")
        tau_text.Draw("same")
        chi2_lifetime_text.Draw("same")
        integral_yield.Draw("same")
        canv.SetLogy()
        canv.Write() # write to file
        canv.Print(f'plots/cpt_and_lifetime/{h_corrected_yields[i_split].GetName()}.pdf')

    if SPLIT:
        # ratios
        ROOT.gStyle.SetOptStat(0)
        h_ratio = ROOT.TH1D(h_corrected_yields[0])
        h_ratio.SetName(f'fRatio_{cent_bins[0]}_{cent_bins[1]}')
        h_ratio.SetTitle(f'{cent_bins[0]}-{cent_bins[1]}%')
        h_ratio.Divide(h_corrected_yields[0], h_corrected_yields[1], 1, 1)
        h_ratio.GetYaxis().SetTitle("Ratio ^{3}_{#bar{#Lambda}}#overline{H} / ^{3}_{#Lambda}H")
        h_ratio.GetYaxis().SetRangeUser(0., 1.8)
        h_ratio.GetXaxis().SetRangeUser(0., 35.)
        h_ratio.SetMarkerStyle(20)
        h_ratio.SetMarkerSize(0.8)
        h_ratio.Fit("pol0")
        h_ratio.Write()

        # plot ratios
        c = ROOT.TCanvas("c", "c")
        c.SetTicks(1, 1)
        c.cd()
        h_ratio.Draw()
        formatted_ratio = "{:.2f}".format(h_ratio.GetFunction("pol0").GetParameter(0))
        formatted_ratio_error = "{:.2f}".format(h_ratio.GetFunction("pol0").GetParError(0))
        text_x_position = 20
        ratio_text = ROOT.TLatex(text_x_position, 1.6, f"R = {formatted_ratio} #pm {formatted_ratio_error}")
        ratio_text.SetTextFont(44)
        ratio_text.SetTextSize(28)
        ratio_text.Draw("same")
        formatted_chi2 = "{:.2f}".format(h_ratio.GetFunction("pol0").GetChisquare())
        chi2_text = ROOT.TLatex(text_x_position, 1.45, "#chi^{2}/NDF = "+formatted_chi2+"/"+str(h_ratio.GetFunction("pol0").GetNDF()))
        chi2_text.SetTextFont(44)
        chi2_text.SetTextSize(28)
        chi2_text.Draw("same")
        c.Print(f"plots/{h_ratio.GetName()}.pdf")

        del h_corrected_yields
        del h_ratio

ratio_file.Close()
