#!/usr/bin/env python3

import ROOT
import numpy as np

n_pt_bins = 8
pt_bins = [0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90]
cent_bins = [[0,5],[5,10],[30,50]]
pt_bins_array = np.asarray(pt_bins)

input_file = ROOT.TFile.Open("out/SignalProtonGausDExpSignal1_LongMCTracks_1.root")
input_file_eff = ROOT.TFile.Open("out/EfficiencyProtonMC_21l5_false__.root")
out_file = ROOT.TFile.Open("CheckTPCSignalShape.root","recreate")

for cent_bin in cent_bins:
    h_fraction = ROOT.TH1D(f"hFraction_{cent_bin[0]}_{cent_bin[1]}",";#it{p}_{T} (GeV/#it{c});fraction",n_pt_bins,pt_bins_array)
    h_eff_a = input_file_eff.Get(f"_/fAEff_TPC_{cent_bin[0]}_{cent_bin[1]}")
    h_eff_m = input_file_eff.Get(f"_/fMEff_TPC_{cent_bin[0]}_{cent_bin[1]}")
    for pt_bin in zip(pt_bins[:-1],pt_bins[1:]):
        format_pt_bin_0 = "{:.2f}".format(pt_bin[0])
        format_pt_bin_1 = "{:.2f}".format(pt_bin[1])
        h_anti = input_file.Get(f"_1_1_1/fATPCSignal_{cent_bin[0]}_{cent_bin[1]}_{format_pt_bin_0}_{format_pt_bin_1}")
        h_matt = input_file.Get(f"_1_1_1/fMTPCSignal_{cent_bin[0]}_{cent_bin[1]}_{format_pt_bin_0}_{format_pt_bin_1}")
        anti_counts_3_sigma = h_anti.Integral(h_anti.FindBin(-3.),h_anti.FindBin(2.999))
        matt_counts_3_sigma = h_matt.Integral(h_anti.FindBin(-3.),h_anti.FindBin(2.999))
        eff_a = h_eff_a.GetBinContent(h_eff_a.FindBin(pt_bin[0]+0.001))
        eff_m = h_eff_m.GetBinContent(h_eff_m.FindBin(pt_bin[0]+0.001))
        counts_diff = matt_counts_3_sigma/eff_m-anti_counts_3_sigma/eff_a
        tot = matt_counts_3_sigma/eff_m
        h_fraction.SetBinContent(h_fraction.FindBin(pt_bin[0]+0.001),counts_diff/tot)
    h_fraction.Write()

out_file.Close()