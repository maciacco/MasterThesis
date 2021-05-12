#!/usr/bin/env python3
import ROOT

path_he3 = './He3_PbPb/out'
path_hyp = './Hypertriton_PbPb'
centrality_classes = [[0, 5], [5, 10], [30, 50]]

file_he3 = ROOT.TFile.Open(path_he3 + '/SpectraHe3.root')
file_hyp = ROOT.TFile.Open(path_hyp + '/Ratio.root')

file_out = ROOT.TFile.Open('FinalPlot.root', 'recreate')

for cent in centrality_classes:

    # get histograms
    ratio_he3 = file_he3.Get(f'1.0_89_1_1_1/fRatio_{cent[0]}_{cent[1]}')
    ratio_hyp = file_hyp.Get(f'fRatio_{cent[0]}_{cent[1]}')

    # get fit functions
    fit_he3 = ratio_he3.GetFunction("pol0")
    fit_hyp = ratio_hyp.GetFunction("pol0")

    # get ratios and errors
    ratio_he3 = fit_he3.GetParameter(0)
    ratio_he3_err = fit_he3.GetParError(0)

    ratio_hyp = fit_hyp.GetParameter(0)
    ratio_hyp_err = fit_hyp.GetParError(0)

    # final plot
    ROOT.gStyle.SetOptStat(0)
    ratios_vs_b = ROOT.TH1D(f'fRatio_vs_b_{cent[0]}_{cent[1]}',';B+S/3; Antimatter / Matter', 10, -0.5, 9.5)
    ratios_vs_b.SetBinContent(10, ratio_he3)
    ratios_vs_b.SetBinError(10, ratio_he3_err)
    ratios_vs_b.SetBinContent(9, ratio_hyp)
    ratios_vs_b.SetBinError(9, ratio_hyp_err)

    # set labels
    a = ratios_vs_b.GetXaxis()
    a.SetBit(ROOT.TAxis.kLabelsHori)
    a.SetBinLabel(1,'0')
    a.SetBinLabel(4,'1')
    a.SetBinLabel(7,'2')
    a.SetBinLabel(10,'3')
    a.SetLabelSize(0.05)

    # write to file
    ratios_vs_b.Write()

file_out.Close()
