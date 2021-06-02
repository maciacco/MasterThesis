#!/usr/bin/env python3
import ROOT

path_he3 = './He3_PbPb/out'
path_hyp = './Hypertriton_PbPb'
centrality_classes = [[0, 5], [5, 10], [30, 50]]
centrality_colors = [ROOT.kOrange+7, ROOT.kAzure+7, ROOT.kTeal+2]

ROOT.gStyle.SetPadTickY(1)
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetTextFont(4)

file_he3 = ROOT.TFile.Open(path_he3 + '/SpectraHe3.root')
file_hyp = ROOT.TFile.Open(path_hyp + '/Ratio.root')
file_he3_syst = ROOT.TFile.Open(path_he3 + '/SystematicsAll.root')
file_hyp_syst = ROOT.TFile.Open(path_hyp + '/Systematics.root')

file_out = ROOT.TFile.Open('FinalPlot.root', 'recreate')

for i_cent, cent in enumerate(centrality_classes):

    # get histograms
    ratio_he3 = file_he3.Get(f'1.0_89_1_1_1/fRatio_{cent[0]}_{cent[1]}')
    ratio_hyp = file_hyp.Get(f'fRatio_{cent[0]}_{cent[1]}')
    ratio_he3_distribution = file_he3_syst.Get(f'hist/fFitPar_{cent[0]}_{cent[1]}')
    ratio_hyp_distribution = file_hyp_syst.Get(f'fParameterDistribution_{cent[0]}_{cent[1]}')

    # get fit functions
    fit_he3 = ratio_he3.GetFunction("pol0")
    fit_hyp = ratio_hyp.GetFunction("pol0")

    # get ratios and errors
    ratio_he3 = fit_he3.GetParameter(0)
    ratio_he3_err = fit_he3.GetParError(0)

    ratio_hyp = fit_hyp.GetParameter(0)
    ratio_hyp_err = fit_hyp.GetParError(0)

    # systematic error
    syst_he3 = ratio_he3_distribution.GetRMS()
    syst_hyp = ratio_hyp_distribution.GetRMS()

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
    a.SetLabelSize(0.062)

    # set histogram style
    ratios_vs_b.SetTitle(f'{cent[0]}-{cent[1]}%')
    ratios_vs_b.SetMarkerColor(centrality_colors[i_cent])
    ratios_vs_b.SetLineColor(centrality_colors[i_cent])
    ratios_vs_b.SetMarkerStyle(20)
    ratios_vs_b.SetMarkerSize(0.8)
    ratios_vs_b.GetYaxis().SetRangeUser(0., 1.2)

    ratios_vs_b_graph = ROOT.TGraphErrors(ratios_vs_b)
    for _ in range(8):
        ratios_vs_b_graph.RemovePoint(0)
    
    # set systematic errors
    ratios_vs_b_graph.SetPointError(0, 0.5, syst_hyp)
    ratios_vs_b_graph.SetPointError(1, 0.5, syst_he3)
    
    ratios_vs_b_graph.SetFillColor(ROOT.kWhite)
    ratios_vs_b_graph.SetFillStyle(3000)
    ratios_vs_b_graph.SetLineColor(centrality_colors[i_cent])
    ratios_vs_b_graph.SetDrawOption("AP5")

    # write to file
    ratios_vs_b.Write()
    
    # save plot
    c = ROOT.TCanvas("c", "c")
    c.cd()
    ratios_vs_b.Draw("pe")
    ratios_vs_b_graph.Draw("P5 same")
    c.Print(f"Ratios_{cent[0]}_{cent[1]}.png")

file_out.Close()
