#!/usr/bin/env python3
import os
import pickle
import warnings
import argparse

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import ROOT
import uproot
import yaml
from helpers import significance_error, ndarray2roo

SPLIT = False
use_crystalball = True
MAX_EFF = 1.00
mass_PDG = 1.115683 # GeV/c^2

ROOT.gInterpreter.ProcessLine("#include \"../utils/RooDSCBShape.h\"")
ROOT.gInterpreter.ProcessLine(".L ../utils/RooDSCBShape.cxx+")

# avoid pandas warning
warnings.simplefilter(action='ignore', category=FutureWarning)
ROOT.gROOT.SetBatch()

parser = argparse.ArgumentParser(prog='signal_extraction', allow_abbrev=True)
parser.add_argument('-bkgExpo', action='store_true')
args = parser.parse_args()

BKG_EXPO = args.bkgExpo

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
CT_BINS = params['CT_BINS']
CT_BINS_CENT = params['CT_BINS_CENT']
CENTRALITY_LIST = params['CENTRALITY_LIST']
RANDOM_STATE = params['RANDOM_STATE']
##################################################################

# split matter/antimatter
SPLIT_LIST = ['all']
if SPLIT:
    SPLIT_LIST = ['antimatter', 'matter']

bkg_shape = 'pol1'
if BKG_EXPO:
    bkg_shape = 'expo'

score_eff_arrays_dict = pickle.load(open("file_score_eff_dict", "rb"))
eff_array = np.arange(0.10, MAX_EFF, 0.01)

for split in SPLIT_LIST:
    for i_cent_bins in range(len(CENTRALITY_LIST)):
        cent_bins = CENTRALITY_LIST[i_cent_bins]
        for ct_bins in zip(CT_BINS_CENT[i_cent_bins][:-1], CT_BINS_CENT[i_cent_bins][1:]):

            # ct_bins_df_index = int(CT_BINS_CENT[i_cent_bins][0]/5 -1) # where 5 is the bin size
            # print(f'Index = {ct_bins_df_index}')
            bin = f'{split}_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}'
            # bin_df = f'{split}_{cent_bins[0]}_{cent_bins[1]}_{CT_BINS[ct_bins_df_index][0]}_{CT_BINS[ct_bins_df_index][1]}'
            df_data = pd.read_parquet(f'df/{bin}.parquet.gzip')
            df_signal = pd.read_parquet(f'df/mc_{bin}')

            # ROOT.Math.MinimizerOptions.SetDefaultTolerance(1e-2)
            root_file_signal_extraction = ROOT.TFile("SignalExtraction.root", "update")
            root_file_signal_extraction.mkdir(f'{bin}_{bkg_shape}')

            # raw yields histogram
            h_raw_yields = ROOT.TH1D("fRawYields", "fRawYields", 101, -0.005, 1.005)

            # mass histogram
            h_mass = ROOT.TH1D("fMass","fMass", 101, -0.005, 1.005)

            # significance histogram
            h_significance = ROOT.TH1D("fSignificance", "fSignificance", 101, -0.005, 1.005)

            for eff_score in zip(eff_array, score_eff_arrays_dict[bin]):
                if (ct_bins[0] > 0.5) and (eff_score[0] < 0.69 or eff_score[0] > 0.99):
                    continue
                formatted_eff = "{:.2f}".format(eff_score[0])
                print(f'processing {bin}: eff = {eff_score[0]:.2f}, score = {eff_score[1]:.2f}...')

                df_data_sel = df_data.query(f'model_output > {eff_score[1]}')
                df_signal_sel = df_signal.query(f'model_output > {eff_score[1]} and y_true == 1')

                # get invariant mass distribution (data and mc)
                roo_m = ROOT.RooRealVar("m", "#it{M} (p + #pi^{-})", 1.09, 1.14, "GeV/#it{c}^{2}")
                roo_data_unbinned = ndarray2roo(np.array(df_data_sel['mass']), roo_m)
                roo_mc_signal = ndarray2roo(np.array(df_signal_sel['mass']), roo_m)
                roo_m.setBins(200)
                roo_data = ROOT.RooDataHist('data','data',ROOT.RooArgSet(roo_m),roo_data_unbinned)
                roo_mc_signal = ROOT.RooDataHist('data','data',ROOT.RooArgSet(roo_m),roo_mc_signal)

                # declare fit model
                # kde
                roo_n_signal = ROOT.RooRealVar('N_{signal}', 'Nsignal', 1., 1.e6)
                delta_mass = ROOT.RooRealVar("#deltam", 'deltaM', -0.004, 0.004, 'GeV/c^{2}')
                shifted_mass = ROOT.RooAddition("mPrime", "m + #Deltam", ROOT.RooArgList(roo_m, delta_mass))
                #roo_signal = ROOT.RooKeysPdf("signal", "signal", shifted_mass, roo_m,
                           #                  roo_mc_signal, ROOT.RooKeysPdf.NoMirror, 2)

                # cb signal
                mass = ROOT.RooRealVar('mass','mass',1.11,1.12)
                sigma_left = ROOT.RooRealVar('sigma_left','sigma_left',0.,0.005)
                sigma_right = ROOT.RooRealVar('sigma_right','sigma_right',0.,0.005)
                alpha_left = ROOT.RooRealVar('alpha_left','alpha_left',0.,2.)
                alpha_right = ROOT.RooRealVar('alpha_right','alpha_right',0.,2.)
                n_left = ROOT.RooRealVar('n_left','n_left',0.,10.)
                n_right = ROOT.RooRealVar('n_right','n_right',0.,10.)
                roo_signal = ROOT.RooDSCBShape('signal','signal',roo_m,mass,sigma_left,alpha_left,n_left,alpha_right,n_right)    
                roo_signal_copy = ROOT.RooDSCBShape(roo_signal)

                # background
                roo_n_background = ROOT.RooRealVar('N_{bkg}', 'Nbackground', 1., 1.e7)
                roo_slope = ROOT.RooRealVar('slope', 'slope', -20., 20.)
                roo_bkg = ROOT.RooRealVar()

                if not BKG_EXPO:
                    roo_bkg = ROOT.RooPolynomial('background', 'background', roo_m, ROOT.RooArgList(roo_slope))
                else:
                    roo_bkg = ROOT.RooExponential('background', 'background', roo_m, roo_slope)

                # model
                roo_model = ROOT.RooAddPdf(
                    'model', 'model', ROOT.RooArgList(roo_signal, roo_bkg),
                    ROOT.RooArgList(roo_n_signal, roo_n_background))

                # fit
                ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)
                ROOT.RooMsgService.instance().setSilentMode(ROOT.kTRUE)
                ROOT.gErrorIgnoreLevel = ROOT.kError
                roo_model.fitTo(roo_data, ROOT.RooFit.Save(), ROOT.RooFit.Extended(ROOT.kTRUE))
                r = roo_model.fitTo(roo_data, ROOT.RooFit.Save(), ROOT.RooFit.Extended(ROOT.kTRUE))

                print(f'fit status: {r.status()}')
                if r.status() == 0: # and delta_mass.getError() > 1.e-6:

                    # plot
                    nBins = 200
                    xframe = roo_m.frame(1.09, 1.14, nBins)
                    xframe.SetTitle(
                        str(ct_bins[0]) + '#leq #it{c}t<' + str(ct_bins[1]) + ' cm, ' + str(cent_bins[0]) + '-' +
                        str(cent_bins[1]) + '%, BDT efficiency = ' + str(formatted_eff))
                    xframe.SetName(f'fInvMass_{formatted_eff}')
                    roo_data.plotOn(xframe, ROOT.RooFit.Name('data'))
                    roo_model.plotOn(
                        xframe, ROOT.RooFit.Components('background'),
                        ROOT.RooFit.Name('background'),
                        ROOT.RooFit.LineStyle(ROOT.kDashed),
                        ROOT.RooFit.LineColor(ROOT.kGreen))
                    roo_model.plotOn(xframe, ROOT.RooFit.Components('signal'), ROOT.RooFit.Name('signal'),
                                     ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kRed))
                    roo_model.plotOn(xframe, ROOT.RooFit.Name('model'), ROOT.RooFit.LineColor(ROOT.kBlue))

                    formatted_chi2 = "{:.2f}".format(xframe.chiSquare('model', 'data'))
                    roo_model.paramOn(xframe, ROOT.RooFit.Label(
                        '#chi^{2}/NDF = '+formatted_chi2),
                        ROOT.RooFit.Layout(0.57, 0.85, 0.88))
                    xframe.getAttText().SetTextFont(44)
                    xframe.getAttText().SetTextSize(20)
                    xframe.getAttLine().SetLineWidth(0)

                    print(f'chi2/NDF: {formatted_chi2}, edm: {r.edm()}')
                    if float(formatted_chi2) < 10: # and r.edm() < 1:

                        # fit mc distribution to get sigma and mass
                        roo_mean_mc = ROOT.RooRealVar("mean", "mean", 1.11, 1.12)
                        roo_sigma_mc = ROOT.RooRealVar("sigma", "sigma", 0.0001, 0.0040)
                        gaus = ROOT.RooGaussian('gaus', 'gaus', roo_m, roo_mean_mc, roo_sigma_mc)
                        gaus.fitTo(roo_mc_signal)

                        # fit mc distribution with dscb
                        mass_mc = ROOT.RooRealVar('mass','mass',1.11,1.12)
                        sigma_left_mc = ROOT.RooRealVar('sigma_left','sigma_left',0.,0.005)
                        sigma_right_mc = ROOT.RooRealVar('sigma_right','sigma_right',0.,0.005)
                        alpha_left_mc = ROOT.RooRealVar('alpha_left','alpha_left',0.,2.)
                        alpha_right_mc = ROOT.RooRealVar('alpha_right','alpha_right',0.,2.)
                        n_left_mc = ROOT.RooRealVar('n_left','n_left',0.,10.)
                        n_right_mc = ROOT.RooRealVar('n_right','n_right',0.,10.)
                        roo_signal_mc = ROOT.RooDSCBShape('signal','signal',roo_m,mass_mc,sigma_left_mc,alpha_left_mc,n_left_mc,alpha_right_mc,n_right_mc) 
                        roo_signal_plot = ROOT.RooDSCBShape(roo_signal_mc)
                        for _ in range(2):
                            roo_signal_mc.fitTo(roo_mc_signal)

                        # mass
                        mass_val = mass_mc.getVal()

                        # significance
                        m_set = ROOT.RooArgSet(roo_m)
                        normSet = ROOT.RooFit.NormSet(m_set)
                        roo_m.setRange(
                            'signalRange', mass_val - 3 * roo_sigma_mc.getVal(),
                            mass_val + 3 * roo_sigma_mc.getVal())
                        signal_int = (roo_model.pdfList().at(0).createIntegral(
                            m_set, normSet, ROOT.RooFit.Range("signalRange"))).getVal()
                        print(f'signal integral = {signal_int}')
                        bkg_int = (roo_model.pdfList().at(1).createIntegral(
                            m_set, normSet, ROOT.RooFit.Range("signalRange"))).getVal()
                        print(f'background integral = {bkg_int}')
                        sig = signal_int*roo_n_signal.getVal()
                        bkg = bkg_int*roo_n_background.getVal()
                        significance_val = sig/np.sqrt(sig+bkg)
                        significance_err = significance_error(sig, bkg)

                        # fill significance histogram
                        eff_index = h_significance.FindBin(float(formatted_eff))
                        h_significance.SetBinContent(eff_index, significance_val)
                        h_significance.SetBinError(eff_index, significance_err)

                        if significance_val > 2.95:
                            # fill raw yields histogram
                            h_raw_yields.SetBinContent(eff_index, roo_n_signal.getVal())
                            h_raw_yields.SetBinError(eff_index, roo_n_signal.getError())

                            # fill mass histogram
                            h_mass.SetBinContent(eff_index, mass.getVal()-mass_mc.getVal()+mass_PDG)
                            h_mass.SetBinError(eff_index, np.sqrt(mass.getError()*mass.getError()+mass_mc.getError()*mass_mc.getError()))

                            # write to file
                            root_file_signal_extraction.cd(f'{bin}_{bkg_shape}')
                            xframe.Write()

                            # draw on canvas and save plots
                            canv = ROOT.TCanvas()
                            canv.cd()
                            text_mass = ROOT.TLatex(
                                1.092, 0.74 * xframe.GetMaximum(),
                                "#it{m}_{#Lambda} = " + "{:.6f}".format(mass_val) + " GeV/#it{c^{2}}")
                            text_mass.SetTextFont(44)
                            text_mass.SetTextSize(20)
                            text_signif = ROOT.TLatex(1.092, 0.91 * xframe.GetMaximum(),
                                                    "S/#sqrt{S+B} (3#sigma) = " + "{:.3f}".format(significance_val) + " #pm " +
                                                    "{:.3f}".format(significance_err))
                            text_signif.SetTextFont(44)
                            text_signif.SetTextSize(20)
                            text_sig = ROOT.TLatex(1.092, 0.84 * xframe.GetMaximum(), "S (3#sigma) = " + "{:.1f}".format(sig) + " #pm " + "{:.1f}".format(signal_int*roo_n_signal.getError()))
                            text_sig.SetTextFont(44)
                            text_sig.SetTextSize(20)
                            text_bkg = ROOT.TLatex(1.092, 0.77 * xframe.GetMaximum(), "B (3#sigma) = " + "{:.1f}".format(bkg) + " #pm" + "{:.1f}".format(bkg_int*roo_n_background.getError()))
                            text_bkg.SetTextFont(44)
                            text_bkg.SetTextSize(20)
                            xframe.Draw("")
                            # text_mass.Draw("same")
                            text_signif.Draw("same")
                            text_sig.Draw("same")
                            text_bkg.Draw("same")
                            print(
                                f'significance = {"{:.3f}".format(significance_val)} +/- {"{:.3f}".format(significance_err)}')
                            if not os.path.isdir('plots/signal_extraction'):
                                os.mkdir('plots/signal_extraction')
                            if not os.path.isdir(f'plots/signal_extraction/{bin}_{bkg_shape}'):
                                os.mkdir(f'plots/signal_extraction/{bin}_{bkg_shape}')
                            canv.Print(f'plots/signal_extraction/{bin}_{bkg_shape}/{eff_score[0]:.2f}_{bin}.pdf')

                            # plot kde and mc
                            frame = roo_m.frame(1.09, 1.14, 130)
                            frame.SetTitle(str(cent_bins[0])+"-"+str(cent_bins[1])+"%, "+str(ct_bins[0])+"#leq #it{c}t<"+str(ct_bins[1])+" cm, BDT efficiency = "+str(formatted_eff))
                            roo_mc_signal.plotOn(frame)
                            roo_signal_plot.plotOn(frame, ROOT.RooFit.Name("DSCB"))
                            gaus.plotOn(frame, ROOT.RooFit.Name("gaussian"), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.LineStyle(ROOT.kDashed))
                            cc = ROOT.TCanvas("cc", "cc")
                            if not os.path.isdir('plots/kde_signal'):
                                os.mkdir('plots/kde_signal')
                            if not os.path.isdir(f'plots/kde_signal/{bin}'):
                                os.mkdir(f'plots/kde_signal/{bin}')
                            frame.Draw()
                            leg_mc = ROOT.TLegend(0.7, 0.8, 0.85, 0.7)
                            leg_mc.AddEntry(frame.findObject("DSCB"), "DSCB")
                            leg_mc.AddEntry(frame.findObject("gaussian"), "Gaussian")
                            leg_mc.SetTextFont(44)
                            leg_mc.SetTextSize(20)
                            leg_mc.SetBorderSize(0)
                            leg_mc.Draw("same")
                            cc.SetLogy(ROOT.kTRUE)
                            cc.Write()
                            cc.Print(f'plots/kde_signal/{bin}/{formatted_eff}_{bin}.pdf')

            h_raw_yields.GetXaxis().SetTitle("BDT efficiency")
            h_raw_yields.GetYaxis().SetTitle("#it{N_{raw}}")
            h_raw_yields.Write()

            h_mass.GetXaxis().SetTitle("BDT efficiency")
            h_mass.GetYaxis().SetTitle("#it{m} (GeV/#it{c^2})")
            h_mass.Write()

            h_significance.GetXaxis().SetTitle("BDT efficiency")
            h_significance.GetYaxis().SetTitle("S / #sqrt{S + B}")
            h_significance.Write()

            root_file_signal_extraction.Close()
