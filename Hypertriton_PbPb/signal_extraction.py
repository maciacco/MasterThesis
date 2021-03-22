#!/usr/bin/env python3
import os
import pickle
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import ROOT
import uproot
import xgboost as xgb
import yaml
from hipe4ml.model_handler import ModelHandler
from hipe4ml.tree_handler import TreeHandler


def ndarray2roo(ndarray, var):
    if isinstance(ndarray, ROOT.RooDataSet):
        print('Already a RooDataSet')
        return ndarray

    assert isinstance(ndarray, np.ndarray), 'Did not receive NumPy array'
    assert len(ndarray.shape) == 1, 'Can only handle 1d array'

    name = var.GetName()
    x = np.zeros(1, dtype=np.float64)

    tree = ROOT.TTree('tree', 'tree')
    tree.Branch(f'{name}', x, f'{name}/D')

    for i in ndarray:
        x[0] = i
        tree.Fill()

    array_roo = ROOT.RooDataSet(
        'data', 'dataset from tree', tree, ROOT.RooArgSet(var))
    return array_roo


SPLIT = True

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

DATA_PATH = params['DATA_PATH']
CT_BINS = params['CT_BINS']
CENTRALITY_LIST = params['CENTRALITY_LIST']
RANDOM_STATE = params['RANDOM_STATE']
##################################################################

# split matter/antimatter
SPLIT_LIST = ['']
if SPLIT:
    SPLIT_LIST = ['antimatter', 'matter']

score_eff_arrays_dict = pickle.load(open("file_score_eff_dict", "rb"))
eff_array = np.arange(0.10, 0.51, 0.01)

for split in SPLIT_LIST:
    for cent_bins in CENTRALITY_LIST:
        for ct_bins in zip(CT_BINS[:-1], CT_BINS[1:]):

            bin = f'{split}_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}'
            df_data = pd.read_parquet(f'df/{bin}')

            for eff_score in zip(eff_array, score_eff_arrays_dict[bin]):
                formatted_eff = "{:.2f}".format(1-eff_score[0])
                print(f'eff = {1-eff_score[0]:.2f}, score = {eff_score[1]:.2f}')
                df_data_sel = df_data.query(f'model_output > {eff_score[1]}')

                # get invariant mass distribution
                root_file_signal_extraction = ROOT.TFile("SignalExtraction.root", "update")
                root_file_signal_extraction.mkdir(f'Efficiency_{1-eff_score[0]:.2f}')
                roo_m = ROOT.RooRealVar("m", "#it{M} (^{3}He + #pi^{-})", 2.96, 3.04, "GeV/#it{c}^{2}")
                roo_data = ndarray2roo(np.array(df_data_sel['m']), roo_m)

                # declare fit model (gaus + pol2)
                roo_n_signal = ROOT.RooRealVar('Nsignal', 'N_{signal}', 1., 50.)
                roo_mean = ROOT.RooRealVar('mean', '#mu', 2.991, 2.9880, 2.9950)
                roo_sigma = ROOT.RooRealVar('sigma', '#sigma', 0.0010, 0.0060)
                roo_signal = ROOT.RooGaussian('signal', 'signal', roo_m, roo_mean, roo_sigma)

                roo_n_background = ROOT.RooRealVar('Nbackground', 'N_{bkg}', 1500., 1.e3, 1.e4)
                roo_a = ROOT.RooRealVar('a', 'a', -0.01, 0.01)
                roo_b = ROOT.RooRealVar('b', 'b', -.5, -1.e-5)
                roo_bkg = ROOT.RooPolynomial('background', 'background', roo_m, ROOT.RooArgList(roo_b, roo_a))

                roo_model = ROOT.RooAddPdf(
                    'model', 'model', ROOT.RooArgList(roo_signal, roo_bkg),
                    ROOT.RooArgList(roo_n_signal, roo_n_background))

                # fit
                ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)
                ROOT.RooMsgService.instance().setSilentMode(ROOT.kTRUE)
                ROOT.gErrorIgnoreLevel = ROOT.kError
                roo_bkg.fitTo(roo_data, ROOT.RooFit.Save())

                # plot
                xframe = roo_m.frame(2.96, 3.04, 32)
                xframe.SetTitle(
                    str(ct_bins[0]) + '#leq #it{c}t<' + str(ct_bins[1]) + ' cm, ' + str(cent_bins[0]) + '-' +
                    str(cent_bins[1]) + '%, ' + str(formatted_eff))
                xframe.SetName(f'fInvMass_{bin}')
                roo_data.plotOn(xframe, ROOT.RooFit.Name('data'))
                roo_model.plotOn(
                    xframe, ROOT.RooFit.Components('background'),
                    ROOT.RooFit.Name('background'),
                    ROOT.RooFit.LineStyle(ROOT.kDashed),
                    ROOT.RooFit.LineColor(ROOT.kGreen))
                roo_model.plotOn(xframe, ROOT.RooFit.Components('signal'), ROOT.RooFit.Name('signal'),
                                 ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kRed))
                roo_model.plotOn(xframe, ROOT.RooFit.Name('model'), ROOT.RooFit.LineColor(ROOT.kBlue))
                roo_model.paramOn(xframe, ROOT.RooFit.Label(
                    '#chi^{2}/NDF = '+str(xframe.chiSquare('model', 'data'))),
                    ROOT.RooFit.Layout(0.68, 0.96, 0.96))

                # write to file
                root_file_signal_extraction.cd(f'Efficiency_{1-eff_score[0]:.2f}')
                xframe.Write()

                # save plots
                canv = ROOT.TCanvas()
                canv.cd()
                xframe.Draw("")
                if not os.path.isdir('plots/signal_extraction'):
                    os.mkdir('plots/signal_extraction')
                canv.Print(f'plots/signal_extraction/{1-eff_score[0]:.2f}_{bin}.png')

                root_file_signal_extraction.Close()
