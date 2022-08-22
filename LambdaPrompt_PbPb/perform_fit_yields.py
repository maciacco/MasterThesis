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
from hipe4ml import analysis_utils, plot_utils

SPLIT = True
MAX_EFF = 1

N_BINS_INV_MASS = 100
N_BINS_OUTPUT_SCORE = 20
SHM_FRAC = [0.505269,0.505189,0.505109,0.505031,0.504954,0.504877,0.504802,0.504727,0.504653,0.50458,0.504508,0.504436,0.504365,0.504295,0.504225,0.504156,0.504088,0.50402,0.503953,0.503886,0.50382]
ROOT.gInterpreter.ProcessLine("#include \"../utils/RooDSCBShape.h\"")
ROOT.gInterpreter.ProcessLine(".L ../utils/RooDSCBShape.cxx+")
warnings.simplefilter(action='ignore', category=FutureWarning)
ROOT.EnableImplicitMT(40)
ROOT.gROOT.SetBatch()

ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)
ROOT.RooMsgService.instance().setSilentMode(ROOT.kTRUE)
ROOT.gErrorIgnoreLevel = ROOT.kError

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

score_eff_dict = pickle.load(open('file_score_eff_dict','rb'))
eff_array = np.arange(0.10, MAX_EFF, 0.01)

if not os.path.isdir("plots/signal_extraction_"):
    os.mkdir("plots/signal_extraction_")

#df_data = uproot.open(os.path.expandvars(f"/data/mciacco/LambdaPrompt_PbPb/AnalysisResults_data_filtered_dcaV0PV.root"))['LambdaTree'].arrays(library="pd")

for i_cent_bins in range(len(CENTRALITY_LIST)):

    cent_bins = CENTRALITY_LIST[i_cent_bins]
    if not os.path.isdir(f"plots/signal_extraction_/{cent_bins[0]}_{cent_bins[1]}"):
        os.mkdir(f"plots/signal_extraction_/{cent_bins[0]}_{cent_bins[1]}")
    for ct_ in CT_BINS:
        if ct_[0] < 5 or ct_[1] > 40:
            continue
        for ct_bins in zip(CT_BINS_CENT[i_cent_bins][:-1], CT_BINS_CENT[i_cent_bins][1:]):
            if ct_bins[0] < ct_[0] or ct_bins[1] > ct_[1]:
                continue
        
            h_raw_yields = [[],[]] # 0 -> antim; 1 -> m
            h_raw_yields[0] = [ROOT.TH1D("fRawYields_pol","fRawYields_pol",100,0,1),ROOT.TH1D("fRawYields_expo","fRawYields_expo",100,0,1)]

            h_raw_yields[1] = [ROOT.TH1D("fRawYields_pol","fRawYields_pol",100,0,1),ROOT.TH1D("fRawYields_expo","fRawYields_expo",100,0,1)]

            bin = f'all_{cent_bins[0]}_{cent_bins[1]}_{ct_[0]}_{ct_[1]}_data_apply'
            df_data = pd.read_parquet(f'df_test/{bin}.parquet.gzip')
            #df_data = uproot.open(os.path.expandvars(f"/data/mciacco/LambdaPrompt_PbPb/df_data_{cent_bins[0]}_{cent_bins[1]}/AnalysisResults_lambda_{cent_bins[0]}_{cent_bins[1]}_ct_{ct_bins[0]}_{ct_bins[1]}.root"))['LambdaTreeBDTOut'].arrays(library="pd")
            df_data_ = df_data.query(f"mass > 1.09 and mass < 1.15 and ct >= {ct_bins[0]} and ct < {ct_bins[1]} and centrality > {cent_bins[0]} and centrality < {cent_bins[0]}")
            df_data_sidebands = df_data.query(f"(mass < 1.102 or mass > 1.13) and ct >= {ct_bins[0]} and ct < {ct_bins[1]}").sample(frac=0.5)

            for bdt_score in np.arange(0,1,0.05):
                print(f'processing {bin}: bkg cut = {bdt_score:.4f}')

                for i_split, split in enumerate(SPLIT_LIST):
                    split_ineq_sign = '> -0.1'
                    if SPLIT:
                        split_ineq_sign = '> 0.5'
                        if split == 'antimatter':
                            split_ineq_sign = '< 0.5'

                    # get invariant mass
                    df_data_cut = df_data_.query(f"model_output_0 < .2 and model_output_1 > {bdt_score} and pt < 3.5")
                    data_inv_mass_all = df_data_cut["mass"]
                    df_data_cut = df_data_cut.query(f"matter {split_ineq_sign}")
                    data_inv_mass = df_data_cut["mass"]

                    for i_bkg_func, bkg_function in enumerate(['pol','expo']):
                        # fit to invariant mass
                        inv_mass = ROOT.RooRealVar("m","#it{M} (p + #pi^{-})",1.09,1.15,"GeV/#it{c}^{2}")
                        inv_mass.setBins(N_BINS_INV_MASS)
                        inv_mass_roo_data = helpers.ndarray2roo(data_inv_mass.to_numpy(),inv_mass)
                        inv_mass_roo_data_all = helpers.ndarray2roo(data_inv_mass_all.to_numpy(),inv_mass)
                        inv_mass_data = ROOT.RooDataHist("db_m","db_m",ROOT.RooArgList(inv_mass),inv_mass_roo_data)
                        inv_mass_data_all = ROOT.RooDataHist("db_m_all","db_m_all",ROOT.RooArgList(inv_mass),inv_mass_roo_data_all)
                        #inv_mass_roo_mc = helpers.ndarray2roo(mc_inv_mass_sig.to_numpy(),inv_mass)
                        #inv_mass_mc = ROOT.RooDataHist("db_m_mc","db_m_mc",ROOT.RooArgList(inv_mass),inv_mass_roo_mc)
                        mass = ROOT.RooRealVar('#it{m}_{#Lambda}','mass',1.102,1.13,"GeV/#it{c}^{2}")
                        sigma = ROOT.RooRealVar('#sigma','sigma',0.001,0.004,"GeV/#it{c}^{2}")
                        alpha_left = ROOT.RooRealVar('#alpha_{left}','alpha_left',0.8,5.)
                        alpha_right = ROOT.RooRealVar('#alpha_{right}','alpha_right',.8,5.)
                        n_left = ROOT.RooRealVar('n_{left}','n_left',4.,30.)
                        n_right = ROOT.RooRealVar('n_{right}','n_right',4.,30.)
                        par_a = ROOT.RooRealVar('a','a',-20.,20.,"#it{c}^{2}/GeV")
                        par_b = ROOT.RooRealVar('b','b',-20.,20.,"#it{c}^{4}/GeV^{2}")
                        slope = ROOT.RooRealVar('slope','slope',-50.,50.,"#it{c}^{2}/GeV")
                        signal_pdf = ROOT.RooDSCBShape('signal','signal',inv_mass,mass,sigma,alpha_left,n_left,alpha_right,n_right)
                        #signal_pdf = ROOT.RooHistPdf("signal", "signal", ROOT.RooArgList(inv_mass), inv_mass_mc)
                        bkg_pdf = ROOT.RooPolynomial('bkg','bkg',inv_mass,ROOT.RooArgList(par_a,par_b))
                        if bkg_function == 'expo':
                            bkg_pdf = ROOT.RooExponential('bkg','bkg',inv_mass,slope)
                        w_signal = ROOT.RooRealVar("#it{f}_{signal}","w_signal",0.,1.)
                        n_tot = ROOT.RooRealVar("#it{N}_{tot}","n_tot",0.,1.e8)
                        model_mass_ = ROOT.RooAddPdf("model_mass_","model_mass_",signal_pdf,bkg_pdf,w_signal)
                        model_mass = ROOT.RooAddPdf("model_mass","model_mass",ROOT.RooArgList(model_mass_),ROOT.RooArgList(n_tot))
                        model_mass.fitTo(inv_mass_data_all,ROOT.RooFit.Save())
                        r=model_mass.fitTo(inv_mass_data_all,ROOT.RooFit.Save())

                        if (r.covQual() < 2.5):
                            continue

                        frame_mass = inv_mass.frame(ROOT.RooFit.Name(f"fMass_{ct_bins[0]}_{ct_bins[1]}_{bdt_score:.2f}"))
                        frame_mass.SetTitle(f"BDT score = {bdt_score:.2f}")
                        inv_mass_data.plotOn(frame_mass,ROOT.RooFit.Name("inv_mass_data"))
                        model_mass.plotOn(frame_mass,ROOT.RooFit.Components("bkg"),ROOT.RooFit.Name('background'),ROOT.RooFit.LineStyle(ROOT.kDashed),ROOT.RooFit.LineColor(ROOT.kGreen))
                        model_mass.plotOn(frame_mass,ROOT.RooFit.Components("signal"),ROOT.RooFit.Name('signal'),ROOT.RooFit.LineStyle(ROOT.kDashed),ROOT.RooFit.LineColor(ROOT.kOrange+7))
                        model_mass.plotOn(frame_mass,ROOT.RooFit.Name('model_mass'))
                        model_mass.paramOn(frame_mass, ROOT.RooFit.Parameters((mass, sigma)), ROOT.RooFit.Label("#chi^{2}/NDF = " + "{:.2f}".format(frame_mass.chiSquare("model_mass", "inv_mass_data"))), ROOT.RooFit.Layout(0.509565,0.829542,0.872206))

                        f = ROOT.TFile(f"fit_yield_noRadius_nodcaV0tracksCut_lambda_splitCentInMC_saveHistos_{cent_bins[0]}_{cent_bins[1]}.root","update")
                        f.cd()
                        f.mkdir(f'{split}_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}_{bkg_function}')
                        f.cd(f'{split}_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}_{bkg_function}')
                        frame_mass.getAttText().SetTextFont(44)
                        frame_mass.getAttText().SetTextSize(15)
                        frame_mass.getAttLine().SetLineWidth(0)
                        frame_mass.Write()
                        c_sim_fit = ROOT.TCanvas(f"cFit_{ct_bins[0]}_{ct_bins[1]}_{bdt_score:.2f}",f"cFit_{ct_bins[0]}_{ct_bins[1]}_{bdt_score:.2f}",1200,550)
                        c_sim_fit.cd()
                        frame_mass.Draw()
                        ROOT.gPad.Modified()
                        ROOT.gPad.Update()

                        if not os.path.isdir(f"plots/signal_extraction_/{cent_bins[0]}_{cent_bins[1]}/{ct_bins[0]}_{ct_bins[1]}"):
                            os.mkdir(f"plots/signal_extraction_/{cent_bins[0]}_{cent_bins[1]}/{ct_bins[0]}_{ct_bins[1]}")
                        c_sim_fit.Print(f"plots/signal_extraction_/{cent_bins[0]}_{cent_bins[1]}/{ct_bins[0]}_{ct_bins[1]}/{split}_{bkg_function}_{bdt_score:.2f}.pdf")
                        c_sim_fit.Write()

                        f.Close()

                        h_raw_yields[i_split][i_bkg_func].SetBinContent(h_raw_yields[i_split][i_bkg_func].FindBin(bdt_score+0.0001),n_tot.getVal()*w_signal.getVal())
                        dN_dn_tot = w_signal.getVal()
                        dN_dw_sig = n_tot.getVal()
                        sig_n_tot = n_tot.getError()
                        sig_w_signal = w_signal.getError()
                        cov_n_w_signal = r.correlation(w_signal,n_tot)*sig_n_tot*sig_w_signal
                        var_N = dN_dn_tot*dN_dn_tot*sig_n_tot*sig_n_tot + dN_dw_sig*dN_dw_sig*sig_w_signal*sig_w_signal + 2*dN_dn_tot*dN_dw_sig*cov_n_w_signal
                        h_raw_yields[i_split][i_bkg_func].SetBinError(h_raw_yields[i_split][i_bkg_func].FindBin(bdt_score+0.0001),np.sqrt(var_N))
                        print(f'************* fraction result = {w_signal.getVal()} *************')

            f = ROOT.TFile(f"fit_yield_noRadius_nodcaV0tracksCut_lambda_splitCentInMC_saveHistos_{cent_bins[0]}_{cent_bins[1]}.root","update")
            for i_split, split in enumerate(SPLIT_LIST):
                split_ineq_sign = '> -0.1'
                if SPLIT:
                    split_ineq_sign = '> 0.5'
                    if split == 'antimatter':
                        split_ineq_sign = '< 0.5'
                for i_bkg_func, bkg_function in enumerate(['pol','expo']):
                    f.cd(f'{split}_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}_{bkg_function}')
                    h_raw_yields[i_split][i_bkg_func].Write()
            f.Close()
