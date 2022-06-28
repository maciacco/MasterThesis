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

SPLIT = False
MAX_EFF = 1

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

for split in SPLIT_LIST:
    split_ineq_sign = '> -0.1'
    if SPLIT:
        split_ineq_sign = '> 0.5'
        if split == 'antimatter':
            split_ineq_sign = '< 0.5'

    for i_cent_bins in range(len(CENTRALITY_LIST)):
        # df_data = uproot.open(os.path.expandvars("/data/mciacco/LambdaPrompt_PbPb/AnalysisResults.root"))['LambdaTreeBDTOut'].arrays(library="pd")

        cent_bins = CENTRALITY_LIST[i_cent_bins]
        for ct_ in CT_BINS:
            if ct_[0] < 10 or ct_[1] > 40:
                continue
            for ct_bins in zip(CT_BINS_CENT[i_cent_bins][:-1], CT_BINS_CENT[i_cent_bins][1:]):
                if ct_bins[0] < ct_[0] or ct_bins[1] > ct_[1]:
                    continue
                # if ct_bins[0] < 32 or ct_bins[1] > ct_[1]:
                #     continue
            
                for bkg_function in ['pol','expo']:
                    h_raw_yields = ROOT.TH1D("fRawYields","fRawYields",100,0,1)

                    bin = f'{split}_{cent_bins[0]}_{cent_bins[1]}_{ct_[0]}_{ct_[1]}'
                    df_data = pd.read_parquet(f'df/{bin}.parquet.gzip')

                    bin_mc_sig = f'mc_sig_all_{cent_bins[0]}_{cent_bins[1]}_{ct_[0]}_{ct_[1]}'
                    bin_mc_bkg = f'mc_bkg_all_{cent_bins[0]}_{cent_bins[1]}_{ct_[0]}_{ct_[1]}'
                    df_mc_sig = pd.read_parquet(f'df/{bin_mc_sig}.parquet.gzip')
                    df_mc_sig = df_mc_sig.query(f"ct > {ct_bins[0]} and ct < {ct_bins[1]} and mass > 1.09 and mass < 1.15")
                    df_mc_bkg = pd.read_parquet(f'df/{bin_mc_bkg}.parquet.gzip')
                    df_mc_ = [df_mc_sig, df_mc_bkg]
                    df_mc = pd.concat(df_mc_)
                    df_mc = df_mc.query(f"ct > {ct_bins[0]} and ct < {ct_bins[1]} and mass > 1.09 and mass < 1.15")

                    eff_selected = np.arange(0.1, MAX_EFF, 0.04)
                    y_truth_tmp = df_mc_sig['flag'].apply(lambda x : x if x == 1 else 0)
                    score = analysis_utils.score_from_efficiency_array(
                        y_truth_tmp, df_mc_sig['model_output_0'], eff_selected, keep_lower=True)
                    #print(y_truth_tmp)
                    #print(score)

                    del df_mc_sig, df_mc_bkg

                    for bdt_score, bdt_eff in zip(score, eff_selected):
                        if bdt_eff < 0.85 or bdt_eff > 0.87:
                            continue
                        print(f'processing {bin}: bkg cut = {bdt_score:.4f}')

                        # apply cut
                        df_data_cut = df_data.query(f"model_output_0 < {bdt_score} and ct > {ct_bins[0]} and ct < {ct_bins[1]} and mass > 1.09 and mass < 1.15")
                        df_mc_cut = df_mc.query(f"model_output_0 < {bdt_score} and ct > {ct_bins[0]} and ct < {ct_bins[1]} and mass > 1.09 and mass < 1.15")
                        # df_mc_no_bdt = df_mc.query(f"ct > {ct_bins[0]} and ct < {ct_bins[1]} and mass > 1.09 and mass < 1.15")

                        # get invariant mass
                        data_inv_mass = df_data_cut["mass"]
                        mc_inv_mass_bkg = df_mc_cut.query("y_true_from_flags==0")["mass"]
                        mc_inv_mass_sig = df_mc_cut.query("y_true_from_flags==2 or y_true_from_flags==1")["mass"]

                        plt.hist(data_inv_mass, bins=1000)
                        plt.savefig("plt.pdf")

                        # get non-prompt bdt score
                        data_bdt_non_prompt_out = df_data_cut["model_output_2"]
                        mc_prompt_bdt_non_prompt_out = df_mc_cut.query("y_true_from_flags==2")["model_output_2"]
                        # mc_prompt_bdt_non_prompt_out_uncut = df_mc_no_bdt.query("y_true_from_flags==2")["model_output_2"]
                        # bdt_eff = mc_prompt_bdt_non_prompt_out.count()/mc_prompt_bdt_non_prompt_out_uncut.count()
                        # print(f"bdt_eff = {bdt_eff}")
                        mc_non_prompt_bdt_non_prompt_out = df_mc_cut.query("y_true_from_flags==1")["model_output_2"]
                        mc_background_bdt_non_prompt_out = df_mc_cut.query("y_true_from_flags==0")["model_output_2"]

                        # mc_count_p = mc_prompt_bdt_non_prompt_out.count()
                        # mc_count_np = mc_non_prompt_bdt_non_prompt_out.count()

                        # fit to invariant mass
                        inv_mass = ROOT.RooRealVar("m","m",1.09,1.15)
                        inv_mass.setBins(50)
                        inv_mass_roo_data = helpers.ndarray2roo(data_inv_mass.to_numpy(),inv_mass)
                        inv_mass_data = ROOT.RooDataHist("db_m","db_m",ROOT.RooArgList(inv_mass),inv_mass_roo_data)
                        inv_mass_roo_mc = helpers.ndarray2roo(mc_inv_mass_sig.to_numpy(),inv_mass)
                        inv_mass_mc = ROOT.RooDataHist("db_m_mc","db_m_mc",ROOT.RooArgList(inv_mass),inv_mass_roo_mc)
                        mass = ROOT.RooRealVar('mass','mass',1.11,1.12)
                        sigma = ROOT.RooRealVar('sigma','sigma',0.001,0.005)
                        alpha_left = ROOT.RooRealVar('alpha_left','alpha_left',0.,2.)
                        alpha_right = ROOT.RooRealVar('alpha_right','alpha_right',0.,2.)
                        n_left = ROOT.RooRealVar('n_left','n_left',0.,15.)
                        n_right = ROOT.RooRealVar('n_right','n_right',0.,15.)
                        par_a = ROOT.RooRealVar('a','a',-20.,20.)
                        par_b = ROOT.RooRealVar('b','b',-20.,20.)
                        slope = ROOT.RooRealVar('slope','slope',-20.,20.)
                        #signal_pdf = ROOT.RooGaussian('signal','signal',inv_mass,mass,sigma)
                        signal_pdf = ROOT.RooDSCBShape('signal','signal',inv_mass,mass,sigma,alpha_left,n_left,alpha_right,n_right)
                        bkg_pdf = ROOT.RooPolynomial('bkg','bkg',inv_mass,ROOT.RooArgList(par_a,par_b))
                        if bkg_function == 'expo':
                            bkg_pdf = ROOT.RooExponential('bkg','bkg',inv_mass,slope)
                        w_signal = ROOT.RooRealVar("w_signal","w_signal",0.,1.)
                        n_tot = ROOT.RooRealVar("n_tot","n_tot",0.,1.e8)
                        model_mass_ = ROOT.RooAddPdf("model_mass_","model_mass_",signal_pdf,bkg_pdf,w_signal)
                        model_mass = ROOT.RooAddPdf("model_mass","model_mass",ROOT.RooArgList(model_mass_),ROOT.RooArgList(n_tot))
                        model_mass.fitTo(inv_mass_data)
                        # sigma.setConstant()
                        # mass.setConstant()
                        # alpha_left.setConstant()
                        # alpha_right.setConstant()
                        # n_left.setConstant()
                        # n_right.setConstant()
                        # par_a.setConstant()
                        # par_b.setConstant()
                        # n_tot.setConstant()
                        # w_signal.setConstant()

                        # fit to bdt output
                        bdt_out = ROOT.RooRealVar("BDT out","BDT out",0.,1.)
                        bdt_out.setBins(10)
                        bdt_roo_data = helpers.ndarray2roo(data_bdt_non_prompt_out.to_numpy(), bdt_out)
                        bdt_data = ROOT.RooDataHist("dh", "dh", ROOT.RooArgList(bdt_out), bdt_roo_data)
                        bdt_roo_mc_prompt = helpers.ndarray2roo(mc_prompt_bdt_non_prompt_out.to_numpy(), bdt_out)
                        bdt_mc_prompt = ROOT.RooDataHist("dp", "dp", ROOT.RooArgList(bdt_out), bdt_roo_mc_prompt)
                        # bdt_mc_prompt_pdf = ROOT.RooHistPdf("dppdf", "dppdf", ROOT.RooArgList(bdt_out), bdt_mc_prompt)
                        bdt_mc_prompt_pdf = ROOT.RooParamHistFunc("dppdf", "dppdf", bdt_mc_prompt)
                        bdt_roo_mc_non_prompt = helpers.ndarray2roo(mc_non_prompt_bdt_non_prompt_out.to_numpy(), bdt_out)
                        bdt_mc_non_prompt = ROOT.RooDataHist("dnp", "dnp", ROOT.RooArgList(bdt_out), bdt_roo_mc_non_prompt)
                        # bdt_mc_non_prompt_pdf = ROOT.RooHistPdf("dnppdf", "dnppdf", ROOT.RooArgList(bdt_out), bdt_mc_non_prompt)
                        bdt_mc_non_prompt_pdf = ROOT.RooParamHistFunc("dnppdf", "dnppdf", bdt_mc_non_prompt)
                        bdt_roo_mc_bkg = helpers.ndarray2roo(mc_background_bdt_non_prompt_out.to_numpy(), bdt_out)
                        bdt_mc_bkg = ROOT.RooDataHist("db", "db", ROOT.RooArgList(bdt_out), bdt_roo_mc_bkg)
                        # bdt_mc_bkg_pdf = ROOT.RooHistPdf("dbpdf", "dbpdf", ROOT.RooArgList(bdt_out), bdt_mc_bkg)
                        bdt_mc_bkg_pdf = ROOT.RooParamHistFunc("dbpdf", "dbpdf", bdt_mc_bkg)
                        frame = bdt_out.frame(ROOT.RooFit.Name(f"fBDTOutNonPrompt_{ct_bins[0]}_{ct_bins[1]}_{bdt_eff:.2f}"))
                        w_non_prompt = ROOT.RooRealVar("w_non_prompt","w_non_prompt",0.,1.)
                        #model_signal = ROOT.RooAddPdf("model_signal","model_signal",bdt_mc_non_prompt_pdf,bdt_mc_prompt_pdf,w_non_prompt)
                        #model = ROOT.RooAddPdf("model","model",model_signal,bdt_mc_bkg_pdf,w_signal)
                        
                        hc_prompt = ROOT.RooHistConstraint("hc_prompt","hc_prompt",bdt_mc_prompt_pdf)
                        hc_non_prompt = ROOT.RooHistConstraint("hc_non_prompt","hc_non_prompt",bdt_mc_non_prompt_pdf)
                        hc_background = ROOT.RooHistConstraint("hc_background","hc_background",bdt_mc_bkg_pdf)

                        h_p = bdt_mc_prompt.createHistogram("BDT out")
                        h_np = bdt_mc_non_prompt.createHistogram("BDT out")
                        h_b = bdt_mc_bkg.createHistogram("BDT out")
                        h_p.Scale(1/h_p.GetSumOfWeights())
                        h_np.Scale(1/h_np.GetSumOfWeights())
                        h_b.Scale(1/h_b.GetSumOfWeights())

                        bdt_mc_prompt_ = ROOT.RooDataHist("dp_", "dp_", ROOT.RooArgList(bdt_out), h_p)
                        bdt_mc_non_prompt_ = ROOT.RooDataHist("dnp_", "dnp_", ROOT.RooArgList(bdt_out), h_np)
                        bdt_mc_bkg_ = ROOT.RooDataHist("db_", "db_", ROOT.RooArgList(bdt_out), h_b)
                        bdt_mc_prompt_pdf_ = ROOT.RooParamHistFunc("dppdf_", "dppdf_", bdt_mc_prompt_, bdt_mc_prompt_pdf)
                        bdt_mc_non_prompt_pdf_ = ROOT.RooParamHistFunc("dnppdf_", "dnppdf_", bdt_mc_non_prompt_, bdt_mc_non_prompt_pdf)
                        bdt_mc_bkg_pdf_ = ROOT.RooParamHistFunc("dbpdf_", "dbpdf_", bdt_mc_bkg_, bdt_mc_bkg_pdf)

                        model_signal = ROOT.RooRealSumPdf("model_signal","model_signal",bdt_mc_non_prompt_pdf_,bdt_mc_prompt_pdf_,w_non_prompt)
                        model_tmp = ROOT.RooRealSumPdf("model_tmp","model_tmp",model_signal,bdt_mc_bkg_pdf_,w_signal)
                        model = ROOT.RooProdPdf("model","model",ROOT.RooArgSet(hc_background,hc_non_prompt,hc_prompt),ROOT.RooFit.Conditional(ROOT.RooArgSet(model_tmp),ROOT.RooArgSet(bdt_out)))
                        list_pdf = model.getConnectedParameters(bdt_out)

                        # simultaneous fit
                        sample = ROOT.RooCategory("sample","sample")
                        sample.defineType("bdt")
                        sample.defineType("mass")
                        combData = ROOT.RooDataHist("combData","combData",ROOT.RooArgSet(bdt_out,inv_mass),ROOT.RooFit.Index(sample),ROOT.RooFit.Import("bdt",bdt_data),ROOT.RooFit.Import("mass",inv_mass_data))
                        simPdf = ROOT.RooSimultaneous("simPdf","simPdf",sample)
                        simPdf.addPdf(model,"bdt")
                        simPdf.addPdf(model_mass,"mass")
                        simPdf.fitTo(combData,ROOT.RooFit.Save())
                        r = simPdf.fitTo(combData,ROOT.RooFit.Save())
                        if r.status() != 0:
                            continue

                        # mass plot
                        frame_mass = inv_mass.frame(ROOT.RooFit.Name(f"fMass_{ct_bins[0]}_{ct_bins[1]}_{bdt_eff:.2f}"))
                        inv_mass_data.plotOn(frame_mass,ROOT.RooFit.Name("inv_mass_data"))
                        model_mass.plotOn(frame_mass,ROOT.RooFit.Components("bkg"),ROOT.RooFit.Name('background'),ROOT.RooFit.LineStyle(ROOT.kDashed),ROOT.RooFit.LineColor(ROOT.kGreen))
                        model_mass.plotOn(frame_mass,ROOT.RooFit.Components("signal"),ROOT.RooFit.Name('signal'),ROOT.RooFit.LineStyle(ROOT.kDashed),ROOT.RooFit.LineColor(ROOT.kRed))
                        model_mass.plotOn(frame_mass,ROOT.RooFit.Name('model_mass'))
                        model_mass.paramOn(frame_mass, ROOT.RooFit.Label("#chi^{2}/NDF = " + "{:.2f}".format(frame_mass.chiSquare("model_mass", "inv_mass_data"))), ROOT.RooFit.Layout(0.58812,0.911028,0.861955))
                
                        # bdt plot
                        bdt_data.plotOn(frame)
                        model.plotOn(frame,ROOT.RooFit.Components('dbpdf_'),ROOT.RooFit.Name('background'),ROOT.RooFit.LineStyle(ROOT.kDashed),ROOT.RooFit.LineColor(ROOT.kGreen))
                        model.plotOn(frame,ROOT.RooFit.Components('dnppdf_'),ROOT.RooFit.Name('non-prompt'),ROOT.RooFit.LineStyle(ROOT.kDashed),ROOT.RooFit.LineColor(ROOT.kOrange))
                        model.plotOn(frame,ROOT.RooFit.Components('dppdf_'),ROOT.RooFit.Name('prompt'),ROOT.RooFit.LineStyle(ROOT.kDashed),ROOT.RooFit.LineColor(ROOT.kRed))
                        model.plotOn(frame,ROOT.RooFit.Name('model'),ROOT.RooFit.LineColor(ROOT.kBlue))

                        f = ROOT.TFile("out_1_bb_.root","update")
                        f.cd()
                        f.mkdir(f'{split}_{cent_bins[0]}_{cent_bins[0]}_{ct_bins[0]}_{ct_bins[1]}_{bkg_function}')
                        f.cd(f'{split}_{cent_bins[0]}_{cent_bins[0]}_{ct_bins[0]}_{ct_bins[1]}_{bkg_function}')
                        frame.Write()
                        frame_mass.Write()
                        # h_p.Write()
                        # h_np.Write()
                        # h_b.Write()
                        f.Close()

                        h_raw_yields.SetBinContent(h_raw_yields.FindBin(bdt_eff+0.0001),n_tot.getVal()*w_signal.getVal()*(1-w_non_prompt.getVal())/bdt_eff)
                        dN_dn_tot = w_signal.getVal()*(1-w_non_prompt.getVal())
                        dN_dw_sig = n_tot.getVal()*(1-w_non_prompt.getVal())
                        dN_dw_non_prompt = -n_tot.getVal()*w_signal.getVal()
                        sig_n_tot = n_tot.getError()
                        sig_w_signal = w_signal.getError()
                        sig_w_non_prompt = w_non_prompt.getError()
                        cov_w_w = 0
                        cov_n_w_signal = 0
                        cov_n_w_np = 0
                        cov_w_w = r.correlation(w_signal,w_non_prompt)*sig_w_non_prompt*sig_w_signal
                        cov_n_w_signal = r.correlation(w_signal,n_tot)*sig_n_tot*sig_w_signal
                        cov_n_w_np = r.correlation(w_non_prompt,n_tot)*sig_n_tot*sig_w_non_prompt
                        var_N = dN_dn_tot*dN_dn_tot*sig_n_tot*sig_n_tot + dN_dw_sig*dN_dw_sig*sig_w_signal*sig_w_signal + dN_dw_non_prompt*dN_dw_non_prompt*sig_w_non_prompt*sig_w_non_prompt + 2*dN_dn_tot*dN_dw_sig*cov_n_w_signal + 2*dN_dn_tot*dN_dw_non_prompt*cov_n_w_np + 2*dN_dw_sig*dN_dw_non_prompt*cov_w_w
                        h_raw_yields.SetBinError(h_raw_yields.FindBin(bdt_eff+0.0001),np.sqrt(var_N)/bdt_eff)
                        #print(f'Non-prompt / prompt = {mc_count_np/(mc_count_p+mc_count_np)}')

                    f = ROOT.TFile("out_1_bb_.root","update")
                    f.cd(f'{split}_{cent_bins[0]}_{cent_bins[0]}_{ct_bins[0]}_{ct_bins[1]}_{bkg_function}')
                    h_raw_yields.Write()
                    f.Close()