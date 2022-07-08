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

if not os.path.isdir("plots/signal_extraction"):
    os.mkdir("plots/signal_extraction")

for i_cent_bins in range(len(CENTRALITY_LIST)):

    cent_bins = CENTRALITY_LIST[i_cent_bins]
    if not os.path.isdir(f"plots/signal_extraction/{cent_bins[0]}_{cent_bins[1]}"):
        os.mkdir(f"plots/signal_extraction/{cent_bins[0]}_{cent_bins[1]}")
    for ct_ in CT_BINS:
        if ct_[0] < 10 or ct_[1] > 40:
            continue
        for ct_bins in zip(CT_BINS_CENT[i_cent_bins][:-1], CT_BINS_CENT[i_cent_bins][1:]):
            if ct_bins[0] < ct_[0] or ct_bins[1] > ct_[1]:
                continue
            # if ct_bins[0] < 32 or ct_bins[1] > ct_[1]:
            #     continue
        
            h_raw_yields = [[],[]] # 0 -> antim; 1 -> m
            h_raw_yields[0] = [ROOT.TH1D("fRawYields_pol","fRawYields_pol",100,0,1),ROOT.TH1D("fRawYields_expo","fRawYields_expo",100,0,1)]

            h_raw_yields[1] = [ROOT.TH1D("fRawYields_pol","fRawYields_pol",100,0,1),ROOT.TH1D("fRawYields_expo","fRawYields_expo",100,0,1)]

            bin = f'all_0_90_{ct_[0]}_{ct_[1]}'
            # df_data = pd.read_parquet(f'df/{bin}.parquet.gzip')
            df_data = uproot.open(os.path.expandvars(f"/data/mciacco/LambdaPrompt_PbPb/AnalysisResults_lambda_{cent_bins[0]}_{cent_bins[1]}.root"))['LambdaTreeBDTOut'].arrays(library="pd")
            df_data_sidebands = df_data.query(f"ct > {ct_bins[0]} and ct < {ct_bins[1]} and (mass < 1.095 or mass > 1.14)")
            df_data = df_data.query(f"ct > {ct_bins[0]} and ct < {ct_bins[1]} and mass > 1.095 and mass < 1.14")
            bin_mc_sig = f'all_mc_sig_0_90_{ct_[0]}_{ct_[1]}_allCandidates' # _reweight_{cent_bins[0]}_{cent_bins[1]}'
            bin_mc_bkg = f'mc_bkg_all_0_90_{ct_[0]}_{ct_[1]}'
            df_mc_sig = pd.read_parquet(f'df/{bin_mc_sig}.parquet.gzip')
            df_mc_bkg = pd.read_parquet(f'df/{bin_mc_bkg}.parquet.gzip')
            low_cent = cent_bins[0]
            up_cent = cent_bins[1]
            df_mc_bkg = df_mc_bkg.query(f'(mass < 1.105 or mass > 1.13) and ct > {ct_bins[0]} and ct < {ct_bins[1]} and mass > 1.095 and mass < 1.14 and centrality > {low_cent} and centrality < {up_cent}')
            df_mc_bkg_ = pd.read_parquet(f'df/{bin}.parquet.gzip')
            df_mc_bkg_ = df_mc_bkg_.query(f'(mass < 1.105 or mass > 1.13) and ct > {ct_bins[0]} and ct < {ct_bins[1]} and mass > 1.095 and mass < 1.14 and centrality > {low_cent} and centrality < {up_cent}')
            df_mc_sig = df_mc_sig.query(f"ct > {ct_bins[0]} and ct < {ct_bins[1]} and mass > 1.095 and mass < 1.14")
            # df_mc_no_bdt = df_mc_sig.copy()
            # n_row = df_mc_bkg_.shape[0]
            # zero_list = []
            # for _ in np.arange(n_row):
            #     zero_list.append(0)
            n_row = df_data_sidebands.shape[0]
            zero_list = []
            for _ in np.arange(n_row):
                zero_list.append(0)
            df_data_sidebands.loc[:,'y_true_from_flags'] = zero_list
            # print(df_mc_sig)
            df_mc_ = [df_data_sidebands, df_mc_sig] #, df_mc_bkg, df_mc_bkg_]
            df_mc = pd.concat(df_mc_)

            eff_selected = np.arange(0.1, MAX_EFF, 0.01)
            y_truth_tmp = df_mc_sig['flag'].apply(lambda x : x if x == 1 else 0)
            score = analysis_utils.score_from_efficiency_array(
                y_truth_tmp, df_mc_sig['model_output_0'], eff_selected, keep_lower=True)
            # print(y_truth_tmp)
            # print(score)

            del df_mc_sig , df_data_sidebands, df_mc_bkg_, df_mc_bkg

            for bdt_score, bdt_eff in zip(score, eff_selected):
                if bdt_eff < 0.75 or bdt_eff > 0.85: #or bdt_score > 0.15:
                    continue
                print(f'processing {bin}: bkg cut = {bdt_score:.4f}')

                # apply cut
                df_mc_cut = df_mc.query(f"(model_output_0 < {bdt_score} or bdtOutputBackground < {bdt_score})")
                
                for i_split, split in enumerate(SPLIT_LIST):
                    split_ineq_sign = '> -0.1'
                    if SPLIT:
                        split_ineq_sign = '> 0.5'
                        if split == 'antimatter':
                            split_ineq_sign = '< 0.5'
                    df_data_cut = df_data.query(f"bdtOutputBackground < {bdt_score} and matter {split_ineq_sign}")

                    # df_mc_no_bdt = df_mc.query(f"ct > {ct_bins[0]} and ct < {ct_bins[1]} and mass > 1.095 and mass < 1.14")

                    # get invariant mass
                    data_inv_mass = df_data_cut["mass"]
                    # mc_inv_mass_bkg = df_mc_cut.query("y_true_from_flags==0")["mass"]
                    # mc_inv_mass_sig = df_mc_cut.query("y_true_from_flags==2 or y_true_from_flags==1")["mass"]

                    # plt.hist(data_inv_mass, bins=1000)
                    # plt.savefig("plt.pdf")

                    # get non-prompt bdt score
                    data_bdt_non_prompt_out = df_data_cut["bdtOutputPrompt"]
                    mc_prompt_bdt_non_prompt_out = df_mc_cut.query("y_true_from_flags==2")["model_output_2"]
                    # mc_prompt_bdt_non_prompt_out_uncut = df_mc_no_bdt.query("y_true_from_flags==2")["model_output_2"]
                    # bdt_eff_corrected = mc_prompt_bdt_non_prompt_out.shape[0]/mc_prompt_bdt_non_prompt_out_uncut.shape[0]
                    # print(f"bdt_eff = {bdt_eff}")
                    mc_non_prompt_bdt_non_prompt_out = df_mc_cut.query("y_true_from_flags==1")["model_output_2"]
                    mc_background_bdt_non_prompt_out = df_mc_cut.query("y_true_from_flags==0")["bdtOutputPrompt"]

                    # mc_count_p = mc_prompt_bdt_non_prompt_out.count()
                    # mc_count_np = mc_non_prompt_bdt_non_prompt_out.count()

                    for i_bkg_func, bkg_function in enumerate(['pol','expo']):
                        # fit to invariant mass
                        inv_mass = ROOT.RooRealVar("m","#it{M} (p + #pi^{-})",1.095,1.14,"GeV/#it{c}^{2}")
                        inv_mass.setBins(100)
                        inv_mass_roo_data = helpers.ndarray2roo(data_inv_mass.to_numpy(),inv_mass)
                        inv_mass_data = ROOT.RooDataHist("db_m","db_m",ROOT.RooArgList(inv_mass),inv_mass_roo_data)
                        # inv_mass_roo_mc = helpers.ndarray2roo(mc_inv_mass_sig.to_numpy(),inv_mass)
                        # inv_mass_mc = ROOT.RooDataHist("db_m_mc","db_m_mc",ROOT.RooArgList(inv_mass),inv_mass_roo_mc)
                        mass = ROOT.RooRealVar('#it{m}_{#Lambda}','mass',1.11,1.12,"GeV/#it{c}^{2}")
                        sigma = ROOT.RooRealVar('#sigma','sigma',0.001,0.005,"GeV/#it{c}^{2}")
                        alpha_left = ROOT.RooRealVar('#alpha_{left}','alpha_left',1.,5.)
                        alpha_right = ROOT.RooRealVar('#alpha_{right}','alpha_right',1.,5.)
                        n_left = ROOT.RooRealVar('n_{left}','n_left',5.,20.)
                        n_right = ROOT.RooRealVar('n_{right}','n_right',5.,20.)
                        par_a = ROOT.RooRealVar('a','a',-20.,20.,"#it{c}^{2}/GeV")
                        par_b = ROOT.RooRealVar('b','b',-20.,20.,"#it{c}^{4}/GeV^{2}")
                        slope = ROOT.RooRealVar('slope','slope',-20.,20.,"#it{c}^{2}/GeV")
                        #signal_pdf = ROOT.RooGaussian('signal','signal',inv_mass,mass,sigma)
                        signal_pdf = ROOT.RooDSCBShape('signal','signal',inv_mass,mass,sigma,alpha_left,n_left,alpha_right,n_right)
                        bkg_pdf = ROOT.RooPolynomial('bkg','bkg',inv_mass,ROOT.RooArgList(par_a,par_b))
                        if bkg_function == 'expo':
                            bkg_pdf = ROOT.RooExponential('bkg','bkg',inv_mass,slope)
                        w_signal = ROOT.RooRealVar("#it{f}_{signal}","w_signal",0.,1.)
                        n_tot = ROOT.RooRealVar("#it{N}_{tot}","n_tot",0.,1.e8)
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
                        bdt_out = ROOT.RooRealVar("BDT out","BDT output for prompt",0.,1.)
                        bdt_out.setBins(10)
                        bdt_roo_data = helpers.ndarray2roo(data_bdt_non_prompt_out.to_numpy(), bdt_out)
                        bdt_data = ROOT.RooDataHist("dh", "dh", ROOT.RooArgList(bdt_out), bdt_roo_data)
                        bdt_roo_mc_prompt = helpers.ndarray2roo(mc_prompt_bdt_non_prompt_out.to_numpy(), bdt_out)
                        bdt_mc_prompt = ROOT.RooDataHist("dp", "dp", ROOT.RooArgList(bdt_out), bdt_roo_mc_prompt)
                        bdt_mc_prompt_pdf = ROOT.RooHistPdf("dppdf", "dppdf", ROOT.RooArgList(bdt_out), bdt_mc_prompt)
                        #bdt_mc_prompt_pdf = ROOT.RooParamHistFunc("dppdf", "dppdf", bdt_mc_prompt)
                        bdt_roo_mc_non_prompt = helpers.ndarray2roo(mc_non_prompt_bdt_non_prompt_out.to_numpy(), bdt_out)
                        bdt_mc_non_prompt = ROOT.RooDataHist("dnp", "dnp", ROOT.RooArgList(bdt_out), bdt_roo_mc_non_prompt)
                        bdt_mc_non_prompt_pdf = ROOT.RooHistPdf("dnppdf", "dnppdf", ROOT.RooArgList(bdt_out), bdt_mc_non_prompt)
                        #bdt_mc_non_prompt_pdf = ROOT.RooParamHistFunc("dnppdf", "dnppdf", bdt_mc_non_prompt)
                        bdt_roo_mc_bkg = helpers.ndarray2roo(mc_background_bdt_non_prompt_out.to_numpy(), bdt_out)
                        bdt_mc_bkg = ROOT.RooDataHist("db", "db", ROOT.RooArgList(bdt_out), bdt_roo_mc_bkg)
                        bdt_mc_bkg_pdf = ROOT.RooHistPdf("dbpdf", "dbpdf", ROOT.RooArgList(bdt_out), bdt_mc_bkg)
                        #bdt_mc_bkg_pdf = ROOT.RooParamHistFunc("dbpdf", "dbpdf", bdt_mc_bkg)
                        frame = bdt_out.frame(ROOT.RooFit.Name(f"fBDTOutNonPrompt_{ct_bins[0]}_{ct_bins[1]}_{bdt_eff:.2f}"))
                        frame.SetTitle(str(ct_bins[0]) + " #leq #it{c}t < " + str(ct_bins[1]) + " cm")
                        w_non_prompt = ROOT.RooRealVar("#it{f}_{non-prompt}","w_non_prompt",0.,1.)
                        model_signal = ROOT.RooAddPdf("model_signal","model_signal",bdt_mc_non_prompt_pdf,bdt_mc_prompt_pdf,w_non_prompt)
                        model__ = ROOT.RooAddPdf("model__","model__",model_signal,bdt_mc_bkg_pdf,w_signal)
                        model = ROOT.RooAddPdf("model","model",ROOT.RooArgList(model__),ROOT.RooArgList(n_tot))
                        # model_signal = ROOT.RooRealSumPdf("model_signal","model_signal",bdt_mc_non_prompt_pdf,bdt_mc_prompt_pdf,w_non_prompt)
                        # model_tmp = ROOT.RooRealSumPdf("model_tmp","model_tmp",model_signal,bdt_mc_bkg_pdf,w_signal)
                        # hc_prompt = ROOT.RooHistConstraint("hc_prompt","hc_prompt",bdt_mc_prompt_pdf)
                        # hc_non_prompt = ROOT.RooHistConstraint("hc_non_prompt","hc_non_prompt",bdt_mc_non_prompt_pdf)
                        # hc_background = ROOT.RooHistConstraint("hc_background","hc_background",bdt_mc_bkg_pdf)
                        # model = ROOT.RooProdPdf("model","model",ROOT.RooArgSet(hc_prompt,hc_non_prompt,hc_background),ROOT.RooFit.Conditional(model_tmp,bdt_out))

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
                        if (r.status() != 0) or (r.covQual() < 3):
                            continue

                        # mass plot
                        frame_mass = inv_mass.frame(ROOT.RooFit.Name(f"fMass_{ct_bins[0]}_{ct_bins[1]}_{bdt_eff:.2f}"))
                        frame_mass.SetTitle(f"BDT eff = {bdt_eff:.2f}")
                        inv_mass_data.plotOn(frame_mass,ROOT.RooFit.Name("inv_mass_data"))
                        model_mass.plotOn(frame_mass,ROOT.RooFit.Components("bkg"),ROOT.RooFit.Name('background'),ROOT.RooFit.LineStyle(ROOT.kDashed),ROOT.RooFit.LineColor(ROOT.kGreen))
                        model_mass.plotOn(frame_mass,ROOT.RooFit.Components("signal"),ROOT.RooFit.Name('signal'),ROOT.RooFit.LineStyle(ROOT.kDashed),ROOT.RooFit.LineColor(ROOT.kRed))
                        model_mass.plotOn(frame_mass,ROOT.RooFit.Name('model_mass'))
                        model_mass.paramOn(frame_mass, ROOT.RooFit.Label("#chi^{2}/NDF = " + "{:.2f}".format(frame_mass.chiSquare("model_mass", "inv_mass_data"))), ROOT.RooFit.Layout(0.546084,0.867801,0.883868))
                
                        # bdt plot
                        bdt_data.plotOn(frame,ROOT.RooFit.Name("bdt_out_data"))
                        model.plotOn(frame,ROOT.RooFit.Components('dbpdf'),ROOT.RooFit.Name('background'),ROOT.RooFit.LineStyle(ROOT.kDashed),ROOT.RooFit.LineColor(ROOT.kGreen))
                        model.plotOn(frame,ROOT.RooFit.Components('dnppdf'),ROOT.RooFit.Name('non-prompt'),ROOT.RooFit.LineStyle(ROOT.kDashed),ROOT.RooFit.LineColor(ROOT.kOrange))
                        model.plotOn(frame,ROOT.RooFit.Components('dppdf'),ROOT.RooFit.Name('prompt'),ROOT.RooFit.LineStyle(ROOT.kDashed),ROOT.RooFit.LineColor(ROOT.kRed))
                        model.plotOn(frame,ROOT.RooFit.Name('model'),ROOT.RooFit.LineColor(ROOT.kBlue))
                        model.paramOn(frame, ROOT.RooFit.Label("#chi^{2}/NDF = " + "{:.2f}".format(frame.chiSquare("model", "bdt_out_data"))), ROOT.RooFit.Layout(0.14959,0.47130,0.860544))

                        h_p = bdt_mc_prompt.createHistogram("BDT out")
                        h_np = bdt_mc_non_prompt.createHistogram("BDT out")
                        h_b = bdt_mc_bkg.createHistogram("BDT out")

                        f = ROOT.TFile(f"fit_lambda_splitCentInMC_saveHistos_{cent_bins[0]}_{cent_bins[1]}.root","update")
                        f.cd()
                        f.mkdir(f'{split}_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}_{bkg_function}')
                        f.cd(f'{split}_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}_{bkg_function}')
                        frame.getAttText().SetTextFont(44)
                        frame.getAttText().SetTextSize(15)
                        frame.getAttLine().SetLineWidth(0)
                        frame.Write()
                        frame_mass.getAttText().SetTextFont(44)
                        frame_mass.getAttText().SetTextSize(15)
                        frame_mass.getAttLine().SetLineWidth(0)
                        frame_mass.Write()
                        h_p.Write()
                        h_np.Write()
                        h_b.Write()
                        c_sim_fit = ROOT.TCanvas(f"cSimFit_{ct_bins[0]}_{ct_bins[1]}_{bdt_eff:.2f}",f"cSimFit_{ct_bins[0]}_{ct_bins[1]}_{bdt_eff:.2f}",1200,550)
                        c_sim_fit.Divide(2,1)
                        c_sim_fit.cd(1)
                        frame.Draw()
                        c_sim_fit.cd(2)
                        frame_mass.Draw()

                        if not os.path.isdir(f"plots/signal_extraction/{cent_bins[0]}_{cent_bins[1]}/{ct_bins[0]}_{ct_bins[1]}"):
                            os.mkdir(f"plots/signal_extraction/{cent_bins[0]}_{cent_bins[1]}/{ct_bins[0]}_{ct_bins[1]}")
                        c_sim_fit.Print(f"plots/signal_extraction/{cent_bins[0]}_{cent_bins[1]}/{ct_bins[0]}_{ct_bins[1]}/{split}_{bkg_function}_{bdt_eff:.2f}.pdf")
                        c_sim_fit.Write()

                        f.Close()

                        h_raw_yields[i_split][i_bkg_func].SetBinContent(h_raw_yields[i_split][i_bkg_func].FindBin(bdt_eff+0.0001),n_tot.getVal()*w_signal.getVal()*(1-w_non_prompt.getVal())/bdt_eff)
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
                        h_raw_yields[i_split][i_bkg_func].SetBinError(h_raw_yields[i_split][i_bkg_func].FindBin(bdt_eff+0.0001),np.sqrt(var_N)/bdt_eff)
                        #print(f'Non-prompt / prompt = {mc_count_np/(mc_count_p+mc_count_np)}')

            f = ROOT.TFile(f"fit_lambda_splitCentInMC_saveHistos_{cent_bins[0]}_{cent_bins[1]}.root","update")
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