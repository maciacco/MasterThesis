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
SHM_XI_MINUS = [0.000305003,0.000329881,0.000356453,0.000384808,0.000415043,0.000447256,0.000481548,0.000518026,0.000556799,0.000597981,0.000641689,0.000688045,0.000737174,0.000789207,0.000844277,0.000902522,0.000964086,0.00102912,0.00109776,0.00117018,0.00124654]
SHM_XI_ZERO = [0.000311499,0.000336799,0.000363813,0.000392631,0.00042335,0.000456068,0.000490888,0.000527915,0.00056726,0.000609038,0.000653365,0.000700363,0.000750159,0.000802882,0.000858667,0.000917653,0.000979981,0.0010458,0.00111526,0.00118852,0.00126573]
SHM_OM = [4.5159e-05,4.94286e-05,5.40398e-05,5.9015e-05,6.43773e-05,7.0151e-05,7.63618e-05,8.30361e-05,9.02019e-05,9.78884e-05,0.000106126,0.000114946,0.000124381,0.000134467,0.000145239,0.000156733,0.000168989,0.000182047,0.000195948,0.000210736,0.000226455]
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

for i_cent_bins in range(len(CENTRALITY_LIST)):

    cent_bins = CENTRALITY_LIST[i_cent_bins]
    df_mc_tot = uproot.open(os.path.expandvars(f"/data/mciacco/LambdaPrompt_PbPb/AnalysisResults_reweight_BW_{cent_bins[0]}_{cent_bins[1]}.root"))['LambdaTree'].arrays(library="pd")
    df_mc_tot_om = df_mc_tot.query("ptMC > 0.5 and ptMC < 3.5 and flag==4")
    n_gen_om = df_mc_tot_om.shape[0]
    df_mc_xi0_tot = uproot.open(os.path.expandvars(f"/data/mciacco/LambdaPrompt_PbPb/AnalysisResults_Xi0_reweight_BW_{cent_bins[0]}_{cent_bins[1]}.root"))['LambdaTree'].arrays(library="pd")
    df_mc_xi0_tot_filter = df_mc_xi0_tot.query("ptMC > 0.5 and ptMC < 3.5 and flag==2")
    n_gen_xi0 = df_mc_xi0_tot_filter["weightMC"].sum()
    n_gen_xi = df_mc_xi0_tot_filter.shape[0]

    if not os.path.isdir(f"plots/signal_extraction_/{cent_bins[0]}_{cent_bins[1]}"):
        os.mkdir(f"plots/signal_extraction_/{cent_bins[0]}_{cent_bins[1]}")
    for ct_ in CT_BINS:
        if ct_[0] < 10 or ct_[1] > 40:
            continue
        for ct_bins in zip(CT_BINS_CENT[i_cent_bins][:-1], CT_BINS_CENT[i_cent_bins][1:]):
            if ct_bins[0] < ct_[0] or ct_bins[1] > ct_[1]:
                continue
        
            df_data = uproot.open(os.path.expandvars(f"/data/mciacco/LambdaPrompt_PbPb/df_data_{cent_bins[0]}_{cent_bins[1]}/AnalysisResults_lambda_{cent_bins[0]}_{cent_bins[1]}_ct_{ct_bins[0]}_{ct_bins[1]}.root"))['LambdaTreeBDTOut'].arrays(library="pd")
            df_data__ = df_data.query(f"mass > 1.09 and mass < 1.15 and ct >= {ct_bins[0]} and ct < {ct_bins[1]}")
            h_raw_yields = [[],[]] # 0 -> antim; 1 -> m
            h_raw_yields[0] = [ROOT.TH1D("fRawYields_pol","fRawYields_pol",100,0,1),ROOT.TH1D("fRawYields_expo","fRawYields_expo",100,0,1)]

            h_raw_yields[1] = [ROOT.TH1D("fRawYields_pol","fRawYields_pol",100,0,1),ROOT.TH1D("fRawYields_expo","fRawYields_expo",100,0,1)]

            bin = f'all_{cent_bins[0]}_{cent_bins[1]}_{ct_[0]}_{ct_[1]}_data_apply'
            #df_data = pd.read_parquet(f'df_test/{bin}.parquet.gzip')
            #df_data = uproot.open(os.path.expandvars(f"/data/mciacco/LambdaPrompt_PbPb/df_data_{cent_bins[0]}_{cent_bins[1]}/AnalysisResults_lambda_{cent_bins[0]}_{cent_bins[1]}_ct_{ct_bins[0]}_{ct_bins[1]}.root"))['LambdaTreeBDTOut'].arrays(library="pd")
            df_data_ = df_data__.query(f"mass > 1.09 and mass < 1.15 and ct >= {ct_bins[0]} and ct < {ct_bins[1]}")
            df_data_sidebands = df_data__.query(f"(mass < 1.102 or mass > 1.13) and ct >= {ct_bins[0]} and ct < {ct_bins[1]}").sample(frac=0.5)

            n_row = df_data_sidebands.shape[0]
            zero_list = []
            for _ in np.arange(n_row):
                zero_list.append(0)
            df_data_sidebands.loc[:,'y_true_from_flags'] = zero_list

            df_mc = pd.read_parquet(f"df_test/test_data_{cent_bins[0]}_{cent_bins[1]}_{ct_[0]}_{ct_[1]}_predict.parquet.gzip")
            
            # compute weights of Xi- and Xi0 feeddown components
            df_mc_prompt = df_mc.query('flag==1')
            df_mc_np_xi = df_mc.query('flag==2')
            df_mc_np_om = df_mc.query('flag==4')
            df_mc_np_xi0 = df_mc_np_xi.sample(frac=0.5)
            df_mc_np_xi.drop(df_mc_np_xi0.index)

            # merge mc data frames
            df_mc_np_xi0['flagXi0'] = df_mc_np_xi0['flag']
            df_mc__ = [df_mc_np_xi0,df_mc_np_xi,df_mc_np_om,df_mc_prompt]
            df_mc_ = pd.concat(df_mc__)
            df_mc_['bdtOutputBackground'] = df_mc_['model_output_background']
            df_mc_ = df_mc_.query(f'ct > {ct_bins[0]} and ct < {ct_bins[1]} and mass > 1.09 and mass < 1.15')

            eff_selected = np.arange(0.1, MAX_EFF, 0.01)
            y_truth_tmp = df_mc_['flag'].apply(lambda x : x if x == 1 else 0)
            score = analysis_utils.score_from_efficiency_array(
                y_truth_tmp, df_mc_['model_output_background'], eff_selected, keep_lower=True)
            # print(y_truth_tmp)
            # print(score)

            # del df_data_sidebands, df_mc_

            for bdt_score, bdt_eff in zip(score, eff_selected):
                if bdt_score < 0.2 or bdt_score > 0.30:
                    continue
                print(f'processing {bin}: bkg cut = {bdt_score:.4f}')

                
                for i_split, split in enumerate(SPLIT_LIST):
                    split_ineq_sign = '> -0.1'
                    if SPLIT:
                        split_ineq_sign = '> 0.5'
                        if split == 'antimatter':
                            split_ineq_sign = '< 0.5'
                    # apply cut
                    df_mc_cut = df_mc_.query(f"model_output_background < {bdt_score} and pt < 3.5")
                    print(f'size of df = {df_mc_cut.shape[0]}')
                    df_mc_cut_bkg = df_data_sidebands.query(f"matter {split_ineq_sign} and bdtOutputBackground < {bdt_score} and pt < 3.5")
                    eff = (df_mc_cut.query(f"flag==1").shape[0])/(df_mc_.query(f"flag==1").shape[0])

                    # get invariant mass
                    df_mc_cut["bdtOutputPrompt"] = df_mc_cut["model_output_prompt"]
                    df_data_cut = df_data_.query(f"bdtOutputBackground < {bdt_score} and pt < 3.5")
                    data_inv_mass_all = df_data_cut["mass"]
                    df_data_cut = df_data_cut.query(f"matter {split_ineq_sign}")
                    data_inv_mass = df_data_cut["mass"]
                    mc_inv_mass_sig = df_mc_cut.query("flag==1")["mass"]

                    # get prompt bdt prompt
                    data_bdt_prompt_out = df_data_cut["bdtOutputPrompt"]
                    mc_prompt_bdt_prompt_out = df_mc_cut.query("flag == 1")["model_output_prompt"]
                    mc_non_prompt_bdt_prompt_out_xiMinus = df_mc_cut.query("flag==2")["model_output_prompt"]
                    mc_non_prompt_bdt_prompt_out_Omega = df_mc_cut.query("flag==4")["model_output_prompt"]
                    mc_non_prompt_bdt_prompt_out_xiZero = df_mc_cut.query("flag==2 and flagXi0==2")["model_output_prompt"]
                    mc_non_prompt_weights_out_xiZero = df_mc_cut.query("flag==2 and flagXi0==2")["weightMC"]
                    mc_background_bdt_prompt_out = df_mc_cut_bkg["bdtOutputPrompt"]

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
                        mass = ROOT.RooRealVar('#it{m}_{#Lambda}','mass',1.11,1.12,"GeV/#it{c}^{2}")
                        sigma = ROOT.RooRealVar('#sigma','sigma',0.001,0.003,"GeV/#it{c}^{2}")
                        alpha_left = ROOT.RooRealVar('#alpha_{left}','alpha_left',1.,5.)
                        alpha_right = ROOT.RooRealVar('#alpha_{right}','alpha_right',1.,5.)
                        n_left = ROOT.RooRealVar('n_{left}','n_left',5.,30.)
                        n_right = ROOT.RooRealVar('n_{right}','n_right',5.,30.)
                        par_a = ROOT.RooRealVar('a','a',-20.,20.,"#it{c}^{2}/GeV")
                        par_b = ROOT.RooRealVar('b','b',-20.,20.,"#it{c}^{4}/GeV^{2}")
                        slope = ROOT.RooRealVar('slope','slope',-50.,50.,"#it{c}^{2}/GeV")
                        signal_pdf = ROOT.RooDSCBShape('signal','signal',inv_mass,mass,sigma,alpha_left,n_left,alpha_right,n_right)
                        #signal_pdf = ROOT.RooHistPdf("signal", "signal", ROOT.RooArgList(inv_mass), inv_mass_mc)
                        bkg_pdf = ROOT.RooChebychev('bkg','bkg',inv_mass,ROOT.RooArgList(par_a,par_b))
                        if bkg_function == 'expo':
                            bkg_pdf = ROOT.RooExponential('bkg','bkg',inv_mass,slope)
                        w_signal = ROOT.RooRealVar("#it{f}_{signal}","w_signal",0.,1.)
                        n_tot = ROOT.RooRealVar("#it{N}_{tot}","n_tot",0.,1.e9)
                        model_mass_ = ROOT.RooAddPdf("model_mass_","model_mass_",signal_pdf,bkg_pdf,w_signal)
                        model_mass = ROOT.RooAddPdf("model_mass","model_mass",ROOT.RooArgList(model_mass_),ROOT.RooArgList(n_tot))
                        for _ in range(2):
                            model_mass.fitTo(inv_mass_data_all)
                        sigma.setConstant()
                        mass.setConstant()
                        alpha_left.setConstant()
                        alpha_right.setConstant()
                        n_left.setConstant()
                        n_right.setConstant()

                        # fit to bdt output (2)
                        bdt_out = ROOT.RooRealVar("BDT out","BDT output for prompt",0.,1.)
                        bdt_out.setBins(N_BINS_OUTPUT_SCORE)
                        bdt_roo_data = helpers.ndarray2roo(data_bdt_prompt_out.to_numpy(), bdt_out)
                        bdt_data = ROOT.RooDataHist("dh", "dh", ROOT.RooArgList(bdt_out), bdt_roo_data)
                        bdt_roo_mc_prompt = helpers.ndarray2roo(mc_prompt_bdt_prompt_out.to_numpy(), bdt_out)
                        bdt_mc_prompt = ROOT.RooDataHist("dp", "dp", ROOT.RooArgList(bdt_out), bdt_roo_mc_prompt)
                        bdt_mc_prompt_pdf = ROOT.RooHistPdf("dppdf", "dppdf", ROOT.RooArgList(bdt_out), bdt_mc_prompt)
                        h_np_xi = ROOT.TH1D("h_np_xi","h_np_xi",N_BINS_OUTPUT_SCORE,0,1)
                        h_np_om = ROOT.TH1D("h_np_om","h_np_om",N_BINS_OUTPUT_SCORE,0,1)
                        h_np_xi_ = ROOT.TH1D("h_np_xi_zero","h_np_xi_zero",N_BINS_OUTPUT_SCORE,0,1)
                        h_np_xi_w = np.full((1,mc_non_prompt_bdt_prompt_out_xiMinus.shape[0]),1,dtype='float')
                        h_np_om_w = np.full((1,mc_non_prompt_bdt_prompt_out_Omega.shape[0]),1,dtype='float')
                        h_np_xi.FillN(mc_non_prompt_bdt_prompt_out_xiMinus.shape[0],np.asarray(mc_non_prompt_bdt_prompt_out_xiMinus.to_numpy(),dtype='float'),h_np_xi_w)
                        eff_xi_minus = h_np_xi.Integral()/n_gen_xi
                        h_np_xi_.FillN(mc_non_prompt_bdt_prompt_out_xiZero.shape[0],np.asarray(mc_non_prompt_bdt_prompt_out_xiZero.to_numpy(),dtype='float'),np.asarray(mc_non_prompt_weights_out_xiZero.to_numpy(),dtype='float'))
                        eff_xi_zero = h_np_xi_.Integral()/n_gen_xi0
                        h_np_om.FillN(mc_non_prompt_bdt_prompt_out_Omega.shape[0],np.asarray(mc_non_prompt_bdt_prompt_out_Omega.to_numpy(),dtype='float'),h_np_om_w)
                        eff_om = h_np_om.Integral()/n_gen_om
                        print(f"eff_xi_minus = {eff_xi_minus}, eff_xi_zero = {eff_xi_zero}, eff_om = {eff_om}")
                        f_xi_minus = SHM_XI_MINUS[10]*eff_xi_minus/(SHM_XI_MINUS[10]*eff_xi_minus+SHM_XI_ZERO[10]*eff_xi_zero+SHM_OM[10]*0.678*eff_om)
                        f_xi_zero = SHM_XI_ZERO[10]*eff_xi_zero/(SHM_XI_MINUS[10]*eff_xi_minus+SHM_XI_ZERO[10]*eff_xi_zero+SHM_OM[10]*0.678*eff_om)
                        f_om = 0.678*SHM_OM[10]*eff_om/(SHM_XI_MINUS[10]*eff_xi_minus+SHM_XI_ZERO[10]*eff_xi_zero+SHM_OM[10]*0.678*eff_om)
                        h_np_xi.Scale(f_xi_minus/h_np_xi.Integral())
                        h_np_xi_.Scale(f_xi_zero/h_np_xi_.Integral())
                        h_np_om.Scale(f_om/h_np_om.Integral())
                        h_np_xi_minus = ROOT.TH1D(h_np_xi)
                        h_np_xi_minus.SetName("h_np_xi_minus")
                        h_np_xi.Add(h_np_xi_)
                        h_np_xi.Add(h_np_om)
                        
                        bdt_mc_non_prompt = ROOT.RooDataHist("dnp", "dnp", ROOT.RooArgList(bdt_out), h_np_xi)
                        bdt_mc_non_prompt_pdf = ROOT.RooHistPdf("dnppdf", "dnppdf", ROOT.RooArgList(bdt_out), bdt_mc_non_prompt)
                        bdt_roo_mc_bkg = helpers.ndarray2roo(mc_background_bdt_prompt_out.to_numpy(), bdt_out)
                        bdt_mc_bkg = ROOT.RooDataHist("db", "db", ROOT.RooArgList(bdt_out), bdt_roo_mc_bkg)
                        bdt_mc_bkg_pdf = ROOT.RooHistPdf("dbpdf", "dbpdf", ROOT.RooArgList(bdt_out), bdt_mc_bkg)
                        frame = bdt_out.frame(ROOT.RooFit.Name(f"fBDTOutPrompt_{ct_bins[0]}_{ct_bins[1]}_{bdt_eff:.2f}"))
                        frame.SetTitle(str(ct_bins[0]) + " #leq #it{c}t < " + str(ct_bins[1]) + " cm")
                        w_non_prompt = ROOT.RooRealVar("#it{f}_{non-prompt}","w_non_prompt",0.,.5)
                        model_signal = ROOT.RooAddPdf("model_signal","model_signal",bdt_mc_non_prompt_pdf,bdt_mc_prompt_pdf,w_non_prompt)
                        model__ = ROOT.RooAddPdf("model__","model__",model_signal,bdt_mc_bkg_pdf,w_signal)
                        model = ROOT.RooAddPdf("model","model",ROOT.RooArgList(model__),ROOT.RooArgList(n_tot))

                        # simultaneous fit
                        sample = ROOT.RooCategory("sample","sample")
                        sample.defineType("bdt")
                        sample.defineType("mass")
                        combData = ROOT.RooDataHist("combData","combData",ROOT.RooArgSet(bdt_out,inv_mass),ROOT.RooFit.Index(sample),ROOT.RooFit.Import("bdt",bdt_data),ROOT.RooFit.Import("mass",inv_mass_data)) # ,ROOT.RooFit.Import("bdt_1",bdt_data_1),ROOT.RooFit.Import("bdt_0",bdt_data_0)
                        simPdf = ROOT.RooSimultaneous("simPdf","simPdf",sample)
                        simPdf.addPdf(model,"bdt")
                        simPdf.addPdf(model_mass,"mass")
                        simPdf.fitTo(combData,ROOT.RooFit.Save())
                        r = simPdf.fitTo(combData,ROOT.RooFit.Save())
                        if (r.covQual() < 1.5):
                            continue

                        # # # # # # # # PLOT FEATURES (REWEIGHTING STILL TO BE INTRODUCED!!!) # # # # # # # #
                        # splot features
                        # t_col_range = [[0.992,1],[0,1.2550588235294118],[0,10.120156862745098],[0,10.120156862745098],[0,2.626633712308453],[-5,5],[0,101.6],[0,10],[0,40],[0,1]]
                        # t_col_n_bins = [290,127,128,128,524,256,256,101,401,201]
                        # for score_regions_cuts in [[0,0.2],[0.8,1],[0,1]]:
                        #     continue
                        #     for i_col, t_col in enumerate(TRAINING_COLUMNS_LIST):
                        #         bdt_out.setRange("cutRange",score_regions_cuts[0],score_regions_cuts[1])
                        #         bdt_set = ROOT.RooArgSet(bdt_out)
                        #         norm_set = ROOT.RooFit.NormSet(bdt_set)
                        #         int_background_cut = (model__.pdfList().at(1).createIntegral(bdt_set, norm_set, ROOT.RooFit.Range("cutRange"))).getVal()
                        #         int_non_prompt_cut = (model_signal.pdfList().at(0).createIntegral(bdt_set, norm_set, ROOT.RooFit.Range("cutRange"))).getVal()
                        #         int_prompt_cut = (model_signal.pdfList().at(1).createIntegral(bdt_set, norm_set, ROOT.RooFit.Range("cutRange"))).getVal()
                        #         int_cut = (model.createIntegral(bdt_set, norm_set, ROOT.RooFit.Range("cutRange"))).getVal()
                        #         w_non_prompt_ = ROOT.RooRealVar("#it{N}_{non-prompt}_","w_non_prompt_",n_tot.getVal()*w_signal.getVal()*w_non_prompt.getVal()*int_non_prompt_cut,0.,1.e9)
                        #         w_prompt_ = ROOT.RooRealVar("#it{N}_{prompt}_","w_prompt_",n_tot.getVal()*w_signal.getVal()*(1-w_non_prompt.getVal())*int_prompt_cut,0.,1.e9)
                        #         w_bkg_ = ROOT.RooRealVar("#it{N}_{bkg}_","w_bkg_",n_tot.getVal()*(1-w_signal.getVal())*int_background_cut,0.,1.e9)
                        #         # if t_col != 'pt':
                        #         #     continue
                        #         print(t_col)
                        #         # get non-prompt bdt score
                        #         data_t_col = df_data_cut.query(f"bdtOutputPrompt > {score_regions_cuts[0]} and bdtOutputPrompt < {score_regions_cuts[1]}")[t_col]
                        #         min_t_col = data_t_col.min()
                        #         max_t_col = data_t_col.max()
                        #         mc_prompt_t_col = df_mc_cut.query(f"flag == 1 and model_output_prompt > {score_regions_cuts[0]} and model_output_prompt < {score_regions_cuts[1]}")[t_col]
                        #         mc_non_prompt_t_col = df_mc_cut.query(f"(flag == 2 or flag == 4) and (flag != 1 or flag != 0) and model_output_prompt > {score_regions_cuts[0]} and model_output_prompt < {score_regions_cuts[1]}")[t_col]
                        #         mc_background_t_col = df_mc_cut_bkg.query(f"y_true_from_flags==0 and bdtOutputPrompt > {score_regions_cuts[0]} and bdtOutputPrompt < {score_regions_cuts[1]}")[t_col]

                        #         t_col_var = ROOT.RooRealVar(f"{t_col}",f"{t_col}",t_col_range[i_col][0],t_col_range[i_col][1])
                        #         t_col_var.setBins(t_col_n_bins[i_col]-1)
                        #         t_col_var.setRange("full",t_col_range[i_col][0],t_col_range[i_col][1])
                        #         t_col_roo_data = helpers.ndarray2roo(data_t_col.to_numpy(), t_col_var)
                        #         t_col_var.setRange("norm",min_t_col,max_t_col)
                        #         t_col_data = ROOT.RooDataHist(f"dh_{t_col}", f"dh_{t_col}", ROOT.RooArgList(t_col_var), t_col_roo_data)
                        #         t_col_roo_mc_prompt = helpers.ndarray2roo(mc_prompt_t_col.to_numpy(), t_col_var)
                        #         t_col_mc_prompt = ROOT.RooDataHist(f"dp_{t_col}", f"dp_{t_col}", ROOT.RooArgList(t_col_var), t_col_roo_mc_prompt)
                        #         t_col_mc_prompt_pdf = ROOT.RooHistPdf(f"dp_{t_col}_pdf", f"dp_{t_col}_pdf", ROOT.RooArgList(t_col_var), t_col_mc_prompt)
                        #         t_col_roo_mc_non_prompt = helpers.ndarray2roo(mc_non_prompt_t_col.to_numpy(), t_col_var)
                        #         t_col_mc_non_prompt = ROOT.RooDataHist(f"dnp_{t_col}", f"dnp_{t_col}", ROOT.RooArgList(t_col_var), t_col_roo_mc_non_prompt)
                        #         t_col_mc_non_prompt_pdf = ROOT.RooHistPdf(f"dnp_{t_col}_pdf", f"dnp_{t_col}_pdf", ROOT.RooArgList(t_col_var), t_col_mc_non_prompt)
                        #         t_col_roo_mc_bkg = helpers.ndarray2roo(mc_background_t_col.to_numpy(), t_col_var)
                        #         t_col_mc_bkg = ROOT.RooDataHist(f"db_{t_col}_", f"db_{t_col}_", ROOT.RooArgList(t_col_var), t_col_roo_mc_bkg)
                        #         t_col_mc_bkg_pdf = ROOT.RooHistPdf(f"db_{t_col}_pdf", f"db_{t_col}_pdf", ROOT.RooArgList(t_col_var), t_col_mc_bkg)
                        #         frame_t_col = t_col_var.frame(ROOT.RooFit.Name(f"{t_col}_{ct_bins[0]}_{ct_bins[1]}"))
                        #         frame_t_col.SetTitle(str(ct_bins[0]) + " #leq #it{c}t < " + str(ct_bins[1]) + " cm")
                        #         model_t_col = ROOT.RooAddPdf(f"model_{t_col}",f"model_{t_col}",ROOT.RooArgList(t_col_mc_prompt_pdf,t_col_mc_non_prompt_pdf,t_col_mc_bkg_pdf),ROOT.RooArgList(w_prompt_,w_non_prompt_,w_bkg_))

                        #         t_col_data_ = t_col_data.createHistogram(t_col)
                        #         t_col_mc_prompt_ = t_col_mc_prompt.createHistogram(t_col)
                        #         t_col_mc_prompt_.Scale(w_prompt_.getVal()/t_col_mc_prompt_.Integral())
                        #         t_col_mc_non_prompt_ = t_col_mc_non_prompt.createHistogram(t_col)
                        #         t_col_mc_non_prompt_.Scale(w_non_prompt_.getVal()/t_col_mc_non_prompt_.Integral())
                        #         t_col_mc_bkg_ = t_col_mc_bkg.createHistogram(t_col)
                        #         t_col_mc_bkg_.Scale(w_bkg_.getVal()/t_col_mc_bkg_.Integral())
                        #         t_col_tot_hist = ROOT.TH1D("h","h",t_col_n_bins[i_col]-1,t_col_range[i_col][0],t_col_range[i_col][1])
                        #         t_col_tot_hist.Add(t_col_mc_prompt_)
                        #         t_col_tot_hist.Add(t_col_mc_non_prompt_)
                        #         t_col_tot_hist.Add(t_col_mc_bkg_)
                        #         c = ROOT.TCanvas(f"c_{t_col}_{score_regions_cuts[0]}_{score_regions_cuts[1]}",f"c_{t_col}_{score_regions_cuts[0]}_{score_regions_cuts[1]}")
                        #         c.cd()
                        #         c.SetLogy()
                        #         t_col_tot_hist.SetLineColor(ROOT.kBlue)
                        #         t_col_mc_prompt_.SetLineColor(ROOT.kRed)
                        #         t_col_mc_non_prompt_.SetLineColor(ROOT.kOrange)
                        #         t_col_mc_bkg_.SetLineColor(ROOT.kGreen)
                        #         t_col_tot_hist.GetXaxis().SetTitle(f"{t_col}")
                        #         t_col_tot_hist.GetYaxis().SetTitle(f"Entries/{t_col_tot_hist.GetXaxis().GetBinWidth(1)}")
                        #         t_col_tot_hist.Draw("histo")
                        #         t_col_data_.Draw("pesame")
                        #         t_col_mc_prompt_.Draw("histosame")
                        #         t_col_mc_non_prompt_.Draw("histosame")
                        #         t_col_mc_bkg_.Draw("histosame")

                        #         # t_col plot
                        #         t_col_data.plotOn(frame_t_col,ROOT.RooFit.Name(f"{t_col}_data"))
                        #         model_t_col.plotOn(frame_t_col,ROOT.RooFit.Components(f'db_{t_col}_pdf'),ROOT.RooFit.Name('background'),ROOT.RooFit.LineStyle(ROOT.kDashed),ROOT.RooFit.LineColor(ROOT.kGreen))
                        #         model_t_col.plotOn(frame_t_col,ROOT.RooFit.Components(f'dnp_{t_col}_pdf'),ROOT.RooFit.Name('non-prompt'),ROOT.RooFit.LineStyle(ROOT.kDashed),ROOT.RooFit.LineColor(ROOT.kOrange))
                        #         model_t_col.plotOn(frame_t_col,ROOT.RooFit.Components(f'dp_{t_col}_pdf'),ROOT.RooFit.Name('prompt'),ROOT.RooFit.LineStyle(ROOT.kDashed),ROOT.RooFit.LineColor(ROOT.kRed))
                        #         model_t_col.plotOn(frame_t_col,ROOT.RooFit.Name(f'model_{t_col}'),ROOT.RooFit.LineColor(ROOT.kBlue))
                        #         model_t_col.paramOn(frame_t_col, ROOT.RooFit.Label("#chi^{2}/NDF = " + "{:.2f}".format(frame_t_col.chiSquare(f"model_{t_col}", f"{t_col}_data"))), ROOT.RooFit.Layout(0.14959,0.47130,0.860544))

                        #         f = ROOT.TFile(f"fit_{cent_bins[0]}_{cent_bins[1]}.root","update")
                        #         f.cd()
                        #         f.mkdir(f'{split}_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}_{bkg_function}')
                        #         f.cd(f'{split}_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}_{bkg_function}')
                        #         #frame_t_col.Write()
                        #         c.Write()
                        #         f.Close()

                        # mass plot
                        mass.setConstant(False)
                        sigma.setConstant(False)
                        alpha_left.setConstant(False)
                        alpha_right.setConstant(False)
                        n_left.setConstant(False)
                        n_right.setConstant(False)
                        frame_mass = inv_mass.frame(ROOT.RooFit.Name(f"fMass_{ct_bins[0]}_{ct_bins[1]}_{bdt_eff:.2f}"))
                        frame_mass.SetTitle(f"BDT eff = {eff:.2f}")
                        inv_mass_data.plotOn(frame_mass,ROOT.RooFit.Name("inv_mass_data"))
                        model_mass.plotOn(frame_mass,ROOT.RooFit.Components("bkg"),ROOT.RooFit.Name('background'),ROOT.RooFit.LineStyle(ROOT.kDashed),ROOT.RooFit.LineColor(ROOT.kGreen))
                        model_mass.plotOn(frame_mass,ROOT.RooFit.Components("signal"),ROOT.RooFit.Name('signal'),ROOT.RooFit.LineStyle(ROOT.kDashed),ROOT.RooFit.LineColor(ROOT.kOrange+7))
                        model_mass.plotOn(frame_mass,ROOT.RooFit.Name('model_mass'))
                        
                        # bdt plot
                        bdt_data.plotOn(frame,ROOT.RooFit.Name("bdt_out_data"))
                        model.plotOn(frame,ROOT.RooFit.Components('dbpdf'),ROOT.RooFit.Name('background'),ROOT.RooFit.LineStyle(ROOT.kDashed),ROOT.RooFit.LineColor(ROOT.kGreen))
                        model.plotOn(frame,ROOT.RooFit.Components('dnppdf'),ROOT.RooFit.Name('non-prompt'),ROOT.RooFit.LineStyle(ROOT.kDashed),ROOT.RooFit.LineColor(ROOT.kOrange))
                        model.plotOn(frame,ROOT.RooFit.Components('dppdf'),ROOT.RooFit.Name('prompt'),ROOT.RooFit.LineStyle(ROOT.kDashed),ROOT.RooFit.LineColor(ROOT.kRed))
                        model.plotOn(frame,ROOT.RooFit.Name('model'),ROOT.RooFit.LineColor(ROOT.kBlue))
                        fit_param_bkg = 2
                        if bkg_function=='expo':
                            fit_param_bkg = 1
                        model_mass.paramOn(frame_mass, ROOT.RooFit.Parameters((mass, sigma)), ROOT.RooFit.Label("#chi^{2}/NDF = " + "{:.2f}".format((frame.chiSquare("model", "bdt_out_data")*frame.GetNbinsX()+frame_mass.chiSquare("model_mass", "inv_mass_data")*frame_mass.GetNbinsX())/(frame.GetNbinsX()+frame_mass.GetNbinsX()-fit_param_bkg-4))), ROOT.RooFit.Layout(0.509565,0.829542,0.872206))
                        model.paramOn(frame, ROOT.RooFit.Parameters((w_signal,w_non_prompt,n_tot)), ROOT.RooFit.Layout(0.160024,0.478262,0.33965))

                        h_model = model.createHistogram("BDT out")
                        h_data = bdt_data.createHistogram("BDT out")
                        h_model.Scale(h_data.Integral()/h_model.Integral())
                        c_ratio_data_MC = ROOT.TCanvas("pull","pull")
                        h_ratio = ROOT.TH1D()
                        c_ratio_data_MC.cd()
                        h_data.Divide(h_model)
                        h_data.Draw()

                        h_p = bdt_mc_prompt.createHistogram("BDT out")
                        h_np = bdt_mc_non_prompt.createHistogram("BDT out")
                        h_b = bdt_mc_bkg.createHistogram("BDT out")

                        f = ROOT.TFile(f"fit_{cent_bins[0]}_{cent_bins[1]}.root","update")
                        f.cd()
                        f.mkdir(f'{split}_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}_{bkg_function}')
                        f.cd(f'{split}_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}_{bkg_function}')
                        frame.getAttText().SetTextFont(44)
                        frame.getAttText().SetTextSize(15)
                        frame.getAttLine().SetLineWidth(0)
                        frame.Write()
                        #frame_1.Write()
                        #frame_0.Write()
                        frame_mass.getAttText().SetTextFont(44)
                        frame_mass.getAttText().SetTextSize(15)
                        frame_mass.getAttLine().SetLineWidth(0)
                        frame_mass.Write()
                        h_p.Write()
                        h_np.Write()
                        h_b.Write()
                        h_np_xi_minus.Write()
                        h_np_om.Write()
                        h_np_xi_.Write()
                        c_ratio_data_MC.Write()
                        c_sim_fit = ROOT.TCanvas(f"cSimFit_{ct_bins[0]}_{ct_bins[1]}_{bdt_eff:.2f}",f"cSimFit_{ct_bins[0]}_{ct_bins[1]}_{bdt_eff:.2f}",1200,550)
                        c_sim_fit.Divide(2,1)
                        c_sim_fit.cd(1)
                        c_sim_fit.GetPad(1).SetLogy()
                        frame.Draw()
                        c_sim_fit.cd(2)
                        frame_mass.Draw()
                        ROOT.gPad.Modified()
                        ROOT.gPad.Update()

                        if not os.path.isdir(f"plots/signal_extraction_/{cent_bins[0]}_{cent_bins[1]}/{ct_bins[0]}_{ct_bins[1]}"):
                            os.mkdir(f"plots/signal_extraction_/{cent_bins[0]}_{cent_bins[1]}/{ct_bins[0]}_{ct_bins[1]}")
                        c_sim_fit.Print(f"plots/signal_extraction_/{cent_bins[0]}_{cent_bins[1]}/{ct_bins[0]}_{ct_bins[1]}/{split}_{bkg_function}_{bdt_eff:.2f}.png")
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
                        print(f'************* fraction result = {w_signal.getVal()*(1-w_non_prompt.getVal())} *************')

            f = ROOT.TFile(f"fit_{cent_bins[0]}_{cent_bins[1]}.root","update")
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
