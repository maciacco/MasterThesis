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
from ctypes import c_double

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

colors = [ROOT.kRed,ROOT.kOrange,ROOT.kGreen]

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
            # if ct_bins[0] < 12 or ct_bins[1] > ct_[1]:
            #     continue
        
            h_raw_yields = [[],[]] # 0 -> antim; 1 -> m
            h_raw_yields[0] = [ROOT.TH1D("fRawYields_pol","fRawYields_pol",100,0,1),ROOT.TH1D("fRawYields_expo","fRawYields_expo",100,0,1)]

            h_raw_yields[1] = [ROOT.TH1D("fRawYields_pol","fRawYields_pol",100,0,1),ROOT.TH1D("fRawYields_expo","fRawYields_expo",100,0,1)]

            bin = f'all_{cent_bins[0]}_{cent_bins[1]}_{ct_[0]}_{ct_[1]}_data_apply'
            df_data = pd.read_parquet(f'df_test/{bin}.parquet.gzip')
            #df_data = uproot.open(os.path.expandvars(f"/data/mciacco/LambdaPrompt_PbPb/df_data_{cent_bins[0]}_{cent_bins[1]}/AnalysisResults_lambda_{cent_bins[0]}_{cent_bins[1]}_ct_{ct_bins[0]}_{ct_bins[1]}.root"))['LambdaTreeBDTOut'].arrays(library="pd")
            df_data_ = df_data.query(f"mass > 1.09 and mass < 1.15 and ct >= {ct_bins[0]} and ct < {ct_bins[1]} and pt > 0.5 and pt < 3.5") #and centrality >= {cent_bins[0]} and centrality < {cent_bins[1]} and pt > 0.5 and pt < 3.5 and ct > {ct_bins[0]} and ct < {ct_bins[1]} and tpcClV0Pi > 69 and tpcClV0Pr > 69 and radius > 3 and radius < 100 and dcaPrPV < 20 and dcaPiPV < 20 and dcaV0PV < 10 and eta < 0.8 and eta > -0.8")
            df_data_sidebands = df_data.query(f"(mass < 1.102 or mass > 1.13) and ct >= {ct_bins[0]} and ct < {ct_bins[1]} and pt > 0.5 and pt < 3.5").sample(frac=0.5) # and ct >= {ct_bins[0]} and ct < {ct_bins[1]} and centrality >= {cent_bins[0]} and centrality < {cent_bins[1]} and pt > 0.5 and pt < 3.5 and ct > {ct_bins[0]} and ct < {ct_bins[1]} and tpcClV0Pi > 69 and tpcClV0Pr > 69 and radius > 3 and radius < 100 and dcaPrPV < 20 and dcaPiPV < 20 and dcaV0PV < 10 and eta < 0.8 and eta > -0.8").sample(frac=0.5)
            # bin_mc_sig = f'all_mc_sig_0_90_{ct_[0]}_{ct_[1]}_allCandidates' # _reweight_{cent_bins[0]}_{cent_bins[1]}'
            # bin_mc_bkg = f'mc_bkg_all_0_90_{ct_[0]}_{ct_[1]}'
            # df_mc_sig = pd.read_parquet(f'df_test/{bin_mc_sig}.parquet.gzip')
            # df_mc_bkg = pd.read_parquet(f'df_test/{bin_mc_bkg}.parquet.gzip')
            # df_mc_bkg = df_mc_bkg.query(f'(mass < 1.105 or mass > 1.13) and ct > {ct_bins[0]} and ct < {ct_bins[1]} and mass > 1.09 and mass < 1.15 and centrality > {low_cent} and centrality < {up_cent}')
            # df_mc_bkg_ = pd.read_parquet(f'df_test/{bin}.parquet.gzip')
            # df_mc_bkg_ = df_mc_bkg_.query(f'(mass < 1.105 or mass > 1.13) and ct > {ct_bins[0]} and ct < {ct_bins[1]} and mass > 1.09 and mass < 1.15 and centrality > {low_cent} and centrality < {up_cent}')
            # df_mc_sig = df_mc_sig.query(f"ct > {ct_bins[0]} and ct < {ct_bins[1]} and mass > 1.09 and mass < 1.15")
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
            # df_mc_ = [df_data_sidebands, df_mc_sig] #, df_mc_bkg, df_mc_bkg_]
            # df_mc = pd.concat(df_mc_)

            df_mc_1 = pd.read_parquet(f"df_test/test_data_{cent_bins[0]}_{cent_bins[1]}_{ct_[0]}_{ct_[1]}_predict.parquet.gzip")
            # df_mc_2 = pd.read_parquet(f"df_test/train_data_{cent_bins[0]}_{cent_bins[1]}_{ct_[0]}_{ct_[1]}_predict.parquet.gzip")
            # df_data_sidebands.rename(columns={'model_output_0':'model_output_background'})
            df_mc__ = [df_mc_1]
            df_mc_ = pd.concat(df_mc__)
            df_mc_ = df_mc_.query(f'ct > {ct_bins[0]} and ct < {ct_bins[1]} and mass > 1.09 and mass < 1.15 and pt > 0.5 and pt < 3.5')

            eff_selected = np.arange(0.1, MAX_EFF, 0.1)
            y_truth_tmp = df_mc_['flag'].apply(lambda x : x if x == 1 else 0)
            score = analysis_utils.score_from_efficiency_array(
                y_truth_tmp, df_mc_['model_output_background'], eff_selected, keep_lower=True)
            # print(y_truth_tmp)
            # print(score)

            # del df_data_sidebands, df_mc_

            for bdt_score, bdt_eff in zip(score, eff_selected):
                if bdt_eff < 0.49 or bdt_eff > 0.51: #or bdt_score > 0.15:
                    continue
                print(f'processing {bin}: bkg cut = {bdt_score:.4f}')

                
                for i_split, split in enumerate(SPLIT_LIST):
                    split_ineq_sign = '> -0.1'
                    if SPLIT:
                        split_ineq_sign = '> 0.5'
                        if split == 'antimatter':
                            split_ineq_sign = '< 0.5'
                    # apply cut
                    df_mc_cut = df_mc_.query(f"model_output_background < 0.3")
                    print(f'size of df = {df_mc_cut.shape[0]}')
                    df_mc_cut_bkg = df_data_sidebands.query(f"model_output_0 < 0.3")
                    df_data_cut = df_data_.query(f"matter {split_ineq_sign} and model_output_0 < 0.3")
                    eff = (df_mc_cut.query(f"flag==1").shape[0])/(df_mc_.query(f"flag==1").shape[0])
                    # plt.hist(df_data_cut["model_output_1"],bins=1000)
                    # plt.hist(df_mc_cut_bkg["model_output_1"],bins=1000)
                    # plt.hist(df_mc_cut["model_output_non_prompt"],bins=1000)
                    # plt.yscale("log")
                    # plt.savefig("cosPA_data.png")
                    # plt.close()

                    # df_mc_no_bdt = df_mc.query(f"ct > {ct_bins[0]} and ct < {ct_bins[1]} and mass > 1.09 and mass < 1.15")

                    # get invariant mass
                    df_mc_cut["model_output_2"] = df_mc_cut["model_output_prompt"]
                    # df_data_cut = df_mc_cut.sample(frac=0.5)
                    # df_data_cut = df_data_cut.query(f"pdg {split_ineq_sign}")
                    data_inv_mass = df_data_cut["mass"]
                    # mc_inv_mass_bkg = df_mc_cut.query("y_true_from_flags==0")["mass"]
                    # mc_inv_mass_sig = df_mc_cut.query("y_true_from_flags==2 or y_true_from_flags==1")["mass"]

                    # plt.hist(data_inv_mass, bins=1000)
                    # plt.savefig("plt.pdf")

                    # get non-prompt bdt score
                    data_bdt_non_prompt_out = df_data_cut["model_output_2"]#["model_output_2"]
                    mc_prompt_bdt_non_prompt_out = df_mc_cut.query("flag == 1")["model_output_prompt"] # .drop(df_data_cut.index)
                    # mc_prompt_bdt_non_prompt_out_uncut = df_mc_no_bdt.query("y_true_from_flags==2")["model_output_prompt"]
                    # bdt_eff_corrected = mc_prompt_bdt_non_prompt_out.shape[0]/mc_prompt_bdt_non_prompt_out_uncut.shape[0]
                    # print(f"bdt_eff = {bdt_eff}")

                    # mc_non_prompt_bdt_non_prompt_out = df_mc_cut.drop(df_data_cut.index).query("(flag == 2) and (flag != 1 or flag != 0)")["model_output_prompt"]
                    mc_non_prompt_bdt_non_prompt_out = df_mc_cut.query("flag == 2 or flag == 4")["model_output_prompt"]
                    mc_background_bdt_non_prompt_out = df_mc_cut_bkg["model_output_2"]#["model_output_2"]

                    # mc_count_p = mc_prompt_bdt_non_prompt_out.count()
                    # mc_count_np = mc_non_prompt_bdt_non_prompt_out.count()

                    for i_bkg_func, bkg_function in enumerate(['pol','expo']):
                        # fit to invariant mass
                        inv_mass = ROOT.RooRealVar("m","#it{M} (p + #pi^{-})",1.09,1.15,"GeV/#it{c}^{2}")
                        inv_mass.setBins(100)
                        inv_mass_roo_data = helpers.ndarray2roo(data_inv_mass.to_numpy(),inv_mass)
                        inv_mass_data = ROOT.RooDataHist("db_m","db_m",ROOT.RooArgList(inv_mass),inv_mass_roo_data)
                        # inv_mass_roo_mc = helpers.ndarray2roo(mc_inv_mass_sig.to_numpy(),inv_mass)
                        # inv_mass_mc = ROOT.RooDataHist("db_m_mc","db_m_mc",ROOT.RooArgList(inv_mass),inv_mass_roo_mc)
                        mass = ROOT.RooRealVar('#it{m}_{#Lambda}','mass',1.113,1.12,"GeV/#it{c}^{2}")
                        sigma = ROOT.RooRealVar('#sigma','sigma',0.001,0.003,"GeV/#it{c}^{2}")
                        alpha_left = ROOT.RooRealVar('#alpha_{left}','alpha_left',0.9,5.)
                        alpha_right = ROOT.RooRealVar('#alpha_{right}','alpha_right',0.9,5.)
                        n_left = ROOT.RooRealVar('n_{left}','n_left',4.,50.)
                        n_right = ROOT.RooRealVar('n_{right}','n_right',4.,50.)
                        par_a = ROOT.RooRealVar('a','a',-10.,10.,"#it{c}^{2}/GeV")
                        par_b = ROOT.RooRealVar('b','b',-10.,10.,"#it{c}^{4}/GeV^{2}")
                        slope = ROOT.RooRealVar('slope','slope',-30.,0.,"#it{c}^{2}/GeV")
                        #signal_pdf = ROOT.RooGaussian('signal','signal',inv_mass,mass,sigma)
                        signal_pdf = ROOT.RooDSCBShape('signal','signal',inv_mass,mass,sigma,alpha_left,n_left,alpha_right,n_right)
                        bkg_pdf = ROOT.RooPolynomial('bkg','bkg',inv_mass,ROOT.RooArgList(par_a,par_b))
                        if bkg_function == 'expo':
                            bkg_pdf = ROOT.RooExponential('bkg','bkg',inv_mass,slope)
                        f_signal = ROOT.RooRealVar("#it{f}_{signal}","f_signal",0.,1.)
                        n_tot = ROOT.RooRealVar("#it{N}_{tot}","n_tot",0.,1.e9)
                        model_mass_ = ROOT.RooAddPdf("model_mass_","model_mass_",signal_pdf,bkg_pdf,f_signal)
                        model_mass = ROOT.RooAddPdf("model_mass","model_mass",ROOT.RooArgList(model_mass_),ROOT.RooArgList(n_tot))
                        model_mass.fitTo(inv_mass_data)
                        r = model_mass.fitTo(inv_mass_data,ROOT.RooFit.Save())
                        if r.covQual() < 1.5:
                            continue

                        # fit to bdt output
                        bdt_out = ROOT.RooRealVar("BDT out","BDT out",0.,1)
                        bdt_out.setBins(50)
                        bdt_roo_data = helpers.ndarray2roo(data_bdt_non_prompt_out.to_numpy(), bdt_out)
                        bdt_data = ROOT.RooDataHist("dh", "dh", ROOT.RooArgList(bdt_out), bdt_roo_data)
                        bdt_roo_mc_prompt = helpers.ndarray2roo(mc_prompt_bdt_non_prompt_out.to_numpy(), bdt_out)
                        bdt_mc_prompt = ROOT.RooDataHist("dp", "dp", ROOT.RooArgList(bdt_out), bdt_roo_mc_prompt)
                        bdt_roo_mc_non_prompt = helpers.ndarray2roo(mc_non_prompt_bdt_non_prompt_out.to_numpy(), bdt_out)
                        bdt_mc_non_prompt = ROOT.RooDataHist("dnp", "dnp", ROOT.RooArgList(bdt_out), bdt_roo_mc_non_prompt)
                        bdt_roo_mc_bkg = helpers.ndarray2roo(mc_background_bdt_non_prompt_out.to_numpy(), bdt_out)
                        bdt_mc_bkg = ROOT.RooDataHist("db", "db", ROOT.RooArgList(bdt_out), bdt_roo_mc_bkg)

                        h_data = bdt_data.createHistogram("BDT out")
                        h_p = bdt_mc_prompt.createHistogram("BDT out")
                        h_np = bdt_mc_non_prompt.createHistogram("BDT out")
                        h_b = bdt_mc_bkg.createHistogram("BDT out")

                        # TFF fit
                        mc = ROOT.TObjArray(3)
                        mc.Add(h_p)
                        mc.Add(h_np)
                        mc.Add(h_b)
                        fit = ROOT.TFractionFitter(h_data,mc)
                        fit.Constrain(0,0,1)
                        fit.Constrain(1,0,1)
                        fit.Constrain(2,0,1)
                        fitter = fit.GetFitter()
                        fitter.Config().ParSettings(2).SetValue(1-f_signal.getVal())
                        fitter.Config().ParSettings(2).Fix()
                        status = fit.Fit()
                        result = ROOT.TH1F()
                        h_fit_pred = []
                        h_fit_resl = []
                        int_data = h_data.Integral()
                        if status==0:
                            h_result = fit.GetPlot()
                            result = ROOT.TH1F(h_result)
                            for i_mc in range(3):
                                h_fit_resl.append([c_double(),c_double()])
                                fit.GetResult(i_mc,h_fit_resl[i_mc][0],h_fit_resl[i_mc][1])
                                h_fit_pred.append(fit.GetMCPrediction(i_mc))
                                h_fit_pred[i_mc].SetLineColor(colors[i_mc])
                                h_fit_pred[i_mc].Scale(h_fit_resl[i_mc][0].value*int_data/h_fit_pred[i_mc].Integral())
                        else:
                            continue

                        # mass plot
                        frame_mass = inv_mass.frame(ROOT.RooFit.Name(f"fMass_{ct_bins[0]}_{ct_bins[1]}_{bdt_eff:.2f}"))
                        frame_mass.SetTitle(f"BDT eff = {eff:.2f}")
                        inv_mass_data.plotOn(frame_mass,ROOT.RooFit.Name("inv_mass_data"))
                        model_mass.plotOn(frame_mass,ROOT.RooFit.Components("bkg"),ROOT.RooFit.Name('background'),ROOT.RooFit.LineStyle(ROOT.kDashed),ROOT.RooFit.LineColor(ROOT.kGreen))
                        model_mass.plotOn(frame_mass,ROOT.RooFit.Components("signal"),ROOT.RooFit.Name('signal'),ROOT.RooFit.LineStyle(ROOT.kDashed),ROOT.RooFit.LineColor(ROOT.kRed))
                        model_mass.plotOn(frame_mass,ROOT.RooFit.Name('model_mass'))
                        model_mass.paramOn(frame_mass, ROOT.RooFit.Label("#chi^{2}/NDF = " + "{:.2f}".format(frame_mass.chiSquare("model_mass", "inv_mass_data"))), ROOT.RooFit.Layout(0.546084,0.867801,0.883868))
                
                        # bdt plot


                        f = ROOT.TFile(f"fit_9_noRadius_nodcaV0tracksCut_lambda_splitCentInMC_saveHistos_{cent_bins[0]}_{cent_bins[1]}.root","update")
                        f.cd()
                        f.mkdir(f'{split}_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}_{bkg_function}')
                        f.cd(f'{split}_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}_{bkg_function}')
                        frame_mass.getAttText().SetTextFont(44)
                        frame_mass.getAttText().SetTextSize(15)
                        frame_mass.getAttLine().SetLineWidth(0)
                        frame_mass.Write()
                        h_p.Write()
                        h_np.Write()
                        h_b.Write()
                        h_data.Write()
                        c_sim_fit = ROOT.TCanvas(f"cSimFit_{ct_bins[0]}_{ct_bins[1]}_{bdt_eff:.2f}",f"cSimFit_{ct_bins[0]}_{ct_bins[1]}_{bdt_eff:.2f}",1200,550)
                        c_sim_fit.Divide(2,1)
                        c_sim_fit.cd(1)
                        c_sim_fit.GetPad(1).SetLogy()
                        h_data.Draw()
                        result.Draw("same")
                        for i_mc in range(3):
                            h_fit_pred[i_mc].Draw("same")
                        c_sim_fit.cd(2)
                        frame_mass.Draw()
                        ROOT.gPad.Modified()
                        ROOT.gPad.Update()

                        if not os.path.isdir(f"plots/signal_extraction_/{cent_bins[0]}_{cent_bins[1]}/{ct_bins[0]}_{ct_bins[1]}"):
                            os.mkdir(f"plots/signal_extraction_/{cent_bins[0]}_{cent_bins[1]}/{ct_bins[0]}_{ct_bins[1]}")
                        c_sim_fit.Print(f"plots/signal_extraction_/{cent_bins[0]}_{cent_bins[1]}/{ct_bins[0]}_{ct_bins[1]}/{split}_{bkg_function}_{bdt_eff:.2f}.pdf")
                        c_sim_fit.Write()

                        f.Close()
                        scale_frac = h_fit_res[0][0].value
                        h_raw_yields[i_split][i_bkg_func].SetBinContent(h_raw_yields[i_split][i_bkg_func].FindBin(bdt_eff+0.0001),n_tot.getVal()*scale_frac/eff)
                        h_raw_yields[i_split][i_bkg_func].SetBinError(h_raw_yields[i_split][i_bkg_func].FindBin(bdt_eff+0.0001),n_tot.getVal()*h_fit_resl[0][1].value/eff)
                        print(f'************* fraction result = {h_fit_resl[0][0]} *************')

            f = ROOT.TFile(f"fit_9_noRadius_nodcaV0tracksCut_lambda_splitCentInMC_saveHistos_{cent_bins[0]}_{cent_bins[1]}.root","update")
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