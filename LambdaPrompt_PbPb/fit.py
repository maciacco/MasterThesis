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

SPLIT = False

ROOT.gInterpreter.ProcessLine("#include \"../utils/RooDSCBShape.h\"")
ROOT.gInterpreter.ProcessLine(".L ../utils/RooDSCBShape.cxx+")
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

f = ROOT.TFile("out_1.root","recreate")

score_eff_dict = pickle.load(open('file_score_eff_dict','rb'))

for split in SPLIT_LIST:
    split_ineq_sign = '> -0.1'
    if SPLIT:
        split_ineq_sign = '> 0.5'
        if split == 'antimatter':
            split_ineq_sign = '< 0.5'

    for i_cent_bins in range(len(CENTRALITY_LIST)):
        cent_bins = CENTRALITY_LIST[i_cent_bins]
        for ct_bins in zip(CT_BINS_CENT[i_cent_bins][:-1], CT_BINS_CENT[i_cent_bins][1:]):
            if ct_bins[0] < 5 or ct_bins[1] > 5.5:
                continue
            h_raw_yields = ROOT.TH1D("fRawYields","fRawYields",100,0,100)
            for bdt_eff in range(61,82):
                bdt_cut = 1-bdt_eff/100.
                bin = f'{split}_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}'
                bin_mc = f'test_data_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}'
                df_data = uproot.open(os.path.expandvars("../data/LambdaPrompt_PbPb/outTree.root"))['LambdaTreeBDTOut'].arrays(library="pd")
                df_mc = pd.read_parquet(f'df/df_alidbz/df/{bin_mc}.parquet.gzip')

                # apply cut
                df_data_cut = df_data.query(f"bdtOutputPrompt > {bdt_cut} and ct > {ct_bins[0]} and ct < {ct_bins[1]} and mass > 1.09 and mass < 1.15")
                df_mc_cut = df_mc.query(f"model_output_prompt > {bdt_cut} and mass > 1.09 and mass < 1.15")

                # get invariant mass
                data_inv_mass = df_data_cut["mass"]
                mc_inv_mass_bkg = df_mc_cut.query("y_true==0")["mass"]
                mc_inv_mass_sig = df_mc_cut.query("y_true==2 or y_true==1")["mass"]

                plt.hist(data_inv_mass, bins=1000)
                plt.savefig("plt.pdf")

                # get non-prompt bdt score
                data_bdt_non_prompt_out = df_data_cut["bdtOutputNonPrompt"]
                mc_prompt_bdt_non_prompt_out = df_mc_cut.query("y_true==2")["model_output_non_prompt"]
                mc_non_prompt_bdt_non_prompt_out = df_mc_cut.query("y_true==1")["model_output_non_prompt"]
                mc_background_bdt_non_prompt_out = df_mc_cut.query("y_true==0")["model_output_non_prompt"]

                mc_count_p = mc_prompt_bdt_non_prompt_out.count()
                mc_count_np = mc_non_prompt_bdt_non_prompt_out.count()

                # fit to invariant mass
                inv_mass = ROOT.RooRealVar("m","m",1.09,1.15)
                inv_mass.setBins(200)
                inv_mass_roo_data = helpers.ndarray2roo(data_inv_mass.to_numpy(),inv_mass)
                inv_mass_data = ROOT.RooDataHist("db_m","db_m",ROOT.RooArgList(inv_mass),inv_mass_roo_data)
                inv_mass_roo_mc = helpers.ndarray2roo(mc_inv_mass_sig.to_numpy(),inv_mass)
                inv_mass_mc = ROOT.RooDataHist("db_m_mc","db_m_mc",ROOT.RooArgList(inv_mass),inv_mass_roo_mc)
                mass = ROOT.RooRealVar('mass','mass',1.112,1.118)
                sigma = ROOT.RooRealVar('sigma','sigma',0.,0.005)
                # alpha_left = ROOT.RooRealVar('alpha_left','alpha_left',0.,2.)
                # alpha_right = ROOT.RooRealVar('alpha_right','alpha_right',0.,2.)
                # n_left = ROOT.RooRealVar('n_left','n_left',0.,15.)
                # n_right = ROOT.RooRealVar('n_right','n_right',0.,15.)
                par_a = ROOT.RooRealVar('a','a',-10.,10.)
                par_b = ROOT.RooRealVar('b','b',-10.,10.)
                signal_pdf = ROOT.RooGaussian('signal','signal',inv_mass,mass,sigma)
                # signal_pdf = ROOT.RooDSCBShape('signal','signal',inv_mass,mass,sigma,alpha_left,n_left,alpha_right,n_right)
                signal_pdf.fitTo(inv_mass_mc)
                sigma.setConstant()
                # mass.setConstant()
                # alpha_left.setConstant()
                # alpha_right.setConstant()
                # n_left.setConstant()
                # n_right.setConstant()
                bkg_pdf = ROOT.RooPolynomial('bkg','bkg',inv_mass,ROOT.RooArgList(par_a,par_b))
                w_signal = ROOT.RooRealVar("w_signal","w_signal",0.,1.)
                n_tot = ROOT.RooRealVar("n_tot","n_tot",0.,5.e2)
                model_mass_ = ROOT.RooAddPdf("model_mass_","model_mass_",signal_pdf,bkg_pdf,w_signal)
                model_mass = ROOT.RooAddPdf("model_mass","model_mass",ROOT.RooArgList(model_mass_),ROOT.RooArgList(n_tot))

                # fit to bdt output
                bdt_out = ROOT.RooRealVar("BDT out","BDT out",0.,1.)
                bdt_out.setBins(100)
                bdt_roo_data = helpers.ndarray2roo(data_bdt_non_prompt_out.to_numpy(), bdt_out)
                bdt_data = ROOT.RooDataHist("dh", "dh", ROOT.RooArgList(bdt_out), bdt_roo_data)
                bdt_roo_mc_prompt = helpers.ndarray2roo(mc_prompt_bdt_non_prompt_out.to_numpy(), bdt_out)
                bdt_mc_prompt = ROOT.RooDataHist("dp", "dp", ROOT.RooArgList(bdt_out), bdt_roo_mc_prompt)
                bdt_mc_prompt_pdf = ROOT.RooHistPdf("dppdf", "dppdf", ROOT.RooArgList(bdt_out), bdt_mc_prompt)
                bdt_roo_mc_non_prompt = helpers.ndarray2roo(mc_non_prompt_bdt_non_prompt_out.to_numpy(), bdt_out)
                bdt_mc_non_prompt = ROOT.RooDataHist("dnp", "dnp", ROOT.RooArgList(bdt_out), bdt_roo_mc_non_prompt)
                bdt_mc_non_prompt_pdf = ROOT.RooHistPdf("dnppdf", "dnppdf", ROOT.RooArgList(bdt_out), bdt_mc_non_prompt)
                bdt_roo_mc_bkg = helpers.ndarray2roo(mc_background_bdt_non_prompt_out.to_numpy(), bdt_out)
                bdt_mc_bkg = ROOT.RooDataHist("db", "db", ROOT.RooArgList(bdt_out), bdt_roo_mc_bkg)
                bdt_mc_bkg_pdf = ROOT.RooHistPdf("dbpdf", "dbpdf", ROOT.RooArgList(bdt_out), bdt_mc_bkg)
                frame = bdt_out.frame(ROOT.RooFit.Name(f"fBDTOutNonPrompt_{ct_bins[0]}_{ct_bins[1]}"))
                w_non_prompt = ROOT.RooRealVar("w_non_prompt","w_non_prompt",0.,1.)
                model_signal = ROOT.RooAddPdf("model_signal","model_signal",bdt_mc_non_prompt_pdf,bdt_mc_prompt_pdf,w_non_prompt)
                model = ROOT.RooAddPdf("model","model",model_signal,bdt_mc_bkg_pdf,w_signal)

                # simultaneous fit
                ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)
                ROOT.RooMsgService.instance().setSilentMode(ROOT.kTRUE)
                ROOT.gErrorIgnoreLevel = ROOT.kError
                sample = ROOT.RooCategory("sample","sample")
                sample.defineType("bdt")
                sample.defineType("mass")
                combData = ROOT.RooDataHist("combData","combData",ROOT.RooArgSet(bdt_out,inv_mass),ROOT.RooFit.Index(sample),ROOT.RooFit.Import("bdt",bdt_data),ROOT.RooFit.Import("mass",inv_mass_data))
                simPdf = ROOT.RooSimultaneous("simPdf","simPdf",sample)
                simPdf.addPdf(model,"bdt")
                simPdf.addPdf(model_mass,"mass")
                simPdf.fitTo(combData)

                # mass plot
                frame_mass = inv_mass.frame(ROOT.RooFit.Name(f"fMass_{ct_bins[0]}_{ct_bins[1]}"))
                inv_mass_data.plotOn(frame_mass)
                model_mass.plotOn(frame_mass,ROOT.RooFit.Components("bkg"),ROOT.RooFit.Name('background'),ROOT.RooFit.LineStyle(ROOT.kDashed),ROOT.RooFit.LineColor(ROOT.kGreen))
                model_mass.plotOn(frame_mass,ROOT.RooFit.Components("signal"),ROOT.RooFit.Name('signal'),ROOT.RooFit.LineStyle(ROOT.kDashed),ROOT.RooFit.LineColor(ROOT.kRed))
                model_mass.plotOn(frame_mass)

                # bdt plot
                bdt_data.plotOn(frame)
                model.plotOn(frame,ROOT.RooFit.Components('dbpdf'),ROOT.RooFit.Name('background'),ROOT.RooFit.LineStyle(ROOT.kDashed),ROOT.RooFit.LineColor(ROOT.kGreen))
                model.plotOn(frame,ROOT.RooFit.Components('dnppdf'),ROOT.RooFit.Name('non-prompt'),ROOT.RooFit.LineStyle(ROOT.kDashed),ROOT.RooFit.LineColor(ROOT.kOrange))
                model.plotOn(frame,ROOT.RooFit.Components('dppdf'),ROOT.RooFit.Name('prompt'),ROOT.RooFit.LineStyle(ROOT.kDashed),ROOT.RooFit.LineColor(ROOT.kRed))
                model.plotOn(frame,ROOT.RooFit.Name('model'),ROOT.RooFit.LineColor(ROOT.kBlue))

                f.cd()
                f.mkdir(f'{split}_{cent_bins[0]}_{cent_bins[0]}_{ct_bins[0]}_{ct_bins[1]}')
                f.cd(f'{split}_{cent_bins[0]}_{cent_bins[0]}_{ct_bins[0]}_{ct_bins[1]}')
                frame.Write()
                frame_mass.Write()

                h_raw_yields.SetBinContent(int(bdt_eff),n_tot.getVal()*w_signal.getVal()*(1-w_non_prompt.getVal()))
                print(f'Non-prompt / prompt = {mc_count_np/(mc_count_p+mc_count_np)}')

            h_raw_yields.Write()