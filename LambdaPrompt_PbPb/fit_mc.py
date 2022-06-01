import sys
sys.argv.append( '-b-' )
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

SPLIT = False
use_root = True

colors = [ROOT.kBlue, ROOT.kRed, ROOT.kBlack,ROOT.kOrange, ROOT.kGreen, ROOT.kGray+2]
ROOT.gInterpreter.ProcessLine("#include \"../utils/RooDSCBShape.h\"")
ROOT.gInterpreter.ProcessLine(".L ../utils/RooDSCBShape.cxx+")

##################################################################
# read configuration file
##################################################################
config = 'config.yaml'
with open(os.path.expandvars(config), 'r') as stream:
    try:
        params = yaml.full_load(stream)
    except yaml.YAMLError as exc:
        print(exc)

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

f = ROOT.TFile("out.root","recreate")

for split in SPLIT_LIST:
    split_ineq_sign = '> -0.1'
    if SPLIT:
        split_ineq_sign = '> 0.5'
        if split == 'antimatter':
            split_ineq_sign = '< 0.5'

    for i_cent_bins in range(len(CENTRALITY_LIST)):
        cent_bins = CENTRALITY_LIST[i_cent_bins]
        for ct_bins in zip(CT_BINS_CENT[i_cent_bins][:-1], CT_BINS_CENT[i_cent_bins][1:]):
            if ct_bins[0] < 1 or ct_bins[1] > 20:
                 continue
            bdt_cut = 1
            bin = f'{split}_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}'
            print(bin)
            bin_mc = f'test_data_{cent_bins[0]}_{cent_bins[1]}_{ct_bins[0]}_{ct_bins[1]}'
            df_data = pd.read_parquet(f'df/{bin_mc}.parquet.gzip')
            df_mc = pd.read_parquet(f'df/{bin_mc}.parquet.gzip')

            # apply cut
            df_data_cut = df_data.query(f"model_output_background < {bdt_cut} and model_output_prompt > 0 and mass > 1.09 and mass < 1.15")
            df_mc_cut = df_mc.query(f"model_output_background < {bdt_cut} and model_output_prompt > 0 and mass > 1.09 and mass < 1.15")

            # get invariant mass
            data_inv_mass = df_data_cut["mass"]
            mc_inv_mass_bkg = df_mc_cut.query("y_true==0")["mass"]

            plt.hist(data_inv_mass, bins=1000)
            plt.savefig("plt.pdf")
            plt.close()

            # get non-prompt bdt score
            data_bdt_non_prompt_out = df_data_cut["model_output_non_prompt"]
            mc_prompt_bdt_non_prompt_out = df_mc_cut.query("y_true==2")["model_output_non_prompt"]
            mc_non_prompt_bdt_non_prompt_out = df_mc_cut.query("y_true==1")["model_output_non_prompt"]
            mc_background_bdt_non_prompt_out = df_mc_cut.query("y_true==0")["model_output_non_prompt"]

            mc_count_p = mc_prompt_bdt_non_prompt_out.count()
            mc_count_np = mc_non_prompt_bdt_non_prompt_out.count()

            # fit to invariant mass
            inv_mass = ROOT.RooRealVar("m","m",1.09,1.15)
            inv_mass_roo_data = helpers.ndarray2roo(data_inv_mass.to_numpy(),inv_mass)
            inv_mass_data = ROOT.RooDataHist("db_m","db_m",ROOT.RooArgList(inv_mass),inv_mass_roo_data)
            inv_mass_roo_bkg = helpers.ndarray2roo(mc_inv_mass_bkg.to_numpy(),inv_mass)
            inv_mass_bkg = ROOT.RooDataHist("db_m","db_m",ROOT.RooArgList(inv_mass),inv_mass_roo_bkg)
            inv_mass_bkg_pdf = ROOT.RooHistPdf("db_m_pdf","db_m_pdf",ROOT.RooArgList(inv_mass),inv_mass_bkg)
            mass = ROOT.RooRealVar('mass','mass',1.11,1.12)
            sigma = ROOT.RooRealVar('sigma','sigma',0.,0.005)
            alpha_left = ROOT.RooRealVar('alpha_left','alpha_left',0.,2.)
            alpha_right = ROOT.RooRealVar('alpha_right','alpha_right',0.,2.)
            n_left = ROOT.RooRealVar('n_left','n_left',0.,15.)
            n_right = ROOT.RooRealVar('n_right','n_right',0.,15.)
            signal_pdf = ROOT.RooDSCBShape('signal','signal',inv_mass,mass,sigma,alpha_left,n_left,alpha_right,n_right)
            w_signal = ROOT.RooRealVar("w_signal","w_signal",0.,1.)
            model_mass = ROOT.RooAddPdf("model_mass","model_mass",signal_pdf,inv_mass_bkg_pdf,w_signal)

            # fit to bdt output
            bdt_out = ROOT.RooRealVar("BDT out","BDT out",0.,1.)
            bdt_out.setBins(800)
            bdt_roo_data = helpers.ndarray2roo(data_bdt_non_prompt_out.to_numpy(), bdt_out)
            bdt_data = ROOT.RooDataHist("dh", "dh", ROOT.RooArgList(bdt_out), bdt_roo_data)
            bdt_roo_mc_prompt = helpers.ndarray2roo(mc_prompt_bdt_non_prompt_out.to_numpy(), bdt_out)
            bdt_mc_prompt = ROOT.RooDataHist("dp", "dp", ROOT.RooArgList(bdt_out), bdt_roo_mc_prompt)
            bdt_mc_prompt_pdf = ROOT.RooHistPdf("dppdf", "dppdf", ROOT.RooArgList(bdt_out), bdt_mc_prompt)
            
            h = bdt_mc_prompt_pdf.createHistogram("BDT out")
            h.GetXaxis().SetRangeUser(0.,0.25)
            max_bdt_prompt = h.GetBinCenter(h.GetMaximumBin())
            h.GetXaxis().SetRangeUser(0.,1.)
            gaus_1 = ROOT.TF1("gaus_1","gaus",0,1)
            gaus_2 = ROOT.TF1("gaus_2","gaus",0,1)
            h.Fit("gaus_1","QR","",max_bdt_prompt-0.005,max_bdt_prompt+0.005)
            mu = h.GetFunction("gaus_1").GetParameter(1)
            sig = h.GetFunction("gaus_1").GetParameter(2)
            print(f"mu_1 = {mu}; sig_1 = {sig}")
            bdt_score_region = [[mu-1*sig,mu+1*sig],[mu+1*sig,mu+3*sig],[mu+3*sig,mu+5*sig],[mu+5*sig,mu+7*sig],[mu+7*sig,mu+9*sig],[mu+9*sig,mu+11*sig]]

            bdt_roo_mc_non_prompt = helpers.ndarray2roo(mc_non_prompt_bdt_non_prompt_out.to_numpy(), bdt_out)
            bdt_mc_non_prompt = ROOT.RooDataHist("dnp", "dnp", ROOT.RooArgList(bdt_out), bdt_roo_mc_non_prompt)
            bdt_mc_non_prompt_pdf = ROOT.RooHistPdf("dnppdf", "dnppdf", ROOT.RooArgList(bdt_out), bdt_mc_non_prompt)
            bdt_roo_mc_bkg = helpers.ndarray2roo(mc_background_bdt_non_prompt_out.to_numpy(), bdt_out)
            bdt_mc_bkg = ROOT.RooDataHist("db", "db", ROOT.RooArgList(bdt_out), bdt_roo_mc_bkg)
            bdt_mc_bkg_pdf = ROOT.RooHistPdf("dbpdf", "dbpdf", ROOT.RooArgList(bdt_out), bdt_mc_bkg)
            frame = bdt_out.frame()
            w_non_prompt = ROOT.RooRealVar("w_non_prompt","w_non_prompt",0.,1.)
            model_signal = ROOT.RooAddPdf("model_signal","model_signal",bdt_mc_non_prompt_pdf,bdt_mc_prompt_pdf,w_non_prompt)
            model = ROOT.RooAddPdf("model","model",model_signal,bdt_mc_bkg_pdf,w_signal)
            
            # simultaneous fit
            sample = ROOT.RooCategory("sample","sample")
            sample.defineType("bdt")
            sample.defineType("mass")
            combData = ROOT.RooDataHist("combData","combData",ROOT.RooArgSet(bdt_out,inv_mass),ROOT.RooFit.Index(sample),ROOT.RooFit.Import("bdt",bdt_data),ROOT.RooFit.Import("mass",inv_mass_data))
            simPdf = ROOT.RooSimultaneous("simPdf","simPdf",sample)
            simPdf.addPdf(model,"bdt")
            simPdf.addPdf(model_mass,"mass")
            simPdf.fitTo(combData)

            # mass plot
            frame_mass = inv_mass.frame()
            inv_mass_data.plotOn(frame_mass)
            model_mass.plotOn(frame_mass,ROOT.RooFit.Components("db_m_pdf"),ROOT.RooFit.Name('background'),ROOT.RooFit.LineStyle(ROOT.kDashed),ROOT.RooFit.LineColor(ROOT.kGreen))
            model_mass.plotOn(frame_mass,ROOT.RooFit.Components("signal"),ROOT.RooFit.Name('signal'),ROOT.RooFit.LineStyle(ROOT.kDashed),ROOT.RooFit.LineColor(ROOT.kRed))
            model_mass.plotOn(frame_mass)

            # bdt plot
            bdt_data.plotOn(frame)
            model.plotOn(frame,ROOT.RooFit.Components('dbpdf'),ROOT.RooFit.Name('background'),ROOT.RooFit.LineStyle(ROOT.kDashed),ROOT.RooFit.LineColor(ROOT.kGreen))
            model.plotOn(frame,ROOT.RooFit.Components('dnppdf'),ROOT.RooFit.Name('non-prompt'),ROOT.RooFit.LineStyle(ROOT.kDashed),ROOT.RooFit.LineColor(ROOT.kOrange))
            model.plotOn(frame,ROOT.RooFit.Components('dppdf'),ROOT.RooFit.Name('prompt'),ROOT.RooFit.LineStyle(ROOT.kDashed),ROOT.RooFit.LineColor(ROOT.kRed))
            model.plotOn(frame,ROOT.RooFit.Name('model'),ROOT.RooFit.LineColor(ROOT.kBlue))

            # features
            f.mkdir(bin)
            features = ['cosPA','dcaV0tracks','dcaPiPV','dcaPrPV','dcaV0PV','tpcNsigmaPr','radius']
            regions_labels = ["main","second"]
            for i_f, ft in enumerate(features):
                d = []
                plot = ROOT.TCanvas(f"fCanv_{ft}",f"fCanv_{ft}")
                h_features = []
                for i_r, ranges in enumerate(bdt_score_region):
                    d.append(df_mc_cut.query(f"y_true==2 and model_output_non_prompt > {ranges[0]} and model_output_non_prompt < {ranges[1]}")[ft].to_numpy())
                    mn = d[i_r].min()
                    mx = d[i_r].max()
                    h_features.append(ROOT.TH1D(f"fHist_{ft}",f";{ft};Entries",400,mn,mx))
                    for i_l in np.arange(d[i_r].size):
                        h_features[i_r].Fill(d[i_r][i_l])
                    f.cd(bin)
                    h_features[i_r].SetLineColor(colors[i_r])
                    h_features[i_r].SetMarkerColor(colors[i_r])
                    h_features[i_r].SetMarkerStyle(20)
                    h_features[i_r].SetMarkerSize(1)
                    h_features[i_r].Scale(1/h_features[i_r].Integral())
                    h_features[i_r].SetDrawOption("LP")
                    h_features[i_r].Write()
                    plot.cd()
                    if i_r == 0:
                        h_features[i_r].Draw("lpex0")
                    else:
                        h_features[i_r].Draw("lpex0same")
                f.cd(bin)
                plot.SetLogy()
                plot.Write()
                #plot.Print(f"{plot.GetName()}.png")
                plt.hist(d,bins=40,alpha=0.5,density=True)
                plt.yscale('log')
                plt.savefig(f"h_{ft}")
                plt.close()

            f.cd(bin)
            frame.Write()
            h.Write()
            plot_2 = ROOT.TCanvas("c_2","c_2")
            h.Draw("")
            l = []
            for i_r, r in enumerate(bdt_score_region):
                l.append(ROOT.TLine(r[0],0,r[0],1))
                l[i_r].Draw("same")
            plot_2.Write()
            gaus_1.Write()
            gaus_2.Write()
            frame_mass.Write()

            print(f'Non-prompt / prompt = {mc_count_np/(mc_count_p+mc_count_np)}')
            print(f'N_prompt = {mc_count_p}')

f.Close()