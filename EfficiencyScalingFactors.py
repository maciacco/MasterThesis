#!/usr/bin/env python3

import numpy as np
import ROOT

data_path_proton = "data/Proton_PbPb"
data_path_he3 = "data/He3_PbPb"

input_file_eff_default_proton = ROOT.TFile.Open(f"{data_path_proton}/EfficiencyProtonMC_21l5_false_XS.root")
input_file_eff_increased_proton = ROOT.TFile.Open(f"{data_path_proton}/EfficiencyProtonMC_21l5_false_XS+.root")
input_file_eff_default_he3 = ROOT.TFile.Open(f"{data_path_he3}/EfficiencyHe3_XSDefault.root")
input_file_eff_increased_he3 = ROOT.TFile.Open(f"{data_path_he3}/EfficiencyHe3_XSPlus.root")

out_file = ROOT.TFile.Open("scaling_factors.root","recreate")

par_fit2g4_antip = 1.02829
par_fit2g4_antip_err = 0.0073243 
par_fit2g4_p = 0.993
par_fit2g4_p_err = 0.035
par_fit2g4_antihe3 = 0.83
par_fit2g4_antihe3_err = 0.07
par_fit2g4_he3 = 0.86836
par_fit2g4_he3_err = 0.0580437

out_file.cd() 
eff_default_antip = input_file_eff_default_proton.Get("_/fAEff_TOF_-10_0;1")
eff_increased_antip = input_file_eff_increased_proton.Get("_/fAEff_TOF_-10_0;1")
eff_difference_antip = ROOT.TH1D(eff_default_antip)
eff_difference_antip.Add(eff_increased_antip,-1.)
eff_difference_antip.Divide(eff_difference_antip,eff_default_antip,1,-0.5)
scaling_factor = ROOT.TH1D(eff_difference_antip)
for iB in range(scaling_factor.GetNbinsX()+1):
    if iB == 0:
        continue
    tmp = scaling_factor.GetBinContent(iB)
    tmp_err = scaling_factor.GetBinError(iB)
    scaling_factor.SetBinContent(iB,1.-tmp*(1.-par_fit2g4_antip))
    scaling_factor.SetBinError(iB,tmp_err*(1.-par_fit2g4_antip))
fit_func = ROOT.TF1("fit_func","[0]*TMath::Power(x,[1])",1.,2.,2)
fit_func.SetParLimits(0,0.,1.1)
scaling_factor.Fit(fit_func,"MR","",1.,2.)
scaling_factor.Write()

out_file.cd() 
eff_default_p = input_file_eff_default_proton.Get("_/fMEff_TOF_-10_0;1")
eff_increased_p = input_file_eff_increased_proton.Get("_/fMEff_TOF_-10_0;1")
eff_difference_p = ROOT.TH1D(eff_default_p)
eff_difference_p.Add(eff_increased_p,-1.)
eff_difference_p.Divide(eff_difference_p,eff_default_p,1,-0.5)
scaling_factor = ROOT.TH1D(eff_difference_p)
for iB in range(scaling_factor.GetNbinsX()+1):
    if iB == 0:
        continue
    tmp = scaling_factor.GetBinContent(iB)
    tmp_err = scaling_factor.GetBinError(iB)
    scaling_factor.SetBinContent(iB,1.-tmp*(1.-par_fit2g4_p))
    scaling_factor.SetBinError(iB,tmp_err*(1.-par_fit2g4_p))
fit_func = ROOT.TF1("fit_func_1","[0]*TMath::Power(x,[1])",1.,2.,2)
fit_func.SetParLimits(0,0.,1.1)
scaling_factor.Fit(fit_func,"MR","",1.,2.)
scaling_factor.Write()

out_file.cd() 
eff_default_he3 = input_file_eff_default_he3.Get("fMEff_TPC_0_5;1")
eff_increased_he3 = input_file_eff_increased_he3.Get("fMEff_TPC_0_5;1")
eff_difference_he3 = ROOT.TH1D(eff_default_he3)
eff_difference_he3.Add(eff_increased_he3,-1.)
eff_difference_he3.Divide(eff_difference_he3,eff_default_he3,1,-0.5)
scaling_factor = ROOT.TH1D(eff_difference_he3)
for iB in range(scaling_factor.GetNbinsX()+1):
    if iB == 0:
        continue
    tmp = scaling_factor.GetBinContent(iB)
    tmp_err = scaling_factor.GetBinError(iB)
    scaling_factor.SetBinContent(iB,1.-tmp*(1.-par_fit2g4_he3))
    scaling_factor.SetBinError(iB,tmp_err*(1.-par_fit2g4_he3))
fit_func = ROOT.TF1("fit_func_2","[0]*TMath::Power(x,[1])",1.,10.,2)
fit_func.SetParLimits(0,0.,1.1)
scaling_factor.Fit(fit_func,"MR","",1.,10.)
scaling_factor.Write()

out_file.cd() 
eff_default_antihe3 = input_file_eff_default_he3.Get("fAEff_TPC_0_5;1")
eff_increased_antihe3 = input_file_eff_increased_he3.Get("fAEff_TPC_0_5;1")
eff_difference_antihe3 = ROOT.TH1D(eff_default_antihe3)
eff_difference_antihe3.Add(eff_increased_antihe3,-1.)
eff_difference_antihe3.Divide(eff_difference_antihe3,eff_default_antihe3,1,-0.5)
scaling_factor = ROOT.TH1D(eff_difference_antihe3)
for iB in range(scaling_factor.GetNbinsX()+1):
    if iB == 0:
        continue
    tmp = scaling_factor.GetBinContent(iB)
    tmp_err = scaling_factor.GetBinError(iB)
    scaling_factor.SetBinContent(iB,1.-tmp*(1.-par_fit2g4_antihe3))
    scaling_factor.SetBinError(iB,tmp_err*(1.-par_fit2g4_antihe3))
fit_func = ROOT.TF1("fit_func_3","[0]*TMath::Power(x,[1])",1.,10.,2)
fit_func.SetParLimits(0,0.,1.1)
scaling_factor.Fit(fit_func,"MR","",1.,10.)
scaling_factor.Write()

out_file.Close()