#!/usr/bin/env python3

import numpy as np
import ROOT

data_path_proton = "data/Proton_PbPb"
data_path_he3 = "data/He3_PbPb"
data_path_triton = "Triton_PbPb/out"

ROOT.gStyle.SetOptStat(00000000000)

input_file_eff_default_proton = ROOT.TFile.Open(f"{data_path_proton}/EfficiencyProtonMC_21l5_false_XS.root")
input_file_eff_increased_proton = ROOT.TFile.Open(f"{data_path_proton}/EfficiencyProtonMC_21l5_false_XS+.root")
input_file_eff_default_he3 = ROOT.TFile.Open(f"{data_path_he3}/EfficiencyHe3_XSDefault.root")
input_file_eff_increased_he3 = ROOT.TFile.Open(f"{data_path_he3}/EfficiencyHe3_XSPlus.root")
input_file_eff_default_triton = ROOT.TFile.Open(f"{data_path_triton}/EfficiencyTriton_XSDefault.root")
input_file_eff_increased_triton = ROOT.TFile.Open(f"{data_path_triton}/EfficiencyTriton_XSPlus.root")

out_file = ROOT.TFile.Open("scaling_factors.root","recreate")

par_fit2g4_antip = 1.02829
par_fit2g4_antip_err = 0.0073243 
par_fit2g4_p = 0.993
par_fit2g4_p_err = 0.035
par_fit2g4_antihe3 = 0.83
par_fit2g4_antihe3_err = 0.07
par_fit2g4_he3 = 0.86836
par_fit2g4_he3_err = 0.0580437
par_fit2g4_antitriton = 0.931474
par_fit2g4_antitriton_err = 0.07
par_fit2g4_triton = 0.86836
par_fit2g4_triton_err =  0.186728

# proton
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
    #scaling_factor.SetBinContent(iB,1.-tmp*(1.-par_fit2g4_antip))
    #scaling_factor.SetBinError(iB,tmp_err*(1.-par_fit2g4_antip))
    scaling_factor.SetBinContent(iB,-tmp*0.0856785)
    scaling_factor.SetBinError(iB,tmp_err*0.0856785)
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

# he3
out_file.cd() 
cHe3 = ROOT.TCanvas("cHe3","cHe3")
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
scaling_factor.SetTitle("^{3}He")
scaling_factor.GetYaxis().SetTitle("Scaling factor")
scaling_factor.Write()
scaling_factor.Draw("")
cHe3.Print("ScalingFactorHe3.pdf")

out_file.cd() 
cAntiHe3 = ROOT.TCanvas("cAntiHe3","cAntiHe3")
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
scaling_factor.SetTitle("^{3}#bar{He}")
scaling_factor.GetYaxis().SetTitle("Scaling factor")
scaling_factor.Write()
scaling_factor.Draw("")
cAntiHe3.Print("ScalingFactorAntiHe3.pdf")

# triton
out_file.cd() 
ctriton = ROOT.TCanvas("ctriton","ctriton")
eff_default_triton = input_file_eff_default_triton.Get("fMEff_TPC_-10_0;1")
eff_default_triton.Scale(1.0348)
eff_increased_triton = input_file_eff_increased_triton.Get("fMEff_TPC_-10_0;1")
eff_difference_triton = ROOT.TH1D(eff_default_triton)
eff_difference_triton.Add(eff_increased_triton,-1.)
eff_difference_triton.Divide(eff_difference_triton,eff_default_triton,1,(1-(1.5))/0.88)
scaling_factor = ROOT.TH1D(eff_difference_triton)
for iB in range(scaling_factor.GetNbinsX()+1):
    if iB == 0:
        continue
    tmp = scaling_factor.GetBinContent(iB)
    tmp_err = scaling_factor.GetBinError(iB)
    scaling_factor.SetBinContent(iB,1.-tmp*par_fit2g4_triton*(1.-par_fit2g4_triton))
    scaling_factor.SetBinError(iB,tmp_err*par_fit2g4_triton*(1.-par_fit2g4_triton))
    if iB == 1:
        scaling_factor.SetBinContent(iB,0)
        scaling_factor.SetBinError(iB,0)
fit_func = ROOT.TF1("fit_func_2","[0]*TMath::Power(x,[1])",1.,10.,2)
fit_func.SetParLimits(0,0.,1.1)
scaling_factor.Fit(fit_func,"MR","",1.6,3.)
scaling_factor.SetTitle("^{3}H")
scaling_factor.GetYaxis().SetTitle("Scaling factor")
scaling_factor.GetYaxis().SetRangeUser(1.,1.1)
scaling_factor.GetXaxis().SetRangeUser(1.6,3)
scaling_factor.Write()
scaling_factor.Draw("")
ctriton.Print("ScalingFactortriton.pdf")

out_file.cd() 
cAntitriton = ROOT.TCanvas("cAntitriton","cAntitriton")
eff_default_antitriton = input_file_eff_default_triton.Get("fAEff_TPC_-10_0;1")
eff_default_antitriton.Scale(1.06862)
eff_increased_antitriton = input_file_eff_increased_triton.Get("fAEff_TPC_-10_0;1")
eff_difference_antitriton = ROOT.TH1D(eff_default_antitriton)
eff_difference_antitriton.Add(eff_increased_antitriton,-1.)
eff_difference_antitriton.Divide(eff_difference_antitriton,eff_default_antitriton,1,(1-(1.5))/0.88)
scaling_factor = ROOT.TH1D(eff_difference_antitriton)
for iB in range(scaling_factor.GetNbinsX()+1):
    if iB == 0:
        continue
    tmp = scaling_factor.GetBinContent(iB)
    tmp_err = scaling_factor.GetBinError(iB)
    scaling_factor.SetBinContent(iB,1.-tmp*par_fit2g4_antitriton*(1.-par_fit2g4_antitriton))
    scaling_factor.SetBinError(iB,tmp_err*par_fit2g4_antitriton*(1.-par_fit2g4_antitriton))
    if iB == 1:
        scaling_factor.SetBinContent(iB,0)
        scaling_factor.SetBinError(iB,0)
fit_func = ROOT.TF1("fit_func_3","[0]*TMath::Power(x,[1])",1.,10.,2)
#fit_func.SetParLimits(0,0.,1.1)
scaling_factor.Fit(fit_func,"MR","",1.6,3.)
scaling_factor.SetTitle("^{3}#bar{H}")
scaling_factor.GetYaxis().SetTitle("Scaling factor")
scaling_factor.GetYaxis().SetRangeUser(1.,1.2)
scaling_factor.GetXaxis().SetRangeUser(1.6,3)
scaling_factor.Write()
scaling_factor.Draw("")
cAntitriton.Print("ScalingFactorAntitriton.pdf")

out_file.Close()