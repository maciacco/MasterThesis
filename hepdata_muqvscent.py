import ROOT
import numpy as np

ifile = ROOT.TFile("muQvsCent.root")

chi2 = ifile.Get("h;1")
covMat = ifile.Get("covM_0")

with open("Chi2MuQ.yaml", "w") as ofile:
    ofile.write("dependent_variables:\n")
    ofile.write("- header: {name: '$\\chi^{2}$'}\n")
    ofile.write("  values:\n")
    for i in range(0, chi2.GetNbinsX()):
        ofile.write("  - {value: " + "{:.6f}".format(chi2.GetBinContent(i + 1)) + "}\n")
    ofile.write("independent_variables:\n")
    ofile.write("- header: {name: '$\\mu_Q$ (MeV)'}\n")
    ofile.write("  values:\n")
    for i in range(0, chi2.GetNbinsX()):
        ofile.write("  - {value: " + "{:.6f}".format(chi2.GetBinCenter(i + 1)) + "}\n")

with open("CovarianceMatrixMuQCent.yaml", "w") as ofile:
    ofile.write("independent_variables:\n")
    ofile.write("- header: {name: Centrality}\n")
    ofile.write("  values:\n")
    for i in range(0, 5):
        ofile.write("  - {value: \"0-5%\"}\n")
        ofile.write("  - {value: \"5-10%\"}\n")
        ofile.write("  - {value: \"10-30%\"}\n")
        ofile.write("  - {value: \"30-50%\"}\n")
        ofile.write("  - {value: \"50-90%\"}\n")
    ofile.write("- header: {name: Centrality}\n")
    ofile.write("  values:\n")
    for i in range(0, 5):
        ofile.write("  - {value: \"0-5%\"}\n")
    for i in range(0, 5):
        ofile.write("  - {value: \"5-10%\"}\n")
    for i in range(0, 5):
        ofile.write("  - {value: \"10-30%\"}\n")
    for i in range(0, 5):
        ofile.write("  - {value: \"30-50%\"}\n")
    for i in range(0, 5):
        ofile.write("  - {value: \"50-90%\"}\n")
    ofile.write("dependent_variables:\n")
    ofile.write("- header: {name: total covariance matrix}\n")
    ofile.write("  values:\n")
    for i in range(1, 6):
        for j in range(1, 6):
            if j < i + 1:
                ofile.write("  - {value: " + "{:.6f}".format(covMat.GetBinContent(i, 6 - j)) + "}\n")
            else:
                ofile.write("  - {value: " + "{:.6f}".format(covMat.GetBinContent(j, 6 - i)) + "}\n")
