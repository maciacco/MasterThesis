import ROOT
import numpy as np

ROOT.gStyle.SetStripDecimals(0)
marker = [34, 45, 47, 33, 29]

chi2 = []
chi2_fixmuQ = []
chi2_fixmuQ_nopions = []
fit_res = []
fit_res_fixMuQ = []
fit_res_fixMuQ_nopions = []

fit_res.append([0.732082, 0.375965, 0.157341, 0.12004, -0.783795])
chi2.append(0.853811)
fit_res_fixMuQ.append([1.12657, 0.0999665, -0.0283236])
chi2_fixmuQ.append(11.6751)
fit_res_fixMuQ_nopions.append([1.13196, 0.100004, -0.0284593])
chi2_fixmuQ_nopions.append(0.515368)
fit_res.append([0.263754, 1.26683, 0.158628, 0.12032, -0.785797])
chi2.append(3.71891)
fit_res_fixMuQ.append([1.57095, 0.0997854, -0.0394962])
chi2_fixmuQ.append(120.671)
fit_res_fixMuQ_nopions.append([1.59169, 0.0999803, -0.0400177])
chi2_fixmuQ_nopions.append(1.89386)
fit_res.append([1.19837, -0.518325, 0.155481, 0.119549, -0.780839])
chi2.append(1.08175)
fit_res_fixMuQ.append([0.682777, 0.100033, -0.017166])
chi2_fixmuQ.append(19.358)
fit_res_fixMuQ_nopions.append([0.672583, 0.100017, -0.0169098])
chi2_fixmuQ_nopions.append(1.26066)
fit_res.append([0.590557, 0.741128, 0.158119, 0.120289, -0.785134])
chi2.append(0.716383)
fit_res_fixMuQ.append([1.36059, 0.0999093, -0.0342073])
chi2_fixmuQ.append(41.3203)
fit_res_fixMuQ_nopions.append([1.37137, 0.100009, -0.0344785])
chi2_fixmuQ_nopions.append(0.165534)
fit_res.append([0.873264, 0.0102989, 0.156495, 0.119764, -0.782333])
chi2.append(1.18213)
fit_res_fixMuQ.append([0.892101, 0.0999822, -0.0224288])
chi2_fixmuQ.append(1.19981)
fit_res_fixMuQ_nopions.append([0.892611, 0.100003, -0.0224416])
chi2_fixmuQ_nopions.append(1.13344)

fit_res.append([0.730364, 0.302259, 0.150009, 0.11095, -0.760148])
chi2.append(4.34376)
fit_res_fixMuQ.append([1.04772, 0.0997644, -0.0263413])
chi2_fixmuQ.append(12.6284)
fit_res_fixMuQ_nopions.append([1.05344, 0.0997956, -0.026485])
chi2_fixmuQ_nopions.append(3.96268)
fit_res.append([0.268336, 1.19585, 0.15133, 0.111216, -0.762433])
chi2.append(4.7174)
fit_res_fixMuQ.append([1.50084, 0.0995951, -0.0377335])
chi2_fixmuQ.append(126.55)
fit_res_fixMuQ_nopions.append([1.5235, 0.0997715, -0.0383032])
chi2_fixmuQ_nopions.append(2.67027)
fit_res.append([1.19, -0.594735, 0.148423, 0.110616, -0.757425])
chi2.append(7.53303)
fit_res_fixMuQ.append([0.595799, 0.0998555, -0.0149793])
chi2_fixmuQ.append(35.7424)
fit_res_fixMuQ_nopions.append([0.585925, 0.0989956, -0.014731])
chi2_fixmuQ_nopions.append(5)
fit_res.append([0.589086, 0.66753, 0.150689, 0.111135, -0.761433])
chi2.append(4.78684)
fit_res_fixMuQ.append([1.28152, 0.0997129, -0.0322195])
chi2_fixmuQ.append(43.4794)
fit_res_fixMuQ_nopions.append([1.29367, 0.099786, -0.032525])
chi2_fixmuQ_nopions.append(4.07446)
fit_res.append([0.871353, -0.0635727, 0.149266, 0.110743, -0.758742])
chi2.append(4.07776)
fit_res_fixMuQ.append([0.81406, 0.099781, -0.0204667])
chi2_fixmuQ.append(4.23337)
fit_res_fixMuQ_nopions.append([0.813155, 0.0998048, -0.0204439])
chi2_fixmuQ_nopions.append(4.07014)

fit_res.append([0.759228, -0.238615, 0.138515, 0.115424, -0.857544])
chi2.append(3.40606)
fit_res_fixMuQ.append([0.515186, 0.0732047, -0.0129526])
chi2_fixmuQ.append(7.37409)
fit_res_fixMuQ_nopions.append([0.512595, 0.0732049, -0.0128874])
chi2_fixmuQ_nopions.append(3.13772)
fit_res.append([0.317167, 0.659052, 0.140228, 0.116181, -0.860022])
chi2.append(2.15119)
fit_res_fixMuQ.append([1.00787, 0.0731753, -0.0253395])
chi2_fixmuQ.append(35.8421)
fit_res_fixMuQ_nopions.append([1.01276, 0.0731883, -0.0254624])
chi2_fixmuQ_nopions.append(1.54537)
fit_res.append([1.19678, -1.13738, 0.136317, 0.114356, -0.854184])
chi2.append(9)
fit_res_fixMuQ.append([0.0219876, 0.0733615, -0.000552801])
chi2_fixmuQ.append(108.149)
fit_res_fixMuQ_nopions.append([0.0142132, 0.0728957, -0.00035734])
chi2_fixmuQ_nopions.append(5)
fit_res.append([0.759228, -0.238615, 0.138515, 0.115424, -0.857544])
chi2.append(3.40606)
fit_res_fixMuQ.append([0.515186, 0.0732047, -0.0129526])
chi2_fixmuQ.append(7.37409)
fit_res_fixMuQ_nopions.append([0.512595, 0.0732049, -0.0128874])
chi2_fixmuQ_nopions.append(3.13772)
fit_res.append([0.759228, -0.238615, 0.138515, 0.115424, -0.857544])
chi2.append(3.40606)
fit_res_fixMuQ.append([0.515186, 0.0732047, -0.0129526])
chi2_fixmuQ.append(7.37409)
fit_res_fixMuQ_nopions.append([0.512595, 0.0732049, -0.0128874])
chi2_fixmuQ_nopions.append(3.13772)

fit_res.append([0.825807, 0.136206, 0.0942823, 0.0723968, -0.795644])
chi2.append(2.93335)
fit_res_fixMuQ.append([0.976131, 0.0585533, -0.0245414])
chi2_fixmuQ.append(7.30645)
fit_res_fixMuQ_nopions.append([0.978234, 0.0585657, -0.0245943])
chi2_fixmuQ_nopions.append(2.98952)
fit_res.append([0.329507, 1.0359, 0.0952184, 0.0726401, -0.797982])
chi2.append(6.34478)
fit_res_fixMuQ.append([1.41389, 0.058485, -0.0355474])
chi2_fixmuQ.append(220.935)
fit_res_fixMuQ_nopions.append([1.42864, 0.0584718, -0.0359184])
chi2_fixmuQ_nopions.append(5)
fit_res.append([1.31811, -0.766036, 0.0931977, 0.0721053, -0.792898])
chi2.append(1.946)
fit_res_fixMuQ.append([0.541527, 0.0586082, -0.0136148])
chi2_fixmuQ.append(113.489)
fit_res_fixMuQ_nopions.append([0.527268, 0.0585711, -0.0132563])
chi2_fixmuQ_nopions.append(1.09386)
fit_res.append([0.730795, 0.380368, 0.0946166, 0.0725128, -0.796577])
chi2.append(2.56288)
fit_res_fixMuQ.append([1.13538, 0.0585432, -0.0285452])
chi2_fixmuQ.append(32.922)
fit_res_fixMuQ_nopions.append([1.1406, 0.0585632, -0.0286765])
chi2_fixmuQ_nopions.append(2.79656)
fit_res.append([0.920577, -0.108123, 0.0939439, 0.07228, -0.794696])
chi2.append(3.37534)
fit_res_fixMuQ.append([0.817574, 0.0585581, -0.0205551])
chi2_fixmuQ.append(5.09777)
fit_res_fixMuQ_nopions.append([0.815846, 0.0585684, -0.0205116])
chi2_fixmuQ_nopions.append(3.23815)

fit_res.append([0.714549, -0.271576, 0.188885, 0.160705, -0.878581])
chi2.append(1.62153)
fit_res_fixMuQ.append([0.433449, 0.0926987, -0.0108976])
chi2_fixmuQ.append(4.33548)
fit_res_fixMuQ_nopions.append([0.430914, 0.0926928, -0.0108338])
chi2_fixmuQ_nopions.append(1.34692)
fit_res.append([0.228141, 0.632279, 0.191322, 0.161876, -0.880854])
chi2.append(0.281467)
fit_res_fixMuQ.append([0.892337, 0.0926727, -0.0224347])
chi2_fixmuQ.append(16.1368)
fit_res_fixMuQ_nopions.append([0.895879, 0.0926761, -0.0225238])
chi2_fixmuQ_nopions.append(0.416871)
fit_res.append([1.19689, -1.17718, 0.187196, 0.160096, -0.87708])
chi2.append(4.01259)
fit_res_fixMuQ.append([-0.0262955, 0.0928863, 0.000661107])
chi2_fixmuQ.append(57.7456)
fit_res_fixMuQ_nopions.append([-0.034709, 0.0927288, 0.000872637])
chi2_fixmuQ_nopions.append(2.9769)
fit_res.append([0.714549, -0.271576, 0.188885, 0.160705, -0.878581])
chi2.append(1.62153)
fit_res_fixMuQ.append([0.433449, 0.0926987, -0.0108976])
chi2_fixmuQ.append(4.33548)
fit_res_fixMuQ_nopions.append([0.430914, 0.0926928, -0.0108338])
chi2_fixmuQ_nopions.append(1.34692)
fit_res.append([0.714549, -0.271576, 0.188885, 0.160705, -0.878581])
chi2.append(1.62153)
fit_res_fixMuQ.append([0.433449, 0.0926987, -0.0108976])
chi2_fixmuQ.append(4.33548)
fit_res_fixMuQ_nopions.append([0.430914, 0.0926928, -0.0108338])
chi2_fixmuQ_nopions.append(1.34692)


CENTRALITY_CLASSES = [[0,5], [5,10], [30,50], [10, 30], [50,90]]
centrality_colors = [ROOT.kRed, ROOT.kOrange-3, ROOT.kAzure+4,ROOT.kGreen+2,ROOT.kMagenta+2]
centrality_colors_area = [ROOT.kRed, ROOT.kOrange-3, ROOT.kAzure+4,ROOT.kGreen+2,ROOT.kMagenta+2]

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)

ellipses = []
ellipses_corr = []
ellipses_corr_ = []
ellipses_polarity = []
points = []

cov_files = [open(f"CovMuBMuQ_{i}.yaml", "w") for i in range(5)]
for cov_file in cov_files:
    cov_file.write("independent_variables:\n")

    cov_file.write("- header: {name: Chemical potentials (MeV)}\n")
    cov_file.write("  values:\n")
    for _ in range(2):
        cov_file.write("  - {value: \'$\mu_Q$\'}\n")
        cov_file.write("  - {value: \'$\mu_B$\'}\n")

    cov_file.write("- header: {name: Chemical potentials (MeV)}\n")
    cov_file.write("  values:\n")
    for _ in range(2):
        cov_file.write("  - {value: \'$\mu_Q$\'}\n")
    for _ in range(2):
        cov_file.write("  - {value: \'$\mu_B$\'}\n")

    cov_file.write("dependent_variables:\n")
    cov_file.write("- header: {name: Covariance matrix}\n")
    cov_file.write("  values:\n")

with open("MuBMuQ.yaml", "w") as ofile:
    ofile.write("independent_variables:\n")
    ofile.write("- header: {name: Centrality}\n")
    ofile.write("  values:\n")
    ofile.write("  - {value: \"0-5%\"}\n")
    ofile.write("  - {value: \"5-10%\"}\n")
    ofile.write("  - {value: \"10-30%\"}\n")
    ofile.write("  - {value: \"30-50%\"}\n")
    ofile.write("  - {value: \"50-90%\"}\n")

    ofile.write("dependent_variables:\n")
    ofile.write("- header: {name: '$\mu_Q$ (MeV)'}\n")
    ofile.write("  values:\n")
    for ii_cent, cent in enumerate(CENTRALITY_CLASSES):
        i_cent = ii_cent
        if (ii_cent == 2):
            i_cent = 3
        elif (ii_cent == 3):
            i_cent = 2
        muI = fit_res[int(5*i_cent)][1]
        muI_err_corr = np.abs(fit_res[int(5*i_cent+2)][1]-fit_res[int(5*i_cent+1)][1])*0.5
        muI_err_polarity = np.abs(fit_res[int(5*i_cent+4)][1]-fit_res[int(5*i_cent+3)][1])*0.5
        muI_err = fit_res[int(5*i_cent)][3]
        muB_err = fit_res[int(5*i_cent)][2]
        corr_ = fit_res[int(5*i_cent)][4]

        ofile.write("  - errors:\n")
        ofile.write("    - {label: Uncorrelated uncertainty, symerror: " + "{:.6f}".format(muI_err) + "}\n")
        ofile.write("    - {label: Correlated uncertainty, symerror: " + "{:.6f}".format(np.sqrt(muI_err_corr**2+muI_err_polarity**2)) + "}\n")
        ofile.write("    value: " + "{:.6f}".format(muI) + "\n")

        cov_files[ii_cent].write("  - {value: " + "{:.6f}".format(muI_err ** 2) + "}\n")
        cov_files[ii_cent].write("  - {value: " + "{:.6f}".format(corr_ * muI_err * muB_err) + "}\n")
        cov_files[ii_cent].write("  - {value: " + "{:.6f}".format(corr_ * muI_err * muB_err) + "}\n")
        cov_files[ii_cent].write("  - {value: " + "{:.6f}".format(muB_err ** 2) + "}\n")
        # print(np.sqrt(muI_err_corr**2+muI_err_polarity**2))

    ofile.write("- header: {name: '$\mu_B$ (MeV)'}\n")
    ofile.write("  values:\n")
    for ii_cent, cent in enumerate(CENTRALITY_CLASSES):
        i_cent = ii_cent
        if (ii_cent == 2):
            i_cent = 3
        elif (ii_cent == 3):
            i_cent = 2
        muB = fit_res[int(5*i_cent)][0]
        muB_err_corr = np.abs(fit_res[int(5*i_cent+2)][0]-fit_res[int(5*i_cent+1)][0])*0.5
        muB_err_polarity = np.abs(fit_res[int(5*i_cent+4)][0]-fit_res[int(5*i_cent+3)][0])*0.5
        # print(np.sqrt(muB_err_corr**2+muB_err_polarity**2))
        muB_err = fit_res[int(5*i_cent)][2]

        ofile.write("  - errors:\n")
        ofile.write("    - {label: Uncorrelated uncertainty, symerror: " + "{:.6f}".format(muB_err) + "}\n")
        ofile.write("    - {label: Correlated uncertainty, symerror: " + "{:.6f}".format(np.sqrt(muB_err_corr**2+muB_err_polarity**2)) + "}\n")
        ofile.write("    value: " + "{:.6f}".format(muB) + "\n")
