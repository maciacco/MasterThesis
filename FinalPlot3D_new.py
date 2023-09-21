#!/usr/bin/env python3

import ROOT
import numpy as np
ROOT.gStyle.SetHistTopMargin(0)
path_he3 = './He3_PbPb/out'
path_triton = './Triton_PbPb/out'
path_hyp = './Hypertriton_PbPb'
path_proton = './Proton_PbPb/out'
path_pion = './Pion_PbPb/out'
path_omega = './Omega_PbPb'
centrality_classes = [[0, 5], [5, 10], [30, 50], [10, 30], [50,90]]
centrality_colors = [ROOT.kRed, ROOT.kOrange-3,ROOT.kAzure+4,ROOT.kGreen+2, ROOT.kMagenta+2]
particle_ratios = ["#bar{#Omega}^{+} / #Omega^{-}","#pi^{-} / #pi^{+}","#bar{p} / p","{}_{#bar{#Lambda}}^{3}#bar{H} / ^{3}_{#Lambda}H","^{3}#bar{H} / ^{3}H","^{3}#bar{He} / ^{3}He"]
n_part = [2.5,7.5,40,20,70]
n_part_err = [1.,1.,1.,1.,1.]

fixmuQ = False
nopions = False

# chi2 = [1.0596,4.77823,3.71918,3.40207,1.77357]
# chi2_fixmuQ = [10.7356,12.7981,10.4624,7.41971,5.14968]
# chi2_fixmuQ_nopions = [0.695672,4.0161,3.03759,2.90952,1.26367]
chi2 = []
chi2_fixmuQ = []
chi2_fixmuQ_nopions = []

ndf = [4,4,4,4,3]

fit_res = []
fit_res_fixMuQ = []
fit_res_fixMuQ_nopions = []

# fit_res.append([0.646566, 0.459121, 0.160271, 0.123529, -0.792434])
# chi2.append(0.970221)
# fit_res_fixMuQ.append([1.12586, 0.0999647, -0.0283059])
# chi2_fixmuQ.append(15.9348)
# fit_res_fixMuQ_nopions.append([1.13196, 0.100004, -0.0284593])
# chi2_fixmuQ_nopions.append(0.515368)
# fit_res.append([0.179215, 1.3487, 0.161486, 0.123803, -0.794266])
# chi2.append(3.95279)
# fit_res_fixMuQ.append([1.57088, 0.099782, -0.0394945])
# chi2_fixmuQ.append(128.816)
# fit_res_fixMuQ_nopions.append([1.59169, 0.0999803, -0.0400177])
# chi2_fixmuQ_nopions.append(1.89386)
# fit_res.append([1.11276, -0.434212, 0.158374, 0.122998, -0.789521])
# chi2.append(1.10316)
# fit_res_fixMuQ.append([0.680615, 0.100026, -0.0171117])
# chi2_fixmuQ.append(13.1512)
# fit_res_fixMuQ_nopions.append([0.672583, 0.100017, -0.0169098])
# chi2_fixmuQ_nopions.append(1.26066)
# fit_res.append([0.50497, 0.824122, 0.161033, 0.123783, -0.793697])
# chi2.append(0.846032)
# fit_res_fixMuQ.append([1.36011, 0.099905, -0.0341952])
# chi2_fixmuQ.append(47.9625)
# fit_res_fixMuQ_nopions.append([1.37137, 0.100009, -0.0344785])
# chi2_fixmuQ_nopions.append(0.165534)
# fit_res.append([0.787958, 0.0935532, 0.159421, 0.123239, -0.791019])
# chi2.append(1.28637)
# fit_res_fixMuQ.append([0.890967, 0.0999835, -0.0224003])
# chi2_fixmuQ.append(2.00803)
# fit_res_fixMuQ_nopions.append([0.892611, 0.100003, -0.0224416])
# chi2_fixmuQ_nopions.append(1.13344)

# fit_res.append([0.781426, 0.172088, 0.129462, 0.120502, -0.575297])
# chi2.append(0.132713)
# fit_res_fixMuQ.append([1.04704, 0.0997645, -0.0263242])
# chi2_fixmuQ.append(15.8296)
# fit_res_fixMuQ_nopions.append([1.05344, 0.0997956, -0.026485])
# chi2_fixmuQ_nopions.append(3.96268)
# fit_res.append([0.20045, 1.2615, 0.154049, 0.114629, -0.771964])
# chi2.append(4.91523)
# fit_res_fixMuQ.append([1.50106, 0.0995939, -0.0377391])
# chi2_fixmuQ.append(132.243)
# fit_res_fixMuQ_nopions.append([1.5235, 0.0997715, -0.0383032])
# chi2_fixmuQ_nopions.append(2.67027)
# fit_res.append([1.12138, -0.527379, 0.15115, 0.113992, -0.767127])
# chi2.append(7.56368)
# fit_res_fixMuQ.append([0.593652, 0.0998487, -0.0149253])
# chi2_fixmuQ.append(28.3695)
# fit_res_fixMuQ_nopions.append([0.585925, 0.0989956, -0.014731])
# chi2_fixmuQ_nopions.append(5)
# fit_res.append([0.520498, 0.733978, 0.153444, 0.114551, -0.771052])
# chi2.append(4.90904)
# fit_res_fixMuQ.append([1.28119, 0.0997111, -0.0322112])
# chi2_fixmuQ.append(48.6687)
# fit_res_fixMuQ_nopions.append([1.29367, 0.099786, -0.032525])
# chi2_fixmuQ_nopions.append(4.07446)
# fit_res.append([0.80296, 0.00311398, 0.152021, 0.114138, -0.768449])
# chi2.append(4.1706)
# fit_res_fixMuQ.append([0.812729, 0.0997815, -0.0204332])
# chi2_fixmuQ.append(4.11037)
# fit_res_fixMuQ_nopions.append([0.813155, 0.0998048, -0.0204439])
# chi2_fixmuQ_nopions.append(4.07014)

# fit_res.append([0.697399, -0.178802, 0.141653, 0.118879, -0.864176])
# chi2.append(3.39079)
# fit_res_fixMuQ.append([0.51438, 0.0732017, -0.0129323])
# chi2_fixmuQ.append(5.36281)
# fit_res_fixMuQ_nopions.append([0.512595, 0.0732049, -0.0128874])
# chi2_fixmuQ_nopions.append(3.13772)
# fit_res.append([0.256138, 0.717743, 0.143405, 0.119689, -0.866585])
# chi2.append(2.2324)
# fit_res_fixMuQ.append([1.00781, 0.0731757, -0.0253379])
# chi2_fixmuQ.append(39.7421)
# fit_res_fixMuQ_nopions.append([1.01276, 0.0731883, -0.0254624])
# chi2_fixmuQ_nopions.append(1.54537)
# fit_res.append([1.13366, -1.0757, 0.139298, 0.117667, -0.860728])
# chi2.append(9)
# fit_res_fixMuQ.append([0.0208937, 0.0733509, -0.000525299])
# chi2_fixmuQ.append(92.8789)
# fit_res_fixMuQ_nopions.append([0.0142132, 0.0728957, -0.00035734])
# chi2_fixmuQ_nopions.append(5)
# fit_res.append([0.697399, -0.178802, 0.141653, 0.118879, -0.864176])
# chi2.append(3.39079)
# fit_res_fixMuQ.append([0.51438, 0.0732017, -0.0129323])
# chi2_fixmuQ.append(5.36281)
# fit_res_fixMuQ_nopions.append([0.512595, 0.0732049, -0.0128874])
# chi2_fixmuQ_nopions.append(3.13772)
# fit_res.append([0.697399, -0.178802, 0.141653, 0.118879, -0.864176])
# chi2.append(3.39079)
# fit_res_fixMuQ.append([0.51438, 0.0732017, -0.0129323])
# chi2_fixmuQ.append(5.36281)
# fit_res_fixMuQ_nopions.append([0.512595, 0.0732049, -0.0128874])
# chi2_fixmuQ_nopions.append(3.13772)

# fit_res.append([0.800627, 0.160439, 0.0958831, 0.0743065, -0.80318])
# chi2.append(2.93535)
# fit_res_fixMuQ.append([0.975953, 0.0585538, -0.024537])
# chi2_fixmuQ.append(8.5143)
# fit_res_fixMuQ_nopions.append([0.978234, 0.0585657, -0.0245943])
# chi2_fixmuQ_nopions.append(2.98952)
# fit_res.append([0.304551, 1.05982, 0.0968229, 0.074562, -0.805448])
# chi2.append(6.38844)
# fit_res_fixMuQ.append([1.41439, 0.0584859, -0.03556])
# chi2_fixmuQ.append(219.35)
# fit_res_fixMuQ_nopions.append([1.42864, 0.0584718, -0.0359184])
# chi2_fixmuQ_nopions.append(5)
# fit_res.append([1.29288, -0.741484, 0.0947708, 0.0739918, -0.800457])
# chi2.append(1.91356)
# fit_res_fixMuQ.append([0.540445, 0.058606, -0.0135876])
# chi2_fixmuQ.append(101.153)
# fit_res_fixMuQ_nopions.append([0.527268, 0.0585711, -0.0132563])
# chi2_fixmuQ_nopions.append(1.09386)
# fit_res.append([0.705583, 0.404583, 0.0962231, 0.0744291, -0.804094])
# chi2.append(2.57278)
# fit_res_fixMuQ.append([1.13538, 0.0585435, -0.0285453])
# chi2_fixmuQ.append(34.9773)
# fit_res_fixMuQ_nopions.append([1.1406, 0.0585632, -0.0286765])
# chi2_fixmuQ_nopions.append(2.79656)
# fit_res.append([0.895439, -0.0838715, 0.0955375, 0.0741821, -0.802246])
# chi2.append(3.36936)
# fit_res_fixMuQ.append([0.817058, 0.0585581, -0.0205421])
# chi2_fixmuQ.append(4.23893)
# fit_res_fixMuQ_nopions.append([0.815846, 0.0585684, -0.0205116])
# chi2_fixmuQ_nopions.append(3.23815)

# fit_res.append([0.649342, -0.208717, 0.19268, 0.164762, -0.883546])
# chi2.append(1.59046)
# fit_res_fixMuQ.append([0.432747, 0.0926938, -0.0108799])
# chi2_fixmuQ.append(3.04905)
# fit_res_fixMuQ_nopions.append([0.430914, 0.0926928, -0.0108338])
# chi2_fixmuQ_nopions.append(1.34692)
# fit_res.append([0.163009, 0.694674, 0.195181, 0.166011, -0.885776])
# chi2.append(0.272443)
# fit_res_fixMuQ.append([0.892224, 0.0926737, -0.0224319])
# chi2_fixmuQ.append(18.4048)
# fit_res_fixMuQ_nopions.append([0.895879, 0.0926761, -0.0225238])
# chi2_fixmuQ_nopions.append(0.416871)
# fit_res.append([1.13111, -1.11327, 0.190841, 0.16404, -0.881969])
# chi2.append(3.97182)
# fit_res_fixMuQ.append([-0.0271187, 0.0928731, 0.000681804])
# chi2_fixmuQ.append(49.7348)
# fit_res_fixMuQ_nopions.append([-0.034709, 0.0927288, 0.000872637])
# chi2_fixmuQ_nopions.append(2.9769)
# fit_res.append([0.649342, -0.208717, 0.19268, 0.164762, -0.883546])
# chi2.append(1.59046)
# fit_res_fixMuQ.append([0.432747, 0.0926938, -0.0108799])
# chi2_fixmuQ.append(3.04905)
# fit_res_fixMuQ_nopions.append([0.430914, 0.0926928, -0.0108338])
# chi2_fixmuQ_nopions.append(1.34692)
# fit_res.append([0.649342, -0.208717, 0.19268, 0.164762, -0.883546])
# chi2.append(1.59046)
# fit_res_fixMuQ.append([0.432747, 0.0926938, -0.0108799])
# chi2_fixmuQ.append(3.04905)
# fit_res_fixMuQ_nopions.append([0.430914, 0.0926928, -0.0108338])
# chi2_fixmuQ_nopions.append(1.34692)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# fit_res.append([0.611881, 0.492826, 0.161671, 0.125195, -0.796389])
# chi2.append(1.02021)
# fit_res_fixMuQ.append([1.12563, 0.0999642, -0.0283])
# chi2_fixmuQ.append(17.714)
# fit_res_fixMuQ_nopions.append([1.13196, 0.100004, -0.0284593])
# chi2_fixmuQ_nopions.append(0.515368)
# fit_res.append([0.145081, 1.38176, 0.162856, 0.125467, -0.798152])
# chi2.append(4.05147)
# fit_res_fixMuQ.append([1.57095, 0.099781, -0.0394963])
# chi2_fixmuQ.append(131.532)
# fit_res_fixMuQ_nopions.append([1.59169, 0.0999803, -0.0400177])
# chi2_fixmuQ_nopions.append(1.89386)
# fit_res.append([1.07792, -0.400021, 0.159759, 0.124645, -0.793503])
# chi2.append(1.11432)
# fit_res_fixMuQ.append([0.679763, 0.100023, -0.0170903])
# chi2_fixmuQ.append(11.0382)
# fit_res_fixMuQ_nopions.append([0.672583, 0.100017, -0.0169098])
# chi2_fixmuQ_nopions.append(1.26066)
# fit_res.append([0.470279, 0.857752, 0.162426, 0.125451, -0.79762])
# chi2.append(0.901731)
# fit_res_fixMuQ.append([1.35998, 0.0999037, -0.0341922])
# chi2_fixmuQ.append(50.4854)
# fit_res_fixMuQ_nopions.append([1.37137, 0.100009, -0.0344785])
# chi2_fixmuQ_nopions.append(0.165534)
# fit_res.append([0.753342, 0.127303, 0.160819, 0.124898, -0.794996])
# chi2.append(1.33128)
# fit_res_fixMuQ.append([0.890589, 0.0999843, -0.0223908])
# chi2_fixmuQ.append(2.56582)
# fit_res_fixMuQ_nopions.append([0.892611, 0.100003, -0.0224416])
# chi2_fixmuQ_nopions.append(1.13344)

# fit_res.append([0.634577, 0.395302, 0.154117, 0.116017, -0.774337])
# chi2.append(4.49366)
# fit_res_fixMuQ.append([1.04683, 0.0997648, -0.026319])
# chi2_fixmuQ.append(17.1163)
# fit_res_fixMuQ_nopions.append([1.05344, 0.0997956, -0.026485])
# chi2_fixmuQ_nopions.append(3.96268)
# fit_res.append([0.173597, 1.28747, 0.15538, 0.116291, -0.776422])
# chi2.append(4.99626)
# fit_res_fixMuQ.append([1.50128, 0.0995939, -0.0377446])
# chi2_fixmuQ.append(133.749)
# fit_res_fixMuQ_nopions.append([1.5235, 0.0997715, -0.0383032])
# chi2_fixmuQ_nopions.append(2.67027)
# fit_res.append([1.09399, -0.500516, 0.152482, 0.115636, -0.771664])
# chi2.append(7.57676)
# fit_res_fixMuQ.append([0.592808, 0.0998458, -0.0149041])
# chi2_fixmuQ.append(25.7504)
# fit_res_fixMuQ_nopions.append([0.585925, 0.0989956, -0.014731])
# chi2_fixmuQ_nopions.append(5)
# fit_res.append([0.493244, 0.760375, 0.154789, 0.116215, -0.775545])
# chi2.append(4.95974)
# fit_res_fixMuQ.append([1.28115, 0.0997107, -0.0322101])
# chi2_fixmuQ.append(50.4917)
# fit_res_fixMuQ_nopions.append([1.29367, 0.099786, -0.032525])
# chi2_fixmuQ_nopions.append(4.07446)
# fit_res.append([0.775753, 0.0296193, 0.153364, 0.11579, -0.772982])
# chi2.append(4.20893)
# fit_res_fixMuQ.append([0.812271, 0.099782, -0.0204217])
# chi2_fixmuQ.append(4.25541)
# fit_res_fixMuQ_nopions.append([0.813155, 0.0998048, -0.0204439])
# chi2_fixmuQ_nopions.append(4.07014)

# fit_res.append([0.677893, -0.159951, 0.143174, 0.120557, -0.867234])
# chi2.append(3.38539)
# fit_res_fixMuQ.append([0.514132, 0.0732009, -0.012926])
# chi2_fixmuQ.append(4.86164)
# fit_res_fixMuQ_nopions.append([0.512595, 0.0732049, -0.0128874])
# chi2_fixmuQ_nopions.append(3.13772)
# fit_res.append([0.237029, 0.736122, 0.144947, 0.121393, -0.869613])
# chi2.append(2.25918)
# fit_res_fixMuQ.append([1.00785, 0.0731759, -0.025339])
# chi2_fixmuQ.append(40.5747)
# fit_res_fixMuQ_nopions.append([1.01276, 0.0731883, -0.0254624])
# chi2_fixmuQ_nopions.append(1.54537)
# fit_res.append([1.11364, -1.05615, 0.140754, 0.119284, -0.863774])
# chi2.append(9)
# fit_res_fixMuQ.append([0.0204826, 0.0733471, -0.000514964])
# chi2_fixmuQ.append(87.7127)
# fit_res_fixMuQ_nopions.append([0.0142132, 0.0728957, -0.00035734])
# chi2_fixmuQ_nopions.append(5)
# fit_res.append([0.677893, -0.159951, 0.143174, 0.120557, -0.867234])
# chi2.append(3.38539)
# fit_res_fixMuQ.append([0.514132, 0.0732009, -0.012926])
# chi2_fixmuQ.append(4.86164)
# fit_res_fixMuQ_nopions.append([0.512595, 0.0732049, -0.0128874])
# chi2_fixmuQ_nopions.append(3.13772)
# fit_res.append([0.677893, -0.159951, 0.143174, 0.120557, -0.867234])
# chi2.append(3.38539)
# fit_res_fixMuQ.append([0.514132, 0.0732009, -0.012926])
# chi2_fixmuQ.append(4.86164)
# fit_res_fixMuQ_nopions.append([0.512595, 0.0732049, -0.0128874])
# chi2_fixmuQ_nopions.append(3.13772)

# fit_res.append([0.80553, 0.15569, 0.0967162, 0.0753125, -0.80698])
# chi2.append(2.9354)
# fit_res_fixMuQ.append([0.976071, 0.0585541, -0.0245399])
# chi2_fixmuQ.append(8.07804)
# fit_res_fixMuQ_nopions.append([0.978234, 0.0585657, -0.0245943])
# chi2_fixmuQ_nopions.append(2.98952)
# fit_res.append([0.309528, 1.05507, 0.0976674, 0.0755773, -0.809232])
# chi2.append(6.38841)
# fit_res_fixMuQ.append([1.41488, 0.0584872, -0.0355722])
# chi2_fixmuQ.append(211.865)
# fit_res_fixMuQ_nopions.append([1.42864, 0.0584718, -0.0359184])
# chi2_fixmuQ_nopions.append(5)
# fit_res.append([1.29772, -0.746214, 0.095595, 0.0749882, -0.804287])
# chi2.append(1.91157)
# fit_res_fixMuQ.append([0.540206, 0.058606, -0.0135816])
# chi2_fixmuQ.append(99.7734)
# fit_res_fixMuQ_nopions.append([0.527268, 0.0585711, -0.0132563])
# chi2_fixmuQ_nopions.append(1.09386)
# fit_res.append([0.710507, 0.39983, 0.0970604, 0.0754391, -0.807888])
# chi2.append(2.57464)
# fit_res_fixMuQ.append([1.13559, 0.058544, -0.0285505])
# chi2_fixmuQ.append(33.4162)
# fit_res_fixMuQ_nopions.append([1.1406, 0.0585632, -0.0286765])
# chi2_fixmuQ_nopions.append(2.79656)
# fit_res.append([0.900319, -0.0886135, 0.0963667, 0.0751842, -0.806055])
# chi2.append(3.36734)
# fit_res_fixMuQ.append([0.817107, 0.0585585, -0.0205433])
# chi2_fixmuQ.append(4.34208)
# fit_res_fixMuQ_nopions.append([0.815846, 0.0585684, -0.0205116])
# chi2_fixmuQ_nopions.append(3.23815)

# fit_res.append([0.656659, -0.215802, 0.194585, 0.166852, -0.885967])
# chi2.append(1.59099)
# fit_res_fixMuQ.append([0.432767, 0.0926943, -0.0108804])
# chi2_fixmuQ.append(3.12024)
# fit_res_fixMuQ_nopions.append([0.430914, 0.0926928, -0.0108338])
# chi2_fixmuQ_nopions.append(1.34692)
# fit_res.append([0.170469, 0.687535, 0.197126, 0.168137, -0.888176])
# chi2.append(0.273679)
# fit_res_fixMuQ.append([0.892376, 0.0926738, -0.0224357])
# chi2_fixmuQ.append(17.5955)
# fit_res_fixMuQ_nopions.append([0.895879, 0.0926761, -0.0225238])
# chi2_fixmuQ_nopions.append(0.416871)
# fit_res.append([1.1383, -1.12026, 0.192732, 0.166113, -0.884417])
# chi2.append(3.97179)
# fit_res_fixMuQ.append([-0.0272752, 0.0928723, 0.000685739])
# chi2_fixmuQ.append(49.1472)
# fit_res_fixMuQ_nopions.append([-0.034709, 0.0927288, 0.000872637])
# chi2_fixmuQ_nopions.append(2.9769)
# fit_res.append([0.656659, -0.215802, 0.194585, 0.166852, -0.885967])
# chi2.append(1.59099)
# fit_res_fixMuQ.append([0.432767, 0.0926943, -0.0108804])
# chi2_fixmuQ.append(3.12024)
# fit_res_fixMuQ_nopions.append([0.430914, 0.0926928, -0.0108338])
# chi2_fixmuQ_nopions.append(1.34692)
# fit_res.append([0.656659, -0.215802, 0.194585, 0.166852, -0.885967])
# chi2.append(1.59099)
# fit_res_fixMuQ.append([0.432767, 0.0926943, -0.0108804])
# chi2_fixmuQ.append(3.12024)
# fit_res_fixMuQ_nopions.append([0.430914, 0.0926928, -0.0108338])
# chi2_fixmuQ_nopions.append(1.34692)



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
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


fit_val=[[0.996927,0.985795,0.962759,0.96744,0.967332,0.99511],
[0.997528,0.986755,0.964582,0.968351,0.968615,0.995346],
[1.0019,0.993512,0.97713,0.974125,0.977209,0.996642],
[0.998884,0.98768,0.965026,0.966724,0.968067,0.995239],
[1.00222,0.994464,0.979506,0.97608,-1,0.996947]]

fit_val_fixmuQ=[[1.00022,0.986473,0.959204,0.958864,0.961662,0.994283],
[1.00021,0.987157,0.961242,0.960918,0.963579,0.994573],
[1.00011,0.993497,0.980246,0.980079,0.981448,0.997257],
[1.0002,0.987818,0.963211,0.962903,0.965431,0.994854],
[1.00009,0.994564,0.983469,0.983329,-1,0.997708]]

fit_val_fixmuQ_nopions=[[1.00023,0.98604,0.957918,0.957567,0.960451,0.994099],
[1.00021,0.986979,0.960711,0.960383,0.96308,0.994498],
[1.0001,0.993672,0.980775,0.980612,0.981945,0.997331],
[1.0002,0.98786,0.963337,0.96303,0.96555,0.994872],
[1.00009,0.994633,0.983678,0.98354,-1,0.997737]]

def muBmuQ(fixmuQ,nopions,i_cent):
    muQ = 0
    muB = 0
    muB_err = 0
    muQ_err = 0
    muQ_corr = 0
    muB_corr = 0
    muB_polarity = 0
    muQ_polarity = 0
    fit_r = fit_res
    if fixmuQ:
        fit_r = fit_res_fixMuQ
        if nopions:
            fit_r = fit_res_fixMuQ_nopions
    if not fixmuQ and not nopions:
        muB = fit_r[i_cent*5][0]
        muB_err = fit_r[i_cent*5][2]
        muB_corr = np.abs(fit_r[(i_cent)*5+1][0]-fit_r[(i_cent)*5+2][0])*0.5
        muB_polarity = np.abs(fit_r[(i_cent)*5+3][0]-fit_r[(i_cent)*5+4][0])*0.5
        muQ = fit_r[i_cent*5][1]
        muQ_err = fit_r[i_cent*5][3]
        muQ_corr = np.abs((fit_r[(i_cent)*5+1][1]-fit_r[(i_cent)*5+2][1])*0.5)
        muQ_polarity = np.abs((fit_r[(i_cent)*5+3][1]-fit_r[(i_cent)*5+4][1])*0.5)
    else:
        muB = fit_r[i_cent*5][0]
        muB_err = fit_r[i_cent*5][1]
        muB_corr = np.abs(fit_r[(i_cent)*5+1][0]-fit_r[(i_cent)*5+2][0])*0.5
        muB_polarity = np.abs(fit_r[(i_cent)*5+3][0]-fit_r[(i_cent)*5+4][0])*0.5
        muQ = fit_r[i_cent*5][2]
        muQ_err = fit_r[i_cent*5][1]/fit_r[i_cent*5][0]*muQ
        muQ_corr = 0
        muQ_polarity = 0
    return [muB,muB_err,muB_corr,muQ,muQ_err,muQ_corr,muB_polarity,muQ_polarity]

apply_scaling_radius = False

TLATEX_TEXT_SIZE = 28
CENT_LIMIT_NUCLEI = 4

ROOT.gStyle.SetPadTickY(1)
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)
ROOT.gStyle.SetTextFont(44)

file_out = ROOT.TFile.Open('FinalPlot3D_new_2.root', 'recreate') # fixmuQ_nopions
if apply_scaling_radius:
    file_out = ROOT.TFile.Open('FinalPlot3D_LHC22b9.root', 'recreate')
file_hijing = ROOT.TFile.Open('HIJINGRatios.root')

hMuBCentFrame = ROOT.TH1D("frame",";Centrality (%);#it{#mu}_{#it{B}} (MeV)",1,0,90)
hMuQCentFrame = ROOT.TH1D("frame",";#LT#it{N}_{part}#GT;#it{#mu}_{#it{Q}} (MeV)",1,0,450)
gMuBCent = ROOT.TGraphErrors()
gMuBCentCorrUnc = ROOT.TGraphErrors()
gMuBCentPolarityUnc = ROOT.TGraphErrors()
gMuQCent = ROOT.TGraphErrors()
gMuQCentCorrUnc = ROOT.TGraphErrors()
gMuQCentPolarityUnc = ROOT.TGraphErrors()

for i_cent, cent in enumerate(centrality_classes):
    fixmuQ_nopions = ''
    if fixmuQ:
        fixmuQ_nopions = 'fixmuQ'
    if fixmuQ and nopions:
        fixmuQ_nopions = 'fixmuQ_nopions'

    # get histograms
    file_material_budget = ROOT.TFile.Open("MaterialBudgetUncertainty.root")
    file_material_budget_hyp = ROOT.TFile.Open("MaterialBudgetUncertaintyHyp.root")
    file_he3 = ROOT.TFile.Open(path_he3 + '/SpectraHe3_kINT7.root')
    file_triton = ROOT.TFile.Open(path_triton + '/SpectraHe3.root')
    file_hyp = ROOT.TFile.Open(path_hyp + '/Ratio.root')
    file_proton = ROOT.TFile.Open(path_proton + '/out_tof_extend_4_old/SystematicsAllEPtNotCombinedTOF_extend_4.root')
    file_pion = ROOT.TFile.Open(path_pion + '/SystematicsAllEPtNotCombined_extend2.root')
    file_omega = ROOT.TFile.Open(path_omega + '/ratio_cutCompetingMass-3.root')
    file_he3_syst = ROOT.TFile.Open(path_he3 + '/SystematicsAll_extend2.root')
    file_triton_syst = ROOT.TFile.Open(path_triton + '/SystematicsAll_extend.root')
    file_he3_syst_eff_prim = ROOT.TFile.Open(path_he3 + '/SystematicsEfficiencyPrimary.root')
    file_hyp_syst = ROOT.TFile.Open(path_hyp + '/SystematicsRatio.root')
    # file_he3_scaling = ROOT.TFile.Open('He3_PbPb/MaterialRadiusCheck.root')
    # file_proton_scaling = ROOT.TFile.Open('Proton_PbPb/MaterialRadiusCheck.root')
    # file_pion_scaling = ROOT.TFile.Open('Pion_PbPb/MaterialRadiusCheck.root')
    file_out.cd()
    ratio_he3 = file_he3.Get(f'1.0_89_0.1_2.5_1_1_1/fRatio_{cent[0]}_{cent[1]}')
    ratio_proton = file_proton.Get(f'fRatio_{cent[0]}_{cent[1]}')
    ratio_pion = file_pion.Get(f'fRatio_{cent[0]}_{cent[1]}')
    ratio_omega = file_omega.Get(f'h_ratio_{cent[0]}_{cent[1]}')
    ratio_triton = ROOT.TH1D()
    ratio_hyp = ROOT.TH1D()
    ratio_triton_distribution = ROOT.TH1D()
    ratio_hyp_distribution = ROOT.TH1D()
    if i_cent < CENT_LIMIT_NUCLEI:
        ratio_hyp = file_hyp.Get(f'fRatio_{cent[0]}_{cent[1]}')
    ratio_triton = file_triton.Get(f'1.0_69_0.1_2.5_1_1_1/fRatio_{cent[0]}_{cent[1]}')
    ratio_triton_distribution = file_triton_syst.Get(f'hist/fFitPar_{cent[0]}_{cent[1]}')
    ratio_hyp_distribution = file_hyp_syst.Get(f'fParameterDistribution_{cent[0]}_{cent[1]}')
    ratio_he3_distribution = file_he3_syst.Get(f'hist/fFitPar_{cent[0]}_{cent[1]}')
    ratio_he3_distribution_eff_prim = file_he3_syst_eff_prim.Get(f'fRatioDistribution_{cent[0]}_{cent[1]}')
    ratio_proton_pt_correlated = file_proton.Get(f'fRatioDistributionTrials_{cent[0]}_{cent[1]}')
    ratio_pion_pt_correlated = file_pion.Get(f'fRatioDistributionTrials_{cent[0]}_{cent[1]}')
    ratio_omega_distribution = file_omega.Get(f'h_sys;{1+2*i_cent}')

    ratio_he3.Fit("pol0","QR","",2.,10.)
    ratio_he3.Write()

    # triton
    ratio_triton.Fit("pol0","QR","",1.6,3.)
    ratio_triton.Write()

    if i_cent < CENT_LIMIT_NUCLEI:
        ratio_hyp.Fit("pol0","QR","",0.,35.)
        ratio_hyp.Write()

    # proton
    ratio_proton.Fit("pol0","QR","",0.5,3.)
    ratio_proton.Write()

    # pion
    ratio_pion.Fit("pol0","QR","",0.7,1.6)
    ratio_pion.Write()

    # omega
    ratio_omega.Fit("pol0","QR","",1,10)
    ratio_omega.Write()

    # get fit functions
    fit_he3 = ratio_he3.GetFunction("pol0")
    fit_proton = ratio_proton.GetFunction("pol0")
    fit_pion = ratio_pion.GetFunction("pol0")
    fit_omega = ratio_omega.GetFunction("pol0")
    fit_triton = ratio_triton.GetFunction("pol0")
    fit_hyp = ROOT.TF1()
    if i_cent < CENT_LIMIT_NUCLEI:
        fit_hyp = ratio_hyp.GetFunction("pol0")

    # get ratios and errors
    ratio_he3 = fit_he3.GetParameter(0)
    ratio_he3_err = fit_he3.GetParError(0)/ratio_he3

    ratio_triton_err = 0.
    ratio_hyp_err = 0.
    ratio_triton = fit_triton.GetParameter(0)
    ratio_triton_err = fit_triton.GetParError(0)/ratio_triton
    if i_cent < CENT_LIMIT_NUCLEI:
        ratio_hyp = fit_hyp.GetParameter(0)
        ratio_hyp_err = fit_hyp.GetParError(0)/ratio_hyp
    
    ratio_proton = fit_proton.GetParameter(0)
    ratio_proton_err = 0#fit_proton.GetParError(0)

    ratio_pion = fit_pion.GetParameter(0)
    ratio_pion_err = 0#fit_pion.GetParError(0)

    ratio_omega = fit_omega.GetParameter(0)
    ratio_omega_err = fit_omega.GetParError(0)

    # systematic error
    syst_he3 = ratio_he3_distribution.GetRMS()/ratio_he3
    syst_he3_eff_prim = ratio_he3_distribution_eff_prim.GetRMS()/ratio_he3
    syst_omega = ratio_omega_distribution.GetRMS()
    syst_triton = 0.
    syst_hyp = 0.
    syst_triton = ratio_triton_distribution.GetRMS()/ratio_triton
    if i_cent < CENT_LIMIT_NUCLEI:
        syst_hyp = ratio_hyp_distribution.GetRMS()/ratio_hyp

    # uncorrelated sys
    syst_he3 = syst_he3*ratio_he3
    syst_he3_eff_prim = syst_he3_eff_prim*ratio_he3
    syst_he3 = np.sqrt(syst_he3*syst_he3+syst_he3_eff_prim*syst_he3_eff_prim)
    syst_triton = syst_triton*ratio_triton
    syst_triton = np.sqrt(syst_triton*syst_triton)
    if i_cent < CENT_LIMIT_NUCLEI:
        syst_hyp = syst_hyp*ratio_hyp
        syst_hyp = np.sqrt(syst_hyp*syst_hyp)
    syst_proton = fit_proton.GetParError(0)
    syst_pion = fit_pion.GetParError(0)
    
    # final plot
    stat_proton = ratio_proton_err
    stat_pion = ratio_pion_err
    stat_he3 = ratio_he3_err
    stat_hyp = 0
    stat_triton = ratio_triton_err
    if i_cent < CENT_LIMIT_NUCLEI:
        stat_hyp = ratio_hyp_err
    ratios_vs_b = ROOT.TH2D(f'fRatio_vs_b_{cent[0]}_{cent[1]}', ';#it{B}+#it{S}/3;#it{Q}-#it{S}/3;Antimatter / Matter', 10, -0.5, 9.5, 7, -0.5, 6.5)
    ratios_vs_b.SetBinContent(10, 7, ratio_he3)
    ratios_vs_b.SetBinError(10, 7, np.sqrt(syst_he3*syst_he3+stat_he3*stat_he3))
    ratios_vs_b.SetBinContent(4, 4, ratio_proton)
    ratios_vs_b.SetBinError(4, 4, np.sqrt(syst_proton*syst_proton+stat_proton*stat_proton))
    ratios_vs_b.SetBinContent(1, 4, ratio_pion)
    ratios_vs_b.SetBinError(1, 4, np.sqrt(syst_pion*syst_pion+stat_pion*stat_pion))
    ratios_vs_b.SetBinContent(1, 1, ratio_omega)
    ratios_vs_b.SetBinError(1, 1, np.sqrt(syst_omega*syst_omega+ratio_omega_err*ratio_omega_err))
    ratios_vs_b.SetBinContent(10, 4, ratio_triton)
    ratios_vs_b.SetBinError(10, 4, np.sqrt(syst_triton*syst_triton+stat_triton*stat_triton))
    if i_cent < CENT_LIMIT_NUCLEI:
        ratios_vs_b.SetBinContent(9, 5, ratio_hyp)
        ratios_vs_b.SetBinError(9, 5, np.sqrt(syst_hyp*syst_hyp+stat_hyp*stat_hyp))

    # set labels
    a = ratios_vs_b.GetXaxis()
    a.SetBit(ROOT.TAxis.kLabelsHori)
    a.SetBinLabel(1, '0')
    a.SetBinLabel(4, '1')
    a.SetBinLabel(7, '2')
    a.SetBinLabel(10, '3')
    a.SetLabelSize(0.062)

    aa = ratios_vs_b.GetYaxis()
    aa.SetBit(ROOT.TAxis.kLabelsHori)
    aa.SetBinLabel(1, '0')
    aa.SetBinLabel(4, '1')
    aa.SetBinLabel(7, '2')
    aa.SetLabelSize(0.062)

    # chi2 text
    formatted_chi2 = "{:.2f}".format(chi2[i_cent*5])
    if fixmuQ:
        formatted_chi2 = "{:.2f}".format(chi2_fixmuQ[i_cent*5])
        if nopions:
            formatted_chi2 = "{:.2f}".format(chi2_fixmuQ_nopions[i_cent*5])
    dof = ndf[i_cent]
    if fixmuQ and not nopions:
        dof = ndf[i_cent] + 1
    text_chi2 = ROOT.TLatex(-0.27, -0.44, "#chi^{2}/NDF = "+formatted_chi2+"/"+str(dof))
    text_chi2.SetTextSize(TLATEX_TEXT_SIZE)
    text_chi2.SetTextColor(ROOT.kBlack)

    # write to file
    ratios_vs_b.Write()

    # make final plot for approval (preliminary)
    leg_ratios_particle = ROOT.TLegend(0.164659,0.291729,0.359438,0.44812)
    leg_ratios_particle.SetTextFont(43)
    leg_ratios_particle.SetTextSize(22)
    cRatiosParticle = ROOT.TCanvas(f"cRatiosParticle_{cent[0]}_{cent[1]}",f"cRatiosParticle_{cent[0]}_{cent[1]}",500,500)
    pad1 = ROOT.TPad("pad1","pad1",0.0,0.3,1.0,1.0,0)
    pad1.SetFillColor(0)
    pad1.SetFrameBorderMode(0)
    pad1.SetBottomMargin(0.0)
    pad1.Draw()
    pad2 = ROOT.TPad("pad2", "pad2", 0.0, 0.0, 1.0, 0.3, 0)
    pad2.SetFillColor(0)
    pad2.SetFrameBorderMode(0)
    pad2.SetFillColor(0)
    pad2.SetFrameBorderMode(0)
    pad2.SetTopMargin(0.)
    pad2.SetBottomMargin(0.27)
    pad2.SetGridy()
    pad2.Draw()
    pad1.cd()
    hRatiosParticle = ROOT.TH1D(f"hRatiosParticle_{cent[0]}_{cent[1]}",f" ",6,0,6)
    hNSigmaRatioFitParticle = ROOT.TH1D(f"hNSigmaRatioFitParticle_{cent[0]}_{cent[1]}",f"{cent[0]}-{cent[1]}%",6,0,6)
    hRatiosParticleFit = ROOT.TH1D(f"hRatiosParticleFit_{cent[0]}_{cent[1]}",f"{cent[0]}-{cent[1]}%",6,0,6)
    hRatiosParticle.GetYaxis().SetTitle("Ratio")
    hRatiosParticle.GetYaxis().CenterTitle()
    hRatiosParticle.GetXaxis().SetTickSize(0.07*0.4)
    hRatiosParticle.GetYaxis().SetTitleSize(0.13*0.45)
    hRatiosParticle.GetYaxis().SetLabelSize(0.12*0.42)
    hRatiosParticle.GetYaxis().SetTitleOffset(0.9)
    for i_part in range(0,6):
        hRatiosParticle.GetXaxis().SetBinLabel(i_part+1,particle_ratios[i_part])
        hRatiosParticle.GetXaxis().SetLabelSize(0.06)
        ratio = -9999
        ratio_err = 0
        fit = -999
        if i_part == 0:
            ratio = ratio_omega
            ratio_err = np.sqrt(syst_omega*syst_omega+ratio_omega_err*ratio_omega_err)
            fit = fit_val[i_cent][5]
            if fixmuQ:
                fit = fit_val_fixmuQ[i_cent][5]
            if fixmuQ and nopions:
                fit = fit_val_fixmuQ_nopions[i_cent][5]
        if i_part == 1:
            ratio = ratio_pion
            ratio_err = np.sqrt(syst_pion*syst_pion+stat_pion*stat_pion)
            fit = fit_val[i_cent][0]
            if fixmuQ:
                fit = fit_val_fixmuQ[i_cent][0]
            if fixmuQ and nopions:
                fit = fit_val_fixmuQ_nopions[i_cent][0]
        if i_part == 2:
            ratio = ratio_proton
            ratio_err = np.sqrt(syst_proton*syst_proton+stat_proton*stat_proton)
            fit = fit_val[i_cent][1]
            if fixmuQ:
                fit = fit_val_fixmuQ[i_cent][1]
            if fixmuQ and nopions:
                fit = fit_val_fixmuQ_nopions[i_cent][1]
        if i_part == 3 and i_cent < CENT_LIMIT_NUCLEI:
            ratio = ratio_hyp
            ratio_err = np.sqrt(syst_hyp*syst_hyp+stat_hyp*stat_hyp)
            fit = fit_val[i_cent][4]
            if fixmuQ:
                fit = fit_val_fixmuQ[i_cent][4]
            if fixmuQ and nopions:
                fit = fit_val_fixmuQ_nopions[i_cent][4]
        if i_part == 4:
            ratio = ratio_triton
            ratio_err = np.sqrt(syst_triton*syst_triton+stat_triton*stat_triton)
            fit = fit_val[i_cent][3]
            if fixmuQ:
                fit = fit_val_fixmuQ[i_cent][3]
            if fixmuQ and nopions:
                fit = fit_val_fixmuQ_nopions[i_cent][3]
        if i_part == 5:
            ratio = ratio_he3
            ratio_err = np.sqrt(syst_he3*syst_he3+stat_he3*stat_he3)
            fit = fit_val[i_cent][2]
            if fixmuQ:
                fit = fit_val_fixmuQ[i_cent][2]
            if fixmuQ and nopions:
                fit = fit_val_fixmuQ_nopions[i_cent][2]
        hRatiosParticle.SetBinContent(i_part+1,ratio)
        hRatiosParticle.SetBinError(i_part+1,ratio_err)
        hRatiosParticleFit.SetBinContent(i_part+1,fit)
        hRatiosParticleFit.SetBinError(i_part+1,0)
        print(f"fit = {fit}")
    hRatiosParticle.GetYaxis().SetNdivisions(10)
    hRatiosParticle.GetYaxis().SetRangeUser(0.67,1.18)
    gRatiosParticle = ROOT.TGraphErrors(hRatiosParticle)
    gRatiosParticle.SetName(f"gRatiosParticle_{cent[0]}_{cent[1]}")
    gRatiosParticleFit = ROOT.TGraphErrors(hRatiosParticleFit)
    gRatiosParticleFit.SetName(f"gRatiosParticleFit_{cent[0]}_{cent[1]}")
    for i_part in range(0,6):
        gRatiosParticle.SetPointError(i_part,0,hRatiosParticle.GetBinError(i_part+1))
        gRatiosParticleFit.SetPointError(i_part,0.3,0)
    
    # data to fit ratio
    for i_ticks in range(6):
        hNSigmaRatioFitParticle.GetYaxis().SetNdivisions(5)
    for i_part in range(0,6):
        hNSigmaRatioFitParticle.GetXaxis().SetBinLabel(i_part+1,particle_ratios[i_part])
        hNSigmaRatioFitParticle.GetXaxis().SetLabelSize(0.2)
        ratio = -999.
        ratio_err = -1.
        fit = -9999.
        if i_part == 0:
            ratio = ratio_omega
            ratio_err = np.sqrt(syst_omega*syst_omega+ratio_omega_err*ratio_omega_err)
            fit = fit_val[i_cent][5]
            if fixmuQ:
                fit = fit_val_fixmuQ[i_cent][5]
            if fixmuQ and nopions:
                fit = fit_val_fixmuQ_nopions[i_cent][5]
        if i_part == 1:
            ratio = ratio_pion
            ratio_err = np.sqrt(syst_pion*syst_pion+stat_pion*stat_pion)
            fit = fit_val[i_cent][0]
            if fixmuQ:
                fit = fit_val_fixmuQ[i_cent][0]
            if fixmuQ and nopions:
                fit = fit_val_fixmuQ_nopions[i_cent][0]
        if i_part == 2:
            ratio = ratio_proton
            ratio_err = np.sqrt(syst_proton*syst_proton+stat_proton*stat_proton)
            fit = fit_val[i_cent][1]
            if fixmuQ:
                fit = fit_val_fixmuQ[i_cent][1]
            if fixmuQ and nopions:
                fit = fit_val_fixmuQ_nopions[i_cent][1]
        if i_part == 3 and i_cent < CENT_LIMIT_NUCLEI:
            ratio = ratio_hyp
            ratio_err = np.sqrt(syst_hyp*syst_hyp+stat_hyp*stat_hyp)
            fit = fit_val[i_cent][4]
            if fixmuQ:
                fit = fit_val_fixmuQ[i_cent][4]
            if fixmuQ and nopions:
                fit = fit_val_fixmuQ_nopions[i_cent][4]
        if i_part == 4:
            ratio = ratio_triton
            ratio_err = np.sqrt(syst_triton*syst_triton+stat_triton*stat_triton)
            fit = fit_val[i_cent][3]
            if fixmuQ:
                fit = fit_val_fixmuQ[i_cent][3]
            if fixmuQ and nopions:
                fit = fit_val_fixmuQ_nopions[i_cent][3]
        if i_part == 5:
            ratio = ratio_he3
            ratio_err = np.sqrt(syst_he3*syst_he3+stat_he3*stat_he3)
            fit = fit_val[i_cent][2]
            if fixmuQ:
                fit = fit_val_fixmuQ[i_cent][2]
            if fixmuQ and nopions:
                fit = fit_val_fixmuQ_nopions[i_cent][2]
        if ratio > 0:
            hNSigmaRatioFitParticle.SetBinContent(i_part+1,(ratio-fit)/ratio_err)
            hNSigmaRatioFitParticle.SetBinError(i_part+1,0)
        else:
            hNSigmaRatioFitParticle.SetBinContent(i_part+1,-99999.)
            hNSigmaRatioFitParticle.SetBinError(i_part+1,0)
    hNSigmaRatioFitParticle.SetLineColor(centrality_colors[i_cent])
    hNSigmaRatioFitParticle.SetMarkerColor(centrality_colors[i_cent])
    hNSigmaRatioFitParticle.SetTitle("")
    hNSigmaRatioFitParticle.GetXaxis().SetTickSize(0.07)
    hNSigmaRatioFitParticle.GetYaxis().SetTickSize(0.035)
    hNSigmaRatioFitParticle.GetYaxis().SetTitleSize(0.13)
    hNSigmaRatioFitParticle.GetYaxis().SetLabelSize(0.12)
    hNSigmaRatioFitParticle.GetYaxis().SetTitleOffset(0.4)
    hNSigmaRatioFitParticle.GetXaxis().SetLabelOffset(0.01)
    # hNSigmaRatioFitParticle.SetLineWidth(1)
    hNSigmaRatioFitParticle.GetYaxis().SetTitle("pull")#"#frac{data - fit}{#sigma_{data}}")
    hNSigmaRatioFitParticle.GetYaxis().CenterTitle()
    hNSigmaRatioFitParticle.GetYaxis().SetRangeUser(-4.7,4.7)
    # gHijingNsigma = ROOT.TGraphErrors(gHijingPredictions)
    # gHijingNsigma.SetPointY(0,(gHijingNsigma.GetPointY(0)-fit_expo.Eval(0,1))/np.sqrt(syst_pion*syst_pion+stat_pion*stat_pion))
    # gHijingNsigma.SetPointY(1,(gHijingNsigma.GetPointY(1)-fit_expo.Eval(3,0.5))/np.sqrt(syst_proton*syst_proton+stat_proton*stat_proton))
    # hLineZero = ROOT.TLine(0.,0.,6.,0.)
    # hLineZero.SetLineColor(ROOT.kBlack)
    # hLineZero.SetLineStyle(ROOT.kDashed)
    # hLineOnePos = ROOT.TLine(0.,1.,6.,1.)
    # hLineOnePos.SetLineColor(ROOT.kBlack)
    # hLineOnePos.SetLineStyle(ROOT.kDashed)
    # hLineOneNeg = ROOT.TLine(0.,-1.,6.,-1.)
    # hLineOneNeg.SetLineColor(ROOT.kBlack)
    # hLineOneNeg.SetLineStyle(ROOT.kDashed)
    # hLineTwoPos = ROOT.TLine(0.,2.,6.,2.)
    # hLineTwoPos.SetLineColor(ROOT.kBlack)
    # hLineTwoPos.SetLineStyle(ROOT.kDashed)
    # hLineTwoNeg = ROOT.TLine(0.,-2.,6.,-2.)
    # hLineTwoNeg.SetLineColor(ROOT.kBlack)
    # hLineTwoNeg.SetLineStyle(ROOT.kDashed)
    hLineZero = ROOT.TLine(0.,0.,6.,0.)
    hLineZero.SetLineColor(ROOT.kBlack)
    hLineZero.SetLineStyle(ROOT.kDashed)
    hLineOnePos = ROOT.TLine(0.,2.,6.,2.)
    hLineOnePos.SetLineColor(ROOT.kBlack)
    hLineOnePos.SetLineStyle(ROOT.kDashed)
    hLineOneNeg = ROOT.TLine(0.,-2.,6.,-2.)
    hLineOneNeg.SetLineColor(ROOT.kBlack)
    hLineOneNeg.SetLineStyle(ROOT.kDashed)
    hLineTwoPos = ROOT.TLine(0.,4.,6.,4.)
    hLineTwoPos.SetLineColor(ROOT.kBlack)
    hLineTwoPos.SetLineStyle(ROOT.kDashed)
    hLineTwoNeg = ROOT.TLine(0.,-4.,6.,-4.)
    hLineTwoNeg.SetLineColor(ROOT.kBlack)
    hLineTwoNeg.SetLineStyle(ROOT.kDashed)

    hRatiosParticle.GetYaxis().SetRangeUser(0.45,1.3)
    text_alice = ROOT.TLatex(0.174185,0.17764,"ALICE")
    text_alice.SetNDC(True)
    text_alice.SetTextFont(43)
    text_alice.SetTextSize(28)
    text_energy = ROOT.TLatex(0.175439,0.0732919,"Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV")
    text_energy.SetNDC(True)
    text_energy.SetTextFont(43)
    text_energy.SetTextSize(25)
    text_centrality = ROOT.TLatex(0.18,0.8,f"{cent[0]}-{cent[1]}%")
    text_centrality.SetNDC(True)
    text_centrality.SetTextFont(43)
    text_centrality.SetTextSize(25)
    gRatiosParticle = ROOT.TGraphErrors(hRatiosParticle)
    gRatiosParticleFit = ROOT.TGraphErrors(hRatiosParticleFit)
    for i_part in range(0,6):
        gRatiosParticle.SetPointError(i_part,0,hRatiosParticle.GetBinError(i_part+1))
        gRatiosParticleFit.SetPointError(i_part,0.3,0)
    hRatiosParticle.SetLineColor(ROOT.kWhite)
    hRatiosParticle.SetMarkerColor(ROOT.kWhite)
    hRatiosParticle.Draw("")
    hRatiosParticle.GetYaxis().SetDecimals()
    hRatiosParticle.GetXaxis().SetDecimals()
    gRatiosParticle.SetMarkerStyle(20)
    gRatiosParticle.SetMarkerSize(1.4)
    gRatiosParticle.SetLineWidth(1)
    gRatiosParticle.SetLineColor(centrality_colors[i_cent])
    gRatiosParticle.SetMarkerColor(centrality_colors[i_cent])
    gRatiosParticleFit.SetMarkerSize(0)
    gRatiosParticleFit.SetLineWidth(2)
    gRatiosParticleFit.SetLineColor(ROOT.kBlack)
    # gHijingPredictions.SetLineColor(ROOT.kGreen+1)
    # gHijingPredictions.SetLineWidth(2)
    gRatiosParticleFit.Draw("e same")
    # gHijingPredictions.Draw("pe same")
    gRatiosParticle.Draw("pe same")
    leg_ratios_particle.AddEntry(gRatiosParticle,"data","pe")
    leg_ratios_particle.AddEntry(gRatiosParticleFit,"fit")

    muBmuQ_list = muBmuQ(fixmuQ,nopions,i_cent)
    text_mu_b = ROOT.TLatex(1,0,"#it{#mu}_{B} = "+"{:.2f}".format(muBmuQ_list[0])+" #pm "+"{:.2f}".format(muBmuQ_list[1])+"(unc.) #pm "+"{:.2f}".format(muBmuQ_list[2])+"(corr.) MeV")
    text_mu_b.SetName("mu_b")
    text_mu_q = ROOT.TLatex(1.2,0,"#it{#mu}_{Q} = "+"{:.2f}".format(muBmuQ_list[3])+" #pm "+"{:.2f}".format(muBmuQ_list[4])+"(unc.) #pm "+"{:.2f}".format(muBmuQ_list[5])+"(corr.) MeV")
    if fixmuQ:
        text_mu_q = ROOT.TLatex(1.2,0,"#it{#mu}_{Q} = "+"{:.2f}".format(muBmuQ_list[3])+" MeV")
    text_mu_q.SetName("mu_q")
    text_chi2.SetName("chi2")
    text_centrality.SetName("centrality")
    # leg_ratios_particle.AddEntry(gHijingPredictions,"HIJING")
    text_mu_b.SetTextFont(43)
    text_mu_q.SetTextFont(43)
    text_chi2.SetTextFont(43)
    text_mu_b.SetTextSize(10)
    text_mu_q.SetTextSize(10)
    text_chi2.SetTextSize(10)
    text_mu_b.SetNDC(True)
    text_mu_q.SetNDC(True)
    text_chi2.SetNDC(True)
    text_mu_b.SetX(0.18)
    text_mu_b.SetY(0.16)
    text_mu_q.SetX(0.18)
    text_mu_q.SetY(0.07)
    text_chi2.SetX(0.18)
    text_chi2.SetY(0.24)
    text_mu_b.Draw("same")
    text_mu_q.Draw("same")
    text_chi2.Draw("same")
    #text_alice.Draw("same")
    #text_energy.Draw("same")
    text_centrality.Draw("same")
    leg_ratios_particle.Draw("same")
    pad2.cd()
    hNSigmaRatioFitParticle.SetMarkerStyle(20)
    hNSigmaRatioFitParticle.SetMarkerSize(1.4)
    hNSigmaRatioFitParticle.Draw("axis")
    hLineZero.Draw("same")
    hLineOnePos.Draw("same")
    hLineOneNeg.Draw("same")
    hLineTwoPos.Draw("same")
    hLineTwoNeg.Draw("same")
    # gHijingNsigma.SetLineColor(ROOT.kGreen+1)
    # gHijingNsigma.SetLineWidth(2)
    hNSigmaRatioFitParticle.Draw("psame")
    # gHijingNsigma.Draw("pe same")
    cRatiosParticle.Write()
    cRatiosParticle.Print(f"SHMfitToRatios_{cent[0]}_{cent[1]}{fixmuQ_nopions}.pdf")
    gRatiosParticle.Write()
    gRatiosParticleFit.Write()
    hNSigmaRatioFitParticle.Write()

    gMuBCent.AddPoint(n_part[i_cent],muBmuQ_list[0])
    gMuBCent.SetPointError(gMuBCent.GetN()-1,n_part_err[i_cent],muBmuQ_list[1])
    gMuBCentCorrUnc.AddPoint(n_part[i_cent],muBmuQ_list[0])
    gMuBCentCorrUnc.SetPointError(gMuBCent.GetN()-1,.5,muBmuQ_list[2])
    gMuBCentPolarityUnc.AddPoint(n_part[i_cent],muBmuQ_list[0])
    gMuBCentPolarityUnc.SetPointError(gMuBCent.GetN()-1,.5,muBmuQ_list[6])
    gMuQCent.AddPoint(n_part[i_cent],muBmuQ_list[3])
    gMuQCent.SetPointError(gMuBCent.GetN()-1,0.,muBmuQ_list[4])
    gMuQCentCorrUnc.AddPoint(n_part[i_cent],muBmuQ_list[3])
    gMuQCentCorrUnc.SetPointError(gMuQCent.GetN()-1,.5,muBmuQ_list[5])
    gMuQCentPolarityUnc.AddPoint(n_part[i_cent],muBmuQ_list[3])
    gMuQCentPolarityUnc.SetPointError(gMuQCent.GetN()-1,.5,muBmuQ_list[7])
    # print(f'mu_I_error = {mu_I_error}')

gMuBCent.SetName("MuBCent")
gMuBCentCorrUnc.SetName("gMuBCentCorrUnc")
gMuBCentPolarityUnc.SetName("gMuBCentPolarityUnc")

legMuBCent = ROOT.TLegend(0.18,0.15,0.60,0.28)
gPublishedResult = ROOT.TGraphErrors()
gPublishedResult.SetName("NaturePoint")
gPublishedResult.AddPoint(5,0.7)
gPublishedResult.SetPointError(0,5,3.8)
gPublishedResult.SetMarkerStyle(0)
gPublishedResult.SetMarkerSize(0)
gPublishedResult.SetLineWidth(0)
gPublishedResult.SetFillStyle(3345)
gPublishedResult.SetFillColor(ROOT.kBlack)
text_alice = ROOT.TLatex(0.18,0.8,"ALICE")
text_alice.SetTextFont(43)
text_alice.SetTextSize(28)
text_alice.SetNDC()
text_energy = ROOT.TLatex(0.18,0.73,"Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV")
text_energy.SetTextFont(43)
text_energy.SetTextSize(25)
text_energy.SetNDC()
text_uncorr_uncertainty = ROOT.TLatex(0.18,0.66,"Uncorrelated uncertainties only")
text_uncorr_uncertainty.SetTextFont(43)
text_uncorr_uncertainty.SetNDC()
cMuBCent = ROOT.TCanvas("cMuBCent","cMuBCent")
cMuBCent.cd()
gMuBCent.SetMarkerStyle(20)
gMuBCent.SetMarkerSize(1.2)
gMuBCent.SetMarkerColor(ROOT.kRed)
gMuBCent.SetLineColor(ROOT.kRed)
gMuBCent.SetFillColor(0)
gMuBCent.SetFillStyle(0)
gMuBCentCorrUnc.SetMarkerColor(ROOT.kRed-9)
gMuBCentCorrUnc.SetLineWidth(0)
gMuBCentCorrUnc.SetMarkerSize(0)
gMuBCentCorrUnc.SetLineColor(0)
gMuBCentCorrUnc.SetFillColor(ROOT.kRed-9)
hMuBCentFrame.GetYaxis().SetRangeUser(-3.6,5.0)
hMuBCentFrame.GetYaxis().SetNdivisions(510)
hMuBCentFrame.GetYaxis().SetDecimals()
hMuBCentFrame.GetXaxis().SetDecimals()
hMuBCentFrame.Draw("axis")
gPublishedResult.Draw("pe5same")
gMuBCentCorrUnc.Draw("pe5same")
gMuBCentPolarityUnc.Draw("pesame")
gMuBCent.Draw("pe5same")
legMuBCent.SetNColumns(2)
legMuBCent.SetTextFont(43)
legMuBCent.SetTextSize(22)
legMuBCent.AddEntry(gMuBCent,"Uncorr. uncert.", "pf")
legMuBCent.AddEntry(gMuBCentCorrUnc,"Corr. uncert.")
legMuBCent.AddEntry(gPublishedResult,"SHM fit, #it{Nature} #bf{561}, 321-330 (2018)")
legMuBCent.Draw("same")
text_alice.Draw("same")
text_energy.Draw("same")
cMuBCent.Write()
cMuBCent.Print("MuBvsNpart.eps")

cB = ROOT.TCanvas("cMuBuncorr","cMuBuncorr")
cB.cd()
gMuBCent.Fit("pol0")
gMuBCent.GetYaxis().SetTitle("#it{#mu}_{B} (MeV)")
gMuBCent.GetXaxis().SetTitle("Centrality (%)")
gMuBCent.Draw("ape")


cQ = ROOT.TCanvas("cMuQuncorr","cMuQuncorr")
cQ.cd()
gMuQCent.SetName("MuQCent")
gMuQCentCorrUnc.SetName("MuQCentCorrUnc")
gMuQCentPolarityUnc.SetName("MuQCentPolarity")
gMuQCent.SetMarkerStyle(20)
gMuQCent.SetMarkerSize(1.2)
gMuQCent.SetMarkerColor(ROOT.kRed)
gMuQCent.SetLineColor(ROOT.kRed)
gMuQCent.SetFillColor(0)
gMuQCent.SetFillStyle(0)
gMuQCent.Fit("pol0")
gMuQCent.GetYaxis().SetTitle("#it{#mu}_{Q} (MeV)")
gMuQCent.GetXaxis().SetTitle("Centrality (%)")
gMuQCentPolarityUnc.Draw("ape5")
gMuQCentCorrUnc.Draw("pesame")
gMuQCent.Draw("pesame")

cB.Write()
cQ.Write()
gMuBCent.Write()


file_out.Close()