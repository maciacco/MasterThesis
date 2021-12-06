import ROOT
import numpy

df_mc = ROOT.RDataFrame("LambdaTree","../data/Lambda_PbPb/mc_reduced_tenth.root")
df_mc.Filter("isReconstructed").Snapshot("LambdaTree","../data/Lambda_PbPb/mc_reduced_tenth_rec.root")