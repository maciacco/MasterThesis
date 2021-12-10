import ROOT

df_mc = ROOT.RDataFrame("LambdaTree","../data/Lambda_PbPb/mc_cols.root")
df_mc.Filter("isReconstructed").Snapshot("LambdaTree","../data/Lambda_PbPb/mc_cols_rec.root")