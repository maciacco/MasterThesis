void reweight_events(){
  ROOT::EnableImplicitMT(32);
  ROOT::RDataFrame df("XiOmegaTree","/data/mciacco/Omega_PbPb/merge_omegas/AnalysisResults_omega.root");
  TFile inFile("/data/mciacco/Omega_PbPb/AnalysisResults_18.root");
  TH1D *cent=(TH1D*)inFile.Get("Centrality_selected");
  auto reweight = [cent](unsigned char centrality){
    int i_bin = cent->FindBin((double)centrality+0.0001);
    double h_val = cent->GetBinContent(i_bin);
    double threshold = 0;
    if (centrality>=10&&centrality<30) threshold = cent->GetBinContent(cent->FindBin(20));
    if (centrality>=50) threshold = cent->GetBinContent(cent->FindBin(60));
    bool cut = gRandom->Rndm()*h_val < threshold;
    if (centrality<10||(centrality>=30&&centrality<50)) cut=true;
    return cut;
  };
  df.Filter(reweight,{"centrality"}).Snapshot("XiOmegaTree","/data/mciacco/Omega_PbPb/merge_omegas/AnalysisResults_omega_reweight.root");
}