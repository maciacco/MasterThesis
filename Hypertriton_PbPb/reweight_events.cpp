void reweight_events(){
  ROOT::EnableImplicitMT(32);
  ROOT::RDataFrame df("DataTable","/data/mciacco/Hypertriton_PbPb/DataTable_18.root");
  TFile out_cent("out_cent.root","recreate");
  TF1 fit_10_30("fit_10_30","pol3");
  TF1 fit_50_90("fit_50_90","pol3");
  auto h_cent_ = df.Histo1D({"h","h",100,0,100},"centrality");
  TH1D* h_cent = new TH1D(*h_cent_);
  h_cent->Write();
  h_cent->Fit("fit_10_30","IRM+","",12,28);
  h_cent->Fit("fit_50_90","IRM+","",52,70);
  TFile inFile("/data/mciacco/Hypertriton_PbPb/AnalysisResults_18.root");
  TH1D *cent=(TH1D*)inFile.Get("Centrality_selected");
  auto reweight = [h_cent](float centrality_){
    int centrality = (int)centrality_;
    int i_bin = h_cent->FindBin(centrality_+0.0001);
    double h_val = h_cent->GetBinContent(i_bin);
    double threshold = 0;
    if (centrality>=10&&centrality<30) threshold = h_cent->GetFunction("fit_10_30")->Eval(h_cent->GetBinCenter(i_bin));
    if (centrality>=50&&centrality<70) threshold = h_cent->GetFunction("fit_50_90")->Eval(h_cent->GetBinCenter(i_bin));
    bool cut = (gRandom->Rndm()*h_val) < threshold;
    if (centrality<10||(centrality>=30&&centrality<50)||centrality>=70) cut=true;
    return cut;
  };
  df.Filter(reweight,{"centrality"}).Snapshot("DataTable","/data/mciacco/Hypertriton_PbPb/DataTable_18_reweight.root");

  out_cent.Close();
}