const char* cent_labels[] = {"0_5","5_10","10_30","30_50","50_90"};
const char* ax_labels[] = {"#Omega","#pi","p","{}^{3}_{#Lambda}H","{}^{3}H","{}^{3}He"};


void plot_fits_theFist(){
  TFile in_file("FinalPlot3D_new.root");

  TCanvas c("c","c",1500,1000);
  c.Divide(3,2);

  std::vector<TCanvas*> cc;
  for (auto label : cent_labels){
    cc.push_back((TCanvas*)in_file.Get(Form("cRatiosParticle_%s",label)));
  }

  for (int i_c{0}; i_c < sizeof(cent_labels)/sizeof(cent_labels[0]); ++i_c){
    c.cd(i_c+1);
    if (i_c > 1) c.cd(i_c+2);
    auto h_ratio = (TH1D*)(cc[i_c]->FindObject(Form("hRatiosParticle_%s",cent_labels[i_c])));
    h_ratio->GetYaxis()->SetNdivisions(10);
    h_ratio->GetYaxis()->SetLabelFont(44);
    h_ratio->GetYaxis()->SetLabelSize(20);
    h_ratio->GetYaxis()->SetTitleFont(44);
    h_ratio->GetYaxis()->SetTitle("antimatter / matter");
    h_ratio->GetYaxis()->SetTitleSize(25);
    h_ratio->GetYaxis()->SetTitleOffset(1.2);
    auto g_ratio = (TGraphErrors*)(cc[i_c]->FindObject(Form("Graph_from_hRatiosParticle_%s",cent_labels[i_c])));
    g_ratio->SetMarkerSize(2.);
    auto h_pull = (TH1D*)(cc[i_c]->FindObject(Form("hNSigmaRatioFitParticle_%s",cent_labels[i_c])));
    h_pull->GetYaxis()->SetRangeUser(-4.9,4.9);
    for (int i_b{0}; i_b < h_pull->GetNbinsX(); ++i_b){
      h_pull->GetXaxis()->SetBinLabel(i_b+1,ax_labels[i_b]);
    }
    h_pull->GetXaxis()->SetLabelFont(44);
    h_pull->GetXaxis()->SetLabelSize(40);
    h_pull->GetXaxis()->SetLabelOffset(.02);
    h_pull->GetYaxis()->SetLabelFont(44);
    h_pull->GetYaxis()->SetLabelSize(20);
    h_pull->GetYaxis()->SetTitleFont(44);
    h_pull->GetYaxis()->SetTitleSize(25);
    h_pull->GetYaxis()->SetTitleOffset(1.1);
    h_pull->SetMarkerSize(2.);
    cc[i_c]->DrawClonePad();
    auto leg = (TLegend*)c.FindObject("TPave");
    leg->Delete();
  }

  TGraphErrors g;
  g.AddPoint(-100,-100);
  TF1 f("f","pol0");
  f.SetParameter(0,-1000);
  c.cd(3);
  g.SetMarkerStyle(20);
  g.SetMarkerSize(2.);
  g.SetLineColor(kBlack);
  g.SetMarkerColor(kBlack);
  g.Draw("pe");
  f.SetLineColor(kBlack);
  f.Draw("same");
  TLatex t;
  t.SetTextFont(44);
  t.SetTextSize(40);
  t.DrawLatexNDC(0.15,0.8,"ALICE");
  t.SetTextSize(35);
  t.DrawLatexNDC(0.15,0.7,"Pb-Pb #sqrt{#it{s}_{NN}}=5.02 TeV");
  //t.DrawLatexNDC(0.15,0.6,"|y| < 0.5");
  TLegend l(0.15,0.65,0.7,0.35);
  l.AddEntry(&g,"data","pe");
  l.AddEntry(&f,"fit (Thermal-FIST)","l");
  l.SetTextFont(44);
  l.SetTextSize(30);
  l.Draw("same");

  TFile out("fit_out_prova.root","recreate");
  c.Write();
  out.Close();
  c.Print("fits.pdf");
}