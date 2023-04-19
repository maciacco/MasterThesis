constexpr float minpt = 0;
constexpr float maxpt = 6;
constexpr float miny[] = {0.3,-2.3};
constexpr float maxy[] = {1.3,2.3};

constexpr double sx[2]{0.355,1.-2*sx[0]};
constexpr double sy[2]{0.35,0.15};
constexpr double sy_top[2]{0.,-0.0525};
constexpr double sy_bot[2]{0.,+0.0525};
constexpr double fx = sx[1];
constexpr double fy = 0.85;
double global_y = 1.;
double global_x = 0.;
const char *titles[2]={"; ;antimatter / matter", "; ;pull"};

const char *particle_ratios[] = {"#Omega","#pi","p","{}^{3}_{#Lambda}H","{}^{3}H","{}^{3}He"};
constexpr int nRow=4;
constexpr int nCol=3;

std::array<TPad*,12> CreatePads(TCanvas* &cv)
{
  if (!cv) cv = new TCanvas;
  std::array<TPad*,12> pads{nullptr};

  for (int iP = 0; iP < 12; ++iP) {
    const int col = iP % nCol;
    const int row = iP / nCol;

    int top = (row == 0);
    int bot_global = (row == 3);
    int top_low = (row == 2);
    int bot = (row == (nRow - 1) || row == (nRow-3));
    int left = col==0;
    int right = col==(nCol-1);
    int center = !(left || right);

    std::cout << sx[0] * col << "\t" << global_y - sy[bot] - sy_top[top_low] - sy_bot[bot_global] << "\t" << sx[0] + col * (1. - sx[0]) << "\t" << global_y << std::endl;
    if (!((row== 1 && col==2)||(row==0 && col== 2))) pads[iP] = new TPad(Form("ratio%i",iP),"",global_x, std::abs(global_y - sy[bot] - sy_top[top_low] - sy_bot[bot_global]), global_x + sx[center], global_y);
    else if (row== 0 && col==2) pads[iP] = new TPad(Form("ratio%i",iP),"",global_x, std::abs(global_y - sy[0] - sy[1]), global_x + sx[center], global_y);
    else pads[iP] = new TPad(Form("ratio%i",iP),"",0., 0., 0.5, 0.5);
    if (col == nCol - 1) global_y -= (sy[bot] + sy_top[top_low]);
    global_x = (right) ? 0. : global_x + sx[center];

    pads[iP]->SetRightMargin(right * (1 - fx / sx[center]));
    pads[iP]->SetLeftMargin(left * (1 - fx / sx[center]));
    pads[iP]->SetTopMargin(0.0);
    if (row == 0) {
      pads[iP]->SetTopMargin(top * (1 - fy));
      if (col == 2)pads[iP]->SetTopMargin(top * (1 - 0.895));
    }
    pads[iP]->SetBottomMargin((bot_global) * (1 - 0.74074074));

    cv->cd();
    pads[iP]->Draw();
    pads[iP]->cd();
    if (!(row== 1 && col==2)){
      if (row == 1 || row == 3)pads[iP]->SetGridy();

      TH2F *rframe = new TH2F(Form("rframe%i",iP),titles[row%2],6,minpt,maxpt,100,miny[bot],maxy[bot]);
      for (int i_part=0;i_part<6;i_part++){
        rframe->GetXaxis()->SetBinLabel(i_part+1,particle_ratios[i_part]);
        rframe->GetXaxis()->ChangeLabel(i_part+1,-1.,-1.,21);
      }

      rframe->GetYaxis()->CenterTitle();
      rframe->GetYaxis()->SetTickLength(0.01 / sx[center] / (1+bot*fy));
      if (row == 2)rframe->GetYaxis()->SetTickLength(0.01 / sx[center] / (1+0.2));
      if (row == 3)rframe->GetYaxis()->SetTickLength(0.01 / sx[center] / (1+0.2));
      if (row==0 && col == 2) rframe->GetYaxis()->SetTickLength(0);
      rframe->GetYaxis()->SetTitleSize(27);
      rframe->GetYaxis()->SetTitleFont((!col) * 43);
      rframe->GetYaxis()->SetTitleOffset(1.2);
      rframe->GetYaxis()->SetLabelOffset(0.01);
      rframe->GetYaxis()->SetNdivisions(505);
      rframe->GetYaxis()->SetDecimals(1);
      if (row==1 || row==3) rframe->GetYaxis()->SetNdivisions(5);
      rframe->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
      rframe->GetYaxis()->SetLabelSize((!col) * 15);

      /* rframe->GetXaxis()->CenterTitle(); */
      rframe->GetXaxis()->SetTickLength(0.004 / sx[1-center] / (sy[bot] + sy_top[top_low] + sy_bot[bot_global]));
      if (row==0 && col == 2) rframe->GetXaxis()->SetTickLength(0);
      rframe->GetXaxis()->SetLabelSize(27);
      rframe->GetXaxis()->SetLabelFont(43);
      rframe->GetXaxis()->SetLabelOffset(0.07);
      rframe->GetXaxis()->SetNdivisions(505);
      rframe->GetXaxis()->SetDecimals(1);
      std::cout<<rframe->GetXaxis()->GetLabelSize()<<" "<<rframe->GetXaxis()->GetTickLength()<<" "<<rframe->GetYaxis()->GetLabelSize()<<" "<<rframe->GetYaxis()->GetTickLength()<<std::endl;

      rframe->Draw("col");
    }
  }
  return pads;
}

Color_t colors[]={kOrange+7, kAzure+4, kTeal+4};
double x_limits[3][2]={{0.7,1.6},{1.0,3.0},{2.,8.}};
int i_particle=2;
const char *format_fit_results[] = {"R=%.3f #pm %.3f","R=%.3f #pm %.3f","R=%.2f #pm%.2f (stat.) #pm%.2f (syst.) #pm%.2f (mat.)"};
const char *format_fit_results_more[] = {"R=%.4f #pm %.4f"};
double text_position_up_right[3][2]={{1.2,1.023},{2.12,1.063},{1.9,1.4}};
double text_position_low_right[3][2]={{1.2,1.017},{2.12,1.050},{1.9,1.3}};
double text_position_up_left[3][2]={{0.71,1.023},{1.06,1.063},{8.2,1.4}};
double text_position_low_left[3][2]={{0.71,1.017},{1.06,1.050},{1.9,1.2}};
const char* cent_labels[] = {"0_5","5_10","10_30","30_50","50_90"};
const char* ax_labels[] = {"#Omega","#pi","p","{}^{3}_{#Lambda}H","{}^{3}H","{}^{3}He"};


void plot_fits_try(){
  gStyle->SetOptStat(0);
  TFile in_file("FinalPlot3D_new_2.root");

  TCanvas *c=new TCanvas("c","c",1300,1000);
  auto pads=CreatePads(c);

  std::vector<TCanvas*> cc;
  for (auto label : cent_labels){
    cc.push_back((TCanvas*)in_file.Get(Form("cRatiosParticle_%s",label)));
  }

  for (int i_c{0}; i_c < sizeof(cent_labels)/sizeof(cent_labels[0]); ++i_c){
    pads[i_c]->cd();
    if (i_c > 1) pads[i_c+4]->cd();
    auto g_ratio = (TGraphErrors*)(cc[i_c]->FindObject(Form("Graph_from_hRatiosParticle_%s",cent_labels[i_c])));
    auto g_fit = (TGraphErrors*)(cc[i_c]->FindObject(Form("Graph_from_hRatiosParticleFit_%s",cent_labels[i_c])));
    auto t_cent = (TLatex*)cc[i_c]->FindObject("centrality");
    auto t_mub = (TLatex*)cc[i_c]->FindObject("mu_b");
    auto t_muq = (TLatex*)cc[i_c]->FindObject("mu_q");
    auto t_chi2 = (TLatex*)cc[i_c]->FindObject("chi2");
    g_ratio->SetMarkerSize(2.);
    g_fit->Draw("esame");
    g_ratio->Draw("samepe");
    t_cent->SetNDC(false);
    t_mub->SetNDC(false);
    t_muq->SetNDC(false);
    t_chi2->SetNDC(false);
    t_cent->SetTextFont(44);
    t_mub->SetTextFont(44);
    t_muq->SetTextFont(44);
    t_chi2->SetTextFont(44);
    t_cent->SetTextSize(23);
    t_mub->SetTextSize(23);
    t_muq->SetTextSize(23);
    t_chi2->SetTextSize(23);
    t_cent->SetX(0.3);
    t_mub->SetX(0.3);
    t_muq->SetX(0.3);
    t_chi2->SetX(0.3);
    t_cent->SetY(1.15);
    t_mub->SetY(0.5);
    t_muq->SetY(0.4);
    t_chi2->SetY(0.6);
    t_cent->Draw("same");
    t_mub->Draw("same");
    t_muq->Draw("same");
    t_chi2->Draw("same");
    auto h_pull = (TH1D*)(cc[i_c]->FindObject(Form("hNSigmaRatioFitParticle_%s",cent_labels[i_c])));
    h_pull->SetMarkerSize(2.);
    pads[i_c+3]->cd();
    std::cout<<"pad = "<<1+(i_c)*2<<std::endl;
    if (i_c > 1) pads[i_c+7]->cd();
    TGraph *g_pull = new TGraph(h_pull);
    g_pull->Draw("pesame");
    //cc[i_c]->DrawClonePad();
    auto leg = (TLegend*)c->FindObject("TPave");
    //leg->Delete();
  }

  TGraphErrors g;
  g.AddPoint(-100,-100);
  TF1 f("f","pol0");
  f.SetParameter(0,-1000);
  pads[2]->cd();
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
  t.DrawLatexNDC(0.03,0.7,"ALICE");
  t.SetTextSize(35);
  t.DrawLatexNDC(0.03,0.6,"Pb-Pb #sqrt{#it{s}_{NN}}=5.02 TeV");
  //t.DrawLatexNDC(0.15,0.6,"|y| < 0.5");
  TLegend l(0.03,0.55,0.5,0.35);
  l.SetTextSize(35);
  l.AddEntry(&g,"data","pe");
  l.AddEntry(&f,"Thermal-FIST","l");
  t.SetTextSize(25);
  t.DrawLatex(1.05, 0.64, "#it{T}_{ch}=155#pm2 MeV");
  t.DrawLatex(1.05, 0.56, "#mu_{#it{S}} constrained");
  l.SetTextFont(44);
  l.SetTextSize(25);
  l.Draw("same");

  TFile out("fit_out_prova.root","recreate");
  c->Write();
  out.Close();
  c->Print("fits.pdf");
}