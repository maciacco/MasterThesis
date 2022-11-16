constexpr float minpt = 0;
constexpr float maxpt = 6;
constexpr float miny[] = {0.65,-2.3};
constexpr float maxy[] = {1.28,2.3};

constexpr double sx[2]{0.355,1.-2.*sx[0]};
constexpr double sy[2]{0.7,1.-sy[0]};
constexpr double fx = sx[1];
constexpr double fy = 0.6;
double global_y = 1.;
double global_x = 0.;
const char *titles[2]={"; ;Ratio", "; ;pull"};
const double margin_val = 0.02;

//const char *particle_ratios[] = {"#frac{#pi^{-}}{#pi^{+}}","#frac{#bar{#Omega}^{+}}{#Omega^{-}}","#frac{#bar{p}}{p}","#frac{^{3}#bar{He}}{^{3}He}","#frac{{}_{#bar{#Lambda}}^{3}#bar{H}}{^{3}_{#Lambda}H}","#frac{^{3}#bar{H}}{^{3}H}"};
//const char *particle_ratios[] = {"#pi^{-} / #pi^{+}","#bar{#Omega}^{+} / #Omega^{-}","#bar{p} / p","^{3}#bar{He} / ^{3}He","{}_{#bar{#Lambda}}^{3}#bar{H} / ^{3}_{#Lambda}H","^{3}#bar{H} / ^{3}H"};
const char *particle_ratios[] = {"#pi","p","#Omega","{}^{3}He","{}^{3}_{#Lambda}H","{}^{3}H"};
constexpr int nRow=2;
constexpr int nCol=3;

std::array<TPad*,6> CreatePads(TCanvas* &cv)
{
  if (!cv) cv = new TCanvas;
  std::array<TPad*,6> pads{nullptr};

  for (int iP = 0; iP < 6; ++iP) {
    const int col = iP % nCol;
    const int row = iP / nCol;

    const int top = (row == 0);
    const int bot = (row == (nRow - 1));
    const int left = col==0;
    const int right = col==(nCol-1);
    const int center = !(left || right);

    //std::cout << sx[0] * col << "\t" << global_y - sy[bot] << "\t" << sx[0] + col * (1. - sx[0]) << "\t" << global_y << std::endl;
    pads[iP] = new TPad(Form("ratio%i",iP),"",global_x, std::abs(global_y - sy[bot]), global_x + sx[center], global_y);
    if (col == nCol - 1) global_y -= sy[bot];
    global_x = (right) ? 0. : global_x + sx[center];

    pads[iP]->SetRightMargin(right * (1 - fx / sx[center]));
    pads[iP]->SetLeftMargin(left * (1 - fx / sx[center]));
    // if (row == 0){
    //   double f[2] = {0.25, 0.25};
    //   if (left) {
    //     pads[iP]->SetRightMargin(f[0] * (1 - fx / sx_0[center]) * 1.10);
    //     pads[iP]->SetLeftMargin(1.1 * (1 - fx / sx_0[center]) * 1.10);
    //   }
    //   else if (right){
    //     pads[iP]->SetRightMargin(1.1 * (1 - fx / sx_0[center]) * 1.10);
    //     pads[iP]->SetLeftMargin(f[1] * (1 - fx / sx_0[center]) * 1.10);
    //   }
    //   else {
    //     pads[iP]->SetRightMargin(f[0] * (1 - fx / sx_0[center]) * 1.10);
    //     pads[iP]->SetLeftMargin(f[1] * (1 - fx / sx_0[center]) * 1.10);
    //   }
    // }
    if (row==0){
      pads[iP]->SetRightMargin(1 - (fx - margin_val) / sx[center]);
      pads[iP]->SetLeftMargin(1 - (fx - margin_val) / sx[center]);
      if (left)
        pads[iP]->SetRightMargin(margin_val / sx[center]);
      if (right)
        pads[iP]->SetLeftMargin(margin_val / sx[center]);
    }
    pads[iP]->SetTopMargin(0.0);
    if (row == 0) pads[iP]->SetTopMargin(0.1);
    pads[iP]->SetBottomMargin(bot * (1 - fy));
    if (row == 0) pads[iP]->SetBottomMargin(0.05);

    cv->cd();
    pads[iP]->Draw();
    pads[iP]->cd();
    if (row == 1)pads[iP]->SetGridy();

    if (row==1){
      TH2F *rframe = new TH2F(Form("rframe%i",iP),titles[row],6,minpt,maxpt,100,miny[row],maxy[row]);
      for (int i_part=0;i_part<6;i_part++)
        rframe->GetXaxis()->SetBinLabel(i_part+1,particle_ratios[i_part]);

      rframe->GetYaxis()->CenterTitle();
      rframe->GetYaxis()->SetTickLength(0.01 / sx[center] * (1+bot*fy));
      rframe->GetYaxis()->SetTitleSize(30);
      rframe->GetYaxis()->SetTitleFont((!col) * 43);
      rframe->GetYaxis()->SetTitleOffset(1.5);
      rframe->GetYaxis()->SetLabelOffset(0.01);
      rframe->GetYaxis()->SetNdivisions(505);
      rframe->GetYaxis()->SetDecimals(1);
      if (row==1) rframe->GetYaxis()->SetNdivisions(5);
      rframe->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
      rframe->GetYaxis()->SetLabelSize((!col) * 20);

      /* rframe->GetXaxis()->CenterTitle(); */
      rframe->GetXaxis()->SetTickLength(0.004 / sx[1-center] / sy[bot]);
      rframe->GetXaxis()->SetLabelSize(30);
      rframe->GetXaxis()->SetLabelFont(43);
      rframe->GetXaxis()->SetLabelOffset(0.04);
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
double point_bins[6][2]={{1,11},{1,21},{4,16},{9,11},{10,16},{10,6}};

void plot_fits_3d() {
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kDarkBodyRadiator);
  const float leftmargin = 0.16;
  const float rigthmargin = 0.04;

  const float deltaRatio = 0.89;
  const float kMarkerSize = 0.9;

  const float sizeRatioPadY = 0.15;
  const float sizeLastPad = (1. - 3 * sizeRatioPadY) * 0.5f;
  const float titleYoffset = 2.3;

  TLine l;
  l.SetLineStyle(kDashed);
  l.SetLineColor(kBlack);

  TCanvas *cv = new TCanvas("cv","cv",2100,700);
  auto pads = CreatePads(cv);

  TLatex text;
  text.SetTextFont(63);
  pads[0]->cd();
  double mean_y = 0.5*(miny[0]+maxy[0]);
  double half_width_y = 0.5*(maxy[0]-miny[0]);
  double mean_x = 0.5*(minpt+maxpt);
  double half_width_x = 0.5*(maxpt-minpt);

  const string labels[3]{"0-5%","5-10%","30-50%"};
  const string names[3]{"0_5","5_10","30_50"};
  const string chi2[3]{"1.72","3.70","2.55"};
  TFile input("FinalPlot3D.root");
  TH2D *h3D[3];
  TGraph2DErrors *g3D[3];
  TF2 *fFit[3];
  TH1D* hRatio_[3];
  TLegend *leg;
  TH1D* hRatio[3];
  for (int iP = 0; iP < 3; ++iP) {
    pads[iP]->cd();
    pads[iP]->SetTheta(40);
    pads[iP]->SetPhi(30);
    text.SetTextSize(18);
    text.DrawText(mean_x-0.85*half_width_x,mean_y+0.8*half_width_y,labels[iP].data());
    auto ff = (TF2*)input.Get(Form("fit_expo_%s",names[iP].data()));
    h3D[iP] = (TH2D*)input.Get(Form("fRatio_vs_b_%s",names[iP].data()));
    fFit[iP] = (TF2*)input.Get(Form("fit_expo_%s",names[iP].data()));
    hRatio_[iP] = (TH1D*)input.Get(Form("hNSigmaRatioFitParticle_%s",names[iP].data()));
    hRatio[iP] = new TH1D(*hRatio_[iP]);
    int map[]={1,2,0,5,3,4};
    for (int ip=0; ip<6; ++ip){
      hRatio[iP]->SetBinContent(ip+1,hRatio_[iP]->GetBinContent(map[ip]+1));
    }
    fFit[iP]->SetLineWidth(2);
    fFit[iP]->SetMinimum(0.67);
    fFit[iP]->SetMaximum(1.13);
    g3D[iP]=new TGraph2DErrors();
    for (int ip=0; ip<6; ++ip){
      double x_bin = h3D[iP]->GetXaxis()->GetBinCenter(point_bins[ip][0]);
      double y_bin = h3D[iP]->GetYaxis()->GetBinCenter(point_bins[ip][1]);
      g3D[iP]->AddPoint(x_bin,y_bin,h3D[iP]->GetBinContent(point_bins[ip][0],point_bins[ip][1]));
      g3D[iP]->SetPointError(ip,0.,0.,h3D[iP]->GetBinError(point_bins[ip][0],point_bins[ip][1]));
      std::cout<<"x = "<<x_bin<<"; y = "<<y_bin<<"; ERROR = "<<h3D[iP]->GetBinError(point_bins[ip][0],point_bins[ip][1])<<std::endl;
    }
    h3D[iP]->SetMarkerStyle(20);
    h3D[iP]->SetMarkerSize(2);
    h3D[iP]->SetMarkerColor(colors[iP]);
    h3D[iP]->SetLineColor(colors[iP]);
    h3D[iP]->GetXaxis()->SetRangeUser(-0.5,9.5);
    h3D[iP]->GetYaxis()->SetRangeUser(-0.5,6.5);
    h3D[iP]->GetZaxis()->SetRangeUser(0.67,1.13);
    h3D[iP]->GetZaxis()->SetNdivisions(3);
    h3D[iP]->GetYaxis()->SetNdivisions(5);
    if (iP>0) h3D[iP]->GetZaxis()->SetTitle("");
    h3D[iP]->Draw("pe");
    if (iP==0){
      text.SetTextFont(43);
      text.SetTextSize(32);
      text.DrawLatexNDC(0.1,0.92,"ALICE, Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
    }
    fFit[iP]->SetLineColor(kOrange+1);
    fFit[iP]->SetLineWidth(1);
    fFit[iP]->Draw("surf1 same a fb bb");
    h3D[iP]->Draw("pe same");
    text.SetTextSize(25);
    text.DrawLatex(-0.15,-0.5,Form("#chi^{2}/NDF = %.1f/%d",ff->GetChisquare(),ff->GetNDF()));
    text.DrawLatex(0.5,0.8,labels[iP].data());
    // if (iP==0){
    //   leg=new TLegend(0.22808,0.172684,0.611348,0.273155);
    //   leg->AddEntry(fFit[iP],"fit");
    //   leg->SetTextFont(43);
    //   leg->SetTextSize(18);
    //   leg->Draw("same");
    // }
    pads[iP+3]->cd();
    hRatio[iP]->Draw("psame");
    hRatio[iP]->SetMarkerColor(colors[iP]);
    hRatio[iP]->SetMarkerStyle(20);
    hRatio[iP]->SetMarkerSize(2);
  }
  cv->SaveAs("SHMfitToRatios_all.eps");

}
