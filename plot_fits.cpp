constexpr float minpt = 0;
constexpr float maxpt = 4;
constexpr float miny[] = {0.66,-2.3};
constexpr float maxy[] = {1.14,2.3};

constexpr double sx[2]{0.355,1.-2*sx[0]};
constexpr double sy[2]{0.7,1.-sy[0]};
constexpr double fx = sx[1];
constexpr double fy = 0.6;
double global_y = 1.;
double global_x = 0.;
const char *titles[2]={"; ;Ratio", "; ;#frac{data - fit}{#sigma_{data}}"};

const char *particle_ratios[] = {"#pi^{-} / #pi^{+}","#bar{p} / p","{}_{#bar{#Lambda}}^{3}#bar{H} / ^{3}_{#Lambda}H","^{3}#bar{He} / ^{3}He"};
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
    pads[iP]->SetTopMargin(0.0);
    if (row == 0) pads[iP]->SetTopMargin(0.01);
    pads[iP]->SetBottomMargin(bot * (1 - fy));

    cv->cd();
    pads[iP]->Draw();
    pads[iP]->cd();
    if (row == 1)pads[iP]->SetGridy();

    TH2F *rframe = new TH2F(Form("rframe%i",iP),titles[row],4,minpt,maxpt,100,miny[row],maxy[row]);
    for (int i_part=0;i_part<4;i_part++)
      rframe->GetXaxis()->SetBinLabel(i_part+1,particle_ratios[i_part]);

    rframe->GetYaxis()->CenterTitle();
    rframe->GetYaxis()->SetTickLength(0.01 / sx[center] * (1+bot*fy));
    rframe->GetYaxis()->SetTitleSize(20);
    rframe->GetYaxis()->SetTitleFont((!col) * 43);
    rframe->GetYaxis()->SetTitleOffset(1.5);
    rframe->GetYaxis()->SetLabelOffset(0.01);
    rframe->GetYaxis()->SetNdivisions(505);
    rframe->GetYaxis()->SetDecimals(1);
    if (row==1) rframe->GetYaxis()->SetNdivisions(5);
    rframe->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    rframe->GetYaxis()->SetLabelSize((!col) * 15);

    /* rframe->GetXaxis()->CenterTitle(); */
    rframe->GetXaxis()->SetTickLength(0.004 / sx[1-center] / sy[bot]);
    rframe->GetXaxis()->SetLabelSize(20);
    rframe->GetXaxis()->SetLabelFont(43);
    rframe->GetXaxis()->SetLabelOffset(0.04);
    rframe->GetXaxis()->SetNdivisions(505);
    rframe->GetXaxis()->SetDecimals(1);
    std::cout<<rframe->GetXaxis()->GetLabelSize()<<" "<<rframe->GetXaxis()->GetTickLength()<<" "<<rframe->GetYaxis()->GetLabelSize()<<" "<<rframe->GetYaxis()->GetTickLength()<<std::endl;

    rframe->Draw("col");
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

void plot_fits() {
  gStyle->SetOptStat(0);
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

  TCanvas *cv = new TCanvas("cv","cv",1200,480);
  auto pads = CreatePads(cv);

  TLatex text;
  text.SetTextFont(63);
  text.SetTextSize(18);
  pads[0]->cd();
  double mean_y = 0.5*(miny[0]+maxy[0]);
  double half_width_y = 0.5*(maxy[0]-miny[0]);
  double mean_x = 0.5*(minpt+maxpt);
  double half_width_x = 0.5*(maxpt-minpt);
  text.DrawText(mean_x-0.85*half_width_x,mean_y-0.4*half_width_y,"ALICE Preliminary");

  pads[0]->cd();
  text.SetTextFont(43);
  text.DrawLatex(mean_x-0.85*half_width_x,mean_y-0.6*half_width_y,"Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");

  const string labels[3]{"0-5%","5-10%","30-50%"};
  const string names[3]{"0_5","5_10","30_50"};
  const string chi2[3]{"1.72","3.70","2.55"};
  TFile input("FinalPlot3D.root");
  TGraphErrors *g[3],*gFit[3];
  TH1D* hRatio[3];
  for (int iP = 0; iP < 3; ++iP) {
    pads[iP]->cd();
    text.SetTextSize(18);
    text.DrawText(mean_x-0.85*half_width_x,mean_y+0.8*half_width_y,labels[iP].data());
    auto ff = (TF2*)input.Get(Form("fit_expo_%s",names[iP].data()));
    text.DrawLatex(mean_x-0.85*half_width_x,mean_y-0.8*half_width_y,Form("#chi^{2}/NDF = %.1f/%d",ff->GetChisquare(),ff->GetNDF()));
    g[iP] = (TGraphErrors*)input.Get(Form("Graph_from_hRatiosParticle_%s",names[iP].data()));
    g[iP]->SetMarkerSize(1.0);
    gFit[iP] = (TGraphErrors*)input.Get(Form("Graph_from_hRatiosParticleFit_%s",names[iP].data()));
    hRatio[iP] = (TH1D*)input.Get(Form("hNSigmaRatioFitParticle_%s",names[iP].data()));
    gFit[iP]->SetLineWidth(2);
    gFit[iP]->Draw("esame");
    g[iP]->Draw("pesame");  
    pads[iP+3]->cd();
    hRatio[iP]->Draw("psame");
    hRatio[iP]->SetMarkerSize(1.0);
  }
  cv->SaveAs("SHMfitToRatios_all.eps");

}
