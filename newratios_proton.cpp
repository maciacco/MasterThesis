constexpr float minpt = 0.95;
constexpr float maxpt = 3.05;
constexpr float miny = 0.92;
constexpr float maxy = 1.08;
std::array<TPad*,3> CreatePads(TCanvas* &cv)
{
  if (!cv) cv = new TCanvas;
  std::array<TPad*,3> pads{nullptr};

  constexpr double sx[2]{1,1.-sx[0]};
  constexpr double sy[2]{0.365, 0.27};
  constexpr double fx = 0.9;
  constexpr double fy = 0.3;

  double global_y = 1.;
  for (int iP = 0; iP < 3; ++iP) {
    const int col = 0;
    const int row = iP;

    const int top = (row == 0);
    const int bot = row == 2;
    const int mid = !(top || bot);

    std::cout << sx[0] * col << "\t" << global_y - sy[mid] << "\t" << sx[0] + col * (1. - sx[0]) << "\t" << global_y << std::endl;
    pads[iP] = new TPad(Form("ratio%i",iP),"", 0, std::abs(global_y - sy[mid]), 1, global_y);
    global_y -= sy[mid];

    pads[iP]->SetRightMargin(0.01);
    pads[iP]->SetLeftMargin(0.15);
    pads[iP]->SetTopMargin(top * (1 - fy / sy[mid]));
    pads[iP]->SetBottomMargin(bot * (1 - fy / sy[mid]));

    cv->cd();
    pads[iP]->Draw();
    pads[iP]->cd();
    TH2F *rframe = new TH2F(Form("rframe%i",iP),";#it{p}_{T} (GeV/#it{c});#bar{p}/p;",100,minpt,maxpt,100,miny,maxy);

    rframe->GetYaxis()->CenterTitle();
    rframe->GetYaxis()->SetTickLength(0.03 * sx[0] / sx[col]);
    rframe->GetYaxis()->SetTitleSize(20);
    rframe->GetYaxis()->SetTitleFont((!col) * 43);
    rframe->GetYaxis()->SetTitleOffset(2.);
    rframe->GetYaxis()->SetNdivisions(505);
    rframe->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    rframe->GetYaxis()->SetLabelSize((!col) * 15);

    rframe->GetXaxis()->SetTickLength(0.03 * sy[0] / sy[mid]);
    rframe->GetXaxis()->SetTitleSize(bot * 20);
    rframe->GetXaxis()->SetTitleFont(43);
    rframe->GetXaxis()->SetTitleOffset(3);
    rframe->GetXaxis()->SetNdivisions(510);
    rframe->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    rframe->GetXaxis()->SetLabelSize(bot * 15);

    rframe->Draw("col");
  }
  return pads;
}

Color_t colors[]={kOrange+7, kAzure+4, kTeal+4};
double x_limits[2][2]={{0.7,1.6},{1.0,3.0}};
int i_particle=1;
const char *format_fit_results[] = {"Ratio = %.3f #pm %.3f","Ratio = %.3f #pm %.3f"};
const char *format_fit_results_more[] = {"Ratio = %.4f #pm %.4f"};
double text_position_up_right[2][2]={{1.2,1.023},{2.12,1.063}};
double text_position_low_right[2][2]={{1.2,1.017},{2.12,1.050}};
double text_position_up_left[2][2]={{0.71,1.023},{1.06,1.063}};
double text_position_low_left[2][2]={{0.71,1.017},{1.06,1.050}};

void newratios_proton() {
  gStyle->SetOptStat(0);
  const float leftmargin = 0.16;
  const float rigthmargin = 0.04;

  const float minpt = 0.4;
  const float maxpt = 8.4;
  const float deltaRatio = 0.89;
  const float kMarkerSize = 0.9;

  const float sizeRatioPadY = 0.15;
  const float sizeLastPad = (1. - 3 * sizeRatioPadY) * 0.5f;
  const float titleYoffset = 2.3;

  TLine l;
  l.SetLineStyle(kDashed);
  l.SetLineColor(kBlack);

  TCanvas *cv = new TCanvas("cv","cv",480,800);
  auto pads = CreatePads(cv);

  TLatex text;
  text.SetTextFont(63);
  text.SetTextSize(18);
  pads[0]->cd();
  text.DrawText(text_position_up_right[i_particle][0],text_position_up_right[i_particle][1],"ALICE Preliminary");

  pads[0]->cd();
  text.SetTextFont(43);
  text.DrawLatex(text_position_low_right[i_particle][0],text_position_low_right[i_particle][1],"Pb-Pb #sqrt{#it{s}_{NN}}=5.02 TeV");

  const string labels[3]{"0-5%","5-10%","30-50%"};
  const string names[3]{"0_5","5_10","30_50"};
  //text.SetTextAlign(31);
  TFile input("Proton_PbPb/out/SystematicsAllEPtNotCombined.root");
  TGraphErrors *g[3];
  TH1D* h[3];
  for (int iP = 0; iP < 3; ++iP) {
    pads[iP]->cd();
    text.DrawText(text_position_up_left[i_particle][0],text_position_up_left[i_particle][1],labels[iP].data());
    h[iP] = (TH1D*)input.Get(Form("fRatio_%s",names[iP].data()));
    g[iP]=new TGraphErrors(h[iP]);
    g[iP]->SetMarkerStyle(20);
    g[iP]->SetMarkerSize(1.);
    g[iP]->SetMarkerColor(colors[iP]);
    g[iP]->SetLineColor(colors[iP]);
    g[iP]->SetFillStyle(0);
    int iPoint{0};
    while (iPoint<g[iP]->GetN()){
      double x=g[iP]->GetPointX(iPoint);
      if (x<x_limits[i_particle][0]||x>x_limits[i_particle][1]){
        g[iP]->RemovePoint(iPoint);
      }
      else ++iPoint;
    }
    h[iP]->GetFunction("pol0")->Draw("same");
    g[iP]->Draw("pe5same");
    if (i_particle==0 && iP==1)
      text.DrawLatex(text_position_low_left[i_particle][0],text_position_low_left[i_particle][1],Form(format_fit_results_more[i_particle],h[iP]->GetFunction("pol0")->GetParameter(0),h[iP]->GetFunction("pol0")->GetParError(0)));
    else text.DrawLatex(text_position_low_left[i_particle][0],text_position_low_left[i_particle][1],Form(format_fit_results[i_particle],h[iP]->GetFunction("pol0")->GetParameter(0),h[iP]->GetFunction("pol0")->GetParError(0)));
    // TH1* syst = (TH1*)input.Get(Form("ratio/%i/syst",iP));
    //stat->Draw("esamex0");
    // syst->Draw("e2same");
    // l.DrawLine(0.7,1.,6.3,1.);
  }
  cv->SaveAs("RatioRun2_proton.pdf");

}
