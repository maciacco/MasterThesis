// particle sequence: pions, protons, hypertriton ,he3

constexpr float minpt[] = {0.66, 0.4, 0., 1.6};
constexpr float maxpt[] = {1.64, 3.1, 37., 8.4};
constexpr float miny[] = {0.965, 0.94, 0.3, 0.55};
constexpr float maxy[] = {1.035, 1.06, 1.7, 1.45};
const char *histoTitles[] = {";#it{p}_{T} (GeV/#it{c});#pi^{-}/#pi^{+};", ";#it{p}_{T} (GeV/#it{c});#bar{p}/p;", ";#it{ct} (cm);{}^{3}_{#bar{#Lambda}}#bar{H}/^{3}_{#Lambda}H;", ";#it{p}_{T} (GeV/#it{c});^{3}#bar{He}/^{3}He;"};
const char *inputFiles[] = {"Pion_PbPb/out/SystematicsAllEPtNotCombined.root", "Proton_PbPb/out/SystematicsAllEPtNotCombinedTOFTEST.root", "Hypertriton_PbPb/Ratio.root", "He3_PbPb/out/SpectraHe3.root"};
const char *outputFiles[] = {"pion","proton","hyp","he3"};
const char *inputFilesSys[] = {"", "", "Hypertriton_PbPb/Systematics.root", "He3_PbPb/out/SystematicsAll.root"};
const char *inputHistoSys[] = {"", "", "fParameterDistribution_%s", "hist/fFitPar_%s"};
const char *inputFilesSysMC[] = {"", "", "", "He3_PbPb/out/SystematicsEfficiencyPrimary.root"};
const double xsec_err_sq[] = {0., 0.00115085*0.00115085+0.000172469*0.000172469, 0.00861352*0.00861352+0.00473124*0.00473124, 0.00861352*0.00861352+0.00473124*0.00473124};
const char *dir[] = {"", "", "", "1.0_89_0.1_2.5_1_1_1/"};
const int n_part = 4;
const bool draw_text = true;

std::array<TPad*,3> CreatePads(TCanvas* &cv, int i_part=0)
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
    TH2F *rframe = new TH2F(Form("rframe%i",iP),histoTitles[i_part],100,minpt[i_part],maxpt[i_part],100,miny[i_part],maxy[i_part]);

    rframe->GetYaxis()->CenterTitle();
    rframe->GetYaxis()->SetTickLength(0.03 * sx[0] / sx[col]);
    rframe->GetYaxis()->SetTitleSize(20);
    rframe->GetYaxis()->SetTitleFont((!col) * 43);
    rframe->GetYaxis()->SetTitleOffset(2.);
    rframe->GetYaxis()->SetNdivisions(505);
    rframe->GetYaxis()->SetDecimals(1);
    rframe->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    rframe->GetYaxis()->SetLabelSize((!col) * 15);

    rframe->GetXaxis()->SetTickLength(0.03 * sy[0] / sy[mid]);
    rframe->GetXaxis()->SetTitleSize(bot * 20);
    rframe->GetXaxis()->SetTitleFont(43);
    rframe->GetXaxis()->SetTitleOffset(3);
    rframe->GetXaxis()->SetNdivisions(510);
    rframe->GetXaxis()->SetDecimals(1);
    rframe->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    rframe->GetXaxis()->SetLabelSize(bot * 15);

    rframe->Draw("col");
  }
  return pads;
}

Color_t colors[]={kOrange+7, kAzure+4, kTeal+4};
double x_limits[][2]={{0.7,1.6},{0.5, 3.},{2., 35.},{2., 8.}};
double x_limits_30_50[][2]={{0.7,1.6},{.5, 3.},{2., 14.},{2., 7.}};
double material_ratio_instance[] = {4,3,2,1};
const char *format_fit_results[] = {"R=%.3f #pm0.000 (stat.) #pm%.3f (syst.) #pm%.3f (mat.)","R=%.3f #pm0.000 (stat.) #pm%.3f (syst.) #pm%.3f (mat.)","R=%.2f #pm%.2f (stat.) #pm%.2f (syst.) #pm%.2f (mat.)","R=%.2f #pm%.2f (stat.) #pm%.2f (syst.) #pm%.2f (mat.)"};

void newratios() {
  gStyle->SetOptStat(0);
  const float leftmargin = 0.16;
  const float rigthmargin = 0.04;
  const float deltaRatio = 0.89;
  const float kMarkerSize = 0.9;

  const float sizeRatioPadY = 0.15;
  const float sizeLastPad = (1. - 3 * sizeRatioPadY) * 0.5f;
  const float titleYoffset = 2.3;

  for (int i_part =0; i_part<n_part; ++i_part){
    TCanvas *cv = new TCanvas("cv","cv",480,800);
    auto pads = CreatePads(cv, i_part);

    TLatex text;
    text.SetTextFont(63);
    text.SetTextSize(18);
    pads[0]->cd();
    double mean_y = 0.5*(miny[i_part]+maxy[i_part]);
    double half_width_y = 0.5*(maxy[i_part]-miny[i_part]);
    double mean_x = 0.5*(minpt[i_part]+maxpt[i_part]);
    double half_width_x = 0.5*(maxpt[i_part]-minpt[i_part]);
    text.DrawText(mean_x+0.08*half_width_x,mean_y+0.75*half_width_y,"ALICE Preliminary");

    pads[0]->cd();
    text.SetTextFont(43);
    text.DrawLatex(mean_x+0.08*half_width_x,mean_y+0.55*half_width_y,"Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");

    const string labels[3]{"0-5%","5-10%","30-50%"};
    const string names[3]{"0_5","5_10","30_50"};
    TFile input(inputFiles[i_part]);
    TFile input_material("FinalPlot3D.root");
    TFile input_sys(inputFilesSys[i_part]);
    TFile input_sys_MC(inputFilesSysMC[i_part]);
    TGraphErrors *g[3],*gSys[3];
    TH1D* h[3];
    TF1 *f[3];
    for (int iP = 0; iP < 3; ++iP) {
      pads[iP]->cd();
      text.SetTextSize(18);
      text.DrawText(mean_x-0.89*half_width_x,mean_y+0.75*half_width_y,labels[iP].data());
      h[iP] = (TH1D*)input.Get(Form("%sfRatio_%s",dir[i_part],names[iP].data()));
      g[iP]=new TGraphErrors(h[iP]);
      g[iP]->SetMarkerStyle(20);
      g[iP]->SetMarkerSize(1.);
      g[iP]->SetMarkerColor(colors[iP]);
      g[iP]->SetLineColor(colors[iP]);
      g[iP]->SetFillStyle(0);
      double ratio_mat[3];
      for (int iMat=0;iMat<3;++iMat){
        auto hMaterial=(TH1D*)input_material.Get(Form("fRatio_%s;%.0f",names[iP].data(),iMat*4+material_ratio_instance[i_part]));
        ratio_mat[iMat] = hMaterial->GetFunction("pol0")->GetParameter(0);
      }
      double material_error = std::abs(ratio_mat[1]-ratio_mat[0])*0.5;
      int iPoint{0};
      while (iPoint<g[iP]->GetN()){
        double x=g[iP]->GetPointX(iPoint);
        if ( ((x<x_limits[i_part][0]||x>x_limits[i_part][1])&&iP<2) || ((x<x_limits_30_50[i_part][0]||x>x_limits_30_50[i_part][1])&&iP==2) ){
          g[iP]->RemovePoint(iPoint);
        }
        else ++iPoint;
      }
      double sys_tot = 0.;
      if (i_part>1){
      TH1D *sys=(TH1D*)input_sys.Get(Form(inputHistoSys[i_part],names[iP].data()));
        double sysMC_val = 0;
        if (i_part>2) {
          TH1D *sys_MC=(TH1D*)input_sys_MC.Get(Form("fRatioDistribution_%s",names[iP].data()));
          sysMC_val=sys_MC->GetRMS();
        }
        gSys[iP]=new TGraphErrors(*g[iP]);
        sys_tot = TMath::Sqrt(sys->GetRMS()*sys->GetRMS()+sysMC_val*sysMC_val);
        for (int iPoint=0;iPoint<g[iP]->GetN();++iPoint){
          gSys[iP]->SetPointError(iPoint,h[iP]->GetBinWidth(h[iP]->FindBin(g[iP]->GetPointX(iPoint)))*0.5,sys_tot);
          g[iP]->SetPointError(iPoint,0.,g[iP]->GetErrorY(iPoint));
        }
      }
      f[iP]=new TF1("fit","pol0",minpt[i_part],maxpt[i_part]);
      f[iP]->SetParameter(0,h[iP]->GetFunction("pol0")->GetParameter(0));
      f[iP]->SetLineColor(kBlack);
      f[iP]->Draw("same");
      i_part<2 ? g[iP]->Draw("pe5same") : g[iP]->Draw("pesame");
      if (i_part>1) gSys[iP]->Draw("pe5same");
      text.SetTextSize(15);
      double sys_err = sqrt(h[iP]->GetFunction("pol0")->GetParameter(0)*h[iP]->GetFunction("pol0")->GetParameter(0)*(xsec_err_sq[i_part])+(1-i_part/2)*h[iP]->GetFunction("pol0")->GetParError(0)*h[iP]->GetFunction("pol0")->GetParError(0)+(i_part/2)*sys_tot*sys_tot);
      if (draw_text) {
        if (i_part<2)text.DrawLatex(mean_x-0.91*half_width_x,mean_y-0.8*half_width_y,Form(format_fit_results[i_part],h[iP]->GetFunction("pol0")->GetParameter(0),sys_err,material_error));
        else text.DrawLatex(mean_x-0.91*half_width_x,mean_y-0.8*half_width_y,Form(format_fit_results[i_part],h[iP]->GetFunction("pol0")->GetParameter(0),h[iP]->GetFunction("pol0")->GetParError(0),sys_err,material_error));
      }
    }
    cv->SaveAs(Form("RatioRun2_%s.eps",outputFiles[i_part]));
  }
}
