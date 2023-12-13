// particle sequence: pions, protons, hypertriton ,he3

constexpr double pol_err[] = {0.003, 0.003, 0.002, 0., 0.};
constexpr float minpt[] = {0.66, 0.4, 0., 1.6, 0.5, 1.5};
constexpr float maxpt[] = {1.64, 3.1, 37., 8.4, 10.5, 3.1};
constexpr float miny[] = {0.94, 0.93, -0.4, -0.4, 0.5, -0.6};
constexpr float maxy[] = {1.06, 1.07, 2.4, 2.4, 1.5, 2.6};
const char *histoTitles[] = {";#it{p}_{T} (GeV/#it{c});#pi^{-}/#pi^{+};", ";#it{p}_{T} (GeV/#it{c});#bar{p}/p;", ";#it{ct} (cm);{}^{3}_{#bar{#Lambda}}#bar{H}/^{3}_{#Lambda}H;", ";#it{p}_{T} (GeV/#it{c});^{3}#bar{He}/^{3}He;", ";#it{ct} (cm);#bar{#Omega}^{+}/#Omega^{-};", ";#it{p}_{T} (GeV/#it{c});{}^{3}#bar{H}/^{3}H;"};
const char *inputFiles[] = {"Pion_PbPb/out/SystematicsAllEPtNotCombined_extend2.root", "Proton_PbPb/out/SystematicsAllEPtNotCombinedTOF_extend_4.root", "Hypertriton_PbPb/Ratio.root", "He3_PbPb/out/SpectraHe3_kINT7.root", "Omega_PbPb/ratio_cutCompetingMass-3.root", "Triton_PbPb/out/SpectraHe3.root"};
const char *outputFiles[] = {"pion","proton","hyp","he3","omega","triton"};
const char *outputFormat[] = {"ratio\tstat\tuncorr_sys\tcorr_sys","ratio\tstat\tuncorr_sys\tcorr_sys","ratio\tstat\tuncorr_sys\tcorr_sys","ratio\tstat\tuncorr_sys\tcorr_sys","ratio\tstat\tuncorr_sys","ratio\tstat\tuncorr_sys\tcorr_sys"};
const char *inputFilesSys[] = {"", "", "Hypertriton_PbPb/SystematicsRatio.root", "He3_PbPb/out/SystematicsAll_extend2.root", "Omega_PbPb/ratio_cutCompetingMass-3.root", "Triton_PbPb/out/SystematicsAll_extend.root"};
const char *inputHistoSys[] = {"fRatioDistributionTrials_%s", "fRatioDistributionTrials_%s", "fParameterDistribution_%s", "hist/fFitPar_%s", "h_sys;%d", "hist/fFitPar_%s"};
const char *inputFilesSysMC[] = {"", "", "", "He3_PbPb/out/SystematicsEfficiencyPrimary.root", "", ""};
const double xsec_err_sq[] = {0.0073795806*0.0073795806, 0.00129523*0.00129523+0.000173974*0.000173974, 0.00905074*0.00905074+0.00487767*0.00487767, 0.00905074*0.00905074+0.00487767*0.00487767, 0, 0.107459*0.107459+0.0177524*0.0177524};
const char *dir[] = {"", "", "", "1.0_89_0.1_2.5_1_1_1/", "", "1.0_69_0.1_2.5_1_1_1/"};
const int n_part = 6;
const bool draw_text = true;

std::array<TPad*,5> CreatePads(TCanvas* &cv, int i_part=0)
{
  if (!cv) cv = new TCanvas;
  std::array<TPad*,5> pads{nullptr};

  constexpr double sx[2]{1,1.-sx[0]};
  constexpr double sy[2]{0.23, 0.18};
  constexpr double fx = 0.9;
  constexpr double fy = 0.17;

  double global_y = 1.;
  for (int iP = 0; iP < 5; ++iP) {
    const int col = 0;
    const int row = iP;

    const int top = (row == 0);
    const int bot = row == 4;
    const int mid = !(top || bot);

    std::cout << sx[0] * col << "\t" << global_y - sy[mid] << "\t" << sx[0] + col * (1. - sx[0]) << "\t" << global_y << std::endl;
    pads[iP] = new TPad(Form("ratio%i",iP),"", 0, std::abs(global_y - sy[mid]), 1, global_y);
    global_y -= sy[mid];

    pads[iP]->SetRightMargin(0.01);
    pads[iP]->SetLeftMargin(0.13);
    pads[iP]->SetTopMargin(top * (1 - fy / sy[mid]));
    pads[iP]->SetBottomMargin(bot * (1 - fy / sy[mid]));

    cv->cd();
    pads[iP]->Draw();
    pads[iP]->cd();
    TH2F *rframe = new TH2F(Form("rframe%i",iP),histoTitles[i_part],100,minpt[i_part],maxpt[i_part],100,miny[i_part],maxy[i_part]);

    rframe->GetYaxis()->CenterTitle();
    rframe->GetYaxis()->SetTickLength(0.03 * (iP > 0 && iP < 4 ? .75 : 1.));
    rframe->GetYaxis()->SetTitleSize(26);
    rframe->GetYaxis()->SetTitleFont((!col) * 43);
    rframe->GetYaxis()->SetTitleOffset((i_part == 3 || i_part == 2 || i_part == 5) ? 1.08 : ( (i_part < 2) ? 1.18 : 1.15));
    rframe->GetYaxis()->SetNdivisions(505);
    rframe->GetYaxis()->SetDecimals(1);
    rframe->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    rframe->GetYaxis()->SetLabelSize((!col) * 17);

    rframe->GetXaxis()->SetTickLength(0.03 * sy[0] / sy[mid]);
    rframe->GetXaxis()->SetTitleSize(bot * 26);
    rframe->GetXaxis()->SetTitleFont(43);
    rframe->GetXaxis()->SetTitleOffset(0.85);
    rframe->GetXaxis()->SetNdivisions(510);
    rframe->GetXaxis()->SetDecimals(1);
    rframe->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    rframe->GetXaxis()->SetLabelSize(bot * 17);

    rframe->Draw("col");
  }
  return pads;
}

Color_t colors[]={kRed, kOrange-3, kAzure+4,kGreen+2,kMagenta+2};
double x_limits[][2]={{0.7,1.6},{0.5, 3.},{2., 35.},{2., 8.},{1.,10.},{1.6, 3.}};
double x_limits_30_50[][2]={{0.7,1.6},{.5, 3.},{2., 14.},{2., 7.},{1.,10.},{1.6, 3.}};
double material_ratio_instance[] = {5,4,3,1,-10,2};
const char *format_fit_results[] = {"#it{R} = %.3f #pm %.3f (uncorr.) #pm %.3f (corr.)","#it{R} = %.3f #pm %.3f (uncorr.) #pm %.3f (corr.)","#it{R} = %.2f #pm %.2f (uncorr.) #pm %.2f (corr.)","#it{R} = %.2f #pm %.2f (uncorr.) #pm %.2f (corr.)","#it{R} = %.3f #pm %.3f (uncorr.)","#it{R} = %.2f #pm %.2f (uncorr.) #pm %.2f (corr.)"};
const char *format_out_results[] = {"%s\t%.12f\t0.000000000000\t%.12f\t%.12f","%s\t%.12f\t0.000000000000\t%.12f\t%.12f","%s\t%.12f\t%.12f\t%.12f\t%.12f","%s\t%.12f\t%.12f\t%.12f\t%.12f","%s\t%.12f\t%.12f\t%.12f","%s\t%.12f\t%.12f\t%.12f\t%.12f"};

void newratios() {
  std::ofstream out_ratios_file;
  out_ratios_file.open("run2_ratios.dat");
  gStyle->SetOptStat(0);
  const float leftmargin = 0.14;
  const float rigthmargin = 0.04;
  const float deltaRatio = 0.89;
  const float kMarkerSize = 0.9;

  const float sizeRatioPadY = 0.15;
  const float sizeLastPad = (1. - 3 * sizeRatioPadY) * 0.5f;
  const float titleYoffset = 2.3;

  for (int i_part =0; i_part<n_part; ++i_part){
    out_ratios_file << outputFiles[i_part] << std::endl;
    out_ratios_file << outputFormat[i_part] << std::endl;
    TCanvas *cv = new TCanvas("cv","cv",480,900);
    auto pads = CreatePads(cv, i_part);

    TLatex text;
    text.SetTextFont(43);
    text.SetTextSize(24);
    pads[0]->cd();
    double mean_y = 0.5*(miny[i_part]+maxy[i_part]);
    double half_width_y = 0.5*(maxy[i_part]-miny[i_part]);
    double mean_x = 0.5*(minpt[i_part]+maxpt[i_part]);
    double half_width_x = 0.5*(maxpt[i_part]-minpt[i_part]);
    text.DrawText(mean_x+0.6*half_width_x,mean_y+0.68*half_width_y,"ALICE");

    pads[0]->cd();
    text.SetTextFont(43);
    text.SetTextSize(24);
    text.DrawLatex(mean_x-0.22*half_width_x,mean_y+0.39*half_width_y,"Pb#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV");

    const string labels[5]{"0-5%","5-10%","30-50%","10-30%","50-90%"};
    const string names[5]{"0_5","5_10","30_50","10_30","50_90"};
    TFile input(inputFiles[i_part]);
    TFile input_material("FinalPlot3D.root");
    TFile input_sys(inputFilesSys[i_part]);
    TFile input_sys_MC(inputFilesSysMC[i_part]);
    TGraphErrors *g[5],*gSys[5];
    TH1D* h[5];
    TF1 *f[5];
    for (int iPP = 0; iPP < 5; ++iPP) {
      int iP = iPP;
      if (i_part==2&&iPP==4)continue;
      if (iPP==2)iP=3;
      if (iPP==3)iP=2;
      pads[iPP]->cd();
      text.SetTextSize(24);
      text.DrawText(mean_x-0.93*half_width_x,mean_y+0.68*half_width_y,labels[iP].data());
      if (i_part!=4) h[iP] = (TH1D*)input.Get(Form("%sfRatio_%s",dir[i_part],names[iP].data()));
      else h[iP] = (TH1D*)input.Get(Form("%sh_ratio_%s",dir[i_part],names[iP].data()));
      g[iP]=new TGraphErrors(h[iP]);
      g[iP]->SetMarkerStyle(20);
      g[iP]->SetMarkerSize(1.);
      g[iP]->SetMarkerColor(colors[iP]);
      g[iP]->SetLineColor(colors[iP]);
      g[iP]->SetFillStyle(0);
      double ratio_mat[3];
      double material_error=0;
      if (i_part!=4){
        for (int iMat=0;iMat<3;++iMat){
          auto hMaterial=(TH1D*)input_material.Get(Form("fRatio_%s;%.0f",names[iP].data(),iMat*(5-(int)(iP>3))+(material_ratio_instance[i_part])-((int)(iP>3 && i_part<2))));
          ratio_mat[iMat] = hMaterial->GetFunction("pol0")->GetParameter(0);
        }
        material_error = std::abs(ratio_mat[1]-ratio_mat[0])*0.5;
      }

      int iPoint{0};
      while (iPoint<g[iP]->GetN()){
        double x=g[iP]->GetPointX(iPoint);
        if ( ((x<x_limits[i_part][0]||x>x_limits[i_part][1])&&((iP!=2&&i_part!=3)||((iP!=2&&iP!=3&&i_part==3)))) || ((x<x_limits_30_50[i_part][0]||x>x_limits_30_50[i_part][1])&&((iP==2&&i_part!=3)||((iP>=2&&iP<=3&&i_part==3))||((iP>=2&&iP<=3&&i_part==2)))) ){
          g[iP]->RemovePoint(iPoint);
        }
        else ++iPoint;
      }
      g[iP]->SetLineColor(colors[iP]);
      g[iP]->SetMarkerColor(colors[iP]);

      double sys_tot = 0.;
      if (i_part>1){
        TH1D *sys=nullptr;
        if (i_part!=4){
          sys=(TH1D*)input_sys.Get(Form(inputHistoSys[i_part],names[iP].data()));
        }
        else sys=(TH1D*)input_sys.Get(Form(inputHistoSys[i_part],(iP*2)+1));
        double sysMC_val = 0;
        if (i_part==3) {
          TH1D *sys_MC=(TH1D*)input_sys_MC.Get(Form("fRatioDistribution_%s",names[iP].data()));
          sysMC_val=sys_MC->GetRMS();
        }
        gSys[iP]=new TGraphErrors(*g[iP]);
        gSys[iP]->SetLineColor(colors[iP]);
        gSys[iP]->SetMarkerColor(colors[iP]);
        sys_tot = TMath::Sqrt(sys->GetRMS()*sys->GetRMS()+sysMC_val*sysMC_val);
        for (int iPoint=0;iPoint<g[iP]->GetN();++iPoint){
          gSys[iP]->SetPointError(iPoint,h[iP]->GetBinWidth(h[iP]->FindBin(g[iP]->GetPointX(iPoint)))*0.5,sys_tot);
          g[iP]->SetPointError(iPoint,0.,g[iP]->GetErrorY(iPoint));


          std::cout<<"i_part = "<<i_part<<std::endl;
        }
        //std::cout<<"i_part = "<<i_part<<std::endl;
      }
      f[iP]=new TF1("fit","pol0",minpt[i_part],maxpt[i_part]);
      f[iP]->SetParameter(0,h[iP]->GetFunction("pol0")->GetParameter(0));
      f[iP]->SetLineColor(kBlack);
      f[iP]->Draw("same");
      i_part<2 ? g[iP]->Draw("pe5same") : g[iP]->Draw("pesame");
      if (i_part>1) gSys[iP]->Draw("pe5same");
      text.SetTextSize(18);
      double sys_err = sqrt((double)(i_part<2)*h[iP]->GetFunction("pol0")->GetParError(0)*h[iP]->GetFunction("pol0")->GetParError(0)+(double)(i_part>1)*sys_tot*sys_tot);
      if (i_part<2){
        TH1D* sys=(TH1D*)input.Get(Form(inputHistoSys[i_part],names[iP].data()));
        sys_tot = TMath::Sqrt(sys->GetRMS()*sys->GetRMS());
      }
      //else sys_tot = 0.;
      material_error = sqrt(material_error*material_error+h[iP]->GetFunction("pol0")->GetParameter(0)*h[iP]->GetFunction("pol0")->GetParameter(0)*(xsec_err_sq[i_part])+(double)(i_part<2)*sys_tot*sys_tot);
      if (draw_text) {
        if (i_part<2){
          text.DrawLatex(mean_x-0.91*half_width_x,mean_y-0.86*half_width_y,Form(format_fit_results[i_part],h[iP]->GetFunction("pol0")->GetParameter(0),sys_err,sqrt(material_error*material_error + pol_err[iPP]*pol_err[iPP])));
          out_ratios_file <<Form(format_out_results[i_part],labels[iP].data(),h[iP]->GetFunction("pol0")->GetParameter(0),sys_err,sqrt(material_error*material_error + pol_err[iPP]*pol_err[iPP]))<< std::endl;
        }
        else if (i_part==4){
          text.DrawLatex(mean_x-0.91*half_width_x,mean_y-0.86*half_width_y,Form(format_fit_results[i_part],h[iP]->GetFunction("pol0")->GetParameter(0),std::hypot(h[iP]->GetFunction("pol0")->GetParError(0),sys_tot)));
          out_ratios_file <<Form(format_out_results[i_part],labels[iP].data(),h[iP]->GetFunction("pol0")->GetParameter(0),std::hypot(h[iP]->GetFunction("pol0")->GetParError(0),sys_tot))<< std::endl;
        }
        else {
          text.DrawLatex(mean_x-0.91*half_width_x,mean_y-0.86*half_width_y,Form(format_fit_results[i_part],h[iP]->GetFunction("pol0")->GetParameter(0),std::hypot(h[iP]->GetFunction("pol0")->GetParError(0),sys_err),material_error));
          out_ratios_file <<Form(format_out_results[i_part],labels[iP].data(),h[iP]->GetFunction("pol0")->GetParameter(0),std::hypot(h[iP]->GetFunction("pol0")->GetParError(0),sys_err),material_error)<< std::endl;
        }
        text.DrawLatex(mean_x-0.5*half_width_x,mean_y+0.7*half_width_y,Form("#chi^{2}/NDF = %.2f/%d",h[iP]->GetFunction("pol0")->GetChisquare(),h[iP]->GetFunction("pol0")->GetNDF()));
      }
    }
    cv->SaveAs(Form("RatioRun2_%s.pdf",outputFiles[i_part]));
    out_ratios_file << "\n" << std::endl;
  }
  out_ratios_file.close();
}
