// particle sequence: pions, protons, hypertriton ,he3

constexpr double pol_err[] = {0.003, 0.003, 0.002, 0., 0.};
constexpr float minpt[] = {0.66, 0.4, 0., 1.6, 0.5, 1.5};
constexpr float maxpt[] = {1.64, 3.1, 37., 8.4, 10.5, 3.1};
constexpr float miny[] = {0.94, 0.93, -0.4, -0.4, 0.5, -0.6};
constexpr float maxy[] = {1.06, 1.07, 2.4, 2.4, 1.5, 2.6};
const char *histoTitles[] = {";#it{p}_{T} (GeV/#it{c});#pi^{-}/#pi^{+};", ";#it{p}_{T} (GeV/#it{c});#bar{p}/p;", ";#it{ct} (cm);{}^{3}_{#bar{#Lambda}}#bar{H}/^{3}_{#Lambda}H;", ";#it{p}_{T} (GeV/#it{c});^{3}#bar{He}/^{3}He;", ";#it{ct} (cm);#bar{#Omega}^{+}/#Omega^{-};", ";#it{p}_{T} (GeV/#it{c});{}^{3}#bar{H}/^{3}H;"};
const char *inputFiles[] = {"Pion_PbPb/out/SystematicsAllEPtNotCombined_extend2.root", "Proton_PbPb/out/SystematicsAllEPtNotCombinedTOF_extend_4.root", "Hypertriton_PbPb/Ratio.root", "He3_PbPb/out/SpectraHe3_kINT7.root", "Omega_PbPb/ratio_cutCompetingMass-3.root", "Triton_PbPb/out/SpectraHe3.root"};
const char *outputFiles[] = {"pion","proton","hyp","he3","omega","triton"};
const char *outputFiles2[] = {"$\\pi^+$","$\\mathrm{p}$","^{3}_{\\Lambda}\\mathrm{H}","^{3}\\mathrm{He}","\\Omega^-","^{3}\\mathrm{H}"};
const char *outputFormat[] = {"ratio\tstat\tuncorr_sys\tcorr_sys","ratio\tstat\tuncorr_sys\tcorr_sys","ratio\tstat\tuncorr_sys\tcorr_sys","ratio\tstat\tuncorr_sys\tcorr_sys","ratio\tstat\tuncorr_sys","ratio\tstat\tuncorr_sys\tcorr_sys"};
const char *inputFilesSys[] = {"", "", "Hypertriton_PbPb/SystematicsRatio.root", "He3_PbPb/out/SystematicsAll_extend2.root", "Omega_PbPb/ratio_cutCompetingMass-3.root", "Triton_PbPb/out/SystematicsAll_extend.root"};
const char *inputHistoSys[] = {"fRatioDistributionTrials_%s", "fRatioDistributionTrials_%s", "fParameterDistribution_%s", "hist/fFitPar_%s", "h_sys;%d", "hist/fFitPar_%s"};
const char *inputFilesSysMC[] = {"", "", "", "He3_PbPb/out/SystematicsEfficiencyPrimary.root", "", ""};
const double xsec_err_sq[] = {0.0073795806*0.0073795806, 0.00129523*0.00129523+0.000173974*0.000173974, 0.00905074*0.00905074+0.00487767*0.00487767, 0.00905074*0.00905074+0.00487767*0.00487767, 0, 0.107459*0.107459+0.0177524*0.0177524};
const char *dir[] = {"", "", "", "1.0_89_0.1_2.5_1_1_1/", "", "1.0_69_0.1_2.5_1_1_1/"};
const int n_part = 6;

const char *particle_ratios[] = {"\\pi^{-} / \\pi^{+}","\\bar{\\mathrm{p}} / \\mathrm{p}","{}_{\\overline{\\Lambda}}^{3}\\overline{\\mathrm{H}} / ^{3}_{\\Lambda}\\mathrm{H}","^{3}\\overline{\\mathrm{He}} / ^{3}\\mathrm{He}", "\\overline{\\Omega}^{+} / \\Omega^{-}","^{3}\\overline{\\mathrm{H}} / ^{3}\\mathrm{H}"};

double x_limits[][2]={{0.7,1.6},{0.5, 3.},{2., 35.},{2., 8.},{1.,10.},{1.6, 3.}};
double x_limits_30_50[][2]={{0.7,1.6},{.5, 3.},{2., 14.},{2., 7.},{1.,10.},{1.6, 3.}};
double material_ratio_instance[] = {5,4,3,1,-10,2};
const char *format_fit_results[] = {"#it{R} = %.3f #pm %.3f (uncorr.) #pm %.3f (corr.)","#it{R} = %.3f #pm %.3f (uncorr.) #pm %.3f (corr.)","#it{R} = %.2f #pm %.2f (uncorr.) #pm %.2f (corr.)","#it{R} = %.2f #pm %.2f (uncorr.) #pm %.2f (corr.)","#it{R} = %.3f #pm %.3f (uncorr.)","#it{R} = %.2f #pm %.2f (uncorr.) #pm %.2f (corr.)"};
const char *format_out_results[] = {"%s\t%.12f\t0.000000000000\t%.12f\t%.12f","%s\t%.12f\t0.000000000000\t%.12f\t%.12f","%s\t%.12f\t%.12f\t%.12f\t%.12f","%s\t%.12f\t%.12f\t%.12f\t%.12f","%s\t%.12f\t%.12f\t%.12f","%s\t%.12f\t%.12f\t%.12f\t%.12f"};

void hepdata_ratios() {
  std::ofstream ofile_avg[5];
  for (int ii{0}; ii < 5; ++ii){
    int i = ii;
    if (ii == 2) i = 3;
    else if (ii == 3) i = 2;
    ofile_avg[i].open(Form("Ratios_%d.yaml", ii));
    ofile_avg[i] << "independent_variables:" << "\n";
    ofile_avg[i] << "- header: {name: Species}" << "\n";
    ofile_avg[i] << "  values:" << "\n";
    for (int iP{0}; iP < 6; ++iP){
      if (i == 4 && iP == 2) continue;
      ofile_avg[i] << "  - {value: \'$" << outputFiles2[iP] << "$\'}" << "\n";
    }
    ofile_avg[i] << "dependent_variables:" << "\n";
    ofile_avg[i] << "- header: {name: Antimatter / matter}" << "\n";
    ofile_avg[i] << "  values:" << "\n";
  }

  for (int i_part =0; i_part<n_part; ++i_part){
    TLatex text;

    const string labels[5]{"0-5","5-10","30-50","10-30","50-90"};
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

      std::ofstream ofile(Form("Ratios_%s_%d.yaml", outputFiles[i_part], iPP));
      ofile << "dependent_variables:" << "\n";
      ofile << "- header: {name: '$" << particle_ratios[i_part] << "$'}" << "\n";

      if (i_part!=4) h[iP] = (TH1D*)input.Get(Form("%sfRatio_%s",dir[i_part],names[iP].data()));
      else h[iP] = (TH1D*)input.Get(Form("%sh_ratio_%s",dir[i_part],names[iP].data()));
      g[iP]=new TGraphErrors(h[iP]);

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
        sys_tot = TMath::Sqrt(sys->GetRMS()*sys->GetRMS()+sysMC_val*sysMC_val);
        // if (i_part == 3) std::cout << sys_tot << std::endl;
        for (int iPoint=0;iPoint<g[iP]->GetN();++iPoint){
          gSys[iP]->SetPointError(iPoint,h[iP]->GetBinWidth(h[iP]->FindBin(g[iP]->GetPointX(iPoint)))*0.5,sys_tot);
          g[iP]->SetPointError(iPoint,0.,g[iP]->GetErrorY(iPoint));
        }
        //std::cout<<"i_part = "<<i_part<<std::endl;
      }
      double sys_err = sqrt((double)(i_part<2)*h[iP]->GetFunction("pol0")->GetParError(0)*h[iP]->GetFunction("pol0")->GetParError(0)+(double)(i_part>1)*sys_tot*sys_tot);
      if (i_part<2){
        TH1D* sys=(TH1D*)input.Get(Form(inputHistoSys[i_part],names[iP].data()));
        sys_tot = TMath::Sqrt(sys->GetRMS()*sys->GetRMS());
      }
      //else sys_tot = 0.;
      material_error = sqrt(material_error*material_error+h[iP]->GetFunction("pol0")->GetParameter(0)*h[iP]->GetFunction("pol0")->GetParameter(0)*(xsec_err_sq[i_part])+(double)(i_part<2)*sys_tot*sys_tot);
      if (true) {
        if (i_part<2){
          ofile << "  values:" << "\n";
          for (int i{0}; i < g[iP]->GetN(); ++i){
            ofile << "  - errors:" << "\n";
            ofile << "    - {label: Uncorrelated uncertainty, symerror: " << setprecision(6) << g[iP]->GetErrorY(i) << "}" << "\n";
            ofile << "    - {label: Total correlated systematic, symerror: " << setprecision(6) << sqrt(material_error*material_error + pol_err[iPP]*pol_err[iPP]) << "}" << "\n";
            ofile << "    value: " << setprecision(6) << g[iP]->GetPointY(i) << "\n";
          }
          ofile << "independent_variables:" << "\n";
          ofile << "- header: {name: '$p_{\\mathrm{T}}$ (GeV/$c$)'}" << "\n";
          ofile << "  values:" << "\n";
          for (int i{0}; i < g[iP]->GetN(); ++i){
            ofile << "  - {high: " << setprecision(6) << g[iP]->GetPointX(i) + g[iP]->GetErrorX(i) << ", low: " << g[iP]->GetPointX(i) - g[iP]->GetErrorX(i) << "}" << "\n";
          }

          ofile_avg[iP] << "  - errors:" << "\n";
          ofile_avg[iP] << "    - {label: Total uncorrelated uncertainty, symerror: " << setprecision(6) << sys_err << "}\n";
          ofile_avg[iP] << "    - {label: Total correlated uncertainty, symerror: " << setprecision(6) << sqrt(material_error*material_error + pol_err[iPP]*pol_err[iPP]) << "}\n";
          ofile_avg[iP] << "    value: " << setprecision(6) << h[iP]->GetFunction("pol0")->GetParameter(0) << "\n";
        }
        else if (i_part==4){
          ofile << "  values:" << "\n";
          for (int i{0}; i < g[iP]->GetN() - 1 * (iPP == 4); ++i){
            ofile << "  - errors:" << "\n";
            ofile << "    - {label: Statistical uncertainty, symerror: " << setprecision(6) << g[iP]->GetErrorY(i) << "}" << "\n";
            ofile << "    - {label: Total uncorrelated systematic, symerror: " << setprecision(6) << gSys[iP]->GetErrorY(i) << "}" << "\n";
            ofile << "    value: " << setprecision(6) << g[iP]->GetPointY(i) << "\n";
          }
          ofile << "independent_variables:" << "\n";
          ofile << "- header: {name: '$c\\mathrm{t}$ (cm)'}" << "\n";
          ofile << "  values:" << "\n";
          for (int i{0}; i < g[iP]->GetN() - 1 * (iPP == 4); ++i){
            ofile << "  - {high: " << setprecision(6) << gSys[iP]->GetPointX(i) + gSys[iP]->GetErrorX(i) << ", low: " << gSys[iP]->GetPointX(i) - gSys[iP]->GetErrorX(i) << "}" << "\n";
          }

          ofile_avg[iP] << "  - errors:" << "\n";
          ofile_avg[iP] << "    - {label: Total uncorrelated uncertainty, symerror: " << setprecision(6) << std::hypot(h[iP]->GetFunction("pol0")->GetParError(0),sys_tot) << "}\n";
          ofile_avg[iP] << "    value: " << setprecision(6) << h[iP]->GetFunction("pol0")->GetParameter(0) << "\n";
          //ofile <<Form(format_out_results[i_part],labels[iP].data(),h[iP]->GetFunction("pol0")->GetParameter(0),std::hypot(h[iP]->GetFunction("pol0")->GetParError(0),sys_tot))<< std::endl;
        }
        else {
          ofile << "  values:" << "\n";
          for (int i{0}; i < g[iP]->GetN(); ++i){
            ofile << "  - errors:" << "\n";
            ofile << "    - {label: Statistical uncertainty, symerror: " << setprecision(6) << g[iP]->GetErrorY(i) << "}" << "\n";
            ofile << "    - {label: Total uncorrelated systematic, symerror: " << setprecision(6) << gSys[iP]->GetErrorY(i) << "}" << "\n";
            ofile << "    - {label: Total correlated systematic, symerror: " << setprecision(6) << material_error << "}" << "\n";
            ofile << "    value: " << setprecision(6) << g[iP]->GetPointY(i) << "\n";
          }
          ofile << "independent_variables:" << "\n";
          ofile << "- header: {name: " << ((i_part == 2) ? "'$c\\mathrm{t}$ (cm)'" : "'$p_{\\mathrm{T}}$ (GeV/$c$)'") << "}" << "\n";
          ofile << "  values:" << "\n";
          for (int i{0}; i < g[iP]->GetN(); ++i){
            ofile << "  - {high: " << setprecision(6) << gSys[iP]->GetPointX(i) + gSys[iP]->GetErrorX(i) << ", low: " << gSys[iP]->GetPointX(i) - gSys[iP]->GetErrorX(i) << "}" << "\n";
          }

          ofile_avg[iP] << "  - errors:" << "\n";
          ofile_avg[iP] << "    - {label: Total uncorrelated uncertainty, symerror: " << setprecision(6) << std::hypot(h[iP]->GetFunction("pol0")->GetParError(0),sys_err) << "}\n";
          ofile_avg[iP] << "    - {label: Total correlated uncertainty, symerror: " << setprecision(6) << material_error << "}\n";
          ofile_avg[iP] << "    value: " << setprecision(6) << h[iP]->GetFunction("pol0")->GetParameter(0) << "\n";
          //ofile <<Form(format_out_results[i_part],labels[iP].data(),h[iP]->GetFunction("pol0")->GetParameter(0),std::hypot(h[iP]->GetFunction("pol0")->GetParError(0),sys_err),material_error)<< std::endl;
        }
      }
    }
  }
}
