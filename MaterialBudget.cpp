const char *species[]={"Pion","Proton","He3"};
const char *antiMatter[]={"A","M"};
const char *var[]={"Minus","Plus"};
const double cent[3][2]={{0,5},{5,10},{30,50}};

const bool VERBOSE=false;

void MaterialBudget(){
  TFile outFile("MaterialBudgetUncertainty.root","recreate");
  TFile fitSHM3D("FinalPlot3D.root");
  TFile *f[3][3];
  TH1D *effRatio[3][2][2];
  for (int iSp=0; iSp<3; ++iSp){
    for (int iC=0;iC<3;++iC){
      for (int iM=0; iM<2; ++iM){
        TH1D *eff[3];
        for (int iProd=0; iProd<3; ++iProd){
          f[iSp][iProd]=new TFile(Form("%s_PbPb/out/Efficiency%s_LHC22b9_%d.root",species[iSp],species[iSp],iProd+1));
          if (!f[iSp][iProd]) continue;
          TString detector("TOF");
          TString directory("_/");
          if (iSp>1) {
            detector=Form("TPC");
            directory=Form("");
          }
          if (VERBOSE) {
            std::cout<<"File Name = "<<Form("%s_PbPb/out/Efficiency%s_LHC22b9_%d.root",species[iSp],species[iSp],iProd+1)<<std::endl;
            std::cout<<"Efficiency name = "<<Form("%sf%sEff_%s_%.0f_%.0f",directory.Data(),antiMatter[iM],detector.Data(),cent[iC][0],cent[iC][1])<<std::endl;
          }
          eff[iProd]=(TH1D*)f[iSp][iProd]->Get(Form("%sf%sEff_%s_%.0f_%.0f",directory.Data(),antiMatter[iM],detector.Data(),cent[iC][0],cent[iC][1]));
          if (!eff[iProd]) continue;
          //outFile.cd();
          //eff[iProd]->Write();
        }
        for (int iVar=0; iVar<2; ++iVar){
          effRatio[iSp][iM][iVar]=new TH1D(*eff[0]);
          effRatio[iSp][iM][iVar]->SetName(Form("f%s%sEffRatio%s_%.0f_%.0f",antiMatter[iM],species[iSp],var[iVar],cent[iC][0],cent[iC][1]));
          iVar == 0 ? effRatio[iSp][iM][iVar]->Divide(eff[0],eff[1]) : effRatio[iSp][iM][iVar]->Divide(eff[2],eff[1]);
          effRatio[iSp][iM][iVar]->Fit("pol0");
          outFile.cd();
          effRatio[iSp][iM][iVar]->Write();
        }
      }
    }
  }
  for (int iC=0;iC<3;++iC){
    for (int iVar=0;iVar<2;++iVar){
      TH2D *fRatios=(TH2D*)fitSHM3D.Get(Form("fRatio_vs_b_%.0f_%.0f",cent[iC][0],cent[iC][1]));
      fRatios->SetBinContent(10, 6, fRatios->GetBinContent(10,6)*effRatio[2][1][iVar]->GetFunction("pol0")->GetParameter(0)/effRatio[2][0][iVar]->GetFunction("pol0")->GetParameter(0));
      fRatios->SetBinContent(9, 1, fRatios->GetBinContent(9,1)*effRatio[2][1][iVar]->GetFunction("pol0")->GetParameter(0)/effRatio[2][0][iVar]->GetFunction("pol0")->GetParameter(0));
      fRatios->SetBinContent(4, 6, fRatios->GetBinContent(4,6)*effRatio[1][1][iVar]->GetFunction("pol0")->GetParameter(0)/effRatio[1][0][iVar]->GetFunction("pol0")->GetParameter(0));
      fRatios->SetBinContent(1, 11, fRatios->GetBinContent(1,11)*effRatio[0][1][iVar]->GetFunction("pol0")->GetParameter(0)/effRatio[0][0][iVar]->GetFunction("pol0")->GetParameter(0));
      fRatios->SetBinError(10, 6, fRatios->GetBinError(10,6)*effRatio[2][1][iVar]->GetFunction("pol0")->GetParameter(0)/effRatio[2][0][iVar]->GetFunction("pol0")->GetParameter(0));
      fRatios->SetBinError(9, 1, fRatios->GetBinError(9,1)*effRatio[2][1][iVar]->GetFunction("pol0")->GetParameter(0)/effRatio[2][0][iVar]->GetFunction("pol0")->GetParameter(0));
      fRatios->SetBinError(4, 6, fRatios->GetBinError(4,6)*effRatio[1][1][iVar]->GetFunction("pol0")->GetParameter(0)/effRatio[1][0][iVar]->GetFunction("pol0")->GetParameter(0));
      fRatios->SetBinError(1, 11, fRatios->GetBinError(1,11)*effRatio[0][1][iVar]->GetFunction("pol0")->GetParameter(0)/effRatio[0][0][iVar]->GetFunction("pol0")->GetParameter(0));
      TF2 fitExpo(Form("fitExpo_%.0f_%.0f",cent[iC][0],cent[iC][1]), "TMath::Exp(-2./3.*[0]*x-2.*[1]*y)", -0.5, 9.5,-0.05,1.05,"");
      fRatios->Fit(Form("fitExpo_%.0f_%.0f",cent[iC][0],cent[iC][1]),"SN");
      //chi2 and fit parameter
      std::cout<<"chi2 = "<<fitExpo.GetChisquare()<<std::endl;
      std::cout<<"prob = "<<fitExpo.GetProb()<<std::endl;
      double fit_parameter_0 = fitExpo.GetParameter(0);
      double fit_parameterError_0 = fitExpo.GetParError(0);
      double fit_parameter_1 = fitExpo.GetParameter(1);
      double fit_parameterError_1 = fitExpo.GetParError(1);
      std::cout<<"mu_B / T = "<<fit_parameter_0<<" +/- "<<fit_parameterError_0<<"; mu_I3 / T = "<<fit_parameter_1<<" +/- "<<fit_parameterError_1<<std::endl;
      double temperature = 155.;
      std::cout<<"mu_B (T = 155 MeV) = "<<fit_parameter_0*155<<" +/- "<<fit_parameterError_0*155<<" +/- "<<fit_parameter_0*2<<" MeV; mu_I3 (T = 155 MeV) = "<<fit_parameter_1*155<<" +/- "<<fit_parameterError_1*155<<" +/- "<<fit_parameter_1*2<<" MeV"<<std::endl;
    }
  }
  outFile.Close();
}