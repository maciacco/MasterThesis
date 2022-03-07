const char *species[]={"Pion","Proton","He3"};
const char *antiMatter[]={"A","M"};
const char *var[]={"Minus","Plus"};
const double cent[3][2]={{0,5},{5,10},{30,50}};
const double fitRange[3][2]={{0.5,1.6},{0.8,3.},{1.5,10}};
const bool VERBOSE=false;

void MaterialBudget(){
  TFile fitSHM3D("FinalPlot3D.root");
  TFile outFile("MaterialBudgetUncertainty.root","recreate");
  TFile *f[3][3];
  TH1D *effRatio[3][2][2];
  TH1D *effError[3][2];
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
          for (int iP=1;iP<effRatio[iSp][iM][iVar]->GetNbinsX();++iP){
            if (effRatio[iSp][iM][iVar]->GetBinContent(iP)<1.e-9) effRatio[iSp][iM][iVar]->SetBinError(iP,0.);
          }
          outFile.cd();
          effRatio[iSp][iM][iVar]->Write();
        }
        effError[iSp][iM]=new TH1D(*eff[0]);
        effError[iSp][iM]->SetName(Form("f%s%sEffError_%.0f_%.0f",antiMatter[iM],species[iSp],cent[iC][0],cent[iC][1]));
        effError[iSp][iM]->Add(effRatio[iSp][iM][1],effRatio[iSp][iM][0],-1./sqrt(12),1./sqrt(12));
        if (iC==0 && iM==1) {
          effError[iSp][iM]->SetBinContent(effError[iSp][iM]->FindBin(1.98),0.);
          effError[iSp][iM]->SetBinError(effError[iSp][iM]->FindBin(1.98),0.);
        }
        TF1 fitFunction("fitFunction","[0]*TMath::Power(x,[1])",fitRange[iSp][0],fitRange[iSp][1]);
        effError[iSp][iM]->Fit("fitFunction","QR","",fitRange[iSp][0],fitRange[iSp][1]);
        effError[iSp][iM]->Write();
      }
    }
  }
  outFile.Close();
}