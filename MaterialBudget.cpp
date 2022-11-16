const char *species[]={"Pion","Proton","He3","Triton","Proton"};
const char *antiMatter[]={"A","M"};
const char *var[]={"Minus","Plus"};
const char *var_title[]={"-4.5%","+4.5%"};
const double cent[3][2]={{0,5},{5,10},{30,50}};
const double fitRange[5][2]={{0.5,1.5},{0.8,3.},{1.5,10},{1.6,3},{0.5,.85}};
const bool VERBOSE=false;

void MaterialBudget(){
  gStyle->SetOptStat(0);
  TFile fitSHM3D("FinalPlot3D.root");
  TFile outFile("MaterialBudgetUncertainty.root","recreate");
  TFile *f[5][3];
  TH1D *effRatio[5][2][2];
  TH1D *effError[5][2];
  for (int iSp=0; iSp<5; ++iSp){
    for (int iC=0;iC<3;++iC){
      for (int iM=0; iM<2; ++iM){
        TH1D *eff[3];
        for (int iProd=0; iProd<3; ++iProd){
          if (iSp<4)f[iSp][iProd]=new TFile(Form("%s_PbPb/out/Efficiency%s_LHC22b9_%d.root",species[iSp],species[iSp],iProd+1));
          else f[iSp][iProd]=new TFile(Form("%s_PbPb/out/Efficiency%s_LowPt_LHC22b9_%d.root",species[iSp],species[iSp],iProd+1));
          if (!f[iSp][iProd]) continue;
          TString detector("TOF");
          TString directory("_/");
          if (iSp>0&&iSp<4&&iSp!=1) directory=Form("");
          if (iSp>1&&iSp<4) {
            detector=Form("TPC");
            directory=Form("");
          }
          if (VERBOSE) {
            std::cout<<"File Name = "<<f[iSp][iProd]->GetName()<<std::endl;
            std::cout<<"Efficiency name = "<<Form("%sf%sEff_%s_%.0f_%.0f",directory.Data(),antiMatter[iM],detector.Data(),cent[iC][0],cent[iC][1])<<std::endl;
          }
          eff[iProd]=(TH1D*)f[iSp][iProd]->Get(Form("%sf%sEff_%s_%.0f_%.0f",directory.Data(),antiMatter[iM],detector.Data(),0.,90./* ,cent[iC][0],cent[iC][1] */));
          //if (iSp==0)eff[iProd]=(TH1D*)f[iSp][iProd]->Get(Form("%sf%sEff_%s_0_90",directory.Data(),antiMatter[iM],detector.Data()));
          if (!eff[iProd]) continue;
          //outFile.cd();
          //eff[iProd]->Write();
        }
        for (int iVar=0; iVar<2; ++iVar){
          effRatio[iSp][iM][iVar]=new TH1D(*eff[0]);
          if (iSp<4) effRatio[iSp][iM][iVar]->SetName(Form("f%s%sEffRatio%s_%.0f_%.0f",antiMatter[iM],species[iSp],var[iVar],cent[iC][0],cent[iC][1]));
          else if (iSp==4) effRatio[iSp][iM][iVar]->SetName(Form("f%s%sLowPTEffRatio%s_%.0f_%.0f",antiMatter[iM],species[iSp],var[iVar],cent[iC][0],cent[iC][1]));
          iVar == 0 ? effRatio[iSp][iM][iVar]->Divide(eff[0],eff[1]) : effRatio[iSp][iM][iVar]->Divide(eff[2],eff[1]);
          for (int iP=1;iP<effRatio[iSp][iM][iVar]->GetNbinsX();++iP){
            if (effRatio[iSp][iM][iVar]->GetBinContent(iP)<1.e-9) effRatio[iSp][iM][iVar]->SetBinError(iP,0.);
          }
          outFile.cd();
          effRatio[iSp][iM][iVar]->SetTitle(Form("Efficiency ratio (%s) %s, Material %s, 0-90%%",antiMatter[iM],species[iSp],var_title[iVar]));
          effRatio[iSp][iM][iVar]->GetXaxis()->SetRangeUser(fitRange[iSp][0],fitRange[iSp][1]);
          TCanvas ratio(effRatio[iSp][iM][iVar]->GetName(),effRatio[iSp][iM][iVar]->GetTitle());
          effRatio[iSp][iM][iVar]->GetYaxis()->SetTitle("(#epsilon#times A)^{var.}/(#epsilon#times A)^{std.}");
          effRatio[iSp][iM][iVar]->Draw();
          ratio.Print(Form("%s.pdf",effRatio[iSp][iM][iVar]->GetName()));
          effRatio[iSp][iM][iVar]->Write();
        }
        effError[iSp][iM]=new TH1D(*eff[0]);
        if (iSp < 4) effError[iSp][iM]->SetName(Form("f%s%sEffError_%.0f_%.0f",antiMatter[iM],species[iSp],cent[iC][0],cent[iC][1]));
        else if (iSp==4) effError[iSp][iM]->SetName(Form("f%s%sLowPTEffError_%.0f_%.0f",antiMatter[iM],species[iSp],cent[iC][0],cent[iC][1]));
        effError[iSp][iM]->Add(effRatio[iSp][iM][1],effRatio[iSp][iM][0],-1./sqrt(12),1./sqrt(12));
        for (int iB=0;iB<effError[iSp][iM]->GetNbinsX()+1;++iB){
          effError[iSp][iM]->SetBinContent(iB,std::abs(effError[iSp][iM]->GetBinContent(iB)));
        }
        if (iC==0 && iM==1) {
          effError[iSp][iM]->SetBinContent(effError[iSp][iM]->FindBin(1.98),0.);
          effError[iSp][iM]->SetBinError(effError[iSp][iM]->FindBin(1.98),0.);
          if (iSp==0){
            effError[iSp][iM]->SetBinContent(effError[iSp][iM]->FindBin(1.12),0.);
            effError[iSp][iM]->SetBinError(effError[iSp][iM]->FindBin(1.12),0.);
          }
        }
        effError[iSp][iM]->GetXaxis()->SetRangeUser(fitRange[iSp][0],fitRange[iSp][1]);
        TF1 fitFunction("fitFunction","[0]*TMath::Power(x,[1])",fitRange[iSp][0],fitRange[iSp][1]);
        effError[iSp][iM]->Fit("fitFunction","QR","",fitRange[iSp][0],fitRange[iSp][1]);
        effError[iSp][iM]->Write();
        TCanvas c(Form("c%s",effError[iSp][iM]->GetName()),Form("c%s",effError[iSp][iM]->GetName()));
        effError[iSp][iM]->SetTitle(" ");
        effError[iSp][iM]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        effError[iSp][iM]->GetYaxis()->SetTitle("#epsilon #times A relative uncertainty");
        effError[iSp][iM]->Draw();
        //c.Write();
        c.Print(Form("%s.pdf",effError[iSp][iM]->GetName()));
      }
    }
  }
  outFile.Close();
}