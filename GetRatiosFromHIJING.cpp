#include "utils/Utils.h"

using namespace utils;

const int kNCentClasses = 3;
const int kNParticleSpecies = 8;
const double kCentralityClasses[3][2]={{0,5},{5,10},{30,50}};
const char* kAntiMatter[2]={"A","M"};
const char *kParticleSpecies[]={"Pion","Kaon","Proton","Neutron","Deuteron","Lambda","Xi","Omega"};
const double kQuantumNumbers[kNParticleSpecies][2]={{0,1},{1.,0.5},{3.,0.5},{3.,-0.5},{6.,0.},{2.,0},{1.,0.5},{0.,0.}}; // 3*B+S,I

void GetRatiosFromHIJING(){
  TFile fileInput("data/AnalysisResultsWithHIJINGGen.root");
  TFile fileOutput("HIJINGRatios.root","recreate");

  TGraphErrors muBCent, muICent;

  TTList *list= (TTList*)fileInput.Get("nuclei_pion_mcTrue_");
  TH3D *fParticleSpeciesProduction3D = (TH3D*)list->Get("fParticleProd");

  for (int iC=0;iC<kNCentClasses;++iC){
    TH1D fRatios1D(Form("fRatios1D_%.0f_%0.f",kCentralityClasses[iC][0],kCentralityClasses[iC][1]),Form("%.0f-%.0f%%",kCentralityClasses[iC][0],kCentralityClasses[iC][1]),10, -0.5, 9.5);
    TH2D fRatios(Form("fRatios_%.0f_%0.f",kCentralityClasses[iC][0],kCentralityClasses[iC][1]),Form("%.0f-%.0f%%",kCentralityClasses[iC][0],kCentralityClasses[iC][1]),10, -0.5, 9.5, 21, -1.05, 1.05);
    TH1D *fSpeciesSpectrum[2];
    TH1D *fRatio;
    for (int iSpecies=1;iSpecies<kNParticleSpecies*2+1;++iSpecies){
      TString hRatioName=Form("f%sRatio_%.0f_%.0f",kParticleSpecies[(iSpecies-1)/2],kCentralityClasses[iC][0],kCentralityClasses[iC][1]);
      int centBinMin=fParticleSpeciesProduction3D->GetXaxis()->FindBin(kCentralityClasses[iC][0]+0.001);
      int centBinMax=fParticleSpeciesProduction3D->GetXaxis()->FindBin(kCentralityClasses[iC][1]-0.001);
      double yield[2];
      int tmpSpecies = iSpecies;
      for (int iM=1;iM>-1;--iM){
        TString hProjName=Form("f%s%sSpectrum_%.0f_%.0f",kAntiMatter[iM],kParticleSpecies[(iSpecies-1)/2],kCentralityClasses[iC][0],kCentralityClasses[iC][1]);
        if (iM==0) iSpecies+=1;
        fSpeciesSpectrum[iM]=fParticleSpeciesProduction3D->ProjectionZ(hProjName,centBinMin,centBinMax,iSpecies,iSpecies);
        std::cout<<"iSpecies="<<iSpecies<<", CentBinMin="<<centBinMin<<"; CentBinMax="<<centBinMax<<std::endl;
        fSpeciesSpectrum[iM]->Sumw2();
        fSpeciesSpectrum[iM]->Write();
        yield[iM]=fSpeciesSpectrum[iM]->Integral();
      }
      fRatio=new TH1D(*fSpeciesSpectrum[0]);
      fRatio->SetName(hRatioName);
      std::cout<<"Divide: "<<fSpeciesSpectrum[0]->GetName()<<", "<<fSpeciesSpectrum[1]->GetName()<<std::endl;
      fRatio->Divide(fSpeciesSpectrum[0],fSpeciesSpectrum[1],1,1);
      fRatio->SetFillColor(kRed);
      fRatio->SetFillStyle(1001);
      fRatio->Write();
      if (/* (tmpSpecies/2==0)|| */(tmpSpecies/2==2)||(tmpSpecies/2==3)){
        fRatios.SetBinContent(fRatios.GetXaxis()->FindBin(kQuantumNumbers[(tmpSpecies)/2][0]+0.001),fRatios.GetYaxis()->FindBin(kQuantumNumbers[(tmpSpecies)/2][1]+0.001),yield[0]/yield[1]);
        fRatios.SetBinError(fRatios.GetXaxis()->FindBin(kQuantumNumbers[(tmpSpecies)/2][0]+0.001),fRatios.GetYaxis()->FindBin(kQuantumNumbers[(tmpSpecies)/2][1]+0.001),yield[0]/yield[1]*TMath::Sqrt(1/yield[0]+1/yield[1]));
        /* if (kQuantumNumbers[(tmpSpecies)/2][1]>-0.1||kQuantumNumbers[(tmpSpecies)/2][1]<0.6){
          fRatios1D.SetBinContent(fRatios1D.GetXaxis()->FindBin(kQuantumNumbers[(tmpSpecies)/2][0]+0.001),yield[0]/yield[1]);
          fRatios1D.SetBinError(fRatios1D.GetXaxis()->FindBin(kQuantumNumbers[(tmpSpecies)/2][0]+0.001),yield[0]/yield[1]*TMath::Sqrt(1/yield[0]+1/yield[1]));
        } */
      }
    }
    TF2 fitExpo(Form("fitExpo_%.0f_%.0f",kCentralityClasses[iC][0],kCentralityClasses[iC][1]), "TMath::Exp(-2./3.*[0]*x-2.*[1]*y)", -0.5, 9.5,-0.05,1.05,"");
    fRatios.Fit(Form("fitExpo_%.0f_%.0f",kCentralityClasses[iC][0],kCentralityClasses[iC][1]),"SN");
    fRatios.GetZaxis()->SetRangeUser(0.8,1.0);
    fitExpo.SetMinimum(0.8);
    fitExpo.SetMaximum(1.0);
    fitExpo.Write();
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
    fRatios.Write();
    fRatios1D.GetXaxis()->SetTitle("3B+S");
    fRatios1D.GetYaxis()->SetTitle("Ratio");
    fRatios1D.Write();
    muBCent.AddPoint((kCentralityClasses[iC][0]+kCentralityClasses[iC][1])/2,fit_parameter_0*155);
    muBCent.SetPointError(muBCent.GetN()-1,0,fit_parameterError_0*155);
    muICent.AddPoint((kCentralityClasses[iC][0]+kCentralityClasses[iC][1])/2,fit_parameter_1*155);
    muICent.SetPointError(muICent.GetN()-1,0,fit_parameterError_1*155);
  }
  muBCent.SetFillStyle(3145);
  muICent.SetFillStyle(3145);
  muBCent.SetFillColor(kGreen+2);
  muICent.SetFillColor(kGreen+2);
  muBCent.Write("gMuBCent");
  muICent.Write("gMuICent");
  fileOutput.Close();
}