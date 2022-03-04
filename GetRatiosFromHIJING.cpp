#include "utils/Utils.h"

using namespace utils;

const int kNCentClasses = 3;
const int kNParticleSpecies = 8;
const double kCentralityClasses[3][2]={{0,5},{5,10},{30,50}};
const char* kAntiMatter[2]={"A","M"};
const char *kParticleSpecies[]={"Pion","Kaon","Proton","Neutron","Deuteron","Lambda","Xi","Omega"};

void GetRatiosFromHIJING(){
  TFile fileInput("data/AnalysisResultsWithHIJINGGen.root");
  TFile fileOutput("HIJINGRatios.root","recreate");

  TTList *list= (TTList*)fileInput.Get("nuclei_pion_mcTrue_");
  TH3D *fParticleSpeciesProduction3D = (TH3D*)list->Get("fParticleProd");

  for (int iC=0;iC<kNCentClasses;++iC){
    TH1D *fSpeciesSpectrum[2];
    TH1D *fRatio;
    for (int iSpecies=1;iSpecies<kNParticleSpecies*2+1;++iSpecies){
      TString hRatioName=Form("f%sRatio_%.0f_%.0f",kParticleSpecies[iSpecies/2],kCentralityClasses[iC][0],kCentralityClasses[iC][1]);
      int centBinMin=fParticleSpeciesProduction3D->GetXaxis()->FindBin(kCentralityClasses[iC][0]+0.001);
      int centBinMax=fParticleSpeciesProduction3D->GetXaxis()->FindBin(kCentralityClasses[iC][1]-0.001);
      for (int iM=1;iM>-1;--iM){
        TString hProjName=Form("f%s%sSpectrum_%.0f_%.0f",kAntiMatter[iM],kParticleSpecies[(iSpecies-1)/2],kCentralityClasses[iC][0],kCentralityClasses[iC][1]);
        if (iM==0) iSpecies+=1;
        fSpeciesSpectrum[iM]=fParticleSpeciesProduction3D->ProjectionZ(hProjName,centBinMin,centBinMax,iSpecies,iSpecies);
        std::cout<<"iSpecies="<<iSpecies<<", CentBinMin="<<centBinMin<<"; CentBinMax="<<centBinMax<<std::endl;
        fSpeciesSpectrum[iM]->Sumw2();
        fSpeciesSpectrum[iM]->Write();
      }
      fRatio=new TH1D(*fSpeciesSpectrum[0]);
      fRatio->SetName(hRatioName);
      std::cout<<"Divide: "<<fSpeciesSpectrum[0]->GetName()<<", "<<fSpeciesSpectrum[1]->GetName()<<std::endl;
      fRatio->Divide(fSpeciesSpectrum[0],fSpeciesSpectrum[1],1,1);
      fRatio->SetFillColor(kRed);
      fRatio->SetFillStyle(1001);
      fRatio->Write();
    }
  }

  fileOutput.Close();
}