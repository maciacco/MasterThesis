// ITS recalibration

#include <TFile.h>
#include <Riostream.h>
#include <TH3F.h>
#include <TH1D.h>
#include <TF1.h>

#include "../utils/Config.h"
#include "../utils/Utils.h"

using namespace proton;
using namespace utils;

void ITSrecalibration(){
  TFile *inFile=TFile::Open("../data/Proton_PbPb/AnalysisResults_ITS.root");
  if (!inFile) std::cout << "Input file does not exist!" << std::endl;
  TFile outFile("out/ITScalib.root","recreate");
  TTList *ls_proton=(TTList*)inFile->Get("nuclei_proton_");
  for (int iM=0;iM<2;++iM){
    TH3F *fNSigmaITSpTCent=(TH3F*)ls_proton->Get(Form("f%sITSnSigma",kAntimatterMatter[iM]));
    TH1D fITSCalibMean(Form("f%sITSCalibMean",kAntimatterMatter[iM]),";#it{p}_{T} (GeV/#it{c});Mean (a.u.)",kNPtBins,kPtBins);
    TH1D fITSCalibSigma(Form("f%sITSCalibSigma",kAntimatterMatter[iM]),";#it{p}_{T} (GeV/#it{c});Sigma (a.u.)",kNPtBins,kPtBins);
    outFile.mkdir(Form("fits_%s",kAntimatterMatter[iM]));
    for (int iB=1;iB<kNPtBins+1;iB++){
      outFile.cd(Form("fits_%s",kAntimatterMatter[iM]));
      TH1D* fNSigmaITSpT=fNSigmaITSpTCent->ProjectionZ(Form("proj_%.2f_%.2f",kPtBins[iB-1],kPtBins[iB]),1,10,iB,iB);
      fNSigmaITSpT->Fit("gaus","RQL","",-1.5,1.5);
      if (fNSigmaITSpT->GetFunction("gaus")->GetChisquare()/fNSigmaITSpT->GetFunction("gaus")->GetNDF() > 3.) continue;
      fITSCalibMean.SetBinContent(iB,fNSigmaITSpT->GetFunction("gaus")->GetParameter(1));
      fITSCalibMean.SetBinError(iB,fNSigmaITSpT->GetFunction("gaus")->GetParError(1));
      fITSCalibSigma.SetBinContent(iB,fNSigmaITSpT->GetFunction("gaus")->GetParameter(2));
      fITSCalibSigma.SetBinError(iB,fNSigmaITSpT->GetFunction("gaus")->GetParError(2));
      fNSigmaITSpT->Write();
    }
    outFile.cd("");
    fITSCalibMean.Write();
    fITSCalibSigma.Write();
  }
}