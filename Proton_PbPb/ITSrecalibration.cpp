// ITS recalibration

#include <TFile.h>
#include <Riostream.h>
#include <TH3F.h>
#include <TH1D.h>
#include <TF1.h>
#include <RooRealVar.h>
#include <RooPlot.h>
#include <RooExponential.h>
#include <RooAddPdf.h>
#include <RooAbsPdf.h>

#include "../utils/Config.h"
#include "../utils/Utils.h"
#include "../utils/RooGausExp.h"

using namespace proton;
using namespace utils;

void ITSrecalibration(){
  TFile *inFile=TFile::Open("../data/Proton_PbPb/AnalysisResults_.root");
  if (!inFile) std::cout << "Input file does not exist!" << std::endl;
  TFile outFile("out/ITScalib.root","recreate");
  TTList *ls_proton=(TTList*)inFile->Get("nuclei_proton_");
  TH3F *fNSigmaITSpTCentA=(TH3F*)ls_proton->Get(Form("f%sITSnSigma",kAntimatterMatter[0]));
  TH3F *fNSigmaITSpTCentM=(TH3F*)ls_proton->Get(Form("f%sITSnSigma",kAntimatterMatter[1]));
  TH1D fITSCalibMean("fITSCalibMean",";#it{p}_{T} (GeV/#it{c});Mean (a.u.)",kNPtBins,kPtBins);
  TH1D fITSCalibSigma("fITSCalibSigma",";#it{p}_{T} (GeV/#it{c});Sigma (a.u.)",kNPtBins,kPtBins);
  TH1D fITSCalibMeanRoo("fITSCalibMeanRoo",";#it{p}_{T} (GeV/#it{c});Mean (a.u.)",kNPtBins,kPtBins);
  TH1D fITSCalibSigmaRoo("fITSCalibSigmaRoo",";#it{p}_{T} (GeV/#it{c});Sigma (a.u.)",kNPtBins,kPtBins);
  outFile.mkdir("fits");
  for (int iB=1;iB<kNPtBins+1;iB++){
    if(kPtBins[iB-1]>0.89) continue;
    outFile.cd("fits");
    TH1D* fNSigmaITSpT=fNSigmaITSpTCentA->ProjectionZ(Form("proj_%.2f_%.2f",kPtBins[iB-1],kPtBins[iB]),1,10,iB,iB);
    TH1D* fNSigmaITSpTM=fNSigmaITSpTCentM->ProjectionZ(Form("proj_%.2f_%.2f",kPtBins[iB-1],kPtBins[iB]),1,10,iB,iB);
    fNSigmaITSpT->Add(fNSigmaITSpTM);
    fNSigmaITSpT->Rebin(1);
    double lowNsigma = -3.8;
    if (kPtBins[iB-1]>0.39 && kPtBins[iB-1]<0.6) lowNsigma = -3.2;
    else if (kPtBins[iB-1]>0.6) lowNsigma = -2.7;
    RooRealVar nSigma("nSigma","nSigma",lowNsigma,4.);
    RooDataHist data("data","data",RooArgList(nSigma),fNSigmaITSpT);
    RooPlot *xframe = nSigma.frame();
    RooRealVar mu("mu","mu",-0.5,-1.,1.);
    RooRealVar sig("sig","sig",0.5,0.,1.);
    RooRealVar tau("tau","tau",0.5,0.,2.);
    RooRealVar slope("slope","slope",-10.,10.);
    RooGausExp signal("signal","signal",nSigma,mu,sig,tau);
    RooExponential bkg("bkg","bkg",nSigma,slope);
    RooRealVar w("w","w",0.,1.);
    RooAddPdf model("model","model",signal,bkg,w);
    model.fitTo(data);
    data.plotOn(xframe,RooFit::Name("data"));
    model.plotOn(xframe,RooFit::Name("model"));
    model.plotOn(xframe,RooFit::Components("signal"),RooFit::Name("signal"),RooFit::LineColor(kRed),RooFit::LineStyle(kDashed));
    model.plotOn(xframe,RooFit::Components("bkg"),RooFit::Name("bkg"),RooFit::LineColor(kGreen),RooFit::LineStyle(kDashed));
    TF1 func("func","gaus(0)+expo(3)");
    func.SetParLimits(0,0,1.e9);
    func.SetParLimits(1,-1.5,0.);
    func.SetParLimits(2,0.5,2.);
    func.SetParameter(1,-0.7);
    func.SetParameter(2,0.7);
    func.SetParLimits(4,-20.,10.);
    if (kPtBins[iB-1]<0.5){for (int i=0;i<2;++i)fNSigmaITSpT->Fit("func","RQL","",-3.5,0.2);}
    else if (kPtBins[iB-1]<0.6){for (int i=0;i<2;++i)fNSigmaITSpT->Fit("func","RQL","",-3.,0.2);}
    else{for (int i=0;i<2;++i)fNSigmaITSpT->Fit("func","RQL","",-2.5,0.2);}
    //if (fNSigmaITSpT->GetFunction("func")->GetChisquare()/fNSigmaITSpT->GetFunction("func")->GetNDF() > 3.) continue;
    fITSCalibMean.SetBinContent(iB,fNSigmaITSpT->GetFunction("func")->GetParameter(1));
    fITSCalibMean.SetBinError(iB,fNSigmaITSpT->GetFunction("func")->GetParError(1));
    fITSCalibSigma.SetBinContent(iB,fNSigmaITSpT->GetFunction("func")->GetParameter(2));
    fITSCalibSigma.SetBinError(iB,fNSigmaITSpT->GetFunction("func")->GetParError(2));
    fITSCalibMeanRoo.SetBinContent(iB,mu.getVal());
    fITSCalibMeanRoo.SetBinError(iB,mu.getError());
    fITSCalibSigmaRoo.SetBinContent(iB,sig.getVal());
    fITSCalibSigmaRoo.SetBinError(iB,sig.getError());
    fNSigmaITSpT->Write();
    xframe->Write();
  }
  outFile.cd("");
  fITSCalibMean.Write();
  fITSCalibSigma.Write();
  fITSCalibMeanRoo.Write();
  fITSCalibSigmaRoo.Write();
}