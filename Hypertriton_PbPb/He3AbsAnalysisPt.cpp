#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <Riostream.h>

#include "AliMCEvent.h"
#include "AliVParticle.h"
#include "AliMCEventHandler.h"
#endif

# define DEBUG_ABSORPTION

double kHe3Mass = 2.809230089;

constexpr double Sq(double x) { return x * x; }

double Dist(const double a[3], const double b[3]) { return std::sqrt(Sq(a[0] - b[0]) + Sq(a[1] - b[1]) + Sq(a[2] - b[2])); }

double ComputeHe3Ct(AliVParticle *he3Part, AliVParticle *dauPart)
{
    double primVertex[3];
    double secVertex[3];
    he3Part->XvYvZv(primVertex);
    dauPart->XvYvZv(secVertex);
    double decLength = Dist(secVertex, primVertex);
    // std::cout<<kHe3Mass * decLength / he3Part->P()<<std::endl;
    return kHe3Mass * decLength / he3Part->P();
}


void He3AbsAnalysisPt(std::string pathToSimulation = "./")
{
    int nCent = 3;
    int nCentMB = 3;
    double centClasses[][2] = {{0, 5}, {5, 10}, {30, 50}};
    int centIndex[3] = {1, 2, 4};

    int nFun = 1;
    int nFunMB = 1;

    TF1* fun[3][5];
    TF1* funMB[3][4];

    std::string names[5]{"BGBW","Boltzmann","Mt-exp","Pt-exp","LevyTsallis"};
    std::string namesMB[4]{"BlastWave","Boltzmann","LevyTsallis","Mt-exp"};
    std::string antimattMatt[2]{"anti","matt"};

    TFile inputFun("Anti_fits.root");
    TFile inputFunMB("BlastWaveFits.root");
    TFile inputCent("AnalysisResults_18.root"); // centrality distribution (for MB analysis)

    // get functions/histograms from files
    // analysis in centrality classes
    for (int iCent{0}; iCent < nCent; ++iCent) {
        for (int iFun{0}; iFun < nFun; ++iFun) {
            fun[iCent][iFun] = (TF1*)inputFun.Get(Form("%s/%d/%s%d",names[iFun].c_str(),centIndex[iCent],names[iFun].c_str(),centIndex[iCent]));
        }
    }
    // MB analysis
    for (int iCent{0}; iCent < 3; ++iCent) {
        for (int iFun{0}; iFun < nFunMB; ++iFun) {
            funMB[iCent][iFun] = (TF1*)inputFunMB.Get(Form("%s/%s%d",namesMB[iFun].c_str(),namesMB[iFun].c_str(),iCent));
        }
    }
    TH1D *centHist = (TH1D*)inputCent.Get("Centrality_selected");

    // declare histograms
    TH1D *recPt[3][2][5], *recCt[3][2][5], *genPt[3][2][5], *genCt[3][2][5];
    TH1D *absorptionCt[3][2][5];
    TH1D *daughterCharge[3][2][5];
    TH2F *recPtCt[3][2][5], *genPtCt[3][2][5];
    TH1D *recPtNotAbsorbed[3][2][5], *nDaughters[3][2][5], *daughterPDGCode[3][2][5];
    TH1D *recInt[3][2][5], *genInt[3][2][5];                             // pt and ct integrated analysis
    TH1D *recPtMB[2][4], *recCtMB[2][4], *genPtMB[2][4], *genCtMB[2][4]; // MB analysis

    double maxPt[3][5];
    double maxPtMB[3][4];
    double maxEv; // event centrality distribution from data

    double ctBins[] = {0, 2, 4, 8, 14, 35};
    double ptBins[51];

    for (int i=0; i<51; ++i) ptBins[i] = 0.2*i;

    // book histograms
    // analysis in centrality classes
    for (int iCent{0}; iCent < nCent; ++iCent) {
        for (int iFun{0}; iFun < nFun; ++iFun) {
            for (int iMatt{0}; iMatt < 2; ++iMatt) {
                recPt[iCent][iMatt][iFun] = new TH1D(Form("recPt_%.0f_%.0f_%s_%s",centClasses[iCent][0],centClasses[iCent][1],antimattMatt[iMatt].data(),names[iFun].data()), "; #it{p}_{T} (GeV/#it{c}); Counts", 50, 0, 10);
                recCt[iCent][iMatt][iFun] = new TH1D(Form("recCt_%.0f_%.0f_%s_%s",centClasses[iCent][0],centClasses[iCent][1],antimattMatt[iMatt].data(),names[iFun].data()), "; #it{c}t (cm); Counts", 5, ctBins);
                genPt[iCent][iMatt][iFun] = new TH1D(Form("genPt_%.0f_%.0f_%s_%s",centClasses[iCent][0],centClasses[iCent][1],antimattMatt[iMatt].data(),names[iFun].data()), "; #it{p}_{T} (GeV/#it{c}); Counts", 50, 0, 10);
                genCt[iCent][iMatt][iFun] = new TH1D(Form("genCt_%.0f_%.0f_%s_%s",centClasses[iCent][0],centClasses[iCent][1],antimattMatt[iMatt].data(),names[iFun].data()), "; #it{c}t (cm);Counts;", 5, ctBins);

                daughterCharge[iCent][iMatt][iFun] = new TH1D(Form("daughterCharge_%.0f_%.0f_%s_%s",centClasses[iCent][0],centClasses[iCent][1],antimattMatt[iMatt].data(),names[iFun].data()), "; daughterCharge (#it{e});Counts;", 11, -5, 5);

                recInt[iCent][iMatt][iFun] = new TH1D(Form("recInt_%.0f_%.0f_%s_%s",centClasses[iCent][0],centClasses[iCent][1],antimattMatt[iMatt].data(),names[iFun].data()), "; #it{c}t (cm);Counts;", 1, 0, 35);
                genInt[iCent][iMatt][iFun] = new TH1D(Form("genInt_%.0f_%.0f_%s_%s",centClasses[iCent][0],centClasses[iCent][1],antimattMatt[iMatt].data(),names[iFun].data()), "; #it{c}t (cm);Counts;", 1, 0, 35);

                recPtCt[iCent][iMatt][iFun] = new TH2F(Form("recPtCt_%.0f_%.0f_%s_%s",centClasses[iCent][0],centClasses[iCent][1],antimattMatt[iMatt].data(),names[iFun].data()), "; #it{p}_{T} (GeV/#it{c}); #it{c}t (cm);", 50, ptBins, 5, ctBins);
                genPtCt[iCent][iMatt][iFun] = new TH2F(Form("genPtCt_%.0f_%.0f_%s_%s",centClasses[iCent][0],centClasses[iCent][1],antimattMatt[iMatt].data(),names[iFun].data()), "; #it{p}_{T} (GeV/#it{c}); #it{c}t (cm);", 50, ptBins, 5, ctBins);

                absorptionCt[iCent][iMatt][iFun] = new TH1D(Form("absorptionCt_%.0f_%.0f_%s_%s",centClasses[iCent][0],centClasses[iCent][1],antimattMatt[iMatt].data(),names[iFun].data()), "; #it{c}t (cm);Counts;", 301, -1, 300);
                recPtNotAbsorbed[iCent][iMatt][iFun] = new TH1D(Form("recPtNotAbsorbed_%.0f_%.0f_%s_%s",centClasses[iCent][0],centClasses[iCent][1],antimattMatt[iMatt].data(),names[iFun].data()), "; #it{p}_{T} (GeV/#it{c}); Counts", 50, 0, 10);
                nDaughters[iCent][iMatt][iFun] = new TH1D(Form("nDaughters_%.0f_%.0f_%s_%s",centClasses[iCent][0],centClasses[iCent][1],antimattMatt[iMatt].data(),names[iFun].data()), "; nDaughters; Counts;", 1001, 0, 1000);
                daughterPDGCode[iCent][iMatt][iFun] = new TH1D(Form("daughterPDGCode_%.0f_%.0f_%s_%s",centClasses[iCent][0],centClasses[iCent][1],antimattMatt[iMatt].data(),names[iFun].data()), "; nDaughters; Counts;", 100, 0, 99);
            }
            maxPt[iCent][iFun] = fun[iCent][iFun]->GetMaximum();
        }
    }
    // MB analysis
    for (int iFun{0}; iFun < nFunMB; ++iFun) {
        for (int iMatt{0}; iMatt < 2; ++iMatt) {
            recPtMB[iMatt][iFun] = new TH1D(Form("recPtMB_%s_%s",antimattMatt[iMatt].data(),namesMB[iFun].data()), "; #it{p}_{T} (GeV/#it{c}); Counts", 50, 0, 10);
            recCtMB[iMatt][iFun] = new TH1D(Form("recCtMB_%s_%s",antimattMatt[iMatt].data(),namesMB[iFun].data()), "; #it{c}t (cm);Counts;", 5, ctBins);
            genPtMB[iMatt][iFun] = new TH1D(Form("genPtMB_%s_%s",antimattMatt[iMatt].data(),namesMB[iFun].data()), "; #it{p}_{T} (GeV/#it{c}); Counts", 50, 0, 10);
            genCtMB[iMatt][iFun] = new TH1D(Form("genCtMB_%s_%s",antimattMatt[iMatt].data(),namesMB[iFun].data()), "; #it{c}t (cm);Counts;", 5, ctBins);
        }
        for (int iCent{0}; iCent < nCentMB; ++iCent) {
            maxPtMB[iCent][iFun] = funMB[iCent][iFun]->GetMaximum();
        }
    }

    int kNorm = 0;
    int dirs[] = {295585, 295586};
    int subdirs[] = {38, 160};

    for (Int_t iDir = 0; iDir <= 1; iDir++)
    {
        for (Int_t subdir = 1; subdir <= subdirs[iDir]; ++subdir)
        {
            AliMCEventHandler mcEvHandler("mcEvHandler", "MC Event Handler");
            std::string pathToDir = pathToSimulation + "/" + Form("%d/%03d/", dirs[iDir], subdir);
            std::cout << "Analysing " << pathToDir << std::endl;
            mcEvHandler.SetInputPath(pathToDir.data());
            mcEvHandler.SetReadTR(false);
            mcEvHandler.Init("");
            Int_t iEvent = 0;

            while (mcEvHandler.LoadEvent(iEvent++))
            {
                AliMCEvent *mcEv = mcEvHandler.MCEvent();

                for (Int_t i = 0; i < mcEv->GetNumberOfTracks(); ++i)
                {
                    AliVParticle *part = mcEv->GetTrack(i);
                    if (part->Pt() < 2) continue; // pt cut

                    if (part->IsPhysicalPrimary() && std::abs(part->PdgCode()) == 1000020030)
                    {
                        int isMatter = 0;
                        if (part->PdgCode() > 0) isMatter = 1;

                        // analysis performed in centrality classes
                        for (int iCent{0}; iCent < nCent; ++iCent) {
                            for (int iFun{0}; iFun < nFun; ++iFun) {
                                float hypPtShapeNum = fun[iCent][iFun]->Eval(part->Pt());
                                if (gRandom->Rndm() * maxPt[iCent][iFun] > hypPtShapeNum) {
                                    continue;
                                }
                                kNorm++;

                                double decCt = gRandom->Exp(7.6);
                                genPt[iCent][isMatter][iFun]->Fill(part->Pt());
                                genCt[iCent][isMatter][iFun]->Fill(decCt);
                                genInt[iCent][isMatter][iFun]->Fill(decCt);
                                genPtCt[iCent][isMatter][iFun]->Fill(part->Pt(),decCt);
                                int checkedDaughters=0;
                                bool isAbsorbed = false;
                                double absCt = -1;
                                if(part->GetNDaughters()>0){
                                    for (int c = part->GetDaughterFirst(); c <= part->GetDaughterLast(); ++c)
                                    {
                                        AliVParticle *dPart = mcEv->GetTrack(c);
                                        int dPartPDG = std::abs(dPart->PdgCode());
                                        if (dPartPDG != 22 && dPartPDG != 11)
                                        {
                                            absCt = ComputeHe3Ct(part, dPart);
                                            isAbsorbed = absCt < decCt;
                                            daughterCharge[iCent][isMatter][iFun]->Fill(dPart->Charge());
                                            if (!isAbsorbed) {
                                                recPt[iCent][isMatter][iFun]->Fill(part->Pt());
                                                recCt[iCent][isMatter][iFun]->Fill(decCt);
                                                recInt[iCent][isMatter][iFun]->Fill(decCt);
                                                recPtCt[iCent][isMatter][iFun]->Fill(part->Pt(),decCt);
                                            }
                                            break;
                                        }
                                        checkedDaughters++;
                                    }
                                }

                                absorptionCt[iCent][isMatter][iFun]->Fill(absCt);
                                // not-absorbed particles
                                if (absCt < -0.5) {
                                    recPt[iCent][isMatter][iFun]->Fill(part->Pt());
                                    recCt[iCent][isMatter][iFun]->Fill(decCt);
                                    recInt[iCent][isMatter][iFun]->Fill(decCt);
                                    recPtCt[iCent][isMatter][iFun]->Fill(part->Pt(),decCt);

                                    #ifdef DEBUG_ABSORPTION
                                        nDaughters[iCent][isMatter][iFun]->Fill(part->GetNDaughters());
                                        if(part->GetNDaughters()>0){
                                            for (int c = part->GetDaughterFirst(); c <= part->GetDaughterLast(); ++c)
                                            {
                                                AliVParticle *dPart = mcEv->GetTrack(c);
                                                int dPartPDG = std::abs(dPart->PdgCode());
                                                daughterPDGCode[iCent][isMatter][iFun]->Fill(dPartPDG);
                                            }
                                        }
                                        recPtNotAbsorbed[iCent][isMatter][iFun]->Fill(part->Pt());
                                    #endif // DEBUG_ABSORPTION
                                }
                            }
                        }

                        // MB analysis
                        for (int iFun{0}; iFun < nFunMB; ++iFun) {
                            // sample centrality
                            float cent = (gRandom->Rndm() * 90.);
                            float evFromUniformCent = centHist->GetBinContent(centHist->FindBin(cent));
                            if (gRandom->Rndm() * maxEv > evFromUniformCent) {
                                continue;
                            }
                            int iCent = 0;
                            if (cent > 10 && cent < 40) iCent = 1;
                            else if (cent > 40) iCent = 2; 

                            float hypPtShapeNum = funMB[iCent][iFun]->Eval(part->Pt());
                            if (gRandom->Rndm() * maxPtMB[iCent][iFun] > hypPtShapeNum) {
                                continue;
                            }
                            kNorm++;

                            double decCt = gRandom->Exp(7.6);
                            genPtMB[isMatter][iFun]->Fill(part->Pt());
                            genCtMB[isMatter][iFun]->Fill(decCt);
                            int checkedDaughters=0;
                            if(part->GetNDaughters()>0){
                                for (int c = part->GetDaughterFirst(); c <= part->GetDaughterLast(); ++c)
                                {
                                    AliVParticle *dPart = mcEv->GetTrack(c);
                                    int dPartPDG = std::abs(dPart->PdgCode());
                                    if (dPartPDG != 22 && dPartPDG != 11)
                                    {
                                        double absCt = ComputeHe3Ct(part, dPart);
                                        bool isAbsorbed = absCt < decCt;
                                        if (!isAbsorbed) {
                                            recPtMB[isMatter][iFun]->Fill(part->Pt());
                                            recCtMB[isMatter][iFun]->Fill(decCt);
                                        }
                                        break;
                                    }
                                    checkedDaughters++;
                                }
                            }

                            // not-absorbed particles
                            if (checkedDaughters==part->GetNDaughters()) {
                                recPtMB[isMatter][iFun]->Fill(part->Pt());
                                recCtMB[isMatter][iFun]->Fill(decCt);
                            }
                        }

                    }
                }
            }
        }
    }

    TFile fFile((pathToSimulation + "recPtHe3_1.5.root").data(), "recreate");
    std::cout << "Total number passing the selections " << kNorm << std::endl;

    for (int iCent{0}; iCent < nCent; ++iCent) {
        fFile.mkdir(Form("%.0f_%.0f",centClasses[iCent][0],centClasses[iCent][1]));
        fFile.cd(Form("%.0f_%.0f",centClasses[iCent][0],centClasses[iCent][1]));
        for (int iMatt{0}; iMatt < 2; ++iMatt) {
            for (int iFun{0}; iFun < nFun; ++iFun) {
                recPt[iCent][iMatt][iFun]->Write();
                recCt[iCent][iMatt][iFun]->Write();
                genPt[iCent][iMatt][iFun]->Write();
                genCt[iCent][iMatt][iFun]->Write();
                TGraphAsymmErrors effPt(recPt[iCent][iMatt][iFun],genPt[iCent][iMatt][iFun]);
                effPt.Write(Form("effPt_%.0f_%.0f_%s_%s",centClasses[iCent][0],centClasses[iCent][1],antimattMatt[iMatt].data(),names[iFun].data()));
                TGraphAsymmErrors effCt(recCt[iCent][iMatt][iFun],genCt[iCent][iMatt][iFun]);
                effCt.Write(Form("effCt_%.0f_%.0f_%s_%s",centClasses[iCent][0],centClasses[iCent][1],antimattMatt[iMatt].data(),names[iFun].data()));
                absorptionCt[iCent][iMatt][iFun]->Write();
                recPtCt[iCent][iMatt][iFun]->Write();
                recPtNotAbsorbed[iCent][iMatt][iFun]->Write();
                genPtCt[iCent][iMatt][iFun]->Write();
                daughterCharge[iCent][iMatt][iFun]->Write();
                recInt[iCent][iMatt][iFun]->Write();
                genInt[iCent][iMatt][iFun]->Write();
                TGraphAsymmErrors effInt(recInt[iCent][iMatt][iFun],genInt[iCent][iMatt][iFun]);
                effInt.Write(Form("effInt_%.0f_%.0f_%s_%s",centClasses[iCent][0],centClasses[iCent][1],antimattMatt[iMatt].data(),names[iFun].data()));

                nDaughters[iCent][iMatt][iFun]->Write();
                daughterPDGCode[iCent][iMatt][iFun]->Write();
            }
        }
    }

    fFile.mkdir("MB");
    fFile.cd("MB");
    for (int iMatt{0}; iMatt < 2; ++iMatt) {
        for (int iFun{0}; iFun < nFunMB; ++iFun) {
            recPtMB[iMatt][iFun]->Write();
            recCtMB[iMatt][iFun]->Write();
            genPtMB[iMatt][iFun]->Write();
            genCtMB[iMatt][iFun]->Write();
            TGraphAsymmErrors effPtMB(recPtMB[iMatt][iFun],genPtMB[iMatt][iFun]);
            effPtMB.Write(Form("effPtMB_%s_%s",antimattMatt[iMatt].data(),namesMB[iFun].data()));
            TGraphAsymmErrors effCtMB(recCtMB[iMatt][iFun],genCtMB[iMatt][iFun]);
            effCtMB.Write(Form("effCtMB_%s_%s",antimattMatt[iMatt].data(),namesMB[iFun].data()));
        }
    }

    fFile.Close();
}
