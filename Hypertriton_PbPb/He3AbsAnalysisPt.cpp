#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <Riostream.h>

#include "AliMCEvent.h"
#include "AliVParticle.h"
#include "AliMCEventHandler.h"
#endif

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
    double centClasses[][2] = {{0, 5}, {5, 10}, {30, 50}};
    int centIndex[3] = {0, 1, 3};
    int nFun = 1;
    TF1* fun[3][5];
    std::string names[5]{"BGBW","Boltzmann","Mt-exp","Pt-exp","LevyTsallis"};
    std::string antimattMatt[2]{"anti","matt"};
    TFile inputFun("Anti_fits.root");

    for (int iCent{0}; iCent < nCent; ++iCent) {
        for (int iFun{0}; iFun < nFun; ++iFun) {
            fun[iCent][iFun] = (TF1*)inputFun.Get(Form("%s/%d/%s%d",names[iFun].c_str(),centIndex[iCent],names[iFun].c_str(),centIndex[iCent]));
        }
    }

    TH1D *recPt[3][2][4], *recCt[3][2][4], *genPt[3][2][4], *genCt[3][2][4];
    double maxPt[3][4];

    double ctBins[] = {2, 4, 8, 14, 35};

    for (int iCent{0}; iCent < nCent; ++iCent) {
        for (int iFun{0}; iFun < nFun; ++iFun) {
            for (int iMatt{0}; iMatt < 2; ++iMatt) {
                recPt[iCent][iMatt][iFun] = new TH1D(Form("recPt_%.0f_%.0f_%s_%s",centClasses[iCent][0],centClasses[iCent][1],antimattMatt[iMatt].data(),names[iFun].data()), "; #it{p}_{T} (GeV/#it{c}); Counts", 50, 0, 10);
                recCt[iCent][iMatt][iFun] = new TH1D(Form("recCt_%.0f_%.0f_%s_%s",centClasses[iCent][0],centClasses[iCent][1],antimattMatt[iMatt].data(),names[iFun].data()), "; #it{c}t (cm); Counts", 4, ctBins);
                genPt[iCent][iMatt][iFun] = new TH1D(Form("genPt_%.0f_%.0f_%s_%s",centClasses[iCent][0],centClasses[iCent][1],antimattMatt[iMatt].data(),names[iFun].data()), "; #it{p}_{T} (GeV/#it{c}); Counts", 50, 0, 10);
                genCt[iCent][iMatt][iFun] = new TH1D(Form("genCt_%.0f_%.0f_%s_%s",centClasses[iCent][0],centClasses[iCent][1],antimattMatt[iMatt].data(),names[iFun].data()), "; #it{c}t (cm);Counts;", 4, ctBins);
            }
            maxPt[iCent][iFun] = fun[iCent][iFun]->GetMaximum();
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
                    if (part->IsPhysicalPrimary() && std::abs(part->PdgCode()) == 1000020030)
                    {
                        if (part->Pt() < 2) continue; // pt cut

                        for (int iCent{0}; iCent < nCent; ++iCent) {
                            for (int iFun{0}; iFun < nFun; ++iFun) {
                                float hypPtShapeNum = fun[iCent][iFun]->Eval(part->Pt());
                                if (gRandom->Rndm() * maxPt[iCent][iFun] > hypPtShapeNum) {
                                    continue;
                                }
                                kNorm++;

                                int isMatter = 0;
                                if (part->PdgCode() > 0) isMatter = 1;

                                genPt[iCent][isMatter][iFun]->Fill(part->Pt());
                                for (int c = part->GetDaughterFirst(); c < part->GetDaughterLast(); c++)
                                {
                                    AliVParticle *dPart = mcEv->GetTrack(c);
                                    int dPartPDG = std::abs(dPart->PdgCode());
                                    if (dPartPDG != 22 && dPartPDG != 11)
                                    {
                                        double absCt = ComputeHe3Ct(part, dPart);
                                        double decCt = gRandom->Exp(7.6);
                                        genCt[iCent][isMatter][iFun]->Fill(decCt);
                                        bool isAbsorbed = absCt < decCt;
                                        if (!isAbsorbed) {
                                            recPt[iCent][isMatter][iFun]->Fill(part->Pt());
                                            recCt[iCent][isMatter][iFun]->Fill(decCt);
                                        }
                                        break;
                                    }
                                }
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
            }
        }
    }
    fFile.Close();
}
