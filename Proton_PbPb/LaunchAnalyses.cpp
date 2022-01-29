// LaunchAnalyses.cpp
// This macro launches all the analyses needed for systematics computation

#include <stdlib.h>
#include <string.h>

#include <Riostream.h>
#include <TStopwatch.h>
#include <TROOT.h>
#include <TSystem.h>

#include "../utils/Config.h"

using namespace proton;

double roi_n_sigma[] = {7.5, 8., 8.5};

void LaunchAnalyses(const bool analyse = false)
{
  TStopwatch swatch;
  swatch.Start(true);
  std::ofstream outFile;
  outFile.open("ProcessedFiles.dat");

  auto tmpNCutDCAz = kNCutDCAz-1;
  auto tmpNTPCPidSigmas = kNTPCPidSigmas-1;
  auto tmpNCutTPCClusters = kNCutTPCClusters-1;
  auto tmpNCutDCAxy = kNCutDCAxy-1;

  for (int iCutSettings = -1; iCutSettings < tmpNCutDCAz + tmpNTPCPidSigmas + tmpNCutTPCClusters/*  + tmpNCutDCAxy */; ++iCutSettings)
  {
    char hname[100];

    // bool binCountingFlag = 1 - iBin;
    int cutVariable = 0;
    int cutIndex = iCutSettings;
    char *fullCutSettings = Form("999");
    if (iCutSettings > -1) {
      if (iCutSettings >= tmpNCutDCAz && iCutSettings < (tmpNCutDCAz + tmpNTPCPidSigmas))
      {
        cutVariable = 1;
        cutIndex -= tmpNCutDCAz;
      }
      else if (iCutSettings >= (tmpNCutDCAz + tmpNTPCPidSigmas) && iCutSettings < (tmpNCutDCAz + tmpNTPCPidSigmas + tmpNCutTPCClusters))
      {
        cutVariable = 2;
        cutIndex -= (tmpNCutDCAz + tmpNTPCPidSigmas);
      }
      else if (iCutSettings >= (tmpNCutDCAz + tmpNTPCPidSigmas + tmpNCutTPCClusters))
      {
        cutVariable = 3;
        cutIndex -= (tmpNCutDCAz + tmpNTPCPidSigmas + tmpNCutTPCClusters);
      }
      fullCutSettings = Form("%s%d", cutSettings[cutVariable], cutIndex);
    }
    std::cout << "fullCutSettings = " << fullCutSettings << std::endl;

    double DCAxyCut = 0.12;
    if (cutVariable == 3) DCAxyCut = kCutDCAxyVariations[cutIndex];

    for (int iBkg = 1; iBkg < 2; ++iBkg)
    {
      for(int iNsigma = 0; iNsigma < 3; ++iNsigma) {
        for (int iG3G4Prim = 0; iG3G4Prim < 2; ++iG3G4Prim){
          std::cout << "bkg selection = " << iBkg << "; roiNsigma = " << roi_n_sigma[iNsigma] << "; dcaxycut = " << DCAxyCut << "; G3G4Prim = " << iG3G4Prim << std::endl;
          if (analyse)
          {
            gSystem->Exec(Form("bash ~/Code/MasterThesis/Proton_PbPb/scripts/LaunchAnalysisSignEffPrim.sh %s 1 1 %f %f %d", fullCutSettings, roi_n_sigma[iNsigma], DCAxyCut, iG3G4Prim));
          }

          for (int iSgm = 0; iSgm < 2; ++iSgm)
          {
            // bool binCountingFlag = 1 - iBin;
            bool sigmoidFlag = 1 - iSgm;
            auto tmpFullCutSettings = fullCutSettings;
            if (iCutSettings == -1) tmpFullCutSettings = Form("");
            auto spectraNameId = Form("%s_%d_%d_%d_%d",tmpFullCutSettings, iBkg, sigmoidFlag, iNsigma, iG3G4Prim);
            std::cout << "SigmoidCorrection = " << kBoolString[sigmoidFlag] << "; cutSettings = " << fullCutSettings << "..." << std::endl;
            outFile << "SigmoidCorrection = " << kBoolString[sigmoidFlag] << "; cutSettings = " << fullCutSettings << "..."
                    << "\n";

            if (analyse)
            {
              gSystem->Exec(Form("bash ~/Code/MasterThesis/Proton_PbPb/scripts/LaunchAnalysisSpec.sh %s 1 %d %d %s %f %d", fullCutSettings, iBkg, sigmoidFlag, spectraNameId, roi_n_sigma[iNsigma], iG3G4Prim));
            }
          }
        }
      }
    }

  }

  swatch.Stop();
  std::cout << "Elapsed (real) time = " << swatch.RealTime() << " s" << std::endl;

  outFile.close();
}