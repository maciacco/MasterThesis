// LaunchAnalyses.cpp
// This macro launches all the analyses needed for systematics computation

#include <stdlib.h>
#include <string.h>

#include <Riostream.h>
#include <TStopwatch.h>
#include <TROOT.h>
#include <TSystem.h>

#include "../utils/Config.h"

using namespace pion;

double roi_n_sigma_down[] = {1.0,1.5,2.0};
double roi_n_sigma_up[] = {10.,11.,12.};

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

  for (int iCutSettings = -1; iCutSettings < tmpNCutDCAz + tmpNTPCPidSigmas + tmpNCutTPCClusters /* + tmpNCutDCAxy */; ++iCutSettings)
  {
    char hname[100];

    for (int iBkg = 1; iBkg < 2; ++iBkg)
    {
      for(int iNsigmaUp = 0; iNsigmaUp < 3; ++iNsigmaUp) {
        for(int iNsigmaDown = 0; iNsigmaDown < 3; ++iNsigmaDown) {
          
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
          //if (cutVariable == 3) DCAxyCut = kCutDCAxyVariations[cutIndex];

          std::cout << "bkg selection = " << iBkg << "; roiNsigmaDown = " << roi_n_sigma_down[iNsigmaDown] << "; roiNsigmaUp = " << roi_n_sigma_up[iNsigmaUp] << "; dcaxycut = " << DCAxyCut << std::endl;
          if (analyse)
          {
            gSystem->Exec(Form("bash ~/Code/MasterThesis/Pion_PbPb/scripts/LaunchAnalysisSignEffPrim.sh %s 1 1 %f %f %f", fullCutSettings, roi_n_sigma_down[iNsigmaDown], roi_n_sigma_up[iNsigmaUp], DCAxyCut));
          }

          for (int iSgm = 0; iSgm < 1; ++iSgm)
          {
            // bool binCountingFlag = 1 - iBin;
            bool sigmoidFlag = 1 - iSgm;
            auto tmpFullCutSettings = fullCutSettings;
            if (iCutSettings == -1) tmpFullCutSettings = Form("");
            auto spectraNameId = Form("%s_%d_%d_%d_%d",tmpFullCutSettings, iBkg, sigmoidFlag, iNsigmaDown, iNsigmaUp);
            std::cout << "SigmoidCorrection = " << kBoolString[sigmoidFlag] << "; cutSettings = " << fullCutSettings << "..." << std::endl;
            outFile << "SigmoidCorrection = " << kBoolString[sigmoidFlag] << "; cutSettings = " << fullCutSettings << "..."
                    << "\n";

            if (analyse)
            {
              gSystem->Exec(Form("bash ~/Code/MasterThesis/Pion_PbPb/scripts/LaunchAnalysisSpec.sh %s 1 %d %d %s %f %f", fullCutSettings, iBkg, sigmoidFlag, spectraNameId, roi_n_sigma_down[iNsigmaDown], roi_n_sigma_up[iNsigmaUp]));
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