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
double mismatch_n_sigma_down[] = {8.0,8.5,9.0};
double mismatch_n_sigma_up[] = {12.5,13.5,14.5};

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
  auto tmpNCutChi2TPC = kNCutChi2TPC;

  for (int iCutSettings = -1; iCutSettings < tmpNCutDCAz + tmpNTPCPidSigmas + tmpNCutTPCClusters + tmpNCutDCAxy + tmpNCutChi2TPC; ++iCutSettings)
  {
    char hname[100];

    for (int iBkg = 1; iBkg < 2; ++iBkg)
    {
      for(int iNsigmaUp = 0; iNsigmaUp < 6; ++iNsigmaUp) {
        for(int iNsigmaDown = 0; iNsigmaDown < 6; ++iNsigmaDown) {
          // bool binCountingFlag = 1 - iBin;
          int iNsigmaROIDown=1, iNsigmaROIUp=1, iNsigmaMismatchDown=1, iNsigmaMismatchUp=1;
          if (iNsigmaDown<3 && iNsigmaUp<3) {
            iNsigmaROIDown = iNsigmaDown;
            iNsigmaROIUp = iNsigmaUp;
            iNsigmaMismatchDown=1;
            iNsigmaMismatchUp=1;
          }
          else if (iNsigmaDown>2 && iNsigmaUp>2){
            iNsigmaMismatchDown = iNsigmaDown-3;
            iNsigmaMismatchUp = iNsigmaUp-3;
            iNsigmaROIUp=1;
            iNsigmaROIDown=1;
          }
          else continue;
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
            else if (iCutSettings >= (tmpNCutDCAz + tmpNTPCPidSigmas + tmpNCutTPCClusters) && iCutSettings < (tmpNCutDCAz + tmpNTPCPidSigmas + tmpNCutTPCClusters + tmpNCutDCAxy))
            {
              cutVariable = 3;
              cutIndex -= (tmpNCutDCAz + tmpNTPCPidSigmas + tmpNCutTPCClusters);
            }
            else if (iCutSettings >= (tmpNCutDCAz + tmpNTPCPidSigmas + tmpNCutTPCClusters + tmpNCutDCAxy))
            {
              cutVariable = 4;
              cutIndex -= (tmpNCutDCAz + tmpNTPCPidSigmas + tmpNCutTPCClusters + tmpNCutDCAxy);
            }
            fullCutSettings = Form("%s%d", cutSettings[cutVariable], cutIndex);
          }
          std::cout << "fullCutSettings = " << fullCutSettings << std::endl;

          double DCAxyCut = 0.12;
          if (cutVariable == 3) DCAxyCut = kCutDCAxyVariations[cutIndex];

          std::cout << "bkg selection = " << iBkg << "; roiNsigmaROIDown = " << roi_n_sigma_down[iNsigmaROIDown] << "; roiNsigmaROIUp = " << roi_n_sigma_up[iNsigmaROIUp] << "; roiNsigmaDownMismatch = " << mismatch_n_sigma_down[iNsigmaMismatchDown] << "; roiNsigmaMismatchUp = " << mismatch_n_sigma_up[iNsigmaMismatchUp] << "; dcaxycut = " << DCAxyCut << std::endl;
          if (analyse)
          {
            gSystem->Exec(Form("bash ~/Code/MuBFromRatios_PbPb/Pion_PbPb/scripts/LaunchAnalysisSignEffPrim.sh %s 1 1 %f %f %f %f %f", fullCutSettings, roi_n_sigma_down[iNsigmaROIDown], roi_n_sigma_up[iNsigmaROIUp], mismatch_n_sigma_down[iNsigmaMismatchDown], mismatch_n_sigma_up[iNsigmaMismatchUp], DCAxyCut));
          }

          for (int iSgm = 1; iSgm < 2; ++iSgm)
          {
            // bool binCountingFlag = 1 - iBin;
            bool sigmoidFlag = 1 - iSgm;
            auto tmpFullCutSettings = fullCutSettings;
            if (iCutSettings == -1) tmpFullCutSettings = Form("");
            auto spectraNameId = Form("%s_%d_%d_%d_%d_%d_%d",tmpFullCutSettings, iBkg, sigmoidFlag, iNsigmaROIDown, iNsigmaROIUp, iNsigmaMismatchDown, iNsigmaMismatchUp);
            std::cout << "SigmoidCorrection = " << kBoolString[sigmoidFlag] << "; cutSettings = " << fullCutSettings << "..." << std::endl;
            outFile << "SigmoidCorrection = " << kBoolString[sigmoidFlag] << "; cutSettings = " << fullCutSettings << "..."
                    << "\n";

            if (analyse)
            {
              gSystem->Exec(Form("bash ~/Code/MuBFromRatios_PbPb/Pion_PbPb/scripts/LaunchAnalysisSpec.sh %s 1 %d %d %s %f %f %f %f", fullCutSettings, iBkg, sigmoidFlag, spectraNameId, roi_n_sigma_down[iNsigmaROIDown], roi_n_sigma_up[iNsigmaROIUp], mismatch_n_sigma_down[iNsigmaMismatchDown], mismatch_n_sigma_up[iNsigmaMismatchUp]));
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