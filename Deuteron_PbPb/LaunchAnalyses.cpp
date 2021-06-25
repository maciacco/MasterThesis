// LaunchAnalyses.cpp
// This macro launches all the analyses needed for systematics computation

#include <stdlib.h>
#include <string.h>

#include <Riostream.h>
#include <TStopwatch.h>
#include <TROOT.h>
#include <TSystem.h>

#include "../utils/Config.h"

using namespace deuteron;

void LaunchAnalyses(const bool analyse = false)
{
  TStopwatch swatch;
  swatch.Start(true);
  std::ofstream outFile;
  outFile.open("ProcessedFiles.dat");

  for (int iCutSettings = 0; iCutSettings < kNCutDCAz + kNTPCPidSigmas + kNCutTPCClusters; ++iCutSettings)
  {
    char hname[100];

    // bool binCountingFlag = 1 - iBin;
    int cutVariable = 0;
    int cutIndex = iCutSettings;
    if (iCutSettings >= kNCutDCAz && iCutSettings < (kNCutDCAz + kNTPCPidSigmas))
    {
      cutVariable = 1;
      cutIndex -= kNCutDCAz;
    }
    else if (iCutSettings >= (kNCutDCAz + kNTPCPidSigmas))
    {
      cutVariable = 2;
      cutIndex -= (kNCutDCAz + kNTPCPidSigmas);
    }
    auto fullCutSettings = Form("%s%d", cutSettings[cutVariable], cutIndex);
    std::cout << "fullCutSettings = " << fullCutSettings << std::endl;

    for (int iBkg = 0; iBkg < 2; ++iBkg)
    {
      if (analyse)
      {
        gSystem->Exec(Form("bash ~/Code/MasterThesis/Deuteron_PbPb/scripts/LaunchAnalysisSignEffPrim.sh %s 0 %d", fullCutSettings, iBkg));
      }

      for (int iSgm = 0; iSgm < 2; ++iSgm)
      {
        // bool binCountingFlag = 1 - iBin;
        bool sigmoidFlag = 1 - iSgm;
        auto spectraNameId = Form("%s_%d_%d",fullCutSettings, iBkg, sigmoidFlag);
        std::cout << "SigmoidCorrection = " << kBoolString[sigmoidFlag] << "; cutSettings = " << fullCutSettings << "..." << std::endl;
        outFile << "SigmoidCorrection = " << kBoolString[sigmoidFlag] << "; cutSettings = " << fullCutSettings << "..."
                << "\n";

        if (analyse)
        {
          gSystem->Exec(Form("bash ~/Code/MasterThesis/Deuteron_PbPb/scripts/LaunchAnalysisSpec.sh %s 0 %d %d %s", fullCutSettings, iBkg, sigmoidFlag, spectraNameId));
        }
      }
    }
  }

  swatch.Stop();
  std::cout << "Elapsed (real) time = " << swatch.RealTime() << " s" << std::endl;

  outFile.close();
}