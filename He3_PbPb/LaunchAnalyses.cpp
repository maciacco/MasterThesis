// LaunchAnalyses.cpp
// This macro launches all the analyses needed for systematics computation

#include <stdlib.h>
#include <string.h>

#include <Riostream.h>
#include <TStopwatch.h>

#include "../utils/Config.h"

void LaunchAnalyses(const bool analyse = false)
{
  TStopwatch swatch;
  swatch.Start(true);
  std::ofstream outFile;
  outFile.open("ProcessedFiles.dat");

  for (int iDCA = 0; iDCA < kNCutDCAz; ++iDCA)
  {
    for (int iCls = 0; iCls < kNCutTPCClusters; ++iCls)
    {
      for (int iBin = 0; iBin < 2; ++iBin)
      {
        for (int iBkg = 0; iBkg < 2; ++iBkg)
        {
          for (int iSgm = 0; iSgm < 2; ++iSgm)
          {
            char hname[100];
            char par[100];

            bool binCountingFlag = 1 - iBin;
            bool expFlag = 1 - iBkg;
            bool sigmoidFlag = 1 - iSgm;
            if (kCutTPCClusters[iCls] < 100)
              snprintf(hname, 14, "%1.1f_0%d_%d_%d_%d", kCutDCAz[iDCA], kCutTPCClusters[iCls], binCountingFlag, expFlag, sigmoidFlag);
            else
              snprintf(hname, 14, "%1.1f_%d_%d_%d_%d", kCutDCAz[iDCA], kCutTPCClusters[iCls], binCountingFlag, expFlag, sigmoidFlag);
            sprintf(par, "bash ./scripts/AnalysisSys.sh %1.1f %d %s %s %s %s", kCutDCAz[iDCA], kCutTPCClusters[iCls], kBoolString[binCountingFlag], kBoolString[expFlag], kBoolString[sigmoidFlag], hname);

            std::cout << "|DCAz| < " << kCutDCAz[iDCA] << ", nClsTPC > " << kCutTPCClusters[iCls] << ", binCounting = " << kBoolString[binCountingFlag] << ", expBackground = " << kBoolString[expFlag] << ", sigmoidCorrection = " << kBoolString[sigmoidFlag] << "; processing treeHe3_" << hname << "..." << std::endl;
            outFile << "|DCAz| < " << kCutDCAz[iDCA] << ", nClsTPC > " << kCutTPCClusters[iCls] << ", binCounting = " << kBoolString[binCountingFlag] << ", expBackground = " << kBoolString[expFlag] << ", sigmoidCorrection = " << kBoolString[sigmoidFlag] << "; processing treeHe3_" << hname << "..."
                    << "\n";

            if (analyse)
              system(par);
          }
        }
      }
    }
  }

  swatch.Stop();
  std::cout << "Elapsed (real) time = " << swatch.RealTime() << " s" << std::endl;

  outFile.close();
}