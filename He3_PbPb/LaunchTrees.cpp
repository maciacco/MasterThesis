// LaunchTrees.cpp
// This macro launches all the trees needed for systematics computation

#include <stdlib.h>
#include <string.h>

#include <Riostream.h>
#include <TStopwatch.h>

#include "../utils/Config.h"

using namespace he3;

void LaunchTrees(const bool analyse = false)
{
  TStopwatch swatch;
  swatch.Start(true);
  std::ofstream outFile;
  outFile.open("ProcessedTrees.dat");

  for (int iDCA = 0; iDCA < kNCutDCAz; ++iDCA)
  {
    for (int iCls = 0; iCls < kNCutTPCClusters; ++iCls)
    {
      for (int iDCAxy = 0; iDCAxy < kNCutDCAxy; ++iDCAxy)
      {
        char hname[100];
        char par[100];

        if (kCutTPCClusters[iCls] < 100)
          snprintf(hname, 12, "%1.1f_0%d_%1.2f", kCutDCAz[iDCA], kCutTPCClusters[iCls], kCutDCAxy[iDCAxy]);
        else
          snprintf(hname, 12, "%1.1f_%d_%1.2f", kCutDCAz[iDCA], kCutTPCClusters[iCls], kCutDCAxy[iDCAxy]);
        sprintf(par, "bash ./scripts/ReadTreeSys.sh %1.1f %d %1.2f %s", kCutDCAz[iDCA], kCutTPCClusters[iCls], kCutDCAxy[iDCAxy], hname);

        std::cout << "|DCAz| < " << kCutDCAz[iDCA] << ", nClsTPC > " << kCutTPCClusters[iCls] << ", |DCAxy| < " << kCutDCAxy[iDCAxy] << "; processing treeHe3_" << hname << "..." << std::endl;
        outFile << "|DCAz| < " << kCutDCAz[iDCA] << ", nClsTPC > " << kCutTPCClusters[iCls] << ", |DCAxy| < " << kCutDCAxy[iDCAxy] << "; processing treeHe3_" << hname << "..."
                << "\n";

        if (analyse)
          system(par);
      }
    }
  }

  swatch.Stop();
  std::cout << "Elapsed (real) time = " << swatch.RealTime() << " s" << std::endl;

  outFile.close();
}