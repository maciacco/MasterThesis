// LaunchTrees.cpp
// This macro launches all the trees needed for systematics computation

#include <stdlib.h>
#include <string.h>

#include <Riostream.h>
#include <TStopwatch.h>

#include "../Config.h"

using namespace triton;

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
        for (int iChi2TPC = 0; iChi2TPC < kNCutChi2TPC; ++iChi2TPC)
        {
          char hname[100];
          char par[100];

          if (kCutTPCClusters[iCls] < 100)
            snprintf(hname, 17, "%1.1f_0%d_%1.2f_%1.2f", kCutDCAz[iDCA], kCutTPCClusters[iCls], kCutDCAxy[iDCAxy], kCutChi2TPC[iChi2TPC]);
          else
            snprintf(hname, 17, "%1.1f_%d_%1.2f_%1.2f", kCutDCAz[iDCA], kCutTPCClusters[iCls], kCutDCAxy[iDCAxy], kCutChi2TPC[iChi2TPC]);
          sprintf(par, "bash ./scripts/ReadTreeSys.sh %1.1f %d %1.2f %1.2f %s", kCutDCAz[iDCA], kCutTPCClusters[iCls], kCutDCAxy[iDCAxy], kCutChi2TPC[iChi2TPC], hname);

          std::cout << "|DCAz| < " << kCutDCAz[iDCA] << ", nClsTPC > " << kCutTPCClusters[iCls] << ", |DCAxy| < " << kCutDCAxy[iDCAxy] << ", chi2TPC < " << kCutChi2TPC[iChi2TPC] << "; processing treeHe3_" << hname << "..." << std::endl;
          outFile << "|DCAz| < " << kCutDCAz[iDCA] << ", nClsTPC > " << kCutTPCClusters[iCls] << ", |DCAxy| < " << kCutDCAxy[iDCAxy] << ", chi2TPC < " << kCutChi2TPC[iChi2TPC] << "; processing treeHe3_" << hname << "..."
                  << "\n";

          if (analyse)
            system(par);
        }
      }
    }
  }

  swatch.Stop();
  std::cout << "Elapsed (real) time = " << swatch.RealTime() << " s" << std::endl;

  outFile.close();
}