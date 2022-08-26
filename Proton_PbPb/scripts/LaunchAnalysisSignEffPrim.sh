# LaunchAnalysis.sh

# parameters
cutSettings="$1"
binCountingFlag="$2"
expFlag="$3" # 1->sum of 2 exp
roiNsigma="$4"
dcaxycut="$5"
G3G4Prim="$6"
extractRatios=1

fileData="AnalysisResults"
fileDataEff="AnalysisResults_LHC21l5_full_largeDCA_cutChi2" # uncomment for high pt proton analysis (this is also the MCinj file in that case)
fileMCInj="AnalysisResults_LHC21l5_full_largeDCA_cutChi2" # uncomment for high pt proton analysis (this is also the MCinj file in that case)
#fileDataEff="AnalysisResults_LHC21l5_lowPtProton" # uncomment for low pt proton analysis
fileMC="mc"
fileMCInj="AnalysisResults_LHC21l5_LambdaCtMotherAndProtonsITSPID"
signalNameEff="SignalProtonSysEffTOF"
spectraNameEff="SpectraProtonSysEffTOF"
signalName="SignalProtonSysTOF"
spectraName="SpectraProtonSysTOF"
EfficiencyHe3="EfficiencyProtonSysTOF"
PrimaryHe3="PrimaryProtonSysTOF"
PrimaryHe3TPC="PrimaryProtonSysTOFTPC"

# create output directories
DIR_OUT=out
DIR_PLOT=plots

if [ -d "$DIR_OUT" ]; then
    echo "Directory '$DIR_OUT' already exists"
else
    mkdir $DIR_OUT
fi

if [ -d "$DIR_PLOT" ]; then
    echo "Directory '$DIR_PLOT' already exists"
else
    mkdir $DIR_PLOT
fi

# launch analysis
if (($cutSettings==999)); then
    cutSettings=""
    dcaxycut=0.12
fi

echo $cutSettings
argumentSignal="$binCountingFlag,$expFlag"

if [ $extractRatios -eq 1 ]; then
    root -b -l <<EOF
.L ../utils/RooGausExp.cxx+
.L ../utils/RooDSCBShape.cxx+
.L ../utils/RooGausDExp.cxx+
.L SignalBinned.cpp+
.L Efficiency.cpp+
.L Secondary.cpp+
.L Spectra.cpp+
SignalBinned("$cutSettings",$roiNsigma,$argumentSignal,"$fileData","$signalName","update")
Secondary("$cutSettings",$dcaxycut,"$fileData","$fileMC","$PrimaryHe3",$G3G4Prim)
.q
EOF
fi

root -b -l <<EOF
//.L SecondaryTPC.cpp+
//SecondaryTPC("$cutSettings",$dcaxycut,"$fileData","$fileMC","$PrimaryHe3TPC",$G3G4Prim)
EOF

if [ $extractRatios -eq 1 ]; then
    root -b -l <<EOF
.L ../utils/RooGausExp.cxx+
.L ../utils/RooDSCBShape.cxx+
.L ../utils/RooGausDExp.cxx+
.L SignalBinnedMC.cpp+
.L EfficiencyNew.cpp+
.L SecondaryMC.cpp+
SignalBinnedMC("$cutSettings",$roiNsigma,$argumentSignal,"$fileDataEff","$signalNameEff","recreate")
SecondaryMC("$cutSettings",$dcaxycut,"$fileDataEff","$fileMCInj","$PrimaryHe3Eff")
EfficiencyNew("$cutSettings","$fileMCInj","$EfficiencyHe3","$signalNameEff","$PrimaryHe3Eff")//,"update")
.q
EOF
fi