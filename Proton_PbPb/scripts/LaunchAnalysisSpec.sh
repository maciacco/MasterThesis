# LaunchAnalysis.sh

# parameters
cutSettings="$1"
binCountingFlag="$2"
expFlag="$3" # 1->sum of 2 exp
sigmoidFlag="$4"
spectraHistNameId="$5"
roiNsigma="$6"
G3G4Prim="$7"
extractRatios=1

fileData="AnalysisResults"
fileMC="mc_20g7_20210929"
signalName="SignalProtonSysTOF_extend_4"
spectraName="SpectraProtonSysTOF_extend_4"
spectraNameTPC="SpectraProtonSysTOF_extend_4TPC"
EfficiencyHe3="EfficiencyProtonSysTOF_extend_4"
PrimaryHe3="PrimaryProtonSysTOF_extend_4"
PrimaryHe3TPC="PrimaryProtonSysTOF_extend_4TPC"

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
fi

echo $cutSettings
argumentSignal="$binCountingFlag,$expFlag"

if [ $extractRatios -eq 1 ]; then
    root -b -l <<EOF
.L Spectra.cpp+
Spectra("$cutSettings",$roiNsigma,$G3G4Prim,$argumentSignal,$sigmoidFlag,"$spectraHistNameId","$spectraName","update","AnalysisResults","$signalName","$EfficiencyHe3","$PrimaryHe3")
.q
EOF
fi

if [ $extractRatios -eq 1 ]; then
    root -b -l <<EOF
//.L SpectraTPC.cpp+
//SpectraTPC("$cutSettings",$roiNsigma,$G3G4Prim,$argumentSignal,$sigmoidFlag,"$spectraHistNameId","$spectraNameTPC","update","AnalysisResults","$signalName","$EfficiencyHe3","$PrimaryHe3TPC")
.q
EOF
fi
