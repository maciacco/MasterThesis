# LaunchAnalysis.sh

# parameters
cutSettings="$1"
binCountingFlag="$2"
expFlag="$3" # 1->sum of 2 exp
sigmoidFlag="$4"
spectraHistNameId="$5"
roiNsigmaDown="$6"
roiNsigmaUp="$7"
mismatchNsigmaDown="$8"
mismatchNsigmaUp="$9"
extractRatios=1

fileData="AnalysisResults"
fileMC="mc_20g7_20210929"
signalName="SignalPionSys"
spectraName="SpectraPionSys"
EfficiencyHe3="EfficiencyPionSys"
PrimaryHe3="PrimaryPionSys"

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
Spectra("$cutSettings",$roiNsigmaDown,$roiNsigmaUp,$mismatchNsigmaDown,$mismatchNsigmaUp,$argumentSignal,$sigmoidFlag,"$spectraHistNameId","$spectraName","update","AnalysisResults","$signalName","$EfficiencyHe3","$PrimaryHe3")
.q
EOF
fi
