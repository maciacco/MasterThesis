# LaunchAnalysis.sh

# parameters
cutSettings="$1"
binCountingFlag="$2"
expFlag="$3" # 1->sum of 2 exp
roiNsigmaMin="$4"
roiNsigmaMax="$5"
mismatchNsigmaMin="$6"
mismatchNsigmaMax="$7"
dcaxycut="$8"
extractRatios=1

fileData="AnalysisResults"
fileDataEff="AnalysisResults_LHC21l5_full_largeDCA_cutChi2"
fileMC="mc"
signalNameEff="SignalPionSysEff_extend3"
spectraNameEff="SpectraPionSysEff_extend3"
signalName="SignalPionSys_extend3"
spectraName="SpectraPionSys_extend3"
EfficiencyHe3="EfficiencyPionSys_extend3"
PrimaryHe3="PrimaryPionSys_extend3"

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
SignalBinned("$cutSettings",$roiNsigmaMin,$roiNsigmaMax,$mismatchNsigmaMin,$mismatchNsigmaMax,$argumentSignal,"$fileData","$signalName","update")
Secondary("$cutSettings",$dcaxycut,"$fileData","$fileMC","$PrimaryHe3")//,true) // uncomment to use roofit
.q
EOF
fi

if [ $extractRatios -eq 1 ]; then
    root -b -l <<EOF
.L ../utils/RooGausExp.cxx+
.L ../utils/RooDSCBShape.cxx+
.L ../utils/RooGausDExp.cxx+
.L SignalBinnedMCv2.cpp+
.L EfficiencyNew.cpp+
.L SecondaryMC.cpp+
SignalBinnedMCv2("$cutSettings",$roiNsigmaMin,$roiNsigmaMax,$mismatchNsigmaMin,$mismatchNsigmaMax,$argumentSignal,"$fileDataEff","$signalNameEff","recreate")
SecondaryMC("$cutSettings",$dcaxycut,"$fileDataEff","$fileDataEff","$PrimaryHe3Eff")
EfficiencyNew("$cutSettings","$fileDataEff","$EfficiencyHe3","$signalNameEff","$PrimaryHe3Eff")//,"update")
.q
EOF
fi