# LaunchAnalysis.sh

# parameters
cutSettings=""
binCountingFlag=1
expFlag=1 # 1->sum of 2 exp, 0 -> sum of exp and pol
sigmoidFlag=0
spectraHistNameId=""
extractRatios=1

#fileData="AnalysisResults_LHC22b9_lowPT_MB"
#fileMC="AnalysisResults_LHC22b9_lowPT_MB"
#fileData="AnalysisResults_LHC21l5_full_largeDCA_cutChi2"
#fileMC="AnalysisResults_LHC21l5_full_largeDCA_cutChi2"
#fileData="AnalysisResults_LHC21l5_lowPtProton_ITSrecalibrated"
#fileMC="AnalysisResults_LHC21l5_lowPtProton_ITSrecalibrated"
fileData="../AnalysisResults_LHC22b9_3"
fileMC="../AnalysisResults_LHC22b9_3"
#signalName="SignalProtonMC_22b9_lowPT_MB"
#spectraName="SpectraProtonMCMC_22b9_lowPT_MB"
#EfficiencyHe3="EfficiencyProton_lowPt_LHC22b9_3"
#PrimaryHe3="PrimaryProtonMC_22b9_lowPT_MB"
signalName="SignalProtonTPC_MC21l5_raw_primaryTPC_try"
spectraName="SpectraProtonTPC_MC21l5_raw_primaryTPC_try"
#EfficiencyHe3="EfficiencyProtonMC_21l5_LOWPT_TESTTPC" #"EfficiencyProtonMC_21l5_LOWPT_TESTTPC_try"
EfficiencyHe3="EfficiencyProton_LHC22b9_3"
PrimaryHe3="PrimaryProtonTPC_MC21l5_raw_primaryTPC_try"

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

echo $cutSettings
argumentSignal="$binCountingFlag,$expFlag"

if [ $extractRatios -eq 1 ]; then
    root -b -l <<EOF
.L ../utils/RooGausExp.cxx+
.L ../utils/RooDSCBShape.cxx+
.L ../utils/RooGausDExp.cxx+
.L SignalBinnedMC.cpp+
.L SecondaryMC.cpp+
.L EfficiencyNew.cpp+
.L Efficiency.cpp+
//.L Spectra.cpp+
SignalBinnedMC("$cutSettings",8.,$argumentSignal,"$fileData","$signalName","recreate")
SecondaryMC("$cutSettings",0.12,"$fileData","$fileMC","$PrimaryHe3")
EfficiencyNew("$cutSettings","$fileMC","$EfficiencyHe3","$signalName","$PrimaryHe3")
//Efficiency("$cutSettings","$fileMC","$EfficiencyHe3")
.q
EOF
fi
