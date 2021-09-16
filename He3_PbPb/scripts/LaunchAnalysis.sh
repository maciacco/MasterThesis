# LaunchAnalysis.sh
# This scripts launchs the whole analysis with tracking PID cut (select He3)

# parameters
cutDCAz=1.
cutTPCcls=89
binCountingFlag=1
expFlag=1 # 0->pol1, 1->expo, 2->pol2
sigmoidFlag=1
spectraHistNameId="1.0_89_1_1_1"
readTree=0
extractRatios=1

treeData="TreeOutData"
treeMC="TreeOutMC"
signalName="SignalHe3"
spectraName="SpectraHe3"
EfficiencyHe3="EfficiencyHe3"
EfficiencyHe3SecWD="EfficiencyHe3SecWD"
PrimaryHe3="PrimaryHe3"

# create output directories
DIR_OUT=out
DIR_RES=results
DIR_PLOT=plots

if [ -d "$DIR_OUT" ]; then
    echo "Directory '$DIR_OUT' already exists"
else
    mkdir $DIR_OUT
fi

if [ -d "$DIR_RES" ]; then
    echo "Directory '$DIR_RES' already exists"
else
    mkdir $DIR_RES
fi

if [ -d "$DIR_PLOT" ]; then
    echo "Directory '$DIR_PLOT' already exists"
else
    mkdir $DIR_PLOT
fi

# launch analysis

argumentCuts="$cutDCAz,$cutTPCcls"
argumentSignal="$binCountingFlag,$expFlag"

if [ $readTree -eq 1 ]; then
    root -b -l <<EOF
.L ReadTreeData.cpp+
.L ReadTreeMC.cpp+

// ReadTreeData($argumentCuts,"TreeOutData_NoPID","recreate","true") // no PID
// ReadTreeData($argumentCuts,"TreeOutData_He3PID","recreate","trackingPID==7") // kHe3
// ReadTreeData($argumentCuts,"TreeOutData_AlphaPID","recreate","trackingPID==8") // kHe4
ReadTreeData($argumentCuts,"$treeData","recreate","( ( (std::abs(pt)<2.5f) && (trackingPID==7) ) || !(std::abs(pt)<2.5f) )")
ReadTreeMC($argumentCuts,"$treeMC","recreate","( ( (std::abs(pt)<2.5f) && (trackingPID==7) ) || !(std::abs(pt)<2.5f) )")
EOF
fi

# delete unused files
rm results/TreeOutData_NoPID_anti.root
rm results/TreeOutData_NoPID_matt.root

rm results/TreeOutData_He3PID_anti.root
rm results/TreeOutData_He3PID_matt.root

rm results/TreeOutData_anti.root
rm results/TreeOutData_matt.root

rm results/TreeOutData_AlphaPID_anti.root
rm results/TreeOutData_AlphaPID_matt.root

if [ $extractRatios -eq 1 ]; then
    root -b -l <<EOF
.L SignalUnbinned.cpp+
.L Efficiency.cpp+
.L EfficiencySec.cpp+
.L Secondary.cpp+
.L Spectra.cpp+
.L SignalLoss.cpp+
// SignalUnbinned($argumentCuts,$argumentSignal,"TreeOutData_NoPID","SignalHe3_NoPID","recreate",true,false)
// SignalUnbinned($argumentCuts,$argumentSignal,"TreeOutData_He3PID","SignalHe3_He3PID","recreate",true,false)
// SignalUnbinned($argumentCuts,$argumentSignal,"TreeOutData_AlphaPID","SignalHe3_AlphaPID","recreate",true,true)
SignalUnbinned($argumentCuts,$argumentSignal,"$treeData","$signalName","recreate")
Efficiency($argumentCuts,"$treeMC","$EfficiencyHe3")
EfficiencySec($argumentCuts,"$treeMC","$EfficiencyHe3SecWD",$HYPER_TO_HE3_RATIO)
Secondary($argumentCuts,"$treeData","$treeMC","$EfficiencyHe3SecWD","$PrimaryHe3")
Spectra($argumentCuts,$argumentSignal,$sigmoidFlag,"$spectraHistNameId","$spectraName","recreate","AnalysisResults","$signalName","$EfficiencyHe3","$PrimaryHe3")
// SignalLoss()
.q
EOF
fi
