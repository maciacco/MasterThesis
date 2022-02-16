# LaunchAnalysis.sh
# This scripts launchs the whole analysis with tracking PID cut (select He3)

# parameters
cutDCAz=1.
cutTPCcls=89
binCountingFlag=1
expFlag=1 # 0->pol1, 1->expo, 2->pol2
sigmoidFlag=1
spectraHistNameId="1.0_89_0.1_1_1_1"
readTree=0
extractRatios=1

treeData="TreeOutData"
# treeMC="TreeOutMC_XSPlus"
treeMC="TreeOutMC"
signalName="SignalHe3"
spectraName="SpectraHe3"
# EfficiencyHe3="EfficiencyHe3_XSPlus"
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
ReadTreeData($argumentCuts,0.1,"$treeData","recreate","( ( (std::abs(pt)<2.5f) && (trackingPID==7) ) || !(std::abs(pt)<2.5f) )")
ReadTreeMC($argumentCuts,0.1,"$treeMC","recreate","( ( (std::abs(pt)<2.5f) && (trackingPID==7) ) || !(std::abs(pt)<2.5f) )")
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
.L AbsorptionError.cpp+
// SignalUnbinned($argumentCuts,$argumentSignal,"TreeOutData_NoPID","SignalHe3_NoPID","recreate",true,false)
// SignalUnbinned($argumentCuts,$argumentSignal,"TreeOutData_He3PID","SignalHe3_He3PID","recreate",true,false)
// SignalUnbinned($argumentCuts,$argumentSignal,"TreeOutData_AlphaPID","SignalHe3_AlphaPID","recreate",true,true)
SignalUnbinned($argumentCuts,0.1f,$argumentSignal,"$treeData","$signalName","recreate")
Efficiency($argumentCuts,0.1f,"$treeMC","$EfficiencyHe3")
EfficiencySec($argumentCuts,0.1f,"$treeMC","$EfficiencyHe3SecWD",0.3365047128558935)
Secondary($argumentCuts,0.1f,"$treeData","$treeMC","$EfficiencyHe3SecWD","$PrimaryHe3")
Spectra($argumentCuts,0.1f,$argumentSignal,$sigmoidFlag,"$spectraHistNameId","$spectraName","recreate","AnalysisResults","$signalName","$EfficiencyHe3","$PrimaryHe3")
// SignalLoss()
//AbsorptionError("AbsError","recreate","$spectraName")
.q
EOF
fi
