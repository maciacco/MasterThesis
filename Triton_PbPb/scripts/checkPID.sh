# checkPID.sh
# study trackingPID cut

readTreeDat=1
plotProjDat=1
readTreeMc=1
plotProjMc=1

if [ $readTreeDat -eq 1 ]; then
    root -b -l <<EOF
.L ReadTreeData.cpp+
ReadTreeData(1.,89,"TreeOutData_H3PID_all","recreate","(trackingPID==6)") // select H3
ReadTreeData(1.,89,"TreeOutData_He3PID_all","recreate","(trackingPID==7)") // select He3
ReadTreeData(1.,89,"TreeOutData_He4PID_all","recreate","(trackingPID==8)") // select He4
EOF
fi

if [ $plotProjDat -eq 1 ]; then
    root -b -l <<EOF
.L SignalUnbinned.cpp+
SignalUnbinned(1.,89,1,1,"TreeOutData_H3PID_all","Projection_H3PID_all","recreate",false) // projection H3
SignalUnbinned(1.,89,1,1,"TreeOutData_He3PID_all","Projection_He3PID_all","recreate",false) // projection He3
SignalUnbinned(1.,89,1,1,"TreeOutData_He4PID_all","Projection_He4PID_all","recreate",false) // projection He4
EOF
fi

if [ $readTreeMc -eq 1 ]; then
    root -b -l <<EOF
.L ReadTreeData.cpp+
ReadTreeData(1.,89,"TreeOutMc_H3PID_all","recreate","(trackingPID==6)",false) // select H3
ReadTreeData(1.,89,"TreeOutMc_He3PID_all","recreate","(trackingPID==7)",false) // select He3
ReadTreeData(1.,89,"TreeOutMc_He4PID_all","recreate","(trackingPID==8)",false) // select He4
EOF
fi

if [ $plotProjMc -eq 1 ]; then
    root -b -l <<EOF
.L SignalUnbinned.cpp+
SignalUnbinned(1.,89,1,1,"TreeOutMc_H3PID_all","Projection_Mc_H3PID_all","recreate",false) // projection H3
SignalUnbinned(1.,89,1,1,"TreeOutMc_He3PID_all","Projection_Mc_He3PID_all","recreate",false) // projection He3
SignalUnbinned(1.,89,1,1,"TreeOutMc_He4PID_all","Projection_Mc_He4PID_all","recreate",false) // projection He4
EOF
fi

rm results/TreeOutData_H3PID_all_anti.root
rm results/TreeOutData_He3PID_all_anti.root
rm results/TreeOutData_He4PID_all_anti.root
rm results/TreeOutMc_H3PID_all_anti.root
rm results/TreeOutMc_He3PID_all_anti.root
rm results/TreeOutMc_He4PID_all_anti.root

rm results/TreeOutData_H3PID_all_matt.root
rm results/TreeOutData_He3PID_all_matt.root
rm results/TreeOutData_He4PID_all_matt.root
rm results/TreeOutMc_H3PID_all_matt.root
rm results/TreeOutMc_He3PID_all_matt.root
rm results/TreeOutMc_He4PID_all_matt.root
