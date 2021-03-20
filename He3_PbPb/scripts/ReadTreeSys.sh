# ReadTreeSys.sh
# Launch ReadTree macros for systematic uncertainties computation (cut variation)

cutDCAz="$1"
cutTPCcls="$2"

# create output directories
DIR_OUT=out
DIR_RES=results

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

# launch analysis

argumentTree="$cutDCAz,$cutTPCcls"
argumentSignal="$binCountingFlag,$expFlag"

root -b -l <<EOF
.L ReadTreeData.cpp+
.L ReadTreeMC.cpp+
ReadTreeData($argumentTree,"TreeOutDataSys","update")
ReadTreeMC($argumentTree,"TreeOutMCSys","update")
.q
EOF

# delete unused files
rm results/TreeOutDataSys_anti.root results/TreeOutDataSys_matt.root
