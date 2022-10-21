# Systematics.sh
# Launch all systematics analyses

# all variations
root -b -l << EOF
.x Systematics.cpp(kNPoints,true,true,true,true,"SystematicsAll")
EOF

# # only cut variation
# root -b -l << EOF
# .x Systematics.cpp(kNPoints,true,false,false,false,"SystematicsCutVar")
# EOF

# # only variation of signal extraction (bin counting vs. fit)
# root -b -l << EOF
# .x Systematics.cpp(kNPoints,false,true,false,false,"SystematicsBinCounting")
# EOF

# # only variation of background shape
# root -b -l << EOF
# .x Systematics.cpp(kNPoints,false,false,true,false,"SystematicsExp")
# EOF

# # only variation of primary fraction computation
# root -b -l << EOF
# .x Systematics.cpp(kNPoints,false,false,false,true,"SystematicsSigmoid")
# EOF

# systematics due to signal extraction
root -b -l << EOF
.x Systematics.cpp(kNPoints,true,false,true,false,"SystematicsPol1AndCutVar")
EOF

# systematics due to signal extraction
root -b -l << EOF
.x Systematics.cpp(kNPoints,true,false,false,false,"SystematicsExpoAndCutVar")
EOF
