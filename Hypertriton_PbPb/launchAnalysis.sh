# preselection efficiency
#python3 plot_efficiencies_correlations.py -m -s -eff

# plot feature distributions and correlations
#python3 plot_efficiencies_correlations.py -s -m

# train/test models and compute scores
#python3 ml_analysis.py -m -s -t -c

# train/test plots
#python3 ml_analysis.py -m -c

# application
#python3 ml_analysis.py -s -m -a

# significance scan
python3 significance_scan.py

# signal extraction (pol1)
# python3 signal_extraction.py

# signal extraction (expo)
# python3 signal_extraction.py -b

# absorption correction
# python3 he3_absorption_analysis.py

# lifetime fits + ratio
#python3 ratio.py

# systematics (ratio)
#python3 systematics_ratio.py

# systematics (cpt)
#python3 systematics_cpt.py
