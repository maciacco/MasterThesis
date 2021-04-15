# preselection efficiency
python3 ml_analysis.py -m -s -eff

# plot feature distributions and correlations
# python3 ml_analysis.py -s -m

# train/test models and compute scores
python3 ml_analysis.py -m -s -t -c

# application
python3 ml_analysis.py -s -m -a

# significance scan
python3 significance_scan.py

# signal extraction
python3 signal_extraction.py

# lifetime fits + ratio
python3 ratio.py
