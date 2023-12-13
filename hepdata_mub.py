import ROOT
import numpy as np



with open("MuB.yaml", "w") as ofile:
    ofile.write("independent_variables:\n")
    ofile.write("- header: {name: \'$\sqrt{s_{\\rm{NN}}}$ (TeV)\'}\n")
    ofile.write("  values:\n")
    ofile.write("  - {value: 5.020}\n")

    ofile.write("dependent_variables:\n")
    ofile.write("- header: {name: '$\mu_B$ (MeV)'}\n")
    ofile.write("  values:\n")
    ofile.write("    - errors:\n")
    ofile.write("      - {label: Total uncertainty, symerror: 0.45}\n")
    ofile.write("      value: 0.71")