import numpy as np
import matplotlib.pyplot as plt
import pickle

bdt_score_dict = pickle.load(open("./file_score_eff_dict","rb"))
bdt_score = bdt_score_dict["all_0_5_2_3"]
bdt_eff = np.arange(0.1,1.0,0.01)
plt.scatter(bdt_score,bdt_eff,color="red")
plt.xlabel("BDT output score")
plt.ylabel("BDT signal efficiency")
plt.title(r"$0-5\%, 2\leq c\mathrm{t} < 3\ \mathrm{cm}$")
plt.ylim([0.,1.1])
plt.grid()
plt.savefig("BDT_eff_score_0_5_2_ct_3_cm.pdf")