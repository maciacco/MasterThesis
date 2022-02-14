import numpy as np
import matplotlib.pyplot as plt
import pickle

bdt_score_dict = pickle.load(open("file_score_eff_dict","rb"))
bdt_score = bdt_score_dict["matter_0_5_8_14"]
bdt_eff = np.arange(0.1,1.0,0.01)
plt.scatter(bdt_score,bdt_eff,color="red")
bdt_score_dict_2 = pickle.load(open("file_score_eff_dict_2","rb"))
bdt_score_2 = bdt_score_dict_2["matter_0_5_8_14"]
bdt_eff_2 = np.arange(0.1,1.0,0.01)
print(f"eff_test = {bdt_score[64]}, eff = {bdt_eff[64]}, eff_train = {bdt_score_2[67]}, eff = {bdt_eff_2[67]}")
plt.scatter(bdt_score_2,bdt_eff_2,color="blue")
plt.xlabel("BDT output score")
plt.ylabel("BDT signal efficiency")
plt.title(r"$0-5\%, 8\leq c\mathrm{t} < 14\ \mathrm{cm}$")
plt.ylim([0.,1.1])
plt.grid()
plt.savefig("BDT_eff_score_0_5_8_ct_14_cm.pdf")