import numpy as np
import matplotlib.pyplot as plt


one_big_exp = np.loadtxt("../data/oneBigexperiment.txt")
radii = one_big_exp[0,:]
radii_square = one_big_exp[1,:]

fig,ax = plt.subplots(1,2,figsize=(10,10))

ax[0].hist(radii)
ax[0].set_title("Histogram of Radii")
ax[1].hist(radii_square)
ax[1].set_title("Histogram of Radius Squared")

plt.savefig("One_Big_Experiment")