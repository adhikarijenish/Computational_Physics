import numpy as np
import matplotlib.pyplot as plt


one_big_exp = np.loadtxt("../data/oneBigexperiment.txt")
radii = one_big_exp[:,0]
radii_square = one_big_exp[:,1]

fig,ax = plt.subplots(1,2,figsize=(10,10))

ax[0].hist(radii)
ax[0].set_xlabel(r"Radius")
ax[0].set_ylabel(r"#",loc='top')
ax[0].set_title("Histogram of Radii")
ax[1].hist(radii_square)
ax[1].set_ylabel(r"#",loc='top')
ax[1].set_xlabel(r"Radius Squared")
ax[1].set_title("Histogram of Radius Squared")

plt.savefig("One_Big_Experiment_Radii_histogram")
plt.close()

indicator = 4*(radii_square <=1)
mean = np.mean(indicator)
std = np.std(indicator)
n, bins, patches = plt.hist(indicator,label='Histgrom of indicator')
plt.vlines(mean,0,n.max(),colors='r',label='Mean of samples')
plt.vlines(np.pi,0,n.max(),colors='b',label=r'$\pi$')
plt.fill_between([mean-std,mean+std],0,8000,alpha=0.2,label='Standard Deviation')
plt.xlabel(r'Indicatior variable')
plt.ylabel(r'#',loc='top')
plt.legend()
plt.savefig('Pi_comparison')
