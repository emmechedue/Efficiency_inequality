# Here I perform several analysis on the Inequality_efficieny results

import numpy
import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.colors as mc
import matplotlib.mlab as mlab
from matplotlib.ticker import NullFormatter
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
import matplotlib.path as mpath
import matplotlib.lines as mlines
from pylab import *
from numpy import array
from configobj import ConfigObj
from subprocess import call
import os

#****************************************Let's first of all get all the parameters:     **********************
config = ConfigObj("./parameters.txt")

dummy = config["N"]
N=int(dummy)
dummy = config["T"]
T=int(dummy)
dummy = config["S"]
S=int(dummy)
dummy = config["Q"]
Q=float(dummy)
dummy = config["mu"]
mu=float(dummy)
dummy = config["sigmag"]
sigmag=float(dummy)
dummy = config["W0"]
W0=float(dummy)
dummy = config["Choice"]
cho=int(dummy)
dummy = config["monincr"]
monincr=int(dummy)
dummy = config["beta"]
beta=double(dummy)

#***************************************************************************************

#**************************** Now I load everything: *****************

a=numpy.loadtxt("./time1.txt")
data1=a.transpose()

a=numpy.loadtxt("./time2.txt")
data2=a.transpose()

a=numpy.loadtxt("./time3.txt")
data3=a.transpose()

a=numpy.loadtxt("./time4.txt")
data4=a.transpose()

#******************************************************

#*************** Now let's make the plot for the Gini coefficient ************************
plt.plot(data1[0], data1[1])
plt.plot(data1[0], data2[1])
plt.plot(data1[0], data3[1])
plt.plot(data1[0], data4[1])
plt.title("Gini coefficient for different ranking systems")
legend1=r'$\alpha_i  w_i  r_i = O_i$'
legend2=r'$\alpha_i  w_i$'
legend3=r'$\alpha_i  r_i$'
legend4=r'$\alpha_i$'
plt.legend((legend1,legend2,legend3,legend4),loc=1)
plt.xlabel('time')
plt.ylabel('Gini coefficient')
plt.tight_layout()
plt.savefig('gini.pdf',dpi=100)
close()
#******************************************************************************

#********** Now let's make the plot for the Efficiency, 

# First we have to define how many time points I have:
timepoints=T + 1


# Now let's plot the relative growth for each of the 4 rules
plt.plot(data1[0], data1[2])
plt.plot(data1[0], data2[2])
plt.plot(data1[0], data3[2])
plt.plot(data1[0], data4[2])
plt.title("Total wealth (?) for different ranking systems")
legend1=r'$\alpha_i  w_i  r_i = O_i$'
legend2=r'$\alpha_i  w_i$'
legend3=r'$\alpha_i  r_i$'
legend4=r'$\alpha_i$'
plt.legend((legend1,legend2,legend3,legend4),loc=4)
plt.xlabel('time')
#plt.ylabel(r'Relative Growth, $\left(w[t]-w[t-1]\right)/w[t-1]$')
plt.ylabel(r'Total wealth')
plt.tight_layout()
plt.savefig('growth.pdf',dpi=100)
close()

#***********************************************************************************************

#*********** Now let's plot the average cooperation level for each time step.
plt.plot(data1[0], data1[3])
plt.plot(data1[0], data2[3])
plt.plot(data1[0], data3[3])
plt.plot(data1[0], data4[3])
plt.title("Average cooperation level for different ranking systems")
legend1=r'$\alpha_i  w_i  r_i = O_i$'
legend2=r'$\alpha_i  w_i$'
legend3=r'$\alpha_i  r_i$'
legend4=r'$\alpha_i$'
plt.legend((legend1,legend2,legend3,legend4),loc=1)
plt.xlabel('time')
plt.ylabel('Average cooperation level')
plt.tight_layout()
plt.savefig('cooperation.pdf',dpi=100)
close()
