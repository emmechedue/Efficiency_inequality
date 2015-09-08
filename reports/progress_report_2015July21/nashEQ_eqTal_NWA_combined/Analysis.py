# Here I perform several analysis on the Inequality_efficieny results for ensemble

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

dummy = config["NE"]
NE=int(dummy)
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

#***************   LET'S START WITH ALL THE DIFFERENT ANALYSIS *******************************

#**************************** LOADING EVERYTHING *****************
# Loading 1
rankingmeth= str(1)

a=numpy.loadtxt("./cooperation"+rankingmeth+".txt")
cooperation1=a.transpose()

a=numpy.loadtxt("./efficiency"+rankingmeth+".txt")
efficiency1=a.transpose()

a=numpy.loadtxt("./giniCoef"+rankingmeth+".txt")
ginicoef1=a.transpose()

a=numpy.loadtxt("./wealth"+rankingmeth+".txt")
wealth1=a.transpose()

# Loading 2
rankingmeth= str(2)

a=numpy.loadtxt("./cooperation"+rankingmeth+".txt")
cooperation2=a.transpose()

a=numpy.loadtxt("./efficiency"+rankingmeth+".txt")
efficiency2=a.transpose()

a=numpy.loadtxt("./giniCoef"+rankingmeth+".txt")
ginicoef2=a.transpose()

a=numpy.loadtxt("./wealth"+rankingmeth+".txt")
wealth2=a.transpose()

# Loading 3
rankingmeth= str(3)

a=numpy.loadtxt("./cooperation"+rankingmeth+".txt")
cooperation3=a.transpose()

a=numpy.loadtxt("./efficiency"+rankingmeth+".txt")
efficiency3=a.transpose()

a=numpy.loadtxt("./giniCoef"+rankingmeth+".txt")
ginicoef3=a.transpose()

a=numpy.loadtxt("./wealth"+rankingmeth+".txt")
wealth3=a.transpose()

# Loading 4
rankingmeth= str(4)

a=numpy.loadtxt("./cooperation"+rankingmeth+".txt")
cooperation4=a.transpose()

a=numpy.loadtxt("./efficiency"+rankingmeth+".txt")
efficiency4=a.transpose()

a=numpy.loadtxt("./giniCoef"+rankingmeth+".txt")
ginicoef4=a.transpose()

a=numpy.loadtxt("./wealth"+rankingmeth+".txt")
wealth4=a.transpose()

#***********************************************************

#*************** Now let's compute all the averages ******************

TMAX = T+1

#Here I define all the arrays i will need
#1
COOP1=range(TMAX)
EFF1=range(TMAX)
GINI1=range(TMAX)
WEALTH1=range(TMAX)

COOPerrorstdplus1=range(TMAX)
EFFerrorstdplus1=range(TMAX)
GINIerrorstdplus1=range(TMAX)
WEALTHerrorstdplus1=range(TMAX)

COOPerrorstdminus1=range(TMAX)
EFFerrorstdminus1=range(TMAX)
GINIerrorstdminus1=range(TMAX)
WEALTHerrorstdminus1=range(TMAX)

#2
COOP2=range(TMAX)
EFF2=range(TMAX)
GINI2=range(TMAX)
WEALTH2=range(TMAX)

COOPerrorstdplus2=range(TMAX)
EFFerrorstdplus2=range(TMAX)
GINIerrorstdplus2=range(TMAX)
WEALTHerrorstdplus2=range(TMAX)

COOPerrorstdminus2=range(TMAX)
EFFerrorstdminus2=range(TMAX)
GINIerrorstdminus2=range(TMAX)
WEALTHerrorstdminus2=range(TMAX)

#3
COOP3=range(TMAX)
EFF3=range(TMAX)
GINI3=range(TMAX)
WEALTH3=range(TMAX)

COOPerrorstdplus3=range(TMAX)
EFFerrorstdplus3=range(TMAX)
GINIerrorstdplus3=range(TMAX)
WEALTHerrorstdplus3=range(TMAX)

COOPerrorstdminus3=range(TMAX)
EFFerrorstdminus3=range(TMAX)
GINIerrorstdminus3=range(TMAX)
WEALTHerrorstdminus3=range(TMAX)

#4
COOP4=range(TMAX)
EFF4=range(TMAX)
GINI4=range(TMAX)
WEALTH4=range(TMAX)

COOPerrorstdplus4=range(TMAX)
EFFerrorstdplus4=range(TMAX)
GINIerrorstdplus4=range(TMAX)
WEALTHerrorstdplus4=range(TMAX)

COOPerrorstdminus4=range(TMAX)
EFFerrorstdminus4=range(TMAX)
GINIerrorstdminus4=range(TMAX)
WEALTHerrorstdminus4=range(TMAX)

# For all
Root=math.sqrt(NE)

COOPstd=range(TMAX)
EFFstd=range(TMAX)
GINIstd=range(TMAX)
WEALTHstd=range(TMAX)

# Computing all the arrays
for i in range(TMAX):
	#1
	COOP1[i]=cooperation1[i].mean()
	EFF1[i]=efficiency1[i].mean()
	GINI1[i]=ginicoef1[i].mean()
	WEALTH1[i]=wealth1[i].mean() 
	#
	COOPstd[i]=cooperation1[i].std()  #Here I don't need to number them, because I just use them in the computation of the errors
	EFFstd[i]=efficiency1[i].std()
	GINIstd[i]=ginicoef1[i].std()
	WEALTHstd[i]=wealth1[i].std() 
	#
	COOPerrorstdplus1[i]=COOP1[i]+1.96*COOPstd[i]/Root
	COOPerrorstdminus1[i]=COOP1[i]-1.96*COOPstd[i]/Root
	EFFerrorstdplus1[i]=EFF1[i]+1.96*EFFstd[i]/Root
	EFFerrorstdminus1[i]=EFF1[i]-1.96*EFFstd[i]/Root
	GINIerrorstdplus1[i]=GINI1[i]+1.96*GINIstd[i]/Root
	GINIerrorstdminus1[i]=GINI1[i]-1.96*GINIstd[i]/Root
	WEALTHerrorstdplus1[i]=WEALTH1[i]+1.96*WEALTHstd[i]/Root
	WEALTHerrorstdminus1[i]=WEALTH1[i]-1.96*WEALTHstd[i]/Root
	#2
	COOP2[i]=cooperation2[i].mean()
	EFF2[i]=efficiency2[i].mean()
	GINI2[i]=ginicoef2[i].mean()
	WEALTH2[i]=wealth2[i].mean() 
	#
	COOPstd[i]=cooperation2[i].std()  #Here I don't need to number them, because I just use them in the computation of the errors
	EFFstd[i]=efficiency2[i].std()
	GINIstd[i]=ginicoef2[i].std()
	WEALTHstd[i]=wealth2[i].std() 
	#
	COOPerrorstdplus2[i]=COOP2[i]+1.96*COOPstd[i]/Root
	COOPerrorstdminus2[i]=COOP2[i]-1.96*COOPstd[i]/Root
	EFFerrorstdplus2[i]=EFF2[i]+1.96*EFFstd[i]/Root
	EFFerrorstdminus2[i]=EFF2[i]-1.96*EFFstd[i]/Root
	GINIerrorstdplus2[i]=GINI2[i]+1.96*GINIstd[i]/Root
	GINIerrorstdminus2[i]=GINI2[i]-1.96*GINIstd[i]/Root
	WEALTHerrorstdplus2[i]=WEALTH2[i]+1.96*WEALTHstd[i]/Root
	WEALTHerrorstdminus2[i]=WEALTH2[i]-1.96*WEALTHstd[i]/Root
	#3
	COOP3[i]=cooperation3[i].mean()
	EFF3[i]=efficiency3[i].mean()
	GINI3[i]=ginicoef3[i].mean()
	WEALTH3[i]=wealth3[i].mean() 
	#
	COOPstd[i]=cooperation3[i].std()  #Here I don't need to number them, because I just use them in the computation of the errors
	EFFstd[i]=efficiency3[i].std()
	GINIstd[i]=ginicoef3[i].std()
	WEALTHstd[i]=wealth3[i].std() 
	#
	COOPerrorstdplus3[i]=COOP3[i]+1.96*COOPstd[i]/Root
	COOPerrorstdminus3[i]=COOP3[i]-1.96*COOPstd[i]/Root
	EFFerrorstdplus3[i]=EFF3[i]+1.96*EFFstd[i]/Root
	EFFerrorstdminus3[i]=EFF3[i]-1.96*EFFstd[i]/Root
	GINIerrorstdplus3[i]=GINI3[i]+1.96*GINIstd[i]/Root
	GINIerrorstdminus3[i]=GINI3[i]-1.96*GINIstd[i]/Root
	WEALTHerrorstdplus3[i]=WEALTH3[i]+1.96*WEALTHstd[i]/Root
	WEALTHerrorstdminus3[i]=WEALTH3[i]-1.96*WEALTHstd[i]/Root
	#4
	COOP4[i]=cooperation4[i].mean()
	EFF4[i]=efficiency4[i].mean()
	GINI4[i]=ginicoef4[i].mean()
	WEALTH4[i]=wealth4[i].mean() 
	#
	COOPstd[i]=cooperation4[i].std()  #Here I don't need to number them, because I just use them in the computation of the errors
	EFFstd[i]=efficiency4[i].std()
	GINIstd[i]=ginicoef4[i].std()
	WEALTHstd[i]=wealth4[i].std() 
	#
	COOPerrorstdplus4[i]=COOP4[i]+1.96*COOPstd[i]/Root
	COOPerrorstdminus4[i]=COOP4[i]-1.96*COOPstd[i]/Root
	EFFerrorstdplus4[i]=EFF4[i]+1.96*EFFstd[i]/Root
	EFFerrorstdminus4[i]=EFF4[i]-1.96*EFFstd[i]/Root
	GINIerrorstdplus4[i]=GINI4[i]+1.96*GINIstd[i]/Root
	GINIerrorstdminus4[i]=GINI4[i]-1.96*GINIstd[i]/Root
	WEALTHerrorstdplus4[i]=WEALTH4[i]+1.96*WEALTHstd[i]/Root
	WEALTHerrorstdminus4[i]=WEALTH4[i]-1.96*WEALTHstd[i]/Root
	

# *************************** NOW WE ARE DONE AND WE CAN DO ALL THE PLOTTING !!! *********
legend1=r'$\alpha_i  w_i  r_i = O_i$'
legend2=r'$\alpha_i  w_i$'
legend3=r'$\alpha_i  r_i$'
legend4=r'$\alpha_i$'

time = range(TMAX)

# Plot for the cooperation
figure(num=None, figsize=(16, 12), dpi=160, facecolor='w', edgecolor='k')

plot(time,COOP1,'b',label=legend1)
plot(time,COOP2,'g',label=legend2)
plot(time,COOP3,'r',label=legend3)
plot(time,COOP4,'c',label=legend4)
#1
plot(time,COOPerrorstdplus1,'b--')
plot(time,COOPerrorstdminus1,'b--')
#2
plot(time,COOPerrorstdplus2,'g--')
plot(time,COOPerrorstdminus2,'g--')
#3
plot(time,COOPerrorstdplus3,'r--')
plot(time,COOPerrorstdminus3,'r--')
#4
plot(time,COOPerrorstdplus4,'c--')
plot(time,COOPerrorstdminus4,'c--')

axis([0, T, 0, 1])
title("Ensemble average of cooperation level vs. t")
ylabel("Average cooperation")
xlabel("t")
plt.legend((legend1,legend2,legend3,legend4),loc=1)
plt.tight_layout()
plt.savefig("cooperation.pdf",dpi=100)
plt.close()

# Plot for the efficiency
figure(num=None, figsize=(16, 12), dpi=160, facecolor='w', edgecolor='k')

plot(time,EFF1,'b',label=legend1)
plot(time,EFF2,'g',label=legend2)
plot(time,EFF3,'r',label=legend3)
plot(time,EFF4,'c',label=legend4)
#1
plot(time,EFFerrorstdplus1,'b--')
plot(time,EFFerrorstdminus1,'b--')
#2
plot(time,EFFerrorstdplus2,'g--')
plot(time,EFFerrorstdminus2,'g--')
#3
plot(time,EFFerrorstdplus3,'r--')
plot(time,EFFerrorstdminus3,'r--')
#4
plot(time,EFFerrorstdplus4,'c--')
plot(time,EFFerrorstdminus4,'c--')

title("Ensemble average of the efficiency vs. t")
ylabel("Average efficiency")
xlabel("t")
plt.legend((legend1,legend2,legend3,legend4),loc=1)
plt.tight_layout()
plt.savefig("efficiency.pdf",dpi=100)
plt.close()

# Plot for the gini coefficient
figure(num=None, figsize=(16, 12), dpi=160, facecolor='w', edgecolor='k')

plot(time,GINI1,'b',label=legend1)
plot(time,GINI2,'g',label=legend2)
plot(time,GINI3,'r',label=legend3)
plot(time,GINI4,'c',label=legend4)
#1
plot(time,GINIerrorstdplus1,'b--')
plot(time,GINIerrorstdminus1,'b--')
#2
plot(time,GINIerrorstdplus2,'g--')
plot(time,GINIerrorstdminus2,'g--')
#3
plot(time,GINIerrorstdplus3,'r--')
plot(time,GINIerrorstdminus3,'r--')
#4
plot(time,GINIerrorstdplus4,'c--')
plot(time,GINIerrorstdminus4,'c--')

title("Ensemble average of the gini coefficient vs. t")
ylabel("Average gini coefficient")
xlabel("t")
plt.legend((legend1,legend2,legend3,legend4),loc=1)
plt.tight_layout()
plt.savefig("gini.pdf",dpi=100)
plt.close()

# Plot for the wealth
figure(num=None, figsize=(16, 12), dpi=160, facecolor='w', edgecolor='k')
'''
plot(time,WEALTH1,'b',label=legend1)
plot(time,WEALTH2,'g',label=legend2)
plot(time,WEALTH3,'r',label=legend3)
plot(time,WEALTH4,'c',label=legend4)
#1
plot(time,WEALTHerrorstdplus1,'b--')
plot(time,WEALTHerrorstdminus1,'b--')
#2
plot(time,WEALTHerrorstdplus2,'g--')
plot(time,WEALTHerrorstdminus2,'g--')
#3
plot(time,WEALTHerrorstdplus3,'r--')
plot(time,WEALTHerrorstdminus3,'r--')
#4
plot(time,WEALTHerrorstdplus4,'c--')
plot(time,WEALTHerrorstdminus4,'c--')
'''
# Log scale
semilogy(time,WEALTH1,'b',label=legend1)
semilogy(time,WEALTH2,'g',label=legend2)
semilogy(time,WEALTH3,'r',label=legend3)
semilogy(time,WEALTH4,'c',label=legend4)
#1
semilogy(time,WEALTHerrorstdplus1,'b--')
semilogy(time,WEALTHerrorstdminus1,'b--')
#2
semilogy(time,WEALTHerrorstdplus2,'g--')
semilogy(time,WEALTHerrorstdminus2,'g--')
#3
semilogy(time,WEALTHerrorstdplus3,'r--')
semilogy(time,WEALTHerrorstdminus3,'r--')
#4
semilogy(time,WEALTHerrorstdplus4,'c--')
semilogy(time,WEALTHerrorstdminus4,'c--')

# Log scale
#plt.set_yscale('log', nonposy='clip', linthreshy=[YM,-YM] )

title("Ensemble average of the wealth vs. t")
ylabel("Average wealth")
xlabel("t")
plt.legend((legend1,legend2,legend3,legend4),loc=1)
plt.tight_layout()
plt.savefig("wealth.pdf",dpi=100)
plt.close()

