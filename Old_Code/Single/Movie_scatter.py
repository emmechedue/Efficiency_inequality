#*****************HERE I MAKE A MOVIE THAT SHOWS HOW THE WEALTH DISTRIBUTION EVOLVES DURING TIME ****************************

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

#********************** This is the nomenclature ********************
legend1=r'$\alpha_i  w_i  r_i = O_i$'
legend2=r'$\alpha_i  w_i$'
legend3=r'$\alpha_i  r_i$'
legend4=r'$\alpha_i$'

#*********************************************************

timepoints=T  # Here I compute how many steps I have


#************* NOW LET'S MAKE THE VIDEOS *********************

#*********************************** FIRST VIDEO *********************************************************************

#********************* As first thing, let's load everything: ********************************

t, c_average= np.loadtxt("./time1.txt", usecols=(0,3), unpack=True) #Here I load the time and the average cooperation
w_table=np.loadtxt("./wealth1.txt")
c_table=np.loadtxt("./cooperation1.txt")
r_table=np.loadtxt("./talent1.txt", usecols=(1,), unpack=True)
#**************************************************************************************


call(["mkdir", "videoscatter1"]) #I am creating a new directory to not create confusion
name = "./videoscatter1/scatter_1_"
filetype = '.png'

# Let's make all the pictures
for i in range(timepoints+1): #It is +1 because I have the 0th time step
	
	average = np.mean(w_table[i]) #This is the average wealth at time i
		
	# Coloring schema info
	color = r_table
	rmax = np.amax(r_table)
	rmin = np.amin(r_table)
	
	
	# The actual plot
	fig = plt.figure()
	ax = fig.add_subplot(111)
	lll=ax.scatter(c_table[i], w_table[i], c= color, s = 20,marker='o', vmin=rmin, vmax=rmax)
	plt.title("Scatter plot of wealth and cooperation in case "+legend1+" is rewarded")
	plt.ylabel('Wealth')
	plt.xlabel('Percentage of contribution')
	plt.plot(c_average[i],average, 'k_',ms=16,mew=3)
	plt.plot(c_average[i],average, 'k|',ms=16,mew=3)
	cbar=plt.colorbar(lll)
	cbar.set_label('Talent', rotation=270)
	#ax.set_ylim([0,1])
	ax.set_xlim([0,1])
	plt.tight_layout()
	
	# save plot as png file
	if i<10: ending = '00' + str(i)
	if 10 <= i and i<100: ending = '0' + str(i)
	if 100<= i and i<1000: ending = str(i)
	#if 1000<= i and i<10000: ending = str(i)
	filename = name + ending + filetype
	plt.savefig(filename)
	plt.close()
	print i
	

# Now we make the first movie
os.system("mencoder -really-quiet -mc 0 -noskip -skiplimit 0 -ovc lavc -lavcopts vcodec=mpeg4:vhq:trell:mbd=2:v4mv:vb_strategy=0:vlelim=0:vcelim=0:cmp=6:subcmp=6:precmp=6:predia=3:dia=3:vme=4:vqscale=1 \"mf://./videoscatter1/*.png\" -mf type=png:fps=5 -o scatter_1.avi")




#*********************************** SECOND VIDEO *********************************************************************

#********************* As first thing, let's load everything: ********************************

t, c_average= np.loadtxt("./time2.txt", usecols=(0,3), unpack=True) #Here I load the time and the average cooperation
w_table=np.loadtxt("./wealth2.txt")
c_table=np.loadtxt("./cooperation2.txt")
r_table=np.loadtxt("./talent2.txt", usecols=(1,), unpack=True)
#**************************************************************************************


call(["mkdir", "videoscatter2"]) #I am creating a new directory to not create confusion
name = "./videoscatter2/scatter_2_"
filetype = '.png'

# Let's make all the pictures
for i in range(timepoints+1): #It is +1 because I have the 0th time step
	
	average = np.mean(w_table[i]) #This is the average wealth at time i
		
	# Coloring schema info
	color = r_table
	rmax = np.amax(r_table)
	rmin = np.amin(r_table)
	
	
	# The actual plot
	fig = plt.figure()
	ax = fig.add_subplot(111)
	lll=ax.scatter(c_table[i], w_table[i], c= color, s = 20,marker='o', vmin=rmin, vmax=rmax)
	plt.title("Scatter plot of wealth and cooperation in case "+legend2+" is rewarded")
	plt.ylabel('Wealth')
	plt.xlabel('Percentage of contribution')
	plt.plot(c_average[i],average, 'k_',ms=16,mew=3)
	plt.plot(c_average[i],average, 'k|',ms=16,mew=3)
	cbar=plt.colorbar(lll)
	cbar.set_label('Talent', rotation=270)
	#ax.set_ylim([0,1])
	ax.set_xlim([0,1])
	plt.tight_layout()
	
	# save plot as png file
	if i<10: ending = '00' + str(i)
	if 10 <= i and i<100: ending = '0' + str(i)
	if 100<= i and i<1000: ending = str(i)
	#if 1000<= i and i<10000: ending = str(i)
	filename = name + ending + filetype
	plt.savefig(filename)
	plt.close()
	print i
	

# Now we make the first movie
os.system("mencoder -really-quiet -mc 0 -noskip -skiplimit 0 -ovc lavc -lavcopts vcodec=mpeg4:vhq:trell:mbd=2:v4mv:vb_strategy=0:vlelim=0:vcelim=0:cmp=6:subcmp=6:precmp=6:predia=3:dia=3:vme=4:vqscale=1 \"mf://./videoscatter2/*.png\" -mf type=png:fps=5 -o scatter_2.avi")



#*********************************** THIRD VIDEO *********************************************************************

#********************* As first thing, let's load everything: ********************************

t, c_average= np.loadtxt("./time3.txt", usecols=(0,3), unpack=True) #Here I load the time and the average cooperation
w_table=np.loadtxt("./wealth3.txt")
c_table=np.loadtxt("./cooperation3.txt")
r_table=np.loadtxt("./talent3.txt", usecols=(1,), unpack=True)
#**************************************************************************************


call(["mkdir", "videoscatter3"]) #I am creating a new directory to not create confusion
name = "./videoscatter3/scatter_3_"
filetype = '.png'

# Let's make all the pictures
for i in range(timepoints+1): #It is +1 because I have the 0th time step
	
	average = np.mean(w_table[i]) #This is the average wealth at time i
		
	# Coloring schema info
	color = r_table
	rmax = np.amax(r_table)
	rmin = np.amin(r_table)
	
	
	# The actual plot
	fig = plt.figure()
	ax = fig.add_subplot(111)
	lll=ax.scatter(c_table[i], w_table[i], c= color, s = 20,marker='o', vmin=rmin, vmax=rmax)
	plt.title("Scatter plot of wealth and cooperation in case "+legend3+" is rewarded")
	plt.ylabel('Wealth')
	plt.xlabel('Percentage of contribution')
	plt.plot(c_average[i],average, 'k_',ms=16,mew=3)
	plt.plot(c_average[i],average, 'k|',ms=16,mew=3)
	cbar=plt.colorbar(lll)
	cbar.set_label('Talent', rotation=270)
	#ax.set_ylim([0,1])
	ax.set_xlim([0,1])
	plt.tight_layout()	
	
	# save plot as png file
	if i<10: ending = '00' + str(i)
	if 10 <= i and i<100: ending = '0' + str(i)
	if 100<= i and i<1000: ending = str(i)
	#if 1000<= i and i<10000: ending = str(i)
	filename = name + ending + filetype
	plt.savefig(filename)
	plt.close()
	print i
	

# Now we make the first movie
os.system("mencoder -really-quiet -mc 0 -noskip -skiplimit 0 -ovc lavc -lavcopts vcodec=mpeg4:vhq:trell:mbd=2:v4mv:vb_strategy=0:vlelim=0:vcelim=0:cmp=6:subcmp=6:precmp=6:predia=3:dia=3:vme=4:vqscale=1 \"mf://./videoscatter3/*.png\" -mf type=png:fps=5 -o scatter_3.avi")





#*********************************** FOURTH VIDEO *********************************************************************

#********************* As first thing, let's load everything: ********************************

t, c_average= np.loadtxt("./time4.txt", usecols=(0,3), unpack=True) #Here I load the time and the average cooperation
w_table=np.loadtxt("./wealth4.txt")
c_table=np.loadtxt("./cooperation4.txt")
r_table=np.loadtxt("./talent4.txt", usecols=(1,), unpack=True)
#**************************************************************************************


call(["mkdir", "videoscatter4"]) #I am creating a new directory to not create confusion
name = "./videoscatter4/scatter_4_"
filetype = '.png'

# Let's make all the pictures
for i in range(timepoints+1): #It is +1 because I have the 0th time step
	
	average = np.mean(w_table[i]) #This is the average wealth at time i
		
	# Coloring schema info
	color = r_table
	rmax = np.amax(r_table)
	rmin = np.amin(r_table)
	
	
	# The actual plot
	fig = plt.figure()
	ax = fig.add_subplot(111)
	lll=ax.scatter(c_table[i], w_table[i], c= color, s = 20,marker='o', vmin=rmin, vmax=rmax)
	plt.title("Scatter plot of wealth and cooperation in case "+legend4+" is rewarded")
	plt.ylabel('Wealth')
	plt.xlabel('Percentage of contribution')
	plt.plot(c_average[i],average, 'k_',ms=16,mew=3)
	plt.plot(c_average[i],average, 'k|',ms=16,mew=3)
	cbar=plt.colorbar(lll)
	cbar.set_label('Talent', rotation=270)
	#ax.set_ylim([0,1])
	ax.set_xlim([0,1])
	plt.tight_layout()
	
	# save plot as png file
	if i<10: ending = '00' + str(i)
	if 10 <= i and i<100: ending = '0' + str(i)
	if 100<= i and i<1000: ending = str(i)
	#if 1000<= i and i<10000: ending = str(i)
	filename = name + ending + filetype
	plt.savefig(filename)
	plt.close()
	print i
	

# Now we make the first movie
os.system("mencoder -really-quiet -mc 0 -noskip -skiplimit 0 -ovc lavc -lavcopts vcodec=mpeg4:vhq:trell:mbd=2:v4mv:vb_strategy=0:vlelim=0:vcelim=0:cmp=6:subcmp=6:precmp=6:predia=3:dia=3:vme=4:vqscale=1 \"mf://./videoscatter4/*.png\" -mf type=png:fps=5 -o scatter_4.avi")









