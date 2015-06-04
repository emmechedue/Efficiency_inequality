#*****************HERE I MAKE A MOVIE THAT SHOWS HOW THE WEALTH DISTRIBUTION EVOLVES DURING TIME ****************************

import matplotlib
import numpy
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
dummy = config["interval"]
interval = float(dummy)
dummy = config["S"]
S=int(dummy)
dummy = config["Q"]
Q=float(dummy)
dummy = config["lambda"]
lam=float(dummy)
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

data1=numpy.loadtxt("./wealth1.txt")

data2=numpy.loadtxt("./wealth2.txt")

data3=numpy.loadtxt("./wealth3.txt")

data4=numpy.loadtxt("./wealth4.txt")

#******************************************************

#********************** This is the nomenclature ********************
legend1=r'$\alpha_i  w_i  r_i = O_i$'
legend2=r'$\alpha_i  w_i$'
legend3=r'$\alpha_i  r_i$'
legend4=r'$\alpha_i$'

#*********************************************************

#************************* NOW LET'S MAKE THE VIDEOS ********************

#Let's compute the number of steps:
timepoints=int(T/interval)

#This is in how many intervals I want to divide my histograms
KK=100

#let's prepare the x axis	
raggio=numpy.zeros(KK)
for j in range(KK):
	raggio[j]=(1./KK) * j


#************** Let's make the first video ****************************
call(["mkdir", "video1"]) #I am creating a new directory to not create confusion
name = "./video1/wealth_1_"
filetype = '.png'


# Let's make all the pictures
for i in range(timepoints+1):
	
	#Let's compute the histogram
	fuffo=zeros(KK) #Here I just fill the histofram with zeros
	test=data1[i]
	
	#Here I count how many hits per bin
	
	for j in range(N): 
		dummy1=test[j]*KK
		dummy2=floor(dummy1)
		fuffo[dummy2]=fuffo[dummy2] + 1
    
   	
   	# The actual plot
	plt.plot(raggio,fuffo)
	plt.axis([0, 1, 0, N])
	plt.title("Wealth distribution in case "+legend1+" is rewarded")
	plt.xlabel('normalized wealth w')
	plt.ylabel('Amount of people with wealth w')
	
	# The clock, this is an inset axes over the main axes where I plot the clock
	time=double(i)/timepoints #Here I compute the fraction of the pie that has to be filled, i.e. the percentage of the time
	a = axes([.65, .6, .2, .2], axisbg='w',aspect=1.0)
	fraction=[1-time,time]
	mycolors=["white","black"]
	plt.pie(fraction,colors=mycolors,startangle=90)

	
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
os.system("mencoder -really-quiet -mc 0 -noskip -skiplimit 0 -ovc lavc -lavcopts vcodec=mpeg4:vhq:trell:mbd=2:v4mv:vb_strategy=0:vlelim=0:vcelim=0:cmp=6:subcmp=6:precmp=6:predia=3:dia=3:vme=4:vqscale=1 \"mf://./video1/*.png\" -mf type=png:fps=5 -o video_1.avi")


#************** Let's make the second video ****************************
call(["mkdir", "video2"]) #I am creating a new directory to not create confusion
name = "./video2/wealth_2_"
filetype = '.png'

# Let's make all the pictures
for i in range(timepoints+1):

	#Let's compute the histogram
	fuffo=zeros(KK) #Here I just fill the histofram with zeros
	test=data2[i]
	
	#Here I count how many hits per bin
	
	for j in range(N): 
		dummy1=test[j]*KK
		dummy2=floor(dummy1)
		fuffo[dummy2]=fuffo[dummy2] + 1
    
   	
   	# The actual plot
	plt.plot(raggio,fuffo)
	plt.axis([0, 1, 0, N])
	plt.title("Wealth distribution in case "+legend2+" is rewarded")
	plt.xlabel('normalized wealth w')
	plt.ylabel('Amount of people with wealth w')
	
	# The clock, this is an inset axes over the main axes where I plot the clock
	time=double(i)/timepoints #Here I compute the fraction of the pie that has to be filled, i.e. the percentage of the time
	a = axes([.65, .6, .2, .2], axisbg='w',aspect=1.0)
	fraction=[1-time,time]
	mycolors=["white","black"]
	plt.pie(fraction,colors=mycolors,startangle=90)

	
	# save plot as png file
	if i<10: ending = '00' + str(i)
	if 10 <= i and i<100: ending = '0' + str(i)
	if 100<= i and i<1000: ending = str(i)
	#if 1000<= i and i<10000: ending = str(i)
	filename = name + ending + filetype
	plt.savefig(filename)
	plt.close()
	print i
	

# Now we make the second movie
os.system("mencoder -really-quiet -mc 0 -noskip -skiplimit 0 -ovc lavc -lavcopts vcodec=mpeg4:vhq:trell:mbd=2:v4mv:vb_strategy=0:vlelim=0:vcelim=0:cmp=6:subcmp=6:precmp=6:predia=3:dia=3:vme=4:vqscale=1 \"mf://./video2/*.png\" -mf type=png:fps=5 -o video_2.avi")



#************** Let's make the third video ****************************
call(["mkdir", "video3"]) #I am creating a new directory to not create confusion
name = "./video3/wealth_3_"
filetype = '.png'

# Let's make all the pictures
for i in range(timepoints+1):

	#Let's compute the histogram
	fuffo=zeros(KK) #Here I just fill the histofram with zeros
	test=data3[i]
	
	#Here I count how many hits per bin
	
	for j in range(N): 
		dummy1=test[j]*KK
		dummy2=floor(dummy1)
		fuffo[dummy2]=fuffo[dummy2] + 1
    
   	
   	# The actual plot
	plt.plot(raggio,fuffo)
	plt.axis([0, 1, 0, N])
	plt.title("Wealth distribution in case "+legend3+" is rewarded")
	plt.xlabel('normalized wealth w')
	plt.ylabel('Amount of people with wealth w')
	
	# The clock, this is an inset axes over the main axes where I plot the clock
	time=double(i)/timepoints #Here I compute the fraction of the pie that has to be filled, i.e. the percentage of the time
	a = axes([.65, .6, .2, .2], axisbg='w',aspect=1.0)
	fraction=[1-time,time]
	mycolors=["white","black"]
	plt.pie(fraction,colors=mycolors,startangle=90)

	
	# save plot as png file
	if i<10: ending = '00' + str(i)
	if 10 <= i and i<100: ending = '0' + str(i)
	if 100<= i and i<1000: ending = str(i)
	#if 1000<= i and i<10000: ending = str(i)
	filename = name + ending + filetype
	plt.savefig(filename)
	plt.close()
	print i
	

# Now we make the third movie
os.system("mencoder -really-quiet -mc 0 -noskip -skiplimit 0 -ovc lavc -lavcopts vcodec=mpeg4:vhq:trell:mbd=2:v4mv:vb_strategy=0:vlelim=0:vcelim=0:cmp=6:subcmp=6:precmp=6:predia=3:dia=3:vme=4:vqscale=1 \"mf://./video3/*.png\" -mf type=png:fps=5 -o video_3.avi")



#************** Let's make the fourth video ****************************
call(["mkdir", "video4"]) #I am creating a new directory to not create confusion
name = "./video4/wealth_4_"
filetype = '.png'

# Let's make all the pictures
for i in range(timepoints+1):

	#Let's compute the histogram
	fuffo=zeros(KK) #Here I just fill the histofram with zeros
	test=data4[i]
	
	#Here I count how many hits per bin
	
	for j in range(N): 
		dummy1=test[j]*KK
		dummy2=floor(dummy1)
		fuffo[dummy2]=fuffo[dummy2] + 1
    
   	
   	# The actual plot
	plt.plot(raggio,fuffo)
	plt.axis([0, 1, 0, N])
	plt.title("Wealth distribution in case "+legend4+" is rewarded")
	plt.xlabel('normalized wealth w')
	plt.ylabel('Amount of people with wealth w')
	
	# The clock, this is an inset axes over the main axes where I plot the clock
	time=double(i)/timepoints #Here I compute the fraction of the pie that has to be filled, i.e. the percentage of the time
	a = axes([.65, .6, .2, .2], axisbg='w',aspect=1.0)
	fraction=[1-time,time]
	mycolors=["white","black"]
	plt.pie(fraction,colors=mycolors,startangle=90)

	
	# save plot as png file
	if i<10: ending = '00' + str(i)
	if 10 <= i and i<100: ending = '0' + str(i)
	if 100<= i and i<1000: ending = str(i)
	#if 1000<= i and i<10000: ending = str(i)
	filename = name + ending + filetype
	plt.savefig(filename)
	plt.close()
	print i
	

# Now we make the fourth movie
os.system("mencoder -really-quiet -mc 0 -noskip -skiplimit 0 -ovc lavc -lavcopts vcodec=mpeg4:vhq:trell:mbd=2:v4mv:vb_strategy=0:vlelim=0:vcelim=0:cmp=6:subcmp=6:precmp=6:predia=3:dia=3:vme=4:vqscale=1 \"mf://./video4/*.png\" -mf type=png:fps=5 -o video_4.avi")



