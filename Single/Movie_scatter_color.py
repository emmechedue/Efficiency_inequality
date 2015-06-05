#*****************HERE I MAKE A MOVIE THAT SHOWS HOW THE WEALTH DISTRIBUTION EVOLVES DURING TIME ****************************

import numpy
import matplotlib
import matplotlib.pyplot as plt
from pylab import *
from numpy import array
from configobj import ConfigObj
from subprocess import call
import os
from textwrap import wrap
import sys, getopt


def main(argv):
	# Update the default values with commandline options if necessary
	simuNR = ""					# Simulation specifier
	rankSys = 1					# Ranking system
	paramFile = "./parameters.txt"			# Parameter file
	outDIR = "./videoscatter"					# Output DIR
	outBase = "scatter_ranking_" + str(rankSys)+"_" # Output file name base

	try:
		opts, args = getopt.getopt(argv,"hs:r:p:o:n:")
	except getopt.GetoptError:
		print 'Wrong usage, see -h'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
	        	print 'python Movie_scatter_color.py \n'
	        	print '\n' 
	        	print '-s <simulation specifier> \n' 
	        	print '\twealth$s.txt\n'
	        	print '\tcooperation$s.txt\n'
	        	print '\ttalent$s.txt\n'
	        	print '\tDEFAULT="" \n'
	        	print '\n'
	        	print '-r <Ranking type>\n'
	        	print '\t1 -> Everything: O[i]=r[i]*alpha[i]*w[i]\n'
	        	print '\t2 -> Money contributed: O[i]=alpha[i]*w[i]\n'
	        	print '\t3 -> Western meritocracy: O[i]=r[i]*alpha[i]\n'
	        	print '\t4 -> Willingness to contribute: O[i]=alpha[i]\n'
	        	print '\tDEFAULT=1 \n'		
	        	print '\n'
	        	print '-p <parameter file>\n'
	        	print '\tDEFAULT="./parameters.txt" \n'
	        	print '-o <output directory>\n'
	        	print '\tDEFAULT="./videoscatter" \n'
	        	print '-n <output file names>\n'
	        	print '\tDEFAULT="scatter_ranking_$r_*" \n'
			sys.exit()
		elif opt in ("-s"):
			simuNR = arg		
		elif opt in ("-r"):
			print(arg)
			rankSys = int(arg)
			outBase = "scatter_ranking_" + str(rankSys)+"_"
		elif opt in ("-p"):
			paramFile = arg
		elif opt in ("-o"):
			outDIR = arg
		elif opt in ("-n"):
			outBase = arg

        # Now setup the plotting info 

	# Time info
        config = ConfigObj(paramFile)
        dummy = config["T"]
        T=int(dummy)
        
	# Ranking system string
	leg_rank = ""
	if (rankSys == 1):
		leg_rank=r'$\alpha_i  w_i  r_i = O_i$'
	if (rankSys == 2):
		leg_rank=r'$\alpha_i  w_i$'
	if (rankSys == 3):
		print 'test'
		leg_rank=r'$\alpha_i  r_i$'
	if (rankSys == 4):
		leg_rank=r'$\alpha_i$' 	

	# Read data from input files
        w_table=np.loadtxt("./wealth" + simuNR + ".txt")
        c_table=np.loadtxt("./cooperation" + simuNR + ".txt")
        r_table=np.loadtxt("./talent" + simuNR + ".txt", usecols=(1,), unpack=True)

	# Make output dir
        call(["mkdir", outDIR]) 
        name = outDIR + "/" + outBase
        filetype = '.png'
        
        # Coloring schema info
        color = r_table
        rmax = np.amax(r_table)
        rmin = np.amin(r_table)
        
        # The limiting values of y-axis
        YL = np.amax(w_table)
        #Rmax = 5*np.amax(color) + 4
        
        # Let's make all the pictures
        for i in range(T+1): #It is +1 because I have the 0th time step
        	
        	average = np.mean(w_table[i]) #This is the average wealth at time i
        	area = np.pi * ((5/YL) * w_table[i] + 4)**2
        	# The actual plot
        	fig = plt.figure()
        
        	ax = fig.add_subplot(111)
        	
        	l = ax.scatter(c_table[i], w_table[i], c=color, s = area, marker='o', vmin=rmin, vmax=rmax)
        	fig.colorbar(l)
        	plt.title("Scatter plot of wealth and cooperation \n in case "+leg_rank+" is rewarded")
        	#title = ax.set_title("\n".join(wrap("Scatter plot of wealth and cooperation in case "+legend1+" is rewarded", 60)))
        
        	#plt.ylabel('Normalized wealth')
        	plt.ylabel('Wealth')
        	plt.xlabel('Percentage of contribution')
        
        	#plt.plot(c_average[i],average, 'k_',ms=16,mew=3)
        	#plt.plot(c_average[i],average, 'k|',ms=16,mew=3)
        	ax.set_ylim([0,YL])
        	ax.set_xlim([0,1])	
        
        	# save plot as png file
        	if i<10: ending = '00' + str(i)
        	if 10 <= i and i<100: ending = '0' + str(i)
        	if 100<= i and i<1000: ending = str(i)
        	#if 1000<= i and i<10000: ending = str(i)
        	filename = name + ending + filetype
        	plt.savefig(filename)
        	plt.close()
        	print i
        	
        
  
		os.system("mencoder -really-quiet -mc 0 -noskip -skiplimit 0 -ovc lavc -lavcopts vcodec=mpeg4:vhq:trell:mbd=2:v4mv:vb_strategy=0:vlelim=0:vcelim=0:cmp=6:subcmp=6:precmp=6:predia=3:dia=3:vme=4:vqscale=1 \"mf://"+outDIR+"/*.png\" -mf type=png:fps=5 -o "+outBase+".avi")

if __name__ == "__main__":
   main(sys.argv[1:])



















