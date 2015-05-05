/* 
 * File:   main.cpp
 * Author: sduca
 *
 * Created on February 4, 2015, 4:51 PM
 */

// Standard C++ library includes 
#include<cstdio>
#include<cstdlib>
#include<fstream>
#include<math.h>
#include<iostream>
#include <new>
#include <algorithm>    
#include <vector> 

// Gnu Scientific Library (gsl) includes
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h> //Needed for the beta pdf

// My includes
#include "Constants.h"
#include "Fcts.h"


using namespace std;

int main() {
	// ********************************************************** //
	// ********************* INITIALISATION ********************* //		
	// ********************************************************** //

	// *********** Initialising fields and variables *********** //
	
	// RNG Stuff 	
	gsl_rng *gslpointer; 	// Pointer to the type of rng
	FILE *pfile; 		// File to read from /usr/urandom
	unsigned int seed; 	// Seed of the random number generator

	// My Constants instance for reading/holding all the input varibles 	
	Constants cons;

	// Variables for describing the agents ang groups
	int M; 			// Number of groups
	int k; 			// Index of agent that at that particular time is changing its mind
	double r[cons.N];	// Array to hold talent of each agent
	double w[cons.N];	// Array to hold wealth of each agent
	double alpha[cons.N];	// Array to hold contribution coefficient of each agent 
	int rank[cons.N];	// Ranking array: shows in which order the agents are. 
				// IMPORTANT CONVENTION: The HIGHER the ORDER of agent, the BETTER the PLACEMENT !!!
	double O[cons.N];	// Array to hold the output for each turn 

	// Supporting variables for different purposes 
	double oldt, t; 	// The time and an time index variable to keep track of time
	int enne; 		// Int variable that I need in a time thing
	int i; 			// Support variable
	double ranvar;		// A random variable
	double eff, oldeff;	// These are to print the efficiency and then renormalize the wealth. 
				// Here the efficiency is defined as percentage growth
	int counteff;		// This I need to make proper averages of the efficiency

	// ofstream for outputing and output files
	ofstream filep;		// Parameter output stream 
	ofstream filet;		// Main/time output stream
	ofstream filew;		// Wealth output stream
	ofstream filec;		// Cooperation output stream
	
	const char filenamep[]="parameters.txt";	// Parameter output filename	
	const char filenamet[]="time.txt";		// Main/time output stream
	const char filenamew[]="wealth.txt"; 		// Wealth output filename
	const char filenamec[]="cooperation.txt";	// Cooperation output filename



	// *********** Initialising the GSL Mersenne twister RNG & a C srand RNG *********** //

    	// GSL Mersenne Twister 19937 RNG
	pfile = fopen ("/dev/urandom", "r");
	i=fread (&seed, sizeof (seed), 1, pfile);	// I added the rand= ... just to not be bothered anymore by the warnings!
	fclose(pfile);
	gslpointer = gsl_rng_alloc(gsl_rng_mt19937); 	// I'm using the "Mersenne Twister" generator!
	gsl_rng_set(gslpointer,seed); 			// Starting the generator
	
	// C srand RNG
	srand(seed); 					// It's unlikely, but I might need it to break ties in the ordering

	
	// *********** Initialising the OUTPUT FILES with HEADERS *********** //

	// The WEALTH file 	
	filew.open(filenamew,ios::out|ios::trunc); 
	if(filew.is_open()){
		filew << "#Wealth at each time step for the simulation Efficiency_inequality with:"<<endl;
		filew << "#N=" << cons.N << " T="<<cons.T << " interval=" << cons.interval << " S=";
		filew << cons.S << " Q=" << cons.Q << " lambda=" << cons.lambda << " mu="<<cons.mu;
		filew << " sigmag=" << cons.sigmag << " W0=" << cons.W0 << " 	choice=" << cons.Choice;
		filew << " monincr=" << cons.monincr << " beta=" << cons.beta << " seed=" << seed << endl;
	}

	// The COOPERATION file
	filec.open(filenamec,ios::out|ios::trunc);
	if(filec.is_open()){ 
		filec << "#Strategy at each time step for the simulation Efficiency_inequality with:"<<endl;
		filec << "#N=" << cons.N << " T=" << cons.T << " interval=" << cons.interval << " S=";
		filec << cons.S << " Q=" << cons.Q << " lambda=" << cons.lambda << " mu=" << cons.mu;
		filec << " sigmag=" << cons.sigmag << " W0=" << cons.W0 << " choice=" << cons.Choice;
		filec << " monincr=" << cons.monincr << " beta=" << cons.beta << " seed=" << seed << endl;
	}

	// The TIME file
	filet.open(filenamet,ios::out|ios::trunc); //Open the time file
	if(filet.is_open()){
		filet << "#Results for the simulation Efficiency_inequality with:"<<endl;
		filet << "#N=" << cons.N << " T=" << cons.T << " interval=" << cons.interval << " S=";
		filet << cons.S << " Q=" << cons.Q << " lambda=" << cons.lambda << " mu=" << cons.mu;
		filet << " sigmag=" << cons.sigmag << " W0=" << cons.W0 << " choice=" << cons.Choice;
		filet << " monincr=" << cons.monincr << " beta=" << cons.beta << " seed=" << seed << endl;
		filet << "#This is in the form of t, Gini coefficient, growth percentage and average cooperation" << endl;
	}



	// *********** Initialise description of agents and groups (START GENERATING ALL THE STUFF) *********** //

	// Check that number of agents N is exactly divisible by size of groups S 
	if((cons.N%cons.S) != 0){ 
		// cout<<"The number of agents is not exactly divisible by the number of groups"<<endl;
		cout<<"The number of agents is not exactly divisible by the size of groups"<<endl;
		exit(1);
	}
	// Compute number of groups M = N/S
	else{M=cons.N/cons.S;} 

	// FILL THE VECTORS OF TALENT , WEALTH AND RANKING AND INITIALIZE EVERYTHING
	t=0;
	eff=0;
	oldeff=0; 		// I don't really need to initialize this, but it's better anyway

	// The arrays describing agents	
	for(i=0;i<cons.N;i++){
	
		// ***** WHY NOT INITIALISE WEALTH OF ALL AGENTS TO BE 1 !? ***** //
		// Initial wealth W0 is equal for all
    		//w[i]=cons.W0;		
		w[i] = 1.0;		
		
		// Intial rank is efectivly random
    		rank[i]=i;		
		
		// Gaussian distributed talent. Have to add the mean of the gaussian because the generator has mean zero.  
    		r[i]=gsl_ran_gaussian(gslpointer,cons.sigmag)+cons.mu;

		// All initial trategies are full defection = contribution of all agents is 0
    		alpha[i]=0; 
	}
	// ***** WHY NOT INITIALISE WEALTH OF ALL AGENTS TO BE 1 !?  --> DON'T RENORMALIZE! ***** //
	// Renormalise wealth w[]
	//renormalizearray(w,cons.N);

	// Now I fill up the output files at t=0
	printstuffsingleloop(filet, filew, filec, t, eff, w,alpha, cons);	// Print time, gini coefficient, efficiency 
									 	// and all the wealth of the initial state.
	// RESET the counter for the AVERAGE of the EFFICIENTCY	
	counteff = -1; 


	// ********************* END OF INITIALISATION ********************* //		


	// ************************************************************* //
	// ********************* THE BIG TIME LOOP ********************* //		
	// ************************************************************* //
	
	do {

		// *** Sample and update the time (with an exponentially distributed RNG) *** //

		// Generate a random timestep lenght from an exponetially distributed RNG
		ranvar=gsl_ran_exponential(gslpointer,cons.lambda);

		// If the new time exceeds maximum simulation time const.T --> Simply print everything and then break the loop
		if (t + ranvar > cons.T) { 
			oldt = cons.T - t; 			// Note that here oldt has a different use, 
								// I'm simply using this variable since I don't need it anymore

			// Calculate how many printing intervals were there bitween the last timestep and the final timestep
			enne = floor(oldt / cons.interval);

			// Print the output (time, Gini coefficient, efficiency and all the wealth)
			// for all the printing intervals between last and the final timesteps.
			for (i = 1; i < (enne + 1); i++) { 

				// Note that here I have "<" of (enne + 1) and not "<=", since here I need the time explicitely!
				// The (enne + 1) is because I need it to multply. Meaning that all this is because I start from 1. 
				// And I have to do it since I still have to update the time t.

				// Print time, Gini coefficient, efficiency and all the wealth for the (t + i*cinst.interval) timestep 
				printstuffsingleloop(filet, filew, filec, t + i * cons.interval, eff, w,alpha, cons); 
				// RESET the counter for the AVERAGE of the EFFICIENTCY					
				counteff = -1; 
			}
			
			// Print the output (time, Gini coefficient, efficiency and all the wealth) for the last timestep (at time T)
			printstuffsingleloop(filet, filew,filec, cons.T, eff, w, alpha, cons);
			// RESET the counter for the AVERAGE of the EFFICIENTCY	 
			counteff = -1; 
			
			// EXIT from the DO LOOP !
			break; 
		}
		

		// Else check if the timestep exceeds the printout interval and if necessary
		// reprint the old situation before updating the system!!
		if (ranvar > cons.interval) { 
	
			// Calculate how many printing intervals were there bitween the last timestep and the new timestep
			enne = floor(ranvar / cons.interval);

			// Print the output (time, Gini coefficient, efficiency and all the wealth)
			// for all the printing intervals between last and the new timesteps.
			for (i = 1; i <= enne; i++) { 

				// Note that "<=" is used because I start from 1

				// Print time, Gini coefficient, efficiency and all the wealth for each if the printout intervals 
				printstuffsingleloop(filet, filew,filec, t + i * cons.interval, eff, w, alpha, cons); 
				// RESET the counter for the AVERAGE of the EFFICIENTCY
				counteff = -1; 

				// Print the time to commandline for checking the runtime. 
				cout << "The time is " << t + i * cons.interval << endl; 
			}
			
			// Shift the time by the printout intervals just done (const.interval*enne) and correct the timestep length (ranvar)					
			ranvar = ranvar - cons.interval*enne;
			t = t + enne * cons.interval;
		}

		// Update the time. I have to do it after I do the check for the time.
		t = t + ranvar; 
		// Update oldt.
		oldt = oldt + ranvar; 

		// *** Update the stragedy of a random agent *** //

		// Randomly select an agent (uniformly distributed RNG)
		k=gsl_rng_uniform_int(gslpointer,cons.N);

		// Here I update the strategy of agent k according to a logit distribution.
		updatestrategy(k,alpha,w,r,cons, gslpointer);
 
		// Here I would update the strategy completely at random!
		//alpha[k]=gsl_ran_beta(gslpointer,2,2); 
		
		// I generate the ranking and compute the output for each player.
		makeorder(r,w,alpha,rank,cons,O); 

		// Here I use the ranking to generate the groups and update the wealth of the players.
		formgroups(w, alpha, O,rank, cons, M); 

		// Here I compute the total amount of wealth. Since the total wealth was normalized to one before,
		// this is the increase in percentage of wealth.
		eff = sumofvector(w, cons.N) - 1. ; 

		// These 3 lines here are to compute proper averages of eff. 
		// This is in case I perform many iterations without printing.
		counteff++; 
		eff = updateaverages(oldeff, eff, counteff);
		oldeff = eff;

		// ***** DON'T RENORMALISE THE WEALTH ***** // 
		// Finally, Here I renormalize the array.
		//renormalizearray(w,cons.N); 

		// Check whether I have exceeded the printout interval.
		if(oldt>=cons.interval){ 
   
			// Print time, Gini coefficient, efficiency and all the wealth for the aldt timestep.
			printstuffsingleloop(filet, filew,filec, t, eff, w, alpha,cons); 
			// RESET the counter for the AVERAGE of the EFFICIENTCY	 
			counteff = -1;
	
			 // Subtract by oldt the value of interval to start counting again.
			oldt=oldt -cons.interval;

			// Print the time to commandline for checking the runtime.
			cout<<"The time is "<<t<<endl; 
		}
	
	}while(t<=cons.T);


	// ********************* END OF BIG TIME LOOP ********************* //	

	// *********************************************************************** //
	// ********************* PRINT PARAMETERS & FINALISE ********************* //		
	// *********************************************************************** //

	filep.open(filenamep,ios::out|ios::trunc);
	printparamsingleloop(filep,cons);
	filep.close();
	filet.close();
	filew.close();
	filec.close();
	return 0;
}
