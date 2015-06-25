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
#include "sml_Constants.h"
#include "sml_Fcts.h"


using namespace std;

int main() {
	// ********************************************************** //
	// ********************* INITIALISATION ********************* //		
	// ********************************************************** //

	// Print for testing
	cout << "Starting initialising" << endl;


	// *********** Initialising fields and variables *********** //
	
	// RNG Stuff 	
	gsl_rng *gslpointer; 	// Pointer to the type of rng
	FILE *pfile; 		// File to read from /usr/urandom
	unsigned int seed; 	// Seed of the random number generator

	// My Constants instance for reading/holding all the input varibles 	
	Constants cons;

	// Variables for describing the agents ang groups
	int M; 			// Number of groups
	//int k; 		// Index of agent that at that particular time is changing its mind // !!! NOT NECESSARY IN THIS NEW SCHEMA !!! 
	double r[cons.N];	// Array to hold talent of each agent
	double w[cons.N];	// Array to hold wealth of each agent
	double OLDw[cons.N];	// Array to hold OLDwealth of each agent
	double alpha[cons.N];	// Array to hold contribution coefficient of each agent (contribution == stragety)
	int rank[cons.N];	// Ranking array: shows in which order the agents are. 
				// IMPORTANT CONVENTION: The HIGHER the ORDER of agent, the BETTER the PLACEMENT !!!
	double O[cons.N];	// Array to hold the output for each turn 

	// Variables for the memory thingy
	int N_alpha = cons.monincr + 1; 		// The number of different strategies possible
	double alpha_interval = 1.0/cons.monincr;	// The interval lenght of the contribution discritasation
	double *memory = new double [cons.N*N_alpha]; 	// The array for the memory
	int current_alpha[cons.N];			// The array to hold the index of the currently used strategies.
	double mUpDate=cons.mUpDate;			// The coef. for memory update 

	// Supporting variables for different purposes 
	//double oldt, t; 	// The time and an time index variable to keep track of time
	//double oldt;		// The oldtime // !!! NOT NECESSARY IN THIS NEW SCHEMA !!! 
	int t;
	//int enne; 		// Int variable that I need in a time thing // !!! NOT NECESSARY IN THIS NEW SCHEMA !!! 
	int i; 			// Support variable 
	//double ranvar;	// A random variable // !!! NOT NECESSARY IN THIS NEW SCHEMA !!! 
	double eff;	// These is to print the efficiency and then renormalize the wealth. Here the efficiency is defined as percentage growth
	double wealth, oldwealth; //Wealth is the wealth at every time step. Oldwealth is a support variable I need to compute the efficiency 

	// ofstream for outputing and output files
	ofstream filep;		// Parameter output stream 
	ofstream filet;		// Main/time output stream
	ofstream filew;		// Wealth output stream
	ofstream filec;		// Cooperation output stream
	ofstream filer;		// Talent output stream
	ofstream filem;		// Memory output stream
	
	const char filenamep[]="parameters.txt";	// Parameter output filename	
	const char filenamet[]="time.txt";		// Main/time output stream
	const char filenamew[]="wealth.txt"; 		// Wealth output filename
	const char filenamec[]="cooperation.txt";	// Cooperation output filename
	const char filenamer[]="talent.txt";		// Talent output filename
	const char filenamem[]="memory.txt";		// Memory output filename



	// *********** Initialising the GSL Mersenne twister RNG & a C srand RNG *********** //

	// Print for testing
	cout << "Initialising the RNG" << endl;


    	// GSL Mersenne Twister 19937 RNG
	pfile = fopen ("/dev/urandom", "r");
	i=fread (&seed, sizeof (seed), 1, pfile);	// I added the rand= ... just to not be bothered anymore by the warnings!
	fclose(pfile);
	gslpointer = gsl_rng_alloc(gsl_rng_mt19937); 	// I'm using the "Mersenne Twister" generator!
	gsl_rng_set(gslpointer,seed); 			// Starting the generator
	
	// C srand RNG
	srand(seed); 					// It's unlikely, but I might need it to break ties in the ordering

	// *********** Initialising the OUTPUT FILES with HEADERS *********** //

	// Print for testing
	cout << "Initialising the output files" << endl;


	// The WEALTH file 	
	filew.open(filenamew,ios::out|ios::trunc); 
	if(filew.is_open()){
		filew << "# Wealth at each time step for the simulation Efficiency_inequality with:"<<endl;
		filew << "#N=" << cons.N << " T="<<cons.T  << " S=";
		filew << cons.S << " Q=" << cons.Q << " mu="<<cons.mu;
		filew << " sigmag=" << cons.sigmag << " W0=" << cons.W0 << " 	choice=" << cons.Choice;
		filew << " monincr=" << cons.monincr << " beta=" << cons.beta << " seed=" << seed << endl;
	}

	// The MEMORY file 	
	filem.open(filenamem,ios::out|ios::trunc); 
	if(filem.is_open()){
		filem << "# Memory at each time step for the simulation Efficiency_inequality with:"<<endl;
		filem << "#N=" << cons.N << " T="<<cons.T  << " S=";
		filem << cons.S << " Q=" << cons.Q << " mu="<<cons.mu;
		filem << " sigmag=" << cons.sigmag << " W0=" << cons.W0 << " 	choice=" << cons.Choice;
		filem << " monincr=" << cons.monincr << " beta=" << cons.beta << " seed=" << seed << endl;
	}

	// The COOPERATION file
	filec.open(filenamec,ios::out|ios::trunc);
	if(filec.is_open()){ 
		filec << "#Strategy at each time step for the simulation Efficiency_inequality with:"<<endl;
		filec << "#N=" << cons.N << " T=" << cons.T  << " S=";
		filec << cons.S << " Q=" << cons.Q  << " mu=" << cons.mu;
		filec << " sigmag=" << cons.sigmag << " W0=" << cons.W0 << " choice=" << cons.Choice;
		filec << " monincr=" << cons.monincr << " beta=" << cons.beta << " seed=" << seed << endl;
	}

	// The TIME file
	filet.open(filenamet,ios::out|ios::trunc); //Open the time file
	if(filet.is_open()){
		filet << "#Results for the simulation Efficiency_inequality with:"<<endl;
		filet << "#N=" << cons.N << " T=" << cons.T  << " S=";
		filet << cons.S << " Q=" << cons.Q <<  " mu=" << cons.mu;
		filet << " sigmag=" << cons.sigmag << " W0=" << cons.W0 << " choice=" << cons.Choice;
		filet << " monincr=" << cons.monincr << " beta=" << cons.beta << " seed=" << seed << endl;
		filet << "#This is in the form of t, Gini coefficient, growth percentage and average cooperation" << endl;
	}

	// The TALENT file 
	filer.open(filenamer,ios::out|ios::trunc); //Open the talent file
	if(filer.is_open()){
		filer << "# Talent for each player, with simulation configurations:"<<endl;
		filer << "#N=" << cons.N << " T=" << cons.T << " S=";
		filer << cons.S << " Q=" << cons.Q <<  " mu=" << cons.mu;
		filer << " sigmag=" << cons.sigmag << " W0=" << cons.W0 << " choice=" << cons.Choice;
		filer << " monincr=" << cons.monincr << " beta=" << cons.beta << " seed=" << seed << endl;
		filer << "#" << endl;
		filer << "# Player\tTalent" << endl;
	}

	// The PARAMETER file	
	filep.open(filenamep,ios::out|ios::trunc);
	printparamsingleloop(filep,cons);
	filep.close();


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

	// The arrays describing agents	
	for(i=0;i<cons.N;i++){
	
		// Initial wealth W0 is equal for all
    		w[i]=cons.W0;		
		// Intial rank is effectivly random
    		rank[i]=i;		
		
		// Gaussian distributed talent. Have to add the mean of the gaussian because the generator has mean zero.  
    		r[i]=gsl_ran_gaussian(gslpointer,cons.sigmag)+cons.mu;

		// Print out the talent for each player 
		filer << i << "\t" << r[i] << endl;

		// All initial strategies are full defection => contribution of all agents is 0
    		alpha[i]=0; 

		// The memory array 
		for (int j=0; j<N_alpha; j++){
			int ID = i*N_alpha + j;
			memory[ID] = 1.0;
		}
	}

	// Close the talent file as it is not used anymore
	filer.close();

	// Renormalise wealth w[]
	//renormalizearray(w,cons.N);

	//Total initial wealth. I don't really need to initialize this, but it's better anyway
	wealth = sumofvector(w, cons.N);

	//It is important to initialize this!
	oldwealth=wealth; 

	// Now I fill up the output files at t=0
	printstuffsingleloop(filet, filew, filec, t, eff, w,alpha, cons); 	// Print time, gini coefficient, efficiency and all the wealth of the initial state.

	// Print memory for testing purposes
	printMemory(filem, t, cons.N, N_alpha, memory);

	// Print for testing
	cout << "Initialisation complete" << endl;



	// ********************* END OF INITIALISATION ********************* //		


	// ************************************************************* //
	// ********************* THE BIG TIME LOOP ********************* //		
	// ************************************************************* //

	for ( t=1; t<(cons.T+1); t++ ){

		// ** Save the previous wealth vector --> for memory updating ** //
		for (int j=0; j<cons.N; j++){
			OLDw[j] = w[j];
		}

		// *** Update strategy for all players *** //
		for( i=0; i<cons.N; i++){

			// The strategy is updated based on the memory of the agent
			//updatestrategy(i,alpha,oldalpha,w,r,cons, gslpointer);
			int tmp = newStrategySML(i, memory, N_alpha, cons, gslpointer);
			current_alpha[i] = tmp;
			alpha[i] = tmp*alpha_interval;
		}

		// *** "Play the game" with the new strategys *** //

		// I generate the ranking and compute the output for each player.
		makeorder(r,w,alpha,rank,cons,O); 

		// Here I use the ranking to generate the groups and update the wealth of the players.
		formgroups(w, alpha, O,rank, cons, M); 

		// Here I compute the total amount of wealth. And the efficiency, which is defined as the percentage increase of wealth.
		wealth = sumofvector(w, cons.N);
		eff = wealth - oldwealth;
		eff = eff/oldwealth; // I separate those 2 because in this way it should be better for numerical errors, right?
		oldwealth = wealth;		 

		// ** Update the memory, using the growth of wealth ** //
		for (int j=0; j<cons.N; j++){

			// Calculate the growth for the agent
			double tmp_growth = (w[j] - OLDw[j])/OLDw[j];

			// Get the memory ID
			int ID = j*N_alpha + current_alpha[j];

			// Find the current memory value
			float tmp_mem = memory[ID];

			// Calculate the new memory
			memory[ID] = tmp_mem*(1.0 + mUpDate*tmp_growth);			
		}


		// Print time, Gini coefficient, efficiency and all the wealth for the aldt timestep.
		printstuffsingleloop(filet, filew, filec, t, eff, w, alpha,cons); 

		// Print memory for testing purposes
		printMemory(filem, t, cons.N, N_alpha, memory);


		// Print the time to commandline for checking the runtime.
		cout<<"Timestep "<< t << endl; 		


	}



	// ********************* END OF BIG TIME LOOP ********************* //	

	// ******************************************************************* //
	// ********************* CLOSE ALL THE OFSTRESMS ********************* //		
	// ******************************************************************* //

	filet.close();
	filew.close();
	filec.close();
	filem.close();
	return 0;
}
