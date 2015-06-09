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
#include <iomanip>
#include <new>
#include <algorithm>    
#include <vector> 

// Gnu Scientific Library (gsl) includes
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h> //Needed for the beta pdf

// My includes
#include "ensemble_Constants.h"
#include "ensemble_Fcts.h"


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

	// Variables for describing the agents and groups
	int NE = cons.NE;	// Number of "societies" in the ensemble
	int N = cons.N;		// Number of agents in one "society" 
	int M; 			// Number of groups
/*	double r[NE*N];		// Array to hold talent of each agent
	double w[NE*N];		// Array to hold wealth of each agent
	double alpha[NE*N];	// Array to hold contribution coefficient of each agent (contribution == strategy == co-op)
	//double oldalpha[NE*N];	// Array to hold the contribution coefficient (strategy) for each agent from the last timestep
	int rank[NE*N];		// Ranking array: shows in which order the agents are. 
				// IMPORTANT CONVENTION: The HIGHER the ORDER of agent, the BETTER the PLACEMENT !!!
	double O[N];		// Array to hold the output for each turn !!! THIS DOES NOT NEED MORE THEN N ENTRIES !!!
*/
	double r[N];		// Array to hold talent of each agent
	double w[N];		// Array to hold wealth of each agent
	double alpha[N];	// Array to hold contribution coefficient of each agent (contribution == strategy == co-op)
	double oldalpha[N];	// Array to hold the contribution coefficient (strategy) for each agent from the last timestep
	int rank[N];		// Ranking array: shows in which order the agents are. 
				// IMPORTANT CONVENTION: The HIGHER the ORDER of agent, the BETTER the PLACEMENT !!!
	double O[N];		// Array to hold the output for each turn !!! THIS DOES NOT NEED MORE THEN N ENTRIES !!!



	// Supporting variables for different purposes 
	//double oldt, t; 	// The time and an time index variable to keep track of time
	//double oldt;		// The oldtime // !!! NOT NECESSARY IN THIS NEW SCHEMA !!! 
	int t;
	//int enne; 		// Int variable that I need in a time thing // !!! NOT NECESSARY IN THIS NEW SCHEMA !!! 
	//int i; 			// Support variable 
	//double ranvar;	// A random variable // !!! NOT NECESSARY IN THIS NEW SCHEMA !!! 
	//double eff[NE];		// These are to print the efficiency and then renormalize the wealth.
	double eff;		// These are to print the efficiency and then renormalize the wealth.
	//double oldeff[NE];	//  
				// Here the efficiency is defined as percentage growth
/*	double wealth[NE];  	// Wealth is the total wealth at every time step for each society
	double oldwealth[NE];	// Oldwealth is a support variable I need to compute the efficiency	
*/	double wealth;  	// Wealth is the total wealth at every time step for each society
	double oldwealth;	// Oldwealth is a support variable I need to compute the efficiency
		
	//int counteff;		// This I need to make proper averages of the efficiency	// !!! NOT NECESSARY IN THIS NEW SCHEMA !!! 

	// ofstream for outputing and output files
	ofstream filep;		// Parameter output stream 
//	ofstream filet;		// Main/time output stream
	ofstream filew;		// Wealth output stream
	ofstream filec;		// Cooperation output stream
	ofstream filer;		// Talent output stream
	ofstream fileg;		// Gini coef. output stream
	ofstream fileEf;		// Efficientcy output stream
	
	const char filenamep[]="parameters.txt";	// File to hold the simulation input parameters	
	//const char filenamet[]="time.txt";		// File to hold the ensemble averages for each timestep: en.avg TotWealth, Growth%, Gini and avg.co-op
	const char filenamew[]="wealth.txt"; 		// File to hold the total wealth for each ensemble for each timestep
	const char filenamec[]="cooperation.txt";	// File to hold the average co-op in each "society" for each timestep
	const char filenamer[]="talent.txt";		// File to hold the talent for each individual in each society (rows: individuals, columns: societies)
	const char filenameg[]="giniCoef.txt";		// File to hold the Gini coef. for each society on each time step.
	const char filenameEf[]="efficiency.txt";		// File to hold the efficientcy for each society on each time step.


	// *********** Initialising the GSL Mersenne twister RNG & a C srand RNG *********** //

    // GSL Mersenne Twister 19937 RNG
	pfile = fopen ("/dev/urandom", "r");
	//int ir;
	fread (&seed, sizeof (seed), 1, pfile);	// I added the rand= ... just to not be bothered anymore by the warnings!
	//cout << "Why do we use this ir = fread (&seed, sizeof (seed), 1, pfile) ?" << endl;
	//cout << "The above print-out is just to get rid of the warning" << endl;
	fclose(pfile);
	gslpointer = gsl_rng_alloc(gsl_rng_mt19937); 	// I'm using the "Mersenne Twister" generator!
	gsl_rng_set(gslpointer,seed); 			// Starting the generator
	
	// C srand RNG
	srand(seed); 					// It's unlikely, but I might need it to break ties in the ordering

	
	// *********** Initialising the OUTPUT FILES with HEADERS *********** //

	// The WEALTH file 	
	filew.open(filenamew,ios::out|ios::trunc); 
	if(filew.is_open()){
		filew << "# Total wealth of each 'society' at each time step for the simulation Efficiency_inequality with:"<<endl;
		filew << "#N=" << cons.N << " T="<<cons.T  << " S=";
		filew << cons.S << " Q=" << cons.Q << " mu="<<cons.mu;
		filew << " sigmag=" << cons.sigmag << " W0=" << cons.W0 << " 	choice=" << cons.Choice;
		filew << " monincr=" << cons.monincr << " beta=" << cons.beta << " seed=" << seed << endl;
		//filew << "Time step\tTotal wealth for each 'society'" << endl;
		filew << "# Total wealth -> Rows: 'societies'\tColumns: timesteps" << endl;	
		filew << "#" << endl;
	}

	// The WEALTH file 	
	fileEf.open(filenameEf,ios::out|ios::trunc); 
	if(filew.is_open()){
		fileEf << "# Efficiency of each 'society' at each time step for the simulation Efficiency_inequality with:"<<endl;
		fileEf << "#N=" << cons.N << " T="<<cons.T  << " S=";
		fileEf << cons.S << " Q=" << cons.Q << " mu="<<cons.mu;
		fileEf << " sigmag=" << cons.sigmag << " W0=" << cons.W0 << " 	choice=" << cons.Choice;
		fileEf << " monincr=" << cons.monincr << " beta=" << cons.beta << " seed=" << seed << endl;
		//filew << "Time step\tTotal wealth for each 'society'" << endl;
		fileEf << "# Efficiency -> Rows: 'societies'\tColumns: timesteps" << endl;	
		filew << "#" << endl;
	}

	// The COOPERATION file
	filec.open(filenamec,ios::out|ios::trunc);
	if(filec.is_open()){ 
		filec << "# Average strategy/co-op for each 'society' at each time step for the simulation Efficiency_inequality with:"<<endl;
		filec << "#N=" << cons.N << " T=" << cons.T  << " S=";
		filec << cons.S << " Q=" << cons.Q  << " mu=" << cons.mu;
		filec << " sigmag=" << cons.sigmag << " W0=" << cons.W0 << " choice=" << cons.Choice;
		filec << " monincr=" << cons.monincr << " beta=" << cons.beta << " seed=" << seed << endl;
		//filec << "Time step\tAverage co-op/strategy for each 'society'" << endl;	
		filec << "# Average strategy/co-op -> Rows: 'societies'\tColumns: timesteps" << endl;		
		filec << "#" << endl;
	}

	// The GINI COEF. file
	fileg.open(filenameg,ios::out|ios::trunc);
	if(filec.is_open()){ 
		fileg << "# Gini coef. for each 'society' at each time step for the simulation Efficiency_inequality with:"<<endl;
		fileg << "#N=" << cons.N << " T=" << cons.T  << " S=";
		fileg << cons.S << " Q=" << cons.Q  << " mu=" << cons.mu;
		fileg << " sigmag=" << cons.sigmag << " W0=" << cons.W0 << " choice=" << cons.Choice;
		fileg << " monincr=" << cons.monincr << " beta=" << cons.beta << " seed=" << seed << endl;
		//fileg << "Time step\tGini coef. for each 'society'" << endl;
		fileg << "# Gini coef -> Rows: 'societies'\tColumns: timesteps" << endl;	
		fileg << "#" << endl;
	}

	/*// The TIME file
	filet.open(filenamet,ios::out|ios::trunc); //Open the time file
	if(filet.is_open()){
		filet << "# Ensemble averages at each time step for the simulation Efficiency_inequality with:"<<endl;
		filet << "#N=" << cons.N << " T=" << cons.T  << " S=";
		filet << cons.S << " Q=" << cons.Q <<  " mu=" << cons.mu;
		filet << " sigmag=" << cons.sigmag << " W0=" << cons.W0 << " choice=" << cons.Choice;
		filet << " monincr=" << cons.monincr << " beta=" << cons.beta << " seed=" << seed << endl;
		filet << "# Time step\tGini coefficient\tGrowth percentage\tAverage cooperation" << endl;
		filet << "#" << endl;
	}*/

	// The TALENT file 
	filer.open(filenamer,ios::out|ios::trunc); //Open the talent file
	if(filer.is_open()){
		filer << "# Talent for each player in each society with simulation configurations:"<<endl;
		filer << "#N=" << cons.N << " T=" << cons.T << " S=";
		filer << cons.S << " Q=" << cons.Q <<  " mu=" << cons.mu;
		filer << " sigmag=" << cons.sigmag << " W0=" << cons.W0 << " choice=" << cons.Choice;
		filer << " monincr=" << cons.monincr << " beta=" << cons.beta << " seed=" << seed << endl;
		filer << "#" << endl;
		filer << "# Rows: 'societies'\tColumns: individuals/players" << endl;
		filer << "#" << endl;
	}

	// The PARAMETER file	
	filep.open(filenamep,ios::out|ios::trunc);
	ensemblePrintParams(filep,cons);
	filep.close();


	// ---- Get the number of GROUPS ---- //
	// Check that number of agents N is exactly divisible by size of groups S 
	if((cons.N%cons.S) != 0){ 
		// cout<<"The number of agents is not exactly divisible by the number of groups"<<endl;
		cout<<"The number of agents is not exactly divisible by the size of groups"<<endl;
		exit(1);
	}
	// Compute number of groups M = N/S
	else{M=cons.N/cons.S;} 



	// ----------- LOOP OVER THE ENSEMBLE ------------ //
	for (int k=0; k<NE; k++){

		// *********** Initialise description of agents and groups (START GENERATING ALL THE STUFF) *********** //
	

		// FILL THE VECTORS OF TALENT , WEALTH AND RANKING AND INITIALIZE EVERYTHING
		t=0;			
	
	
		// The arrays describing agents	& eff
		//for(int k=0; k<NE; k++){
		
			eff = 0.0;
			//oldeff[k] = 0.0;
	
			for(int i=0; i<N; i++){
		
				// Index of the agent
				//int ID = k*N + i;
	
				// Initial wealth W0 is equal for all
		    		//w[ID]=cons.W0;
				w[i]=cons.W0;		
	
				// Intial rank is efectivly random
		    		//rank[ID]=i;
				rank[i]=i;		
				
				// Gaussian distributed talent. Have to add the mean of the gaussian because the generator has mean zero.  
		    		//r[ID]=gsl_ran_gaussian(gslpointer,cons.sigmag)+cons.mu;
				r[i]=gsl_ran_gaussian(gslpointer,cons.sigmag)+cons.mu;
	
				// All initial trategies are full defection = contribution of all agents is 0
		    		//alpha[ID]=0; 
				alpha[i]=0;
			}
		//}
	
		// Print out the talent for each player 
		for(int i=0; i<N; i++){
			//for(int k=0; k<NE; k++){
				//int ID = k*NE + i;		
				//filer << r[ID] << "\t" ;
				//filer << left << setw(8) << r[ID] << "\t" ;
				filer << left << setw(8) << r[i] << "\t" ;
			//}
			
		}
		filer << endl;	
	
		// Close the talent file as it is not used anymore
		//filer.close();
	
		// Renormalise wealth w[]
		//renormalizearray(w,cons.N);
	
		//for(int k=0; k<NE; k++){
			
		/*	// Total initial wealth. I don't really need to initialize this, but it's better anyway
			int begin = N*k;
			int end = N*(k+1) - 1;
			wealth[k] = sumofSubVector(w, begin, end);
	
			//It is important to initialize this!
			oldwealth[k] = wealth[k]; 
		//}
*/
		// Total initial wealth. I don't really need to initialize this, but it's better anyway
		wealth = sumofvector(w, N);

		//It is important to initialize this!
		oldwealth = wealth;
	
	
		// Print the initial state // 
		// First column: Time Step = 0, Following Columns: 
		// filet/time.txt --> Ensemble averages for total wealth, growth%/eff, Gini coef., avg Co-op/strategy
		// filew/wealth.txt --> Total wealth for each 'society'
		// filec/cooperation.txt --> Average co-op/strategy for each 'society'
		// fileg/GiniCoef.txt --> Gini coef. for each 'society' 
		//ensemblePrint(filet, filew, filec, fileg, t, eff, wealth, w, alpha, cons, N, NE);
		ensemblePrint(fileEf, filew, filec, fileg, t, eff, wealth, w, alpha, cons, N);
	
	
	
		// ********************* END OF INITIALISATION ********************* //		
	
	
		// ************************************************************* //
		// ********************* THE BIG TIME LOOP ********************* //		
		// ************************************************************* //
	
		for ( t=1; t<(cons.T+1); t++ ){
	
			// *** Save the old strategy *** // 		
			for(int i=0; i<(N); i++){
				oldalpha[i] = alpha[i];
			}
			
			// *** Loop over all NE societies *** //
			//for (int k=0; k<NE; k++){
				/*
				// ** Subarrays corresponding to k-th society ** //
				double alpha_k[N];
				FillSubArray(alpha, alpha_k, k, N);			
				double w_k[N];
				FillSubArray(w, w_k, k, N);
				double r_k[N];
				FillSubArray(r, r_k, k, N);
				int rank_k[N];
				FillSubArray(rank, rank_k, k, N);
			*/	
			
				// ** Update strategy for each player in each in k-th society ** //
				for (int i=0; i<N; i++){
					//alpha[k*N + i] = NewStrategy(i, alpha_k, w_k, r_k, cons, gslpointer);
					alpha[i] = NewStrategy(i, oldalpha, w, r, cons, gslpointer);
				}
	
				// ** Now !!! CHANGE/UPDATE !!! the values in !!! alpha_k !!! ** //
				/*for (int i=0; i<N; i++){
					alpha_k[i] = alpha[k*N + i];
				}*/
			
				// ** "Play the game" with the new strategys ** //
	
				// I generate the ranking and compute the output for each player.
				//makeorder(r_k, w_k, alpha_k, rank_k, cons, O);
				makeorder(r, w, alpha, rank, cons, O); 
	
				// Here I use the ranking to generate the groups and update the wealth of the players.
				//formgroups(w_k, alpha_k, O, rank_k, cons, M);
				formgroups(w, alpha, O, rank, cons, M); 
	
				// Here I compute the total amount of wealth. And the efficiency, which is defined as the percentage increase of wealth.
/*				wealth[k] = sumofvector(w_k, N);
				eff[k] = wealth[k] - oldwealth[k];
				eff[k] = eff[k]/oldwealth[k]; // I separate those 2 because in this way it should be better for numerical errors, right?
				oldwealth[k] = wealth[k];		 		
				//oldeff[k] = eff[k];
*/
				wealth = sumofvector(w, N);
				eff = wealth - oldwealth;
				eff = eff/oldwealth; // I separate those 2 because in this way it should be better for numerical errors, right?
				oldwealth = wealth;		 		
				//oldeff[k] = eff[k];
	
				/*// ** UPDATE the !! BigArrays !! ** //
				for (int i=0; i<N; i++){
					rank[k*N + i] = rank_k[i];
					w[k*N + i] = w_k[i];
				}*/
			
			//}
	
			// *** Outputting *** // 
	
	
			// Print the state at time step t // 
			// First column: Time Step = t, Following Columns: 
			// filet/time.txt --> Ensemble averages for total wealth, growth%/eff, Gini coef., avg Co-op/strategy
			// filew/wealth.txt --> Total wealth for each 'society'
			// filec/cooperation.txt --> Average co-op/strategy for each 'society'
			// fileg/GiniCoef.txt --> Gini coef. for each 'society' 
			//ensemblePrint(filet, filew, filec, fileg, t, eff, wealth, w, alpha, cons, N, NE);
			ensemblePrint(fileEf, filew, filec, fileg, t, eff, wealth, w, alpha, cons, N);
				
	
	
			// Print the time to commandline for checking the runtime //
			cout << "Society nr " << k + 1 << ", timestep "<< t << endl; 		
	
	
		}
		
		// ********************* END OF BIG TIME LOOP ********************* //

		// --- End the line for k-th society in output files --- //
		filew << endl;
		filec << endl;
		fileg << endl;
		fileEf << endl;
	
	
	}

	

	// ******************************************************************* //
	// ********************* CLOSE ALL THE OFSTRESMS ********************* //		
	// ******************************************************************* //

	//filet.close();
	filew.close();
	filec.close();
	fileg.close();
	filer.close();
	fileEf.close();
	return 0;
}
