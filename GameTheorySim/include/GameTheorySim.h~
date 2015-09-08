// ------------------------------------------------ //
// -------- GameTheorySim class header file --------//
// ------------------------------------------------ //

#ifndef GAMETHEORYSIM_H_INCLUDED
#define GAMETHEORYSIM_H_INCLUDED


// C++ Std includes 
#include<cstdio>
#include<cstdlib>
#include<iostream>
#include<fstream>
#include<cmath>
#include<math.h>
#include<iomanip>
#include <algorithm>    
#include <vector>       
#include <new>
#include <string>
#include <sstream>

// 3rd party library includes
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h> //Needed for the beta pdf


//#include <boost/multiprecision/cpp_dec_float.hpp>
//#include <boost/multiprecision/gmp.hpp>
#include <boost/multiprecision/mpfr.hpp>



// For convinience use boost::multiprecision as namespace
//using namespace boost::multiprecision;

// *** MAKE CUSTOM MULTIPRECISION TYPE *** //

// Uncomment this for the usual boost cpp multiprecision (only need boost) --> ALSO CHANGE LIB INCLUDE 
//typedef number<cpp_dec_float<256, long long> > myMP_float;

// Umcomment this for the GMP multiprecision --> ALSO CHANGE LIB INCLUDE
//typedef number<gmp_float<100> > myMP_float;

// Uncoment this for the mpfr multiprecision (The fastest and most recommended by Boost and GNU) --> ALSO CHANGE LIB INCLUDE
typedef boost::multiprecision::number< boost::multiprecision::backends::mpfr_float_backend<100> > myMP_float;


// My includes 
//#include "mpEXP.h" 
#include "timeEngine.h"


using namespace std;

////////////////////////////////////////////////
// ----- GameTheorySim Class definition ----- // 
////////////////////////////////////////////////
class GameTheorySim {

	
	public:

	// Disallow copying and assignement operations
	GameTheorySim(const GameTheorySim&) = delete;
    	GameTheorySim& operator=(const GameTheorySim&) = delete; 

	// --- Constructor and destructor --- //
	GameTheorySim(string SimName, bool SML);	// IF SML == TRUE -->> Simple memory learning simulation ELSE Nash eq.
	//GameTheorySim();	
	~GameTheorySim();

	// The society initialization method for simple memory learning (SML) schmea //
	void initSociety_SML_NWA(int k, bool MemTest);

	// The society initialization method for Nash Equlibrium schmea //
	void initSociety_NashEQ_NWA(int k);

	// Run simulation for one full society for SML //
	void runSociety_SML_NWA(int k, bool MemTest);

	// Run simulation for one full society for NashEQ //
	void runSociety_NashEQ_NWA(int k);

	// Run Full simulation: either SML or NashEQ //
	void runSimulation_NWA(string type, bool MemTest);




	private:
	
	// ----------------------------------------------------- //
	// ----- The fields and varibles of the simulation ----- // 
	// ----------------------------------------------------- //

	// RNG Stuff 	
	gsl_rng *gslpointer; 	// Pointer to the type of rng
	FILE *pfile; 		// File to read from /usr/urandom
	unsigned int seed; 	// Seed of the random number generator

	// TimeEngine for timestamp
	TimeEngine *myTimeEngine;

	// Simulation size variables //
	int NE;			// Number of "societies" in the ensemble
	int N;			// Number of agents in one "society" 
	int M; 			// Number of groups
	int S;			// Size of groups

	// Arrays describing the agents //
	double *r;		// Array to hold INVESTMENT TALENT of each agent
	double *w;		// Array to hold WEALTH of each agent
	double *w0;		// Array to hold the INVESTMENT CAP of each agent
	double *xi;		// Array to hold the LEARNIGN SKILL/TALENT of each agent
	double *alpha;	 	// Array to hold contribution coefficient of each agent (contribution == strategy == co-op)
	double *oldalpha;	// Array to hold the contribution coefficient (strategy) for each agent from the last timestep
	double *OLDw;		// Array to hold OLDwealth of each agent
	int *rank;		// Ranking array: shows in which order the agents are. 
				// IMPORTANT CONVENTION: The HIGHER the ORDER of agent, the BETTER the PLACEMENT !!!
	double *O;		// Array to hold the output for each turn !!! THIS DOES NOT NEED MORE THEN N ENTRIES !!!


	// Variables for the INVESTMENT TALENT distribution //
	double muR; 		// Mean of the Gaussian distribution for the INVESTMENT TALENT --> If mean is negtive: EUQAL INVESTMENT TALENT
	double sigmagR;		// Variance sigma of the Gaussian distribution for the talent

	// Variables for the INVESTMENT CAP distribution //
	double muW0; 		// Mean of the Gaussian distribution for the INVESTMENT CAP --> If mean is negtive: EUQAL INVESTMENT CAP
	double sigmagW0;	// Variance sigma of the Gaussian distribution for the talent


	// Variables for the LEARNING SKILL/TALENT distribution //
	double muXI; 		// Mean of the Gaussian distribution for the LEARNING SKILL/TALENT --> If mean is negtive: EUQAL LEARNING SKILL/TALENT
	double sigmagXI;	// Variance sigma of the Gaussian distribution for the talent


	// Simulation settings 
	int Choice;		// The choice of grouping schema
	int monincr; 		// Indicates in how many intervals I want to divide the interval [0,1] for the contributions.
				// The number of strategies available to each player is this number +1
        double beta; 		// This is the noise parameter in the logit
	double Q; 		// Return rate
	int T;			// Time at which I want to stop my simulation
	

	// Variables for memory of the learning schema //
	int N_alpha; 					// The number of different strategies possible
	double alpha_interval;				// The interval lenght of the contribution discritasation
	double *memory; 				// The array for the memory

	// Supporting variables for different purposes //
	int t;
	int i; 			// Support variable 
	double eff;		// Holds the production each turn ('normalized')
	double wealth;  	// Wealth is the total wealth at every time step for each society
	double oldwealth;	// Oldwealth is a support variable I need to compute the efficiency
	double totInvestCap; 	// The total investment cap
		

	// ofstream for outputing and output files //
	ofstream filep;		// Parameter output stream 
	ofstream filew;		// Wealth output stream
	ofstream filec;		// Cooperation output stream
	ofstream filer;		// Investment Talent output stream
	ofstream fileW0;	// Investment Cap output stream
	ofstream fileXI;	// Learning Skill/Talent output stream
	ofstream fileg;		// Gini coef. output stream
	ofstream fileEf;	// Efficientcy output stream
	ofstream MemOut;	// Memory testing output stream

	const string filenamep = string ("parameters.txt");		// File to hold the simulation input parameters	
	const string filenamew = string ("wealth.txt"); 		// File to hold the total wealth for each ensemble for each timestep
	const string filenamec = string ("cooperation.txt");		// File to hold the average co-op in each "society" for each timestep
	const string filenamer = string ("talent.txt");			// File to hold the talent for each individual in each society (rows: societies, columns: individuals)
	const string filenameW0 = string ("investment_cap.txt");	// File to hold the investment cap for each individual in each society (rows: societies, columns: individuals)
	const string filenameXI = string ("learning_skill.txt");	// File to hold the learning skill/talent for each individual in each society (rows: societies, columns: individuals)
	const string filenameg = string ("giniCoef.txt");		// File to hold the Gini coef. for each society on each time step.
	const string filenameEf = string ("efficiency.txt");		// File to hold the efficientcy for each society on each time step.
	const string FNMemOut = string ("memory_test.txt");		// File to hold the memory testing output at each timestep.
	



	// ----- Class methods necessary for the simulation ----- //


	// --- The main methods for the simulation --- //
	
	// - Only NO Wealth Accumulation (NWA) methods at the moment - //

	// The makeorder mehtod, which ALSO fills the OUTPUT array //
	void makeorder_NWA(double *alpha, int *rank, double *O);		

	// The GROUPING methods, which ALSO claculate the UTILITY //
	void formgroups_NWA(double *w, double *alpha, double *O, int *rank);			

	// Methods for calculating the utility of q-th agent without forming groups //
	double getutility_NWA(int q, double *w, double *alpha, double *O, int *rank);	

	// Methods for gettign a new strategy for the k-th player //
	double NewStrategyNashEQ_NWA(int k);						// Nash equlibrium
	double newStrategySML_NWA(int k); 						// Simple learning schema 

	// Memory update methods for the simple learning schema //
	void memoryUpdate_NWA(); 


	

	// --- Support functions --- //

	// --- Binary search --- //
	int binaryprobsearch(double *Gamma, int M, double x);

	// Method to return the rank of the agent k. It's the inverse of rank //
	int inverserank(int k, int *rank);

	// Method for computing the Gini Coef. as defined in Wikipedia //
	double computegini(double *w);

	// Method for summing over vector of length n //
	double sumofvector(double *a, int n);

	// Method for calculating exponential function using Boost Multiprecision lib
	myMP_float my_mpEXP(double x);



	// --- Printing/Outputting functions --- //

	// Ensemble print -->> The standart output //
	void ensemblePrint();

	// ensemblePrintParams -->> prints the parameters.txt file //
	void ensemblePrintParams();

	// ensemblePrintTalent -->> Prints the INVESTMENT TALENT to the talent.txt file //
	void ensemblePrintTalent();

	// ensemblePrintW0 -->> Prints the INVESTMENT CAP to the investment_cap.txt file //
	void ensemblePrintW0();

	// ensemblePrintXI -->> Prints the LEARNING SKILL/TALENT to the learning_skill.txt file //
	void ensemblePrintXI();

	// printMemory -->> Prints the memory of the learning schema || for testing purposes -->> have to give a opened ofstream //
	void printMemory(int t);






}; 

#endif // GAMETHEORYSIM_H_INCLUDED
