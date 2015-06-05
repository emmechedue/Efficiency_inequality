/* 
 * File:   Constants.cpp
 * Author: sduca
 * 
 * Created on February 4, 2015, 2:06 PM
 */


#include<cstdio>
#include<cstdlib>
#include<fstream>
#include "ensemble_Constants.h"

using namespace std;
/*

Constants::Constants(){
	NE= 300; // Number of "societies" in the ensemble
	N= 100; //Number of agents playing the game
	T= 35;//Time at which I want to stop my simulation
	S= 4; //Size of the groups
	Q= 0.3; //Return rate
	mu= 1; //Mean of the Gaussian distribution for the talent
	sigmag= 0.22; //Variance sigma of the Gaussian distribution for the talent
	W0= 2; //Initial wealth value
	Choice= 1; //Indicates the choice of what are we using to rank agents in the groups
    monincr= 100; //Indicates in how many intervals I want to divide the interval [0,1] for the contributions. The number of strategies available to each player is this number +1
    beta= 2; //This is the noise parameter in the logit
}
Constants::~Constants(){}

*/

//This reads from a file:

Constants::Constants(){ //Note that name must be the entire path; i.e. "./config.conf"
	char line[256];
	int linenum=0;
	int count=0, M=11; //M is the amount of parameters I have to give, count will range from 0 to M-1
	double vector[M]; //will store the M parameters
	FILE *pfile;

	pfile = fopen ("./config.conf" , "r");
	//pfile= fopen("/project/theorie/s/Stefano.Duca/Analysis/Prog/config.conf", "r"); //Here I have to put the folder where the config file will be!

	while(fgets(line, 256, pfile) != NULL)
	{
		linenum++;
		if(line[0] == '#') {continue;} //I'm going to the next line without reading and incrementing count
		sscanf(line, "%*s %lf", &vector[count]);
		count ++;
	}
	
	//Here I cast when I actually have integers
	NE=int(vector[0]); // Number of "societies" in the ensemble
	N=int(vector[1]); //Number of agents playing the game
	T=vector[2];//Time at which I want to stop my simulation
	S=int(vector[3]); //Size of the groups
	Q=vector[4]; //Return rate
	mu=vector[5]; //Mean of the Gaussian distribution for the talent
	sigmag=vector[6]; //Variance sigma of the Gaussian distribution for the talent
	W0=vector[7]; //Initial wealth value
	Choice=int(vector[8]); //Indicates the choice of what are we using to rank agents in the groups
    monincr=int(vector[9]); //Indicates in how many intervals I want to divide the interval [0,1] for the contributions. The number of strategies available to each player is this number +1
    beta= vector[10]; //This is the noise parameter in the logit
}
Constants::~Constants(){}
