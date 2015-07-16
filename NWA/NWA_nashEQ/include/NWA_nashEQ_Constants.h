/* 
 * File:   Constants.h
 * Author: sduca
 *
 * Created on February 4, 2015, 2:06 PM
 */

#pragma once
#include<cstdio>
#include<cstdlib>


class Constants{
	public:
	Constants();	//Default constructor and destructor
	~Constants();
	int NE;		// Number of "societies" in the ensemble
	int N;		//Number of agents playing the game
	int T;		//Time at which I want to stop my simulation 
	int S; 		//Size of the groups
	double Q; 	//Return rate
	double mu; 	//Mean of the Gaussian distribution for the talent --> If mean is negtive: EUQAL TALEN
	double sigmag; //Variance sigma of the Gaussian distribution for the talent
	double W0; 	//Initial wealth value
	int Choice; 	//Indicates the choice of what are we using to rank agents in the groups
    int monincr; 	//Indicates in how many intervals I want to divide the interval [0,1] for the contributions. The number of strategies available to each player is this number +1
    double beta; 	//This is the noise parameter in the logit
};
