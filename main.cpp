/* 
 * File:   main.cpp
 * Author: sduca
 *
 * Created on February 4, 2015, 4:51 PM
 */

#include<cstdio>
#include<cstdlib>
#include<fstream>
#include<math.h>
#include<iostream>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h> //Needed for the beta pdf
#include "Constants.h"
#include "Fcts.h"
#include <new>
#include <algorithm>    
#include <vector> 


using namespace std;

int main() {
    gsl_rng *gslpointer; //Pointer to the type of rng
    FILE *pfile; //file to read from /usr/urandom
    unsigned int seed; //Seed of the random number generator
    Constants cons; //The class containing all the constants
    int M; //Number of groups
    double oldt, t; //The time and a i variable to keep track of time
    int enne; //int variable that I need in a time thing
    int i; //Support variable
    double r[cons.N], w[cons.N], alpha[cons.N]; //talent, wealth and contribution for each agent
    int rank[cons.N]; //Ranking vector. It shows in which order the agents are. IMPORTANT CONVENTION: The higher the order, the better the placement!!
    double O[cons.N]; //The output vector at each turn
    ofstream filep, filet,filew, filec;//Output files. Filep is for the parameters file, filet for the main output one and filew for the wealth
    const char filenamew[]="wealth.txt"; //Tese are all the output files that I need to print
    const char filenamep[]="parameters.txt";
    const char filenamet[]="time.txt";
    const char filenamec[]="cooperation.txt";
    double ranvar; //random variable
    int k; //agent that at that particular time is changing its mind
    double eff, oldeff; //These are to print the efficiency and then renormalize the wealth. Here the efficiency is defined as percentage  growth
    int counteff; //This I need to make proper averages of the efficiency

    /***************** HERE I INITIALIZE THE RANDOM NUMBER GENERATOR ***********************/
	
    pfile = fopen ("/dev/urandom", "r");
    i=fread (&seed, sizeof (seed), 1, pfile); //I added the rand= ... just to not be bothered anymore by the warnings!
    fclose(pfile);
    gslpointer = gsl_rng_alloc(gsl_rng_mt19937); //I'm using the "Mersenne Twister" generator!
    gsl_rng_set(gslpointer,seed); // Starting the generator
    
    srand(seed); //I also start the c random number generator. It's unlikely, but I might need it to break ties in the ordering
    /***************************************************************************************************/
    
    /**************************** HERE I PRINT OUT THE HEADER OF THE OUTPUT FILES ******************************************/
    filew.open(filenamew,ios::out|ios::trunc); //Open the wealth file
    filew<<"#Wealth at each time step for the simulation Efficiency_inequality with:"<<endl;
    filew<<"#N="<<cons.N<<" T="<<cons.T<<" interval="<<cons.interval<<" S="<<cons.S<<" Q="<<cons.Q<<" lambda="<<cons.lambda<<" mu="<<cons.mu<<" sigmag="<<cons.sigmag<<" W0="<<cons.W0<<" choice="<<cons.Choice<<" monincr="<<cons.monincr<<" beta="<<cons.beta<<" seed="<<seed<<endl;
    
    filec.open(filenamec,ios::out|ios::trunc); //Open the wealth file
    filec<<"#Strategy at each time step for the simulation Efficiency_inequality with:"<<endl;
    filec<<"#N="<<cons.N<<" T="<<cons.T<<" interval="<<cons.interval<<" S="<<cons.S<<" Q="<<cons.Q<<" lambda="<<cons.lambda<<" mu="<<cons.mu<<" sigmag="<<cons.sigmag<<" W0="<<cons.W0<<" choice="<<cons.Choice<<" monincr="<<cons.monincr<<" beta="<<cons.beta<<" seed="<<seed<<endl;
    
    filet.open(filenamet,ios::out|ios::trunc); //Open the time file
    filet<<"#Results for the simulation Efficiency_inequality with:"<<endl;
    filet<<"#N="<<cons.N<<" T="<<cons.T<<" interval="<<cons.interval<<" S="<<cons.S<<" Q="<<cons.Q<<" lambda="<<cons.lambda<<" mu="<<cons.mu<<" sigmag="<<cons.sigmag<<" W0="<<cons.W0<<" choice="<<cons.Choice<<" monincr="<<cons.monincr<<" beta="<<cons.beta<<" seed="<<seed<<endl;
    filet<<"#This is in the form of t, Gini coefficient, growth percentage and average cooperation"<<endl;
    /*****************************************************************************************************************/
    
    /***************** START GENERATING ALL THE STUFF ***********************************************************/
    if((cons.N%cons.S) != 0){ //I check that N is exactly divisible by S
    	cout<<"The number of agents is not exactly divisible by the number of groups"<<endl;
    	exit(1);
    }
    else{M=cons.N/cons.S;} //I compute the number of groups defined as N/S
    
    /* FILL THE VECTORS OF TALENT , WEALTH AND RANKING AND INITIALIZE EVERYTHING*/
    t=0;
    eff=0;
    oldeff=0; //I don't really need to initialize this, but it's better anyway
    for(i=0;i<cons.N;i++){
        w[i]=cons.W0;
	rank[i]=i; 
	r[i]=gsl_ran_gaussian(gslpointer,cons.sigmag)+cons.mu ; //I have to add the mean of the gaussian because the generator has mean zero
        alpha[i]=0; //Initially, I set all the strategies to full defection
    }
    renormalizearray(w,cons.N);

    // Now I fill up the output files at t=0
    printstuffsingleloop(filet, filew, filec, t, eff, w,alpha, cons); //Print time, gini coefficient, efficiency and all the wealth for each time step
    counteff = -1; // This here is to reset the counter for the average of the efficiency

    /***************** END OF INITIALIZATION  **********************************************/
    
    /************************************** HERE THE BIG TIME LOOP STARTS ********************************/
    do {
        
        /*Here I sample and update the time*/
        ranvar=gsl_ran_exponential(gslpointer,cons.lambda);
        if (t + ranvar > cons.T) { //If the new time would be bigger than T, I simply print out everything I have to print and then break the loop
            oldt = cons.T - t; //Note that here oldt has a different use, I'm simply using this variable since I don't need it anymore
            enne = floor(oldt / cons.interval);
            for (i = 1; i < (enne + 1); i++) { //Note that here I have less of enne + 1 and not <=, since here I need the time explicitely! The enne+1 is because I need it to multply. Meaning that all this is because I start from 1. And I haveo to do it since I still have to update the time t
                printstuffsingleloop(filet, filew, filec, t + i * cons.interval, eff, w,alpha, cons); //Print time, Gini coefficient, efficiency and all the wealth for each time step
                counteff = -1; // This here is to reset the counter for the average of the efficiency
            }
            printstuffsingleloop(filet, filew,filec, cons.T, eff, w, alpha, cons); //Print time, Gini coefficient, efficiency and all the wealth for each time step
            counteff = -1; // This here is to reset the counter for the average of the efficiency
            break; //Exit from the do loop
        }
        if (ranvar > cons.interval) { //Here is to check if I have to reprint the old situation before update the system!
            enne = floor(ranvar / cons.interval);
            for (i = 1; i <= enne; i++) { //<=because I start from 1
                printstuffsingleloop(filet, filew,filec, t + i * cons.interval, eff, w, alpha, cons); //Print time, Gini coefficient, efficiency and all the wealth for each time step
                counteff = -1; // This here is to reset the counter for the average of the efficiency
                cout << "The time is " << t + i * cons.interval << endl; //Just to check
            }
            ranvar = ranvar - cons.interval*enne;
            t = t + enne * cons.interval;
        }
        t = t + ranvar; //Update the time. I have to do it after I do the check for the time
        oldt = oldt + ranvar; //Update oldt
        
        k=gsl_rng_uniform_int(gslpointer,cons.N); //Here I uniformly at random select an agent
        updatestrategy(k,alpha,w,r,cons, gslpointer); // Here I update the strategy of agent k according to a logit distribution.
        //alpha[k]=gsl_ran_beta(gslpointer,2,2); //Here I would update the strategy completely at random!
        makeorder(r,w,alpha,rank,cons,O); //I generate the ranking and compute the output for each player
        formgroups(w, alpha, O,rank, cons, M); // Here I use the ranking to generate the groups and update the wealth of the players
        eff = sumofvector(w, cons.N) - 1. ; //Here I compute the total amount of wealth. Since the total wealth was normalized to one before, this is the increase in percentage of wealth
        counteff++; //These 3 lines here are to compute proper averages of eff. This is in case I perform many iterations without printing
        eff = updateaverages(oldeff, eff, counteff);
        oldeff = eff;
        renormalizearray(w,cons.N); //Finally, Here I renormalize the array.
        
        if(oldt>=cons.interval){ //Checks whether I have to print or not
           
            printstuffsingleloop(filet, filew,filec, t, eff, w, alpha,cons); //Print time, Gini coefficient, efficiency and all the wealth for each time step
            counteff = -1; // This here is to reset the counter for the average of the efficiency 
            oldt=oldt -cons.interval; //Subtract by oldt the value of interval to start counting again
            cout<<"The time is "<<t<<endl; //Just to check
        }
		
    }while(t<=cons.T);
    /*****************************************************************************************************/
    /***********Now I just print all the parameters to a file! Then I close everything */
    filep.open(filenamep,ios::out|ios::trunc);
    printparamsingleloop(filep,cons);
    filep.close();
    filet.close();
    filew.close();
    filec.close();
    return 0;
}