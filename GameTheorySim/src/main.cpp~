/* 
 * File:   main.cpp
 * Author: sduca
 *
 * Created on February 4, 2015, 4:51 PM
 */

#include "GameTheorySim.h"
#include "readcmd.h"


using namespace std;

int main (int argc, char* argv[]){

	// Define variables
	bool SML;
	string simType;
	int control=1;

	// Check which kind of simulation I wish to conduct //
	if(cmdOptionExists(argv, argv+argc, "-SML"))
	{
		SML=true;
		simType="SML";
	}
	else if(cmdOptionExists(argv, argv+argc, "-NashEQ"))
	{
		SML=false;
		simType="NashEQ";
	}
	else{
		cout << "Incorrect simulation type! Please choose either: -SML or -NashEQ" << endl;
		control=-1;
		exit(1);
	}


	if (control > 0){
	
		GameTheorySim myGTsim ("GameTheorySim Test", SML);
		myGTsim.runSimulation_NWA(simType);
	
	}

}
