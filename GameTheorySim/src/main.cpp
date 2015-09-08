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
	string simName("GameTheorySim");
	int control=1;
	bool MemTest = false;

	// Check which kind of simulation I wish to conduct //
	if(cmdOptionExists(argv, argv+argc, "-SML"))
	{
		SML=true;
		simType="SML";

		// Check for memory test //
		if(cmdOptionExists(argv, argv+argc, "-MT"))
		{
			MemTest=true;
		}
	
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

	// Can change the simulation name with commandline parameters if I wish //
	if(cmdOptionExists(argv, argv+argc, "-N"))
	{
		char * tmp = getCmdOption(argv, argv + argc, "-N");
		simName = string(tmp); 
	}


	if (control > 0){
	
		GameTheorySim myGTsim (simName, SML);
		myGTsim.runSimulation_NWA(simType, MemTest);
	
	}

}
