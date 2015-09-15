// --------------------------------------------------------------------------- //
// -------- GameTheorySim derived class --> GTS_quadPayoff source file --------//
// --------------------------------------------------------------------------- //


#include "GTS_quadPayoff.h"


using namespace std;


// -------------------------------------- //
// ----- Constructor and destructor ----- //
// -------------------------------------- //

// Constructor //
GTS_quadPayoff::GTS_quadPayoff(string SimName, bool SML) : GameTheorySim (SimName, SML) {}



// ----------------------------------------------------------------------------------------------- //
// --- The Payoff function -->> this is where the derived class is different from the original --- //
// ----------------------------------------------------------------------------------------------- //


// --- The payoff function -->> calculates the payoff --- //

double GTS_quadPayoff::payoff(double w, double w0, double alpha, double gain){

	return (w - w0*alpha*alpha + gain);

}



// -------------------------------------------------------------- //
// --- Output file -->> Slight difference on the header files --- //
// -------------------------------------------------------------- //


// Initialization of the output files //
void GTS_quadPayoff::initOutput(string SimName, bool SML){
	// --- Open the output files and write the headers --- //

	// Define the Header template
	string HEAD ("# ");
	HEAD += myTimeEngine->get_timestamp();
	HEAD += "\tOutput of Simulation ";
	HEAD += SimName;
	HEAD += ":\t";

	// Define the parameters
	string PARAM_LIST;
	ostringstream StrStream;   // string stream used for the conversion
	StrStream << "# Parameters:\tNE=" << NE << " N=" << N << " T=" << T  << " S=" << S;
	StrStream << " Q=" << Q << " monincr=" << monincr << " beta=" << beta << " choice=" << Choice; 
	StrStream << " muR="<< muR << " sigmagR=" << sigmagR << " muW0="<< muW0 << " sigmagW0=" << sigmagW0 ;
	StrStream << " muXI="<< muXI << " sigmagXI=" << sigmagXI <<  " seed=" << seed;
	if(SML){
		StrStream << " type=SML/quadraticPayoff";
	}
	else {
		StrStream << " type=NashEq/quadraticPayoff";
	}
	PARAM_LIST = StrStream.str();

	// The WEALTH file 	
	filew.open(filenamew,ios::out|ios::trunc); 
	if(filew.is_open()){
		filew << HEAD << "Total Wealth" << endl;
		filew << PARAM_LIST << endl; 
		filew << "# Total wealth -> Rows: 'societies' \tColumns: timesteps" << endl;	
		filew << "#" << endl;
	}

	// The EFFICIENCY file 	
	fileEf.open(filenameEf,ios::out|ios::trunc); 
	if(fileEf.is_open()){
		fileEf << HEAD << "Efficiency" << endl;
		fileEf << PARAM_LIST << endl;
		fileEf << "# Efficiency -> Rows: 'societies' \tColumns: timesteps" << endl;	
		fileEf << "#" << endl;
	}

	// The COOPERATION file
	filec.open(filenamec,ios::out|ios::trunc);
	if(filec.is_open()){ 
		filec << HEAD << "Average Co-op (strategy)" << endl;
		filec << PARAM_LIST << endl;	
		filec << "# Average Co-op (strategy) -> Rows: 'societies' \tColumns: timesteps" << endl;		
		filec << "#" << endl;
	}

	// The GINI COEF. file
	fileg.open(filenameg,ios::out|ios::trunc);
	if(filec.is_open()){ 
		fileg << HEAD << "Gini coef." << endl;
		fileg << PARAM_LIST << endl;
		fileg << "# Gini coef -> Rows: 'societies' \tColumns: timesteps" << endl;	
		fileg << "#" << endl;
	}

	// The INVESTMENT TALENT file
	filer.open(filenamer,ios::out|ios::trunc); //Open the talent file
	if(filer.is_open()){
		filer << HEAD << "Investment Talent of each player in 'society'" << endl;
		filer << PARAM_LIST << endl;
		filer << "# Rows: 'societies'\tColumns: individuals/players" << endl;
		filer << "#" << endl;
	}

	// The INVESTMENT CAP file
	fileW0.open(filenameW0,ios::out|ios::trunc); //Open the talent file
	if(fileW0.is_open()){
		fileW0 << HEAD << "Investment Cap of each player in 'society'" << endl;
		fileW0 << PARAM_LIST << endl;
		fileW0 << "# Rows: 'societies'\tColumns: individuals/players" << endl;
		fileW0 << "#" << endl;
	}

	// The LEARNING SKILL/TALENT file
	if(SML){
		fileXI.open(filenameXI,ios::out|ios::trunc); //Open the talent file
		if(fileXI.is_open()){
			fileXI << HEAD << "Learning skill/talent of each player in 'society'" << endl;
			fileXI << PARAM_LIST << endl;
			fileXI << "# Rows: 'societies'\tColumns: individuals/players" << endl;
			fileXI << "#" << endl;
		}

		MemOut.open(FNMemOut,ios::out|ios::trunc); //Open the memory test file
		if(MemOut.is_open()){
			MemOut << HEAD << "Memory of each player in 'society'" << endl;
			MemOut << PARAM_LIST << endl;
			MemOut << "# Rows: 'Players'\tColumns: strategies" << endl;
			MemOut << "#" << endl;
		}
	}


}


/*
// ----------------------------------------------------------------------------------------------- //





// --- The method, which don't differ from the original GTS class --- //
	
	// The society initialization method for simple memory learning (SML) schmea //
	void GTS_quadPayoff::initSociety_SML_NWA(int k, bool MemTest) { GameTheorySim::initSociety_SML_NWA(k, MemTest); } 

	// The society initialization method for Nash Equlibrium schmea //
	void GTS_quadPayoff::initSociety_NashEQ_NWA(int k) { GameTheorySim::initSociety_NashEQ_NWA(k); } 

	// Run simulation for one full society for SML //
	void GTS_quadPayoff::runSociety_SML_NWA(int k, bool MemTest) { GameTheorySim::runSociety_SML_NWA(k, MemTest); } 

	// Run simulation for one full society for NashEQ //
	void GTS_quadPayoff::runSociety_NashEQ_NWA(int k) { GameTheorySim::runSociety_NashEQ_NWA(k); } 

	// Run Full simulation: either SML or NashEQ //
	void GTS_quadPayoff::runSimulation_NWA(string type, bool MemTest) { GameTheorySim::runSimulation_NWA(type, MemTest); } 


	// The makeorder mehtod, which ALSO fills the OUTPUT array //
	void GTS_quadPayoff::makeorder_NWA(double *alpha, int *rank, double *O) { GameTheorySim::makeorder_NWA(alpha, rank, O); } 		

	// The GROUPING methods, which ALSO claculate the UTILITY //
	void GTS_quadPayoff::formgroups_NWA(double *w, double *alpha, double *O, int *rank) { GameTheorySim::formgroups_NWA(w, alpha, O, rank); } 	

	// Methods for calculating the utility of q-th agent without forming groups //
	double GTS_quadPayoff::getutility_NWA(int q, double *w, double *alpha, double *O, int *rank) { return GameTheorySim::getutility_NWA(q, w, alpha, O, rank); } 


	// Methods for gettign a new strategy for the k-th player //
	double GTS_quadPayoff::NewStrategyNashEQ_NWA(int k) { return GameTheorySim::NewStrategyNashEQ_NWA(k); } 			// Nash equlibrium
	double GTS_quadPayoff::newStrategySML_NWA(int k) { return GameTheorySim::newStrategySML_NWA(k); }  			// Simple learning schema 

	// Memory update methods for the simple learning schema //
	void GTS_quadPayoff::memoryUpdate_NWA() { GameTheorySim::memoryUpdate_NWA(); } 




	// --- Printing/Outputting functions --- //

	// Ensemble print -->> The standart output //
	void ensemblePrint() { GameTheorySim::ensemblePrint() ;}

	// ensemblePrintParams -->> prints the parameters.txt file //
	void ensemblePrintParams(){ GameTheorySim::ensemblePrintParams() ;}

	// ensemblePrintTalent -->> Prints the INVESTMENT TALENT to the talent.txt file //
	void ensemblePrintTalent(){ GameTheorySim::ensemblePrintTalent() ;}

	// ensemblePrintW0 -->> Prints the INVESTMENT CAP to the investment_cap.txt file //
	void ensemblePrintW0(){ GameTheorySim::ensemblePrintW0() ;}

	// ensemblePrintXI -->> Prints the LEARNING SKILL/TALENT to the learning_skill.txt file //
	void ensemblePrintXI(){ GameTheorySim::ensemblePrintXI() ;}

	// printMemory -->> Prints the memory of the learning schema || for testing purposes -->> have to give a opened ofstream //
	void printMemory(int t){ GameTheorySim::printMemory(t) ;}

*/

















