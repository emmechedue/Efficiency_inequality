// --------------------------------------------------------------------------- //
// -------- GameTheorySim derived class --> GTS_quadPayoff header file --------//
// --------------------------------------------------------------------------- //

#ifndef GTS_QUADPAYOFF_H_INCLUDED
#define GTS_QUADPAYOFF_H_INCLUDED


// --- Include the parent class GameTheorySim --- //
#include "GameTheorySim.h"


using namespace std;



///////////////////////////////////////////////////////////////////////////
// ----- GameTheorySim derived class --> GTS_quadPayoff definition ----- // 
///////////////////////////////////////////////////////////////////////////
class GTS_quadPayoff : public GameTheorySim {

	public:

	
	// --- Constructor and destructor --- //
	GTS_quadPayoff(string SimName, bool SML);	// IF SML == TRUE -->> Simple memory learning simulation ELSE Nash eq.
	//GameTheorySim();	
	//~GTS_quadPayoff();

	// ----------------------------------------------------------------------------------------------- //
	// --- The Payoff function -->> this is where the derived class is different from the original --- //
	// ----------------------------------------------------------------------------------------------- //

	double payoff(double w, double w0, double alpha, double gain);


	// -------------------------------------------------------------- //
	// --- Output file -->> Slight difference on the header files --- //
	// -------------------------------------------------------------- //

	void initOutput(string SimName, bool SML);

	// ----------------------------------------------------------------------------------------------- //

/*

	// --- The method, which don't differ from the original GTS class --- //
	
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


	// The makeorder mehtod, which ALSO fills the OUTPUT array //
	void makeorder_NWA(double *alpha, int *rank, double *O);		

	// The GROUPING methods, which ALSO claculate the UTILITY //
	void formgroups_NWA(double *w, double *alpha, double *O, int *rank);	
	//void formgroups_NLP_NWA(double *w, double *alpha, double *O, int *rank);		// Non-linear payoff schema

	// Methods for calculating the utility of q-th agent without forming groups //
	double getutility_NWA(int q, double *w, double *alpha, double *O, int *rank);
	//double getutility_NLP_NWA(int q, double *w, double *alpha, double *O, int *rank);	// Non-linear payoff schema	


	// Methods for gettign a new strategy for the k-th player //
	double NewStrategyNashEQ_NWA(int k);						// Nash equlibrium
	//double NewStrategyNashEQ_NWA(int k);						// Nash equlibrium
	double newStrategySML_NWA(int k); 						// Simple learning schema 

	// Memory update methods for the simple learning schema //
	void memoryUpdate_NWA(); 


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

*/

};


#endif // GTS_QUADPAYOFF_H_INCLUDED
