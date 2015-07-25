// ------------------------------------------------ //
// -------- GameTheorySim class source file --------//
// ------------------------------------------------ //

#include "GameTheorySim.h"

using namespace std;


// -------------------------------------- //
// ----- Constructor and destructor ----- //
// -------------------------------------- //

// Constructor //
GameTheorySim::GameTheorySim(string SimName, bool SML){
//GameTheorySim::GameTheorySim(){

	// --- Read input parameters from config.conf --- //

	// Define parameters for parsing the input file
	char line[256];
	int linenum=0;
	int NP=14; 		// NP --> amount of parameters I have to give 
	int count=0;		// count will range from 0 to M-1
	double ParamVec[NP]; 	//will store the M parameters
	FILE *pfile;

	// Open the config.conf file
	pfile = fopen ("./config.conf" , "r");
	
	// Read the parameter into ParamVec
	while(fgets(line, 256, pfile) != NULL)
	{
		linenum++;
		if(line[0] == '#') {continue;} //I'm going to the next line without reading and incrementing count
		sscanf(line, "%*s %lf", &ParamVec[count]);
		count ++;
	}

	// Close the config.conf file
	fclose (pfile);
	
	// Read parameter from the ParamVec
	NE = int(ParamVec[0]); 		// Number of "societies" in the ensemble
	N = int(ParamVec[1]); 		// Number of agents in a "society"
	T = ParamVec[2];		// Time at which I want to stop my simulation (Simulation lenght)
	S = int(ParamVec[3]); 		// Size of the groups
	Q = ParamVec[4]; 		// Return rate	
	monincr = int(ParamVec[5]); 	// Indicates in how many intervals I want to divide the interval [0,1] for the contributions. 
					// The number of strategies available to each player is this number +1
        beta = ParamVec[6]; 		// This is the noise parameter in the logit
	Choice = int(ParamVec[7]); 	// Indicates the choice of what are we using to rank agents in the groups
	muR = ParamVec[8]; 		// Mean of the Gaussian distribution for the INVESTMENT TALENT --> Negative == const. as abs(muR)
	sigmagR = ParamVec[9]; 		// Variance sigma of the Gaussian distribution for the INVESTMENT TALENT
	muW0 = ParamVec[10]; 		// Mean of the Gaussian distribution for the IINVESTMENT CAP --> Negative == const. as abs(muW0)
	sigmagW0 = ParamVec[11]; 	// Variance sigma of the Gaussian distribution for the INVESTMENT CAP
	muXI = ParamVec[12]; 		// Mean of the Gaussian distribution for the LEARNING SKILL/TALENT --> Negative == const. as abs(muXI)
	sigmagXI = ParamVec[13]; 	// Variance sigma of the Gaussian distribution for the LEARNING SKILL/TALENT



	// --- Get the number of GROUPS --- //
	// Check that number of agents N is exactly divisible by size of groups S 
	if((N%S) == 0){ 
		// IF is excatly devisable --> Compute number of groups M = N/S
		M=N/S;		
	}
	else{
		// IF NOT excatly devisable -->	Give and ERROR and exit
		cout<<"The number of agents is not exactly divisible by the size of groups"<<endl;
		exit(1);
	} 

	//cout << "Testing class contructor: M = " << M << " N = " << N << " S = " << S << endl;
	
	// --- Now dymanically initializing the arrays for describing the agents --- //
	r = new double[N];		// Array to hold INVESTMENT TALENT of each agent
	w = new double[N];		// Array to hold WEALTH of each agent
	w0 = new double[N];		// Array to hold the INVESTMENT CAP of each agent
	xi = new double[N];		// Array to hold the LEARNIGN SKILL/TALENT of each agent
	alpha = new double[N];	 	// Array to hold contribution coefficient of each agent (contribution == strategy == co-op)
	oldalpha = new double[N];	// Array to hold the contribution coefficient (strategy) for each agent from the last timestep
	OLDw = new double[N];		// Array to hold OLDwealth of each agent
	rank = new int[N];		// Ranking array: shows in which order the agents are. 
					// IMPORTANT CONVENTION: The HIGHER the ORDER of agent, the BETTER the PLACEMENT !!!
	O = new double[N];		// Array to hold the output for each turn !!! THIS DOES NOT NEED MORE THEN N ENTRIES !!!




	// --- Initialising the GSL Mersenne twister RNG & a C srand RNG --- //
	// GSL Mersenne Twister 19937 RNG
	pfile = fopen ("/dev/urandom", "r");
	fread (&seed, sizeof (seed), 1, pfile);		// I added the rand= ... just to not be bothered anymore by the warnings!
	fclose(pfile);
	gslpointer = gsl_rng_alloc(gsl_rng_mt19937); 	// I'm using the "Mersenne Twister" generator!
	gsl_rng_set(gslpointer,seed); 			// Starting the generator
	
	// C srand RNG
	srand(seed); 					// It's unlikely, but I might need it to break ties in the ordering


	// --- Initialise the memory and stuff in connection with it --- //
	N_alpha = monincr + 1; 		// The number of different strategies possible
	alpha_interval = 1.0/monincr;		// The interval lenght of the contribution discritasation
	memory = new double [N*N_alpha]; 	// The array for the memory 


	// --- Initialise the TimeEngine --- //
	myTimeEngine = new TimeEngine();


	// --- Initialise the ofstreams --- //
	// ofstream for outputing and output files //
/*	filep = new ofstream(filenamep,ios::out|ios::trunc));		// Parameter output stream 
	filew = new ofstream(filenamew,ios::out|ios::trunc));		// Wealth output stream
	filec= new ofstream(filenamec,ios::out|ios::trunc));		// Cooperation output stream
	filer= new ofstream(filenamer,ios::out|ios::trunc));		// Investment Talent output stream
	fileW0= new ofstream(filenameW0,ios::out|ios::trunc));	// Investment Cap output stream
	fileXI= new ofstream(filenameXI,ios::out|ios::trunc));	// Learning Skill/Talent output stream
	fileg= new ofstream(filenameg,ios::out|ios::trunc));		// Gini coef. output stream
	fileEf= new ofstream(filenameEf,ios::out|ios::trunc));	// Efficientcy output stream	
*/

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
		StrStream << " type=SML";
	}
	else {
		StrStream << " type=NashEq";
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
	}


}

// Destructor //
GameTheorySim::~GameTheorySim(){

	// Clean up the dynamically allocated memory
	// --- Now dymanically initializing the arrays for describing the agents --- //
	delete[] r;		// Array to hold INVESTMENT TALENT of each agent
	delete[] w;		// Array to hold WEALTH of each agent
	delete[] w0;		// Array to hold the INVESTMENT CAP of each agent
	delete[] xi;		// Array to hold the LEARNIGN SKILL/TALENT of each agent
	delete[] alpha;	 	// Array to hold contribution coefficient of each agent (contribution == strategy == co-op)
	delete[] oldalpha;	// Array to hold the contribution coefficient (strategy) for each agent from the last timestep
	delete[] OLDw;		// Array to hold OLDwealth of each agent
	delete[] rank;		// Ranking array: shows in which order the agents are. 
					// IMPORTANT CONVENTION: The HIGHER the ORDER of agent, the BETTER the PLACEMENT !!!
	delete[] O;		// Array to hold the output for each turn !!! THIS DOES NOT NEED MORE THEN N ENTRIES !!!
	delete[] memory;


	// Close the output files 
	filew.close();
	filec.close();
	fileg.close();
	fileEf.close();
	filer.close();
	fileW0.close();
	if(fileXI.is_open()){
		fileXI.close();
	}


}


// ------------------------------- //
// ----- Full simulation run ----- //
// ------------------------------- //


// --- Run Full simulation: either SML or NashEQ --- //
void GameTheorySim::runSimulation_NWA(string type){

	// Check which type of simulation I wish to conduct. 	
	if (type == "SML"){

		// Loop over all of the societies
		for (int scn = 0; scn < NE; scn++){

			// Print parameters
			ensemblePrintParams();

			// Initialise society
			initSociety_SML_NWA(scn);

			// Run the simulation for that society
			runSociety_SML_NWA(scn);
		}
	}
	else if (type == "NashEQ"){

		// Loop over all of the societies
		for (int scn = 0; scn < NE; scn++){

			// Print parameters
			ensemblePrintParams();

			// Initialise society
			initSociety_NashEQ_NWA(scn);

			// Run the simulation for that society
			runSociety_NashEQ_NWA(scn);
		}

	}
	else {
			cout << "Incorrect simulation type! Please choose either: SML or NashEQ" << endl;
			exit(1);
	}

	
}



// ---------------------------------------- //
// ----- Initialisations of Societies ----- //
// ---------------------------------------- //


// --- The society initialization method for simple memory learning (SML) schmea --- //
void GameTheorySim::initSociety_SML_NWA(int k){

	// Print for testing
	cout << "Initialising SML society nr " << k << endl;
	
	// FILL THE VECTORS OF TALENT , WEALTH AND RANKING AND INITIALIZE EVERYTHING

	// Initial timestep is 0
	t=0;			

	// Efficiency at t=0 is 0
	eff = 0.0;

			
	// The arrays describing agents
	for(int i=0; i<N; i++){
		

		// Intial rank is efectivly random
		rank[i]=i;		
				
		// INVESTMENT TALENT //
		// IF mu positive -->> Gussian distributed 
		if(muR > 0 ){  
  			r[i]=gsl_ran_gaussian(gslpointer,sigmagR) + muR;
		}
		// ELSE all eqaul to abs(mu)
		else{ 
			r[i] = abs(muR);
		}

		// INVESTMENT CAP & INITIAL WEALTH //
		// IF mu positive -->> Gussian distributed 
		if(muW0 > 0 ){  
  			w0[i]=gsl_ran_gaussian(gslpointer,sigmagW0) + muW0;
			w[i] = w0[i];
		}
		// ELSE all eqaul to abs(mu)
		else{ 
			w0[i] = abs(muW0);
			w[i] = w0[i];
		}

		// LEARNING SKILL/TALENT //
		// IF mu positive -->> Gussian distributed 
		if(muXI > 0 ){  
  			xi[i]=gsl_ran_gaussian(gslpointer,sigmagXI) + muXI;
		}
		// ELSE all eqaul to abs(mu)
		else{ 
			xi[i] = abs(muXI);
		}

		// All initial trategies are full defection = contribution of all agents is 0
		alpha[i]=0;
	
		// The memory array 
		for (int j=0; j<N_alpha; j++){
			int ID = i*N_alpha + j;
			memory[ID] = 1.0;
		}
	}
	
	// Print out the INVESTMENT TALENT for each player 
	ensemblePrintTalent();

	// Print out the INVESTMENT CAP for each player 
	ensemblePrintW0();

	// Print out the LEARNING SKILL/TALENT for each player 
	ensemblePrintXI();
	
	// Total initial wealth. I don't really need to initialize this, but it's better anyway
	wealth = sumofvector(w, N);

	//It is important to initialize this!
	oldwealth = wealth;
	
	// Total investment cap
	totInvestCap = sumofvector(w0, N);
	
	
	// Print the initial state // 
	// First column: Time Step = 0, Following Columns: 
	// filet/time.txt --> Ensemble averages for total wealth, growth%/eff, Gini coef., avg Co-op/strategy
	// filew/wealth.txt --> Total wealth for each 'society'
	// filec/cooperation.txt --> Average co-op/strategy for each 'society'
	// fileg/GiniCoef.txt --> Gini coef. for each 'society' 
	//ensemblePrint(filet, filew, filec, fileg, t, eff, wealth, w, alpha, cons, N, NE);
	ensemblePrint();
	
	
	// Print for testing
	cout << "Initialisation of SML society nr " << k << " complete" << endl;
	cout << endl;	

}


// --- The society initialization method for Nash Equalibrium schmea --- //
void GameTheorySim::initSociety_NashEQ_NWA(int k){

	// Print for testing
	cout << "Initialising NashEQ society nr " << k << endl;
	
	// FILL THE VECTORS OF TALENT , WEALTH AND RANKING AND INITIALIZE EVERYTHING

	// Initial timestep is 0
	t=0;			

	// Efficiency at t=0 is 0
	eff = 0.0;

			
	// The arrays describing agents
	for(int i=0; i<N; i++){
		

		// Intial rank is efectivly random
		rank[i]=i;		
				
		// INVESTMENT TALENT //
		// IF mu positive -->> Gussian distributed 
		if(muR > 0 ){  
  			r[i]=gsl_ran_gaussian(gslpointer,sigmagR) + muR;
		}
		// ELSE all eqaul to abs(mu)
		else{ 
			r[i] = abs(muR);
		}

		// INVESTMENT CAP & INITIAL WEALTH //
		// IF mu positive -->> Gussian distributed 
		if(muW0 > 0 ){  
  			w0[i]=gsl_ran_gaussian(gslpointer,sigmagW0) + muW0;
			w[i] = w0[i];
		}
		// ELSE all eqaul to abs(mu)
		else{ 
			w0[i] = abs(muW0);
			w[i] = w0[i];
		}

		// All initial trategies are full defection = contribution of all agents is 0
		alpha[i]=0;
	
	}
	
	// Print out the INVESTMENT TALENT for each player 
	ensemblePrintTalent();

	// Print out the INVESTMENT CAP for each player 
	ensemblePrintW0();

	// Total initial wealth. I don't really need to initialize this, but it's better anyway
	wealth = sumofvector(w, N);

	//It is important to initialize this!
	oldwealth = wealth;

	// Total investment cap
	totInvestCap = sumofvector(w0, N);
	
	
	// Print the initial state // 
	// First column: Time Step = 0, Following Columns: 
	// filet/time.txt --> Ensemble averages for total wealth, growth%/eff, Gini coef., avg Co-op/strategy
	// filew/wealth.txt --> Total wealth for each 'society'
	// filec/cooperation.txt --> Average co-op/strategy for each 'society'
	// fileg/GiniCoef.txt --> Gini coef. for each 'society' 
	//ensemblePrint(filet, filew, filec, fileg, t, eff, wealth, w, alpha, cons, N, NE);
	ensemblePrint();
	
	
	// Print for testing
	cout << "Initialisation of NashEQ society nr " << k << " complete" << endl;
	cout << endl;	


}



// --------------------------------- //
// ----- Running the Societies ----- //
// --------------------------------- //


// --- Run simulation for one full society for SML -- //
void GameTheorySim::runSociety_SML_NWA(int k){


	// ************************************************************* //
	// ********************* THE BIG TIME LOOP ********************* //		
	// ************************************************************* //

	for ( t=1; t<(T+1); t++ ){

		// ** Save the previous wealth vector --> for memory updating ** //
		for (int j=0; j<N; j++){
			OLDw[j] = w[j];
		}

		// Testing 
		// cout << "Test1" << endl;

		// *** Update strategy for all players *** //
		for( i=0; i<N; i++){

			// The strategy is updated based on the memory of the agent
			//updatestrategy(i,alpha,oldalpha,w,r,cons, gslpointer);
			alpha[i] = newStrategySML_NWA(i);

		}

		// Testing 
		// cout << "Test2: Number of agent in a society is N = " << N << endl;

		// *** "Play the game" with the new strategys *** //

		// I generate the ranking and compute the output for each player.
		makeorder_NWA(alpha, rank, O); 
		
		// Testing 
		// cout << "Test3: Number of groups M = " << M << endl;

		// Here I use the ranking to generate the groups and update the wealth of the players.
		//formgroups(w, alpha, O,rank, cons, M); 
		formgroups_NWA(w, alpha, O, rank); 

		// Testing 
		// cout << "Test4" << endl;

		// Here I compute the total amount of wealth. And the efficiency, which is the production per timestep.
		wealth = sumofvector(w, N);
		eff = wealth - oldwealth;
		//eff = eff/oldwealth; // I separate those 2 because in this way it should be better for numerical errors, right?
		
		// --- Because there is a maximum limit on amount of investment, it would be correct to calculate efficiency compared to the max limit --- //  	
		eff = eff/(totInvestCap);
		oldwealth = wealth;		 

		// Testing 
		// cout << "Test5" << endl;

		// Each player updates their memory for each possible strategy
		// they could have taken, by fixing the strategies of other 
		// player to the ones they picked and seeing how much they would
		// have won. 
		this->memoryUpdate_NWA();

		// Testing 
		// cout << "Test6" << endl;


		// *** Outputting *** // 


		// Print the state at time step t // 
		// First column: Time Step = t, Following Columns: 
		// filet/time.txt --> Ensemble averages for total wealth, growth%/eff, Gini coef., avg Co-op/strategy
		// filew/wealth.txt --> Total wealth for each 'society'
		// filec/cooperation.txt --> Average co-op/strategy for each 'society'
		// fileg/GiniCoef.txt --> Gini coef. for each 'society' 
		//ensemblePrint(filet, filew, filec, fileg, t, eff, wealth, w, alpha, cons, N, NE);
		ensemblePrint();
			


		// Print the time to commandline for checking the runtime //
		cout << "SML Society nr " << k + 1 << ", timestep "<< t << endl; 		


	}
	
	// ********************* END OF BIG TIME LOOP ********************* //

	// --- End the line for k-th society in output files --- //
	filew << endl;
	filec << endl;
	fileg << endl;
	fileEf << endl;

	
	// Print for testing
	cout << endl;		
	cout << "Simulation for SML society nr " << k << " complete" << endl;
	cout << endl;
	cout << endl;



}





// --- Run simulation for one full society for NashEQ -- //
void GameTheorySim::runSociety_NashEQ_NWA(int k){

	// ************************************************************* //
	// ********************* THE BIG TIME LOOP ********************* //		
	// ************************************************************* //

	for ( t=1; t<(T+1); t++ ){

		// *** Save the old strategy *** // 		
		for(int i=0; i<(N); i++){
			oldalpha[i] = alpha[i];
		}


		
		// ** Update strategy for each player in each in k-th society ** //
		for (int i=0; i<N; i++){
			//alpha[k*N + i] = NewStrategy(i, alpha_k, w_k, r_k, cons, gslpointer);
			alpha[i] = NewStrategyNashEQ_NWA(i);
		}
		
		// *** "Play the game" with the new strategys *** //

		// I generate the ranking and compute the output for each player.
		makeorder_NWA(alpha, rank, O); 
		
		// Testing 
		// cout << "Test3" << endl;

		// Here I use the ranking to generate the groups and update the wealth of the players.
		//formgroups(w, alpha, O,rank, cons, M); 
		formgroups_NWA(w, alpha, O,rank); 

		// Testing 
		// cout << "Test4" << endl;

		// Here I compute the total amount of wealth. And the efficiency, which is the production per timestep.
		wealth = sumofvector(w, N);
		eff = wealth - oldwealth;
		//eff = eff/oldwealth; // I separate those 2 because in this way it should be better for numerical errors, right?
		
		// --- Because there is a maximum limit on amount of investment, it would be correct to calculate efficiency compared to the max limit --- //  	
		eff = eff/(totInvestCap);
		oldwealth = wealth;		 

		// Testing 
		// cout << "Test5" << endl;	 		


		// *** Outputting *** // 


		// Print the state at time step t // 
		// First column: Time Step = t, Following Columns: 
		// filet/time.txt --> Ensemble averages for total wealth, growth%/eff, Gini coef., avg Co-op/strategy
		// filew/wealth.txt --> Total wealth for each 'society'
		// filec/cooperation.txt --> Average co-op/strategy for each 'society'
		// fileg/GiniCoef.txt --> Gini coef. for each 'society' 
		//ensemblePrint(filet, filew, filec, fileg, t, eff, wealth, w, alpha, cons, N, NE);
		ensemblePrint();
			


		// Print the time to commandline for checking the runtime //
		cout << "NashEQ Society nr " << k + 1 << ", timestep "<< t << endl; 		


	}
	
	// ********************* END OF BIG TIME LOOP ********************* //

	// --- End the line for k-th society in output files --- //
	filew << endl;
	filec << endl;
	fileg << endl;
	fileEf << endl;

	// Print for testing
	cout << endl;		
	cout << "Simulation for NashEQ society nr " << k << " complete" << endl;
	cout << endl;
	cout << endl;


}








// ----------------------------------------------- //
// ----- The main methods for the simulation ----- //
// ----------------------------------------------- //

// --- makeorder_NWA -->> No wealth accumulation --- //

void GameTheorySim::makeorder_NWA(double *alpha,int *rank, double *O){

	int enne=N; 			// Dummy variable that is equal to N
	double X[enne]; 		// This is the quantity that will decide the ordering
        vector<int> Y (enne,0);		// I have to create the array due to how sort works
        vector<int> Shuf (enne,0);	// Vector that I will shuffle, to provide a random ranking in case of ties
	int i; 				// Dummy variable

	// ****************** Generate the output vector //
	for(i=0;i<enne;i++){ 
		O[i]=r[i]*alpha[i]*w0[i];
		//cout << "Testing makeorder_NWA. Output of agent " << i << " is " << O[i] << endl;
                Y[i]=i;
	}
	// ********************* COMPUTE THE ORDERING ARRAY **************** //
	switch(Choice){
		case 1:
			for(i=0;i<enne;i++){
				X[i]=O[i];
			}	
			//cout << "Testing makeorder_NWA. Choice is 1" << endl;		
			break;
		case 2:
			for(i=0;i<enne;i++){
				X[i]=alpha[i]*w0[i];
			}
			//cout << "Testing makeorder_NWA. Choice is 2" << endl;
			break;
		case 3:
			for(i=0;i<enne;i++){
				X[i]=r[i]*alpha[i];
			}			
			//cout << "Testing makeorder_NWA. Choice is 3" << endl;
			break;
		case 4:
			for(i=0;i<enne;i++){
				X[i]=alpha[i];
			}
			//cout << "Testing makeorder_NWA. Choice is 4" << endl;
			break;
		default:
			cout<<"ERROR IN THE SWITCH IN THE MAKEORDER FUNCTION";
			break;
	}

	// *************** DEFINITION OF THE COMPARING FCT AND THEN SORTING *************** //
        
	// Here I create a random order that I can use to break ties. Here i fill in the vector with increasing numbers
        for(i=0;i<enne;i++){ 
            Shuf[i]=i;
        }
        random_shuffle(begin(Shuf), end(Shuf)); 	// Here I shuffle the vector
        
	auto fancycompare = [&](int i, int j){ 		// This function order the agents in terms of their output. 
							// With this particular configuration it puts the agents 
							// that contribute more in the beginning and the others at the end.
		double trial;
	
		trial= X[i]-X[j];
		if(trial==0){
			trial=Shuf[i] - Shuf[j]; 	// If the 2 contributions are equal, I resort to the ordering of Shuf
		}
		//Hence trial can never be ==0
                if(trial<0){return false;} 		// If you want ordering s.t. the lowest contributors are in the beginning, you have to switch that.
		else{return true;}	
	};	
	
	sort(Y.begin(),Y.end(),fancycompare); 		// Here I sort the array
	
	// ******************** COPYING THE RESULTS BACK IN THE ORIGINAL ARRAY **************** //

	// Here I'm just copying the result back in the original array
	for(i=0;i<enne;i++){ 
		rank[i]=Y[i];
	}
	return;
}




// --- formgroups_NWA -->> No wealth accumulation --- //

void GameTheorySim::formgroups_NWA(double *w, double *alpha, double *O, int *rank){

	// i and j are dummy variables, k needs to be updated by hand and will be used for the rank
	int i,j,k; 
	double gain; // The gain in the group

	//cout << "Testing formgroups_NWA: M = " << M << endl;
	k=0; // Initialize k to zero
	for(i=0;i<M;i++){ // Here i sum over all the groups
		
		gain=0; // Set the gain of the group to 0.
		
		for(j=0;j<S;j++){
			gain = gain + O[rank[k+j]]; //Here I am adding to the common pool the contribution due to rank[k]+j
		}
		
		gain = gain*Q; // Correct for the return rate
		
		for(j=0; j<S;j++){
			// !!!!! THIS IS DIFFERENT THEN THE USUAL CASE !!!!!! //
			w[rank[k]] = w[rank[k]] - w0[rank[k]]*alpha[rank[k]] + gain; //Here I update the wealth of all the players
			//cout << "Testing formgroups. Wealth of agent " << rank[k] << " is " << w[rank[k]] <<endl; 
			k++; //Here I update k so that it keeps track for the next round
		}
		
	}
	return;
}



// --- getutility_NWA -->> No wealth accumulation --- //

double GameTheorySim::getutility_NWA(int q, double *w, double *alpha, double *O, int *rank){
	int j,k; //i and j are dummy variables, k needs to be updated by hand and will be used for the rank
	double gain, utility; //The gain in the group and the utlity of q
        int g; //The group in which q is
        int place; //The rank of q
        
        // As first thing we have to compute what is the ranking of //
        place=inverserank(q,rank);
        // As second thing we have to compute in which group the q-th agent is ending up //
        g = place/S;
        
        // Then I skip directly to that group and compute the gain in that group //
        k=S * g; //I start from the k-th ranked agent
        gain=0; //Set the gain of the group to 0.
	
        for(j=0;j<S;j++){
            gain = gain + O[rank[k+j]]; //Here I am adding to the common pool the contribution due to rank[k+j]. Due to how the systyem is designed, one of those people is q
	}
		
	gain = gain*Q; //Correct for the return rate
		
        utility = w[q] -w0[q]*alpha[q] + gain; //This is the utility
	//cout << "Testing getutility! Utility of agent " << q << " is " << utility << endl;	
	
	return utility;
}





// --- newStrategyNashEQ_NWA -->> Nash Equlibrium updating && No Wealth Accumulation--- // 

double GameTheorySim::NewStrategyNashEQ_NWA(int k){
    
    double utility, oldutility; // they contain the possible utility for player k for a strategy and the previous one
    double dummyalpha[N]; //A dummy vector where I copy all the strategies of the other players and then I update the one of player k
    int oldplace=0; //Dummy vector to store the old ranking (to check if the ranking changes or not according to the new strategy). Initialize to a negative for the first time it checks the if
    int place; //To store the ranking of the k-th agent
    int i;
    int dummyrank[N]; //The vector that will do all the ranking
    double dummyO[N]; //The dummy vector of outputs
    
	vector<double> cumprob (N_alpha,0); //The array with the cumulative probabilities

	// ***** FOR probarr[] and sum I'm going to use myMP_float from myEXP.h ****** //
	myMP_float probarr[N_alpha];
	myMP_float sum;
    


    for(i=0; i<N ; i++){ //Here I copy the vector of strategy into this fake one. Note that I don't directly pass the pointer because I want to manipulate this one. Hence I have to copy all the array
        dummyalpha[i]=oldalpha[i];
	//dummyalpha[i]=alpha[i];
    }

    sum=0; //This is to normalize the probability at the end
    
    // This is the loop over all the possible strategies for player k //
    // Here I also compute the probability of a certain strategy //
    for(i=0; i < N_alpha ; i++){
        dummyalpha[k]  = i*alpha_interval; //Set the strategy of k 

        //makeorder(r,w,dummyalpha,dummyrank,cons,dummyO); //Make the rankings
	// NWA ordering        
	makeorder_NWA(dummyalpha, dummyrank, dummyO);

	place=inverserank(k,dummyrank); //The position where k is placed
        if(place==oldplace){ //If the rank of k is the same, I don't have to reform the groups. I just update the utility
            utility=oldutility + alpha_interval * w0[k] * (Q * r[k] -1); //It's fine because for i=0, I'm sure that this if will never take place because oldplace is initialized to a negative number
        }
        else{ //In this other case I actually have to reform the group
            utility=getutility_NWA(k, w, dummyalpha, dummyO, dummyrank); //Actually form the group and compute the utility
        }
        
	
        //probarr[i] = logitprob(beta , utility); //Here I compute the probability according to the logit and the sum
	
	// ***** I use my_mpEXP from myEXP.h for taking the exponent ***** //
	probarr[i] = my_mpEXP(beta*utility);
        sum = sum + probarr[i];

	// ***** Now print the exponent for control ***** //
	//cout << "mpEXP(" << beta << "*" << utility << "): " << probarr[i] << endl;
        
        oldutility=utility;
        oldplace=place;
    }
    
    // Now let's renormalize the probabilities and create the array of cumulative probabilities //
    
    probarr[0] =  probarr[0]/sum;

	// ***** Here I declear the myMP_float type probarr[i] to be double (just in case there is no automatic conversion ***** //
	cumprob[0] =  double(probarr[0]); //I have to do the first by end.
    for(i=1;i < N_alpha; i++){
        probarr[i]=  probarr[i]/sum;

	// ***** Here I declear the myMP_float type probarr[i] to be double (just in case there is no automatic conversion ***** //
        //cumprob[i]=cumprob[i-1] + probarr[i]; //Here I sum the cumulative probability
	cumprob[i]=cumprob[i-1] + double(probarr[i]); //Here I sum the cumulative probability
    }
    // Now let's sample the strategy //
    
	// ***** AS sum now is of type myMP_float, I'm going to use an other variable for the RNG ****** //

//    sum= gsl_ran_flat(gslpointer,0,1); //Generate a random number btw 0 and 1. Note that I use sum here just to not use another variable
//    i=binaryprobsearch(&cumprob[0],p,sum); //Here I compute which strategy is the agent k using. Note that I use i here just to not use another variable


    double rng_tmp = gsl_ran_flat(gslpointer,0,1); //Generate a random number btw 0 and 1. Note that I use sum here just to not use another variable
    i=binaryprobsearch(&cumprob[0],N_alpha,rng_tmp); //Here I compute which strategy is the agent k using. Note that I use i here just to not use another variable
   // Note that if I had to sample many times from the same distribution, it would have been smarter to use the GSL tool for "General Discrete Distributions"
   // Here I need to sample only one, hence I sample it by hand, using a binary search method 
    
    //alpha[k]=i*alpha_interval; //Here I set the value for alpha[k]. I.e. I set the strategy for agent k
    
	// Return the new strategy
    return i*alpha_interval;
}





// --- newStrategySML_NWA -->> Simple learning algorithm using memory && No Wealth Accumulation--- //

double GameTheorySim::newStrategySML_NWA(int k){

	vector<double> cumprob (N_alpha,0); //The array with the cumulative probabilities

	// ***** FOR probarr[] and sum I'm going to use myMP_float from myEXP.h ****** //
	myMP_float probarr[N_alpha];
	myMP_float sum;
    
	sum=0; //This is to normalize the probability at the end
    
	// ************* This is the loop over all the possible strategies for player k *****************
	//  Here I also compute the probability of a certain strategy 
	for(int i=0; i < N_alpha ; i++){
	
		// Find the appropriate ID
		int ID = k*N_alpha + i;
		
		// Testing
		//cout << memory[ID] << endl;

		//Here I compute the probability according to the logit and the sum
	    	// ***** I use my_mpEXP from myEXP.h for taking the exponent ***** //
	    	probarr[i] = my_mpEXP(beta*memory[ID]);
	    	sum = sum + probarr[i];
	
	    	// ***** Now print the exponent for control ***** //
	    	//cout << "mpEXP(" << beta << "*" << utility << "): " << probarr[i] << endl;
	}
	
	// Now let's renormalize the probabilities and create the array of cumulative probabilities //
	
	probarr[0] =  probarr[0]/sum;
	
	// ***** Here I declare the myMP_float type probarr[i] to be double (just in case there is no automatic conversion ***** //
	cumprob[0] =  double(probarr[0]); //I have to do the first by end.
	for(int i=1; i < N_alpha; i++){
	    probarr[i]=  probarr[i]/sum;
	    cumprob[i]=cumprob[i-1] + double(probarr[i]); //Here I sum the cumulative probability
	}
	// Now let's sample the strategy //
	
	// ***** AS sum now is of type myMP_float, I'm going to use an other variable for the RNG ****** //

	// Generate a random number btw 0 and 1.
	double RN = gsl_ran_flat(gslpointer,0,1);  
	// Here I compute which strategy is the agent k using. Note that I use i here just to not use another variable
	int tmp_out = binaryprobsearch(&cumprob[0],N_alpha,RN); 
	
	double out = tmp_out*alpha_interval; //Here I set the value for alpha[k]. I.e. I set the strategy for agent k
    
	return out;
}





// --- memoryUpdate_NWA -->> SML memory update && No Wealth Accumulation--- //

// Each player updates their memory for each possible strategy
// they could have taken, by fixing the strategies of other 
// player to the ones they picked and seeing how much they would
// have won.

void GameTheorySim::memoryUpdate_NWA(){

	// --- Dummy variables --- //
	int dummyrank[N];		//The vector that will do all the ranking
	double dummyO[N];		//The dummy vector of outputs
	double dummyalpha[N];	//A dummy vector where I copy all the strategies of the other players and then I update the one of player k
	int oldplace=-1;		//Dummy vector to store the old ranking (to check if the ranking changes or not according to the new strategy). Initialize to a negative
	int place;			//To store the ranking of the k-th agent
	double utility, oldutility;	//they contain the possible utility for player k for a strategy and the previous one


	// --- Fill the dummyalpha --- //
	for (int k=0; k<N; k++){
		dummyalpha[k] = alpha[k];
	}

	// ----- Loop over players ------ //
	for (int k=0; k<N; k++){

		// ---- Loop over all strategies for k-th player ---- //
		for (int i=0; i<N_alpha; i++){


			// -- Play the game for the k-th strategy -- //

			// Set the strategy of k
			dummyalpha[k]  = i*alpha_interval;

			// Make the rankings 
			//makeorder(r,OLDw,dummyalpha,dummyrank,cons,dummyO);
			// !!! This is where NO_WEALTH_ACCUMULATION is hidden !!! //
			//makeorder_NWA(r, W0, alpha, rank, cons, O); 
			makeorder_NWA(dummyalpha, dummyrank, dummyO);

			// The position where k is placed	
			place=inverserank(k,dummyrank);			
		
			// Calculate the utility for the k-th player //
			if(!(place==oldplace)){ 

				// Form the group and compute the utility
            			utility=getutility_NWA(k, OLDw, dummyalpha, dummyO, dummyrank); 
        		}
			// If the rank of k is the same, I don't have to reform the groups. I just update the utility
		        else{ 
		             utility=oldutility + alpha_interval * w0[k] * (Q * r[k] -1);
        		}



			// -- Now that I have the utility I can update the memory -- //
			
			// Calculate the growth for the k-th agent
			//double tmp_growth = utility/OLDw[k];
			double tmp_growth = (utility - OLDw[k])/(w0[k]);

			// Get the memory ID
			int ID = k*N_alpha + i;

			// Find the current memory value
			float tmp_mem = memory[ID];

			// Calculate the new memory
			memory[ID] = tmp_mem*(1.0 + xi[k]*tmp_growth);





			// -- Save the oldutility and oldplace parameters -- //	
			oldutility=utility;
			oldplace=place;
		}


		// ---- Reset the correct strategy for k-th player ---- //
		dummyalpha[k] = alpha[k];

	}
}






// ------------------------------- //
// ----- The support methods ----- //
// ------------------------------- //


// --- Binary search --- //
// Gamma --> the array of cumulative prob
// M --> the lenght of the array 
// x --> the random number between zero and one
int GameTheorySim::binaryprobsearch(double *Gamma, int M, double x) { 
    int a, b, l, result;
    bool check;
    a = 0;
    b = M - 1;
    
    do {
        l = (a + b) / 2;
        if (x <= Gamma[l]) {
            if ((x >= Gamma[l - 1])&&(l > 0)) {
                result = l;
                check = true;
            } else {
                if (l > 0) {
                    b = l;
                    check = false;
                } else {
                    result = 0;
                    check = true;
                }
            }
        } else {
            if (x <= Gamma[l + 1]) {
                result = l + 1;
                check = true;
            } else {
                a = l;
                check = false;
            }
        }
    } while (check == false);
    return result;
}



// --- Inversrank -->> returns the rank of the agent k. It's the inverse of rank //
int GameTheorySim::inverserank(int k, int *rank){
    int i; //dummy variable
    int out; //the output
    
    for(i=0;i< N; i++){
        if(rank[i]==k){
            out=i;
            break;
        }
    }
    return out;
}

// --- computegini -->> Compute the Gini coef. //
// The Gini coefficient is computed as defined from Wikipedia, in particular,
// given that y[i] is the wealth of an individual i, and y[i] is ordered in increasing order (i.e. y[i]<=y[i+1]) we have:
// G=1/n(n+1-2((Sum{(n+1-i)*y[i]})/(Sum{y[i]})))

double GameTheorySim::computegini(double *w){
	
	int enne=N; 
	vector<double> y(w,w+enne); //Vector of wealth. I have to define it like that because I have to order it
	double G;
	double dummy1, dummy2; //dummy doubles
	int i; //dummy integer
	
	sort(y.begin(),y.end()); //Here I sort the wealth vector increasing order
	
	dummy1=0;
	for(i=0; i<enne;i++){ //Here I compute the sum at the numerator
		dummy1=dummy1+(enne-i)*y[i]; //Here it should be enne +1 -i for i that goes between 1 and enne => enne +1 -(i +1) for i that goes from 0 to enne-1
	}
	
	dummy2=0;
	for(i=0; i<enne;i++){ //Here I compute the sum at the denominator
		dummy2=dummy2+y[i];
	}
	
	G=enne + 1. -2.*(dummy1/dummy2); //Here I finish the computation
	G=G/enne;
	
	return G;
}


// --- sumofvector -->> This function just returns the sum over the vector of length n //
double GameTheorySim::sumofvector(double *a, int n){
    double res;
    int i;
    
    res=0;
    
    for(i=0;i<n;i++){
        res = res + a[i];
    }
    
    return res;
}



// mpEXPcpp - Exponent function using Boost::multiprecision cpp_dec_float 

myMP_float GameTheorySim::my_mpEXP(double x){
	myMP_float tmpX = myMP_float(x);
	
	return boost::multiprecision::exp(tmpX); 
}






// ------------------------------------------- //
// ----- The printing/outputting methods ----- //
// ------------------------------------------- //

// --- Ensemble print -->> The standart output --- //
void GameTheorySim::ensemblePrint(){

	double Gini = 0.0;		// Gini for each society
	double Coop = 0.0;		// Average co-op for each society     

	// Calculate the Gini coef & Co-op and their and total wealth & growth % ensemble averages //

	// The arithmetic average co-op for k-th 'society'
	for (int i=0; i<N; i++){
		Coop += alpha[i];
	}
	Coop = Coop/N;

	// Gini coef
	Gini = computegini(w);

	
	// Printing //

	filew  << left << setw(12) << wealth << "\t";
	filec  << left << setw(12) << Coop << "\t";
	fileg  << left << setw(12) << Gini << "\t";
	fileEf  << left << setw(12) << eff << "\t";

	// End the print function
	return;
}



// --- ensemblePrintParams -->> prints the parameters.txt file  --- //
void GameTheorySim::ensemblePrintParams(){
	filep.open(filenamep,ios::out|ios::trunc);
	if(filep.is_open()){
		filep << "# Number of 'societies' in the ensemble" << endl;
		filep << "NE= " << NE << endl << endl;
		filep<<"#Number of agents playing the game:"<<endl;
		filep<<"N= "<<N<<endl<<endl;
		filep<<"#Time at which I stop the simulation:"<<endl;
		filep<<"T= "<<T<<endl<<endl;
		filep<<"#Size of the groups:"<<endl;
		filep<<"S= "<<S<<endl<<endl;
		filep<<"#Return rate:"<<endl;
		filep<<"Q= "<<Q<<endl<<endl;
		filep<<"#Indicates in how many intervals I want to divide the interval [0,1] for the contributions. The number of strategies available to each player is this number +1"<<endl;
		filep<<"monincr= "<<monincr<<endl<<endl;
		filep<<"#This is the noise parameter in the logit:"<<endl;
		filep<<"beta= "<<beta<<endl<<endl;
		filep<<"#Number that indicates the choice of what are we using to rank agents in the groups:"<<endl;
		filep<<"Choice= "<<Choice<<endl<<endl;
		filep<<"#Mean of the Gaussian distribution for the INVESTMENT TALENT --> Negative == const. as abs(muR):"<<endl;
		filep<<"muR= "<<muR<<endl<<endl;
		filep<<"#Variance sigma of the Gaussian distribution for the IVESTMENT TALENT:"<<endl;
		filep<<"sigmagR= "<<sigmagR<<endl<<endl;
		filep<<"#Mean of the Gaussian distribution for the INVESTMENT CAP --> Negative == const. as abs(muW0):"<<endl;
		filep<<"muW0= "<<muW0<<endl<<endl;
		filep<<"#Variance sigma of the Gaussian distribution for the IVESTMENT CAP:"<<endl;
		filep<<"sigmagW0= "<<sigmagW0<<endl<<endl;
		filep<<"#Mean of the Gaussian distribution for the LEARNING SKILL/TALENT --> Negative == const. as abs(muXI):"<<endl;
		filep<<"muXI= "<<muXI<<endl<<endl;
		filep<<"#Variance sigma of the Gaussian distribution for the LEARNING SKILL/TALENT:"<<endl;
		filep<<"sigmagXI= "<<sigmagXI<<endl<<endl;
	}
	filep.close();
	return;
}


// --- ensemblePrintTalent -->> Prints the INVESTMENT TALENT to the talent.txt file --- //
void GameTheorySim::ensemblePrintTalent(){

	if(filer.is_open()){
		// Print out the talent for each player 
		for(int i=0; i<N; i++){

			filer << left << setw(8) << r[i] << "\t" ;
				
		}	
	}
}



// --- ensemblePrintW0 -->> Prints the INVESTMENT CAP to the talent.txt file --- //
void GameTheorySim::ensemblePrintW0(){

	if(fileW0.is_open()){
		// Print out the talent for each player 
		for(int i=0; i<N; i++){

			fileW0 << left << setw(8) << w0[i] << "\t" ;
				
		}	
	}
}


// --- ensemblePrintXI -->> Prints the LEARNING SKILL/TALENT to the talent.txt file --- //
void GameTheorySim::ensemblePrintXI(){

	if(fileXI.is_open()){
		// Print out the talent for each player 
		for(int i=0; i<N; i++){

			fileXI << left << setw(8) << xi[i] << "\t" ;
				
		}	
	}
}



