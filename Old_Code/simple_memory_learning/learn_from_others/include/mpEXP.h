// -------------------------------------- //
// --- INCLUDE THE MULTIPRECISION LIB --- //

//#include <boost/multiprecision/cpp_dec_float.hpp>
//#include <boost/multiprecision/gmp.hpp>
#include <boost/multiprecision/mpfr.hpp>



// For convinience use boost::multiprecision as namespace
using namespace boost::multiprecision;

// *** MAKE CUSTOM MULTIPRECISION TYPE *** //

// Uncomment this for the usual boost cpp multiprecision (only need boost) --> ALSO CHANGE LIB INCLUDE 
//typedef number<cpp_dec_float<256, long long> > myMP_float;

// Umcomment this for the GMP multiprecision --> ALSO CHANGE LIB INCLUDE
//typedef number<gmp_float<100> > myMP_float;

// Uncoment this for the mpfr multiprecision (The fastest and most recommended by Boost and GNU) --> ALSO CHANGE LIB INCLUDE
typedef number<mpfr_float_backend<100> > myMP_float;

// mpEXPcpp - Exponent function using Boost::multiprecision cpp_dec_float 

myMP_float mpEXP(double x){
	myMP_float tmpX = myMP_float(x);
	
	return boost::multiprecision::exp(tmpX); 
}




// ----- Stuff for testing purposes ----- //

/*
// Multiprecision (cpp_dec_float) method for calculating probability 
// S - utility array; N - lenght of S; beta - parameter beta;
// i - index of the utility for which we wish to calulate the probability
cpp_dec_float_100 mpPRcppMP(double S[], int N, int i, double beta){
	// Calculate exp(beta*S[i])
	cpp_dec_float_100 eSi = mpEXPcpp(beta*S[i]);

	// Calculate the sum of exp(beta*S[k])
	cpp_dec_float_100 sum = 0;
	for (int k=0; k<N; k++){
		sum += mpEXPcpp(beta*S[k]);
	}
	
	// Calculate the division
	cpp_dec_float_100 tmp = eSi/sum;

	// Output
	return tmp;
}

double mpPRcppDB(double S[], int N, int i, double beta){
	// Calculate exp(beta*S[i])
	cpp_dec_float_100 eSi = mpEXPcpp(beta*S[i]);

	// Calculate the sum of exp(beta*S[k])
	cpp_dec_float_100 sum = 0;
	for (int k=0; k<N; k++){
		sum += mpEXPcpp(beta*S[k]);
	}
	
	// Calculate the division
	cpp_dec_float_100 tmp = eSi/sum;

	// Output
	return double(tmp);
	
}

*/
