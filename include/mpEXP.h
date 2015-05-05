// include BOOST multiprecision cpp_dec_float
#include <boost/multiprecision/cpp_dec_float.hpp>



// For convinience use boost::multiprecision as namespace
using namespace boost::multiprecision;


// mpEXPcpp - Exponent function using Boost::multiprecision cpp_dec_float 

cpp_dec_float_100 mpEXPcpp(double x){
	cpp_dec_float_100 tmpX = cpp_dec_float_100(x);
	
	return boost::multiprecision::exp(tmpX); 
}


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
