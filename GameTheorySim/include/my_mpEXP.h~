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

myMP_float my_mpEXP(double x){
	myMP_float tmpX = myMP_float(x);
	
	return boost::multiprecision::exp(tmpX); 
}



