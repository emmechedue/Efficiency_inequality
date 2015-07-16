#####################################################
# You may need to change the parameters under here. #

CXX=g++

# Flags for the libraries etc.
# gsl for the RNG
# gmp for multiprecision --> needed for GMP & MPFR 
# mpfr for multiprecision --> needed for MPFR

# FOR LINKING THE LIBRARIES ON EULER 
GMP=$(HOME)/bin/GMP/include
MPFR=$(HOME)/bin/MPFR/include
LD_GMP=$(HOME)/bin/GMP/lib 
LD_MPFR=$(HOME)/bin/MPFR/lib 
#LDFLAGS = -lgsl -lgslcblas -L$(LD_GMP) -lgmp -lgmpxx -L$(LD_MPFR) -lmpfr
LDFLAGS = -lgsl -lgslcblas -L$(LD_GMP) -lgmp  -L$(LD_MPFR) -lmpfr


# Flag for your Boost --> needed for multiprecision (all three versions)
#BOOST=/usr/local/boost_1_58_0
#BOOST=/usr/include/boost
# BOOST ON EULER
BOOST=/cluster/apps/boost/1.55.0/x86_64/serial/gcc_4.8.2/include

# Flags for include folder (header files)
INCFLAGS = -Iinclude -I$(BOOST) -I$(MPFR) -I$(GMP)

# Set default compiler parameters
# -Wall 	shows all warnings when compiling, always use this!
# -std=c++11 	enables the C++11 standard mode
CXXFLAGS = -Wall -std=c++11 $(LDFLAGS) $(INCFLAGS) 



#####################################################

# The simu type
TYPE = ensemble_sml

# Declare the name of our program 
#PROGRAM = inequality_inefficiency_logit_single_loop
PROGRAM = run_$(TYPE)

# The needed object files (we make these from the source code)
OBJ = obj/$(TYPE)_main.o obj/$(TYPE)_Constants.o 

# This is the first target. It will be built when you run 'make' or 'make all'
all: $(PROGRAM) 

# Target for cleaning up. 
clean:	
		rm $(OBJ)
		rm $(PROGRAM)

# Rule for linking IMPORTANT! The space in front of $(CXX) is a TABULATOR!
$(PROGRAM): $(OBJ)
	$(CXX) $(OBJ) $(INCFLAGS) $(LDFLAGS) -o $(PROGRAM)

# Rules for compiling
obj/$(TYPE)_main.o: src/$(TYPE)_main.cpp
	$(CXX) -c src/$(TYPE)_main.cpp -o obj/$(TYPE)_main.o $(CXXFLAGS) 

obj/$(TYPE)_Constants.o: src/$(TYPE)_Constants.cpp
	$(CXX) -c src/$(TYPE)_Constants.cpp -o obj/$(TYPE)_Constants.o $(CXXFLAGS) 
