#####################################################
# You may need to change the parameters under here. #

CXX=g++

# Flags for the libraries etc.
# gsl for the RNG
# gmp for multiprecision --> needed for GMP & MPFR 
# mpfr for multiprecision --> needed for MPFR
#LDFLAGS = -lgsl -lgslcblas -lgmp -lgmpxx -lmpfr
LDFLAGS = -lgsl -lgslcblas -lgmp -lmpfr

# Flag for your Boost --> needed for multiprecision (all three versions)
BOOST=/usr/local/boost_1_58_0
#BOOST=/usr/include/boost

# Flags for include folder (header files)
INCFLAGS = -Iinclude -I$(BOOST)


# Set default compiler parameters
# -Wall 	shows all warnings when compiling, always use this!
# -std=c++11 	enables the C++11 standard mode
CXXFLAGS = -Wall -std=c++11  $(LDFLAGS) $(INCFLAGS) 



#####################################################

# The simu type


# Declare the name of our program 
#PROGRAM = inequality_inefficiency_logit_single_loop
PROGRAM = run

# The needed object files (we make these from the source code)
OBJ = obj/main.o obj/GameTheorySim.o obj/GTS_quadPayoff.o obj/TimeEngine.o

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
obj/main.o: src/main.cpp
	$(CXX) -c src/main.cpp -o obj/main.o $(CXXFLAGS) 

obj/GameTheorySim.o: src/GameTheorySim.cpp
	$(CXX) -c src/GameTheorySim.cpp -o obj/GameTheorySim.o $(CXXFLAGS) 

obj/GTS_quadPayoff.o: src/GTS_quadPayoff.cpp
	$(CXX) -c src/GTS_quadPayoff.cpp -o obj/GTS_quadPayoff.o $(CXXFLAGS) 

obj/TimeEngine.o: src/timeEngine.cpp
	$(CXX) -c src/timeEngine.cpp -o obj/TimeEngine.o $(CXXFLAGS) 
