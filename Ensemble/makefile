#####################################################
# You may need to change the parameters under here. #

CXX=g++

# Flags for the libraries etc.
# gsl for the RNG
# gmp for multiprecision nr. 
# mpfr for multiprecision nr.
LDFLAGS = -lgsl -lgslcblas -lgmp -lgmpxx -lmpfr

# Flag for your Boost 
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
TYPE = ensemble

# Declare the name of our program 
#PROGRAM = inequality_inefficiency_logit_single_loop
PROGRAM = $(TYPE)_run

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
