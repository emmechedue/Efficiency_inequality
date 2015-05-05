#####################################################
# You may need to change the parameters under here. #

CXX=g++

# Flags for the libraries etc.
LDFLAGS = -lgsl -lgslcblas

# Flag for your Boost 
BOOST=/usr/local/boost_1_58_0

# Flags for include folder (header files)
INCFLAGS = -Iinclude -I$(BOOST)

# Set default compiler parameters
# -Wall 	shows all warnings when compiling, always use this!
# -std=c++11 	enables the C++11 standard mode
CXXFLAGS = -Wall -std=c++11 $(LDFLAGS) $(INCFLAGS) 



#####################################################

# Declare the name of our program 
PROGRAM = inequality_inefficiency_logit_single_loop

# The needed object files (we make these from the source code)
OBJ = obj/main.o obj/Constants.o 

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

obj/Constants.o: src/Constants.cpp
	$(CXX) -c src/Constants.cpp -o obj/Constants.o $(CXXFLAGS) 
