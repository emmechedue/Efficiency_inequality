// A small class for getting printing the writing
// time into the output file.


#ifndef TIMEENGINE_H_INCLUDED
#define TIMEENGINE_H_INCLUDED

#include <stdio.h>
#include <time.h>
#include <cstdlib>
#include <string>
#include <vector>

using namespace std;

class TimeEngine{

public:

    // Variables for storing information
    // about the time
    time_t rawtime;
    struct tm * timeinfo;

    // Constructor.
    TimeEngine();

    // Destructor.
    ~TimeEngine();

    // Mathod to get the time.
    string get_timestamp();


};


#endif // TIMEENGINE_H_INCLUDED
