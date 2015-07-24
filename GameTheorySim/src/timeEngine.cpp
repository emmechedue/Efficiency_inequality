#include "timeEngine.h"


// Constructor.
TimeEngine::TimeEngine(){}

// Destructor.
TimeEngine::~TimeEngine(){}

// Mathod to get the time.
string TimeEngine::get_timestamp(){
    char tmp [19];

    time (&rawtime);
    timeinfo = localtime (&rawtime);

    strftime (tmp,19,"%d.%m.%y  %X",timeinfo);

    string rtn(tmp);

    return rtn;


}
