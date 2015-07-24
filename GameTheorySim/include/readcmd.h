#include <stdlib.h>
#include <algorithm>
#include <fstream>


using namespace std;

// Tools for reading commandline option
bool cmdOptionExists(char** begin, char** end, const std::string& option)
{
	return std::find(begin, end, option) != end;
}

char* getCmdOption(char ** begin, char ** end, const std::string & option)
{
	char ** itr = std::find(begin, end, option);
	if (itr != end && ++itr != end)
	{
	    return *itr;
	}
	return 0;
}

