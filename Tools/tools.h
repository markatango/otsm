#pragma once

#define WINDOWS_PLATFORM

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <cmath>
#include <utility>
#include <limits>
#include <queue>
#include <ctime>
#include <utility>

#ifdef WINDOWS_PLATFORM
#include <windows.h>
#endif

#define TOOLS_BEGIN namespace tools{
#define TOOLS_END };

using std::vector;
using std::string;

TOOLS_BEGIN

vector<string> Split(const string& strLine, char delim);

//	Thread Timer
class CThreadTime
{
public:
	void Begin();
	__int64 End();

private:
	__int64 FileTimeToQuadWord(PFILETIME pft);

private:  
	FILETIME ftKernelTimeStart;  
	FILETIME ftKernelTimeEnd;  
	FILETIME ftUserTimeStart;  
	FILETIME ftUserTimeEnd;  
	FILETIME ftDummy;   

};

TOOLS_END