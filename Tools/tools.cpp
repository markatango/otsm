#include "tools.h"

TOOLS_BEGIN


vector<string> Split(const string& strLine, char delim)
{
	vector<string> vecTokens;
	string cur;
	for(size_t i = 0; i < strLine.size(); ++i)
	{
		if(strLine[i] == delim)
		{
			if(0 != cur.size())
			{
				vecTokens.push_back(cur);
				cur = "";
			}
		}
		else if(isspace(strLine[i]))	//	skip space
			continue;
		else 
			cur.push_back(strLine[i]);
	}
	if(0 != cur.size())
		vecTokens.push_back(cur);

	return vecTokens;
}

void CThreadTime::Begin()
{
	GetThreadTimes(GetCurrentThread(), &ftDummy, &ftDummy, &ftKernelTimeStart, &ftUserTimeStart);
}

__int64 CThreadTime::End()
{
	GetThreadTimes(GetCurrentThread(), &ftDummy, &ftDummy, &ftKernelTimeEnd, &ftUserTimeEnd);  

	__int64 qwKernelTimeElapsed = FileTimeToQuadWord(&ftKernelTimeEnd)  - FileTimeToQuadWord(&ftKernelTimeStart);  
	__int64 qwUserTimeElapsed = FileTimeToQuadWord(&ftUserTimeEnd) - FileTimeToQuadWord(&ftUserTimeStart);  

	//   Get total time duration by adding the kernel and user times.  
	//   the default is 100ns, so we convert it to ms  
	return (qwKernelTimeElapsed + qwUserTimeElapsed) / 10000;   
}

__int64 CThreadTime::FileTimeToQuadWord(PFILETIME pft)
{
	return (Int64ShllMod32(pft->dwHighDateTime, 32) | pft->dwLowDateTime);
}

TOOLS_END