#include "utility.h"

UTILITY_BEGIN

void Generator(double min, double max, long long amount, string path)
{
	srand(static_cast<unsigned>(time(0)));

	if(min > max)
	{
		double tmp = min;
		min = max;
		max = tmp;
	}

	ofstream out(path.c_str(), std::ios::app);

	int precision = 100;

	int range = max;

	int jump = 1000000;

	cout << "Generate " << amount << " data." << endl;

	for(long long i = 0; i < amount; ++i)
	{
		int value = rand() % range;

		out << value << endl; 
		if(0 == i % jump)
			cout << i << " data generated." << endl;
	}
	
	out.close();
}

void SineGenerator(long long amount, string path, double a, double w, double b)
{
	ofstream out(path.c_str(), std::ios::app);

	for(long long i = 0; i < amount; ++i)
	{
		double value = a * sin(w * i) + b;
		out << value << endl;
	}

	out.close();

}

void RandomDataGenerator(long long amount, string path, double min, double max, double sharpness)
{
	if(min > max)
	{
		double tmp = min;
		min = max;
		max = tmp;
	}

	ofstream out(path.c_str());

	srand(time(static_cast<time_t*>(0)));
	
	int range = max - min;
	double top = abs(max);
	double bottom = abs(min);
	double preValue = rand() % range - bottom;

	int jump = 1000000;

	clock_t start, end;

	start = clock();

	for(long long i = 0; i < amount; ++i)
	{
		double signRate = rand() % range;

		double posRate = (max - preValue) / range;
		if(signRate >= preValue)
		{
			double add = (rand() % (int)(max - preValue)) * sharpness;
			//cout << preValue << "+" << add << endl;
			preValue += add;
		}
		else
		{
			double sub = (rand() % (int)(preValue - min)) * sharpness;
			//cout << preValue << "-" << sub << endl;
			preValue -= sub;
		}

		out << preValue << endl;

		/*
		if(0 == i % jump)
		{
			end = clock();
			cout << i << " records generated, time elapse:" << (end - start) / CLOCKS_PER_SEC << " seconds." << endl;
		}
		*/
	}

	out.close();
}

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


bool ReadDataFromFile(vector<double>& vecDatum, const string& strFilename, char delim, int iStartPos)
{
	ifstream in(strFilename.c_str());
	istringstream iss;

	int j = 0;

	while(false == in.eof())
	{
		if(++j < iStartPos)
			continue;

		string line;
		getline(in, line);

		if(0 == line.size() || (false == isspace(line[0]) && false == isdigit(line[0]) && line[0] != '-'))
			continue;
		else
		{
			vector<string> vecLine = Split(line, delim);

			for(size_t i = 0; i < vecLine.size(); ++i)
			{
				iss.str(vecLine[i]);
				double dCur;
				iss >> dCur;
				iss.clear();
				vecDatum.push_back(dCur);
			}
		}
	}	

	in.close();

	return true;
}

//	computing utilities
double PointDistance(double x1, double y1, double x2, double y2)
{
	double sqx = (x1 - x2) * (x1 - x2);
	double sqy = (y1 - y2) * (y1 - y2);

	return sqrt(sqx + sqy);
}

double Point2LineDistance(const Point &p, const Point &lineStart, const Point &lineEnd)
{
	//	Calculate line y = kx + m
	double k = (lineStart.second - lineEnd.second) / (lineStart.first - lineEnd.first);
	double m = lineStart.second - k * lineStart.first;

	//	Distance of point to line
	double dist = abs(k * p.first - p.second + m);
	dist /= sqrt(k * k + 1);

	return dist;
}

double Line2LinePerpendicularDistance(const Point &p1, const Point &p2, const Point &p3, const Point &p4)
{
	// Calculate perpendicular distance
	double perDist1 = Point2LineDistance(p1, p3, p4);
	double perDist2 = Point2LineDistance(p2, p3, p4);
	double perDist;
	if(0 == perDist1 && 0 == perDist2)
		perDist = 0;
	else
		perDist = ((perDist1 * perDist1) + (perDist2 + perDist2)) / (perDist1 + perDist2);

	return perDist;
}

double Line2LineAngleDistance(double perDist1, double perDist2, const Point &p1, const Point &p2, const Point &p3, const Point &p4)
{
	if(perDist1 == perDist2)
		return 0;

	double line1 = perDist1 > perDist2 ? (perDist1 - perDist2) : (perDist2 - perDist1);
	double lineOrg = PointDistance(p1.first, p1.second, p2.first, p2.second);
	double lineSeg = PointDistance(p3.first, p3.second, p4.first, p4.second);

	return lineSeg * line1/lineOrg;
}

double Line2LineDistance(const Point &p1, const Point &p2, const Point &p3, const Point &p4, double w1, double w2)
{
	double dPer = Line2LinePerpendicularDistance(p1, p2, p3, p4);
	double dAng = Line2LineAngleDistance(Point2LineDistance(p1, p3, p4), Point2LineDistance(p2, p3, p4), p1, p2, p3, p4);
	double l2lDist = 0.0;
	
	if(0 != dPer)
		l2lDist += w1 * log10(dPer) / log10(2.0);
	if(0 != dAng)
		l2lDist += w2 * log10(dAng) / log10(2.0);

	return  l2lDist;
}

double MDLseq(const deque<Point> &slideWindow, int startIndex, int curIndex)
{
	double LH = log10(PointDistance(slideWindow[startIndex].first, slideWindow[startIndex].second, 
						slideWindow[curIndex].first, slideWindow[curIndex].second)) / log10(2.0);
		
	double LDH = 0.0;
	
	size_t l = curIndex - startIndex;
	/*
	for(size_t i = startIndex; i < l; ++i)
	{
		LDH += Point2LineDistance(slideWindow[i], slideWindow[startIndex], slideWindow[curIndex]);
	}
	*/

	double angleDist = 0.0;
	double perpendicularDist = 0.0;
	Point p3 = slideWindow[startIndex];
	Point p4 = slideWindow[curIndex];
	for(size_t i = startIndex; i < l - 1; ++i)
	{
		Point p1 = slideWindow[i];
		Point p2 = slideWindow[i + 1];
		
		double perDist1 = Point2LineDistance(p1, p3, p4);
		double perDist2 = Point2LineDistance(p2, p3, p4);

		angleDist += Line2LineAngleDistance(perDist1, perDist2, p1, p2, p3, p4);
		perpendicularDist += Line2LinePerpendicularDistance(p1, p2, p3, p4);
	}



	LDH = log10(angleDist) / log10(2.0) + log10(perpendicularDist) / log10(2.0);



	return LH + LDH;
}

double MDLnoseq(const deque<Point> &slideWindow, int startIndex, int curIndex)		//	LDH = 0;
{
	return log10(PointDistance(slideWindow[startIndex].first, slideWindow[startIndex].second, 
						slideWindow[curIndex].first, slideWindow[curIndex].second)) / log10(2.0);
}

double MDLseq4Line(const deque<Point> &slideWindow, int startIndex, int curIndex)
{
	double LH = log10(PointDistance(slideWindow[startIndex].first, slideWindow[startIndex].second, 
						slideWindow[curIndex].first, slideWindow[curIndex].second)) / log10(2.0);
	
	double LDH = 0.0;

	


	return LH + LDH;
}

double MDLnoseq4Line(const deque<Point> &slideWindow, int startIndex, int curIndex)
{
	return MDLnoseq(slideWindow, startIndex, curIndex);
}

UTILITY_END