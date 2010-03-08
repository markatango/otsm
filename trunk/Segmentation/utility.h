#pragma once

#include "headers.h"

UTILITY_BEGIN

using std::vector;
using std::deque;
using std::string;
using std::ifstream;
using std::istringstream;
using std::ofstream;
using std::endl;
using std::cout;

//	Data generator
void Generator(double min, double max, long long amount, string path);
void SineGenerator(long long amount, string path, double a, double w, double b);
void RandomDataGenerator(long long amount, string path, double min, double max, double sharpness = 0.001);

vector<string> Split(const string& strLine, char delim);
bool ReadDataFromFile(vector<double>& vecDatum, const string& strFilename, char delim = ' ', int iStartPos = 0);

//	computing utilities
double PointDistance(double x1, double y1, double x2, double y2);


//	Distance from point p to line (lineStart, lineEnd)
double Point2LineDistance(const Point &p, const Point &lineStart, const Point &lineEnd);

//	Angle distance from line(p1,p2) to line(p3,p4)
double Line2LineAngleDistance(double perDist1, double perDist2, const Point &p1, const Point &p2, const Point &p3, const Point &p4);

//	Perpendicular distance
double Line2LinePerpendicularDistance(const Point &p1, const Point &p2, const Point &p3, const Point &p4);

//	Line to line distance
double Line2LineDistance(const Point &p1, const Point &p2, const Point &p3, const Point &p4, double w1 = 1.0, double w2 = 1.0);


//	MDL related computing
double MDLseq(const deque<Point> &slideWindow, int startIndex, int curIndex);
double MDLnoseq(const deque<Point> &slideWindow, int startIndex, int curIndex);

double MDLseq4Line(const deque<Point> &slideWindow, int startIndex, int curIndex);
double MDLnoseq4Line(const deque<Point> &slideWindow, int startIndex, int curIndex);

UTILITY_END