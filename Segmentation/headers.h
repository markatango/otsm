#pragma once

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

using std::pair;
using std::vector;
using std::queue;
using std::deque;


//	define namespace
#define UTILITY_BEGIN namespace utility{ 
#define UTILITY_END } 

#define SEGMENTING_BEGIN namespace segmenting{ 
#define SEGMENTING_END } 

#define MOTIF_BEGIN namespace motif{
#define MOTIF_END }

typedef long Index;
typedef pair<Index, double> Point;		//	first represents x, second represents y.
typedef pair<Point, Point> Partition;


