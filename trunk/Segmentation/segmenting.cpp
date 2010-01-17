#include "segmenting.h"

SEGMENTING_BEGIN

using namespace std;
using namespace utility;


/////////////////////////////////
//	Class Segmenting
/////////////////////////////////

void Segmenting::Merge(vector<Partition> &vecPartition, int index1, int index2)
{
	//	Head of first partition and Tail of second partition becomes the new partition
	Partition partition(vecPartition[index1].first, vecPartition[index2].second);

	vecPartition[index1] = partition;
	vector<Partition>::iterator itr = vecPartition.begin() + index2;

	//	delete the second partition
	vecPartition.erase(remove(itr, itr + 1, vecPartition[index2]));
}

long long Segmenting::GetSegmentSize() const
{
	return m_SegmentSize;
}


/////////////////////////////////
//	Class Parameterable
/////////////////////////////////

/////////////////////////////////
//	Class BottomUp
/////////////////////////////////

BottomUp::BottomUp(double dMaxError, const std::vector<double> vecDatum) 
	: ParameterMemorySegmenting(dMaxError, vecDatum)
{

}

double BottomUp::MergeCost(const Partition &p1, const Partition &p2)
{
	//	Calculate line y = kx + m
	double k = (p2.second.second - p1.first.second) / (p2.second.first - p1.first.first);
	double m = p1.first.second - k * p1.first.first;
	
	//	Distance of point 2 to line
	double dist1 = abs(k * p1.second.first - p1.second.second + m) / sqrt(k * k + 1);

	//	Distance of point 3 to line
	double dist2 = abs(k * p2.first.first - p2.first.second + m) / sqrt(k * k + 1);

	double cost = sqrt(dist1 * dist1 + dist2 * dist2);

	return cost;
}

vector<Point> BottomUp::Approximate()
{
	vector<Point> vecAppr;			//	approximate result
	vector<Partition> vecPartition;	//	first step partition
	vector<double> vecCost;			//	cost for each partition pair

	size_t minSize = this->m_vecDatum.size() / 100;
	minSize = minSize < 1 ? 1 : minSize;

	//	Generate fine approximation
	for(size_t i = 0; i < this->m_vecDatum.size() - 1; i += 2)
	{
		Point point1(i, this->m_vecDatum[i]);
		Point point2(i + 1, this->m_vecDatum[i + 1]);
		Partition partition(point1, point2);
		vecPartition.push_back(partition);
	}

	cout << "Fine partition size:" << vecPartition.size() << endl;

	//	Get cost for fine approximation
	int indexMin = 0;
	double dMinError = std::numeric_limits<double>::max();
	double avgCost = 0;
	for(size_t i = 0; i < vecPartition.size() - 1; ++i)
	{
		double cost = this->MergeCost(vecPartition[i], vecPartition[i + 1]);
		if(dMinError > cost)
		{
			dMinError = cost;
			indexMin = i;
		}
		vecCost.push_back(cost);
		avgCost += cost;
	}

	avgCost /= vecCost.size();

	cout << "AvgCost:" << avgCost << endl;

	//	Iterative till too coarse
	while(dMinError < this->m_MaxError)
	{
		if(0 == indexMin)		//	second partition been merged
		{
			this->Merge(vecPartition, indexMin, indexMin + 1);

			if(vecCost.size() != 1)
				vecCost[indexMin + 1] = this->MergeCost(vecPartition[indexMin], vecPartition[indexMin + 1]);

			vecCost.erase(remove(vecCost.begin() + indexMin, vecCost.end(), vecCost[indexMin]));
		}
		else if(vecCost.size() == indexMin)	//	last partition been merged
		{
			this->Merge(vecPartition, indexMin, indexMin + 1);

			if(vecCost.size() != 1)
				vecCost[indexMin - 1] = this->MergeCost(vecPartition[indexMin - 1], vecPartition[indexMin]);
			
			vecCost.erase(remove(vecCost.begin() + indexMin, vecCost.end(), vecCost[indexMin]));
		}
		else	//	middle partitions been merged
		{
			this->Merge(vecPartition, indexMin, indexMin + 1);
			vecCost[indexMin - 1] = this->MergeCost(vecPartition[indexMin - 1], vecPartition[indexMin]);

			if(vecCost.size() - 1 != indexMin)
				vecCost[indexMin + 1] = this->MergeCost(vecPartition[indexMin], vecPartition[indexMin + 1]);
			
			vecCost.erase(remove(vecCost.begin() + indexMin, vecCost.end(), vecCost[indexMin]));
		}
		//cout << "vecCost.size():" << vecCost.size() << ", cost:" << vecCost[indexMin] << endl;
		//vecCost[indexMin - 1] = this->MergeCost(vecPartition[indexMin - 1], vecPartition[indexMin]);

		int indexMinTmp = 0;
		double dMinErrorTmp = std::numeric_limits<double>::max();;
		for(size_t i = 0; i < vecCost.size(); ++i)
		{
			//double cost = this->MergeCost(vecPartition[i], vecPartition[i + 1]);
			if(dMinErrorTmp >= vecCost[i])
			{
				dMinErrorTmp = vecCost[i];
				indexMinTmp = i;
			}
		}
		dMinError = dMinErrorTmp;
		indexMin = indexMinTmp;

		if(minSize > vecPartition.size())
			break;
	}

	cout << "Coarse partition size:" << vecPartition.size() << endl; 

	ofstream out1("G:\\fine.txt");
	for(size_t i = 0; i < this->m_vecDatum.size(); ++i)
	{
		out1 << this->m_vecDatum[i] << endl;
	}
	out1.close();

	ofstream out2("G:\\coarse1.txt");
	for(size_t i = 0; i < vecPartition.size(); ++i)
	{
		out2 << vecPartition[i].first.first << " " << vecPartition[i].first.second << endl;
		out2 << vecPartition[i].second.first << " " << vecPartition[i].second.second << endl;
		Point p1(vecPartition[i].first.first, vecPartition[i].first.second);
		Point p2(vecPartition[i].second.first, vecPartition[i].second.second);
		vecAppr.push_back(p1);
		vecAppr.push_back(p2);
	}
	out2.close();

	return vecAppr;
}




/////////////////////////////////
//	Class ParameterSegmenting
/////////////////////////////////
ContinuousBottomUp::ContinuousBottomUp(double dMaxError, const std::vector<double> vecDatum) 
	: ParameterMemorySegmenting(dMaxError, vecDatum)
{
	m_Outputfile = "G:\\dataset\\TimeSeries\\output";
}

void ContinuousBottomUp::SetOutputFilepath(string filepath)
{
	m_Outputfile = filepath;
}

vector<Point> ContinuousBottomUp::Approximate()
{
	vector<Point> vecAppr;			//	approximate result
	vector<Partition> vecPartition;	//	first step partition
	vector<double> vecCost;			//	cost for each partition pair

	size_t minSize = this->m_vecDatum.size() / 100;
	minSize = minSize < 1 ? 1 : minSize;

	//	Generate fine approximation
	for(size_t i = 0; i < this->m_vecDatum.size() - 1; ++i)
	{
		Point point1(i, this->m_vecDatum[i]);
		Point point2(i + 1, this->m_vecDatum[i + 1]);
		Partition partition(point1, point2);
		vecPartition.push_back(partition);
	}

	cout << "Fine partition size:" << vecPartition.size() << endl;

	//	Get cost for fine approximation
	double avgCost = 0;
	int indexMin = 0;
	double dMinError = std::numeric_limits<double>::max();
	for(size_t i = 0; i < vecPartition.size() - 1; ++i)
	{
		double cost = this->MergeCost(vecPartition[i], vecPartition[i + 1]);
		if(dMinError > cost)
		{
			dMinError = cost;
			indexMin = i;
		}
		vecCost.push_back(cost);
		avgCost += cost;
	}

	avgCost /= vecCost.size();

	cout << "Average cost:" << avgCost << endl;

	//	Iterative till too coarse
	while(dMinError < this->m_MaxError)
	{
		//cout << "Merged " << vecPartition[indexMin + 1].first.first << " with min cost:" << vecCost[indexMin] << endl;

		this->Merge(vecPartition, indexMin, indexMin + 1);

		if(0 == indexMin)		//	second partition been merged
		{
			if(vecCost.size() != 1)	//	if it is also the last one
				vecCost[indexMin + 1] = this->MergeCost(vecPartition[indexMin], vecPartition[indexMin + 1]);

			vecCost.erase(remove(vecCost.begin() + indexMin, vecCost.end(), vecCost[indexMin]));
		}
		else if(vecCost.size() == indexMin)	//	last partition been merged
		{
			if(vecCost.size() != 1)	//	if it is also the first one
				vecCost[indexMin - 1] = this->MergeCost(vecPartition[indexMin - 1], vecPartition[indexMin]);
			
			vecCost.erase(remove(vecCost.begin() + indexMin, vecCost.end(), vecCost[indexMin]));
		}
		else	//	middle partitions been merged
		{
			vecCost[indexMin - 1] = this->MergeCost(vecPartition[indexMin - 1], vecPartition[indexMin]);

			if(vecCost.size() - 1 != indexMin)
				vecCost[indexMin + 1] = this->MergeCost(vecPartition[indexMin], vecPartition[indexMin + 1]);
			
			vecCost.erase(remove(vecCost.begin() + indexMin, vecCost.end(), vecCost[indexMin]));
		}

		//cout << "vecCost.size():" << vecCost.size() << ", cost:" << vecCost[indexMin] << endl;

		int indexMinTmp = 0;
		double dMinErrorTmp = std::numeric_limits<double>::max();
		for(size_t i = 0; i < vecCost.size(); ++i)
		{
			//double cost = this->MergeCost(vecPartition[i], vecPartition[i + 1]);
			if(dMinErrorTmp >= vecCost[i])
			{
				dMinErrorTmp = vecCost[i];
				indexMinTmp = i;
			}
		}
		dMinError = dMinErrorTmp;
		indexMin = indexMinTmp;

		if(minSize > vecPartition.size())
			break;
	}

	cout << "Coarse partition size:" << vecPartition.size() << endl; 

	/*
	ofstream out1("G:\\fine2.txt");
	for(size_t i = 0; i < this->m_vecDatum.size(); ++i)
	{
		out1 << this->m_vecDatum[i] << endl;
	}
	out1.close();
	*/


	if("NULL" != m_Outputfile)
	{
		ofstream out2(this->m_Outputfile.c_str());
		for(size_t i = 0; i < vecPartition.size(); ++i)
		{
			out2 << vecPartition[i].first.first << " " << vecPartition[i].first.second << endl;
			Point p(vecPartition[i].first.first, vecPartition[i].first.second);
			vecAppr.push_back(p);
		}
		out2 << vecPartition[vecPartition.size() - 1].second.first << " " << vecPartition[vecPartition.size() - 1].second.second << endl;
		Point p(vecPartition[vecPartition.size() - 1].second.first, vecPartition[vecPartition.size() - 1].second.second);
		vecAppr.push_back(p);
		out2.close();

		m_SegmentSize = vecAppr.size();
	}
	else
	{
		for(size_t i = 0; i < vecPartition.size(); ++i)
		{
			Point p(vecPartition[i].first.first, vecPartition[i].first.second);
			vecAppr.push_back(p);
		}

		Point p(vecPartition[vecPartition.size() - 1].second.first, vecPartition[vecPartition.size() - 1].second.second);
		vecAppr.push_back(p);

		m_SegmentSize = vecAppr.size();
	}

	return vecAppr;
}

double ContinuousBottomUp::MergeCost(const Partition &p1, const Partition &p2)
{
	//	Calculate line y = kx + m
	double k = (p2.second.second - p1.first.second) / (p2.second.first - p1.first.first);
	double m = p1.first.second - k * p1.first.first;
	
	//	Distance of point 2 to line
	double dist1 = abs(k * p1.second.first - p1.second.second + m) / sqrt(k * k + 1);

	double cost = sqrt(dist1 * dist1);

	return cost;
}


/////////////////////////////////
//	Class SlideWindowButtomUp
/////////////////////////////////
SlideWindowBottomUp::SlideWindowBottomUp(double dMaxError, size_t windowSize, istream *in, ofstream *out) 
	: Parameterable(dMaxError), m_WindowSize(windowSize), m_PartitionsPerTime(1), m_ApprSize(0),
		m_BatchSize(1), m_In(in), m_Out(out), m_CurrentCount(0)
{
	m_LowerBound = windowSize / 2;
	m_UpperBound = windowSize * 2;
	m_ApprSize = 0;
};

void SlideWindowBottomUp::SetBatchSize(size_t batchSize)
{
	m_BatchSize = batchSize;
}

void SlideWindowBottomUp::SetPartitionsPerTime(size_t partitionsPerTime)
{
	m_PartitionsPerTime = partitionsPerTime;
}

vector<Point> SlideWindowBottomUp::SlideBottomUp(vector<Partition> &vecPartition)
{
	//cout << "In bottom up-------------------------." << endl;
	/*
	for(size_t i = 0; i < vecPartition.size(); ++i)
	{
		cout << vecPartition[i].first.first << "," << vecPartition[i].first.second << endl;
		cout << vecPartition[i].second.first << "," << vecPartition[i].second.second << endl;
	}
	*/
	//cout << "---------------" << endl;

	vector<Point> vecAppr;
	vector<double> vecCost;
	size_t minSize = 2;

	//	Get cost for fine approximation
	double avgCost = 0;
	int indexMin = 0;
	double dMinError = std::numeric_limits<double>::max();
	
	for(size_t i = 0; i < vecPartition.size() - 1; ++i)
	{
		double cost = this->MergeCost(vecPartition[i], vecPartition[i + 1]);
		
		if(dMinError > cost)
		{
			dMinError = cost;
			indexMin = i;
		}
		
		vecCost.push_back(cost);
		avgCost += cost;
	}

	avgCost /= vecCost.size();
	
	//	Iterative till too coarse
	while(dMinError < this->m_MaxError)
	{
		//cout << "Merged " << vecPartition[indexMin + 1].first.first << " with min cost:" << vecCost[indexMin] << endl;

		this->Merge(vecPartition, indexMin, indexMin + 1);

		if(0 == indexMin)		//	second partition been merged
		{
			if(vecCost.size() != 1)	//	if it is also the last one
				vecCost[indexMin + 1] = this->MergeCost(vecPartition[indexMin], vecPartition[indexMin + 1]);

			vecCost.erase(remove(vecCost.begin() + indexMin, vecCost.end(), vecCost[indexMin]));
		}
		else if(vecCost.size() == indexMin)	//	last partition been merged
		{
			if(vecCost.size() != 1)	//	if it is also the first one
				vecCost[indexMin - 1] = this->MergeCost(vecPartition[indexMin - 1], vecPartition[indexMin]);
			
			vecCost.erase(remove(vecCost.begin() + indexMin, vecCost.end(), vecCost[indexMin]));
		}
		else	//	middle partitions been merged
		{
			vecCost[indexMin - 1] = this->MergeCost(vecPartition[indexMin - 1], vecPartition[indexMin]);

			if(vecCost.size() - 1 != indexMin)
				vecCost[indexMin + 1] = this->MergeCost(vecPartition[indexMin], vecPartition[indexMin + 1]);
			
			vecCost.erase(remove(vecCost.begin() + indexMin, vecCost.end(), vecCost[indexMin]));
		}

		//cout << "vecCost.size():" << vecCost.size() << ", cost:" << vecCost[indexMin] << endl;

		int indexMinTmp = 0;
		double dMinErrorTmp = std::numeric_limits<double>::max();
		for(size_t i = 0; i < vecCost.size(); ++i)
		{
			//double cost = this->MergeCost(vecPartition[i], vecPartition[i + 1]);
			if(dMinErrorTmp >= vecCost[i])
			{
				dMinErrorTmp = vecCost[i];
				indexMinTmp = i;
			}
		}
		dMinError = dMinErrorTmp;
		indexMin = indexMinTmp;

		if(minSize > vecPartition.size())
			break;
	}

	for(size_t i = 0; i < vecPartition.size(); ++i)
	{
		Point p(vecPartition[i].first.first, vecPartition[i].first.second);
		vecAppr.push_back(p);
	}

	//Point p(vecPartition[vecPartition.size() - 1].second.first, vecPartition[vecPartition.size() - 1].second.second);
	//vecAppr.push_back(p);

	//cout << "In bottom up:vecAppr.size()" << vecAppr.size() << endl;

	/*
	cout << "[";
	for(size_t i = 0; i < vecAppr.size(); ++i)
	{
		cout << "(" << vecAppr[i].first  << " " << vecAppr[i].second << ")";
	}
	cout << "]";
	*/

	return vecAppr;
}

vector<Point> SlideWindowBottomUp::Approximate()
{
	vector<Point> vecAppr;				//	store approximated point
	deque<Point> slideWindow;			//	slide window that stores points
	vector<Partition> slidePartition;	//	store partition in slide window

	istringstream iss;
	size_t size = 0;		//	size of current slide window

	//	initial slide window
	while((!m_In->eof()) && (size < this->m_WindowSize))
	{
		string line;
		getline(*m_In, line);

		if(string::npos != line.find(' '))	//	multiple data for each line
		{
			vector<string> vecDatum = utility::Split(line, ' ');
			for(size_t i = 0; i < vecDatum.size(); ++i)
			{
				double y;
				iss.str(vecDatum[i]);
				iss >> y;
				iss.clear();
				Point p(this->m_CurrentCount++, y);
				slideWindow.push_back(p);
				++size;
			}
		}
		else	//	single data for each line
		{
			double y;
			iss.str(line);
			iss >> y;
			iss.clear();
			Point p(this->m_CurrentCount++, y);
			slideWindow.push_back(p);
			++size;
		}
	}

	//	create partitions from slide window
	for(size_t i = 0; i < slideWindow.size() - 1; ++i)
	{
		Partition prt(slideWindow[i], slideWindow[i + 1]);
		slidePartition.push_back(prt);
	}

	//	process continuous comming data
	vector<Point> buffer;
	Index lastIndex = 0;
	long apprSize = 0;
	while(true == this->ReadIn(buffer, slideWindow.size()))
	{
		/*
		cout << "[";
		for(size_t i = 0; i < slideWindow.size(); ++i)
		{
			cout << "(" << slideWindow[i].first << "," << slideWindow[i].second << "),"; 
		}
		cout << "]" << endl;
		*/

		//cout << slidePartition.size() << endl;

		//	append new approximate partition to vecAppr
		vector<Point> vecCoarsePartition = this->SlideBottomUp(slidePartition);
		//cout << vecCoarsePartition.size() << endl;

		vector<Point>::iterator itr = vecCoarsePartition.begin();
		for(size_t i = 0; i < this->m_PartitionsPerTime + 1; ++i)
		{
			if(i < vecCoarsePartition.size())
			{
				*m_Out << itr->first << " " << itr->second << endl;
				//vecAppr.push_back(*itr++);
				itr++;
				apprSize++;
			}
		}


		//cout << "vecAppr.size():" << vecAppr.size() << endl;

		//cout << "before take out. size:" << slideWindow.size() << endl;
		this->TakeOut(slideWindow, vecCoarsePartition);
		//cout  << "after takeout. size:" << slideWindow.size() << endl;
		vecCoarsePartition.clear();	//	clear partition of slide window

		//cout << buffer.size() << endl;
		//	add new data into slide window
		for(size_t i = 0; i < buffer.size(); ++i)
		{
			slideWindow.push_back(buffer[i]);
		}

		//	update slidePartition according to slide window
		slidePartition.clear();
		for(size_t i = 0; i < slideWindow.size() - 1; ++i)
		{
			Partition prt(slideWindow[i], slideWindow[i + 1]);
			slidePartition.push_back(prt);
		}

		buffer.clear();	// clear buffer
	}

	(*m_Out).close();

	/*
	for(size_t i = 0; i < vecAppr.size(); ++i)
	{
		Point p = vecAppr[i];
		cout << p.first << "," << p.second << endl;
	}
	*/

	cout << "Approximate size:" << apprSize << endl;

	/*
	ofstream out("G:\\coarse3.txt");

	for(size_t i = 0; i < vecAppr.size(); ++i)
	{
		out << vecAppr[i].first << " " << vecAppr[i].second << endl;
	}
	*/

	m_ApprSize = apprSize;

	return vecAppr;
}

bool SlideWindowBottomUp::ReadIn(vector<Point> &buffer, int curWindowSize)
{
	size_t size = 0;
	istringstream iss;



	while(curWindowSize + size++ < m_WindowSize)
	{
		if(m_In->eof())
			return false;

		string line;
		getline(*m_In, line);

		if(string::npos != line.find(' '))	//	multiple data for each line
		{
			vector<string> vecDatum = utility::Split(line, ' ');
			for(size_t i = 0; i < vecDatum.size(); ++i)
			{
				double y;
				iss.str(vecDatum[i]);
				iss >> y;
				iss.clear();
				Point p(this->m_CurrentCount++, y);
				buffer.push_back(p);
				++size;
			}
		}
		else	//	single data for each line
		{
			double y;
			iss.str(line);
			iss >> y;
			iss.clear();
			Point p(this->m_CurrentCount++, y);
			buffer.push_back(p);
			++size;
		}
	}

	return true;
}

void SlideWindowBottomUp::TakeOut(deque<Point> &slideWindow, const vector<Point> &slidePartition)
{
	int offset = this->m_PartitionsPerTime < slidePartition.size() - 1 ? this->m_PartitionsPerTime : slidePartition.size() - 1;

	int x = slidePartition[offset].first;

	int beforeSize = slideWindow.size();
	

	//	iteratively delete the front of slideWindow until meet condition
	int swx = slideWindow.front().first;
	//cout << "x is:" << x  << " swx is:" << swx << endl;
	while(swx <= x)
	{
		//cout << "x is:" << x  << " swx is:" << swx << endl;
		slideWindow.erase(remove(slideWindow.begin(), slideWindow.end(), slideWindow[0]));
		swx = slideWindow.front().first;
	}

	//cout << "remove " << beforeSize - slideWindow.size() << " points from slide window." << endl;

}

double SlideWindowBottomUp::MergeCost(const Partition &p1, const Partition &p2)
{
	//	Calculate line y = kx + m
	double k = (p2.second.second - p1.first.second) / (p2.second.first - p1.first.first);
	double m = p1.first.second - k * p1.first.first;
	
	//	Distance of point 2 to line
	double dist1 = abs(k * p1.second.first - p1.second.second + m) / sqrt(k * k + 1);

	double cost = sqrt(dist1 * dist1);

	return cost;
}

long long SlideWindowBottomUp::GetSize() const
{
	return m_CurrentCount;
}

long long SlideWindowBottomUp::GetApprSize() const
{
	return m_ApprSize;
}


/////////////////////////////////
//	Class MDLSlideWindow
/////////////////////////////////
MDLSlideWindow::MDLSlideWindow(size_t windowSize, std::istream *in, std::ofstream *out) 
	: m_WindowSize(windowSize), m_In(in), m_Out(out), m_BatchSize(1), m_CurrentSize(0)
{
	m_BatchSize = 1;
}

void MDLSlideWindow::SetBatchSize(size_t batchSize)
{
	if(batchSize > 0 && batchSize < m_WindowSize)
		m_BatchSize = batchSize;
}

vector<Point> MDLSlideWindow::MDLSegmenting()
{
	vector<Point> appr;

	if(!m_SlideWindow.empty())
		appr.push_back(m_SlideWindow.front());

	size_t startIndex = 0;
	size_t length = 1;
	size_t curIndex = 0;

	size_t segFound = 0;
	size_t windowSize = m_SlideWindow.size();

	while(startIndex + length < windowSize)
	{
		curIndex = startIndex + length;
		double costseq = MDLseq(m_SlideWindow, startIndex, curIndex);
		double costnoseq = MDLnoseq(m_SlideWindow, startIndex, curIndex);
		//cout << startIndex << " to " << curIndex << ", seq:\t" << costseq << endl;
		//cout << startIndex << " to " << curIndex << ", noseq:\t" << costnoseq << endl;
		if(costseq > costnoseq)	//	find the first segment
		{
			appr.push_back(m_SlideWindow[curIndex]);
			++segFound;
			
			if(segFound == m_BatchSize)
				return appr;
			
			startIndex = curIndex;
			length = 1;
		}
		else
		{
			++length;
		}
	}

	return appr;
}

void MDLSlideWindow::Approximate()
{
	cout << "Begin to approximate..." << endl;
	vector<Point> appr;
	long long apprSize = 0;
	istringstream iss;
	
	//	fill in slide window
	while((!m_In->eof()) && (m_SlideWindow.size() < m_WindowSize))	//	continuous reading in while new data arrive
	{
		string line;
		getline(*m_In, line);

		if(line.find(' '))	//	multiple data per line
		{
			vector<string> vecValues = Split(line, ' ');

			for(size_t i = 0; i < vecValues.size(); ++i)
			{
				iss.str(vecValues[i]);
				double value;
				iss >> value;
				iss.clear();
				Point p(m_CurrentSize, value);
				m_SlideWindow.push_back(p);
				++m_CurrentSize;
			}
		}
		else	//	one data per line
		{
			iss.str(line);
			double value;
			iss >> value;
			iss.clear();
			
			Point p(m_CurrentSize, value);
			m_SlideWindow.push_back(p);
			++m_CurrentSize;
		}
	}


	//	stop until no data in
	//ofstream out("G:\\coarsemdl.txt");
	Point lastPoint;
	int jump = 1000000;
	int batch = 0;
	clock_t start, end;
	start = clock();
	int j = 0;
	while(true == ReadIn())
	{
		vector<Point> segPoints = MDLSegmenting();


		/*
		for(size_t i = 0; i < segPoints.size(); ++i)
		{
			cout << segPoints[i].first << "," << segPoints[i].second << endl;
		}
		cout << "end " << j++ << endl;
		*/
		
		//cout << appr.size() << endl;

		/*
		for(size_t i = 0; i < m_SlideWindow.size(); ++i)
		{	
			cout << "[" << m_SlideWindow[i].first << "," << m_SlideWindow[i].second << "],";
			//if(appr.empty() || m_SlideWindow[i] != appr[appr.size() - 1])
			//	appr.push_back(m_SlideWindow[i]);
		}
		cout << endl;
		*/

		if(2 <= segPoints.size())
		{
			Point last = segPoints[segPoints.size() - 1];

			for(size_t i = 0; i < segPoints.size(); ++i)
			{
				if(segPoints[i] != lastPoint)
					appr.push_back(segPoints[i]);
			}

			TakeOut(last);
			//cout << "take out" << endl;
		}
		else if(1 == segPoints.size())	//	data cannot be compressed
		{
			for(size_t i = 0; i < m_SlideWindow.size(); ++i)
			{	
				//cout << "[" << m_SlideWindow[i].first << "," << m_SlideWindow[i].second << "],";
				if(m_SlideWindow[i] != lastPoint)
					appr.push_back(m_SlideWindow[i]);
			}
			//cout << endl;

			m_SlideWindow.clear();
		}

		for(size_t i = 0; i < appr.size(); ++i)
		{
			*m_Out << appr[i].first << " " << appr[i].second << endl;
		}	

		lastPoint = appr[appr.size() - 1];
		apprSize += appr.size();
		appr.clear();

		if(batch < m_CurrentSize / jump)
		{
			++batch;
			end = clock();
			cout << m_CurrentSize << " records has been segmented for " << apprSize << " segments, time elapse:" << (end - start) / CLOCKS_PER_SEC << endl;
		}

	}

	cout << "Total size:" << m_CurrentSize << endl;
	cout << "Approximate size:" << apprSize << endl;

	m_ApprSize = apprSize;

	m_Out->close();
}

bool MDLSlideWindow::ReadIn()	//	fill the slide window
{
	if(m_SlideWindow.size() >= m_WindowSize)
		return true;

	istringstream iss;

	while(m_SlideWindow.size() < m_WindowSize)
	{
		if(m_In->eof())
			return false;

		string line;
		getline(*m_In, line);

		if(line.find(' '))	//	multiple data per line
		{
			vector<string> vecValues = Split(line, ' ');

			for(size_t i = 0; i < vecValues.size(); ++i)
			{
				iss.str(vecValues[i]);
				double value;
				iss >> value;
				iss.clear();
				Point p(m_CurrentSize, value);
				m_SlideWindow.push_back(p);
				++m_CurrentSize;
			}
		}
		else	//	one data per line
		{
			iss.str(line);
			double value;
			iss >> value;
			iss.clear();
			
			Point p(m_CurrentSize, value);
			m_SlideWindow.push_back(p);
			++m_CurrentSize;
		}
	}

	return true;
}

void MDLSlideWindow::TakeOut(const Point &last)
{
	Index index = last.first;

	while(true)
	{
		Point p = m_SlideWindow.front();
		Index idx = p.first;
		if(idx >= index)
		{
			return;
		}
		else
		{
			m_SlideWindow.pop_front();
		}
		;
	}
}

long long MDLSlideWindow::GetSize() const
{
	return m_CurrentSize;
}

long long MDLSlideWindow::GetApprSize() const
{
	return m_ApprSize;
}

SEGMENTING_END