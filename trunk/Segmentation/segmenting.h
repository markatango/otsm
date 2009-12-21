#pragma once

#include "headers.h"
#include "utility.h"

SEGMENTING_BEGIN


/**
*	Base class Segmenting.
*/
class Segmenting
{
public:
	Segmenting(){}; 
	void Merge(vector<Partition> &vecPartition, int index1, int index2);
	virtual vector<Point> Approximate() = 0;
	
};

class Parameterable
{
public:
	Parameterable(double dMaxError) : m_MaxError(dMaxError){};

protected:
	virtual double MergeCost(const Partition &p1, const Partition &p2) = 0;

protected:
	double m_MaxError;
};

class MemoryBased
{
public:
	MemoryBased(vector<double> vecDatum) : m_vecDatum(vecDatum){};
	
protected:
	vector<double> m_vecDatum;
};


class ParameterMemorySegmenting : public Parameterable, public MemoryBased, public Segmenting
{
public:
	ParameterMemorySegmenting(double dMaxError, vector<double> vecDatum) : Parameterable(dMaxError), MemoryBased(vecDatum){};

};


/**
/	Class BottomUp. Approximate line every two points.
*/
class BottomUp : public ParameterMemorySegmenting
{
public:
	BottomUp(double dMaxError, const vector<double> vecDatum);
	vector<Point> Approximate();

protected:
	double MergeCost(const Partition &p1, const Partition &p2);
};

/**
*	Class ContinuousButtomUp. Appcoximate line with continuous points.
*/
class ContinuousBottomUp : public ParameterMemorySegmenting
{
public:
	ContinuousBottomUp(double dMaxError, const vector<double> vecDatum);
	vector<Point> Approximate();

protected:
	double MergeCost(const Partition &p1, const Partition &p2);
};

/**
*	Class SlideWindowBottomUp. Appoximate line with slide window.
*/
class SlideWindowBottomUp : public Segmenting, public Parameterable
{
public:
	SlideWindowBottomUp(double dMaxError, size_t windowSize, std::istream *in, std::ofstream *out);
	vector<Point> Approximate();
	void SetBatchSize(size_t batchSize);
	void SetPartitionsPerTime(size_t partitionsPerTime); 
	long long GetSize() const;
	long long GetApprSize() const;


protected:
	double MergeCost(const Partition &p1, const Partition &p2);
	vector<Point> SlideBottomUp(vector<Partition> &vecPartition);

private:
	bool ReadIn(vector<Point> &buffer, int curWindowSize);
	void TakeOut(deque<Point> &slideWindow, const vector<Point> &slidePartition);

protected:
	size_t m_BatchSize;			//	the size of data load in per time

	size_t m_WindowSize;
	size_t m_PartitionsPerTime;	//	add how many partition during each calculationg of slide window
	std::istream *m_In;
	std::ofstream *m_Out;

	long m_LowerBound;
	long m_UpperBound;
	Index m_CurrentCount;
	long long m_ApprSize;

	vector<Point> m_Seg;

};

/**
*	Class MDLSlideWindow.
*
*/
class MDLSlideWindow
{
public:
	MDLSlideWindow(size_t windowSize, std::istream *in, std::ofstream *out);
	void SetBatchSize(size_t batchSize);

	void Approximate();
	long long GetSize() const;
	long long GetApprSize() const;

protected:
	vector<Point> MDLSegmenting();
	bool ReadIn();
	void TakeOut(const Point &last);

private:
	size_t m_WindowSize;		//	size of slide window
	size_t m_BatchSize;			//	number of segment generate per time at most
	std::istream *m_In;
	std::ofstream *m_Out;
	deque<Point> m_SlideWindow;

	Index m_CurrentSize;			//	size of data has been read
	long long m_ApprSize;
};

SEGMENTING_END

