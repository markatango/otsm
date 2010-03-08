#pragma once

#include "headers.h"

MOTIF_BEGIN

/**
*	Definition of Candidate Motif Buffer
*/
typedef std::vector<double> TimeSeries;
typedef std::pair<long long, TimeSeries> TimeSeriesCandidate;
typedef std::vector<TimeSeriesCandidate> CMB;

/**
*	Distance measurement.
*/
inline double EuclideanDistance(TimeSeries &t1, TimeSeries &t2, std::string name = "");

struct Parameters
{
	unsigned long slideWindowSize;
	unsigned long motifLength;
	double R;
};

struct ParametersNaive : public Parameters
{
	unsigned long step;
};

struct ParametersSIMMA : public Parameters
{
	unsigned long k;
	double sigma;
};

/**
*	Definition of Brute Force Find Motif Algorithm.
*/
class BFFindMotif
{
public:
	BFFindMotif(unsigned long windowSize, unsigned long len, double R, unsigned long step, std::istream *in);
	~BFFindMotif();
	
	//	setters
	void SetParameters(const ParametersNaive &params);
	void SetSlideWindowSize(unsigned long windowSize);
	void SetMotifLength(unsigned long len);
	void SetRadius(double R);
	void SetStep(unsigned long step);
	void SetIStream(std::istream *in);

	std::priority_queue<TimeSeriesCandidate> GetTopKMotifs(size_t k);

	void FindMotif();
	Index GetLastIndex() const;
	size_t GetCandidateMotifSize() const;


protected:
	void FillWindow();
	void FindMotifSub(deque<Point> &window);


protected:

	//	properties
	unsigned long m_SlideWindowSize;	//	The size of slide window, default 1024
	unsigned long m_MotifLength;		//	The size of motif length, should be less than size of slide window, default 1/4 of slide window
	double m_Radius;		//	The distance radius to judge whether two time series are different. 2R < dist(t1, t2) means new motif appear, R > dist(t1, t2) means t1, t2 are similar
	unsigned long m_Step;			//	The step size for time series stream, if new motif or similar motifs are found , the next search potision should jump forward R location for fear of trivial match

	Index m_CurIndex;
	std::istream *m_In;
	CMB m_CandidateMotif;
	std::deque<Point> m_SlideWindow;
};





/**
*	Definition of SIMMA algorithm
*/
class SIMMA
{
public:
	SIMMA(unsigned long windowSize, unsigned long len, double R, std::istream *in, double sigma, long k);
	~SIMMA();

	//	setters
	void SetParameters(const ParametersSIMMA &params);
	void SetSlideWindowSize(unsigned long windowSize);
	void SetMotifLength(unsigned long len);
	void SetRadius(double R);
	void SetK(unsigned long k);
	void SetIStream(std::istream *in);
	void SetBufferCheck(bool check);

	std::priority_queue<TimeSeriesCandidate> GetTopKMotifs(size_t k);
	void FindMotif();
	size_t GetCandidateMotifSize() const;
	Index GetLastIndex() const;
	size_t GetShrinkCount() const;

protected:
	void FillWindow();
	void FindMotifSub(long &bufferCount);
	long NonTrivialStep(const long curPos);
	void BufferCheck();

protected:
	//	properties
	unsigned long m_SlideWindowSize;	//	The size of slide window, default 1024
	unsigned long m_MotifLength;		//	The size of motif length, should be less than size of slide window, default 1/4 of slide window
	double m_Radius;		//	The distance radius to judge whether two time series are different. 2R < dist(t1, t2) means new motif appear, R > dist(t1, t2) means t1, t2 are similar
	double m_Sigma;			//	Used to check whether buffer should be shrinked
	long m_K;				//	Number of top k motifs
	bool m_bBufferCheck;		//	Whether check buffer periodically
	size_t m_ShrinkCount;

	Index m_CurIndex;
	std::istream *m_In;
	CMB m_CandidateMotif;
	std::deque<Point> m_SlideWindow;
};

MOTIF_END