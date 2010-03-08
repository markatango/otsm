#include "motif.h"

using namespace std;

MOTIF_BEGIN

inline double EuclideanDistance(TimeSeries &t1, TimeSeries &t2, string name)
{
	double dist = 0.0;

	if(t1.size() != t2.size())
	{
		cout << t1.size() << "," << t2.size() << endl;
		cout << "Called by " << name << endl;
		system("pause");
	}

	for(size_t i = 0; i < t1.size(); ++i)
	{
		dist += pow(t1[i] - t2[i], 2.0);
	}

	return sqrt(dist);
}

BFFindMotif::BFFindMotif(unsigned long windowSize, unsigned long len, double R, unsigned long step, istream *in) 
	: m_In(in), m_CurIndex(0)
{
	SetSlideWindowSize(windowSize);
	SetMotifLength(len);
	SetRadius(R);
	SetStep(step);
}

BFFindMotif::~BFFindMotif()
{
	m_SlideWindow.clear();
	m_CandidateMotif.clear();
}

void BFFindMotif::SetSlideWindowSize(unsigned long windowSize)
{
	m_SlideWindowSize = windowSize;
}

void BFFindMotif::SetMotifLength(unsigned long len)
{
	if(0 == len || len >= 0.5 * m_SlideWindowSize)
		m_MotifLength = static_cast<long>(0.25 * m_SlideWindowSize);
	else
		m_MotifLength = len;
}

void BFFindMotif::SetRadius(double R)
{
	if(0.0 >= R)
		m_Radius = 0.1;
	else
		m_Radius = R;
}

void BFFindMotif::SetStep(unsigned long step)
{
	if(0.125 * m_MotifLength >= step || m_MotifLength < step)
		m_Step = m_MotifLength;
	else
		m_Step = step;
}

void BFFindMotif::SetParameters(const ParametersNaive &params)
{
	SetSlideWindowSize(params.slideWindowSize);
	SetMotifLength(params.motifLength);
	SetRadius(params.R);
	SetStep(params.step);
}

void BFFindMotif::SetIStream(istream *in)
{
	m_In = in;
}

priority_queue<TimeSeriesCandidate> BFFindMotif::GetTopKMotifs(size_t k)
{
	priority_queue<TimeSeriesCandidate> topkTS;
	priority_queue<TimeSeriesCandidate> pq;

	for(size_t i = 0; i < m_CandidateMotif.size(); ++i)
	{
		pq.push(m_CandidateMotif[i]);
	}

	for(size_t i = 0; i < k; ++i)
	{
		if(pq.empty())
			break;

		topkTS.push(pq.top());
		pq.pop();
	}
	
	return topkTS;
}

void BFFindMotif::FindMotif()
{
	//	Fill up slide window
	FillWindow();

	//	Find motifs through FindMotifSub
	while(false == m_In->eof())
	{
		this->FindMotifSub(m_SlideWindow);
		this->FillWindow();
	}
}

void BFFindMotif::FillWindow()
{
	istringstream iss;

	while((false == m_In->eof()) && (m_SlideWindow.size() < m_SlideWindowSize))
	{
		string line;
		double value;
		getline(*m_In, line);
		iss.str(line);
		iss >> value;
		iss.clear();

		m_SlideWindow.push_back(make_pair<Index, double>(++m_CurIndex, value));

		++m_CurIndex;

		
		if(m_CurIndex % 5000000 == 0)
		{
			cout << m_CurIndex << " records read." << endl;
		}
		
	}
}

void BFFindMotif::FindMotifSub(std::deque<Point> &window)
{
	size_t i = 0;
	for(i = 0; i < window.size() - m_MotifLength; ++i)
	{
		double distance = 0.0;
		bool newMotif = true;
		TimeSeries ts;
		ts.reserve(m_MotifLength);
		if(m_SlideWindow.size() >= m_MotifLength)		//	Only process slide window larger than motif length
		{
			//	Get time series
			for(size_t j = i; j < i + m_MotifLength; ++j)
			{
				ts.push_back(window[j].second);
			}

			//	Compare with candidate motif
			for(size_t j = 0; j < m_CandidateMotif.size(); ++j)
			{
				distance = EuclideanDistance(m_CandidateMotif[j].second, ts);

				if((2 * m_Radius > distance) && (m_Radius < distance))	//	Neither new motif nor similar motif
				{
					newMotif = false;
				}
				else if(m_Radius > distance)	//	Similar motif
				{
					m_CandidateMotif[j].first++;
					i += m_Step;
					newMotif = false;
					break;		//	Impossible to be similar with other candidates
				}
			}

			if(true == newMotif)	//	New motif
			{
				m_CandidateMotif.push_back(make_pair<long long, TimeSeries>(1, ts));
				i += m_Step;
			}
			
		}
		else
		{
			cerr << "Window size:" << m_SlideWindow.size() << endl;
		}
	
	}
	
	for(size_t k = 0; k < i; ++k)
	{
		window.pop_front();
	}
	
}

Index BFFindMotif::GetLastIndex() const
{
	return m_CurIndex;
}

size_t BFFindMotif::GetCandidateMotifSize() const
{
	return m_CandidateMotif.size();
}





/**
*	SIMMA
*/
SIMMA::SIMMA(unsigned long windowSize, unsigned long len, double R, std::istream *in, double sigma, long k)
	: m_SlideWindowSize(windowSize), m_MotifLength(len), m_Radius(R), m_In(in), 
	m_Sigma(sigma), m_K(k), m_CurIndex(0), m_CandidateMotif(), m_SlideWindow(), m_bBufferCheck(true), m_ShrinkCount(0)
{

}

SIMMA::~SIMMA()
{
	m_SlideWindow.clear();
	m_CandidateMotif.clear();
}

void SIMMA::SetSlideWindowSize(unsigned long windowSize)
{
	m_SlideWindowSize = windowSize;
}

void SIMMA::SetMotifLength(unsigned long len)
{
	if(0 == len || len >= 0.5 * m_SlideWindowSize)
		m_MotifLength = static_cast<long>(0.25 * m_SlideWindowSize);
	else
		m_MotifLength = len;
}

void SIMMA::SetRadius(double R)
{
	if(0.0 >= R)
		m_Radius = 0.1;
	else
		m_Radius = R;
}

void SIMMA::SetK(unsigned long k)
{
	m_K = k;
}

void SIMMA::SetParameters(const ParametersSIMMA &params)
{
	SetSlideWindowSize(params.slideWindowSize);
	SetMotifLength(params.motifLength);
	SetRadius(params.R);
	SetK(params.k);
}

void SIMMA::SetIStream(istream *in)
{
	m_In = in;
}

void SIMMA::SetBufferCheck(bool check)
{
	m_bBufferCheck = check;
}

size_t SIMMA::GetCandidateMotifSize() const
{
	return m_CandidateMotif.size();
}

Index SIMMA::GetLastIndex() const
{
	return m_CurIndex;
}

size_t SIMMA::GetShrinkCount() const
{
	return m_ShrinkCount;
}

std::priority_queue<TimeSeriesCandidate> SIMMA::GetTopKMotifs(size_t k)
{
	if(k >= m_CandidateMotif.size() * 0.00000002)
	{
		priority_queue<TimeSeriesCandidate> topkTS;
		priority_queue<TimeSeriesCandidate> pq;

		for(size_t i = 0; i < m_CandidateMotif.size(); ++i)
		{
			pq.push(m_CandidateMotif[i]);
		}

		for(size_t i = 0; i < k; ++i)
		{
			if(pq.empty())
				break;

			topkTS.push(pq.top());
			pq.pop();
		}

		return topkTS;
	}
	else
	{
		//	use k-means to cluster top 20% motifs
		
	}
}

void SIMMA::FindMotif()
{
	//	Fill up slide window
	FillWindow();

	long bufferCount = 0;
	//	Find motifs through FindMotifSub
	while(false == m_In->eof())
	{
		this->FindMotifSub(bufferCount);
		this->FillWindow();
	}
}

void SIMMA:: FindMotifSub(long &bufferCount)
{
	long lLastSize = m_K;			//	last size of candidate motif buffer
	long step = 32;

	size_t i = 0;
	for(i = 0; i < m_SlideWindow.size() - m_MotifLength; ++i)
	{
		double distance = 0.0;
		bool newMotif = true;
		TimeSeries ts;
		ts.reserve(m_MotifLength);

		++bufferCount;

		//	Check buffer when size big enough
		if(true == m_bBufferCheck)
		{
			if((2 * m_K < bufferCount) && (5 * m_K < (int)m_CandidateMotif.size()))
			{
				this->BufferCheck();
				lLastSize = m_SlideWindow.size();
				bufferCount = 0;
			}
		}

		if(m_SlideWindow.size() >= m_MotifLength)	//	Only process slide window larger than motif length
		{
			//	Get time series
			for(size_t j = i; j < i + m_MotifLength; ++j)
			{
				ts.push_back(m_SlideWindow[j].second);
			}

			//	Compare with candidate motif
			for(size_t j = 0; j < m_CandidateMotif.size(); ++j)
			{
				distance = EuclideanDistance(m_CandidateMotif[j].second, ts);

				if((2 * m_Radius > distance) && (m_Radius < distance))	//	Neither new motif nor similar motif
				{
					newMotif = false;
				}
				else if(m_Radius > distance)	//	Similar motif
				{
					m_CandidateMotif[j].first++;
					long jump = this->NonTrivialStep(i);

					if(0 != jump)
						i += jump;
					else
						i = m_SlideWindow.size() - m_MotifLength;

					newMotif = false;
					break;		//	Impossible to be similar with other candidates
				}
			}
			
			//	Check whether current time series is new motif
			if(true == newMotif)
			{
				m_CandidateMotif.push_back(make_pair<long long, TimeSeries>(1, ts));
				long jump = this->NonTrivialStep(i);

				if(0 != jump)
					i += jump;
				else
					i = m_SlideWindow.size() - m_MotifLength;
			}
		}
	}
	
	//	Remove used elements
	for(size_t k = 0; k < i; ++k)
	{
		m_SlideWindow.pop_front();
	}
}

void SIMMA::FillWindow()
{
	istringstream iss;

	while((false == m_In->eof()) && (m_SlideWindow.size() < m_SlideWindowSize))
	{
		string line;
		double value;
		getline(*m_In, line);
		iss.str(line);
		iss >> value;
		iss.clear();

		m_SlideWindow.push_back(make_pair<Index, double>(++m_CurIndex, value));

		++m_CurIndex;

		
		if(m_CurIndex % 5000000 == 0)
		{
			cout << m_CurIndex << " records read." << endl;
		}
		
	}
}

long SIMMA::NonTrivialStep(const long curPos)
{
	long lStep = 2;							//	default step distance
	long lCur = curPos;						//	current position of sub time series sequence
	long lNext = lCur + lStep / 2;
	long motifLen = m_MotifLength;			//	length of motif to be discovered
	long windowLen = m_SlideWindow.size();	//	length of slide window

	TimeSeries ts1;	//	current sub time series sequence
	for(int j = lCur; j < lCur + motifLen; ++j)
	{
		ts1.push_back(m_SlideWindow[j].second);
	}

	TimeSeries ts2;	//	next sub time series sequence
	for(int pos = lNext; pos < lNext + motifLen; ++pos)	//	next sub time series sequence can be get by copy from current sub time series sequence
	{
		ts2.push_back(m_SlideWindow[pos].second);
	}

	double distance = EuclideanDistance(ts1, ts2);
	if(distance > m_Radius)
		return lStep;

	long lTmp = lStep;
	lStep *= 2;
	lStep += lTmp;
	lNext = lCur + lStep / 2;

	while(lCur + lStep + motifLen < windowLen)
	{
		
		if(lNext + 1 + lStep > windowLen)	//	out of range
			return 0;
		else
		{
			//	update ts2, pop first step points and push new step point into st2
			for(int cand = 0; cand < lStep; ++cand)
			{
				ts2.erase(ts2.begin());
				int index = lNext + 1 + cand;
				ts2.push_back(m_SlideWindow[index].second);
			}
		}
		
		double distance = EuclideanDistance(ts1, ts2);
		if(distance > m_Radius)
			return lStep;

		lTmp = lStep;
		lStep *= 2;
		lStep += lTmp;
		lNext = lStep / 2;
	}

	return 0;

}

void SIMMA::BufferCheck()
{
	bool shrink = true;

	long double frac = 0.0;	//	get sum of frequency
	long long sum = 0;

	size_t before = m_CandidateMotif.size();
	
	for(size_t i = 0; i < m_CandidateMotif.size(); ++i)
	{
		frac += 1 / m_CandidateMotif[i].first;
		sum += m_CandidateMotif[i].first;
	}

	for(int i = 0; i < m_K; ++i)
	{
		double f = 1 / (m_K * frac);
		double p = double(m_CandidateMotif[i].first / sum);

		if(f - p > m_Sigma)	//	if any one of the probability exceeds the error threshold, don't shrink
		{
			shrink = false;
			break;
		}
	}

	CMB tmp;
	if(true == shrink)
	{
		++m_ShrinkCount;
		for(int i = 0; i < m_K; ++i)
		{
			tmp.push_back(m_CandidateMotif[i]);
		}
		m_CandidateMotif = tmp;

		cout << "shrink from " << before << " to " << m_CandidateMotif.size() << endl;
	}

}


MOTIF_END