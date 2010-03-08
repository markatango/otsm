#include "experiments.h"

using namespace std;
using namespace motif;

void Test()
{
	cout << "Begin to test finding motif." << endl;

	string infile("G:\\dataset\\TimeSeries\\sythetic\\sythetic2");
	ifstream *in = new ifstream(infile.c_str());

	if(true == in->bad())
	{
		cerr << "Input file is bad!" << endl;
		return;
	}


	unsigned long windowSize = 512;
	unsigned long motifLength = 64;
	unsigned long step = 16;
	double R = 20;
	size_t k = 126;

	BFFindMotif *findMotif = new BFFindMotif(windowSize, motifLength, R, step, in);
	clock_t start, end;
	start = clock();
	findMotif->FindMotif();
	end = clock();
	priority_queue<TimeSeriesCandidate> topk = findMotif->GetTopKMotifs(k);

	cout << "There are totally " << findMotif->GetLastIndex() << " sample points." << endl;
	cout << "Find out " << findMotif->GetCandidateMotifSize() << " candidate motifs." << endl;
	cout << "Time cost:" << (end - start) / CLOCKS_PER_SEC << " seconds." << endl;

	ofstream *out = new ofstream("G:\\dataset\\TimeSeries\\real\\motif");

	size_t motifSize = k < topk.size()? k : topk.size();

	/*
	for(size_t i = 0; i < motifSize; ++i)
	{
	TimeSeriesCandidate ts = topk.top();
	topk.pop();
	cout << "[" << ts.first << "]:";
	for(size_t j = 0; j < ts.second.size(); ++j)
	{
	//cout << (ts.second)[j] << "\t";
	*out << (ts.second)[j] << "\t";
	}
	*out << endl; 
	//cout << endl;
	}
	*/
	for(size_t i = 0; i < motifSize; ++i)
	{
		TimeSeriesCandidate ts = topk.top();
		topk.pop();
		*out << ts.first << "\t";
	}
	out->close();

	cout << "End of test finding motif." << endl;
}

/**
*	Naive related experiments
*/

void NaiveLongTailExperiments()
{
	cout << "-------Naive long tail test----------" << endl;
	//	Test for real datasets
	string inFilePrefix = "G:\\dataset\\TimeSeries\\real";
	string outFilePrefix = "G:\\";
	ParametersNaive params;
	params.slideWindowSize = 512;
	params.motifLength = 64;
	params.R = 30;
	params.step = 16;


	//	IBM data
	NaiveLongTailExperiment(inFilePrefix + "\\IBM_data", outFilePrefix + "\\IBM_lt", params);
	//	Dodger
	NaiveLongTailExperiment(inFilePrefix + "\\Traffic_data", outFilePrefix + "\\Traffic_lt", params);
	//	ICU
	NaiveLongTailExperiment(inFilePrefix + "\\ICU_data", outFilePrefix + "\\ICU_lt", params);


	/*
	//	synthetic
	string inSyFilePrefix = "G:\\dataset\\TimeSeries\\sythetic";
	for(size_t i = 9; i < 10; ++i)
	{
	stringstream ss;
	ss << i;
	string si = ss.str();
	ss.clear();
	NaiveLongTailExperiment(inSyFilePrefix + "\\sythetic" + si, outFilePrefix + "\\sythetic_lt" + si, params);
	}
	*/


}

void NaiveLongTailExperiment(string infile, string outfile, const ParametersNaive &params)
{
	cout << "\t-----------" << infile << "------------" << endl;
	ifstream *in = new ifstream(infile.c_str());
	ofstream *out = new ofstream(outfile.c_str());

	/*
	unsigned long windowSize = 512;
	unsigned long motifLength = 64;
	unsigned long step = 16;
	double R = 20;
	*/

	size_t k = 20000;

	BFFindMotif *findMotif = new BFFindMotif(params.slideWindowSize, params.motifLength, params.R, params.step, in);
	clock_t start, end;
	start = clock();
	findMotif->FindMotif();
	end = clock();
	priority_queue<TimeSeriesCandidate> topk = findMotif->GetTopKMotifs(k);

	for(size_t i = 0; i < topk.size(); ++i)
	{
		TimeSeriesCandidate ts = topk.top();
		topk.pop();
		*out << ts.first << "\t";
	}
	out->close();

	cout << "There are totally " << findMotif->GetLastIndex() << " sample points." << endl;
	cout << "Find out " << findMotif->GetCandidateMotifSize() << " candidate motifs." << endl;
	cout << "Time cost:" << (end - start) / CLOCKS_PER_SEC << " seconds." << endl;

}

/**
*	SIMMA related experiments
*/
void SIMMALongTailExperiments()
{
	cout << "-------SIMMA long tail test----------" << endl;

	//	Test for real datasets
	string inFilePrefix = "G:\\dataset\\TimeSeries\\real";
	string outFilePrefix = "G:\\";
	ParametersSIMMA params;
	params.slideWindowSize = 512;
	params.motifLength = 64;
	params.R = 120;
	params.k = 10;
	params.sigma = 0.0002;


	//	IBM data
	SIMMALongTailExperiment(inFilePrefix + "\\IBM_data", outFilePrefix + "\\IBM_SIMMA_lt", params);
	//	Dodger
	SIMMALongTailExperiment(inFilePrefix + "\\Traffic_data", outFilePrefix + "\\Traffic_SIMMA_lt", params);
	//	ICU
	SIMMALongTailExperiment(inFilePrefix + "\\ICU_data", outFilePrefix + "\\ICU_SIMMA_lt", params);


	/*	
	//	synthetic
	string inSyFilePrefix = "G:\\dataset\\TimeSeries\\sythetic";
	for(size_t i = 9; i < 10; ++i)
	{
	stringstream ss;
	ss << i;
	string si = ss.str();
	ss.clear();
	SIMMALongTailExperiment(inSyFilePrefix + "\\sythetic" + si, outFilePrefix + "\\sythetic_lt" + si, params);
	}
	*/
}

void SIMMALongTailExperiment(std::string infile, std::string  outfile, const motif::ParametersSIMMA &params)
{
	cout << "\t-----------" << infile << "------------" << endl;
	ifstream *in = new ifstream(infile.c_str());
	ofstream *out = new ofstream(outfile.c_str());
	size_t k = 1000;

	SIMMA *simma = new SIMMA(params.slideWindowSize, params.motifLength, params.R, in, params.sigma, params.k);
	clock_t start, end;
	start = clock();
	simma->FindMotif();
	end = clock();
	priority_queue<TimeSeriesCandidate> topk = simma->GetTopKMotifs(k);

	for(size_t i = 0; i < topk.size(); ++i)
	{
		TimeSeriesCandidate ts = topk.top();
		topk.pop();
		*out << ts.first << "\t";
	}
	out->close();

	cout << "There are totally " << simma->GetLastIndex() << " sample points." << endl;
	cout << "Find out " << simma->GetCandidateMotifSize() << " candidate motifs." << endl;
	cout << "Time cost:" << (end - start) / CLOCKS_PER_SEC << " seconds." << endl;
}

void SIMMARvsMotifSize()
{
	cout << "---SIMMA efficiency test---" << endl;
	ParametersSIMMA paramsSIMMA;
	paramsSIMMA.slideWindowSize = 512;
	paramsSIMMA.motifLength = 64;
	paramsSIMMA.k = 10;
	paramsSIMMA.sigma = 0.0002;

	long k = 10000;

	double R = 5.0;
	string inFile = "G:\\dataset\\TimeSeries\\real\\IBM_data";
	string outLog = "G:\\log\\SIMMARlog";


	ofstream *out = new ofstream(outLog.c_str());

	tools::CThreadTime cThreadTime;
	cThreadTime.Begin();
	for(size_t i = 0; i < 200; ++i)
	{
		ifstream *in = new ifstream(inFile.c_str());
		paramsSIMMA.R = R + i * 1;
		cout << "\tRadius = " << paramsSIMMA.R << endl;
		SIMMA *simma = new SIMMA(paramsSIMMA.slideWindowSize, paramsSIMMA.motifLength, paramsSIMMA.R, in, paramsSIMMA.sigma, paramsSIMMA.k);
		simma->SetBufferCheck(false);
		simma->FindMotif();
		priority_queue<TimeSeriesCandidate> topk = simma->GetTopKMotifs(k);
		size_t size = topk.size();
		*out << paramsSIMMA.R << "\t" << size << endl;
		in->close();
	}
	*out << endl;
	out->close();

	__int64 elapse = cThreadTime.End();
	cout << "Time cost: " << elapse << endl;
}

void NaiveRvsMotifSize()
{
	cout << "---Naive efficiency test---" << endl;
	ParametersNaive params;
	params.slideWindowSize = 512;
	params.motifLength = 64;
	params.step = 16;

	long k = 10000;

	double R = 5.0;
	string inFile = "G:\\dataset\\TimeSeries\\real\\IBM_data";
	string outLog = "G:\\log\\NaiveRlog";


	ofstream *out = new ofstream(outLog.c_str());

	tools::CThreadTime cThreadTime;
	cThreadTime.Begin();
	for(size_t i = 0; i < 200; ++i)
	{
		ifstream *in = new ifstream(inFile.c_str());
		params.R = R + i * 1;
		cout << "\tRadius = " << params.R << endl;
		BFFindMotif *findMotif = new BFFindMotif(params.slideWindowSize, params.motifLength, params.R, params.step, in);
		findMotif->FindMotif();
		priority_queue<TimeSeriesCandidate> topk = findMotif->GetTopKMotifs(k);
		size_t size = topk.size();
		*out << params.R << "\t" << size << endl;
		in->close();
	}
	*out << endl;
	out->close();

	__int64 elapse = cThreadTime.End();
	cout << "Time cost: " << elapse << endl;
}

void SIMMASlideWindowSizevsMotifSize()
{
	cout << "---SIMMA efficiency test for slide window---" << endl;
	ParametersSIMMA paramsSIMMA;
	paramsSIMMA.motifLength = 64;
	paramsSIMMA.R = 120;
	paramsSIMMA.k = 10;
	paramsSIMMA.sigma = 0.0002;

	long k = 10000;

	unsigned long slideWindowSize = 128;
	string inFile = "G:\\dataset\\TimeSeries\\real\\IBM_data";
	string outLog = "G:\\log\\SIMMASlideWindowSizeLog";


	ofstream *out = new ofstream(outLog.c_str());

	tools::CThreadTime cThreadTime;
	cThreadTime.Begin();
	for(size_t i = 0; i < 50; ++i)
	{
		ifstream *in = new ifstream(inFile.c_str());
		paramsSIMMA.slideWindowSize = slideWindowSize + i * 32;
		cout << "\tSlide window size = " << paramsSIMMA.slideWindowSize << endl;
		SIMMA *simma = new SIMMA(paramsSIMMA.slideWindowSize, paramsSIMMA.motifLength, paramsSIMMA.R, in, paramsSIMMA.sigma, paramsSIMMA.k);
		simma->FindMotif();
		priority_queue<TimeSeriesCandidate> topk = simma->GetTopKMotifs(k);
		size_t size = topk.size();
		*out << paramsSIMMA.slideWindowSize << "\t" << size << endl;
		in->close();
	}
	*out << endl;
	out->close();

	__int64 elapse = cThreadTime.End();
	cout << "Time cost: " << elapse << endl;
}

void NaiveSlideWindowSizevsMotifSize()
{
	cout << "---Naive efficiency test for slide window---" << endl;
	ParametersNaive params;
	params.motifLength = 64;
	params.R = 120;
	params.step = 16;

	long k = 10000;

	unsigned long slideWindowSize = 128;
	string inFile = "G:\\dataset\\TimeSeries\\real\\IBM_data";
	string outLog = "G:\\log\\NaiveSlideWindowSizeLog";


	ofstream *out = new ofstream(outLog.c_str());

	tools::CThreadTime cThreadTime;
	cThreadTime.Begin();
	for(size_t i = 0; i < 50; ++i)
	{
		ifstream *in = new ifstream(inFile.c_str());
		params.slideWindowSize = slideWindowSize + i * 32;
		cout << "\tSlide window size = " << params.slideWindowSize << endl;
		BFFindMotif *findMotif = new BFFindMotif(params.slideWindowSize, params.motifLength, params.R, params.step, in);
		findMotif->FindMotif();
		priority_queue<TimeSeriesCandidate> topk = findMotif->GetTopKMotifs(k);
		size_t size = topk.size();
		*out << params.slideWindowSize << "\t" << size << endl;
		in->close();
	}
	*out << endl;
	out->close();

	__int64 elapse = cThreadTime.End();
	cout << "Time cost: " << elapse << endl;
}

void SIMMAEfficiencyMotifLength()
{
	cout <<"---SIMMA efficiency test for motif length---" << endl;

	ParametersSIMMA paramsSIMMA;
	paramsSIMMA.slideWindowSize = 4096;
	paramsSIMMA.motifLength = 32;
	paramsSIMMA.R = 120;
	paramsSIMMA.k = 10;
	paramsSIMMA.sigma = 0.0002;

	long k = 10000;

	unsigned long motifLength = 32;
	unsigned long motifStep = 4;
	string inFile = "G:\\dataset\\TimeSeries\\real\\IBM_data";
	string outLog = "G:\\log\\SIMMAMotifLengthLog";

	ofstream *out = new ofstream(outLog.c_str());

	for(size_t i = 0; i < 500; ++i)
	{
		ifstream *in = new ifstream(inFile.c_str());
		paramsSIMMA.motifLength = motifLength + i * motifStep;
		cout << "\tMotif length = " << paramsSIMMA.motifLength << endl;

		clock_t start, end;
		start = clock();
		SIMMA *simma = new SIMMA(paramsSIMMA.slideWindowSize, paramsSIMMA.motifLength, paramsSIMMA.R, in, paramsSIMMA.sigma, paramsSIMMA.k);
		simma->FindMotif();
		priority_queue<TimeSeriesCandidate> topk = simma->GetTopKMotifs(k);
		size_t size = topk.size();
		end = clock();
		*out << paramsSIMMA.motifLength << "\t" << (end - start) / CLOCKS_PER_SEC << endl;
		in->close();
	}
	*out << endl;
	out->close();

}

void NaiveEfficiencyMotifLength()
{
	cout <<"---Naive efficiency test for motif length---" << endl;

	ParametersNaive params;
	params.slideWindowSize = 4096;
	params.motifLength = 32;
	params.R = 120;
	params.step = 16;

	long k = 10000;

	unsigned long motifLength = 32;
	unsigned long motifStep = 4;
	string inFile = "G:\\dataset\\TimeSeries\\real\\IBM_data";
	string outLog = "G:\\log\\NaiveMotifLengthLog";

	ofstream *out = new ofstream(outLog.c_str());

	for(size_t i = 0; i < 500; ++i)
	{
		ifstream *in = new ifstream(inFile.c_str());
		params.motifLength = motifLength + i * 16;
		cout << "\tMotif length = " << params.motifLength << endl;

		clock_t start, end;
		start = clock();
		BFFindMotif *findMotif = new BFFindMotif(params.slideWindowSize, params.motifLength, params.R, params.step, in);
		findMotif->FindMotif();
		priority_queue<TimeSeriesCandidate> topk = findMotif->GetTopKMotifs(k);
		size_t size = topk.size();
		end = clock();
		*out << params.motifLength << "\t" << (end - start) / CLOCKS_PER_SEC << endl;
		in->close();
	}
	*out << endl;
	out->close();
}

void SIMMAMotifNumber()
{

	string inFilePrefix = "G:\\dataset\\TimeSeries\\sythetic\\sythetic";
	string outLogMotifNumber = "G:\\log\\SIMMAMotifNumber";
	ofstream *outMotif = new ofstream(outLogMotifNumber.c_str());

	//	SIMMA motif number
	//	SIMMA efficiency
	cout << "--SIMMA motif number--" << endl;
	ParametersSIMMA paramsSIMMA;
	paramsSIMMA.slideWindowSize = 512;
	paramsSIMMA.motifLength = 64;
	paramsSIMMA.R = 120;
	paramsSIMMA.k = 10;
	paramsSIMMA.sigma = 0.0002;

	for(size_t i = 0; i < 10; ++i)
	{
		if(i > 3 && i < 8)
			continue;

		stringstream ss;
		ss << (i + 1);
		string inFile = inFilePrefix + ss.str();
		ifstream *in = new ifstream(inFile.c_str());
		cout << "\nProcess file:" << inFile << endl;

		size_t k = 1000;

		SIMMA *simma = new SIMMA(paramsSIMMA.slideWindowSize, paramsSIMMA.motifLength, paramsSIMMA.R, in, paramsSIMMA.sigma, paramsSIMMA.k);
		simma->SetBufferCheck(false);
		simma->FindMotif();
		priority_queue<TimeSeriesCandidate> topk = simma->GetTopKMotifs(k);
		size_t size = topk.size();

		*outMotif << size << "\t";	
		in->close();
	}
	outMotif->close();
}

void NaiveMotifNumber()
{
	string inFilePrefix = "G:\\dataset\\TimeSeries\\sythetic\\sythetic";
	string outLogMotifNumber = "G:\\log\\NaiveMotifNumber";
	ofstream *outMotif = new ofstream(outLogMotifNumber.c_str());

	cout << "--SIMMA motif number--" << endl;
	ParametersNaive paramsNaive;
	paramsNaive.slideWindowSize = 512;
	paramsNaive.motifLength = 64;
	paramsNaive.R = 120;
	paramsNaive.step = 16;

	for(size_t i = 0; i < 10; ++i)
	{
		if(i > 3 && i < 8)
			continue;

		stringstream ss;
		ss << (i + 1);
		string inFile = inFilePrefix + ss.str();
		ifstream *in = new ifstream(inFile.c_str());
		cout << " \nProcessed file:" << inFile << endl;

		size_t k = 1000000;
		BFFindMotif *findMotif = new BFFindMotif(paramsNaive.slideWindowSize, paramsNaive.motifLength, paramsNaive.R, paramsNaive.step, in);
		findMotif->FindMotif();

		priority_queue<TimeSeriesCandidate> topk = findMotif->GetTopKMotifs(k);
		size_t size = topk.size();

		*outMotif << size << endl;
		in->close();
	}

	*outMotif << endl;
	outMotif->close();
}


void SIMMAEfficiency()
{

	string outLogEfficiency = "G:\\log\\simmaEfficiency";
	ofstream *out = new ofstream(outLogEfficiency.c_str(), ios::app);
	string inFilePrefix = "G:\\dataset\\TimeSeries\\sythetic\\sythetic";

	//	SIMMA efficiency
	cout << "--SIMMA Efficiency--" << endl;
	ParametersSIMMA paramsSIMMA;
	paramsSIMMA.slideWindowSize = 512;
	paramsSIMMA.motifLength = 64;
	paramsSIMMA.R = 50;
	paramsSIMMA.k = 10;
	paramsSIMMA.sigma = 0.005;

	for(size_t i = 0; i < 10; ++i)
	{
		if(i > 3 && i < 8)
			continue;

		stringstream ss;
		ss << (i + 1);
		string inFile = inFilePrefix + ss.str();
		ifstream *in = new ifstream(inFile.c_str());
		cout << "\nProcess file:" << inFile << endl;

		size_t k = 1000;

		SIMMA *simma = new SIMMA(paramsSIMMA.slideWindowSize, paramsSIMMA.motifLength, paramsSIMMA.R, in, paramsSIMMA.sigma, paramsSIMMA.k);
		tools::CThreadTime cThreadTime;
		cThreadTime.Begin();
		simma->FindMotif();
		long elapse = cThreadTime.End();

		*out << elapse << "\t";	
		in->close();
	}

	*out << endl;
	out->close();
}

void NaiveEfficiency()
{
	string outLogEfficiency = "G:\\log\\naiveEfficiency";
	ofstream *out = new ofstream(outLogEfficiency.c_str(), ios::app);
	string inFilePrefix = "G:\\dataset\\TimeSeries\\sythetic\\sythetic";

	//	naive efficiency
	cout << "--naive Efficiency--" << endl; 
	ParametersNaive paramsNaive;
	paramsNaive.slideWindowSize = 512;
	paramsNaive.motifLength = 64;
	paramsNaive.R = 50;
	paramsNaive.step = 16;

	for(size_t i = 0; i < 10; ++i)
	{
		if(i > 3 && i < 8)
			continue;

		stringstream ss;
		ss << (i + 1);
		string inFile = inFilePrefix + ss.str();
		ifstream *in = new ifstream(inFile.c_str());
		cout << " \nProcess file:" << inFile << endl;

		size_t k = 1000000;
		tools::CThreadTime cThreadTime;
		BFFindMotif *findMotif = new BFFindMotif(paramsNaive.slideWindowSize, paramsNaive.motifLength, paramsNaive.R, paramsNaive.step, in);
		cThreadTime.Begin();
		findMotif->FindMotif();
		long lasts = cThreadTime.End();

		*out << lasts << "\t";	
		in->close();
	}

	*out << endl;
	out->close();
}

void SIMMAREfficiency()
{

	cout << "---SIMMA efficiency test---" << endl;
	ParametersSIMMA paramsSIMMA;
	paramsSIMMA.slideWindowSize = 512;
	paramsSIMMA.motifLength = 64;
	paramsSIMMA.k = 10;
	paramsSIMMA.sigma = 0.005;

	long k = 10000;

	double R = 10.0;
	string inFile = "G:\\dataset\\TimeSeries\\sythetic\\sythetic1";
	string outLog = "G:\\log\\SIMMAREfficiencylog";


	ofstream *out = new ofstream(outLog.c_str());

	tools::CThreadTime cThreadTime;
	cThreadTime.Begin();
	for(size_t i = 0; i < 100; ++i)
	{
		ifstream *in = new ifstream(inFile.c_str());
		paramsSIMMA.R = R + i * 1;
		cout << "\tRadius = " << paramsSIMMA.R << ", time costs:";
		SIMMA *simma = new SIMMA(paramsSIMMA.slideWindowSize, paramsSIMMA.motifLength, paramsSIMMA.R, in, paramsSIMMA.sigma, paramsSIMMA.k);
		simma->SetBufferCheck(false);
		clock_t start, end;
		start = clock();
		simma->FindMotif();
		priority_queue<TimeSeriesCandidate> topk = simma->GetTopKMotifs(k);
		size_t size = topk.size();
		end = clock();
		long lasts = (end - start) / CLOCKS_PER_SEC; 
		cout << lasts << " seconds." << endl;
		*out << paramsSIMMA.R << "\t" << lasts << endl;
		in->close();
	}
	*out << endl;
	out->close();

	__int64 elapse = cThreadTime.End();
	cout << "Time cost: " << elapse << endl;
}

void NaiveREfficiency()
{

	cout << "---Naive R vs Efficiency test---" << endl;
	ParametersNaive params;
	params.slideWindowSize = 512;
	params.motifLength = 64;
	params.step = 16;

	long k = 10000;

	double R = 10.0;
	string inFile = "G:\\dataset\\TimeSeries\\sythetic\\sythetic1";
	string outLog = "G:\\log\\NaiveREfficiencylog";


	ofstream *out = new ofstream(outLog.c_str());

	for(size_t i = 0; i < 100; ++i)
	{
		ifstream *in = new ifstream(inFile.c_str());
		params.R = R + i * 1;
		cout << "\tRadius = " << params.R << ", time costs:";

		BFFindMotif *findMotif = new BFFindMotif(params.slideWindowSize, params.motifLength, params.R, params.step, in);
		clock_t start, end;
		start = clock();
		findMotif->FindMotif();
		priority_queue<TimeSeriesCandidate> topk = findMotif->GetTopKMotifs(k);
		size_t size = topk.size();
		end = clock();
		long lasts = (end - start) / CLOCKS_PER_SEC;
		cout << lasts << " seconds." << endl;
		*out << params.R << "\t" << lasts << endl;
		in->close();
	}
	*out << endl;
	out->close();

}

void SIMMAShrinkTest()
{

	string inFilePrefix = "G:\\dataset\\TimeSeries\\sythetic\\sythetic";
	string outLogMotifNumber = "G:\\log\\SIMMShrinkTest";
	ofstream *outMotif = new ofstream(outLogMotifNumber.c_str());

	//	SIMMA motif number
	//	SIMMA efficiency
	cout << "--SIMMA motif number--" << endl;
	ParametersSIMMA paramsSIMMA;
	paramsSIMMA.slideWindowSize = 512;
	paramsSIMMA.motifLength = 64;
	paramsSIMMA.R = 30;
	paramsSIMMA.k = 10;
	paramsSIMMA.sigma = 0.005;

	for(size_t i = 0; i < 10; ++i)
	{
		if(i > 3 && i < 8)
			continue;

		stringstream ss;
		ss << (i + 1);
		string inFile = inFilePrefix + ss.str();
		ifstream *in = new ifstream(inFile.c_str());
		cout << "\nProcess file:" << inFile << endl;

		size_t k = 1000;

		SIMMA *simma = new SIMMA(paramsSIMMA.slideWindowSize, paramsSIMMA.motifLength, paramsSIMMA.R, in, paramsSIMMA.sigma, paramsSIMMA.k);
		simma->FindMotif();
		size_t shrinkCount = simma->GetShrinkCount();
		*outMotif << shrinkCount << endl;	
		in->close();
	}
	outMotif->close();
}