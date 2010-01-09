#include "experiments.h"
#include "utility.h"

using namespace std;
using namespace utility;
using namespace segmenting;


void Experiment()
{

}

//	Generate 11 datasets
void ExperimentGenerateData(string pathPrefix)
{
	cout << "Begin to generate data." << endl;

	string outfile = pathPrefix + "\\DataGenerateLog.txt";
	ofstream outlog(outfile.c_str());

	long long amount[] = {10 * 1000, 100 * 1000, 1000 * 1000, 
							10 * 1000 * 1000, 10 * 1000 * 1000, 10 * 1000 * 1000, 10 * 1000 * 1000, 10 * 1000 * 1000, 
							100 * 1000 * 1000,			//	100M
							1000 * 1000 * 1000,			//	1G
							(long long)10 * 1000 * 1000 * 1000};	//	10G

	double min = 0;
	double max[] = {3000, 3000, 3000, 3000, 3000, 3000, 3000, 3000, 3000, 5000, 5000};
	double sharpness[] = {0.001, 0.001, 0.001, 0.001, 0.002, 0.0005, 0.0002, 0.0001, 0.00005, 0.00001, 0.00001};


	//	generate dataset with 10k, 100k, 1M, 10M, 100M, 1G,  10G
	for(size_t i = 10; i < 11; ++i)
	{
		clock_t start, end;
		start = clock();
		cout << "Generate dataset" << i << " with size " << amount[i] << endl;
		stringstream ss;
		ss << (i + 1);
		string suffix;
		ss >> suffix;

		string path = pathPrefix + "\\sythetic" + suffix;
		RandomDataGenerator(amount[i], path, min, max[i]);

		end = clock();
		double lasts = (end -start) / CLOCKS_PER_SEC;
		outlog << "dataset" << i << endl << "Time Series Size:" 
			<< amount[i] << endl << "Processing time:" << lasts << endl << endl;
	}

	cout << "End of generating data." << endl;
}

void ExperimentEfficient()
{
	string pathPrefix = "G:\\dataset\\TimeSeries\\sythetic";
	ExperimentEfficientGlobal(pathPrefix);
	ExperimentEfficientSlideWindow(pathPrefix);
	ExperimentEfficientMDLSlideWindow(pathPrefix);
}

void ExperimentEfficientGlobal(std::string pathPrefix)
{

}

void ExperimentEfficientSlideWindow(std::string pathPrefix)
{
	string logFile = pathPrefix + "\\SlideWindowEfficientLog.txt";
	ofstream outlog(logFile.c_str());

	int windowSize = 500;
	double threshold = 1.1;

	for(size_t i = 1; i <=8; ++i)
	{
		clock_t start, end;
		start = clock();

		stringstream ss;
		ss << i;
		string suffix;
		ss >> suffix;
		string infile = pathPrefix + "\\sythetic" + suffix;
		string outfile = pathPrefix + "\\result\\sw" + suffix;

		cout << "Begining of processing dataset" << i << endl;
		cout << "Input filepath:" << infile << endl;
		cout << "Output filepath:" << outfile << endl;

		ifstream *in = new ifstream(infile.c_str());
		ofstream *out = new ofstream(outfile.c_str());
		SlideWindowBottomUp *swbu = new SlideWindowBottomUp(threshold, windowSize, in, out);
		swbu->SetPartitionsPerTime(5);
		swbu->Approximate();
		end = clock();
		double lasts = (end - start) / CLOCKS_PER_SEC;
		cout << "End of processing dataset" << i << ", time elapse:" << lasts << endl << endl;
		outlog << "dataset" << i << endl << "Time Series Size:" << swbu->GetSize() << endl
			<< "Approximate Size:" << swbu->GetApprSize() << endl << "Processing time:" << lasts << endl << endl;
	}

	outlog.close();
}

void ExperimentEfficientMDLSlideWindow(std::string pathPrefix)
{
	string logFile = pathPrefix + "\\MDLSlideWindowEfficientLog.txt";
	ofstream outlog(logFile.c_str());


	int windowSize = 50;

	for(size_t i = 3; i <=3; ++i)
	{
		clock_t start, end;
		start = clock();

		stringstream ss;
		ss << i;
		string suffix;
		ss >> suffix;
		string infile = pathPrefix + "\\sythetic" + suffix;
		string outfile = pathPrefix + "\\result\\mdl" + suffix;

		cout << "Begining of processing dataset" << i << endl;
		cout << "Input filepath:" << infile << endl;
		cout << "Output filepath:" << outfile << endl;

		ifstream *in = new ifstream(infile.c_str());
		ofstream *out = new ofstream(outfile.c_str());
		MDLSlideWindow *mdlsw = new MDLSlideWindow(windowSize, in, out);
		//mdlsw->SetBatchSize(20);
		mdlsw->Approximate();
		end = clock();
		double lasts = (end - start) / CLOCKS_PER_SEC;
		cout << "End of processing dataset" << i << ", time elapse:" << lasts << endl << endl;
		outlog << "dataset" << i << endl << "Time Series Size:" << mdlsw->GetSize() << endl
			<< "Approximate Size:" << mdlsw->GetApprSize() << endl << "Processing time:" << lasts << endl << endl;
	}
	outlog.close();
}

void ExperimentPrecision(const string &pathPrefix)
{
	string IBM = pathPrefix + "\\IBM_data";
	string ICU = pathPrefix + "\\ICU_data";
	string Traffic = pathPrefix + "\\Traffic_data";

	string IBMoutput_bu = IBM + "_result_bu";
	string ICUoutput_bu = ICU + "_result_bu";
	string Trafficoutput_bu = Traffic + "_result_bu";

	string IBMoutput_mdl = IBM + "_result_mdl";
	string ICUoutput_mdl = ICU + "_result_mdl";
	string Trafficoutput_mdl = Traffic + "_result_mdl";

	string seglog = pathPrefix + "\\SegmentingLog.txt";
	string measurelog = pathPrefix + "\\MeasureLog.txt";

	/*	Segmenting	*/
	double threshold1 = 1.4;
	ExperimentBottomUp(IBM, IBMoutput_bu, seglog, threshold1);
	ExperimentMDLSlideWindow(IBM, IBMoutput_mdl, seglog);

	double threshold2 = 3.55;
	ExperimentBottomUp(ICU, ICUoutput_bu, seglog, threshold2);
	ExperimentMDLSlideWindow(ICU, ICUoutput_mdl, seglog);

	double threshold3 = 1.005;
	ExperimentBottomUp(Traffic, Trafficoutput_bu, seglog, threshold3);
	ExperimentMDLSlideWindow(Traffic, Trafficoutput_mdl, seglog);
	

	/*	Measure	*/
	/*
	ExperimentMeasure(IBM, IBMoutput_bu, measurelog);
	ExperimentMeasure(IBM, IBMoutput_mdl, measurelog);

	ExperimentMeasure(ICU, ICUoutput_bu, measurelog);
	ExperimentMeasure(ICU, ICUoutput_mdl, measurelog);

	ExperimentMeasure(Traffic, Trafficoutput_bu, measurelog);
	ExperimentMeasure(Traffic, Trafficoutput_mdl, measurelog);
	*/
}

void ExperimentBottomUp(const string &filename, const string &resultFilename, const string &log, double threshold)
{
	cout << "Begining of processing dataset." << endl;
	cout << "Input filepath:" << filename << endl;
	cout << "Output filepath:" << resultFilename << endl;

	ofstream outlog(log.c_str(), ios::app);

	clock_t start, end;
	start = clock();

	vector<double> vec;
	bool success = ReadDataFromFile(vec, filename);
	
	ContinuousBottomUp *bottomUp = new ContinuousBottomUp(threshold, vec);
	bottomUp->SetOutputFilepath(resultFilename);
	bottomUp->Approximate();

	end = clock();
	double lasts = (end - start) / CLOCKS_PER_SEC;

	cout << "End of processing dataset, time elapse:" << lasts << endl << endl;
	outlog << "Time Series Size:" << vec.size() << endl
		<< "Approximate Size:" << bottomUp->GetSegmentSize() << endl << "Processing time:" << lasts << endl << endl;
	
}

void ExperimentMDLSlideWindow(const string &filename, const string &resultFilename, const string &log)
{
	cout << "Begining of processing dataset." << endl;
	cout << "Input filepath:" << filename << endl;
	cout << "Output filepath:" << resultFilename << endl;

	ifstream *in = new ifstream(filename.c_str());
	ofstream *out = new ofstream(resultFilename.c_str());
	ofstream outlog(log.c_str(), ios::app);

	clock_t start, end;
	start = clock();

	size_t windowSize = 500;
	MDLSlideWindow *mdlsw = new MDLSlideWindow(windowSize, in, out);
	mdlsw->Approximate();

	end = clock();
	double lasts = (end - start) / CLOCKS_PER_SEC;

	cout << "End of processing dataset, time elapse:" << lasts << endl << endl;
	outlog << "Time Series Size:" << mdlsw->GetSize() << endl
		<< "Approximate Size:" << mdlsw->GetApprSize() << endl << "Processing time:" << lasts << endl << endl;

}

void ExperimentMeasure(const string &filename, const string &resultFilename, const string &log)
{
	cout << "Begining of measureing dataset." << endl;
	cout << "Input filepath:" << filename << endl;
	cout << "Output filepath:" << resultFilename << endl;

	ifstream inResult(resultFilename.c_str());
	ifstream inOrg(filename.c_str());
	ofstream outlog(log.c_str(), ios::app);

	clock_t start, end;
	start = clock();

	istringstream iss;

	double error_rate = 0.0;
	long long segcount = 0;

	while(false == inResult.eof())
	{
		string line;
		getline(inResult, line, '\n');

		size_t pos = line.find_first_of(" ");
		string strIndexResult = line.substr(0, pos);
		string strValueResult = line.substr(pos, line.size() - pos);

		long long indexResult;
		double valueResult;

		iss.str(strIndexResult);
		iss >> indexResult;
		iss.clear();
		iss.str(strValueResult);
		iss >> valueResult;
		iss.clear();

		long long indexOrg;
		double valueOrg;


		long orgcount = 0;
		double seg_error = 0.0;

		do
		{
			string lineOrg;
			getline(inOrg, lineOrg, '\n');

			size_t pos = lineOrg.find_first_of(" ");
			string strIndexOrg = line.substr(0, pos);
			string strValueOrg = line.substr(pos, line.size() - pos);

			long long indexOrg;
			double valueOrg;

			iss.str(strIndexOrg);
			iss >> indexOrg;
			iss.clear();
			iss.str(strValueOrg);
			iss >> valueOrg;
			iss.clear();

			++orgcount;

			

		}while(indexOrg < indexResult || inOrg.eof());

	}

	inResult.close();
	inOrg.close();


	end = clock();
	double lasts = (end - start) / CLOCKS_PER_SEC;

	cout << "End of processing dataset, time elapse:" << lasts << endl << endl;

}
