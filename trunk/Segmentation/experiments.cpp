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
	double sharpness[] = {0.001, 0.001, 0.001, 0.001, 0.001, 0.0005, 0.0002, 0.0001, 0.00005, 0.00001, 0.00001};


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
			<< amount[i] << endl << "Professing time:" << lasts << endl << endl;
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
	double threshold = 0.3;

	for(size_t i = 1; i <=1; ++i)
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
			<< "Approximate Size:" << swbu->GetApprSize() << endl << "Professing time:" << lasts << endl << endl;
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
		mdlsw->SetBatchSize(20);
		mdlsw->Approximate();
		end = clock();
		double lasts = (end - start) / CLOCKS_PER_SEC;
		cout << "End of processing dataset" << i << ", time elapse:" << lasts << endl << endl;
		outlog << "dataset" << i << endl << "Time Series Size:" << mdlsw->GetSize() << endl
			<< "Approximate Size:" << mdlsw->GetApprSize() << endl << "Professing time:" << lasts << endl << endl;
	}
	outlog.close();
}

void ExperimentPrecision()
{

}