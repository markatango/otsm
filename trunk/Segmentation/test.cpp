#include "utility.h"
#include "segmenting.h"
#include "experiments.h"

using namespace std;
using namespace utility;
using namespace segmenting;

void GenerateData(int argc, char **argv)
{
	cout << "Begin to generate data." << endl;

	double min = 0;
	double max = 1000;

	string path = "G:\\dataset\\Time Series\\sythetic_10k";

	long long amount = 2500 * 1000;	// unit gigabytes
	amount *= 1000;	//	totally 2.5 gigabytes

	if(2 == argc)
	{
		istringstream iss;
		iss.str(argv[1]);
		iss >> amount;
		cout << "Amount set to " << amount << endl;
	}

	if(4 <= argc)
	{
		istringstream iss;
		iss.str(argv[2]);
		iss >> min;
		iss.clear();

		iss.str(argv[3]);
		iss >> max;
		iss.clear();

	}

	clock_t start, end;

	start = clock();

	Generator(min, max, amount, path);

	end = clock();

	long long last = (end - start) / CLOCKS_PER_SEC;

	cout << "Time cost: " << last << " seconds." << endl;
}

void SineDataGenerator(int argc, char **argv)
{
	cout << "Begin to generate data." << endl;
	cout << "Param size:" << argc << endl;

	string path = "G:\\dataset\\Time Series\\sythetic_10k";

	double a = 100.0;
	double w = 1.0;
	double b = 0.0;

	long long amount = 10000;	// unit gigabytes

	if(2 <= argc)
	{
		istringstream iss;
		iss.str(argv[1]);
		iss >> amount;
		cout << "Amount set to: " << amount << endl;
	}


	if(3 <= argc)
	{
		istringstream iss;
		iss.str(argv[2]);
		iss >> a;
		cout << "a set to:" << a << endl;
	}

	if(4 <= argc)
	{
		istringstream iss;
		iss.str(argv[3]);
		iss >> w;
		cout << "w set to:" << w << endl;
	}

	if(5 <= argc)
	{
		istringstream iss;
		iss.str(argv[4]);
		iss >> b;
		cout << "b set to:" << b << endl;
	}


	clock_t start, end;

	start = clock();

	SineGenerator(amount, path, a, w, b);

	end = clock();

	long long last = (end - start) / CLOCKS_PER_SEC;

	cout << "Time cost: " << last << " seconds." << endl;
}

void RandomGenerate(int argc, char **argv)
{
	cout << "Begin to generate data." << endl;

	long long amount = 10000;
	string path = "G:\\dataset\\Time Series\\sythetic";
	double min = 100;
	double max = 2000;
	double sharpness = 0.001;

	if(2 <= argc)
	{
		istringstream iss;
		iss.str(argv[1]);
		iss >> amount;
		cout << "Datasize:" << amount << endl;
	}

	if(3 <= argc)
	{
		path = argv[2];
		cout << "Filepath:" << path << endl;
	}

	if(4 <= argc)
	{
		istringstream iss;
		iss.str(argv[3]);
		iss >> min;
		cout << "Min:" << min << endl;
	}

	if(5 <= argc)
	{
		istringstream iss;
		iss.str(argv[4]);
		iss >> max;
		cout << "Max:" << max << endl;
	}

	if(6 <= argc)
	{
		istringstream iss;
		iss.str(argv[5]);
		iss >> sharpness;
		cout << "Sharpness:" << sharpness << endl;
	}

	clock_t start, end;
	start = clock();

	RandomDataGenerator(amount, path, min, max);

	end = clock();

	double lasts = (end - start) / CLOCKS_PER_SEC;
	cout << "Time cost:" << lasts << endl;
}

void main(int argc, char **argv)
{
	vector<double> vec;
	//string filename = "G:\\dataset\\TimeSeries\\MHSETS\\ECOLOGY\\ARCTIC.1";
	//string filename = "G:\\dataset\\TimeSeries\\MHSETS\\ASKEW\\ASKEW23.1";
	//string filename = "G:\\dataset\\TimeSeries\\MHSETS\\ASTATKIE\\FISHERT.1";
	string filename = "G:\\dataset\\TimeSeries\\sythetic\\sythetic1";
	//string outfilename = "G:\\dataset\\TimeSeries\\result\\sythetic";
	string outfilename;// = "g:\\dataset\\TimeSeries\\result\\sythetic";

	if(2 <= argc)
	{
		filename = argv[1];
	}

	size_t pos = filename.find_last_of("\\");
	outfilename = "g:\\dataset\\TimeSeries\\sythetic\\result\\sw" + filename.substr(pos + 1, filename.size() - pos);


	bool success = ReadDataFromFile(vec, filename);
	
	/*
	for(size_t i = 0; i < vec.size(); ++i)
	{
		cout << vec[i] << endl;
	}
	*/
	

	double threshold = 2.2;
	//double threshold = 13.2;

	/*
	ParameterMemorySegmenting *bottomUp1 = new BottomUp(threshold, vec);
	bottomUp1->Approximate();

	ParameterMemorySegmenting *bottomUp2 = new ContinuousBottomUp(threshold, vec);
	bottomUp2->Approximate();
	*/	

	size_t windowSize = 500;
	ifstream *in = new ifstream(filename.c_str());
	ofstream *out = new ofstream(outfilename.c_str());

	/*
	SlideWindowBottomUp *swbp = new SlideWindowBottomUp(threshold, windowSize, in, out);
	swbp->SetBatchSize(3);
	swbp->SetPartitionsPerTime(1);
	swbp->Approximate();
	*/

	clock_t start, end;
	start = clock();

	//MDLSlideWindow *mdlsw = new MDLSlideWindow(windowSize, in, out);
	//mdlsw->SetBatchSize(3);
	//mdlsw->Approximate();

	//SineDataGenerator(argc, argv);
	//RandomGenerate(argc, argv);

	string pathPrefix = "G:\\dataset\\TimeSeries\\sythetic";
	//ExperimentGenerateData(pathPrefix);
	//ExperimentEfficientSlideWindow(pathPrefix);
	//ExperimentEfficientMDLSlideWindow(pathPrefix);

	string realDataPathPrefix = "G:\\dataset\\TimeSeries\\real";
	//ExperimentPrecision(realDataPathPrefix);
	
	ExperimentThresholdTest();

	end = clock();
	long long lasts = (end - start) / CLOCKS_PER_SEC;
	cout << "Cost time:" << lasts << " seconds." << endl;

	cout << endl << endl;

}