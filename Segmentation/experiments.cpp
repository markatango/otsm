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
	//ExperimentEfficientGlobal(pathPrefix);
	ExperimentEfficientSlideWindow(pathPrefix);
	//ExperimentEfficientMDLSlideWindow(pathPrefix);
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

	for(size_t i = 9; i <=9; ++i)
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
	//	Precision of real data

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
	/*
	double threshold1 = 0.91;
	ExperimentBottomUp(IBM, IBMoutput_bu, seglog, threshold1);
	ExperimentMDLSlideWindow(IBM, IBMoutput_mdl, seglog);

	
	double threshold2 = 3.0;
	ExperimentBottomUp(ICU, ICUoutput_bu, seglog, threshold2);
	ExperimentMDLSlideWindow(ICU, ICUoutput_mdl, seglog);

	double threshold3 = 1.005;
	ExperimentBottomUp(Traffic, Trafficoutput_bu, seglog, threshold3);
	ExperimentMDLSlideWindow(Traffic, Trafficoutput_mdl, seglog);
	*/

	/*	Measure	*/
	/*
	ExperimentMeasure(IBM, IBMoutput_bu, measurelog);
	
	ExperimentMeasure(IBM, IBMoutput_mdl, measurelog);

	ExperimentMeasure(ICU, ICUoutput_bu, measurelog);
	ExperimentMeasure(ICU, ICUoutput_mdl, measurelog);
	
	
	ExperimentMeasure(Traffic, Trafficoutput_bu, measurelog);
	ExperimentMeasure(Traffic, Trafficoutput_mdl, measurelog);
	*/
	
	//	Precision of sythetic data
	char buffer[10];
	for(size_t i = 8; i < 10; ++i)
	{
		string orgFilepath = "G:\\dataset\\TimeSeries\\sythetic\\";
		orgFilepath += "\\sythetic";
		orgFilepath += _itoa_s(i, buffer, 10);;
		string swFilepath = "G:\\dataset\\TimeSeries\\sythetic\\result\\";
		swFilepath += "\\sw";
		swFilepath += _itoa_s(i, buffer, 10);;
		string mdlFilepath = "G:\\dataset\\TimeSeries\\sythetic\\result\\";
		mdlFilepath += "\\mdl";
		mdlFilepath += _itoa_s(i, buffer, 10);

		ExperimentMeasure(orgFilepath, swFilepath, measurelog);
		ExperimentMeasure(orgFilepath, mdlFilepath, measurelog);
	}
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

	size_t windowSize = 50;
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
	cout << "Result filepath:" << resultFilename << endl;

	ifstream inResult(resultFilename.c_str());
	ifstream inOrg(filename.c_str());
	ofstream outlog(log.c_str(), ios::app);

	clock_t start, end;
	start = clock();

	istringstream iss;

	long double error_rate = 0.0;
	long long segcount = 0;		//	number of segments

	Point segPoint1;
	Point segPoint2;

	//	Read first point of segment
	string line;
	getline(inResult, line, '\n');
	++segcount;

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
	segPoint1 = std::make_pair(indexResult, valueResult);

	long long org = 0;		//	total number of original points

	while(false == inResult.eof())
	{
		++segcount;

		string line;
		getline(inResult, line, '\n');
		
		if(line == "")	//	eof
			break;
		
		long long indexResult;
		double valueResult;
		size_t posSeg = line.find_first_of(" ");
		string strIndexResult = line.substr(0, posSeg);
		
		//printf("%d, Total size %d, From %d with size %d.\n", ++seg, line.size(), posSeg, line.size() - posSeg - 1);
		string strValueResult = line.substr(posSeg, line.size() - posSeg);
		
		iss.str(strIndexResult);
		iss >> indexResult;
		iss.clear();
		iss.str(strValueResult);
		iss >> valueResult;
		iss.clear();

		segPoint2 = std::make_pair(indexResult, valueResult);
		//cout << "[(" << segPoint1.first << "," << segPoint1.second << "),(" << segPoint2.first << "," << segPoint2.second << ")]" << endl;
		
		long orgcount = 0;		//	original point number in segment
		double seg_error = 0.0;	//	error rate in segment
		
		
		//	Read first point of original point
		long long indexOrg;
		double valueOrg;

		Point orgPoint1;
		Point orgPoint2;


		string lineOrg;
		getline(inOrg, lineOrg, '\n');

		++orgcount;

		//cout << "\tin seg[";
		//cout << "(" << org << "," << lineOrg << "),";

		indexOrg = org;
		iss.str(lineOrg);
		iss >> valueOrg;
		iss.clear();

		orgPoint1 = std::make_pair(indexOrg, valueOrg);
		
		do
		{
			++orgcount;

			string lineOrg;
			getline(inOrg, lineOrg, '\n');

			if(lineOrg == "")
				break;

			indexOrg = ++org;			
			iss.str(lineOrg);
			iss >> valueOrg;
			iss.clear();
			//cout << "(" << indexOrg << "," << valueOrg << "),";

			orgPoint2 = std::make_pair(indexOrg, valueOrg);

			//	Measure error
			if((segPoint1.first == orgPoint1.first) && (segPoint2.first == orgPoint2.first))
				seg_error += 0;
			else
			{
				double v = Line2LineDistance(orgPoint1, orgPoint2, segPoint1, segPoint2);
				//cout << endl << seg_error << "+=" << v << endl;
				seg_error += v;
			}
			//cout << "segerror:" << seg_error << endl;

			orgPoint1 = orgPoint2;
		
		}while(indexOrg < indexResult || inOrg.eof());
		//cout << "]" << endl;
		

		//error_rate += seg_error / orgcount;

		//cout << "Orgcount:" << orgcount << endl;
		
		if(0 != orgcount)
			error_rate += (seg_error / orgcount);

		segPoint1 = segPoint2;
		
		//system("pause");
	}
	

	//error_rate = error_rate / segcount;

	cout << "Segcount:" << segcount << endl;
	cout << (error_rate /= segcount) << endl;


	end = clock();
	double lasts = (end - start) / CLOCKS_PER_SEC;

	outlog << "Original data file:" << filename.c_str() << endl
		<< "Result data file:" << resultFilename.c_str() << endl
		<< "Error rate:" << error_rate << endl;
	cout << "End of processing dataset, time elapse:" << lasts << endl << endl;

	inResult.close();
	inOrg.close();
	outlog.close();

}


void ExperimentThresholdTest()
{
	const size_t sizefile = 3;
	string filenames[sizefile] = {"G:\\dataset\\TimeSeries\\real\\IBM_data", "G:\\dataset\\TimeSeries\\real\\ICU_data", "G:\\dataset\\TimeSeries\\real\\Traffic_data"};
	string resultFilename = "NULL";
	string logFile = "G:\\dataset\\TimeSeries\\real\\logSegSizeTest.txt";

	vector<double> thresholds;

	const size_t size = 100;

	double threshold = 0.5;
	double step = 0.1;

	ofstream out(logFile.c_str(), ios::app);

	for(size_t i = 0; i < sizefile; ++i)
	{
		vector<double> vec;
		ReadDataFromFile(vec, filenames[i]);
		cout << "Size:" << vec.size() << endl;

		out << "Filename:" << filenames[i] << endl;
		out << "size = [";
		for(size_t j = 0; j < size; ++j)
		{
			double thr = threshold + j * step;

			ContinuousBottomUp *bottomUp = new ContinuousBottomUp(thr, vec);
			bottomUp->SetOutputFilepath(resultFilename);
			bottomUp->Approximate();
			long long segSize = bottomUp->GetSegmentSize();

			out << segSize << ","; 
		}
		out << "];" << endl << endl;
	}

	out.close();

}

void ExperimentWindowsizeBatchsize()
{
	int windowSize = 100;
	int windowSizeStep = 50;
	int batchSize = 10;
	int batchSizeStep = 5;

	string filename = "G:\\dataset\\TimeSeries\\sythetic\\sythetic2";
	string log = "G:\\dataset\\TimeSeries\\sythetic\\windowSizeBatchSizelog.txt";
	string resultFilename = "NULL";

	ofstream outlog(log.c_str());

	outlog << "Y = [";
	for(size_t i = 0; i < 10; ++i)
	{
		int curWindowSize = windowSize + i * windowSizeStep;
		
		for(size_t j = 0; j < 10; ++j)
		{
			int curBatchSize = batchSize + i * batchSize;

			ifstream *in = new ifstream(filename.c_str());
			ofstream *out = new ofstream(resultFilename.c_str());
			

			clock_t start, end;
			start = clock();

			MDLSlideWindow *mdlsw = new MDLSlideWindow(curWindowSize, in, out);
			mdlsw->SetOutput(false);
			mdlsw->SetBatchSize(curBatchSize);
			mdlsw->Approximate();

			end = clock();

			long double lasts = ((double)(end - start)) / (double)CLOCKS_PER_SEC;
			
			outlog << lasts << " ";
			
		}

		outlog << ";\n";
	}
	outlog << "];";

	outlog.close();

}