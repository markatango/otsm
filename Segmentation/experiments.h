#include "headers.h"
#include "segmenting.h"

void Experiment();

void ExperimentGenerateData(std::string pathPrefix);

void ExperimentEfficient();
void ExperimentEfficientGlobal(std::string pathPrefix);
void ExperimentEfficientSlideWindow(std::string pathPrefix);
void ExperimentEfficientMDLSlideWindow(std::string pathPrefix);



void ExperimentPrecision(const std::string &pathPrefix);

//	Run bottom-up on datasets
void ExperimentBottomUp(const std::string &filename, const std::string &resultFilename, const std::string &log, double threshold);
//	Run PRESEE on datasets
void ExperimentMDLSlideWindow(const std::string &filename, const std::string &resultFilename, const std::string &log);

//	Measure precision on datasets
void ExperimentMeasure(const std::string &filename, const std::string &resultFilename, const std::string &log);

//	Mearure the relationship between error threshold and number of segments
void ExperimentThresholdTest();

//	Test the efficiency between window size and batch size, window size from 100 to 1000 step 100, batch size from 1 to 5 step 1
void ExperimentWindowsizeBatchsize();