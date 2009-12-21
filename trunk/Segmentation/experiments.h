#include "headers.h"
#include "segmenting.h"

void Experiment();

void ExperimentGenerateData(std::string pathPrefix);

void ExperimentEfficient();
void ExperimentEfficientGlobal(std::string pathPrefix);
void ExperimentEfficientSlideWindow(std::string pathPrefix);
void ExperimentEfficientMDLSlideWindow(std::string pathPrefix);

void ExperimentPrecision();