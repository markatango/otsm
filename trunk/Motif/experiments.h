#include "motif.h"
#include "tools.h"

void Test();

//	Test the long tail phenomeno of motif frequency in time series stream
void NaiveLongTailExperiments();
void NaiveLongTailExperiment(std::string infile, std::string outfile, const motif::ParametersNaive &params);

//	Test the long tail phenomeno of motif frequency through SIMMA
void SIMMALongTailExperiments();
void SIMMALongTailExperiment(std::string infile, std::string  outfile, const motif::ParametersSIMMA &params);

//	Test the relationship between R and efficiency
void SIMMARvsMotifSize();
void NaiveRvsMotifSize();

//	Test the relationship between slide window size and efficiency
void SIMMASlideWindowSizevsMotifSize();
void NaiveSlideWindowSizevsMotifSize();

//	Test the relationship between motif length and efficiency
void SIMMAEfficiencyMotifLength();
void NaiveEfficiencyMotifLength();


//	Test the relationship between motif number and dataset size
void SIMMAMotifNumber();
void NaiveMotifNumber();

//	Compare efficiency between naive and SIMMA
void SIMMAEfficiency();
void NaiveEfficiency();

//	Test the relationship between R and efficiency
void SIMMAREfficiency();
void NaiveREfficiency();

//	Test SIMMA shrink time
void SIMMAShrinkTest();