#include "experiments.h"

using namespace std;

int main()
{

	cout << "Motif." << endl;
	
	//Test();
	//NaiveLongTailExperiments();

	//SIMMALongTailExperiments();

	//	Test the relationship between R and efficiency
	//SIMMARvsMotifSize();
	//NaiveRvsMotifSize();

	//	Test the relationship between slide window size and efficiency
	//SIMMASlideWindowSizevsMotifSize();
	//NaiveSlideWindowSizevsMotifSize();

	//	Test the relationship between motif length and efficiency
	//SIMMAEfficiencyMotifLength();
	//NaiveEfficiencyMotifLength();

	//	Test the efficiency of naive and SIMMA
	NaiveEfficiency();
	SIMMAEfficiency();

	//	Test the relationship between R and efficiency
	NaiveREfficiency();
	SIMMAREfficiency();

	//	Test SIMMA shrink time
	SIMMAShrinkTest();

	return 0;
}