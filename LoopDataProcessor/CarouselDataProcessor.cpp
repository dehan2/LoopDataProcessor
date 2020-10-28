#include <iostream>
#include "ConstForCarouselDataProcessor.h"
#include "CarouselDataManager.h"

using namespace std;

int main(int argc, char** argv)
{
	CarouselDataManager manager;

	//Argument 1 : COOP Mode - E/C/M
	char* modeArg = argv[1];
	CD_OPERATION_MODE COOPMode;
	switch (modeArg[0])
	{
	case 'E':
	default:
		COOPMode = CD_OPERATION_MODE::EVENT_SEQ_GENERATION;
		cout << "Mode: Event Sequence Generation" << endl;
		break;
	case 'C':
		COOPMode = CD_OPERATION_MODE::CLOSE_APPROACH_SCAN;
		cout << "Mode: Close Approach Scan" << endl;
		break;
	}



	//Argument 2 : Start Time
	double predictionStartTime = stof(argv[2]);

	//Argument 3 : Number of threads to be used in COOP
	int maxThreads = stoi(argv[3]);

	//Argument 4 : Close Approach Scan Threshold 
	/*double cutoffDistance = 0;
	if (argc > 4)
		cutoffDistance = stof(argv[4]);*/

	double cutoffDistance = 0.2;

	switch (COOPMode)
	{
	case CD_OPERATION_MODE::EVENT_SEQ_GENERATION:
	default:
		manager.propagate_EVENTSEQ_through_parallel_computation(maxThreads, INITIAL_TIME_FRAG_SIZE);
		break;
	case CD_OPERATION_MODE::CLOSE_APPROACH_SCAN:
		manager.generate_PPDB_through_parallel_computation(maxThreads, cutoffDistance);
		break;
	}

	cout << "Process finish" << endl;

	return 0;
}