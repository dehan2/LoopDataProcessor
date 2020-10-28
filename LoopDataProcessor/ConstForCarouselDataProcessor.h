#pragma once

#include <string>

using namespace std;

const float INITIAL_TIME_FRAG_SIZE = 100;
const float TCA_INTEGRATION_THRESHOLD = 1.0e-1;


enum class CD_OPERATION_MODE
{
	EVENT_SEQ_GENERATION,
	CLOSE_APPROACH_SCAN
};




struct TCAReport
{
	int primaryID;
	int secondaryID;
	double distanceOfClosestApproach;
	double timeOfClosestApproach;
	double closeApproachEnteringTime;
	double closeApproachLeavingTime;
};



struct FailSegmentInfo
{
	float startTime;
	float timeWindow;
	double processedTime;
	string message;
};


const int TIME_FRAGMENT_SCALE = 1;
const double TIME_FRAGMENT_SCALE_REVERSE = 1;