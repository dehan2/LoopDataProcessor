#pragma once


class CarouselObject;
struct SegmentTransitionEvent
{
	double time;
	CarouselObject* corrObj;
};


struct compare_segment_transition_events_in_ascending_order {
	bool operator()(const SegmentTransitionEvent& lhs, const SegmentTransitionEvent& rhs)
	{
		return lhs.time > rhs.time;
	}
};


const double SECONDS_PER_TCA = 1.0;
const double OBJECT_SIZE = 1.0E-10;
const int START_ID_OF_OBJECT = 1;


struct CarouselObjectData
{
	int ID;
	double radius;
	double period;
};