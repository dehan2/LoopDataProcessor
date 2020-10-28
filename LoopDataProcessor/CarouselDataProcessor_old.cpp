#include <iostream>
#include <chrono>

#include "CarouselObject.h"
#include "CarouselDataEngine.h"

int main()
{
	int approxLevel = 100;

	int numSecondaryObjects = 1000;
	double missDistance = 0.1;
	double secondaryCarouselRadius = 10.0;
	double secondaryCarouselPeriod = 10.0;
	double distanceBetweenSecondaryCarousels = 1;

    CarouselDataEngine engine;

	/*const list<CarouselObject>& objects = engine.get_carousel_objects();

	const CarouselObject& primary = objects.front();

	for (int i = 0; i < 10 * numSecondaryObjects; i++)
	{
		double time = i * 1.0e-1;
		cout << "Time: " << time << endl;
		rg_Point3D primaryCoord = primary.calculate_position_of_original_at_time(time);

		for (auto& obj : objects)
		{
			rg_Point3D coord = obj.calculate_position_of_original_at_time(time);
			double distanceFromCenter = coord.magnitude();
			double distanceToPrimary = coord.distance(primaryCoord);
			cout << "Obj[" << obj.get_ID() << "]: distance: " << distanceFromCenter << ", distance to primary: "<<distanceToPrimary <<", center: [" << coord.getX() << ", " << coord.getY() << ", " << coord.getZ() << "]" << endl;
		}
	}*/


	//Process A. Event Seq Generation
	/*engine.set_approximation_level(approxLevel);
	engine.generate_carousel_data(numSecondaryObjects, missDistance, secondaryCarouselRadius, secondaryCarouselPeriod, distanceBetweenSecondaryCarousels);
	
	cout << "Step1: Initialization" << endl;
	Clock::time_point begin = Clock::now();
	engine.initialize_carousel_engine(numSecondaryObjects);
	Clock::time_point end = Clock::now();
	double timeElapsed = chrono::duration_cast<chrono::nanoseconds>(end - begin).count() * 1.0e-6;
	cout << "Initialization Finish: " << timeElapsed << "ms" << endl;

	cout << "Step2: EvSeq propagation" << endl;
	begin = Clock::now();
	engine.propagate_EVENTSEQ(0);
	end = Clock::now();
	timeElapsed = chrono::duration_cast<chrono::nanoseconds>(end - begin).count() * 1.0e-6;
	cout << "Event Seq Propagation Finish: " << timeElapsed << "ms" << endl;

	cout << "Step3: Validation" << endl;
	begin = Clock::now();
	engine.go_to_past_using_event_sequence(0.0);
	end = Clock::now();
	timeElapsed = chrono::duration_cast<chrono::nanoseconds>(end - begin).count() * 1.0e-6;
	cout << "Validation Finish: " << timeElapsed << "ms" << endl;

	string eventSeqFileName = engine.make_event_sequence_file_name(numSecondaryObjects, approxLevel);
	string carouselDataFileName = engine.make_carousel_data_file_name(numSecondaryObjects, approxLevel);
	engine.save_carousel_data_file(carouselDataFileName);
	engine.save_event_sequence(eventSeqFileName);*/


	//Process B. Generate PPDB
	double terminateCondition = 1.0e-6;
	double cutoffDistance = 0.2;
	
	string carouselDataFileName = engine.make_carousel_data_file_name(numSecondaryObjects, approxLevel);
	
	cout << "Step1: Initialization" << endl;
	Clock::time_point begin = Clock::now();
	
	engine.load_carousel_data_file(carouselDataFileName);
	engine.initialize_carousel_engine(numSecondaryObjects);
	string eventSeqFileName = engine.make_event_sequence_file_name(numSecondaryObjects, approxLevel);
	engine.load_event_sequence_file(eventSeqFileName);
	
	Clock::time_point end = Clock::now();
	double timeElapsed = chrono::duration_cast<chrono::nanoseconds>(end - begin).count() * 1.0e-6;
	cout << "Initialization Finish: " << timeElapsed << "ms" << endl;
	
	cout << "Step2: Generate PPDB" << endl;
	begin = Clock::now();

	engine.set_cutoff_value_for_CA(terminateCondition);
	engine.generate_pairwise_proximity_database(cutoffDistance);

	end = Clock::now();
	timeElapsed = chrono::duration_cast<chrono::nanoseconds>(end - begin).count() * 1.0e-6;
	cout << "Generate PPDB Finish: " << timeElapsed << "ms" << endl;

	const list<SearchInterval>& searchIntervals = engine.get_search_intervals();

	int index = 0;
	for (auto& searchInterval : searchIntervals)
	{
		int primaryID = searchInterval.primary->get_ID();
		int secondaryID = searchInterval.secondary->get_ID();
		if (primaryID > secondaryID)
			swap(primaryID, secondaryID);

		cout << "CA[ " << ++index << " ]: pair [ " << primaryID << ", " << secondaryID << " ], TCA: " << searchInterval.timeOfClosestApproach << ", distance: " << searchInterval.minDistance;
		cout << ", SI: [ " << searchInterval.startTime << ", " << searchInterval.endTime << " ]" << endl;
	}
}





/*const list<CarouselObject>& objects = engine.get_carousel_objects();

const CarouselObject& primary = objects.front();

vector<pair<double, double>> minDistance(numSecondaryObjects);
for (int i = 0; i < numSecondaryObjects; i++)
{
	minDistance.at(i) = { DBL_MAX, 0 };
}

for (int i = 0; i < 1e3 * numSecondaryObjects; i++)
{
	double time = i * 1.0e-3;
	//cout << "Time: " << time << endl;
	rg_Point3D primaryCoord = primary.calculate_position_of_original_at_time(time);

	for (auto& obj : objects)
	{
		rg_Point3D coord = obj.calculate_position_of_original_at_time(time);
		double distanceFromCenter = coord.magnitude();
		double distanceToPrimary = coord.distance(primaryCoord);
		//cout << "Obj[" << obj.get_ID() << "]: distance: " << distanceFromCenter << ", distance to primary: " << distanceToPrimary << ", center: [" << coord.getX() << ", " << coord.getY() << ", " << coord.getZ() << "]" << endl;

		if (obj.get_ID() > 1 && distanceToPrimary < minDistance.at(obj.get_ID() - 2).first)
		{
			minDistance.at(obj.get_ID() - 2) = { distanceToPrimary, time };
		}
	}
}

for (int i = 0; i < numSecondaryObjects; i++)
{
	cout << "Min " << i + 2 << " : " << minDistance.at(i).first << " , " << minDistance.at(i).second << endl;
}*/