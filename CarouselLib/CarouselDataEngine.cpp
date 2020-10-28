#include "CarouselDataEngine.h"
#include "ConstForCarouselLib.h"
#include <random>

#define _USE_MATH_DEFINES
#include <math.h>

CarouselDataEngine::CarouselDataEngine()
	:DynamicVDEngine()
{
	m_dynamicVD.set_process_collision_option(false);
}



CarouselDataEngine::CarouselDataEngine(const int& ID)
	:CarouselDataEngine()
{
	m_ID = ID;
}



CarouselDataEngine::~CarouselDataEngine()
{
	clear();
}



void CarouselDataEngine::clear()
{
	DynamicVDEngine::clear();
}



void CarouselDataEngine::load_carousel_data_file(const string& filePath)
{
	/*ifstream fin;
	fin.open(filePath);

	if (fin.is_open())
	{
		char lineData[256];
		fin.getline(lineData, 256);
		while (lineData[0] == '%')
		{
			fin.getline(lineData, 256);
		}

		//Parse command
		char* context;
		string delimiter = " \t";

		string token = strtok_s(lineData, delimiter.c_str(), &context);
		m_approxLevel = stoi(token);

		token = strtok_s(NULL, delimiter.c_str(), &context);
		m_numSecondaryObjects = stoi(token);

		token = strtok_s(NULL, delimiter.c_str(), &context);
		m_missDistance = stof(token);

		token = strtok_s(NULL, delimiter.c_str(), &context);
		m_secondaryCarouselRadius = stof(token);

		token = strtok_s(NULL, delimiter.c_str(), &context);
		m_minSecondaryCarouselPeriod = stof(token);

		token = strtok_s(NULL, delimiter.c_str(), &context);
		m_maxSecondaryCarouselPeriod = stof(token);

		token = strtok_s(NULL, delimiter.c_str(), &context);
		m_distanceBetweenSecondaryCarousels = stof(token);

		list<CarouselObjectData> secondaryInfos;
		for (int i = 0; i < m_numSecondaryObjects; i++)
		{
			CarouselObjectData secondaryInfo;
			fin.getline(lineData, 256);
			
			string token = strtok_s(lineData, delimiter.c_str(), &context);
			secondaryInfo.ID = stoi(token);

			token = strtok_s(NULL, delimiter.c_str(), &context);
			secondaryInfo.radius = stof(token);

			token = strtok_s(NULL, delimiter.c_str(), &context);
			secondaryInfo.period = stof(token);

			secondaryInfos.push_back(secondaryInfo);
		}

		reproduce_carousel_data(secondaryInfos);
	}

	fin.close();*/
}



void CarouselDataEngine::save_carousel_data_file(const string& filePath)
{
	/*ofstream fout(filePath);
	fout << "% Carousel Data\n";
	fout << "% ApproxLevel\tnumSecondaryObjects\tmissDistance\tSecCarRadius\tMinSecCarPeriod\tMaxSecCarPeriod\tDistanceBtwSecCars\n";
	fout << m_approxLevel << "\t" << m_numSecondaryObjects << "\t" <<m_missDistance<<"\t"<< m_secondaryCarouselRadius << "\t" << m_minSecondaryCarouselPeriod<<"\t"<<m_maxSecondaryCarouselPeriod << "\t" << m_distanceBetweenSecondaryCarousels<<"\n";
	for (auto& object : m_carouselObjects)
	{
		if(object.get_ID() > 1)
			fout << object.get_ID() << "\t" << object.get_carousel_radius() << "\t" << object.get_period() << "\n";
	}
	fout.close();*/
}



string CarouselDataEngine::make_carousel_data_file_name(const int& numSecondaryObjects, const int& approxLevel)
{
	string carouselDataFileName = string("CD_") + to_string(numSecondaryObjects) + string("_") + to_string(approxLevel) + (".cdf");
	return carouselDataFileName;
}



void CarouselDataEngine::propagate_EVENTSEQ(const int& warehouseSize)
{
	add_approx_segment_transition_events_within_time_window();
	try
	{
		DynamicVDEngine::propagate_EVENTSEQ(warehouseSize);
	}
	catch (DVD_ERROR error)
	{
		throw error;
	}
}



void CarouselDataEngine::save_event_sequence_for_carousel_data(const int& numSecondaryObjects, const int& startTime, const int& predictionTimeWindow) const
{
	string fileName = make_event_sequence_file_name_for_carousel_data(numSecondaryObjects, startTime, predictionTimeWindow);
	save_event_sequence(fileName);
}



string CarouselDataEngine::make_event_sequence_file_name_for_carousel_data(const int& numSecondaryObjects, const int& startTime, const int& predictionTimeWindow) const
{
	string carouselDataFileName = string("EvSeq_CD_") + to_string(numSecondaryObjects) + string("_") + to_string(startTime) + string("_") + to_string (predictionTimeWindow) + string(".evseq");
	return carouselDataFileName;
}



void CarouselDataEngine::process_next_custom_event_for_event_sequence_propagation()
{
	const SegmentTransitionEvent& segmentationTransitionEvent = m_approxSegmentTransitionQ.top();
	double segmentationTransitionEventTime = segmentationTransitionEvent.time;

	CarouselObject* object = segmentationTransitionEvent.corrObj;
	//cout << "Segment transition : " << segmentChangeBall->get_ID() << " at " << nextSegmentChangeEventInfo.second << endl;
	m_dynamicVD.add_velocity_change_event(object->get_ID(), segmentationTransitionEvent.time);

	double secondsPerApproxSegment = object->calculate_time_per_segment();

	m_approxSegmentTransitionQ.pop();
	m_approxSegmentTransitionQ.push({ segmentationTransitionEventTime + secondsPerApproxSegment, object });
}



double CarouselDataEngine::find_imminent_custom_event_time_for_event_sequence_propagation() const
{
	if (m_approxSegmentTransitionQ.empty() == false)
	{
		const SegmentTransitionEvent& segmentationTransitionEvent = m_approxSegmentTransitionQ.top();
		return segmentationTransitionEvent.time;
	}
	else
	{
		return DVD_MAX;
	}
}



void CarouselDataEngine::initialize_carousel_engine(const double& predictionTimeWindow)
{
	m_status = DVDOS_STATUS::INITIALIZE;
	
	enlist_carousel_objects_to_DVD();

	add_center_ball();
	double carouselSpaceSize = 4*m_carouselObjects.front().get_carousel_radius();
	add_bounding_static_objects(carouselSpaceSize, OBJECT_SIZE, m_carouselObjects.size()+ START_ID_OF_OBJECT+1);

	m_dynamicVD.set_prediction_time_window(predictionTimeWindow);
	m_dynamicVD.construct_initial_VD(m_objects);

	/*string qtfFilePath = string("QTF_")+to_string(m_numSecondaryObjects)+string(".qtf");
	bool isQTFAvailable = construct_initial_VD_from_QTF_file(qtfFilePath);
	if (isQTFAvailable == false)
	{
		m_dynamicVD.construct_initial_VD(m_objects);
		save_quasi_triangulation_file(qtfFilePath);
	}*/
}



void CarouselDataEngine::initialize_carousel_engine(const list<CarouselObjectData>& secondaryInfos, int numApproxSegments, const double& startTime, const double& predictionTimeWindow)
{
	m_status = DVDOS_STATUS::INITIALIZE;

	initialize_carousel_objects(secondaryInfos, startTime);
	enlist_carousel_objects_to_DVD();

	add_center_ball();
	double carouselSpaceSize = 4 * m_carouselObjects.front().get_carousel_radius();
	add_bounding_static_objects(carouselSpaceSize, OBJECT_SIZE, m_carouselObjects.size() + START_ID_OF_OBJECT + 1);

	m_dynamicVD.set_prediction_time_window(predictionTimeWindow);
	//m_dynamicVD.construct_initial_VD(m_objects);
}



void CarouselDataEngine::initialize_carousel_objects(const list<CarouselObjectData>& secondaryInfos, const double& startTime)
{
	double primaryCarouselPeriod = m_numSecondaryObjects * SECONDS_PER_TCA;
	double primaryCarouselRadius = m_distanceBetweenSecondaryCarousels * m_numSecondaryObjects * 0.5 / M_PI;
	double primaryAngularPosition = 2 * M_PI / primaryCarouselPeriod * startTime;

	m_carouselObjects.push_back(CarouselObject(START_ID_OF_OBJECT, primaryAngularPosition, rg_Point3D(0, 0, 0), primaryCarouselRadius, primaryCarouselPeriod, rg_Point3D(1, 0, 0), rg_Point3D(0, 1, 0)));

	//cout << "Primary radius: " << primaryCarouselRadius << endl;

	double distanceToSecondaryCarouselCenter = primaryCarouselRadius + m_missDistance + m_secondaryCarouselRadius;

	for (auto& secondaryInfo : secondaryInfos)
	{
		double angleToSecondaryCarouselCenter = 2 * M_PI * (secondaryInfo.ID - 1) / m_numSecondaryObjects;
		rg_Point3D secondaryCarouselCenter(cos(angleToSecondaryCarouselCenter), sin(angleToSecondaryCarouselCenter), 0);
		secondaryCarouselCenter = secondaryCarouselCenter * distanceToSecondaryCarouselCenter;

		double secondaryPeriod = secondaryInfo.period;
		//cout << "Period of " << i + 2 << " object: " << period << endl;
		double secondaryAngularVelocity = 2 * M_PI / secondaryPeriod;

		rg_Point3D secondaryCarouselUVector = -secondaryCarouselCenter;
		secondaryCarouselUVector.normalize();
		rg_Point3D secondaryCarouselVVector(0, 0, 1);
		double secondaryAngularPosition = 2 * M_PI - (secondaryInfo.ID-1) * secondaryAngularVelocity; //angular position at 0 time.
		secondaryAngularPosition += startTime * secondaryAngularVelocity;

		m_carouselObjects.push_back(CarouselObject(secondaryInfo.ID, secondaryAngularPosition, secondaryCarouselCenter, m_secondaryCarouselRadius,
			secondaryPeriod, secondaryCarouselUVector, secondaryCarouselVVector));
	}

	for (auto& carouselObject : m_carouselObjects)
	{
		carouselObject.initialize_carousel_object(m_approxLevel, 0.0);

		//for debug
		//double targetTime = carouselObject.get_ID() - 1 - startTime;
		//rg_Point3D primaryAtTime = m_carouselObjects.front().calculate_position_of_original_at_time(targetTime);
		//rg_Point3D seoncdaryAtTime = carouselObject.calculate_position_of_original_at_time(targetTime);
		//cout << "CO[" << carouselObject.get_ID() << "] - distance at "<<targetTime<<": "<<primaryAtTime.distance(seoncdaryAtTime)<< endl;


		//rg_Point3D center = carouselObject.get_sphere_center();
		//rg_Point3D carouselCenter = carouselObject.get_carousel_center();
		//cout << "CO[" << carouselObject.get_ID() << "] - center[" << center.getX() << ", " << center.getY() << ", " << center.getZ() << "], CC: [" << carouselCenter.getX() << ", " << carouselCenter.getY() << ", " << carouselCenter.getZ() << "], distance: [" <<carouselCenter.magnitude()<<", "<<carouselCenter.distance(center)<<"]"<< endl;


		//rg_Point3D center = carouselObject.get_sphere_center();
		//rg_Point3D velocity = carouselObject.get_velocity();
		//cout << "CO[" << carouselObject.get_ID() << "] - center[" << center.getX() << ", " << center.getY() << ", " << center.getZ() << "], Velocity: [" << velocity.getX() << ", " << velocity.getY() << ", " << velocity.getZ() << "]" << endl;
	}
}



void CarouselDataEngine::add_center_ball()
{
	int centerBallID = m_carouselObjects.size() + START_ID_OF_OBJECT + m_staticObjects.size();
	m_staticObjects.push_back(DynamicObject(centerBallID, DYNAMIC_OBJECT_TYPE::STATIC_OBJECT, Sphere(rg_Point3D(), OBJECT_SIZE), rg_Point3D()));
	m_objects.push_back(&m_staticObjects.back());
}



void CarouselDataEngine::enlist_carousel_objects_to_DVD()
{
	for (auto& carouselObject : m_carouselObjects)
	{
		m_objects.push_back(&carouselObject);
		m_mapFromIDToDynamicObjects[carouselObject.get_ID()] = &carouselObject;
	}
}



void CarouselDataEngine::construct_initial_VD_and_save_as_qtf_file(const string& qtfFileName)
{
	m_dynamicVD.construct_initial_VD(m_objects);
	save_quasi_triangulation_file(qtfFileName);
}



string CarouselDataEngine::make_qtf_file_name_for_carousel_data(const int& numObject, const int& timeFragmentStart, const int& timeFragmentSize) const
{
	string eventSeqFileName = make_event_sequence_file_name_for_carousel_data(numObject, timeFragmentStart, timeFragmentSize);
	string qtfFileName = eventSeqFileName + "_qt.qtf";
	return qtfFileName;
}



void CarouselDataEngine::add_approx_segment_transition_events_within_time_window()
{
	for (auto& carouselObject : m_carouselObjects)
	{
		double secondsPerSegment = carouselObject.calculate_time_per_segment();
		if (secondsPerSegment < m_dynamicVD.get_prediction_time_window())
		{
			m_approxSegmentTransitionQ.push({ secondsPerSegment, &carouselObject });
		}
	}
}



double CarouselDataEngine::calculate_real_distance_for_time(const CarouselObject* primary, const CarouselObject* secondary, const double& targetTime)
{
	rg_Point3D primaryCoord = primary->calculate_position_of_original_at_time(targetTime);
	rg_Point3D secondaryCoord = secondary->calculate_position_of_original_at_time(targetTime);
	double distance = primaryCoord.distance(secondaryCoord);
	return distance;
}



void CarouselDataEngine::generate_carousel_data(int numSecondaryObjects, double missDistance, double secondaryCarouselRadius, double secondaryCarouselPeriod, double distanceBetweenSecondaryCarousels)
{
	cout << "Generate Carousel data" << endl;

	m_numSecondaryObjects = numSecondaryObjects;
	m_missDistance = missDistance;
	m_secondaryCarouselRadius = secondaryCarouselRadius;
	//m_secondaryCarouselPeriod = secondaryCarouselPeriod;
	m_distanceBetweenSecondaryCarousels = distanceBetweenSecondaryCarousels;


	double primaryCarouselPeriod = numSecondaryObjects * SECONDS_PER_TCA;
	double primaryCarouselRadius = distanceBetweenSecondaryCarousels * numSecondaryObjects * 0.5 / M_PI;
	m_carouselObjects.push_back(CarouselObject(START_ID_OF_OBJECT, 0, rg_Point3D(0, 0, 0), primaryCarouselRadius, primaryCarouselPeriod, rg_Point3D(1, 0, 0), rg_Point3D(0, 1, 0)));

	cout << "Primary radius: " << primaryCarouselRadius << endl;

	//Randomize secondary period

	//m_minSecondaryCarouselPeriod = 500;
	//m_maxSecondaryCarouselPeriod = 1000;

	random_device rd;
	mt19937 generator(rd());
	uniform_real_distribution<double> distribution(500, 1000);

	double distanceToSecondaryCarouselCenter = primaryCarouselRadius + missDistance + secondaryCarouselRadius;
	cout << "Distance to Sec Car Center: " << distanceToSecondaryCarouselCenter << endl;

	for (int i = 0; i < numSecondaryObjects; i++)
	{
		double angleToSecondaryCarouselCenter = 2 * M_PI * (i + 1) / numSecondaryObjects;
		rg_Point3D secondaryCarouselCenter(cos(angleToSecondaryCarouselCenter), sin(angleToSecondaryCarouselCenter), 0);
		secondaryCarouselCenter = secondaryCarouselCenter * distanceToSecondaryCarouselCenter;

		double period = distribution(generator);
		//cout << "Period of " << i + 2 << " object: " << period << endl;
		double secondaryAngularVelocity = 2 * M_PI / period;

		rg_Point3D secondaryCarouselUVector = -secondaryCarouselCenter;
		secondaryCarouselUVector.normalize();
		rg_Point3D secondaryCarouselVVector(0, 0, 1);
		double secondaryAngularPosition = 2 * M_PI - (i + 1) * secondaryAngularVelocity;
		
		m_carouselObjects.push_back(CarouselObject(START_ID_OF_OBJECT+i+1, secondaryAngularPosition, secondaryCarouselCenter, secondaryCarouselRadius,
			period, secondaryCarouselUVector, secondaryCarouselVVector));
	}

	for (auto& carouselObject : m_carouselObjects)
	{
		carouselObject.initialize_carousel_object(m_approxLevel);

		//for debug
		//rg_Point3D center = carouselObject.get_sphere_center();
		//rg_Point3D velocity = carouselObject.get_velocity();
		//cout << "CO[" << carouselObject.get_ID() << "] - center[" << center.getX() << ", " << center.getY() << ", " << center.getZ() << "], Velocity: [" << velocity.getX() << ", " << velocity.getY() << ", " << velocity.getZ() << "]" << endl;
	}
}



void CarouselDataEngine::reproduce_carousel_data(const list<CarouselObjectData>& carouselObjectInfos)
{
	cout << "Reproduce Carousel data" << endl;

	double primaryCarouselPeriod = m_numSecondaryObjects * SECONDS_PER_TCA;
	double primaryCarouselRadius = m_distanceBetweenSecondaryCarousels * m_numSecondaryObjects * 0.5 / M_PI;
	m_carouselObjects.push_back(CarouselObject(START_ID_OF_OBJECT, 0, rg_Point3D(0, 0, 0), primaryCarouselRadius, primaryCarouselPeriod, rg_Point3D(1, 0, 0), rg_Point3D(0, 1, 0)));

	cout << "Primary radius: " << primaryCarouselRadius << endl;

	for (auto& objectInfo : carouselObjectInfos)
	{
		int ID = objectInfo.ID;
		double radius = objectInfo.radius;
		double distanceToSecondaryCarouselCenter = primaryCarouselRadius + m_missDistance + radius;

		double angleToSecondaryCarouselCenter = 2 * M_PI * (ID - 1) / m_numSecondaryObjects;
		rg_Point3D secondaryCarouselCenter(cos(angleToSecondaryCarouselCenter), sin(angleToSecondaryCarouselCenter), 0);
		secondaryCarouselCenter = secondaryCarouselCenter * distanceToSecondaryCarouselCenter;

		double period = objectInfo.period;
		//cout << "Period of " << i + 2 << " object: " << period << endl;
		double secondaryAngularVelocity = 2 * M_PI / period;

		rg_Point3D secondaryCarouselUVector = -secondaryCarouselCenter;
		secondaryCarouselUVector.normalize();
		rg_Point3D secondaryCarouselVVector(0, 0, 1);
		double secondaryAngularPosition = 2 * M_PI - (ID -1) * secondaryAngularVelocity;

		m_carouselObjects.push_back(CarouselObject(ID, secondaryAngularPosition, secondaryCarouselCenter, radius,
			period, secondaryCarouselUVector, secondaryCarouselVVector));
	}

	for (auto& carouselObject : m_carouselObjects)
	{
		carouselObject.initialize_carousel_object(m_approxLevel);

		//for debug
		//rg_Point3D center = carouselObject.get_sphere_center();
		//rg_Point3D velocity = carouselObject.get_velocity();
		//cout << "CO[" << carouselObject.get_ID() << "] - center[" << center.getX() << ", " << center.getY() << ", " << center.getZ() << "], Velocity: [" << velocity.getX() << ", " << velocity.getY() << ", " << velocity.getZ() << "]" << endl;
	}
}
