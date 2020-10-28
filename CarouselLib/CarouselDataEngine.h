#pragma once

#include "DynamicVDEngine.h"
#include "CarouselObject.h"
#include <list>

class CarouselDataEngine : public DynamicVDEngine
{
protected:
	list<CarouselObject>	m_carouselObjects;
	int m_approxLevel = 10;

	int m_numSecondaryObjects = 0;
	double m_missDistance = 0.0;
	double m_secondaryCarouselRadius = 0.0;
	double m_distanceBetweenSecondaryCarousels = 0.0;
	//double m_minSecondaryCarouselPeriod = 0.0;
	//double m_maxSecondaryCarouselPeriod = 0.0;

	priority_queue<SegmentTransitionEvent, vector<SegmentTransitionEvent>, compare_segment_transition_events_in_ascending_order> m_approxSegmentTransitionQ;

public:
	CarouselDataEngine();
	CarouselDataEngine(const int& ID);
	virtual ~CarouselDataEngine();

	virtual void clear();

	inline const list<CarouselObject>& get_carousel_objects() const { return m_carouselObjects; }
	inline const int&	get_approximation_level() const { return m_approxLevel; }
	inline const int& get_num_secondary_objects() const { return m_numSecondaryObjects; }
	inline const double& get_secondary_carousel_radius() const { return m_secondaryCarouselRadius; }
	inline const double& get_miss_distance() const { return m_missDistance; }
	inline const double& get_distance_between_secondary_carousels() const { return m_distanceBetweenSecondaryCarousels; }

	inline void set_approximation_level(const int& approxLevel) { m_approxLevel = approxLevel; }
	inline void set_num_secondary_objects(const int& numSecondaryObjects) { m_numSecondaryObjects = numSecondaryObjects; }
	inline void set_secondary_carousel_radius(const double& secondaryCarouselRadius) { m_secondaryCarouselRadius = secondaryCarouselRadius; }
	inline void set_miss_distance(const double& missDistance) { m_missDistance = missDistance; }
	inline void set_distance_between_secondary_carousels(const double& distanceBetweenSecondaryCarousels) { m_distanceBetweenSecondaryCarousels = distanceBetweenSecondaryCarousels; }


	//////////////////////////////////////////////////////////////////////////
	// PART 1. File Load/Save - Inherited
	//////////////////////////////////////////////////////////////////////////

	virtual void	load_carousel_data_file(const string& filePath);
	virtual void	save_carousel_data_file(const string& filePath);
	string make_carousel_data_file_name(const int& numSecondaryObjects, const int& approxLevel);

	//////////////////////////////////////////////////////////////////////////
	// PART 2. Simulation Control - Inherited
	//////////////////////////////////////////////////////////////////////////

	//2-1. Event Seq. Generation
	virtual void	propagate_EVENTSEQ(const int& warehouseSize);
	void save_event_sequence_for_carousel_data(const int& numSecondaryObjects, const int& startTime, const int& predictionTimeWindow) const;
	string make_event_sequence_file_name_for_carousel_data(const int& numSecondaryObjects, const int& startTime, const int& predictionTimeWindow) const;

	virtual void	process_next_custom_event_for_event_sequence_propagation();
	virtual double	find_imminent_custom_event_time_for_event_sequence_propagation() const;


	//////////////////////////////////////////////////////////////////////////
	// PART 3. Carousel Engine
	//////////////////////////////////////////////////////////////////////////

	//3-1. Initialization

	void initialize_carousel_engine(const double& predictionTimeWindow);
	void initialize_carousel_engine(const list<CarouselObjectData>& secondaryInfos, int numApproxSegments, const double& startTime, const double& predictionTimeWindow);

	void initialize_carousel_objects(const list<CarouselObjectData>& secondaryInfos, const double& startTime);
	
	void add_center_ball();
	void enlist_carousel_objects_to_DVD();
	
	void construct_initial_VD_and_save_as_qtf_file(const string& qtfFileName);
	string make_qtf_file_name_for_carousel_data(const int& numObject, const int& timeFragmentStart, const int& timeFragmentSize) const;


	void add_approx_segment_transition_events_within_time_window();

	//3-2. Utilities
	static double		calculate_real_distance_for_time(const CarouselObject* primary, const CarouselObject* secondary, const double& targetTime);

	//3-3. Generate Carousel Data
	void generate_carousel_data(int numSecondaryObjects, double missDistance, double secondaryCarouselRadius, double secondaryCarouselPeriod, double distanceBetweenSecondaryCarousels);
	void reproduce_carousel_data(const list<CarouselObjectData>& carouselObjectInfos);
};

