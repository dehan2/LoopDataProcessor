#pragma once

#include <thread>
#include <mutex>
#include "CarouselDataEngine.h"
#include "ConstForCarouselDataProcessor.h"

using namespace std;

class CarouselDataManager
{
private:

	//Infos for Carousel Data Generation
	int m_numSecondaryObjects = 1000;
	int m_approxLevel = 100;
	
	double m_missDistance = 0.1;
	double m_secondaryCarouselRadius = 10.0;
	//double m_secondaryCarouselPeriod = 10.0;
	double m_distanceBetweenSecondaryCarousels = 1;
	double m_minSecondaryCarouselPeriod = 90.0;
	double m_maxSecondaryCarouselPeriod = 110.0;
	
	string m_directory;

	list<CarouselObjectData> m_secondaryInfos;

	list<thread*> m_workers;
	list<pair<int, int>> m_timeFragmentsQ;
	mutex m_mutexForTimeFragQ;

	list<TCAReport> m_TCAReports;
	mutex m_mutexForTCAReports;

	list<pair<int, int>> m_succeedTimeFragments;
	mutex m_mutexForSucceedTimeFragments;

	list<pair<int, int>> m_failedTimeFragments;
	mutex m_mutexForFailedTimeFragments;

	list<FailSegmentInfo> m_failReports;
	mutex m_mutexForFailReport;

public:

	CarouselDataManager();
	virtual ~CarouselDataManager();

	CarouselDataManager& operator=(const CarouselDataManager& rhs);

	void copy(const CarouselDataManager& rhs);
	void clear();

	//////////////////////////////////////////////////////////////////////////
	// MAIN FUNCTIONS
	//////////////////////////////////////////////////////////////////////////

	void propagate_EVENTSEQ_through_parallel_computation(const int& maxThreads, const double& fragmentSize);
	void generate_PPDB_through_parallel_computation(const int& maxThreads, const double& cutoffDistance);

	//////////////////////////////////////////////////////////////////////////
	// SUB FUNCTIONS
	//////////////////////////////////////////////////////////////////////////

	void read_data_generation_command(const string& filePath);
	void generate_carousel_data();

	virtual void	load_carousel_data_file(const string& filePath);
	virtual void	save_carousel_data_file(const string& filePath);
	string make_carousel_data_file_name(const int& numSecondaryObjects, const int& approxLevel);

	void make_initial_time_fragments(const double& totalPredictionTimeWindow, const double& startTime, const double& fragmentSize);
	void read_time_fragment_file(const string& fileName, const double& predictionStartTime, const double& predictionEndTime);

	void make_thread_for_EVENTSEQ_propagation(CarouselDataEngine* engine, const pair<int, int>& timeFragment);
	void process_failed_EVSEQ_propagation(const CarouselDataEngine* engine, const pair<int, int>& timeFragment);

	void save_succeed_time_fragment(const pair<int, int>& timeFragment);
	void save_failed_engine_info(const CarouselDataEngine* engine, const double& timeFragmentStart, const double& timeFragmentWindow);
	void add_new_time_fragments_to_Q(const list<pair<int, int>>& newTimeFragments);

	void make_thread_for_PPDB_generation(CarouselDataEngine* engine, const pair<int, int>& timeFragment, const double& cutoffDistance);
	void save_TCA_reports_for_engine(const CarouselDataEngine* engine, const double& timeFragmentStart);
	void process_failed_PPDB_generation(const CarouselDataEngine* engine, const pair<int, int>& timeFragment);

	void integrate_TCA_reports(const double& integrationThreshold);
	void make_map_from_pair_to_TCA_reports(map<int, map<int, list<TCAReport>>>& mapFromPairsToTCAReports);
	void integrate_TCA_reports_for_pair(list<TCAReport>& TCAReportsForPair, const double& integrationThreshold);

	void write_result_files_for_EVSEQ_propagation(const double& timeElapsed);
	void write_result_files_for_PPDB_generation(const double& timeElapsed);
	void write_PPDB(const list<TCAReport>& TCAReports, const string& fileName, double timeElapsed);

	bool report_engines_status(const list<pair<CarouselDataEngine*, pair<int, int>>>& engines, const double& timeElapsed) const;

};

