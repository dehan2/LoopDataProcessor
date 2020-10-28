#include "CarouselDataManager.h"

#include <random>

#define _USE_MATH_DEFINES
#include <math.h>

#ifdef _WIN32
#include <Windows.h>
#else
#include <unistd.h>
#endif

#include <iostream>
#include <iomanip>

using namespace std;


CarouselDataManager::CarouselDataManager()
{

}



CarouselDataManager::~CarouselDataManager()
{
	clear();
}



CarouselDataManager& CarouselDataManager::operator=(const CarouselDataManager& rhs)
{
	if (this == &rhs)
		return *this;

	copy(rhs);
	return *this;
}



void CarouselDataManager::copy(const CarouselDataManager& rhs)
{

}



void CarouselDataManager::clear()
{

}



void CarouselDataManager::propagate_EVENTSEQ_through_parallel_computation(const int& maxThreads, const double& fragmentSize)
{
	generate_carousel_data();
	double timeWindow = m_numSecondaryObjects;

	make_initial_time_fragments(timeWindow, 0, fragmentSize);

	int engineID = 0;

	Clock::time_point begin = Clock::now();
	Clock::time_point end;
	while (!m_timeFragmentsQ.empty())
	{
		list<pair<CarouselDataEngine*, pair<int, int>>> enginesOnWorking;

		while (enginesOnWorking.size() < maxThreads
			&& !m_timeFragmentsQ.empty())
		{
			CarouselDataEngine* engine = new CarouselDataEngine(engineID++);
			enginesOnWorking.push_back({ engine, m_timeFragmentsQ.front() });
			m_timeFragmentsQ.pop_front();
		}


		for (auto& enginePair : enginesOnWorking)
		{
			thread* initThread = new thread(&CarouselDataManager::make_thread_for_EVENTSEQ_propagation, this, enginePair.first, enginePair.second);
			m_workers.emplace_back(initThread);
		}

		end = Clock::now();
		double timeElapsed = chrono::duration_cast<chrono::milliseconds>(end - begin).count();

		bool isAllFinish = false;
		while (isAllFinish == false)
		{
			isAllFinish = report_engines_status(enginesOnWorking, timeElapsed);

#ifdef _WIN32
			Sleep(1000);
#else
			sleep(1);
#endif
		}


		for (auto& worker : m_workers)
		{
			worker->join();
			delete worker;
		}
		m_workers.clear();

		for (auto& enginePair : enginesOnWorking)
		{
			enginePair.first->clear();
			delete enginePair.first;
		}

		cout << "Threads are joined - Initialization finish" << endl;
	}

	end = Clock::now();
	double timeElapsed = chrono::duration_cast<chrono::milliseconds>(end - begin).count();
	write_result_files_for_EVSEQ_propagation(timeElapsed);
}



void CarouselDataManager::generate_PPDB_through_parallel_computation(const int& maxThreads, const double& cutoffDistance)
{
	string carouselDataFileName = make_carousel_data_file_name(m_numSecondaryObjects, m_approxLevel);
	load_carousel_data_file(carouselDataFileName);

	double timeWindow = m_numSecondaryObjects;
	read_time_fragment_file("SucceedTimeFragments.txt", 0.0, timeWindow);

	int engineID = 0;

	Clock::time_point begin = Clock::now();
	Clock::time_point end;
	while (!m_timeFragmentsQ.empty())
	{
		list<pair<CarouselDataEngine*, pair<int, int>>> enginesOnWorking;

		while (enginesOnWorking.size() < maxThreads
			&& !m_timeFragmentsQ.empty())
		{
			CarouselDataEngine* engine = new CarouselDataEngine(engineID++);
			enginesOnWorking.push_back({ engine, m_timeFragmentsQ.front() });
			m_timeFragmentsQ.pop_front();
		}


		for (auto& enginePair : enginesOnWorking)
		{
			thread* initThread = new thread(&CarouselDataManager::make_thread_for_PPDB_generation, this, enginePair.first, enginePair.second, cutoffDistance);
			m_workers.emplace_back(initThread);
		}

		end = Clock::now();
		double timeElapsed = chrono::duration_cast<chrono::milliseconds>(end - begin).count();

		bool isAllFinish = false;
		while (isAllFinish == false)
		{
			isAllFinish = report_engines_status(enginesOnWorking, timeElapsed);

#ifdef _WIN32
			Sleep(1000);
#else
			sleep(1);
#endif
		}


		for (auto& worker : m_workers)
		{
			worker->join();
			delete worker;
		}
		m_workers.clear();

		for (auto& enginePair : enginesOnWorking)
		{
			enginePair.first->clear();
			delete enginePair.first;
		}

		cout << "Threads are joined - Initialization finish" << endl;
	}

	end = Clock::now();
	double timeElapsed = chrono::duration_cast<chrono::milliseconds>(end - begin).count();
	write_result_files_for_PPDB_generation(timeElapsed);
}



void CarouselDataManager::generate_carousel_data()
{
	cout << "Generate Carousel data" << endl;

	random_device rd;
	mt19937 generator(rd());
	uniform_real_distribution<double> distribution(m_minSecondaryCarouselPeriod, m_maxSecondaryCarouselPeriod);

	double primaryCarouselRadius = m_distanceBetweenSecondaryCarousels * m_numSecondaryObjects * 0.5 / M_PI;

	double distanceToSecondaryCarouselCenter = primaryCarouselRadius + m_missDistance + m_secondaryCarouselRadius;
	cout << "Distance to Sec Car Center: " << distanceToSecondaryCarouselCenter << endl;

	for (int i = 0; i < m_numSecondaryObjects; i++)
	{
		double secondaryPeriod = distribution(generator);
		
		CarouselObjectData secondaryInfo;
		secondaryInfo.ID = i + START_ID_OF_OBJECT + 1;
		secondaryInfo.period = secondaryPeriod;
		secondaryInfo.radius = m_secondaryCarouselRadius;

		m_secondaryInfos.push_back(secondaryInfo);
	}
}



void CarouselDataManager::load_carousel_data_file(const string & filePath)
{
	ifstream fin;
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

			m_secondaryInfos.push_back(secondaryInfo);
		}
	}

	fin.close();
}



void CarouselDataManager::save_carousel_data_file(const string& filePath)
{
	ofstream fout(filePath);
	fout << "% Carousel Data\n";
	fout << "% ApproxLevel\tnumSecondaryObjects\tmissDistance\tSecCarRadius\tMinSecCarPeriod\tMaxSecCarPeriod\tDistanceBtwSecCars\n";
	fout << m_approxLevel << "\t" << m_numSecondaryObjects << "\t" << m_missDistance << "\t" << m_secondaryCarouselRadius << "\t" << m_minSecondaryCarouselPeriod << "\t" << m_maxSecondaryCarouselPeriod << "\t" << m_distanceBetweenSecondaryCarousels << "\n";
	for (auto& secondaryInfo : m_secondaryInfos)
	{
		fout << secondaryInfo.ID << "\t" << secondaryInfo.radius << "\t" << secondaryInfo.period << "\n";
	}
	fout.close();
}



string CarouselDataManager::make_carousel_data_file_name(const int& numSecondaryObjects, const int& approxLevel)
{
	string carouselDataFileName = string("CD_") + to_string(numSecondaryObjects) + string("_") + to_string(approxLevel) + (".cdf");
	return carouselDataFileName;
}



void CarouselDataManager::make_initial_time_fragments(const double& totalPredictionTimeWindow, const double& startTime, const double& fragmentSize)
{
	int scaledCurrTime = startTime * TIME_FRAGMENT_SCALE;
	int scaledFragmentSize = fragmentSize * TIME_FRAGMENT_SCALE;
	int scaledPTWindow = totalPredictionTimeWindow * TIME_FRAGMENT_SCALE;
	while (scaledCurrTime < scaledPTWindow)
	{
		m_timeFragmentsQ.push_back({ scaledCurrTime, scaledFragmentSize });
		scaledCurrTime += scaledFragmentSize;
	}
}



void CarouselDataManager::read_time_fragment_file(const string& fileName, const double& predictionStartTime, const double& predictionEndTime)
{
	int scaledStartTime = predictionStartTime * TIME_FRAGMENT_SCALE;
	int scaledEndTime = predictionEndTime * TIME_FRAGMENT_SCALE;

	ifstream fin;
	fin.open(fileName);

	if (fin.is_open())
	{
		while (!fin.eof())
		{
			char lineData[256];
			fin.getline(lineData, 256);

			if (lineData[0] != '%')
			{
				char* context;
				string delimiter = "\t";

				string token = strtok_s(lineData, delimiter.c_str(), &context);
				int timeFragmentStart = stoi(token);

				token = strtok_s(NULL, delimiter.c_str(), &context);
				int timeFragmentSize = stoi(token);

				if (timeFragmentStart >= scaledStartTime && timeFragmentStart <= scaledEndTime)
					m_timeFragmentsQ.push_back({ timeFragmentStart, timeFragmentSize });
			}
		}
	}

	fin.close();
}



void CarouselDataManager::make_thread_for_EVENTSEQ_propagation(CarouselDataEngine* engine, const pair<int, int>& timeFragment)
{
	engine->set_approximation_level(m_approxLevel);
	engine->set_miss_distance(m_missDistance);
	engine->set_num_secondary_objects(m_numSecondaryObjects);
	engine->set_distance_between_secondary_carousels(m_distanceBetweenSecondaryCarousels);
	engine->set_secondary_carousel_radius(m_secondaryCarouselRadius);

	const double currTime = timeFragment.first * TIME_FRAGMENT_SCALE_REVERSE;
	const double timeFragmentWindow = timeFragment.second * TIME_FRAGMENT_SCALE_REVERSE;
	engine->initialize_carousel_engine(m_secondaryInfos, m_approxLevel, currTime, timeFragmentWindow);
	
	string qtfFilePath = engine->make_qtf_file_name_for_carousel_data(m_numSecondaryObjects, timeFragment.first, timeFragment.second);
	engine->construct_initial_VD_and_save_as_qtf_file(qtfFilePath);
	
	/*bool isQTFAvailable = engine->construct_initial_VD_from_QTF_file(qtfFilePath);
	if (isQTFAvailable == false)
	{
		engine->construct_initial_VD_and_save_as_qtf_file(qtfFilePath);
	}*/

	if (PRINT_PROGRESS_LEVEL >= 1)
		cout << "Start EVSEQ propagation: [" << timeFragment.first << ", " << timeFragment.second << "]" << endl;

	engine->propagate_EVENTSEQ(0);

	if (engine->get_status() == DVDOS_STATUS::FAIL)
	{
		if (PRINT_PROGRESS_LEVEL >= 1)
			cout << "EVSEQ propagation Failed: " << engine->make_error_message() << endl;

		process_failed_EVSEQ_propagation(engine, timeFragment);
		return;
	}

	if (PRINT_PROGRESS_LEVEL >= 1)
		cout << "Start Validation" << endl;

	engine->validate_EVENTSEQ_by_go_back_to_start();

	if (engine->get_status() == DVDOS_STATUS::COMPLETE)
	{
		if (PRINT_PROGRESS_LEVEL >= 1)
			cout << "Validation Finish" << endl;

		engine->save_event_sequence_for_carousel_data(m_numSecondaryObjects, timeFragment.first, timeFragment.second);
		save_succeed_time_fragment(timeFragment);
	}
	else
	{
		if (PRINT_PROGRESS_LEVEL >= 1)
			cout << "Validation Failed: " << engine->make_error_message() << endl;

		process_failed_EVSEQ_propagation(engine, timeFragment);
	}
}



void CarouselDataManager::process_failed_EVSEQ_propagation(const CarouselDataEngine* engine, const pair<int, int>& timeFragment)
{
	save_failed_engine_info(engine, timeFragment.first * TIME_FRAGMENT_SCALE_REVERSE, timeFragment.second * TIME_FRAGMENT_SCALE_REVERSE);

	int timeFragmentWindow = timeFragment.second;

	//Divide by 2, 3, 5
	int divider = 0;
	if (timeFragmentWindow % 2 == 0)
	{
		divider = 2;
	}
	else if (timeFragmentWindow % 3 == 0)
	{
		divider = 3;
	}
	else if (timeFragmentWindow % 5 == 0)
	{
		divider = 5;
	}
	else
	{
		cout << "Fail to divide: " << timeFragmentWindow * TIME_FRAGMENT_SCALE_REVERSE << endl;
	}


	if (divider != 0)
	{
		list<pair<int, int>> newTimeFragments;
		int newTimeFragmentWindow = timeFragmentWindow / divider;
		for (int i = 0; i < divider; i++)
		{
			int newTimeFragmentStart = timeFragment.first + newTimeFragmentWindow * i;
			newTimeFragments.push_back(make_pair(newTimeFragmentStart, newTimeFragmentWindow));
		}

		add_new_time_fragments_to_Q(newTimeFragments);
	}
}



void CarouselDataManager::save_succeed_time_fragment(const pair<int, int>& timeFragment)
{
	lock_guard<mutex> lock(m_mutexForSucceedTimeFragments);
	m_succeedTimeFragments.push_back(timeFragment);
}



void CarouselDataManager::save_failed_engine_info(const CarouselDataEngine* engine, const double& timeFragmentStart, const double& timeFragmentWindow)
{
	lock_guard<mutex> lock(m_mutexForFailReport);

	FailSegmentInfo failInfo;
	failInfo.startTime = timeFragmentStart;
	failInfo.timeWindow = timeFragmentWindow;
	failInfo.processedTime = engine->get_error_report().time;
	failInfo.message = engine->make_error_message();

	m_failReports.push_back(failInfo);
}



void CarouselDataManager::add_new_time_fragments_to_Q(const list<pair<int, int>>& newTimeFragments)
{
	lock_guard<mutex> lock(m_mutexForTimeFragQ);
	for (auto& timeFrag : newTimeFragments)
	{
		m_timeFragmentsQ.push_back(timeFrag);
	}
}



void CarouselDataManager::make_thread_for_PPDB_generation(CarouselDataEngine* engine, const pair<int, int>& timeFragment, const double& cutoffDistance)
{
	engine->set_approximation_level(m_approxLevel);
	engine->set_miss_distance(m_missDistance);
	engine->set_num_secondary_objects(m_numSecondaryObjects);
	engine->set_distance_between_secondary_carousels(m_distanceBetweenSecondaryCarousels);
	engine->set_secondary_carousel_radius(m_secondaryCarouselRadius);

	cout << "Step1: Initialization" << endl;
	Clock::time_point begin = Clock::now();

	const double currTime = timeFragment.first * TIME_FRAGMENT_SCALE_REVERSE;
	const double timeFragmentWindow = timeFragment.second * TIME_FRAGMENT_SCALE_REVERSE;

	engine->initialize_carousel_engine(m_secondaryInfos, m_approxLevel, currTime, timeFragmentWindow);
	
	string eventSeqFileName = engine->make_event_sequence_file_name_for_carousel_data(m_numSecondaryObjects, timeFragment.first, timeFragment.second);
	engine->load_event_sequence_file(eventSeqFileName);

	string qtfFilePath = engine->make_qtf_file_name_for_carousel_data(m_numSecondaryObjects, timeFragment.first, timeFragment.second);
	bool isQTFAvailable = engine->construct_initial_VD_from_QTF_file(qtfFilePath);
	if (isQTFAvailable == false)
	{
		cout << "No QTF is available: " << qtfFilePath << endl;
		process_failed_PPDB_generation(engine, timeFragment);
		return;
	}

	if (PRINT_PROGRESS_LEVEL >= 1)
		cout << "Start PPDB generation" << endl;

	double terminateCondition = 1.0e-6;
	engine->set_terminating_condition_for_DCA_scan(terminateCondition);
	engine->generate_pairwise_proximity_database(cutoffDistance);

	Clock::time_point end = Clock::now();
	double timeElapsed = chrono::duration_cast<chrono::nanoseconds>(end - begin).count() * 1.0e-6;
	cout << "Generate PPDB Finish: " << timeElapsed << "ms" << endl;

	/*const list<SearchInterval>& searchIntervals = engine->get_search_intervals();

	int index = 0;
	for (auto& searchInterval : searchIntervals)
	{
		int primaryID = searchInterval.primary->get_ID();
		int secondaryID = searchInterval.secondary->get_ID();
		if (primaryID > secondaryID)
			swap(primaryID, secondaryID);

		cout << "CA[ " << ++index << " ]: pair [ " << primaryID << ", " << secondaryID << " ], TCA: " << searchInterval.timeOfClosestApproach << ", distance: " << searchInterval.minDistance;
		cout << ", SI: [ " << searchInterval.startTime << ", " << searchInterval.endTime << " ]" << endl;
	}*/



	if (engine->get_status() == DVDOS_STATUS::COMPLETE)
	{
		if (PRINT_PROGRESS_LEVEL >= 1)
			cout << "PPDB Generation Finish: " << timeElapsed << endl;

		save_TCA_reports_for_engine(engine, timeFragment.first * TIME_FRAGMENT_SCALE_REVERSE);
		save_succeed_time_fragment(timeFragment);
	}
	else
	{
		if (PRINT_PROGRESS_LEVEL >= 1)
			cout << "PPDB Generation Failed: " << engine->make_error_message() << endl;

		process_failed_PPDB_generation(engine, timeFragment);
	}
}



void CarouselDataManager::save_TCA_reports_for_engine(const CarouselDataEngine* engine, const double& timeFragmentStart)
{
	lock_guard<mutex> lock(m_mutexForTCAReports);
	const list<SearchInterval>& foundSearchIntervals = engine->get_search_intervals();
	for (auto& searchInterval : foundSearchIntervals)
	{
		TCAReport report;

		report.primaryID = searchInterval.primary->get_ID();
		report.secondaryID = searchInterval.secondary->get_ID();
		report.distanceOfClosestApproach = searchInterval.minDistance;
		report.timeOfClosestApproach = searchInterval.timeOfClosestApproach + timeFragmentStart;
		report.closeApproachEnteringTime = searchInterval.startTime + timeFragmentStart;
		report.closeApproachLeavingTime = searchInterval.endTime + timeFragmentStart;
		m_TCAReports.push_back(report);
	}
}



void CarouselDataManager::process_failed_PPDB_generation(const CarouselDataEngine* engine, const pair<int, int>& timeFragment)
{
	save_failed_engine_info(engine, timeFragment.first * TIME_FRAGMENT_SCALE_REVERSE, timeFragment.second * TIME_FRAGMENT_SCALE_REVERSE);
}



bool compare_TCA_reports_by_time_order(const TCAReport& lhs, const TCAReport& rhs)
{
	return lhs.timeOfClosestApproach < rhs.timeOfClosestApproach;
}



void CarouselDataManager::integrate_TCA_reports(const double& integrationThreshold)
{
	map<int, map<int, list<TCAReport>>> mapFromPairsToTCAReports;
	make_map_from_pair_to_TCA_reports(mapFromPairsToTCAReports);

	m_TCAReports.clear();

	for (auto& primary : mapFromPairsToTCAReports)
	{
		for (auto& secondary : primary.second)
		{
			list<TCAReport>& reportsForPair = secondary.second;
			integrate_TCA_reports_for_pair(reportsForPair, integrationThreshold);
		}
	}

	m_TCAReports.sort(compare_TCA_reports_by_time_order);
}



void CarouselDataManager::make_map_from_pair_to_TCA_reports(map<int, map<int, list<TCAReport>>>& mapFromPairsToTCAReports)
{
	for (auto& report : m_TCAReports)
	{
		auto itForPrimary = mapFromPairsToTCAReports.find(report.primaryID);
		if (itForPrimary == mapFromPairsToTCAReports.end())
		{
			mapFromPairsToTCAReports[report.primaryID] = map<int, list<TCAReport>>();
			mapFromPairsToTCAReports.at(report.primaryID)[report.secondaryID] = list<TCAReport>();
			mapFromPairsToTCAReports.at(report.primaryID).at(report.secondaryID).push_back(report);
		}
		else
		{
			map<int, list<TCAReport>>& mapFromSecondaryToTCAReports = (*itForPrimary).second;
			auto itForSecondary = mapFromSecondaryToTCAReports.find(report.secondaryID);
			if (itForSecondary == mapFromSecondaryToTCAReports.end())
			{
				mapFromSecondaryToTCAReports[report.secondaryID] = list<TCAReport>();
				mapFromSecondaryToTCAReports.at(report.secondaryID).push_back(report);
			}
			else
			{
				list<TCAReport>& TCAReportsForPair = (*itForSecondary).second;
				TCAReportsForPair.push_back(report);
			}
		}
	}
}



void CarouselDataManager::integrate_TCA_reports_for_pair(list<TCAReport>& TCAReportsForPair, const double& integrationThreshold)
{
	list<TCAReport> integratedReports;

	TCAReportsForPair.sort(compare_TCA_reports_by_time_order);

	auto itForTCAReport = TCAReportsForPair.begin();
	while (itForTCAReport != TCAReportsForPair.end())
	{
		TCAReport integratedReport = (*itForTCAReport);
		auto itForCandidate = next(itForTCAReport);
		while (itForCandidate != TCAReportsForPair.end())
		{
			if (((*itForCandidate).closeApproachEnteringTime - integratedReport.closeApproachLeavingTime) < integrationThreshold)
			{
				//Integrate two report
				integratedReport.closeApproachLeavingTime = (*itForCandidate).closeApproachLeavingTime;
				if ((*itForCandidate).distanceOfClosestApproach < integratedReport.distanceOfClosestApproach)
				{
					integratedReport.distanceOfClosestApproach = (*itForCandidate).distanceOfClosestApproach;
					integratedReport.timeOfClosestApproach = (*itForCandidate).timeOfClosestApproach;
				}
				itForCandidate = TCAReportsForPair.erase(itForCandidate);
			}
			else
			{
				break;
			}
		}

		m_TCAReports.push_back(integratedReport);
		itForTCAReport++;
	}
}



bool compare_time_fragments_in_ascending_order(const pair<int, int>& lhs, const pair<int, int>& rhs)
{
	return lhs.first < rhs.first;
}



void CarouselDataManager::write_result_files_for_EVSEQ_propagation(const double& timeElapsed)
{
	cout << "EVSEQ propagation finish: " << timeElapsed << endl;

	string carouselDataFileName = make_carousel_data_file_name(m_numSecondaryObjects, m_approxLevel);
	save_carousel_data_file(carouselDataFileName);


	m_succeedTimeFragments.sort(compare_time_fragments_in_ascending_order);

	string summaryFileName = "Summary-EventSeqGen_" + to_string(m_numSecondaryObjects) + string("_") + to_string(m_approxLevel) + string(".txt");

	ofstream summaryFile(summaryFileName);
	summaryFile << fixed << setprecision(2);
	summaryFile << "Total Elapsed time: " << timeElapsed << "ms" << "\n\n";

	summaryFile << "Processed Time Fragments" << "\n\n";

	for (auto& timeFragments : m_succeedTimeFragments)
		summaryFile << timeFragments.first << "\t" << timeFragments.second << "\n";

	summaryFile << "\n\n Failed Time Fragments" << "\n\n";
	for (auto& failReport : m_failReports)
		summaryFile << failReport.startTime << "\t" << failReport.timeWindow << "\t" << failReport.processedTime << "\t" << failReport.message << "\n";
	summaryFile.close();


	ofstream timeFragInfoFile("SucceedTimeFragments.txt");
	timeFragInfoFile << "% Total Time Fragments: " << m_succeedTimeFragments.size();
	timeFragInfoFile << fixed << setprecision(2);
	for (auto& timeFragments : m_succeedTimeFragments)
	{
		timeFragInfoFile << "\n" << timeFragments.first << "\t" << timeFragments.second;
	}
	timeFragInfoFile.close();
}



void CarouselDataManager::write_result_files_for_PPDB_generation(const double& timeElapsed)
{
	cout << "Scan Close Approach Event Finish: " << timeElapsed << endl;
	string summaryFileName = "Summary-ScanCA_" + to_string(m_numSecondaryObjects) + string("_") + to_string(m_approxLevel) + string(".txt");

	ofstream summaryFile(summaryFileName);
	summaryFile << fixed << setprecision(2);
	summaryFile << "Total Elapsed time: " << timeElapsed << "ms" << "\n\n";

	summaryFile << "Processed Time Fragments" << "\n\n";

	for (auto& timeFragments : m_succeedTimeFragments)
		summaryFile << timeFragments.first << "\t" << timeFragments.second << "\n";

	summaryFile << "\n\n Failed Time Fragments" << "\n\n";
	for (auto& failReport : m_failReports)
		summaryFile << failReport.startTime << "\t" << failReport.timeWindow << "\t" << failReport.processedTime << "\t" << failReport.message << "\n";
	summaryFile.close();

	integrate_TCA_reports(TCA_INTEGRATION_THRESHOLD);

	string PPDBFileName = string("PPDB_") + to_string(m_numSecondaryObjects) + string("_") + to_string(m_approxLevel) + string(".txt");
	write_PPDB(m_TCAReports, PPDBFileName, timeElapsed);
}



void CarouselDataManager::write_PPDB(const list<TCAReport>& TCAReports, const string& fileName, double timeElapsed)
{
	ofstream resultFile(fileName);
	resultFile << fixed << setprecision(3);

	resultFile << "% Close Approach Under " << 2*m_missDistance << "km:\n%\n";
	resultFile << "% Time Elapsed: " << timeElapsed * 1.0e-3 << " sec\n%\n";

	resultFile << "% PrimaryID\tSecondaryID\tMinDistance(km)\tYear\tMon\tDay\tHour\tMin\tSec\n%\n";

	for (auto& report : TCAReports)
	{
		resultFile << "\n";
		resultFile << report.primaryID << "\t" << report.secondaryID << "\t";
		resultFile << report.distanceOfClosestApproach << "\t" << report.timeOfClosestApproach << "\t" << report.closeApproachEnteringTime << "\t" << report.closeApproachLeavingTime;
	}

	resultFile.close();
}



bool CarouselDataManager::report_engines_status(const list<pair<CarouselDataEngine*, pair<int, int>>>& engines, const double& timeElapsed) const
{
	cout << "Time elapsed: " << timeElapsed << "ms" << endl;
	bool isAllFinish = true;
	for (auto& enginePair : engines)
	{
		CarouselDataEngine* engine = enginePair.first;
		const pair<int, int>& timeFragment = enginePair.second;
		cout << "[" << timeFragment.first * TIME_FRAGMENT_SCALE_REVERSE << ", " << timeFragment.second * TIME_FRAGMENT_SCALE_REVERSE << "]: ";
		cout << engine->make_status_message() << endl;

		if (engine->get_status() != DVDOS_STATUS::COMPLETE && engine->get_status() != DVDOS_STATUS::FAIL)
		{
			isAllFinish = false;
		}
	}
	cout << endl;

	return isAllFinish;
}

