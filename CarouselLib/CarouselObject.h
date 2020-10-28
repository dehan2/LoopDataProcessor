#pragma once

#include "rg_Point3D.h"
#include "DynamicObject.h"
#include "constForDynamicVD3DLib.h"
#include "ConstForCarouselLib.h"
#include <vector>

using namespace std;

class CarouselObject : public DynamicObject
{
protected:
	double m_initialAngularPosition = 0.0;

	rg_Point3D m_carouselCenter;
	double m_carouselRadius = 0.0;
	double m_period = 0.0;

	rg_Point3D m_uVector;
	rg_Point3D m_vVector;

	int m_trajectoryIndex = 0;
	vector<rg_Point3D> m_trajectories;


public:
	CarouselObject();
	CarouselObject(const int& ID, const double& angularPosition, const rg_Point3D& center, const double& radius, const double& period, const rg_Point3D& uVector, const rg_Point3D& vVector);
	CarouselObject(const CarouselObject& rhs);

	virtual ~CarouselObject();

	CarouselObject& operator=(const CarouselObject& rhs);

	virtual void copy(const CarouselObject& rhs);
	virtual void clear();

	inline const double& get_initial_angular_position() const { return m_initialAngularPosition; }
	inline const rg_Point3D& get_carousel_center()	const { return m_carouselCenter; }
	inline const double& get_carousel_radius() const { return m_carouselRadius; }
	inline const double& get_period() const { return m_period; }
	inline const rg_Point3D& get_uVector() const { return m_uVector; }
	inline const rg_Point3D& get_vVector() const { return m_vVector; }
	inline const int& get_trajectory_index() const { return m_trajectoryIndex; }
	inline const vector<rg_Point3D>& get_trajectories() const { return m_trajectories; }

	inline void set_initial_angular_position(const double& angularPosition) { m_initialAngularPosition = angularPosition; }
	inline void set_carousel_center(const rg_Point3D& carouselCenter)	{ m_carouselCenter = carouselCenter; }
	inline void set_carousel_radius(const double& carouselRadius) { m_carouselRadius = carouselRadius; }
	inline void set_period(const double& period) { m_period = period; }
	inline void set_uVector(const rg_Point3D& uVector) { m_uVector = uVector; }
	inline void set_vVector(const rg_Point3D& vVector) { m_vVector = vVector; }
	inline void set_trajectory_index(const int& trajectoryIndex) { m_trajectoryIndex = trajectoryIndex; }

	void initialize_carousel_object(const int& numSegments, const double& targetTime = 0.0);
	
	//inherited from DynamicObject
	virtual rg_Point3D	calculate_position_of_original_at_time(const double& targetTime) const;
	virtual double	calculate_approximation_error() const;

	virtual void change_velocity(const bool& isForward);
	inline double calculate_time_per_segment() const { return m_period / m_trajectories.size(); }
};

