#include "CarouselObject.h"

#define _USE_MATH_DEFINES
#include <math.h>

CarouselObject::CarouselObject()
	: DynamicObject()
{

}



CarouselObject::CarouselObject(const int& ID, const double& angularPosition, const rg_Point3D& center, const double& radius, const double& period, const rg_Point3D& uVector, const rg_Point3D& vVector)
	: DynamicObject(ID, DYNAMIC_OBJECT_TYPE::NORMAL_OBJECT, Sphere(rg_Point3D(), OBJECT_SIZE), rg_Point3D())
{
	m_initialAngularPosition = angularPosition;
	m_carouselCenter = center;
	m_carouselRadius = radius;
	m_period = period;
	m_uVector = uVector;
	m_vVector = vVector;
}



CarouselObject::CarouselObject(const CarouselObject& rhs)
{
	copy(rhs);
}



CarouselObject::~CarouselObject()
{
	clear();
}



CarouselObject& CarouselObject::operator=(const CarouselObject& rhs)
{
	if (this == &rhs)
		return *this;

	copy(rhs);
	return *this;
}



void CarouselObject::copy(const CarouselObject& rhs)
{
	DynamicObject::copy(rhs);
	m_initialAngularPosition = rhs.m_initialAngularPosition;
	m_carouselCenter = rhs.m_carouselCenter;
	m_carouselRadius = rhs.m_carouselRadius;
	m_period = rhs.m_period;
	m_uVector = rhs.m_uVector;
	m_vVector = rhs.m_vVector;
}



void CarouselObject::clear()
{
	//Do nothing
}



void CarouselObject::initialize_carousel_object(const int& numSegments, const double& targetTime)
{
	m_trajectoryIndex = 0;
	m_trajectories.clear();

	double anglePerSegment = 2 * M_PI / numSegments;
	for (int i = 0; i < numSegments; i++)
	{
		double time = targetTime + i * m_period / numSegments;
		rg_Point3D coord = calculate_position_of_original_at_time(time);
		m_trajectories.push_back(coord);
	}

	set_sphere_center(m_trajectories.front());
	
	rg_Point3D velocity = (m_trajectories.at(1) - get_sphere_center()) / calculate_time_per_segment();
	set_velocity(velocity);
}



rg_Point3D CarouselObject::calculate_position_of_original_at_time(const double& targetTime) const
{
	double angularVelocity = 2 * M_PI / m_period;
	double angleAtTime = m_initialAngularPosition + angularVelocity * targetTime;
	rg_Point3D positionAtTime = m_carouselCenter + m_carouselRadius*(m_uVector * cos(angleAtTime) + m_vVector * sin(angleAtTime));
	return positionAtTime;
}



double CarouselObject::calculate_approximation_error() const
{
	int approximationLevel = m_trajectories.size();
	double approxError = m_carouselRadius * (1 - cos(M_PI / approximationLevel));
	return approxError;
}



void CarouselObject::change_velocity(const bool& isForward)
{
	++m_trajectoryIndex;
	if (m_trajectoryIndex >= m_trajectories.size())
		m_trajectoryIndex = 0;

	int nextTrajectoryIndex = m_trajectoryIndex + 1;
	if (nextTrajectoryIndex >= m_trajectories.size())
		nextTrajectoryIndex = 0;

	const rg_Point3D& nextTrajectory = m_trajectories.at(nextTrajectoryIndex);
	rg_Point3D velocity = (nextTrajectory - get_sphere_center())/calculate_time_per_segment();
	set_velocity(velocity);
}
