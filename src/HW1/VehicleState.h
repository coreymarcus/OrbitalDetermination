#ifndef OD_SRC_HW1_VEHICLESTATE_H_
#define OD_SRC_HW1_VEHICLESTATE_H_

#include <Eigen/Dense>

namespace VehicleState{

	//class to store and propagate the state of a vehicle
	class Propagator{

	public:

		//vehicle position and velocity (ECI)
		Eigen::Vector3d pos_;
		Eigen::Vector3d vel_;

		//orbital elements
		double a_; //semi-major axis
		double e_; //eccentricity
		double i_; //inclination
		double ascend_; //longitude of ascending node
		double periap_; //argument of periapsis
		double nu_; //true anomaly

		//physics constants
		double mu_; //gravity constant [km^3/sec]

		//time
		double t_;

		// Constructor
		Propagator();

		//Use orbital elements to update the state vector
		void OE2State();

		//Use state vector to update orbital elements
		void State2OE();

		//call the numerical propagtor from the current time, t_, to t_ + dt
		void Propagate(double dt);


	}; // class Propagator



} // namespace VehicleState


#endif