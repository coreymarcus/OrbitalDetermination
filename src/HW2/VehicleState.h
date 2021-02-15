#ifndef OD_SRC_HW1_VEHICLESTATE_H_
#define OD_SRC_HW1_VEHICLESTATE_H_

#include <Eigen/Dense>
#include <boost/numeric/odeint.hpp> //integrator
#include <vector> // basic cpp vectors

namespace VehicleState{

	// type definitions for the propagator
	namespace ode = boost::numeric::odeint;
	typedef std::vector< double > state_type;
	typedef ode::runge_kutta_dopri5< state_type > error_stepper_type;
	typedef ode::controlled_runge_kutta< error_stepper_type > controlled_stepper_type;

	//class to store and propagate the state of a vehicle
	class Propagator{

	public:

		//parameters
		bool useJ2_;

		//vehicle position and velocity (ECI) [km] and [km/sec]
		Eigen::Vector3d pos_;
		Eigen::Vector3d vel_;

		//orbital elements
		double a_; //semi-major axis
		double e_; //eccentricity
		double i_; //inclination
		double ascend_; //longitude of ascending node
		double periap_; //argument of periapsis
		double nu_; //true anomaly
		double P_; //orbital period
		double T_p_; // time of perigee passage 

		//physics constants
		double mu_; //gravity constant [km^3/sec]
		double J2_; //higher order gravity parameters
		double Rearth_; //radius of the earth [km]

		//time
		double t_;

		//ode integrator elements
		double dt_var_; //initial guess for a variable step size
		double reltol_; //relative tolerance
		double abstol_; //absolute tolerance		

		// Constructor
		Propagator();

		//Use orbital elements to update the state vector
		void OE2State();

		//Use state vector to update orbital elements
		void State2OE();

		//call the numerical propagator from the current time, t_, to t_ + dt
		void Propagate(double dt);

		//accessor for the acceleration vector
		Eigen::Vector3d GetAccelVector();

		//accessor for the angular momentum vector
		Eigen::Vector3d GetAngMomVector();

		//accessor for specific kinetic energy
		double GetKEsp();

		//accessor for specific potential energy
		double GetPEsp();

		// () operator function reserved for orbit propagation with odeint
		void operator()( const state_type &x , state_type &dxdt , const double /* t */ );

	private:




	}; // class Propagator



} // namespace VehicleState


#endif