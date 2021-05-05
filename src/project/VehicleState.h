#ifndef OD_SRC_PROJECT_VEHICLESTATE_H_
#define OD_SRC_PROJECT_VEHICLESTATE_H_

#include <Eigen/Dense>
#include <boost/numeric/odeint.hpp> //integrator
#include <vector> // basic cpp vectors
#include "Util.h" //utility functions

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
		bool usedrag_;
		bool useSRP_; //SRP in propagation?
		bool useLuniSolar_; //should we include the effects of sun and moon
		bool use20x20_; //should we use a spherical harmonics 20x20 gravity model?

		//vehicle position and velocity (ECI) [km] and [km/sec]
		Eigen::Vector3d pos_;
		Eigen::Vector3d vel_;

		//storage for most recent STM calculation by Propagate()
		Eigen::MatrixXd STM_;

		//orbital elements (all angles should be degrees)
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
		double mu_sun_; // sun gravity constant
		double mu_moon_; // moon gravity constant
		double AU_; //length of one atstronomical unit [km]
		double J2_; //higher order gravity parameters
		double J3_;
		double Rearth_; //radius of the earth [km]
		double earthrotationspeed_; //rotation [rad/sec]
		double C_D_; // coefficient of drag
		double A_; // effective area, m^2
		double m_; // vehicle mass, kg
		double rho_0_; // standard air density, kg/m^3
		double r0_; // param for air density calc, km
		double H_; //param for air density calc, km

		//higher order gravity model pointer
		Util::EGM96Grav* gravmodel_;

		//time
		double t_; //seconds since initial time
		double t_JD_; // julian date (UTC) at vehicle state initialization! Should never change

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

		//call the numerical propagator from the current time, t_, to t_ + dt. Option to integrate the STM as well
		void Propagate(double dt, bool intSTM);

		//accessor for the acceleration vector
		Eigen::Vector3d GetAccelVector();

		//accessor for the angular momentum vector
		Eigen::Vector3d GetAngMomVector();

		//accessor for specific kinetic energy
		double GetKEsp();

		//accessor for specific potential energy
		double GetPEsp();

		// () operator function reserved for orbit propagation with odeint
		void operator()( const state_type &x , state_type &dxdt , const double t );

		//obtains the expected range and range rate given the ECEF location of an observation
		Eigen::Vector2d GetRangeAndRate(Eigen::Vector3d pos_station_ecef, double tof);

		//obtain the jacobian of range and range rate measurement given ECEF location of observation
		Eigen::MatrixXd GetRangeAndRateJac(Eigen::Vector3d pos_station_ecef);

	private:

		//private flag for conveniently passing intSTM argument of Propagate() to operator()
		bool intSTM_; 




	}; // class Propagator



} // namespace VehicleState


#endif