#include <iostream> //outputs
#include <boost/numeric/odeint.hpp> //integrator
#include <math.h> // trig functions
#include <Eigen/Dense> // vectors
#include "VehicleState.h" // headerfile containing propagator class
// #include "Estimators.h" // headerfile containing estimator classes
#include "Util.h" // utility functions
#include <iomanip>      // std::setprecision
#include <unsupported/Eigen/MatrixFunctions>


int main() {

	//precision
	std::cout << std::setprecision(17);

	//form a propagator object
	VehicleState::Propagator propobj;

	//set physics constants
	propobj.mu_ = 398600.4415; // km^3/sec^2
	propobj.mu_sun_ = 132712440018.0; //km^3/sec^2
	propobj.mu_moon_ = 4902.800066; //km^3/sec^2
	propobj.AU_ = 149597870.7;
	propobj.J2_ = 0.00108248;
	propobj.Rearth_ = 6378.1363; //km
	propobj.earthrotationspeed_ = 7.292115146706979 * pow(10.0,-5.0); // rad/sec
	propobj.C_D_ = 1.88;
	propobj.A_ = 3.6; // m^2
	propobj.m_ = 2000.0; //kg
	propobj.rho_0_ = 3.614*pow(10.0,-13.0); //kg/m^3
	propobj.r0_ = 700.0 + propobj.Rearth_; //km
	propobj.H_ = 88.6670; //km

	//parameters
	propobj.useJ2_ = false;
	propobj.usedrag_ = false;
	propobj.useSRP_ = false;
	propobj.useLuniSolar_ = false; //gravity of sun and moon

	//set objects position and velocity
	Eigen::Vector3d pos0;
	pos0[0] = 6990.077798814194;
	pos0[1] = 1617.465311978378;
	pos0[2] = 22.679810569245355;

	Eigen::Vector3d vel0;
	vel0[0] = -1.67513972506056;
	vel0[1] = 7.27372441330686;
	vel0[2] = 0.252688512916741;

	propobj.pos_ = pos0;
	propobj.vel_ = vel0;
	propobj.t_ = 0.0;
	propobj.t_JD_ = Util::JulianDateNatural2JD(2018.0, 2.0, 1.0, 5.0, 0.0, 0.0);

	//convert pos and vel into orbital elements
	propobj.State2OE();

	//set tolerance options
	propobj.abstol_ = 1.0*pow(10.0,-16.0);
	propobj.reltol_ = 3.0*pow(10.0,-14.0);
	propobj.dt_var_ = 0.1;

	//propagate the orbit
	double dt = 60; //seconds for propagation
	// int N = 21600/dt; // propagate for desired time
	int N = 6*60;
	for (int ii = 0; ii < N; ++ii){

		//propagate
		propobj.Propagate(dt, false);

		//update orbital elements
		propobj.State2OE();

	}

	std::cout << "Propagation Result: \n";
	std::cout << "final position:\n" << propobj.pos_ << std::endl;
	std::cout << "final velocity:\n" << propobj.vel_ << std::endl;


	return 0;
	
	
} // main