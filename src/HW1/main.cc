#include <iostream> //outputs
#include <boost/numeric/odeint.hpp> //integrator
#include <math.h> // trig functions
#include <Eigen/Dense> // vectors
#include "matplotlibcpp.h" // plotting
#include "VehicleState.h" // headerfile containing propagator class

/* The type of container used to hold the state vector */
typedef std::vector< double > state_type;

//namespaces, don't ever use name space, simply call them aliases
namespace plt = matplotlibcpp;
namespace ode = boost::numeric::odeint;

int main() {

	//form a propagator object
	VehicleState::Propagator propobj;

	//set grav constant
	propobj.mu_ = 398600.5; // km^3/sec^2

	//set objects position and velocity
	Eigen::Vector3d pos0;
	pos0[0] = -2436.45;
	pos0[1] = -2346.45;
	pos0[2] = 6891.037;

	Eigen::Vector3d vel0;
	vel0[0] = 5.088611;
	vel0[1] = -5.088611;
	vel0[2] = 0;

	propobj.pos_ = pos0;
	propobj.vel_ = vel0;
	propobj.t_ = 0.0;

	//convert pos and vel into orbital elements
	propobj.State2OE();

	
	
} // main