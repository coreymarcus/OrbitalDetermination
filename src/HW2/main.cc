#include <iostream> //outputs
#include <boost/numeric/odeint.hpp> //integrator
#include <math.h> // trig functions
#include <Eigen/Dense> // vectors
#include "VehicleState.h" // headerfile containing propagator class
#include "Util.h" // utility functions

/* The type of container used to hold the state vector */
typedef std::vector< double > state_type;

//namespaces, don't ever use name space, simply call them aliases
// namespace mat = matlab::engine;

int main() {

	//form a propagator object
	VehicleState::Propagator propobj;

	//set physics constants
	propobj.mu_ = 398600.4; // km^3/sec^2
	propobj.J2_ = 0.00108248;
	propobj.Rearth_ = 6378.145; //km

	//set objects position and velocity
	Eigen::Vector3d pos0;
	pos0[0] = -2436.45;
	pos0[1] = -2436.45;
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

	//convert OE back to state to verify things working correctly
	propobj.OE2State();

	//set tolerance options
	propobj.abstol_ = 1*pow(10,-20);
	propobj.reltol_ = 1*pow(10,-12);
	propobj.dt_var_ = 0.1;

	//propagate the orbit for two revolutions
	double dt = 20; //seconds for propagation
	const int N = 700; // aproximate number for two orbits
	Eigen::Matrix<double, 12, N> xhist; //state and ang momentum
	Eigen::Matrix<double, 1, N> thist; //time
	Eigen::Matrix<double, 2, N> ehist; //energy
	for (int ii = 0; ii < N; ++ii){

		//propagate
		propobj.Propagate(dt);

		//store state
		xhist.block(0,ii,3,1) = propobj.pos_;
		xhist.block(3,ii,3,1) = propobj.vel_;
		xhist.block(6,ii,3,1) = propobj.GetAccelVector();
		xhist.block(9,ii,3,1) = propobj.GetAngMomVector();
		thist(0,ii) = propobj.t_;
		ehist(0,ii) = propobj.GetKEsp();
		ehist(1,ii) = propobj.GetPEsp();
	}

	//write data to csv
	// Util::Eigen2csv("../data/xhist_HW2.csv",xhist);
	// Util::Eigen2csv("../data/thist_HW2.csv",thist);

	std::cout << "done!" << std::endl;
	
	
} // main