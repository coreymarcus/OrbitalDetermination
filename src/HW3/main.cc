#include <iostream> //outputs
#include <boost/numeric/odeint.hpp> //integrator
#include <math.h> // trig functions
#include <Eigen/Dense> // vectors
#include "VehicleState.h" // headerfile containing propagator class
#include "Util.h" // utility functions
#include <iomanip>      // std::setprecision
#include <unsupported/Eigen/MatrixFunctions>

/* The type of container used to hold the state vector */
typedef std::vector< double > state_type;

//namespaces, don't ever use name space, simply call them aliases
// namespace mat = matlab::engine;

int main() {

	//form a propagator object
	VehicleState::Propagator propobj;

	//set physics constants
	propobj.mu_ = 1.0; // dimensionless
	propobj.J2_ = 0.0;
	propobj.Rearth_ = 0.0;
	propobj.earthrotationspeed_ = 0.0;
	propobj.C_D_ = 0.0;
	propobj.A_ = 0.0;
	propobj.m_ = 0.0;
	propobj.rho_0_ = 0.0;
	propobj.r0_ = 0.0;
	propobj.H_ = 0.0;

	//parameters
	propobj.useJ2_ = false;
	propobj.usedrag_ = false;

	//set objects position and velocity
	Eigen::Vector3d pos0;
	pos0[0] = 1.0;
	pos0[1] = 0.0;
	pos0[2] = 0.0;

	Eigen::Vector3d vel0;
	vel0[0] = 0.0;
	vel0[1] = 1.0;
	vel0[2] = 0.0;

	propobj.pos_ = pos0;
	propobj.vel_ = vel0;
	propobj.t_ = 0.0;

	//convert pos and vel into orbital elements
	propobj.State2OE();

	//set tolerance options
	propobj.abstol_ = 1*pow(10,-16);
	propobj.reltol_ = 3*pow(10,-14);
	propobj.dt_var_ = 0.1;

	//propagate the orbit
	double dt = 10.0; //TU for propagation
	const int N = 1; // propagate for 11 steps
	Eigen::MatrixXd xhist(6,N); //state history
	Eigen::MatrixXd thist(1,N); //time
	for (int ii = 0; ii < N; ++ii){

		//propagate
		propobj.Propagate(dt, true);

		//update orbital elements
		propobj.State2OE();

		//store variables
		xhist.block(0,ii,3,1) = propobj.pos_;
		xhist.block(3,ii,3,1) = propobj.pos_;

	}

	std::cout << std::setprecision(17);
	std::cout << "Problem 1: \n";
	std::cout << propobj.STM_ << std::endl;

	// std::cout << std::setprecision(17);
	// std::cout << "Problem 1: \n";
	// std::cout << "i = 1 position:\n" << xhist.block(0,1,3,1) << std::endl;
	// std::cout << "i = 1 velocity:\n" << xhist.block(3,1,3,1)  << std::endl;
	// std::cout << "i = 10 position:\n" << propobj.pos_ << std::endl;
	// std::cout << "i = 10 velocity:\n" << propobj.vel_ << std::endl;

	// //write data to csv
	// Util::Eigen2csv("../data/xhist_HW3.csv",xhist);
	// Util::Eigen2csv("../data/thist_HW3.csv",thist);

	std::cout << "done!" << std::endl;
	
	
} // main