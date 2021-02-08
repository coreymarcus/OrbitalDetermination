#include <iostream> //outputs
#include <boost/numeric/odeint.hpp> //integrator
#include <math.h> // trig functions
#include <Eigen/Dense> // vectors
#include "VehicleState.h" // headerfile containing propagator class
#include "Util.h" // utility functions
#include "MatlabEngine.hpp" // matlab API
#include "MatlabDataArray.hpp" // data things

/* The type of container used to hold the state vector */
typedef std::vector< double > state_type;

//namespaces, don't ever use name space, simply call them aliases
namespace mat = matlab::engine;

int main() {

	//look for matlab sessions
	std::vector<mat::String> matnames = mat::findMATLAB();
	if(matnames.size() == 0){
		std::cout << "Error: MATLAB Shared Session not Found!" << "\n";
		std::cout << "Open MATLAB and run command matlab.engine.shareEngine" << "\n";
	}
	std::cout << mat::convertUTF16StringToUTF8String(matnames[0]) << std::endl;
    std::unique_ptr<mat::MATLABEngine> matlabPtr = mat::connectMATLAB(matnames[0]); 

	//form a propagator object
	VehicleState::Propagator propobj;

	//set grav constant
	propobj.mu_ = 398600.5; // km^3/sec^2

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
	const int N = 238; // aproximate number for two orbits
	Eigen::Matrix<double, 9, N> xhist;
	for (int ii = 0; ii < N; ++ii){

		//propagate
		propobj.Propagate(dt);

		//store state
		xhist.block(0,ii,3,1) = propobj.pos_;
		xhist.block(3,ii,3,1) = propobj.vel_;
		xhist.block(6,ii,3,1) = propobj.GetAccelVector();
	}

	//write data to csv
	Util::Eigen2csv("../data/xhist_HW1.csv",xhist);

	//run the matlab plotting script
	matlabPtr->eval(u"run(\'/home/cm58349/Documents/OrbitalDetermination/src/HW1/plotting.m\')");

	
	
} // main