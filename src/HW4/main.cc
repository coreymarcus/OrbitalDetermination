#include <iostream> //outputs
// #include <boost/numeric/odeint.hpp> //integrator
#include <math.h> // trig functions
#include <Eigen/Dense> // vectors
// #include "VehicleState.h" // headerfile containing propagator class
// #include "Estimators.h" // headerfile containing estimator classes
#include "Util.h" // utility functions
#include <iomanip>      // std::setprecision
#include <unsupported/Eigen/MatrixFunctions>


int main() {

	//initialize position and time
	Eigen::Vector3d pos_ECEF(-28738.3218400000, -30844.0723200000, -6.71800000000000);
	double time_JD = 2458088.500555556;

	//convert position
	Eigen::Vector3d pos_ECI = Util::ECEF2ECI(pos_ECEF, time_JD);


	return 0;
	
	
} // main