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

	//precision
	std::cout << std::setprecision(17);

	//initialize position and time
	Eigen::Vector3d pos_ECEF(-28738.3218400000, -30844.0723200000, -6.71800000000000);
	double time_JD = 2458088.50055556;

	//convert position
	Eigen::Vector3d pos_ECI = Util::ECEF2ECI(pos_ECEF, time_JD);

	//final position
	std::cout << "ECI Position: \n" << pos_ECI << "\n";

	//target position
	Eigen::Vector3d targ_pos_ECI(19165.44514777874, -37549.06140374086, -41.043609948282580);

	//error
	std::cout << "Error [m]: \n" << pow(10,3)*(targ_pos_ECI - pos_ECI) << "\n";


	return 0;
	
	
} // main