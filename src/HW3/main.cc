#include <iostream> //outputs
#include <boost/numeric/odeint.hpp> //integrator
#include <math.h> // trig functions
#include <Eigen/Dense> // vectors
#include "VehicleState.h" // headerfile containing propagator class
#include "Estimators.h" // headerfile containing estimator classes
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

	//form the peturbed object
	Eigen::Vector3d petpos(1.0*pow(10,-6), -1.0*pow(10,-6), 0.0);
	Eigen::Vector3d petvel(1.0*pow(10,-6), 1.0*pow(10,-6), 0.0);
	VehicleState::Propagator peturbed;
	peturbed = propobj;
	peturbed.pos_ = pos0 - petpos;
	peturbed.vel_ = vel0 - petvel;

	//initialize matrix storing STM from t = 0 and peturbation vector
	Eigen::VectorXd dev(6);
	dev.block(0,0,3,1) = petpos;
	dev.block(3,0,3,1) = petvel;
	Eigen::MatrixXd STM = Eigen::MatrixXd::Identity(6,6);

	//propagate the orbit
	double dt = 10.0; //TU for propagation
	const int N = 10; // propagate for N steps
	Eigen::MatrixXd truediff(6,N); // x - x_pet
	Eigen::MatrixXd estdiff(6,N); // STM*pet
	Eigen::MatrixXd thist(1,N); //time
	for (int ii = 0; ii < N; ++ii){

		//propagate
		propobj.Propagate(dt, true);
		peturbed.Propagate(dt, false);

		//update orbital elements
		propobj.State2OE();
		peturbed.State2OE();

		//estimate deviation from nominal with STM
		STM = propobj.STM_ * STM;
		estdiff.block(0,ii,6,1) = STM*dev;

		//store variables
		truediff.block(0,ii,3,1) = propobj.pos_ - peturbed.pos_;
		truediff.block(3,ii,3,1) = propobj.vel_ - peturbed.vel_;


	}

	std::cout << std::setprecision(17);
	std::cout << "Problem 1: \n";
	std::cout << truediff.block(0,N-1,3,1) << std::endl;
	std::cout << estdiff.block(0,N-1,3,1) << std::endl;

	// std::cout << std::setprecision(17);
	// std::cout << "Problem 1: \n";
	// std::cout << "i = 1 position:\n" << xhist.block(0,1,3,1) << std::endl;
	// std::cout << "i = 1 velocity:\n" << xhist.block(3,1,3,1)  << std::endl;
	// std::cout << "i = 10 position:\n" << propobj.pos_ << std::endl;
	// std::cout << "i = 10 velocity:\n" << propobj.vel_ << std::endl;

	// //write data to csv
	// Util::Eigen2csv("../data/xhist_HW3.csv",xhist);
	// Util::Eigen2csv("../data/thist_HW3.csv",thist);

	
	// Do some problem 2 stuff

	//initialize estimator
	Estimate::LS estLS;

	//set properties
	Eigen::MatrixXd H(3,1);
	Eigen::MatrixXd W = Eigen::MatrixXd::Identity(3,3);
	Eigen::MatrixXd Pbar(1,1);
	Eigen::VectorXd mu(1);
	H(0,0) = 1.0;
	H(1,0) = 1.0;
	H(2,0) = 1.0;
	W(0,0) = 2.0;
	Pbar(0,0) = 1/2.0;
	mu(0) = 2.0;
	estLS.H_ = H;
	estLS.R_ = W.inverse();
	estLS.Pbar_ = Pbar;
	estLS.mu_ = mu;

	//measurement
	Eigen::VectorXd y(3);
	y(0) = 1.0;
	y(1) = 2.0;
	y(2) = 1.0;

	//perform estimate
	estLS.CalcEstimate(y);

	//output
	std::cout << "\n Problem 2: \n";
	std::cout << "x_hat: " << estLS.xhat_ << std::endl;
	std::cout << "e_hat: " << y - H*estLS.xhat_ << std::endl;
	
	
} // main