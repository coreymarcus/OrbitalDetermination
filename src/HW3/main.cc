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
	Eigen::MatrixXd truediff(4,N); // x - x_pet
	Eigen::MatrixXd estdiff(4,N); // STM*pet
	Eigen::MatrixXd nomtraj(4,N);
	Eigen::MatrixXd perttraj(4,N);
	Eigen::MatrixXd thist(1,N); //time
	Eigen::MatrixXd STM1;
	for (int ii = 0; ii < N; ++ii){

		//propagate
		propobj.Propagate(dt, false);
		peturbed.Propagate(dt, true);

		//update orbital elements
		propobj.State2OE();
		peturbed.State2OE();

		//estimate deviation from nominal with STM
		STM = peturbed.STM_ * STM;
		estdiff.block(0,ii,6,1) = STM*dev;

		//store first STM
		if (ii == 0) {
			STM1 = STM;
		}

		//store variables
		truediff.block(0,ii,2,1) = propobj.pos_.block(0,0,2,1) - peturbed.pos_.block(0,0,2,1);
		truediff.block(2,ii,2,1) = propobj.vel_.block(0,0,2,1) - peturbed.vel_.block(0,0,2,1);
		nomtraj.block(0,ii,2,1) = propobj.pos_.block(0,0,2,1);
		nomtraj.block(2,ii,2,1) = propobj.vel_.block(0,0,2,1);
		perttraj.block(0,ii,2,1) = peturbed.pos_.block(0,0,2,1);
		perttraj.block(2,ii,2,1) = peturbed.vel_.block(0,0,2,1);


	}

	std::cout << std::setprecision(5);
	std::cout << "Problem 1a: \n" << "i=1\n";
	std::cout << nomtraj.block(0,0,4,1) << "\n" << "i=10\n";
	std::cout << nomtraj.block(0,N-1,4,1) << "\n" << "\n";

	std::cout << "Problem 1b: \n" << "i=1\n";
	std::cout << perttraj.block(0,0,4,1) << "\n" << "i=10\n";
	std::cout << perttraj.block(0,N-1,4,1) << "\n" << "\n";

	std::cout << "STM from 0 to 1: \n" << STM1 << "\n \n";
	std::cout << "STM from 0 to 10: \n" << STM << "\n \n";

	std::cout << "Problem 1-c: \n" << "STM Inverse:\n" << STM.inverse() << "\n";
	std::cout << "STM times Inverse:\n" << STM*STM.inverse() << "\n";

	std::cout << "Problem 1-d: \n" << "(1)\n" << truediff.block(0,0,4,1) << "\n \n" << truediff.block(0,N-1,4,1) << "\n";
	std::cout << "(2)\n" << estdiff.block(0,0,4,1) << "\n \n" << estdiff.block(0,N-1,4,1) << "\n";
	std::cout << "(3)\n" << truediff.block(0,0,4,1) - estdiff.block(0,0,4,1) << "\n \n" << truediff.block(0,N-1,4,1) - estdiff.block(0,N-1,4,1) << "\n";
	
	// Do some problem 2 stuff

	//initialize estimator
	Estimate::KF estKF;

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
	estKF.H_ = H;
	estKF.R_ = W.inverse();
	estKF.Pbar_ = Pbar;
	estKF.mu_ = mu;

	//measurement
	Eigen::VectorXd y(3);
	y(0) = 1.0;
	y(1) = 2.0;
	y(2) = 1.0;

	//perform estimate
	estKF.CalcEstimate(y);
	Eigen::VectorXd errhat = y - H*estKF.xhat_;

	//output
	std::cout << "Problem 2: \n";
	std::cout << "x_hat: " << estKF.xhat_ << std::endl;
	std::cout << "e_hat: " << errhat << std::endl;

	std::cout << "done!\n";
	return 0;
	
	
} // main