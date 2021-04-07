#include <iostream> //outputs
#include <boost/numeric/odeint.hpp> //integrator
#include <math.h> // trig functions
#include <Eigen/Dense> // vectors
#include "VehicleState.h" // headerfile containing propagator class
#include "Estimators.h" // headerfile containing estimator classes
#include "Util.h" // utility functions
#include <iomanip>      // std::setprecision
#include <unsupported/Eigen/MatrixFunctions>


int main() {

	//precision
	std::cout << std::setprecision(5);

	//form a propagator object
	VehicleState::Propagator propobj;

	//set physics constants
	propobj.mu_ = 398600.4415; // km^3/sec^2
	propobj.mu_sun_ = 132712440018.0; //km^3/sec^2
	propobj.mu_moon_ = 4902.800066; //km^3/sec^2
	propobj.AU_ = 149597870.7;
	propobj.J2_ = 0.00108248;
	propobj.J3_ = 0.0000025327;
	propobj.Rearth_ = 6378.1363; //km
	propobj.earthrotationspeed_ = 7.292115146706979 * pow(10.0,-5.0); // rad/sec
	propobj.C_D_ = 1.88;
	propobj.A_ = 3.6; // m^2
	propobj.m_ = 2000.0; //kg
	propobj.rho_0_ = 3.614*pow(10.0,-13.0); //kg/m^3
	propobj.r0_ = 700.0 + propobj.Rearth_; //km
	propobj.H_ = 88.6670; //km

	double c = 299792.0; //speed of light km/sec

	//parameters
	propobj.useJ2_ = true;
	propobj.use20x20_ = true;
	propobj.usedrag_ = true;
	propobj.useSRP_ = true;
	propobj.useLuniSolar_ = true; //gravity of sun and moon

	//form a gravity object
	Util::EGM96Grav gravmodel;
	gravmodel.LoadNormCoeffs("../data/egm96_C_normalized.csv", "../data/egm96_S_normalized.csv");
	gravmodel.NormCoeffs2Reg();
	gravmodel.mu_ = propobj.mu_;
	gravmodel.Rearth_ = propobj.Rearth_;
	Eigen::MatrixXd nut80 = Util::LoadDatFile("../data/nut80.csv", 106, 10);
	Eigen::MatrixXd iau1980 = Util::LoadDatFile("../data/iau1980modifiedHW6.csv",15, 4);
	gravmodel.nut80ptr_ = &nut80; 
	gravmodel.iau1980ptr_ = &iau1980;
	propobj.gravmodel_ = &gravmodel;

	//set objects position and velocity
	Eigen::Vector3d pos0;
	pos0[0] = 6984.464593016789876855909824371337890625;
	pos0[1] = 1612.2223725952417225926183164119720458984375;
	pos0[2] = 13.0914313353482452129128432716242969036102294921875;

	Eigen::Vector3d vel0;
	vel0[0] = -1.676644766682913623156991889118216931819915771484375;
	vel0[1] = 7.26144494619244529332036108826287090778350830078125;
	vel0[2] = 0.25988992108511244083501878776587545871734619140625;

	//set observation station ECEF location
	Eigen::Vector3d obs_station1{-6143.584, 1364.250, 1033.743}; //atoll
	Eigen::Vector3d obs_station2{1907.295, 6030.810, -817.119}; //diego garcia
	Eigen::Vector3d obs_station3{2390.310, -5564.341, 1994.578}; //arecibo

	propobj.pos_ = pos0;
	propobj.vel_ = vel0;
	propobj.t_ = 0.0;
	propobj.t_JD_ = Util::JulianDateNatural2JD(2018.0, 2.0, 1.0, 5.0, 0.0, 0.0);

	//initialize STM
	Eigen::MatrixXd STM = Eigen::MatrixXd::Identity(6,6);

	//convert pos and vel into orbital elements
	propobj.State2OE();

	//set tolerance options
	propobj.abstol_ = 1.0*pow(10.0,-16.0);
	propobj.reltol_ = 3.0*pow(10.0,-14.0);
	propobj.dt_var_ = 0.1;

	//create an estimator
	Estimate::UKF ukf;

	//test propagate
	Eigen::MatrixXd Phat0 = Eigen::MatrixXd::Zero(3,3);
	Phat0(0,0) = 3.0;
	Phat0(0,1) = 1.0;
	Phat0(1,0) = 1.0;
	Phat0(1,1) = 10.0;
	Phat0(2,2) = 5.0;
	Phat0(1,2) = 1.5;
	Phat0(2,1) = 1.5;
	Eigen::VectorXd xhat0 = Eigen::VectorXd::Zero(3);
	xhat0(0) = 1.0;
	xhat0(1) = 2.3;
	xhat0(2) = -1.0;
	
	ukf.Phat_ = Phat0;
	ukf.xhat_ = xhat0;
	ukf.k_ = 0.5;
	ukf.n_ = 3;
	ukf.GetSigmaPoints();

	std::cout << "Xi: \n" << ukf.Xi_ << "\n";

	ukf.SigmaPts2Estimate();
	std::cout << xhat0 << "\n" << Phat0 << "\n";
	std::cout << ukf.xhat_ << "\n" << ukf.Phat_ << "\n";

	exit(0);

	// timing
	double dt; //seconds for propagation
	int N = 113; // number of measurements

	//measurements
	Eigen::MatrixXd zbar = Eigen::MatrixXd::Zero(2,N);
	Eigen::VectorXd ziter = propobj.GetRangeAndRate(obs_station1);
	zbar.block(0,0,2,1) = ziter; //first measurement

	//load the true measurements
	Eigen::MatrixXd z_true = Util::LoadDatFile("../data/meas_HW5.csv",113,4);

	// propagate starting with the second measurement
	for (int ii = 1; ii < N; ++ii){

		//get this measurement time of flight
		double tof = z_true(ii,2)/c;

		//get the time we need to propagate to get to this measurement
		dt = z_true(ii,1) - tof - z_true(ii-1,1);

		//propagate
		propobj.Propagate(dt, true);

		//update STM
		STM = propobj.STM_ * STM;

		//detirmine which tracking station was used
		int stationID = (int) z_true(ii, 0);

		switch(stationID) {
			case 1:
				ziter = propobj.GetRangeAndRate(obs_station1);
				break;

			case 2:
				ziter = propobj.GetRangeAndRate(obs_station2);
				break;

			case 3:
				ziter = propobj.GetRangeAndRate(obs_station3);
				break;

			default: std::cout << "Error: bad case in measurement \n";
		}

		//write
		zbar.block(0,ii,2,1) = ziter;

		//propagate forward to get to the measurement's time of arrival
		propobj.Propagate(tof, true);

		//update STM
		STM = propobj.STM_ * STM;

	}

	//true propagation result for 20x20 grav and 6 hours
	Eigen::Vector3d pos_true;

	//20x20 only
	// pos_true[0] = -5153.77192464701;
	// pos_true[1] = -4954.43582198601;
	// pos_true[2] = -144.832817220038;

	//20x20 plus lunisolar
	// pos_true[0] = -5153.7904826402;
	// pos_true[1] = -4954.42147166823;
	// pos_true[2] = -144.825029304757;

	//20x20 lunisolar drag
	// pos_true[0] = -5153.78037301524;
	// pos_true[1] = -4954.43083112784;
	// pos_true[2] = -144.825403424162;

	//20x20 lunisolar drag SRP
	// pos_true[0] = -5153.78219138574;
	// pos_true[1] = -4954.43108360912;
	// pos_true[2] = -144.82549360719;


	return 0;
	
	
} // main