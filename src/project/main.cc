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
	pos0[0] = 6984.45711518852;
	pos0[1] = 1612.2547582643;
	pos0[2] = 13.0925904314402;

	Eigen::Vector3d vel0;
	vel0[0] = -1.67667852227336;
	vel0[1] = 7.26143715396544;
	vel0[2] = 0.259889857225218;

	//set observation station ECEF location
	Eigen::Vector3d obs_station1{-6143.584, 1364.250, 1033.743}; //atoll / Kwajalein
	Eigen::Vector3d obs_station2{1907.295, 6030.810, -817.119}; //diego garcia
	Eigen::Vector3d obs_station3{2390.310, -5564.341, 1994.578}; //arecibo

	//observation station measurement noise
	Eigen::MatrixXd R1 = Eigen::MatrixXd::Zero(2,2);
	Eigen::MatrixXd R2 = Eigen::MatrixXd::Zero(2,2);
	Eigen::MatrixXd R3 = Eigen::MatrixXd::Zero(2,2);
	R1(0,0) = pow(10.0/1000.0,2);
	R2(0,0) = pow(5.0/1000.0,2);
	R3(0,0) = pow(10.0/1000.0,2);
	R1(1,1) = pow(0.5/1000000.0,2);
	R2(1,1) = pow(1.0/1000000.0,2);
	R3(1,1) = pow(0.5/1000000.0,2);

	propobj.pos_ = pos0;
	propobj.vel_ = vel0;
	propobj.t_JD_ = Util::JulianDateNatural2JD(2018.0, 3.0, 30.0, 8.0, 55.0, 3.0);

	//convert pos and vel into orbital elements
	propobj.State2OE();

	//set tolerance options
	propobj.abstol_ = 1.0*pow(10.0,-16.0);
	propobj.reltol_ = 3.0*pow(10.0,-14.0);
	propobj.dt_var_ = 0.1;

	//create an estimator
	Estimate::UKF ukf;

	//initial estimate
	Eigen::MatrixXd Phat0 = Eigen::MatrixXd::Zero(6,6);
	Phat0.block(0,0,3,3) = 3.0*Eigen::MatrixXd::Identity(3,3);
	Phat0.block(3,3,3,3) = 0.1*Eigen::MatrixXd::Identity(3,3);
	Eigen::VectorXd xhat0(6);
	xhat0.segment(0,3) = pos0;
	xhat0.segment(3,3) = vel0;
	
	ukf.Phat_ = Phat0;
	ukf.xhat_ = xhat0;
	ukf.k_ = 0.5;
	ukf.n_ = 6;
	ukf.m_ = 2;

	std::cout << "xhat: \n" << ukf.xhat_ << "\n";
	std::cout << "Phat: \n" << ukf.Phat_ << "\n";

	//number of sigma points
	int Nsig = 2*6 + 1;

	// timing
	double dt; //seconds for propagation
	int N = 113; // number of measurements

	//load the measurements
	Eigen::MatrixXd z = Util::LoadDatFile("../data/meas_proj_set1.csv",435,4);

	//process the first measurement outside the loop
	double tof = z(0,2)/c;

	//start at the predicted time
	propobj.t_ = -1.0*tof;

	//initialize a vector of propagation objects for each sigma point
	std::vector<VehicleState::Propagator> propobj_vec(Nsig, propobj);

	//generate sigma points
	ukf.GetSigmaPoints();

	//detirmine which tracking station was used
	int stationID = (int) z(0, 0);

	//get the station position
	Eigen::Vector3d obs_station_iter;
	switch(stationID) {
		case 1:
			obs_station_iter = obs_station1;
			ukf.R_ = R1;
			break;

		case 2:
			obs_station_iter = obs_station2;
			ukf.R_ = R2;
			break;

		case 3:
			obs_station_iter = obs_station3;
			ukf.R_ = R3;
			break;

		default: std::cout << "Error: bad case in measurement \n";
	}

	//initialze predicted measurements for each sigma point
	Eigen::MatrixXd Y = Eigen::MatrixXd::Zero(2,Nsig);

	//cycle through sigma points, writing them to each element of the list
	//and getting a predicted measurement
	for (int i = 0; i < Nsig; ++i)	{
		
		//extract sig state
		Eigen::VectorXd xi = ukf.Xi_.block(0,i,6,1);

		//assign
		propobj_vec[i].pos_ = xi.segment(0,3);
		propobj_vec[i].vel_ = xi.segment(3,3);

		//predicted measurement
		Y.block(0,i,2,1) = propobj_vec[i].GetRangeAndRate(obs_station_iter);
	}

	//use these sigma points to find an estimate
	Eigen::VectorXd ziter = z.block(0,2,1,2).transpose();
	ukf.CalcEstimate(ziter, Y);

	std::cout << "xhat: \n" << ukf.xhat_ << "\n";
	std::cout << "Phat: \n" << ukf.Phat_ << "\n";

	exit(0);

	// propagate starting with the second measurement
	for (int ii = 1; ii < N; ++ii){

		//get this measurement time of flight
		tof = z(ii,2)/c;

		//get the time we need to propagate to get to this measurement
		dt = z(ii,1) - tof - z(ii-1,1);

		//propagate
		propobj.Propagate(dt, false);

		// //detirmine which tracking station was used
		// stationID = (int) z(ii, 0);

		// switch(stationID) {
		// 	case 1:
		// 		ziter = propobj.GetRangeAndRate(obs_station1);
		// 		break;

		// 	case 2:
		// 		ziter = propobj.GetRangeAndRate(obs_station2);
		// 		break;

		// 	case 3:
		// 		ziter = propobj.GetRangeAndRate(obs_station3);
		// 		break;

		// 	default: std::cout << "Error: bad case in measurement \n";
		// }

		//propagate forward to get to the measurement's time of arrival
		propobj.Propagate(tof, false);

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