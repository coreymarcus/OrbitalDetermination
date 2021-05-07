#include <iostream> //outputs
#include <boost/numeric/odeint.hpp> //integrator
#include <math.h> // trig functions
#include <Eigen/Dense> // vectors
#include "VehicleState.h" // headerfile containing propagator class
#include "Estimators.h" // headerfile containing estimator classes
#include "Util.h" // utility functions
#include <iomanip>      // std::setprecision
#include <unsupported/Eigen/MatrixFunctions>


int main(int argc, char** argv) {

	//process inputs
	std::string project_case = argv[1];

	//precision
	// std::cout << std::setprecision(5);
	std::cout << std::setprecision(17);

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
	// propobj.C_D_ = 2.0;
	propobj.C_D_ = 1.88;
	// propobj.C_D_ = 1.80;
	// propobj.C_D_ = 1.70;
	// propobj.C_D_ = 1.60;
	// propobj.C_D_ = 1.25;
	// propobj.C_D_ = 0.0;
	// propobj.C_D_ = 18.8;
	// propobj.C_D_ = 188.0;
	propobj.A_ = 22; // m^2
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
	Eigen::MatrixXd iau1980 = Util::LoadDatFile("../data/iau1980modified.csv",17563, 5);
	gravmodel.nut80ptr_ = &nut80; 
	gravmodel.iau1980ptr_ = &iau1980;
	propobj.gravmodel_ = &gravmodel;

	//set objects position and velocity (from assignment document)
	// Eigen::Vector3d pos0;
	// pos0[0] = 6984.45711518852;
	// pos0[1] = 1612.2547582643;
	// pos0[2] = 13.0925904314402;

	// Eigen::Vector3d vel0;
	// vel0[0] = -1.67667852227336;
	// vel0[1] = 7.26143715396544;
	// vel0[2] = 0.259889857225218;

	//set objects position and velocity (from optimization)
	// Eigen::Vector3d pos0;
	// pos0[0] = 6978.83947078333;
	// pos0[1] = 1617.08566078009;
	// pos0[2] = 19.5045324160835;

	// Eigen::Vector3d vel0;
	// vel0[0] = -1.66346314624123;
	// vel0[1] = 7.26036443567597;
	// vel0[2] = 0.270402425183416;

	//set objects position and velocity (from backwards prop)
	Eigen::Vector3d pos0;
	pos0[0] = 6978.6333198140128;
	pos0[1] = 1616.5877041718923;
	pos0[2] = 19.521992036638412;

	Eigen::Vector3d vel0;
	vel0[0] = -1.6630417308551808;
	vel0[1] = 7.2607803252380343;
	vel0[2] = 0.27056529506678861;

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

	//underweight range measurements
	R1(1,1) = 2.0*R1(1,1);
	R2(1,1) = 2.0*R2(1,1);
	R3(1,1) = 2.0*R3(1,1);

	//observation station biases
	Eigen::Vector2d bias1(0.0,0.0);
	Eigen::Vector2d bias2(0.0,0.0);
	Eigen::Vector2d bias3(0.020,0.0);
	// Eigen::Vector2d bias3(0.0,0.0);

	//different cases
	std::string xhat_NAG_filename;
	std::string Phat_NAG_filename;

	//write residuals logic
	bool writeresiduals = false;

	// case A
	if(project_case.compare("A") == 0) {
		std::cout << "Case A\n";
		R1(1,1) = 10000000;
		R2(1,1) = 10000000;
		R3(1,1) = 10000000;
		xhat_NAG_filename = "../data/xhat_A_NAG.csv";
		Phat_NAG_filename = "../data/Phat_A_NAG.csv";
		// var_i = 5.0*pow(10.0,-15.0); //used for NAG1
		// var_i = 1.0*pow(10.0,-19.0); //used for finding bias
	}

	// case B
	if(project_case.compare("B") == 0) {
		std::cout << "Case B\n";
		R1(0,0) = 10000000;
		R2(0,0) = 10000000;
		R3(0,0) = 10000000;
		xhat_NAG_filename = "../data/xhat_B_NAG.csv";
		Phat_NAG_filename = "../data/Phat_B_NAG.csv";
		// var_i = 5.0*pow(10.0,-15.0); //used for NAG1
		// var_i = 1.0*pow(10.0,-19.0); //used for finding bias
	}

	// case C
	if(project_case.compare("C") == 0) {
		std::cout << "Case C\n";
		R2 = 1000000*R2;
		R3 = 1000000*R3;
		xhat_NAG_filename = "../data/xhat_C_NAG.csv";
		Phat_NAG_filename = "../data/Phat_C_NAG.csv";
		// var_i = 5.0*pow(10.0,-15.0); //used for NAG1
		// var_i = 1.0*pow(10.0,-19.0); //used for finding bias
	}

	// case D
	if(project_case.compare("D") == 0) {
		std::cout << "Case D\n";
		R1 = 1000000*R2;
		R3 = 1000000*R3;
		xhat_NAG_filename = "../data/xhat_D_NAG.csv";
		Phat_NAG_filename = "../data/Phat_D_NAG.csv";
		// var_i = 5.0*pow(10.0,-14.0); //used for NAG1 - case D
		// var_i = 1.0*pow(10.0,-19.0); //used for finding bias
	}

	// case E
	if(project_case.compare("E") == 0) {
		std::cout << "Case E\n";
		R1 = 1000000*R2;
		R2 = 1000000*R3;
		xhat_NAG_filename = "../data/xhat_E_NAG.csv";
		Phat_NAG_filename = "../data/Phat_E_NAG.csv";
		// var_i = 5.0*pow(10.0,-15.0); //used for NAG1
		// var_i = 1.0*pow(10.0,-19.0); //used for finding bias
	}

	// case F
	if(project_case.compare("F") == 0) {
		std::cout << "Case F\n";
		xhat_NAG_filename = "../data/xhat_F_NAG.csv";
		Phat_NAG_filename = "../data/Phat_F_NAG.csv";
		// var_i = 5.0*pow(10.0,-15.0); //used for NAG1
		// var_i = 1.0*pow(10.0,-19.0); //used for finding bias
		writeresiduals = true;
	}

	// stations 1 and 2 only
	if(project_case.compare("12") == 0){
		std::cout << "Case 1 and 2 only\n";
		writeresiduals = true;
		xhat_NAG_filename = "../data/xhat_only12.csv";
		Phat_NAG_filename = "../data/Phat_only12.csv";
		// var_i = 1.0*pow(10.0,-19.0);
		R3 = 10000000*Eigen::Matrix2d::Identity(2,2);
		// R1(1,1) = 10000000;
		// R2(1,1) = 10000000;
		// R3(1,1) = 10000000;
	}

	// stations 1 and 3 only
	if(project_case.compare("13") == 0){
		std::cout << "Case 1 and 3 only\n";
		writeresiduals = true;
		xhat_NAG_filename = "../data/xhat_only13.csv";
		Phat_NAG_filename = "../data/Phat_only13.csv";
		// var_i = 1.0*pow(10.0,-19.0);
		R2 = 10000000*Eigen::Matrix2d::Identity(2,2);
		// R1(1,1) = 10000000;
		// R2(1,1) = 10000000;
		// R3(1,1) = 10000000;
	}

	// stations 2 and 3 only
	if(project_case.compare("23") == 0){
		std::cout << "Case 2 and 3 only\n";
		writeresiduals = true;
		xhat_NAG_filename = "../data/xhat_only23.csv";
		Phat_NAG_filename = "../data/Phat_only23.csv";
		// var_i = 1.0*pow(10.0,-19.0);
		R1 = 10000000*Eigen::Matrix2d::Identity(2,2);
		// R1(1,1) = 10000000;
		// R2(1,1) = 10000000;
		// R3(1,1) = 10000000;
	}


	//process noise matrix
	Eigen::MatrixXd Q_sub = Eigen::MatrixXd::Identity(3,3);
	// double var_i = sqrt((1.0/3.0)*pow(25.0/(21600.0*21600.0),2));
	// double var_i = pow(2.0*5.0/(21600.0*21600.0),2.0);
	// double var_i = 5.0*pow(10.0,-15.0); //used for NAG1
	// double var_i = 5.0*pow(10.0,-14.0); //used for NAG1 - case D
	// double var_rad = 1.0*pow(10.0,-16.0); //try 1E-12 * 1E-12
	// double var_in = 1.0*pow(10.0,-17.0);
	// double var_cross = 1.0*pow(10.0,-17.0);
	// double var_rad = 1.0*pow(10.0,-13.0)*pow(10.0,-13.0); //try 1E-12 * 1E-12
	// double var_in = 1.0*pow(10.0,-13.5)*pow(10.0,-13.5);
	// double var_cross = 1.0*pow(10.0,-13.5)*pow(10.0,-13.5);
	// double var_rad = 1.0*pow(10.0,-12.0)*pow(10.0,-12.0); //try 1E-12 * 1E-12
	// double var_in = 1.0*pow(10.0,-12.5)*pow(10.0,-12.5);
	// double var_cross = 1.0*pow(10.0,-12.5)*pow(10.0,-12.5);
	// double var_rad = 1.0*pow(10.0,-11.0)*pow(10.0,-11.0);
	// double var_in = 1.0*pow(10.0,-11.5)*pow(10.0,-11.5);
	// double var_cross = 1.0*pow(10.0,-11.5)*pow(10.0,-11.5);
	double var_rad = 1.0*pow(10.0,-10.0)*pow(10.0,-10.0);
	double var_in = 1.0*pow(10.0,-10.5)*pow(10.0,-10.5);
	double var_cross = 1.0*pow(10.0,-10.5)*pow(10.0,-10.5);
	// double var_rad = 1.0*pow(10.0,-8.0)*pow(10.0,-8.0);
	// double var_in = 1.0*pow(10.0,-8.5)*pow(10.0,-8.5);
	// double var_cross = 1.0*pow(10.0,-8.5)*pow(10.0,-8.5);
	// double var_rad = 1.0*pow(10.0,-7.0)*pow(10.0,-7.0);
	// double var_in = 1.0*pow(10.0,-7.5)*pow(10.0,-7.5);
	// double var_cross = 1.0*pow(10.0,-7.5)*pow(10.0,-7.5);

	//process noise for Cd estimation
	double var_Cd = 1.0*pow(10.0,-10.0);
	std::cout << "var_Cd: " << var_Cd << "\n";


	//construct rest of Q
	Q_sub(0,0) = var_rad;
	Q_sub(1,1) = var_in;
	Q_sub(2,2) = var_cross;
	Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(7,7);

	// timing
	double dt; //seconds for propagation
	int N = 435; // number of measurements for set 1
	// int N = 1289; //number of measurements for set 2
	// int N = 2570; //number of measurements for set 3

	//initialize state for object
	propobj.pos_ = pos0;
	propobj.vel_ = vel0;
	propobj.t_JD_ = Util::JulianDateNatural2JD(2018.0, 3.0, 23.0, 8.0, 55.0, 3.0); //initial epoch
	// double t_dV1 = Util::JulianDateNatural2JD(2018.0, 3.0, 30.0, 8.0, 55.0, 3.0); //dV1
	double t_dV1 = Util::JulianDateNatural2JD(2018.0, 3.0, 24.0, 8.0, 55.0, 3.0); //one day only

	// std::cout << "Natural Julian Date: " << propobj.t_JD_ << "\n";

	//convert pos and vel into orbital elements
	propobj.State2OE();

	//set tolerance options
	propobj.abstol_ = 1.0*pow(10.0,-12.0);
	propobj.reltol_ = 3.0*pow(10.0,-10.0);
	propobj.dt_var_ = 0.1;

	//create an estimator
	Estimate::UKF ukf;

	//initial estimate
	Eigen::MatrixXd Phat0 = Eigen::MatrixXd::Zero(7,7);
	// Phat0.block(0,0,3,3) = 100.0*Eigen::MatrixXd::Identity(3,3);
	// Phat0.block(3,3,3,3) = 0.01*Eigen::MatrixXd::Identity(3,3);

	//values from backwards prop
	Phat0.block(0,0,3,3) = 5.0*pow(10.0,-7.0)*Eigen::MatrixXd::Identity(3,3);
	Phat0.block(3,3,3,3) = 1.0*pow(10.0,-11.0)*Eigen::MatrixXd::Identity(3,3);
	Phat0(6,6) = 0.3;

	Eigen::VectorXd xhat0(7);
	xhat0.segment(0,3) = pos0;
	xhat0.segment(3,3) = vel0;
	xhat0[6] = propobj.C_D_;
	
	ukf.Phat_ = Phat0;
	ukf.xhat_ = xhat0;
	ukf.k_ = 0.5;
	// ukf.k_ = 0.9;
	ukf.n_ = xhat0.size();
	ukf.m_ = 2;

	//number of sigma points
	int Nsig = 2*xhat0.size() + 1;


	//load the measurements
	Eigen::MatrixXd z = Util::LoadDatFile("../data/meas_proj_set1.csv",N,4);
	// Eigen::MatrixXd z = Util::LoadDatFile("../data/meas_proj_set2.csv",N,4);
	// Eigen::MatrixXd z = Util::LoadDatFile("../data/meas_proj_set3.csv",N,4);

	//process the first measurement outside the loop
	double tof = z(0,2)/c;

	//start at the predicted time
	propobj.t_ = -1.0*tof;

	//initialize a vector of propagation objects for each sigma point
	std::vector<VehicleState::Propagator> propobj_vec(Nsig, propobj);

	//initialize data storage
	Eigen::MatrixXd xhat_mat = Eigen::MatrixXd::Zero(7,N);
	Eigen::MatrixXd Phat_mat = Eigen::MatrixXd::Zero(49,N);
	Eigen::MatrixXd Pyy_mat = Eigen::MatrixXd::Zero(4,N);
	Eigen::MatrixXd prefit_res = Eigen::MatrixXd::Zero(2,N);
	Eigen::MatrixXd postfit_res = Eigen::MatrixXd::Zero(2,N);

	//generate sigma points
	ukf.GetSigmaPoints();

	//detirmine which tracking station was used
	int stationID = (int) z(0, 0);

	//get the station position
	Eigen::Vector3d obs_station_iter;
	Eigen::Vector2d bias_iter;
	switch(stationID) {
		case 1:
			obs_station_iter = obs_station1;
			ukf.R_ = R1;
			bias_iter = bias1;

			// std::cout << "Station 1 \n";
			break;

		case 2:
			obs_station_iter = obs_station2;
			ukf.R_ = R2;
			bias_iter = bias2;

			// std::cout << "Station 2 \n";
			break;

		case 3:
			obs_station_iter = obs_station3;
			ukf.R_ = R3;
			bias_iter = bias3;

			// std::cout << "Station 3 \n";
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
		propobj_vec[i].C_D_ = xi[6];

		//predicted measurement
		Y.block(0,i,2,1) = propobj_vec[i].GetRangeAndRate(obs_station_iter, tof) + bias_iter;

		//now that we have the predicted measurement, propagate forward to the current time
		propobj_vec[i].Propagate(tof, false);
	}

	//assign estimate to propobj for residual calculation
	propobj.Propagate(tof,false);
	Eigen::Vector2d prefit_pred = propobj.GetRangeAndRate(obs_station_iter, tof) + bias_iter;

	//use these sigma points to find an estimate
	Eigen::VectorXd ziter = z.block(0,2,1,2).transpose();

	std::cout << "Measurement 1: \n" << ziter << "\n";
	std::cout << "Prediction 1: \n" << prefit_pred << "\n";

	Eigen::MatrixXd Pyy = ukf.CalcEstimate(ziter, Y);

	//assign estimate to propobj for residual calculation
	propobj.pos_ = ukf.xhat_.segment(0,3);
	propobj.vel_ = ukf.xhat_.segment(3,3);
	propobj.C_D_ = ukf.xhat_[6];
	Eigen::Vector2d postfit_pred = propobj.GetRangeAndRate(obs_station_iter, tof) + bias_iter;

	//store data
	xhat_mat.block(0,0,7,1) = ukf.xhat_;
	Eigen::Map<Eigen::VectorXd> Phat_vec(ukf.Phat_.data(), ukf.Phat_.size());
	Eigen::Map<Eigen::VectorXd> Pyy_vec(Pyy.data(), Pyy.size());
	Phat_mat.block(0,0,49,1) = Phat_vec;
	Pyy_mat.block(0,0,4,1) = Pyy_vec;
	prefit_res.block(0,0,2,1) = ziter - prefit_pred - bias_iter;
	postfit_res.block(0,0,2,1) = ziter - postfit_pred - bias_iter;

	std::cout << ziter - prefit_pred << "\n";
	// exit(0);

	//maximum time for propagation
	double maxproptime = 89.3; 

	// propagate starting with the second measurement
	for (int ii = 1; ii < N; ++ii){

		//////////////////// Propagate /////////////////////

		//get this measurement time of flight
		tof = z(ii,2)/c;

		//get the time we need to propagate to get to this measurement
		dt = z(ii,1) - tof - propobj_vec[0].t_;

		//we will propagate no more than maxproptime seconds at a time to avoid issues with process noise
		double N_prop = floor(dt/maxproptime);
		double rem = dt - N_prop*maxproptime;

		std::cout << "dt: " << dt << " N: " << N_prop << " rem: " << rem << "\n";

		//create UKF sigma points
		ukf.GetSigmaPoints();

		//propagate through the intervals
		if(N_prop > 0) {

			for (int k = 0; k < N_prop; ++k) {

				//propagate each sigma point
				for (int j = 0; j < Nsig; ++j){

					//extract sig state
					Eigen::VectorXd xi = ukf.Xi_.block(0,j,6,1);

					//assign
					propobj_vec[j].pos_ = xi.segment(0,3);
					propobj_vec[j].vel_ = xi.segment(3,3);
					propobj_vec[j].C_D_ = xi[6];

					//propagate
					propobj_vec[j].Propagate(maxproptime,false);

					//update sigma point
					ukf.Xi_.block(0,j,3,1) = propobj_vec[j].pos_;
					ukf.Xi_.block(3,j,3,1) = propobj_vec[j].vel_;
				}

				//Update the estimate
				ukf.SigmaPts2Estimate();

				//get orbital elements
				propobj_vec[Nsig - 1].State2OE();

				//create rotation matricies
				double d2r = M_PI/180.0;
				double Ohm = propobj_vec[Nsig - 1].ascend_;
				double w = propobj_vec[Nsig - 1].periap_;
				double inc = propobj_vec[Nsig - 1].i_;
				double theta = propobj_vec[Nsig - 1].nu_;

				std::vector<int> order{3,1,3};
				Eigen::Vector3d angles1(Ohm*d2r, inc*d2r, w*d2r);
				Eigen::Vector3d angles2(-1.0*theta*d2r, 0.0, 0.0);

				Eigen::Matrix3d R_ECI2PQW = Util::Angle2RotM(angles1, order);
				Eigen::Matrix3d R_PQW2RSW = Util::Angle2RotM(angles2, order);
				Eigen::Matrix3d R_total = R_PQW2RSW*R_ECI2PQW;

				Eigen::Matrix3d Q_iter = R_total*Q_sub*R_total.transpose();

				//add process noise
				Q.block(0,0,3,3) = 0.25*pow(maxproptime,4.0)*Q_iter;
				Q.block(0,3,3,3) = 0.5*pow(maxproptime,3.0)*Q_iter;
				Q.block(3,0,3,3) = 0.5*pow(maxproptime,3.0)*Q_iter;
				Q.block(3,3,3,3) = 0.25*pow(maxproptime,2.0)*Q_iter;
				Q(6,6) = pow(maxproptime,2.0)*var_Cd;
				ukf.Phat_ = ukf.Phat_ + Q;

				//Get new sigma points
				ukf.GetSigmaPoints();

			}

		}

		//propagate through the remainder
		for (int j = 0; j < Nsig; ++j){

			//extract sig state
			Eigen::VectorXd xi = ukf.Xi_.block(0,j,6,1);

			//assign
			propobj_vec[j].pos_ = xi.segment(0,3);
			propobj_vec[j].vel_ = xi.segment(3,3);
			propobj_vec[j].C_D_ = xi[6];

			//propagate
			propobj_vec[j].Propagate(rem,false);

			//update sigma point
			ukf.Xi_.block(0,j,3,1) = propobj_vec[j].pos_;
			ukf.Xi_.block(3,j,3,1) = propobj_vec[j].vel_;
		}

		//Update the estimate
		ukf.SigmaPts2Estimate();

		//get orbital elements
		propobj_vec[Nsig - 1].State2OE();

		//create rotation matricies
		double d2r = M_PI/180.0;
		double Ohm = propobj_vec[Nsig - 1].ascend_;
		double w = propobj_vec[Nsig - 1].periap_;
		double inc = propobj_vec[Nsig - 1].i_;
		double theta = propobj_vec[Nsig - 1].nu_;

		std::vector<int> order{3,1,3};
		Eigen::Vector3d angles1(Ohm*d2r, inc*d2r, w*d2r);
		Eigen::Vector3d angles2(-1.0*theta*d2r, 0.0, 0.0);

		Eigen::Matrix3d R_ECI2PQW = Util::Angle2RotM(angles1, order);
		Eigen::Matrix3d R_PQW2RSW = Util::Angle2RotM(angles2, order);
		Eigen::Matrix3d R_total = R_PQW2RSW*R_ECI2PQW;

		Eigen::Matrix3d Q_iter = R_total*Q_sub*R_total.transpose();

		//add process noise
		Q.block(0,0,3,3) = 0.25*pow(rem,4.0)*Q_iter;
		Q.block(0,3,3,3) = 0.5*pow(rem,3.0)*Q_iter;
		Q.block(3,0,3,3) = 0.5*pow(rem,3.0)*Q_iter;
		Q.block(3,3,3,3) = 0.25*pow(rem,2.0)*Q_iter;
		Q(6,6) = pow(rem,2.0)*var_Cd;
		ukf.Phat_ = ukf.Phat_ + Q;

		//Get new sigma points
		ukf.GetSigmaPoints();

		//update time on the propobj
		propobj.t_ = propobj_vec[0].t_;

		//////////////////// Update /////////////////////


		//determine which tracking station was used
		stationID = (int) z(ii, 0);

		switch(stationID) {
			case 1:
				obs_station_iter = obs_station1;
				ukf.R_ = R1;
				bias_iter = bias1;
				break;

			case 2:
				obs_station_iter = obs_station2;
				ukf.R_ = R2;
				bias_iter = bias2;
				break;

			case 3:
				obs_station_iter = obs_station3;
				ukf.R_ = R3;
				bias_iter = bias3;
				break;

			default: std::cout << "Error: bad case in measurement \n";
		}

		//set time appropriately for dopplar shift
		propobj.t_ += tof;

		//assign estimate to propobj for residual calculation
		propobj.pos_ = ukf.xhat_.segment(0,3);
		propobj.vel_ = ukf.xhat_.segment(3,3);
		propobj.C_D_ = ukf.xhat_[6];
		// prefit_pred = propobj.GetRangeAndRate(obs_station_iter, tof) + bias_iter;

		//reset prefit predicted measurement
		prefit_pred[0] = 0.0;
		prefit_pred[1] = 0.0;

		//cycle through sigma points, writing them to each element of the list
		//and getting a predicted measurement
		for (int j = 0; j < Nsig; ++j)	{
			
			//extract sig state
			Eigen::VectorXd xi = ukf.Xi_.block(0,j,6,1);

			//assign
			propobj_vec[j].pos_ = xi.segment(0,3);
			propobj_vec[j].vel_ = xi.segment(3,3);
			propobj_vec[j].C_D_ = xi[6];

			//set time properly for dopplar shift
			propobj_vec[j].t_ += tof;

			//predicted measurement
			Y.block(0,j,2,1) = propobj_vec[j].GetRangeAndRate(obs_station_iter, tof) + bias_iter;

			//reset dopplar correction
			propobj_vec[j].t_ -= tof;

			//prefit predicted measuremnet
			prefit_pred += ukf.w_[j]*Y.block(0,j,2,1);
		}

		// std::cout << "Measurement Block: \n" << Y << "\n";

		//measurement
		ziter = z.block(ii,2,1,2).transpose();

		//perform update
		Pyy = ukf.CalcEstimate(ziter, Y);
		ukf.GetSigmaPoints();

		//get postfit predicted measurement
		postfit_pred[0] = 0.0;
		postfit_pred[1] = 0.0;
		for (int j = 0; j < Nsig; ++j)	{
			
			//extract sig state
			Eigen::VectorXd xi = ukf.Xi_.block(0,j,6,1);

			//assign
			propobj_vec[j].pos_ = xi.segment(0,3);
			propobj_vec[j].vel_ = xi.segment(3,3);
			propobj_vec[j].C_D_ = xi[6];

			//set time properly for dopplar shift
			propobj_vec[j].t_ += tof;

			//predicted measurement
			Y.block(0,j,2,1) = propobj_vec[j].GetRangeAndRate(obs_station_iter, tof) + bias_iter;

			//reset dopplar correction
			propobj_vec[j].t_ -= tof;

			//prefit predicted measuremnet
			postfit_pred += ukf.w_[j]*Y.block(0,j,2,1);
		}

		//assign estimate to propobj for residual calculation
		propobj.pos_ = ukf.xhat_.segment(0,3);
		propobj.vel_ = ukf.xhat_.segment(3,3);
		propobj.C_D_ = ukf.xhat_[6];
		// postfit_pred = propobj.GetRangeAndRate(obs_station_iter, tof) + bias_iter;

		// std::cout << "weights: \n" << ukf.w_ << "\n";
		// std::cout << "postfit_pred: \n" << postfit_pred << "\n";

		//undo timeshift for dopplar
		propobj.t_ -= tof;

		//store data
		xhat_mat.block(0,ii,7,1) = ukf.xhat_;
		Eigen::Map<Eigen::VectorXd> Phat_vec_iter(ukf.Phat_.data(), ukf.Phat_.size());
		Eigen::Map<Eigen::VectorXd> Pyy_vec_iter(Pyy.data(), Pyy.size());
		Phat_mat.block(0,ii,49,1) = Phat_vec_iter;
		Pyy_mat.block(0,ii,4,1) = Pyy_vec_iter;
		prefit_res.block(0,ii,2,1) = ziter - prefit_pred;
		postfit_res.block(0,ii,2,1) = ziter - postfit_pred;

		//////////// Propagate through the TOF ///////////////

		//propagate each sigma point
		for (int j = 0; j < Nsig; ++j){

			//extract sig state
			Eigen::VectorXd xi = ukf.Xi_.block(0,j,6,1);

			//assign
			propobj_vec[j].pos_ = xi.segment(0,3);
			propobj_vec[j].vel_ = xi.segment(3,3);
			propobj_vec[j].C_D_ = xi[6];

			//propagate
			propobj_vec[j].Propagate(tof,false);

			//update sigma point
			ukf.Xi_.block(0,j,3,1) = propobj_vec[j].pos_;
			ukf.Xi_.block(3,j,3,1) = propobj_vec[j].vel_;
		}

		//use sigma points to update estimate
		ukf.SigmaPts2Estimate();

		//////////////////////////////////////////////////////////////////

		std::cout << "Project Case: " << project_case << " Station ID: "<< stationID << "\n";
		std::cout << "postfit: \n" << ziter - postfit_pred << "\n";
		std::cout << "Cd Estimate: " << ukf.xhat_[6] << " StdDev: " << sqrt(ukf.Phat_(6,6)) << "\n";
		// std::cout << "Phat: \n" << ukf.Phat_ << "\n";
		// std::cout << "Q: \n" << Q << "\n";

	}

	//write out data (to avoid dumb mistakes, only write residuals for case F)
	Util::Eigen2csv("../data/xhat_proj.csv", xhat_mat);
	Util::Eigen2csv("../data/Phat_proj.csv", Phat_mat);
	if (writeresiduals) {
		Util::Eigen2csv("../data/Pyy_proj.csv", Pyy_mat);
		Util::Eigen2csv("../data/prefit_res_proj.csv", prefit_res);
		Util::Eigen2csv("../data/postfit_res_proj.csv", postfit_res);
	}
	
	//do the final propagation for the NAG
	double t_total = 24.0*60.0*60.0*(t_dV1 - propobj_vec[0].t_JD_);
	double t_remain = t_total - propobj_vec[0].t_;

	std::cout << "Propating X seconds to delivery time: " << t_remain << "\n";

	//create UKF sigma points
	ukf.GetSigmaPoints();

	//we will propagate no more than 90 seconds at a time to avoid issues with process noise
	double N_prop = floor(t_remain/maxproptime);
	double rem = t_remain - N_prop*maxproptime;

	std::cout << "N: " << N_prop << " rem: " << rem << "\n";

	//create UKF sigma points
	ukf.GetSigmaPoints();

	//propagate through the intervals
	if(N_prop > 0) {

		for (int k = 0; k < N_prop; ++k) {

			//propagate each sigma point
			for (int j = 0; j < Nsig; ++j){

				//extract sig state
				Eigen::VectorXd xi = ukf.Xi_.block(0,j,6,1);

				//assign
				propobj_vec[j].pos_ = xi.segment(0,3);
				propobj_vec[j].vel_ = xi.segment(3,3);
				propobj_vec[j].C_D_ = xi[6];

				//propagate
				propobj_vec[j].Propagate(maxproptime,false);

				//update sigma point
				ukf.Xi_.block(0,j,3,1) = propobj_vec[j].pos_;
				ukf.Xi_.block(3,j,3,1) = propobj_vec[j].vel_;
			}

			//Update the estimate
			ukf.SigmaPts2Estimate();

			//get orbital elements
			propobj_vec[Nsig - 1].State2OE();

			//create rotation matricies
			double d2r = M_PI/180.0;
			double Ohm = propobj_vec[Nsig - 1].ascend_;
			double w = propobj_vec[Nsig - 1].periap_;
			double inc = propobj_vec[Nsig - 1].i_;
			double theta = propobj_vec[Nsig - 1].nu_;

			std::vector<int> order{3,1,3};
			Eigen::Vector3d angles1(Ohm*d2r, inc*d2r, w*d2r);
			Eigen::Vector3d angles2(-1.0*theta*d2r, 0.0, 0.0);

			Eigen::Matrix3d R_ECI2PQW = Util::Angle2RotM(angles1, order);
			Eigen::Matrix3d R_PQW2RSW = Util::Angle2RotM(angles2, order);
			Eigen::Matrix3d R_total = R_PQW2RSW*R_ECI2PQW;

			Eigen::Matrix3d Q_iter = R_total*Q_sub*R_total.transpose();

			//add process noise
			Q.block(0,0,3,3) = 0.25*pow(maxproptime,4.0)*Q_iter;
			Q.block(0,3,3,3) = 0.5*pow(maxproptime,3.0)*Q_iter;
			Q.block(3,0,3,3) = 0.5*pow(maxproptime,3.0)*Q_iter;
			Q.block(3,3,3,3) = 0.25*pow(maxproptime,2.0)*Q_iter;
			Q(6,6) = pow(maxproptime,2.0)*var_Cd;
			ukf.Phat_ = ukf.Phat_ + Q;

			//Get new sigma points
			ukf.GetSigmaPoints();

		}

	}

	//propagate through the remainder
	for (int j = 0; j < Nsig; ++j){

		//extract sig state
		Eigen::VectorXd xi = ukf.Xi_.block(0,j,6,1);

		//assign
		propobj_vec[j].pos_ = xi.segment(0,3);
		propobj_vec[j].vel_ = xi.segment(3,3);
		propobj_vec[j].C_D_ = xi[6];

		//propagate
		propobj_vec[j].Propagate(rem,false);

		//update sigma point
		ukf.Xi_.block(0,j,3,1) = propobj_vec[j].pos_;
		ukf.Xi_.block(3,j,3,1) = propobj_vec[j].vel_;
	}

	//Update the estimate
	ukf.SigmaPts2Estimate();

	//get orbital elements
	propobj_vec[Nsig - 1].State2OE();

	//create rotation matricies
	double d2r = M_PI/180.0;
	double Ohm = propobj_vec[Nsig - 1].ascend_;
	double w = propobj_vec[Nsig - 1].periap_;
	double inc = propobj_vec[Nsig - 1].i_;
	double theta = propobj_vec[Nsig - 1].nu_;

	std::vector<int> order{3,1,3};
	Eigen::Vector3d angles1(Ohm*d2r, inc*d2r, w*d2r);
	Eigen::Vector3d angles2(-1.0*theta*d2r, 0.0, 0.0);

	Eigen::Matrix3d R_ECI2PQW = Util::Angle2RotM(angles1, order);
	Eigen::Matrix3d R_PQW2RSW = Util::Angle2RotM(angles2, order);
	Eigen::Matrix3d R_total = R_PQW2RSW*R_ECI2PQW;

	Eigen::Matrix3d Q_iter = R_total*Q_sub*R_total.transpose();

	//add process noise
	Q.block(0,0,3,3) = 0.25*pow(rem,4.0)*Q_iter;
	Q.block(0,3,3,3) = 0.5*pow(rem,3.0)*Q_iter;
	Q.block(3,0,3,3) = 0.5*pow(rem,3.0)*Q_iter;
	Q.block(3,3,3,3) = 0.25*pow(rem,2.0)*Q_iter;
	Q(6,6) = pow(rem,2.0)*var_Cd;
	ukf.Phat_ = ukf.Phat_ + Q;

	std::cout << "xhat: \n" << ukf.xhat_.segment(0,3) << "\n";
	std::cout << "Phat: \n" << ukf.Phat_.block(0,0,3,3) << "\n";

	Util::Eigen2csv(xhat_NAG_filename, ukf.xhat_.segment(0,6));
	Util::Eigen2csv(Phat_NAG_filename, ukf.Phat_.block(0,0,3,3));

	//update tolerances for backwards propagation
	for (int i = 0; i < Nsig; ++i){
		propobj_vec[i].abstol_ = 1.0*pow(10.0,-16.0);
		propobj_vec[i].reltol_ = 3.0*pow(10.0,-14.0);
		propobj_vec[i].dt_var_ = -0.1;
	}

	//for case F, propagate back to t0 to get an initial estimate	
	// if(project_case.compare("F") == 0 && false){

	// 	std::cout << "Estimating in reverse to get initial conditions... \n";

	// 	for (int ii = N-1; ii >= 0; --ii){

	// 		//////////////////// Propagate /////////////////////

	// 		//get this measurement time of flight
	// 		tof = z(ii,2)/c;

	// 		//difference between current time and time of measurement arrival
	// 		dt = propobj_vec[0].t_ - z(ii,1);

	// 		//magnitude of time to propagate
	// 		double tprop = dt + tof;

	// 		//add process noise
	// 		Q.block(0,0,3,3) = 0.25*pow(tprop,4.0)*Q_sub;
	// 		Q.block(0,3,3,3) = 0.5*pow(tprop,3.0)*Q_sub;
	// 		Q.block(3,0,3,3) = 0.5*pow(tprop,3.0)*Q_sub;
	// 		Q.block(3,3,3,3) = 0.25*pow(tprop,2.0)*Q_sub;
	// 		ukf.Phat_ = ukf.Phat_ + Q;

	// 		//create UKF sigma points
	// 		ukf.GetSigmaPoints();

	// 		//propagate each sigma point
	// 		for (int j = 0; j < Nsig; ++j){

	// 			//extract sig state
	// 			Eigen::VectorXd xi = ukf.Xi_.block(0,j,6,1);

	// 			//assign
	// 			propobj_vec[j].pos_ = xi.segment(0,3);
	// 			propobj_vec[j].vel_ = xi.segment(3,3);

	// 			//propagate backwards in time
	// 			propobj_vec[j].Propagate(-1.0*tprop,false);

	// 			//update sigma point
	// 			ukf.Xi_.block(0,j,3,1) = propobj_vec[j].pos_;
	// 			ukf.Xi_.block(3,j,3,1) = propobj_vec[j].vel_;
	// 		}

	// 		//Update the estimate
	// 		ukf.SigmaPts2Estimate();

	// 		//Get new sigma points
	// 		ukf.GetSigmaPoints();

	// 		//////////////////// Update /////////////////////


	// 		//determine which tracking station was used
	// 		stationID = (int) z(ii, 0);

	// 		switch(stationID) {
	// 			case 1:
	// 				obs_station_iter = obs_station1;
	// 				ukf.R_ = R1;
	// 				bias_iter = bias1;
	// 				break;

	// 			case 2:
	// 				obs_station_iter = obs_station2;
	// 				ukf.R_ = R2;
	// 				bias_iter = bias2;
	// 				break;

	// 			case 3:
	// 				obs_station_iter = obs_station3;
	// 				ukf.R_ = R3;
	// 				bias_iter = bias3;
	// 				break;

	// 			default: std::cout << "Error: bad case in measurement \n";
	// 		}

	// 		//cycle through sigma points, writing them to each element of the list
	// 		//and getting a predicted measurement
	// 		for (int j = 0; j < Nsig; ++j)	{
				
	// 			//extract sig state
	// 			Eigen::VectorXd xi = ukf.Xi_.block(0,j,6,1);

	// 			//assign
	// 			propobj_vec[j].pos_ = xi.segment(0,3);
	// 			propobj_vec[j].vel_ = xi.segment(3,3);

	// 			//predicted measurement
	// 			Y.block(0,j,2,1) = propobj_vec[j].GetRangeAndRate(obs_station_iter, tof) + bias_iter;
	// 		}

	// 		//measurement
	// 		ziter = z.block(ii,2,1,2).transpose();

	// 		//perform update
	// 		Pyy = ukf.CalcEstimate(ziter, Y);
	// 		ukf.GetSigmaPoints();

	// 		//assign estimate to propobj for residual calculation
	// 		propobj.pos_ = ukf.xhat_.segment(0,3);
	// 		propobj.vel_ = ukf.xhat_.segment(3,3);
	// 		propobj.t_ = propobj_vec[0].t_;
	// 		postfit_pred = propobj.GetRangeAndRate(obs_station_iter, tof) + bias_iter;

	// 		//////////////////////////////////////////////////////////////////

	// 		std::cout << "Backwards Prop time: " << propobj_vec[0].t_ << "\n";
	// 		std::cout << "postfit: \n" << ziter - postfit_pred << "\n";
	// 		// std::cout << "Phat: \n" << ukf.Phat_ << "\n";
	// 		// std::cout << "Q: \n" << Q << "\n";

	// 	}

	// 	//prop to t=0

	// 	//propagate each sigma point
	// 	for (int j = 0; j < Nsig; ++j){

	// 		//extract sig state
	// 		Eigen::VectorXd xi = ukf.Xi_.block(0,j,6,1);

	// 		//assign
	// 		propobj_vec[j].pos_ = xi.segment(0,3);
	// 		propobj_vec[j].vel_ = xi.segment(3,3);
	// 		propobj_vec[j].dt_var_ = 0.1;

	// 		//propagate backwards in time
	// 		propobj_vec[j].Propagate(tof,false);

	// 		//update sigma point
	// 		ukf.Xi_.block(0,j,3,1) = propobj_vec[j].pos_;
	// 		ukf.Xi_.block(3,j,3,1) = propobj_vec[j].vel_;
	// 	}

	// 	//Update the estimate
	// 	ukf.SigmaPts2Estimate();

	// 	std::cout << "Initial Time: \n" << propobj_vec[0].t_ << "\n";
	// 	std::cout << "Initial Estimate after Backwards Prop: \n" << ukf.xhat_ << "\n";
	// 	std::cout << "Initial Estimate Covariance after Backwards Prop: \n" << ukf.Phat_ << "\n";
	// }			

	return 0;
	
	
} // main
