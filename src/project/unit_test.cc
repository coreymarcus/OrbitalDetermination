#include <Eigen/Dense> // vectors
#include "VehicleState.h" // headerfile containing propagator class
#include "Util.h" // utility functions
// #include "orbital.h" //Enrico's code for comparison



int main() {

	//////////////// Test ECI to ECEF //////////////////

	//precision
	std::cout << std::setprecision(17);

	//initialize position and time
	Eigen::Vector3d pos_ECEF(-28738.3218400000, -30844.0723200000, -6.71800000000000);

	//target position
	Eigen::Vector3d targ_pos_ECI(19165.44514777874, -37549.06140374086, -41.043609948282580);

	//load some data files
	Eigen::MatrixXd nut80 = Util::LoadDatFile("../data/nut80.csv", 106, 10);
	Eigen::MatrixXd iau1980 = Util::LoadDatFile("../data/iau1980modified.csv",17563, 5);

	//check Julian date algorithm
	double t_JD = Util::JulianDateNatural2JD(2017.0, 12.0, 1.0, 0.0, 0.0, 48.0000384);

	std::cout << "Julian Date: " << t_JD << "\n";

	// get rotation matrix
	std::vector<Eigen::Matrix3d> matvec(5, Eigen::Matrix3d::Identity(3,3));
	Eigen::Matrix3d R = Util::ECEF2ECI(t_JD, &nut80, &iau1980, &matvec);

	//rotate
	Eigen::Vector3d pos_ECI = R*pos_ECEF;

	//error
	std::cout << "Error [cm]: \n" << 100.0*1000.0*(targ_pos_ECI - pos_ECI) << "\n";

	//load the test times
	Eigen::MatrixXd t_R_test = Util::LoadDatFile("../data/t_ECF_test.csv",865,1);

	//load test position
	// Eigen::MatrixXd r_ECF_test = Util::LoadDatFile("../data/r_ECF_test.csv",1,3);
	Eigen::Vector3d pos_ECEF2(6378.136300000000, 6378.136300000000, 6378.136300000000);

	std::cout << "ECEF Test Position: \n" << pos_ECEF2 << "\n";

	//reset julian date to project epoch
	t_JD = Util::JulianDateNatural2JD(2018.0, 3.0, 23.0, 8.0, 55.0, 3.0);

	//initialize ECI matrix
	Eigen::MatrixXd pos_ECI2 = Eigen::MatrixXd::Zero(3,865);

	//perform rotations
	for (int i = 0; i < 865; ++i)	{
		R = Util::ECEF2ECI(t_JD + t_R_test(i,0)/(24.0*60.0*60.0), &nut80, &iau1980, &matvec);
		pos_ECI2.block(0,i,3,1) = R*pos_ECEF2;
	}

	//save data
	Util::Eigen2csv("../data/pred_ECIpos.csv", pos_ECI2);


	// exit(0);s
	////////////////////////////////////////////////////


	///////////////// Test Propagator /////////////////////

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
	propobj.A_ = 15; // m^2
	propobj.m_ = 2000.0; //kg
	propobj.rho_0_ = 3.614*pow(10.0,-13.0); //kg/m^3
	propobj.r0_ = 700.0 + propobj.Rearth_; //km
	propobj.H_ = 88.6670; //km

	//parameters
	propobj.useJ2_ = false;
	propobj.use20x20_ = true;
	propobj.usedrag_ = false;
	propobj.useSRP_ = false;
	propobj.useLuniSolar_ = false; //gravity of sun and moon

	//form a gravity object
	Util::EGM96Grav gravmodel;
	gravmodel.LoadNormCoeffs("../data/egm96_C_normalized.csv", "../data/egm96_S_normalized.csv");
	gravmodel.NormCoeffs2Reg();
	gravmodel.mu_ = propobj.mu_;
	gravmodel.Rearth_ = propobj.Rearth_;
	gravmodel.nut80ptr_ = &nut80; 
	gravmodel.iau1980ptr_ = &iau1980;
	propobj.gravmodel_ = &gravmodel;

	//set objects position and velocity (from lighttime_truth.mat)
	Eigen::Vector3d pos0;
	pos0[0] = 6984.46459301679;
	pos0[1] = 1612.22237259524;
	pos0[2] = 13.0914313353482;

	Eigen::Vector3d vel0;
	vel0[0] = -1.67664476668291;
	vel0[1] = 7.26144494619245;
	vel0[2] =  0.259889921085112;

	//set objects position and velocity (from HW5 assignment)
	// Eigen::Vector3d pos0;
	// pos0[0] = 6990077.798814194;
	// pos0[1] = 1617465.311978378;
	// pos0[2] = 22679.810569245355;
	// pos0 = pos0/1000.0; //convert units

	// Eigen::Vector3d vel0;
	// vel0[0] = -1675.13972506056;
	// vel0[1] = 7273.72441330686;
	// vel0[2] =  252.688512916741;
	// vel0 = vel0/1000.0; //convert units

	propobj.pos_ = pos0;
	propobj.vel_ = vel0;

	propobj.t_JD_ = Util::JulianDateNatural2JD(2018.0, 2.0, 1.0, 5.0, 0.0, 0.0);

	//set tolerance options
	propobj.abstol_ = 1.0*pow(10.0,-16.0);
	propobj.reltol_ = 3.0*pow(10.0,-14.0);
	propobj.dt_var_ = 0.1;

	// propagate
	// propobj.Propagate(21600.0, false);

	//true final position and velocity

	//lighttime_truth.mat
	// Eigen::Vector3d posf_true(-5153.82535324096, -4954.3861758395, -144.823727065817);
	// Eigen::Vector3d velf_true(5.17802136652852, -5.38752135424995, -0.211929233793885);

	// .xlsx file (20x20 only)
	Eigen::Vector3d posf_true(-5153.77192464701, -4954.43582198601, -144.832817220038);
	Eigen::Vector3d velf_true(5.17808108556007, -5.38747088174872, -0.211928868806419);

	//output
	// std::cout << "Final Position:\n" << propobj.pos_ << "\n" << "Error [m]:\n" << 1000.0*(propobj.pos_ - posf_true) << "\n";
	// std::cout << "Final Velocity:\n" << propobj.vel_ << "\n" << "Error [m/sec]:\n" << 1000.0*(propobj.vel_ - velf_true) << "\n";

	//////////////////////////////////////////////////////////////

	//////////////// Test Measurement Predictions /////////////////////

	//simple propagation
	propobj.useJ2_ = true;
	propobj.use20x20_ = true;
	propobj.usedrag_ = true;
	propobj.useSRP_ = true;
	propobj.useLuniSolar_ = true; //gravity of sun and moon

	// load true state and measurements
	Eigen::MatrixXd lighttime_truth = Util::LoadDatFile("../data/lighttime_truth.csv",113, 7);
	Eigen::MatrixXd true_meas = Util::LoadDatFile("../data/LEO_Data_Apparent.csv",113, 4);
	//initialize predicted measurements
	Eigen::MatrixXd pred_meas = Eigen::MatrixXd::Zero(2,113);

	//set observation station ECEF location
	Eigen::Vector3d obs_station1{-6143.584, 1364.250, 1033.743}; //atoll / Kwajalein
	Eigen::Vector3d obs_station2{1907.295, 6030.810, -817.119}; //diego garcia
	Eigen::Vector3d obs_station3{2390.310, -5564.341, 1994.578}; //arecibo

	double c = 299792.0; //speed of light km/sec

	propobj.dt_var_ = -0.1; //only doing backwards props

	//vector for use with choosing station
	Eigen::Vector3d obs_station_iter;

	// cycle through all the measurements, create a predicted measurement at each
	for (int i = 0; i < 113; ++i){
		
		//set time of object
		propobj.t_ = lighttime_truth(i,6);

		//set true state of object
		propobj.pos_ = lighttime_truth.block(i,0,1,3).transpose();
		propobj.vel_ = lighttime_truth.block(i,3,1,3).transpose();

		Eigen::Vector3d pos1;
		Eigen::Vector3d vel1;
		Eigen::Vector3d pos2 = propobj.pos_;
		Eigen::Vector3d vel2 = propobj.vel_;

		//approximate tof
		double tof = true_meas(i,2)/c;

		//propagate backwards in time to measurement
		// propobj.Propagate(-1.0*tof,false);

		pos1 = propobj.pos_;
		vel1 = propobj.vel_;

		std::cout << "Pos1: \n" << pos1 << "\n";
		std::cout << "Pos2: \n" << pos2 << "\n";
		std::cout << "Pos1 + tof*Vel1: \n" << pos1 + tof*vel1 << "\n";
		std::cout << "tof: " << tof << "\n";
		std::cout << "Pos1 + tof*Vel1 - Pos2: \n" << pos1 + tof*vel1 - pos2 << "\n";	

		//determine which tracking station was used
		int stationID = (int) true_meas(i, 0);

		switch(stationID) {
			case 1:
				obs_station_iter = obs_station1;
				break;

			case 2:
				obs_station_iter = obs_station2;
				break;

			case 3:
				obs_station_iter = obs_station3;
				break;

			default: std::cout << "Error: bad case in measurement \n";
		}

		//reset time
		propobj.t_ = lighttime_truth(i,6);

		//approximate measurement
		pred_meas.block(0,i,2,1) = propobj.GetRangeAndRate(obs_station_iter, 1.0*tof);

		// for (int j = 0; j < 5; ++j){

		// 	//approximate measurement
		// 	pred_meas.block(0,i,2,1) = propobj.GetRangeAndRate(obs_station_iter, 1.0*tof);

		// 	std::cout << "Change in ToF: " << tof - pred_meas(0,i)/c << "\n";

		// 	//reset time of flight
		// 	tof = pred_meas(0,i)/c;
		// }
		

		// exit(0);

	}

	//write out the predicted measurements
	Util::Eigen2csv("../data/pred_meas.csv", pred_meas);

	return 0;
}