#include <iostream> //outputs
#include <boost/numeric/odeint.hpp> //integrator
#include <math.h> // trig functions
#include <Eigen/Dense> // vectors
#include "VehicleState.h" // headerfile containing propagator class
// #include "Estimators.h" // headerfile containing estimator classes
#include "Util.h" // utility functions
#include <iomanip>      // std::setprecision
#include <unsupported/Eigen/MatrixFunctions>


int main() {

	//precision
	std::cout << std::setprecision(17);

	//form a propagator object
	VehicleState::Propagator propobj;

	//set physics constants
	propobj.mu_ = 398600.4415; // km^3/sec^2
	propobj.mu_sun_ = 132712440018.0; //km^3/sec^2
	propobj.mu_moon_ = 4902.800066; //km^3/sec^2
	propobj.AU_ = 149597870.7;
	propobj.J2_ = 0.00108248;
	propobj.Rearth_ = 6378.1363; //km
	propobj.earthrotationspeed_ = 7.292115146706979 * pow(10.0,-5.0); // rad/sec
	propobj.C_D_ = 1.88;
	propobj.A_ = 3.6; // m^2
	propobj.m_ = 2000.0; //kg
	propobj.rho_0_ = 3.614*pow(10.0,-13.0); //kg/m^3
	propobj.r0_ = 700.0 + propobj.Rearth_; //km
	propobj.H_ = 88.6670; //km

	//parameters
	propobj.useJ2_ = false;
	propobj.use20x20_ = false;
	propobj.usedrag_ = false;
	propobj.useSRP_ = false;
	propobj.useLuniSolar_ = false; //gravity of sun and moon

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

	propobj.pos_ = pos0;
	propobj.vel_ = vel0;
	propobj.t_ = 0.0;
	propobj.t_JD_ = Util::JulianDateNatural2JD(2018.0, 2.0, 1.0, 5.0, 0.0, 0.0);

	Eigen::Vector3d accel1 = gravmodel.GetGravAccel(pos0, propobj.t_JD_);

	std::cout << "20x20 accel: \n" << accel1 << "\n";

	//convert pos and vel into orbital elements
	propobj.State2OE();

	//set tolerance options
	propobj.abstol_ = 1.0*pow(10.0,-16.0);
	propobj.reltol_ = 3.0*pow(10.0,-14.0);
	propobj.dt_var_ = 0.1;

	//propagate the orbit
	double dt = 60; //seconds for propagation
	int N = 21600/dt; // propagate for desired time
	// int N = 6*60;
	for (int ii = 0; ii < N; ++ii){

		//propagate
		propobj.Propagate(dt, false);

		//update orbital elements
		propobj.State2OE();

	}

	std::cout << "Propagation Result: \n";
	std::cout << "final position:\n" << propobj.pos_ << std::endl;
	std::cout << "final velocity:\n" << propobj.vel_ << std::endl;


	return 0;
	
	
} // main